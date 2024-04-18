!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
module comm_mh_specind_mod
  use comm_signal_mod
  implicit none



contains


  subroutine sample_template_mh(outdir, cpar, handle, handle_noise)
    implicit none
    character(len=*),               intent(in)    :: outdir
    type(planck_rng),               intent(inout) :: handle, handle_noise
    type(comm_params) :: cpar

    ! Given a component, propose an amplitude to shift the global amplitude of any component.

  end subroutine sample_template_mh

  subroutine sample_mbbtab_mh(outdir, cpar, handle, handle_noise)
    implicit none
    character(len=*),               intent(in)    :: outdir
    type(planck_rng),               intent(inout) :: handle, handle_noise
    type(comm_params) :: cpar

    real(dp)     :: chisq, my_chisq, chisq_old, chisq_new, chisq_prop, mval, mval_0
    integer(i4b) :: band, ierr, i, j, k, l, pol, pix
    logical(lgt)  :: include_comp, reject, todo
    character(len=512) :: tokens(10), str_buff, operation
    class(comm_comp),   pointer           :: c => null()
    class(comm_map), pointer              :: invN_res => null(), map => null(), sig => null(), res => null()



    ! Loop over sampling groups

    do l = 1, cpar%mcmc_num_user_samp_groups


       mval_0 = -1d0
       k = 0

       c => compList
       do while (associated(c))
          k = k + 1
          select type (c)
          type is (comm_MBBtab_comp)
            if (maxval(c%theta_steplen(c%npar+1:,l)) > 0) then
               mval = maxval(c%theta_steplen(3:,l))
               mval_0 = max(mval, mval_0)
               exit
            end if
          end select
          c => c%nextComp()
       end do
       call mpi_bcast(mval_0, 1, MPI_DOUBLE_PRECISION, &
         & 0, data(1)%info%comm, ierr)

       if (mval_0 <= 0d0) then
         cycle
       end if

       call output_FITS_sample(cpar, 1050, .true.)


       ! Calculate initial chisq
       chisq_old = 0d0
       call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_old, &
                          & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

       if (cpar%myid_chain .eq. 0) then
         write(*,*) '| Old chisq is ', chisq_old
       end if

       c => compList
       do while (associated(c))
                       
          select type (c)
          type is (comm_MBBtab_comp)

            if (c%id == k) then

               if (cpar%myid_chain .eq. 0) then
                 write(*,*) trim(c%label)
                 c%SEDtab_buff = c%SEDtab
                 do i = 1, c%ntab
                    c%SEDtab(3,i) = c%SEDtab(3,i) + rand_gauss(handle) * c%theta_steplen(c%npar+i,l)
                 end do
                 write(*,*) 'MBBtab original', c%SEDtab_buff(3,:)
                 write(*,*) 'MBBtab proposal', c%SEDtab(3,:)
               end if


            end if
            call mpi_bcast(c%SEDtab, size(c%SEDtab), MPI_DOUBLE_PRECISION, &
              & 0, data(1)%info%comm, ierr)
            call mpi_bcast(c%SEDtab_buff, size(c%SEDtab), MPI_DOUBLE_PRECISION, &
              & 0, data(1)%info%comm, ierr)

          end select
          
          !go to next compedt
          c => c%nextComp()
       end do


       ! Update mixing matrices
       call update_mixing_matrices(update_F_int=.true.)


       ! Perform component separation
       if (cpar%myid == 0) write(*,*) trim(cpar%mcmc_update_cg_groups(l))
       if (trim(cpar%mcmc_update_cg_groups(l)) == 'none') then
          if (cpar%myid == 0) write(*,*) 'No groups to sample'
       else
          call sample_all_amps_by_CG(cpar, handle, handle_noise, store_buff=.true., cg_groups=cpar%mcmc_update_cg_groups(l))
       end if

       call output_FITS_sample(cpar, 2050, .true.)

       chisq_prop = 0d0
       call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_prop, &
                          & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

       if (cpar%myid_chain .eq. 0) then
         write(*,*) "|    Proposal chisq is ", chisq_prop
         write(*,*) "|    Delta chi^2 is    ", chisq_prop - chisq_old
       end if

       ! Check MH statistic
       reject = log(rand_uni(handle)) > (chisq_old - chisq_prop)/2
       call mpi_bcast(reject, 1, MPI_LOGICAL, 0, data(1)%info%comm, ierr)


       if (reject) then
         if (cpar%myid_chain == 0) then
           write(*,*) '| '
           write(*,*) '| MH step rejected, returning to original tabulated values.'
           write(*,*) '| '
         end if


         c => compList
         do while (associated(c))
                         
            select type (c)
            type is (comm_MBBtab_comp)

              if (c%id == k) then
                 if (cpar%myid_chain .eq. 0) then
                   do i = 1, c%ntab
                      c%SEDtab(3,i) = c%SEDtab_buff(3,i)
                   end do
                 end if

                 call mpi_bcast(c%SEDtab, size(c%SEDtab), MPI_DOUBLE_PRECISION, &
                   & 0, data(1)%info%comm, ierr)
              end if

            end select
            
            !go to next component
            c => c%nextComp()
         end do

         ! Update mixing matrices
         call update_mixing_matrices(update_F_int=.true.)


         ! Instead of doing compsep, revert the amplitudes here
         call revert_CG_amps(cpar)

       else
         if (cpar%myid_chain == 0) then
           write(*,*) '| '
           write(*,*) '| MH step accepted'
           write(*,*) '| '
         end if
       end if


     end do


     if (cpar%myid == 0) then
       write(*,*) 'Finished sampling mbbtab'
     end if

  end subroutine sample_mbbtab_mh

  subroutine sample_specind_mh(outdir, cpar, handle, handle_noise)
    implicit none
    character(len=*),               intent(in)    :: outdir
    type(planck_rng),               intent(inout) :: handle, handle_noise
    type(comm_params) :: cpar

    real(dp)     :: chisq, my_chisq, chisq_old, chisq_new, chisq_prop, mval, mval_0
    integer(i4b) :: band, ierr, i, j, k, pol, pix, l
    logical(lgt)  :: include_comp, reject, todo, has_alms
    character(len=512) :: tokens(10), str_buff, operation
    class(comm_comp),   pointer           :: c => null()
    class(comm_map), pointer              :: invN_res => null(), map => null(), sig => null(), res => null()



    ! Loop over sampling groups

    do l = 1, cpar%mcmc_num_user_samp_groups

       mval_0 = -1d0
       k = 0

       c => compList
       do while (associated(c))
          k = k + 1
          if (.not. allocated(c%theta_steplen)) then
            c => c%nextComp()
          else
            select type (c)
            class is (comm_diffuse_comp)
              if (maxval(c%theta_steplen(1:c%npar,l)) > 0) then
                 mval = maxval(c%theta_steplen(1:c%npar,l))
                 mval_0 = max(mval, mval_0)
                 exit
              end if
            end select
            c => c%nextComp()
          end if
       end do
       call mpi_bcast(mval_0, 1, MPI_DOUBLE_PRECISION, &
         & 0, data(1)%info%comm, ierr)

       if (mval_0 <= 0d0) then
         cycle
       end if



       ! Calculate initial chisq
       chisq_old = 0d0
       call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_old, &
                          & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

       if (cpar%myid_chain .eq. 0) then
         write(*,*) '| Old chisq is ', chisq_old
       end if



       c => compList
       do while (associated(c))
         if (.not. allocated(c%theta_steplen)) then
             c => c%nextComp()
             cycle
          end if
                       
          do j = 1, c%npar
             if (c%theta_steplen(j,l) == 0.d0) cycle
             select type (c)
             class is (comm_diffuse_comp)

               if (c%id == k .and. cpar%myid_chain == 0) then

                   write(*,*) trim(c%label), j
                   do pol = 1, c%theta(j)%p%info%nmaps
                      c%theta_pixreg_buff(1,pol,j) = c%theta_pixreg(1,pol,j)
                      c%theta_pixreg(1,pol,j) = c%theta_pixreg(1,pol,j) + rand_gauss(handle) * c%theta_steplen(j,l)
                   end do
                   write(*,*) 'theta_pixreg original', c%theta_pixreg_buff(1,:,j)
                   write(*,*) 'theta_pixreg proposal', c%theta_pixreg(1,:,j)

                   do pol = 1,c%theta(j)%p%info%nmaps
                      do pix = 0, c%theta(j)%p%info%np-1
                         c%theta(j)%p%map(pix,pol) = c%theta_pixreg(pol,c%ind_pixreg_arr(pix,pol,j),j)
                      end do
                   end do

               end if



               call mpi_bcast(c%theta_pixreg_buff, size(c%theta_pixreg_buff), MPI_DOUBLE_PRECISION, &
                 & 0, data(1)%info%comm, ierr)
               call mpi_bcast(c%theta_pixreg,      size(c%theta_pixreg),      MPI_DOUBLE_PRECISION, &
                 & 0, data(1)%info%comm, ierr)

               if (cpar%myid_chain .eq. 0 .and. c%theta(j)%p%info%nalm .ne. 0) then
                 has_alms = .true.
               else
                 has_alms = .false.
               end if
               call mpi_bcast(has_alms, 1,  MPI_LOGICAL, 0, data(1)%info%comm, ierr)
               if (has_alms) call c%theta(j)%p%YtW_scalar()

             end select

          end do
          
          !go to next component
          c => c%nextComp()
              
       end do

       ! Update mixing matrices
       call update_mixing_matrices(update_F_int=.true.)


       ! Perform component separation
       if (cpar%myid == 0) write(*,*) 'Should only sample groups ', trim(cpar%mcmc_update_cg_groups(l))
       if (trim(cpar%mcmc_update_cg_groups(l)) == 'none') then
          if (cpar%myid == 0) write(*,*) 'No groups to sample'
       else
          call sample_all_amps_by_CG(cpar, handle, handle_noise, store_buff=.true., cg_groups=cpar%mcmc_update_cg_groups(l))
       end if



       chisq_prop = 0d0
       call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_prop, &
                          & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

       if (cpar%myid_chain .eq. 0) then
         write(*,*) "|    Proposal chisq is ", chisq_prop
         write(*,*) "|    Delta chi^2 is    ", chisq_prop - chisq_old
       end if

       ! Check MH statistic
       reject = log(rand_uni(handle)) > (chisq_old - chisq_prop)/2
       call mpi_bcast(reject, 1, MPI_LOGICAL, 0, data(1)%info%comm, ierr)


       if (reject) then
         if (cpar%myid_chain == 0) then
           write(*,*) '| '
           write(*,*) '| MH step rejected, returning to original tabulated values.'
           write(*,*) '| '
         end if




         ! Instead of doing compsep, revert the amplitudes here
         c => compList
         do while (associated(c))
            if (c%npar == 0) then
               c => c%nextComp()
               cycle
            end if
                         
            do j = 1, c%npar
               if (c%theta_steplen(j,l) == 0.d0) cycle
               select type (c)
               class is (comm_diffuse_comp)

                 if (c%id == k) then

                    if (cpar%myid_chain .eq. 0) then
                      do pol = 1, c%theta(j)%p%info%nmaps
                         c%theta_pixreg(1,pol,j) = c%theta_pixreg_buff(1,pol,j)
                      end do
                    end if

                    call mpi_bcast(c%theta_pixreg,      size(c%theta_pixreg),      MPI_DOUBLE_PRECISION, &
                      & 0, data(1)%info%comm, ierr)

                    do pol = 1,c%theta(j)%p%info%nmaps
                       do pix = 0, c%theta(j)%p%info%np-1
                          c%theta(j)%p%map(pix,pol) = c%theta_pixreg(pol,c%ind_pixreg_arr(pix,pol,j),j)
                       end do
                    end do


                    if (cpar%myid_chain .eq. 0 .and. c%theta(j)%p%info%nalm .ne. 0) then
                      has_alms = .true.
                    else
                      has_alms = .false.
                    end if
                    call mpi_bcast(has_alms, 1,  MPI_LOGICAL, 0, data(1)%info%comm, ierr)
                    if (has_alms) call c%theta(j)%p%YtW_scalar()
                  end if
               end select

            end do
            
            !go to next component
            c => c%nextComp()
                
         end do

         ! Update mixing matrices
         call update_mixing_matrices(update_F_int=.true.)

         call revert_CG_amps(cpar)

       else
         if (cpar%myid_chain == 0) then
           write(*,*) '| '
           write(*,*) '| MH step accepted'
           write(*,*) '| '
         end if
       end if


     end do


     if (cpar%myid == 0) then
       write(*,*) 'Finished sampling specind'
     end if

  end subroutine sample_specind_mh

  subroutine initialize_mh_mod(cpar)
    !
    !   Initializes parameters for the MH proposals, including
    !    1. The masks
    !    2. The bands for which chi^2 is evaluated
    !    3. The components to sample and their proposals
    !    4. The specific CG groups to do amplitudes sampling on
    !

    ! Need to create a new set of variables that hold the step length
    implicit none
    type(comm_params) :: cpar


    type(comm_mapinfo), pointer :: info => null(), info_def => null(), info_ud
    class(comm_map), pointer :: indmask, mask_ud


    integer(i4b) :: i, j, k, nside, lmax, nmaps, n, m, n_tokens, n_in_group, ierr
    integer(i4b), allocatable, dimension(:,:)  :: bands_to_sample
    real(dp) :: sigma
    logical(lgt) :: pol
    class(comm_comp),   pointer           :: c => null()


    character(len=128) :: tokens(100), comp_tokens(2), comp_names(2), comp_bands(2), wire_to(2), wire_from(2)

    allocate(cpar%mcmc_group_bands_indices(cpar%mcmc_num_user_samp_groups, numband))

    cpar%mcmc_group_bands_indices = 0


    ! Need to add an argument to sample_all_amps_by_CG that includes this list
    ! if (cpar%myid == 0) then
    !     do i = 1, cpar%mcmc_num_user_samp_groups
    !        write(*,*) i, trim(cpar%mcmc_update_cg_groups(i)), trim(cpar%mcmc_samp_groups(i))
    !        if (trim(cpar%mcmc_update_cg_groups(i)) .ne. 'none') then
    !          call get_tokens(cpar%mcmc_update_cg_groups(i), ',', tokens, n)
    !          if (n == 0) then
    !            write(*,*) "Something is wrong with your CG groups"
    !            write(*,*) trim(cpar%mcmc_update_cg_groups(i))
    !          end if
    !        end if
    !     end do
    ! end if


    if (cpar%myid == 0) then
      do i = 1, cpar%mcmc_num_user_samp_groups
        call get_tokens(cpar%mcmc_samp_group_bands(i), ',', tokens, n)
        if (n == 0) then
          write(*,*) 'You have misformatted MCMC_SAMPLING_GROUP_CHISQ_BANDS'
          stop
        else
          n_in_group = 0
          do j = 1, n
            k = findloc(cpar%ds_label, tokens(j), dim=1)
            if (k .ne. 0) then
              n_in_group = n_in_group + 1
              cpar%mcmc_group_bands_indices(i, n_in_group) = k
            end if
          end do

        end if
      end do
    end if

    call mpi_bcast(cpar%mcmc_group_bands_indices, size(cpar%mcmc_group_bands_indices), &
                   & MPI_INTEGER, 0, data(1)%info%comm, ierr)




    !if (cpar%myid == 0) then
        do i = 1, cpar%mcmc_num_user_samp_groups


          if (cpar%myid == 0) write(*,*) trim(cpar%mcmc_samp_groups(i))


           call get_tokens(cpar%mcmc_samp_groups(i), ',', tokens, n_tokens)
           do j = 1, n_tokens
             if (trim(tokens(1)) == 'none') exit

             call get_tokens(tokens(j), '>', comp_tokens, n)
             if (n == 2) then
               call get_tokens(comp_tokens(1), ':', wire_from, num=n)
               call get_tokens(comp_tokens(2), ':', wire_to, num=n)


             else if (n == 1) then
                 !write(*,*) ''
                 !write(*,*) 'No wiring'
                 call get_tokens(tokens(j), '%', comp_tokens)
                 read(comp_tokens(2), *) sigma

                 call get_tokens(comp_tokens(1), ':', comp_names)


                 if (trim(comp_names(1)) == 'gain') then
                   if (cpar%myid == 0) write(*,*) 'We are sampling gain', sigma

                 else if (trim(comp_names(2)) == 'scale') then
                   if (cpar%myid == 0) write(*,*) 'We are scaling the amplitude'

                 else if (comp_names(2)(1:3) == 'tab') then
                   ! write(*,*) 'We are dealing with tabulated values: ', trim(comp_names(2)), sigma
                   ! write(*,*) trim(comp_names(2)(4:))
                   ! Currently, there is no distinguishing between 05a and 05b.
                   ! We need to make the SEDtab more extendable, like having a dictionary for each band.
                   !c%theta_steplen(cind, i) = 
                   ! Parse the comp_names, get the index
                   call get_tokens(comp_names(2), '@', comp_bands)

                   ! We currently have all bands being 5a/b, 05a/b, etc. In principle we need to do this better, but for now we will
                   ! just get everything except the last index.

                   m = len_trim(comp_bands(2))
                   read(comp_bands(2)(1:m-1), *) m
                  

                   !write(*,*) 'comp names ', trim(comp_names(1))

                   c => compList
                   do while (associated(c))
                      if (trim(c%label) == trim(comp_names(1))) then
                        !       (beta+T+ntab, n_mcmc_samp_groups)
                        !write(*,*) 'set sigma_SEDtab  to ',sigma
                        !write(*,*) i, 2+m, c%theta_steplen(2+m,i)
                        c%theta_steplen(2+m,i) = sigma
                        !write(*,*) i, 2+m, c%theta_steplen(2+m,i)
                      end if
                      c => c%nextComp()
                   end do

                 else if (comp_names(2)(1:4) == 'beta') then
                   c => compList
                   do while (associated(c))
                     if (trim(c%label) == trim(comp_names(1))) then
                        !       (beta+T+ntab, n_mcmc_samp_groups)
                        ! or    (beta+T,      n_mcmc_samp_groups)
                        c%theta_steplen(1,i) = sigma
                        if (cpar%myid == 0) write(*,*) 'set sigma_beta  to ',sigma, i
                      end if
                      c => c%nextComp()
                   end do
                 else if (comp_names(2)(1:1) == 'T') then
                   c => compList
                   do while (associated(c))
                     if (trim(c%label) == trim(comp_names(1))) then
                        !       (beta+T+ntab, n_mcmc_samp_groups)
                        ! or    (beta+T,      n_mcmc_samp_groups)
                        c%theta_steplen(2,i) = sigma
                        if (cpar%myid == 0) write(*,*) 'set sigma_T', sigma, i
                      end if
                      c => c%nextComp()
                   end do
                 else
                   if (cpar%myid == 0) write(*,*) 'Potentially poorly formatted ', trim(comp_tokens(1))
                 end if
             else
               write(*,*) 'Error in comp tokens', comp_tokens
               stop
             end if


           end do
        end do
    !end if

  end subroutine initialize_mh_mod

end module comm_mh_specind_mod
