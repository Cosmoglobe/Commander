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



  subroutine sample_gain_firas(outdir, cpar, handle, handle_noise, l)
    implicit none
    character(len=*),               intent(in)    :: outdir
    type(planck_rng),               intent(inout) :: handle, handle_noise
    integer(i4b),                   intent(in)    :: l
    type(comm_params) :: cpar


    real(dp)     :: chisq, my_chisq, chisq_old, chisq_new, chisq_prop, mval, mval_0
    real(dp), allocatable, dimension(:) :: sigmas, scales
    integer(i4b) :: band, ierr, i, j, k, m, pol, pix, n_scales, num_to_samp
    logical(lgt)  :: include_comp, reject, todo
    character(len=512) :: tokens(10), str_buff, operation
    class(comm_comp),   pointer           :: c => null()


    !do l = 1, cpar%mcmc_num_user_samp_groups
       ! Check if there are any gains to sample
    do m = 1, cpar%mcmc_samp_group_numstep(l)
       if (cpar%myid == 0) write(*,*) '|   Running MH sample ', m, ' out of ', cpar%mcmc_samp_group_numstep(l)
       num_to_samp = 0
       do i = 1, numband
         if (data(i)%gain_sigmas(l) > 0d0) then
           num_to_samp = num_to_samp + 1
         end if
       end do

       if (num_to_samp == 0) return
       if (cpar%myid == 0) then
         write(*,*) '| MH sampling group ', l
       end if

       call update_mixing_matrices(update_F_int=.true.)
       ! Calculate initial chisq
       chisq_old = 0d0
       call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_old, &
                          & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

       if (cpar%myid_chain .eq. 0) then
         write(*,*) '| Old chisq is ', chisq_old
       end if

       call store_buffer()

       ! Sample the gains
       do i = 1, numband
         data(i)%gain_tmp = data(i)%gain
         if (data(i)%gain_sigmas(l) > 0d0) then
           if (cpar%myid_chain == 0) then
             data(i)%gain = data(i)%gain + rand_gauss(handle)*data(i)%gain_sigmas(l)
             write(*,*) '| Gain ', trim(data(i)%label), ' sampled from ', &
                       & data(i)%gain_tmp, ' to ', data(i)%gain
           end if
           call mpi_bcast(data(i)%gain, 1, MPI_DOUBLE_PRECISION, 0, data(1)%info%comm, ierr)
         end if

       end do

       ! If there is any partner band, update that as well
       do i = 1, numband
         if (data(i)%gain_stat(l) > 0) then
           data(data(i)%gain_stat(l))%gain = data(i)%gain
           call mpi_bcast(data(data(i)%gain_stat(l))%gain, 1, MPI_DOUBLE_PRECISION, &
             & 0, data(1)%info%comm, ierr)
           if (cpar%myid_chain == 0) then
             write(*,*) '| Wired gain ', trim(data(data(i)%gain_stat(l))%label), & 
               & ' from ', data(data(i)%gain_stat(l))%gain_tmp, &
               & ' to ', data(data(i)%gain_stat(l))%gain
           end if
         end if
       end do

       ! Update mixing matrices
       call update_mixing_matrices(update_F_int=.true.)

       ! Perform component separation
       if (trim(cpar%mcmc_update_cg_groups(l)) == 'none') then
          if (cpar%myid == 0) write(*,*) '| No groups to sample'
       else
          if (cpar%myid == 0) write(*,*) '| Sampling CG groups ',trim(cpar%mcmc_update_cg_groups(l))
          call sample_all_amps_by_CG(cpar, handle, handle_noise, cg_groups=cpar%mcmc_update_cg_groups(l))
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

         ! Reset gains
         do i = 1, numband
           data(i)%gain = data(i)%gain_tmp
         end do

         ! Update mixing matrices
         call update_mixing_matrices(update_F_int=.true.)


         ! Instead of doing compsep, revert the amplitudes here
         if (trim(cpar%mcmc_update_cg_groups(l)) .ne. 'none') call revert_CG_amps(cpar)

         chisq_prop = 0d0
         call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_prop, &
                            & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

         if (cpar%myid_chain == 0) then
           write(*,*) '| '
           write(*,*) '| Chisq reverted back to ', chisq_prop, ' should be ', chisq_old
           write(*,*) '| '
         end if

       else
         if (cpar%myid_chain == 0) then
           write(*,*) '| '
           write(*,*) '| MH step accepted'
           write(*,*) '| '
         end if
       end if

  end do


    if (cpar%myid == 0) then
      write(*,*) '|   Finished sampling map gains'
    end if


  end subroutine sample_gain_firas

  subroutine sample_template_mh(outdir, cpar, handle, handle_noise, l)
    implicit none
    character(len=*),               intent(in)    :: outdir
    type(planck_rng),               intent(inout) :: handle, handle_noise
    integer(i4b),                   intent(in)    :: l
    type(comm_params) :: cpar


    real(dp)     :: chisq, my_chisq, chisq_old, chisq_new, chisq_prop, mval, mval_0
    real(dp), allocatable, dimension(:) :: scales
    integer(i4b) :: band, ierr, i, j, k, m, pol, pix, n_scales
    logical(lgt)  :: include_comp, reject, todo
    character(len=512) :: tokens(10), str_buff, operation
    class(comm_comp),   pointer           :: c => null()

    ! Given a component, propose an amplitude to shift the global amplitude of any component.


    !do l = 1, cpar%mcmc_num_user_samp_groups
    do m = 1, cpar%mcmc_samp_group_numstep(l)
       if (cpar%myid == 0) write(*,*) '|   Running MH sample ', m, ' out of ', cpar%mcmc_samp_group_numstep(l)

       n_scales = 0
       c => compList
       do while (associated(c))
          if (c%scale_sigma(l) > 0d0) then
            n_scales = n_scales + 1
          end if
          c => c%nextComp()
       end do

       if (n_scales == 0) return

       if (cpar%myid == 0) then
         write(*,*) '| MH sampling group ', l
       end if

       allocate(scales(n_scales))



       ! Calculate initial chisq
       chisq_old = 0d0
       call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_old, &
                          & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

       if (cpar%myid_chain .eq. 0) then
         write(*,*) '| Old chisq is ', chisq_old
       end if


       call store_buffer()



       ! Scale parameters

       i = 0
       c => compList
       do while (associated(c))
          if (c%scale_sigma(l) > 0d0) then
            i = i + 1
            if (cpar%myid == 0) then
              if (cpar%myid == 0) write(*,*) '|    ', trim(c%label)
              scales(i) = 1 + rand_gauss(handle)*c%scale_sigma(l)
              write(*,*) '|   Scaling by ', scales(i)
            end if
            call mpi_bcast(scales(i), 1, MPI_DOUBLE_PRECISION, 0, data(1)%info%comm, ierr)
            select type(c)
            class is (comm_diffuse_comp)
              c%x%alm = c%x%alm*scales(i)
              c%x_scale = c%x_scale * scales(i)
              !call c%x%Y
            class is (comm_template_comp)
              c%T%map = c%T%map*scales(i)
              c%T_scale = c%T_scale * scales(i)
            class default
              write(*,*) "You have not set behavior for class ", trim(c%class)
              stop
            end select
          end if
          c => c%nextComp()
       end do

       ! Perform component separation
       if (trim(cpar%mcmc_update_cg_groups(l)) == 'none') then
          if (cpar%myid == 0) write(*,*) '| No groups to sample'
       else
          if (cpar%myid == 0) write(*,*) '| Sampling CG groups ',trim(cpar%mcmc_update_cg_groups(l))
          call sample_all_amps_by_CG(cpar, handle, handle_noise, cg_groups=cpar%mcmc_update_cg_groups(l))
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

         i = 0
         c => compList
         do while (associated(c))
            if (c%scale_sigma(l) > 0d0) then
               i = i + 1
               select type(c)
               class is (comm_diffuse_comp)
                  c%x%alm = c%x%alm/scales(i)
                  !call c%x%Y
               class is (comm_template_comp)
                  c%T%map = c%T%map/scales(i)
               class default
                  write(*,*) "You have not set behavior for class ", trim(c%class)
                  stop
               end select
            end if
            c => c%nextComp()
         end do

         
         ! Instead of doing compsep, revert the amplitudes here
         if (trim(cpar%mcmc_update_cg_groups(l)) .ne. 'none') call revert_CG_amps(cpar)

         chisq_prop = 0d0
         call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_prop, &
                            & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

         if (cpar%myid_chain == 0) then
           write(*,*) '| '
           write(*,*) '| Chisq reverted back to ', chisq_prop, ' should be ', chisq_old
           write(*,*) '| '
         end if

       else
         if (cpar%myid_chain == 0) then
           write(*,*) '| '
           write(*,*) '| MH step accepted'
           write(*,*) '| '
         end if
       end if




       deallocate(scales)

    end do

    if (cpar%myid == 0) then
      write(*,*) '|   Finished sampling amplitude parameter'
    end if

  end subroutine sample_template_mh

  subroutine sample_mbbtab_mh(outdir, cpar, handle, handle_noise, l)
    implicit none
    character(len=*),               intent(in)    :: outdir
    type(planck_rng),               intent(inout) :: handle, handle_noise
    integer(i4b),                   intent(in)    :: l
    type(comm_params) :: cpar

    real(dp)     :: chisq, my_chisq, chisq_old, chisq_new, chisq_prop, mval, mval_0
    integer(i4b) :: band, ierr, i, j, k, m, pol, pix
    logical(lgt)  :: include_comp, reject, todo
    character(len=512) :: tokens(10), str_buff, operation
    class(comm_comp),   pointer           :: c => null()
    class(comm_map), pointer              :: invN_res => null(), map => null(), sig => null(), res => null()



    ! Loop over sampling groups

    !do l = 1, cpar%mcmc_num_user_samp_groups
    do m = 1, cpar%mcmc_samp_group_numstep(l)
       if (cpar%myid == 0) write(*,*) '|   Running MH sample ', m, ' out of ', cpar%mcmc_samp_group_numstep(l)

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
         return
       end if
       if (cpar%myid == 0) then
         write(*,*) '| MH sampling group ', l
       end if

       ! Calculate initial chisq
       chisq_old = 0d0
       call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_old, &
                          & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

       if (cpar%myid_chain .eq. 0) then
         write(*,*) '| Old chisq is ', chisq_old
       end if

       call store_buffer()

       c => compList
       do while (associated(c))
                       
          select type (c)
          type is (comm_MBBtab_comp)

            !if (c%id == k) then
            if (maxval(c%theta_steplen(c%npar+1:,l)) > 0) then

               if (cpar%myid_chain .eq. 0) then
                 write(*,*) '| ', trim(c%label)
                 c%SEDtab_buff = c%SEDtab
                 write(*,*) '| MBBtab original'
                 do i = 1, c%ntab
                   if (c%theta_steplen(c%npar+i,l) > 0) then
                     write(*,*) '|         ', c%SEDtab(3,i)
                   end if
                 end do
                 write(*,*) '| MBBtab proposal'
                 do i = 1, c%ntab
                   if (c%theta_steplen(c%npar+i,l) > 0) then
                     c%SEDtab(3,i) = c%SEDtab(3,i) + rand_gauss(handle) * c%theta_steplen(c%npar+i,l)
                     write(*,*) '|         ',  c%SEDtab(3,i)
                   end if
                 end do
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
       if (trim(cpar%mcmc_update_cg_groups(l)) == 'none') then
          if (cpar%myid == 0) write(*,*) '| No groups to sample'
       else
          if (cpar%myid == 0) write(*,*) '| Sampling CG groups ',trim(cpar%mcmc_update_cg_groups(l))
          call sample_all_amps_by_CG(cpar, handle, handle_noise, cg_groups=cpar%mcmc_update_cg_groups(l))
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


         c => compList
         do while (associated(c))
                         
            select type (c)
            type is (comm_MBBtab_comp)

              !if (c%id == k) then
              if (maxval(c%theta_steplen(c%npar+1:,l)) > 0) then
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
         if (trim(cpar%mcmc_update_cg_groups(l)) .ne. 'none') call revert_CG_amps(cpar)

         chisq_prop = 0d0
         call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_prop, &
                            & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

         if (cpar%myid_chain == 0) then
           write(*,*) '| '
           write(*,*) '| Chisq reverted back to ', chisq_prop, ' should be ', chisq_old
           write(*,*) '| '
         end if

       else
         if (cpar%myid_chain == 0) then
           write(*,*) '| '
           write(*,*) '| MH step accepted'
           write(*,*) '| '
         end if
       end if


     end do


     if (cpar%myid == 0) then
       write(*,*) '|   Finished sampling mbbtab'
     end if

  end subroutine sample_mbbtab_mh

  subroutine sample_specind_mh(outdir, cpar, handle, handle_noise, l)
    implicit none
    character(len=*),               intent(in)    :: outdir
    type(planck_rng),               intent(inout) :: handle, handle_noise
    integer(i4b),                   intent(in)    :: l
    type(comm_params) :: cpar

    real(dp)     :: chisq, my_chisq, chisq_old, chisq_new, chisq_prop, mval, mval_0
    integer(i4b) :: band, ierr, i, j, k, m, pol, pix
    logical(lgt)  :: include_comp, reject, todo, has_alms
    character(len=512) :: tokens(10), str_buff, operation
    class(comm_comp),   pointer           :: c => null()
    class(comm_map), pointer              :: invN_res => null(), map => null(), sig => null(), res => null()



    ! Loop over sampling groups

    !do l = 1, cpar%mcmc_num_user_samp_groups
    do m = 1, cpar%mcmc_samp_group_numstep(l)
       if (cpar%myid == 0) write(*,*) '|   Running MH sample ', m, ' out of ', cpar%mcmc_samp_group_numstep(l)

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
         return
       end if
       if (cpar%myid == 0) then
         write(*,*) '| MH sampling group ', l
       end if



       ! Calculate initial chisq
       chisq_old = 0d0
       call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_old, &
                          & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

       if (cpar%myid_chain .eq. 0) then
         write(*,*) '| Old chisq is ', chisq_old
       end if

       call store_buffer()


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

                do pol = 1,c%theta(j)%p%info%nmaps
                   if (c%lmax_ind_pol(pol,j) >= 0) cycle
                   do pix = 0, c%theta(j)%p%info%np-1
                      c%theta_pixreg(c%ind_pixreg_arr(pix,pol,j),pol,j) = c%theta(j)%p%map(pix,pol) 
                   end do
                end do

                if (any(c%lmax_ind_pol(:,j) >= 0)) call c%theta(j)%p%Y_scalar()
                do pol = 1,c%theta(j)%p%info%nmaps
                   if (c%lmax_ind_pol(pol,j) == -1) cycle
                  do pix = 0, c%theta(j)%p%info%np-1
                     c%theta_pixreg(c%ind_pixreg_arr(pix,pol,j),pol,j) = c%theta(j)%p%map(pix,pol) 
                  end do
               end do
               c%theta_pixreg_buff(:,:,j) = c%theta_pixreg(:,:,j)

               if (cpar%myid_chain == 0) then
                   do pol = 1, c%theta(j)%p%info%nmaps
                      c%theta_pixreg(1,pol,j) = c%theta_pixreg(1,pol,j) + rand_gauss(handle) * c%theta_steplen(j,l)
                   end do

                   write(*,*) '| ', trim(c%label), j
                   write(*,*) '| theta_pixreg original', c%theta_pixreg_buff(:,:,j)
                   write(*,*) '| theta_pixreg proposal', c%theta_pixreg(:,:,j)
               end if               
               call mpi_bcast(c%theta_pixreg(:,:,j),      size(c%theta_pixreg(:,:,j)),      MPI_DOUBLE_PRECISION, &
                 & 0, data(1)%info%comm, ierr)

               do pol = 1,c%theta(j)%p%info%nmaps
                  do pix = 0, c%theta(j)%p%info%np-1
                     c%theta(j)%p%map(pix,pol) = c%theta_pixreg(c%ind_pixreg_arr(pix,pol,j),pol,j)
                  end do
               end do
               if (any(c%lmax_ind_pol(:,j) >= 0)) call c%theta(j)%p%YtW_scalar()
 
             end select

          end do
          
          !go to next component
          c => c%nextComp()
              
       end do

       ! Update mixing matrices
       call update_mixing_matrices(update_F_int=.true.)

       ! Perform component separation
       if (trim(cpar%mcmc_update_cg_groups(l)) == 'none') then
          if (cpar%myid == 0) write(*,*) '| No groups to sample'
       else
          if (cpar%myid == 0) write(*,*) '| Sampling CG groups ',trim(cpar%mcmc_update_cg_groups(l))
          call sample_all_amps_by_CG(cpar, handle, handle_noise, cg_groups=cpar%mcmc_update_cg_groups(l))
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
            if (.not. allocated(c%theta_steplen)) then
               c => c%nextComp()
               cycle
            end if
                         
            do j = 1, c%npar
               if (c%theta_steplen(j,l) == 0.d0) cycle
               select type (c)
               class is (comm_diffuse_comp)
                  do pol = 1,c%theta(j)%p%info%nmaps
                     do pix = 0, c%theta(j)%p%info%np-1
                        c%theta(j)%p%map(pix,pol) = c%theta_pixreg_buff(c%ind_pixreg_arr(pix,pol,j),pol,j)
                     end do
                  end do
                  c%theta_pixreg(:,:,j) = c%theta_pixreg_buff(:,:,j)
                  if (any(c%lmax_ind_pol(:,j) >= 0)) call c%theta(j)%p%YtW_scalar()
               end select
            end do
            
            !go to next component
            c => c%nextComp()
                
         end do

         ! Update mixing matrices
         call update_mixing_matrices(update_F_int=.true.)

         if (trim(cpar%mcmc_update_cg_groups(l)) .ne. 'none') call revert_CG_amps(cpar)

         chisq_prop = 0d0
         call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_prop, &
                            & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

         if (cpar%myid_chain == 0) then
           write(*,*) '| '
           write(*,*) '| Chisq reverted back to ', chisq_prop, ' should be ', chisq_old
           write(*,*) '| '
         end if

       else
         if (cpar%myid_chain == 0) then
           write(*,*) '| '
           write(*,*) '| MH step accepted'
           write(*,*) '| '
         end if
       end if


     end do


     if (cpar%myid == 0) then
       write(*,*) '|   Finished sampling specind'
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

    implicit none
    type(comm_params) :: cpar


    type(comm_mapinfo), pointer :: info => null(), info_def => null(), info_ud
    class(comm_map), pointer :: indmask, mask_ud


    integer(i4b) :: i, j, k, nside, lmax, nmaps, n, m, n_tokens, n_in_group, ierr, num_gains
    integer(i4b) :: wire_to_ind, wire_from_ind
    integer(i4b), allocatable, dimension(:,:)  :: bands_to_sample
    real(dp) :: sigma
    logical(lgt) :: pol
    class(comm_comp),   pointer           :: c => null()


    character(len=128) :: tokens(100), comp_tokens(2), comp_names(2), comp_bands(2), wire_to(2), wire_from(2)
    character(len=2) :: itext

    allocate(cpar%mcmc_group_bands_indices(cpar%mcmc_num_user_samp_groups, numband))

    cpar%mcmc_group_bands_indices = 0


    ! Need to add an argument to sample_all_amps_by_CG that includes this list
    if (cpar%myid == 0) then
        do i = 1, cpar%mcmc_num_user_samp_groups
           !write(*,*) i, trim(cpar%mcmc_update_cg_groups(i)), trim(cpar%mcmc_samp_groups(i))
           if (trim(cpar%mcmc_update_cg_groups(i)) .ne. 'none') then
             call get_tokens(cpar%mcmc_update_cg_groups(i), ',', tokens, n)
             if (n == 0) then
               write(*,*) "Something is wrong with your CG groups"
               write(*,*) trim(cpar%mcmc_update_cg_groups(i))
             end if
           end if
        end do
    end if


    do i = 1, cpar%mcmc_num_user_samp_groups
      call get_tokens(cpar%mcmc_samp_group_bands(i), ',', tokens, n)
      if (n == 0) then
        write(*,*) 'You have misformatted MCMC_SAMPLING_GROUP_CHISQ_BANDS'
        stop
      else

      n_in_group = 0
      do j = 1, n
        k = findloc(cpar%ds_label, tokens(j), dim=1)
        do k = 1, numband
          if (trim(tokens(j)) == trim(data(k)%label)) then
            n_in_group = n_in_group + 1
            cpar%mcmc_group_bands_indices(i, n_in_group) = k
            exit
          end if
        end do
      end do

      if (maxval(cpar%mcmc_group_bands_indices(i,:)) == 0) then
        if (cpar%myid == 0) then
          call int2string(i,itext)
          write(*,*) 'MCMC_SAMPLING_GROUP_CHISQ_BANDS'//itext, ' does not reference any band that is loaded'
        end if
        stop
      end if

      end if
    end do



    ! Initializing, allocating all gain proposal lengths
    do i = 1, numband
      allocate(data(i)%gain_sigmas(cpar%mcmc_num_user_samp_groups))
      data(i)%gain_sigmas = 0d0
    end do

    do i = 1, cpar%mcmc_num_user_samp_groups


      call get_tokens(cpar%mcmc_samp_groups(i), ',', tokens, n_tokens)
      do j = 1, n_tokens
        if (trim(tokens(1)) == 'none') exit

        call get_tokens(tokens(j), '>', comp_tokens, n)
        if (n == 2) then
          call get_tokens(comp_tokens(1), ':', wire_from, num=n)
          call get_tokens(comp_tokens(2), ':', wire_to, num=n)
          wire_from_ind = 0
          wire_to_ind   = 0
          do k = 1, numband
            if (trim(data(k)%label) == trim(wire_from(2))) then
              wire_from_ind = k
            end if
            if (trim(data(k)%label) == trim(wire_to(2))) then
              wire_to_ind = k
            end if
          end do
          if (wire_from_ind == 0 .or. wire_to_ind == 0) then
            write(*,*) "You wiring points to a band that is not being used. ", trim(tokens(j))
            stop
          end if

          data(wire_to_ind)%gain_stat(i) = wire_from_ind


        else if (n == 1) then
            call get_tokens(tokens(j), '%', comp_tokens)
            read(comp_tokens(2), *) sigma

            call get_tokens(comp_tokens(1), ':', comp_names)


            if (trim(comp_names(1)) == 'gain') then
                do k = 1, numband
                  if (trim(comp_names(2)) .eq. trim(data(k)%label)) then
                    data(k)%gain_sigmas(i) = sigma
                  end if
                end do

            else if (trim(comp_names(2)) == 'scale') then

              c => compList
              do while (associated(c))
                 if (.not. allocated(c%scale_sigma)) then
                   allocate(c%scale_sigma(cpar%mcmc_num_user_samp_groups))
                   c%scale_sigma = 0d0
                 end if
                 if (trim(c%label) == trim(comp_names(1))) then
                   c%scale_sigma(i) = sigma
                 end if
                 c => c%nextComp()
              end do



            else if (comp_names(2)(1:3) == 'tab') then
              ! Get bin index
              call get_tokens(comp_names(2), '@', comp_bands)
              if (cpar%myid == 0 .and. .not. is_numeric(comp_bands(2))) then
                call int2string(i,itext)
                write(*,*) 'MCMC_SAMPLING_GROUP_PARAMS'//itext, '   ', trim(comp_bands(2)), ' should refer to an integer bin'
                write(*,*) trim(cpar%mcmc_samp_groups(i))
                stop
              else
                call mpi_barrier(cpar%comm_chain, ierr)
              end if
              read(comp_bands(2), *) m
             

              c => compList
              do while (associated(c))
                 if (trim(c%label) == trim(comp_names(1))) then
                   !       (beta+T+ntab, n_mcmc_samp_groups)
                   c%theta_steplen(2+m,i) = sigma
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

  end subroutine initialize_mh_mod


  subroutine store_buffer()
    implicit none
    integer(i4b) :: i, ind
    class(comm_comp), pointer :: c => null()

    ind = 1
    c   => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          c%x%alm_buff = c%x%alm
       class is (comm_ptsrc_comp)
          c%x_buff = c%x
       class is (comm_template_comp)
          c%x_buff = c%x
       end select
       c => c%nextComp()
    end do
    
  end subroutine store_buffer

  subroutine revert_CG_amps(cpar)
    implicit none

    type(comm_params), intent(in)    :: cpar

    integer(i4b) :: i, ind
    class(comm_comp), pointer :: c => null()

    if (cpar%myid == 0) then
      write(*,*) '|   Reverting to buffer values. Did you run sample_maps_with_CG with '
      write(*,*) '|   store_buff = .true.?'
    end if

    c   => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
         c%x%alm = c%x%alm_buff
       class is (comm_ptsrc_comp)
         c%x = c%x_buff
       class is (comm_template_comp)
         c%x = c%x_buff
       class default
         write(*,*) "What is this? no buffer exists", trim(c%type)
       end select
       c => c%nextComp()
    end do

  end subroutine revert_CG_amps

end module comm_mh_specind_mod
