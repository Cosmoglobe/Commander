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

  subroutine sample_mbbtab_mh_sample(outdir, cpar, handle, handle_noise)
    implicit none
    character(len=*),               intent(in)    :: outdir
    type(planck_rng),               intent(inout) :: handle, handle_noise
    type(comm_params) :: cpar

    real(dp)     :: chisq, my_chisq, chisq_old, chisq_new, chisq_prop
    integer(i4b) :: band, samp_group, ierr, i, j, k, pol, pix
    logical(lgt)  :: include_comp, reject, todo
    character(len=512) :: tokens(10), str_buff, operation
    class(comm_comp),   pointer           :: c => null()
    class(comm_map), pointer              :: invN_res => null(), map => null(), sig => null(), res => null()



    ! Calculate initial chisq
    todo = .true.
    c => compList
    chisq_old = 0d0
    do while (associated(c))
       select type (c)
       type is (comm_md_comp)
         continue
       class is (comm_diffuse_comp)
         if (todo) then
           if (allocated(c%indmask)) then
             call compute_chisq(c%comm, chisq_fullsky=chisq_old, mask=c%indmask)
           end if
           todo = .false.
         end if
       end select
       c => c%nextComp()
    end do

    if (cpar%myid_chain .eq. 0) then
      write(*,*) '| Old chisq is ', chisq_old
    end if

    c => compList
    do while (associated(c))
                    
       select type (c)
       type is (comm_MBBtab_comp)

         if (cpar%myid_chain .eq. 0) then
           write(*,*) trim(c%label)
           c%SEDtab_buff = c%SEDtab
           do i = 1, c%ntab
              c%SEDtab(3,i) = c%SEDtab(3,i) + rand_gauss(handle) * 0.1
           end do
           write(*,*) 'MBBtab original', c%SEDtab_buff(3,:)
           write(*,*) 'MBBtab proposal', c%SEDtab(3,:)
         end if


         call mpi_bcast(c%SEDtab, size(c%SEDtab), MPI_DOUBLE_PRECISION, &
           & 0, data(1)%info%comm, ierr)
         call mpi_bcast(c%SEDtab_buff, size(c%SEDtab), MPI_DOUBLE_PRECISION, &
           & 0, data(1)%info%comm, ierr)

       end select
       
       !go to next component
       c => c%nextComp()
    end do

    call update_mixing_matrices(update_F_int=.true.)

    ! Perform component separation
    call timer%start(TOT_AMPSAMP)
    do samp_group = 1, cpar%cg_num_user_samp_groups
       if (cpar%myid_chain == 0) then
          write(*,fmt='(a,i4,a,i4,a,i4)') ' |  Chain = ', cpar%mychain, ' -- CG sample group = ', &
               & samp_group, ' of ', cpar%cg_num_user_samp_groups
       end if
       call sample_amps_by_CG(cpar, samp_group, handle, handle_noise)

       if (trim(cpar%cmb_dipole_prior_mask) /= 'none') call apply_cmb_dipole_prior(cpar, handle)

    end do
    call timer%stop(TOT_AMPSAMP)


    todo = .true.
    chisq_prop = 0d0
    c => compList
    do while (associated(c))
       select type (c)
       type is (comm_md_comp)
         continue
       class is (comm_diffuse_comp)
         if (todo) then
           if (allocated(c%indmask)) then
             call compute_chisq(c%comm, chisq_fullsky=chisq_prop, mask=c%indmask)
           end if
           todo = .false.
         end if
       end select
       c => c%nextComp()
    end do

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
        write(*,*) '| Rejected chisq is ', chisq_old
        write(*,*) '| '
      end if


      c => compList
      do while (associated(c))
                      
         select type (c)
         type is (comm_MBBtab_comp)

           if (cpar%myid_chain .eq. 0) then
             c%SEDtab_buff = c%SEDtab
             do i = 1, c%ntab
                c%SEDtab(3,i) = c%SEDtab_buff(3,i)
             end do
           end if


           call mpi_bcast(c%SEDtab, size(c%SEDtab), MPI_DOUBLE_PRECISION, &
             & 0, data(1)%info%comm, ierr)
           call mpi_bcast(c%SEDtab_buff, size(c%SEDtab), MPI_DOUBLE_PRECISION, &
             & 0, data(1)%info%comm, ierr)

         end select
         
         !go to next component
         c => c%nextComp()
      end do

      ! Update mixing matrices
      call update_mixing_matrices(update_F_int=.true.)

      ! Perform component separation
      call timer%start(TOT_AMPSAMP)
      do samp_group = 1, cpar%cg_num_user_samp_groups
         if (cpar%myid_chain == 0) then
            write(*,fmt='(a,i4,a,i4,a,i4)') ' |  Chain = ', cpar%mychain, ' -- CG sample group = ', &
                 & samp_group, ' of ', cpar%cg_num_user_samp_groups
         end if
         call sample_amps_by_CG(cpar, samp_group, handle, handle_noise)

         if (trim(cpar%cmb_dipole_prior_mask) /= 'none') call apply_cmb_dipole_prior(cpar, handle)

      end do
      call timer%stop(TOT_AMPSAMP)

      todo = .true.
      c => compList
      chisq_old = 0d0
      do while (associated(c))
         select type (c)
         type is (comm_md_comp)
           continue
         class is (comm_diffuse_comp)
           if (todo) then
             if (allocated(c%indmask)) then
               call compute_chisq(c%comm, chisq_fullsky=chisq_old, mask=c%indmask)
             end if
             todo = .false.
           end if
         end select
         c => c%nextComp()
      end do

      if (cpar%myid_chain == 0) then
        write(*,*) '| '
        write(*,*) '| Current MH chisq is ', chisq_old
        write(*,*) '| '
      end if

    else
      if (cpar%myid_chain == 0) then
        write(*,*) '| '
        write(*,*) '| MH step accepted'
        write(*,*) '| Current MH chisq is ', chisq_prop
      end if
    end if

  end subroutine sample_mbbtab_mh_sample

  subroutine sample_specind_mh_sample(outdir, cpar, handle, handle_noise)
    implicit none
    character(len=*),               intent(in)    :: outdir
    type(planck_rng),               intent(inout) :: handle, handle_noise
    type(comm_params) :: cpar

    real(dp)     :: chisq, my_chisq, chisq_old, chisq_new, chisq_prop
    integer(i4b) :: band, samp_group, ierr, i, j, k, pol, pix
    logical(lgt)  :: include_comp, reject, todo, has_alms
    character(len=512) :: tokens(10), str_buff, operation
    class(comm_comp),   pointer           :: c => null()
    class(comm_map), pointer              :: invN_res => null(), map => null(), sig => null(), res => null()



    ! Calculate initial chisq
    todo = .true.
    c => compList
    chisq_old = 0d0
    do while (associated(c))
       select type (c)
       type is (comm_md_comp)
         continue
       class is (comm_diffuse_comp)
         if (todo) then
           if (allocated(c%indmask)) then
             ! This is a hack, and should be replaced with a real mask call.
             call compute_chisq(c%comm, chisq_fullsky=chisq_old, mask=c%indmask)
             todo = .false.
           end if
         end if
       end select
       c => c%nextComp()
    end do

    if (todo) then
      write(*,*) 'Big oops'
      stop
    end if

    if (cpar%myid_chain .eq. 0) then
      write(*,*) '| Old chisq is ', nint(chisq_old, i8b)
    end if

    ! Propose spectral index MH step - for a given component or all components?

    ! From the sample_nonlin_params section
    c => compList
    do while (associated(c))
       if (c%npar == 0) then
          c => c%nextComp()
          cycle
       end if
                    
       do j = 1, c%npar
          if (c%p_gauss(2,j) == 0.d0) then
            cycle
          end if
          select type (c)
          class is (comm_diffuse_comp)


            if (cpar%myid_chain .eq. 0) then
              write(*,*) trim(c%label), j
              do pol = 1, c%theta(j)%p%info%nmaps
                 c%theta_pixreg_buff(1,pol,j) = c%theta_pixreg(1,pol,j)
                 c%theta_pixreg(1,pol,j) = c%theta_pixreg(1,pol,j) + rand_gauss(handle) * c%p_gauss(2,j)
              end do
              write(*,*) 'theta_pixreg original', c%theta_pixreg_buff(1,:,j)
              write(*,*) 'theta_pixreg proposal', c%theta_pixreg(1,:,j)
            end if


            call mpi_bcast(c%theta_pixreg_buff, size(c%theta_pixreg_buff), MPI_DOUBLE_PRECISION, &
              & 0, data(1)%info%comm, ierr)
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


          end select

       end do
       
       !go to next component
       c => c%nextComp()
           
    end do

    ! Update mixing matrices
    c => compList
    do while (associated(c))
       call c%updateMixmat
       c => c%nextComp()
    end do

    ! Perform component separation
    call timer%start(TOT_AMPSAMP)
    do samp_group = 1, cpar%cg_num_user_samp_groups
       if (cpar%myid_chain == 0) then
          write(*,fmt='(a,i4,a,i4,a,i4)') ' |  Chain = ', cpar%mychain, ' -- CG sample group = ', &
               & samp_group, ' of ', cpar%cg_num_user_samp_groups
       end if
       call sample_amps_by_CG(cpar, samp_group, handle, handle_noise)

       if (trim(cpar%cmb_dipole_prior_mask) /= 'none') call apply_cmb_dipole_prior(cpar, handle)

    end do
    call timer%stop(TOT_AMPSAMP)


    todo = .true.
    chisq_prop = 0d0
    c => compList
    do while (associated(c))
       select type (c)
       type is (comm_md_comp)
         continue
       class is (comm_diffuse_comp)
         if (todo) then
           if (allocated(c%indmask)) then
             call compute_chisq(c%comm, chisq_fullsky=chisq_prop, mask=c%indmask)
             todo = .false.
           end if
         end if
       end select
       c => c%nextComp()
    end do

    if (cpar%myid_chain .eq. 0) then
      write(*,*) "|    Proposal chisq is ", nint(chisq_prop, i8b)
      write(*,*) "|    Delta chi^2 is    ", nint(chisq_prop - chisq_old, i8b)
    end if

    ! Check MH statistic
    reject = log(rand_uni(handle)) > (chisq_old - chisq_prop)/2
    call mpi_bcast(reject, 1, MPI_LOGICAL, 0, data(1)%info%comm, ierr)


    if (reject) then
      if (cpar%myid_chain == 0) then
        write(*,*) '| '
        write(*,*) '| MH step rejected, returning to original spectral indices.'
        write(*,*) '| Rejected chisq is ', nint(chisq_old, i8b)
        write(*,*) '| '
      end if

      c => compList
      do while (associated(c))
         if (c%npar == 0) then
            c => c%nextComp()
            cycle
         end if
                      
         do j = 1, c%npar
            if (c%p_gauss(2,j) == 0.d0) then
              cycle
            end if
            select type (c)
            class is (comm_diffuse_comp)

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
            end select

         end do
         
         !go to next component
         c => c%nextComp()
             
      end do

      ! Update mixing matrices
      c => compList
      do while (associated(c))
         call c%updateMixmat
         c => c%nextComp()
      end do

      ! Perform component separation
      call timer%start(TOT_AMPSAMP)
      do samp_group = 1, cpar%cg_num_user_samp_groups
         if (cpar%myid_chain == 0) then
            write(*,fmt='(a,i4,a,i4,a,i4)') ' |  Chain = ', cpar%mychain, ' -- CG sample group = ', &
                 & samp_group, ' of ', cpar%cg_num_user_samp_groups
         end if
         call sample_amps_by_CG(cpar, samp_group, handle, handle_noise)

         if (trim(cpar%cmb_dipole_prior_mask) /= 'none') call apply_cmb_dipole_prior(cpar, handle)

      end do
      call timer%stop(TOT_AMPSAMP)

      todo = .true.
      c => compList
      chisq_old = 0d0
      do while (associated(c))
         select type (c)
         type is (comm_md_comp)
           continue
         class is (comm_diffuse_comp)
           if (todo) then
             if (allocated(c%indmask)) then
               call compute_chisq(c%comm, chisq_fullsky=chisq_old, mask=c%indmask)
               todo = .false.
             end if
           end if
         end select
         c => c%nextComp()
      end do

      if (cpar%myid_chain == 0) then
        write(*,*) '| '
        write(*,*) '| Current MH chisq is ', nint(chisq_old, i8b)
        write(*,*) '| '
      end if

    else
      if (cpar%myid_chain == 0) then
        write(*,*) '| '
        write(*,*) '| MH step accepted'
        write(*,*) '| Current MH chisq is ', nint(chisq_prop, i8b)
      end if
    end if

  end subroutine sample_specind_mh_sample

end module comm_mh_specind_mod
