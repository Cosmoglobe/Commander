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



  subroutine sample_gain_firas(outdir, cpar, handle, handle_noise)
    implicit none
    character(len=*),               intent(in)    :: outdir
    type(planck_rng),               intent(inout) :: handle, handle_noise
    type(comm_params) :: cpar


    !integer(i4b)  :: i, l, n_firas, n_sample, band, ntok, root, ierr, samp_group
    !real(dp)      :: chisq, my_chisq, sigma, chisq_old, chisq_new, chisq_prop
    !real(dp)      :: MAX_DELTA_G = 0.3d0
    !logical(lgt)  :: include_comp, reject
    !character(len=4) :: chain_text
    !character(len=6) :: iter_text
    !character(len=512) :: tokens(10), str_buff, operation
    !integer(i4b), allocatable,  dimension(:) :: bands_sample, bands_firas
    !real(dp), allocatable, dimension(:) :: gains_prop, gains_old, chisqs_old, chisqs_prop
    !class(comm_comp),   pointer           :: c => null()
    !class(comm_map), pointer              :: invN_res => null(), map => null(), sig => null(), res => null()


    real(dp)     :: chisq, my_chisq, chisq_old, chisq_new, chisq_prop, mval, mval_0
    real(dp), allocatable, dimension(:) :: sigmas, scales
    integer(i4b) :: band, ierr, i, j, k, l, pol, pix, n_scales, num_to_samp
    logical(lgt)  :: include_comp, reject, todo
    character(len=512) :: tokens(10), str_buff, operation
    class(comm_comp),   pointer           :: c => null()


    do l = 1, cpar%mcmc_num_user_samp_groups
       ! Check if there are any gains to sample
       num_to_samp = 0
       do i = 1, numband
         if (data(i)%gain_sigmas(l) > 0d0) then
           num_to_samp = num_to_samp + 1
         end if
       end do

       if (num_to_samp == 0) cycle

       ! Calculate initial chisq
       chisq_old = 0d0
       call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_old, &
                          & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

       if (cpar%myid_chain .eq. 0) then
         write(*,*) '| Old chisq is ', chisq_old
       end if



       ! Sample the gains
       do i = 1, numband
         data(i)%gain_tmp = data(i)%gain
         if (cpar%myid == 0 .and. data(i)%gain_sigmas(l) > 0d0) then
           data(i)%gain = data(i)%gain + rand_gauss(handle)*data(i)%gain_sigmas(l)
           write(*,*) 'Gain sampled from ', data(i)%gain_tmp, data(i)%gain
         end if
         call mpi_bcast(data(i)%gain, 1, MPI_DOUBLE_PRECISION, 0, data(1)%info%comm, ierr)
         call mpi_bcast(data(i)%gain_tmp, 1, MPI_DOUBLE_PRECISION, 0, data(1)%info%comm, ierr)
       end do

       ! If there is any partner band, update that as well


       ! Update mixing matrices
       call update_mixing_matrices(update_F_int=.true.)


       ! Perform component separation
       if (trim(cpar%mcmc_update_cg_groups(l)) == 'none') then
          if (cpar%myid == 0) write(*,*) 'No groups to sample'
       else
          call sample_all_amps_by_CG(cpar, handle, handle_noise, store_buff=.true., cg_groups=cpar%mcmc_update_cg_groups(l))
       end if

       chisq_prop = 0d0
       call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_prop, &
                          & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

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
         call revert_CG_amps(cpar)

       else
         if (cpar%myid_chain == 0) then
           write(*,*) '| '
           write(*,*) '| MH step accepted'
           write(*,*) '| '
         end if
       end if

    !root = 0
  end do



    ! Gain sampling for 545, 857, and 100-240 um, using FIRAS as the calibrator. 
    ! That is, we make a MH proposal to change the gain of all of these channels, 
    ! compute one step of the compsep amplitude sampling, and make an accept/reject decision based on the FIRAS chisq.

    !operation = cpar%operation
    !cpar%operation = 'optimize'


    ! n_sample = 0
    ! n_firas = 0
    ! do i = 1, numband
    !   ! Finds bands that we want to calibrate against FIRAS
    !   str_buff = data(i)%gain_comp
    !   call toupper(str_buff)
    !   if (index(str_buff, 'FIRAS') .ne. 0 .and. data(i)%sample_gain) then
    !     n_sample = n_sample + 1
    !   end if


    !   ! Identifies FIRAS bands
    !   str_buff = data(i)%label
    !   call toupper(str_buff)
    !   if (index(str_buff, 'FIRAS') .ne. 0) then
    !     n_firas = n_firas + 1
    !   end if

    ! end do

    ! if (n_sample .eq. 0) then
    !   return
    ! else if (n_firas .eq. 0 .and. n_sample > 0) then
    !   write(*,*) 'No FIRAS bands loaded, cannot calibrate as asked'
    !   stop
    ! else
    !   if (cpar%myid == root) then
    !      write(*,*) '|'
    !      write(*,*) '| MH sampling gain based on FIRAS'
    !      write(*,*) '|'
    !   end if
    ! end if
    !    
    ! allocate(bands_firas(n_firas), bands_sample(n_sample), gains_prop(n_sample), gains_old(n_sample))
    ! allocate(chisqs_old(n_firas), chisqs_prop(n_firas))

    ! n_sample = 0
    ! n_firas = 0
    ! do i = 1, numband
    !   str_buff = data(i)%gain_comp
    !   call toupper(str_buff)
    !   if (index(str_buff, 'FIRAS') .ne. 0 .and. (data(i)%sample_gain)) then
    !     n_sample = n_sample + 1
    !     bands_sample(n_sample) = i
    !   end if

    !   str_buff = data(i)%label
    !   call toupper(str_buff)
    !   if (index(str_buff, 'FIRAS') .ne. 0) then
    !     n_firas = n_firas + 1
    !     bands_firas(n_firas) = i
    !   end if

    ! end do

    ! chisq_old = 0d0

    ! do band = 1, n_firas
    !     ! Build reference signal
    !     sig => comm_map(data(bands_firas(band))%info)
    !     res => comm_map(data(bands_firas(band))%info)
    !     c => compList
    !     do while (associated(c))
    !        call get_tokens(trim(adjustl(data(bands_firas(band))%gain_comp)), ',', tokens, num=ntok)
    !        include_comp = .false.
    !        do i = 1, ntok
    !           if (trim(c%label) == trim(tokens(i)) .or. trim(tokens(i)) == 'all') then
    !              include_comp = .true.
    !              exit
    !           end if
    !        end do
    !        c => c%nextComp()
    !     end do

    !     ! Compute residual
    !     res                => compute_residual(bands_firas(band))
    !     data(bands_firas(band))%res%map =  res%map

    !     invN_res     => comm_map(res)
    !     call data(bands_firas(band))%N%invN(invN_res)! Multiply with (invN)

    !     if (associated(data(bands_firas(band))%gainmask)) then
    !        res%map      = res%map      * data(bands_firas(band))%gainmask%map
    !     end if

    !     my_chisq    = sum(res%map * invN_res%map)
    !     call mpi_reduce(my_chisq,    chisq,    1, MPI_DOUBLE_PRECISION, MPI_SUM, root, data(bands_firas(band))%info%comm, ierr)
    !     chisq_old = chisq_old + chisq
    !     chisqs_old(band) = chisq
    ! end do


    ! ! MH Step
    ! sigma = 0.005
    ! !sigma = 1e-6

    ! do i = 1, n_sample
    !   if (data(bands_sample(i))%info%myid == root) then
    !     gains_old(i) = data(bands_sample(i))%gain
    !     gains_prop(i) = gains_old(i) + rand_gauss(handle)*sigma
    !   end if

    !   call mpi_bcast(data(bands_sample(i))%gain, 1, MPI_DOUBLE_PRECISION, root, data(bands_sample(i))%info%comm, ierr)

    ! end do

    ! call mpi_bcast(gains_old, size(gains_old), MPI_DOUBLE_PRECISION, root, data(bands_sample(1))%info%comm, ierr)
    ! call mpi_bcast(gains_prop, size(gains_prop), MPI_DOUBLE_PRECISION, root, data(bands_sample(1))%info%comm, ierr)

    ! do i = 1, n_sample
    !    data(bands_sample(i))%gain = gains_prop(i)
    ! end do


    ! ! Update mixing matrices
    ! c => compList
    ! do while (associated(c))
    !    call c%updateMixmat
    !    c => c%nextComp()
    ! end do


    ! if (cpar%myid_chain == 0) then
    !   do i = 1, n_sample
    !      write(*,*) '| ', trim(data(bands_sample(i))%label), gains_old(i), gains_prop(i)
    !   end do
    !   write(*,*) '|  Generating chisq for proposal gains'
    ! end if

    ! ! Do component separation
    ! call sample_all_amps_by_CG(cpar, handle, handle_noise)

    ! chisq_prop = 0d0

    ! if (data(bands_firas(1))%info%myid == root) then
    !    write(*,*) '|'
    !    do band = 1, n_firas
    !       write(*,*) '|  chisq old ', trim(data(bands_firas(band))%label), chisqs_old(band)
    !    end do
    !    write(*,*) '|'
    ! end if

    ! do band = 1, n_firas
    !     ! Build reference signal
    !     sig => comm_map(data(bands_firas(band))%info)
    !     res => comm_map(data(bands_firas(band))%info)
    !     c => compList
    !     do while (associated(c))
    !        call get_tokens(trim(adjustl(data(bands_firas(band))%gain_comp)), ',', tokens, num=ntok)
    !        include_comp = .false.
    !        do i = 1, ntok
    !           if (trim(c%label) == trim(tokens(i)) .or. trim(tokens(i)) == 'all') then
    !              include_comp = .true.
    !              exit
    !           end if
    !        end do
    !        c => c%nextComp()
    !     end do

    !     ! Compute residual
    !     res                => compute_residual(bands_firas(band))
    !     data(bands_firas(band))%res%map =  res%map

    !     invN_res     => comm_map(res)
    !     call data(bands_firas(band))%N%invN(invN_res)! Multiply with (invN)

    !     if (associated(data(bands_firas(band))%gainmask)) then
    !        res%map      = res%map      * data(bands_firas(band))%gainmask%map
    !     end if

    !     my_chisq    = sum(res%map * invN_res%map)
    !     call mpi_reduce(my_chisq,    chisq,    1, MPI_DOUBLE_PRECISION, MPI_SUM, root, data(band)%info%comm, ierr)
    !     chisq_prop = chisq_prop + chisq
    !     chisqs_prop(band) = chisq
    !     if (cpar%myid_chain == root) then
    !       write(*,*) '|  chisq prop ', trim(data(bands_firas(band))%label), chisqs_prop(band)
    !     end if
    ! end do
    ! if (cpar%myid_chain == root) then
    !    write(*,*) '|'
    !    do i = 1, n_sample
    !      write(*,*) '| ', trim(data(bands_sample(i))%label), ' old:', gains_old(i), ' prop:', gains_prop(i)
    !    end do
    !    write(*,*) '|'
    !    write(*,*) '| chisq diffs per band:'
    !    do i = 1, n_firas
    !      write(*,*) '|  ', trim(data(bands_firas(i))%label), chisqs_prop(i) - chisqs_old(i)
    !    end do
    !    write(*,*) '| chisq diff: ', chisq_prop-chisq_old
    ! end if


    ! reject = log(rand_uni(handle)) > (chisq_old - chisq_prop)/2
    ! call mpi_bcast(reject, 1, MPI_LOGICAL, root, data(bands_sample(1))%info%comm, ierr)


    ! if (reject) then
    !   if (cpar%myid_chain == 0) then
    !     write(*,*) '| '
    !     write(*,*) '| MH step rejected, returning to original gains.'
    !     write(*,*) '| '
    !   end if
    !   do i = 1, n_sample
    !     data(bands_sample(i))%gain = gains_old(i)
    !   end do

    !   ! Update mixing matrices
    !   c => compList
    !   do while (associated(c))
    !      call c%updateMixmat
    !      c => c%nextComp()
    !   end do

    ! else
    !   if (data(bands_firas(1))%info%myid == root) then
    !     write(*,*) '| '
    !     write(*,*) '| MH step accepted'
    !     write(*,*) '| New gains are'
    !     do i = 1, n_sample
    !       write(*,*) '| ', trim(data(bands_sample(i))%label), ':', gains_prop(i)
    !     end do
    !     write(*,*) '| '
    !   end if
    ! end if

    ! deallocate(bands_firas, bands_sample, gains_prop, gains_old)
    ! call invN_res%dealloc(); deallocate(invN_res)
    !cpar%operation = operation

  end subroutine sample_gain_firas

  subroutine sample_template_mh(outdir, cpar, handle, handle_noise)
    implicit none
    character(len=*),               intent(in)    :: outdir
    type(planck_rng),               intent(inout) :: handle, handle_noise
    type(comm_params) :: cpar


    real(dp)     :: chisq, my_chisq, chisq_old, chisq_new, chisq_prop, mval, mval_0
    real(dp), allocatable, dimension(:) :: sigmas, scales
    integer(i4b) :: band, ierr, i, j, k, l, pol, pix, n_scales
    logical(lgt)  :: include_comp, reject, todo
    character(len=512) :: tokens(10), str_buff, operation
    class(comm_comp),   pointer           :: c => null()

    ! Given a component, propose an amplitude to shift the global amplitude of any component.


    do l = 1, cpar%mcmc_num_user_samp_groups



       n_scales = 0
       c => compList
       do while (associated(c))
          if (c%scale_sigma(l) > 0d0) then
            n_scales = n_scales + 1
          end if
          c => c%nextComp()
       end do

       if (n_scales == 0) cycle

       allocate(sigmas(n_scales), scales(n_scales))



       ! Calculate initial chisq
       chisq_old = 0d0
       call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_old, &
                          & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

       if (cpar%myid_chain .eq. 0) then
         write(*,*) '| Old chisq is ', chisq_old
       end if





       ! Scale parameters
       i = 0
       c => compList
       do while (associated(c))
          if (c%scale_sigma(l) > 0d0) then
            i = i + 1
            sigmas(i) = c%scale_sigma(l)
            if (cpar%myid == 0) then
              scales(i) = 1 + rand_gauss(handle)*sigmas(i)
              write(*,*) 'Scaling by ', scales(i)
            end if
            call mpi_bcast(scales(i), 1, MPI_DOUBLE_PRECISION, 0, data(1)%info%comm, ierr)
            select type(c)
            class is (comm_diffuse_comp)
              c%x%map = c%x%map*scales(i)
              call c%x%YtW
            class is (comm_template_comp)
              c%T%map = c%T%map*scales(i)
            class default
              write(*,*) "You have not set behavior for class ", trim(c%class)
              stop
            end select
          end if
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

       chisq_prop = 0d0
       call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_prop, &
                          & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

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
                c%x%map = c%x%map/scales(i)
                call c%x%YtW
              class is (comm_template_comp)
                c%T%map = c%T%map/scales(i)
              class default
                write(*,*) "You have not set behavior for class ", trim(c%class)
                stop
              end select
            end if
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




       deallocate(sigmas, scales)

    end do

    if (cpar%myid == 0) then
      write(*,*) 'Amplitude parameter'
    end if

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

       chisq_prop = 0d0
       call compute_chisq(data(1)%info%comm, chisq_fullsky=chisq_prop, &
                          & maskpath=cpar%mcmc_samp_group_mask(l), band_list=cpar%mcmc_group_bands_indices(l,:))

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
         write(*,*) "|    Proposal chisq is ", nint(chisq_prop, i8b)
         write(*,*) "|    Delta chi^2 is    ", nint(chisq_prop - chisq_old, i8b)
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


    integer(i4b) :: i, j, k, nside, lmax, nmaps, n, m, n_tokens, n_in_group, ierr, num_gains
    integer(i4b), allocatable, dimension(:,:)  :: bands_to_sample
    real(dp) :: sigma
    logical(lgt) :: pol
    class(comm_comp),   pointer           :: c => null()


    character(len=128) :: tokens(100), comp_tokens(2), comp_names(2), comp_bands(2), wire_to(2), wire_from(2)

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

      end if
    end do


    ! Initializing, allocating all gain proposal lengths
    do i = 1, numband
      allocate(data(i)%gain_sigmas(cpar%mcmc_num_user_samp_groups))
      data(i)%gain_sigmas = 0d0
    end do

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
              if (cpar%myid == 0) then
                write(*,*) 'We are sampling gain', sigma
                write(*,*) trim(comp_names(2)), 'band gain to be sampled'
                do k = 1, numband
                  if (trim(comp_names(2)) .eq. trim(data(k)%label)) then
                    data(k)%gain_sigmas(i) = sigma
                  end if
                end do
              end if

            else if (trim(comp_names(2)) == 'scale') then
              if (cpar%myid == 0) write(*,*) 'We are scaling the amplitude'

              c => compList
              do while (associated(c))
                 if (.not. allocated(c%scale_sigma)) then
                   allocate(c%scale_sigma(cpar%mcmc_num_user_samp_groups))
                   c%scale_sigma = 0d0
                 end if
                 if (trim(c%label) == trim(comp_names(1))) then
                   c%scale_sigma(i) = sigma
                 else
                   c%scale_sigma(i) = 0d0
                 end if
                 c => c%nextComp()
              end do



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

  end subroutine initialize_mh_mod

end module comm_mh_specind_mod
