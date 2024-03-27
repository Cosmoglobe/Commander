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
module comm_gain_mod
  use comm_signal_mod
  implicit none

contains
  
  subroutine sample_gain(operation, band, outdir, chain, iter, resamp_hard_prior, handle)
    implicit none
    integer(i4b),     intent(in)    :: band
    character(len=*), intent(in)    :: operation, outdir
    integer(i4b),     intent(in)    :: chain, iter
    logical(lgt),     intent(in)    :: resamp_hard_prior
    type(planck_rng), intent(inout) :: handle

    integer(i4b)  :: i, l, lmin, lmax, ierr, root, ncomp, dl_low, dl_high, ntok
    real(dp)      :: my_sigma, my_mu, mu, sigma, gain_new, chisq, mychisq
    real(dp)      :: MAX_DELTA_G = 0.3d0
    logical(lgt)  :: include_comp
    character(len=4) :: chain_text
    character(len=6) :: iter_text
    character(len=512) :: tokens(10)
    real(dp), allocatable, dimension(:,:) :: m, cls1, cls2
    class(comm_comp),   pointer           :: c => null()
    class(comm_map), pointer              :: invN_sig => null(), map => null(), sig => null(), res => null()


    ierr    = 0
    root    = 0
    dl_low  = 10
    dl_high = 25

    ! Handle bands with hard gain prior
    if (data(band)%gain_prior(2) < 0.d0) then
!!$       data(band)%gain = sum(data(6:9)%gain)/4.d0
!!$       return
       if (resamp_hard_prior) then
          if (data(band)%info%myid == root) then
             data(band)%gain = data(band)%gain_prior(1) + rand_gauss(handle) * abs(data(band)%gain_prior(2))
          end if
          call mpi_bcast(data(band)%gain, 1, MPI_DOUBLE_PRECISION, 0, data(band)%info%comm, ierr)  
       end if
       return
    end if

    ! Build reference signal
    sig => comm_map(data(band)%info)
    res => comm_map(data(band)%info)
    c => compList
    do while (associated(c))
       call get_tokens(trim(adjustl(data(band)%gain_comp)), ',', tokens, num=ntok)
       include_comp = .false.
       do i = 1, ntok
          if (trim(c%label) == trim(tokens(i)) .or. trim(tokens(i)) == 'all') then
             include_comp = .true.
             exit
          end if
       end do

       if (include_comp) then
          ! Add current component to calibration signal
          allocate(m(0:data(band)%info%np-1,data(band)%info%nmaps))
          m       = c%getBand(band)
          sig%map = sig%map + m
          deallocate(m)
       end if
       c => c%nextComp()
    end do

    ! Compute residual
    res                => compute_residual(band)
    data(band)%res%map =  res%map

    ! Add reference signal to residual
    res%map = data(band)%res%map + sig%map

    ! Divide by old gain
    sig%map = sig%map / data(band)%gain

    lmin = data(band)%gain_lmin
    lmax = data(band)%gain_lmax
    if (lmin > 0 .and. lmax > 0) then

       ! Apply mask
       if (associated(data(band)%gainmask)) then
          sig%map = sig%map * data(band)%gainmask%map
          res%map = res%map * data(band)%gainmask%map
       end if

       ! Compute cross-spectrum
       allocate(cls1(0:sig%info%lmax,sig%info%nspec))
       allocate(cls2(0:res%info%lmax,res%info%nspec))
       call sig%YtW
       call res%YtW

       call sig%getSigmaL(cls1)
       call sig%getCrossSigmaL(res, cls2)

       data(band)%gain = mean(cls2(lmin:lmax,1)/cls1(lmin:lmax,1))

       if (data(band)%info%myid == root) then
          call int2string(chain, chain_text)
          call int2string(iter,  iter_text)
          open(58,file=trim(outdir)//'/gain_cl_'//trim(data(band)%label)//'_c'//chain_text//'_k'//iter_text//'.dat', recl=1024)
          write(58,*) '#  ell     Ratio    C_l_cross     C_l_sig'
          do l = lmin, lmax
             write(58,*) l, cls2(l,1)/cls1(l,1), cls2(l,1), cls1(l,1)
          end do
          close(58)
       end if

       deallocate(cls1, cls2)

    else
       ! Correlate in pixel space with standard likelihood fit
       invN_sig     => comm_map(sig)
       call data(band)%N%invN(invN_sig)! Multiply with (invN)
       if (associated(data(band)%gainmask)) then
          invN_sig%map = invN_sig%map * data(band)%gainmask%map
          sig%map      = sig%map      * data(band)%gainmask%map
          res%map      = res%map      * data(band)%gainmask%map
       end if

       !call invN_sig%writeFITS('invN_sig_'//trim(data(band)%label)//'.fits')

       my_sigma = sum(sig%map * invN_sig%map)
       my_mu    = sum(res%map * invN_sig%map)
       call mpi_reduce(my_mu,    mu,    1, MPI_DOUBLE_PRECISION, MPI_SUM, root, data(band)%info%comm, ierr)
       call mpi_reduce(my_sigma, sigma, 1, MPI_DOUBLE_PRECISION, MPI_SUM, root, data(band)%info%comm, ierr)
       if (data(band)%info%myid == root) then
          ! Compute mu and sigma from likelihood term
          mu       = mu / sigma
          sigma    = sqrt(1.d0 / sigma)
          if (trim(operation) == 'optimize') then ! Optimize
             gain_new = mu
          else
             gain_new = mu + sigma * rand_gauss(handle)
          end if
          ! Only allow relatively small changes between steps, and not outside the range from 0.01 to 0.01
          data(band)%gain = min(max(gain_new, data(band)%gain-MAX_DELTA_G), data(band)%gain+MAX_DELTA_G)
          write(*,*) ' Posterior mean gain = ', mu, sigma, data(band)%gain
       end if

       ! Distribute new gains
       call mpi_bcast(data(band)%gain, 1, MPI_DOUBLE_PRECISION, 0, data(band)%info%comm, ierr)

       call invN_sig%dealloc(); deallocate(invN_sig)
    end if

    ! Subtract scaled reference signal to residual
    data(band)%res%map = res%map - data(band)%gain * sig%map

    ! Output residual signal and residual for debugging purposes
    if (.true.) then
       call sig%writeFITS(trim(outdir)//'/gain_sig_'//trim(data(band)%label)//'.fits')
       call res%writeFITS(trim(outdir)//'/gain_out'//trim(data(band)%label)//'.fits')
       call data(band)%res%writeFITS(trim(outdir)//'/gain_inp_'//trim(data(band)%label)//'.fits')
    end if

    call sig%dealloc(); deallocate(sig)
    call res%dealloc(); deallocate(res)

  end subroutine sample_gain


  subroutine sample_gain_firas(outdir, cpar, handle, handle_noise)
    implicit none
    character(len=*),               intent(in)    :: outdir
    type(planck_rng),               intent(inout) :: handle, handle_noise
    type(comm_params) :: cpar


    integer(i4b)  :: i, l, n_firas, n_sample, band, ntok, root, ierr, samp_group
    real(dp)      :: chisq, my_chisq, sigma, chisq_old, chisq_new, chisq_prop
    real(dp)      :: MAX_DELTA_G = 0.3d0
    logical(lgt)  :: include_comp, reject
    character(len=4) :: chain_text
    character(len=6) :: iter_text
    character(len=512) :: tokens(10), str_buff, operation
    integer(i4b), allocatable,  dimension(:) :: bands_sample, bands_firas
    real(dp), allocatable, dimension(:) :: gains_prop, gains_old, chisqs_old, chisqs_prop
    class(comm_comp),   pointer           :: c => null()
    class(comm_map), pointer              :: invN_res => null(), map => null(), sig => null(), res => null()



    root = 0



    ! Gain sampling for 545, 857, and 100-240 um, using FIRAS as the calibrator. 
    ! That is, we make a MH proposal to change the gain of all of these channels, 
    ! compute one step of the compsep amplitude sampling, and make an accept/reject decision based on the FIRAS chisq.

    !operation = cpar%operation
    !cpar%operation = 'optimize'


    n_sample = 0
    n_firas = 0
    do i = 1, numband
      ! Finds bands that we want to calibrate against FIRAS
      str_buff = data(i)%gain_comp
      call toupper(str_buff)
      if (index(str_buff, 'FIRAS') .ne. 0 .and. data(i)%sample_gain) then
        n_sample = n_sample + 1
      end if


      ! Identifies FIRAS bands
      str_buff = data(i)%label
      call toupper(str_buff)
      if (index(str_buff, 'FIRAS') .ne. 0) then
        n_firas = n_firas + 1
      end if

    end do

    if (n_sample .eq. 0) then
      return
    else if (n_firas .eq. 0 .and. n_sample > 0) then
      write(*,*) 'No FIRAS bands loaded, cannot calibrate as asked'
      stop
    else
      if (cpar%myid == root) then
         write(*,*) '|'
         write(*,*) '| MH sampling gain based on FIRAS'
         write(*,*) '|'
      end if
    end if
       
    allocate(bands_firas(n_firas), bands_sample(n_sample), gains_prop(n_sample), gains_old(n_sample))
    allocate(chisqs_old(n_firas), chisqs_prop(n_firas))

    n_sample = 0
    n_firas = 0
    do i = 1, numband
      str_buff = data(i)%gain_comp
      call toupper(str_buff)
      if (index(str_buff, 'FIRAS') .ne. 0 .and. (data(i)%sample_gain)) then
        n_sample = n_sample + 1
        bands_sample(n_sample) = i
      end if

      str_buff = data(i)%label
      call toupper(str_buff)
      if (index(str_buff, 'FIRAS') .ne. 0) then
        n_firas = n_firas + 1
        bands_firas(n_firas) = i
      end if

    end do

    chisq_old = 0d0

    do band = 1, n_firas
        ! Build reference signal
        sig => comm_map(data(bands_firas(band))%info)
        res => comm_map(data(bands_firas(band))%info)
        c => compList
        do while (associated(c))
           call get_tokens(trim(adjustl(data(bands_firas(band))%gain_comp)), ',', tokens, num=ntok)
           include_comp = .false.
           do i = 1, ntok
              if (trim(c%label) == trim(tokens(i)) .or. trim(tokens(i)) == 'all') then
                 include_comp = .true.
                 exit
              end if
           end do
           c => c%nextComp()
        end do

        ! Compute residual
        res                => compute_residual(bands_firas(band))
        data(bands_firas(band))%res%map =  res%map

        invN_res     => comm_map(res)
        call data(bands_firas(band))%N%invN(invN_res)! Multiply with (invN)

        if (associated(data(bands_firas(band))%gainmask)) then
           res%map      = res%map      * data(bands_firas(band))%gainmask%map
        end if

        my_chisq    = sum(res%map * invN_res%map)
        call mpi_reduce(my_chisq,    chisq,    1, MPI_DOUBLE_PRECISION, MPI_SUM, root, data(bands_firas(band))%info%comm, ierr)
        chisq_old = chisq_old + chisq
        chisqs_old(band) = chisq
    end do


    ! MH Step
    sigma = 0.005
    !sigma = 1e-6

    do i = 1, n_sample
      if (data(bands_sample(i))%info%myid == root) then
        gains_old(i) = data(bands_sample(i))%gain
        gains_prop(i) = gains_old(i) + rand_gauss(handle)*sigma
      end if

      call mpi_bcast(data(bands_sample(i))%gain, 1, MPI_DOUBLE_PRECISION, root, data(bands_sample(i))%info%comm, ierr)

    end do

    call mpi_bcast(gains_old, size(gains_old), MPI_DOUBLE_PRECISION, root, data(bands_sample(1))%info%comm, ierr)
    call mpi_bcast(gains_prop, size(gains_prop), MPI_DOUBLE_PRECISION, root, data(bands_sample(1))%info%comm, ierr)

    do i = 1, n_sample
       data(bands_sample(i))%gain = gains_prop(i)
    end do


    ! Update mixing matrices
    c => compList
    do while (associated(c))
       call c%updateMixmat
       c => c%nextComp()
    end do

    ! Do component separation

    if (cpar%myid_chain == 0) then
      do i = 1, n_sample
         write(*,*) '| ', trim(data(bands_sample(i))%label), gains_old(i), gains_prop(i)
      end do
      write(*,*) '|  Generating chisq for proposal gains'
    end if

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

    chisq_prop = 0d0

    if (data(bands_firas(1))%info%myid == root) then
       write(*,*) '|'
       do band = 1, n_firas
          write(*,*) '|  chisq old ', trim(data(bands_firas(band))%label), chisqs_old(band)
       end do
       write(*,*) '|'
    end if

    do band = 1, n_firas
        ! Build reference signal
        sig => comm_map(data(bands_firas(band))%info)
        res => comm_map(data(bands_firas(band))%info)
        c => compList
        do while (associated(c))
           call get_tokens(trim(adjustl(data(bands_firas(band))%gain_comp)), ',', tokens, num=ntok)
           include_comp = .false.
           do i = 1, ntok
              if (trim(c%label) == trim(tokens(i)) .or. trim(tokens(i)) == 'all') then
                 include_comp = .true.
                 exit
              end if
           end do
           c => c%nextComp()
        end do

        ! Compute residual
        res                => compute_residual(bands_firas(band))
        data(bands_firas(band))%res%map =  res%map

        invN_res     => comm_map(res)
        call data(bands_firas(band))%N%invN(invN_res)! Multiply with (invN)

        if (associated(data(bands_firas(band))%gainmask)) then
           res%map      = res%map      * data(bands_firas(band))%gainmask%map
        end if

        my_chisq    = sum(res%map * invN_res%map)
        call mpi_reduce(my_chisq,    chisq,    1, MPI_DOUBLE_PRECISION, MPI_SUM, root, data(band)%info%comm, ierr)
        chisq_prop = chisq_prop + chisq
        chisqs_prop(band) = chisq
        if (cpar%myid_chain == root) then
          write(*,*) '|  chisq prop ', trim(data(bands_firas(band))%label), chisqs_prop(band)
        end if
    end do
    if (cpar%myid_chain == root) then
       write(*,*) '|'
       do i = 1, n_sample
         write(*,*) '| ', trim(data(bands_sample(i))%label), ' old:', gains_old(i), ' prop:', gains_prop(i)
       end do
       write(*,*) '|'
       write(*,*) '| chisq diffs per band:'
       do i = 1, n_firas
         write(*,*) '|  ', trim(data(bands_firas(i))%label), chisqs_prop(i) - chisqs_old(i)
       end do
       write(*,*) '| chisq diff: ', chisq_prop-chisq_old
    end if


    reject = log(rand_uni(handle)) > (chisq_old - chisq_prop)/2
    call mpi_bcast(reject, 1, MPI_LOGICAL, root, data(bands_sample(1))%info%comm, ierr)


    if (reject) then
      if (cpar%myid_chain == 0) then
        write(*,*) '| '
        write(*,*) '| MH step rejected, returning to original gains.'
        write(*,*) '| '
      end if
      do i = 1, n_sample
        data(bands_sample(i))%gain = gains_old(i)
      end do

      ! Update mixing matrices
      c => compList
      do while (associated(c))
         call c%updateMixmat
         c => c%nextComp()
      end do

    else
      if (data(bands_firas(1))%info%myid == root) then
        write(*,*) '| '
        write(*,*) '| MH step accepted'
        write(*,*) '| New gains are'
        do i = 1, n_sample
          write(*,*) '| ', trim(data(bands_sample(i))%label), ':', gains_prop(i)
        end do
        write(*,*) '| '
      end if
    end if

    deallocate(bands_firas, bands_sample, gains_prop, gains_old)
    call invN_res%dealloc(); deallocate(invN_res)
    !cpar%operation = operation

  end subroutine sample_gain_firas

end module comm_gain_mod
