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



end module comm_gain_mod
