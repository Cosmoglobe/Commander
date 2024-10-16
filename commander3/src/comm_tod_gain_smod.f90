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
submodule (comm_tod_gain_mod) comm_tod_gain_smod
  use omp_lib
contains

  module subroutine calculate_gain_mean_std_per_scan(tod, scan_id, s_invsqrtN, mask, s_tot, handle, mask_lowres, tod_arr)
      ! 
      ! Calculate the scan-specific contributions to the gain estimation. 
      !
      ! This subroutine is called during the processing of each scan, and the
      ! results are stored in the TOD object. Then, after each scan has been
      ! processed in the main loop, these data are used to calculate the total
      ! gain estimate for each detector, since these gain estimates depend on
      ! knowing all the scan-specific gain contributions. The data will be
      ! downsampled before calculating the estimate.
      !
      ! Arguments:
      ! ----------
      ! tod:           derived class (comm_tod)
      !                TOD object containing all scan-relevant data. Will be
      !                modified.
      ! scan_id:       integer(i4b)
      !                The ID of the current scan.
      ! s_invsqrtN:    real(sp) array
      !                The product of the reference signal we want to calibrate
      !                on multiplied by the square root of the inverse noise
      !                matrix for this scan.
      ! handle:        derived class (planck_rng)
      !                Random number generator handle. Will be modified.
      ! mask_lowres:   real(sp) array
      !                The mask for the low-resolution downsampled data.
      !
      ! tod_arr:       To be deprecated soon.
      !
      ! Returns:
      ! --------
      ! tod:           derived class (comm_tod)
      !                Will update the fields tod%scans(scan_id)%d(:)%dgain and
      !                tod%scans(scan_id)%d(:)%gain_invsigma, which contain
      !                incremental estimates to be used for global gain
      !                estimation. 
      !                dgain = s_ref * n_invN * residual where residual is
      !                     d - (g_0 + g_i) * s_tot (g_0 being the absolute gain, and
      !                     g_i the absolute gain per detector)
      !                gain_invsigma =  s_ref * n_invN * s_ref

      
    implicit none
    class(comm_tod),                      intent(inout) :: tod
    real(sp),             dimension(:,:), intent(in)    :: tod_arr
    real(sp),             dimension(:,:), intent(in)    :: s_invsqrtN, mask, s_tot
    integer(i4b),                         intent(in)    :: scan_id
    type(planck_rng),                     intent(inout) :: handle
    real(sp),             dimension(:,:), intent(in), optional :: mask_lowres


    real(sp), allocatable, dimension(:,:) :: residual
    real(sp), allocatable, dimension(:)   :: r_fill
    integer(i4b) :: ext(2), j, i, ndet, ntod
    character(len=5) :: itext

    ndet = tod%ndet
    ntod = size(s_tot,1)

    allocate(r_fill(size(s_tot,1)))
    call tod%downsample_tod(s_tot(:,1), ext)    
    allocate(residual(ext(1):ext(2),ndet))
    do j = 1, ndet
       if (.not. tod%scans(scan_id)%d(j)%accept) then
          residual(:,j) = 0.d0
          cycle
       end if
       r_fill = tod_arr(:, j) - (tod%gain0(0) + tod%gain0(j)) * s_tot(:,j)
       call fill_all_masked(r_fill, mask(:,j), ntod, trim(tod%operation) == 'sample', real(tod%scans(scan_id)%d(j)%N_psd%sigma0, sp), handle, tod%scans(scan_id)%chunk_num)
       call tod%downsample_tod(r_fill, ext, residual(:,j))
    end do
    call multiply_inv_N(tod, scan_id, residual, sampfreq=tod%samprate_lowres, pow=0.5d0)

    do j = 1, ndet
       if (.not. tod%scans(scan_id)%d(j)%accept) then
          tod%scans(scan_id)%d(j)%gain  = 0.d0
          tod%scans(scan_id)%d(j)%dgain = 0.d0
       else
          if (present(mask_lowres)) then
             tod%scans(scan_id)%d(j)%dgain         = sum(s_invsqrtN(:,j) * residual(:,j) * mask_lowres(:,j))
             tod%scans(scan_id)%d(j)%gain_invsigma = sum(s_invsqrtN(:,j) ** 2  * mask_lowres(:,j))
          else
             tod%scans(scan_id)%d(j)%dgain         = sum(s_invsqrtN(:,j) * residual(:,j))
             tod%scans(scan_id)%d(j)%gain_invsigma = sum(s_invsqrtN(:,j) ** 2)
          end if
          if (tod%scans(scan_id)%d(j)%gain_invsigma < 0.d0) then
             write(*,*) 'Warning: Not positive definite invN = ', tod%scanid(scan_id), j, tod%scans(scan_id)%d(j)%gain_invsigma
          end if
       end if
    end do
!    tod%scans(scan_id)%d(det)%dgain      = sum(stot_invN * residual)
!    tod%scans(scan_id)%d(det)%gain_invsigma = sum(stot_invN * s_tot)
!    if (tod%scans(scan_id)%d(det)%gain_invsigma < 0) then
!       write(*,*) 's', sum(mask * stot_invN * s_tot), sum(stot_invN * s_tot)
!    end if

    !write(*,*) det, scan_id, real(tod%scans(scan_id)%d(det)%dgain/tod%scans(scan_id)%d(det)%gain_invsigma,sp), real(tod%gain0(0),sp), real(tod%gain0(det),sp), real(g_old,sp)

   ! write(*,*) tod%scanid(scan_id), real(tod%scans(scan_id)%d(1)%dgain/tod%scans(scan_id)%d(3)%gain_invsigma,sp), real(tod%gain0(0) + tod%gain0(3) + tod%scans(scan_id)%d(3)%dgain/tod%scans(scan_id)%d(3)%gain_invsigma,sp), '# deltagain'

    !if (.false. .and. trim(tod%freq) == '030' .and. mod(tod%scanid(scan_id),100) == 0) then
    if (.false.) then
       call int2string(tod%scanid(scan_id), itext)
       !write(*,*) 'gain'//itext//'   = ', tod%gain0(0) + tod%gain0(1), tod%scans(scan_id)%d(1)%dgain/tod%scans(scan_id)%d(1)%gain_invsigma
       open(58,file='gain_delta_'//itext//'.dat')
       do i = ext(1), ext(2)
          write(58,*) i, residual(i,1)
       end do
       write(58,*)
       write(58,*)
       do i = 1, size(s_invsqrtN,1)
          write(58,*) i, s_invsqrtN(i,1)
       end do
       write(58,*)
       do i = 1, size(s_tot,1)
          write(58,*) i, tod_arr(i, 1) - (tod%gain0(0) +  tod%gain0(1)) * s_tot(i,1)
       end do
       write(58,*)
       do i = 1, size(s_tot,1)
          write(58,*) i, tod_arr(i, 1)
       end do
       close(58)
    end if

    deallocate(residual, r_fill)

  end subroutine calculate_gain_mean_std_per_scan

! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  ! Eirik: Update this routine to sample time-dependent gains properly; results should be stored in self%scans(i)%d(j)%gain, with gain0(0) and gain0(i) included
  module subroutine sample_smooth_gain(tod, handle, dipole_mods, smooth)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),   dimension(:, :),       intent(in)     :: dipole_mods
    logical(lgt), optional,            intent(in)     :: smooth


!    real(sp), dimension(:, :) :: inv_gain_covar ! To be replaced by proper matrix, and to be input as an argument
    integer(i4b) :: i, j, k, l, ndet, nscan_tot, ierr, ind(1)
    integer(i4b) :: currstart, currend, window, i1, i2, pid_id, range_end
    real(dp)     :: mu, denom, sum_inv_sigma_squared, sum_weighted_gain, g_tot, g_curr, sigma_curr, fknee, sigma_0, alpha, invvar, minvar, var, x, delta
    real(dp), allocatable, dimension(:)     :: lhs, rhs, g_smooth, a, xx, yy
    real(dp), allocatable, dimension(:)     :: g_over_s, temp_invsigsquared
    real(dp), allocatable, dimension(:)     :: summed_invsigsquared, smoothed_gain
    real(dp), allocatable, dimension(:,:,:) :: g
    real(dp), allocatable, dimension(:,:)   :: pidrange_gainarr
    integer(i4b),   allocatable, dimension(:, :) :: window_sizes
    integer(i4b), save :: cnt = 0
    character(len=128)  :: kernel_type
    logical(lgt) :: smooth_, sample_per_jump

    smooth_ = .true.
    if (present(smooth) ) smooth_ = smooth

    ndet       = tod%ndet
    nscan_tot  = tod%nscan_tot

    ! Collect all gain estimates on the root processor
    allocate(g(nscan_tot,ndet,2))
    allocate(g_over_s(nscan_tot))
    allocate(lhs(nscan_tot), rhs(nscan_tot))
    g = 0.d0
    do j = 1, ndet
       do i = 1, tod%nscan
          k        = tod%scanid(i)
          if (.not. tod%scans(i)%d(j)%accept) cycle
          g(k,j,1) = tod%scans(i)%d(j)%dgain
          g(k,j,2) = tod%scans(i)%d(j)%gain_invsigma
       end do
    end do

    if (.not. smooth) then
       ! If the time-dependent gain has high signal-to-noise per scan, then we
       ! don't need to Wiener filter, and instead implement equation (28) of
       ! Gjerlow et al. 2020.
       do j = 1, ndet
          do i = 1, tod%nscan
             k        = tod%scanid(i)
             if (.not. tod%scans(i)%d(j)%accept) cycle
             tod%scans(i)%d(j)%dgain = g(k,j,1)/g(k,j,2)
             if (trim(tod%operation)=='sample') tod%scans(i)%d(j)%dgain = tod%scans(i)%d(j)%dgain + rand_gauss(handle)/sqrt(g(k,j,2))
             tod%scans(i)%d(j)%gain  = tod%gain0(0) + tod%gain0(j) + tod%scans(i)%d(j)%dgain
          end do
       end do

       ! Perform poly-fit
       allocate(xx(tod%nscan_tot), yy(tod%nscan_tot))
       allocate(a(0:16))
       do j = 1, tod%ndet
          xx = 0.d0
          do i = 1, tod%nscan
             xx(tod%scanid(i)) = tod%scans(i)%d(j)%dgain
          end do
          call mpi_allreduce(xx, yy, size(xx), MPI_DOUBLE_PRECISION, &
               & MPI_SUM, tod%comm, ierr)

!!$          do i = 1, tod%nscan_tot
!!$             xx(i) = (i-1.d0)/tod%nscan_tot
!!$          end do

!!$          call fit_polynomial(xx, yy, a)
!!$          do i = 1, tod%nscan
!!$             x = (tod%scanid(i)-1.d0)/tod%nscan_tot
!!$             tod%scans(i)%d(j)%dgain = 0.d0
!!$             do k = 1, size(a)-1 
!!$                tod%scans(i)%d(j)%dgain = tod%scans(i)%d(j)%dgain + a(k)*x**k
!!$             end do
!!$             tod%scans(i)%d(j)%gain  = tod%gain0(0) + tod%gain0(j) + tod%scans(i)%d(j)%dgain
!!$          end do

          delta = sum(yy)/size(yy)
          do i = 1, tod%nscan
             tod%scans(i)%d(j)%dgain = tod%scans(i)%d(j)%dgain - delta
             tod%scans(i)%d(j)%gain  = tod%scans(i)%d(j)%gain  - delta
          end do

       end do
       deallocate(a, xx, yy)

    else
       call mpi_allreduce(mpi_in_place, g, size(g), MPI_DOUBLE_PRECISION, MPI_SUM, &
            & tod%comm, ierr)

       do j = 1+tod%myid, ndet, tod%numprocs
         if (all(g(:, j, 1) == 0)) continue
      
          ! Tune uncertainties to allow for proper compromise between smoothing and stiffness; set gain sigma_0 to minimum of the empirical variance
          if (count(g(:,j,2)>0) > 500) call compute_minimum_sigma0(g(:,j,2), 100, tod%gain_sigma_0(j))
          tod%gain_alpha(j) = -1.d0             ! Physically motivated value
          tod%gain_fknee(j) = tod%gain_samprate ! makes sigma_0 = true standard devation per sample

          !if (j == 1) then
          !   open(58,file='gain_in.dat')
          !   do k = 1, size(g,1)
          !      if (g(k,j,2) > 0) then
          !         write(58,*) k, g(k,j,1)/g(k,j,2), 1/sqrt(g(k,j,2))
          !      end if
          !   end do
          !   close(58)
             !write(*,*) '|  psd = ', tod%gain_sigma_0(j), tod%gain_alpha(j), tod%gain_fknee(j)

          !   open(68,file='g.unf', form='unformatted')
          !   write(68) size(g,1)
          !   write(68) g(:,j,1)
          !   write(68) g(:,j,2)
          !   write(68) tod%gain_sigma_0(j)
          !   write(68) tod%gain_alpha(j)
          !   write(68) tod%gain_fknee(j)
          !   close(68)
          !end if

          sample_per_jump = .false. .and. (size(tod%jumplist(j, :)) > 2)
          if (sample_per_jump) then
             if (tod%myid == 0) then
                write(*, *) 'Estimating per gain jump'
             end if
             do k = 1, size(tod%jumplist(j, :)) - 2
                call wiener_filtered_gain(g(tod%jumplist(j, k):tod%jumplist(j, k+1), j, 1), g(tod%jumplist(j, k):tod%jumplist(j, k+1), j, 2), &
                     & tod%gain_samprate, tod%gain_sigma_0(j), tod%gain_alpha(j), tod%gain_fknee(j), trim(tod%operation)=='sample', handle)
             end do
               ! Final chunk
             call wiener_filtered_gain(g(tod%jumplist(j, k):, j, 1), g(tod%jumplist(j, k):, j, 2), &
                  & tod%gain_samprate, tod%gain_sigma_0(j), tod%gain_alpha(j), tod%gain_fknee(j), trim(tod%operation)=='sample', handle)
          else
             call wiener_filtered_gain(g(:, j, 1), g(:, j, 2), tod%gain_samprate, tod%gain_sigma_0(j), tod%gain_alpha(j), &
                & tod%gain_fknee(j), trim(tod%operation)=='sample', handle)
          end if

         ! Force flat average to zero
         mu = 0.d0
         denom = 0.d0
         do k = 1, nscan_tot
            if (g(k, j, 2) > 0.d0) then
               mu    = mu    + g(k,j,1)
               denom = denom + 1.d0
            end if
         end do
         mu = mu / denom
         !write(*,*) 'g = ', mu
         where(g(:,j,2) > 0.d0) 
            g(:,j,1) = g(:,j,1) - mu
         end where

          !if (j == 1) then
          !   open(58,file='gain_out.dat')
          !   do k = 1, size(g,1)
          !      if (g(k,j,2) > 0) then
          !         write(58,*) k, g(k,j,1), 1/sqrt(g(k,j,2))
          !      end if
          !   end do
          !   close(58)
          !end if
       end do
!!$       call mpi_finalize(ierr)
!!$       stop

       ! Distribute and update results
       do j = 1, ndet
         call mpi_bcast(tod%gain_sigma_0(j),              1, MPI_DOUBLE_PRECISION, mod(j-1,tod%numprocs), tod%comm, ierr)
         call mpi_bcast(g(:,j,:),            size(g(:,j,:)), MPI_DOUBLE_PRECISION, mod(j-1,tod%numprocs), tod%comm, ierr)    
       end do
       do j = 1, ndet
          do i = 1, tod%nscan
             k        = tod%scanid(i)
             !if (g(k, j, 2) <= 0.d0) cycle
             tod%scans(i)%d(j)%dgain = g(k,j,1)
             tod%scans(i)%d(j)%gain  = tod%gain0(0) + tod%gain0(j) + g(k,j,1)
             !write(*,*) j, k,  tod%gain0(0), tod%gain0(j), g(k,j,1), tod%scans(i)%d(j)%gain 
          end do
       end do
    end if

    deallocate(g, lhs, rhs, g_over_s)

  end subroutine sample_smooth_gain

  module subroutine compute_minimum_sigma0(g2, window, sigma0)
    implicit none
    real(dp), dimension(:), intent(in)    :: g2
    integer(i4b),           intent(in)    :: window
    real(dp),               intent(out)   :: sigma0

    integer(i4b) :: k, l, nscan_tot, n
    real(dp)     :: var, minvar, invvar

    nscan_tot = size(g2)

    ! Compute minimum variance 
    minvar = 0.d0
    do k = window+1, nscan_tot-window ! First, compute the minimum variance, as estimated over a gliding window
       invvar = 0.d0
       do l = -window, window
          if (g2(k+l) > 0) then
             invvar = invvar + g2(k+1)
          end if
       end do
       if (count(g2(k-window:k+window) > 0) > 0) then
          invvar = invvar / count(g2(k-window:k+window) > 0)
          if (invvar > minvar) minvar = invvar
       end if
    end do
    sigma0 = 1d0 * 1/sqrt(minvar)

  end subroutine compute_minimum_sigma0


   ! This is implementing equation 16, adding up all the terms over all the sums
   ! the sum i is over the detector.
  module subroutine accumulate_abscal(tod, scan, mask, s_sub, s_invsqrtN, A_abs, b_abs, handle, out, s_highres, mask_lowres, tod_arr)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    integer(i4b),                      intent(in)     :: scan
    real(sp),          dimension(:,:), intent(in)     :: mask, s_sub
    real(sp),          dimension(:,:), intent(in)     :: s_invsqrtN
    real(dp),          dimension(:),   intent(inout)  :: A_abs, b_abs
    type(planck_rng),                  intent(inout)  :: handle
    logical(lgt), intent(in) :: out
    real(sp),          dimension(:,:), intent(in), optional :: s_highres
    real(sp),          dimension(:,:), intent(in), optional :: mask_lowres
    real(sp),          dimension(:,:), intent(in), optional :: tod_arr
 
    real(sp), allocatable, dimension(:,:)     :: residual
    real(sp), allocatable, dimension(:)       :: r_fill
    real(dp)     :: A, b, scale, dA, db
    integer(i4b) :: i, j, ext(2), ndet, ntod
    character(len=5) :: itext

    ndet = tod%ndet
    ntod = size(s_sub,1)
!    p = 5000 ! Buffer width
!    q = size(mask)-p
    call tod%downsample_tod(s_sub(:,1), ext)    
    allocate(residual(ext(1):ext(2),ndet), r_fill(size(s_sub,1)))
    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) then
          residual(:,j) = 0.
          cycle
       end if
       r_fill = tod_arr(:,j) - s_sub(:,j)
       call fill_all_masked(r_fill, mask(:,j), ntod, trim(tod%operation) == 'sample', abs(real(tod%scans(scan)%d(j)%N_psd%sigma0, sp)), handle, tod%scans(scan)%chunk_num)
       call tod%downsample_tod(r_fill, ext, residual(:,j))
    end do

    call multiply_inv_N(tod, scan, residual, sampfreq=tod%samprate_lowres, pow=0.5d0)

    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (present(mask_lowres)) then
          dA = sum(s_invsqrtN(:,j) ** 2    * mask_lowres(:,j))
          db = sum(s_invsqrtN(:,j) * residual(:,j) * mask_lowres(:,j))
       else
          dA = sum(s_invsqrtN(:,j) ** 2)
          db = sum(s_invsqrtN(:,j) * residual(:,j))
       end if

       if (dA == 0.) then
          tod%scans(scan)%d(j)%accept = .false.
          write(*,*) '| Rejecting scan in gain due to no unmasked samples: ', tod%scanid(scan), j
!       else if (abs(db/sqrt(dA)) > 1000. .or. 1/sqrt(dA) < 1d-6) then
!          tod%scans(scan)%d(j)%accept = .false.
!          write(*,*) '| Rejecting scan in abs gain due to dubious uncertainty: ', tod%scanid(scan), j, db/dA, 1/sqrt(dA)
       else 
          A_abs(j) = A_abs(j) + dA
          b_abs(j) = b_abs(j) + db
          !write(*,*) scan, j, db/dA
       end if

       !if (tod%scanid(scan) == 30 .and. out) then
!       if (out .and. dA > 0.d0) then
!          write(*,*) tod%scanid(scan), real(db/dA,sp), real(1/sqrt(dA),sp), '  # abs', j
!         !write(*,*) tod%scanid(scan), sum(abs(s_invN(:,j))), sum(abs(residual(:,j))), sum(abs(s_ref(:,j))), '  # absK', j
!       else if (dA > 0.d0) then
!          write(*,*) tod%scanid(scan), real(db/dA,sp), real(1/sqrt(dA),sp), '  # rel', j, tod%gain0(0), tod%scans(scan)%d(j)%gain
!       end if
    end do

    deallocate(residual, r_fill)

  end subroutine accumulate_abscal

  ! Sample absolute gain from orbital dipole alone 
  module subroutine sample_abscal_from_orbital(tod, handle, A_abs, b_abs)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs

    integer(i4b) :: i, j, ierr
    real(dp), allocatable, dimension(:) :: A, b

    ! Collect contributions from all cores
    allocate(A(tod%ndet), b(tod%ndet))
    call mpi_reduce(A_abs, A, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)
    call mpi_reduce(b_abs, b, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)

    ! Compute gain update and distribute to all cores
    if (tod%myid == 0) then
       tod%gain0(0) = sum(b)/sum(A)
       if (trim(tod%operation) == 'sample') then
          ! Add fluctuation term if requested
          tod%gain0(0) = tod%gain0(0) + 1.d0/sqrt(sum(A)) * rand_gauss(handle)
       end if
       if (tod%verbosity > 1) then
         write(*,fmt='(a,f12.8)') ' |      abscal = ', tod%gain0(0)
         !write(*,*) 'sum(b), sum(A) = ', sum(b), sum(A)
       end if
    end if
    call mpi_bcast(tod%gain0(0), 1,  MPI_DOUBLE_PRECISION, 0, &
         & tod%info%comm, ierr)

    do j = 1, tod%nscan
       do i = 1, tod%ndet
          tod%scans(j)%d(i)%gain = tod%gain0(0) + tod%gain0(i) + tod%scans(j)%d(i)%dgain 
       end do
    end do

    deallocate(A, b)

  end subroutine sample_abscal_from_orbital

  ! Sample absolute gain from orbital dipole alone 
  module subroutine sample_relcal(tod, handle, A_abs, b_abs)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs

    integer(i4b) :: i, j, k, ierr
    integer(i4b), allocatable, dimension(:) :: ind
    real(dp),     allocatable, dimension(:) :: A, b, rhs, x, tmp
    real(dp),     allocatable, dimension(:, :) :: coeff_matrix

    ! Collect contributions from all cores
    allocate(A(tod%ndet), b(tod%ndet), rhs(tod%ndet+1), x(tod%ndet+1), ind(tod%ndet+1), tmp(tod%ndet+1))
    allocate(coeff_matrix(tod%ndet+1, tod%ndet+1))
    call mpi_reduce(A_abs, A, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)
    call mpi_reduce(b_abs, b, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)

    coeff_matrix = 0.d0
    rhs          = 0.d0
    x            = 0.d0
    ind          = 0     ! List of indices for active (non-zero) elements in linear system
    k            = 0
    if (tod%myid == 0) then
       do j = 1, tod%ndet
         coeff_matrix(j, j) = A(j)
         rhs(j) = b(j) 
         if (trim(tod%operation) == 'sample') rhs(j) = rhs(j) + sqrt(A(j)) * rand_gauss(handle)
         coeff_matrix(j, tod%ndet+1) = 0.5d0
         coeff_matrix(tod%ndet+1, j) = 1
         if (coeff_matrix(j,j) > 0.d0) then
            k      = k+1
            ind(k) = j
         end if
       end do
       k = k+1
       ind(k) = tod%ndet+1
       coeff_matrix(tod%ndet+1, tod%ndet+1) = 0.d0
       rhs(tod%ndet+1) = 0.d0
       call solve_system_real(coeff_matrix(ind(1:k),ind(1:k)), tmp(1:k), rhs(ind(1:k)))
       x(ind(1:k)) = tmp(1:k)
       if (tod%verbosity > 1) then
         write(*,*) '|      relcal = ', real(x(1:tod%ndet),sp)
       end if
    end if
    call mpi_bcast(x, tod%ndet+1, MPI_DOUBLE_PRECISION, 0, &
       & tod%info%comm, ierr)

    tod%gain0(1:tod%ndet) = x(1:tod%ndet)
    deallocate(coeff_matrix, rhs, A, b, x, ind, tmp)

    do j = 1, tod%nscan
       do i = 1, tod%ndet
          tod%scans(j)%d(i)%gain = tod%gain0(0) + tod%gain0(i) + tod%scans(j)%d(i)%dgain 
       end do
    end do

  end subroutine sample_relcal

  module subroutine sample_imbal_cal(tod, handle, A_abs, b_abs)
    !
    !  Subroutine to sample the transmission imbalance parameters, defined in
    !  the WMAP data model as the terms x_im; given the definition
    !  d_{A/B} = T_{A/B} \pm Q_{A/B} cos(2 gamma_{A/B}) \pm U_{A/B} sin(2 gamma_{A/B})
    !  we have
    !  d = g[(1+x_im)*d_A - (1-x_im)*d_B]
    !    = g(d_A - d_B) + g*x_im*(d_A + d_B)
    !  Returns x_{im,1} for detectors 13/14, and x_{im,2} for detectors 23/24.
    !
    !
    !  Arguments (fixed):
    !  ------------------
    !  A_abs: real(dp)
    !     Accumulated A_abs = s_ref^T N^-1 s_ref for all scans
    !  b_abs: real(dp)
    !     Accumulated b_abs = s_ref^T N^-1 s_sub for all scans
    !
    !  
    !  Arguments (modified):
    !  ---------------------
    !  tod: comm_WMAP_tod
    !     The entire tod object. tod%x_im estimated and optionally sampled
    !  handle: planck_rng derived type 
    !     Healpix definition for random number generation
    !     so that the same sequence can be resumed later on from that same point
    !
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs

    integer(i4b) :: i, j, ierr
    real(dp), allocatable, dimension(:) :: A, b

    ! Collect contributions from all cores
    allocate(A(tod%ndet), b(tod%ndet))
    call mpi_reduce(A_abs, A, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)
    call mpi_reduce(b_abs, b, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)

    ! Compute gain update and distribute to all cores
    if (tod%myid == 0) then
       tod%x_im(1) = sum(b(1:2))/sum(A(1:2))
       tod%x_im(3) = sum(b(3:4))/sum(A(3:4))
       if (trim(tod%operation) == 'sample') then
          ! Add fluctuation term if requested
          tod%x_im(1) = tod%x_im(1) + 1.d0/sqrt(sum(A(1:2))) * rand_gauss(handle)
          tod%x_im(3) = tod%x_im(3) + 1.d0/sqrt(sum(A(3:4))) * rand_gauss(handle)
       end if
       !tod%x_im = 0d0
       tod%x_im(2) = tod%x_im(1)
       tod%x_im(4) = tod%x_im(3)
       if (tod%verbosity > 1) then
         write(*,*) '|      imbal =', real(tod%x_im(1),sp), real(tod%x_im(3),sp)
       end if
    end if
    call mpi_bcast(tod%x_im, 4,  MPI_DOUBLE_PRECISION, 0, &
         & tod%info%comm, ierr)


    deallocate(A, b)

  end subroutine sample_imbal_cal

  module subroutine get_smoothing_windows(tod, windows, dipole_mods)
     implicit none

     integer(i4b), dimension(:, :), intent(out)         :: windows
     class(comm_tod),               intent(in)          :: tod
     real(dp), dimension(:, :),     intent(in)          :: dipole_mods

     real(sp), allocatable, dimension(:)                :: nonzero_mask
     integer(i4b)  :: window_size_dipole_minimum, window_size_dipole_maximum
     real(dp)   :: low_dipole_thresh, curr_dipole
     real(dp)   :: high_dipole_thresh

     real(dp)   :: reduce_fac 

     integer(i4b)   :: i, j, window_size


     high_dipole_thresh = 5d-6
     low_dipole_thresh = 1d-6
     reduce_fac = 1.d0

     select case (trim(tod%freq))
         case ('030')
            window_size_dipole_minimum = 1200
            window_size_dipole_maximum = 400
         case ('044')
            window_size_dipole_minimum = 1500
            window_size_dipole_maximum = 300
         case ('070')
            window_size_dipole_minimum = 1880
            window_size_dipole_maximum = 400
         case ('023-WMAP_K')
             window_size_dipole_minimum = 1 !120
             window_size_dipole_maximum = 1 ! 40
         case default
            window_size_dipole_minimum = 1 !1800
            window_size_dipole_maximum = 1 !400
      end select

      do i = 1, tod%nscan_tot
         do j = 1, tod%ndet
            curr_dipole = dipole_mods(i, j)
            if (curr_dipole == 0) then
               windows(j, i) = 0
            else if (curr_dipole < low_dipole_thresh) then
               windows(j, i) = window_size_dipole_minimum
            else if (curr_dipole > high_dipole_thresh) then
               windows(j, i) = window_size_dipole_maximum
            else
               ! Interpolate
               windows(j, i) = window_size_dipole_minimum + &
                  & int(real(window_size_dipole_maximum - & 
                  & window_size_dipole_minimum, dp) / (high_dipole_thresh - &
                  & low_dipole_thresh) * (curr_dipole - low_dipole_thresh))
            end if
            windows(j, i) = windows(j, i) * reduce_fac
         end do
      end do

  end subroutine get_smoothing_windows


  module subroutine wiener_filtered_gain(b, inv_N_wn, samprate, sigma_0, alpha, fknee, sample, &
     & handle)
     !
     ! Given a spectral model for the gain, samples a wiener filtered (smoothed)
     ! estimate of the gain given data. The basic equation to be solved for is
     !  Ax = b where
     !  A = inv_G + s_ref * inv_N_wn * s_ref, where inv_G is the covariance
     !                                        matrix of the gain, given by the
     !                                        spectral model.
     !  and 
     !  b = s_ref * inv_N_wn * residual.
     !
     ! If a sample is to be drawn, a fluctuation term will also be added.
     !
     ! The spectral model is given by P(f) = sigma_0**2 * (f/f_knee) ** alpha
     !
     ! Arguments:
     ! ----------
     ! b:           real(dp) array
     !              The "right hand side" of the equation Ax = b, where x is the
     !              wiener filtered solution. b = dgain from
     !              calculate_gain_mean_std_per_scan. Will contain the solution.
     ! inv_N_wn:    real(dp) array
     !              The inverse white noise diagonal covariance matrix
     ! sigma_0:     real(dp)
     !              The current estimate of sigma_0 in the spectral model
     ! alpha:       real(dp)
     !              The current estimate of alpha in the spectral model
     ! fknee:       real(dp)
     !              The current estimate of fknee in the spectral model
     ! sample:      logical(lgt)
     !              Whether to draw a sample from the distribution or just
     !              return the best-fit solution.
     ! handle:      derived class (planck_rng)
     !              Random number generator handle. Will be modified.
     !
     ! Returns:
     ! --------
     ! b:           real(dp) array
     !              At exit, will contain the solution to Ax = b, including
     !              fluctuations if sampling is turned on.

     implicit none

     real(dp), dimension(:), intent(inout)   :: b
     real(dp), dimension(:), intent(in)      :: inv_N_wn
     real(dp), intent(in)                    :: samprate, sigma_0, alpha, fknee
     logical(lgt), intent(in)                :: sample
     type(planck_rng)                        :: handle

     real(dp), allocatable, dimension(:)     :: freqs, dt, inv_N_corr, b2, inv_N_wn2
     complex(dpc), allocatable, dimension(:) :: dv
     real(dp), allocatable, dimension(:)     :: fluctuations, temp
     real(dp), allocatable, dimension(:)     :: precond
     complex(dpc), allocatable, dimension(:) :: fourier_fluctuations
     integer*8          :: plan_fwd, plan_back
     integer(i4b)       :: nscan, nfft, n, nomp, err
     real(dp)           :: sigma0_wn

     integer(i4b)   :: i

     nscan = size(b)
     nfft = 2*nscan
     n = nfft / 2 + 1
     !samprate = 1.d0 / (24.d0*60.d0 * 60.d0) ! Just assuming a pid per hour for now
     allocate(freqs(0:n-1))
     do i = 0, n - 1
        freqs(i) = i * (samprate * 0.5) / (n - 1)
     end do

     allocate(b2(nfft), inv_N_wn2(nfft))
     b2(1:nscan)                   = b
     b2(2*nscan:nscan+1:-1)        = b
     inv_N_wn2(1:nscan)            = inv_N_wn
     inv_N_wn2(2*nscan:nscan+1:-1) = inv_N_wn

     nomp = OMP_GET_THREAD_NUM()
     call dfftw_init_threads(err)
     call dfftw_plan_with_nthreads(nomp)
     allocate(dt(nfft), dv(0:n-1))
     call dfftw_plan_dft_r2c_1d(plan_fwd, nfft, dt, dv, fftw_estimate + fftw_unaligned)
     call dfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)
     ! We use temporary arrays inside the transformation function, so no need
     ! for these anymore
     deallocate(dt, dv)

     allocate(inv_N_corr(0:n-1))
     allocate(fluctuations(nfft))
     allocate(precond(0:n-1))
     allocate(fourier_fluctuations(0:n-1))

     !write(*, *) 'Sigma_0: ', sigma_0
     !write(*, *) 'alpha: ', alpha
     !write(*, *) 'fknee: ', fknee

     call calculate_invcov(sigma_0, alpha, fknee, freqs, inv_N_corr)
     if (sample) then
        fourier_fluctuations = 0.d0
        do i = 1, n-1
           fourier_fluctuations(i) = cmplx(rand_gauss(handle),rand_gauss(handle))/sqrt(2.d0) * sqrt(inv_N_corr(i))
        end do
        call fft_back(fourier_fluctuations, fluctuations, plan_back)

        do i = 1, nfft
           if (inv_N_wn2(i) > 0) then
              fluctuations(i) = fluctuations(i) + sqrt(inv_N_wn2(i)) * rand_gauss(handle)
           end if
           b2(i) = b2(i) + fluctuations(i)
         end do
      end if

      !write(*,*) 'precond = ', maxval(inv_N_wn), median(inv_N_wn)
      do i = 0, n-1
         precond(i) = 1.d0/(inv_N_corr(i) + maxval(inv_N_wn2))
         !write(*,*) i, inv_N_corr(i), maxval(inv_N_wn), precond(i)
         !precond(i) = 1.d0/inv_N_wn(i)
      end do
      !write(*,*) i, inv_N_corr(i), maxval(inv_N_wn), precond(i)

      b2 = solve_cg_gain(inv_N_wn2, inv_N_corr, b2, precond, plan_fwd, plan_back)
      deallocate(inv_N_corr, freqs, fluctuations, fourier_fluctuations, precond)
        
      b = b2(1:nscan)
  
      call dfftw_destroy_plan(plan_fwd)                                           
      call dfftw_destroy_plan(plan_back) 

     deallocate(b2, inv_N_wn2)

  end subroutine wiener_filtered_gain

  module subroutine fft_fwd(time, fourier, plan_fwd)
     !
     ! Given a fft plan, transforms a vector from time domain to Fourier domain.
     !
     ! Arguments:
     ! ----------
     ! time:            real(dp) array
     !                  The original time-domain vector.
     ! fourier:         complex(dpc) array
     !                  The array that will contain the Fourier vector.
     ! plan_fwd:        integer*8
     !                  The fft plan to carry out the transformation
     !
     ! Returns:
     ! --------
     ! fourier_vector:  complex(dpc) array
     !                  At exit will contain the Fourier domain vector.

     implicit none
     real(dp),     dimension(1:), intent(in)     :: time
     complex(dpc), dimension(0:), intent(out)    :: fourier
     integer*8,                   intent(in)     :: plan_fwd
     real(dp),     dimension(1:size(time))       :: dt
     integer(i4b)                                :: ntime

     ntime = size(time)
     dt(1:ntime) = time
     !dt(2*ntime:ntime+1:-1) = dt(1:ntime)
     call timer%start(TOT_FFT)
     call dfftw_execute_dft_r2c(plan_fwd, dt, fourier)
     call timer%stop(TOT_FFT)
     fourier = fourier/sqrt(real(size(dt),dp))

   end subroutine fft_fwd

  module subroutine fft_back(fourier, time, plan_back)
     !
     ! Given a fft plan, transforms a vector from Fourier domain to time domain.
     !
     ! Arguments:
     ! ----------
     ! fourier:              complex(dpc) array
     !                  The original Fourier domain vector.
     ! time:            real(dp) array
     !                  The vector that will contain the transformed vector.
     ! plan_back:       integer*8
     !                  The fft plan to carry out the transformation.
     !
     ! Returns:
     ! --------
     ! vector:          real(dp) array
     !                  At exit will contain the time-domain vector.
     implicit none
     complex(dpc), dimension(0:), intent(in)     :: fourier
     real(dp),     dimension(1:), intent(out)    :: time
     real(dp),     dimension(size(time))         :: dt
     integer*8,                   intent(in)     :: plan_back

     call timer%start(TOT_FFT)
     call dfftw_execute_dft_c2r(plan_back, fourier, dt)
     call timer%stop(TOT_FFT)
     dt = dt/sqrt(real(size(dt),dp))
     time = dt(1:size(time))

   end subroutine fft_back

  module function solve_cg_gain(inv_N_wn, inv_N_corr, b, precond, plan_fwd, plan_back) result(solution) !, &
!     & with_precond)
     !
     ! Specialized function for solving the gain Wiener filter equation using
     ! the CG algorithm.
     !
     ! Arguments:
     ! ----------
     ! inv_N_wn:    real(dp) array
     !              The inverse white noise covariance matrix in time domain.
     ! inv_N_corr:  real(dp) array
     !              The inverse gain covariance matrix in Fourier domain.
     ! b:           real(dp) array
     !              The right hand side of the Wiener filter equation.
     ! precond:     real(dp) array
     !              The preconditioner matrix in time domain.
     ! plan_fwd:    integer*8
     !              The FFT forward plan.
     ! plan_back:   integer*8
     !              The FFT backward plan.
     !
     ! Returns:
     ! --------
     ! solve_cg_gain real(dp) array
     !               The solution to the Wiener filter equation.
     implicit none
     real(dp), dimension(1:), intent(in) :: inv_N_wn, b
     real(dp), dimension(0:), intent(in) :: inv_N_corr, precond
     integer*8             ,  intent(in) :: plan_fwd, plan_back
     real(dp), dimension(size(b))       :: solution
!     logical(lgt)                       :: with_precond


     real(dp), allocatable, dimension(:)    :: initial_guess, prop_sol, residual
     real(dp), allocatable, dimension(:)    :: p, Abyp, new_residual, z, new_z
     real(dp)       :: alpha, beta, orig_residual
     integer(i4b)   :: iterations, nscan, i
     logical(lgt)   :: converged
     character(len=8)   :: itext

!      open(58, file='gain_cg_invcov.dat')
!      do i = 1, size(inv_N_corr)
!         write(58, *) inv_N_corr(i)
!      end do
!      close(58)
!      open(58, file='gain_cg_invnwn.dat')
!      do i = 1, size(inv_N_wn)
!         write(58, *) inv_N_wn(i)
!      end do
!      close(58)


     nscan = size(b)
     allocate(initial_guess(nscan), prop_sol(nscan), residual(nscan), p(nscan))
     allocate(Abyp(nscan), new_residual(nscan), z(nscan), new_z(nscan))

     initial_guess = 0.d0
     prop_sol = initial_guess
     residual = b - tot_mat_mul_by_vector(inv_N_wn, inv_N_corr, prop_sol, &
        & plan_fwd, plan_back)
     orig_residual = sum(abs(residual))
     !z = residual / precond
     call apply_cg_precond(residual, z, precond, plan_fwd, plan_back)
!      open(58, file='gain_cg_residual.dat')
!      do i = 1, nscan
!         write(58, *) residual(i)
!      end do
!      close(58)

     p = z
     iterations = 0
     converged = .false.
     do while ((.not. converged) .and. (iterations < 100000))
         Abyp = tot_mat_mul_by_vector(inv_N_wn, inv_N_corr, p, plan_fwd, &
            & plan_back, filewrite=(iterations == 0))
!         if (iterations == 0) then
!            open(58, file='gain_cg_abyp.dat')
!            do i = 1, nscan
!               write(58, *) Abyp(i)
!            end do
!            close(58)
!         end if

         alpha = sum(residual * z) / sum(p * Abyp) 
         prop_sol = prop_sol + alpha * p
!         if (iterations == 0) then
!            write(*,*) 'Alpha:', alpha
!            open(58, file='gain_cg_prop_sol.dat')
!            do i = 1, nscan
!               write(58, *) prop_sol(i)
!            end do
!            close(58)
!         end if

         new_residual = residual - alpha * Abyp
!         if (iterations == 0) then
!            open(58, file='gain_cg_new_residual.dat')
!            do i = 1, nscan
!               write(58, *) new_residual(i)
!            end do
!            close(58)
!         end if

         !if (mod(iterations,1000) == 0) write(*,*) iterations, sum(abs(new_residual))/orig_residual
         if (sum(abs(new_residual))/orig_residual < 1e-6) then 
            converged = .true.
            exit
         end if
         !new_z = new_residual / precond
         call apply_cg_precond(new_residual, new_z, precond, plan_fwd, plan_back)
         beta = sum(new_residual * new_z) / sum(residual * z)
         p = new_z + beta * p
!         if (iterations == 0) then
!            write(*,*) 'Beta:', beta
!            open(58, file='gain_cg_p.dat')
!            do i = 1, nscan
!               write(58, *) p(i)
!            end do
!            close(58)
!         end if

         residual = new_residual
         z = new_z
         iterations = iterations + 1
         if (.true. .and. mod(iterations, 1) == 0) then
!            write(*, *) "Gain CG search res: ", sum(abs(new_residual))/orig_residual, sum(abs(new_residual))
!            call int2string(iterations, itext)
!            open(58, file='gain_cg_' // itext // '.dat')
!            do i = 1, nscan
!               write(58, *) prop_sol(i)
!            end do
!            close(58)
         end if
!         if (iterations == 1) then
!            open(58, file='gain_cg_' // itext // '.dat')
!            call int2string(iterations, itext)
!            do i = 1, nscan
!               write(58, *) prop_sol(i)
!            end do
!            close(58)
!         end if
     end do
     if (.not. converged) then
        write(*, *) "Gain CG search did not converge."
     end if
!     if (with_precond) then
!        write(*, *) "With preconditioner"
!     else
!        write(*, *) "Without preconditioner"
!     end if
     !write(*, *) "Gain CG iterations: ", iterations
     solution = prop_sol

     deallocate(initial_guess, prop_sol, residual, p, Abyp, new_residual, z, new_z)

  end function solve_cg_gain

  module subroutine apply_cg_precond(vec_in, vec_out, precond, plan_fwd, plan_back)
    implicit none
    real(dp), dimension(1:), intent(in)  :: vec_in
    real(dp), dimension(0:), intent(in)  :: precond
    real(dp), dimension(1:), intent(out) :: vec_out
    integer*8             , intent(in)  :: plan_fwd, plan_back    

    integer(i4b) :: i
    complex(dpc), dimension(size(precond))        :: fourier_vector

    call fft_fwd(vec_in, fourier_vector, plan_fwd)
    fourier_vector = fourier_vector * precond
    call fft_back(fourier_vector, vec_out, plan_back)

    !vec_out = vec_in !* precond

  end subroutine apply_cg_precond

  module function tot_mat_mul_by_vector(time_mat, fourier_mat, vector, plan_fwd, &
     & plan_back, filewrite) result(res)
     !
     ! Multiplies a sum of a time-domain diagonal matrix and a Fourier-domain
     ! diagonal matrix by a time-domain vector.
     !
     ! Arguments:
     ! ----------
     ! time_mat:    real(dp) array
     !              The matrix that is diagonal in time domain.
     ! fourier_mat: real(dp) array
     !              The matrix that is diagonal in Fourier domain.
     ! vector:      real(dp) array
     !              The vector to multiply by.
     ! plan_fwd:    integer*8
     !              The FFT forward plan.
     ! plan_back:   integer*8
     !              The FFT backward plan.
     ! filewrite:   logical(lgt), optional
     !              If present, writes various debug info into files.
     !
     ! Returns:
     ! --------
     ! tot_mat_mul_by_vector:   real(dp) array
     !                          The time-domain result of multiplying the
     !                          matrices by the input vector.
     implicit none

     real(dp), dimension(1:), intent(in)               :: vector
     real(dp), dimension(size(vector))                 :: res
     real(dp), dimension(size(vector)), intent(in)     :: time_mat
     real(dp), dimension(0:) , intent(in)     :: fourier_mat
     integer*8              , intent(in)     :: plan_fwd, plan_back
     logical(lgt)           , optional       :: filewrite

     logical(lgt)       :: write_file
     integer(i4b)       :: i
     real(dp), dimension(size(vector))       :: temp_vector
     complex(dpc), dimension(0:size(fourier_mat)-1)        :: fourier_vector

     write_file = .false.
     if (present(filewrite)) then
         write_file = filewrite
      end if

     if (all(vector .eq. 0)) then
        res = 0.d0
        return
     end if
!     if (write_file) then
!         open(58, file='gain_cg_vector_in.dat')
!         do i = 1, size(vector)
!            write(58, *) vector(i)
!         end do
!         close(58)
!      end if

     call fft_fwd(vector, fourier_vector, plan_fwd)
!     if (write_file) then
!         open(58, file='gain_cg_fourier_vector.dat')
!         do i = 1, size(fourier_vector)
!            write(58, *) abs(fourier_vector(i))
!         end do
!         close(58)
!      end if

     fourier_vector = fourier_vector * fourier_mat
!     if (write_file) then
!         open(58, file='gain_cg_matmul_fourier_vector.dat')
!         do i = 1, size(fourier_vector)
!            write(58, *) abs(fourier_vector(i))
!         end do
!         close(58)
!      end if


     call fft_back(fourier_vector, temp_vector, plan_back)
     !call fft_back(fourier_vector, temp_vector, plan_back, norm='transpose')
!     if (write_file) then
!         open(58, file='gain_cg_vector_after_fourier.dat')
!         do i = 1, size(temp_vector)
!            write(58, *) temp_vector(i)
!         end do
!         close(58)
!      end if

!     if (write_file) then
!         open(58, file='gain_cg_vector_time_matmul.dat')
!         do i = 1, size(vector)
!            write(58, *) vector(i) * time_mat(i)
!         end do
!         close(58)
!      end if

     !res = size(vector)*vector * time_mat + size(vector)**4*temp_vector
     !res = size(vector)*vector * time_mat + temp_vector
     !write(*,*) 'a', real(vector(1:4) * time_mat(1:4),sp)
     !write(*,*) 'b', real(temp_vector(1:4),sp)

     res = vector * time_mat + temp_vector
     !res = vector * time_mat

  end function tot_mat_mul_by_vector

  module subroutine calculate_invcov(sigma_0, alpha, fknee, freqs, invcov)
     !
     ! Calculate the inverse covariance matrix of the gain PSD given the PSD
     ! parameters as a diagonal matrix in Fourier domain.
     !
     ! The spectral model is given by P(f) = sigma_0**2 * (f/f_knee) ** alpha
     !
     ! Arguments:
     ! ----------
     ! sigma_0:     real(dp)
     !              The sigma_0 parameter of the spectral model. 
     ! alpha:       real(dp)
     !              The alpha parameter of the spectral model. 
     ! fknee:       real(dp)
     !              The fknee parameter of the spectral model. 
     ! freqs:       real(dp) array
     !              The frequencies at which to calculate the covariance matrix.
     !
     ! Returns:
     ! --------
     ! calculate_invcov:    real(dp) array
     !                      The covariance matrix in Fourier space, evaluated at
     !                      the frequencies given by 'freqs'.

     implicit none
     real(dp), dimension(0:), intent(in)     :: freqs
     real(dp), intent(in)                   :: sigma_0, fknee, alpha
     real(dp), dimension(0:), intent(out)    :: invcov
     real(dp)                               :: target_apod_freq
     integer(i4b)       :: apod_ind

     integer(i4b) :: i
     real(dp) :: apod

     invcov(0) = 0.d0
     target_apod_freq = 1d-6
      do i = 1, size(freqs) -1
         apod = (1 + freqs(i)/target_apod_freq)**5
         invcov(i) = min(1.d0 / (sigma_0 ** 2 * (1+(freqs(i)/fknee) ** alpha)) * apod, 1d12)
!        else
!           invcov(i) = 1.d12
!        end if
     end do

  end subroutine calculate_invcov

  module function calc_sigma_0(gain)
     !
     ! Function to estimate sigma_0 in the gain PSD model. Sigma_0 is estimated
     ! as the standard deviation of
     !      res(i) = gain(i) - gain(i-1)
     ! over all scans.
     !
     ! Arguments:
     ! ----------
     ! gain:        real(dp) array
     !              The gains per scan.
     !
     ! Returns:
     ! --------
     ! calc_sigma_0: real(dp)
     !               The estimated sigma_0.
     implicit none
     real(dp)       :: calc_sigma_0
     real(dp), dimension(:) :: gain

     real(dp), dimension(size(gain))    :: res
     real(dp)   :: std, mean
     integer(i4b)   :: i

     mean = gain(1) 
     res = 0.d0
     do i = 2, size(gain)
         res(i) = gain(i) - gain(i-1)
         mean = mean + res(i)
      end do
      mean = mean / size(gain)
      std = 0.d0
      do i = 2, size(gain)
         std = std + (res(i) - mean) ** 2
      end do
      if (std == 0) then
         open(58, file='gain_cg_std_gain.dat')
         do i = 1, size(gain)
            write(58, *) gain(i)
         end do
         close(58)
      end if
         
      calc_sigma_0 = sqrt(std / (2 * (size(gain) - 1)))
  end function calc_sigma_0

  module subroutine sample_gain_psd(tod, handle)
     !
     ! "Master" subroutine for gain PSD sampling. Uses the current gain solution
     ! and current values for the PSD parameters to sample new values for these
     ! parameters.
     !
     ! Arguments:
     ! ----------
     ! tod:         derived class (comm_tod)
     !              The TOD object that contains all relevant data. The fields
     !              containing the gain PSD parameters will be updated.
     ! handle:      derived class (planck_rng)
     !              Random number generator handle. Will be modified.
     !
     ! Returns:
     ! --------
     ! tod will be updated with a new set of values for the gain PSD model,
     ! drawn from their conditional distribution given the current gain
     ! solution.
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle

    real(dp), allocatable, dimension(:,:) :: g
    real(dp), allocatable, dimension(:)     :: freqs, dt, currgain, gain_ps
    complex(dpc), allocatable, dimension(:) :: dv, gain_fourier
    real(dp)    :: samprate
    integer(i4b)    :: nscan, nfft, n, ndet, nscan_tot, ierr
    integer(i4b)    :: i, j, k
    integer*8       :: plan_fwd, plan_back

    ndet       = tod%ndet
    nscan_tot  = tod%nscan_tot
    ! Collect gains on all processors
    allocate(g(nscan_tot,ndet))
    g = 0.d0
    do j = 1, ndet
       do i = 1, tod%nscan
          k        = tod%scanid(i)
          if (.not. tod%scans(i)%d(j)%accept) cycle
          g(k,j) = tod%scans(i)%d(j)%gain
       end do
    end do
    call mpi_allreduce(mpi_in_place, g, size(g), MPI_DOUBLE_PRECISION, MPI_SUM, &
         & tod%comm, ierr)

!          tod%scans(j)%d(i)%gain = tod%gain0(0) + tod%gain0(i) + tod%scans(j)%d(i)%dgain 
    nfft = nscan_tot !* 2
    n = nfft / 2 + 1
    samprate = 1.d0 / (60.d0 * 60.d0) ! Just assuming a pid per hour for now
    allocate(dt(nfft), dv(0:n-1))
    call dfftw_plan_dft_r2c_1d(plan_fwd, nfft, dt, dv, fftw_estimate + fftw_unaligned)
    call dfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)
    ! We use temporary arrays inside the transformation function, so no need
    ! for these anymore
    deallocate(dt, dv)

    allocate(freqs(0:n-1))
    do i = 0, n - 1
       freqs(i) = i * (samprate * 0.5) / (n - 1)
    end do
    allocate(currgain(nscan_tot))
    allocate(gain_ps(0:n-1), gain_fourier(0:n-1))
    do j = 1 + tod%myid, tod%ndet, tod%numprocs
!      currgain = 0.d0
!      do i = 1, tod%nscan
!         k        = tod%scanid(i)
!         currgain(k) = tod%scans(i)%d(j)%gain
!      end do
      call fft_fwd(g(:, j), gain_fourier, plan_fwd)
      gain_ps = abs(gain_fourier) ** 2
!       inv_N_corr = calculate_invcov(tod%gain_sigma_0(j), tod%gain_alpha(j), tod%gain_fknee(j), freqs)
      call sample_psd_params_by_mh(gain_ps, freqs, &
         & tod%gain_sigma_0(j), tod%gain_fknee(j), tod%gain_alpha(j), &
         & tod%gain_sigma0_std, tod%gain_fknee_std, tod%gain_alpha_std, handle)
    end do
    do j = 1, ndet
      call mpi_bcast(tod%gain_sigma_0(j), 1, MPI_DOUBLE_PRECISION, &
           & mod(j-1,tod%numprocs), tod%comm, ierr)
      call mpi_bcast(tod%gain_fknee(j), 1, MPI_DOUBLE_PRECISION, &
           & mod(j-1,tod%numprocs), tod%comm, ierr)
      call mpi_bcast(tod%gain_alpha(j), 1, MPI_DOUBLE_PRECISION, &
           & mod(j-1,tod%numprocs), tod%comm, ierr)
    end do
    deallocate(freqs, currgain, gain_ps)

  end subroutine sample_gain_psd

  module subroutine sample_psd_params_by_mh(gain_ps, freqs, sigma_0, fknee, alpha, &
     & sigma0_std, fknee_std, alpha_std, handle)
    ! 
    ! Uses the Metropolis-Hastings algorithm to draw a sample of the gain PSD
    ! parameters given the current gain solution.
    !
    ! The subroutine will do an initial run to estimate a proposal covariance
    ! matrix, then use that covariance matrix in the final run. The last sample
    ! of the chain is used as the drawn sample.
    !
    ! Arguments:
    ! ----------
    ! gain_ps:      real(dp) array
    !               The current gain power spectrum in fourier space.
    ! freqs:        real(dp) array
    !               The frequencies at which the power spectrum is calculated.
    ! sigma_0:      real(dp)
    !               The sigma_0 parameter of the spectral model. 
    ! fknee:        real(dp)
    !               The fknee parameter of the spectral model. Will be modified.
    ! alpha:        real(dp)
    !               The alpha parameter of the spectral model. Will be modified.
    ! fknee_std:    real(dp)
    !               The initial proposal standard deviation for fknee - will be
    !               re-estimated during the sampling.
    ! alpha_std:    real(dp)
    !               The initial proposal standard deviation for alpha - will be
    !               re-estimated during the sampling.
    ! handle:       derived class (planck_rng)
    !               Random number generator handle. Will be modified.
    !
    ! Returns:
    !---------
    ! fknee and alpha will contain new values sampled from the posterior given
    ! by the current gain power spectrum.
    implicit none

    real(dp), dimension(0:), intent(in)      :: gain_ps
    real(dp), dimension(0:), intent(in)      :: freqs
    real(dp), intent(inout)                 :: sigma_0, fknee, alpha
    real(dp), intent(in)                    :: sigma0_std, fknee_std, alpha_std
    type(planck_rng),                  intent(inout)  :: handle

    real(dp), dimension(3, 3)               :: propcov
    real(dp), allocatable, dimension(:, :)  :: samples
    integer(i4b)        :: i, nsamp_tot, nsamp_tune, nsamp_burnin

    nsamp_tot  = 5000
    nsamp_tune = 2000
    nsamp_burnin = 1000

    propcov = 0.d0
    propcov(1, 1) = sigma0_std ** 2
    propcov(2, 2) = fknee_std ** 2
    propcov(3, 3) = alpha_std ** 2

    allocate(samples(nsamp_tot, 3))
    samples = 0.d0
    call run_mh(sigma_0, fknee, alpha, propcov, nsamp_tune, samples, gain_ps, freqs, .true., handle)
    open(58, file='samples_first.dat')
    do i = 1, nsamp_tune
      write(58, *) real(samples(i, :),sp)
    end do
    close(58)
    call compute_covariance_matrix(samples(nsamp_burnin+1:nsamp_tune, :), propcov)
    sigma_0 = samples(nsamp_tune, 1)
    fknee = samples(nsamp_tune, 2)
    alpha = samples(nsamp_tune, 3)
    call run_mh(sigma_0, fknee, alpha, propcov, nsamp_tot, samples, gain_ps, freqs, .false., handle)
    open(58, file='samples_second.dat')
    do i = 1, nsamp_tot
      write(58, *) real(samples(i, :),sp)
    end do
    close(58)
    sigma_0 = samples(nsamp_tot, 1)
    fknee = samples(nsamp_tot, 2)
    alpha = samples(nsamp_tot, 3)
    write(*,fmt='(a,2f18.12,f8.3)') "Gain PSD {sigma0,fknee,alpha}:", sigma_0, fknee, alpha

  end subroutine sample_psd_params_by_mh

  module subroutine run_mh(sigma_0, fknee, alpha, propcov, num_samples, samples, &
     & gain_ps, freqs, adjust_scaling_full, handle)
  !
  ! The core Metropolis-Hastings routine used for gain PSD sampling.
  !
  ! Arguments:
  ! ----------
  ! fknee:                  real(dp)
  !                         The current value of fknee.
  ! alpha:                  real(dp)
  !                         The current value of alpha.
  ! propcov:                real(dp) 2x2 matrix
  !                         The covariance matrix of the proposal density.
  ! sigma_0:                real(dp)
  !                         The sigma_0 parameter.
  ! num_samples:            integer(i4b)
  !                         The number of samples to draw.
  ! samples:                real(dp) array
  !                         Output array for the samples.
  ! gain_ps:                real(dp) array
  !                         The current gain power spectrum.
  ! freqs:                  real(dp) array
  !                         The frequencies at which the gain power spectrum is
  !                         evaluated.
  ! adjust_scaling_full:    logical(lgt)
  !                         If true, will continue to adjust the scaling
  !                         parameter (the step length) throughout. If false,
  !                         will only adjust the step length during the first
  !                         half of the sampling.
  ! handle:                 derived class (planck_rng)
  !                         Random number generator handle. Will be modified.
  !
  ! Returns:
  ! --------
  ! the 'samples' array will contain num_samples samples drawn from the
  ! posterior distribution of the gain PSD given the current gain power
  ! spectrum. Only the last half of the samples are strictly speaking
  ! statistically proper samples, and only if adjust_scaling_full is false.
    implicit none

    real(dp), intent(in)       :: fknee, alpha, sigma_0
    real(dp), dimension(3, 3), intent(in)  :: propcov
    integer(i4b), intent(in)       :: num_samples
    real(dp), dimension(:, :), intent(out) :: samples
    real(dp), dimension(0:),    intent(in)  :: freqs  
    real(dp), dimension(0:), intent(in)  :: gain_ps
    logical(lgt)          ,    intent(in)  :: adjust_scaling_full
    type(planck_rng),                  intent(inout)  :: handle


    integer(i4b)       :: i, accepted
    real(dp)                :: scaling_factor
    real(dp), dimension(3, 3)   :: sqrt_cov
    real(dp), dimension(3)  :: curr_vec, prop_vec, eta
    real(dp)                :: curr_lnl, prop_lnl, acc_rate

    samples = 0.d0
    curr_vec = 0.d0
    prop_vec = 0.d0
    accepted = 0
    scaling_factor = 1.d0
    
    curr_vec(1) = sigma_0
    curr_vec(2) = fknee
    curr_vec(3) = alpha
    curr_lnL = psd_loglike(curr_vec(1), curr_vec(2), curr_vec(3), gain_ps, freqs)
!    write(*, *) "Curr_lnl: ", curr_lnL
    sqrt_cov = propcov
    call compute_hermitian_root_with_mask(sqrt_cov, 0.5d0)
    do i = 1, num_samples
      eta = [rand_gauss(handle), rand_gauss(handle), rand_gauss(handle)]
      prop_vec = curr_vec + scaling_factor * matmul(sqrt_cov, eta)
      prop_lnL = psd_loglike(prop_vec(1), prop_vec(2), prop_vec(3), gain_ps, freqs)
      !write(*, *) "prop_lnL", prop_lnL
!      if (prop_lnL /= -1d30 .and. rand_uni(handle) <= exp(prop_lnL - curr_lnL)) then
      !write(*,*) prop_vec
      if (log(rand_uni(handle)) <= prop_lnL - curr_lnL) then
         !write(*,*) prop_lnL, curr_lnL, 'accept'
         curr_vec = prop_vec
         curr_lnL = prop_lnL
         accepted = accepted + 1
      else
         !write(*,*) prop_lnL, curr_lnL, 'reject'
      end if
      samples(i, :) = curr_vec
      if (mod(i, 100) == 0) then
         acc_rate = real(accepted, dp) / real(i, dp)
         if ((adjust_scaling_full) .or. (i < int(num_samples / 2))) then
            if (acc_rate > 0.5) then
               scaling_factor = scaling_factor * 2
            else if (acc_rate < 0.2) then
               scaling_factor = scaling_factor * 0.5
            end if
         end if
      end if
    end do
!    write(*, *) "Final acceptance rate: ", acc_rate
!    write(*, *) "Final scaling factor: ", scaling_factor

  end subroutine run_mh

  module function psd_loglike(sigma_0, fknee, alpha, gain_ps, freqs)
     !
     ! Calculates the log-likelihood given the gain and psd parameters.
     !
     ! The spectral model is given by P(f) = sigma_0**2 * (f/f_knee) ** alpha,
     ! and the log-likelihood is given by taking the sum over all frequencies of
     ! the gain power spectrum times the inverse of the covariance matrix given
     ! by the spectral model, minus the logarithm of this covariance matrix.
     !
     ! Arguments:
     ! ----------
     ! fknee:       real(dp)
     !              The value of the fknee parameter.
     ! alpha:       real(dp)
     !              The value of the alpha parameter.
     ! sigma_0:     real(dp)
     !              The value of sigma_0
     ! gain_ps:     real(dp) array
     !              The power spectrum of the current gain solution.
     ! freqs:       real(dp) array
     !              The frequencies at which the gain power spectrum is
     !              evaluated.
     !
     ! Returns:
     ! --------
     ! psd_loglike      real(dp)
     !                  The log-likelihood value of the parameter combination
     !                  given the current gain power spectrum.
     implicit none

     real(dp)        :: psd_loglike
     real(dp)        :: fknee, alpha, sigma_0
     real(dp), dimension(0:)  :: gain_ps, freqs

     real(dp) :: lambda
     real(dp), dimension(size(gain_ps))  :: inv_N_corr

     if (sigma_0 <= 0.d0 .or. fknee <= 0.d0 .or. alpha < -3.d0) then
        psd_loglike = -1d30
        return
     end if
     !lambda = 1e8
     call calculate_invcov(sigma_0, alpha, fknee, freqs, inv_N_corr)
     
     psd_loglike = -sum(gain_ps(1:2000) * inv_N_corr(1:2000) - log(inv_N_corr(1:2000))) !- lambda*sigma_0
     !write(*,*) sigma_0, sum(gain_ps(2:) * inv_N_corr(2:) - log(inv_N_corr(2:)))!,  lambda*sigma_0

  end function psd_loglike
    

end submodule comm_tod_gain_smod
