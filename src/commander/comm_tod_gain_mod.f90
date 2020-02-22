module comm_tod_gain_mod
  use comm_tod_mod
  use comm_tod_noise_mod
  use comm_utils
  implicit none



contains


  ! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  ! Haavard: Get rid of explicit n_corr, and replace 1/sigma**2 with proper invN multiplication
!  subroutine sample_gain_per_scan(self, handle, det, scan_id, n_corr, mask, s_ref)
!   subroutine calculate_gain_mean_std_per_scan(self, det, scan_id, s_tot, invn, mask)
   subroutine calculate_gain_mean_std_per_scan(tod, scan_id, s_invN, mask, s_ref, s_tot)
    implicit none
    class(comm_tod),                      intent(inout) :: tod
    real(sp),             dimension(:,:), intent(in)    :: s_invN, mask, s_ref, s_tot
    integer(i4b),                         intent(in)    :: scan_id

    real(sp), allocatable, dimension(:,:) :: residual
    real(sp), allocatable, dimension(:)   :: r_fill
    integer(i4b) :: ext(2), j, i, ndet, ntod
    character(len=5) :: itext

    ndet = tod%ndet
    ntod = size(s_tot,1)

!    allocate(residual(size(stot_invN)))
      !g_old = tod%scans(scan_id)%d(det)%gain 

!   residual = (tod%scans(scan_id)%d(det)%tod - (tod%gain0(0) + &
!         & tod%gain0(det)) * s_tot) * mask

    allocate(r_fill(size(s_tot,1)))
    call tod%downsample_tod(s_tot(:,1), ext)    
    allocate(residual(ext(1):ext(2),ndet))
    do j = 1, ndet
       if (.not. tod%scans(scan_id)%d(j)%accept) then
          residual(:,j) = 0.d0
          cycle
       end if
       r_fill = tod%scans(scan_id)%d(j)%tod - (tod%gain0(0) + &
            & tod%gain0(j)) * s_tot(:,j)
       call fill_all_masked(r_fill, mask(:,j), ntod, .false.)
       call tod%downsample_tod(r_fill, ext, residual(:,j))
    end do
    call multiply_inv_N(tod, scan_id, residual, sampfreq=tod%samprate_lowres, pow=0.5d0)

    ! Get a proper invn multiplication from Haavard's routine here
    do j = 1, ndet
       if (.not. tod%scans(scan_id)%d(j)%accept) then
          tod%scans(scan_id)%d(j)%gain  = 0.d0
          tod%scans(scan_id)%d(j)%dgain = 0.d0
       else
          tod%scans(scan_id)%d(j)%dgain      = sum(s_invN(:,j) * residual(:,j))
          tod%scans(scan_id)%d(j)%gain_sigma = sum(s_invN(:,j) * s_ref(:,j))
          if (tod%scans(scan_id)%d(j)%gain_sigma < 0.d0) then
             write(*,*) 'Warning: Not positive definite invN = ', tod%scanid(scan_id), j, tod%scans(scan_id)%d(j)%gain_sigma
          end if
       end if
    end do
!    tod%scans(scan_id)%d(det)%dgain      = sum(stot_invN * residual)
!    tod%scans(scan_id)%d(det)%gain_sigma = sum(stot_invN * s_tot)
!    if (tod%scans(scan_id)%d(det)%gain_sigma < 0) then
!       write(*,*) 's', sum(mask * stot_invN * s_tot), sum(stot_invN * s_tot)
!    end if

!!$    if (tod%scans(scan_id)%d(det)%dgain/tod%scans(scan_id)%d(det)%gain_sigma > 5.) then
!!$       open(58, file='tod.dat')
!!$       do i = p, q
!!$          if (mask(i) > 0.5) write(58,*) i, residual(i), s_tot(i)*0.001 
!!$       end do
!!$       close(58)
!!$       stop
!!$    end if

    !write(*,*) det, scan_id, real(tod%scans(scan_id)%d(det)%dgain/tod%scans(scan_id)%d(det)%gain_sigma,sp), real(tod%gain0(0),sp), real(tod%gain0(det),sp), real(g_old,sp)

   ! write(*,*) tod%scanid(scan_id), real(tod%scans(scan_id)%d(1)%dgain/tod%scans(scan_id)%d(3)%gain_sigma,sp), real(tod%gain0(0) + tod%gain0(3) + tod%scans(scan_id)%d(3)%dgain/tod%scans(scan_id)%d(3)%gain_sigma,sp), '# deltagain'

    if (.false. .and. trim(tod%freq) == '030' .and. mod(tod%scanid(scan_id),100) == 0) then
       call int2string(tod%scanid(scan_id), itext)
       !write(*,*) 'gain'//itext//'   = ', tod%gain0(0) + tod%gain0(1), tod%scans(scan_id)%d(1)%dgain/tod%scans(scan_id)%d(1)%gain_sigma
       open(58,file='gain_delta_'//itext//'.dat')
       do i = ext(1), ext(2)
          write(58,*) i, residual(i,1)
       end do
       write(58,*)
       do i = 1, size(s_ref,1)
          write(58,*) i, s_ref(i,1)
       end do
       write(58,*)
       do i = 1, size(s_ref,1)
          write(58,*) i, s_invN(i,1)
       end do
       write(58,*)
       do i = 1, size(s_tot,1)
          write(58,*) i, tod%scans(scan_id)%d(1)%tod(i) - (tod%gain0(0) + &
            & tod%gain0(1)) * s_tot(i,1)
       end do
       write(58,*)
       do i = 1, size(s_tot,1)
          write(58,*) i, tod%scans(scan_id)%d(1)%tod(i)
       end do
       close(58)
    end if


    deallocate(residual, r_fill)

  end subroutine calculate_gain_mean_std_per_scan


  ! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  ! Eirik: Update this routine to sample time-dependent gains properly; results should be stored in self%scans(i)%d(j)%gain, with gain0(0) and gain0(i) included
  subroutine sample_smooth_gain(tod, handle)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle

!    real(sp), dimension(:, :) :: inv_gain_covar ! To be replaced by proper matrix, and to be input as an argument
    integer(i4b) :: i, j, k, ndet, nscan_tot, ierr, ind(1)
    integer(i4b) :: currstart, currend, window, i1, i2
    real(dp)     :: mu, denom, sum_inv_sigma_squared, sum_weighted_gain, g_tot, g_curr, sigma_curr
    real(dp), allocatable, dimension(:)     :: lhs, rhs, g_smooth
    real(dp), allocatable, dimension(:,:,:) :: g

    ndet       = tod%ndet
    nscan_tot  = tod%nscan_tot

    ! Collect all gain estimates on the root processor
    allocate(g(nscan_tot,ndet,2))
    allocate(lhs(nscan_tot), rhs(nscan_tot))
    g = 0.d0
    do j = 1, ndet
       do i = 1, tod%nscan
          k        = tod%scanid(i)
          if (.not. tod%scans(i)%d(j)%accept) cycle
          g(k,j,1) = tod%scans(i)%d(j)%dgain
          g(k,j,2) = tod%scans(i)%d(j)%gain_sigma
       end do
    end do
    if (tod%myid == 0) then
       call mpi_reduce(mpi_in_place, g, size(g), MPI_DOUBLE_PRECISION, MPI_SUM, &
            & 0, tod%comm, ierr)
    else
       call mpi_reduce(g,            g, size(g), MPI_DOUBLE_PRECISION, MPI_SUM, &
            & 0, tod%comm, ierr)
    end if

    if (tod%myid == 0) then
!       nbin = nscan_tot / binsize + 1

!!$       open(58,file='gain.dat', recl=1024)
!!$       do j = 1, ndet
!!$          do k = 1, nscan_tot
!!$             !if (g(k,j,2) /= 0) then
!!$             if (g(k,j,2) > 0) then
!!$                write(58,*) j, k, real(g(k,j,1)/g(k,j,2),sp), real(g(k,j,1),sp), real(g(k,j,2),sp)
!!$             else
!!$                write(58,*) j, k, 0.
!!$             end if
!!$          end do
!!$          write(58,*)
!!$       end do
!!$       close(58)

       do j = 1, ndet
         lhs = 0.d0
         rhs = 0.d0
         k = 0
         do while (k < nscan_tot)
            k         = k + 1
            currstart = k
            !if (g(k, j, 2) <= 0.d0) cycle
            if (g(k,j,2) > 0.d0) then
               sum_inv_sigma_squared = g(k, j, 2)
               sum_weighted_gain     = g(k, j, 1) !/ g(k, j, 2)
            else
               sum_inv_sigma_squared = 0.d0
               sum_weighted_gain     = 0.d0
            end if
            g_curr                = 0.d0
            sigma_curr            = 1.d0
            !do while (sqrt(sum_inv_sigma_squared) < 10000. .and. k < nscan_tot)
            do while (g_curr / sigma_curr < 200.d0 .and. k < nscan_tot)
               k = k + 1
               if (g(k,j,2) > 0.d0) then
                  sum_weighted_gain     = sum_weighted_gain     + g(k, j, 1) !/ g(k, j, 2)
                  sum_inv_sigma_squared = sum_inv_sigma_squared + g(k, j, 2)
               end if
               if (sum_inv_sigma_squared > 0.d0) then
                  g_curr     = tod%gain0(0) + tod%gain0(j) + sum_weighted_gain / sum_inv_sigma_squared
                  sigma_curr = 1.d0 / sqrt(sum_inv_sigma_squared)
                  !write(*,*) j, ', S/N = ', real(g_curr/sigma_curr,sp)
               else
                  g_curr     = 0.d0
                  sigma_curr = 1.d0
               end if
            end do
            currend = k
            if (sum_inv_sigma_squared > 0.d0) then
               g_tot = sum_weighted_gain / sum_inv_sigma_squared
               if (trim(tod%operation) == 'sample') then
                  ! Add fluctuation term if requested
                  g_tot = g_tot + rand_gauss(handle) / sqrt(sum_inv_sigma_squared)
               end if
            else
               g_tot = 0.d0
            end if
            !write(*,*) currstart, currend, g_tot, 1 / sqrt(sum_inv_sigma_squared)
            !write(*,*) currstart, currend, nscan_tot, g_tot 
            g(currstart:currend, j, 1) = g_tot 
            g(currstart:currend, j, 2) = 1.d0 / sqrt(sum_inv_sigma_squared)
         end do

         ! Doing binning for now, but should use proper covariance sampling
!         do i = 1, nbin
!            b1 = (i-1) * binsize + 1
!            b2 = min(i * binsize, nscan_tot)
!            if (count(g(b1:b2, j, 2) > 0) <= 1) cycle
!            n = 0
!            do k = b1, b2
!               if (g(k,j,2) <= 0.d0) then
!                  g(k,j,1) = 0.d0
!                  cycle
!               end if
!               if (g(k,j,1) == 0) cycle
!               n = n + 1
!               ! To be replaced by proper matrix multiplications
!               lhs(k) = g(k, j, 2)
!   !            rhs(k) = rhs(k) + g(k, j, 1) + sqrt(inv_gain_covar(j, k)) * rand_gauss(handle) + sqrt(g(k, j, 2)) * rand_gauss(handle) this is proper when we have the covariance matrix
!               ! HKE: Disabling fluctuations for now 
!               rhs(k) = g(k, j, 1) !+ sqrt(g(k, j, 2)) * rand_gauss(handle)
!!               g(k, j, 1) = rhs(k) / lhs(k)
!               vals(n) = rhs(k) / lhs(k)
!               if (abs(vals(n)) > 0.05d0) vals(n) = 0.d0
!            end do
!            where (g(b1:b2, j, 2) > 0)
!               g(b1:b2, j, 1) = median(vals(1:n))
!            end where
!         end do

         ! Apply running average smoothing
         allocate(g_smooth(nscan_tot))
         window = 500
         do k = 1, nscan_tot
            i1 = max(k-window,1)
            i2 = min(k+window,nscan_tot)
            g_smooth(k) = sum(g(i1:i2,j,1)) / (i2-i1)
            !g_smooth(k) = median(g(i1:i2,j,1)) 
         end do
         g(:,j,1) = g_smooth - mean(g_smooth)
         deallocate(g_smooth)
 
!!$         ! Apply running average smoothing
!!$         allocate(g_smooth(nscan_tot))
!!$         window = 10000
!!$         do k = 1, nscan_tot
!!$            i1 = max(k-window,1)
!!$            i2 = min(k+window,nscan_tot)
!!$            ind = minloc(g(i1:i2,j,2))
!!$            g_smooth(k) = g(ind(1),j,1)
!!$            !g_smooth(k) = median(g(i1:i2,j,1)) 
!!$         end do
!!$         g(:,j,1) = g_smooth - mean(g_smooth)
!!$         deallocate(g_smooth)

!!$         mu  = 0.d0
!!$         denom = 0.d0
!!$         do k = 1, nscan_tot
!!$            if (g(k, j, 2) <= 0.d0) cycle
!!$            mu         = mu + g(k, j, 1)! * g(k,j,2)
!!$            denom      = denom + 1.d0! * g(k,j,2)
!!$!            write(*,*) j, k, g(k,j,1), rhs(k)/lhs(k)
!!$         end do
!!$         mu = mu / denom
!!$
!!$         ! Make sure fluctuations sum up to zero
!!$         !write(*,*) 'mu = ', mu
!!$         g(:,j,1) = g(:,j,1) - mu
       end do
    end if

    ! Distribute and update results
    call mpi_bcast(g, size(g),  MPI_DOUBLE_PRECISION, 0, tod%comm, ierr)    
    do j = 1, ndet
       do i = 1, tod%nscan
          k        = tod%scanid(i)
          !if (g(k, j, 2) <= 0.d0) cycle
          tod%scans(i)%d(j)%dgain = g(k,j,1)
          tod%scans(i)%d(j)%gain  = tod%gain0(0) + tod%gain0(j) + g(k,j,1)
          !write(*,*) j, k,  tod%scans(i)%d(j)%gain 
       end do
    end do

!!$    call mpi_finalize(ierr)
!!$    stop

    deallocate(g, lhs, rhs)

  end subroutine sample_smooth_gain


  ! Haavard: Remove current monopole fit, and replace inv_sigmasq with a proper invN(alpha,fknee) multiplication
!  subroutine accumulate_abscal(self, scan, det, mask, s_sub, &
!       & s_orb, A_abs, b_abs)
   subroutine accumulate_abscal(tod, scan, mask, s_sub, s_ref, s_invN, A_abs, b_abs)
    implicit none
    class(comm_tod),                   intent(in)     :: tod
    integer(i4b),                      intent(in)     :: scan
    real(sp),          dimension(:,:), intent(in)     :: mask, s_sub, s_ref
    real(sp),          dimension(:,:), intent(in)     :: s_invN
    real(dp),          dimension(:),   intent(inout)  :: A_abs, b_abs

    real(sp), allocatable, dimension(:,:)     :: residual
    real(sp), allocatable, dimension(:)       :: r_fill
    real(dp)     :: A, b
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
       r_fill = tod%scans(scan)%d(j)%tod-s_sub(:,j)
       call fill_all_masked(r_fill, mask(:,j), ntod, .false.)
       call tod%downsample_tod(r_fill, ext, residual(:,j))
    end do

    call multiply_inv_N(tod, scan, residual, sampfreq=tod%samprate_lowres, pow=0.5d0)

    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       A_abs(j) = A_abs(j) + sum(s_invN(:,j) * s_ref(:,j))
       b_abs(j) = b_abs(j) + sum(s_invN(:,j) * residual(:,j))
    end do
    !write(*,*) sum(abs(s_sub)), sum(abs(tod%scans(scan)%d(det)%tod))

!!$    if (trim(tod%freq) == '070') then
!!$       write(*,*) tod%scanid(scan), real(b/A,sp), real(1/sqrt(A),sp), '  # abs70', det
!!$    end if
!!$

!!$    if (mod(tod%scanid(scan),100) == 0) then
!!$       call int2string(tod%scanid(scan), itext)
!!$       !write(*,*) 'gain'//itext//'   = ', tod%gain0(0) + tod%gain0(1), tod%gain0(0), tod%gain0(1)
!!$       open(58,file='gainfit3_'//itext//'.dat')
!!$       do i = 1, size(s_ref,1)
!!$          write(58,*) i, residual(i,1)
!!$       end do
!!$       write(58,*)
!!$       do i = 1, size(s_ref,1)
!!$          write(58,*) i, s_ref(i,1)
!!$       end do
!!$       write(58,*)
!!$       do i = 1, size(s_ref,1)
!!$          write(58,*) i, s_invN(i,1)
!!$       end do
!!$       write(58,*)
!!$       do i = 1, size(s_sub,1)
!!$          write(58,*) i, s_sub(i,1)
!!$       end do
!!$       write(58,*)
!!$       do i = 1, size(s_sub,1)
!!$          write(58,*) i, tod%scans(scan)%d(1)%tod(i)
!!$       end do
!!$       close(58)
!!$    end if



    !if (det == 1) write(*,*) tod%scanid(scan), b/A, '  # absgain', real(A,sp), real(b,sp)

    deallocate(residual, r_fill)

  end subroutine accumulate_abscal

  ! Sample absolute gain from orbital dipole alone 
  subroutine sample_abscal_from_orbital(tod, handle, A_abs, b_abs)
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
       write(*,*) 'abscal = ', tod%gain0(0)
    end if
    call mpi_bcast(tod%gain0(0), tod%ndet,  MPI_DOUBLE_PRECISION, 0, &
         & tod%info%comm, ierr)

    do j = 1, tod%nscan
       do i = 1, tod%ndet
          tod%scans(j)%d(i)%gain = tod%gain0(0) + tod%gain0(i) + tod%scans(j)%d(i)%dgain 
       end do
    end do

!!$    call mpi_finalize(ierr)
!!$    stop

    deallocate(A, b)

  end subroutine sample_abscal_from_orbital

  ! Sample absolute gain from orbital dipole alone 
  ! Eirik: Change this routine to sample the constant offset per radiometer; put the result in tod%gain0(i)
  subroutine sample_relcal(tod, handle, A_abs, b_abs)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs

    integer(i4b) :: i, j, ierr
    real(dp), allocatable, dimension(:) :: A, b, rhs, x
    real(dp), allocatable, dimension(:, :) :: coeff_matrix

    ! Collect contributions from all cores
    allocate(A(tod%ndet), b(tod%ndet), rhs(tod%ndet+1), x(tod%ndet+1))
    allocate(coeff_matrix(tod%ndet+1, tod%ndet+1))
    call mpi_reduce(A_abs, A, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)
    call mpi_reduce(b_abs, b, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)

    coeff_matrix = 0.d0
    rhs = 0.d0
    if (tod%myid == 0) then
       do j = 1, tod%ndet
         coeff_matrix(j, j) = A(j)
         rhs(j) = b(j) + sqrt(A(j)) * rand_gauss(handle)
         coeff_matrix(j, tod%ndet+1) = 0.5d0
         coeff_matrix(tod%ndet+1, j) = 1
       end do
       coeff_matrix(tod%ndet+1, tod%ndet+1) = 0.d0
       rhs(tod%ndet+1) = 0.d0
       call solve_system_real(coeff_matrix, x, rhs)
       
       write(*,*) 'relcal = ', real(x,sp)
    end if
    call mpi_bcast(x, tod%ndet+1, MPI_DOUBLE_PRECISION, 0, &
       & tod%info%comm, ierr)

    tod%gain0(1:tod%ndet) = x(1:tod%ndet)
    deallocate(coeff_matrix, rhs, A, b, x)

    do j = 1, tod%nscan
       do i = 1, tod%ndet
          tod%scans(j)%d(i)%gain = tod%gain0(0) + tod%gain0(i) + tod%scans(j)%d(i)%dgain 
       end do
    end do

  end subroutine sample_relcal



end module comm_tod_gain_mod
