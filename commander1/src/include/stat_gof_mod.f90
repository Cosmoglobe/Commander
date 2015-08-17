module stat_gof_mod
  use healpix_types
  use sort_utils
  use math_tools
  use num_rec
  implicit none


contains

  subroutine compute_anderson_darling_norm_dist(X, A_sq)
    implicit none

    real(dp), dimension(1:), intent(in)  :: X
    real(dp),                intent(out) :: A_sq

    integer(i4b) :: i, j, n
    real(dp)     :: mu, sigma_sq, t1, t2
    real(dp), allocatable, dimension(:) :: Y

    n = size(X)

    ! Compute mean and variance
    mu       = sum(X) / real(n,dp)
    sigma_sq = sum((X-mu)**2) / real(n-1,dp)

    ! Sort data
    allocate(Y(n))
    Y = X
    call cpu_time(t1)
    call QuickSort_real(Y)    
    call cpu_time(t2)
!    write(*,*) 'sort = ', t2-t1

    ! Normalize data
    Y = (Y-mu) / sqrt(sigma_sq)

    ! Compute Anderson-Darling statistic
    A_sq = -real(n,dp)
    do i = 1, n
       A_sq = A_sq - 1.d0 / real(n,dp) * (2.d0*real(i,dp)-1) * (log(Phi_gauss(Y(i))) + log(1.d0 - Phi_gauss(Y(n+1-i))))
    end do
    call cpu_time(t1)
!    write(*,*) 'comp = ', t1-t2
    
    ! Correct for bias from estimating mu and sigma
    A_sq = A_sq * (1.d0 + 4.d0 / real(n,dp) - 25.d0 / real(n*n,dp))

    deallocate(Y)

!!$!    if (A_sq > 10.3d0) then
!!$       open(68,file='X.dat')
!!$       do i = 1, n
!!$          write(68,*) X(i)
!!$       end do
!!$       close(68)
!!$       stop
!!$!    end if

  end subroutine compute_anderson_darling_norm_dist


  subroutine check_gof_gaussian(mu, sigma_sq, x_low, x_high, samples, p)
    implicit none

    real(dp),                intent(in)  :: mu, sigma_sq, x_low, x_high
    real(dp), dimension(1:), intent(in)  :: samples
    real(dp),                intent(out) :: p

    integer(i4b) :: i, j, bin, nbin, numsamp, n
    real(dp)     :: dx, x, chisq, C, gln
    real(dp), allocatable, dimension(:) :: hist, Y, bins, sorted

    numsamp = size(samples)
    nbin    = nint(sqrt(real(numsamp,dp)))
    dx      = (x_high-x_low) / real(nbin,dp)

    ! Find number of valid samples
    n = count(samples > x_low .and. samples < x_high)
!    do i = 1, numsamp
!       if (samples(i) > x_low .and. samples(i) < x_high) n = n+1
!    end do

    ! Sort samples
    allocate(sorted(n))
    j = 1
    do i = 1, numsamp
       if (samples(i) > x_low .and. samples(i) < x_high) then
          sorted(j) = samples(i)
          j         = j+1
       end if
    end do
    call QuickSort_real(sorted)

    ! Set up bins
    allocate(bins(0:nbin))
    call get_bins(sorted, nbin, bins)

    ! Compute expected distribution
    allocate(Y(nbin))
    do i = 1, nbin
       x    = 0.5d0 * (bins(i-1) + bins(i))
       Y(i) = exp(-0.5d0 * (x-mu)**2 / sigma_sq) * (bins(i)-bins(i-1))
    end do

    ! Normalize distribution
    Y = Y / sum(Y) * real(n,dp)


    do while (minval(Y) < 5.d0)
       nbin = nbin/2
       deallocate(bins)
       deallocate(Y)
       allocate(bins(0:nbin))
       allocate(Y(nbin))
       call get_bins(sorted, nbin, bins)
       do i = 1, nbin
          x    = 0.5d0 * (bins(i-1) + bins(i))
          Y(i) = exp(-0.5d0 * (x-mu)**2 / sigma_sq) * (bins(i)-bins(i-1))
       end do
       Y = Y / sum(Y) * real(n,dp)
    end do


    ! Compute histogram
    allocate(hist(nbin))
    hist = real(n,dp) / real(nbin,dp)
!    do i = 1, nbin
!       hist(i) = hist(i) / (bins(i)-bins(i-1))
!    end do
!    hist = 0.d0
!    do i = 1, numsamp
!       if (samples(i) > x_low .and. samples(i) < x_high) then
!          bin       = locate_bin(bins, samples(i))
!          hist(bin) = hist(bin) + 1.d0
!       end if
!    end do
    
    ! Compute C statistic
    C = 0.d0
    do i = 1, nbin
       C = C + (hist(i) - Y(i))**2 / Y(i)
    end do

!    write(*,*) 'nbin = ', nbin
!    write(*,*) bins
!    open(58,file='gof.dat')
!    do i = 1, nbin
!       write(58,*) 0.5d0*(bins(i-1)+bins(i)), hist(i) / (bins(i)-bins(i-1))
!    end do
!    write(58,*)
!    do i = 1, nbin
!       write(58,*) 0.5d0*(bins(i-1)+bins(i)), Y(i) / (bins(i)-bins(i-1))
!    end do
!    close(58)

!    open(58,file='samples.dat')
!    do i = 1, n
!       write(58,*) samples(i)
!    end do
!    close(58)
!    stop

    if (nbin < 4) then
!       write(*,*) 'Error: Too few samples to compute Gaussian goodness-of-fit'
!       stop
       p = 1.d0
    else
       ! Compute probability
       write(*,*) 'ERROR!!'
!       p = gammp(0.5d0*real(nbin-3,dp), 0.5d0*C)
    end if
    
!    write(*,*) nbin, C, p
!    write(*,*) 0.5d0*real(nbin-3,dp),0.5d0*C, p
!    if (p > 0.9d0 .and. p < 0.99d0) stop
    
    deallocate(hist)
    deallocate(Y)
    if (allocated(bins))   deallocate(bins)
    if (allocated(sorted)) deallocate(sorted)

  end subroutine check_gof_gaussian

  subroutine get_bins(x, nbin, bins)
    implicit none

    real(dp),     dimension(1:), intent(in)  :: x
    integer(i4b),                intent(in)  :: nbin
    real(dp),     dimension(0:), intent(out) :: bins

    integer(i4b) :: i, n

    n = size(x)

    bins(0) = x(1)
    do i = 1, nbin-1
       bins(i) = x(nint(n*real(i,dp)/real(nbin,dp)))
    end do
    bins(nbin) = x(n)

  end subroutine get_bins


  subroutine compute_kolmogorov_smirnow_2side(X, Y, K)
    implicit none

    real(dp), dimension(1:), intent(in)  :: X, Y
    real(dp),                intent(out) :: K

    integer(i4b) :: n, np, i, j
    real(dp)     :: F_X, F_Y, dX, dY, X_arg, Y_arg, D
    real(dp), allocatable, dimension(:) :: Xs, Ys

    n  = size(X)
    np = size(Y)
    dX = 1.d0 / real(n,dp)
    dY = 1.d0 / real(np,dp)

    allocate(Xs(n))
    allocate(Ys(np))
    Xs = X
    Ys = Y
    call QuickSort_real(Xs)
    call QuickSort_real(Ys)

    ! Compute the KS statistic, without binning the samples
    i     = n
    j     = np
    X_arg = Xs(i)
    Y_arg = Ys(j)
    F_x   = 1.d0
    F_y   = 1.d0
    D     = 0.d0

    do while (i > 1 .and. j > 1)
       
    end do



    deallocate(Xs)
    deallocate(Ys)

  end subroutine compute_kolmogorov_smirnow_2side


  function Phi_gauss(X)
    implicit none

    real(dp), intent(in) :: X
    real(dp)             :: Phi_gauss

    if (X < -6.d0) then
       Phi_gauss = 1.d-10
    else if (X > 6.d0) then
       Phi_gauss = 1.d0 - 1.d-10
    else
       Phi_gauss = 0.5d0 * (1.d0 + corr_erf(X/sqrt(2.d0)))
    end if

  end function Phi_gauss


  function locate_bin(xx, x)
    implicit none

    real(dp),               intent(in) :: x
    real(dp), dimension(:), intent(in) :: xx
    integer(i4b)                       :: locate_bin

    integer(i4b) :: n, jl, jm, ju
    logical(lgt) :: ascnd

    n     = size(xx)
    ascnd = (xx(n) >= xx(1))
    jl    = 0
    ju    = n+1

    do
       if (ju-jl <= 1) exit
       jm = (ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       end if
    end do

    if (x == xx(1)) then
       locate_bin = 1
    else if (x == xx(n)) then
       locate_bin = n-1
    else
       locate_bin = jl
    end if

  end function locate_bin


end module stat_gof_mod
