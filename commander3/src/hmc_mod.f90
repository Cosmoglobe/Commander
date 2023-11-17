module hmc_mod
  use healpix_types
  use rngmod
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Implementation of Hamiltonian Monte Carlo, specifically the No-U-Turn
  ! Sampler (NUTS) and standard HMC. The algorithms follow the implementation
  ! given in Hoffman & Gelman (2014), with additional changes referenced in
  ! Betancourt (2018).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  subroutine Leapfrog(x, p, eps, grad_func)
    implicit none
    real(dp), dimension(:), intent(inout) :: x, p
    real(dp),               intent(in)    :: eps

    interface
      function grad_func(x)
        use healpix_types
        implicit none
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(size(x))       :: grad_func
      end function grad_func
    end interface

    ! Leapfrog integrator, performs sympletic integration of Hamilton's equation

    p = p + grad_func(x)*eps/2
    x = x + eps*p
    p = p + grad_func(x)*eps/2

  end subroutine Leapfrog

  subroutine hmc(theta, lnlike, grad_lnlike, n_steps, eps, handle, length, M)
    implicit none
    real(dp), dimension(:),          intent(inout) :: theta
    integer(i4b),                       intent(in) :: n_steps
    real(dp),                           intent(in) :: eps
    type(planck_rng),                intent(inout) :: handle
    integer(i4b), optional,             intent(in) :: length
    real(dp), optional, dimension(:),   intent(in) :: M
    interface
       function lnlike(theta)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in) :: theta
         real(dp)                           :: lnlike
       end function lnlike

       function grad_lnlike(theta)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in) :: theta
         real(dp), dimension(size(theta))   :: grad_lnlike
       end function grad_lnlike
    end interface


    integer(i4b) :: i, j, npar, L
    real(dp) :: alpha
    real(dp), allocatable, dimension(:) :: p0, p, theta_prop, mass

    npar = size(theta)

    allocate(p0(npar), p(npar), theta_prop(npar), mass(npar))

    if (present(M)) then
      mass = M
    else
      mass = 1.d0
    end if
    if (present(length)) then
      L = length
    else
      L = int(1/eps)
    end if

    theta_prop = theta

    do i = 1, n_steps
      do j = 1, npar
        p0(j) = rand_gauss(handle)*mass(j)
      end do
      theta_prop = theta
      p = p0

      do j = 1, L
        call Leapfrog(theta_prop, p, eps, grad_lnlike)
      end do
      alpha = min(1.d0, exp(lnlike(theta_prop) - lnlike(theta) &
                     &- 0.5*(sum(p**2/mass) - sum(p0**2/mass)))  )
      if (alpha > rand_uni(handle)) then
        theta = theta_prop
      end if
      write(*,*) theta(1)
    end do

    deallocate(p0, p, theta_prop, mass)


  end subroutine hmc


  function lnlike_hmc_test(theta)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in)  :: theta
    real(dp)                            :: lnlike_hmc_test
    !
    ! Using a multivariate standard normal Gaussian, and I want to estimate the means.
    !

    lnlike_hmc_test = -sum(theta**2)/2

  end function

  function grad_lnlike_hmc_test(theta)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in)  :: theta
    real(dp), dimension(size(theta))    :: grad_lnlike_hmc_test

    grad_lnlike_hmc_test = -theta

  end function


  ! Example run
  !  ! HMC testing
  !  real(dp), dimension(5) :: param_test
  !  real(dp) :: time_step

  !  param_test = 5.
  !  time_step = 0.01

  !  if (cpar%myid == cpar%root) then
  !      write(*,*) "first", param_test(1)
  !      call hmc(param_test, lnlike_hmc_test, grad_lnlike_hmc_test, 100, time_step, handle)
  !      write(*,*) "last", param_test(1)
  !  end if



end module hmc_mod
