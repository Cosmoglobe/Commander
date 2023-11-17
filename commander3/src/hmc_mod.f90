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

  subroutine hmc(theta, lnlike, grad_lnlike, n_steps, eps, handle, length, M)
    ! Algorithm 1 from Hoffman & Gelman (2014)
    ! Requires step size eps. n_steps is the number of independent draws,
    ! while "length" how many steps the Leapfrog integrator should be run. H&G
    ! recommend eps*length = 1, so it defaults to this unless specified
    ! otherwise. Could require some tuning to figure out the best epsilon value.

    !
    ! As implemented, this only returns the last of n_steps samples.
    !

    !
    ! Optional parameter M is the "mass matrix", which is theoretically the
    ! covariance matrix that the momenta are drawn from, but in practice this is
    ! diagonal. Gelman's Bayesian Data Analysis (3ed) suggests M could be the
    ! inverse covariance matrix of the parameters theta. "...better scaling of M
    ! will merely make HMC more efficient."
    !
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
    real(dp), allocatable, dimension(:) :: p0, p, theta_prop, mass, theta_new, p_new

    npar = size(theta)

    allocate(p0(npar), p(npar), theta_prop(npar), mass(npar), theta_new(npar), p_new(npar))

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
        call Leapfrog(theta_prop, p, eps, theta_new, p_new, grad_lnlike, mass)
        theta_prop = theta_new
        p = p_new
      end do
      alpha = min(1.d0, exp(lnlike(theta_prop) - lnlike(theta) &
                     &- 0.5*(sum(p**2/mass) - sum(p0**2/mass)))  )
      if (alpha > rand_uni(handle)) then
        theta = theta_prop
      end if
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

  subroutine Leapfrog(x, p, x_new, p_new, eps, grad_func, mass)
    implicit none
    real(dp), dimension(:), intent(inout) :: x, p
    real(dp), dimension(:), intent(out)   :: x_new, p_new
    real(dp),               intent(in)    :: eps
    real(dp), dimension(:), intent(in)    :: mass

    interface
      function grad_func(x)
        use healpix_types
        implicit none
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(size(x))       :: grad_func
      end function grad_func
    end interface

    ! Leapfrog integrator, performs sympletic integration of Hamilton's equation

    if (minval(mass) .le. 0d0) then
      write(*,*) "Warning: Mass matrix in HMC has nonpositive values"
    end if

    p_new = p + grad_func(x)*eps/2
    x_new = x + eps*p/mass
    p_new = p_new + grad_func(x_new)*eps/2

  end subroutine Leapfrog


  function FindReasonableEpsilon(theta, lnlike, grad_lnlike, handle, M)
       ! def FindReasonableEpsilon(theta, lnlike, grad_lnlike):
       !     eps = 1
       !     r = np.random.randn(theta.size)
       !     thetap, rp = Leapfrog(theta, r, eps, grad_lnlike)
       !     # p = exp(lnlike(theta) - 0.5*r**2)
       !     pp_over_p = np.exp(lnlike(thetap) - lnlike(theta)- 0.5*(rp.dot(rp) - r.dot(r)))
       !     if pp_over_p > 0.5:
       !         a = 1
       !     else:
       !         a = -1
       !     while pp_over_p**a > 2**-a:
       !         eps = 2**a*eps
       !         thetap, rp = Leapfrog(theta, r, eps, grad_lnlike)
       !         pp_over_p = np.exp(lnlike(thetap) - lnlike(theta)- 0.5*(rp.dot(rp) - r.dot(r)))
       !     return eps
    implicit none
    real(dp), dimension(:),          intent(inout) :: theta
    type(planck_rng),                intent(inout) :: handle
    real(dp), optional, dimension(:),   intent(in) :: M
    real(dp)                                       :: eps, a, npar, pp_over_p
    real(dp), dimension(size(theta))               :: p, p0, theta_prop, mass, theta_new, p_new
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

    npar = size(theta)

    if (present(M)) then
      mass = M
    else
      mass = 1.d0
    end if

    eps = 1d0
    do i = 1, npar
      p(i) = rand_gauss(handle)
    end do


    call Leapfrog(theta, p, theta_new, p_new, eps, grad_lnlike, mass)

    pp_over_p =  exp(lnlike(theta_new) - lnlike(theta) - 0.5*(sum(p_new**2/mass) - sum(p**2/mass)))

    if (pp_over_p > 0.5) then
      a = 1.d0
    else
      a = -1.d0
    end if

    do
      if ((pp_over_p)**a > 2**(-a)) exit
      eps = 2**a*eps
      call Leapfrog(theta, p, theta_new, p_new, eps, grad_lnlike, mass)
      pp_over_p =  exp(lnlike(theta_new) - lnlike(theta) - 0.5*(sum(p_new**2/mass) - sum(p**2/mass)))
    end do

    FindReasonableEpsilon = eps


  end function FindReasonableEpsilon





end module hmc_mod
