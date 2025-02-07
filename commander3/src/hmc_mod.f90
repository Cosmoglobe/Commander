! Example run
!  ! HMC testing
!  real(dp), dimension(5) :: param_test
!  real(dp) :: time_step
!  param_test = 1.d0
!  time_step = 1d-1
!
!  if (cpar%myid == cpar%root) then
!      write(*,*) "first", param_test(1)
!      call hmc(param_test, lnlike_hmc_test, grad_lnlike_hmc_test, 10000, time_step, handle)
!      call nuts(param_test, lnlike_hmc_test, grad_lnlike_hmc_test, 10000, time_step, handle)
!      write(*,*) "last", param_test(1)
!  end if
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
    real(dp),                        intent(inout) :: eps
    type(planck_rng),                intent(inout) :: handle
    integer(i4b), optional,             intent(in) :: length
    real(dp), optional, dimension(:),   intent(in) :: M


  end subroutine hmc

  subroutine nuts(theta, lnlike, grad_lnlike, n_steps, eps, handle, M)
    ! Algorithm 3 from Hoffman & Gelman (2014), efficient NUTS with slice sampler
    ! Does not take length as an argument, as it samples until the U-Turn
    ! condition is reached (momentum and (theta_plus - theta_minus) are
    ! perpendicular).

    ! Takes theta, lnlike, grad_lnlike, n_steps, eps, adn handle as arguments,
    ! with optional mass matrix M.
    ! After performing n_steps samples, the subroutine will return the latest
    ! sample. 
    implicit none
    real(dp), dimension(:),          intent(inout) :: theta
    integer(i4b),                       intent(in) :: n_steps
    real(dp),                           intent(inout) :: eps
    type(planck_rng),                intent(inout) :: handle
    real(dp), optional, dimension(:),   intent(in) :: M


  end subroutine nuts

  recursive subroutine BuildTree(theta, p, logu, v, j, eps, theta0, p0, lnlike, grad_lnlike, mass, &
                      & theta_minus, p_minus, theta_plus, p_plus, theta_p, p_p, n_p, s_p, &
                      & alpha, n_alpha, handle)
    implicit none
    real(dp), dimension(:),          intent(inout) :: theta, p, theta_minus, p_minus, theta_plus, p_plus, theta_p, p_p, mass, theta0, p0
    real(dp),                           intent(in) :: logu, eps
    integer(i4b),                       intent(in) :: v, j
    integer(i4b),                       intent(out) :: n_p, s_p
    integer(i4b),                       intent(out) :: n_alpha
    real(dp),                           intent(out) :: alpha
    type(planck_rng),                intent(inout) :: handle


  end subroutine BuildTree


  function lnlike_hmc_test(theta)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in)  :: theta
    real(dp)                            :: lnlike_hmc_test

    lnlike_hmc_test = -sum(theta**2)/2

  end function

  function grad_lnlike_hmc_test(theta)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in)  :: theta
    real(dp), dimension(size(theta))    :: grad_lnlike_hmc_test

    grad_lnlike_hmc_test = -theta

  end function



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


  end subroutine Leapfrog


  function FindReasonableEpsilon(theta, lnlike, grad_lnlike, handle, M)
    !
    ! If you have no idea what the timestep should be, this will give the value
    ! where the change in energy will be roughly 0.5 -- generally not a great
    ! estimate, but good for starting a tuning run.
    !
    implicit none
    real(dp), dimension(:),          intent(inout) :: theta
    type(planck_rng),                intent(inout) :: handle
    real(dp), optional, dimension(:),   intent(in) :: M
    real(dp)                                       :: FindReasonableEpsilon
    real(dp)                                       :: eps, npar, pp_over_p
    real(dp), dimension(size(theta))               :: p, p0, theta_prop, mass, theta_new, p_new


    FindReasonableEpsilon = 0d0


  end function FindReasonableEpsilon





end module hmc_mod
