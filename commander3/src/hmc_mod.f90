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
  use comm_utils
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Implementation of Hamiltonian Monte Carlo, specifically the No-U-Turn
  ! Sampler (NUTS) and standard HMC. The algorithms follow the implementation
  ! given in Hoffman & Gelman (2014), with additional changes referenced in
  ! Betancourt (2018).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  subroutine hmc(theta, lnlike, grad_lnlike, n_steps, eps, handle, length, M)
    ! Algorithm 5 (standard hmc plus dual averaging to adjust epsilon) from Hoffman & Gelman (2014) 
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

    integer(i4b) :: i, j, npar, L, t0, M_adapt, unit
    real(dp) :: alpha, H_m, logeps, logepsbar, gamm, kappa, mu, delta, acc_rate
    real(dp), dimension(size(theta)) :: p0, p, theta_prop, mass, theta_new, p_new

    unit = getlun()
    open(unit,file='angela_outputs/hmc/gaussian_tests/hmc_state.dat', recl=2048)
    write(unit, *) '#ITERATION  | POSITION_X  |  MOMENTUM_X'

    acc_rate = 0.0

    H_m = 0
    mu = log(10*eps)
    logepsbar = 0d0
    logeps = log(eps)
    gamm = 0.05
    t0 = 10
    kappa = 0.75
    delta = 0.65

    npar = size(theta)

    if (present(M)) then
      mass = M
    else
      mass = 1.d0
    end if

    if (present(length)) then
      L = length
    else
      L = max(1, int(1/eps))  !assumpion eps*L=1 this requires further tests
    end if

    theta_prop = theta

    M_adapt = n_steps/10

    do i = 1, n_steps + M_adapt
      do j = 1, npar
        p0(j) = rand_gauss(handle)*mass(j)
      end do
      theta_prop = theta
      p = p0

      do j = 1, L
        call Leapfrog(theta_prop, p, theta_new, p_new, eps, grad_lnlike, mass)
        theta_prop = theta_new
        p = p_new
        !write(*,*) '1', lnlike(theta_prop) - 0.5*sum(p**2/mass), theta_prop(1), p(1)
      end do
      alpha = min(1.d0, exp(lnlike(theta_prop) - lnlike(theta) &
                     &- 0.d5*(sum(p**2/mass) + 0.d5*sum(p0**2/mass)))  )
      if (alpha > rand_uni(handle)) then
        theta = theta_prop
        p = -p_new
        acc_rate = acc_rate + 1
      else
        theta = theta
        p = p0
      end if

      if (i .le. M_adapt) then
        H_m = (1d0 - 1d0/real(i+t0,dp))*H_m + (delta - alpha)/real(i+t0,dp)
        logeps = mu - sqrt(real(i,dp))/gamm*H_m
        logepsbar = i**(-kappa)*logeps + (1-i**(-kappa))*logepsbar
        eps = exp(logeps)
      !!  write(*,*) '2', theta(1), eps, exp(logepsbar)
      !!else
      !!  write(*,*) '2', theta(1)
      !end if
      !if (i .eq. M_adapt) then
      else
        eps = exp(logepsbar)
      end if

      if (mod(i, 1)==0) then
        write(unit, *)  i, '  | ', theta(1), '  | ', p(1)
      end if

    end do

    write(*,*) "final acceptance rate ", acc_rate/real(n_steps + M_adapt)

  end subroutine hmc

  subroutine nuts(theta, lnlike, grad_lnlike, n_steps, handle, eps, M)
    ! Algorithm 6 (efficient nuts as in alg 3 plus dual averaging) from Hoffman & Gelman (2014), efficient NUTS with slice sampler
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
    real(dp), optional,                           intent(inout) :: eps
    type(planck_rng),                intent(inout) :: handle
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


    integer(i4b) :: i, j, k,  npar, L, n, s, vel, n_p, s_p, n_pp, s_pp
    integer(i4b) :: n_alpha, t0, M_adapt
    real(dp) :: alpha, logu, H_m, logeps, logepsbar, gamm, kappa, mu, delta
    real(dp), dimension(size(theta)) :: theta_plus, theta_minus, theta_p
    real(dp), dimension(size(theta)) :: p_plus, p_minus, p_p, p, buff1, buff2, mass

    integer(i4b) :: unit

    npar = size(theta)

    if (.not. present(eps)) then
      eps = FindReasonableEpsilon(theta, lnlike, grad_lnlike, handle) 
    end if

    H_m = 0
    mu = log(10*eps)
    logepsbar = 0d0
    logeps = log(eps)
    gamm = 0.05
    t0 = 10
    kappa = 0.75
    delta = 0.65

    unit = getlun()
    open(unit,file='angela_outputs/hmc/gaussian_tests/two_hmc_state.dat', recl=2048)
    write(unit, *) '#ITERATION  | POSITION_X  | MOMENTUM_X  | POSITION_Y  | MOMENTUM_Y | LOGLIKE  |  GRADLOGLIKE'

    if (present(M)) then
      mass = M
    else
      mass = 1.d0
    end if

    M_adapt = n_steps/10

    do k = 1, n_steps + M_adapt
      do j = 1, npar
        p(j) = mass(j)*rand_gauss(handle)
      end do
      logu = lnlike(theta) - 0.5*sum(p**2/mass) + log(rand_uni(handle))
      theta_minus = theta
      theta_plus  = theta
      p_minus = p
      p_plus  = p

      j = 0
      n = 1
      s = 1

      do 
        if (rand_uni(handle) < 0.5) then
          vel = -1
          call BuildTree(theta_minus, p_minus, logu, vel, j, eps, theta, p, lnlike, grad_lnlike, mass, &
              & theta_minus, p_minus, buff1, buff2, theta_p, p_p, n_p, s_p, alpha, n_alpha, handle)
        else
          vel = 1
          call BuildTree(theta_plus,  p_plus,  logu, vel, j, eps, theta, p, lnlike, grad_lnlike, mass, &
              & buff1, buff2, theta_plus, p_plus,   theta_p, p_p, n_p, s_p, alpha, n_alpha, handle)
        end if
        if (s_p == 1) then
          if (rand_uni(handle) < min(1, n_p/n)) then
            theta = theta_p
            p = p_p
          end if
        end if
        n = n + n_p
        if ((dot_product(theta_plus - theta_minus, p_minus) < 0) .or. &
            (dot_product(theta_plus - theta_minus, p_plus) < 0)) then
            s = 0
        else
            s = s_p
        end if
        j = j + 1

        if (s .ne. 1) exit
      end do

      if (k .le. M_adapt) then
        H_m = (1d0 - 1d0/real(k+t0,dp))*H_m + (delta - alpha/n_alpha)/real(k+t0,dp)
        logeps = mu - sqrt(real(k,dp))/gamm*H_m
        logepsbar = k**(-kappa)*logeps + (1-k**(-kappa))*logepsbar
        eps = exp(logeps)
!     !   write(*,*) '4', theta(1), eps, exp(logepsbar)
!     ! else
!     !   write(*,*) '4', theta(1)
      !end if
      !if (k .eq. M_adapt) then
      else
        eps = exp(logepsbar)
      end if

      if ((mod(k, 10)==0) .or. (k < 10)) then
        write(unit, *)  k, '  | ', theta(1), '  | ', p(1), '  | ', theta(2), '  | ', p(2), '  | ', lnlike_hmc_test(theta), '  | ', grad_lnlike_hmc_test(theta)
      end if

    end do
  end subroutine nuts

  recursive subroutine BuildTree(theta, p, logu, v, j, eps, theta0, p0, lnlike, grad_lnlike, mass, &
                      & theta_minus, p_minus, theta_plus, p_plus, theta_p, p_p, n_p, s_p, &
                      & alpha, n_alpha, handle)
  ! BuildTree function as in algorithm 6 of Hoffman and Gelman
    implicit none
    real(dp), dimension(:),          intent(inout) :: theta, p, theta_minus, p_minus, theta_plus, p_plus, theta_p, p_p, mass, theta0, p0
    real(dp),                           intent(in) :: logu, eps
    integer(i4b),                       intent(in) :: v, j
    integer(i4b),                       intent(out) :: n_p, s_p
    integer(i4b),                       intent(out) :: n_alpha
    real(dp),                           intent(out) :: alpha
    type(planck_rng),                intent(inout) :: handle
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


    integer(i4b) :: s, s_pp, n_pp
    integer(i4b) :: n_ap, n_app
    real(dp) :: alpha_p, alpha_pp, E0, Ep
    real(dp), dimension(size(theta)) :: buff1, buff2, theta_pp, p_pp

    real(dp) :: deltamax = 1000

    if (j == 0) then
      call Leapfrog(theta, p, theta_p, p_p, v*eps, grad_lnlike, mass)
      !write(*,*) '3', lnlike(theta_p) - 0.5*sum(p_p**2/mass), theta_p(1), p_p(1)
      if (logu < lnlike(theta_p) - 0.5*sum(p_p**2/mass)) then
        n_p = 1
      else
        n_p = 0
      end if

      if (logu < deltamax + lnlike(theta_p) - 0.5*sum(p_p**2/mass)) then
        s_p = 1
      else
        s_p = 0
      end if

      theta_minus = theta_p
      theta_plus  = theta_p

      p_plus  = p_p
      p_minus = p_p

      n_alpha = 1
      E0 = lnlike(theta0)  - 0.5*sum(p0**2/mass)
      Ep = lnlike(theta_p) - 0.5*sum(p_p**2/mass)
      alpha = min(1d0, exp(Ep - E0))

    else
      call BuildTree(theta, p, logu, v, j-1, eps, theta0, p0, lnlike, grad_lnlike, mass, &
          & theta_minus, p_minus, theta_plus, p_plus, theta_p, p_p, n_p, s_p, alpha_p, n_ap, handle)

      if (s_p == 1) then
        if (v == -1) then
            call BuildTree(theta_minus, p_minus, logu, v, j-1, eps, theta0, p0, lnlike, grad_lnlike, mass, &
                & theta_minus, p_minus, buff1, buff2, theta_pp, p_pp, n_pp, s_pp, &
                & alpha_pp, n_app, handle)
        else
            call BuildTree(theta_plus, p_plus, logu, v, j-1, eps, theta0, p0, lnlike, grad_lnlike, mass, &
                & buff1, buff2, theta_plus,  p_plus,  theta_pp, p_pp, n_pp, s_pp, &
                & alpha_pp, n_app, handle)
        end if

        n_alpha = n_ap + n_app
        !alpha = log(exp(alpha_p) + exp(alpha_pp))
        alpha = alpha_p + alpha_pp
        !if (alpha_pp < 0) then
        !  write(*,*) alpha_pp, n_app, j, n_pp, s_pp, logu, theta_plus(1), p_plus(1), theta_minus(1), p_minus(1)
        !  stop
        !end if

        if (rand_uni(handle)*(n_p + n_pp) < n_pp) then
          theta_p = theta_pp
          p_p     = p_pp
        end if

        if ((dot_product(theta_plus - theta_minus, p_minus) < 0) .or. &
            (dot_product(theta_plus - theta_minus, p_plus) < 0)) then
            s_p = 0
        else
            s_p = s_p*s_pp
        end if

        n_p = n_p + n_pp
        !if ((n_ap < 0) .or. (n_app < 0)) then
        !    write(*,*) n_ap, n_app, alpha_p, alpha_pp
        !    stop
        !end if
      else
        n_alpha = n_ap
        alpha = alpha_p
      end if

    end if

  end subroutine BuildTree


  function lnlike_hmc_test(theta) result(lnlike)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in)  :: theta
    real(dp)                            :: lnlike

    lnlike = -(1.0_dp/2.0_dp)*(theta(1)**2+theta(2)**2)

  end function

  function grad_lnlike_hmc_test(theta) result(grad_lnlike)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in)  :: theta
    real(dp), dimension(size(theta))    :: grad_lnlike

    grad_lnlike(1) = -theta(1)
    grad_lnlike(2) = -theta(2)

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

    ! Leapfrog integrator, performs sympletic integration of Hamilton's equation

    if (minval(mass) .le. 0d0) then
      write(*,*) "Warning: Mass matrix in HMC has nonpositive values"
    end if

    p_new = p +     grad_func(x)    *eps/2
    x_new = x +     p_new/mass      *eps
    p_new = p_new + grad_func(x_new)*eps/2

  end subroutine Leapfrog


  function FindReasonableEpsilon(theta, lnlike, grad_lnlike, handle, M) result(final_eps)
    !
    ! If you have no idea what the timestep should be, this will give the value
    ! where the change in energy will be roughly 0.5 -- generally not a great
    ! estimate, but good for starting a tuning run.
    !
    implicit none
    real(dp), dimension(:),          intent(inout) :: theta
    type(planck_rng),                intent(inout) :: handle
    real(dp), optional, dimension(:),   intent(in) :: M    
    real(dp)                                       :: final_eps
    real(dp)                                       :: eps, npar, pp_over_p
    real(dp), dimension(size(theta))               :: p, p0, theta_prop, mass, theta_new, p_new
    integer(i4b)                                   :: i, a
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

    theta_new = theta
    p_new = p

    call Leapfrog(theta, p, theta_new, p_new, eps, grad_lnlike, mass)

    pp_over_p =  exp(lnlike(theta_new) - lnlike(theta) - 0.d5*(sum(p_new**2/mass) + 0.d5*sum(p**2/mass)))

    if (pp_over_p > 0.5) then
      a = 1
    else
      a = -1
    end if

    do 
      eps = eps*2.d0**a
      call Leapfrog(theta, p, theta_new, p_new, eps, grad_lnlike, mass)
      pp_over_p =  exp(lnlike(theta_new) - lnlike(theta) - 0.d5*(sum(p_new**2/mass) + 0.d5*sum(p**2/mass)))
      !write(*,*) eps, pp_over_p
      if (pp_over_p**a < 2.d0**(-a)) exit
    end do

    final_eps = eps

  end function FindReasonableEpsilon





end module hmc_mod
