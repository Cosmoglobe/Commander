module ARS_mod
  use healpix_types
  use rngmod
  use sort_utils
  implicit none

  integer(i4b), parameter, private :: MAX_ARS_TRIES       = 100
  integer(i4b), parameter, private :: MAX_ARS_BOUND_TRIES = 10
  real(dp),     parameter, private :: EXP_OVERFLOW_LIMIT = 20.d0

contains

  function sample_ARS(handle, lnL, x_init, prior, status, neval, x_n_out, S_n_out, bounds_ok)
    implicit none

    type(planck_rng)                                           :: handle
    real(dp),              dimension(:), intent(in)            :: x_init
    real(dp)                                                   :: sample_ARS
    real(dp),              dimension(2),              optional :: prior
    integer(i4b),                                     optional :: neval, status
    real(dp), allocatable, dimension(:), intent(out), optional :: x_n_out, S_n_out
    logical(lgt),                        intent(in),  optional :: bounds_ok
    interface
       function lnL(x)
         use healpix_types
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: lnL
       end function lnL
    end interface

    integer(i4b) :: i, j, k, n, iter, stat
    real(dp)     :: p(2), x_new, y_new, h, u
    real(dp), dimension(MAX_ARS_TRIES) :: x_n
    real(dp), dimension(MAX_ARS_TRIES) :: S_n
    real(dp), allocatable, dimension(:) :: xv, yv

    if (present(prior)) then
       p = prior
    else
       p = [-1.d100, 1.d100]
    end if

    ! Initialize search
    stat = 0
    n    = 0
    do i = 1, size(x_init)
       if (x_init(i) >= p(1) .and. x_init(i) <= p(2)) then
          n      = n+1
          x_n(n) = x_init(n)
          S_n(n) = lnL(x_n(n))
       end if
    end do

    ! Check if distribution is uniform; if so, sample directly from a uniform distribution
    if (all(S_n(1:n) == S_n(1))) then
       sample_ARS = prior(1) + (prior(2)-prior(1)) * rand_uni(handle)
       return
    end if

    ! Check that peak is bounded; if not do a golden ratio search
    do while (S_n(2) < S_n(1) .and. x_n(1) > p(1))
       x_new = max(x_n(1) - 1.61803d0*(x_n(2)-x_n(1)), p(1))
       y_new = lnL(x_new)
       call update_sample_set(x_new, y_new, x_n, S_n, n, stat)
    end do

    do while (S_n(n) > S_n(n-1) .and. x_n(n) < p(2))
       x_new = min(x_n(n) + 1.61803d0*(x_n(n)-x_n(n-1)), p(2))
       y_new = lnL(x_new)
       call update_sample_set(x_new, y_new, x_n, S_n, n, stat)
    end do
    if (present(bounds_ok)) then
       if (stat /= 0) then
          stat = stat+1
          goto 99
       end if
    else
       if (stat /= 0 .or. x_n(1) == p(1) .or. x_n(n) == p(2)) then
          stat = stat+1
          goto 99
       end if
    end if

    call QuickSort_dp_dist(S_n(1:n), x_n(1:n))
    call get_vertices(x_n, S_n, n, xv, yv, p, stat)
    if (stat /= 0) goto 99

    ! Loop until a sample has been drawn successfully
    k = 0
    do while (k < MAX_ARS_TRIES)
       k = k+1

       ! Draw from piecewise exponential envelope
       x_new = sample_piecewise_exponential(handle, xv, yv, n, stat)
       if (stat /= 0) goto 99

       ! Draw uniform random number
       u = rand_uni(handle)
       
       ! Find likelihood and envelope value at x_new
       y_new = lnL(x_new)
       h     = envelope_function(x_new, xv, yv, n, stat)
       if (stat /= 0) goto 99

       if (u > exp(y_new - h)) then
          ! Reject; adapt the envelope
          call update_sample_set(x_new, y_new, x_n, S_n, n, stat)
          if (stat /= 0) goto 99

          call get_vertices(x_n, S_n, n, xv, yv, p, stat)
          if (stat /= 0) goto 99
       else
          ! Accept; return sample
          sample_ARS = x_new
          exit
       end if

    end do

    if (k == MAX_ARS_TRIES) then
       write(*,*) 'ARS_mod -- Too many attempts without successful draw'
       sample_ARS = 1.d30
       stat       = stat+10
    end if

99  if (allocated(xv)) deallocate(xv)
    if (allocated(yv)) deallocate(yv)

    if (stat /= 0) then
       sample_ARS = 1.d30
    else
       sample_ARS = max(min(sample_ARS, p(2)), p(1))
    end if
    if (present(status)) status = stat
    if (present(neval)) neval = n
    if (present(x_n_out) .and. n > 0) then
       allocate(x_n_out(n), S_n_out(n))
       x_n_out = x_n(1:n)
       S_n_out = S_n(1:n)
    end if
    if (present(status)) then
       status = stat
    else
       if (stat /= 0) then
          write(*,*) 'ARS sampler failed; stat = ', stat
          stop
       end if
    end if

  end function sample_ARS

  subroutine get_vertices(x, y, n, xv, yv, xbounds, stat)
    implicit none

    real(dp),                  dimension(:), intent(in)    :: x, y
    integer(i4b),                            intent(in)    :: n
    real(dp),     allocatable, dimension(:), intent(out)   :: xv, yv
    real(dp),                  dimension(2), intent(in)    :: xbounds
    integer(i4b),                            intent(inout) :: stat

    integer(i4b) :: i, m
    real(dp)     :: xrange(2), a0, a1, b0, b1

!    write(*,*) 'x = ', real(x(1:n),sp)

    if (n < 3) then
       write(*,*) 'ARS_mod -- n must be larger than 3'
       stat = stat+1
       return
       !stop
    end if

    ! Allocate internal data structures
    m      = 2*n+1
    if (allocated(xv)) deallocate(xv); if (allocated(yv)) deallocate(yv)
    allocate(xv(m), yv(m))
    xv(1) = xbounds(1); yv(1) = y(1) + (y(2)-y(1))*(xbounds(1)-x(1))/(x(2)-x(1))
    xv(2) = x(1); yv(2) = y(1)
    xv(3) = x(1); yv(3) = y(2) + (y(3)-y(2))*(x(1)-x(2))/(x(3)-x(2))
    xv(4) = x(2); yv(4) = y(2)

!    write(*,*) xv(1:4)
!    write(*,*) yv(1:4)
!    stop

    do i = 2, n-2
       b0 = y(i) - y(i-1); b1 = y(i+2) - y(i+1)
       a0 = x(i) - x(i-1); a1 = x(i+2) - x(i+1)
       xv(2*i+1) = (b0*a1*x(i-1) - b1*a0*x(i+1) + a0*a1*(y(i+1)-y(i-1))) / &
            & (b0*a1 - b1*a0)
       yv(2*i+1) = y(i-1) + b0*(xv(2*i+1)-x(i-1))/a0
       xv(2*i+2) = x(i+1)
       yv(2*i+2) = y(i+1)
    end do

    xv(m-2) = x(n); yv(m-2) = y(n-2) + (y(n-1)-y(n-2))*(x(n)-x(n-2))/(x(n-1)-x(n-2))
    xv(m-1) = x(n); yv(m-1) = y(n)
    xv(m)   = xbounds(2); yv(m) = y(n-1) + (y(n)-y(n-1))*(xbounds(2)-x(n-1))/(x(n)-x(n-1))

    ! Check that envelope is consistent
    do i = 1, m-1
       if (xv(i+1) < xv(i)) then
          stat = stat + 1
!!$          write(*,*) 'x_n = ', x(1:n)
!!$          write(*,*) 'S_n = ', y(1:n)
!!$          write(*,*) 'error in get_vertices', real(xv,sp)
          !stop
       end if
    end do

!    write(*,*) 'inni', real(xv,sp)
!    stop
  end subroutine get_vertices

  subroutine update_sample_set(x_new, y_new, x, y, n, stat)
    implicit none

    real(dp),                   intent(in)    :: x_new, y_new
    real(dp),     dimension(:), intent(inout) :: x, y
    integer(i4b),               intent(inout) :: n, stat

    integer(i4b) :: i

    if (n == MAX_ARS_TRIES) then
       !write(*,*) 'ARS_mod -- sample set full. Increase MAX_ARS_TRIES.'
       !stop
       stat = stat+1
       write(*,*) 'error max_arx_tries'
       return
    end if

    do i = 1, n
       if (x_new < x(i)) exit
    end do
    x(i+1:n+1) = x(i:n); x(i) = x_new
    y(i+1:n+1) = y(i:n); y(i) = y_new
    n      = n+1

  end subroutine update_sample_set

  function envelope_function(x, xv, yv, n, stat)
    implicit none

    real(dp),                    intent(in) :: x
    real(dp),     dimension(1:), intent(in) :: xv, yv
    integer(i4b),                intent(in) :: n
    real(dp)                                :: envelope_function
    integer(i4b),                intent(inout) :: stat

    real(dp)     :: dx, this_dx, s, c
    integer(i4b) :: i, j, m

    m       = 2*n+1
    if (x < xv(1) .or. x > xv(m)) then
       write(*,*) 'ARS_mod -- x out of prior range; x = ', real(x,sp)
       stat = stat+1
       return
       !stop
    end if

    ! Find which piece of the envelope function x corresponds to
    j       = m
    dx      = 2.d0 * (xv(m)-xv(1))
    this_dx = 0.d0
    do i = 1, m
       this_dx = xv(i) - x
       if (this_dx > 0.d0 .and. this_dx < dx) then
          j  = i
          dx = this_dx
       end if
    end do
    
    ! Compute slope and offset
    s = (yv(j) - yv(j-1)) / (xv(j) - xv(j-1))
    c = yv(j) - s*xv(j)
    
!    write(*,*) j
!    write(*,*) 's ', s, c
!    write(*,*) 'xv = ', xv(j-1), xv(j)
!    write(*,*) 'yv = ', yv(j-1), yv(j)

    envelope_function = s*x + c

  end function envelope_function

  function sample_piecewise_exponential(handle, xv, yv, n, stat)
    implicit none

    type(planck_rng)                      :: handle
    real(dp), dimension(:), intent(in)    :: xv, yv
    integer(i4b),           intent(in)    :: n
    real(dp)                              :: sample_piecewise_exponential
    integer(i4b),           intent(inout) :: stat

    real(dp)     :: y0, u, du, this_du, f, x
    integer(i4b) :: i, j, m
    real(dp), allocatable, dimension(:) :: s, cdf

    m = 2*n+1

    ! Rescale yv to avoid overruns
    y0 = maxval(yv)

    ! Calculate cdf and gradient between vertices
    allocate(s(m-1), cdf(m))
    cdf(1) = 0.d0
    do i = 1, m-1
       if (xv(i+1) == xv(i)) then
          cdf(i+1) = cdf(i)
       else
          s(i) = (yv(i+1)-yv(i)) / (xv(i+1)-xv(i))
          cdf(i+1) = cdf(i) + (exp(yv(i+1)-y0) - exp(yv(i)-y0)) / s(i)
       end if
    end do

    ! Test CDF for validity (infinite or too big), and positive semi-finiteness
    if (any(cdf > 1.d100)) then
       write(*,*) 'ARS_mod -- CDF is invalid; integrates to infinity'
       stat = stat+1
       deallocate(s, cdf)
       return
       !write(*,*) 'ARS_mod -- CDF is invalid; integrates to infinity'
       !stop
    end if
    if (any(cdf < 0.d0)) then
!!$       write(*,*) 'ARS_mod -- CDF is invalid; negative values found'
!!$       write(*,*) 'xv', real(xv,sp)
!!$       write(*,*) 'yv', real(yv,sp)
!!$       write(*,*) 'cdf', real(cdf,sp)
!!$       write(*,*) 's', real(s,sp)
!!$       do i = 1, m-1
!!$          if (xv(i+1) == xv(i)) then
!!$             cdf(i+1) = cdf(i)
!!$          else
!!$             s(i) = (yv(i+1)-yv(i)) / (xv(i+1)-xv(i))
!!$             cdf(i+1) = cdf(i) + (exp(yv(i+1)-y0) - exp(yv(i)-y0)) / s(i)
!!$             if ((exp(yv(i+1)-y0) - exp(yv(i)-y0)) / s(i) < 0.d0) then
!!$                write(*,*) s(i)
!!$                write(*,*) exp(yv(i+1)-y0), exp(yv(i)-y0), exp(yv(i+1)-y0)-exp(yv(i)-y0)
!!$             end if
!!$          end if
!!$       end do

       stat = stat+1
       deallocate(s, cdf)
       return
    end if

    ! Draw a uniform random number
    u = rand_uni(handle)

    ! Find which piece of the envelope it corresponds to
    j = n
    du = 2.d0*cdf(m)
    do i = 1, m
       this_du = u*cdf(m) - cdf(i)
       if (this_du > 0.d0 .and. this_du <= du) then
          j  = i
          du = this_du
       end if
    end do

    ! Get sample of using analytic equation for this piece
    if (-(yv(j)-y0) >  EXP_OVERFLOW_LIMIT) then
       x = xv(j) + (log((u*cdf(m) - cdf(j)) * s(j)) - (yv(j)-y0)) / s(j)
    else
       f = (u*cdf(m) - cdf(j)) * s(j) * exp(-(yv(j) - y0))
       if (f < -1.d0) then
          write(*,*) 'ARS_mod -- error: f < -1.d0'
          stat = stat+1
          deallocate(s, cdf)
          return
          !write(*,*) 'ARS_mod -- error: f < -1.d0'
          !stop
       end if
       x = xv(j) + log(1.d0 + f) / s(j)
    end if

!    write(*,*) real(xv,sp)
!    write(*,*) real(cdf,sp)
!    stop

    sample_piecewise_exponential = x

    deallocate(s, cdf)

  end function sample_piecewise_exponential

end module ARS_mod
