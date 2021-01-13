module interpolation_mod
  use healpix_types
  implicit none


contains

  ! These routines are taken from Numerical Recipes. 

  subroutine spline(x, y, yp1, ypn, y2)
    implicit none

    real(dp),               intent(in)  :: yp1, ypn
    real(dp), dimension(:), intent(in)  :: x, y
    real(dp), dimension(:), intent(out) :: y2

    integer(i4b) :: n
    real(dp), dimension(size(x)) :: a, b, c, r

    n = size(x)

    c(1:n-1) = x(2:n) - x(1:n-1)
    r(1:n-1) = 6.0 * ((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1) = r(2:n-1) - r(1:n-2)
    a(2:n-1) = c(1:n-2)
    b(2:n-1) = 2.0 * (c(2:n-1)+a(2:n-1))
    b(1)     = 1.0
    b(n)     = 1.0

    if (yp1 > 0.99e30) then
       r(1) = 0.0
       c(1) = 0.0
    else
       r(1) = (3.0 / (x(2)-x(1))) * ((y(2)-y(1))/(x(2)-x(1)) - yp1)
       c(1) = 0.5
    end if

    if (ypn > 0.99e30) then
       r(n) = 0.
       a(n) = 0.
    else
       r(n) = (-3.0 / (x(n)-x(n-1))) * ((y(n)-y(n-1))/(x(n)-x(n-1)) - ypn)
       a(n) = 0.5
    end if

    call tridag(a(2:n), b(1:n), c(1:n-1), r(1:n), y2(1:n))

  end subroutine spline

  function splint(xa, ya, y2a, x)
    implicit none

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: xa, ya, y2a
    real(dp)                            :: splint

    integer(i4b) :: khi, klo, n
    real(dp)     :: a, b, h

    n = size(xa)

    klo = max(min(locate(xa,x),n-1),1)
    khi = klo+1
    h   = xa(khi) - xa(klo)
    a   = (xa(khi) - x) / h
    b   = (x - xa(klo)) / h
    
    splint = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi))*(h**2)/6.0

  end function splint


  subroutine splie2(x1a, x2a, d1, d2, ya, y2a)
    implicit none

    real(dp), dimension(:),   intent(in)  :: x1a, x2a, d1, d2
    real(dp), dimension(:,:), intent(in)  :: ya
    real(dp), dimension(:,:), intent(out) :: y2a

    integer(i4b) :: j, m, ndum

    m    = size(x1a)
    ndum = size(x2a)

    do j = 1, ndum
       call spline(x1a, ya(:,j), d1(j), d2(j), y2a(:,j))
!       call spline(x2a, ya(j,:), 1.0e30, 1.0e30, y2a(j,:))
    end do

  end subroutine splie2

  function splin2(x1a, x2a, ya, y2a, x1,x2)
    implicit none

    real(dp),                 intent(in) :: x1, x2
    real(dp), dimension(:),   intent(in) :: x1a, x2a
    real(dp), dimension(:,:), intent(in) :: ya, y2a
    real(dp)                             :: splin2

    integer(i4b) :: j, m, ndum
    real(dp), dimension(size(x2a)) :: yytmp, y2tmp2

    m    = size(x1a)
    ndum = size(x2a)

    do j = 1, ndum
       yytmp(j) = splint(x1a, ya(:,j), y2a(:,j), x1)
    end do

!    call spline(x2a, yytmp, 0.d0, 0., y2tmp2)
    call spline(x1a, yytmp, 1.d30, 1.d30, y2tmp2)
    
    splin2 = splint(x2a, yytmp, y2tmp2, x2)

  end function splin2

  
  function locate(xx, x)
    implicit none

    real(dp),               intent(in) :: x
    real(dp), dimension(:), intent(in) :: xx
    integer(i4b)                       :: locate

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
       locate = 1
    else if (x == xx(n)) then
       locate = n-1
    else
       locate = jl
    end if

  end function locate

  subroutine tridag(a, b, c, r, u)
    implicit none

    real(dp), dimension(:), intent(in)  :: a, b, c, r
    real(dp), dimension(:), intent(out) :: u

    real(dp)     :: bet
    integer(i4b) :: n, j
    real(dp), dimension(size(b)) :: gam
    
    n   = size(b)
    bet = b(1)

    u(1) = r(1) / bet
    do j = 2, n
       gam(j) = c(j-1)/bet
       bet    = b(j)-a(j-1)*gam(j)
       u(j)   = (r(j)-a(j-1)*u(j-1))/bet
    end do

    do j = n-1, 1, -1
       u(j) = u(j) - gam(j+1)*u(j+1)
    end do

  end subroutine tridag

end module interpolation_mod
