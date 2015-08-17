module spline_1D_mod
  use healpix_types
  use locate_mod
  use math_tools
  implicit none

  type spline_type
     real(dp), dimension(:), allocatable :: x, y, y2
     real(dp)                            :: boundary(2)
     logical(lgt)                        :: regular
     logical(lgt)                        :: linear
  end type

  interface spline
     module procedure spline_simple, spline_plain
  end interface

  interface splint
     module procedure splint_simple, splint_plain
  end interface

  interface splint_multi
     module procedure splint_simple_multi
  end interface

contains

  ! Simple spline interface routines. Uses the spline
  ! type which makes internal copies, so uses slightly
  ! more memory than necessary, but probably not a big
  ! deal. And saves you from having to pass in 3 arguments
  ! all the time.
  subroutine spline_simple(s, x, y, boundary, regular, linear)
    implicit none
    type(spline_type),       intent(inout) :: s
    real(dp),                intent(in)    :: x(:), y(:)
    real(dp),      optional, intent(in)    :: boundary(2)
    logical(lgt),  optional, intent(in)    :: regular, linear
    call free_spline(s)
    s%boundary = 1d30
    s%regular  = .false.
    s%linear   = .false.
    if(present(boundary)) s%boundary = boundary
    if(present(regular))  s%regular  = regular
    if(present(linear))   s%linear   = linear
    allocate(s%x(size(x)),s%y(size(x)),s%y2(size(x)))
    s%x = x
    s%y = y
    if(.not. s%linear) call spline(s%x, s%y, s%boundary(1), s%boundary(2), s%y2)
  end subroutine

  function splint_simple(s, x) result(y)
    implicit none
    type(spline_type), intent(in) :: s
    real(dp),          intent(in) :: x
    real(dp)                      :: y
    integer(i4b)                  :: klo, khi
    if(s%linear) then
       klo = max(min(locate(s%x,x),size(s%x)-1),1)
       khi = klo+1
       y = s%y(klo) + (s%y(khi)-s%y(klo))*(x-s%x(klo))/(s%x(khi)-s%x(klo))
    elseif(s%regular) then
       y = splint_uniform_grid(s%x, s%y, s%y2, x)
    else
       y = splint(s%x, s%y, s%y2, x)
    end if
  end function

  subroutine splint_simple_multi(s, x, y)
    implicit none
    type(spline_type), intent(in) :: s
    real(dp),          intent(in) :: x(:)
    real(dp),          intent(out):: y(:)
    integer(i4b)                  :: i
    do i = 1, size(x)
       y(i) = splint(s, x(i))
    end do
  end subroutine

  subroutine free_spline(s)
    implicit none
    type(spline_type) :: s
    if(allocated(s%x))  deallocate(s%x)
    if(allocated(s%y))  deallocate(s%y)
    if(allocated(s%y2)) deallocate(s%y2)
  end subroutine

  ! Routines from Numerical Recipes
  subroutine spline_plain(x, y, yp1, ypn, y2)
    implicit none

    real(dp),               intent(in)  :: yp1, ypn
    real(dp), dimension(:), intent(in)  :: x, y
    real(dp), dimension(:), intent(out) :: y2

    integer(i4b) :: i, n, m
    real(dp), dimension(:), allocatable :: a, b, c, r

    n = size(x)
    allocate(a(n),b(n),c(n),r(n))

    c(1:n-1) = x(2:n)-x(1:n-1)
    r(1:n-1) = 6.d0*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1) = r(2:n-1)-r(1:n-2)
    a(2:n-1) = c(1:n-2)
    b(2:n-1) = 2.d0*(c(2:n-1)+a(2:n-1))
    b(1)     = 1.d0
    b(n)     = 1.d0

    if (yp1 > 0.99d30) then
       r(1) = 0.d0
       c(1) = 0.d0
    else
       r(1) = (3.d0 / (x(2)-x(1))) * ((y(2)-y(1))/(x(2)-x(1)) - yp1)
       c(1) = 0.5d0
    end if

    if (ypn > 0.99d30) then
       r(n) = 0.d0
       a(n) = 0.d0
    else
       r(n) = (-3.d0 / (x(n)-x(n-1))) * ((y(n)-y(n-1))/(x(n)-x(n-1)) - ypn)
       a(n) = 0.5d0
    end if

    call tridag(a(2:n), b(1:n), c(1:n-1), r(1:n), y2(1:n))

    deallocate(a,b,c,r)
  end subroutine spline_plain

  function splint_plain(xa, ya, y2a, x)
    implicit none

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: xa, ya, y2a
    real(dp)                            :: splint_plain

    integer(i4b) :: khi, klo, n
    real(dp)     :: a, b, h

    n = size(xa)

    klo = max(min(locate(xa,x),n-1),1)
    khi = klo+1

    h   = xa(khi) - xa(klo)
    a   = (xa(khi) - x) / h
    b   = (x - xa(klo)) / h
    
    splint_plain = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi))*(h**2)/6.d0

  end function splint_plain


  function splint_uniform_grid(xa, ya, y2a, x)
    implicit none

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: xa, ya, y2a
    real(dp)                            :: splint_uniform_grid

    integer(i4b) :: khi, klo, n
    real(dp)     :: a, b, h

    n   = size(xa)
    h   = xa(2) - xa(1)
    klo = min(max(int((x-xa(1))/h)+1,1),n-1)
    khi = klo+1
    a   = (xa(khi) - x) / h
    b   = (x - xa(klo)) / h
    
    splint_uniform_grid = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi))*(h**2)/6.d0

  end function splint_uniform_grid


  function splint_deriv(xa, ya, y2a, x)
    implicit none

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: xa, ya, y2a
    real(dp)                            :: splint_deriv

    integer(i4b) :: khi, klo, n
    real(dp)     :: a, b, h

    n = size(xa)

    klo = max(min(locate(xa,x),n-1),1)
    khi = klo+1
    h   = xa(khi) - xa(klo)
    a   = (xa(khi) - x) / h
    b   = (x - xa(klo)) / h

    splint_deriv = (ya(khi) - ya(klo)) / h - (3.d0 * a**2 - 1.d0) / 6.d0 * h * y2a(klo) + &
         & (3.d0 * b**2 - 1.d0) / 6.d0 * h * y2a(khi)

  end function splint_deriv


  subroutine splint_deriv_all_nodes(xa, ya, y2a, deriv)
    implicit none

    real(dp), dimension(:), intent(in)  :: xa, ya, y2a
    real(dp), dimension(:), intent(out) :: deriv

    integer(i4b) :: i, n

    n = size(xa)
    do i = 1, n-1
       deriv(i) = (ya(i+1)-ya(i)) / (xa(i+1)-xa(i)) - (xa(i+1)-xa(i)) * y2a(i) / 3.d0
    end do
    deriv(n) = (ya(n)-ya(n-1)) / (xa(n)-xa(n-1)) + (xa(n)-xa(n-1)) * y2a(n-1) / 3.d0

  end subroutine splint_deriv_all_nodes


  subroutine smooth_spline(weight, alpha, x, y, yp1, ypn, y2, variance)
    implicit none

    character(len=*),       intent(in)    :: weight
    real(dp),               intent(in)    :: alpha, yp1, ypn
    real(dp), dimension(:), intent(in)    :: x
    real(dp), dimension(:), intent(in), optional    :: variance
    real(dp), dimension(:), intent(inout) :: y
    real(dp), dimension(:), intent(out)   :: y2

    integer(i4b) :: i, j, n, row, col, kd, ldab, ldb, info, nrhs
    character(len=1) :: uplo
    real(dp), allocatable, dimension(:)   :: h, QtY, W, WQy2
    real(dp), allocatable, dimension(:,:) :: Q, R, M, QtQ
    real(dp),              dimension(3)   :: Q_row, Q_col

    n = size(x)

    allocate(QtY(2:n-1))
    allocate(h(1:n))
    allocate(W(n))
    allocate(Q(n,-1:1))
    allocate(WQy2(1:n))
    allocate(R(n,0:2))
    allocate(M(0:2,2:n-1))
    allocate(QtQ(n,0:2))
    
    ! Set up step length array h
    h = 0.d0
    do i = 1, n-1
       h(i) = x(i+1)-x(i)
    end do

    ! Compute QtY
    QtY = 0.d0
    do i = 2, n-1
       QtY(i) = (y(i+1)-y(i)) / h(i) - (y(i)-y(i-1)) / h(i-1)
    end do


    ! Set up weight array W
    if (trim(weight) == 'inv_var') then
 !      W(1) = (y(2)-y(1)) / (x(2)-x(1))
 !      do i = 2, n
 !         W(i) = (y(i)-y(i-1)) / (x(i)-x(i-1))
 !      end do

       do i = 1, n
          if (present(variance)) then
             W = variance
          else
             W(i) = y(i)**2
             !if (y(i) < 0.5d0) then
             !   W(i) = y(i)**2
             !else
             !   W(i) = (1.d0 - y(i))**2
             !end if
          end if
       end do
!       W(1) = 0.d0
!       W(n) = 0.d0

!       W = 0.d0
!       do i = 2, n-1
!          W(i) = y(i)
!       end do
    else if (trim(weight) == 'uniform') then
       W = 0.d0
       do i = 1, n
          W(i) = 1.d0
       end do
    else if (trim(weight) == 'N_weights') then
       W = 0.d0
       do i = 1, n
          if (variance(i) <= 1.d0) then
             W(i) = 1.d10
          else if (variance(i) < 6.d4) then
             W(i) = exp(-(variance(i)/1000.d0))  !1.d0/variance(i)**2
          end if
       end do
    else
       write(*,*) 'smooth_spline: Unknown weighting scheme = ', trim(weight)
       stop
    end if

!    write(*,*) W
!    stop

    ! Set up Q
    Q = 0.d0
    do j = 2, n-1
       Q(j,-1) =  1.d0 / h(j-1)               
       Q(j,0)  = -1.d0 / h(j-1) - 1.d0 / h(j) 
       Q(j,1)  =  1.d0 / h(j)                 
    end do

    ! Compute Q^T W Q
    QtQ = 0.d0
    do j = 1, n
       col = j

       Q_col = 0.d0
       if (col > 1) Q_col(1) = W(col-1) * Q(col,-1)
                    Q_col(2) = W(col)   * Q(col, 0)
       if (col < n) Q_col(3) = W(col+1) * Q(col, 1)


       do i = 0, 2
          row = j+i
          if (row > n) cycle

          Q_row = 0.d0
          if (row > 1 .and. -i-1 > -2) Q_row(1) = Q(row,-1-i)
          if (              -i   > -2) Q_row(2) = Q(row, 0-i)  
          if (row < n .and. -i+1 > -2) Q_row(3) = Q(row, 1-i)

          QtQ(j,i) = dot_product(Q_row,Q_col)
       end do
    end do


    ! Set up R
    R = 0.d0
    do i = 2, n-1
       R(i,0) = (h(i-1) + h(i)) / 3.d0
       if (i < n-1)  R(i,1) = h(i)/6.d0
    end do


    ! Compute R + alpha Q^t W Q
    M = 0.d0
    do i = 2, n-1
       do j = 0, 2
          M(j,i) = R(i,j) + alpha * QtQ(i,j)
       end do
    end do

    ! Cholesky decompose M
    uplo = 'l'
    kd   = 2
    ldab = kd+1
    ldb  = n-2
    nrhs = 1
    call dpbtrf(uplo, n-2, kd, M, ldab, info)

    if (info /= 0) then
       write(*,*) 'smooth_spline: info = ', info

       open(58,file='spline.dat')
       do i = 1, n
          write(58,*) x(i), y(i), W(i)
       end do
       close(58)
       stop
    end if

    ! Solve for second derivatives
    call dpbtrs(uplo, n-2, kd, nrhs, M, ldab, QtY, ldb, info)
    y2(1)     = 0.d0
    y2(2:n-1) = QtY
    y2(n)     = 0.d0

 
    ! Find the spline knot values
    WQy2 = 0.d0
    do i = 1, n
       if (i > 1) Q_row(1) = Q(i-1,1)
                  Q_row(2) = Q(i,  0)
       if (i < n) Q_row(3) = Q(i+1,-1)

       if (i > 1) Q_col(1) = y2(i-1)
                  Q_col(2) = y2(i)
       if (i < n) Q_col(3) = y2(i+1)

       WQy2(i) = W(i) * dot_product(Q_row,Q_col)
    end do

    y = y - alpha * WQy2

    deallocate(QtY)
    deallocate(h)
    deallocate(W)
    deallocate(Q)
    deallocate(R)
    deallocate(M)
    deallocate(QtQ)

  end subroutine smooth_spline

  subroutine qsimp(x_0, y_0, y2_0, a, b, s)
    implicit none

    real(dp), dimension(1:), intent(in)  :: x_0, y_0, y2_0
    real(dp),                intent(in)  :: a, b
    real(dp),                intent(out) :: s
    
    integer(i4b) :: j, JMAX = 20
    real(dp)     :: os, ost, st, eps = 1.d-6

    ost = -1.d30
    os  = -1.d30

    do j = 1, JMAX
       call trapzd(x_0, y_0, y2_0, a, b, st, j)
       s = (4.d0*st - ost)/3.d0
       if (j >= 5) then
          if ((abs(s-os) < eps*abs(os)) .or. ((s == 0.d0) .and. (os == 0.d0))) return
       end if
       os  = s
       ost = st
    end do

  end subroutine qsimp

  subroutine trapzd(x_0, y_0, y2_0, a, b, s, n)
    implicit none

    real(dp),     dimension(1:), intent(in)    :: x_0, y_0, y2_0
    real(dp),                    intent(in)    :: a, b
    real(dp),                    intent(inout) :: s
    integer(i4b),                intent(in)    :: n

    integer(i4b) :: it, j
    real(dp)     :: del, tot, tnm, x

    if (n == 1) then
       s = 0.5d0 * (b-a)*(splint(x_0, y_0, y2_0, a) + splint(x_0, y_0, y2_0, b))
    else
       it  = 2**(n-2)
       tnm = real(it,dp)
       del = (b-a)/tnm
       x   = a + 0.5d0*del
       tot = 0.d0
       do j = 1, it
          tot = tot + splint(x_0, y_0, y2_0, x)
          x   = x+del
       end do
       s = 0.5d0*(s + (b-a)*tot/tnm)
    end if
    return

  end subroutine trapzd


  function zriddr(x_0, y_0, y2_0, x1, x2, zeropt, xacc)
    implicit none

    real(dp), dimension(1:), intent(in) :: x_0, y_0, y2_0
    real(dp),                intent(in) :: x1, x2, xacc, zeropt
    real(dp)                            :: zriddr

    integer(i4b) :: maxit=60
    real(dp)     :: UNUSED=-1.11d30

    integer(i4b) :: j
    real(dp)     :: fh, fl, fm, fnew, s, xh, xl, xm, xnew

    fl = splint(x_0, y_0, y2_0, x1)-zeropt
    fh = splint(x_0, y_0, y2_0, x2)-zeropt

    if (((fl > 0.d0) .and. (fh < 0.d0)) .or. ((fl < 0.d0) .and. (fh > 0.d0))) then
       xl = x1
       xh = x2
       zriddr = UNUSED
       do j = 1, MAXIT
          xm = 0.5d0*(xl+xh)
          fm = splint(x_0, y_0, y2_0, xm)-zeropt
          s  = sqrt(fm**2 - fl*fh)
          if (s == 0.d0) return
          xnew = xm+(xm-xl)*(sign(1.d0,fl-fh)*fm/s)
          if (abs(xnew-zriddr) < xacc) return
          zriddr = xnew
          fnew   = splint(x_0, y_0, y2_0, zriddr)-zeropt
          if (fnew == 0.d0) return
          if (sign(fm,fnew) /= fm) then
             xl = xm
             fl = fm
             xh = zriddr
             fh = fnew
          else if (sign(fl,fnew) /= fl) then
             xh = zriddr
             fh = fnew
          else if (sign(fh,fnew) /= fh) then
             xl = zriddr
             fl = fnew
          else
             write(*,*) 'Should never be here in zriddr'
          end if
          if (abs(xh-xl) < xacc) return
       end do
    else

       open(68,file='func.dat')
       do j = 1, size(y_0)
          write(68,*) x_0(j), y_0(j)
       end do
       close(68)

       write(*,*) 'Root not bracketed in zriddr'
       write(*,*) 'xlow   = ', x1
       write(*,*) 'xhigh  = ', x2
       write(*,*) 'zeropt = ', zeropt
       stop

    end if

    return

  end function zriddr






end module spline_1D_mod
