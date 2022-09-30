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
     logical(lgt)                        :: verbose
  end type spline_type
  
  interface spline
     module procedure spline_simple, spline_plain
  end interface spline
  
  interface splint
     module procedure splint_simple, splint_plain
  end interface splint
  
  interface splint_multi
     module procedure splint_simple_multi
  end interface splint_multi
  
contains
  
  ! Simple spline interface routines. Uses the spline
  ! type which makes internal copies, so uses slightly
  ! more memory than necessary, but probably not a big
  ! deal. And saves you from having to pass in 3 arguments
  ! all the time.
  subroutine spline_simple(s, x, y, boundary, regular, linear, verbose)
    implicit none
    type(spline_type),       intent(inout) :: s
    real(dp),                intent(in)    :: x(:), y(:)
    real(dp),      optional, intent(in)    :: boundary(2)
    logical(lgt),  optional, intent(in)    :: regular, linear, verbose
    call free_spline(s)
    s%boundary = 1d30
    s%regular  = .false.
    s%linear   = .false.
    s%verbose  = .false.
    if(present(boundary)) s%boundary = boundary
    if(present(regular))  s%regular  = regular
    if(present(linear))   s%linear   = linear
    if(present(verbose))  s%verbose  = verbose
    allocate(s%x(size(x)),s%y(size(x)),s%y2(size(x)))
    s%x = x
    s%y = y
    if(.not. s%linear) call spline(s%x, s%y, s%boundary(1), s%boundary(2), s%y2)
  end subroutine spline_simple
  
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
       y = splint(s%x, s%y, s%y2, x, s%verbose)
    end if
    if(s%verbose) write(*,*) "splint_simple:", x, y
  end function splint_simple

  subroutine splint_simple_multi(s, x, y)
    implicit none
    type(spline_type), intent(in) :: s
    real(dp),          intent(in) :: x(:)
    real(dp),          intent(out):: y(:)
    integer(i4b)                  :: i
    do i = 1, size(x)
       y(i) = splint(s, x(i))
    end do
  end subroutine splint_simple_multi
  
  subroutine free_spline(s)
    implicit none
    type(spline_type) :: s
    if(allocated(s%x))  deallocate(s%x)
    if(allocated(s%y))  deallocate(s%y)
    if(allocated(s%y2)) deallocate(s%y2)
  end subroutine free_spline
  
  ! Routines from Numerical Recipes
  subroutine spline_plain(x, y, yp1, ypn, y2)
    implicit none
    
    real(dp),               intent(in)  :: yp1, ypn
    real(dp), dimension(:), intent(in)  :: x, y
    real(dp), dimension(:), intent(out) :: y2
    
    integer(i4b) :: n
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

  function splint_plain(xa, ya, y2a, x, verbose)
    implicit none

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: xa, ya, y2a
    real(dp)                            :: splint_plain
    logical(lgt), optional, intent(in)  :: verbose

    integer(i4b) :: khi, klo, n
    real(dp)     :: a, b, h
    logical(lgt) :: verb

    if(.not. present(verbose)) then
      verb = .false.
    else
      verb = verbose
    end if

    n = size(xa)

    klo = max(min(locate(xa,x),n-1),1)
    khi = klo+1

    h   = xa(khi) - xa(klo)
    a   = (xa(khi) - x) / h
    b   = (x - xa(klo)) / h
    
    splint_plain = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi))*(h**2)/6.d0

    if(verb) then 
      write(*,*) "splint_plain", x, splint_plain, xa(klo), xa(khi), ya(klo), ya(khi), locate(xa, x)
    end if

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
