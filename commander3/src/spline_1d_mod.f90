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
  subroutine spline_simple(s, x, y, boundary, regular, linear, verbose)
    implicit none
    type(spline_type),       intent(inout) :: s
    real(dp),                intent(in)    :: x(:), y(:)
    real(dp),      optional, intent(in)    :: boundary(2)
    logical(lgt),  optional, intent(in)    :: regular, linear, verbose
  end subroutine

  function splint_simple(s, x) result(y)
    implicit none
    type(spline_type), intent(in) :: s
    real(dp),          intent(in) :: x
    real(dp)                      :: y
    integer(i4b)                  :: klo, khi

    y = 0d0

  end function

  subroutine splint_simple_multi(s, x, y)
    implicit none
    type(spline_type), intent(in) :: s
    real(dp),          intent(in) :: x(:)
    real(dp),          intent(out):: y(:)
    integer(i4b)                  :: i
    y = 0d0
  end subroutine

  subroutine free_spline(s)
    implicit none
    type(spline_type) :: s
  end subroutine

  ! Routines from Numerical Recipes
  subroutine spline_plain(x, y, yp1, ypn, y2)
    implicit none

    real(dp),               intent(in)  :: yp1, ypn
    real(dp), dimension(:), intent(in)  :: x, y
    real(dp), dimension(:), intent(out) :: y2

    integer(i4b) :: n
    real(dp), dimension(:), allocatable :: a, b, c, r

    y2 = 0d0

  end subroutine spline_plain

  function splint_plain(xa, ya, y2a, x, verbose)
    implicit none

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: xa, ya, y2a
    real(dp)                            :: splint_plain
    logical(lgt), optional, intent(in)  :: verbose

    splint_plain = 0d0


  end function splint_plain


  function splint_uniform_grid(xa, ya, y2a, x)
    implicit none

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: xa, ya, y2a
    real(dp)                            :: splint_uniform_grid

    integer(i4b) :: khi, klo, n
    real(dp)     :: a, b, h

    
    splint_uniform_grid = 0d0

  end function splint_uniform_grid


  function splint_deriv(xa, ya, y2a, x)
    implicit none

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: xa, ya, y2a
    real(dp)                            :: splint_deriv


    splint_deriv = 0d0

  end function splint_deriv


  subroutine splint_deriv_all_nodes(xa, ya, y2a, deriv)
    implicit none

    real(dp), dimension(:), intent(in)  :: xa, ya, y2a
    real(dp), dimension(:), intent(out) :: deriv

    integer(i4b) :: i, n

    deriv = 0d0
  end subroutine splint_deriv_all_nodes



  subroutine qsimp(x_0, y_0, y2_0, a, b, s)
    implicit none

    real(dp), dimension(1:), intent(in)  :: x_0, y_0, y2_0
    real(dp),                intent(in)  :: a, b
    real(dp),                intent(out) :: s
  
    s = 0d0
  end subroutine qsimp

  subroutine trapzd(x_0, y_0, y2_0, a, b, s, n)
    implicit none

    real(dp),     dimension(1:), intent(in)    :: x_0, y_0, y2_0
    real(dp),                    intent(in)    :: a, b
    real(dp),                    intent(inout) :: s
    integer(i4b),                intent(in)    :: n

    integer(i4b) :: it, j
    real(dp)     :: del, tot, tnm, x

    s = 0d0
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

    zriddr = 0d0

  end function zriddr






end module spline_1D_mod
