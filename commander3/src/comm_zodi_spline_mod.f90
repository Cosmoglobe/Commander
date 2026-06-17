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
module comm_zodi_spline_mod
   use spline_1d_mod
   use comm_defs
   implicit none

   private
   public zodi_spline

   type :: zodi_spline
      integer(i4b) :: knots  !number of spline knots inside the interval, should depend on how wide is the spline interval
      real(dp) :: nu_left, nu_right !frequency where MBB becomes spline and spline becomes BB
      real(dp) :: sed_ampl, sed_b, T_zodi !parameters of zodi MBB and BB
      type(spline_type) :: spl

   contains
      procedure :: build
      procedure, private :: update
      procedure :: eval
   end type zodi_spline

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
! TODO: Add error checks

   subroutine build(self, knots, nu_left, nu_right, sed_ampl, sed_b, T_zodi)
      implicit none
      class(zodi_spline), intent(inout) :: self
      integer(i4b), intent(in) :: knots
      real(dp), intent(in) :: nu_left, nu_right, sed_ampl, sed_b, T_zodi

      self%knots = knots
      self%nu_left = nu_left
      self%nu_right = nu_right
      self%sed_ampl = sed_ampl
      self%sed_b = sed_b
      self%T_zodi = T_zodi

      call self%update()
    
   end subroutine build

   subroutine update(self)
      implicit none
      class(zodi_spline), intent(inout) :: self
      integer(i4b) :: n_k, i !number of knots+limits of the integer, variable used in do loop
      real(dp), allocatable :: x(:), y(:) !they'll be zs%spl%x and y
      real(dp) :: xbb, delta_nu, half_nu !exponential in BB expression, width of spline interval portion, central frequency
      real(dp) :: mbb_prime, bb_prime !derivatives at boundaries

      if (self%knots < 2) then 
         n_k = 3 !at least knots=[left, right, midpoint]
      else 
         n_k = self%knots+2 
      end if 

      allocate(x(n_k), y(n_k))

      !The function is MBB -> sed_left -> spline -> sed_right -> BB

      x(1) = self%nu_left
      xbb = (h*x(1))/(k_b*self%T_zodi)
      if (xbb > EXP_OVERFLOW) then
         y(1) = 0.d0
         mbb_prime = 0.d0
      else 
         y(1) = ((2.d0*h*x(1)**3)/c**2)*(1.d0/(exp(xbb)-1.d0))*self%sed_ampl*x(1)**self%sed_b
         mbb_prime = ((2.d0*h*self%sed_ampl*x(1)**(self%sed_b+2))/c**2)*(1.d0/(exp(xbb)-1.d0)**2)*((self%sed_b+3)*(exp(xbb)-1.d0)-(xbb*exp(xbb)))
      end if 
      
      x(n_k) = self%nu_right
      xbb = (h*x(n_k))/(k_b*self%T_zodi)
      if (xbb > EXP_OVERFLOW) then
         y(n_k) = 0.d0
         bb_prime = 0.d0
      else 
         y(n_k) = ((2.d0*h*x(n_k)**3)/c**2)*(1.d0/(exp(xbb)-1.d0))
         bb_prime = ((2.d0*h*x(n_k)**2)/c**2)*(1.d0/(exp(xbb)-1.d0)**2)*(3.d0*(exp(xbb)-1.d0)-(xbb*exp(xbb)))
      end if
      
      delta_nu = (self%nu_right - self%nu_left)/real(n_k-1, dp)
      half_nu = self%nu_left+((self%nu_right - self%nu_left)/2.d0)
      do i = 2, n_k-1
         x(i) = x(i-1)+delta_nu
         xbb = (h*x(i))/(k_b*self%T_zodi)
         if (xbb > EXP_OVERFLOW) then
            y(i) = 0.d0
         else
            !Here I assume that the joint between MBB and BB is the midpoint of the spline interval, maybe this is wrong
            if (x(i) < half_nu) then 
               y(i) = ((2.d0*h*x(i)**3)/c**2)*(1.d0/(exp(xbb)-1.d0))*self%sed_ampl*x(i)**self%sed_b
            else 
               y(i) = ((2.d0*h*x(i)**3)/c**2)*(1.d0/(exp(xbb)-1.d0))
            end if
         end if
      end do

      call spline(self%spl, x, y, boundary=[mbb_prime, bb_prime], regular=.false., linear=.false.)
      deallocate(x,y)

   end subroutine  update

   function eval(self, new_nu)
      implicit none
      class(zodi_spline), intent(in) :: self
      real(dp), intent(in) :: new_nu
      real(dp) :: eval

      !check that new_nu is actually in spline interval
      if (.not.((self%nu_left < new_nu) .and. (new_nu < self%nu_right))) then 
         write(*,*) "WARNING! zodi_spline called for a frequency not in zodi_spline_interval."
         eval = 0.d0 !CHECK THIS
         return
      end if 

      eval = splint(self%spl, new_nu)

   end function eval

end module comm_zodi_spline_mod