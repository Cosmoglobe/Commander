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
module comm_F_int_0D_mod
  use comm_F_int_mod
  use comm_comp_mod
  use comm_bp_mod
  use comm_utils
  implicit none

  private
  public comm_F_int_0D

  type, extends (comm_F_int) :: comm_F_int_0D
     real(dp)                  :: f_precomp
     class(comm_comp), pointer :: comp => null()
   contains
     ! Data procedures
     procedure :: eval         => evalIntF
     procedure :: eval_deriv   => evalIntdF
     procedure :: update       => updateIntF
  end type comm_F_int_0D

  interface comm_F_int_0D
     procedure constructor
  end interface comm_F_int_0D

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(comp, bp, pol, f_precomp)
    implicit none
    class(comm_comp),     intent(in), target   :: comp
    class(comm_bp),       intent(in), target   :: bp
    integer(i4b),         intent(in)           :: pol
    real(dp),             intent(in), optional :: f_precomp
    class(comm_F_int_0D), pointer              :: constructor

    integer(i4b) :: j, m
    real(dp), allocatable, dimension(:) :: s
    
    allocate(constructor)
    constructor%comp => comp
    constructor%bp   => bp

    if (present(f_precomp)) then
       call constructor%update(f_precomp)
    else
       call constructor%update
    end if

  end function constructor

  
  ! Evaluate SED in brightness temperature normalized to reference frequency
  function evalIntF(self, theta)
    class(comm_F_int_0D),                intent(in) :: self
    real(dp),             dimension(1:), intent(in) :: theta
    real(dp)                                        :: evalIntF
    evalIntF = self%f_precomp
  end function evalIntF

  ! Evaluate partial derivative of SED in brightness temperature normalized to reference frequency
  function evalIntdF(self, theta, par)
    class(comm_F_int_0D),                intent(in) :: self
    real(dp),             dimension(1:), intent(in) :: theta
    integer(i4b),                        intent(in) :: par
    real(dp)                                        :: evalIntdF
    evalIntdF = 0.d0
  end function evalIntdF

  ! Compute/update integration look-up tables
  subroutine updateIntF(self, f_precomp, pol)
    class(comm_F_int_0D), intent(inout)           :: self
    real(dp),             intent(in),    optional :: f_precomp
    integer(i4b),         intent(in),    optional :: pol

    integer(i4b) :: m, j
    real(dp), allocatable, dimension(:) :: s

    if (present(f_precomp)) then
       self%f_precomp = f_precomp
    else
       ! Evaluate the bandpass integrated SED
       m = self%bp%n
       allocate(s(m))
       do j = 1, m
          s(j) = self%comp%S(nu=self%bp%nu(j))
       end do
       self%f_precomp = self%bp%SED2F(s)
       deallocate(s)
    end if
  end subroutine updateIntF

end module comm_F_int_0D_mod
