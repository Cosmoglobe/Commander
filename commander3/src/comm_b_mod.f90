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
module comm_B_mod
  use comm_map_mod
  implicit none

  public comm_B, B_ptr

  type, abstract :: comm_B
     ! Data variables
     character(len=512)           :: type
     real(dp)                     :: r_max
     real(dp)                     :: mb_eff
     logical(lgt)                 :: almFromConv 
     class(comm_mapinfo), pointer :: info => null()
     real(dp),          allocatable, dimension(:,:) :: b_l
   contains
     ! Data procedures
     procedure(matmulB),     deferred :: conv
     procedure(matmulInvB),  deferred :: deconv
     procedure                     :: getBTheta
     procedure                     :: initBTheta
  end type comm_B

  abstract interface
     subroutine matmulB(self, trans, map)
       import comm_map, comm_B, dp, i4b, lgt
       implicit none
       class(comm_B),   intent(in)    :: self
       logical(lgt),    intent(in)    :: trans
       class(comm_map), intent(inout) :: map
     end subroutine matmulB

     subroutine matmulInvB(self, trans, map)
       import comm_map, comm_B, dp, i4b, lgt
       implicit none
       class(comm_B),   intent(in)    :: self
       logical(lgt),    intent(in)    :: trans
       class(comm_map), intent(inout) :: map
     end subroutine matmulInvB
  end interface

  type B_ptr
     class(comm_B), pointer :: p => null()
  end type B_ptr

  ! Local variables
  integer(i4b), parameter :: n_beam = 1000  ! Number of sample point in beam spline

contains

  ! Note: Assume b_l(P) = b_l(T) for now.
  subroutine initBTheta(self, b_l, filename)
    implicit none
    class(comm_B),                      intent(inout) :: self
    real(dp),         dimension(0:,1:), intent(in), optional  :: b_l
    character(len=*),                   intent(in), optional  :: filename

    
  end subroutine initBTheta

  function getBTheta(self, r, pol)
    implicit none
    class(comm_B), intent(in) :: self
    real(dp),      intent(in) :: r
    integer(i4b),  intent(in) :: pol
    real(dp)                  :: getBTheta

    getBTheta = 0d0

  end function getBTheta
    
  
end module comm_B_mod
