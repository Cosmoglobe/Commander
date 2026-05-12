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

module comm_tod_Tbol_mod
  use comm_utils
  use comm_status_mod
  implicit none
  private

  public comm_Tbol, Tbol_ptr

  type :: comm_Tbol
  contains
    procedure :: estimate_params
    procedure :: convolve
  end type comm_Tbol   

  interface comm_Tbol
    procedure Tbol_constructor
  end interface comm_Tbol

  type Tbol_ptr
    class(comm_Tbol), pointer :: p => null()
  end type Tbol_ptr

contains
 
  !initializes a comm_Tbol class
  function Tbol_constructor() result(c)
    implicit none
    class(comm_Tbol), pointer                  :: c

    allocate(c)
    
  end function Tbol_constructor
 
  ! estimates Tbol parameters for each detector
  subroutine estimate_params(self)
    implicit none
    class(comm_Tbol), intent(inout)      :: self
  end subroutine estimate_params

  ! Routine for convolving with bolometer transfer function, tod = T * tod
  subroutine convolve(self, tod)
    implicit none
    class(comm_Tbol),               intent(in)    :: self
    real(sp),         dimension(:), intent(inout) :: tod
  end subroutine convolve

end module comm_tod_Tbol_mod
