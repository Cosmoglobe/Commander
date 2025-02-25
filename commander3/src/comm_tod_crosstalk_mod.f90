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

module comm_tod_crosstalk_mod
  use comm_tod_mod
  use comm_status_mod
  use comm_tod_driver_mod
  implicit none
  private

  public comm_crosstalk

  type :: comm_crosstalk
    real(dp), dimension (:,:), allocatable      :: crosstalk_matrix
  contains
    procedure :: estimate_crosstalk_matrix
    procedure :: remove_crosstalk_signal
  end type comm_crosstalk   

  interface comm_crosstalk
    procedure xtalk_constructor
  end interface comm_crosstalk

  type crosstalk_pointer
    class(comm_crosstalk), pointer :: p => null()
  end type crosstalk_pointer

contains
 
  !initializes a comm_crosstalk class
  function xtalk_constructor(correlations) result(c)
    implicit none
    logical(lgt), dimension(:,:),    intent(in)     :: correlations
    class(comm_crosstalk), pointer                  :: c

    allocate(c)
    allocate(c%crosstalk_matrix(size(correlations(1,:)), size(correlations(:,1))))

  end function xtalk_constructor
 
  ! estimates the crosstalk coeficients between each detector
  subroutine estimate_crosstalk_matrix(self)
    implicit none
    class(comm_crosstalk),                            intent(in)      :: self


  end subroutine estimate_crosstalk_matrix

  subroutine remove_crosstalk_signal(self, sd, i)
    implicit none
    class(comm_crosstalk),                      intent(in)      :: self
    class(comm_scandata),                       intent(inout)   :: sd
    integer(i4b),                               intent(in)      :: i

  end subroutine remove_crosstalk_signal

end module comm_tod_crosstalk_mod
