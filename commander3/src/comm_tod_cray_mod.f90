!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University
! of Oslo.
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
module comm_tod_cray_mod
  use comm_param_mod
  implicit none

  private
  public comm_cray, cray_ptr

  type :: comm_cray
    real(dp), allocatable, dimension(:,:)     :: templates
    integer(i4b) :: ntemps
  contains
    procedure :: build_cray_templates
    procedure :: fit_cray_amplitudes
  end type comm_cray

  interface comm_cray 
    procedure constructor
  end interface comm_cray

  type cray_ptr
    class(comm_cray), pointer :: p => null()
  end type cray_ptr

contains

    ! Constructor

    function constructor(ntemps)
      implicit none
      integer(i4b),            intent(in) :: ntemps


      class(comm_cray), pointer           :: constructor

      allocate(constructor)

      constructor%ntemps = ntemps

    end function constructor

    ! builds a set of templates for cosmic ray response from the TOD
    subroutine build_cray_templates(self)
      implicit none
      class(comm_cray),                          intent(inout) :: self


    end subroutine build_cray_templates

    ! fits the constructed templates to the cosmic rays in the timestreams
    subroutine fit_cray_amplitudes(self)
      implicit none
      class(comm_cray),                          intent(inout) :: self

    end subroutine fit_cray_amplitudes

end module comm_tod_cray_mod
