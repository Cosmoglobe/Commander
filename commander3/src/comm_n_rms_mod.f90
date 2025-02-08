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
module comm_N_rms_mod
  use comm_N_mod
  implicit none

  private
  public comm_N_rms, comm_N_rms_ptr
  
  type, extends (comm_N) :: comm_N_rms
     class(comm_map), pointer :: siN        => null()
     class(comm_map), pointer :: rms0       => null()
   contains
     procedure :: update_N    => update_N_rms
  end type comm_N_rms

  interface comm_N_rms
     procedure constructor
  end interface comm_N_rms

  type comm_N_rms_ptr
     type(comm_N_rms), pointer :: p => null()
  end type comm_N_rms_ptr
  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(info)
    implicit none
    class(comm_N_rms),                  pointer       :: constructor
    type(comm_mapinfo), target,         intent(in)    :: info

    ! General parameters
    allocate(constructor)

    call constructor%update_N(info)

  end function constructor

  subroutine update_N_rms(self, info)
    implicit none
    class(comm_N_rms),                   intent(inout)          :: self
    class(comm_mapinfo),                 intent(in)             :: info

    self%rms0     => comm_map(info)
    self%siN     => comm_map(self%rms0)

    self%invN_diag => comm_map(info)
    self%invN_diag%map = 0d0
    call compute_invN_lm(self%invN_diag)
    write(*,*) "Miracle"

  end subroutine update_N_rms





end module comm_N_rms_mod
