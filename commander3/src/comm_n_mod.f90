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
module comm_N_mod
  use comm_map_mod
  implicit none

  type :: comm_N
     ! Data variables
     character(len=512)       :: type
     integer(i4b)             :: nside, nmaps, np, npix, comm
     class(comm_map),     pointer :: invN_diag => null()
     class(comm_map),     pointer :: rms_reg   => null()
     class(comm_mapinfo), pointer :: info      => null()
  end type comm_N

  type comm_N_ptr
     class(comm_N), pointer :: p => null()
  end type comm_N_ptr


contains

  subroutine compute_invN_lm(invN_diag)
    implicit none

    class(comm_map),  intent(inout) :: invN_diag

    call invN_diag%YtW_scalar
    
  end subroutine compute_invN_lm

  
end module comm_N_mod
