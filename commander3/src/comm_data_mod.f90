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
module comm_data_mod
  use comm_N_mod
  use comm_N_rms_mod
  implicit none

  type comm_data_set
     class(comm_mapinfo), pointer :: rmsinfo   => null()
     class(comm_map),     pointer :: map       => null()
     class(comm_map),     pointer :: mask      => null()
     class(comm_N),       pointer :: N         => null()
  end type comm_data_set

  type(comm_data_set), allocatable, dimension(:) :: data
  
contains

  subroutine initialize_data_mod(cpar)
    !
    ! Routine to initialise Commander3 data
    !
    implicit none
    type(comm_params), intent(in)    :: cpar

    allocate(data(1))

    data(1)%rmsinfo  => comm_mapinfo(cpar%comm_chain, 256, 750, 1)
    data(1)%mask     => comm_map(data(1)%rmsinfo)
    data(1)%N        => comm_N_rms(data(1)%rmsinfo)

  end subroutine initialize_data_mod

end module comm_data_mod
