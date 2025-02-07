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
! Author: Sigurd Kirkevold Næss
!
!================================================================================
! This module implements a handy type of shared output file called a status
! file, which prints time, id, memory and lun information at given checkpoints in
! the file. Very useful for debugging resource leaks and timing the program.
! If, for some reason, you want a status file that does not actually do anything,
! pass an empty string for the file name.

module comm_status_mod
  use comm_shared_output_mod
  use comm_timing_mod
  implicit none

contains
end module
