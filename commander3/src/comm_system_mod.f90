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

module comm_system_mod
  use healpix_types
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer(i8b), external :: get_mem_use, get_max_mem_use, get_mem_use2
  integer(i4b), external :: get_pid, open_atomic_file, nfork, popen, svn_revision, ishift, count_set_bits
  logical(lgt), external :: is_dir, mkdir_c, mv_c, rm_c
  integer(i4b), parameter :: SIGINT = 2, SIGSEGV = 11, SIGTERM = 15, SIGBUS = 7
  real(c_double),     bind(C, name="nan")          :: nan
  real(c_double),     bind(C, name="snan")         :: snan
  real(c_double)                                   :: infinity 
  integer(c_int64_t), bind(C, name="fe_divbyzero") :: fe_divbyzero
  integer(c_int64_t), bind(C, name="fe_inexact")   :: fe_inexact
  integer(c_int64_t), bind(C, name="fe_nan")       :: fe_nan
  integer(c_int64_t), bind(C, name="fe_overflow")  :: fe_overflow
  integer(c_int64_t), bind(C, name="fe_underflow") :: fe_underflow

contains

  subroutine fsleep(t)
    implicit none
    real(dp)     :: t
    integer(i4b) :: i
    i = nint(t*1d6)
    call usleep(i)
  end subroutine


end module
