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
module comm_bp_utils
  use comm_utils
  use comm_hdf_mod
  implicit none 

  interface comp_sz_thermo
     module procedure compute_sz_thermo_single, compute_sz_thermo_array
  end interface comp_sz_thermo

  interface comp_a2t
     module procedure compute_ant2thermo_single
  end interface comp_a2t

  interface comp_bnu_prime
     module procedure compute_bnu_prime_single, compute_bnu_prime_array
  end interface comp_bnu_prime

  interface comp_bnu_prime_RJ
     module procedure compute_bnu_prime_RJ_single, compute_bnu_prime_RJ_array
  end interface comp_bnu_prime_RJ
  
contains

  function compute_sz_thermo_single(nu)
    implicit none
    real(dp), intent(in) :: nu
    real(dp)             :: compute_sz_thermo_single

    real(dp) :: x

    compute_sz_thermo_single = 0d0
    
  end function compute_sz_thermo_single

  function compute_sz_thermo_array(nu) result (y)
    implicit none
    real(dp),              dimension(:), intent(in)  :: nu
    real(dp), allocatable, dimension(:)              :: y

    real(dp), allocatable, dimension(:) :: x

    y = 0d0

  end function compute_sz_thermo_array
  

  function compute_ant2thermo_single(nu)
    implicit none
    real(dp), intent(in)  :: nu
    real(dp)              :: compute_ant2thermo_single

    real(dp)     :: x

    compute_ant2thermo_single = 0d0
    
  end function compute_ant2thermo_single

  function compute_bnu_prime_array(nu)
    implicit none

    real(dp), dimension(:),        intent(in)  :: nu
    real(dp), dimension(size(nu))              :: compute_bnu_prime_array

    integer(i4b) :: i
    real(dp)     :: x

    compute_bnu_prime_array = 0d0

    
  end function compute_bnu_prime_array

  function compute_bnu_prime_single(nu)
    implicit none

    real(dp), intent(in)  :: nu
    real(dp)              :: compute_bnu_prime_single

    real(dp)     :: x

    compute_bnu_prime_single = 0d0
    
  end function compute_bnu_prime_single

  function compute_bnu_prime_RJ_array(nu)
    implicit none

    real(dp), dimension(:),        intent(in)  :: nu
    real(dp), dimension(size(nu))              :: compute_bnu_prime_RJ_array

    compute_bnu_prime_RJ_array = 2.d0*k_B*nu**2/c**2
    
  end function compute_bnu_prime_RJ_array

  function compute_bnu_prime_RJ_single(nu)
    implicit none

    real(dp), intent(in)  :: nu
    real(dp)              :: compute_bnu_prime_RJ_single

    compute_bnu_prime_RJ_single = 0d0
    
  end function compute_bnu_prime_RJ_single

  function dBdnu(nu)
    implicit none

    real(dp), intent(in)  :: nu
    real(dp)              :: dBdnu

    real(dp)     :: x

    dBdnu = 0d0
    
  end function dBdnu

  function dB_rj_dnu(nu)
    implicit none

    real(dp), intent(in)  :: nu
    real(dp)              :: dB_rj_dnu

    dB_rj_dnu = 2*nu**2*k_b / c**2
    
  end function dB_rj_dnu

  ! Routine for reading bandpass files for one detecor with threshold


end module comm_bp_utils
