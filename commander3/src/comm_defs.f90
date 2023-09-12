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
module comm_defs
  use healpix_types
  implicit none

  !**************************************************
  !               Global variables
  !**************************************************
  real(dp), parameter :: k_B       = 1.3806503d-23
  real(dp), parameter :: h         = 1.0545726691251021d-34 * 2.d0*pi 
  real(dp), parameter :: c         = 2.99792458d8
  real(dp)            :: T_CMB     = 2.72548d0
  real(dp)            :: T_CMB_DIP = 3359.5d-6
  real(dp)            :: v_solar(3)= [0.d0, 0.d0, 0.d0]
  real(dp)            :: EXP_OVERFLOW = 700.d0
  integer(i4b)        :: zodi_nside = 64
  !**************************************************
  !               Counters
  !**************************************************
  integer(i4b), parameter :: GAIN     = 1
  integer(i4b), parameter :: NOISEAMP = 2
  
end module comm_defs
