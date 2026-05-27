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
  real(dp)            :: EXP_OVERFLOW = 700.d0
  real(dp)            :: SECOND_TO_DAY = 1.1574074074074073e-05
  integer(i4b)        :: zodi_nside = 512

  !**************************************************
  !               Counters
  !**************************************************
  integer(i4b), parameter :: GAIN     = 1
  integer(i4b), parameter :: NOISEAMP = 2

  !**************************************************
  !               Constants
  !**************************************************

  integer(i4b), parameter :: NBIN_EARTH_ELON = 180


  !**************************************************
  !               TOD mask bit definitions
  !   bit 0 = good pixel, bit 1 = bad pixel
  !**************************************************  
  integer(i4b), parameter, public :: TODMASK_MAP        = 0 ! Defines map footprint
  integer(i4b), parameter, public :: TODMASK_NCORR      = 1 
  integer(i4b), parameter, public :: TODMASK_GAIN       = 2
  integer(i4b), parameter, public :: TODMASK_BANDPASS   = 3
  integer(i4b), parameter, public :: TODMASK_ADC        = 4
  integer(i4b), parameter, public :: TODMASK_ZODI       = 5
  integer(i4b), parameter, public :: TODMASK_PROC       = 6  

  !**************************************************
  !           Scan data bit definitions
  !**************************************************
  integer(i4b), parameter :: SD_TOT      =  0
  integer(i4b), parameter :: SD_HORNS    =  1 ! Compute individual horn streams 
  integer(i4b), parameter :: SD_BASE     =  2 ! Pointing and flags
  integer(i4b), parameter :: SD_IND      =  3
  integer(i4b), parameter :: SD_SPUR     =  4
  integer(i4b), parameter :: SD_MASK     =  5
  integer(i4b), parameter :: SD_TOD      =  6
  integer(i4b), parameter :: SD_SKY      =  7
  integer(i4b), parameter :: SD_NCORR    =  9
  integer(i4b), parameter :: SD_BP       = 10
  integer(i4b), parameter :: SD_BP_PROP  = 11
  integer(i4b), parameter :: SD_ZODI     = 12
  integer(i4b), parameter :: SD_DARK     = 13
  integer(i4b), parameter :: SD_SL       = 14
  integer(i4b), parameter :: SD_SKY_PROP = 15
  integer(i4b), parameter :: SD_ORB      = 16
  integer(i4b), parameter :: SD_INST     = 17
  integer(i4b), parameter :: SD_GAIN     = 18
  integer(i4b), parameter :: SD_CALIB    = 19
  integer(i4b), parameter :: SD_MONO     = 20
  integer(i4b), parameter :: SD_OBJCTR   = 21
  integer(i4b), parameter :: SD_JUMP     = 22
  integer(i4b), parameter :: SD_SPIKE    = 23

  !**************************************************
  !           Ephemeris definitions
  !**************************************************
  integer(i4b), parameter :: EPH_NUM_OBJECTS = 2
  integer(i4b), parameter :: EPH_JUPITER     = 1
  integer(i4b), parameter :: EPH_VENUS       = 2

  !**************************************************
  !           TOD data bit definitions
  !**************************************************
  integer(i4b), parameter :: TOD_RAMP_RESET      =  5
  integer(i4b), parameter :: TOD_CALLAMP1        =  15
  integer(i4b), parameter :: TOD_CALLAMP2        =  16


end module comm_defs
