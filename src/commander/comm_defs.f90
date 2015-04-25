module comm_defs
  use healpix_types
  implicit none

  !**************************************************
  !               Global variables
  !**************************************************
  real(dp), parameter :: k_B      = 1.3806503d-23
  real(dp), parameter :: h        = 1.0545726691251021d-34 * 2.d0*pi 
  real(dp), parameter :: c        = 2.99792458d8
  real(dp)            :: T_CMB    = 2.7255d0

  !**************************************************
  !               Counters
  !**************************************************
  integer(i4b), parameter :: GAIN     = 1
  integer(i4b), parameter :: NOISEAMP = 2
  integer(i4b), parameter :: MONOPOLE = 3
  integer(i4b), parameter :: DIPOLE_X = 4
  integer(i4b), parameter :: DIPOLE_Y = 5
  integer(i4b), parameter :: DIPOLE_Z = 6
  

  
end module comm_defs
