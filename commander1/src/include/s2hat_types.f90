
module s2hat_types

  ! definitions needed for s2hat package
  ! Radek Stompor (APC) February, 2007

  use healpix_types
  
  implicit none

  include "mpif.h"

  ! define pixel choices here - HEALPIX the only one so far ...
  integer(i4b), parameter :: PIXCHOICE_HEALPIX = 0  ! for HEALPIX pixelization
  integer(i4b), parameter :: PIXCHOICE_GLESP = 1    ! for GLESP pixelization
  integer(i4b), parameter :: PIXCHOICE_ECP = 2      ! for ECP pixelization
  integer(i4b), parameter :: PIXCHOICE_GLCP = 3     ! for GLCP pixelization

  ! spin convention 
  integer(i4b), parameter :: SPIN_CONV_SIGN = -1    ! for Healpix like choice

  type pixparameters
    integer(i4b) :: par1
    integer(i4b) :: par2
  end type pixparameters

  type pixeltype                        ! pixelization is assumed to be iso-latitudal with pixels evenly spaced for each latitude
					! and a symmetry present between northern and southern hemispheres.
    ! pixel types defined as above - 8 byte ints to avoid padding problems in the C-interface
    integer(i8b) :: type
    ! a total number of pixels
    integer(i8b) :: npixsall
    ! a total number of iso-latitude rings in the north hemisphere (including equator)
    integer(i8b) :: nringsall
    ! a total maximum number of pixels per iso-ring
    integer(i8b) :: nphmx
    ! a number of the first pixel for each ring
    integer(i8b), dimension(:), pointer :: fpix  ! ranges [1:nringsall]
    ! a number of pixel for each iso-ring
    integer(i8b), dimension(:), pointer :: nph   ! ranges [1:nringsall]
    ! an offset of the 1st pixel for each ring with respect to a meridian 0 (radians)
    real(dp), dimension(:), pointer :: kphi  ! ranges [1:nringsall]
    ! quadrature weights for each ring
    real(dp), dimension(:), pointer :: qwght  ! ranges [1:nringsall]    
    ! pixel center separations for each iso-ring
    real(dp), dimension(:), pointer :: pixphi    ! ranges [1:nringsall]
    ! pixel area (assumed to be constant) for each iso-ring
    real(dp), dimension(:), pointer :: parea     ! ranges [1:nringsall]
    ! cosines of the polar angle for each iso-ring
    real(dp), dimension(:), pointer :: cth       ! ranges [1:nringsall]
    ! sines of the polar angle for each iso-ring (redundant)
    real(dp), dimension(:), pointer :: sth       ! ranges [1:nringsall]
  end type pixeltype
  
  type scandef     ! defines the scan parameters needed for the s2hat transforms
     ! a total number of observed pixels ! all 8byte int to simplify the C-interfacing
     integer(i8b) :: npixsobs
     ! a total number of observed rings
     integer(i8b) :: nringsobs
     ! observed flags:
     ! northern hemisphere (includes equator)
     integer(i8b), dimension(:), pointer :: nfl ! ranges [1:nringsall]
     ! southern hemisphere 
     integer(i8b), dimension(:), pointer :: sfl ! ranges [1:nringsall]
     ! either northern or southern (redundant)
     integer(i8b), dimension(:), pointer :: fl  ! ranges [1:nringsall]
  end type scandef

  type( pixeltype), target, public :: globalPixelization
  type( scandef), target, public :: globalScan

end module s2hat_types
