module comm_tod_QUIET_mod
  !
  !   Module which contains all the QUIET time ordered data processing and routines
  !   for a given frequency band
  !
  !   Main Methods
  !   ------------
  !   constructor(cpar, id_abs, info, tod_type)
  !       Initialization routine that reads in, allocates and associates 
  !       all data needed for TOD processing
  !   process_QUIET_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
  !       Routine which processes the time ordered data
  !
  use comm_tod_mod
  use comm_param_mod
  use comm_map_mod
  use comm_conviqt_mod
  use pix_tools
  use healpix_types
  use comm_huffman_mod
  use comm_hdf_mod
  use comm_fft_mod
  use spline_1D_mod
  use comm_4D_map_mod
  use comm_zodi_mod
  use comm_tod_mapmaking_mod
  use comm_tod_pointing_mod
  use comm_tod_gain_mod
  use comm_tod_bandpass_mod
  use comm_tod_orbdipole_mod
  use comm_tod_driver_mod
  use comm_utils
  implicit none

  private
  public comm_QUIET_tod

  type, extends(comm_tod) :: comm_QUIET_tod
    character(len=20), allocatable, dimension(:) :: labels ! names of fields
  contains
    procedure     :: process_tod           => process_QUIET_tod
  end type comm_QUIET_tod

contains


  !**************************************************
  !             Constructor
  !**************************************************
  function constructor(cpar, id_abs, info, tod_type)
    !
    ! Constructor function that gathers all the instrument parameters in a pointer
    ! and constructs the objects
    !
    ! Arguments:
    ! ----------
    ! cpar:     derived type
    !           Object containing parameters from the parameterfile.
    ! id_abs:   integer
    !           The index of the current band within the parameters, related to cpar
    ! info:     map_info structure
    !           Information about the maps for this band, e.g.,
    !           how the maps are distributed in memory
    ! tod_type: string
    !           Instrument specific tod type
    !
    ! Returns
    ! ----------
    ! constructor: pointer
    !              Pointer that contains all instrument data
    !
    implicit none
    type(comm_params),      intent(in) :: cpar
    integer(i4b),           intent(in) :: id_abs
    class(comm_mapinfo),    target     :: info
    character(len=128),     intent(in) :: tod_type
    class(comm_QUIET_tod),  pointer    :: constructor


  end function constructor

  !**************************************************
  !             Driver routine
  !**************************************************
  subroutine process_QUIET_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
    !
    ! Routine that processes the QUIET time ordered data.
    ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
    ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms.
    ! Writes maps to disc in fits format
    !
    ! Arguments:
    ! ----------
    ! self:     pointer of comm_QUIET_tod class
    !           Points to output of the constructor
    ! chaindir: string
    !           Directory for output files
    ! chain:    integer
    !           Index number of the chain being run
    ! iter:     integer
    !           Gibbs iteration number
    ! handle:   planck_rng derived type
    !           Healpix definition for random number generation
    !           so that the same sequence can be resumed later on from that same point
    ! map_in:   array
    !           Array of dimension (ndet,ndelta) with pointer to maps,
    !           with both access to maps and changing them.
    !           ndet is the number of detectors and
    !           ndelta is the number of bandpass deltas being considered
    ! delta:    array
    !           Array of bandpass corrections with dimensions (0:ndet,npar,ndelta)
    !           where ndet is number of detectors, npar is number of parameters
    !           and ndelta is the number of bandpass deltas being considered
    !
    ! Returns:
    ! ----------
    ! map_out: comm_map class
    !          Final output map after TOD processing combined for all detectors
    ! rms_out: comm_map class
    !          Final output rms map after TOD processing combined for all detectors
    !
    implicit none
    class(comm_QUIET_tod),            intent(inout) :: self
    character(len=*),                    intent(in) :: chaindir
    integer(i4b),                        intent(in) :: chain, iter
    type(planck_rng),                 intent(inout) :: handle
    type(map_ptr), dimension(1:, 1:), intent(inout) :: map_in       ! (ndet,ndelta)
    real(dp),  dimension(0:, 1:, 1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
    class(comm_map),                  intent(inout) :: map_out      ! Combined output map
    class(comm_map),                  intent(inout) :: rms_out      ! Combined output rms


  end subroutine process_QUIET_tod

end module comm_tod_QUIET_mod
