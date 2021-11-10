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

  ! comm_tod is an abstarct type
  type, extends(comm_tod) :: comm_QUIET_tod
    character(len=20), allocatable, dimension(:) :: labels ! names of fields
  contains
    procedure :: process_tod    => process_QUIET_tod
    procedure :: read_tod_inst  => read_tod_inst_QUIET
    procedure :: read_scan_inst => read_scan_inst_QUIET
    procedure :: initHDF_inst   => initHDF_QUIET
    procedure :: dumpToHDF_inst => dumpToHDF_QUIET
  end type comm_QUIET_tod

  ! Making it as interface to implement in submodule
  ! Note: Why? Because I want to :) Jokes aside -- to 
  ! illustrate the use of submodules. It should potentially
  ! make recompillation faster and code easier to maintain
  ! as we split one big module into separate instances.
  interface comm_QUIET_tod
    module function constructor(cpar, id_abs, info, tod_type)
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
      implicit none
      type(comm_params),      intent(in) :: cpar
      integer(i4b),           intent(in) :: id_abs
      class(comm_mapinfo),    target     :: info
      character(len=128),     intent(in) :: tod_type
      class(comm_QUIET_tod),  pointer    :: constructor

    end function constructor
  end interface comm_QUIET_tod

  interface
    module subroutine process_QUIET_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
      !
      ! Routine that processes the QUIET time ordered data.
      ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
      ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps 
      ! and rms.
      ! 
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
      !
      implicit none
      class(comm_QUIET_tod),            intent(inout) :: self
      character(len=*),                 intent(in)    :: chaindir
      integer(i4b),                     intent(in)    :: chain, iter
      type(planck_rng),                 intent(inout) :: handle
      type(map_ptr), dimension(1:, 1:), intent(inout) :: map_in    ! (ndet,ndelta)
      real(dp),  dimension(0:, 1:, 1:), intent(inout) :: delta     ! (0:ndet,npar,ndelta) BP corrections
      class(comm_map), intent(inout) :: map_out      ! Combined output map
      class(comm_map), intent(inout) :: rms_out      ! Combined output rms
      type(map_ptr),   dimension(1:,1:),   intent(inout), optional :: map_gain       ! (ndet,1)

    end subroutine process_QUIET_tod
  end interface 

contains

  subroutine read_tod_inst_QUIET(self, file)
    ! 
    ! Reads QUIET-specific common fields from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_QUIET_tod)
    !           QUIET-specific TOD object
    ! file:     derived type (hdf_file)
    !           Already open HDF file handle; only root includes this
    !
    ! Returns
    ! ----------
    ! None, but updates self
    !
    implicit none
    class(comm_QUIET_tod),               intent(inout)          :: self
    type(hdf_file),                      intent(in),   optional :: file
  end subroutine read_tod_inst_QUIET
  
  subroutine read_scan_inst_QUIET(self, file, slabel, detlabels, scan)
    ! 
    ! Reads QUIET-specific scan information from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_QUIET_tod)
    !           QUIET-specific TOD object
    ! file:     derived type (hdf_file)
    !           Already open HDF file handle
    ! slabel:   string
    !           Scan label, e.g., "000001/"
    ! detlabels: string (array)
    !           Array of detector labels, e.g., ["27M", "27S"]
    ! scan:     derived class (comm_scan)
    !           Scan object
    !
    ! Returns
    ! ----------
    ! None, but updates scan object
    !
    implicit none
    class(comm_QUIET_tod),               intent(in)    :: self
    type(hdf_file),                      intent(in)    :: file
    character(len=*),                    intent(in)    :: slabel
    character(len=*), dimension(:),      intent(in)    :: detlabels
    class(comm_scan),                    intent(inout) :: scan
  end subroutine read_scan_inst_QUIET

  subroutine initHDF_QUIET(self, chainfile, path)
    ! 
    ! Initializes QUIET-specific TOD parameters from existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_QUIET_tod)
    !           QUIET-specific TOD object
    ! chainfile: derived type (hdf_file)
    !           Already open HDF file handle to existing chainfile
    ! path:   string
    !           HDF path to current dataset, e.g., "000001/tod/030"
    !
    ! Returns
    ! ----------
    ! None
    !
    implicit none
    class(comm_QUIET_tod),               intent(inout)  :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
  end subroutine initHDF_QUIET
  
  subroutine dumpToHDF_QUIET(self, chainfile, path)
    ! 
    ! Writes QUIET-specific TOD parameters to existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_QUIET_tod)
    !           QUIET-specific TOD object
    ! chainfile: derived type (hdf_file)
    !           Already open HDF file handle to existing chainfile
    ! path:   string
    !           HDF path to current dataset, e.g., "000001/tod/030"
    !
    ! Returns
    ! ----------
    ! None
    !
    implicit none
    class(comm_QUIET_tod),               intent(in)     :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
  end subroutine dumpToHDF_QUIET

end module comm_tod_QUIET_mod
