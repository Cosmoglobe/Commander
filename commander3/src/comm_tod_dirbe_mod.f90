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
module comm_tod_dirbe_mod
  !
  !   Module which contains all the LFI time ordered data processing and routines
  !   for a given frequency band
  !
  !   Main Methods
  !   ------------
  !   constructor(cpar, id_abs, info, tod_type)
  !       Initialization routine that reads in, allocates and associates
  !       all data needed for TOD processing
  !   process_LFI_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
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
  use comm_shared_arr_mod
  use spline_1D_mod
  use comm_4D_map_mod
  use comm_tod_driver_mod
  use comm_utils
  use comm_tod_adc_mod
  implicit none

  private
  public comm_dirbe_tod

  type, extends(comm_tod) :: comm_dirbe_tod
     integer(i4b) :: nbin_spike
     integer(i4b) :: nbin_adc
     logical(lgt) :: use_dpc_adc
     logical(lgt) :: use_dpc_gain_modulation
     real(dp),             allocatable, dimension(:)       :: mb_eff
     real(dp),             allocatable, dimension(:,:)     :: diode_weights
     type(spline_type),    allocatable, dimension(:,:)     :: ref_splint ! ndet, ndiode/2
     type(adc_pointer),    allocatable, dimension(:,:)     :: adc_corrections ! ndet, n_diode
     real(dp),             allocatable, dimension(:,:)     :: spike_templates ! nbin, ndet
     real(dp),             allocatable, dimension(:,:)     :: spike_amplitude ! nscan, ndet
     real(dp),             allocatable, dimension(:,:,:)   :: R               ! nscan, ndet, ndiode/2
     type(double_pointer), allocatable, dimension(:)       :: gmf_splits      ! ndet
     logical(lgt),         allocatable, dimension(:,:)     :: apply_adc       ! ndet, n_diode
     character(len=10),    allocatable, dimension(:,:)     :: adc_mode
   contains
     procedure     :: process_tod             => process_dirbe_tod
     procedure     :: diode2tod_inst          => diode2tod_dirbe
     procedure     :: load_instrument_inst    => load_instrument_dirbe
     procedure     :: dumpToHDF_inst          => dumpToHDF_dirbe
     procedure     :: construct_corrtemp_inst => construct_corrtemp_dirbe
     procedure     :: initHDF_inst            => initHDF_dirbe
     procedure     :: remove_fixed_scans      => remove_fixed_scans_dirbe
     procedure     :: filter_reference_load
     procedure     :: compute_ref_load_filter
     procedure     :: get_nsmooth
     procedure     :: get_freq_bins
     procedure     :: preprocess_L1_to_L2
  end type comm_dirbe_tod

  interface comm_dirbe_tod
     procedure constructor
  end interface comm_dirbe_tod

  type double_pointer
    real(dp), pointer, dimension(:) :: p => null() 
  end type double_pointer

interface

  !**************************************************
  !             Constructor
  !**************************************************
  module function constructor(handle, cpar, id_abs, info, tod_type) result(res)
    !
    ! Constructor function that gathers all the instrument parameters in a pointer
    ! and constructs the objects
    !
    ! Arguments:
    ! ----------
    ! handle:   type(planck_rng)
    !           Healpix random number type
    ! cpar:     derived type
    !           Object containing parameters from the parameterfile.
    ! id_abs:   integer
    !           The index of the current band within the parameters, related to cpar
    ! info:     map_info structure
    !           Information about the maps for this band, like how the maps are distributed in memory
    ! tod_type: string
    !           Instrument specific tod type
    !
    ! Returns
    ! ----------
    ! constructor: pointer
    !              Pointer that contains all instrument data
    !
    implicit none
    type(planck_rng),          intent(inout) :: handle
    type(comm_params),         intent(in)    :: cpar
    integer(i4b),              intent(in)    :: id_abs
    class(comm_mapinfo),       target        :: info
    character(len=128),        intent(in)    :: tod_type
    class(comm_dirbe_tod),     pointer       :: res
  end function constructor

  !**************************************************
  !             Driver routine
  !**************************************************
  module subroutine process_dirbe_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
    !
    ! Routine that processes the DIRBE time ordered data.
    ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
    ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms.
    ! Writes maps to disc in fits format
    !
    ! Arguments:
    ! ----------
    ! self:     pointer of comm_DIRBE_tod class
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

    implicit none
    class(comm_dirbe_tod),                    intent(inout) :: self
    character(len=*),                         intent(in)    :: chaindir
    integer(i4b),                             intent(in)    :: chain, iter
    type(planck_rng),                         intent(inout) :: handle
    type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
    real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
    class(comm_map),                          intent(inout) :: map_out      ! Combined output map
    class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms
    type(map_ptr),       dimension(1:,1:),    intent(inout), optional :: map_gain       ! (ndet,1)
  end subroutine process_dirbe_tod
  
  
  module subroutine load_instrument_dirbe(self, instfile, band)
    !
    ! Reads the DIRBE specific fields from the instrument file
    ! Implements comm_tod_mod::load_instrument_inst
    !
    ! Arguments:
    !
    ! self : comm_DIRBE_tod
    !    the DIRBE tod object (this class)
    ! file : hdf_file
    !    the open file handle for the instrument file
    ! band : int
    !    the index of the current detector
    ! 
    ! Returns : None
    implicit none
    class(comm_dirbe_tod),               intent(inout) :: self
    type(hdf_file),                      intent(in)    :: instfile
    integer(i4b),                        intent(in)    :: band
  end subroutine load_instrument_dirbe
 
  module subroutine initHDF_dirbe(self, chainfile, path)
    ! 
    ! Initializes instrument-specific TOD parameters from existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
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
    class(comm_dirbe_tod),               intent(inout)  :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
  end subroutine initHDF_dirbe

 
  module subroutine diode2tod_dirbe(self, scan, map_sky, procmask, tod)
    ! 
    ! Generates detector-coadded TOD from low-level diode data
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
    ! scan:     int
    !           Scan ID number
    ! procmask: array of sp
    !           processing mask that cuts out the galaxy
    !
    ! Returns
    ! ----------
    ! tod:      ntod x ndet sp array
    !           Output detector TOD generated from raw diode data
    !
    implicit none
    class(comm_dirbe_tod),                     intent(inout) :: self
    integer(i4b),                              intent(in)    :: scan
    real(sp),          dimension(0:,1:,1:,1:), intent(in)    :: map_sky
    real(sp),          dimension(0:),          intent(in)    :: procmask
    real(sp),          dimension(:,:),         intent(out)   :: tod
  end subroutine diode2tod_dirbe

  module function get_nsmooth(self)
    implicit none
    class(comm_dirbe_tod),  intent(in) :: self
    integer(i4b)                       :: get_nsmooth  
  end function get_nsmooth

  module subroutine compute_ref_load_filter(self, data_in, binned_out, nu_out, err)
    ! 
    ! Computes the binned weiner filter for the reference load
    !
    ! Arguments:
    ! ----------
    ! 
    ! self:     comm_tod_DIRBE object
    !           TOD processing class
    ! data_in:  float array (ntod, ndiode)
    !           input diode timestreams
    !
    ! Returns:
    ! --------
    !
    ! binned_out : float array
    !              array of filter transfer function for ref load
    ! nu_out     : float_array
    !              frequencies that index binned_out
    ! err        : error flag; 0 if OK, 1 if no data
    implicit none
    class(comm_dirbe_tod),        intent(in)    :: self
    real(sp),     dimension(:,:), intent(in)    :: data_in
    real(dp),     dimension(:,:), intent(inout) :: binned_out
    real(dp),     dimension(:),   intent(in)    :: nu_out
    integer(i4b),                 intent(out)   :: err
  end subroutine compute_ref_load_filter

  module subroutine get_freq_bins(self, freqs)
    implicit none
    class(comm_dirbe_tod),  intent(in)    :: self
    real(dp), dimension(:), intent(inout) :: freqs
  end subroutine get_freq_bins


  module subroutine filter_reference_load(self, det, data)
    class(comm_dirbe_tod),             intent(in)      :: self
    integer(i4b),                      intent(in)      :: det
    real(sp), dimension(:,:),          intent(inout)   :: data
  end subroutine filter_reference_load

  module subroutine dumpToHDF_dirbe(self, chainfile, path)
    ! 
    ! Writes instrument-specific TOD parameters to existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
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
    class(comm_dirbe_tod),               intent(in)     :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
  end subroutine dumpToHDF_dirbe

  module subroutine sample_1Hz_spikes(tod, handle, map_sky, procmask, procmask2)
    !   Sample DIRBE specific 1Hz spikes shapes and amplitudes
    !
    !   Arguments:
    !   ----------
    !   tod:      comm_tod derived type
    !             contains TOD-specific information
    !   handle:   planck_rng derived type
    !             Healpix definition for random number generation
    !             so that the same sequence can be resumed later on from that same point
    !   map_sky:
    implicit none
    class(comm_dirbe_tod),                        intent(inout) :: tod
    type(planck_rng),                             intent(inout) :: handle
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2
  end subroutine sample_1Hz_spikes

  module subroutine construct_corrtemp_dirbe(self, scan, pix, psi, s)
    !  Construct an DIRBE instrument-specific correction template; for now contains 1Hz template only
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  scan: int
    !       scan number
    !  pix: int
    !       index for pixel
    !  psi: int
    !       integer label for polarization angle
    !
    !  Returns:
    !  --------
    !  s:   real (sp)
    !       output template timestream
    implicit none
    class(comm_dirbe_tod),                 intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    integer(i4b),        dimension(:,:),   intent(in)    :: pix, psi
    real(sp),            dimension(:,:),   intent(out)   :: s
  end subroutine construct_corrtemp_dirbe


  module subroutine preprocess_L1_to_L2(self, map_sky, procmask)
    implicit none
    class(comm_dirbe_tod),                        intent(inout) :: self
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:),           intent(in)    :: procmask
  end subroutine preprocess_L1_to_L2

  module subroutine remove_fixed_scans_dirbe(self)
    ! 
    ! Sets accept = .false. for known bad scans
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
    !
    ! Returns
    ! ----------
    ! None
    !
    implicit none
    class(comm_dirbe_tod), intent(inout)  :: self
  end subroutine remove_fixed_scans_dirbe

end interface

end module comm_tod_dirbe_mod

