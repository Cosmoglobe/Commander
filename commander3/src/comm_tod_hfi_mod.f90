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
module comm_tod_hfi_mod
  !   Module which contains all the HFI time ordered data processing and routines
  !   for a given frequency band
  !
  !   Main Methods
  !   ------------
  !   constructor(cpar, id_abs, info, tod_type)
  !       Initialization routine that reads in, allocates and associates
  !       all data needed for TOD processing
  !   process_HFI_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
  !       Routine which processes the time ordered data
  !
  use comm_tod_driver_mod
  use comm_tod_cray_mod
  use comm_conviqt_mod
  use comm_tod_crosstalk_mod
  use comm_tod_pixhist_mod
  use comm_tod_adc_binfit_mod
  implicit none

  private
  public comm_hfi_tod

  type, extends(comm_tod) :: comm_hfi_tod
     real(sp) :: f_spin
     integer(i4b), allocatable, dimension(:,:) :: adu_range   ! (ndet,min/max)
     class(comm_crosstalk),    pointer :: xtalk
     type(adc_binfit_pointer), allocatable, dimension(:) :: adc ! (ndet)
     real(sp), allocatable, dimension(:) :: pol_eff ! (ndet)
   contains
     procedure     :: process_tod             => process_hfi_tod
     procedure     :: read_tod_inst           => read_tod_inst_hfi
     procedure     :: read_scan_inst          => read_scan_inst_hfi
     procedure     :: load_instrument_inst    => load_instrument_hfi
     procedure     :: initHDF_inst            => initHDF_hfi
     procedure     :: dumpToHDF_inst          => dumpToHDF_hfi
     procedure     :: construct_corrtemp_inst => construct_corrtemp_hfi
     procedure     :: apply_nonlin_corr_inst  => apply_nonlin_corr_hfi

     procedure, private     :: stitch_hfi_dc_level
     procedure, private     :: hfi_dark_correction
     procedure, private     :: estimate_hfi_4k_lines
     procedure, private     :: deconvolve_rolloff
     procedure, private     :: fill_gaps
     procedure, private     :: sample_adc_and_baselines
     procedure, private     :: compute_adu_range
  end type comm_hfi_tod

  interface comm_hfi_tod
     procedure constructor_hfi
  end interface comm_hfi_tod

interface

  !**************************************************
  !             Constructor
  !**************************************************
  module function constructor_hfi(cpar, id, id_abs, info, tod_type) result(c)
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
    !           Information about the maps for this band, like how the maps are distributed in memory
    ! tod_type: string
    !           Instrument specific tod type
    !
    ! Returns
    ! ----------
    ! constructor: pointer
    !              Pointer that contains all instrument data

    implicit none
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs
    class(comm_mapinfo),     target     :: info
    character(len=128),      intent(in) :: tod_type
    class(comm_hfi_tod),     pointer    :: c

  end function constructor_hfi

  !**************************************************
  !             Driver routine
  !**************************************************
  module subroutine process_hfi_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
    !
    ! Routine that processes the HFI time ordered data.
    ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
    ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms.
    ! Writes maps to disc in fits format
    !
    ! Arguments:
    ! ----------
    ! self:     pointer of comm_HFI_tod class
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
    class(comm_hfi_tod),                      intent(inout) :: self
    character(len=*),                         intent(in)    :: chaindir
    integer(i4b),                             intent(in)    :: chain, iter
    type(planck_rng),                         intent(inout) :: handle
    type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
    real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
    class(comm_map),                          intent(inout) :: map_out      ! Combined output map
    class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms
    type(map_ptr),       dimension(1:,1:),    intent(inout), optional :: map_gain       ! (ndet,1)
  end subroutine process_hfi_tod

  module subroutine load_instrument_hfi(self, instfile, band)
    !
    ! Reads the HFI specific fields from the instrument file
    ! Implements comm_tod_mod::load_instrument_inst
    !
    ! Arguments:
    !
    ! self : comm_HFI_tod
    !    the HFI tod object (this class)
    ! file : hdf_file
    !    the open file handle for the instrument file
    ! band : int
    !    the index of the current detector
    ! 
    ! Returns : None
    implicit none
    class(comm_hfi_tod),                 intent(inout) :: self
    type(hdf_file),                      intent(in)    :: instfile
    integer(i4b),                        intent(in)    :: band
  end subroutine load_instrument_hfi


  module subroutine sample_hfi_baselines(self, tod, scan, handle, subtract_s_tot)

    ! 
    ! Estimates baselines for MODULATED data, separate for odd and even samples
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_scandata)
    !           HFI-specific TOD object
    ! tod:      comm_tod derived type
    !             contains TOD-specific information         
    ! scan:     scan ID
    ! handle:   planck_rng derived type
    !           Healpix definition for random number generation
    !           so that the same sequence can be resumed later on from that same
    !           point
    !           
    !
    ! Returns
    ! ----------
    !   None, but updates tod%scans(scan)%d(:)%baseline  (for odd samples)
    !                     tod%scans(scan)%d(:)%baseline2 (for even samples)
    !                     tod%scans(scan)%d(:)%gain (temporary solution)
    !
    implicit none
    class(comm_scandata),                 intent(in)    :: self
    class(comm_hfi_tod),                  intent(inout) :: tod
    integer(i4b),                         intent(in)    :: scan
    type(planck_rng),                     intent(inout) :: handle
    logical(lgt),                         intent(in), optional :: subtract_s_tot
  end subroutine sample_hfi_baselines

  module subroutine set_modulation_phase(self, tod, scan)
    ! 
    ! Estimates baselines for MODULATED data, separate for odd and even samples
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_scandata)
    !           HFI-specific TOD object
    ! tod:      comm_tod derived type
    !             contains TOD-specific information         
    ! scan:     scan ID
    ! handle:   planck_rng derived type
    !           Healpix definition for random number generation
    !           so that the same sequence can be resumed later on from that same
    !           point
    !           
    !
    ! Returns
    ! ----------
    !   None, but updates tod%scans(scan)%d(:)%baseline1  (for odd samples)
    !                     tod%scans(scan)%d(:)%baseline2 (for even samples)
    !                     tod%scans(scan)%d(:)%gain (temporary solution)
    !
    implicit none
    class(comm_scandata),                 intent(in)    :: self
    class(comm_hfi_tod),                  intent(inout) :: tod
    integer(i4b),                         intent(in)    :: scan
  end subroutine set_modulation_phase
  

  module subroutine demodulate_tod(self, tod, scan)
    ! 
    ! Demodulate HFI TOD
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_scandata)
    !           HFI-specific TOD object
    ! tod:      comm_tod derived type
    !             contains TOD-specific information         
    ! scan:     Scan ID
    !
    ! Returns
    ! ----------
    !   None, but updates self%tod
    !
    implicit none
    class(comm_scandata),                         intent(inout) :: self
    class(comm_hfi_tod),                          intent(in)    :: tod
    integer(i4b),                                 intent(in)    :: scan
  end subroutine demodulate_tod

  
  module subroutine read_tod_inst_hfi(self, file)
    ! 
    ! Reads HFI-specific common fields from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_HFI_tod)
    !           HFI-specific TOD object
    ! file:     derived type (hdf_file)
    !           Already open HDF file handle; only root includes this
    !
    ! Returns
    ! ----------
    ! None, but updates self
    !
    implicit none
    class(comm_HFI_tod),                 intent(inout)          :: self
    type(hdf_file),                      intent(in),   optional :: file
  end subroutine read_tod_inst_hfi
  
  module subroutine read_scan_inst_hfi(self, file, slabel, detlabels, scan)
    ! 
    ! Reads HFI-specific scan information from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_HFI_tod)
    !           HFI-specific TOD object
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
    class(comm_HFI_tod),                 intent(in)    :: self
    type(hdf_file),                      intent(in)    :: file
    character(len=*),                    intent(in)    :: slabel
    character(len=*), dimension(:),      intent(in)    :: detlabels
    class(comm_scan),                    intent(inout) :: scan
  end subroutine read_scan_inst_hfi

  module subroutine initHDF_hfi(self, chainfile, path)
    ! 
    ! Initializes HFI-specific TOD parameters from existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_HFI_tod)
    !           HFI-specific TOD object
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
    class(comm_HFI_tod),                 intent(inout)  :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
  end subroutine initHDF_hfi
  
  module subroutine dumpToHDF_hfi(self, chainfile, path)
    ! 
    ! Writes HFI-specific TOD parameters to existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_HFI_tod)
    !           HFI-specific TOD object
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
    class(comm_hfi_tod),                 intent(in)     :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
  end subroutine dumpToHDF_hfi

  module subroutine construct_corrtemp_hfi(self, scan, pix, psi, s, det)
    !  Construct an LFI instrument-specific correction template; for now contains 1Hz template only
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
    class(comm_hfi_tod),                   intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    integer(i4b),        dimension(:,:),   intent(in)    :: pix, psi
    real(sp),            dimension(:,:),   intent(out)   :: s
    integer(i4b),                          intent(in), optional :: det
  end subroutine construct_corrtemp_hfi

  module subroutine apply_nonlin_corr_hfi(self, scan, sd, skip_nonlin, handle, det)
    !  Construct and apply HFI instrument-specific non-linear corrections
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  scan: int
    !       scan number
    !
    !  sd: comm_scandata
    !       structure holding the scan data
    !
    !  Returns:
    !  --------
    !  s:   real (sp)
    !       output template timestream
    implicit none
    class(comm_hfi_tod),                   intent(inout) :: self
    integer(i4b),                          intent(in)    :: scan
    class(comm_scandata),                  intent(inout) :: sd
    integer(i4b),                          intent(in)    :: skip_nonlin
    type(planck_rng),            optional, intent(inout) :: handle
    integer(i4b),                optional, intent(in)    :: det
  end subroutine apply_nonlin_corr_hfi

  module subroutine stitch_hfi_dc_level(self, scan, sd)
    !  Construct and apply HFI instrument-specific non-linear corrections
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  scan: int
    !       scan number
    !  sd: comm_scandata object
    !       structure holding the data for each scan
    !
    implicit none
    class(comm_hfi_tod),                   intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    class(comm_scandata),                  intent(inout) :: sd
  end subroutine stitch_hfi_dc_level

  module subroutine hfi_dark_correction(self, scan, sd)
    !  Construct and apply HFI instrument-specific corrections
    !  from dark bolometer timestreams
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  scan: int
    !       scan number
    !  sd: comm_scandata object
    !       structure holding the data for each scan
    !
    implicit none
    class(comm_hfi_tod),                   intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    class(comm_scandata),                  intent(inout) :: sd
  end subroutine hfi_dark_correction

  module subroutine sample_adc_and_baselines(self, handle, det, map_sky, procmask)
    !  Sample ADC parameters
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    implicit none
    class(comm_hfi_tod),                 intent(inout) :: self
    type(planck_rng),                    intent(inout) :: handle
    integer(i4b),                        intent(in)    :: det
    real(sp),          dimension(1:,1:), intent(in)    :: map_sky
    real(sp),          dimension(0:),    intent(in)    :: procmask
  end subroutine sample_adc_and_baselines

  module subroutine estimate_hfi_4k_lines(self, scan, sd)
    !  Construct and apply HFI instrument-specific corrections
    !  from 4k lines
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  scan: int
    !       scan number
    !  sd: comm_scandata object
    !       structure holding the data for each scan
    !
    implicit none
    class(comm_hfi_tod),                   intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    class(comm_scandata),                  intent(inout) :: sd
  end subroutine estimate_hfi_4k_lines

  module subroutine deconvolve_rolloff(self, tod, scan, i_det, s_sub, mask, flag, handle, ps_output)
    ! Deconvolves high frequency rolloff in noise spectrum
    !
    ! Arguments:
    ! ----------
    ! self: comm_tod object
    !
    ! tod: real(sp) array
    !      residual tod
    ! scan: int
    !       scan number
    ! i_det: int
    !        detector id
    ! s_sub: real(sp) array
    !        sky signal template
    ! mask: real(sp) array
    !       mask array
    ! flag:  integer array
    !        quality flags
    ! handle: planck_rng
    !         rng handle
    implicit none
    class(comm_hfi_tod),                       intent(inout) :: self
    real(sp),                   dimension(1:), intent(inout) :: tod
    integer(i4b),                              intent(in)    :: scan, i_det
    real(sp),                   dimension(1:), intent(in)    :: s_sub
    real(sp),                   dimension(1:), intent(in)    :: mask
    integer(i4b),               dimension(1:), intent(inout) :: flag
    type(planck_rng),                          intent(inout) :: handle
    character(len=*),                optional, intent(in)    :: ps_output
  end subroutine deconvolve_rolloff

  module subroutine fill_gaps(self, tod, handle, scan, i_det, mask, s_sub, pix, nomono, dospike, ps_output, filling)
    ! Fill gaps in flagged samples
    !
    ! Arguments:
    ! ----------
    ! self: comm_tod object
    !
    ! tod: real(sp) array
    !      residual tod
    ! handle: planck_rng
    !         rng handle
    ! scan: int
    !       scan number
    ! i_det: int
    !        detector id
    ! s_sub: real(sp) array
    !        sky signal template
    ! mask:  real(sp) array
    !        mask
    ! pix:   int array
    !         
    ! nomono: logical
    !         option to remove monopole
    ! dospike: logical
    !          option to flag spikes
    ! ps_output: string
    !            filename to output corrected power spectrum
    ! filling: string
    !          filling mod (white noise, chunks, zero)
    implicit none
    class(comm_hfi_tod),                      intent(in)    :: self
    real(sp),               dimension(1:),    intent(inout) :: tod
    type(planck_rng),                         intent(inout) :: handle
    integer(i4b),                             intent(in)    :: scan, i_det
    real(sp),               dimension(1:),    intent(in)    :: mask, s_sub
    integer(i4b), optional, dimension(1:,1:), intent(in)    :: pix
    logical(lgt),                   optional, intent(in)    :: nomono
    logical(lgt),                   optional, intent(in)    :: dospike
    character(len=*),               optional, intent(in)    :: ps_output
    character(len=*),               optional, intent(in)    :: filling
  end subroutine fill_gaps

  module subroutine compute_adu_range(self)
    ! 
    ! Computes ADU range over all scans and unmasked samples; for ADC correction
    ! Must be called after dynamic mask definition
    !
    ! Arguments:
    ! ----------
    ! self:      comm_tod derived type
    !             contains TOD-specific information         
    ! Returns
    ! ----------
    !   None, but updates tod%adu_range
    !
    implicit none
    class(comm_hfi_tod),                  intent(inout) :: self
  end subroutine compute_adu_range

  
end interface

end module comm_tod_hfi_mod
