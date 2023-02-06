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
module comm_tod_SPIDER_mod
   !   Module which contains all the SPIDER time ordered data processing and routines
   !   for a given frequency band
   !
   !   Main Methods
   !   ------------
   !   constructor(cpar, id_abs, info, tod_type)
   !       Initialization routine that reads in, allocates and associates
   !       all data needed for TOD processing
   !   process_SPIDER_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
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
   use comm_tod_jump_mod
   implicit none
 
   private
   public comm_SPIDER_tod
 
   type, extends(comm_tod) :: comm_SPIDER_tod
    contains
      procedure     :: process_tod        => process_SPIDER_tod
      procedure     :: read_tod_inst      => read_tod_inst_SPIDER
      procedure     :: read_scan_inst     => read_scan_inst_SPIDER
      procedure     :: initHDF_inst       => initHDF_SPIDER
      procedure     :: dumpToHDF_inst     => dumpToHDF_SPIDER
   end type comm_SPIDER_tod
 
   interface comm_SPIDER_tod
      procedure constructor
   end interface comm_SPIDER_tod
 
 interface
 
   !**************************************************
   !             Constructor
   !**************************************************
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
     integer(i4b),            intent(in) :: id_abs
     class(comm_mapinfo),     target     :: info
     character(len=128),      intent(in) :: tod_type
     class(comm_SPIDER_tod),  pointer    :: constructor 
   end function constructor
 
   !**************************************************
   !             Driver routine
   !**************************************************
   module subroutine process_SPIDER_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
     !
     ! Routine that processes the SPIDER time ordered data.
     ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
     ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms.
     ! Writes maps to disc in fits format
     !
     ! Arguments:
     ! ----------
     ! self:     pointer of comm_SPIDER_tod class
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
     class(comm_SPIDER_tod),                   intent(inout) :: self
     character(len=*),                         intent(in)    :: chaindir
     integer(i4b),                             intent(in)    :: chain, iter
     type(planck_rng),                         intent(inout) :: handle
     type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
     real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
     class(comm_map),                          intent(inout) :: map_out      ! Combined output map
     class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms 
     type(map_ptr),     dimension(:,:),   intent(inout), optional :: map_gain

   end subroutine process_SPIDER_tod
 
   
   module subroutine read_tod_inst_SPIDER(self, file)
     ! 
     ! Reads SPIDER-specific common fields from TOD fileset
     ! 
     ! Arguments:
     ! ----------
     ! self:     derived class (comm_SPIDER_tod)
     !           SPIDER-specific TOD object
     ! file:     derived type (hdf_file)
     !           Already open HDF file handle; only root includes this
     !
     ! Returns
     ! ----------
     ! None, but updates self
     !
     implicit none
     class(comm_SPIDER_tod),                 intent(inout)          :: self
     type(hdf_file),                      intent(in),   optional :: file
   end subroutine read_tod_inst_SPIDER
   
   module subroutine read_scan_inst_SPIDER(self, file, slabel, detlabels, scan)
     ! 
     ! Reads SPIDER-specific scan information from TOD fileset
     ! 
     ! Arguments:
     ! ----------
     ! self:     derived class (comm_SPIDER_tod)
     !           SPIDER-specific TOD object
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
     class(comm_SPIDER_tod),              intent(in)    :: self
     type(hdf_file),                      intent(in)    :: file
     character(len=*),                    intent(in)    :: slabel
     character(len=*), dimension(:),      intent(in)    :: detlabels
     class(comm_scan),                    intent(inout) :: scan
   end subroutine read_scan_inst_SPIDER
 
   module subroutine initHDF_SPIDER(self, chainfile, path)
     ! 
     ! Initializes SPIDER-specific TOD parameters from existing chain file
     ! 
     ! Arguments:
     ! ----------
     ! self:     derived class (comm_SPIDER_tod)
     !           SPIDER-specific TOD object
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
     class(comm_SPIDER_tod),              intent(inout)  :: self
     type(hdf_file),                      intent(in)     :: chainfile
     character(len=*),                    intent(in)     :: path
   end subroutine initHDF_SPIDER
   
   module subroutine dumpToHDF_SPIDER(self, chainfile, path)
     ! 
     ! Writes SPIDER-specific TOD parameters to existing chain file
     ! 
     ! Arguments:
     ! ----------
     ! self:     derived class (comm_SPIDER_tod)
     !           SPIDER-specific TOD object
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
     class(comm_SPIDER_tod),              intent(in)     :: self
     type(hdf_file),                      intent(in)     :: chainfile
     character(len=*),                    intent(in)     :: path
   end subroutine dumpToHDF_SPIDER

   module subroutine write2file(filename, iter, param)
      implicit none
   
      character(len=*), intent(in)         :: filename
      real(dp), intent(in)                 :: param
      integer(i4b), intent(in)             :: iter
    end subroutine write2file
  end interface
 end module comm_tod_SPIDER_mod

 
