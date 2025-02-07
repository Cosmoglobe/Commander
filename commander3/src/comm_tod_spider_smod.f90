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
submodule (comm_tod_spider_mod) comm_tod_SPIDER_mod
contains
 
   !**************************************************
   !             Constructor
   !**************************************************
   module function constructor_spider(cpar, id, id_abs, info, tod_type) result(c)
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
     class(comm_SPIDER_tod),  pointer    :: c
 
     integer(i4b) :: i, nside_beam, lmax_beam, nmaps_beam, ierr
     logical(lgt) :: pol_beam
     character(len=8)     :: det, det_partner
     character(len=2)     :: row_str, row_partner_str
     integer(i4b)         :: row_int, row_partner_int
     integer(i4b)         :: det_index(1)


     ! Allocate object
     allocate(c)
 

 
   end function constructor_spider
 
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
 

     real(dp)            :: t1, t2
     integer(i4b)        :: i, j, k, l, ierr, ndelta, nside, npix, nmaps
     logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, output_scanlist
     type(comm_binmap)   :: binmap
     type(comm_scandata) :: sd
     character(len=4)    :: ctext, myid_text
     character(len=6)    :: samptext, scantext
     character(len=512)  :: prefix, postfix, prefix4D, filename
     character(len=512), allocatable, dimension(:) :: slist
     real(sp), allocatable, dimension(:)       :: procmask, procmask2, sigma0
     real(sp), allocatable, dimension(:,:)     :: s_buf
     real(sp), allocatable, dimension(:,:,:)   :: d_calib
     real(sp), allocatable, dimension(:,:,:,:) :: map_sky, m_gain
     real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf

     real(sp),     allocatable, dimension(:,:)   :: s_jump, tod_gapfill
     real(sp),     allocatable, dimension(:,:,:) :: jump_calib
     integer(i4b), allocatable, dimension(:,:)   :: jumps, offset_range, jumpflag_range
     real(sp),     allocatable, dimension(:)     :: offset_level
     type(comm_binmap)                           :: jump_map
     character(len=4)                            :: it_label
     logical(lgt)                                :: debug
     real(sp),    allocatable, dimension(:)      :: test_array




     call int2string(iter, ctext)
     call int2string(iter, it_label)
 
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
   
      integer(i4b)                           :: unit, io_error
      logical                                :: existing
   
      unit = 22
   
      inquire(file=trim(filename),exist=existing)
      if (existing) then
         open(unit,file=trim(filename),status='old',position='append',action='write',iostat=io_error)
      else
         open(unit,file=trim(filename),status='replace',action='write',iostat=io_error)
      end if
   
      write(unit,*), iter, param
   
      close(unit)
    end subroutine write2file
 
 end submodule comm_tod_SPIDER_mod
 
