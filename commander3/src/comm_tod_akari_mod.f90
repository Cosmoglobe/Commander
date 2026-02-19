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
module comm_tod_akari_mod
   !   Module which contains all the LiteBIRD time ordered data processing and routines
   !   for a given frequency band
   !
   !   Main Methods
   !   ------------
   !   constructor(cpar, id_abs, info, tod_type)
   !       Initialization routine that reads in, allocates and associates 
   !       all data needed for TOD processing
   !   process_akari_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
   !       Routine which processes the time ordered data
   !
  use comm_tod_driver_mod
  use comm_tod_pixhist_mod
  use comm_tod_mapmaking_mod
  use comm_tod_cgmap_mod
   implicit none

   private
   public comm_akari_tod

   type, extends(comm_tod) :: comm_akari_tod
      integer(i4b)                                  :: ntempl       ! Number of tod correction templates
      integer(i4b), allocatable, dimension(:,:,:)   :: nsamp_templ  ! length of each template [ntempl, ndet, nscan]
      real(dp),     allocatable, dimension(:,:,:,:) :: tod_correction_templ  ! [nsamp,ntempl,ndet,nscan]
      class(comm_dynmask), pointer :: dynmask
   contains
     procedure     :: process_tod             => process_akari_tod
     procedure     :: apply_fast_flags_inst   => apply_fast_flags_akari
     procedure     :: construct_corrtemp_inst => construct_corrtemp_akari
     procedure     :: sample_ramp, init_sample_ramp
   end type comm_akari_tod

   interface comm_akari_tod
      procedure constructor_akari
   end interface comm_akari_tod

interface
      
   !**************************************************
   !             Constructor
   !**************************************************
   module function constructor_akari(cpar, id, id_abs, info, tod_type) result(c)
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
      type(comm_params),       intent(in) :: cpar          !comm_param structure, list of all the input parameters
      integer(i4b),            intent(in) :: id, id_abs        !index of the current band within the parameters 
      class(comm_mapinfo),     target     :: info
      character(len=128),      intent(in) :: tod_type      !
      class(comm_akari_tod),      pointer    :: c
    end function constructor_akari

   !**************************************************
   !             Driver routine
   !**************************************************
   module subroutine process_akari_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
      ! 
      ! Routine that processes the akari Calibrated Individual Observations. 
      ! Samples absolute and relative bandpass, gain and correlated noise in time domain, 
      ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms. 
      ! Writes maps to disc in fits format
      ! 
      ! Arguments:
      ! ----------
      ! self:     pointer of comm_akari_tod class
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
      class(comm_akari_tod),                    intent(inout) :: self
      character(len=*),                         intent(in)    :: chaindir
      integer(i4b),                             intent(in)    :: chain, iter
      type(planck_rng),                         intent(inout) :: handle
      type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
      real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
      class(comm_map),                          intent(inout) :: map_out      ! Combined output map
      class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms
      type(map_ptr),       dimension(1:),       intent(inout), optional :: map_gain       ! (ndet,1)
   end subroutine process_akari_tod   
   
   module subroutine apply_fast_flags_akari(self, sd)
     !  Apply fast flags to sd%flag; should only depend on time, pix or flag arrays, not TOD
     !  Expensive operations should instead be added to the dynamic mask
     !
     !  Arguments:
     !  ----------
     !  self: comm_tod object
     !
     implicit none
     class(comm_akari_tod),                 intent(inout)    :: self
     class(comm_scandata),                  intent(inout)    :: sd
   end subroutine apply_fast_flags_akari

   
   module subroutine construct_corrtemp_akari(self, sd, det)
     ! construct AKARI instrument-specific correction template from ramp template
     ! puts this into sd%s_inst
    implicit none
    class(comm_akari_tod), intent(in)             :: self
    class(comm_scandata),  intent(inout)          :: sd
    integer(i4b),          intent(in),   optional :: det
  end subroutine construct_corrtemp_akari

  module subroutine init_sample_ramp(self, sd)
    ! scans tod flags to find length of tod segment between ramp flags (= length of correction template)
    ! puts this into self%nsamp_templ
    implicit none
    class(comm_akari_tod), intent(inout)       ::self
    class(comm_scandata),  intent(in)          :: sd

  end subroutine init_sample_ramp

  module subroutine sample_ramp(self, sd)
    ! find template of baseline beween to ramp events, by averaging over all events in a scan
    ! put this into self%tod_correction_templ
    implicit none
    class(comm_akari_tod), intent(inout)       ::self
    class(comm_scandata),  intent(in)          :: sd

  end subroutine sample_ramp

   
end interface
   
end module comm_tod_akari_mod
