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
module comm_tod_DIRBE_mod
   !   Module which contains all the LiteBIRD time ordered data processing and routines
   !   for a given frequency band
   !
   !   Main Methods
   !   ------------
   !   constructor(cpar, id_abs, info, tod_type)
   !       Initialization routine that reads in, allocates and associates 
   !       all data needed for TOD processing
   !   process_DIRBE_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
   !       Routine which processes the time ordered data
   !
   use comm_tod_driver_mod
   implicit none

   private
   public comm_dirbe_tod

   type, extends(comm_tod) :: comm_dirbe_tod
   contains
      procedure     :: process_tod          => process_DIRBE_tod
      procedure     :: dumpToHDF_inst       => dumpToHDF_DIRBE
   end type comm_dirbe_tod

   interface comm_dirbe_tod
      procedure constructor_dirbe
   end interface comm_dirbe_tod

contains
   !**************************************************
   !             Constructor
   !**************************************************
   function constructor_dirbe(cpar, id, id_abs, info, tod_type) result(c)
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
      class(comm_dirbe_tod),      pointer    :: c

      integer(i4b) :: i, j, nside_beam, lmax_beam, nmaps_beam, ierr
      logical(lgt) :: pol_beam


      ! Allocate object
      allocate(c)

    end function constructor_dirbe


    subroutine dumpToHDF_DIRBE(self, chainfile, path)
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


      if (self%myid == 0) then
         if (self%map_solar_allocated == .true.) then
           call write_hdf(chainfile, trim(adjustl(path))//'map_solar',  self%map_solar)
         end if
      end if

    end subroutine dumpToHDF_DIRBE

   !**************************************************
   !             Driver routine
   !**************************************************
   subroutine process_DIRBE_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
      ! 
      ! Routine that processes the DIRBE Calibrated Individual Observations. 
      ! Samples absolute and relative bandpass, gain and correlated noise in time domain, 
      ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms. 
      ! Writes maps to disc in fits format
      ! 
      ! Arguments:
      ! ----------
      ! self:     pointer of comm_dirbe_tod class
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
   end subroutine process_DIRBE_tod   


end module comm_tod_DIRBE_mod
