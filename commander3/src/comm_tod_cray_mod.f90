!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University
! of Oslo.
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
module comm_tod_cray_mod
  use comm_param_mod
  implicit none

  private
  public constructor_cray, comm_tod_cray, cray_ptr

  ! Cosmic ray types
  integer(i4b), parameter :: CR_SHORT = 1
  integer(i4b), parameter :: CR_LONG  = 2
  integer(i4b), parameter :: CR_SLOW  = 3
  
  ! Class for single cosmic ray event
  type :: cray_event
     integer(i4b) :: type    ! CR type, defined above
     integer(i4b) :: nsamp   ! Maximum TOD segment length
     real(dp), allocatable, dimension(:) :: p   ! Free shape parameters
     real(sp), allocatable, dimension(:) :: T   ! CR shape template; normalized to 1
   contains
     procedure :: build_tod_template
  end type cray_event
  
  type :: comm_tod_cray
     integer(i4b)      :: n        ! Number of events in current scan
     integer(i4b)      :: ntypes   ! Number of active base types
     character(len=64) :: freq     ! Frequency label
     integer(i4b)      :: det      ! Detector ID
     integer(i4b)      :: scanid   ! Absolute scan ID
     type(cray_event), pointer,     dimension(:) :: T_base       ! Default templates
     type(cray_event), pointer,     dimension(:) :: event_list   ! List of individual event objects; 
     real(sp),         allocatable, dimension(:) :: amp          ! CR amplitudes
     integer(i4b),     allocatable, dimension(:) :: event_type   ! >0 = pointers to base events; < 0 = fit individually
     integer(i4b),     allocatable, dimension(:) :: event_length ! TOD truncation length of each event
   contains
     procedure :: update_base_templates
     procedure :: detect_events
     procedure :: select_cr_type
     procedure :: fit_local_shape_params
     procedure :: fit_amps
     procedure :: generate => generate_cray_correction
     procedure :: params2cache                            ! IO communication routine
  end type comm_tod_cray

  interface comm_tod_cray 
    procedure constructor_cray
  end interface comm_tod_cray

  type cray_ptr
    class(cray_event), pointer :: p => null()
 end type cray_ptr

contains

    ! Constructor
    function constructor_cray(freq, det, scanid, active_cr_types) result(c)
      implicit none
      character(len=*),                   intent(in) :: freq             ! Frequency label
      integer(i4b),                       intent(in) :: det              ! Detector ID
      integer(i4b),                       intent(in) :: scanid             ! Scan ID (absolute)
      integer(i4b),         dimension(:), intent(in) :: active_cr_types  ! List of active types
      class(comm_tod_cray), pointer          :: c

      allocate(c)
      c%freq    = freq
      c%det     = det
      c%scanid  = scanid
      c%ntypes  = size(active_cr_types)
      
    end function constructor_cray

    ! Routine for constructing TOD shape template for current event/base model
    subroutine build_tod_template(self)
      implicit none
      class(cray_event),              intent(in)  :: self
    end subroutine build_tod_template
    
    ! builds a set of templates for cosmic ray response from the TOD
    subroutine update_base_templates(self)
      implicit none
      class(comm_tod_cray),                          intent(inout) :: self

    end subroutine update_base_templates

    ! Search for events; store list of impact samples; allocate event arrays; initialize defaults; define local length
    ! This routine should only be called (at most) once per main run
    subroutine detect_events(self, tod)
      implicit none
      class(comm_tod_cray),                      intent(inout) :: self
      real(sp), dimension(:),                    intent(in)    :: tod

      integer(i4b) :: i
      
      self%n = 0
      do i = 1, self%n
         self%event_length(i) = 0 ! Number of samples in current template; scales with estimated amplitude
      end do

    end subroutine detect_events

    ! Search for events; store list of impact samples; allocate event arrays; initialize defaults
    subroutine select_cr_type(self, tod)
      implicit none
      class(comm_tod_cray),                      intent(inout) :: self
      real(sp), dimension(:),                    intent(in)    :: tod

      integer(i4b) :: i

      do i = 1, self%n
         self%event_type(i) = 0
      end do

    end subroutine select_cr_type

    
    ! fits the constructed templates to the cosmic rays in the timestreams
    subroutine fit_local_shape_params(self, tod)
      implicit none
      class(comm_tod_cray),                      intent(inout) :: self
      real(sp), dimension(:),                    intent(in)    :: tod

      integer(i4b) :: i
               
      ! Fit parameters for individual events
      do i = 1, self%n
         if (self%event_type(i) > 0) cycle  ! Skip elements modelled by global templates
      end do
         
    end subroutine fit_local_shape_params

    ! fits the constructed templates to the cosmic rays in the timestreams
    subroutine fit_amps(self, tod)
      implicit none
      class(comm_tod_cray),                      intent(inout) :: self
      real(sp), dimension(:),                    intent(in)    :: tod
      self%amp = 0.
    end subroutine fit_amps
    
    ! Routine for constructing total CR model for a given scan
    subroutine generate_cray_correction(self, s_cray)
      implicit none
      class(comm_tod_cray),               intent(in)  :: self
      real(sp),             dimension(:), intent(out) :: s_cray

      integer(i4b) :: i

      s_cray = 0.
      do i = 1, self%n
         ! Add contribution for current event to s_cray
      end do
      
    end subroutine generate_cray_correction

    ! IO communication routine
    subroutine params2cache(self, cache)
      implicit none
      class(comm_tod_cray),                         intent(in)    :: self
      real(sp),          allocatable, dimension(:), intent(inout) :: cache

      if (allocated(cache)) then
         ! Initialize model parameters from cache
      else
         ! Store model parameters in cache
      end if
      
    end subroutine params2cache

end module comm_tod_cray_mod
