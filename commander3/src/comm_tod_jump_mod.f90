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
module comm_tod_jump_mod
  use comm_param_mod
  implicit none

  private
  public comm_tod_jump

  type :: comm_tod_jump
    integer(i4b) :: njump
    integer(dp), allocatable, dimension(:)     :: jump_pos
    integer(dp), allocatable, dimension(:,:)   :: mask       ! (njump,2), start/stop
    real(sp),    allocatable, dimension(:)     :: jump_amp
    integer(i4b) :: maxseg
  contains
    procedure :: find_jumps
    procedure :: sample_jump_amplitudes
    procedure :: get_jump_tod
 end type comm_tod_jump

  interface comm_tod_jump 
    procedure constructor_jump
  end interface comm_tod_jump

contains

    ! Constructor
    function constructor_jump() result(c)
      implicit none
      class(comm_tod_jump), pointer           :: c

      allocate(c)

    end function constructor_jump

    subroutine find_jumps(self, tod, signal, flag, samprate)
      implicit none
      
      ! Detects jumps in the residual and stores their positions and amplitudes
      class(comm_tod_jump),                intent(inout) :: self
      real(sp),              dimension(:), intent(inout) :: tod
      real(sp),              dimension(:), intent(in)    :: signal
      integer(i4b),          dimension(:), intent(inout) :: flag
      real(dp),                            intent(in)    :: samprate !KG NOTE: make sure to include this in the call for find_jumps.

      ! KG: these follow the param names of Simons prototype
      ! TODO: check that these are the names we want
      real(sp), allocatable, dimension(:) :: resdidual
      real(dp), allocatable, dimension(:) :: detection_work, bmed, filter_score, step_diff
      logical,  allocatable, dimension(:) :: good
      integer(i4b)                        :: n, nblock, iteration, scale_idx, block_idx
      integer(i4b)                        :: sample_idx, window_blocks
      real(dp)                            :: best_score, score, noise_scale, step_amp

      ! Perform matched filer

      ! Count number of detections
      self%njump = 0

      ! Store jump inidices in self%jump_pos

      ! Identidy bad range; store in self%mask

      ! Correct tod

      ! Update flag array
      
    end subroutine find_jumps

    subroutine sample_jump_amplitudes(self, tod)
      implicit none

      ! Refits the jump ampl. from the residual levels around each jump
      class(comm_tod_jump),               intent(inout) :: self
      real(sp),             dimension(:), intent(in)    :: tod

      real(dp), allocatable, dimension(:) :: tod_dp, segment_levels, segment_values
      logical,  allocatable, dimension(:) :: good
      integer(i4b)                        :: jump_idx, segment_start, segment_stop, n, maxseg

      ! TODO: Refit function goes here

    end subroutine sample_jump_amplitudes

    subroutine get_jump_tod(self, tod)
      implicit none

      ! Builds the staircase correction
      class(comm_tod_jump),               intent(in)   :: self
      real(sp),             dimension(:), intent(out)  :: tod

      integer(i4b)                                     :: jump_idx, sample_idx

      ! TODO: staircase stuff goes here

    end subroutine get_jump_tod

    subroutine block_medians(samples, good, block_median)
      implicit none

      ! Computes a median for each block of samples
      real(dp),             dimension(:), intent(in)  :: samples
      logical,              dimension(:), intent(in)  :: good
      real(dp),             dimension(:), intent(out) :: block_median

      ! TODO: fill in.
    end subroutine block_medians

    subroutine step_filter(block_median, window_blocks, filter_score, step_difference)
      implicit none

      ! Assigns a score to the jump-candidates
      real(dp),             dimension(:), intent(in)  :: block_median
      integer(i4b),                       intent(in)  :: window_blocks
      real(dp),             dimension(:), intent(out) :: filter_score
      real(dp),             dimension(:), intent(out) :: step_difference

      ! TODO: fill in
    end subroutine step_filter

    integer(i4b) function refine(samples, good, first_sample, last_sample, window_samples) result(best_sample)
      implicit none

      ! Refines (=determines) the accurate jump position
      real(dp),             dimension(:), intent(in)  :: samples
      logical,              dimension(:), intent(in)  :: good
      integer(i4b),                       intent(in)  :: first_sample, last_sample, window_samples

      ! TODO: fill in
    end function refine

    real(dp) function robust_noise_scale(values) result(noise_scale)
      implicit none

      ! Estimates the noise scale via median absolute deviation (mad)
      real(dp),             dimension(:), intent(in)  :: values

      logical,  allocatable, dimension(:) :: finite
      real(dp), allocatable, dimension(:) :: finite_values, deviations
      real(dp) :: center

      ! TODO: Fill out
    end function robust_noise_scale

    real(dp) function median_values(values) result(median)
      implicit none

      real(dp),            dimension(:), intent(in)  :: values

      real(dp), allocatable, dimension(:)  :: sorted_values
      integer(i4b)                         :: n, middle

      ! TODO: write equivalent to np.median here OR resuse the one in comm_utils and delete. 
    end function median_values


end module comm_tod_jump_mod
