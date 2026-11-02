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
  use comm_utils, only: median
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none

  private
  public comm_tod_jump

  real(dp), parameter :: min_fraction = 0.4_dp

  type :: comm_tod_jump
    integer(i4b) :: njump
    integer(dp), allocatable, dimension(:)     :: jump_pos
    integer(dp), allocatable, dimension(:,:)   :: mask       ! (njump,2), start/stop
    real(sp),    allocatable, dimension(:)     :: jump_amp
    integer(i4b) :: block_size
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
      c%block_size = 64_i4b !TODO: I've asked Simon if this should be changeable

    end function constructor_jump

    subroutine find_jumps(self, tod, signal, flag, samprate)
      implicit none
      
      ! Detects jumps in the residual and stores their positions and amplitudes
      class(comm_tod_jump),                intent(inout) :: self
      real(sp),              dimension(:), intent(inout) :: tod
      real(sp),              dimension(:), intent(in)    :: signal
      integer(i4b),          dimension(:), intent(inout) :: flag
      real(dp),                            intent(in)    :: samprate !KG NOTE: make sure to include this in the call for find_jumps.

      ! TODO: check that these are the names we want
      real(sp), allocatable, dimension(:) :: residual
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

    subroutine get_jump_tod(self, jump_tod)
      implicit none

      ! Builds the staircase correction (jump_tod)
      class(comm_tod_jump),               intent(in)   :: self
      real(sp),             dimension(:), intent(out)  :: jump_tod !TODO "stair" in Simons code, ok to change?

      integer(i4b)                                     :: jump_idx
      integer(dp)                                      :: sample_idx

      jump_tod = 0.0_sp

      if (.not. allocated(self%jump_pos)) return
      if (.not. allocated(self%jump_amp)) return

      do jump_idx = 1, min(self%njump, size(self%jump_pos), size(self%jump_amp))
        sample_idx = self%jump_pos(jump_idx)
        if (sample_idx < 1 .or. sample_idx > size(jump_tod)) cycle

        jump_tod(sample_idx:) = jump_tod(sample_idx:) + self%jump_amp(jump_idx)
      end do

    end subroutine get_jump_tod

    subroutine block_medians(samples, good, block_size, block_median)
      implicit none

      ! Computes a median for each block with enough good samples
      real(dp),             dimension(:), intent(in)  :: samples
      logical,              dimension(:), intent(in)  :: good
      integer(i4b),                       intent(in)  :: block_size
      real(dp),             dimension(:), intent(out) :: block_median
      
      real(dp), allocatable, dimension(:) :: block_values
      integer(i4b)                        :: block_idx, first_sample, last_sample

      block_median = ieee_value(0.0_dp, ieee_quiet_nan)
      if (block_size <= 0) return

      do block_idx = 1, size(block_median)
        first_sample = (block_idx - 1) * block_size + 1
        last_sample = block_idx * block_size
        block_values = pack(samples(first_sample:last_sample), good(first_sample:last_sample))

        if (size(block_values) > 0 .and. size(block_values) >= block_size / 2) then
          block_median(block_idx) = median(block_values)
        end if
      end do

    end subroutine block_medians

    subroutine step_filter(block_median, window_blocks, filter_score, step_diff)
      implicit none

      ! Assigns a score to the jump-candidates via a matched filter with a step func. 
      real(dp),             dimension(:), intent(in)  :: block_median
      integer(i4b),                       intent(in)  :: window_blocks
      real(dp),             dimension(:), intent(out) :: filter_score, step_diff
      
      real(dp),     allocatable, dimension(:) :: cumulative_sum
      integer(i4b), allocatable, dimension(:) :: cumulative_count
      integer(i4b)                            :: nblock, block_idx, lo, hi, before_count, after_count
      real(dp)                                :: before_mean, after_mean

      nblock = size(block_median)
      filter_score = ieee_value(0.0_dp, ieee_quiet_nan)
      step_diff = ieee_value(0.0_dp, ieee_quiet_nan)
      if (window_blocks <= 0) return

      allocate(cumulative_sum(nblock + 1), cumulative_count(nblock + 1))
      cumulative_sum(1) = 0.0_dp
      cumulative_count(1) = 0

      do block_idx = 1, nblock
        cumulative_sum(block_idx + 1) = cumulative_sum(block_idx)
        cumulative_count(block_idx + 1) = cumulative_count(block_idx)
        if (ieee_is_finite(block_median(block_idx))) then
          cumulative_sum(block_idx + 1) = cumulative_sum(block_idx + 1) + block_median(block_idx)
          cumulative_count(block_idx + 1) = cumulative_count(block_idx + 1) + 1
        end if
      end do
      
      do block_idx = 1, nblock
        lo = max(1_i4b, block_idx - window_blocks)
        hi = min(nblock, block_idx + window_blocks - 1)
        before_count = cumulative_count(block_idx) - cumulative_count(lo)
        after_count = cumulative_count(hi + 1) - cumulative_count(block_idx)
        if (real(before_count, dp) < min_fraction * window_blocks .or. &
            real(after_count, dp) < min_fraction * window_blocks) cycle
        
        before_mean = (cumulative_sum(block_idx) - cumulative_sum(lo)) / real(before_count, dp)
        after_mean = (cumulative_sum(hi + 1) - cumulative_sum(block_idx)) / real(after_count, dp)
        step_diff(block_idx) = after_mean - before_mean

        filter_score(block_idx) = step_diff(block_idx)/sqrt(1.0_dp/after_count + 1.0_dp/before_count)
      end do

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

      ! Estimates the noise scale using the scaled median absolute deviation (mad)
      real(dp),              dimension(:), intent(in)  :: values

      logical,  allocatable, dimension(:) :: finite
      real(dp), allocatable, dimension(:) :: finite_values

      noise_scale = 0.0_dp
      finite = ieee_is_finite(values)
      finite_values = pack(values, finite)
      if (size(finite_values) == 0) return

      noise_scale = 1.4826_dp * median(abs(finite_values - median(finite_values)))

    end function robust_noise_scale


end module comm_tod_jump_mod
