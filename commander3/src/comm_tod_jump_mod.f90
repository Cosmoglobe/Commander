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

  real(dp),     parameter :: min_fraction = 0.4_dp
  real(dp),     parameter :: detection_threshold = 5.0_dp
  integer(i4b), parameter :: max_jumps = 10_i4b ! Safe-guard, more than 10 jumps can indicate a drift instead

  type :: comm_tod_jump
    integer(i4b) :: njump
    integer(dp), allocatable, dimension(:)     :: jump_pos
    integer(dp), allocatable, dimension(:,:)   :: mask       ! (njump,2), start/stop
    real(sp),    allocatable, dimension(:)     :: jump_amp
    real(dp),    allocatable, dimension(:)     :: scales_s
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
      c%njump = 0_i4b
      c%maxseg = 0_i4b ! this will get set in-code by the samprate
      allocate(c%jump_pos(0), c%jump_amp(0), c%mask(0,2))
      c%block_size = 64_i4b            !TODO: I've asked Simon if this should be changeable
      c%scales_s = [240.0_dp, 60.0_dp] !TODO: Same for this!

    end function constructor_jump

    ! Detects jumps in the residual and stores their positions and amplitudes
    subroutine find_jumps(self, tod, signal, flag, samprate)
      implicit none
      
      class(comm_tod_jump),                intent(inout) :: self
      real(sp),              dimension(:), intent(inout) :: tod
      real(sp),              dimension(:), intent(in)    :: signal
      integer(i4b),          dimension(:), intent(inout) :: flag
      real(dp),                            intent(in)    :: samprate

      ! TODO: check that these are the names we want
      real(sp), allocatable, dimension(:) :: residual, residual_work
      real(dp), allocatable, dimension(:) :: block_median, filter_score, step_diff
      logical,  allocatable, dimension(:) :: good
      integer(i4b)                        :: n, nblock, iteration, scale_idx, block_idx
      integer(i4b)                        :: sample_idx, best_block_idx, window_blocks
      integer(i4b)                        :: hi, lo, jump_idx, other_idx
      integer(dp)                         :: pos_temp
      real(dp)                            :: best_score, score, noise_scale, step_amp

      ! Clear the previous jump solution
      self%njump = 0
      deallocate(self%jump_pos, self%jump_amp, self%mask)
      allocate(self%jump_pos(0), self%jump_amp(0), self%mask(0,2))

      self%maxseg = int(maxval(self%scales_s) * samprate, i4b)
      
      n = size(tod)
      allocate(residual(n), residual_work(n), good(n))
      residual = tod - signal
      good = ieee_is_finite(residual) ! ignores non-finite residual values
      if (.not. any(good)) return

      ! Centers the residual and creates a work-copy 
      residual_work  = residual - median(pack(residual, good))
      
      !! Perform matched filter
      if (self%block_size <= 0_i4b) return
      nblock = n / self%block_size
      if (nblock == 0) return
      allocate(block_median(nblock), filter_score(nblock), step_diff(nblock))
      
      do iteration = 1, max_jumps
        ! recomputes block-medians after subtracting any previously found jumps
        call block_medians(residual_work, good, self%block_size, block_median)

        best_score = 0.0_dp
        best_block_idx = 0
        sample_idx = 0
        step_amp = 0.0_dp

        do scale_idx = 1, size(self%scales_s)
          window_blocks = max(4_i4b, int(self%scales_s(scale_idx) * samprate, i4b) / self%block_size)
          call step_filter(block_median, window_blocks, filter_score, step_diff)

          ! scores the candidates
          noise_scale = robust_noise_scale(filter_score)
          if (.not. ieee_is_finite(noise_scale)) cycle
          if (noise_scale <= 0.0_dp) cycle
          
          do block_idx = 1, nblock
            if (.not. ieee_is_finite(filter_score(block_idx))) cycle

            score = abs(filter_score(block_idx)) / noise_scale
            if (score > best_score) then
              best_score = score
              best_block_idx = block_idx
              step_amp = step_diff(block_idx)
            end if
          end do
        end do
        
        ! stop searching if none of the candidates reaches the threshold
        if (best_score < detection_threshold) exit
        
        ! each iter finds and refines *one* jump then updates residual_work
        lo = (best_block_idx - 2) * self%block_size + 1 ! TODO: double check that I got hte indexing right
        hi = best_block_idx * self%block_size + 1
        sample_idx = refine(residual_work, good, lo, hi, 2_i4b * self%block_size)
        
        ! Jump indicies are stored in self%jump_pos, amplitudes in self%jump_amp
        self%njump = self%njump + 1_i4b
        self%jump_pos = [self%jump_pos, int(sample_idx, dp)]
        self%jump_amp = [self%jump_amp, real(step_amp, sp)]

        residual_work(sample_idx:) = residual_work(sample_idx:) - real(step_amp, sp)
      
      end do
      
      ! Sorts the detected pos. and refits amplitudes using the full residual 
      do jump_idx = 1, self%njump - 1
        do other_idx = jump_idx + 1, self%njump
          if (self%jump_pos(other_idx) < self%jump_pos(jump_idx)) then
            pos_temp = self%jump_pos(jump_idx)
            self%jump_pos(jump_idx) = self%jump_pos(other_idx)
            self%jump_pos(other_idx) = pos_temp
          end if
        end do
      end do

      call self%sample_jump_amplitudes(residual)


      ! Identify bad range; store in self%mask
      ! TODO: I need to check with HK or someone else about this. 

      ! Update flag array
      ! TODO: Same as above.
      
    end subroutine find_jumps

    ! Refits the jump ampl. from the residual levels around each jump
    subroutine sample_jump_amplitudes(self, residual)
      implicit none

      class(comm_tod_jump),               intent(inout) :: self
      real(sp),             dimension(:), intent(in)    :: residual

      real(dp),     allocatable, dimension(:) :: seg_levels, seg_values
      logical,      allocatable, dimension(:) :: good
      integer(i4b), allocatable, dimension(:) :: edges
      integer(i4b)                            :: seg_idx, jump_idx, seg_start, seg_stop, n, maxseg

      n = size(residual)

      ! Resize the jump_amp array to match the current number of jumps
      if (size(self%jump_amp) /= self%njump) then
        deallocate(self%jump_amp)
        allocate(self%jump_amp(self%njump))
      end if

      self%jump_amp = 0.0_sp
      if (self%njump == 0) return
      
      maxseg = self%maxseg
      allocate(good(n), seg_levels(self%njump + 1), edges(self%njump + 2))
      good = ieee_is_finite(residual)
      seg_levels = ieee_value(0.0_dp, ieee_quiet_nan)

      edges(1) = 1
      edges(2:self%njump + 1) = int(self%jump_pos(1:self%njump), i4b)
      edges(self%njump + 2) = n + 1

      do seg_idx = 1, self%njump + 1
        seg_start = edges(seg_idx)
        seg_stop = edges(seg_idx + 1)

        ! First segment measured at its end, the rest are measured at their start
        if (seg_idx == 1) then
          seg_start = max(seg_start, seg_stop - maxseg)
        else
          seg_stop = min(seg_stop, seg_start + maxseg)
        end if
        
        seg_values = pack(residual(seg_start:seg_stop - 1), &
                          good(seg_start:seg_stop - 1))

        ! At least two samples are needed to estimate the level
        if (size(seg_values) >= 2) then
          seg_levels(seg_idx) = median(seg_values)
        end if
      end do

      do jump_idx = 1, self%njump
        ! The ampl. is zero if either of the neighbouring levels are invalid
        if (ieee_is_finite(seg_levels(jump_idx)) .and. ieee_is_finite(seg_levels(jump_idx + 1))) then
          self%jump_amp(jump_idx) = real(seg_levels(jump_idx + 1) - seg_levels(jump_idx), sp)
        end if
      end do


    ! Builds the staircase correction (jump_tod)
    subroutine get_jump_tod(self, jump_tod)
      implicit none

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

    ! Computes a median for each block with enough good samples
    subroutine block_medians(samples, good, block_size, block_median)
      implicit none

      real(sp),             dimension(:), intent(in)  :: samples
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

    ! Assigns a score to the jump-candidates via a matched filter with a step func. 
    subroutine step_filter(block_median, window_blocks, filter_score, step_diff)
      implicit none

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

    ! Refines (=determines) the accurate jump position
    integer(i4b) function refine(samples, good, lo, hi, window_samples) result(best_sample)
      implicit none

      real(sp),             dimension(:), intent(in)  :: samples
      logical,              dimension(:), intent(in)  :: good
      integer(i4b),                       intent(in)  :: lo, hi, window_samples
      
      real(dp), allocatable, dimension(:) :: before, after
      integer(i4b)                        :: candidate_sample, before_window_start, after_window_end
      real(dp)                            :: step_diff, best_diff
      
      best_sample = max(lo, 2_i4b)
      best_diff = -1.0_dp
      if (window_samples <= 0) return
     
      ! Searches candidate positions in [lo, hi) to match the Python behaviour  
      ! Starts at sample 2 to make sure at least one sample is in 'before'
      do candidate_sample = max(lo, 2_i4b), min(hi - 1, size(samples) - 1)  
        before_window_start = max(1_i4b, candidate_sample - window_samples)
        after_window_end = min(size(samples), candidate_sample + window_samples - 1)

        before = pack(samples(before_window_start:candidate_sample - 1), &
                      good(before_window_start:candidate_sample - 1))
        after = pack(samples(candidate_sample:after_window_end), &
                     good(candidate_sample:after_window_end))
        
        if (size(before) < 2 .or. size(after) < 2) cycle

        step_diff = abs(median(after) - median(before))
        if (step_diff > best_diff) then
          best_sample = candidate_sample
          best_diff = step_diff
        end if
      end do
    
    end function refine

    ! Estimates the noise scale using the scaled median absolute deviation (mad)
    real(dp) function robust_noise_scale(values) result(noise_scale)
      implicit none

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
