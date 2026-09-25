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

    subroutine find_jumps(self, tod, signal, flag)
      implicit none
      class(comm_tod_jump),               intent(inout) :: self
      real(sp),             dimension(:), intent(inout) :: tod
      real(sp),             dimension(:), intent(in)    :: signal
      integer(i4b),         dimension(:), intent(inout) :: flag

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
      class(comm_tod_jump),               intent(inout) :: self
      real(sp),             dimension(:), intent(in)    :: tod

    end subroutine sample_jump_amplitudes

    subroutine get_jump_tod(self, tod)
      implicit none
      class(comm_tod_jump),               intent(in)     :: self
      real(sp),             dimension(:), intent(out)    :: tod

    end subroutine get_jump_tod

end module comm_tod_jump_mod
