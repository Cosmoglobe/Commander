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

module comm_tod_spike_mod
  use comm_utils
  use comm_status_mod
  implicit none
  private

  public comm_tod_spike, tod_spike_ptr

  type :: comm_tod_spike
   contains
     procedure :: estimate => estimate_spike_template
     procedure :: generate => generate_spike_correction
  end type comm_tod_spike

  interface comm_tod_spike
     procedure spike_constructor
  end interface comm_tod_spike

  type tod_spike_ptr
     class(comm_tod_spike), pointer :: p => null()
  end type tod_spike_ptr

contains
 
  !initializes a comm_spike class
  function spike_constructor() result(c)
    implicit none
    class(comm_tod_spike), pointer                  :: c

    allocate(c)
    
  end function spike_constructor
 
  ! estimates spike parameters for each detector
  subroutine estimate_spike_template(self)
    implicit none
    class(comm_tod_spike), intent(inout)      :: self
  end subroutine estimate_spike_template

  ! Routine for convolving with bolometer transfer function, tod = T * tod
  subroutine generate_spike_correction(self, s_spike)
    implicit none
    class(comm_tod_spike),              intent(in)  :: self
    real(sp),             dimension(:), intent(out) :: s_spike

    s_spike = 0.
  end subroutine generate_spike_correction

end module comm_tod_spike_mod
