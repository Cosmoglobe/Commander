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
module comm_B_bl_mod
  use comm_B_mod
  implicit none

  private
  public comm_B_bl, comm_B_bl_ptr

  type, extends (comm_B) :: comm_B_bl
   contains
     ! Data procedures
     procedure :: conv           => matmulB_bl
     procedure :: deconv         => matmulInvB_bl
     procedure :: update         => updateBeam
  end type comm_B_bl

  interface comm_B_bl
     procedure constructor
  end interface comm_B_bl

  type comm_B_bl_ptr
     type(comm_B_bl), pointer :: p => null()
  end type comm_B_bl_ptr

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, id_abs, fwhm, mb_eff, nside, pixwin, init_realspace)
    implicit none
    type(comm_params),                  intent(in)           :: cpar
    type(comm_mapinfo), target,         intent(in)           :: info
    integer(i4b),                       intent(in)           :: id, id_abs
    real(dp),                           intent(in), optional :: fwhm
    real(dp),                           intent(in), optional :: mb_eff
    character(len=*),                   intent(in), optional :: pixwin
    integer(i4b),                       intent(in), optional :: nside
    logical(lgt),                       intent(in), optional :: init_realspace
    class(comm_B_bl),   pointer                              :: constructor

    character(len=4)   :: nside_text
    character(len=512) :: dir
    
    ! General parameters
    allocate(constructor)
    
  end function constructor
  
  subroutine matmulB_bl(self, trans, map)
    implicit none
    class(comm_B_bl), intent(in)    :: self
    logical(lgt),     intent(in)    :: trans
    class(comm_map),  intent(inout) :: map


  end subroutine matmulB_bl

  subroutine matmulInvB_bl(self, trans, map)
    implicit none
    class(comm_B_bl), intent(in)    :: self
    logical(lgt),     intent(in)    :: trans
    class(comm_map),  intent(inout) :: map


  end subroutine matmulInvB_bl

  subroutine updateBeam(self, b_l_norm, mb_eff) 
    implicit none
    class(comm_B_bl),                   intent(inout)           :: self
    real(dp),         dimension(0:,1:), intent(in),    optional :: b_l_norm
    real(dp),                           intent(in),    optional :: mb_eff

    if (present(mb_eff)) self%mb_eff = mb_eff
    
  end subroutine updateBeam
  
  
end module comm_B_bl_mod
