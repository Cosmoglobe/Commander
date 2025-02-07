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
module comm_B_FIRAS_mod
  use comm_B_mod
  implicit none

  private
  public comm_B_FIRAS, comm_B_FIRAS_ptr

  type, extends (comm_B) :: comm_B_FIRAS
     integer(i4b) :: nside_hires, npix_hires
     real(dp)     :: M_ecl(3,3)
     real(dp), allocatable, dimension(:,:) :: vecs
   contains
     ! Data procedures
     procedure :: conv           => matmulB_firas
     procedure :: deconv         => matmulInvB_firas
     procedure :: update         => updateBeam
  end type comm_B_FIRAS

  interface comm_B_FIRAS
     procedure constructor
  end interface comm_B_FIRAS

  type comm_B_FIRAS_ptr
     type(comm_B_FIRAS), pointer :: p => null()
  end type comm_B_FIRAS_ptr

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
    class(comm_B_FIRAS),   pointer                              :: constructor

    integer(i4b)       :: i
    character(len=4)   :: nside_text
    character(len=512) :: dir
    real(dp)           :: phi_ecl, theta_ecl
    
    ! General parameters
    allocate(constructor)
       
  end function constructor
  
  subroutine matmulB_firas(self, trans, map)
    implicit none
    class(comm_B_FIRAS), intent(in)    :: self
    logical(lgt),        intent(in)    :: trans
    class(comm_map),     intent(inout) :: map

    real(dp)           :: vec(3), theta, vec0(3), weight, w, vec_rot(3), vec1(3)
    real(dp)           :: M_rot(3,3), dalpha0, alpha0, dalpha1, alpha1
    real(dp)           :: phi0, phi1, phi_rot, theta_rot, M_sat(3,3), alpha, phi_ecl, theta_ecl

  end subroutine matmulB_firas

  subroutine matmulInvB_firas(self, trans, map)
    implicit none
    class(comm_B_FIRAS), intent(in)    :: self
    logical(lgt),     intent(in)    :: trans
    class(comm_map),  intent(inout) :: map

    integer(i4b) :: i, l, j

    write(*,*) 'matmulInvB not implemented in b_FIRAS'
    call mpi_finalize(i)
    stop

  end subroutine matmulInvB_firas

  subroutine updateBeam(self, b_l_norm, mb_eff) 
    implicit none
    class(comm_B_FIRAS),                intent(inout)           :: self
    real(dp),         dimension(0:,1:), intent(in),    optional :: b_l_norm
    real(dp),                           intent(in),    optional :: mb_eff

    if (present(mb_eff)) self%mb_eff = mb_eff
    
  end subroutine updateBeam

  subroutine crossproduct(vector1, vector2, crossvector)
    implicit none 
    real(dp), dimension(3), intent(in)  :: vector1, vector2
    real(dp), dimension(3), intent(out) :: crossvector
    crossvector=[vector1(2)*vector2(3)-vector1(3)*vector2(2), &
               & vector1(3)*vector2(1)-vector1(1)*vector2(3), &
               & vector1(1)*vector2(2)-vector1(2)*vector2(1) ]
  end subroutine crossproduct
  
end module comm_B_FIRAS_mod
