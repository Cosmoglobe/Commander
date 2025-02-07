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
module comm_pah_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_pah_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_pah_comp
     real(dp)          :: nu_p0
   contains
     procedure :: S    => evalSED_pah
  end type comm_pah_comp

  interface comm_pah_comp
     procedure constructor_pah
  end interface comm_pah_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_pah(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_pah_comp), pointer   :: c

    integer(i4b) :: ind(1), i, j, k
    real(dp), allocatable, dimension(:,:) :: SED


    ! General parameters
    allocate(c)

  end function constructor_pah

  function evalSED_pah(self, nu, band, pol, theta)
    implicit none
    class(comm_pah_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_pah

    evalSED_pah = 0d0

  end function evalSED_pah
  
end module comm_pah_comp_mod
