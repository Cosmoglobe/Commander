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
module comm_curvature_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_curvature_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_curvature_comp
   contains
     procedure :: S    => evalSED_curvature
  end type comm_curvature_comp

  interface comm_curvature_comp
     procedure constructor_curvature
  end interface comm_curvature_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_curvature(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_curvature_comp), pointer   :: c

    integer(i4b) :: i, j, k, l, m, n, p, ierr
    type(comm_mapinfo), pointer :: info => null()
    real(dp)           :: par_dp
    integer(i4b), allocatable, dimension(:) :: sum_pix
    real(dp),    allocatable, dimension(:) :: sum_theta, sum_proplen, sum_nprop
    character(len=512) :: temptxt, partxt
    integer(i4b) :: smooth_scale, p_min, p_max
    class(comm_mapinfo), pointer :: info2 => null()
    class(comm_map),     pointer :: tp => null() 
    class(comm_map),     pointer :: tp_smooth => null() 

    ! General parameters
    allocate(c)


  end function constructor_curvature

  ! Definition:
  !    SED  = (nu/nu_ref)**(beta+0.5*C_s*ln(nu/nu_ref))
  ! where 
  !    beta = theta(1), C_s = theta(2)
  function evalSED_curvature(self, nu, band, pol, theta)
    implicit none
    class(comm_curvature_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_curvature

    evalSED_curvature = (nu/self%nu_ref(pol))**(theta(1)+theta(2)*log(nu/self%nu_ref(pol)))

  end function evalSED_curvature
  
end module comm_curvature_comp_mod
