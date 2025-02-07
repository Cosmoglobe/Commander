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
module comm_exp_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_exp_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_exp_comp
   contains
     procedure :: S    => evalSED_exp
  end type comm_exp_comp

  interface comm_exp_comp
     procedure constructor_exp
  end interface comm_exp_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_exp(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_exp_comp), pointer   :: c

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


  end function constructor_exp

  ! Definition:
  !    SED  = exp(theta * (nu-nu_ref) / nu_ref)
  ! where 
  !    beta = theta(1)
  function evalSED_exp(self, nu, band, pol, theta)
    implicit none
    class(comm_exp_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_exp

    if (nu < 10d9) then
       evalSED_exp = 0.d0
    else
       evalSED_exp = exp(theta(1) * (nu-self%nu_ref(pol))/self%nu_ref(pol))
    end if

  end function evalSED_exp
  
end module comm_exp_comp_mod
