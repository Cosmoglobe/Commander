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
module comm_powlaw_break_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_powlaw_break_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_powlaw_break_comp
     real(dp), allocatable, dimension(:) :: nu_break
   contains
     procedure :: S    => evalSED_powlaw_break
  end type comm_powlaw_break_comp

  interface comm_powlaw_break_comp
     procedure constructor_powlaw_break
  end interface comm_powlaw_break_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_powlaw_break(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),           intent(in) :: cpar
    integer(i4b),                intent(in) :: id, id_abs
    class(comm_powlaw_break_comp), pointer   :: c

    integer(i4b) :: i, j, k, l, m, n, p, ierr
    integer(i4b), allocatable, dimension(:) :: sum_pix
    real(dp),     allocatable, dimension(:) :: sum_theta, sum_proplen, sum_nprop

    real(dp)                     :: par_dp
    character(len=512)           :: temptxt, partxt
    integer(i4b)                 :: smooth_scale, p_min, p_max
    type(comm_mapinfo),  pointer :: info => null()
    class(comm_mapinfo), pointer :: info2 => null()
    class(comm_map),     pointer :: tp => null() 
    class(comm_map),     pointer :: tp_smooth => null() 

    ! General parameters
    allocate(c)


  end function constructor_powlaw_break

  ! Definition:
  !    SED  = (nu/nu_ref)**(beta+dbeta) if nu < nu_break
  !    SED  = (nu/nu_ref)**(beta) if nu => nu_break
  ! where 
  !    beta = theta(1), dbeta = theta(2)
  function evalSED_powlaw_break(self, nu, band, pol, theta)
    implicit none
    class(comm_powlaw_break_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: beta,dbeta,nu_break
    real(dp)                                      :: evalSED_powlaw_break

    beta     = theta(1)
    dbeta    = theta(2)
    nu_break = self%nu_break(pol)

    if (nu < nu_break) then
       evalSED_powlaw_break = (nu/self%nu_ref(pol))**(beta)
    else
       evalSED_powlaw_break = (nu/self%nu_ref(pol))**(beta+dbeta)
    end if

  end function evalSED_powlaw_break
  
end module comm_powlaw_break_comp_mod
