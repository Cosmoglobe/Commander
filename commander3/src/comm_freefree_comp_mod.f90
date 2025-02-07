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
module comm_freefree_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_freefree_comp

  !**************************************************
  !      Free-free component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_freefree_comp
   contains
     procedure :: S    => evalSED_freefree
  end type comm_freefree_comp

  interface comm_freefree_comp
     procedure constructor_freefree
  end interface comm_freefree_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_freefree(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_freefree_comp), pointer   :: c

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


  end function constructor_freefree

  ! Definition:
  !      x  = h*nu/(k_b*T)
  !    SED  = (nu/nu_ref)**(beta+1) * (exp(x_ref)-1)/(exp(x)-1)
  ! where 
  !    beta = theta(1)
  function evalSED_freefree(self, nu, band, pol, theta)
    implicit none
    class(comm_freefree_comp),    intent(in)      :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_freefree
    real(dp)     :: S, S_ref, EM, T_e
    real(dp)     :: g, g_ref, Z_i, tau, tau_ref, EM1, Te

!!$    EM      = theta(1)
!!$    !EM1 = 1.d0 
!!$    Te      = theta(2)
!!$    Z_i     = 1.d0
!!$    g       = log(exp(5.960d0 - sqrt(3.d0)/pi * log(Z_i * nu/1.d9          * (Te/1.d4)**(-1.5d0))) + 2.71828d0)
!!$    !g_ref   = log(exp(5.960d0 - sqrt(3.d0)/pi * log(Z_i * self%nu_ref/1.d9 * (Te/1.d4)**(-1.5d0))) + 2.71828d0)
!!$    tau     = 5.468d-2 * Te**(-1.5d0) * (nu/1.d9)**(-2)          * EM * g
!!$    !tau_ref = 5.468d-2 * Te**(-1.5d0) * (self%nu_ref/1.d9)**(-2) * EM * g_ref
!!$
!!$    evalSED_freefree = 1.d6 * Te * (1.d0 - exp(-tau)) !/ (1.d0 - exp(-tau_ref)) 
!!$    
!!$    return
!!$    !write(*,*) "1:", evalSED_freefree


    !EM    = theta(1) ! Not used
    T_e   = theta(1)
    S     = log(exp(5.960d0 - sqrt(3.d0)/pi * log(1.d0 * nu    /1.d9 * (T_e/1.d4)**(-1.5d0))) + 2.71828d0)
    S_ref = log(exp(5.960d0 - sqrt(3.d0)/pi * log(1.d0 * self%nu_ref(pol)/1.d9 * (T_e/1.d4)**(-1.5d0))) + 2.71828d0)
    !evalSED_freefree = S/S_ref * exp(-h*(nu-self%nu_ref(pol))/k_b/T_e) * (nu/self%nu_ref(pol))**(-2)
    !evalSED_freefree = S/S_ref * (nu/self%nu_ref(pol))**-2
    evalSED_freefree = S/S_ref * (nu/self%nu_ref(pol))**(-2)
    !write(*,*) "2",evalSED_freefree
    
  end function evalSED_freefree
  
end module comm_freefree_comp_mod
