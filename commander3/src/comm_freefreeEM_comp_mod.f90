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
module comm_freefreeEM_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_freefreeEM_comp

  !*******************************************************
  !      Free-free with non-linear emission measure (EM)
  !      Note: Amplitude map must always be fixed to 1
  !*******************************************************
  type, extends (comm_diffuse_comp) :: comm_freefreeEM_comp
   contains
     procedure :: S    => evalSED_ffEM
  end type comm_freefreeEM_comp

  interface comm_freefreeEM_comp
     procedure constructor_ffEM
  end interface comm_freefreeEM_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_ffEM(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_freefreeEM_comp), pointer   :: c

    integer(i4b) :: i, j, k, l, m, n, p, ierr
    type(comm_mapinfo), pointer :: info => null()
    real(dp)           :: par_dp
    integer(i4b), allocatable, dimension(:) :: sum_pix
    real(dp),    allocatable, dimension(:) :: sum_theta, sum_proplen, sum_nprop
    character(len=512) :: temptxt, partxt
    integer(i4b) :: smooth_scale, p_min, p_max
    class(comm_mapinfo), pointer :: info2 => null()

    ! General parameters
    allocate(c)



  end function constructor_ffEM

  ! Definition:
  function evalSED_ffEM(self, nu, band, pol, theta)
    implicit none
    class(comm_freefreeEM_comp),    intent(in)      :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_ffEM
    real(dp)     :: S, S_ref, EM, T_e
    real(dp)     :: g, g_ref, Z_i, tau, tau_ref, EM1, Te

    EM      = exp(theta(1))
    Te      = theta(2)
    Z_i     = 1.d0
    g       = log(exp(5.960d0 - sqrt(3.d0)/pi * log(Z_i * nu/1.d9 * (Te/1.d4)**(-1.5d0))) + 2.71828d0)
    tau     = 5.468d-2 * Te**(-1.5d0) * (nu/1.d9)**(-2) * EM * g

    evalSED_ffEM = 1.d6 * Te * (1.d0 - exp(-tau)) 

  end function evalSED_ffEM

  
end module comm_freefreeEM_comp_mod
