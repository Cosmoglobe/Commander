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
module comm_ame_lognormal_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_ame_lognormal_comp

  !**************************************************
  !        AME - lognormal width component
  !  https://arxiv.org/pdf/2001.07159.pdf (Eqn. 6)
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_ame_lognormal_comp
     contains
       procedure :: s => evalSED_amelognormal
  end type comm_ame_lognormal_comp

  interface comm_ame_lognormal_comp
     procedure constructor_ame_lognormal
  end interface comm_ame_lognormal_comp

contains

  function constructor_ame_lognormal(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),  intent(in) :: cpar
    integer(i4b),       intent(in) :: id, id_abs
    class(comm_ame_lognormal_comp), pointer :: c

    integer(i4b) :: i, j, k
    type(comm_mapinfo), pointer :: info => null()

    allocate(c)


  end function constructor_ame_lognormal

  ! Definition:
  !
  ! S_AME(nu) = A_AME * exp{ -1/2 * [ln(nu/nu_AME)/W_AME]**2}
  !
  ! where nu_AME is the peak frequency and W_AME is the width 
  ! of the spectrum

  function evalSED_amelognormal(self, nu, band, pol, theta)
    implicit none
    class(comm_ame_lognormal_comp),    intent(in) :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_amelognormal
    real(dp)                                      :: peak
    real(dp)                                      :: nu_AME, W_AME

    nu_AME  = theta(1)
    W_AME   = theta(2)

    peak    = nu_AME*1d9

    ! Moving from flux density to T_RJ requires an multiplicative (1/nu^2) 
    
    evalSED_amelognormal = exp(-0.5*(log(nu/peak)/W_AME)**2)*(self%nu_ref(pol)/nu)**2
  end function evalSED_amelognormal
end module comm_ame_lognormal_mod
