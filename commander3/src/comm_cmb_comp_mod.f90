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
module comm_cmb_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_cmb_comp, apply_cmb_dipole_prior

  !**************************************************
  !           CMB component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_cmb_comp
     real(dp) :: cmb_dipole_prior(3)
   contains
     procedure :: S            => evalSED_cmb
     procedure :: update_F_int => updateIntF_cmb
  end type comm_cmb_comp

  interface comm_cmb_comp
     procedure constructor_cmb
  end interface comm_cmb_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_cmb(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_cmb_comp), pointer   :: c

    integer(i4b) :: i, j, k
    real(dp)     :: f
    

    ! General parameters
    allocate(c)

  end function constructor_cmb

  ! Definition:
  !    SED  = conversion between thermodynamic and brightness temperature = 1/a2t
  function evalSED_cmb(self, nu, band, pol, theta)
    implicit none
    class(comm_cmb_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_cmb

    real(dp) :: x
    x           = 0
    evalSED_cmb = (x**2 * exp(x)) / (exp(x)-1.d0)**2

  end function evalSED_cmb

  ! Update band integration lookup tables
  subroutine updateIntF_cmb(self, band)
    implicit none
    class(comm_cmb_comp),                    intent(inout)        :: self
    integer(i4b),                            intent(in), optional :: band

    integer(i4b) :: i, j, k
    real(dp)     :: f

    do k = 1, 3
       do i = 1, numband
          do j = 0, data(i)%ndet
             if (k > 1) then
                if (self%nu_ref(k) == self%nu_ref(k-1)) then
                   self%F_int(k,i,j)%p => self%F_int(k-1,i,j)%p
                   cycle
                end if
             end if
          end do
       end do
    end do

  end subroutine updateIntF_cmb

  subroutine apply_cmb_dipole_prior(cpar, handle)
    implicit none
    type(comm_params),   intent(in)    :: cpar
    type(planck_rng),    intent(inout) :: handle

    integer(i4b) :: i, l, m
    real(dp)     :: md(0:3)
    class(comm_map),  pointer :: map
    class(comm_comp), pointer :: c => null()


  end subroutine apply_cmb_dipole_prior


end module comm_cmb_comp_mod
