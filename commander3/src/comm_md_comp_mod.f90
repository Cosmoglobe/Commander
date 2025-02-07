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
module comm_md_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_md_comp, initialize_md_comps

  !**************************************************
  !           Monopole/dipole component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_md_comp
     !integer(i4b)                            :: ref_band
     logical(lgt) :: mono_from_prior  ! true if the band used as zero-level prior
     real(dp)     :: mono_alm         ! alm value of monopole for when we CG-sample,
                                      ! can revert to pre CG sampling value afterwards
   contains
     procedure :: S    => evalSED_md
  end type comm_md_comp

  interface comm_md_comp
     procedure constructor_md
  end interface comm_md_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_md(cpar, id, id_abs, band, label, mu, rms, def) result(c)
    implicit none
    type(comm_params),              intent(in) :: cpar
    integer(i4b),                   intent(in) :: id, id_abs, band
    character(len=*),               intent(in) :: label
    real(dp),         dimension(4), intent(in) :: mu, def
    real(dp),         dimension(2), intent(in) :: rms
    class(comm_md_comp), pointer               :: c

    integer(i4b) :: i, j, k, l, m, n
    character(len=16), dimension(1000) :: comp_label
    type(comm_mapinfo), pointer :: info => null()

    ! General parameters
    allocate(c)

    
  end function constructor_md

  ! Definition:
  !    SED  = delta_{band,ref_band}
  function evalSED_md(self, nu, band, pol, theta)
    class(comm_md_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_md

    integer(i4b) :: i, ind

    evalSED_md = 1.d0

  end function evalSED_md

  function initialize_md_comps(cpar, id, id_abs, n)
    implicit none
    type(comm_params),   intent(in)  :: cpar
    integer(i4b),        intent(in)  :: id, id_abs
    integer(i4b),        intent(out) :: n
    class(comm_md_comp), pointer   :: initialize_md_comps

    integer(i4b)        :: i, unit
    real(dp)            :: mu(4), rms(2), def(4)
    character(len=1024) :: line, label
    class(comm_comp), pointer :: c => null()

    initialize_md_comps => null()
    n = 0

  end function initialize_md_comps
  
end module comm_md_comp_mod
