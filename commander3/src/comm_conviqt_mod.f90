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
module comm_conviqt_mod
  use iso_c_binding, only : c_ptr, c_double, c_int
  use comm_map_mod
  use comm_shared_arr_mod
  implicit none

  private
  public comm_conviqt, conviqt_ptr

  type :: comm_conviqt
    integer(i4b) :: lmax, mmax, nmaps, bmax, nside, npix, comm, optim, psisteps, win, myid
    real(dp), allocatable, dimension(:)        :: lnorm
    type(shared_2d_sp) :: c
    !type(shared_1d_int) :: pixLookup
    real(sp) :: psires
    class(comm_mapinfo), pointer :: info => null()
    type(shared_2d_spc) :: alm_beam
  contains
    procedure     :: interp
    procedure     :: precompute_sky
    procedure     :: get_alms
    procedure     :: dealloc
  end type comm_conviqt

  interface comm_conviqt
    procedure constructor
  end interface comm_conviqt

  type conviqt_ptr
     class(comm_conviqt), pointer :: p => null()
  end type conviqt_ptr

contains

  !Constructor

  function constructor(myid_shared, comm_shared, myid_inter, comm_inter, nside, lmax, nmaps, bmax, beam, map, optim)
    implicit none
    integer(i4b),                 intent(in) :: nside, lmax, bmax, nmaps
    integer(i4b),                 intent(in) :: myid_shared, comm_shared, myid_inter, comm_inter
    class(comm_map),              intent(in) :: beam, map
    integer(i4b),                 intent(in) :: optim ! desired optimization flags
    class(comm_conviqt), pointer             :: constructor

    integer(i4b) :: i, j, k, l, m, ierr, nalm_tot, nalm
    integer(i4b), allocatable, dimension(:)   :: ind
    complex(spc), allocatable, dimension(:,:) :: alm
    real(dp) :: theta, phi


    allocate(constructor)

  end function constructor


  !interpolates the precomputed map to the psi angle
  function interp(self, pix, psi)
    implicit none
    class(comm_conviqt),       intent(in) :: self
    integer(i4b),              intent(in) :: pix
    real(dp),                  intent(in) :: psi
    real(sp)                              :: interp

    real(sp)       :: x0, x1, unwrap
    integer(i4b)   :: pixnum, psii, psiu, bpsi

    interp = 0


  end function interp

  !precomputes the maps and stores to the allocated array
  !TODO: handle polarized sidelobes?
  subroutine precompute_sky(self, map)
    implicit none
    class(comm_conviqt),                          intent(inout) :: self
    class(comm_map),                              intent(in)    :: map ! Must contain alms

 
  end subroutine precompute_sky

  ! Note: lnorm is set to 1 here. Todo: Figure out why 
  subroutine get_alms(self, m_b, map, alm)
    implicit none
    class(comm_conviqt),        intent(in)  :: self
    integer(i4b),               intent(in)  :: m_b
    class(comm_map),            intent(in)  :: map
    real(dp), dimension(0:,1:), intent(out) :: alm

    integer(i4b) :: i, ineg, j, l, m
    real(dp)     :: spinsign, mfac, sqrt_two
    complex(dpc) :: v1, v2, alm_b(3), alm_s(3), almc

    spinsign = 1; if (m_b /= 0) spinsign = -1
    sqrt_two = sqrt(2.d0)

    alm = 0.d0

  end subroutine get_alms


  !deallocate memory used by the class
  subroutine dealloc(self)
    implicit none
    class(comm_conviqt), intent(inout) :: self

    deallocate(self%lnorm)
    call dealloc_shared_2d_spc(self%alm_beam)
    call dealloc_shared_2d_sp(self%c)
    !call dealloc_shared_1d_int(self%pixLookup)

  end subroutine dealloc

end module comm_conviqt_mod
