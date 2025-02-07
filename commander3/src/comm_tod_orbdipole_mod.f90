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
module comm_tod_orbdipole_mod
  use comm_map_mod
  implicit none

  type :: comm_orbdipole
    integer(i4b) :: ndet!, subsample
    integer(i4b) :: comm
    logical(lgt) :: beam_4pi
    real(dp),       dimension(:,:), allocatable :: orb_dp_s !precomputed s integrals for orbital dipole sidelobe term
    class(map_ptr), dimension(:),   allocatable :: beam
  contains
    procedure :: precompute_orb_dp_s
    procedure :: compute_CMB_dipole
    procedure :: compute_4pi_product
  end type comm_orbdipole

  interface comm_orbdipole
    procedure constructor
  end interface comm_orbdipole

  type orbdipole_pointer
    class(comm_orbdipole), pointer :: p => null()
  end type orbdipole_pointer

contains

  function constructor(beam, comm)
    implicit none
    class(map_ptr),      dimension(:), target, optional :: beam
    integer(i4b),                              optional :: comm
    class(comm_orbdipole), pointer :: constructor
    integer(i4b) :: i

    allocate(constructor)

  end function constructor

  subroutine precompute_orb_dp_s(self) 
    implicit none
    class(comm_orbdipole),              intent(inout) :: self

    integer(i4b) :: i, j, ierr
    real(dp), dimension(3) :: v
    real(dp) :: pixVal



  end subroutine precompute_orb_dp_s

  subroutine compute_CMB_dipole(self, det, v_ref, nu, &
       & relativistic, beam_4pi, P, s_dip, factor, v_ref_next)
    ! Evaluates the CMB dipole as a function of time
    !
    !
    !   Arguments:
    !   ---------
    !   self: comm_orbdipole object
    !
    !   det: int
    !        detector index
    !   v_ref: double, array of length 3
    !        velocity of observer in km/s, Galactic coordinates
    !   relativistic: logical
    !        if True, computes relativistic correction
    !   beam_4pi: logical
    !        if True, uses the full main beam map, else uses pencil beam.
    !   P: double, array
    !        Pointing array
    !        if beam_4pi, array of phi/theta/psi values for TOD
    !        else, array of unit vectors for TOD pointing
    !   factor: double (optional)
    !        multiplicative factor if ad-hoc unit correction is needed.
    !
    !   Returns:
    !   --------
    !   s_dip: real (sp)
    !        Array of dipole template timestreams for given detector
    implicit none
    class(comm_orbdipole),                 intent(inout)  :: self
    integer(i4b),                          intent(in)  :: det
    real(dp),                              intent(in)  :: v_ref(3), nu
    logical(lgt),                          intent(in)  :: relativistic 
    logical(lgt),                          intent(in)  :: beam_4pi
    real(dp),            dimension(1:,1:), intent(in)  :: P
    real(sp),            dimension(:),     intent(out) :: s_dip
    real(dp),                              intent(in), optional :: factor
    real(dp),                              intent(in), optional :: v_ref_next(3)

    s_dip = 0

  end subroutine compute_CMB_dipole


  function compute_4pi_product(self, det, q, Pang, v_ref) result(prod)
    implicit none
    class(comm_orbdipole), intent(in) :: self
    integer(i4b),          intent(in) :: det
    real(dp),              intent(in) :: q, Pang(3), v_ref(3)
    real(dp)                          :: prod

    real(dp), dimension(3,3) :: rot_mat
    real(dp), dimension(3)   :: vnorm


    prod = 0d0

  end function compute_4pi_product

end module comm_tod_orbdipole_mod
