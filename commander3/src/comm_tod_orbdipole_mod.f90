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
  use comm_utils
  use comm_defs
  use comm_map_mod
  use spline_1D_mod
  implicit none

  private 
  public comm_orbdipole, orbdipole_pointer

  type :: comm_orbdipole
    integer(i4b) :: ndet!, subsample
    integer(i4b) :: comm
    logical(lgt) :: beam_4pi
    real(dp),       dimension(:,:), allocatable :: orb_dp_s !precomputed s integrals for orbital dipole sidelobe term
    class(map_ptr), dimension(:),   allocatable :: beam
    type(spline_type) :: s
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
    !constructor%subsample = 20
    if (present(beam)) then
       constructor%ndet      = size(beam)
       constructor%beam_4pi  = .false.
       constructor%comm      = beam(1)%p%info%comm

       allocate(constructor%beam(constructor%ndet))
       do i=1, constructor%ndet
          constructor%beam(i)%p => beam(i)%p
       end do
       
       call constructor%precompute_orb_dp_s()
    else
       constructor%ndet      = 1
       constructor%beam_4pi  = .false.
       constructor%comm      = comm
    end if

  end function constructor

  subroutine precompute_orb_dp_s(self) 
    implicit none
    class(comm_orbdipole),              intent(inout) :: self

    integer(i4b) :: i, j, ierr
    real(dp), dimension(3) :: v
    real(dp) :: pixVal

    allocate(self%orb_dp_s(self%ndet, 10)) 
    do i = 1, self%ndet
      self%orb_dp_s(i, :) = 0
      do j = 0, self%beam(i)%p%info%np-1
         call pix2vec_ring(self%beam(i)%p%info%nside, self%beam(i)%p%info%pix(j+1),v)
         pixVal = self%beam(i)%p%map(j, 1)
         !if(pixVal < 0) then
         !   pixVal = 0.d0
         !end if
         !x
         self%orb_dp_s(i, 1) = self%orb_dp_s(i, 1) + pixVal * v(1)
         !y 
         self%orb_dp_s(i, 2) = self%orb_dp_s(i, 2) + pixVal * v(2)
         !z 
         self%orb_dp_s(i, 3) = self%orb_dp_s(i, 3) + pixVal * v(3)
         !x^2 
         self%orb_dp_s(i, 4) = self%orb_dp_s(i, 4) + pixVal*v(1)*v(1)
         !2xy 
         self%orb_dp_s(i, 5) = self%orb_dp_s(i, 5) + 2.d0 * pixVal* v(1) * v(2)
         !2xz 
         self%orb_dp_s(i, 5) = self%orb_dp_s(i, 5) + 2.d0 * pixVal * v(1) *v(2)
         !2xz 
         self%orb_dp_s(i, 6) = self%orb_dp_s(i, 6) + 2.d0 * pixVal* v(1) * v(3)
         !y^2 
         self%orb_dp_s(i, 7) = self%orb_dp_s(i, 7)+pixVal*v(2)*v(2)
         !2yz 
         self%orb_dp_s(i, 8) = self%orb_dp_s(i, 8) + 2.d0 * pixVal* v(2) * v(3)
         !z^2 
         self%orb_dp_s(i, 9) = self%orb_dp_s(i, 9) + pixVal*v(3)*v(3)
         !full beam integral for normalization
         self%orb_dp_s(i, 10) = self%orb_dp_s(i,10) + pixVal
         !if(self%myid == 38) then
         !  write(*,*) v(1), v(2), v(3), pixVal, self%orb_dp_s(i, 6)
         !end if
       end do
    end do

   call mpi_allreduce(MPI_IN_PLACE, self%orb_dp_s, size(self%orb_dp_s), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)

   do i = 1, self%ndet
      self%orb_dp_s(i,:) = self%orb_dp_s(i,:)*4*pi/real(self%beam(i)%p%info%npix)
   end do

!!$   if (self%tod%myid == 0) then
!!$      do i = 1, 10
!!$        write(*,*) self%orb_dp_s(1, i)
!!$      end do
!!$    end if

    !npipe s factors for 27M
!    self%orb_dp_s(:,1) = 0.005130801850847007
!    self%orb_dp_s(:,2) = -0.000516559072591428
!    self%orb_dp_s(:,3) = 0.995234628256561
!    self%orb_dp_s(:,4) = 0.00483461658765793
!    self%orb_dp_s(:,5) = 2.d0* -0.00043088175007651217
!    self%orb_dp_s(:,6) = 2.d0*0.0007072028589201922
!    self%orb_dp_s(:,7) = 0.0005094291355884364
!    self%orb_dp_s(:,8) = 2.d0* -0.00010373401957447322
!    self%orb_dp_s(:,9) = 0.9946559542767537
!    self%orb_dp_s(:,10) = 0.9927374627686433

    !dpc s factors for 27M
!    self%orb_dp_s(:,1) = 4.6679499857542042e-03
!    self%orb_dp_s(:,2) = 1.5034470382246318e-03  
!    self%orb_dp_s(:,3) = 9.9334531414858640e-01  
!    self%orb_dp_s(:,4) = 4.1262573481758826e-03
!    self%orb_dp_s(:,5) = 1.1183327794935952e-03  
!    self%orb_dp_s(:,6) = 8.0451444124375745e-04  
!    self%orb_dp_s(:,7) = 8.2453898616820286e-04
!    self%orb_dp_s(:,8) = 2.9241722968962171e-04  
!    self%orb_dp_s(:,9) = 9.9277120225676552e-01
!    self%orb_dp_s(:,10)= 9.9772199859110755e-01


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
    !        if True, comoputes relativistic correction
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

    real(dp)     :: b, x, q, b_dot, f, vp_ref(3), xx, v(3)
    integer(i4b) :: i, j, k, s_len, ntod, subsample
    real(dp), dimension(:), allocatable :: x_vec, y_vec

    subsample = 20

    f      = 1.d0; if (present(factor)) f = factor
    ntod   = size(s_dip)
    if ((ntod - 1) == 0) stop '!!!!!!!ntod = 1, orb dipole bug!!!!!!!'
    x      = h * nu/(k_B * T_CMB)
    q      = (x/2.d0)*(exp(x)+1)/(exp(x) -1)
    vp_ref = v_ref; if (present(v_ref_next)) vp_ref = v_ref_next ! Velocity for next PID; for linear interpolation

    if (.not. beam_4pi) then
       do i = 1, ntod !length of the tod
          xx       = real(i-1,dp)/real(ntod-1,dp)
          v        = (1.d0-xx) * v_ref + xx*vp_ref
          b_dot    = dot_product(v, P(:,i))/c
          s_dip(i) = b_dot
          if (relativistic) then
             b         = sqrt(sum(v**2))/c
             s_dip(i)  = s_dip(i) + q*(b_dot**2  - b**2/3. )
          end if
          s_dip(i)  = f * T_CMB * s_dip(i)
       end do
    else 
       s_len = ntod/subsample
       allocate(x_vec(s_len), y_vec(s_len))
       do k = 1, s_len !number of subsampled samples
          xx       = real(k-1,dp)/real(s_len-1,dp)
          v        = (1.d0-xx)*v_ref + xx*vp_ref
          j        = subsample * (k-1) + 1
          X_vec(k) = j
          y_vec(k) = self%compute_4pi_product(det, q, P(:,j), v) * f
       end do
       
       !spline the subsampled dipole to the full resolution
       call spline_simple(self%s, x_vec, y_vec, regular=.true.)

       do i = 1, s_len*subsample
         s_dip(i) = splint_simple(self%s, real(i,dp))
       end do

       !handle the last irregular length bit of each chunk as a special case 
       !so that we can use regular=true in the spline code for speed; use end velocity for this bin
       do j=s_len*subsample+1, ntod
         s_dip(j) = self%compute_4pi_product(det, q, P(:,j), vp_ref) 
       end do

       deallocate(x_vec, y_vec)
       call free_spline(self%s)
    end if

  end subroutine compute_CMB_dipole


  function compute_4pi_product(self, det, q, Pang, v_ref) result(prod)
    implicit none
    class(comm_orbdipole), intent(in) :: self
    integer(i4b),          intent(in) :: det
    real(dp),              intent(in) :: q, Pang(3), v_ref(3)
    real(dp)                          :: prod

    real(dp), dimension(3,3) :: rot_mat
    real(dp), dimension(3)   :: vnorm

    call compute_euler_matrix_zyz(-Pang(3), -Pang(2), -Pang(1), rot_mat)
    vnorm = matmul(rot_mat, v_ref) / c

    ! Equation C.5 in NPIPE paper
    prod = vnorm(1)*self%orb_dp_s(det,1)+ &
         & vnorm(2)*self%orb_dp_s(det,2)+ &
         & vnorm(3)*self%orb_dp_s(det,3)+ &
         &  q*(vnorm(1)*vnorm(1)*self%orb_dp_s(det,4) + &
         &     vnorm(1)*vnorm(2)*self%orb_dp_s(det,5) + &
         &     vnorm(1)*vnorm(3)*self%orb_dp_s(det,6) + &
         &     vnorm(2)*vnorm(2)*self%orb_dp_s(det,7) + &
         &     vnorm(2)*vnorm(3)*self%orb_dp_s(det,8) + &
         &     vnorm(3)*vnorm(3)*self%orb_dp_s(det,9))

    prod = T_CMB*prod/self%orb_dp_s(det,10)

  end function compute_4pi_product

end module comm_tod_orbdipole_mod
