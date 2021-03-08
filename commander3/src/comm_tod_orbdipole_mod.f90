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
  use comm_tod_mod
  use spline_1D_mod
  implicit none

  private 
  public comm_orbdipole, orbdipole_pointer

  type :: comm_orbdipole
    real(dp),           allocatable, dimension(:, :) :: orb_dp_s !precomputed s integrals for orbital dipole sidelobe term
    class(map_ptr), dimension(:), allocatable :: beam
    class(comm_tod), pointer :: tod
    type(spline_type) :: s
    integer(i4b) :: subsample

  contains
    procedure :: precompute_orb_dp_s
    procedure :: compute_orbital_dipole_4pi
    procedure :: compute_orbital_dipole_4pi_diff
    procedure :: compute_orbital_dipole_pencil
    procedure :: compute_4pi_product
    procedure :: compute_solar_dipole_4pi
    procedure :: compute_solar_dipole_pencil
    procedure :: compute_4pi_product_sol

  end type comm_orbdipole

  interface comm_orbdipole
    procedure constructor
  end interface comm_orbdipole

  type orbdipole_pointer
    class(comm_orbdipole), pointer :: p => null()
  end type orbdipole_pointer

contains

  function constructor(tod, beam)
    implicit none
    class(comm_tod), target :: tod
    class(map_ptr), dimension(:), target :: beam
    class(comm_orbdipole), pointer :: constructor
    integer(i4b) :: i
    allocate(constructor)

    constructor%tod => tod
    allocate(constructor%beam(constructor%tod%ndet))
    do i=1, constructor%tod%ndet
      constructor%beam(i)%p => beam(i)%p
    end do

    call constructor%precompute_orb_dp_s()

    constructor%subsample = 20

  end function constructor

  subroutine precompute_orb_dp_s(self) 
    implicit none
    class(comm_orbdipole),              intent(inout) :: self

    integer(i4b) :: i, j, ierr
    real(dp), dimension(3) :: v
    real(dp) :: pixVal

    allocate(self%orb_dp_s(self%tod%ndet, 10)) 
    do i = 1, self%tod%ndet
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
         self%orb_dp_s(i, 7) =self%orb_dp_s(i, 7)+pixVal*v(2)*v(2)
         !2yz 
         self%orb_dp_s(i, 8) = self%orb_dp_s(i, 8) + 2.d0 * pixVal* v(2) * v(3)
         !z^2 
         self%orb_dp_s(i, 9) =self%orb_dp_s(i, 9) + pixVal*v(3)*v(3)
         !full beam integral for normalization
         self%orb_dp_s(i, 10) = self%orb_dp_s(i,10) + pixVal
         !if(self%myid == 38) then
         !  write(*,*) v(1), v(2), v(3), pixVal, self%orb_dp_s(i, 6)
         !end if
       end do
    end do

   call mpi_allreduce(MPI_IN_PLACE, self%orb_dp_s, size(self%orb_dp_s), MPI_DOUBLE_PRECISION, MPI_SUM, self%tod%info%comm, ierr)

   do i = 1, self%tod%ndet
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

  subroutine compute_orbital_dipole_pencil(self, ind, pix, psi, s_orb, factor)
    implicit none
    class(comm_orbdipole),               intent(in)  :: self
    integer(i4b),                        intent(in)  :: ind !scan nr/index
    integer(i4b),        dimension(:),   intent(in)  :: pix, psi
    real(sp),            dimension(:,:), intent(out) :: s_orb
    real(dp),               intent(in), optional     :: factor
    real(dp) :: b, x, q, b_dot, f
    integer(i4b) :: i, j

    if (present(factor)) then 
      f = factor
    else
      f = 1
    end if

    b = sqrt(sum(self%tod%scans(ind)%v_sun**2))/c   ! beta for the given scan
    x = h * self%tod%central_freq/(k_B * T_CMB)
    q = (x/2.d0)*(exp(x)+1)/(exp(x) -1)


    do i = 1,self%tod%ndet
       if (.not. self%tod%scans(ind)%d(i)%accept) then
         s_orb(:,i) = 0
         cycle
       end if
       do j=1,self%tod%scans(ind)%ntod !length of the tod
          b_dot = dot_product(self%tod%scans(ind)%v_sun, self%tod%pix2vec(:,pix(j)))/c
          s_orb(j,i) = f*T_CMB * (b_dot + q*(b_dot**2 - b**2/3.))
       end do
    end do

  end subroutine compute_orbital_dipole_pencil

  subroutine compute_solar_dipole_pencil(self, ind, pix, s_sol, factor)
    implicit none
    class(comm_orbdipole),               intent(in)  :: self
    integer(i4b),                        intent(in)  :: ind !scan nr/index
    integer(i4b),        dimension(:),   intent(in)  :: pix
    real(sp),            dimension(:),   intent(out) :: s_sol
    real(sp),                            intent(in), optional :: factor
    real(dp) :: b, x, q, b_dot, phi, theta
    real(sp) :: f
    real(dp), dimension(3) :: vnorm
    integer(i4b) :: i, j

    if (.not. self%tod%scans(ind)%d(i)%accept) then
       s_sol = 0
       return
    end if

    f = 1.; if (present(factor)) f = factor
    phi   = 4.607145626489432  ! 263.97*pi/180
    theta = 0.7278022980816355 ! (90-48.3)*pi/180
    vnorm = (/ sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) /)

    do j=1,self%tod%scans(ind)%ntod !length of the tod
       b_dot = dot_product(vnorm, self%tod%pix2vec(:,pix(j)))
       s_sol(j) = f * T_CMB_DIP * b_dot
    end do
  end subroutine compute_solar_dipole_pencil

  function compute_4pi_product(self, p, psiInd, chunkInd, i, q) result(prod)
    implicit none
    class(comm_orbdipole), intent(in) :: self
    integer(i4b),          intent(in) :: p, psiInd, chunkInd, i
    real(dp),              intent(in) :: q
    real(dp)                          :: prod

    real(dp) :: theta, psi_d, phi
    real(dp), dimension(3,3) :: rot_mat
    real(dp), dimension(3)   :: vnorm

    theta = self%tod%ind2ang(1,p)
    phi   = self%tod%ind2ang(2,p)
    psi_d = self%tod%psi(psiInd)
    !write(*,*), j, phi, theta, psi_d, rot_mat 
    call compute_euler_matrix_zyz(-psi_d, -theta, -phi, rot_mat)
    vnorm = matmul(rot_mat, self%tod%scans(chunkInd)%v_sun)
    vnorm = vnorm / c
    ! Equation C.5 in NPIPE paper
    prod = vnorm(1)*self%orb_dp_s(i,1)+ &
          &vnorm(2)*self%orb_dp_s(i,2)+ &
          &vnorm(3)*self%orb_dp_s(i,3)+ &
          & q*(vnorm(1)*vnorm(1)*self%orb_dp_s(i,4) + &
          &    vnorm(1)*vnorm(2)*self%orb_dp_s(i,5) + &
          &    vnorm(1)*vnorm(3)*self%orb_dp_s(i,6) + &
          &    vnorm(2)*vnorm(2)*self%orb_dp_s(i,7) + &
          &    vnorm(2)*vnorm(3)*self%orb_dp_s(i,8) + &
          &    vnorm(3)*vnorm(3)*self%orb_dp_s(i,9))

    prod = T_CMB*prod/self%orb_dp_s(i,10)

  end function compute_4pi_product

  function compute_4pi_product_sol(self, p, psiInd, chunkInd, i, q) result(prod)
    implicit none
    class(comm_orbdipole), intent(in) :: self
    integer(i4b),          intent(in) :: p, psiInd, chunkInd, i
    real(dp),              intent(in) :: q
    real(dp)                          :: prod

    real(dp) :: theta, psi_d, phi
    real(dp), dimension(3,3) :: rot_mat
    real(dp), dimension(3)   :: vnorm
    real(dp), dimension(3)   :: T_dip

    phi   = 4.607145626489432  ! 263.97*pi/180
    theta = 0.7278022980816355 ! (90-48.3)*pi/180
    vnorm = (/ sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) /)

    theta = self%tod%ind2ang(1,p)
    phi   = self%tod%ind2ang(2,p)
    psi_d = self%tod%psi(psiInd)
    !write(*,*), j, phi, theta, psi_d, rot_mat 
    call compute_euler_matrix_zyz(-psi_d, -theta, -phi, rot_mat)
    vnorm = matmul(rot_mat, vnorm)
    ! Equation C.5 in NPIPE paper
    prod = vnorm(1)*self%orb_dp_s(i,1)+ &
          &vnorm(2)*self%orb_dp_s(i,2)+ &
          &vnorm(3)*self%orb_dp_s(i,3)+ &
          & q*(vnorm(1)*vnorm(1)*self%orb_dp_s(i,4) + &
          &    vnorm(1)*vnorm(2)*self%orb_dp_s(i,5) + &
          &    vnorm(1)*vnorm(3)*self%orb_dp_s(i,6) + &
          &    vnorm(2)*vnorm(2)*self%orb_dp_s(i,7) + &
          &    vnorm(2)*vnorm(3)*self%orb_dp_s(i,8) + &
          &    vnorm(3)*vnorm(3)*self%orb_dp_s(i,9))

    prod = T_CMB_DIP*prod/self%orb_dp_s(i,10)

  end function compute_4pi_product_sol

  ! Compute map with white noise assumption from correlated noise 
  ! corrected and calibrated data, d' = (d-n_corr-n_temp)/gain 
  subroutine compute_orbital_dipole_4pi(self, ind, pix, psi, s_orb, factor)
    implicit none
    class(comm_orbdipole),               intent(inout)  :: self
    integer(i4b),                        intent(in)  :: ind !scan nr/index
    integer(i4b),        dimension(:,:), intent(in)  :: pix, psi
    real(sp),            dimension(:,:), intent(out) :: s_orb
    real(dp),               intent(in), optional     :: factor
    integer(i4b)         :: i, j, p, s_len, k, psiInd
    real(dp)             :: x, T_0, q, pix_dir(3), b, b_dot, summation, j_real, f
    real(dp), parameter  :: h = 6.62607015d-34   ! Planck's constant [Js]
    real(dp), dimension(:), allocatable :: x_vec, y_vec

    f = 1.; if (present(factor)) f = factor
    !these are the npipe paper definitions
    !TODO: maybe also use the bandpass shift to modify the central frequency? 
    !will that matter?
    x = h * self%tod%central_freq/(k_B * T_CMB)
    q = (x/2.d0)*(exp(x)+1)/(exp(x) -1)

    do i = 1,self%tod%ndet
       if (.not. self%tod%scans(ind)%d(i)%accept) cycle
       s_len = self%tod%scans(ind)%ntod/self%subsample
       allocate(x_vec(s_len), y_vec(s_len))
       do k=0,s_len-1 !number of subsampled samples
          j = (self%subsample * k) + 1

          p      = self%tod%pix2ind(pix(j,i))
          psiInd = psi(j,i)

          summation = f*self%compute_4pi_product(p, psiInd, ind, i, q)
          x_vec(k+1) = j
          y_vec(k+1) = summation
       end do
       
       !spline the subsampled dipole to the full resolution
       call spline_simple(self%s, x_vec, y_vec, regular=.true.)

       do j=1, s_len*self%subsample
         j_real = j
         s_orb(j,i) = splint_simple(self%s, j_real)
       end do

       !handle the last irregular length bit of each chunk as a special case 
       !so that we can use regular=true in the spline code for speed
       do j=s_len*self%subsample, self%tod%scans(ind)%ntod
         p      = self%tod%pix2ind(pix(j,i))
         psiInd = psi(j,i)
         
         s_orb(j,i) = self%compute_4pi_product(p, psiInd, ind, i, q)

       end do

       deallocate(x_vec, y_vec)
   end do

   call free_spline(self%s)
   
   !if (trim(self%tod%label(1)) == '27M') then
   ! close(58)
   ! close(59)
   !end if
  end subroutine compute_orbital_dipole_4pi

  subroutine compute_orbital_dipole_4pi_diff(self, ind, pix, psi, s_orb, i, factor)
    implicit none
    class(comm_orbdipole),               intent(inout)  :: self
    integer(i4b),                        intent(in)  :: ind !scan nr/index
    integer(i4b),                        intent(in)  :: i   !detector
    integer(i4b),        dimension(:), intent(in)  :: pix, psi
    real(sp),            dimension(:,:), intent(out) :: s_orb
    real(dp),               intent(in), optional     :: factor
    integer(i4b)         :: j, p, s_len, k, psiInd
    real(dp)             :: x, T_0, q, pix_dir(3), b, b_dot, summation, j_real, f
    real(dp), parameter  :: h = 6.62607015d-34   ! Planck's constant [Js]
    real(dp), dimension(:), allocatable :: x_vec, y_vec

    f = 1.d0; if (present(factor)) f = factor
    !these are the npipe paper definitions
    !TODO: maybe also use the bandpass shift to modify the central frequency? 
    !will that matter?
    x = h * self%tod%central_freq/(k_B * T_CMB)
    q = (x/2.d0)*(exp(x)+1.d0)/(exp(x) -1.d0)

    s_len = self%tod%scans(ind)%ntod/self%subsample
    allocate(x_vec(s_len), y_vec(s_len))
    do k=0,s_len-1 !number of subsampled samples
       j = (self%subsample * k) + 1

       p      = self%tod%pix2ind(pix(j))
       psiInd = psi(j)

       summation = f*self%compute_4pi_product(p, psiInd, ind, i, q)
       x_vec(k+1) = j
       y_vec(k+1) = summation
    end do
    
    !spline the subsampled dipole to the full resolution
    call spline_simple(self%s, x_vec, y_vec, regular=.true.)

    do j=1, s_len*self%subsample
      j_real = j
      s_orb(j,i) = splint_simple(self%s, j_real)
    end do

    !handle the last irregular length bit of each chunk as a special case 
    !so that we can use regular=true in the spline code for speed
    do j=s_len*self%subsample, self%tod%scans(ind)%ntod
      p      = self%tod%pix2ind(pix(j))
      psiInd = psi(j)
      
      s_orb(j,i) = self%compute_4pi_product(p, psiInd, ind, i, q)

    end do

    deallocate(x_vec, y_vec)

    call free_spline(self%s)
   
  end subroutine compute_orbital_dipole_4pi_diff

  subroutine compute_solar_dipole_4pi(self, ind, pix, psi, s_sol)
    implicit none
    class(comm_orbdipole),               intent(inout)  :: self
    integer(i4b),                        intent(in)  :: ind !scan nr/index
    integer(i4b),        dimension(:,:), intent(in)  :: pix, psi
    real(sp),            dimension(:,:), intent(out) :: s_sol
    integer(i4b)         :: i, j, p, s_len, k, psiInd
    real(dp)             :: x, T_0, q, pix_dir(3), b, b_dot, summation, j_real
    real(dp), parameter  :: h = 6.62607015d-34   ! Planck's constant [Js]
    real(dp), dimension(:), allocatable :: x_vec, y_vec

    !these are the npipe paper definitions
    !TODO: maybe also use the bandpass shift to modify the central frequency? 
    !will that matter?
    x = h * self%tod%central_freq/(k_B * T_CMB)
    q = (x/2.d0)*(exp(x)+1)/(exp(x) -1)

    do i = 1,self%tod%ndet
       if (.not. self%tod%scans(ind)%d(i)%accept) cycle
       s_len = self%tod%scans(ind)%ntod/self%subsample
       allocate(x_vec(s_len), y_vec(s_len))
       do k=0,s_len-1 !number of subsampled samples
          j = (self%subsample * k) + 1

          p      = self%tod%pix2ind(pix(j,i))
          psiInd = psi(j,i)

          summation = self%compute_4pi_product_sol(p, psiInd, ind, i, q)
          x_vec(k+1) = j
          y_vec(k+1) = summation
       end do
       
       !spline the subsampled dipole to the full resolution
       call spline_simple(self%s, x_vec, y_vec, regular=.true.)

       do j=1, s_len*self%subsample
         j_real = j
         s_sol(j,i) = splint_simple(self%s, j_real)
       end do

       !handle the last irregular length bit of each chunk as a special case 
       !so that we can use regular=true in the spline code for speed
       do j=s_len*self%subsample, self%tod%scans(ind)%ntod
         p      = self%tod%pix2ind(pix(j,i))
         psiInd = psi(j,i)
         
         s_sol(j,i) = self%compute_4pi_product_sol(p, psiInd, ind, i, q)

       end do

       deallocate(x_vec, y_vec)
   end do

   call free_spline(self%s)
   
  end subroutine compute_solar_dipole_4pi
end module comm_tod_orbdipole_mod
