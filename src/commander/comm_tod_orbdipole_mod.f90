module comm_tod_orbdipole_mod
  use comm_utils
  use comm_defs
  use comm_map_mod
  use comm_tod_mod
  implicit none

  private 
  public comm_orbdipole, orbdipole_pointer

  type :: comm_orbdipole
    real(dp),           allocatable, dimension(:, :) :: orb_dp_s !precomputed s integrals for orbital dipole sidelobe term
    class(map_ptr), dimension(:), allocatable :: beam
    class(comm_tod), pointer :: tod

  contains
    procedure precompute_orb_dp_s
    procedure compute_orbital_dipole_4pi
    procedure compute_orbital_dipole_pencil

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
      do j = 0, self%tod%info%np-1
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

   if (self%tod%myid == 0) then
      do i = 1, 10
        write(*,*) self%orb_dp_s(1, i)
      end do
    end if

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

  subroutine compute_orbital_dipole_pencil(self, ind, pix, psi, s_orb)
    implicit none
    class(comm_orbdipole),               intent(in)  :: self
    integer(i4b),                        intent(in)  :: ind !scan nr/index
    integer(i4b),        dimension(:,:), intent(in)  :: pix, psi
    real(sp),            dimension(:,:), intent(out) :: s_orb
    real(dp) :: b, x, q, b_dot
    integer(i4b) :: i, j

    b = sqrt(sum(self%tod%scans(ind)%v_sun**2))/c   ! beta for the given scan
    x = h * self%tod%central_freq/(k_B * T_CMB)
    q = (x/2.d0)*(exp(x)+1)/(exp(x) -1)


    do i = 1,self%tod%ndet
       if (.not. self%tod%scans(ind)%d(i)%accept) cycle
       do j=1,self%tod%scans(ind)%ntod !length of the tod
          b_dot = dot_product(self%tod%scans(ind)%v_sun, self%tod%pix2vec(:,pix(j,i)))/c
          !s_orb(j,i) = real(T_CMB  * b_dot,sp) !* self%tod%mb_eff(i) !only dipole,
          !1.d6 to make it uK, as [T_CMB] = K
          !s_orb(j,i) = T_CMB  * 1.d6 * b_dot !only dipole, 1.d6 to make it uK,
          !as [T_CMB] = K
          !s_orb(j,i) = T_CMB * 1.d6 * (b_dot + q*b_dot**2) ! with quadrupole
          s_orb(j,i) = T_CMB * (b_dot + q*(b_dot**2 - b**2/3.)) ! net zero monopole
        end do
    end do

  end subroutine compute_orbital_dipole_pencil


  ! Compute map with white noise assumption from correlated noise 
  ! corrected and calibrated data, d' = (d-n_corr-n_temp)/gain 
  subroutine compute_orbital_dipole_4pi(self, ind, pix, psi, s_orb)
    implicit none
    class(comm_orbdipole),               intent(in)  :: self
    integer(i4b),                        intent(in)  :: ind !scan nr/index
    integer(i4b),        dimension(:,:), intent(in)  :: pix, psi
    real(sp),            dimension(:,:), intent(out) :: s_orb
    integer(i4b)         :: i, j, p
    real(dp)             :: x, T_0, q, pix_dir(3), b, b_dot, summation
    real(dp)             :: theta, phi, psi_d
    real(dp), dimension(3,3) :: rot_mat
    real(dp), dimension(3)   :: vnorm
    real(dp), parameter  :: h = 6.62607015d-34   ! Planck's constant [Js]

    !these are the npipe paper definitions
    !TODO: maybe also use the bandpass shift to modify the central frequency? 
    !will that matter?
    x = h * self%tod%central_freq/(k_B * T_CMB)
    q = (x/2.d0)*(exp(x)+1)/(exp(x) -1)

    !if (trim(self%tod%label(1)) == '27M') then 
    !  open(58,file='orb_27M_npipe.dat')
    !  if (self%tod%scans(ind)%chunk_num == 27) write(58,*) " SCET    THETA    PHI    PSI    TOD    ORB_DP    ORB_FSL"
    !  open(59,file='orb_27M_debug.dat')
    !end if
    do i = 1,self%tod%ndet 
       if (.not. self%tod%scans(ind)%d(i)%accept) cycle
       do j=1,self%tod%scans(ind)%ntod !length of the tod

          call pix2ang_ring(self%tod%info%nside, pix(j,i), theta, phi)
          !rotate v_sun into frame where pointing is along z axis
          !write(*,*) -phi, -theta, -self%tod%psi(psi(j,i)), psi(j,i)
          p     = self%tod%pix2ind(pix(j,i))
          theta = self%tod%ind2ang(1,p)
          phi   = self%tod%ind2ang(2,p)
          psi_d = self%tod%psi(psi(j,i))
          !write(*,*), j, phi, theta, psi_d, rot_mat 
          !call compute_euler_matrix_zyz(-phi, -theta, -psi_d, rot_mat)
          call compute_euler_matrix_zyz(-psi_d, -theta, -phi, rot_mat)
          vnorm = matmul(rot_mat, self%tod%scans(ind)%v_sun)
          vnorm = vnorm / c
          summation = vnorm(1)*self%orb_dp_s(i,1)+vnorm(2)*self%orb_dp_s(i,2)+& 
          &vnorm(3)*self%orb_dp_s(i,3)+q*(vnorm(1)*vnorm(1)*self%orb_dp_s(i,4)+&
          &vnorm(1)*vnorm(2)*self%orb_dp_s(i,5) + vnorm(1)*vnorm(3)* &
          &self%orb_dp_s(i,6) + vnorm(2)*vnorm(2)*self%orb_dp_s(i,7) + &
          &vnorm(2)*vnorm(3)*self%orb_dp_s(i,8) + vnorm(3)*vnorm(3)*&
          &self%orb_dp_s(i,9))
          !if (trim(self%tod%label(i)) == '27M' .and. self%tod%scans(ind)%chunk_num == 27) write(59,*) self%tod%scans(ind)%t0(3)/1000000.d0 + real(j-1)/(real(self%tod%samprate)), theta, phi, psi_d, self%tod%scans(ind)%v_sun, vnorm, self%tod%scans(ind)%d(i)%tod(j), s_orb(j,i), T_CMB*summation
          !if (trim(self%tod%label(i)) == '27M' .and. self%tod%scans(ind)%chunk_num == 27) write(58,"(ES25.18,  ES15.7,  ES15.7,  ES15.7,  ES15.7,  ES15.7,  ES15.7)") self%tod%scans(ind)%t0(2)/2**16 + real(j-1)/(real(self%tod%samprate)), theta, phi, psi_d, self%tod%scans(ind)%d(i)%tod(j)/self%tod%scans(ind)%d(i)%gain, s_orb(j,i), T_CMB*summation/self%orb_dp_s(i, 10)

          !s_orb(j,i) = s_orb(j,i) + T_CMB *summation
          s_orb(j,i) = T_CMB*summation/self%orb_dp_s(i,10)

       end do
   end do
   !if (trim(self%tod%label(1)) == '27M') then
   ! close(58)
   ! close(59)
   !end if
  end subroutine compute_orbital_dipole_4pi
end module comm_tod_orbdipole_mod
