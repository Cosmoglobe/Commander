module comm_F_int_2D_mod
  use comm_F_int_mod
  use comm_comp_mod
  use comm_bp_mod
  use spline_2D_mod
  use comm_utils
  implicit none

  private
  public comm_F_int_2D

  type, extends (comm_F_int) :: comm_F_int_2D
     real(dp), allocatable, dimension(:)       :: x, y
     real(dp), allocatable, dimension(:,:,:,:) :: coeff
   contains
     ! Data procedures
     procedure :: eval => evalIntF
  end type comm_F_int_2D

  interface comm_F_int_2D
     procedure constructor
  end interface comm_F_int_2D

  integer(i4b) :: n = 100
  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(comp, bp, pol)
    implicit none
    class(comm_comp),     intent(in) :: comp
    class(comm_bp),       intent(in) :: bp
    integer(i4b),         intent(in) :: pol
    class(comm_F_int_2D), pointer    :: constructor

    integer(i4b) :: i, j, k, m, ierr
    real(dp)     :: t1, t2, t3, t4
    real(dp), allocatable, dimension(:)   :: s
    real(dp), allocatable, dimension(:,:) :: F, Fp
    
    allocate(constructor)

    call wall_time(t1)

    m = bp%n
    allocate(constructor%x(n), constructor%y(n), F(n,n), Fp(n,n), s(m))

    do i = 1, n
       constructor%x(i) = comp%p_uni(1,1) + (comp%p_uni(2,1)-comp%p_uni(1,1))/(n-1) * (i-1)
       constructor%y(i) = comp%p_uni(1,2) + (comp%p_uni(2,2)-comp%p_uni(1,2))/(n-1) * (i-1)
    end do

    ! Evaluate the bandpass integrated SED over all relevant parameters
    F = 0.d0
    call wall_time(t3)
    do i = 1+comp%myid, n, comp%numprocs
       do j = 1, n
          do k = 1, m
             s(k) = comp%S(nu=bp%nu(k), pol=pol, theta=[constructor%x(i),constructor%y(j)])
          end do
          F(i,j) = bp%SED2F(s)
       end do
    end do
    call wall_time(t4)
    !if (comp%myid == 0) write(*,*) 'a = ', real(t4-t3,sp)
    call wall_time(t3)
    call mpi_allreduce(F, Fp,            n*n, MPI_DOUBLE_PRECISION, MPI_SUM, comp%comm, ierr)
    call wall_time(t4)
    !if (comp%myid == 0) write(*,*) 'b = ', real(t4-t3,sp)

    ! Precompute spline object
    call wall_time(t3)
    allocate(constructor%coeff(4,4,n,n))
    call splie2_full_precomp_mpi(comp%comm, constructor%x, constructor%y, Fp, constructor%coeff)
    call wall_time(t4)
    !if (comp%myid == 0) write(*,*) 'c = ', real(t4-t3,sp)
    
    call wall_time(t2)
    !if (comp%myid == 0) write(*,*) 'sum = ', sum(abs(constructor%coeff)), real(t2-t1,sp)
    !call mpi_finalize(i)
    !stop

    deallocate(F, Fp, s)
    
  end function constructor

  
  ! Evaluate SED in brightness temperature normalized to reference frequency
  function evalIntF(self, theta)
    class(comm_F_int_2D),                intent(in) :: self
    real(dp),             dimension(1:), intent(in) :: theta
    real(dp)                                        :: evalIntF
    evalIntF = splin2_full_precomp(self%x, self%y, self%coeff, theta(1), theta(2)) 
  end function evalIntF

end module comm_F_int_2D_mod
