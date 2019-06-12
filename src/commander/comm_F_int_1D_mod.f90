module comm_F_int_1D_mod
  use comm_F_int_mod
  use comm_comp_mod
  use comm_bp_mod
  use spline_1D_mod
  use comm_utils
  implicit none

  private
  public comm_F_int_1D

  type, extends (comm_F_int) :: comm_F_int_1D
     type(spline_type) :: s
     class(comm_comp), pointer :: comp
   contains
     ! Data procedures
     procedure :: eval => evalIntF
     procedure :: update => updateIntF
  end type comm_F_int_1D

  interface comm_F_int_1D
     procedure constructor
  end interface comm_F_int_1D

  integer(i4b) :: n = 1000
  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(comp, bp, pol)
    implicit none
    class(comm_comp),     intent(in), target  :: comp
    class(comm_bp),       intent(in), target  :: bp
    integer(i4b),         intent(in)          :: pol
    class(comm_F_int_1D),             pointer :: constructor

    integer(i4b) :: i, j, m, ierr
    real(dp), allocatable, dimension(:) :: theta, F, s
    
    allocate(constructor)
    constructor%comp => comp
    constructor%bp   => bp

    call constructor%update(pol=pol)

  end function constructor

  
  ! Evaluate SED in brightness temperature normalized to reference frequency
  function evalIntF(self, theta)
    class(comm_F_int_1D),                intent(in) :: self
    real(dp),             dimension(1:), intent(in) :: theta
    real(dp)                                        :: evalIntF
    evalIntF = splint_simple(self%s, theta(1))
  end function evalIntF

  ! Compute/update integration look-up tables
  subroutine updateIntF(self, f_precomp, pol)
    class(comm_F_int_1D), intent(inout)           :: self
    real(dp),             intent(in),    optional :: f_precomp
    integer(i4b),         intent(in),    optional :: pol

    integer(i4b) :: i, j, m, ierr
    real(dp), allocatable, dimension(:) :: theta, F, s

    m = self%bp%n
    allocate(theta(n), F(n), s(m))

    do i = 1, n
       theta(i) = self%comp%p_uni(1,1) + (self%comp%p_uni(2,1)-self%comp%p_uni(1,1))/(n-1) * (i-1)
    end do

    ! Evaluate the bandpass integrated SED over all relevant parameters
    F = 0.d0
    do i = 1+self%comp%myid, n, self%comp%numprocs
       do j = 1, m
          s(j) = self%comp%S(nu=self%bp%nu(j), pol=pol, theta=theta(i:i))
       end do
       F(i) = self%bp%SED2F(s)
    end do
    call mpi_allreduce(MPI_IN_PLACE, F, n, MPI_DOUBLE_PRECISION, MPI_SUM, self%comp%comm, ierr)

    ! Precompute spline object
    call spline_simple(self%s, theta, F, regular=.true.)

    deallocate(theta, F, s)

  end subroutine updateIntF


end module comm_F_int_1D_mod
