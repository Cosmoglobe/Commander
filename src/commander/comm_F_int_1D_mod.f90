module comm_F_int_1D_mod
  use comm_F_int_mod
  use comm_comp_mod
  use comm_bp_mod
  use spline_1D_mod
  implicit none

  private
  public comm_F_int_1D

  type, extends (comm_F_int) :: comm_F_int_1D
     type(spline_type) :: s
   contains
     ! Data procedures
     procedure :: eval => evalIntF
  end type comm_F_int_1D

  interface comm_F_int_1D
     procedure constructor
  end interface comm_F_int_1D

  integer(i4b) :: n = 1000
  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(comp, bp)
    implicit none
    class(comm_comp),     intent(in) :: comp
    class(comm_bp),       intent(in) :: bp
    class(comm_F_int_1D), pointer    :: constructor

    integer(i4b) :: i, j, m
    real(dp), allocatable, dimension(:) :: theta, F, s
    
    allocate(constructor)

    m = bp%n
    allocate(theta(n), F(n), s(m))

    ! Evaluate the bandpass integrated SED over all relevant parameters
    do i = 1, n
       theta(i) = comp%p_uni(1,1) + (comp%p_uni(2,1)-comp%p_uni(1,1))/(n-1) * (i-1)
       do j = 1, m
          s(j) = comp%S(nu=bp%nu(j), theta=theta(i:i))
       end do
       F(i) = bp%SED2F(s)
    end do

    ! Precompute spline object
    call spline_simple(constructor%s, theta, F, regular=.true.)
    
    deallocate(theta, F, s)
    
  end function constructor

  
  ! Evaluate SED in brightness temperature normalized to reference frequency
  function evalIntF(self, theta)
    class(comm_F_int_1D),                intent(in) :: self
    real(dp),             dimension(1:), intent(in) :: theta
    real(dp)                                        :: evalIntF
    evalIntF = splint_simple(self%s, theta(1))
  end function evalIntF

end module comm_F_int_1D_mod
