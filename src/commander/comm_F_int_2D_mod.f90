module comm_F_int_2D_mod
  use comm_F_int_mod
  use comm_comp_mod
  use comm_bp_mod
  use spline_2D_mod
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
  function constructor(comp, bp)
    implicit none
    class(comm_comp),     intent(in) :: comp
    class(comm_bp),       intent(in) :: bp
    class(comm_F_int_2D), pointer    :: constructor

    integer(i4b) :: i, j, k, m
    real(dp), allocatable, dimension(:)   :: s
    real(dp), allocatable, dimension(:,:) :: F
    
    allocate(constructor)

    m = bp%n
    allocate(constructor%x(n), constructor%y(n), F(n,n), s(m))

    ! Evaluate the bandpass integrated SED over all relevant parameters
    do i = 1, n
       constructor%x(i) = comp%p_uni(1,1) + (comp%p_uni(2,1)-comp%p_uni(1,1))/(n-1) * (i-1)
       do j = 1, n
          constructor%y(j) = comp%p_uni(1,2) + (comp%p_uni(2,2)-comp%p_uni(1,2))/(n-1) * (j-1)
          do k = 1, m
             s(k) = comp%S(nu=bp%nu(k), theta=[constructor%x(i),constructor%y(j)])
          end do
          F(i,j) = bp%SED2F(s)
       end do
    end do

    ! Precompute spline object
    allocate(constructor%coeff(4,4,n,n))
    call splie2_full_precomp(constructor%x, constructor%y, F, constructor%coeff)
    
    deallocate(F, s)
    
  end function constructor

  
  ! Evaluate SED in brightness temperature normalized to reference frequency
  function evalIntF(self, theta)
    class(comm_F_int_2D),                intent(in) :: self
    real(dp),             dimension(1:), intent(in) :: theta
    real(dp)                                        :: evalIntF
    evalIntF = splin2_full_precomp(self%x, self%y, self%coeff, theta(1), theta(2)) 
  end function evalIntF

end module comm_F_int_2D_mod
