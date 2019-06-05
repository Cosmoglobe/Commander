module comm_F_int_0D_mod
  use comm_F_int_mod
  use comm_comp_mod
  use comm_bp_mod
  use comm_utils
  implicit none

  private
  public comm_F_int_0D

  type, extends (comm_F_int) :: comm_F_int_0D
     real(dp) :: f_precomp
   contains
     ! Data procedures
     procedure :: eval => evalIntF
  end type comm_F_int_0D

  interface comm_F_int_0D
     procedure constructor
  end interface comm_F_int_0D

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(comp, bp, pol, f_precomp)
    implicit none
    class(comm_comp),     intent(in)           :: comp
    class(comm_bp),       intent(in)           :: bp
    integer(i4b),         intent(in)           :: pol
    real(dp),             intent(in), optional :: f_precomp
    class(comm_F_int_0D), pointer              :: constructor

    integer(i4b) :: j, m
    real(dp), allocatable, dimension(:) :: s
    
    allocate(constructor)

    if (present(f_precomp)) then
       constructor%f_precomp = f_precomp
    else
       ! Evaluate the bandpass integrated SED
       m = bp%n
       allocate(s(m))
       do j = 1, m
          s(j) = comp%S(nu=bp%nu(j))
       end do
       constructor%f_precomp = bp%SED2F(s)
       deallocate(s)
    end if

  end function constructor

  
  ! Evaluate SED in brightness temperature normalized to reference frequency
  function evalIntF(self, theta)
    class(comm_F_int_0D),                intent(in) :: self
    real(dp),             dimension(1:), intent(in) :: theta
    real(dp)                                        :: evalIntF
    evalIntF = self%f_precomp
  end function evalIntF

end module comm_F_int_0D_mod
