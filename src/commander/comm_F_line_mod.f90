module comm_F_line_mod
  use comm_F_int_mod
  use comm_comp_mod
  use comm_bp_mod
  use comm_utils
  implicit none

  private
  public comm_F_line

  type, extends (comm_F_int) :: comm_F_line
     logical(lgt) :: active
     real(dp)     :: f_precomp
     integer(i4b) :: ind
   contains
     ! Data procedures
     procedure :: eval => evalIntF
  end type comm_F_line

  interface comm_F_line
     procedure constructor
  end interface comm_F_line

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(comp, bp, active, scale, ind)
    implicit none
    class(comm_comp),                 intent(in) :: comp
    class(comm_bp),                   intent(in) :: bp
    logical(lgt),                     intent(in) :: active
    real(dp),                         intent(in) :: scale
    integer(i4b),                     intent(in) :: ind
    class(comm_F_line), pointer                  :: constructor

    allocate(constructor)
    constructor%active = active
    if (active) then
       constructor%f_precomp = scale
       constructor%ind       = ind
    else
       constructor%f_precomp = 0.d0
    end if

  end function constructor

  
  ! Evaluate SED in brightness temperature normalized to reference frequency
  function evalIntF(self, theta)
    class(comm_F_line),                  intent(in) :: self
    real(dp),             dimension(1:), intent(in) :: theta
    real(dp)                                        :: evalIntF
    if (self%active) then
       if (self%ind > 0) then
          evalIntF = self%f_precomp * theta(self%ind)
       else
          evalIntF = self%f_precomp
       end if
    else 
       evalIntF = 0.d0
    end if
  end function evalIntF

end module comm_F_line_mod
