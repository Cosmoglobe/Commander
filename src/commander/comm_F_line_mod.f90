module comm_F_line_mod
  use comm_F_int_mod
  use comm_comp_mod
  use comm_bp_mod
  use comm_utils
  implicit none

  private
  public comm_F_line

  type, extends (comm_F_int) :: comm_F_line
     class(comm_comp), pointer :: comp
     logical(lgt) :: active
     real(dp)     :: f_precomp
     integer(i4b) :: ind
   contains
     ! Data procedures
     procedure :: eval => evalIntF
     procedure :: update => updateIntF
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
    class(comm_comp),                 intent(in), target :: comp
    class(comm_bp),                   intent(in), target :: bp
    logical(lgt),                     intent(in)         :: active
    real(dp),                         intent(in)         :: scale 
    integer(i4b),                     intent(in)         :: ind
    class(comm_F_line), pointer                          :: constructor

    allocate(constructor)
    constructor%comp   => comp
    constructor%bp     => bp
    constructor%active =  active

    if (constructor%active) then
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

  ! Compute/update integration look-up tables
  subroutine updateIntF(self, f_precomp, pol)
    class(comm_F_line),   intent(inout)           :: self
    real(dp),             intent(in),    optional :: f_precomp
    integer(i4b),         intent(in),    optional :: pol

  end subroutine updateIntF


end module comm_F_line_mod
