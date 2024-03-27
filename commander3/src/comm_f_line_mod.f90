module comm_F_line_mod
  use comm_F_int_mod
  implicit none

  private
  public comm_F_line

  type, extends (comm_F_int) :: comm_F_line
     logical(lgt) :: active
     real(dp)     :: f_precomp
     integer(i4b) :: ind
   contains
     ! Data procedures
     procedure :: eval       => evalIntF_line
     procedure :: eval_deriv => evalIntdF_line
     procedure :: update     => updateIntF_line
  end type comm_F_line

  interface comm_F_line
     procedure constructor_line
  end interface comm_F_line

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_line(comp, bp, active, scale, ind) result(c)
    implicit none
    class(comm_comp),                 intent(in), target :: comp
    class(comm_bp),                   intent(in), target :: bp
    logical(lgt),                     intent(in)         :: active
    real(dp),                         intent(in)         :: scale 
    integer(i4b),                     intent(in)         :: ind
    class(comm_F_line), pointer                          :: c

    allocate(c)
    c%comp   => comp
    c%bp     => bp
    c%active =  active

    if (c%active) then
       c%f_precomp = scale
       c%ind       = ind
    else
       c%f_precomp = 0.d0
    end if


  end function constructor_line

  
  ! Evaluate SED in brightness temperature normalized to reference frequency
  function evalIntF_line(self, theta)
    class(comm_F_line),                  intent(in) :: self
    real(dp),             dimension(1:), intent(in) :: theta
    real(dp)                                        :: evalIntF_line
    if (self%active) then
       if (self%ind > 0) then
          evalIntF_line = self%f_precomp * theta(self%ind)
       else
          evalIntF_line = self%f_precomp
       end if
    else 
       evalIntF_line = 0.d0
    end if
  end function evalIntF_line
  
  ! Evaluate partial derivative of SED in brightness temperature normalized to reference frequency
  function evalIntdF_line(self, theta, par)
    class(comm_F_line),               intent(in) :: self
    real(dp),          dimension(1:), intent(in) :: theta
    integer(i4b),                     intent(in) :: par
    real(dp)                                     :: evalIntdF_line
    if (self%active) then
       if (self%ind == par) then
          evalIntdF_line = self%f_precomp
       else
          evalIntdF_line = 0.d0
       end if
    else 
       evalIntdF_line = 0.d0
    end if
  end function evalIntdF_line

  ! Compute/update integration look-up tables
  subroutine updateIntF_line(self, f_precomp, pol)
    class(comm_F_line),   intent(inout)           :: self
    real(dp),             intent(in),    optional :: f_precomp
    integer(i4b),         intent(in),    optional :: pol

  end subroutine updateIntF_line


end module comm_F_line_mod
