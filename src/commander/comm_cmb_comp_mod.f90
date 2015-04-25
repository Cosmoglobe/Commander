module comm_cmb_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  implicit none

  private
  public comm_cmb_comp

  !**************************************************
  !           CMB component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_cmb_comp
   contains
     procedure :: S    => evalSED
  end type comm_cmb_comp

  interface comm_cmb_comp
     procedure constructor
  end interface comm_cmb_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id
    class(comm_cmb_comp), pointer   :: constructor

    ! General parameters
    allocate(constructor)
    call constructor%initDiffuse(cpar, id)

    ! Component specific parameters
    constructor%npar = 0
    
  end function constructor

  ! Definition:
  !    SED  = conversion between thermodynamic and brightness temperature = 1/a2t
  function evalSED(self, nu, band, theta)
    class(comm_cmb_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    real(dp) :: x
    x       = h*nu / (k_B*T_CMB)
    evalSED = (x**2 * exp(x)) / (exp(x)-1.d0)**2
    
  end function evalSED
  
end module comm_cmb_comp_mod
