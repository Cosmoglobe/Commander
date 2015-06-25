module comm_powlaw_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  implicit none

  private
  public comm_powlaw_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_powlaw_comp
   contains
     procedure :: S    => evalSED
  end type comm_powlaw_comp

  interface comm_powlaw_comp
     procedure constructor
  end interface comm_powlaw_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id
    class(comm_powlaw_comp), pointer   :: constructor

    ! General parameters
    allocate(constructor)
    call constructor%initDiffuse(cpar, id)

    ! Component specific parameters
    constructor%npar = 1
    allocate(constructor%theta_def(1), constructor%p_gauss(2,1), constructor%p_uni(2,1))
    constructor%theta_def(1) = cpar%cs_theta_def(id,1)
    constructor%p_uni(:,1)   = cpar%cs_p_uni(id,:,1)
    constructor%p_gauss(:,1) = cpar%cs_p_gauss(id,:,1)
  end function constructor

  ! Definition:
  !    SED  = (nu/nu_ref)**beta
  ! where 
  !    beta = theta(1)
  function evalSED(self, nu, band, theta)
    class(comm_powlaw_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    evalSED = (nu/self%nu_ref)**theta(1)

  end function evalSED
  
end module comm_powlaw_comp_mod
