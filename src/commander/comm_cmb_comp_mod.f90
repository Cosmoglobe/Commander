module comm_cmb_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_F_int_0D_mod
  use comm_data_mod
  use comm_bp_utils
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
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_cmb_comp), pointer   :: constructor

    integer(i4b) :: i
    real(dp)     :: f
    
    ! General parameters
    allocate(constructor)
    call constructor%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    constructor%npar         = 0

    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(numband))
    do i = 1, numband
       f = comp_a2t(constructor%nu_ref) / data(i)%bp%a2t * data(i)%RJ2data()
       constructor%F_int(i)%p => comm_F_int_0D(constructor, data(i)%bp, f_precomp=f)
    end do
    
    ! Initialize mixing matrix
    call constructor%updateMixmat

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
