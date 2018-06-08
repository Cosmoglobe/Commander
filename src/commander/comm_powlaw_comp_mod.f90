module comm_powlaw_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_map_mod
  use comm_F_int_1D_mod
  use comm_data_mod
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
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_powlaw_comp), pointer   :: constructor

    integer(i4b) :: i
    type(comm_mapinfo), pointer :: info

    ! General parameters
    allocate(constructor)
    call constructor%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    constructor%npar         = 1
    allocate(constructor%theta_def(1), constructor%p_gauss(2,1), constructor%p_uni(2,1))
    allocate(constructor%poltype(1), constructor%indlabel(1))
    allocate(constructor%nu_min_ind(1), constructor%nu_max_ind(1))
    constructor%poltype(1)   = cpar%cs_poltype(1,id_abs)
    constructor%theta_def(1) = cpar%cs_theta_def(1,id_abs)
    constructor%p_uni(:,1)   = cpar%cs_p_uni(id_abs,:,1)
    constructor%p_gauss(:,1) = cpar%cs_p_gauss(id_abs,:,1)
    constructor%nu_min_ind(1) = cpar%cs_nu_min(id_abs,1)
    constructor%nu_max_ind(1) = cpar%cs_nu_max(id_abs,1)
    constructor%indlabel(1)  = 'beta'

    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & constructor%nmaps, constructor%pol)

    allocate(constructor%theta(1))
    if (trim(cpar%cs_input_ind(1,id_abs)) == 'default') then
       constructor%theta(1)%p => comm_map(info)
       constructor%theta(1)%p%map = constructor%theta_def(1)
    else
       ! Read map from FITS file, and convert to alms
       constructor%theta(1)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_input_ind(1,id_abs)))
    end if
    if (constructor%lmax_ind >= 0) call constructor%theta(1)%p%YtW_scalar

    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(numband))
    do i = 1, numband
       constructor%F_int(i)%p => comm_F_int_1D(constructor, data(i)%bp)
    end do

    ! Initialize mixing matrix
    call constructor%updateMixmat

    ! Set up smoothing scale information
    allocate(constructor%smooth_scale(constructor%npar))
    constructor%smooth_scale = cpar%cs_smooth_scale(id_abs,1:constructor%npar)

  end function constructor

  ! Definition:
  !    SED  = (nu/nu_ref)**beta
  ! where 
  !    beta = theta(1)
  function evalSED(self, nu, band, theta)
    implicit none
    class(comm_powlaw_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    evalSED = (nu/self%nu_ref)**theta(1)

  end function evalSED
  
end module comm_powlaw_comp_mod
