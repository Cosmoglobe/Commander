module comm_spindust_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_map_mod
  use comm_F_int_1D_mod
  use comm_data_mod
  use spline_1D_mod
  implicit none

  private
  public comm_spindust_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_spindust_comp
     real(dp)          :: nu_p0, nu_min, nu_max
     type(spline_type) :: SED_spline
   contains
     procedure :: S    => evalSED
  end type comm_spindust_comp

  interface comm_spindust_comp
     procedure constructor
  end interface comm_spindust_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_spindust_comp), pointer   :: constructor

    integer(i4b) :: i, ind(1)
    real(dp), allocatable, dimension(:,:) :: SED
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
    constructor%indlabel(1)  = 'nu_p'
    constructor%nu_min_ind(1) = cpar%cs_nu_min(id_abs,1)
    constructor%nu_max_ind(1) = cpar%cs_nu_max(id_abs,1)

    ! Component specific parameters for 2 parameter model
    !constructor%npar         = 2
    !allocate(constructor%theta_def(2), constructor%p_gauss(2,2), constructor%p_uni(2,2))
    !allocate(constructor%poltype(2), constructor%indlabel(2))
    !allocate(constructor%nu_min_ind(2), constructor%nu_max_ind(2))
    !do i = 1, 2
    !   constructor%poltype(i)   = cpar%cs_poltype(i,id_abs)
    !   constructor%theta_def(i) = cpar%cs_theta_def(i,id_abs)
    !   constructor%p_uni(:,i)   = cpar%cs_p_uni(id_abs,:,i)
    !   constructor%p_gauss(:,i) = cpar%cs_p_gauss(id_abs,:,i)
    !   constructor%nu_min_ind(i) = cpar%cs_nu_min(id_abs,i)
    !   constructor%nu_max_ind(i) = cpar%cs_nu_max(id_abs,i)
    !end do
    !constructor%indlabel  = ['nu_p','alpha']
    
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

    ! Initialize spectral index map for 2 parameter model
    !allocate(constructor%theta(constructor%npar))
    !do i = 1, constructor%npar
    !   if (trim(cpar%cs_input_ind(i,id_abs)) == 'default') then
    !      constructor%theta(i)%p => comm_map(info)
    !      constructor%theta(i)%p%map = constructor%theta_def(i)
    !   else
    !      ! Read map from FITS file, and convert to alms
    !      constructor%theta(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_input_ind(i,id_$
    !   end if
    !   if (constructor%lmax_ind >= 0) call constructor%theta(1)%p%YtW_scalar
    !end do

    ! Initialize spectral template !CHANGE??
    call read_spectrum(trim(cpar%datadir)//'/'//trim(cpar%cs_SED_template(1,id_abs)), SED)
    ind                = maxloc(SED(:,2))
    constructor%nu_p0  = SED(ind(1),1)
    constructor%nu_min = minval(SED(:,1))
    constructor%nu_max = maxval(SED(:,1))
    SED                = log(SED)
    call spline(constructor%SED_spline, SED(:,1), SED(:,2))
    deallocate(SED)

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
    class(comm_spindust_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    real(dp) :: scale, nu_p

    nu_p    = theta(1)
    !alpha   = theta(2)
    scale   = self%nu_p0 / (nu_p*1.d9) ! nu_p is in GHz

    if (scale*nu < self%nu_min .or. scale*nu > self%nu_max) then
       evalSED = 0.d0
    else
       evalSED = exp(splint(self%SED_spline, log(scale*nu))) / &
               & exp(splint(self%SED_spline, log(scale*self%nu_ref))) * (self%nu_ref/nu)**2
    !           & exp(splint(self%SED_spline, log(scale*self%nu_ref))) * (self%nu_ref/nu)**(2.d0-alpha)
    end if

  end function evalSED
  
end module comm_spindust_comp_mod
