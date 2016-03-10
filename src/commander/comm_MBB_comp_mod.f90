module comm_MBB_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_map_mod
  use comm_F_int_2D_mod
  use comm_data_mod
  implicit none

  private
  public comm_MBB_comp

  !**************************************************
  !      Modified Black Body (MBB) component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_MBB_comp
   contains
     procedure :: S    => evalSED
  end type comm_MBB_comp

  interface comm_MBB_comp
     procedure constructor
  end interface comm_MBB_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_MBB_comp), pointer   :: constructor

    integer(i4b) :: i
    type(comm_mapinfo), pointer :: info

    ! General parameters
    allocate(constructor)
    call constructor%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    constructor%npar         = 2
    allocate(constructor%theta_def(2), constructor%p_gauss(2,2), constructor%p_uni(2,2))
    allocate(constructor%poltype(2), constructor%indlabel(2))
    do i = 1, 2
       constructor%theta_def(i) = cpar%cs_theta_def(i,id_abs)
       constructor%p_uni(:,i)   = cpar%cs_p_uni(id_abs,:,i)
       constructor%p_gauss(:,i) = cpar%cs_p_gauss(id_abs,:,i)
    end do
    constructor%indlabel  = ['beta', 'T']

    ! Initialize spectral parameter maps
    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & constructor%nmaps, constructor%pol)

    allocate(constructor%theta(2))
    do i = 1, 2
       if (trim(cpar%cs_input_ind(i,id_abs)) == 'default') then
          constructor%theta(i)%p => comm_map(info)
          constructor%theta(i)%p%map = constructor%theta_def(i)
       else
          ! Read map from FITS file, and convert to alms
          constructor%theta(i)%p => comm_map(info, cpar%cs_input_ind(i,id_abs))
       end if
       call constructor%theta(i)%p%YtW
    end do

    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(numband))
    do i = 1, numband
       constructor%F_int(i)%p => comm_F_int_2D(constructor, data(i)%bp)
    end do

    ! Initialize mixing matrix
    call constructor%updateMixmat

  end function constructor

  ! Definition:
  !      x  = h*nu/(k_b*T)
  !    SED  = (nu/nu_ref)**(beta+1) * (exp(x_ref)-1)/(exp(x)-1)
  ! where 
  !    beta = theta(1)
  function evalSED(self, nu, band, theta)
    implicit none
    class(comm_MBB_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    real(dp) :: x, x_ref, beta, T
    
    beta    = theta(1)
    T       = theta(2)
    x       = h*nu          / (k_b*T)
    x_ref   = h*self%nu_ref / (k_b*T)
    evalSED = (nu/self%nu_ref)**(beta+1.d0) * (exp(x_ref)-1.d0)/(exp(x)-1.d0)

  end function evalSED
  
end module comm_MBB_comp_mod
