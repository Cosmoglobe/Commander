module comm_spindust_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_spindust_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_spindust_comp
     real(dp)          :: nu_p0, nu_min_SED, nu_max_SED
     type(spline_type) :: SED_spline
   contains
     procedure :: S    => evalSED_spindust
  end type comm_spindust_comp

  interface comm_spindust_comp
     procedure constructor_spindust
  end interface comm_spindust_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_spindust(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_spindust_comp), pointer   :: c

    integer(i4b) :: ind(1)
    real(dp), allocatable, dimension(:,:) :: SED
    integer(i4b) :: i, j, k, l, m, n, p, ierr
    type(comm_mapinfo), pointer :: info => null()
    real(dp)           :: par_dp
    integer(i4b), allocatable, dimension(:) :: sum_pix
    real(dp),    allocatable, dimension(:) :: sum_theta, sum_proplen, sum_nprop
    character(len=512) :: temptxt, partxt
    integer(i4b) :: smooth_scale, p_min, p_max
    class(comm_mapinfo), pointer :: info2 => null()
    class(comm_map),     pointer :: tp => null() 
    class(comm_map),     pointer :: tp_smooth => null() 

    ! General parameters
    allocate(c)

  end function constructor_spindust

  ! Definition:
  !    SED  = (nu/nu_ref)**beta
  ! where 
  !    beta = theta(1)
  function evalSED_spindust(self, nu, band, pol, theta)
    implicit none
    class(comm_spindust_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band    
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_spindust

    real(dp) :: scale, nu_p

    nu_p    = theta(1)
    !alpha   = theta(2)
    scale   = self%nu_p0 / (nu_p*1.d9) ! nu_p is in GHz

    if (scale*nu < self%nu_min_SED .or. scale*nu > self%nu_max_SED) then
       evalSED_spindust = 0.d0
    else
       evalSED_spindust = exp(splint(self%SED_spline, log(scale*nu))) / &
               & exp(splint(self%SED_spline, log(scale*self%nu_ref(pol)))) * (self%nu_ref(pol)/nu)**2
    !           & exp(splint(self%SED_spline, log(scale*self%nu_ref))) * (self%nu_ref/nu)**(2.d0-alpha)
    end if

  end function evalSED_spindust
  
end module comm_spindust_comp_mod
