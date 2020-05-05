module comm_freefree_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_map_mod
  use comm_F_int_1D_mod
  use comm_data_mod
  implicit none

  private
  public comm_freefree_comp

  !**************************************************
  !      Free-free component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_freefree_comp
   contains
     procedure :: S    => evalSED
  end type comm_freefree_comp

  interface comm_freefree_comp
     procedure constructor
  end interface comm_freefree_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_freefree_comp), pointer   :: constructor

    integer(i4b) :: i, j, k
    type(comm_mapinfo), pointer :: info => null()

    ! General parameters
    allocate(constructor)
    call constructor%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    
    constructor%npar         = 1
    allocate(constructor%theta_def(1), constructor%p_gauss(2,1), constructor%p_uni(2,1))
    allocate(constructor%poltype(1), constructor%indlabel(1))
    allocate(constructor%nu_min_ind(1), constructor%nu_max_ind(1))
    do i = 1, 1
       constructor%poltype(i)   = cpar%cs_poltype(i,id_abs)
       constructor%theta_def(i) = cpar%cs_theta_def(i,id_abs)
       constructor%p_uni(:,i)   = cpar%cs_p_uni(id_abs,:,i)
       constructor%p_gauss(:,i) = cpar%cs_p_gauss(id_abs,:,i)
       constructor%nu_min_ind(i) = cpar%cs_nu_min(id_abs,i)
       constructor%nu_max_ind(i) = cpar%cs_nu_max(id_abs,i)
    end do
    constructor%indlabel  = ['Te']

    ! Init alm 
    if (constructor%lmax_ind >= 0) call constructor%initSpecindProp(cpar, id, id_abs)

    !constructor%npar         = 1
    !allocate(constructor%theta_def(1), constructor%p_gauss(1,1), constructor%p_uni(1,1))
    !allocate(constructor%poltype(1), constructor%indlabel(1))
    !allocate(constructor%nu_min_ind(1), constructor%nu_max_ind(1))
    !i = 1
    !constructor%poltype(i)   = cpar%cs_poltype(i,id_abs)
    !constructor%theta_def(i) = cpar%cs_theta_def(i,id_abs)
    !constructor%p_uni(:,i)   = cpar%cs_p_uni(id_abs,:,i)
    !constructor%p_gauss(:,i) = cpar%cs_p_gauss(id_abs,:,i)
    !constructor%nu_min_ind(i) = cpar%cs_nu_min(id_abs,i)
    !constructor%nu_max_ind(i) = cpar%cs_nu_max(id_abs,i)
    
    !constructor%indlabel  = ['Te']

    ! Initialize spectral parameter maps
    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & constructor%nmaps, constructor%pol)

    allocate(constructor%theta(constructor%npar))
    do i = 1, constructor%npar
       if (trim(cpar%cs_input_ind(i,id_abs)) == 'default') then
          constructor%theta(i)%p => comm_map(info)
          constructor%theta(i)%p%map = constructor%theta_def(i)
       else
          ! Read map from FITS file, and convert to alms
          constructor%theta(i)%p => comm_map(info, trim(cpar%datadir)//'/'//trim(cpar%cs_input_ind(i,id_abs)))
       end if
       if (constructor%lmax_ind >= 0) call constructor%theta(i)%p%YtW_scalar
    end do

    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(3,numband,0:constructor%ndet))
    do k = 1, 3
       do i = 1, numband
          do j = 0, data(i)%ndet
             if (k > 1) then
                if (constructor%nu_ref(k) == constructor%nu_ref(k-1)) then
                   constructor%F_int(k,i,j)%p => constructor%F_int(k-1,i,j)%p
                   cycle
                end if
             end if
             constructor%F_int(k,i,j)%p => comm_F_int_1D(constructor, data(i)%bp(j)%p, k)
          end do
       end do
    end do

    ! Initialize mixing matrix
    call constructor%updateMixmat

    ! Set up smoothing scale information
    allocate(constructor%smooth_scale(constructor%npar))
    constructor%smooth_scale = cpar%cs_smooth_scale(id_abs,1:constructor%npar)

  end function constructor

  ! Definition:
  !      x  = h*nu/(k_b*T)
  !    SED  = (nu/nu_ref)**(beta+1) * (exp(x_ref)-1)/(exp(x)-1)
  ! where 
  !    beta = theta(1)
  function evalSED(self, nu, band, pol, theta)
    implicit none
    class(comm_freefree_comp),    intent(in)      :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED
    real(dp)     :: S, S_ref, EM, T_e
    real(dp)     :: g, g_ref, Z_i, tau, tau_ref, EM1, Te

!!$    EM      = theta(1)
!!$    !EM1 = 1.d0 
!!$    Te      = theta(2)
!!$    Z_i     = 1.d0
!!$    g       = log(exp(5.960d0 - sqrt(3.d0)/pi * log(Z_i * nu/1.d9          * (Te/1.d4)**(-1.5d0))) + 2.71828d0)
!!$    !g_ref   = log(exp(5.960d0 - sqrt(3.d0)/pi * log(Z_i * self%nu_ref/1.d9 * (Te/1.d4)**(-1.5d0))) + 2.71828d0)
!!$    tau     = 5.468d-2 * Te**(-1.5d0) * (nu/1.d9)**(-2)          * EM * g
!!$    !tau_ref = 5.468d-2 * Te**(-1.5d0) * (self%nu_ref/1.d9)**(-2) * EM * g_ref
!!$
!!$    evalSED = 1.d6 * Te * (1.d0 - exp(-tau)) !/ (1.d0 - exp(-tau_ref)) 
!!$    
!!$    return
!!$    !write(*,*) "1:", evalSED


    !EM    = theta(1) ! Not used
    T_e   = theta(1)
    S     = log(exp(5.960d0 - sqrt(3.d0)/pi * log(1.d0 * nu    /1.d9 * (T_e/1.d4)**(-1.5d0))) + 2.71828d0)
    S_ref = log(exp(5.960d0 - sqrt(3.d0)/pi * log(1.d0 * self%nu_ref(pol)/1.d9 * (T_e/1.d4)**(-1.5d0))) + 2.71828d0)
    !evalSED = S/S_ref * exp(-h*(nu-self%nu_ref(pol))/k_b/T_e) * (nu/self%nu_ref(pol))**(-2)
    evalSED = S/S_ref * (nu/self%nu_ref(pol))**-2
    !write(*,*) "2",evalSED
    
  end function evalSED
  
end module comm_freefree_comp_mod
