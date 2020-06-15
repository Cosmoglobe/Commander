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
    allocate(constructor)


    constructor%npar         = 2
    allocate(constructor%poltype(constructor%npar))
    do i = 1, constructor%npar
       constructor%poltype(i)   = cpar%cs_poltype(i,id_abs)
    end do
    call constructor%initLmaxSpecind(cpar, id, id_abs)

    call constructor%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    allocate(constructor%theta_def(2), constructor%p_gauss(2,2), constructor%p_uni(2,2))
    allocate(constructor%indlabel(2))
    allocate(constructor%nu_min_ind(2), constructor%nu_max_ind(2))
    do i = 1, 2
       constructor%theta_def(i)  = cpar%cs_theta_def(i,id_abs)
       constructor%p_uni(:,i)    = cpar%cs_p_uni(id_abs,:,i)
       constructor%p_gauss(:,i)  = cpar%cs_p_gauss(id_abs,:,i)
       constructor%nu_min_ind(i) = cpar%cs_nu_min(id_abs,i)
       constructor%nu_max_ind(i) = cpar%cs_nu_max(id_abs,i)
    end do
    constructor%indlabel  = ['beta', 'T   ']

    ! Init alm 
    if (constructor%lmax_ind >= 0) call constructor%initSpecindProp(cpar, id, id_abs)

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
             constructor%F_int(k,i,j)%p => comm_F_int_2D(constructor, data(i)%bp(j)%p, k)
          end do
       end do
    end do

    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & constructor%nmaps, constructor%pol)

    allocate(constructor%theta(constructor%npar))
    do i = 1, constructor%npar
       if (trim(cpar%cs_input_ind(i,id_abs)) == 'default') then
          constructor%theta(i)%p => comm_map(info)
          constructor%theta(i)%p%map = constructor%theta_def(i)
       else
          ! Read map from FITS file, and convert to alms
          constructor%theta(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_input_ind(i,id_abs)))
       end if

       !convert spec. ind. pixel map to alms if lmax_ind >= 0
       if (constructor%lmax_ind >= 0) then
          ! if lmax >= 0 we can get alm values for the theta map
          call constructor%theta(i)%p%YtW_scalar
       end if
    end do

    call constructor%initPixregSampling(cpar, id, id_abs)

    ! Initialize mixing matrix
    call constructor%updateMixmat

  end function constructor

  ! Definition:
  !      x  = h*nu/(k_b*T)
  !    SED  = (nu/nu_ref)**(beta+1) * (exp(x_ref)-1)/(exp(x)-1)
  ! where 
  !    beta = theta(1)
  function evalSED(self, nu, band, pol, theta)
    implicit none
    class(comm_MBB_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    real(dp) :: x, x_ref, beta, T
    
    beta    = theta(1)
    T       = theta(2)
    x       = h*nu               / (k_b*T)
    x_ref   = h*self%nu_ref(pol) / (k_b*T)
    evalSED = (nu/self%nu_ref(pol))**(beta+1.d0) * (exp(x_ref)-1.d0)/(exp(x)-1.d0)

  end function evalSED
  
end module comm_MBB_comp_mod
