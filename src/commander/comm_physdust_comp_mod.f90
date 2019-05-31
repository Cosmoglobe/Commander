module comm_physdust_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_map_mod
  use comm_F_int_1D_mod
  use comm_data_mod
  use spline_2D_mod
  implicit none

  private
  public comm_physdust_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_physdust_comp
     integer(i4b) :: num_nu, num_u, num_comp, num_u_int
     real(dp)     :: log_umax, gamma, alpha
     real(dp)     :: u_arr_min, u_arr_max, du_arr
     real(dp), allocatable, dimension(:)         :: log_dust_wav, dust_u, extcrv, extpol, scat, amps
     real(dp), allocatable, dimension(:,:,:,:,:) :: coeff_arr
   contains
     procedure :: S    => evalSED
  end type comm_physdust_comp

  interface comm_physdust_comp
     procedure constructor
  end interface comm_physdust_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_physdust_comp), pointer   :: constructor

    integer(i4b) :: i, j, k, num_nu, unit
    real(dp)     :: dummy
    type(comm_mapinfo), pointer :: info

    real(dp), allocatable, dimension(:,:,:) :: comp_mat


    ! General parameters
    allocate(constructor)
    call constructor%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    constructor%npar         = 1
    allocate(constructor%theta_def(1), constructor%p_gauss(2,1), constructor%p_uni(2,1))
    allocate(constructor%poltype(1), constructor%indlabel(1))
    constructor%poltype(1)   = cpar%cs_poltype(1,id_abs)
    constructor%theta_def(1) = cpar%cs_theta_def(1,id_abs)
    constructor%p_uni(:,1)   = cpar%cs_p_uni(id_abs,:,1)
    constructor%p_gauss(:,1) = cpar%cs_p_gauss(id_abs,:,1)
    constructor%indlabel(1)  = 'Umin'

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
    if (constructor%lmax_ind >= 0) call constructor%theta(1)%p%YtW

    ! Initialize SED templates
    unit                  = getlun()
    num_nu                = 1000    ! Number of frequencies in dust files
    constructor%num_nu    = num_nu
    constructor%num_u     = 11      ! Number of radiation field strengths
    constructor%u_arr_min = -0.5d0
    constructor%u_arr_max =  0.5d0
    constructor%du_arr    = (constructor%u_arr_max - constructor%u_arr_min)/(constructor%num_u - 1.d0)
    constructor%num_comp  = 4 ! Number of dust components

    allocate(constructor%coeff_arr(4,4,num_nu,constructor%num_u,constructor%num_comp))
    allocate(constructor%log_dust_wav(num_nu), constructor%dust_u(constructor%num_u))
    allocate(constructor%extcrv(num_nu), constructor%extpol(num_nu), constructor%scat(num_nu))
    allocate(constructor%amps(constructor%num_comp))

    constructor%log_umax  = cpar%cs_auxpar(1,id_abs)
    constructor%gamma     = cpar%cs_auxpar(2,id_abs)
    constructor%alpha     = cpar%cs_auxpar(3,id_abs)
    constructor%amps      = cpar%cs_auxpar(4:3+constructor%num_comp,id_abs)

    ! Make array of log U
    do i = 1, constructor%num_u
       constructor%dust_u(i) = -0.5d0 + constructor%du_arr*(i-1)
    enddo

    allocate(comp_mat(num_nu,constructor%num_u,constructor%num_comp))    
    do k = 1, constructor%num_comp
       open(unit,file=trim(cpar%datadir)//'/'//trim(cpar%cs_SED_template(k,id_abs)),recl=4096)
       do i = 1, num_nu
          read(unit,fmt='(1P26E11.3)') constructor%log_dust_wav(i), constructor%extcrv(i),&
               & constructor%extpol(i), constructor%scat(i), (comp_mat(i,j,k),j=1,constructor%num_u), &
                (dummy,j=1,constructor%num_u)
          constructor%log_dust_wav(i) = log(constructor%log_dust_wav(i))
       enddo
       call splie2_full_precomp(constructor%log_dust_wav,constructor%dust_u,&
            & log(comp_mat(:,:,k)), constructor%coeff_arr(:,:,:,:,k))
        close(59)
     enddo
     deallocate(comp_mat)


    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(numband))
    do i = 1, numband
       constructor%F_int(i)%p => comm_F_int_1D(constructor, data(i)%bp(0)%p)
    end do

    ! Initialize mixing matrix
    call constructor%updateMixmat

  end function constructor

  ! Definition:
  !    SED  = (nu/nu_ref)**beta
  ! where 
  !    beta = theta(1)
  function evalSED(self, nu, band, theta)
    implicit none
    class(comm_physdust_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    integer(i4b) :: i, j, num_u_int
    real(dp) :: SED_physdust, SED_norm, wav_in, wav_ref, c
    real(dp) :: log_umin, umin, log_umax, umax, uval, du, fdu

    if (nu < 2d9) then
       evalSED = 0.d0
       return
    end if
     
    num_u_int = 100 ! How many dU steps to perform
    c         = 2.99792458d8
     
    log_umin     = theta(1)
    log_umax     = self%log_umax
    SED_physdust = 0.d0
    SED_norm     = 0.d0
    wav_in       = c/nu          * 1d6 ! In um
    wav_ref      = c/self%nu_ref * 1d6 ! In um
    umin         = 10**log_umin
    umax         = 10**log_umax
     do i = 1, self%num_comp
        ! Delta function component
        SED_physdust = SED_physdust + (1.d0-self%gamma)* &
             & (self%amps(i)*exp(splin2_full_precomp_irreg(self%log_dust_wav,self%dust_u,&
             & self%coeff_arr(:,:,:,:,i), log(wav_in), log_umin)))
        SED_norm = SED_norm + (1.d0-self%gamma)*&
             & (self%amps(i)*exp(splin2_full_precomp_irreg(self%log_dust_wav,self%dust_u,self%coeff_arr(:,:,:,:,i), &
             log(wav_ref),log_umin)))

        ! Distribution of U values
        ! See, e.g., Aniano et al 2012, Eqs. 9 & 10
        if (self%gamma /= 0.d0) then
           du = umin*((umax/umin)**(1.d0/(num_u_int-1.d0)) - 1d0) ! Formally dU/U
           do j = 1, num_u_int ! integrate over U distribution
              uval = umin*(umax/umin)**((j-1.d0)/(num_u_int-1.d0))
              if (self%alpha /= 1.d0) then
                 fdu = uval**(1.d0-self%alpha)*du*self%gamma*&
                      & (self%alpha-1.d0)/(umin**(1.d0-self%alpha) - umax**(1.d0-self%alpha))
              else
                 ! Special case of alpha = 1 to ensure integrability
                 fdu = du*self%gamma/log(umax/umin)
              endif
              SED_physdust = SED_physdust + self%amps(i)* &
                   & exp(splin2_full_precomp_irreg(self%log_dust_wav,self%dust_u,self%coeff_arr(:,:,:,:,i), &
                   & log(wav_in),log10(uval)))*fdu
              SED_norm = SED_norm + self%amps(i)* &
                   & exp(splin2_full_precomp_irreg(self%log_dust_wav,self%dust_u,self%coeff_arr(:,:,:,:,i), &
                   log(wav_ref),log10(uval)))*fdu
           enddo
        endif
     enddo
     evalSED = (SED_physdust / SED_norm) * (self%nu_ref/nu)**3 ! Normalize to reference in T units

  end function evalSED
  
end module comm_physdust_comp_mod
