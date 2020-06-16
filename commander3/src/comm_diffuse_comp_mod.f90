module comm_diffuse_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_map_mod
  use comm_data_mod
  use comm_F_int_mod
  use comm_Cl_mod
  use math_tools
  use comm_cr_utils
  use comm_cr_precond_mod
  use comm_hdf_mod
  use InvSamp_mod
  use powell_mod
  use comm_beam_mod
  implicit none

  private
  public comm_diffuse_comp, add_to_npre, updateDiffPrecond, initDiffPrecond, applyDiffPrecond, &
       & res_smooth, rms_smooth, print_precond_mat, nullify_monopole_amp, recompute_diffuse_precond
  
  !**************************************************
  !            Diffuse component class
  !**************************************************
  type, abstract, extends (comm_comp) :: comm_diffuse_comp
     character(len=512) :: cltype
     integer(i4b)       :: nside, nx, x0, ndet
     logical(lgt)       :: pol, output_mixmat, output_EB, apply_jeffreys, almsamp_pixreg
     integer(i4b)       :: lmax_amp, lmax_ind, lmax_prior, lpiv, l_apod, lmax_pre_lowl
     integer(i4b)       :: lmax_def, nside_def, ndef, nalm_tot, sample_first_niter

     real(dp),     allocatable, dimension(:,:)   :: sigma_priors, steplen, chisq_min
     real(dp),     allocatable, dimension(:,:,:,:)   :: L
     integer(i4b), allocatable, dimension(:,:)   :: corrlen     
     logical(lgt),    dimension(:), allocatable :: L_read

     real(dp)           :: cg_scale, latmask, fwhm_def, test
     real(dp),           allocatable, dimension(:,:)   :: cls
     real(dp),           allocatable, dimension(:,:,:) :: F_mean
     character(len=128), allocatable, dimension(:,:)   :: pol_lnLtype     ! {'chisq', 'ridge', 'marginal'}
     integer(i4b),       allocatable, dimension(:,:)   :: lmax_ind_pol    ! lmax per poltype per spec. ind  
     integer(i4b),       allocatable, dimension(:,:)   :: lmax_ind_mix    ! equal to lmax_ind_pol, but 0 where lmax=-1 and fullsky pixreg
     integer(i4b),       allocatable, dimension(:,:)   :: pol_pixreg_type ! {1=fullsky, 2=single_pix, 3=pixel_regions}
     integer(i4b),       allocatable, dimension(:,:)   :: nprop_uni       ! limits on nprop
     integer(i4b),       allocatable, dimension(:,:)   :: npixreg         ! number of pixel regions
     integer(i4b),       allocatable, dimension(:,:,:) :: ind_pixreg_arr  ! number of pixel regions
     real(dp),           allocatable, dimension(:,:,:) :: theta_pixreg    ! thetas for pixregs, per poltype, per ind.
     real(dp),           allocatable, dimension(:,:,:) :: proplen_pixreg  ! proposal length for pixregs
     integer(i4b),       allocatable, dimension(:,:,:) :: nprop_pixreg    ! number of proposals for pixregs
     integer(i4b),       allocatable, dimension(:,:,:) :: npix_pixreg     ! number of pixels per pixel region
     logical(lgt),       allocatable, dimension(:,:)   :: pol_sample_nprop   ! sample the corr. length in first iteration
     logical(lgt),       allocatable, dimension(:,:)   :: pol_sample_proplen ! sample the prop. length in first iteration
     logical(lgt),       allocatable, dimension(:,:)   :: first_ind_sample
     class(map_ptr),     allocatable, dimension(:)     :: pol_ind_mask
     class(map_ptr),     allocatable, dimension(:)     :: pol_proplen
     class(map_ptr),     allocatable, dimension(:)     :: ind_pixreg_map   !map with defined pixelregions
     class(map_ptr),     allocatable, dimension(:)     :: pol_nprop
     class(comm_B_bl_ptr), allocatable, dimension(:)   :: B_pp_fr

     character(len=512) :: mono_prior_type
     class(comm_map),               pointer     :: mono_prior_map => null()
     class(comm_map),               pointer     :: mask => null()
     class(comm_map),               pointer     :: procmask => null()
     class(map_ptr),  dimension(:), allocatable :: indmask
     class(comm_map),               pointer     :: defmask => null()
     class(comm_map),               pointer     :: priormask => null()
     class(comm_map),               pointer     :: x => null()           ! Spatial parameters
     class(comm_map),               pointer     :: x_smooth => null()    ! Spatial parameters
     class(comm_map),               pointer     :: mu => null()          ! Spatial prior mean
     class(comm_B),                 pointer     :: B_out => null()       ! Output beam
     class(comm_Cl),                pointer     :: Cl => null()          ! Power spectrum
     class(map_ptr),  dimension(:), allocatable :: theta        ! Spectral parameters
     class(map_ptr),  dimension(:), allocatable :: theta_smooth ! Spectral parameters
     type(map_ptr),   dimension(:,:), allocatable :: F          ! Mixing matrix
     type(map_ptr),   dimension(:,:), allocatable :: dF         ! Derivative of mixing matrix
     real(dp),        dimension(:,:), allocatable :: invM_lowl  ! (0:nalm-1,0:nalm-1)
     real(dp),        dimension(:,:), allocatable :: Z_def      ! (0:nalm-1,ndef)
     real(dp),        dimension(:,:), allocatable :: invM_def   ! (0:nalm-1,0:nalm-1)
     logical(lgt),    dimension(:,:), allocatable :: F_null     ! Don't allocate space for null mixmat's
     type(F_int_ptr), dimension(:,:,:), allocatable :: F_int        ! SED integrator
   contains
     procedure :: initDiffuse
     procedure :: initPixregSampling
     procedure :: initSpecindProp
     procedure :: initLmaxSpecind
     procedure :: updateMixmat  => updateDiffuseMixmat
     procedure :: update_F_int  => updateDiffuseFInt
!!$     procedure :: dumpHDF    => dumpDiffuseToHDF
     procedure :: getBand       => evalDiffuseBand
     procedure :: projectBand   => projectDiffuseBand
     procedure :: dumpFITS      => dumpDiffuseToFITS
     procedure :: initHDF       => initDiffuseHDF
     procedure :: sampleSpecInd => sampleDiffuseSpecInd
     procedure :: sampleDiffuseSpecIndFullsky
     procedure :: sampleDiffuseSpecIndSinglePix
     procedure :: sampleDiffuseSpecIndPixReg
     procedure :: updateLowlPrecond
     procedure :: applyLowlPrecond
     procedure :: updateDeflatePrecond
     procedure :: applyDeflatePrecond
     procedure :: applyMonoDipolePrior
  end type comm_diffuse_comp

  type diff_ptr
     class(comm_diffuse_comp), pointer :: p => null()
  end type diff_ptr
  
  integer(i4b) :: npre      =  0
  integer(i4b) :: lmax_pre  = -1
  integer(i4b) :: nside_pre = 1000000
  integer(i4b) :: nmaps_pre = -1
  logical(lgt) :: recompute_diffuse_precond = .true.
  logical(lgt) :: output_cg_eigenvals
  logical(lgt), private :: only_pol
  character(len=512) :: outdir, precond_type
  integer(i4b),        allocatable, dimension(:) :: ind_pre
  class(comm_mapinfo), pointer                   :: info_pre => null()
  class(diff_ptr),     allocatable, dimension(:) :: diffComps

  character(len=24), private :: operation

  ! Variables for non-linear search
  class(comm_diffuse_comp), pointer,       private :: c_lnL => null()
  integer(i4b),                            private :: k_lnL, p_lnL, id_lnL
  real(dp), allocatable, dimension(:),     private :: a_lnL
  real(dp), allocatable, dimension(:),     private :: theta_lnL        
  real(dp), allocatable, dimension(:,:),   private :: buffer_lnL        
  logical(lgt),                            private :: apply_mixmat = .true.
  type(map_ptr),        allocatable, dimension(:) :: res_smooth
  type(comm_N_rms_ptr), allocatable, dimension(:) :: rms_smooth
  
contains

  subroutine initDiffuse(self, cpar, id, id_abs)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    character(len=512) :: filename
    integer(i4b) :: i, j, l, m, ntot, nloc, p
    real(dp) :: fwhm_prior, sigma_prior
    logical(lgt) :: exist
    type(comm_mapinfo), pointer :: info => null(), info_def => null(), info_ud
    class(comm_map), pointer :: indmask, mask_ud
    
    call self%initComp(cpar, id, id_abs)

    call update_status(status, "init_diffuse_start")

    ! Initialize variables specific to diffuse source type
    self%pol           = cpar%cs_polarization(id_abs)
    self%nside         = cpar%cs_nside(id_abs)
    self%lmax_amp      = cpar%cs_lmax_amp(id_abs)
    self%lmax_prior    = cpar%cs_lmax_amp_prior(id_abs)
    self%l_apod        = cpar%cs_l_apod(id_abs)

    if(self%npar == 0) then
       self%lmax_ind = 0 !default
       allocate(self%lmax_ind_mix(3,1))
       self%lmax_ind_mix = 0
    end if

    self%cltype        = cpar%cs_cltype(id_abs)
    self%cg_scale      = cpar%cs_cg_scale(id_abs)
    self%nmaps         = 1; if (self%pol) self%nmaps = 3
    self%output_mixmat = cpar%output_mixmat
    self%latmask       = cpar%cs_latmask(id_abs)
    self%apply_jeffreys = .false.
    self%sample_first_niter = cpar%cs_local_burn_in

    only_pol           = cpar%only_pol
    output_cg_eigenvals = cpar%output_cg_eigenvals
    outdir              = cpar%outdir
    precond_type        = cpar%cg_precond
    info               => comm_mapinfo(cpar%comm_chain, self%nside, self%lmax_amp, self%nmaps, self%pol)
    if (self%pol) then
       self%output_EB     = cpar%cs_output_EB(id_abs)
    else
       self%output_EB     = .false.
    end if
    operation          = cpar%operation


    ! Diffuse preconditioner variables
    npre      = npre+1
    lmax_pre  = max(lmax_pre,  self%lmax_amp)
    nside_pre = min(nside_pre, self%nside)
    nmaps_pre = max(nmaps_pre, self%nmaps)

    ! CMB low-l preconditioner
    if (trim(self%type) == 'cmb' .and. trim(precond_type) == 'diagonal' &
         & .and. self%lmax_pre_lowl > -1) then
       self%lmax_pre_lowl = cpar%cg_lmax_precond
       ntot               = (self%lmax_pre_lowl+1)**2
       nloc               = count(info%lm(1,:) <= self%lmax_pre_lowl)
       if (nloc > 0) allocate(self%invM_lowl(0:ntot-1,0:nloc-1))
    else
       self%lmax_pre_lowl = -1
    end if

    ! Initialize amplitude map
    if (trim(cpar%cs_input_amp(id_abs)) == 'zero' .or. trim(cpar%cs_input_amp(id_abs)) == 'none') then
       self%x => comm_map(info)
    else
       ! Read map from FITS file, and convert to alms
       self%x => comm_map(info, trim(cpar%datadir)//'/'//trim(cpar%cs_input_amp(id_abs)))
       do i = 1, self%x%info%nmaps
          self%x%map(:,i) = self%x%map(:,i) / (self%RJ2unit_(i)*self%cg_scale)
       end do
       call self%x%YtW
    end if
    self%ncr = size(self%x%alm)

    ! Initialize output beam
    self%B_out => comm_B_bl(cpar, self%x%info, 0, 0, fwhm=cpar%cs_fwhm(id_abs), nside=self%nside,&
         & init_realspace=.false.)

!!$    if (trim(cpar%cs_input_amp(id_abs)) /= 'zero' .and. trim(cpar%cs_input_amp(id_abs)) /= 'none') then
!!$       call self%B_out%deconv(.false., self%x)
!!$    end if

    ! Read component mask
    if (trim(cpar%cs_mask(id_abs)) /= 'fullsky' .and. self%latmask < 0.d0) then
       self%mask => comm_map(self%x%info, trim(cpar%datadir)//'/'//trim(cpar%cs_mask(id_abs)), &
            & udgrade=.true.)
    end if

    ! Read processing mask
    if (trim(cpar%ds_procmask) /= 'none') then
       self%procmask => comm_map(self%x%info, trim(cpar%datadir)//'/'//trim(cpar%ds_procmask), &
            & udgrade=.true.)
    end if

    ! Read index sampling mask; downgrade to each channel; re-use for equal Nsides; skip for dense matrices
    if (trim(cpar%cs_indmask(id_abs)) /= 'fullsky') then
       indmask => comm_map(self%x%info, trim(cpar%datadir)//'/'//trim(cpar%cs_indmask(id_abs)), &
            & udgrade=.true.)
       allocate(self%indmask(numband))
       do i = 1, numband
          if (trim(data(i)%N%type) == 'QUmat') cycle
          exist = .false.
          do j = 1, i-1
             if (data(i)%info%nside == data(j)%info%nside .and. &
                  & data(i)%info%nmaps == data(j)%info%nmaps) then
                self%indmask(i)%p => self%indmask(j)%p
                exist = .true.
                exit
             end if
          end do
          if (.not. exist) then
             info_ud => comm_mapinfo(cpar%comm_chain, indmask%info%nside, 0, data(i)%info%nmaps, data(i)%info%nmaps==3)
             mask_ud => comm_map(info_ud)
             do j = 1, min(data(i)%info%nmaps, indmask%info%nmaps)
                mask_ud%map(:,j) = indmask%map(:,j)
             end do
             self%indmask(i)%p => comm_map(data(i)%info)
             call mask_ud%udgrade(self%indmask(i)%p)
             where (self%indmask(i)%p%map > 0.5d0)
                self%indmask(i)%p%map = 1.d0
             elsewhere
                self%indmask(i)%p%map = 0.d0
             end where
             call mask_ud%dealloc()
          end if
!!$          call self%indmask(i)%p%writeFITS('TEST.fits')
!!$          call mpi_finalize(j)
!!$          stop
       end do
       call indmask%dealloc()
    end if

    ! Read deflation mask
    if (trim(cpar%cs_defmask(id_abs)) /= 'fullsky') then
       self%lmax_def      = 512
       self%nside_def     = 32
       self%fwhm_def      = 90.d0
       info_def => comm_mapinfo(self%comm, self%nside_def, self%lmax_def, 1, .false.)
       self%defmask => comm_map(info_def, trim(cpar%datadir)//'/'//trim(cpar%cs_defmask(id_abs)), udgrade=.true.)
       where (self%defmask%map < 0.5)
          self%defmask%map = 0.d0
       elsewhere
          self%defmask%map = 1.d0
       end where
    else
       self%lmax_def      = -1
       self%nside_def     = 0
       self%fwhm_def      = 0.d0
    end if

    ! Initialize prior mean
    if (trim(cpar%cs_prior_amp(id_abs)) /= 'none') then
       self%mu => comm_map(info, trim(cpar%datadir)//'/'//trim(cpar%cs_prior_amp(id_abs)))
       call self%mu%YtW
    end if

    ! Allocate mixing matrix
    call update_status(status, "init_pre_mix")
    self%ndet = maxval(data%ndet)
    allocate(self%F(numband,0:self%ndet), self%F_mean(numband,0:self%ndet,self%nmaps), self%F_null(numband,0:self%ndet))
    self%F_mean = 0.d0
    do i = 1, numband
       info      => comm_mapinfo(cpar%comm_chain, data(i)%info%nside, &
            & self%lmax_ind, data(i)%info%nmaps, data(i)%info%pol)
!       write(*,*) i, 'ndet = ', data(i)%ndet, shape(self%F), info%nside
       do j = 0, data(i)%ndet
          self%F(i,j)%p    => comm_map(info)
          self%F_null(i,j) =  .false.
       end do
    end do
    call update_status(status, "init_postmix")


    ! Initialize power spectrum
    self%Cl => comm_Cl(cpar, self%x%info, id, id_abs)

    ! Initialize pointers for non-linear search
    if (.not. allocated(res_smooth)) allocate(res_smooth(numband))
    if (.not. allocated(rms_smooth)) allocate(rms_smooth(numband))

    ! Set up monopole prior
    self%mono_prior_type = get_token(cpar%cs_mono_prior(id_abs), ":", 1)
    if (trim(self%mono_prior_type) /= 'none') then
       filename = get_token(cpar%cs_mono_prior(id_abs), ":", 2)
       self%mono_prior_map => comm_map(self%x%info, trim(cpar%datadir)//'/'//trim(filename))
    end if

    call update_status(status, "init_diffuse_end")
    
  end subroutine initDiffuse

  subroutine initLmaxSpecind(self, cpar, id, id_abs)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    integer(i4b) :: i, j, k, l, m, nmaps, p, ierr

    call update_status(status, "lmax_specind_start")

    !check nmaps of the component
    nmaps = 1
    if (cpar%cs_polarization(id_abs)) nmaps=3

    if (self%npar==0) return !do not go further, lmax_ind is set in initDiffuse 

    allocate(self%lmax_ind_pol(3,self%npar))  ! {integer}: lmax per. polarization (poltype index) per spec. ind.

    !self%lmax_ind      = cpar%cs_lmax_ind(id_abs) !not to be used anymore

    self%lmax_ind = -1 !default
    do i = 1,self%npar
       do p = 1, self%poltype(i)
          l = cpar%cs_lmax_ind_pol(p,i,id_abs)

          !assign lmax per spec ind per polarization (poltype)
          if (self%poltype(i)==1) then !all polarizations have the same lmax
             self%lmax_ind_pol(:,i) = l 
          else if (self%poltype(i)==2) then 
             if (p == 1) then 
                if (cpar%only_pol) cycle
                self%lmax_ind_pol(1,i) = l
             else
                !if only_pol, we don't want a possible Stokes I to slow down a potential lmax_ind == 0 speed up in SHT
                if (cpar%only_pol) self%lmax_ind_pol(1,i) = l 
                if (p > nmaps) then
                   self%lmax_ind_pol(2:3,i) = self%lmax_ind_pol(nmaps,i)
                   cycle
                else
                   self%lmax_ind_pol(2:3,i) = l
                end if
             end if
          else if (self%poltype(i)==3) then 
             if (p == 1) then 
                if (cpar%only_pol) cycle
                self%lmax_ind_pol(1,i) = l
             else if (p == 2) then 
                !if only_pol, we don't want a possible Stokes I to slow down a potential lmax_ind == 0 speed up in SHT
                if (cpar%only_pol) self%lmax_ind_pol(1,i) = l 
                if (p > nmaps) then
                   self%lmax_ind_pol(2,i) = self%lmax_ind_pol(nmaps,i)
                   cycle
                else
                   self%lmax_ind_pol(2,i) = l
                end if
             else
                if (p > nmaps) then
                   self%lmax_ind_pol(3,i) = self%lmax_ind_pol(nmaps,i)
                   cycle
                else
                   self%lmax_ind_pol(3,i) = l
                end if
             end if
          else
             write(*,fmt='(a,i1,a)') 'Incompatiple poltype ',p,' for '//trim(self%label)//' '//trim(self%indlabel(i))
             stop
          end if
          if (l > self%lmax_ind) self%lmax_ind = l
       end do
    end do


    allocate(self%lmax_ind_mix(3,self%npar))  ! {integer}: lmax per. polarization (poltype index) per spec. ind.
                                              ! Set to zero if lmax = -1 and fullsky pixregions are used

    self%lmax_ind_mix = self%lmax_ind_pol
    ! Comment: We are defining lmax per poltype index (or polarization) so that only the valid polarizations contribute
    !          to the overall lmax_ind. It also makes sure that if all valid lmax_ind parameters per poltype index and 
    !          parameter is 0, then all non-active are set to zero as well, i.e. we save some time on SHT.
    !          
    ! Note: This subroutine must be called for all diffuse components with npar > 0, before initializing diffuse comp,
    !       i.e. initDiffuse(), but after the number of spec. ind. parameters (npar) and poltype has been set.
    call update_status(status, "lmax_specind_end")

  end subroutine initLmaxSpecind

  subroutine initPixregSampling(self, cpar, id, id_abs)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    integer(i4b) :: i, j, k, l, m, n, p, ierr
    type(comm_mapinfo), pointer :: info => null()
    real(dp)           :: par_dp
    integer(i4b), allocatable, dimension(:) :: sum_pix
    real(dp),     allocatable, dimension(:) :: sum_theta, sum_proplen, sum_nprop
    real(dp),     allocatable, dimension(:,:) :: m_in, m_out, buffer
    character(len=512) :: temptxt, partxt
    integer(i4b) :: smooth_scale, p_min, p_max
    class(comm_mapinfo), pointer :: info2 => null(), info_ud => null()
    class(comm_map),     pointer :: tp => null() 
    class(comm_map),     pointer :: tp_smooth => null() 


    call update_status(status, "initPixreg_specind_start")
    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, self%nside, self%lmax_ind, &
         & self%nmaps, self%pol)



    ! Set up smoothing scale information
    allocate(self%smooth_scale(self%npar))
    self%smooth_scale = cpar%cs_smooth_scale(id_abs,1:self%npar)


    !!! (local) sampling specific parameters!!!

    allocate(self%pol_lnLtype(3,self%npar))        ! {chisq, ridge, marginal}: evaluation type for lnL
    allocate(self%pol_sample_nprop(3,self%npar))   ! {.true., .false.}: sample nprop on first iteration
    allocate(self%pol_sample_proplen(3,self%npar)) ! {.true., .false.}: sample proplen on first iteration
    allocate(self%pol_pixreg_type(3,self%npar))    ! {1=fullsky, 2=single_pix, 3=pixel_regions}
    allocate(self%nprop_uni(2,self%npar))          ! {integer}: upper and lower limits on nprop
    allocate(self%npixreg(3,self%npar))            ! {integer}: number of pixel regions per poltye per spec ind
    allocate(self%first_ind_sample(3,self%npar)) !used for pixelregion sampling
    self%first_ind_sample=.true.

    call update_status(status, "initPixreg_specind_pixreg_type")

    self%npixreg = 0
    do i = 1,self%npar
       self%nprop_uni(:,i)=cpar%cs_spec_uni_nprop(:,i,id_abs)
       do j = 1,self%poltype(i)
          self%pol_lnLtype(j,i)  = cpar%cs_spec_lnLtype(j,i,id_abs)
          self%pol_sample_nprop(j,i) = cpar%cs_spec_samp_nprop(j,i,id_abs)
          self%pol_sample_proplen(j,i) = cpar%cs_spec_samp_proplen(j,i,id_abs)

          if (self%lmax_ind_pol(j,i) < 0) then
             if (trim(cpar%cs_spec_pixreg(j,i,id_abs))=='fullsky') then
                self%pol_pixreg_type(j,i) = 1
                self%npixreg(j,i) = 1
             else if (trim(cpar%cs_spec_pixreg(j,i,id_abs))=='single_pix') then
                self%pol_pixreg_type(j,i) = 2
                if (cpar%nside_smooth(smooth_scale) < self%theta(i)%p%info%nside) then
                   self%npixreg(j,i) = 12*(cpar%nside_smooth(self%smooth_scale(i))**2)
                else
                   self%npixreg(j,i) = 12*(self%theta(i)%p%info%nside**2)
                end if
             else if (trim(cpar%cs_spec_pixreg(j,i,id_abs))=='pixreg') then
                self%pol_pixreg_type(j,i) = 3
                self%npixreg(j,i) = cpar%cs_spec_npixreg(j,i,id_abs) 
             else
                write(*,*) 'Unspecified pixel region type for poltype',j,'of spectral index',i,'in component',id_abs
                stop
             end if
          else
             !check if alm using defined pixel regions to suggest theta (then to smooth and get alms from map)
             !if so, pixreg values needs to be added. 
             if (cpar%almsamp_pixreg) then
                self%pol_pixreg_type(j,i) = 3
                self%npixreg(j,i) = cpar%cs_spec_npixreg(j,i,id_abs) 
             else
                self%pol_pixreg_type(j,i) = 0
             end if
          end if
       end do
    end do

    call update_status(status, "initPixreg_specind_allocate")

    if (any(self%pol_pixreg_type(:,:) > 0)) then
       !find highest number of pixel regions
       k = 0
       do i = 1,self%npar
          do j = 1,self%poltype(i)
             if (j > self%nmaps) cycle
             if (self%npixreg(j,i) > k) k = self%npixreg(j,i)
          end do
       end do
       allocate(self%nprop_pixreg(k,3,self%npar))
       allocate(self%npix_pixreg(k,3,self%npar))
       allocate(self%proplen_pixreg(k,3,self%npar))
       allocate(self%B_pp_fr(self%npar))
    end if

    ! Set up spectral index sampling masks, proposal length maps and nprop maps
    allocate(self%pol_ind_mask(self%npar)) ! masks per spectral index (all poltypes)
    allocate(self%pol_nprop(self%npar))    ! nprop map per spectral index (all poltypes
    allocate(self%pol_proplen(self%npar))  ! proplen map per spectral index (all poltypes)
    allocate(self%ind_pixreg_map(self%npar))   ! pixel region map per spectral index (all poltypes)
    if (any(self%pol_pixreg_type(:,:) > 0)) then
       k=0
       do i = 1,self%npar
          do j = 1,self%poltype(i)
             if (self%pol_pixreg_type(j,i) > 0) then
                if (self%npixreg(j,i) > k) k = self%npixreg(j,i)
             end if
          end do
       end do
       allocate(self%theta_pixreg(0:k,3,self%npar))
    end if

    info => comm_mapinfo(cpar%comm_chain, self%nside, self%lmax_ind, &
         & self%nmaps, self%pol)

    if (any(self%pol_pixreg_type > 0)) then
       if (any(self%poltype > 1)) then
          allocate(self%ind_pixreg_arr(0:info%np-1,3,self%npar)) 
       else
          allocate(self%ind_pixreg_arr(0:info%np-1,1,self%npar)) 
       end if
       self%ind_pixreg_arr = 0 !all pixels assigned to pixelregion 0 (not to be sampled), will read in pixreg later
    end if


    do i = 1,self%npar
       call update_status(status, "initPixreg_specind_mask")
       ! spec. ind. mask
       if (trim(cpar%cs_spec_mask(i,id_abs)) == 'fullsky') then
          self%pol_ind_mask(i)%p => comm_map(info)
          self%pol_ind_mask(i)%p%map = 1.d0
       else
          ! Read map from FITS file
          self%pol_ind_mask(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_spec_mask(i,id_abs)))

          if (min(self%poltype(i),self%nmaps) > &
               & self%pol_ind_mask(i)%p%info%nmaps) then
             write(*,fmt='(a,i2,a,i2,a,i2)') trim(self%indlabel(i))//' mask has fewer maps (', & 
                  & self%pol_ind_mask(i)%p%info%nmaps,') than poltype (',self%poltype(i), &
                  & ') for component nr. ',id_abs
             stop
          else if (self%theta(i)%p%info%nside /= self%pol_ind_mask(i)%p%info%nside) then
             write(*,fmt='(a,i4,a,i4,a,i2)') trim(self%indlabel(i))//' mask has different nside (', & 
                  & self%pol_ind_mask(i)%p%info%nside,') than the spectral index map (', &
                  & self%theta(i)%p%info%nside, ') for component nr. ',id_abs
             stop
          end if
       end if

       call update_status(status, "initPixreg_specind_proplen")

       ! prop. len. map
       if (trim(cpar%cs_spec_proplen(i,id_abs)) == 'fullsky') then
          self%pol_proplen(i)%p => comm_map(info)
          self%pol_proplen(i)%p%map = 1.d0
       else
          ! Read map from FITS file
          self%pol_proplen(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_spec_proplen(i,id_abs)))
          if (min(self%poltype(i),self%nmaps) > &
               & self%pol_proplen(i)%p%info%nmaps) then
             write(*,fmt='(a,i2,a,i2,a,i2)') trim(self%indlabel(i))//' proplen map has fewer maps (', & 
                  & self%pol_ind_mask(i)%p%info%nmaps,') than poltype (',self%poltype(i), &
                  & ') for component nr. ',id_abs
             stop
          else if (self%theta(i)%p%info%nside /= self%pol_proplen(i)%p%info%nside) then
             write(*,fmt='(a,i4,a,i4,a,i2)') trim(self%indlabel(i))//' proplen map has different nside (', & 
                  & self%pol_ind_mask(i)%p%info%nside,') than the spectral index map (', &
                  & self%theta(i)%p%info%nside, ') for component nr. ',id_abs
             stop
          end if
       end if
       
       ! replace proplen of input map for given poltype with user specified value, if given
       do j = 1,self%poltype(i)
          if (j > self%nmaps) cycle
          if (cpar%cs_spec_proplen_init(j,i,id_abs) > 0.d0) then
             self%pol_proplen(i)%p%map(:,j)=cpar%cs_spec_proplen_init(j,i,id_abs)
          end if
       end do

       call update_status(status, "initPixreg_specind_nprop")
    
       ! nprop map
       if (trim(cpar%cs_spec_nprop(i,id_abs)) == 'fullsky') then
          self%pol_nprop(i)%p => comm_map(info)
          self%pol_nprop(i)%p%map = 1.d0
       else
          ! Read map from FITS file
          self%pol_nprop(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_spec_nprop(i,id_abs)))
          if (min(self%poltype(i),self%nmaps) > &
               & self%pol_nprop(i)%p%info%nmaps) then
             write(*,fmt='(a,i2,a,i2,a,i2)') trim(self%indlabel(i))//' nprop map has fewer maps (', & 
                  & self%pol_ind_mask(i)%p%info%nmaps,') than poltype (',self%poltype(i), &
                  & ') for component nr. ',id_abs
             stop
          else if (self%theta(i)%p%info%nside /= self%pol_nprop(i)%p%info%nside) then
             write(*,fmt='(a,i4,a,i4,a,i2)') trim(self%indlabel(i))//' nprop map has different nside (', & 
                  & self%pol_ind_mask(i)%p%info%nside,') than the spectral index map (', &
                  & self%theta(i)%p%info%nside, ') for component nr. ',id_abs
             stop
          end if
       end if
       ! replace nprop of input map for given poltype with user specified value, if given       
       do j = 1,self%poltype(i)
          if (j > self%nmaps) cycle
          if (cpar%cs_spec_nprop_init(j,i,id_abs) > 0) then
             self%pol_nprop(i)%p%map(:,j)=cpar%cs_spec_nprop_init(j,i,id_abs)*1.d0
          end if
       end do

       ! limit nprop
       self%pol_nprop(i)%p%map = min(max(self%pol_nprop(i)%p%map,self%nprop_uni(1,i)*1.d0), &
            & self%nprop_uni(2,i)*1.d0)

       call update_status(status, "initPixreg_specind_pixel_regions")

       ! initialize pixel regions if relevant 
       if (any(self%pol_pixreg_type(1:self%poltype(i),i) > 0)) then

          info2  => comm_mapinfo(self%theta(i)%p%info%comm, self%nside, &
               & 2*self%nside, 1, .false.) 

          smooth_scale = self%smooth_scale(i)
          if (cpar%num_smooth_scales > 0 .and. smooth_scale > 0) then
             if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then
                !create beam for smoothing
                self%B_pp_fr(i)%p => comm_B_bl(cpar, info2, 1, 1, fwhm=cpar%fwhm_smooth(smooth_scale), &
                     & init_realspace=.false.)
             end if
          end if

          if (trim(cpar%cs_spec_pixreg_map(i,id_abs)) == 'fullsky') then
             self%ind_pixreg_map(i)%p => comm_map(info)
             self%ind_pixreg_map(i)%p%map = 1.d0
          else if (trim(cpar%cs_spec_pixreg_map(i,id_abs)) == 'none') then
             self%ind_pixreg_map(i)%p => comm_map(info)
             self%ind_pixreg_map(i)%p%map = 0.d0
          else
             ! Read map from FITS file
             self%ind_pixreg_map(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_spec_pixreg_map(i,id_abs)))
             if (min(self%poltype(i),self%nmaps) > &
                  & self%ind_pixreg_map(i)%p%info%nmaps) then
                write(*,fmt='(a,i2,a,i2,a,i2)') trim(self%indlabel(i))//' pixreg map has fewer maps (', & 
                     & self%pol_ind_mask(i)%p%info%nmaps,') than poltype (',self%poltype(i), &
                     & ') for component nr. ',id_abs
                stop
             else if (self%theta(i)%p%info%nside /= self%ind_pixreg_map(i)%p%info%nside) then
                write(*,fmt='(a,i4,a,i4,a,i2)') trim(self%indlabel(i))//' pixreg map has different nside (', & 
                     & self%pol_ind_mask(i)%p%info%nside,') than the spectral index map (', &
                     & self%theta(i)%p%info%nside, ') for component nr. ',id_abs
                stop
             end if
          end if

          ! Check if pixel region type = 'fullsky' og 'single_pix' for given poltype index  
          do j = 1,self%poltype(i)
             if (j > self%nmaps) cycle
             if (self%pol_pixreg_type(j,i) == 1) then !fullsky
                self%ind_pixreg_map(i)%p%map(:,j) = 1.d0
             else if (self%pol_pixreg_type(j,i) == 2) then 
                !single pix at smoothing scale pix size
                !make a map at smooth scale nside, give pixels value from 1 to npix
                !udgrade to full resolution
                if (cpar%nside_smooth(smooth_scale) > self%nside) then
                   !only allow nside up to full resolution
                   info_ud  => comm_mapinfo(self%theta(i)%p%info%comm, &
                        & self%theta(i)%p%info%nside, -1, 1, .false.) 
                else
                   info_ud  => comm_mapinfo(self%theta(i)%p%info%comm, &
                        & cpar%nside_smooth(smooth_scale), -1, 1, .false.) 
                end if

                if (cpar%nside_smooth(smooth_scale) < self%theta(i)%p%info%nside) then
                   tp => comm_map(info_ud)
                   tp%map(:,1) = tp%info%pix*1.d0 + 1.d0 !set pixel value equal to pixel number+1
                   allocate(m_in(0:info_ud%npix-1,1))
                   allocate(m_out(0:self%ind_pixreg_map(i)%p%info%npix-1, 1))
                   allocate(buffer(0:self%ind_pixreg_map(i)%p%info%npix-1, 1))
                   m_in(tp%info%pix,1) = tp%map(:,1)
                   call udgrade_ring(m_in, info_ud%nside, m_out, &
                        & self%ind_pixreg_map(i)%p%info%nside)
                   call mpi_allreduce(m_out, buffer, size(m_out), MPI_DOUBLE_PRECISION, &
                        & MPI_SUM, self%ind_pixreg_map(i)%p%info%comm, ierr)

                   self%ind_pixreg_map(i)%p%map(:,j) = buffer(self%ind_pixreg_map(i)%p%info%pix,1)
                   deallocate(m_in, m_out, buffer)
                   call tp%dealloc()
                   tp => null()
                else
                   self%ind_pixreg_map(i)%p%map(:,j) = self%ind_pixreg_map(i)%p%info%pix*1.d0 +1.d0
                end if
             end if
          end do

          ! check if there is a non-smoothed init map for theta of pixelregions
          if (cpar%cs_pixreg_init_theta(i,id_abs) == 'none') then
             tp => comm_map(self%theta(i)%p%info)
             tp%map = self%theta(i)%p%map !take avrage from existing theta map
          else
             !read map from init map (non-smoothed theta map)
             tp => comm_map(self%theta(i)%p%info, trim(cpar%datadir) // '/' // trim(cpar%cs_pixreg_init_theta(i,id_abs)))
          end if

          !compute the average theta in each pixel region for the poltype indices that sample theta using pixel regions
          do j = 1,self%poltype(i)
             if (j > self%nmaps) cycle
             self%theta_pixreg(:,j,i)=self%p_gauss(1,i) !prior
             if (self%pol_pixreg_type(j,i) < 1) cycle

             n=self%npixreg(j,i)
             allocate(sum_pix(n),sum_theta(n),sum_nprop(n),sum_proplen(n))
             sum_theta=0.d0
             sum_pix=0
             sum_nprop=0.d0
             sum_proplen=0.d0

             do k = 0,self%theta(i)%p%info%np-1
                do m = 1,n
                   if ( self%ind_pixreg_map(i)%p%map(k,j) > (m-0.5d0) .and. &
                        & self%ind_pixreg_map(i)%p%map(k,j) < (m+0.5d0) ) then
                      !sum_theta(m)=sum_theta(m)+self%theta(i)%p%map(k,j)
                      sum_theta(m)=sum_theta(m)+tp%map(k,j)
                      sum_proplen(m)=sum_proplen(m)+self%pol_proplen(i)%p%map(k,j)
                      sum_nprop(m)=sum_nprop(m)+self%pol_nprop(i)%p%map(k,j)
                      sum_pix(m)=sum_pix(m)+1
                      self%ind_pixreg_arr(k,j,i)=m !assign pixel region index 
                      exit
                   end if
                end do
             end do

             !allreduce
             call mpi_allreduce(MPI_IN_PLACE, sum_theta, n, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
             call mpi_allreduce(MPI_IN_PLACE, sum_proplen, n, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
             call mpi_allreduce(MPI_IN_PLACE, sum_nprop, n, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
             call mpi_allreduce(MPI_IN_PLACE, sum_pix, n, MPI_INTEGER, MPI_SUM, info%comm, ierr)

             do k = 1,n
                !if (cpar%myid == cpar%root) write(*,*) 'pixreg',k,'  -- numbe of pixels',sum_pix(k)
                if (sum_pix(k) > 0) then
                   self%theta_pixreg(k,j,i) = sum_theta(k)/(1.d0*sum_pix(k))
                   self%nprop_pixreg(k,j,i) = IDINT(sum_nprop(k)/(1.d0*sum_pix(k)))
                   self%proplen_pixreg(k,j,i) = sum_proplen(k)/(1.d0*sum_pix(k))
                else
                   self%theta_pixreg(k,j,i) = self%p_gauss(1,i) ! the prior as theta
                   self%nprop_pixreg(k,j,i) = 0
                   self%proplen_pixreg(k,j,i) = 1.d0
                end if
                self%npix_pixreg(k,j,i)=sum_pix(k)
             end do
             self%theta_pixreg(0,j,i)=self%p_gauss(1,i) !all pixels in region 0 has the prior as theta
             
             deallocate(sum_pix,sum_theta,sum_proplen,sum_nprop)
             
             !Should also assign (and smooth) theta map. This might be dangerous so we don't assign for now
             ! might put in a flag for this if wanted
             if (.false.) then
                if (self%poltype(i) == 1) then
                   p_min = 1; p_max = self%nmaps
                   if (cpar%only_pol) p_min = 2
                else if (self%poltype(i) == 2) then
                   if (j == 1) then
                      p_min = 1; p_max = 1
                   else
                      p_min = 2; p_max = self%nmaps
                   end if
                else if (self%poltype(i) == 3) then
                   p_min = j
                   p_max = j
                else
                   write(*,*) 'Unsupported polarization type'
                   stop
                end if

                smooth_scale = self%smooth_scale(i)
                if (cpar%num_smooth_scales > 0 .and. smooth_scale > 0) then
                   if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then
                      if (p_min<=p_max) then !just a guaranty that we dont smooth for nothing

                         !spec. ind. map with 1 map (will be smoothed like zero spin map using the existing code)
                         tp => comm_map(info2)

                         do k = 0,info2%np-1
                            tp%map(k,1) = self%theta_pixreg(self%ind_pixreg_arr(k,j,i),j,i)
                         end do
                         !smooth single map as intensity (i.e. zero spin)
                         call smooth_map(info2, .false., &
                              & self%B_pp_fr(i)%p%b_l*0.d0+1.d0, tp, &  
                              & self%B_pp_fr(i)%p%b_l, tp_smooth)

                         do k = p_min,p_max
                            self%theta(i)%p%map(:,k) = tp_smooth%map(:,1)
                         end do
                         call tp_smooth%dealloc()
                         call tp%dealloc()

                      end if
                   end if
                end if !num smooth scales > 0
             end if
          end do !poltype

          call tp%dealloc()

       end if !any pixreg_type > 0

       call update_status(status, "initPixreg_specind_lmax_ind_mix")

       !final check to see if lmax < 0 and pixel region is fullsky
       do j = 1,self%poltype(i)
          if (j > self%nmaps) then
             self%lmax_ind_mix(j:,i) = self%lmax_ind_mix(self%nmaps,i)
             cycle
          end if

          if (self%poltype(i) == 1) then
             p_min = 1; p_max = 3
          else if (self%poltype(i) == 2) then
             if (j == 1) then
                if (cpar%only_pol) cycle
                p_min = 1; p_max = 1
             else
                p_min = 2; p_max = 3
                if (cpar%only_pol) p_min = 1 !add the first index to this case
             end if
          else if (self%poltype(i) == 3) then
             if (cpar%only_pol .and. j == 1) cycle
             p_min = j
             p_max = j
             if (cpar%only_pol .and. j == 2) p_min = 1 !add the first index to this case
          else
             write(*,*) 'Unsupported polarization type'
             stop
          end if
             
          if (self%lmax_ind_pol(j,i) >= 0) then
             self%lmax_ind_mix(p_min:p_max,i) = self%lmax_ind_pol(j,i) !in case only_pol and poltype = 2 has lmax > 0
          else if (self%pol_pixreg_type(j,i)==1) then !pixel region is defined fullsky
             self%lmax_ind_mix(p_min:p_max,i) = 0
          else if (self%pol_pixreg_type(j,i) == 3) then
             !not to sure about these checks, omitting them for now, these could mess up debug
             if (.false.) then
                if (trim(cpar%cs_spec_pixreg_map(i,id_abs)) == 'fullsky') then !pixel region map is a fullsky map
                   self%lmax_ind_mix(p_min:p_max,i) = 0
                else
                   !check if any pixel region has as many pixels as the full theta map (i.e. only one pixel region is defined)
                   if (any(self%npix_pixreg(:,j,i)==self%theta(i)%p%info%npix)) then 
                      self%lmax_ind_mix(p_min:p_max,i) = 0
                   end if
                end if
             end if
          end if
       end do

    end do !npar

    if (self%theta(1)%p%info%myid == 0 .and. .false.) then !debug

       write(*,*) ''
       write(*,*) 'lmax_ind_pol'
       do i = 1,self%npar
          write(*,*) i, self%lmax_ind_pol(:,i)
       end do
       write(*,*) ''
       write(*,*) 'lmax_ind_mix'
       do i = 1,self%npar
          write(*,*) i, self%lmax_ind_mix(:,i)
       end do
       write(*,*) ''
       if (all(self%lmax_ind_mix(:,:) == 0)) write(*,*) '  all lmax_ind_mix zero'

    end if

    call update_status(status, "initPixreg_specind_end")

  end subroutine initPixregSampling

  subroutine initSpecindProp(self,cpar, id, id_abs)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    character(len=512) :: filename
    integer(i4b) :: i, j, l, m, ntot, nloc, p
    real(dp) :: fwhm_prior, sigma_prior
    logical(lgt) :: exist



    ! Init alm sampling params (Trygve)
    allocate(self%corrlen(self%npar, self%nmaps))
    self%corrlen    = 0     ! Init correlation length
    
    ! Init bool for L_read
    allocate(self%L_read(self%npar))
    self%L_read    = .false.

    ! save minimum chisq per iteration
    allocate(self%chisq_min(self%npar, self%nmaps))

    self%nalm_tot = (self%lmax_ind + 1)**2

    ! Init smooth prior
    allocate(self%sigma_priors(0:self%nalm_tot-1,self%npar)) !a_00 is given by different one

    fwhm_prior = cpar%almsamp_prior_fwhm   !600.d0 ! 1200.d0
    do j = 1, self%npar
       self%sigma_priors(0,j) = 0.05 !p_gauss(2,j)*0.1
       if (self%nalm_tot > 1) then
          ! Saving and smoothing priors
          i = 1
          do l = 1, self%lmax_ind ! Skip a_00 - m-major ordering (0,0)(1,0)(2,0)(1,1)(1,-1)(2,1)(2,-1)(2,2)(2,-2)
             self%sigma_priors(i,j) = 0.01*self%sigma_priors(0,j)*exp(-0.5d0*l*(l+1)*(fwhm_prior * pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
             i = i + 1
          end do
          do l = 1, self%lmax_ind
             do m = 1, l
                self%sigma_priors(i,j) = 0.01*self%sigma_priors(0,j)*exp(-0.5d0*l*(l+1)*(fwhm_prior * pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
                i = i + 1
                self%sigma_priors(i,j) = 0.01*self%sigma_priors(0,j)*exp(-0.5d0*l*(l+1)*(fwhm_prior * pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
                i = i + 1
             end do
          end do
       end if
    end do

    ! Initialize cholesky matrix
    if (cpar%almsamp_pixreg)  then
       allocate(self%L(0:50, 0:50, self%nmaps, self%npar)) ! Set arbitrary number of max regions
    else
       allocate(self%L(0:self%nalm_tot-1, 0:self%nalm_tot-1, self%nmaps, self%npar))
    end if

    allocate(self%steplen(self%nmaps, self%npar)) !a_00 is given by different one              
    self%L = 0.d0
    self%steplen = 1 !0.3d0

    ! Filename formatting
    do j = 1, self%npar
       filename = trim(cpar%datadir)//'/init_alm_'//trim(self%label)//'_'//trim(self%indlabel(j))//'.dat'

       inquire(file=filename, exist=exist)
       if (exist) then ! If present cholesky file
          self%L_read(j) = .true.
          if (self%myid == 0) write(*,*) " - ALM init file found for "//trim(self%label)//" "//trim(self%indlabel(j))
          open(unit=11, file=filename, recl=10000)
          read(11,*) self%corrlen(j,:)
          read(11,*) self%L(:,:,:,j)
          close(11)
       else
          if (self%myid == 0) write(*,*) " - ALM init file NOT found for "//trim(self%label)//" "//trim(self%indlabel(j))
          if (cpar%almsamp_pixreg) then
             self%L(:,:,:,j) = self%sigma_priors(0,j)
          else
             do p = 0, self%nalm_tot-1
                self%L(p,p,:,j) = self%sigma_priors(p,j)
             end do
          end if
       end if
    end do
    
  end subroutine initSpecindProp


  subroutine initDiffPrecond(comm)
    implicit none
    integer(i4b),                intent(in) :: comm

    select case (trim(precond_type))
    case ("diagonal")
       call initDiffPrecond_diagonal(comm)
    case ("pseudoinv")
       call initDiffPrecond_pseudoinv(comm)
    case default
       call report_error("Preconditioner type not supported")
    end select

  end subroutine initDiffPrecond

  subroutine initDiffPrecond_diagonal(comm)
    implicit none
    integer(i4b),                intent(in) :: comm

    integer(i4b) :: i, i1, i2, j, k1, k2, q, l, m, n
    real(dp)     :: t1, t2
    integer(i4b), allocatable, dimension(:) :: ind
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: p1 => null(), p2 => null()
    real(dp),     allocatable, dimension(:,:) :: mat

    call update_status(status, "init_diffpre1")
    if (npre == 0) return
    if (allocated(P_cr%invM_diff)) return
    
    if (.not. allocated(diffComps)) then
       ! Set up an array of all the diffuse components
       allocate(diffComps(npre))
       c => compList
       i =  1
       do while (associated(c))
          select type (c)
          class is (comm_diffuse_comp)
             diffComps(i)%p => c
             i              =  i+1
          end select
          c => c%next()
       end do
       info_pre => comm_mapinfo(comm, nside_pre, lmax_pre, nmaps_pre, nmaps_pre==3)
    end if
    call update_status(status, "init_diffpre2")
    
    ! Build frequency-dependent part of preconditioner
    call wall_time(t1)
    allocate(P_cr%invM_diff(0:info_pre%nalm-1,info_pre%nmaps))
    !!$OMP PARALLEL PRIVATE(mat, ind, j, i1, l, m, q, i2, k1, p1, k2, n)
    allocate(mat(npre,npre), ind(npre))
    do j = 1, info_pre%nmaps
       !!$OMP DO SCHEDULE(guided)
       do i1 = 0, info_pre%nalm-1
          call info_pre%i2lm(i1, l, m)
          mat = 0.d0
          do q = 1, numband
             call data(q)%info%lm2i(l,m,i2)
             if (i2 == -1) cycle
             if (j > data(q)%info%nmaps) cycle
             do k1 = 1, npre
                p1 => diffComps(k1)%p
                if (l > p1%lmax_amp) cycle
                if (j > p1%nmaps) cycle
                do k2 = 1, npre
                   p2 => diffComps(k2)%p
                   if (l > p2%lmax_amp) cycle
                   if (j > p2%nmaps) cycle
                   mat(k1,k2) = mat(k1,k2) + &
                        & data(q)%N%invN_diag%alm(i2,j) * &  ! invN_{lm,lm}
                        & data(q)%B(0)%p%b_l(l,j)**2 * &          ! b_l^2
                        & p1%F_mean(q,0,j) * p2%F_mean(q,0,j)    ! F(c1)*F(c2)
                end do
             end do
          end do

          n = 0
          allocate(P_cr%invM_diff(i1,j)%comp2ind(npre))
          P_cr%invM_diff(i1,j)%comp2ind = -1
          do k1 = 1, npre
             if (mat(k1,k1) > 0.d0) then
                n = n+1
                ind(n) = k1
                P_cr%invM_diff(i1,j)%comp2ind(k1) = n
             end if
          end do
          P_cr%invM_diff(i1,j)%n = n
          allocate(P_cr%invM_diff(i1,j)%ind(n))
          allocate(P_cr%invM_diff(i1,j)%M0(n,n), P_cr%invM_diff(i1,j)%M(n,n))
          P_cr%invM_diff(i1,j)%ind = ind(1:n)
          P_cr%invM_diff(i1,j)%M0   = mat(ind(1:n),ind(1:n))
       end do
       !!$OMP END DO
    end do
    deallocate(ind, mat)
    !!$OMP END PARALLEL
    call wall_time(t2)
    call update_status(status, "init_diffpre3")

  end subroutine initDiffPrecond_diagonal


  subroutine initDiffPrecond_pseudoinv(comm)
    implicit none
    integer(i4b),                intent(in) :: comm

    integer(i4b) :: i, i1, i2, j, k1, k2, q, l, m, n
    real(dp)     :: t1, t2
    integer(i4b), allocatable, dimension(:) :: ind
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: p1 => null(), p2 => null()
    real(dp),     allocatable, dimension(:,:) :: mat

    if (npre == 0) return
    if (allocated(P_cr%invM_diff)) return
    
    if (.not. allocated(diffComps)) then
       ! Set up an array of all the diffuse components
       allocate(diffComps(npre))
       c => compList
       i =  1
       do while (associated(c))
          select type (c)
          class is (comm_diffuse_comp)
             diffComps(i)%p => c
             i              =  i+1
          end select
          c => c%next()
       end do
       info_pre => comm_mapinfo(comm, nside_pre, lmax_pre, nmaps_pre, nmaps_pre==3)
    end if
    
    ! Allocate space for pseudo-inverse of U
    call wall_time(t1)
    allocate(P_cr%invM_diff(0:lmax_pre,info_pre%nmaps))
    do j = 1, info_pre%nmaps
       do l = 0, lmax_pre
          allocate(P_cr%invM_diff(l,j)%M(npre,numband+npre))
       end do
    end do

  end subroutine initDiffPrecond_pseudoinv


  subroutine updateDiffPrecond(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update

    select case (trim(precond_type))
    case ("diagonal")
       call updateDiffPrecond_diagonal(samp_group, force_update)
    case ("pseudoinv")
       call updateDiffPrecond_pseudoinv(samp_group, force_update)
    case default
       call report_error("Preconditioner type not supported")
    end select

  end subroutine updateDiffPrecond

  subroutine updateDiffPrecond_diagonal(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update

    integer(i4b) :: i, j, k, k1, k2, l, m, n, ierr, unit, p, q
    real(dp)     :: W_ref, t1, t2
    logical(lgt), save :: first_call = .true.
    integer(i4b), allocatable, dimension(:) :: ind
    real(dp),     allocatable, dimension(:) :: W
    real(dp),     allocatable, dimension(:) :: cond, W_min
    real(dp),     allocatable, dimension(:,:) :: alm

    if (npre == 0) return
    if (.not. recompute_diffuse_precond .and. .not. force_update) return
    
    ! Initialize current preconditioner to F^t * B^t * invN * B * F
    !self%invM    = self%invM0
    do j = 1, info_pre%nmaps
       do i = 0, info_pre%nalm-1
          if (P_cr%invM_diff(i,j)%n > 0) P_cr%invM_diff(i,j)%M = P_cr%invM_diff(i,j)%M0
       end do
    end do

!!$    if (info_pre%myid == 0) then
!!$       n = self%invM_(0,1)%n
!!$       allocate(W(n))
!!$       call get_eigenvalues(self%invM_(0,1)%M, W)
!!$       write(*,*) 'W0 = ', real(W,sp)
!!$       deallocate(W)
!!$       write(*,*) 
!!$       do i = 1, 3
!!$          write(*,*) real(self%invM_(0,1)%M(i,:),sp)
!!$       end do
!!$       write(*,*)
!!$    end if

    ! Right-multiply with sqrt(Cl)
    call wall_time(t1)
    do k1 = 1, npre
       if (trim(diffComps(k1)%p%cltype) == 'none') cycle
       !$OMP PARALLEL PRIVATE(alm, k2, j, i, p, q)
       allocate(alm(0:info_pre%nalm-1,info_pre%nmaps))
       !$OMP DO SCHEDULE(guided)
       do k2 = 1, npre
          do j = 1, info_pre%nmaps
             do i = 0, info_pre%nalm-1
                if (P_cr%invM_diff(i,j)%n == 0) cycle
                p = P_cr%invM_diff(i,j)%comp2ind(k1)
                q = P_cr%invM_diff(i,j)%comp2ind(k2)
                if (p /= -1 .and. q /= -1) then
                   alm(i,j) = P_cr%invM_diff(i,j)%M(q,p)
                else
                   alm(i,j) = 0.d0
                end if
             end do
          end do
          !if (info_pre%myid == 0) write(*,*) 'a', k1, k2, alm(4,1)
          call diffComps(k1)%p%Cl%sqrtS(alm=alm, info=info_pre, diag=.true.)
          !if (info_pre%myid == 0) write(*,*) 'b', k1, k2, alm(4,1)
          !call mpi_finalize(j)
          !stop
          do j = 1, info_pre%nmaps
             do i = 0, info_pre%nalm-1
                if (P_cr%invM_diff(i,j)%n == 0) cycle                
                p = P_cr%invM_diff(i,j)%comp2ind(k1)
                q = P_cr%invM_diff(i,j)%comp2ind(k2)
                if (p /= -1 .and. q /= -1) P_cr%invM_diff(i,j)%M(q,p) = alm(i,j)
             end do
          end do
       end do
       !$OMP END DO
       deallocate(alm)
       !$OMP END PARALLEL
    end do

    ! Left-multiply with sqrt(Cl)
    do k1 = 1, npre
       if (trim(diffComps(k1)%p%cltype) == 'none') cycle
       !$OMP PARALLEL PRIVATE(alm, k2, j, i, p, q)
       allocate(alm(0:info_pre%nalm-1,info_pre%nmaps))
       !$OMP DO SCHEDULE(guided)
       do k2 = 1, npre
          do j = 1, info_pre%nmaps
             do i = 0, info_pre%nalm-1
                if (P_cr%invM_diff(i,j)%n == 0) cycle                
                p = P_cr%invM_diff(i,j)%comp2ind(k1)
                q = P_cr%invM_diff(i,j)%comp2ind(k2)
                if (p /= -1 .and. q /= -1) then
                   alm(i,j) = P_cr%invM_diff(i,j)%M(p,q)
                else
                   alm(i,j) = 0.d0
                end if
             end do
          end do
!          if (info_pre%myid == 0) write(*,*) 'c', k1, k2, alm(4,1)
          call diffComps(k1)%p%Cl%sqrtS(alm=alm, info=info_pre, diag=.true.)
!          if (info_pre%myid == 0) write(*,*) 'd', k1, k2, alm(4,1)
          do j = 1, info_pre%nmaps
             do i = 0, info_pre%nalm-1
                if (P_cr%invM_diff(i,j)%n == 0) cycle                
                p = P_cr%invM_diff(i,j)%comp2ind(k1)
                q = P_cr%invM_diff(i,j)%comp2ind(k2)
                if (p /= -1 .and. q /= -1) P_cr%invM_diff(i,j)%M(p,q) = alm(i,j)
             end do
          end do
       end do
       !$OMP END DO
       deallocate(alm)
       !$OMP END PARALLEL
    end do
    !call wall_time(t2)
    !write(*,*) 'sqrtS = ', t2-t1

    ! Nullify temperature block if only polarization
    if (only_pol) then
       do i = 0, info_pre%nalm-1
          if (P_cr%invM_diff(i,1)%n == 0) cycle                
          P_cr%invM_diff(i,1)%M = 0.d0
       end do
    end if

!!$    if (info_pre%myid == 0) then
!!$       n = self%invM_(0,1)%n
!!$       allocate(W(n))
!!$       call get_eigenvalues(self%invM_(0,1)%M, W)
!!$       write(*,*) 'W1 = ', real(W,sp)
!!$       deallocate(W)
!!$       write(*,*) 
!!$       do i = 1, 3
!!$          write(*,*) real(self%invM_(0,1)%M(i,:),sp)
!!$       end do
!!$       write(*,*)
!!$    end if


    ! Add unity 
    do k1 = 1, npre
       if (trim(diffComps(k1)%p%cltype) == 'none') cycle
       !!$OMP PARALLEL PRIVATE(i,l,m,j,p)
       !!$OMP DO SCHEDULE(guided)
       do i = 0, info_pre%nalm-1
          call info_pre%i2lm(i, l, m)
          !if (info_pre%myid == 1) then
          !write(*,*) info_pre%myid, i, l, m
          !end if
          if (l <= diffComps(k1)%p%lmax_amp) then
             do j = 1, info_pre%nmaps
                if (P_cr%invM_diff(i,j)%n == 0) cycle                
                p = P_cr%invM_diff(i,j)%comp2ind(k1)
                if (p > 0) P_cr%invM_diff(i,j)%M(p,p) = P_cr%invM_diff(i,j)%M(p,p) + 1.d0
             end do
          end if
       end do
       !!$OMP END DO
       !!$OMP END PARALLEL
    end do

!!$    if (info_pre%myid == 0) then
!!$       n = self%invM_(0,1)%n
!!$       allocate(W(n))
!!$       call get_eigenvalues(self%invM_(0,1)%M, W)
!!$       write(*,*) 'W2 = ', real(W,sp)
!!$       deallocate(W)
!!$    end if

!!$    call mpi_finalize(i)
!!$    stop


    ! Nullify elements that are not involved in current sample group
    do j = 1, info_pre%nmaps
       do i = 0, info_pre%nalm-1
          if (P_cr%invM_diff(i,j)%n == 0) cycle                
          do k1 = 1, npre
             if (diffComps(k1)%p%active_samp_group(samp_group)) cycle
             p = P_cr%invM_diff(i,j)%comp2ind(k1)
             if (p == -1) cycle
             P_cr%invM_diff(i,j)%M(p,:) = 0.d0
             P_cr%invM_diff(i,j)%M(:,p) = 0.d0
          end do
       end do
    end do

    ! Print out worst condition number in first call
    call wall_time(t1)
    if (first_call .and. output_cg_eigenvals) then
       allocate(cond(0:lmax_pre), W_min(0:lmax_pre))
       cond  = 0.d0
       W_min = 1.d30
       do j = 1, nmaps_pre
          do i = 0, info_pre%nalm-1
             call info_pre%i2lm(i, l, m)
             n = P_cr%invM_diff(i,j)%n
             if (n > 0) then
                allocate(W(n), ind(n))
                call get_eigenvalues(P_cr%invM_diff(i,j)%M, W)
                W_ref    = minval(abs(W))
                cond(l)  = max(cond(l), maxval(abs(W/W_ref)))
                W_min(l) = min(W_min(l), minval(W))
                deallocate(W, ind)
             end if
          end do
       end do
       call mpi_allreduce(MPI_IN_PLACE, cond,  lmax_pre+1, MPI_DOUBLE_PRECISION, MPI_MAX, info_pre%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, W_min, lmax_pre+1, MPI_DOUBLE_PRECISION, MPI_MIN, info_pre%comm, ierr)
       if (info_pre%myid == 0) then
          write(*,*) 'Precond -- largest condition number = ', real(maxval(cond),sp)
          write(*,*) 'Precond -- smallest eigenvalue      = ', real(minval(W_min),sp)
          unit = getlun()
          open(unit, file=trim(outdir)//'/precond_eigenvals.dat', recl=1024)
          do l = 0, lmax_pre
             !write(unit,fmt='(i8,2e16.8)') l, cond(l), W_min(l)
             write(unit,*) l, cond(l), W_min(l)
          end do
          close(unit)
       end if

       first_call = .false.
       deallocate(cond, W_min)
    end if
    call wall_time(t2)
    !write(*,*) 'eigen = ', t2-t1

    ! Invert matrix
    call wall_time(t1)
    do j = 1, nmaps_pre
       do i = 0, info_pre%nalm-1
!!$          if (info_pre%myid == 0) then
!!$             do k = 1, P_cr%invM_diff(i,j)%n
!!$                write(*,*) real(P_cr%invM_diff(i,j)%M(k,:),sp)
!!$             end do
!!$          end if
          if (P_cr%invM_diff(i,j)%n > 0) then
             if (any(P_cr%invM_diff(i,j)%M /= 0.d0)) call invert_matrix_with_mask(P_cr%invM_diff(i,j)%M)
          end if
       end do
    end do
    call wall_time(t2)
    !write(*,*) 'invert = ', t2-t1

    ! Disable preconditioner update
    recompute_diffuse_precond = .false.
       
  end subroutine updateDiffPrecond_diagonal


  subroutine updateDiffPrecond_pseudoinv(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update

    integer(i4b) :: i, j, k, l, n, ierr, p, q
    real(dp)     :: t1, t2, Cl
    real(dp),     allocatable, dimension(:,:) :: mat

    if (npre == 0) return
    if (.not. recompute_diffuse_precond .and. .not. force_update) return

    ! Build pinv_U
    call wall_time(t1)
    allocate(mat(numband+npre,npre))
    do j = 1, info_pre%nmaps
       do l = 0, lmax_pre

          ! Data section
          mat = 0.d0
          do q = 1, numband
             if (l > data(q)%info%lmax) cycle
             if (j > data(q)%info%nmaps) cycle
             do k = 1, npre
                if (l > diffComps(k)%p%lmax_amp) cycle
                if (j > diffComps(k)%p%nmaps) cycle
                if (.not. diffComps(k)%p%active_samp_group(samp_group)) cycle
                mat(q,k) = data(q)%N%alpha_nu(j) * &
                          & data(q)%B(0)%p%b_l(l,j) * &
                          & diffComps(k)%p%F_mean(q,0,j)

                if (trim(diffComps(k)%p%Cl%type) /= 'none') then
                   mat(q,k) = mat(q,k)*sqrt(diffComps(k)%p%Cl%getCl(l,j))
                end if
                if (mat(q,k) /= mat(q,k)) then
                   if (trim(diffComps(k)%p%Cl%type) /= 'none') then
                      write(*,*) q, j, l, k, data(q)%N%alpha_nu(j), data(q)%B(0)%p%b_l(l,j), &
                           & diffComps(k)%p%F_mean(q,0,j),diffComps(k)%p%Cl%getCl(l,j)
                   else
                      write(*,*) q, j, l, k, data(q)%N%alpha_nu(j), data(q)%B(0)%p%b_l(l,j),&
                           & diffComps(k)%p%F_mean(q,0,j)
                   end if

                end if
             end do
          end do

          ! Prior section
          do k = 1, npre
             if (trim(diffComps(k)%p%Cl%type) == 'none' .or. l > diffComps(k)%p%lmax_amp .or. &
                  & .not. diffComps(k)%p%active_samp_group(samp_group)) cycle
             mat(numband+k,k) = 1.d0
          end do

          ! Store pseudo-inverse of U
          call compute_pseudo_inverse(mat, P_cr%invM_diff(l,j)%M)
          !P_cr%invM_diff(l,j)%M = transpose(mat)

!!$          !if (info_pre%myid == 0 .and. l > 4498 .and. l < 4503) then
!!$          if (info_pre%myid == 0) then
!!$             write(*,*)
!!$             do k = 1, numband+npre
!!$                write(*,*) mat(k,:)
!!$             end do
!!$             write(*,*) shape(mat)
!!$             do k = 1, npre
!!$                write(*,*) P_cr%invM_diff(l,j)%M(k,:)
!!$             end do
!!$          end if

       end do
    end do
    deallocate(mat)
    call wall_time(t2)

    ! Set up array of active components
    if (allocated(ind_pre)) deallocate(ind_pre)
    n = 0
    do k = 1, npre
       if (diffComps(k)%p%active_samp_group(samp_group)) n = n+1
    end do
    if (n > 0) then
       allocate(ind_pre(n))
       i = 0
       do k = 1, npre
          if (diffComps(k)%p%active_samp_group(samp_group)) then
             i = i+1
             ind_pre(i) = k
          end if
       end do
    end if

!!$    call mpi_finalize(k)
!!$    stop

    ! Disable preconditioner update
    recompute_diffuse_precond = .false.

  end subroutine updateDiffPrecond_pseudoinv
    
  
  ! Evaluate amplitude map in brightness temperature at reference frequency
  subroutine updateDiffuseMixmat(self, theta, beta, band, df, par)
    implicit none
    class(comm_diffuse_comp),                  intent(inout)           :: self
    class(comm_map),           dimension(:),   intent(in),    optional :: theta
    real(dp),  dimension(:,:,:),               intent(in),    optional :: beta  ! Not used here
    integer(i4b),                              intent(in),    optional :: band
    class(map_ptr), dimension(:),              intent(inout), optional :: df    ! Derivative of mixmat with respect to parameter par; for Jeffreys prior
    integer(i4b),                              intent(in),    optional :: par   ! Parameter ID for derivative

    integer(i4b) :: i, j, k, l, n, p, p_min, p_max, nmaps, ierr
    real(dp)     :: lat, lon, t1, t2
    logical(lgt) :: precomp, mixmatnull, bad ! NEW
    character(len=2) :: ctext
    real(dp),        allocatable, dimension(:,:,:) :: theta_p
    real(dp),        allocatable, dimension(:)     :: nu, s, buffer, buff2
    class(comm_mapinfo),          pointer          :: info => null(), info_tp => null()
    class(comm_map),              pointer          :: t => null(), tp => null(), td => null()
    class(map_ptr),  allocatable, dimension(:)     :: theta_prev


    if (trim(self%type) == 'md') return

    call update_status(status, "mixupdate1 " // trim(self%label))
    
    ! Copy over alms from input structure, and compute pixel-space parameter maps
    if (present(theta)) then
       do i = 1, self%npar
          self%theta(i)%p%alm = theta(i)%alm          
       end do
    end if
    
    ! Compute mixing matrix
!!$    allocate(theta_prev(self%npar))
!!$    do j = 1, self%npar
!!$       nullify(theta_prev(j)%p)
!!$    end do

    do i = 1, numband
       
       ! Only update requested band if present
       if (present(band)) then
          if (i /= band) cycle
       end if

       ! Compute spectral parameters at the correct resolution for this channel
       if (self%npar > 0) then
          nmaps = min(data(i)%info%nmaps, self%theta(1)%p%info%nmaps)
          allocate(theta_p(0:data(i)%info%np-1,nmaps,self%npar))

          do j = 1,self%npar
             info => comm_mapinfo(data(i)%info%comm, data(i)%info%nside, &
                  & self%theta(j)%p%info%lmax, nmaps, data(i)%info%pol)
             td => comm_map(info)
             
             ! if any polarization is alm sampled. Only use alms to set polarizations with alm sampling
             if (any(self%lmax_ind_pol(1:self%poltype(j),j) >= 0)) then
                t => comm_map(info)
                t%alm(:,1:nmaps) = self%theta(j)%p%alm(:,1:nmaps)
                call t%Y_scalar()
                do p = 1,self%poltype(j)
                   if (self%lmax_ind_pol(p,j) < 0) cycle
                   if (self%poltype(j) == 1) then
                      p_min=1
                      p_max=info%nmaps
                      if (only_pol) p_min = 2
                   else if (self%poltype(j)==2) then
                      if (p == 1) then
                         p_min = 1
                         p_max = 1
                      else
                         p_min = 2
                         p_max = info%nmaps
                      end if
                   else if (self%poltype(j)==3) then
                      p_min = p
                      p_max = p
                   else
                      write(*,*) '  Unknown poltype in component ',self%label,', parameter ',self%indlabel(j) 
                      stop
                   end if

                   do k = p_min,p_max
                      td%map(:,k) = t%map(:,k)
                   end do
                end do
                call t%dealloc()
             end if

             ! if any polarization is local sampled. Only set theta using polarizations with local sampling
             if (any(self%lmax_ind_pol(1:self%poltype(j),j) < 0)) then
                info_tp => comm_mapinfo(self%theta(j)%p%info%comm, self%theta(j)%p%info%nside, &
                     & 2*self%theta(j)%p%info%nside, nmaps, .false.)
                tp    => comm_map(info_tp)
                tp%map(:,1:nmaps) = self%theta(j)%p%map(:,1:nmaps)
                call tp%YtW_scalar()
                info => comm_mapinfo(data(i)%info%comm, data(i)%info%nside, &
                     & data(i)%info%lmax, nmaps, .false.)
                t    => comm_map(info)
                call tp%alm_equal(t)
                call t%Y_scalar()
                where (t%map < self%p_uni(1,j))
                   t%map = self%p_uni(1,j)
                end where
                where (t%map > self%p_uni(2,j))
                   t%map = self%p_uni(2,j)
                end where

!!$                call tp%udgrade(t)
!!$
                bad = any(t%map == 0.d0)
                call mpi_allreduce(mpi_in_place, bad, 1, &
                     & MPI_LOGICAL, MPI_LOR, self%x%info%comm, ierr)
                if (bad) then
                   write(*,*) trim(self%label), i, j
                   call int2string(j, ctext)
                   call t%writeFITS("beta1.fits")
                   call tp%writeFITS("beta2.fits")
                   call mpi_finalize(k)
                   stop
                end if

                call tp%dealloc()

                
                call wall_time(t2)

                do p = 1,self%poltype(j)
                   if (self%lmax_ind_pol(p,j) >= 0) cycle
                   if (self%poltype(j) == 1) then
                      p_min=1
                      p_max=info%nmaps
                      if (only_pol) p_min = 2
                   else if (self%poltype(j)==2) then
                      if (p == 1) then
                         p_min = 1
                         p_max = 1
                      else
                         p_min = 2
                         p_max = info%nmaps
                      end if
                   else if (self%poltype(j)==3) then
                      p_min = p
                      p_max = p
                   else
                      write(*,*) '  Unknown poltype in component ',self%label,', parameter ',self%indlabel(j) 
                      stop
                   end if

                   do k = p_min,p_max
                      td%map(:,k) = t%map(:,k)
                   end do
                end do

                call t%dealloc()

                !if (info%myid == 0) write(*,*) 'udgrade = ', t2-t1
             end if
             theta_p(:,:,j) = td%map
             !write(*,*) 'q1', minval(theta_p(:,:,j)), maxval(theta_p(:,:,j))
!!$             if (associated(theta_prev(j)%p)) call theta_prev(j)%p%dealloc()
!!$             theta_prev(j)%p => comm_map(t)
             call td%dealloc()
          end do
       end if


       do l = 0, data(i)%ndet
          
          ! Don't update null mixing matrices
          if (self%F_null(i,l)) then
             if (present(df)) df(i)%p%map = 0.d0
             cycle
          end if
          
!!$          ! Compute spectral parameters at the correct resolution for this channel
!!$          if (self%npar > 0) then
!!$             nmaps = min(data(i)%info%nmaps, self%theta(1)%p%info%nmaps)
!!$             allocate(theta_p(0:data(i)%info%np-1,nmaps,self%npar))
!!$             
!!$             do j = 1, self%npar
!!$                if (.false. .and. associated(theta_prev(j)%p)) then
!!$                   if (data(i)%info%nside == theta_prev(j)%p%info%nside .and. &
!!$                        & self%theta(j)%p%info%lmax == theta_prev(j)%p%info%lmax  .and. &
!!$                        & nmaps == theta_prev(j)%p%info%nmaps .and. &
!!$                        & data(i)%info%pol == theta_prev(j)%p%info%pol) then !
!!$                      theta_p(:,:,j) = theta_prev(j)%p%map
!!$                      cycle
!!$                   end if
!!$                end if
!!$                info => comm_mapinfo(data(i)%info%comm, data(i)%info%nside, self%theta(j)%p%info%lmax, nmaps, data(i)%info%pol)
!!$                t    => comm_map(info)
!!$                if (self%lmax_ind >= 0) then
!!$                   t%alm(:,1:nmaps) = self%theta(j)%p%alm(:,1:nmaps)
!!$                   call t%Y_scalar
!!$                else
!!$                   call wall_time(t1)
!!$                   call self%theta(j)%p%udgrade(t)
!!$                   call wall_time(t2)
!!$                   !if (info%myid == 0) write(*,*) 'udgrade = ', t2-t1
!!$                end if
!!$                theta_p(:,:,j) = t%map
!!$                if (associated(theta_prev(j)%p)) call theta_prev(j)%p%dealloc()
!!$                theta_prev(j)%p => comm_map(t)
!!$                call t%dealloc()
!!$             end do
!!$          end if
          
          
          ! Loop over all pixels, computing mixing matrix for each
          !allocate(theta_p(self%npar,self%nmaps))
          call wall_time(t1)
          do j = 0, self%F(i,l)%p%info%np-1
             if (self%latmask >= 0.d0) then
                call pix2ang_ring(data(i)%info%nside, data(i)%info%pix(j+1), lat, lon)
                lat = 0.5d0*pi-lat
                if (abs(lat) < self%latmask) then
                   self%F(i,l)%p%map(j,:) = 0.d0
                   cycle
                end if
             end if

!          if (all(self%lmax_ind_mix(1:min(self%nmaps,data(i)%info%nmaps)) == 0)) then  !if (self%lmax_ind == 0) then
!             cycle
!          end if

!!$          if (self%npar > 0) then
!!$             ! Collect all parameters
!!$             t0 => t
!!$             theta_p(1,:) = t0%map(j,:)
!!$             do k = 2, self%npar
!!$                t0 => t0%next()
!!$                theta_p(k,:) = t0%map(j,:)
!!$             end do
!!$
!!$             ! Check polarization type
!!$             if (self%nmaps == 3) then
!!$                do k = 1, self%npar
!!$                   if (self%poltype(k) < 2) theta_p(k,2) = theta_p(k,1)
!!$                   if (self%poltype(k) < 3) theta_p(k,3) = theta_p(k,2)
!!$                end do
!!$             end if
!!$          end if

             ! NEW ! Check band sensitivity before mixing matrix update
             ! Possible labels are "broadband", "cmb", "synch", "dust", "co10", "co21", "co32", "ff", "ame"
             if (data(i)%comp_sens == "broadband") then
                ! If broadband, calculate mixing matrix
                mixmatnull = .false.
             else
                ! If component sensitivity, only calculate mixmat on that component.
                mixmatnull = .true.
                If (data(i)%comp_sens == self%label) then
                   mixmatnull = .false.
                end if
             end if
             
             ! Temperature
             if (self%npar > 0) then
                if (mixmatnull) then
                   self%F(i,l)%p%map(j,1) = 0.0
                else
                   if (trim(self%label) == 'dust' .and. any(theta_p(j,1,:)==0.d0)) then
                      write(*,*) i, l, j, real(theta_p(j,1,:),sp)
                      stop
                   end if
                   self%F(i,l)%p%map(j,1) = self%F_int(1,i,l)%p%eval(theta_p(j,1,:)) * data(i)%gain * self%cg_scale
                end if
             else
                if (mixmatnull) then 
                   self%F(i,l)%p%map(j,1) = 0.0
                else
                   self%F(i,l)%p%map(j,1) = self%F_int(1,i,l)%p%eval([0.d0]) * data(i)%gain * self%cg_scale
                end if
             end if
             
             ! Polarization
             if (self%nmaps == 3 .and. data(i)%info%nmaps == 3) then
                ! Stokes Q
                if (self%npar == 0) then
                   self%F(i,l)%p%map(j,2) = self%F(i,l)%p%map(j,1) 
                else if (all(self%poltype < 2)) then
                   self%F(i,l)%p%map(j,2) = self%F(i,l)%p%map(j,1) 
                else
                   if (trim(self%label) == 'dust' .and. any(theta_p(j,1,:)==0.d0)) then
                      write(*,*) i, l, j, real(theta_p(j,1,:),sp)
                   end if
                   if (self%npar > 0) then
                      self%F(i,l)%p%map(j,2) = self%F_int(2,i,l)%p%eval(theta_p(j,2,:)) * data(i)%gain * self%cg_scale
                   else
                      self%F(i,l)%p%map(j,2) = self%F_int(2,i,l)%p%eval([0.d0]) * data(i)%gain * self%cg_scale
                   end if
                end if
                
                ! Stokes U
                if (self%npar == 0) then
                   self%F(i,l)%p%map(j,3) = self%F(i,l)%p%map(j,2) 
                else if (all(self%poltype < 3)) then
                   self%F(i,l)%p%map(j,3) = self%F(i,l)%p%map(j,2) 
                else
                   if (trim(self%label) == 'dust' .and. any(theta_p(j,1,:)==0.d0)) then
                      write(*,*) i, l, j, real(theta_p(j,1,:),sp)
                   end if
                   if (self%npar > 0) then
                      self%F(i,l)%p%map(j,3) = self%F_int(3,i,l)%p%eval(theta_p(j,3,:)) * data(i)%gain * self%cg_scale
                   else
                      self%F(i,l)%p%map(j,3) = self%F_int(3,i,l)%p%eval([0.d0]) * data(i)%gain * self%cg_scale
                   end if
                end if
             end if
          
             ! Compute derivative if requested
             if (present(df) .and. l == 0) then
                if (self%npar > 0) then
                   do k = 1, nmaps
                      if (k <= self%poltype(par)) then
                         df(i)%p%map(j,k) = self%F_int(k,i,l)%p%eval_deriv(theta_p(j,k,:),par) * data(i)%gain * self%cg_scale
                      end if
                   end do
                else
                   df(i)%p%map = 0.d0
                end if
             end if

          end do

          call wall_time(t2)
          !if (self%x%info%myid == 0) write(*,*) 'eval = ', t2-t1
                
          ! Compute mixing matrix average; for preconditioning only
          allocate(buffer(self%nmaps), buff2(self%nmaps))
          do j = 1, min(self%nmaps, data(i)%info%nmaps)
             self%F_mean(i,l,j) = sum(self%F(i,l)%p%map(:,j))
          end do
          buff2 = self%F_mean(i,l,:)
          call mpi_allreduce(buff2, buffer, self%nmaps, &
               & MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
          self%F_mean(i,l,:) = buffer / self%F(i,l)%p%info%npix
          deallocate(buffer,buff2)

!!$       if (self%npar > 0) then
!!$          call t%dealloc(clean_info=.true.)
!!$          call t0%dealloc(clean_info=.true.)
!!$          nullify(t)
!!$          nullify(t0)
!!$       end if
    
       end do
       if (allocated(theta_p)) deallocate(theta_p)
    end do
!!$    do j = 1, self%npar
!!$       if (associated(theta_prev(j)%p)) call theta_prev(j)%p%dealloc()
!!$    end do
!!$    deallocate(theta_prev)


    call update_status(status, "mixupdate2 " // trim(self%label))

    ! Request preconditioner update
    recompute_diffuse_precond = .true.

  end subroutine updateDiffuseMixmat



  function evalDiffuseBand(self, band, amp_in, pix, alm_out, det)
    implicit none
    class(comm_diffuse_comp),                     intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
    real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
    logical(lgt),                                 intent(in),  optional :: alm_out
    integer(i4b),                                 intent(in),  optional :: det
    real(dp),        dimension(:,:), allocatable                        :: evalDiffuseBand

    integer(i4b) :: i, j, np, nmaps, lmax, nmaps_comp, d
    logical(lgt) :: alm_out_
    real(dp)     :: t1, t2
    class(comm_mapinfo), pointer :: info 
    class(comm_map),     pointer :: m

    alm_out_ = .false.; if (present(alm_out)) alm_out_ = alm_out
    d = 0; if (present(det)) d = det

    if (self%F_null(band,0)) then
       if (alm_out_) then
          if (.not. allocated(evalDiffuseBand)) allocate(evalDiffuseBand(0:data(band)%info%nalm-1,data(band)%info%nmaps))
       else
          if (.not. allocated(evalDiffuseBand)) allocate(evalDiffuseBand(0:data(band)%info%np-1,data(band)%info%nmaps))
       end if
       evalDiffuseBand = 0.d0
       return
    end if

    ! Initialize amplitude map
    nmaps =  min(data(band)%info%nmaps, self%nmaps)
    info  => comm_mapinfo(data(band)%info%comm, data(band)%info%nside, data(band)%info%lmax, nmaps, nmaps==3)
    m     => comm_map(info)
    if (present(amp_in)) then
       m%alm(:,1:nmaps) = amp_in(:,1:nmaps)
       !m%alm(:,1:nmaps) = amp_in
    else
       call self%x%alm_equal(m)
       !m%alm(:,1:nmaps) = self%x%alm(:,1:nmaps)
    end if
!!$    if (self%lmax_amp > data(band)%map%info%lmax) then
!!$       ! Nullify elements above band-specific lmax to avoid aliasing during projection
!!$       do i = 0, m%info%nalm-1
!!$          if (m%info%lm(1,i) > data(band)%info%lmax) m%alm(i,:) = 0.d0
!!$       end do
!!$    end if

    if (apply_mixmat) then
       ! Scale to correct frequency through multiplication with mixing matrix

       if (all(self%lmax_ind_mix(1:nmaps,:) == 0) .and. self%latmask < 0.d0) then
          do i = 1, m%info%nmaps
             m%alm(:,i) = m%alm(:,i) * self%F_mean(band,d,i)
          end do
       else
          call m%Y()
          m%map(:,1:nmaps) = m%map(:,1:nmaps) * self%F(band,d)%p%map(:,1:nmaps)
          call m%YtW()
       end if
    end if
       
    ! Convolve with band-specific beam
    call data(band)%B(d)%p%conv(trans=.false., map=m)
    if (.not. alm_out_) call m%Y()

    ! Return correct data product
    if (alm_out_) then
       !if (.not. allocated(evalDiffuseBand)) allocate(evalDiffuseBand(0:self%x%info%nalm-1,self%x%info%nmaps))
       if (.not. allocated(evalDiffuseBand)) allocate(evalDiffuseBand(0:data(band)%info%nalm-1,data(band)%info%nmaps))
       if (nmaps /= data(band)%info%nmaps) evalDiffuseBand = 0.d0
       evalDiffuseBand(:,1:nmaps) = m%alm(:,1:nmaps)
    else
       if (.not. allocated(evalDiffuseBand)) allocate(evalDiffuseBand(0:data(band)%info%np-1,data(band)%info%nmaps))
       if (nmaps /= data(band)%info%nmaps) evalDiffuseBand = 0.d0
       evalDiffuseBand(:,1:nmaps) = m%map(:,1:nmaps)
    end if
       

    ! Clean up
    call m%dealloc()
    nullify(info)

  end function evalDiffuseBand

  ! Return component projected from map
  function projectDiffuseBand(self, band, map, alm_in, det)
    implicit none
    class(comm_diffuse_comp),                     intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    class(comm_map),                              intent(in)            :: map
    logical(lgt),                                 intent(in), optional  :: alm_in
    integer(i4b),                                 intent(in), optional  :: det
    real(dp),        dimension(:,:), allocatable                        :: projectDiffuseBand

    integer(i4b) :: i, nmaps, d
    logical(lgt) :: alm_in_
    class(comm_mapinfo), pointer :: info_in => null(), info_out => null()
    class(comm_map),     pointer :: m       => null(), m_out    => null()

    if (self%F_null(band,0)) then
       if (.not. allocated(projectDiffuseBand)) allocate(projectDiffuseBand(0:self%x%info%nalm-1,self%x%info%nmaps))
       projectDiffuseBand = 0.d0
       return
    end if

    d = 0; if (present(det)) d = det

    nmaps     =  min(self%x%info%nmaps, map%info%nmaps)
    info_in   => comm_mapinfo(self%x%info%comm, map%info%nside, map%info%lmax, nmaps, nmaps==3)
    info_out  => comm_mapinfo(self%x%info%comm, map%info%nside, self%lmax_amp, nmaps, nmaps==3)
    m         => comm_map(info_in)
    m_out     => comm_map(info_out)
    alm_in_ = .false.; if (present(alm_in)) alm_in_ = alm_in

    ! Convolve with band-specific beam
    if (alm_in) then
       m%alm = map%alm(:,1:nmaps)
    else
       m%map = map%map(:,1:nmaps)
       call m%Yt()
    end if
    call data(band)%B(d)%p%conv(trans=.true., map=m)
    
    if (all(self%lmax_ind_mix(1:nmaps,:) == 0) .and. self%latmask < 0.d0) then
       do i = 1, nmaps
          m%alm(:,i) = m%alm(:,i) * self%F_mean(band,d,i)
       end do
    else
       call m%Y()
       m%map(:,1:nmaps) = m%map(:,1:nmaps) * self%F(band,d)%p%map(:,1:nmaps)
       call m%YtW()
    end if
    call m%alm_equal(m_out)

    if (.not. allocated(projectDiffuseBand)) allocate(projectDiffuseBand(0:self%x%info%nalm-1,self%x%info%nmaps))
    projectDiffuseBand = m_out%alm

    call m%dealloc()
    call m_out%dealloc()

  end function projectDiffuseBand


  subroutine applyDiffPrecond(x)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x

    select case (trim(precond_type))
    case ("diagonal")
       call applyDiffPrecond_diagonal(x)
    case ("pseudoinv")
       call applyDiffPrecond_pseudoinv(x)
    case default
       call report_error("Preconditioner type not supported")
    end select

  end subroutine applyDiffPrecond


  subroutine applyDiffPrecond_diagonal(x)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x

    integer(i4b)              :: i, j, k, l, m, nmaps
    real(dp), allocatable, dimension(:,:)   :: alm
    real(dp), allocatable, dimension(:,:,:) :: y

    if (npre == 0) return
    
    ! Reformat linear array into y(npre,nalm,nmaps) structure
    allocate(y(npre,0:info_pre%nalm-1,info_pre%nmaps))
    y = 0.d0
    do i = 1, npre
       nmaps = diffComps(i)%p%x%info%nmaps
       call cr_extract_comp(diffComps(i)%p%id, x, alm)
       do j = 0, diffComps(i)%p%x%info%nalm-1
          call diffComps(i)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          y(i,k,1:nmaps) = alm(j,1:nmaps)
       end do
       deallocate(alm)
    end do

    ! Multiply with preconditioner
    do j = 1, nmaps_pre
       do i = 0, info_pre%nalm-1
          if (P_cr%invM_diff(i,j)%n == 0) cycle
          y(P_cr%invM_diff(i,j)%ind,i,j) = &
               & matmul(P_cr%invM_diff(i,j)%M, y(P_cr%invM_diff(i,j)%ind,i,j))
       end do
    end do

    ! Reformat y(npre,nalm,nmaps) structure into linear array
    do i = 1, npre
       nmaps = diffComps(i)%p%x%info%nmaps
       allocate(alm(0:diffComps(i)%p%x%info%nalm-1,nmaps))
       alm = 0.d0
       do j = 0, diffComps(i)%p%x%info%nalm-1
          call diffComps(i)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          alm(j,1:nmaps) = y(i,k,1:nmaps)
       end do
       call cr_insert_comp(diffComps(i)%p%id, .false., alm, x)
       deallocate(alm)
    end do
    
    deallocate(y)

  end subroutine applyDiffPrecond_diagonal


  subroutine applyDiffPrecond_pseudoinv(x)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x

    integer(i4b)              :: i, ii, j, k, l, m, p, q, qq, nmaps, npre_int
    real(dp)                  :: t1, t2
    real(dp), allocatable, dimension(:)     :: w, w2
    real(dp), allocatable, dimension(:,:)   :: alm
    real(dp), allocatable, dimension(:,:,:) :: y, z
    class(comm_map), pointer                :: invN_x => null()

    if (.not. allocated(ind_pre)) return

    npre_int = size(ind_pre)
    if (npre_int <= 0) return
    
    ! Reformat linear array into y(npre,nalm,nmaps) structure
!    call update_status(status, "pseudo1")
    allocate(y(npre_int,0:info_pre%nalm-1,info_pre%nmaps))
    allocate(z(npre_int,0:info_pre%nalm-1,info_pre%nmaps))
    y = 0.d0
    do i = 1, npre_int
       ii = ind_pre(i)
       nmaps = diffComps(ii)%p%x%info%nmaps
       call cr_extract_comp(diffComps(ii)%p%id, x, alm)
       do j = 0, diffComps(ii)%p%x%info%nalm-1
          call diffComps(ii)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          y(i,k,1:nmaps) = alm(j,1:nmaps)
       end do
       deallocate(alm)
    end do

    ! Frequency-dependent terms
    z = 0.d0
    do k = 1, numband
!       call update_status(status, "pseudo2")
       invN_x => comm_map(data(k)%info)
       nmaps  =  data(k)%info%nmaps
       
       ! Sum over (U^plus)^t
       !!$OMP PARALLEL DEFAULT(shared) PRIVATE(i,l,m,j,q,qq,p)
       !!$OMP DO SCHEDULE(guided)
       do i = 0, data(k)%info%nalm-1
          call data(k)%info%i2lm(i, l, m)
          if (l > info_pre%lmax) cycle
          call info_pre%lm2i(l,m,j)
          do q = 1, npre_int
             qq = ind_pre(q)
             do p = 1, nmaps
                invN_x%alm(i,p) = invN_x%alm(i,p) + P_cr%invM_diff(l,p)%M(qq,k) * y(q,j,p)
             end do
          end do
       end do
       !!$OMP END DO
       !!$OMP END PARALLEL
!       call update_status(status, "pseudo3")

       ! Multiply by T
       call invN_x%WY
!       call update_status(status, "pseudo4")
       call data(k)%N%N(invN_x)
!       call update_status(status, "pseudo5")
       call invN_x%YtW
!       call update_status(status, "pseudo6")
       do i = 1, nmaps
          invN_x%alm(:,i) = invN_x%alm(:,i) * data(k)%N%alpha_nu(i)**2
       end do

       ! Sum over U^plus
       !!$OMP PARALLEL DEFAULT(shared) PRIVATE(i,l,m,j,q,qq,p)
       !!$OMP DO SCHEDULE(guided)
       do i = 0, data(k)%info%nalm-1
          call data(k)%info%i2lm(i, l, m)
          if (l > info_pre%lmax) cycle
          call info_pre%lm2i(l,m,j)
          do q = 1, npre_int
             qq = ind_pre(q)
             do p = 1, nmaps
                z(q,j,p) = z(q,j,p) + P_cr%invM_diff(l,p)%M(qq,k) * invN_x%alm(i,p)
             end do
          end do
       end do
       !!$OMP END DO
       !!$OMP END PARALLEL
!       call update_status(status, "pseudo7")

       call invN_x%dealloc()
    end do

    ! Prior terms
!    call update_status(status, "pseudo7.1")
    call wall_time(t1)
    !!$OMP PARALLEL DEFAULT(shared) PRIVATE(i,l,k,j,m,w,w2,p)
    allocate(w(npre_int), w2(npre_int))
    !!$OMP DO SCHEDULE(guided)
    do i = 0, info_pre%nalm-1
       do p = 1, info_pre%nmaps
          !call info_pre%i2lm(i, l, m)
          l = info_pre%lm(1,i)
          m = info_pre%lm(2,i)
          w        = y(:,i,p)
          w2       = 0.d0
          do j = 1, npre_int
             do k = 1, npre_int
                w2(j) = w2(j) + P_cr%invM_diff(l,p)%M(ind_pre(k),numband+ind_pre(j))*w(k)
             end do
          end do
          w       = 0.d0
          do j = 1, npre_int
             do k = 1, npre_int
                w(j) = w(j) + P_cr%invM_diff(l,p)%M(ind_pre(j),numband+ind_pre(k))*w2(k)
             end do
          end do
          z(:,i,p) = z(:,i,p) + w
       end do
    end do
    deallocate(w, w2)
    !!$OMP END DO
    !!$OMP END PARALLEL
    call wall_time(t2)
    !if (info_pre%myid == 0 .or. info_pre%myid == 25) write(*,*) info_pre%myid, ', nalm = ', info_pre%nalm, real(t2-t1,sp)
!    call update_status(status, "pseudo8")

    ! Reformat z(npre,nalm,nmaps) structure into linear array
    do i = 1, npre_int
       ii = ind_pre(i)
       nmaps = diffComps(ii)%p%x%info%nmaps
       allocate(alm(0:diffComps(ii)%p%x%info%nalm-1,nmaps))
       alm = 0.d0
       do j = 0, diffComps(ii)%p%x%info%nalm-1
          call diffComps(ii)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          alm(j,1:nmaps) = z(i,k,1:nmaps)
       end do
       call cr_insert_comp(diffComps(ii)%p%id, .false., alm, x)
       deallocate(alm)
    end do
!    call update_status(status, "pseudo9")
    
    deallocate(y, z)

  end subroutine applyDiffPrecond_pseudoinv


  
  ! Dump current sample to HEALPix FITS file
  subroutine dumpDiffuseToFITS(self, iter, chainfile, output_hdf, postfix, dir)
    implicit none
    class(comm_diffuse_comp),                intent(inout)        :: self
    integer(i4b),                            intent(in)           :: iter
    type(hdf_file),                          intent(in)           :: chainfile
    logical(lgt),                            intent(in)           :: output_hdf
    character(len=*),                        intent(in)           :: postfix
    character(len=*),                        intent(in)           :: dir

    integer(i4b)       :: i, l, j, k, m, ierr, unit
    integer(i4b)       :: p, p_min, p_max, npr, npol
    real(dp)           :: vals(10)
    logical(lgt)       :: exist, first_call = .true.
    character(len=6)   :: itext
    character(len=512) :: filename, path
    class(comm_mapinfo), pointer :: info => null()
    class(comm_map), pointer :: map => null(), tp => null()
    real(dp), allocatable, dimension(:,:) :: sigma_l
    real(dp),     allocatable, dimension(:,:) :: dp_pixreg
    integer(i4b), allocatable, dimension(:,:) :: int_pixreg

    if (trim(self%type) == 'md') then
       filename = trim(dir)//'/md_' // trim(postfix) // '.dat'
       vals = 0.d0
       do i = 0, self%x%info%nalm-1
          call self%x%info%i2lm(i,l,m)
          if (l == 0) then                 ! Monopole
             vals(1)  = 1.d0/sqrt(4.d0*pi) * self%x%alm(i,1)             * self%RJ2unit_(1)
             vals(5)  = 1.d0/sqrt(4.d0*pi) * self%mu%alm(i,1)            * self%RJ2unit_(1)
             vals(9)  = 1.d0/sqrt(4.d0*pi) * sqrt(self%Cl%Dl(0,1))       * self%RJ2unit_(1)
          end if
          if (l == 1 .and. m == -1) then   ! Y dipole
             vals(3)  = 1.d0/sqrt(4.d0*pi/3.d0) * self%x%alm(i,1)        * self%RJ2unit_(1)
             vals(7)  = 1.d0/sqrt(4.d0*pi/3.d0) * self%mu%alm(i,1)       * self%RJ2unit_(1)
          end if
          if (l == 1 .and. m ==  0) then   ! Z dipole
             vals(4)  = 1.d0/sqrt(4.d0*pi/3.d0) * self%x%alm(i,1)        * self%RJ2unit_(1)
             vals(8)  = 1.d0/sqrt(4.d0*pi/3.d0) * self%mu%alm(i,1)       * self%RJ2unit_(1)
          end if
          if (l == 1 .and. m ==  1) then   ! X dipole
             vals(2)  = -1.d0/sqrt(4.d0*pi/3.d0) * self%x%alm(i,1)        * self%RJ2unit_(1)
             vals(6)  = -1.d0/sqrt(4.d0*pi/3.d0) * self%mu%alm(i,1)       * self%RJ2unit_(1)
             vals(10) =  1.d0/sqrt(4.d0*pi/3.d0) * sqrt(self%Cl%Dl(1,1))  * self%RJ2unit_(1)
          end if
       end do
       call mpi_allreduce(MPI_IN_PLACE, vals, 10, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       if (self%x%info%myid == 0) then
          if (output_hdf) then
             call int2string(iter, itext)
             path = trim(adjustl(itext))//'/md/'//trim(adjustl(self%label))
             call write_hdf(chainfile, trim(adjustl(path)), vals(1:4))
          end if
          
          inquire(file=trim(filename), exist=exist)
          unit = getlun()
          if (first_call) then
             open(unit, file=trim(filename), recl=1024)
             write(unit,*) '# Band          M_def       X_def       Y_def'//&
                  & '      Z_def     M_mean      X_mean      Y_mean      Z_mean'// &
                  & '      M_rms       D_rms'
             first_call = .false.
          else
             open(unit, file=trim(filename), recl=1024, position='append')
          end if
          write(unit,'(a12,10e12.3)') trim(self%label), vals
          close(unit)
       end if
    else

       call update_status(status, "writeFITS_1")

       ! Write amplitude
       map => comm_map(self%x)
       do i = 1, map%info%nmaps
          map%alm(:,i) = map%alm(:,i) * self%RJ2unit_(i) * self%cg_scale  ! Output in requested units
       end do

       !call update_status(status, "writeFITS_2")

       if (output_hdf) then
          call int2string(iter, itext)
          path = '/'//trim(adjustl(itext))//'/'//trim(adjustl(self%label))
          !write(*,*) 'path', trim(path)
          if (self%x%info%myid == 0) call create_hdf_group(chainfile, trim(adjustl(path)))
       end if

       !call update_status(status, "writeFITS_3")

       filename = trim(self%label) // '_' // trim(postfix) // '.fits'
       call self%B_out%conv(trans=.false., map=map)
       call map%Y
       do i = 1, map%info%nmaps
          map%alm(:,i) = self%x%alm(:,i) * self%RJ2unit_(i) * self%cg_scale  ! Replace convolved with original alms
       end do
       !call update_status(status, "writeFITS_4")

       !call self%apply_proc_mask(map)

       if (output_hdf) then
          !write(*,*) 'path2', trim(path)//'/amp_'
          call map%writeFITS(trim(dir)//'/'//trim(filename), &
               & hdffile=chainfile, hdfpath=trim(path)//'/amp_', output_hdf_map=.false.)
       else
          call map%writeFITS(trim(dir)//'/'//trim(filename))
       end if
       call map%dealloc()
       call update_status(status, "writeFITS_5")

       if (self%output_EB) then
          map => comm_map(self%x)

          do i = 1, map%info%nmaps
             map%alm(:,i) = map%alm(:,i) * self%RJ2unit_(i) * self%cg_scale  ! Output in requested units
          end do
          
          filename = trim(self%label) // '_' // trim(postfix) // '_TEB.fits'
          call self%B_out%conv(trans=.false., map=map)
          call map%Y_EB
          !call self%apply_proc_mask(map)
          call map%writeFITS(trim(dir)//'/'//trim(filename))
          call map%dealloc()
       end if
       !call update_status(status, "writeFITS_6")
       
       allocate(sigma_l(0:self%x%info%lmax,self%x%info%nspec))
       call self%x%getSigmaL(sigma_l)
       k = 1
       do i = 1, self%x%info%nmaps
          do j = i, self%x%info%nmaps
             sigma_l(:,k) = sigma_l(:,k) * self%RJ2unit(j)
             k = k+1
          end do
       end do
       if (self%x%info%myid == 0) then
          filename = trim(dir)//'/sigma_l_'//trim(self%label) // '_' // trim(postfix) // '.dat'
          call write_sigma_l(filename, sigma_l)
          if (output_hdf) call write_hdf(chainfile, trim(adjustl(path))//'/sigma_l', sigma_l)             
       end if
       deallocate(sigma_l)
       !call update_status(status, "writeFITS_7")

       ! Write spectral index maps
       do i = 1, self%npar
          filename = trim(self%label) // '_' // trim(self%indlabel(i)) // '_' // &
               & trim(postfix) // '.fits'

          ! need to only take alms from polarizations with lmax > 0
          if (self%lmax_ind >= 0) then
             info => comm_mapinfo(self%theta(i)%p%info%comm, self%theta(i)%p%info%nside, &
                  & self%theta(i)%p%info%lmax, self%theta(i)%p%info%nmaps, self%theta(i)%p%info%pol)
             
             ! if any polarization is alm sampled. Only use alms to set polarizations with alm sampling
             if (all(self%lmax_ind_pol(1:self%poltype(i),i) >= 0)) then
                call self%theta(i)%p%Y_scalar()
             else if (any(self%lmax_ind_pol(1:self%poltype(i),i) >= 0)) then
                tp => comm_map(info)
                tp%alm = self%theta(i)%p%alm
                call tp%Y_scalar
                do p = 1,self%poltype(i)
                   if (self%lmax_ind_pol(p,i) < 0) cycle
                   if (self%poltype(i) == 1) then
                      p_min=1
                      p_max=info%nmaps
                      if (only_pol) p_min = 2
                   else if (self%poltype(i)==2) then
                      if (p == 1) then
                         p_min = 1
                         p_max = 1
                      else
                         p_min = 2
                         p_max = info%nmaps
                      end if
                   else if (self%poltype(i)==3) then
                      p_min = p
                      p_max = p
                   else
                      write(*,*) '  Unknown poltype in component ',self%label,', parameter ',self%indlabel(i) 
                      stop
                   end if

                   do k = p_min,p_max
                      self%theta(i)%p%map(:,k) = tp%map(:,k)
                   end do
                end do
                call tp%dealloc()
             end if
          end if


          if (output_hdf) then
             call self%theta(i)%p%writeFITS(trim(dir)//'/'//trim(filename), &
                  & hdffile=chainfile, hdfpath=trim(path)//'/'//trim(adjustl(self%indlabel(i)))//&
                  & '_')
          else
             call self%theta(i)%p%writeFITS(trim(dir)//'/'//trim(filename))
          end if
          
          !write proposal length and number of proposals maps if local sampling was used
          if (any(self%lmax_ind_pol(:self%poltype(i),i)<0)) then
             filename = trim(self%label) // '_' // trim(self%indlabel(i)) // &
                  & '_proplen_'  // trim(postfix) // '.fits'
             call self%pol_proplen(i)%p%writeFITS(trim(dir)//'/'//trim(filename))

             filename = trim(self%label) // '_' // trim(self%indlabel(i)) // &
                  & '_nprop_'  // trim(postfix) // '.fits'
             call self%pol_nprop(i)%p%writeFITS(trim(dir)//'/'//trim(filename))

          end if

          !if pixelregions, create map without smoothed thetas (for input in new runs)
          if (any(self%pol_pixreg_type(1:self%poltype(i),i) > 0)) then
             
             info => comm_mapinfo(self%theta(i)%p%info%comm, self%theta(i)%p%info%nside, &
                  & self%theta(i)%p%info%lmax, self%theta(i)%p%info%nmaps, self%theta(i)%p%info%pol)
             tp => comm_map(info)
             tp%map = self%theta(i)%p%map
             do p = 1,self%poltype(i)
                if (self%pol_pixreg_type(p,i) /=3) cycle
                if (self%poltype(i) == 1) then
                   p_min=1
                   p_max=info%nmaps
                   if (only_pol) p_min = 2
                else if (self%poltype(i)==2) then
                   if (p == 1) then
                      p_min = 1
                      p_max = 1
                   else
                      p_min = 2
                      p_max = info%nmaps
                   end if
                else if (self%poltype(i)==3) then
                   p_min = p
                   p_max = p
                else
                   write(*,*) '  Unknown poltype in component ',self%label,', parameter ',self%indlabel(i) 
                   stop
                end if

                do j = 0,info%np-1
                   tp%map(j,p_min:p_max) = self%theta_pixreg(self%ind_pixreg_arr(j,p,i),p,i)
                end do
             end do
             filename = trim(self%label) // '_' // trim(self%indlabel(i)) // &
                  & '_noSmooth_'  // trim(postfix) // '.fits'
             call tp%writeFITS(trim(dir)//'/'//trim(filename))
             call tp%dealloc()

          end if

          !output theta, proposal length and number of proposals per pixel region to HDF
          if (output_hdf) then
             npol=min(self%nmaps,self%poltype(i))!only concerned about the maps/poltypes in use
             if (any(self%pol_pixreg_type(:npol,i) > 0)) then
                npr=0
                do j = 1,npol
                   if (self%npixreg(j,i)>npr) npr = self%npixreg(j,i)
                end do
                if (npr == 0) then !no pixelregions, theta = prior
                   if (self%theta(i)%p%info%myid == 0) write(*,*) 'No defined pixel regions for ',trim(self%label)//'_'//&
                        & trim(self%indlabel(i))//' in writing to HDF'
                else
                   allocate(dp_pixreg(npr,npol),int_pixreg(npr,npol))
                   !pixel region values for theta
                   dp_pixreg=self%theta_pixreg(1:npr,1:npol,i)
                   if (self%theta(i)%p%info%myid == 0) call write_hdf(chainfile, trim(path)//'/'//&
                        & trim(adjustl(self%indlabel(i)))//'_pixreg_val', real(dp_pixreg,sp))
                   !pixel region values for proposal length
                   dp_pixreg=self%proplen_pixreg(1:npr,1:npol,i)
                   if (self%theta(i)%p%info%myid == 0) call write_hdf(chainfile, trim(path)//'/'//&
                        & trim(adjustl(self%indlabel(i)))//'_pixreg_proplen', real(dp_pixreg,sp))
                   !pixel region values for number of proposals
                   int_pixreg=self%nprop_pixreg(1:npr,1:npol,i)
                   if (self%theta(i)%p%info%myid == 0) call write_hdf(chainfile, trim(path)//'/'//&
                        & trim(adjustl(self%indlabel(i)))//'_pixreg_nprop', int_pixreg)

                   deallocate(dp_pixreg,int_pixreg)
                end if
             end if

          end if

       end do
       !call update_status(status, "writeFITS_8")
       
       ! Write mixing matrices
       if (self%output_mixmat) then
          do i = 1, numband
             if (self%F_null(i,0)) cycle
             filename = 'mixmat_' // trim(self%label) // '_' // trim(data(i)%label) // '_' // &
                  & trim(postfix) // '.fits'
             call self%F(i,0)%p%writeFITS(trim(dir)//'/'//trim(filename))
          end do
       end if
       call update_status(status, "writeFITS_9")
    end if

  end subroutine dumpDiffuseToFITS

  ! Dump current sample to HEALPix FITS file
  subroutine initDiffuseHDF(self, cpar, hdffile, hdfpath)
    implicit none
    class(comm_diffuse_comp),  intent(inout) :: self
    type(comm_params),         intent(in)    :: cpar    
    type(hdf_file),            intent(in)    :: hdffile
    character(len=*),          intent(in)    :: hdfpath

    integer(i4b)       :: i, j, l, m, ierr
    integer(i4b)       :: p, p_min, p_max, npr, npol
    real(dp)           :: md(4)
    character(len=512) :: path
    class(comm_mapinfo), pointer :: info => null()
    class(comm_map), pointer     :: tp => null()
    real(dp),     allocatable, dimension(:,:) :: dp_pixreg
    integer(i4b), allocatable, dimension(:,:) :: int_pixreg
    
    if (trim(self%type) == 'md') then
       call read_hdf(hdffile, trim(adjustl(hdfpath))//'md/'//trim(adjustl(self%label)), md)
       do i = 0, self%x%info%nalm-1
          call self%x%info%i2lm(i,l,m)
          if (l == 0) then                 
             self%x%alm(i,1) =  md(1) * sqrt(4.d0*pi)      / self%RJ2unit_(1)  ! Monopole
          else if (l == 1 .and. m == -1) then   
             self%x%alm(i,1) =  md(3) * sqrt(4.d0*pi/3.d0) / self%RJ2unit_(1)  ! Y dipole
          else if (l == 1 .and. m ==  0) then   
             self%x%alm(i,1) =  md(4) * sqrt(4.d0*pi/3.d0) / self%RJ2unit_(1)  ! Z dipole
          else if (l == 1 .and. m ==  1) then   
             self%x%alm(i,1) = -md(2) * sqrt(4.d0*pi/3.d0) / self%RJ2unit_(1)  ! X dipole
          end if
       end do
    else
       path = trim(adjustl(hdfpath))//trim(adjustl(self%label))
       call self%Cl%initHDF(hdffile, path)
       call self%x%readHDF(hdffile, trim(adjustl(path))//'/amp_alm', .false.)
       !call self%x%readHDF(hdffile, trim(adjustl(path))//'/amp_map', .true.)    ! Read amplitudes
       do i = 1, self%x%info%nmaps
         self%x%alm(:,i) = self%x%alm(:,i) / (self%RJ2unit_(i) * self%cg_scale)
       end do

!!$       if (trim(self%label) == 'cmb') then
!!$          do i = 0, self%x%info%nalm-1
!!$             if (self%x%info%lm(1,i) == 0 .and. self%x%info%lm(2,i) == 0) then
!!$                write(*,*) 'Adding monopole'
!!$                self%x%alm(i,1) = self%x%alm(i,1) + 6 *sqrt(4*pi)
!!$             end if
!!$             if (self%x%info%lm(1,i) == 1 .and. self%x%info%lm(2,i) == 1) then
!!$                write(*,*) 'Adding x dipole'
!!$                self%x%alm(i,1) = self%x%alm(i,1) - 6.77d0 * sqrt(4.d0*pi/3.d0)
!!$             end if
!!$          end do
!!$       end if

!!$       if (trim(self%label) == 'ff') then
!!$          do i = 0, self%x%info%nalm-1
!!$             if (self%x%info%lm(1,i) == 0 .and. self%x%info%lm(2,i) == 0) then
!!$                write(*,*) 'Adding monopole'
!!$                self%x%alm(i,1) = self%x%alm(i,1) + 16 *sqrt(4*pi)
!!$             end if
!!$             if (self%x%info%lm(1,i) == 1 .and. self%x%info%lm(2,i) == 1) then
!!$                write(*,*) 'Adding x dipole'
!!$                self%x%alm(i,1) = self%x%alm(i,1) - 6.77d0 * sqrt(4.d0*pi/3.d0)
!!$             end if
!!$          end do
!!$       end if

!!$       if (trim(self%label) == 'ame') then
!!$          do i = 0, self%x%info%nalm-1
!!$             if (self%x%info%lm(1,i) == 0 .and. self%x%info%lm(2,i) == 0) then
!!$                write(*,*) 'Adding monopole'
!!$                self%x%alm(i,1) = self%x%alm(i,1) + 21.6 *sqrt(4*pi)
!!$             end if
!!$          end do
!!$       end if
!!$
!!$       if (trim(self%label) == 'ame2') then
!!$          do i = 0, self%x%info%nalm-1
!!$             if (self%x%info%lm(1,i) == 0 .and. self%x%info%lm(2,i) == 0) then
!!$                write(*,*) 'Adding monopole'
!!$                self%x%alm(i,1) = self%x%alm(i,1) + 11 *sqrt(4*pi)
!!$             end if
!!$          end do
!!$       end if
!!$
!!$       if (trim(self%label) == 'ame') then
!!$          do i = 0, self%x%info%nalm-1
!!$             if (self%x%info%lm(1,i) == 1 .and. self%x%info%lm(2,i) == 1) then
!!$                write(*,*) 'Adding synch dipole'
!!$                self%x%alm(i,1) = self%x%alm(i,1) + 25.d0
!!$             end if
!!$          end do
!!$       end if
!!$
!!$       if (trim(self%label) == 'ame2') then
!!$          do i = 0, self%x%info%nalm-1
!!$             if (self%x%info%lm(1,i) == 1 .and. self%x%info%lm(2,i) == 1) then
!!$                write(*,*) 'Adding synch dipole'
!!$                self%x%alm(i,1) = self%x%alm(i,1) + 14.d0
!!$             end if
!!$          end do
!!$       end if

       do i = 1, self%npar
          call self%theta(i)%p%readHDF(hdffile, trim(path)//'/'//trim(adjustl(self%indlabel(i)))//&
               & '_map', .true.)
          if (self%lmax_ind >= 0) then
             call self%theta(i)%p%readHDF(hdffile, trim(path)//'/'//trim(adjustl(self%indlabel(i)))//&
                  & '_alm', .false.)
             ! need to only take alms from polarizations with lmax > 0
             ! if any polarization is alm sampled. Only use alms to set polarizations with alm sampling
             if (all(self%lmax_ind_pol(1:self%poltype(i),i) >= 0)) then
                call self%theta(i)%p%Y_scalar
             else if (any(self%lmax_ind_pol(1:self%poltype(i),i) >= 0)) then
                info => comm_mapinfo(self%theta(i)%p%info%comm, self%theta(i)%p%info%nside, &
                     & self%theta(i)%p%info%lmax, self%theta(i)%p%info%nmaps, self%theta(i)%p%info%pol)
                tp => comm_map(info)
                tp%alm = self%theta(i)%p%alm
                call tp%Y_scalar
                do p = 1,self%poltype(i)
                   if (self%lmax_ind_pol(p,i) < 0) cycle
                   if (self%poltype(i) == 1) then
                      p_min=1
                      p_max=info%nmaps
                      if (only_pol) p_min = 2
                   else if (self%poltype(i)==2) then
                      if (p == 1) then
                         p_min = 1
                         p_max = 1
                      else
                         p_min = 2
                         p_max = info%nmaps
                      end if
                   else if (self%poltype(i)==3) then
                      p_min = p
                      p_max = p
                   else
                      write(*,*) '  Unknown poltype in component ',self%label,', parameter ',self%indlabel(i) 
                      stop
                   end if

                   do j = p_min,p_max
                      self%theta(i)%p%map(:,j) = tp%map(:,j)
                   end do
                end do
                call tp%dealloc()
             end if
          end if !lmax_ind > 0
          !if (trim(self%label) == 'dust' .and. i == 1) self%theta(i)%p%map(:,1) = 1.65d0 
          !if (trim(self%label) == 'dust' .and. i == 2) self%theta(i)%p%map(:,1) = 18.d0 
          !if (trim(self%label) == 'dust' .and. i > 1) self%theta(i)%p%map(:,1) = 1.6d0 
          !if (trim(self%label) == 'dust' .and. i == 1) self%theta(i)%p%alm(:,:) = 1.65d0 * sqrt(4*pi)
!!$          if (trim(self%label) == 'dust' .and. i == 2) self%theta(i)%p%alm(:,:) = 18 * sqrt(4*pi)
!!$          if (trim(self%label) == 'synch' .and. i == 1) self%theta(i)%p%alm(:,:) = -3.0d0 * sqrt(4*pi)
          !if (trim(self%label) == 'ame' .and. i == 1) self%theta(i)%p%alm(:,1) = self%theta(i)%p%alm(:,1) + 0.5d0*sqrt(4*pi)

          !Need to initialize pixelregions and local sampler from chain as well (where relevant)
          npol=min(self%nmaps,self%poltype(i))!only concerned about the maps/poltypes in use
          if (any(self%pol_pixreg_type(:npol,i) > 0)) then
             npr=0
             do j = 1,npol
                if (self%npixreg(j,i)>npr) npr = self%npixreg(j,i)
             end do
             if (npr == 0) then !no pixelregions, theta = prior
                if (self%theta(i)%p%info%myid == 0) write(*,*) 'No defined pixel regions for ',trim(self%label)//'_'//&
                     & trim(self%indlabel(i))
             else
                allocate(dp_pixreg(npr,npol),int_pixreg(npr,npol))
                !pixel region values for theta
                if (self%theta(i)%p%info%myid == 0) call read_hdf(hdffile, trim(path)//'/'//&
                     & trim(adjustl(self%indlabel(i)))//'_pixreg_val', dp_pixreg)
                call mpi_bcast(dp_pixreg, size(dp_pixreg),  MPI_DOUBLE_PRECISION, 0, self%theta(i)%p%info%comm, ierr)
                self%theta_pixreg(1:npr,1:npol,i)=dp_pixreg
                !pixel region values for proposal length
                if (self%theta(i)%p%info%myid == 0) call read_hdf(hdffile, trim(path)//'/'//&
                     & trim(adjustl(self%indlabel(i)))//'_pixreg_proplen', dp_pixreg)
                call mpi_bcast(dp_pixreg, size(dp_pixreg),  MPI_DOUBLE_PRECISION, 0, self%theta(i)%p%info%comm, ierr)
                self%proplen_pixreg(1:npr,1:npol,i)=dp_pixreg
                !pixel region values for number of proposals
                if (self%theta(i)%p%info%myid == 0) call read_hdf(hdffile, trim(path)//'/'//&
                     & trim(adjustl(self%indlabel(i)))//'_pixreg_nprop', int_pixreg)
                call mpi_bcast(int_pixreg, size(int_pixreg),  MPI_INTEGER, 0, self%theta(i)%p%info%comm, ierr)
                self%nprop_pixreg(1:npr,1:npol,i)=int_pixreg

                deallocate(dp_pixreg,int_pixreg)
             end if
          end if
       end do !i = 1,npar
    end if


!!$       if (trim(self%label)=='dust') then
!!$          call self%theta(2)%p%writeFITS("test.fits")
!!$          call mpi_finalize(i)
!!$          stop
!!$       end if

    !if (trim(self%label) == 'dust') write(*,*) 'range beta = ', minval(self%theta(1)%p%map), maxval(self%theta(1)%p%map)
    call self%updateMixmat

  end subroutine initDiffuseHDF
  
  subroutine add_to_npre(n, nside, lmax, nmaps)
    implicit none
    integer(i4b), intent(in) :: n, nside, lmax, nmaps
    npre      = npre + n
    nside_pre = min(nside_pre, nside)
    lmax_pre  = max(lmax_pre, lmax)
    nmaps_pre = max(nmaps_pre, nmaps)
  end subroutine add_to_npre

  ! Sample spectral parameters
  subroutine sampleDiffuseSpecInd(self, cpar, handle, id, iter)
    implicit none
    class(comm_diffuse_comp),                intent(inout)        :: self
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: id
    integer(i4b),                            intent(in)           :: iter

    integer(i4b) :: i, j, k, l, n, p, q, pix, ierr, ind(1), counter, n_ok
    integer(i4b) :: i_min, i_max, status, n_gibbs, n_pix, n_pix_tot, flag, npar, np, nmaps
    real(dp)     :: a, b, a_tot, b_tot, s, t1, t2, x_min, x_max, delta_lnL_threshold
    real(dp)     :: mu, sigma, par, w, mu_p, sigma_p, a_old, chisq, chisq_old, chisq_tot, unitconv
    real(dp)     :: x(1), theta_min, theta_max
    logical(lgt) :: ok
    logical(lgt), save :: first_call = .true.
    class(comm_comp), pointer :: c => null()
    real(dp),     allocatable, dimension(:)   :: lnL, P_tot, F, theta, a_curr
    real(dp),     allocatable, dimension(:,:) :: amp, buffer, alm_old
    !Following is for the local sampler
    real(dp)     :: old_theta, new_theta, mixing_old, mixing_new, lnL_new, lnL_old, res_lnL, delta_lnL, accept_rate
    integer(i4b) :: i_s, p_min, p_max, pixreg_nprop
    integer(i4b) :: n_spec_prop, n_accept, n_corr_prop, n_prop_limit, n_corr_limit, corr_len
    logical(lgt) :: first_sample
    class(comm_mapinfo),            pointer :: info => null()

    delta_lnL_threshold = 25.d0
    n                   = 101
    n_ok                = 50
    first_call          = .false.

    !c_lnL       => self
    id_lnL      = id
    c           => compList     ! Extremely ugly hack...
    do while (self%id /= c%id)
       c => c%next()
    end do
    select type (c)
    class is (comm_diffuse_comp)
       c_lnL => c
    end select

    info  => comm_mapinfo(c_lnL%x%info%comm, c_lnL%x%info%nside, &
         & c_lnL%x%info%lmax, c_lnL%x%info%nmaps, c_lnL%x%info%pol)


    npar      = c%npar
    np        = self%x_smooth%info%np
    nmaps     = self%x_smooth%info%nmaps

    theta_min = c_lnL%p_uni(1,id_lnL)
    theta_max = c_lnL%p_uni(2,id_lnL)

    !needs rewriting to only sample polarizations covered by polt_id
    if (trim(operation) == 'optimize') then
       allocate(buffer(0:self%x_smooth%info%np-1,self%x_smooth%info%nmaps))
       buffer = max(min(self%theta(id)%p%map,theta_max),theta_min)
       !do p = 1, nmaps
       do p = 1,self%poltype(id)  !I think this is the only change needed (now it only samples the poltype index polt_id)
          if (self%poltype(id) > 1 .and. only_pol .and. p == 1) cycle
          if (p > self%poltype(id)) cycle
          p_lnL = p
          !!$OMP PARALLEL DEFAULT(shared) PRIVATE(k,k_lnL,x,theta_lnL,a_lnL,i,ierr)
          allocate(theta_lnL(npar),a_lnL(self%nmaps))

          !!$OMP DO SCHEDULE(guided)
          do k = 0, np-1
             ! Perform non-linear search
             k_lnL     = k
             x(1)      = buffer(k,p)
             a_lnL     = self%x_smooth%map(k,:)
             do i = 1, npar
                if (i == id) cycle
                theta_lnL(i) = self%theta_smooth(i)%p%map(k,p)
             end do
             call powell(x, lnL_diffuse_multi, ierr)
             !write(*,*) k, x, ierr
             if (ierr == 0) then
                if (self%poltype(id) == 1) then
                   buffer(k,:) = x(1)
                else if (self%poltype(id) == 2) then
                   if (p == 1) then
                      buffer(k,1)   = x(1)
                   else
                      buffer(k,2:3) = x(1)
                   end if
                else if (self%poltype(id) == 3) then
                   buffer(k,p) = x(1)
                end if
             end if

          end do
          !!!$OMP END DO
          deallocate(theta_lnl, a_lnl)
          !!!$OMP END PARALLEL
       end do
       self%theta(id)%p%map = buffer

       ! Update mixing matrix
       call self%updateMixmat
       
       ! Ask for CG preconditioner update
       recompute_diffuse_precond = .true.

       deallocate(buffer)

   
    else if (trim(operation) == 'sample') then

       !allocate a buffer "map" with the same number of pixels as the processor is handling
       !if I understand correctly, each comm_map%info has info for the given processor. np = npix for that processor
       ! when a map is made, only the pixels that the processor is to handle is put in the map.
       ! that is why buffer goes to x_smooth%info%np and not to 12*nside**2
       ! Note: pix goes from 0 to np-1, so pix is not equal to the true pixel number. But, as all comm_maps of 
       ! the same nside have their pixels divided equally, then we still look up the correct pixels indep. of proc.
       
       allocate(buffer_lnL(0:self%x_smooth%info%np-1,self%x_smooth%info%nmaps))
       buffer_lnL=max(min(self%theta(id)%p%map,theta_max),theta_min) 

       do p = 1,self%poltype(id)
          if (self%lmax_ind_pol(p,id) >= 0) cycle !this set of polarizations are not to be local sampled (is checked before this point)
          if (self%poltype(id) > 1 .and. only_pol .and. p == 1) cycle !only polarization (poltype > 1)

          if (info%myid == 0) write(*,*) 'Sampling poltype index',p,'of ',self%poltype(id) !Needed?

          
          call wall_time(t1)
          p_lnL = p !global parameter for poltype
          if (self%pol_pixreg_type(p,id)==1) then
             call self%sampleDiffuseSpecIndFullsky(handle, id, p, iter)
          else if (self%pol_pixreg_type(p,id)==2) then
             call self%sampleDiffuseSpecIndSinglePix(handle, id, p, iter)
          else if (self%pol_pixreg_type(p,id)==3) then
             call self%sampleDiffuseSpecIndPixReg(cpar, handle, id, p, iter) ! not yet written, needs some thought on execution
          else
             write(*,*) 'Undefined spectral index sample region'
             write(*,*) 'Component:',trim(self%label),'ind:',trim(self%indlabel(id))
             stop
          end if
          call wall_time(t2)
          if (info%myid == 0) write(*,*) 'poltype:',self%poltype(id),' pol:', &
               & p,'CPU time specind = ', real(t2-t1,sp)



       end do

       !after sampling is done we assign the spectral index its new value(s)
       self%theta(id)%p%map = buffer_lnL

       ! Update mixing matrix
       call self%updateMixmat

       ! Ask for CG preconditioner update
       if (self%cg_unique_sampgroup > 0) recompute_diffuse_precond = .true.

       ! deallocate

       deallocate(buffer_lnL)



    end if !sample

  end subroutine sampleDiffuseSpecInd

  subroutine sampleDiffuseSpecIndSinglePix(self, handle, id, p, iter)
    implicit none
    class(comm_diffuse_comp),                intent(inout)        :: self
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: id
    integer(i4b),                            intent(in)           :: p       !incoming polarization
    integer(i4b),                            intent(in)           :: iter    !Gibbs iteration

    integer(i4b) :: i, j, k, l, n, q, pix, ierr, ind(1), counter, n_ok
    integer(i4b) :: i_min, i_max, status, n_gibbs, n_pix, n_pix_tot, flag, npar, np, nmaps
    real(dp)     :: a, b, a_tot, b_tot, s, t1, t2, x_min, x_max, delta_lnL_threshold
    real(dp)     :: mu, sigma, par, w, mu_p, sigma_p, a_old, chisq, chisq_old, chisq_tot, unitconv
    real(dp)     :: x(1), theta_min, theta_max
    logical(lgt) :: ok
    logical(lgt), save :: first_call = .true.
    class(comm_comp), pointer :: c => null()
    real(dp),     allocatable, dimension(:)   :: lnL, P_tot, F, theta, a_curr
    real(dp),     allocatable, dimension(:,:) :: amp, buffer, alm_old
    !Following is for the local sampler
    real(dp)     :: old_theta, new_theta, mixing_old, mixing_new, lnL_new, lnL_old, res_lnL, delta_lnL
    real(dp)     :: accept_rate, accept_scale, lnL_sum, sum_proplen, sum_accept
    integer(i4b) :: i_s, p_min, p_max, pixreg_nprop, band_count, burn_in, count_accept
    integer(i4b) :: n_spec_prop, n_accept, n_corr_prop, n_prop_limit, n_corr_limit, corr_len
    logical(lgt) :: first_sample, use_det, burned_in
    real(dp),     allocatable, dimension(:) :: all_thetas, data_arr, invN_arr, mixing_old_arr, mixing_new_arr
    real(dp),     allocatable, dimension(:) :: theta_corr_arr
    integer(i4b), allocatable, dimension(:) :: band_i, pol_j
    class(comm_mapinfo),            pointer :: info => null()


    info  => comm_mapinfo(c_lnL%x%info%comm, c_lnL%x%info%nside, &
         & c_lnL%x%info%lmax, c_lnL%x%info%nmaps, c_lnL%x%info%pol)
                
                
    delta_lnL_threshold = 25.d0
    n                   = 101
    n_ok                = 50
    burn_in             = 20
    burned_in           = .false.
    first_call          = .false.
    
    n_prop_limit = 40
    n_corr_limit = 40
    count_accept = 0
    sum_accept   = 0.d0
    sum_proplen   = 0.d0
    !p_lnL = p

    ! we have allready got globally from sampleDiffuseSpecInd: c_lnL, p_lnL (also in input), id_lnL (also in input),
    !     buffer_lnL (theta, limited by min/max)

    if (c_lnL%poltype(id) == 1) then
       p_min = 1; p_max = c_lnL%nmaps
       if (only_pol) p_min = 2
    else if (c_lnL%poltype(id) == 2) then
       if (p == 1) then
          p_min = 1; p_max = 1
       else
          p_min = 2; p_max = c_lnL%nmaps
       end if
    else if (c_lnL%poltype(id) == 3) then
       p_min = p
       p_max = p
    else
       write(*,*) 'Unsupported polarization type'
       stop
    end if

    npar      = c_lnL%npar
    np        = self%x_smooth%info%np
    nmaps     = self%x_smooth%info%nmaps

    theta_min = c_lnL%p_uni(1,id)
    theta_max = c_lnL%p_uni(2,id)

    !set up which bands and polarizations to include
    !this is to exclude the need of checking this for each time the mixing matrix needs to be computed per pixel
    allocate(band_i(1000),pol_j(1000))
    band_count=0
    do k = 1,numband !run over all active bands
       !check if the band is associated with the component, i.e. frequencies are within 
       !freq. limits for the component
       if (.not. associated(rms_smooth(k)%p)) cycle
       if (data(k)%bp(0)%p%nu_c < self%nu_min_ind(id) .or. data(k)%bp(0)%p%nu_c > self%nu_max_ind(id)) cycle

       do l = p_min,p_max !loop all maps of band k associated with p (given poltype)
          if (l > data(k)%info%nmaps) cycle !no map in band k for polarization l
          band_count = band_count + 1
          band_i(band_count)=k
          pol_j(band_count)=l
       end do
    end do

    if (band_count==0) then
       buffer_lnL(:,p_min:p_max)=c_lnL%p_gauss(1,id) !set all thata to prior, as no bands are valid, no data
       deallocate(band_i,pol_j)
       if (info%myid == 0)  write(*,*) 'no data bands available for sampling of spec ind'
       return
    else
       if (info%myid == 0 .and. .true.) then
             write(*,fmt='(a)') '  Active bands'
          do k = 1,band_count
             write(*,fmt='(a,i1)') '   band: '//trim(data(band_i(k))%label)//', -- polarization: ',pol_j(k)
          end do
       end if
    end if

    allocate(all_thetas(npar))
    allocate(mixing_old_arr(band_count),mixing_new_arr(band_count),invN_arr(band_count),data_arr(band_count))
    do pix = 0,np-1 !loop over the number of pixels the processor is handling
       
       if (self%pol_ind_mask(id)%p%map(pix,p) < 0.5d0) then     ! if pixel is masked out
          buffer_lnL(pix,p_min:p_max) = c_lnL%p_gauss(1,id) ! set spec. ind to prior
          self%pol_nprop(id)%p%map(pix,p)=0.d0 !to mark that the pixel is masked out, we set nprop to zero
          cycle
       end if
       
       if (self%pol_sample_nprop(p,id) .or. self%pol_sample_proplen(p,id)) then
          pixreg_nprop = 10000 !should be enough to find correlation length (and proposal length if prompted) 
          allocate(theta_corr_arr(pixreg_nprop))
          n_spec_prop = 0
          n_corr_prop = 0 
          !self%pol_sample_nprop(p,id) = boolean array of size (poltype,npar)
          !self%pol_sample_proplen(p,id) = boolean array of size (poltype,npar)
       else
          pixreg_nprop = IDINT(self%pol_nprop(id)%p%map(pix,p)) ! comm_map is real(dp), use IDINT to get the integer value
       end if

       n_accept = 0
       first_sample = .true.

       do j = 1,pixreg_nprop !propose new sample n times (for faster mixing in total)
          
          if (first_sample) then
             old_theta = buffer_lnL(pix,p)

             !get the values of the remaining spec inds of the component
             do i = 1, npar
                if (i == id) cycle
                all_thetas(i) = self%theta_smooth(i)%p%map(pix,p)
             end do
          end if

          !draw new spec ind
          new_theta = old_theta + self%pol_proplen(id)%p%map(pix,p)*rand_gauss(handle)

          !init log-like for new sample
          if (first_sample) lnL_old = 0.d0
          lnL_new = 0.d0
          
          if (new_theta > theta_max .or. new_theta < theta_min) then !new proposal is outside limits
             lnL_new = -1.d30 
             ! skip the true caclulation of lnL, we reject the sample ~100%
          else
             ! do a clever way of computing the lnL for the different evaluation types
             if (trim(self%pol_lnLtype(p,id))=='chisq') then
                !write the chisq sampling here

                do k = 1,band_count !run over all active bands
                   ! get conversion factor from amplitude to data (i.e. mixing matrix element)           
                   ! both for old and new spec. ind.
                   if (first_sample) then
                      all_thetas(id)=old_theta
                      mixing_old = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                           & data(band_i(k))%gain * c_lnL%cg_scale
                   end if
                   all_thetas(id)=new_theta
                   mixing_new = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                        & data(band_i(k))%gain * c_lnL%cg_scale

                   !compute chisq for pixel 'pix' of the associated band
                   if (first_sample) then
                      res_lnL = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_old*self%x_smooth%map(pix,pol_j(k))
                      lnL_old = lnL_old -0.5d0*(res_lnL*rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k)))**2
                   end if
                   res_lnL = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_new*self%x_smooth%map(pix,pol_j(k))
                   lnL_new = lnL_new -0.5d0*(res_lnL*rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k)))**2
                end do

             else if ((trim(self%pol_lnLtype(p,id))=='ridge') .or. &
                  & (trim(self%pol_lnLtype(p,id))=='marginal')) then
                if (trim(self%pol_lnLtype(p,id))=='ridge') then
                   use_det = .false.
                else
                   use_det = .true.
                end if

                !! build mixing matrix, invN and data for pixel
                do k = 1,band_count !run over all active bands
                   ! get conversion factor from amplitude to data (i.e. mixing matrix element)           
                   ! both for old and new spec. ind.
                   if (first_sample) then
                      all_thetas(id)=old_theta
                      mixing_old_arr(k) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                           & data(band_i(k))%gain * c_lnL%cg_scale
                   end if
                   all_thetas(id)=new_theta
                   mixing_new_arr(k) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                        & data(band_i(k))%gain * c_lnL%cg_scale

                   data_arr(k)=res_smooth(band_i(k))%p%map(pix,pol_j(k))
                   invN_arr(k)=rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k))**2 !assumed diagonal and uncorrelated 
                end do

                !compute the marginal/ridge log-likelihood for the pixel
                if (first_sample) then
                   lnL_old = lnL_old + comp_lnL_marginal_diagonal(mixing_old_arr, invN_arr, data_arr, &
                        & use_det, band_count)
                end if
                lnL_new = lnL_new + comp_lnL_marginal_diagonal(mixing_new_arr, invN_arr, data_arr, &
                     & use_det, band_count)

             else 
                write(*,*) 'invalid polarized lnL sampler type'
                stop

             end if !chisq/marginal/ridge

             !add spec. ind. priors
             if (first_sample) then
                all_thetas(id)=old_theta
                do l = 1, npar
                   if (c_lnL%p_gauss(2,l) > 0.d0) then
                      lnL_old = lnL_old - 0.5d0 * (all_thetas(l)-c_lnL%p_gauss(1,l))**2 / c_lnL%p_gauss(2,l)**2 
                   end if
                   all_thetas(id)=new_theta
                end do
             end if
             do l = 1, npar
                if (c_lnL%p_gauss(2,l) > 0.d0) then
                   lnL_new = lnL_new - 0.5d0 * (all_thetas(l)-c_lnL%p_gauss(1,l))**2 / c_lnL%p_gauss(2,l)**2 
                end if
             end do
             !first sample done
             first_sample = .false.

          end if !outside spec limits

          !accept/reject new spec ind
          delta_lnL = lnL_new-lnL_old

          if (delta_lnL < 0.d0) then
             !if delta_lnL is more negative than -25.d0, limit to -25.d0
             if (abs(delta_lnL) > delta_lnL_threshold) delta_lnL = -delta_lnL_threshold 

             !draw random uniform number
             a = rand_uni(handle) !draw uniform number from 0 to 1
             if (exp(delta_lnL) > a) then
                !accept
                old_theta = new_theta
                lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
                n_accept = n_accept + 1
             end if
          else
             !accept new sample, higher likelihood
             old_theta = new_theta
             lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
             n_accept = n_accept + 1
          end if
          
          if (j > burn_in) burned_in = .true.
          ! evaluate proposal_length, then correlation length. 
          !This should only be done the first gibbs iteration, if prompted to do so from parameter file.
          if (self%pol_sample_proplen(p,id)) then
             n_spec_prop = n_spec_prop + 1
             accept_rate = n_accept*1.d0/n_spec_prop
             if (.not. burned_in) then 
                !giving som time to let the first steps "burn in" if one starts far away from max lnL
                n_spec_prop = 0
                n_accept = 0
             else if (n_spec_prop > n_prop_limit) then
                if (accept_rate > 0.7d0) then
                   accept_scale = 1.5d0
                else if (accept_rate < 0.3d0) then
                   accept_scale = 0.5d0
                else 
                   accept_scale = -1.d0
                end if

                if (accept_scale > 0.d0) then
                   self%pol_proplen(id)%p%map(pix,p) = self%pol_proplen(id)%p%map(pix,p)* &
                        & accept_scale !set new proposal length
                   n_accept = 0 !reset with new prop length
                   n_spec_prop = 0 !reset with new prop length
                else
                   sum_accept = sum_accept + accept_rate
                   count_accept = count_accept + 1
                   sum_proplen = sum_proplen + self%pol_proplen(id)%p%map(pix,p)
                   if (.not. self%pol_sample_nprop(p,id)) then
                      exit 
                      !ugly, but we only end up here in the first gibbs sample if we are only fitting proposal
                      !lengths and not correlation length
                   end if
                end if
             end if

          else if (self%pol_sample_nprop(p,id)) then
             ! if we are keeping track of spec inds for corr. length, then save current spec ind
             n_corr_prop = n_corr_prop + 1
             theta_corr_arr(n_corr_prop) = old_theta
             if (.not. burned_in) then 
                !giving som time to let the first steps "burn in" if one starts far away from max lnL
                n_corr_prop = 0
             else if (n_corr_prop > n_corr_limit) then
                corr_len = calc_corr_len(theta_corr_arr(1:n_corr_prop),n_corr_prop) !returns -1 if no corr.len. is found
                if (corr_len > 0) then ! .or. n_corr_prop > 4*self%nprop_uni(2,id)) then
                   !set nprop (for pix region) to x*corr_len, x=2, but inside limits from parameter file
                   self%pol_nprop(id)%p%map(pix,p) = 1.d0*min(max(2*corr_len,self%nprop_uni(1,id)),self%nprop_uni(2,id))
                   exit 
                   !ugly, but we only end up here in the first gibbs sample if we are fitting correlation length

                end if
             end if
          else
             accept_rate = n_accept*1.d0/j
             sum_accept = sum_accept + accept_rate
             count_accept = count_accept + 1
             sum_proplen = sum_proplen + self%pol_proplen(id)%p%map(pix,p)
          end if
          
       end do !pixreg_nprop


       if (allocated(theta_corr_arr)) deallocate(theta_corr_arr)
       
       !assign last accepted theta value (old_theta) to the pixels of the current pixel region, 
       !in the polarizations given by poltype and p
       buffer_lnL(pix,p_min:p_max) = old_theta
       
    end do !pix = 1,np


    call mpi_allreduce(MPI_IN_PLACE, sum_accept, 1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
    call mpi_allreduce(MPI_IN_PLACE, count_accept, 1, MPI_INTEGER, MPI_SUM, info%comm, ierr)
    call mpi_allreduce(MPI_IN_PLACE, sum_proplen, 1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
    if (info%myid==0) then
       write(*,fmt='(a,f6.4)') '   avg. accept rate     = ', sum_accept/count_accept
       write(*,fmt='(a,e14.4)') '   avg. proposal length = ', sum_proplen/count_accept
    end if
    if (iter >= 2) then !only stop sampling after 2nd iteration, in case first iter is far off.
       self%pol_sample_nprop(p,id) = .false. !set sampling of correlation length (number of proposals and
       self%pol_sample_proplen(p,id) = .false. !proposal length to false (this is only to be done first gibbs iteration)
    end if
    !deallocate
    deallocate(mixing_old_arr,mixing_new_arr,data_arr,invN_arr,all_thetas)
    deallocate(band_i,pol_j)


  end subroutine sampleDiffuseSpecIndSinglePix


  subroutine sampleDiffuseSpecIndFullsky(self, handle, id, p, iter)
    implicit none
    class(comm_diffuse_comp),                intent(inout)        :: self
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: id
    integer(i4b),                            intent(in)           :: p       !incoming polarization
    integer(i4b),                            intent(in)           :: iter    !Gibbs iteration

    integer(i4b) :: i, j, k, l, n, q, pix, ierr, ind(1), counter, n_ok
    integer(i4b) :: i_min, i_max, status, n_gibbs, n_pix, n_pix_tot, flag, npar, np, nmaps
    real(dp)     :: a, b, a_tot, b_tot, s, t1, t2, x_min, x_max, delta_lnL_threshold
    real(dp)     :: mu, sigma, par, w, mu_p, sigma_p, a_old, chisq, chisq_old, chisq_tot, unitconv
    real(dp)     :: buff1_r(1), buff2_r(1), theta_min, theta_max
    logical(lgt) :: ok
    logical(lgt), save :: first_call = .true.
    class(comm_comp), pointer :: c => null()
    !Following is for the local sampler
    real(dp)     :: old_theta, new_theta, mixing_old, mixing_new, lnL_new, lnL_old, res_lnL, delta_lnL
    real(dp)     :: accept_rate, accept_scale, lnL_sum, init_theta
    integer(i4b) :: i_s, p_min, p_max, pixreg_nprop, band_count, pix_count, buff1_i(1), buff2_i(1), burn_in
    integer(i4b) :: n_spec_prop, n_accept, n_corr_prop, n_prop_limit, n_corr_limit, corr_len, out_every
    logical(lgt) :: first_sample, loop_exit, use_det, burned_in
    real(dp),     allocatable, dimension(:)   :: all_thetas, data_arr, invN_arr, mixing_old_arr, mixing_new_arr
    real(dp),     allocatable, dimension(:)   :: theta_corr_arr
    integer(i4b), allocatable, dimension(:)   :: band_i, pol_j
    class(comm_mapinfo),            pointer   :: info => null()


    info  => comm_mapinfo(c_lnL%x%info%comm, c_lnL%x%info%nside, &
         & c_lnL%x%info%lmax, c_lnL%x%info%nmaps, c_lnL%x%info%pol)
                

    delta_lnL_threshold = 25.d0
    n                   = 101
    n_ok                = 50
    burn_in             = 20
    burned_in           = .false.
    first_call          = .false.
    
    n_prop_limit = 20
    n_corr_limit = 20
    out_every    = 10
    !p_lnL = p

    ! we have allready got globally from sampleDiffuseSpecInd: c_lnL, p_lnL (also in input), id_lnL (also in input),
    !     buffer_lnL (theta, limited by min/max)

    if (c_lnL%poltype(id) == 1) then
       p_min = 1; p_max = c_lnL%nmaps
       if (only_pol) p_min = 2
    else if (c_lnL%poltype(id) == 2) then
       if (p == 1) then
          p_min = 1; p_max = 1
       else
          p_min = 2; p_max = c_lnL%nmaps
       end if
    else if (c_lnL%poltype(id) == 3) then
       p_min = p
       p_max = p
    else
       write(*,*) 'Unsupported polarization type'
       stop
    end if

    npar      = c_lnL%npar
    np        = self%x_smooth%info%np
    nmaps     = self%x_smooth%info%nmaps

    theta_min = c_lnL%p_uni(1,id)
    theta_max = c_lnL%p_uni(2,id)

    !set up which bands and polarizations to include

    allocate(band_i(3*numband),pol_j(3*numband))

    band_count=0
    do k = 1,numband !run over all active bands
       !check if the band is associated with the smoothed component, and if band frequencies are within 
       !freq. limits for the component
       if (.not. associated(rms_smooth(k)%p)) cycle
       if (data(k)%bp(0)%p%nu_c < self%nu_min_ind(id) .or. data(k)%bp(0)%p%nu_c > self%nu_max_ind(id)) cycle
                   
       do l = p_min,p_max !loop all maps of band k associated with p (given poltype)
          if (l > data(k)%info%nmaps) cycle !no map in band k for polarization l
          band_count = band_count + 1
          band_i(band_count)=k
          pol_j(band_count)=l
       end do
    end do

    if (band_count==0) then
       buffer_lnL(:,p_min:p_max)=c_lnL%p_gauss(1,id) !set theta to prior, as no bands are valid, no data
       deallocate(band_i,pol_j)
       if (info%myid == 0)  write(*,*) 'no data bands available for sampling of spec ind'
       return
    else
       if (info%myid == 0 .and. .true.) then
             write(*,fmt='(a)') '  Active bands'
          do k = 1,band_count
             write(*,fmt='(a,i1)') '   band: '//trim(data(band_i(k))%label)//', -- polarization: ',pol_j(k)
          end do
       end if
    end if

    allocate(all_thetas(npar))
    allocate(mixing_old_arr(band_count),mixing_new_arr(band_count),invN_arr(band_count),data_arr(band_count))

    n_spec_prop = 0
    n_accept = 0
    n_corr_prop = 0 
    if (self%pol_sample_nprop(p,id) .or. self%pol_sample_proplen(p,id)) then
       pixreg_nprop = 10000 !should be enough to find correlation length (and proposal length if prompted) 
       allocate(theta_corr_arr(pixreg_nprop))
       !self%pol_sample_nprop(j,p,id) = boolean array of size (n_pixreg,poltype)
       !self%pol_sample_proplen(j,p,id) = boolean array of size (n_pixreg,poltype)
    else
       pixreg_nprop = IDINT(self%pol_nprop(id)%p%map(0,p)) ! comm_map is real(dp), use IDINT() to convert to integer
    end if

    first_sample  = .true.
    loop_exit     = .false.
    
    call wall_time(t1)
    do j = 1,pixreg_nprop !propose new sample n times (for faster mixing in total)
       if (info%myid == 0) then
          if (first_sample) then
             old_theta = buffer_lnL(0,p_min) !taking the first pixel value (should be full sky) 
                ! that the root processor operates on
             init_theta = old_theta
             write(*,fmt='(a, f10.5)') "  initial sepc. ind. value: ", init_theta
          end if

          !draw new spec ind, 
          new_theta = old_theta + self%pol_proplen(id)%p%map(0,p)*rand_gauss(handle)

       end if


       !broadcast new_theta, and new old_theta, and all_thetas to the other processors
       if (first_sample) call mpi_bcast(old_theta, 1, MPI_DOUBLE_PRECISION, &
            & 0, c_lnL%comm, ierr)
       call mpi_bcast(new_theta, 1, MPI_DOUBLE_PRECISION, 0, c_lnL%comm, ierr)

       if (info%myid==0 .and. mod(j,out_every)==0) then
          write(*,fmt='(a, i6, a, f10.5, a, f10.5)') "  proposal: ", j," -- Current ind: ", &
               & old_theta, " -- proposed ind: ", new_theta
       end if


       !init log-like for new sample
       if (first_sample) lnL_old = 0.d0
       lnL_new = 0.d0       

       if (new_theta > theta_max .or. new_theta < theta_min) then !new proposal is outside limits
          lnL_new = -1.d30 
          ! skip the true caclulation of lnL, we reject the sample ~100%
          if (info%myid==0) write(*,fmt='(a, f10.5, a, f10.5)') "    Proposed ind outside limits.  min: ", &
               & theta_min," -- max: ", theta_max
       else
          pix_count = 0
          all_thetas(id)=new_theta
          !lnL type should split here
          if (trim(self%pol_lnLtype(p,id))=='chisq') then
             !loop bands and then pixels

             do pix = 0,np-1 !loop over pixels covered by the processor
                if (self%pol_ind_mask(id)%p%map(pix,p) < 0.5d0) cycle     ! if pixel is masked out

                !get the values of the remaining spec inds of the component for the given pixel
                do i = 1, npar
                   if (i == id) cycle
                   all_thetas(i) = self%theta_smooth(i)%p%map(pix,p) 
                end do

                do k = 1,band_count !run over all active bands
                   pix_count = pix_count + 1

                   ! get conversion factor from amplitude to data (i.e. mixing matrix element)           
                   ! both for old and new spec. ind. and calc. log likelihood (chisq)
                   if (first_sample) then
                      all_thetas(id)=old_theta
                      mixing_old = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                           & data(band_i(k))%gain * c_lnL%cg_scale
                      all_thetas(id)=new_theta

                      res_lnL = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_old* &
                           & self%x_smooth%map(pix,pol_j(k))
                      lnL_old = lnL_old -0.5d0*(res_lnL*rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k)))**2

                   end if
                   mixing_new = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                        & data(band_i(k))%gain * c_lnL%cg_scale

                   res_lnL = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_new* &
                        & self%x_smooth%map(pix,pol_j(k))
                   lnL_new = lnL_new -0.5d0*(res_lnL*rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k)))**2
                end do
             end do

          else if ((trim(self%pol_lnLtype(p,id))=='ridge') .or. &
               & (trim(self%pol_lnLtype(p,id))=='marginal')) then
             if (trim(self%pol_lnLtype(p,id))=='ridge') then
                use_det = .false.
             else
                use_det = .true.
             end if

             do pix = 0,np-1 !More precise, we need to loop over pixels covered by the processor
                
                if (self%pol_ind_mask(id)%p%map(pix,p) < 0.5d0) cycle     ! if pixel is masked out

                !get the values of the remaining spec inds of the component for the given pixel
                do i = 1, npar
                   if (i == id) cycle
                   all_thetas(i) = self%theta_smooth(i)%p%map(pix,p) 
                end do

                !build mixing matrix
                do k = 1,band_count !run over all active bands
                   pix_count = pix_count + 1

                   if (first_sample) then
                      all_thetas(id)=old_theta
                      mixing_old_arr(k) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                           & data(band_i(k))%gain * c_lnL%cg_scale
                      all_thetas(id)=new_theta
                   end if
                   mixing_new_arr(k) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                        & data(band_i(k))%gain * c_lnL%cg_scale
                   pix_count = pix_count + 1
                   data_arr(k)=res_smooth(band_i(k))%p%map(pix,pol_j(k))
                   invN_arr(k)=rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k))**2 !assumed diagonal and uncorrelated 
                end do

                !compute the marginal/ridge log-like for the pixel
                if (first_sample) then
                   lnL_old = lnL_old + comp_lnL_marginal_diagonal(mixing_old_arr, invN_arr, data_arr, &
                        & use_det, band_count)
                end if
                lnL_new = lnL_new + comp_lnL_marginal_diagonal(mixing_new_arr, invN_arr, data_arr, &
                     & use_det, band_count)
             end do
          else 
             if (info%myid==0) write(*,*) 'invalid polarized lnL sampler type'
             stop

          end if !chisq/marginal/ridge


          !Sum lnL_new (and lnL_old) among all processors
          if (first_sample) then
             !there is a problem with the mpi_allreduce call. Commander3 crashes
             call mpi_allreduce(lnL_old, buff2_r, 1, MPI_DOUBLE_PRECISION, & 
                  & MPI_SUM, info%comm, ierr)
             lnL_old = buff2_r(1)

             call mpi_allreduce(pix_count, k, 1, MPI_INTEGER, & 
                  & MPI_SUM, info%comm, ierr)
             pix_count = k

             if (pix_count == 0) then! there are no data that is valid
                old_theta = c_lnL%p_gauss(1,id) !set theta to prior
                !set nprop to 0, write error message, then exit sampling
                self%pol_nprop(id)%p%map(:,p) = 0.d0
                if (info%myid==0) write(*,*) "    No valid pixels to sample spectral index from"
                exit
             end if

          end if

          call mpi_allreduce(lnL_new, lnl_sum, 1, MPI_DOUBLE_PRECISION, & 
               & MPI_SUM, info%comm, ierr)
          lnL_new = lnL_sum

          !add spec. ind. priors
          if (info%myid == 0) then
             if (first_sample) then
                all_thetas(id)=old_theta
                do l = 1, npar
                   if (c_lnL%p_gauss(2,l) > 0.d0) then
                      lnL_old = lnL_old - 0.5d0 * (all_thetas(l)-c_lnL%p_gauss(1,l))**2 / c_lnL%p_gauss(2,l)**2 
                   end if
                   all_thetas(id)=new_theta
                end do
             end if
             do l = 1, npar
                if (c_lnL%p_gauss(2,l) > 0.d0) then
                   lnL_new = lnL_new - 0.5d0 * (all_thetas(l)-c_lnL%p_gauss(1,l))**2 / c_lnL%p_gauss(2,l)**2 
                end if
             end do
          end if

          !first sample done
          first_sample = .false.

          if (pix_count==0) then
             !set theta to prior
             old_theta = c_lnL%p_gauss(1,id)
             !set number of proposals to zero, i.e. spectral index will not be sampled in later iterations 
             !as there are no unmasked data
             self%pol_nprop(id)%p%map(:,p)=0.d0
             exit   !exit sampling
          end if

       end if !new_theta outside spec limits

       if (info%myid == 0) then
          !accept/reject new spec ind
          delta_lnL = lnL_new-lnL_old

          if (mod(j,out_every)==0) write(*,fmt='(a, e14.5)') "    lnL_new - lnL_old = ", delta_lnL


          if (delta_lnL < 0.d0) then
             !if delta_lnL is more negative than -25.d0, limit to -25.d0
             if (abs(delta_lnL) > delta_lnL_threshold) delta_lnL = -delta_lnL_threshold 

             !draw random uniform number
             a = rand_uni(handle) !draw uniform number from 0 to 1
             if (exp(delta_lnL) > a) then
                !accept
                old_theta = new_theta
                lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
                n_accept = n_accept + 1
             end if
          else
             !accept new sample, higher likelihood
             old_theta = new_theta
             lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
             n_accept = n_accept + 1
          end if

          if (j > burn_in) burned_in = .true.
          ! evaluate proposal_length, then correlation length. 
          !This should only be done the first gibbs iteration, if prompted to do so from parameter file.
          if (self%pol_sample_proplen(p,id)) then
             n_spec_prop = n_spec_prop + 1
             accept_rate = n_accept*1.d0/n_spec_prop
             if (mod(n_spec_prop,out_every)==0) write(*,fmt='(a, f6.4)') "   accept rate = ", accept_rate
             if (.not. burned_in) then 
                !giving som time to let the first steps "burn in" if one starts far away from max lnL
                n_spec_prop = 0
                n_accept = 0
             else if (n_spec_prop > n_prop_limit) then
                if (accept_rate > 0.7d0) then
                   accept_scale = 1.5d0
                else if (accept_rate < 0.3d0) then
                   accept_scale = 0.5d0
                else 
                   accept_scale = -1.d0
                end if

                if (accept_scale > 0.d0) then
                   !prop_len only ever used by root processor
                   self%pol_proplen(id)%p%map(0,p) = self%pol_proplen(id)%p%map(0,p)*accept_scale !set new proposal length
                   n_accept = 0 !reset with new prop length
                   n_spec_prop = 0 !reset with new prop length
                   write(*,fmt='(a, e14.5)') "     New prop. len. = ", self%pol_proplen(id)%p%map(0,p)
                else
                   self%pol_sample_proplen(p,id)=.false. !making sure we dont go back into prop. len. sampling
                   !only necessary to update for myid==0
                   if (.not. self%pol_sample_nprop(p,id)) then
                      loop_exit=.true. !we have found the ideal proposal length, not looking for correlation length
                   end if
                end if
             end if

          else if (self%pol_sample_nprop(p,id)) then
             ! if we are keeping track of spec inds for corr. length, then save current spec ind
             n_corr_prop = n_corr_prop + 1
             theta_corr_arr(n_corr_prop) = old_theta
             if (.not. burned_in) then
                n_corr_prop = 0
             else if (n_corr_prop > n_corr_limit) then
                corr_len = calc_corr_len(theta_corr_arr(1:n_corr_prop),n_corr_prop) !returns -1 if no corr.len. is found
                if (corr_len > 0) then ! .or. n_corr_prop > 4*self%nprop_uni(2,id)) then
                   !set pixreg_nprop (for pix region) to x*corr_len, x=2, but inside limits from parameter file
                   self%pol_nprop(id)%p%map(:,p) = 1.d0*min(max(2*corr_len,self%nprop_uni(1,id)),self%nprop_uni(2,id))
                   loop_exit=.true. !we have found the correlation length
                   self%pol_sample_nprop(p,id)=.false. !making sure we dont go back into corr. len. sampling
                end if
             end if
          else
             accept_rate = n_accept*1.d0/j
             if (mod(j,out_every)==0) write(*,fmt='(a, f6.4)') "   accept rate = ", accept_rate
          end if
          
       end if

       call wall_time(t2)

       if (info%myid == 0 .and. mod(j,out_every)==0) write(*,*) '   Sample:',j,'  Wall time per sample:',real((t2-t1)/j,sp)

       call mpi_bcast(loop_exit, 1, MPI_LOGICAL, 0, c_lnL%comm, ierr)

       if (loop_exit) exit
       
    end do !nprop

    !bcast proposal length
    call mpi_bcast(self%pol_proplen(id)%p%map(0,p), 1, MPI_DOUBLE_PRECISION, 0, &
         & c_lnL%comm, ierr)
    self%pol_proplen(id)%p%map(:,p) = self%pol_proplen(id)%p%map(0,p)

    !bcast number of proposals
    call mpi_bcast(self%pol_nprop(id)%p%map(0,p), 1, MPI_DOUBLE_PRECISION, 0, &
         & c_lnL%comm, ierr)
    self%pol_nprop(id)%p%map(:,p) = self%pol_nprop(id)%p%map(0,p)

    !bcast last valid theta to all procs to update theta map
    call mpi_bcast(old_theta, 1, MPI_DOUBLE_PRECISION, 0, c_lnL%comm, ierr)

    buffer_lnL(:,p_min:p_max) = old_theta

    if (info%myid==0) then
       write(*,fmt='(a, f10.5)') "    Final spec. ind. value: ", old_theta
       write(*,fmt='(a, i5, a, e14.5)') "    Number of proposals: ",IDINT(self%pol_nprop(id)%p%map(0,p)), &
            & "  --  Proposal length: ", self%pol_proplen(id)%p%map(0,p)
       write(*,fmt='(a, e14.5)') "    Difference in spec. ind., new - old: ", old_theta-init_theta
       write(*,*) ''
    end if

    if (iter >= 2) then !only stop sampling after 2nd iteration, in case first iter is far off.
       self%pol_sample_nprop(p,id) = .false. !set sampling of correlation length (number of proposals) and
       self%pol_sample_proplen(p,id) = .false. !proposal length to false (this is only to be done in the first 
                                               !gibbs iteration anyways)
    end if

    !deallocate
    deallocate(mixing_old_arr,mixing_new_arr,data_arr,invN_arr,all_thetas)
    deallocate(band_i,pol_j)
    if (allocated(theta_corr_arr)) deallocate(theta_corr_arr)

  end subroutine sampleDiffuseSpecIndFullsky


  subroutine sampleDiffuseSpecIndPixReg(self, cpar, handle, id, p, iter)
    implicit none
    class(comm_diffuse_comp),                intent(inout)        :: self
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: id
    integer(i4b),                            intent(in)           :: p       !incoming polarization
    integer(i4b),                            intent(in)           :: iter    !Gibbs iteration

    integer(i4b) :: i, j, k, l, n, q, pix, ierr, ind(1), counter, n_ok
    integer(i4b) :: i_min, i_max, status, n_gibbs, n_pix, n_pix_tot, flag, npar, np, nmaps
    real(dp)     :: a, b, a_tot, b_tot, s, t1, t2, x_min, x_max, delta_lnL_threshold
    real(dp)     :: mu, sigma, par, w, mu_p, sigma_p, a_old, chisq, chisq_old, chisq_tot, unitconv
    real(dp)     :: buff1_r(1), buff2_r(1), theta_min, theta_max
    logical(lgt) :: ok
    logical(lgt), save :: first_call = .true.
    class(comm_comp), pointer :: c => null()
    !Following is for the local sampler
    real(dp)     :: mixing_old, mixing_new, lnL_new, lnL_old, res_lnL, delta_lnL, lnL_prior, lnL_init
    real(dp)     :: accept_rate, accept_scale, lnL_sum, proplen, chisq_jeffreys
    integer(i4b) :: i_s, p_min, p_max, pixreg_nprop, band_count, pix_count, buff1_i(1), buff2_i(1), burn_in
    integer(i4b) :: n_spec_prop, n_accept, n_corr_prop, n_prop_limit, n_corr_limit, corr_len, out_every
    integer(i4b) :: npixreg, smooth_scale
    logical(lgt) :: first_sample, loop_exit, use_det, burned_in, sampled_nprop, sampled_proplen
    real(dp),      allocatable, dimension(:) :: all_thetas, data_arr, invN_arr, mixing_old_arr, mixing_new_arr
    real(dp),      allocatable, dimension(:) :: theta_corr_arr, old_thetas, new_thetas, init_thetas
    real(dp),      allocatable, dimension(:) :: old_theta_smooth, new_theta_smooth
    integer(i4b),  allocatable, dimension(:) :: band_i, pol_j
    class(comm_mapinfo),             pointer :: info => null()
    class(comm_map),                 pointer :: theta_single_pol => null() ! Spectral parameter of polt_id index
    class(comm_map),                 pointer :: theta_single_pol_smooth => null() ! smoothed spec. ind. of polt_id index
    type(map_ptr), allocatable, dimension(:) :: df


    info  => comm_mapinfo(c_lnL%x%info%comm, c_lnL%x%info%nside, &
         & c_lnL%x%info%lmax, c_lnL%x%info%nmaps, c_lnL%x%info%pol)
                

    delta_lnL_threshold = 25.d0
    n                   = 101
    n_ok                = 50
    burn_in             = 20
    burned_in           = .false.
    first_call          = .false.
    
    n_prop_limit = 20
    n_corr_limit = 20
    out_every    = 10
    !p_lnL = p

    ! we have allready got globally from sampleDiffuseSpecInd: c_lnL, p_lnL (also in input), id_lnL (also in input),
    !     buffer_lnL (theta, limited by min/max)

    if (c_lnL%poltype(id) == 1) then
       p_min = 1; p_max = c_lnL%nmaps
       if (only_pol) p_min = 2
    else if (c_lnL%poltype(id) == 2) then
       if (p == 1) then
          p_min = 1; p_max = 1
       else
          p_min = 2; p_max = c_lnL%nmaps
       end if
    else if (c_lnL%poltype(id) == 3) then
       p_min = p
       p_max = p
    else
       write(*,*) 'Unsupported polarization type'
       stop
    end if

    npar      = c_lnL%npar
    np        = self%x_smooth%info%np
    nmaps     = self%x_smooth%info%nmaps
    npixreg   = self%npixreg(p,id)

    theta_min = c_lnL%p_uni(1,id)
    theta_max = c_lnL%p_uni(2,id)

    !Pixregs are going to have a slightly different structure than fullsky (fullsky will be changed later too)

    !1) we need active bands for evaluation.

    !2) We need to set up sampling of theta (nprop, proplen)

    !2.1) draw new theta for each pixel reg, bcast and smooth theta map

    !2.2) evaluate loglike. If chisq use normal chisq code. If lowres_chisq use trygves new code. 
    !if marg/ridge use underlying code (for now)

    !step 1) 
    !set up which bands and polarizations to include, make sure to set masks of unincluded bands/polarizations to zero

    allocate(band_i(3*numband),pol_j(3*numband))

    band_count=0
    do k = 1,numband !run over all active bands
       !check if the band is associated with the smoothed component, and if band frequencies are within 
       !freq. limits for the component
       if (.not. associated(rms_smooth(k)%p)) cycle
       if (data(k)%bp(0)%p%nu_c < self%nu_min_ind(id) .or. data(k)%bp(0)%p%nu_c > self%nu_max_ind(id)) cycle
                   
       do l = p_min,p_max !loop all maps of band k associated with p (given poltype)
          if (l > data(k)%info%nmaps) cycle !no map in band k for polarization l
          band_count = band_count + 1
          band_i(band_count)=k
          pol_j(band_count)=l
       end do
    end do

    if (band_count==0) then
       buffer_lnL(:,p_min:p_max)=c_lnL%p_gauss(1,id) !set theta to prior, as no bands are valid, no data
       deallocate(band_i,pol_j)
       if (info%myid == 0)  write(*,*) 'no data bands available for sampling of spec ind'
       return
    else
       if (info%myid == 0 .and. .true.) then
             write(*,fmt='(a)') '  Active bands'
          do k = 1,band_count
             write(*,fmt='(a,i1)') '   band: '//trim(data(band_i(k))%label)//', -- polarization: ',pol_j(k)
          end do
       end if
    end if

    ! This is used for marginal/ridge sampling
    allocate(all_thetas(npar))
    allocate(mixing_old_arr(band_count),mixing_new_arr(band_count),invN_arr(band_count),data_arr(band_count))

    !This is used for all
    allocate(old_thetas(0:npixreg),new_thetas(0:npixreg),init_thetas(0:npixreg))
    allocate(old_theta_smooth(0:np-1), new_theta_smooth(0:np-1))

    !This is used for chisq
    if (self%apply_jeffreys) then
       allocate(df(numband))
       do k = 1, numband
          df(k)%p => comm_map(data(k)%info)
       end do
    end if

    n_spec_prop = 0
    n_accept = 0
    n_corr_prop = 0 
    sampled_nprop = .true.
    sampled_proplen = .true.

    if (self%pol_sample_nprop(p,id) .or. self%pol_sample_proplen(p,id)) then
       pixreg_nprop = 10000 !should be enough to find correlation length (and proposal length if prompted) 
       allocate(theta_corr_arr(pixreg_nprop))
       !self%pol_sample_nprop(j,p,id) = boolean array of size (n_pixreg,poltype)
       !self%pol_sample_proplen(j,p,id) = boolean array of size (n_pixreg,poltype)
       if (self%pol_sample_nprop(p,id)) sampled_nprop = .true.
       if (self%pol_sample_proplen(p,id)) sampled_proplen = .true.
    else
       if (np > 0) then
          pixreg_nprop = IDINT(sum(self%pol_nprop(id)%p%map(:,p))/np) ! comm_map is real(dp), use IDINT() to convert to integer
       else
          pixreg_nprop = 0
       end if
    end if

    first_sample  = .true.
    loop_exit     = .false.

    !Init a smoothing map with nmaps = 1
    
    call wall_time(t1)
    lnl_init=0.d0

    !Step 2)

    do j = 1,pixreg_nprop !propose new sample n times (for faster mixing in total)
       if (info%myid == 0) then

          if (first_sample) then
             if (np > 0) then
                proplen=sum(self%pol_proplen(id)%p%map(:,p))/np
             else
                proplen=1.d0
             end if

             old_thetas = self%theta_pixreg(:,p,id)
             old_thetas = min(max(old_thetas,theta_min),theta_max)
             ! that the root processor operates on
             init_thetas = old_thetas
             write(*,fmt='(a, f10.5)') "  initial (avg) sepc. ind. value: ", sum(init_thetas)/npixreg
          end if

          !draw new spec ind, using same proposal length per pixreg
          !may add code later to support sampling of one pixelregion at the time later

          new_thetas(0)=old_thetas(0) !prior value, not to be sampled
          do i = 1,npixreg
             new_thetas(i) = old_thetas(i) + proplen*rand_gauss(handle)
          end do

       end if


       !broadcast new_theta, and new old_theta, and all_thetas to the other processors
       if (first_sample) call mpi_bcast(old_thetas(0:npixreg), npixreg+1, MPI_DOUBLE_PRECISION, &
            & 0, c_lnL%comm, ierr)
       call mpi_bcast(new_thetas(0:npixreg), npixreg+1, MPI_DOUBLE_PRECISION, 0, c_lnL%comm, ierr)

       if (info%myid==0 .and. mod(j,out_every)==0) then
          write(*,fmt='(a, i6, a, f10.5, a, f10.5)') "  proposal: ", j," -- Current (avg) ind (per pixreg): ", &
               & sum(old_thetas(1:npixreg))/npixreg, " -- proposed (avg) ind (per pixreg): ",  &
               & sum(new_thetas(1:npixreg))/npixreg
       end if


       !init log-like for new sample
       if (first_sample) lnL_old = 0.d0
       lnL_new = 0.d0       

       if (any(new_thetas > theta_max) .or. any(new_thetas < theta_min)) then !new proposal is outside limits
          lnL_new = -1.d30 
          ! skip the true caclulation of lnL, we reject the sample ~100%
          if (info%myid==0) write(*,fmt='(a, f10.5, a, f10.5)') "    Proposed ind (partially) outside limits.  min: ", &
               & theta_min," -- max: ", theta_max
       else
          !Step 2.1) set up theta maps and smooth, using beam equal to post-processing beam of the given smoothing scale

          if (first_sample) then
             !set up the old theta map
             do pix=0,np-1
                old_theta_smooth(pix)=old_thetas(self%ind_pixreg_arr(pix,p,id))
             end do
          end if
          !set up the new theta map
          do pix=0,np-1
             new_theta_smooth(pix)=new_thetas(self%ind_pixreg_arr(pix,p,id))
          end do

          smooth_scale = c_lnL%smooth_scale(id)
          if (cpar%num_smooth_scales > 0) then
             if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then
                if (p_min<=p_max) then !just a guaranty that we dont smooth for nothing
                   ! Smooth index map with a postprocessing beam
                   !deallocate(c%theta_smooth)             
                   
                   info  => comm_mapinfo(c_lnL%theta(id)%p%info%comm, cpar%nside_smooth(smooth_scale), &
                        & cpar%lmax_smooth(smooth_scale), 1, c_lnL%theta(id)%p%info%pol) !only want 1 map

                   !spec. ind. map with 1 map (will be smoothed like zero spin map using the existing code)
                   theta_single_pol => comm_map(info)
                   if (first_sample) then
                      theta_single_pol%map(:,1) = old_theta_smooth
                      !smooth single map as intensity (i.e. zero spin)
                      call smooth_map(info, .false., &
                           & data(1)%B_postproc(smooth_scale)%p%b_l*0.d0+1.d0, theta_single_pol, &  
                           & data(1)%B_postproc(smooth_scale)%p%b_l, theta_single_pol_smooth)

                      ! assign smoothed theta map to theta buffers
                      old_theta_smooth = theta_single_pol_smooth%map(:,1)

                      call theta_single_pol_smooth%dealloc()
                      theta_single_pol_smooth => null()
                   end if
                   theta_single_pol%map(:,1) = new_theta_smooth
                   !smooth single map as intensity (i.e. zero spin)
                   call smooth_map(info, .false., &
                        & data(1)%B_postproc(smooth_scale)%p%b_l*0.d0+1.d0, theta_single_pol, &  
                        & data(1)%B_postproc(smooth_scale)%p%b_l, theta_single_pol_smooth)

                   ! assign smoothed theta map to theta buffers
                   new_theta_smooth = theta_single_pol_smooth%map(:,1)

                   call theta_single_pol_smooth%dealloc()
                   theta_single_pol_smooth => null()


                   call theta_single_pol%dealloc()
                   theta_single_pol => null()

                end if
             end if
          end if

          pix_count = 0
          !lnL type should split here
          if (trim(self%pol_lnLtype(p,id))=='chisq') then

             !!due to circular dependencies, we have to comment out this for now
             !!use normal, fullsky chisq
             !if (first_sample) then
             !   do i = p_min,p_max
             !      self%theta(id)%p%map(:,i)=old_theta_smooth
             !   end do
             !   
             !   ! Update mixing matrix with new alms
             !   if (c_lnL%apply_jeffreys) then
             !      call self%updateMixmat(df=df, par=id)
             !      call compute_jeffreys_prior(c_lnL, df, p, id, chisq_jeffreys)
             !   else
             !      call self%updateMixmat
             !   end if
             !
             !   ! Calculate proposed chisq
             !   if (allocated(c_lnL%indmask)) then
             !      call compute_chisq(c_lnL%comm, chisq_fullsky=lnL_old, mask=c_lnL%indmask)
             !   else
             !      call compute_chisq(c_lnL%comm, chisq_fullsky=lnl_old)
             !   end if
             !   if (c_lnL%apply_jeffreys) lnL_old = lnL_old + chisq_jeffreys
             !   lnL_old = -0.5d0*lnL_old
             !end if
             !
             !do i = p_min,p_max
             !   self%theta(id)%p%map(:,i)=new_theta_smooth
             !end do
             !
             !! Update mixing matrix with new alms
             !if (c_lnL%apply_jeffreys) then
             !   call self%updateMixmat(df=df, par=id)
             !   call compute_jeffreys_prior(c_lnL, df, p, id, chisq_jeffreys)
             !else
             !   call self%updateMixmat
             !end if
             !
             !! Calculate proposed chisq
             !if (allocated(c_lnL%indmask)) then
             !   call compute_chisq(c_lnL%comm, chisq_fullsky=lnL_new, mask=c_lnL%indmask)
             !else
             !   call compute_chisq(c_lnL%comm, chisq_fullsky=lnl_new)
             !end if
             !if (c_lnL%apply_jeffreys) lnL_new = lnL_new + chisq_jeffreys
             !lnL_new = -0.5d0*lnL_new
             !
             lnL_old=0.d0
             lnL_new=0.d0
             
          else if (trim(self%pol_lnLtype(p,id))=='chisq_lowres') then
             !use low resolution (nside 16?) chisq evaluation
             !need to add code for support later
             lnL_old=0.d0
             lnL_new=0.d0

          else if (trim(self%pol_lnLtype(p,id))=='chisq_smooth') then

             !The old, evaluate all chisq at smoothing scale, chisq evaluation
             do pix = 0,np-1 !loop over pixels covered by the processor
                if (self%pol_ind_mask(id)%p%map(pix,p) < 0.5d0) cycle     ! if pixel is masked out

                all_thetas(id)=new_theta_smooth(pix)
                !get the values of the remaining spec inds of the component for the given pixel
                do i = 1, npar
                   if (i == id) cycle
                   all_thetas(i) = self%theta_smooth(i)%p%map(pix,p) 
                end do

                do k = 1,band_count !run over all active bands
                   pix_count = pix_count + 1

                   ! get conversion factor from amplitude to data (i.e. mixing matrix element)           
                   ! both for old and new spec. ind. and calc. log likelihood (chisq)
                   if (first_sample) then
                      all_thetas(id)=old_theta_smooth(pix)
                      mixing_old = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                           & data(band_i(k))%gain * c_lnL%cg_scale
                      all_thetas(id)=new_theta_smooth(pix)

                      res_lnL = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_old* &
                           & self%x_smooth%map(pix,pol_j(k))
                      lnL_old = lnL_old -0.5d0*(res_lnL*rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k)))**2

                   end if
                   mixing_new = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                        & data(band_i(k))%gain * c_lnL%cg_scale

                   res_lnL = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_new* &
                        & self%x_smooth%map(pix,pol_j(k))
                   lnL_new = lnL_new -0.5d0*(res_lnL*rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k)))**2
                end do
             end do

          else if ((trim(self%pol_lnLtype(p,id))=='ridge') .or. &
               & (trim(self%pol_lnLtype(p,id))=='marginal')) then
             if (trim(self%pol_lnLtype(p,id))=='ridge') then
                use_det = .false.
             else
                use_det = .true.
             end if

             do pix = 0,np-1 !More precise, we need to loop over pixels covered by the processor
                
                if (self%pol_ind_mask(id)%p%map(pix,p) < 0.5d0) cycle     ! if pixel is masked out

                all_thetas(id)=new_theta_smooth(pix)
                !get the values of the remaining spec inds of the component for the given pixel
                do i = 1, npar
                   if (i == id) cycle
                   all_thetas(i) = self%theta_smooth(i)%p%map(pix,p) 
                end do

                !build mixing matrix
                do k = 1,band_count !run over all active bands
                   pix_count = pix_count + 1

                   if (first_sample) then
                      all_thetas(id)=old_theta_smooth(pix)
                      mixing_old_arr(k) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                           & data(band_i(k))%gain * c_lnL%cg_scale
                      all_thetas(id)=new_theta_smooth(pix)
                   end if
                   mixing_new_arr(k) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                        & data(band_i(k))%gain * c_lnL%cg_scale
                   pix_count = pix_count + 1
                   data_arr(k)=res_smooth(band_i(k))%p%map(pix,pol_j(k))
                   invN_arr(k)=rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k))**2 !assumed diagonal and uncorrelated 
                end do

                !compute the marginal/ridge log-like for the pixel
                if (first_sample) then
                   lnL_old = lnL_old + comp_lnL_marginal_diagonal(mixing_old_arr, invN_arr, data_arr, &
                        & use_det, band_count)
                end if
                lnL_new = lnL_new + comp_lnL_marginal_diagonal(mixing_new_arr, invN_arr, data_arr, &
                     & use_det, band_count)
             end do
          else 
             if (info%myid==0) write(*,*) 'invalid polarized lnL sampler type'
             stop

          end if !chisq/marginal/ridge


          !Sum lnL_new (and lnL_old) among all processors
          if (first_sample) then
             !there is a problem with the mpi_allreduce call. Commander3 crashes
             call mpi_allreduce(MPI_IN_PLACE, lnL_old, 1, MPI_DOUBLE_PRECISION, & 
                  & MPI_SUM, info%comm, ierr)

             if (.not. ( (trim(self%pol_lnLtype(p,id))=='chisq') .or. (trim(self%pol_lnLtype(p,id))=='chisq_lowres'))) then
                call mpi_allreduce(MPI_IN_PLACE, pix_count, 1, MPI_INTEGER, & 
                     & MPI_SUM, info%comm, ierr)
                
                if (pix_count == 0) then! there are no data that is valid
                   old_thetas = c_lnL%p_gauss(1,id) !set theta to prior
                   !set nprop to 0, write error message, then exit sampling
                   self%pol_nprop(id)%p%map(:,p) = 0.d0
                   if (info%myid==0) write(*,*) "    No valid pixels to sample spectral index from"
                   exit
                end if
             end if
          end if

          call mpi_allreduce(MPI_IN_PLACE, lnL_new, 1, MPI_DOUBLE_PRECISION, & 
               & MPI_SUM, info%comm, ierr)

          !add spec. ind. prior
          if (c_lnL%p_gauss(2,id) > 0.d0) then
             !Find average prior "chisq" and add it to lnL
             if (first_sample) then
                lnL_prior=0.d0
                do pix = 0,np-1
                   lnL_prior = (old_theta_smooth(pix)-c_lnL%p_gauss(1,id))**2
                end do
                if (np > 0) lnL_prior = -0.5d0 * lnl_prior/c_lnL%p_gauss(2,id)**2
                pix_count=np
                call mpi_allreduce(MPI_IN_PLACE, lnL_prior, 1, MPI_DOUBLE_PRECISION, & 
                     & MPI_SUM, info%comm, ierr)
                call mpi_allreduce(MPI_IN_PLACE, pix_count, 1, MPI_INTEGER, & 
                     & MPI_SUM, info%comm, ierr)
                lnL_old = lnL_old + lnL_prior/pix_count
                lnl_init = lnL_old
             end if
             lnL_prior=0.d0
             do pix = 0,np-1
                lnL_prior = (new_theta_smooth(pix)-c_lnL%p_gauss(1,id))**2
             end do
             if (np > 0) lnL_prior = -0.5d0 * lnl_prior/c_lnL%p_gauss(2,id)**2
             pix_count=np
             call mpi_allreduce(MPI_IN_PLACE, lnL_prior, 1, MPI_DOUBLE_PRECISION, & 
                  & MPI_SUM, info%comm, ierr)
             call mpi_allreduce(MPI_IN_PLACE, pix_count, 1, MPI_INTEGER, & 
                  & MPI_SUM, info%comm, ierr)
             lnL_new = lnL_new + lnL_prior/pix_count
          end if

          !first sample done
          first_sample = .false.

       end if !new_theta outside spec limits

       if (info%myid == 0) then
          !accept/reject new spec ind
          delta_lnL = lnL_new-lnL_old

          if (mod(j,out_every)==0) write(*,fmt='(a, e14.5)') "    lnL_new - lnL_old = ", delta_lnL


          if (delta_lnL < 0.d0) then
             !if delta_lnL is more negative than -25.d0, limit to -25.d0
             if (abs(delta_lnL) > delta_lnL_threshold) delta_lnL = -delta_lnL_threshold 

             !draw random uniform number
             a = rand_uni(handle) !draw uniform number from 0 to 1
             if (exp(delta_lnL) > a) then
                !accept
                old_thetas = new_thetas
                lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
                n_accept = n_accept + 1
             end if
          else
             !accept new sample, higher likelihood
             old_thetas = new_thetas
             lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
             n_accept = n_accept + 1
          end if

          if (j > burn_in) burned_in = .true.
          ! evaluate proposal_length, then correlation length. 
          !This should only be done the first gibbs iteration, if prompted to do so from parameter file.
          if (self%pol_sample_proplen(p,id)) then
             n_spec_prop = n_spec_prop + 1
             accept_rate = n_accept*1.d0/n_spec_prop
             if (mod(n_spec_prop,out_every)==0) write(*,fmt='(a, f6.4)') "   accept rate = ", accept_rate
             if (.not. burned_in) then 
                !giving som time to let the first steps "burn in" if one starts far away from max lnL
                n_spec_prop = 0
                n_accept = 0
             else if (n_spec_prop > n_prop_limit) then
                if (accept_rate > 0.7d0) then
                   accept_scale = 1.5d0
                else if (accept_rate < 0.3d0) then
                   accept_scale = 0.5d0
                else 
                   accept_scale = -1.d0
                end if

                if (accept_scale > 0.d0) then
                   !prop_len only ever used by root processor
                   proplen = proplen*accept_scale !set new proposal length
                   n_accept = 0 !reset with new prop length
                   n_spec_prop = 0 !reset with new prop length
                   write(*,fmt='(a, e14.5)') "     New prop. len. = ", proplen
                else
                   self%pol_sample_proplen(p,id)=.false. !making sure we dont go back into prop. len. sampling
                   !only necessary to update for myid==0
                   if (.not. self%pol_sample_nprop(p,id)) then
                      loop_exit=.true. !we have found the ideal proposal length, not looking for correlation length
                   end if
                end if
             end if

          else if (self%pol_sample_nprop(p,id)) then
             ! if we are keeping track of spec inds for corr. length, then save current spec ind
             n_corr_prop = n_corr_prop + 1
             theta_corr_arr(n_corr_prop) = sum(old_thetas(1:npixreg))/npixreg
             if (.not. burned_in) then
                n_corr_prop = 0
             else if (n_corr_prop > n_corr_limit) then
                corr_len = calc_corr_len(theta_corr_arr(1:n_corr_prop),n_corr_prop) !returns -1 if no corr.len. is found
                if (corr_len > 0) then ! .or. n_corr_prop > 4*self%nprop_uni(2,id)) then
                   !set pixreg_nprop (for pix region) to x*corr_len, x=2, but inside limits from parameter file
                   self%pol_nprop(id)%p%map(:,p) = 1.d0*min(max(2*corr_len,self%nprop_uni(1,id)),self%nprop_uni(2,id))
                   loop_exit=.true. !we have found the correlation length
                   self%pol_sample_nprop(p,id)=.false. !making sure we dont go back into corr. len. sampling
                end if
             end if
          else
             accept_rate = n_accept*1.d0/j
             if (mod(j,out_every)==0) write(*,fmt='(a, f6.4)') "   accept rate = ", accept_rate
          end if
          
       end if

       call wall_time(t2)

       if (info%myid == 0 .and. mod(j,out_every)==0) write(*,*) '   Sample:',j,'  Wall time per sample:',real((t2-t1)/j,sp)

       call mpi_bcast(loop_exit, 1, MPI_LOGICAL, 0, c_lnL%comm, ierr)

       if (loop_exit) exit
       
    end do !nprop

    !bcast proposal length
    call mpi_bcast(proplen, 1, MPI_DOUBLE_PRECISION, 0, &
         & c_lnL%comm, ierr)
    self%pol_proplen(id)%p%map(:,p) = proplen

    !bcast number of proposals
    call mpi_bcast(self%pol_nprop(id)%p%map(0,p), 1, MPI_DOUBLE_PRECISION, 0, &
         & c_lnL%comm, ierr)
    self%pol_nprop(id)%p%map(:,p) = self%pol_nprop(id)%p%map(0,p)

    !bcast last valid theta to all procs to update theta map
    call mpi_bcast(old_thetas(0:npixreg), npixreg+1, MPI_DOUBLE_PRECISION, 0, c_lnL%comm, ierr)
    !assign thatas to pixel regions
    do pix=0,np-1
       old_theta_smooth=old_thetas(self%ind_pixreg_arr(pix,p,id))
    end do

    do i = p_min,p_max
       buffer_lnL(:,i) = old_theta_smooth
    end do
    self%theta_pixreg(0:npixreg,p,id)=old_thetas

    if (info%myid==0) then
       write(*,fmt='(a, f10.5)') "    Final (avg) spec. ind. value: ", sum(old_thetas(1:npixreg)/npixreg)
       write(*,fmt='(a, i5, a, e14.5)') "    Number of proposals: ",IDINT(self%pol_nprop(id)%p%map(0,p)), &
            & "  --  Proposal length: ", proplen
       write(*,fmt='(a, e14.5)') "    Difference in (avg) spec. ind., new - old: ", &
            & sum((old_thetas(1:npixreg)-init_thetas(1:npixreg))/npixreg)
       write(*,fmt='(a, e14.5)') "    New Log-Likelihood                      : ", &
            & lnl_old
       write(*,fmt='(a, e14.5)') "    Difference in Log-Likelihood (new - old): ", &
            & lnl_old-lnl_init
       write(*,*) ''
    end if

    if (iter >= self%sample_first_niter) then !only stop sampling after n'th iteration, in case first iter is far off.
       self%pol_sample_nprop(p,id) = .false. !set sampling of correlation length (number of proposals) and
       self%pol_sample_proplen(p,id) = .false. !proposal length to false (this is only to be done in the first 
                                               !gibbs iteration anyways)
    else !reset sampling to true if we have sampled
       if (sampled_nprop) self%pol_sample_nprop(p,id) = .true. 
       if (sampled_proplen) self%pol_sample_proplen(p,id) = .true. 
    end if

    !deallocate
    deallocate(mixing_old_arr,mixing_new_arr,data_arr,invN_arr,all_thetas)
    deallocate(band_i,pol_j)
    if (allocated(theta_corr_arr)) deallocate(theta_corr_arr)
    deallocate(old_thetas,new_thetas,init_thetas,old_theta_smooth, new_theta_smooth)
    if (allocated(df)) deallocate(df)


  end subroutine sampleDiffuseSpecIndPixReg

  function calc_corr_len(spec_arr,n_spec) 
    implicit none
    integer(i4b),               intent(in) :: n_spec
    real(dp),     dimension(:), intent(in) :: spec_arr
    integer(i4b)                           :: calc_corr_len

    integer(i4b) :: i, j, maxcorr, ns
    real(dp)     :: x_mean, y_mean, sig_x, sig_y
    real(dp), dimension(:), allocatable :: xarr, yarr,covarr

    calc_corr_len = -1 !default if no corr length is found

    maxcorr = n_spec/2
    allocate(xarr(maxcorr),yarr(maxcorr),covarr(maxcorr))

    do i = 1,maxcorr
       ns=maxcorr-i
       xarr=spec_arr(1:ns)
       yarr=spec_arr(1+i:ns+i)
       x_mean=sum(xarr)/ns
       y_mean=sum(yarr)/ns
       sig_x = sqrt(sum((xarr-x_mean)**2)/ns)
       sig_y = sqrt(sum((yarr-y_mean)**2)/ns)
       covarr(i) = sum((xarr-x_mean)*(yarr-y_mean))/(ns*sig_x*sig_y)
       if (covarr(i) < 0.d0) then
          calc_corr_len = i
          exit
       end if
    end do

    deallocate(xarr,yarr,covarr)

  end function calc_corr_len

  function comp_lnL_marginal_diagonal(mixing,invN_arr,data,use_det,arr_len)
    implicit none
    logical(lgt),               intent(in)           :: use_det
    real(dp),     dimension(:), intent(in)           :: mixing
    real(dp),     dimension(:), intent(in)           :: invN_arr
    real(dp),     dimension(:), intent(in)           :: data
    integer(i4b),               intent(in), optional :: arr_len
    real(dp)                                         :: comp_lnL_marginal_diagonal

    integer(i4b) :: i, j, mat_len
    real(dp)     :: MNd,MNM
    real(dp), dimension(:), allocatable :: MN

    if (present(arr_len)) then
       allocate(MN(arr_len))
       MN=mixing(1:arr_len)*invN_arr(1:arr_len)
       MNd=sum(MN*data(1:arr_len))
       MNM=sum(MN*mixing(1:arr_len))
    else
       allocate(MN(size(mixing)))
       MN=mixing*invN_arr
       MNd=sum(MN*data)
       MNM=sum(MN*mixing)
    end if

    comp_lnL_marginal_diagonal = 0.d0

    if (MNM /= 0.d0) then 
       MNM=1.d0/MNM !invert 1x1 matrix
    else
       comp_lnL_marginal_diagonal=1.d30 !MNM = 0.d0, i.e. no inversion possible 
       deallocate(MN)
       return
    end if

    comp_lnL_marginal_diagonal = 0.5d0*MNd*MNM*MNd

    !determinant of 1x1 matrix is the value of the matrix itself
    if (use_det) comp_lnL_marginal_diagonal = comp_lnL_marginal_diagonal - 0.5d0*log(MNM) 

    deallocate(MN)


  end function comp_lnL_marginal_diagonal

  ! Sample spectral parameters
  subroutine computeMarginalLogLikelihood(self, handle, id, use_det, marg_fullsky, marg_map, theta_map, theta_mask, p_ind)
    implicit none
    class(comm_diffuse_comp),                intent(inout)           :: self
    type(planck_rng),                        intent(inout)           :: handle
    integer(i4b),                            intent(in)              :: id
    logical(lgt),                            intent(in)              :: use_det
    real(dp),                                intent(out),   optional :: marg_fullsky
    class(comm_map),                         intent(inout), optional :: marg_map
    class(comm_map),                         intent(in),    optional :: theta_map
    class(comm_map),                         intent(in),    optional :: theta_mask
    !type(map_ptr),            dimension(1:), intent(in),    optional :: data_mask !data masks, if they are to be used they must have been ud_graded to correct nside
    integer(i4b),                            intent(in),    optional :: p_ind !poltype index, only compute for the given poltype ind

    integer(i4b) :: i, j, k, l, n, p, q, pix, ierr, ind(1), counter, n_ok
    integer(i4b) :: i_min, i_max, status, n_gibbs, iter, n_pix, n_pix_tot, flag, npar, np, nmaps
    real(dp)     :: a, b, a_tot, b_tot, s, t1, t2, x_min, x_max, delta_lnL_threshold
    real(dp)     :: mu, sigma, par, w, mu_p, sigma_p, a_old, chisq, chisq_old, chisq_tot, unitconv
    real(dp)     :: x(1), theta_min, theta_max
    logical(lgt) :: ok
    logical(lgt), save :: first_call = .true.
    class(comm_comp), pointer :: c => null()
    !Following is for the local sampler
    real(dp)     :: lnL_new
    integer(i4b) :: i_s, p_min, p_max, band_count
    real(dp),     allocatable, dimension(:) :: all_thetas, data_arr, invN_arr, mixing_arr
    real(dp),     allocatable, dimension(:) :: theta_corr_arr
    integer(i4b), allocatable, dimension(:) :: band_i, pol_j
    class(comm_mapinfo),            pointer :: info => null()


    if (.not. (present(marg_fullsky) .or. present(marg_map))) return !no log-like output specified

    delta_lnL_threshold = 25.d0
    n                   = 101
    n_ok                = 50

    !c_lnL       => self
    id_lnL      = id
    c           => compList     ! Extremely ugly hack...
    do while (self%id /= c%id)
       c => c%next()
    end do
    select type (c)
    class is (comm_diffuse_comp)
       c_lnL => c
    end select

    info  => comm_mapinfo(c_lnL%x%info%comm, c_lnL%x%info%nside, &
         & c_lnL%x%info%lmax, c_lnL%x%info%nmaps, c_lnL%x%info%pol)


    npar      = c%npar
    np        = self%x_smooth%info%np
    nmaps     = self%x_smooth%info%nmaps

    theta_min = c_lnL%p_uni(1,id_lnL)
    theta_max = c_lnL%p_uni(2,id_lnL)
       
    allocate(buffer_lnL(0:self%x_smooth%info%np-1,self%x_smooth%info%nmaps))
    if (present(theta_map)) then
       buffer_lnL=max(min(theta_map%map,theta_max),theta_min) 
    else
       buffer_lnL=max(min(self%theta(id)%p%map,theta_max),theta_min) 
    end if

    ! Choose which polarization fields to include in evaluation
    if (present(p_ind)) then
       if (p_ind == 1) then
          if (c_lnL%poltype(id_lnL) == 1) then
             p_min = 1; p_max = c_lnL%nmaps
             if (only_pol) p_min = 2
          else 
             p_min = 1; p_max = 1
          end if
       else if (p_ind == 2) then
          if (c_lnL%poltype(id_lnL) == 1) then
             write(*,*) 'Trying to evaluate poltype index (', p_ind, ') higher than defined poltype'
             write(*,*) 'Component '//trim(self%label)//', spec. ind. '//trim(self%indlabel(id))
             stop
          else if (c_lnL%poltype(id_lnL) == 2) then
             p_min = 2; p_max = c_lnL%nmaps
          else if (c_lnL%poltype(id_lnL) == 3) then
             p_min = 2; p_max = 2
          else
             write(*,*) 'Unsupported polarization type'
             write(*,*) 'Component '//trim(self%label)//', spec. ind. '//trim(self%indlabel(id))
             stop
          end if
       else if  (p_ind == 3) then
          if (c_lnL%poltype(id_lnL) < 3) then
             write(*,*) 'Trying to evaluate poltype index (', p_ind, ') higher than defined poltype in'
             write(*,*) 'Component '//trim(self%label)//', spec. ind. '//trim(self%indlabel(id))
             stop
          else if (c_lnL%poltype(id_lnL) == 3) then
             p_min = 3; p_max = 3
          else
             write(*,*) 'Unsupported polarization type'
             write(*,*) 'Component '//trim(self%label)//', spec. ind. '//trim(self%indlabel(id))
             stop
          end if
       else
          write(*,*) 'Unsupported polarization index', p_ind, 'in marginal log-likelihood evaluation' 
          stop
       end if
    else
       p_min = 1; p_max = c_lnL%nmaps
       if (only_pol) p_min = 2
    end if

    allocate(band_i(3*numband),pol_j(3*numband))

    band_count=0
    do k = 1,numband !run over all active bands
       !check if the band is associated with the smoothed component, and if band frequencies are within 
       !freq. limits for the component
       if (.not. associated(rms_smooth(k)%p)) cycle
       if (data(k)%bp(0)%p%nu_c < self%nu_min_ind(id) .or. data(k)%bp(0)%p%nu_c > self%nu_max_ind(id)) cycle
                   
       do l = p_min,p_max !loop all maps of band k associated with p (given poltype)
          if (l > data(k)%info%nmaps) cycle !no map in band k for polarization l
          band_count = band_count + 1
          band_i(band_count)=k
          pol_j(band_count)=l
       end do
    end do

    if (band_count==0) then
       if (present(marg_fullsky)) marg_fullsky = 0.d0
       if (present(marg_map))     marg_map%map = 0.d0
       deallocate(band_i,pol_j)
       if (info%myid == 0) then
          write(*,*) 'no data bands available for evaluating marginal log-likelihood'
          write(*,*) 'Component '//trim(self%label)//', spec. ind. '//trim(self%indlabel(id))
       end if
       return
    end if

    allocate(all_thetas(npar))
    allocate(mixing_arr(band_count),invN_arr(band_count),data_arr(band_count))

    if (present(marg_fullsky)) marg_fullsky = 0.d0
    if (present(marg_map)) marg_map%map = 0.d0

    do pix = 0,np-1 !loop over the number of pixels the processor is handling
       
       if (present(theta_mask)) then 
          if (theta_mask%map(pix,1) < 0.5d0) cycle
       end if

       
       !! build mixing matrix, invN and data for pixel
       do k = 1,band_count !run over all active bands

          ! get conversion factor from amplitude to data (i.e. mixing matrix element)           
          ! both for old and new spec. ind.
          do i = 1, npar
             if (i == id) then
                all_thetas(i) = buffer_lnL(pix,pol_j(k))
             else
                all_thetas(i) = self%theta_smooth(i)%p%map(pix,pol_j(k))
             end if
          end do

          mixing_arr(k) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
               & data(band_i(k))%gain * c_lnL%cg_scale

          data_arr(k)=res_smooth(band_i(k))%p%map(pix,pol_j(k))
          invN_arr(k)=rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k))**2 !assumed diagonal and uncorrelated 
       end do

       lnL_new = comp_lnL_marginal_diagonal(mixing_arr, invN_arr, data_arr, &
            & use_det, band_count)

       if (present(marg_fullsky)) marg_fullsky = marg_fullsky + lnL_new
       if (present(marg_map)) marg_map%map(pix,1)  = lnL_new
    end do !pix = 1,np

    if (present(marg_fullsky)) then
       call mpi_allreduce(MPI_IN_PLACE, marg_fullsky, 1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
    end if

    deallocate(all_thetas)
    deallocate(invN_arr)
    deallocate(data_arr)
    deallocate(mixing_arr)    
    deallocate(buffer_lnL)

  end subroutine computeMarginalLogLikelihood


  function lnL_diffuse_multi(p)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp)                                     :: lnL_diffuse_multi
    
    integer(i4b) :: i, l, k, q, pix, ierr, flag, p_min, p_max
    real(dp)     :: lnL, amp, s, a
    real(dp), allocatable, dimension(:) :: theta

    allocate(theta(c_lnL%npar))
    theta         = theta_lnL
    theta(id_lnL) = p(1)
    
    ! Check spectral index priors
    do l = 1, c_lnL%npar
       if (theta(l) < c_lnL%p_uni(1,l) .or. theta(l) > c_lnL%p_uni(2,l)) then
          lnL_diffuse_multi = 1.d30
          deallocate(theta)
          return
       end if
    end do

    ! Choose which polarization fields to include
    if (c_lnL%poltype(id_lnL) == 1) then
       p_min = 1; p_max = c_lnL%nmaps
       if (only_pol) p_min = 2
    else if (c_lnL%poltype(id_lnL) == 2) then
       if (p_lnL == 1) then
          p_min = 1; p_max = 1
       else
          p_min = 2; p_max = c_lnL%nmaps
       end if
    else if (c_lnL%poltype(id_lnL) == 3) then
       p_min = p_lnL
       p_max = p_lnL
    else
       write(*,*) 'Unsupported polarization type'
       stop
    end if

!!$    if (c_lnL%x%info%pix(k_lnL) == 10000) then
!!$       write(*,*)
!!$       write(*,*) 'Theta = ', p(1)
!!$       write(*,*) 'amp   = ', a_lnL(2:3)
!!$       write(*,*) 'p     = ', p_lnL, p_min, p_max
!!$    end if


    lnL = 0.d0
    do k = p_min, p_max

       do l = 1, numband
       !if (c_lnL%x%info%myid == 0) write(*,*) l, numband
          if (.not. associated(rms_smooth(l)%p)) cycle
          if (k > data(l)%info%nmaps) cycle
          
       ! Compute predicted source amplitude for current band
          s = a_lnL(k) * c_lnL%F_int(k,l,0)%p%eval(theta) * data(l)%gain * c_lnL%cg_scale
          
          ! Compute likelihood 
          lnL = lnL - 0.5d0 * (res_smooth(l)%p%map(k_lnL,k)-s)**2 * rms_smooth(l)%p%siN%map(k_lnL,k)**2

!          if (c_lnL%x%info%myid == 0) write(*,fmt='(2i4,3f8.2,f12.2)') k, l, real(res_smooth(l)%p%map(k_lnL,k),sp), real(s,sp), real(1.d0/rms_smooth(l)%p%siN%map(k_lnL,k),sp), real(lnL,sp)

!!$       if (c_lnL%x%info%pix(k_lnL) == 10000) then
!!$          write(*,fmt='(5f10.3,f16.3)') data(l)%bp%nu_c/1.d9, res_smooth(l)%p%map(k_lnL,k), s, res_smooth(l)%p%map(k_lnL,k)-a_lnL(k) * c_lnL%F_int(l)%p%eval(theta) * data(l)%gain * c_lnL%cg_scale, rms_smooth(l)%p%siN%map(k_lnL,k), (res_smooth(l)%p%map(k_lnL,k)-s)**2 * rms_smooth(l)%p%siN%map(k_lnL,k)**2
!!$       end if

       end do
    end do

!!$    call mpi_finalize(i)
!!$    stop

!    if (c_lnL%x%info%myid == 0) write(*,*)
!    if (c_lnL%x%info%myid == 0) write(*,fmt='(i4,f8.4,f12.2)') k_lnL, theta, lnL

    ! Apply index priors
    do l = 1, c_lnL%npar
       if (c_lnL%p_gauss(2,l) > 0.d0) then
          lnL = lnL - 0.5d0 * (theta(l)-c_lnL%p_gauss(1,l))**2 / c_lnL%p_gauss(2,l)**2 
       end if
    end do

!    if (c_lnL%x%info%myid == 0) write(*,fmt='(i4,f8.4,f12.2)') k_lnL, theta, lnL
    
    ! Return chi-square
    lnL_diffuse_multi = -2.d0*lnL

!    if (c_lnL%x%info%myid == 0) write(*,fmt='(i4,f8.4,f12.2)') k_lnL, theta, lnL_diffuse_multi

!!$    if (c_lnL%x%info%pix(k_lnL) == 10000) then
!!$       write(*,*) 'lnL = ', lnL_diffuse_multi
!!$    end if

!!$    if (c_lnL%x%info%pix(k_lnL) == 534044) then
!!$       write(*,*) 'lnL = ', lnL, theta(id_lnL)
!!$       write(*,*)
!!$    end if


    deallocate(theta)

  end function lnL_diffuse_multi




!!$  ! Sample spectral parameters
!!$  subroutine sampleDiffuseSpecInd(self, handle, id)
!!$    implicit none
!!$    class(comm_diffuse_comp),                intent(inout)        :: self
!!$    type(planck_rng),                        intent(inout)        :: handle
!!$    integer(i4b),                            intent(in)           :: id
!!$
!!$    integer(i4b) :: i, j, k, l, n, p, q, pix, ierr, ind(1), counter, n_ok, i_min, i_max, status, n_gibbs, iter, n_pix, n_pix_tot, flag
!!$    real(dp)     :: a, b, a_tot, b_tot, s, t1, t2, x_min, x_max, delta_lnL_threshold, mu, sigma, w, mu_p, sigma_p, a_old, chisq, chisq_tot, unitconv
!!$    logical(lgt) :: ok
!!$    logical(lgt), save :: first_call = .true.
!!$    class(comm_comp), pointer :: c
!!$    real(dp),     allocatable, dimension(:)   :: x, lnL, P_tot, F, theta, a_curr
!!$    real(dp),     allocatable, dimension(:,:) :: amp
!!$
!!$    delta_lnL_threshold = 25.d0
!!$    n                   = 101
!!$    n_ok                = 50
!!$    first_call          = .false.
!!$
!!$    if (trim(operation) == 'optimize') then
!!$
!!$       allocate(amps_lnL(0:data(6)%info%np-1,data(6)%info%nmaps,6:9))
!!$       !allocate(amps_lnL(0:data(1)%info%np-1,data(1)%info%nmaps,numband))
!!$       apply_mixmat = .false.
!!$       do i = 6, 9
!!$          amps_lnL(:,:,i) = self%getBand(i)
!!$       end do
!!$       apply_mixmat = .true.
!!$
!!$
!!$       do p = 1, data(6)%info%nmaps
!!$          do k = 0, data(6)%info%np-1
!!$             p_lnL       = p
!!$             k_lnL       = k
!!$             c           => compList     ! Extremely ugly hack...
!!$             do while (self%id /= c%id)
!!$                c => c%next()
!!$             end do
!!$             select type (c)
!!$             class is (comm_diffuse_comp)
!!$                c_lnL => c
!!$             end select
!!$             
!!$             ! Perform non-linear search
!!$             if (self%p_gauss(2,1) == 0.d0 .and. self%p_gauss(2,2) > 0.d0) then
!!$                ! Dust T
!!$                allocate(x(1))
!!$                x(1)         = self%theta(2)%p%map(k,p)
!!$                a_old_lnL    = self%x%map(k,p)
!!$                beta_old_lnL = self%theta(1)%p%map(k,p)
!!$                call powell(x, lnL_diffuse_multi, ierr)
!!$                self%theta(2)%p%map(k,p) = x(1)
!!$                deallocate(x)
!!$             end if
!!$
!!$             if (self%p_gauss(2,1) > 0.d0 .and. self%p_gauss(2,2) == 0.d0) then
!!$                ! Dust T
!!$                allocate(x(1))
!!$                x(1)         = self%theta(1)%p%map(k,p)
!!$                a_old_lnL    = self%x%map(k,p)
!!$                T_old_lnL    = self%theta(2)%p%map(k,p)
!!$                call powell(x, lnL_diffuse_multi, ierr)
!!$                self%theta(1)%p%map(k,p) = x(1)
!!$                deallocate(x)
!!$             end if
!!$
!!$
!!$             if (self%p_gauss(2,1) > 0.d0 .and. self%p_gauss(2,2) > 0.d0) then
!!$                ! Dust beta and T
!!$                allocate(x(2))
!!$                x(1)         = self%theta(1)%p%map(k,p)
!!$                x(2)         = self%theta(2)%p%map(k,p)
!!$                a_old_lnL    = self%x%map(k,p)
!!$                call powell(x, lnL_diffuse_multi, ierr)
!!$                self%theta(1)%p%map(k,p) = x(1)
!!$                self%theta(2)%p%map(k,p) = x(2)
!!$                deallocate(x)
!!$             end if
!!$
!!$
!!$
!!$
!!$
!!$          end do
!!$       end do
!!$
!!$       ! Update mixing matrix
!!$       call self%updateMixmat
!!$       
!!$       ! Ask for CG preconditioner update
!!$       if (self%cg_unique_sampgroup > 0) recompute_diffuse_precond = .true.
!!$
!!$       deallocate(amps_lnL)
!!$       return
!!$    end if
!!$
!!$  end subroutine sampleDiffuseSpecInd
!!$
!!$
!!$  function lnL_diffuse_multi(p)
!!$    use healpix_types
!!$    implicit none
!!$    real(dp), dimension(:), intent(in), optional :: p
!!$    real(dp)                                     :: lnL_diffuse_multi
!!$    
!!$    integer(i4b) :: i, l, k, q, pix, ierr, flag
!!$    real(dp)     :: lnL, amp, s, a
!!$    real(dp), allocatable, dimension(:) :: theta
!!$    
!!$    allocate(theta(c_lnL%npar))
!!$    if (c_lnL%p_gauss(2,1) > 0.d0 .and. c_lnL%p_gauss(2,2) > 0.d0) then
!!$       theta(1) = p(1)
!!$       theta(2) = p(2)
!!$    else if (c_lnL%p_gauss(2,1) > 0.d0 .and. c_lnL%p_gauss(2,2) == 0.d0) then
!!$       theta(1) = p(1)
!!$       theta(2) = T_old_lnL
!!$    else if (c_lnL%p_gauss(2,1) == 0.d0 .and. c_lnL%p_gauss(2,2) > 0.d0) then
!!$       theta(1) = beta_old_lnL
!!$       theta(2) = p(1)
!!$    end if
!!$    
!!$    ! Check spectral index priors
!!$    do l = 1, c_lnL%npar
!!$       if (theta(l) < c_lnL%p_uni(1,l) .or. theta(l) > c_lnL%p_uni(2,l)) then
!!$          lnL_diffuse_multi = 1.d30
!!$          deallocate(theta)
!!$          return
!!$       end if
!!$    end do
!!$    
!!$    lnL = 0.d0
!!$    do l = 6, 9
!!$    !do l = 1, numband
!!$       
!!$       ! Compute mixing matrix
!!$       s = c_lnL%F_int(l)%p%eval(theta) * data(l)%gain * c_lnL%cg_scale
!!$       
!!$       ! Compute predicted source amplitude for current band
!!$       a = s * amps_lnL(k_lnL,p_lnL,l)
!!$       
!!$       ! Compute likelihood 
!!$       lnL = lnL - 0.5d0 * (data(l)%res%map(k_lnL,p_lnL)-a)**2 / data(l)%N%rms_pix(k_lnL,p_lnL)**2
!!$
!!$    end do
!!$    
!!$    ! Apply index priors
!!$    do l = 1, c_lnL%npar
!!$       if (c_lnL%p_gauss(2,l) > 0.d0) then
!!$          lnL = lnL - 0.5d0 * (theta(l)-c_lnL%p_gauss(1,l))**2 / c_lnL%p_gauss(2,l)**2 
!!$       end if
!!$    end do
!!$    
!!$    ! Return chi-square
!!$    lnL_diffuse_multi = -2.d0*lnL
!!$
!!$    deallocate(theta)
!!$
!!$  end function lnL_diffuse_multi

  
  subroutine print_precond_mat
    implicit none

    integer(i4b) :: l, m, i, j, p
    real(dp), allocatable, dimension(:)   :: W
    real(dp), allocatable, dimension(:,:) :: mat

    if (info_pre%myid /= 0) return
    p = 2

    open(58,file=trim(outdir)//'/precond_W.dat',recl=1024)
    if (trim(precond_type) == 'diagonal') then
       do l = 0, info_pre%lmax
          call info_pre%lm2i(l, 0, i)
          write(*,*) 
          write(*,*) l 
          do j = 1, size(P_cr%invM_diff(i,p)%M(j,:),1)
             write(*,*) real(P_cr%invM_diff(i,p)%M(j,:),sp)
          end do
          allocate(W(P_cr%invM_diff(i,p)%n))
          call get_eigenvalues(P_cr%invM_diff(i,p)%M, W)
          write(58,*) l, real(W,sp)
          deallocate(W)
       end do
    else
       allocate(mat(npre,npre))
       do l = 0, info_pre%lmax
          mat = matmul(P_cr%invM_diff(l,p)%M, transpose(P_cr%invM_diff(l,p)%M))
          write(*,*) 
          write(*,*) l 
          do j = 1, npre
             write(*,*) real(mat(j,:),sp)
          end do
          allocate(W(npre))
          call get_eigenvalues(mat, W)
          write(58,*) l, real(W,sp)
          deallocate(W)
       end do
       deallocate(mat)
    end if
    close(58)


  end subroutine print_precond_mat

  subroutine updateDiffuseFInt(self, band)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self
    integer(i4b),             intent(in),   optional :: band

    integer(i4b) :: i, j, k

    if (present(band)) then
       do i = 1, data(band)%info%nmaps
          do j = 0, data(band)%ndet
             call self%F_int(i,band,j)%p%update(pol=i)
          end do
       end do
    else
       do k = 1, numband
          do i = 1, data(k)%info%nmaps
             do j = 0, data(k)%ndet
                call self%F_int(i,k,j)%p%update(pol=i)
             end do
          end do
       end do
    end if

  end subroutine updateDiffuseFInt

  subroutine updateLowlPrecond(self)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self

    real(dp)                  :: t1, t2
    integer(i4b)              :: i, j, k, l, m, q, lp, mp, myid, nalm, ierr, nmaps
    class(comm_map),     pointer :: map => null(), map2 => null(), tot => null()
    class(comm_mapinfo), pointer :: info => null()
    real(dp),        allocatable, dimension(:,:) :: invM, buffer

    if (self%lmax_pre_lowl < 0) return

    call wall_time(t1)

    ! Initialize arrays
    nalm = (self%lmax_pre_lowl+1)**2
    allocate(invM(0:nalm-1,0:nalm-1), buffer(0:nalm-1,0:nalm-1))
    invM = 0.d0

    info => comm_mapinfo(self%comm, 2, 2*self%lmax_pre_lowl, self%nmaps, self%nmaps==3)
    map  => comm_map(info)
    tot  => comm_map(info)
    do l = 0, self%lmax_pre_lowl
       !if (info%myid == 0) write(*,*) '  Low-ell init l = ', l, self%lmax_pre_lowl, trim(self%label)
       do m = -l, l
          ! Set up unit vector
          call info%lm2i(l,m,j)
          map%alm = 0.d0
          if (j >= 0) map%alm(j,1) = 1.d0

          call self%Cl%sqrtS(alm=map%alm, info=map%info) ! Multiply with sqrt(Cl)

          ! Add frequency dependent terms
          tot%alm = 0.d0
          do i = 1, numband

             ! Compute component-summed map, ie., column-wise matrix elements
             nmaps    = min(self%nmaps, data(i)%info%nmaps)
             info     => comm_mapinfo(self%comm, data(i)%N%nside_chisq_lowres, 2*self%lmax_pre_lowl, &
                  & nmaps, nmaps==3)
             map2     => comm_map(info)
             do k = 1, nmaps
                map2%alm(:,k) = map%alm(:,k) * self%F_mean(i,0,k)
             end do
             if (any(map2%alm /= map2%alm)) then
                write(*,*) 'a', l, m, i
             end if
             call data(i)%B(0)%p%conv(trans=.false., map=map2)
             if (any(map2%alm /= map2%alm)) then
                write(*,*) 'b', l, m, i
             end if
             call map2%Y()
             if (any(map2%map /= map2%map)) then
                write(*,*) 'c', l, m, i
             end if

             ! Multiply with invN
             call data(i)%N%InvN_lowres(map2)
             if (any(map2%map /= map2%map)) then
                write(*,*) 'd', l, m, i
             end if

             ! Project summed map into components, ie., row-wise matrix elements
             call map2%Yt()
             if (any(map2%alm /= map2%alm)) then
                write(*,*) 'e', l, m, i
             end if
             call data(i)%B(0)%p%conv(trans=.true., map=map2)
             if (any(map2%alm /= map2%alm)) then
                write(*,*) 'f', l, m, i
             end if

             map2%alm(:,1) = map2%alm(:,1) * self%F_mean(i,0,1)

             if (any(map2%alm /= map2%alm)) then
                write(*,*) 'g/', l, m, i
             end if


             ! Add up alms
             tot%alm(:,1) = tot%alm(:,1) + map2%alm(:,1)

             call map2%dealloc()
          end do

             if (any(tot%alm /= tot%alm)) then
                write(*,*) 'h', l, m
             end if

          ! Add prior term and multiply with sqrt(S) for relevant components
          call self%Cl%sqrtS(alm=tot%alm, info=tot%info)

             if (any(tot%alm /= tot%alm)) then
                write(*,*) 'i', l, m
             end if


          ! Add (unity) prior term, and store column
          i = l**2 + l + m
          do lp = 0, self%lmax_pre_lowl
             do mp = -lp, lp
                call tot%info%lm2i(lp,mp,k)
                if (k >= 0) then
                   j = lp**2 + lp + mp
                   invM(i,j) = tot%alm(k,1)
!!$                   if (i == 0 .or. j == 0) then
!!$                      write(*,*) i, j, l, m, lp, mp, k, tot%alm(k,1)
!!$                   end if
                   if (i == j) invM(i,j) = invM(i,j) + 1.d0
                end if
             end do
          end do

       end do
    end do

    call mpi_allreduce(invM, buffer, size(invM), MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
    invM = buffer

    if (self%x%info%myid == 0) then
       do i = 0, 3
          write(*,*) real(invM(i,0:3),sp)
       end do
    end if

    call invert_matrix(invM, cholesky=.true.)

    if (self%x%info%myid == 0) then
       do i = 0, 3
!          write(*,*) real(invM(i,0:3),sp)
       end do
    end if

!!$    call mpi_finalize(ierr)
!!$    stop

    ! Store matrix rows
    q = 0
    do l = 0, self%lmax_pre_lowl
       do m = -l, l
          call self%x%info%lm2i(l,m,j)
          if (j >= 0) then
             k = l**2 + l + m
             do i = 0, nalm-1
                self%invM_lowl(i,q) = invM(i,k)
             end do
             q = q+1
          end if
       end do
    end do

    deallocate(invM, buffer)
    call map%dealloc()
    call tot%dealloc()

    call wall_time(t2)
    if (info%myid == 0) write(*,*) '  Low-ell init = ', t2-t1


  end subroutine updateLowlPrecond

  subroutine applyLowlPrecond(self, alm, alm0)
    implicit none
    class(comm_diffuse_comp),                   intent(in)     :: self
    real(dp),                 dimension(0:,1:), intent(inout)  :: alm
    real(dp),                 dimension(0:,1:), intent(in)     :: alm0

    real(dp)                  :: t1, t2
    integer(i4b)              :: i, j, k, l, m, q, myid, nalm, ntot, ierr
    real(dp), allocatable, dimension(:) :: y, yloc

    ntot = (self%lmax_pre_lowl+1)**2
    allocate(y(0:ntot-1))
    if (allocated(self%invM_lowl)) then
       nalm = size(self%invM_lowl,2)
       allocate(yloc(0:nalm-1))
       q   = 0
       do l = 0, self%lmax_pre_lowl
          do m = -l, l
             call self%x%info%lm2i(l,m,j)
             if (j >= 0) then
                yloc(q) = alm(j,1)
                q       = q+1
             end if
          end do
       end do
       y    = matmul(self%invM_lowl,yloc)
       deallocate(yloc)
    else
       nalm = 0
       y    = 0.d0
    end if
    call mpi_allreduce(MPI_IN_PLACE, y, ntot, MPI_DOUBLE_PRECISION, MPI_SUM, &
         & self%x%info%comm, ierr)

    alm = alm0
    do l = 0, self%lmax_pre_lowl
       do m = -l, l
          call self%x%info%lm2i(l,m,j)
          if (j >= 0) then
             k = l**2 + l + m
             alm(j,1) = y(k)
          end if
       end do
    end do


    deallocate(y)

  end subroutine applyLowlPrecond


  subroutine updateDeflatePrecond(self)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self

    real(dp)                  :: t1, t2, norm
    integer(i4b)              :: i, j, k, l, m, p, q, lp, mp, myid, nalm, ierr, nmaps, nmax, ndef, nside
    logical(lgt)              :: add
    class(comm_map),     pointer :: map => null(), map2 => null(), tot => null()
    class(comm_mapinfo), pointer :: info => null()
    real(dp),        allocatable, dimension(:,:) :: invM, buffer, Z
    real(dp),        allocatable, dimension(:) :: W

    if (self%lmax_def < 0) return

    call wall_time(t1)

    if (.not. allocated(self%Z_def)) then
       nside     = min(self%nside,128)
       self%ndef = 0
       q         = (nside/self%nside_def)**2
       info => comm_mapinfo(self%comm, nside, self%lmax_def, 1, .false.)
       map  => comm_map(info)

       nalm      = info%nalm
       nmax      = min(10000, info%npix)
       allocate(Z(0:nalm-1,nmax))

       call setup_needlets(info, self%nside_def, self%defmask, Z, self%ndef)

!!$       do i = 0, 12*self%nside_def**2-1
!!$          if (self%myid == 0 .and. mod(i,100) == 0) write(*,*) i, 12*self%nside_def**2-1
!!$          add = .false.
!!$          do j = 1, self%defmask%info%np
!!$             if (self%defmask%info%pix(j) == i) add = (self%defmask%map(j-1,1) == 0.d0)  
!!$          end do
!!$          call mpi_allreduce(MPI_IN_PLACE, add, 1, MPI_LOGICAL, MPI_LOR, info%comm, ierr)
!!$          if (.not. add) cycle
!!$
!!$          map%map = 0.d0
!!$          p       = i*q 
!!$          call nest2ring(nside, p, j)
!!$          do p = 1, info%np
!!$             if (info%pix(p) == j) map%map(p-1,1) = 1.d0 
!!$          end do
!!$          call map%smooth(self%fwhm_def)
!!$          norm = maxval(map%map(:,1))
!!$          call mpi_allreduce(MPI_IN_PLACE, norm, 1, MPI_DOUBLE_PRECISION, MPI_MAX, info%comm, ierr)
!!$          map%map(:,1) = map%map(:,1) / norm
!!$
!!$          self%ndef = self%ndef + 1
!!$          if (self%ndef > nmax) then
!!$             write(*,*) ' ERROR: Too many deflated basis functions'
!!$             stop
!!$          end if
!!$          Z(:,self%ndef) = map%alm(:,1)
!!$         
!!$       end do

!!$       Z = 0.d0
!!$       do l = 0, 30
!!$          do m = -l, l
!!$             self%ndef = self%ndef + 1
!!$             call map%info%lm2i(l,m,i)
!!$             if (i > -1) Z(i,self%ndef) = 1.d0
!!$          end do
!!$       end do


!!$       self%ndef = 4
!!$       Z = 0.d0
!!$       do p = 0, info%nalm-1
!!$          if (info%lm(1,p) == 0) Z(p,1) = 1.d0
!!$       end do
!!$       do p = 0, info%nalm-1
!!$          if (info%lm(1,p) == 1 .and. info%lm(2,p) == -1) Z(p,2) = 1.d0
!!$       end do
!!$       do p = 0, info%nalm-1
!!$          if (info%lm(1,p) == 1 .and. info%lm(2,p) == 0) Z(p,3) = 1.d0
!!$       end do
!!$       do p = 0, info%nalm-1
!!$          if (info%lm(1,p) == 1 .and. info%lm(2,p) == 1) Z(p,4) = 1.d0
!!$       end do
!!$       do p = 0, info%nalm-1
!!$          if (info%lm(1,p) == 2 .and. info%lm(2,p) == -2) Z(p,5) = 1.d0
!!$       end do
!!$       do p = 0, info%nalm-1
!!$          if (info%lm(1,p) == 2 .and. info%lm(2,p) == -1) Z(p,6) = 1.d0
!!$       end do
!!$       do p = 0, info%nalm-1
!!$          if (info%lm(1,p) == 2 .and. info%lm(2,p) == 0) Z(p,7) = 1.d0
!!$       end do
!!$       do p = 0, info%nalm-1
!!$          if (info%lm(1,p) == 2 .and. info%lm(2,p) == 1) Z(p,8) = 1.d0
!!$       end do
!!$       do p = 0, info%nalm-1
!!$          if (info%lm(1,p) == 2 .and. info%lm(2,p) == 2) Z(p,9) = 1.d0
!!$       end do
 
       if (self%myid == 0) write(*,*) 'nz = ', self%ndef

       allocate(self%Z_def(0:nalm-1,self%ndef))
       self%Z_def = Z(:,1:self%ndef)

       if (self%myid == 0) allocate(self%invM_def(self%ndef,self%ndef))

       call map%dealloc()
       deallocate(Z)
    end if

    ndef = size(self%Z_def,2)
    info => comm_mapinfo(self%comm, 2, self%lmax_def, self%nmaps, self%nmaps==3)
    map  => comm_map(info)
    tot  => comm_map(info)
    allocate(invM(ndef,ndef), buffer(ndef,ndef))
    do l = 1, ndef

       if (info%myid == 0 .and. mod(l,100) == 1) write(*,*) '  Precomputing deflate mode = ', l, ndef, trim(self%label)
       ! Set up basis vector
       map%alm = 0.d0
       map%alm(:,1) = self%Z_def(:,l)

       call self%Cl%sqrtS(alm=map%alm, info=map%info) ! Multiply with sqrt(Cl)

       ! Add frequency dependent terms
       tot%alm = 0.d0
       do i = 1, numband

          ! Compute component-summed map, ie., column-wise matrix elements
          nmaps  = data(i)%info%nmaps
          info     => comm_mapinfo(self%comm, data(i)%N%nside_chisq_lowres, self%lmax_def, &
               & nmaps, nmaps==3)
          map2     => comm_map(info)
          do k = 1, min(self%nmaps, nmaps)
             map2%alm(:,k) = map%alm(:,k) * self%F_mean(i,0,k)
          end do
          call data(i)%B(0)%p%conv(trans=.false., map=map2)
          call map2%Y()

          ! Multiply with invN
          call data(i)%N%InvN_lowres(map2)

             ! Project summed map into components, ie., row-wise matrix elements
          call map2%Yt()
          call data(i)%B(0)%p%conv(trans=.true., map=map2)
          map2%alm(:,1) = map2%alm(:,1) * self%F_mean(i,0,1)

          ! Add up alms
          do k = 1, min(self%nmaps, nmaps)
             tot%alm(:,k) = tot%alm(:,k) + map2%alm(:,k)
          end do

          call map2%dealloc()
       end do

       ! Add prior term and multiply with sqrt(S) for relevant components
       call self%Cl%sqrtS(alm=tot%alm, info=tot%info)

       ! Add (unity) prior term, and store column
       tot%alm(:,1) = tot%alm(:,1) + self%Z_def(:,l)

       if (info%nalm > 0) then
          do j = 1, ndef
             invM(l,j) = sum(self%Z_def(:,j)*tot%alm(:,1))
          end do
!          write(*,*) self%myid, l, sum(invM(l,:)), sum(tot%alm(:,1))
       else
          invM(l,:) = 0.d0
       end if
          
    end do

    call mpi_reduce(invM, buffer, size(invM), MPI_DOUBLE_PRECISION, MPI_SUM, 0, info%comm, ierr)

    if (self%myid == 0) then
       self%invM_def = buffer
       allocate(W(ndef))
       call get_eigenvalues(buffer, W)
       write(*,*) 'W =', real(W,sp)
       deallocate(W)
       call invert_matrix(self%invM_def, cholesky=.true.)

!!$       do i= 1, self%ndef
!!$          write(*,*) real(self%invM_def(i,:),sp)
!!$       end do
    end if

!!$call mpi_finalize(ierr)
!!$stop

    deallocate(invM, buffer)
    call map%dealloc()
    call tot%dealloc()

    call wall_time(t2)
    if (info%myid == 0) write(*,*) '  Deflate init = ', t2-t1

  end subroutine updateDeflatePrecond



  subroutine applyDeflatePrecond(self, alm, Qalm)
    implicit none
    class(comm_diffuse_comp),                   intent(in)     :: self
    real(dp),                 dimension(0:,1:), intent(in)     :: alm
    real(dp),                 dimension(0:,1:), intent(out)    :: Qalm

    real(dp)                  :: t1, t2
    integer(i4b)              :: i, j, k, l, m, myid, nalm, ntot, ierr, ndef
    real(dp), allocatable, dimension(:) :: y, ytot
    class(comm_map),     pointer :: map => null(), map2 => null(), tot => null()
    class(comm_mapinfo), pointer :: info => null()

    Qalm= 0.d0
    if (self%lmax_def < 0) return

    ndef = size(self%Z_def,2)
    allocate(y(ndef), ytot(ndef))

    ! Multiply with Z
    info => comm_mapinfo(self%comm, 2, self%lmax_def, 1, .false.)
    map  => comm_map(info)
    map%alm = 0.d0
    do i = 0, info%nalm-1
       call info%i2lm(i,l,m)
       call self%x%info%lm2i(l,m,j)
       if (j > -1) map%alm(i,1) = alm(j,1)
    end do
    if (info%nalm > 0) then
       y = matmul(transpose(self%Z_def), map%alm(:,1))
    else
       y = 0.d0
    end if
    call mpi_reduce(y, ytot, ndef, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
         & self%x%info%comm, ierr)

    ! Multiply with invE
    if (self%myid == 0) then
       y = matmul(self%invM_def, ytot)
    end if
    call mpi_bcast(y, ndef, MPI_DOUBLE_PRECISION, 0, self%x%info%comm, ierr)

    ! Multiply with Z^t
    map%alm(:,1) = matmul(self%Z_def, y)
    
    do i = 0, info%nalm-1
       call info%i2lm(i,l,m)
       call self%x%info%lm2i(l,m,j)
       if (j > -1) Qalm(j,1) = map%alm(i,1)
    end do
    
    call map%dealloc()
    deallocate(y, ytot)

  end subroutine applyDeflatePrecond

  subroutine setup_needlets(info, nside_def, defmask, Z, ndef)
    implicit none
    class(comm_mapinfo), pointer,          intent(in)  :: info
    integer(i4b),                          intent(in)  :: nside_def
    class(comm_map),                       intent(in)  :: defmask
    real(dp),            dimension(0:,1:), intent(out) :: Z
    integer(i4b),                          intent(out) :: ndef

    integer(i4b) :: i, j, l, p, q, nside, m, res, ierr
    real(dp)     :: B, bl, Bj
    character(len=4) :: ntext
    logical(lgt) :: add
    class(comm_map),     pointer :: map => null()
    real(dp), allocatable, dimension(:) :: x, t, f, psi, phi, b0
    type(spline_type) :: sb, spsi, sphi

    ! Set up needlet kernel
    m    = 10000
    B    = 2.72d0
    allocate(x(m), t(m), f(m), psi(m), phi(m), b0(m))
    do i = 1, m
       x(i) = -1.d0 + 2.d0*(i-1.d0)/(m-1.d0)
    end do
    f = exp(-1.d0/(1.d0-x**2))
!!$    if (info%myid == 0) then
!!$       open(58,file='f.dat')
!!$       do i = 1, m
!!$          write(58,*) x(i), f(i)
!!$       end do
!!$       close(58)
!!$    end if

    psi(1) = 0.d0
    do i = 2, m
       psi(i) = psi(i-1) + 0.5d0*(f(i)+f(i-1))*(x(i)-x(i-1))
    end do
    psi = psi / psi(m)
    call spline(spsi, x, psi)
!!$    if (info%myid == 0) then
!!$       open(58,file='psi.dat')
!!$       do i = 1, m
!!$          write(58,*) x(i), psi(i)
!!$       end do
!!$       close(58)
!!$    end if

    do i = 1, m
       t(i) = (i-1.d0)/(m-1.d0)
       if (t(i) < 1.d0/B) then
          phi(i) = 1.d0
       else
          phi(i) = splint(spsi, 1.d0-2.d0*B/(B-1.d0)*(t(i)-1.d0/B))
       end if
    end do
!!$    if (info%myid == 0) then
!!$       open(58,file='phi.dat')
!!$       do i = 1, m
!!$          write(58,*) t(i), phi(i)
!!$       end do
!!$       close(58)
!!$    end if
    call spline(sphi, t, phi)

    do i = 1, m
       t(i)  = B*(i-1.d0)/(m-1.d0)
       if (t(i) < 1.d0) then
          b0(i) = sqrt(splint(sphi, t(i)/B) - splint(sphi, t(i)))
          !b0(i) = (splint(sphi, t(i)/B) - splint(sphi, t(i)))
       else
          b0(i) = sqrt(splint(sphi, t(i)/B))
          !b0(i) = (splint(sphi, t(i)/B))
       end if
       if (b0(i) < 1d-10 .or. b0(i) /= b0(i)) b0(i) = 0.d0
    end do
    call spline(sb, t, b0)
!!$    if (info%myid == 0) then
!!$       open(58,file='b.dat')
!!$       do i = 1, m
!!$          write(58,*) t(i), b0(i)
!!$       end do
!!$       close(58)
!!$    end if


!!$    if (info%myid == 0)then
!!$       do i = 1, m
!!$          write(*,*) i, t(i), b0(i)
!!$       end do

!!$       open(58,file='needlet.dat')
!!$       do i = 0, info%lmax
!!$          write(*,*)  i, splint(sb, i/B**5)
!!$          write(58,*) i, splint(sb, i/B**5)
!!$       end do
!!$       close(58)
!!$    end if
!!$    call mpi_finalize(ierr)
!!$    stop

    map   => comm_map(info)
    ndef  = 0.d0
    nside = 1
    Bj    = B
    do while(nside <= nside_def)
       do i = 0, 12*nside**2-1
          p = i*(defmask%info%nside/nside)**2
          call nest2ring(defmask%info%nside, p, j)
          add = .false.
          do p = 1, defmask%info%np
             if (defmask%info%pix(p) == j) add = (defmask%map(p-1,1) == 0.d0)  
          end do
          call mpi_allreduce(MPI_IN_PLACE, add, 1, MPI_LOGICAL, MPI_LOR, info%comm, ierr)
          if (.not. add) cycle
          
          ndef    = ndef+1
          if (ndef > size(Z,2)) then
             write(*,*) ' ERROR: Too many deflated basis functions'
             stop
          end if

          map%map = 0.d0
          p       = i*(info%nside/nside)**2
          call nest2ring(info%nside, p, j)
          do p = 1, info%np
             if (info%pix(p) == j) map%map(p,1) = 1.d0
          end do
          call map%YtW

          do j = 0, info%nalm-1
             l = info%lm(1,j)
             if (l > Bj) then
                bl = 0.d0
             else
                bl = splint(sb, l/Bj)
             end if
             map%alm(j,1) = bl * map%alm(j,1) 
          end do

          if (mod(ndef,10) == 0) then
             call int2string(ndef, ntext)
             call map%Y
             call map%writeFITS('need'//ntext//'.fits')
          end if

          Z(:,ndef) = map%alm(:,1)
       end do
       if (info%myid == 0) write(*,*) nside, ndef

       nside = 2*nside
       Bj    = B*Bj
    end do

    ! Add low-l modes
!!$    do l = 0, 10
!!$       do m = -l, l
!!$          ndef = ndef + 1
!!$          call map%info%lm2i(l,m,i)
!!$          if (i > -1) Z(i,ndef) = 1.d0
!!$       end do
!!$    end do


    call free_spline(sb)
    call free_spline(spsi)
    call free_spline(sphi)
    call map%dealloc()
    deallocate(x, t, f, psi, phi, b0)

  end subroutine setup_needlets

  subroutine applyMonoDipolePrior(self)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self

    integer(i4b) :: i, j, k, l, m, ierr
    real(dp)     :: mu(0:3), a, b, Amat(0:3,0:3), bmat(0:3), v(0:3)
    class(comm_map), pointer :: map 

    if (trim(self%mono_prior_type) == 'none') then ! No active monopole prior
       return
    end if

    ! Compute monopole offset given the specified prior
    mu = 0.d0

    ! Generate real-space component map
    map => comm_map(self%x)
    call self%B_out%conv(trans=.false., map=map)
    call map%Y
!!$    call map%writeFITS('ff.fits')
!!$    call mpi_finalize(ierr)
!!$    stop

    if (trim(self%mono_prior_type) == 'monopole') then        ! Set monopole to zero outside user-specified mask

       ! Compute mean outside mask; no noise weighting for now at least
       a = sum(map%map(:,1) * self%mono_prior_map%map(:,1))
       b = sum(self%mono_prior_map%map(:,1))
       call mpi_allreduce(MPI_IN_PLACE, a, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, b, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       mu(0) = a / b

       ! Subtract mean in real space
       self%x%map(:,1) = self%x%map(:,1) - mu(0)

       if (self%x%info%myid == 0) write(*,fmt='(a,f10.3)') '   Monopole prior correction for '//trim(self%label)//': ', mu(0)

    else if (trim(self%mono_prior_type) == 'monopole+dipole') then        ! Set monopole and dipole to zero outside user-specified mask
       ! Generate real-space component map

       ! Compute mean outside mask; no noise weighting for now at least
       Amat = 0.d0; bmat = 0.d0
       do i = 0, self%x%info%np-1
          if (self%mono_prior_map%map(i,1) < 0.5d0) cycle
          v(0) = 1.d0
          call pix2vec_ring(self%x%info%nside, self%x%info%pix(i+1), v(1:3))
          do j = 0, 3
             do k = 0, 3
                Amat(j,k) = Amat(j,k) + v(j)*v(k)
             end do
             bmat(j) = bmat(j) + v(j)*map%map(i,1)
          end do
       end do
       call mpi_allreduce(MPI_IN_PLACE, Amat, 16, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, bmat,  4, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       call solve_system_real(Amat, mu, bmat)

       ! Subtract mean and dipole in real space
       do i = 0, self%x%info%np-1
          v(0) = 1.d0
          call pix2vec_ring(self%x%info%nside, self%x%info%pix(i+1), v(1:3))
          self%x%map(i,1) = self%x%map(i,1) - sum(v*mu)
       end do

       if (self%x%info%myid == 0) write(*,fmt='(a,4f10.3)') '   Monopole prior correction for '//trim(self%label)//': ', mu

    else if (trim(self%mono_prior_type) == 'crosscorr') then ! Enforce zero intercept in correlation with specified map
       write(*,*) 'Error: Cross-correlation monopole prior not implemented yet'
       stop
    end if
    
    ! Subtract mean in harmonic space
    do i = 0, self%x%info%nalm-1
       if (self%x%info%lm(1,i) == 0 .and. self%x%info%lm(2,i) == 0) then
          self%x%alm(i,1) = self%x%alm(i,1) - mu(0) * sqrt(4.d0*pi)
       end if
       if (self%x%info%lm(1,i) == 1 .and. self%x%info%lm(2,i) == -1) then
          self%x%alm(i,1) = self%x%alm(i,1) - mu(2) * sqrt(4.d0*pi/3.d0)
       end if
       if (self%x%info%lm(1,i) == 1 .and. self%x%info%lm(2,i) == 0) then
          self%x%alm(i,1) = self%x%alm(i,1) - mu(3) * sqrt(4.d0*pi/3.d0)
       end if
       if (self%x%info%lm(1,i) == 1 .and. self%x%info%lm(2,i) == 1) then
          self%x%alm(i,1) = self%x%alm(i,1) + mu(1) * sqrt(4.d0*pi/3.d0)
       end if
    end do

    call map%dealloc()
    
  end subroutine applyMonoDipolePrior

  subroutine nullify_monopole_amp(band)
    implicit none
    character(len=*), intent(in) :: band

    integer(i4b) :: i
    class(comm_comp), pointer :: c => null()

    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%label) == trim(band)) then
             do i = 0, c%x%info%nalm-1
                if (c%x%info%lm(1,i) == 0 .and. c%x%info%lm(2,i) == 0) then
                   c%x%alm(i,1) = 0.d0
                   return
                end if
             end do
          end if
       end select
       c => c%next()
    end do

  end subroutine nullify_monopole_amp


end module comm_diffuse_comp_mod
