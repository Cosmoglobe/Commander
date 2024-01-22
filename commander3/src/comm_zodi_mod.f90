module comm_zodi_mod
   use comm_utils
   use comm_param_mod
   use spline_1D_mod
   use comm_hdf_mod
   use hashtbl

   implicit none

   private
   public initialize_zodi_mod, get_s_zodi, ZodiModel, ZodiComponent, zodi_model


   type, abstract :: ZodiComponent
      real(dp) :: x_0, y_0, z_0, incl, Omega, n_0
      real(dp), allocatable :: sin_omega, cos_omega, sin_incl, cos_incl
   contains
      procedure :: init_base_comp
      procedure(init_interface),    deferred :: init
      procedure(density_interface), deferred :: get_density
      procedure(prior_interface),   deferred :: init_priors_and_scales
      procedure(p2m_interface),     deferred :: param2model
      procedure(m2p_interface),     deferred :: model2param
   end type ZodiComponent

   type :: ZodiComponentContainer
      integer(i4b) :: npar, start_ind
      class(ZodiComponent), allocatable :: c
      character(len=32), allocatable :: labels(:)
   end type ZodiComponentContainer

   abstract interface
      subroutine init_interface(self)
         import dp, ZodiComponent
         class(ZodiComponent)  :: self
      end subroutine init_interface

      subroutine density_interface(self, X_vec, theta, n_out)
         ! Returns the dust density (n) of the component at heliocentric
         ! coordinates (x, y, z) and the earths longitude (theta).
         import i4b, dp, ZodiComponent
         class(ZodiComponent) :: self
         real(dp), intent(in) :: X_vec(:, :)
         real(dp), intent(in) :: theta
         real(dp), intent(out):: n_out(:)
         real(dp) :: x_prime, y_prime, z_prime, R, Z_midplane
      end subroutine density_interface

      subroutine prior_interface(self, start_ind, prior, scale)
        import dp, i4b, ZodiComponent
        class(ZodiComponent),       intent(in)    :: self
        integer(i4b),               intent(in)    :: start_ind
        real(dp), dimension(1:,1:), intent(inout) :: prior
        real(dp), dimension(1:),    intent(inout) :: scale
      end subroutine prior_interface

      subroutine p2m_interface(self, x)
        import dp,ZodiComponent
        class(ZodiComponent),       intent(inout) :: self
        real(dp), dimension(1:),    intent(in)    :: x
      end subroutine p2m_interface

      subroutine m2p_interface(self, x)
        import dp,ZodiComponent
        class(ZodiComponent),       intent(in)  :: self
        real(dp), dimension(1:),    intent(out) :: x
      end subroutine m2p_interface

   end interface

   type, extends(ZodiComponent) :: ZodiCloud
      real(dp) :: alpha, beta, gamma, mu
   contains
      procedure :: init => init_cloud
      procedure :: get_density => get_density_cloud
      procedure :: init_priors_and_scales => init_cloud_priors_and_scales
      procedure :: param2model => param2model_cloud
      procedure :: model2param => model2param_cloud
   end type ZodiCloud

   type, extends(ZodiComponent) :: ZodiBand
      real(dp) :: delta_zeta, delta_r, v, p
      real(dp) :: delta_zeta_rad = 0.d0
   contains
      procedure :: init => init_band
      procedure :: get_density => get_density_band
      procedure :: init_priors_and_scales => init_band_priors_and_scales
      procedure :: param2model => param2model_band
      procedure :: model2param => model2param_band
   end type ZodiBand

   type, extends(ZodiComponent) :: ZodiRing
      real(dp) :: R_0, sigma_r, sigma_z, theta_0, sigma_theta
      real(dp) :: theta_0_rad = 0.d0
      real(dp) :: sigma_theta_rad = 0.d0
   contains
      procedure :: init => init_ring
      procedure :: get_density => get_density_ring
      procedure :: init_priors_and_scales => init_ring_priors_and_scales
      procedure :: param2model => param2model_ring
      procedure :: model2param => model2param_ring
   end type ZodiRing

   type, extends(ZodiComponent) :: ZodiFeature
      real(dp) :: R_0, sigma_r, sigma_z, theta_0, sigma_theta
      real(dp) :: theta_0_rad = 0.d0
      real(dp) :: sigma_theta_rad = 0.d0
   contains
      procedure :: init => init_feature
      procedure :: get_density => get_density_feature
      procedure :: init_priors_and_scales => init_feature_priors_and_scales
      procedure :: param2model => param2model_feature
      procedure :: model2param => model2param_feature
   end type ZodiFeature

   type, extends(ZodiComponent) :: ZodiInterstellar
      real(dp)  :: R, alpha
   contains
      procedure :: init => init_interstellar
      procedure :: get_density => get_density_interstellar
      procedure :: init_priors_and_scales => init_interstellar_priors_and_scales
      procedure :: param2model => param2model_interstellar
      procedure :: model2param => model2param_interstellar
   end type ZodiInterstellar

   type, extends(ZodiComponent) :: ZodiFan
      real(dp) :: Q, P, gamma, Z_midplane_0, R_outer
   contains
      procedure :: init => init_fan
      procedure :: get_density => get_density_fan
      procedure :: init_priors_and_scales => init_fan_priors_and_scales
      procedure :: param2model => param2model_fan
      procedure :: model2param => model2param_fan
   end type ZodiFan

   type, extends(ZodiComponent) :: ZodiComet
      real(dp) :: P, Z_midplane_0, R_inner, R_outer
   contains
      procedure :: init => init_comet
      procedure :: get_density => get_density_comet
      procedure :: init_priors_and_scales => init_comet_priors_and_scales
      procedure :: param2model => param2model_comet
      procedure :: model2param => model2param_comet
   end type ZodiComet

   type :: ZodiModel
      class(ZodiComponentContainer), allocatable :: comps(:)
      character(len=24), allocatable :: comp_labels(:), general_labels(:), par_labels(:)
      integer(i4b) :: n_comps, n_params, n_common_params, n_general_params
      logical(lgt) :: joint_mono
      real(dp)     :: min_solar_elong, max_solar_elong
      real(dp)     :: T_0, delta
      real(dp), dimension(10) :: F_sun = [2.3405606d8, 1.2309874d8, 64292872d0, 35733824d0, 5763843d0, 1327989.4d0, 230553.73d0, 82999.336d0, 42346.605d0, 14409.608d0] * 1d-20 ! convert to specific intensity units
      real(dp), dimension(10) :: C0 = [-0.94209999, -0.52670002, -0.4312, 0., 0., 0., 0., 0., 0., 0.]
      real(dp), dimension(10) :: C1 = [0.1214, 0.18719999, 0.1715, 0., 0., 0., 0., 0., 0., 0.]
      real(dp), dimension(10) :: C2 = [-0.1648, -0.59829998, -0.63330001, 0., 0., 0., 0., 0., 0., 0.]

      integer(i4b) :: npar_tot
      integer(i4b), allocatable, dimension(:,:) :: theta_stat
      logical(lgt), allocatable, dimension(:,:) :: sampgroup_active_band
      integer(i4b), allocatable, dimension(:)   :: theta2band
      real(dp),     allocatable, dimension(:,:) :: theta_prior
      real(dp),     allocatable, dimension(:)   :: theta_scale
    contains
      procedure :: init_comps, init_general_params, model_to_chain, params_to_model2, model_to_params2, comp_from_chain, get_par_ind, init_general_priors_and_scales
   end type ZodiModel

   type(ZodiModel), target :: zodi_model
   integer(i4b)            :: numband
   character(len=128)      :: zodi_refband
   character(len=128), allocatable, dimension(:) :: band_labels
   character(len=128), allocatable, dimension(:) :: band_instlabels
   character(len=128), allocatable, dimension(:) :: band_todtype
   real(dp),           allocatable, dimension(:) :: band_nu_c
   
contains
   subroutine initialize_zodi_mod(cpar)
      type(comm_params), intent(in) :: cpar
      integer(i4b) :: i, j, k, ierr, ind, npar, ntok
      character(len=256) :: file_path
      real(dp), allocatable :: comp_params(:, :)
      character(len=128), allocatable :: comp_labels(:)
      character(len=128) :: tokens(100)

      ! Find number of bands and labels
      numband         = count(cpar%ds_active)
      allocate(band_labels(numband),band_instlabels(numband),band_todtype(numband),band_nu_c(numband))
      band_labels     = pack(cpar%ds_label, cpar%ds_active)
      band_instlabels = pack(cpar%ds_instlabel, cpar%ds_active)
      band_todtype    = pack(cpar%ds_tod_type, cpar%ds_active)
      band_nu_c       = pack(cpar%ds_nu_c, cpar%ds_active)
      zodi_refband    = cpar%zs_refband
      
      ! Set model and zodi_mod parameters from cpar
      zodi_model%n_comps = cpar%zs_ncomps
      allocate(comp_params(zodi_model%n_comps, size(cpar%zs_comp_params, dim=1)))
      zodi_model%general_labels   = cpar%zodi_param_labels%general
      zodi_model%comp_labels      = cpar%zs_comp_labels(1:zodi_model%n_comps)
      zodi_model%n_common_params  = size(cpar%zodi_param_labels%common)
      zodi_model%n_general_params = size(cpar%zodi_param_labels%general)
      zodi_model%joint_mono       = cpar%zs_joint_mono
      
      comp_params = cpar%zs_comp_params(:, :, 1)
      do i = 1, zodi_model%n_comps
         if (trim(adjustl(cpar%zs_init_hdf(i))) /= 'none') then 
            call zodi_model%comp_from_chain(cpar, comp_params, i)
         end if
      end do
      call zodi_model%init_comps(comp_params, cpar%zs_comp_types, cpar%zodi_param_labels)
      call zodi_model%init_general_params(cpar%zs_general_params(:, 1))

      ! Find total number of free parameters, and set up parameter mapping
      zodi_model%npar_tot = zodi_model%n_general_params + 2*zodi_model%n_comps*numband + numband
      zodi_model%comps(1)%start_ind = zodi_model%n_general_params + 1
      do i = 1, zodi_model%n_comps
         zodi_model%npar_tot = zodi_model%npar_tot + zodi_model%comps(i)%npar
         if (i < zodi_model%n_comps) then
            zodi_model%comps(i+1)%start_ind = zodi_model%comps(i)%start_ind + &
                 & zodi_model%comps(i)%npar + 2*numband
         end if
      end do

      ! stat =  0  -> sample freely
      ! stat = -1  -> fix to input
      ! stat = -2  -> fix to zero
      ! stat = -3  -> fix to unity
      ! stat >  0  -> set equal to parameter stat
      allocate(zodi_model%theta_stat(zodi_model%npar_tot,0:cpar%zs_num_samp_groups))
      allocate(zodi_model%theta2band(zodi_model%npar_tot))
      allocate(zodi_model%theta_prior(4,zodi_model%npar_tot)) ! [min,max,mean,rms]
      allocate(zodi_model%theta_scale(zodi_model%npar_tot))
      allocate(zodi_model%par_labels(zodi_model%npar_tot))
      
      ! Set up sampling groups
      allocate(zodi_model%sampgroup_active_band(numband,cpar%zs_num_samp_groups))
      zodi_model%sampgroup_active_band = .false.
      do i = 1, cpar%zs_num_samp_groups
         call get_tokens(cpar%zs_samp_group_bands(i), ',', tokens, ntok) 
         do j = 1, ntok
            k = get_string_index(band_labels, tokens(j))
            zodi_model%sampgroup_active_band(k,i) = .true.
         end do
         call samp_group2stat(cpar, i, zodi_model%sampgroup_active_band(:,i), zodi_model%theta_stat(:,i))
      end do
      do i = 1, zodi_model%npar_tot
         zodi_model%theta_stat(i,0) = maxval(zodi_model%theta_stat(i,1:cpar%zs_num_samp_groups))
      end do

      ! Initialize parameter-band mapping
      zodi_model%theta2band(1:zodi_model%n_general_params) = 0 ! General paraneters affect all bands
      do i = 1, zodi_model%n_comps
         ind  = zodi_model%comps(i)%start_ind
         npar = zodi_model%comps(i)%npar
         zodi_model%theta2band(ind:ind+npar-1) = 0 ! Shape paraneters affect all band
         do j = 1, numband
            zodi_model%theta2band(ind+npar+j-1) = j   ! Emissivity only affect band j
            zodi_model%theta2band(ind+npar+numband+j-1) = j ! The same for albedo
         end do
      end do

      ! Initialize priors and scale factors
      call zodi_model%init_general_priors_and_scales(zodi_model%theta_prior, &
           & zodi_model%theta_scale)
      zodi_model%par_labels(1:zodi_model%n_general_params) = &
           & zodi_model%general_labels
      do i = 1, zodi_model%n_comps
         ! Shape parameters
         ind = zodi_model%comps(i)%start_ind
         call zodi_model%comps(i)%c%init_priors_and_scales(ind, &
              & zodi_model%theta_prior, zodi_model%theta_scale)
         zodi_model%par_labels(ind:ind+zodi_model%comps(i)%npar-1) = &
              & zodi_model%comps(i)%labels

         ! Emissivity and albedo
         ind = zodi_model%comps(i)%start_ind + zodi_model%comps(i)%npar-1
         do j = 1, numband
            zodi_model%theta_prior(:,ind+j) = [0.d0, 5.d0, 1.d0, -1.d0] ! Emissivity
            zodi_model%theta_scale(ind+j)   = 1.d0
            zodi_model%par_labels(ind+j)    = 'em@'//trim(band_labels(j))

            zodi_model%theta_prior(:,ind+numband+j) = [0.d0, 1.d0, 0.3d0, -1.d0] ! Albedo
            zodi_model%theta_scale(ind+numband+j)   = 1.d0
            zodi_model%par_labels(ind+numband+j)    = 'al@'//trim(band_labels(j))
         end do
      end do
      ! Monopoles
      do j = 1, numband
         ind = zodi_model%npar_tot - numband + j
         zodi_model%theta_prior(:,ind) = [0.d0, 1d30, 0.d0, -1.d0] ! Priors
         zodi_model%theta_scale(ind)   = 1.d0
         zodi_model%par_labels(ind)    = 'm@'//trim(band_labels(j))
      end do
      
      if (cpar%myid_chain == 0) then
         write(*,*) ' Total number of free zodi parameters = ', count(zodi_model%theta_stat(:,0)==0)
!!$         do i = 1, zodi_model%npar_tot
!!$            write(*,*) i,  ', stat=', zodi_model%theta_stat(i,:)
!!$         end do
      end if

!!$      call mpi_finalize(ierr)
!!$      stop
!!$      
    end subroutine initialize_zodi_mod

   subroutine init_general_params(self, general_params)
      class(ZodiModel), intent(inout) :: self
      real(dp), intent(in) :: general_params(:)
      self%T_0 = general_params(1)
      self%delta = general_params(2)
   end subroutine

   subroutine init_base_comp(self)
     class(ZodiComponent) :: self
      self%sin_omega = sin(self%Omega*deg2rad)
      self%cos_omega = cos(self%Omega*deg2rad)
      self%sin_incl  = sin(self%incl*deg2rad)
      self%cos_incl  = cos(self%incl*deg2rad)
   end subroutine init_base_comp

   subroutine init_cloud(self)
      class(ZodiCloud) :: self
      call self%init_base_comp()
   end subroutine init_cloud

   subroutine init_band(self)
     class(ZodiBand) :: self
      self%delta_zeta_rad = self%delta_zeta*deg2rad
      call self%init_base_comp()
   end subroutine init_band

   subroutine init_ring(self)
      class(ZodiRing) :: self
      self%theta_0_rad = self%theta_0*deg2rad
      self%sigma_theta_rad = self%sigma_theta*deg2rad
      call self%init_base_comp()
   end subroutine init_ring

   subroutine init_feature(self)
      class(ZodiFeature) :: self
      self%theta_0_rad = self%theta_0*deg2rad
      self%sigma_theta_rad = self%sigma_theta*deg2rad
      call self%init_base_comp()
   end subroutine init_feature

   subroutine init_interstellar(self)
     class(ZodiInterstellar) :: self
     call self%init_base_comp()
   end subroutine init_interstellar

   subroutine init_fan(self)
      class(ZodiFan) :: self
      call self%init_base_comp()
    end subroutine init_fan

    subroutine init_comet(self)
      class(ZodiComet) :: self
      call self%init_base_comp()
    end subroutine init_comet

    subroutine init_general_priors_and_scales(self, prior, scale)
      implicit none
      class(ZodiModel),           intent(in)    :: self
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:),    intent(inout) :: scale

      prior(:,1) = [250.d0, 300.d0, 286.d0, 5.d0] ! T_0
      scale(1)   = 286.d0
      prior(:,2) = [0.4d0, 0.5d0, 0.467d0, 0.004d0] ! delta
      scale(2)   = 0.4d0      
    end subroutine init_general_priors_and_scales

    subroutine init_cloud_priors_and_scales(self, start_ind, prior, scale)
      implicit none
      class(ZodiCloud),           intent(in)    :: self
      integer(i4b),               intent(in)    :: start_ind
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:),    intent(inout) :: scale
      
      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0)   = 1.d-9
      prior(:,start_ind+1) = [-30.d0, 30.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1)   = 1.d0      
      prior(:,start_ind+2) = [-720.d0, 720.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2)   = 1.d0      
      prior(:,start_ind+3) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3)   = 1.d0      
      prior(:,start_ind+4) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4)   = 1.d0      
      prior(:,start_ind+5) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5)   = 1.d0      
      ! Component-specific parameters
      prior(:,start_ind+6) = [1.d0, 1.5d0, 1.34d0, -1.d0] ! alpha
      scale(start_ind+6)   = 1.d0      
      prior(:,start_ind+7) = [3.d0, 5d0, 4.14d0, -1.d0] ! beta
      scale(start_ind+7)   = 1.d0      
      prior(:,start_ind+8) = [0.3d0, 1.1d0, 0.942d0, -1.d0] ! gamma
      scale(start_ind+8)   = 1.d0      
      prior(:,start_ind+9) = [0.1d0, 0.4d0, 0.189d0, -1.d0] ! mu
      scale(start_ind+9)   = 1.d0      
    end subroutine init_cloud_priors_and_scales

    subroutine init_band_priors_and_scales(self, start_ind, prior, scale)
      implicit none
      class(ZodiBand),            intent(in)    :: self
      integer(i4b),               intent(in)    :: start_ind
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:),    intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0)   = 1.d-9
      prior(:,start_ind+1) = [-30.d0, 30.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1)   = 1.d0      
      prior(:,start_ind+2) = [-720.d0, 720.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2)   = 1.d0      
      prior(:,start_ind+3) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3)   = 1.d0      
      prior(:,start_ind+4) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4)   = 1.d0      
      prior(:,start_ind+5) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5)   = 1.d0      
      ! Component-specific parameters
      prior(:,start_ind+6) = [0.d0, 30d0, 0d0, -1.d0] ! delta_zeta
      scale(start_ind+6)   = 1.d0      
      prior(:,start_ind+7) = [0.8d0, 5.4d0, 4.14d0, -1.d0] ! delta_r
      scale(start_ind+7)   = 1.d0      
      prior(:,start_ind+8) = [0.01d0, 1.1d0, 0.942d0, -1.d0] ! v
      scale(start_ind+8)   = 1.d0      
      prior(:,start_ind+9) = [0.1d0, 5d0, 0.189d0, -1.d0] ! p
      scale(start_ind+9)   = 1.d0      
    end subroutine init_band_priors_and_scales

    subroutine init_ring_priors_and_scales(self, start_ind, prior, scale)
      implicit none
      class(ZodiRing),            intent(in)    :: self
      integer(i4b),               intent(in)    :: start_ind
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:),    intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0)   = 1.d-9
      prior(:,start_ind+1) = [-30.d0, 30.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1)   = 1.d0      
      prior(:,start_ind+2) = [-720.d0, 720.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2)   = 1.d0      
      prior(:,start_ind+3) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3)   = 1.d0      
      prior(:,start_ind+4) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4)   = 1.d0      
      prior(:,start_ind+5) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5)   = 1.d0      
      ! Component-specific parameters
      prior(:,start_ind+6) = [0.9d0, 1.1d0, 0d0, -1.d0] ! r
      scale(start_ind+6)   = 1.d0      
      prior(:,start_ind+7) = [0.d0, 0.3d0, 0.2d0, -1.d0] ! delta_r
      scale(start_ind+7)   = 1.d0      
      prior(:,start_ind+8) = [0.0d0, 0.2d0, 0.1d0, -1.d0] ! delta_z
      scale(start_ind+8)   = 1.d0      
      prior(:,start_ind+9) = [-60.d-3, 60.d-3, 0.d0, -1.d0] ! theta
      scale(start_ind+9)   = 1.d0      
      prior(:,start_ind+10) = [0.d0, 0.001d0, 0.d0, -1.d0] ! sigma_theta
      scale(start_ind+10)   = 1.d0      
    end subroutine init_ring_priors_and_scales

    subroutine init_feature_priors_and_scales(self, start_ind, prior, scale)
      implicit none
      class(ZodiFeature),         intent(in)    :: self
      integer(i4b),               intent(in)    :: start_ind
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:),    intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0)   = 1.d-9
      prior(:,start_ind+1) = [-30.d0, 30.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1)   = 1.d0      
      prior(:,start_ind+2) = [-720.d0, 720.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2)   = 1.d0      
      prior(:,start_ind+3) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3)   = 1.d0      
      prior(:,start_ind+4) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4)   = 1.d0      
      prior(:,start_ind+5) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5)   = 1.d0      
      ! Component-specific parameters
      prior(:,start_ind+6) = [0.9d0, 1.1d0, 0d0, -1.d0] ! r
      scale(start_ind+6)   = 1.d0      
      prior(:,start_ind+7) = [0.d0, 0.3d0, 0.2d0, -1.d0] ! delta_r
      scale(start_ind+7)   = 1.d0      
      prior(:,start_ind+8) = [0.0d0, 0.2d0, 0.1d0, -1.d0] ! delta_z
      scale(start_ind+8)   = 1.d0      
      prior(:,start_ind+9) = [-20.d0, 20.d0, 0.d0, -1.d0] ! theta
      scale(start_ind+9)   = 1.d0      
      prior(:,start_ind+10) = [0.d0, 30.d0, 0.d0, -1.d0] ! sigma_theta
      scale(start_ind+10)   = 1.d0      
    end subroutine init_feature_priors_and_scales

    subroutine init_interstellar_priors_and_scales(self, start_ind, prior, scale)
      implicit none
      class(ZodiInterstellar),    intent(in)    :: self
      integer(i4b),               intent(in)    :: start_ind
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:),    intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0)   = 1.d-9
      prior(:,start_ind+1) = [0.d0, 00.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1)   = 1.d0      
      prior(:,start_ind+2) = [0.d0, 0.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2)   = 1.d0      
      prior(:,start_ind+3) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3)   = 1.d0      
      prior(:,start_ind+4) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4)   = 1.d0      
      prior(:,start_ind+5) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5)   = 1.d0      
      ! Component-specific parameters
!!$      prior(:,start_ind+6) = [0.9d0, 1.1d0, 0d0, -1.d0] ! r
!!$      scale(start_ind+6)   = 1.d0      
!!$      prior(:,start_ind+7) = [0.d0, 0.3d0, 0.2d0, -1.d0] ! delta_r
!!$      scale(start_ind+7)   = 1.d0      
!!$      prior(:,start_ind+8) = [0.0d0, 0.2d0, 0.1d0, -1.d0] ! delta_z
!!$      scale(start_ind+8)   = 1.d0      
!!$      prior(:,start_ind+9) = [-60.d0, 60.d0, 0.d0, -1.d0] ! theta
!!$      scale(start_ind+9)   = 1.d0      
!!$      prior(:,start_ind+10) = [0.d0, 60.d0, 0.d0, -1.d0] ! sigma_theta
!!$      scale(start_ind+10)   = 1.d0      
    end subroutine init_interstellar_priors_and_scales

    subroutine init_fan_priors_and_scales(self, start_ind, prior, scale)
      implicit none
      class(ZodiFan),             intent(in)    :: self
      integer(i4b),               intent(in)    :: start_ind
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:),    intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0)   = 1.d-9
      prior(:,start_ind+1) = [0.d0, 00.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1)   = 1.d0      
      prior(:,start_ind+2) = [0.d0, 0.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2)   = 1.d0      
      prior(:,start_ind+3) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3)   = 1.d0      
      prior(:,start_ind+4) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4)   = 1.d0      
      prior(:,start_ind+5) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5)   = 1.d0      
      ! Component-specific parameters
!!$      prior(:,start_ind+6) = [0.9d0, 1.1d0, 0d0, -1.d0] ! r
!!$      scale(start_ind+6)   = 1.d0      
!!$      prior(:,start_ind+7) = [0.d0, 0.3d0, 0.2d0, -1.d0] ! delta_r
!!$      scale(start_ind+7)   = 1.d0      
!!$      prior(:,start_ind+8) = [0.0d0, 0.2d0, 0.1d0, -1.d0] ! delta_z
!!$      scale(start_ind+8)   = 1.d0      
!!$      prior(:,start_ind+9) = [-60.d0, 60.d0, 0.d0, -1.d0] ! theta
!!$      scale(start_ind+9)   = 1.d0      
!!$      prior(:,start_ind+10) = [0.d0, 60.d0, 0.d0, -1.d0] ! sigma_theta
!!$      scale(start_ind+10)   = 1.d0      
    end subroutine init_fan_priors_and_scales

    subroutine init_comet_priors_and_scales(self, start_ind, prior, scale)
      implicit none
      class(ZodiComet),           intent(in)    :: self
      integer(i4b),               intent(in)    :: start_ind
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:),    intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0)   = 1.d-9
      prior(:,start_ind+1) = [0.d0, 00.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1)   = 1.d0      
      prior(:,start_ind+2) = [0.d0, 0.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2)   = 1.d0      
      prior(:,start_ind+3) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3)   = 1.d0      
      prior(:,start_ind+4) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4)   = 1.d0      
      prior(:,start_ind+5) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5)   = 1.d0      
      ! Component-specific parameters
!!$      prior(:,start_ind+6) = [0.9d0, 1.1d0, 0d0, -1.d0] ! r
!!$      scale(start_ind+6)   = 1.d0      
!!$      prior(:,start_ind+7) = [0.d0, 0.3d0, 0.2d0, -1.d0] ! delta_r
!!$      scale(start_ind+7)   = 1.d0      
!!$      prior(:,start_ind+8) = [0.0d0, 0.2d0, 0.1d0, -1.d0] ! delta_z
!!$      scale(start_ind+8)   = 1.d0      
!!$      prior(:,start_ind+9) = [-60.d0, 60.d0, 0.d0, -1.d0] ! theta
!!$      scale(start_ind+9)   = 1.d0      
!!$      prior(:,start_ind+10) = [0.d0, 60.d0, 0.d0, -1.d0] ! sigma_theta
!!$      scale(start_ind+10)   = 1.d0      
    end subroutine init_comet_priors_and_scales

    subroutine param2model_cloud(self, x)
      implicit none
      class(ZodiCloud),                intent(inout) :: self
      real(dp),         dimension(1:), intent(in)    :: x
      self%n_0   = x(1)
      self%incl  = x(2)
      self%Omega = x(3)
      self%x_0   = x(4)
      self%y_0   = x(5)
      self%z_0   = x(6)
      self%alpha = x(7)
      self%beta  = x(8)
      self%gamma = x(9)
      self%mu    = x(10)      
    end subroutine param2model_cloud

    subroutine model2param_cloud(self, x)
      implicit none
      class(ZodiCloud),                intent(in)  :: self
      real(dp),         dimension(1:), intent(out) :: x
      x(1)  = self%n_0  
      x(2)  = self%incl 
      x(3)  = self%Omega 
      x(4)  = self%x_0   
      x(5)  = self%y_0   
      x(6)  = self%z_0   
      x(7)  = self%alpha 
      x(8)  = self%beta  
      x(9)  = self%gamma 
      x(10) = self%mu    
    end subroutine model2param_cloud

    subroutine param2model_band(self, x)
      implicit none
      class(ZodiBand),                 intent(inout) :: self
      real(dp),         dimension(1:), intent(in)    :: x
      self%n_0        = x(1)
      self%incl       = x(2)
      self%Omega      = x(3)
      self%x_0        = x(4)
      self%y_0        = x(5)
      self%z_0        = x(6)
      self%delta_zeta = x(7)
      self%delta_r    = x(8)
      self%v          = x(9)
      self%p          = x(10)      
    end subroutine param2model_band

    subroutine model2param_band(self, x)
      implicit none
      class(ZodiBand),                 intent(in)  :: self
      real(dp),         dimension(1:), intent(out) :: x
      x(1)  = self%n_0  
      x(2)  = self%incl 
      x(3)  = self%Omega 
      x(4)  = self%x_0   
      x(5)  = self%y_0   
      x(6)  = self%z_0   
      x(7)  = self%delta_zeta 
      x(8)  = self%delta_r
      x(9)  = self%v
      x(10) = self%p    
    end subroutine model2param_band

    subroutine param2model_ring(self, x)
      implicit none
      class(ZodiRing),                 intent(inout) :: self
      real(dp),         dimension(1:), intent(in)    :: x
      self%n_0         = x(1)
      self%incl        = x(2)
      self%Omega       = x(3)
      self%x_0         = x(4)
      self%y_0         = x(5)
      self%z_0         = x(6)
      self%R_0         = x(7)
      self%sigma_r     = x(8)
      self%sigma_z     = x(9)
      self%theta_0     = x(10)
      self%sigma_theta = x(11)      
    end subroutine param2model_ring

    subroutine model2param_ring(self, x)
      implicit none
      class(ZodiRing),                 intent(in)  :: self
      real(dp),         dimension(1:), intent(out) :: x
      x(1)  = self%n_0  
      x(2)  = self%incl 
      x(3)  = self%Omega 
      x(4)  = self%x_0   
      x(5)  = self%y_0   
      x(6)  = self%z_0   
      x(7)  = self%R_0 
      x(8)  = self%sigma_r  
      x(9)  = self%sigma_z 
      x(10) = self%theta_0
      x(11) = self%sigma_theta
    end subroutine model2param_ring

    subroutine param2model_feature(self, x)
      implicit none
      class(ZodiFeature),                intent(inout) :: self
      real(dp),           dimension(1:), intent(in)    :: x
      self%n_0         = x(1)
      self%incl        = x(2)
      self%Omega       = x(3)
      self%x_0         = x(4)
      self%y_0         = x(5)
      self%z_0         = x(6)
      self%R_0         = x(7)
      self%sigma_r     = x(8)
      self%sigma_z     = x(9)
      self%theta_0     = x(10)
      self%sigma_theta = x(11)      
    end subroutine param2model_feature

    subroutine model2param_feature(self, x)
      implicit none
      class(ZodiFeature),                intent(in)  :: self
      real(dp),           dimension(1:), intent(out) :: x
      x(1)  = self%n_0  
      x(2)  = self%incl 
      x(3)  = self%Omega 
      x(4)  = self%x_0   
      x(5)  = self%y_0   
      x(6)  = self%z_0   
      x(7)  = self%R_0 
      x(8)  = self%sigma_r  
      x(9)  = self%sigma_z 
      x(10) = self%theta_0
      x(11) = self%sigma_theta
    end subroutine model2param_feature

    subroutine param2model_interstellar(self, x)
      implicit none
      class(ZodiInterstellar),                intent(inout) :: self
      real(dp),                dimension(1:), intent(in)    :: x
      self%n_0   = x(1)
      self%incl  = x(2)
      self%Omega = x(3)
      self%x_0   = x(4)
      self%y_0   = x(5)
      self%z_0   = x(6)
    end subroutine param2model_interstellar

    subroutine model2param_interstellar(self, x)
      implicit none
      class(ZodiInterstellar),                intent(in)  :: self
      real(dp),                dimension(1:), intent(out) :: x
      x(1)  = self%n_0  
      x(2)  = self%incl 
      x(3)  = self%Omega 
      x(4)  = self%x_0   
      x(5)  = self%y_0   
      x(6)  = self%z_0   
    end subroutine model2param_interstellar

    subroutine param2model_fan(self, x)
      implicit none
      class(ZodiFan),                intent(inout) :: self
      real(dp),       dimension(1:), intent(in)    :: x
      self%n_0          = x(1)
      self%incl         = x(2)
      self%Omega        = x(3)
      self%x_0          = x(4)
      self%y_0          = x(5)
      self%z_0          = x(6)
      self%Q            = x(7)
      self%P            = x(8)
      self%gamma        = x(9)
      self%Z_midplane_0 = x(10)
      self%R_outer      = x(11)
    end subroutine param2model_fan

    subroutine model2param_fan(self, x)
      implicit none
      class(ZodiFan),                intent(in)  :: self
      real(dp),       dimension(1:), intent(out) :: x
      x(1)  = self%n_0  
      x(2)  = self%incl 
      x(3)  = self%Omega 
      x(4)  = self%x_0   
      x(5)  = self%y_0   
      x(6)  = self%z_0   
      x(7)  = self%Q 
      x(8)  = self%P  
      x(9)  = self%gamma 
      x(10) = self%Z_midplane_0
      x(11) = self%R_outer
    end subroutine model2param_fan

    subroutine param2model_comet(self, x)
      implicit none
      class(ZodiComet),                intent(inout) :: self
      real(dp),         dimension(1:), intent(in)    :: x
      self%n_0          = x(1)
      self%incl         = x(2)
      self%Omega        = x(3)
      self%x_0          = x(4)
      self%y_0          = x(5)
      self%z_0          = x(6)
      self%P            = x(7)
      self%Z_midplane_0 = x(8)
      self%R_inner      = x(9)
      self%R_outer      = x(10)      
    end subroutine param2model_comet

    subroutine model2param_comet(self, x)
      implicit none
      class(ZodiComet),                intent(in)  :: self
      real(dp),         dimension(1:), intent(out) :: x
      x(1)  = self%n_0  
      x(2)  = self%incl 
      x(3)  = self%Omega 
      x(4)  = self%x_0   
      x(5)  = self%y_0   
      x(6)  = self%z_0   
      x(7)  = self%P 
      x(8)  = self%Z_midplane_0  
      x(9)  = self%R_inner
      x(10) = self%R_outer    
    end subroutine model2param_comet

    
   subroutine get_density_cloud(self, X_vec, theta, n_out)
      class(ZodiCloud) :: self
      real(dp), dimension(:, :), intent(in) :: X_vec
      real(dp), intent(in) :: theta
      real(dp), dimension(:), intent(out) :: n_out
      integer(i4b) :: i
      real(dp) :: R, Z_midplane, zeta, g, x_prime, y_prime, z_prime

      do i = 1, size(n_out)
         x_prime = X_vec(1, i) - self%x_0
         y_prime = X_vec(2, i) - self%y_0
         z_prime = X_vec(3, i) - self%z_0

         R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
         Z_midplane = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl
         zeta = abs(Z_midplane/R)

         if (zeta < self%mu) then
            g = (zeta*zeta)/(2.d0*self%mu)
         else
            g = zeta - (0.5d0*self%mu)
         end if

         n_out(i) = self%n_0*R**(-self%alpha)*exp(-self%beta*g**self%gamma)
      end do
   end subroutine get_density_cloud

   subroutine get_density_band(self, X_vec, theta, n_out)
      class(ZodiBand) :: self
      real(dp), dimension(:, :), intent(in)  :: X_vec
      real(dp), intent(in) :: theta
      real(dp), dimension(:), intent(out) :: n_out
      integer(i4b) :: i
      real(dp) :: x_prime, y_prime, z_prime, R, Z_midplane, zeta, zeta_over_delta_zeta, term1, term2, term3, term4, R_ratio

      do i = 1, size(n_out)
         x_prime = X_vec(1, i) - self%x_0
         y_prime = X_vec(2, i) - self%y_0
         z_prime = X_vec(3, i) - self%z_0

         R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
         Z_midplane = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl
         zeta = abs(Z_midplane/R)
         zeta_over_delta_zeta = zeta/self%delta_zeta_rad
         term1 = (3.d0*self%n_0)/R
         term2 = exp(-(zeta_over_delta_zeta**6))

         ! Differs from eq 8 in K98 by a factor of 1/self.v. See Planck XIV
         ! section 4.1.2.
         ! term3 = self%v + (zeta_over_delta_zeta**self%p)
         term3 = 1.d0 + (zeta_over_delta_zeta**self%p)/self%v
         R_ratio = R/self%delta_r
         if (abs(R_ratio) > 1d12) then ! overflow
            term4 = 1.
         else 
            term4 = 1. - exp(-((R/self%delta_r)**20))
         end if

         n_out(i) = term1*term2*term3*term4
      end do
   end subroutine get_density_band

   subroutine get_density_ring(self, X_vec, theta, n_out)
      class(ZodiRing) :: self
      real(dp), dimension(:, :), intent(in)  :: X_vec
      real(dp), intent(in) :: theta
      real(dp), dimension(:), intent(out) :: n_out
      integer(i4b) :: i
      real(dp) :: x_prime, y_prime, z_prime, R, Z_midplane, term1, term2, theta_prime

      do i = 1, size(n_out)
         x_prime = X_vec(1, i) - self%x_0
         y_prime = X_vec(2, i) - self%y_0
         z_prime = X_vec(3, i) - self%z_0
         theta_prime = atan2(y_prime, x_prime) - theta - self%theta_0_rad

         ! Constraining the angle to the limit [-pi, pi]
         do while (theta_prime < -pi)
            theta_prime = theta_prime + 2.d0*pi
         end do
         do while (theta_prime > pi)
            theta_prime = theta_prime - 2.d0*pi
         end do

         R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
         Z_midplane = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl
         term1 = -((R - self%R_0)**2)/self.sigma_r**2
         term2 = abs(Z_midplane/self.sigma_z)

         if (self%sigma_theta_rad <= 0.d0 .or. abs(theta_prime) > self%sigma_theta_rad) then
            n_out(i) = self%n_0*exp(term1 - term2)
         else 
            n_out(i) = 0.d0
         end if
      end do
   end subroutine get_density_ring

   subroutine get_density_feature(self, X_vec, theta, n_out)
      class(ZodiFeature) :: self
      real(dp), dimension(:, :), intent(in) :: X_vec
      real(dp), intent(in) :: theta
      real(dp), dimension(:), intent(out) :: n_out
      integer(i4b) :: i
      real(dp) :: x_prime, y_prime, z_prime, R, Z_midplane, theta_prime, exp_term

      do i = 1, size(n_out)
         x_prime = X_vec(1, i) - self%x_0
         y_prime = X_vec(2, i) - self%y_0
         z_prime = X_vec(3, i) - self%z_0
         theta_prime = atan2(y_prime, x_prime) - theta - self%theta_0_rad

         ! Constraining the angle to the limit [-pi, pi]
         do while (theta_prime < -pi)
            theta_prime = theta_prime + 2.d0*pi
         end do
         do while (theta_prime > pi)
            theta_prime = theta_prime - 2.d0*pi
         end do

         R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
         Z_midplane = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl

         exp_term = ((R - self%R_0)**2/self%sigma_r**2) + (abs(Z_midplane)/self%sigma_z) + (theta_prime**2/self%sigma_theta_rad**2)
         n_out(i) = self%n_0*exp(-exp_term)
      end do
   end subroutine get_density_feature

   subroutine get_density_interstellar(self, X_vec, theta, n_out)
      class(ZodiInterstellar) :: self
      real(dp), dimension(:, :), intent(in) :: X_vec
      real(dp), intent(in) :: theta
      real(dp), dimension(:), intent(out) :: n_out
      integer(i4b) :: i
      real(dp) :: R, Z_midplane, zeta, g, x_prime, y_prime, z_prime

      do i = 1, size(n_out)
         n_out(i) = self%n_0
      end do
   end subroutine get_density_interstellar

   subroutine get_density_fan(self, X_vec, theta, n_out)
      class(ZodiFan) :: self
      real(dp), dimension(:, :), intent(in) :: X_vec
      real(dp), intent(in) :: theta
      real(dp), dimension(:), intent(out) :: n_out
      integer(i4b) :: i
      real(dp) :: R, Z_midplane, x_prime, y_prime, z_prime, beta, sin_beta, Z_midplane_abs, epsilon, f

      do i = 1, size(n_out)
         x_prime = X_vec(1, i) - self%x_0
         y_prime = X_vec(2, i) - self%y_0
         z_prime = X_vec(3, i) - self%z_0

         R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
         if (R > self%R_outer) then
            n_out(i) = 0.d0
            cycle
         end if

         Z_midplane = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl

         sin_beta = Z_midplane / R
         beta = asin(sin_beta)
         Z_midplane_abs = abs(Z_midplane)

         if (Z_midplane_abs < self%Z_midplane_0) then
            epsilon = 2. - (Z_midplane_abs / self%Z_midplane_0)
         else
            epsilon = 1.
         end if
         f = cos(beta) ** self%Q * exp(-self%P * sin(abs(beta) ** epsilon))
         n_out(i) = self%n_0 * R ** (-self%gamma) * f
      end do
   end subroutine get_density_fan

   subroutine get_density_comet(self, X_vec, theta, n_out)
      class(ZodiComet) :: self
      real(dp), dimension(:, :), intent(in) :: X_vec
      real(dp), intent(in) :: theta
      real(dp), dimension(:), intent(out) :: n_out
      integer(i4b) :: i
      real(dp) :: R, Z_midplane, x_prime, y_prime, z_prime, beta, sin_beta, Z_midplane_abs, epsilon, f

      do i = 1, size(n_out)
         x_prime = X_vec(1, i) - self%x_0
         y_prime = X_vec(2, i) - self%y_0
         z_prime = X_vec(3, i) - self%z_0

         R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
         if ((R > self%R_outer) .or. (R < self%R_inner)) then
            n_out(i) = 0.d0
            cycle
         end if
         Z_midplane = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl

         sin_beta = Z_midplane / R
         beta = asin(sin_beta)
         Z_midplane_abs = abs(Z_midplane)

         if (Z_midplane_abs < self%Z_midplane_0) then
            epsilon = 2. - (Z_midplane_abs / self%Z_midplane_0)
         else
            epsilon = 1.
         end if
         f = exp(-self%P * sin(abs(beta) ** epsilon))
         n_out(i) = 0.37 * self%n_0 * f / R
      end do
   end subroutine get_density_comet

   subroutine init_comps(self, params, comp_types, param_labels)
      ! Initializes the components in the zodi model and computes the number of parameters in the model.
      class(ZodiModel), target, intent(inout) :: self
      real(dp), intent(in) :: params(:, :)
      character(len=*), intent(in) :: comp_types(:)
      class(InterplanetaryDustParamLabels), intent(in) :: param_labels
      integer(i4b) :: i, ierr
      allocate (self%comps(self%n_comps))
      self%n_params = self%n_general_params
      ! NOTE: The order of the parameters in the `params` array must match the below order of readin
      do i = 1, self%n_comps
         !self%comps(i)%labels = [param_labels%common]
         select case (trim(adjustl(comp_types(i))))
         case ('cloud')
            allocate (ZodiCloud::self%comps(i)%c)
            self%comps(i)%c = ZodiCloud(&
                & n_0=params(i, 1), &
                & incl=params(i, 2), &
                & Omega=params(i, 3), &
                & x_0=params(i, 4), &
                & y_0=params(i, 5), &
                & z_0=params(i, 6), &
                & alpha=params(i, 7), &
                & beta=params(i, 8), &
                & gamma=params(i, 9), &
                & mu=params(i, 10) &
                &)
            allocate(self%comps(i)%labels(10))
            self%comps(i)%labels = [param_labels%common, param_labels%cloud]
         case ('band')
            allocate (ZodiBand::self%comps(i)%c)
            self%comps(i)%c = ZodiBand(&
                & n_0=params(i, 1), &
                & incl=params(i, 2), &
                & Omega=params(i, 3), &
                & x_0=params(i, 4), &
                & y_0=params(i, 5), &
                & z_0=params(i, 6), &
                & delta_zeta=params(i, 7), &
                & delta_r=params(i, 8), &
                & v=params(i, 9), &
                & p=params(i, 10) &
            &)
            allocate(self%comps(i)%labels(10))
            self%comps(i)%labels = [param_labels%common, param_labels%band]
         case ('ring')
            allocate (ZodiRing::self%comps(i)%c)
            self%comps(i)%c = ZodiRing(&
               & n_0=params(i, 1), &
               & incl=params(i, 2), &
               & Omega=params(i, 3), &
               & x_0=params(i, 4), &
               & y_0=params(i, 5), &
               & z_0=params(i, 6), &
               & R_0=params(i, 7), &
               & sigma_r=params(i, 8), &
               & sigma_z=params(i, 9), &
               & theta_0=params(i, 10), &
               & sigma_theta=params(i, 11) &
            &)
            allocate(self%comps(i)%labels(11))
            self%comps(i)%labels = [param_labels%common, param_labels%ring]
         case ('feature')
            allocate (ZodiFeature::self%comps(i)%c)
            self%comps(i)%c = ZodiFeature(&
               & n_0=params(i, 1), &
               & incl=params(i, 2), &
               & Omega=params(i, 3), &
               & x_0=params(i, 4), &
               & y_0=params(i, 5), &
               & z_0=params(i, 6), &
               & R_0=params(i, 7), &
               & sigma_r=params(i, 8), &
               & sigma_z=params(i, 9), &
               & theta_0=params(i, 10), &
               & sigma_theta=params(i, 11) &
               &)
            allocate(self%comps(i)%labels(11))
            self%comps(i)%labels = [param_labels%common, param_labels%feature]
         case ('interstellar')
            allocate (ZodiInterstellar::self%comps(i)%c)
            self%comps(i)%c = ZodiInterstellar(&
                 & n_0=params(i, 1), &
                 & incl=params(i, 2), &
                 & Omega=params(i, 3), &
                 & x_0=params(i, 4), &
                 & y_0=params(i, 5), &
                 & z_0=params(i, 6), &
                 & R=params(i, 7), &
                 & alpha=params(i, 8) &
                 &)
            allocate(self%comps(i)%labels(8))
            self%comps(i)%labels = [param_labels%common, param_labels%interstellar]
         case ('fan')
            allocate (ZodiFan::self%comps(i)%c)
            self%comps(i)%c = ZodiFan(&
                 & n_0=params(i, 1), &
                 & incl=params(i, 2), &
                 & Omega=params(i, 3), &
                 & x_0=params(i, 4), &
                 & y_0=params(i, 5), &
                 & z_0=params(i, 6), &
                 & Q=params(i, 7), &
                 & P=params(i, 8), &
                 & gamma=params(i, 9), &
                 & Z_midplane_0=params(i, 10), &
                 & R_outer=params(i, 11) &
            &)
            allocate(self%comps(i)%labels(11))
            self%comps(i)%labels = [param_labels%common, param_labels%fan]
         case ('comet')
            allocate (ZodiFan::self%comps(i)%c)
            self%comps(i)%c = ZodiComet(&
                 & n_0=params(i, 1), &
                 & incl=params(i, 2), &
                 & Omega=params(i, 3), &
                 & x_0=params(i, 4), &
                 & y_0=params(i, 5), &
                 & z_0=params(i, 6), &
                 & P=params(i, 7), &
                 & Z_midplane_0=params(i, 8), &
                 & R_inner=params(i, 9), &
                 & R_outer=params(i, 10) &
                 &)
            allocate(self%comps(i)%labels(10))
            self%comps(i)%labels = [param_labels%common, param_labels%comet]
         case default
            print *, 'Invalid zodi component type in zodi `init_from_params`:', trim(adjustl(comp_types(i)))
            stop
         end select
         self%comps(i)%npar = size(self%comps(i)%labels)
         call self%comps(i)%c%init()
         self%n_params = self%n_params + size(self%comps(i)%labels)
      end do
   end subroutine



    subroutine model_to_params2(self, x, labels)
      ! Dumps a zodi model to a parameter vector `x`. If `labels` is present, it is populated with
      ! the corresponding parameter labels.
      class(ZodiModel), intent(in) :: self
      real(dp), intent(out) :: x(:)
      character(len=*), allocatable, optional, intent(inout) :: labels(:)
      character(len=128), allocatable :: labels_copy(:), comp_label_upper(:)
      integer(i4b) :: i, j, running_idx

      if (size(x) /= self%n_params) stop "Error: argument 'x' has the wrong size. must be `size(zodi_model%n_params)`"
      if (present(labels)) then
         if (allocated(labels)) stop "`labels` must not be allocated at the time of passing it in to `model_to_params`"
         allocate(comp_label_upper(self%n_comps))
      end if

      running_idx = 0
      do i = 1, self%n_comps
         x(running_idx + 1) = self%comps(i)%c%n_0
         x(running_idx + 2) = self%comps(i)%c%incl
         x(running_idx + 3) = self%comps(i)%c%Omega
         x(running_idx + 4) = self%comps(i)%c%x_0
         x(running_idx + 5) = self%comps(i)%c%y_0
         x(running_idx + 6) = self%comps(i)%c%z_0
         running_idx = running_idx + self%n_common_params
         select type (comp => self%comps(i)%c)
         class is (ZodiCloud)
            x(running_idx + 1) = comp%alpha
            x(running_idx + 2) = comp%beta
            x(running_idx + 3) = comp%gamma
            x(running_idx + 4) = comp%mu
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiBand)
            x(running_idx + 1) = comp%delta_zeta
            x(running_idx + 2) = comp%delta_r
            x(running_idx + 3) = comp%v
            x(running_idx + 4) = comp%p
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiRing)
            x(running_idx + 1) = comp%R_0
            x(running_idx + 2) = comp%sigma_r
            x(running_idx + 3) = comp%sigma_z
            x(running_idx + 4) = comp%theta_0
            x(running_idx + 5) = comp%sigma_theta
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiFeature)
            x(running_idx + 1) = comp%R_0
            x(running_idx + 2) = comp%sigma_r
            x(running_idx + 3) = comp%sigma_z
            x(running_idx + 4) = comp%theta_0
            x(running_idx + 5) = comp%sigma_theta
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiInterstellar)
            x(running_idx + 1) = comp%R
            x(running_idx + 2) = comp%alpha
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiFan)
            x(running_idx + 1) = comp%Q
            x(running_idx + 2) = comp%P
            x(running_idx + 3) = comp%gamma
            x(running_idx + 4) = comp%Z_midplane_0
            x(running_idx + 5) = comp%R_outer
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
            class is (ZodiComet)
            x(running_idx + 1) = comp%P
            x(running_idx + 2) = comp%Z_midplane_0
            x(running_idx + 3) = comp%R_inner
            x(running_idx + 4) = comp%R_outer
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         end select
         if (present(labels)) then
            labels_copy = self%comps(i)%labels
            comp_label_upper(i) = self%comp_labels(i)
            call toupper(comp_label_upper(i))
            do j = 1, size(labels_copy)
               labels_copy(j) = trim(adjustl(comp_label_upper(i)))//'_'//trim(adjustl(labels_copy(j))) 
            end do
               labels = [labels, labels_copy]
         end if
      end do
      x(running_idx + 1) = self%T_0
      x(running_idx + 2) = self%delta
      if (present(labels)) then
         labels = [labels, self%general_labels]
      end if
    end subroutine model_to_params2

   subroutine params_to_model2(self, x)
      ! Dumps a zodi model to a parameter vector `x`.
      class(ZodiModel), intent(inout) :: self
      real(dp), intent(in) :: x(:)
      integer(i4b) :: i, running_idx

      if (size(x) /= self%n_params) then
         write(*,*) "Error: argument 'x' has the wrong size. must be `size(zodi_model%n_params)`", size(x), self%n_params
         stop
      end if

      running_idx = 0
      do i = 1, self%n_comps
         self%comps(i)%c%n_0 = x(running_idx + 1)
         self%comps(i)%c%incl = mod(x(running_idx + 2), 360.) ! degree prior
         self%comps(i)%c%Omega = mod(x(running_idx + 3), 360.) ! degree prior
         self%comps(i)%c%x_0 = x(running_idx + 4)
         self%comps(i)%c%y_0 = x(running_idx + 5)
         self%comps(i)%c%z_0 = x(running_idx + 6)
         running_idx = running_idx + self%n_common_params
         
         ! The order of these operations much match the order tabulated in the labels in `InterplanetaryDustParamLabels`
         select type (comp => self%comps(i)%c)
         class is (ZodiCloud)
            comp%alpha = x(running_idx + 1)
            comp%beta = x(running_idx + 2)
            comp%gamma = x(running_idx + 3)
            comp%mu = x(running_idx + 4)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiBand)
            comp%delta_zeta = mod(x(running_idx + 1), 360.) ! degree prior
            comp%delta_r = x(running_idx + 2)
            comp%v = x(running_idx + 3)
            comp%p = x(running_idx + 4)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiRing)
            comp%R_0 = x(running_idx + 1)
            comp%sigma_r = x(running_idx + 2)
            comp%sigma_z = x(running_idx + 3)
            comp%theta_0 = mod(x(running_idx + 4), 360.) ! degree prior
            comp%sigma_theta = mod(x(running_idx + 5), 360.) ! degree prior
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiFeature)
            comp%R_0 = x(running_idx + 1)
            comp%sigma_r = x(running_idx + 2)
            comp%sigma_z = x(running_idx + 3)
            comp%theta_0 = mod(x(running_idx + 4), 360.) ! degree prior
            comp%sigma_theta = mod(x(running_idx + 5), 360.) ! degree prior
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiInterstellar)
            comp%R = x(running_idx + 1)
            comp%alpha = x(running_idx + 2)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiFan)
            comp%Q = x(running_idx + 1)
            comp%P = x(running_idx + 2)
            comp%gamma = x(running_idx + 3)
            comp%Z_midplane_0 = x(running_idx + 4)
            comp%R_outer = x(running_idx + 5)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiComet)
            comp%P = x(running_idx + 1)
            comp%Z_midplane_0 = x(running_idx + 2)
            comp%R_inner = x(running_idx + 3)
            comp%R_outer = x(running_idx + 4)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         end select
         call self%comps(i)%c%init()
      end do
      self%T_0 = x(running_idx + 1)
      self%delta = x(running_idx + 2)
    end subroutine params_to_model2

   subroutine model_to_chain(self, cpar, iter)
      ! Dumps the zodi model to the chain file
      class(ZodiModel), intent(in) :: self
      type(comm_params), intent(in) :: cpar
      integer(i4b), intent(in) :: iter

      integer(i4b) :: i, j, hdferr, ierr, unit, n_comp_params, param_idx
      logical(lgt) :: exist, init, new_header
      character(len=6) :: itext
      character(len=4) :: ctext
      character(len=512) :: zodi_path, comp_path, comp_group_path, param_path, chainfile, hdfpath, param_label, general_group_path
      character(len=32), allocatable :: labels(:)
      real(dp), allocatable :: params(:)
      type(hdf_file) :: file

      if (.not. cpar%myid_chain == 0) return

      call int2string(cpar%mychain, ctext)
      chainfile = trim(adjustl(cpar%outdir))//'/chain'// &
          & '_c'//trim(adjustl(ctext))//'.h5'

      inquire (file=trim(chainfile), exist=exist)
      call open_hdf_file(chainfile, file, 'b')

      call int2string(iter, itext)
      zodi_path = trim(adjustl(itext))//'/zodi'
      call create_hdf_group(file, trim(adjustl(zodi_path)))


      general_group_path = trim(adjustl(zodi_path))//'/general'
      call create_hdf_group(file, trim(adjustl(general_group_path)))

      comp_group_path = trim(adjustl(zodi_path))//'/comps'
      call create_hdf_group(file, trim(adjustl(comp_group_path)))

      allocate(params(self%n_params))
      call self%model_to_params2(params, labels)
      
      param_idx = 0 
      do i = 1, self%n_comps
         comp_path = trim(adjustl(comp_group_path))//'/'//trim(adjustl(self%comp_labels(i)))//'/'
         call create_hdf_group(file, trim(adjustl(comp_path)))
         n_comp_params = size(self%comps(i)%labels)
         do j = 1, n_comp_params
            param_label = trim(adjustl(comp_path))//'/'//trim(adjustl(self%comps(i)%labels(j)))
            call write_hdf(file, trim(adjustl(param_label)), params(param_idx + j))
         end do
         param_idx = param_idx + n_comp_params
      end do
      if (param_idx + self%n_general_params /= self%n_params) stop "Error: param_idx + self%n_general_params /= self%n_params"
      do i = 1, self%n_general_params
         param_label = trim(adjustl(general_group_path))//'/'//trim(adjustl(self%general_labels(i)))
         call write_hdf(file, trim(adjustl(param_label)), params(param_idx + i))
      end do
      call close_hdf_file(file)
      deallocate(params,labels)
   end subroutine

   subroutine comp_from_chain(self, cpar, params, comp_idx)
      ! Initialize a component from a chain
      class(ZodiModel), target, intent(inout) :: self
      type(comm_params), intent(in) :: cpar
      real(dp), intent(inout) :: params(:, :)
      integer(i4b), intent(in) :: comp_idx

      logical(lgt) :: exist
      integer(i4b) :: i, j, l, ierr, initsamp
      character(len=6) :: itext

      type(hdf_file) :: file

      character(len=32), allocatable :: param_labels(:)
      character(len=512) :: chainfile, group_name

      if (cpar%myid == cpar%root) then
         if (trim(cpar%zs_init_hdf(comp_idx)) == 'default') then
            call get_chainfile_and_samp(trim(cpar%init_chain_prefixes(1)), chainfile, initsamp)
         else
            call get_chainfile_and_samp(trim(cpar%zs_init_hdf(comp_idx)), chainfile, initsamp)
         end if
         inquire (file=trim(chainfile), exist=exist)
         if (.not. exist) call report_error('Zodi init chain does not exist = '//trim(chainfile))
         l = len(trim(chainfile))
         if (.not. ((trim(chainfile(l-2:l)) == '.h5') .or. (trim(chainfile(l-3:l)) == '.hd5'))) call report_error('Zodi init chain must be a .h5 file')
         
         call open_hdf_file(trim(chainfile), file, "r")
         
         call int2string(initsamp, itext)
         group_name = trim(adjustl(itext)//'/zodi/comps/'//trim(adjustl(cpar%zs_comp_labels(comp_idx))))
         if (.not. hdf_group_exists(file, group_name)) then
            print *, "zodi component: ", trim(adjustl(cpar%zs_comp_labels(comp_idx))), "not found in chain:", trim(chainfile)
            stop
         end if 

         param_labels = cpar%zodi_param_labels%get_labels(trim(adjustl(cpar%zs_comp_types(comp_idx))), add_common=.true.)
         do j = 1, size(param_labels)
            call read_hdf(file, trim(adjustl(itext)//'/zodi/comps/'//trim(adjustl(cpar%zs_comp_labels(comp_idx)))// &
                  & '/'//trim(adjustl(param_labels(j)))), params(comp_idx, j))
         end do
         deallocate(param_labels)
      end if
      ! call mpi_bcast(params, sum(shape(params)), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
      call mpi_bcast(params, size(params, dim=1) * size(params, dim=2), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)

   end subroutine
   ! subroutine model_from_chain(self, cpar)
   !    ! Initializes parts of the zodi model from the chain file
   !    class(ZodiModel), target, intent(inout) :: self
   !    type(comm_params), intent(in) :: cpar
   !    logical(lgt) :: exist
   !    integer(i4b) :: i, j, l, unit, ierr, initsamp, hdferr
   !    character(len=6) :: itext

   !    type(hdf_file) :: file
   !    real(dp) :: comp_params(100), params(100, 100), general_params(100)
   !    character(len=32), allocatable :: common_param_labels(:), param_labels(:)
   !    character(len=512) :: chainfile, group_name
   !    TYPE(h5o_info_t) :: object_info

   !    params = 0.
   !    if (cpar%myid == cpar%root) then
   !       call get_chainfile_and_samp(trim(cpar%zs_init_chain), chainfile, initsamp)
   !       inquire (file=trim(chainfile), exist=exist)
   !       if (.not. exist) call report_error('Zodi init chain does not exist = '//trim(chainfile))
   !       l = len(trim(chainfile))
   !       if (.not. ((trim(chainfile(l-2:l)) == '.h5') .or. (trim(chainfile(l-3:l)) == '.hd5'))) call report_error('Zodi init chain must be a .h5 file')
         
   !       call open_hdf_file(trim(chainfile), file, "r")
         
   !       call int2string(initsamp, itext)
   !       do i = 1, cpar%zs_ncomps
   !          group_name = trim(adjustl(itext)//'/zodi/comps/'//trim(adjustl(cpar%zs_comp_labels(i))))
   !          if (.not. hdf_group_exists(file, group_name)) then
   !             params(i, :) = cpar%zs_comp_params(i, :, 1)
   !             cycle
   !          end if 
   !          param_labels = cpar%zodi_param_labels%get_labels(trim(adjustl(cpar%zs_comp_types(i))), add_common=.true.)
   !          do j = 1, size(param_labels)
   !             call read_hdf(file, trim(adjustl(itext)//'/zodi/comps/'//trim(adjustl(cpar%zs_comp_labels(i)))// &
   !                 & '/'//trim(adjustl(param_labels(j)))), comp_params(j))
   !          end do
   !             params(i, :) = comp_params
   !       end do 
   !       do i = 1, self%n_general_params
   !          call read_hdf(file, trim(adjustl(itext)//'/zodi/general/'//trim(adjustl(self%general_labels(i)))), general_params(i))
   !       end do
   !       call close_hdf_file(file)
   !    end if
   !    call mpi_bcast(params, size(params, dim=1) * size(params, dim=2), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
   !    call mpi_bcast(general_params, size(general_params), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
   !    call self%init_general_params(general_params)
   !    call self%init_comps(params, cpar%zs_comp_types, cpar%zodi_param_labels)

   ! end subroutine model_from_chain

   subroutine get_s_zodi(s_therm, s_scat, s_zodi, emissivity, albedo)
      ! Evaluates the zodiacal signal (eq. 20 in ZodiPy paper [k98 model]) given
      ! integrated thermal zodiacal emission and scattered zodiacal light.
      !
      ! Parameters:
      ! -----------
      ! s_scat :
      !     Integrated contribution from scattered sunlight light.
      ! s_therm :
      !     Integrated contribution from thermal interplanetary dust emission.
      ! s_zodi :
      !     Zodiacal signal.
      ! emissivity :
      !     Emissivity of the zodiacal components.
      ! albedo :
      !     Albedo of the zodiacal components.
      real(sp), dimension(:, :), intent(in) :: s_scat, s_therm
      real(sp), dimension(:), intent(out)   :: s_zodi
      real(dp), dimension(:), intent(in) :: emissivity, albedo

      integer(i4b) :: i

      s_zodi = 0.
      do i = 1, size(emissivity)
         s_zodi = s_zodi + ((s_scat(:,i) * albedo(i)) + (1. - albedo(i)) * emissivity(i) * s_therm(:,i))
      end do
   end subroutine get_s_zodi

   function get_par_ind(self, comp, comp_str, param, em_band, al_band, em_string, al_string, mono_band, mono_string)
     implicit none
     class(ZodiModel),     intent(in)           :: self
     class(ZodiComponentContainer), intent(in), target, optional :: comp
     character(len=*),     intent(in), optional :: comp_str
     character(len=*),     intent(in), optional :: param
     integer(i4b),         intent(in), optional :: em_band, al_band, mono_band
     character(len=*),     intent(in), optional :: em_string, al_string, mono_string
     integer(i4b)                               :: get_par_ind
     
     integer(i4b) :: i, band
     class(ZodiComponentContainer), pointer :: c
     character(len=128) :: str1, str2
     
     if (present(comp) .or. present(comp_str)) then
        if (present(comp_str)) then
           i =  get_string_index(self%comp_labels, comp_str)
           c => zodi_model%comps(i)
        else
           c => comp
        end if

        if (present(param)) then
           get_par_ind = c%start_ind + get_string_index(c%labels, param)-1
        else if (present(em_band)) then
           get_par_ind = c%start_ind + c%npar + em_band - 1
        else if (present(al_band)) then
           get_par_ind = c%start_ind + c%npar + numband + al_band - 1
        else if (present(em_string)) then
           get_par_ind = c%start_ind + c%npar + get_string_index(band_labels, em_string) - 1
        else if (present(al_string)) then
           get_par_ind = c%start_ind + c%npar + numband + get_string_index(band_labels, al_string) - 1
        else
           write(*,*) 'get_par_ind error: Need parameter specification'
           stop
        end if
     else
        if (present(param)) then
           ! General parameter
           get_par_ind = get_string_index(self%general_labels, param)
        else if (present(mono_string) .or. present(mono_band)) then
           ! Monopoles are always last
           if (present(mono_string)) then
              band = get_string_index(band_labels, mono_string)
           else
              band = mono_band
           end if
           get_par_ind = self%npar_tot - numband + band
        else
           write(*,*) 'get_par_ind error: Need parameter specification'
           stop
        end if
     end if
     
   end function get_par_ind

   subroutine samp_group2stat(cpar, samp_group, active, stat)
     implicit none     
     type(comm_params), intent(in) :: cpar
     integer(i4b),                   intent(in)    :: samp_group
     logical(lgt),     dimension(:), intent(in)    :: active
     integer(i4b),     dimension(:), intent(inout) :: stat

     integer(i4b) :: i, c, j, k, first, last, n_params, n, m, ind, em_global, al_global, c_to, c_from, band
     character(len=128) :: tokens(100), comp_param(2), wire_from(2), wire_to(2), label, param_label_tokens(10), em_from(2), em_to(2)
     character(len=2048) :: str
     
     ! Default: Fix everything at input
     str  = cpar%zs_samp_groups(samp_group)
     stat = -1
     
     ! Parse user directives
     call get_tokens(str, ',', tokens, n_params) 

     write(*,*) 'a'
     em_global = 0; al_global = 0
     if (cpar%zs_em_global /= 'none') em_global = get_string_index(zodi_model%comp_labels, cpar%zs_em_global)
     if (cpar%zs_al_global /= 'none') al_global = get_string_index(zodi_model%comp_labels, cpar%zs_al_global)
     do i = 1, n_params
        call get_tokens(tokens(i), ':', comp_param, num=n)
        write(*,*) 'a1', tokens(i)
        if (n == 1) then
           ! General parameter
           ind = zodi_model%get_par_ind(param=comp_param(1))
           stat(ind) = 0
        else if (n == 2) then
           if (trim(comp_param(1)) == 'em') then
              write(*,*) 'a2', tokens(i)
              band = get_string_index(band_labels, comp_param(2))
              if (.not. active(band)) cycle
              do j = 1, zodi_model%n_comps
                 ind       = zodi_model%get_par_ind(comp=zodi_model%comps(j), em_string=comp_param(2))
                 stat(ind) = 0
              end do
              cycle
           else if (trim(comp_param(1)) == 'al') then
              write(*,*) 'a3', tokens(i)
              band = get_string_index(band_labels, comp_param(2))
              if (.not. active(band)) cycle
              do j = 1, zodi_model%n_comps
                 ind       = zodi_model%get_par_ind(comp=zodi_model%comps(j), al_string=comp_param(2))
                 stat(ind) = 0
              end do
              cycle
           end if

           ! Component parameters
           label = comp_param(2)
           c     = get_string_index(zodi_model%comp_labels, comp_param(1))
           first = zodi_model%comps(c)%start_ind
           if (trim(label) == 'all') then
              write(*,*) 'a4', label
              last  = first + zodi_model%comps(c)%npar - 1
              stat(first:last) = 0 ! Activate all
              do j = 1, numband
                 if (.not. active(j)) cycle
                 stat(last+j)         = 0
                 stat(last+numband+j) = 0
              end do
              last = last + 2*numband
           else if (trim(label(1:2)) == 'em') then
              write(*,*) 'a5', label
              ! Emissivity
              call get_tokens(label, '@', comp_param, num=n)
              if (n == 1) then
                 first = first + zodi_model%comps(c)%npar
                 do j = 1, numband
                    if (active(j)) stat(first+j-1) = 0 
                 end do
              else
                 band = get_string_index(band_labels, comp_param(2))
                 if (.not. active(band)) cycle
                 ind = zodi_model%get_par_ind(comp=zodi_model%comps(c), em_string=comp_param(2))
                 stat(ind) = 0
              end if
           else if (trim(label(1:2)) == 'al') then
              write(*,*) 'a6', label
              ! Albedo
              call get_tokens(label, '@', comp_param, num=n)
              if (n == 1) then
                 first = first + zodi_model%comps(c)%npar + numband
                 do j = 1, numband
                    if (active(j)) stat(first+j-1) = 0 
                 end do
              else
                 band = get_string_index(band_labels, comp_param(2))
                 if (.not. active(band)) cycle
                 ind = zodi_model%get_par_ind(comp=zodi_model%comps(c), al_string=comp_param(2))
                 stat(ind) = 0
              end if
           else
              write(*,*) 'a7', label
              ! Shape parameter
              ind = zodi_model%get_par_ind(comp=zodi_model%comps(c), param=comp_param(2))
              stat(ind) = 0
           end if
        else
           write(*,*) 'Invalid zodi samp group element = ', trim(tokens(i))
           stop
        end if
     end do

          write(*,*) 'b'
     ! Set up monopoles
     do i = 1, numband
        ind = zodi_model%get_par_ind(mono_band=i)
        if (active(i) .and. cpar%zs_joint_mono .and. band_todtype(i) /= 'none') stat(ind) = 0
     end do

          write(*,*) 'c'
     ! Apply explicit parameter wiring
     call get_tokens(cpar%zs_wiring, ',', tokens, n_params) 
     do i = 1, n_params
        if (trim(tokens(1)) == 'none') exit
        
        call get_tokens(tokens(i), '>', comp_param, num=n)
        if (n /= 2) then
           write(*,*) 'Error: Zodi wiring must contain two tokens: ', tokens(i)
           stop
        end if

        call get_tokens(comp_param(1), ':', wire_from, num=n)
        call get_tokens(comp_param(2), ':', wire_to,   num=m)
        if (n /= 2) then
           write(*,*) 'Error: Zodi wiring must contain two tokens: ', trim(tokens(i))!, trim(comp_param(1)), trim(wire_from(1)), trim(wire_from(2))
           stop
        end if
        if (m /= 2) then
           write(*,*) 'Error: Zodi wiring must contain two tokens: ', trim(tokens(i))!, trim(comp_param(2)), trim(wire_to(1)), trim(wire_to(2))
           write(*,*) 'Error: Zodi wiring must contain two tokens: ', trim(comp_param(2))!, trim(wire_to(1)), trim(wire_to(2))
           stop
        end if

        if (trim(wire_from(2)(1:2)) == 'em') then
           call get_tokens(wire_from(2), '@', em_from, num=n)
           call get_tokens(wire_to(2), '@', em_to, num=n)
           if (n == 1) then
              ! Attach all bands
              do c = 1, numband
                 if (.not. active(c)) cycle
                 c_from = zodi_model%get_par_ind(comp_str=wire_from(1), em_band=c)
                 c_to   = zodi_model%get_par_ind(comp_str=wire_to(1),   em_band=c)
                 stat(c_from) = c_to
              end do
           else
              ! Attach specified band
              j = get_string_index(band_labels, em_from(2))
              k = get_string_index(band_labels, em_to(2))
              if (.not. active(j) .or. .not. active(k)) cycle
              c_from = zodi_model%get_par_ind(comp_str=wire_from(1), em_string=em_from(2))
              c_to   = zodi_model%get_par_ind(comp_str=wire_to(1),   em_string=em_to(2))
              stat(c_from) = c_to
           end if
        else
           write(*,*) 'Unsupported zodi wire = ', trim(tokens(i))
           stop
        end if
     end do

          write(*,*) 'd'
     ! Apply global directives
     if (em_global > 0) then
        do i = 1, zodi_model%n_comps
           if (i == em_global) cycle
           do j = 1, numband
              if (.not. active(j)) cycle
              ind   = zodi_model%get_par_ind(comp=zodi_model%comps(em_global), em_band=j)
              first = zodi_model%get_par_ind(comp=zodi_model%comps(i), em_band=j)
              stat(first) = ind
           end do
        end do
     end if

     if (al_global > 0) then
        do i = 1, zodi_model%n_comps
           if (i == al_global) cycle
           do j = 1, numband
              if (.not. active(j)) cycle
              ind   = zodi_model%get_par_ind(comp=zodi_model%comps(al_global), al_band=j)
              first = zodi_model%get_par_ind(comp=zodi_model%comps(i), al_band=j)
              stat(first) = ind
           end do
        end do
     end if

          write(*,*) 'e'
     ! Match emissivity and albedo for bands with identical instruments
     do i = 1, numband
        if (.not. active(i)) cycle
        do j = i+1, numband
           if (.not. active(j)) cycle
           if (trim(band_instlabels(j)) == trim(band_instlabels(i))) then
              do c = 1, zodi_model%n_comps
                 ind   = zodi_model%get_par_ind(comp=zodi_model%comps(c), em_band=i)
                 first = zodi_model%get_par_ind(comp=zodi_model%comps(c), em_band=j)
                 stat(first) = ind
                 
                 ind   = zodi_model%get_par_ind(comp=zodi_model%comps(c), al_band=i)
                 first = zodi_model%get_par_ind(comp=zodi_model%comps(c), al_band=j)
                 stat(first) = ind
              end do
           end if
        end do
     end do

          write(*,*) 'f'
     ! Apply absolute constraints
     do j = 1, numband
        if (band_todtype(j) == 'none') then
           do i = 1, zodi_model%n_comps
              ind = zodi_model%get_par_ind(comp=zodi_model%comps(i), em_band=j)
              stat(ind)         = -2 ! Fix emissivity to zero
              stat(ind+numband) = -2 ! Fix albedo to zero
           end do
        end if

        if (band_todtype(j) /= 'none' .and. band_nu_c(j) > 1d14) then
           do i = 1, zodi_model%n_comps
              ind = zodi_model%get_par_ind(comp=zodi_model%comps(i), em_band=j)
              stat(ind)         = -3 ! Fix emissivity to unity
           end do
        end if

        if (band_todtype(j) /= 'none' .and. band_nu_c(j) < 70000d9) then
           do i = 1, zodi_model%n_comps
              ind = zodi_model%get_par_ind(comp=zodi_model%comps(i), al_band=j)
              stat(ind)         = -2 ! Fix albedo to zero
           end do
        end if

        if (trim(zodi_refband) == trim(band_labels(j))) then
           do i = 1, zodi_model%n_comps
              ind = zodi_model%get_par_ind(comp=zodi_model%comps(i), em_band=j)
              stat(ind)         = -3 ! Fix emissivity to unity
           end do
        end if
     end do

          write(*,*) 'g'
     ! Short-cut multi-leg wires
     do i = 1, size(stat)
        if (stat(i) > 0) then
           j = stat(i)
           do while (stat(j) > 0)
              stat(i) = stat(j)
              j       = stat(i)
              if (j <= 0) exit
           end do
        end if
     end do

          write(*,*) 'h'
     
   end subroutine samp_group2stat

   function get_string_index(arr, str)
     implicit none
     character(len=*), dimension(:), intent(in) :: arr
     character(len=*),               intent(in) :: str
     integer(i4b)                               :: get_string_index

     integer(i4b) :: i
     character(len=128) :: str1, str2

     str1 = str
     call toupper(str1)
     do i = 1, size(arr)
        str2 = arr(i)
        call toupper(str2)
        if (trim(str1) == trim(str2)) then
           get_string_index = i
           exit
        end if
     end do
     if (i > size(arr)) then
        write(*,*) 'get_string_index: String not found = ', trim(str)
        stop
     end if

   end function get_string_index

 end module comm_zodi_mod
