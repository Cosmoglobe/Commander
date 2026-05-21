module comm_zodi_comp_mod
   use comm_param_mod
   use comm_hdf_mod
   implicit none

   type, abstract :: ZodiComponent
      real(dp) :: x_0, y_0, z_0, incl, Omega, n_0, sed_ampl, sed_cutoff, sed_b
      real(dp), allocatable :: sin_omega, cos_omega, sin_incl, cos_incl
      real(dp), allocatable :: emissivity(:), albedo(:)
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
        real(dp), dimension(1:,1:), intent(inout) :: scale
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

   type, extends(ZodiComponent) :: ZodiWrightCloudRing
      real(dp) :: p1, p3, p4, p5, p6, p7, p8, p9, p10, p13, p14, p15, p16, p17, p18, p19
   contains
      procedure :: init => init_WrightCloudRing
      procedure :: get_density => get_density_WrightCloudRing
      procedure :: init_priors_and_scales => init_WrightCloudRing_priors_and_scales
      procedure :: param2model => param2model_WrightCloudRing
      procedure :: model2param => model2param_WrightCloudRing
   end type ZodiWrightCloudRing

   type, extends(ZodiComponent) :: ZodiWrightBand
      real(dp) :: q1, q5, q6, q7, q8, R_1
   contains
      procedure :: init => init_WrightBand
      procedure :: get_density => get_density_WrightBand
      procedure :: init_priors_and_scales => init_WrightBand_priors_and_scales
      procedure :: param2model => param2model_WrightBand
      procedure :: model2param => model2param_WrightBand
   end type ZodiWrightBand

contains

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
   
   subroutine init_WrightCloudRing(self)
     class(ZodiWrightCloudRing) :: self
     call self%init_base_comp()
   end subroutine init_WrightCloudRing

   subroutine init_WrightBand(self)
     class(ZodiWrightBand) :: self
     call self%init_base_comp()
   end subroutine init_WrightBand

    
   subroutine init_cloud_priors_and_scales(self, start_ind, prior, scale)
     implicit none
     class(ZodiCloud),           intent(in)    :: self
     integer(i4b),               intent(in)    :: start_ind
     real(dp), dimension(1:,1:), intent(inout) :: prior
     real(dp), dimension(1:,1:), intent(inout) :: scale
      
      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0,:) = [1.d-9, 4.d-9]
      prior(:,start_ind+1) = [-30.d0, 30.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1,:) = [1.d0, 0.03d0]
      prior(:,start_ind+2) = [-720.d0, 720.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2,:) = [1.d0, 0.3d0]
      !prior(:,start_ind+3) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! X_0
      prior(:,start_ind+3) = [-0.04d0, 0.04d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3,:) = [1.d0, 1d-3]
      prior(:,start_ind+4) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4,:) = [1.d0, 0.8d-3]
      prior(:,start_ind+5) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5,:) = [1.d0, 0.3d-3]
      prior(:,start_ind+6) = [-1.0d0, 10.0d0, 2.0652d-4, -1.d0] ! sed_ampl
      scale(start_ind+6,:) = [1.d-1, 1.d-4]
      prior(:,start_ind+7) = [1.1d3, 2.4d5, 1.05072985d4, -1.d0] ! sed_cutoff
      scale(start_ind+7,:) = [1.d3, 0.3d5]
      prior(:,start_ind+8) = [0.0d0, 4d0, 2.82996d-1, -1.d0] ! sed_b
      scale(start_ind+8,:) = [1.d-1, 1.d0]
      ! Component-specific parameters
      prior(:,start_ind+9) = [1.d0, 2.d0, 1.34d0, -1.d0] ! alpha
      scale(start_ind+9,:) = [1.d0, 0.02d0]
      prior(:,start_ind+10) = [3.d0, 5d0, 4.14d0, -1.d0] ! beta
      scale(start_ind+10,:) = [1.d0, 0.05d0]
      prior(:,start_ind+11) = [0.3d0, 1.1d0, 0.942d0, -1.d0] ! gamma
      scale(start_ind+11,:) = [1.d0, 0.03d0]
      prior(:,start_ind+12) = [0.1d0, 0.4d0, 0.189d0, -1.d0] ! mu
      scale(start_ind+12,:) = [1.d0, 0.013d0]
    end subroutine init_cloud_priors_and_scales

    subroutine init_band_priors_and_scales(self, start_ind, prior, scale)
      implicit none
      class(ZodiBand),            intent(in)    :: self
      integer(i4b),               intent(in)    :: start_ind
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:,1:), intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0,:) = [1.d-9, 0.2d-9]
      prior(:,start_ind+1) = [-30.d0, 30.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1,:) = [1.d0, 0.05d0]
      prior(:,start_ind+2) = [-720.d0, 720.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2,:) = [1.d0, 2.d0]
      prior(:,start_ind+3) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3,:) = [1.d0, 0.5d-3]
      prior(:,start_ind+4) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4,:) = [1.d0, 0.2d-3]
      prior(:,start_ind+5) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5,:) = [1.d0, 0.3d-3]
      prior(:,start_ind+6) = [-1.0d0, 10.0d0, 2.0652d-4, -1.d0] ! sed_ampl
      scale(start_ind+6,:) = [1.d-1, 1.d-4]
      prior(:,start_ind+7) = [1.1d3, 2.4d5, 1.05072985d4, -1.d0] ! sed_cutoff
      scale(start_ind+7,:) = [1.d3, 0.3d5]
      prior(:,start_ind+8) = [0.0d0, 4d0, 2.82996d-1, -1.d0] ! sed_b
      scale(start_ind+8,:) = [1.d-1, 1.d0]
      ! Component-specific parameters
      prior(:,start_ind+9) = [0.d0, 30d0, 0d0, -1.d0] ! delta_zeta
      scale(start_ind+9,:) = [1.d0, 0.14d0]
      prior(:,start_ind+10) = [0.8d0, 5.4d0, 4.14d0, -1.d0] ! delta_r
      scale(start_ind+10,:) = [1.d0, 0.005d0]
      !prior(:,start_ind+8) = [0.01d0, 1.5d0, 0.942d0, -1.d0] ! v
      prior(:,start_ind+11) = [0.01d0, 3.d0, 0.942d0, -1.d0] ! v
      scale(start_ind+11,:) = [1.d0, 0.1d0]      
      !prior(:,start_ind+9) = [3.99999d0, 4.000001d0, 0.189d0, -1.d0] ! p
      prior(:,start_ind+12) = [2.d0, 6.d0, 0.189d0, -1.d0] ! p
      scale(start_ind+12,:) = [1.d0, 1d-6]      
    end subroutine init_band_priors_and_scales

    subroutine init_ring_priors_and_scales(self, start_ind, prior, scale)
      implicit none
      class(ZodiRing),            intent(in)    :: self
      integer(i4b),               intent(in)    :: start_ind
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:,1:), intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0,:) = [1.d-9, 1.d-11]
      prior(:,start_ind+1) = [-30.d0, 30.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1,:) = [1.d0, 0.1d0]
      prior(:,start_ind+2) = [-720.d0, 720.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2,:) = [1.d0, 1.d0]
      prior(:,start_ind+3) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3,:) = [1.d0, 1d-3]
      prior(:,start_ind+4) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4,:) = [1.d0, 1d-3]
      prior(:,start_ind+5) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5,:) = [1.d0, 1d-3]
      prior(:,start_ind+6) = [-1.0d0, 10.0d0, 2.0652d-4, -1.d0] ! sed_ampl
      scale(start_ind+6,:) = [1.d-1, 1.d-4]
      prior(:,start_ind+7) = [1.1d3, 2.4d5, 1.05072985d4, -1.d0] ! sed_cutoff
      scale(start_ind+7,:) = [1.d3, 0.3d5]
      prior(:,start_ind+8) = [0.0d0, 4d0, 2.82996d-1, -1.d0] ! sed_b
      scale(start_ind+8,:) = [1.d-1, 1.d0]
      ! Component-specific parameters
      prior(:,start_ind+9) = [0.9d0, 1.1d0, 0d0, -1.d0] ! r
      scale(start_ind+9,:) = [1.d0, 0.01d0]
      prior(:,start_ind+10) = [0.d0, 0.3d0, 0.2d0, -1.d0] ! delta_r
      scale(start_ind+10,:) = [1.d0, 0.01d0]
      prior(:,start_ind+11) = [0.0d0, 0.2d0, 0.1d0, -1.d0] ! delta_z
      scale(start_ind+11,:) = [1.d0, 0.01d0]
      prior(:,start_ind+12) = [-60.d-3, 60.d-3, 0.d0, -1.d0] ! theta
      scale(start_ind+12,:) = [1.d0, 0.01d0]
      prior(:,start_ind+13) = [0.d0, 30.d0, 0.d0, -1.d0] ! sigma_theta
      scale(start_ind+13,:) = [1.d0, 0.01d0]
    end subroutine init_ring_priors_and_scales

    subroutine init_feature_priors_and_scales(self, start_ind, prior, scale)
      implicit none
      class(ZodiFeature),         intent(in)    :: self
      integer(i4b),               intent(in)    :: start_ind
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:,1:), intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0,:) = [1.d-9, 1.d-11]
      prior(:,start_ind+1) = [-30.d0, 30.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1,:) = [1.d0, 0.1d0]
      prior(:,start_ind+2) = [-720.d0, 720.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2,:) = [1.d0, 1.d0]
      prior(:,start_ind+3) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3,:) = [1.d0, 1d-3]
      prior(:,start_ind+4) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4,:) = [1.d0, 1d-3]      
      prior(:,start_ind+5) = [-0.001d0, 0.001d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5,:) = [1.d0, 1d-3]
      prior(:,start_ind+6) = [-1.0d0, 10.0d0, 2.0652d-4, -1.d0] ! sed_ampl
      scale(start_ind+6,:) = [1.d-1, 1.d-4]
      prior(:,start_ind+7) = [1.1d3, 2.4d5, 1.05072985d4, -1.d0] ! sed_cutoff
      scale(start_ind+7,:) = [1.d3, 0.3d5]
      prior(:,start_ind+8) = [0.0d0, 4d0, 2.82996d-1, -1.d0] ! sed_b
      scale(start_ind+8,:) = [1.d-1, 1.d0]
      ! Component-specific parameters
      prior(:,start_ind+9) = [0.9d0, 1.1d0, 0d0, -1.d0] ! r
      scale(start_ind+9,:) = [1.d0, 0.01d0]
      prior(:,start_ind+10) = [0.d0, 0.3d0, 0.2d0, -1.d0] ! delta_r
      scale(start_ind+10,:) = [1.d0, 0.01d0]
      prior(:,start_ind+11) = [0.0d0, 0.2d0, 0.1d0, -1.d0] ! delta_z
      scale(start_ind+11,:) = [1.d0, 0.01d0]
      prior(:,start_ind+12) = [-20.d0, 20.d0, 0.d0, -1.d0] ! theta
      scale(start_ind+12,:) = [1.d0, 0.01d0]
      prior(:,start_ind+13) = [0.d0, 30.d0, 0.d0, -1.d0] ! sigma_theta
      scale(start_ind+13,:) = [1.d0, 0.01d0]
    end subroutine init_feature_priors_and_scales

    subroutine init_interstellar_priors_and_scales(self, start_ind, prior, scale)
      implicit none
      class(ZodiInterstellar),    intent(in)    :: self
      integer(i4b),               intent(in)    :: start_ind
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:,1:),    intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0,:) = [1.d-9, 1.d-11]
      prior(:,start_ind+1) = [0.d0, 00.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1,:) = [1.d0, 0.d0]
      prior(:,start_ind+2) = [0.d0, 0.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2,:) = [1.d0, 0.d0]
      prior(:,start_ind+3) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3,:) = [1.d0, 0.d0]
      prior(:,start_ind+4) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4,:) = [1.d0, 0.d0]
      prior(:,start_ind+5) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5,:) = [1.d0, 0.d0]
      prior(:,start_ind+6) = [-1.0d0, 10.0d0, 2.0652d-4, -1.d0] ! sed_ampl
      scale(start_ind+6,:) = [1.d-1, 1.d-4]
      prior(:,start_ind+7) = [1.1d3, 2.4d5, 1.05072985d4, -1.d0] ! sed_cutoff
      scale(start_ind+7,:) = [1.d3, 0.3d5]
      prior(:,start_ind+8) = [0.0d0, 4d0, 2.82996d-1, -1.d0] ! sed_b
      scale(start_ind+8,:) = [1.d-1, 1.d0]
      ! Component-specific parameters
      prior(:,start_ind+9) = [0.d0, 0.0d0, 0d0, -1.d0] ! R, inactive
      scale(start_ind+9,:) = [1.d0, 0.d0]
      prior(:,start_ind+10) = [0.d0, 0.0d0, 0.2d0, -1.d0] ! alpha, inactive
      scale(start_ind+10,:) = [1.d0, 0.d0]
    end subroutine init_interstellar_priors_and_scales

    subroutine init_fan_priors_and_scales(self, start_ind, prior, scale)
      implicit none
      class(ZodiFan),             intent(in)    :: self
      integer(i4b),               intent(in)    :: start_ind
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:,1:), intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0,:) = [1.d-9, 1.d-11]
      prior(:,start_ind+1) = [0.d0, 10.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1,:) = [1.d0, 0.d0]
      prior(:,start_ind+2) = [-720.d0, 720.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2,:) = [1.d0, 0.d0]
      prior(:,start_ind+3) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3,:) = [1.d0, 0.d0]
      prior(:,start_ind+4) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4,:) = [1.d0, 0.d0]
      prior(:,start_ind+5) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5,:) = [1.d0, 0.d0]
      prior(:,start_ind+6) = [-1.0d0, 10.0d0, 2.0652d-4, -1.d0] ! sed_ampl
      scale(start_ind+6,:) = [1.d-1, 1.d-4]
      prior(:,start_ind+7) = [1.1d3, 2.4d5, 1.05072985d4, -1.d0] ! sed_cutoff
      scale(start_ind+7,:) = [1.d3, 0.3d5]
      prior(:,start_ind+8) = [0.0d0, 4d0, 2.82996d-1, -1.d0] ! sed_b
      scale(start_ind+8,:) = [1.d-1, 1.d0]
      ! Component-specific parameters
      prior(:,start_ind+9) = [5d0, 15d0, 0d0, -1.d0] ! Q
      scale(start_ind+9,:) = [1.d0, 0.d0]
      prior(:,start_ind+10) = [1.d0, 3d0, 0.2d0, -1.d0] ! P
      scale(start_ind+10,:) = [1.d0, 0.d0]
      prior(:,start_ind+11) = [0.5d0, 2d0, 0.1d0, -1.d0] ! Gamma
      scale(start_ind+11,:) = [1.d0, 0.d0]
      prior(:,start_ind+12) = [0.d0, 0.3d0, 0.d0, -1.d0] ! Z
      scale(start_ind+12,:) = [1.d0, 0.d0]
      prior(:,start_ind+13) = [1.d0, 5.d0, 0.d0, -1.d0] ! R_max
      scale(start_ind+13,:) = [1.d0, 0.d0]
    end subroutine init_fan_priors_and_scales

    subroutine init_comet_priors_and_scales(self, start_ind, prior, scale)
      implicit none
      class(ZodiComet),           intent(in)    :: self
      integer(i4b),               intent(in)    :: start_ind
      real(dp), dimension(1:,1:), intent(inout) :: prior
      real(dp), dimension(1:,1:), intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0,:) = [1.d-9, 0.d0]
      prior(:,start_ind+1) = [0.d0, 0.d0, 0.d0, -1.d0] ! Incl
      scale(start_ind+1,:) = [1.d0, 0.d0]
      prior(:,start_ind+2) = [0.d0, 0.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2,:) = [1.d0, 0.d0]
      prior(:,start_ind+3) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3,:) = [1.d0, 0.d0]
      prior(:,start_ind+4) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4,:) = [1.d0, 0.d0]
      prior(:,start_ind+5) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5,:) = [1.d0, 0.d0]
      prior(:,start_ind+6) = [-1.0d0, 10.0d0, 2.0652d-4, -1.d0] ! sed_ampl
      scale(start_ind+6,:) = [1.d-1, 1.d-4]
      prior(:,start_ind+7) = [1.1d3, 2.4d5, 1.05072985d4, -1.d0] ! sed_cutoff
      scale(start_ind+7,:) = [1.d3, 0.3d5]
      prior(:,start_ind+8) = [0.0d0, 4d0, 2.82996d-1, -1.d0] ! sed_b
      scale(start_ind+8,:) = [1.d-1, 1.d0]
      ! Component-specific parameters
      prior(:,start_ind+9) = [1d0, 5d0, 0d0, -1.d0] ! P
      scale(start_ind+9,:) = [1.d0, 0.d0]
      prior(:,start_ind+10) = [0.d0, 0.3d0, 0.2d0, -1.d0] ! z_mid
      scale(start_ind+10,:) = [1.d0, 0.d0]
      prior(:,start_ind+11) = [0.5d0, 1.5d0, 0.1d0, -1.d0] ! R_inner
      scale(start_ind+11,:) = [1.d0, 0.d0]
      prior(:,start_ind+12) = [1.5d0, 5.d0, 0.d0, -1.d0] ! R_outer
      scale(start_ind+12,:) = [1.d0, 0.d0]
    end subroutine init_comet_priors_and_scales

    ! See Appendix in Wright (1998) for details; https://iopscience.iop.org/article/10.1086/305345/pdf
   subroutine init_WrightCloudRing_priors_and_scales(self, start_ind, prior, scale)
     implicit none
     class(ZodiWrightCloudRing),           intent(in)    :: self
     integer(i4b),               intent(in)    :: start_ind
     real(dp), dimension(1:,1:), intent(inout) :: prior
     real(dp), dimension(1:,1:), intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-8, -1.d0] ! n_0
      scale(start_ind+0,:) = [1.d-9, 4.d-9]
      prior(:,start_ind+1) = [0.d0,0.d0, 0.d0, -1.d0] ! Incl -- don't fit in Wright model, it's part of the internal parameterization
      scale(start_ind+1,:) = [1.d0, 0.d0]
      prior(:,start_ind+2) = [0.d0, 0.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2,:) = [1.d0, 0.d0]
      prior(:,start_ind+3) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3,:) = [1.d0, 0.d0]
      prior(:,start_ind+4) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4,:) = [1.d0, 0.d0]
      prior(:,start_ind+5) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5,:) = [1.d0, 0.d0]
      prior(:,start_ind+6) = [-1.0d0, 10.0d0, 2.0652d-4, -1.d0] ! sed_ampl
      scale(start_ind+6,:) = [1.d-1, 1.d-4]
      prior(:,start_ind+7) = [1.1d3, 2.4d5, 1.05072985d4, -1.d0] ! sed_cutoff
      scale(start_ind+7,:) = [1.d3, 0.3d5]
      prior(:,start_ind+8) = [0.0d0, 4d0, 2.82996d-1, -1.d0] ! sed_b
      scale(start_ind+8,:) = [1.d-1, 1.d0]
      ! Component-specific parameters
      prior(:,start_ind+9) = [ 1.d0,   1.5d0,  1.2186d0, -1.d0] ! p1 - radial density exponent
      scale(start_ind+9,:) = [1.d0, 0.01d0]
      prior(:,start_ind+10) = [ 3.d0,   4d0,    3.6122d0, -1.d0] ! p3 - vertical "scale height"
      scale(start_ind+10,:) = [1.d0, 0.01d0]
      prior(:,start_ind+11) = [ 0.7d0,  1.1d0,  0.9285d0, -1.d0] ! p4 - vertical density exponent
      scale(start_ind+11,:) = [1.d0, 0.01d0]
      prior(:,start_ind+12) = [-1.8d0, -1.2d0, -1.4766d0, -1.d0] ! p5 - ln(sin i) at break
      scale(start_ind+12,:) = [1.d0, 0.01d0] 
      prior(:,start_ind+13) = [ -1.d0,  1.d0,  0.3705d0, -1.d0] ! p6 - 10 x cloud pole x component
      scale(start_ind+13,:) = [1.d0, 0.01d0]
      prior(:,start_ind+14) = [ -1.d0,  1.d0, -0.0736d0, -1.d0] ! p7 - 10 x cloud pole y component
      scale(start_ind+14,:) = [1.d0, 0.01d0]
      prior(:,start_ind+15) = [ -1.d0,  1.d0, -0.0235d0, -1.d0] ! p8 - 10 x cloud offset x component
      scale(start_ind+15,:) = [1.d0, 0.01d0]
      prior(:,start_ind+16) = [-1.d0,   1.d0, -0.0081d0, -1.d0] ! p9 - 10 x cloud offset y component
      scale(start_ind+16,:) = [1.d0, 0.01d0]
      prior(:,start_ind+17) = [ 0.1d0,  3.d0,  0.7548d0, -1.d0] ! p10 - 10 x density contrast of Dermott ring
      scale(start_ind+17,:) = [1.d0, 0.01d0]
      prior(:,start_ind+18) = [ 0.1d0,  1.d0,  0.4284d0, -1.d0] ! p13 - "dimple" in Dermott ring
      scale(start_ind+18,:) = [1.d0, 0.01d0]
      prior(:,start_ind+19) = [ 0.d0,  50.d0, 27.7741d0, -1.d0] ! p14 - vertical scale for Dermott ring
      scale(start_ind+19,:) = [1.d0, 1.d0]
      prior(:,start_ind+20) = [-0.1d0, 0.1d0, -0.0251d0, -1.d0] ! p15 - spherical term in vertical density
      scale(start_ind+20,:) = [1.d0, 0.01d0]
      prior(:,start_ind+21) = [-0.2d0, 0.2d0,  0.0249d0, -1.d0] ! p16 - (sin i)**2 term in vertical density
      scale(start_ind+21,:) = [1.d0, 0.01d0]
      prior(:,start_ind+22) = [-0.2d0, 0.2d0, -0.0456d0, -1.d0] ! p17 - Additional density at sin i ~ 0.5
      scale(start_ind+22,:) = [1.d0, 0.01d0]
      prior(:,start_ind+23) = [-0.2d0, 0.2d0, -0.1276d0, -1.d0] ! p18 - Additional density at sin i ~ 0.25
      scale(start_ind+23,:) = [1.d0, 0.01d0]
      prior(:,start_ind+24) = [-0.2d0, 0.2d0, -0.0103d0, -1.d0] ! p19 - Additional density at sin i ~ 0.17
      scale(start_ind+24,:) = [1.d0, 0.01d0]
    end subroutine init_WrightCloudRing_priors_and_scales

    ! See Appendix in Wright (1998) for details; https://iopscience.iop.org/article/10.1086/305345/pdf
   subroutine init_WrightBand_priors_and_scales(self, start_ind, prior, scale)
     implicit none
     class(ZodiWrightBand),     intent(in)    :: self
     integer(i4b),               intent(in)    :: start_ind
     real(dp), dimension(1:,1:), intent(inout) :: prior
     real(dp), dimension(1:,1:), intent(inout) :: scale

      ! Common parameters
      prior(:,start_ind+0) = [1.d-11, 1.d-5, 1.d-9, -1.d0] ! n_0
      scale(start_ind+0,:) = [1.d-9, 4.d-9]
      prior(:,start_ind+1) = [0.d0,0.d0, 0.d0, -1.d0] ! Incl -- don't fit in Wright model, these are part of the internal parameterization
      scale(start_ind+1,:) = [1.d0, 0.d0]
      prior(:,start_ind+2) = [0.d0, 0.d0, 0.d0, -1.d0] ! Omega
      scale(start_ind+2,:) = [1.d0, 0.d0]
      prior(:,start_ind+3) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3,:) = [1.d0, 0.d0]
      prior(:,start_ind+4) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4,:) = [1.d0, 0.d0]
      prior(:,start_ind+5) = [0.d0, 0.d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5,:) = [1.d0, 0.d0]
      prior(:,start_ind+6) = [-1.0d0, 10.0d0, 2.0652d-4, -1.d0] ! sed_ampl
      scale(start_ind+6,:) = [1.d-1, 1.d-4]
      prior(:,start_ind+7) = [1.1d3, 2.4d5, 1.05072985d4, -1.d0] ! sed_cutoff
      scale(start_ind+7,:) = [1.d3, 0.3d5]
      prior(:,start_ind+8) = [0.0d0, 4d0, 2.82996d-1, -1.d0] ! sed_b
      scale(start_ind+8,:) = [1.d-1, 1.d0]
      ! Component-specific parameters
      prior(:,start_ind+9) = [ 1.d0,   2.0d0,  1.3849d0, -1.d0] ! q1 - 10 x (sin i)_max for band 1
      scale(start_ind+9,:) = [1.d0, 0.01d0]
      prior(:,start_ind+10) = [-0.3d0, 0.3d0,   0.1735d0, -1.d0] ! q5 - 10 x band pole x component
      scale(start_ind+10,:) = [1.d0, 0.01d0]
      prior(:,start_ind+11) = [-0.3d0, 0.3d0,  -0.2088d0, -1.d0] ! q6 - 10 x band pole y component 
      scale(start_ind+11,:) = [1.d0, 0.01d0]
      prior(:,start_ind+12) = [-2.d0, 2.d0,    -1.5723d0, -1.d0] ! q7 - 10 x band offset x component
      scale(start_ind+12,:) = [1.d0, 0.1d0] 
      prior(:,start_ind+13) = [ -2.d0, 2.d0,  -0.2225d0, -1.d0] ! q8 - 10 x band offset y component
      scale(start_ind+13,:) = [1.d0, 0.01d0]
      prior(:,start_ind+14) = [ 2.d0,  4.d0,   3.14d0,   -1.d0] ! R_1 - Outer radius
      scale(start_ind+14,:) = [1.d0, 0.01d0]
    end subroutine init_WrightBand_priors_and_scales
    
    
    subroutine param2model_cloud(self, x)
      implicit none
      class(ZodiCloud),                intent(inout) :: self
      real(dp),         dimension(1:), intent(in)    :: x
      self%n_0          = x(1)
      self%incl         = x(2)
      self%Omega        = x(3)
      self%x_0          = x(4)
      self%y_0          = x(5)
      self%z_0          = x(6)
      self%sed_ampl     = x(7)
      self%sed_cutoff   = x(8)
      self%sed_b        = x(9)
      self%alpha        = x(10)
      self%beta         = x(11)
      self%gamma        = x(12)
      self%mu           = x(13)       
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
      x(7)  = self%sed_ampl 
      x(8)  = self%sed_cutoff  
      x(9)  = self%sed_b 
      x(10) = self%alpha  
      x(11) = self%beta      
      x(12) = self%gamma    
      x(13) = self%mu           
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
      self%sed_ampl   = x(7)
      self%sed_cutoff = x(8)
      self%sed_b      = x(9)
      self%delta_zeta = x(10)
      self%delta_r    = x(11)
      self%v          = x(12)
      self%p          = x(13)     
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
      x(7)  = self%sed_ampl   
      x(8)  = self%sed_cutoff 
      x(9)  = self%sed_b      
      x(10) = self%delta_zeta 
      x(11) = self%delta_r    
      x(12) = self%v          
      x(13) = self%p             
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
      self%sed_ampl    = x(7)
      self%sed_cutoff  = x(8)
      self%sed_b       = x(9)
      self%R_0         = x(10)
      self%sigma_r     = x(11)
      self%sigma_z     = x(12)
      self%theta_0     = x(13)
      self%sigma_theta = x(14)      
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
      x(7)  = self%sed_ampl    
      x(8)  = self%sed_cutoff  
      x(9)  = self%sed_b       
      x(10) = self%R_0         
      x(11) = self%sigma_r     
      x(12) = self%sigma_z      
      x(13) = self%theta_0     
      x(14) = self%sigma_theta 
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
      self%sed_ampl    = x(7)
      self%sed_cutoff  = x(8)
      self%sed_b       = x(9)
      self%R_0         = x(10)
      self%sigma_r     = x(11)
      self%sigma_z     = x(12)
      self%theta_0     = x(13)
      self%sigma_theta = x(14)   
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
      x(7)  = self%sed_ampl    
      x(8)  = self%sed_cutoff  
      x(9)  = self%sed_b       
      x(10) = self%R_0         
      x(11) = self%sigma_r     
      x(12) = self%sigma_z      
      x(13) = self%theta_0     
      x(14) = self%sigma_theta 
    end subroutine model2param_feature

    subroutine param2model_interstellar(self, x)
      implicit none
      class(ZodiInterstellar),                intent(inout) :: self
      real(dp),                dimension(1:), intent(in)    :: x
      self%n_0          = x(1)
      self%incl         = x(2)
      self%Omega        = x(3)
      self%x_0          = x(4)
      self%y_0          = x(5)
      self%z_0          = x(6)
      self%sed_ampl     = x(7)
      self%sed_cutoff   = x(8)
      self%sed_b        = x(9)
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
      x(7)  = self%sed_ampl    
      x(8)  = self%sed_cutoff  
      x(9)  = self%sed_b  
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
      self%sed_ampl     = x(7)
      self%sed_cutoff   = x(8)
      self%sed_b        = x(9)
      self%Q            = x(10)
      self%P            = x(11)
      self%gamma        = x(12)
      self%Z_midplane_0 = x(13)
      self%R_outer      = x(14)
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
      x(7)  = self%sed_ampl    
      x(8)  = self%sed_cutoff  
      x(9)  = self%sed_b  
      x(10)  = self%Q 
      x(11)  = self%P  
      x(12)  = self%gamma 
      x(13) = self%Z_midplane_0
      x(14) = self%R_outer
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
      self%sed_ampl     = x(7)
      self%sed_cutoff   = x(8)
      self%sed_b        = x(9)
      self%P            = x(10)
      self%Z_midplane_0 = x(11)
      self%R_inner      = x(12)
      self%R_outer      = x(13)      
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
      x(7)  = self%sed_ampl    
      x(8)  = self%sed_cutoff  
      x(9)  = self%sed_b   
      x(10)  = self%P 
      x(11)  = self%Z_midplane_0  
      x(12)  = self%R_inner
      x(13) = self%R_outer    
    end subroutine model2param_comet

    subroutine param2model_WrightCloudRing(self, x)
      implicit none
      class(ZodiWrightCloudRing),                intent(inout) :: self
      real(dp),                   dimension(1:), intent(in)    :: x
      self%n_0          = x(1)
      self%incl         = x(2)
      self%Omega        = x(3)
      self%x_0          = x(4)
      self%y_0          = x(5)
      self%z_0          = x(6)
      self%sed_ampl     = x(7)
      self%sed_cutoff   = x(8)
      self%sed_b        = x(9)
      self%p1           = x(10)
      self%p3           = x(11)
      self%p4           = x(12)
      self%p5           = x(13)
      self%p6           = x(14)
      self%p7           = x(15)
      self%p8           = x(16)
      self%p9           = x(17)
      self%p10          = x(18)
      self%p13          = x(19)
      self%p14          = x(20)
      self%p15          = x(21)
      self%p16          = x(22)
      self%p17          = x(23)
      self%p18          = x(24)
      self%p19          = x(25)
    end subroutine param2model_WrightCloudRing

    subroutine model2param_WrightCloudRing(self, x)
      implicit none
      class(ZodiWrightCloudRing),                intent(in)  :: self
      real(dp),                   dimension(1:), intent(out) :: x
      x(1)  = self%n_0  
      x(2)  = self%incl 
      x(3)  = self%Omega 
      x(4)  = self%x_0   
      x(5)  = self%y_0   
      x(6)  = self%z_0   
      x(7)  = self%sed_ampl    
      x(8)  = self%sed_cutoff  
      x(9)  = self%sed_b   
      x(10)  = self%p1
      x(11)  = self%p3
      x(12)  = self%p4
      x(13) = self%p5
      x(14) = self%p6
      x(15) = self%p7
      x(16) = self%p8
      x(17) = self%p9
      x(18) = self%p10
      x(19) = self%p13
      x(20) = self%p14
      x(21) = self%p15
      x(22) = self%p16
      x(23) = self%p17
      x(24) = self%p18
      x(25) = self%p19
    end subroutine model2param_WrightCloudRing

    subroutine param2model_WrightBand(self, x)
      implicit none
      class(ZodiWrightBand),                intent(inout) :: self
      real(dp),              dimension(1:), intent(in)    :: x
      self%n_0          = x(1)
      self%incl         = x(2)
      self%Omega        = x(3)
      self%x_0          = x(4)
      self%y_0          = x(5)
      self%z_0          = x(6)
      self%sed_ampl     = x(7)
      self%sed_cutoff   = x(8)
      self%sed_b        = x(9)
      self%q1           = x(10)
      self%q5           = x(11)
      self%q6           = x(12)
      self%q7           = x(13)
      self%q8           = x(14)
      self%R_1          = x(15)
    end subroutine param2model_WrightBand

    subroutine model2param_WrightBand(self, x)
      implicit none
      class(ZodiWrightBand),                intent(in)  :: self
      real(dp),              dimension(1:), intent(out) :: x
      x(1)  = self%n_0  
      x(2)  = self%incl 
      x(3)  = self%Omega 
      x(4)  = self%x_0   
      x(5)  = self%y_0   
      x(6)  = self%z_0
      x(7)  = self%sed_ampl    
      x(8)  = self%sed_cutoff  
      x(9)  = self%sed_b   
      x(10)  = self%q1
      x(11)  = self%q5
      x(12)  = self%q6
      x(13) = self%q7
      x(14) = self%q8
      x(15) = self%R_1
    end subroutine model2param_WrightBand
    
    
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

    subroutine get_density_WrightCloudRing(self, X_vec, theta, n_out)
      class(ZodiWrightCloudRing) :: self
      real(dp), dimension(:, :), intent(in) :: X_vec
      real(dp), intent(in) :: theta
      real(dp), dimension(:), intent(out) :: n_out
      integer(i4b) :: i
      real(dp) :: R, z_c(3), o_c(3), R_c
      real(dp) :: f, sin_i, Z, S
      real(dp) :: x_D, y_D, z_D, L_D, A, D

      ! X_vec = r = position to evaluate model in heliocentric coordinates
!!$      open(58,file='n.dat', recl=1024)
!!$      write(*,*) 'p3', self%p3
!!$      write(*,*) 'p4', self%p4
!!$      write(*,*) 'p5', self%p5
!!$      write(*,*) 'p6', self%p6
!!$      write(*,*) 'p7', self%p7
!!$      write(*,*) 'p8', self%p8
!!$      write(*,*) 'p9', self%p9
!!$      write(*,*) 'p10', self%p10
!!$      write(*,*) 'p13', self%p13
!!$      write(*,*) 'p14', self%p14
!!$      write(*,*) 'p15', self%p15
!!$      write(*,*) 'p16', self%p16
!!$      write(*,*) 'p17', self%p17
!!$      write(*,*) 'p18', self%p18
!!$      write(*,*) 'p19', self%p19
      do i = 1, size(n_out)
         R     = sqrt(sum(X_vec(:,i)**2))
         z_c   = [self%p6, self%p7, 10.d0] / sqrt(100.d0 + self%p6**2 + self%p7**2)
         o_c   = [self%p8/10.d0, self%p9/10.d0, 0.d0]
         sin_i = sum(z_c*X_vec(:,i)) / R
         R_c   = R + sum(o_c * X_vec(:,i))

         ! Cloud
         S     = exp(self%p5)
         if (abs(sin_i) > S) then
            Z = abs(sin_i) - 0.5d0*S
         else
            Z = 0.5d0*sin_i**2/S
         end if
         f = exp(-self%p3 * Z**self%p4) + self%p15 + self%p16 * sin_i**2 &
              & + self%p17 *  4.d0*sin_i**2*exp( -4.d0*sin_i**2) &
              & + self%p18 * 16.d0*sin_i**2*exp(-16.d0*sin_i**2) &
              & + self%p19 * 36.d0*sin_i**2*exp(-36.d0*sin_i**2)         
         
         ! Ring
         x_D = X_vec(1,i)*cos(theta) + X_vec(2,i)*sin(theta)
         y_D = X_vec(2,i)*cos(theta) - X_vec(1,i)*sin(theta)
         z_D = X_vec(3,i)
         L_D = abs(atan2(y_D, x_D) + 0.25d0)
         if (L_D < 0.375d0) then
            A = cos(8.d0*pi*L_d/3.d0)
         else if (L_D < 0.75d0) then
            A  = 0.5d0*(cos(8.d0*pi*L_D/3)-1.d0)
         else
            A = 0.d0
         end if
         D   = exp(-56.5d0*(sqrt(x_D**2 + y_D**2) - 1.133d0 + 0.133d0 * self%p13 * exp(-4.d0*L_D**2))**2 - self%p14 * z_D**2/R**2)
         
         ! Total density
         n_out(i) = self%n_0 * R/R_c * f * R_c**(-self%p1) * (1.d0 + 0.1d0 * self%p10 * D*(1.d0+A))
         !write(58,*) R, n_out(i), X_vec(:,i), sin_i, Z, f, self%n_0, R_c, self%p1, self%p10, D, A
      end do
!      close(58)
!      stop
    end subroutine get_density_WrightCloudRing

    subroutine get_density_WrightBand(self, X_vec, theta, n_out)
      class(ZodiWrightBand) :: self
      real(dp), dimension(:, :), intent(in) :: X_vec
      real(dp), intent(in) :: theta
      real(dp), dimension(:), intent(out) :: n_out
      integer(i4b) :: i
      real(dp) :: R, z_b(3), o_b(3), R_b, sin_i

      do i = 1, size(n_out)
         R     = sqrt(sum(X_vec(:,i)**2))
         z_b   = [self%q5, self%q6, 10.d0] / sqrt(100.d0 + self%q5**2 + self%q6**2)
         sin_i = sum(z_b*X_vec(:,i)) / R
         o_b   = [self%q7/10.d0, self%q8/10.d0, 0.d0]
         R_b   = R + sum(o_b * X_vec(:,i))

         if (abs(sin_i) < 0.1d0*self%q1 .and. R_b < self%R_1) then
            !n_out(i) = self%n_0 * R/R_b**2 * cosh(1.72d0*abs(sin_i)/(0.1d0*self%q1))
            n_out(i) = self%n_0 * R/R_b**2 * cosh(1.72d0*abs(sin_i)/self%q1) / pi
         else
            n_out(i) = 0.d0
         end if
      end do
    end subroutine get_density_WrightBand

  end module comm_zodi_comp_mod
