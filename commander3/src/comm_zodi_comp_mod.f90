module comm_zodi_comp_mod
   use comm_param_mod
   use comm_hdf_mod
   implicit none

   type, abstract :: ZodiComponent
      real(dp) :: x_0, y_0, z_0, incl, Omega, n_0
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
      prior(:,start_ind+3) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! X_0
      scale(start_ind+3,:) = [1.d0, 1d-3]
      prior(:,start_ind+4) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! Y_0
      scale(start_ind+4,:) = [1.d0, 0.8d-3]
      prior(:,start_ind+5) = [-0.02d0, 0.02d0, 0.d0, -1.d0] ! ! Z_0
      scale(start_ind+5,:) = [1.d0, 0.3d-3]
      ! Component-specific parameters
      prior(:,start_ind+6) = [1.d0, 2.d0, 1.34d0, -1.d0] ! alpha
      scale(start_ind+6,:) = [1.d0, 0.02d0]
      prior(:,start_ind+7) = [3.d0, 5d0, 4.14d0, -1.d0] ! beta
      scale(start_ind+7,:) = [1.d0, 0.05d0]
      prior(:,start_ind+8) = [0.3d0, 1.1d0, 0.942d0, -1.d0] ! gamma
      scale(start_ind+8,:) = [1.d0, 0.03d0]
      prior(:,start_ind+9) = [0.1d0, 0.4d0, 0.189d0, -1.d0] ! mu
      scale(start_ind+9,:) = [1.d0, 0.013d0]
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
      ! Component-specific parameters
      prior(:,start_ind+6) = [0.d0, 30d0, 0d0, -1.d0] ! delta_zeta
      scale(start_ind+6,:) = [1.d0, 0.14d0]
      prior(:,start_ind+7) = [0.8d0, 5.4d0, 4.14d0, -1.d0] ! delta_r
      scale(start_ind+7,:) = [1.d0, 0.005d0]
      prior(:,start_ind+8) = [0.01d0, 1.5d0, 0.942d0, -1.d0] ! v
      scale(start_ind+8,:) = [1.d0, 0.1d0]      
      prior(:,start_ind+9) = [3.99999d0, 4.000001d0, 0.189d0, -1.d0] ! p
      scale(start_ind+9,:) = [1.d0, 1d-6]      
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
      ! Component-specific parameters
      prior(:,start_ind+6) = [0.9d0, 1.1d0, 0d0, -1.d0] ! r
      scale(start_ind+6,:) = [1.d0, 0.01d0]
      prior(:,start_ind+7) = [0.d0, 0.3d0, 0.2d0, -1.d0] ! delta_r
      scale(start_ind+7,:) = [1.d0, 0.01d0]
      prior(:,start_ind+8) = [0.0d0, 0.2d0, 0.1d0, -1.d0] ! delta_z
      scale(start_ind+8,:) = [1.d0, 0.01d0]
      prior(:,start_ind+9) = [-60.d-3, 60.d-3, 0.d0, -1.d0] ! theta
      scale(start_ind+9,:) = [1.d0, 0.01d0]
      prior(:,start_ind+10) = [0.d0, 30.d0, 0.d0, -1.d0] ! sigma_theta
      scale(start_ind+10,:) = [1.d0, 0.01d0]
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
      ! Component-specific parameters
      prior(:,start_ind+6) = [0.9d0, 1.1d0, 0d0, -1.d0] ! r
      scale(start_ind+6,:) = [1.d0, 0.01d0]
      prior(:,start_ind+7) = [0.d0, 0.3d0, 0.2d0, -1.d0] ! delta_r
      scale(start_ind+7,:) = [1.d0, 0.01d0]
      prior(:,start_ind+8) = [0.0d0, 0.2d0, 0.1d0, -1.d0] ! delta_z
      scale(start_ind+8,:) = [1.d0, 0.01d0]
      prior(:,start_ind+9) = [-20.d0, 20.d0, 0.d0, -1.d0] ! theta
      scale(start_ind+9,:) = [1.d0, 0.01d0]
      prior(:,start_ind+10) = [0.d0, 30.d0, 0.d0, -1.d0] ! sigma_theta
      scale(start_ind+10,:) = [1.d0, 0.01d0]
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
      ! Component-specific parameters
      prior(:,start_ind+6) = [0.d0, 0.0d0, 0d0, -1.d0] ! R, inactive
      scale(start_ind+6,:) = [1.d0, 0.d0]
      prior(:,start_ind+7) = [0.d0, 0.0d0, 0.2d0, -1.d0] ! alpha, inactive
      scale(start_ind+7,:) = [1.d0, 0.d0]
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
      ! Component-specific parameters
      prior(:,start_ind+6) = [5d0, 15d0, 0d0, -1.d0] ! Q
      scale(start_ind+6,:) = [1.d0, 0.d0]
      prior(:,start_ind+7) = [1.d0, 3d0, 0.2d0, -1.d0] ! P
      scale(start_ind+7,:) = [1.d0, 0.d0]
      prior(:,start_ind+8) = [0.5d0, 2d0, 0.1d0, -1.d0] ! Gamma
      scale(start_ind+8,:) = [1.d0, 0.d0]
      prior(:,start_ind+9) = [0.d0, 0.3d0, 0.d0, -1.d0] ! Z
      scale(start_ind+9,:) = [1.d0, 0.d0]
      prior(:,start_ind+10) = [1.d0, 5.d0, 0.d0, -1.d0] ! R_max
      scale(start_ind+10,:) = [1.d0, 0.d0]
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
      ! Component-specific parameters
      prior(:,start_ind+6) = [1d0, 5d0, 0d0, -1.d0] ! P
      scale(start_ind+6,:) = [1.d0, 0.d0]
      prior(:,start_ind+7) = [0.d0, 0.3d0, 0.2d0, -1.d0] ! z_mid
      scale(start_ind+7,:) = [1.d0, 0.d0]
      prior(:,start_ind+8) = [0.5d0, 1.5d0, 0.1d0, -1.d0] ! R_inner
      scale(start_ind+8,:) = [1.d0, 0.d0]
      prior(:,start_ind+9) = [1.5d0, 5.d0, 0.d0, -1.d0] ! R_outer
      scale(start_ind+9,:) = [1.d0, 0.d0]
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
    
  end module comm_zodi_comp_mod
