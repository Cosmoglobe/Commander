module comm_zodi_mod
   use comm_utils
   use comm_param_mod
   use spline_1D_mod
   use comm_hdf_mod
   use hashtbl

   implicit none

   private
   public initialize_zodi_mod, get_s_zodi, get_s_zodi_comp, ZodiModel, ZodiComponent, zodi_model


   type, abstract :: ZodiComponent
      real(dp) :: x_0, y_0, z_0, incl, Omega, n_0
      real(dp), allocatable :: sin_omega, cos_omega, sin_incl, cos_incl
      !priors
   contains
      procedure :: init_base_comp
      procedure(init_interface), deferred :: init
      procedure(density_interface), deferred :: get_density
   end type ZodiComponent

   type :: ZodiComponentContainer
      class(ZodiComponent), allocatable :: c
      character(len=32), allocatable :: labels(:)
   end type

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

   end interface

   type, extends(ZodiComponent) :: ZodiCloud
      real(dp) :: alpha, beta, gamma, mu
   contains
      procedure :: init => init_cloud
      procedure :: get_density => get_density_cloud
   end type ZodiCloud

   type, extends(ZodiComponent) :: ZodiBand
      real(dp) :: delta_zeta, delta_r, v, p
      real(dp) :: delta_zeta_rad = 0.d0
   contains
      procedure :: init => init_band
      procedure :: get_density => get_density_band
   end type ZodiBand

   type, extends(ZodiComponent) :: ZodiRing
      real(dp) :: R_0, sigma_r, sigma_z, theta_0, sigma_theta
      real(dp) :: theta_0_rad = 0.d0
      real(dp) :: sigma_theta_rad = 0.d0
   contains
      procedure :: init => init_ring
      procedure :: get_density => get_density_ring
   end type ZodiRing

   type, extends(ZodiComponent) :: ZodiFeature
      real(dp) :: R_0, sigma_r, sigma_z, theta_0, sigma_theta
      real(dp) :: theta_0_rad = 0.d0
      real(dp) :: sigma_theta_rad = 0.d0
   contains
      procedure :: init => init_feature
      procedure :: get_density => get_density_feature
   end type ZodiFeature

   type, extends(ZodiComponent) :: ZodiInterstellar
      real(dp)  :: R, alpha
   contains
      procedure :: init => init_interstellar
      procedure :: get_density => get_density_interstellar
   end type ZodiInterstellar

   type :: ZodiModel
      class(ZodiComponentContainer), allocatable :: comps(:)
      character(len=128), allocatable :: comp_labels(:), general_labels(:)
      integer(i4b) :: n_comps, n_params, n_common_params, n_general_params
      real(dp) :: T_0, delta
      real(dp), dimension(10) :: F_sun = [2.3405606d8, 1.2309874d8, 64292872d0, 35733824d0, 5763843d0, 1327989.4d0, 230553.73d0, 82999.336d0, 42346.605d0, 14409.608d0] * 1d-20 ! convert to specific intensity units
      real(dp), dimension(10) :: C0 = [-0.94209999, -0.52670002, -0.4312, 0., 0., 0., 0., 0., 0., 0.]
      real(dp), dimension(10) :: C1 = [0.1214, 0.18719999, 0.1715, 0., 0., 0., 0., 0., 0., 0.]
      real(dp), dimension(10) :: C2 = [-0.1648, -0.59829998, -0.63330001, 0., 0., 0., 0., 0., 0., 0.]
   contains
      procedure :: init_comps, init_general_params, params_to_model, model_to_params, model_to_chain, comp_from_chain
   end type ZodiModel
   type(ZodiModel), target :: zodi_model

contains
   subroutine initialize_zodi_mod(cpar)
      type(comm_params), intent(in) :: cpar
      integer(i4b) :: i, ierr
      character(len=256) :: file_path
      real(dp), allocatable :: comp_params(:, :)
      character(len=128), allocatable :: comp_labels(:)

      ! Set model and zodi_mod parameters from cpar
      zodi_model%n_comps = cpar%zs_ncomps
      allocate(comp_params(zodi_model%n_comps, size(cpar%zs_comp_params, dim=1)))
      zodi_model%general_labels = cpar%zodi_param_labels%general
      zodi_model%comp_labels = cpar%zs_comp_labels(1:zodi_model%n_comps)
      zodi_model%n_common_params = size(cpar%zodi_param_labels%common)
      zodi_model%n_general_params = size(cpar%zodi_param_labels%general)
      
      comp_params = cpar%zs_comp_params(:, :, 1)
      do i = 1, zodi_model%n_comps
         if (trim(adjustl(cpar%zs_init_hdf(i))) /= 'none') then 
            call zodi_model%comp_from_chain(cpar, comp_params, i)
         end if
      end do
      call zodi_model%init_comps(comp_params, cpar%zs_comp_types, cpar%zodi_param_labels)
      call zodi_model%init_general_params(cpar%zs_general_params(:, 1))
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
      self%sin_incl = sin(self%incl*deg2rad)
      self%cos_incl = cos(self%incl*deg2rad)
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

         if (abs(theta_prime) < self%sigma_theta_rad) then
            n_out(i) = 0.d0
         else
            n_out(i) = self%n_0*exp(term1 - term2)
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
         self%comps(i)%labels = [param_labels%common]
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
            self%comps(i)%labels = [self%comps(i)%labels, param_labels%cloud]
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
            self%comps(i)%labels = [self%comps(i)%labels, param_labels%band]
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
            self%comps(i)%labels = [self%comps(i)%labels, param_labels%ring]
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
            self%comps(i)%labels = [self%comps(i)%labels, param_labels%feature]
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
            self%comps(i)%labels = [self%comps(i)%labels, param_labels%interstellar]
         case default
            print *, 'Invalid zodi component type in zodi `init_from_params`:', trim(adjustl(comp_types(i)))
            stop
         end select
         call self%comps(i)%c%init()
         self%n_params = self%n_params + size(self%comps(i)%labels)
      end do
   end subroutine

   subroutine model_to_params(self, x, labels)
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
   end subroutine model_to_params

   subroutine params_to_model(self, x)
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
         end select
         call self%comps(i)%c%init()
      end do
      self%T_0 = x(running_idx + 1)
      self%delta = x(running_idx + 2)
   end subroutine params_to_model

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
      call self%model_to_params(params, labels)
      
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
   end subroutine

   subroutine comp_from_chain(self, cpar, params, comp_idx)
      ! Initialize a component from a chain
      class(ZodiModel), target, intent(inout) :: self
      type(comm_params), intent(in) :: cpar
      real(dp), intent(inout) :: params(:, :)
      integer(i4b) :: comp_idx

      logical(lgt) :: exist
      integer(i4b) :: i, j, l, ierr, initsamp
      character(len=6) :: itext

      type(hdf_file) :: file

      character(len=32), allocatable :: param_labels(:)
      character(len=512) :: chainfile, group_name

      if (cpar%myid == cpar%root) then
         call get_chainfile_and_samp(trim(cpar%zs_init_hdf(comp_idx)), chainfile, initsamp)
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

   subroutine get_s_zodi(s_therm, s_scat, s_zodi, emissivity, albedo, alpha)
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
      ! alpha : optional
      !     Scale factor per component
      real(sp), dimension(:, :), intent(in) :: s_scat, s_therm
      real(sp), dimension(:), intent(out)   :: s_zodi
      real(dp), dimension(:), intent(in) :: emissivity, albedo
      real(dp), dimension(:), intent(in), optional :: alpha
      integer(i4b) :: i, n_comps

      n_comps = size(emissivity)
      s_zodi = 0.
      do i = 1, n_comps
         if (present(alpha)) then
            call get_s_zodi_comp(s_therm(:, i), s_scat(:, i), s_zodi, emissivity(i), albedo(i), alpha(i))
         else
            call get_s_zodi_comp(s_therm(:, i), s_scat(:, i), s_zodi, emissivity(i), albedo(i))
         end if
      end do
   end subroutine get_s_zodi

   subroutine get_s_zodi_comp(s_therm_comp, s_scat_comp, s_zodi_comp, emissivity_comp, albedo_comp, alpha_comp)
      ! Evaluates the zodiacal signal (eq. 20 in ZodiPy paper [k98 model]) given
      ! integrated thermal zodiacal emission and scattered zodiacal light for a single
      ! component.
      !
      ! Parameters:
      ! -----------
      ! s_scat_comp :
      !     Integrated contribution from scattered sunlight light.
      ! s_therm_comp :
      !     Integrated contribution from thermal interplanetary dust emission.
      ! s_zodi :
      !     Zodiacal signal.
      ! albedo_comp :
      !     Albedo of the zodiacal component.
      ! emissivity_comp :
      !     Emissivity of the zodiacal component.
      ! alpha_comp : optional
      !     Scale factor for a component
      real(sp), dimension(:), intent(in) :: s_scat_comp, s_therm_comp
      real(sp), dimension(:), intent(inout) :: s_zodi_comp
      real(dp), intent(in) :: emissivity_comp, albedo_comp
      real(dp), intent(in), optional :: alpha_comp
      integer(i4b) :: i, n_comps
      if (present(alpha_comp)) then 
         s_zodi_comp = s_zodi_comp + ((s_scat_comp * albedo_comp) + (1. - albedo_comp) * emissivity_comp * s_therm_comp) * alpha_comp
      else 
         s_zodi_comp = s_zodi_comp + ((s_scat_comp * albedo_comp) + (1. - albedo_comp) * emissivity_comp * s_therm_comp)
      end if
   end subroutine get_s_zodi_comp
end module comm_zodi_mod
