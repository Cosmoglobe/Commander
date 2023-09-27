module comm_zodi_mod
   use comm_utils
   use comm_param_mod
   use spline_1D_mod
   use comm_hdf_mod
   use hashtbl

   implicit none

   private
   public initialize_zodi_mod, get_s_zodi, ZodiModel, ZodiComponent, zodi_model

   ! Global variables
   real(dp) :: gauss_degree

   type, abstract :: ZodiComponent
      real(dp) :: x_0, y_0, z_0, incl, Omega, n_0, t_0, delta
      real(dp), allocatable :: sin_omega, cos_omega, sin_incl, cos_incl
      !priors
   contains
      procedure :: initialize_base
      procedure(initialize_interface), deferred :: initialize
      procedure(density_interface), deferred :: get_density
   end type ZodiComponent

   type :: ZodiComponentContainer
      class(ZodiComponent), allocatable :: c
      character(len=32), allocatable :: labels(:)
   end type

   abstract interface
      subroutine initialize_interface(self)
         import dp, ZodiComponent
         class(ZodiComponent)  :: self
      end subroutine initialize_interface

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
      procedure :: initialize => initialize_cloud
      procedure :: get_density => get_density_cloud
   end type ZodiCloud

   type, extends(ZodiComponent) :: ZodiBand
      real(dp) :: delta_zeta, delta_r, v, p
      real(dp) :: delta_zeta_rad = 0.d0
   contains
      procedure :: initialize => initialize_band
      procedure :: get_density => get_density_band
   end type ZodiBand

   type, extends(ZodiComponent) :: ZodiRing
      real(dp) :: R_0, sigma_r, sigma_z
   contains
      procedure :: initialize => initialize_ring
      procedure :: get_density => get_density_ring
   end type ZodiRing

   type, extends(ZodiComponent) :: ZodiFeature
      real(dp) :: R_0, sigma_r, sigma_z, theta_0, sigma_theta
      real(dp) :: theta_0_rad = 0.d0
      real(dp) :: sigma_theta_rad = 0.d0
   contains
      procedure :: initialize => initialize_feature
      procedure :: get_density => get_density_feature
   end type ZodiFeature

   ! Stores the interplanetary dust model parameters
   type :: ZodiModel
      integer(i4b) :: n_comps, n_params, n_common_params
      real(dp) :: T_0, delta
      class(ZodiComponentContainer), allocatable :: comps(:)
      character(len=128), allocatable :: comp_labels(:)
      real(dp), dimension(10) :: F_sun = [2.3405606d8, 1.2309874d8, 64292872d0, 35733824d0, 5763843d0, 1327989.4d0, 230553.73d0, 82999.336d0, 42346.605d0, 14409.608d0]
      real(dp), dimension(10) :: C0 = [-0.94209999, -0.52670002, -0.4312, 0., 0., 0., 0., 0., 0., 0.]
      real(dp), dimension(10) :: C1 = [0.1214, 0.18719999, 0.1715, 0., 0., 0., 0., 0., 0., 0.]
      real(dp), dimension(10) :: C2 = [-0.1648, -0.59829998, -0.63330001, 0., 0., 0., 0., 0., 0., 0.]
   contains
      procedure :: init_comps, params_to_model, model_to_params, model_to_chain, model_from_chain, model_to_ascii, ascii_to_model
   end type ZodiModel

   ! Global zodi parameter object
   type(ZodiModel), target :: zodi_model

contains
   subroutine initialize_zodi_mod(cpar)
      ! Initialize the zodi module.
      !
      ! Parameters
      ! ----------
      ! cpar: comm_params
      !    Parameter file variables.

      type(comm_params), intent(in) :: cpar
      integer(i4b) :: i, ierr
      character(len=256) :: file_path
      real(dp), allocatable :: param_vec(:)
      character(len=128), allocatable :: comp_labels(:)

      zodi_model%n_comps = cpar%zs_ncomps
      allocate (zodi_model%comp_labels(zodi_model%n_comps))
      zodi_model%comp_labels = cpar%zs_comp_labels(1:zodi_model%n_comps)

      if (cpar%zs_init_chain /= 'none') then
         ! init from model as defined in the paramterfile/defaults
         call zodi_model%model_from_chain(cpar)
      else
         ! init from model structure as defined in the paramterfile/defaults with values from chain
         call zodi_model%init_comps(cpar%zs_comp_params(:, :, 1), cpar%zs_comp_types, cpar%zodi_param_labels)
      end if

      gauss_degree = cpar%zs_los_steps
      ! Currently only support single values temperature
      zodi_model%T_0 = zodi_model%comps(1)%c%t_0
      zodi_model%delta = zodi_model%comps(1)%c%delta

      zodi_model%n_common_params = size(cpar%zodi_param_labels%common)
      zodi_model%n_params = zodi_model%n_comps*size(cpar%zodi_param_labels%common)
      do i = 1, zodi_model%n_comps
         select case (cpar%zs_comp_types(i))
         case ('cloud')
            zodi_model%n_params = zodi_model%n_params + size(cpar%zodi_param_labels%cloud)
         case ('band')
            zodi_model%n_params = zodi_model%n_params + size(cpar%zodi_param_labels%band)
         case ('ring')
            zodi_model%n_params = zodi_model%n_params + size(cpar%zodi_param_labels%ring)
         case ('feature')
            zodi_model%n_params = zodi_model%n_params + size(cpar%zodi_param_labels%feature)
         end select
      end do
      ! call zodi_model%ascii_to_model(cpar, "/mn/stornext/u3/metins/dirbe/chains/chains_testing/init_zodi.dat")
      ! call zodi_model%model_to_ascii(cpar, "/mn/stornext/u3/metins/dirbe/chains/chains_testing/init_zodi.dat")
   end subroutine initialize_zodi_mod

   subroutine get_s_zodi(s_therm, s_scat, s_zodi, emissivity, albedo, comp)
      ! Evaluates the zodiacal signal (eq. 20 in zodipy paper [k98 model]) given
      ! integrated thermal zodiacal emission and scattered zodiacal light.
      !
      ! Parameters:
      ! -----------
      ! emissivity :
      !     Emissivity of the zodiacal components.
      ! albedo :
      !     Albedo of the zodiacal components.
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
      ! comp :
      !     If present, only the component with index comp is used to compute the zodiacal signal.

      real(sp), dimension(:, :), intent(in) :: s_scat, s_therm
      real(sp), dimension(:), intent(inout) :: s_zodi
      real(dp), dimension(:), intent(in) :: emissivity, albedo
      integer(i4b), intent(in), optional :: comp

      integer(i4b) :: i, n_comps

      n_comps = size(emissivity)
      if (present(comp)) then
         s_zodi = s_scat(:, comp)*albedo(comp) + (1.-albedo(comp))*emissivity(comp)*s_therm(:, comp)
      else
         s_zodi = 0.
         do i = 1, n_comps
            s_zodi = s_zodi + s_scat(:, i)*albedo(i) + (1.-albedo(i))*emissivity(i)*s_therm(:, i)
         end do
      end if
   end subroutine get_s_zodi

   subroutine initialize_base(self)
      class(ZodiComponent) :: self
      self%sin_omega = sin(self%Omega*deg2rad)
      self%cos_omega = cos(self%Omega*deg2rad)
      self%sin_incl = sin(self%incl*deg2rad)
      self%cos_incl = cos(self%incl*deg2rad)
   end subroutine initialize_base

   subroutine initialize_cloud(self)
      class(ZodiCloud) :: self
      call self%initialize_base()
   end subroutine initialize_cloud

   subroutine initialize_band(self)
      class(ZodiBand) :: self
      self%delta_zeta_rad = self%delta_zeta*deg2rad
      call self%initialize_base()
   end subroutine initialize_band

   subroutine initialize_ring(self)
      class(ZodiRing) :: self
      call self%initialize_base()
   end subroutine initialize_ring

   subroutine initialize_feature(self)
      class(ZodiFeature) :: self
      self%theta_0_rad = self%theta_0*deg2rad
      self%sigma_theta_rad = self%sigma_theta*deg2rad
      call self%initialize_base()
   end subroutine initialize_feature

   subroutine get_density_cloud(self, X_vec, theta, n_out)
      class(ZodiCloud) :: self
      real(dp), dimension(:, :), intent(in) :: X_vec
      real(dp), intent(in) :: theta
      real(dp), dimension(:), intent(out) :: n_out
      integer(i4b) :: i
      real(dp) :: R, Z_midplane, zeta, g, x_prime, y_prime, z_prime

      do i = 1, gauss_degree
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
      real(dp) :: x_prime, y_prime, z_prime, R, Z_midplane, zeta, zeta_over_delta_zeta, term1, term2, term3, term4

      do i = 1, gauss_degree
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
         term3 = 1.d0 + (zeta_over_delta_zeta**self%p)/self%v
         term4 = 1.d0 - exp(-((R/self%delta_r)**20))

         n_out(i) = term1*term2*term3*term4
      end do
   end subroutine get_density_band

   subroutine get_density_ring(self, X_vec, theta, n_out)
      class(ZodiRing) :: self
      real(dp), dimension(:, :), intent(in)  :: X_vec
      real(dp), intent(in) :: theta
      real(dp), dimension(:), intent(out) :: n_out
      integer(i4b) :: i
      real(dp) :: x_prime, y_prime, z_prime, R, Z_midplane, term1, term2

      do i = 1, gauss_degree
         x_prime = X_vec(1, i) - self%x_0
         y_prime = X_vec(2, i) - self%y_0
         z_prime = X_vec(3, i) - self%z_0

         R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
         Z_midplane = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl
         term1 = -((R - self%R_0)**2)/self.sigma_r**2
         term2 = abs(Z_midplane/self.sigma_z)

         n_out(i) = self%n_0*exp(term1 - term2)
      end do
   end subroutine get_density_ring

   subroutine get_density_feature(self, X_vec, theta, n_out)
      class(ZodiFeature) :: self
      real(dp), dimension(:, :), intent(in) :: X_vec
      real(dp), intent(in) :: theta
      real(dp), dimension(:), intent(out) :: n_out
      integer(i4b) :: i
      real(dp) :: x_prime, y_prime, z_prime, R, Z_midplane, theta_prime, exp_term

      do i = 1, gauss_degree
         x_prime = X_vec(1, i) - self%x_0
         y_prime = X_vec(2, i) - self%y_0
         z_prime = X_vec(3, i) - self%z_0
         theta_prime = atan2(y_prime, x_prime) - (theta + self%theta_0_rad)

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

   subroutine init_comps(self, params, comp_types, param_labels)
      class(ZodiModel), target, intent(inout) :: self
      character(len=*), intent(in), optional :: comp_types(:)
      class(InterplanetaryDustParamLabels) :: param_labels
      real(dp) :: params(:, :)
      integer(i4b) :: i, ierr

      allocate (self%comps(self%n_comps))
      
      ! NOTE: order of `zs_comp_params` is important and is given by the lable arrays in `InterplanetaryDustParamLabels` in param_mod
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
                & t_0=params(i, 7), &
                & delta=params(i, 8), &
                & alpha=params(i, 9), &
                & beta=params(i, 10), &
                & gamma=params(i, 11), &
                & mu=params(i, 12) &
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
                & t_0=params(i, 7), &
                & delta=params(i, 8), &
                & delta_zeta=params(i, 9), &
                & delta_r=params(i, 10), &
                & v=params(i, 11), &
                & p=params(i, 12) &
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
               & t_0=params(i, 7), &
               & delta=params(i, 8), &
               & R_0=params(i, 9), &
               & sigma_r=params(i, 10), &
               & sigma_z=params(i, 11) &
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
               & t_0=params(i, 7), &
               & delta=params(i, 8), &
               & R_0=params(i, 9), &
               & sigma_r=params(i, 10), &
               & sigma_z=params(i, 11), &
               & theta_0=params(i, 12), &
               & sigma_theta=params(i, 13) &
            &)
            self%comps(i)%labels = [self%comps(i)%labels, param_labels%feature]
         case default
            print *, 'Invalid zodi component type in zodi `init_from_params`:', trim(adjustl(comp_types(i)))
            stop
         end select
         call self%comps(i)%c%initialize()
      end do
   end subroutine

   subroutine model_to_params(self, x, labels)
      ! Dumps a zodi model to a parameter vector `x`.
      ! If `labels` is present, it is populated with the corresponding labels.
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
         x(running_idx + 7) = self%comps(i)%c%t_0
         x(running_idx + 8) = self%comps(i)%c%delta
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
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiFeature)
            x(running_idx + 1) = comp%R_0
            x(running_idx + 2) = comp%sigma_r
            x(running_idx + 3) = comp%sigma_z
            x(running_idx + 4) = comp%theta_0
            x(running_idx + 5) = comp%sigma_theta
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         end select
         if (present(labels)) then
            labels_copy = self%comps(i)%labels
            call upcase(self%comp_labels(i), comp_label_upper(i))
            do j = 1, size(labels_copy)
               labels_copy(j) = trim(adjustl(comp_label_upper(i)))//'_'//trim(adjustl(labels_copy(j))) 
            end do
               labels = [labels, labels_copy]
         end if
      end do
   end subroutine model_to_params

   subroutine params_to_model(self, x)
      ! Dumps a zodi model to a parameter vector
      class(ZodiModel), intent(inout) :: self
      real(dp), intent(in) :: x(:)
      integer(i4b) :: i, running_idx

      if (size(x) /= self%n_params) stop "Error: argument 'x' has the wrong size. must be `size(zodi_model%n_params)`"

      running_idx = 0
      do i = 1, self%n_comps
         ! The order of these operations much match the order tabulated in the labels in `InterplanetaryDustParamLabels`
         self%comps(i)%c%n_0 = x(running_idx + 1)
         self%comps(i)%c%incl = x(running_idx + 2)
         self%comps(i)%c%Omega = x(running_idx + 3)
         self%comps(i)%c%x_0 = x(running_idx + 4)
         self%comps(i)%c%y_0 = x(running_idx + 5)
         self%comps(i)%c%z_0 = x(running_idx + 6)
         self%comps(i)%c%t_0 = x(running_idx + 7)
         self%comps(i)%c%delta = x(running_idx + 8)
         running_idx = running_idx + self%n_common_params
         select type (comp => self%comps(i)%c)
         class is (ZodiCloud)
            comp%alpha = x(running_idx + 1)
            comp%beta = x(running_idx + 2)
            comp%gamma = x(running_idx + 3)
            comp%mu = x(running_idx + 4)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiBand)
            comp%delta_zeta = x(running_idx + 1)
            comp%delta_r = x(running_idx + 2)
            comp%v = x(running_idx + 3)
            comp%p = x(running_idx + 4)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiRing)
            comp%R_0 = x(running_idx + 1)
            comp%sigma_r = x(running_idx + 2)
            comp%sigma_z = x(running_idx + 3)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         class is (ZodiFeature)
            comp%R_0 = x(running_idx + 1)
            comp%sigma_r = x(running_idx + 2)
            comp%sigma_z = x(running_idx + 3)
            comp%theta_0 = x(running_idx + 4)
            comp%sigma_theta = x(running_idx + 5)
            running_idx = running_idx + size(self%comps(i)%labels) - self%n_common_params
         end select
         call self%comps(i)%c%initialize()
      end do
   end subroutine params_to_model


   subroutine model_to_chain(self, cpar, iter)
      ! Writes the zodi model to an hdf file
      class(ZodiModel), intent(in) :: self
      type(comm_params), intent(in) :: cpar
      integer(i4b), intent(in) :: iter

      integer(i4b) :: i, j, hdferr, ierr, unit, n_comp_params, param_idx
      logical(lgt) :: exist, init, new_header
      character(len=6) :: itext
      character(len=4) :: ctext
      character(len=512) :: zodi_path, comp_path, comp_group_path, param_path, chainfile, hdfpath, param_label
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

      comp_group_path = trim(adjustl(zodi_path))//'/params'
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
      call close_hdf_file(file)
   end subroutine

   subroutine model_from_chain(self, cpar)
      class(ZodiModel), target, intent(inout) :: self
      type(comm_params), intent(in) :: cpar
      logical(lgt) :: exist
      integer(i4b) :: i, j, l, unit, ierr, initsamp, hdferr
      character(len=6) :: itext

      type(hdf_file) :: file
      real(dp) :: comp_params(100), params(100, 100)
      character(len=32), allocatable :: common_param_labels(:), param_labels(:)
      character(len=512) :: chainfile
      TYPE(h5o_info_t) :: object_info

      params = 0.
      if (cpar%myid == cpar%root) then
         call get_chainfile_and_samp(trim(cpar%zs_init_chain), chainfile, initsamp)
         inquire (file=trim(chainfile), exist=exist)
         if (.not. exist) call report_error('Zodi init chain does not exist = '//trim(chainfile))
         l = len(trim(chainfile))
         if (.not. ((trim(chainfile(l-2:l)) == '.h5') .or. (trim(chainfile(l-3:l)) == '.hd5'))) call report_error('Zodi init chain must be a .h5 file')
         
         call open_hdf_file(trim(chainfile), file, "r")
         
         call int2string(initsamp, itext)
         do i = 1, cpar%zs_ncomps
            call h5eset_auto_f(0, hdferr)
            call h5oget_info_by_name_f(file%filehandle, trim(adjustl(itext)//'/zodi/params/'//trim(adjustl(cpar%zs_comp_labels(i)))), object_info, hdferr)
            if (hdferr /= 0) cycle
            param_labels = cpar%zodi_param_labels%get_labels(trim(adjustl(cpar%zs_comp_types(i))), add_common=.true.)
            do j = 1, size(param_labels)
               call read_hdf(file, trim(adjustl(itext)//'/zodi/params/'//trim(adjustl(cpar%zs_comp_labels(i)))// &
                   & '/'//trim(adjustl(param_labels(j)))), comp_params(j))
            end do
               params(i, :) = comp_params
         end do 
         call close_hdf_file(file)
      end if
      call mpi_bcast(params, size(params, dim=1) * size(params, dim=2), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
      call self%init_comps(params, cpar%zs_comp_types, cpar%zodi_param_labels)
   end subroutine model_from_chain


   subroutine model_to_ascii(self, cpar, filename)
      class(ZodiModel), target, intent(in) :: self
      type(comm_params), intent(in) :: cpar
      character(len=*), intent(in) :: filename

      integer(i4b) :: io, i, j, running_idx
      logical(lgt) :: exists
      real(dp), allocatable :: params(:)
      integer(i4b), allocatable :: comp_switch_indices(:)
      character(len=128), allocatable :: labels(:)

      if (cpar%myid_chain /= cpar%root) return
      inquire(file=trim(adjustl(filename)), exist=exists)
      if (.not. exists) then
         print *, "zodi asciifile: " // trim(adjustl(filename)) // " does not exist"
         stop
      end if

      open(newunit=io, file=trim(adjustl(filename)), action="write")

      allocate(params(self%n_params))
      call self%model_to_params(params, labels=labels)

      allocate(comp_switch_indices(self%n_comps))

      running_idx = 0
      do i = 1, self%n_comps
         running_idx = running_idx + size(self%comps(i)%labels)
         comp_switch_indices(i) = running_idx
      end do

      do i = 1, self%n_params

         if (any(comp_switch_indices == i)) then
               write(io, fmt='(a, T25, ES12.5, a)') trim(adjustl(labels(i))), params(i), new_line('a')
            else
               write(io, fmt='(a, T25, ES12.5)') trim(adjustl(labels(i))), params(i)
         end if
      end do

      close(io)
   end subroutine

   subroutine ascii_to_model(self, cpar, filename)
      class(ZodiModel), target, intent(inout) :: self
      type(comm_params), intent(in) :: cpar
      character(len=*), intent(in) :: filename
      type(hash_tbl_sll) :: htbl

      integer(i4b) :: i, io, io_status, ierr
      logical(lgt) :: exists
      character(len=128) :: key, val, line
      character(len=128), allocatable :: labels(:)
      real(dp), allocatable :: params(:)

      allocate(params(self%n_params))
      if (cpar%myid_chain == cpar%root) then
         inquire(file=trim(adjustl(filename)), exist=exists)
         if (.not. exists) then
            print *, "zodi asciifile: " // trim(adjustl(filename)) // " does not exist"
            stop
         end if
         
         call init_hash_tbl_sll(htbl, tbl_len=500)

         open(newunit=io, file=trim(adjustl(filename)), action="read")
         io_status = 0
         do while (io_status == 0)
            read(io, "(a)", iostat=io_status) line
            if (io_status == 0 .and. line /= "") then
               read(line, *) key, val
               call tolower(key)
               call put_hash_tbl_sll(htbl, trim(key), trim(val)) 
            end if
         end do
         close(io)

         call self%model_to_params(params, labels)
         params = 0.
         if (size(labels) /= size(params)) stop "Error: size of labels and params do not match"
         do i = 1, size(labels)
            call get_parameter_hashtable(htbl, labels(i), par_dp=params(i))
         end do
      end if

      call mpi_bcast(params, size(params), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
      call self%params_to_model(params)
   end subroutine

end module comm_zodi_mod
