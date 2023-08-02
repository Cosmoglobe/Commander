module comm_zodi_mod
    use comm_utils
    use comm_param_mod
    use spline_1D_mod
    use comm_hdf_mod
    implicit none

    private
    public initialize_zodi_mod, get_s_zodi, ZodiModel, ZodiComponent, zodi_model

    ! Global variables
    integer(i4b) :: gauss_degree, n_interp_points
    real(dp) :: R_cutoff
    logical(lgt) :: use_cloud, use_band1, use_band2, use_band3, use_ring, use_feature
    
    type, abstract :: ZodiComponent
        real(dp) :: x_0, y_0, z_0, incl, Omega, n_0
        real(dp), allocatable :: sin_omega, cos_omega, sin_incl, cos_incl
    contains
        procedure(initialize_interface), deferred :: initialize
        procedure(density_interface), deferred :: get_density
    end type ZodiComponent

    type :: ZodiComponentContainer
        class(ZodiComponent), allocatable :: c
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
        integer(i4b) :: n_comps, n_bands
        integer(i4b) :: n_parameters
        real(dp) :: T_0, delta ! solar system temperature parameters
        real(dp), allocatable, dimension(:) :: nu_ref, solar_irradiance ! spectral parameters
        real(dp), allocatable, dimension(:, :) :: phase_coeffs, emissivity, albedo ! spectral parameters
        class(ZodiComponentContainer), allocatable :: comps(:)
        character(len=10), allocatable :: comp_labels(:)
        type(spline_type) :: solar_irradiance_spl! spline interpolators
        type(spline_type), allocatable, dimension(:) :: phase_coeff_spl
    contains
        procedure :: init_from_defaults, init_from_chain, build_splines, param_vec_to_model, model_to_param_vec, output_to_hd5
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

        call initialize_hyper_parameters(cpar)
        if (cpar%zs_init_chain == 'none') then
            call zodi_model%init_from_defaults(cpar)
        else
            call zodi_model%init_from_chain(cpar)
        end if
    end subroutine initialize_zodi_mod

    subroutine get_s_zodi(emissivity, albedo, s_therm, s_scat, s_zodi)
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
        real(dp), dimension(:), intent(in) :: emissivity, albedo     
        real(sp), dimension(:, :), intent(in) :: s_scat, s_therm
        real(sp), dimension(:), intent(out) :: s_zodi
        integer(i4b) :: i, n_comps

        s_zodi = 0.
        n_comps = size(s_scat, dim=2)
        do i = 1, n_comps
            s_zodi = s_zodi +  s_scat(:, i) * albedo(i) + (1. - albedo(i)) * emissivity(i) * s_therm(:, i)
        end do
    end subroutine get_s_zodi

   ! Methods for initizializing the zodiacal components
    ! -----------------------------------------------------------------------------------
    subroutine initialize_cloud(self)
        class(ZodiCloud) :: self
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_cloud

    subroutine initialize_band(self)
        class(ZodiBand) :: self
        self%delta_zeta_rad = self%delta_zeta * deg2rad
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_band

    subroutine initialize_ring(self)
        class(ZodiRing) :: self
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_ring

    subroutine initialize_feature(self)
        class(ZodiFeature) :: self
        self%theta_0_rad = self%theta_0 * deg2rad
        self%sigma_theta_rad = self%sigma_theta * deg2rad
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_feature

    ! Methods describing the allocatable, densitry distribution of the zodiacal allocatable, components
    ! -----------------------------------------------------------------------------------
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
                g = (zeta * zeta) / (2.d0 * self%mu)
            else
                g = zeta - (0.5d0 * self%mu)
            end if

            n_out(i) = self%n_0 * R**(-self%alpha) * exp(-self%beta * g**self%gamma)
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

            zeta_over_delta_zeta = zeta / self%delta_zeta_rad
            term1 = (3.d0 * self%n_0) / R
            term2 = exp(-(zeta_over_delta_zeta**6))

            ! Differs from eq 8 in K98 by a factor of 1/self.v. See Planck XIV
            ! section 4.1.2.
            term3 = 1.d0 + (zeta_over_delta_zeta**self%p) / self%v
            term4 = 1.d0 - exp(-((R / self%delta_r) ** 20))

            n_out(i) = term1 * term2 * term3 * term4
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
            term1 = -((R - self%R_0) ** 2) / self.sigma_r**2
            term2 = abs(Z_midplane/self.sigma_z)

            n_out(i) = self%n_0 * exp(term1 - term2)
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

            exp_term = ((R - self%R_0) ** 2 / self%sigma_r**2) + (abs(Z_midplane) / self%sigma_z) + (theta_prime**2 / self%sigma_theta_rad**2)

            n_out(i) = self%n_0 * exp(-exp_term)
        end do
    end subroutine get_density_feature

    subroutine initialize_hyper_parameters(cpar)
        ! Initialize global hyper parameters.
        !
        ! Parameters
        ! ----------
        ! cpar: comm_params
        !    Parameter file variables.

        type(comm_params), intent(in) :: cpar

        use_cloud = cpar%zs_use_cloud
        use_band1 = cpar%zs_use_band1
        use_band2 = cpar%zs_use_band2
        use_band3 = cpar%zs_use_band3
        use_ring = cpar%zs_use_ring
        use_feature = cpar%zs_use_feature
        R_cutoff = cpar%zs_los_cut
        gauss_degree = cpar%zs_gauss_quad_order
    end subroutine initialize_hyper_parameters

    subroutine init_from_defaults(self, cpar)
        ! Initialize the kelsal et al 1998 model.
        !
        ! Parameters
        ! ----------
        ! self: zodi_params_k98
        !    Zodi model object.
        ! cpar: comm_params
        !    Parameter file variables.

        class(ZodiModel), target, intent(inout) :: self
        type(comm_params), intent(in) :: cpar
        integer(i4b) :: i

        ! aux
        self%n_bands = cpar%zs_nbands
        self%n_comps = cpar%zs_ncomps
        self%n_parameters = 62 + 2*(self%n_comps * self%n_bands)

        allocate(self%comps(self%n_comps))
        allocate(self%comp_labels(self%n_comps))
        allocate(self%emissivity(self%n_comps, self%n_bands))
        allocate(self%albedo(self%n_comps, self%n_bands))
        ! Tempereature parameters
        self%T_0 = cpar%zs_t_0
        self%delta = cpar%zs_delta

        ! Spectral parameters and splines
        allocate(self%nu_ref, mold=cpar%zs_nu_ref)
        allocate(self%phase_coeffs, mold=cpar%zs_phase_coeff)
        allocate(self%solar_irradiance, mold=cpar%zs_solar_irradiance)
        self%nu_ref = cpar%zs_nu_ref
        self%phase_coeffs = cpar%zs_phase_coeff
        self%solar_irradiance = cpar%zs_solar_irradiance
        self%albedo = cpar%zs_albedo
        self%emissivity = cpar%zs_emissivity
        allocate(self%phase_coeff_spl(3))
        call self%build_splines()

        i = 1
        ! Initialize class instances
        ! if the allocate(ZodiCloud::self%comps(i)%c) notation is confusing, see this
        ! https://fortran-lang.discourse.group/t/class-array-with-different-types-at-each-index/2627
        if (use_cloud) then
            allocate(ZodiCloud::self%comps(i)%c)
            self%comps(i)%c = ZodiCloud(x_0=cpar%zs_common(1, 1), y_0=cpar%zs_common(1, 2), z_0=cpar%zs_common(1, 3), &
                                    incl=cpar%zs_common(1, 4), Omega=cpar%zs_common(1, 5), n_0=cpar%zs_common(1, 6), & 
                                    alpha=cpar%zs_cloud_alpha, beta=cpar%zs_cloud_beta, gamma=cpar%zs_cloud_gamma, &
                                    mu=cpar%zs_cloud_mu)
            self%comp_labels(i) = "cloud"
            i = i + 1
        end if
        if (use_band1) then
            allocate(ZodiBand::self%comps(i)%c)
            self%comps(i)%c = ZodiBand(x_0=cpar%zs_common(2, 1), y_0=cpar%zs_common(2, 2), z_0=cpar%zs_common(2, 3), &
                              incl=cpar%zs_common(2, 4), Omega=cpar%zs_common(2, 5), n_0=cpar%zs_common(2, 6), & 
                              delta_zeta=cpar%zs_bands_delta_zeta(1), delta_r=cpar%zs_bands_delta_r(1), & 
                              v=cpar%zs_bands_v(1), p=cpar%zs_bands_p(1))
            self%comp_labels(i) = "band1"
            i = i + 1
        end if
        if (use_band2) then
            allocate(ZodiBand::self%comps(i)%c)
            self%comps(i)%c = ZodiBand(x_0=cpar%zs_common(3, 1), y_0=cpar%zs_common(3, 2), z_0=cpar%zs_common(3, 3), &
                              incl=cpar%zs_common(3, 4), Omega=cpar%zs_common(3, 5), n_0=cpar%zs_common(3, 6), & 
                              delta_zeta=cpar%zs_bands_delta_zeta(2), delta_r=cpar%zs_bands_delta_r(2), & 
                              v=cpar%zs_bands_v(2), p=cpar%zs_bands_p(2))
            self%comp_labels(i) = "band2"
            i = i + 1
        end if
        if (use_band3) then
            allocate(ZodiBand::self%comps(i)%c)
            self%comps(i)%c = ZodiBand(x_0=cpar%zs_common(4, 1), y_0=cpar%zs_common(4, 2), z_0=cpar%zs_common(4, 3), &
                              incl=cpar%zs_common(4, 4), Omega=cpar%zs_common(4, 5), n_0=cpar%zs_common(4, 6), & 
                              delta_zeta=cpar%zs_bands_delta_zeta(3), delta_r=cpar%zs_bands_delta_r(3), & 
                              v=cpar%zs_bands_v(3), p=cpar%zs_bands_p(3))
            self%comp_labels(i) = "band3"
            i = i + 1
        end if
        if (use_ring) then
            allocate(ZodiRing::self%comps(i)%c)
            self%comps(i)%c = ZodiRing(x_0=cpar%zs_common(5, 1), y_0=cpar%zs_common(5, 2), z_0=cpar%zs_common(5, 3), &
                             incl=cpar%zs_common(5, 4), Omega=cpar%zs_common(5, 5), n_0=cpar%zs_common(5, 6), &
                             R_0=cpar%zs_ring_r, sigma_r=cpar%zs_ring_sigma_r, sigma_z=cpar%zs_ring_sigma_z)
            self%comp_labels(i) = "ring"
            i = i + 1
        end if
        if (use_feature) then
            allocate(ZodiFeature::self%comps(i)%c)
            self%comps(i)%c = ZodiFeature(x_0=cpar%zs_common(6, 1), y_0=cpar%zs_common(6, 2), z_0=cpar%zs_common(6, 3), &
                                   incl=cpar%zs_common(6, 4), Omega=cpar%zs_common(6, 5), n_0=cpar%zs_common(6, 6), &
                                   R_0=cpar%zs_feature_r, sigma_r=cpar%zs_feature_sigma_r, sigma_z=cpar%zs_feature_sigma_z, &
                                   theta_0=cpar%zs_feature_theta, sigma_theta=cpar%zs_feature_sigma_theta)
            self%comp_labels(i) = "feature"
            i = i + 1
        end if

        ! Run initializer on comps
        do i = 1, self%n_comps
            call self%comps(i)%c%initialize()
        end do        
    end subroutine init_from_defaults

    subroutine init_from_chain(self, cpar)
        ! Initialize the kelsal et al 1998 model.
        !
        ! Parameters
        ! ----------
        ! self: zodi_params_k98
        !    Zodi model object.
        ! cpar: comm_params
        !    Parameter file variables.

        class(ZodiModel), target, intent(inout) :: self
        type(comm_params), intent(in) :: cpar
        logical(lgt)        :: exist
        integer(i4b) :: i, l,  unit, ierr, initsamp
            character(len=6)          :: itext, itext2

        type(hdf_file) :: file
        real(dp) :: param_vec(self%n_parameters)
        character(len=512) :: chainfile


        unit = getlun()
        
        ! validate that the chain file exists, and get path to file and which sample to use
        call get_chainfile_and_samp(trim(cpar%zs_init_chain), chainfile, initsamp)
        inquire(file=trim(chainfile), exist=exist)
        if (.not. exist) call report_error('Zodi init chain does not exist = ' // trim(chainfile))
        l = len(trim(chainfile))
        if (.not. ((trim(chainfile(l-2:l)) == '.h5') .or. (trim(chainfile(l-3:l)) == '.hd5'))) call report_error('Zodi init chain must be a .h5 file')
        
        call open_hdf_file(trim(chainfile), file, "r")

        call int2string(initsamp, itext)

        call read_hdf(file, itext//"/zodi/params" , param_vec)
        call close_hdf_file(file)

        call self%init_from_defaults(cpar)
        call self%param_vec_to_model(param_vec)

    end subroutine init_from_chain


    subroutine build_splines(self)
        ! Build splines for the phase function, and solar irradiance
        ! for each component in the model
        class(ZodiModel), intent(inout) :: self
        integer(i4b) :: i

        do i = 1, 3
            call spline_simple(self%phase_coeff_spl(i), self%nu_ref, self%phase_coeffs(:, i), regular=.false.)
        end do
        call spline_simple(self%solar_irradiance_spl, self%nu_ref, self%solar_irradiance, regular=.false.)
    end subroutine build_splines

    subroutine param_vec_to_model(self, x)
        ! Updates model parameters given a parametor vector x
        class(ZodiModel), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        integer(i4b) :: i, j, k

        if (size(x) /= self%n_parameters) stop "param_vec has the wrong size"
        do i = 1, self%n_comps
            select type (a => self%comps(i)%c)
            class is (ZodiCloud)
                a%n_0 = x(1)
                a%incl = x(2)
                a%Omega = x(3)
                a%x_0 = x(4)
                a%y_0 = x(5)
                a%z_0 = x(6)
                a%alpha = x(7)
                a%beta = x(8)
                a%gamma = x(9)
                a%mu = x(10)
            class is (ZodiBand)
                select case (self%comp_labels(i))
                case ('band1')
                    j = 10
                case ('band2')
                    j = 20
                case ('band3')
                    j = 30
                case default
                    stop 'Invalid component label'
                end select
                a%n_0 = x(j + 1)
                a%incl = x(j + 2)
                a%Omega = x(j + 3)
                a%x_0 = x(j + 4)
                a%y_0 = x(j + 5)
                a%z_0 = x(j + 6)
                a%delta_zeta = x(j + 7)
                a%delta_r = x(j + 8)
                a%v = x(j + 9)
                a%p = x(j + 10)
            class is (ZodiRing)
                a%n_0 = x(41)
                a%incl = x(42)
                a%Omega = x(43)
                a%x_0 = x(44)
                a%y_0 = x(45)
                a%z_0 = x(46)
                a%R_0 = x(47)
                a%sigma_r = x(48)
                a%sigma_z = x(49)
            class is (ZodiFeature)
                a%n_0 = x(50)
                a%incl = x(51)
                a%Omega = x(52)
                a%x_0 = x(53)
                a%y_0 = x(54)
                a%z_0 = x(55)
                a%R_0 = x(56)
                a%sigma_r = x(57)
                a%sigma_z = x(58)
                a%theta_0 = x(59)
                a%sigma_theta = x(60)
            class default 
                stop 'Invalid zodi class'
            end select
            call self%comps(i)%c%initialize()
        end do
        self%T_0 = x(61)
        self%delta = x(62)

        !emissivities
        k = 63
        do i = 1, self%n_comps
            do j = 1, self%n_bands
                self%emissivity(i, j) = x(k)
                k = k + 1
            end do
        end do

        !albedos
        k = 63 + self%n_comps * self%n_bands
        do i = 1, self%n_comps
            do j = 1, self%n_bands
                self%albedo(i, j) = x(k)
                k = k + 1
            end do
        end do
    end subroutine

    function model_to_param_vec(self) result(x) 
        ! Dumps a zodi model to a parameter vector
        class(ZodiModel), intent(in) :: self
        real(dp) :: x(self%n_parameters)
        integer(i4b) :: i, j, k
        do i = 1, self%n_comps
            select type (a => self%comps(i)%c)
            class is (ZodiCloud)
                x(1) = a%n_0
                x(2) = a%incl
                x(3) = a%Omega
                x(4) = a%x_0
                x(5) = a%y_0
                x(6) = a%z_0
                x(7) = a%alpha
                x(8) = a%beta
                x(9) = a%gamma
                x(10) = a%mu
            class is (ZodiBand)
                select case (self%comp_labels(i))
                case ('band1')
                    j = 10
                case ('band2')
                    j = 20
                case ('band3')
                    j = 30
                case default
                    stop 'Invalid component label'
                end select                
                x(j + 1) = a%n_0
                x(j + 2) = a%incl
                x(j + 3) = a%Omega
                x(j + 4) = a%x_0
                x(j + 5) = a%y_0
                x(j + 6) = a%z_0
                x(j + 7) = a%delta_zeta
                x(j + 8) = a%delta_r
                x(j + 9) = a%v
                x(j + 10) = a%p
            class is (ZodiRing)
                x(41) = a%n_0
                x(42) = a%incl
                x(43) = a%Omega
                x(44) = a%x_0
                x(45) = a%y_0
                x(46) = a%z_0
                x(47) = a%R_0
                x(48) = a%sigma_r
                x(49) = a%sigma_z
            class is (ZodiFeature)
                x(50) = a%n_0
                x(51) = a%incl
                x(52) = a%Omega
                x(53) = a%x_0
                x(54) = a%y_0
                x(55) = a%z_0
                x(56) = a%R_0
                x(57) = a%sigma_r
                x(58) = a%sigma_z
                x(59) = a%theta_0
                x(60) = a%sigma_theta
            end select
        end do
        x(61) = self%T_0 
        x(62) = self%delta
        
        !emissivities
        k = 63
        do i = 1, self%n_comps
            do j = 1, self%n_bands
                x(k) = self%emissivity(i, j)
                k = k + 1
            end do
        end do

        !albedos
        k = 63 + self%n_comps * self%n_bands
        do i = 1, self%n_comps
            do j = 1, self%n_bands
                x(k) = self%albedo(i, j)
                k = k + 1
            end do
        end do
    end function model_to_param_vec

    subroutine output_to_hd5(self, cpar, iter)
        ! Writes the zodi model to an hdf file
        class(ZodiModel), intent(in) :: self
        type(comm_params), intent(in) :: cpar
        integer(i4b), intent(in) :: iter

        integer(i4b) :: i, j, hdferr, ierr, unit
        logical(lgt) :: exist, init, new_header
        character(len=6) :: itext
        character(len=4) :: ctext
        character(len=512) :: zodi_path, comp_path, param_path, chainfile, hdfpath
        character(len=10), allocatable :: param_names(:)
        real(dp), allocatable :: param_values(:)
        type(hdf_file) :: file
        type(h5o_info_t) :: object_info

        if (.not. cpar%myid_chain == 0) return

        call int2string(cpar%mychain, ctext)
        chainfile = trim(adjustl(cpar%outdir)) // '/chain' // &
            & '_c' // trim(adjustl(ctext)) // '.h5'

        inquire(file=trim(chainfile), exist=exist)
        call open_hdf_file(chainfile, file, 'b')

        call int2string(iter, itext)
        zodi_path = trim(adjustl(itext))//'/zodi'
        call create_hdf_group(file, trim(adjustl(zodi_path)))
        comp_path = trim(adjustl(zodi_path))//'/params/'
        
        call write_hdf(file, trim(adjustl(comp_path)), self%model_to_param_vec())

        call close_hdf_file(file)
    end subroutine

end module comm_zodi_mod