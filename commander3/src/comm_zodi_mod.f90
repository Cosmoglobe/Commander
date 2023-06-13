module comm_zodi_mod
    use comm_utils
    use comm_param_mod
    use spline_1D_mod

    implicit none

    private
    public initialize_zodi_mod, get_s_zodi, zodi_params_k98, ZodiComponent, comp_list, zodi

    ! Global variables
    integer(i4b) :: gauss_degree, n_interp_points
    real(dp) :: R_cutoff
    logical(lgt) :: use_cloud, use_band1, use_band2, use_band3, use_ring, use_feature

    type, abstract :: ZodiComponent
        ! Abstract base class for a zodiacal component which also acts as a 
        ! linked list container.

        ! Pointers to the next/prev links in the linked list
        class(ZodiComponent), pointer :: next_link => null()
        class(ZodiComponent), pointer :: prev_link => null()

        ! Shared component variables
        real(dp) :: x_0, y_0, z_0, incl, Omega, n_0
        real(dp), allocatable :: sin_omega, cos_omega, sin_incl, cos_incl

    contains
        ! Shared component procedures
        procedure(initialize_interface), deferred :: initialize
        procedure(density_interface), deferred :: get_density

        ! Linked list procedures
        procedure :: next
        procedure :: add
    end type ZodiComponent

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

    ! Linked list of zodi components for iteration
    class(ZodiComponent), pointer :: comp_list => null()

    ! Stores the interplanetary dust model parameters
    type :: zodi_params_k98
        integer(i4b) :: n_comps, n_bands
        real(dp) :: T_0, delta ! solar system temperature parameters
        real(dp), allocatable, dimension(:) :: nu_ref, solar_irradiance ! spectral parameters
        real(dp), allocatable, dimension(:, :) :: phase_coeffs ! spectral parameters
        type(ZodiCloud) :: cloud ! Zodiacal components
        type(ZodiBand) :: band1, band2, band3
        type(ZodiRing) :: ring
        type(ZodiFeature) :: feature
        type(spline_type) :: solar_irradiance_spl! spline interpolators
        type(spline_type), allocatable, dimension(:) :: phase_coeff_spl
    contains
        procedure :: initialize_k98_model, build_splines
    end type zodi_params_k98

    ! Global zodi parameter object
    type(zodi_params_k98), target :: zodi

contains
    subroutine initialize_zodi_mod(cpar)
        ! Initialize the zodi module.
        !
        ! Parameters
        ! ----------
        ! cpar: comm_params
        !    Parameter file variables.

        type(comm_params), intent(in) :: cpar
        class(ZodiComponent), pointer :: comp

        call initialize_hyper_parameters(cpar)
        call zodi%initialize_k98_model(cpar)
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


    ! Methods describing the densitry distribution of the zodiacal components
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

    subroutine initialize_k98_model(self, cpar)
        ! Initialize the kelsal et al 1998 model.
        !
        ! Parameters
        ! ----------
        ! self: zodi_params_k98
        !    Zodi model object.
        ! cpar: comm_params
        !    Parameter file variables.

        class(zodi_params_k98), target, intent(inout) :: self
        type(comm_params), intent(in) :: cpar
        class(ZodiComponent), pointer :: comp
        integer(i4b) :: i

        ! aux
        self%n_bands = cpar%zs_nbands
        self%n_comps = cpar%zs_ncomps

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

        allocate(self%phase_coeff_spl(3))
        call self%build_splines()

        ! Initialize class instances
        if (use_cloud) then
            self%cloud = ZodiCloud(x_0=cpar%zs_common(1, 1), y_0=cpar%zs_common(1, 2), z_0=cpar%zs_common(1, 3), &
                                    incl=cpar%zs_common(1, 4), Omega=cpar%zs_common(1, 5), n_0=cpar%zs_common(1, 6), & 
                                    alpha=cpar%zs_cloud_alpha, beta=cpar%zs_cloud_beta, gamma=cpar%zs_cloud_gamma, &
                                    mu=cpar%zs_cloud_mu)
            comp => self%cloud
            call add_component_to_list(comp)
        end if
        if (use_band1) then
            self%band1 = ZodiBand(x_0=cpar%zs_common(2, 1), y_0=cpar%zs_common(2, 2), z_0=cpar%zs_common(2, 3), &
                              incl=cpar%zs_common(2, 4), Omega=cpar%zs_common(2, 5), n_0=cpar%zs_common(2, 6), & 
                              delta_zeta=cpar%zs_bands_delta_zeta(1), delta_r=cpar%zs_bands_delta_r(1), & 
                              v=cpar%zs_bands_v(1), p=cpar%zs_bands_p(1))
            comp => self%band1
            call add_component_to_list(comp)
        end if
        if (use_band2) then
            self%band2 = ZodiBand(x_0=cpar%zs_common(3, 1), y_0=cpar%zs_common(3, 2), z_0=cpar%zs_common(3, 3), &
                              incl=cpar%zs_common(3, 4), Omega=cpar%zs_common(3, 5), n_0=cpar%zs_common(3, 6), & 
                              delta_zeta=cpar%zs_bands_delta_zeta(2), delta_r=cpar%zs_bands_delta_r(2), & 
                              v=cpar%zs_bands_v(2), p=cpar%zs_bands_p(2))
            comp => self%band2
            call add_component_to_list(comp)
        end if
        if (use_band3) then
            self%band3 = ZodiBand(x_0=cpar%zs_common(4, 1), y_0=cpar%zs_common(4, 2), z_0=cpar%zs_common(4, 3), &
                              incl=cpar%zs_common(4, 4), Omega=cpar%zs_common(4, 5), n_0=cpar%zs_common(4, 6), & 
                              delta_zeta=cpar%zs_bands_delta_zeta(3), delta_r=cpar%zs_bands_delta_r(3), & 
                              v=cpar%zs_bands_v(3), p=cpar%zs_bands_p(3))
            comp => self%band3
            call add_component_to_list(comp)
        end if
        if (use_ring) then
            self%ring = ZodiRing(x_0=cpar%zs_common(5, 1), y_0=cpar%zs_common(5, 2), z_0=cpar%zs_common(5, 3), &
                             incl=cpar%zs_common(5, 4), Omega=cpar%zs_common(5, 5), n_0=cpar%zs_common(5, 6), &
                             R_0=cpar%zs_ring_r, sigma_r=cpar%zs_ring_sigma_r, sigma_z=cpar%zs_ring_sigma_z)
            comp => self%ring
            call add_component_to_list(comp)
        end if
        if (use_feature) then
            self%feature = ZodiFeature(x_0=cpar%zs_common(6, 1), y_0=cpar%zs_common(6, 2), z_0=cpar%zs_common(6, 3), &
                                   incl=cpar%zs_common(6, 4), Omega=cpar%zs_common(6, 5), n_0=cpar%zs_common(6, 6), &
                                   R_0=cpar%zs_feature_r, sigma_r=cpar%zs_feature_sigma_r, sigma_z=cpar%zs_feature_sigma_z, &
                                   theta_0=cpar%zs_feature_theta, sigma_theta=cpar%zs_feature_sigma_theta)
            comp => self%feature
            call add_component_to_list(comp)
        end if

        ! Call initialize method on each component
        comp => comp_list
        do while (associated(comp))
            call comp%initialize()
            comp => comp%next()
        end do        
    end subroutine initialize_k98_model

    subroutine build_splines(self)
        ! Build splines for the phase function, and solar irradiance
        ! for each component in the model
        class(zodi_params_k98), intent(inout) :: self
        integer(i4b) :: i

        do i = 1, 3
            call spline_simple(self%phase_coeff_spl(i), self%nu_ref, self%phase_coeffs(:, i), regular=.false.)
        end do
        call spline_simple(self%solar_irradiance_spl, self%nu_ref, self%solar_irradiance, regular=.false.)
    end subroutine build_splines

    ! Functions for constructing and iterating the linked list of zodiacal components
    ! -----------------------------------------------------------------------------------
    function next(self)
        class(ZodiComponent) :: self
        class(ZodiComponent), pointer :: next
        next => self%next_link
    end function next

    subroutine add(self, link)
        class(ZodiComponent), target  :: self
        class(ZodiComponent), pointer :: link
        class(ZodiComponent), pointer :: comp

        comp => self
        do while (associated(comp%next_link))
            comp => comp%next_link
        end do
        link%prev_link => comp
        comp%next_link => link
    end subroutine add

    subroutine add_component_to_list(comp)
        class(ZodiComponent), pointer :: comp

        if (.not. associated(comp_list)) then
            comp_list => comp
        else
            call comp_list%add(comp)
        end if
    end subroutine add_component_to_list

end module comm_zodi_mod