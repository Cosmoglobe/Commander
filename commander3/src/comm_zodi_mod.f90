!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================

module comm_zodi_mod    
    use comm_utils
    use comm_param_mod
    use comm_bp_mod
    use spline_1D_mod
    use comm_tod_mod
    
    implicit none

    private
    public :: initialize_zodi_mod, get_zodi_emission, update_zodi_splines

    ! TODO: probably stick this into an object, for instance comm_zodi_params
    ! which will make it easier to pass around for parameter estimation.
    ! Global parameters
    integer(i4b) :: gauss_degree, n_interp_points
    real(dp) :: EPS = 3.d-14
    real(dp) :: R_cutoff, delta_t_reset, min_ipd_temp, max_ipd_temp, SECOND_TO_DAY

    real(dp), allocatable :: tabulated_earth_time(:), tabulated_earth_pos(:, :), temperature_grid(:)
    logical(lgt) :: use_cloud, use_band1, use_band2, use_band3, use_ring, use_feature, &
                    apply_color_correction, use_unit_emissivity, use_unit_albedo

    ! K98 Model parameter variables
    real(dp) :: T_0, delta
    real(dp) :: cloud_x_0, cloud_y_0, cloud_z_0, cloud_incl, cloud_Omega, cloud_n_0, &
                cloud_alpha, cloud_beta, cloud_gamma, cloud_mu
    real(dp) :: band1_x_0, band1_y_0, band1_z_0, band1_incl, band1_Omega, band1_n_0, &
                band1_delta_zeta, band1_delta_R, band1_v, band1_p
    real(dp) :: band2_x_0, band2_y_0, band2_z_0, band2_incl, band2_Omega, band2_n_0, &
                band2_delta_zeta, band2_delta_R, band2_v, band2_p
    real(dp) :: band3_x_0, band3_y_0, band3_z_0, band3_incl, band3_Omega, band3_n_0, &
                band3_delta_zeta, band3_delta_R, band3_v, band3_p
    real(dp) :: ring_x_0, ring_y_0, ring_z_0, ring_incl, ring_Omega, ring_n_0, &
                ring_R, ring_sigma_R, ring_sigma_z
    real(dp) :: feature_x_0, feature_y_0, feature_z_0, feature_incl, feature_Omega, feature_n_0, &
                feature_R, feature_sigma_R, feature_sigma_z, feature_theta, feature_sigma_theta
    type(spline_type) :: earth_pos_spl_obj(3)

    type, abstract :: ZodiComponent
        ! Abstract base class for a zodiacal component.

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

    type, extends(ZodiComponent) :: Cloud
        real(dp) :: alpha, beta, gamma, mu
        contains
            procedure :: initialize => initialize_cloud
            procedure :: get_density => get_density_cloud
    end type Cloud

    type, extends(ZodiComponent) :: Band
        real(dp) :: delta_zeta, delta_r, v, p
        contains
            procedure :: initialize => initialize_band
            procedure :: get_density => get_density_band
    end type Band

    type, extends(ZodiComponent) :: Ring
        real(dp) :: R_0, sigma_r, sigma_z
        contains
            procedure :: initialize => initialize_ring
            procedure :: get_density => get_density_ring
    end type Ring

    type, extends(ZodiComponent) :: Feature
        real(dp) :: R_0, sigma_r, sigma_z, theta_0, sigma_theta
        contains
            procedure :: initialize => initialize_feature
            procedure :: get_density => get_density_feature
    end type Feature

    ! Initializing global ZodiComponent list and component instances
    class(ZodiComponent), pointer :: comp_list => null()
    type(Cloud), target :: cloud_comp
    type(Band), target :: band1_comp, band2_comp, band3_comp
    type(Ring), target :: ring_comp
    type(Feature), target :: feature_comp

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

        SECOND_TO_DAY = 1.d0 / (60*60*24)

        call initialize_zodi_parameters_from_cpar(cpar)
        call initialize_earth_pos_spline_object(cpar, earth_pos_spl_obj)

        ! Set up interpolation grid for evaluating line of sight b_nu
        allocate(temperature_grid(n_interp_points))
        call linspace(min_ipd_temp, max_ipd_temp, temperature_grid)
        
        ! Initialize zodi components
        if (use_cloud) then
            cloud_comp = Cloud(x_0=cloud_x_0, y_0=cloud_y_0, z_0=cloud_z_0, incl=cloud_incl, &
                               Omega=cloud_Omega, n_0=cloud_n_0, alpha=cloud_alpha, &
                               beta=cloud_beta, gamma=cloud_gamma, mu=cloud_mu)
            comp => cloud_comp
            call add_component_to_list(comp)
        end if
        if (use_band1) then
            band1_comp = Band(x_0=band1_x_0, y_0=band1_y_0, z_0=band1_z_0, incl=band1_incl, &
                              Omega=band1_Omega, n_0=band1_n_0, delta_zeta=band1_delta_zeta, &
                              delta_r=band1_delta_R, v=band1_v, p=band1_p)
            comp => band1_comp
            call add_component_to_list(comp)
        end if
        if (use_band2) then
            band2_comp = Band(x_0=band2_x_0, y_0=band2_y_0, z_0=band2_z_0, incl=band2_incl, &
                              Omega=band2_Omega, n_0=band2_n_0, delta_zeta=band2_delta_zeta, &
                              delta_r=band2_delta_R, v=band2_v, p=band2_p)
            comp => band2_comp
            call add_component_to_list(comp)
        end if
        if (use_band3) then
            band3_comp = Band(x_0=band3_x_0, y_0=band3_y_0, z_0=band3_z_0, incl=band3_incl, &
                              Omega=band3_Omega, n_0=band3_n_0, delta_zeta=band3_delta_zeta, &
                              delta_r=band3_delta_R, v=band3_v, p=band3_p)
            comp => band3_comp
            call add_component_to_list(comp)
        end if
        if (use_ring) then
            ring_comp = Ring(x_0=ring_x_0, y_0=ring_y_0, z_0=ring_z_0, incl=ring_incl, &
                             Omega=ring_Omega, n_0=ring_n_0, R_0=ring_R, sigma_r=ring_sigma_R, &
                             sigma_z=ring_sigma_z)
            comp => ring_comp
            call add_component_to_list(comp)
        end if
        if (use_feature) then
            feature_comp = Feature(x_0=feature_x_0, y_0=feature_y_0, z_0=feature_z_0, &
                                   incl=feature_incl, Omega=feature_Omega, n_0=feature_n_0, & 
                                   R_0=feature_R, sigma_r=feature_sigma_R, sigma_z=feature_sigma_z, &
                                   theta_0=feature_theta, sigma_theta=feature_sigma_theta)
            comp => feature_comp
            call add_component_to_list(comp)
        end if

        ! Add components to list for easy iteration
        comp => comp_list
        do while (associated(comp))
            call comp%initialize()
            comp => comp%next()
        end do

    end subroutine initialize_zodi_mod

    subroutine get_zodi_emission(tod, pix, scan, s_zodi)
        ! Returns the predicted zodiacal emission for a scan (chunk of time-ordered data).
        !
        ! Parameters
        ! ----------
        ! tod : class(comm_tod)
        !     The TOD object holding the spline objects to update.
        ! pix : integer(i4b), dimension(ntod, ndet, nhorns)
        !     The pixel indices of each time-ordered observation.
        ! scan : integer(i4b)
        !     The scan number.
        !
        ! Returns
        ! -------
        ! s_zodi : real(sp), dimension(ntod, ndet)
        !     The predicted zodiacal emission for each time-ordered observation.

        class(ZodiComponent), pointer :: comp

        class(comm_tod), intent(inout) :: tod
        integer(i4b), intent(in) :: pix(:, :), scan
        real(sp), intent(inout) :: s_zodi(:, :)

        logical(lgt) :: scattering
        integer(i4b) :: i, j, k, pixel, lookup_idx, n_det, n_tod
        real(dp) :: earth_lon, R_obs, R_max, dt_tod, obs_time
        real(dp) :: unit_vector(3), X_unit_LOS(3, gauss_degree), X_LOS(3, gauss_degree), obs_pos(3), earth_pos(3)
        real(dp), dimension(gauss_degree) :: R_LOS, T_LOS, density_LOS, gauss_nodes, gauss_weights, &
                                             comp_emission_LOS, solar_flux_LOS, scattering_angle, phase_function, b_nu_LOS

        if (.not. tod%zodi_tod_params_are_initialized) then 
            print *, "Error: Zodi parameters have not been initialized."
            print *, "make sure that `initialize_zodi_tod_parameters` is called in the constructor of `tod_your_experiment_mod`"
            stop
        else if (.not. allocated(tod%zodi_b_nu_spl_obj)) then
            print *, "Error: Zodi splines have not been initialized."
            print *, "make sure that `initialize_zodi_splines` is called at least once before this function is ran, and then again everytime the bandpasses are updated."
            stop
        end if

        s_zodi = 0.d0
        n_tod = size(pix, dim=1)
        n_det = size(pix, dim=2)

        dt_tod = tod%samprate * SECOND_TO_DAY ! dt between two samples in units of days (assumes equispaced tods)
        obs_pos = tod%scans(scan)%satpos
        R_obs = norm2(obs_pos)

        do i = 1, 3
            earth_pos(i) = splint_simple(earth_pos_spl_obj(i), tod%scans(scan)%t0(1))
        end do        
        earth_lon = atan(earth_pos(2), earth_pos(1))

        if (count(tod%zodi_spl_albedos /= 0.d0) > 0) then
            scattering = .true.
        else
            scattering = .false.
        end if

        do i = 1, n_tod
            ! Reset cache if time between last cache update and current time is larger than `delta_t_reset`.
            ! NOTE: this cache is only effective if the scans a core handles are in chronological order.
            obs_time = tod%scans(scan)%t0(1) + (i - 1) * dt_tod
            if (((obs_time - tod%zodi_cache_time) >= delta_t_reset) &
                .and. &
                (obs_time > tod%zodi_min_obs_time .and. obs_time < tod%zodi_max_obs_time) &
            ) then 
                do j = 1, 3 
                    earth_pos(j) = splint_simple(earth_pos_spl_obj(j), obs_time)
                    obs_pos(j) = splint_simple(tod%zodi_obs_pos_spl_obj(j), obs_time)
                end do            
                R_obs = norm2(obs_pos)
                earth_lon = atan(earth_pos(2), earth_pos(1))
                tod%zodi_cache = -1.d0
                tod%zodi_cache_time = obs_time
            end if

            do j = 1, n_det   
                lookup_idx = tod%pix2ind(pix(i, j))
                if (tod%zodi_cache(lookup_idx, j) >= 0.d0) then
                    s_zodi(i, j) = tod%zodi_cache(lookup_idx, j)
                    cycle
                end if

                unit_vector = tod%ind2vec_ecl(:, lookup_idx)
                call get_R_max(unit_vector, obs_pos, R_obs, R_max)
                call gauss_legendre_quadrature(x1=EPS, x2=R_max, n=gauss_degree, x=gauss_nodes, w=gauss_weights)

                do k = 1, 3
                    X_unit_LOS(k, :) = gauss_nodes * unit_vector(k)
                    X_LOS(k, :) = X_unit_LOS(k, :) + obs_pos(k)
                end do
                R_LOS = norm2(X_LOS, dim=1)           

                if (scattering) then
                    solar_flux_LOS = tod%zodi_spl_solar_irradiance(j) / R_LOS**2
                    call get_scattering_angle(X_LOS, X_unit_LOS, R_LOS, scattering_angle)
                    call get_phase_function(scattering_angle, tod%zodi_spl_phase_coeffs(j, :), tod%zodi_phase_func_normalization(j), phase_function)
                end if

                call get_dust_grain_temperature(R_LOS, T_LOS)
                call splint_simple_multi(tod%zodi_b_nu_spl_obj(j), T_LOS, b_nu_LOS)

                comp => comp_list
                k = 1
                do while (associated(comp))
                    call comp%get_density(X_LOS, earth_lon, density_LOS)
                    comp_emission_LOS = (1.d0 - tod%zodi_spl_albedos(j, k)) * (tod%zodi_spl_emissivities(j, k) * b_nu_LOS)
                    if (scattering) comp_emission_LOS = comp_emission_LOS + (tod%zodi_spl_albedos(j, k) * solar_flux_LOS * phase_function)
                    comp_emission_LOS = comp_emission_LOS * density_LOS

                    s_zodi(i, j) = s_zodi(i, j) + sum(comp_emission_LOS * gauss_weights)
                    comp => comp%next()
                    k = k + 1
                end do
                tod%zodi_cache(lookup_idx, j) = s_zodi(i, j)
            end do
        end do
    end subroutine get_zodi_emission

    subroutine update_zodi_splines(tod, bandpass, det)
        ! Updates the spectral spline objects in the TOD object.
        !
        ! In the K98 model, several spectral parameters are tabulated at individual frequencies, 
        ! which we need to evaluate over the bandpass. In a future version, we may want to fit 
        ! a modified blackbody which would allow us to drop using some of these spline objects.
        !
        ! Parameters
        ! ----------
        ! tod : class(comm_tod)
        !     The TOD object holding the spline objects to update.
        ! bandpass : class(comm_bp_ptr)
        !   The bandpass object holding the bandpass frequencies, and the SED2F function 
        !   (bandpass integration).
        ! det : integer(i4b)
        !   The detector to update the spline objects for.

        class(comm_tod),   intent(inout) :: tod
        class(comm_bp_ptr), intent(in) :: bandpass
        integer(i4b), intent(in) :: det

        real(dp), allocatable :: b_nu(:), integrals(:)
        integer(i4b) :: i, j

        allocate(integrals(n_interp_points))
        allocate(b_nu(bandpass%p%n))
        do i = 1, n_interp_points    
            call get_blackbody_emission(bandpass%p%nu, temperature_grid(i), b_nu)
            integrals(i) = bandpass%p%SED2F(b_nu)
        end do
        call spline_simple(tod%zodi_b_nu_spl_obj(det), temperature_grid, integrals, regular=.true.)

        do j = 1, size(tod%zodi_spl_emissivities, dim=2)
            tod%zodi_spl_emissivities(det, j) = splint_simple(tod%zodi_emissivity_spl_obj(j), bandpass%p%nu_c)
        end do

        do j = 1, size(tod%zodi_spl_albedos, dim=2)
            tod%zodi_spl_albedos(det, j) = splint_simple(tod%zodi_albedo_spl_obj(j), bandpass%p%nu_c)
        end do
        do j = 1, size(tod%zodi_spl_phase_coeffs, dim=2) 
            tod%zodi_spl_phase_coeffs(det, j) = splint_simple(tod%zodi_phase_coeff_spl_obj(j), bandpass%p%nu_c)
        end do

        tod%zodi_spl_solar_irradiance = splint_simple(tod%zodi_solar_irradiance_spl_obj, bandpass%p%nu_c)
        call get_phase_normalization(tod%zodi_spl_phase_coeffs(det, :), tod%zodi_phase_func_normalization(det))
    end subroutine update_zodi_splines

    ! Functions for evaluating the zodiacal emission
    ! -----------------------------------------------------------------------------------
    subroutine get_R_max(unit_vector, obs_pos, R_obs, R_max)
        ! Computes R_max (the length of the LOS such that it stops exactly at los_cutoff_radius).

        real(dp), intent(in), dimension(:) :: unit_vector, obs_pos
        real(dp), intent(in) :: R_obs
        real(dp), intent(out) :: R_max
        real(dp) :: lon, lat, cos_lat, b, d, q

        lon = atan(unit_vector(2), unit_vector(1))
        lat = asin(unit_vector(3))
        cos_lat = cos(lat)
        b = 2.d0 * (obs_pos(1) * cos_lat * cos(lon) + obs_pos(2) * cos_lat * sin(lon))
        d = R_obs**2 - R_cutoff**2
        q = -0.5d0 * b * (1.d0 + sqrt(b**2 - (4.d0 * d)) / abs(b))
        R_max = max(q, d / q)
    end subroutine get_R_max

    subroutine get_dust_grain_temperature(R, T_out)
        real(dp), dimension(:), intent(in) :: R
        real(dp), dimension(:), intent(out) :: T_out
        T_out = T_0 * R ** (-delta)
    end subroutine get_dust_grain_temperature

    subroutine get_blackbody_emission(nus, T, b_nu)
        real(dp), intent(in) :: nus(:), T
        real(dp), dimension(:), intent(out) :: b_nu
        b_nu = 1d20 * ((2 * h * nus**3) / (c*c)) / (exp((h * nus) / (k_B * T)) - 1.d0)
    end subroutine get_blackbody_emission

    subroutine get_scattering_angle(X_helio_vec_LOS, X_vec_LOS, R_helio_LOS, scattering_angle)
        real(dp), intent(in) :: X_helio_vec_LOS(:, :), X_vec_LOS(:, :), R_helio_LOS(:)
        real(dp), dimension(:), intent(out) :: scattering_angle
        real(dp), dimension(gauss_degree) :: cos_theta, R_LOS
        
        R_LOS = norm2(X_vec_LOS, dim=1)
        cos_theta = sum(X_helio_vec_LOS * X_vec_LOS, dim=1) / (R_LOS * R_helio_LOS)
        ! clip cos(theta) to [-1, 1]
        where (cos_theta > 1)
            cos_theta = 1
        elsewhere (cos_theta < -1)
            cos_theta = -1
        endwhere         

        scattering_angle = acos(-cos_theta)
    end subroutine get_scattering_angle

    subroutine get_phase_function(scattering_angle, phase_coefficients, normalization_factor, phase_function)
        real(dp), intent(in) :: scattering_angle(:), phase_coefficients(:), normalization_factor
        real(dp), intent(out) :: phase_function(:)

        phase_function = normalization_factor * (phase_coefficients(1) + phase_coefficients(2) &
                         * scattering_angle + exp(phase_coefficients(3) * scattering_angle))
    end subroutine

    subroutine get_phase_normalization(phase_coefficients, normalization_factor)
        real(dp), intent(in) :: phase_coefficients(:)
        real(dp), intent(out) :: normalization_factor
        real(dp) :: term1, term2, term3, term4

        term1 = 2.d0 * pi
        term2 = 2.d0 * phase_coefficients(1)
        term3 = pi * phase_coefficients(2)
        term4 = (exp(phase_coefficients(3) * pi) + 1.d0) / (phase_coefficients(3)**2 + 1.d0)
        normalization_factor = 1.d0 / (term1 * (term2 + term3 + term4))
    end subroutine

    ! Methods for initizializing the zodiacal components
    ! -----------------------------------------------------------------------------------
    subroutine initialize_cloud(self)
        class(Cloud) :: self
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_cloud

    subroutine initialize_band(self)
        class(Band) :: self
        self%delta_zeta = self%delta_zeta * deg2rad
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_band

    subroutine initialize_ring(self)
        class(Ring) :: self
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_ring

    subroutine initialize_feature(self)
        class(Feature) :: self
        self%theta_0 = self%theta_0 * deg2rad
        self%sigma_theta = self%sigma_theta * deg2rad
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_feature


    ! Methods describing the densitry distribution of the zodiacal components
    ! -----------------------------------------------------------------------------------
    subroutine get_density_cloud(self, X_vec, theta, n_out)
        class(Cloud) :: self
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
        class(Band) :: self
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

            zeta_over_delta_zeta = zeta / self%delta_zeta
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
        class(Ring) :: self
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
        class(Feature) :: self
        real(dp), dimension(:, :), intent(in) :: X_vec
        real(dp), intent(in) :: theta
        real(dp), dimension(:), intent(out) :: n_out
        integer(i4b) :: i
        real(dp) :: x_prime, y_prime, z_prime, R, Z_midplane, theta_prime, exp_term

        do i = 1, gauss_degree
            x_prime = X_vec(1, i) - self%x_0
            y_prime = X_vec(2, i) - self%y_0
            z_prime = X_vec(3, i) - self%z_0
            theta_prime = atan2(y_prime, x_prime) - (theta + self%theta_0)

            ! Constraining the angle to the limit [-pi, pi]
            do while (theta_prime < -pi)
                theta_prime = theta_prime + 2.d0*pi
            end do
            do while (theta_prime > pi)
                theta_prime = theta_prime - 2.d0*pi
            end do

            R = sqrt(x_prime*x_prime + y_prime*y_prime + z_prime*z_prime)
            Z_midplane = (x_prime*self%sin_omega - y_prime*self%cos_omega)*self%sin_incl + z_prime*self%cos_incl

            exp_term = ((R - self%R_0) ** 2 / self%sigma_r**2) + (abs(Z_midplane) / self%sigma_z) + (theta_prime**2 / self%sigma_theta**2)

            n_out(i) = self%n_0 * exp(-exp_term)
        end do
    end subroutine get_density_feature


    ! Functions for constructing and iterating the linked list of zodiacal components
    ! -----------------------------------------------------------------------------------
    function next(self)
        class(ZodiComponent) :: self
        class(ZodiComponent), pointer :: next
        next => self%next_link
    end function next

    subroutine add(self,link)
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

    ! Functions for initializing zodi_mod
    ! -----------------------------------------------------------------------------------
    subroutine initialize_earth_pos_spline_object(cpar, earth_spl_obj)
        ! Returns the spline object which is used to evaluate the earth position

        type(comm_params), intent(in) :: cpar
        type(spline_type), intent(inout) :: earth_spl_obj(3)
        integer :: i, n_earthpos, unit

        unit = getlun()
        open(unit, file=trim(trim(cpar%datadir)//'/'//trim("earth_pos_1980-2050_ephem_de432s.txt")))
        read(unit, *) n_earthpos
        read(unit, *) ! skip header

        allocate(tabulated_earth_pos(3, n_earthpos))
        allocate(tabulated_earth_time(n_earthpos))
        do i = 1, n_earthpos
        read(unit,*) tabulated_earth_time(i), tabulated_earth_pos(1, i), tabulated_earth_pos(2, i), tabulated_earth_pos(3, i)
        end do
        close(unit)
        do i = 1, 3
            call spline_simple(earth_spl_obj(i), tabulated_earth_time, tabulated_earth_pos(i, :), regular=.true.)
        end do
    end subroutine initialize_earth_pos_spline_object

    subroutine initialize_zodi_parameters_from_cpar(cpar)
        ! Initialize hyper parameters for commander run
        type(comm_params), intent(in) :: cpar

        ! Hypper parameters
        use_cloud = cpar%zs_use_cloud
        use_band1 = cpar%zs_use_band1
        use_band2 = cpar%zs_use_band2
        use_band3 = cpar%zs_use_band3
        use_ring = cpar%zs_use_ring
        use_feature = cpar%zs_use_feature
        use_unit_emissivity = cpar%zs_use_unit_emissivity
        R_cutoff = cpar%zs_los_cut
        gauss_degree = cpar%zs_gauss_quad_order
        delta_t_reset = cpar%zs_delta_t_reset
        min_ipd_temp = cpar%zs_min_ipd_temp
        max_ipd_temp = cpar%zs_max_ipd_temp
        n_interp_points = cpar%zs_n_interp_points

        ! Source parameters
        T_0 = cpar%zs_t_0
        delta = cpar%zs_delta

        ! Shape parameters
        cloud_x_0 = cpar%zs_common(1, 1)
        cloud_y_0 = cpar%zs_common(1, 2)
        cloud_z_0 = cpar%zs_common(1, 3)
        cloud_incl = cpar%zs_common(1, 4)
        cloud_Omega = cpar%zs_common(1, 5)
        cloud_n_0 = cpar%zs_common(1, 6)
        cloud_alpha = cpar%zs_cloud_alpha
        cloud_beta = cpar%zs_cloud_beta
        cloud_gamma = cpar%zs_cloud_gamma
        cloud_mu = cpar%zs_cloud_mu

        band1_x_0 = cpar%zs_common(2, 1)
        band1_y_0 = cpar%zs_common(2, 2)
        band1_z_0 = cpar%zs_common(2, 3)
        band1_incl = cpar%zs_common(2, 4)
        band1_Omega = cpar%zs_common(2, 5)
        band1_n_0 = cpar%zs_common(2, 6)
        band1_delta_zeta = cpar%zs_bands_delta_zeta(1)
        band1_delta_R = cpar%zs_bands_delta_r(1)
        band1_v = cpar%zs_bands_v(1)
        band1_p = cpar%zs_bands_p(1)

        band2_x_0 = cpar%zs_common(3, 1)
        band2_y_0 = cpar%zs_common(3, 2)
        band2_z_0 = cpar%zs_common(3, 3)
        band2_incl = cpar%zs_common(3, 4)
        band2_Omega = cpar%zs_common(3, 5)
        band2_n_0 = cpar%zs_common(3, 6)
        band2_delta_zeta = cpar%zs_bands_delta_zeta(2)
        band2_delta_R = cpar%zs_bands_delta_r(2)
        band2_v = cpar%zs_bands_v(2)
        band2_p = cpar%zs_bands_p(2)

        band3_x_0 = cpar%zs_common(4, 1)
        band3_y_0 = cpar%zs_common(4, 2)
        band3_z_0 = cpar%zs_common(4, 3)
        band3_incl = cpar%zs_common(4, 4)
        band3_Omega = cpar%zs_common(4, 5)
        band3_n_0 = cpar%zs_common(4, 6)
        band3_delta_zeta = cpar%zs_bands_delta_zeta(3)
        band3_delta_R = cpar%zs_bands_delta_r(3)
        band3_v = cpar%zs_bands_v(3)
        band3_p = cpar%zs_bands_p(3)

        ring_x_0 = cpar%zs_common(5, 1)
        ring_y_0 = cpar%zs_common(5, 2)
        ring_z_0 = cpar%zs_common(5, 3)
        ring_incl = cpar%zs_common(5, 4)
        ring_Omega = cpar%zs_common(5, 5)
        ring_n_0 = cpar%zs_common(5, 6)
        ring_R = cpar%zs_ring_r
        ring_sigma_R = cpar%zs_ring_sigma_r
        ring_sigma_z = cpar%zs_ring_sigma_z

        feature_x_0 = cpar%zs_common(6, 1)
        feature_y_0 = cpar%zs_common(6, 2)
        feature_z_0 = cpar%zs_common(6, 3)
        feature_incl = cpar%zs_common(6, 4)
        feature_Omega = cpar%zs_common(6, 5)
        feature_n_0 = cpar%zs_common(6, 6)
        feature_R = cpar%zs_feature_r
        feature_sigma_R = cpar%zs_feature_sigma_r
        feature_sigma_z = cpar%zs_feature_sigma_z
        feature_theta = cpar%zs_feature_theta
        feature_sigma_theta = cpar%zs_feature_sigma_theta
    end subroutine initialize_zodi_parameters_from_cpar
end module comm_zodi_mod