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
    !   The zodi module handles simulating and fitting zodiacal emission at tod level.
    !
    !   Methods
    !   -------
    !   initialize_zodi_mod
    !       Initializes the zodi_mod.

    !   get_zodi_emission
    !       Method that simulates the zodiacal emission observed by the instrument for
    !       a given pixel.
    use comm_utils
    use comm_param_mod
    use comm_bp_mod
    use spline_1D_mod
    
    implicit none

    private
    public :: initialize_zodi_mod, get_zodi_emission, update_zodi_spline_obj, init_zodi_spline_objects

    ! Global parameters
    integer(i4b) :: gauss_degree
    integer(i4b) :: n_interpolation_points=100
    real(dp) :: min_ipd_temperature=40, max_ipd_temperature=550

    real(dp) :: EPS = 3.d-14
    real(dp) :: R_cutoff
    real(dp), allocatable :: unique_nsides(:), tabulated_earth_time(:), temperature_grid(:)
    real(dp), allocatable :: tabulated_earth_pos(:, :)
    type(spline_type) :: solar_irradiance_spline_obj
    type(spline_type) :: spline_earth_pos_obj(3), spline_obs_pos_obj(3)
    type(spline_type), allocatable :: emissivity_spline_obj(:), albedo_spline_obj(:), phase_coeff_spline_obj(:), b_nu_spline_obj(:)

    logical(lgt) :: use_cloud, use_band1, use_band2, use_band3, use_ring, use_feature, &
                    apply_color_correction, use_unit_emissivity, use_unit_albedo, scattering

    ! Model parameters
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
    real(dp), allocatable :: emissivities(:, :), albedos(:, :), phase_coeffs(:, :)
    real(dp), allocatable :: splined_emissivities(:, :), splined_albedos(:, :), splined_phase_coeffs(:, :), solar_irradiances(:), splined_solar_irradiance(:), phase_function_normalization(:)



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
    type(Feature), target  :: feature_comp

    ! Initializing container for precomputed pixel to unit vectors for each unique 
    ! nside in the Commander run.
    type UnitVector
        real(dp), allocatable :: elements(:, :)
    end type UnitVector

    type UnitVectorList
        type(UnitVector), allocatable :: vectors(:)
    end type UnitVectorList
    type(UnitVectorList) :: unit_vector_list


contains
    subroutine initialize_zodi_mod(cpar)
        !   Initialize `comm_zodi_mod` specific quantities and zodiacal components.
        !
        !   Parameters
        !   ----------
        !   cpar: comm_params
        !      Parameter file variables.
        implicit none

        type(comm_params), intent(in) :: cpar
        class(ZodiComponent), pointer :: comp

        integer(i4b) :: i, j, npix, nside, unit, n_earthpos
        character(len=1024) :: tabulated_earth_pos_filename
        real(dp) :: vec(3), ecliptic_to_galactic_matrix(3, 3)
        integer(i4b), allocatable :: sorted_unique_nsides(:)
        real(dp), allocatable :: galactic_vec(:, :)

        EPS = 3.d-14

        ! Initialize hyper, shape, and source parameters for ipd model from files
        call init_hyper_parameters(cpar)
        call init_source_parameters(cpar)
        call init_shape_parameters(cpar)

        ! Initialize Zodi components
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

        comp => comp_list
        do while (associated(comp))
            call comp%initialize()
            comp => comp%next()
        end do

        ! Precompute unit vector coordinates for galactic to ecliptic coordinates
        call gal_to_ecl_conversion_matrix(ecliptic_to_galactic_matrix)

        sorted_unique_nsides = unique_sort(pack(cpar%ds_nside, cpar%ds_nside /= 0))
        allocate(unit_vector_list%vectors(size(sorted_unique_nsides)))
        do i = 1, size(sorted_unique_nsides)
            nside = sorted_unique_nsides(i)
            npix = nside2npix(nside)
            allocate(unit_vector_list%vectors(i)%elements(0:npix-1, 3))
            allocate(galactic_vec(0:npix-1, 3))
        
            do j = 0, npix - 1
                call pix2vec_ring(nside, j, vec)
                galactic_vec(j, 1) = vec(1)
                galactic_vec(j, 2) = vec(2)
                galactic_vec(j, 3) = vec(3)
            end do
        
            unit_vector_list%vectors(i)%elements = matmul(galactic_vec, ecliptic_to_galactic_matrix)
            deallocate(galactic_vec)
        end do

        ! Read in tabulated earth position
        unit = getlun()
        tabulated_earth_pos_filename = trim(cpar%datadir)//'/'//trim("earth_pos_1980-2050_ephem_de432s.txt")
        open(unit, file=trim(tabulated_earth_pos_filename))
        read(unit, *) n_earthpos
        read(unit, *) ! skip header
        allocate(tabulated_earth_pos(3, n_earthpos))
        allocate(tabulated_earth_time(n_earthpos))
        do i = 1, n_earthpos
        read(unit,*) tabulated_earth_time(i), tabulated_earth_pos(1, i), tabulated_earth_pos(2, i), tabulated_earth_pos(3, i)
        end do
        close(unit)

        ! Set up spline objects
        allocate(emissivity_spline_obj(cpar%zs_ncomps))
        allocate(albedo_spline_obj(cpar%zs_ncomps))
        allocate(phase_coeff_spline_obj(3))

        ! Earths position
        do i = 1, 3
            call spline_simple(spline_earth_pos_obj(i), tabulated_earth_time, tabulated_earth_pos(i, :), regular=.true.)
        end do


        do i = 1, cpar%zs_ncomps
            call spline_simple(emissivity_spline_obj(i), cpar%zs_nu_ref, emissivities(:, i), regular=.false.)
            call spline_simple(albedo_spline_obj(i), cpar%zs_nu_ref, albedos(:, i), regular=.false.)
        end do
        do i = 1, 3
            call spline_simple(phase_coeff_spline_obj(i), cpar%zs_nu_ref, phase_coeffs(:, i), regular=.false.)
        end do
        call spline_simple(solar_irradiance_spline_obj, cpar%zs_nu_ref, solar_irradiances, regular=.false.)
        
        allocate(temperature_grid(n_interpolation_points))

        call linspace(min_ipd_temperature, max_ipd_temperature, temperature_grid)

    end subroutine initialize_zodi_mod

    subroutine init_zodi_spline_objects(obs_time, obs_pos)
        ! Initializes various spline objects for the zodi mod. This function is required to be called 
        ! in the tod_x_mod for an experiment .
        implicit none
        real(dp), intent(in) :: obs_time(:), obs_pos(:, :)
        integer(i4b) :: i
        do i = 1, 3
            call spline_simple(spline_obs_pos_obj(i), obs_time, obs_pos(i, :))
        end do
    end subroutine init_zodi_spline_objects

    subroutine get_zodi_emission(nside, pix, obs_time_at_chunk_start, satpos, s_zodi)
        !   Compute simulated zodiacal emission.
        !
        !   Given a set of observations, the observer position and time of 
        !   observation is used to perform line-of-sight integrations through 
        !   the interplanetary dust distribution in the solar system.
        !
        !   Args:
        !       nside (int): HEALPix resolution parameter for pixels in `pix`.
        !       pix (1D array): Pixel indices representing observations in a chunk of time-ordered 
        !           data. Dimensions are (`n_tod`, `n_det`).
        !       obs_time_at_chunk_start (float): Time of the first observation in the chunk of
        !
        !   Returns:
        !       s_zodi (1D array): Simulated zodiacal emission over the chunk of time-ordered data
        !           pixel pointing. Dimensions are (`n_tod`, `n_det`).


        implicit none
        class(ZodiComponent), pointer :: comp

        integer(i4b), intent(in) :: nside, pix(1:, 1:)
        real(dp), intent(in) :: obs_time_at_chunk_start, satpos(3)
        real(sp), intent(out) :: s_zodi(1:, 1:)
        
        integer(i4b) :: i, j, k, l, pix_idx, n_detectors, n_tods, npix
        real(dp) :: earth_lon, R_obs, R_max, samp_rate, dt_tod, time_tod, delta_t_reset, SECOND_TO_DAY, current_time
        real(dp) :: unit_vector(3), X_vec_LOS(3, gauss_degree), X_helio_vec_LOS(3, gauss_degree), tod_obs_pos(3), tod_earth_pos(3)
        real(dp), allocatable :: tabulated_unit_vectors(:, :), cached_s_zodi(:)

        ! Line of sight arrays and quantities
        real(dp), dimension(gauss_degree) :: R_helio_LOS, R_LOS, T_LOS, density_LOS, gauss_nodes, gauss_weights, &
                                             comp_emission_LOS, solar_flux_LOS, scattering_angle, phase_function, b_nu_LOS
        n_tods = size(pix, dim=1)
        n_detectors = size(pix, dim=2)
        npix = nside2npix(nside)
        s_zodi = 0.d0

        ! HYPERPARAMETERS
        samp_rate = 1.d0 / 8
        SECOND_TO_DAY = 1.d0 / (60*60*24)
        dt_tod = samp_rate * second_to_day ! dt between two samples in units of days
        current_time = obs_time_at_chunk_start

        allocate(cached_s_zodi(0:npix-1))
        cached_s_zodi = 0.d0

        ! Get unit vectors corresponding to the observed pixels (only on first chunk)
        allocate(tabulated_unit_vectors(0:npix-1, 3))
        call get_tabulated_unit_vectors(npix, tabulated_unit_vectors)

        if (.not. allocated(b_nu_spline_obj)) stop "b_nu_spline_obj not allocated"

        do l = 1, 3 ! for x, y, z
            tod_earth_pos(l) = splint_simple(spline_earth_pos_obj(l), current_time)
            tod_obs_pos(l) = splint_simple(spline_obs_pos_obj(l), current_time)
        end do        
        
        delta_t_reset = 0.1 ! time before reseting cache
        tod_obs_pos = satpos
        R_obs = norm2(tod_obs_pos)
        earth_lon = atan(tod_earth_pos(2), tod_earth_pos(1))

        ! Loop over each detectors time-ordered pointing chunks. For each unique pixel
        ! observed (cross dectors) perform a line of sight integral solving Eq. (20) in 
        ! San et al. 2022 (https://www.aanda.org/articles/aa/pdf/forth/aa44133-22.pdf).
        do i = 1, n_tods
            ! After a time, delta_t_reset, update the earth and observer position and reset cached s_zodi.
            time_tod = obs_time_at_chunk_start + (i - 1) * dt_tod
            if ((time_tod - current_time) >= delta_t_reset) then
                do l = 1, 3 
                    tod_earth_pos(l) = splint_simple(spline_earth_pos_obj(l), time_tod)
                    tod_obs_pos(l) = splint_simple(spline_obs_pos_obj(l), time_tod)
                end do            
                R_obs = norm2(tod_obs_pos)
                earth_lon = atan(tod_earth_pos(2), tod_earth_pos(1))
                cached_s_zodi = 0.d0
                current_time = time_tod
            end if

            do j = 1, n_detectors    
                pix_idx = pix(i, j)
                if (cached_s_zodi(pix_idx) /= 0.d0) then
                    s_zodi(i, j) = cached_s_zodi(pix_idx) ! Use cached value from previous chunk ordetector
                    cycle
                end if

                unit_vector = tabulated_unit_vectors(pix_idx, :)
                call get_R_max(unit_vector, tod_obs_pos, R_obs, R_max)
                call gauss_legendre_quadrature(x1=EPS, x2=R_max, n=gauss_degree, x=gauss_nodes, w=gauss_weights)

                do k = 1, 3 ! for x, y, z
                    X_vec_LOS(k, :) = gauss_nodes * unit_vector(k)
                    X_helio_vec_LOS(k, :) = X_vec_LOS(k, :) + tod_obs_pos(k)
                end do
                R_helio_LOS = norm2(X_helio_vec_LOS, dim=1)           

                if (scattering) then
                    solar_flux_LOS = splined_solar_irradiance(j) / R_helio_LOS**2
                    R_LOS = norm2(X_vec_LOS, dim=1)
                    call get_scattering_angle(X_helio_vec_LOS, X_vec_LOS, R_helio_LOS, R_LOS, scattering_angle)
                    call get_phase_function(scattering_angle, splined_phase_coeffs(j, :), phase_function_normalization(j), phase_function)
                end if

                call get_dust_grain_temperature(R_helio_LOS, T_LOS)
                call splint_simple_multi(b_nu_spline_obj(j), T_LOS, b_nu_LOS)

                comp => comp_list
                k = 1
                do while (associated(comp))
                    call comp%get_density(X_helio_vec_LOS, earth_lon, density_LOS)
                    comp_emission_LOS = (1.d0 - splined_albedos(j, k)) * (splined_emissivities(j, k) * b_nu_LOS)
                    if (scattering) comp_emission_LOS = comp_emission_LOS + (splined_albedos(j, k) * solar_flux_LOS * phase_function)
                    comp_emission_LOS = comp_emission_LOS * density_LOS

                    s_zodi(i, j) = s_zodi(i, j) + sum(comp_emission_LOS * gauss_weights)
                    comp => comp%next()
                    k = k + 1
                end do
                cached_s_zodi(pix_idx) = s_zodi(i, j) ! Update cache with newly computed LOS emission
            end do
        end do
    end subroutine get_zodi_emission


    ! Functions used in `get_zodi_emission`
    ! -------------------------------------
    subroutine get_R_max(unit_vector, obs_pos, R_obs, R_max)
        ! Computes R_max (the length of the LOS such that it stops exactly at los_cutoff_radius).
        implicit none
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
        implicit none
        real(dp), dimension(:), intent(in) :: R
        real(dp), dimension(:), intent(out) :: T_out
        T_out = T_0 * R ** (-delta)
    end subroutine get_dust_grain_temperature

    subroutine get_blackbody_emission(nus, T, b_nu)
        implicit none
        real(dp), intent(in) :: nus(:), T
        real(dp), dimension(:), intent(out) :: b_nu
        b_nu = 1d20 * ((2 * h * nus**3) / (c*c)) / (exp((h * nus) / (k_B * T)) - 1.d0)
    end subroutine get_blackbody_emission

    subroutine get_scattering_angle(X_helio_vec_LOS, X_vec_LOS, R_helio_LOS, R_LOS, scattering_angle)
        implicit none
        real(dp), intent(in) :: X_helio_vec_LOS(:, :), X_vec_LOS(:, :), R_helio_LOS(:), R_LOS(:)
        real(dp), dimension(:), intent(out) :: scattering_angle
        real(dp), dimension(gauss_degree) :: cos_theta

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
        implicit none
        real(dp), intent(in) :: scattering_angle(:), phase_coefficients(:), normalization_factor
        real(dp), intent(out) :: phase_function(:)

        phase_function = normalization_factor * (phase_coefficients(1) + phase_coefficients(2) &
                         * scattering_angle + exp(phase_coefficients(3) * scattering_angle))
    end subroutine

    subroutine get_phase_normalization(phase_coefficients, normalization_factor)
        implicit none
        real(dp), intent(in) :: phase_coefficients(:)
        real(dp), intent(out) :: normalization_factor
        real(dp) :: term1, term2, term3, term4

        term1 = 2.d0 * pi
        term2 = 2.d0 * phase_coefficients(1)
        term3 = pi * phase_coefficients(2)
        term4 = (exp(phase_coefficients(3) * pi) + 1.d0) / (phase_coefficients(3)**2 + 1.d0)
        normalization_factor = 1.d0 / (term1 * (term2 + term3 + term4))
    end subroutine

    subroutine update_zodi_spline_obj(bandpass)
        ! Updates the spline object which is used to evaluate b_nu over the bandpass

        implicit none
        class(comm_bp_ptr), intent(in) :: bandpass(:)

        real(dp), allocatable :: b_nu(:)
        real(dp) :: integrals(size(temperature_grid))
        integer(i4b) :: i, j, k, n_det

        n_det = size(bandpass) - 1
        ! Until the zodi stuff is moved into the tod class, we need to check if the spline object are allocated
        ! which should be done the first time the mixing matrix is updated.
        if (.not. allocated(b_nu_spline_obj)) allocate(b_nu_spline_obj(n_det))
        if (.not. allocated(splined_emissivities)) allocate(splined_emissivities(n_det, size(emissivities, dim=2)))
        if (.not. allocated(splined_albedos)) allocate(splined_albedos(n_det, size(albedos, dim=2)))
        if (.not. allocated(splined_phase_coeffs)) allocate(splined_phase_coeffs(n_det, 3))
        if (.not. allocated(splined_solar_irradiance)) allocate(splined_solar_irradiance(n_det))
        if (.not. allocated(phase_function_normalization)) allocate(phase_function_normalization(n_det))

        splined_emissivities = 0.d0
        splined_albedos = 0.d0
        splined_phase_coeffs = 0.d0
        splined_solar_irradiance = 0.d0
        phase_function_normalization = 0.d0
        
        do j = 1, n_det
            allocate(b_nu(bandpass(j)%p%n))
            do i = 1, n_interpolation_points    
                call get_blackbody_emission(bandpass(j)%p%nu, temperature_grid(i), b_nu)
                integrals(i) = bandpass(j)%p%SED2F(b_nu)
            end do
            call spline_simple(b_nu_spline_obj(j), temperature_grid, integrals)

            do k = 1, size(splined_emissivities)
                splined_emissivities(j, k) = splint_simple(emissivity_spline_obj(k), bandpass(j)%p%nu_c)
            end do

            do k = 1, size(splined_albedos)
                splined_albedos(j, k) = splint_simple(albedo_spline_obj(k), bandpass(j)%p%nu_c)
            end do
            if (count(splined_albedos /= 0.d0) > 0) then
                scattering = .true.
                do k = 1, size(splined_phase_coeffs, dim=2) 
                    splined_phase_coeffs(j, k) = splint_simple(phase_coeff_spline_obj(k), bandpass(j)%p%nu_c)
                end do

                splined_solar_irradiance(j) = splint_simple(solar_irradiance_spline_obj, bandpass(j)%p%nu_c)
                call get_phase_normalization(splined_phase_coeffs(j, :), phase_function_normalization(j))
            else 
                scattering = .false.
            end if
            deallocate(b_nu)
        end do

    end subroutine update_zodi_spline_obj

    ! Methods initizializing the zodiacal components
    ! ----------------------------------------------
    subroutine initialize_cloud(self)
        implicit none
        class(Cloud) :: self
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_cloud

    subroutine initialize_band(self)
        implicit none
        class(Band) :: self
        self%delta_zeta = self%delta_zeta * deg2rad
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_band

    subroutine initialize_ring(self)
        implicit none
        class(Ring) :: self
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_ring

    subroutine initialize_feature(self)
        implicit none
        class(Feature) :: self
        self%theta_0 = self%theta_0 * deg2rad
        self%sigma_theta = self%sigma_theta * deg2rad
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_feature


    ! Methods describing the three dimensional parametric density distribution of 
    ! the zodiacal components.
    ! ---------------------------------------------------------------------------
    subroutine get_density_cloud(self, X_vec, theta, n_out)
        implicit none
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
        implicit none
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
        implicit none
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
        implicit none
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


    ! Methods for constructing and iterating the linked list of zodiacal components
    ! -----------------------------------------------------------------------------
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
        implicit none
        class(ZodiComponent), pointer :: comp

        if (.not. associated(comp_list)) then
            comp_list => comp
        else
            call comp_list%add(comp)
        end if
    end subroutine add_component_to_list


    ! Utility functions
    ! ------------------
    function unique_sort(array) result(unique_sorted_array)
        implicit none
        integer :: idx, min_val, max_val
        integer, intent(in), dimension(:) :: array
        integer, dimension(size(array)) :: unique
        integer, dimension(:), allocatable :: unique_sorted_array

        idx = 0
        min_val = minval(array) - 1
        max_val = maxval(array)
        do while (min_val < max_val)
            idx = idx + 1
            min_val = minval(array, mask=array > min_val)
            unique(idx) = min_val
        enddo

        allocate(unique_sorted_array(idx), source=unique(1:idx))
    end function unique_sort

    subroutine gal_to_ecl_conversion_matrix(matrix)
        ! Ecliptic to galactic rotation matrix
        implicit none
        real(dp), dimension(3,3) :: matrix

        matrix(1,1) =  -0.054882486d0
        matrix(1,2) =  -0.993821033d0
        matrix(1,3) =  -0.096476249d0
        matrix(2,1) =   0.494116468d0
        matrix(2,2) =  -0.110993846d0
        matrix(2,3) =   0.862281440d0
        matrix(3,1) =  -0.867661702d0
        matrix(3,2) =  -0.000346354d0
        matrix(3,3) =   0.497154957d0

        ! call invert_matrix_dp(matrix)
    end subroutine gal_to_ecl_conversion_matrix

    subroutine get_tabulated_unit_vectors(npix, unit_vectors)
        ! Routine which selects coordinate transformation map based on resolution
        implicit none
        integer(i4b), intent(in) :: npix
        real(dp), dimension(:,:), intent(inout) :: unit_vectors

        integer(i4b) :: i

        do i = 1, size(unit_vector_list%vectors)
            if (size(unit_vector_list%vectors(i)%elements(:, 1)) == npix) then
                unit_vectors = unit_vector_list%vectors(i)%elements
                exit
            end if
        end do
    end subroutine get_tabulated_unit_vectors



    ! Functions for initializing zodi parameters from cpar
    ! ----------------------------------------------------
    subroutine init_hyper_parameters(cpar)
        ! Initialize hyper parameters for commander run
        implicit none 
        type(comm_params), intent(in) :: cpar
        
        use_cloud = cpar%zs_use_cloud
        use_band1 = cpar%zs_use_band1
        use_band2 = cpar%zs_use_band2
        use_band3 = cpar%zs_use_band3
        use_ring = cpar%zs_use_ring
        use_feature = cpar%zs_use_feature
        use_unit_emissivity = cpar%zs_use_unit_emissivity
        R_cutoff = cpar%zs_los_cut
        gauss_degree = cpar%zs_gauss_quad_order
    end subroutine init_hyper_parameters

    subroutine init_source_parameters(cpar)
        ! Initialize source parameters given interplanetary dust model
        implicit none 
        type(comm_params), intent(in) :: cpar

        T_0 = cpar%zs_t_0
        delta = cpar%zs_delta
        allocate(emissivities(cpar%zs_nbands, cpar%zs_ncomps))
        allocate(albedos(cpar%zs_nbands, cpar%zs_ncomps))
        allocate(phase_coeffs(cpar%zs_nbands, 3))
        allocate(solar_irradiances(cpar%zs_nbands))
        emissivities = cpar%zs_emissivity
        albedos = cpar%zs_albedo
        phase_coeffs = cpar%zs_phase_coeff
        solar_irradiances = cpar%zs_solar_irradiance
    end subroutine init_source_parameters

    subroutine init_shape_parameters(cpar)
        ! Initialize interplanetary dust shape parameters
        implicit none 
        type(comm_params), intent(in) :: cpar

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
    end subroutine init_shape_parameters
end module comm_zodi_mod