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
    !   NOTE: in the current setup, the zodi module only accepts IPD components as
    !   described in the DIRBE IPD model. All model parameters are provided in the 
    !   parameter files. By default, the IPD shape, source and hyper parameters in
    !   "parameter_files/defaults/components/zodi/*" are used, but these can be 
    !   overwritten your own parameterfile.
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
    public :: initialize_zodi_mod, get_zodi_emission
    integer(i4b) :: gauss_quad_order
    real(dp) :: R_LOS_cutoff, EPS, delta_t_reset_cash, previous_chunk_obs_time
    real(dp), allocatable, dimension(:) :: unique_nsides
    real(dp), allocatable, dimension(:) :: tabulated_earth_time
    real(sp), allocatable, dimension(:) :: cached_zodi
    real(dp), allocatable, dimension(:,:) :: tabulated_earth_pos
    character(len=512) :: freq_correction_type
    type(spline_type) :: solar_irradiance_spline_obj
    type(spline_type), dimension(3) :: spline_earth_pos_obj
    type(spline_type), allocatable, dimension(:) :: emissivity_spline_obj, albedo_spline_obj, phase_coeff_spline_obj
    logical(lgt) :: use_cloud, use_band1, use_band2, use_band3, use_ring, use_feature, &
                    apply_color_correction, use_unit_emissivity, use_unit_albedo

    ! model parameters
    real(dp) :: T_0, delta, splined_solar_irradiance
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
    real(dp), allocatable, dimension(:) :: solar_irradiances
    real(dp), allocatable, dimension(:, :) :: emissivities, albedos, phase_functions
    real(dp), allocatable, dimension(:) :: splined_emissivities, splined_albedos, splined_phase_coeffs


    type, abstract :: ZodiComponent
        ! Abstract base Zodical component class.

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

        subroutine density_interface(self, r_vec, theta, n_out)
            ! Returns the dust density (n) of the component at heliocentric 
            ! coordinates (x, y, z) and the earths longitude (theta).

            import i4b, dp, ZodiComponent
            class(ZodiComponent) :: self
            real(dp), intent(in), dimension(:, :) :: r_vec ! shape: (3, gauss_quadorder)
            real(dp), intent(in) :: theta
            real(dp), intent(out), dimension(:) :: n_out
            real(dp) :: r_prime_vec, R, Z_midplane
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

    ! Initializing list to o
    class(ZodiComponent), pointer :: comp_list => null()
    type(Cloud),          target  :: cloud_comp
    type(Band),           target  :: band1_comp, band2_comp, band3_comp
    type(Ring),           target  :: ring_comp
    type(Feature),        target  :: feature_comp

    ! Initializing container for precomputed pixel to unit vectors for each unique 
    ! nside in the Commander run.
    type UnitVector
        real(dp), dimension(:,:), allocatable :: elements
    end type UnitVector

    type UnitVectorList
        type(UnitVector), dimension(:), allocatable :: vectors
    end type UnitVectorList
    type(UnitVectorList) :: unit_vectors


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
        real(dp), dimension(3) :: vec
        real(dp), dimension(3,3) :: ecliptic_to_galactic_matrix
        integer(i4b), dimension(:), allocatable :: sorted_unique_nsides
        real(dp), dimension(:,:), allocatable :: galactic_vec

        EPS = 3.d-14

        ! Initialize hyper, shape, and source parameters for ipd model from cpar
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
        allocate(unit_vectors%vectors(size(sorted_unique_nsides)))

        do i = 1, size(sorted_unique_nsides)
            nside = sorted_unique_nsides(i)
            npix = nside2npix(nside)
            allocate(unit_vectors%vectors(i)%elements(0:npix-1, 3))
            allocate(galactic_vec(0:npix-1, 3))
        
            do j = 0, npix - 1
                call pix2vec_ring(nside, j, vec)
                galactic_vec(j, 1) = vec(1)
                galactic_vec(j, 2) = vec(2)
                galactic_vec(j, 3) = vec(3)
            end do
        
            unit_vectors%vectors(i)%elements = matmul(galactic_vec, ecliptic_to_galactic_matrix)
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
        allocate(splined_emissivities(cpar%zs_ncomps))
        allocate(splined_albedos(cpar%zs_ncomps))
        allocate(splined_phase_coeffs(3))

        ! Earths position
        do i = 1, 3
            call spline_simple(spline_earth_pos_obj(i), tabulated_earth_time, tabulated_earth_pos(i, :), regular=.true.)
        end do

        ! Source parameters
        do i = 1, cpar%zs_ncomps
            if (.not. use_unit_emissivity) then
                call spline_simple(emissivity_spline_obj(i), cpar%zs_nu_ref, cpar%zs_emissivity(:, i), regular=.false.)
            else 
                splined_emissivities = 1.d0
            end if 
            if (.not. use_unit_albedo) then
                call spline_simple(albedo_spline_obj(i), cpar%zs_nu_ref, cpar%zs_albedo(:, i), regular=.false.)
            else
                splined_albedos = 1.d0
            end if
        end do
        do i = 1, 3
            call spline_simple(phase_coeff_spline_obj(i), cpar%zs_nu_ref, cpar%zs_phase_coeff(:, i), regular=.false.)
        end do
        call spline_simple(solar_irradiance_spline_obj, cpar%zs_nu_ref, cpar%zs_solar_irradiance, regular=.false.)


        previous_chunk_obs_time = 0 ! Set initial previous chunk observation time to 0
    end subroutine initialize_zodi_mod

    subroutine get_zodi_emission(nside, pix, sat_pos, obs_time, bandpass, s_zodi)
        !   Simulates the zodiacal emission over a chunk of time-ordered data.
        !
        !   Given a set of observations, the observer position and time of 
        !   observation is used to perform line-of-sight integrations through 
        !   the interplanetary dust distribution in the Solar System, yielding 
        !   a prediction of the zodiacal emission seen by the observer.
        !
        !   Parameters
        !   ----------
        !   nside: int
        !       Resolution of the pixels given in `pix`.
        !   pix: array
        !       Pixel indices representing observations in a chunk of the 
        !       time-ordered data. Dimensions are (`n_tod`, `n_det`)
        !   sat_pos: real
        !       Heliocentric ecliptic cartesian position of the observer in AU.
        !   obs_time: real
        !       Time of observation in MJD. Assumes a fixed time of observeration
        !       for the entire tod chunk. Chunk sizes exceeding 1 day in time will
        !       result in poor zodiacal emission estimates.
        !   bandpass: comm_bp_ptr
        !       Bandpass object used to apply bandpass or color corrections to the
        !       predicted zodiacal emission.
        !
        !   Returns
        !   -------
        !   s_zodi: array
        !       Estimated zodiacal emission over the chunk of time-ordered data.
        !       Dimensions are (`n_tod`, `n_det`).
        implicit none
        class(ZodiComponent), pointer :: comp

        integer(i4b), intent(in) :: nside
        integer(i4b), dimension(1:,1:), intent(in) :: pix
        real(dp), dimension(3), intent(in) :: sat_pos
        real(dp), intent(in) :: obs_time
        class(comm_bp_ptr), dimension(:), intent(in) :: bandpass
        real(sp), dimension(1:,1:), intent(out) :: s_zodi

        integer(i4b) :: i, j, k, pixel_idx, los_step, n_detectors, n_tods, npix
        logical(lgt) :: scattering
        real(dp), dimension(3) :: u_vec, r_obs_vec, r_earth_vec
        real(dp), dimension(gauss_quad_order) :: R_helio_LOS
        real(dp), dimension(3, gauss_quad_order) :: r_vec_LOS, r_helio_vec_LOS
        real(dp) :: theta_earth, R_obs, R_max, nu_det
        real(dp) :: b_nu_delta_term1, b_nu_delta_term2
        real(dp) :: phase_function_normalization = 0.d0
        real(dp), dimension(:), allocatable :: nu_ratio, b_nu_bandpass_term1, b_nu_bandpass_term2
        real(dp), dimension(:,:), allocatable :: unit_vector_map, b_nu_bandpass_LOS, b_nu_ratio_LOS
        real(dp), dimension(gauss_quad_order) :: gauss_grid, gauss_weights
        real(dp), dimension(gauss_quad_order) :: T_LOS, density_LOS                              
        real(dp), dimension(gauss_quad_order) :: comp_emission_LOS, b_nu_bandpass_integrated_LOS, &
                                                 b_nu_center_LOS, b_nu_colorcorr_LOS, b_nu_freq_corrected_LOS
        real(dp), dimension(gauss_quad_order) :: solar_flux_LOS, scattering_angle, phase_function

        n_tods = size(pix,1)
        n_detectors = size(pix,2)
        npix = nside2npix(nside)

        ! Get unit vectors corresponding to the observed pixels (only on first chunk)
        if (.not. allocated(unit_vector_map)) then 
            allocate(unit_vector_map(0:npix-1, 3))
            call get_unit_vector_map(nside, unit_vector_map)
        end if

        ! Allocate cached zodi array (only on first chunk)
        if (.not. allocated(cached_zodi)) then
            allocate(cached_zodi(0:npix-1))
            cached_zodi = 0.d0
        end if

        ! Reset quantites from previous chunk
        comp_emission_LOS = 0.d0
        gauss_grid = 0.d0
        gauss_weights = 0.d0
        s_zodi = 0.d0

        ! Reset cached zodi if time since last chunk > DELTA_T_ZODI ~ 1day. 
        if ((obs_time - previous_chunk_obs_time) > delta_t_reset_cash) cached_zodi = 0.d0

        ! Get Earths position given `obs_time`
        do i = 1, 3 
            r_earth_vec(i) = splint_simple(spline_earth_pos_obj(i), obs_time)
        end do

        r_obs_vec = sat_pos
        R_obs = norm2(r_obs_vec)

        ! Get longitude of the earth
        theta_earth = atan(r_earth_vec(2), r_earth_vec(1))


        ! Loop over each detectors time-ordered pointing chunks. For each unique pixel
        ! observed (cross dectors) perform a line of sight integral solving Eq. (20) in 
        ! San et al. 2022 (https://www.aanda.org/articles/aa/pdf/forth/aa44133-22.pdf).
        do j = 1, n_detectors
            ! Allocate detector specific blackbody quantities
            allocate(b_nu_bandpass_term1(bandpass(j)%p%n))
            allocate(b_nu_bandpass_term2(bandpass(j)%p%n))
            allocate(b_nu_bandpass_LOS(gauss_quad_order, bandpass(j)%p%n))
            allocate(b_nu_ratio_LOS(gauss_quad_order, bandpass(j)%p%n))
            allocate(nu_ratio(bandpass(j)%p%n))
            b_nu_delta_term1 = (2 * h * bandpass(j)%p%nu_c**3) / (c*c)
            b_nu_delta_term2 = (h * bandpass(j)%p%nu_c)/ k_B
            b_nu_bandpass_term1 = (2 * h * bandpass(j)%p%nu**3) / (c*c)
            b_nu_bandpass_term2 = (h * bandpass(j)%p%nu)/ k_B

            ! Get interpolated source parameters
            if (.not. use_unit_emissivity) then
                do k = 1, size(splined_emissivities)
                    splined_emissivities(k) = splint_simple(emissivity_spline_obj(k), bandpass(j)%p%nu_c)
                end do
            end if
            if (.not. use_unit_albedo) then
                do k = 1, size(splined_albedos)
                    splined_albedos(k) = splint_simple(albedo_spline_obj(k), bandpass(j)%p%nu_c)
                end do
            end if
            do k = 1, size(splined_phase_coeffs)
                splined_phase_coeffs(k) = splint_simple(phase_coeff_spline_obj(k), bandpass(j)%p%nu_c)
            end do
            splined_solar_irradiance = splint_simple(solar_irradiance_spline_obj, bandpass(j)%p%nu_c)

            ! If any of the interpolated albedos are non zero, we must take contributions from 
            ! scattered sunlight into account when computing the zodiacal emission
            if (count(splined_albedos /= 0.d0) > 0) then
                scattering = .true.
                call get_phase_normalization(splined_phase_coeffs, phase_function_normalization)
            else
                scattering = .false.
            end if

            do i = 1, n_tods
                pixel_idx = pix(i, j)
                if (cached_zodi(pixel_idx) /= 0.d0) then
                    s_zodi(i, j) = cached_zodi(pixel_idx)
                else                  
                    u_vec = unit_vector_map(pixel_idx, :)      

                    call get_R_max(u_vec, r_obs_vec, R_obs, R_LOS_cutoff, R_max)
                    call gauss_legendre_quadrature(x1=EPS, x2=R_max, n=gauss_quad_order, x=gauss_grid, w=gauss_weights)

                    ! Line of sight grid points from solar system origin with shape (3, gauss_quad_order)
                    r_vec_LOS(1, :) = gauss_grid * u_vec(1)
                    r_vec_LOS(2, :) = gauss_grid * u_vec(2)
                    r_vec_LOS(3, :) = gauss_grid * u_vec(3)

                    ! Line of sight grid points from observer in heliocentric coordinates with shape (3, gauss_quad_order)
                    r_helio_vec_LOS(1, :) = r_vec_LOS(1, :) + r_obs_vec(1)
                    r_helio_vec_LOS(2, :) = r_vec_LOS(2, :) + r_obs_vec(2)
                    r_helio_vec_LOS(3, :) = r_vec_LOS(3, :) + r_obs_vec(3)
                    R_helio_LOS = norm2(r_helio_vec_LOS, dim=1)
    
                    call get_dust_grain_temperature(R_helio_LOS, T_LOS)

                    ! Compute blackbody LOS terms which depend on the frequency correction specified in the hyper parameters
                    select case (trim(freq_correction_type))
                        case ("delta")        
                            call get_blackbody_emission_delta(T_LOS, b_nu_delta_term1, b_nu_delta_term2, b_nu_center_LOS)
                            b_nu_freq_corrected_LOS = b_nu_center_LOS
                        case ("bandpass")
                            call get_blackbody_emission_bp(T_LOS, b_nu_bandpass_term1, b_nu_bandpass_term2, b_nu_bandpass_LOS)
                            do los_step = 1, gauss_quad_order
                                b_nu_bandpass_integrated_LOS(los_step) = bandpass(j)%p%SED2F(b_nu_bandpass_LOS(los_step, :))
                            end do
                            b_nu_freq_corrected_LOS = b_nu_bandpass_integrated_LOS
                        case("color")
                            call get_blackbody_emission_delta(T_LOS, b_nu_delta_term1, b_nu_delta_term2, b_nu_center_LOS)
                            call get_blackbody_emission_bp(T_LOS, b_nu_bandpass_term1, b_nu_bandpass_term2, b_nu_bandpass_LOS)
                            do los_step = 1, gauss_quad_order
                                b_nu_ratio_LOS(los_step, :) = b_nu_bandpass_LOS(los_step, :) / b_nu_center_LOS(los_step)
                                b_nu_colorcorr_LOS(los_step) = tsum(bandpass(j)%p%nu, b_nu_ratio_LOS(los_step, :) * bandpass(j)%p%tau)
                            end do
                            nu_ratio = bandpass(j)%p%nu_c / bandpass(j)%p%nu
                            b_nu_colorcorr_LOS = b_nu_colorcorr_LOS / tsum(bandpass(j)%p%nu, nu_ratio * bandpass(j)%p%tau)
                            b_nu_freq_corrected_LOS = b_nu_center_LOS * b_nu_colorcorr_LOS
                    end select

                    if (scattering) then
                        solar_flux_LOS = splined_solar_irradiance / R_helio_LOS**2
                        call get_scattering_angle(r_helio_vec_LOS, R_helio_LOS, r_vec_LOS, scattering_angle)
                        call get_phase_function(scattering_angle, splined_phase_coeffs, phase_function_normalization, phase_function)
                    end if

                    comp => comp_list
                    k = 1
                    do while (associated(comp))
                        call comp%get_density(r_helio_vec_LOS, theta_earth, density_LOS)
                        comp_emission_LOS = (1.d0 - splined_albedos(k)) * (splined_emissivities(k) * b_nu_freq_corrected_LOS)
                        if (scattering) comp_emission_LOS = comp_emission_LOS + (splined_albedos(k) * solar_flux_LOS * phase_function)
                        comp_emission_LOS = comp_emission_LOS * density_LOS
                        s_zodi(i, j) = s_zodi(i, j) + sum(comp_emission_LOS * gauss_weights)

                        comp => comp%next()
                        k = k + 1
                    end do
                    cached_zodi(pixel_idx) = s_zodi(i, j) ! Update cache with newly computed LOS emission
                end if
            end do
            deallocate(b_nu_bandpass_term1)
            deallocate(b_nu_bandpass_term2)
            deallocate(b_nu_bandpass_LOS)
            deallocate(b_nu_ratio_LOS)
            deallocate(nu_ratio)
        end do
        
        previous_chunk_obs_time = obs_time ! Store prevous chunks obs time
    end subroutine get_zodi_emission


    ! Functions used in `get_zodi_emission`
    ! -------------------------------------
    subroutine get_R_max(u_vec, r_obs_vec, R_obs, R_cutoff, R_max)
        ! Computes R_max (the length of the LOS such that it stops exactly at R_cutoff).
        implicit none
        real(dp), dimension(:), intent(in) :: u_vec, r_obs_vec
        real(dp), intent(in) :: R_obs, R_cutoff
        real(dp), intent(out) :: R_max
        real(dp) :: lon, lat, cos_lat, b, d, q

        lon = atan(u_vec(2), u_vec(1))
        lat = asin(u_vec(3))
        cos_lat = cos(lat)

        b = 2.d0 * (r_obs_vec(1) * cos_lat * cos(lon) + r_obs_vec(2) * cos_lat * sin(lon))
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

    subroutine get_blackbody_emission_bp(T, b_nu_bandpass_term1, b_nu_bandpass_term2, b_nu_out)
        implicit none
        real(dp), dimension(:), intent(in) :: T, b_nu_bandpass_term1, b_nu_bandpass_term2
        real(dp), dimension(:, :), intent(out) :: b_nu_out
        integer(i4b) :: i
        do i = 1, gauss_quad_order
            b_nu_out(i, :) = b_nu_bandpass_term1/(exp(b_nu_bandpass_term2/T(i)) - 1.d0)
        end do
        b_nu_out = b_nu_out * 1d20 !Convert from W/s/m^2/sr to MJy/sr
    end subroutine get_blackbody_emission_bp

    subroutine get_blackbody_emission_delta(T, b_nu_delta_term1, b_nu_delta_term2, b_nu_out)
        implicit none
        real(dp), dimension(:), intent(in) :: T
        real(dp), intent(in) :: b_nu_delta_term1, b_nu_delta_term2
        real(dp), dimension(:), intent(out) :: b_nu_out
        b_nu_out = (b_nu_delta_term1/(exp(b_nu_delta_term2/T) - 1.d0)) * 1d20 !Convert from W/s/m^2/sr to MJy/sr
    end subroutine get_blackbody_emission_delta

    subroutine get_scattering_angle(r_helio_vec_LOS, R_helio_LOS, r_vec_LOS, scattering_angle)
        implicit none
        real(dp), dimension(:, :), intent(in) :: r_helio_vec_LOS, r_vec_LOS
        real(dp), dimension(:), intent(in) :: R_helio_LOS
        real(dp), dimension(:), intent(out) :: scattering_angle
        real(dp), dimension(gauss_quad_order) :: R_LOS, cos_theta

        R_LOS = norm2(r_vec_LOS, dim=1)
        cos_theta = (sum(r_helio_vec_LOS * r_vec_LOS, dim=1))  / (R_helio_LOS * R_LOS)

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
        real(dp), intent(in) :: scattering_angle(:), phase_coefficients(:)
        real(dp), intent(in) :: normalization_factor
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
        freq_correction_type = cpar%zs_freq_correction_type
        R_LOS_cutoff = cpar%zs_los_cut
        gauss_quad_order = cpar%zs_gauss_quad_order
        delta_t_reset_cash = cpar%zs_delta_t
    end subroutine init_hyper_parameters

    subroutine init_source_parameters(cpar)
        ! Initialize source parameters given interplanetary dust model
        implicit none 
        type(comm_params), intent(in) :: cpar

        T_0 = cpar%zs_t_0
        delta = cpar%zs_delta
        allocate(emissivities(cpar%zs_nbands, cpar%zs_ncomps))
        allocate(albedos(cpar%zs_nbands, cpar%zs_ncomps))
        allocate(phase_functions(cpar%zs_nbands, 3))
        allocate(solar_irradiances(cpar%zs_nbands))
        emissivities = cpar%zs_emissivity
        albedos = cpar%zs_emissivity
        phase_functions = cpar%zs_emissivity
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

        ring_x_0 = cpar%zs_common(4, 1)
        ring_y_0 = cpar%zs_common(4, 2)
        ring_z_0 = cpar%zs_common(4, 3)
        ring_incl = cpar%zs_common(4, 4)
        ring_Omega = cpar%zs_common(4, 5)
        ring_n_0 = cpar%zs_common(4, 6)
        ring_R = cpar%zs_ring_r
        ring_sigma_R = cpar%zs_ring_sigma_r
        ring_sigma_z = cpar%zs_ring_sigma_z

        feature_x_0 = cpar%zs_common(5, 1)
        feature_y_0 = cpar%zs_common(5, 2)
        feature_z_0 = cpar%zs_common(5, 3)
        feature_incl = cpar%zs_common(5, 4)
        feature_Omega = cpar%zs_common(5, 5)
        feature_n_0 = cpar%zs_common(5, 6)
        feature_R = cpar%zs_feature_r
        feature_sigma_R = cpar%zs_feature_sigma_r
        feature_sigma_z = cpar%zs_feature_sigma_z
        feature_theta = cpar%zs_feature_theta
        feature_sigma_theta = cpar%zs_feature_sigma_theta
    end subroutine init_shape_parameters


    ! Methods for the zodiacal components
    ! -----------------------------------
    subroutine initialize_cloud(self)
        implicit none
        class(Cloud) :: self
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_cloud

    subroutine get_density_cloud(self, r_vec, theta, n_out)
        implicit none
        class(Cloud) :: self
        real(dp), dimension(:, :), intent(in) :: r_vec
        real(dp), intent(in) :: theta
        real(dp), dimension(:), intent(out) :: n_out
        integer(i4b) :: i
        real(dp) :: R_prime, Z_prime, g, zeta
        real(dp), dimension(3) :: r_prime_vec

        do i = 1, gauss_quad_order
            r_prime_vec(1) = r_vec(1, i) - self%x_0
            r_prime_vec(2) = r_vec(2, i) - self%y_0
            r_prime_vec(3) = r_vec(3, i) - self%z_0

            R_prime = norm2(r_prime_vec, dim=1)
            Z_prime = (r_prime_vec(1)*self%sin_omega - r_prime_vec(2)*self%cos_omega)*self%sin_incl &
                      + r_prime_vec(3)*self%cos_incl
            zeta = abs(Z_prime/R_prime)

            if (zeta < self%mu) then
                g = (zeta * zeta) / (2.d0 * self%mu)
            else
                g = zeta - (0.5d0 * self%mu)
            end if

            n_out(i) = self%n_0 * R_prime**(-self%alpha) * exp(-self%beta * g**self%gamma)
        end do
    end subroutine get_density_cloud

    subroutine initialize_band(self)
        implicit none
        class(Band) :: self
        self%delta_zeta = self%delta_zeta * deg2rad
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_band

    subroutine get_density_band(self, r_vec, theta, n_out)
        implicit none
        class(Band) :: self
        real(dp), dimension(:, :), intent(in) :: r_vec
        real(dp), intent(in) :: theta
        real(dp), dimension(:), intent(out) :: n_out
        integer(i4b) :: i
        real(dp) :: R_prime, Z_prime, zeta, zeta_over_delta_zeta, term1, term2, term3, term4
        real(dp), dimension(3) :: r_prime_vec

        do i = 1, gauss_quad_order
            r_prime_vec(1) = r_vec(1, i) - self%x_0
            r_prime_vec(2) = r_vec(2, i) - self%y_0
            r_prime_vec(3) = r_vec(3, i) - self%z_0

            R_prime = norm2(r_prime_vec, dim=1)
            Z_prime = (r_prime_vec(1)*self%sin_omega - r_prime_vec(2)*self%cos_omega)*self%sin_incl &
                      + r_prime_vec(3)*self%cos_incl
            zeta = abs(Z_prime/R_prime)

            zeta_over_delta_zeta = zeta / self%delta_zeta
            term1 = (3.d0 * self%n_0) / R_prime
            term2 = exp(-(zeta_over_delta_zeta**6))

            ! Differs from eq 8 in K98 by a factor of 1/self.v. See Planck XIV
            ! section 4.1.2.
            term3 = 1.d0 + (zeta_over_delta_zeta**self%p) / self%v
            term4 = 1.d0 - exp(-((R_prime / self%delta_r) ** 20))

            n_out(i) = term1 * term2 * term3 * term4
        end do
    end subroutine get_density_band

    subroutine initialize_ring(self)
        implicit none
        class(Ring) :: self
        self%sin_omega = sin(self%Omega * deg2rad)
        self%cos_omega = cos(self%Omega * deg2rad)
        self%sin_incl = sin(self%incl * deg2rad)
        self%cos_incl = cos(self%incl * deg2rad)
    end subroutine initialize_ring

    subroutine get_density_ring(self, r_vec, theta, n_out)
        implicit none
        class(Ring) :: self
        real(dp), dimension(:, :), intent(in) :: r_vec
        real(dp), intent(in) :: theta
        real(dp), dimension(:), intent(out) :: n_out
        integer(i4b) :: i
        real(dp) :: R_prime, Z_prime, zeta, term1, term2
        real(dp), dimension(3) :: r_prime_vec

        do i = 1, gauss_quad_order
            r_prime_vec(1) = r_vec(1, i) - self%x_0
            r_prime_vec(2) = r_vec(2, i) - self%y_0
            r_prime_vec(3) = r_vec(3, i) - self%z_0

            R_prime = norm2(r_prime_vec, dim=1)
            Z_prime = (r_prime_vec(1)*self%sin_omega - r_prime_vec(2)*self%cos_omega)*self%sin_incl &
                      + r_prime_vec(3)*self%cos_incl
            zeta = abs(Z_prime/R_prime)

            term1 = -((R_prime - self%R_0) ** 2) / self.sigma_r**2
            term2 = abs(Z_prime/self.sigma_z)

            n_out(i) = self%n_0 * exp(term1 - term2)
        end do
    end subroutine get_density_ring

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

    subroutine get_density_feature(self, r_vec, theta, n_out)
        implicit none
        class(Feature) :: self
        real(dp), dimension(:, :), intent(in) :: r_vec
        real(dp), intent(in) :: theta
        real(dp), dimension(:), intent(out) :: n_out
        integer(i4b) :: i
        real(dp) :: R_prime, Z_prime, zeta, theta_prime, exp_term
        real(dp), dimension(3) :: r_prime_vec

        do i = 1, gauss_quad_order
            r_prime_vec(1) = r_vec(1, i) - self%x_0
            r_prime_vec(2) = r_vec(2, i) - self%y_0
            r_prime_vec(3) = r_vec(3, i) - self%z_0

            R_prime = norm2(r_prime_vec, dim=1)
            Z_prime = (r_prime_vec(1)*self%sin_omega - r_prime_vec(2)*self%cos_omega)*self%sin_incl &
                      + r_prime_vec(3)*self%cos_incl
            zeta = abs(Z_prime/R_prime)

            theta_prime = atan2(r_prime_vec(2), r_prime_vec(1)) - theta - self%theta_0

            ! clip theta_prime to [-pi, pi]
            do while (theta_prime < -pi)
                theta_prime = theta_prime + 2.d0*pi
            end do
            do while (theta_prime > pi)
                theta_prime = theta_prime - 2.d0*pi
            end do
            
            exp_term = ((R_prime - self%R_0) ** 2 / self%sigma_r**2) + (abs(Z_prime) / self%sigma_z) + (theta_prime**2 / self%sigma_theta**2)

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


    ! Auxilary functions
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

    subroutine get_unit_vector_map(band_nside, unit_vector_map)
        ! Routine which selects coordinate transformation map based on resolution
        implicit none
        integer(i4b), intent(in) :: band_nside
        real(dp), dimension(:,:), intent(inout) :: unit_vector_map
        integer(i4b) :: i, npix, nside_idx

        npix = nside2npix(band_nside)
        do i = 1, size(unit_vectors%vectors)
            if (size(unit_vectors%vectors(i)%elements(:, 1)) == npix) then
                unit_vector_map = unit_vectors%vectors(i)%elements
            end if
        end do
    end subroutine get_unit_vector_map
end module comm_zodi_mod