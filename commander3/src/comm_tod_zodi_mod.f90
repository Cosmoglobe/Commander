module comm_tod_zodi_mod
    use comm_param_mod
    use comm_utils
    use comm_tod_mod
    use comm_zodi_mod
    use spline_1D_mod
    use comm_bp_mod
    implicit none

    private
    public initialize_tod_zodi_mod, get_zodi_emission, update_zodi_splines

    ! Constants
    real(dp) :: R_MIN = 3.d-14
    real(dp) :: SECOND_TO_DAY = 1.1574074074074073e-05

    ! Shared global parameters
    type(spline_type) :: earth_pos_spl_obj(3)
    integer(i4b) :: n_interp_points, gauss_degree
    real(dp), allocatable :: temperature_grid(:), gauss_nodes(:), gauss_weights(:)
    real(dp) :: R_cutoff, delta_t_reset, min_ipd_temp, max_ipd_temp

contains 
    subroutine initialize_tod_zodi_mod(cpar)
        ! Initializes the zodiacal emission module.
        !
        ! Parameters
        ! ----------
        ! cpar : class(comm_param)
        !     The parameter object.

        type(comm_params), intent(in) :: cpar

        ! hyper parameters
        delta_t_reset = cpar%zs_delta_t_reset
        min_ipd_temp = cpar%zs_min_ipd_temp
        max_ipd_temp = cpar%zs_max_ipd_temp
        n_interp_points = cpar%zs_n_interp_points
        R_cutoff = cpar%zs_los_cut
        gauss_degree = cpar%zs_gauss_quad_order

        ! Set up Gauss-Legendre quadrature
        allocate(gauss_nodes(gauss_degree), gauss_weights(gauss_degree))
        call leggaus(gauss_degree, gauss_nodes, gauss_weights)

        ! Set up interpolation grid for evaluating temperature
        allocate(temperature_grid(n_interp_points))
        call linspace(min_ipd_temp, max_ipd_temp, temperature_grid)

        ! Read earth position from file and set up spline object
        call initialize_earth_pos_spline(cpar)

    end subroutine initialize_tod_zodi_mod


    subroutine get_zodi_emission(tod, pix, scan, det, s_zodi_scat, s_zodi_therm, model)
        ! Returns the predicted zodiacal emission for a scan (chunk of time-ordered data).
        !
        ! Parameters
        ! ----------
        ! tod : class(comm_tod)
        !     The TOD object holding the spline objects to update.
        ! pix : integer(i4b), dimension(ntod)
        !     The pixel indices of each time-ordered observation.
        ! det : integer(i4b)
        !     The detector index.
        ! scan : integer(i4b)
        !     The scan number.
        ! s_zodi_scat : real(sp), dimension(ntod, ncomps)
        !     Contribution from scattered sunlight light.
        ! s_zodi_therm : real(sp), dimension(ntod, ncomps)
        !     Contribution from thermal interplanetary dust emission.

        !
        ! Returns
        ! -------
        ! s_zodi_scat : real(sp), dimension(ntod, ncomps, ndet)
        !     Contribution from scattered sunlight light.
        ! s_zodi_therm : real(sp), dimension(ntod, ncomps, ndet)
        !     Contribution from thermal interplanetary dust emission.
        !
        class(comm_tod), intent(inout) :: tod
        integer(i4b), intent(in) :: pix(:), scan, det
        real(sp), dimension(:, :), intent(inout) :: s_zodi_scat, s_zodi_therm
        type(zodi_model), intent(in) :: model
        integer(i4b) :: i, j, k, pixel, lookup_idx, n_tod, ierr
        real(dp) :: earth_lon, R_obs, R_max, dt_tod, obs_time
        real(dp) :: unit_vector(3), X_unit_LOS(3, gauss_degree), X_LOS(3, gauss_degree), obs_pos(3), earth_pos(3)
        real(dp), dimension(gauss_degree) :: R_LOS, T_LOS, density_LOS, solar_flux_LOS, scattering_angle, phase_function, b_nu_LOS

        s_zodi_scat = 0.
        s_zodi_therm = 0.
        n_tod = size(pix, dim=1)

        dt_tod = (1./tod%samprate) * SECOND_TO_DAY ! dt between two samples in units of days (assumes equispaced tods)
        obs_pos = tod%scans(scan)%x0_obs
        earth_pos = tod%scans(scan)%x0_earth
        R_obs = norm2(obs_pos)
        obs_time = tod%scans(scan)%t0(1)   
        earth_lon = atan(earth_pos(2), earth_pos(1))

        do i = 1, n_tod
            ! Reset cache if time between last cache update and current time is larger than `delta_t_reset`.
            ! NOTE: this cache is only effective if the scans a core handles are in chronological order.
            obs_time = obs_time + dt_tod
            if ((obs_time - tod%zodi_cache_time) >= delta_t_reset) then 
                do j = 1, 3 
                    earth_pos(j) = splint_simple(tod%x_earth_spline(j), obs_time)
                    obs_pos(j) = splint_simple(tod%x_obs_spline(j), obs_time)
                end do            
                R_obs = norm2(obs_pos)
                earth_lon = atan(earth_pos(2), earth_pos(1))
                call tod%clear_zodi_cache(obs_time)
            end if

            lookup_idx = tod%pix2ind(pix(i))
            if (tod%zodi_therm_cache(lookup_idx, 1, det) > 0.d0) then
                s_zodi_scat(i, :) = tod%zodi_scat_cache(lookup_idx, :, det)
                s_zodi_therm(i, :) = tod%zodi_therm_cache(lookup_idx, :, det)
                cycle
            end if

            unit_vector = tod%ind2vec_ecl(:, lookup_idx)
            call get_R_max(unit_vector, obs_pos, R_obs, R_max)

            do k = 1, 3
                ! Convert quadrature range from [-1, 1] to [R_min, R_max]
                X_unit_LOS(k, :) = (0.5 * (R_max - R_MIN)) * gauss_nodes + (0.5 * (R_max + R_MIN))
                X_unit_LOS(k, :) = X_unit_LOS(k, :) * unit_vector(k)
                X_LOS(k, :) = X_unit_LOS(k, :) + obs_pos(k)
            end do
            R_LOS = norm2(X_LOS, dim=1)           

            if (tod%zodi_scattering) then
                solar_flux_LOS = tod%zodi_spl_solar_irradiance(det) / R_LOS**2
                call get_scattering_angle(X_LOS, X_unit_LOS, R_LOS, scattering_angle)
                call get_phase_function(scattering_angle, tod%zodi_spl_phase_coeffs(det, :), tod%zodi_phase_func_normalization(det), phase_function)
            end if

            call get_dust_grain_temperature(R_LOS, T_LOS, model%T_0, model%delta)
            call splint_simple_multi(tod%zodi_b_nu_spl_obj(det), T_LOS, b_nu_LOS)

            do k = 1, model%n_comps
                call model%comps(k)%c%get_density(X_LOS, earth_lon, density_LOS)
                if (tod%zodi_scattering) then 
                    s_zodi_scat(i, k) = sum(density_LOS * solar_flux_LOS * phase_function * gauss_weights)
                    tod%zodi_scat_cache(lookup_idx, k, det) = s_zodi_scat(i, k)
                end if
                s_zodi_therm(i, k) = sum(density_LOS * b_nu_LOS * gauss_weights) * 0.5 * (R_max - R_MIN)
                tod%zodi_therm_cache(lookup_idx, k, det) = s_zodi_therm(i, k)
            end do
        end do
    end subroutine get_zodi_emission

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
        b = 2. * (obs_pos(1) * cos_lat * cos(lon) + obs_pos(2) * cos_lat * sin(lon))
        d = R_obs**2 - R_cutoff**2
        q = -0.5 * b * (1. + sqrt(b**2 - (4. * d)) / abs(b))
        R_max = max(q, d / q)
    end subroutine get_R_max

    subroutine get_dust_grain_temperature(R, T_out, T_0, delta)
        real(dp), dimension(:), intent(in) :: R
        real(dp), dimension(:), intent(out) :: T_out
        real(dp), intent(in) :: T_0, delta
        T_out = T_0 * R ** (-delta)
    end subroutine get_dust_grain_temperature

    subroutine get_blackbody_emission(nus, T, b_nu)
        real(dp), intent(in) :: nus(:), T
        real(dp), dimension(:), intent(out) :: b_nu
        b_nu = 1d20 * ((2. * h * nus**3) / (c*c)) / (exp((h * nus) / (k_B * T)) - 1.)
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

        term1 = 2. * pi
        term2 = 2. * phase_coefficients(1)
        term3 = pi * phase_coefficients(2)
        term4 = (exp(phase_coefficients(3) * pi) + 1.) / (phase_coefficients(3)**2 + 1.)
        normalization_factor = 1. / (term1 * (term2 + term3 + term4))
    end subroutine


    subroutine initialize_earth_pos_spline(cpar)
        ! Returns the spline object which is used to evaluate the earth position

        type(comm_params), intent(in) :: cpar

        integer :: i, n_earthpos, unit
        real(dp), allocatable :: tabulated_earth_time(:), tabulated_earth_pos(:, :)
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
            call spline_simple(earth_pos_spl_obj(i), tabulated_earth_time, tabulated_earth_pos(i, :), regular=.true.)
        end do
    end subroutine initialize_earth_pos_spline


    subroutine update_zodi_splines(tod, bandpass, det, model)
        ! Updates the spectral spline objects in the TOD object.
        !
        ! In the K98 model, several spectral parameters are tabulated at individual frequencies, 
        ! which we need to evaluate over the bandpass. In a future version, we may want to fit 
        ! a modified blackbody which would allow us to drop using some of these spline objects.
        !
        !  -----------------------------------------------------------------------------------------
        ! | The difficulty with this functino is that it needs the access to the bandpass, so is is |
        ! | very limited in where it can be excecuted in commander.                                 |
        !  -----------------------------------------------------------------------------------------
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
        type(zodi_model), intent(inout) :: model

        real(dp), allocatable :: b_nu(:), integrals(:)
        integer(i4b) :: i, j
    
        allocate(integrals(n_interp_points))
        allocate(b_nu(bandpass%p%n))
        do i = 1, n_interp_points    
           !write(*,*) 'alloc2', allocated(bandpass%p%nu)
            call get_blackbody_emission(bandpass%p%nu, temperature_grid(i), b_nu)
            integrals(i) = tsum(bandpass%p%nu, bandpass%p%tau * b_nu)
        end do
        call spline_simple(tod%zodi_b_nu_spl_obj(det), temperature_grid, integrals, regular=.true.)
        do j = 1, size(tod%zodi_spl_phase_coeffs, dim=2) 
            tod%zodi_spl_phase_coeffs(det, j) = splint_simple(model%phase_coeff_spl(j), bandpass%p%nu_c)
        end do

        tod%zodi_spl_solar_irradiance = splint_simple(model%solar_irradiance_spl, bandpass%p%nu_c)
        call get_phase_normalization(tod%zodi_spl_phase_coeffs(det, :), tod%zodi_phase_func_normalization(det))
    end subroutine update_zodi_splines

end module comm_tod_zodi_mod
