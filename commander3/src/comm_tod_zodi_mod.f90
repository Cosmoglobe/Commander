module comm_tod_zodi_mod
    use comm_param_mod
    use comm_utils
    use comm_tod_mod
    use comm_zodi_mod
    use spline_1D_mod
    use comm_bp_mod

    implicit none

    private
    public get_zodi_emission, update_zodi_splines, initialize_tod_zodi_mod, solve_Ax_zodi, sample_dynamic_zodi_parameters, accumulate_zodi_emissivities, accumulate_zodi_albedos, downsample_zodi_res_and_pointing, sample_zodi_params

    ! Constants
    real(dp) :: EPS = 3.d-14
    real(dp) :: SECOND_TO_DAY = 1.1574074074074073e-05

    ! Shared global parameters
    type(spline_type) :: earth_pos_spl_obj(3)
    real(dp), allocatable :: temperature_grid(:)
    real(dp) :: R_cutoff, delta_t_reset, min_ipd_temp, max_ipd_temp
    integer(i4b) :: n_interp_points, gauss_degree

contains 
    subroutine initialize_tod_zodi_mod(cpar, tod)
        ! Initializes the zodiacal emission module.
        !
        ! Parameters
        ! ----------
        ! cpar : class(comm_param)
        !     The parameter object.

        type(comm_params), intent(in) :: cpar
        class(comm_tod), intent(inout) :: tod
        integer(i4b) :: i, ierr
        real(dp), allocatable :: obs_time(:), obs_pos(:, :), r

        ! hyper parameters
        delta_t_reset = cpar%zs_delta_t_reset
        min_ipd_temp = cpar%zs_min_ipd_temp
        max_ipd_temp = cpar%zs_max_ipd_temp
        n_interp_points = cpar%zs_n_interp_points
        gauss_degree = cpar%zs_gauss_quad_order
        R_cutoff = cpar%zs_los_cut

        ! Set up interpolation grid for evaluating line of sight b_nu
        if (.not. allocated(temperature_grid)) allocate(temperature_grid(n_interp_points))
        call linspace(min_ipd_temp, max_ipd_temp, temperature_grid)

        call initialize_earth_pos_spline(cpar)

        allocate(obs_time(tod%nscan_tot))
        allocate(obs_pos(3, tod%nscan_tot))

        ! Set up spline objects for observer position (requires knowing the full scan list ahead of time)
        obs_time = 0.d0
        obs_pos = 0.d0
        do i = 1, tod%nscan
            obs_time(tod%scanid(i)) = tod%scans(i)%t0(1)
        end do
        call mpi_allreduce(MPI_IN_PLACE, obs_time, size(obs_time), &
                & MPI_DOUBLE_PRECISION, MPI_SUM, tod%comm, ierr)

        do i = 1, tod%nscan
            obs_pos(:, tod%scanid(i)) = tod%scans(i)%satpos
        end do
        call mpi_allreduce(MPI_IN_PLACE, obs_pos, size(obs_pos), &
                & MPI_DOUBLE_PRECISION, MPI_SUM, tod%comm, ierr)
        
        do i = 1, 3
            call spline_simple(tod%zodi_obs_pos_spl_obj(i), obs_time, obs_pos(i, :))
        end do
        tod%zodi_init_cache_time = tod%scans(1)%t0(1)
        call tod%reset_zodi_cache()
        tod%zodi_min_obs_time = minval(obs_time)
        tod%zodi_max_obs_time = maxval(obs_time)
        
        allocate(tod%zodi_emissivity(zodi%n_comps))
        allocate(tod%zodi_albedo(zodi%n_comps))
        tod%zodi_emissivity(:) = cpar%ds_zodi_emissivity(tod%band, :)
        tod%zodi_albedo(:) = cpar%ds_zodi_albedo(tod%band, :)
        
        if (count(tod%zodi_albedo /= 0.d0) > 0) then
            tod%zodi_scattering = .true.
        else
            tod%zodi_scattering = .false.
        end if
        !allocate spectral quantities
        allocate(tod%zodi_spl_phase_coeffs(tod%ndet, 3))
        allocate(tod%zodi_spl_solar_irradiance(tod%ndet))
        allocate(tod%zodi_phase_func_normalization(tod%ndet))
        allocate(tod%zodi_b_nu_spl_obj(tod%ndet))

        ! Inform `get_zodi_emission` that the relevant tod zodi parameters have been sucessfully initialized
        tod%zodi_tod_params_are_initialized = .true.
    end subroutine initialize_tod_zodi_mod


    subroutine get_zodi_emission(tod, pix, scan, det, s_zodi_scat, s_zodi_therm)
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
        integer(i4b) :: i, j, k, pixel, lookup_idx, n_tod, ierr
        real(dp) :: earth_lon, R_obs, R_max, dt_tod, obs_time
        real(dp) :: unit_vector(3), X_unit_LOS(3, gauss_degree), X_LOS(3, gauss_degree), obs_pos(3), earth_pos(3)
        real(dp), dimension(gauss_degree) :: R_LOS, T_LOS, density_LOS, gauss_nodes, gauss_weights, &
                                            solar_flux_LOS, scattering_angle, phase_function, b_nu_LOS
        class(ZodiComponent), pointer :: comp

        if (.not. tod%zodi_tod_params_are_initialized) then 
            print *, "Error: Zodi parameters have not been initialized."
            print *, "make sure that `initialize_zodi_tod_parameters` is called in the constructor of `tod_your_experiment_mod`"
            stop
        else if (.not. allocated(tod%zodi_b_nu_spl_obj)) then
            print *, "Error: Zodi splines have not been initialized."
            print *, "make sure that `initialize_zodi_splines` is called at least once before this function is ran, and then again everytime the bandpasses are updated."
            stop
        end if

        s_zodi_scat = 0.d0
        s_zodi_therm = 0.d0
        n_tod = size(pix, dim=1)

        dt_tod = (1.d0/tod%samprate) * SECOND_TO_DAY ! dt between two samples in units of days (assumes equispaced tods)
        obs_pos = tod%scans(scan)%satpos
        R_obs = norm2(obs_pos)
        obs_time = tod%scans(scan)%t0(1)
        do i = 1, 3
            earth_pos(i) = splint_simple(earth_pos_spl_obj(i), tod%scans(scan)%t0(1))
        end do        

        earth_lon = atan(earth_pos(2), earth_pos(1))

        do i = 1, n_tod
            ! Reset cache if time between last cache update and current time is larger than `delta_t_reset`.
            ! NOTE: this cache is only effective if the scans a core handles are in chronological order.
            obs_time = obs_time + dt_tod
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
                call tod%reset_zodi_cache(obs_time)
            end if

            lookup_idx = tod%pix2ind(pix(i))
            if (tod%zodi_therm_cache(lookup_idx, 1, det) > 0.d0) then
                s_zodi_scat(i, :) = tod%zodi_scat_cache(lookup_idx, :, det)
                s_zodi_therm(i, :) = tod%zodi_therm_cache(lookup_idx, :, det)
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

            if (tod%zodi_scattering) then
                solar_flux_LOS = tod%zodi_spl_solar_irradiance(det) / R_LOS**2
                call get_scattering_angle(X_LOS, X_unit_LOS, R_LOS, scattering_angle)
                call get_phase_function(scattering_angle, tod%zodi_spl_phase_coeffs(det, :), tod%zodi_phase_func_normalization(det), phase_function)
            end if

            call get_dust_grain_temperature(R_LOS, T_LOS)
            call splint_simple_multi(tod%zodi_b_nu_spl_obj(det), T_LOS, b_nu_LOS)

            comp => comp_list
            k = 1
            do while (associated(comp))
                call comp%get_density(X_LOS, earth_lon, density_LOS)
                if (tod%zodi_scattering) then 
                    s_zodi_scat(i, k) = sum(density_LOS * solar_flux_LOS * phase_function * gauss_weights)
                    tod%zodi_scat_cache(lookup_idx, k, det) = s_zodi_scat(i, k)
                end if
                s_zodi_therm(i, k) = sum(density_LOS * b_nu_LOS * gauss_weights)
                tod%zodi_therm_cache(lookup_idx, k, det) = s_zodi_therm(i, k)
                comp => comp%next()
                k = k + 1
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
        b = 2.d0 * (obs_pos(1) * cos_lat * cos(lon) + obs_pos(2) * cos_lat * sin(lon))
        d = R_obs**2 - R_cutoff**2
        q = -0.5d0 * b * (1.d0 + sqrt(b**2 - (4.d0 * d)) / abs(b))
        R_max = max(q, d / q)
    end subroutine get_R_max

    subroutine get_dust_grain_temperature(R, T_out)
        real(dp), dimension(:), intent(in) :: R
        real(dp), dimension(:), intent(out) :: T_out
        T_out = zodi%T_0 * R ** (-zodi%delta)
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


    subroutine update_zodi_splines(tod, bandpass, det)
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

        real(dp), allocatable :: b_nu(:), integrals(:)
        integer(i4b) :: i, j
    
        allocate(integrals(n_interp_points))
        allocate(b_nu(bandpass%p%n))
        do i = 1, n_interp_points    
            call get_blackbody_emission(bandpass%p%nu, temperature_grid(i), b_nu)
            integrals(i) = tsum(bandpass%p%nu, bandpass%p%tau * b_nu)
        end do
        call spline_simple(tod%zodi_b_nu_spl_obj(det), temperature_grid, integrals, regular=.true.)
        do j = 1, size(tod%zodi_spl_phase_coeffs, dim=2) 
            tod%zodi_spl_phase_coeffs(det, j) = splint_simple(zodi%phase_coeff_spl(j), bandpass%p%nu_c)
        end do

        tod%zodi_spl_solar_irradiance = splint_simple(zodi%solar_irradiance_spl, bandpass%p%nu_c)
        call get_phase_normalization(tod%zodi_spl_phase_coeffs(det, :), tod%zodi_phase_func_normalization(det))
    end subroutine update_zodi_splines

    subroutine accumulate_zodi_emissivities(tod, s_therm, s_scat, res, A_T_A, AY, group_comps)
        ! Returns the A^T A and A Y matrices when solving the normal equations 
        ! (AX = Y, where X is the emissivity vector).
        !
        ! TODO: Add the covariance matrix
        !
        ! Parameters
        ! ----------
        ! tod :
        !     The TOD object holding the component emissivities and albedos.
        ! s_therm : (ntod, ncomps)
        !     The thermal zodiacal emission integral (without being scaled by emissivities).
        ! s_scat : (ntod, ncomps)
        !     The scattered zodiacal light integral (without being scaled by albedos).
        ! res : (ntod)
        !     The residual (data - sky model).
        ! A_T_A : (ncomps, ncomps)
        !     The A^T A matrix.
        ! AY : (ncomps)
        !     The A Y matrix.
        ! group_comps :
        !     Whether to group the components (bands into one group, and feature + ring into another).
        class(comm_tod), intent(in) :: tod
        real(dp), intent(in) :: s_therm(:, :), s_scat(:, :), res(:)
        real(dp), intent(inout) :: A_T_A(:, :), AY(:)
        logical(lgt), intent(in) :: group_comps
        integer(i4b) :: i, j, k, det, ierr, ndet, ntod, ncomp
        real(dp) :: term1, term2, residual
        
        ntod = size(s_therm, dim=1)
        if (group_comps) then
            ncomp = 3
        else 
            ncomp = size(s_therm, dim=2)
        end if
        do i = 1, ntod
            do j = 1, ncomp
                if (group_comps .and. j == 2) then
                    term1 = sum(s_therm(i, 2:4), dim=1) * (1.d0 - tod%zodi_albedo(2))
                else if (group_comps .and. j == 3) then
                    term1 = sum(s_therm(i, 5:6), dim=1) * (1.d0 - tod%zodi_albedo(5))
                else
                    term1 = s_therm(i, j) * (1.d0 - tod%zodi_albedo(j))
                end if
                residual = res(i) - sum(s_scat(i, :) * tod%zodi_albedo(:), dim=1)
                AY(j) = AY(j) + residual * term1
                do k = 1, ncomp
                    if (group_comps .and. k == 2) then
                        term2 = sum(s_therm(i, 2:4), dim=1) * (1.d0 - tod%zodi_albedo(2))
                    else if (group_comps .and. k == 3) then
                        term2 = sum(s_therm(i, 5:6), dim=1) * (1.d0 - tod%zodi_albedo(5))
                    else
                        term2 = s_therm(i, k) * (1.d0 - tod%zodi_albedo(k))
                    end if
                    A_T_A(j, k) = A_T_A(j, k) + term1 * term2
                end do
            end do
        end do
    end subroutine accumulate_zodi_emissivities

    subroutine accumulate_zodi_albedos(tod, s_therm, s_scat, res, A_T_A, AY, group_comps)
        ! Returns the A^T A and A Y matrices when solving the normal equations 
        ! (AX = Y, where X is the albedo vector).
        !
        ! TODO: Add the covariance matrix
        !
        ! Parameters
        ! ----------
        ! tod :
        !     The TOD object holding the component emissivities and albedos.
        ! s_therm : (ntod, ncomps)
        !     The thermal zodiacal emission integral (without being scaled by emissivities).
        ! s_scat : (ntod, ncomps)
        !     The scattered zodiacal light integral (without being scaled by albedos).
        ! res : (ntod)
        !     The residual (data - sky model).
        ! A_T_A : (ncomps, ncomps)
        !     The A^T A matrix.
        ! AY : (ncomps)
        !     The A Y matrix.
        ! group_comps :
        !     Whether to group the components (bands into one group, and feature + ring into another).
        class(comm_tod), intent(in) :: tod
        real(dp), intent(in):: s_therm(:, :), s_scat(:, :), res(:)
        real(dp), intent(inout) :: A_T_A(:, :), AY(:)
        logical(lgt), intent(in) :: group_comps
        integer(i4b) :: i, j, k, det, ierr, ndet, ntod, ncomp
        real(dp) :: term1, term2, residual
        
        ntod = size(s_scat, dim=1)
        if (group_comps) then
            ncomp = 3
        else 
            ncomp = size(s_scat, dim=2)
        end if

        do i = 1, ntod
            do j = 1, ncomp
                if (group_comps .and. j == 2) then
                    term1 = sum(s_scat(i, 2:4), dim=1) - (tod%zodi_emissivity(2) * sum(s_therm(i, 2:4), dim=1))
                else if (group_comps .and. j == 3) then
                    term1 = sum(s_scat(i, 5:6), dim=1) - (tod%zodi_emissivity(5) * sum(s_therm(i, 5:6), dim=1))
                else
                    term1 = s_scat(i, j) - (tod%zodi_emissivity(j) * s_therm(i, j))
                end if
                residual = res(i) - sum(s_therm(i, :) * tod%zodi_emissivity(:), dim=1) 
                AY(j) = AY(j) + residual * term1
                do k = 1, ncomp
                    if (group_comps .and. k == 2) then
                        term2 = sum(s_scat(i, 2:4), dim=1) - (tod%zodi_emissivity(2) * sum(s_therm(i, 2:4), dim=1))
                    else if (group_comps .and. k == 3) then
                        term2 = sum(s_scat(i, 5:6), dim=1) - (tod%zodi_emissivity(5) * sum(s_therm(i, 5:6), dim=1))
                    else
                        term2 = s_scat(i, k) - (tod%zodi_emissivity(k) * s_therm(i, k))
                    end if
                    A_T_A(j, k) = A_T_A(j, k) + term1 * term2
                end do
            end do
        end do
    end subroutine accumulate_zodi_albedos

    subroutine solve_Ax_zodi(A_T_A, AY, handle, X)
        ! Solve the normal equations and return the parameter vector X (albedos or emissivities).
        ! X = (A^T A)^{-1} (A Y)
        !
        ! TODO: Add the covariance matrix
        !
        ! Parameters
        ! ----------
        ! A_T_A:
        !   (A^T A) matrix.
        ! AY:
        !   (A Y) vector.
        ! handle:
        !   The random number generator handle.
        real(dp), dimension(:, :), intent(in):: A_T_A
        real(dp), dimension(:), intent(in) :: AY
        type(planck_rng), intent(inout) :: handle
        real(dp), dimension(:), intent(inout) :: X

        real(dp), allocatable, dimension(:, :) :: A_T_A_inv, A_T_A_inv_sqrt
        real(dp), allocatable, dimension(:) :: eta
        integer(i4b) :: i, ierr

        allocate(A_T_A_inv, mold=A_T_A)
        allocate(A_T_A_inv_sqrt, mold=A_T_A)
        allocate(eta, mold=AY)
        A_T_A_inv = A_T_A

        call invert_matrix(A_T_A_inv)
        call cholesky_decompose(A_T_A_inv, A_T_A_inv_sqrt)
        do i = 1, size(AY)
            eta(i) = rand_gauss(handle)
        end do
        X = matmul(A_T_A_inv, AY) + matmul(A_T_A_inv_sqrt, eta)
    end subroutine solve_Ax_zodi

    subroutine sample_dynamic_zodi_parameters(tod)
        class(comm_tod), intent(inout) :: tod
        integer(i4b) :: i, ierr
    end subroutine sample_dynamic_zodi_parameters

    subroutine downsample_zodi_res_and_pointing(tod, scan_id, s_tod, s_tot, s_zodi, pix, flag, mask)
        ! Downsamples the residual and pointing to allow for fast zodi parameter estimation.
        ! The downsampled arrays are stored in the detscan object (e.g tod%scans(scan_id)%d(j)%downsamp_res)
        !
        ! Parameters
        ! ----------
        ! tod:
        !     The TOD object.
        ! scan_id:
        !     The current scan to process.
        ! s_tod:
        !     The TOD signal.
        ! s_tot:
        !     The total sky signal.
        ! s_zodi:
        !     The zodi signal.
        ! pix:
        !     The pixel indices.
        ! flag:
        !     Flagged tods.
        ! mask:
        !     Galaxy mask.
        !
        class(comm_tod), intent(inout) :: tod
        integer(i4b), intent(in):: scan_id, pix(:, :, :), flag(:, :)
        real(sp), intent(in) :: s_tod(:, :), s_tot(:, :), s_zodi(:, :), mask(:)

        real(dp) :: box_width
        real(sp), allocatable, dimension(:) :: res, res_truncated
        integer(i4b), allocatable, dimension(:) :: pix_truncated
        logical(lgt), allocatable, dimension(:) :: indices
        integer(i4b) :: i, j, ext(2), ndet, ntod, box_halfwidth, upper_bound

        ndet = tod%ndet
        ntod = size(s_tot, dim=1)

        ! Get box car size and make sure it has odd numbers so that we can pick out a center pixel.
        box_width = tod%samprate / tod%samprate_lowres
        if (box_width < 1.) stop "Cannot downsample zodi tods if box car width is less than 1 sample!"
        box_halfwidth = box_width / 2
        if (mod(box_halfwidth, 2) == 0) box_width = box_width + 1.

        allocate(res(ntod), indices(ntod))
        do j = 1, ndet
            ! Add zodi back to the residual
            res = s_tod(:, j) - s_tot(:,j) + s_zodi(:, j)

            ! Get indices of flagged data + masked pixels and remove these from the fitted values
            indices = .true.
            where ((iand(flag(:, j), tod%flag0) .ne. 0.) .or. (mask .eq. 0.)) indices = .false.
            res_truncated = pack(res, indices)
            pix_truncated = pack(pix(:, j, 1), indices)

            ! Get downsampled shape (ext), and allocate the downsampled arrays
            call tod%downsample_tod(res_truncated, ext, step=box_width)
            allocate(tod%scans(scan_id)%d(j)%downsamp_res(ext(1):ext(2)))
            allocate(tod%scans(scan_id)%d(j)%downsamp_pointing(ext(1):ext(2)))
            tod%scans(scan_id)%d(j)%downsamp_res = 0.

            ! Downsample the residual
            call tod%downsample_tod(res_truncated, ext, tod%scans(scan_id)%d(j)%downsamp_res, step=box_width)

            ! Construct the downsampled pointing array (pick central pixel in each box)
            ! This is a bit cluncky, but we have essentially copied the code in downsample_tod and 
            ! picked out the center pixel
            tod%scans(scan_id)%d(j)%downsamp_pointing = 0
            do i = 1, int(size(res_truncated) / box_width)
                tod%scans(scan_id)%d(j)%downsamp_pointing(i) = pix_truncated(floor(i * box_width))
            end do

            ! Remove padded data (including the first and last sample as these might have boundary effects)
            upper_bound = ubound(tod%scans(scan_id)%d(j)%downsamp_res, dim=1) + ext(1) - 1 
            tod%scans(scan_id)%d(j)%downsamp_pointing = tod%scans(scan_id)%d(j)%downsamp_pointing(1:upper_bound)
            tod%scans(scan_id)%d(j)%downsamp_res = tod%scans(scan_id)%d(j)%downsamp_res(1:upper_bound)

        end do
    end subroutine

    subroutine sample_zodi_params(tod, handle, group_comps)
        class(comm_tod), intent(inout) :: tod
        type(planck_rng), intent(inout) :: handle
        logical(lgt), intent(in) :: group_comps
        real(sp), allocatable, dimension(:, :) :: s_scat, s_therm
        integer(i4b) :: i, j, ierr, scan, ndet, nhorn, nscan, ntod, n_comps_to_fit
        real(dp), allocatable, dimension(:, :) :: A_T_A_emiss, A_T_A_albedo, A_T_A_emiss_reduced, A_T_A_albedo_reduced
        real(dp), allocatable, dimension(:) :: AY_emiss, AY_albedo, AY_emiss_reduced, AY_albedo_reduced, emissivities, albedos

        ndet = tod%ndet
        nscan = tod%nscan

        ! Get the number of components to fit (depends on if we want to sample the k98 groups or all components individually)
        ! and initialize Ax = Y matrices
        n_comps_to_fit = zodi%n_comps
        if (group_comps) n_comps_to_fit = 3
        allocate(A_T_A_emiss(n_comps_to_fit, n_comps_to_fit), A_T_A_albedo(n_comps_to_fit, n_comps_to_fit))
        allocate(A_T_A_emiss_reduced(n_comps_to_fit, n_comps_to_fit), A_T_A_albedo_reduced(n_comps_to_fit, n_comps_to_fit))
        allocate(AY_emiss(n_comps_to_fit), AY_albedo(n_comps_to_fit), AY_emiss_reduced(n_comps_to_fit), AY_albedo_reduced(n_comps_to_fit))
        allocate(emissivities(n_comps_to_fit), albedos((n_comps_to_fit)))
        A_T_A_emiss = 0.
        A_T_A_albedo = 0.
        A_T_A_albedo_reduced = 0.
        A_T_A_emiss_reduced = 0.
        AY_emiss = 0.
        AY_albedo = 0.
        AY_emiss_reduced = 0.
        AY_albedo_reduced = 0.

        ! Loop over downsampled data, and evaluate emissivity
        do scan = 1, nscan
            do j = 1, tod%ndet
                ntod = size(tod%scans(scan)%d(j)%downsamp_res)
                allocate(s_scat(ntod, zodi%n_comps), s_therm(ntod, zodi%n_comps))
                call get_zodi_emission(tod, tod%scans(scan)%d(j)%downsamp_pointing, scan, j, s_scat, s_therm)
                call accumulate_zodi_emissivities(tod, real(s_therm, dp), real(s_scat, dp), real(tod%scans(scan)%d(j)%downsamp_res, dp), A_T_A_emiss, AY_emiss, group_comps)
                deallocate(s_scat, s_therm)
            end do
        end do
        call mpi_reduce(A_T_A_emiss, A_T_A_emiss_reduced, size(A_T_A_emiss), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(AY_emiss, AY_emiss_reduced, size(AY_emiss), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (tod%myid == 0) then
            call solve_Ax_zodi(A_T_A_emiss_reduced, AY_emiss_reduced, handle, emissivities)
            ! Prior on emissivity (eps > 0)
            ! where (X < 0)
            !    X = 0
            ! endwhere
            if (group_comps) then
               tod%zodi_emissivity(1) = emissivities(1)
               tod%zodi_emissivity(2:4) = emissivities(2)
               tod%zodi_emissivity(5:6) = emissivities(3)
            else 
               tod%zodi_emissivity = emissivities
            end if
            print *, "Sampled emissivity: ", emissivities
        end if
        ! Loop over downsampled data, and evaluate albedo
        do scan = 1, nscan
            do j = 1, tod%ndet
                ntod = size(tod%scans(scan)%d(j)%downsamp_res)
                allocate(s_scat(ntod, zodi%n_comps), s_therm(ntod, zodi%n_comps))
                call get_zodi_emission(tod, tod%scans(scan)%d(j)%downsamp_pointing, scan, j, s_scat, s_therm)
                call accumulate_zodi_albedos(tod, real(s_therm, dp), real(s_scat, dp), real(tod%scans(scan)%d(j)%downsamp_res, dp), A_T_A_albedo, AY_albedo, group_comps)
                deallocate(s_scat, s_therm)
            end do
        end do
        call mpi_reduce(A_T_A_albedo, A_T_A_albedo_reduced, size(A_T_A_albedo), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(AY_albedo, AY_albedo_reduced, size(AY_albedo), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (tod%myid == 0) then
            call solve_Ax_zodi(A_T_A_albedo_reduced, AY_albedo_reduced, handle, albedos)
            if (group_comps) then
               tod%zodi_albedo(1) = albedos(1)
               tod%zodi_albedo(2:4) = albedos(2)
               tod%zodi_albedo(5:6) = albedos(3)
            else 
               tod%zodi_albedo = albedos
            end if
            print *, "Sampled albedo: ", albedos
        end if
    end subroutine

end module comm_tod_zodi_mod