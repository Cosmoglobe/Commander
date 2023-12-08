module comm_tod_zodi_mod
   use comm_param_mod
   use comm_utils
   use comm_tod_mod
   use comm_zodi_mod
   use spline_1D_mod
   use comm_bp_mod
   use comm_hdf_mod
   implicit none

   private
   public initialize_tod_zodi_mod, get_zodi_emission, update_zodi_splines, output_tod_params_to_hd5, read_tod_zodi_params, get_instantaneous_zodi_emission, get_interp_zodi_emission

   type :: ZodiCompLOS
      real(dp) :: R_min, R_max
      real(dp), allocatable, dimension(:) :: gauss_nodes, gauss_weights
      real(dp), allocatable, dimension(:) :: R, T, n, F_sol, Theta, Phi, B_nu
      real(dp), allocatable, dimension(:, :) :: X, X_unit
   end type

   real(dp) :: R_MIN = 3.d-14, R_CUTOFF = 5.2, EPS = TINY(1.0_dp), delta_t_reset
   real(dp), allocatable :: T_grid(:), B_nu_integrals(:)
   type(ZodiCompLOS), allocatable, dimension(:) :: comp_LOS

contains
   subroutine initialize_tod_zodi_mod(cpar)
      ! Initializes the zodiacal emission module.
      !
      ! Parameters
      ! ----------
      ! cpar : class(comm_param)
      !     The parameter object.

      type(comm_params), intent(in) :: cpar
      real(dp) :: min_temp = 40.0, max_temp = 550.0
      integer(i4b) :: n_interp_points = 100
      integer(i4b) :: i, gauss_degree

      allocate(comp_LOS(zodi_model%n_comps))
      allocate(B_nu_integrals(n_interp_points))

      ! hyper parameters
      delta_t_reset = cpar%zs_delta_t_reset

      ! Set up Gauss-Legendre quadrature
      do i = 1, zodi_model%n_comps
         gauss_degree = cpar%zs_los_steps(i)
         allocate (comp_LOS(i)%gauss_nodes(gauss_degree), comp_LOS(i)%gauss_weights(gauss_degree))
         if (cpar%zs_r_min(i) == 0.) then
            comp_LOS(i)%R_min = R_MIN
         else
            comp_LOS(i)%R_min = cpar%zs_r_min(i)
         end if
         comp_LOS(i)%R_max = cpar%zs_r_max(i)
         call leggaus(gauss_degree, comp_LOS(i)%gauss_nodes, comp_LOS(i)%gauss_weights)
         allocate(comp_LOS(i)%X_unit(3, gauss_degree))
         allocate(comp_LOS(i)%X(3, gauss_degree))
         allocate(comp_LOS(i)%R(gauss_degree))
         allocate(comp_LOS(i)%T(gauss_degree))
         allocate(comp_LOS(i)%n(gauss_degree))
         allocate(comp_LOS(i)%F_sol(gauss_degree))
         allocate(comp_LOS(i)%Theta(gauss_degree))
         allocate(comp_LOS(i)%Phi(gauss_degree))
         allocate(comp_LOS(i)%B_nu(gauss_degree))
      end do

      ! Set up interpolation grid for evaluating temperature
      allocate (T_grid(n_interp_points))
      call linspace(min_temp, max_temp, T_grid)

   end subroutine initialize_tod_zodi_mod

   subroutine get_zodi_emission(tod, pix, scan, det, s_zodi_scat, s_zodi_therm, model, always_scattering, use_lowres_pointing, comp)
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
      ! model : type(ZodiModel)
      !     The zodiacal emission model.
      ! always_scattering : logical(lgt), optional
      !     If present, this overrides the default behavior of only including scattering when the albedo is non-zero.
      ! use_lowres_pointing : logical(lgt), optional
      !     If present, the input pixels are converted to low resolution pixels before evaluating the zodiacal emission.
      ! comp : integer(i4b), optional
      !     If present, only evaluate the zodiacal emission for this component.
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
      type(ZodiModel), intent(in) :: model
      logical(lgt), intent(in), optional :: always_scattering, use_lowres_pointing
      integer(i4b), intent(in), optional :: comp
      integer(i4b) :: i, j, k, l, pix_at_zodi_nside, lookup_idx, n_tod, ierr, cache_hits
      logical(lgt) :: scattering, use_lowres
      real(dp) :: earth_lon, R_obs, R_min, R_max, dt_tod, obs_time, phase_normalization, C0, C1, C2
      real(dp) :: unit_vector(3), obs_pos(3), earth_pos(3)
      !real(dp), dimension(gauss_degree) :: R_LOS, T_LOS, density_LOS, solar_flux_LOS, scattering_angle, phase_function, b_nu_LOS
      
      s_zodi_scat = 0.
      s_zodi_therm = 0.
      n_tod = size(pix, dim=1)

      dt_tod = (1./tod%samprate)*SECOND_TO_DAY ! dt between two samples in units of days (assumes equispaced tods)
      obs_time = tod%scans(scan)%t0(1)
      do i = 1, 3
         earth_pos(i) = splint_simple(model%earth_pos_interpolator(i), obs_time)
      end do
      obs_pos = earth_pos!tod%scans(scan)%x0_obs
      R_obs = norm2(obs_pos)
      earth_lon = atan(earth_pos(2), earth_pos(1))

      C0 = zodi_model%C0(tod%band)
      C1 = zodi_model%C1(tod%band)
      C2 = zodi_model%C2(tod%band)
      phase_normalization = get_phase_normalization(C0, C1, C2)
      if (present(always_scattering)) then
         scattering = always_scattering
      else
         scattering = any(tod%zodi_albedo > EPS)
      end if
      ! select the correct cache
      if (present(use_lowres_pointing)) then
         if (tod%nside == zodi_nside) then
            use_lowres = .false.
         else
            if (.not. allocated(tod%zodi_therm_cache_lowres)) stop "zodi cache not allocated. `use_lowres_pointing` should only be true when sampling zodi."
            if (.not. allocated(tod%scans(scan)%downsamp_obs_time)) then
               print *, tod%band, scan, "lowres obs_time not allocated"
               stop
            end if
            use_lowres = .true.
         end if
      else
         use_lowres = .false.
      end if

      cache_hits = 0
      do i = 1, n_tod
         ! Reset cache if time between last cache update and current time is larger than `delta_t_reset`.
         ! NOTE: this cache is only effective if the scans a core handles are in chronological order.
         if (use_lowres) then
            obs_time = tod%scans(scan)%downsamp_obs_time(i)
         else
            obs_time = obs_time + dt_tod ! assumes a time continuous TOD
         end if
         if ((obs_time - tod%zodi_cache_time) >= delta_t_reset) then
            do j = 1, 3
               earth_pos(j) = splint_simple(model%earth_pos_interpolator(j), obs_time)
               obs_pos(j) = earth_pos(j) !splint_simple(tod%x_obs_spline(j), obs_time)
            end do
            R_obs = norm2(obs_pos)
            earth_lon = atan(earth_pos(2), earth_pos(1))
            call tod%clear_zodi_cache(obs_time)
         end if

         ! Get lookup index for cache. If the pixel is already cached, used that value.
         if (use_lowres) then
            lookup_idx = tod%pix2ind_lowres(tod%udgrade_pix_zodi(pix(i)))
            if (tod%zodi_therm_cache_lowres(lookup_idx, 1, det) > 0.d0) then
               if (scattering) s_zodi_scat(i, :) = tod%zodi_scat_cache_lowres(lookup_idx, :, det)
               s_zodi_therm(i, :) = tod%zodi_therm_cache_lowres(lookup_idx, :, det)
               cache_hits = cache_hits + 1
               cycle
            end if
            unit_vector = tod%ind2vec_ecl_lowres(:, lookup_idx)
         else
            lookup_idx = tod%pix2ind(pix(i))
            if (tod%zodi_therm_cache(lookup_idx, 1, det) > 0.d0) then
               if (scattering) s_zodi_scat(i, :) = tod%zodi_scat_cache(lookup_idx, :, det)
               s_zodi_therm(i, :) = tod%zodi_therm_cache(lookup_idx, :, det)
               cache_hits = cache_hits + 1
               cycle
            end if
            unit_vector = tod%ind2vec_ecl(:, lookup_idx)
         end if

         do k = 1, model%n_comps
            ! If comp is present we only evaluate the zodi emission for that component.
            ! If comp == 0 then we evaluate the zodi emission for all components.
            if (present(comp)) then
               if (k /= comp .and. comp /= 0) cycle
            end if

            ! Get line of sight integration range
            call get_sphere_intersection(unit_vector, obs_pos, R_obs, comp_LOS(k)%R_min, R_min)
            call get_sphere_intersection(unit_vector, obs_pos, R_obs, comp_LOS(k)%R_max, R_max)

            do l = 1, 3
               ! Convert quadrature range from [-1, 1] to [R_min, R_max]
               comp_LOS(k)%X_unit(l, :) = (0.5 * (R_max - R_MIN)) * comp_LOS(k)%gauss_nodes + (0.5 * (R_max + R_MIN))
               comp_LOS(k)%X_unit(l, :) = comp_LOS(k)%X_unit(l, :) * unit_vector(l)
               comp_LOS(k)%X(l, :) = comp_LOS(k)%X_unit(l, :) + obs_pos(l)
            end do
            comp_LOS(k)%R = norm2(comp_LOS(k)%X, dim=1)

            if (scattering) then
               comp_LOS(k)%F_sol = model%F_sun(tod%band)/comp_LOS(k)%R**2
               call get_scattering_angle(comp_LOS(k)%X, comp_LOS(k)%X_unit, comp_LOS(k)%R, comp_LOS(k)%Theta)
               call get_phase_function(comp_LOS(k)%Theta, C0, C1, C2, phase_normalization, comp_LOS(k)%Phi)
            end if

            call get_dust_grain_temperature(comp_LOS(k)%R, comp_LOS(k)%T, model%T_0, model%delta)
            call splint_simple_multi(tod%zodi_b_nu_spl_obj(det), comp_LOS(k)%T, comp_LOS(k)%B_nu)

            call model%comps(k)%c%get_density(comp_LOS(k)%X, earth_lon, comp_LOS(k)%n)
            if (scattering) then
               s_zodi_scat(i, k) = sum(comp_LOS(k)%n*comp_LOS(k)%F_sol*comp_LOS(k)%Phi*comp_LOS(k)%gauss_weights) * 0.5*(R_max - R_MIN) * 1d20
               if (use_lowres) then
                  tod%zodi_scat_cache_lowres(lookup_idx, k, det) = s_zodi_scat(i, k)
               else
                  tod%zodi_scat_cache(lookup_idx, k, det) = s_zodi_scat(i, k)
               end if
            end if
            s_zodi_therm(i, k) = sum(comp_LOS(k)%n*comp_LOS(k)%B_nu*comp_LOS(k)%gauss_weights) * 0.5 * (R_max - R_MIN) * 1d20
            if (use_lowres) then
               tod%zodi_therm_cache_lowres(lookup_idx, k, det) = s_zodi_therm(i, k)
            else
               tod%zodi_therm_cache(lookup_idx, k, det) = s_zodi_therm(i, k)
            end if
         end do
      end do
   end subroutine get_zodi_emission


   subroutine get_instantaneous_zodi_emission(tod, pix, obs_time, obs_pos, det, s_zodi, model)
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
      ! model : type(ZodiModel)
      !     The zodiacal emission model.
      ! always_scattering : logical(lgt), optional
      !     If present, this overrides the default behavior of only including scattering when the albedo is non-zero.
      ! use_lowres_pointing : logical(lgt), optional
      !     If present, the input pixels are converted to low resolution pixels before evaluating the zodiacal emission.
      ! comp : integer(i4b), optional
      !     If present, only evaluate the zodiacal emission for this component.
      !
      ! Returns
      ! -------
      ! s_zodi_scat : real(sp), dimension(ntod, ncomps, ndet)
      !     Contribution from scattered sunlight light.
      ! s_zodi_therm : real(sp), dimension(ntod, ncomps, ndet)
      !     Contribution from thermal interplanetary dust emission.
      !
      class(comm_tod), intent(inout) :: tod
      integer(i4b), intent(in) :: pix(:), det
      real(dp) :: obs_time, obs_pos(3)
      real(sp), dimension(:, :), intent(inout) :: s_zodi
      type(ZodiModel), intent(in) :: model

      integer(i4b) :: i, j, k, l, pix_at_zodi_nside, lookup_idx, n_tod, ierr, cache_hits
      logical(lgt) :: scattering, use_lowres
      real(dp) :: earth_lon, R_obs, R_min, R_max, dt_tod, phase_normalization, C0, C1, C2
      real(dp) :: unit_vector(3), earth_pos(3), rotation_matrix(3, 3)
      !real(dp), dimension(gauss_degree) :: R_LOS, T_LOS, density_LOS, solar_flux_LOS, scattering_angle, phase_function, b_nu_LOS

      s_zodi = 0.
      n_tod = size(pix)

      R_obs = norm2(obs_pos)
      obs_pos = obs_pos
      do i = 1, 3
         earth_pos(i) = splint_simple(model%earth_pos_interpolator(i), obs_time)
      end do
      earth_lon = atan(earth_pos(2), earth_pos(1))

      call ecl_to_gal_rot_mat(rotation_matrix)

      C0 = zodi_model%C0(tod%band)
      C1 = zodi_model%C1(tod%band)
      C2 = zodi_model%C2(tod%band)
      phase_normalization = get_phase_normalization(C0, C1, C2)

      scattering = any(tod%zodi_albedo > EPS)

      R_obs = norm2(obs_pos)
      earth_lon = atan(earth_pos(2), earth_pos(1))
      
      do i = 1, n_tod
         ! Reset cache if time between last cache update and current time is larger than `delta_t_reset`.
         
         ! Get lookup index for cache. If the pixel is already cached, used that value.
         call pix2vec_ring(zodi_nside, pix(i), unit_vector)
         unit_vector = matmul(unit_vector, rotation_matrix)
         ! print*, unit_vector
         do k = 1, model%n_comps
            ! Get line of sight integration range
            call get_sphere_intersection(unit_vector, obs_pos, R_obs, comp_LOS(k)%R_min, R_min)
            call get_sphere_intersection(unit_vector, obs_pos, R_obs, comp_LOS(k)%R_max, R_max)

            do l = 1, 3
               ! Convert quadrature range from [-1, 1] to [R_min, R_max]
               comp_LOS(k)%X_unit(l, :) = (0.5 * (R_max - R_MIN)) * comp_LOS(k)%gauss_nodes + (0.5 * (R_max + R_MIN))
               comp_LOS(k)%X_unit(l, :) = comp_LOS(k)%X_unit(l, :) * unit_vector(l)
               comp_LOS(k)%X(l, :) = comp_LOS(k)%X_unit(l, :) + obs_pos(l)
            end do
            comp_LOS(k)%R = norm2(comp_LOS(k)%X, dim=1)
            if (scattering) then
               comp_LOS(k)%F_sol = model%F_sun(tod%band)/comp_LOS(k)%R**2
               call get_scattering_angle(comp_LOS(k)%X, comp_LOS(k)%X_unit, comp_LOS(k)%R, comp_LOS(k)%Theta)
               call get_phase_function(comp_LOS(k)%Theta, C0, C1, C2, phase_normalization, comp_LOS(k)%Phi)
            end if

            call get_dust_grain_temperature(comp_LOS(k)%R, comp_LOS(k)%T, model%T_0, model%delta)
            call splint_simple_multi(tod%zodi_b_nu_spl_obj(det), comp_LOS(k)%T, comp_LOS(k)%B_nu)

            call model%comps(k)%c%get_density(comp_LOS(k)%X, earth_lon, comp_LOS(k)%n)
            if (scattering) then
               s_zodi(i, k) = s_zodi(i, k) + tod%zodi_albedo(k) * sum(comp_LOS(k)%n*comp_LOS(k)%F_sol*comp_LOS(k)%Phi*comp_LOS(k)%gauss_weights) * 0.5*(R_max - R_MIN) * 1d20
            end if
            s_zodi(i, k) = s_zodi(i, k) +  (1. - tod%zodi_albedo(k)) * tod%zodi_emissivity(k) * sum(comp_LOS(k)%n*comp_LOS(k)%B_nu*comp_LOS(k)%gauss_weights) * 0.5 * (R_max - R_MIN) * 1d20
         end do
      end do
   end subroutine get_instantaneous_zodi_emission


   subroutine get_interp_zodi_emission(tod, pix, scan, s_zodi)
      class(comm_tod), intent(inout) :: tod
      integer(i4b), intent(in) :: pix(:)
      integer(i4b), intent(in) :: scan
      real(sp), dimension(:), intent(inout) :: s_zodi

      integer(i4b) :: i, nobs, k, ncoeffs
      real(dp), allocatable :: obs_time_samp(:), freqs(:)
      real(dp) :: dt_tod, obs_time, period, dt, tmp

      
      nobs = size(tod%zodi_fourier_cube, dim=1)
      ncoeffs = nobs/2 + 1

      period = 365.242199

      obs_time_samp = [(i * (period/real(nobs)), i = 1, nobs)]
      obs_time_samp = obs_time_samp +  47871.0100490336
      dt = obs_time_samp(2) - obs_time_samp(1)

      allocate(freqs(0:nobs-1))

      freqs = 0.
      do k = 1, ncoeffs
         if (k < nobs/2.) then
            freqs(k) = k / (dt * nobs)
         else
            freqs(k) = (k - nobs) / (dt * nobs)
         end if
         if (k < ncoeffs - 1) freqs(nobs - k) = - freqs(k)
      end do

      dt_tod = (1./tod%samprate)*SECOND_TO_DAY ! dt between two samples in units of days (assumes equispaced tods)
      obs_time = tod%scans(scan)%t0(1)
      
      ! print *,"pix bounds:", lbound(pix), ubound(pix)
      ! print *,"pix min max:", minval(pix), maxval(pix)
      ! print *,"udgrade_pix_zodi bounds:", lbound(tod%udgrade_pix_zodi), ubound(tod%udgrade_pix_zodi)
      ! print *,"udgrade_pix_zodi min max:", minval(tod%udgrade_pix_zodi), maxval(tod%udgrade_pix_zodi)
      ! stop
      do i = 1, size(pix)
         ! print *, "coeffs:", tod%zodi_fourier_cube(:, pix(i))
         ! print *, "freqs:", freqs
         ! print *, "obs_time:", obs_time
         ! print *, "obs_time_samp(1):", obs_time_samp(1)
         ! stop
         call fft_interp_zodi(tod%zodi_fourier_cube(:, tod%udgrade_pix_zodi(pix(i))), freqs, obs_time, obs_time_samp(1), tmp)
         s_zodi(i) = real(tmp, sp)
         obs_time = obs_time + dt_tod
      end do

   end subroutine get_interp_zodi_emission

   subroutine fft_interp_zodi(fft_coeffs, fft_freqs, t, t_0, s_zodi_pix)
      complex(spc), intent(in) :: fft_coeffs(:)
      real(dp), intent(in) :: fft_freqs(:), t, t_0
      real(dp), intent(out) :: s_zodi_pix

      complex(spc) :: i = cmplx(0., 1.)
      integer(i4b) :: k

      s_zodi_pix = 0
      do k = 1, size(fft_coeffs)
         s_zodi_pix = s_zodi_pix + fft_coeffs(k) * exp(i * 2. * pi  * fft_freqs(k) * (t-t_0))
         ! s_zodi_pix = s_zodi_pix + exp(i * 2. * pi  * fft_freqs(k) * mod(t, t_0))
      end do
      s_zodi_pix = s_zodi_pix / size(fft_coeffs)
   end subroutine fft_interp_zodi

   ! subroutine fft_interp_zodi(fft_coeffs, fft_freqs, t, t_0, s_zodi_pix)
   !    complex(spc), intent(in) :: fft_coeffs(:)
   !    real(dp), intent(in) :: fft_freqs(:), t, t_0
   !    real(sp), intent(out) :: s_zodi_pix

   !    complex(spc) :: s_zodi_pix_complex
   !    complex(spc) :: i = cmplx(0., 1.)
   !    integer(i4b) :: k

   !    s_zodi_pix_complex = cmplx(0., 0.)
   !    do k = 1, size(fft_coeffs)
   !       ! print *, "s_zodi_pix_complex", s_zodi_pix_complex
   !       ! print *, "fft_coeffs(k)", fft_coeffs(k)
   !       ! print *, "fft_freqs(k)", fft_freqs(k)
   !       ! print *, "mod(t, t_0)", mod(t, t_0)
   !       ! print *, "i", i
   !       ! print *, "exp", exp(i * 2. * pi  * fft_freqs(k) * mod(t, t_0))

   !       ! s_zodi_pix_complex = s_zodi_pix_complex + fft_coeffs(k) * exp(i * 2. * pi  * fft_freqs(k) * mod(t, t_0))
   !       s_zodi_pix_complex = s_zodi_pix_complex + fft_coeffs(k)
   !    end do
   !    s_zodi_pix = real(s_zodi_pix_complex, sp) / size(fft_coeffs)
   ! end subroutine fft_interp_zodi

   ! Functions for evaluating the zodiacal emission
   ! -----------------------------------------------------------------------------------
   subroutine get_sphere_intersection(unit_vector, obs_pos, R_obs, R_cutoff, R_intersection)
      ! Computes R_max (the length of the LOS such that it stops exactly at los_cutoff_radius).

      real(dp), intent(in), dimension(:) :: unit_vector, obs_pos
      real(dp), intent(in) :: R_obs, R_cutoff
      real(dp), intent(out) :: R_intersection
      real(dp) :: lon, lat, cos_lat, b, d, q

      if (R_obs > R_cutoff) then
         R_intersection = EPS
         return
      end if

      lon = atan(unit_vector(2), unit_vector(1))
      lat = asin(unit_vector(3))
      cos_lat = cos(lat)
      b = 2.*(obs_pos(1)*cos_lat*cos(lon) + obs_pos(2)*cos_lat*sin(lon))
      d = R_obs**2 - R_cutoff**2
      q = -0.5*b*(1.+sqrt(b**2 - (4.*d))/abs(b))
      R_intersection = max(q, d/q)
   end subroutine get_sphere_intersection

   subroutine get_dust_grain_temperature(R, T_out, T_0, delta)
      real(dp), dimension(:), intent(in) :: R
      real(dp), dimension(:), intent(out) :: T_out
      real(dp), intent(in) :: T_0, delta
      T_out = T_0*R**(-delta)
   end subroutine get_dust_grain_temperature

   subroutine get_blackbody_emission(nus, T, b_nu)
      real(dp), intent(in) :: nus(:), T
      real(dp), dimension(:), intent(out) :: b_nu
      b_nu = ((2.*h*nus**3)/(c*c))/(exp((h*nus)/(k_B*T)) - 1.)
   end subroutine get_blackbody_emission

   subroutine get_scattering_angle(X_helio_vec_LOS, X_vec_LOS, R_helio_LOS, scattering_angle)
      real(dp), intent(in) :: X_helio_vec_LOS(:, :), X_vec_LOS(:, :), R_helio_LOS(:)
      real(dp), dimension(:), intent(out) :: scattering_angle
      real(dp), allocatable, dimension(:) :: cos_theta, R_LOS

      allocate(cos_theta(size(X_vec_LOS, dim=1)))
      allocate(R_LOS(size(X_vec_LOS, dim=1)))


      R_LOS = norm2(X_vec_LOS, dim=1)
      cos_theta = sum(X_helio_vec_LOS*X_vec_LOS, dim=1)/(R_LOS*R_helio_LOS)
      ! clip cos(theta) to [-1, 1]
      where (cos_theta > 1)
         cos_theta = 1
      elsewhere(cos_theta < -1)
         cos_theta = -1
      end where
      scattering_angle = acos(-cos_theta)
   end subroutine get_scattering_angle

   subroutine get_phase_function(Theta, C0 , C1 , C2, N, phase_function)
      real(dp), intent(in) :: Theta(:), C0, C1, C2, N
      real(dp), intent(out) :: phase_function(:)
      phase_function = N *  (C0 + (C1 * Theta) + exp(C2 * Theta))
   end subroutine

   function get_phase_normalization(C0, C1, C2) result(N)
      real(dp), intent(in) :: C0, C1, C2
      real(dp) :: term1, term2, term3, term4, N

      term1 = 2.*pi
      term2 = 2.*C0
      term3 = pi*C1
      term4 = (exp(C2 * pi) + 1.)/(C2**2 + 1.)
      N = 1. / (term1 * (term2 + term3 + term4))
   end function

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

      class(comm_tod), intent(inout) :: tod
      class(comm_bp_ptr), intent(in) :: bandpass
      integer(i4b), intent(in) :: det
      type(ZodiModel), intent(inout) :: model

      real(dp), allocatable :: b_nu(:)
      integer(i4b) :: i, j

      allocate (b_nu(bandpass%p%n))
      do i = 1, size(B_nu_integrals)
         call get_blackbody_emission(bandpass%p%nu, T_grid(i), b_nu)
         B_nu_integrals(i) = tsum(bandpass%p%nu, bandpass%p%tau*b_nu)
      end do
      call spline_simple(tod%zodi_b_nu_spl_obj(det), T_grid, B_nu_integrals, regular=.true.)
   end subroutine update_zodi_splines

   subroutine output_tod_params_to_hd5(cpar, model, tod, iter)
      ! Writes the zodi model to an hdf file
      type(comm_params), intent(in) :: cpar
      class(comm_tod), intent(inout) :: tod
      type(ZodiModel), intent(in) :: model
      integer(i4b), intent(in) :: iter

      integer(i4b) :: i, j, hdferr, ierr, unit
      logical(lgt) :: exist, init, new_header
      character(len=6) :: itext
      character(len=4) :: ctext
      character(len=512) :: zodi_path, tod_path, band_path, det_path, comp_path, chainfile, hdfpath
      character(len=10), allocatable :: param_names(:)
      real(dp), allocatable :: param_values(:)
      type(hdf_file) :: file
      type(h5o_info_t) :: object_info

      if (.not. cpar%myid_chain == 0) return

      call int2string(cpar%mychain, ctext)
      chainfile = trim(adjustl(cpar%outdir))//'/chain'// &
      & '_c'//trim(adjustl(ctext))//'.h5'

      inquire (file=trim(chainfile), exist=exist)
      call open_hdf_file(chainfile, file, 'b')

      call int2string(iter, itext)
      zodi_path = trim(adjustl(itext))//'/zodi'
      tod_path = trim(adjustl(zodi_path))//'/tod'

      call create_hdf_group(file, trim(adjustl(tod_path)))
      band_path = trim(adjustl(tod_path))//'/'//trim(adjustl(tod%freq))
      call create_hdf_group(file, trim(adjustl(band_path)))
      do i = 1, model%n_comps
         comp_path = trim(adjustl(band_path))//'/'//trim(adjustl(model%comp_labels(i)))//'/'
         call create_hdf_group(file, trim(adjustl(comp_path)))
         call write_hdf(file, trim(adjustl(comp_path))//'/emissivity', tod%zodi_emissivity(i))
         call write_hdf(file, trim(adjustl(comp_path))//'/albedo', tod%zodi_albedo(i))
      end do
      call close_hdf_file(file)
   end subroutine

   subroutine read_tod_zodi_params(cpar, model, tod)
      type(comm_params), intent(in) :: cpar
      type(ZodiModel), intent(in) :: model
      class(comm_tod), intent(inout) :: tod

      logical(lgt) :: exist
      integer(i4b) :: i, l, comp, unit, ierr, initsamp
      character(len=6) :: itext, itext2

      type(hdf_file) :: file
      real(dp) :: lambda, lambda_min, lambda_max
      real(dp), allocatable :: emissivity(:), albedo(:)
      character(len=512) :: chainfile, emissivity_path, albedo_path, band_path, comp_path, tod_path, group_name

      allocate(tod%zodi_emissivity(tod%zodi_n_comps))
      allocate(tod%zodi_albedo(tod%zodi_n_comps))
      lambda_min = 0.1
      lambda_max = 4.
      do i = 1, cpar%zs_ncomps
         if (trim(adjustl(cpar%zs_init_hdf(i))) == 'none') then
            tod%zodi_emissivity(i) = 1.
            lambda = (c / tod%nu_c(1)) * 1e6 ! in microns
            if ((lambda_min < lambda) .and. (lambda_max > lambda)) then
               tod%zodi_albedo(i) = 0.5
            else
               tod%zodi_albedo(i) = 0.
            end if
            cycle
         end if

         if (cpar%myid == cpar%root) then
            unit = getlun()
            call get_chainfile_and_samp(trim(cpar%zs_init_hdf(i)), chainfile, initsamp)
            inquire(file=trim(chainfile), exist=exist)
            if (.not. exist) call report_error('Zodi init chain does not exist = ' // trim(chainfile))
            l = len(trim(chainfile))
            if (.not. ((trim(chainfile(l-2:l)) == '.h5') .or. (trim(chainfile(l-3:l)) == '.hd5'))) call report_error('Zodi init chain must be a .h5 file')
            call open_hdf_file(trim(chainfile), file, "r")

            call int2string(initsamp, itext)

            tod_path = trim(adjustl(itext//"/zodi/tod/"))
            band_path = trim(adjustl(tod_path))//trim(adjustl(tod%freq))

            if (.not. hdf_group_exists(file, band_path)) then
               print *, "Zodi init chain does not contain emissivities or albedos for band: " // trim(adjustl(tod%freq))
               stop
            end if

            comp_path = trim(adjustl(band_path))//'/'//trim(adjustl(model%comp_labels(i)))//'/'
            if (hdf_group_exists(file, comp_path)) then
               call read_hdf(file, trim(adjustl(comp_path))//'/emissivity', tod%zodi_emissivity(i))
               call read_hdf(file, trim(adjustl(comp_path))//'/albedo', tod%zodi_albedo(i))
            else
               tod%zodi_emissivity(i) = 1.
               tod%zodi_albedo(i) = 0.
            end if
         end if
      end do
      call mpi_bcast(tod%zodi_emissivity, size(tod%zodi_emissivity), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
      call mpi_bcast(tod%zodi_albedo, size(tod%zodi_albedo), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
   end subroutine read_tod_zodi_params
end module comm_tod_zodi_mod
