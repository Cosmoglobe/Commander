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
   public initialize_tod_zodi_mod, get_zodi_emission, update_zodi_splines, output_tod_params_to_hd5, read_tod_zodi_params

   ! Constants
   real(dp) :: R_MIN = 3.d-14, R_CUTOFF = 5.2, EPS = TINY(1.0_dp)

   ! Shared global parameters
   type(spline_type) :: earth_pos_spl_obj(3)
   real(dp), allocatable :: temperature_grid(:), gauss_nodes(:), gauss_weights(:), b_nu_integrals(:)
   real(dp) :: delta_t_reset, min_ipd_temp, max_ipd_temp
   integer(i4b) :: gauss_degree

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

      allocate (b_nu_integrals(n_interp_points))

      ! hyper parameters
      delta_t_reset = cpar%zs_delta_t_reset
      gauss_degree = cpar%zs_los_steps

      ! Set up Gauss-Legendre quadrature
      allocate (gauss_nodes(gauss_degree), gauss_weights(gauss_degree))
      call leggaus(gauss_degree, gauss_nodes, gauss_weights)

      ! Set up interpolation grid for evaluating temperature
      allocate (temperature_grid(n_interp_points))
      call linspace(min_temp, max_temp, temperature_grid)

      ! Read earth position from file and set up spline object
      call initialize_earth_pos_spline(cpar)

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

      integer(i4b) :: i, j, k, pix_at_zodi_nside, lookup_idx, n_tod, ierr, cache_hits
      logical(lgt) :: scattering, use_lowres
      real(dp) :: earth_lon, R_obs, R_max, dt_tod, obs_time, phase_normalization, C0, C1, C2
      real(dp) :: unit_vector(3), X_unit_LOS(3, gauss_degree), X_LOS(3, gauss_degree), obs_pos(3), earth_pos(3)
      real(dp), dimension(gauss_degree) :: R_LOS, T_LOS, density_LOS, solar_flux_LOS, scattering_angle, phase_function, b_nu_LOS

      s_zodi_scat = 0.
      s_zodi_therm = 0.
      n_tod = size(pix, dim=1)

      dt_tod = (1./tod%samprate)*SECOND_TO_DAY ! dt between two samples in units of days (assumes equispaced tods)
      obs_pos = tod%scans(scan)%x0_obs
      earth_pos = tod%scans(scan)%x0_earth
      R_obs = norm2(obs_pos)
      obs_time = tod%scans(scan)%t0(1)
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
               earth_pos(j) = splint_simple(tod%x_earth_spline(j), obs_time)
               obs_pos(j) = splint_simple(tod%x_obs_spline(j), obs_time)
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

         call get_R_max(unit_vector, obs_pos, R_obs, R_max)

         do k = 1, 3
            ! Convert quadrature range from [-1, 1] to [R_min, R_max]
            X_unit_LOS(k, :) = (0.5 * (R_max - R_MIN)) * gauss_nodes + (0.5 * (R_max + R_MIN))
            X_unit_LOS(k, :) = X_unit_LOS(k, :) * unit_vector(k)
            X_LOS(k, :) = X_unit_LOS(k, :) + obs_pos(k)
         end do
         R_LOS = norm2(X_LOS, dim=1)

         if (scattering) then
            solar_flux_LOS = model%F_sun(tod%band)/R_LOS**2
            call get_scattering_angle(X_LOS, X_unit_LOS, R_LOS, scattering_angle)
            call get_phase_function(scattering_angle, C0, C1, C2, phase_normalization, phase_function)
         end if

         call get_dust_grain_temperature(R_LOS, T_LOS, model%T_0, model%delta)
         call splint_simple_multi(tod%zodi_b_nu_spl_obj(det), T_LOS, b_nu_LOS)

         do k = 1, model%n_comps
            ! If comp is present we only evaluate the zodi emission for that component.
            ! If comp == 0 then we evaluate the zodi emission for all components.
            if (present(comp)) then
               if (k /= comp .and. comp /= 0) cycle
            end if
            call model%comps(k)%c%get_density(X_LOS, earth_lon, density_LOS)
            if (scattering) then
               s_zodi_scat(i, k) = sum(density_LOS*solar_flux_LOS*phase_function*gauss_weights) * 0.5*(R_max - R_MIN) * 1d20
               if (use_lowres) then
                  tod%zodi_scat_cache_lowres(lookup_idx, k, det) = s_zodi_scat(i, k)
               else
                  tod%zodi_scat_cache(lookup_idx, k, det) = s_zodi_scat(i, k)
               end if
            end if
            s_zodi_therm(i, k) = sum(density_LOS*b_nu_LOS*gauss_weights) * 0.5 * (R_max - R_MIN) * 1d20
            if (use_lowres) then
               tod%zodi_therm_cache_lowres(lookup_idx, k, det) = s_zodi_therm(i, k)
            else
               tod%zodi_therm_cache(lookup_idx, k, det) = s_zodi_therm(i, k)
            end if
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
      b = 2.*(obs_pos(1)*cos_lat*cos(lon) + obs_pos(2)*cos_lat*sin(lon))
      d = R_obs**2 - R_CUTOFF**2
      q = -0.5*b*(1.+sqrt(b**2 - (4.*d))/abs(b))
      R_max = max(q, d/q)
   end subroutine get_R_max

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
      real(dp), dimension(gauss_degree) :: cos_theta, R_LOS

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

   subroutine initialize_earth_pos_spline(cpar)
      ! Returns the spline object which is used to evaluate the earth position

      type(comm_params), intent(in) :: cpar

      integer :: i, n_earthpos, unit
      real(dp), allocatable :: tabulated_earth_time(:), tabulated_earth_pos(:, :)
      unit = getlun()
      open (unit, file=trim(trim(cpar%datadir)//'/'//trim("earth_pos_1980-2050_ephem_de432s.txt")))
      read (unit, *) n_earthpos
      read (unit, *) ! skip header

      allocate (tabulated_earth_pos(3, n_earthpos))
      allocate (tabulated_earth_time(n_earthpos))
      do i = 1, n_earthpos
         read (unit, *) tabulated_earth_time(i), tabulated_earth_pos(1, i), tabulated_earth_pos(2, i), tabulated_earth_pos(3, i)
      end do
      close (unit)
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

      class(comm_tod), intent(inout) :: tod
      class(comm_bp_ptr), intent(in) :: bandpass
      integer(i4b), intent(in) :: det
      type(ZodiModel), intent(inout) :: model

      real(dp), allocatable :: b_nu(:)
      integer(i4b) :: i, j

      allocate (b_nu(bandpass%p%n))
      do i = 1, size(b_nu_integrals)
         !write(*,*) 'alloc2', allocated(bandpass%p%nu)
         call get_blackbody_emission(bandpass%p%nu, temperature_grid(i), b_nu)
         b_nu_integrals(i) = tsum(bandpass%p%nu, bandpass%p%tau*b_nu)
      end do
      call spline_simple(tod%zodi_b_nu_spl_obj(det), temperature_grid, b_nu_integrals, regular=.true.)
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
