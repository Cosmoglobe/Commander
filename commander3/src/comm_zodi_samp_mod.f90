module comm_zodi_samp_mod
   use comm_tod_mod
   use comm_utils
   use comm_data_mod
   use comm_tod_zodi_mod
   use comm_zodi_mod
   use comm_tod_driver_mod
   use comm_chisq_mod
   use powell_mod
   implicit none

   private
   public initialize_zodi_samp_mod, downsamp_invariant_structs, project_and_downsamp_sky, compute_downsamp_zodi, sample_zodi_group, sample_linear_zodi, zodi_model_to_ascii, ascii_to_zodi_model

   real(dp), allocatable :: chisq_previous, step_size, priors(:, :), step_sizes_emissivity(:, :), step_sizes_albedo(:, :), step_sizes_ipd(:), step_sizes_n0(:)
   integer(i4b) :: n_samp_bands, reference_emissivity_band
   real(dp) :: NO_PRIOR = HUGE(1.0_dp)
   real(dp), dimension(2) :: emissivity_prior, albedo_prior

contains
   subroutine initialize_zodi_samp_mod(cpar)
      ! Initialize the zodi sampling module.
      !
      ! Parameters
      ! ----------
      ! cpar: comm_params
      !    Parameter file variables.

      type(comm_params), intent(inout) :: cpar
      integer(i4b) :: i, idx_start, idx_stop, ref_band_count
      real(dp), allocatable :: param_vec(:)

      ! Figure out how many sampling bands there are and initialize the tod step sizes
      n_samp_bands = 0
      do i = 1, numband
         if (data(i)%tod_type == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle
         n_samp_bands = n_samp_bands + 1
      end do

      ref_band_count = count(cpar%ds_zodi_reference_band == .true.)
      if (ref_band_count > 1) then
         stop "Error: cannot have more than one reference band for zodi emissivity."
      else if (ref_band_count == 0) then
         stop "Error: cannot sample zodi without the reference band being active."
      end if
      

      allocate (step_sizes_albedo(n_samp_bands, zodi_model%n_comps))
      allocate (step_sizes_emissivity(n_samp_bands, zodi_model%n_comps))
      allocate (step_sizes_n0(zodi_model%n_comps))

      do i = 1, zodi_model%n_comps
         step_sizes_n0(i) = 0.1*zodi_model%comps(i)%c%n_0
      end do

      where (cpar%zs_comp_params == 0)
         cpar%zs_comp_params = NO_PRIOR
      end where
      where (cpar%zs_general_params == 0)
         cpar%zs_general_params = NO_PRIOR
      end where

      allocate(param_vec(zodi_model%n_params))
      call zodi_model%model_to_params(param_vec)
      allocate (step_sizes_ipd(zodi_model%n_params))
      step_sizes_ipd = 0.1*param_vec

      allocate (priors(zodi_model%n_params, 2))
      idx_start = 1
      do i = 1, zodi_model%n_comps
         idx_stop = idx_start + size(zodi_model%comps(i)%labels) - 1
         priors(idx_start:idx_stop, :) = cpar%zs_comp_params(i, :size(zodi_model%comps(i)%labels), 2:)
         idx_start = idx_stop + 1
      end do 
      do i = 1, zodi_model%n_general_params
         priors(idx_start, :) = cpar%zs_general_params(i, 2:)
         idx_start = idx_start + 1
      end do

      emissivity_prior = [0., 5.]
      albedo_prior = [0., 1.]
   end subroutine initialize_zodi_samp_mod

   function get_boxwidth(samprate_lowres, samprate) result(box_width)
      ! Returns the boxcar width for downsampling the zodi tods
      real(dp), intent(in) :: samprate_lowres, samprate
      real(dp) :: box_width

      if (samprate_lowres >= samprate) stop "Cannot downsample zodi tods if lowres samprate is greater than highres samprate!"
      box_width = samprate/samprate_lowres
      box_width = real(nint(box_width))
      if (box_width < 1.) stop "Cannot downsample zodi tods if box car width is less than 1 sample!"
      if (mod(int(box_width), 2) == 0) box_width = box_width + 1.
   end function

   subroutine sample_linear_zodi(cpar, handle, gibbs_iter, model, verbose)
      type(comm_params), intent(in) :: cpar
      type(planck_rng), intent(inout) :: handle
      integer(i4b), intent(in) :: gibbs_iter
      type(ZodiModel), intent(inout) :: model
      logical(lgt), intent(in), optional :: verbose

      integer(i4b) :: i, j, k, ndet, nscan, ntod, nprop, scan, ierr, n_accepted
      integer(i4b) :: albedo_idx, emissivity_idx
      real(sp), allocatable :: s_zodi_scan(:)
      real(dp) :: chisq_tod, chisq_current, chisq_diff, ln_acceptance_probability, accept_rate, chisq_lnL, box_width
      real(dp), allocatable :: param_vec(:), powell_vec(:), emissivity_new(:, :), albedo_new(:, :), emissivity_prev(:, :), albedo_prev(:, :), n0_new(:), n0_prev(:), n0_initial(:)
      logical(lgt) :: accepted, verbose_, tuning, skip
      type(hdf_file) :: tod_file

      if (present(verbose)) then
         verbose_ = verbose
      else
         verbose_ = .false.
      end if
      if (verbose_ .and. (cpar%myid == cpar%root)) print *, "sampling emissivities and albeods"

      nprop = 250
      tuning = gibbs_iter <= 25

      ! Update the cached downsampled zodi with the new model
      call compute_downsamp_zodi(cpar, model)

      ! Start by sampling n_0 at reference band
      allocate (n0_new(model%n_comps))
      allocate (n0_prev(model%n_comps))
      allocate (n0_initial(model%n_comps))

      do j = 1, model%n_comps
         n0_initial(j) = model%comps(j)%c%n_0
      end do
      n0_new = n0_initial
      n0_prev = n0_new

      n_accepted = 0
      accept_rate = 0.

      do k = 0, nprop
         skip = .false.
         chisq_tod = 0.
         chisq_current = 0.

         if (k > 0) then
            if (cpar%myid == cpar%root) then
               do j = 1, model%n_comps
                  n0_new(j) = n0_new(j) + (step_sizes_n0(j) * rand_gauss(handle))
               end do
            end if
            call mpi_bcast(n0_new, size(n0_new), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
         end if
         do i = 1, numband
            if (data(i)%tod_type == "none") cycle
            if (.not. data(i)%tod%sample_zodi) cycle

            if (chisq_tod >= 1.d30) exit

            ndet = data(i)%tod%ndet
            nscan = data(i)%tod%nscan

            box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)

            ! Update current_model with new parameters
            do scan = 1, nscan
               do j = 1, ndet
                  if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle

                  ! Evaluate zodi model with newly proposed values for each band and calculate chisq
                  call get_s_zodi_with_n0(&
                      & s_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                      & s_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                      & s_zodi=data(i)%tod%scans(scan)%d(j)%downsamp_zodi, &
                      & emissivity=data(i)%tod%zodi_emissivity, &
                      & albedo=data(i)%tod%zodi_albedo, &
                      & n_ratio=n0_new/n0_initial &
                  &)
                  chisq_tod = chisq_tod + sum( &
                  & ((data(i)%tod%scans(scan)%d(j)%downsamp_tod &
                  &   - data(i)%tod%scans(scan)%d(j)%downsamp_sky &
                  &   - data(i)%tod%scans(scan)%d(j)%downsamp_zodi &
                  & )/(data(i)%tod%scans(scan)%d(j)%N_psd%sigma0/sqrt(box_width)))**2 &
                  &)
               end do
            end do
         end do
         ! Reduce chisq to root process
         call mpi_reduce(chisq_tod, chisq_current, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cpar%root, MPI_COMM_WORLD, ierr)

         ! Use chisq from the first iteration where we dont draw new parameters as the base chisq
         if (k == 0) chisq_previous = chisq_current

         if (k > 0) then
            ! Root checks if the new proposal is accepted
            if (cpar%myid == cpar%root) then
               chisq_diff = max(chisq_current - chisq_previous, 0.)
               ln_acceptance_probability = -0.5*chisq_diff
               accepted = ln_acceptance_probability > log(rand_uni(handle))
               if (accepted) print *, chisq_current
            end if

            call mpi_bcast(accepted, 1, MPI_LOGICAL, cpar%root, cpar%comm_chain, ierr)
            ! Update model if proposal is accepted
            if (accepted) then
               n_accepted = n_accepted + 1
               chisq_previous = chisq_current
               n0_prev = n0_new
            else
               n0_new = n0_prev
            end if
            accept_rate = real(n_accepted)/real(k)
         end if
      end do

      if (verbose_ .and. (cpar%myid == cpar%root)) print *, "n0 accept rate:", (real(n_accepted)/real(nprop))*100

      if (tuning) then
         if (accept_rate < 0.1) then
            step_sizes_n0 = step_sizes_n0 / 2.
         else if (accept_rate > 0.4) then
            step_sizes_n0 = step_sizes_n0 * 2.
         end if
      end if

      ! Update model with new n_0 values
      do j = 1, model%n_comps
         model%comps(j)%c%n_0 = n0_new(j)
      end do

      ! Sample emissivity and albedo
      allocate (emissivity_new(n_samp_bands, model%n_comps))
      allocate (albedo_new(n_samp_bands, model%n_comps))
      allocate (emissivity_prev(n_samp_bands, model%n_comps))
      allocate (albedo_prev(n_samp_bands, model%n_comps))

      accept_rate = 0
      n_accepted = 0
      do i = 1, numband
         if (data(i)%tod_type == "none") cycle
         if (.not. data(i)%tod%sample_zodi) cycle
         accept_rate = 0.

         ndet = data(i)%tod%ndet
         nscan = data(i)%tod%nscan

         box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)

         emissivity_new(i, :) = data(i)%tod%zodi_emissivity
         albedo_new(i, :) = data(i)%tod%zodi_albedo
         emissivity_prev = emissivity_new
         albedo_prev = albedo_new
         n_accepted = 0
         do k = 0, nprop

            skip = .false.
            chisq_tod = 0.
            chisq_current = 0.

            if (k > 0) then
               if (cpar%myid == cpar%root) then
                  sample_loop: do j = 1, model%n_comps
                     if (cpar%ds_zodi_reference_band(i)) then
                        emissivity_new(i, j) = 1.
                     else
                        emissivity_new(i, j) = emissivity_new(i, j) + (step_sizes_emissivity(i, j)*rand_gauss(handle))
                     end if
                     if (emissivity_new(i, j) < emissivity_prior(1) .or. emissivity_new(i, j) > emissivity_prior(2)) then
                        if (tuning) step_sizes_emissivity(i, j) = step_sizes_emissivity(i, j)/2.
                        chisq_tod = 1.d30
                        skip = .true.
                        emissivity_new(i, :) = emissivity_prev(i, :)
                        exit sample_loop
                     end if
                     albedo_new(i, j) = albedo_new(i, j) + (step_sizes_albedo(i, j)*rand_gauss(handle))
                     if (albedo_new(i, j) < albedo_prior(1) .or. albedo_new(i, j) > albedo_prior(2)) then
                        if (tuning) step_sizes_albedo(i, j) = step_sizes_albedo(i, j)/2.
                        chisq_tod = 1.d30
                        skip = .true.
                        albedo_new(i, :) = albedo_prev(i, :)
                        exit sample_loop
                     end if
                  end do sample_loop
               end if
               call mpi_bcast(emissivity_new, size(emissivity_new), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
               call mpi_bcast(albedo_new, size(albedo_new), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
            end if
            if (.not. skip) then
               ! Update current_model with new parameters
               scan_loop: do scan = 1, nscan
                  do j = 1, ndet
                     if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle

                     ! Evaluate zodi model with newly proposed values for each band and calculate chisq
                     call get_s_zodi_with_n0(&
                         & s_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                         & s_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                         & s_zodi=data(i)%tod%scans(scan)%d(j)%downsamp_zodi, &
                         & emissivity=emissivity_new(i, :), &
                         & albedo=albedo_new(i, :), &
                         & n_ratio=n0_new/n0_initial &
                     &)
                     chisq_tod = chisq_tod + sum( &
                        & ((data(i)%tod%scans(scan)%d(j)%downsamp_tod &
                        &   - data(i)%tod%scans(scan)%d(j)%downsamp_sky &
                        &   - data(i)%tod%scans(scan)%d(j)%downsamp_zodi &
                        & )/(data(i)%tod%scans(scan)%d(j)%N_psd%sigma0/sqrt(box_width)))**2 &
                     &)
                     ! If chisq is already too large, skip rest of the evaluation and go directly to rejection
                     if (chisq_tod >= 1.d30) exit scan_loop
                  end do
               end do scan_loop
            end if

            ! Reduce chisq to root process
            call mpi_reduce(chisq_tod, chisq_current, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cpar%root, MPI_COMM_WORLD, ierr)
            ! Use chisq from the first iteration where we dont draw new parameters as the base chisq
            if (k == 0) then
               chisq_previous = chisq_current
            else if (k > 0) then
               ! Root checks if the new proposal is accepted
               if (cpar%myid == cpar%root) then
                  chisq_diff = max(chisq_current - chisq_previous, 0.)
                  ln_acceptance_probability = -0.5*chisq_diff
                  accepted = ln_acceptance_probability > log(rand_uni(handle))
               end if

               call mpi_bcast(accepted, 1, MPI_LOGICAL, cpar%root, cpar%comm_chain, ierr)
               ! Update model if proposal is accepted
               if (accepted) then
                  n_accepted = n_accepted + 1
                  chisq_previous = chisq_current
                  emissivity_prev(i, :) = emissivity_new(i, :)
                  albedo_prev(i, :) = albedo_new(i, :)
               else
                  emissivity_new(i, :) = emissivity_prev(i, :)
               end if
               accept_rate = real(n_accepted)/real(k)
            end if
         end do
         if (accept_rate < 0.1) then
            step_sizes_emissivity(i, :) = step_sizes_emissivity(i, :)/2.
            step_sizes_albedo(i, :) = step_sizes_albedo(i, :)/2.
         else if (accept_rate > 0.4) then
            step_sizes_emissivity(i, :) = step_sizes_emissivity(i, :)*2.
            step_sizes_albedo(i, :) = step_sizes_albedo(i, :)*2.
         end if
            if (verbose_ .and. (cpar%myid == cpar%root)) print *, "emissivity and albedo accept rate:", (real(n_accepted) / real(k))*100., "band: ", trim(adjustl(data(i)%tod%freq))
      end do
      do i = 1, numband
         data(i)%tod%zodi_emissivity = emissivity_new(i, :)
         data(i)%tod%zodi_albedo = albedo_new(i, :)
      end do
   end subroutine

   subroutine sample_zodi_group(cpar, handle, gibbs_iter, model, verbose)
      type(comm_params), intent(in) :: cpar
      type(planck_rng), intent(inout) :: handle
      integer(i4b), intent(in) :: gibbs_iter
      type(ZodiModel), intent(inout) :: model
      logical(lgt), intent(in), optional :: verbose

      integer(i4b) :: i, j, k, prop, group_idx, flag, ndet, nscan, ntod, n_proposals, scan, ierr, n_accepted, param_idx
      real(dp) :: chisq_tod, chisq_current, chisq_diff, ln_acceptance_probability, accept_rate, chisq_lnL, box_width
      real(dp), allocatable :: param_vec(:)
      integer(i4b), allocatable :: group_indices(:)
      type(ZodiModel) :: current_model, previous_model
      logical(lgt) :: accepted, verbose_, tuning, skip
      type(hdf_file) :: tod_file
      character(len=128), allocatable :: param_labels(:)

      if (present(verbose)) then
         verbose_ = verbose
      else
         verbose_ = .false.
      end if
      if (verbose_ .and. (cpar%myid == cpar%root)) print *, "sampling interplanetary dust parameters"

      current_model = model
      previous_model = model

      tuning = gibbs_iter <= 25
      n_proposals = 25

      allocate(param_vec(current_model%n_params))
      call current_model%model_to_params(param_vec, param_labels)

      ! sample all parameters in a group jointly
      do group_idx = 1, cpar%zs_num_samp_groups
         if (cpar%myid == cpar%root) print *, "sampling zodi group: ", group_idx, " of ", cpar%zs_num_samp_groups

         ! Get indices of the parameters to sample in the current sampling group
         call parse_samp_group_strings(cpar%zs_samp_groups(group_idx), param_labels, group_indices)

         n_accepted = 0
         do prop = 0, n_proposals
            ! if propsed parameter is outside of priors we dont want to immediately reject
            skip = .false.

            ! Reset chisq for current proposal
            chisq_tod = 0.
            chisq_current = 0.

            ! If first iteration we dont want to draw, just compute the chisq with the base model
            if (prop > 0) then
               call current_model%model_to_params(param_vec)
               ! Root process draws new set of zodi parameters and broadcasts to all processes
               if (cpar%myid == cpar%root) then
                  do i = 1, size(group_indices)
                     param_idx = group_indices(i)
                     param_vec(param_idx) = param_vec(param_idx) + (step_sizes_ipd(param_idx)*rand_gauss(handle))
                     skip = prior_is_violated(param_vec(param_idx), priors(param_idx, :))
                     if (skip) then
                        chisq_tod = 1.d30
                        exit
                     end if
                  end do
               end if
               call mpi_bcast(param_vec, size(param_vec), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
               call current_model%params_to_model(param_vec)
            end if

            do i = 1, numband
               if (skip) exit ! drop evaluation if prior is violated
               if (data(i)%tod_type == "none") cycle
               if (.not. data(i)%tod%sample_zodi) cycle

               ! If chisq is already too large, skip rest of the evaluation and go directly to rejection
               if (chisq_tod >= 1.d30) exit

               ndet = data(i)%tod%ndet
               nscan = data(i)%tod%nscan

               box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)

               ! Make sure that the zodi cache is cleared before each new band
               call data(i)%tod%clear_zodi_cache()

               ! Evaluate zodi model with newly proposed values for each band and calculate chisq
               do scan = 1, nscan
                  ! Skip scan if no accepted data
                  do j = 1, ndet
                     if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle

                     call get_zodi_emission(&
                         & tod=data(i)%tod, &
                         & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pix, &
                         & scan=scan, &
                         & det=j, &
                         & s_zodi_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                         & s_zodi_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                         & model=current_model, &
                         & use_lowres_pointing=.true. &
                         ! & comp=comp &
                     &)
                     call get_s_zodi(&
                         & s_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                         & s_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                         & s_zodi=data(i)%tod%scans(scan)%d(j)%downsamp_zodi, &
                         & emissivity=data(i)%tod%zodi_emissivity, &
                         & albedo=data(i)%tod%zodi_albedo &
                     &)

                     chisq_tod = chisq_tod + sum( &
                        & ((data(i)%tod%scans(scan)%d(j)%downsamp_tod &
                        &   - data(i)%tod%scans(scan)%d(j)%downsamp_sky &
                        &   - data(i)%tod%scans(scan)%d(j)%downsamp_zodi &
                        & )/(data(i)%tod%scans(scan)%d(j)%N_psd%sigma0/sqrt(box_width)))**2 &
                     &)
                  end do
               end do
            end do

            ! Reduce chisq to root process
            call mpi_reduce(chisq_tod, chisq_current, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cpar%root, MPI_COMM_WORLD, ierr)

            ! Use chisq from the first iteration where we dont draw new parameters as the base chisq
            if (prop == 0) chisq_previous = chisq_current

            if (prop > 0) then
               ! Root checks if the new proposal is accepted
               if (cpar%myid == cpar%root) then
                  chisq_diff = max(chisq_current - chisq_previous, 0.)
                  ln_acceptance_probability = -0.5*chisq_diff
                  select case (cpar%operation)
                  case ("sample")
                     accepted = ln_acceptance_probability > log(rand_uni(handle))
                  case ("optimize")
                     accepted = chisq_diff < 0.
                  case default
                     stop "Error: invalid operation"
                  end select

                  if (verbose_ .and. accepted) then
                     print *, "proposal:", prop
                     print *, "chisq", chisq_current
                     print *, "accept rate:", (real(n_accepted)/real(prop))*100.
                     print *, " "
                  end if
               end if

               call mpi_bcast(accepted, 1, MPI_LOGICAL, cpar%root, cpar%comm_chain, ierr)
               ! Update model if proposal is accepted
               if (accepted) then
                  n_accepted = n_accepted + 1
                  chisq_previous = chisq_current
                  previous_model = current_model
               else
                  current_model = previous_model
               end if
               accept_rate = real(n_accepted)/real(prop)
            end if
         end do
         if (accept_rate < 0.1) then
            do i = 1, size(group_indices)
               param_idx = group_indices(i)
               step_sizes_ipd(param_idx) = step_sizes_ipd(param_idx)/2.
            end do
         else if (accept_rate > 0.6) then
            do i = 1, size(group_indices)
               param_idx = group_indices(i)
               step_sizes_ipd(param_idx) = step_sizes_ipd(param_idx)*2.
            end do
         end if
         model = current_model
      end do
   end subroutine

   function prior_is_violated(value, prior) result(violated)
      real(dp), intent(in) :: value
      real(dp), intent(in) :: prior(2)
      logical(lgt) :: violated

      if (value < prior(1)) then
         violated = prior(1) /= NO_PRIOR
      else if (value > prior(2)) then
         violated = prior(2) /= NO_PRIOR
      end if
   end function

   subroutine downsamp_invariant_structs(cpar)
      ! Downsamples pointing, tod and caches the zodi mask (sky mask + flags)
      type(comm_params), intent(in) :: cpar

      integer(i4b) :: i, j, k, scan, npix, nmaps, ndelta, ext(2), upper_bound, padding, ierr, ntod, ndet, nhorn, ndownsamp
      real(dp) :: box_width, dt_tod
      real(sp), allocatable :: tod(:), mask(:)
      integer(i4b), allocatable :: pix(:, :), psi(:, :), flag(:)
      real(dp), allocatable, dimension(:, :) :: m_buf
      integer(i4b), allocatable, dimension(:) :: downsamp_pix
      logical(lgt), allocatable, dimension(:) :: downsamp_mask_idx
      real(sp), allocatable, dimension(:) :: downsamp_mask, downsamp_tod, downsamp_obs_time, obs_time
      real(sp), allocatable, dimension(:) :: procmask_zodi
      type(hdf_file) :: tod_file

      padding = 5
      if (cpar%myid == cpar%root) print *, "downsampling tod and pointing"
      ! For each zodi band, create downsampled residual time stream

      do i = 1, numband
         ! Only generate downsampled arrays for tod bands and bands where we have zodi
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         ndelta = 1
         npix = 12*data(i)%map%info%nside**2
         nmaps = data(i)%map%info%nmaps
         ndet = data(i)%tod%ndet
         nhorn = data(i)%tod%nhorn

         ! Build zodi mask
         allocate (m_buf(0:npix - 1, nmaps))
         allocate (procmask_zodi(0:npix - 1))
         call data(i)%tod%procmask_zodi%bcast_fullsky_map(m_buf); procmask_zodi = m_buf(:, 1)

         box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)
         do scan = 1, data(i)%tod%nscan
            if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
            ntod = data(i)%tod%scans(scan)%ntod
            ! Allocate temporary variables to be downsampled
            allocate (pix(ntod, nhorn))
            allocate (psi(ntod, nhorn))
            allocate (flag(ntod))
            allocate (tod(ntod))
            allocate (mask(ntod))

            ! Get shape of downsampled structs
            call data(i)%tod%downsample_tod(tod, ext, step=box_width)

            ! Allocate downsampled strucs
            allocate (downsamp_tod(ext(1):ext(2)))
            allocate (downsamp_mask(ext(1):ext(2)))
            allocate (downsamp_pix(ext(1):ext(2)))
            allocate (downsamp_mask_idx(ext(1):ext(2)))

            ! downsample obs_time
            allocate (obs_time(ntod))
            dt_tod = (1./data(i)%tod%samprate)*SECOND_TO_DAY ! dt between two samples in units of days (assumes equispaced tods)
            do k = 1, ntod
               obs_time(k) = data(i)%tod%scans(scan)%t0(1) + (k - 1)*dt_tod
            end do
            allocate (downsamp_obs_time(ext(1):ext(2)))
            call data(i)%tod%downsample_tod(obs_time, ext, downsamp_obs_time, step=box_width)

            do j = 1, data(i)%tod%ndet
               pix = 0.
               psi = 0.
               flag = 0
               tod = 0.
               mask = 0.

               downsamp_tod = 0.
               downsamp_mask = 0.
               downsamp_pix = 0
               downsamp_mask_idx = .true.

               ! Decompress pointing and flags
               call data(i)%tod%decompress_pointing_and_flags(scan, j, pix, psi, flag)

               ! Decompress tod if necessary
               if (data(i)%tod%compressed_tod) then
                  call data(i)%tod%decompress_tod(scan, j, tod)
               else
                  tod = data(i)%tod%scans(scan)%d(j)%tod
               end if

               ! Cache zodi sky mask + flags
               do k = 1, data(i)%tod%scans(scan)%ntod
                  mask(k) = procmask_zodi(pix(k, 1))
                  if (iand(flag(k), data(i)%tod%flag0) .ne. 0) mask(k) = 0.
               end do

               ! Downsample the mask
               call data(i)%tod%downsample_tod(mask, ext, downsamp_mask, step=box_width)

               ! Downsample the tods
               call data(i)%tod%downsample_tod(tod, ext, downsamp_tod, step=box_width)

               ! Construct the downsampled pix array (pick central pixel in each box)
               do k = 1, int(size(tod)/box_width)
                  downsamp_pix(k) = pix((k*box_width), 1)
               end do
               where (downsamp_mask < 1.) downsamp_mask_idx = .false.

               ! Find upper bound and truncate the downsampled arrays by where they were padded and
               ! filter out the masked values. Store the result in the tod objects for use in sampling.
               upper_bound = ubound(downsamp_mask, dim=1)

                    data(i)%tod%scans(scan)%d(j)%downsamp_pix = pack(downsamp_pix(padding:upper_bound-padding), downsamp_mask_idx(padding:upper_bound-padding))
                    data(i)%tod%scans(scan)%d(j)%downsamp_tod = pack(downsamp_tod(padding:upper_bound-padding), downsamp_mask_idx(padding:upper_bound-padding))
               if (j == 1) then
                        data(i)%tod%scans(scan)%downsamp_obs_time = pack(downsamp_obs_time(padding:upper_bound-padding), downsamp_mask_idx(padding:upper_bound-padding))
               end if

               ! Allocate other downsampled quantities with same shape
               ndownsamp = size(data(i)%tod%scans(scan)%d(j)%downsamp_pix)
               allocate (data(i)%tod%scans(scan)%d(j)%downsamp_sky(ndownsamp))
               allocate (data(i)%tod%scans(scan)%d(j)%downsamp_zodi(ndownsamp))
               allocate (data(i)%tod%scans(scan)%d(j)%downsamp_scat(ndownsamp, zodi_model%n_comps))
               allocate (data(i)%tod%scans(scan)%d(j)%downsamp_therm(ndownsamp, zodi_model%n_comps))
               deallocate (pix, psi, flag, tod, mask, downsamp_tod, downsamp_mask, downsamp_pix, downsamp_mask_idx)
               deallocate (obs_time, downsamp_obs_time)
            end do
         end do

         deallocate (m_buf, procmask_zodi)
      end do
   end subroutine

   subroutine project_and_downsamp_sky(cpar)
      type(comm_params), intent(in) :: cpar

      integer(i4b) :: i, j, k, pix, scan, ext(2), upper_bound, padding, ierr, ntod, ndet, nmaps, ndelta
      real(dp) :: box_width
      type(map_ptr), allocatable, dimension(:, :) :: sky_signal
      real(sp), allocatable, dimension(:) :: downsamp_s_sky
      real(sp), allocatable, dimension(:, :, :, :) :: map_sky
      type(hdf_file) :: tod_file

      if (cpar%myid == cpar%root) print *, "projecting downsampled sky"

      padding = 5
      ! For each zodi band, create downsampled residual time stream
      do i = 1, numband
         ! Only generate downsampled arrays for tod bands and bands where we have zodi
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         ndelta = 1
         nmaps = data(i)%map%info%nmaps
         ndet = data(i)%tod%ndet
         allocate (sky_signal(data(i)%tod%ndet, ndelta))
         do j = 1, data(i)%tod%ndet
            call get_sky_signal(i, j, sky_signal(j, ndelta)%p, mono=.false.)
         end do

         allocate (map_sky(nmaps, data(i)%tod%nobs, 0:data(i)%tod%ndet, 1))
         call distribute_sky_maps(data(i)%tod, sky_signal, 1.e0, map_sky)

         do scan = 1, data(i)%tod%nscan
            if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle

            ! Project sky map down to downsampled timestream and cache
            do j = 1, ndet
               do k = 1, size(data(i)%tod%scans(scan)%d(1)%downsamp_pix)
                  pix = data(i)%tod%pix2ind(data(i)%tod%scans(scan)%d(j)%downsamp_pix(k))
                  data(i)%tod%scans(scan)%d(j)%downsamp_sky(k) = map_sky(1, pix, j, 1)  !zodi is only temperature (for now)
               end do
            end do
         end do
         deallocate (map_sky, sky_signal)
      end do
   end subroutine

   subroutine compute_downsamp_zodi(cpar, model)
      type(comm_params), intent(in) :: cpar
      type(ZodiModel), intent(inout) :: model
      integer(i4b) :: i, j, scan, ierr, ndet

      if (cpar%myid == cpar%root) print *, "computing downsampled zodi"

      do i = 1, numband
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         do scan = 1, data(i)%tod%nscan
            if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle

            ndet = data(i)%tod%ndet
            do j = 1, ndet
               call get_zodi_emission(&
                   & tod=data(i)%tod, &
                   & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pix, &
                   & scan=scan, &
                   & det=j, &
                   & s_zodi_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                   & s_zodi_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                   & model=model, &
                   & use_lowres_pointing=.true. &
               &)
            end do
         end do
      end do
   end subroutine

   subroutine get_s_zodi_with_n0(s_therm, s_scat, s_zodi, emissivity, albedo, n_ratio)
      real(sp), dimension(:, :), intent(in) :: s_scat, s_therm
      real(sp), dimension(:), intent(inout) :: s_zodi
      real(dp), dimension(:), intent(in) :: emissivity, albedo, n_ratio
      integer(i4b) :: i, n_comps

      n_comps = size(emissivity)
      s_zodi = 0.
      do i = 1, n_comps
         s_zodi = s_zodi + (s_scat(:, i)*albedo(i) + (1.-albedo(i))*emissivity(i)*s_therm(:, i))*n_ratio(i)
      end do
   end subroutine get_s_zodi_with_n0

   subroutine parse_samp_group_strings(samp_group_str, param_labels, param_indices)
      character(len=*), intent(in) :: samp_group_str
      character(len=*), intent(in) :: param_labels(:)
      integer(i4b), allocatable, intent(inout) :: param_indices(:)
      character(len=128) :: tokens(100), comp_param(2), label, param_label_tokens(10)
      character(len=128), allocatable :: tokens_trunc(:)
      integer(i4b) :: i, j, n_params
      logical(lgt) :: found


      if (allocated(param_indices)) deallocate(param_indices)

      call get_tokens(samp_group_str, ',', tokens, n_params) 
      tokens_trunc = tokens(1:n_params)
      allocate(param_indices(n_params))
      param_indices = -1
      do i = 1, size(tokens_trunc)
         call get_tokens(tokens_trunc(i), ':', comp_param) 
         call toupper(comp_param(1))
         call toupper(comp_param(2))
         if (trim(adjustl(comp_param(2))) == "ALL") then
            do j = 1, size(param_labels)
               call get_tokens(param_labels(j), "_", param_label_tokens) 
               call toupper(param_label_tokens(1))
               if (trim(adjustl(param_label_tokens(2))) == "N_0") then ! dont add N_0 to sampling groups with comp:all
                  print *, param_labels(j)
                  cycle
               end if
               if (trim(adjustl(comp_param(1))) == param_label_tokens(1)) then
                  if (.not. any(param_indices == j)) param_indices = [param_indices , j]
               end if
            end do
         else
            found = .false.
            label = trim(adjustl(comp_param(1)))//"_"//trim(adjustl(comp_param(2)))
            do j = 1, size(param_labels)
               if (label == param_labels(j)) then
                  if (.not. any(param_indices == j)) param_indices(i) = j
                  found = .true.
                  exit
               end if
            end do
            if (.not. found) then
               print *, "Error: invalid zodi sampling group parameter label :" // trim(adjustl(label)) 
               stop
            end if 
         end if
      end do
      param_indices = pack(param_indices, param_indices > 0)
   end subroutine


   ! function get_powell_vec(param_vec) result(powell_vec)
   !     real(dp), allocatable :: param_vec(:)
   !     real(dp), allocatable :: powell_vec(:)
   !     powell_vec = pack(param_vec, active_params == .true.)
   ! end function get_powell_vec

   ! subroutine update_param_vec_from_powell_vec(param_vec, powell_vec)
   !     real(dp), intent(inout) :: param_vec(:)
   !     real(dp), intent(in) :: powell_vec(:)
   !     integer :: i, j
   !     j = 1
   !     do i = 1, size(param_vec)
   !         if (active_params(i) == .true.) then
   !             param_vec(i) = powell_vec(j)
   !             j = j + 1
   !         end if
   !     end do
   ! end subroutine update_param_vec_from_powell_vec

   ! function lnL_zodi(p)
   !     use healpix_types
   !     implicit none
   !     real(dp), dimension(:), intent(in), optional :: p
   !     real(dp) :: lnL_zodi

   !     real(dp), allocatable :: theta(:)
   !     real(sp), allocatable :: s_scat(:, :), s_therm(:, :), s_zodi(:), res(:)
   !     type(ZodiModel) :: model

   !     integer(i4b) :: i, j, ntod, ndet, nscan, scan, ierr, flag
   !     real(dp) :: chisq_tod, chisq, chisq_buff
   !     real(dp), allocatable :: p_copy(:), param_vec(:)
   !     allocate(theta(n_active_params))

   !     if (data(1)%tod%myid == 0) then
   !         flag = 1
   !         call mpi_bcast(flag, 1, MPI_INTEGER, 0, data(1)%tod%comm, ierr)
   !         theta = p
   !     end if
   !     call mpi_bcast(theta, size(theta), MPI_DOUBLE_PRECISION, 0, data(1)%tod%comm, ierr)

   !     ! Get param vector from zodi model
   !     model = zodi_model
   !     allocate(param_vec(model%n_parameters))
   !     param_vec = model%model_to_param_vec()

   !     ! Check priors
   !     j = 1
   !     do i = 1, size(param_vec)
   !         if (active_params(i) == .true.) then
   !             if (theta(j) < priors(i, 1) .or. theta(j) > priors(i, 2)) then
   !                 lnL_zodi = 1.d30
   !                 return
   !             end if
   !             j = j + 1
   !         end if
   !     end do

   !     ! Update param_vec with new values from powell_vec (theta)
   !     call update_param_vec_from_powell_vec(param_vec, theta)

   !     ! update zodi model with new param_vec
   !     call model%param_vec_to_model(param_vec)

   !     ! Calculate chisq for each band
   !     chisq_tod = 0.
   !     chisq = 0.
   !     lnL_zodi = 0.

   !     do i = 1, numband
   !         ! Skip none tod bands
   !         if (data(i)%tod_type == "none") cycle
   !         ! Skip tod bands where we dont want to sample zodi
   !         if (.not. data(i)%tod%sample_zodi) cycle

   !         ndet = data(i)%tod%ndet
   !         nscan = data(i)%tod%nscan

   !         ! Make sure that the zodi cache is cleared before each new band
   !         call data(i)%tod%clear_zodi_cache()

   !         ! Evaluate zodi model with newly proposed values for each band and calculate chisq
   !         do scan = 1, nscan
   !             ! Skip scan if no accepted data
   !             do j = 1, ndet
   !                 if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
   !                 ntod = size(data(i)%tod%scans(scan)%d(j)%downsamp_res)
   !                 allocate(s_scat(ntod, model%n_comps), s_therm(ntod, model%n_comps), s_zodi(ntod))
   !                 call get_zodi_emission(&
   !                     & tod=data(i)%tod, &
   !                     & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pix, &
   !                     & scan=scan, &
   !                     & det=j, &
   !                     & s_zodi_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
   !                     & s_zodi_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
   !                     & model=model &
   !                 &)
   !                 call get_s_zodi(&
   !                     & s_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
   !                     & s_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
   !                     & s_zodi=s_zodi, &
   !                     & emissivity=data(i)%tod%zodi_emissivity, &
   !                     & albedo=data(i)%tod%zodi_albedo &
   !                 &)
   !                 chisq_buff = sum(((data(i)%tod%scans(scan)%d(j)%downsamp_res - s_zodi)/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2)

   !                 chisq_tod = chisq_tod + chisq_buff
   !                 deallocate(s_scat, s_therm, s_zodi)
   !             end do
   !         end do
   !     end do

   !     ! Reduce chisq to root process
   !     call mpi_reduce(chisq_tod, chisq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, data(1)%tod%comm, ierr)

   !     if (data(1)%tod%myid == 0) then
   !         ! print *, chisq, model%comps(1)%c%n_0
   !         lnL_zodi = chisq
   !     end if

   ! end function lnL_zodi


   subroutine zodi_model_to_ascii(cpar, model, filename)
      ! Dumps the zodi model to an ascii file on the format {COMP}_{PARAM} = {VALUE}.
      class(ZodiModel), target, intent(in) :: model
      type(comm_params), intent(in) :: cpar
      character(len=*), intent(in) :: filename

      integer(i4b) :: io, i, j, running_idx
      logical(lgt) :: exists
      real(dp), allocatable :: params(:)
      integer(i4b), allocatable :: comp_switch_indices(:)
      character(len=128), allocatable :: labels(:)
      character(len=512) :: concatenated_string, val

      if (cpar%myid_chain /= cpar%root) return
      inquire(file=trim(adjustl(filename)), exist=exists)
      if (.not. exists) then
         print *, "zodi asciifile: " // trim(adjustl(filename)) // " does not exist"
         stop
      end if

      open(newunit=io, file=trim(adjustl(filename)), action="write")

      allocate(params(model%n_params))
      call model%model_to_params(params, labels=labels)

      allocate(comp_switch_indices(model%n_comps))

      running_idx = 0
      do i = 1, model%n_comps
         running_idx = running_idx + size(model%comps(i)%labels)
         comp_switch_indices(i) = running_idx
      end do

      do i = 1, model%n_params
         if (any(comp_switch_indices == i)) then
               write(io, fmt='(a, T25, a, ES12.5, a)') trim(adjustl(labels(i))), "= ", params(i), new_line('a')
            else
               write(io, fmt='(a, T25, a, ES12.5)') trim(adjustl(labels(i))), "= ", params(i)
         end if
      end do

      write(io, fmt='(a)') ''
      do i = 1, numband
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle
         concatenated_string = ""
         do j = 1, model%n_comps
            write(val, fmt='(ES12.5)') data(i)%tod%zodi_emissivity(j)
            concatenated_string = trim(adjustl(concatenated_string)) // "," // trim(adjustl(val))
         end do
         write(io, fmt='(a, T25, a, a)') trim(adjustl("EMISSIVITY_"//trim(adjustl(data(i)%tod%freq)))), "= ", trim(adjustl(concatenated_string(2:)))
      end do

      write(io, fmt='(a)') ''
      do i = 1, numband
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle
         concatenated_string = ""
         do j = 1, model%n_comps
            write(val, fmt='(ES12.5)') data(i)%tod%zodi_albedo(j)
            concatenated_string = trim(adjustl(concatenated_string)) // "," // trim(adjustl(val))
         end do
         write(io, fmt='(a, T25, a, a)') trim(adjustl("ALBEDO_"//trim(adjustl(data(i)%tod%freq)))), "= ", trim(adjustl(concatenated_string(2:)))
      end do

      close(io)
   end subroutine

   subroutine ascii_to_zodi_model(cpar, model, filename)
      ! Reads in and updates the zodi model from an ascii file on the format {COMP}_{PARAM} = {VALUE}.
      class(ZodiModel), target, intent(inout) :: model
      type(comm_params), intent(in) :: cpar
      character(len=*), intent(in) :: filename
      type(hash_tbl_sll) :: htbl

      integer(i4b) :: i, j, io, io_status, ierr, n_comps
      logical(lgt) :: exists
      character(len=512) :: key, val, line
      character(len=128), allocatable :: labels(:)
      characteR(len=128) :: toks(100)
      characteR(len=512) :: concatenated_string
      real(dp), allocatable :: params(:)

      allocate(params(model%n_params))
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
               j = index(line, "=")
               if (j == 0) then
                  print *, "Error: invalid line in ascii file: ", trim(adjustl(line))
                  close(io)
                  stop
               end if

               key = trim(adjustl(line(:j-1)))
               val = trim(adjustl(line(j+1:)))
               call tolower(key)
               call put_hash_tbl_sll(htbl, trim(adjustl(key)), trim(adjustl(val))) 
            end if
         end do
         close(io)

         call model%model_to_params(params, labels)
         params = 0.
         if (size(labels) /= size(params)) stop "Error: size of labels and params do not match"
         do i = 1, size(labels)
            call get_parameter_hashtable(htbl, labels(i), par_dp=params(i))
         end do
      end if

      call mpi_bcast(params, size(params), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
      call model%params_to_model(params)

      do i = 1, numband
         if (cpar%myid == 0) then
            if (trim(data(i)%tod_type) == 'none') cycle
            if (.not. data(i)%tod%subtract_zodi) cycle

            call get_parameter_hashtable(htbl, trim(adjustl("EMISSIVITY_"//trim(adjustl(data(i)%tod%freq)))), par_string=concatenated_string)
            call get_tokens(trim(adjustl(concatenated_string)), ',', toks, n_comps)
            if (n_comps /= model%n_comps) stop "Error: number of components in ascii file does not match model emissivity"
            do j = 1, n_comps
               read(toks(j), *) data(i)%tod%zodi_emissivity(j)
            end do

            call get_parameter_hashtable(htbl, trim(adjustl("ALBEDO_"//trim(adjustl(data(i)%tod%freq)))), par_string=concatenated_string)
            call get_tokens(trim(adjustl(concatenated_string)), ',', toks, n_comps)
            if (n_comps /= model%n_comps) stop "Error: number of components in ascii file does not match model albedo"
            do j = 1, n_comps
               read(toks(j), *) data(i)%tod%zodi_albedo(j)
            end do
         end if
         call mpi_bcast(data(i)%tod%zodi_emissivity, size(data(i)%tod%zodi_emissivity), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
         call mpi_bcast(data(i)%tod%zodi_albedo, size(data(i)%tod%zodi_albedo), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
      end do
   end subroutine

end module
