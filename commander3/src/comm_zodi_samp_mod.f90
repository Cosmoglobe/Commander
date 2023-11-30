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
   public initialize_zodi_samp_mod, downsamp_invariant_structs, project_and_downsamp_sky, compute_downsamp_zodi, sample_zodi_group, sample_linear_zodi, zodi_model_to_ascii, ascii_to_zodi_model, minimize_zodi_with_powell, get_chisq_priors, precompute_lowres_zodi_lookups, remove_glitches_from_downsamped_zodi_quantities

   real(dp), allocatable :: chisq_previous, step_size, prior_vec(:, :), prior_vec_powell(:, :), step_sizes_emissivity(:, :), step_sizes_albedo(:, :), step_sizes_ipd(:), step_sizes_n0(:), theta_0(:)
   real(dp), allocatable :: powell_emissivity(:, :), powell_albedo(:, :)
   integer(i4b) :: n_samp_bands, reference_emissivity_band
   real(dp) :: EPS = TINY(1.0_dp)
   real(dp), dimension(2) :: emissivity_prior, albedo_prior
   logical(lgt), allocatable :: powell_included_params(:), ref_band(:)
   character(len=128), allocatable, dimension(:) :: implemented_sampling_algorithms

contains
   subroutine initialize_zodi_samp_mod(cpar)
      ! Initialize the zodi sampling module.
      !
      ! Parameters
      ! ----------
      ! cpar: comm_params
      !    Parameter file variables.

      type(comm_params), intent(inout) :: cpar
      integer(i4b) :: i, j, idx_start, idx_stop, ref_band_count, ierr, group_idx
      real(dp), allocatable :: param_vec(:)
      character(len=128), allocatable :: labels(:)
      integer(i4b), allocatable :: indices(:)
      ! Figure out how many sampling bands there are and initialize the tod step sizes

      implemented_sampling_algorithms = ["powell", "mh"]
      if (.not. any(implemented_sampling_algorithms == cpar%zs_sample_method)) then
         if (cpar%myid == 0) then 
            print *, "Error: invalid sampling method for zodi, must be one of: ", [(trim(adjustl(implemented_sampling_algorithms(i)))//", ", i=1, size(implemented_sampling_algorithms))]
            stop
         end if
      end if
      n_samp_bands = 0
      do i = 1, numband
         if (data(i)%tod_type == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle
         n_samp_bands = n_samp_bands + 1
      end do
      ref_band = cpar%ds_zodi_reference_band
      ref_band_count = count(cpar%ds_zodi_reference_band == .true.)
      if (trim(adjustl(cpar%zs_sample_method)) == "mh") then
         if (ref_band_count > 1) then
            stop "Error: cannot have more than one reference band for zodi emissivity."
         else if (ref_band_count == 0) then
            stop "Error: cannot sample zodi without the reference band being active."
         end if
      end if
      if (trim(adjustl(cpar%zs_sample_method)) == "powell") then
         allocate(powell_albedo(n_samp_bands, zodi_model%n_comps))
         allocate(powell_emissivity(n_samp_bands, zodi_model%n_comps))
      end if

      ! crappy way of picking initial stepsizes for the zodi parameters
      allocate (step_sizes_albedo(n_samp_bands, zodi_model%n_comps))
      allocate (step_sizes_emissivity(n_samp_bands, zodi_model%n_comps))
      allocate (step_sizes_n0(zodi_model%n_comps))

      do i = 1, zodi_model%n_comps
         step_sizes_n0(i) = 0.01*zodi_model%comps(i)%c%n_0
      end do
      do i = 1, numband
         step_sizes_emissivity(i,:) =  0.01 * data(i)%tod%zodi_emissivity
         step_sizes_albedo(i, :) =  0.01 * data(i)%tod%zodi_albedo
      end do

      allocate(param_vec(zodi_model%n_params))
      call zodi_model%model_to_params(param_vec, labels)

      allocate (step_sizes_ipd(zodi_model%n_params))
      step_sizes_ipd = 0.005*param_vec

      allocate (prior_vec(zodi_model%n_params, 3))
      idx_start = 1
      do i = 1, zodi_model%n_comps
         idx_stop = idx_start + size(zodi_model%comps(i)%labels) - 1
         prior_vec(idx_start:idx_stop, :) = cpar%zs_comp_params(i, :size(zodi_model%comps(i)%labels), 2:)
         idx_start = idx_stop + 1
      end do 
      do i = 1, zodi_model%n_general_params
         prior_vec(idx_start, :) = cpar%zs_general_params(i, 2:)
         idx_start = idx_start + 1
      end do

      emissivity_prior = [0., 5.]
      albedo_prior = [0., 1.]

      ! validate sampling group parameters
      do group_idx = 1, cpar%zs_num_samp_groups
         ! Get indices of the parameters to sample in the current sampling group
         call parse_samp_group_strings(cpar%zs_samp_groups(group_idx), labels, indices)
      end do
      ! do i = 1, size(labels)
      !    if (index(trim(adjustl(labels(i))), "X_0") /= 0) then
      !       step_sizes_ipd(i) = step_sizes_ipd(i) * 10.
      !    else if (index(trim(adjustl(labels(i))), "Y_0") /= 0) then
      !       step_sizes_ipd(i) = step_sizes_ipd(i) * 10.
      !    else if (index(trim(adjustl(labels(i))), "Z_0") /= 0) then
      !       step_sizes_ipd(i) = step_sizes_ipd(i) * 10.
      !    end if
      ! end do

      ! if (trim(adjustl(cpar%zs_covar_chain)) /= "none") then
      !    call read_zodi_covariance(cpar)
      ! end if 
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

      !call compute_downsamp_zodi(cpar, model)
      !if (cpar%myid == cpar%root) print *, new_line('A'), "sampling and subtracting monopole"
      !call sample_and_subtract_monopole(cpar, handle)
      if (cpar%myid == cpar%root) print *, new_line('A'), "sampling n0"
      call sample_n0_emissivity_and_albedo(cpar, handle, gibbs_iter, model)
   end subroutine

   subroutine sample_and_subtract_monopole(cpar, handle)
      type(comm_params), intent(in) :: cpar
      type(planck_rng), intent(inout) :: handle
      
      type(hdf_file) :: tod_file
      integer(i4b) :: i, j, prop, nprop, ndet, scan, nscan, ierr, n_accepted
      integer(i4b), allocatable, dimension(:) :: n_downsamp_tod
      real(dp) :: box_width, chisq_scan, chisq_new, chisq_prev, chisq_diff, ln_accept_prob, accept_rate
      real(dp), allocatable, dimension(:) :: monopole_prev, monopole_new, monopole_init, monopole_step, eta
      logical(lgt) :: is_root, accepted, first_prop
      character(len=3) :: itext
      character(len=512) :: path

      nprop = 200
      is_root = cpar%myid == cpar%root
      allocate(monopole_init(n_samp_bands), monopole_new(n_samp_bands), monopole_prev(n_samp_bands), eta(n_samp_bands), n_downsamp_tod(n_samp_bands))
      n_downsamp_tod = 0
      monopole_init = 0.
      monopole_new = 0.
      monopole_prev = monopole_new

      n_accepted = 0
      do prop = 0, nprop
         first_prop = prop == 0
         if (.not. first_prop .and. is_root) then
            eta = [(rand_gauss(handle), i = 1, n_samp_bands)]
            monopole_new = monopole_prev + monopole_step * eta
         else
            monopole_new = monopole_prev
         end if
         call mpi_bcast(monopole_new, 1, MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)

         chisq_scan = 0.
         chisq_new = 0.

         ! compute chisq for new model
         do i = 1, numband
            if ((data(i)%tod_type == "none") .or. (.not. data(i)%tod%sample_zodi)) cycle
            ndet = data(i)%tod%ndet
            nscan = data(i)%tod%nscan
            box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)
            do scan = 1, nscan
               do j = 1, ndet
                  if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
                  
                  n_downsamp_tod(i) = n_downsamp_tod(i) + size(data(i)%tod%scans(scan)%d(j)%downsamp_pix)
                  ! rescale downsampled zodi comp-wise emission with newly proposed n0s
                  call get_s_zodi(&
                     & s_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                     & s_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                     & s_zodi=data(i)%tod%scans(scan)%d(j)%downsamp_zodi, &
                     & emissivity=data(i)%tod%zodi_emissivity, &
                     & albedo=data(i)%tod%zodi_albedo &
                  &)
                  ! compute monopole step size
                  if (first_prop) then 
                     ! print *, "tod", maxval(data(i)%tod%scans(scan)%d(j)%downsamp_tod), minval(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                     ! print *, "sky", maxval(data(i)%tod%scans(scan)%d(j)%downsamp_sky), minval(data(i)%tod%scans(scan)%d(j)%downsamp_sky)
                     ! print *, "zodi", maxval(data(i)%tod%scans(scan)%d(j)%downsamp_zodi), minval(data(i)%tod%scans(scan)%d(j)%downsamp_zodi)
                     ! print *, "pix", maxval(data(i)%tod%scans(scan)%d(j)%downsamp_pix), minval(data(i)%tod%scans(scan)%d(j)%downsamp_pix)
                     ! stop
                     monopole_init(i) = monopole_init(i) + sum( &
                        &  (data(i)%tod%scans(scan)%d(j)%downsamp_tod &
                        &  - data(i)%tod%scans(scan)%d(j)%downsamp_sky &
                        &  - data(i)%tod%scans(scan)%d(j)%downsamp_zodi) &
                     & )
                  end if 

                  chisq_scan = chisq_scan + sum( &
                     & ((data(i)%tod%scans(scan)%d(j)%downsamp_tod &
                     &   - data(i)%tod%scans(scan)%d(j)%downsamp_sky &
                     &   - data(i)%tod%scans(scan)%d(j)%downsamp_zodi &
                     &   - monopole_new(i) &
                     & )/(data(i)%tod%scans(scan)%d(j)%N_psd%sigma0/sqrt(box_width)))**2 &
                  &)
               end do
            end do
         end do
         call mpi_reduce(chisq_scan, chisq_new, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cpar%root, MPI_COMM_WORLD, ierr)

         ! Use chisq from the first iteration where we dont draw new parameters as the base chisq
         if (first_prop) then 
            call mpi_allreduce(MPI_IN_PLACE, monopole_init, size(monopole_init), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            call mpi_allreduce(MPI_IN_PLACE, n_downsamp_tod, size(n_downsamp_tod), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
            monopole_prev = monopole_init/n_downsamp_tod
            monopole_step = 0.1*monopole_prev
            chisq_prev = chisq_new
         else 
            if (is_root) then
               chisq_diff = max(chisq_new - chisq_prev, 0.)
               ln_accept_prob = -0.5 * chisq_diff
               select case (cpar%operation)
                  case ("sample")
                     accepted = ln_accept_prob > log(rand_uni(handle))
                  case ("optimize")
                     accepted = chisq_diff < 0.
                  case default
                     stop "Error: invalid cpar%operation in comm_zodi_samp_mod sample_monopole routine"
                  end select
            end if
            call mpi_bcast(accepted, 1, MPI_LOGICAL, cpar%root, cpar%comm_chain, ierr)

            if (accepted) then
               n_accepted = n_accepted + 1
               chisq_prev = chisq_new
               monopole_prev = monopole_new
               if (is_root) then
                  print *, "proposal:", prop
                  write(*, '(a, ES12.5)') "chisq: ", chisq_new
                  print *, "monopole: ", monopole_new
                  print *, " "
               end if
            else
               monopole_new = monopole_prev
            end if
            accept_rate = real(n_accepted)/real(prop)
         end if
      end do

      if (is_root) print *, "monopole accept rate: ", (real(n_accepted)/real(nprop))*100.

      do i = 1, numband
         if ((data(i)%tod_type == "none") .or. (.not. data(i)%tod%sample_zodi)) cycle
         do scan = 1, nscan
            do j = 1, ndet
               if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
               data(i)%tod%scans(scan)%d(j)%downsamp_tod = data(i)%tod%scans(scan)%d(j)%downsamp_tod - monopole_prev(i)
            end do
         end do
      end do
   end subroutine

   subroutine sample_n0_emissivity_and_albedo(cpar, handle, gibbs_iter, model)
      type(comm_params), intent(in) :: cpar
      type(planck_rng), intent(inout) :: handle
      integer(i4b), intent(in) :: gibbs_iter
      type(ZodiModel), intent(inout) :: model

      integer(i4b) :: i, j, prop, nprop, ndet, scan, nscan, ierr, n_accepted, comp
      real(dp) :: box_width, chisq_scan, chisq_new, chisq_prev, chisq_diff, ln_accept_prob, accept_rate
      real(dp), allocatable :: n0_new(:), n0_prev(:), eta_n0(:), eta_emissivity(:, :), eta_albedo(:, :), emissivity_new(:, :), emissivity_prev(:, :), albedo_new(:, :), albedo_prev(:, :)
      logical(lgt) :: is_root, accepted, tuning, first_prop, scattering

      nprop = 1000
      tuning = gibbs_iter <= 25
      is_root = cpar%myid == cpar%root

      n0_prev = [(model%comps(i)%c%n_0, i = 1, model%n_comps)]
      n0_new = n0_prev

      allocate(emissivity_prev(n_samp_bands, model%n_comps))
      allocate(albedo_prev(n_samp_bands, model%n_comps))
      allocate(eta_emissivity(n_samp_bands, model%n_comps))
      allocate(eta_albedo(n_samp_bands, model%n_comps))
      do i = 1, numband
         emissivity_prev(i, :) = data(i)%tod%zodi_emissivity
         albedo_prev(i, :) = data(i)%tod%zodi_albedo
      end do
      emissivity_new = emissivity_prev
      albedo_new = albedo_prev

      n_accepted = 0
      do prop = 0, nprop
         first_prop = prop == 0
         ! propose new n0s
         if (.not. first_prop) then
            if (is_root) then
               eta_n0 = [(rand_gauss(handle), i = 1, model%n_comps)]
               n0_new = n0_prev + (step_sizes_n0 * eta_n0)
               do i = 1, numband
                  if ((data(i)%tod_type == "none") .or. (.not. data(i)%tod%sample_zodi)) cycle
                  if (cpar%ds_zodi_reference_band(data(i)%id_abs)) then
                     eta_emissivity(i, :) = 0.
                  else 
                     eta_emissivity(i, :) = [(rand_gauss(handle), j = 1, model%n_comps)]
                  end if
                  eta_albedo(i, :) = [(rand_gauss(handle), j = 1, model%n_comps)]
               end do
               emissivity_new = emissivity_prev + (step_sizes_emissivity * eta_emissivity)
               albedo_new = albedo_prev + (step_sizes_albedo * eta_albedo)
            end if
            call mpi_bcast(n0_new, size(n0_new), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
            call mpi_bcast(emissivity_new, size(emissivity_new), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
            call mpi_bcast(albedo_new, size(albedo_new), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
         end if

         chisq_scan = 0.
         chisq_new = 0.
         ! compute chisq for new model
         do i = 1, numband
            if ((data(i)%tod_type == "none") .or. (.not. data(i)%tod%sample_zodi)) cycle
            ndet = data(i)%tod%ndet
            nscan = data(i)%tod%nscan
            box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)
            scattering = any(data(i)%tod%zodi_albedo > EPS)
            do scan = 1, nscan
               do j = 1, ndet
                  if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
                  ! rescale downsampled zodi comp-wise emission with newly proposed n0s
                  call get_s_zodi(&
                     & s_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                     & s_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                     & s_zodi=data(i)%tod%scans(scan)%d(j)%downsamp_zodi, &
                     & emissivity=emissivity_new(i, :), &
                     & albedo=albedo_new(i, :), &
                     & alpha=n0_new/n0_prev &
                  &)
                  chisq_scan = chisq_scan + sum( &
                  & ((data(i)%tod%scans(scan)%d(j)%downsamp_tod &
                  &   - data(i)%tod%scans(scan)%d(j)%downsamp_sky &
                  &   - data(i)%tod%scans(scan)%d(j)%downsamp_zodi &
                  & )/(data(i)%tod%scans(scan)%d(j)%N_psd%sigma0/sqrt(box_width)))**2 &
                  &)
               end do
            end do
         end do
         call mpi_reduce(chisq_scan, chisq_new, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cpar%root, MPI_COMM_WORLD, ierr)

         ! Use chisq from the first iteration where we dont draw new parameters as the base chisq
         if (first_prop) then 
            chisq_prev = chisq_new
         else 
            if (is_root) then
               chisq_diff = max(chisq_new - chisq_prev, 0.)
               ln_accept_prob = -0.5 * chisq_diff
               select case (cpar%operation)
                  case ("sample")
                     accepted = ln_accept_prob > log(rand_uni(handle))
                  case ("optimize")
                     accepted = chisq_diff < 0.
                  case default
                     stop "Error: invalid cpar%operation in comm_zodi_samp_mod sample_n0 routine"
               end select
            end if
            call mpi_bcast(accepted, 1, MPI_LOGICAL, cpar%root, cpar%comm_chain, ierr)

            if (accepted) then
               n_accepted = n_accepted + 1
               chisq_prev = chisq_new
               n0_prev = n0_new
               emissivity_prev = emissivity_new
               albedo_prev = albedo_new
               if (is_root) then
                  print *, "proposal:", prop
                  write(*, '(a, ES12.5)') "chisq: ", chisq_new
                  print *, " "
               end if
            else
               n0_new = n0_prev
               emissivity_new = emissivity_prev
               albedo_new = albedo_prev
            end if
            accept_rate = real(n_accepted)/real(prop)
         end if
      end do
      if (tuning) then
         if (accept_rate < 0.1) then
            step_sizes_n0 = step_sizes_n0 / 2.
            step_sizes_emissivity = step_sizes_emissivity / 2.
            step_sizes_albedo = step_sizes_albedo / 2.
         else if (accept_rate > 0.4) then
            step_sizes_n0 = step_sizes_n0 * 2.
            step_sizes_emissivity = step_sizes_emissivity * 2.
            step_sizes_albedo = step_sizes_albedo * 2.
         end if
      end if
      if (is_root) print *, "n0 accept rate: ", (real(n_accepted)/real(nprop))*100.
      do i = 1, model%n_comps
         model%comps(i)%c%n_0 = n0_prev(i)
      end do
      do i = 1, numband
         data(i)%tod%zodi_albedo = albedo_prev(i, :)
         data(i)%tod%zodi_emissivity = emissivity_prev(i, :)
      end do
   end subroutine

   subroutine sample_zodi_group(cpar, handle, gibbs_iter, model, verbose)
      type(comm_params), intent(in) :: cpar
      type(planck_rng), intent(inout) :: handle
      integer(i4b), intent(in) :: gibbs_iter
      type(ZodiModel), intent(inout) :: model
      logical(lgt), intent(in), optional :: verbose

      integer(i4b) :: i, j, k, prop, group_idx, flag, ndet, nscan, ntod, n_proposals, scan, ierr, n_accepted, param_idx
      real(dp) :: chisq_tod, chisq_current, chisq_prior, chisq_diff, ln_acceptance_probability, accept_rate, chisq_lnL, box_width
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

      tuning = .true.
      n_proposals = 500

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
                     param_vec(param_idx) = param_vec(param_idx) + (step_sizes_ipd(param_idx) * rand_gauss(handle))
                  end do
                  chisq_prior = get_chisq_priors(param_vec, prior_vec)
                  if (chisq_prior >= 1.d30) skip = .true.
               end if
               call mpi_bcast(chisq_prior, 1, MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
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

                     if (chisq_tod >= 1.d30) exit
                  end do
               end do
            end do

            ! Reduce chisq to root process
            call mpi_reduce(chisq_tod, chisq_current, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cpar%root, MPI_COMM_WORLD, ierr)

            ! Add prior penalty to chisq
            if (prop > 0 ) chisq_current = chisq_current + chisq_prior

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
                     write(*, '(a, ES12.5)') "chisq: ", chisq_current
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
         if (accept_rate < 0.05) then
            do i = 1, size(group_indices)
               param_idx = group_indices(i)
               step_sizes_ipd(param_idx) = step_sizes_ipd(param_idx) / 2.
            end do
         else if (accept_rate > 0.95) then
            do i = 1, size(group_indices)
               param_idx = group_indices(i)
               step_sizes_ipd(param_idx) = step_sizes_ipd(param_idx) * 2.
            end do
         end if
         model = current_model
      end do
   end subroutine

   function get_chisq_priors(params, priors) result(chisq_prior)
      real(dp), intent(in) :: params(:)
      real(dp), intent(in) :: priors(:, :)
      real(dp) :: chisq_prior
      integer(i4b) :: i
      logical(lgt) :: prior_is_violated

      chisq_prior = 0.
      do i = 1, size(params)
         if (priors(i, 3) == 0) then
            prior_is_violated = params(i) <= priors(i, 1) .or. params(i) >= priors(i, 2)
            if (prior_is_violated) then
               chisq_prior = 1.d30
               return
            end if
         else if (priors(i, 3) == 1) then
            chisq_prior = chisq_prior + (params(i) - priors(i, 1))**2 / priors(i, 2)**2
         end if
      end do
   end function 

   subroutine downsamp_invariant_structs(cpar)
      ! Downsamples pointing, tod and caches the zodi mask (sky mask + flags)
      type(comm_params), intent(in) :: cpar

      integer(i4b) :: i, j, k, l, g, scan, npix, nmaps, ndelta, ext(2), padding, ierr, ntod, ndet, nhorn, ndownsamp, box_halfwidth
      real(dp) :: box_width, dt_tod
      real(sp), allocatable :: tod(:), mask(:), downsamp_vec(:, :)
      integer(i4b), allocatable :: pix(:, :), psi(:, :), flag(:)
      real(dp), allocatable, dimension(:, :) :: m_buf, vec
      integer(i4b), allocatable, dimension(:) :: downsamp_pix
      logical(lgt), allocatable, dimension(:) :: downsamp_mask_idx
      real(sp), allocatable, dimension(:) :: downsamp_mask, downsamp_tod, downsamp_obs_time, obs_time
      real(sp), allocatable, dimension(:) :: procmask_zodi
      type(hdf_file) :: tod_file
      character(len=4) :: scan_str

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
         box_halfwidth = anint(box_width / 2.)
         do scan = 1, data(i)%tod%nscan
            if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
            ntod = data(i)%tod%scans(scan)%ntod
            ! Allocate temporary variables to be downsampled
            allocate (pix(ntod, nhorn))
            allocate (psi(ntod, nhorn))
            allocate (flag(ntod))
            allocate (tod(ntod))
            allocate (mask(ntod))
            allocate (vec(3, ntod))

            ! Get shape of downsampled structs
            call data(i)%tod%downsample_tod(tod, ext, step=box_width)

            ! Allocate downsampled strucs
            allocate (downsamp_tod(ext(1):ext(2)))
            allocate (downsamp_mask(ext(1):ext(2)))
            allocate (downsamp_pix(ext(1):ext(2)))
            allocate (downsamp_mask_idx(ext(1):ext(2)))
            allocate (downsamp_vec(3, ext(1):ext(2)))

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
               downsamp_vec = 0. 

               call data(i)%tod%decompress_pointing_and_flags(scan, j, pix, psi, flag)

               if (data(i)%tod%compressed_tod) then
                  call data(i)%tod%decompress_tod(scan, j, tod)
               else
                  tod = data(i)%tod%scans(scan)%d(j)%tod
               end if

               do k = 1, data(i)%tod%scans(scan)%ntod
                  mask(k) = procmask_zodi(pix(k, 1))
                  if (iand(flag(k), data(i)%tod%flag0) .ne. 0) mask(k) = 0.
                  vec(:, k) = data(i)%tod%ind2vec(:, data(i)%tod%pix2ind(pix(k, 1)))
               end do

               where (mask > 0.5) 
                  mask = 1.
               elsewhere
                  mask = 0.
               end where

               do k = 1, 3
                  call data(i)%tod%downsample_tod(real(vec(k, :), sp), ext, downsamp_vec(k, :), step=box_width)
               end do

               do k = 0, ext(2)-padding
                  call vec2pix_ring(data(i)%tod%nside, real(downsamp_vec(:, k), dp), downsamp_pix(k))
               end do

               ! Downsample the mask
               call data(i)%tod%downsample_tod(mask, ext, downsamp_mask, step=box_width)

               ! Downsample the tods
               call data(i)%tod%downsample_tod(tod, ext, downsamp_tod, step=box_width)

               ! Construct the downsampled pix array (pick central pixel in each box)

               ! do k = 1, int(size(tod)/box_width) + 1
               !    l = floor(min(((k-1)*box_width + box_halfwidth), real(ntod, dp)))
               !    downsamp_pix(k) = pix(l, 1)
               ! end do

               where (downsamp_mask < 1.) downsamp_mask_idx = .false.

               ! Find upper bound and truncate the downsampled arrays by where they were padded and
               ! filter out the masked values. Store the result in the tod objects for use in sampling.
               data(i)%tod%scans(scan)%d(j)%downsamp_pix = pack(downsamp_pix(0:ext(2)-padding), downsamp_mask_idx(0:ext(2)-padding))
               data(i)%tod%scans(scan)%d(j)%downsamp_tod = pack(downsamp_tod(0:ext(2)-padding), downsamp_mask_idx(0:ext(2)-padding))
               if (j == 1) data(i)%tod%scans(scan)%downsamp_obs_time = pack(downsamp_obs_time(0:ext(2)-padding), downsamp_mask_idx(0:ext(2)-padding))
               ! write timestreams to files
               ! call int2string(data(i)%tod%scanid(scan), scan_str)
               ! call open_hdf_file(trim(adjustl("/mn/stornext/u3/metins/dirbe/chains/chains_downsamp/dtod_"//scan_str//".h5")), tod_file, 'w')
               ! call write_hdf(tod_file, '/dtod', data(i)%tod%scans(scan)%d(j)%downsamp_tod)
               ! call write_hdf(tod_file, '/vec', vec)
               ! call write_hdf(tod_file, '/dvec', downsamp_vec)
               ! call write_hdf(tod_file, '/pix', pix(:, 1))
               ! call write_hdf(tod_file, '/dpix', downsamp_pix)
               ! call close_hdf_file(tod_file)

               ! Allocate other downsampled quantities with same shape
               ndownsamp = size(data(i)%tod%scans(scan)%d(j)%downsamp_pix)
               allocate (data(i)%tod%scans(scan)%d(j)%downsamp_zodi(ndownsamp))
               allocate (data(i)%tod%scans(scan)%d(j)%downsamp_scat(ndownsamp, zodi_model%n_comps))
               allocate (data(i)%tod%scans(scan)%d(j)%downsamp_therm(ndownsamp, zodi_model%n_comps))
            end do
            deallocate (obs_time, downsamp_obs_time)
            deallocate (pix, psi, vec, flag, tod, mask, downsamp_tod, downsamp_mask, downsamp_pix, downsamp_mask_idx, downsamp_vec)
         end do

         deallocate (m_buf, procmask_zodi)
      end do
   end subroutine
   

   subroutine precompute_lowres_zodi_lookups(cpar)
      ! Loop over each band with zodi and precompute lookup tables for lowres zodi caching
      type(comm_params), intent(in) :: cpar
      integer(i4b) :: i, j, k, l, scan, n_lowres_obs, nobs_downsamp, nobs_lowres, pix_high, pix_low, ierr
      integer(i4b), allocatable :: pix2ind_highres(:), ind2pix_highres(:)
      real(dp), allocatable :: ind2vec_zodi_temp(:, :)
      real(dp) :: rotation_matrix(3, 3)

      call ecl_to_gal_rot_mat(rotation_matrix)

      do i = 1, numband
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         allocate(pix2ind_highres(0:12*data(i)%tod%nside**2-1))         
         pix2ind_highres = 0
         
         do scan = 1, data(i)%tod%nscan
            if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
            do j = 1, data(i)%tod%ndet
               do k = 1, size(data(i)%tod%scans(scan)%d(j)%downsamp_pix)
                  pix2ind_highres(data(i)%tod%scans(scan)%d(j)%downsamp_pix(k)) = 1
               end do
            end do
         end do
         
         nobs_downsamp = count(pix2ind_highres == 1)
         
         allocate(ind2vec_zodi_temp(3, nobs_downsamp))
         allocate(ind2pix_highres(nobs_downsamp))

         j = 1
         do k = 0, 12*data(i)%tod%nside**2-1
            if (pix2ind_highres(k) == 1) then
               ind2pix_highres(j) = k
               j = j+1
            end if
         end do

         allocate(data(i)%tod%pix2ind_lowres(0:12*zodi_nside**2-1))
         data(i)%tod%pix2ind_lowres = 0
         ind2vec_zodi_temp = 0.

         j = 1
         do k = 1, nobs_downsamp
            pix_high = ind2pix_highres(k)
            pix_low = data(i)%tod%udgrade_pix_zodi(pix_high)
            if (data(i)%tod%pix2ind_lowres(pix_low) == 0) then
               data(i)%tod%pix2ind_lowres(pix_low) = j
               call pix2vec_ring(zodi_nside, pix_low, ind2vec_zodi_temp(:, j))
            end if
            j =  j + 1
         end do

         nobs_lowres = j - 1
         allocate(data(i)%tod%ind2vec_ecl_lowres(3, nobs_lowres))
         data(i)%tod%ind2vec_ecl_lowres = ind2vec_zodi_temp(:, 1:nobs_lowres)
         
         allocate(data(i)%tod%zodi_scat_cache_lowres(nobs_lowres, data(i)%tod%zodi_n_comps, data(i)%tod%ndet))
         allocate(data(i)%tod%zodi_therm_cache_lowres(nobs_lowres, data(i)%tod%zodi_n_comps, data(i)%tod%ndet))
         data(i)%tod%zodi_scat_cache_lowres = -1.d0
         data(i)%tod%zodi_therm_cache_lowres = -1.d0
         
         do k = 1, nobs_lowres
            data(i)%tod%ind2vec_ecl_lowres(:, k) = matmul(data(i)%tod%ind2vec_ecl_lowres(:, k), rotation_matrix)
         end do   
         deallocate(ind2vec_zodi_temp, pix2ind_highres, ind2pix_highres)
      end do
   end subroutine

   subroutine project_and_downsamp_sky(cpar)
      type(comm_params), intent(in) :: cpar
      integer(i4b) :: i, j, k, scan, ext(2), upper_bound, padding, ierr, ntod, nhorn, npix, ndet, nmaps, ndelta
      real(dp) :: box_width
      type(map_ptr), allocatable, dimension(:, :) :: sky_signal
      real(sp), allocatable, dimension(:) :: downsamp_sky, downsamp_mask
      logical(lgt), allocatable, dimension(:) :: downsamp_mask_idx
      real(sp), allocatable, dimension(:, :, :, :) :: map_sky
      real(dp), allocatable, dimension(:, :) :: m_buf
      type(hdf_file) :: tod_file
      integer(i4b), allocatable :: pix(:, :), psi(:, :), flag(:)
      real(sp), allocatable :: mask(:), sky(:)
      real(sp), allocatable, dimension(:) :: procmask_zodi

      if (cpar%myid == cpar%root) print *, "projecting downsampled sky"

      padding = 5
      ! For each zodi band, create downsampled residual time stream
      do i = 1, numband
         ! Only generate downsampled arrays for tod bands and bands where we have zodi
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         npix = 12*data(i)%map%info%nside**2
         ndelta = 1
         nmaps = data(i)%map%info%nmaps
         ndet = data(i)%tod%ndet
         nhorn = data(i)%tod%nhorn
         box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)

         allocate (sky_signal(data(i)%tod%ndet, ndelta))
         do j = 1, data(i)%tod%ndet
            call get_sky_signal(i, j, sky_signal(j, ndelta)%p, mono=.true.)
         end do

         allocate (map_sky(nmaps, data(i)%tod%nobs, 0:data(i)%tod%ndet, 1))
         call distribute_sky_maps(data(i)%tod, sky_signal, 1.e0, map_sky)

         allocate (m_buf(0:npix - 1, nmaps))
         allocate (procmask_zodi(0:npix - 1))
         call data(i)%tod%procmask_zodi%bcast_fullsky_map(m_buf)
         procmask_zodi = m_buf(:, 1)

         do scan = 1, data(i)%tod%nscan
            if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
            ntod = data(i)%tod%scans(scan)%ntod

            do j = 1, ndet
               allocate (pix(ntod, nhorn))
               allocate (psi(ntod, nhorn))
               allocate (flag(ntod))
               allocate (mask(ntod))
               allocate (sky(ntod))
               pix = 0
               psi = 0
               flag = 0
               mask = 0.
               sky = 0.

               ! Project sky map down to downsampled timestream and cache
               call data(i)%tod%decompress_pointing_and_flags(scan, j, pix, psi, flag)

               do k = 1, ntod
                  mask(k) = procmask_zodi(pix(k, 1))
                  if (iand(flag(k), data(i)%tod%flag0) .ne. 0) mask(k) = 0.
               end do
               where (mask > 0.) mask = 1. ! make sure mask is binary before downsampling
               
               do k = 1, ntod
                  sky(k) = map_sky(1, data(i)%tod%pix2ind(pix(k, 1)), j, 1)  !zodi is only temperature (for now)
               end do
               
               ! Get downsample shape
               call data(i)%tod%downsample_tod(sky, ext, step=box_width)

               allocate (downsamp_sky(ext(1):ext(2)))
               allocate (downsamp_mask(ext(1):ext(2)))
               allocate (downsamp_mask_idx(ext(1):ext(2)))
               downsamp_mask_idx = .true.

               call data(i)%tod%downsample_tod(sky, ext, downsamp_sky, step=box_width)

               call data(i)%tod%downsample_tod(mask, ext, downsamp_mask, step=box_width)
               where (downsamp_mask < 1.) downsamp_mask_idx = .false.

               data(i)%tod%scans(scan)%d(j)%downsamp_sky = pack(downsamp_sky(0:ext(2)-padding), downsamp_mask_idx(0:ext(2)-padding))
               deallocate(pix, psi, flag, mask, sky)
               deallocate(downsamp_sky, downsamp_mask_idx, downsamp_mask)
            end do

         end do
         deallocate (map_sky, sky_signal, m_buf, procmask_zodi)
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

   subroutine remove_glitches_from_downsamped_zodi_quantities(cpar)
      type(comm_params), intent(in) :: cpar
      integer(i4b) :: i, j, k, scan, ierr, non_glitch_size, thinstep
      real(dp) :: box_width, rms
      real(sp), allocatable :: res(:)
      logical(lgt), allocatable :: glitch_mask(:)
      real(sp), allocatable :: downsamp_scat_comp(:, :), downsamp_therm_comp(:, :)

      thinstep = nint(1.d0/cpar%zs_tod_thin_factor)
      do i = 1, numband
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         do scan = 1, data(i)%tod%nscan
            box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)
            do j = 1, data(i)%tod%ndet
               if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
               
               ! Search for strong outliers
               res = data(i)%tod%scans(scan)%d(j)%downsamp_tod - data(i)%tod%scans(scan)%d(j)%downsamp_sky - data(i)%tod%scans(scan)%d(j)%downsamp_zodi
               !write(*,*) 'a', sqrt(mean(real(data(i)%tod%scans(scan)%d(j)%downsamp_tod**2, dp)))
               !write(*,*) 'b', sqrt(mean(real(data(i)%tod%scans(scan)%d(j)%downsamp_sky**2, dp)))
               !write(*,*) 'c', sqrt(mean(real(data(i)%tod%scans(scan)%d(j)%downsamp_zodi**2, dp)))
               rms = sqrt(mean(real(res**2, dp)))
               glitch_mask = abs(res) > 5. * real(rms, sp)

               ! Apply TOD thinning
               do k = 1, size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                  !if (mod(k,thinstep) /= 0 .and. abs(res(k))/(data(i)%tod%scans(scan)%d(j)%N_psd%sigma0/sqrt(box_width)) < cpar%zs_tod_thin_threshold) then
                  if (mod(k,thinstep) /= 0) then
                     glitch_mask(k) = .true.
                  end if
               end do
              
               ! Remove flagged samples
               data(i)%tod%scans(scan)%d(j)%downsamp_tod = pack(data(i)%tod%scans(scan)%d(j)%downsamp_tod, .not. glitch_mask)
               data(i)%tod%scans(scan)%d(j)%downsamp_sky = pack(data(i)%tod%scans(scan)%d(j)%downsamp_sky, .not. glitch_mask)
               data(i)%tod%scans(scan)%d(j)%downsamp_zodi = pack(data(i)%tod%scans(scan)%d(j)%downsamp_zodi, .not. glitch_mask)
               data(i)%tod%scans(scan)%d(j)%downsamp_pix = pack(data(i)%tod%scans(scan)%d(j)%downsamp_pix, .not. glitch_mask)
               non_glitch_size = count(.not. glitch_mask)

               ! pack doesnt work on multidimensional arrays so here we manually reallocate the zodi caches to the new sizes
               allocate(downsamp_scat_comp(non_glitch_size, zodi_model%n_comps))
               allocate(downsamp_therm_comp(non_glitch_size, zodi_model%n_comps))
               do k = 1, zodi_model%n_comps
                  downsamp_scat_comp(:, k) = pack(data(i)%tod%scans(scan)%d(j)%downsamp_scat(:, k), .not. glitch_mask)
                  downsamp_therm_comp(:, k) = pack(data(i)%tod%scans(scan)%d(j)%downsamp_therm(:, k), .not. glitch_mask)
               end do
               deallocate(data(i)%tod%scans(scan)%d(j)%downsamp_scat)
               deallocate(data(i)%tod%scans(scan)%d(j)%downsamp_therm)
               allocate(data(i)%tod%scans(scan)%d(j)%downsamp_scat(non_glitch_size, zodi_model%n_comps))
               allocate(data(i)%tod%scans(scan)%d(j)%downsamp_therm(non_glitch_size, zodi_model%n_comps))

               data(i)%tod%scans(scan)%d(j)%downsamp_therm = downsamp_therm_comp
               data(i)%tod%scans(scan)%d(j)%downsamp_scat = downsamp_scat_comp
               deallocate(downsamp_scat_comp, downsamp_therm_comp)
            end do
         end do
      end do
   end subroutine

   subroutine minimize_zodi_with_powell(cpar)
      type(comm_params), intent(in) :: cpar
      logical(lgt), save :: first_call = .true.
      real(dp), allocatable :: theta(:), theta_phys(:), theta_new(:), theta_final(:)
      character(len=128), allocatable :: labels(:)
      integer(i4b) :: i, j, ierr, flag, group_idx, end_idx, n_bands
      real(dp) :: chisq_lnL
      real(dp), allocatable, dimension(:) :: prior_vec_powell_min, prior_vec_powell_max, prior_vec_powell_type 
      integer(i4b), allocatable :: indices(:)

      allocate(theta(zodi_model%n_params))
      call zodi_model%model_to_params(theta, labels)
      if (cpar%myid == 0) print *, "minimizing zodi parameters with powell"
      
      ! filter out parameters for powell search
      allocate(powell_included_params(zodi_model%n_params))
      powell_included_params = .false.

      ! Get parameters to sample from sampling groups
      do group_idx = 1, cpar%zs_num_samp_groups
         call parse_samp_group_strings(cpar%zs_samp_groups(group_idx), labels, indices)
         do i = 1, size(indices)
            powell_included_params(indices(i)) = .true.
         end do
      end do
      prior_vec_powell_min = prior_vec(:, 1)
      prior_vec_powell_max = prior_vec(:, 2)
      prior_vec_powell_type = prior_vec(:, 3)

      ! ! append emissivities and albedoes to the vector given to powell.
      do i = 1, numband
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle
         if (cpar%ds_zodi_reference_band(data(i)%id_abs)) cycle
         theta = [theta, data(i)%tod%zodi_emissivity]
         powell_included_params = [powell_included_params, [(.true. , j=1, size(data(i)%tod%zodi_emissivity))]]
         prior_vec_powell_min = [prior_vec_powell_min, [(emissivity_prior(1) , j=1, size(data(i)%tod%zodi_emissivity))]]
         prior_vec_powell_max = [prior_vec_powell_max, [(emissivity_prior(2) , j=1, size(data(i)%tod%zodi_emissivity))]]
         prior_vec_powell_type = [prior_vec_powell_type, [(real(0, dp), j=1, size(data(i)%tod%zodi_emissivity))]]
      end do
      
      do i = 1, numband
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle
         if (.not. any(data(i)%tod%zodi_albedo > EPS)) cycle
         theta = [theta, data(i)%tod%zodi_albedo]
         powell_included_params = [powell_included_params, [(.true. , j=1, size(data(i)%tod%zodi_albedo))]]
         prior_vec_powell_min = [prior_vec_powell_min, [(albedo_prior(1) , j=1, size(data(i)%tod%zodi_albedo))]]
         prior_vec_powell_max = [prior_vec_powell_max, [(albedo_prior(2) , j=1, size(data(i)%tod%zodi_albedo))]]
         prior_vec_powell_type = [prior_vec_powell_type, [(real(0, dp), j=1, size(data(i)%tod%zodi_albedo))]]
      end do
      
      allocate(prior_vec_powell(size(prior_vec_powell_min), 3))
      prior_vec_powell(:, 1) = prior_vec_powell_min
      prior_vec_powell(:, 2) = prior_vec_powell_max
      prior_vec_powell(:, 3) = prior_vec_powell_type

      if (any(pack(theta, powell_included_params) == 0.)) then
         do i = 1, size(theta)
            write(*,*) i, theta(i)
         end do
         stop "theta_0 contains zeros in zodi powell sampling. Cant compute physical units due to theta/theta_0"
      end if

      allocate(theta_0(size(theta)))
      theta_0 = theta

      allocate(theta_phys(count(powell_included_params)))
      
      if (cpar%myid == cpar%root) then
         !filter out N_0 parameters and scale to physical units
         theta_phys = pack(theta / theta_0, powell_included_params)
         !call powell(theta_phys, lnL_zodi, ierr)
         call powell(theta_phys, lnL_zodi, ierr, tolerance=1d-3)
         flag = 0
         call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
      else
         do while (.true.)
             call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
             if (flag == 1) then 
                 chisq_lnL = lnL_zodi()
             else
                 exit
             end if
         end do
      end if

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_bcast(theta_phys, size(theta_phys), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
      
      allocate(theta_final(zodi_model%n_params))
      theta_final = theta(1:zodi_model%n_params)

      j = 1
      do i = 1, size(theta)
         if (powell_included_params(i)) then
            if (i <= zodi_model%n_params) then 
               theta_final(i) = theta_phys(j) * theta_0(i)
            end if
            theta(i) = theta_phys(j) * theta_0(i)
            j = j + 1
         end if
      end do

      ! update model with final parameters
      call zodi_model%params_to_model(theta_final)

      ! update emissivities and albedos with final parameters
      j = 0
      do i = numband, 1, -1
         if (data(i)%tod_type == "none") cycle
         if (.not. data(i)%tod%sample_zodi) cycle
         if (any(data(i)%tod%zodi_albedo > EPS)) then
            j = j + 1
            data(i)%tod%zodi_albedo = theta(size(theta)-(j*zodi_model%n_comps)+1:size(theta)-(j-1)*zodi_model%n_comps)
         else 
            data(i)%tod%zodi_albedo = data(i)%tod%zodi_albedo
         end if 
      end do
      do i = numband, 1, -1
         if (data(i)%tod_type == "none") cycle
         if (.not. data(i)%tod%sample_zodi) cycle
         if (.not. ref_band(data(i)%id_abs)) then
            j = j + 1
            data(i)%tod%zodi_emissivity = theta(size(theta)-(j*zodi_model%n_comps)+1:size(theta)-(j-1)*zodi_model%n_comps)
         else
            data(i)%tod%zodi_emissivity = 1.
         end if
      end do
      deallocate(prior_vec_powell, theta_0, powell_included_params)
   end subroutine

   function lnL_zodi(p)
      use healpix_types
      implicit none
      real(dp), dimension(:), intent(in), optional :: p
      real(dp) :: lnL_zodi
      type(ZodiModel) :: model

      real(dp), allocatable :: theta(:), theta_phys(:), theta_full(:)
      character(len=128), allocatable :: labels(:)
      real(dp) :: chisq, chisq_tod, chisq_prior, box_width, t1, t2, t3, t4
      integer(i4b) :: i, j, k, scan, ntod, ndet, nscan, flag, ierr, n_bands, ndof, ndof_tot
      character(len=4) :: scan_str
      type(hdf_file) :: tod_file

      call wall_time(t1)
      model = zodi_model
      allocate(theta_phys(count(powell_included_params)))
      
      if (data(1)%tod%myid == 0) then
         flag = 1
         call mpi_bcast(flag, 1, MPI_INTEGER, 0, data(1)%tod%comm, ierr)
         theta_phys = p
      end if
      call mpi_bcast(theta_phys, size(theta_phys), MPI_DOUBLE_PRECISION, 0, data(1)%tod%comm, ierr)
      

      allocate(theta(model%n_params))
      call model%model_to_params(theta, labels)

      allocate(theta_full(size(powell_included_params)))
      theta_full(1:model%n_params) = theta
      
      j = 1
      do i = 1, size(theta_full)
         if (powell_included_params(i)) then
            if (i <= model%n_params) then 
               theta(i) = theta_phys(j) * theta_0(i)
            end if
            theta_full(i) = theta_phys(j) * theta_0(i)
            j = j + 1
         end if
      end do

      call model%params_to_model(theta)

      ! if (data(1)%tod%myid == 0) call print_zodi_model(theta, labels)


      j = 0
      do i = numband, 1, -1
         if (data(i)%tod_type == "none") cycle
         if (.not. data(i)%tod%sample_zodi) cycle
         if (any(data(i)%tod%zodi_albedo > EPS)) then
            j = j + 1
            powell_albedo(i, :) = theta_full(size(theta_full)-(j*zodi_model%n_comps)+1:size(theta_full)-(j-1)*zodi_model%n_comps)
         else 
            powell_albedo(i, :) = data(i)%tod%zodi_albedo
         end if 
      end do
      do i = numband, 1, -1
         if (data(i)%tod_type == "none") cycle
         if (.not. data(i)%tod%sample_zodi) cycle
         if (.not. ref_band(data(i)%id_abs)) then
            j = j + 1
            powell_emissivity(i, :) = theta_full(size(theta_full)-(j*zodi_model%n_comps)+1:size(theta_full)-(j-1)*zodi_model%n_comps)
         else
            powell_emissivity(i, :) = 1.
         end if
      end do

      ! rescale parameters 
      call model%model_to_params(theta)
      theta_full(1:model%n_params) = theta
      
      if (data(1)%tod%myid == 0) then
         call print_zodi_model(theta, labels, chisq)
      end if

      chisq_tod = get_chisq_priors(theta_full, prior_vec_powell)
      chisq = 0.
      lnL_zodi = 0.
      ndof = 0

      do i = 1, numband
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

               call wall_time(t3)
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
               call wall_time(t4)
               !if (data(1)%tod%myid == 0) write(*,*) ' CPU1 = ', t4-t3
               call get_s_zodi(&
                   & s_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                   & s_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                   & s_zodi=data(i)%tod%scans(scan)%d(j)%downsamp_zodi, &
                   & emissivity=powell_emissivity(i, :), &
                   & albedo=powell_albedo(i, :) &
                   &)
               call wall_time(t3)
               !if (data(1)%tod%myid == 0) write(*,*) ' CPU2 = ', t3-t4
               
               chisq = sum( &
                  & ((data(i)%tod%scans(scan)%d(j)%downsamp_tod &
                  &   - data(i)%tod%scans(scan)%d(j)%downsamp_sky &
                  &   - data(i)%tod%scans(scan)%d(j)%downsamp_zodi &
                  & )/(data(i)%tod%scans(scan)%d(j)%N_psd%sigma0/sqrt(box_width)))**2 &
                  &)
               chisq_tod = chisq_tod + chisq
               call wall_time(t4)
               !if (data(1)%tod%myid == 0) write(*,*) ' CPU3 = ', t4-t3


               
               ndof = ndof + size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
               if (chisq_tod >= 1.d30) exit
               ! call int2string(data(i)%tod%scanid(scan), scan_str)
               ! call open_hdf_file(trim(adjustl("/mn/stornext/u3/metins/dirbe/chains/chains_downsamp/dtodlnl_"//scan_str//".h5")), tod_file, 'w')
               ! call write_hdf(tod_file, '/dtod', data(i)%tod%scans(scan)%d(j)%downsamp_tod)
               ! call write_hdf(tod_file, '/dzodi', data(i)%tod%scans(scan)%d(j)%downsamp_zodi)
               ! call write_hdf(tod_file, '/dsky', data(i)%tod%scans(scan)%d(j)%downsamp_sky)
               ! call write_hdf(tod_file, '/dpix', data(i)%tod%scans(scan)%d(j)%downsamp_pix)
               ! call close_hdf_file(tod_file)

               if (data(1)%tod%myid == 0 .and. scan == 1) then
                  !write(*,*) "scan = ", data(i)%tod%scanid(scan), sum(abs(data(i)%tod%scans(scan)%d(j)%downsamp_tod)), sum(abs(data(i)%tod%scans(scan)%d(j)%downsamp_sky)), sum(abs(data(i)%tod%scans(scan)%d(j)%downsamp_zodi)), data(i)%tod%scans(scan)%d(j)%N_psd%sigma0
                  open(58,file='res.dat')
                  do k = 1, size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                     write(58,*) data(i)%tod%scans(scan)%d(j)%downsamp_tod(k) &
                          &   - data(i)%tod%scans(scan)%d(j)%downsamp_sky(k) &
                          &   - data(i)%tod%scans(scan)%d(j)%downsamp_zodi(k)
                  end do
                  close(58)
               end if
               call wall_time(t3)
               !if (data(1)%tod%myid == 0) write(*,*) ' CPU4 = ', t3-t4

            end do
         end do
      end do
      ! call mpi_barrier(MPI_COMM_WORLD, ierr)
      ! stop

      ! Reduce chisq to root process
      call mpi_reduce(chisq_tod, chisq,    1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, data(1)%tod%comm, ierr)
      call mpi_reduce(ndof,      ndof_tot, 1, MPI_INTEGER, MPI_SUM, 0, data(1)%tod%comm, ierr)

      call wall_time(t4)
      !if (data(1)%tod%myid == 0) write(*,*) ' CPU5 = ', t4-t3

      if (data(1)%tod%myid == 0) then
         lnL_zodi = chisq
         call wall_time(t2)
         if (ndof_tot > 0) write(*,fmt='(a,e16.8,a,f10.4,a,f8.3)') "chisq_zodi = ", chisq, ", chisq_red = ", chisq/ndof_tot, ", time = ", t2-t1
         write(*,*)
         ! call print_zodi_model(theta, labels, chisq)
      end if
   end function

   subroutine get_s_zodi_with_n0(s_therm, s_scat, s_zodi, emissivity, albedo, n_0_comp_ratio, comp)
      real(sp), dimension(:, :), intent(in) :: s_scat, s_therm
      real(sp), dimension(:), intent(inout) :: s_zodi
      real(dp), dimension(:), intent(in) :: emissivity, albedo, n_0_comp_ratio
      integer(i4b) :: comp
      integer(i4b) :: i

      s_zodi = 0.
      do i = 1, size(s_therm, dim= 2)
         s_zodi = s_zodi + (s_scat(:, i) * albedo(i) + (1. - albedo(i)) * emissivity(i) * s_therm(:, i)) * n_0_comp_ratio(i)
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


   subroutine zodi_model_to_ascii(cpar, model, filename, overwrite)
      ! Dumps the zodi model to an ascii file on the format {COMP}_{PARAM} = {VALUE}.
      class(ZodiModel), target, intent(in) :: model
      type(comm_params), intent(in) :: cpar
      character(len=*), intent(in) :: filename
      logical(lgt), intent(in), optional :: overwrite

      integer(i4b) :: io, i, j, running_idx
      logical(lgt) :: exists, overwrite_
      real(dp), allocatable :: params(:)
      integer(i4b), allocatable :: comp_switch_indices(:)
      character(len=128), allocatable :: labels(:)
      character(len=512) :: concatenated_string, val

      if (present(overwrite)) then
         overwrite_ = overwrite
      else
         overwrite_ = .false.
      end if

      if (cpar%myid_chain /= cpar%root) return
      inquire(file=trim(adjustl(filename)), exist=exists)
      if (exists .and. (.not. overwrite_)) then
         print *, "zodi asciifile: " // trim(adjustl(filename)) // " exists and overwrite = .false."
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


   subroutine print_zodi_model(theta, labels, chisq)
      real(dp), allocatable, intent(in) :: theta(:)
      character(len=128), allocatable, intent(in) :: labels(:)
      real(dp), intent(in), optional :: chisq
      integer(i4b) :: i, j, k, l, idx, n_cols, comp, n_params
      
      n_params = size(theta)
      n_cols = 5
      comp = 1
      do i = 1, n_params
         idx = index(trim(adjustl(labels(i))), "_")
         if (i < n_params) then
            if (i == 1) then 
               write(*, "(a, a)", advance="no") trim(adjustl(labels(i)(:idx-1))), ":  "
            else if (trim(adjustl(labels(i)(:idx))) /= trim(adjustl(labels(i-1)(:idx)))) then
               comp = comp + 1
               ! write(*, *), ""
               if (comp > zodi_model%n_comps) then
                  write(*, "(a)", advance="no") "OTHER:  "
               else 
                  write(*, "(a, a)", advance="no") trim(adjustl(labels(i)(:idx-1))), ":  "
               end if
            end if
         end if
         if (mod(i, n_cols) == 0 .or. (i<n_params .and. (trim(adjustl(labels(i)(:idx))) /= trim(adjustl(labels(i+1)(:idx)))))) then
            write(*, "(a,a,g0.4,a)") trim(adjustl(labels(i)(idx+1:))), "=", theta(i), ",  "
         else
            write(*, "(a,a,g0.4,a)", advance="no") trim(adjustl(labels(i)(idx+1:))), "=", theta(i), ",  "
         end if
      end do
      print *, ""
   end subroutine
 
end module
