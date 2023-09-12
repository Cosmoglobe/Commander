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
    public initialize_zodi_samp_mod, downsamp_invariant_structs, project_and_downsamp_sky, compute_downsamp_zodi, sample_zodi_emissivity_and_albedo, active_params, n_active_geom_params
    ! public initialize_zodi_samp_mod, sample_zodi_emissivity_and_albedo, sample_zodi_parameter, downsamp_invariant_structs, project_and_downsamp_sky, cleanup_zodi, active_params, n_active_geom_params

    ! globals
    real(dp), allocatable :: chisq_previous, step_size, step_sizes(:), priors(:, :), step_sizes_emissivity(:, :), step_sizes_albedo(:, :)
    logical(lgt), allocatable :: active_params(:)
    integer(i4b) :: n_samp_bands, n_active_params, n_active_geom_params, cloud_start, cloud_stop, band1_start, band1_stop, band2_start, band2_stop, band3_start, band3_stop, ring_start, ring_stop, feature_start, feature_stop, emissivity_start, emissivity_stop, albedo_start, albedo_stop
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

        type(comm_params), intent(in) :: cpar
        integer(i4b) :: i, j
        

        ! Figure out how many sampling bands there are and initialize the tod step sizes
        n_samp_bands = 0
        do i = 1, numband
            if (data(i)%tod_type == 'none') cycle
            if (.not. data(i)%tod%subtract_zodi) cycle
            n_samp_bands = n_samp_bands + 1
        end do

        allocate(step_sizes_albedo(n_samp_bands, zodi_model%n_comps))
        allocate(step_sizes_emissivity(n_samp_bands, zodi_model%n_comps))
        j = 0
        do i = 1, numband
            if (data(i)%tod_type == 'none') cycle
            if (.not. data(i)%tod%subtract_zodi) cycle
            j = j + 1
            step_sizes_albedo(j, :) = 0.1 * data(i)%tod%zodi_albedo
            step_sizes_emissivity(j, :) = 0.1 * data(i)%tod%zodi_emissivity
        end do
                
        allocate(priors(zodi_model%n_parameters, 2))

        emissivity_prior = [0., 5.]
        albedo_prior = [0., 1.]

    end subroutine initialize_zodi_samp_mod

    subroutine initialize_zodi_samp_mod_old(cpar)
        ! Initialize the zodi sampling module.
        !
        ! Parameters
        ! ----------
        ! cpar: comm_params
        !    Parameter file variables.

        type(comm_params), intent(in) :: cpar
        integer(i4b) :: i, j, n_spec_params
        real(dp), allocatable :: param_vec(:)

        chisq_previous = 1d20
        step_size = 0.01

        allocate(priors(zodi_model%n_parameters, 2))
        allocate(active_params(zodi_model%n_parameters))
        allocate(param_vec(zodi_model%n_parameters))
        param_vec = zodi_model%model_to_param_vec()

        allocate(step_sizes(size(param_vec)))
        step_sizes = 0.1 * param_vec

        n_spec_params = zodi_model%n_bands * zodi_model%n_comps


        ! Set active parameters (Those we want to fit)
        ! Cloud
        active_params(1) = .false. !n_0
        active_params(2) = .true. !incl
        active_params(3) = .true. !omega
        active_params(4) = .false. !x_0
        active_params(5) = .false. !y_0
        active_params(6) = .false. !z_0
        active_params(7) = .false. !alpha
        active_params(8) = .false. !beta
        active_params(9) = .false. !gamma
        active_params(10) = .false. !mu
        cloud_start = 1 ; cloud_stop = 10

        ! Band1
        active_params(11) = .false. !n_0
        active_params(12) = .false. !incl
        active_params(13) = .false. !omega
        active_params(14) = .false. !x_0
        active_params(15) = .false. !y_0
        active_params(16) = .false. !z_0
        active_params(17) = .false. !delta_zeta
        active_params(18) = .false. !delta_r
        active_params(19) = .false. !v
        active_params(20) = .false. !p
        band1_start = 11 ; band1_stop = 20

        ! Band2
        active_params(21) = .false. !n_0
        active_params(22) = .false. !incl
        active_params(23) = .false. !omega
        active_params(24) = .false. !x_0
        active_params(25) = .false. !y_0
        active_params(26) = .false. !z_0
        active_params(27) = .false. !delta_zeta
        active_params(28) = .false. !delta_r
        active_params(29) = .false. !v
        active_params(30) = .false. !p
        band2_start = 21 ; band2_stop = 30

        ! Band3
        active_params(31) = .false. !n_0
        active_params(32) = .false. !incl
        active_params(33) = .false. !omega
        active_params(34) = .false. !x_0
        active_params(35) = .false. !y_0
        active_params(36) = .false. !z_0
        active_params(37) = .false. !delta_zeta
        active_params(38) = .false. !delta_r
        active_params(39) = .false. !v
        active_params(40) = .false. !p
        band3_start = 31 ; band3_stop = 40

        ! Ring
        active_params(41) = .false. !n_0
        active_params(42) = .false. !incl
        active_params(43) = .false. !omega
        active_params(44) = .false. !x_0
        active_params(45) = .false. !y_0
        active_params(46) = .false. !z_0
        active_params(47) = .false. !R_0
        active_params(48) = .false. !sigma_r
        active_params(49) = .false. !sigma_z
        ring_start = 41 ; ring_stop = 49

        ! Feature
        active_params(50) = .false. !n_0
        active_params(51) = .false. !incl
        active_params(52) = .false. !omega
        active_params(53) = .false. !x_0
        active_params(54) = .false. !y_0
        active_params(55) = .false. !z_0
        active_params(56) = .false. !R_0
        active_params(57) = .false. !sigma_r
        active_params(58) = .false. !sigma_z
        active_params(59) = .false. !theta_0
        active_params(60) = .false. !sigma_theta
        feature_start = 50 ; feature_stop = 60

        ! Other
        active_params(61) = .false. !T_0
        active_params(62) = .false. !delta
        
        n_active_geom_params = count(active_params == .true.)

        ! Deactivate emissivities for band 6 (reference emissivity of 1.0)
        do i = 1, zodi_model%n_comps
            active_params(get_emissivity_idx(i, 6, zodi_model%n_bands, zodi_model%n_comps)) = .false.
        end do

        n_active_params = count(active_params == .true.)

        ! Priors for zodi parameters
        ! Cloud
        priors(1, :) = [1.13d-8, 1.13d-6] !n_0
        priors(2, :) = NO_PRIOR !incl
        priors(3, :) = NO_PRIOR !omega
        priors(4, :) = [-0.25, 0.25] !x_0
        priors(5, :) = [-0.25, 0.25] !y_0
        priors(6, :) = [-0.25, 0.25] !z_0
        priors(7, :) = NO_PRIOR !alpha
        priors(8, :) = NO_PRIOR !beta
        priors(9, :) = NO_PRIOR !gamma
        priors(10, :) = NO_PRIOR !mu

        ! Band1
        priors(11, :) = [5.59d-12, 5.59d-8] !n_0
        priors(12, :) = NO_PRIOR !incl
        priors(13, :) = NO_PRIOR !omega
        priors(14, :) = [-0.25, 0.25] !x_0
        priors(15, :) = [-0.25, 0.25] !y_0
        priors(16, :) = [-0.25, 0.25] !z_0
        priors(17, :) = NO_PRIOR !delta_zeta
        priors(18, :) = [0.5, 6.] !delta_r
        priors(19, :) = NO_PRIOR !v
        priors(20, :) = NO_PRIOR !p

        ! Band2
        priors(21, :) = [1.99d-11, 1.99d-7] !n_0
        priors(22, :) = NO_PRIOR !incl
        priors(23, :) = NO_PRIOR !omega
        priors(24, :) = [-0.25, 0.25] !x_0
        priors(25, :) = [-0.25, 0.25] !y_0
        priors(26, :) = [-0.25, 0.25] !z_0
        priors(27, :) = NO_PRIOR !delta_zeta
        priors(28, :) = [0.5, 6.] !delta_r
        priors(29, :) = NO_PRIOR !v
        priors(30, :) = NO_PRIOR !p

        ! Band3
        priors(31, :) = [1.44d-12, 1.44d-8] !n_0
        priors(32, :) = NO_PRIOR !incl
        priors(33, :) = NO_PRIOR !omega
        priors(34, :) = [-0.25, 0.25] !x_0
        priors(35, :) = [-0.25, 0.25] !y_0
        priors(36, :) = [-0.25, 0.25] !z_0
        priors(37, :) = NO_PRIOR !delta_zeta
        priors(38, :) = [0., 6.] !delta_r
        priors(39, :) = NO_PRIOR !v
        priors(40, :) = NO_PRIOR !p

        ! Ring
        priors(41, :) = [1.83d-10, 1.83d-6] !n_0
        priors(42, :) = NO_PRIOR !incl
        priors(43, :) = NO_PRIOR !omega
        priors(44, :) = [-0.25, 0.25] !x_0
        priors(45, :) = [-0.25, 0.25] !y_0
        priors(46, :) = [-0.25, 0.25] !z_0
        priors(47, :) = [0.5, 2.] !R_0
        priors(48, :) = NO_PRIOR !sigma_r
        priors(49, :) = NO_PRIOR !sigma_z

        ! Feature
        priors(50, :) = [1.9d-10, 1.9d-6] !n_0
        priors(51, :) = NO_PRIOR !incl
        priors(52, :) = NO_PRIOR !omega
        priors(53, :) = [-0.25, 0.25] !x_0
        priors(54, :) = [-0.25, 0.25] !y_0
        priors(55, :) = [-0.25, 0.25] !z_0
        priors(56, :) = [0.5, 2.] !R_0
        priors(57, :) = NO_PRIOR !sigma_r
        priors(58, :) = NO_PRIOR !sigma_z
        priors(59, :) = NO_PRIOR !theta_0
        priors(60, :) = NO_PRIOR !sigma_theta

        ! Other
        priors(61, :) = [200., 400.] !T_0
        priors(62, :) = [0.1, 1.] !delta


        do i = 1, numband
            if (data(i)%tod_type == 'none') cycle
            if (.not. data(i)%tod%subtract_zodi) cycle
            
        end do
        ! Emissivity
        priors(emissivity_start:emissivity_stop, 1) = 0.
        priors(emissivity_start:emissivity_stop, 2) = 5.

        ! Albedo
        priors(albedo_start:albedo_stop, 1) = 0.
        priors(albedo_start:albedo_stop, 2) = 1.


    end subroutine initialize_zodi_samp_mod_old

    function get_boxwidth(samprate_lowres, samprate) result(box_width)
        ! Returns the boxcar width for downsampling the zodi tods
        real(dp), intent(in) :: samprate_lowres, samprate
        real(dp) :: box_width
        
        if (samprate_lowres >= samprate) stop "Cannot downsample zodi tods if lowres samprate is greater than highres samprate!"
        box_width = samprate / samprate_lowres
        box_width = real(nint(box_width))
        if (box_width < 1.) stop "Cannot downsample zodi tods if box car width is less than 1 sample!"
        if (mod(int(box_width), 2) == 0) box_width = box_width + 1.
    end function

    function get_emissivity_idx(comp_idx, band, n_bands, n_comps) result(idx)
        ! Returns the index of the emissivity for a given component and band
        integer(i4b), intent(in) :: comp_idx, band, n_bands, n_comps
        integer(i4b) :: idx
        idx = 62 + (band - 1) * n_comps + comp_idx
    end function get_emissivity_idx

    function get_albedo_idx(comp_idx, band, n_bands, n_comps) result(idx)
        ! Returns the index of the albedo for a given component and band
        integer(i4b), intent(in) :: comp_idx, band, n_bands, n_comps
        integer(i4b) :: idx
        idx = 62 + (n_comps * n_bands) + (band - 1) * n_comps + comp_idx
    end function get_albedo_idx

    function get_comp_from_param_idx(param_idx) result(comp)
        ! Returns the component index for a given parameter index
        integer(i4b), intent(in) :: param_idx

        integer(i4b) :: comp

        if (cloud_start <= param_idx <= cloud_stop) then 
            comp = 1
        else if (band1_start <= param_idx <= band1_stop) then
            comp = 2
        else if (band2_start <= param_idx <= band2_stop) then
            comp = 3
        else if (band3_start <= param_idx <= band3_stop) then
            comp = 4
        else if (ring_start <= param_idx <= ring_stop) then
            comp = 5
        else if (feature_start <= param_idx <= feature_stop) then
            comp = 6
        else
            comp = 0
        end if

    end function get_comp_from_param_idx

    function get_powell_vec(param_vec) result(powell_vec)
        real(dp), allocatable :: param_vec(:)
        real(dp), allocatable :: powell_vec(:)
        powell_vec = pack(param_vec, active_params == .true.)
    end function get_powell_vec

    subroutine update_param_vec_from_powell_vec(param_vec, powell_vec)
        real(dp), intent(inout) :: param_vec(:)
        real(dp), intent(in) :: powell_vec(:)
        integer :: i, j
        j = 1
        do i = 1, size(param_vec)
            if (active_params(i) == .true.) then
                param_vec(i) = powell_vec(j)
                j = j + 1
            end if
        end do
    end subroutine update_param_vec_from_powell_vec


    ! subroutine sample_zodi_model(cpar, handle)
    !     ! Metropolis-Hastings for nproposals of new sets of zodi parameters
    !     ! Todo: Split into gibbs steps for each component
    !     type(comm_params), intent(in) :: cpar
    !     type(planck_rng), intent(inout) :: handle
    !     integer(i4b) :: i, j, k, flag, ndet, nscan, ntod, nprop, scan, ierr, n_accepted, n_tot_tod, n_tot_tod_reduced, n_tot_tod_current
    !     real(sp), allocatable :: s_scat(:, :), s_therm(:, :), s_zodi(:)
    !     real(dp) :: chisq_tod, chisq_current, chisq_diff, ln_acceptance_probability, accept_rate, chisq_lnL
    !     real(dp), allocatable :: param_vec(:), powell_vec(:)
    !     type(ZodiModel) :: current_model, previous_model
    !     logical(lgt) :: accepted, operation
    !     ! real(dp) :: n0_grid(100), chisq_grid(100)

    !     !debugging stuff
    !     character(len=256) :: chaindir
    !     character(len=512) :: path, param_path, param_name
    !     character(len=6) :: itext
    !     type(hdf_file) :: tod_file
    !     integer(i4b) :: l

    !     ! Metropolis-Hastings for nproposals of new sets of zodi parameters
    !     nprop = cpar%zs_nprop
    !     n_accepted = 0
    !     current_model = zodi_model
    !     previous_model = zodi_model

    !     allocate(param_vec(current_model%n_parameters))
    !     allocate(powell_vec(current_model%n_comps))  

    !     ! param_vec = current_model%model_to_param_vec()
    !     ! powell_vec = get_powell_vec(param_vec)
        
    !     ! Start initial parameter minimization. Root will aquire new set of parameters. While this happens
    !     ! the other processes will wait for the new parameters to be broadcasted. Then they will all compute
    !     ! the their own lnL.
    !     ! if (cpar%myid == cpar%root) then
    !     !     print *, "before", powell_vec
    !     !     call powell(powell_vec, lnL_zodi, ierr)
    !     !     flag = 0
    !     !     call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
    !     ! else
    !     !     do while (.true.)
    !     !         call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
    !     !         if (flag == 1) then 
    !     !             ! Here we call the function without a powell vector, because it will recieve it in internally from MPI
    !     !             chisq_lnL = lnL_zodi()
    !     !             ! print *, chisq_lnL
    !     !         else
    !     !             exit
    !     !         end if
    !     !     end do
    !     ! end if

    !     ! ! Wait for every process to exit
    !     ! call mpi_barrier(MPI_COMM_WORLD, ierr)

    !     ! call mpi_bcast(powell_vec, size(powell_vec), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
    !     ! if (cpar%myid == cpar%root) then
    !     !     print *, "after", powell_vec, "ierr", ierr
    !     ! end if
    !     ! stop

    !     ! ! update param_vec with new parameters from powell
    !     ! call update_param_vec_from_powell_vec(param_vec, powell_vec)

    !     ! update model with new param_vec
    !     ! call current_model%param_vec_to_model(param_vec)

    !     do k = 0, nprop
    !         ! Reset chisq for current proposal
    !         chisq_tod = 0.
    !         chisq_current = 0.
            
    !         ! Root process draws new set of zodi parameters and broadcasts to all processes
    !         ! If first iteration we dont want to draw, just compute the chisq with the base model
    !         if (k > 0) then
    !             if (cpar%myid == cpar%root) then
    !                 param_vec = current_model%model_to_param_vec()
    !                 do i = 1, size(param_vec)
    !                     if (active_params(i) == .true.) then
    !                         param_vec(i) = param_vec(i) + (step_sizes(i) * param_vec(i) * rand_gauss(handle))
    !                         if (priors(i, 1) == NO_PRIOR .or. priors(i, 2) == NO_PRIOR) then 
    !                             continue
    !                         else if (param_vec(i) < priors(i, 1) .or. param_vec(i) > priors(i, 2)) then 
    !                             step_sizes(i) = step_sizes(i) / 2.
    !                             chisq_tod = 1.d30
    !                             exit
    !                         end if
    !                     end if
    !                 end do
    !             end if
    !             call mpi_bcast(param_vec, size(param_vec), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
    !             call current_model%param_vec_to_model(param_vec)
    !         end if

    !         do i = 1, numband
    !             ! Skip none tod bands
    !             if (data(i)%tod_type == "none") cycle
    !             ! Skip tod bands where we dont want to sample zodi
    !             if (.not. data(i)%tod%sample_zodi) cycle

    !             if (chisq_tod >= 1.d30) exit

    !             ndet = data(i)%tod%ndet
    !             nscan = data(i)%tod%nscan
                
    !             ! Make sure that the zodi cache is cleared before each new band
    !             call data(i)%tod%clear_zodi_cache()

    !             ! Evaluate zodi model with newly proposed values for each band and calculate chisq
    !             do scan = 1, nscan
    !                 ! Skip scan if no accepted data
    !                 do j = 1, ndet
    !                     if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
    !                     ntod = size(data(i)%tod%scans(scan)%d(j)%downsamp_res)

    !                     allocate(s_scat(ntod, zodi_model%n_comps), s_therm(ntod, zodi_model%n_comps), s_zodi(ntod))
    !                     call get_zodi_emission(&
    !                         & tod=data(i)%tod, &
    !                         & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pix, &
    !                         & scan=scan, &
    !                         & det=j, &
    !                         & s_zodi_scat=s_scat, &
    !                         & s_zodi_therm=s_therm, &
    !                         & model=current_model &
    !                     &)
    !                     call get_s_zodi(&
    !                         & s_therm=s_therm, &
    !                         & s_scat=s_scat, &
    !                         & s_zodi=s_zodi, &
    !                         & emissivity=data(i)%tod%zodi_emissivity, &
    !                         & albedo=data(i)%tod%zodi_albedo &
    !                     &)
                    
    !                     chisq_tod = chisq_tod + sum(((data(i)%tod%scans(scan)%d(j)%downsamp_res - s_zodi)/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2)
    !                     deallocate(s_scat, s_therm, s_zodi)
    !                 end do
    !             end do
    !         end do

    !         ! Reduce chisq to root process
    !         call mpi_reduce(chisq_tod, chisq_current, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cpar%root, MPI_COMM_WORLD, ierr)

    !         ! Use chisq from the first iteration where we dont draw new parameters as the base chisq
    !         if (k == 0) chisq_previous = chisq_current

    !         if (k > 0) then
    !             ! Root checks if the new proposal is accepted
    !             if (cpar%myid == cpar%root) then
    !                 chisq_diff = max(chisq_current - chisq_previous, 0.)
    !                 ln_acceptance_probability = -0.5 * chisq_diff
    !                 accepted = ln_acceptance_probability > log(rand_uni(handle))
                    
    !                 print *, "proposal:", k
    !                 print *, "chisq_current:", chisq_current
    !                 print *, "chisq_diff:", chisq_diff
    !                 print *, "ln a:", ln_acceptance_probability 
    !                 print *, "accepted:", accepted
    !                 print *, "accept rate:", (real(n_accepted) / real(k))*100.
    !                 print *, " "
    !             end if

    !             call mpi_bcast(accepted, 1, MPI_LOGICAL, cpar%root, cpar%comm_chain, ierr)
    !             ! Update model if proposal is accepted
    !             if (accepted) then
    !                 n_accepted = n_accepted + 1
    !                 chisq_previous = chisq_current
    !                 previous_model = current_model
    !             else
    !                 current_model = previous_model
    !             end if
    !             accept_rate = n_accepted / k
    !         end if
    !     end do  
    !     if (accept_rate < 10.) then
    !         step_sizes = step_sizes / 2.
    !         print *, "lowering stepsize"
    !     else if (accept_rate > 60.) then 
    !         step_sizes = step_sizes * 2.
    !         print *, "increasing stepsize"
    !     end if
        
    !     zodi_model = current_model
    ! end subroutine


    ! subroutine init_scandata_and_downsamp_zodi(cpar, comp_idx)
    !     !     ! This routine mimics parts of the process_TOD where we calibrate the timestreams 
    !     !     ! put them into th scan_data type, and then downsample them for use in zodi fitting.
    !     !     ! See the process_TOD routines for more details on what happens here.
    !     !     type(comm_params), intent(in) :: cpar
    !     !     integer(i4b), intent(in), optional :: comp_idx

    !     !     integer(i4b) :: i, j , k, nside, npix, nmaps, ndelta, box_width
    !     !     real(sp), allocatable, dimension(:, :, :, :) :: map_sky
    !     !     real(sp), allocatable, dimension(:) :: procmask, procmask2, procmask_zodi
    !     !     real(dp), allocatable, dimension(:, :) :: m_buf
    !     !     integer(i4b), allocatable, dimension(:) :: downsamp_width
    !     !     logical(lgt) :: downsamp_all_equally
    !     !     type(comm_scandata) :: sd
    !     !     type(map_ptr), allocatable, dimension(:, :) :: s_sky


    !     !     if (present(comp_idx)) then
    !     !         allocate(downsamp_width(zodi_model%n_comps))
    !     !         downsamp_all_equally = .false.
    !     !     else 
    !     !         downsamp_all_equally = .true.
    !     !     end if

    !     !     ! For each zodi band, create downsampled residual time stream
    !     !     do i = 1, numband

    !     !         ! Only generate downsampled arrays for tod bands and bands where we have zodi
    !     !         if (trim(data(i)%tod_type) == 'none') cycle
    !     !         if (.not. data(i)%tod%subtract_zodi) cycle

    !     !         ! Clear previously allocated arrays if they exist
    !     !         call data(i)%tod%deallocate_downsampled_zodi()
    !     !         call data(i)%tod%clear_zodi_cache()

    !     !         ndelta = 1
    !     !         nside = data(i)%map%info%nside
    !     !         npix = 12 * nside**2
    !     !         nmaps = data(i)%map%info%nmaps

    !     !         allocate(m_buf(0:npix-1, nmaps), procmask(0:npix-1), procmask2(0:npix-1), procmask_zodi(0:npix-1))
    !     !         call data(i)%tod%procmask%bcast_fullsky_map(m_buf); procmask  = m_buf(:,1)
    !     !         call data(i)%tod%procmask2%bcast_fullsky_map(m_buf); procmask2 = m_buf(:,1)
    !     !         call data(i)%tod%procmask_zodi%bcast_fullsky_map(m_buf); procmask_zodi = m_buf(:,1)
                
    !     !         allocate(s_sky(data(i)%tod%ndet, ndelta))
    !     !         do j = 1, data(i)%tod%ndet
    !     !             call get_sky_signal(i, j, s_sky(j, ndelta)%p, mono=.false.)
    !     !         end do

    !     !         allocate(map_sky(nmaps, data(i)%tod%nobs, 0:data(i)%tod%ndet, 1))
    !     !         call distribute_sky_maps(data(i)%tod, s_sky, 1.e0, map_sky)

    !     !         ! Start with same downsamp box width for all comps (change in future)
    !     !         if (data(i)%tod%samprate_lowres >= data(i)%tod%samprate) stop "Cannot downsample zodi tods if lowres samprate is greater than highres samprate!"
    !     !         box_width = data(i)%tod%samprate / data(i)%tod%samprate_lowres
    !     !         if (box_width < 1.) stop "Cannot downsample zodi tods if box car width is less than 1 sample!"

    !     !     ! Make box have odd size so that we can pick out a center pixel
    !     !         if (mod(box_width, 2) == 0) then 
    !     !             box_fullwidth = box_width + 1.
    !     !         else 
    !     !             box_fullwidth = box_width
    !     !         end if
    !     !         box_halfwidth = box_fullwidth / 2
    !     !         do j = 1, data(i)%tod%nscan
    !     !             if (.not. any(data(i)%tod%scans(j)%d%accept)) cycle

    !     !             ! not calibrating the data so we use substitute map_gain with map_sky
    !     !             call sd%init_singlehorn(data(i)%tod, j, map_sky, map_sky, procmask, procmask2, procmask_zodi)
    !     !             call boxcar_downsamp_zodi_res_and_pointing(data(i)%tod, sd, j, padding=5, box_width=box_width)
    !     !             ! prepare scandata object for next scan
    !     !             call sd%dealloc
    !     !         end do

    !     !         deallocate(m_buf, procmask, procmask2, procmask_zodi)
    !     !         deallocate(s_sky)
    !     !         deallocate(map_sky)
    !     !     end do
    ! end subroutine


    ! subroutine boxcar_downsamppointing(tod, sd, scan_id, padding, box_width)
    !     ! Downsamples the tod residual and pointing given a box car width for use in zodi sampling.
    !     class(comm_tod), intent(inout) :: tod
    !     !     type(comm_scandata), intent(in) :: sd
    !     !     integer(i4b), intent(in) :: scan_id, padding
    !     !     integer(i4b), intent(in) :: box_width

    !     !     integer(i4b) :: i, j, ext(2), upper_bound, box_halfwidth
    !     !     real(dp) :: box_fullwidth
    !     !     real(sp), allocatable, dimension(:) :: downsamp_mask, downsamp_pointing
    !     !     logical(lgt), allocatable, dimension(:) :: downsamp_masked_indices

    !     !     if (box_width < 1.) stop "Cannot downsample zodi tods if box car width is less than 1 sample!"

    !     !    ! Make box have odd size so that we can pick out a center pixel
    !     !     if (mod(box_width, 2) == 0) then 
    !     !         box_fullwidth = box_width + 1.
    !     !     else 
    !     !         box_fullwidth = box_width
    !     !     end if
    !     !     box_halfwidth = box_fullwidth / 2

    !     !     do j = 1, tod%ndet
    !     !         ! Get downsampled shape (ext), and allocate the downsampled arrays
    !     !         call tod%downsample_tod(sd%tod(:, j), ext, step=box_fullwidth)

    !     !         ! Allocate these the first gibbs iter
    !     !         allocate(downsamp_pointing(ext(1):ext(2)))
    !     !         allocate(downsamp_mask(ext(1):ext(2)))
    !     !         allocate(downsamp_masked_indices(ext(1):ext(2)))

    !     !         downsamp_mask = 0.
    !     !         downsamp_pointing = 0

    !     !         ! Downsample the mask
    !     !         call tod%downsample_tod(sd%mask_zodi(:, j), ext, downsamp_mask, step=box_fullwidth)

    !     !         ! Construct the downsampled pointing array (pick central pixel in each box)
    !     !         do i = 1, int(size(downsamp_mask) / box_fullwidth)
    !     !             downsamp_pointing(i) = sd%pix((i * box_fullwidth), j, 1)
    !     !         end do

    !     !         downsamp_masked_indices = .true.
    !     !         where (downsamp_mask < 1.) downsamp_masked_indices = .false.
                
    !     !         ! Find upper bound and truncate the downsampled arrays by where they were padded and
    !     !         ! filter out the masked values. Store the result in the tod objects for use in sampling.
    !     !         upper_bound = ubound(downsamp_mask, dim=1)
    !     !         tod%scans(scan_id)%d(j)%downsamp_pointing = pack(downsamp_pointing(padding:upper_bound-padding), downsamp_masked_indices(padding:upper_bound-padding))
    !     !         deallocate(downsamp_pointing, downsamp_mask, downsamp_masked_indices)
    !     !     end do
    ! end subroutine boxcar_downsamp_zodi_res_and_pointing


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


    ! subroutine gibbs_sample_zodi_comp(cpar, handle, comp_idx, zodi_model_samp, verbose)
    !     !     type(comm_params), intent(in) :: cpar
    !     !     type(planck_rng), intent(inout) :: handle
    !     !     integer(i4b), intent(in) :: comp_idx
    !     !     type(ZodiModel), intent(inout) :: zodi_model_samp
    !     !     logical(lgt), intent(in), optional :: verbose

    !     !     integer(i4b) :: i, j, k, param_start, param_stop, flag, ndet, nscan, ntod, nprop, scan, ierr, n_accepted, n_tot_tod, n_tot_tod_reduced, n_tot_tod_current
    !     !     real(sp), allocatable :: s_scat(:, :), s_therm(:, :), s_zodi(:)
    !     !     real(dp) :: chisq_tod, chisq_current, chisq_diff, ln_acceptance_probability, accept_rate, chisq_lnL
    !     !     real(dp), allocatable :: param_vec(:), powell_vec(:)
    !     !     type(ZodiModel) :: current_model, previous_model
    !     !     logical(lgt) :: accepted, verbose_
    !     !     type(hdf_file) :: tod_file

    !     !     if (present(verbose)) then
    !     !         verbose_ = verbose
    !     !     else
    !     !         verbose_ = .false.
    !     !     end if

    !     !     ! Get parameter indicies for this component
    !     !     select case (trim(zodi_model%comp_labels(comp_idx))) 
    !     !     case ("cloud")
    !     !         param_start = cloud_start
    !     !         param_stop = cloud_stop
    !     !     case ("band1")
    !     !         param_start = band1_start
    !     !         param_stop = band1_stop
    !     !     case ("band2")
    !     !         param_start = band2_start
    !     !         param_stop = band2_stop
    !     !     case ("band3")
    !     !         param_start = band3_start
    !     !         param_stop = band3_stop
    !     !     case ("ring")
    !     !         param_start = ring_start
    !     !         param_stop = ring_stop
    !     !     case ("feature")
    !     !         param_start = feature_start
    !     !         param_stop = feature_stop
    !     !     case default 
    !     !         stop "Error: invalid component label"
    !     !     end select

    !     !     if (verbose_ .and. (cpar%myid == cpar%root)) then
    !     !         print *, "==== Zodi gibbs sampling step: ", trim(zodi_model%comp_labels(comp_idx)), " ===="
    !     !         print *, " "
    !     !     end if 

    !     !     ! Metropolis-Hastings for nproposals of new sets of zodi parameters
    !     !     nprop = cpar%zs_nprop

    !     !     ! Make two copies of the zodi model, one for the current proposal and one for the previous proposal
    !     !     current_model = zodi_model_samp
    !     !     previous_model = zodi_model_samp

    !     !     allocate(param_vec(current_model%n_parameters))

    !     !     n_accepted = 0
    !     !     do k = 0, nprop
    !     !         ! Reset chisq for current proposal
    !     !         chisq_tod = 0.
    !     !         chisq_current = 0.
                
    !     !         ! If first iteration we dont want to draw, just compute the chisq with the base model
    !     !         if (k > 0) then

    !     !             ! Root process draws new set of zodi parameters and broadcasts to all processes
    !     !             if (cpar%myid == cpar%root) then
    !     !                 param_vec = current_model%model_to_param_vec()

    !     !                 ! propse new set of parameters for component
    !     !                 do i = param_start, param_stop
    !     !                     if (.not. active_params(i)) cycle
    !     !                     param_vec(i) = param_vec(i) + (step_sizes(i) * rand_gauss(handle))

    !     !                     ! If parameter is outside of prior range, immediatly reject proposal and halve step size
    !     !                     if (priors(i, 1) == NO_PRIOR .or. priors(i, 2) == NO_PRIOR) then 
    !     !                         continue
    !     !                     else if (param_vec(i) < priors(i, 1) .or. param_vec(i) > priors(i, 2)) then 
    !     !                         step_sizes(i) = step_sizes(i) / 2.
    !     !                         chisq_tod = 1.d30
    !     !                         exit
    !     !                     end if
    !     !                 end do
    !     !             end if
    !     !             call mpi_bcast(param_vec, size(param_vec), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
    !     !             ! Update current_model with new parameters
    !     !             call current_model%param_vec_to_model(param_vec)
    !     !         end if

    !     !         ! loop over downsampled residuals for each tod band with zodi
    !     !         do i = 1, numband
    !     !             if (data(i)%tod_type == "none") cycle
    !     !             if (.not. data(i)%tod%sample_zodi) cycle

    !     !             ! If chisq is already too large, skip rest of the evaluation and go directly to rejection
    !     !             if (chisq_tod >= 1.d30) exit

    !     !             ndet = data(i)%tod%ndet
    !     !             nscan = data(i)%tod%nscan
                    
    !     !             ! Make sure that the zodi cache is cleared before each new band
    !     !             call data(i)%tod%clear_zodi_cache()

    !     !             ! Evaluate zodi model with newly proposed values for each band and calculate chisq
    !     !             do scan = 1, nscan
    !     !                 ! Skip scan if no accepted data
    !     !                 do j = 1, ndet
    !     !                     if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
    !     !                     ntod = size(data(i)%tod%scans(scan)%d(j)%downsamp_res)

    !     !                     allocate(s_scat(ntod, zodi_model%n_comps), s_therm(ntod, zodi_model%n_comps), s_zodi(ntod))
    !     !                     call get_zodi_emission(&
    !     !                         & tod=data(i)%tod, &
    !     !                         & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pointing, &
    !     !                         & scan=scan, &
    !     !                         & det=j, &
    !     !                         & s_zodi_scat=s_scat, &
    !     !                         & s_zodi_therm=s_therm, &
    !     !                         & model=current_model &
    !     !                     &)
    !     !                     call get_s_zodi(&
    !     !                         & s_therm=s_therm, &
    !     !                         & s_scat=s_scat, &
    !     !                         & s_zodi=s_zodi, &
    !     !                         & emissivity=data(i)%tod%zodi_emissivity, &
    !     !                        & albedo=data(i)%tod%zodi_albedo &
    !     !                     &)
    !     !                     ! if (.true. .and. cpar%myid == cpar%root) then
    !     !                     !     call open_hdf_file(trim("/mn/stornext/d5/data/metins/dirbe/chains/chains_testing/downsamp_timestream.h5"), tod_file, 'w')
    !     !                     !     call write_hdf(tod_file, '/r', data(i)%tod%scans(scan)%d(j)%downsamp_res)
    !     !                     !     call write_hdf(tod_file, '/p', data(i)%tod%scans(scan)%d(j)%downsamp_pointing)
    !     !                     !     call write_hdf(tod_file, '/z', s_zodi)
    !     !                     !     call close_hdf_file(tod_file)
    !     !                     !     stop
    !     !                     ! end if
    !     !                     chisq_tod = chisq_tod + sum(((data(i)%tod%scans(scan)%d(j)%downsamp_res - s_zodi)/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2)
    !     !                     deallocate(s_scat, s_therm, s_zodi)
    !     !                 end do
    !     !             end do
    !     !         end do

    !     !         ! Reduce chisq to root process
    !     !         call mpi_reduce(chisq_tod, chisq_current, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cpar%root, MPI_COMM_WORLD, ierr)

    !     !         ! Use chisq from the first iteration where we dont draw new parameters as the base chisq
    !     !         if (k == 0) chisq_previous = chisq_current

    !     !         if (k > 0) then
    !     !             ! Root checks if the new proposal is accepted
    !     !             if (cpar%myid == cpar%root) then
    !     !                 chisq_diff = max(chisq_current - chisq_previous, 0.)
    !     !                 ln_acceptance_probability = -0.5 * chisq_diff
    !     !                 accepted = ln_acceptance_probability > log(rand_uni(handle))
    !     !                 if (verbose_) then
    !     !                     print *, "proposal:", k
    !     !                     print *, "chisq_current:", chisq_current
    !     !                     print *, "chisq_diff:", chisq_diff
    !     !                     print *, "ln a:", ln_acceptance_probability 
    !     !                     print *, "accepted:", accepted
    !     !                     print *, "accept rate:", (real(n_accepted) / real(k))*100.
    !     !                     print *, " "
    !     !                 end if
    !     !             end if

    !     !             call mpi_bcast(accepted, 1, MPI_LOGICAL, cpar%root, cpar%comm_chain, ierr)
    !     !             ! Update model if proposal is accepted
    !     !             if (accepted) then
    !     !                 n_accepted = n_accepted + 1
    !     !                 chisq_previous = chisq_current
    !     !                 previous_model = current_model
    !     !             else
    !     !                 current_model = previous_model
    !     !             end if
    !     !             accept_rate = real(n_accepted) / real(k)
    !     !         end if
    !     !     end do  
    !     !     if (accept_rate < 0.1) then
    !     !         step_sizes(param_start:param_stop) = step_sizes(param_start:param_stop) / 2.
    !     !         if (verbose_ .and. (cpar%myid == cpar%root)) print *, "lowering stepsize", accept_rate
    !     !     else if (accept_rate > 0.6) then 
    !     !         step_sizes(param_start:param_stop) = step_sizes(param_start:param_stop) * 2.
    !     !         if (verbose_ .and. (cpar%myid == cpar%root)) print *, "increasing stepsize", accept_rate
    !     !     end if
            
    !     !     zodi_model_samp = current_model
    ! end subroutine


    subroutine sample_zodi_emissivity_and_albedo(cpar, handle, gibbs_iter, zodi_model_samp, verbose)
        type(comm_params), intent(in) :: cpar
        type(planck_rng), intent(inout) :: handle
        integer(i4b), intent(in) :: gibbs_iter
        type(ZodiModel), intent(inout) :: zodi_model_samp
        logical(lgt), intent(in), optional :: verbose

        integer(i4b) :: i, j, k, ndet, nscan, ntod, nprop, scan, ierr, n_accepted
        integer(i4b) :: albedo_idx, emissivity_idx
        real(sp), allocatable :: s_zodi_scan(:)
        real(dp) :: chisq_tod, chisq_current, chisq_diff, ln_acceptance_probability, accept_rate, chisq_lnL
        real(dp), allocatable :: param_vec(:), powell_vec(:), emissivity_new(:, :), albedo_new(:, :), emissivity_prev(:, :), albedo_prev(:, :)
        logical(lgt) :: accepted, verbose_, tuning, skip
        type(hdf_file) :: tod_file

        if (present(verbose)) then
            verbose_ = verbose
        else
            verbose_ = .false.
        end if


        if (verbose_ .and. (cpar%myid == cpar%root)) then
            print *, "==== Zodi emissivity and albedo sampling step ===="
            print *, " "
        end if 

        allocate(emissivity_new(n_samp_bands, zodi_model%n_comps))
        allocate(albedo_new(n_samp_bands, zodi_model%n_comps))
        allocate(emissivity_prev(n_samp_bands, zodi_model%n_comps))
        allocate(albedo_prev(n_samp_bands, zodi_model%n_comps))

        ! Metropolis-Hastings for nproposals of new sets of zodi parameters
        nprop = 200
        ! nprop = cpar%zs_nprop

        ! Make two copies of the zodi model, one for the current proposal and one for the previous proposal
        tuning = gibbs_iter <= 25
        do i = 1, numband
            if (data(i)%tod_type == "none") cycle
            if (.not. data(i)%tod%sample_zodi) cycle
            accept_rate = 0.

            ndet = data(i)%tod%ndet
            nscan = data(i)%tod%nscan

            emissivity_new(i, :) = data(i)%tod%zodi_emissivity
            albedo_new(i, :) = data(i)%tod%zodi_albedo
            emissivity_prev = emissivity_new
            albedo_prev = albedo_new

            n_accepted = 0
            do k = 0, nprop
                
                skip = .false.
                chisq_tod = 0.
                chisq_current = 0.

                ! TODO: code in sampling of tod parameters emissivity and albedo
                if (k > 0) then
                    ! Root process draws new set of zodi parameters and broadcasts to all processes
                    if (cpar%myid == cpar%root) then
                        sample_loop: do j = 1, zodi_model%n_comps
                            emissivity_new(i, j) = emissivity_new(i, j) + (step_sizes_emissivity(i, j) * rand_gauss(handle)) 
                            if (emissivity_new(i, j) < emissivity_prior(1) .or. emissivity_new(i, j) > emissivity_prior(2)) then 
                                if (tuning) step_sizes_emissivity(i, j) = step_sizes_emissivity(i, j) / 2.
                                chisq_tod = 1.d30
                                skip = .true.
                                emissivity_new(i, :) = emissivity_prev(i, :)
                                exit sample_loop
                            end if
                            albedo_new(i, j) = albedo_new(i, j) + (step_sizes_albedo(i, j) * rand_gauss(handle)) 
                            if (albedo_new(i, j) < albedo_prior(1) .or. albedo_new(i, j) > albedo_prior(2)) then 
                                if (tuning) step_sizes_albedo(i, j) = step_sizes_albedo(i, j) / 2.
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
                            call get_s_zodi(&
                                & s_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                                & s_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                                & s_zodi=data(i)%tod%scans(scan)%d(j)%downsamp_zodi, &
                                & emissivity=emissivity_new(i, :), &
                                & albedo=albedo_new(i, :) &
                            &)
                            chisq_tod = chisq_tod + sum( &
                               & (  (data(i)%tod%scans(scan)%d(j)%downsamp_tod &
                               &   - data(i)%tod%scans(scan)%d(j)%downsamp_sky &
                               &   - data(i)%tod%scans(scan)%d(j)%downsamp_zodi &
                               & ) / data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2 &
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
                        ln_acceptance_probability = -0.5 * chisq_diff
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
                    accept_rate = real(n_accepted) / real(k)
                end if
            end do  
            if (accept_rate < 0.1) then
                step_sizes_emissivity(i, :) = step_sizes_emissivity(i, :) / 2.
                step_sizes_albedo(i, :) = step_sizes_albedo(i, :) / 2.
            else if (accept_rate > 0.6) then 
                step_sizes_emissivity(i, :) = step_sizes_emissivity(i, :) * 2.
                step_sizes_albedo(i, :) = step_sizes_albedo(i, :) * 2.
            end if
            if (verbose_ .and. (cpar%myid == cpar%root)) print *, "emissivity and albedo accept rate:", (real(n_accepted) / real(k))*100., "band: ", trim(adjustl(data(i)%tod%freq))
        end do
        do i = 1, numband
            data(i)%tod%zodi_emissivity = emissivity_new(i, :)
            data(i)%tod%zodi_albedo = albedo_new(i, :)
        end do
    end subroutine

    subroutine sample_zodi_group(cpar, handle, gibbs_iter, zodi_model_samp, group, verbose)
        type(comm_params), intent(in) :: cpar
        type(planck_rng), intent(inout) :: handle
        integer(i4b), intent(in) :: gibbs_iter
        type(ZodiModel), intent(inout) :: zodi_model_samp
        integer(i4b), intent(in) :: group(:)
        logical(lgt), intent(in), optional :: verbose

        integer(i4b) :: i, j, k, flag, ndet, nscan, ntod, nprop, scan, ierr, n_accepted
        real(dp) :: chisq_tod, chisq_current, chisq_diff, ln_acceptance_probability, accept_rate, chisq_lnL
        real(dp), allocatable :: param_vec(:)
        type(ZodiModel) :: current_model, previous_model
        logical(lgt) :: accepted, verbose_, tuning, skip
        type(hdf_file) :: tod_file

        if (present(verbose)) then
            verbose_ = verbose
        else
            verbose_ = .false.
        end if

        

    end subroutine

    ! subroutine sample_zodi_parameter(cpar, handle, gibbs_iter, param_idx, zodi_model_samp, verbose)
    !     ! Gibbs sample one by one zodi parameter
    !     type(comm_params), intent(in) :: cpar
    !     type(planck_rng), intent(inout) :: handle
    !     integer(i4b), intent(in) :: gibbs_iter
    !     integer(i4b), intent(in) :: param_idx
    !     type(ZodiModel), intent(inout) :: zodi_model_samp
    !     logical(lgt), intent(in), optional :: verbose

    !     integer(i4b) :: i, j, k, param_start, param_stop, flag, ndet, nscan, ntod, nprop, scan, ierr, n_accepted, n_tot_tod, n_tot_tod_reduced, n_tot_tod_current, comp
    !     real(sp), allocatable :: s_scat(:, :), s_therm(:, :), s_zodi(:)
    !     real(dp) :: chisq_tod, chisq_current, chisq_diff, ln_acceptance_probability, accept_rate, chisq_lnL
    !     real(dp), allocatable :: param_vec(:), powell_vec(:)
    !     type(ZodiModel) :: current_model, previous_model
    !     logical(lgt) :: accepted, verbose_, tuning, skip
    !     type(hdf_file) :: tod_file

    !     if (present(verbose)) then
    !         verbose_ = verbose
    !     else
    !         verbose_ = .false.
    !     end if

    !     if (verbose_ .and. (cpar%myid == cpar%root)) then
    !         print *, "==== Zodi gibbs sampling step for param number: ", param_idx, " ===="
    !         print *, " "
    !     end if 

    !     comp = get_comp_from_param_idx(param_idx)
    !     print *, "comp idx: ", comp
    !     stop
    !     ! Metropolis-Hastings for nproposals of new sets of zodi parameters
    !     nprop = cpar%zs_nprop

    !     ! Make two copies of the zodi model, one for the current proposal and one for the previous proposal
    !     current_model = zodi_model_samp
    !     previous_model = zodi_model_samp

    !     allocate(param_vec(current_model%n_parameters))

    !     if (gibbs_iter >= 25) then
    !         tuning = .false.
    !     else
    !         tuning = .true.
    !     end if
    !     n_accepted = 0

    !     do k = 0, nprop
    !         ! if propsed parameter is outside of priors we dont want to immediately reject
    !         skip = .false.

    !         ! Reset chisq for current proposal
    !         chisq_tod = 0.
    !         chisq_current = 0.
    !         ! If first iteration we dont want to draw, just compute the chisq with the base model
    !         if (k > 0) then

    !             ! Root process draws new set of zodi parameters and broadcasts to all processes
    !             param_vec = current_model%model_to_param_vec()
    !             if (cpar%myid == cpar%root) then
    !                 param_vec(param_idx) = param_vec(param_idx) + (step_sizes(param_idx) * rand_gauss(handle))

    !                 ! If parameter is outside of prior range, immediatly reject proposal and halve step size
    !                 if (priors(param_idx, 1) == NO_PRIOR .or. priors(param_idx, 2) == NO_PRIOR) then 
    !                     continue
    !                 else if (param_vec(param_idx) < priors(param_idx, 1) .or. param_vec(param_idx) > priors(param_idx, 2)) then 
    !                     if (tuning) step_sizes(param_idx) = step_sizes(param_idx) / 2.
    !                     chisq_tod = 1.d30
    !                     skip = .true.
    !                 end if
    !             end if
    !             call mpi_bcast(param_vec, size(param_vec), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
    !             ! Update current_model with new parameters
    !             call current_model%param_vec_to_model(param_vec)
    !         end if

    !         ! loop over downsampled residuals for each tod band with zodi
    !         do i = 1, numband
    !             if (skip) exit ! prior was violated, skip rest of evaluation
    !             if (data(i)%tod_type == "none") cycle
    !             if (.not. data(i)%tod%sample_zodi) cycle

    !             ! If chisq is already too large, skip rest of the evaluation and go directly to rejection
    !             if (chisq_tod >= 1.d30) exit

    !             ndet = data(i)%tod%ndet
    !             nscan = data(i)%tod%nscan
                
    !             ! Make sure that the zodi cache is cleared before each new band
    !             call data(i)%tod%clear_zodi_cache()

    !             ! Evaluate zodi model with newly proposed values for each band and calculate chisq
    !             do scan = 1, nscan
    !                 ! Skip scan if no accepted data
    !                 do j = 1, ndet
    !                     if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
    !                     ntod = size(data(i)%tod%scans(scan)%d(j)%downsamp_res)

    !                     allocate(s_scat(ntod, zodi_model%n_comps), s_therm(ntod, zodi_model%n_comps), s_zodi(ntod))
    !                     call get_zodi_emission(&
    !                         & tod=data(i)%tod, &
    !                         & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pix, &
    !                         & scan=scan, &
    !                         & det=j, &
    !                         & s_zodi_scat=s_scat, &
    !                         & s_zodi_therm=s_therm, &
    !                         & model=current_model, &
    !                         & use_lowres_pointing=.true., &
    !                         & comp=comp &
    !                     &)
    !                     call get_s_zodi(&
    !                         & s_therm=s_therm, &
    !                         & s_scat=s_scat, &
    !                         & s_zodi=s_zodi, &
    !                         & emissivity=data(i)%tod%zodi_emissivity, &
    !                         & albedo=data(i)%tod%zodi_albedo &
    !                         & comp=comp &
    !                     &)

    !                     if (.false.) then
    !                         if (cpar%myid == cpar%root) then
    !                             call open_hdf_file(trim("/mn/stornext/d5/data/metins/dirbe/chains/chains_testing/downsamp_timestream.h5"), tod_file, 'w')
    !                             call write_hdf(tod_file, '/r', data(i)%tod%scans(scan)%d(j)%downsamp_res)
    !                             call write_hdf(tod_file, '/z', s_zodi)
    !                             call close_hdf_file(tod_file)
    !                             stop
    !                         end if
    !                         call mpi_barrier(cpar%comm_chain, ierr)
    !                         stop
    !                     end if
                        
    !                     chisq_tod = chisq_tod + sum(((data(i)%tod%scans(scan)%d(j)%downsamp_res - s_zodi)/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2)
    !                     deallocate(s_scat, s_therm, s_zodi)
    !                 end do
    !             end do
    !         end do

    !         ! Reduce chisq to root process
    !         call mpi_reduce(chisq_tod, chisq_current, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cpar%root, MPI_COMM_WORLD, ierr)

    !         ! Use chisq from the first iteration where we dont draw new parameters as the base chisq
    !         if (k == 0) chisq_previous = chisq_current

    !         if (k > 0) then
    !             ! Root checks if the new proposal is accepted
    !             if (cpar%myid == cpar%root) then
    !                 chisq_diff = max(chisq_current - chisq_previous, 0.)
    !                 ln_acceptance_probability = -0.5 * chisq_diff
    !                 select case (cpar%operation)
    !                 case ("sample")
    !                     accepted = ln_acceptance_probability > log(rand_uni(handle))
    !                 case ("optimize")
    !                     accepted = chisq_diff < 0.
    !                 case default
    !                     stop "Error: invalid operation"
    !                 end select

    !                 if (verbose_) then
    !                     print *, "proposal:", k
    !                     print *, "chisq_current:", chisq_current
    !                     print *, "chisq_diff:", chisq_diff
    !                     print *, "ln a:", ln_acceptance_probability 
    !                     print *, "accepted:", accepted
    !                     print *, "accept rate:", (real(n_accepted) / real(k))*100.
    !                     print *, " "
    !                 end if
    !             end if

    !             call mpi_bcast(accepted, 1, MPI_LOGICAL, cpar%root, cpar%comm_chain, ierr)
    !             ! Update model if proposal is accepted
    !             if (accepted) then
    !                 n_accepted = n_accepted + 1
    !                 chisq_previous = chisq_current
    !                 previous_model = current_model
    !             else
    !                 current_model = previous_model
    !             end if
    !             accept_rate = real(n_accepted) / real(k)
    !         end if
    !     end do  
    !     if (accept_rate < 0.1) then
    !         step_sizes(param_idx) = step_sizes(param_idx) / 2.
    !         if (verbose_ .and. (cpar%myid == cpar%root)) print *, "lowering stepsize", accept_rate
    !     else if (accept_rate > 0.6) then 
    !         step_sizes(param_idx) = step_sizes(param_idx) * 2.
    !         if (verbose_ .and. (cpar%myid == cpar%root)) print *, "increasing stepsize", accept_rate
    !     end if
        
    !     zodi_model_samp = current_model
    ! end subroutine



    subroutine downsamp_invariant_structs(cpar)
        ! Downsamples pointing, tod and caches the zodi mask (sky mask + flags)
        type(comm_params), intent(in) :: cpar

        integer(i4b) :: i, j, k, scan, npix, nmaps, ndelta, ext(2), upper_bound, padding, ierr, ntod, ndet, nhorn, ndownsamp
        real(dp) :: box_width
        real(sp), allocatable :: tod(:), mask(:)
        integer(i4b), allocatable :: pix(:, :), psi(:, :), flag(:)
        real(dp), allocatable, dimension(:, :) :: m_buf
        integer(i4b), allocatable, dimension(:) :: downsamp_pix
        logical(lgt), allocatable, dimension(:) :: downsamp_mask_idx
        real(sp), allocatable, dimension(:) :: downsamp_mask, downsamp_tod
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
            npix = 12 * data(i)%map%info%nside**2
            nmaps = data(i)%map%info%nmaps
            ndet = data(i)%tod%ndet
            nhorn = data(i)%tod%nhorn

            ! Build zodi mask
            allocate(m_buf(0:npix-1, nmaps))
            allocate(procmask_zodi(0:npix-1))
            call data(i)%tod%procmask_zodi%bcast_fullsky_map(m_buf); procmask_zodi = m_buf(:,1)

            box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)
            do scan = 1, data(i)%tod%nscan
                if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle

                ntod = data(i)%tod%scans(scan)%ntod
                ! Allocate temporary variables to be downsampled
                allocate(pix(ntod, nhorn))
                allocate(psi(ntod, nhorn))
                allocate(flag(ntod))
                allocate(tod(ntod))
                allocate(mask(ntod))

                ! Get shape of downsampled structs
                call data(i)%tod%downsample_tod(tod, ext, step=box_width)

                ! Allocate downsampled strucs
                allocate(downsamp_tod(ext(1):ext(2)))
                allocate(downsamp_mask(ext(1):ext(2)))
                allocate(downsamp_pix(ext(1):ext(2)))
                allocate(downsamp_mask_idx(ext(1):ext(2)))

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
                    do k = 1, int(size(tod) / box_width)
                        downsamp_pix(k) = pix((k * box_width), 1)
                    end do
                    where (downsamp_mask < 1.) downsamp_mask_idx = .false.

                    ! Find upper bound and truncate the downsampled arrays by where they were padded and
                    ! filter out the masked values. Store the result in the tod objects for use in sampling.
                    upper_bound = ubound(downsamp_mask, dim=1)
                    
                    data(i)%tod%scans(scan)%d(j)%downsamp_pix = pack(downsamp_pix(padding:upper_bound-padding), downsamp_mask_idx(padding:upper_bound-padding))
                    data(i)%tod%scans(scan)%d(j)%downsamp_tod = pack(downsamp_tod(padding:upper_bound-padding), downsamp_mask_idx(padding:upper_bound-padding))

                    ! Allocate other downsampled quantities with same shape 
                    ndownsamp = size(data(i)%tod%scans(scan)%d(j)%downsamp_pix)
                    allocate(data(i)%tod%scans(scan)%d(j)%downsamp_sky(ndownsamp))
                    allocate(data(i)%tod%scans(scan)%d(j)%downsamp_zodi(ndownsamp))
                    allocate(data(i)%tod%scans(scan)%d(j)%downsamp_scat(ndownsamp, zodi_model%n_comps))
                    allocate(data(i)%tod%scans(scan)%d(j)%downsamp_therm(ndownsamp, zodi_model%n_comps))
                end do
                deallocate(pix, psi, flag, tod, mask, downsamp_tod, downsamp_mask, downsamp_pix, downsamp_mask_idx)
            end do

            deallocate(m_buf, procmask_zodi)
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
            allocate(sky_signal(data(i)%tod%ndet, ndelta))
            do j = 1, data(i)%tod%ndet
                call get_sky_signal(i, j, sky_signal(j, ndelta)%p, mono=.false.)
            end do

            allocate(map_sky(nmaps, data(i)%tod%nobs, 0:data(i)%tod%ndet, 1))
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
            deallocate(map_sky, sky_signal)
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
end module
