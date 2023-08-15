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
    public initialize_zodi_samp_mod, sample_zodi_model, init_scandata_and_downsamp_zodi, gibbs_sample_zodi_comp, sample_zodi_emissivity_and_albedo, sample_zodi_parameter, active_params

    ! globals
    real(dp), allocatable :: chisq_previous, step_size, step_sizes(:), priors(:, :)
    logical(lgt), allocatable :: active_params(:)
    integer(i4b) :: n_active_params, cloud_start, cloud_stop, band1_start, band1_stop, band2_start, band2_stop, band3_start, band3_stop, ring_start, ring_stop, feature_start, feature_stop, emissivity_start, emissivity_stop, albedo_start, albedo_stop
    real(dp) :: NO_PRIOR = HUGE(1.0_dp)

    ! Custom types for albedo and emissivity fitting
    type :: ZodiTimestreams
        real(sp), allocatable, dimension(:, :) :: s_therm   
        real(sp), allocatable, dimension(:, :) :: s_scat    
    end type ZodiTimestreams

    type :: ZodiDetector
        type(ZodiTimestreams), allocatable, dimension(:) :: d
    end type ZodiDetector

    type :: ZodiScan
        type(ZodiDetector), allocatable, dimension(:) :: scans
    end type ZodiScan

contains
    subroutine initialize_zodi_samp_mod(cpar)
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
        active_params(1) = .true. !n_0
        active_params(2) = .true. !incl
        active_params(3) = .true. !omega
        active_params(4) = .true. !x_0
        active_params(5) = .true. !y_0
        active_params(6) = .true. !z_0
        active_params(7) = .true. !alpha
        active_params(8) = .true. !beta
        active_params(9) = .true. !gamma
        active_params(10) = .true. !mu
        cloud_start = 1 ; cloud_stop = 10

        ! Band1
        active_params(11) = .true. !n_0
        active_params(12) = .true. !incl
        active_params(13) = .true. !omega
        active_params(14) = .true. !x_0
        active_params(15) = .true. !y_0
        active_params(16) = .true. !z_0
        active_params(17) = .true. !delta_zeta
        active_params(18) = .true. !delta_r
        active_params(19) = .true. !v
        active_params(20) = .true. !p
        band1_start = 11 ; band1_stop = 20

        ! Band2
        active_params(21) = .true. !n_0
        active_params(22) = .true. !incl
        active_params(23) = .true. !omega
        active_params(24) = .true. !x_0
        active_params(25) = .true. !y_0
        active_params(26) = .true. !z_0
        active_params(27) = .true. !delta_zeta
        active_params(28) = .true. !delta_r
        active_params(29) = .true. !v
        active_params(30) = .true. !p
        band2_start = 21 ; band2_stop = 30

        ! Band3
        active_params(31) = .true. !n_0
        active_params(32) = .true. !incl
        active_params(33) = .true. !omega
        active_params(34) = .true. !x_0
        active_params(35) = .true. !y_0
        active_params(36) = .true. !z_0
        active_params(37) = .true. !delta_zeta
        active_params(38) = .true. !delta_r
        active_params(39) = .true. !v
        active_params(40) = .true. !p
        band3_start = 31 ; band3_stop = 40

        ! Ring
        active_params(41) = .true. !n_0
        active_params(42) = .true. !incl
        active_params(43) = .true. !omega
        active_params(44) = .true. !x_0
        active_params(45) = .true. !y_0
        active_params(46) = .true. !z_0
        active_params(47) = .true. !R_0
        active_params(48) = .true. !sigma_r
        active_params(49) = .true. !sigma_z
        ring_start = 41 ; ring_stop = 49

        ! Feature
        active_params(50) = .true. !n_0
        active_params(51) = .true. !incl
        active_params(52) = .true. !omega
        active_params(53) = .true. !x_0
        active_params(54) = .true. !y_0
        active_params(55) = .true. !z_0
        active_params(56) = .true. !R_0
        active_params(57) = .true. !sigma_r
        active_params(58) = .true. !sigma_z
        active_params(59) = .true. !theta_0
        active_params(60) = .true. !sigma_theta
        feature_start = 50 ; feature_stop = 60

        ! Other
        active_params(61) = .true. !T_0
        active_params(62) = .true. !delta

        ! emissivities and albedos
        emissivity_start = 63 ; emissivity_stop = 63 + n_spec_params - 1
        albedo_start = 63 + n_spec_params ; albedo_stop = 63 + 2 * n_spec_params - 1

        active_params(emissivity_start:emissivity_stop) = .true.
        active_params(albedo_start:albedo_stop) = .true.

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
        priors(7, :) = [0., 10.] !alpha
        priors(8, :) = [0., 10.] !beta
        priors(9, :) = [0., 10.] !gamma
        priors(10, :) = [0., 10.] !mu

        ! Band1
        priors(11, :) = [5.59d-12, 5.59d-8] !n_0
        priors(12, :) = NO_PRIOR !incl
        priors(13, :) = NO_PRIOR !omega
        priors(14, :) = [-0.25, 0.25] !x_0
        priors(15, :) = [-0.25, 0.25] !y_0
        priors(16, :) = [-0.25, 0.25] !z_0
        priors(17, :) = [0., 50.] !delta_zeta
        priors(18, :) = [0.5, 6.] !delta_r
        priors(19, :) = [0., 10.] !v
        priors(20, :) = [0., 10.] !p

        ! Band2
        priors(21, :) = [1.99d-11, 1.99d-7] !n_0
        priors(22, :) = NO_PRIOR !incl
        priors(23, :) = NO_PRIOR !omega
        priors(24, :) = [-0.25, 0.25] !x_0
        priors(25, :) = [-0.25, 0.25] !y_0
        priors(26, :) = [-0.25, 0.25] !z_0
        priors(27, :) = [0., 50.] !delta_zeta
        priors(28, :) = [0.5, 6.] !delta_r
        priors(29, :) = [0., 10.] !v
        priors(30, :) = [0., 10.] !p

        ! Band3
        priors(31, :) = [1.44d-12, 1.44d-8] !n_0
        priors(32, :) = NO_PRIOR !incl
        priors(33, :) = NO_PRIOR !omega
        priors(34, :) = [-0.25, 0.25] !x_0
        priors(35, :) = [-0.25, 0.25] !y_0
        priors(36, :) = [-0.25, 0.25] !z_0
        priors(37, :) = [0., 50.] !delta_zeta
        priors(38, :) = [0., 6.] !delta_r
        priors(39, :) = [0., 10.] !v
        priors(40, :) = [0., 10.] !p

        ! Ring
        priors(41, :) = [1.83d-10, 1.83d-6] !n_0
        priors(42, :) = NO_PRIOR !incl
        priors(43, :) = NO_PRIOR !omega
        priors(44, :) = [-0.25, 0.25] !x_0
        priors(45, :) = [-0.25, 0.25] !y_0
        priors(46, :) = [-0.25, 0.25] !z_0
        priors(47, :) = [0.5, 2.] !R_0
        priors(48, :) = [0., 5.] !sigma_r
        priors(49, :) = [0., 5.] !sigma_z

        ! Feature
        priors(50, :) = [1.9d-10, 1.9d-6] !n_0
        priors(51, :) = NO_PRIOR !incl
        priors(52, :) = NO_PRIOR !omega
        priors(53, :) = [-0.25, 0.25] !x_0
        priors(54, :) = [-0.25, 0.25] !y_0
        priors(55, :) = [-0.25, 0.25] !z_0
        priors(56, :) = [0.5, 2.] !R_0
        priors(57, :) = [0., 5.] !sigma_r
        priors(58, :) = [0., 5.] !sigma_z
        priors(59, :) = NO_PRIOR !theta_0
        priors(60, :) = [0., 100.] !sigma_theta

        ! Other
        priors(61, :) = [200., 400.] !T_0
        priors(62, :) = [0.1, 1.] !delta

        ! Emissivity
        priors(emissivity_start:emissivity_stop, 1) = 0.
        priors(emissivity_start:emissivity_stop, 2) = 5.

        ! Albedo
        priors(albedo_start:albedo_stop, 1) = 0.
        priors(albedo_start:albedo_stop, 2) = 1.


    end subroutine initialize_zodi_samp_mod

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


    subroutine sample_zodi_model(cpar, handle)
        ! Metropolis-Hastings for nproposals of new sets of zodi parameters
        ! Todo: Split into gibbs steps for each component
        type(comm_params), intent(in) :: cpar
        type(planck_rng), intent(inout) :: handle
        integer(i4b) :: i, j, k, flag, ndet, nscan, ntod, nprop, scan, ierr, n_accepted, n_tot_tod, n_tot_tod_reduced, n_tot_tod_current
        real(sp), allocatable :: s_scat(:, :), s_therm(:, :), s_zodi(:)
        real(dp) :: chisq_tod, chisq_current, chisq_diff, ln_acceptance_probability, accept_rate, chisq_lnL
        real(dp), allocatable :: param_vec(:), powell_vec(:)
        type(ZodiModel) :: current_model, previous_model
        logical(lgt) :: accepted, operation
        ! real(dp) :: n0_grid(100), chisq_grid(100)

        !debugging stuff
        character(len=256) :: chaindir
        character(len=512) :: path, param_path, param_name
        character(len=6) :: itext
        type(hdf_file) :: tod_file
        integer(i4b) :: l

        ! Metropolis-Hastings for nproposals of new sets of zodi parameters
        nprop = cpar%zs_nprop
        n_accepted = 0
        current_model = zodi_model
        previous_model = zodi_model

        allocate(param_vec(current_model%n_parameters))
        allocate(powell_vec(current_model%n_comps))  

        ! param_vec = current_model%model_to_param_vec()
        ! powell_vec = get_powell_vec(param_vec)
        
        ! Start initial parameter minimization. Root will aquire new set of parameters. While this happens
        ! the other processes will wait for the new parameters to be broadcasted. Then they will all compute
        ! the their own lnL.
        ! if (cpar%myid == cpar%root) then
        !     print *, "before", powell_vec
        !     call powell(powell_vec, lnL_zodi, ierr)
        !     flag = 0
        !     call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
        ! else
        !     do while (.true.)
        !         call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
        !         if (flag == 1) then 
        !             ! Here we call the function without a powell vector, because it will recieve it in internally from MPI
        !             chisq_lnL = lnL_zodi()
        !             ! print *, chisq_lnL
        !         else
        !             exit
        !         end if
        !     end do
        ! end if

        ! ! Wait for every process to exit
        ! call mpi_barrier(MPI_COMM_WORLD, ierr)

        ! call mpi_bcast(powell_vec, size(powell_vec), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
        ! if (cpar%myid == cpar%root) then
        !     print *, "after", powell_vec, "ierr", ierr
        ! end if
        ! stop

        ! ! update param_vec with new parameters from powell
        ! call update_param_vec_from_powell_vec(param_vec, powell_vec)

        ! update model with new param_vec
        ! call current_model%param_vec_to_model(param_vec)

        do k = 0, nprop
            ! Reset chisq for current proposal
            chisq_tod = 0.
            chisq_current = 0.
            
            ! Root process draws new set of zodi parameters and broadcasts to all processes
            ! If first iteration we dont want to draw, just compute the chisq with the base model
            if (k > 0) then
                if (cpar%myid == cpar%root) then
                    param_vec = current_model%model_to_param_vec()
                    do i = 1, size(param_vec)
                        if (active_params(i) == .true.) then
                            param_vec(i) = param_vec(i) + (step_sizes(i) * param_vec(i) * rand_gauss(handle))
                            if (priors(i, 1) == NO_PRIOR .or. priors(i, 2) == NO_PRIOR) then 
                                continue
                            else if (param_vec(i) < priors(i, 1) .or. param_vec(i) > priors(i, 2)) then 
                                step_sizes(i) = step_sizes(i) / 2.
                                chisq_tod = 1.d30
                                exit
                            end if
                        end if
                    end do
                end if
                call mpi_bcast(param_vec, size(param_vec), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
                call current_model%param_vec_to_model(param_vec)
            end if

            do i = 1, numband
                ! Skip none tod bands
                if (data(i)%tod_type == "none") cycle
                ! Skip tod bands where we dont want to sample zodi
                if (.not. data(i)%tod%sample_zodi) cycle

                if (chisq_tod >= 1.d30) exit

                ndet = data(i)%tod%ndet
                nscan = data(i)%tod%nscan
                
                ! Make sure that the zodi cache is cleared before each new band
                call data(i)%tod%clear_zodi_cache()

                ! Evaluate zodi model with newly proposed values for each band and calculate chisq
                do scan = 1, nscan
                    ! Skip scan if no accepted data
                    do j = 1, ndet
                        if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
                        ntod = size(data(i)%tod%scans(scan)%d(j)%downsamp_res)

                        allocate(s_scat(ntod, zodi_model%n_comps), s_therm(ntod, zodi_model%n_comps), s_zodi(ntod))
                        call get_zodi_emission(&
                            & tod=data(i)%tod, &
                            & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pointing, &
                            & scan=scan, &
                            & det=j, &
                            & s_zodi_scat=s_scat, &
                            & s_zodi_therm=s_therm, &
                            & model=current_model &
                        &)
                        call get_s_zodi(&
                            & s_therm=s_therm, &
                            & s_scat=s_scat, &
                            & s_zodi=s_zodi, &
                            & band=data(i)%tod%band, &
                            & model=current_model &
                        &)
                    
                        chisq_tod = chisq_tod + sum(((data(i)%tod%scans(scan)%d(j)%downsamp_res - s_zodi)/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2)
                        deallocate(s_scat, s_therm, s_zodi)
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
                    ln_acceptance_probability = -0.5 * chisq_diff
                    accepted = ln_acceptance_probability > log(rand_uni(handle))
                    
                    print *, "proposal:", k
                    print *, "chisq_current:", chisq_current
                    print *, "chisq_diff:", chisq_diff
                    print *, "ln a:", ln_acceptance_probability 
                    print *, "accepted:", accepted
                    print *, "accept rate:", (real(n_accepted) / real(k))*100.
                    print *, " "
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
                accept_rate = n_accepted / k
            end if
        end do  
        if (accept_rate < 10.) then
            step_sizes = step_sizes / 2.
            print *, "lowering stepsize"
        else if (accept_rate > 60.) then 
            step_sizes = step_sizes * 2.
            print *, "increasing stepsize"
        end if
        
        zodi_model = current_model
    end subroutine


    subroutine init_scandata_and_downsamp_zodi(cpar, comp_idx)
        ! This routine mimics parts of the process_TOD where we calibrate the timestreams 
        ! put them into th scan_data type, and then downsample them for use in zodi fitting.
        ! See the process_TOD routines for more details on what happens here.
        type(comm_params), intent(in) :: cpar
        integer(i4b), intent(in), optional :: comp_idx

        integer(i4b) :: i, j , k, nside, npix, nmaps, ndelta, box_width
        real(sp), allocatable, dimension(:, :, :, :) :: map_sky
        real(sp), allocatable, dimension(:) :: procmask, procmask2, procmask_zodi
        real(dp), allocatable, dimension(:, :) :: m_buf
        integer(i4b), allocatable, dimension(:) :: downsamp_width
        logical(lgt) :: downsamp_all_equally
        type(comm_scandata) :: sd
        type(map_ptr), allocatable, dimension(:, :) :: s_sky


        if (present(comp_idx)) then
            allocate(downsamp_width(zodi_model%n_comps))
            downsamp_all_equally = .false.
        else 
            downsamp_all_equally = .true.
        end if

        ! For each zodi band, create downsampled residual time stream
        do i = 1, numband

            ! Only generate downsampled arrays for tod bands and bands where we have zodi
            if (trim(data(i)%tod_type) == 'none') cycle
            if (.not. data(i)%tod%subtract_zodi) cycle

            ! Clear previously allocated arrays if they exist
            call data(i)%tod%deallocate_downsampled_zodi()
            call data(i)%tod%clear_zodi_cache()

            ndelta = 1
            nside = data(i)%map%info%nside
            npix = 12 * nside**2
            nmaps = data(i)%map%info%nmaps

            allocate(m_buf(0:npix-1, nmaps), procmask(0:npix-1), procmask2(0:npix-1), procmask_zodi(0:npix-1))
            call data(i)%tod%procmask%bcast_fullsky_map(m_buf); procmask  = m_buf(:,1)
            call data(i)%tod%procmask2%bcast_fullsky_map(m_buf); procmask2 = m_buf(:,1)
            call data(i)%tod%procmask_zodi%bcast_fullsky_map(m_buf); procmask_zodi = m_buf(:,1)
            
            allocate(s_sky(data(i)%tod%ndet, ndelta))
            do j = 1, data(i)%tod%ndet
                call get_sky_signal(i, j, s_sky(j, ndelta)%p, mono=.false.)
            end do

            allocate(map_sky(nmaps, data(i)%tod%nobs, 0:data(i)%tod%ndet, 1))
            call distribute_sky_maps(data(i)%tod, s_sky, 1.e0, map_sky)

            ! Start with same downsamp box width for all comps (change in future)
            if (downsamp_all_equally) then
                box_width = data(i)%tod%samprate / data(i)%tod%samprate_lowres
            else 
                stop "not implemented"
            end if

            do j = 1, data(i)%tod%nscan
                if (.not. any(data(i)%tod%scans(j)%d%accept)) cycle

                ! not calibrating the data so we use substitute map_gain with map_sky
                call sd%init_singlehorn(data(i)%tod, j, map_sky, map_sky, procmask, procmask2, procmask_zodi)
                call boxcar_downsamp_zodi_res_and_pointing(data(i)%tod, sd, j, padding=5, box_width=box_width)
                ! prepare scandata object for next scan
                call sd%dealloc
            end do

            deallocate(m_buf, procmask, procmask2, procmask_zodi)
            deallocate(s_sky)
            deallocate(map_sky)
        end do
    end subroutine


    subroutine boxcar_downsamp_zodi_res_and_pointing(tod, sd, scan_id, padding, box_width)
        ! Downsamples the tod residual and pointing given a box car width for use in zodi sampling.
        class(comm_tod), intent(inout) :: tod
        type(comm_scandata), intent(in) :: sd
        integer(i4b), intent(in) :: scan_id, padding
        integer(i4b), intent(in) :: box_width

        integer(i4b) :: i, j, ext(2), upper_bound, box_halfwidth
        real(dp) :: box_fullwidth
        real(sp), allocatable, dimension(:) :: res, res_truncated, downsamp_mask, downsamp_res, downsamp_pointing
        logical(lgt), allocatable, dimension(:) :: downsamp_masked_indices

        if (box_width < 1.) stop "Cannot downsample zodi tods if box car width is less than 1 sample!"

       ! Make box have odd size so that we can pick out a center pixel
        if (mod(box_width, 2) == 0) then 
            box_fullwidth = box_width + 1.
        else 
            box_fullwidth = box_width
        end if
        box_halfwidth = box_fullwidth / 2

        allocate(res(size(sd%s_tot, dim=1)))

        do j = 1, tod%ndet
            ! Add zodi back to the residual and apply mask
            res = (sd%tod(:, j) - sd%s_tot(:,j) + sd%s_zodi(:, j))

            ! Get downsampled shape (ext), and allocate the downsampled arrays
            call tod%downsample_tod(res, ext, step=box_fullwidth)

            ! Allocate these the first gibbs iter
            allocate(downsamp_res(ext(1):ext(2)))
            allocate(downsamp_pointing(ext(1):ext(2)))
            allocate(downsamp_mask(ext(1):ext(2)))
            allocate(downsamp_masked_indices(ext(1):ext(2)))

            downsamp_mask = 0.
            downsamp_res = 0.
            downsamp_pointing = 0

            ! Downsample the mask
            call tod%downsample_tod(sd%mask_zodi(:, j), ext, downsamp_mask, step=box_fullwidth)

            ! Downsample the residual
            call tod%downsample_tod(res, ext, downsamp_res, step=box_fullwidth)

            ! Construct the downsampled pointing array (pick central pixel in each box)
            do i = 1, int(size(res) / box_fullwidth)
                downsamp_pointing(i) = sd%pix((i * box_fullwidth), j, 1)
            end do

            ! Remove flagged values in the downsampled residual
            downsamp_masked_indices = .true.
            where (downsamp_mask < 1.) downsamp_masked_indices = .false.
            
            ! Find upper bound and truncate the downsampled arrays by where they were padded and
            ! filter out the masked values. Store the result in the tod objects for use in sampling.
            upper_bound = ubound(downsamp_res, dim=1)
            tod%scans(scan_id)%d(j)%downsamp_res = pack(downsamp_res(padding:upper_bound-padding), downsamp_masked_indices(padding:upper_bound-padding))
            tod%scans(scan_id)%d(j)%downsamp_pointing = pack(downsamp_pointing(padding:upper_bound-padding), downsamp_masked_indices(padding:upper_bound-padding))
            deallocate(downsamp_res, downsamp_pointing, downsamp_mask, downsamp_masked_indices)
        end do
    end subroutine boxcar_downsamp_zodi_res_and_pointing


    function lnL_zodi(p)
        use healpix_types
        implicit none
        real(dp), dimension(:), intent(in), optional :: p
        real(dp) :: lnL_zodi

        real(dp), allocatable :: theta(:)
        real(sp), allocatable :: s_scat(:, :), s_therm(:, :), s_zodi(:)
        type(ZodiModel) :: model

        integer(i4b) :: i, j, ntod, ndet, nscan, scan, ierr, flag
        real(dp) :: chisq_tod, chisq, chisq_buff
        real(dp), allocatable :: p_copy(:), param_vec(:)
        allocate(theta(n_active_params))

        if (data(1)%tod%myid == 0) then
            flag = 1
            call mpi_bcast(flag, 1, MPI_INTEGER, 0, data(1)%tod%comm, ierr)
            theta = p
        end if
        call mpi_bcast(theta, size(theta), MPI_DOUBLE_PRECISION, 0, data(1)%tod%comm, ierr)

        ! Get param vector from zodi model
        model = zodi_model
        allocate(param_vec(model%n_parameters))
        param_vec = model%model_to_param_vec()
        
        ! Check priors
        j = 1
        do i = 1, size(param_vec)
            if (active_params(i) == .true.) then
                if (theta(j) < priors(i, 1) .or. theta(j) > priors(i, 2)) then
                    lnL_zodi = 1.d30
                    return
                end if
                j = j + 1
            end if
        end do

        ! Update param_vec with new values from powell_vec (theta)
        call update_param_vec_from_powell_vec(param_vec, theta)

        ! update zodi model with new param_vec
        call model%param_vec_to_model(param_vec)


        ! Calculate chisq for each band
        chisq_tod = 0.
        chisq = 0.
        lnL_zodi = 0.

        do i = 1, numband
            ! Skip none tod bands
            if (data(i)%tod_type == "none") cycle
            ! Skip tod bands where we dont want to sample zodi
            if (.not. data(i)%tod%sample_zodi) cycle

            ndet = data(i)%tod%ndet
            nscan = data(i)%tod%nscan
            
            ! Make sure that the zodi cache is cleared before each new band
            call data(i)%tod%clear_zodi_cache()

            ! Evaluate zodi model with newly proposed values for each band and calculate chisq
            do scan = 1, nscan
                ! Skip scan if no accepted data
                do j = 1, ndet
                    if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
                    ntod = size(data(i)%tod%scans(scan)%d(j)%downsamp_res)
                    allocate(s_scat(ntod, model%n_comps), s_therm(ntod, model%n_comps), s_zodi(ntod))
                    call get_zodi_emission(&
                        & tod=data(i)%tod, &
                        & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pointing, &
                        & scan=scan, &
                        & det=j, &
                        & s_zodi_scat=s_scat, &
                        & s_zodi_therm=s_therm, &
                        & model=model &
                    &)
                    call get_s_zodi(&
                        & s_therm=s_therm, &
                        & s_scat=s_scat, &
                        & s_zodi=s_zodi, &
                        & band=data(i)%tod%band, &
                        & model=model &
                    &)
                    chisq_buff = sum(((data(i)%tod%scans(scan)%d(j)%downsamp_res - s_zodi)/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2)

                    chisq_tod = chisq_tod + chisq_buff
                    deallocate(s_scat, s_therm, s_zodi)
                end do
            end do
        end do

        ! Reduce chisq to root process
        call mpi_reduce(chisq_tod, chisq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, data(1)%tod%comm, ierr)

        if (data(1)%tod%myid == 0) then
            ! print *, chisq, model%comps(1)%c%n_0
            lnL_zodi = chisq
        end if

    end function lnL_zodi


    subroutine sample_zodi_parameter(cpar, handle, gibbs_iter, param_idx, zodi_model_samp, verbose)
        ! Gibbs sample one by one zodi parameter
        type(comm_params), intent(in) :: cpar
        type(planck_rng), intent(inout) :: handle
        integer(i4b), intent(in) :: gibbs_iter
        integer(i4b), intent(in) :: param_idx
        type(ZodiModel), intent(inout) :: zodi_model_samp
        logical(lgt), intent(in), optional :: verbose

        integer(i4b) :: i, j, k, param_start, param_stop, flag, ndet, nscan, ntod, nprop, scan, ierr, n_accepted, n_tot_tod, n_tot_tod_reduced, n_tot_tod_current
        real(sp), allocatable :: s_scat(:, :), s_therm(:, :), s_zodi(:)
        real(dp) :: chisq_tod, chisq_current, chisq_diff, ln_acceptance_probability, accept_rate, chisq_lnL
        real(dp), allocatable :: param_vec(:), powell_vec(:)
        type(ZodiModel) :: current_model, previous_model
        logical(lgt) :: accepted, verbose_, tuning, skip
        type(hdf_file) :: tod_file

        if (present(verbose)) then
            verbose_ = verbose
        else
            verbose_ = .false.
        end if

        if (verbose_ .and. (cpar%myid == cpar%root)) then
            print *, "==== Zodi gibbs sampling step for param number: ", param_idx, " ===="
            print *, " "
        end if 

        ! Metropolis-Hastings for nproposals of new sets of zodi parameters
        nprop = cpar%zs_nprop

        ! Make two copies of the zodi model, one for the current proposal and one for the previous proposal
        current_model = zodi_model_samp
        previous_model = zodi_model_samp

        allocate(param_vec(current_model%n_parameters))

        if (gibbs_iter >= 25) then
            tuning = .false.
        else
            tuning = .true.
        end if
        n_accepted = 0

        do k = 0, nprop
            ! if propsed parameter is outside of priors we dont want to immediately reject
            skip = .false.

            ! Reset chisq for current proposal
            chisq_tod = 0.
            chisq_current = 0.
            ! If first iteration we dont want to draw, just compute the chisq with the base model
            if (k > 0) then

                ! Root process draws new set of zodi parameters and broadcasts to all processes
                param_vec = current_model%model_to_param_vec()
                if (cpar%myid == cpar%root) then
                    param_vec(param_idx) = param_vec(param_idx) + (step_sizes(param_idx) * rand_gauss(handle))

                    ! If parameter is outside of prior range, immediatly reject proposal and halve step size
                    if (priors(param_idx, 1) == NO_PRIOR .or. priors(param_idx, 2) == NO_PRIOR) then 
                        continue
                    else if (param_vec(param_idx) < priors(param_idx, 1) .or. param_vec(param_idx) > priors(param_idx, 2)) then 
                        if (tuning) step_sizes(param_idx) = step_sizes(param_idx) / 2.
                        chisq_tod = 1.d30
                        skip = .true.
                    end if
                end if
                call mpi_bcast(param_vec, size(param_vec), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
                ! Update current_model with new parameters
                call current_model%param_vec_to_model(param_vec)
            end if

            ! loop over downsampled residuals for each tod band with zodi
            do i = 1, numband
                if (skip) exit ! prior was violated, skip rest of evaluation
                if (data(i)%tod_type == "none") cycle
                if (.not. data(i)%tod%sample_zodi) cycle

                ! If chisq is already too large, skip rest of the evaluation and go directly to rejection
                if (chisq_tod >= 1.d30) exit

                ndet = data(i)%tod%ndet
                nscan = data(i)%tod%nscan
                
                ! Make sure that the zodi cache is cleared before each new band
                call data(i)%tod%clear_zodi_cache()

                ! Evaluate zodi model with newly proposed values for each band and calculate chisq
                do scan = 1, nscan
                    ! Skip scan if no accepted data
                    do j = 1, ndet
                        if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
                        ntod = size(data(i)%tod%scans(scan)%d(j)%downsamp_res)

                        allocate(s_scat(ntod, zodi_model%n_comps), s_therm(ntod, zodi_model%n_comps), s_zodi(ntod))
                        call get_zodi_emission(&
                            & tod=data(i)%tod, &
                            & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pointing, &
                            & scan=scan, &
                            & det=j, &
                            & s_zodi_scat=s_scat, &
                            & s_zodi_therm=s_therm, &
                            & model=current_model &
                        &)
                        call get_s_zodi(&
                            & s_therm=s_therm, &
                            & s_scat=s_scat, &
                            & s_zodi=s_zodi, &
                            & band=data(i)%tod%band, &
                            & model=current_model &
                        &)

                        if (.false.) then
                            if (cpar%myid == cpar%root) then
                                call open_hdf_file(trim("/mn/stornext/d5/data/metins/dirbe/chains/chains_testing/downsamp_timestream.h5"), tod_file, 'w')
                                call write_hdf(tod_file, '/r', data(i)%tod%scans(scan)%d(j)%downsamp_res)
                                call write_hdf(tod_file, '/z', s_zodi)
                                call close_hdf_file(tod_file)
                                stop
                            end if
                            call mpi_barrier(cpar%comm_chain, ierr)
                            stop
                        end if
                        
                        chisq_tod = chisq_tod + sum(((data(i)%tod%scans(scan)%d(j)%downsamp_res - s_zodi)/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2)
                        deallocate(s_scat, s_therm, s_zodi)
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
                    ln_acceptance_probability = -0.5 * chisq_diff
                    select case (cpar%operation)
                    case ("sample")
                        accepted = ln_acceptance_probability > log(rand_uni(handle))
                    case ("optimize")
                        accepted = chisq_diff < 0.
                    case default
                        stop "Error: invalid operation"
                    end select

                    if (verbose_) then
                        print *, "proposal:", k
                        print *, "chisq_current:", chisq_current
                        print *, "chisq_diff:", chisq_diff
                        print *, "ln a:", ln_acceptance_probability 
                        print *, "accepted:", accepted
                        print *, "accept rate:", (real(n_accepted) / real(k))*100.
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
                accept_rate = real(n_accepted) / real(k)
            end if
        end do  
        if (accept_rate < 0.1) then
            step_sizes(param_idx) = step_sizes(param_idx) / 2.
            if (verbose_ .and. (cpar%myid == cpar%root)) print *, "lowering stepsize", accept_rate
        else if (accept_rate > 0.6) then 
            step_sizes(param_idx) = step_sizes(param_idx) * 2.
            if (verbose_ .and. (cpar%myid == cpar%root)) print *, "increasing stepsize", accept_rate
        end if
        
        zodi_model_samp = current_model
    end subroutine


    subroutine gibbs_sample_zodi_comp(cpar, handle, comp_idx, zodi_model_samp, verbose)
        type(comm_params), intent(in) :: cpar
        type(planck_rng), intent(inout) :: handle
        integer(i4b), intent(in) :: comp_idx
        type(ZodiModel), intent(inout) :: zodi_model_samp
        logical(lgt), intent(in), optional :: verbose

        integer(i4b) :: i, j, k, param_start, param_stop, flag, ndet, nscan, ntod, nprop, scan, ierr, n_accepted, n_tot_tod, n_tot_tod_reduced, n_tot_tod_current
        real(sp), allocatable :: s_scat(:, :), s_therm(:, :), s_zodi(:)
        real(dp) :: chisq_tod, chisq_current, chisq_diff, ln_acceptance_probability, accept_rate, chisq_lnL
        real(dp), allocatable :: param_vec(:), powell_vec(:)
        type(ZodiModel) :: current_model, previous_model
        logical(lgt) :: accepted, verbose_
        type(hdf_file) :: tod_file

        if (present(verbose)) then
            verbose_ = verbose
        else
            verbose_ = .false.
        end if

        ! Get parameter indicies for this component
        select case (trim(zodi_model%comp_labels(comp_idx))) 
        case ("cloud")
            param_start = cloud_start
            param_stop = cloud_stop
        case ("band1")
            param_start = band1_start
            param_stop = band1_stop
        case ("band2")
            param_start = band2_start
            param_stop = band2_stop
        case ("band3")
            param_start = band3_start
            param_stop = band3_stop
        case ("ring")
            param_start = ring_start
            param_stop = ring_stop
        case ("feature")
            param_start = feature_start
            param_stop = feature_stop
        case default 
            stop "Error: invalid component label"
        end select

        if (verbose_ .and. (cpar%myid == cpar%root)) then
            print *, "==== Zodi gibbs sampling step: ", trim(zodi_model%comp_labels(comp_idx)), " ===="
            print *, " "
        end if 

        ! Metropolis-Hastings for nproposals of new sets of zodi parameters
        nprop = cpar%zs_nprop

        ! Make two copies of the zodi model, one for the current proposal and one for the previous proposal
        current_model = zodi_model_samp
        previous_model = zodi_model_samp

        allocate(param_vec(current_model%n_parameters))

        n_accepted = 0
        do k = 0, nprop
            ! Reset chisq for current proposal
            chisq_tod = 0.
            chisq_current = 0.
            
            ! If first iteration we dont want to draw, just compute the chisq with the base model
            if (k > 0) then

                ! Root process draws new set of zodi parameters and broadcasts to all processes
                if (cpar%myid == cpar%root) then
                    param_vec = current_model%model_to_param_vec()

                    ! propse new set of parameters for component
                    do i = param_start, param_stop
                        if (.not. active_params(i)) cycle
                        param_vec(i) = param_vec(i) + (step_sizes(i) * rand_gauss(handle))

                        ! If parameter is outside of prior range, immediatly reject proposal and halve step size
                        if (priors(i, 1) == NO_PRIOR .or. priors(i, 2) == NO_PRIOR) then 
                            continue
                        else if (param_vec(i) < priors(i, 1) .or. param_vec(i) > priors(i, 2)) then 
                            step_sizes(i) = step_sizes(i) / 2.
                            chisq_tod = 1.d30
                            exit
                        end if
                    end do
                end if
                call mpi_bcast(param_vec, size(param_vec), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
                ! Update current_model with new parameters
                call current_model%param_vec_to_model(param_vec)
            end if

            ! loop over downsampled residuals for each tod band with zodi
            do i = 1, numband
                if (data(i)%tod_type == "none") cycle
                if (.not. data(i)%tod%sample_zodi) cycle

                ! If chisq is already too large, skip rest of the evaluation and go directly to rejection
                if (chisq_tod >= 1.d30) exit

                ndet = data(i)%tod%ndet
                nscan = data(i)%tod%nscan
                
                ! Make sure that the zodi cache is cleared before each new band
                call data(i)%tod%clear_zodi_cache()

                ! Evaluate zodi model with newly proposed values for each band and calculate chisq
                do scan = 1, nscan
                    ! Skip scan if no accepted data
                    do j = 1, ndet
                        if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
                        ntod = size(data(i)%tod%scans(scan)%d(j)%downsamp_res)

                        allocate(s_scat(ntod, zodi_model%n_comps), s_therm(ntod, zodi_model%n_comps), s_zodi(ntod))
                        call get_zodi_emission(&
                            & tod=data(i)%tod, &
                            & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pointing, &
                            & scan=scan, &
                            & det=j, &
                            & s_zodi_scat=s_scat, &
                            & s_zodi_therm=s_therm, &
                            & model=current_model &
                        &)
                        call get_s_zodi(&
                            & s_therm=s_therm, &
                            & s_scat=s_scat, &
                            & s_zodi=s_zodi, &
                            & band=data(i)%tod%band, &
                            & model=current_model &
                        &)
                        ! if (.true. .and. cpar%myid == cpar%root) then
                        !     call open_hdf_file(trim("/mn/stornext/d5/data/metins/dirbe/chains/chains_testing/downsamp_timestream.h5"), tod_file, 'w')
                        !     call write_hdf(tod_file, '/r', data(i)%tod%scans(scan)%d(j)%downsamp_res)
                        !     call write_hdf(tod_file, '/p', data(i)%tod%scans(scan)%d(j)%downsamp_pointing)
                        !     call write_hdf(tod_file, '/z', s_zodi)
                        !     call close_hdf_file(tod_file)
                        !     stop
                        ! end if
                        chisq_tod = chisq_tod + sum(((data(i)%tod%scans(scan)%d(j)%downsamp_res - s_zodi)/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2)
                        deallocate(s_scat, s_therm, s_zodi)
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
                    ln_acceptance_probability = -0.5 * chisq_diff
                    accepted = ln_acceptance_probability > log(rand_uni(handle))
                    if (verbose_) then
                        print *, "proposal:", k
                        print *, "chisq_current:", chisq_current
                        print *, "chisq_diff:", chisq_diff
                        print *, "ln a:", ln_acceptance_probability 
                        print *, "accepted:", accepted
                        print *, "accept rate:", (real(n_accepted) / real(k))*100.
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
                accept_rate = real(n_accepted) / real(k)
            end if
        end do  
        if (accept_rate < 0.1) then
            step_sizes(param_start:param_stop) = step_sizes(param_start:param_stop) / 2.
            if (verbose_ .and. (cpar%myid == cpar%root)) print *, "lowering stepsize", accept_rate
        else if (accept_rate > 0.6) then 
            step_sizes(param_start:param_stop) = step_sizes(param_start:param_stop) * 2.
            if (verbose_ .and. (cpar%myid == cpar%root)) print *, "increasing stepsize", accept_rate
        end if
        
        zodi_model_samp = current_model
    end subroutine


    subroutine sample_zodi_emissivity_and_albedo(cpar, handle, gibbs_iter, zodi_model_samp, verbose)
        type(comm_params), intent(in) :: cpar
        type(planck_rng), intent(inout) :: handle
        integer(i4b), intent(in) :: gibbs_iter
        type(ZodiModel), intent(inout) :: zodi_model_samp
        logical(lgt), intent(in), optional :: verbose

        integer(i4b) :: i, j, k, param_start, param_stop, flag, ndet, nscan, ntod, nprop, scan, ierr, n_accepted, n_tot_tod, n_tot_tod_reduced, n_tot_tod_current
        integer(i4b) :: albedo_idx, emissivity_idx
        real(sp), allocatable :: s_zodi_scan(:)
        real(dp) :: chisq_tod, chisq_current, chisq_diff, ln_acceptance_probability, accept_rate, chisq_lnL
        real(dp), allocatable :: param_vec(:), powell_vec(:)
        type(ZodiModel) :: current_model, previous_model
        logical(lgt) :: accepted, verbose_, tuning, skip
        type(hdf_file) :: tod_file
        type(ZodiScan) :: zodi_timestreams

        if (present(verbose)) then
            verbose_ = verbose
        else
            verbose_ = .false.
        end if


        if (verbose_ .and. (cpar%myid == cpar%root)) then
            print *, "==== Zodi emissivity and albedo sampling step ===="
            print *, " "
        end if 

        ! Metropolis-Hastings for nproposals of new sets of zodi parameters
        nprop = 50
        ! nprop = cpar%zs_nprop

        ! Make two copies of the zodi model, one for the current proposal and one for the previous proposal

        if (gibbs_iter >= 25) then
            tuning = .false.
        else
            tuning = .true.
        end if
        tuning = .false.
        current_model = zodi_model_samp
        previous_model = zodi_model_samp
        allocate(param_vec(current_model%n_parameters))
        do i = 1, numband

            if (data(i)%tod_type == "none") cycle
            if (.not. data(i)%tod%sample_zodi) cycle

            ndet = data(i)%tod%ndet
            nscan = data(i)%tod%nscan
            allocate(zodi_timestreams%scans(nscan))

            ! Make sure that the zodi cache is cleared before each new band
            call data(i)%tod%clear_zodi_cache()

            ! Loop over over scans to generate zodi emission (just needed once per band)
            do scan = 1, nscan
                allocate(zodi_timestreams%scans(scan)%d(ndet))

                do j = 1, ndet
                    ! Skip scan if no accepted data
                    if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle

                    ntod = size(data(i)%tod%scans(scan)%d(j)%downsamp_res)
                    allocate(zodi_timestreams%scans(scan)%d(j)%s_scat(ntod, zodi_model%n_comps))
                    allocate(zodi_timestreams%scans(scan)%d(j)%s_therm(ntod, zodi_model%n_comps))

                    call get_zodi_emission(&
                        & tod=data(i)%tod, &
                        & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pointing, &
                        & scan=scan, &
                        & det=j, &
                        & s_zodi_scat=zodi_timestreams%scans(scan)%d(j)%s_scat, &
                        & s_zodi_therm=zodi_timestreams%scans(scan)%d(j)%s_therm, &
                        & model=current_model, &
                        & always_scattering=.false. &
                    &)
                end do
            end do

            n_accepted = 0
            do k = 0, nprop
                
                chisq_tod = 0.
                chisq_current = 0.

                ! if propsed parameter is outside of priors we dont want to immediately reject
                skip = .false.
                ! Root process draws new set of zodi parameters and broadcasts to all processes

                param_vec = current_model%model_to_param_vec()
                if (cpar%myid == cpar%root) then
                    ! Update all comp emissivites and albedos jointly (linear parameters)
                    param_loop: do j = 1, current_model%n_comps
                        ! Emissivities
                        emissivity_idx = get_emissivity_idx(j, data(i)%tod%band, current_model%n_bands, current_model%n_comps)
                        if (.not. active_params(emissivity_idx)) cycle

                        if (param_vec(emissivity_idx) /= 0.) then
                            param_vec(emissivity_idx) = param_vec(emissivity_idx) + (step_sizes(emissivity_idx) * rand_gauss(handle)) 
                            if (priors(emissivity_idx, 1) == NO_PRIOR .or. priors(emissivity_idx, 2) == NO_PRIOR) then 
                                continue
                            else if (param_vec(emissivity_idx) < priors(emissivity_idx, 1) .or. param_vec(emissivity_idx) > priors(emissivity_idx, 2)) then 
                                if (tuning) step_sizes(emissivity_idx) = step_sizes(emissivity_idx) / 2.
                                chisq_tod = 1.d30
                                skip = .true.
                                exit param_loop
                            end if
                        end if

                        ! Albedos
                        albedo_idx = get_albedo_idx(j, data(i)%tod%band, current_model%n_bands, current_model%n_comps)
                        if (.not. active_params(albedo_idx)) cycle

                        if (param_vec(albedo_idx) /= 0.) then
                            param_vec(albedo_idx) = param_vec(albedo_idx) + (step_sizes(albedo_idx) * rand_gauss(handle)) 

                            if (priors(albedo_idx, 1) == NO_PRIOR .or. priors(albedo_idx, 2) == NO_PRIOR) then 
                                continue
                            else if (param_vec(albedo_idx) < priors(albedo_idx, 1) .or. param_vec(albedo_idx) > priors(albedo_idx, 2)) then 
                                if (tuning) step_sizes(albedo_idx) = step_sizes(albedo_idx) / 2.
                                chisq_tod = 1.d30
                                skip = .true.
                                exit param_loop
                            end if
                        end if
                    end do param_loop
                end if

                if (tuning) call mpi_bcast(step_sizes, size(step_sizes), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
                call mpi_bcast(param_vec, size(param_vec), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
                call current_model%param_vec_to_model(param_vec)

                if (.not. skip) then
                    ! Update current_model with new parameters
                    scan_loop: do scan = 1, nscan
                        do j = 1, ndet
                            if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle

                            ntod = size(data(i)%tod%scans(scan)%d(j)%downsamp_res)
                            allocate(s_zodi_scan(ntod))
                            ! Evaluate zodi model with newly proposed values for each band and calculate chisq
                            call get_s_zodi(&
                                & s_therm=zodi_timestreams%scans(scan)%d(j)%s_therm, &
                                & s_scat=zodi_timestreams%scans(scan)%d(j)%s_scat, &
                                & s_zodi=s_zodi_scan, &
                                & band=data(i)%tod%band, &
                                & model=current_model &
                            &)
                            chisq_tod = chisq_tod + sum(((data(i)%tod%scans(scan)%d(j)%downsamp_res - s_zodi_scan)/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2)
                            deallocate(s_zodi_scan)

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
                        if (verbose_) then
                            print *, "proposal:", k
                            print *, "chisq_current:", chisq_current
                            print *, "chisq_diff:", chisq_diff
                            print *, "ln a:", ln_acceptance_probability 
                            print *, "accepted:", accepted
                            print *, "accept rate:", (real(n_accepted) / real(k))*100.
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
                    accept_rate = real(n_accepted) / real(k)
                end if
            end do  
            if (accept_rate < 0.1) then
                step_sizes(emissivity_start:albedo_stop) = step_sizes(emissivity_start:albedo_stop) / 2.
                if (verbose_ .and. (cpar%myid == cpar%root)) print *, "lowering stepsize", step_sizes(emissivity_start)
            else if (accept_rate > 0.6) then 
                step_sizes(emissivity_start:albedo_stop) = step_sizes(emissivity_start:albedo_stop) * 2.
                if (verbose_ .and. (cpar%myid == cpar%root)) print *, "increasing stepsize", step_sizes(emissivity_start)
            end if
            deallocate(zodi_timestreams%scans)
        end do
        zodi_model_samp = current_model
    end subroutine
       
end module
