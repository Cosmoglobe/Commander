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
    public initialize_zodi_samp_mod, sample_zodi_model, init_scandata_and_downsamp_zodi, gibbs_sample_zodi_comp

    ! globals
    real(dp), allocatable :: chisq_previous, step_size, step_sizes(:), priors(:, :)
    logical(lgt), allocatable :: active_params(:)
    integer(i4b) :: n_active_params, cloud_start, cloud_stop, band1_start, band1_stop, band2_start, band2_stop, band3_start, band3_stop, ring_start, ring_stop, feature_start, feature_stop
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
        allocate(step_sizes(zodi_model%n_parameters))
        allocate(active_params(zodi_model%n_parameters))
        allocate(param_vec(zodi_model%n_parameters))
        param_vec = zodi_model%model_to_param_vec()
        step_sizes = 0.1
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
        active_params(61) = .false. !T_0
        active_params(62) = .false. !delta

        ! emissivities
        active_params(63 : 63 + n_spec_params) = .false.

        ! albedos
        active_params(63 + n_spec_params : 63 + 2 * n_spec_params) = .false.

        n_active_params = count(active_params == .true.)

        ! Priors for zodi parameters
        ! Cloud
        priors(1, :) = [1.0d-7, 1.4d-7] !n_0
        priors(2, :) = [0., 4.] !incl
        priors(3, :) = [67., 87.] !omega
        priors(4, :) = [-0.05, 0.05] !x_0
        priors(5, :) = [-0.05, 0.05] !y_0
        priors(6, :) = [-0.05, 0.05] !z_0
        priors(7, :) = [1., 2.] !alpha
        priors(8, :) = [3., 5.] !beta
        priors(9, :) = [0.7, 1.2] !gamma
        priors(10, :) = [0.0, 0.5] !mu

        ! Band1
        priors(11, :) = [5.59d-11, 5.59d-9] !n_0
        priors(12, :) = [0., 2.] !incl
        priors(13, :) = [0., 180.] !omega
        priors(14, :) = [-0.1, 0.1] !x_0
        priors(15, :) = [-0.1, 0.1] !y_0
        priors(16, :) = [-0.1, 0.1] !z_0
        priors(17, :) = [6.78, 10.78] !delta_zeta
        priors(18, :) = [1., 3.] !delta_r
        priors(19, :) = [0., 1.] !v
        priors(20, :) = [2., 6.] !p

        ! Band2
        priors(21, :) = [1.99d-10, 1.99d-8] !n_0
        priors(22, :) = [0., 10.] !incl
        priors(23, :) = [10., 50.] !omega
        priors(24, :) = [-0.1, 0.1] !x_0
        priors(25, :) = [-0.1, 0.1] !y_0
        priors(26, :) = [-0.1, 0.1] !z_0
        priors(27, :) = [0., 5.] !delta_zeta
        priors(28, :) = [0.6, 3.] !delta_r
        priors(29, :) = [0., 3.] !v
        priors(30, :) = [2., 6.] !p

        ! Band3
        priors(31, :) = [1.44d-11, 1.44d-9] !n_0
        priors(32, :) = [0., 5.] !incl
        priors(33, :) = [40., 120.] !omega
        priors(34, :) = [-0.1, 0.1] !x_0
        priors(35, :) = [-0.1, 0.1] !y_0
        priors(36, :) = [-0.1, 0.1] !z_0
        priors(37, :) = [10., 20.] !delta_zeta
        priors(38, :) = [1., 4.] !delta_r
        priors(39, :) = [0., 1.] !v
        priors(40, :) = [2., 6.] !p

        ! Ring
        priors(41, :) = [1.83d-9, 1.83d-7] !n_0
        priors(42, :) = [0., 2.] !incl
        priors(43, :) = [0., 40.] !omega
        priors(44, :) = [-0.1, 0.1] !x_0
        priors(45, :) = [-0.1, 0.1] !y_0
        priors(46, :) = [-0.1, 0.1] !z_0
        priors(47, :) = [0.9, 1.1] !R_0
        priors(48, :) = [0., 0.05] !sigma_r
        priors(49, :) = [0., 0.1] !sigma_z

        ! Feature
        priors(50, :) = [1.9d-9, 1.9d-7] !n_0
        priors(51, :) = [0., 2.] !incl
        priors(52, :) = [0., 40.] !omega
        priors(53, :) = [-0.1, 0.1] !x_0
        priors(54, :) = [-0.1, 0.1] !y_0
        priors(55, :) = [-0.1, 0.1] !z_0
        priors(56, :) = [0.85, 1.15] !R_0
        priors(57, :) = [0., 0.2] !sigma_r
        priors(58, :) = [0., 0.2] !sigma_z
        priors(59, :) = [-20., 5.] !theta_0
        priors(60, :) = [0., 20.] !sigma_theta

        ! Other
        priors(61, :) = [250., 315.] !T_0
        priors(62, :) = [0.2, 0.7] !delta

        ! Emissivity
        priors(63 : 63 + n_spec_params, 0) = 0.
        priors(63 : 63 + n_spec_params, 1) = 5.

        ! Albedo
        priors(63 + n_spec_params : 63 + 2 * n_spec_params, 0) = 0.
        priors(63 + n_spec_params : 63 + 2 * n_spec_params, 1) = 1.
    end subroutine initialize_zodi_samp_mod

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
                            print *, "before", param_vec(i)
                            param_vec(i) = param_vec(i) + (step_sizes(i) * param_vec(i) * rand_gauss(handle))
                            print *, "after", param_vec(i)
                            if (param_vec(i) < priors(i, 1) .or. param_vec(i) > priors(i, 2)) then 
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
                    if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
                    do j = 1, ndet
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
                            & emissivity=data(i)%tod%zodi_emissivity, &
                            & albedo=data(i)%tod%zodi_albedo, &
                            & s_therm=s_therm, &
                            & s_scat=s_scat, &
                            & s_zodi=s_zodi &
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
        integer(i4b), intent(in) :: comp_idx

        integer(i4b) :: i, j , k, nside, npix, nmaps, ndelta
        real(sp), allocatable, dimension(:, :, :, :) :: map_sky
        real(sp), allocatable, dimension(:) :: procmask, procmask2, procmask_zodi
        real(dp), allocatable, dimension(:, :) :: m_buf
        integer(i4b), allocatable, dimension(:) :: downsamp_width
        type(comm_scandata) :: sd
        type(map_ptr), allocatable, dimension(:, :) :: s_sky

        allocate(downsamp_width(zodi_model%n_comps))

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
            downsamp_width = data(i)%tod%samprate / data(i)%tod%samprate_lowres

            do j = 1, data(i)%tod%nscan
                if (.not. any(data(i)%tod%scans(j)%d%accept)) cycle

                ! not calibrating the data so we use substitute map_gain with map_sky
                call sd%init_singlehorn(data(i)%tod, j, map_sky, map_sky, procmask, procmask2, procmask_zodi)
                call boxcar_downsamp_zodi_res_and_pointing(data(i)%tod, sd, j, padding=5, box_width=downsamp_width(comp_idx))

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
        if (mod(box_halfwidth, 2) == 0) then 
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

    ! function for Powell()
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
                if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
                do j = 1, ndet
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
                        & emissivity=data(i)%tod%zodi_emissivity, &
                        & albedo=data(i)%tod%zodi_albedo, &
                        & s_therm=s_therm, &
                        & s_scat=s_scat, &
                        & s_zodi=s_zodi &
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




    subroutine gibbs_sample_zodi_comp(cpar, handle, comp_label, verbose)
        type(comm_params), intent(in) :: cpar
        type(planck_rng), intent(inout) :: handle
        character(len=*), intent(in) :: comp_label
        logical(lgt), intent(in), optional :: verbose
        integer(i4b) :: i, j, k, param_start, param_stop, flag, ndet, nscan, ntod, nprop, scan, ierr, n_accepted, n_tot_tod, n_tot_tod_reduced, n_tot_tod_current
        real(sp), allocatable :: s_scat(:, :), s_therm(:, :), s_zodi(:)
        real(dp) :: chisq_tod, chisq_current, chisq_diff, ln_acceptance_probability, accept_rate, chisq_lnL
        real(dp), allocatable :: param_vec(:), powell_vec(:)
        type(ZodiModel) :: current_model, previous_model
        logical(lgt) :: accepted, verbose_

        if (present(verbose)) then
            verbose_ = verbose
        else
            verbose_ = .false.
        end if

        if (verbose_ .and. (cpar%myid == cpar%root)) then
            print *, "==== Zodi gibbs sampling step: ", trim(comp_label), " ===="
            print *, " "
        end if 

        ! Get parameter indicies for this component
        select case (trim(comp_label)) 
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
        end select

        ! Metropolis-Hastings for nproposals of new sets of zodi parameters
        nprop = cpar%zs_nprop

        ! Make two copies of the zodi model, one for the current proposal and one for the previous proposal
        current_model = zodi_model
        previous_model = zodi_model

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
                        param_vec(i) = param_vec(i) + (step_sizes(i) * param_vec(i) * rand_gauss(handle))

                        ! If parameter is outside of prior range, immediatly reject proposal and halve step size
                        if (param_vec(i) < priors(i, 1) .or. param_vec(i) > priors(i, 2)) then 
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
                    if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
                    do j = 1, ndet
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
                            & emissivity=data(i)%tod%zodi_emissivity, &
                            & albedo=data(i)%tod%zodi_albedo, &
                            & s_therm=s_therm, &
                            & s_scat=s_scat, &
                            & s_zodi=s_zodi &
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
                accept_rate = n_accepted / k
            end if
        end do  
        if (accept_rate < 10.) then
            step_sizes(param_start:param_stop) = step_sizes(param_start:param_stop) / 2.
            if (verbose_ .and. (cpar%myid == cpar%root)) print *, "lowering stepsize"
        else if (accept_rate > 60.) then 
            step_sizes(param_start:param_stop) = step_sizes(param_start:param_stop) * 2.
            if (verbose_ .and. (cpar%myid == cpar%root)) print *, "increasing stepsize"
        end if
        
        zodi_model = current_model
    end subroutine

end module
