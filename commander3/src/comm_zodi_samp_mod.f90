module comm_zodi_samp_mod
    use comm_tod_mod
    use comm_utils
    use comm_data_mod
    use comm_tod_zodi_mod
    use comm_zodi_mod
    use comm_output_mod
    implicit none

    private
    public sample_zodi_model, output_zodi_model_to_hdf, initialize_zodi_samp_mod

    ! globals
    real(dp), allocatable :: chisq_iter(:), param_vec_iter(:, :)

contains
    subroutine initialize_zodi_samp_mod(cpar)
        ! Initialize the zodi sampling module.
        !
        ! Parameters
        ! ----------
        ! cpar: comm_params
        !    Parameter file variables.

        type(comm_params), intent(in) :: cpar
        allocate(chisq_iter(cpar%zs_nprop))
        allocate(param_vec_iter(cpar%zs_nprop, sampled_zodi_model%N_PARAMETERS))
    end subroutine initialize_zodi_samp_mod

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

    subroutine sample_linear_zodi_params(tod, handle, group_comps)
        class(comm_tod), intent(inout) :: tod
        type(planck_rng), intent(inout) :: handle
        logical(lgt), intent(in) :: group_comps
        real(sp), allocatable, dimension(:, :) :: s_scat, s_therm
        integer(i4b) :: i, j, ierr, scan, ndet, nhorn, nscan, ntod, n_comps_to_fit, padded_ubound
        real(dp), allocatable, dimension(:, :) :: A_T_A_emiss, A_T_A_albedo, A_T_A_emiss_reduced, A_T_A_albedo_reduced
        real(dp), allocatable, dimension(:) :: AY_emiss, AY_albedo, AY_emiss_reduced, AY_albedo_reduced, emissivities, albedos
        real(dp) :: n_used

        ndet = tod%ndet
        nscan = tod%nscan

        ! Get the number of components to fit (depends on if we want to sample the k98 groups or all components individually)
        ! and initialize Ax = Y matrices
        n_comps_to_fit = base_zodi_model%n_comps
        if (group_comps) then
            n_comps_to_fit = 3
        end if
        allocate(A_T_A_emiss(n_comps_to_fit, n_comps_to_fit))
        allocate(A_T_A_emiss_reduced(n_comps_to_fit, n_comps_to_fit))
        allocate(AY_emiss(n_comps_to_fit))
        allocate(AY_emiss_reduced(n_comps_to_fit))
        allocate(emissivities(n_comps_to_fit))
        A_T_A_emiss = 0.
        A_T_A_emiss_reduced = 0.
        AY_emiss = 0.
        AY_emiss_reduced = 0.

        ! Loop over downsampled data, and evaluate emissivity
        do scan = 1, nscan
            do j = 1, tod%ndet
                ntod = size(tod%scans(scan)%d(j)%downsamp_res)
                padded_ubound = ubound(tod%scans(scan)%d(j)%downsamp_res, dim=1) - 5 - 1 ! - 5 for padding
                allocate(s_scat(ntod, base_zodi_model%n_comps), s_therm(ntod, base_zodi_model%n_comps))
                call get_zodi_emission(tod, tod%scans(scan)%d(j)%downsamp_pointing(1:padded_ubound), scan, j, s_scat, s_therm, base_zodi_model)
                call accumulate_zodi_emissivities(tod, real(s_therm, dp), real(s_scat, dp), real(tod%scans(scan)%d(j)%downsamp_res(1:padded_ubound), dp), A_T_A_emiss, AY_emiss, group_comps)
                deallocate(s_scat, s_therm)
            end do
        end do
        call mpi_reduce(A_T_A_emiss, A_T_A_emiss_reduced, size(A_T_A_emiss), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call mpi_reduce(AY_emiss, AY_emiss_reduced, size(AY_emiss), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (tod%myid == 0) then
            call solve_Ax_zodi(A_T_A_emiss_reduced, AY_emiss_reduced, handle, emissivities)
            if (group_comps) then
               tod%zodi_emissivity(1) = emissivities(1)
               tod%zodi_emissivity(2:4) = emissivities(2)
               tod%zodi_emissivity(5:6) = emissivities(3)
            else 
               tod%zodi_emissivity = emissivities
            end if
            print *, "Sampled emissivity: ", emissivities
        end if
        call mpi_bcast(tod%zodi_emissivity, size(tod%zodi_emissivity), MPI_DOUBLE_PRECISION, 0, tod%comm, ierr)

        ! ! ! Loop over downsampled data, and evaluate albedo
        ! allocate(A_T_A_albedo(n_comps_to_fit, n_comps_to_fit))
        ! allocate(A_T_A_albedo_reduced(n_comps_to_fit, n_comps_to_fit))
        ! allocate(AY_albedo(n_comps_to_fit))
        ! allocate(AY_albedo_reduced(n_comps_to_fit))
        ! A_T_A_albedo = 0.
        ! A_T_A_albedo_reduced = 0.
        ! AY_albedo = 0.
        ! AY_albedo_reduced = 0.
        ! allocate(albedos(n_comps_to_fit))

        ! do scan = 1, nscan
        !     do j = 1, tod%ndet
        !         ntod = size(tod%scans(scan)%d(j)%downsamp_res)
        !         padded_ubound = ubound(tod%scans(scan)%d(j)%downsamp_res, dim=1) - 5 - 1 ! - 5 for padding
        !         allocate(s_scat(ntod, zodi%n_comps), s_therm(ntod, zodi%n_comps))
        !         call get_zodi_emission(tod, tod%scans(scan)%d(j)%downsamp_pointing(1:padded_ubound), scan, j, s_scat, s_therm)
        !         call accumulate_zodi_albedos(tod, real(s_therm, dp), real(s_scat, dp), real(tod%scans(scan)%d(j)%downsamp_res(1:padded_ubound), dp), A_T_A_albedo, AY_albedo, group_comps)
        !         deallocate(s_scat, s_therm)
        !     end do
        ! end do
        ! call mpi_reduce(A_T_A_albedo, A_T_A_albedo_reduced, size(A_T_A_albedo), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        ! call mpi_reduce(AY_albedo, AY_albedo_reduced, size(AY_albedo), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        ! if (tod%myid == 0) then
        !     call solve_Ax_zodi(A_T_A_albedo_reduced, AY_albedo_reduced, handle, albedos)
        !     if (group_comps) then
        !        tod%zodi_albedo(1) = albedos(1)
        !        tod%zodi_albedo(2:4) = albedos(2)
        !        tod%zodi_albedo(5:6) = albedos(3)
        !     else 
        !        tod%zodi_albedo = albedos
        !     end if
        !     print *, "Sampled albedo: ", albedos
        ! end if
        ! call mpi_bcast(tod%zodi_albedo, size(tod%zodi_albedo), MPI_DOUBLE_PRECISION, 0, tod%comm, ierr)
    end subroutine

    subroutine draw_n0(model, cpar, handle)
        type(zodi_model), intent(inout) :: model
        type(comm_params), intent(in) :: cpar
        type(planck_rng), intent(inout) :: handle
        integer(i4b) :: i, ierr
        real(dp) :: n0_std(6), i_std(6), c(6), eta

        n0_std = [6.4d-10, 7.2d-11, 1.28d-10, 2.34d-11, 1.27d-9, 1.42d-9]
        i_std = [6.4d-10, 7.2d-11, 1.28d-10, 2.34d-11, 1.27d-9, 1.42d-9]
        c = [1.d-2, 1.d-2, 1.d22, 1.d-2, 1.d-2, 1.d-2]

        ! do i = 1, 1
        do i = 1, model%n_comps
            if (cpar%myid == cpar%root) then
                ! n0
                eta = rand_gauss(handle)
                model%comps(i)%c%n_0 = model%comps(i)%c%n_0 + (n0_std(i) * eta)

                ! incl
                ! eta = rand_gauss(handle)
                ! model%comps(i)%c%i = model%comps(i)%c%i + (i_std(i) * eta)
            end if
            call mpi_bcast(model%comps(i)%c%n_0, 1, MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
        end do
    end subroutine 


    subroutine sample_zodi_model(cpar, handle)
        ! Metropolis-Hastings for nproposals of new sets of zodi parameters
        ! Todo: Split into gibbs steps for each component
        type(comm_params), intent(in) :: cpar
        type(planck_rng), intent(inout) :: handle
        integer(i4b) :: i, j, k, ndet, nscan, ntod, nprop, scan, ierr, n_accepted, n_tot_tod, n_tot_tod_reduced, n_tot_tod_current
        real(sp), allocatable :: s_scat(:, :), s_therm(:, :), s_zodi(:)
        real(dp) :: chisq_tod, chisq_current, chisq_previous, chisq_diff, ln_acceptance_probability
        type(zodi_model) :: current_model, previous_model
        logical(lgt) :: accepted
        ! real(dp) :: n0_grid(100), chisq_grid(100)

        !debugging stuff
        character(len=256) :: chaindir
        character(len=512) :: path, param_path, param_name
        character(len=6) :: itext
        real(dp) :: chisq_buf
        type(hdf_file) :: tod_file
        integer(i4b) :: l


        nprop = 25

        ! Metropolis-Hastings for nproposals of new sets of zodi parameters
        chisq_previous = 1d20
        n_accepted = 0

        current_model = sampled_zodi_model
        previous_model = sampled_zodi_model

        do k = 0, nprop
            ! Reset chisq for current proposal
            chisq_tod = 0.
            chisq_current = 0.

            ! Root process draws new set of zodi parameters and broadcasts to all processes
            ! If first iteration we dont want to draw, just compute the chisq with the base model
            if (k > 0) call draw_n0(current_model, cpar, handle)
   
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

                        ! Add noise to simulations
                        ! do l = 1, ntod
                        !     data(i)%tod%scans(scan)%d(j)%downsamp_res(l) = data(i)%tod%scans(scan)%d(j)%downsamp_res(l) + rand_gauss(handle) * data(i)%tod%scans(scan)%d(j)%N_psd%sigma0
                        ! end do

                        allocate(s_scat(ntod, base_zodi_model%n_comps), s_therm(ntod, base_zodi_model%n_comps), s_zodi(ntod))
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
                        
                        if (.false. .and. cpar%myid == cpar%root .and. scan == 1) then
                            chaindir = "/mn/stornext/d5/data/metins/dirbe/chains/chains_testing"
                            call open_hdf_file(trim(chaindir)//'/timestream.h5', tod_file, 'w')
                            call int2string(k, itext)
                            path = trim(adjustl(itext))//'/'
                            call create_hdf_group(tod_file, trim(adjustl(path)))
                            param_name = "res"
                            param_path = trim(adjustl(path))//trim(adjustl(param_name))
                            call write_hdf(tod_file, trim(adjustl(param_path)), data(i)%tod%scans(scan)%d(j)%downsamp_res)
                            param_name = "pix"
                            param_path = trim(adjustl(path))//trim(adjustl(param_name))
                            call write_hdf(tod_file, trim(adjustl(param_path)), data(i)%tod%scans(scan)%d(j)%downsamp_pointing)
                            param_name = "zodi"
                            param_path = trim(adjustl(path))//trim(adjustl(param_name))
                            call write_hdf(tod_file, trim(adjustl(param_path)), s_zodi)
                            param_name = "n0"
                            param_path = trim(adjustl(path))//trim(adjustl(param_name))
                            call write_hdf(tod_file, trim(adjustl(param_path)), current_model%comps(1)%c%n_0)
                            if (k > 0) then
                                call close_hdf_file(tod_file)
                                stop                        
                            endif
                        end if 

                        ! if (cpar%myid == 0) print *, "chisq_before: ", chisq_tod
                        chisq_buf = sum(((data(i)%tod%scans(scan)%d(j)%downsamp_res - s_zodi)/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2)
                        ! if (cpar%myid == 0) print *, "chisq_buf: ", chisq_buf
                        chisq_tod = chisq_tod + chisq_buf
                        ! if (cpar%myid == 0) print *, "chisq_after: ", chisq_tod

                        ! Normalize chisq_tod
                        ! chisq_tod = chisq_tod - ntod/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0**2
                        deallocate(s_scat, s_therm, s_zodi)
                    end do
                end do
            end do

            ! Reduce chisq to root process
            ! call mpi_reduce(n_tot_tod, n_tot_tod_reduced, 1, MPI_INTEGER, MPI_SUM, cpar%root, MPI_COMM_WORLD, ierr)
            call mpi_reduce(chisq_tod, chisq_current, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cpar%root, MPI_COMM_WORLD, ierr)

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
                chisq_iter(k) = chisq_current
                param_vec_iter(k, :) = current_model%model_to_param_vec()
            end if
        end do  
        sampled_zodi_model = current_model
        ! if (cpar%myid == cpar%root) then
        !     chaindir = "/mn/stornext/d5/data/metins/dirbe/chains/chains_testing"
        !     call open_hdf_file(trim(chaindir)//'/chisq_r.h5', tod_file, 'w')
        !     call write_hdf(tod_file, '/chisq', chisq_grid)
        !     call close_hdf_file(tod_file)
        !     stop
        ! end if 
        

    end subroutine

    subroutine output_zodi_model_to_hdf(cpar, iter, model)
        ! Writes the zodi model to an hdf file
        type(comm_params), intent(in) :: cpar
        integer(i4b), intent(in) :: iter
        type(zodi_model), intent(in) :: model


        integer(i4b) :: i, j, hdferr, ierr, unit
        logical(lgt) :: exist, init, new_header
        character(len=6) :: itext
        character(len=4) :: ctext
        character(len=512) :: zodi_path, comp_path, param_path, chainfile, hdfpath
        character(len=10), allocatable :: param_names(:)
        real(dp), allocatable :: param_values(:)
        type(hdf_file) :: file
        type(h5o_info_t) :: object_info

        if (.not. cpar%myid_chain == 0) return

        call int2string(cpar%mychain, ctext)
        call int2string(iter,         itext)
        chainfile = trim(adjustl(cpar%outdir)) // '/chain' // &
            & '_c' // trim(adjustl(ctext)) // '.h5'

        inquire(file=trim(chainfile), exist=exist)
        call open_hdf_file(chainfile, file, 'b')

        call int2string(iter, itext)
        zodi_path = trim(adjustl(itext))//'/zodi'
        call create_hdf_group(file, trim(adjustl(zodi_path)))
        comp_path = trim(adjustl(zodi_path))//'/params/'
        call create_hdf_group(file, trim(adjustl(comp_path)))

        call write_hdf(file, trim(adjustl(comp_path)), param_vec_iter)

        ! Write chisq
        comp_path = trim(adjustl(zodi_path))//'/chisq/'
        call create_hdf_group(file, trim(adjustl(comp_path)))
        call write_hdf(file, trim(adjustl(comp_path))//'/chisq', chisq_iter)

        call close_hdf_file(file)
    end subroutine



end module
