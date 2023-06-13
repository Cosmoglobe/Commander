module comm_zodi_samp_mod
    use comm_tod_mod
    use comm_utils
    use comm_data_mod
    use comm_tod_zodi_mod
    use comm_zodi_mod
    implicit none

    private
    public sample_zodi_model

contains
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
        n_comps_to_fit = zodi%n_comps
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
                allocate(s_scat(ntod, zodi%n_comps), s_therm(ntod, zodi%n_comps))
                call get_zodi_emission(tod, tod%scans(scan)%d(j)%downsamp_pointing(1:padded_ubound), scan, j, s_scat, s_therm)
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

    subroutine update_xyz(handle)
        type(planck_rng), intent(inout) :: handle
        integer(i4b) :: i
        class(ZodiComponent), pointer :: comp

        comp => comp_list
        do while (associated(comp))
            print *, comp%x_0
            comp => comp%next()
        end do    
        stop
    end subroutine 


    subroutine sample_zodi_model(handle)
        ! Sample zodi model parameters
        type(planck_rng), intent(inout) :: handle
        integer(i4b) :: i, j, k, ndet, nscan, ntod, nprop, scan, ierr
        real(sp), allocatable :: s_scat(:, :), s_therm(:, :), s_zodi(:)
        real(dp) :: chisq, chisq_prev, reduced_chisq

        nprop = 1
        
        ! Metropolis-Hastings for nproposals of new sets of zodi parameters
        do k = 1, nprop
            chisq = 0.

            ! TODO: Loop over categories of shape parameters and propose new values (to all components)
            if (data(i)%tod%myid) call update_xyz(handle)
            do i = 1, numband
                ! Skip bands where we dont want to sample zodi
                if (.not. data(i)%tod%sample_zodi) cycle
                if (.not. allocated(data(i)%tod%scans(1)%d(1)%downsamp_res)) stop "cannot sample zodi parameters because downsamp_res and downsamp_pointing isnt allocated in the bands respective tod module. "

                ndet = data(i)%tod%ndet
                nscan = data(i)%tod%nscan

                ! Evaluate zodi model with newly proposed values for each band and calculate chisq
                do scan = 1, nscan
                    do j = 1, ndet
                        ! Allocate arrays
                        ntod = size(data(i)%tod%scans(scan)%d(j)%downsamp_res)
                        allocate(s_scat(ntod, zodi%n_comps), s_therm(ntod, zodi%n_comps), s_zodi(ntod))
                        call get_zodi_emission(&
                            & tod=data(i)%tod, &
                            & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pointing, &
                            & scan=scan, &
                            & det=j, &
                            & s_zodi_scat=s_scat, &
                            & s_zodi_therm=s_therm &
                        &)
                        s_zodi = 0.
                        call get_s_zodi(&
                            & emissivity=data(i)%tod%zodi_emissivity, &
                            & albedo=data(i)%tod%zodi_albedo, &
                            & s_therm=s_therm, &
                            & s_scat=s_scat, &
                            & s_zodi=s_zodi &
                        &)
                        chisq = chisq + sum((data(i)%tod%scans(scan)%d(j)%downsamp_res - s_zodi)**2)
                    end do
                    deallocate(s_scat, s_therm, s_zodi)
                end do

                ! Reduce chisq to root process
                call mpi_reduce(chisq, reduced_chisq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

                ! TODO: Accept or reject new parameter


                if (data(i)%tod%myid == 0) then
                    ! Sample new zodi parameters
                    print *, "chisq", reduced_chisq
                    stop
                end if

                ! Deallocate downsampled data and pointing
            end do
            chisq_prev = chisq
        end do  
    end subroutine
end module
