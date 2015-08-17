module comm_mp_mod
  use healpix_types
  use comm_genvec_mod
  use comm_utils
  use comm_cgd_matmul_mod
  use comm_cgd_precond_mod
  use comm_data_mod
  use comm_fg_mod
  use comm_chisq_mod
  !use comm_checkpoint_mod
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2006, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************

  integer(i4b),                             private :: comm_chain, myid_chain
  integer(i4b),                             private :: ierr, root = 0
  integer(i4b), dimension(MPI_STATUS_SIZE), private :: status
  type(planck_rng),                         private :: handle

contains

  subroutine initialize_mp_mod(comm_chain_in, handle_in)
    implicit none

    integer(i4b),     intent(in)    :: comm_chain_in
    type(planck_rng), intent(inout) :: handle_in
    
    integer(i4b) :: seed

    comm_chain = comm_chain_in
    call mpi_comm_rank(comm_chain, myid_chain, ierr)

    seed = nint(100000000*rand_uni(handle_in))
    call rand_init(handle, seed)

  end subroutine initialize_mp_mod


  subroutine slave_carry_out_operation(slave_work)
    implicit none

    logical(lgt), intent(inout) :: slave_work

    integer(i4b) :: operation, comp, par
    real(dp)     :: rms


    call mpi_bcast(operation, 1, MPI_INTEGER, root, comm_chain, ierr)

    if (operation == 1) then

       ! Compute a (A^T N^-1 A)*x product
       call compute_At_invN_A_x

    else if (operation == 2) then

       ! Set up the right-hand side of equations for the signal component
       call compute_signal_rhs
       
    else if (operation == 3) then

       ! Gain sampler
       call sample_gain(handle, map_id)

    else if (operation == 4) then
       
       ! Empty slot

    else if (operation == 5) then
       
       ! Empty slot

    else if (operation == 6) then
       
       ! Initialize preconditioner
       call initialize_preconditioner

    else if (operation == 7) then

       ! Compute residuals
       call compute_residual

    else if (operation == 8) then

       ! Compute signal map
       call compute_signal_map

    else if (operation == 9) then

       ! Empty slot

    else if (operation == 10) then

       ! Empty slot

    else if (operation == 11) then

       ! Initialize preconditioner
       call initialize_preconditioner

    else if (operation == 12) then

       ! Sample spectral index map
       call update_fg_pix_response_map(map_id, pixels, fg_pix_spec_response)

    else if (operation == 13) then

       ! Sample spectral index map
       call sample_fg_pix_response_map

    else if (operation == 14) then

       ! Update par_rms
       call mpi_bcast(comp, 1, MPI_INTEGER,          root, comm_chain, ierr)
       call mpi_bcast(par,  1, MPI_INTEGER,          root, comm_chain, ierr)
       call mpi_bcast(rms,  1, MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
       fg_components(comp)%p_rms(par) = rms

    else if (operation == 15) then

       ! Compute chisquares
       call compute_total_chisq(map_id)

    else if (operation == 16) then

       ! Compute chisquares
       call set_mask_state

    else if (operation == 17) then

       ! Set checkpoint
       !call checkpoint_all_modules

    else if (operation == 18) then

       call set_noiseamp

    else if (operation == 21) then

!       call sample_tab_fg_spectra(map_id, residual, inv_N_covar_lowres(:,:,:,1), fg_pix_spec_response, pixels)

    else if (operation == 26) then

       call output_component_map(map_id)

    else if (operation == 27) then

       call update_eff_fg_spectrum

    else if (operation == 28) then

       call output_mixing_matrix(map_id)

    else if (operation == 29) then

       ! Sample spectral index map
!       call fit_all_fixed_mono_dipoles(map_id, cmbmap, cmbmaps_lowres, fg_temp_free, fg_temp_free_lowres)

    else if (operation == 30) then

       !call subtract_fg_mono_dipole

    else if (operation == 31) then

       call init_ind

    else if (operation == 32) then

       !call sample_global_fg_params

    else if (operation == 33) then

       call enforce_pos_amp

    else if (operation == 34) then

       call mpi_bcast(exclude_pos_fg_amp, 1, MPI_LOGICAL, root, comm_chain, ierr)

    else if (operation == 35) then

       call mpi_bcast(sample_temp_coeffs, 1, MPI_LOGICAL, root, comm_chain, ierr)

    else if (operation == 36) then

       !call est_const_ind

    else if (operation == 40) then

       ! Gain sampler
       call sample_bandpass(handle, map_id)

    else if (operation == 50) then

       call output_ml_map_engine(myid_chain)

    else if (operation == 0) then

       ! Quit
       slave_work = .false.

    end if

  end subroutine slave_carry_out_operation


  subroutine compute_rhs_product(add_wf, add_S_omega, add_N_omega, rhs, fluct)
    implicit none

    logical(lgt), intent(in)    :: add_wf, add_S_omega, add_N_omega
    type(genvec), intent(inout) :: rhs
    type(genvec), intent(inout), optional :: fluct

    call mpi_bcast(2, 1, MPI_INTEGER, root, comm_chain, ierr)

    call compute_signal_rhs(add_wf, add_S_omega, add_N_omega, rhs, fluct)

  end subroutine compute_rhs_product


  subroutine compute_residuals(s, subtract_signal)
    implicit none

    type(genvec), intent(in)      :: s
    logical(lgt), intent(in)      :: subtract_signal

    call mpi_bcast(7, 1, MPI_INTEGER, root, comm_chain, ierr)
    call compute_residual(s, subtract_signal)       

  end subroutine compute_residuals

  subroutine compute_signal_maps(coeff)
    implicit none

    type(genvec), intent(inout) :: coeff

    call mpi_bcast(8, 1, MPI_INTEGER, root, comm_chain, ierr)
    call compute_signal_map(coeff)

  end subroutine compute_signal_maps


  ! Chi-square calculation
  subroutine compute_chisq(map_id, coeff, output_stats, chisq_fullsky, chisq_highlat, chisq_map, &
       & chisq_rms, chisq_band, nu_band, chain, iter)
    implicit none

    integer(i4b),                   intent(in)               :: map_id
    type(genvec),                   intent(in)               :: coeff
    logical(lgt),                   intent(in),  optional    :: output_stats
    real(dp),                       intent(out), optional    :: chisq_highlat, chisq_fullsky
    real(dp),                       intent(out), optional    :: chisq_rms
    real(dp),     dimension(1:),    intent(out), optional    :: chisq_band, nu_band
    real(dp),     dimension(0:,1:), intent(out), optional    :: chisq_map
    integer(i4b),                   intent(in),  optional    :: chain, iter

    ! Multiply with beam and inverse noise matrix
    call mpi_bcast(15, 1, MPI_INTEGER, root, comm_chain, ierr)
    call compute_total_chisq(map_id, coeff, output_stats=output_stats, chisq_map=chisq_map, &
         & chisq_fullsky=chisq_fullsky, chisq_highlat=chisq_highlat, chisq_rms=chisq_rms, &
         & chisq_band=chisq_band, nu_band=nu_band, chain=chain, iter=iter)    

  end subroutine compute_chisq


  ! Matrix multiplication routines
  subroutine init_precond(precond_type)
    implicit none

    integer(i4b), intent(in) :: precond_type

    call mpi_bcast(11, 1, MPI_INTEGER, root, comm_chain, ierr)
    call initialize_preconditioner(precond_type)
       
  end subroutine init_precond

  ! Update response maps
  subroutine update_fg_pix_response_maps(index_map)
    implicit none

    real(dp), dimension(0:,1:,1:), intent(in) :: index_map

    call mpi_bcast(12, 1, MPI_INTEGER, root, comm_chain, ierr)
    call update_fg_pix_response_map(map_id, pixels, fg_pix_spec_response, index_map)
    
  end subroutine update_fg_pix_response_maps

  ! Sample foreground spectral indices
  subroutine sample_spectral_param_map(s, residuals, inv_N_noise, fg_amp, index_map, stat)
    implicit none

    type(genvec),                     intent(in)    :: s
    real(dp), dimension(0:,1:,1:),    intent(in)    :: residuals, fg_amp
    real(dp), dimension(0:,1:,1:),    intent(in)    :: inv_N_noise
    real(dp), dimension(0:,1:,1:),    intent(out)   :: index_map
    integer(i4b),                     intent(inout) :: stat

    call mpi_bcast(13, 1, MPI_INTEGER, root, comm_chain, ierr)
    call sample_fg_pix_response_map(s, residuals, inv_N_noise, fg_amp, index_map, stat)
    
  end subroutine sample_spectral_param_map

  ! Sample foreground spectral indices
  subroutine update_par_rms(comp, par, rms)
    implicit none

    integer(i4b), intent(in) :: comp, par
    real(dp),     intent(in) :: rms

    call mpi_bcast(14,   1, MPI_INTEGER,          root, comm_chain, ierr)
    call mpi_bcast(comp, 1, MPI_INTEGER,          root, comm_chain, ierr)
    call mpi_bcast(par,  1, MPI_INTEGER,          root, comm_chain, ierr)
    call mpi_bcast(rms,  1, MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
    fg_components(comp)%p_rms(par) = rms
    
  end subroutine update_par_rms

  subroutine update_all_eff_fg_spectrum(comp)
    implicit none

    integer(i4b), intent(in), optional :: comp

    call mpi_bcast(27,   1, MPI_INTEGER,          root, comm_chain, ierr)
    call update_eff_fg_spectrum(comp)

  end subroutine update_all_eff_fg_spectrum

  ! Update CG constraints
!!$  subroutine update_cgd_constraint_module
!!$    implicit none
!!$
!!$    call mpi_bcast(14, 1, MPI_INTEGER, root, comm_chain, ierr)
!!$    call update_cgd_constraints
!!$    
!!$  end subroutine update_cgd_constraint_module

  ! Set mask_state
  subroutine update_mask_state(mask_state)
    implicit none

    integer(i4b), intent(in) :: mask_state

    call mpi_bcast(16, 1, MPI_INTEGER, root, comm_chain, ierr)
    call set_mask_state(mask_state)
    
  end subroutine update_mask_state

  ! Checkpoint session
!!$  subroutine set_checkpoint(operation)
!!$    implicit none
!!$
!!$    character(len=128), intent(in) :: operation
!!$
!!$    call mpi_bcast(17, 1, MPI_INTEGER, root, comm_chain, ierr)
!!$    call checkpoint_all_modules(operation)
!!$    
!!$  end subroutine set_checkpoint

  ! Checkpoint session
  subroutine set_noiseamps(amps)
    implicit none

    real(dp), dimension(1:), intent(in) :: amps

    call mpi_bcast(18,    1, MPI_INTEGER,          root, comm_chain, ierr)
    call set_noiseamp(amps)
    
  end subroutine set_noiseamps

!  subroutine sample_tabulated_fg_spectra(fg_amp)
!    implicit none

!    real(dp), dimension(0:,1:,1:),    intent(in) :: fg_amp

!    call mpi_bcast(21,    1, MPI_INTEGER,          root, comm_chain, ierr)
!    call sample_tab_fg_spectra(map_id, residual, inv_N_covar_lowres(:,:,:,1), fg_pix_spec_response, pixels, fg_amp)

!  end subroutine sample_tabulated_fg_spectra

  subroutine output_component_maps(map_id, s, chain, iter, chain_dir)
    implicit none

    integer(i4b),     intent(in) :: map_id, chain, iter
    character(len=*), intent(in) :: chain_dir
    type(genvec),     intent(in) :: s

    call mpi_bcast(26,    1, MPI_INTEGER,          root, comm_chain, ierr)
    call output_component_map(map_id, s, chain, iter, chain_dir)

  end subroutine output_component_maps

  subroutine output_mixing_matrices(map_id, chain, iter, chain_dir)
    implicit none

    integer(i4b),     intent(in) :: map_id, chain, iter
    character(len=*), intent(in) :: chain_dir

    call mpi_bcast(28,    1, MPI_INTEGER,          root, comm_chain, ierr)
    call output_mixing_matrix(map_id, chain, iter, chain_dir)

  end subroutine output_mixing_matrices

  ! Sample foreground spectral indices
!!$  subroutine fit_fixed_mono_dipoles(chain_dir, chain, iter, residuals, inv_N_noise, fg_amp, index_map, md_out)
!!$    implicit none
!!$
!!$    character(len=*),                 intent(in)  :: chain_dir
!!$    integer(i4b),                     intent(in)  :: chain, iter
!!$    real(dp), dimension(0:,1:,1:),    intent(in)  :: residuals, fg_amp
!!$    real(dp), dimension(0:,1:,1:),    intent(in)  :: inv_N_noise
!!$    real(dp), dimension(0:,1:,1:),    intent(out) :: index_map
!!$    real(dp), dimension(1:,1:),       intent(out) :: md_out
!!$
!!$    call mpi_bcast(29, 1, MPI_INTEGER, root, comm_chain, ierr)
!!$!    call fit_all_fixed_mono_dipoles(map_id, cmbmap, cmbmaps_lowres, fg_temp_free, fg_temp_free_lowres, &
!!$!         & chain_dir, chain, iter, residuals, inv_N_noise, fg_amp, index_map, md_out)
!!$
!!$  end subroutine fit_fixed_mono_dipoles

!!$  subroutine sample_a2t_factors(s)
!!$    implicit none
!!$
!!$    type(genvec), intent(in) :: s
!!$
!!$    call mpi_bcast(27,    1, MPI_INTEGER,          root, comm_chain, ierr)
!!$    call sample_a2t(map_id, s)
!!$
!!$  end subroutine sample_a2t_factors

!!$  subroutine subtract_fg_mono_dipoles(md)
!!$    implicit none
!!$
!!$    real(dp), dimension(1:,1:), intent(in) :: md
!!$
!!$    call mpi_bcast(30,    1, MPI_INTEGER,          root, comm_chain, ierr)
!!$    !call subtract_fg_mono_dipole(md)
!!$
!!$  end subroutine subtract_fg_mono_dipoles

  subroutine init_ind_by_nonlin_search_co(residuals, fg_amp, inv_N_noise, index_map)
    implicit none

    real(dp), dimension(0:,1:,1:),    intent(in)    :: residuals, fg_amp, inv_N_noise
    real(dp), dimension(0:,1:,1:),    intent(inout) :: index_map

    call mpi_bcast(35,    1, MPI_INTEGER,          root, comm_chain, ierr)
!    call init_ind_co(residuals, inv_N_noise, index_map)

  end subroutine init_ind_by_nonlin_search_co

  subroutine init_ind_by_nonlin_search(residuals, inv_N_noise, index_map, fg_amp)
    implicit none

    real(dp), dimension(0:,1:,1:),    intent(in)    :: residuals
    real(dp), dimension(0:,1:,1:),    intent(in)    :: inv_N_noise
    real(dp), dimension(0:,1:,1:),    intent(inout) :: index_map
    real(dp), dimension(0:,1:,1:),    intent(inout) :: fg_amp

    call mpi_bcast(31,    1, MPI_INTEGER,          root, comm_chain, ierr)
    call init_ind(residuals, inv_N_noise, index_map, fg_amp)

  end subroutine init_ind_by_nonlin_search

!!$  subroutine init_ind_by_nonlin_search_ame(residuals, inv_N_noise, index_map)
!!$    implicit none
!!$
!!$    real(dp), dimension(0:,1:,1:),    intent(in)    :: residuals
!!$    real(dp), dimension(0:,1:,1:),    intent(in)    :: inv_N_noise
!!$    real(dp), dimension(0:,1:,1:),    intent(inout) :: index_map
!!$
!!$    call mpi_bcast(34,    1, MPI_INTEGER,          root, comm_chain, ierr)
!!$    !call init_ind_ame(residuals, inv_N_noise, index_map)
!!$
!!$  end subroutine init_ind_by_nonlin_search_ame

!!$  subroutine estimate_constant_indices(residuals, inv_N_noise, index_map)
!!$    implicit none
!!$
!!$    real(dp), dimension(0:,1:,1:),    intent(in)    :: residuals
!!$    real(dp), dimension(0:,1:,1:),    intent(in)    :: inv_N_noise
!!$    real(dp), dimension(0:,1:,1:),    intent(inout) :: index_map
!!$
!!$    call mpi_bcast(36,    1, MPI_INTEGER,          root, comm_chain, ierr)
!!$    !call est_const_ind(residuals, inv_N_noise, index_map)
!!$
!!$  end subroutine estimate_constant_indices

!!$  subroutine get_masked_region_alms(alms_in, alms_out)
!!$    implicit none
!!$
!!$    real(dp), dimension(1:,1:), intent(in)  :: alms_in
!!$    real(dp), dimension(1:,1:), intent(out) :: alms_out
!!$
!!$    call mpi_bcast(32,    1, MPI_INTEGER,          root, comm_chain, ierr)
!!$    call compute_masked_region_alms(alms_in, alms_out)
!!$
!!$  end subroutine get_masked_region_alms

  subroutine enforce_pos_amps(chaindir, residuals, inv_N_noise, fg_amp, index_map, burnin) 
    implicit none
    
    character(len=*),                 intent(in)    :: chaindir
    real(dp), dimension(0:,1:,1:),    intent(in)    :: residuals
    real(dp), dimension(0:,1:,1:),    intent(inout) :: fg_amp
    real(dp), dimension(0:,1:,1:),    intent(in)    :: inv_N_noise
    real(dp), dimension(0:,1:,1:),    intent(in)    :: index_map
    logical(lgt),                     intent(in)    :: burnin

    call mpi_bcast(33,    1, MPI_INTEGER,          root, comm_chain, ierr)
    call enforce_pos_amp(chaindir, residuals, inv_N_noise, fg_amp, index_map, burnin)

  end subroutine enforce_pos_amps

  subroutine set_exclude_fg_amp(exclude_amp)
    implicit none

    logical(lgt), intent(in) :: exclude_amp

    exclude_pos_fg_amp = exclude_amp
    call mpi_bcast(34,                 1, MPI_INTEGER, root, comm_chain, ierr)
    call mpi_bcast(exclude_pos_fg_amp, 1, MPI_LOGICAL, root, comm_chain, ierr)

  end subroutine set_exclude_fg_amp

  subroutine set_sample_temp_coeffs(sample_coeffs)
    implicit none

    logical(lgt), intent(in) :: sample_coeffs

    sample_temp_coeffs = sample_coeffs
    call mpi_bcast(35,                 1, MPI_INTEGER, root, comm_chain, ierr)
    call mpi_bcast(sample_temp_coeffs, 1, MPI_LOGICAL, root, comm_chain, ierr)

  end subroutine set_sample_temp_coeffs


  subroutine sample_gains(s, chain, iter)
    implicit none

    type(genvec), intent(in) :: s
    integer(i4b), intent(in) :: chain, iter

    call mpi_bcast(3,                  1, MPI_INTEGER, root, comm_chain, ierr)
    call sample_gain(handle, map_id, s, chain, iter)

  end subroutine sample_gains

  subroutine sample_bandpasses(s, index_map, chain, iter)
    implicit none

    type(genvec),                      intent(in) :: s
    real(dp),     dimension(0:,1:,1:), intent(in) :: index_map
    integer(i4b),                      intent(in) :: chain, iter

    call mpi_bcast(40,                  1, MPI_INTEGER, root, comm_chain, ierr)
    call sample_bandpass(handle, map_id, s, index_map, chain, iter)

  end subroutine sample_bandpasses


  subroutine output_ml_map(paramfile)
    implicit none

    character(len=*),                  intent(in) :: paramfile

    call mpi_bcast(50,                  1, MPI_INTEGER, root, comm_chain, ierr)
    call output_ml_map_engine(myid_chain, paramfile)

  end subroutine output_ml_map

  subroutine free_slaves
    implicit none

    call mpi_bcast(0, 1, MPI_INTEGER, root, comm_chain, ierr)

  end subroutine free_slaves


end module comm_mp_mod
