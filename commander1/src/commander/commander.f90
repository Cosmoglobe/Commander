program commander
  use comm_signal_mod
  use comm_mp_mod
  use comm_chisq_mod
  use comm_mcmc_mod
  use comm_Cl_sampling_mod
  use comm_cgd_mod
  use comm_noisesamp_mod
  use comm_global_par_mod
  implicit none


  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2013, all rights reserved                *
  ! *                                                                   *
  ! *                                                                   *
  ! *   NB! The code is provided as is, and *no* guarantees are given   *
  ! *       as far as either accuracy or correctness goes. Even though  *
  ! *       it is fairly well tested, there may be (and likely are)     *
  ! *       bugs in this code.                                          *
  ! *                                                                   *
  ! *  If used for published results, please cite these papers:         *
  ! *                                                                   *
  ! *      - Jewell et al. 2004, ApJ, 609, 1                            *
  ! *      - Wandelt et al. 2004, Phys. Rev. D, 70, 083511              *
  ! *      - Eriksen et al. 2004, ApJS, 155, 227 (Commander)            *
  ! *      - Eriksen et al. 2008, ApJ, 676, 10  (Joint FG + CMB)        *
  ! *                                                                   *
  ! *********************************************************************


  ! ******************************************************
  ! *                     Conventions                    *
  ! *                                                    *
  ! *  1) All maps are stored internally in RING format  *
  ! *                                                    *
  ! ******************************************************

  integer(i4b)        :: iargc, ierr, numprocs, myid, unit, root, num_chain, num_chain_per_realization
  integer(i4b)        :: i, j, k, l, m, s, lmax_sht, nside_sht, nmaps_sht, base_see
  integer(i4b)        :: mychain, num_realizations, chain
  integer(i4b)        :: precompute_plms, nband, num_proc_per_band
  integer(i4b)        :: irow, icol, verbosity, comm_master, band_id, base_seed
  integer(i4b)        :: comm_chain, comm_alms, comm_data, myid_chain, myid_alms, myid_data
  logical(lgt)        :: slave_work, master, polar
  type(planck_rng)    :: rng_handle
  character(len=128)  :: paramfile, filename, solution_method, chain_dir, operation

  real(dp),     allocatable, dimension(:,:)         :: cl_i
  real(dp),     allocatable, dimension(:,:)         :: weights
  real(dp),     allocatable, dimension(:,:,:)       :: fg_param_map
  integer(i4b), dimension(MPI_STATUS_SIZE)          :: status

  real(dp) :: nu

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)

  root         = 0
  unit         = getlun()

  if (iargc() == 0) then

     if (myid == root) then
        write(*,*) 'Usage: commander [infofile]'
     end if
 
     call mpi_finalize(ierr)
     stop
  end if


  ! **************************************************************
  ! *          Get parameters and set up working groups          *
  ! **************************************************************

  ! Construct working groups out of all working groups
  call getarg(1,paramfile)
  call get_parameter(paramfile, 'NUM_GROUPS',                par_int=num_chain)
  call get_parameter(paramfile, 'NUM_CHAIN_PER_REALIZATION', par_int=num_chain_per_realization)
  call get_parameter(paramfile, 'NUMBAND',                   par_int=nband)
  call get_parameter(paramfile, 'NUM_PROC_PER_BAND',         par_int=num_proc_per_band)
  call get_parameter(paramfile, 'VERBOSITY',                 par_int=verbosity)
  call get_parameter(paramfile, 'SOLUTION_METHOD',           par_string=solution_method)
  call get_parameter(paramfile, 'NUM_REALIZATIONS',          par_int=num_realizations)
  call get_parameter(paramfile, 'OPERATION',                 par_string=operation)

  ! Check that the number of processors match the number of chains and maps
  if (num_chain*nband*num_proc_per_band /= numprocs) then
     if (myid == root .and. verbosity > 0) then
        write(*,*) ''
        write(*,*) 'ERROR: Number of processors does not match number of chains and maps'
        write(*,*) '       num_groups        = ', num_chain
        write(*,*) '       numband           = ', nband
        write(*,*) '       num_proc_per_band = ', num_proc_per_band
        write(*,*) '       numprocs          = ', numprocs, ' /= ', &
             & num_chain*nband*num_proc_per_band
        write(*,*) ''
     end if
     
     call mpi_finalize(ierr)
     stop
  end if

  ! Divide the processors into working groups
  icol = myid/(nband*num_proc_per_band)
  irow = mod(myid,nband*num_proc_per_band)
  call mpi_comm_split(MPI_COMM_WORLD, icol, irow, comm_chain,  ierr) ! communicator for each chain
  call mpi_comm_split(MPI_COMM_WORLD, irow, icol, comm_master, ierr) ! communicator for masters
  master = (mod(myid,nband*num_proc_per_band) == 0)

  call mpi_comm_rank(comm_chain, myid_chain, ierr)
  irow = myid_chain/num_proc_per_band
  icol = mod(myid_chain, num_proc_per_band)
  call mpi_comm_split(comm_chain, irow, icol, comm_alms, ierr) ! communicator for mpi_alm_tools
  call mpi_comm_split(comm_chain, icol, irow, comm_data, ierr) ! communicator for data nodes
  call mpi_comm_rank(comm_alms, myid_alms, ierr)
  call mpi_comm_rank(comm_data, myid_data, ierr)
  mychain = myid/(nband*num_proc_per_band) + 1

  ! Read input parameters from external file
  if (myid == root .and. verbosity > 2) write(*,*) 'Reading parameters from ', trim(paramfile)

  call get_parameter(paramfile, 'LMAX',               par_int=lmax_sht)
  call get_parameter(paramfile, 'NSIDE',              par_int=nside_sht)
  call get_parameter(paramfile, 'BASE_SEED',          par_int=base_seed)
  call get_parameter(paramfile, 'PRECOMPUTE_PLMS',    par_int=precompute_plms)
  call get_parameter(paramfile, 'POLARIZATION',       par_lgt=polar)
  call get_parameter(paramfile, 'CHAIN_DIRECTORY',    par_string=chain_dir)
  if (polar) then
     nmaps_sht = 3
  else
     nmaps_sht = 1
  end if

  ! Create chain directory if it does not exist
  if (myid == root) then
     call mkdir(trim(chain_dir))
     call mkdir(trim(chain_dir)//'/objects')
  end if

  ! Output a little information to notify the user that something is happening
  if (myid == root .and. verbosity > 0) then
     write(*,*) ''
     write(*,*) '       **********   Commander   *************'
     write(*,*) ''
     write(*,*) '          Verbosity                         = ', verbosity
     write(*,*) ''
     write(*,*) '          Number of processor groups        = ', num_chain
     write(*,*) '          Number of chains per realization  = ', num_chain_per_realization
     write(*,*) '          Number of individual channels     = ', nband
     write(*,*) '          Number of processors per band     = ', num_proc_per_band
     write(*,*) '          Number of processors              = ', numprocs
     write(*,*) '          Nside                             = ', nside_sht
     write(*,*) '          Lmax                              = ', lmax_sht
     write(*,*) '          Polarization                      = ', polar
     write(*,*) ''
  end if

  ! Distribute essential information to nodes
  if (myid == root .and. verbosity > 1) write(*,*) '     Initializing data structures'

  ! Initialize random number generator
  if (myid == root) then
     call rand_init(rng_handle, base_seed)
     do i = 1, numprocs-1
        base_seed = nint(rand_uni(rng_handle)*1000000.d0)
        call mpi_send(base_seed, 1, MPI_INTEGER, i, 98, MPI_COMM_WORLD, ierr)
     end do
  else 
     call mpi_recv(base_seed, 1, MPI_INTEGER, root, 98, MPI_COMM_WORLD, status, ierr)
     call rand_init(rng_handle, base_seed)
  end if


  ! ************************************************
  ! *               Initialize modules             *
  ! ************************************************

  band_id = myid_data+1

  call initialize_genvec_mod(paramfile)

  ! Initialize foregrounds module
  if (myid_chain == root .and. verbosity > 2) write(*,*) 'Initializing foreground modules'
  call initialize_bp_mod(myid, chain, comm_chain, paramfile, rng_handle)

!!$  if (myid == root) then
!!$     write(*,*) '1-0'
!!$     do i = 1, numband
!!$        write(*,*) trim(bp(i)%label), get_bp_line_ant(i, 115.27d9) *  compute_ant2thermo_single(115.27d9)
!!$     end do
!!$     write(*,*) '2-1'
!!$     do i = 1, numband
!!$        write(*,*) trim(bp(i)%label), get_bp_line_ant(i, 230.54d9) *  compute_ant2thermo_single(230.54d9)
!!$     end do
!!$     write(*,*) '3-2'
!!$     do i = 1, numband
!!$        write(*,*) trim(bp(i)%label), get_bp_line_ant(i, 345.80d9) *  compute_ant2thermo_single(345.80d9)
!!$     end do
!!$  end if

  call initialize_fg_mod(mychain, comm_chain, comm_alms, rng_handle, paramfile)
  if (.false.) then
     if (myid_chain == root) then
        write(*,*) 'Writing debug spectra'
        call output_debug_fg_spectra
     else
!!$        fg_components(1)%p_rms(1) = 0.1d0
!!$        fg_components(1)%p_rms(2) = 1.d0
!!$        call update_eff_fg_spectrum(frequencies)
!!$        fg_components(1)%p_rms(1) = 0.2d0
!!$        fg_components(1)%p_rms(2) = 2.d0
!!$        call update_eff_fg_spectrum(frequencies)
!!$        fg_components(1)%p_rms(1) = 0.3d0
!!$        fg_components(1)%p_rms(2) = 3.d0
!!$        call update_eff_fg_spectrum(frequencies)
     end if
     call mpi_finalize(ierr)
     stop
  end if

  if (myid_alms == root) then

     ! Initialize the Cl utility module (binning, priors etc.)
     call initialize_Cl_util_mod(paramfile)

     ! Initialize the Cl sampler
     call initialize_Cl_sampling_mod(paramfile)

     ! Initialize the beam convolution module
     if (myid_chain == root .and. verbosity > 2) write(*,*) 'Initializing beam module'
     call initialize_beam_mod(paramfile, comm_alms, comm_chain, band_id)

     ! Initialize the mpi_alm_tools module
     if (myid_chain == root .and. verbosity > 2) write(*,*) 'Initializing mpi_alm_tools'
     allocate(weights(1:2*nside_sht,nmaps_sht))
     weights = 1.d0
     call mpi_initialize_alm_tools(comm_alms, nside_sht, lmax_sht, &
          & lmax_sht, [0.d0, 0.d0], polar, precompute_plms, weights)
     deallocate(weights)

     ! Initialize the invN multiplication module
     if (myid_chain == root .and. verbosity > 2) write(*,*) 'Initializing noise module'
     call initialize_N_mult_mod(paramfile, comm_alms, comm_chain, band_id)

     ! Initialize the data module
     if (myid_chain == root .and. verbosity > 2) write(*,*) 'Initializing data module'
     base_seed = nint(rand_uni(rng_handle)*1000000.d0)
     call initialize_data_mod(paramfile, comm_chain, comm_data, comm_alms, base_seed, band_id)
     ! Initialize signal module
     if (myid_chain == root .and. verbosity > 2) write(*,*) 'Initializing signal module'
     if (master) call initialize_signal_mod(paramfile)

     if (size(pixels) == 12*nside_sht**2 .and. precompute_plms == 1) then
        allocate(plms(0:nside_sht*(lmax_sht+1)*(lmax_sht+2)-1,nmaps_sht))
        call plm_gen(nside_sht, lmax_sht, lmax_sht, plms)
     end if

  else

     ! Initialize the beam convolution module
     call initialize_beam_mod(paramfile, comm_alms, comm_chain, band_id)

     ! Initialize the mpi_alm_tools module
     call mpi_initialize_alm_tools(comm_alms, nside, lmax, &
          & lmax, [0.d0, 0.d0], polar, precompute_plms)

     ! Initialize the invN multiplication module
     call initialize_N_mult_mod(paramfile, comm_alms, comm_chain, band_id)

     ! Initialize data module
     base_seed = nint(rand_uni(rng_handle)*1000000.d0)
     call initialize_data_mod(paramfile, comm_chain, comm_data, comm_alms, base_seed)

  end if

  ! Initialize CG constraint module 
!  call initialize_cgd_constraint_mod(paramfile, rng_handle, comm_chain, comm_data, comm_alms)

  ! Initialize Conjugate Gradients module
  if (myid_chain == root .and. verbosity > 2) write(*,*) 'Initializing CG module'
  call initialize_cgd_mod(mychain, comm_chain, comm_data, comm_alms, paramfile)

  ! Initialize Direct Solution module
  if (trim(solution_method) == 'brute_force' .and. myid_chain == root) then
     call initialize_direct_solution_mod(paramfile)
  end if

  ! Initialize the MP communication module
  call initialize_mp_mod(comm_chain, rng_handle)

  ! Initialize signal multiplication module
  call initialize_S_mult_mod(lmax, nmaps)
  
  ! Initialize noise RMS amplitude module
  call initialize_noisesamp_mod(paramfile)

  ! Initialize the CGD matrix multiplication module
  base_seed = nint(rand_uni(rng_handle)*1000000.d0)
  if (myid_chain == root .and. verbosity > 2) write(*,*) 'Initializing CG multiplication module'
  call initialize_cgd_matmul_mod(paramfile, comm_alms, &
       & comm_chain, base_seed, npix, nmaps, lmax)

  ! Initialize chi-square module
  call initialize_chisq_mod(paramfile, comm_chain, comm_alms, comm_data)

  ! Set up initial power spectrum for preconditioner evaluation only
  if (master) then
     allocate(cl_i(0:lmax, nspec))
     cl_i = 0.d0
     cl_i(2:lmax,1) = 1.d0
     if (nspec > 1) then
        cl_i(2:lmax,4) = 1.d0
        cl_i(2:lmax,6) = 1.d0
     end if
     call update_S_mult_mod(cl_i)
  end if

  ! Set up foreground spectral index maps if required
  if (myid_chain == root .and. verbosity > 2) write(*,*) 'Initializing spectral response maps'
  if (sample_fg_pix) then
     if (master) then
        allocate(fg_param_map(0:npix-1, nmaps, num_fg_par))
        call initialize_index_map(paramfile, rng_handle, fg_param_map)
        call update_fg_pix_response_map(map_id, pixels, fg_pix_spec_response, fg_param_map)
        deallocate(fg_param_map)
     else
        call update_fg_pix_response_map(map_id, pixels, fg_pix_spec_response)
     end if
  end if
  !call update_cgd_constraints

  if (solution_method == 'CG') then
     ! Initialize CG preconditioner module
     if (myid_chain == root .and. verbosity > 2) write(*,*) 'Initializing CG preconditioner module'
     call pre_initialize_preconditioner(chain_dir, comm_chain, comm_data, comm_alms, comm_master, &
          & paramfile)
     if (master) then
        call free_slaves
     else
        slave_work = .true.
        do while (slave_work)
           call slave_carry_out_operation(slave_work)
        end do
     end if
  end if

  ! Set up random power spectrum for sampling
  if (master) then
     call initialize_ran_powspec(paramfile, rng_handle, cl_i)
     call update_S_mult_mod(cl_i)     
  end if

  if (myid_chain == root) then
     ! Initialize MCMC module
     if (myid_chain == root .and. verbosity > 2) write(*,*) 'Initializing MCMC module'
     call initialize_mcmc_mod(paramfile, comm_master)
     !call initialize_cl_spline_mod(paramfile)
  end if




  ! **************************************************************
  ! *                   Carry out computations                   *
  ! **************************************************************

  if (myid == root .and. verbosity > 0) write(*,*) '     Starting Gibbs sampling'

  do i = 1, num_realizations

     if ((mychain-1)/num_chain_per_realization /= mod(i-1, num_chain / num_chain_per_realization)) cycle
     chain = (i-1)*num_chain_per_realization + mod(mychain-1,num_chain_per_realization)+1

     ! Read maps
     call read_map_realization(paramfile, i)

     if (master) then
        
        ! Set up random power spectrum
        call initialize_ran_powspec(paramfile, rng_handle, cl_i)
        call update_S_mult_mod(cl_i)
        
        ! Compute chain
        call compute_single_chain(paramfile, chain, rng_handle, cl_i)
        
        call free_slaves
     
     else
        
        slave_work = .true.
        do while (slave_work)
           call slave_carry_out_operation(slave_work)
        end do
        
     end if
     
  end do

  ! Wait for everybody to exit
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! Clean up
  if (myid == root .and. verbosity > 1) write(*,*) '     Cleaning up and finalizing'
  call clean_up_N_mult_mod
  call clean_up_S_mult_mod
  call mpi_cleanup_alm_tools
  !call cleanup_cgd_constraint_mod

  if (master) then
     call cleanup_signal_mod
     call cleanup_mcmc_mod
  end if


  ! Clean up arrays
  if (allocated(cl_i))  deallocate(cl_i)


  ! And exit
  call mpi_finalize(ierr)


contains


  ! Main computational routine
  subroutine compute_single_chain(paramfile, chain, rng_handle, cl_i)
    implicit none

    integer(i4b),                         intent(in)    :: chain
    type(planck_rng),                     intent(inout) :: rng_handle
    real(dp),           dimension(0:,1:), intent(inout) :: cl_i 
    character(len=128),                   intent(in)    :: paramfile

    integer(i4b)       :: ierr, numprocs, myid, root, precond_type
    integer(i4b)       :: iter, i, j, k, l, m, num_gibbs_iter, num_step_ml_search
    integer(i4b)       :: first_iteration, ind, iter_thin, thinning_factor, stat, num_bp_step, skip_freq
    logical(lgt)       :: condition_on_cl, enforce_zero_cl, sample_inside_mask, exist, sample_bp
    logical(lgt)       :: output_ml_map_and_covmat
    real(dp)           :: t1, t2, chisq_fullsky, chisq_highlat
    character(len=256) :: filename, chain_dir, operation, cmb_md_maskfile
    type(genvec)       :: s_i, s_i_fg, s_prev

    real(dp),     allocatable, dimension(:,:)   :: cl_prev, cmb_md_mask
    real(dp),     allocatable, dimension(:,:,:) :: fg_amp_lowres, fg_param_map, fg_param_map_prev
    integer(i4b),              dimension(3)     :: l_sub, num_sub_samples
    character(len=256), allocatable, dimension(:,:) :: det_defs

    real(dp), allocatable, dimension(:,:,:) :: debug

    ! Read parameters
    call get_parameter(paramfile, 'OPERATION',                     par_string=operation)
    call get_parameter(paramfile, 'NUM_GIBBS_ITER',                par_int=num_gibbs_iter)
    call get_parameter(paramfile, 'THINNING_FACTOR',               par_int=thinning_factor)
    call get_parameter(paramfile, 'ENFORCE_ZERO_CL',               par_lgt=enforce_zero_cl)
    call get_parameter(paramfile, 'CONDITION_ON_CL',               par_lgt=condition_on_cl)
    call get_parameter(paramfile, 'SAMPLE_INSIDE_MASK',            par_lgt=sample_inside_mask)
    call get_parameter(paramfile, 'CHAIN_DIRECTORY',               par_string=chain_dir)
    call get_parameter(paramfile, 'VERBOSITY',                     par_int=verbosity)
    call get_parameter(paramfile, 'CMB_MONO_DIPOLE_MASK',          par_string=cmb_md_maskfile)
    call get_parameter(paramfile, 'NUM_STEP_WITH_ML_SEARCH',       par_int=num_step_ml_search)
    call get_parameter(paramfile, 'NUM_BP_MCMC_SUBSTEPS',          par_int=num_bp_step)
    call get_parameter(paramfile, 'NUM_FAST_SUBSTEPS',             par_int=skip_freq)
    call get_parameter(paramfile, 'OUTPUT_ML_MAP_AND_COVMAT',      par_lgt=output_ml_map_and_covmat)

    ! Check if there are masked pixels; if not, do not sample "inside mask"
    if (sample_inside_mask .and. nu_highlat == nu_fullsky) then
       write(*,*) 'Warning: Turning off sampling inside mask because there are no masked pixels'
       sample_inside_mask = .false.
    end if

    if (trim(cmb_md_maskfile) == 'fullsky') then
       allocate(cmb_md_mask(0:npix-1,nmaps))
       cmb_md_mask = 1.d0
    else if (trim(cmb_md_maskfile) == 'none') then
       
    else
       allocate(cmb_md_mask(0:npix-1,nmaps))
       call read_map(cmb_md_maskfile, cmb_md_mask)
    end if
 

    ! Allocate data structures
    allocate(cl_prev(0:lmax, nspec))
    call allocate_genvec(s_i)
    call allocate_genvec(s_i_fg)
    call allocate_genvec(s_prev)
    cl_prev     = cl_i

    ! Initialize template amplitudes
    do i = 1, numband
       do j = 1, num_fg_temp
          s_i%temp_amp(j,i) = tempamp(j,i)
       end do
    end do

    ! Initialize foreground amplitudes
    call init_fg_amps(paramfile, s_i%fg_amp)

    if (sample_fg_pix) then
       allocate(fg_amp_lowres(0:npix-1, nmaps, num_fg_comp))
       allocate(fg_param_map(0:npix-1, nmaps, num_fg_par))
       allocate(fg_param_map_prev(0:npix-1, nmaps, num_fg_par))
       call initialize_index_map(paramfile, rng_handle, fg_param_map)
       fg_param_map_prev = fg_param_map
       call update_fg_pix_response_maps(fg_param_map)
       if (.not. all(fg_pix_spec_response == fg_pix_spec_response)) then
          do j = 1, num_fg_comp
             do k = 1, nmaps
                do i = 0, map_size-1 
                   if (is_nan(fg_pix_spec_response(i,k,j))) then
                      write(*,*) i, k, j, fg_pix_spec_response(i,k,j)
                   end if
                end do
             end do
          end do
          write(*,*) 'nan'
          stop
       end if
    end if

    ! Set up template amplitude constraints
    !call update_cgd_constraint_module

    ! Initialize chains
    first_iteration = 1
    call initialize_chain_files(chain_dir, chain, num_gibbs_iter)

    !call output_component_maps(map_id, s_i, chain, 10, chain_dir)
    !call output_sample(paramfile, 2, 10, s_i, skip_freq, cl_i, fg_param_map, noiseamp, bp%gain, bp%delta)
!    stop
    if (output_ml_map_and_covmat) then
       call compute_residuals(s_i, .false.)
       call output_ml_map(paramfile)
       return
    end if

    ! Do the Gibbs sampling
    iter = first_iteration
    do while (iter <= num_gibbs_iter)

       if (mod(iter,100) == 0) write(*,*) 'Chain no. ', chain, ' -- generating sample no. ', iter
       if (verbosity > 0) write(*,*) 'Chain no. ', chain, ' -- generating sample no. ', iter
       call cpu_time(t1)

       !call set_sample_temp_coeffs(iter >= 2)
       call set_sample_temp_coeffs(.true.)
       
       if (.false. .and. num_step_ml_search > 0 .and. iter == 1) then
          call compute_lowres_residual(s_i)
          call init_ind_by_nonlin_search(residuals_lowres, inv_N_scaled, fg_param_map, s_i%fg_amp)
          call update_fg_pix_response_maps(fg_param_map)
          !call update_cgd_constraint_module
       end if
       call genvec_set_equal(s_i, s_prev)
       
       stat = 0
       
!!$       write(*,*) 'a'
!!$       call output_sample(paramfile, 2, 0, s_i, skip_freq, cl_i, fg_param_map, noiseamp, bp%gain, bp%delta)
!!$       stop

       ! ********************************************
       !             Sample amplitudes
       ! ********************************************
       
       if (verbosity > 1) write(*,*) 'Chain no. ', chain, ' -- sampling signal with hard mask'
       
       if (enforce_zero_cl) then
          cl_i = 0.d0
          call update_S_mult_mod(cl_i)
       end if

       ! Outside mask
       !if (.not. enforce_zero_cl) then
       if (.true.) then
          precond_type       = 2
          call update_mask_state(OUTSIDE_MASK)
          if (num_fg_comp > 0) then
             call set_exclude_fg_amp(&
                  & (any(fg_components%enforce_positive_amplitude) .and. any(s_i%fg_amp /= 0.d0)) .or. &
                  & (iter <= num_step_ml_search))
          else
             call set_exclude_fg_amp(.true.)
          end if

          call compute_residuals(s_i, .false.)
          call sample_signal_component(iter, precond_type, s_i, stat) ! BUG HERE

          do j = 1, num_fg_comp
             where (mask_lowres < 0.5d0) s_i%fg_amp(:,:,j) = s_prev%fg_amp(:,:,j)
             if (fg_components(j)%enforce_positive_amplitude) s_i%fg_amp(:,:,j) = s_prev%fg_amp(:,:,j)
             
             ! NEW - Overwrite sampled (but frozen) amplitude with previous value
             if (.not. fg_components(j)%sample_amplitudes) s_i%fg_amp(:,:,j) = s_prev%fg_amp(:,:,j)
          end do

          if (num_fg_temp > 0) then
             where (fix_temp)
                s_i%temp_amp = s_prev%temp_amp
             end where
          end if
             
             ! Enforce zero CMB monopole and dipole if requested; transfer CMB offsets to template coefficients
          if (.not. enforce_zero_cl .and. sample_fg_pix) call set_pix_cmb_equal_to_Cl_cmb(s_i)
          call set_exclude_fg_amp(.false.)
       end if
       

!!$       write(*,*) 'b2'
!!$       call output_sample(paramfile, 2, 2, s_i, skip_freq, cl_i, fg_param_map, noiseamp, bp%gain, bp%delta)    

       ! Do components with positivity prior
       if (sample_fg_pix .and. stat == 0 .and. sample_T_modes) then
          if (any(fg_components%enforce_positive_amplitude)) then
             call compute_lowres_residual(s_i)
             call enforce_pos_amps(chain_dir, residuals_lowres, inv_N_scaled, &
                  & s_i%fg_amp, fg_param_map, iter==1) 
          end if
       end if
!!$       write(*,*) 'c'
!!$       call output_sample(paramfile, 2, 3, s_i, skip_freq, cl_i, fg_param_map, noiseamp, bp%gain, bp%delta)
       if (allocated(cmb_md_mask) .and. sample_T_modes) call enforce_zero_cmb_md(s_i, cmb_md_mask(:,1))


!!$       write(*,*) 'd'
!!$       call output_sample(paramfile, 2, 4, s_i, skip_freq, cl_i, fg_param_map, noiseamp, bp%gain, bp%delta)
    
       do j = 1, nmaps
          do i = 1, numcomp
             if (is_nan(s_i%cmb_amp(i,j))) stat = stat+1
          end do
       end do

       ! Various sanity checks
       do k = 1, num_fg_comp
          do j = 1, nmaps
             do i = 0, npix-1
                if (is_nan(s_i%fg_amp(i,j,k))) stat = stat+1
             end do
          end do
       end do

       ! **********************************************
       !             Sample spectral indices
       ! **********************************************
       if (num_fg_par > 0 .and. stat == 0) then
          if (verbosity > 1) write(*,*) 'Chain no. ', chain, ' -- sampling spectral indices'
          call compute_lowres_residual(s_i)
          call compute_lowres_fg_amp(s_i%fg_amp, fg_amp_lowres)
          call sample_spectral_param_map(s_i, residuals_lowres, inv_N_scaled, &
               & fg_amp_lowres, fg_param_map, stat)
          if (stat == 0) call update_fg_pix_response_maps(fg_param_map)
       end if
!!$       write(*,*) 'e'
!!$       call output_sample(paramfile, 2, 5, s_i, skip_freq, cl_i, fg_param_map, noiseamp, bp%gain, bp%delta)
    
!       stop

       ! **********************************************
       !             Sample global parameters
       ! **********************************************
       if (.false. .and. num_fg_par > 0 .and. stat == 0) then
          if (verbosity > 1) write(*,*) 'Chain no. ', chain, ' -- sampling global parameters'
          call compute_lowres_residual(s_i)
          call compute_lowres_fg_amp(s_i%fg_amp, fg_amp_lowres)
          call sample_global_fg_par(rng_handle, s_i, residuals_lowres, inv_N_scaled, fg_param_map, stat)
          if (stat == 0) call update_fg_pix_response_maps(fg_param_map)
       end if

!!$       write(*,*) 'f'
!!$       call output_sample(paramfile, 2, 6, s_i, skip_freq, cl_i, fg_param_map, noiseamp, bp%gain, bp%delta)
!       stop


       ! **********************************************
       !             Sample CMB power spectrum
       ! **********************************************
       if (.not. enforce_zero_cl .and. .not. condition_on_cl .and. stat == 0) then
          if (verbosity > 1) write(*,*) 'Chain no. ', chain, ' -- sampling power spectrum'
          if (trim(operation) == 'sample') then
             call sample_Cls(s_i%cmb_amp, rng_handle, cl_i)
          else if (trim(operation) == 'optimize') then
             call compute_sigma_l(s_i%cmb_amp, cl_i)
             call bin_spectrum(cl_i)
          else
             write(*,*) 'Unknown operation: ', trim(operation)
             stop
          end if
          call update_S_mult_mod(cl_i)
          
          ! Sample high-l CMB modes quickly
          call sample_cls_and_alms_by_mcmc(rng_handle, cl_i, s_i)
          call set_pix_cmb_equal_to_Cl_cmb(s_i)
          call update_S_mult_mod(cl_i)
       end if
       !stop
       ! **********************************************
       !             Sample instrumental parameters
       ! **********************************************
       call update_mask_state(OUTSIDE_MASK)
       if (.true.) then
          ! Sample noise RMS amplitude
          if (any(N_prior(:,2) /= 0.d0)) then
             if (verbosity > 1) write(*,*) 'Chain no. ', chain, ' -- sampling noise RMS amplitudes'
             call sample_noiseamps(rng_handle, s_i, noiseamp)
             call set_noiseamps(noiseamp)
             do j = 1, numband
                inv_N_scaled(:,:,j) = inv_N_lowres(:,:,j) / noiseamp(j)**2
             end do
          end if
!!$          write(*,*) 'f'
!!$          call output_sample(paramfile, 2, 6, s_i, skip_freq, cl_i, fg_param_map, noiseamp, bp%gain, bp%delta)

          ! Sample gain
          if (any(g_gauss(:,2) /= 0.d0)) then
             if (verbosity > 1) write(*,*) 'Chain no. ', chain, ' -- sampling gains'
             call sample_gains(s_i, chain, iter)
             if (stat == 0) call update_fg_pix_response_maps(fg_param_map)
          end if

!!$          write(*,*) 'g'
!!$          call output_sample(paramfile, 2, 7, s_i, skip_freq, cl_i, fg_param_map, noiseamp, bp%gain, bp%delta)
!!$          stop

          ! Sample bandpass errors
          if (any(bp%delta_rms > 0.d0) .and. num_bp_step /= 0 .and. skip_freq > 0) then
             if (mod(iter,skip_freq) == 0) then
                if (verbosity > 1) write(*,*) 'Chain no. ', chain, ' -- sampling bandpass errors'
                call sample_bandpasses(s_i, fg_param_map, chain, iter)
             end if
          end if

!!$          write(*,*) 'g'
!!$          call output_sample(paramfile, 2, 7, s_i, skip_freq, cl_i, fg_param_map, noiseamp, bp%gain, bp%delta)
       end if

       if (iter <= num_step_ml_search .and. skip_freq > 0) then
          if (mod(iter,skip_freq) == 0) then
             call compute_lowres_residual(s_i)
             call init_ind_by_nonlin_search(residuals_lowres, inv_N_scaled, fg_param_map, s_i%fg_amp)
             call update_fg_pix_response_maps(fg_param_map)
          end if
       end if

!!$       write(*,*) 'i'
!!$       call output_sample(paramfile, 2, 7, s_i, skip_freq, cl_i, fg_param_map, noiseamp, bp%gain, bp%delta)

       ! Various sanity checks
       do k = 1, num_fg_par
          do j = 1, nmaps
             do i = 0, npix-1
                if (is_nan(fg_param_map(i,j,k))) stat = stat+1
             end do
          end do
       end do
       
       do j = 1, nspec
          do l = 0, lmax
             if (is_nan(cl_i(l,j))) stat = stat+1
          end do
       end do


       if (stat == 0) then

          ! Store sample to disk
          if (mod(iter,thinning_factor) == 0) then
             call output_sample(paramfile, chain, iter, s_i, skip_freq, cl_i, fg_param_map, noiseamp, bp%gain, bp%delta)
          end if

          ! Move to next sample
          call genvec_set_equal(s_i, s_prev)
          cl_prev           = cl_i
          iter              = iter + 1
          if (allocated(fg_param_map)) fg_param_map_prev = fg_param_map

       else

          if (verbosity > 0) write(*,*) 'Chain no. ', chain, ' -- sample rejected = ', iter
          call genvec_set_equal(s_prev, s_i)
          cl_i         = cl_prev
          if (allocated(fg_param_map)) fg_param_map = fg_param_map_prev
          call update_S_mult_mod(cl_i)
          call update_fg_pix_response_maps(fg_param_map)
          !call update_cgd_constraint_module

       end if

       if (output_ml_map_and_covmat) then
          call compute_residuals(s_i, .false.)
          call output_ml_map(paramfile)
       end if
          
       call cpu_time(t2)
       if (verbosity > 0) write(*,*) 'Chain no. ', chain, ' -- total CPU time = ', real(t2-t1,dp)

    end do

    call deallocate_genvec(s_i)

    if (sample_fg_pix) then
       deallocate(fg_amp_lowres)
       deallocate(fg_param_map)
    end if

  end subroutine compute_single_chain


  subroutine initialize_chain_files(chaindir, chain, num_gibbs_iter)
    implicit none

    character(len=*), intent(in) :: chaindir
    integer(i4b),     intent(in) :: chain, num_gibbs_iter

    integer(i4b)       :: i, unit
    character(len=4)   :: chain_text
    character(len=256) :: filename

    call int2string(chain, chain_text)
    unit = getlun()

    ! Initialize result files
    filename = trim(chain_dir) // '/' // 'chain_cls_no' // chain_text // '.unf'
    open(unit,file=trim(filename), form='unformatted', status='replace')
    write(unit) num_gibbs_iter
    close(unit)
    
    filename = trim(chain_dir) // '/' // 'chain_sigma_no' // chain_text // '.unf'
    open(unit,file=trim(filename), form='unformatted', status='replace')
    write(unit) num_gibbs_iter
    close(unit)
    
    filename = trim(chain_dir) // '/' // 'chain_alms_no' // chain_text // '.unf'
    open(unit,file=trim(filename), form='unformatted', status='replace')
    write(unit) num_gibbs_iter
    close(unit)
    
    filename = trim(chain_dir) // '/' // 'chain_fgtemp_no' // chain_text // '.unf'
    open(unit,file=trim(filename), form='unformatted', status='replace')
    write(unit) num_gibbs_iter
    close(unit)
    
    filename = trim(chain_dir) // '/' // 'chain_chisq_no' // chain_text // '.unf'
    open(unit,file=trim(filename), form='unformatted', status='replace')
    write(unit) num_gibbs_iter
    close(unit)

    filename = trim(chain_dir) // '/gain_no' // chain_text // '.dat'
    open(unit,file=trim(filename), status='replace', recl=2048)
    write(unit,fmt='(a)',advance='no') '#    Sample    '
    do i = 1, numband
       write(unit,fmt='(a)',advance='no') bp(i)%label 
    end do
    write(unit,*) 
    close(unit)

    filename = trim(chain_dir) // '/bp_no' // chain_text // '.dat'
    open(unit,file=trim(filename), status='replace', recl=2048)
    write(unit,fmt='(a)',advance='no') '#    Sample    '
    do i = 1, numband
       write(unit,fmt='(a)',advance='no') bp(i)%label 
    end do
    write(unit,*) 
    close(unit)

    filename = trim(chain_dir) // '/noiseamp_no' // chain_text // '.dat'
    open(unit,file=trim(filename), status='replace', recl=2048)
    write(unit,fmt='(a)',advance='no') '#    Sample    '
    do i = 1, numband
       write(unit,fmt='(a)',advance='no') bp(i)%label 
    end do
    write(unit,*) 
    close(unit)

    if (sample_fg_pix) then
       filename = trim(chain_dir) // '/' // 'chain_fg_amps_no' // chain_text // '.unf'
       open(unit,file=trim(filename), form='unformatted', status='replace')
       write(unit) num_gibbs_iter
       close(unit)
       
       filename = trim(chain_dir) // '/' // 'chain_fg_ind_no' // chain_text // '.unf'
       open(unit,file=trim(filename), form='unformatted', status='replace')
       write(unit) num_gibbs_iter
       close(unit)

       filename = trim(chain_dir) // '/fg_ind_mean_no' // chain_text // '.dat'
       open(unit,file=trim(filename), status='replace', recl=2048)
       write(unit,fmt='(a)',advance='no') '#    Sample    '
       do i = 1, num_fg_comp
          do j = 1, fg_components(i)%npar
             write(unit,fmt='(a,a)',advance='no') '   ', fg_components(i)%indlabel(j)
          end do
       end do
       write(unit,fmt='(a)',advance='no') 'Chi_full   Chi_hilat   nu_mask'
       write(unit,*) 
       close(unit)

       filename = trim(chain_dir) // '/par_rms_no' // chain_text // '.dat'
       open(unit,file=trim(filename), status='replace', recl=2048)
       write(unit,fmt='(a)',advance='no') '#              '
       do i = 1, num_fg_comp
          do j = 1, fg_components(i)%npar
             if (fg_components(i)%p_rms_gauss(j,2) <= 0 .or. &
                  & fg_components(i)%p_rms_uni(j,2) <= fg_components(i)%p_rms_uni(j,1)) cycle
             write(unit,fmt='(a,a)',advance='no') ' ', fg_components(i)%label
          end do
       end do
       write(unit,*)
       write(unit,fmt='(a)',advance='no') '#    Sample    '
       do i = 1, num_fg_comp
          do j = 1, fg_components(i)%npar
             if (fg_components(i)%p_rms_gauss(j,2) <= 0 .or. &
                  & fg_components(i)%p_rms_uni(j,2) <= fg_components(i)%p_rms_uni(j,1)) cycle
             write(unit,fmt='(a,a)',advance='no') ' ', fg_components(i)%indlabel(j)
          end do
       end do
       write(unit,fmt='(a)',advance='no') 'Chi_full   Chi_hilat   nu_mask'
       write(unit,*) 
       close(unit)
       
    end if

  end subroutine initialize_chain_files

  subroutine output_sample(paramfile, chain, iter, s_i, skip_freq, cl_i, fg_param_map, &
       & noiseamps, gain, delta)
    implicit none

    integer(i4b),                      intent(in) :: chain, iter, skip_freq
    character(len=*),                  intent(in) :: paramfile
    type(genvec),                      intent(in) :: s_i
    real(dp),     dimension(0:,1:),    intent(in) :: cl_i
    real(dp),     dimension(0:,1:,1:), intent(in) :: fg_param_map
    real(dp),     dimension(1:),       intent(in) :: noiseamps, gain, delta

    real(dp)           :: lat, lon, chisq_fullsky, chisq_highlat, scale
    integer(i4b)       :: i, j, l, unit
    character(len=256) :: filename
    character(len=4)   :: i_text, chain_text
    character(len=5)   :: iter_text
    integer(i4b), allocatable, dimension(:), save     :: pix
    real(dp),     allocatable, dimension(:,:), save   :: cmb_pix
    real(dp),     allocatable, dimension(:,:)   :: sigma_out, chisq_map, map, tempamps
    real(dp),     allocatable, dimension(:,:,:) :: par_smooth
    real(sp),     allocatable, dimension(:,:,:) :: fg_scaled
    integer(i4b),       save :: output_thinning_factor = 0, num_objects
    logical(lgt),       save :: output_cross_corr
    logical(lgt),       save :: output_data(5)
    character(len=256), save :: chain_dir, object_dir
    character(len=256), allocatable, dimension(:), save :: objname

    if (output_thinning_factor == 0) then
       call get_parameter(paramfile, 'OUTPUT_THINNING_FACTOR',          par_int=output_thinning_factor)
       call get_parameter(paramfile, 'OUTPUT_CROSS_CORRELATION_STATS',  par_lgt=output_cross_corr)
       call get_parameter(paramfile, 'OUTPUT_SPECTRA',                  par_lgt=output_data(1))
       call get_parameter(paramfile, 'OUTPUT_ALMS',                     par_lgt=output_data(2))
       call get_parameter(paramfile, 'OUTPUT_MAPS',                     par_lgt=output_data(3))
       call get_parameter(paramfile, 'OUTPUT_FOREGROUNDS',              par_lgt=output_data(4))
       call get_parameter(paramfile, 'OUTPUT_CHISQ',                    par_lgt=output_data(5))
       call get_parameter(paramfile, 'CHAIN_DIRECTORY',                 par_string=chain_dir)
       call get_parameter(paramfile, 'NUM_OBJECTS',                     par_int=num_objects)
       object_dir = trim(chain_dir) // '/objects'

       allocate(pix(num_objects), cmb_pix(num_objects,nmaps), objname(num_objects))
       do i = 1, num_objects
          call int2string(i, i_text)
          call get_parameter(paramfile, 'OBJECT_NAME'//i_text,      par_string=objname(i))
          call get_parameter(paramfile, 'OBJECT_LATITUDE'//i_text,  par_dp=lat)
          call get_parameter(paramfile, 'OBJECT_LONGITUDE'//i_text, par_dp=lon)
          call ang2pix_ring(nside,(90.d0-lat)*DEG2RAD, lon*DEG2RAD, pix(i))
       end do
    end if

    if (num_fg_comp > 0) then
       allocate(par_smooth(0:npix-1,nmaps,num_fg_par))
       call get_smooth_par_map(fg_param_map, par_smooth)
    end if

    if (num_fg_temp > 0) then
       allocate(tempamps(num_fg_temp,numband))
       tempamps = s_i%temp_amp
    end if

    unit = getlun()

    allocate(sigma_out(0:lmax, nspec), chisq_map(0:npix-1,nmaps))
    call int2string(iter, iter_text)
    call int2string(chain, chain_text)

    call compute_sigma_l(s_i%cmb_amp, sigma_out)
    do l = 0, lmax
       sigma_out(l,:) = sigma_out(l,:) * l*(l+1)/(2.d0*pi)
    end do

    filename = trim(chain_dir) // '/chain_cls_no' // chain_text // '.unf'
    open(unit,file=trim(filename), form='unformatted', position='append')
    write(unit) real(cl_i,sp)
    close(unit)
          
    filename = trim(chain_dir) // '/chain_sigma_no' // chain_text // '.unf'
    open(unit,file=trim(filename), form='unformatted', position='append')
    write(unit) real(sigma_out,sp)
    close(unit)

    filename = trim(chain_dir) // '/noiseamp_no' // chain_text // '.dat'
    open(unit,file=trim(filename), position='append', recl=2048)
    write(unit,fmt='(i8)',advance='no') iter
    do i = 1, numband
       write(unit,fmt='(f14.6)',advance='no') noiseamps(i)
    end do
    write(unit,*)
    close(unit)

    filename = trim(chain_dir) // '/gain_no' // chain_text // '.dat'
    open(unit,file=trim(filename), position='append', recl=2048)
    write(unit,fmt='(i8)',advance='no') iter
    do i = 1, numband
       write(unit,fmt='(f14.6)',advance='no') gain(i)
    end do
    write(unit,*)
    close(unit)

    filename = trim(chain_dir) // '/bp_no' // chain_text // '.dat'
    open(unit,file=trim(filename), position='append', recl=2048)
    write(unit,fmt='(i8)',advance='no') iter
    do i = 1, numband
       write(unit,fmt='(f14.6)',advance='no') delta(i)
    end do
    write(unit,*)
    close(unit)

    if (num_fg_temp > 0) then
       filename = trim(chain_dir) // '/chain_fgtemp_no' // chain_text // '.unf'
       open(unit,file=trim(filename), form='unformatted', position='append')
       write(unit) real(tempamps,sp)
       close(unit)
    end if
          
    call compute_chisq(map_id, s_i, chisq_map=chisq_map, chain=chain, iter=iter, &
         & output_stats=mod(iter,output_thinning_factor)==0, chisq_fullsky=chisq_fullsky, &
         & chisq_highlat=chisq_highlat)
    if (verbosity > 1) write(*,fmt='(a,i4,a,e16.8,a,e16.8)') 'Chain no. ', &
         & chain, ' -- chi^2     = ', chisq_fullsky, ', highlat =', chisq_highlat
    filename = trim(chain_dir) // '/' // 'chain_chisq_no' // chain_text // '.unf'
    open(unit,file=trim(filename), form='unformatted', position='append')
    write(unit) real(sum(chisq_map),sp)
    close(unit)

    if (sample_fg_pix) then
       filename = trim(chain_dir) // '/fg_ind_mean_no' // chain_text // '.dat'
       open(unit,file=trim(filename), position='append', recl=2048)
       k = 1
       write(unit,fmt='(i11)',advance='no') iter
       do i = 1, num_fg_comp
          do j = 1, fg_components(i)%npar
             write(unit,fmt='(e15.3)',advance='no') sum(par_smooth(:,:,k)*fg_components(i)%mask) / &
                  & sum(fg_components(i)%mask)
             k = k+1
          end do
       end do
       write(unit,fmt='(f11.3)',advance='no') chisq_fullsky / nu_fullsky
       write(unit,fmt='(f11.3)',advance='no') chisq_highlat / nu_highlat
       write(unit,fmt='(e11.3)',advance='no') real(nu_highlat,sp)
       write(unit,*) 
       close(unit)

       filename = trim(chain_dir) // '/par_rms_no' // chain_text // '.dat'
       open(unit,file=trim(filename), position='append', recl=2048)
       k = 1
       write(unit,fmt='(i11)',advance='no') iter
       do i = 1, num_fg_comp
          do j = 1, fg_components(i)%npar
             if (fg_components(i)%p_rms_gauss(j,2) <= 0 .or. &
                  & fg_components(i)%p_rms_uni(j,2) <= fg_components(i)%p_rms_uni(j,1)) cycle
             write(unit,fmt='(e11.3)',advance='no') fg_components(i)%p_rms(j)
             k = k+1
          end do
       end do
       write(unit,fmt='(f11.3)',advance='no') chisq_fullsky / nu_fullsky
       write(unit,fmt='(f11.3)',advance='no') chisq_highlat / nu_highlat
       write(unit,fmt='(e11.3)',advance='no') real(nu_highlat,sp)
       write(unit,*) 
       close(unit)
    end if

    if (mod(iter, output_thinning_factor) == 0) then
       filename = trim(chain_dir) // '/' // 'chain_alms_no' // chain_text // '.unf'
       open(unit,file=trim(filename), form='unformatted', position='append')
       write(unit) real(s_i%cmb_amp,sp)
       close(unit)
       
       if (sample_fg_pix) then
          allocate(fg_scaled(0:npix-1,nmaps,num_fg_comp))
          fg_scaled = s_i%fg_amp
          do i = 1, num_fg_comp
             fg_scaled(:,:,i) = fg_scaled(:,:,i) * ant2unit(fg_components(i)%amp_unit, &
                  & fg_components(i)%nu_ref, band_ref=fg_components(i)%ref_band)
          end do

          filename = trim(chain_dir) // '/' // 'chain_fg_amps_no' // chain_text // '.unf'
          open(unit,file=trim(filename), form='unformatted', position='append')
          write(unit) fg_scaled
          close(unit)
          deallocate(fg_scaled)
          
          allocate(fg_scaled(0:npix-1,nmaps,num_fg_par))
          fg_scaled = fg_param_map
          k = 1
          do i = 1, num_fg_comp
             do j = 1, fg_components(i)%npar
                if (trim(fg_components(i)%type) == 'freefree_EM' .and. j == 1) then
                   fg_scaled(:,:,k) = exp(fg_scaled(:,:,k))
                end if
                k = k+1
             end do
          end do

          filename = trim(chain_dir) // '/' // 'chain_fg_ind_no' // chain_text // '.unf'
          open(unit,file=trim(filename), form='unformatted', position='append')
          write(unit) fg_scaled
          close(unit)
          deallocate(fg_scaled)
       end if

       if (output_data(1)) then
          call output_spectra_from_iteration(chain, iter, chain_dir, cl_i, sigma_out)
       end if
       
       if (output_data(2)) then
          call output_alms_from_iteration(nside, lmax, chain, iter, chain_dir, s_i%cmb_amp)
       end if

       cmb_pix = 0.d0
       if (output_data(3)) then
          call output_maps_from_iteration(nside, lmax, chain, iter, chain_dir, s_i%cmb_amp, &
               & fwhm_lowres, pix, cmb_pix)
       end if

       if (num_objects > 0) then
          do i = 1, num_objects
             call output_pixel_to_file(object_dir, chain, objname(i), pix(i), iter, &
                  & cmbmaps_lowres(pix(i),:,:), inv_N_scaled(pix(i),:,:), cmb_pix(i,:), &
                  & tempamps, fg_temp_lowres(pix(i),:,:,:), &
                  & s_i%fg_amp(pix(i),:,:), par_smooth(pix(i),:,:), chisq_map(pix(i),1))
          end do
       end if

       if (output_data(4)) then
          if (num_fg_temp > 0) then 
             call output_foregrounds(chain, iter, chain_dir, tempamps)
          end if
          if (sample_fg_pix) then
             j = 1
             do i = 1, num_fg_signal
                filename = trim(chain_dir) // '/' // trim(fg_components(i)%label) // &
                     & '_c' // chain_text // '_k' // iter_text // '.fits'
                scale = ant2unit(fg_components(i)%amp_unit, fg_components(i)%nu_ref, &
                     & band_ref=fg_components(i)%ref_band)
                call write_map(filename, s_i%fg_amp(:,:,i)*scale, unit=fg_components(i)%amp_unit, &
                     & comptype=fg_components(i)%label, nu_ref=fg_components(i)%nu_ref)
                do k = 1, fg_components(i)%npar
                   filename = trim(chain_dir) // '/' // trim(fg_components(i)%label) // '_' // &
                        & trim(fg_components(i)%indlabel(k)) // '_c' // chain_text // '_k' // &
                        & iter_text // '.fits'
                   if (trim(fg_components(i)%type) == 'freefree_EM' .and. k == 1) then
                      call write_map(filename, exp(par_smooth(:,:,j)), unit=fg_components(i)%ind_unit(k), &
                           & ttype=fg_components(i)%ttype(k))                      
                   else
                      call write_map(filename, par_smooth(:,:,j), unit=fg_components(i)%ind_unit(k), &
                           & ttype=fg_components(i)%ttype(k))
                   end if
                   if (fg_components(i)%fwhm_p(k) > 0.d0) then
                      filename = trim(chain_dir) // '/' // trim(fg_components(i)%label) // '_' // &
                           & trim(fg_components(i)%indlabel(k)) // '_c' // chain_text // '_k' // &
                           & iter_text // '_nosmooth.fits'
                      if (trim(fg_components(i)%type) == 'freefree_EM' .and. k == 1) then
                         call write_map(filename, exp(fg_param_map(:,:,j)), unit=fg_components(i)%ind_unit(k), &
                              & ttype=fg_components(i)%ttype(k))                      
                      else
                         call write_map(filename, fg_param_map(:,:,j), unit=fg_components(i)%ind_unit(k), &
                              & ttype=fg_components(i)%ttype(k))
                      end if
                   end if
                   j = j+1
                end do
             end do
          end if
          
       end if
       
       if (output_data(5)) then
          call output_chisq_from_iteration(chain, iter, chain_dir, chisq_map)
       end if

    end if

    call output_component_maps(map_id, s_i, chain, iter, chain_dir)
    call output_mixing_matrices(map_id, chain, iter, chain_dir)
    if (skip_freq > 0) then
       !if(mod(iter,skip_freq) == 0) call output_fg_overview(paramfile, chain, iter, chain_dir, s_i, par_smooth)
    end if

    if (output_cross_corr .and. num_fg_comp > 0) call output_signal_correlations(s_i, chisq_map, &
         & chain, iter, chain_dir)
          
    deallocate(sigma_out, chisq_map)
    if (allocated(par_smooth)) deallocate(par_smooth)
    if (allocated(tempamps))   deallocate(tempamps)

  end subroutine output_sample

  subroutine output_debug_fg_spectra
    implicit none

    integer(i4b) :: i, j, k, q, n
    real(dp)     :: s1, s2, nu, nu_min, nu_max

    real(dp),     allocatable, dimension(:) :: spec
    real(dp), dimension(17) :: bands = [0.408d0, 23.d0, 31.d0, 41d0, 61d0, 94d0, 30d0, 44d0, 70d0, 100d0, &
         & 143d0, 217d0, 353d0, 545d0, 857d0, 1249d0, 2998d0]

    allocate(spec(numband))
    if (.true.) then

       n      = 1000
       nu_min = 10d0
       nu_max = 3000d0
!!$       open(50,file='spec_test.dat')
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          write(50,*) nu, compute_physical_dust_spectrum(nu,0.d0,fg_components(1)%us,545.d9,fg_components(1)%S_phys_dust, fg_components(1)%S_dust_coeff, [0.d0,0.d0])
!!$       end do
!!$       close(50)
!!$       stop
!!$       open(58,file='bands.dat')
!!$       do i = 1, 15
!!$          write(58,*) bands(i)*0.98, -1e12
!!$          write(58,*) bands(i)*0.98, -100
!!$          write(58,*) bands(i)*0.98, -1e-6
!!$          write(58,*) bands(i)*0.98, 1e-6
!!$          write(58,*) bands(i)*0.98, 100
!!$          write(58,*) bands(i)*0.98, 1e12
!!$          write(58,*) bands(i)*1.02, 1e12
!!$          write(58,*) bands(i)*1.02, 100
!!$          write(58,*) bands(i)*1.02, 1e-6
!!$          write(58,*) bands(i)*1.02, -1e-6
!!$          write(58,*) bands(i)*1.02, -100
!!$          write(58,*) bands(i)*1.02, -1e12
!!$          write(58,*)
!!$       end do
!!$       close(58)
       

       open(58,file='s.dat')
!!$       do j = 1, num_fg_comp
!!$          do q = 1, numband
!!$             k = i2f(q)
!!$             s1 = get_effective_fg_spectrum(fg_components(j), k, fg_components(j)%priors(:,3))
!!$             if (trim(bp(k)%unit) == 'uK_cmb') then
!!$                write(*,*) real(bp(k)%nu_c/1.d9,sp), s1 / bp(k)%a2t
!!$                write(58,*) bp(k)%nu_c/1.d9, s1 / bp(k)%a2t
!!$             else
!!$                write(*,*) real(bp(k)%nu_c/1.d9,sp), s1 * bp(k)%f2t / bp(k)%a2t
!!$                write(58,*) bp(k)%nu_c/1.d9, s1 * bp(k)%f2t / bp(k)%a2t
!!$             end if
!!$          end do
!!$          write(58,*)
!!$       end do
!!$       close(58)
!!$       stop


       n      = 1000
       nu_min = 0.1d0
       nu_max = 3000d0
       ! CMB
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = 1.d0 / compute_ant2thermo_single(nu)
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       ! TSZ
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = 1.d0/((2.d0*nu**2*k_b/c**2 / &
!!$               & (compute_bnu_prime_single(nu) * sz_thermo(nu)))) * 1d6
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       ! Synchrotron 
!!$       nu_min = 0.001d0
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_power_law_break_spectrum(nu, 0.408d9, 0.0d0, -3.0d0, 3.d0, [0.d0,0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_power_law_break_spectrum(nu, 0.408d9, 0.4d0, -3.d0, 3.d0, [0.d0,0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_power_law_break_spectrum(nu, 0.408d9, 0.4d0, -3.2d0, 3.d0, [0.d0,0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_power_law_break_spectrum(nu, 0.408d9, 0.4d0, -2.8d0, 3.d0, [0.d0,0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)

!       open(58,file='s.dat')
!       do j = 1, num_fg_comp
!          do q = 1, numband
!             k = i2f(q)
!             s1 = get_effective_fg_spectrum(fg_components(j), k, fg_components(j)%priors(:,3))
!             if (trim(bp(k)%unit) == 'uK_cmb') then
!                write(*,*) real(bp(k)%nu_c/1.d9,sp), s1 / bp(k)%a2t
!                write(58,*) bp(k)%nu_c/1.d9, s1 / bp(k)%a2t
!             else
!                write(*,*) real(bp(k)%nu_c/1.d9,sp), s1 * bp(k)%f2t / bp(k)%a2t
!                write(58,*) bp(k)%nu_c/1.d9, s1 * bp(k)%f2t / bp(k)%a2t
!             end if
!          end do
!          write(58,*)
!       end do
!       close(58)
!       stop


       n      = 1000
       nu_min = 0.1d0
       nu_max = 3000d0
       ! CMB
 !      do i = 1, n
 !         nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
 !         s1 = 1.d0 / compute_ant2thermo_single(nu)
 !         write(58,*) nu/1.d9, s1
 !      end do
 !      write(58,*)

       ! TSZ
!       do i = 1, n
!          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!          s1 = 1.d0/((2.d0*nu**2*k_b/c**2 / &
!               & (compute_bnu_prime_single(nu) * sz_thermo(nu)))) * 1d6
!          write(58,*) nu/1.d9, s1
!       end do
!       write(58,*)

       ! Synchrotron 
!       nu_min = 0.001d0
!       do i = 1, n
!          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!          s1 = compute_power_law_break_spectrum(nu, 0.408d9, 0.0d0, -3.0d0, 3.d0, [0.d0,0.d0])
!          write(58,*) nu/1.d9, s1
!       end do
!       write(58,*)

!       do i = 1, n
!          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!          s1 = compute_power_law_break_spectrum(nu, 0.408d9, 0.4d0, -3.d0, 3.d0, [0.d0,0.d0])
!          write(58,*) nu/1.d9, s1
!       end do
!       write(58,*)

!       do i = 1, n
!          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!          s1 = compute_power_law_break_spectrum(nu, 0.408d9, 0.4d0, -3.2d0, 3.d0, [0.d0,0.d0])
!          write(58,*) nu/1.d9, s1
!       end do
!       write(58,*)

!       do i = 1, n
!          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!          s1 = compute_power_law_break_spectrum(nu, 0.408d9, 0.4d0, -2.8d0, 3.d0, [0.d0,0.d0])
!          write(58,*) nu/1.d9, s1
!       end do
!       write(58,*)

!!$       close(58)
!!$       stop


       n      = 1000
       nu_min = 0.1d0
       nu_max = 3000d0
!!$       ! CMB
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = 1.d0 / compute_ant2thermo_single(nu)
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       ! TSZ
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = 1.d0/((2.d0*nu**2*k_b/c**2 / &
!!$               & (compute_bnu_prime_single(nu) * sz_thermo(nu)))) * 1d6
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
       ! Synchrotron 
       !nu_min = 0.001d0
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          !s1 = compute_power_law_break_spectrum(nu, 0.408d9, 0.0d0, -3.0d0, 3.d0, [0.d0,0.d0])
!!$          s1 = compute_power_law_spectrum(nu, 21.d9, -3.0d0, -0.053d0, [0.d0,0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_power_law_break_spectrum(nu, 0.408d9, 0.4d0, -3.d0, 3.d0, [0.d0,0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_power_law_break_spectrum(nu, 0.408d9, 0.4d0, -3.2d0, 3.d0, [0.d0,0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_power_law_break_spectrum(nu, 0.408d9, 0.4d0, -2.8d0, 3.d0, [0.d0,0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       ! Thermal dust
!!$       nu_min = 0.1d0
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_one_comp_dust_spectrum(nu, 545.d9, 1.5d0, 21d0, 100.d0, 0.d0, [0.d0, 0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_one_comp_dust_spectrum(nu, 545.d9, 2.d0, 21d0, 100.d0, 0.d0, [0.d0, 0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_one_comp_dust_spectrum(nu, 545.d9, 1.d0, 21d0, 100.d0, 0.d0, [0.d0, 0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_one_comp_dust_spectrum(nu, 545.d9, 1.5d0, 16.d0, 100.d0, 0.d0, [0.d0, 0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_one_comp_dust_spectrum(nu, 545.d9, 1.5d0, 26.d0, 100.d0, 0.d0, [0.d0, 0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       ! Spinning dust
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_AME_freq_shift_spectrum(nu, 22.8d9, 21.d0, 30.d9, fg_components(5)%S_nu_ref, [0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_AME_freq_shift_spectrum(nu, 22.8d9, 11.d0, 30.d9, fg_components(5)%S_nu_ref, [0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_AME_freq_shift_spectrum(nu, 22.8d9, 31.d0, 30.d9, fg_components(5)%S_nu_ref, [0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)


       ! Spinning dust
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_AME_freq_shift_spectrum(nu, 22.8d9, 21.d0, 30.d9, fg_components(5)%S_nu_ref, [0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_AME_freq_shift_spectrum(nu, 22.8d9, 11.d0, 30.d9, fg_components(5)%S_nu_ref, [0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_AME_freq_shift_spectrum(nu, 22.8d9, 31.d0, 30.d9, fg_components(5)%S_nu_ref, [0.d0])
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)

       ! Free-free
       nu_min = 0.00001d0
       do i = 1, n
          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
          s1 = compute_freefree_EM_spectrum(nu, log(1.d0), 7000d0)
          !s2 = compute_freefree_spectrum(nu, 23.d9, 7000d0)
          write(58,*) nu/1.d9, s1, s2
       end do
       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_freefree_EM_spectrum(nu, log(1.d-2), 7000d0)
!!$          s2 = compute_freefree_spectrum(nu, 23.d9, 7000d0)
!!$          write(58,*) nu/1.d9, s1, s2
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_freefree_EM_spectrum(nu, log(1.d2), 7000d0)
!!$          s2 = compute_freefree_spectrum(nu, 23d9, 7000d0)
!!$          write(58,*) nu/1.d9, s1, s2
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_freefree_EM_spectrum(nu, log(1.d0), 500d0)
!!$          s2 = compute_freefree_spectrum(nu, 23.d9, 500d0)
!!$          write(58,*) nu/1.d9, s1, s2
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_freefree_EM_spectrum(nu, log(1.d0), 20000d0)
!!$          s2 = compute_freefree_spectrum(nu, 23.d9, 20000d0)
!!$          write(58,*) nu/1.d9, s1, s2
!!$       end do
!!$       write(58,*)

!!$       ! CO
!!$       nu = 100d9
!!$       s1 = compute_CO_multiline_spectrum(4, fg_components(3)%co_band, 100d0, [0.6d0, 0.3d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)
!!$
!!$       nu = 217d9*0.99
!!$       s1 = compute_CO_multiline_spectrum(6, fg_components(3)%co_band, 100d0, [0.5d0, 0.3d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)
!!$
!!$       nu = 217d9
!!$       s1 = compute_CO_multiline_spectrum(6, fg_components(3)%co_band, 100d0, [0.6d0, 0.3d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)
!!$
!!$       nu = 217d9*1.01
!!$       s1 = compute_CO_multiline_spectrum(6, fg_components(3)%co_band, 100d0, [0.7d0, 0.3d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)
!!$
!!$       nu = 353d9*0.99
!!$       s1 = compute_CO_multiline_spectrum(7, fg_components(3)%co_band, 100d0, [0.6d0, 0.2d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)
!!$
!!$       nu = 353d9
!!$       s1 = compute_CO_multiline_spectrum(7, fg_components(3)%co_band, 100d0, [0.6d0, 0.25d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)
!!$
!!$       nu = 353d9*1.01
!!$       s1 = compute_CO_multiline_spectrum(7, fg_components(3)%co_band, 100d0, [0.6d0, 0.3d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)


!!$       nu_min = 0.00001d0
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_freefree_EM_spectrum(nu, log(1.d0), 7000d0)
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_freefree_EM_spectrum(nu, log(1.d-2), 7000d0)
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_freefree_EM_spectrum(nu, log(1.d2), 7000d0)
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_freefree_EM_spectrum(nu, log(1.d0), 500d0)
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       do i = 1, n
!!$          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!!$          s1 = compute_freefree_EM_spectrum(nu, log(1.d0), 20000d0)
!!$          write(58,*) nu/1.d9, s1
!!$       end do
!!$       write(58,*)
!!$
!!$       ! CO
!!$       nu = 100d9
!!$       s1 = compute_CO_multiline_spectrum(4, fg_components(3)%co_band, 100d0, [0.6d0, 0.3d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)
!!$
!!$       nu = 217d9*0.99
!!$       s1 = compute_CO_multiline_spectrum(6, fg_components(3)%co_band, 100d0, [0.5d0, 0.3d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)
!!$
!!$       nu = 217d9
!!$       s1 = compute_CO_multiline_spectrum(6, fg_components(3)%co_band, 100d0, [0.6d0, 0.3d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)
!!$
!!$       nu = 217d9*1.01
!!$       s1 = compute_CO_multiline_spectrum(6, fg_components(3)%co_band, 100d0, [0.7d0, 0.3d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)
!!$
!!$       nu = 353d9*0.99
!!$       s1 = compute_CO_multiline_spectrum(7, fg_components(3)%co_band, 100d0, [0.6d0, 0.2d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)
!!$
!!$       nu = 353d9
!!$       s1 = compute_CO_multiline_spectrum(7, fg_components(3)%co_band, 100d0, [0.6d0, 0.25d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)
!!$
!!$       nu = 353d9*1.01
!!$       s1 = compute_CO_multiline_spectrum(7, fg_components(3)%co_band, 100d0, [0.6d0, 0.3d0])
!!$       write(58,*) nu/1.d9, 1.d-6
!!$       write(58,*) nu/1.d9, s1
!!$       write(58,*)


!       nu_min = 0.00001d0
!       do i = 1, n
!          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!          s1 = compute_freefree_EM_spectrum(nu, log(1.d0), 7000d0)
!          write(58,*) nu/1.d9, s1
!       end do
!       write(58,*)

!       do i = 1, n
!          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!          s1 = compute_freefree_EM_spectrum(nu, log(1.d-2), 7000d0)
!          write(58,*) nu/1.d9, s1
!       end do
!       write(58,*)

!       do i = 1, n
!          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!          s1 = compute_freefree_EM_spectrum(nu, log(1.d2), 7000d0)
!          write(58,*) nu/1.d9, s1
!       end do
!       write(58,*)

!       do i = 1, n
!          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!          s1 = compute_freefree_EM_spectrum(nu, log(1.d0), 500d0)
!          write(58,*) nu/1.d9, s1
!       end do
!       write(58,*)

!       do i = 1, n
!          nu = nu_min * (nu_max/nu_min)**((i-1.d0)/(n-1.d0)) * 1d9
!          s1 = compute_freefree_EM_spectrum(nu, log(1.d0), 20000d0)
!          write(58,*) nu/1.d9, s1
!       end do
!       write(58,*)

       ! CO
!       nu = 100d9
!       s1 = compute_CO_multiline_spectrum(4, fg_components(3)%co_band, 100d0, [0.6d0, 0.3d0])
!       write(58,*) nu/1.d9, 1.d-6
!       write(58,*) nu/1.d9, s1
!       write(58,*)

!       nu = 217d9*0.99
!       s1 = compute_CO_multiline_spectrum(6, fg_components(3)%co_band, 100d0, [0.5d0, 0.3d0])
!       write(58,*) nu/1.d9, 1.d-6
!       write(58,*) nu/1.d9, s1
!       write(58,*)

!       nu = 217d9
!       s1 = compute_CO_multiline_spectrum(6, fg_components(3)%co_band, 100d0, [0.6d0, 0.3d0])
!       write(58,*) nu/1.d9, 1.d-6
!       write(58,*) nu/1.d9, s1
!       write(58,*)

!       nu = 217d9*1.01
!       s1 = compute_CO_multiline_spectrum(6, fg_components(3)%co_band, 100d0, [0.7d0, 0.3d0])
!       write(58,*) nu/1.d9, 1.d-6
!       write(58,*) nu/1.d9, s1
!       write(58,*)

!       nu = 353d9*0.99
!       s1 = compute_CO_multiline_spectrum(7, fg_components(3)%co_band, 100d0, [0.6d0, 0.2d0])
!       write(58,*) nu/1.d9, 1.d-6
!       write(58,*) nu/1.d9, s1
!       write(58,*)

!       nu = 353d9
!       s1 = compute_CO_multiline_spectrum(7, fg_components(3)%co_band, 100d0, [0.6d0, 0.25d0])
!       write(58,*) nu/1.d9, 1.d-6
!       write(58,*) nu/1.d9, s1
!       write(58,*)

!       nu = 353d9*1.01
!       s1 = compute_CO_multiline_spectrum(7, fg_components(3)%co_band, 100d0, [0.6d0, 0.3d0])
!       write(58,*) nu/1.d9, 1.d-6
!       write(58,*) nu/1.d9, s1
!       write(58,*)

       close(58)

    end if
     
    call mpi_finalize(ierr)
    stop

  end subroutine output_debug_fg_spectra

  subroutine output_pixel_to_file(chain_dir, chain, objname, pixel, iter, data, invN, cmb_pix, &
       & temp_amp, fg_temp, fg_amp, fg_param, chisq)
    implicit none
    
    character(len=*),                   intent(in) :: chain_dir, objname
    integer(i4b),                       intent(in) :: pixel, chain, iter
    real(dp),         dimension(:,:),   intent(in) :: data, invN, fg_param, fg_amp, temp_amp
    real(dp),         dimension(:,:,:), intent(in) :: fg_temp
    real(dp),         dimension(:),     intent(in) :: cmb_pix
    real(dp),                           intent(in) :: chisq

    real(dp)     :: f, scale
    integer(i4b) :: i, j, k, q, unit
    real(dp),     allocatable, dimension(:) :: spec_tot
    character(len=1), dimension(3) :: ind2label = ['T','Q','U']
    character(len=256) :: filename
    character(len=5)   :: iter_text
    character(len=3)   :: chain_text
    type(fg_params)    :: fg_par
    
    call int2string(chain, chain_text)
    call int2string(iter,  iter_text)
    unit = getlun()
    
    allocate(spec_tot(numband))

    do j = 1, nmaps
       if (all(invN(j,:) == 0.d0)) cycle

       filename = trim(chain_dir) // '/' // trim(objname) // '_' // &
            & ind2label(j) // '_c' // chain_text // '_k' // iter_text // '.dat'
       open(unit,file=trim(filename), recl=1024)
       write(unit,*) '# Pixel = ', pixel, ', object = ', trim(objname), ', chisq = ', real(chisq,sp)

       ! Print out data
       write(unit,*) '# Data with template corrections:'
       do i = 1, numband
          q = i2f(i)
          if (invN(j,q) > 0.d0) then
             f = data(j,q)
             do k = 1, size(temp_amp,1)
                f = f - temp_amp(k,q) * fg_temp(j,q,k)
             end do
             scale = spec2data(q, intype='uK_ant', outtype=bp(q)%unit)
             write(unit,fmt='(3e16.8)') bp(q)%nu_c/1d9, f/scale, &
                  & sqrt(1.d0/invN(j,q))/scale
          end if
       end do
       write(unit,*)

       write(unit,*) '# Data without template corrections:'
       do i = 1, numband
          q = i2f(i)
          if (invN(j,q) > 0.d0) then
             f = data(j,q)
             scale = spec2data(q, intype='uK_ant', outtype=bp(q)%unit)
             write(unit,fmt='(3e16.8)') bp(q)%nu_c/1d9, f/scale, &
                  & sqrt(1.d0/invN(j,q))/scale
          end if
       end do
       write(unit,*)

       ! Print out cmb amplitude
       spec_tot = 0.d0
       if (cmb_pix(j) /= 0.d0) then
          write(unit,*) '# CMB with Cls'
          do i = 1, numband
             q = i2f(i)
             if (invN(j,q) > 0.d0) then
                f = cmb_pix(j)/bp(q)%a2t
                spec_tot(q) = spec_tot(q) + f
                if (f > 0.d0) write(unit,fmt='(3e16.8)') bp(q)%nu_c/1d9, f, 0.d0
             end if
          end do
       end if
       write(unit,*)

       ! Print out signal components
       call reorder_fg_params(fg_param, fg_par)
       do k = 1, num_fg_comp
          scale = ant2unit(fg_components(k)%amp_unit, fg_components(k)%nu_ref, &
               & band_ref=fg_components(k)%ref_band)
          write(unit,*) '# Component = ', trim(fg_components(k)%label)
          write(unit,*) '# Amplitude = ', fg_amp(j,k), ' ', trim(fg_components(k)%amp_unit)
          write(unit,*) '# Ref freq  = ', fg_components(k)%nu_ref/1d9
          do i = 1, numband
             q = i2f(i)
             if (invN(j,q) > 0.d0) then
                scale = spec2data(q, intype='uK_ant', outtype=bp(q)%unit)
                f = get_effective_fg_spectrum(fg_components(k), q, fg_par%comp(k)%p(j,:), &
                     & pixel=pixel, pol=j) * fg_amp(j,k) / scale
                spec_tot(q) = spec_tot(q) + f
                if (f /= 0.d0) write(unit,fmt='(3e16.8)') bp(q)%nu_c/1d9, f, 0.d0
             end if
          end do
          write(unit,*)
       end do
       call deallocate_fg_params(fg_par)

       ! Print out total signal
       write(unit,*) '# Total signal component'
       do i = 1, numband
          q = i2f(i)
          if (invN(j,q) > 0.d0) then
             write(unit,fmt='(3e16.8)') bp(q)%nu_c/1d9, spec_tot(q), 0.d0
          end if
       end do
       close(unit)

    end do

    deallocate(spec_tot)

  end subroutine output_pixel_to_file

!  subroutine output_foregrounds(unit, chain, iteration, coeff, map)
  subroutine output_foregrounds(chain, iteration, chain_dir, coeff)
    implicit none

    integer(i4b),                   intent(in) :: chain, iteration
    character(len=128),             intent(in) :: chain_dir
    real(dp),     dimension(1:,1:), intent(in) :: coeff

    integer(i4b)       :: k, j, l, numband, unit
    character(len=4)   :: chain_text
    character(len=5)   :: i_text
    character(len=128) :: filename

    unit    = getlun()
    numband = size(coeff(1,:))
    
    call int2string(chain,     chain_text)
    call int2string(iteration, i_text)
    
    filename = trim(chain_dir) // '/temp_amp_c' // chain_text // '_k' // i_text // '.dat'
    open(unit,file=trim(filename), recl=1024)
    do k = 1, numband
       write(unit,*) bp(k)%label, real(coeff(:,k),sp)
    end do
    close(unit)

  end subroutine output_foregrounds

  subroutine init_fg_amps(paramfile, fg_amp)
    implicit none

    character(len=*),                      intent(in)  :: paramfile
    real(dp),         dimension(0:,1:,1:), intent(out) :: fg_amp

    integer(i4b)       :: i, j, unit
    logical(lgt)       :: exist
    real(dp)           :: scale
    character(len=512) :: filename
    character(len=2)   :: i_text

    if (num_fg_comp == 0) return

    unit   = getlun()
    fg_amp = 0.d0
    do i = 1, num_fg_comp
       if (trim(fg_components(i)%type) == 'freefree_EM') then
          fg_amp(:,:,i) = 1.d0
          cycle
       end if

       call int2string(i, i_text)
       call get_parameter(paramfile, 'INITIAL_AMPLITUDE_MAP'//i_text, par_string=filename)
       inquire(file=trim(filename), exist=exist)
       if (exist) then
          call read_map(filename, fg_amp(:,:,i))
       end if

       ! Convert from external native units to internal antenna temperature units
       scale = ant2unit(fg_components(i)%amp_unit, fg_components(i)%nu_ref, &
            & band_ref=fg_components(i)%ref_band)
       fg_amp(:,:,i) = fg_amp(:,:,i) / scale

    end do

  end subroutine init_fg_amps
  
!!$  subroutine read_temp_amp_init(filename, s_i)
!!$    implicit none
!!$
!!$    character(len=*), intent(in)    :: filename
!!$    type(genvec),     intent(inout) :: s_i
!!$
!!$    logical(lgt) :: exist
!!$    integer(i4b) :: unit
!!$    character(len=32) :: label
!!$
!!$    inquire(file=trim(filename), exist=exist)
!!$    if (exist) then
!!$       unit = getlun()
!!$       open(unit,file=trim(filename))
!!$       do i = 1, numband
!!$          read(unit,*) label, s_i%temp_amp(:,i)
!!$          do j = 1, num_fg_temp
!!$             if (fix_temp(j,i)) s_i%temp_amp(j,i) = 0.d0
!!$          end do
!!$       end do
!!$       close(unit)
!!$    end if
!!$
!!$  end subroutine read_temp_amp_init


  subroutine output_fg_overview(parfile,chain, iter, chain_dir, s_i, par_smooth)
    implicit none

    integer(i4b),                   intent(in) :: chain, iter
    character(len=128),             intent(in) :: parfile, chain_dir
    type(genvec),                   intent(in) :: s_i
    real(dp), dimension(0:,1:,1:),  intent(in) :: par_smooth

    real(dp)           :: dnu, mu, sigma, nu0, t1, t2
    integer(i4b)       :: c, i, k, j, l, p, numband, unit, n, nspline, k_min, k_max
    character(len=4)   :: chain_text
    character(len=5)   :: i_text
    character(len=128) :: filename

    real(dp), allocatable, dimension(:)           :: nu, f2
    real(dp), allocatable, dimension(:,:)         :: f, map

    real(dp),                                save :: nu_min, nu_max
    real(dp), allocatable, dimension(:,:,:), save :: mask

    integer(i4b), parameter :: nbands = 9
    real(dp), dimension(nbands) :: nuc       = [30, 44, 70, 100, 143, 217, 353, 545, 857]
    real(dp), dimension(nbands) :: bandwidth = [ 6,  9, 14,  33,  47,  72, 116, 181, 285]
    !real(dp), dimension(nbands) :: nuc       = [23, 33, 41, 61, 94, 30, 44, 70, 100, 143, 217, 353, 545, 857]
    !real(dp), dimension(nbands) :: bandwidth = [6,   7,  8, 14, 25,  6,  9, 14,  33,  47,  72, 116, 181, 285]

    unit    = getlun()
    n       = 50
    nspline = 1000
    
    if (.not. (allocated(mask))) then
       allocate(mask(0:npix-1,nmaps,2))
       call get_parameter(paramfile, 'OVERVIEW_MASK_HIGHFG', par_string=filename)
       call read_map(filename, mask(:,:,1))
       call get_parameter(paramfile, 'OVERVIEW_MASK_LOWFG',  par_string=filename)
       call read_map(filename, mask(:,:,2))
       call get_parameter(paramfile, 'OVERVIEW_FREQ_MIN',    par_dp=nu_min)
       call get_parameter(paramfile, 'OVERVIEW_FREQ_MAX',    par_dp=nu_max)
       nu_min = nu_min * 1d9
       nu_max = nu_max * 1d9
    end if

    call int2string(chain,     chain_text)
    call int2string(iter,      i_text)

    allocate(nu(n), f(n,2), f2(n), map(0:npix-1,nmaps))
    do k = 1, n
       nu(k) = nu_min * (nu_max/nu_min)**((k-1.d0)/(n-1.d0))
    end do

    ! Temperature
    if (any(s_i%fg_amp(:,1,:) /= 0.d0)) then
       
       call wall_time(t1)
       filename = trim(chain_dir) // '/overview_T_c' // chain_text // '_k' // i_text // '.dat'
       open(unit,file=trim(filename), recl=1024)
       
       do i = 1, nbands
          write(unit,*) nuc(i)-0.5*bandwidth(i), 1.d-6
          write(unit,*) nuc(i)-0.5*bandwidth(i), 1.d3
          write(unit,*) nuc(i)+0.5*bandwidth(i), 1.d3
          write(unit,*) nuc(i)+0.5*bandwidth(i), 1.d-6
          write(unit,*)
       end do
       
       j = 1
       do i = 1, num_fg_comp
          if (trim(fg_components(i)%type) == 'CO_multiline') then
             j = j + fg_components(i)%npar
             cycle
          end if
          do k = 1, n
             map = 0.d0
             do p = 0, npix-1
                if (any(mask(p,1,:)==1.d0)) then
                   map(p,1) = s_i%fg_amp(p,1,i) * get_ideal_fg_spectrum(fg_components(i), &
                        & par_smooth(p,1,j:j+fg_components(i)%npar-1), nu(k))
                end if
             end do
             mu     = sum(map(:,1)*mask(:,1,1)) / sum(mask(:,1,1))
             f(k,1) = log(max(sqrt(sum((map(:,1)*mask(:,1,1)-mu)**2) / (sum(mask(:,1,1))-1.d0)),1.d-30))
             mu     = sum(map(:,1)*mask(:,1,2)) / sum(mask(:,1,2))
             f(k,2) = log(max(sqrt(sum((map(:,1)*mask(:,1,2)-mu)**2) / (sum(mask(:,1,2))-1.d0)),1.d-30))
          end do
          
          k_min = 1
          do while (f(k_min,1) == log(1.d-30) .or. f(k_min,2) == log(1.d-30) .and. k_min < n)
             k_min = k_min+1
          end do
          
          k_max = n
          do while (f(k_max,1) == log(1.d-30) .or. f(k_max,2) == log(1.d-30) .and. k_max > k_min)
             k_max = k_max-1
          end do
          
          call spline(nu(k_min:k_max), f(k_min:k_max,1), 1.d30, 1.d30, f2(k_min:k_max))
          do k = 1, nspline
             nu0 = nu_min + (nu_max-nu_min) / (nspline-1.d0) * (k-1.d0)
             if (nu0 < nu(k_min) .or. nu0 > nu(k_max)) cycle
             write(unit,*) nu0*1d-9, exp(splint(nu(k_min:k_max), f(k_min:k_max,1), f2(k_min:k_max), nu0))
          end do
          call spline(nu(k_min:k_max), f(k_min:k_max,2), 1.d30, 1.d30, f2(k_min:k_max))
          do k = nspline, 1, -1
             nu0 = nu_min + (nu_max-nu_min) / (nspline-1.d0) * (k-1.d0)
             if (nu0 < nu(k_min) .or. nu0 > nu(k_max)) cycle
             write(unit,*) nu0*1d-9, exp(splint(nu(k_min:k_max), f(k_min:k_max,2), f2(k_min:k_max), nu0))
          end do
          write(unit,*)
          
          j = j + fg_components(i)%npar
          
       end do
       close(unit)
    end if

    ! Polarization
    if (polarization .and. any(s_i%fg_amp(:,2:3,:) /= 0.d0)) then
       
       call wall_time(t1)
       filename = trim(chain_dir) // '/overview_P_c' // chain_text // '_k' // i_text // '.dat'
       open(unit,file=trim(filename), recl=1024)
       
       do i = 1, nbands
          write(unit,*) nuc(i)-0.5*bandwidth(i), 1.d-6
          write(unit,*) nuc(i)-0.5*bandwidth(i), 1.d3
          write(unit,*) nuc(i)+0.5*bandwidth(i), 1.d3
          write(unit,*) nuc(i)+0.5*bandwidth(i), 1.d-6
          write(unit,*)
       end do
       
       j = 1
       do i = 1, num_fg_comp
          if (trim(fg_components(i)%type) == 'CO_multiline') then
             j = j + fg_components(i)%npar
             cycle
          end if
          do k = 1, n
             map = 0.d0
             do p = 0, npix-1
                if (any(mask(p,2:3,:)==1.d0)) then
                   map(p,1) = sqrt(s_i%fg_amp(p,2,i)**2+s_i%fg_amp(p,3,i)**2) * &
                        & get_ideal_fg_spectrum(fg_components(i), &
                        & par_smooth(p,1,j:j+fg_components(i)%npar-1), nu(k))
                end if
             end do
             mu     = sum(map(:,1)*mask(:,2,1)) / sum(mask(:,2,1))
             f(k,1) = log(max(sqrt(sum((map(:,1)*mask(:,2,1)-mu)**2) / (sum(mask(:,2,1))-1.d0)),1.d-30))
             mu     = sum(map(:,1)*mask(:,2,2)) / sum(mask(:,2,2))
             f(k,2) = log(max(sqrt(sum((map(:,1)*mask(:,2,2)-mu)**2) / (sum(mask(:,2,2))-1.d0)),1.d-30))
          end do
          
          k_min = 1
          do while (f(k_min,1) == log(1.d-30) .or. f(k_min,2) == log(1.d-30) .and. k_min < n)
             k_min = k_min+1
          end do
          
          k_max = n
          do while (f(k_max,1) == log(1.d-30) .or. f(k_max,2) == log(1.d-30) .and. k_max > k_min)
             k_max = k_max-1
          end do
          
          call spline(nu(k_min:k_max), f(k_min:k_max,1), 1.d30, 1.d30, f2(k_min:k_max))
          do k = 1, nspline
             nu0 = nu_min + (nu_max-nu_min) / (nspline-1.d0) * (k-1.d0)
             if (nu0 < nu(k_min) .or. nu0 > nu(k_max)) cycle
             write(unit,*) nu0*1d-9, exp(splint(nu(k_min:k_max), f(k_min:k_max,1), f2(k_min:k_max), nu0))
          end do
          call spline(nu(k_min:k_max), f(k_min:k_max,2), 1.d30, 1.d30, f2(k_min:k_max))
          do k = nspline, 1, -1
             nu0 = nu_min + (nu_max-nu_min) / (nspline-1.d0) * (k-1.d0)
             if (nu0 < nu(k_min) .or. nu0 > nu(k_max)) cycle
             write(unit,*) nu0*1d-9, exp(splint(nu(k_min:k_max), f(k_min:k_max,2), f2(k_min:k_max), nu0))
          end do
          write(unit,*)
          
          j = j + fg_components(i)%npar
          
       end do
       close(unit)
    end if
    call wall_time(t2)
    write(*,*) 'Wall time overview = ', t2-t1
    
    deallocate(nu, f, f2, map)


  end subroutine output_fg_overview


end program commander
