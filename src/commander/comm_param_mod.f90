!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Commander parameter module !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module comm_param_mod
  use comm_utils
  use comm_status_mod
  use hashtbl
  implicit none

  ! Note: This module reads in the Commander parameter file as the first operation
  !       in the program. This is primarily intended to avoid crashes after hours
  !       of running because of user errors; catch these early, and report back
  !       problems. Then, copy parameters over to module structures if convenient
  !       at a later stage. 

  integer(i4b), parameter, private :: MAXPAR       = 10
  integer(i4b), parameter, private :: MAXAUXPAR    = 10
  integer(i4b), parameter, private :: MAXSAMPGROUP = 100
  type(status_file)                :: status
  
  type comm_params

     ! MPI info
     integer(i4b) :: myid, numprocs, root = 0
     integer(i4b) :: myid_chain, numprocs_chain, comm_chain, mychain
     integer(i4b) :: myid_shared, comm_shared, myid_inter, comm_inter
     integer(i4b), dimension(MPI_STATUS_SIZE)          :: status

     ! Global parameters
     character(len=24)  :: operation
     logical(lgt)       :: resamp_CMB
     integer(i4b)       :: first_samp_resamp, last_samp_resamp, numsamp_per_resamp
     integer(i4b)       :: verbosity, base_seed, base_seed_noise, numchain, num_smooth_scales
     integer(i4b)       :: num_gibbs_iter
     character(len=512) :: chain_status, init_chain_prefix
     real(dp)           :: T_CMB
     character(len=512) :: MJysr_convention
     character(len=512) :: fft_magic_number_file
     logical(lgt)       :: only_pol, optimize_alm
     logical(lgt)       :: enable_TOD_analysis
     integer(i4b)       :: tod_freq
     integer(i4b)       :: nsamp_alm, nside_chisq_lowres, prior_fwhm, burnin !alm sampler
     integer(i4b)       :: resamp_hard_gain_prior_nth_iter
     integer(i4b)       :: output_4D_map_nth_iter
     logical(lgt)       :: include_tod_zodi
     real(dp),           allocatable, dimension(:)     :: fwhm_smooth
     real(dp),           allocatable, dimension(:)     :: fwhm_postproc_smooth
     integer(i4b),       allocatable, dimension(:)     :: lmax_smooth
     integer(i4b),       allocatable, dimension(:)     :: nside_smooth
     character(len=512), allocatable, dimension(:)     :: pixwin_smooth

     ! Output parameters
     character(len=512) :: outdir
     integer(i4b)       :: nside_chisq, nmaps_chisq
     logical(lgt)       :: pol_chisq, output_mixmat, output_residuals, output_chisq, output_cg_eigenvals
     integer(i4b)       :: output_cg_freq
     logical(lgt)       :: output_input_model, ignore_gain_bp, output_debug_seds, output_sig_per_band
     logical(lgt)       :: sample_signal_amplitudes, sample_specind
     
     ! Numerical parameters
     character(len=512) :: cg_conv_crit, cg_precond
     integer(i4b)       :: cg_lmax_precond, cg_maxiter, cg_num_samp_groups, cg_num_user_samp_groups, cg_miniter, cg_check_conv_freq
     logical(lgt)       :: cg_init_zero, set_noise_to_mean
     real(dp)           :: cg_tol
     integer(i4b)       :: num_bp_prop
     character(len=512), dimension(MAXSAMPGROUP) :: cg_samp_group
     character(len=512), dimension(MAXSAMPGROUP) :: cg_samp_group_mask

     ! Data parameters
     integer(i4b)       :: numband
     character(len=512) :: datadir, ds_sourcemask, ds_procmask
     logical(lgt),       allocatable, dimension(:)   :: ds_active
     integer(i4b),       allocatable, dimension(:)   :: ds_period
     logical(lgt),       allocatable, dimension(:)   :: ds_polarization
     integer(i4b),       allocatable, dimension(:)   :: ds_nside
     integer(i4b),       allocatable, dimension(:)   :: ds_lmax
     character(len=512), allocatable, dimension(:)   :: ds_label
     character(len=512), allocatable, dimension(:)   :: ds_unit
     character(len=512), allocatable, dimension(:)   :: ds_noise_format
     character(len=512), allocatable, dimension(:)   :: ds_mapfile
     character(len=512), allocatable, dimension(:)   :: ds_noisefile
     character(len=512), allocatable, dimension(:)   :: ds_regnoise
     character(len=512), allocatable, dimension(:,:) :: ds_noise_rms_smooth
     real(dp),           allocatable, dimension(:)   :: ds_noise_uni_fsky
     character(len=512), allocatable, dimension(:)   :: ds_maskfile
     character(len=512), allocatable, dimension(:)   :: ds_maskfile_calib
     character(len=512), allocatable, dimension(:)   :: ds_beamtype
     character(len=512), allocatable, dimension(:)   :: ds_blfile
     character(len=512), allocatable, dimension(:)   :: ds_btheta_file
     character(len=512), allocatable, dimension(:)   :: ds_pixwin
     logical(lgt),       allocatable, dimension(:)   :: ds_samp_noiseamp
     character(len=512), allocatable, dimension(:)   :: ds_bptype
     character(len=512), allocatable, dimension(:)   :: ds_bpfile
     character(len=512), allocatable, dimension(:)   :: ds_bpmodel
     real(dp),           allocatable, dimension(:)   :: ds_nu_c
     logical(lgt),       allocatable, dimension(:)   :: ds_sample_gain
     real(dp),           allocatable, dimension(:,:) :: ds_gain_prior
     character(len=512), allocatable, dimension(:)   :: ds_gain_calib_comp
     integer(i4b),       allocatable, dimension(:)   :: ds_gain_lmin
     integer(i4b),       allocatable, dimension(:)   :: ds_gain_lmax
     character(len=512), allocatable, dimension(:)   :: ds_gain_apodmask
     character(len=512), allocatable, dimension(:)   :: ds_gain_fwhm
     real(dp),           allocatable, dimension(:,:) :: ds_defaults
     character(len=512), allocatable, dimension(:)   :: ds_component_sensitivity
     !TOD data parameters
     character(len=512), allocatable, dimension(:)   :: ds_tod_type
     character(len=512), allocatable, dimension(:)   :: ds_tod_procmask1
     character(len=512), allocatable, dimension(:)   :: ds_tod_procmask2
     character(len=512), allocatable, dimension(:)   :: ds_tod_filelist
     character(len=512), allocatable, dimension(:)   :: ds_tod_instfile
     character(len=512), allocatable, dimension(:)   :: ds_tod_dets
     character(len=512), allocatable, dimension(:)   :: ds_tod_bp_init
     character(len=512), allocatable, dimension(:)   :: ds_tod_initHDF
     integer(i4b),       allocatable, dimension(:,:) :: ds_tod_scanrange
     integer(i4b),       allocatable, dimension(:)   :: ds_tod_tot_numscan
     integer(i4b),       allocatable, dimension(:)   :: ds_tod_flag
     logical(lgt),       allocatable, dimension(:)   :: ds_tod_orb_abscal

     ! Component parameters
     character(len=512) :: cs_inst_parfile
     character(len=512) :: cs_init_inst_hdf
     integer(i4b)       :: cs_ncomp, cs_ncomp_tot
     real(dp)           :: cmb_dipole_prior(3)
     character(len=512) :: cmb_dipole_prior_mask
     logical(lgt),       allocatable, dimension(:)     :: cs_include
     character(len=512), allocatable, dimension(:)     :: cs_initHDF
     character(len=512), allocatable, dimension(:)     :: cs_label
     character(len=512), allocatable, dimension(:)     :: cs_type
     character(len=512), allocatable, dimension(:)     :: cs_class
     logical(lgt),       allocatable, dimension(:)     :: cs_polarization
     real(dp),           allocatable, dimension(:)     :: cs_cg_scale
     !integer(i4b),       allocatable, dimension(:)     :: cs_cg_samp_group
     integer(i4b),       allocatable, dimension(:)     :: cs_nside
     integer(i4b),       allocatable, dimension(:,:)   :: cs_poltype
     integer(i4b),       allocatable, dimension(:)     :: cs_lmax_amp
     integer(i4b),       allocatable, dimension(:)     :: cs_l_apod
     integer(i4b),       allocatable, dimension(:)     :: cs_lmax_ind
     character(len=512), allocatable, dimension(:)     :: cs_unit
     real(dp),           allocatable, dimension(:,:)   :: cs_nu_ref
     character(len=512), allocatable, dimension(:)     :: cs_band_ref
     real(dp),           allocatable, dimension(:)     :: cs_fwhm
     character(len=512), allocatable, dimension(:)     :: cs_cltype
     character(len=512), allocatable, dimension(:)     :: cs_clfile
     character(len=512), allocatable, dimension(:)     :: cs_binfile
     integer(i4b),       allocatable, dimension(:)     :: cs_lpivot
     character(len=512), allocatable, dimension(:)     :: cs_mask
     character(len=512), allocatable, dimension(:)     :: cs_mono_prior
     real(dp),           allocatable, dimension(:)     :: cs_latmask
     character(len=512), allocatable, dimension(:)     :: cs_indmask
     character(len=512), allocatable, dimension(:)     :: cs_defmask
     real(dp),           allocatable, dimension(:,:)   :: cs_cl_prior
     real(dp),           allocatable, dimension(:,:)   :: cs_cl_amp_def
     real(dp),           allocatable, dimension(:,:)   :: cs_cl_beta_def
     integer(i4b),       allocatable, dimension(:)     :: cs_cl_poltype
     logical(lgt),       allocatable, dimension(:)     :: cs_output_EB
     character(len=512), allocatable, dimension(:)     :: cs_input_amp
     character(len=512), allocatable, dimension(:)     :: cs_prior_amp
     character(len=512), allocatable, dimension(:,:)   :: cs_input_ind
     character(len=512), allocatable, dimension(:,:)   :: cs_SED_template
     real(dp),           allocatable, dimension(:,:)   :: cs_theta_def
     integer(i4b),       allocatable, dimension(:,:)   :: cs_smooth_scale
     real(dp),           allocatable, dimension(:,:,:) :: cs_p_gauss
     real(dp),           allocatable, dimension(:,:,:) :: cs_p_uni
     character(len=512), allocatable, dimension(:)     :: cs_catalog
     character(len=512), allocatable, dimension(:)     :: cs_ptsrc_template
     real(dp),           allocatable, dimension(:,:)   :: cs_nu_min
     real(dp),           allocatable, dimension(:,:)   :: cs_nu_max
     logical(lgt),       allocatable, dimension(:)     :: cs_burn_in
     logical(lgt),       allocatable, dimension(:)     :: cs_output_ptsrc_beam
     logical(lgt),       allocatable, dimension(:)     :: cs_apply_pos_prior
     real(dp),           allocatable, dimension(:)     :: cs_min_src_dist
     real(dp),           allocatable, dimension(:)     :: cs_amp_rms_scale
     real(dp),           allocatable, dimension(:,:)   :: cs_auxpar
     logical(lgt),       allocatable, dimension(:)     :: cs_apply_jeffreys
  end type comm_params

contains

  ! ********************************************************
  !                     Driver routines
  ! ********************************************************
  subroutine read_comm_params(cpar)
    implicit none
    type(hash_tbl_sll) :: htable
    type(comm_params), intent(inout) :: cpar
    
    integer(i4b)       :: paramfile_len, ierr, i
    character(len=512) :: paramfile
    character(len=512), allocatable, dimension(:) :: paramfile_cache

    call getarg(1, paramfile)
    ! read parameter file once, save to ascii array
    ! Need to know how long the file is to allocate ascii array
    if (cpar%myid == cpar%root) then
       call get_file_length(paramfile,paramfile_len)
    end if
    
    call mpi_bcast(paramfile_len, 1, MPI_INTEGER, cpar%root, MPI_COMM_WORLD, ierr)
    allocate(paramfile_cache(paramfile_len))

    if (cpar%myid == cpar%root) then
       call read_paramfile_to_ascii(paramfile,paramfile_cache)
    end if
    do i=1,paramfile_len
       call mpi_bcast(paramfile_cache(i), 512, MPI_CHAR, cpar%root, MPI_COMM_WORLD, ierr)
    end do
    
    !Initialize a hash table
    call init_hash_tbl_sll(htable,tbl_len=10*paramfile_len)
    ! Put the parameter file into the hash table
    call put_ascii_into_hashtable(paramfile_cache,htable)
    deallocate(paramfile_cache)

    ! Read parameters from the hash table
    call read_global_params_hash(htable,cpar)
    call read_data_params_hash(htable,cpar)
    call read_component_params_hash(htable,cpar)

    ! Override parameter choices if user if RESAMPLE_CMB = .true.
    if (cpar%resamp_CMB) then
       cpar%operation           = 'sample'     ! Force sampling
       cpar%cg_precond          = 'diagonal'   ! Use diagonal precond to fill in the mask
       cpar%enable_TOD_analysis = .false.      ! Disable TOD analysis
       cpar%sample_specind      = .false.      ! Disable non-linear parameter fitting
       !cpar%num_gibbs_iter      = (cpar%last_samp_resamp-cpar%first_samp_resamp+1)*cpar%numsamp_per_resamp
!!$       do i = 1, cpar%cs_ncomp_tot
!!$          if (trim(cpar%cs_type(i)) /= 'cmb') then
!!$             cpar%cs_cg_samp_group(i) = 0   ! Disable all other components than CMB from CG search
!!$          end if
!!$       end do
    end if
    
    !Deallocate hash table
    call free_hash_tbl_sll(htable)
  end subroutine read_comm_params

  subroutine initialize_mpi_struct(cpar, handle, handle_noise, reinit_rng)
    implicit none
    type(comm_params), intent(inout) :: cpar
    type(planck_rng),  intent(out)   :: handle, handle_noise
    integer(i4b),      intent(in), optional :: reinit_rng

    integer(i4b) :: i, j, m, n, ierr, mpistat(MPI_STATUS_SIZE)
    integer(i4b), allocatable, dimension(:,:) :: ind
    
    ! !the following commented lines are moved to commander.f90 in order to initialize
    ! !the mpi structure so that one only reads parameter file with one processor
    !call mpi_init(ierr)
    !call mpi_comm_rank(MPI_COMM_WORLD, cpar%myid, ierr)
    !call mpi_comm_size(MPI_COMM_WORLD, cpar%numprocs, ierr)
    !cpar%root = 0
    if (.not. (present(reinit_rng))) then
       cpar%numchain = min(cpar%numchain, cpar%numprocs)

       allocate(ind(0:cpar%numprocs-1,2))
       n = 0
       do i = 1, cpar%numchain
          m = cpar%numprocs / cpar%numchain
          if ((cpar%numprocs-(cpar%numprocs/cpar%numchain)*cpar%numchain) >= i) m = m+1
          ind(n:n+m-1,1) = i
          do j = 0, m-1
             ind(n+j,2) = j
          end do
          n = n+m
       end do

       cpar%mychain    = ind(cpar%myid,1)
       cpar%myid_chain = ind(cpar%myid,2)

       call mpi_comm_split(MPI_COMM_WORLD, cpar%mychain, cpar%myid_chain, cpar%comm_chain,  ierr) 
       call mpi_comm_size(cpar%comm_chain, cpar%numprocs_chain, ierr)

       !Communicators for shared memory access
       call mpi_comm_split_type(cpar%comm_chain, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, cpar%comm_shared, ierr) 
       call mpi_comm_rank(cpar%comm_shared, cpar%myid_shared, ierr)
       call mpi_comm_split(cpar%comm_chain, cpar%myid_shared, 0, cpar%comm_inter, ierr)
       call mpi_comm_rank(cpar%comm_inter, cpar%myid_inter, ierr)

       deallocate(ind)
    end if

    ! Initialize random number generator
    if (cpar%myid == cpar%root) then
       call rand_init(handle, cpar%base_seed)
       if (present(reinit_rng)) then
          do i = 1, reinit_rng
             j = nint(rand_uni(handle)*1000000.d0)
             call rand_init(handle, j)
          end do
       end if
       do i = 1, cpar%numprocs-1
          j = nint(rand_uni(handle)*1000000.d0)
          call mpi_send(j, 1, MPI_INTEGER, i, 98, MPI_COMM_WORLD, ierr)
       end do
       call rand_init(handle_noise, cpar%base_seed_noise)
       do i = 1, cpar%numprocs-1
          j = nint(rand_uni(handle_noise)*1000000.d0)
          call mpi_send(j, 1, MPI_INTEGER, i, 98, MPI_COMM_WORLD, ierr)
       end do
    else 
       call mpi_recv(j, 1, MPI_INTEGER, cpar%root, 98, MPI_COMM_WORLD, mpistat, ierr)
       call rand_init(handle, j)
       call mpi_recv(j, 1, MPI_INTEGER, cpar%root, 98, MPI_COMM_WORLD, mpistat, ierr)
       call rand_init(handle_noise, j)
    end if

  end subroutine initialize_mpi_struct
  
  ! ********************************************************
  !              Specialized routines; one per module
  ! ********************************************************

  subroutine read_global_params_hash(htbl, cpar)
    implicit none

    type(hash_tbl_sll), intent(in) :: htbl
    type(comm_params),  intent(inout) :: cpar

    integer(i4b)     :: i
    character(len=2) :: itext
    
    
    call get_parameter_hashtable(htbl, 'VERBOSITY',                par_int=cpar%verbosity)
    call get_parameter_hashtable(htbl, 'OPERATION',                par_string=cpar%operation)
    call get_parameter_hashtable(htbl, 'RESAMPLE_CMB',             par_lgt=cpar%resamp_CMB)

    call get_parameter_hashtable(htbl, 'BASE_SEED',                par_int=cpar%base_seed)
    !call get_parameter_hashtable(htbl, 'BASE_SEED_NOISE',          par_int=cpar%base_seed_noise)
    cpar%base_seed_noise = 0  ! Not currently in use
    call get_parameter_hashtable(htbl, 'NUMCHAIN',                 par_int=cpar%numchain)
    call get_parameter_hashtable(htbl, 'NUM_GIBBS_ITER',           par_int=cpar%num_gibbs_iter)
    call get_parameter_hashtable(htbl, 'CHAIN_STATUS',             par_string=cpar%chain_status)
    call get_parameter_hashtable(htbl, 'INIT_CHAIN',               par_string=cpar%init_chain_prefix)
    call get_parameter_hashtable(htbl, 'SAMPLE_ONLY_POLARIZATION', par_lgt=cpar%only_pol)

    call get_parameter_hashtable(htbl, 'CG_CONVERGENCE_CRITERION', par_string=cpar%cg_conv_crit)
    call get_parameter_hashtable(htbl, 'CG_PRECOND_TYPE',          par_string=cpar%cg_precond)
    call get_parameter_hashtable(htbl, 'CG_LMAX_PRECOND',          par_int=cpar%cg_lmax_precond)
    call get_parameter_hashtable(htbl, 'CG_MINITER',               par_int=cpar%cg_miniter)
    call get_parameter_hashtable(htbl, 'CG_MAXITER',               par_int=cpar%cg_maxiter)
    call get_parameter_hashtable(htbl, 'CG_TOLERANCE',             par_dp=cpar%cg_tol)
    call get_parameter_hashtable(htbl, 'CG_CONV_CHECK_FREQUENCY',  par_int=cpar%cg_check_conv_freq)
    call get_parameter_hashtable(htbl, 'CG_INIT_AMPS_ON_ZERO',     par_lgt=cpar%cg_init_zero)
    call get_parameter_hashtable(htbl, 'SET_ALL_NOISE_MAPS_TO_MEAN',     par_lgt=cpar%set_noise_to_mean)

    call get_parameter_hashtable(htbl, 'T_CMB',                    par_dp=cpar%T_cmb)
    call get_parameter_hashtable(htbl, 'MJYSR_CONVENTION',         par_string=cpar%MJysr_convention)

    call get_parameter_hashtable(htbl, 'OUTPUT_DIRECTORY',         par_string=cpar%outdir)

    call get_parameter_hashtable(htbl, 'NSIDE_CHISQ',              par_int=cpar%nside_chisq)
    call get_parameter_hashtable(htbl, 'POLARIZATION_CHISQ',       par_lgt=cpar%pol_chisq)
    cpar%nmaps_chisq = 1; if (cpar%pol_chisq) cpar%nmaps_chisq = 3

    call get_parameter_hashtable(htbl, 'OUTPUT_MIXING_MATRIX',     par_lgt=cpar%output_mixmat)
    call get_parameter_hashtable(htbl, 'OUTPUT_RESIDUAL_MAPS',     par_lgt=cpar%output_residuals)
    call get_parameter_hashtable(htbl, 'OUTPUT_CHISQ_MAP',         par_lgt=cpar%output_chisq)
    call get_parameter_hashtable(htbl, 'OUTPUT_EVERY_NTH_CG_ITERATION', par_int=cpar%output_cg_freq)
    call get_parameter_hashtable(htbl, 'OUTPUT_CG_PRECOND_EIGENVALS', par_lgt=cpar%output_cg_eigenvals)
    call get_parameter_hashtable(htbl, 'OUTPUT_INPUT_MODEL',       par_lgt=cpar%output_input_model)
    call get_parameter_hashtable(htbl, 'IGNORE_GAIN_AND_BANDPASS_CORR', par_lgt=cpar%ignore_gain_bp)
    call get_parameter_hashtable(htbl, 'OUTPUT_DEBUG_SEDS',        par_lgt=cpar%output_debug_seds)
    call get_parameter_hashtable(htbl, 'OUTPUT_SIGNALS_PER_BAND',  par_lgt=cpar%output_sig_per_band)

    call get_parameter_hashtable(htbl, 'SAMPLE_SIGNAL_AMPLITUDES', par_lgt=cpar%sample_signal_amplitudes)
    call get_parameter_hashtable(htbl, 'SAMPLE_SPECTRAL_INDICES',  par_lgt=cpar%sample_specind)

    call get_parameter_hashtable(htbl, 'NUM_SMOOTHING_SCALES',     par_int=cpar%num_smooth_scales)

    call get_parameter_hashtable(htbl, 'ENABLE_TOD_ANALYSIS',      par_lgt=cpar%enable_TOD_analysis)

    call get_parameter_hashtable(htbl, 'NUMITER_RESAMPLE_HARD_GAIN_PRIORS', par_int=cpar%resamp_hard_gain_prior_nth_iter)

    if (cpar%enable_TOD_analysis) then
       call get_parameter_hashtable(htbl, 'FFTW3_MAGIC_NUMBERS',   par_string=cpar%fft_magic_number_file)
       call get_parameter_hashtable(htbl, 'TOD_NUM_BP_PROPOSALS_PER_ITER', par_int=cpar%num_bp_prop)
       call get_parameter_hashtable(htbl, 'NUM_GIBBS_STEPS_PER_TOD_SAMPLE', par_int=cpar%tod_freq)
       call get_parameter_hashtable(htbl, 'TOD_OUTPUT_4D_MAP_EVERY_NTH_ITER', par_int=cpar%output_4D_map_nth_iter)
       call get_parameter_hashtable(htbl, 'TOD_INCLUDE_ZODI',      par_lgt=cpar%include_TOD_zodi)
    end if

    if (cpar%resamp_CMB) then
       call get_parameter_hashtable(htbl, 'FIRST_SAMPLE_FOR_CMB_RESAMP',   par_int=cpar%first_samp_resamp)
       call get_parameter_hashtable(htbl, 'LAST_SAMPLE_FOR_CMB_RESAMP',    par_int=cpar%last_samp_resamp)
       call get_parameter_hashtable(htbl, 'NUM_SUBSAMP_PER_MAIN_SAMPLE',   par_int=cpar%numsamp_per_resamp)
    end if

    allocate(cpar%fwhm_smooth(cpar%num_smooth_scales))
    allocate(cpar%fwhm_postproc_smooth(cpar%num_smooth_scales))
    allocate(cpar%lmax_smooth(cpar%num_smooth_scales))
    allocate(cpar%nside_smooth(cpar%num_smooth_scales))
    allocate(cpar%pixwin_smooth(cpar%num_smooth_scales))
    do i = 1, cpar%num_smooth_scales
       call int2string(i, itext)
       call get_parameter_hashtable(htbl, 'SMOOTHING_SCALE_FWHM'//itext, len_itext=len(trim(itext)), par_dp=cpar%fwhm_smooth(i))
       call get_parameter_hashtable(htbl, 'SMOOTHING_SCALE_FWHM_POSTPROC'//itext, &
            & len_itext=len(trim(itext)), par_dp=cpar%fwhm_postproc_smooth(i))
       call get_parameter_hashtable(htbl, 'SMOOTHING_SCALE_LMAX'//itext, len_itext=len(trim(itext)), par_int=cpar%lmax_smooth(i))
       call get_parameter_hashtable(htbl, 'SMOOTHING_SCALE_NSIDE'//itext, len_itext=len(trim(itext)), par_int=cpar%nside_smooth(i))
       call get_parameter_hashtable(htbl, 'SMOOTHING_SCALE_PIXWIN'//itext, &
            & len_itext=len(trim(itext)), par_string=cpar%pixwin_smooth(i))
    end do
    
    call get_parameter_hashtable(htbl, 'ALMSAMP_NSAMP_ALM',          par_int=cpar%nsamp_alm)
    call get_parameter_hashtable(htbl, 'ALMSAMP_BURN_IN',            par_int=cpar%burnin)
    call get_parameter_hashtable(htbl, 'ALMSAMP_PRIOR_FWHM',         par_int=cpar%prior_fwhm)
    call get_parameter_hashtable(htbl, 'ALMSAMP_NSIDE_CHISQ_LOWRES', par_int=cpar%nside_chisq_lowres)
    call get_parameter_hashtable(htbl, 'ALMSAMP_OPTIMIZE_ALM',       par_lgt=cpar%optimize_alm)

  end subroutine read_global_params_hash


  subroutine read_data_params_hash(htbl, cpar)
    implicit none

    type(hash_tbl_sll), intent(in) :: htbl
    type(comm_params),  intent(inout) :: cpar

    integer(i4b)     :: i, j, n,len_itext
    character(len=3) :: itext
    character(len=2) :: jtext

    len_itext=len(trim(itext))
    call get_parameter_hashtable(htbl, 'NUMBAND',             par_int=cpar%numband)
    call get_parameter_hashtable(htbl, 'DATA_DIRECTORY',      par_string=cpar%datadir)
    call get_parameter_hashtable(htbl, 'SOURCE_MASKFILE',     par_string=cpar%ds_sourcemask)
    call get_parameter_hashtable(htbl, 'PROCESSING_MASKFILE', par_string=cpar%ds_procmask)

    n = cpar%numband
    allocate(cpar%ds_active(n), cpar%ds_label(n))
    allocate(cpar%ds_polarization(n), cpar%ds_nside(n), cpar%ds_lmax(n))
    allocate(cpar%ds_unit(n), cpar%ds_noise_format(n), cpar%ds_mapfile(n))
    allocate(cpar%ds_noisefile(n), cpar%ds_maskfile(n), cpar%ds_maskfile_calib(n))
    allocate(cpar%ds_regnoise(n))
    allocate(cpar%ds_noise_rms_smooth(n,cpar%num_smooth_scales))
    allocate(cpar%ds_samp_noiseamp(n), cpar%ds_noise_uni_fsky(n))
    allocate(cpar%ds_bptype(n), cpar%ds_nu_c(n), cpar%ds_bpfile(n), cpar%ds_bpmodel(n))
    allocate(cpar%ds_period(n), cpar%ds_beamtype(n), cpar%ds_blfile(n))
    allocate(cpar%ds_pixwin(n), cpar%ds_btheta_file(n))
    allocate(cpar%ds_sample_gain(n), cpar%ds_gain_prior(n,2), cpar%ds_gain_calib_comp(n), cpar%ds_gain_lmax(n))
    allocate(cpar%ds_gain_lmin(n), cpar%ds_gain_apodmask(n), cpar%ds_gain_fwhm(n))
    allocate(cpar%ds_defaults(n,2))
    allocate(cpar%ds_component_sensitivity(n))
    allocate(cpar%ds_tod_type(n), cpar%ds_tod_filelist(n), cpar%ds_tod_initHDF(n))
    allocate(cpar%ds_tod_procmask1(n), cpar%ds_tod_procmask2(n), cpar%ds_tod_bp_init(n))
    allocate(cpar%ds_tod_instfile(n), cpar%ds_tod_dets(n), cpar%ds_tod_scanrange(n,2))
    allocate(cpar%ds_tod_tot_numscan(n), cpar%ds_tod_flag(n), cpar%ds_tod_orb_abscal(n))

    do i = 1, n
       call int2string(i, itext)
       call get_parameter_hashtable(htbl, 'INCLUDE_BAND'//itext, len_itext=len_itext, par_lgt=cpar%ds_active(i))
       if (.not. cpar%ds_active(i)) cycle
       call get_parameter_hashtable(htbl, 'BAND_OBS_PERIOD'//itext, len_itext=len_itext, par_int=cpar%ds_period(i))
       call get_parameter_hashtable(htbl, 'BAND_LABEL'//itext, len_itext=len_itext, par_string=cpar%ds_label(i))
       call get_parameter_hashtable(htbl, 'BAND_POLARIZATION'//itext, len_itext=len_itext, par_lgt=cpar%ds_polarization(i))
       call get_parameter_hashtable(htbl, 'BAND_NSIDE'//itext, len_itext=len_itext, par_int=cpar%ds_nside(i))
       call get_parameter_hashtable(htbl, 'BAND_LMAX'//itext, len_itext=len_itext, par_int=cpar%ds_lmax(i))
       call get_parameter_hashtable(htbl, 'BAND_UNIT'//itext, len_itext=len_itext, par_string=cpar%ds_unit(i))
       call get_parameter_hashtable(htbl, 'BAND_NOISE_FORMAT'//itext, len_itext=len_itext, par_string=cpar%ds_noise_format(i))
       call get_parameter_hashtable(htbl, 'BAND_MAPFILE'//itext, len_itext=len_itext, par_string=cpar%ds_mapfile(i))
       call get_parameter_hashtable(htbl, 'BAND_NOISEFILE'//itext, len_itext=len_itext, par_string=cpar%ds_noisefile(i))
       call get_parameter_hashtable(htbl, 'BAND_REG_NOISEFILE'//itext, len_itext=len_itext, par_string=cpar%ds_regnoise(i))
       call get_parameter_hashtable(htbl, 'BAND_NOISE_UNIFORMIZE_FSKY'//itext, len_itext=len_itext, &
            & par_dp=cpar%ds_noise_uni_fsky(i))
       call get_parameter_hashtable(htbl, 'BAND_MASKFILE'//itext, len_itext=len_itext,        par_string=cpar%ds_maskfile(i))
       call get_parameter_hashtable(htbl, 'BAND_MASKFILE_CALIB'//itext, len_itext=len_itext,  par_string=cpar%ds_maskfile_calib(i))
       call get_parameter_hashtable(htbl, 'BAND_BEAMTYPE'//itext, len_itext=len_itext,        par_string=cpar%ds_beamtype(i))
       call get_parameter_hashtable(htbl, 'BAND_BEAM_B_L_FILE'//itext, len_itext=len_itext,   par_string=cpar%ds_blfile(i))
       call get_parameter_hashtable(htbl, 'BAND_BEAM_B_PTSRC_FILE'//itext, len_itext=len_itext, par_string=cpar%ds_btheta_file(i))
       call get_parameter_hashtable(htbl, 'BAND_PIXEL_WINDOW'//itext, len_itext=len_itext,    par_string=cpar%ds_pixwin(i))
       call get_parameter_hashtable(htbl, 'BAND_SAMP_NOISE_AMP'//itext, len_itext=len_itext,  par_lgt=cpar%ds_samp_noiseamp(i))
       call get_parameter_hashtable(htbl, 'BAND_BANDPASS_TYPE'//itext, len_itext=len_itext,   par_string=cpar%ds_bptype(i))
       call get_parameter_hashtable(htbl, 'BAND_NOMINAL_FREQ'//itext, len_itext=len_itext,    par_dp=cpar%ds_nu_c(i))
       call get_parameter_hashtable(htbl, 'BAND_BANDPASSFILE'//itext, len_itext=len_itext,    par_string=cpar%ds_bpfile(i))
       call get_parameter_hashtable(htbl, 'BAND_BANDPASS_MODEL'//itext, len_itext=len_itext,  par_string=cpar%ds_bpmodel(i))
       call get_parameter_hashtable(htbl, 'BAND_SAMP_GAIN'//itext, len_itext=len_itext,       par_lgt=cpar%ds_sample_gain(i))
       call get_parameter_hashtable(htbl, 'BAND_GAIN_PRIOR_MEAN'//itext, len_itext=len_itext, par_dp=cpar%ds_gain_prior(i,1))
       call get_parameter_hashtable(htbl, 'BAND_GAIN_PRIOR_RMS'//itext, len_itext=len_itext, par_dp=cpar%ds_gain_prior(i,2))
       call get_parameter_hashtable(htbl, 'BAND_GAIN_CALIB_COMP'//itext, len_itext=len_itext, par_string=cpar%ds_gain_calib_comp(i))
       call get_parameter_hashtable(htbl, 'BAND_GAIN_LMAX'//itext, len_itext=len_itext,       par_int=cpar%ds_gain_lmax(i))
       call get_parameter_hashtable(htbl, 'BAND_GAIN_LMIN'//itext, len_itext=len_itext,       par_int=cpar%ds_gain_lmin(i))
       call get_parameter_hashtable(htbl, 'BAND_GAIN_APOD_MASK'//itext, len_itext=len_itext,  par_string=cpar%ds_gain_apodmask(i))
       call get_parameter_hashtable(htbl, 'BAND_GAIN_APOD_FWHM'//itext, len_itext=len_itext,  par_string=cpar%ds_gain_fwhm(i))
       call get_parameter_hashtable(htbl, 'BAND_DEFAULT_GAIN'//itext, len_itext=len_itext,    par_dp=cpar%ds_defaults(i,GAIN))
       call get_parameter_hashtable(htbl, 'BAND_DEFAULT_NOISEAMP'//itext, len_itext=len_itext,par_dp=cpar%ds_defaults(i,NOISEAMP))
       call get_parameter_hashtable(htbl, 'BAND_COMPONENT_SENSITIVITY'//itext, len_itext=len_itext, &
            & par_string=cpar%ds_component_sensitivity(i))

       !read in all TOD parameters
       call get_parameter_hashtable(htbl, 'BAND_TOD_TYPE'//itext, len_itext=len_itext, par_string=cpar%ds_tod_type(i))
       if (cpar%enable_TOD_analysis .or. cpar%resamp_CMB) then
          if (trim(cpar%ds_tod_type(i)) /= 'none') then
             call get_parameter_hashtable(htbl, 'BAND_TOD_INIT_FROM_HDF'//itext, len_itext=len_itext, &
                  & par_string=cpar%ds_tod_initHDF(i))
          end if
       end if
       if (trim(cpar%ds_tod_type(i)) /= 'none') then
          call get_parameter_hashtable(htbl, 'BAND_TOD_DETECTOR_LIST'//itext, len_itext=len_itext, &
               & par_string=cpar%ds_tod_dets(i))
       end if

       if (cpar%enable_TOD_analysis) then
          if (trim(cpar%ds_tod_type(i)) /= 'none') then
             !all other tod things
             call get_parameter_hashtable(htbl, 'BAND_TOD_MAIN_PROCMASK'//itext, len_itext=len_itext, &
                  & par_string=cpar%ds_tod_procmask1(i))
             call get_parameter_hashtable(htbl, 'BAND_TOD_SMALL_PROCMASK'//itext, len_itext=len_itext, &
                  & par_string=cpar%ds_tod_procmask2(i))
             call get_parameter_hashtable(htbl, 'BAND_TOD_FILELIST'//itext, len_itext=len_itext, &
                  & par_string=cpar%ds_tod_filelist(i))
             call get_parameter_hashtable(htbl, 'BAND_TOD_START_SCANID'//itext, len_itext=len_itext, &
                  & par_int=cpar%ds_tod_scanrange(i,1))
             call get_parameter_hashtable(htbl, 'BAND_TOD_END_SCANID'//itext, len_itext=len_itext, &
                  & par_int=cpar%ds_tod_scanrange(i,2))
             call get_parameter_hashtable(htbl, 'BAND_TOD_TOT_NUMSCAN'//itext, len_itext=len_itext, &
                  & par_int=cpar%ds_tod_tot_numscan(i))
             call get_parameter_hashtable(htbl, 'BAND_TOD_FLAG'//itext, len_itext=len_itext, &
                  & par_int=cpar%ds_tod_flag(i))
             call get_parameter_hashtable(htbl, 'BAND_TOD_ORBITAL_ONLY_ABSCAL'//itext, len_itext=len_itext, &
                  & par_lgt=cpar%ds_tod_orb_abscal(i))
             call get_parameter_hashtable(htbl, 'BAND_TOD_RIMO'//itext, len_itext=len_itext, &
                  & par_string=cpar%ds_tod_instfile(i))
             call get_parameter_hashtable(htbl, 'BAND_TOD_BP_INIT_PROP'//itext, len_itext=len_itext, &
                  & par_string=cpar%ds_tod_bp_init(i))
          end if
       end if

       do j = 1, cpar%num_smooth_scales
          call int2string(j, jtext)          
          call get_parameter_hashtable(htbl, 'BAND_NOISE_RMS'//itext//'_SMOOTH'//jtext, &
               & par_string=cpar%ds_noise_rms_smooth(i,j))
       end do

       if (cpar%enable_TOD_analysis) then
          call get_parameter_hashtable(htbl, 'BAND_TOD_TYPE'//itext, len_itext=len_itext, &
               & par_string=cpar%ds_tod_type(i))
       end if

    end do

    ! Convert to proper internal units where necessary
    cpar%ds_nu_c = cpar%ds_nu_c * 1d9                           ! From GHz to Hz
    
  end subroutine read_data_params_hash


  subroutine read_component_params_hash(htbl, cpar)
    implicit none

    type(hash_tbl_sll), intent(in) :: htbl
    type(comm_params),  intent(inout) :: cpar

    integer(i4b)       :: i, n, len_itext
    real(dp)           :: amp, lat, lon
    character(len=2)   :: itext
    character(len=512) :: maskfile, tokens(4)
    
    len_itext=len(trim(itext))
    call get_parameter_hashtable(htbl, 'INSTRUMENT_PARAM_FILE', par_string=cpar%cs_inst_parfile)
    call get_parameter_hashtable(htbl, 'INIT_INSTRUMENT_FROM_HDF', par_string=cpar%cs_init_inst_hdf)
    call get_parameter_hashtable(htbl, 'NUM_SIGNAL_COMPONENTS', par_int=cpar%cs_ncomp_tot)
    call get_parameter_hashtable(htbl, 'NUM_CG_SAMPLING_GROUPS', par_int=cpar%cg_num_user_samp_groups)

    do i = 1, cpar%cg_num_user_samp_groups
       call int2string(i, itext)
       call get_parameter_hashtable(htbl, 'SAMPLING_GROUP'//itext, par_string=cpar%cg_samp_group(i))
       call get_parameter_hashtable(htbl, 'SAMPLING_GROUP_MASK'//itext, par_string=cpar%cg_samp_group_mask(i))
    end do

    n = cpar%cs_ncomp_tot
    allocate(cpar%cs_include(n), cpar%cs_label(n), cpar%cs_type(n), cpar%cs_class(n))
    allocate(cpar%cs_polarization(n), cpar%cs_nside(n), cpar%cs_lmax_amp(n), cpar%cs_lmax_ind(n))
    allocate(cpar%cs_l_apod(n), cpar%cs_output_EB(n), cpar%cs_initHDF(n))
    allocate(cpar%cs_unit(n), cpar%cs_nu_ref(n,3), cpar%cs_cltype(n), cpar%cs_cl_poltype(n))
    allocate(cpar%cs_clfile(n), cpar%cs_binfile(n), cpar%cs_band_ref(n))
    allocate(cpar%cs_lpivot(n), cpar%cs_mask(n), cpar%cs_mono_prior(n), cpar%cs_fwhm(n), cpar%cs_poltype(MAXPAR,n))
    allocate(cpar%cs_latmask(n), cpar%cs_defmask(n))
    allocate(cpar%cs_indmask(n), cpar%cs_amp_rms_scale(n))
    allocate(cpar%cs_cl_amp_def(n,3), cpar%cs_cl_beta_def(n,3), cpar%cs_cl_prior(n,2))
    allocate(cpar%cs_input_amp(n), cpar%cs_prior_amp(n), cpar%cs_input_ind(MAXPAR,n))
    allocate(cpar%cs_theta_def(MAXPAR,n), cpar%cs_p_uni(n,2,MAXPAR), cpar%cs_p_gauss(n,2,MAXPAR))
    allocate(cpar%cs_catalog(n), cpar%cs_SED_template(4,n), cpar%cs_cg_scale(n))!, cpar%cs_cg_samp_group(n))
    allocate(cpar%cs_ptsrc_template(n), cpar%cs_output_ptsrc_beam(n), cpar%cs_min_src_dist(n))
    allocate(cpar%cs_auxpar(MAXAUXPAR,n), cpar%cs_apply_pos_prior(n))
    allocate(cpar%cs_nu_min(n,MAXPAR), cpar%cs_nu_max(n,MAXPAR), cpar%cs_burn_in(n))
    allocate(cpar%cs_smooth_scale(n,MAXPAR), cpar%cs_apply_jeffreys(n))
    do i = 1, n
       call int2string(i, itext)
       call get_parameter_hashtable(htbl, 'INCLUDE_COMP'//itext, len_itext=len_itext,         par_lgt=cpar%cs_include(i))
       call get_parameter_hashtable(htbl, 'COMP_INIT_FROM_HDF'//itext, len_itext=len_itext,   par_string=cpar%cs_initHDF(i))
       call get_parameter_hashtable(htbl, 'COMP_LABEL'//itext, len_itext=len_itext,           par_string=cpar%cs_label(i))
       call get_parameter_hashtable(htbl, 'COMP_TYPE'//itext, len_itext=len_itext,            par_string=cpar%cs_type(i))
       call get_parameter_hashtable(htbl, 'COMP_CLASS'//itext, len_itext=len_itext,           par_string=cpar%cs_class(i))
       !call get_parameter_hashtable(htbl, 'COMP_CG_SAMPLE_GROUP'//itext, len_itext=len_itext, par_int=cpar%cs_cg_samp_group(i))
       !if (.not. cpar%cs_include(i)) cpar%cs_cg_samp_group(i) = 0
       if (trim(cpar%cs_type(i)) == 'md') then
          call get_parameter_hashtable(htbl, 'COMP_POLARIZATION'//itext, len_itext=len_itext, &
               & par_lgt=cpar%cs_polarization(i))
          call get_parameter_hashtable(htbl, 'COMP_MD_DEFINITION_FILE'//itext, len_itext=len_itext, &
               & par_string=cpar%cs_SED_template(1,i))
       else if (trim(cpar%cs_class(i)) == 'template') then
          call get_parameter_hashtable(htbl, 'COMP_TEMPLATE_DEFINITION_FILE'//itext, len_itext=len_itext, &
               & par_string=cpar%cs_SED_template(1,i))
          call get_parameter_hashtable(htbl, 'COMP_DEFAULT_AMPLITUDE'//itext, len_itext=len_itext, &
               & par_dp=cpar%cs_theta_def(1,i))
          call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_MEAN'//itext, len_itext=len_itext, &
               & par_dp=cpar%cs_p_gauss(i,1,1))
          call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_RMS'//itext, len_itext=len_itext,  &
               & par_dp=cpar%cs_p_gauss(i,2,1))
       else if (trim(cpar%cs_class(i)) == 'diffuse') then
          call get_parameter_hashtable(htbl, 'COMP_POLARIZATION'//itext, len_itext=len_itext,    par_lgt=cpar%cs_polarization(i))
          if (cpar%cs_polarization(i)) &
               & call get_parameter_hashtable(htbl, 'COMP_OUTPUT_EB_MAP'//itext, len_itext=len_itext, par_lgt=cpar%cs_output_EB(i))
          call get_parameter_hashtable(htbl, 'COMP_CG_SCALE'//itext, len_itext=len_itext,        par_dp=cpar%cs_cg_scale(i))
          call get_parameter_hashtable(htbl, 'COMP_NSIDE'//itext, len_itext=len_itext,           par_int=cpar%cs_nside(i))
          call get_parameter_hashtable(htbl, 'COMP_LMAX_AMP'//itext, len_itext=len_itext,        par_int=cpar%cs_lmax_amp(i))
          call get_parameter_hashtable(htbl, 'COMP_L_APOD'//itext, len_itext=len_itext,          par_int=cpar%cs_l_apod(i))
          call get_parameter_hashtable(htbl, 'COMP_LMAX_IND'//itext, len_itext=len_itext,        par_int=cpar%cs_lmax_ind(i))
          call get_parameter_hashtable(htbl, 'COMP_UNIT'//itext, len_itext=len_itext,            par_string=cpar%cs_unit(i))
          call get_parameter_hashtable(htbl, 'COMP_NU_REF_T'//itext, len_itext=len_itext,          par_dp=cpar%cs_nu_ref(i,1))
          call get_parameter_hashtable(htbl, 'COMP_NU_REF_P'//itext, len_itext=len_itext,          par_dp=cpar%cs_nu_ref(i,2))
          cpar%cs_nu_ref(i,3) = cpar%cs_nu_ref(i,2)
          call get_parameter_hashtable(htbl, 'COMP_CL_TYPE'//itext, len_itext=len_itext,         par_string=cpar%cs_cltype(i))
          call get_parameter_hashtable(htbl, 'COMP_INPUT_AMP_MAP'//itext, len_itext=len_itext,   par_string=cpar%cs_input_amp(i))
          call get_parameter_hashtable(htbl, 'COMP_PRIOR_AMP_MAP'//itext, len_itext=len_itext,   par_string=cpar%cs_prior_amp(i))
          call get_parameter_hashtable(htbl, 'COMP_OUTPUT_FWHM'//itext, len_itext=len_itext,     par_dp=cpar%cs_fwhm(i))
          if (trim(cpar%cs_cltype(i)) == 'binned') then
             call get_parameter_hashtable(htbl, 'COMP_CL_BIN_FILE'//itext, len_itext=len_itext,     par_string=cpar%cs_binfile(i))
             call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_FILE'//itext, len_itext=len_itext, par_string=cpar%cs_clfile(i))
          else if (trim(cpar%cs_cltype(i)) == 'power_law' .or. &
               & trim(cpar%cs_cltype(i)) == 'exp' .or. trim(cpar%cs_cltype(i))=='gauss') then
             call get_parameter_hashtable(htbl, 'COMP_CL_POLTYPE'//itext, len_itext=len_itext,      par_int=cpar%cs_cl_poltype(i))
             call get_parameter_hashtable(htbl, 'COMP_CL_L_PIVOT'//itext, len_itext=len_itext,      par_int=cpar%cs_lpivot(i))
             call get_parameter_hashtable(htbl, 'COMP_CL_BETA_PRIOR_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_cl_prior(i,1))
             call get_parameter_hashtable(htbl, 'COMP_CL_BETA_PRIOR_RMS'//itext, len_itext=len_itext,&
                  & par_dp=cpar%cs_cl_prior(i,2))
             call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_AMP_T'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_cl_amp_def(i,1))
             call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_BETA_T'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_cl_beta_def(i,1))
             if (cpar%cs_polarization(i)) then
                call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_AMP_E'//itext, len_itext=len_itext, &
                     & par_dp=cpar%cs_cl_amp_def(i,2))
                call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_BETA_E'//itext, len_itext=len_itext,&
                     & par_dp=cpar%cs_cl_beta_def(i,2))
                call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_AMP_B'//itext, len_itext=len_itext, &
                     & par_dp=cpar%cs_cl_amp_def(i,3))
                call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_BETA_B'//itext, len_itext=len_itext, &  
                & par_dp=cpar%cs_cl_beta_def(i,3))
             end if
             cpar%cs_cl_amp_def(i,:) = cpar%cs_cl_amp_def(i,:) / cpar%cs_cg_scale(i)**2
          end if
          call get_parameter_hashtable(htbl, 'COMP_MONOPOLE_PRIOR'//itext, len_itext=len_itext, par_string=cpar%cs_mono_prior(i))
          call get_parameter_hashtable(htbl, 'COMP_MASK'//itext, len_itext=len_itext,            par_string=cpar%cs_mask(i))
          maskfile = adjustl(trim(cpar%cs_mask(i)))
          if (maskfile(1:4) == '|b|<') then
             read(maskfile(5:),*) cpar%cs_latmask(i)
             cpar%cs_latmask(i) = cpar%cs_latmask(i) * DEG2RAD
          else
             cpar%cs_latmask(i) = -1.d0
          end if
          cpar%cs_indmask(i) = 'fullsky'

          if (cpar%resamp_CMB .and. trim(cpar%cs_type(i)) == 'cmb') then
             call get_parameter_hashtable(htbl, 'COMP_DEFLATION_MASK'//itext, len_itext=len_itext, &
                  & par_string=cpar%cs_defmask(i))
          else
             cpar%cs_defmask(i) = 'fullsky'
          end if

          select case (trim(cpar%cs_type(i)))
          case ('cmb')
             call get_parameter_hashtable(htbl, 'CMB_DIPOLE_PRIOR', par_string=cpar%cmb_dipole_prior_mask)
             if (trim(cpar%cmb_dipole_prior_mask) /= 'none') then
                call get_tokens(trim(adjustl(cpar%cmb_dipole_prior_mask)), ';', tokens)
                cpar%cmb_dipole_prior_mask = tokens(1)
                read(tokens(2),*) amp
                read(tokens(3),*) lon
                read(tokens(4),*) lat
                cpar%cmb_dipole_prior(1) = -sqrt(4.d0*pi/3.d0) * amp * cos(lon*pi/180.d0) * sin((90.d0-lat)*pi/180.d0)
                cpar%cmb_dipole_prior(2) =  sqrt(4.d0*pi/3.d0) * amp * sin(lon*pi/180.d0) * sin((90.d0-lat)*pi/180.d0)
                cpar%cmb_dipole_prior(3) =  sqrt(4.d0*pi/3.d0) * amp *                      cos((90.d0-lat)*pi/180.d0)
             end if
          case ('power_law')
             call get_parameter_hashtable(htbl, 'COMP_BETA_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
             call get_parameter_hashtable(htbl, 'COMP_INPUT_BETA_MAP'//itext, len_itext=len_itext,        &
                  & par_string=cpar%cs_input_ind(1,i))
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_BETA'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_BETA_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_BETA_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_BETA_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_BETA_RMS'//itext, len_itext=len_itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_INDMASK'//itext, len_itext=len_itext,               &
                  & par_string=cpar%cs_indmask(i))
             call get_parameter_hashtable(htbl, 'COMP_BETA_SMOOTHING_SCALE'//itext, len_itext=len_itext,  &
                  & par_int=cpar%cs_smooth_scale(i,1))
             call get_parameter_hashtable(htbl, 'COMP_BETA_NU_MIN'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_min(i,1))
             call get_parameter_hashtable(htbl, 'COMP_BETA_NU_MAX'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_max(i,1))
             call get_parameter_hashtable(htbl, 'COMP_APPLY_JEFFREYS_PRIOR'//itext, len_itext=len_itext,   par_lgt=cpar%cs_apply_jeffreys(i))
             cpar%cs_apply_jeffreys(i) = .false. ! Disabled until properly debugged and validated
          case ('physdust')
             call get_parameter_hashtable(htbl, 'COMP_UMIN_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
             call get_parameter_hashtable(htbl, 'COMP_INPUT_UMIN_MAP'//itext, len_itext=len_itext,        &
                  & par_string=cpar%cs_input_ind(1,i))
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_UMIN'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_UMIN_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_UMIN_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_UMIN_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_UMIN_RMS'//itext, len_itext=len_itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_UMAX'//itext, len_itext=len_itext,  par_dp=cpar%cs_auxpar(1,i))
             call get_parameter_hashtable(htbl, 'COMP_GAMMA'//itext, len_itext=len_itext,  par_dp=cpar%cs_auxpar(2,i))
             call get_parameter_hashtable(htbl, 'COMP_ALPHA'//itext, len_itext=len_itext,  par_dp=cpar%cs_auxpar(3,i))
             call get_parameter_hashtable(htbl, 'COMP_SIL_AMP1_'//itext, len_itext=len_itext,  par_dp=cpar%cs_auxpar(4,i))
             call get_parameter_hashtable(htbl, 'COMP_SIL_AMP2_'//itext, len_itext=len_itext,  par_dp=cpar%cs_auxpar(5,i))
             call get_parameter_hashtable(htbl, 'COMP_CARB_AMP1_'//itext, len_itext=len_itext,  par_dp=cpar%cs_auxpar(6,i))
             call get_parameter_hashtable(htbl, 'COMP_CARB_AMP2_'//itext, len_itext=len_itext,  par_dp=cpar%cs_auxpar(7,i))
             call get_parameter_hashtable(htbl, 'COMP_SIL_FILE1_'//itext, len_itext=len_itext,  &
                  & par_string=cpar%cs_SED_template(1,i))
             call get_parameter_hashtable(htbl, 'COMP_SIL_FILE2_'//itext, len_itext=len_itext,  &
                  & par_string=cpar%cs_SED_template(2,i))
             call get_parameter_hashtable(htbl, 'COMP_CARB_FILE1_'//itext, len_itext=len_itext, &
                  & par_string=cpar%cs_SED_template(3,i))
             call get_parameter_hashtable(htbl, 'COMP_CARB_FILE2_'//itext, len_itext=len_itext, &
                  & par_string=cpar%cs_SED_template(4,i))
             call get_parameter_hashtable(htbl, 'COMP_INDMASK'//itext, len_itext=len_itext, &
                  & par_string=cpar%cs_indmask(i))
             call get_parameter_hashtable(htbl, 'COMP_APPLY_JEFFREYS_PRIOR'//itext, len_itext=len_itext,   par_lgt=cpar%cs_apply_jeffreys(i))
             write(*,*) 'ERROR -- smoothing not yet supported.'
             stop
          case ('spindust')
             call get_parameter_hashtable(htbl, 'COMP_NU_P_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
             call get_parameter_hashtable(htbl, 'COMP_INPUT_NU_P_MAP'//itext, len_itext=len_itext,        &
                  & par_string=cpar%cs_input_ind(1,i))
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_NU_P'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_NU_P_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_NU_P_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_NU_P_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_NU_P_RMS'//itext, len_itext=len_itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_SED_TEMPLATE'//itext, len_itext=len_itext,  &
                  & par_string=cpar%cs_SED_template(1,i))
             call get_parameter_hashtable(htbl, 'COMP_INDMASK'//itext, len_itext=len_itext,         par_string=cpar%cs_indmask(i))
             call get_parameter_hashtable(htbl, 'COMP_NU_P_SMOOTHING_SCALE'//itext, len_itext=len_itext,  &
                  & par_int=cpar%cs_smooth_scale(i,1))
             call get_parameter_hashtable(htbl, 'COMP_NU_P_NU_MIN'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_min(i,1))
             call get_parameter_hashtable(htbl, 'COMP_NU_P_NU_MAX'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_max(i,1))
             call get_parameter_hashtable(htbl, 'COMP_APPLY_JEFFREYS_PRIOR'//itext, len_itext=len_itext,   par_lgt=cpar%cs_apply_jeffreys(i))
          case ('spindust2')
             call get_parameter_hashtable(htbl, 'COMP_NU_P_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
             call get_parameter_hashtable(htbl, 'COMP_INPUT_NU_P_MAP'//itext, len_itext=len_itext,        &
                  & par_string=cpar%cs_input_ind(1,i))
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_NU_P'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_NU_P_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_NU_P_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_NU_P_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_NU_P_RMS'//itext, len_itext=len_itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_ALPHA_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(2,i))
             call get_parameter_hashtable(htbl, 'COMP_INPUT_ALPHA_MAP'//itext, len_itext=len_itext,        &
                  & par_string=cpar%cs_input_ind(2,i))
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_ALPHA'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(2,i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_ALPHA_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_ALPHA_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_ALPHA_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_ALPHA_RMS'//itext, len_itext=len_itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,2))
             call get_parameter_hashtable(htbl, 'COMP_SED_TEMPLATE'//itext, len_itext=len_itext,  &
                  & par_string=cpar%cs_SED_template(1,i))
             call get_parameter_hashtable(htbl, 'COMP_INDMASK'//itext, len_itext=len_itext,         par_string=cpar%cs_indmask(i))
             call get_parameter_hashtable(htbl, 'COMP_NU_P_SMOOTHING_SCALE'//itext, len_itext=len_itext,  &
                  & par_int=cpar%cs_smooth_scale(i,1))
             call get_parameter_hashtable(htbl, 'COMP_ALPHA_SMOOTHING_SCALE'//itext, len_itext=len_itext,  &
                  & par_int=cpar%cs_smooth_scale(i,2))
             call get_parameter_hashtable(htbl, 'COMP_NU_P_NU_MIN'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_min(i,1))
             call get_parameter_hashtable(htbl, 'COMP_NU_P_NU_MAX'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_max(i,1))
             call get_parameter_hashtable(htbl, 'COMP_ALPHA_NU_MIN'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_min(i,2))          
             call get_parameter_hashtable(htbl, 'COMP_ALPHA_NU_MAX'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_max(i,2))
             call get_parameter_hashtable(htbl, 'COMP_APPLY_JEFFREYS_PRIOR'//itext, len_itext=len_itext,   par_lgt=cpar%cs_apply_jeffreys(i))
          case ('MBB')
             call get_parameter_hashtable(htbl, 'COMP_BETA_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
             call get_parameter_hashtable(htbl, 'COMP_INPUT_BETA_MAP'//itext, len_itext=len_itext,        &
                  & par_string=cpar%cs_input_ind(1,i))
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_BETA'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_BETA_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_BETA_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_BETA_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_BETA_RMS'//itext, len_itext=len_itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_T_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(2,i))
             call get_parameter_hashtable(htbl, 'COMP_INPUT_T_MAP'//itext, len_itext=len_itext,        &
                  & par_string=cpar%cs_input_ind(2,i))
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_T'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(2,i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_T_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_T_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_T_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_T_RMS'//itext, len_itext=len_itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,2))
             call get_parameter_hashtable(htbl, 'COMP_INDMASK'//itext, len_itext=len_itext,         par_string=cpar%cs_indmask(i))
             call get_parameter_hashtable(htbl, 'COMP_BETA_SMOOTHING_SCALE'//itext, len_itext=len_itext,  &
                  & par_int=cpar%cs_smooth_scale(i,1))
             call get_parameter_hashtable(htbl, 'COMP_T_SMOOTHING_SCALE'//itext, len_itext=len_itext,  &
                  & par_int=cpar%cs_smooth_scale(i,2))
             call get_parameter_hashtable(htbl, 'COMP_BETA_NU_MIN'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_min(i,1))
             call get_parameter_hashtable(htbl, 'COMP_BETA_NU_MAX'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_max(i,1))
             call get_parameter_hashtable(htbl, 'COMP_T_NU_MIN'//itext, len_itext=len_itext,      par_dp=cpar%cs_nu_min(i,2))
             call get_parameter_hashtable(htbl, 'COMP_T_NU_MAX'//itext, len_itext=len_itext,      par_dp=cpar%cs_nu_max(i,2))
             call get_parameter_hashtable(htbl, 'COMP_APPLY_JEFFREYS_PRIOR'//itext, len_itext=len_itext,   par_lgt=cpar%cs_apply_jeffreys(i))
          case ('freefree')
!!$             call get_parameter_hashtable(htbl, 'COMP_EM_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
!!$             call get_parameter_hashtable(htbl, 'COMP_INPUT_EM_MAP'//itext, len_itext=len_itext,        &
!!$                  & par_string=cpar%cs_input_ind(1,i))
!!$             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_EM'//itext, len_itext=len_itext,          &
!!$                  & par_dp=cpar%cs_theta_def(1,i))
!!$             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_EM_LOW'//itext, len_itext=len_itext,    &
!!$                  & par_dp=cpar%cs_p_uni(i,1,1))
!!$             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_EM_HIGH'//itext, len_itext=len_itext,   &
!!$                  & par_dp=cpar%cs_p_uni(i,2,1))
!!$             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_EM_MEAN'//itext, len_itext=len_itext, &
!!$                  & par_dp=cpar%cs_p_gauss(i,1,1))
!!$             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_EM_RMS'//itext, len_itext=len_itext,  &
!!$                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_TE_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
             call get_parameter_hashtable(htbl, 'COMP_INPUT_TE_MAP'//itext, len_itext=len_itext,        &
                  & par_string=cpar%cs_input_ind(1,i))
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_TE'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_TE_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_TE_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_TE_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_TE_RMS'//itext, len_itext=len_itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_INDMASK'//itext, len_itext=len_itext,         par_string=cpar%cs_indmask(i))
!!$             call get_parameter_hashtable(htbl, 'COMP_EM_SMOOTHING_SCALE'//itext, len_itext=len_itext,  &
!!$                  & par_int=cpar%cs_smooth_scale(i,1))
             call get_parameter_hashtable(htbl, 'COMP_T_E_SMOOTHING_SCALE'//itext, len_itext=len_itext,  &
                  & par_int=cpar%cs_smooth_scale(i,1))
!!$             call get_parameter_hashtable(htbl, 'COMP_EM_NU_MIN'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_min(i,1))
!!$             call get_parameter_hashtable(htbl, 'COMP_EM_NU_MAX'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_max(i,1))
             call get_parameter_hashtable(htbl, 'COMP_T_E_NU_MIN'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_min(i,1))
             call get_parameter_hashtable(htbl, 'COMP_T_E_NU_MAX'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_max(i,1))
             call get_parameter_hashtable(htbl, 'COMP_APPLY_JEFFREYS_PRIOR'//itext, len_itext=len_itext,   par_lgt=cpar%cs_apply_jeffreys(i))
          case ('line')
             call get_parameter_hashtable(htbl, 'COMP_LINE_TEMPLATE'//itext, len_itext=len_itext,  &
                  & par_string=cpar%cs_SED_template(1,i))
             call get_parameter_hashtable(htbl, 'COMP_BAND_REF'//itext, len_itext=len_itext, &
                  & par_string=cpar%cs_band_ref(i))
             call get_parameter_hashtable(htbl, 'COMP_INDMASK'//itext, len_itext=len_itext,         par_string=cpar%cs_indmask(i))
             call get_parameter_hashtable(htbl, 'COMP_APPLY_JEFFREYS_PRIOR'//itext, len_itext=len_itext,   par_lgt=cpar%cs_apply_jeffreys(i))
          end select

       else if (trim(cpar%cs_class(i)) == 'ptsrc') then
          call get_parameter_hashtable(htbl, 'COMP_POLARIZATION'//itext, len_itext=len_itext,    par_lgt=cpar%cs_polarization(i))
          call get_parameter_hashtable(htbl, 'COMP_CATALOG'//itext, len_itext=len_itext,  par_string=cpar%cs_catalog(i))
          call get_parameter_hashtable(htbl, 'COMP_PTSRC_TEMPLATE'//itext, len_itext=len_itext,  &
               & par_string=cpar%cs_ptsrc_template(i))
          call get_parameter_hashtable(htbl, 'COMP_BURN_IN_ON_FIRST_SAMPLE'//itext, &
               & len_itext=len_itext,  par_lgt=cpar%cs_burn_in(i))
          call get_parameter_hashtable(htbl, 'COMP_AMP_RMS_SCALE_FACTOR'//itext, len_itext=len_itext,  &
               & par_dp=cpar%cs_amp_rms_scale(i))
          call get_parameter_hashtable(htbl, 'COMP_OUTPUT_PTSRC_TEMPLATE'//itext, len_itext=len_itext,  &
               & par_lgt=cpar%cs_output_ptsrc_beam(i))
          call get_parameter_hashtable(htbl, 'COMP_APPLY_POSITIVITY_PRIOR'//itext, len_itext=len_itext, &
               & par_lgt=cpar%cs_apply_pos_prior(i))
          !if (cpar%cs_apply_pos_prior(i)) cpar%cs_cg_samp_group(i) = 0
          call get_parameter_hashtable(htbl, 'COMP_MIN_DIST_BETWEEN_SRC'//itext, &
               & len_itext=len_itext, par_dp=cpar%cs_min_src_dist(i))
          call get_parameter_hashtable(htbl, 'COMP_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
          call get_parameter_hashtable(htbl, 'COMP_NSIDE'//itext, len_itext=len_itext,    par_int=cpar%cs_nside(i))
          call get_parameter_hashtable(htbl, 'COMP_NU_REF'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_ref(i,1))
          cpar%cs_nu_ref(i,2:3) = cpar%cs_nu_ref(i,1)
          call get_parameter_hashtable(htbl, 'COMP_CG_SCALE'//itext, len_itext=len_itext, par_dp=cpar%cs_cg_scale(i))
          select case (trim(cpar%cs_type(i)))
          case ('radio')
             call get_parameter_hashtable(htbl, 'COMP_APPLY_JEFFREYS_PRIOR'//itext, len_itext=len_itext,   par_lgt=cpar%cs_apply_jeffreys(i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_ALPHA_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_ALPHA_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_ALPHA_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_ALPHA_RMS'//itext, len_itext=len_itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_ALPHA'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_BETA_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_BETA_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_BETA_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_BETA_RMS'//itext, len_itext=len_itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,2))
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_BETA'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(2,i))
             call get_parameter_hashtable(htbl, 'COMP_ALPHA_NU_MIN'//itext, len_itext=len_itext,  par_dp=cpar%cs_nu_min(i,1))
             call get_parameter_hashtable(htbl, 'COMP_ALPHA_NU_MAX'//itext, len_itext=len_itext,  par_dp=cpar%cs_nu_max(i,1))
             call get_parameter_hashtable(htbl, 'COMP_BETA_NU_MIN'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_min(i,2))
             call get_parameter_hashtable(htbl, 'COMP_BETA_NU_MAX'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_max(i,2))
          case ('fir')
             call get_parameter_hashtable(htbl, 'COMP_APPLY_JEFFREYS_PRIOR'//itext, len_itext=len_itext,   par_lgt=cpar%cs_apply_jeffreys(i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_BETA_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_BETA_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_BETA_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_BETA_RMS'//itext, len_itext=len_itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))             
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_BETA'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_T_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_T_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_T_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_T_RMS'//itext, len_itext=len_itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,2))
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_T'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(2,i))             
             call get_parameter_hashtable(htbl, 'COMP_BETA_NU_MIN'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_min(i,1))
             call get_parameter_hashtable(htbl, 'COMP_BETA_NU_MAX'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_max(i,1))
             call get_parameter_hashtable(htbl, 'COMP_T_NU_MIN'//itext, len_itext=len_itext,      par_dp=cpar%cs_nu_min(i,2))
             call get_parameter_hashtable(htbl, 'COMP_T_NU_MAX'//itext, len_itext=len_itext,      par_dp=cpar%cs_nu_max(i,2))
          end select
       end if
       
    end do
    cpar%cs_ncomp           = count(cpar%cs_include)
    !cpar%cg_num_samp_groups = maxval(cpar%cs_cg_samp_group)

    ! Convert to proper units
    cpar%cs_nu_ref = 1d9 * cpar%cs_nu_ref
    cpar%cs_nu_min = 1d9 * cpar%cs_nu_min
    cpar%cs_nu_max = 1d9 * cpar%cs_nu_max


  end subroutine read_component_params_hash




  ! ********************************************************
  !                     Utility routines
  ! ********************************************************


  subroutine parse_parameter(line, parname, found, par_int, par_char, par_string, par_sp, par_dp, par_lgt)
    implicit none
    character(len=*)           :: line, parname
    character(len=256)         :: toks(2), key, value, par
    logical(lgt)               :: found
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt

    integer(i4b) :: n

    call get_tokens(trim(line), "=", group="''" // '""', maxnum=2, toks=toks, num=n)
    if(n < 2) then
       found = .false.
       return
    end if
    key = get_token(toks(1), " ", 1, group="''" // '""')
    value = get_token(toks(2), " ", 1, group="''" // '""')
    par = parname
    call tolower(key)
    call tolower(par)
    if (trim(key) == trim(par)) then
       if (present(par_int)) then
          read(value,*) par_int
       elseif (present(par_char)) then
          read(value,*) par_char
       elseif (present(par_string)) then
          read(value,*) par_string
       elseif (present(par_sp)) then
          read(value,*) par_sp
       elseif (present(par_dp)) then
          read(value,*) par_dp
       elseif (present(par_lgt)) then
          read(value,*) par_lgt
       else
          write(*,*) "get_parameter: Reached unreachable point!"
       end if
       found = .true.
    else
       found = .false.
    end if
  end subroutine parse_parameter

  !gets parameters from input arguments in Commander call
  subroutine get_parameter_arg(parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present, desc)
    implicit none
    character(len=*)           :: parname
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present
    character(len=*), optional :: desc

    character(len=512) :: line
    integer(i4b)       :: i
    logical(lgt)       :: found
    do i = 1, iargc()
       call getarg(i, line)
       if(line(1:2) /= "--") cycle
       call parse_parameter(line(3:), parname, found, par_int, par_char, par_string, par_sp, par_dp, par_lgt)
       if(found) then
          if(present(par_present)) par_present = .true.
          return
       end if
    end do
    if(present(par_present)) then
       par_present = .false.
    else
       write(*,*) "get_parameter_arg: Fatal error: Cannot find " // trim(parname) // " in argument list!"
       if(present(desc)) write(*,*) trim(desc)
       stop
    end if
  end subroutine

  ! Loops through the parameter files and children, counting lines.
  ! No error reporting.
  subroutine dump_expanded_paramfile(parfile, outfile)
    implicit none
    character(len=*)           :: parfile, outfile
    integer(i4b), parameter    :: maxdepth = 256
    integer(i4b)               :: depth, units(maxdepth), i, num, ounit
    character(len=1024)        :: key, value, arg

    num = 0
    depth = 1
    ounit = getlun()
    open(ounit,file=outfile,action="write")
    write(ounit,fmt="(a)",advance="no") '# Arguments:'
    do i = 1, iargc()
       call getarg(i, arg)
       write(ounit,fmt="(a)",advance="no") " '" // trim(arg) // "'"
    end do
    write(ounit,*)

    units(depth) = getlun()
    open(units(depth),file=parfile,status="old",action="read")
    do while(depth >= 1)
       read(units(depth),*,end=1) key
       backspace(units(depth))
       if (key=='@INCLUDE') then
          ! Recurse into the new file
          read(units(depth),*,end=1) key, value
          write(ounit,fmt='(a)') pad("",depth-1," ") // "# File: " // trim(value)
          depth=depth+1
          units(depth) = getlun()
          open(units(depth),file=value,status="old")
       else
          read(units(depth),fmt="(a)") value
          write(ounit,fmt='(a)') pad("",depth-1," ") // trim(value)
       end if
       cycle
1      close(units(depth))
       depth = depth-1
    end do
    close(ounit)
  end subroutine dump_expanded_paramfile


  function get_token(string, sep, num, group, allow_empty) result(res)
    implicit none
    character(len=*)           :: string, sep
    character(len=len(string)) :: res
    character(len=*), optional :: group
    logical(lgt),     optional :: allow_empty
    integer(i4b)               :: i, num, ext(2)
    ext = -1
    do i = 1, num
       call tokenize(string, sep, ext, group, allow_empty)
    end do
    res = string(ext(1):ext(2))
  end function get_token

  ! Fill all tokens into toks, and the num filled into num
  subroutine get_tokens(string, sep, toks, num, group, maxnum, allow_empty)
    implicit none
    character(len=*) :: string, sep
    character(len=*) :: toks(:)
    character(len=*), optional :: group
    integer(i4b),     optional :: num, maxnum
    logical(lgt),     optional :: allow_empty
    integer(i4b) :: n, ext(2), nmax
    ext = -1
    n = 0
    nmax = size(toks); if(present(maxnum)) nmax = maxnum
    call tokenize(string, sep, ext, group, allow_empty)
    do while(ext(1) > 0 .and. n < nmax)
       n = n+1
       toks(n) = string(ext(1):ext(2))
       call tokenize(string, sep, ext, group, allow_empty)
    end do
    if(present(num)) num = n
  end subroutine get_tokens

  function has_token(token, string, sep, group, allow_empty) result(res)
    implicit none
    character(len=*) :: token, string, sep
    character(len=*), optional :: group
    logical(lgt),     optional :: allow_empty
    logical(lgt)     :: res
    integer(i4b)     :: ext(2)
    res = .true.
    ext = -1
    call tokenize(string, sep, ext, group, allow_empty)
    do while(ext(1) > 0)
       if(string(ext(1):ext(2)) == trim(token)) return
       call tokenize(string, sep, ext, group, allow_empty)
    end do
    res = .false.
  end function has_token

  function num_tokens(string, sep, group, allow_empty) result(res)
    implicit none
    character(len=*) :: string, sep
    character(len=*), optional :: group
    logical(lgt),     optional :: allow_empty
    integer(i4b)     :: res, ext(2)
    ext = -1
    res = 0
    call tokenize(string, sep, ext, group, allow_empty)
    do while(ext(1) > 0)
       res = res+1
       call tokenize(string, sep, ext, group, allow_empty)
    end do
  end function num_tokens

  subroutine tokenize(string, sep, ext, group, allow_empty)
    implicit none
    character(len=*) :: string, sep
    character(len=*), optional   :: group
    character(len=256)  :: op, cl
    integer(i4b), save           :: level(256), nl
    integer(i4b), intent(inout)  :: ext(2)
    logical(lgt), optional       :: allow_empty

    integer(i4b) :: i, j, o, c, ng
    logical(lgt) :: intok, hit, empty

    empty = .false.; if(present(allow_empty)) empty = allow_empty

    if(ext(2) >= len(string)) then
       ext = (/ 0, -1 /)
       return
    end if
    ng = 0
    if(present(group)) then
       ng = len_trim(group)/2
       do i = 1, ng
          op(i:i) = group(2*i-1:2*i-1)
          cl(i:i) = group(2*i:2*i)
       end do
    end if
    if(ext(2) <= 0) then
       level = 0
       nl = 0
    end if
    intok = .false.
    j     = 1
    do i = ext(2)+2, len(string)
       hit = .false.
       c = index(cl(1:ng), string(i:i))
       if(c /= 0) then; if(level(c) > 0) then
          level(c) = level(c) - 1
          if(level(c) == 0) nl = nl - 1
          hit = .true.
       end if; end if
       if(nl == 0) then
          ! Are we in a separator or not
          if(index(sep, string(i:i)) == 0) then
             ! Nope, so we must be in a token. Register start of token.
             if(.not. intok) then
                j = i
                intok = .true.
             end if
          else
             ! Yes. This either means that a token is done, and we should
             ! return it, or that we are waiting for a new token, in
             ! which case do nothing.
             if(intok) then
                ext = (/ j, i-1 /)
                return
             elseif(empty) then
                ext = (/ i, i-1 /)
                return
             end if
          end if
       end if
       o = index(op(1:ng), string(i:i))
       if(o /= 0 .and. .not. hit) then
          if(level(o) == 0) nl = nl + 1
          level(o) = level(o) + 1
       end if
    end do
    ! Handle last token
    if(intok) then
       ext = (/ j, i-1 /)
    elseif(empty) then
       ext = (/ i, i-1 /)
    else
       ext = (/ 0, -1 /)
    end if
  end subroutine tokenize

  subroutine validate_params(cpar)
    implicit none
    type(comm_params), intent(inout) :: cpar

    integer(i4b) :: i, j
    character(len=512) :: datadir, chaindir

    datadir  = trim(cpar%datadir) // '/'
    chaindir = trim(cpar%outdir) // '/'
    
    do i = 1, cpar%cg_num_user_samp_groups
       if (trim(cpar%cg_samp_group_mask(i)) /= 'fullsky') call validate_file(trim(datadir)//trim(cpar%cg_samp_group_mask(i)))
    end do

    ! Check that all dataset files exist
    do i = 1, cpar%numband
       if (.not. cpar%ds_active(i)) cycle

       call validate_file(trim(datadir)//trim(cpar%ds_mapfile(i)))           ! Map file
       call validate_file(trim(datadir)//trim(cpar%ds_noisefile(i)))         ! Noise file
       if (trim(cpar%ds_maskfile(i)) /= 'fullsky') &
            & call validate_file(trim(datadir)//trim(cpar%ds_maskfile(i)))   ! Mask file
       if (trim(cpar%ds_bptype(i)) /= 'delta') &
            & call validate_file(trim(datadir)//trim(cpar%ds_bpfile(i)))     ! Bandpass
       call validate_file(trim(datadir)//trim(cpar%ds_pixwin(i)))            ! Pixel window
       call validate_file(trim(datadir)//trim(cpar%ds_blfile(i)))            ! Beam b_l file
       if (trim(cpar%ds_btheta_file(i)) /= 'none') &
            & call validate_file(trim(datadir)//trim(cpar%ds_btheta_file(i))) ! Point source file
       do j = 1, cpar%num_smooth_scales
          if (trim(cpar%ds_noise_rms_smooth(i,j)) /= 'none' .and. trim(cpar%ds_noise_rms_smooth(i,j))/= 'native') &
               & call validate_file(trim(datadir)//trim(cpar%ds_noise_rms_smooth(i,j)))  ! Smoothed RMS file
       end do

       if (cpar%enable_TOD_analysis .and. trim(cpar%ds_tod_type(i)) /= 'none') then
          call validate_file(trim(datadir)//trim(cpar%ds_tod_procmask1(i)))  ! Procmask1
          call validate_file(trim(datadir)//trim(cpar%ds_tod_procmask2(i)))  ! Procmask2
          call validate_file(trim(datadir)//trim(cpar%ds_tod_filelist(i)))   ! Filelist
          call validate_file(trim(datadir)//trim(cpar%ds_tod_instfile(i)))   ! Instrument file, RIMO
          call validate_file(trim(datadir)//trim(cpar%ds_tod_bp_init(i)))    ! BP prop and init
       end if

    end do

    ! Instrument data base
    if (trim(cpar%cs_inst_parfile) /= 'none') &
         & call validate_file(trim(datadir)//trim(cpar%cs_inst_parfile))  ! Instrument data base
    if (trim(cpar%ds_sourcemask) /= 'none') &
         & call validate_file(trim(datadir)//trim(cpar%ds_sourcemask))    ! Source mask
    if (trim(cpar%ds_procmask) /= 'none') &
         & call validate_file(trim(datadir)//trim(cpar%ds_procmask))      ! Source mask

    ! Check component files
    do i = 1, cpar%cs_ncomp_tot
       if (.not. cpar%cs_include(i)) cycle

       if (trim(cpar%cs_type(i)) == 'md') then
          call validate_file(trim(datadir)//trim(cpar%cs_SED_template(1,i)))
       else if (trim(cpar%cs_class(i)) == 'diffuse') then
          if (trim(cpar%cs_input_amp(i)) /= 'none') &
               call validate_file(trim(datadir)//trim(cpar%cs_input_amp(i)))
          if (trim(cpar%cs_prior_amp(i)) /= 'none') &
               call validate_file(trim(datadir)//trim(cpar%cs_prior_amp(i)))
          if (trim(cpar%cs_cltype(i)) == 'binned') then
             call validate_file(trim(datadir)//trim(cpar%cs_binfile(i)))
             call validate_file(trim(datadir)//trim(cpar%cs_clfile(i)))             
          end if
          if (trim(cpar%cs_mask(i)) /= 'fullsky') &
               call validate_file(trim(datadir)//trim(cpar%cs_mask(i)))          
          if (trim(cpar%cs_indmask(i)) /= 'fullsky') &
               call validate_file(trim(datadir)//trim(cpar%cs_indmask(i)))          
          
          select case (trim(cpar%cs_type(i)))
          case ('power_law')
             if (trim(cpar%cs_input_ind(1,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(1,i)))
          case ('physdust')
             if (trim(cpar%cs_input_ind(1,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(1,i)))
             call validate_file(trim(datadir)//trim(cpar%cs_SED_template(1,i)))
             call validate_file(trim(datadir)//trim(cpar%cs_SED_template(2,i)))
             call validate_file(trim(datadir)//trim(cpar%cs_SED_template(3,i)))
             call validate_file(trim(datadir)//trim(cpar%cs_SED_template(4,i)))
          case ('spindust')
             if (trim(cpar%cs_input_ind(1,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(1,i)))
             call validate_file(trim(datadir)//trim(cpar%cs_SED_template(1,i)))             
          case ('spindust2')
             if (trim(cpar%cs_input_ind(1,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(1,i)))
             if (trim(cpar%cs_input_ind(2,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(2,i)))
             call validate_file(trim(datadir)//trim(cpar%cs_SED_template(1,i)))
          case ('MBB')
             if (trim(cpar%cs_input_ind(1,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(1,i)))
             if (trim(cpar%cs_input_ind(2,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(2,i)))
          case ('freefree')
!!$             if (trim(cpar%cs_input_ind(1,i)) /= 'default') &
!!$                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(1,i)))
             if (trim(cpar%cs_input_ind(1,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(1,i)))             
          case ('line')
             call validate_file(trim(datadir)//trim(cpar%cs_SED_template(1,i)))
          end select

       else if (trim(cpar%cs_class(i)) == 'ptsrc') then
          call validate_file(trim(datadir)//trim(cpar%cs_catalog(i)))
          call validate_file(trim(datadir)//trim(cpar%cs_ptsrc_template(i)), &
               & should_exist=.not. cpar%cs_output_ptsrc_beam(i))
       else if (trim(cpar%cs_type(i)) == 'template' .or. trim(cpar%cs_type(i)) == 'cmb_relquad') then
          call validate_file(trim(datadir)//trim(cpar%cs_SED_template(1,i)))
       end if
       
    end do
    
  end subroutine validate_params

  subroutine validate_file(filename, should_exist)
    implicit none
    character(len=*), intent(in)           :: filename
    logical(lgt),     intent(in), optional :: should_exist
    logical(lgt) :: exist, default
    default = .true.; if (present(should_exist)) default = should_exist
    inquire(file=trim(filename), exist=exist)
    if (exist .neqv. default) then
       if (default) then
          call report_error('Error: File does not exist = '//trim(filename))
       else
          call report_error('Error: File already exists = '//trim(filename))
       end if
    else
    end if
  end subroutine validate_file

  subroutine read_paramfile_to_ascii(paramfile,paramfile_cache)
    implicit none
    character(len=512),                            intent(in)  :: paramfile
    character(len=512), allocatable, dimension(:), intent(inout) :: paramfile_cache
    integer(i4b), parameter    :: maxdepth = 256
    integer(i4b)               :: depth, units(maxdepth), line_nr, paramfile_len, i
    character(len=512)         :: key, value, filenames(maxdepth), line
    ! read file to ascii array

    line_nr = 0
    depth = 1
    units(depth) = getlun()
    filenames(depth) = paramfile
    open(units(depth),file=trim(paramfile),status="old",err=4)
    do while(depth >= 1)
       read(units(depth),*,end=1) key
       if (key(1:1)=='#') cycle
       backspace(units(depth))

       if (key(1:1)=='@') then
          if(key == '@INCLUDE') then
             ! Recurse into the new file
             read(units(depth),*,end=1) key, value
             depth=depth+1
             units(depth) = getlun()
             filenames(depth) = value
             open(units(depth),file=value,status="old",err=2)
          else
             goto 3
          end if
       else
          read(units(depth),fmt="(a)") line
          !if we get here we have read a new line from the parameter file(s)
          line_nr = line_nr + 1
          write(paramfile_cache(line_nr),fmt="(a)") line
       end if
       cycle
       ! We get here if we reached the end of a file. Close it and
       ! return to the file above.
1      close(units(depth))
       !write(*,*) "Exiting file " // filenames(depth)
       depth = depth-1
    end do
    return

    ! ===== Error handling section ======

    ! Case 1: Include file error
2   write(*,*) "Error: Cannot open include file '" // trim(value) // "'"
    write(*,*) " in file " // trim(filenames(depth-1))
    do i = depth-2, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
    do i = depth-1, 1, -1; close(units(i)); end do
    stop

    ! Case 2: Directive error
3   write(*,*) "Error: Unrecognized directive '" // trim(key) //"'"
    write(*,*) " in file " // trim(filenames(depth))
    do i = depth-1, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
    do i = depth, 1, -1; close(units(i)); end do
    stop

    ! Case 3: Top level parameter file unreadable
4   write(*,*) "Error: Cannot open parameter file '" // trim(paramfile) // "'"
    stop

  end subroutine read_paramfile_to_ascii
    
  subroutine get_file_length(filename,length)
    implicit none 
    character(len=512), intent(in)  :: filename
    integer(i4b),       intent(out) :: length
    integer(i4b), parameter    :: maxdepth = 256
    integer(i4b)               :: depth, units(maxdepth), i
    character(len=512)         :: key, value, filenames(maxdepth), line


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    length = 0
    depth = 1
    units(depth) = getlun()
    filenames(depth) = filename
    open(units(depth),file=trim(filename),status="old",err=4)
    do while(depth >= 1)
       read(units(depth),*,end=1) key
       if (key(1:1)=='#') cycle
       backspace(units(depth))

       if (key(1:1)=='@') then
          if(key == '@INCLUDE') then
             ! Recurse into the new file
             read(units(depth),*,end=1) key, value
             depth=depth+1
             units(depth) = getlun()
             filenames(depth) = value
             open(units(depth),file=value,status="old",err=2)
          else
             goto 3
          end if
       else
          read(units(depth),fmt="(a)") line
          !if we get here we have read a new line from the parameter file(s)
          length = length + 1
       end if
       cycle
       ! We get here if we reached the end of a file. Close it and
       ! return to the file above.
1      close(units(depth))
       !write(*,*) "Exiting file " // filenames(depth)
       depth = depth-1
    end do
    return
    ! ===== Error handling section ======

    ! Case 1: Include file error
2   write(*,*) "Error: Cannot open include file '" // trim(value) // "'"
    write(*,*) " in file " // trim(filenames(depth-1))
    do i = depth-2, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
    do i = depth-1, 1, -1; close(units(i)); end do
    stop

    ! Case 2: Directive error
3   write(*,*) "Error: Unrecognized directive '" // trim(key) //"'"
    write(*,*) " in file " // trim(filenames(depth))
    do i = depth-1, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
    do i = depth, 1, -1; close(units(i)); end do
    stop

    ! Case 3: Top level parameter file unreadable
4   write(*,*) "Error: Cannot open parameter file '" // trim(filename) // "'"
    stop

  end subroutine get_file_length

  ! A subroutine for debugging
  subroutine print_ascii_parameter_file  !(paramfile_cache)
    implicit none
    !character(len=512), allocatable, dimension(:), intent(in) :: paramfile_cache
    integer(i4b)      :: unit
    character(len=32) :: paramfile_out

    unit = getlun()
    paramfile_out="read_ascii_paramfile.txt"
    open(unit, file=trim(paramfile_out),err=1)

    !If the opening of the output parameter file fails
1   write(*,*) "Error: Cannot open ascii file '" // trim(paramfile_out) // "'"
    stop
  end subroutine print_ascii_parameter_file

  ! filling the hash table with elements from the parameter file (ascii array) 
  subroutine put_ascii_into_hashtable(asciitbl,htbl)
    implicit none
    character(len=512), allocatable, dimension(:), intent(in) :: asciitbl
    type(hash_tbl_sll), intent(inout) :: htbl
    character(len=512) :: key, val
    character(len=256) :: toks(2)
    integer            :: i, n
    do i = 1,size(asciitbl)
       call get_tokens(trim(asciitbl(i)), "=", group="''" // '""', maxnum=2, toks=toks, num=n)
       if(n < 2) then ! only need the lines where one has 'key'='value'
          cycle
       end if
       key = get_token(toks(1), " ", 1, group="''" // '""')
       val = get_token(toks(2), " ", 1, group="''" // '""')
       call tolower(key)  ! we don't differentiate btw. upper and lower case
       if (key=="") cycle !we don't need blank lines
       call put_hash_tbl_sll(htbl,trim(key),trim(val)) 
    end do
    return

    write(*,*) "Error: Cannot read ascii line:", i, "line = '" // trim(asciitbl(i)) // "'"
    stop

  end subroutine put_ascii_into_hashtable

  ! read parameter from input argument or hash table
  subroutine get_parameter_hashtable(htbl, parname, len_itext, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present, desc)
    implicit none
    type(hash_tbl_sll), intent(in) :: htbl 
    character(len=*),   intent(in) :: parname
    integer(i4b),     optional :: len_itext
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present
    character(len=*), optional :: desc

    logical(lgt)               :: found
   
    found = .false.
    call get_parameter_arg(parname, par_int, par_char, par_string, par_sp, par_dp, par_lgt, found, desc)
    if(found) then
       if(present(par_present)) par_present = .true.
    else
       call get_parameter_from_hash(htbl, parname, len_itext, par_int, &
            & par_char, par_string, par_sp, par_dp, par_lgt, par_present, desc)
    end if
  end subroutine get_parameter_hashtable

  ! getting parameter value from hash table
  subroutine get_parameter_from_hash(htbl, parname, len_itext, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present, desc)
    implicit none
    type(hash_tbl_sll), intent(in) :: htbl
    character(len=*),   intent(in) :: parname
    integer(i4b),     optional :: len_itext
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present
    character(len=*), optional :: desc
    character(len=256)         :: key
    character(len=:), ALLOCATABLE   :: itext,jtext
    CHARACTER(len=:), ALLOCATABLE   :: val,val2,val3
    integer(i4b)                    :: i,j

    key=trim(parname)
    call tolower(key)
    call get_hash_tbl_sll(htbl,trim(key),val)
    if (.not. allocated(val)) then
       if (.not. present(len_itext)) goto 1
       allocate(character(len=len_itext) :: itext,jtext)
       itext=key(len(trim(key))-(len_itext-1):len(trim(key)))
       call get_hash_tbl_sll(htbl,'band_default_params'//trim(itext),val2)
       if (allocated(val2)) then
          read(val2,*) j
          if (j /= 0) then
             call int2string(j, jtext)
             call get_hash_tbl_sll(htbl,'band_default_params'//trim(jtext),val3)
             if (allocated(val3)) then
                read(val3,*) i
                if (i /= 0) goto 2
             end if
             call get_hash_tbl_sll(htbl,key(1:len(trim(key))-len_itext)//trim(jtext),val)
             if (.not. allocated(val)) goto 3
          else
             goto 1
          end if
       else
          goto 1
       end if
       deallocate(itext,jtext)
    end if
    
    if (present(par_int)) then
       read(val,*) par_int
    elseif (present(par_char)) then
       read(val,*) par_char
    elseif (present(par_string)) then
       read(val,*) par_string
    elseif (present(par_sp)) then
       read(val,*) par_sp
    elseif (present(par_dp)) then
       read(val,*) par_dp
    elseif (present(par_lgt)) then
       read(val,*) par_lgt
    else
       write(*,*) "get_parameter: Reached unreachable point!"
    end if

    deallocate(val)
    return

    !if (cpar%myid == cpar%root) then

1   write(*,*) "Error: Could not find parameter '" // trim(parname) // "'"
    write(*,*) ""
    stop


2   write(*,*) "Error: Recursive default parameters, bands " // &
         & trim(jtext) // " and " //trim(itext)
    write(*,*) ""
    stop

3   write(*,*) "Error: Could not find parameter '" // trim(parname) // &
         & "' from default '"//key(1:len(trim(key))-len_itext)//trim(jtext)//"'"
    write(*,*) ""
    stop
     
  end subroutine get_parameter_from_hash

  subroutine get_chainfile_and_samp(string, chainfile, initsamp)
    implicit none
    character(len=*),   intent(in)  :: string
    character(len=512), intent(out) :: chainfile
    integer(i4b),       intent(out) :: initsamp

    integer(i4b) :: i, num
    character(len=512), dimension(2) :: toks

    call get_tokens(string, ":", toks, num)    
    chainfile = toks(1)
    read(toks(2),*) initsamp

  end subroutine get_chainfile_and_samp

  subroutine define_cg_samp_groups(cpar)
    implicit none
    type(comm_params), intent(inout) :: cpar

    integer(i4b) :: i, j, k, n
    character(len=16), dimension(1000) :: comp_label

    ! Add user specified sample groups
    cpar%cg_num_samp_groups = cpar%cg_num_user_samp_groups 

    ! Add one sample group per component
    do i = 1, cpar%cs_ncomp_tot
       if (cpar%cs_include(i)) then
          cpar%cg_num_samp_groups                     = cpar%cg_num_samp_groups + 1
          cpar%cg_samp_group(cpar%cg_num_samp_groups) = trim(cpar%cs_label(i))
       end if
    end do

    ! Expand md type if present
    do i = 1, cpar%cg_num_samp_groups
       call get_tokens(cpar%cg_samp_group(i), ",", comp_label, n)
       do j = 1, n
          if (trim(comp_label(j)) == 'md') then
             do k = 1, cpar%numband
                if (cpar%ds_active(k)) cpar%cg_samp_group(i) = trim(cpar%cg_samp_group(i))//','//trim(cpar%ds_label(k))
             end do
          end if
       end do
    end do

    ! More groups may be defined here

    
    if (cpar%cg_num_samp_groups > MAXSAMPGROUP) then
       write(*,*) 'Error -- too many CG sampling groups defined. Increase MAXSAMPGROUP'
       stop
    end if

  end subroutine define_cg_samp_groups

end module comm_param_mod
