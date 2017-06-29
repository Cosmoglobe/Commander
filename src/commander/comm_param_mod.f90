module comm_param_mod
  use comm_utils
  use comm_status_mod
  implicit none

  ! Note: This module reads in the Commander parameter file as the first operation
  !       in the program. This is primarily intended to avoid crashes after hours
  !       of running because of user errors; catch these early, and report back
  !       problems. Then, copy parameters over to module structures if convenient
  !       at a later stage. 

  integer(i4b), parameter, private :: MAXPAR    = 10
  integer(i4b), parameter, private :: MAXAUXPAR = 10
  type(status_file)                :: status
  
  type comm_params

     ! MPI info
     integer(i4b) :: myid, numprocs, root = 0
     integer(i4b) :: myid_chain, numprocs_chain, comm_chain, mychain
     integer(i4b), dimension(MPI_STATUS_SIZE)          :: status

     ! Global parameters
     character(len=24)  :: operation
     integer(i4b)       :: verbosity, base_seed, numchain, num_smooth_scales
     integer(i4b)       :: num_gibbs_iter, num_ml_iter, init_samp
     character(len=512) :: chain_prefix, init_chain_prefix
     real(dp)           :: T_CMB
     character(len=512) :: MJysr_convention
     logical(lgt)       :: only_pol
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
     logical(lgt)       :: sample_signal_amplitudes
     
     ! Numerical parameters
     integer(i4b)      :: cg_lmax_precond, cg_maxiter, cg_num_samp_groups, cg_miniter
     real(dp)          :: cg_tol

     ! Data parameters
     integer(i4b)       :: numband
     character(len=512) :: datadir, ds_sourcemask, ds_procmask
     real(dp)           :: ds_fwhm_proc
     logical(lgt),       allocatable, dimension(:)   :: ds_active
     integer(i4b),       allocatable, dimension(:)   :: ds_period
     logical(lgt),       allocatable, dimension(:)   :: ds_polarization
     integer(i4b),       allocatable, dimension(:)   :: ds_nside
     integer(i4b),       allocatable, dimension(:)   :: ds_lmax
     character(len=512), allocatable, dimension(:)   :: ds_label
     character(len=512), allocatable, dimension(:)   :: ds_unit
     character(len=512), allocatable, dimension(:)   :: ds_noise_format
     character(len=512), allocatable, dimension(:)   :: ds_mapfile
     character(len=512), allocatable, dimension(:)   :: ds_noise_rms
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
     character(len=512), allocatable, dimension(:)   :: ds_gain_calib_comp
     integer(i4b),       allocatable, dimension(:)   :: ds_gain_lmin
     integer(i4b),       allocatable, dimension(:)   :: ds_gain_lmax
     character(len=512), allocatable, dimension(:)   :: ds_gain_apodmask
     character(len=512), allocatable, dimension(:)   :: ds_gain_fwhm
     real(dp),           allocatable, dimension(:,:) :: ds_defaults

     ! Component parameters
     character(len=512) :: cs_inst_parfile
     integer(i4b)       :: cs_ncomp, cs_ncomp_tot
     logical(lgt),       allocatable, dimension(:)     :: cs_include
     logical(lgt),       allocatable, dimension(:)     :: cs_initHDF
     character(len=512), allocatable, dimension(:)     :: cs_label
     character(len=512), allocatable, dimension(:)     :: cs_type
     character(len=512), allocatable, dimension(:)     :: cs_class
     logical(lgt),       allocatable, dimension(:)     :: cs_polarization
     real(dp),           allocatable, dimension(:)     :: cs_cg_scale
     integer(i4b),       allocatable, dimension(:)     :: cs_cg_samp_group
     integer(i4b),       allocatable, dimension(:)     :: cs_nside
     integer(i4b),       allocatable, dimension(:,:)   :: cs_poltype
     integer(i4b),       allocatable, dimension(:)     :: cs_lmax_amp
     integer(i4b),       allocatable, dimension(:)     :: cs_l_apod
     integer(i4b),       allocatable, dimension(:)     :: cs_lmax_ind
     character(len=512), allocatable, dimension(:)     :: cs_unit
     real(dp),           allocatable, dimension(:)     :: cs_nu_ref
     character(len=512), allocatable, dimension(:)     :: cs_band_ref
     real(dp),           allocatable, dimension(:)     :: cs_fwhm
     character(len=512), allocatable, dimension(:)     :: cs_cltype
     character(len=512), allocatable, dimension(:)     :: cs_clfile
     character(len=512), allocatable, dimension(:)     :: cs_binfile
     integer(i4b),       allocatable, dimension(:)     :: cs_lpivot
     character(len=512), allocatable, dimension(:)     :: cs_mask
     character(len=512), allocatable, dimension(:)     :: cs_indmask
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
     
  end type comm_params

contains

  ! ********************************************************
  !                     Driver routines
  ! ********************************************************
  subroutine read_comm_params(cpar)
    implicit none

    type(comm_params), intent(inout) :: cpar

    character(len=512) :: paramfile
    
    call getarg(1, paramfile)
    call read_global_params(paramfile, cpar)
    call read_data_params(paramfile, cpar)
    call read_component_params(paramfile, cpar)
    
  end subroutine read_comm_params

  subroutine initialize_mpi_struct(cpar, handle)
    implicit none
    type(comm_params), intent(inout) :: cpar
    type(planck_rng),  intent(out)   :: handle

    integer(i4b) :: i, j, m, n, ierr
    integer(i4b), allocatable, dimension(:,:) :: ind

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, cpar%myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, cpar%numprocs, ierr)
    cpar%root = 0
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

    ! Initialize random number generator
    if (cpar%myid == cpar%root) then
       call rand_init(handle, cpar%base_seed)
       do i = 1, cpar%numprocs-1
          j = nint(rand_uni(handle)*1000000.d0)
          call mpi_send(j, 1, MPI_INTEGER, i, 98, MPI_COMM_WORLD, ierr)
       end do
    else 
       call mpi_recv(j, 1, MPI_INTEGER, cpar%root, 98, MPI_COMM_WORLD, status, ierr)
       call rand_init(handle, j)
    end if
    
    deallocate(ind)

  end subroutine initialize_mpi_struct
  
  ! ********************************************************
  !              Specialized routines; one per module
  ! ********************************************************

  subroutine read_global_params(paramfile, cpar)
    implicit none

    character(len=*),  intent(in)    :: paramfile
    type(comm_params), intent(inout) :: cpar

    integer(i4b)     :: i
    character(len=2) :: itext

    call get_parameter(paramfile, 'VERBOSITY',                par_int=cpar%verbosity)
    call get_parameter(paramfile, 'OPERATION',                par_string=cpar%operation)

    call get_parameter(paramfile, 'BASE_SEED',                par_int=cpar%base_seed)
    call get_parameter(paramfile, 'NUMCHAIN',                 par_int=cpar%numchain)
    call get_parameter(paramfile, 'NUM_GIBBS_ITER',           par_int=cpar%num_gibbs_iter)
    call get_parameter(paramfile, 'NUM_ITER_WITH_ML_SEARCH',  par_int=cpar%num_ml_iter)
    call get_parameter(paramfile, 'CHAIN_PREFIX',             par_string=cpar%chain_prefix)
    call get_parameter(paramfile, 'INIT_CHAIN',               par_string=cpar%init_chain_prefix)
    call get_parameter(paramfile, 'INIT_SAMPLE_NUMBER',       par_int=cpar%init_samp)
    call get_parameter(paramfile, 'SAMPLE_ONLY_POLARIZATION', par_lgt=cpar%only_pol)

    call get_parameter(paramfile, 'CG_LMAX_PRECOND',          par_int=cpar%cg_lmax_precond)
    call get_parameter(paramfile, 'CG_MINITER',               par_int=cpar%cg_miniter)
    call get_parameter(paramfile, 'CG_MAXITER',               par_int=cpar%cg_maxiter)
    call get_parameter(paramfile, 'CG_TOLERANCE',             par_dp=cpar%cg_tol)

    call get_parameter(paramfile, 'T_CMB',                    par_dp=cpar%T_cmb)
    call get_parameter(paramfile, 'MJYSR_CONVENTION',         par_string=cpar%MJysr_convention)

    call get_parameter(paramfile, 'OUTPUT_DIRECTORY',         par_string=cpar%outdir)

    call get_parameter(paramfile, 'NSIDE_CHISQ',              par_int=cpar%nside_chisq)
    call get_parameter(paramfile, 'POLARIZATION_CHISQ',       par_lgt=cpar%pol_chisq)
    cpar%nmaps_chisq = 1; if (cpar%pol_chisq) cpar%nmaps_chisq = 3

    call get_parameter(paramfile, 'OUTPUT_MIXING_MATRIX',     par_lgt=cpar%output_mixmat)
    call get_parameter(paramfile, 'OUTPUT_RESIDUAL_MAPS',     par_lgt=cpar%output_residuals)
    call get_parameter(paramfile, 'OUTPUT_CHISQ_MAP',         par_lgt=cpar%output_chisq)
    call get_parameter(paramfile, 'OUTPUT_EVERY_NTH_CG_ITERATION', par_int=cpar%output_cg_freq)
    call get_parameter(paramfile, 'OUTPUT_CG_PRECOND_EIGENVALS', par_lgt=cpar%output_cg_eigenvals)
    call get_parameter(paramfile, 'OUTPUT_INPUT_MODEL',       par_lgt=cpar%output_input_model)
    call get_parameter(paramfile, 'IGNORE_GAIN_AND_BANDPASS_CORR', par_lgt=cpar%ignore_gain_bp)
    call get_parameter(paramfile, 'OUTPUT_DEBUG_SEDS',        par_lgt=cpar%output_debug_seds)
    call get_parameter(paramfile, 'OUTPUT_SIGNALS_PER_BAND',  par_lgt=cpar%output_sig_per_band)

    call get_parameter(paramfile, 'SAMPLE_SIGNAL_AMPLITUDES', par_lgt=cpar%sample_signal_amplitudes)

    call get_parameter(paramfile, 'NUM_SMOOTHING_SCALES',     par_int=cpar%num_smooth_scales)

    allocate(cpar%fwhm_smooth(cpar%num_smooth_scales))
    allocate(cpar%fwhm_postproc_smooth(cpar%num_smooth_scales))
    allocate(cpar%lmax_smooth(cpar%num_smooth_scales))
    allocate(cpar%nside_smooth(cpar%num_smooth_scales))
    allocate(cpar%pixwin_smooth(cpar%num_smooth_scales))
    do i = 1, cpar%num_smooth_scales
       call int2string(i, itext)
       call get_parameter(paramfile, 'FWHM_SMOOTHING_SCALE'//itext, par_dp=cpar%fwhm_smooth(i))
       call get_parameter(paramfile, 'FWHM_POSTPROC_SMOOTHING_SCALE'//itext, par_dp=cpar%fwhm_postproc_smooth(i))
       call get_parameter(paramfile, 'LMAX_SMOOTHING_SCALE'//itext, par_int=cpar%lmax_smooth(i))
       call get_parameter(paramfile, 'NSIDE_SMOOTHING_SCALE'//itext, par_int=cpar%nside_smooth(i))
       call get_parameter(paramfile, 'PIXWIN_SMOOTHING_SCALE'//itext, par_string=cpar%pixwin_smooth(i))
    end do

  end subroutine read_global_params


  subroutine read_data_params(paramfile, cpar)
    implicit none

    character(len=*),  intent(in)    :: paramfile
    type(comm_params), intent(inout) :: cpar

    integer(i4b)     :: i, j, n
    character(len=3) :: itext
    character(len=2) :: jtext
    
    call get_parameter(paramfile, 'NUMBAND',             par_int=cpar%numband)
    call get_parameter(paramfile, 'DATA_DIRECTORY',      par_string=cpar%datadir)
    call get_parameter(paramfile, 'SOURCE_MASKFILE',     par_string=cpar%ds_sourcemask)
    call get_parameter(paramfile, 'PROCESSING_MASKFILE', par_string=cpar%ds_procmask)
    call get_parameter(paramfile, 'PROC_SMOOTH_SCALE',   par_dp=cpar%ds_fwhm_proc)

    n = cpar%numband
    allocate(cpar%ds_active(n), cpar%ds_label(n))
    allocate(cpar%ds_polarization(n), cpar%ds_nside(n), cpar%ds_lmax(n))
    allocate(cpar%ds_unit(n), cpar%ds_noise_format(n), cpar%ds_mapfile(n))
    allocate(cpar%ds_noise_rms(n), cpar%ds_maskfile(n), cpar%ds_maskfile_calib(n))
    allocate(cpar%ds_noise_rms_smooth(n,cpar%num_smooth_scales))
    allocate(cpar%ds_samp_noiseamp(n), cpar%ds_noise_uni_fsky(n))
    allocate(cpar%ds_bptype(n), cpar%ds_nu_c(n), cpar%ds_bpfile(n), cpar%ds_bpmodel(n))
    allocate(cpar%ds_period(n), cpar%ds_beamtype(n), cpar%ds_blfile(n))
    allocate(cpar%ds_pixwin(n), cpar%ds_btheta_file(n))
    allocate(cpar%ds_sample_gain(n), cpar%ds_gain_calib_comp(n), cpar%ds_gain_lmax(n))
    allocate(cpar%ds_gain_lmin(n), cpar%ds_gain_apodmask(n), cpar%ds_gain_fwhm(n))
    allocate(cpar%ds_defaults(n,2))
    do i = 1, n
       call int2string(i, itext)
       call get_parameter(paramfile, 'INCLUDE_BAND'//itext,         par_lgt=cpar%ds_active(i))
       call get_parameter(paramfile, 'BAND_OBS_PERIOD'//itext,      par_int=cpar%ds_period(i))
       call get_parameter(paramfile, 'BAND_LABEL'//itext,           par_string=cpar%ds_label(i))
       call get_parameter(paramfile, 'BAND_POLARIZATION'//itext,    par_lgt=cpar%ds_polarization(i))
       call get_parameter(paramfile, 'BAND_NSIDE'//itext,           par_int=cpar%ds_nside(i))
       call get_parameter(paramfile, 'BAND_LMAX'//itext,            par_int=cpar%ds_lmax(i))
       call get_parameter(paramfile, 'BAND_UNIT'//itext,            par_string=cpar%ds_unit(i))
       call get_parameter(paramfile, 'BAND_NOISE_FORMAT'//itext,    par_string=cpar%ds_noise_format(i))
       call get_parameter(paramfile, 'BAND_MAPFILE'//itext,         par_string=cpar%ds_mapfile(i))
       call get_parameter(paramfile, 'BAND_NOISE_RMS'//itext,       par_string=cpar%ds_noise_rms(i))
       call get_parameter(paramfile, 'BAND_NOISE_UNIFORMIZE_FSKY'//itext, par_dp=cpar%ds_noise_uni_fsky(i))
       call get_parameter(paramfile, 'BAND_MASKFILE'//itext,        par_string=cpar%ds_maskfile(i))
       call get_parameter(paramfile, 'BAND_MASKFILE_CALIB'//itext,  par_string=cpar%ds_maskfile_calib(i))
       call get_parameter(paramfile, 'BAND_BEAMTYPE'//itext,        par_string=cpar%ds_beamtype(i))
       call get_parameter(paramfile, 'BAND_BEAM_B_L_FILE'//itext,   par_string=cpar%ds_blfile(i))
       call get_parameter(paramfile, 'BAND_BEAM_B_PTSRC_FILE'//itext, par_string=cpar%ds_btheta_file(i))
       call get_parameter(paramfile, 'BAND_PIXEL_WINDOW'//itext,    par_string=cpar%ds_pixwin(i))
       call get_parameter(paramfile, 'BAND_SAMP_NOISE_AMP'//itext,  par_lgt=cpar%ds_samp_noiseamp(i))
       call get_parameter(paramfile, 'BAND_BANDPASS_TYPE'//itext,   par_string=cpar%ds_bptype(i))
       call get_parameter(paramfile, 'BAND_NOMINAL_FREQ'//itext,    par_dp=cpar%ds_nu_c(i))
       call get_parameter(paramfile, 'BAND_BANDPASSFILE'//itext,    par_string=cpar%ds_bpfile(i))
       call get_parameter(paramfile, 'BAND_BANDPASS_MODEL'//itext,  par_string=cpar%ds_bpmodel(i))
       call get_parameter(paramfile, 'BAND_SAMP_GAIN'//itext,       par_lgt=cpar%ds_sample_gain(i))
       call get_parameter(paramfile, 'BAND_GAIN_CALIB_COMP'//itext, par_string=cpar%ds_gain_calib_comp(i))
       call get_parameter(paramfile, 'BAND_GAIN_LMAX'//itext,       par_int=cpar%ds_gain_lmax(i))
       call get_parameter(paramfile, 'BAND_GAIN_LMIN'//itext,       par_int=cpar%ds_gain_lmin(i))
       call get_parameter(paramfile, 'BAND_GAIN_APOD_MASK'//itext,  par_string=cpar%ds_gain_apodmask(i))
       call get_parameter(paramfile, 'BAND_GAIN_APOD_FWHM'//itext,  par_string=cpar%ds_gain_fwhm(i))
       call get_parameter(paramfile, 'BAND_DEFAULT_GAIN'//itext,    par_dp=cpar%ds_defaults(i,GAIN))
       call get_parameter(paramfile, 'BAND_DEFAULT_NOISEAMP'//itext,par_dp=cpar%ds_defaults(i,NOISEAMP))
       do j = 1, cpar%num_smooth_scales
          call int2string(j, jtext)          
          call get_parameter(paramfile, 'BAND_NOISE_RMS'//itext//'_SMOOTH'//jtext, &
               & par_string=cpar%ds_noise_rms_smooth(i,j))
       end do
    end do

    ! Convert to proper internal units where necessary
    cpar%ds_nu_c = cpar%ds_nu_c * 1d9                           ! From GHz to Hz
    
  end subroutine read_data_params


  subroutine read_component_params(paramfile, cpar)
    implicit none

    character(len=*),  intent(in)    :: paramfile
    type(comm_params), intent(inout) :: cpar

    integer(i4b)     :: i, n
    character(len=2) :: itext
    
    call get_parameter(paramfile, 'INSTRUMENT_PARAM_FILE', par_string=cpar%cs_inst_parfile)
    call get_parameter(paramfile, 'NUM_SIGNAL_COMPONENTS', par_int=cpar%cs_ncomp_tot)

    n = cpar%cs_ncomp_tot
    allocate(cpar%cs_include(n), cpar%cs_label(n), cpar%cs_type(n), cpar%cs_class(n))
    allocate(cpar%cs_polarization(n), cpar%cs_nside(n), cpar%cs_lmax_amp(n), cpar%cs_lmax_ind(n))
    allocate(cpar%cs_l_apod(n), cpar%cs_output_EB(n), cpar%cs_initHDF(n))
    allocate(cpar%cs_unit(n), cpar%cs_nu_ref(n), cpar%cs_cltype(n), cpar%cs_cl_poltype(n))
    allocate(cpar%cs_clfile(n), cpar%cs_binfile(n), cpar%cs_band_ref(n))
    allocate(cpar%cs_lpivot(n), cpar%cs_mask(n), cpar%cs_fwhm(n), cpar%cs_poltype(MAXPAR,n))
    allocate(cpar%cs_indmask(n), cpar%cs_amp_rms_scale(n))
    allocate(cpar%cs_cl_amp_def(n,3), cpar%cs_cl_beta_def(n,3), cpar%cs_cl_prior(n,2))
    allocate(cpar%cs_input_amp(n), cpar%cs_prior_amp(n), cpar%cs_input_ind(MAXPAR,n))
    allocate(cpar%cs_theta_def(MAXPAR,n), cpar%cs_p_uni(n,2,MAXPAR), cpar%cs_p_gauss(n,2,MAXPAR))
    allocate(cpar%cs_catalog(n), cpar%cs_SED_template(4,n), cpar%cs_cg_scale(n), cpar%cs_cg_samp_group(n))
    allocate(cpar%cs_ptsrc_template(n), cpar%cs_output_ptsrc_beam(n), cpar%cs_min_src_dist(n))
    allocate(cpar%cs_auxpar(MAXAUXPAR,n), cpar%cs_apply_pos_prior(n))
    allocate(cpar%cs_nu_min(n,MAXPAR), cpar%cs_nu_max(n,MAXPAR), cpar%cs_burn_in(n))
    allocate(cpar%cs_smooth_scale(n,MAXPAR))
    do i = 1, n
       call int2string(i, itext)
       call get_parameter(paramfile, 'INCLUDE_COMP'//itext,         par_lgt=cpar%cs_include(i))
       call get_parameter(paramfile, 'COMP_INIT_FROM_HDF'//itext,   par_lgt=cpar%cs_initHDF(i))
       call get_parameter(paramfile, 'COMP_LABEL'//itext,           par_string=cpar%cs_label(i))
       call get_parameter(paramfile, 'COMP_TYPE'//itext,            par_string=cpar%cs_type(i))
       call get_parameter(paramfile, 'COMP_CLASS'//itext,           par_string=cpar%cs_class(i))
       call get_parameter(paramfile, 'COMP_CG_SAMPLE_GROUP'//itext, par_int=cpar%cs_cg_samp_group(i))
       if (trim(cpar%cs_type(i)) == 'md') then
          call get_parameter(paramfile, 'COMP_POLARIZATION'//itext,    par_lgt=cpar%cs_polarization(i))
          call get_parameter(paramfile, 'COMP_MD_DEFINITION_FILE'//itext, par_string=cpar%cs_SED_template(1,i))
       else if (trim(cpar%cs_type(i)) == 'template') then
          call get_parameter(paramfile, 'COMP_TEMPLATE_DEFINITION_FILE'//itext, &
               & par_string=cpar%cs_SED_template(1,i))
       else if (trim(cpar%cs_class(i)) == 'diffuse') then
          call get_parameter(paramfile, 'COMP_POLARIZATION'//itext,    par_lgt=cpar%cs_polarization(i))
          if (cpar%cs_polarization(i)) &
               & call get_parameter(paramfile, 'COMP_OUTPUT_EB_MAP'//itext,        par_lgt=cpar%cs_output_EB(i))
          call get_parameter(paramfile, 'COMP_CG_SCALE'//itext,        par_dp=cpar%cs_cg_scale(i))
          call get_parameter(paramfile, 'COMP_NSIDE'//itext,           par_int=cpar%cs_nside(i))
          call get_parameter(paramfile, 'COMP_LMAX_AMP'//itext,        par_int=cpar%cs_lmax_amp(i))
          call get_parameter(paramfile, 'COMP_L_APOD'//itext,          par_int=cpar%cs_l_apod(i))
          call get_parameter(paramfile, 'COMP_LMAX_IND'//itext,        par_int=cpar%cs_lmax_ind(i))
          call get_parameter(paramfile, 'COMP_UNIT'//itext,            par_string=cpar%cs_unit(i))
          call get_parameter(paramfile, 'COMP_NU_REF'//itext,          par_dp=cpar%cs_nu_ref(i))
          call get_parameter(paramfile, 'COMP_CL_TYPE'//itext,         par_string=cpar%cs_cltype(i))
          call get_parameter(paramfile, 'COMP_INPUT_AMP_MAP'//itext,   par_string=cpar%cs_input_amp(i))
          call get_parameter(paramfile, 'COMP_PRIOR_AMP_MAP'//itext,   par_string=cpar%cs_prior_amp(i))
          call get_parameter(paramfile, 'COMP_OUTPUT_FWHM'//itext,     par_dp=cpar%cs_fwhm(i))
          if (trim(cpar%cs_cltype(i)) == 'binned') then
             call get_parameter(paramfile, 'COMP_CL_BIN_FILE'//itext,     par_string=cpar%cs_binfile(i))
             call get_parameter(paramfile, 'COMP_CL_DEFAULT_FILE'//itext,     par_string=cpar%cs_clfile(i))
          else if (trim(cpar%cs_cltype(i)) == 'power_law' .or. trim(cpar%cs_cltype(i)) == 'exp' .or. trim(cpar%cs_cltype(i))=='gauss') then
             call get_parameter(paramfile, 'COMP_CL_POLTYPE'//itext,      par_int=cpar%cs_cl_poltype(i))
             call get_parameter(paramfile, 'COMP_CL_L_PIVOT'//itext,      par_int=cpar%cs_lpivot(i))
             call get_parameter(paramfile, 'COMP_CL_BETA_PRIOR_MEAN'//itext, par_dp=cpar%cs_cl_prior(i,1))
             call get_parameter(paramfile, 'COMP_CL_BETA_PRIOR_RMS'//itext, par_dp=cpar%cs_cl_prior(i,2))
             call get_parameter(paramfile, 'COMP_CL_DEFAULT_AMP_T'//itext,  par_dp=cpar%cs_cl_amp_def(i,1))
             call get_parameter(paramfile, 'COMP_CL_DEFAULT_BETA_T'//itext,  par_dp=cpar%cs_cl_beta_def(i,1))
             if (cpar%cs_polarization(i)) then
                call get_parameter(paramfile, 'COMP_CL_DEFAULT_AMP_E'//itext,  par_dp=cpar%cs_cl_amp_def(i,2))
                call get_parameter(paramfile, 'COMP_CL_DEFAULT_BETA_E'//itext,  par_dp=cpar%cs_cl_beta_def(i,2))
                call get_parameter(paramfile, 'COMP_CL_DEFAULT_AMP_B'//itext,  par_dp=cpar%cs_cl_amp_def(i,3))
                call get_parameter(paramfile, 'COMP_CL_DEFAULT_BETA_B'//itext,  par_dp=cpar%cs_cl_beta_def(i,3))
             end if
             cpar%cs_cl_amp_def(i,:) = cpar%cs_cl_amp_def(i,:) / cpar%cs_cg_scale(i)**2
          end if
          call get_parameter(paramfile, 'COMP_MASK'//itext,            par_string=cpar%cs_mask(i))
          cpar%cs_indmask(i) = 'fullsky'

          select case (trim(cpar%cs_type(i)))
          case ('power_law')
             call get_parameter(paramfile, 'COMP_BETA_POLTYPE'//itext,  par_int=cpar%cs_poltype(1,i))
             call get_parameter(paramfile, 'COMP_INPUT_BETA_MAP'//itext,        &
                  & par_string=cpar%cs_input_ind(1,i))
             call get_parameter(paramfile, 'COMP_DEFAULT_BETA'//itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_BETA_LOW'//itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_BETA_HIGH'//itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_BETA_MEAN'//itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_BETA_RMS'//itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter(paramfile, 'COMP_INDMASK'//itext,               &
                  & par_string=cpar%cs_indmask(i))
             call get_parameter(paramfile, 'COMP_BETA_SMOOTHING_SCALE'//itext,  &
                  & par_int=cpar%cs_smooth_scale(i,1))
             call get_parameter(paramfile, 'COMP_BETA_NU_MIN'//itext,   par_dp=cpar%cs_nu_min(i,1))
             call get_parameter(paramfile, 'COMP_BETA_NU_MAX'//itext,   par_dp=cpar%cs_nu_max(i,1))
          case ('physdust')
             call get_parameter(paramfile, 'COMP_UMIN_POLTYPE'//itext,  par_int=cpar%cs_poltype(1,i))
             call get_parameter(paramfile, 'COMP_INPUT_UMIN_MAP'//itext,        &
                  & par_string=cpar%cs_input_ind(1,i))
             call get_parameter(paramfile, 'COMP_DEFAULT_UMIN'//itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_UMIN_LOW'//itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_UMIN_HIGH'//itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_UMIN_MEAN'//itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_UMIN_RMS'//itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter(paramfile, 'COMP_UMAX'//itext,  par_dp=cpar%cs_auxpar(1,i))
             call get_parameter(paramfile, 'COMP_GAMMA'//itext,  par_dp=cpar%cs_auxpar(2,i))
             call get_parameter(paramfile, 'COMP_ALPHA'//itext,  par_dp=cpar%cs_auxpar(3,i))
             call get_parameter(paramfile, 'COMP_SIL_AMP1_'//itext,  par_dp=cpar%cs_auxpar(4,i))
             call get_parameter(paramfile, 'COMP_SIL_AMP2_'//itext,  par_dp=cpar%cs_auxpar(5,i))
             call get_parameter(paramfile, 'COMP_CARB_AMP1_'//itext,  par_dp=cpar%cs_auxpar(6,i))
             call get_parameter(paramfile, 'COMP_CARB_AMP2_'//itext,  par_dp=cpar%cs_auxpar(7,i))
             call get_parameter(paramfile, 'COMP_SIL_FILE1_'//itext,  par_string=cpar%cs_SED_template(1,i))
             call get_parameter(paramfile, 'COMP_SIL_FILE2_'//itext,  par_string=cpar%cs_SED_template(2,i))
             call get_parameter(paramfile, 'COMP_CARB_FILE1_'//itext,  par_string=cpar%cs_SED_template(3,i))
             call get_parameter(paramfile, 'COMP_CARB_FILE2_'//itext,  par_string=cpar%cs_SED_template(4,i))
             call get_parameter(paramfile, 'COMP_INDMASK'//itext,         par_string=cpar%cs_indmask(i))
             write(*,*) 'ERROR -- smoothing not yet supported.'
             stop
          case ('spindust')
             call get_parameter(paramfile, 'COMP_NU_P_POLTYPE'//itext,  par_int=cpar%cs_poltype(1,i))
             call get_parameter(paramfile, 'COMP_INPUT_NU_P_MAP'//itext,        &
                  & par_string=cpar%cs_input_ind(1,i))
             call get_parameter(paramfile, 'COMP_DEFAULT_NU_P'//itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_NU_P_LOW'//itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_NU_P_HIGH'//itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_NU_P_MEAN'//itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_NU_P_RMS'//itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter(paramfile, 'COMP_SED_TEMPLATE'//itext,  &
                  & par_string=cpar%cs_SED_template(1,i))
             call get_parameter(paramfile, 'COMP_INDMASK'//itext,         par_string=cpar%cs_indmask(i))
             call get_parameter(paramfile, 'COMP_NU_P_SMOOTHING_SCALE'//itext,  &
                  & par_int=cpar%cs_smooth_scale(i,1))
             call get_parameter(paramfile, 'COMP_NU_P_NU_MIN'//itext,   par_dp=cpar%cs_nu_min(i,1))
             call get_parameter(paramfile, 'COMP_NU_P_NU_MAX'//itext,   par_dp=cpar%cs_nu_max(i,1))
          case ('MBB')
             call get_parameter(paramfile, 'COMP_BETA_POLTYPE'//itext,  par_int=cpar%cs_poltype(1,i))
             call get_parameter(paramfile, 'COMP_INPUT_BETA_MAP'//itext,        &
                  & par_string=cpar%cs_input_ind(1,i))
             call get_parameter(paramfile, 'COMP_DEFAULT_BETA'//itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_BETA_LOW'//itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_BETA_HIGH'//itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_BETA_MEAN'//itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_BETA_RMS'//itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter(paramfile, 'COMP_T_POLTYPE'//itext,  par_int=cpar%cs_poltype(2,i))
             call get_parameter(paramfile, 'COMP_INPUT_T_MAP'//itext,        &
                  & par_string=cpar%cs_input_ind(2,i))
             call get_parameter(paramfile, 'COMP_DEFAULT_T'//itext,          &
                  & par_dp=cpar%cs_theta_def(2,i))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_T_LOW'//itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,2))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_T_HIGH'//itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,2))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_T_MEAN'//itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,2))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_T_RMS'//itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,2))
             call get_parameter(paramfile, 'COMP_INDMASK'//itext,         par_string=cpar%cs_indmask(i))
             call get_parameter(paramfile, 'COMP_BETA_SMOOTHING_SCALE'//itext,  &
                  & par_int=cpar%cs_smooth_scale(i,1))
             call get_parameter(paramfile, 'COMP_T_SMOOTHING_SCALE'//itext,  &
                  & par_int=cpar%cs_smooth_scale(i,2))
             call get_parameter(paramfile, 'COMP_BETA_NU_MIN'//itext,   par_dp=cpar%cs_nu_min(i,1))
             call get_parameter(paramfile, 'COMP_BETA_NU_MAX'//itext,   par_dp=cpar%cs_nu_max(i,1))
             call get_parameter(paramfile, 'COMP_T_NU_MIN'//itext,      par_dp=cpar%cs_nu_min(i,2))
             call get_parameter(paramfile, 'COMP_T_NU_MAX'//itext,      par_dp=cpar%cs_nu_max(i,2))
          case ('freefree')
             call get_parameter(paramfile, 'COMP_EM_POLTYPE'//itext,  par_int=cpar%cs_poltype(1,i))
             call get_parameter(paramfile, 'COMP_INPUT_EM_MAP'//itext,        &
                  & par_string=cpar%cs_input_ind(1,i))
             call get_parameter(paramfile, 'COMP_DEFAULT_EM'//itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_EM_LOW'//itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_EM_HIGH'//itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_EM_MEAN'//itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_EM_RMS'//itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter(paramfile, 'COMP_TE_POLTYPE'//itext,  par_int=cpar%cs_poltype(2,i))
             call get_parameter(paramfile, 'COMP_INPUT_TE_MAP'//itext,        &
                  & par_string=cpar%cs_input_ind(2,i))
             call get_parameter(paramfile, 'COMP_DEFAULT_TE'//itext,          &
                  & par_dp=cpar%cs_theta_def(2,i))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_TE_LOW'//itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,2))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_TE_HIGH'//itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,2))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_TE_MEAN'//itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,2))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_TE_RMS'//itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,2))
             call get_parameter(paramfile, 'COMP_INDMASK'//itext,         par_string=cpar%cs_indmask(i))
             call get_parameter(paramfile, 'COMP_EM_SMOOTHING_SCALE'//itext,  &
                  & par_int=cpar%cs_smooth_scale(i,1))
             call get_parameter(paramfile, 'COMP_T_E_SMOOTHING_SCALE'//itext,  &
                  & par_int=cpar%cs_smooth_scale(i,2))
             call get_parameter(paramfile, 'COMP_EM_NU_MIN'//itext,   par_dp=cpar%cs_nu_min(i,1))
             call get_parameter(paramfile, 'COMP_EM_NU_MAX'//itext,   par_dp=cpar%cs_nu_max(i,1))
             call get_parameter(paramfile, 'COMP_T_E_NU_MIN'//itext,   par_dp=cpar%cs_nu_min(i,2))
             call get_parameter(paramfile, 'COMP_T_E_NU_MAX'//itext,   par_dp=cpar%cs_nu_max(i,2))
          case ('line')
             call get_parameter(paramfile, 'COMP_LINE_TEMPLATE'//itext,  &
                  & par_string=cpar%cs_SED_template(1,i))
             call get_parameter(paramfile, 'COMP_BAND_REF'//itext, &
                  & par_string=cpar%cs_band_ref(i))
             call get_parameter(paramfile, 'COMP_INDMASK'//itext,         par_string=cpar%cs_indmask(i))
          end select

       else if (trim(cpar%cs_class(i)) == 'ptsrc') then
          call get_parameter(paramfile, 'COMP_POLARIZATION'//itext,    par_lgt=cpar%cs_polarization(i))
          call get_parameter(paramfile, 'COMP_CATALOG'//itext,  par_string=cpar%cs_catalog(i))
          call get_parameter(paramfile, 'COMP_PTSRC_TEMPLATE'//itext,  &
               & par_string=cpar%cs_ptsrc_template(i))
          call get_parameter(paramfile, 'COMP_BURN_IN_ON_FIRST_SAMPLE'//itext,  par_lgt=cpar%cs_burn_in(i))
          call get_parameter(paramfile, 'COMP_AMP_RMS_SCALE_FACTOR'//itext,  par_dp=cpar%cs_amp_rms_scale(i))
          call get_parameter(paramfile, 'COMP_OUTPUT_PTSRC_TEMPLATE'//itext,  &
               & par_lgt=cpar%cs_output_ptsrc_beam(i))
          call get_parameter(paramfile, 'COMP_APPLY_POSITIVITY_PRIOR'//itext, par_lgt=cpar%cs_apply_pos_prior(i))
          if (cpar%cs_apply_pos_prior(i)) cpar%cs_cg_samp_group(i) = 0
          call get_parameter(paramfile, 'COMP_MIN_DIST_BETWEEN_SRC'//itext, par_dp=cpar%cs_min_src_dist(i))
          call get_parameter(paramfile, 'COMP_POLTYPE'//itext,  par_int=cpar%cs_poltype(1,i))
          call get_parameter(paramfile, 'COMP_NSIDE'//itext,    par_int=cpar%cs_nside(i))
          call get_parameter(paramfile, 'COMP_NU_REF'//itext,   par_dp=cpar%cs_nu_ref(i))
          call get_parameter(paramfile, 'COMP_CG_SCALE'//itext, par_dp=cpar%cs_cg_scale(i))
          select case (trim(cpar%cs_type(i)))
          case ('radio')
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_ALPHA_LOW'//itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_ALPHA_HIGH'//itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_ALPHA_MEAN'//itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_ALPHA_RMS'//itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))
             call get_parameter(paramfile, 'COMP_DEFAULT_ALPHA'//itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_BETA_LOW'//itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,2))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_BETA_HIGH'//itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,2))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_BETA_MEAN'//itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,2))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_BETA_RMS'//itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,2))
             call get_parameter(paramfile, 'COMP_DEFAULT_BETA'//itext,          &
                  & par_dp=cpar%cs_theta_def(2,i))
             call get_parameter(paramfile, 'COMP_ALPHA_NU_MIN'//itext,  par_dp=cpar%cs_nu_min(i,1))
             call get_parameter(paramfile, 'COMP_ALPHA_NU_MAX'//itext,  par_dp=cpar%cs_nu_max(i,1))
             call get_parameter(paramfile, 'COMP_BETA_NU_MIN'//itext,   par_dp=cpar%cs_nu_min(i,2))
             call get_parameter(paramfile, 'COMP_BETA_NU_MAX'//itext,   par_dp=cpar%cs_nu_max(i,2))
          case ('fir')
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_BETA_LOW'//itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_BETA_HIGH'//itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_BETA_MEAN'//itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_BETA_RMS'//itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,1))             
             call get_parameter(paramfile, 'COMP_DEFAULT_BETA'//itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_T_LOW'//itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,2))
             call get_parameter(paramfile, 'COMP_PRIOR_UNI_T_HIGH'//itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,2))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_T_MEAN'//itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,2))
             call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_T_RMS'//itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,2))
             call get_parameter(paramfile, 'COMP_DEFAULT_T'//itext,          &
                  & par_dp=cpar%cs_theta_def(2,i))             
             call get_parameter(paramfile, 'COMP_BETA_NU_MIN'//itext,   par_dp=cpar%cs_nu_min(i,1))
             call get_parameter(paramfile, 'COMP_BETA_NU_MAX'//itext,   par_dp=cpar%cs_nu_max(i,1))
             call get_parameter(paramfile, 'COMP_T_NU_MIN'//itext,      par_dp=cpar%cs_nu_min(i,2))
             call get_parameter(paramfile, 'COMP_T_NU_MAX'//itext,      par_dp=cpar%cs_nu_max(i,2))
          end select
       end if
       
    end do
    cpar%cs_ncomp           = count(cpar%cs_include)
    cpar%cg_num_samp_groups = maxval(cpar%cs_cg_samp_group)

    ! Convert to proper units
    cpar%cs_nu_ref = 1d9 * cpar%cs_nu_ref
    cpar%cs_nu_min = 1d9 * cpar%cs_nu_min
    cpar%cs_nu_max = 1d9 * cpar%cs_nu_max


  end subroutine read_component_params



  ! ********************************************************
  !                     Utility routines
  ! ********************************************************
  subroutine get_parameter(parfile, parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present, desc)
    implicit none
    character(len=*)           :: parfile, parname
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present
    character(len=*), optional :: desc

    logical(lgt)               :: found
    integer(i4b)               :: unit

    unit = getlun()
    found = .false.
    call get_parameter_arg(parname, par_int, par_char, par_string, par_sp, par_dp, par_lgt, found, desc)
    if(found) then
       if(present(par_present)) par_present = .true.
    else
       call get_parameter_parfile(parfile, parname, par_int, par_char, par_string, par_sp, par_dp, par_lgt, par_present, desc)
    end if
  end subroutine

  subroutine get_parameter_parfile(parfile, parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present, desc)
    implicit none
    character(len=*)           :: parfile, parname
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present
    character(len=*), optional :: desc

    logical(lgt)               :: found
    integer(i4b), parameter    :: maxdepth = 256
    integer(i4b)               :: depth, units(maxdepth), i
    character(len=512)         :: key, value, filenames(maxdepth), line

    depth = 1
    units(depth) = getlun()
    !write(*,*) "Entering file " // trim(parfile)
    filenames(depth) = parfile
    open(units(depth),file=trim(parfile),status="old",err=4)
    do while(depth >= 1)
       read(units(depth),*,end=1) key
       if (key(1:1)=='#') cycle
       backspace(units(depth))

       if (key(1:1)=='@') then
          if(key == '@INCLUDE') then
             ! Recurse into the new file
             read(units(depth),*,end=1) key, value
             !write(*,*) "Entering file " // trim(value)
             depth=depth+1
             units(depth) = getlun()
             filenames(depth) = value
             open(units(depth),file=value,status="old",err=2)
          else
             goto 3
          end if
       else
          read(units(depth),fmt="(a)") line
          call parse_parameter(line, parname, found, par_int, par_char, par_string, par_sp, par_dp, par_lgt)
          if(found) then
             ! Match found, so clean up and return.
             do i = depth, 1, -1; close(units(i)); end do
             if(present(par_present)) par_present = .true.
             return
          end if
       end if
       cycle
       ! We get here if we reached the end of a file. Close it and
       ! return to the file above.
1      close(units(depth))
       !write(*,*) "Exiting file " // filenames(depth)
       depth = depth-1
    end do

    ! ===== Error handling section ======
    ! Case 1: Failed to find matching parameter in file
    if(present(par_present)) then
       par_present = .false.
       return
    else
       write(*,*) "Error: Cannot find " // trim(parname) // &
            & " in " // trim(parfile) // " or included files."
       if(present(desc)) write(*,*) trim(desc)
       stop
    end if

    ! Case 2: Include file error
2   write(*,*) "Error: Cannot open include file '" // trim(value) // "'"
    write(*,*) " in file " // trim(filenames(depth))
    do i = depth-1, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
    do i = depth-1, 1, -1; close(units(i)); end do
    stop

    ! Case 3: Directive error
3   write(*,*) "Error: Unrecognized directive '" // trim(key) //"'"
    write(*,*) " in file " // trim(filenames(depth))
    do i = depth-1, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
    do i = depth, 1, -1; close(units(i)); end do
    stop

    ! Case 4: Top level parameter file unreadable
4   write(*,*) "Error: Cannot open parameter file '" // trim(parfile) // "'"
    stop
  end subroutine

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

  subroutine get_parameter_arr(arr, parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present)
    implicit none
    character(len=*)           :: arr(:)
    character(len=*)           :: parname
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present

    integer(i4b)       :: i
    logical(lgt)       :: found
    do i = 1, size(arr)
       call parse_parameter(arr(i), parname, found, par_int, par_char, par_string, par_sp, par_dp, par_lgt)
       if(found) then
          if(present(par_present)) par_present = .true.
          return
       end if
    end do
    if(present(par_present)) then
       par_present = .false.
    else
       write(*,*) "get_parameter_arr: Fatal error: Cannot find " // trim(parname) // " in argument list!"
       stop
    end if
  end subroutine


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
  end subroutine


  function get_token(string, sep, num, group, allow_empty) result(res)
    implicit none
    character(len=*)           :: string, sep
    character(len=len(string)) :: res
    character(len=*), optional :: group
    logical(lgt),     optional :: allow_empty
    integer(i4b)               :: i, num, ext(2)
    ext = -1
    do i = 1, num; call tokenize(string, sep, ext, group, allow_empty); end do
    res = string(ext(1):ext(2))
  end function

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
  end subroutine

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
  end function

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
  end function

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
  end subroutine

  subroutine validate_params(cpar)
    implicit none
    type(comm_params), intent(inout) :: cpar

    integer(i4b) :: i, j
    logical(lgt) :: exist
    character(len=512) :: datadir, chaindir, filename

    datadir  = trim(cpar%datadir) // '/'
    chaindir = trim(cpar%outdir) // '/'
    
    ! Check that all dataset files exist
    do i = 1, cpar%numband
       call validate_file(trim(datadir)//trim(cpar%ds_mapfile(i)))           ! Map file
       if (trim(cpar%ds_maskfile(i)) /= 'fullsky') &
            & call validate_file(trim(datadir)//trim(cpar%ds_maskfile(i)))   ! Mask file
       if (trim(cpar%ds_noise_format(i)) == 'rms') &
            & call validate_file(trim(datadir)//trim(cpar%ds_noise_rms(i)))  ! RMS file
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
    end do

    ! Instrument data base
    if (trim(cpar%cs_inst_parfile) /= 'none') &
         & call validate_file(trim(datadir)//trim(cpar%cs_inst_parfile))   ! Instrument data base
    if (trim(cpar%ds_sourcemask) /= 'none') &
         & call validate_file(trim(cpar%datadir)//'/'//trim(cpar%ds_sourcemask))   ! Source mask
    if (trim(cpar%ds_procmask) /= 'none') &
         & call validate_file(trim(cpar%datadir)//'/'//trim(cpar%ds_procmask))   ! Source mask

    ! Check component files
    do i = 1, cpar%cs_ncomp_tot
       if (.not. cpar%cs_include(i)) cycle

       if (trim(cpar%cs_type(i)) == 'md') then
          call validate_file(trim(cpar%cs_SED_template(1,i)))
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
          case ('MBB')
             if (trim(cpar%cs_input_ind(1,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(1,i)))
             if (trim(cpar%cs_input_ind(2,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(2,i)))
          case ('freefree')
             if (trim(cpar%cs_input_ind(1,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(1,i)))
             if (trim(cpar%cs_input_ind(2,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(2,i)))             
          case ('line')
             call validate_file(trim(datadir)//trim(cpar%cs_SED_template(1,i)))
          end select

       else if (trim(cpar%cs_class(i)) == 'ptsrc') then
          call validate_file(trim(datadir)//trim(cpar%cs_catalog(i)))
          call validate_file(trim(datadir)//trim(cpar%cs_ptsrc_template(i)), &
               & should_exist=.not. cpar%cs_output_ptsrc_beam(i))
       else if (trim(cpar%cs_class(i)) == 'template') then
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

end module comm_param_mod
