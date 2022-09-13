!================================================================================
! 
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
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
     integer(i4b)       :: num_gibbs_iter, thinning, num_init_chains
     character(len=512) :: chain_status, init_chain_prefix
     real(dp)           :: T_CMB
     character(len=512) :: MJysr_convention
     character(len=512) :: fft_magic_number_file
     character(len=512) :: output_comps
     logical(lgt)       :: only_pol, only_I
     logical(lgt)       :: enable_TOD_analysis
     logical(lgt)       :: enable_TOD_simulations !< start commander in simulation regime
     integer(i4b)       :: tod_freq
     integer(i4b)       :: resamp_hard_gain_prior_nth_iter
     integer(i4b)       :: output_4D_map_nth_iter, output_aux_maps
     logical(lgt)       :: include_tod_zodi
     real(dp),           allocatable, dimension(:)     :: fwhm_smooth
     real(dp),           allocatable, dimension(:)     :: fwhm_postproc_smooth
     integer(i4b),       allocatable, dimension(:)     :: lmax_smooth
     integer(i4b),       allocatable, dimension(:)     :: nside_smooth
     character(len=512), allocatable, dimension(:)     :: pixwin_smooth
     character(len=512), allocatable, dimension(:)     :: init_chain_prefixes
     character(len=512)                                :: sims_output_dir !< simulations directory

     ! alm-sampler
     integer(i4b)       :: almsamp_nsamp, almsamp_nside_chisq_lowres, almsamp_prior_fwhm, almsamp_burnin
     logical(lgt)       :: almsamp_optimize, almsamp_apply_prior, almsamp_pixreg, almsamp_priorsamp_frozen

     ! Output parameters
     character(len=512) :: outdir
     integer(i4b)       :: nside_chisq, nmaps_chisq
     logical(lgt)       :: pol_chisq, output_mixmat, output_residuals, output_chisq, output_cg_eigenvals
     integer(i4b)       :: output_cg_freq
     logical(lgt)       :: output_input_model, ignore_gain_bp, output_debug_seds, output_sig_per_band
     logical(lgt)       :: sample_signal_amplitudes, sample_specind, sample_powspec
     
     ! Numerical parameters
     character(len=512) :: cg_conv_crit, cg_precond
     integer(i4b)       :: cg_lmax_precond, cg_maxiter, cg_num_samp_groups, cg_num_user_samp_groups, cg_miniter, cg_check_conv_freq, cg_samp_group_md
     logical(lgt)       :: cg_init_zero, set_noise_to_mean
     real(dp)           :: cg_tol
     integer(i4b)       :: num_bp_prop
     character(len=512), dimension(MAXSAMPGROUP) :: cg_samp_group
     character(len=512), dimension(MAXSAMPGROUP) :: cg_samp_group_mask
     integer(i4b),       dimension(MAXSAMPGROUP) :: cg_samp_group_maxiter

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
     integer(i4b),       allocatable, dimension(:)   :: ds_noise_lcut
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
     character(len=512), allocatable, dimension(:)   :: ds_tod_jumplist
     character(len=512), allocatable, dimension(:)   :: ds_tod_instfile
     character(len=512), allocatable, dimension(:)   :: ds_tod_dets
     character(len=512), allocatable, dimension(:)   :: ds_tod_bp_init
     character(len=512), allocatable, dimension(:)   :: ds_tod_initHDF
     character(len=512), allocatable, dimension(:)   :: ds_tod_level
     integer(i4b),       allocatable, dimension(:,:) :: ds_tod_scanrange
     integer(i4b),       allocatable, dimension(:)   :: ds_tod_tot_numscan
     integer(i4b),       allocatable, dimension(:)   :: ds_tod_flag
     integer(i4b),       allocatable, dimension(:)   :: ds_tod_halfring
     logical(lgt),       allocatable, dimension(:)   :: ds_tod_orb_abscal
     logical(lgt),       allocatable, dimension(:)   :: ds_tod_subtract_zodi

     ! Component parameters
     character(len=512) :: cs_inst_parfile
     character(len=512) :: cs_init_inst_hdf
     integer(i4b)       :: cs_ncomp, cs_ncomp_tot, cs_local_burn_in
     logical(lgt)       :: cs_output_localsamp_maps
     real(dp)           :: cmb_dipole_prior(3)
     character(len=512) :: cmb_dipole_prior_mask
     logical(lgt),       allocatable, dimension(:)     :: cs_include
     character(len=512), allocatable, dimension(:)     :: cs_initHDF
     character(len=512), allocatable, dimension(:)     :: cs_label
     character(len=512), allocatable, dimension(:)     :: cs_type
     character(len=512), allocatable, dimension(:)     :: cs_class
     logical(lgt),       allocatable, dimension(:)     :: cs_polarization
     real(dp),           allocatable, dimension(:,:)   :: cs_cg_scale
     integer(i4b),       allocatable, dimension(:)     :: cs_nside
     integer(i4b),       allocatable, dimension(:,:)   :: cs_poltype
     integer(i4b),       allocatable, dimension(:)     :: cs_cg_samp_group_maxiter
     character(len=512), allocatable, dimension(:,:,:) :: cs_spec_lnLtype
     character(len=512), allocatable, dimension(:,:,:) :: cs_spec_pixreg
     character(len=512), allocatable, dimension(:,:)   :: cs_spec_pixreg_map
     character(len=512), allocatable, dimension(:,:,:) :: cs_spec_fix_pixreg
     character(len=512), allocatable, dimension(:,:,:) :: cs_spec_pixreg_priors
     character(len=512), allocatable, dimension(:,:)   :: cs_spec_mask
     character(len=512), allocatable, dimension(:,:)   :: cs_spec_nprop
     character(len=512), allocatable, dimension(:,:)   :: cs_spec_proplen
     character(len=512), allocatable, dimension(:,:)   :: cs_spec_mono_mask
     character(len=512), allocatable, dimension(:,:)   :: cs_spec_mono_freeze
     character(len=512), allocatable, dimension(:,:)   :: cs_spec_mono_type
     character(len=512), allocatable, dimension(:,:)   :: cs_almsamp_init
     character(len=512), allocatable, dimension(:,:)   :: cs_pixreg_init_theta
     integer(i4b),       allocatable, dimension(:,:,:) :: cs_spec_nprop_init
     real(dp),           allocatable, dimension(:,:,:) :: cs_spec_proplen_init
     real(dp),           allocatable, dimension(:,:)   :: cs_spec_corr_limit
     real(dp),           allocatable, dimension(:,:,:,:) :: cs_theta_prior
     integer(i4b),       allocatable, dimension(:,:,:) :: cs_spec_uni_nprop
     logical(lgt),       allocatable, dimension(:,:,:) :: cs_spec_samp_nprop
     logical(lgt),       allocatable, dimension(:,:,:) :: cs_spec_samp_proplen
     logical(lgt),       allocatable, dimension(:,:)   :: cs_spec_mono_combined
     logical(lgt),       allocatable, dimension(:,:)   :: cs_spec_corr_convergence
     integer(i4b),       allocatable, dimension(:,:,:) :: cs_spec_npixreg
     integer(i4b),       allocatable, dimension(:)     :: cs_samp_samp_params_niter
     integer(i4b),       allocatable, dimension(:,:,:) :: cs_lmax_ind_pol
     integer(i4b),       allocatable, dimension(:)     :: cs_lmin_amp
     integer(i4b),       allocatable, dimension(:)     :: cs_lmax_amp
     integer(i4b),       allocatable, dimension(:)     :: cs_lmax_amp_prior
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
     real(dp),           allocatable, dimension(:,:)   :: cs_cl_theta_def
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
     character(len=512), allocatable, dimension(:)     :: cs_init_catalog
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

     ! Zodi parameters (IPD model)
     ! Hyper parameters
     integer(14b)       :: zs_gauss_quad_order
     real(dp)           :: zs_los_cut, zs_delta_t
     logical(lgt)       :: zs_use_cloud, zs_use_band1, zs_use_band2, zs_use_band3, zs_use_ring, &
                           zs_use_feature, zs_use_unit_emissivity
     character(len=512) :: zs_freq_correction_type

     ! Component parameters
     real(dp), allocatable, dimension(:, :) :: zs_common ! shape: (n_comps, 6)
     real(dp)                               :: zs_cloud_alpha, zs_cloud_beta, zs_cloud_gamma, zs_cloud_mu
     real(dp)                               :: zs_band1_delta_zeta, zs_band1_v, zs_band1_p, zs_band1_delta_r
     real(dp)                               :: zs_band2_delta_zeta, zs_band2_v, zs_band2_p, zs_band2_delta_r
     real(dp)                               :: zs_band3_delta_zeta, zs_band3_v, zs_band3_p, zs_band3_delta_r
     real(dp)                               :: zs_ring_r, zs_ring_sigma_r, zs_ring_sigma_z
     real(dp)                               :: zs_feature_r, zs_feature_sigma_r, zs_feature_sigma_z, &
                                               zs_feature_theta, zs_feature_sigma_theta
     ! Source parameters
     real(dp), allocatable, dimension(:, :) :: zs_emissivity, zs_albedo, zs_phase_function ! (n_band, n_comp)
     real(dp), allocatable, dimension(:)    :: zs_solar_irradiance ! (n_band)

  end type comm_params

contains

  ! ********************************************************
  !                     Driver routines
  ! ********************************************************
  subroutine read_comm_params(cpar)
    implicit none
    type(hash_tbl_sll) :: htable
    type(comm_params), intent(inout) :: cpar
    
    integer(i4b)       :: paramfile_len, ierr, i, idx
    character(len=512) :: paramfile, paramfile_name
    character(len=512), allocatable, dimension(:) :: paramfile_cache

    call getarg(1, paramfile)
    ! read parameter file once, save to ascii array
    ! Guess at file size, will be resized later if needed

    if (cpar%myid == cpar%root) then
      paramfile_len = 512
      allocate(paramfile_cache(paramfile_len))
      call read_paramfile_to_ascii(paramfile,paramfile_cache,paramfile_len)
    end if
    call mpi_bcast(paramfile_len, 1, MPI_INTEGER, cpar%root, MPI_COMM_WORLD, ierr)
    if (cpar%myid /= cpar%root) then 
      allocate(paramfile_cache(paramfile_len))
    end if
    do i=1,paramfile_len
       call mpi_bcast(paramfile_cache(i), 512, MPI_CHAR, cpar%root, MPI_COMM_WORLD, ierr)
    end do
    !Initialize a hash table
    call init_hash_tbl_sll(htable,tbl_len=10*paramfile_len)
    ! Put the parameter file into the hash table
    call put_ascii_into_hashtable(paramfile_cache,htable)

    ! Read parameters from the hash table
    call read_global_params_hash(htable,cpar)
    call read_data_params_hash(htable,cpar)
    call read_component_params_hash(htable,cpar)
    call read_zodi_params_hash(htable, cpar) 
    !output parameter file to output directory
    if (cpar%myid == cpar%root) then
      idx = index(paramfile, '/', back=.true.) 
      paramfile_name = trim(cpar%outdir)//'/'//paramfile(idx+1:len(paramfile))
      call save_ascii_parameter_file(paramfile_name, paramfile_cache) 
    end if 
    !deallocate ascii cache
    deallocate(paramfile_cache)

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
       cpar%init_chain_prefix = cpar%init_chain_prefixes(mod(cpar%mychain-1,cpar%num_init_chains)+1)

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
    call get_parameter_hashtable(htbl, 'THINNING_FACTOR',           par_int=cpar%thinning)
    call get_parameter_hashtable(htbl, 'CHAIN_STATUS',             par_string=cpar%chain_status)
    !call get_parameter_hashtable(htbl, 'INIT_CHAIN',               par_string=cpar%init_chain_prefix)
    call get_parameter_hashtable(htbl, 'NUM_INIT_CHAINS',          par_int=cpar%num_init_chains)
    allocate(cpar%init_chain_prefixes(cpar%num_init_chains))
    do i = 1, cpar%num_init_chains
       call int2string(i,itext)
       call get_parameter_hashtable(htbl, 'INIT_CHAIN'//itext,     par_string=cpar%init_chain_prefixes(i))
    end do
    call get_parameter_hashtable(htbl, 'SAMPLE_ONLY_POLARIZATION', par_lgt=cpar%only_pol)
    call get_parameter_hashtable(htbl, 'SAMPLE_ONLY_TEMPERATURE', par_lgt=cpar%only_I)
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
    call get_parameter_hashtable(htbl, 'OUTPUT_COMPS_TO_CHAINDIR', par_string=cpar%output_comps)

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
    call get_parameter_hashtable(htbl, 'SAMPLE_POWSPEC',           par_lgt=cpar%sample_powspec)

    call get_parameter_hashtable(htbl, 'NUM_SMOOTHING_SCALES',     par_int=cpar%num_smooth_scales)

    call get_parameter_hashtable(htbl, 'ENABLE_TOD_ANALYSIS',      par_lgt=cpar%enable_TOD_analysis)
    !----------------------------------------------------------------------------------
    ! Commander3 simulations parameters
    call get_parameter_hashtable(htbl, 'ENABLE_TOD_SIMULATIONS',   par_lgt=cpar%enable_TOD_simulations)
    if (cpar%enable_TOD_simulations) cpar%num_gibbs_iter = 1
    call get_parameter_hashtable(htbl, 'SIMS_OUTPUT_DIRECTORY',    par_string=cpar%sims_output_dir)
    !----------------------------------------------------------------------------------

    call get_parameter_hashtable(htbl, 'NUMITER_RESAMPLE_HARD_GAIN_PRIORS', par_int=cpar%resamp_hard_gain_prior_nth_iter)

    if (cpar%enable_TOD_analysis) then
       call get_parameter_hashtable(htbl, 'FFTW3_MAGIC_NUMBERS',   par_string=cpar%fft_magic_number_file)
       call get_parameter_hashtable(htbl, 'TOD_NUM_BP_PROPOSALS_PER_ITER', par_int=cpar%num_bp_prop)
       call get_parameter_hashtable(htbl, 'NUM_GIBBS_STEPS_PER_TOD_SAMPLE', par_int=cpar%tod_freq)
       call get_parameter_hashtable(htbl, 'TOD_OUTPUT_4D_MAP_EVERY_NTH_ITER', par_int=cpar%output_4D_map_nth_iter)
       call get_parameter_hashtable(htbl, 'TOD_OUTPUT_AUXILIARY_MAPS_EVERY_NTH_ITER', par_int=cpar%output_aux_maps)
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
    
    call get_parameter_hashtable(htbl, 'ALMSAMP_NSAMP_ALM',          par_int=cpar%almsamp_nsamp)
    call get_parameter_hashtable(htbl, 'ALMSAMP_BURN_IN',            par_int=cpar%almsamp_burnin)
    call get_parameter_hashtable(htbl, 'ALMSAMP_PRIOR_FWHM',         par_int=cpar%almsamp_prior_fwhm)
    call get_parameter_hashtable(htbl, 'ALMSAMP_NSIDE_CHISQ_LOWRES', par_int=cpar%almsamp_nside_chisq_lowres)
    call get_parameter_hashtable(htbl, 'ALMSAMP_OPTIMIZE_ALM',       par_lgt=cpar%almsamp_optimize)
    call get_parameter_hashtable(htbl, 'ALMSAMP_APPLY_PRIOR',        par_lgt=cpar%almsamp_apply_prior)
    call get_parameter_hashtable(htbl, 'ALMSAMP_PIXREG',             par_lgt=cpar%almsamp_pixreg)
    if (cpar%almsamp_pixreg) &
         & call get_parameter_hashtable(htbl, 'ALMSAMP_PRIORSAMP_FROZEN_REGIONS', par_lgt=cpar%almsamp_priorsamp_frozen)

  end subroutine read_global_params_hash


  subroutine read_data_params_hash(htbl, cpar)
    implicit none

    type(hash_tbl_sll), intent(in) :: htbl
    type(comm_params),  intent(inout) :: cpar

    integer(i4b)     :: i, j, n,len_itext
    character(len=3) :: itext
    character(len=2) :: jtext

    len_itext=len(trim(itext)) !! FIXME
    call get_parameter_hashtable(htbl, 'NUMBAND',             par_int=cpar%numband)
    call get_parameter_hashtable(htbl, 'DATA_DIRECTORY',      par_string=cpar%datadir)
    call get_parameter_hashtable(htbl, 'SOURCE_MASKFILE',     par_string=cpar%ds_sourcemask)
    call get_parameter_hashtable(htbl, 'PROCESSING_MASKFILE', par_string=cpar%ds_procmask)

    n = cpar%numband
    allocate(cpar%ds_active(n), cpar%ds_label(n))
    allocate(cpar%ds_polarization(n), cpar%ds_nside(n), cpar%ds_lmax(n))
    allocate(cpar%ds_unit(n), cpar%ds_noise_format(n), cpar%ds_mapfile(n))
    allocate(cpar%ds_noisefile(n), cpar%ds_maskfile(n), cpar%ds_maskfile_calib(n))
    allocate(cpar%ds_regnoise(n), cpar%ds_noise_lcut(n))
    allocate(cpar%ds_noise_rms_smooth(n,cpar%num_smooth_scales))
    allocate(cpar%ds_samp_noiseamp(n), cpar%ds_noise_uni_fsky(n))
    allocate(cpar%ds_bptype(n), cpar%ds_nu_c(n), cpar%ds_bpfile(n), cpar%ds_bpmodel(n))
    allocate(cpar%ds_period(n), cpar%ds_beamtype(n), cpar%ds_blfile(n))
    allocate(cpar%ds_pixwin(n), cpar%ds_btheta_file(n))
    allocate(cpar%ds_sample_gain(n), cpar%ds_gain_prior(n,2), cpar%ds_gain_calib_comp(n), cpar%ds_gain_lmax(n))
    allocate(cpar%ds_gain_lmin(n), cpar%ds_gain_apodmask(n), cpar%ds_gain_fwhm(n))
    allocate(cpar%ds_defaults(n,2))
    allocate(cpar%ds_component_sensitivity(n))
    allocate(cpar%ds_tod_type(n), cpar%ds_tod_filelist(n), cpar%ds_tod_jumplist(n), cpar%ds_tod_initHDF(n), cpar%ds_tod_level(n))
    allocate(cpar%ds_tod_procmask1(n), cpar%ds_tod_procmask2(n), cpar%ds_tod_bp_init(n))
    allocate(cpar%ds_tod_instfile(n), cpar%ds_tod_dets(n), cpar%ds_tod_scanrange(n,2))
    allocate(cpar%ds_tod_tot_numscan(n), cpar%ds_tod_flag(n), cpar%ds_tod_orb_abscal(n), cpar%ds_tod_halfring(n), cpar%ds_tod_subtract_zodi(n))

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
       if (trim(cpar%ds_noise_format(i)) == 'lcut') then
          call get_parameter_hashtable(htbl, 'BAND_NOISE_LCUT'//itext, len_itext=len_itext, par_int=cpar%ds_noise_lcut(i))
       end if
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
             call get_parameter_hashtable(htbl, 'BAND_TOD_JUMPLIST'//itext, len_itext=len_itext, &
                  & par_string=cpar%ds_tod_jumplist(i))
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
             call get_parameter_hashtable(htbl, 'BAND_TOD_SUBTRACT_ZODI'//itext, len_itext=len_itext, &
                  & par_lgt=cpar%ds_tod_subtract_zodi(i))         
             call get_parameter_hashtable(htbl, 'BAND_TOD_RIMO'//itext, len_itext=len_itext, &
                  & par_string=cpar%ds_tod_instfile(i))
             call get_parameter_hashtable(htbl, 'BAND_TOD_BP_INIT_PROP'//itext, len_itext=len_itext, &
                  & par_string=cpar%ds_tod_bp_init(i))
             call get_parameter_hashtable(htbl, 'BAND_TOD_HALFRING'//itext, len_itext=len_itext, par_int=cpar%ds_tod_halfring(i))
             call get_parameter_hashtable(htbl, "BAND_TOD_LEVEL"//itext, len_itext=len_itext, par_string=cpar%ds_tod_level(i))
             if(cpar%ds_tod_level(i) .ne. 'L1' .and. cpar%ds_tod_level(i) .ne. 'L2') then
               write(*,*) "Unrecognized BAND_TOD_LEVEL"//itext//" parameter", trim(cpar%ds_tod_level(i))
               stop
             end if
          end if
       end if

       do j = 1, cpar%num_smooth_scales
          call int2string(j, jtext)          
          call get_parameter_hashtable(htbl, 'BAND_NOISE_RMS'//itext//'_SMOOTH'//jtext, &
               & par_string=cpar%ds_noise_rms_smooth(i,j))
          if (trim(cpar%ds_noise_rms_smooth(i,j)) == 'native') then
             if (cpar%ds_noise_format(i) == 'QUcov') then
                cycle !we allow this, as residuals are udgraded to nside of QUcov
             else if (cpar%ds_nside(i) /= cpar%nside_smooth(j)) then
                write(*,fmt='(a,i3,a,i2)') "nside of band ",i," doesn't match the nside of smoothing scale ",j
                stop
             end if
          end if
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

    integer(i4b)       :: i, j, k, n, len_itext
    real(dp)           :: amp, lat, lon
    character(len=2)   :: itext
    character(len=512) :: maskfile, tokens(4)
    character(len=512) :: pol_labels(3)
    logical(lgt)       :: bool_flag

    pol_labels(1)='INT'
    pol_labels(2)='POL'
    pol_labels(3)='POL3'

    len_itext=len(trim(itext)) !FIXME!!
    call get_parameter_hashtable(htbl, 'INSTRUMENT_PARAM_FILE', par_string=cpar%cs_inst_parfile)
    call get_parameter_hashtable(htbl, 'INIT_INSTRUMENT_FROM_HDF', par_string=cpar%cs_init_inst_hdf)
    call get_parameter_hashtable(htbl, 'NUM_SIGNAL_COMPONENTS', par_int=cpar%cs_ncomp_tot)
    call get_parameter_hashtable(htbl, 'NUM_CG_SAMPLING_GROUPS', par_int=cpar%cg_num_user_samp_groups)

    do i = 1, cpar%cg_num_user_samp_groups
       call int2string(i, itext)
       call get_parameter_hashtable(htbl, 'CG_SAMPLING_GROUP'//itext, par_string=cpar%cg_samp_group(i))
       call get_parameter_hashtable(htbl, 'CG_SAMPLING_GROUP_MASK'//itext, par_string=cpar%cg_samp_group_mask(i))
       call get_parameter_hashtable(htbl, 'CG_SAMPLING_GROUP_MAXITER'//itext, par_int=cpar%cg_samp_group_maxiter(i))
    end do
    call get_parameter_hashtable(htbl, 'LOCALSAMP_BURN_IN', par_int=cpar%cs_local_burn_in)
    call get_parameter_hashtable(htbl, 'LOCALSAMP_OUTPUT_MAPS', par_lgt=cpar%cs_output_localsamp_maps)


    n = cpar%cs_ncomp_tot
    allocate(cpar%cs_include(n), cpar%cs_label(n), cpar%cs_type(n), cpar%cs_class(n))
    allocate(cpar%cs_spec_lnLtype(3,MAXPAR,n))
    allocate(cpar%cs_pixreg_init_theta(MAXPAR,n))
    allocate(cpar%cs_almsamp_init(MAXPAR,n),cpar%cs_theta_prior(2,3,MAXPAR,n))
    allocate(cpar%cs_spec_pixreg(3,MAXPAR,n),cpar%cs_spec_mask(MAXPAR,n))
    allocate(cpar%cs_spec_nprop(MAXPAR,n),cpar%cs_spec_uni_nprop(2,MAXPAR,n))
    allocate(cpar%cs_spec_proplen(MAXPAR,n))
    allocate(cpar%cs_spec_nprop_init(3,MAXPAR,n),cpar%cs_spec_proplen_init(3,MAXPAR,n))
    allocate(cpar%cs_spec_samp_nprop(3,MAXPAR,n),cpar%cs_spec_samp_proplen(3,MAXPAR,n))
    allocate(cpar%cs_spec_npixreg(3,MAXPAR,n),cpar%cs_spec_pixreg_map(MAXPAR,n))
    allocate(cpar%cs_spec_fix_pixreg(3,MAXPAR,n))
    allocate(cpar%cs_spec_pixreg_priors(3,MAXPAR,n))
    allocate(cpar%cs_lmax_ind(n), cpar%cs_lmax_ind_pol(3,MAXPAR,n))
    allocate(cpar%cs_polarization(n), cpar%cs_nside(n), cpar%cs_lmax_amp(n), cpar%cs_lmax_amp_prior(n))
    allocate(cpar%cs_l_apod(n), cpar%cs_output_EB(n), cpar%cs_initHDF(n), cpar%cs_lmin_amp(n))
    allocate(cpar%cs_unit(n), cpar%cs_nu_ref(n,3), cpar%cs_cltype(n), cpar%cs_cl_poltype(n))
    allocate(cpar%cs_clfile(n), cpar%cs_binfile(n), cpar%cs_band_ref(n))
    allocate(cpar%cs_lpivot(n), cpar%cs_mask(n), cpar%cs_mono_prior(n), cpar%cs_fwhm(n), cpar%cs_poltype(MAXPAR,n))
    allocate(cpar%cs_latmask(n), cpar%cs_defmask(n), cpar%cs_cg_samp_group_maxiter(n))
    allocate(cpar%cs_indmask(n), cpar%cs_amp_rms_scale(n))
    allocate(cpar%cs_cl_amp_def(n,3), cpar%cs_cl_beta_def(n,3), cpar%cs_cl_theta_def(n,3), cpar%cs_cl_prior(n,2))
    allocate(cpar%cs_input_amp(n), cpar%cs_prior_amp(n), cpar%cs_input_ind(MAXPAR,n))
    allocate(cpar%cs_theta_def(MAXPAR,n), cpar%cs_p_uni(n,2,MAXPAR), cpar%cs_p_gauss(n,2,MAXPAR))
    allocate(cpar%cs_catalog(n), cpar%cs_init_catalog(n), cpar%cs_SED_template(4,n), cpar%cs_cg_scale(3,n))
    allocate(cpar%cs_ptsrc_template(n), cpar%cs_output_ptsrc_beam(n), cpar%cs_min_src_dist(n))
    allocate(cpar%cs_auxpar(MAXAUXPAR,n), cpar%cs_apply_pos_prior(n))
    allocate(cpar%cs_nu_min(n,MAXPAR), cpar%cs_nu_max(n,MAXPAR), cpar%cs_burn_in(n))
    allocate(cpar%cs_smooth_scale(n,MAXPAR), cpar%cs_apply_jeffreys(n))
    allocate(cpar%cs_spec_mono_combined(n,MAXPAR),cpar%cs_spec_mono_mask(n,MAXPAR),cpar%cs_spec_mono_type(n,MAXPAR))
    allocate(cpar%cs_spec_corr_convergence(MAXPAR,n),cpar%cs_spec_corr_limit(MAXPAR,n))
    allocate(cpar%cs_spec_mono_freeze(n,MAXPAR))
    cpar%cs_spec_mono_combined=.false. !by default
    cpar%cs_spec_corr_convergence=.false. !by default

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
          call get_parameter_hashtable(htbl, 'COMP_CG_SCALE_T'//itext, len_itext=len_itext,        par_dp=cpar%cs_cg_scale(1,i))
          call get_parameter_hashtable(htbl, 'COMP_CG_SCALE_P'//itext, len_itext=len_itext,        par_dp=cpar%cs_cg_scale(2,i))
          cpar%cs_cg_scale(3,i)=cpar%cs_cg_scale(2,i)
          if (.not. trim(cpar%cs_type(i)) == 'cmb') call get_parameter_hashtable(htbl, & 
               &'COMP_CG_SAMP_GROUP_MAXITER'//itext, len_itext=len_itext, par_int=cpar%cs_cg_samp_group_maxiter(i))
          !call get_parameter_hashtable(htbl, 'COMP_CG_SAMP_GROUP_MAXITER'//itext, len_itext=len_itext,        par_int=cpar%cs_cg_samp_group_maxiter(i)) !put it in the components with npar > 0
          call get_parameter_hashtable(htbl, 'COMP_NSIDE'//itext, len_itext=len_itext,           par_int=cpar%cs_nside(i))
          call get_parameter_hashtable(htbl, 'COMP_LMIN_AMP'//itext, len_itext=len_itext,        par_int=cpar%cs_lmin_amp(i))
          call get_parameter_hashtable(htbl, 'COMP_LMAX_AMP'//itext, len_itext=len_itext,        par_int=cpar%cs_lmax_amp(i))
          call get_parameter_hashtable(htbl, 'COMP_L_APOD'//itext, len_itext=len_itext,          par_int=cpar%cs_l_apod(i))
          !call get_parameter_hashtable(htbl, 'COMP_LMAX_IND'//itext, len_itext=len_itext,        par_int=cpar%cs_lmax_ind(i)) ! to be fased out
          call get_parameter_hashtable(htbl, 'COMP_UNIT'//itext, len_itext=len_itext,            par_string=cpar%cs_unit(i))
          call get_parameter_hashtable(htbl, 'COMP_NU_REF_T'//itext, len_itext=len_itext,          par_dp=cpar%cs_nu_ref(i,1))
          call get_parameter_hashtable(htbl, 'COMP_NU_REF_P'//itext, len_itext=len_itext,          par_dp=cpar%cs_nu_ref(i,2))
          cpar%cs_nu_ref(i,3) = cpar%cs_nu_ref(i,2)
          call get_parameter_hashtable(htbl, 'COMP_CL_TYPE'//itext, len_itext=len_itext,         par_string=cpar%cs_cltype(i))
          call get_parameter_hashtable(htbl, 'COMP_INPUT_AMP_MAP'//itext, len_itext=len_itext,   par_string=cpar%cs_input_amp(i))
          call get_parameter_hashtable(htbl, 'COMP_PRIOR_AMP_MAP'//itext, len_itext=len_itext,   par_string=cpar%cs_prior_amp(i))
          if (trim(cpar%cs_prior_amp(i)) /= 'none') then
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_AMP_LMAX'//itext, len_itext=len_itext,        par_int=cpar%cs_lmax_amp_prior(i))
          else
             cpar%cs_lmax_amp_prior(i) = -1
          end if
          call get_parameter_hashtable(htbl, 'COMP_OUTPUT_FWHM'//itext, len_itext=len_itext,     par_dp=cpar%cs_fwhm(i))

          if (trim(cpar%cs_cltype(i)) == 'binned') then
             call get_parameter_hashtable(htbl, 'COMP_CL_BIN_FILE'//itext, len_itext=len_itext,     par_string=cpar%cs_binfile(i))
             call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_FILE'//itext, len_itext=len_itext, par_string=cpar%cs_clfile(i))
          else if (trim(cpar%cs_cltype(i)) == 'power_law' .or. &
               & trim(cpar%cs_cltype(i)) == 'exp' .or. trim(cpar%cs_cltype(i))=='gauss' .or. trim(cpar%cs_cltype(i))=='power_law_gauss') then
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
             if (trim(cpar%cs_cltype(i))=='power_law_gauss') then
                call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_THETA_T'//itext, len_itext=len_itext, &
                     & par_dp=cpar%cs_cl_theta_def(i,1))
             end if
             if (cpar%cs_polarization(i)) then
                call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_AMP_E'//itext, len_itext=len_itext, &
                     & par_dp=cpar%cs_cl_amp_def(i,2))
                call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_BETA_E'//itext, len_itext=len_itext,&
                     & par_dp=cpar%cs_cl_beta_def(i,2))
                call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_AMP_B'//itext, len_itext=len_itext, &
                     & par_dp=cpar%cs_cl_amp_def(i,3))
                call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_BETA_B'//itext, len_itext=len_itext, &  
                & par_dp=cpar%cs_cl_beta_def(i,3))
                if (trim(cpar%cs_cltype(i))=='power_law_gauss') then
                   call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_THETA_E'//itext, len_itext=len_itext, &
                     & par_dp=cpar%cs_cl_theta_def(i,2))
                   call get_parameter_hashtable(htbl, 'COMP_CL_DEFAULT_THETA_B'//itext, len_itext=len_itext, &
                        & par_dp=cpar%cs_cl_theta_def(i,3))
                end if
             end if
             cpar%cs_cl_amp_def(i,:) = cpar%cs_cl_amp_def(i,:) / cpar%cs_cg_scale(:,i)**2
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
             call get_parameter_hashtable(htbl, 'COMP_BETA_POLTYPE'//itext, len_itext=len_itext,  &
                  & par_int=cpar%cs_poltype(1,i))
             k = cpar%cs_poltype(1,i)
             if (.not. cpar%cs_polarization(i)) k = 1 
             do j = 1,k
                call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_LMAX'//itext, &
                     & len_itext=len_itext,        par_int=cpar%cs_lmax_ind_pol(j,1,i))
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0) then
                   call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_LNLTYPE'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_lnLtype(j,1,i))
                   if (trim(cpar%cs_spec_lnLtype(j,1,i)) == 'prior') then
                      call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_PRIOR_MEAN'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(1,j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_PRIOR_RMS'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(2,j,1,i))
                   else
                      call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_SAMPLE_NPROP'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_nprop(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_SAMPLE_PROPLEN'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_proplen(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_NPROP_INIT'//itext, &
                           & len_itext=len_itext, par_int=cpar%cs_spec_nprop_init(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_PROPLEN_INIT'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_spec_proplen_init(j,1,i))
                   end if
                end if
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_NUM_PIXREG'//itext, &
                        & len_itext=len_itext, par_int=cpar%cs_spec_npixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_FIX_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_fix_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_PIXREG_PRIORS'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_priors(j,1,i))
                end if
             end do
             do j = 1,k
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_BETA_PIXREG_MAP'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_map(1,i))
                   exit
                end if
             end do
             bool_flag=.false.
             do j = 1,k
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0 ) bool_flag=.true.
             end do
             if (bool_flag .or. cpar%almsamp_pixreg) &
                  & call get_parameter_hashtable(htbl, 'COMP_BETA_PIXREG_INITVALUE_MAP'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_pixreg_init_theta(1,i))
             if (any(cpar%cs_lmax_ind_pol(:k,1,i) >= 0)) &
                  & call get_parameter_hashtable(htbl, 'COMP_BETA_ALMSAMP_INIT'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_almsamp_init(1,i))
             if (bool_flag) then
                call get_parameter_hashtable(htbl, 'COMP_BETA_UNI_NPROP_LOW'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(1,1,i))
                call get_parameter_hashtable(htbl, 'COMP_BETA_UNI_NPROP_HIGH'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(2,1,i))
                call get_parameter_hashtable(htbl, 'COMP_BETA_MASK'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_mask(1,i))
                call get_parameter_hashtable(htbl, 'COMP_BETA_NPROP'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_nprop(1,i))
                call get_parameter_hashtable(htbl, 'COMP_BETA_PROPLEN'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_proplen(1,i))
                call get_parameter_hashtable(htbl, 'COMP_BETA_CORRELATION_CONVERGENCE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_corr_convergence(1,i))
                if (cpar%cs_spec_corr_convergence(1,i))  call get_parameter_hashtable(htbl, &
                     & 'COMP_BETA_CORRELATION_CONVERGENCE_LIMIT'//itext, &
                     & len_itext=len_itext, par_dp=cpar%cs_spec_corr_limit(1,i))
                call get_parameter_hashtable(htbl, 'COMP_BETA_COMBINED_MONOPOLE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_mono_combined(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_BETA_COMBINED_MONOPOLE_MASK'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_mask(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_BETA_COMBINED_MONOPOLE_TYPE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_type(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_BETA_COMBINED_MONOPOLE_FREEZE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_freeze(i,1))
             end if
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

             do j=1,1
                if (cpar%cs_smooth_scale(i,j) > cpar%num_smooth_scales) then
                   write(*,fmt='(a,i2,a,i2,a,i2,a,i2)') 'Smoothing scale ',cpar%cs_smooth_scale(i,j), &
                        & ' for index nr. ',j,' in component nr. ', i,' is larger than the number of smoothing scales: ', &
                        & cpar%num_smooth_scales
                   stop
                end if
             end do

          case ('physdust')
             call get_parameter_hashtable(htbl, 'COMP_UMIN_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
             k = cpar%cs_poltype(1,i)
             if (.not. cpar%cs_polarization(i)) k = 1 
             do j = 1,k
                call get_parameter_hashtable(htbl, 'COMP_UMIN_'//trim(pol_labels(j))//'_LMAX'//itext, &
                     & len_itext=len_itext,        par_int=cpar%cs_lmax_ind_pol(j,1,i))
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0) then
                   call get_parameter_hashtable(htbl, 'COMP_UMIN_'//trim(pol_labels(j))//'_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_UMIN_'//trim(pol_labels(j))//'_LNLTYPE'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_lnLtype(j,1,i))
                   if (trim(cpar%cs_spec_lnLtype(j,1,i)) == 'prior') then
                      call get_parameter_hashtable(htbl, 'COMP_UMIN_'//trim(pol_labels(j))//'_PRIOR_MEAN'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(1,j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_UMIN_'//trim(pol_labels(j))//'_PRIOR_RMS'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(2,j,1,i))
                   else
                      call get_parameter_hashtable(htbl, 'COMP_UMIN_'//trim(pol_labels(j))//'_SAMPLE_NPROP'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_nprop(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_UMIN_'//trim(pol_labels(j))//'_SAMPLE_PROPLEN'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_proplen(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_UMIN_'//trim(pol_labels(j))//'_NPROP_INIT'//itext, &
                           & len_itext=len_itext, par_int=cpar%cs_spec_nprop_init(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_UMIN_'//trim(pol_labels(j))//'_PROPLEN_INIT'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_spec_proplen_init(j,1,i))
                   end if
                end if
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_UMIN_'//trim(pol_labels(j))//'_NUM_PIXREG'//itext, &
                        & len_itext=len_itext, par_int=cpar%cs_spec_npixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_UMIN_'//trim(pol_labels(j))//'_FIX_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_fix_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_UMIN_'//trim(pol_labels(j))//'_PIXREG_PRIORS'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_priors(j,1,i))
                end if
             end do
             do j = 1,k
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_UMIN_PIXREG_MAP'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_map(1,i))
                   exit
                end if
             end do
             bool_flag=.false.
             do j = 1,k
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0 ) bool_flag=.true.
             end do
             if (bool_flag .or. cpar%almsamp_pixreg) &
                  & call get_parameter_hashtable(htbl, 'COMP_UMIN_PIXREG_INITVALUE_MAP'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_pixreg_init_theta(1,i))
             if (any(cpar%cs_lmax_ind_pol(:k,1,i) >= 0)) &
                  & call get_parameter_hashtable(htbl, 'COMP_UMIN_ALMSAMP_INIT'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_almsamp_init(1,i))
             if (bool_flag) then
                call get_parameter_hashtable(htbl, 'COMP_UMIN_UNI_NPROP_LOW'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(1,1,i))
                call get_parameter_hashtable(htbl, 'COMP_UMIN_UNI_NPROP_HIGH'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(2,1,i))
                call get_parameter_hashtable(htbl, 'COMP_UMIN_MASK'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_mask(1,i))
                call get_parameter_hashtable(htbl, 'COMP_UMIN_NPROP'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_nprop(1,i))
                call get_parameter_hashtable(htbl, 'COMP_UMIN_PROPLEN'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_proplen(1,i))
                call get_parameter_hashtable(htbl, 'COMP_UMIN_CORRELATION_CONVERGENCE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_corr_convergence(1,i))
                if (cpar%cs_spec_corr_convergence(1,i))  call get_parameter_hashtable(htbl, &
                     & 'COMP_UNIM_CORRELATION_CONVERGENCE_LIMIT'//itext, &
                     & len_itext=len_itext, par_dp=cpar%cs_spec_corr_limit(1,i))
                call get_parameter_hashtable(htbl, 'COMP_UMIN_COMBINED_MONOPOLE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_mono_combined(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_UMIN_COMBINED_MONOPOLE_MASK'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_mask(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_UMIN_COMBINED_MONOPOLE_TYPE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_type(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_UMIN_COMBINED_MONOPOLE_FREEZE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_freeze(i,1))
             end if
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

             write(*,*) 'ERROR -- smoothing not yet supported for physdust.'
             stop

          case ('spindust')
             call get_parameter_hashtable(htbl, 'COMP_NU_P_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
             k = cpar%cs_poltype(1,i)
             if (.not. cpar%cs_polarization(i)) k = 1 
             do j = 1,k
                call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_LMAX'//itext, &
                     & len_itext=len_itext,        par_int=cpar%cs_lmax_ind_pol(j,1,i))
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0) then
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_LNLTYPE'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_lnLtype(j,1,i))
                   if (trim(cpar%cs_spec_lnLtype(j,1,i)) == 'prior') then
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PRIOR_MEAN'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(1,j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PRIOR_RMS'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(2,j,1,i))
                   else
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_SAMPLE_NPROP'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_nprop(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_SAMPLE_PROPLEN'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_proplen(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_NPROP_INIT'//itext, &
                           & len_itext=len_itext, par_int=cpar%cs_spec_nprop_init(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PROPLEN_INIT'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_spec_proplen_init(j,1,i))
                   end if
                end if
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_NUM_PIXREG'//itext, &
                        & len_itext=len_itext, par_int=cpar%cs_spec_npixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_FIX_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_fix_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PIXREG_PRIORS'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_priors(j,1,i))
                end if
             end do
             do j = 1,k
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_PIXREG_MAP'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_map(1,i))
                   exit
                end if
             end do
             bool_flag=.false.
             do j = 1,k
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0 ) bool_flag=.true.
             end do
             if (bool_flag .or. cpar%almsamp_pixreg) &
                  & call get_parameter_hashtable(htbl, 'COMP_NU_P_PIXREG_INITVALUE_MAP'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_pixreg_init_theta(1,i))
             if (any(cpar%cs_lmax_ind_pol(:k,1,i) >= 0)) &
                  & call get_parameter_hashtable(htbl, 'COMP_NU_P_ALMSAMP_INIT'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_almsamp_init(1,i))
             if (bool_flag) then
                call get_parameter_hashtable(htbl, 'COMP_NU_P_UNI_NPROP_LOW'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(1,1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_UNI_NPROP_HIGH'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(2,1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_MASK'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_mask(1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_NPROP'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_nprop(1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_PROPLEN'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_proplen(1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_CORRELATION_CONVERGENCE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_corr_convergence(1,i))
                if (cpar%cs_spec_corr_convergence(1,i))  call get_parameter_hashtable(htbl, &
                     & 'COMP_NU_P_CORRELATION_CONVERGENCE_LIMIT'//itext, &
                     & len_itext=len_itext, par_dp=cpar%cs_spec_corr_limit(1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_COMBINED_MONOPOLE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_mono_combined(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_NU_P_COMBINED_MONOPOLE_MASK'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_mask(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_NU_P_COMBINED_MONOPOLE_TYPE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_type(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_NU_P_COMBINED_MONOPOLE_FREEZE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_freeze(i,1))
             end if
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
             do j=1,1
                if (cpar%cs_smooth_scale(i,1) > cpar%num_smooth_scales) then
                   write(*,fmt='(a,i2,a,i2,a,i2,a,i2)') 'Smoothing scale ',cpar%cs_smooth_scale(i,j), &
                        & ' for index nr. ',j,' in component nr. ', i,' is larger than the number of smoothing scales: ', &
                        & cpar%num_smooth_scales
                   stop
                end if
             end do
          ! Spinning dust 2 parameter model
          case ('spindust2')
             call get_parameter_hashtable(htbl, 'COMP_NU_P_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
             k = cpar%cs_poltype(1,i)
             if (.not. cpar%cs_polarization(i)) k = 1 
             do j = 1,k
                call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_LMAX'//itext, &
                     & len_itext=len_itext,        par_int=cpar%cs_lmax_ind_pol(j,1,i))
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0) then
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_LNLTYPE'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_lnLtype(j,1,i))
                   if (trim(cpar%cs_spec_lnLtype(j,1,i)) == 'prior') then
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PRIOR_MEAN'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(1,j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PRIOR_RMS'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(2,j,1,i))
                   else
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_SAMPLE_NPROP'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_nprop(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_SAMPLE_PROPLEN'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_proplen(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_NPROP_INIT'//itext, &
                           & len_itext=len_itext, par_int=cpar%cs_spec_nprop_init(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PROPLEN_INIT'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_spec_proplen_init(j,1,i))
                   end if
                end if
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_NUM_PIXREG'//itext, &
                        & len_itext=len_itext, par_int=cpar%cs_spec_npixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_FIX_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_fix_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PIXREG_PRIORS'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_priors(j,1,i))
                end if
             end do
             do j = 1,k
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_PIXREG_MAP'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_map(1,i))
                   exit
                end if
             end do
             bool_flag=.false.
             do j = 1,k
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0 ) bool_flag=.true.
             end do
             if (bool_flag .or. cpar%almsamp_pixreg) &
                  & call get_parameter_hashtable(htbl, 'COMP_NU_P_PIXREG_INITVALUE_MAP'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_pixreg_init_theta(1,i))
             if (any(cpar%cs_lmax_ind_pol(:k,1,i) >= 0)) &
                  & call get_parameter_hashtable(htbl, 'COMP_NU_P_ALMSAMP_INIT'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_almsamp_init(1,i))
             if (bool_flag) then
                call get_parameter_hashtable(htbl, 'COMP_NU_P_UNI_NPROP_LOW'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(1,1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_UNI_NPROP_HIGH'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(2,1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_MASK'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_mask(1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_NPROP'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_nprop(1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_PROPLEN'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_proplen(1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_CORRELATION_CONVERGENCE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_corr_convergence(1,i))
                if (cpar%cs_spec_corr_convergence(1,i))  call get_parameter_hashtable(htbl, &
                     & 'COMP_NU_P_CORRELATION_CONVERGENCE_LIMIT'//itext, &
                     & len_itext=len_itext, par_dp=cpar%cs_spec_corr_limit(1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_COMBINED_MONOPOLE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_mono_combined(i,2))
                if (cpar%cs_spec_mono_combined(i,2)) call get_parameter_hashtable(htbl, &
                     & 'COMP_NU_P_COMBINED_MONOPOLE_MASK'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_mask(i,2))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_NU_P_COMBINED_MONOPOLE_TYPE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_type(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_NU_P_COMBINED_MONOPOLE_FREEZE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_freeze(i,1))
             end if
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
             k = cpar%cs_poltype(2,i)
             if (.not. cpar%cs_polarization(i)) k = 1 
             do j = 1,k
                call get_parameter_hashtable(htbl, 'COMP_ALPHA_'//trim(pol_labels(j))//'_LMAX'//itext, &
                     & len_itext=len_itext,        par_int=cpar%cs_lmax_ind_pol(j,2,i))
                if (cpar%cs_lmax_ind_pol(j,2,i) < 0) then
                   call get_parameter_hashtable(htbl, 'COMP_ALPHA_'//trim(pol_labels(j))//'_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg(j,2,i))
                   call get_parameter_hashtable(htbl, 'COMP_ALPHA_'//trim(pol_labels(j))//'_LNLTYPE'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_lnLtype(j,2,i))
                   if (trim(cpar%cs_spec_lnLtype(j,2,i)) == 'prior') then
                      call get_parameter_hashtable(htbl, 'COMP_ALPHA_'//trim(pol_labels(j))//'_PRIOR_MEAN'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(1,j,2,i))
                      call get_parameter_hashtable(htbl, 'COMP_ALPHA_'//trim(pol_labels(j))//'_PRIOR_RMS'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(2,j,2,i))
                   else
                      call get_parameter_hashtable(htbl, 'COMP_ALPHA_'//trim(pol_labels(j))//'_SAMPLE_NPROP'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_nprop(j,2,i))
                      call get_parameter_hashtable(htbl, 'COMP_ALPHA_'//trim(pol_labels(j))//'_SAMPLE_PROPLEN'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_proplen(j,2,i))
                      call get_parameter_hashtable(htbl, 'COMP_ALPHA_'//trim(pol_labels(j))//'_NPROP_INIT'//itext, &
                           & len_itext=len_itext, par_int=cpar%cs_spec_nprop_init(j,2,i))
                      call get_parameter_hashtable(htbl, 'COMP_ALPHA_'//trim(pol_labels(j))//'_PROPLEN_INIT'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_spec_proplen_init(j,2,i))
                   end if
                end if
                if (trim(cpar%cs_spec_pixreg(j,2,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_ALPHA_'//trim(pol_labels(j))//'_NUM_PIXREG'//itext, &
                        & len_itext=len_itext, par_int=cpar%cs_spec_npixreg(j,2,i))
                   call get_parameter_hashtable(htbl, 'COMP_ALPHA_'//trim(pol_labels(j))//'_FIX_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_fix_pixreg(j,2,i))
                   call get_parameter_hashtable(htbl, 'COMP_ALPHA_'//trim(pol_labels(j))//'_PIXREG_PRIORS'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_priors(j,2,i))
                end if
             end do
             do j = 1,k
                if (trim(cpar%cs_spec_pixreg(j,2,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_ALPHA_PIXREG_MAP'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_map(2,i))
                   exit
                end if
             end do
             bool_flag=.false.
             do j = 1,k
                if (cpar%cs_lmax_ind_pol(j,2,i) < 0 ) bool_flag=.true.
             end do
             if (bool_flag .or. cpar%almsamp_pixreg) &
                  & call get_parameter_hashtable(htbl, 'COMP_ALPHA_PIXREG_INITVALUE_MAP'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_pixreg_init_theta(2,i))
             if (any(cpar%cs_lmax_ind_pol(:k,2,i) >= 0)) &
                  & call get_parameter_hashtable(htbl, 'COMP_ALPHA_ALMSAMP_INIT'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_almsamp_init(2,i))
             if (bool_flag) then
                call get_parameter_hashtable(htbl, 'COMP_ALPHA_UNI_NPROP_LOW'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(1,2,i))
                call get_parameter_hashtable(htbl, 'COMP_ALPHA_UNI_NPROP_HIGH'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(2,2,i))
                call get_parameter_hashtable(htbl, 'COMP_ALPHA_MASK'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_mask(2,i))
                call get_parameter_hashtable(htbl, 'COMP_ALPHA_NPROP'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_nprop(2,i))
                call get_parameter_hashtable(htbl, 'COMP_ALPHA_PROPLEN'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_proplen(2,i))
                call get_parameter_hashtable(htbl, 'COMP_ALPHA_CORRELATION_CONVERGENCE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_corr_convergence(2,i))
                if (cpar%cs_spec_corr_convergence(2,i))  call get_parameter_hashtable(htbl, &
                     & 'COMP_ALPHA_CORRELATION_CONVERGENCE_LIMIT'//itext, &
                     & len_itext=len_itext, par_dp=cpar%cs_spec_corr_limit(2,i))
                call get_parameter_hashtable(htbl, 'COMP_ALPHA_COMBINED_MONOPOLE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_mono_combined(i,2))
                if (cpar%cs_spec_mono_combined(i,2)) call get_parameter_hashtable(htbl, &
                     & 'COMP_ALPHA_COMBINED_MONOPOLE_MASK'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_mask(i,2))
                if (cpar%cs_spec_mono_combined(i,2)) call get_parameter_hashtable(htbl, &
                     & 'COMP_ALPHA_COMBINED_MONOPOLE_TYPE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_type(i,2))
                if (cpar%cs_spec_mono_combined(i,2)) call get_parameter_hashtable(htbl, &
                     & 'COMP_ALPHA_COMBINED_MONOPOLE_FREEZE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_freeze(i,2))
             end if
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
             do j=1,2
                if (cpar%cs_smooth_scale(i,1) > cpar%num_smooth_scales) then
                   write(*,fmt='(a,i2,a,i2,a,i2,a,i2)') 'Smoothing scale ',cpar%cs_smooth_scale(i,j), &
                        & ' for index nr. ',j,' in component nr. ', i,' is larger than the number of smoothing scales: ', &
                        & cpar%num_smooth_scales
                   stop
                end if
             end do

          case ('lognormal')
             call get_parameter_hashtable(htbl, 'COMP_NU_P_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
             k = cpar%cs_poltype(1,i)
             if (.not. cpar%cs_polarization(i)) k = 1 
             do j = 1,k
                call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_LMAX'//itext, &
                     & len_itext=len_itext,        par_int=cpar%cs_lmax_ind_pol(j,1,i))
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0) then
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_LNLTYPE'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_lnLtype(j,1,i))
                   if (trim(cpar%cs_spec_lnLtype(j,1,i)) == 'prior') then
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PRIOR_MEAN'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(1,j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PRIOR_RMS'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(2,j,1,i))
                   else
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_SAMPLE_NPROP'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_nprop(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_SAMPLE_PROPLEN'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_proplen(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_NPROP_INIT'//itext, &
                           & len_itext=len_itext, par_int=cpar%cs_spec_nprop_init(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PROPLEN_INIT'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_spec_proplen_init(j,1,i))
                   end if
                end if
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_NUM_PIXREG'//itext, &
                        & len_itext=len_itext, par_int=cpar%cs_spec_npixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_FIX_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_fix_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_'//trim(pol_labels(j))//'_PIXREG_PRIORS'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_priors(j,1,i))
                end if
             end do
             do j = 1,k
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_NU_P_PIXREG_MAP'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_map(1,i))
                   exit
                end if
             end do
             bool_flag=.false.
             do j = 1,k
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0 ) bool_flag=.true.
             end do
             if (bool_flag .or. cpar%almsamp_pixreg) &
                  & call get_parameter_hashtable(htbl, 'COMP_NU_P_PIXREG_INITVALUE_MAP'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_pixreg_init_theta(1,i))
             if (any(cpar%cs_lmax_ind_pol(:k,1,i) >= 0)) &
                  & call get_parameter_hashtable(htbl, 'COMP_NU_P_ALMSAMP_INIT'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_almsamp_init(1,i))
             if (bool_flag) then
                call get_parameter_hashtable(htbl, 'COMP_NU_P_UNI_NPROP_LOW'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(1,1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_UNI_NPROP_HIGH'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(2,1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_MASK'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_mask(1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_NPROP'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_nprop(1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_PROPLEN'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_proplen(1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_CORRELATION_CONVERGENCE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_corr_convergence(1,i))
                if (cpar%cs_spec_corr_convergence(1,i))  call get_parameter_hashtable(htbl, &
                     & 'COMP_NU_P_CORRELATION_CONVERGENCE_LIMIT'//itext, &
                     & len_itext=len_itext, par_dp=cpar%cs_spec_corr_limit(1,i))
                call get_parameter_hashtable(htbl, 'COMP_NU_P_COMBINED_MONOPOLE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_mono_combined(i,2))
                if (cpar%cs_spec_mono_combined(i,2)) call get_parameter_hashtable(htbl, &
                     & 'COMP_NU_P_COMBINED_MONOPOLE_MASK'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_mask(i,2))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_NU_P_COMBINED_MONOPOLE_TYPE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_type(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_NU_P_COMBINED_MONOPOLE_FREEZE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_freeze(i,1))
             end if
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
             call get_parameter_hashtable(htbl, 'COMP_W_AME_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(2,i))
             k = cpar%cs_poltype(2,i)
             if (.not. cpar%cs_polarization(i)) k = 1 
             do j = 1,k
                call get_parameter_hashtable(htbl, 'COMP_W_AME_'//trim(pol_labels(j))//'_LMAX'//itext, &
                     & len_itext=len_itext,        par_int=cpar%cs_lmax_ind_pol(j,2,i))
                if (cpar%cs_lmax_ind_pol(j,2,i) < 0) then
                   call get_parameter_hashtable(htbl, 'COMP_W_AME_'//trim(pol_labels(j))//'_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg(j,2,i))
                   call get_parameter_hashtable(htbl, 'COMP_W_AME_'//trim(pol_labels(j))//'_LNLTYPE'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_lnLtype(j,2,i))
                   if (trim(cpar%cs_spec_lnLtype(j,2,i)) == 'prior') then
                      call get_parameter_hashtable(htbl, 'COMP_W_AME_'//trim(pol_labels(j))//'_PRIOR_MEAN'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(1,j,2,i))
                      call get_parameter_hashtable(htbl, 'COMP_W_AME_'//trim(pol_labels(j))//'_PRIOR_RMS'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(2,j,2,i))
                   else
                      call get_parameter_hashtable(htbl, 'COMP_W_AME_'//trim(pol_labels(j))//'_SAMPLE_NPROP'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_nprop(j,2,i))
                      call get_parameter_hashtable(htbl, 'COMP_W_AME_'//trim(pol_labels(j))//'_SAMPLE_PROPLEN'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_proplen(j,2,i))
                      call get_parameter_hashtable(htbl, 'COMP_W_AME_'//trim(pol_labels(j))//'_NPROP_INIT'//itext, &
                           & len_itext=len_itext, par_int=cpar%cs_spec_nprop_init(j,2,i))
                      call get_parameter_hashtable(htbl, 'COMP_W_AME_'//trim(pol_labels(j))//'_PROPLEN_INIT'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_spec_proplen_init(j,2,i))
                   end if
                end if
                if (trim(cpar%cs_spec_pixreg(j,2,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_W_AME_'//trim(pol_labels(j))//'_NUM_PIXREG'//itext, &
                        & len_itext=len_itext, par_int=cpar%cs_spec_npixreg(j,2,i))
                   call get_parameter_hashtable(htbl, 'COMP_W_AME_'//trim(pol_labels(j))//'_FIX_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_fix_pixreg(j,2,i))
                   call get_parameter_hashtable(htbl, 'COMP_W_AME_'//trim(pol_labels(j))//'_PIXREG_PRIORS'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_priors(j,2,i))
                end if
             end do
             do j = 1,k
                if (trim(cpar%cs_spec_pixreg(j,2,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_W_AME_PIXREG_MAP'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_map(2,i))
                   exit
                end if
             end do
             bool_flag=.false.
             do j = 1,k
                if (cpar%cs_lmax_ind_pol(j,2,i) < 0 ) bool_flag=.true.
             end do
             if (bool_flag .or. cpar%almsamp_pixreg) &
                  & call get_parameter_hashtable(htbl, 'COMP_W_AME_PIXREG_INITVALUE_MAP'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_pixreg_init_theta(2,i))
             if (any(cpar%cs_lmax_ind_pol(:k,2,i) >= 0)) &
                  & call get_parameter_hashtable(htbl, 'COMP_W_AME_ALMSAMP_INIT'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_almsamp_init(2,i))
             if (bool_flag) then
                call get_parameter_hashtable(htbl, 'COMP_W_AME_UNI_NPROP_LOW'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(1,2,i))
                call get_parameter_hashtable(htbl, 'COMP_W_AME_UNI_NPROP_HIGH'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(2,2,i))
                call get_parameter_hashtable(htbl, 'COMP_W_AME_MASK'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_mask(2,i))
                call get_parameter_hashtable(htbl, 'COMP_W_AME_NPROP'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_nprop(2,i))
                call get_parameter_hashtable(htbl, 'COMP_W_AME_PROPLEN'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_proplen(2,i))
                call get_parameter_hashtable(htbl, 'COMP_W_AME_CORRELATION_CONVERGENCE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_corr_convergence(2,i))
                if (cpar%cs_spec_corr_convergence(2,i))  call get_parameter_hashtable(htbl, &
                     & 'COMP_W_AME_CORRELATION_CONVERGENCE_LIMIT'//itext, &
                     & len_itext=len_itext, par_dp=cpar%cs_spec_corr_limit(2,i))
                call get_parameter_hashtable(htbl, 'COMP_W_AME_COMBINED_MONOPOLE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_mono_combined(i,2))
                if (cpar%cs_spec_mono_combined(i,2)) call get_parameter_hashtable(htbl, &
                     & 'COMP_W_AME_COMBINED_MONOPOLE_MASK'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_mask(i,2))
                if (cpar%cs_spec_mono_combined(i,2)) call get_parameter_hashtable(htbl, &
                     & 'COMP_W_AME_COMBINED_MONOPOLE_TYPE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_type(i,2))
                if (cpar%cs_spec_mono_combined(i,2)) call get_parameter_hashtable(htbl, &
                     & 'COMP_W_AME_COMBINED_MONOPOLE_FREEZE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_freeze(i,2))
             end if
             call get_parameter_hashtable(htbl, 'COMP_INPUT_W_AME_MAP'//itext, len_itext=len_itext,        &
                  & par_string=cpar%cs_input_ind(2,i))
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_W_AME'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(2,i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_W_AME_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_W_AME_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_W_AME_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,2))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_W_AME_RMS'//itext, len_itext=len_itext,  &
                  & par_dp=cpar%cs_p_gauss(i,2,2))
             call get_parameter_hashtable(htbl, 'COMP_INDMASK'//itext, len_itext=len_itext,         par_string=cpar%cs_indmask(i))
             call get_parameter_hashtable(htbl, 'COMP_NU_P_SMOOTHING_SCALE'//itext, len_itext=len_itext,  &
                  & par_int=cpar%cs_smooth_scale(i,1))
             call get_parameter_hashtable(htbl, 'COMP_W_AME_SMOOTHING_SCALE'//itext, len_itext=len_itext,  &
                  & par_int=cpar%cs_smooth_scale(i,2))
             call get_parameter_hashtable(htbl, 'COMP_NU_P_NU_MIN'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_min(i,1))
             call get_parameter_hashtable(htbl, 'COMP_NU_P_NU_MAX'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_max(i,1))
             call get_parameter_hashtable(htbl, 'COMP_W_AME_NU_MIN'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_min(i,2))          
             call get_parameter_hashtable(htbl, 'COMP_W_AME_NU_MAX'//itext, len_itext=len_itext,   par_dp=cpar%cs_nu_max(i,2))
             call get_parameter_hashtable(htbl, 'COMP_APPLY_JEFFREYS_PRIOR'//itext, len_itext=len_itext,   par_lgt=cpar%cs_apply_jeffreys(i))
             do j=1,2
                if (cpar%cs_smooth_scale(i,1) > cpar%num_smooth_scales) then
                   write(*,fmt='(a,i2,a,i2,a,i2,a,i2)') 'Smoothing scale ',cpar%cs_smooth_scale(i,j), &
                        & ' for index nr. ',j,' in component nr. ', i,' is larger than the number of smoothing scales: ', &
                        & cpar%num_smooth_scales
                   stop
                end if
             end do

          case ('MBB')
             call get_parameter_hashtable(htbl, 'COMP_BETA_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
             k = cpar%cs_poltype(1,i)
             if (.not. cpar%cs_polarization(i)) k = 1 
             do j = 1,k
                call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_LMAX'//itext, &
                     & len_itext=len_itext,        par_int=cpar%cs_lmax_ind_pol(j,1,i))
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0) then
                   call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_LNLTYPE'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_lnLtype(j,1,i))
                   if (trim(cpar%cs_spec_lnLtype(j,1,i)) == 'prior') then
                      call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_PRIOR_MEAN'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(1,j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_PRIOR_RMS'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(2,j,1,i))
                   else
                      call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_SAMPLE_NPROP'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_nprop(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_SAMPLE_PROPLEN'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_proplen(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_NPROP_INIT'//itext, &
                           & len_itext=len_itext, par_int=cpar%cs_spec_nprop_init(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_PROPLEN_INIT'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_spec_proplen_init(j,1,i))
                   end if
                end if
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_NUM_PIXREG'//itext, &
                        & len_itext=len_itext, par_int=cpar%cs_spec_npixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_FIX_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_fix_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_BETA_'//trim(pol_labels(j))//'_PIXREG_PRIORS'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_priors(j,1,i))
                end if
             end do
             do j = 1,k
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_BETA_PIXREG_MAP'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_map(1,i))
                   exit
                end if
             end do
             bool_flag=.false.
             do j = 1,k
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0 ) bool_flag=.true.
             end do
             if (bool_flag .or. cpar%almsamp_pixreg) &
                  & call get_parameter_hashtable(htbl, 'COMP_BETA_PIXREG_INITVALUE_MAP'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_pixreg_init_theta(1,i))
             if (any(cpar%cs_lmax_ind_pol(:k,1,i) >= 0)) &
                  & call get_parameter_hashtable(htbl, 'COMP_BETA_ALMSAMP_INIT'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_almsamp_init(1,i))
             if (bool_flag) then
                call get_parameter_hashtable(htbl, 'COMP_BETA_UNI_NPROP_LOW'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(1,1,i))
                call get_parameter_hashtable(htbl, 'COMP_BETA_UNI_NPROP_HIGH'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(2,1,i))
                call get_parameter_hashtable(htbl, 'COMP_BETA_MASK'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_mask(1,i))
                call get_parameter_hashtable(htbl, 'COMP_BETA_NPROP'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_nprop(1,i))
                call get_parameter_hashtable(htbl, 'COMP_BETA_PROPLEN'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_proplen(1,i))
                call get_parameter_hashtable(htbl, 'COMP_BETA_CORRELATION_CONVERGENCE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_corr_convergence(1,i))
                if (cpar%cs_spec_corr_convergence(1,i))  call get_parameter_hashtable(htbl, &
                     & 'COMP_BETA_CORRELATION_CONVERGENCE_LIMIT'//itext, &
                     & len_itext=len_itext, par_dp=cpar%cs_spec_corr_limit(1,i))
                call get_parameter_hashtable(htbl, 'COMP_BETA_COMBINED_MONOPOLE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_mono_combined(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_BETA_COMBINED_MONOPOLE_MASK'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_mask(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_BETA_COMBINED_MONOPOLE_TYPE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_type(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_BETA_COMBINED_MONOPOLE_FREEZE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_freeze(i,1))
             end if
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
             k = cpar%cs_poltype(2,i)
             if (.not. cpar%cs_polarization(i)) k = 1 
             do j = 1,k
                call get_parameter_hashtable(htbl, 'COMP_T_'//trim(pol_labels(j))//'_LMAX'//itext, &
                     & len_itext=len_itext,        par_int=cpar%cs_lmax_ind_pol(j,2,i))
                if (cpar%cs_lmax_ind_pol(j,2,i) < 0) then
                   call get_parameter_hashtable(htbl, 'COMP_T_'//trim(pol_labels(j))//'_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg(j,2,i))
                   call get_parameter_hashtable(htbl, 'COMP_T_'//trim(pol_labels(j))//'_LNLTYPE'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_lnLtype(j,2,i))
                   if (trim(cpar%cs_spec_lnLtype(j,2,i)) == 'prior') then
                      call get_parameter_hashtable(htbl, 'COMP_T_'//trim(pol_labels(j))//'_PRIOR_MEAN'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(1,j,2,i))
                      call get_parameter_hashtable(htbl, 'COMP_T_'//trim(pol_labels(j))//'_PRIOR_RMS'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(2,j,2,i))
                   else
                      call get_parameter_hashtable(htbl, 'COMP_T_'//trim(pol_labels(j))//'_SAMPLE_NPROP'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_nprop(j,2,i))
                      call get_parameter_hashtable(htbl, 'COMP_T_'//trim(pol_labels(j))//'_SAMPLE_PROPLEN'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_proplen(j,2,i))
                      call get_parameter_hashtable(htbl, 'COMP_T_'//trim(pol_labels(j))//'_NPROP_INIT'//itext, &
                           & len_itext=len_itext, par_int=cpar%cs_spec_nprop_init(j,2,i))
                      call get_parameter_hashtable(htbl, 'COMP_T_'//trim(pol_labels(j))//'_PROPLEN_INIT'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_spec_proplen_init(j,2,i))
                   end if
                end if
                if (trim(cpar%cs_spec_pixreg(j,2,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_T_'//trim(pol_labels(j))//'_NUM_PIXREG'//itext, &
                        & len_itext=len_itext, par_int=cpar%cs_spec_npixreg(j,2,i))
                   call get_parameter_hashtable(htbl, 'COMP_T_'//trim(pol_labels(j))//'_FIX_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_fix_pixreg(j,2,i))
                   call get_parameter_hashtable(htbl, 'COMP_T_'//trim(pol_labels(j))//'_PIXREG_PRIORS'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_priors(j,2,i))
                end if
             end do
             do j = 1,k
                if (trim(cpar%cs_spec_pixreg(j,2,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_T_PIXREG_MAP'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_map(2,i))
                   exit
                end if
             end do
             bool_flag=.false.
             do j = 1,k
                if (cpar%cs_lmax_ind_pol(j,2,i) < 0 ) bool_flag=.true.
             end do
             if (bool_flag .or. cpar%almsamp_pixreg) &
                  & call get_parameter_hashtable(htbl, 'COMP_T_PIXREG_INITVALUE_MAP'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_pixreg_init_theta(2,i))
             if (any(cpar%cs_lmax_ind_pol(:k,2,i) >= 0)) &
                  & call get_parameter_hashtable(htbl, 'COMP_T_ALMSAMP_INIT'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_almsamp_init(2,i))
             if (bool_flag) then
                call get_parameter_hashtable(htbl, 'COMP_T_UNI_NPROP_LOW'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(1,2,i))
                call get_parameter_hashtable(htbl, 'COMP_T_UNI_NPROP_HIGH'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(2,2,i))
                call get_parameter_hashtable(htbl, 'COMP_T_MASK'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_mask(2,i))
                call get_parameter_hashtable(htbl, 'COMP_T_NPROP'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_nprop(2,i))
                call get_parameter_hashtable(htbl, 'COMP_T_PROPLEN'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_proplen(2,i))
                call get_parameter_hashtable(htbl, 'COMP_T_CORRELATION_CONVERGENCE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_corr_convergence(2,i))
                if (cpar%cs_spec_corr_convergence(2,i))  call get_parameter_hashtable(htbl, &
                     & 'COMP_T_CORRELATION_CONVERGENCE_LIMIT'//itext, &
                     & len_itext=len_itext, par_dp=cpar%cs_spec_corr_limit(2,i))
                call get_parameter_hashtable(htbl, 'COMP_T_COMBINED_MONOPOLE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_mono_combined(i,2))
                if (cpar%cs_spec_mono_combined(i,2)) call get_parameter_hashtable(htbl, &
                     & 'COMP_T_COMBINED_MONOPOLE_MASK'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_mask(i,2))
                if (cpar%cs_spec_mono_combined(i,2)) call get_parameter_hashtable(htbl, &
                     & 'COMP_T_COMBINED_MONOPOLE_TYPE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_type(i,2))
                if (cpar%cs_spec_mono_combined(i,2)) call get_parameter_hashtable(htbl, &
                     & 'COMP_T_COMBINED_MONOPOLE_FREEZE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_freeze(i,2))
             end if
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
             do j=1,2
                if (cpar%cs_smooth_scale(i,1) > cpar%num_smooth_scales) then
                   write(*,fmt='(a,i2,a,i2,a,i2,a,i2)') 'Smoothing scale ',cpar%cs_smooth_scale(i,j), &
                        & ' for index nr. ',j,' in component nr. ', i,' is larger than the number of smoothing scales: ', &
                        & cpar%num_smooth_scales
                   stop
                end if
             end do

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
             call get_parameter_hashtable(htbl, 'COMP_T_E_POLTYPE'//itext, len_itext=len_itext,  par_int=cpar%cs_poltype(1,i))
             k = cpar%cs_poltype(1,i)
             if (.not. cpar%cs_polarization(i)) k = 1 
             do j = 1,k
                call get_parameter_hashtable(htbl, 'COMP_T_E_'//trim(pol_labels(j))//'_LMAX'//itext, &
                     & len_itext=len_itext,        par_int=cpar%cs_lmax_ind_pol(j,1,i))
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0) then
                   call get_parameter_hashtable(htbl, 'COMP_T_E_'//trim(pol_labels(j))//'_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_T_E_'//trim(pol_labels(j))//'_LNLTYPE'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_lnLtype(j,1,i))
                   if (trim(cpar%cs_spec_lnLtype(j,1,i)) == 'prior') then
                      call get_parameter_hashtable(htbl, 'COMP_T_E_'//trim(pol_labels(j))//'_PRIOR_MEAN'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(1,j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_T_E_'//trim(pol_labels(j))//'_PRIOR_RMS'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_theta_prior(2,j,1,i))
                   else
                      call get_parameter_hashtable(htbl, 'COMP_T_E_'//trim(pol_labels(j))//'_SAMPLE_NPROP'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_nprop(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_T_E_'//trim(pol_labels(j))//'_SAMPLE_PROPLEN'//itext, &
                           & len_itext=len_itext, par_lgt=cpar%cs_spec_samp_proplen(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_T_E_'//trim(pol_labels(j))//'_NPROP_INIT'//itext, &
                           & len_itext=len_itext, par_int=cpar%cs_spec_nprop_init(j,1,i))
                      call get_parameter_hashtable(htbl, 'COMP_T_E_'//trim(pol_labels(j))//'_PROPLEN_INIT'//itext, &
                           & len_itext=len_itext, par_dp=cpar%cs_spec_proplen_init(j,1,i))
                   end if
                end if
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_T_E_'//trim(pol_labels(j))//'_NUM_PIXREG'//itext, &
                        & len_itext=len_itext, par_int=cpar%cs_spec_npixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_T_E_'//trim(pol_labels(j))//'_FIX_PIXREG'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_fix_pixreg(j,1,i))
                   call get_parameter_hashtable(htbl, 'COMP_T_E_'//trim(pol_labels(j))//'_PIXREG_PRIORS'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_priors(j,1,i))
                end if
             end do
             do j = 1,k
                if (trim(cpar%cs_spec_pixreg(j,1,i)) == 'pixreg' .or. cpar%almsamp_pixreg) then
                   call get_parameter_hashtable(htbl, 'COMP_T_E_PIXREG_MAP'//itext, &
                        & len_itext=len_itext, par_string=cpar%cs_spec_pixreg_map(1,i))
                   exit
                end if
             end do
             bool_flag=.false.
             do j = 1,k
                if (cpar%cs_lmax_ind_pol(j,1,i) < 0 ) bool_flag=.true.
             end do
             if (bool_flag .or. cpar%almsamp_pixreg) &
                  & call get_parameter_hashtable(htbl, 'COMP_T_E_PIXREG_INITVALUE_MAP'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_pixreg_init_theta(1,i))
             if (any(cpar%cs_lmax_ind_pol(:k,1,i) >= 0)) &
                  & call get_parameter_hashtable(htbl, 'COMP_T_E_ALMSAMP_INIT'//itext, &
                  & len_itext=len_itext, par_string=cpar%cs_almsamp_init(1,i))
             if (bool_flag) then
                call get_parameter_hashtable(htbl, 'COMP_T_E_UNI_NPROP_LOW'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(1,1,i))
                call get_parameter_hashtable(htbl, 'COMP_T_E_UNI_NPROP_HIGH'//itext, len_itext=len_itext,  &
                     & par_int=cpar%cs_spec_uni_nprop(2,1,i))
                call get_parameter_hashtable(htbl, 'COMP_T_E_MASK'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_mask(1,i))
                call get_parameter_hashtable(htbl, 'COMP_T_E_NPROP'//itext, & 
                     & len_itext=len_itext, par_string=cpar%cs_spec_nprop(1,i))
                call get_parameter_hashtable(htbl, 'COMP_T_E_PROPLEN'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_proplen(1,i))
                call get_parameter_hashtable(htbl, 'COMP_T_E_CORRELATION_CONVERGENCE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_corr_convergence(1,i))
                if (cpar%cs_spec_corr_convergence(1,i))  call get_parameter_hashtable(htbl, &
                     & 'COMP_T_E_CORRELATION_CONVERGENCE_LIMIT'//itext, &
                     & len_itext=len_itext, par_dp=cpar%cs_spec_corr_limit(1,i))
                call get_parameter_hashtable(htbl, 'COMP_T_E_COMBINED_MONOPOLE_SAMPLING'//itext, &
                     & len_itext=len_itext, par_lgt=cpar%cs_spec_mono_combined(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_T_E_COMBINED_MONOPOLE_MASK'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_mask(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_T_E_COMBINED_MONOPOLE_TYPE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_type(i,1))
                if (cpar%cs_spec_mono_combined(i,1)) call get_parameter_hashtable(htbl, &
                     & 'COMP_T_E_COMBINED_MONOPOLE_FREEZE'//itext, &
                     & len_itext=len_itext, par_string=cpar%cs_spec_mono_freeze(i,1))
             end if
             call get_parameter_hashtable(htbl, 'COMP_INPUT_T_E_MAP'//itext, len_itext=len_itext,        &
                  & par_string=cpar%cs_input_ind(1,i))
             call get_parameter_hashtable(htbl, 'COMP_DEFAULT_T_E'//itext, len_itext=len_itext,          &
                  & par_dp=cpar%cs_theta_def(1,i))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_T_E_LOW'//itext, len_itext=len_itext,    &
                  & par_dp=cpar%cs_p_uni(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_T_E_HIGH'//itext, len_itext=len_itext,   &
                  & par_dp=cpar%cs_p_uni(i,2,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_T_E_MEAN'//itext, len_itext=len_itext, &
                  & par_dp=cpar%cs_p_gauss(i,1,1))
             call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_T_E_RMS'//itext, len_itext=len_itext,  &
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
             do j=1,1
                if (cpar%cs_smooth_scale(i,1) > cpar%num_smooth_scales) then
                   write(*,fmt='(a,i2,a,i2,a,i2,a,i2)') 'Smoothing scale ',cpar%cs_smooth_scale(i,j), &
                        & ' for index nr. ',j,' in component nr. ', i,' is larger than the number of smoothing scales: ', &
                        & cpar%num_smooth_scales
                   stop
                end if
             end do

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
          call get_parameter_hashtable(htbl, 'COMP_INIT_CATALOG'//itext, len_itext=len_itext,  par_string=cpar%cs_init_catalog(i))
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
          call get_parameter_hashtable(htbl, 'COMP_CG_SCALE'//itext, len_itext=len_itext, par_dp=cpar%cs_cg_scale(1,i))
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


    subroutine read_zodi_params_hash(htbl, cpar)
        implicit none

        type(hash_tbl_sll), intent(in) :: htbl
        type(comm_params),  intent(inout) :: cpar

        ! Read in model hyper parameters
        call get_parameter_hashtable(htbl, 'ZODI_GAUSS_QUAD_ORDER', par_int=cpar%zs_gauss_quad_order)
        call get_parameter_hashtable(htbl, 'ZODI_LOS_CUT', par_dp=cpar%zs_los_cut)
        call get_parameter_hashtable(htbl, 'ZODI_DELTA_T', par_dp=cpar%zs_delta_t)
        call get_parameter_hashtable(htbl, 'ZODI_USE_CLOUD', par_lgt=cpar%zs_use_cloud
        call get_parameter_hashtable(htbl, 'ZODI_USE_BAND1', par_lgt=cpar%zs_use_band1)
        call get_parameter_hashtable(htbl, 'ZODI_USE_BAND2', par_lgt=cpar%zs_use_band2)
        call get_parameter_hashtable(htbl, 'ZODI_USE_BAND3', par_lgt=cpar%zs_use_band3)
        call get_parameter_hashtable(htbl, 'ZODI_USE_RING', par_lgt=cpar%zs_use_ring)
        call get_parameter_hashtable(htbl, 'ZODI_USE_FEATURE', par_lgt=cpar%zs_use_feature)
        call get_parameter_hashtable(htbl, 'ZODI_USE_UNIT_EMISSIVITY', par_lgt=cpar%zs_use_unit_emissivity)
        call get_parameter_hashtable(htbl, 'ZODI_FREQ_CORRECTION_TYPE', par_string=cpar%zs_freq_correction_type)

        print *, cpar%zs_freq_correction_type
        stop

        call get_parameter_hashtable(htbl, 'ZODI_CLOUD_MU', par_dp=cpar%zs_cloud_mu)
        stop

    end subroutine read_zodi_params_hash


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
    do i = 1, command_argument_count() !iargc()
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
    integer(i4b)               :: depth, units(maxdepth), i, num, ounit, stat
    character(len=1024)        :: key, value, arg, default_path

    call get_environment_variable("COMMANDER_PARAMS_DEFAULT", default_path, status=stat)


    num = 0
    depth = 1
    ounit = getlun()
    open(ounit,file=outfile,action="write")
    write(ounit,fmt="(a)",advance="no") '# Arguments:'
    do i = 1, command_argument_count() !iargc()
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
       else if (key=='@DEFAULT') then
          ! Recurse into the new file
          if(stat /=0) then
            write(*,*) "@DEFAULT directive present but getting COMMANDER_PARAMS_DEFAULT returned ", stat
            stop
          end if
          read(units(depth),*,end=1) key, value
          write(ounit,fmt='(a)') pad("",depth-1," ") // "# File: " // trim(default_path)//'/'//trim(value)
          depth=depth+1
          units(depth) = getlun()
          open(units(depth),file=trim(default_path)//'/'//trim(value),status="old")
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

  subroutine str2int(str,int,stat)
    implicit none
    ! Arguments
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int
    integer,intent(out)         :: stat

    read(str,*,iostat=stat)  int
  end subroutine str2int

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

  subroutine get_detectors(filename, directory, detectors, num_dets)
     !
     ! Reads detector names from a text file and saves them in a character array.
     !
     ! Arguments:
     ! ----------
     ! filename:  character string
     !            Filename of the file where detector names are stored.
     ! directory: character string
     !            Directory where file is stored.
     ! num_dets:  integer (optional)
     !            Number of detectors
     !
     ! Return:
     ! -------
     ! detectors: character array
     !            Initially empty array is filled with detector names. 
     ! 
     implicit none
     character(len=*), intent(in)           :: filename, directory
     character(len=*), intent(inout)        :: detectors(:)
     integer(i4b),     intent(in), optional :: num_dets

     character(len=500)           :: detector_list_file
     integer(i4b)                 :: unit,io_error,counter, ndet, i
     character(len=8)             :: line

     if (present(num_dets)) then
          ndet = num_dets
     else
          ndet = size(detectors)
     end if

     unit = 20
     detector_list_file = trim(adjustl(directory))//'/'//trim(adjustl(filename))

     open(unit,file=trim(detector_list_file),status='old',action='read',iostat=io_error)
     if (io_error == 0) then
     ! Do nothing
     else
         write(*,*) 'Could not open file: ', trim(adjustl(detector_list_file))
         stop
     end if

     do i=1, ndet
         read(unit,'(a)') line
         if ((line(1:1) == '#') .or. (line(1:1) == '')) then
            cycle
         else
            detectors(i) = line
         end if
     end do
     close(unit)
  end subroutine get_detectors

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

  integer(i4b) function count_detectors(filename, directory)
     ! 
     ! Takes in the filename and directory of a detector list and returns the number of 
     ! detectors in that list. Each detector has to be written on a separate line, as 
     ! the function simply counts the lines of the file that don't start in '#'.
     !
     ! Arguments:
     ! ----------
     ! filename:    character string
     !              Filename of the detector list             
     ! directory:   character string
     !              Directory where file is located
     !
     ! Returns:
     ! --------
     ! count_detectors: integer
     !                  Number of lines in the file that are not commented out using '#'.
     !
     implicit none
     character(len=*) :: filename, directory

     character(len=500)           :: detector_list_file
     integer(i4b)                 :: unit,io_error,counter
     logical                      :: counting
     character(len=8)             :: line

     unit = 20
     detector_list_file = trim(adjustl(directory))//'/'//trim(adjustl(filename))

     open(unit,file=detector_list_file, status='old', action='read', iostat=io_error)
     if (io_error == 0) then
          ! Do nothing
     else
          write(*,*) 'Could not open file: ', trim(adjustl(detector_list_file))
          stop
     end if

     counting = .true.
     counter = 0
     do while(counting)
          read(unit,'(a)',end=1) line
          if ((line(1:1) == '#') .or. (line(1:1) == '')) then
          cycle
          else
          counter = counter + 1
          end if
     end do
     1  close(unit)

     count_detectors = counter
  end function count_detectors

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
    logical(lgt) :: exist

    datadir  = trim(cpar%datadir) // '/'
    chaindir = trim(cpar%outdir) // '/'

#ifdef USE_INTEL   
      !verify that the output directory exists
      inquire(directory=cpar%outdir, exist=exist) 
      if (.not. exist) then
        write(*,*) "Error: the specified output directory ", trim(cpar%outdir), " does not exist"
        stop
      end if 
#endif

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
          if (trim(cpar%ds_noise_rms_smooth(i,j)) /= 'none' .and. trim(cpar%ds_noise_rms_smooth(i,j))/= 'native' .and. trim(cpar%ds_noise_rms_smooth(i,j)) /= "") &
               & call validate_file(trim(datadir)//trim(cpar%ds_noise_rms_smooth(i,j)))  ! Smoothed RMS file
       end do

       if (cpar%enable_TOD_analysis .and. trim(cpar%ds_tod_type(i)) /= 'none') then
          call validate_file(trim(datadir)//trim(cpar%ds_tod_procmask1(i)))  ! Procmask1
          call validate_file(trim(datadir)//trim(cpar%ds_tod_procmask2(i)))  ! Procmask2
          call validate_file(trim(datadir)//trim(cpar%ds_tod_filelist(i)))   ! Filelist
          if (trim(cpar%ds_tod_jumplist(i)) /= 'none') &
               & call validate_file(trim(datadir)//trim(cpar%ds_tod_jumplist(i)))   ! Jumplist
          call validate_file(trim(datadir)//trim(cpar%ds_tod_instfile(i)))   ! Instrument file, RIMO
          if (trim(cpar%ds_tod_bp_init(i)) /= 'none') call validate_file(trim(datadir)//trim(cpar%ds_tod_bp_init(i)))    ! BP prop and init
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
             if (cpar%cs_spec_mono_combined(i,1) .and. trim(cpar%cs_spec_mono_mask(i,1)) /= 'fullsky') &
                  & call validate_file(trim(datadir)//trim(cpar%cs_spec_mono_mask(i,1)))
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
             if (cpar%cs_spec_mono_combined(i,1) .and. trim(cpar%cs_spec_mono_mask(i,1)) /= 'fullsky') &
                  & call validate_file(trim(datadir)//trim(cpar%cs_spec_mono_mask(i,1)))
          case ('spindust2')
             if (trim(cpar%cs_input_ind(1,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(1,i)))
             if (trim(cpar%cs_input_ind(2,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(2,i)))
             call validate_file(trim(datadir)//trim(cpar%cs_SED_template(1,i)))
             if (cpar%cs_spec_mono_combined(i,1) .and. trim(cpar%cs_spec_mono_mask(i,1)) /= 'fullsky') &
                  & call validate_file(trim(datadir)//trim(cpar%cs_spec_mono_mask(i,1)))
             if (cpar%cs_spec_mono_combined(i,2) .and. trim(cpar%cs_spec_mono_mask(i,2)) /= 'fullsky') &
                  & call validate_file(trim(datadir)//trim(cpar%cs_spec_mono_mask(i,2)))
          case ('MBB')
             if (trim(cpar%cs_input_ind(1,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(1,i)))
             if (trim(cpar%cs_input_ind(2,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(2,i)))
             if (cpar%cs_spec_mono_combined(i,1) .and. trim(cpar%cs_spec_mono_mask(i,1)) /= 'fullsky') &
                  & call validate_file(trim(datadir)//trim(cpar%cs_spec_mono_mask(i,1)))
             if (cpar%cs_spec_mono_combined(i,2) .and. trim(cpar%cs_spec_mono_mask(i,2)) /= 'fullsky') &
                  & call validate_file(trim(datadir)//trim(cpar%cs_spec_mono_mask(i,2)))
          case ('freefree')
!!$             if (trim(cpar%cs_input_ind(1,i)) /= 'default') &
!!$                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(1,i)))
             if (trim(cpar%cs_input_ind(1,i)) /= 'default') &
                  call validate_file(trim(datadir)//trim(cpar%cs_input_ind(1,i)))             
             if (cpar%cs_spec_mono_combined(i,1) .and. trim(cpar%cs_spec_mono_mask(i,1)) /= 'fullsky') &
                  & call validate_file(trim(datadir)//trim(cpar%cs_spec_mono_mask(i,1)))
          case ('line')
             call validate_file(trim(datadir)//trim(cpar%cs_SED_template(1,i)))
          end select

       else if (trim(cpar%cs_class(i)) == 'ptsrc') then
          call validate_file(trim(datadir)//trim(cpar%cs_catalog(i)))
          if (trim(cpar%cs_init_catalog(i)) /= 'none') then
             call validate_file(trim(datadir)//trim(cpar%cs_init_catalog(i)))
          end if
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

  subroutine read_paramfile_to_ascii(paramfile,paramfile_cache, paramfile_len)
    implicit none
    character(len=512),                            intent(in)  :: paramfile
    character(len=512), allocatable, dimension(:), intent(inout) :: paramfile_cache
    integer(i4b),intent(inout) :: paramfile_len

    integer(i4b), parameter    :: maxdepth = 256
    integer(i4b)               :: depth, units(maxdepth), line_nr,i, stat, pos
    character(len=512)         :: key, value, filenames(maxdepth), line
    character(len=1024)        :: default_path
    character(len=3)           :: band_num

    character(len=512), allocatable, dimension(:) :: new_cache

    ! read file to ascii array

    
    call get_environment_variable("COMMANDER_PARAMS_DEFAULT", default_path, status=stat)

    band_num = 'XXX'
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
             read(units(depth),fmt="(a)",end=1) line
             pos = index(line, ' ')
             value = adjustl(line(pos:len(line)))
             pos = index(value, ' ')
             value = trim(value(1:pos))

             depth=depth+1
             units(depth) = getlun()
             filenames(depth) = value
             open(units(depth),file=value,status="old",err=2)
          else if(key == '@DEFAULT') then
             if(stat /= 0) then
               write(*,*) "Paramater file uses @DEFAULT command but the environment variable COMMANDER_PARAMS_DEFAULT returns ", stat
               stop
             end if
             ! Recurse to the default new file
             read(units(depth),fmt="(a)",end=1) line
             pos = index(line, ' ')
             value = adjustl(line(pos:len(line)))
             pos = index(value, ' ')
             value = trim(value(1:pos))

             depth = depth+1
             units(depth) = getlun()
             filenames(depth) = trim(default_path)//'/'//trim(value)
             open(units(depth),file=filenames(depth), status="old", err=2)
          else if(key == '@START') then
             read(units(depth),*,end=1) key, value
             if(band_num /= 'XXX') then
               write(*,*) "Error starting band number ", trim(value), ", band ", band_num, " has not ended in file ", trim(filenames(depth))
               stop
             end if
             band_num = value(1:len(value))
          else if(key == '@END') then
             read(units(depth),*,end=1) key, value
             if(value(1:len(value)) /= band_num) then
               write(*,*) "Error ending band ", trim(value), ", current band is ", band_num, " in file ", trim(filenames(depth))
               stop
             end if
             band_num = 'XXX'

          else
             goto 3
          end if
       else
          read(units(depth),fmt="(a)") line
          !if we get here we have read a new line from the parameter file(s)
          line_nr = line_nr + 1
          if(line_nr > paramfile_len) then !we need to resize the cache array
            allocate(new_cache(2*paramfile_len))

            new_cache(1:paramfile_len) = paramfile_cache(1:paramfile_len)
            deallocate(paramfile_cache)
            call move_alloc(new_cache, paramfile_cache)

            paramfile_len = paramfile_len*2
          end if
          if(band_num /= 'XXX') then !active @START directive
            !replace the string &&& with the band number given by @START
            pos = index(line, '&&&') !if this is a band
            if(pos > 0)  line(pos:pos+2)=band_num
            pos = index(line, '&&') !this could be a component
            if(pos > 0) line(pos:pos+1)=band_num
          else
            pos = index(line, '&&') !check for the special chars outside START
            if(pos > 0) then
              write(*,*) "Warning: parameter line ", line, " found outside of a START-END block"
            end if
          end if
          write(paramfile_cache(line_nr),fmt="(a)") line

       end if
       cycle
       ! We get here if we reached the end of a file. Close it and
       ! return to the file above.
1      close(units(depth))
       !write(*,*) "Exiting file " // filenames(depth)
       depth = depth-1
    end do
    !resize the cache to be the exact length
    allocate(new_cache(line_nr))

    new_cache(1:line_nr) = paramfile_cache(1:line_nr)
    deallocate(paramfile_cache)

    call move_alloc(new_cache, paramfile_cache)
    paramfile_len = line_nr
    return

    ! ===== Error handling section ======

    ! Case 1: Include file error
2   write(*,*) "Error: Cannot open include file '" // trim(filenames(depth)) // "'"
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

  ! outputs the parameter file to the path provided
  subroutine save_ascii_parameter_file(outfile, ascii_table)
    implicit none
    character(len=512), intent(in) :: outfile
    character(len=512), dimension(:), intent(in) :: ascii_table

    integer(i4b)      :: unit, i

    write(*,*) "Saving parameter file to ", trim(outfile)

    unit = getlun()
    open(unit, file=trim(outfile),err=1)

    do i=1, size(ascii_table) 
      write(unit, '(a)') trim(ascii_table(i))
    end do
    close(unit)
    return
    !If the opening of the output parameter file fails
1   write(*,*) "Error: Cannot open output file '" // trim(outfile) // "'"
    stop
  end subroutine save_ascii_parameter_file

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
       goto 1
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
       !read(val,*) par_string
       par_string = val
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
          cpar%cg_num_samp_groups                             = cpar%cg_num_samp_groups + 1
          cpar%cg_samp_group(cpar%cg_num_samp_groups)         = trim(cpar%cs_label(i))
          cpar%cg_samp_group_mask(cpar%cg_num_samp_groups)    = 'fullsky'
          if (trim(cpar%cs_class(i)) == 'diffuse') then
             if (trim(cpar%cs_type(i)) == 'cmb') then
                cpar%cg_samp_group_maxiter(cpar%cg_num_samp_groups) = 150
             else
                cpar%cg_samp_group_maxiter(cpar%cg_num_samp_groups) = cpar%cs_cg_samp_group_maxiter(i)
             end if
          else
             cpar%cg_samp_group_maxiter(cpar%cg_num_samp_groups) = 150
          end if
       end if
    end do

    ! Expand md type if present
    cpar%cg_samp_group_md = -1 !no pure mono-/dipole CG sample group exists 
    do i = 1, cpar%cg_num_samp_groups
       call get_tokens(cpar%cg_samp_group(i), ",", comp_label, n)
       do j = 1, n
          if (trim(comp_label(j)) == 'md') then
             if (n==1 .and. cpar%cg_samp_group_md < 0) then
                cpar%cg_samp_group_md = i !a pure mono-/dipole CG sample group exists, used in specific cases 
             else
             end if
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
