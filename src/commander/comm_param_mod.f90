module comm_param_mod
  use comm_utils
  use comm_status_mod
  implicit none

  ! Note: This module reads in the Commander parameter file as the first operation
  !       in the program. This is primarily intended to avoid crashes after hours
  !       of running because of user errors; catch these early, and report back
  !       problems. Then, copy parameters over to module structures if convenient
  !       at a later stage. 

  integer(i4b), parameter, private :: MAXPAR = 10
  type(status_file)                :: status
  
  type comm_params

     ! MPI info
     integer(i4b) :: myid, numprocs, root = 0
     integer(i4b) :: myid_chain, numprocs_chain, comm_chain, mychain
     integer(i4b), dimension(MPI_STATUS_SIZE)          :: status

     ! Global parameters
     character(len=24)  :: operation
     integer(i4b)       :: verbosity, base_seed, numchain
     integer(i4b)       :: num_gibbs_iter, num_ml_iter
     logical(lgt)       :: override_init
     character(len=512) :: initfile
     real(dp)           :: T_CMB
     character(len=512) :: MJysr_convention

     ! Numerical parameters
     integer(i4b)      :: cg_lmax_precond, cg_maxiter
     real(dp)          :: cg_tol

     ! Data parameters
     integer(i4b)       :: numband
     character(len=512) :: datadir
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
     character(len=512), allocatable, dimension(:)   :: ds_maskfile
     character(len=512), allocatable, dimension(:)   :: ds_maskfile_calib
     character(len=512), allocatable, dimension(:)   :: ds_beamtype
     character(len=512), allocatable, dimension(:)   :: ds_blfile
     character(len=512), allocatable, dimension(:)   :: ds_pixwin
     logical(lgt),       allocatable, dimension(:)   :: ds_samp_monopole
     logical(lgt),       allocatable, dimension(:)   :: ds_samp_dipole
     logical(lgt),       allocatable, dimension(:)   :: ds_samp_noiseamp
     character(len=512), allocatable, dimension(:)   :: ds_bptype
     character(len=512), allocatable, dimension(:)   :: ds_bpfile
     character(len=512), allocatable, dimension(:)   :: ds_bpmodel
     real(dp),           allocatable, dimension(:)   :: ds_nu_c
     real(dp),           allocatable, dimension(:)   :: ds_co2t
     logical(lgt),       allocatable, dimension(:)   :: ds_sample_gain
     character(len=512), allocatable, dimension(:)   :: ds_gain_calib_comp
     integer(i4b),       allocatable, dimension(:)   :: ds_gain_lmin
     integer(i4b),       allocatable, dimension(:)   :: ds_gain_lmax
     character(len=512), allocatable, dimension(:)   :: ds_gain_apodmask
     character(len=512), allocatable, dimension(:)   :: ds_gain_fwhm
     integer(i4b),       allocatable, dimension(:)   :: ds_numtemp
     character(len=512), allocatable, dimension(:,:) :: ds_tempname
     logical(lgt),       allocatable, dimension(:,:) :: ds_samptemp
     real(dp),           allocatable, dimension(:,:) :: ds_defaults
     real(dp),           allocatable, dimension(:,:) :: ds_temp_defaults

     ! Component parameters
     integer(i4b) :: cs_ncomp
     logical(lgt),       allocatable, dimension(:)     :: cs_include
     character(len=512), allocatable, dimension(:)     :: cs_label
     character(len=512), allocatable, dimension(:)     :: cs_type
     character(len=512), allocatable, dimension(:)     :: cs_class
     logical(lgt),       allocatable, dimension(:)     :: cs_polarization
     integer(i4b),       allocatable, dimension(:)     :: cs_nside
     integer(i4b),       allocatable, dimension(:)     :: cs_lmax_amp
     integer(i4b),       allocatable, dimension(:)     :: cs_lmax_ind
     character(len=512), allocatable, dimension(:)     :: cs_unit
     real(dp),           allocatable, dimension(:)     :: cs_nu_ref
     character(len=512), allocatable, dimension(:)     :: cs_cltype
     logical(lgt),       allocatable, dimension(:)     :: cs_samp_cls
     character(len=512), allocatable, dimension(:)     :: cs_clfile
     character(len=512), allocatable, dimension(:)     :: cs_binfile
     integer(i4b),       allocatable, dimension(:)     :: cs_lpivot
     character(len=512), allocatable, dimension(:)     :: cs_mask
     real(dp),           allocatable, dimension(:)     :: cs_amp_def
     character(len=512), allocatable, dimension(:)     :: cs_filedef
     real(dp),           allocatable, dimension(:,:)   :: cs_theta_def
     real(dp),           allocatable, dimension(:,:,:) :: cs_p_gauss
     real(dp),           allocatable, dimension(:,:,:) :: cs_p_uni
     
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

  subroutine initialize_mpi_struct(cpar)
    implicit none

    type(comm_params), intent(inout) :: cpar

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
    
    deallocate(ind)

  end subroutine initialize_mpi_struct
  
  ! ********************************************************
  !              Specialized routines; one per module
  ! ********************************************************

  subroutine read_global_params(paramfile, cpar)
    implicit none

    character(len=*),  intent(in)    :: paramfile
    type(comm_params), intent(inout) :: cpar

    call get_parameter(paramfile, 'VERBOSITY',                par_int=cpar%verbosity)
    call get_parameter(paramfile, 'OPERATION',                par_string=cpar%operation)

    call get_parameter(paramfile, 'BASE_SEED',                par_int=cpar%base_seed)
    call get_parameter(paramfile, 'NUMCHAIN',                 par_int=cpar%numchain)
    call get_parameter(paramfile, 'NUM_GIBBS_ITER',           par_int=cpar%num_gibbs_iter)
    call get_parameter(paramfile, 'NUM_ITER_WITH_ML_SEARCH',  par_int=cpar%num_ml_iter)
    call get_parameter(paramfile, 'OVERRIDE_CHAIN_INIT',      par_lgt=cpar%override_init)
    call get_parameter(paramfile, 'INIT_CHAIN_FILE',          par_string=cpar%initfile)

    call get_parameter(paramfile, 'CG_LMAX_PRECOND',          par_int=cpar%cg_lmax_precond)
    call get_parameter(paramfile, 'CG_MAXITER',               par_int=cpar%cg_maxiter)
    call get_parameter(paramfile, 'CG_TOLERANCE',             par_dp=cpar%cg_tol)

    call get_parameter(paramfile, 'T_CMB',                    par_dp=cpar%T_cmb)
    call get_parameter(paramfile, 'MJysr_convention',         par_string=cpar%MJysr_convention)

  end subroutine read_global_params


  subroutine read_data_params(paramfile, cpar)
    implicit none

    character(len=*),  intent(in)    :: paramfile
    type(comm_params), intent(inout) :: cpar

    integer(i4b)     :: i, j, n
    character(len=2) :: itext, jtext
    
    call get_parameter(paramfile, 'NUMBAND',        par_int=cpar%numband)
    call get_parameter(paramfile, 'DATA_DIRECTORY', par_string=cpar%datadir)

    n = cpar%numband
    allocate(cpar%ds_active(n), cpar%ds_label(n))
    allocate(cpar%ds_polarization(n), cpar%ds_nside(n), cpar%ds_lmax(n))
    allocate(cpar%ds_unit(n), cpar%ds_noise_format(n), cpar%ds_mapfile(n))
    allocate(cpar%ds_noise_rms(n), cpar%ds_maskfile(n), cpar%ds_maskfile_calib(n))
    allocate(cpar%ds_samp_monopole(n), cpar%ds_samp_dipole(n), cpar%ds_samp_noiseamp(n))
    allocate(cpar%ds_bptype(n), cpar%ds_nu_c(n), cpar%ds_bpfile(n), cpar%ds_bpmodel(n))
    allocate(cpar%ds_co2t(n), cpar%ds_period(n), cpar%ds_beamtype(n), cpar%ds_blfile(n))
    allocate(cpar%ds_pixwin(n))
    allocate(cpar%ds_sample_gain(n), cpar%ds_gain_calib_comp(n), cpar%ds_gain_lmax(n))
    allocate(cpar%ds_gain_lmin(n), cpar%ds_gain_apodmask(n), cpar%ds_gain_fwhm(n))
    allocate(cpar%ds_numtemp(n), cpar%ds_tempname(n,10), cpar%ds_samptemp(n,10))
    allocate(cpar%ds_defaults(n,6), cpar%ds_temp_defaults(n,10))
    do i = 1, n
       call int2string(i, itext)
       call get_parameter(paramfile, 'BAND_ACTIVE'//itext,          par_lgt=cpar%ds_active(i))
       call get_parameter(paramfile, 'BAND_OBS_PERIOD'//itext,      par_int=cpar%ds_period(i))
       call get_parameter(paramfile, 'BAND_LABEL'//itext,           par_string=cpar%ds_label(i))
       call get_parameter(paramfile, 'BAND_POLARIZATION'//itext,    par_lgt=cpar%ds_polarization(i))
       call get_parameter(paramfile, 'BAND_NSIDE'//itext,           par_int=cpar%ds_nside(i))
       call get_parameter(paramfile, 'BAND_LMAX'//itext,            par_int=cpar%ds_lmax(i))
       call get_parameter(paramfile, 'BAND_UNIT'//itext,            par_string=cpar%ds_unit(i))
       call get_parameter(paramfile, 'BAND_NOISE_FORMAT'//itext,    par_string=cpar%ds_noise_format(i))
       call get_parameter(paramfile, 'BAND_MAPFILE'//itext,         par_string=cpar%ds_mapfile(i))
       call get_parameter(paramfile, 'BAND_NOISE_RMS'//itext,       par_string=cpar%ds_noise_rms(i))
       call get_parameter(paramfile, 'BAND_MASKFILE'//itext,        par_string=cpar%ds_maskfile(i))
       call get_parameter(paramfile, 'BAND_MASKFILE_CALIB'//itext,  par_string=cpar%ds_maskfile_calib(i))
       call get_parameter(paramfile, 'BAND_BEAMTYPE'//itext,        par_string=cpar%ds_beamtype(i))
       call get_parameter(paramfile, 'BAND_BEAM_B_L_FILE'//itext,   par_string=cpar%ds_blfile(i))
       call get_parameter(paramfile, 'BAND_PIXEL_WINDOW'//itext,    par_string=cpar%ds_pixwin(i))
       call get_parameter(paramfile, 'BAND_SAMP_MONOPOLE'//itext,   par_lgt=cpar%ds_samp_monopole(i))
       call get_parameter(paramfile, 'BAND_SAMP_DIPOLE'//itext,     par_lgt=cpar%ds_samp_dipole(i))
       call get_parameter(paramfile, 'BAND_SAMP_NOISE_AMP'//itext,  par_lgt=cpar%ds_samp_noiseamp(i))
       call get_parameter(paramfile, 'BAND_BANDPASS_TYPE'//itext,   par_string=cpar%ds_bptype(i))
       call get_parameter(paramfile, 'BAND_NOMINAL_FREQ'//itext,    par_dp=cpar%ds_nu_c(i))
       call get_parameter(paramfile, 'BAND_BANDPASSFILE'//itext,    par_string=cpar%ds_bpfile(i))
       call get_parameter(paramfile, 'BAND_BANDPASS_MODEL'//itext,  par_string=cpar%ds_bpmodel(i))
       call get_parameter(paramfile, 'BAND_CO2T'//itext,            par_dp=cpar%ds_co2t(i))
       call get_parameter(paramfile, 'BAND_SAMP_GAIN'//itext,       par_lgt=cpar%ds_sample_gain(i))
       call get_parameter(paramfile, 'BAND_GAIN_CALIB_COMP'//itext, par_string=cpar%ds_gain_calib_comp(i))
       call get_parameter(paramfile, 'BAND_GAIN_LMAX'//itext,       par_int=cpar%ds_gain_lmax(i))
       call get_parameter(paramfile, 'BAND_GAIN_LMIN'//itext,       par_int=cpar%ds_gain_lmin(i))
       call get_parameter(paramfile, 'BAND_GAIN_APOD_MASK'//itext,  par_string=cpar%ds_gain_apodmask(i))
       call get_parameter(paramfile, 'BAND_GAIN_APOD_FWHM'//itext,  par_string=cpar%ds_gain_fwhm(i))
       call get_parameter(paramfile, 'BAND_NUM_TEMPLATES'//itext,   par_int=cpar%ds_numtemp(i))
       call get_parameter(paramfile, 'BAND_DEFAULT_GAIN'//itext,    par_dp=cpar%ds_defaults(i,GAIN))
       call get_parameter(paramfile, 'BAND_DEFAULT_NOISEAMP'//itext,par_dp=cpar%ds_defaults(i,NOISEAMP))
       call get_parameter(paramfile, 'BAND_DEFAULT_MONOPOLE'//itext,par_dp=cpar%ds_defaults(i,MONOPOLE))
       call get_parameter(paramfile, 'BAND_DEFAULT_X_DIPOLE'//itext,par_dp=cpar%ds_defaults(i,DIPOLE_X))
       call get_parameter(paramfile, 'BAND_DEFAULT_Y_DIPOLE'//itext,par_dp=cpar%ds_defaults(i,DIPOLE_Y))
       call get_parameter(paramfile, 'BAND_DEFAULT_Z_DIPOLE'//itext,par_dp=cpar%ds_defaults(i,DIPOLE_Z))
       do j = 1, cpar%ds_numtemp(i)
          call int2string(j, jtext)
          call get_parameter(paramfile, 'BAND_TEMPLATE'//itext//'_'//jtext,    par_string=cpar%ds_tempname(i,j))
          call get_parameter(paramfile, 'BAND_SAMP_TEMP'//itext//'_'//jtext,   par_lgt=cpar%ds_samptemp(i,j))
          call get_parameter(paramfile, 'BAND_DEFAULT_TEMP'//itext//'_'//jtext,par_dp=cpar%ds_temp_defaults(i,j))
       end do
    end do

    ! Convert to proper internal units where necessary
    cpar%ds_nu_c = cpar%ds_nu_c * 1d9  ! Input is GHz; internal is Hz
    
  end subroutine read_data_params


  subroutine read_component_params(paramfile, cpar)
    implicit none

    character(len=*),  intent(in)    :: paramfile
    type(comm_params), intent(inout) :: cpar

    integer(i4b)     :: i, n
    character(len=2) :: itext
    
    call get_parameter(paramfile, 'NUM_SIGNAL_COMPONENTS', par_int=cpar%cs_ncomp)

    n = cpar%cs_ncomp
    allocate(cpar%cs_include(n), cpar%cs_label(n), cpar%cs_type(n), cpar%cs_class(n))
    allocate(cpar%cs_polarization(n), cpar%cs_nside(n), cpar%cs_lmax_amp(n), cpar%cs_lmax_ind(n))
    allocate(cpar%cs_unit(n), cpar%cs_nu_ref(n), cpar%cs_cltype(n))
    allocate(cpar%cs_samp_cls(n), cpar%cs_clfile(n), cpar%cs_binfile(n))
    allocate(cpar%cs_lpivot(n), cpar%cs_mask(n), cpar%cs_amp_def(n))
    allocate(cpar%cs_filedef(n))
    allocate(cpar%cs_theta_def(MAXPAR,n), cpar%cs_p_uni(MAXPAR,2,n), cpar%cs_p_gauss(MAXPAR,2,n))
    do i = 1, n
       call int2string(i, itext)
       call get_parameter(paramfile, 'INCLUDE_COMP'//itext,         par_lgt=cpar%cs_include(i))
       call get_parameter(paramfile, 'COMP_LABEL'//itext,           par_string=cpar%cs_label(i))
       call get_parameter(paramfile, 'COMP_TYPE'//itext,            par_string=cpar%cs_type(i))
       call get_parameter(paramfile, 'COMP_CLASS'//itext,           par_string=cpar%cs_class(i))
       call get_parameter(paramfile, 'COMP_POLARIZATION'//itext,    par_lgt=cpar%cs_polarization(i))
       call get_parameter(paramfile, 'COMP_NSIDE'//itext,           par_int=cpar%cs_nside(i))
       call get_parameter(paramfile, 'COMP_LMAX_AMP'//itext,        par_int=cpar%cs_lmax_amp(i))
       call get_parameter(paramfile, 'COMP_LMAX_IND'//itext,        par_int=cpar%cs_lmax_ind(i))
       call get_parameter(paramfile, 'COMP_UNIT'//itext,            par_string=cpar%cs_unit(i))
       call get_parameter(paramfile, 'COMP_NU_REF'//itext,          par_dp=cpar%cs_nu_ref(i))
       call get_parameter(paramfile, 'COMP_CL_TYPE'//itext,         par_string=cpar%cs_cltype(i))
       call get_parameter(paramfile, 'COMP_SAMP_CLS'//itext,        par_lgt=cpar%cs_samp_cls(i))
       if (trim(cpar%cs_cltype(i)) == 'binned') then
          call get_parameter(paramfile, 'COMP_CL_BIN_FILE'//itext,     par_string=cpar%cs_binfile(i))
       end if
       if (trim(cpar%cs_cltype(i)) == 'power_law') then
          call get_parameter(paramfile, 'COMP_CL_L_PIVOT'//itext,      par_int=cpar%cs_lpivot(i))
          call get_parameter(paramfile, 'COMP_DEFAULT_CL_AMP'//itext,  par_dp=cpar%cs_amp_def(i))
       end if
       if (trim(cpar%cs_cltype(i)) == 'single_l' .or. trim(cpar%cs_cltype(i)) == 'binned') then
          call get_parameter(paramfile, 'COMP_DEFAULT_CL_FILE'//itext, par_string=cpar%cs_filedef(i))
       end if
       call get_parameter(paramfile, 'COMP_MASK'//itext,            par_string=cpar%cs_mask(i))

       if (trim(cpar%cs_type(i)) == 'power_law') then
          call get_parameter(paramfile, 'COMP_DEFAULT_BETA'//itext,  par_dp=cpar%cs_theta_def(i,1))
          call get_parameter(paramfile, 'COMP_PRIOR_UNI_BETA_LOW'//itext,  par_dp=cpar%cs_p_uni(i,1,1))
          call get_parameter(paramfile, 'COMP_PRIOR_UNI_BETA_HIGH'//itext,  par_dp=cpar%cs_p_uni(i,2,1))
          call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_BETA_MEAN'//itext,  par_dp=cpar%cs_p_gauss(i,1,1))
          call get_parameter(paramfile, 'COMP_PRIOR_GAUSS_BETA_RMS'//itext,  par_dp=cpar%cs_p_gauss(i,2,1))
       end if
       
    end do

    ! Convert to proper units
    cpar%cs_nu_ref = cpar%cs_nu_ref * 1d9

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
  


end module comm_param_mod
