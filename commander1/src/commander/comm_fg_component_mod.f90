module comm_fg_component_mod
  use healpix_types
  use comm_utils
  use spline_1D_mod
  use spline_2D_mod
  use sort_utils
  use comm_bp_mod
  implicit none

  integer(i4b), parameter          :: num_recent_point = 3
  integer(i4b), parameter          :: N_RMS_MARG       = 7

  type fg_par_comp
     real(dp), allocatable, dimension(:,:) :: p
  end type fg_par_comp

  type fg_params
     type(fg_par_comp), allocatable, dimension(:) :: comp
  end type fg_params

  type fg_region
     integer(i4b) :: n, n_ext
     integer(i4b), allocatable, dimension(:,:) :: pix     ! (pix,spec)
     integer(i4b), allocatable, dimension(:,:) :: pix_ext ! (pix,spec)
     real(dp),     allocatable, dimension(:)   :: w       ! (pix); smooth region weight
     real(dp),                  dimension(num_recent_point):: recent
  end type fg_region

  type fg_reg_set
     integer(i4b) :: n
     integer(i4b), allocatable, dimension(:) :: r
  end type fg_reg_set
  
  type ind_region
     integer(i4b)                                :: nreg, nregset
     logical(lgt)                                :: independent_pixels
     type(fg_region),  allocatable, dimension(:) :: regions
     type(fg_reg_set), allocatable, dimension(:) :: regset
  end type ind_region
  
  type fg_meta_data
     character(len=24)                             :: type, init_mode, label, amp_unit
     character(len=12), dimension(30)              :: indlabel, ind_unit, ttype
     integer(i4b)                                  :: npar, ref_band, fg_id
     logical(lgt)                                  :: apply_jeffreys_prior
     logical(lgt)                                  :: output_freq_comp_maps
     logical(lgt)                                  :: enforce_positive_amplitude
     logical(lgt)                                  :: sample_amplitudes
     real(dp)                                      :: nu_ref, C_def
     integer(i4b), allocatable, dimension(:)       :: co_band
     real(dp), allocatable, dimension(:,:)         :: par
     real(dp), allocatable, dimension(:)           :: fwhm_p, p_rms
     real(dp), allocatable, dimension(:,:)         :: priors, p_rms_gauss, p_rms_uni
     real(dp), allocatable, dimension(:,:)         :: mask, indmask, priormask
     real(dp), allocatable, dimension(:,:)         :: gauss_prior
     real(dp), allocatable, dimension(:,:,:)       :: S_1D
     real(dp), allocatable, dimension(:,:,:,:,:)   :: S_2D
     real(dp)                                      :: nu_peak
     real(dp)                                      :: nu_flat, frac_flat
     real(dp)                                      :: dbeta, nu_break
     real(dp), allocatable, dimension(:)           :: us
     real(dp), allocatable, dimension(:,:)         :: S_nu_ref
     real(dp), allocatable, dimension(:,:)         :: S_phys_dust
     real(dp), allocatable, dimension(:,:,:,:)     :: S_dust_coeff
     real(dp), pointer,     dimension(:)           :: S_tabulated
     logical(lgt), pointer, dimension(:)           :: sample_S_tab   
     logical(lgt), allocatable, dimension(:)       :: optimize_prior
     type(ind_region), allocatable, dimension(:)   :: indregs
  end type fg_meta_data

  integer(i4b)                     :: num_fg_comp, num_fg_par
  integer(i4b), parameter, private :: numgrid = 100
  integer(i4b),            private :: nmaps, nside, npix
  logical(lgt),            private :: polarization, sample_inside_mask, sample_T_modes
  integer(i4b),            private :: comm_chain, myid_chain, numprocs_chain, root, myid, numprocs

  real(dp) :: alpha_flat, gamma_flat

  type(fg_meta_data), allocatable, dimension(:)            :: fg_components

  real(sp), parameter :: missval = -1.6375e30


contains


  subroutine init_comm_fg_component_mod(paramfile, comm_chain_in)
    implicit none

    integer(i4b),                         intent(in) :: comm_chain_in
    character(len=128),                   intent(in) :: paramfile

    real(dp)           :: beta, f_eff, nu, dnu, S_ref, delta, norm
    integer(i4b)       :: i, j, k, l, m, n, ierr, num_samp_point, ind(1), unit
    character(len=2)   :: i_text, j_text
    logical(lgt)       :: exist, apply_par_rms_smooth, optimize_priors
    character(len=128) :: paramname, filename_spectrum, filename, chaindir, maskname
    character(len=128) :: filename_emission, filename_u
    character(len=5000) :: line
    real(dp), allocatable, dimension(:)       :: freq
    real(dp), pointer,     dimension(:)       :: u_array
    real(dp), pointer,     dimension(:,:)     :: spectrum
    real(dp), allocatable, dimension(:,:)     :: mask
    real(dp), allocatable, dimension(:,:,:,:) :: my_grid
    integer(i4b), dimension(MPI_STATUS_SIZE)  :: status

    root       = 0
    comm_chain = comm_chain_in
    unit       = getlun()
    call mpi_comm_rank(comm_chain, myid_chain, ierr)
    call mpi_comm_size(comm_chain, numprocs_chain, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)

    call get_parameter(paramfile,    'NSIDE',                     par_int=nside)
    call get_parameter(paramfile,    'POLARIZATION',              par_lgt=polarization)
    call get_parameter(paramfile,    'NUM_FG_SIGNAL_COMPONENTS',  par_int=num_fg_comp)
    call get_parameter(paramfile,    'CHAIN_DIRECTORY',           par_string=chaindir)
    call get_parameter(paramfile,    'SAMPLE_INSIDE_MASK',        par_lgt=sample_inside_mask)
    call get_parameter(paramfile,    'APPLY_PAR_RMS_SMOOTHING',   par_lgt=apply_par_rms_smooth)
    call get_parameter(paramfile,    'OPTIMIZE_PRIORS',           par_lgt=optimize_priors)
    call get_parameter(paramfile,    'SAMPLE_TT_MODES',           par_lgt=sample_T_modes)
!    call get_parameter(paramfile,    'ALPHA_FLAT',                par_dp=alpha_flat)
!    call get_parameter(paramfile,    'GAMMA_FLAT',                par_dp=gamma_flat)
    npix       = 12*nside**2

    if (polarization) then
       nmaps = 3
    else
       nmaps = 1
    end if

    call get_parameter(paramfile,    'MASKFILE', par_string=filename)
    allocate(mask(0:12*nside**2-1,nmaps))
    call read_map(filename, mask)
    where (mask < 0.5d0)
       mask = 0.d0
    elsewhere 
       mask = 1.d0
    end where

    allocate(fg_components(num_fg_comp))

    ! Read foreground parameters
    num_fg_par = 0
    do i = 1, num_fg_comp
       
       call int2string(i, i_text)

       fg_components(i)%fg_id = i

       paramname = 'COMP_TYPE' // i_text
       call get_parameter(paramfile, paramname, par_string=fg_components(i)%type)
       
       paramname = 'FG_LABEL' // i_text
       call get_parameter(paramfile, paramname, par_string=fg_components(i)%label)       

       paramname = 'REFERENCE_FREQUENCY' // i_text
       call get_parameter(paramfile, paramname, par_dp=fg_components(i)%nu_ref)       
       fg_components(i)%nu_ref = fg_components(i)%nu_ref * 1.d9

       paramname = 'REFERENCE_BAND' // i_text
       call get_parameter(paramfile, paramname, par_int=fg_components(i)%ref_band)       

       paramname = 'FG_UNIT' // i_text
       call get_parameter(paramfile, paramname, par_string=fg_components(i)%amp_unit)       

       paramname = 'OUTPUT_FREQUENCY_COMPONENT_MAPS' // i_text
       call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%output_freq_comp_maps)

       paramname = 'ENFORCE_POSITIVE_AMPLITUDE' // i_text
       call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%enforce_positive_amplitude)       

       paramname = 'SAMPLE_AMPLITUDES' // i_text
       call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%sample_amplitudes)

       paramname = 'COMPONENT_MASK' // i_text
       call get_parameter(paramfile, paramname, par_string=maskname)
       allocate(fg_components(i)%mask(0:npix-1,nmaps))
       if (trim(maskname) == 'fullsky') then
          fg_components(i)%mask = 1.d0
       else
          call read_map(maskname, fg_components(i)%mask)
          where (fg_components(i)%mask < 0.5d0)
             fg_components(i)%mask = 0.d0
          elsewhere 
             fg_components(i)%mask = 1.d0
          end where
       end if
       if (.not. sample_inside_mask) fg_components(i)%mask = fg_components(i)%mask * mask

       paramname = 'INDEX_SAMPLING_MASK' // i_text
       call get_parameter(paramfile, paramname, par_string=maskname)
       allocate(fg_components(i)%indmask(0:npix-1,nmaps))
       if (trim(maskname) == 'fullsky') then
          fg_components(i)%indmask = 1.d0
       else
          call read_map(maskname, fg_components(i)%indmask)
          where (fg_components(i)%indmask < 0.5d0)
             fg_components(i)%indmask = 0.d0
          elsewhere 
             fg_components(i)%indmask = 1.d0
          end where
       end if
       fg_components(i)%indmask = fg_components(i)%indmask * fg_components(i)%mask

       if (optimize_priors) then
          paramname = 'PRIOR_OPTIMIZATION_MASK' // i_text
          call get_parameter(paramfile, paramname, par_string=maskname)
          allocate(fg_components(i)%priormask(0:npix-1,nmaps))
          if (trim(maskname) == 'fullsky') then
             fg_components(i)%priormask = 1.d0
          else
             call read_map(maskname, fg_components(i)%priormask)
             where (fg_components(i)%priormask < 0.5d0)
                fg_components(i)%priormask = 0.d0
             elsewhere 
                fg_components(i)%priormask = 1.d0
             end where
          end if
          if (.not. sample_inside_mask) fg_components(i)%priormask = fg_components(i)%priormask * mask
       end if

       allocate(fg_components(i)%sample_S_tab(numband))
       allocate(fg_components(i)%S_tabulated(numband))
       fg_components(i)%sample_S_tab = .false.
       fg_components(i)%S_tabulated  = 1.d0   ! Set default to 1, so one can always multiply with it


       if (trim(fg_components(i)%type) == 'one-component_dust') then

          paramname = 'APPLY_JEFFREYS_PRIOR' // i_text
          call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%apply_jeffreys_prior)       

          num_fg_par                   = num_fg_par + 2
          fg_components(i)%indlabel(1) = 'beta'
          fg_components(i)%indlabel(2) = 'Td'
          fg_components(i)%ind_unit(1) = ''
          fg_components(i)%ind_unit(2) = 'K'
          fg_components(i)%ttype(1)    = 'Beta'
          fg_components(i)%ttype(2)    = 'Td'
          fg_components(i)%npar        = 2
          allocate(fg_components(i)%priors(fg_components(i)%npar,3))
          allocate(fg_components(i)%gauss_prior(2,2))
          allocate(fg_components(i)%par(numgrid,2))
          allocate(fg_components(i)%S_2D(4,4,numgrid,numgrid,numband))

          paramname = 'INITIALIZATION_MODE' // i_text
          call get_parameter(paramfile, paramname, par_string=fg_components(i)%init_mode)
             
          paramname = 'NU_FLATTENING' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%nu_flat)
          fg_components(i)%nu_flat = 1.d9 * fg_components(i)%nu_flat

          paramname = 'FRAC_FLATTENING' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%frac_flat)

          paramname = 'DEFAULT_EMISSIVITY' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(1,3))
          paramname = 'DEFAULT_DUST_TEMP' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(2,3))
          
          paramname = 'EMISSIVITY_PRIOR_UNIFORM_LOW'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(1,1))
          paramname = 'EMISSIVITY_PRIOR_UNIFORM_HIGH' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(1,2))
          if (fg_components(i)%priors(1,1) >= fg_components(i)%priors(1,2)) then
             write(*,*) 'Error: Lower prior is larger than or equal to upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          paramname = 'EMISSIVITY_PRIOR_GAUSSIAN_MEAN'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(1,1))
          paramname = 'EMISSIVITY_PRIOR_GAUSSIAN_STDDEV' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(1,2))
          
          paramname = 'DUST_TEMP_PRIOR_UNIFORM_LOW'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(2,1))
          paramname = 'DUST_TEMP_PRIOR_UNIFORM_HIGH' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(2,2))
          if (fg_components(i)%priors(2,1) >= fg_components(i)%priors(2,2)) then
             write(*,*) 'Error: Lower prior is larger than or equal to upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          paramname = 'DUST_TEMP_PRIOR_GAUSSIAN_MEAN'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(2,1))
          paramname = 'DUST_TEMP_PRIOR_GAUSSIAN_STDDEV' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(2,2))

       else if (trim(fg_components(i)%type) == 'physical_dust') then

          paramname = 'APPLY_JEFFREYS_PRIOR' // i_text
          call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%apply_jeffreys_prior)       

          num_fg_par                   = num_fg_par + 1
          fg_components(i)%indlabel(1) = 'U'
          fg_components(i)%ind_unit(1) = 'erg/s/sr/cm^2'
          fg_components(i)%ttype(1)    = 'Intensity'
          fg_components(i)%npar        = 1
          allocate(fg_components(i)%priors(1,3))
          allocate(fg_components(i)%gauss_prior(1,2))
          allocate(fg_components(i)%par(numgrid,1))
          allocate(fg_components(i)%S_1D(numgrid,2,numband))

          call get_parameter(paramfile, 'INITIALIZATION_MODE' // i_text, par_string=fg_components(i)%init_mode)
          call get_parameter(paramfile, 'DEFAULT_U' // i_text, par_dp=fg_components(i)%priors(1,3))

          call get_parameter(paramfile, 'U_PRIOR_UNIFORM_LOW' // i_text, &
               & par_dp=fg_components(i)%priors(1,1))
          call get_parameter(paramfile, 'U_PRIOR_UNIFORM_HIGH' // i_text, &
               & par_dp=fg_components(i)%priors(1,2))

          if (fg_components(i)%priors(1,3) .lt. -3 .or. fg_components(i)%priors(1,3) .gt. 6) then
             write(*,*) 'Error: Prior outside of sampling range for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          if (fg_components(i)%priors(1,1) >= fg_components(i)%priors(1,2)) then
             write(*,*) 'Error: Lower prior is larger than or equal to upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if

          call get_parameter(paramfile, 'U_PRIOR_GAUSSIAN_MEAN' // &
               & i_text, par_dp=fg_components(i)%gauss_prior(1,1))
          call get_parameter(paramfile, 'U_PRIOR_GAUSSIAN_STDDEV' // &
               & i_text, par_dp=fg_components(i)%gauss_prior(1,2))

          paramname = 'EMISSION_FILENAME' // i_text
          call get_parameter(paramfile, paramname, par_string=filename_emission)
          paramname = 'U_FILENAME' // i_text
          call get_parameter(paramfile, paramname, par_string=filename_u)
          call read_phys_dust_spectrum(filename_emission, filename_u, spectrum, u_array, fg_components(i)%nu_ref)

          allocate(fg_components(i)%S_phys_dust(size(spectrum,1),size(spectrum,2)))
          allocate(fg_components(i)%S_dust_coeff(4,4,size(spectrum,1),size(spectrum,2)))
          allocate(fg_components(i)%us(size(spectrum,1)-1))

          fg_components(i)%S_phys_dust(:,:)  = spectrum
          fg_components(i)%us(:)             = u_array(:)

          deallocate(spectrum)
          deallocate(u_array)

          call splie2_full_precomp(fg_components(i)%us(:), log(fg_components(i)%S_phys_dust(1,:)), &
               & fg_components(i)%S_phys_dust(2:,:), fg_components(i)%S_dust_coeff)

       else if (trim(fg_components(i)%type) == 'freefree_EM') then

          paramname = 'APPLY_JEFFREYS_PRIOR' // i_text
          call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%apply_jeffreys_prior)       

          num_fg_par                   = num_fg_par + 2
          fg_components(i)%indlabel(1) = 'EM'
          fg_components(i)%indlabel(2) = 'T_e'
          fg_components(i)%ind_unit(1) = 'pc cm^-6'
          fg_components(i)%ind_unit(2) = 'K'
          fg_components(i)%ttype(1)    = 'EM'
          fg_components(i)%ttype(2)    = 'Te'
          fg_components(i)%npar        = 2
          allocate(fg_components(i)%priors(fg_components(i)%npar,3))
          allocate(fg_components(i)%gauss_prior(2,2))
          allocate(fg_components(i)%par(numgrid,2))
          allocate(fg_components(i)%S_2D(4,4,numgrid,numgrid,numband))

          paramname = 'INITIALIZATION_MODE' // i_text
          call get_parameter(paramfile, paramname, par_string=fg_components(i)%init_mode)
             
          paramname = 'DEFAULT_EM' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(1,3))
          paramname = 'DEFAULT_T_E' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(2,3))
          
          paramname = 'EM_PRIOR_UNIFORM_LOW'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(1,1))
          paramname = 'EM_PRIOR_UNIFORM_HIGH' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(1,2))
          if (fg_components(i)%priors(1,1) >= fg_components(i)%priors(1,2)) then
             write(*,*) 'Error: Lower prior is larger than or equal to upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          paramname = 'EM_PRIOR_GAUSSIAN_MEAN'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(1,1))
          paramname = 'EM_PRIOR_GAUSSIAN_STDDEV' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(1,2))
          
          paramname = 'T_E_PRIOR_UNIFORM_LOW'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(2,1))
          paramname = 'T_E_PRIOR_UNIFORM_HIGH' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(2,2))
          if (fg_components(i)%priors(2,1) >= fg_components(i)%priors(2,2)) then
             write(*,*) 'Error: Lower prior is larger than or equal to upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          paramname = 'T_E_PRIOR_GAUSSIAN_MEAN'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(2,1))
          paramname = 'T_E_PRIOR_GAUSSIAN_STDDEV' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(2,2))

          fg_components(i)%priors(1,:)      = log(fg_components(i)%priors(1,:))
          fg_components(i)%gauss_prior(1,:) = log(fg_components(i)%gauss_prior(1,:))

       else if (trim(fg_components(i)%type) == 'power_law') then

          paramname = 'APPLY_JEFFREYS_PRIOR' // i_text
          call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%apply_jeffreys_prior)       

          num_fg_par                   = num_fg_par + 2
          fg_components(i)%npar        = 2
          fg_components(i)%indlabel(1) = 'beta'
          fg_components(i)%indlabel(2) = 'C'
          fg_components(i)%ind_unit(1) = ''
          fg_components(i)%ind_unit(2) = ''
          fg_components(i)%ttype(1)    = 'Beta'
          fg_components(i)%ttype(2)    = 'C'
          allocate(fg_components(i)%priors(fg_components(i)%npar,3))
          allocate(fg_components(i)%gauss_prior(2,2))
          allocate(fg_components(i)%par(numgrid,2))
          allocate(fg_components(i)%S_2D(4,4,numgrid,numgrid,numband))

          paramname = 'INITIALIZATION_MODE' // i_text
          call get_parameter(paramfile, paramname, par_string=fg_components(i)%init_mode)
             
          paramname = 'DEFAULT_SPECTRAL_INDEX' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(1,3))

          paramname = 'DEFAULT_CURVATURE' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(2,3))

          paramname = 'INDEX_PRIOR_UNIFORM_LOW'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(1,1))
          paramname = 'INDEX_PRIOR_UNIFORM_HIGH' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(1,2))
          if (fg_components(i)%priors(1,1) > fg_components(i)%priors(1,2)) then
             write(*,*) 'Error: Lower prior is larger than upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          paramname = 'INDEX_PRIOR_GAUSSIAN_MEAN'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(1,1))
          paramname = 'INDEX_PRIOR_GAUSSIAN_STDDEV' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(1,2))
          
          paramname = 'CURVATURE_PRIOR_UNIFORM_LOW'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(2,1))
          paramname = 'CURVATURE_PRIOR_UNIFORM_HIGH' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(2,2))
          if (fg_components(i)%priors(2,1) > fg_components(i)%priors(2,2)) then
             write(*,*) 'Error: Lower prior is larger than upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          paramname = 'CURVATURE_PRIOR_GAUSSIAN_MEAN'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(2,1))
          paramname = 'CURVATURE_PRIOR_GAUSSIAN_STDDEV' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(2,2))

       else if (trim(fg_components(i)%type) == 'power_law_break') then

          paramname = 'APPLY_JEFFREYS_PRIOR' // i_text
          call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%apply_jeffreys_prior)       

          num_fg_par                   = num_fg_par + 2
          fg_components(i)%npar        = 2
          fg_components(i)%indlabel(1) = 'beta'
          fg_components(i)%indlabel(2) = 'dbeta'
          fg_components(i)%ind_unit(1) = ''
          fg_components(i)%ind_unit(2) = ''
          fg_components(i)%ttype(1)    = 'Beta'
          fg_components(i)%ttype(2)    = 'Beta'
          allocate(fg_components(i)%priors(fg_components(i)%npar,3))
          allocate(fg_components(i)%gauss_prior(2,2))
          allocate(fg_components(i)%par(numgrid,2))
          allocate(fg_components(i)%S_2D(4,4,numgrid,numgrid,numband))

          paramname = 'INITIALIZATION_MODE' // i_text
          call get_parameter(paramfile, paramname, par_string=fg_components(i)%init_mode)
             
          paramname = 'NU_BREAK' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%nu_break)       

          paramname = 'DEFAULT_BETA' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(1,3))

          paramname = 'DEFAULT_DBETA' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(2,3))

          paramname = 'BETA_PRIOR_UNIFORM_LOW'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(1,1))
          paramname = 'BETA_PRIOR_UNIFORM_HIGH' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(1,2))
          if (fg_components(i)%priors(1,1) > fg_components(i)%priors(1,2)) then
             write(*,*) 'Error: Lower prior is larger than upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          paramname = 'BETA_PRIOR_GAUSSIAN_MEAN'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(1,1))
          paramname = 'BETA_PRIOR_GAUSSIAN_STDDEV' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(1,2))
          
          paramname = 'DBETA_PRIOR_UNIFORM_LOW'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(2,1))
          paramname = 'DBETA_PRIOR_UNIFORM_HIGH' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(2,2))
          if (fg_components(i)%priors(2,1) > fg_components(i)%priors(2,2)) then
             write(*,*) 'Error: Lower prior is larger than upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          paramname = 'DBETA_PRIOR_GAUSSIAN_MEAN'  // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(2,1))
          paramname = 'DBETA_PRIOR_GAUSSIAN_STDDEV' // i_text
          call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(2,2))

       else if (trim(fg_components(i)%type) == 'freefree') then

          paramname = 'APPLY_JEFFREYS_PRIOR' // i_text
          call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%apply_jeffreys_prior)       

          num_fg_par                    = num_fg_par + 1
          fg_components(i)%npar         = 1
          fg_components(i)%indlabel(1)  = 'T_e'
          fg_components(i)%ind_unit(1)  = 'K'
          fg_components(i)%ttype(1)     = 'Te'
          allocate(fg_components(i)%priors(1,3))
          allocate(fg_components(i)%gauss_prior(1,2))
          allocate(fg_components(i)%par(numgrid,1))
          allocate(fg_components(i)%S_1D(numgrid,2,numband))

          call get_parameter(paramfile, 'INITIALIZATION_MODE' // i_text, par_string=fg_components(i)%init_mode)
          call get_parameter(paramfile, 'DEFAULT_T_E' // i_text, par_dp=fg_components(i)%priors(1,3))

          call get_parameter(paramfile, 'T_E_PRIOR_UNIFORM_LOW'// i_text, par_dp=fg_components(i)%priors(1,1))
          call get_parameter(paramfile, 'T_E_PRIOR_UNIFORM_HIGH'//i_text, par_dp=fg_components(i)%priors(1,2))
          if (fg_components(i)%priors(1,1) >= fg_components(i)%priors(1,2)) then
             write(*,*) 'Error: Lower prior is larger than or equal to upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          call get_parameter(paramfile, 'T_E_PRIOR_GAUSSIAN_MEAN'  // &
               & i_text, par_dp=fg_components(i)%gauss_prior(1,1))
          call get_parameter(paramfile, 'T_E_PRIOR_GAUSSIAN_STDDEV' // &
               & i_text, par_dp=fg_components(i)%gauss_prior(1,2))

       else if (trim(fg_components(i)%type) == 'CO_multiline') then

          paramname = 'APPLY_JEFFREYS_PRIOR' // i_text
          call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%apply_jeffreys_prior)       

          paramname = 'INITIALIZATION_MODE' // i_text
          call get_parameter(paramfile, paramname, par_string=fg_components(i)%init_mode)

          paramname = 'NUM_CO_HARMONICS' // i_text
          call get_parameter(paramfile, paramname, par_int=fg_components(i)%npar)

          allocate(fg_components(i)%co_band(fg_components(i)%npar+1))
          fg_components(i)%co_band(1) = fg_components(i)%ref_band

          num_fg_par                   = num_fg_par + fg_components(i)%npar 
          allocate(fg_components(i)%priors(fg_components(i)%npar,3))
          allocate(fg_components(i)%gauss_prior(fg_components(i)%npar,2))
          do j = 1, fg_components(i)%npar
             call int2string(j, j_text)
             paramname = 'LINE_LABEL' // i_text // '_' // j_text
             call get_parameter(paramfile, paramname, par_string=fg_components(i)%indlabel(j))       
             paramname = 'BAND' // i_text // '_' // j_text
             call get_parameter(paramfile, paramname, par_int=fg_components(i)%co_band(j+1))
             paramname = 'DEFAULT_CO_LINE_RATIO_' // i_text // '_' // j_text
             call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(j,3))
             paramname = 'CO_PRIOR_UNIFORM_LOW'  // i_text // '_' // j_text
             call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(j,1))
             paramname = 'CO_PRIOR_UNIFORM_HIGH'  // i_text // '_' // j_text
             call get_parameter(paramfile, paramname, par_dp=fg_components(i)%priors(j,2))
             if (fg_components(i)%priors(j,1) > fg_components(i)%priors(j,2)) then
                write(*,*) 'Error: Lower prior is larger than upper prior for parameter no. ', i
                call mpi_finalize(ierr)
                stop
             end if
             paramname = 'CO_PRIOR_GAUSSIAN_MEAN'  // i_text // '_' // j_text
             call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(j,1))
             paramname = 'CO_PRIOR_GAUSSIAN_STDDEV'  // i_text // '_' // j_text
             call get_parameter(paramfile, paramname, par_dp=fg_components(i)%gauss_prior(j,2))
          end do

       else if (trim(fg_components(i)%type) == 'AME_freq_shift') then

          paramname = 'APPLY_JEFFREYS_PRIOR' // i_text
          call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%apply_jeffreys_prior)       

          num_fg_par                    = num_fg_par + 1
          fg_components(i)%npar         = 1
          fg_components(i)%indlabel(1)  = 'nup'
          fg_components(i)%ind_unit(1)  = 'GHz'
          fg_components(i)%ttype(1)     = 'Frequency'
          allocate(fg_components(i)%priors(1,3))
          allocate(fg_components(i)%gauss_prior(1,2))
          allocate(fg_components(i)%par(numgrid,1))
          allocate(fg_components(i)%S_1D(numgrid,2,numband))

          call get_parameter(paramfile, 'INITIALIZATION_MODE' // i_text, par_string=fg_components(i)%init_mode)
          call get_parameter(paramfile, 'DEFAULT_FREQUENCY' // i_text, par_dp=fg_components(i)%priors(1,3))

          call get_parameter(paramfile, 'FREQUENCY_PRIOR_UNIFORM_LOW'  // &
               & i_text, par_dp=fg_components(i)%priors(1,1))
          call get_parameter(paramfile, 'FREQUENCY_PRIOR_UNIFORM_HIGH' // &
               & i_text, par_dp=fg_components(i)%priors(1,2))
          if (fg_components(i)%priors(1,1) >= fg_components(i)%priors(1,2)) then
             write(*,*) 'Error: Lower prior is larger than or equal to upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          call get_parameter(paramfile, 'FREQUENCY_PRIOR_GAUSSIAN_MEAN'  // &
               & i_text, par_dp=fg_components(i)%gauss_prior(1,1))
          call get_parameter(paramfile, 'FREQUENCY_PRIOR_GAUSSIAN_STDDEV' // &
               & i_text, par_dp=fg_components(i)%gauss_prior(1,2))

          paramname = 'SPECTRUM_FILENAME' // i_text
          call get_parameter(paramfile, paramname, par_string=filename_spectrum)
          call read_spectrum(filename_spectrum, spectrum)
          ! Normalize to unity, to get rid of numerical errors
          spectrum(:,2) = spectrum(:,2) / maxval(spectrum(:,2))
          ! Assume input spectrum is in flux units; convert to antenna units by 1/nu^2
!          spectrum(:,2) = spectrum(:,2) / spectrum(:,1)**2
!          if (myid_chain == root) then
!             open(58,file='ame_int.dat')
!             write(*,*) shape(spectrum)
!             do j = 1, size(spectrum,1)
!                write(58,*) spectrum(j,1), spectrum(j,2), spectrum(j,3)
!             end do
!             close(58)
!          end if
!          call mpi_finalize(ierr)
!          stop
          allocate(fg_components(i)%S_nu_ref(size(spectrum,1),3))
          fg_components(i)%S_nu_ref(:,1:2) = spectrum
          deallocate(spectrum)
          fg_components(i)%s_nu_ref(:,2) = log(fg_components(i)%s_nu_ref(:,2))
          call spline(fg_components(i)%S_nu_ref(:,1), fg_components(i)%S_nu_ref(:,2), &
               & 1.d30, 1.d30, fg_components(i)%S_nu_ref(:,3))
          ind = maxloc(fg_components(i)%S_nu_ref(:,2))
          fg_components(i)%nu_peak = fg_components(i)%S_nu_ref(ind(1),1)

       else if (trim(fg_components(i)%type) == 'AME_freq_shift_2par') then

          paramname = 'APPLY_JEFFREYS_PRIOR' // i_text
          call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%apply_jeffreys_prior)       

          num_fg_par                   = num_fg_par + 2
          fg_components(i)%indlabel(1) = 'nup'
          fg_components(i)%indlabel(2) = 'alpha'
          fg_components(i)%ind_unit(1) = 'GHz'
          fg_components(i)%ind_unit(2) = ''
          fg_components(i)%ttype(1)    = 'Frequency'
          fg_components(i)%ttype(2)    = 'Tilt'
          fg_components(i)%npar        = 2
          allocate(fg_components(i)%priors(fg_components(i)%npar,3))
          allocate(fg_components(i)%gauss_prior(2,2))
          allocate(fg_components(i)%par(numgrid,2))
          allocate(fg_components(i)%S_2D(4,4,numgrid,numgrid,numband))

          paramname = 'INITIALIZATION_MODE' // i_text
          call get_parameter(paramfile, paramname, par_string=fg_components(i)%init_mode)

          call get_parameter(paramfile, 'INITIALIZATION_MODE' // i_text, par_string=fg_components(i)%init_mode)
          call get_parameter(paramfile, 'DEFAULT_FREQUENCY' // i_text, par_dp=fg_components(i)%priors(1,3))
          call get_parameter(paramfile, 'DEFAULT_TILT' // i_text, par_dp=fg_components(i)%priors(2,3))

          call get_parameter(paramfile, 'FREQUENCY_PRIOR_UNIFORM_LOW'  // &
               & i_text, par_dp=fg_components(i)%priors(1,1))
          call get_parameter(paramfile, 'FREQUENCY_PRIOR_UNIFORM_HIGH' // &
               & i_text, par_dp=fg_components(i)%priors(1,2))
          if (fg_components(i)%priors(1,1) >= fg_components(i)%priors(1,2)) then
             write(*,*) 'Error: Lower prior is larger than or equal to upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          call get_parameter(paramfile, 'FREQUENCY_PRIOR_GAUSSIAN_MEAN'  // &
               & i_text, par_dp=fg_components(i)%gauss_prior(1,1))
          call get_parameter(paramfile, 'FREQUENCY_PRIOR_GAUSSIAN_STDDEV' // &
               & i_text, par_dp=fg_components(i)%gauss_prior(1,2))

          call get_parameter(paramfile, 'TILT_PRIOR_UNIFORM_LOW'  // &
               & i_text, par_dp=fg_components(i)%priors(2,1))
          call get_parameter(paramfile, 'TILT_PRIOR_UNIFORM_HIGH' // &
               & i_text, par_dp=fg_components(i)%priors(2,2))
          if (fg_components(i)%priors(2,1) >= fg_components(i)%priors(2,2)) then
             write(*,*) 'Error: Lower prior is larger than or equal to upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          call get_parameter(paramfile, 'TILT_PRIOR_GAUSSIAN_MEAN'  // &
               & i_text, par_dp=fg_components(i)%gauss_prior(2,1))
          call get_parameter(paramfile, 'TILT_PRIOR_GAUSSIAN_STDDEV' // &
               & i_text, par_dp=fg_components(i)%gauss_prior(2,2))

          paramname = 'SPECTRUM_FILENAME' // i_text
          call get_parameter(paramfile, paramname, par_string=filename_spectrum)
          call read_spectrum(filename_spectrum, spectrum)
          spectrum(:,2) = spectrum(:,2) / maxval(spectrum(:,2))

          allocate(fg_components(i)%S_nu_ref(size(spectrum,1),3))
          fg_components(i)%S_nu_ref(:,1:2) = spectrum
          deallocate(spectrum)
          fg_components(i)%s_nu_ref(:,2) = log(fg_components(i)%s_nu_ref(:,2))
          call spline(fg_components(i)%S_nu_ref(:,1), fg_components(i)%S_nu_ref(:,2), &
               & 1.d30, 1.d30, fg_components(i)%S_nu_ref(:,3))
          ind = maxloc(fg_components(i)%S_nu_ref(:,2))
          fg_components(i)%nu_peak = fg_components(i)%S_nu_ref(ind(1),1)

       else if (trim(fg_components(i)%type) == 'AME_lognormal_width') then

          ! This log-normal approximation is a simplification of the Stevenson 2014 spectrum.
          ! We fit just two spectral parameters, nu_peak, and the width of the spectrum w_ame.
          ! The spectrum width is essentially a combination of the three free parameters in the
          ! Stevenson 2014 approximatin. http://arxiv.org/abs/2001.07159

          paramname = 'APPLY_JEFFREYS_PRIOR' // i_text
          call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%apply_jeffreys_prior)       

          num_fg_par                    = num_fg_par + 2
          fg_components(i)%indlabel(1)  = 'nup'
          fg_components(i)%indlabel(2)  = 'width'
          fg_components(i)%ind_unit(1)  = 'GHz'
          fg_components(i)%ind_unit(2)  = ''
          fg_components(i)%ttype(1)     = 'Frequency'
          fg_components(i)%ttype(2)     = 'Width'
          fg_components(i)%npar         = 2
          allocate(fg_components(i)%priors(fg_components(i)%npar,3))
          allocate(fg_components(i)%gauss_prior(2,2))
          allocate(fg_components(i)%par(numgrid,2))
          allocate(fg_components(i)%S_2D(4,4,numgrid,numgrid,numband))

          call get_parameter(paramfile, 'INITIALIZATION_MODE' // i_text, par_string=fg_components(i)%init_mode)
          call get_parameter(paramfile, 'DEFAULT_FREQUENCY' // i_text, par_dp=fg_components(i)%priors(1,3))
          call get_parameter(paramfile, 'DEFAULT_WIDTH' // i_text, par_dp=fg_components(i)%priors(2,3))

          call get_parameter(paramfile, 'FREQUENCY_PRIOR_UNIFORM_LOW'  // &
               & i_text, par_dp=fg_components(i)%priors(1,1))
          call get_parameter(paramfile, 'FREQUENCY_PRIOR_UNIFORM_HIGH' // &
               & i_text, par_dp=fg_components(i)%priors(1,2))
          if (fg_components(i)%priors(1,1) >= fg_components(i)%priors(1,2)) then
             write(*,*) 'Error: Lower prior is larger than or equal to upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          call get_parameter(paramfile, 'FREQUENCY_PRIOR_GAUSSIAN_MEAN'  // &
               & i_text, par_dp=fg_components(i)%gauss_prior(1,1))
          call get_parameter(paramfile, 'FREQUENCY_PRIOR_GAUSSIAN_STDDEV' // &
               & i_text, par_dp=fg_components(i)%gauss_prior(1,2))

          call get_parameter(paramfile, 'WIDTH_PRIOR_UNIFORM_LOW'  // &
               & i_text, par_dp=fg_components(i)%priors(2,1))
          call get_parameter(paramfile, 'WIDTH_PRIOR_UNIFORM_HIGH' // &
               & i_text, par_dp=fg_components(i)%priors(2,2))

          if (fg_components(i)%priors(2,3) .lt. 0.1 .or. fg_components(i)%priors(2,3) .gt. 1.0) then
             write(*,*) 'Error: Prior outside of allowed width values.', i
             call mpi_finalize(ierr)
             stop
          end if

          if (fg_components(i)%priors(2,1) >= fg_components(i)%priors(2,2)) then
             write(*,*) 'Error: Lower prior is larger than or equal to upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          call get_parameter(paramfile, 'WIDTH_PRIOR_GAUSSIAN_MEAN'  // &
               & i_text, par_dp=fg_components(i)%gauss_prior(2,1))
          call get_parameter(paramfile, 'WIDTH_PRIOR_GAUSSIAN_STDDEV' // &
               & i_text, par_dp=fg_components(i)%gauss_prior(2,2))

       else if (trim(fg_components(i)%type) == 'magnetic_dust') then

          paramname = 'APPLY_JEFFREYS_PRIOR' // i_text
          call get_parameter(paramfile, paramname, par_lgt=fg_components(i)%apply_jeffreys_prior)       

          num_fg_par                    = num_fg_par + 1
          fg_components(i)%npar         = 1
          fg_components(i)%indlabel(1)  = 'T_m'
          fg_components(i)%ind_unit(1)  = 'K'
          fg_components(i)%ttype(1)     = 'Tm'
          allocate(fg_components(i)%priors(1,3))
          allocate(fg_components(i)%gauss_prior(1,2))
          allocate(fg_components(i)%par(numgrid,1))
          allocate(fg_components(i)%S_1D(numgrid,2,numband))

          call get_parameter(paramfile, 'INITIALIZATION_MODE' // i_text, par_string=fg_components(i)%init_mode)
          call get_parameter(paramfile, 'DEFAULT_T_M' // i_text, par_dp=fg_components(i)%priors(1,3))

          call get_parameter(paramfile, 'T_M_PRIOR_UNIFORM_LOW'  // &
               & i_text, par_dp=fg_components(i)%priors(1,1))
          call get_parameter(paramfile, 'T_M_PRIOR_UNIFORM_HIGH' // &
               & i_text, par_dp=fg_components(i)%priors(1,2))
          if (fg_components(i)%priors(1,1) >= fg_components(i)%priors(1,2)) then
             write(*,*) 'Error: Lower prior is larger than or equal to upper prior for parameter no. ', i
             call mpi_finalize(ierr)
             stop
          end if
          call get_parameter(paramfile, 'T_M_PRIOR_GAUSSIAN_MEAN'  // &
               & i_text, par_dp=fg_components(i)%gauss_prior(1,1))
          call get_parameter(paramfile, 'T_M_PRIOR_GAUSSIAN_STDDEV' // &
               & i_text, par_dp=fg_components(i)%gauss_prior(1,2))

          paramname = 'SPECTRUM_FILENAME' // i_text
          call get_parameter(paramfile, paramname, par_string=filename_spectrum)
          call read_spectrum(filename_spectrum, spectrum)
          ! Normalize to unity, to get rid of numerical errors
          spectrum(:,2) = spectrum(:,2) / maxval(spectrum(:,2))
          ! Assume input spectrum is in flux units; convert to antenna units by 1/nu^2
!          spectrum(:,2) = spectrum(:,2) / spectrum(:,1)**2
!          if (myid_chain == root) then
!             open(58,file='ame_int.dat')
!             write(*,*) shape(spectrum)
!             do j = 1, size(spectrum,1)
!                write(58,*) spectrum(j,1), spectrum(j,2), spectrum(j,3)
!             end do
!             close(58)
!          end if
!          call mpi_finalize(ierr)
!          stop
          allocate(fg_components(i)%S_nu_ref(size(spectrum,1),3))
          fg_components(i)%S_nu_ref(:,1:2) = spectrum
          deallocate(spectrum)
          call spline(fg_components(i)%S_nu_ref(:,1), fg_components(i)%S_nu_ref(:,2), &
               & 1.d30, 1.d30, fg_components(i)%S_nu_ref(:,3))

       else if (trim(fg_components(i)%type) == 'sz') then
             
          fg_components(i)%npar       = 0
          do m = 1, numband
             fg_components(i)%S_tabulated(m) = (2.d0*fg_components(i)%nu_ref**2*k_b/c**2 / &
               & (compute_bnu_prime_single(fg_components(i)%nu_ref) * &
               & sz_thermo(fg_components(i)%nu_ref))) / bp(m)%a2sz * ant2data(m)
          end do

       else if (trim(fg_components(i)%type) == 'cmb') then
             
          fg_components(i)%npar        = 0
          do j = 1, numband
             fg_components(i)%S_tabulated(j) = compute_ant2thermo_single(fg_components(i)%nu_ref) / bp(j)%a2t * &
                  & ant2data(j)
          end do

       else if (trim(fg_components(i)%type) == 'spectral_line') then
             
          fg_components(i)%npar        = 0
          do j = 1, numband
             fg_components(i)%S_tabulated(j) = get_bp_line_ant(j, fg_components(i)%nu_ref) * ant2data(j)
          end do
             
       else
             
          write(*,*) 'comm_fg_component_mod: Unknown spectrum type = ', trim(fg_components(i)%type)
          stop
             
       end if

       allocate(fg_components(i)%fwhm_p(fg_components(i)%npar))
       do j = 1, fg_components(i)%npar
          call int2string(j, j_text)
          call get_parameter(paramfile, 'FWHM_PAR'//i_text//'_'//j_text, par_dp=fg_components(i)%fwhm_p(j))
       end do

       allocate(fg_components(i)%optimize_prior(fg_components(i)%npar))
       if (optimize_priors) then
          do j = 1, fg_components(i)%npar
             call int2string(j, j_text)
             call get_parameter(paramfile, 'OPTIMIZE_PRIOR_MEAN'//i_text//'_'//j_text, &
                  & par_lgt=fg_components(i)%optimize_prior(j))
          end do
       else
          fg_components(i)%optimize_prior = .false.
       end if

       allocate(fg_components(i)%p_rms(fg_components(i)%npar))
       allocate(fg_components(i)%p_rms_gauss(fg_components(i)%npar,2))
       allocate(fg_components(i)%p_rms_uni(fg_components(i)%npar,2))
       do j = 1, fg_components(i)%npar
          if (.not. apply_par_rms_smooth .or. trim(fg_components(i)%type) == 'CO_multiline') then
             fg_components(i)%p_rms       = 0.d0
             fg_components(i)%p_rms_gauss = 0.d0
             fg_components(i)%p_rms_uni   = 0.d0
          else
             call int2string(j, j_text)
             call get_parameter(paramfile, 'PAR_RMS_DEFAULT'//i_text//'_'//j_text, par_dp=fg_components(i)%p_rms(j))
             call get_parameter(paramfile, 'PAR_RMS_GAUSS_MEAN'//i_text//'_'//j_text, par_dp=fg_components(i)%p_rms_gauss(j,1))
             call get_parameter(paramfile, 'PAR_RMS_GAUSS_RMS'//i_text//'_'//j_text, par_dp=fg_components(i)%p_rms_gauss(j,2))
             call get_parameter(paramfile, 'PAR_RMS_UNIFORM_LOW'//i_text//'_'//j_text, par_dp=fg_components(i)%p_rms_uni(j,1))
             call get_parameter(paramfile, 'PAR_RMS_UNIFORM_HIGH'//i_text//'_'//j_text, par_dp=fg_components(i)%p_rms_uni(j,2))
          end if
       end do

       do j = 1, fg_components(i)%npar
          if (trim(fg_components(i)%type) == 'CO_multiline') cycle
          ! Give some extra space to ensure stable splines
          delta = fg_components(i)%priors(j,2)-fg_components(i)%priors(j,1)
          if (.not. trim(fg_components(i)%type) == 'physical_dust') then
             fg_components(i)%priors(j,1) = fg_components(i)%priors(j,1) - 0.1d0*delta
             fg_components(i)%priors(j,2) = fg_components(i)%priors(j,2) + 0.1d0*delta
          end if
          do k = 1, numgrid
             fg_components(i)%par(k,j) = (fg_components(i)%priors(j,2)-fg_components(i)%priors(j,1)) * &
                  & (k-1) / real(numgrid-1,dp) + fg_components(i)%priors(j,1)
          end do
          if (.not. trim(fg_components(i)%type) == 'physical_dust') then
             fg_components(i)%priors(j,1) = fg_components(i)%priors(j,1) + 0.1d0*delta
             fg_components(i)%priors(j,2) = fg_components(i)%priors(j,2) - 0.1d0*delta
          end if
       end do

       ! Set up index regions
       allocate(fg_components(i)%indregs(fg_components(i)%npar))
       do j = 1, fg_components(i)%npar
          call int2string(j, j_text)
          paramname = 'REGION_DEFINITION' // i_text // '_' // j_text
          call get_parameter(paramfile, paramname, par_string=filename)
          call initialize_index_regions(chaindir, i, j, filename, mask, fg_components(i)%mask, nside, nmaps, &
               & fg_components(i)%fwhm_p(j), fg_components(i)%indregs(j), fg_components(i)%priors(j,1:2))
       end do
          
    end do

    call update_eff_fg_spectrum

    deallocate(mask)

!    write(*,*) fg_components%output_freq_comp_maps
!    call mpi_finalize(ierr)
!    stop

  end subroutine init_comm_fg_component_mod

  subroutine cleanup_comm_fg_component_mod
    implicit none


  end subroutine cleanup_comm_fg_component_mod



  ! ***************************
  ! *                         * 
  ! *    Interface routines   *
  ! *                         *
  ! ***************************


  ! External interface
  function get_effective_fg_spectrum(fg_comp, band, fg_params, pixel, pol)
    implicit none

    integer(i4b),                     intent(in)            :: band
    integer(i4b),                     intent(in), optional  :: pixel, pol
    type(fg_meta_data),               intent(in)            :: fg_comp
    real(dp),           dimension(:), intent(in), optional  :: fg_params
    real(dp)                                                :: get_effective_fg_spectrum

    real(dp)     :: p_low, p_high, x, y, epsilon = 1d-6, S_eff, nu, dnu, T_d, beta, C, scale
    integer(i4b) :: i, j, k, n, num_samp_point
    integer(i4b),              dimension(2)   :: bin
    real(dp),     allocatable, dimension(:,:) :: spectrum

    if (present(pixel)) then
       if (fg_comp%mask(pixel,pol) < 0.5d0) then
          get_effective_fg_spectrum = 0.d0
          return
       end if
    end if

    if (.not. trim(fg_comp%type) == 'CO_multiline') then
       do i = 1, fg_comp%npar
          p_low  = fg_comp%priors(i,1)
          p_high = fg_comp%priors(i,2)
          if (fg_params(i) < p_low-epsilon .or. fg_params(i) > p_high+epsilon) then
             write(*,*) 'fg_component_mod: Error -- parameter out of bounds'
             write(*,*) 'fg_component_mod: myid_chain  = ', myid_chain
             write(*,*) 'fg_component_mod: component   = ', trim(fg_comp%label)
             write(*,*) 'fg_component_mod: band        = ', band
             write(*,*) 'fg_component_mod: param       = ', fg_params(i)
             write(*,*) 'fg_component_mod: lower prior = ', p_low
             write(*,*) 'fg_component_mod: upper prior = ', p_high
             if (present(pixel)) write(*,*) 'fg_component_mod: pixel        = ', pixel
             stop
          end if
       end do
    end if

    if (trim(fg_comp%type) == 'CO_multiline') then

       get_effective_fg_spectrum = compute_CO_multiline_spectrum(band, &
            & fg_comp%co_band(1:fg_comp%npar+1), fg_comp%nu_ref, fg_params) * ant2data(band)

    else if (fg_comp%npar == 0) then
          
       get_effective_fg_spectrum = fg_comp%S_tabulated(band) 
       
    else if (fg_comp%npar == 1) then

       get_effective_fg_spectrum = splint(fg_comp%par(:,1), fg_comp%S_1D(:,1,band), fg_comp%S_1D(:,2,band), &
            & fg_params(1)) * fg_comp%S_tabulated(band) 
          
    else if (fg_comp%npar == 2) then
          
       get_effective_fg_spectrum = splin2_full_precomp(fg_comp%par(:,1), fg_comp%par(:,2), &
            & fg_comp%S_2D(:,:,:,:,band), fg_params(1), fg_params(2)) * fg_comp%S_tabulated(band) 

    else

       write(*,*) 'fg_component_mod: Error -- larger number of parameters than two not allowed'
       stop
          
    end if

    !get_effective_fg_spectrum = get_effective_fg_spectrum * spec2data(band, fg_comp%fg_id, outtype='uK_ant')
    get_effective_fg_spectrum = get_effective_fg_spectrum * bp(band)%gain 

  end function get_effective_fg_spectrum


  function get_effective_deriv_fg_spectrum(fg_comp, band, par_id, fg_params, pixel, pol)
    implicit none

    integer(i4b),                     intent(in)           :: band, par_id
    integer(i4b),                     intent(in), optional :: pixel, pol
    type(fg_meta_data),               intent(in)           :: fg_comp
    real(dp),           dimension(:), intent(in)           :: fg_params
    real(dp)                                               :: get_effective_deriv_fg_spectrum

    real(dp)     :: p_low, p_high, x, y, epsilon = 1d-6, delta = 1.d-10, S_eff_1, S_eff_2, nu, dnu
    real(dp)     :: T_d, beta, C, p(2)
    integer(i4b) :: i, j, k, n, num_samp_point
    integer(i4b), dimension(2) :: bin
    real(dp),     allocatable, dimension(:,:) :: spectrum

    if (present(pixel)) then
       if (fg_comp%mask(pixel,pol) < 0.5d0) then
          get_effective_deriv_fg_spectrum = 0.d0
          return
       end if
    end if

    if (trim(fg_comp%type) == 'CO_multiline') then

       if (band == fg_comp%co_band(1)) then
          get_effective_deriv_fg_spectrum = 0.d0          
       else
          get_effective_deriv_fg_spectrum = compute_CO_multiline_spectrum(band, &
               & fg_comp%co_band(1:fg_comp%npar+1), fg_comp%nu_ref, fg_params, deriv=.true.) 
       end if

    else if (fg_comp%npar == 0) then
          
       get_effective_deriv_fg_spectrum = 0.d0
       
    else if (fg_comp%npar == 1) then

       get_effective_deriv_fg_spectrum = splint_deriv(fg_comp%par(:,1), fg_comp%S_1D(:,1,band), &
            & fg_comp%S_1D(:,2,band), fg_params(1)) * fg_comp%S_tabulated(band) 
          
    else if (fg_comp%npar == 2) then

       p         = fg_params
       p(par_id) = p(par_id) + delta
       S_eff_1 = get_effective_fg_spectrum(fg_comp, band, fg_params)
       S_eff_2 = get_effective_fg_spectrum(fg_comp, band, p)
       get_effective_deriv_fg_spectrum = (S_eff_2-S_eff_1) / delta * fg_comp%S_tabulated(band) 

    else

       write(*,*) 'fg_component_mod: Error -- larger number of parameters than two not allowed'
       stop
          
    end if

    !get_effective_deriv_fg_spectrum = get_effective_deriv_fg_spectrum * spec2data(band, fg_comp%fg_id, &
    !     & outtype='uK_ant')
    get_effective_deriv_fg_spectrum = get_effective_deriv_fg_spectrum * bp(band)%gain

  end function get_effective_deriv_fg_spectrum


  function spec2data(band, comp, cmb, intype, outtype, gain)
    implicit none

    integer(i4b),     intent(in)           :: band
    integer(i4b),     intent(in), optional :: comp
    logical(lgt),     intent(in), optional :: cmb
    character(len=*), intent(in), optional :: intype, outtype
    real(dp),         intent(in), optional :: gain
    real(dp)                               :: spec2data

    character(len=128) :: intype_, outtype_

    outtype_ = bp(band)%unit; if (present(outtype)) outtype_ = outtype
    if (present(cmb)) then
       intype_ = 'uK_cmb'
    else if (present(intype)) then
       intype_ = intype
    else 
       intype_ = fg_components(comp)%amp_unit
    end if

    ! Convert from model units to antenna units
    if (trim(intype_) == 'uK_cmb') then
       if (present(cmb)) then
          spec2data = 1.d0 / bp(band)%a2t 
       else
          spec2data = 1.d0 / bp(fg_components(comp)%ref_band)%a2t 
       end if
    else if (trim(intype_) == 'K km/s') then
       spec2data = bp(fg_components(comp)%ref_band)%co2t / bp(fg_components(comp)%ref_band)%a2t
    else if (trim(intype_) == 'MJy/sr') then
       spec2data = bp(fg_components(comp)%ref_band)%f2t / bp(fg_components(comp)%ref_band)%a2t
    else if (trim(intype_) == 'y_SZ') then
       spec2data = 1.d0 / bp(fg_components(comp)%ref_band)%a2sz
    else if (trim(intype_) == 'uK_ant') then
       spec2data = 1.d0
    else
       write(*,*) 'Unsupported unit = ', trim(intype_)
       stop
    end if

    ! Convert from antenna units to map units
    if (trim(outtype_) == 'uK_cmb') then
       spec2data = spec2data * bp(band)%a2t
    else if (trim(outtype_) == 'MJy/sr') then
       spec2data = spec2data * bp(band)%a2t / bp(band)%f2t
    else if (trim(outtype_) == 'K km/s') then
       spec2data = spec2data * bp(band)%a2t / bp(band)%co2t
    else if (trim(outtype_) == 'y_SZ') then
       spec2data = spec2data * bp(band)%a2sz
    else if (trim(outtype_) == 'uK_ant') then
       spec2data = spec2data 
    else 
       write(*,*) 'Unsupported data unit type = ', trim(outtype_)
       stop
    end if

    ! Apply gain
    if (present(gain)) then
       spec2data = spec2data * gain
    else
       spec2data = spec2data * bp(band)%gain
    end if

  end function spec2data

  subroutine initialize_index_map(paramfile, rng_handle, fg_param_map)
    implicit none
    
    character(len=128),                      intent(in)    :: paramfile
    type(planck_rng),                        intent(inout) :: rng_handle
    real(dp),           dimension(0:,1:,1:), intent(out)   :: fg_param_map
    
    real(dp)           :: index
    integer(i4b)       :: i, j, k, l, m, n, pol_ind, counter, npix, nmaps, p, s, ierr
    character(len=2)   :: i_text, j_text
    character(len=128) :: paramtext, filename
    real(dp), allocatable, dimension(:,:) :: map

    npix  = size(fg_param_map,1)
    nmaps = size(fg_param_map,2)

    counter = 1
    do k = 1, num_fg_comp

       if (fg_components(k)%npar == 0) cycle

       do m = 1, fg_components(k)%npar

          if (trim(fg_components(k)%init_mode) == 'input_map') then
             allocate(map(0:npix-1,nmaps))
             call int2string(k,i_text)
             call int2string(m,j_text)
             call get_parameter(paramfile, 'INIT_INDEX_MAP' // i_text // '_' // j_text, par_string=filename)
             call read_map(filename, map)
             if (trim(fg_components(k)%type) == 'freefree_EM' .and. m ==1) then
                map = log(map)
             end if
             map = min(max(map,fg_components(k)%priors(m,1)), fg_components(k)%priors(m,2))
          end if

          do i = 1, fg_components(k)%indregs(m)%nreg
             if (trim(fg_components(k)%init_mode) == 'input_map') then
                index = 0.d0; n = 0
                do j = 1, fg_components(k)%indregs(m)%regions(i)%n
                   p = fg_components(k)%indregs(m)%regions(i)%pix(j,1)
                   s = fg_components(k)%indregs(m)%regions(i)%pix(j,2)
                   index = index + map(p,s)
                   n = n+1
                end do
                index = index / n
             else if (trim(fg_components(k)%init_mode) == 'prior') then
                index = fg_components(k)%gauss_prior(m,1) + &
                     & fg_components(k)%gauss_prior(m,2) * rand_gauss(rng_handle)
                index = max(min(index, fg_components(k)%priors(m,2)), fg_components(k)%priors(m,1))
             else if (trim(fg_components(k)%init_mode) == 'default') then
                index = fg_components(k)%priors(m,3)
                if (trim(fg_components(k)%type) == 'freefree_EM' .and. m ==1) then
                   index = log(index)
                end if
             else
                write(*,*) 'Unknown initialization mode: ', trim(fg_components(k)%init_mode)
                stop
             end if
             do j = 1, fg_components(k)%indregs(m)%regions(i)%n
                p = fg_components(k)%indregs(m)%regions(i)%pix(j,1)
                s = fg_components(k)%indregs(m)%regions(i)%pix(j,2)
                fg_param_map(p,s,counter) = index
             end do
          end do

!          if (fg_components(k)%fwhm_p(m) > 0.d0) then
!             call smooth_par_map(fg_param_map(:,:,counter), fg_components(k)%fwhm_p(m), &
!                  & fg_components(k)%priors(m,1:2))             
!          end if

          if (allocated(map)) deallocate(map)
          counter = counter + 1             

       end do
    end do

  end subroutine initialize_index_map


  subroutine reorder_fg_params(p_in, p_out)
    implicit none

    real(dp),       dimension(:,:), intent(in)  :: p_in
    type(fg_params)                             :: p_out

    integer(i4b) :: i, j, k, n, counter

    if (.not. allocated(p_out%comp)) then
       allocate(p_out%comp(num_fg_comp))
       do i = 1, num_fg_comp
          n = fg_components(i)%npar       
          if (n > 0) then
             allocate(p_out%comp(i)%p(nmaps,n))
          else
             allocate(p_out%comp(i)%p(nmaps,1))
             p_out%comp(i)%p(:,1) = 0.d0
          end if
       end do
    end if

    do k = 1, nmaps
       counter = 1
       do i = 1, num_fg_comp
          n = fg_components(i)%npar
          if (n > 0) then 
             do j = 1, n
                p_out%comp(i)%p(k,j) = p_in(k,counter)
                counter              = counter + 1
             end do
          end if
       end do
    end do

  end subroutine reorder_fg_params
  
  subroutine deallocate_fg_params(p, p_1D, p_2D)
    implicit none

    type(fg_params),                 optional :: p
    type(fg_params), dimension(:),   optional :: p_1D
    type(fg_params), dimension(:,:), optional :: p_2D

    integer(i4b) :: i, j, k, l

    if (present(p)) then
       do i = 1, size(p%comp)
          if (allocated(p%comp(i)%p)) deallocate(p%comp(i)%p)
       end do
       if (allocated(p%comp)) deallocate(p%comp)
    end if

    if (present(p_1D)) then
       do j = 1, size(p_1D)
          do i = 1, size(p_1D(j)%comp)
             if (allocated(p_1D(j)%comp(i)%p)) deallocate(p_1D(j)%comp(i)%p)
          end do
          if (allocated(p_1D(j)%comp)) deallocate(p_1D(j)%comp)
       end do
    end if

    if (present(p_2D)) then
       do k = 1, size(p_2D,2)
          do j = 1, size(p_1D,1)
             do i = 1, size(p_2D(j,k)%comp)
                if (allocated(p_2D(j,k)%comp(i)%p)) deallocate(p_2D(j,k)%comp(i)%p)
             end do
             if (allocated(p_2D(j,k)%comp)) deallocate(p_2D(j,k)%comp)
          end do
       end do
    end if

  end subroutine deallocate_fg_params


  ! ********************************************************************
  ! *                                                                  *
  ! *               Routines for computing spectra                     *
  ! *                                                                  *
  ! ********************************************************************

  function get_ideal_fg_spectrum(fg_comp, p, nu)
    implicit none

    type(fg_meta_data),               intent(in) :: fg_comp
    real(dp),           dimension(:), intent(in) :: p
    real(dp),                         intent(in) :: nu
    real(dp)                                     :: get_ideal_fg_spectrum

    if (trim(fg_comp%type) == 'CO_multiline') then
       !get_ideal_spectrum = compute_CO_multiline_spectrum(band, &
       !     & fg_comp%co_band(1:fg_comp%npar+1), fg_comp%nu_ref, fg_params)
       get_ideal_fg_spectrum = 0.d0
    else if (trim(fg_comp%type) == 'cmb') then
       get_ideal_fg_spectrum = 1.d0 / compute_ant2thermo_single(nu)
    else if (trim(fg_comp%type) == 'sz') then
       get_ideal_fg_spectrum = sz_thermo_single(nu) / compute_ant2thermo_single(nu)
    else if (trim(fg_comp%type) == 'freefree') then
       get_ideal_fg_spectrum = compute_freefree_spectrum(nu, fg_comp%nu_ref, p(1)) 
    else if (trim(fg_comp%type) == 'AME_freq_shift') then                         
       get_ideal_fg_spectrum = compute_AME_freq_shift_spectrum(nu, fg_comp%nu_ref, p(1), &
            & fg_comp%nu_peak, fg_comp%S_nu_ref, fg_comp%p_rms)
    else if (trim(fg_comp%type) == 'AME_freq_shift_2par') then                         
       get_ideal_fg_spectrum = compute_AME_freq_shift_2par_spectrum(nu, fg_comp%nu_ref, p(1), p(2), &
            & fg_comp%nu_peak, fg_comp%S_nu_ref, fg_comp%p_rms)
    else if (trim(fg_comp%type) == 'AME_lognormal_width') then
       get_ideal_fg_spectrum = compute_AME_lognormal_width_spectrum(nu, fg_comp%nu_ref, p(1), p(2), fg_comp%p_rms)
    else if (trim(fg_comp%type) == 'magnetic_dust') then                         
       get_ideal_fg_spectrum = compute_magnetic_dust_spectrum(nu, fg_comp%nu_ref, &
            & p(1),fg_comp%S_nu_ref, fg_comp%p_rms)
    else if (trim(fg_comp%type) == 'one-component_dust') then
       get_ideal_fg_spectrum = compute_one_comp_dust_spectrum(nu, fg_comp%nu_ref, &
            & p(1), p(2), fg_comp%nu_flat, fg_comp%frac_flat, fg_comp%p_rms)
    else if (trim(fg_comp%type) == 'physical_dust') then                         
       get_ideal_fg_spectrum = compute_physical_dust_spectrum(nu, p(1), fg_comp%us, & 
            fg_comp%nu_ref, fg_comp%S_phys_dust, fg_comp%S_dust_coeff, fg_comp%p_rms)
    else if (trim(fg_comp%type) == 'power_law') then
       get_ideal_fg_spectrum = compute_power_law_spectrum(nu, fg_comp%nu_ref, &
            & p(1), p(2), fg_comp%p_rms)
    else if (trim(fg_comp%type) == 'power_law_break') then
       get_ideal_fg_spectrum = compute_power_law_break_spectrum(nu, fg_comp%nu_ref, &
            & p(2), p(1), fg_comp%nu_break, fg_comp%p_rms) 
    else if (trim(fg_comp%type) == 'freefree_EM') then
       get_ideal_fg_spectrum = compute_freefree_EM_spectrum(nu, p(1), p(2)) 
    end if

  end function get_ideal_fg_spectrum

  function compute_power_law_spectrum(nu, nu_ref, beta, C, rms)
    implicit none

    real(dp),               intent(in) :: nu, nu_ref, beta, C
    real(dp), dimension(:), intent(in) :: rms
    real(dp)                           :: compute_power_law_spectrum

    integer(i4b) :: n_beta, n_C, i, j
    real(dp)     :: w, w_tot, betap, Cp, beta_min, beta_max, dbeta, C_min, C_max, dC, beta_rms, C_rms

    n_beta = 1; beta_min = beta; beta_max = beta; dbeta = 0.d0; beta_rms = 1.d0
    n_C    = 1; C_min    = C;    C_max    = C;    dC    = 0.d0; C_rms    = 1.d0 
    if (rms(1) > 0.d0) then
       n_beta   = N_RMS_MARG
       beta_min = beta - 3.d0*rms(1)
       beta_max = beta + 3.d0*rms(1)
       dbeta    = (beta_max - beta_min) / (n_beta-1)
       beta_rms = rms(1)
    end if
    if (rms(2) > 0.d0) then
       n_C     = N_RMS_MARG
       C_min   = C - 3.d0*rms(2)
       C_max   = C + 3.d0*rms(2)
       dC      = (C_max - C_min) / (n_C-1)
       C_rms   = rms(2)
    end if

    compute_power_law_spectrum = 0.d0
    w_tot                      = 0.d0
    do i = 1, n_beta
       betap = beta_min + (i-1)*dbeta
       do j = 1, n_C
          Cp                         = C_min + (j-1)*dC
          w                          = exp(-0.5d0*(((betap-beta)/beta_rms)**2 + ((Cp-C)/C_rms)**2))
          w_tot                      = w_tot + w
          compute_power_law_spectrum = compute_power_law_spectrum + w*(nu/nu_ref)**(betap + Cp*log(nu/nu_ref))
       end do
    end do
    compute_power_law_spectrum = compute_power_law_spectrum / w_tot

!    compute_power_law_spectrum = (nu/nu_ref)**(beta + C*log(nu/nu_ref))
!    write(*,fmt='(2e16.8,2f8.3,e16.8)') nu, nu_ref, beta, C, compute_power_law_spectrum


  end function compute_power_law_spectrum

  function compute_power_law_break_spectrum(nu, nu_ref, dbeta, beta, nu_break, rms)
    implicit none

    real(dp),               intent(in) :: nu, nu_ref, dbeta, beta, nu_break
    real(dp), dimension(:), intent(in) :: rms
    real(dp)                           :: compute_power_law_break_spectrum
   
    real(dp) :: S, nu_b

    nu_b = 1.d9 * nu_break

    if (nu > nu_b) then
       S = (nu/nu_b)**beta
    else
       S = (nu/nu_b)**(beta+dbeta)
    end if

    if (nu_ref > nu_b) then
       S = S * (nu_b/nu_ref)**beta
    else
       S = S * (nu_b/nu_ref)**(beta+dbeta)
    end if

    compute_power_law_break_spectrum = S

  end function compute_power_law_break_spectrum

  function compute_freefree_spectrum(nu, nu_ref, T_e)
    implicit none

    real(dp),               intent(in) :: nu, nu_ref, T_e
    real(dp)                           :: compute_freefree_spectrum

    real(dp)     :: S, S_ref

    S     = log(exp(5.960d0 - sqrt(3.d0)/pi * log(1.d0 * nu    /1.d9 * (T_e/1.d4)**(-1.5d0))) + 2.71828d0)
    S_ref = log(exp(5.960d0 - sqrt(3.d0)/pi * log(1.d0 * nu_ref/1.d9 * (T_e/1.d4)**(-1.5d0))) + 2.71828d0)
    compute_freefree_spectrum = S/S_ref * exp(-h*(nu-nu_ref)/k_b/T_e) * (nu/nu_ref)**(-2)

  end function compute_freefree_spectrum


  function compute_freefree_EM_spectrum(nu, log_EM, T_e)
    implicit none

    real(dp),               intent(in) :: nu, log_EM, T_e
    real(dp)                           :: compute_freefree_EM_spectrum

    real(dp)     :: g, Z_i, tau, EM

!    if (nu > 400d9) then
!       compute_freefree_EM_spectrum = 0.d0
!       return
!    end if

!    compute_freefree_EM_spectrum    = 0.83d0 * exp(log_EM) * &
!         & (1.d0 + 0.23d0*log(T_e/8000.d0) - 0.15d0 * log(nu/53.d9)) / (nu/53d9)**2 / sqrt(T_e/8000.d0)

!    return


    EM   = exp(log_EM)
    Z_i  = 1.d0
    g    = log(exp(5.960d0 - sqrt(3.d0)/pi * log(Z_i * nu/1.d9 * (T_e/1.d4)**(-1.5d0))) + 2.71828d0)
    tau  = 5.468d-2 * T_e**(-1.5d0) * (nu/1.d9)**(-2) * EM * g

    compute_freefree_EM_spectrum = 1.d6 * T_e * (1.d0 - exp(-tau)) ! Output signal in uK

  end function compute_freefree_EM_spectrum


  function compute_CO_multiline_spectrum(band, co_band, nu_ref, ratios, deriv)
    implicit none

    integer(i4b),               intent(in)  :: band
    integer(i4b), dimension(:), intent(in)  :: co_band
    real(dp),                   intent(in)  :: nu_ref
    real(dp),     dimension(:), intent(in)  :: ratios
    real(dp)                                :: compute_CO_multiline_spectrum
    logical(lgt),               intent(in), optional :: deriv

    integer(i4b) :: i

    compute_CO_multiline_spectrum = 0.d0
    if (band == co_band(1)) then
       compute_CO_multiline_spectrum = 1.d0
    else
       do i = 1, size(ratios)
          if (band == co_band(i+1)) then
             compute_CO_multiline_spectrum = ratios(i) * &
                  & (bp(band)%co2t/bp(band)%a2t) / (bp(co_band(1))%co2t/bp(co_band(1))%a2t)
             if (present(deriv)) then
                compute_CO_multiline_spectrum = &
                     & (bp(band)%co2t/bp(band)%a2t) / (bp(co_band(1))%co2t/bp(co_band(1))%a2t)
             end if
          end if
       end do
    end if

  end function compute_CO_multiline_spectrum

  function compute_AME_freq_shift_spectrum(nu, nu_ref, nu_p, nu_p0, S, rms)
    implicit none

    real(dp),                 intent(in)  :: nu, nu_ref, nu_p, nu_p0
    real(dp), dimension(:),   intent(in)  :: rms
    real(dp), dimension(:,:), intent(in)  :: S
    real(dp)                              :: compute_AME_freq_shift_spectrum

    integer(i4b) :: i, n_nu
    real(dp)     :: scale
    real(dp)     :: w, w_tot, nup, nu_min, nu_max, dnu, nu_rms, nu_low, nu_high

    nu_low = minval(S(:,1)); nu_high = maxval(S(:,1))
    if (nu < nu_low .or. nu > nu_high) then
       compute_AME_freq_shift_spectrum = 0.d0
    else

       n_nu = 1; nu_min = nu_p; nu_max = nu_p; dnu = 0.d0; nu_rms = 1.d0
       if (rms(1) > 0.d0) then
          n_nu   = N_RMS_MARG
          nu_min = max(nu_p - 3.d0*rms(1), 10.d0)
          nu_max = nu_p + 3.d0*rms(1)
          dnu    = (nu_max - nu_min) / (n_nu-1)
          nu_rms = rms(1)
       end if
       
       compute_AME_freq_shift_spectrum = 0.d0
       w_tot                           = 0.d0
       do i = 1, n_nu
          nup   = nu_min + (i-1)*dnu
          w     = exp(-0.5d0*((nup-nu_p)/nu_rms)**2)
          scale = nu_p0 / (nup*1.d9) ! nu_p is in GHz
          if (scale*nu > nu_low .and. scale*nu < nu_high) then
             w_tot = w_tot + w
             compute_AME_freq_shift_spectrum = compute_AME_freq_shift_spectrum + &
                  & w * exp(splint(S(:,1), S(:,2), S(:,3), scale*nu)) / &
                  & exp(splint(S(:,1), S(:,2), S(:,3), scale*nu_ref)) * (nu_ref/nu)**2
          end if
       end do
       if (w_tot > 0.d0) then
          compute_AME_freq_shift_spectrum = compute_AME_freq_shift_spectrum / w_tot
       else
          compute_AME_freq_shift_spectrum = 0.d0
       end if

       if (compute_AME_freq_shift_spectrum < 0.d0) then
          write(*,*) real(nu,sp), real(nu_ref,sp), real(nu_p,sp), real(rms,sp)
          write(*,*) real(compute_AME_freq_shift_spectrum,sp)
          stop
       end if

    end if

  end function compute_AME_freq_shift_spectrum

  function compute_AME_freq_shift_2par_spectrum(nu, nu_ref, nu_p, alpha, nu_p0, S, rms)
    implicit none

    real(dp),                 intent(in)  :: nu, nu_ref, nu_p, nu_p0, alpha
    real(dp), dimension(:),   intent(in)  :: rms
    real(dp), dimension(:,:), intent(in)  :: S
    real(dp)                              :: compute_AME_freq_shift_2par_spectrum

    integer(i4b) :: i, n_nu
    real(dp)     :: scale
    real(dp)     :: w, w_tot, nup, nu_min, nu_max, dnu, nu_rms, nu_low, nu_high

    nu_low = minval(S(:,1)); nu_high = maxval(S(:,1))
    if (nu < nu_low .or. nu > nu_high) then
       compute_AME_freq_shift_2par_spectrum = 0.d0
    else
       scale = nu_p0 / (nu_p*1.d9) ! nu_p is in GHz
       if (scale*nu > nu_low .and. scale*nu < nu_high) then
          compute_AME_freq_shift_2par_spectrum = &
               & exp(splint(S(:,1), S(:,2), S(:,3), scale*nu)) / &
               & exp(splint(S(:,1), S(:,2), S(:,3), scale*nu_ref)) * (nu_ref/nu)**(2.d0-alpha)
       else
          compute_AME_freq_shift_2par_spectrum = 0.d0
       end if
       return
    end if

  end function compute_AME_freq_shift_2par_spectrum

  function compute_AME_lognormal_width_spectrum(nu, nu_ref, nu_p, width, rms)
    implicit none
    real(dp),               intent(in) :: nu, nu_ref, nu_p, width
    real(dp), dimension(:), intent(in) :: rms
    real(dp)                           :: peak, compute_AME_lognormal_width_spectrum

    ! Moving from flux density to T_RJ means an multiplicative (1/nu^2) to the spectrum
    peak = nu_p*1d9
    ! write(*,*) peak, width
    compute_AME_lognormal_width_spectrum = exp(-0.5*(log(nu/peak)/width)**2)*(nu_ref/nu)**2

  end function compute_AME_lognormal_width_spectrum

  function compute_magnetic_dust_spectrum(nu, nu_ref, T_m, S, rms)
    implicit none

    real(dp),                 intent(in)  :: nu, nu_ref, T_m
    real(dp), dimension(:),   intent(in)  :: rms
    real(dp), dimension(:,:), intent(in)  :: S
    real(dp)                              :: compute_magnetic_dust_spectrum

    integer(i4b) :: i, n_T
    real(dp)     :: scale
    real(dp)     :: w, w_tot, T, T_min, T_max, dT, T_rms, nu_low, nu_high

    nu_low = minval(S(:,1)); nu_high = maxval(S(:,1))
    if (nu < nu_low .or. nu > nu_high) then
       compute_magnetic_dust_spectrum = 0.d0
    else

       n_T = 1; T_min = T_m; T_max = T_m; dT = 0.d0; T_rms = 1.d0
       if (rms(1) > 0.d0) then
          n_T    = N_RMS_MARG
          T_min  = T_m - 3.d0*rms(1)
          T_max  = T_m + 3.d0*rms(1)
          dT     = (T_max - T_min) / (n_T-1)
          T_rms = rms(1)
       end if
       
       compute_magnetic_dust_spectrum = 0.d0
       w_tot                          = 0.d0
       do i = 1, n_T
          T     = T_min + (i-1)*dT
          w     = exp(-0.5d0*((T-T_m)/T_rms)**2)
          w_tot = w_tot + w
          compute_magnetic_dust_spectrum = compute_magnetic_dust_spectrum + &
               & w * splint(S(:,1), S(:,2), S(:,3), nu) / &
               &     splint(S(:,1), S(:,2), S(:,3), nu_ref) * &
               & (exp(h*nu_ref/(k_b*T))-1.d0) / (exp(h*nu/(k_b*T))-1.d0) * (nu/nu_ref)**4 * & 
               & (nu_ref/nu)**2
       end do
       compute_magnetic_dust_spectrum = compute_magnetic_dust_spectrum / w_tot

       if (compute_magnetic_dust_spectrum < 0.d0) then
          write(*,*) 'Error in magnetic dust spectrum:'
          write(*,*) real(nu,sp), real(nu_ref,sp), real(T_m,sp), real(rms,sp)
          write(*,*) real(compute_magnetic_dust_spectrum,sp)
          stop
       end if

    end if

  end function compute_magnetic_dust_spectrum

  function compute_physical_dust_spectrum(nu, u, u_array, nu_ref, S, S_coeff, rms)
    implicit none

    real(dp),                 intent(in)      :: nu, u, nu_ref
    real(dp), dimension(:),   intent(in)      :: rms, u_array
    real(dp), dimension(:,:), intent(in)      :: S
    real(dp), dimension(:,:,:,:), intent(in)  :: S_coeff
    real(dp)                                  :: compute_physical_dust_spectrum

    real(dp)    :: nu_low, nu_high

    nu_low = minval(S(1,:)); nu_high = maxval(S(1,:))
    if (nu < nu_low .or. nu > nu_high) then
       compute_physical_dust_spectrum = 0.d0
    else
       compute_physical_dust_spectrum = exp(splin2_full_precomp(u_array, log(S(1,:)), S_coeff, u, log(nu)))
       if (myid_chain == root) then
          if (compute_physical_dust_spectrum < 0.d0) then
             write(*,*) compute_physical_dust_spectrum, u, nu
          end if
       end if
    end if

  end function compute_physical_dust_spectrum


  function compute_one_comp_dust_spectrum(nu, nu_ref, beta, T_d, nu_flat, frac_flat, rms)
    implicit none

    real(dp),               intent(in) :: nu, nu_ref, beta, T_d, nu_flat, frac_flat
    real(dp), dimension(:), intent(in) :: rms
    real(dp)                           :: compute_one_comp_dust_spectrum

    integer(i4b) :: n_beta, n_T, i, j
    real(dp)     :: w, w_tot, betap, Tp, beta_min, beta_max, dbeta, T_min, T_max, dT, beta_rms, T_rms, f, S, x
    real(dp)     :: alpha, gamma, delta, nu_steep
!!$    alpha = alpha_flat
!!$    gamma  = gamma_flat
    x      = h / (k_B*T_d)
!!$    if (nu >= nu_flat) then
    !alpha (fraq_flat) between 2-6
    
    compute_one_comp_dust_spectrum = (exp(x*nu_ref)-1.d0) / (exp(x*nu)-1.d0) * (nu/nu_ref)**(beta+1.d0)
    !compute_one_comp_dust_spectrum = (exp(x*nu_ref)-1.d0) / (exp(x*nu)-1.d0) * (nu/nu_ref)**(beta+1.d0)*(tanh(x*frac_flat*nu)/tanh(x*frac_flat*nu_ref))**(2.d0-beta)
    !compute_one_comp_dust_spectrum = (exp(x*nu_ref)-1.d0) / (exp(x*nu)-1.d0) * (nu/nu_ref)**(beta+1.d0)*(tanh(x*frac_flat*nu)/tanh(x*frac_flat*nu_ref))
    return

    nu_steep = nu_flat !100.d9
    delta    = -frac_flat !0.2d0 ! Index steepening
    if (nu > nu_steep) then
       compute_one_comp_dust_spectrum = (exp(x*nu_ref)-1.d0) / (exp(x*nu)-1.d0) * (nu/nu_ref)**(beta+1.d0)
    else
       compute_one_comp_dust_spectrum = (exp(x*nu_steep)-1.d0) / (exp(x*nu)-1.d0) * (nu/nu_steep)**(beta+1.d0+delta) * &
            (exp(x*nu_ref)-1.d0) / (exp(x*nu_steep)-1.d0) * (nu_steep/nu_ref)**(beta+1.d0)
    end if
!!$    else
!!$       compute_one_comp_dust_spectrum = &
!!$            & ((exp(x*nu_ref)-1.d0) / (exp(x*nu_flat)-1.d0) * (nu_flat/nu_ref)**(beta+1.d0)) * &
!!$            & ((exp(x*nu_flat)-1.d0) / (exp(x*nu)-1.d0) * (nu/nu_flat)**(beta+1.d0+alpha*T_d+gamma*T_d**2))
!!$    end if
!!$
    return

    n_beta = 1; beta_min = beta; beta_max = beta; dbeta = 0.d0; beta_rms = 1.d0
    n_T    = 1; T_min    = T_d;  T_max    = T_d;  dT    = 0.d0; T_rms    = 1.d0 
    if (rms(1) > 0.d0) then
       n_beta   = N_RMS_MARG
       beta_min = beta - 3.d0*rms(1)
       beta_max = beta + 3.d0*rms(1)
       dbeta    = (beta_max - beta_min) / (n_beta-1)
       beta_rms = rms(1)
    end if
    if (rms(2) > 0.d0) then
       n_T     = N_RMS_MARG
       T_min   = max(T_d - 3.d0*rms(2), 5.d0)
       T_max   = T_d + 3.d0*rms(2)
       dT      = (T_max - T_min) / (n_T-1)
       T_rms   = rms(2)
    end if

    compute_one_comp_dust_spectrum = 0.d0
    w_tot                          = 0.d0
    do i = 1, n_beta
       betap = beta_min + (i-1)*dbeta
       do j = 1, n_T
          Tp     = T_min + (j-1)*dT
          x      = h / (k_B*Tp)
          w      = exp(-0.5d0*(((betap-beta)/beta_rms)**2 + ((Tp-T_d)/T_rms)**2))
          w_tot  = w_tot + w
          compute_one_comp_dust_spectrum = compute_one_comp_dust_spectrum + &
               & w*((exp(x*nu_ref)-1.d0) / (exp(x*nu)-1.d0) * (nu / nu_ref) * &
                   & ((nu/nu_ref)**betap + frac_flat*(nu_flat/nu_ref)**betap))
       end do
    end do
    compute_one_comp_dust_spectrum = compute_one_comp_dust_spectrum / w_tot
    
    if (compute_one_comp_dust_spectrum < 0.d0) then
       write(*,*) 'Dust error = ', real(compute_one_comp_dust_spectrum,sp)
       stop
    end if

  end function compute_one_comp_dust_spectrum



  ! ***************************************************************************
  ! *                                                                         * 
  ! *    Routines for computing effective spectra by integration over bands   *
  ! *                                                                         *
  ! ***************************************************************************

!!$  subroutine compute_effective_spectrum(nu_a, nu_b, S_nu, S_eff, S_nu_full)
!!$    implicit none
!!$
!!$    real(dp),                   intent(in)           :: nu_a, nu_b
!!$    real(dp), dimension(1:,1:), intent(in)           :: S_nu
!!$    real(dp),                   intent(out)          :: S_eff
!!$    real(dp), dimension(1:,1:), intent(in), optional :: S_nu_full
!!$
!!$    integer(i4b) :: i, j, numpt
!!$    real(dp)     :: integral, S_ref, d, d_bf, nu0
!!$    real(dp), allocatable, dimension(:) :: y2, x, y
!!$
!!$    if (present(S_nu_full)) then
!!$
!!$       ! Do a nearest neighbor search
!!$       nu0  = 0.5d0 * (nu_a + nu_b)
!!$       d_bf = 1.d30
!!$       do i = 1, size(S_nu_full,1)
!!$          d = abs(nu0 - S_nu_full(i,1))
!!$          if (d < d_bf) then
!!$             S_eff = S_nu_full(i,2)
!!$             d_bf  = d
!!$          end if
!!$       end do
!!$
!!$    else if (nu_a == nu_b) then
!!$
!!$       ! Define spectrum at center frequency
!!$       S_eff = S_nu(1,2)
!!$
!!$    else
!!$
!!$       ! Compute spline
!!$       numpt = size(S_nu(:,1))
!!$       allocate(x(numpt))
!!$       allocate(y(numpt))
!!$       allocate(y2(numpt))
!!$       
!!$       x = S_nu(:,1)
!!$       y = S_nu(:,2)
!!$    
!!$       call spline(x, y, 1.d30, 1.d30, y2)
!!$
!!$       y  = y  / (nu_b - nu_a)
!!$       y2 = y2 / (nu_b - nu_a)
!!$
!!$       ! Integrate over frequency band
!!$       call qsimp(x, y, y2, nu_a, nu_b, S_eff)
!!$
!!$       deallocate(x)
!!$       deallocate(y)
!!$       deallocate(y2)
!!$       
!!$    end if
!!$
!!$  end subroutine compute_effective_spectrum

!!$  function get_nearest_neighbor(x, y, x0)
!!$    implicit none
!!$
!!$    real(dp),                intent(in) :: x0
!!$    real(dp), dimension(1:), intent(in) :: x, y
!!$    real(dp)                            :: get_nearest_neighbor
!!$
!!$    real(dp)     :: d, d_bf
!!$    integer(i4b) :: i, n
!!$
!!$    ! Do a nearest neighbor search
!!$    n = size(x)
!!$    d_bf = 1.d30*maxval(x)
!!$    do i = 1, n
!!$       d = abs(x0 - x(i))
!!$       if (d < d_bf) then
!!$          get_nearest_neighbor = y(i)
!!$          d_bf  = d
!!$       end if
!!$    end do
!!$
!!$  end function get_nearest_neighbor
!!$
!!$
!!$  subroutine compute_eff_power_law_spectrum(nu_a, nu_b, nu_ref, beta, S_eff)
!!$    implicit none
!!$
!!$    real(dp),     intent(in)  :: nu_a, nu_b, nu_ref, beta
!!$    real(dp),     intent(out) :: S_eff
!!$
!!$    real(dp) :: nu_eff
!!$
!!$    ! Analytic integration over bandpass by effective frequency
!!$    if (nu_b > nu_a) then
!!$       nu_eff = (1.d0/(beta+1.d0) * (nu_b**(beta+1.d0) - nu_a**(beta+1.d0)) / &
!!$            & (nu_b - nu_a))**(1.d0/beta)
!!$    else
!!$       nu_eff = nu_a
!!$    end if
!!$
!!$    S_eff = (nu_eff/nu_ref)**beta
!!$
!!$  end subroutine compute_eff_power_law_spectrum


  subroutine update_eff_fg_spectrum(comp)
    implicit none
    
    integer(i4b), intent(in), optional :: comp

    integer(i4b) :: i, j, k, l, m, n, num_samp_point, ierr, c
    real(dp)     :: nu, dnu, nu_eff
    real(dp), allocatable, dimension(:)   :: s
    real(dp), allocatable, dimension(:,:,:) :: my_grid

    if (myid_chain == root) then
       c = -1; if (present(comp)) c = comp
    end if
    call mpi_bcast(c, 1,   MPI_INTEGER, root, comm_chain, ierr)

    ! Initialize the effective spectrum grids
    num_samp_point = 100
    do i = 1, num_fg_comp

       if (c > 0 .and. c /= i) cycle

       if (trim(fg_components(i)%type) == 'CO_multiline') then

          ! Do nothing

       else if (fg_components(i)%npar == 1) then

          do l = 1, numband
             n = bp(l)%n
             allocate(s(n))
             do k = 1, numgrid
                do j = 1, n
                   nu = bp(l)%nu(j)
                   if (trim(fg_components(i)%type) == 'freefree') then
                      s(j) = compute_freefree_spectrum(nu, fg_components(i)%nu_ref, &
                           & fg_components(i)%par(k,1)) 
                   else if (trim(fg_components(i)%type) == 'AME_freq_shift') then                         
                      s(j) = compute_AME_freq_shift_spectrum(nu, fg_components(i)%nu_ref, &
                           & fg_components(i)%par(k,1), fg_components(i)%nu_peak, fg_components(i)%S_nu_ref, &
                           & fg_components(i)%p_rms)
                   else if (trim(fg_components(i)%type) == 'magnetic_dust') then                         
                      s(j) = compute_magnetic_dust_spectrum(nu, fg_components(i)%nu_ref, &
                           & fg_components(i)%par(k,1),fg_components(i)%S_nu_ref, &
                           & fg_components(i)%p_rms)
                   else if (trim(fg_components(i)%type) == 'physical_dust') then
                      s(j) = compute_physical_dust_spectrum(nu, fg_components(i)%par(k,1), &
                           & fg_components(i)%us, fg_components(i)%nu_ref, fg_components(i)%S_phys_dust, &
                           & fg_components(i)%S_dust_coeff, fg_components(i)%p_rms)
                   end if
                end do
                ! there was a spurious STOP here
                fg_components(i)%S_1D(k,1,l) = get_bp_avg_spectrum(l, s)
             end do
             deallocate(s)
          end do
          
          do l = 1, numband
             call spline(fg_components(i)%par(:,1), fg_components(i)%S_1D(:,1,l), 1.d30, 1.d30, &
                  & fg_components(i)%S_1D(:,2,l))
          end do
          
       else if (fg_components(i)%npar == 2) then
          
          allocate(my_grid(numgrid,numgrid,numband))
          my_grid = 0.d0
          do m = 1, numband
             n = bp(m)%n
             allocate(s(n))
             do k = 1+myid_chain, numgrid, numprocs_chain
                do l = 1, numgrid
                   do j = 1, n
                      nu = bp(m)%nu(j)
                      
                      if (trim(fg_components(i)%type) == 'one-component_dust') then
                         s(j) = compute_one_comp_dust_spectrum(nu, fg_components(i)%nu_ref, &
                              & fg_components(i)%par(k,1), fg_components(i)%par(l,2), &
                              & fg_components(i)%nu_flat, fg_components(i)%frac_flat, &
                              & fg_components(i)%p_rms)
                      else if (trim(fg_components(i)%type) == 'AME_lognormal_width') then
                         s(j) = compute_AME_lognormal_width_spectrum(nu, fg_components(i)%nu_ref, &
                              & fg_components(i)%par(k,1), fg_components(i)%par(l,2),&
                              & fg_components(i)%p_rms)
                      else if (trim(fg_components(i)%type) == 'power_law') then
                         s(j) = compute_power_law_spectrum(nu, fg_components(i)%nu_ref, &
                              & fg_components(i)%par(k,1), fg_components(i)%par(l,2), &
                              & fg_components(i)%p_rms)
                      else if (trim(fg_components(i)%type) == 'AME_freq_shift_2par') then                         
                         s(j) = compute_AME_freq_shift_2par_spectrum(nu, fg_components(i)%nu_ref, &
                              & fg_components(i)%par(k,1), fg_components(i)%par(l,2), &
                              & fg_components(i)%nu_peak, fg_components(i)%S_nu_ref, &
                              & fg_components(i)%p_rms)
                      else if (trim(fg_components(i)%type) == 'power_law_break') then
                         s(j) = compute_power_law_break_spectrum(nu, fg_components(i)%nu_ref, &
                              & fg_components(i)%par(l,2), fg_components(i)%par(k,1), &
                              & fg_components(i)%nu_break, fg_components(i)%p_rms) 
                      else if (trim(fg_components(i)%type) == 'freefree_EM') then
                         s(j) = compute_freefree_EM_spectrum(nu, &
                              & fg_components(i)%par(k,1), fg_components(i)%par(l,2)) 
                      end if
                   end do
                   my_grid(k,l,m) = get_bp_avg_spectrum(m, s)
                end do
             end do
             deallocate(s)
          end do
          call mpi_allreduce(MPI_IN_PLACE, my_grid, size(my_grid), MPI_DOUBLE_PRECISION, MPI_SUM, &
               & comm_chain, ierr)

          do m = 1, numband
             call splie2_full_precomp(fg_components(i)%par(:,1), fg_components(i)%par(:,2), my_grid(:,:,m), &
                  & fg_components(i)%S_2D(:,:,:,:,m))
          end do
          deallocate(my_grid)
          
       end if
       
    end do

  end subroutine update_eff_fg_spectrum


  subroutine update_fg_pix_response_map(band, pixels, fg_pix_spec_response, fg_param_map_in)
    implicit none

    integer(i4b),                      intent(in)           :: band
    integer(i4b), dimension(0:),       intent(in)           :: pixels
    real(dp),     dimension(0:,1:,1:), intent(out)          :: fg_pix_spec_response
    real(dp),     dimension(0:,1:,1:), intent(in), optional :: fg_param_map_in
 
    integer(i4b) :: i, j, k, p, numpix_highres, pix_nest, pix_ring, ierr
    real(dp), allocatable, dimension(:,:,:) :: fg_param_map
    type(fg_params)  :: fg_par

    allocate(fg_param_map(0:npix-1, nmaps, num_fg_par))
    if (myid_chain == root) then
       fg_param_map = fg_param_map_in
       call get_smooth_par_map(fg_param_map_in, fg_param_map)
    end if
    call mpi_bcast(fg_param_map, size(fg_param_map), MPI_DOUBLE_PRECISION, root, comm_chain, ierr)

    do i = 0, size(pixels)-1
       call reorder_fg_params(fg_param_map(pixels(i),:,:), fg_par)
       do k = 1, num_fg_comp
          do j = 1, nmaps
             if ((j == 1 .and. .not. sample_T_modes)) then
                fg_pix_spec_response(i,j,k) = 0.d0
             else
                fg_pix_spec_response(i,j,k) = &
                     & get_effective_fg_spectrum(fg_components(k), band, fg_par%comp(k)%p(j,:), &
                     & pixel=i, pol=j)
             end if
          end do
       end do
    end do

    call deallocate_fg_params(fg_par)
    deallocate(fg_param_map)

  end subroutine update_fg_pix_response_map



  ! *********************************************
  ! *                Utilities                  *
  ! *********************************************


  ! Routine for reading a spectrum file
  subroutine read_spectrum(filename, spectrum, nu_ref, sample_S_tab)
    implicit none

    character(len=128),                         intent(in)  :: filename
    real(dp),           pointer, dimension(:,:)             :: spectrum
    real(dp),                                   intent(in), optional  :: nu_ref
    logical(lgt),       pointer, dimension(:),  optional    :: sample_S_tab

    real(dp)            :: S_ref
    integer(i4b)        :: i, j, numpoint, unit, first, last
    character(len=128)  :: nu, val, string
    real(dp), allocatable, dimension(:) :: x, y, y2
    real(dp), allocatable, dimension(:,:) :: S_temp

    unit = getlun()
    open(unit, file=trim(filename))

    ! Find the number of entries
    numpoint = 0
    do while (.true.)
       read(unit,*,end=1) string

       if (string(1:1)=='#') cycle

       backspace(unit)
       read(unit,*,end=1) nu, val
       numpoint = numpoint + 1
    end do
1   close(unit)

    if (numpoint == 0) then
       write(*,*) 'No valid data entries in spectrum file ', trim(filename)
       stop
    end if

    allocate(spectrum(numpoint,2))
    if (present(sample_S_tab)) allocate(sample_S_tab(numpoint))
    i = 0
    open(unit, file=trim(filename))
    do while (.true.)
       read(unit,*,end=2) string

       if (string(1:1)=='#') cycle

       backspace(unit)
       if (present(sample_S_tab)) then
          read(unit,*,end=1) nu, val, sample_S_tab(i+1)
       else
          read(unit,*,end=1) nu, val
       end if
       i = i+1

       read(nu,*)  spectrum(i,1)
       read(val,*) spectrum(i,2)

       do j = 1, i-1
          if (spectrum(i,1) == spectrum(j,1)) then
             write(*,*) 'ERROR: Spectrum file ', trim(filename), ' contains double entries; j = ', j, ', spectrum = ', spectrum(i,1)
             stop
          end if
       end do
  
    end do
2   close(unit)

    ! Convert from GHz to Hz
    spectrum(:,1) = spectrum(:,1) * 1.d9

    if (present(nu_ref)) then
       ! Compute a spline through the spectrum spline
       allocate(x(numpoint))
       allocate(y(numpoint))
       allocate(y2(numpoint))
       
       x = spectrum(:,1)
       y = spectrum(:,2)

       call QuickSort_dp_dist(y, x)
       call spline(x, y, 1.d30, 1.d30, y2)
       
       ! Normalize to reference frequency
       S_ref = splint(x, y, y2, nu_ref)
       spectrum(:,2) = spectrum(:,2) / S_ref

       deallocate(x)
       deallocate(y)
       deallocate(y2)
    end if

  end subroutine read_spectrum

  ! Routine for reading a physical dust model file in RJ units (frequency is in Hz)
  subroutine read_phys_dust_spectrum(emission_file, u_file, spectrum, u_array, nu_ref)
    implicit none

    character(len=128),                         intent(in)  :: emission_file, u_file
    real(dp),           pointer, dimension(:,:)             :: spectrum
    real(dp),           pointer, dimension(:)               :: u_array
    real(dp),                                   intent(in)  :: nu_ref

    real(dp)                                  :: S_ref
    integer(i4b)                              :: i, j, unit, unit2, n_nu, n_u
    real(dp), allocatable, dimension(:)       :: u, nu
    real(dp), allocatable, dimension(:,:)     :: em_RJ, emission
    real(dp), allocatable, dimension(:,:,:,:) :: em_coeff

    unit  = getlun()
    open(unit, file=trim(emission_file))
    
    unit2 = getlun()
    open(unit2, file=trim(u_file))
        
    n_nu    = 500
    n_u     = 181

    allocate(spectrum(n_u+1, n_nu))
    allocate(u_array(n_u))
    allocate(u(n_u), nu(n_nu))
    allocate(emission(n_u, n_nu))
    allocate(em_RJ(n_u+1, n_nu), em_coeff(4, 4, n_u, n_nu))

    read(unit,*) ((em_RJ(i,j), i=1,n_u+1), j=1,n_nu)

    do i = 1, n_u
       read(unit2,fmt='(E31.24)') u(i)
    end do

    do i = 1, n_nu
       nu(i) = em_RJ(1,i)
       do j = 1, n_u
          emission(j,i) = log(em_RJ(j+1,i))!/maxval(em_RJ(j+1,:)))
       end do
    end do

    u_array        = u
    spectrum(1,:)  = nu(:)
    spectrum(2:,:) = emission(:,:)

    deallocate(u)
    deallocate(nu)
    deallocate(emission)
    deallocate(em_RJ)
    deallocate(em_coeff)

  end subroutine read_phys_dust_spectrum

  subroutine set_fg_params_equal(fg_par1, fg_par2)
    implicit none

    type(fg_params) :: fg_par1, fg_par2

    integer(i4b) :: i

    do i = 1, size(fg_par1%comp)
       if (fg_components(i)%npar > 0) fg_par2%comp(i)%p = fg_par1%comp(i)%p
    end do

  end subroutine set_fg_params_equal

  subroutine initialize_index_regions(chaindir, comp, par, filename, mask, mask_comp, nside, nmaps, &
       & fwhm, R, prior)
    implicit none

    character(len=*),                   intent(in)    :: chaindir, filename
    real(dp),         dimension(0:,1:), intent(in)    :: mask, mask_comp
    integer(i4b),                       intent(in)    :: comp, par, nside, nmaps
    real(dp),                           intent(in)    :: fwhm
    type(ind_region),                   intent(inout) :: R
    real(dp),         dimension(2),     intent(in)    :: prior

    integer(i4b) :: i, j, k, l, p, q, npix, m, n, nside_lowres, npix_lowres, p_nest, counter, unit, ierr, level
    integer(i4b) :: nlevel, pix, ntype, t
    real(dp)     :: threshold
    character(len=256) :: name, type_text, n_text, precompfile, fname
    character(len=2)   :: comp_text, par_text
    logical(lgt)       :: exist, new_region, ok
    integer(i4b), allocatable, dimension(:)   :: regset, ind
    real(dp),     allocatable, dimension(:)   :: fwhms
    real(dp),     allocatable, dimension(:)   :: numreg
    integer(i4b), allocatable, dimension(:,:) :: regdef
    real(dp),     allocatable, dimension(:,:) :: map, args, w_tot
    real(dp),     allocatable, dimension(:,:) :: order

    call int2string(comp, comp_text)
    call int2string(par,  par_text)
    precompfile = trim(chaindir) // '/regdef_c' // comp_text // '_p' // par_text // '.unf'
    inquire(file=precompfile, exist=exist) 
    unit = getlun()
    npix = 12*nside**2
    threshold = 1.d-3

    R%independent_pixels = .false.
    if (exist) then

       ! Read precomp file from disk
       open(unit, file=trim(precompfile), form='unformatted')
       read(unit) R%independent_pixels
       read(unit) R%nreg
       allocate(R%regions(R%nreg))
       do i = 1, R%nreg
          read(unit) R%regions(i)%n, R%regions(i)%n_ext
          allocate(R%regions(i)%pix(R%regions(i)%n,2), R%regions(i)%pix_ext(R%regions(i)%n_ext,2))
          allocate(R%regions(i)%w(R%regions(i)%n_ext))
          read(unit) R%regions(i)%pix
          read(unit) R%regions(i)%pix_ext
          read(unit) R%regions(i)%w
       end do
       read(unit) R%nregset
       allocate(R%regset(R%nregset))
       do i = 1, R%nregset
          read(unit) R%regset(i)%n
          allocate(R%regset(i)%r(R%regset(i)%n))
          read(unit) R%regset(i)%r
       end do
       close(unit)

       !write(*,*) myid, myid_chain, R%nreg, R%nregset, sum(real(R%regions(1)%pix_ext,dp))

    else

       read(filename,*) name
       if (trim(name) == 'fullsky') then
          if (trim(type_text) == 'TQU') then
             R%nreg = 1
          else if (trim(type_text) == 'T+QU') then
             write(*,*) 'Spectral type T+QU not yet supported'
             stop
          else if (trim(type_text) == 'T+Q+U') then
             write(*,*) 'Spectral type T+Q+U not yet supported'
             stop
          end if

          ! Define single pixel region
          R%nreg = 1
          allocate(R%regions(1))          
          if (sample_inside_mask) then
             n = count(mask_comp > 0.5d0)
          else
             n = count(mask_comp > 0.5d0 .and. mask > 0.5d0)
          end if
          if (n == 0) then
             R%nreg = 0
          else
             R%regions(1)%n     = n
             R%regions(1)%n_ext = n
             allocate(R%regions(1)%pix(n,2))
             k = 1
             do p = 0, npix-1
                do j = 1, nmaps
                   if (sample_inside_mask) then
                      ok = mask_comp(p,j) > 0.5d0
                   else
                      ok = mask_comp(p,j) > 0.5d0 .and. mask(p,j) > 0.5d0
                   end if
                   if (ok) then
                      R%regions(1)%pix(k,1)     = p
                      R%regions(1)%pix(k,2)     = j
                      k                         = k+1
                   end if
                end do
             end do
             allocate(fwhms(R%nreg))
             fwhms = fwhm
          end if

       else if (trim(name) == 'single_pix') then
          
          read(filename,*) name, type_text, n_text
          read(n_text,*) nside_lowres
          npix_lowres = 12*nside_lowres**2
          q           = npix / npix_lowres
          R%independent_pixels = (fwhm == 0) .and. q == 1
          
          if (trim(type_text) == 'TQU') then
             R%nreg = npix_lowres
             ntype  = 1
          else if (trim(type_text) == 'T+QU') then
             R%nreg = 2*npix_lowres
             ntype  = 2
          else if (trim(type_text) == 'T+Q+U') then
             R%nreg = 3*npix_lowres
             ntype  = 3
          end if
          
          
          allocate(R%regions(2*R%nreg))
          counter = 0
          do t = 1, ntype
             do p = 1, R%nreg/ntype
                !Find number of unmasked pixels in current region
                n = 0
                do i = 0, q-1
                   if (trim(type_text) == 'TQU') then
                      do j = 1, nmaps
                         call nest2ring(nside, (p-1)*q+i, k)
                         if (mask(k,j) > 0.5d0 .and. mask_comp(k,j) > 0.5d0) n = n+1
                      end do
                   else if (trim(type_text) == 'T+QU') then
                      if (t == 1) then
                         do j = 1, 1
                            call nest2ring(nside, (p-1)*q+i, k)
                            if (mask(k,j) > 0.5d0 .and. mask_comp(k,j) > 0.5d0) n = n+1
                         end do
                      else
                         do j = 2, nmaps
                            call nest2ring(nside, (p-1)*q+i, k)
                            if (mask(k,j) > 0.5d0 .and. mask_comp(k,j) > 0.5d0) n = n+1
                         end do
                      end if
                   else if (trim(type_text) == 'T+Q+U') then
                      do j = t, t
                         call nest2ring(nside, (p-1)*q+i, k)
                         if (mask(k,j) > 0.5d0 .and. mask_comp(k,j) > 0.5d0) n = n+1
                      end do
                   end if
                end do
                if (n > 0) then
                   counter = counter+1
                   R%regions(counter)%n = n
                   allocate(R%regions(counter)%pix(n,2))
                   n = 1
                   do i = 0, q-1
                      if (trim(type_text) == 'TQU') then
                         do j = 1, nmaps
                            call nest2ring(nside, (p-1)*q+i, k)
                            if (mask(k,j) > 0.5d0 .and. mask_comp(k,j) > 0.5d0) then
                               R%regions(counter)%pix(n,1) = k
                               R%regions(counter)%pix(n,2) = j
                               n = n+1
                            end if
                         end do
                      else if (trim(type_text) == 'T+QU') then
                         if (t == 1) then
                            do j = 1, 1
                               call nest2ring(nside, (p-1)*q+i, k)
                               if (mask(k,j) > 0.5d0 .and. mask_comp(k,j) > 0.5d0) then
                                  R%regions(counter)%pix(n,1) = k
                                  R%regions(counter)%pix(n,2) = j
                                  n = n+1
                               end if
                            end do
                         else
                            do j = 2, nmaps
                               call nest2ring(nside, (p-1)*q+i, k)
                               if (mask(k,j) > 0.5d0 .and. mask_comp(k,j) > 0.5d0) then
                                  R%regions(counter)%pix(n,1) = k
                                  R%regions(counter)%pix(n,2) = j
                                  n = n+1
                               end if
                            end do
                         end if
                      else if (trim(type_text) == 'T+Q+U') then
                         do j = t, t
                            call nest2ring(nside, (p-1)*q+i, k)
                            if (mask(k,j) > 0.5d0 .and. mask_comp(k,j) > 0.5d0) then
                               R%regions(counter)%pix(n,1) = k
                               R%regions(counter)%pix(n,2) = j
                               n = n+1
                            end if
                         end do
                      end if
                   end do
                end if
                
                if (sample_inside_mask) then
                   !Find number of masked pixels in current region
                   n = 0
                   do i = 0, q-1
                      if (trim(type_text) == 'TQU') then
                         do j = 1, nmaps
                            call nest2ring(nside, (p-1)*q+i, k)
                            if (mask(k,j) <= 0.5d0 .and. mask_comp(k,j) > 0.5d0) n = n+1
                         end do
                      else if (trim(type_text) == 'T+QU') then
                         if (t == 1) then
                            do j = 1, 1
                               call nest2ring(nside, (p-1)*q+i, k)
                               if (mask(k,j) <= 0.5d0 .and. mask_comp(k,j) > 0.5d0) n = n+1
                            end do
                         else
                            do j = 2, nmaps
                               call nest2ring(nside, (p-1)*q+i, k)
                               if (mask(k,j) <= 0.5d0 .and. mask_comp(k,j) > 0.5d0) n = n+1
                            end do
                         end if
                      else if (trim(type_text) == 'T+Q+U') then
                         do j = t, t
                            call nest2ring(nside, (p-1)*q+i, k)
                            if (mask(k,j) <= 0.5d0 .and. mask_comp(k,j) > 0.5d0) n = n+1
                         end do
                      end if
                   end do
                   if (n > 0) then
                      counter = counter+1
                      R%regions(counter)%n = n
                      allocate(R%regions(counter)%pix(n,2))
                      n = 1
                      do i = 0, q-1
                         if (trim(type_text) == 'TQU') then
                            do j = 1, nmaps
                               call nest2ring(nside, (p-1)*q+i, k)
                               if (mask(k,j) <= 0.5d0 .and. mask_comp(k,j) > 0.5d0) then
                                  R%regions(counter)%pix(n,1) = k
                                  R%regions(counter)%pix(n,2) = j
                                  n = n+1
                               end if
                            end do
                         else if (trim(type_text) == 'T+QU') then
                            if (t == 1) then
                               do j = 1, 1
                                  call nest2ring(nside, (p-1)*q+i, k)
                                  if (mask(k,j) <= 0.5d0 .and. mask_comp(k,j) > 0.5d0) then
                                     R%regions(counter)%pix(n,1) = k
                                     R%regions(counter)%pix(n,2) = j
                                     n = n+1
                                  end if
                               end do
                            else
                               do j = 2, nmaps
                                  call nest2ring(nside, (p-1)*q+i, k)
                                  if (mask(k,j) <= 0.5d0 .and. mask_comp(k,j) > 0.5d0) then
                                     R%regions(counter)%pix(n,1) = k
                                     R%regions(counter)%pix(n,2) = j
                                     n = n+1
                                  end if
                               end do
                            end if
                         else if (trim(type_text) == 'T+Q+U') then
                            do j = t, t
                               call nest2ring(nside, (p-1)*q+i, k)
                               if (mask(k,j) <= 0.5d0 .and. mask_comp(k,j) > 0.5d0) then
                                  R%regions(counter)%pix(n,1) = k
                                  R%regions(counter)%pix(n,2) = j
                                  n = n+1
                               end if
                            end do
                         end if
                      end do
                   end if

                end if
             end do
          end do


          R%nreg = counter ! Allow for masked/unmasked region in each case
          allocate(fwhms(R%nreg))
          fwhms = fwhm
          
!!$    else if (trim(name) == 'single_pix_T+QU') then
!!$
!!$       R%nreg = 2*npix
!!$       allocate(R%regions(R%nreg))
!!$       do p = 1, npix
!!$          ! Temperature
!!$          R%regions(p)%n = 1
!!$          allocate(R%regions(p)%pix(R%regions(p)%n,2))
!!$          R%regions(p)%pix(1,1) = p-1
!!$          R%regions(p)%pix(1,2) = 1
!!$          ! Polarization
!!$          R%regions(p+npix)%n = 2
!!$          allocate(R%regions(p+npix)%pix(R%regions(p+npix)%n,2))
!!$          R%regions(p+npix)%pix(1,1) = p-1
!!$          R%regions(p+npix)%pix(1,2) = 2
!!$          R%regions(p+npix)%pix(2,1) = p-1
!!$          R%regions(p+npix)%pix(2,2) = 3
!!$       end do
!!$
!!$    else if (trim(name) == 'single_pix_T+Q+U') then
!!$
!!$       R%nreg = 3*npix
!!$       allocate(R%regions(R%nreg))
!!$       do p = 1, npix
!!$          ! Temperature
!!$          R%regions(p)%n = 1
!!$          allocate(R%regions(p)%pix(R%regions(p)%n,2))
!!$          R%regions(p)%pix(1,1) = p-1
!!$          R%regions(p)%pix(1,2) = 1
!!$          ! Stokes Q
!!$          R%regions(p+npix)%n = 1
!!$          allocate(R%regions(p+npix)%pix(R%regions(p+npix)%n,2))
!!$          R%regions(p+npix)%pix(1,1) = p-1
!!$          R%regions(p+npix)%pix(1,2) = 2
!!$          ! Stokes U
!!$          R%regions(p+2*npix)%n = 1
!!$          allocate(R%regions(p+2*npix)%pix(R%regions(p+2*npix)%n,2))
!!$          R%regions(p+2*npix)%pix(1,1) = p-1
!!$          R%regions(p+2*npix)%pix(1,2) = 3
!!$       end do

       else if (trim(name) == 'S2N_map') then
          
          read(filename,*) name, fname, type_text, nlevel
          allocate(args(3,nlevel))
          read(filename,*) name, fname, type_text, nlevel, args

          allocate(map(0:npix-1,nmaps), fwhms(npix))
          call read_map(fname, map)          

          counter = 0
          allocate(R%regions(npix))
          do level = 1, nlevel
             npix_lowres = 12*nint(args(2,level))**2
             if (npix_lowres == 0) then
                ! Put everything into one region
                n = 0
                do p = 0, npix-1
                   if (map(p,1) >= args(1,level) .and. mask_comp(p,1) > 0.5d0) n = n+1
                end do

                if (n > 0) then
                   counter              = counter+1
                   R%regions(counter)%n = n*nmaps
                   fwhms(counter)       = args(3,level)
                   allocate(R%regions(counter)%pix(n*nmaps,2))
                   k = 1
                   do p = 0, npix-1
                      do l = 1, nmaps
                         if (map(p,1) >= args(1,level) .and. mask_comp(p,1) > 0.5d0) then
                            R%regions(counter)%pix(k,1) = p
                            R%regions(counter)%pix(k,2) = l
                            k                           = k+1
                            map(p,1)                    = -1.d0
                         end if
                      end do
                   end do
                end if

             else
                q           = npix / npix_lowres
                do i = 0, npix_lowres-1
                   ! Find number of pixels
                   n = 0
                   do j = 0, q-1
                      pix = i*q+j
                      call nest2ring(nside, pix, p)
                      if (map(p,1) >= args(1,level) .and. mask_comp(p,1) > 0.5d0) n = n+1
                   end do
                   
                   if (n > 0) then
                      counter              = counter+1
                      R%regions(counter)%n = n*nmaps
                      fwhms(counter)       = args(3,level)
                      allocate(R%regions(counter)%pix(n*nmaps,2))
                      k = 1
                      do j = 0, q-1
                         do l = 1, nmaps
                            pix = i*q+j
                            call nest2ring(nside, pix, p)
                            if (map(p,1) >= args(1,level) .and. mask_comp(p,1) > 0.5d0) then
                               R%regions(counter)%pix(k,1) = p
                               R%regions(counter)%pix(k,2) = l
                               k                           = k+1
                               map(p,1)                    = -1.d0
                            end if
                         end do
                      end do
                   end if
                end do
             end if
          end do
          deallocate(map)
          R%nreg = counter
          
       else
          
          allocate(map(0:npix-1,nmaps), regdef(0:npix-1,nmaps))
          call read_map(filename, map)          
          regdef = nint(map)
          
          ! Set up regions
          R%nreg = maxval(map)
          allocate(R%regions(2*R%nreg))
          counter = 0
          do i = 1, R%nreg
             ! Find number of unmasked pixels in this region
             n = 0
             do j = 1, nmaps
                do p = 0, npix-1
!                   if (regdef(p,j) == i .and. mask(p,j) > 0.5d0 .and. mask_comp(p,j) > 0.5d0) n=n+1
                   if (regdef(p,j) == i .and. mask_comp(p,j) > 0.5d0) n=n+1
                end do
             end do
             
             if (n > 0) then
                counter = counter+1
                R%regions(counter)%n = n
                allocate(R%regions(counter)%pix(n,2))
                k = 1
                do p = 0, npix-1
                   do j = 1, nmaps
!                      if (regdef(p,j) == i .and. mask(p,j) > 0.5d0 .and. mask_comp(p,j) > 0.5d0) then
                      if (regdef(p,j) == i .and. mask_comp(p,j) > 0.5d0) then
                         R%regions(counter)%pix(k,1) = p
                         R%regions(counter)%pix(k,2) = j
                         k                           = k+1
                      end if
                   end do
                end do
             end if
             
!!$             if (sample_inside_mask) then
!!$                
!!$                ! Find number of unmasked pixels in this region
!!$                n = 0
!!$                do j = 1, nmaps
!!$                   do p = 0, npix-1
!!$                      if (regdef(p,j) == i .and. mask(p,j) <= 0.5d0 .and. mask_comp(p,j) > 0.5d0) n=n+1
!!$                   end do
!!$                end do
!!$                
!!$                if (n > 0) then
!!$                   counter = counter+1
!!$                   R%regions(counter)%n = n
!!$                   allocate(R%regions(counter)%pix(n,2))
!!$                   k = 1
!!$                   do p = 0, npix-1
!!$                      do j = 1, nmaps
!!$                         if (regdef(p,j) == i .and. mask(p,j) <= 0.5d0 .and. mask_comp(p,j) > 0.5d0) then
!!$                            R%regions(counter)%pix(k,1) = p
!!$                            R%regions(counter)%pix(k,2) = j
!!$                            k                           = k+1
!!$                         end if
!!$                      end do
!!$                   end do
!!$                end if
!!$             end if
             
          end do
          R%nreg = counter
          allocate(fwhms(R%nreg))
          fwhms = fwhm
          
          deallocate(map, regdef)
          
       end if

!!$       write(*,*) myid_chain, R%nreg
!!$       if (myid_chain == root) then
!!$          allocate(map(0:npix-1,1))
!!$          map = -1.6375d30
!!$          do i = 1, R%nreg
!!$             map(R%regions(i)%pix(:,1),1) = R%regions(i)%n
!!$          end do
!!$          call write_map('regions_raw.fits', map)
!!$       end if
!!$       call mpi_finalize(ierr)
!!$       stop

       ! Set up smooth weights
       allocate(map(0:npix-1,nmaps))
       do i = 1+myid, R%nreg, numprocs
          if (fwhms(i) > 0.d0 .and. R%nreg > 1) then
             map = 0.d0
             do j = 1, R%regions(i)%n
                map(R%regions(i)%pix(j,1),R%regions(i)%pix(j,2)) = 1.d0
             end do
             do j = 1, nmaps
                if (any(map(:,j) /= 0.d0)) then
                   call smooth_par_map(map(:,j:j), fwhms(i), [0.d0, 1.d0])
                end if
             end do
             
             n = count(map > threshold)
             R%regions(i)%n_ext = n
             allocate(R%regions(i)%w(n), R%regions(i)%pix_ext(n,2))
             l = 0
             do k = 0, npix-1
                do j = 1, nmaps
                   if (map(k,j) > threshold) then
                      l = l+1
                      R%regions(i)%pix_ext(l,1) = k
                      R%regions(i)%pix_ext(l,2) = j
                      R%regions(i)%w(l)         = map(k,j)
                   end if
                end do
             end do
          else
             n = R%regions(i)%n
             R%regions(i)%n_ext = n
             allocate(R%regions(i)%w(n), R%regions(i)%pix_ext(n,2))
             R%regions(i)%w       = 1.d0
             R%regions(i)%pix_ext = R%regions(i)%pix
          end if
       end do
       deallocate(map)
       
       ! Synchronize
       do i = 1, R%nreg
          call mpi_bcast(R%regions(i)%n_ext, 1,   MPI_INTEGER, mod(i-1,numprocs), MPI_COMM_WORLD, ierr)
          if (.not. allocated(R%regions(i)%w))       allocate(R%regions(i)%w(R%regions(i)%n_ext))
          if (.not. allocated(R%regions(i)%pix_ext)) allocate(R%regions(i)%pix_ext(R%regions(i)%n_ext,2))
          call mpi_bcast(R%regions(i)%w,       size(R%regions(i)%w),       MPI_DOUBLE_PRECISION, &
               & mod(i-1,numprocs), MPI_COMM_WORLD, ierr)
          call mpi_bcast(R%regions(i)%pix_ext, size(R%regions(i)%pix_ext), MPI_INTEGER,          &
               & mod(i-1,numprocs), MPI_COMM_WORLD, ierr)
       end do

       ! Normalize weights
       allocate(w_tot(0:npix-1,nmaps))
       w_tot = 0.d0
       do i = 1, R%nreg
          do j = 1, R%regions(i)%n_ext
             p = R%regions(i)%pix_ext(j,1)
             l = R%regions(i)%pix_ext(j,2)
             w_tot(p,l) = w_tot(p,l) + R%regions(i)%w(j)
          end do
       end do
       do i = 1, R%nreg
          do j = 1, R%regions(i)%n_ext
             p = R%regions(i)%pix_ext(j,1)
             l = R%regions(i)%pix_ext(j,2)
             R%regions(i)%w(j) = R%regions(i)%w(j) / w_tot(p,l)
          end do
       end do
       deallocate(w_tot)

       ! Define region sets, ie., sets of strictly disjoint regions
       if (R%nreg == 1 .or. all(fwhms <= 0.d0)) then
          R%nregset = 1
          allocate(R%regset(1))
          R%regset(1)%n = R%nreg
          allocate(R%regset(1)%r(R%nreg))
          do i = 1, R%nreg
             R%regset(1)%r(i) = i
          end do
       else
          n = 0
          allocate(regset(R%nreg))
          regset = -1
          do i = 1, R%nreg
             ! Sort sets according to number of elements; try to put in the least populated
             allocate(numreg(n), ind(n))
             if (n > 0) then
                do j = 1, n
                   ind(j)    = j
                   numreg(j) = count(regset(1:i-1) == j)
                end do
                call QuickSort(ind, numreg)
             end if

             do j = 1, n
                ! Attempt to put region i in region set j
                do k = 1, i-1
                   if (regset(k) == ind(j)) then
                      ! Look for duplicate pixels
                      l = 1
                      m = 1
                      do while (l <= R%regions(i)%n_ext .and. m <= R%regions(k)%n_ext)
                         if (R%regions(i)%pix_ext(l,1) == R%regions(k)%pix_ext(m,1) .and. &
                              & R%regions(i)%pix_ext(l,2) == R%regions(k)%pix_ext(m,2)) then
                            ! Failed -- duplicates exist
                            goto 100
                         end if
                         if (R%regions(i)%pix_ext(l,1) <= R%regions(k)%pix_ext(m,1)) then
                            l = l+1
                         else
                            m = m+1
                         end if
                      end do
!!$                      do l = 1, R%regions(i)%n_ext
!!$                         do m = 1, R%regions(k)%n_ext
!!$                            if (R%regions(i)%pix_ext(l,1) == R%regions(k)%pix_ext(m,1) .and. &
!!$                                 & R%regions(i)%pix_ext(l,2) == R%regions(k)%pix_ext(m,2)) then
!!$                               ! Failed -- duplicates exist
!!$                               goto 100
!!$                            end if
!!$                         end do
!!$                      end do
                   end if
                end do
                ! No duplicates were found; put in region set j
                regset(i) = ind(j)
                exit
100             continue
             end do
             deallocate(ind, numreg)
             if (j > n) then
                ! Don't fit in any available region sets; create a new one
                n = n+1
                regset(i) = n
             end if
101          continue
          end do
          R%nregset = n
          allocate(R%regset(n))
          do i = 1, n
             m = count(regset == i)
             R%regset(i)%n = m
             allocate(R%regset(i)%r(m))
             j = 0
             do k = 1, R%nreg
                if (regset(k) == i) then
                   j = j+1
                   R%regset(i)%r(j) = k
                end if
             end do
          end do
          deallocate(regset)
       end if

       ! Barrier that actually works...
       call mpi_allreduce(i, j, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
       if (myid == root) then
          ! Write precomp file to disk
          open(unit, file=trim(precompfile), form='unformatted')
          write(unit) R%independent_pixels
          write(unit) R%nreg
          do i = 1, R%nreg
             write(unit) R%regions(i)%n, R%regions(i)%n_ext
             write(unit) R%regions(i)%pix
             write(unit) R%regions(i)%pix_ext
             write(unit) R%regions(i)%w
          end do
          write(unit) R%nregset
          do i = 1, R%nregset
             write(unit) R%regset(i)%n
             write(unit) R%regset(i)%r
          end do
          close(unit)
       end if

    end if
    ! Barrier that actually works...
    call mpi_allreduce(i, j, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    
    if (myid_chain == root) then
       write(*,*) '   Number of index region sets = ', R%nregset
    end if

    if (myid == root) then
       allocate(map(0:npix-1,R%nregset))
       map = -1.6375d30
       do i = 1, R%nregset
          do j = 1, R%regset(i)%n
             !write(*,*) i, j, R%regset(i)%r(j)
             !write(*,*) R%regions(R%regset(i)%r(j))%pix_ext(:,1)
             !write(*,*) R%regions(R%regset(i)%r(j))%w
             map(R%regions(R%regset(i)%r(j))%pix_ext(:,1),i) = R%regions(R%regset(i)%r(j))%w
          end do
       end do
       fname = trim(chaindir)//'/regsets_c' // comp_text // '_p' // par_text // '.fits'
       call write_map(fname, map(:,1:min(R%nregset,3)))
    end if
    !call mpi_finalize(ierr)
    !stop
    
    ! Initialize recent point list over prior
    do k = 1, R%nreg
       do i = 1, num_recent_point
          R%regions(k)%recent(i) = max(min((i-1)*(prior(2)-prior(1))/(num_recent_point-1.d0),prior(2)),prior(1))
       end do
    end do
    
!!$    if (myid_chain == root) then
!!$       write(*,*) R%regset(1)%r(1:2)
!!$       write(*,*) R%regions(R%regset(1)%r(1))%pix_ext(:,1)
!!$       write(*,*)
!!$       write(*,*) R%regions(R%regset(1)%r(2))%pix_ext(:,1)
!!$    end if

!    call mpi_finalize(ierr)
!    stop

  end subroutine initialize_index_regions

  subroutine get_smooth_par_map(par, par_smooth)
    implicit none

    real(dp), dimension(0:,1:,1:), intent(in)   :: par
    real(dp), dimension(0:,1:,1:), intent(out)  :: par_smooth

    integer(i4b) :: i, j, k, l, p, pix, spec
    real(dp)     :: p0

    p = 1
    do i = 1, num_fg_comp
       do j = 1, fg_components(i)%npar
          par_smooth(:,:,p) = 0.d0
          do k = 1, fg_components(i)%indregs(j)%nreg
             pix  = fg_components(i)%indregs(j)%regions(k)%pix(1,1)
             spec = fg_components(i)%indregs(j)%regions(k)%pix(1,2)
             p0   = par(pix,spec,p)
             do l = 1, fg_components(i)%indregs(j)%regions(k)%n_ext
                pix  = fg_components(i)%indregs(j)%regions(k)%pix_ext(l,1)
                spec = fg_components(i)%indregs(j)%regions(k)%pix_ext(l,2)
                par_smooth(pix,spec,p) = par_smooth(pix,spec,p) + &
                     & fg_components(i)%indregs(j)%regions(k)%w(l) * p0
             end do
          end do
          par_smooth(:,:,p) = max(min(par_smooth(:,:,p), fg_components(i)%priors(j,2)), &
               & fg_components(i)%priors(j,1))
          p = p+1
       end do
    end do

  end subroutine get_smooth_par_map

end module comm_fg_component_mod
