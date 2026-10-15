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
submodule (comm_tod_hfi_mod) comm_tod_hfi_smod
contains

  !**************************************************
  !             Constructor
  !**************************************************
  module function constructor_hfi(cpar, id, id_abs, info, tod_type) result(c)
    !
    ! Constructor function that gathers all the instrument parameters in a pointer
    ! and constructs the objects
    !
    ! Arguments:
    ! ----------
    ! cpar:     derived type
    !           Object containing parameters from the parameterfile.
    ! id_abs:   integer
    !           The index of the current band within the parameters, related to cpar
    ! info:     map_info structure
    !           Information about the maps for this band, like how the maps are distributed in memory
    ! tod_type: string
    !           Instrument specific tod type
    !
    ! Returns
    ! ----------
    ! constructor: pointer
    !              Pointer that contains all instrument data

    implicit none
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs
    class(comm_mapinfo),     target     :: info
    character(len=128),      intent(in) :: tod_type
    class(comm_HFI_tod),     pointer    :: c

    integer(i4b) :: i, j, k, nside_beam, lmax_beam, nmaps_beam, ierr
    logical(lgt) :: pol_beam
    character(len=6) :: pstring

    logical(lgt), dimension(:,:), allocatable :: correlations

    call timer%start(TOD_INIT, id_abs)
    
    ! Allocate object
    allocate(c)

    ! Initialize parameters needed for general initialization
    c%nhorn           = 1
    c%samprate_lowres = 18.  ! Lowres samprate in Hz;  10 times lower than the intrinsic HFI rate for now    
    c%nmaps           = info%nmaps
    c%ndet            = num_tokens(cpar%ds_tod_dets(id_abs), "," )
    c%noise_psd_model = 'oof'       ! Not fitted parameters yet
    !c%noise_psd_model = 'spline'

    ! Initialize common parameters
    call c%tod_constructor(cpar, id, id_abs, info, tod_type)
    
    ! Initialize instrument-specific parameters
    c%reload_TOD        = .true.
    c%compressed_tod    = .true.
    c%correct_sl        = .true.
    c%correct_orb       = .true.
    c%correct_S_crosstalk = .false.
    c%correct_N_crosstalk = .false.
    c%apply_inst_corr   = .true.
    c%orb_4pi_beam      = .true.
    c%symm_flags        = .true.
    c%sample_zodi     = cpar%sample_zodi .and. c%subtract_zodi ! Sample zodi parameters
    c%ntime           = 1
    !TODO: set the number of dark bolometers to be correct
    c%ndark           = 1
    c%n_cray_temps    = 3
    c%ndiode          = 1
    nmaps_beam        = 3
    pol_beam          = .true.
    c%nside_beam      = 512
    c%nside_pixhist   = 64
    c%sol_elong_range = cpar%zs_sol_elong
    
    ! Set up noise PSD type and priors
    c%freq            = cpar%ds_label(id_abs)
    if (c%noise_psd_model == 'spline') then
       c%n_xi            = 13 !94 ! {sigma0, spline nodes}
       allocate(c%xi_n_P_uni(c%n_xi,2))
       allocate(c%xi_n_P_rms(c%n_xi))
       allocate(c%xi_n_nu_fit(c%n_xi,2))

       c%xi_n_P_uni(1,:)  = [10.d0, 3000.d0]  ! Sigma0
       c%xi_n_P_rms(1)    = 50.d0
       c%xi_n_nu_fit(1,:) = [10.d0, 90.d0]
       do i = 2, c%n_xi
          c%xi_n_P_uni(i,:)  = [-6.d0, 12.d0]
          c%xi_n_P_rms(i)    = 6.d0
          c%xi_n_nu_fit(i,:) = [0.d0, 90.d0]
       end do
    else
       c%n_xi            = 3
       !c%noise_psd_model = 'white'       ! Using white noise until we get better estimates of the actual noise PSD
       allocate(c%xi_n_P_uni(c%n_xi,2))
       allocate(c%xi_n_P_rms(c%n_xi))
       allocate(c%xi_n_nu_fit(c%n_xi,2))

       c%xi_n_P_uni(1,:)  = [10d0, 300d0]  ! Sigma0
       c%xi_n_P_uni(2,:)  = [0.001d0, 1d0]  ! fknee
       c%xi_n_P_uni(3,:)  = [-2.5d0, -0.5d0]   ! alpha
       !c%xi_n_P_uni(4,:)  = [ 0.5d0,  4.0d0]  ! fknee
       !c%xi_n_P_uni(5,:)  = [-1.5d0, -0.5d0]   ! alpha
       c%xi_n_nu_fit(1,:) = [0.5d0, 5d0] 
       c%xi_n_nu_fit(2,:) = [0.001d0, 0.5d0]
       c%xi_n_nu_fit(3,:) = [0.001d0, 0.5d0]
       !c%xi_n_nu_fit(4,:) = [0.001d0, 10d0]
       !c%xi_n_nu_fit(5,:) = [0.001d0, 10d0]
       c%xi_n_P_rms       = [10.d0, 0.1d0, 0.1d0] ! [sigma0, fknee, alpha]; sigma0 is not used
    end if
    c%f_spin           = 1./60.                 ! Planck spin frequency in Hz
    
    !c%xi_n_P_rms      = [-1.d0] ! [sigma0]; sigma0 is not used
    c%n_cray_temps    = 3
    c%ndiode = 1

    ! Initialize common parameters
    !call c%tod_constructor(cpar, id, id_abs, info, tod_type)

    ! Initialize instrument-specific parameters
    c%samprate_lowres = 1.  ! Lowres samprate in Hz;  10 times lower than the intrinsic HFI rate for now
    c%nhorn           = 1
    c%rawtod_dp       = .true.
    c%compressed_tod  = .true.
    c%correct_sl      = .true.
    c%correct_orb     = .true.
    c%apply_inst_corr = .true.
    c%correct_Tbol    = .true.
    c%orb_4pi_beam    = .true.
    c%symm_flags      = .true.
    c%sample_zodi     = cpar%sample_zodi .and. c%subtract_zodi ! Sample zodi parameters
    c%nmaps           = info%nmaps
    c%ndet            = num_tokens(cpar%ds_tod_dets(id_abs), "," )
    c%ntime           = 1
    c%max_npole_Tbol  = 7
    !TODO: set the number of dark bolometers to be correct
    c%ndark           = 1

    nmaps_beam        = 3
    pol_beam          = .true.
    c%nside_beam      = 512
    c%nside_pixhist   = 64

    c%gain_samprate    = 1.d0 / 3600.d0 * 100.d0 ! Gain smoothing stiffness
        
    ! Channel specific parameters
    if (c%freq(1:3) == "100" .or. c%freq(1:3) == "143".or. c%freq(1:3) == "217" .or. c%freq(1:3) == "353") then
       c%chisq_threshold  = 1.d30  ! 120d0
       c%sigma0_threshold = 1.d30  !100.d0
       c%accept_threshold = 0.8d0
       c%correct_sl       = .true.
    else if (c%freq(1:3) == "545") then
       c%chisq_threshold  = 1.d30  ! 120d0
       c%sigma0_threshold = 1.d30  !100.d0
       c%accept_threshold = 0.8d0
       c%orb_4pi_beam     = .false.
       c%correct_sl       = .false.
    else if (c%freq(1:3) == "857") then
       c%xi_n_P_uni(1,:)  = [10d0, 3000d0]  ! Sigma0
       c%chisq_threshold  = 1.d30  ! 120d0
       c%sigma0_threshold = 1.d30  !100.d0
       c%accept_threshold = 0.8d0
       c%orb_4pi_beam     = .false.
       c%correct_sl       = .false.
    else
       write(*,*) "Specify HFI defaults"
       call mpi_finalize(ierr)
       stop
    end if
    
    ! Get detector labels
    call get_tokens(cpar%ds_tod_dets(id_abs), ",", c%label)

    ! Identify partners
    c%partner = -1
    do i = 1, c%ndet
       if (len(trim(c%label(i))) == 5) cycle ! Spider web
       ! Search for detector with a instead of b
       pstring = trim(adjustl(c%label(i)))
       if (pstring(6:6) == 'a') then
          pstring(6:6) = 'b'
       else
          pstring(6:6) = 'a'
       end if
       do j = 1, c%ndet
          if (pstring == trim(adjustl(c%label(j)))) then
             c%partner(i) = j
             !if (c%myid == 0) write(*,*) 'Partners = ', trim(adjustl(c%label(i))), '<->', trim(adjustl(c%label(c%partner(i))))
             exit
          end if
       end do
    end do
    
    ! Read the actual TOD
    call c%read_tod(c%label)
    
    ! Initialize bandpass mean and proposal matrix
    call c%initialize_bp_covar(cpar%ds_tod_bp_init(id_abs))

    ! Construct lookup tables
    !if (.not. c%reload_TOD) then
       c%pixcache => comm_tod_pixcache(c%nside, c%nside_beam, c%nmaps, .false.)
       !call c%precompute_lookups()
    !end if

    ! Allocate and initialize bolometer transfer functions
    if (c%correct_Tbol) allocate(c%Tbol(c%ndet))
    
    ! Load the instrument file
    call c%load_instrument_file(c%nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

    ! Collect Sun velocities from all scans
    call c%collect_v_sun
    
    ! Allocate sidelobe convolution data structures
    if (c%correct_sl) allocate(c%slconv(c%ndet,c%nhorn), c%orb_dp)
    !allocate(c%orb_dp)
    c%orb_dp => comm_orbdipole(c%mbeam)  
    !c%orb_dp => comm_orbdipole(comm=c%comm)

    ! Initialize all baseline corrections to zero
    do i = 1, c%nscan
       do j = 1, c%ndet
          c%scans(i)%d(j)%baseline = 0.d0
       end do
    end do

    ! Allocate modulation phase
    allocate(c%mod_phase(c%ndet,c%nscan))
    c%mod_phase = 1.0
   
    ! initialize crosstalk class
    if (c%correct_N_crosstalk) then
       allocate(correlations(c%ndet, c%ndet))
       correlations = .true.    

       !  allocate(c%xtalk)
       c%xtalk => comm_crosstalk(correlations)
    end if

    ! Pre-initialize ADC object
    allocate(c%adc(c%ndet))
    allocate(c%adu_range(c%ndet,2))
    
    ! initialize 4K_lines objects
    c%n_4K_lines = 15
    allocate(c%nus_4K_lines(c%n_4K_lines))
    do i = 1, 8
       c%nus_4K_lines(i) = 180.3737d0 / 18.d0 * i  ! Harmonics
    end do
    c%nus_4K_lines(9) = 16.7d0  ! One of the  main lines
    c%nus_4K_lines(10:15) = [16.d0, 25.d0, 43.4d0, 46.2d0, 47.6d0, 56.8d0] ! NPIPE secondary lines
    c%cooler_correction_domain = 'time' ! 'frequency'

    if (trim(c%cooler_correction_domain) == 'time') then ! Estimation in time domain
       allocate(c%cooler_4K_lines_profile(c%ndet,c%nscan))
    elseif (trim(c%cooler_correction_domain) == 'frequency') then ! In frequency domain
       allocate(c%cooler_4K_lines(c%n_4K_lines,c%ndet,c%nscan))
       do i = 1, c%nscan
          do j = 1, c%ndet
             do k = 1, c%n_4K_lines
                c%cooler_4K_lines(k,j,i)%p => comm_4K_lines()
             end do
          end do
       end do
    end if

    ! Initialize dynamic mask
    c%dynmask => comm_dynmask(c, cpar)
    !c%dynmask%apply_pixhist           = .true.
    !c%dynmask%apply_solar_mask        = .true.
    !c%dynmask%remove_isolated_samples = .true.
    if (c%freq(1:3) /= "353") then
       c%dynmask%threshold_extreme      = 20.0  ! in units of white noise sigma
       c%dynmask%threshold_singlesamp   =  5.0  ! in units of residual sigma
       !c%dynmask%threshold_excessRMS(1) =  1.5  ! in units of residual sigma
       !c%dynmask%window_excessRMS(1)    =    5  ! window size in number of samples
       !c%dynmask%threshold_excessRMS(2) =  2.0  ! in units of residual sigma
       !c%dynmask%window_excessRMS(2)    =   50  ! window size in number of samples
    end if
    
    call timer%stop(TOD_INIT, id_abs)

  end function constructor_hfi

  !**************************************************
  !             Driver routine
  !**************************************************
  module subroutine process_HFI_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
    !
    ! Routine that processes the HFI time ordered data.
    ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
    ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms.
    ! Writes maps to disc in fits format
    !
    ! Arguments:
    ! ----------
    ! self:     pointer of comm_HFI_tod class
    !           Points to output of the constructor
    ! chaindir: string
    !           Directory for output files
    ! chain:    integer
    !           Index number of the chain being run
    ! iter:     integer
    !           Gibbs iteration number
    ! handle:   planck_rng derived type
    !           Healpix definition for random number generation
    !           so that the same sequence can be resumed later on from that same point
    ! map_in:   array
    !           Array of dimension (ndet,ndelta) with pointer to maps,
    !           with both access to maps and changing them.
    !           ndet is the number of detectors and
    !           ndelta is the number of bandpass deltas being considered
    ! delta:    array
    !           Array of bandpass corrections with dimensions (0:ndet,npar,ndelta)
    !           where ndet is number of detectors, npar is number of parameters
    !           and ndelta is the number of bandpass deltas being considered
    !
    ! Returns:
    ! ----------
    ! map_out: comm_map class
    !          Final output map after TOD processing combined for all detectors
    ! rms_out: comm_map class
    !          Final output rms map after TOD processing combined for all detectors

    implicit none
    class(comm_HFI_tod),                      intent(inout) :: self
    character(len=*),                         intent(in)    :: chaindir
    integer(i4b),                             intent(in)    :: chain, iter
    type(planck_rng),                         intent(inout) :: handle
    type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
    real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
    class(comm_map),                          intent(inout) :: map_out      ! Combined output map
    class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms
    type(map_ptr),       dimension(1:),       intent(inout), optional :: map_gain       ! (ndet)

    real(dp)            :: t1, t2
    integer(i4b)        :: i, j, k, h, l, ierr, ndelta, nside, npix, nmaps, dec_wn, oper_default, skip_nonlin_, seed
    logical(lgt)        :: select_data, output_scanlist, output_zodi_comps
    logical(lgt)        :: sample_gain, sample_ncorr, sample_abs_bandpass, sample_rel_bandpass, sample_zodi, sample_adc, make_dyn_mask, sample_xi_n
    logical(lgt)        :: fit_4K_lines
    class(comm_binmap), pointer   :: binmap
    type(comm_scandata) :: sd
    !type(comm_detdata)  :: dd
    character(len=4)    :: ctext, myid_text
    character(len=6)    :: samptext, scantext, itertext, dettext
    character(len=512)  :: prefix, postfix, prefix4D, filename, Sfilename
    character(len=512), allocatable, dimension(:) :: slist
    real(sp),              dimension(9)       :: flag_threshold
    real(sp), allocatable, dimension(:)       :: procmask, procmask2, procmask_zodi, sigma0, freqmask
    real(sp), allocatable, dimension(:)       :: d_prime
    real(sp), allocatable, dimension(:,:)     :: s_buf
    real(sp), allocatable, dimension(:,:,:)   :: d_calib
    !real(sp), allocatable, dimension(:,:,:,:) :: map_sky, m_gain
    real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf
    class(comm_cgmap), pointer :: cgmap
    
    ! file for saving tods
    type(hdf_file) :: tod_file

    call int2string(iter, ctext)
    call update_status(status, "tod_start"//ctext)
    call timer%start(TOD_TOT, self%band)
    
    if (self%first_call) then
       seed = rand_uni(handle) * 100000000
       call rand_init(self%handle, seed)
    end if

    if (self%reload_TOD) then
       do i = 1, self%nscan
          call read_hdf_scan_data(self%scans(i), self, self%hdfname(i), self%scanid(i),&
               & self%ndet, self%label, self%nhorn, self%ndiode, self%diode_names)
       end do

       ! Construct lookup tables
       call self%precompute_lookups()
    end if
    
    call map_in(1,1)%p%writeFITS(trim(self%outdir) // "/input_sky_model_"//trim(self%label(1))//".fits")
    
    call update_status(status, "tod_read"//ctext)
    
    ! Toggle optional operations
    sample_rel_bandpass   = size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
    sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses
    if (.false.) then ! Debug
       ! Do data selection, then start sampling
       sample_gain           = .false. !iter  > 0 !.true.                 
       make_dyn_mask         = iter == 1
       sample_ncorr          = iter  > 0 !.true.
       sample_xi_n      = .false.
       select_data           = iter == 1
       sample_adc            = .false. !.false. !iter  > 1 !.true.
    else if (trim(self%init_from_HDF) == 'none') then
       ! Initialize slowly if not HDF init
       sample_gain           = iter  > 2 !.true.                 
       make_dyn_mask         = iter == 20
       sample_ncorr          = iter > 10 !.true.
       sample_xi_n           = iter > 15 
       select_data           = iter == 25 ! self%first_call  
       sample_adc            = .false. !iter  > 0 ! 3 !.true.
    else
       ! Do data selection, then start sampling
       sample_gain           = iter > 1
       make_dyn_mask         = iter == 1
       sample_ncorr          = iter > 1 !.true.
       sample_xi_n           = iter > 1 !.false.
       select_data           = .false. !iter == 1 ! self%first_call  
       sample_adc            = .false. !iter  > 0 !.true.
    end if
    fit_4K_lines          = .false. !iter > 2
    if (self%freq(1:3) == "545" .or. self%freq(1:3) == "857") make_dyn_mask = .false.

    sample_zodi           = self%sample_zodi .and. self%subtract_zodi ! Sample zodi parameters
    output_zodi_comps     = self%output_zodi_comps .and. self%subtract_zodi ! Output zodi components
    output_scanlist       = mod(iter-1,10) == 0    ! only output scanlist every 10th iteration
    dec_wn                = 2 ! Decimation factor for sigma0; 2 corresponds to 45Hz
    skip_nonlin_ = 100

    if (sample_ncorr) then
       if (self%correct_sl) then
          oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
               & SD_SKY,SD_BP,SD_ORB,SD_INST,SD_DARK,SD_NCORR,SD_SL])
       else
          oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
               & SD_SKY,SD_BP,SD_ORB,SD_INST,SD_DARK,SD_NCORR])
       end if
       !oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
       !     & SD_SKY,SD_BP,SD_ORB,SD_INST,SD_DARK,SD_NCORR])
    else
       if (self%correct_sl) then
           oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
               & SD_SKY,SD_BP,SD_ORB,SD_INST,SD_DARK,SD_SL])
       else
           oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
               & SD_SKY,SD_BP,SD_ORB,SD_INST,SD_DARK])
       end if
       !oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
       !     & SD_SKY,SD_BP,SD_ORB,SD_INST,SD_DARK])
    end if
    
    ! Initialize local variables
    ndelta          = size(delta,3)
    self%n_bp_prop  = ndelta-1
    nside           = map_out%info%nside
    nmaps           = map_out%info%nmaps
    npix            = 12*nside**2
    self%output_n_maps = 3
    if (self%output_aux_maps > 0) then
       if (mod(iter-1,self%output_aux_maps) == 0) self%output_n_maps = 7
    end if
    if (output_zodi_comps) self%output_n_maps = 8 + zodi_model%n_comps

    call int2string(chain, ctext)
    call int2string(iter, samptext)
    call int2string(self%myid, myid_text)
    prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
    postfix = '_c' // ctext // '_k' // samptext // '.fits'

    ! Initialize index-based sky map and mask
    call self%pixcache%init_map_mask(map_in, self%bitmask, map_gain=map_gain)
    call update_status(status, "tod_cache"//ctext)

    !init procmask
    if (sample_rel_bandpass .or. sample_abs_bandpass) then
       allocate(m_buf(0:npix-1,1), procmask2(0:npix-1))
       call self%procmask%bcast_fullsky_map(m_buf)
       procmask2 = real(m_buf(:,1),sp)
       deallocate(m_buf)
    end if

    ! Precompute far sidelobe Conviqt structures
    if (self%correct_sl) then
       call timer%start(TOD_SL_PRE, self%band)
       if (self%myid == 0) write(*,*) 'Precomputing sidelobe convolved sky'
       do i = 1, self%ndet
          !TODO: figure out why this is rotated
          call map_in(i,1)%p%YtW()  ! Compute sky a_lms
          self%slconv(i,1)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
               & 100, 3, 100, self%slbeam(i)%p, map_in(i,1)%p, 2)
       end do
       call timer%stop(TOD_SL_PRE, self%band)
    end if

    call update_status(status, "tod_init")

    !------------------------------------
    ! Perform main sampling steps
    !------------------------------------

    ! Sample ADC and baseline parameters jointly
    if (sample_adc) then
       ! Initialize ADC objects before first call
       if (.not. associated(self%adc(1)%p)) then
          call self%compute_adu_range
          do i = 1, self%ndet
             self%adc(i)%p => comm_adc_binfit(self%comm, self%datadir, self%outdir, self%label(i), 16, &
                  & self%adu_range(i,1), self%adu_range(i,2), 40)
          end do
       end if

       ! Sample ADC parameters -- MUST BE FOLLOWED BY BASELINE SAMPLER
       ! HKE: Comment out for now
       !do i = 1, self%ndet
       !   call self%sample_adc_and_baselines(handle, i)
       !end do
       call update_status(status, "tod_adc"//ctext)
    end if
!!$    
    ! Fit per-chunk low-level non-linearity parameters
    do i = 1, self%nscan ! Disable for now
       ! Skip scan if no accepted data
       if (.not. any(self%scans(i)%d%accept)) cycle

!        open(58,file=trim(chaindir) // '/res.dat', recl=1024)
!        do j = 1, sd%ntod
!           write(58,*) j, sd%tod(j,1), sd%mask(j,1), sd%flag(j,1)
!        end do
!        close(58)

       ! Apply only ADC corrections
       call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd, nonlin_level=1)

       ! Estimate modulation baselines; and set modulation phase
       if (self%first_call .and. trim(self%init_from_HDF) == 'none') then
          call sample_hfi_baselines(sd, self, i, handle, subtract_s_tot=.false.)
          call set_modulation_phase(sd, self, i)
       else
          call sample_hfi_baselines(sd, self, i, handle)
       end if
       call demodulate_tod(sd, self, i)
       
       ! Estimate pre-deconvolution white noise rms
       call sample_noise_psd(self, sd, handle, chaindir, only_sigma0=.true., dec_wn=dec_wn, sigma0_preproc=.true.)
          
       call timer%stop(TOD_NONLIN, self%band)

       ! Fix dc level jumps 
       call self%stitch_hfi_dc_level(i, sd)

       ! Dark bolometer drift correction
       call self%hfi_dark_correction(i, sd)       


       if (fit_4K_lines) then
          if (trim(self%cooler_correction_domain) == 'time') then
             ! 4K Line corrections from TOD
             do j = 1, self%ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                call estimate_hfi_4K_lines_tod(self, sd, j)
             end do
          elseif (trim(self%cooler_correction_domain) == 'frequency') then
             ! In-paint flagged samples for 4K lines estimation
             do j = 1, self%ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                do k = 1, sd%ntod
                   if (iand(sd%flag(k,j), self%flag0) .ne. 0) then
                      sd%tod(k,j) = self%scans(i)%d(j)%gain * sd%s_tot(k,j,0,1) + self%scans(i)%d(j)%N_psd%sigma0_preproc * rand_gauss(handle)
                   end if
                end do
             end do   
             ! 4K Line corrections from PSD
             do j = 1, self%ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                call estimate_hfi_4K_lines_harm(self, sd, j)
             end do
          end if
       end if


       ! Initialize spline noise model
       if (self%noise_psd_model == 'spline') then
          do j = 1, self%ndet
             if (self%scans(i)%d(j)%accept) then
                call whitening_filter(self, sd, j)
                call update_spline_noise_psd(self,sd,i,j)
             end if
          end do
         ! if (self%myid==0) write(*,*) '|  Number of spline noise model parameters:', self%scans(i)%d(1)%N_psd%npar
       end if

       ! Clean up
       call dealloc_scan_data(sd)
    end do
    call timer%start(TOD_WAIT, self%band)
    call mpi_barrier(self%comm, ierr) ! Improve timing information
    call timer%stop(TOD_WAIT, self%band)

    call update_status(status, "tod_nonlin"//ctext)

    ! Fit global timestream contaminants 

    ! Subtract cosmic ray contribution
    !do j=1, self%ndet

    !  call init_det_data_singlehorn(dd, self, j)

    !  call self%cray(j)%p%build_cray_templates()

    !  do i=1, self%nscan
    !    call populate_sd_from_dd(sd, dd, i, j)

    !    call self%cray(j)%p%fit_cray_amplitudes(sd%tod(j,:), sd%s_inst(j, :))

    !    call dealloc_scan_data(sd)
    !  end do

    !  call dd%dealloc
    !end do

    ! Estimate ADC corrections
    !    Not implemented yet

    ! Fit bolometer transfer function parameters?
    !    Not implemented yet     

    ! NOW Sample high level tod components that require cleaned data

    ! Sample gain components in separate TOD loops; marginal with respect to n_corr
    if (sample_gain) then
       ! TODO: Also sample non-linear gain response here?
       call sample_calibration(self, 'abscal', oper_default, handle)
       call sample_calibration(self, 'relcal', oper_default, handle)
       !call sample_calibration(self, 'deltaG', oper_default, handle, smooth=.true.)
       !call sample_calibration(self, 'deltaG', oper_default, handle, smooth=.false.)
       !call sample_calibration(self, 'total', oper_default, handle, smooth=.false.)
       call update_status(status, "tod_calib"//ctext)
    end if
    
    ! Create pixel histograms
    !if (self%first_call) call compute_tod_pixhist(self)
    !call update_status(status, "tod_hist"//ctext)
    
    ! Prepare intermediate data structures
    !call binmap%init(self, .true., .false., nplus2=.false.)
    binmap => comm_binmap(self, .true., sample_rel_bandpass, nplus2=.false.)
    if (sample_abs_bandpass .or. sample_rel_bandpass) then
       allocate(chisq_S(self%ndet,size(delta,3)))
       chisq_S = 0.d0
    end if
    if (output_scanlist) then
       allocate(slist(self%nscan))
       slist   = ''
    end if

    ! Initialize CG mapmaker, maptype = 2 = TQU
    !cgmap => comm_cgmap(self, 2)

    ! Fit higher-level corrections
    if (self%myid == 0) write(*,*) '   --> Sampling ncorr, xi_n, maps'
    call update_status(status, "tod_preloop"//ctext)
    do i = 1, self%nscan

       ! Skip scan if no accepted data
       if (.not. any(self%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd, handle=handle)
       
       ! Create dynamic mask
       if (make_dyn_mask) then
          ! Estimate sigma0 for masking
          if (trim(self%init_from_HDF) == 'none') then
             call sample_noise_psd(self, sd, handle, chaindir, only_sigma0=.true., dec_wn=dec_wn)
          end if
          
          ! Create mask
          do j = 1, sd%ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             call self%dynmask%create(sd, j)
          end do

          call dealloc_scan_data(sd)
          if (.not. any(self%scans(i)%d%accept)) cycle

          ! Update scan data with new flagging
          call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd, handle=handle)
       end if

       ! Sample correlated noise
       if (sample_ncorr) then
          call sample_n_corr(self, sd, handle)
          if (sample_xi_n) then
             allocate(freqmask(0:sd%ntod/2))
             freqmask = 1.
             call create_spin_freqmask(real(self%samprate,sp), self%f_spin, self%f_spin/10., 3.0, freqmask) ! Spin harmonics
             !call create_spin_freqmask(real(self%samprate,sp), 10., 1., 85., freqmask) ! 4K harmonics
             
             call sample_noise_psd(self, sd, handle, chaindir, freqmask=freqmask)
             !call sample_noise_psd(self, sd, handle, chaindir)
             deallocate(freqmask)
          else
             call sample_noise_psd(self, sd, handle, chaindir, only_sigma0=.true., dec_wn=dec_wn)
          end if
       else
          call sample_noise_psd(self, sd, handle, chaindir, only_sigma0=.true., dec_wn=dec_wn)
       end if
       
       ! Compute chisquare
       do j = 1, sd%ndet
          if (.not. self%scans(i)%d(j)%accept) cycle
          call self%compute_tod_chisq(sd, j)
       end do

       ! Select data
       if (select_data) then
          call remove_bad_data(self, i, sd%flag)
       end if

       ! Compute chisquare for bandpass fit
       if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)

       ! Compute calibrated TOD for mapmaking
       allocate(d_calib(binmap%nout,sd%ntod, sd%ndet))
       d_calib = 0.d0
       call compute_calibrated_data(self, i, sd, d_calib)

       ! Bin TOD
       call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib, binmap)

       ! Feed CG mapmaker calibrated and cleaned data
       !call cgmap%load_data(i, self, sd, d_calib(1,:,:))

      ! output tod for debugging
       ! for some reason the first iteration is outputing as "tod_"
       !if (mod(self%scanid(i),5000) == 0 .and. allocated(sd%n_corr)) then
       if (mod(self%scanid(i),5000) == 0) then
          do j = 1, sd%ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             sd%tod(:,j)    = sd%tod(:,j) / self%scans(i)%d(j)%gain
             !sd%n_corr(:,j) = sd%n_corr(:,j) / self%scans(i)%d(j)%gain             
          end do


          call int2string(self%scanid(i), scantext)

          write(*,*) '| Writing tod to hdf'
          call open_hdf_file(trim(chaindir)//'/tod_'//scantext//'_samp'//samptext//'.h5', tod_file, 'w')
          call write_hdf(tod_file, '/tod',      sd%tod)
          call write_hdf(tod_file, '/calib',    d_calib(1,:,:))
          call write_hdf(tod_file, '/res',      d_calib(2,:,:))
          call write_hdf(tod_file, '/flag',     sd%flag)
          !call write_hdf(tod_file, '/n_corr',   sd%n_corr)
          call write_hdf(tod_file, '/orb',      sd%s_orb)
          call write_hdf(tod_file, '/psi',      sd%psi)
          call write_hdf(tod_file, '/pix',      sd%pix)
          call write_hdf(tod_file, '/s_tot',    sd%s_tot)
          call write_hdf(tod_file, '/mask',     sd%mask)
          call close_hdf_file(tod_file)
        end if
       
       ! Update scan list
       call wall_time(t2)
       self%scans(i)%proctime   = self%scans(i)%proctime   + t2-t1
       self%scans(i)%n_proctime = self%scans(i)%n_proctime + 1
       if (output_scanlist) then
          write(slist(i),*) self%scanid(i), '"',trim(self%hdfname(i)), &
               & '"', real(self%scans(i)%proctime/self%scans(i)%n_proctime,sp),&
               & real(self%spinaxis(i,:),sp)
       end if

       ! Clean up
       call dealloc_scan_data(sd)
       deallocate(d_calib)
    end do

    call timer%start(TOD_WAIT, self%band)
    call mpi_barrier(self%comm, ierr) ! Improve timing information
    call timer%stop(TOD_WAIT, self%band)
    call update_status(status, "tod_postloop"//ctext)
    
    if (select_data) then
       ! Remove data based on a gliding RMS window cut for each of the listed
       ! criteria

       !if (trim(self%noise_psd_model) /= 'spline') then
          !                           half-window  [chisq, sigma0, fknee, alpha, base, base1, base2]
          !call remove_tod_outliers(self, 100,      [5.,    5.,     5.,    5.,    0.,   5.,    5.   ])
       !else
          !                                  half-window  [chisq, sigma0, xi_n,  base,  base1, base2]
          !call remove_tod_outliers_spline(self, 100,      [5.,    5.,     9.,    0.,    5.,    5.   ])
       !end if
       
       if (self%symm_flags) then
          ! Remove partners for rejected scans
          do j = 1, self%ndet
             if (self%partner(j) == -1) cycle
             do i = 1, self%nscan
                if (.not. self%scans(i)%d(j)%accept) &
                     & self%scans(i)%d(self%partner(j))%accept = .false.
             end do
          end do
       end if
          
       do i = 1, self%nscan
          do j = 1, self%ndet
             if (.not. self%scans(i)%d(j)%accept) self%scans(i)%d(j)%nsamp_unmasked = 0
          end do
       end do
    end if
    call update_status(status, "tod_outlier"//ctext)

    
    if (self%myid == 0) write(*,*) '   --> Finalizing maps, bp'
    if (make_dyn_mask) then
       call self%dynmask%report
       call update_status(status, "tod_dynstats"//ctext)
    end if
    
    ! Output latest scan list with new timing information
    if (output_scanlist) call self%output_scan_list(slist)

    ! Solve for maps
    ! TODO: update mapmaker to n+2 maps
    ! TODO: handle bolometer transfer function in the mapmaking
    call synchronize_binmap(binmap, self)
    call update_status(status, "tod_binmap1"//ctext)
    if (sample_rel_bandpass) then
       if (self%nmaps > 1) then
          !call finalize_binned_map_nplus2(self, binmap, rms_out, 1.d0, chisq_S=chisq_S, mask=procmask2, correct_transfer=.true.)
          Sfilename = trim(prefix) // 'Smap'// trim(postfix)
          call finalize_binned_map(self, binmap, rms_out, 1.d0, chisq_S=chisq_S, Sfilename=Sfilename, mask=procmask2)
       else
         call finalize_binned_map_unpol(self, binmap, rms_out, 1.d0, chisq_S=chisq_s, mask=procmask2, correct_transfer=.true.)
       end if
    else
       if(self%nmaps > 1) then
          call finalize_binned_map(self, binmap, rms_out, 1.d0)
       else 
         call finalize_binned_map_unpol(self, binmap, rms_out, 1.d0)
       end if
    end if
    map_out%map = binmap%outmaps(1)%p%map
    call update_status(status, "tod_binmap2"//ctext)

    ! Sample bandpass parameters
    if (sample_rel_bandpass .or. sample_abs_bandpass) then
       call sample_bp(self, handle, chisq_S, delta)
       self%bp_delta = delta(:,:,1)
    end if
   
    ! Output maps to disk
    call timer%start(TOD_WRITE)
    call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
    call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
    if (self%output_n_maps > 1) call binmap%outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
    if (self%output_n_maps > 2) call binmap%outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
    if (self%output_n_maps > 3) call binmap%outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix))
    if (self%output_n_maps > 4) call binmap%outmaps(5)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
    if (self%output_n_maps > 5) call binmap%outmaps(6)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
    if (self%output_n_maps > 6) call binmap%outmaps(7)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))
    if (self%output_n_maps > 8 .and. self%subtract_zodi .and. output_zodi_comps) then
       do i = 1, zodi_model%n_comps
          call binmap%outmaps(8+i)%p%writeFITS(trim(prefix)//'zodi_'//trim(zodi_model%comp_labels(i))//trim(postfix))
       end do
    endif
    call timer%stop(TOD_WRITE)
    call update_status(status, "tod_binmap3"//ctext)

    ! Solve for CG map
    !call cgmap%solve(map_out)
    !call dealloc_cgmap(cgmap)
    !call map_out%writeFITS(trim(prefix)//'cgmap'//trim(postfix))
    
    ! Clean up
    call deallocate_binmap(binmap)
    call update_status(status, "tod_binmap4"//ctext)
    if (allocated(slist)) deallocate(slist)
    if (self%correct_sl) then
       do i = 1, self%ndet
          do h = 1, self%nhorn
             call self%slconv(i,h)%p%dealloc(); deallocate(self%slconv(i,h)%p)
          end do
       end do
    end if
    
    if (self%reload_TOD) then
       call self%pixcache%dealloc
       do i = self%nscan, 1, -1
          do j = self%ndet, 1, -1
             if (allocated(self%scans(i)%d(j)%flag))  deallocate(self%scans(i)%d(j)%flag)
             do k = self%nhorn, 1, -1
                if (allocated(self%scans(i)%d(j)%pix(k)%p)) deallocate(self%scans(i)%d(j)%pix(k)%p)
                if (allocated(self%scans(i)%d(j)%psi(k)%p)) deallocate(self%scans(i)%d(j)%psi(k)%p)
             end do
             if (allocated(self%scans(i)%d(j)%psi)) deallocate(self%scans(i)%d(j)%psi)
             if (allocated(self%scans(i)%d(j)%pix)) deallocate(self%scans(i)%d(j)%pix)
             if (allocated(self%scans(i)%d(j)%zdiode)) deallocate(self%scans(i)%d(j)%zdiode)
             if (allocated(self%scans(i)%d(j)%diode)) deallocate(self%scans(i)%d(j)%diode)
             if (allocated(self%scans(i)%d(j)%ztod))  deallocate(self%scans(i)%d(j)%ztod)
             if (allocated(self%scans(i)%d(j)%tod))   deallocate(self%scans(i)%d(j)%tod)
          end do
       end do
    end if

    ! Parameter to check if this is first time routine has been called
    if (allocated(procmask2)) deallocate(procmask2)
    if (allocated(chisq_S)) deallocate(chisq_S)
    self%first_call = .false.

    call update_status(status, "tod_end"//ctext)
    call timer%stop(TOD_TOT, self%band)
    
  end subroutine process_HFI_tod

  module subroutine load_instrument_hfi(self, instfile, band)
    !
    ! Reads the HFI specific fields from the instrument file
    ! Implements comm_tod_mod::load_instrument_inst
    !
    ! Arguments:
    !
    ! self : comm_HFI_tod
    !    the HFI tod object (this class)
    ! file : hdf_file
    !    the open file handle for the instrument file
    ! band : int
    !    the index of the current detector
    ! 
    ! Returns : None
    implicit none
    class(comm_hfi_tod),                 intent(inout) :: self
    type(hdf_file),                      intent(in)    :: instfile
    integer(i4b),                        intent(in)    :: band

    integer(i4b) :: j
    real(dp), allocatable, dimension(:) :: par
    
    call read_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'polEff', self%pol_eff(band))
    self%pol_eff(band) = self%pol_eff(band) * 0.01d0 ! Stored as percentage in the instrument file for now

    if (self%correct_Tbol) then
       allocate(par(2 + 2*self%max_npole_Tbol))
       do j = 1, self%ndet
          call read_hdf(instfile, trim(adjustl(self%label(j)))//'/'//'Tbol', par)
          self%Tbol(j)%p => comm_Tbol(self%samprate, par)
       end do
       deallocate(par)
    end if

  end subroutine load_instrument_hfi


  module subroutine sample_hfi_baselines(self, tod, scan, handle, subtract_s_tot)
    ! 
    ! Estimates baselines for MODULATED data, separate for odd and even samples
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_scandata)
    !           HFI-specific TOD object
    ! tod:      comm_tod derived type
    !             contains TOD-specific information         
    ! scan:     scan ID
    ! handle:   planck_rng derived type
    !           Healpix definition for random number generation
    !           so that the same sequence can be resumed later on from that same
    !           point
    !           
    !
    ! Returns
    ! ----------
    !   None, but updates tod%scans(scan)%d(:)%baseline  (for odd samples)
    !                     tod%scans(scan)%d(:)%baseline2 (for even samples)
    !
    implicit none
    class(comm_scandata),                 intent(in)    :: self
    class(comm_hfi_tod),                  intent(inout) :: tod
    integer(i4b),                         intent(in)    :: scan
    type(planck_rng),                     intent(inout) :: handle
    logical(lgt),                         intent(in), optional :: subtract_s_tot
    
    real(dp) :: eta, A1, A2, x, b1, b2, sgn
    integer(i4b) :: i, j
    logical(lgt) :: sub_s

    call timer%start(TOD_BASELINE, self%band)
    
    sub_s = .true.; if (present(subtract_s_tot)) sub_s = subtract_s_tot
    
    
    ! tod%scans(scan)%d(i)%gain - the gain constant over a scan [real number]
    ! sd = self --- self%s_tot - sky signal model
    ! self%s_tot(self%ntod, self%ndet) - how s_tot structured

    !if (tod%scanid(scan) == 1151) open(58,file='baseline_new.dat', recl=1024)
    do i = 1, tod%ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle
       sgn = tod%mod_phase(i,scan)
       
       ! Odd samples
       A1 = 0.d0; b1 = 0
       do j = 1, self%ntod, 2
          if (self%mask(j,i) == 0) cycle
          A1 = A1 + 1.d0
          b1 = b1 + self%tod_dp(j,i)
          if (sub_s) b1 = b1 - sgn*tod%scans(scan)%d(i)%gain * self%s_tot(j,i,0,1)
          !if (tod%scanid(scan) == 1151) write(58,*) j, self%tod(j,i), sgn, tod%scans(scan)%d(i)%gain, self%s_tot(j,i,0,1)
       end do
       A1 = A1 / tod%scans(scan)%d(i)%N_psd%sigma0**2
       b1 = b1 / tod%scans(scan)%d(i)%N_psd%sigma0**2
       if (A1 == 0.d0) then
          tod%scans(scan)%d(i)%accept = .false.
          cycle
       end if
       tod%scans(scan)%d(i)%baseline1  = b1/A1 + rand_gauss(handle)/sqrt(A1)
       !write(*,*) 'Ab1', tod%scanid(scan), i, b1, real(A1,sp), real(rand_gauss(handle)/sqrt(A1),sp)
       
       ! Even samples
       !if (tod%scanid(scan) == 1151) write(58,*)
       A2 = 0.d0; b2 = 0.d0
       do j = 2, self%ntod, 2
          if (self%mask(j,i) == 0) cycle
          A2 = A2 + 1.d0
          b2 = b2 + self%tod_dp(j,i)
          if (sub_s) b2 = b2 + sgn*tod%scans(scan)%d(i)%gain * self%s_tot(j,i,0,1)
          !if (tod%scanid(scan) == 1151) write(58,*) j, self%tod(j,i), sgn, tod%scans(scan)%d(i)%gain, self%s_tot(j,i,0,1)
       end do
       A2 = A2 / tod%scans(scan)%d(i)%N_psd%sigma0**2
       b2 = b2 / tod%scans(scan)%d(i)%N_psd%sigma0**2
       if (A2 == 0.d0) then
          tod%scans(scan)%d(i)%accept = .false.
          cycle
       end if       
       tod%scans(scan)%d(i)%baseline2  = b2/A2 + rand_gauss(handle)/sqrt(A2)
       !write(*,*) 'Ab2', tod%scanid(scan), i, b2, real(A2,sp), real(rand_gauss(handle)/sqrt(A2),sp)
       
       !write(*,'(a,i8,3f16.3)') "baseline =", tod%scanid(scan), tod%scans(scan)%d(i)%baseline1, tod%scans(scan)%d(i)%baseline2, sgn

    end do
    !if (tod%scanid(scan) == 1151) close(58)

    call timer%stop(TOD_BASELINE, self%band)
    
  end subroutine sample_hfi_baselines

  module subroutine set_modulation_phase(self, tod, scan)
    ! 
    ! Estimates baselines for MODULATED data, separate for odd and even samples
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_scandata)
    !           HFI-specific TOD object
    ! tod:      comm_tod derived type
    !             contains TOD-specific information         
    ! scan:     scan ID
    ! handle:   planck_rng derived type
    !           Healpix definition for random number generation
    !           so that the same sequence can be resumed later on from that same
    !           point
    !           
    !
    ! Returns
    ! ----------
    !   None, but updates tod%scans(scan)%d(:)%baseline1  (for odd samples)
    !                     tod%scans(scan)%d(:)%baseline2 (for even samples)
    !                     tod%scans(scan)%d(:)%gain (temporary solution)
    !
    implicit none
    class(comm_scandata),                 intent(in)    :: self
    class(comm_hfi_tod),                  intent(inout) :: tod
    integer(i4b),                         intent(in)    :: scan
    
    real(dp) :: mu1, mu2, n
    integer(i4b) :: i, j
    
    do i = 1, tod%ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle       
       
       mu1 = 0.d0; n = 0.d0
       do j = 1, self%ntod-1, 2
          if (iand(self%flag(j,  i), tod%flag0) .ne. 0) cycle
          if (iand(self%flag(j+1,i), tod%flag0) .ne. 0) cycle
          if (self%pix(j,i,1) > 0.48*tod%info%npix .and. self%pix(j,i,1) < 0.52*tod%info%npix) then
             mu1 = mu1 + (self%tod_dp(j,  i)-tod%scans(scan)%d(i)%baseline1) - &
                       & (self%tod_dp(j+1,i)-tod%scans(scan)%d(i)%baseline2)
             n   = n  + 1.d0
          end if
       end do

       if (n == 0) then
          write(*,*) "set_modulation_phase: no samples crossing the galactic plane. Scan disabled = ", tod%scanid(scan)
          tod%scans(scan)%d(i)%accept = .false.
       else
          mu1 = mu1/n
       end if
       
       ! saving the phase to the tod object
       if (mu1 < 0.) then
          tod%mod_phase(i,scan) = -1
       end if
    end do

  end subroutine set_modulation_phase

  module subroutine demodulate_tod(self, tod, scan)
    ! 
    ! Demodulate HFI TOD
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_scandata)
    !           HFI-specific TOD object
    ! tod:      comm_tod derived type
    !             contains TOD-specific information         
    ! scan:     Scan ID
    !
    ! Returns
    ! ----------
    !   None, but updates self%tod
    !
    implicit none
    class(comm_scandata),                         intent(inout) :: self
    class(comm_hfi_tod),                          intent(in)    :: tod
    integer(i4b),                                 intent(in)    :: scan

    integer(i4b) :: i, j, d
    real(sp)     :: sgn

!!$    open(58,file='tod_adc.dat')
!!$    do j = 1, self%ntod
!!$       write(58,*) j, self%tod(j,1)
!!$    end do
!!$    close(58)
    
    do i = 1, self%ndet
       d = self%det(i)
       if (.not. tod%scans(scan)%d(d)%accept) cycle       
       sgn = tod%mod_phase(d,scan)
       
       ! Subtract baselines and flip sign of every other sample
       do j = 1, self%ntod
           if (mod(j,2) == 1) then
              self%tod(j,i) =  sgn*(self%tod_dp(j,i) - tod%scans(scan)%d(d)%baseline1)
              !self%tod(j,i) =  (self%tod(j,i) - tod%scans(scan)%d(i)%baseline1)
           else
              self%tod(j,i) = -sgn*(self%tod_dp(j,i) - tod%scans(scan)%d(d)%baseline2)
              !self%tod(j,i) = -(self%tod(j,i) - tod%scans(scan)%d(i)%baseline2)
           end if
       end do
       !if (sgn < 0.) self%tod(:,i) = -self%tod(:,i)
       
    end do

!!$    open(58,file='tod.dat')
!!$    do j = 1, self%ntod
!!$       write(58,*) j, self%tod(j,1), self%mask(j,1)
!!$    end do
!!$    close(58)
    
  end subroutine demodulate_tod

  module subroutine compute_adu_range(self)
    ! 
    ! Computes ADU range over all scans and unmasked samples; for ADC correction
    ! Must be called after dynamic mask definition
    !
    ! Arguments:
    ! ----------
    ! self:      comm_tod derived type
    !             contains TOD-specific information         
    ! Returns
    ! ----------
    !   None, but updates tod%adu_range
    !
    implicit none
    class(comm_hfi_tod),                  intent(inout) :: self
    
    integer(i4b) :: i, j, ntod, scan, ierr, oper
    type(comm_scandata) :: sd
    
    oper = get_sd_operation_code([SD_BASE,SD_TOD])

    do i = 1, self%ndet
       self%adu_range(i,:) = [40*2**16,0]
    end do
    
    do scan = 1, self%nscan
       do i = 1, self%ndet
          if (.not. self%scans(scan)%d(i)%accept) cycle
          call init_scan_data(self, scan, oper, -1, sd, nonlin_level=0)
          do j = 1, sd%ntod
             if (iand(sd%flag(j,i), self%flag0) .eq. 0) then
                self%adu_range(i,1) = min(self%adu_range(i,1), nint(sd%tod_dp(j,i)))
                self%adu_range(i,2) = max(self%adu_range(i,2), nint(sd%tod_dp(j,i)))
             end if
          end do
          call dealloc_scan_data(sd)
       end do
    end do

    call mpi_allreduce(MPI_IN_PLACE, self%adu_range(:,1), self%ndet, &
         & MPI_INTEGER, MPI_MIN, self%comm, ierr)
    call mpi_allreduce(MPI_IN_PLACE, self%adu_range(:,2), self%ndet, &
         & MPI_INTEGER, MPI_MAX, self%comm, ierr)
    self%adu_range(:,1) = self%adu_range(:,1) / 40 - 16 ! Add a little buffer
    self%adu_range(:,2) = self%adu_range(:,2) / 40 + 16
    
    if (self%myid == 0) then
       write(*,*) '  ADU range       min        max        n_ADU'
       do i = 1, self%ndet
          write(*,*) '   ', trim(self%label(i)), self%adu_range(i,:), self%adu_range(i,2)-self%adu_range(i,1)
       end do
    end if

!!$    call mpi_finalize(i)
!!$    stop
    
  end subroutine compute_adu_range


  
  module subroutine read_tod_inst_HFI(self, file)
    ! 
    ! Reads HFI-specific common fields from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_HFI_tod)
    !           HFI-specific TOD object
    ! file:     derived type (hdf_file)
    !           Already open HDF file handle; only root includes this
    !
    ! Returns
    ! ----------
    ! None, but updates self
    !
    implicit none
    class(comm_HFI_tod),                 intent(inout)          :: self
    type(hdf_file),                      intent(in),   optional :: file
  end subroutine read_tod_inst_HFI
  
  module subroutine read_scan_inst_HFI(self, file, slabel, detlabels, scan)
    ! 
    ! Reads HFI-specific scan information from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_HFI_tod)
    !           HFI-specific TOD object
    ! file:     derived type (hdf_file)
    !           Already open HDF file handle
    ! slabel:   string
    !           Scan label, e.g., "000001/"
    ! detlabels: string (array)
    !           Array of detector labels, e.g., ["27M", "27S"]
    ! scan:     derived class (comm_scan)
    !           Scan object
    !
    ! Returns
    ! ----------
    ! None, but updates scan object
    !
    implicit none
    class(comm_HFI_tod),                 intent(in)    :: self
    type(hdf_file),                      intent(in)    :: file
    character(len=*),                    intent(in)    :: slabel
    character(len=*), dimension(:),      intent(in)    :: detlabels
    class(comm_scan),                    intent(inout) :: scan
  end subroutine read_scan_inst_HFI

  module subroutine initHDF_HFI(self, chainfile, path)
    ! 
    ! Initializes HFI-specific TOD parameters from existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_HFI_tod)
    !           HFI-specific TOD object
    ! chainfile: derived type (hdf_file)
    !           Already open HDF file handle to existing chainfile
    ! path:   string
    !           HDF path to current dataset, e.g., "000001/tod/030"
    !
    ! Returns
    ! ----------
    ! None
    !
    implicit none
    class(comm_HFI_tod),                 intent(inout)  :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path

    integer(i4b) :: ierr, i, j, k
    real(dp), allocatable, dimension(:,:,:) :: base
    real(sp), allocatable, dimension(:,:)   :: phase
    real(sp), allocatable, dimension(:,:)   :: Q

    allocate(base(self%nscan_tot,self%ndet,2), phase(self%nscan_tot,self%ndet))
    if (self%myid == 0) then
       call read_hdf(chainfile, trim(adjustl(path))//'baseline_adc', base)
       call read_hdf(chainfile, trim(adjustl(path))//'mod_phase',    phase)
       call read_hdf(chainfile, trim(adjustl(path))//'adu_range', self%adu_range)

       ! Read ADC tables
       if (self%adu_range(1,1) > 0) then
          !allocate(Q(minval(self%adu_range(:,1)):maxval(self%adu_range(:,2)),self%ndet))
          ! HKE -- commenting out for now
          !call read_hdf(chainfile, trim(adjustl(path))//'adc_Q', Q)
          !deallocate(Q)
       end if
    end if

    call mpi_bcast(base,  size(base) , MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    call mpi_bcast(phase, size(phase), MPI_REAL, 0, self%comm, ierr)
    call mpi_bcast(self%adu_range, size(self%adu_range), MPI_INTEGER, 0, self%comm, ierr)

    do j = 1, self%ndet
       do i = 1, self%nscan
          k = self%scanid(i)
          if (.not. self%scans(i)%d(j)%accept) cycle
          self%scans(i)%d(j)%baseline1 = base(k,j,1) 
          self%scans(i)%d(j)%baseline2 = base(k,j,2) 
          self%mod_phase(j,i)          = phase(k,j)
       end do
    end do

!!$    if (self%adu_range(1,1) > 0) then
!!$       if (self%myid /= 0) allocate(Q(minval(self%adu_range(:,1)):maxval(self%adu_range(:,2)),self%ndet))
!!$       call mpi_bcast(Q, size(Q), MPI_REAL, 0, self%comm, ierr)
!!$       do j = 1, self%ndet
!!$          self%adc(j)%p => comm_adc_binfit(self%comm, self%label(j), 16, &
!!$                  & self%adu_range(j,1), self%adu_range(j,2), 40)
!!$          !allocate(self%adc(j)%p%Q(self%adu_range(j,1):self%adu_range(j,2)))
!!$          self%adc(j)%p%Q = Q(self%adu_range(j,1):self%adu_range(j,2),j)
!!$       end do
!!$       deallocate(Q)
!!$    end if
    
    deallocate(base, phase)

  end subroutine initHDF_HFI
  
  module subroutine dumpToHDF_hfi(self, chainfile, path)
    ! 
    ! Writes instrument-specific TOD parameters to existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
    ! chainfile: derived type (hdf_file)
    !           Already open HDF file handle to existing chainfile
    ! path:   string
    !           HDF path to current dataset, e.g., "000001/tod/030"
    !
    ! Returns
    ! ----------
    ! None
    !
    implicit none
    class(comm_hfi_tod),                 intent(in)     :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path

    integer(i4b) :: ierr, i, j, k
    real(dp), allocatable, dimension(:,:,:) :: base
    real(sp), allocatable, dimension(:,:)   :: phase
    real(sp), allocatable, dimension(:,:)   :: Q

    allocate(base(self%nscan_tot,self%ndet,2), phase(self%nscan_tot,self%ndet))
    base  = 0.d0
    phase = 0.0
    do j = 1, self%ndet
       do i = 1, self%nscan
          k = self%scanid(i)
          if (.not. self%scans(i)%d(j)%accept) cycle
          base(k,j,1) = self%scans(i)%d(j)%baseline1
          base(k,j,2) = self%scans(i)%d(j)%baseline2
          phase(k,j)  = self%mod_phase(j,i)
       end do
    end do
    if (self%myid == 0) then
       call mpi_reduce(mpi_in_place, base, size(base), &
            & MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
       call write_hdf(chainfile, trim(adjustl(path))//'baseline_adc', base)
       call mpi_reduce(mpi_in_place, phase, size(phase), &
            & MPI_REAL, MPI_SUM, 0, self%comm, ierr)
       call write_hdf(chainfile, trim(adjustl(path))//'mod_phase', phase)
       call write_hdf(chainfile, trim(adjustl(path))//'adu_range', self%adu_range)

       ! Pad ADC tables
       if (associated(self%adc(1)%p)) then
          allocate(Q(minval(self%adu_range(:,1)):maxval(self%adu_range(:,2)),self%ndet))
          Q = 1.
          do j = 1, self%ndet
             Q(self%adu_range(j,1):self%adu_range(j,2),j) = self%adc(j)%p%Q
          end do
          call write_hdf(chainfile, trim(adjustl(path))//'adc_Q', Q)
          deallocate(Q)
       end if
    else
       call mpi_reduce(base,         base, size(base), &
            & MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
       call mpi_reduce(phase,        phase, size(phase), &
            & MPI_REAL, MPI_SUM, 0, self%comm, ierr)
    end if
    deallocate(base, phase)

  end subroutine dumpToHDF_hfi


  module subroutine construct_corrtemp_hfi(self, sd, det)
    !  Construct an HFI instrument-specific correction template; for now contains 1Hz template only
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  scan: int
    !       scan number
    !  pix: int
    !       index for pixel
    !  psi: int
    !       integer label for polarization angle
    !
    !  Returns:
    !  --------
    !  s:   real (sp)
    !       output template timestream
    implicit none
    class(comm_hfi_tod),  intent(in)           :: self
    class(comm_scandata), intent(inout)        :: sd
    integer(i4b),         intent(in), optional :: det

    integer(i4b) :: i, j, k, nbin, b
    real(dp)     :: dt, t_tot, t

    sd%s_inst = 0.

  end subroutine construct_corrtemp_hfi

  module subroutine apply_nonlin_corr_hfi(self, sd, nonlin_lvl, handle, det)
    !  Construct and apply HFI instrument-specific non-linear corrections
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  scan: int
    !       scan number
    !
    !  sd: comm_scandata
    !       structure holding the scan data
    !
    !  Returns:
    !  --------
    !  s:   real (sp)
    !       output template timestream
    implicit none
    class(comm_hfi_tod),                   intent(inout) :: self
    class(comm_scandata),                  intent(inout) :: sd
    integer(i4b),                          intent(in)    :: nonlin_lvl
    type(planck_rng),            optional, intent(inout) :: handle
    integer(i4b),                optional, intent(in)    :: det

    integer(i4b) :: i, j, d, scan, dec_wn

    scan = sd%scan

    ! Apply ADC corrections to raw self%tod_dp
    if (nonlin_lvl > 0 .and. associated(self%adc(1)%p)) then
       do i = 1, self%ndet
          d = i; if (present(det)) d = det
          if (.not. self%scans(scan)%d(d)%accept) cycle
          if (self%myid == 0 .and. sd%scan == 1 .and. i == 1) then
             write(*,*) trim(self%outdir)//'/adc_before.dat'
             open(58,file=trim(self%outdir)//'/adc_before.dat')
             do j = 1, sd%ntod
                write(58,*) j, sd%tod_dp(j,i), and(sd%flag(j,i), self%flag0)
             end do
             close(58)
          end if
          !call self%adc(d)%p%Q2As(92.)
          !call self%adc(d)%p%As2F
          call self%adc(d)%p%adc_correct(sd%tod_dp(:,i), &
               & iand(sd%flag(:,i), self%flag0) .eq.0, self%mod_phase(i,scan))
          if (self%myid == 0 .and. sd%scan == 1 .and. i == 1) then
             open(58,file=trim(self%outdir)//'/adc_after.dat')
             do j = 1, sd%ntod
                write(58,*) j, sd%tod_dp(j,i), and(sd%flag(j,i), self%flag0)
             end do
             close(58)
          end if
       end do
    end if

    ! Demodulate TOD; convert from double precision (sd%tod_dp) to single precision (sd%tod)
    if (nonlin_lvl > 1) then
       !if (present(handle)) call sample_hfi_baselines(sd, self, scan, handle)
       call demodulate_tod(sd, self, scan)
    end if
    

    ! 4K lines removal and/or in-painting
    if (nonlin_lvl > 2) then
       if (trim(self%cooler_correction_domain) == 'time') then
          do i = 1, self%ndet
             ! Correct 4K lines from TOD
             if (self%scans(scan)%d(i)%accept) call remove_hfi_4K_lines_tod(self, sd, i)
             ! In-paint flagged samples with s_tot + white noise
             d = i; if (present(det)) d = det
             if (.not. self%scans(scan)%d(d)%accept) cycle
             do j = 1, sd%ntod
                if (iand(sd%flag(j,i), self%flag0) .ne. 0) then
                   sd%tod(j,i) = self%scans(scan)%d(d)%gain * sd%s_tot(j,i,0,1) + self%scans(scan)%d(d)%N_psd%sigma0_preproc * rand_gauss(handle)
                end if
             end do
          end do
       elseif (trim(self%cooler_correction_domain) == 'frequency') then
          ! In-paint flagged samples with s_tot + white noise
          do i = 1, self%ndet
             d = i; if (present(det)) d = det
             if (self%scans(scan)%d(d)%accept) then
                do j = 1, sd%ntod
                   if (iand(sd%flag(j,i), self%flag0) .ne. 0) then
                      sd%tod(j,i) = self%scans(scan)%d(d)%gain * sd%s_tot(j,i,0,1) + self%scans(scan)%d(d)%N_psd%sigma0_preproc * rand_gauss(handle)
                   end if
                end do
             end if
             ! Correct 4K lines from PSD
             if (self%scans(scan)%d(i)%accept) call remove_hfi_4K_lines_harm(self, sd, i)
          end do
       end if
    end if

    ! Define high-frequency roll-off ! MODIFIED
    if (.false. .and. nonlin_lvl > 3) then
       do i = 1, self%ndet
          if (.not. self%scans(scan)%d(i)%accept) cycle
          call sd%spline_rolloff(i, self%scans(scan)%d(i)%gain, abs(self%scans(scan)%d(i)%N_psd%sigma0), self%samprate)
       end do
    end if

    ! Deconvolve regularized bolometer transfer function
    if (nonlin_lvl > 3) then
       do i = 1, self%ndet
          if (.not. self%scans(scan)%d(i)%accept) cycle
          if (self%myid == 0 .and. sd%scan == 1 .and. i == 1) then
             open(58,file='tod_before.dat')
             do j = 1, sd%ntod
                write(58,*) j, sd%tod(j,i)
             end do
             close(58)
          end if
          call self%Tbol(i)%p%convolve("preproc", "convolve", sd%tod(:,i), sd%plan_fwd, sd%plan_back)!, &
!                                       & sd%rolloff_filter(i), abs(self%scans(scan)%d(i)%N_psd%sigma0)) ! ADDED
          if (self%myid == 0 .and. sd%scan == 1 .and. i == 1) then
             open(58,file='tod_after.dat')
             do j = 1, sd%ntod
                write(58,*) j, sd%tod(j,i)
             end do
             close(58)
          end if
       end do
    end if

   ! At the end, update spline noise model nodes
   if (.true. .and. self%noise_psd_model == 'spline') then
      do i = 1, self%ndet
         if (self%scans(scan)%d(i)%accept) then
            call whitening_filter(self, sd, i)
            call update_spline_noise_psd(self,sd,scan,i)
         end if
      end do
   end if

   ! Estimate post-deconvolution white noise rms
   dec_wn = 2
   call sample_noise_psd(self, sd, handle, only_sigma0=.true., dec_wn=dec_wn)

  end subroutine apply_nonlin_corr_hfi

  module subroutine stitch_hfi_dc_level(self, scan, sd)
    !  Construct and apply HFI instrument-specific non-linear corrections
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  scan: int
    !       scan number
    !  sd: comm_scandata object
    !       structure holding the data for each scan
    !
    implicit none
    class(comm_hfi_tod),                   intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    class(comm_scandata),                  intent(inout) :: sd

    !estimate dc levels between jumps and then update sd%tod to be flat

  end subroutine stitch_hfi_dc_level

  module subroutine hfi_dark_correction(self, scan, sd)
    !  Construct and apply HFI instrument-specific corrections
    !  from dark bolometer timestreams
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  scan: int
    !       scan number
    !  sd: comm_scandata object
    !       structure holding the data for each scan
    !
    implicit none
    class(comm_hfi_tod),                   intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    class(comm_scandata),                  intent(inout) :: sd

    !estimate dark correction then update sd%tod to be flat

  end subroutine hfi_dark_correction

module subroutine estimate_hfi_4K_lines_tod(self, sd, i_det, output)
    !  Construct and apply HFI instrument-specific corrections
    !  from 4K lines
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  sd: comm_scandata object
    !       scan data
    !  i_det: int
    !       detector id
    !  output: lgt
    !          output TOD    
    implicit none
    class(comm_hfi_tod),               intent(inout) :: self
    class(comm_scandata),              intent(inout) :: sd
    integer(i4b),                      intent(in)    :: i_det
    logical(lgt),            optional, intent(in)    :: output

    integer(i4b) :: i, j, k, l, m, ntod, scan
    integer(i4b) :: n, nphi, dt_chunk, n_chunks, ind, i0, i1
    character(len=6) :: scantext, dettext
    logical(lgt) :: output_
    real(sp)     :: samprate, gain!, fade_width
    real(sp),     allocatable, dimension(:) :: d_prime, chunk, templ, template
    real(dp),     allocatable, dimension(:) :: times

    scan     = sd%scan
    ntod     = self%scans(scan)%ntod
    samprate = self%samprate
    n        = self%n_4K_lines
    nphi     = n+1

    output_  = .false.; if (present(output) .and. mod(self%scanid(scan),5000)==0) output_ = output
    if (output_) then
       call int2string(self%scanid(scan), scantext)
       call int2string(i_det, dettext)
    end if

    dt_chunk = 5*60*samprate ! 5 minutes chunks
    n_chunks = max(1, ntod/dt_chunk) ! Avoid small chunk at the end by defining a larger last scan
    !fade_width = 10 ! [s], time width of fading edges between chunks

    self%cooler_4K_lines_profile(i_det,scan)%n_chunks   = n_chunks
    !self%cooler_4K_lines_profile(i_det,scan)%fade_width = fade_width
    if (.not. allocated(self%cooler_4K_lines_profile(i_det,scan)%coeffs)) then
       allocate(self%cooler_4K_lines_profile(i_det,scan)%i0(n_chunks))
       allocate(self%cooler_4K_lines_profile(i_det,scan)%i1(n_chunks))
       allocate(self%cooler_4K_lines_profile(i_det,scan)%coeffs(n_chunks,2*nphi-2))
       self%cooler_4K_lines_profile(i_det,scan)%coeffs = 0.d0
    end if

    gain = self%scans(scan)%d(i_det)%gain
    allocate(d_prime(ntod), template(ntod))
    template = 0.d0
    d_prime = sd%tod(:,i_det) - gain * sd%s_tot(:,i_det,0,1)

    ! Output starting res tod
    if (output_) then
       open(58,file='res_tod_before_4K_removal_' // scantext // '_' // dettext // '.dat', recl=1024)
       do l = 1, ntod
          write(58,*) l, d_prime(l), sd%mask(l,i_det)
       end do
       close(58)
    end if

    d_prime = d_prime * sd%mask(:,i_det)
    do ind = 1, n_chunks
       i0 = (ind-1)*dt_chunk + 1
       i1 = ind*dt_chunk
       if (ind==n_chunks) i1 = ntod
       m = i1 - i0 + 1

       allocate(times(m), chunk(m), templ(m))
       templ = 0.d0
       chunk = d_prime(i0:i1)
       do j = i0, i1
          times(j-i0+1) = real(j-1,dp)/real(samprate,dp)
       end do
       times = times - times(m/2)

       call template_chunk(templ,times, chunk, sd%mask(i0:i1,i_det), ind)
       template(i0:i1) = templ
       deallocate(times,chunk,templ)
       self%cooler_4K_lines_profile(i_det,scan)%i0(ind) = i0
       self%cooler_4K_lines_profile(i_det,scan)%i1(ind) = i1
    end do


    ! Output corrected res tod
    if (output_) then
       d_prime = d_prime - template
       open(58,file='res_tod_after_4K_removal_' // scantext // '_' // dettext // '.dat', recl=1024)
       do l = 1, ntod
          write(58,*) l, d_prime(l), sd%mask(l,i_det), template(l)
       end do
       close(58)
    end if

    sd%tod(:,i_det) = sd%tod(:,i_det) - template
    deallocate(d_prime, template)

    contains

    subroutine template_chunk(templ,times,chunk,mask,ind)
       real(sp), dimension(:), intent(inout) :: templ
       real(sp), dimension(:), intent(in)    :: chunk, mask
       real(dp), dimension(:), intent(in)    :: times
       integer(i4b),           intent(in)    :: ind

       integer(i4b) :: ierr
       real(dp)     :: phi, tr
       real(dp), allocatable, dimension(:)   :: ATd, row, y, x
       real(dp), allocatable, dimension(:,:) :: ATA

       ! Check if enough unmasked samples
       if (count(mask > 0.5) < 10*2*nphi) then
          self%cooler_4K_lines_profile(i_det,scan)%coeffs(ind,:) = 0.d0
          return
       end if

       allocate(ATA(2*nphi, 2*nphi), ATd(2*nphi), row(2*nphi))
       allocate(y(2*nphi), x(2*nphi))
       ATA = 0.d0; ATd = 0.d0; row = 0.d0
       y = 0.d0; x = 0.d0

       do i = 1, m
          if (mask(i) == 0) cycle
          row(1) = 1.d0
          row(2) = times(i)
          do j = 1, n
             phi = 2*PI*self%nus_4K_lines(j)
             row(2*j+1) = cos(phi*times(i))
             row(2*j+2) = sin(phi*times(i))
          end do

          do j = 1, 2*nphi
             ATd(j) = ATd(j) + row(j)*real(chunk(i),dp)
             do k = j, 2*nphi
                ATA(j,k) = ATA(j,k) + row(j) * row(k)
             end do
          end do
       end do

       ! Symmetrization
       do i = 1, 2*nphi
          do j = i+1, 2*nphi
             ATA(j,i) = ATA(i,j)
          end do
       end do

       tr = 0.d0 ! regularization
       do i = 1, 2*nphi
          tr = tr + ATA(i,i)
       end do
       do i = 1, 2*nphi
          ATA(i,i) = ATA(i,i) + 1.d-12 * tr / real(2*nphi, dp)
       end do

       ! Cholesky
       call cholesky_decompose_single(ATA, ierr=ierr)
       if (ierr /= 0) then
          deallocate(ATA, ATd, row, y, x)
          write(*,*) 'Error: 4K cholesky decomposition failed - scan ', self%scanid(scan)
          self%cooler_4K_lines_profile(i_det,scan)%coeffs(ind,:) = 0.
          return
       end if

       do i = 1, 2*nphi
          y(i) = (ATd(i) - dot_product(ATA(i,1:i-1), y(1:i-1))) / ATA(i,i)
       end do
       do i = 2*nphi, 1, -1
          x(i) = (y(i) - dot_product(ATA(i+1:2*nphi,i), x(i+1:2*nphi))) / ATA(i,i)
       end do

       ! Build template
       do i = 1, n
          phi = 2*PI*self%nus_4K_lines(i)
          templ = templ + x(2*i+1) * cos(phi*times) + x(2*i+2) * sin(phi*times)
       end do

       self%cooler_4K_lines_profile(i_det,scan)%coeffs(ind,:) = x(3:2*nphi)

       deallocate(ATA,ATd,row,y,x)
    end subroutine template_chunk

  end subroutine estimate_hfi_4K_lines_tod

  module subroutine remove_hfi_4K_lines_tod(self, sd, i_det)
    implicit none
    class(comm_hfi_tod),               intent(inout) :: self
    class(comm_scandata),              intent(inout) :: sd
    integer(i4b),                      intent(in)    :: i_det

    integer(i4b) :: i, j, n, ntod, scan
    integer(i4b) :: i0, i1, nphi, n_chunks, ind, m
    real(sp) :: samprate 
    real(dp) :: phi
    logical(lgt) :: remove_4K_lines
    real(sp),     allocatable, dimension(:) :: templ, template, x
    real(dp),     allocatable, dimension(:) :: times

    scan     = sd%scan
    remove_4K_lines = .false.
    if (allocated(self%cooler_4K_lines_profile(i_det,scan)%coeffs)) then
       remove_4K_lines = .true.
    end if
    if (.not. remove_4K_lines) return

    ntod = self%scans(scan)%ntod
    samprate = self%samprate
    n        = self%n_4K_lines
    nphi     = n+1

    allocate(template(ntod), x(2*nphi-2))
    template=  0.d0
    n_chunks = self%cooler_4K_lines_profile(i_det,scan)%n_chunks

    do ind = 1, n_chunks
       i0 = self%cooler_4K_lines_profile(i_det,scan)%i0(ind)
       i1 = self%cooler_4K_lines_profile(i_det,scan)%i1(ind)
       m = i1 - i0 + 1

       allocate(times(m), templ(m))
       templ = 0.d0
       do j = i0, i1
          times(j-i0+1) = real(j-1,dp)/real(samprate,dp)
       end do
       times = times - times(m/2)

       ! Build template
       x = self%cooler_4K_lines_profile(i_det,scan)%coeffs(ind,:)
       do i = 1, n
          phi = 2*PI*self%nus_4K_lines(i)
          templ = templ + x(2*i-1) * cos(phi*times) + x(2*i) * sin(phi*times)
       end do
       template(i0:i1) = templ
       deallocate(times, templ)
    end do

    sd%tod(:,i_det) = sd%tod(:,i_det) - template
    deallocate(template, x)

  end subroutine remove_hfi_4K_lines_tod

  module subroutine estimate_hfi_4K_lines_harm(self, sd, i_det, apply_mask, output)
    !  Construct and apply HFI instrument-specific corrections
    !  from 4K lines
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  sd: comm_scandata object
    !       scan data
    !  i_det: int
    !       detector id
    !  apply_mask: logaical
    !              apply mask to residuals
    !  output: lgt
    !          output PSDs    
    implicit none
    class(comm_hfi_tod),               intent(inout) :: self
    class(comm_scandata),              intent(inout) :: sd
    integer(i4b),                      intent(in)    :: i_det
    logical(lgt),            optional, intent(in)    :: apply_mask
    logical(lgt),            optional, intent(in)    :: output

    integer(i4b) :: i, j, k, l, n, ntod, nomp, nfft, err, scan
    integer(i4b) :: i0, i1, nsub, maxiter, n_f0, n_sig
    integer*8    :: plan_fwd, plan_back
    logical(lgt) :: apply_mask_, output_
    character(len=6) :: scantext, dettext
    real(sp)     :: samprate, fmin, fmax, dnu, peak_val, gain, wn
    real(sp)     :: A_fit, f0_fit, sigma_fit
    real(sp),     allocatable, dimension(:)   :: dt, d_prime, profile, W
    complex(spc), allocatable, dimension(:)   :: dv, dv_4K
    real(sp),     allocatable, dimension(:,:) :: ps

    scan     = sd%scan
    ntod     = self%scans(scan)%ntod
    nomp     = 1
    samprate = self%samprate
    nfft     = 2 * ntod
    n        = nfft / 2 + 1
    wn       = abs(self%scans(scan)%d(i_det)%N_psd%sigma0_preproc)**2

    apply_mask_ = .false.; if (present(apply_mask)) apply_mask_ = apply_mask
    output_  = .false.; if (present(output) .and. mod(self%scanid(scan),5000)==0) output_ = output
    if (output_) then
       call int2string(self%scanid(scan), scantext)
       call int2string(i_det, dettext)
    end if


    call sfftw_init_threads(err)
    call sfftw_plan_with_nthreads(nomp)

    allocate(dt(nfft), dv(0:n-1), ps(1:n-1,2))
    call sfftw_plan_dft_r2c_1d(plan_fwd,  nfft, dt, dv, fftw_estimate + fftw_unaligned)
    call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)

    ! FFT
    gain = self%scans(scan)%d(i_det)%gain
    allocate(d_prime(ntod))
    d_prime = sd%tod(:,i_det) - gain * sd%s_tot(:,i_det,0,1)
    if (apply_mask_) d_prime = d_prime * sd%mask(:,i_det)
    dt(1:ntod)           = d_prime
    dt(2*ntod:ntod+1:-1) = dt(1:ntod)
    call timer%start(TOT_FFT)
    call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
    call timer%stop(TOT_FFT)
    do l = 1, n-1
       ps(l,1) = l*(samprate/2)/(n-1)
       ps(l,2) = abs(dv(l))** 2 / ntod
    end do
    deallocate(d_prime)
    
    ! Output starting noise power spectrum
    if (output_) then
       open(58,file='res_ps_before_4K_removal_' // scantext // '_' // dettext // '.dat', recl=1024)
       do l = 1, n-1
          write(58,*) ps(l,1), ps(l,2)
       end do
       close(58)
    
       allocate(W(1:n-1))
       W = 1.d0
    end if

    dnu = 0.03 ! Hz (size of freq window around each spike)
    maxiter = 1 ! Can be extended to iterative correction
    n_f0 = 20; n_sig = 10
    do i = 1, self%n_4K_lines
       fmin = self%nus_4K_lines(i) - dnu
       fmax = self%nus_4K_lines(i) + dnu

       ! Find window
       i0 = 1
       do while (i0 < ntod .and. ps(i0,1) < fmin)
          i0 = i0 + 1
       end do
       i1 = i0 + 1
       do while (i1 < ntod .and. ps(i1,1) < fmax)
          i1 = i1 + 1
       end do
       nsub = i1 - i0 + 1
       if(nsub < 10) cycle ! too small window
       
       if (allocated(self%cooler_4K_lines(i,i_det,scan)%p%spike_profile)) then
          deallocate(self%cooler_4K_lines(i,i_det,scan)%p%spike_profile)
       end if
       
       self%cooler_4K_lines(i,i_det,scan)%p%i0 = i0
       self%cooler_4K_lines(i,i_det,scan)%p%window = nsub
       allocate(self%cooler_4K_lines(i,i_det,scan)%p%spike_profile(nsub))
       allocate(profile(nsub))
       profile = 0.d0
       self%cooler_4K_lines(i,i_det,scan)%p%spike_profile = 0.d0

       if (.true.) then
          ! Gaussian fit
          ! Iterative cleaning loop
          if (.not. allocated(self%cooler_4K_lines(i,i_det,scan)%p%A_fit)) allocate(self%cooler_4K_lines(i,i_det,scan)%p%A_fit(maxiter))
          if (.not. allocated(self%cooler_4K_lines(i,i_det,scan)%p%f0_fit)) allocate(self%cooler_4K_lines(i,i_det,scan)%p%f0_fit(maxiter))
          if (.not. allocated(self%cooler_4K_lines(i,i_det,scan)%p%sigma_fit)) allocate(self%cooler_4K_lines(i,i_det,scan)%p%sigma_fit(maxiter))
          do j = 1, maxiter
             ! 1. baseline_estimation
             if (allocated(self%cooler_4K_lines(i,i_det,scan)%p%baseline)) then 
                deallocate(self%cooler_4K_lines(i,i_det,scan)%p%baseline)
             end if
             call self%cooler_4K_lines(i,i_det,scan)%p%estimate_4K_baseline(ps(i0:i1,:))

             ! 2. residual = sub_ps - baseline
             ps(i0:i1,2) = ps(i0:i1,2) - self%cooler_4K_lines(i,i_det,scan)%p%baseline

             ! 3. amplitude linear fit
             call self%cooler_4K_lines(i,i_det,scan)%p%A_lin_fit(ps(i0:i1,:),j,n_f0,n_sig)
             A_fit = self%cooler_4K_lines(i,i_det,scan)%p%A_fit(j)
             f0_fit = self%cooler_4K_lines(i,i_det,scan)%p%f0_fit(j)
             sigma_fit = self%cooler_4K_lines(i,i_det,scan)%p%sigma_fit(j)
           
             ! 4. subtract gaussian fit
             ps(i0:i1,2) = ps(i0:i1,2) + self%cooler_4K_lines(i,i_det,scan)%p%baseline
             do k = 1, nsub
                profile(k) = profile(k) + A_fit * exp(-0.5 * ((ps(i0+k-1,1) - f0_fit)/sigma_fit)**2)
             end do

             ! 5. tolerance check?
      
          end do
       else
          ! Delta correction
          do j = 1, maxiter
             peak_val = maxval(ps(i0:i1,2))
             do k = 1, nsub
                if (ps(i0+k-1,2)>=peak_val) then
                   profile(k) = profile(k) + peak_val
                   ps(i0+k-1,2) = wn
                end if
             end do
          end do
        end if


       ! Save line profile
       self%cooler_4K_lines(i,i_det,scan)%p%spike_profile = self%cooler_4K_lines(i,i_det,scan)%p%spike_profile + profile
       deallocate(profile)
       if (output_) W(i0:i1) =  self%cooler_4K_lines(i,i_det,scan)%p%baseline / &
                             & (self%cooler_4K_lines(i,i_det,scan)%p%baseline + self%cooler_4K_lines(i,i_det,scan)%p%spike_profile)
    end do


    if (output_) then
       allocate(dv_4K(0:n-1))
       dv_4K = dv
       dv_4K(1:n-1) = dv_4K(1:n-1) * sqrt(W)
       deallocate(W)

       ! Output ending noise power spectrum
       open(58,file='res_ps_after_4K_removal_' // scantext // '_' // dettext // '.dat', recl=1024)
       do l = 1, n-1
          write(58,*) ps(l,1), abs(dv_4K(l))** 2 / ntod
       end do
       close(58)

       ! Output 4K lines PSD and TOD
       dv_4K = dv - dv_4K

       open(58,file='4K_lines_ps_' // scantext // '_' // dettext // '.dat', recl=1024)
       do l = 1, n-1
          write(58,*) ps(l,1), abs(dv_4K(l))** 2 / ntod
       end do
       close(58)

       call timer%start(TOT_FFT)
       call sfftw_execute_dft_c2r(plan_back, dv_4K, dt)
       call timer%stop(TOT_FFT)
       deallocate(dv_4K)
       dt  = dt / nfft

       open(58,file='4K_lines_tod_' // scantext // '_' // dettext // '.dat', recl=1024)
       do l = 1, n-1
          write(58,*) l, dt(l)
       end do
       close(58)
    end if
    

    !! If in the low-level loop there is no other step after this, there's no need to correct the spectrum.
    !!  Otherwise, uncomment the following lines

    ! FFT total TOD
!    dt(1:ntod)           = sd%tod(:,i_det)
!    dt(2*ntod:ntod+1:-1) = dt(1:ntod)
!    call timer%start(TOT_FFT)
!    call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
!    call timer%stop(TOT_FFT)
!
!    dv(1:n-1) = dv(1:n-1) * sqrt(W)
!    deallocate(W)
!
!    ! FFT^-1 back to TOD
!    call timer%start(TOT_FFT)
!    call sfftw_execute_dft_c2r(plan_back, dv, dt)
!    call timer%stop(TOT_FFT)
!
!    dt  = dt / nfft
!    sd%tod(:,i_det) = dt(1:ntod)

    deallocate(dt, dv, ps)
    call sfftw_destroy_plan(plan_fwd)
    call sfftw_destroy_plan(plan_back)

  end subroutine estimate_hfi_4K_lines_harm

  module subroutine remove_hfi_4K_lines_harm(self, sd, i_det)
    !  Apply HFI instrument-specific corrections from 4K lines
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !  
    !  sd: comm_scandata object
    !      scan data
    !  i_det: int
    !       detector id
    implicit none
    class(comm_hfi_tod),               intent(inout) :: self
    class(comm_scandata),              intent(inout) :: sd
    integer(i4b),                      intent(in)    :: i_det

    integer(i4b) :: i, j, k, l, n, ntod, nomp, nfft, err, scan
    integer(i4b) :: i0, i1, nsub
    integer*8    :: plan_fwd, plan_back
    logical(lgt) :: remove_4K_lines
    real(sp)     :: samprate, fmin, fmax, dnu, peak_val, gain
    real(sp),     allocatable, dimension(:)   :: dt, W
    complex(spc), allocatable, dimension(:)   :: dv

    scan     = sd%scan
    remove_4K_lines = .false.
    do i = 1, self%n_4K_lines
       if (allocated(self%cooler_4K_lines(i,i_det,scan)%p%spike_profile)) then
          remove_4K_lines = .true.
       end if
    end do
    if (.not. remove_4K_lines) return


    ntod = self%scans(scan)%ntod
    nomp     = 1
    samprate = self%samprate
    nfft     = 2 * ntod
    n        = nfft / 2 + 1

    call sfftw_init_threads(err)
    call sfftw_plan_with_nthreads(nomp)

    allocate(dt(nfft), dv(0:n-1), W(1:n-1))
    call sfftw_plan_dft_r2c_1d(plan_fwd,  nfft, dt, dv, fftw_estimate + fftw_unaligned)
    call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)

    ! FFT
    dt(1:ntod)           = sd%tod(:,i_det)
    dt(2*ntod:ntod+1:-1) = dt(1:ntod)
    call timer%start(TOT_FFT)
    call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
    call timer%stop(TOT_FFT)

    W = 1.d0
    do i = 1, self%n_4K_lines
       i0 = self%cooler_4K_lines(i,i_det,scan)%p%i0
       nsub = self%cooler_4K_lines(i,i_det,scan)%p%window
       i1 = nsub + i0 - 1
       if(nsub < 10) cycle ! too small window
       if (.not. allocated(self%cooler_4K_lines(i,i_det,scan)%p%spike_profile)) cycle 

       W(i0:i1) =  self%cooler_4K_lines(i,i_det,scan)%p%baseline / &
                & (self%cooler_4K_lines(i,i_det,scan)%p%baseline + self%cooler_4K_lines(i,i_det,scan)%p%spike_profile)
    end do

    dv(1:n-1) = dv(1:n-1) * sqrt(W)
    deallocate(W)

    ! FFT^-1 back to TOD
    call timer%start(TOT_FFT)
    call sfftw_execute_dft_c2r(plan_back, dv, dt)
    call timer%stop(TOT_FFT)

    dt  = dt / nfft
    sd%tod(:,i_det) = dt(1:ntod)

    deallocate(dt, dv)
    call sfftw_destroy_plan(plan_fwd)
    call sfftw_destroy_plan(plan_back)

  end subroutine remove_hfi_4K_lines_harm


  module subroutine fill_gaps(self, tod, handle, scan, i_det, mask, s_sub, pix, nomono, dospike, ps_output, filling)
    implicit none
    class(comm_hfi_tod),                      intent(in)    :: self
    real(sp),               dimension(1:),    intent(inout) :: tod
    type(planck_rng),                         intent(inout) :: handle
    integer(i4b),                             intent(in)    :: scan, i_det
    real(sp),               dimension(1:),    intent(in)    :: mask, s_sub
    integer(i4b), optional, dimension(1:,1:), intent(in)    :: pix
    logical(lgt),                   optional, intent(in)    :: nomono
    logical(lgt),                   optional, intent(in)    :: dospike
    character(len=*),               optional, intent(in)    :: ps_output
    character(len=*),               optional, intent(in)    :: filling

    integer(i4b) :: i, j, k, l, n, ntod, nomp, nfft, err
    integer(i4b) :: j_end, j_start
    character(len=12) :: filling_
    integer*8    :: plan_fwd, plan_back
    logical(lgt) :: init_masked_region, end_masked_region, nomono_
    real(sp)     :: sigma_0, gain, N_wn, samprate
    real(sp),     allocatable, dimension(:)   :: d_prime, dt
    complex(spc), allocatable, dimension(:)   :: dv
    real(sp),     allocatable, dimension(:,:) :: ps

    integer(i4b) :: n_good, n_bad, good_count, bad_count, chunk_size
    integer(i4b) :: region_lgt ! logical
    real(sp)     :: rand_id
    integer(i4b), allocatable, dimension(:) :: good_id
    integer(i4b), allocatable, dimension(:,:) :: good_info, bad_info  ! 3-d arrays: (j_start,j_end,size)

    nomono_ = .false.; if (present(nomono)) nomono_ = nomono
    filling_ = 'white'; if (present(filling)) filling_ = filling

    ntod = self%scans(scan)%ntod
    allocate(d_prime(ntod))
    gain     = self%scans(scan)%d(i_det)%gain  ! Gain in V / K
    sigma_0  = abs(self%scans(scan)%d(i_det)%N_psd%sigma0)
    N_wn     = sigma_0**2  ! white noise power spectrum

    ! Prepare TOD residual
    d_prime = tod - gain * s_sub

    if (trim(filling_)=='white') then ! filling with white noise
       init_masked_region = .true.
       end_masked_region  = .false.
       do j = 1, ntod
          if (mask(j) == 1.) then
             if (end_masked_region) then
                j_end = j - 1
                call fill_masked_region(d_prime, mask, j_start, j_end, ntod, self%scans(scan)%chunk_num)
                ! Add noise to masked region
                if (trim(self%operation) == "sample") then
                   do k = j_start, j_end
                      d_prime(k) = d_prime(k) + sigma_0 * rand_gauss(handle)
                   end do
                end if
                end_masked_region = .false.
                init_masked_region = .true.
             end if
          else
             if (init_masked_region) then
                init_masked_region = .false.
                end_masked_region = .true.
                j_start = j
             end if
          end if
       end do
       ! if the data ends with a masked region
       if (end_masked_region) then
          j_end = ntod
          call fill_masked_region(d_prime, mask, j_start, j_end, ntod, self%scans(scan)%chunk_num)
          if (trim(self%operation) == "sample") then
             do k = j_start, j_end
                d_prime(k) = d_prime(k) + sigma_0 * rand_gauss(handle)
             end do
          end if
       end if
    else if (trim(filling_)=='chunks') then ! filling with good chunks
       ! Separate flagged and unflagged regions
       n_good = 0; n_bad = 0
       do i = 1, ntod
          if (i==1) then
             if (mask(i) == 0) then
                n_bad = n_bad + 1
             else
                n_good = n_good + 1
             end if
             region_lgt = mask(i)
          else
             if (mask(i) /= region_lgt) then
                if (mask(i) == 0) then
                   n_bad = n_bad + 1
                else
                   n_good = n_good + 1
                end if
                region_lgt = mask(i)
             end if
          end if
       end do

       allocate(good_info(3,n_good), bad_info(3,n_bad))
       good_count = 1; bad_count = 1
       do i = 1, ntod
          if (i==1) then
             j_start = 1
             region_lgt = mask(i)
          else
             if (mask(i) /= region_lgt) then
                j_end = i - 1
                n = j_end - j_start + 1
                if (mask(i) == 0) then
                   good_info(1,good_count) = j_start
                   good_info(2,good_count) = j_end
                   good_info(3,good_count) = n
                   good_count = good_count + 1
                else
                   bad_info(1,bad_count) = j_start
                   bad_info(2,bad_count) = j_end
                   bad_info(3,bad_count) = n
                   bad_count = bad_count + 1
                end if
                j_start = i
                region_lgt = mask(i)
             end if
          end if
       end do

       j_end = ntod
       n = j_end - j_start + 1
       if (mask(ntod) == 0) then
          bad_info(1,bad_count) = j_start
          bad_info(2,bad_count) = j_end
          bad_info(3,bad_count) = n
       else
          good_info(1,good_count) = j_start
          good_info(2,good_count) = j_end
          good_info(3,good_count) = n
       end if


       ! Fill bad data with random chunks from good data
       do i = 1, n_bad
          chunk_size = bad_info(3,i)
          ! Find good regions longer than chunk size
          k = 0
          do j = 1, n_good
             if (good_info(3,j) >= chunk_size) k = k + 1
          end do

          do while (k==0) ! Flagged region is too big
             call random_number(rand_id)
             j = rand_id*n_good + 1
             d_prime(bad_info(1,i):bad_info(1,i)+good_info(3,j)-1) = d_prime(good_info(1,j):good_info(2,j))
             bad_info(1,i) = bad_info(1,i) + good_info(3,j)
             bad_info(2,i) = bad_info(2,i) + good_info(3,j)
             bad_info(3,i) = bad_info(3,1) - good_info(3,j)
             chunk_size = bad_info(3,i)
             do j = 1, n_good
                if (good_info(3,j) >= chunk_size) k = k + 1
             end do
          end do

          allocate(good_id(k))
          k = 0
          do j = 1, n_good
             if (good_info(3,j) >= chunk_size) then
                k = k + 1
                good_id(k) = j
             end if
          end do
          ! Choose chunk from random good region
          call random_number(rand_id)
          k = rand_id*k + 1
          j = good_id(k)
          !write(*,*) 'CPU ',self%scanid(scan),': ', i, ntod, bad_info(1,i), bad_info(2,i), size(good_id), j, size(good_info)
          deallocate(good_id)
          call random_number(rand_id)
          j_start = rand_id*(good_info(3,j) - chunk_size) + good_info(1,j)
          d_prime(bad_info(1,i):bad_info(2,i)) = d_prime(j_start:j_start+chunk_size-1)
       end do
       deallocate(good_info, bad_info)
    else if (trim(filling_)=='zero') then ! filling with 0
       do i = 1, ntod
          if (mask(i)==0) d_prime(i) = 0.d0
       end do
    end if

    ! Identify spikes
    !if (self%first_call .and. .not. (present(dospike))) call find_d_prime_spikes(self, scan, i_det, d_prime, pix)

    !alpha    = self%scans(scan)%d(i_det)%N_psd%alpha
    !nu_knee  = self%scans(scan)%d(i_det)%N_psd%fknee

    ! Remove monopole if requested by user
    if (nomono_) d_prime = d_prime -  sum(d_prime*mask)/sum(mask)

    ! Output power spectrum of the residuals
    if (present(ps_output)) then
!       open(58,file='d_prime_tod_' // ps_output // '.dat', recl=1024)
!       do l = 1, ntod
!          write(58,*) l, tod(l), s_sub(l), d_prime(l), gain, sigma_0
!       end do
!       close(58)
       
       nomp     = 1 !omp_get_max_threads()
       samprate = self%samprate
       nfft     = 2 * ntod
       n        = nfft / 2 + 1

       call sfftw_init_threads(err)
       call sfftw_plan_with_nthreads(nomp)

       allocate(dt(nfft), dv(0:n-1), ps(1:n-1,2))
       call sfftw_plan_dft_r2c_1d(plan_fwd,  nfft, dt, dv, fftw_estimate + fftw_unaligned)
       call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)

       ! FFT
       dt(1:ntod)           = d_prime(:)
       dt(2*ntod:ntod+1:-1) = dt(1:ntod)
       call timer%start(TOT_FFT)
       call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
       call timer%stop(TOT_FFT)
       do l = 1, n-1
          ps(l,1) = l*(samprate/2)/(n-1)
          ps(l,2) = abs(dv(l)) ** 2 / ntod
       end do

       ! Output power spectrum of signal-subtracted gap-filled TOD to disk
       open(58,file='d_prime_ps_' // ps_output // '.dat', recl=1024)
       do l = 1, n-1
          write(58,*) ps(l,1), ps(l,2)
       end do
       close(58)
       

       deallocate(dt, dv, ps)
       call sfftw_destroy_plan(plan_fwd)
       call sfftw_destroy_plan(plan_back)
    end if

    tod = d_prime + gain * s_sub
    deallocate(d_prime)

  end subroutine fill_gaps


  module subroutine sample_adc_and_baselines(self, handle, det)
    !  Sample ADC parameters
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    implicit none
    class(comm_hfi_tod),                       intent(inout) :: self
    type(planck_rng),                          intent(inout) :: handle
    integer(i4b),                              intent(in)    :: det

    integer(i4b)        :: i, j, k, flag, ierr, ntod, decimation, offset, ind, oper
    integer(i8b)        :: ntot
    real(sp)            :: gain, phase, base1, base2
    real(dp)            :: chisq
    type(comm_scandata) :: sd
    character(len=4)    :: id
    integer(i8b), allocatable, dimension(:)   :: numsamp
    real(sp),     allocatable, dimension(:)   :: s_volt, tod, tod_corr
    real(sp),     allocatable, dimension(:)   :: sigma0

    decimation = 45 ! Downsampling rate; MUST BE ODD, to explore both parities

    if (self%correct_sl) then
       oper = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
            & SD_SKY,SD_BP,SD_ORB,SD_INST,SD_DARK,SD_SL])
    else
       oper = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
            & SD_SKY,SD_BP,SD_ORB,SD_INST,SD_DARK])
    end if
    !oper = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
    !     & SD_SKY,SD_BP,SD_ORB,SD_INST,SD_DARK])
    
    ! Count number of unmasked and decimated samples
    allocate(numsamp(self%nscan), sigma0(self%nscan))
    do i = 1, self%nscan
       numsamp(i) = self%scans(i)%ntod / decimation
       sigma0(i)  = self%scans(i)%d(det)%N_psd%sigma0
     end do
    ntot = sum(numsamp)

    ! Prepare reduced dataset
    allocate(s_volt(ntot), tod(ntot), tod_corr(maxval(numsamp)))
    ind = 1
    do i = 1, self%nscan
       if (.not. self%scans(i)%d(det)%accept) then
          numsamp(i) = ind-1
          cycle
       end if
       call init_scan_data(self, i, oper, TODMASK_PROC, sd, det=det, nonlin_level=0)

      gain    = self%scans(i)%d(det)%gain
      phase   = self%mod_phase(det,i)
      base1   = self%scans(i)%d(det)%baseline1
      base2   = self%scans(i)%d(det)%baseline2
      
       ! Decimate data
       j = 0; k = decimation*j+1
       do while (k <= sd%ntod .and. j < numsamp(i))
          if (sd%mask(k,1) == 1.) then
             tod(ind)   = sd%tod(k,1)
             if (mod(k,2) == 1) then
                s_volt(ind) =  phase*gain * sd%s_tot(k,1,0,1) + base1
             else
                s_volt(ind) = -phase*gain * sd%s_tot(k,1,0,1) + base2
             end if
             ind = ind+1
          end if
          j = j+1; k = decimation*j+1
       end do

       ! Store total number of accepted samples until now
       numsamp(i) = ind-1

       ! Clean up
       call dealloc_scan_data(sd)
    end do

    ! Compute number of accepted samples for each scan
    do i = self%nscan, 2, -1
       numsamp(i) = numsamp(i)-numsamp(i-1)
    end do
    ntot = sum(numsamp)

!!$    call int2string(self%myid, id)
!!$    open(58,file='adc_data_'//trim(self%adc(det)%p%label)//'_id'//id//'.dat', recl=1024)
!!$    do i = 1, ntot
!!$       write(58,*) i, tod(i), s_volt(i), tod(i)-s_volt(i)
!!$    end do
!!$    close(58)
    
    
    !write(*,*) self%myid, self%nscan, numsamp, ntot
    
    !  Perform sampling
    if (self%myid == 0) then
       if (.true.) then
          call self%adc(det)%p%powell_adc(powell_chisq_adc_hfi)
       else
          call self%adc(det)%p%mcmc_sample_adc(handle, chisq_adc_hfi)
       end if
       ! Release workers
       flag = 0
       call mpi_bcast(flag, 1, MPI_INTEGER, 0, self%comm, ierr)
    else
       do while (.true.)
          call mpi_bcast(flag, 1, MPI_INTEGER, 0, self%comm, ierr)
          if (flag > 0) then
             chisq = chisq_adc_hfi()
          else
             exit
          end if
       end do
    end if

    call mpi_bcast(self%adc(det)%p%p, self%adc(det)%p%npar_adc, MPI_REAL, 0, &
         & self%comm, ierr)    

    ! Update parameters
    call self%adc(det)%p%param2Q

    ! Clean up
    deallocate(s_volt, tod, sigma0, tod_corr)
    
  contains

    function powell_chisq_adc_hfi(x)
      implicit none
      real(dp), dimension(:), intent(in),  optional :: x
      real(dp)                                      :: powell_chisq_adc_hfi

      powell_chisq_adc_hfi = chisq_adc_hfi(x)

    end function powell_chisq_adc_hfi


    function chisq_adc_hfi(x, ndof) result (chisq)
      implicit none
      real(dp), dimension(:), intent(in),  optional :: x
      integer(i8b),           intent(out), optional :: ndof
      real(dp)                                      :: chisq

      integer(i4b) :: i, j, n
      integer(i8b) :: ndof_sub, ndof_tot, k1, k2
      real(dp)     :: chisq_sub, A, b
      real(sp)     :: baseline
      logical(lgt) :: output_ndof
      real(dp), allocatable, dimension(:) :: p

      allocate(p(self%adc(det)%p%npar_adc))
      if (self%myid == 0) then
         flag = 1; if (present(ndof)) flag = 2
         call mpi_bcast(flag, 1, MPI_INTEGER, 0, self%comm, ierr)
         p = x
      end if
      call mpi_bcast(p, size(p), MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
      output_ndof = (flag == 2)

      ! Check priors
      if (any(p < 0.0)) then
         chisq = 1.d30
         deallocate(p)
         return
      end if
      
      ! Update ADC parameters; assume constant sigma0 for now
      call self%adc(det)%p%param2Q(p)
      call self%adc(det)%p%Q2As(92.)
      call self%adc(det)%p%As2F

!!$      if (self%myid == 0 .and. det==1) then
!!$         call int2string(self%myid, id)
!!$         open(58,file='adc_reddata_'//trim(self%adc(det)%p%label)//'_id'//id//'.dat', recl=1024)
!!$      end if
      
      ! Evaluate chisq
      chisq_sub = 0.d0
      k2 = 0 
      do i = 1, self%nscan
         ! Skip scan if no accepted data
         if (.not. self%scans(i)%d(det)%accept) cycle
         n = numsamp(i); k1 = k2+1; k2 = k2+n
         
         ! Apply ADC to voltages
         do j = 1, n
            !tod_corr(j) = splint(self%adc(det)%p%F, real(s_volt(k1+j-1),dp))
            tod_corr(j) = self%adc(det)%p%invF(nint(tod(k1+j-1)))
         end do
         
!!$         if (self%myid == 0 .and. det==1) then
!!$            do j = 1, n
!!$               write(58,*) k1+j-1, tod(k1+j-1), tod_corr(j), s_volt(k1+j-1), tod_corr(j)-s_volt(k1+j-1), (tod_corr(j)-s_volt(k1+j-1))/sigma0(i), sigma0(i)
!!$            end do
!!$         end if
         
         ! Compute chisq
         do j = 1, n
            chisq_sub = chisq_sub + (tod_corr(j)-s_volt(k1+j-1))**2/sigma0(i)**2
         end do
      end do

!!$      if (self%myid == 0 .and. det==1) close(58)
      
      call mpi_reduce(chisq_sub, chisq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
      if (output_ndof) then
         call mpi_reduce(ntot, ndof_tot, 1, MPI_INTEGER8, MPI_SUM, 0, self%comm, ierr)
         if (self%myid == 0) ndof = ndof_tot
      end if

      !if (self%myid == 0) write(*,*) 'adc chisq =', chisq, ', p = ',  real(p,sp)

    end function chisq_adc_hfi

  end subroutine sample_adc_and_baselines
  
end submodule comm_tod_hfi_smod
