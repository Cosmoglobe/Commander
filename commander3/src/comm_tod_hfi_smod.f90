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
submodule (comm_tod_hfi_mod) comm_tod_hfi_mod
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

    ! Initialize common parameters
    call c%tod_constructor(cpar, id, id_abs, info, tod_type)
    
    ! Initialize instrument-specific parameters

    c%compressed_tod  = .true.
    c%correct_sl      = .false.
    c%correct_orb     = .true.
    c%apply_inst_corr = .true.
    c%orb_4pi_beam    = .false.
    c%symm_flags      = .true.
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
    c%n_xi            = 3
    !c%noise_psd_model = 'white'       ! Using white noise until we get better estimates of the actual noise PSD
    allocate(c%xi_n_P_uni(c%n_xi,2))
    allocate(c%xi_n_P_rms(c%n_xi))
    allocate(c%xi_n_nu_fit(c%n_xi,2))

    c%xi_n_P_uni(1,:)  = [10d0, 300d0]  ! Signa0
    c%xi_n_P_uni(2,:)  = [0.001d0, 20d0]  ! fknee
    c%xi_n_P_uni(3,:)  = [-2.5d0, -0.3d0]   ! alpha
    !c%xi_n_P_uni(4,:)  = [ 0.5d0,  4.0d0]  ! fknee
    !c%xi_n_P_uni(5,:)  = [-1.5d0, -0.5d0]   ! alpha
    c%xi_n_nu_fit(1,:) = [0.001d0, 80d0] 
    c%xi_n_nu_fit(2,:) = [0.001d0, 80d0]
    c%xi_n_nu_fit(3,:) = [0.001d0, 80d0]
    !c%xi_n_nu_fit(4,:) = [0.001d0, 10d0]
    !c%xi_n_nu_fit(5,:) = [0.001d0, 10d0]
    c%xi_n_P_rms       = [10.d0, 0.1d0, 0.1d0] ! [sigma0, fknee, alpha]; sigma0 is not used
    c%f_spin           = 1./60.                 ! Planck spin frequency in Hz
    

    !c%xi_n_P_rms      = [-1.d0] ! [sigma0]; sigma0 is not used
    
    ! Channel specific parameters
    if (c%freq(1:3) == "100") then
       c%chisq_threshold  = 120.d0
       c%sigma0_threshold = 100.d0
       c%accept_threshold = 0.8d0
       c%correct_sl       = .false.
    else if (c%freq(1:3) == "353") then
       c%chisq_threshold  = 1000.d0
       c%sigma0_threshold = 1000.d0
       c%accept_threshold = 0.8d0
       c%correct_sl       = .false.
    else
       c%chisq_threshold  = 100.d0 
       c%accept_threshold = 0.5d0
       c%correct_sl       = .false.
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
    c%pixcache => comm_tod_pixcache(c%nside, c%nside_beam, c%nmaps, .false.)
    call c%precompute_lookups()
    
    ! Load the instrument file
    call c%load_instrument_file(c%nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

    ! Collect Sun velocities from all scans
    call c%collect_v_sun
    
    ! Allocate sidelobe convolution data structures
    if(c%correct_sl) allocate(c%slconv(c%ndet,c%nhorn), c%orb_dp)
    
    allocate(c%orb_dp)
    !c%orb_dp => comm_orbdipole(c%mbeam)  ! HKE: Removed mbeam for now due to crash; should be fixed
    c%orb_dp => comm_orbdipole(comm=c%comm)

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
    allocate(correlations(c%ndet, c%ndet))
    correlations = .true.    

  !  allocate(c%xtalk)
    c%xtalk => comm_crosstalk(correlations)

    ! Pre-initialize ADC object
    allocate(c%adc(c%ndet))
    allocate(c%adu_range(c%ndet,2))

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
    integer(i4b)        :: i, j, k, h, l, ierr, ndelta, nside, npix, nmaps, dec_wn, oper_default
    logical(lgt)        :: select_data, output_scanlist, output_zodi_comps
    logical(lgt)        :: sample_gain, sample_ncorr, sample_abs_bandpass, sample_rel_bandpass, sample_zodi, sample_adc, make_dyn_mask, sample_xi_n
    type(comm_binmap)   :: binmap
    type(comm_scandata) :: sd
    !type(comm_detdata)  :: dd
    character(len=4)    :: ctext, myid_text
    character(len=6)    :: samptext, scantext, itertext
    character(len=512)  :: prefix, postfix, prefix4D, Sfilename
    character(len=512), allocatable, dimension(:) :: slist
    real(sp),              dimension(9)       :: flag_threshold
    real(sp), allocatable, dimension(:)       :: procmask, procmask2, procmask_zodi, sigma0, freqmask
    real(sp), allocatable, dimension(:,:)     :: s_buf
    real(sp), allocatable, dimension(:,:,:)   :: d_calib
    real(sp), allocatable, dimension(:,:,:,:) :: map_sky, m_gain
    real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf

    call int2string(iter, ctext)
    call update_status(status, "tod_start"//ctext)
    call timer%start(TOD_TOT, self%band)
    
    ! Toggle optional operations
    sample_rel_bandpass   = size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
    sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses
    if (.false.) then ! Debug
       ! Do data selection, then start sampling
       sample_gain           = iter  > 0 !.true.                 
       make_dyn_mask         = iter == 1
       sample_ncorr          = iter  > 0 !.true.
       sample_xi_n      = .false.
       select_data           = iter == 1
       sample_adc            = .false. !iter  > 1 !.true.
    else if (trim(self%init_from_HDF) == 'none') then
       ! Initialize slowly if not HDF init
       sample_gain           = iter  > 0 !.true.                 
       make_dyn_mask         = iter == 2
       sample_ncorr          = iter > 4 !.true.
       sample_xi_n           = iter > 5 
       select_data           = iter == 3 ! self%first_call  
       sample_adc            = .false. !iter  > 6 ! 3 !.true.
    else
       ! Do data selection, then start sampling
       sample_gain           = iter  > 1 !.true.                 
       make_dyn_mask         = iter == 1
       sample_ncorr          = iter  > 0 !.true.
       sample_xi_n           = iter > 1 !.false.
       select_data           = iter == 1 ! self%first_call  
       sample_adc            = iter  > 1 !.true.
    end if
    sample_zodi           = self%sample_zodi .and. self%subtract_zodi ! Sample zodi parameters
    output_zodi_comps     = self%output_zodi_comps .and. self%subtract_zodi ! Output zodi components
    output_scanlist       = mod(iter-1,10) == 0    ! only output scanlist every 10th iteration
    dec_wn                = 2 ! Decimation factor for sigma0; 2 corresponds to 45Hz

    !                       Pixhist    Single abs/RMS       RMS ranges     Single     Ranges   Pointing
    if (self%freq(1:3) == "100") then
       flag_threshold     = [  1.0,        20.0, 5.0,         1.5, 2.0, -1.0,   -1.0,      -1.0,     -1.0]
    else if (self%freq(1:3) == "353") then
       flag_threshold     = [  1.0,       -200.0, -50.0,    -1.5, -2.0, -1.0,   -1.0,      -1.0,     -1.0]
    else
       flag_threshold     = [  1.0,        20.0, 5.0,         1.5, 2.0, -1.0,   1.0,      -1.0,     -1.0]
    end if


    oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
         & SD_SKY,SD_BP,SD_ORB,SD_INST,SD_DARK,SD_NCORR])
    
    ! Initialize local variables
    ndelta          = size(delta,3)
    self%n_bp_prop  = ndelta-1
    nside           = map_out%info%nside
    nmaps           = map_out%info%nmaps
    npix            = 12*nside**2
    self%output_n_maps = 3
    if (self%output_aux_maps > 0 .or. .true.) then
       if (mod(iter-1,self%output_aux_maps) == 0) self%output_n_maps = 7
    end if
    if (output_zodi_comps) self%output_n_maps = self%output_n_maps + zodi_model%n_comps

    call int2string(chain, ctext)
    call int2string(iter, samptext)
    call int2string(self%myid, myid_text)
    prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
    postfix = '_c' // ctext // '_k' // samptext // '.fits'

    ! Initialize index-based sky map and mask
    call self%pixcache%init_map_mask(map_in, self%bitmask, map_gain=map_gain)
    call update_status(status, "tod_cache"//ctext)

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

    ! Fit per-chunk low-level non-linearity parameters
    do i = 1, self%nscan
       ! Skip scan if no accepted data
       if (.not. any(self%scans(i)%d%accept)) cycle

       call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd, nonlin_level=0)
       
        ! Subtract A/B detector crosstalk
        ! Not implemented yet

       call timer%start(TOD_NONLIN, self%band)
       if (.false.) then
          ! estimate A/B detector crosstalk coeficients
          call self%xtalk%estimate_crosstalk_matrix(sd)
          call self%xtalk%remove_crosstalk_signal(sd)
       end if

       ! Estimate modulation baselines; and set modulation phase
       if (self%first_call .and. trim(self%init_from_HDF) == 'none') then
          call sample_hfi_baselines(sd, self, i, handle, subtract_s_tot=.false.)
          call set_modulation_phase(sd, self, i)
       else
          call sample_hfi_baselines(sd, self, i, handle)
       end if
       call demodulate_tod(sd, self, i)
       
       if (.not. self%first_call) then
          call int2string(iter, itertext)
          call int2string(self%scanid(i), scantext)
          do j = 1, self%ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             ! fill gaps and deconvolve rolloff
             call fill_gaps(self, sd%tod(:,j), handle, i, j, sd%mask(:,j), sd%s_tot(:,j,0,1), sd%pix(:,:,1),nomono=.true.,filling='white')!,&
                            !& ps_output = 'init_' // itertext // '_' // scantext)
             call deconvolve_rolloff(self, sd%tod(:,j), i, j, sd%s_tot(:,j,0,1), sd%mask(:,j), nomono=.true.)!,&
                                     !& ps_output = itertext // '_' // scantext)
          end do
       end if
       call timer%stop(TOD_NONLIN, self%band)

       ! Fix dc level jumps 
       call self%stitch_hfi_dc_level(i, sd)

       ! Dark bolometer drift correction
       call self%hfi_dark_correction(i, sd)       

       ! 4k Line corrections
       call self%estimate_hfi_4k_lines(i, sd)
 
       ! Clean up
       call dealloc_scan_data(sd)
    end do
    call timer%start(TOD_WAIT, self%band)
    call mpi_barrier(self%comm, ierr) ! Improve timing information
    call timer%stop(TOD_WAIT, self%band)

    
!!$    call update_status(status, "tod_nonlin"//ctext)

    ! Fit global timestream contaminants 

    ! Subtract cosmic ray contribution
!!$    do j=1, self%ndet
!!$
!!$      call init_det_data_singlehorn(dd, self, j)
!!$
!!$      call self%cray(j)%p%build_cray_templates()
!!$
!!$      do i=1, self%nscan
!!$        call populate_sd_from_dd(sd, dd, i, j)
!!$
!!$        call self%cray(j)%p%fit_cray_amplitudes(sd%tod(j,:), sd%s_inst(j, :))
!!$
!!$        call dealloc_scan_data(sd)
!!$      end do
!!$
!!$      call dd%dealloc
!!$    end do

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
       call update_status(status, "tod_calib"//ctext)
    end if
    
    ! Sample ADC and baseline parameters jointly
    if (sample_adc) then
       ! Initialize ADC objects before first call
       if (.not. associated(self%adc(1)%p)) then          
          call self%compute_adu_range
          do i = 1, self%ndet
             self%adc(i)%p => comm_adc_binfit(self%comm, self%label(i), 16, &
                  & self%adu_range(i,1), self%adu_range(i,2), 40)
          end do
       end if

       ! Sample ADC parameters
!!$       do i = 1, self%ndet
!!$          call self%sample_adc_and_baselines(handle, i, map_sky(:,:,i,1), procmask)
!!$       end do
       call update_status(status, "tod_adc"//ctext)
    end if

    ! Sample baselines -- MUST IMMEDIATELY FOLLOW ADC SAMPLER
    do i = 1, self%nscan
       if (.not. any(self%scans(i)%d%accept)) cycle
       call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd, nonlin_level=0)
       call sample_hfi_baselines(sd, self, i, handle)
       call dealloc_scan_data(sd)
    end do
    call update_status(status, "tod_base"//ctext)
    
    ! Create pixel histograms
    if (self%first_call) call compute_tod_pixhist(self)
    call update_status(status, "tod_hist"//ctext)
    
    ! Prepare intermediate data structures
    !call binmap%init(self, .true., .false., nplus2=.false.)
    call binmap%init(self, .true., sample_rel_bandpass, nplus2=.false.)
    if (sample_abs_bandpass .or. sample_rel_bandpass) then
       allocate(chisq_S(self%ndet,size(delta,3)))
       chisq_S = 0.d0
    end if
    if (output_scanlist) then
       allocate(slist(self%nscan))
       slist   = ''
    end if

    ! Fit higher-level corrections
    if (self%myid == 0) write(*,*) '   --> Sampling ncorr, xi_n, maps'
    call update_status(status, "tod_preloop"//ctext)
    do i = 1, self%nscan
       
       ! Skip scan if no accepted data
       if (.not. any(self%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd)

       ! Create dynamic mask
       if (make_dyn_mask) then
          ! Estimate sigma0 for masking
          if (trim(self%init_from_HDF) == 'none') then
             call sample_noise_psd(self, sd, handle, only_sigma0=.true., dec_wn=dec_wn)
          end if

          ! Create mask
          do j = 1, sd%ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             call self%create_dynamic_mask(sd, j, flag_threshold)
          end do
          call dealloc_scan_data(sd)
          if (.not. any(self%scans(i)%d%accept)) cycle

          ! Update scan data with new flagging
          call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd)
       end if

       ! Sample correlated noise
       if (sample_ncorr) then
          call sample_n_corr(self, sd, handle)
          if (sample_xi_n) then
             allocate(freqmask(0:sd%ntod/2))
             freqmask = 1.
             call create_spin_freqmask(real(self%samprate,sp), self%f_spin, self%f_spin/10., 3.0, freqmask) ! Spin harmonics
             call create_spin_freqmask(real(self%samprate,sp), 10., 1., 85., freqmask) ! 4K harmonics
             call sample_noise_psd(self, sd, handle, freqmask=freqmask)
             deallocate(freqmask)
          else
             call sample_noise_psd(self, sd, handle, only_sigma0=.true., dec_wn=dec_wn)
          end if
       else
          call sample_noise_psd(self, sd, handle, only_sigma0=.true., dec_wn=dec_wn)
       end if

       
       ! Compute chisquare
       do j = 1, sd%ndet
          if (.not. self%scans(i)%d(j)%accept) cycle
          call self%compute_tod_chisq(sd, j)
       end do

       ! Select data
       if (select_data) then
          call remove_bad_data(self, i, sd%flag)

          ! Count number of unmasked samples outside the processing mask; for ADC sampling
          do j = 1, sd%ndet
             if (self%scans(i)%d(j)%accept) then
                self%scans(i)%d(j)%nsamp_unmasked = sum(sd%mask(:,j))
             else
                self%scans(i)%d(j)%nsamp_unmasked = 0
             end if
          end do
       end if

       ! Compute chisquare for bandpass fit
       if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)

       ! Compute calibrated TOD for mapmaking
       allocate(d_calib(binmap%nout,sd%ntod, sd%ndet))
       d_calib = 0.d0
       call compute_calibrated_data(self, i, sd, d_calib)
       
       ! Bin TOD
       call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib, binmap)
       
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
!!$       call mpi_finalize(ierr)
!!$       stop
    
    if (select_data) then
       ! Remove data based on a gliding RMS window cut for each of the listed
       ! criteria
       !                           half-window  [chisq, sigma0, fknee, alpha, base, base1, base2]
       call remove_tod_outliers(self, 100,      [5.,    5.,     5.,    5.,    0.,   5.,    5.   ])
       
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
       call self%report_dynamic_mask_stats
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

    ! Clean up
    call binmap%dealloc()
    call update_status(status, "tod_binmap4"//ctext)
    if (allocated(slist)) deallocate(slist)
    if (self%correct_sl) then
       do i = 1, self%ndet
          do h = 1, self%nhorn
             call self%slconv(i,h)%p%dealloc(); deallocate(self%slconv(i,h)%p)
          end do
       end do
    end if

    ! Parameter to check if this is first time routine has been
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

    call read_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'polEff', self%pol_eff(band))
    self%pol_eff(band) = self%pol_eff(band) * 0.01d0 ! Stored as percentage in the instrument file for now

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
          b1 = b1 + self%tod(j,i)
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
       if (tod%scanid(scan) == 1151) write(58,*)
       A2 = 0.d0; b2 = 0.d0
       do j = 2, self%ntod, 2
          if (self%mask(j,i) == 0) cycle
          A2 = A2 + 1.d0
          b2 = b2 + self%tod(j,i)
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
             mu1 = mu1 + (self%tod(j,  i)-tod%scans(scan)%d(i)%baseline1) - &
                       & (self%tod(j+1,i)-tod%scans(scan)%d(i)%baseline2)
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
              self%tod(j,i) =  sgn*(self%tod(j,i) - tod%scans(scan)%d(d)%baseline1)
              !self%tod(j,i) =  (self%tod(j,i) - tod%scans(scan)%d(i)%baseline1)
           else
              self%tod(j,i) = -sgn*(self%tod(j,i) - tod%scans(scan)%d(d)%baseline2)
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
                self%adu_range(i,1) = min(self%adu_range(i,1), nint(sd%tod(j,i)))
                self%adu_range(i,2) = max(self%adu_range(i,2), nint(sd%tod(j,i)))
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

    integer(i4b) :: i, scan

    scan = sd%scan
    
    ! Apply ADC corrections to raw self%tod
    if (associated(self%adc(1)%p)) then
       do i = 1, self%ndet
          if (.not. self%scans(scan)%d(i)%accept) cycle
          call self%adc(i)%p%Q2As(92.)
          call self%adc(i)%p%As2F
          call self%adc(i)%p%adc_correct(self%scanid(scan), i, sd%tod(:,i), &
               & self%scans(scan)%d(i)%N_psd%sigma0, &
               & mask=iand(sd%flag(:,i), self%flag0) .eq.0)
       end do
    end if

    ! Demodulate TOD
    call demodulate_tod(sd, self, scan)
    

    ! Fill gaps and deconvolve high frequency rolloff
!!$    if (skip_nonlin > 2 .and. present(handle)) then
!!$       if (.not. self%first_call) then
!!$          do i = 1, self%ndet
!!$             if (.not. self%scans(scan)%d(i)%accept) cycle
!!$             call fill_gaps(self, sd%tod(:,i), handle, scan, i, sd%mask(:,i), sd%s_tot(:,i), sd%pix(:,:,1), nomono=.true.,filling='white')
!!$             call deconvolve_rolloff(self, sd%tod(:,i), scan, i, sd%s_tot(:,i), sd%mask(:,i), nomono=.true.)
!!$          end do
!!$       end if
!!$    end if


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


  module subroutine estimate_hfi_4k_lines(self, scan, sd)
    !  Construct and apply HFI instrument-specific corrections
    !  from 4k lines
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

    !estimate 4k line signal

  end subroutine estimate_hfi_4k_lines

  module subroutine deconvolve_rolloff(self, tod, scan, i_det, s_sub, mask, nomono, ps_output)
    implicit none
    class(comm_hfi_tod),                       intent(inout) :: self
    real(sp),                   dimension(1:), intent(inout) :: tod
    integer(i4b),                              intent(in)    :: scan, i_det
    real(sp),                   dimension(1:), intent(in)    :: s_sub
    real(sp),         optional, dimension(1:), intent(in)    :: mask
    logical(lgt),     optional,                intent(in)    :: nomono
    character(len=*), optional,                intent(in)    :: ps_output

    integer(i4b) :: i, j, k, l, n, nbin, ntod, nomp, nfft, err
    integer*8    :: plan_fwd, plan_back
    logical(lgt) :: nomono_
    real(sp)     :: gain, N_wn, samprate, dnu, rolloff_scale, eval_spline
    type(spline_type) :: rolloff_filter
    integer(i4b), allocatable, dimension(:)   :: bin_count
    real(sp),     allocatable, dimension(:)   :: d_prime, dt, bin_sum
    complex(spc), allocatable, dimension(:)   :: dv
    real(sp),     allocatable, dimension(:,:) :: ps
    real(dp),     allocatable, dimension(:,:) :: bin_spec

    nomono_ = .false.; if (present(nomono)) nomono_ = nomono

    ntod = self%scans(scan)%ntod
    allocate(d_prime(ntod))
    gain     = self%scans(scan)%d(i_det)%gain  ! Gain in V / K

    ! Prepare TOD residual
    d_prime = tod - gain * s_sub

    ! Remove monopole if requested by user
    if (nomono_ .and. present(mask)) d_prime = d_prime -  sum(d_prime*mask)/sum(mask)

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

    ! Binning
    n = size(ps(:,1))
    dnu = 1.d0 ! Hz
    nbin = (ps(n,1) - ps(1,1))/dnu
    allocate(bin_sum(nbin), bin_count(nbin))
    allocate(bin_spec(nbin,2))

    bin_sum = 0.d0; bin_count = 0
    bin_spec(nbin,2) = 0.d0
    do i = 1, n
       j = (ps(i,1) - ps(1,1))/dnu + 1
       if (j >= 1 .and. j <= nbin) then
          bin_sum(j)   = bin_sum(j) + ps(i,2)
          bin_count(j) = bin_count(j) + 1
       end if
    end do

    do j = 1, nbin
       bin_spec(j,1) = ps(1,1) + (j-0.5d0)*dnu
       if (bin_count(j) > 0) bin_spec(j,2) = bin_sum(j) / bin_count(j)
    end do

    deallocate(bin_sum,bin_count)

    ! Deconvolve high-frequency rolloff
    rolloff_scale = 0.d0
    k = 0
    do j = 1, nbin
       if (bin_spec(j,1) >= 61.d0 .and. bin_spec(j,1) <= 69.d0) then
          k = k + 1
          rolloff_scale = rolloff_scale + bin_spec(j,2)
          !write(*,*) 'scan: ',self%scanid(scan),' ; rolloff_scale in loop = ',rolloff_scale
       end if
    end do
    rolloff_scale = rolloff_scale / k
    !write(*,*) 'scan: ',self%scanid(scan),' ; final rolloff_scale = ',rolloff_scale

    call spline(rolloff_filter, bin_spec(:,1), bin_spec(:,2))
    deallocate(bin_spec)

    do l = 1, n
       if (ps(l,1) > 75.d0) then
          eval_spline = real(splint(rolloff_filter,dble(ps(l,1))),sp) / rolloff_scale
          eval_spline = max(eval_spline,0.0001)
          dv(l) = dv(l) / sqrt(eval_spline)
          ps(l,2) = abs(dv(l)) ** 2 / ntod
       end if
    end do

    call timer%start(TOT_FFT)
    call sfftw_execute_dft_c2r(plan_back, dv, dt)
    call timer%stop(TOT_FFT)

    if(present(ps_output)) then
       open(58,file='deconv_ps_' // ps_output // '.dat', recl=1024)
       do l = 1, n-1
          write(58,*) ps(l,1), ps(l,2)
       end do
       close(58)
    end if

    dt          = dt / nfft
    d_prime = dt(1:ntod)
    tod = d_prime + gain * s_sub

    !if(present(ps_output)) then
    !   open(58,file='deconv_tod_' // ps_output // '.dat', recl=1024)
    !   do l = 1, ntod
    !      write(58,*) l/samprate, d_prime(l)
    !   end do
    !   close(58)
    !end if


    deallocate(dt, dv, d_prime)
    call free_spline(rolloff_filter)
    call dfftw_destroy_plan(plan_fwd)
    call dfftw_destroy_plan(plan_back)

    ! Set new wn level
    N_wn = 1.d30
    allocate(bin_sum(nbin), bin_count(nbin))
    allocate(bin_spec(nbin,2))
    bin_sum = 0.d0; bin_count = 0
    bin_spec(nbin,2) = 0.d0
    do i = 1, n
       j = (ps(i,1) - ps(1,1))/dnu + 1
       if (j >= 1 .and. j <= nbin) then
          bin_sum(j)   = bin_sum(j) + ps(i,2)
          bin_count(j) = bin_count(j) + 1
       end if
    end do

    do j = 1, nbin
       bin_spec(j,1) = ps(1,1) + (j-0.5d0)*dnu
       if (bin_count(j) > 0) bin_spec(j,2) = bin_sum(j) / bin_count(j)
    end do

    deallocate(bin_sum,bin_count, ps)

    do i = 1, nbin
       if (bin_spec(i,2) < N_wn) N_wn = bin_spec(i,2)
    end do

    self%scans(scan)%d(i_det)%N_psd%sigma0 = sqrt(abs(N_wn))
    deallocate(bin_spec)

  end subroutine deconvolve_rolloff

  module subroutine fill_gaps(self, tod, handle, scan, i_det, mask, s_sub, pix, nomono, dospike, ps_output, filling)
    implicit none
    class(comm_hfi_tod),                     intent(in)    :: self
    real(sp),              dimension(1:),    intent(inout) :: tod
    type(planck_rng),                        intent(inout) :: handle
    integer(i4b),                            intent(in)    :: scan, i_det
    real(sp),              dimension(1:),    intent(in)    :: mask, s_sub
    integer(i4b),          dimension(1:,1:), intent(in)    :: pix
    logical(lgt),                  optional, intent(in)    :: nomono
    logical(lgt),                  optional, intent(in)    :: dospike
    character(len=*),              optional, intent(in)    :: ps_output
    character(len=*),              optional, intent(in)    :: filling

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
    filling_ = 'none'; if (present(filling)) filling_ = filling

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
       call dfftw_destroy_plan(plan_fwd)
       call dfftw_destroy_plan(plan_back)
    end if

    tod = d_prime + gain * s_sub
    deallocate(d_prime)

  end subroutine fill_gaps


  module subroutine sample_adc_and_baselines(self, handle, det, map_sky, procmask)
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
    real(sp),          dimension(1:,1:),       intent(in)    :: map_sky
    real(sp),          dimension(0:),          intent(in)    :: procmask

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

    oper = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
         & SD_SKY,SD_BP,SD_ORB,SD_INST,SD_DARK])
    
    ! Count number of unmasked and decimated samples
    allocate(numsamp(self%nscan), sigma0(self%nscan))
    do i = 1, self%nscan
       numsamp(i) = self%scans(i)%d(det)%nsamp_unmasked / decimation
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
       call init_scan_data(self, i, oper, TODMASK_PROC, sd, nonlin_level=0)

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

      powell_chisq_adc_hfi = chisq_adc_hfi(real(x,sp))

    end function powell_chisq_adc_hfi


    function chisq_adc_hfi(x, ndof) result (chisq)
      implicit none
      real(sp), dimension(:), intent(in),  optional :: x
      integer(i8b),           intent(out), optional :: ndof
      real(dp)                                      :: chisq

      integer(i4b) :: i, j, n
      integer(i8b) :: ndof_sub, ndof_tot, k1, k2
      real(dp)     :: chisq_sub, A, b
      real(sp)     :: baseline
      logical(lgt) :: output_ndof
      real(sp), allocatable, dimension(:) :: p

      allocate(p(self%adc(det)%p%npar_adc))
      if (self%myid == 0) then
         flag = 1; if (present(ndof)) flag = 2
         call mpi_bcast(flag, 1, MPI_INTEGER, 0, self%comm, ierr)
         p = x
      end if
      call mpi_bcast(p, size(p), MPI_REAL, 0, self%comm, ierr)
      output_ndof = (flag == 2)

      ! Check priors
      if (any(p < 0.0)) then
         chisq = 1.d30
         deallocate(p)
         return
      end if
      
      ! Update ADC parameters; assume constant sigma0 fir niw
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
!!$               write(58,*) k1+j-1, tod(k1+j-1), tod_corr(j), tod(k1+j-1)-tod_corr(j), (tod(k1+j-1)-tod_corr(j))/sigma0(i), sigma0(i)
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

      if (self%myid == 0) write(*,*) 'adc chisq =', chisq, ', p = ',  real(p,sp)

    end function chisq_adc_hfi

  end subroutine sample_adc_and_baselines
  

end submodule comm_tod_hfi_mod
