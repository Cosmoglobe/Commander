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

    logical(lgt), dimension(:,:), allocatable :: correlations

    ! Allocate object
    allocate(c)

    ! Set up noise PSD type and priors
    c%freq            = cpar%ds_label(id_abs)
    c%n_xi            = 3
    !c%noise_psd_model = 'white'       ! Using white noise until we get better estimates of the actual noise PSD
    c%noise_psd_model = 'oof'       ! Not fitted parameters yet
    allocate(c%xi_n_P_uni(c%n_xi,2))
    allocate(c%xi_n_P_rms(c%n_xi))
    allocate(c%xi_n_nu_fit(c%n_xi,2))

    c%xi_n_P_uni(2,:)  = [0.001d0, 5.0d0]  ! fknee
    c%xi_n_P_uni(3,:)  = [-2.5d0, -0.4d0]   ! alpha
    c%xi_n_nu_fit(1,:) = [3.d0, 10.d0] 
    c%xi_n_nu_fit(2,:) = [0.d0, 3d0]
    c%xi_n_nu_fit(3,:) = [0.d0, 3d0]
    c%xi_n_P_rms       = [-1.d0, 0.1d0, 0.1d0] ! [sigma0, fknee, alpha]; sigma0 is not used
    
    

    !c%xi_n_P_rms      = [-1.d0] ! [sigma0]; sigma0 is not used

    c%n_cray_temps    = 3
    c%ndiode = 1

    ! Initialize common parameters
    call c%tod_constructor(cpar, id, id_abs, info, tod_type)

    ! Initialize instrument-specific parameters
    c%samprate_lowres = 18.  ! Lowres samprate in Hz;  10 times lower than the intrinsic HFI rate for now
    c%nhorn           = 1
    c%compressed_tod  = .true.
    c%correct_sl      = .false.
    c%correct_orb     = .true.
    c%apply_inst_corr = .true.
    c%orb_4pi_beam    = .false.
    c%symm_flags      = .false.
    c%sample_zodi     = cpar%sample_zodi .and. c%subtract_zodi ! Sample zodi parameters
    c%nmaps           = info%nmaps
    c%ndet            = num_tokens(cpar%ds_tod_dets(id_abs), "," )
    c%ntime           = 1
    c%partner         = -1
    !TODO: set the number of dark bolometers to be correct
    c%ndark           = 1

    nmaps_beam        = 3
    pol_beam          = .true.
    c%nside_beam      = 512
    c%nside_pixhist   = 64

    ! Channel specific parameters
    if (c%freq(1:3) == "100") then
       c%chisq_threshold  = 100.d0
       c%sigma0_threshold = 100.d0
       c%accept_threshold = 0.5d0
       c%correct_sl       = .false.
    else
       c%chisq_threshold  = 100.d0 
       c%accept_threshold = 0.5d0
       c%correct_sl       = .false.
    end if

    
    ! Get detector labels
    call get_tokens(cpar%ds_tod_dets(id_abs), ",", c%label)
    
    ! Read the actual TOD
    call c%read_tod(c%label)
    
    ! Initialize bandpass mean and proposal matrix
    !call c%initialize_bp_covar(trim(cpar%datadir)//'/'//cpar%ds_tod_bp_init(id_abs))

    ! Construct lookup tables
    call c%precompute_lookups()

    ! Load the instrument file
    call c%load_instrument_file(c%nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

    ! Allocate sidelobe convolution data structures
    if(c%correct_sl) allocate(c%slconv(c%ndet), c%orb_dp)
    
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
    type(map_ptr),       dimension(1:,1:),    intent(inout), optional :: map_gain       ! (ndet,1)

    real(dp)            :: t1, t2
    integer(i4b)        :: i, j, k, l, ierr, ndelta, nside, npix, nmaps
    logical(lgt)        :: select_data, output_scanlist, output_zodi_comps
    logical(lgt)        :: sample_gain, sample_ncorr, sample_abs_bandpass, sample_rel_bandpass, sample_zodi, sample_adc, make_dyn_mask
    type(comm_binmap)   :: binmap
    type(comm_scandata) :: sd
    type(comm_detdata)  :: dd
    character(len=4)    :: ctext, myid_text
    character(len=6)    :: samptext, scantext
    character(len=512)  :: prefix, postfix, prefix4D, filename
    character(len=512), allocatable, dimension(:) :: slist
    real(sp),              dimension(9)       :: flag_threshold
    real(sp), allocatable, dimension(:)       :: procmask, procmask2, procmask_zodi, sigma0
    real(sp), allocatable, dimension(:,:)     :: s_buf
    real(sp), allocatable, dimension(:,:,:)   :: d_calib
    real(sp), allocatable, dimension(:,:,:,:) :: map_sky, m_gain
    real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf

    call int2string(iter, ctext)
    call update_status(status, "tod_start"//ctext)
    
    ! Toggle optional operations
    sample_rel_bandpass   = size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
    sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses
    sample_gain           = iter  > 0 !.true.                 
    make_dyn_mask         = iter == 2
    sample_ncorr          = iter  > 2 !.true.
    select_data           = iter == 3 ! self%first_call  
    sample_adc            = .false. !iter  > 3 !.true.
    sample_zodi           = self%sample_zodi .and. self%subtract_zodi ! Sample zodi parameters
    output_zodi_comps     = self%output_zodi_comps .and. self%subtract_zodi ! Output zodi components
    output_scanlist       = mod(iter-1,10) == 0    ! only output scanlist every 10th iteration

    !                       Pixhist    Single abs/RMS       RMS ranges     Single     Ranges   Pointing
    if (self%freq(1:3) == "100") then
       flag_threshold     = [  1.0,        20.0, 5.0,         1.5, 2.0, -1.0,   -1.0,      -1.0,     -1.0]
    else
       flag_threshold     = [  1.0,        20.0, 5.0,         1.5, 2.0, -1.0,   1.0,      -1.0,     -1.0]
    end if

    
    ! Initialize local variables
    ndelta          = size(delta,3)
    self%n_bp_prop  = ndelta
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

    ! Distribute maps
    allocate(map_sky(nmaps,self%nobs,0:self%ndet,ndelta))
    allocate(m_gain(nmaps,self%nobs,0:self%ndet,1))
    call distribute_sky_maps(self, map_gain, 1.e0, m_gain) ! in K_cmb
    call distribute_sky_maps(self, map_in, 1.e0, map_sky) ! already in K_cmb

    ! Distribute processing masks
    allocate(m_buf(0:npix-1,nmaps), procmask(0:npix-1), procmask2(0:npix-1))
    call self%procmask%bcast_fullsky_map(m_buf);  procmask  = m_buf(:,1)
    call self%procmask2%bcast_fullsky_map(m_buf); procmask2 = m_buf(:,1)
    if (self%sample_zodi .and. self%subtract_zodi) then
       allocate(procmask_zodi(0:npix-1))
       call self%procmask_zodi%bcast_fullsky_map(m_buf); procmask_zodi = m_buf(:,1)
    end if
    deallocate(m_buf)

    call map_in(1,1)%p%writeFITS(trim(self%outdir) // "/input_sky_model_"//trim(self%label(1))//".fits")
    !call self%procmask%writeFITS("mask.fits")
    
    ! Precompute far sidelobe Conviqt structures
    if (self%correct_sl) then
       if (self%myid == 0) write(*,*) 'Precomputing sidelobe convolved sky'
       do i = 1, self%ndet
          !TODO: figure out why this is rotated
          call map_in(i,1)%p%YtW()  ! Compute sky a_lms
          self%slconv(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
               & 100, 3, 100, self%slbeam(i)%p, map_in(i,1)%p, 2)
       end do
    end if
!    write(*,*) 'qqq', self%myid
!    if (.true. .or. self%myid == 78) write(*,*) 'a', self%myid, self%correct_sl, self%ndet, self%slconv(1)%p%psires
!!$    call mpi_finalize(ierr)
!!$    stop

    call update_status(status, "tod_init")

    !------------------------------------
    ! Perform main sampling steps
    !------------------------------------

    ! estimate A/B detector crosstalk coeficients
    call self%xtalk%estimate_crosstalk_matrix()

    ! Fit per-chunk low-level non-linearity parameters
    do i = 1, self%nscan
        ! Skip scan if no accepted data
        if (.not. any(self%scans(i)%d%accept)) cycle
        call init_scan_data_singlehorn(sd, self, i, map_sky, m_gain, procmask, procmask2, procmask_zodi, skip_nonlin=.true., darkdata=.true.)

!!$        open(58,file='res.dat', recl=1024)
!!$        do j = 1, sd%ntod
!!$           write(58,*) j, sd%tod(j,1), sd%pix(j,1,1), sd%flag(j,1)
!!$        end do
!!$        close(58)
        
        ! Subtract A/B detector crosstalk
        ! Not implemented yet

       call self%xtalk%remove_crosstalk_signal(sd, i)

       ! Estimate modulation baselines; and set modulation phase
       if (self%first_call) then
          call sample_hfi_baselines(sd, self, i, handle, subtract_s_tot=.false.)
          call set_modulation_phase(sd, self, i)
       end if

       ! Fix dc level jumps 
       call self%stitch_hfi_dc_level(i, sd)

       ! Dark bolometer drift correction
       call self%hfi_dark_correction(i, sd)       

       ! 4k Line corrections
       call self%estimate_hfi_4k_lines(i, sd)

 
       ! Clean up
       call dealloc_scan_data(sd)
    end do
    
    ! Fit global timestream contaminants 

    ! Subtract cosmic ray contribution
    do j=1, self%ndet

      call init_det_data_singlehorn(dd, self, j)

      call self%cray(j)%p%build_cray_templates()

!!$      do i=1, self%nscan
!!$        call populate_sd_from_dd(sd, dd, i, j)
!!$
!!$        call self%cray(j)%p%fit_cray_amplitudes(sd%tod(j,:), sd%s_inst(j, :))
!!$
!!$        call dealloc_scan_data(sd)
!!$      end do

      call dd%dealloc
    end do

    ! Estimate ADC corrections
    !    Not implemented yet

    ! Fit bolometer transfer function parameters?
    !    Not implemented yet     


    ! NOW Sample high level tod components that require cleaned data

    ! Sample gain components in separate TOD loops; marginal with respect to n_corr
    if (sample_gain) then
       ! TODO: Also sample non-linear gain response here?
       call sample_calibration(self, 'abscal', handle, map_sky, m_gain, procmask2, procmask2)
       call sample_calibration(self, 'relcal', handle, map_sky, m_gain, procmask, procmask2)
       !call sample_calibration(self, 'deltaG', handle, map_sky, m_gain, procmask, procmask2)
    end if

    ! Sample ADC and baseline parameters jointly
    if (sample_adc) then
       do i = 1, self%ndet
          call self%sample_adc_and_baselines(handle, i, map_sky(:,:,i,1), procmask)
       end do
    end if

    ! Sample baselines -- MUST IMMEDIATELY FOLLOW ADC SAMPLER
    do i = 1, self%nscan
       if (.not. any(self%scans(i)%d%accept)) cycle
       call init_scan_data_singlehorn(sd, self, i, map_sky, m_gain, procmask, procmask2, procmask_zodi, skip_nonlin=.true., darkdata=.true.)
       call sample_hfi_baselines(sd, self, i, handle)
       call dealloc_scan_data(sd)
    end do
    
    ! Create pixel histograms
    if (self%first_call) then
       call compute_tod_pixhist(self, map_sky, m_gain, procmask, procmask2)
    end if


    ! Prepare intermediate data structures
    call binmap%init(self, .true., .false.)!, nplus2=.true.)
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
    do i = 1, self%nscan
       
       ! Skip scan if no accepted data
       if (.not. any(self%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       if (sample_rel_bandpass) then
          call init_scan_data_singlehorn(sd, self, i, map_sky, m_gain, procmask, procmask2, procmask_zodi, init_s_bp=.true., init_s_bp_prop=.true.)
       else if (sample_abs_bandpass) then
          call init_scan_data_singlehorn(sd, self, i, map_sky, m_gain, procmask, procmask2, procmask_zodi, init_s_bp=.true., init_s_sky_prop=.true.)
       else
          call init_scan_data_singlehorn(sd, self, i, map_sky, m_gain, procmask, procmask2, procmask_zodi, init_s_bp=.true.)
       end if

       if (self%myid == 0) then
          open(58,file='tod_gain.dat', recl=1024)
          do j = 1, sd%ntod
             write(58,*) j, sd%s_tot(j,:), sd%tod(j,:)
          end do
          close(58)
       end if
       
       ! Create dynamic mask
       if (make_dyn_mask) then
          ! Estimate sigma0 for masking
          sd%n_corr = 0.d0
          call sample_noise_psd(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr, only_sigma0=.true.)

          ! Create mask
          do j = 1, sd%ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             call self%create_dynamic_mask(i, j, sd%pix(:,j,1), sd%tod(:,j), (sd%tod(:,j)-real(self%scans(i)%d(j)%gain,sp)*sd%s_tot(:,j))/self%scans(i)%d(j)%N_psd%sigma0, &
                  & sd%mask(:,j), sd%flag(:,j), flag_threshold)
          end do
          call dealloc_scan_data(sd)
          if (.not. any(self%scans(i)%d%accept)) cycle

          ! Update scan data with new flagging
          if (sample_rel_bandpass) then
             call init_scan_data_singlehorn(sd, self, i, map_sky, m_gain, procmask, procmask2, procmask_zodi, init_s_bp=.true., init_s_bp_prop=.true.)
          else if (sample_abs_bandpass) then
             call init_scan_data_singlehorn(sd, self, i, map_sky, m_gain, procmask, procmask2, procmask_zodi, init_s_bp=.true., init_s_sky_prop=.true.)
          else
             call init_scan_data_singlehorn(sd, self, i, map_sky, m_gain, procmask, procmask2, procmask_zodi, init_s_bp=.true.)
          end if

          ! Count number of unmasked samples outside the processing mask; for ADC sampling
          do j = 1, sd%ndet
             self%scans(i)%d(j)%nsamp_unmasked = sum(sd%mask(:,j))
          end do
       end if

       ! Sample correlated noise
       if (sample_ncorr) then
          !call sample_n_corr(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr, sd%pix(:,:,1), nomono=.true.)
          call sample_n_corr(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr, sd%pix(:,:,1))
          call sample_noise_psd(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr)
       else
          sd%n_corr = 0.d0
          call sample_noise_psd(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr, only_sigma0=.true.)
       end if
       
       ! Compute chisquare
       do j = 1, sd%ndet
          if (.not. self%scans(i)%d(j)%accept) cycle
          call self%compute_tod_chisq(i, j, sd%mask(:,j), sd%s_sky(:,j), sd%s_sl(:,j) + sd%s_orb(:,j), sd%n_corr(:,j), tod=sd%tod(:, j))
       end do

       ! Select data
       if (select_data) then
          call remove_bad_data(self, i, sd%flag)
          do j = 1, sd%ndet
             if (.not. self%scans(i)%d(j)%accept) self%scans(i)%d(j)%nsamp_unmasked = 0
          end do
       end if

       ! Compute chisquare for bandpass fit
       !if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)

       ! Compute calibrated TOD for mapmaking
       allocate(d_calib(self%output_n_maps,sd%ntod, sd%ndet))
       d_calib = 0.d0
       call compute_calibrated_data(self, i, sd, d_calib)

!!$       if (self%scanid(i) == 500) then
!!$          open(58,file='res'//samptext//'.dat', recl=1024)
!!$          do j = 1, sd%ntod
!!$             write(58,*) j, sd%tod(j,1), sd%n_corr(j,1), d_calib(1,j,1), d_calib(2,j,1), 1-(sd%flag(j,1)/maxval(sd%flag(:,1))), self%psi(sd%psi(j,1,1))*RAD2DEG, self%psi(sd%psi(j,2,1))*RAD2DEG, self%psi(sd%psi(j,3,1))*RAD2DEG, self%psi(sd%psi(j,4,1))*RAD2DEG
!!$          end do
!!$          close(58)
!!$       end if

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

        ! Initialize ADC objects
    if (select_data) then
       call self%compute_adu_range
       do i = 1, self%ndet
          self%adc(i)%p => comm_adc_binfit(self%comm, self%label(i), 16, &
               & self%adu_range(i,1), self%adu_range(i,2), 40)
       end do
    end if

    if (self%myid == 0) write(*,*) '   --> Finalizing maps, bp'
    if (make_dyn_mask) then
       call self%report_dynamic_mask_stats
    end if
    
    ! Output latest scan list with new timing information
    if (output_scanlist) call self%output_scan_list(slist)

    ! Solve for maps
    ! TODO: update mapmaker to n+2 maps
    ! TODO: handle bolometer transfer function in the mapmaking
    call synchronize_binmap(binmap, self)
    if (sample_rel_bandpass) then
       if (self%nmaps > 1) then
         !call finalize_binned_map_nplus2(self, binmap, rms_out, 1.d0, chisq_S=chisq_S, mask=procmask2, correct_transfer=.true.)
         call finalize_binned_map(self, binmap, rms_out, 1.d0, chisq_S=chisq_S, mask=procmask2)
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

    ! Sample bandpass parameters
    if (sample_rel_bandpass .or. sample_abs_bandpass) then
       call sample_bp(self, iter, delta, map_sky, handle, chisq_S)
       self%bp_delta = delta(:,:,1)
    end if
   
    ! Output maps to disk
    call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
    call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
    if (self%output_n_maps > 1) call binmap%outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
    if (self%output_n_maps > 2) call binmap%outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
    if (self%output_n_maps > 3) call binmap%outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix))
    if (self%output_n_maps > 4) call binmap%outmaps(5)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
!!$    if (self%output_n_maps > 5) call binmap%outmaps(6)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
!!$    if (self%output_n_maps > 6) call binmap%outmaps(7)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))
!!$    if (self%output_n_maps > 8 .and. self%subtract_zodi .and. output_zodi_comps) then
!!$       do i = 1, zodi_model%n_comps
!!$          call binmap%outmaps(8+i)%p%writeFITS(trim(prefix)//'zodi_'//trim(zodi_model%comp_labels(i))//trim(postfix))
!!$       end do
!!$    endif

    ! Clean up
    call binmap%dealloc()
    if (allocated(slist)) deallocate(slist)
    deallocate(map_sky, procmask, procmask2, m_gain)
    if (self%correct_sl) then
       do i = 1, self%ndet
          call self%slconv(i)%p%dealloc(); deallocate(self%slconv(i)%p)
       end do
    end if

    ! Parameter to check if this is first time routine has been
    self%first_call = .false.

    call update_status(status, "tod_end"//ctext)

  end subroutine process_HFI_tod


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

    sub_s = .true.; if (present(subtract_s_tot)) sub_s = subtract_s_tot
    
    
    ! tod%scans(scan)%d(i)%gain - the gain constant over a scan [real number]
    ! sd = self --- self%s_tot - sky signal model
    ! self%s_tot(self%ntod, self%ndet) - how s_tot structured
    
    if (tod%scanid(scan) == 500) open(58,file='baseline.dat', recl=1024)
    do i = 1, tod%ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle
       sgn = tod%mod_phase(i,scan)
       
       ! Odd samples
       A1 = 0.d0; b1 = 0
       do j = 1, self%ntod, 2
          if (self%mask(j,i) == 0) cycle
          A1 = A1 + 1.d0
          b1 = b1 + self%tod(j,i)
          if (sub_s) b1 = b1 - sgn*tod%scans(scan)%d(i)%gain * self%s_tot(j,i)
          if (tod%scanid(scan) == 500) write(58,*) j, self%tod(j,i), self%tod(j,i) - sgn*tod%scans(scan)%d(i)%gain * self%s_tot(j,i)
       end do
       A1 = A1 / tod%scans(scan)%d(i)%N_psd%sigma0**2
       b1 = b1 / tod%scans(scan)%d(i)%N_psd%sigma0**2
       if (A1 == 0.d0) then
          tod%scans(scan)%d(i)%accept = .false.
          cycle
       end if
       tod%scans(scan)%d(i)%baseline1  = b1/A1 + rand_gauss(handle)/sqrt(A1)

       ! Even samples
       if (tod%scanid(scan) == 500) write(58,*)
       A2 = 0.d0; b2 = 0.d0
       do j = 2, self%ntod, 2
          if (self%mask(j,i) == 0) cycle
          A2 = A2 + 1.d0
          b2 = b2 + self%tod(j,i)
          if (sub_s) b2 = b2 + sgn*tod%scans(scan)%d(i)%gain * self%s_tot(j,i)
          if (tod%scanid(scan) == 500) write(58,*) j, self%tod(j,i), self%tod(j,i) + sgn*tod%scans(scan)%d(i)%gain * self%s_tot(j,i)
       end do
       A2 = A2 / tod%scans(scan)%d(i)%N_psd%sigma0**2
       b2 = b2 / tod%scans(scan)%d(i)%N_psd%sigma0**2
       if (A2 == 0.d0) then
          tod%scans(scan)%d(i)%accept = .false.
          cycle
       end if       
       tod%scans(scan)%d(i)%baseline2  = b2/A2 + rand_gauss(handle)/sqrt(A2)

       !write(*,'(a,i8,3f16.3)') "baseline =", tod%scanid(scan), tod%scans(scan)%d(i)%baseline1, tod%scans(scan)%d(i)%baseline2, sgn

    end do
    if (tod%scanid(scan) == 500) close(58)

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
    
    real(dp) :: mu, n
    integer(i4b) :: i, j
    
    do i = 1, tod%ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle       
       
       mu = 0.d0; n = 0.d0
       do j = 1, self%ntod, 2
          if (iand(self%flag(j,i), tod%flag0) .ne. 0) cycle
          if (self%pix(j,i,1) > 0.48*tod%info%npix .and. self%pix(j,i,1) < 0.52*tod%info%npix) then
             mu = mu + (self%tod(j,i)-tod%scans(scan)%d(i)%baseline1)
             n  = n  + 1.d0
          end if
       end do

       if (n == 0) then
          write(*,*) "set_modulation_phase: no samples crossing the galactic plane. Scan disabled = ", tod%scanid(scan)
          tod%scans(scan)%d(i)%accept = .false.
       else
          mu = mu/n

          ! saving the phase to the tod object
          if (mu < 0.d0) then
             tod%mod_phase(i,scan) = -1
          end if
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

    integer(i4b) :: i, j
    real(sp)     :: sgn

!!$    open(58,file='tod_adc.dat')
!!$    do j = 1, self%ntod
!!$       write(58,*) j, self%tod(j,1)
!!$    end do
!!$    close(58)
    
    do i = 1, tod%ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle       
       sgn = tod%mod_phase(i,scan)
       
       ! Subtract baselines and flip sign of every other sample
       do j = 1, self%ntod
           if (mod(j,2) == 1) then
              self%tod(j,i) =  sgn*(self%tod(j,i) - tod%scans(scan)%d(i)%baseline1)
              !self%tod(j,i) =  (self%tod(j,i) - tod%scans(scan)%d(i)%baseline1)
           else
              self%tod(j,i) = -sgn*(self%tod(j,i) - tod%scans(scan)%d(i)%baseline2)
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
    
    integer(i4b) :: i, j, ntod, scan, ierr
    real(sp),     allocatable, dimension(:)   :: tod
    integer(i4b), allocatable, dimension(:,:) :: pix, psi
    integer(i4b), allocatable, dimension(:)   :: flag

    do i = 1, self%ndet
       self%adu_range(i,:) = [40*2**16,0]
       do scan = 1, self%nscan
          if (.not. self%scans(scan)%d(i)%accept) cycle
          ntod = self%scans(scan)%ntod
          allocate(tod(ntod), pix(ntod,1), psi(ntod,1), flag(ntod))
          call self%decompress_tod(scan, i, tod)
          call self%decompress_pointing_and_flags(scan, i, pix, psi, flag)
          do j = 1, self%scans(scan)%ntod
             if (iand(flag(j), self%flag0) .eq. 0) then
                self%adu_range(i,1) = min(self%adu_range(i,1), nint(tod(j)))
                self%adu_range(i,2) = max(self%adu_range(i,2), nint(tod(j)))
             end if
          end do
          deallocate(tod, pix, psi, flag)
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
  end subroutine initHDF_HFI
  
  module subroutine dumpToHDF_HFI(self, chainfile, path)
    ! 
    ! Writes HFI-specific TOD parameters to existing chain file
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
    class(comm_HFI_tod),                 intent(in)     :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
  end subroutine dumpToHDF_HFI

  module subroutine construct_corrtemp_hfi(self, scan, pix, psi, s, det)
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
    class(comm_hfi_tod),                   intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    integer(i4b),        dimension(:,:),   intent(in)    :: pix, psi
    real(sp),            dimension(:,:),   intent(out)   :: s
    integer(i4b),                          intent(in), optional :: det

    integer(i4b) :: i, j, k, nbin, b
    real(dp)     :: dt, t_tot, t

    s = 0.

  end subroutine construct_corrtemp_hfi

  module subroutine apply_nonlin_corr_hfi(self, scan, sd, det)
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
    class(comm_hfi_tod),                   intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    class(comm_scandata),                  intent(inout) :: sd
    integer(i4b),                          intent(in), optional :: det

    integer(i4b) :: i
    
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

    integer(i4b)        :: i, j, k, flag, ierr, ntod, decimation, offset, ind
    integer(i8b)        :: ntot
    real(dp)            :: chisq
    type(comm_scandata) :: sd
    character(len=4)    :: id
    integer(i8b), allocatable, dimension(:)   :: numsamp
    real(sp),     allocatable, dimension(:)   :: tod, s_tot, phase, tod_corr
    real(sp),     allocatable, dimension(:)   :: sigma0, gain

    decimation = 45 ! Downsampling rate; MUST BE ODD, to explore both parities
    
    ! Count number of unmasked and decimated samples
    allocate(numsamp(self%nscan), sigma0(self%nscan), gain(self%nscan))
    do i = 1, self%nscan
       numsamp(i) = self%scans(i)%d(det)%nsamp_unmasked / decimation
       sigma0(i)  = self%scans(i)%d(det)%N_psd%sigma0
       gain(i)    = self%scans(i)%d(det)%gain
    end do
    ntot = sum(numsamp)

    ! Prepare reduced dataset
    allocate(tod(ntot), s_tot(ntot), phase(ntot), tod_corr(maxval(numsamp)))
    ind = 1
    do i = 1, self%nscan
       if (.not. self%scans(i)%d(det)%accept) then
          numsamp(i) = ind-1
          cycle
       end if
       call init_scan_data_singlehorn_singledet(sd, self, det, i, map_sky, &
            & procmask, skip_nonlin=.true.)
       !call init_scan_data_singlehorn(sd, self, det, i, map_sky, &
       !     & procmask, skip_nonlin=.true.)
       !call init_scan_data_singlehorn(sd, self, i, m_sky, m_sky, procmask, procmask, skip_nonlin=.true.)
       
       ! Decimate data
       j = 0; k = decimation*j+1
       do while (k <= sd%ntod .and. j < numsamp(i))
          if (sd%mask(k,1) == 1.) then
             tod(ind)   = sd%tod(k,1)
             s_tot(ind) = gain(i) * sd%s_tot(k,1)
             if (mod(k,2) == 1) then
                phase(ind) =  self%mod_phase(det,i)
             else
!                if (.not. allocated(self%mod_phase)) then
!                   write(*,*) 'q4', self%myid, k, ind, allocated(self%mod_phase)
!                end if
                phase(ind) = -self%mod_phase(det,i)
             end if
             !if (self%myid == 0) write(*,*) 'tod = ', ind, numsamp(i), tod(ind), s_tot(ind), phase(ind)
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

    call int2string(self%myid, id)
    open(58,file='adc_data_'//trim(self%adc(det)%p%label)//'_id'//id//'.dat', recl=1024)
    do i = 1, ntot
       write(58,*) i, tod(i), s_tot(i), phase(i)
    end do
    close(58)
    
    
    !write(*,*) self%myid, self%nscan, numsamp, ntot
    
    !  Perform sampling
    if (self%myid == 0) then
       call self%adc(det)%p%mcmc_sample_adc(handle, chisq_adc_hfi2)
       ! Release workers
       flag = 0
       call mpi_bcast(flag, 1, MPI_INTEGER, 0, self%comm, ierr)
    else
       do while (.true.)
          call mpi_bcast(flag, 1, MPI_INTEGER, 0, self%comm, ierr)
          if (flag > 0) then
             chisq = chisq_adc_hfi2()
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
    deallocate(tod, s_tot, phase, sigma0, gain, tod_corr)
    
  contains

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

      ! Update ADC parameters; assume constant sigma0 fir niw
      call self%adc(det)%p%param2Q(p)
      call self%adc(det)%p%Q2As(92.)
      call self%adc(det)%p%As2F

!!$      if (self%myid == 134) then
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
         
         ! Apply ADC correction
         tod_corr(1:n) = tod(k1:k2)
         !write(*,*) self%myid, i, n, k1, k2, shape(tod_corr), shape(tod), minval(tod_corr(1:n)), maxval(tod_corr(1:n))
         call self%adc(det)%p%adc_correct(self%scanid(i), det, &
              & tod_corr(1:n), sigma0(i))          

         ! Sample baseline for and demodulate positive parity samples
         A = 0.d0; b = 0.d0
         do j = 1, n
            if (phase(k1+j-1) == -1.) cycle
            A = A + 1.d0
            b = b + tod_corr(j) - s_tot(k1+j-1)
         end do
         baseline = b/A + sigma0(i)*rand_gauss(handle)/sqrt(A)
         do j = 1, n
            if (phase(k1+j-1) == -1.) cycle
            tod_corr(j) = tod_corr(j) - baseline
         end do

         ! Sample baseline for and demodulate negative parity samples
         A = 0.d0; b = 0.d0
         do j = 1, n
            if (phase(k1+j-1) == 1.) cycle
            A = A + 1.d0
            b = b + tod_corr(j) + s_tot(k1+j-1)
         end do
         baseline = b/A + sigma0(i)*rand_gauss(handle)/sqrt(A)
         do j = 1, n
            if (phase(k1+j-1) == 1.) cycle
            tod_corr(j) = -tod_corr(j) + baseline
         end do

!!$         if (self%myid == 134) then
!!$            do j = 1, n
!!$               write(58,*) k1+j-1, tod(k1+j-1), tod_corr(j), s_tot(k1+j-1), tod_corr(j)-s_tot(k1+j-1)
!!$            end do
!!$         end if
         
         ! Compute chisq
         do j = 1, n
            !write(*,*) i, j, tod_corr(j), s_tot(k1+j-1), sigma0(i)
            chisq_sub = chisq_sub + (tod_corr(j)-s_tot(k1+j-1))**2/sigma0(i)**2
         end do
      end do

!!$      if (self%myid == 134) close(58)

!!$      if (ntot > 0) write(*,*) 'adc chisq =', chisq_sub/ntot, ', myid = ', self%myid
      
      call mpi_reduce(chisq_sub, chisq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
      if (output_ndof) then
         call mpi_reduce(ntot, ndof_tot, 1, MPI_INTEGER8, MPI_SUM, 0, self%comm, ierr)
         if (self%myid == 0) ndof = ndof_tot
      end if

    end function chisq_adc_hfi


    function chisq_adc_hfi2(x, ndof) result (chisq)
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

      ! Update ADC parameters; assume constant sigma0 fir niw
      call self%adc(det)%p%param2Q(p)
      call self%adc(det)%p%Q2As(92.)
      call self%adc(det)%p%As2F

      if (self%myid == 0 .and. det==2) then
         call int2string(self%myid, id)
         open(58,file='adc_reddata_'//trim(self%adc(det)%p%label)//'_id'//id//'.dat', recl=1024)
      end if
      
      ! Evaluate chisq
      chisq_sub = 0.d0
      k2 = 0 
      do i = 1, self%nscan
         ! Skip scan if no accepted data
         if (.not. self%scans(i)%d(det)%accept) cycle
         n = numsamp(i); k1 = k2+1; k2 = k2+n
         
         ! Apply ADC correction
         !tod_corr(1:n) = tod(k1:k2)
         !write(*,*) self%myid, i, n, k1, k2, shape(tod_corr), shape(tod), minval(tod_corr(1:n)), maxval(tod_corr(1:n))
         !call self%adc(det)%p%adc_correct(self%scanid(i), det, &
         !     & tod_corr(1:n), sigma0(i))          

         ! Add modulated signal and baseline
         do j = 1, n
            if (gain(i) > 0.) then
               if (phase(k1+j-1) == 1.) then
                  tod_corr(j) =  s_tot(k1+j-1) + self%scans(i)%d(det)%baseline1
               else
                  tod_corr(j) = -s_tot(k1+j-1) + self%scans(i)%d(det)%baseline2
               end if
            else
               if (phase(k1+j-1) == 1.) then
                  tod_corr(j) = -s_tot(k1+j-1) + self%scans(i)%d(det)%baseline2
               else
                  tod_corr(j) =  s_tot(k1+j-1) + self%scans(i)%d(det)%baseline1
               end if
            end if
         end do

         ! Apply ADC
         do j = 1, n
            tod_corr(j) = splint(self%adc(det)%p%F, real(tod_corr(j),dp))
         end do
         
         if (self%myid == 0 .and. det==2) then
            do j = 1, n
               write(58,*) k1+j-1, tod(k1+j-1), tod_corr(j), s_tot(k1+j-1), tod(k1+j-1)-tod_corr(j)
            end do
         end if
         
         ! Compute chisq
         do j = 1, n
            !write(*,*) i, j, tod_corr(j), s_tot(k1+j-1), sigma0(i)
            chisq_sub = chisq_sub + (tod_corr(j)-tod(k1+j-1))**2/sigma0(i)**2
         end do
      end do

      if (self%myid == 0 .and. det==2) close(58)

!!$      if (ntot > 0) write(*,*) 'adc chisq =', chisq_sub/ntot, ', myid = ', self%myid
      
      call mpi_reduce(chisq_sub, chisq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
      if (output_ndof) then
         call mpi_reduce(ntot, ndof_tot, 1, MPI_INTEGER8, MPI_SUM, 0, self%comm, ierr)
         if (self%myid == 0) ndof = ndof_tot
      end if

    end function chisq_adc_hfi2

  end subroutine sample_adc_and_baselines
  

end submodule comm_tod_hfi_mod
