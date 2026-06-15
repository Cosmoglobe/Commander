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
submodule (comm_tod_solat_mod) comm_tod_solat_smod
contains


  !**************************************************
  !             Constructor
  !**************************************************
  module function constructor_solat(cpar, id, id_abs, info, tod_type) result(c)
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
    class(comm_solat_tod),     pointer    :: c

    integer(i4b) :: i, j, k, nside_beam, lmax_beam, nmaps_beam, ierr
    logical(lgt) :: pol_beam
    character(len=6) :: pstring

    logical(lgt), dimension(:,:), allocatable :: correlations

    call timer%start(TOD_INIT, id_abs)
    
    ! Allocate object
    allocate(c)

    ! Initialize parameters needed for general initialization
    c%nhorn           = 1
    c%samprate_lowres = 18.  ! Lowres samprate in Hz;  10 times lower than the intrinsic lat rate for now    
    c%nmaps           = info%nmaps
    if (index(cpar%ds_tod_dets(id_abs), '.txt') /= 0) then
       c%ndet         = count_detectors(cpar%ds_tod_dets(id_abs)) !, cpar%datadir)
    else
       c%ndet         = num_tokens(cpar%ds_tod_dets(id_abs), ",")
    end if

    c%noise_psd_model = 'oof'       ! Not fitted parameters yet

    ! Initialize common parameters
    call c%tod_constructor(cpar, id, id_abs, info, tod_type)
    
    ! Initialize instrument-specific parameters

    c%compressed_tod    = .true.
    c%correct_sl        = .false.
    c%correct_orb       = .false.
    c%correct_S_crosstalk = .false.
    c%correct_N_crosstalk = .false.
    c%apply_inst_corr   = .false.
    c%orb_4pi_beam      = .false.
    c%symm_flags        = .false.
    c%equal_det_bp_beam = .true.
    c%sample_zodi     = cpar%sample_zodi .and. c%subtract_zodi ! Sample zodi parameters
    c%ntime           = 1
    !TODO: set the number of dark bolometers to be correct
    c%ndark           = 1
    c%n_cray_temps    = 3
    c%ndiode          = 1
    nmaps_beam        = 3
    pol_beam          = .true.
    c%nside_beam      = 128
    c%nside_pixhist   = 128
    c%sol_elong_range = cpar%zs_sol_elong
    
    ! Set up noise PSD type and priors
    c%freq            = cpar%ds_label(id_abs)
    c%n_xi            = 3
    !c%noise_psd_model = 'white'       ! Using white noise until we get better estimates of the actual noise PSD
    allocate(c%xi_n_P_uni(c%n_xi,2))
    allocate(c%xi_n_P_rms(c%n_xi))
    allocate(c%xi_n_nu_fit(c%n_xi,2))

    c%xi_n_P_uni(1,:)  = [10d0, 300d0]  ! Sigma0
    c%xi_n_P_uni(2,:)  = [1d0, 60d0]  ! fknee
    c%xi_n_P_uni(3,:)  = [-4.0d0, -1.0d0]   ! alpha
    !c%xi_n_P_uni(4,:)  = [ 0.5d0,  4.0d0]  ! fknee
    !c%xi_n_P_uni(5,:)  = [-1.5d0, -0.5d0]   ! alpha
    c%xi_n_nu_fit(1,:) = [0.1d0, 30d0] 
    c%xi_n_nu_fit(2,:) = [0.1d0, 30d0]
    c%xi_n_nu_fit(3,:) = [0.1d0, 30d0]
    !c%xi_n_nu_fit(4,:) = [0.001d0, 10d0]
    !c%xi_n_nu_fit(5,:) = [0.001d0, 10d0]
    c%xi_n_P_rms       = [10.d0, 0.1d0, 0.1d0] ! [sigma0, fknee, alpha]; sigma0 is not used
    

    !c%xi_n_P_rms      = [-1.d0] ! [sigma0]; sigma0 is not used
    
    ! Channel specific parameters
    c%chisq_threshold  = 1d30
    c%accept_threshold = 0.5d0
    c%correct_sl       = .false.


    ! Used for saving partial-sky maps to h5 file.
    c%cut_sky = .true.

    
    ! Get detector labels
    if (index(cpar%ds_tod_dets(id_abs), '.txt') /= 0) then
       call get_detectors(cpar%ds_tod_dets(id_abs), c%label)
    else
       call get_tokens(trim(adjustl(cpar%ds_tod_dets(id_abs))), ",", c%label)
    end if


    ! Read the actual TOD
    call c%read_tod(c%label)
    
    ! Initialize bandpass mean and proposal matrix
    call c%initialize_bp_covar(cpar%ds_tod_bp_init(id_abs))

    ! Construct lookup tables
    c%pixcache => comm_tod_pixcache(c%nside, c%nside_beam, c%nmaps, .false., c%equal_det_bp_beam)
    call c%precompute_lookups()
    
    ! Load the instrument file
    call c%load_instrument_file(c%nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

    ! Collect Sun velocities from all scans
    ! call c%collect_v_sun
    
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
    if (c%correct_N_crosstalk) then
       allocate(correlations(c%ndet, c%ndet))
       correlations = .true.    

       !  allocate(c%xtalk)
       c%xtalk => comm_crosstalk(correlations)
    end if

    ! Initialize dynamic mask
    c%dynmask => comm_dynmask(c, cpar)
    c%dynmask%apply_pixhist           = .true.
    c%dynmask%apply_solar_mask        = .false.
    c%dynmask%remove_isolated_samples = .true.
    
    call timer%stop(TOD_INIT, id_abs)

  end function constructor_solat

  !**************************************************
  !             Driver routine
  !**************************************************
  module subroutine process_solat_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
    !
    ! Routine that processes the lat time ordered data.
    ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
    ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms.
    ! Writes maps to disc in fits format
    !
    ! Arguments:
    ! ----------
    ! self:     pointer of comm_solat_tod class
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
    class(comm_solat_tod),                      intent(inout) :: self
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
    logical(lgt)        :: select_data, output_scanlist, output_zodi_comps, output_files
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
    type(hdf_file) :: tod_file


    ! Variables for jump detection
    real(sp),     allocatable, dimension(:,:,:) :: jump_calib
    real(sp),     allocatable, dimension(:)     :: tod_gapfill
    integer(i4b), allocatable, dimension(:)     :: jumps
    integer(i4b), allocatable, dimension(:,:)   :: offset_range, jumpflag_range
    real(sp),     allocatable, dimension(:)     :: offset_level
    type(comm_binmap)                           :: jump_map
    character(len=4)                            :: it_label

    call int2string(iter, ctext)
    call update_status(status, "tod_start"//ctext)
    call timer%start(TOD_TOT, self%band)
    
    ! Toggle optional operations
    sample_rel_bandpass   = size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
    sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses
    sample_rel_bandpass   = .false.
    if (.false.) then ! Debug
       ! Do data selection, then start sampling
       sample_gain           = .true.
       make_dyn_mask         = .false.
       sample_ncorr          = .true.
       sample_xi_n           = .false.
       select_data           = .false.
       output_files           = (iter == 10)
    else if (trim(self%init_from_HDF) == 'none') then
       ! Initialize slowly if not HDF init
       sample_gain           = iter  > 0 !.true.                 
       make_dyn_mask         = iter == 2
       sample_ncorr          = iter > 4 !.true.
       sample_xi_n           = iter > 5 
       select_data           = iter == 3 ! self%first_call  
       output_files           = (iter == 10)
    else
       ! Do data selection, then start sampling
       !sample_gain           = iter  > 2 !.true.                 
       sample_gain           = .false.
       make_dyn_mask         = (iter == 7)
       sample_ncorr          = iter  > 2 !.true.
       sample_xi_n           = iter > 3
       select_data           = (iter == 1) .or. (iter == 6)
       output_files           = (iter == 10)
    end if
    sample_zodi           = self%sample_zodi .and. self%subtract_zodi ! Sample zodi parameters
    output_zodi_comps     = self%output_zodi_comps .and. self%subtract_zodi ! Output zodi components
    output_scanlist       = mod(iter-1,10) == 0    ! only output scanlist every 10th iteration
    dec_wn                = 2 ! Decimation factor for sigma0; 2 corresponds to 45Hz

    !oper_default = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
    !    & SD_SKY,SD_BP,SD_ORB,SD_INST,SD_DARK,SD_NCORR])
    oper_default = get_sd_operation_code([SD_TOT, SD_BASE, SD_TOD, SD_IND, SD_NCORR, SD_SKY, SD_MASK, SD_GAIN, SD_JUMP])

    
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

    self%output_n_maps = 3

    call int2string(chain, ctext)
    call int2string(iter, samptext)
    call int2string(self%myid, myid_text)
    prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
    postfix = '_c' // ctext // '_k' // samptext // '.fits'

    ! Initialize index-based sky map and mask
    call self%pixcache%init_map_mask(map_in, self%bitmask, map_gain=map_gain)
    call update_status(status, "tod_cache"//ctext)


    call update_status(status, "tod_init")

    !------------------------------------
    ! Perform main sampling steps
    !------------------------------------

    if (self%first_call) call compute_tod_pixhist(self)

    if (sample_gain) then
       ! 'abscal': the global constant gain factor
       call sample_calibration(self, 'abscal', oper_default, handle)
       ! 'relcal': the gain factor that is constant in time but varying between detectors
       call sample_calibration(self, 'relcal', oper_default, handle)
       ! 'deltaG': the time-variable and detector-variable gain
       call sample_calibration(self, 'deltaG', oper_default, handle)
    end if
    
    call update_status(status, "initializing binmap")
    ! Prepare intermediate data structures
    !call binmap%init(self, .true., .false., nplus2=.false.)
    call binmap%init(self, .true., sample_rel_bandpass, nplus2=.false.)
    call jump_map%init(self, .true., sample_rel_bandpass, nplus2=.false.)  
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
          call init_scan_data(self, i, oper_default, TODMASK_NCORR, sd)
       end if
      call update_status(status, "made dynmask")


       if (sample_ncorr) then
          call sample_n_corr(self, sd, handle)
          if (sample_xi_n) then
             write(*,*) "Sampling noise psd for scan ", self%scanid(i)
             call sample_noise_psd(self, sd, handle, chaindir ,output_files=output_files)
             write(*,*) "Done sampling scan number   ", self%scanid(i)
          else
             call sample_noise_psd(self, sd, handle, chaindir, only_sigma0=.true., output_files=output_files)
          end if
       else
          call sample_n_corr(self, sd, handle, onlymono=.true.)
          call sample_noise_psd(self, sd, handle, chaindir, only_sigma0=.true., output_files=output_files)
       end if
      call update_status(status, "sampled ncorr")
       
       ! REMOVE JUMPS ----------------------------------
       allocate(jumps(sd%ntod))
       allocate(tod_gapfill(sd%ntod))


       sd%s_jump = 0.0

       do j=1, sd%ndet
        ! Throw away detectors that are more than 80% flagged. Also throw away detectors that don't have a partner. 
        if ((sum(sd%flag(:,j)) > 0.8*sd%ntod) .or. (.not. self%scans(i)%d(j)%accept)) then
           self%scans(i)%d(j)%accept = .false.
           cycle
        end if
   
          ! Retrieve offsets from previous run, if they exist
          if (allocated(self%scans(i)%d(j)%offset_range)) then
             call expand_offset_list(              &
                & self%scans(i)%d(j)%offset_range, &
                & self%scans(i)%d(j)%offset_level, &
                & sd%s_jump(:,j))
          else
             sd%s_jump(:,j) = 0.
          end if
          
          ! Retrieve jump flags from previous run, if they exist
          if (allocated(self%scans(i)%d(j)%jumpflag_range)) then
           call add_jumpflags(                     &
           & self%scans(i)%d(j)%jumpflag_range, &
           & sd%flag(:,j))
        end if
        
        ! Scanning for jumps
        if (.true.) then
           call jump_scan(                                 &
           & sd%tod(:,j) - sd%s_sky(:,j,0,1) - sd%s_jump(:,j) - sd%n_corr(:,j), &
           & sd%flag(:,j),                              &
           & jumps,                                     &
           & offset_range,                              &
           & offset_level,                              &
           & handle,                                    &
           & jumpflag_range,                            &
           & it_label,                                  &
           & chaindir,                                  &
           & .false.)
           
             ! Add offsets to persistent list
             if (.not. allocated(self%scans(i)%d(j)%offset_range)) then
                allocate(self%scans(i)%d(j)%offset_range(size(offset_level),2))
                allocate(self%scans(i)%d(j)%offset_level(size(offset_level)))

                self%scans(i)%d(j)%offset_range = offset_range
                self%scans(i)%d(j)%offset_level = offset_level
             else
                call update_offset_list(              &
                   & offset_range,                    &
                   & offset_level,                    &
                   & self%scans(i)%d(j)%offset_range, &
                   & self%scans(i)%d(j)%offset_level)
             end if

             ! Add jump flags to persistent list
             if (allocated(jumpflag_range)) then
                if (.not. allocated(self%scans(i)%d(j)%jumpflag_range)) then
                   allocate(self%scans(i)%d(j)%jumpflag_range(size(jumpflag_range)/2,2))
                   self%scans(i)%d(j)%jumpflag_range = jumpflag_range
                else
                   call update_jumpflag(jumpflag_range, self%scans(i)%d(j)%jumpflag_range)
                end if
             end if

             call expand_offset_list(                &
                 & self%scans(i)%d(j)%offset_range,  &
                 & self%scans(i)%d(j)%offset_level,  & 
                 & sd%s_jump(:,j))
          end if


          call gap_fill_linear(           &
             & sd%tod(:,j) - sd%s_jump(:,j), &
             & sd%flag(:,j),              &
             & tod_gapfill,               &
             & handle,                    &
             & .true.)


          if (allocated(offset_range))   deallocate(offset_range)
          if (allocated(offset_level))   deallocate(offset_level)
          if (allocated(jumpflag_range)) deallocate(jumpflag_range)

       end do

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

       allocate(jump_calib(1, sd%ntod, sd%ndet))
       jump_calib(1,:,:) = sd%s_jump 
       call update_status(status, "computed calibrated")

       
       ! Bin TOD
       call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib, binmap)
       call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, jump_calib, jump_map) 
       call update_status(status, "binned TOD")
       
       ! Update scan list
       call wall_time(t2)
       self%scans(i)%proctime   = self%scans(i)%proctime   + t2-t1
       self%scans(i)%n_proctime = self%scans(i)%n_proctime + 1
       if (output_scanlist) then
          write(slist(i),*) self%scanid(i), '"',trim(self%hdfname(i)), &
               & '"', real(self%scans(i)%proctime/self%scans(i)%n_proctime,sp),&
               & real(self%spinaxis(i,:),sp)
       end if

       ! For debugging: write TOD to hdf
       if (output_files) then
          !if (self%scanid(i) == 915) then 
             !print *, self%scanid(i)
             call int2string(self%scanid(i), scantext)
             call open_hdf_file(trim(chaindir)//'/res_'//trim(self%label(1))//'_'//scantext//'.h5', tod_file, 'w')
             call write_hdf(tod_file, '/tod', sd%tod)
             call write_hdf(tod_file, '/pix', sd%pix(:,:,1))
             call write_hdf(tod_file, '/flag', sd%flag)
             call write_hdf(tod_file, '/n_corr', sd%n_corr)
             call write_hdf(tod_file, '/s_jump', sd%s_jump)
             call close_hdf_file(tod_file)
          !end if
       end if

       ! Clean up
       call dealloc_scan_data(sd)
       deallocate(d_calib)
       deallocate(jumps, tod_gapfill, jump_calib)
       call update_status(status, "cleaned up")

    end do


    call timer%start(TOD_WAIT, self%band)
    call mpi_barrier(self%comm, ierr) ! Improve timing information
    call timer%stop(TOD_WAIT, self%band)
    call update_status(status, "tod_postloop"//ctext)
!!$       call mpi_finalize(ierr)
!!$       stop
    
    !if (select_data) then
    !   ! Remove data based on a gliding RMS window cut for each of the listed
    !   ! criteria
    !   !                           half-window  [chisq, sigma0, fknee, alpha, base, base1, base2]
    !   call remove_tod_outliers(self, 100,      [100000.,    100000.,     100000.,    100000.,    100000.,   100000.,    100000.   ])
    !      
    !   do i = 1, self%nscan
    !      do j = 1, self%ndet
    !         if (.not. self%scans(i)%d(j)%accept) self%scans(i)%d(j)%nsamp_unmasked = 0
    !      end do
    !   end do
    !end if
    !call update_status(status, "tod_outlier"//ctext)

    
    if (self%myid == 0) write(*,*) '   --> Finalizing maps, bp'
    if (make_dyn_mask) then
       call self%dynmask%report
       call update_status(status, "tod_dynstats"//ctext)
    end if
    
    ! Output latest scan list with new timing information
    if (output_scanlist) call self%output_scan_list(slist)

    ! Solve for maps
    call synchronize_binmap(binmap, self)
    call synchronize_binmap(jump_map, self) 
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
    call update_status(status, "tod_binmap_jump"//ctext)
    call finalize_binned_map(self, jump_map, rms_out, 1.d0) 
    call update_status(status, "tod_binmap2"//ctext)

    ! Sample bandpass parameters
    if (sample_rel_bandpass .or. sample_abs_bandpass) then
       call sample_bp(self, handle, chisq_S, delta)
       self%bp_delta = delta(:,:,1)
    end if
   
    ! Output maps to disk
    call timer%start(TOD_WRITE)
    call map_out%writeFITS(trim(prefix)//'map'//trim(postfix), cut_sky=self%cut_sky)
    call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix), cut_sky=self%cut_sky)
    if (self%output_n_maps > 1) call binmap%outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix), cut_sky=self%cut_sky)
    if (self%output_n_maps > 2) call binmap%outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix), cut_sky=self%cut_sky)
    if (self%output_n_maps > 3) call binmap%outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix), cut_sky=self%cut_sky)
    if (self%output_n_maps > 4) call binmap%outmaps(5)%p%writeFITS(trim(prefix)//'orb'//trim(postfix), cut_sky=self%cut_sky)
    if (self%output_n_maps > 5) call binmap%outmaps(6)%p%writeFITS(trim(prefix)//'sl'//trim(postfix), cut_sky=self%cut_sky)
    call jump_map%outmaps(1)%p%writeFITS(trim(prefix)//'jumps'//trim(postfix), cut_sky=self%cut_sky) 
    call timer%stop(TOD_WRITE)
    call update_status(status, "tod_binmap3"//ctext)

    ! Clean up
    call binmap%dealloc()
    call jump_map%dealloc() 
    call update_status(status, "tod_binmap4, deallocation "//ctext)
    if (allocated(slist)) deallocate(slist)
    if (sample_abs_bandpass .or. sample_rel_bandpass) deallocate(chisq_S)
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
    
  end subroutine process_solat_tod

  module subroutine load_instrument_solat(self, instfile, band)
    !
    ! Reads the lat specific fields from the instrument file
    ! Implements comm_tod_mod::load_instrument_inst
    !
    ! Arguments:
    !
    ! self : comm_solat_tod
    !    the lat tod object (this class)
    ! file : hdf_file
    !    the open file handle for the instrument file
    ! band : int
    !    the index of the current detector
    ! 
    ! Returns : None
    implicit none
    class(comm_solat_tod),                 intent(inout) :: self
    type(hdf_file),                      intent(in)    :: instfile
    integer(i4b),                        intent(in)    :: band

    call read_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'polEff', self%pol_eff(band))
    !self%pol_eff(band) = self%pol_eff(band) * 0.01d0 ! Stored as percentage in the instrument file for now
    self%pol_eff(band) = 1.0

  end subroutine load_instrument_solat


  module subroutine read_tod_inst_solat(self, file)
    ! 
    ! Reads lat-specific common fields from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_solat_tod)
    !           lat-specific TOD object
    ! file:     derived type (hdf_file)
    !           Already open HDF file handle; only root includes this
    !
    ! Returns
    ! ----------
    ! None, but updates self
    !
    implicit none
    class(comm_solat_tod),                 intent(inout)          :: self
    type(hdf_file),                      intent(in),   optional :: file
  end subroutine read_tod_inst_solat
  
  module subroutine read_scan_inst_solat(self, file, slabel, detlabels, scan)
    ! 
    ! Reads lat-specific scan information from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_solat_tod)
    !           lat-specific TOD object
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
    class(comm_solat_tod),                 intent(in)    :: self
    type(hdf_file),                      intent(in)    :: file
    character(len=*),                    intent(in)    :: slabel
    character(len=*), dimension(:),      intent(in)    :: detlabels
    class(comm_scan),                    intent(inout) :: scan
  end subroutine read_scan_inst_solat

  module subroutine initHDF_solat(self, chainfile, path)
    ! 
    ! Initializes lat-specific TOD parameters from existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_solat_tod)
    !           lat-specific TOD object
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
    class(comm_solat_tod),                 intent(inout)  :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path


  end subroutine initHDF_solat
  
  module subroutine dumpToHDF_solat(self, chainfile, path)
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
    class(comm_solat_tod),                 intent(in)     :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path


  end subroutine dumpToHDF_solat


  module subroutine construct_corrtemp_solat(self, sd, det)
    !  Construct an lat instrument-specific correction template; for now contains 1Hz template only
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
    class(comm_solat_tod),  intent(in)           :: self
    class(comm_scandata), intent(inout)        :: sd
    integer(i4b),         intent(in), optional :: det

    integer(i4b) :: i, j, k, nbin, b
    real(dp)     :: dt, t_tot, t

    sd%s_inst = 0.

  end subroutine construct_corrtemp_solat


end submodule comm_tod_solat_smod
