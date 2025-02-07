module comm_tod_driver_mod
  use comm_tod_mod
  use comm_tod_gain_mod
  use comm_tod_noise_mod
  use comm_tod_pointing_mod
  use comm_tod_bandpass_mod
  use comm_tod_orbdipole_mod
  use comm_tod_simulations_mod
  use comm_tod_mapmaking_mod
  use comm_tod_jump_mod
  use comm_tod_adc_mod
  use comm_zodi_mod
  use comm_tod_cray_mod
  use comm_shared_arr_mod
  use comm_huffman_mod
  use comm_4d_map_mod
  !use comm_zodi_samp_mod
  use omp_lib
  implicit none

  ! Class for uncompressed data for a given scan
  type :: comm_scandata
     integer(i4b) :: ntod, ndet, nhorn, ndelta
     real(sp),     allocatable, dimension(:,:)     :: tod           ! Raw data
     real(sp),     allocatable, dimension(:,:)     :: n_corr        ! Correlated noise in V
     real(sp),     allocatable, dimension(:,:)     :: s_sl          ! Sidelobe correction
     real(sp),     allocatable, dimension(:,:)     :: s_sky         ! Stationary sky signal
     real(sp),     allocatable, dimension(:,:,:)   :: s_sky_prop    ! Stationary sky signal proposal for bandpass sampling
     real(sp),     allocatable, dimension(:,:)     :: s_orb         ! Orbital dipole
     real(sp),     allocatable, dimension(:,:)     :: s_mono        ! Detector monopole correction 
     real(sp),     allocatable, dimension(:,:)     :: s_calib       ! Custom calibrator
     real(sp),     allocatable, dimension(:,:)     :: s_calibA      ! Custom calibrator
     real(sp),     allocatable, dimension(:,:)     :: s_calibB      ! Custom calibrator
     real(sp),     allocatable, dimension(:,:)     :: s_bp          ! Bandpass correction
     real(sp),     allocatable, dimension(:,:,:)   :: s_bp_prop     ! Bandpass correction proposal     
     real(sp),     allocatable, dimension(:,:)     :: s_zodi        ! Zodiacal emission
     real(sp),     allocatable, dimension(:,:,:)   :: s_zodi_scat   ! Scattered sunlight contribution to zodi  
     real(sp),     allocatable, dimension(:,:,:)   :: s_zodi_therm  ! Thermal zodiacal emission
     real(sp),     allocatable, dimension(:,:)     :: s_inst        ! Instrument-specific correction template
     real(sp),     allocatable, dimension(:,:)     :: s_tot         ! Total signal
     real(sp),     allocatable, dimension(:,:)     :: s_gain        ! Absolute calibrator
     real(sp),     allocatable, dimension(:,:)     :: mask          ! TOD mask (flags + main processing mask)
     real(sp),     allocatable, dimension(:,:)     :: mask2         ! Small TOD mask, for bandpass sampling
     real(sp),     allocatable, dimension(:,:)     :: mask_zodi     ! Mask for sampling zodi
     integer(i4b), allocatable, dimension(:,:,:)   :: pix           ! Discretized pointing 
     integer(i4b), allocatable, dimension(:,:,:)   :: psi           ! Discretized polarization angle
     integer(i4b), allocatable, dimension(:,:)     :: flag          ! Quality flags
     real(sp),     allocatable, dimension(:,:)     :: s_totA        ! Total signal, horn A (differential only)
     real(sp),     allocatable, dimension(:,:)     :: s_totB        ! Total signal, horn B (differential only)
     real(sp),     allocatable, dimension(:,:)     :: s_gainA        ! Total signal, horn A (differential only)
     real(sp),     allocatable, dimension(:,:)     :: s_gainB        ! Total signal, horn B (differential only)
     real(sp),     allocatable, dimension(:,:)     :: s_orbA        ! Orbital signal, horn A (differential only)
     real(sp),     allocatable, dimension(:,:)     :: s_orbB        ! Orbital signal, horn B (differential only)
     integer(i4b) :: band                                           ! Band ID
   contains
     procedure  :: init_singlehorn   => init_scan_data_singlehorn
     procedure  :: init_differential => init_scan_data_differential
     procedure  :: dealloc           => dealloc_scan_data
  end type comm_scandata


contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Scan data routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_scan_data_singlehorn(sd, tod, scan, map_sky, map_gain, procmask, procmask2, procmask_zodi, &
       & init_s_bp, init_s_bp_prop, init_s_sky_prop, skip_nonlin)
    implicit none
    class(comm_scandata),                      intent(inout)          :: sd    
    class(comm_tod),                           intent(inout)          :: tod
    integer(i4b),                              intent(in)             :: scan
    real(sp),          dimension(1:,1:,0:,1:), intent(in)             :: map_sky
    real(sp),          dimension(1:,1:,0:,1:), intent(in)             :: map_gain
    real(sp),          dimension(0:),          intent(in)             :: procmask
    real(sp),          dimension(0:),          intent(in)             :: procmask2
    real(sp),          dimension(0:),          intent(in),   optional :: procmask_zodi
    logical(lgt),                              intent(in),   optional :: init_s_bp
    logical(lgt),                              intent(in),   optional :: init_s_bp_prop
    logical(lgt),                              intent(in),   optional :: init_s_sky_prop
    logical(lgt),                              intent(in),   optional :: skip_nonlin

    integer(i4b) :: i, j, k, ndelta
    logical(lgt) :: init_s_bp_, init_s_bp_prop_, init_s_sky_prop_, skip_nonlin_




  end subroutine init_scan_data_singlehorn


  subroutine init_scan_data_differential(sd, tod, scan, map_sky, map_gain, procmask, procmask2, &
       & init_s_bp, init_s_bp_prop, init_s_sky_prop, polang)
    implicit none
    class(comm_scandata),                      intent(inout)          :: sd    
    class(comm_tod),                           intent(inout)          :: tod
    integer(i4b),                              intent(in)             :: scan
    real(sp),          dimension(1:,1:,0:,1:), intent(in)             :: map_sky
    real(sp),          dimension(1:,1:,0:,1:), intent(in)             :: map_gain
    real(sp),          dimension(0:),          intent(in)             :: procmask
    real(sp),          dimension(0:),          intent(in)             :: procmask2
    logical(lgt),                              intent(in),   optional :: init_s_bp
    logical(lgt),                              intent(in),   optional :: init_s_bp_prop
    logical(lgt),                              intent(in),   optional :: init_s_sky_prop
    real(dp),                                  intent(in),   optional :: polang

    integer(i4b) :: j, k, ndelta
    logical(lgt) :: init_s_bp_, init_s_bp_prop_, init_s_sky_prop_
    real(sp),     allocatable, dimension(:,:)     :: s_bufA, s_bufB, s_buf2A, s_buf2B      ! Buffer
    !
    ! Note that procmask should be larger, cover the residuals in the galactic
    ! plane, so that it is not used for calibration and correlated noise. 
    ! procmask2 should be smaller and be used only for mapmaking.
    !


  end subroutine init_scan_data_differential


  subroutine dealloc_scan_data(sd)
    implicit none
    class(comm_scandata), intent(inout) :: sd    

    sd%ntod = -1; sd%ndet = -1; sd%nhorn = -1

    ! Deallocate data structures
    deallocate(sd%tod, sd%n_corr, sd%s_sl, sd%s_sky)
    deallocate(sd%s_orb, sd%s_tot, sd%mask)
    deallocate(sd%pix, sd%psi, sd%flag)
    if (allocated(sd%s_sky_prop))    deallocate(sd%s_sky_prop)
    if (allocated(sd%s_bp_prop))     deallocate(sd%s_bp_prop)
    if (allocated(sd%s_bp))          deallocate(sd%s_bp)
    if (allocated(sd%s_mono))        deallocate(sd%s_mono)
    if (allocated(sd%mask2))         deallocate(sd%mask2)
    if (allocated(sd%s_zodi))        deallocate(sd%s_zodi)
    if (allocated(sd%s_zodi_scat))   deallocate(sd%s_zodi_scat)
    if (allocated(sd%s_zodi_therm))  deallocate(sd%s_zodi_therm)
    if (allocated(sd%mask_zodi))     deallocate(sd%mask_zodi)
    if (allocated(sd%s_totA))        deallocate(sd%s_totA)
    if (allocated(sd%s_totB))        deallocate(sd%s_totB)
    if (allocated(sd%s_inst))        deallocate(sd%s_inst)
    if (allocated(sd%s_orbA))        deallocate(sd%s_orbA)
    if (allocated(sd%s_orbB))        deallocate(sd%s_orbB)
    if (allocated(sd%s_gain))        deallocate(sd%s_gain)
    if (allocated(sd%s_gainA))       deallocate(sd%s_gainA)
    if (allocated(sd%s_gainB))       deallocate(sd%s_gainB)

  end subroutine dealloc_scan_data


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Sampling drivers etc.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sample_calibration(tod, mode, handle, map_sky, map_gain, procmask, procmask2, &
      & polang, smooth)
    !
    !   Sample calibration modes
    !   Supported modes = {abscal, relcal, deltaG, imbal}
    !
    !   Subroutine that implements gain and horn imbalance sampling, following
    !   the formalism of Gjerlow et al. 2021
    !
    !   Arguments:
    !   ----------
    !   tod:      comm_tod derived type
    !             contains TOD-specific information
    !   mode:     character array
    !             specifies sampling mode. currently supports absolute calibration,
    !             relative calibration, time-variable calibration, and horn
    !             imbalance.
    !   handle:   planck_rng derived type
    !             Healpix definition for random number generation
    !             so that the same sequence can be resumed later on from that same point
    !   map_sky:
    !
    implicit none
    class(comm_tod),                              intent(inout) :: tod
    character(len=*),                             intent(in)    :: mode
    type(planck_rng),                             intent(inout) :: handle
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:,1:,1:,1:),   intent(in)    :: map_gain
    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2
    real(dp),                                  intent(in),   optional :: polang
    logical(lgt),                              intent(in),   optional :: smooth

    integer(i4b) :: i, j, ext(2), ierr, timer_id
    real(dp)     :: t1, t2
    real(dp), allocatable, dimension(:)   :: A, b
    real(sp), allocatable, dimension(:,:) :: s_invsqrtN, mask_lowres, s_buf
    real(dp), allocatable, dimension(:,:) :: dipole_mod
    type(comm_scandata) :: sd
    logical(lgt) :: smooth_


    smooth_ = .true.
    if (present(smooth))  smooth_=smooth

    if (tod%myid == 0) write(*,*) '|    --> Sampling calibration, mode = ', trim(mode)

    if (trim(mode) == 'abscal') then
       timer_id = TOD_ABSCAL
    else if (trim(mode) == 'relcal') then
       timer_id = TOD_RELCAL
    else if (trim(mode) == 'imbal') then
       timer_id = TOD_IMBAL
    else if (trim(mode) == 'deltaG') then
       timer_id = TOD_DELTAG
    end if


    if (trim(mode) == 'abscal' .or. trim(mode) == 'relcal' .or. trim(mode) == 'imbal') then
       allocate(A(tod%ndet), b(tod%ndet))
       A = 0.d0; b = 0.d0
    else if (trim(mode) == 'deltaG') then
       allocate(dipole_mod(tod%nscan_tot, tod%ndet))
       dipole_mod = 0.d0
    else
       write(*,*) 'Unsupported sampling mode!'
       stop
    end if



    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ![Debug] if (tod%myid == 0) write(*,*) '|    --> Preparing data ' !on, mode = ', trim(mode)
       ! Prepare data
       if (tod%nhorn == 1) then
          call sd%init_singlehorn(tod, i, map_sky, map_gain, procmask, procmask2)
       else
          call sd%init_differential(tod, i, map_sky, map_gain, procmask, procmask2, polang=polang)
       end if

       ![Debug] if (tod%myid == 0) write(*,*) '|    --> Setup filtered calibration signal'! m(mode)
       ! Set up filtered calibration signal, conditional contribution and mask
       call tod%downsample_tod(sd%s_orb(:,1), ext)
       allocate(s_invsqrtN(ext(1):ext(2), tod%ndet))      ! s * invN
       allocate(s_buf(sd%ntod, sd%ndet))
       allocate(mask_lowres(ext(1):ext(2), tod%ndet))
       do j = 1, tod%ndet
          if (.not. tod%scans(i)%d(j)%accept) cycle
          call tod%downsample_tod(sd%mask(:,j), ext, mask_lowres(:,j), threshold=0.9)
          if (trim(mode) == 'abscal') then
             if (trim(tod%abscal_comps) == 'orbital') then
               call tod%downsample_tod(sd%s_orb(:,j), ext, s_invsqrtN(:,j))
             else if (trim(tod%abscal_comps) == 'full') then
                ! Calibrator = total signal
                s_buf(:,j) = sd%s_tot(:,j)
                call fill_all_masked(s_buf(:,j), sd%mask(:,j), sd%ntod, .false., real(tod%scans(i)%d(j)%N_psd%sigma0, sp), handle, tod%scans(i)%chunk_num)
                call tod%downsample_tod(s_buf(:,j), ext, s_invsqrtN(:,j))
             else
                s_buf(:,j) = sd%s_gain(:,j)
                call fill_all_masked(s_buf(:,j), sd%mask(:,j), sd%ntod, .false., real(tod%scans(i)%d(j)%N_psd%sigma0, sp), handle, tod%scans(i)%chunk_num)
                call tod%downsample_tod(s_buf(:,j), ext, s_invsqrtN(:,j))
             end if
          else if (trim(mode) == 'imbal' .and. tod%nhorn == 2) then
             ! Calibrator = common mode signal
             s_buf(:,j) = tod%scans(i)%d(j)%gain*(sd%s_totA(:,j) + sd%s_totB(:,j))
             call fill_all_masked(s_buf(:,j), sd%mask(:,j), sd%ntod, .false., &
               & real(tod%scans(i)%d(j)%N_psd%sigma0, sp), handle, tod%scans(i)%chunk_num)
             call tod%downsample_tod(s_buf(:,j), ext, s_invsqrtN(:,j))
          else
             ! Calibrator = total signal
             s_buf(:,j) = sd%s_tot(:,j)
             call fill_all_masked(s_buf(:,j), sd%mask(:,j), sd%ntod, .false., real(tod%scans(i)%d(j)%N_psd%sigma0, sp), handle, tod%scans(i)%chunk_num)
             call tod%downsample_tod(s_buf(:,j), ext, s_invsqrtN(:,j))
          end if
       end do
       ! [Debug] if (tod%myid == 0) write(*,*) '|    --> Passed the loop with downsampel tod'!(mode)
       call multiply_inv_N(tod, i, s_invsqrtN, sampfreq=tod%samprate_lowres, pow=0.5d0)

       if (trim(mode) == 'abscal' .or. trim(mode) == 'relcal' .or. trim(mode) == 'imbal') then
          ! Constant gain terms; accumulate contribution from this scan
          do j = 1, tod%ndet
             if (.not. tod%scans(i)%d(j)%accept) cycle
             if (trim(mode) == 'abscal') then
                if (trim(tod%abscal_comps) == 'orbital') then
                  s_buf(:,j) = real(tod%gain0(0),sp) * (sd%s_tot(:,j) - sd%s_orb(:,j)) + &
                       & real(tod%gain0(j) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j)
                else if (trim(tod%abscal_comps) == 'full') then
                  s_buf(:,j) = real(tod%gain0(j) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j)
                else
                  s_buf(:,j) = real(tod%gain0(0),sp) * (sd%s_tot(:,j) - sd%s_gain(:,j)) + &
                       & real(tod%gain0(j) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j)
                end if
             else if (trim(mode) == 'relcal') then
                s_buf(:,j) = real(tod%gain0(0) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j)
             else if (trim(mode) == 'imbal') then
                s_buf(:,j) = tod%scans(i)%d(j)%gain * (sd%s_totA(:,j) - sd%s_totB(:,j))
             end if
          end do
          call accumulate_abscal(tod, i, sd%mask, s_buf, s_invsqrtN, A, b, handle, &
              & out=.true., mask_lowres=mask_lowres, tod_arr=sd%tod)
       else
          ! Time-variable gain terms
            call calculate_gain_mean_std_per_scan(tod, i, s_invsqrtN, sd%mask, sd%s_tot, &
              & handle, mask_lowres=mask_lowres, tod_arr=sd%tod)
          do j = 1, tod%ndet
             if (.not. tod%scans(i)%d(j)%accept) cycle
          end do
       end if

       ![Debug] if (tod%myid == 0) write(*,*) '|    --> Passed another loop'!(mode)
       ! Clean up
       call wall_time(t2)
       tod%scans(i)%proctime   = tod%scans(i)%proctime   + t2-t1
       tod%scans(i)%n_proctime = tod%scans(i)%n_proctime + 1
       call dealloc_scan_data(sd)
       deallocate(s_invsqrtN, s_buf, mask_lowres)
    end do


    ! Perform sampling operations
    if (trim(mode) == 'abscal') then
       call sample_abscal_from_orbital(tod, handle, A, b)
    else if (trim(mode) == 'relcal') then
       call sample_relcal(tod, handle, A, b)
    else if (trim(mode) == 'deltaG') then
       call mpi_allreduce(mpi_in_place, dipole_mod, size(dipole_mod), MPI_DOUBLE_PRECISION, MPI_SUM, tod%info%comm, ierr)
       call sample_smooth_gain(tod, handle, dipole_mod, smooth_)
    else if (trim(mode) == 'imbal') then
       call sample_imbal_cal(tod, handle, A, b)
    end if

    ! Clean up
    if (allocated(A))          deallocate(A)
    if (allocated(b))          deallocate(b)
    if (allocated(dipole_mod)) deallocate(dipole_mod)


  end subroutine sample_calibration


  ! Sample baseline
  subroutine sample_baseline(tod, handle, map_sky, map_gain, procmask, procmask2)
    implicit none
    class(comm_tod),                              intent(inout) :: tod
    type(planck_rng),                             intent(inout) :: handle
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_gain
    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2

    integer(i4b) :: i, j
    real(dp)     :: t1, t2
    type(comm_scandata) :: sd

    if (tod%myid == 0) write(*,*) '   --> Sampling baseline'

    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       if (tod%nhorn == 1) then
          call sd%init_singlehorn(tod, i, map_sky, map_gain, procmask, procmask2)
       else
          call init_scan_data_differential(sd, tod, i, map_sky, map_gain, procmask, procmask2)
       end if

       do j = 1, tod%ndet
          ! Return the data to its raw state
          sd%tod(:,j) = sd%tod(:,j) + tod%scans(i)%d(j)%baseline

          ! Estimate the baseline and sample it if requested
          tod%scans(i)%d(j)%baseline =sum((sd%tod(:,j) - tod%scans(i)%d(j)%gain*sd%s_tot(:,j)) &
            & *sd%mask(:,j))/sum(sd%mask(:,j))
          if (trim(tod%operation) == 'sample') then
            tod%scans(i)%d(j)%baseline = tod%scans(i)%d(j)%baseline &
             &  + rand_gauss(handle)/sqrt(sum(sd%mask(:,j)*tod%scans(i)%d(j)%N_psd%sigma0**2))
          end if
       end do

       ! Clean up
       call wall_time(t2)
       tod%scans(i)%proctime   = tod%scans(i)%proctime   + t2-t1
       tod%scans(i)%n_proctime = tod%scans(i)%n_proctime + 1
       call dealloc_scan_data(sd)
    end do

  end subroutine sample_baseline

  subroutine remove_bad_data(tod, scan, flag)
    !
    !   Perform data selection on TOD object
    !
    !   Arguments:
    !   ----------
    !   tod:      comm_tod derived type
    !             contains TOD-specific information. Bad data are removed by 
    !             setting scan%det%accept = .false.
    !   scan:     int (scalar)
    !             Local scan ID for the current core 
    !   flag:     int (ntod x ndet array)
    !             Array with data quality flags
    !
    implicit none
    class(comm_tod),                   intent(inout) :: tod
    integer(i4b),    dimension(1:,1:), intent(in)    :: flag
    integer(i4b),                      intent(in)    :: scan

    integer(i4b) :: j, ntod, ndet
    
    ntod = size(flag,1)
    ndet = size(flag,2)
    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (count(iand(flag(:,j),tod%flag0) .ne. 0) > tod%accept_threshold*ntod) then    ! Discard scans with less than 20% good data
          tod%scans(scan)%d(j)%accept = .false.
          write(*, fmt='(a, i4, a, i8, a, i8)') ' | Reject scan = ', &
            & tod%scanid(scan), ': ', count(iand(flag(:,j),tod%flag0) .ne. 0), &
            &  ' flagged data out of', ntod
       else if (abs(tod%scans(scan)%d(j)%chisq) > tod%chisq_threshold .or. &  ! Discard scans with high chisq or NaNs
            & isNaN(tod%scans(scan)%d(j)%chisq)) then
          write(*,fmt='(a,i8,i5,a,f12.1)') ' | Reject scan, det = ', &
               & tod%scanid(scan), j, ', chisq = ', tod%scans(scan)%d(j)%chisq
          tod%scans(scan)%d(j)%accept = .false.
       end if
    end do
    if (tod%symm_flags) then
       do j = 1, ndet
           if (.not. tod%scans(scan)%d(j)%accept .and. tod%partner(j) >= 0) tod%scans(scan)%d(tod%partner(j))%accept = .false.
       end do
    end if
  end subroutine remove_bad_data

  subroutine compute_chisq_abs_bp(tod, scan, sd, chisq)
    implicit none
    class(comm_tod),                       intent(inout) :: tod
    integer(i4b),                          intent(in)    :: scan
    type(comm_scandata),                   intent(in)    :: sd
    real(dp),            dimension(:,:),   intent(inout) :: chisq

    integer(i4b) :: j, k
    real(sp), allocatable, dimension(:,:) :: s_buf

    allocate(s_buf(sd%ntod,sd%ndet))
    do j = 1, tod%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       s_buf(:,j) =  sd%s_sl(:,j) + sd%s_orb(:,j)
       call tod%compute_tod_chisq(scan, j, sd%mask2(:,j), sd%s_sky(:,j), &
            & s_buf(:,j), sd%n_corr(:,j), sd%tod(:,j), absbp=.true.)
       chisq(j,1) = chisq(j,1) + tod%scans(scan)%d(j)%chisq_prop
       do k = 2, tod%n_bp_prop+1
          call tod%compute_tod_chisq(scan, j, sd%mask2(:,j), sd%s_sky_prop(:,j,k), &
               & s_buf(:,j), sd%n_corr(:,j), sd%tod(:,j), absbp=.true.)
          chisq(j,k) = chisq(j,k) + tod%scans(scan)%d(j)%chisq_prop
       end do
    end do
    deallocate(s_buf)

  end subroutine compute_chisq_abs_bp

  subroutine compute_calibrated_data(tod, scan, sd, d_calib, jump_template)
    !
    !  gets calibrated timestreams
    !
    !  Arguments:
    !  ----------
    !  tod: comm_tod object
    !
    !  scan: integer
    !     integer label for scan
    !  sd:  comm_scandata object
    !  jump_template:  baseline that traces jumping tod level
    !
    !  Returns:
    !  --------
    !  d_calib: real(sp) array
    !     nout x ndet x ntod array of calibrated timestreams
    !
    !  d_calib(1,:,:) - best estimate of calibrated data, with all known
    !    calibrations applied
    !  d_calib(2,:,:) - calibrated TOD with expected sky signal subtracted,
    !    i.e., residual
    !  d_calib(3,:,:) - correlated noise, mean subtracted, in temperature
    !    units
    !  d_calib(4,:,:) - bandpass difference contribution
    !  d_calib(5,:,:) - orbital dipole
    !  d_calib(6,:,:) - sidelobe
    !  d_calib(7,:,:) - zodiacal light emission
    !  d_calib(8,:,:) - instrument correction
    !  d_calib(9 - 9 + n_zodi_comps,:,:) - zodiacal light components
    implicit none
    class(comm_tod),                       intent(in)   :: tod
    integer(i4b),                          intent(in)   :: scan
    type(comm_scandata),                   intent(in)   :: sd
    real(sp),            dimension(:,:,:), intent(out)  :: d_calib
    real(sp), dimension(:,:), intent(in), optional      :: jump_template
    integer(i4b) :: i, j, k, nout
    real(dp)     :: inv_gain
   !  write(*, *) "s_bp:", sd%s_sky(:,1)
    nout = size(d_calib,1)
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       inv_gain = 1.0 / tod%scans(scan)%d(j)%gain
       if (tod%compressed_tod) then
        d_calib(1,:,j) = (sd%tod(:,j) - sd%n_corr(:,j)) &
          & * inv_gain - sd%s_tot(:,j) + sd%s_sky(:,j) - sd%s_bp(:,j)
       else
        d_calib(1,:,j) = (tod%scans(scan)%d(j)%tod - sd%n_corr(:,j)) &
          & * inv_gain - sd%s_tot(:,j) + sd%s_sky(:,j) - sd%s_bp(:,j)
       end if

       if (present(jump_template)) d_calib(1,:,j) = d_calib(1,:,j) - jump_template(:,j) * inv_gain
       if (tod%output_n_maps > 1) d_calib(2,:,j) = d_calib(1,:,j) - sd%s_sky(:,j) + sd%s_bp(:,j)              ! residual
       if (tod%output_n_maps > 2) d_calib(3,:,j) = sd%n_corr(:,j) * inv_gain  ! ncorr
       !if (tod%output_n_maps > 2) d_calib(3,:,j) = (sd%n_corr(:,j) - sum(real(sd%n_corr(:,j),dp)/sd%ntod)) * inv_gain  ! ncorr
       if (tod%output_n_maps > 3) d_calib(4,:,j) = sd%s_bp(:,j)                                               ! bandpass
       if (tod%output_n_maps > 4) d_calib(5,:,j) = sd%s_orb(:,j)                                              ! orbital dipole
       if (tod%output_n_maps > 5) d_calib(6,:,j) = sd%s_sl(:,j)          
       if ((tod%output_n_maps > 6) .and. allocated(sd%s_zodi)) d_calib(7,:,j) = sd%s_zodi(:,j) ! zodi
       if ((tod%output_n_maps > 7) .and. allocated(sd%s_inst)) d_calib(8,:,j) = (sd%s_inst(:,j) - sum(real(sd%s_inst(:,j),dp)/sd%ntod)) * inv_gain  ! instrument specific
       if ((tod%output_n_maps > 8) .and. allocated(sd%s_zodi_scat) .and. allocated(sd%s_zodi_therm)) then
          do i = 1, size(sd%s_zodi_therm, dim=2)
             !write(*,*) 'b',j,i,  tod%scanid(scan), any(sd%s_zodi_scat(:,i:i,j)/=sd%s_zodi_scat(:,i:i,j)), any(sd%s_zodi_therm(:,i:i,j)/=sd%s_zodi_therm(:,i:i,j))
             call get_s_zodi(tod%id, sd%s_zodi_therm(:, i:i, j), sd%s_zodi_scat(:, i:i, j), d_calib(8 + i, :, j), comp_id=i)
         end do
       end if
      !  Bandpass proposals
       do i = 1, nout-tod%output_n_maps
          d_calib(tod%output_n_maps+i,:,j) = d_calib(1,:,j) + sd%s_bp(:,j) - sd%s_bp_prop(:,j,i+1)
       end do

    end do

  end subroutine compute_calibrated_data



  subroutine simulate_tod(self, scan_id, s_tot, n_corr, handle)
    !
    ! Commander3 native simulation routine. It simulates  correlated
    ! noise component, adds it to the commander-sampled total sky 
    ! signal (multiplied by gain factor for a given frequency) and 
    ! overwrites the original timestreams inside copied files.
    !
    !  Arguments:
    !  ----------
    !  s_tot:    real(sp), array(:,:)
    !            Total sky signal 
    !  scan_id:  integer(i4b)
    !            Local scan ID for the current core
    !  handle:   planck_rng derived type
    !            Healpix definition for random number generation
    !
    !  Returns:
    !  --------
    !
    implicit none
    class(comm_tod), intent(inout) :: self
    ! Parameter file variables
    !type(comm_params),                     intent(in)    :: cpar
    ! Other input/output variables
    real(sp), allocatable, dimension(:,:), intent(in)    :: s_tot   !< total sky signal
    real(sp),              dimension(:,:), intent(out)   :: n_corr  !< Correlated noise (output)
    integer(i4b),                          intent(in)    :: scan_id !< current PID
    type(planck_rng),                      intent(inout) :: handle
    ! Simulation variables
    real(sp), allocatable, dimension(:,:) :: tod_per_detector !< simulated tods per detector
    real(sp)                              :: gain   !< detector's gain value
    real(sp)                              :: sigma0
    real(sp) :: N_c
    real(sp) :: samprate
    real(sp) :: fft_norm
    real(dp) :: chisq
    integer(i4b)                          :: ntod !< total amount of ODs
    integer(i4b)                          :: ndet !< total amount of detectors
    byte,    allocatable, dimension(:)    :: ztod 

    ! HDF5 variables
    character(len=6)   :: samptext, scantext
    character(len=512) :: mystring, mysubstring !< dummy values for string manipulation
    integer(i4b)       :: myindex     !< dummy value for string manipulation
    character(len=512) :: currentHDFFile !< hdf5 file which stores simulation output
    character(len=6)   :: pidLabel
    character(len=512) :: detectorLabel
    type(hdf_file)     :: hdf5_file   !< hdf5 file to work with
    type(hdf_file)     :: tod_file
    integer(i4b)       :: hdf5_error  !< hdf5 error status
    integer(HID_T)     :: hdf5_file_id !< File identifier
    integer(HID_T)     :: dset_id     !< Dataset identifier
    integer(hid_t)     :: dtype  ! hdf5 datatype
    integer(HSIZE_T), dimension(1) :: dims
    ! Other variables
    integer(i4b)                          :: i, j, k !< loop variables
    integer(i4b)       :: mpi_err, errorcode !< MPI error status
    integer(i4b)       :: nomp !< Number of threads available
    integer(i4b)       :: omp_err !< OpenMP error status
    integer(i4b) :: omp_get_max_threads
    integer(i4b) :: n, nfft
    integer*8    :: plan_back
    real(sp) :: nu
    !real(sp), allocatable, dimension(:,:) :: n_corr
    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv
    character(len=10) :: processor_label   !< to have a nice output to screen
    integer(i4b) :: ntoks
    character(len=512), dimension(100) :: toks

    !write(*,*) 'sim', self%scanid(scan_id), self%scans(scan_id)%d%accept

    ! shortcuts
    ntod = self%scans(scan_id)%ntod
    ndet = self%ndet

    ! Simulating 1/f noise
    !write(*,*) "Simulating correlated noise"
    nfft = 2 * ntod
    n = nfft / 2 + 1
    nomp = omp_get_max_threads()
    call sfftw_init_threads(omp_err)
    call sfftw_plan_with_nthreads(nomp)
    ! planning FFTW - in principle we should do both forward and backward FFTW,
    ! but in this case we can omit forward one and go directly with backward to
    ! save some time on a whole operation.
    allocate(dt(nfft), dv(0:n-1))
    call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)

    !$OMP PARALLEL PRIVATE(i, j, k, dt, dv, sigma0, nu)
    allocate(dt(nfft), dv(0:n-1)) !, n_corr(ntod, ndet))
    !$OMP DO SCHEDULE(guided)
    do j = 1, ndet
      ! skipping iteration if scan was not accepted
      if (.not. self%scans(scan_id)%d(j)%accept) cycle
      ! getting gain for each detector (units, V / K)
      ! (gain is assumed to be CONSTANT for EACH SCAN)
      gain   = self%scans(scan_id)%d(j)%gain
      sigma0 = self%scans(scan_id)%d(j)%N_psd%sigma0
      samprate = self%samprate
      ! used when adding fluctuation terms to Fourier coeffs (depends on Fourier convention)
      fft_norm = sqrt(1.d0 * nfft)
      !
      !dv(0) = dv(0) + fft_norm * sigma0 * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0)
      dv(0) = 0. ! fft_norm * sigma0 * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0) ! HKE: This expression is not correct for the monopole
      do k = 1, (n - 1)
        nu    = k * (samprate / 2) / (n - 1)
        N_c   = self%scans(scan_id)%d(j)%N_psd%eval_corr(nu)
        dv(k) = cmplx(rand_gauss(handle), rand_gauss(handle)) * sqrt(N_c) /sqrt(2.0)
      end do
      ! Executing Backward FFT
      call sfftw_execute_dft_c2r(plan_back, dv, dt)
      dt = dt / sqrt(1.d0*nfft)
      n_corr(:,j) = dt(1:ntod)
      !write(*,*) "n_corr ", n_corr(:, j)
    end do
    !$OMP END DO
    deallocate(dt, dv)
    !$OMP END PARALLEL

    call sfftw_destroy_plan(plan_back)

    ! Allocating main simulations' array
    allocate(tod_per_detector(ntod, ndet))       ! Simulated tod
    tod_per_detector = 1d30

    ! Main simulation loop
    do i = 1, ntod
      do j = 1, ndet
        ! skipping iteration if scan was not accepted
        if (.not. self%scans(scan_id)%d(j)%accept) cycle
        gain   = self%scans(scan_id)%d(j)%gain
        sigma0 = self%scans(scan_id)%d(j)%N_psd%sigma0
        tod_per_detector(i,j) = gain * s_tot(i,j) + n_corr(i, j) + sigma0 * rand_gauss(handle)
      end do
    end do
    ! Digitizes the data to the nearest integer; probably should mimic the
    ! actual ADC conversion process
    if (self%compressed_tod) tod_per_detector = real(nint(tod_per_detector), kind=sp)

    !----------------------------------------------------------------------------------
    ! Saving stuff to hdf file
    ! Getting the full path and name of the current hdf file to overwrite
    !----------------------------------------------------------------------------------
    mystring = trim(self%hdfname(scan_id))
    mysubstring = '/'

    myindex = index(trim(mystring), trim(mysubstring), back=.true.) + 1


    call get_tokens(trim(mystring), "/", toks=toks, num=ntoks)
    currentHDFFile = trim(self%sims_output_dir)//'/'//trim(toks(ntoks))
    call int2string(self%scanid(scan_id), pidLabel)
    call int2string(self%myid, processor_label)
    write(*,*) "!  Process:", self%myid, "started writing PID: "//trim(pidLabel)//", into:"
    write(*,*) "!  "//trim(toks(ntoks))

    dims(1) = ntod
    call h5open_f(hdf5_error)
    call  h5fopen_f(currentHDFFile, H5F_ACC_RDWR_F, hdf5_file_id, hdf5_error)
    if (hdf5_error /= 0) call h5eprint_f(hdf5_error)

    ! Remake huffman, symbols for tod_per_detector
    ! decompress the zipped tods to remake the tod
    !do j = 1, 4
    !   if (.not. self%scans(scan_id)%d(j)%accept) cycle
    !   call self%decompress_tod(scan_id, j, tod_per_detector(:,j))
    !end do
    ! call hufmak(tod_per_detector, self%scans(scan_id)%todkey)
    ! Need to overwrite the keys in the simulated data

    if (self%compressed_tod) then
      call open_hdf_file(trim(self%sims_output_dir)//'/tod_'//pidLabel//'.h5', tod_file, 'w')
      do k = 1, self%ndet
        detectorLabel = self%label(k)

        call write_hdf(tod_file, '/'//trim(detectorLabel), tod_per_detector(:,k))
        call write_hdf(tod_file, '/xi_n_'//trim(detectorLabel), self%scans(scan_id)%d(k)%N_psd%xi_n)
        call write_hdf(tod_file, '/gain_'//trim(detectorLabel), self%scans(scan_id)%d(k)%gain)

      end do
      call write_hdf(tod_file, '/x_im', self%x_im)
      call close_hdf_file(tod_file)
    else
      do j = 1, ndet
        detectorLabel = self%label(j)
        call h5dopen_f(hdf5_file_id, trim(pidLabel)//'/'//trim(detectorLabel)//'/'//'tod', dset_id, hdf5_error)
        call h5dwrite_f(dset_id, H5T_IEEE_F32LE, tod_per_detector(:,j), dims, hdf5_error)
      end do
    end if
    call h5dclose_f(dset_id, hdf5_error)
    call h5fclose_f(hdf5_file_id, hdf5_error)
    call h5close_f(hdf5_error)


    deallocate(tod_per_detector)
    write(*,*) "!  Process:", self%myid, "finished writing PID: "//trim(pidLabel)//"."

  end subroutine simulate_tod


end module comm_tod_driver_mod
