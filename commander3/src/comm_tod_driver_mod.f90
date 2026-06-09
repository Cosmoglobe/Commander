module comm_tod_driver_mod
  use comm_tod_mod
  use comm_tod_gain_mod
  use comm_tod_noise_mod
  use comm_tod_pointing_mod
  use comm_tod_bandpass_mod
  use comm_tod_orbdipole_mod
  !use comm_tod_simulations_mod
  !use comm_tod_mapmaking_mod
  use comm_tod_jump_mod
  use comm_tod_adc_mod
  use comm_zodi_mod
  use comm_tod_objctr_mod
  use comm_tod_cray_mod
  use comm_shared_arr_mod
  use comm_huffman_mod
  use comm_tod_dynmask_mod
  use comm_fft_mod
  use omp_lib
  implicit none

contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Scan data routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! HFI nonlin levels (0) = skip
  !                   (1) = demodulate
  !                   (2) = (1) + adc corrections
  !                   (3) = (2) + fill gaps + rolloff deconvolution
  !                   (4) = (3) + 4K lines correction
  !
  ! spurious levels   (0) = skip
  !                   (1) = monopole
  !                   (2) = (1) + jumps
  !                   (3) = (2) + instrument-specific
  !                   (4) = (3) + spikes
  subroutine init_scan_data(tod, scan, oper, bitmask0, sd, det, nonlin_level, spur_level, handle)
    implicit none
    class(comm_tod),      intent(inout)           :: tod
    integer(i4b),         intent(in)              :: scan, oper, bitmask0
    class(comm_scandata), intent(inout)           :: sd    
    integer(i4b),         intent(in),    optional :: det
    integer(i4b),         intent(in),    optional :: nonlin_level
    integer(i4b),         intent(in),    optional :: spur_level 
    type(planck_rng),     intent(inout), optional :: handle

    integer(i4b) :: i, j, h, k, b, d, ntod, hmax, ndet, nhorn, nbp, nonlin_lvl, spur_lvl
    integer(i4b) :: n, nfft
    real(sp),     allocatable, dimension(:)   :: dt
    complex(spc), allocatable, dimension(:)   :: dv

    call timer%start(TOD_ALLOC, tod%band)

    ! Initialize defaults and check dependencies
    sd%scan     = scan
    sd%oper     = oper
    sd%ntod     = tod%scans(scan)%ntod
    sd%ndet     = tod%ndet; if (present(det)) sd%ndet = 1
    sd%nhorn    = tod%nhorn
    !sd%bitmask0 = bitmask0
    sd%bitmask0 = 0
    sd%nbp   = 1; if (btest(oper, SD_SKY_PROP) .or. btest(oper,SD_BP_PROP)) &
         & sd%nbp = size(tod%pixcache%map_sky,4)
    sd%hmax     = 0; if (tod%nhorn > 1) sd%hmax = tod%nhorn
    nonlin_lvl  = 100; if (present(nonlin_level)) nonlin_lvl = nonlin_level
    spur_lvl    = 100; if (present(spur_level))   spur_lvl = spur_level
    
    ! Allocate data structures
    ntod  = sd%ntod
    ndet  = sd%ndet
    hmax  = sd%hmax
    nhorn = sd%nhorn
    nbp   = sd%nbp 
                                allocate(sd%det     (ndet))
    if (btest(oper,SD_IND))     allocate(sd%ind     (ntod, ndet, nhorn))
    if (btest(oper,SD_BASE))    allocate(sd%flag    (ntod, ndet))
    if (btest(oper,SD_MASK))    allocate(sd%mask    (ntod, ndet))
    if (btest(oper,SD_BASE))    allocate(sd%pix     (ntod, ndet, nhorn))
    if (btest(oper,SD_BASE))    allocate(sd%psi     (ntod, ndet, nhorn))
    if (btest(oper,SD_TOD))     allocate(sd%tod     (ntod, ndet))
    if (btest(oper,SD_ORB))     allocate(sd%s_orb   (ntod, ndet, 0:hmax))
    if (btest(oper,SD_SL))      allocate(sd%s_sl    (ntod, ndet, 0:hmax))
    if (btest(oper,SD_ZODI) .and. tod%subtract_zodi)    allocate(sd%s_zodi  (ntod, ndet, 0:hmax))
    if (btest(oper,SD_OBJCTR))  allocate(sd%s_objctr(ntod, ndet, 0:hmax))
    if (btest(oper,SD_SKY))     allocate(sd%s_sky   (ntod, ndet, 0:hmax, nbp))
    if (btest(oper,SD_BP))      allocate(sd%s_bp    (ntod, ndet, 0:hmax, nbp))
    if (btest(oper,SD_TOT))     allocate(sd%s_tot   (ntod, ndet, 0:hmax, nbp))
    if (btest(oper,SD_GAIN))    allocate(sd%s_gain  (ntod, ndet, 0:hmax))
    if (btest(oper,SD_NCORR))   allocate(sd%n_corr  (ntod, ndet))
    if (btest(oper,SD_MONO))    allocate(sd%s_mono  (ntod, ndet))
    if (btest(oper,SD_INST))    allocate(sd%s_inst  (ntod, ndet))
    if (btest(oper,SD_JUMP))    allocate(sd%s_jump  (ntod, ndet))
    if (btest(oper,SD_SPIKE))   allocate(sd%s_spike (ntod, ndet))
    if (btest(oper,SD_DARK))    allocate(sd%dark    (ntod, tod%ndark))
    if (btest(oper,SD_SPUR))    allocate(sd%s_spur  (ntod, ndet))
    call timer%stop(TOD_ALLOC, tod%band)


    if (btest(oper,SD_NCORR) .or. tod%correct_Tbol) then
       sd%enable_fft = .true.
       nfft     = get_closest_fft_pow2(sd%ntod)
       n        = nfft / 2 + 1
       allocate(dt(nfft), dv(0:n-1))
       sd%plan_fwd  = fftwf_plan_dft_r2c_1d(nfft, dt, dv, fftw_estimate + fftw_unaligned)
       sd%plan_back = fftwf_plan_dft_c2r_1d(nfft, dv, dt, fftw_estimate + fftw_unaligned)
       deallocate(dt, dv)
    end if

    ! Initialize detector list
    if (present(det)) then
       sd%det(1) = det
    else
       do i = 1, ndet
          sd%det(i) = i
       end do
    end if
    
    ! Initialize arrays that are not set in this routine
    if (btest(oper,SD_NCORR)) sd%n_corr = 0.
    
    ! Decompress pointing, psi and flags for current scan
    call timer%start(TOD_DECOMP, tod%band)
    if (btest(oper,SD_BASE)) then
       call tod%decompress_pointing(sd, det)
       call tod%decompress_flags(sd, det)
       call tod%apply_fast_flags_inst(sd)
       if (sd%ndet > 1 .and. tod%symm_flags) call tod%symmetrize_flags(sd%flag)
    end if
       
    ! Decompress dark detectors
    if (btest(oper,SD_DARK)) then
      do j = 1, tod%ndark
        call tod%decompress_dark_data(scan, j, sd%dark(:,j))
      end do
    end if
    call timer%stop(TOD_DECOMP, tod%band)
    
    ! Prepare TOD
    if (btest(oper,SD_TOD)) then
       if (tod%ndiode == 1) then
          call timer%start(TOD_DECOMP, tod%band)
          do j = 1, ndet
             d = j; if (present(det)) d = det
             if (.not. tod%scans(scan)%d(d)%accept) cycle
             if (tod%compressed_tod) then
                call tod%decompress_tod(scan, d, sd%tod(:,j))
             else
                sd%tod(:,j) = tod%scans(scan)%d(d)%tod
             end if
          end do
          call timer%stop(TOD_DECOMP, tod%band)
       else
          call tod%diode2tod_inst(sd)
       end if
    end if

    ! Precompute pix2ind
    if (btest(oper,SD_IND)) then
       call timer%start(TOD_PIX2IND, tod%band)
       call compute_scan_pix2ind(tod, sd)
       call timer%stop(TOD_PIX2IND, tod%band)
    end if

    ! Disable broken detector-scans
    do j = 1, ndet
       d = j; if (present(det)) d = det
       if (.not. tod%scans(scan)%d(d)%accept) cycle
       if (tod%scans(scan)%d(d)%N_psd%sigma0 <= 0.d0) &
            & tod%scans(scan)%d(d)%accept = .false.
    end do

    
    ! Construct mask
    if (btest(oper,SD_MASK)) then
       call timer%start(TOD_PROJECT, tod%band)
       call project_mask(tod, bitmask0, sd)
       call timer%stop(TOD_PROJECT, tod%band)
       sd%mask = 1.0

       ! Disable broken detector-scans
!!$       do j = 1, ndet
!!$          d = sd%det(j)
!!$          if (.not. tod%scans(scan)%d(d)%accept) cycle
!!$          if (all(sd%mask(:,j) == 0)) tod%scans(scan)%d(d)%accept = .false.
!!$       end do
    end if
    
    ! Construct sky signal template
    if (btest(oper,SD_SKY)) then
       call timer%start(TOD_PROJECT, tod%band)
       call project_sky(tod, sd)
       call timer%stop(TOD_PROJECT, tod%band)
    end if
    
    ! Construct orbital dipole template
    if (btest(oper,SD_ORB)) then
       call timer%start(TOD_ORBITAL, tod%band)
       call tod%construct_orbital_dipole(sd, det)
       call timer%stop(TOD_ORBITAL, tod%band)
    end if

    ! Construct zodical light template
    if (btest(oper,SD_ZODI) .and. tod%subtract_zodi) then
       call timer%start(TOD_ZODI, tod%band)
       call compute_zodi(zodi_model, tod, sd, det)
       call timer%stop(TOD_ZODI, tod%band)
    end if

    ! Construct object-centric signals
    if (btest(oper,SD_OBJCTR)) then
       call timer%start(TOD_OBJCTR, tod%band)
       call compute_obj_centric_signal(tod, sd, det)
       call timer%stop(TOD_OBJCTR, tod%band)
    end if

    ! Construct sidelobe template
    if (btest(oper,SD_SL)) then
       call timer%start(TOD_SL_INT, tod%band)
       call tod%construct_sl_template(sd, det)
       sd%s_sl = 2.d0 * sd%s_sl ! Scaling by a factor of 2, by comparison with LevelS. Should be understood
       call timer%stop(TOD_SL_INT, tod%band)
    end if

    ! Apply Tbol convolution to all optical components
    if (tod%correct_Tbol) then
!!$       if (tod%myid == 0 .and. scan == 1 .and. allocated(sd%s_sky)) then
!!$          open(58,file='sky_before_Tbol.dat', recl=1024)
!!$          do i = 1, ntod
!!$             write(58,*) i, sd%s_sky(i,1,0,1), sd%s_orb(i,1,0)
!!$          end do
!!$          close(58)
!!$       end if
       
       do j = 1, ndet
          d = j; if (present(det)) d = det
          if (.not. tod%scans(scan)%d(d)%accept) cycle
          do h = 0, hmax
             do b = 1, nbp
                if (allocated(sd%s_sky)) call tod%Tbol(d)%p%convolve("full", "convolve", sd%s_sky(:,j,h,b), sd%plan_fwd, sd%plan_back)
                if (allocated(sd%s_bp))  call tod%Tbol(d)%p%convolve("full", "convolve", sd%s_bp(:,j,h,b), sd%plan_fwd, sd%plan_back)
             end do
             if (allocated(sd%s_orb))    call tod%Tbol(d)%p%convolve("full", "convolve", sd%s_orb(:,j,h), sd%plan_fwd, sd%plan_back)
             if (allocated(sd%s_sl))     call tod%Tbol(d)%p%convolve("full", "convolve", sd%s_sl(:,j,h), sd%plan_fwd, sd%plan_back)
             if (allocated(sd%s_zodi))   call tod%Tbol(d)%p%convolve("full", "convolve", sd%s_zodi(:,j,h), sd%plan_fwd, sd%plan_back)
             if (allocated(sd%s_objctr)) call tod%Tbol(d)%p%convolve("full", "convolve", sd%s_objctr(:,j,h), sd%plan_fwd, sd%plan_back)
             if (allocated(sd%s_gain))   call tod%Tbol(d)%p%convolve("full", "convolve", sd%s_gain(:,j,h), sd%plan_fwd, sd%plan_back)
          end do
       end do

!!$       if (tod%myid == 0 .and. scan == 1 .and. allocated(sd%s_sky)) then
!!$          open(58,file='sky_after_Tbol.dat', recl=1024)
!!$          do i = 1, ntod
!!$             write(58,*) i, sd%s_sky(i,1,0,1), sd%s_orb(i,1,0)
!!$          end do
!!$          close(58)
!!$          write(*,*) 'Done'
!!$          call mpi_finalize(i)
!!$          stop
!!$       end if
    end if

    ! Construct monopole correction template
    if (btest(oper,SD_MONO)) then
       do j = 1, ndet
          d = j; if (present(det)) d = det
          if (.not. tod%scans(scan)%d(d)%accept) cycle
          !sd%s_mono(:,j) = tod%mono(j)
          sd%s_mono(:,j) = 0.d0 ! Disabled for now
       end do
    end if

    ! Construct jump correction template
    if (btest(oper,SD_JUMP)) then
       call timer%start(TOD_INSTCORR, tod%band)
       !call tod%construct_jump_corr(sd, det)
       call timer%stop(TOD_INSTCORR, tod%band)
    end if

    ! Construct spike correction template
    if (btest(oper,SD_SPIKE)) then
       call timer%start(TOD_INSTCORR, tod%band)
       call tod%construct_spike_corr(sd)
       call timer%stop(TOD_INSTCORR, tod%band)
    end if
    
    ! Coadd optical components of total sky signal
    if (btest(oper,SD_TOT)) then
       sd%s_tot = 0.
       do j = 1, ndet
          d = j; if (present(det)) d = det
          if (.not. tod%scans(scan)%d(d)%accept) cycle
          if (btest(oper,SD_SKY)) sd%s_tot(:,d,:,:) = sd%s_tot(:,d,:,:) + sd%s_sky(:,d,:,:)
          !write(*,*) 'a4 ', j, sd%s_tot(1,j,0,1)
          do k = 1, nbp
             !write(*,*) 'sky', j, sd%s_sky(1,j,0,1), k
             if (btest(oper,SD_SL))     sd%s_tot(:,d,:,k) = sd%s_tot(:,d,:,k) + sd%s_sl(:,d,:)
             if (btest(oper,SD_ORB))    sd%s_tot(:,d,:,k) = sd%s_tot(:,d,:,k) + sd%s_orb(:,d,:)
             if (allocated(sd%s_zodi))  sd%s_tot(:,d,:,k) = sd%s_tot(:,d,:,k) + sd%s_zodi(:,d,:)
             if (btest(oper,SD_OBJCTR)) sd%s_tot(:,d,:,k) = sd%s_tot(:,d,:,k) + sd%s_objctr(:,d,:)
          end do
       end do
    end if

    ! Coadd horn signals into detector signals; store result in 0th horn-row
    if (nhorn > 1) call tod%coadd_horns(sd)

    ! Add non-optical components into a net spurious component
    if (btest(oper,SD_SPUR)) then
       sd%s_spur = 0.
       do j = 1, ndet
          d = j; if (present(det)) d = det
          if (.not. tod%scans(scan)%d(d)%accept) cycle
          if (btest(oper,SD_MONO))  sd%s_spur(:,d) = sd%s_spur(:,d) + sd%s_mono(:,d)
          if (btest(oper,SD_JUMP))  sd%s_spur(:,d) = sd%s_spur(:,d) + sd%s_jump(:,d)
          if (btest(oper,SD_INST))  sd%s_spur(:,d) = sd%s_spur(:,d) + sd%s_inst(:,d)
          if (btest(oper,SD_SPIKE)) sd%s_spur(:,d) = sd%s_spur(:,d) + sd%s_spike(:,d)
       end do
    end if
    
    ! Apply non-linearity corrections
    if (btest(oper,SD_TOD) .and. nonlin_lvl > 0) then
       call timer%start(TOD_NONLIN, tod%band)
       call tod%apply_nonlin_corr_inst(sd, nonlin_lvl, handle)
       call timer%stop(TOD_NONLIN, tod%band)
    end if

    ! Generate instrument-specific correction template
    if (btest(oper,SD_INST) .and. spur_lvl>2) then
       call timer%start(TOD_INSTCORR, tod%band)
       call tod%construct_corrtemp_inst(sd, det)
       call timer%stop(TOD_INSTCORR, tod%band)
    end if
    
    ! Subtract spurious corrections from TOD according to spur_level
    call timer%start(TOD_INSTCORR, tod%band)
    do j = 1, ndet
       d = j; if (present(det)) d = det
       if (.not. tod%scans(scan)%d(d)%accept) cycle
       if (spur_lvl > 0 .and. btest(oper,SD_MONO))  sd%tod(:,d) = sd%tod(:,d) - sd%s_mono(:,d)
       if (spur_lvl > 1 .and. btest(oper,SD_JUMP))  sd%tod(:,d) = sd%tod(:,d) - sd%s_jump(:,d)
       if (spur_lvl > 2 .and. btest(oper,SD_INST))  sd%tod(:,d) = sd%tod(:,d) - sd%s_inst(:,d)
       if (spur_lvl > 3 .and. btest(oper,SD_SPIKE)) sd%tod(:,d) = sd%tod(:,d) - sd%s_spike(:,d)
    end do
    call timer%stop(TOD_INSTCORR, tod%band)
    
  end subroutine init_scan_data
  
  subroutine dealloc_scan_data(sd)
    implicit none
    class(comm_scandata), intent(inout) :: sd    

    sd%ntod = -1; sd%ndet = -1; sd%nhorn = -1; sd%nbp = -1

    ! Deallocate data structures
    if (allocated(sd%det))           deallocate(sd%det)
    if (allocated(sd%ind))           deallocate(sd%ind)
    if (allocated(sd%tod))           deallocate(sd%tod)
    if (allocated(sd%n_corr))        deallocate(sd%n_corr)
    if (allocated(sd%s_sl))          deallocate(sd%s_sl)
    if (allocated(sd%s_sky))         deallocate(sd%s_sky)
    if (allocated(sd%s_orb))         deallocate(sd%s_orb)
    if (allocated(sd%s_tot))         deallocate(sd%s_tot)
    if (allocated(sd%s_spur))        deallocate(sd%s_spur)
    if (allocated(sd%mask))          deallocate(sd%mask)
    if (allocated(sd%pix))           deallocate(sd%pix)
    if (allocated(sd%psi))           deallocate(sd%psi)
    if (allocated(sd%flag))          deallocate(sd%flag)
    if (allocated(sd%s_bp))          deallocate(sd%s_bp)
    if (allocated(sd%s_mono))        deallocate(sd%s_mono)
    if (allocated(sd%s_jump))        deallocate(sd%s_jump)
    if (allocated(sd%s_spike))       deallocate(sd%s_spike)
    if (allocated(sd%s_zodi))        deallocate(sd%s_zodi)
    if (allocated(sd%s_inst))        deallocate(sd%s_inst)
    if (allocated(sd%s_gain))        deallocate(sd%s_gain)
    if (allocated(sd%dark))          deallocate(sd%dark)
    if (allocated(sd%s_objctr))      deallocate(sd%s_objctr)
    if (allocated(sd%s_calib))       deallocate(sd%s_calib)

    if (sd%enable_fft) then
       call fftw_destroy_plan(sd%plan_fwd)
       call fftw_destroy_plan(sd%plan_back)
    end if

  end subroutine dealloc_scan_data


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Detector data routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! initializes a detector data structure for a single detector over the entire flight
  ! - Flagged data are omitted
  ! - nside_pix -> 0 = skip pix; >0 = ring; <0 = nest
  subroutine init_det_data(tod, det, oper, bitmask0, nside_pix, init_mjd, dd, nonlin_level, spur_level)
    implicit none
    class(comm_tod),     intent(inout) :: tod
    integer(i4b),        intent(in)    :: det
    integer(i4b),        intent(in)    :: oper
    integer(i4b),        intent(in)    :: bitmask0
    integer(i4b),        intent(in)    :: nside_pix 
    logical(lgt),        intent(in)    :: init_mjd
    class(comm_detdata), intent(inout) :: dd
    integer(i4b),        intent(in),   optional :: nonlin_level, spur_level
    
    integer(i4b) :: i, j, k, p, ntod, q
    type(comm_scandata) :: sd
    real(sp)     :: g
    real(sp),     allocatable, dimension(:) :: data
    integer(i4b), allocatable, dimension(:) :: pix
    real(dp),     allocatable, dimension(:) :: mjd

    q = (tod%info%nside/nside_pix)**2
    
    ! Find maximum number of samples
    ntod = 0
    do i = 1, tod%nscan
       if (tod%scans(i)%d(det)%accept) ntod = ntod + tod%scans(i)%ntod
    end do

    allocate(data(ntod))
    if (nside_pix /= 0) allocate(pix(ntod))
    if (init_mjd)       allocate(mjd(ntod))
    
    ! Fill in buffer arrays
    k = 0 ! Number of unmasked samples
    do i = 1, tod%nscan
       if (.not. tod%scans(i)%d(det)%accept) cycle
       call init_scan_data(tod, i, oper, bitmask0, sd, det=det, nonlin_level=nonlin_level, spur_level=spur_level)

       g = tod%scans(i)%d(det)%gain
       do j = 1, sd%ntod
          if (sd%mask(j,1) == 1.) then
             k    = k+1

             ! Initialize tod
             data(k) = sd%tod(j,1)
             if (allocated(sd%s_tot))  data(k) = data(k) - g * sd%s_tot(j,1,0,1)
             if (allocated(sd%s_spur)) data(k) = data(k) -     sd%s_spur(j,1)
             if (allocated(sd%n_corr)) data(k) = data(k) -     sd%n_corr(j,1)

             ! Initialize pix
             if (nside_pix /= 0) then
                pix(k) = sd%pix(j,1,1)
                if (nside_pix < 0) then
                   ! Convert to nest
                   call ring2nest(tod%info%nside, pix(k), p)
                   pix(k) = p
                end if

                if (nside_pix < tod%info%nside .and. nside_pix < 0) then 
                   ! Udgrade in nested ordering
                   pix(k) = pix(k) / q
                else
                   write(*,*) 'Init_det_data: Ring udgrade not yet supported'
                   stop
                end if
             end if

             ! Initialize mjd
             if (init_mjd) mjd(k) = tod%scans(i)%t0(1) + real(j-1,dp)/tod%samprate
          end if
       end do
       
       call dealloc_scan_data(sd)
    end do

    ! Store final arrays
    dd%ntod = k
    allocate(dd%tod(dd%ntod))
    dd%tod = data(1:k)
    deallocate(data)
    
    if (nside_pix /= 0) then
       allocate(dd%pix(dd%ntod))
       dd%pix = pix(1:k)
       !if (tod%myid == 0 .and. det == 2) write(*,*) 'd', dd%pix(1:5)
       deallocate(pix)
    end if
    
    if (init_mjd) then
       allocate(dd%mjd(dd%ntod))
       dd%mjd = mjd(1:k)
       deallocate(mjd)
    end if
    
  end subroutine init_det_data

  ! deallocates a det_data struture
  subroutine dealloc_det_data(dd)
    implicit none
    class(comm_detdata), intent(inout) :: dd

    dd%ntod = -1
    if (allocated(dd%tod))   deallocate(dd%tod)
    if (allocated(dd%pix))   deallocate(dd%pix)
    if (allocated(dd%mjd))   deallocate(dd%mjd)
  end subroutine dealloc_det_data
  
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! Scandata 2 detdata
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$  subroutine populate_sd_from_dd(sd, dd, scan, det)
!!$    implicit none
!!$    class(comm_scandata),                      intent(inout)          :: sd
!!$    class(comm_detdata),                       intent(inout)          :: dd
!!$    integer(i4b),                              intent(in)             :: scan
!!$    integer(i4b),                              intent(in)             :: det
!!$
!!$    allocate(sd%tod(det, dd%ntod(scan)))
!!$    allocate(sd%s_inst(det, dd%ntod(scan)))
!!$
!!$  end subroutine populate_sd_from_dd


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Sampling drivers etc.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sample_calibration(tod, mode, oper, handle, polang, smooth, mask_threshold)
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
    integer(i4b),                                 intent(in)    :: oper
    type(planck_rng),                             intent(inout) :: handle
!!$    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
!!$    real(sp),            dimension(0:,1:,1:,1:),   intent(in)    :: map_gain
!!$    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2
    real(dp),                                  intent(in),   optional :: polang
    logical(lgt),                              intent(in),   optional :: smooth
    real(dp),                                  intent(in),   optional :: mask_threshold

    integer(i4b) :: i, j, ext(2), ierr, timer_id
    real(sp)     :: threshold
    real(dp)     :: t1, t2
    real(dp), allocatable, dimension(:)   :: A, b
    real(sp), allocatable, dimension(:,:) :: s_invsqrtN, mask_lowres, s_buf
    real(dp), allocatable, dimension(:,:) :: dipole_mod
    type(comm_scandata) :: sd
    logical(lgt) :: smooth_


    smooth_ = .true.
    threshold = 0.9; if (present(mask_threshold)) threshold=mask_threshold
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

    call timer%start(timer_id, tod%band)

    if (trim(mode) == 'abscal' .or. trim(mode) == 'relcal' .or. trim(mode) == 'imbal') then
       allocate(A(tod%ndet), b(tod%ndet))
       A = 0.d0; b = 0.d0
    else if (trim(mode) == 'deltaG') then
       allocate(dipole_mod(tod%last_scan, tod%ndet))
       dipole_mod = 0.d0
    else
       write(*,*) 'Unsupported sampling mode!'
       stop
    end if

    call timer%stop(timer_id, tod%band)


    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) then
          !write(*,*) '  No accepted samples in scan = ', tod%scanid(i)
          cycle
       end if
       call wall_time(t1)

       ![Debug] if (tod%myid == 0) write(*,*) '|    --> Preparing data ' !on, mode = ', trim(mode)
       ! Prepare data
       call init_scan_data(tod, i, oper, TODMASK_GAIN, sd, handle=handle)

       ![Debug] if (tod%myid == 0) write(*,*) '|    --> Setup filtered calibration signal'! m(mode)
       ! Set up filtered calibration signal, conditional contribution and mask
       call timer%start(timer_id, tod%band)
       call tod%downsample_tod(sd%tod(:,1), ext)
       allocate(s_invsqrtN(ext(1):ext(2), tod%ndet))      ! s * invN
       allocate(s_buf(sd%ntod, sd%ndet))
       allocate(mask_lowres(ext(1):ext(2), tod%ndet))

       !write(*,*) "sample 1: ext(1)= ", ext(1), " ext(2)= ", ext(2)


       do j = 1, tod%ndet

          if (.not. tod%scans(i)%d(j)%accept) then 
            !write(*,*) "sample_calibration is cycling"
            cycle
          end if
          
          call tod%downsample_tod(sd%mask(:,j), ext, mask_lowres(:,j), threshold=threshold)
          !if (size(sd%mask(:,j)) > 0) write(*,*) "fsky", sum(sd%mask(:,j))/size(sd%mask(:,j)), sum(mask_lowres(:,j))/size(mask_lowres(:,j))
          if (trim(mode) == 'abscal') then
             if (trim(tod%abscal_comps) == 'orbital') then
               call tod%downsample_tod(sd%s_orb(:,j,0), ext, s_invsqrtN(:,j))
             else if (trim(tod%abscal_comps) == 'full') then
                ! Calibrator = total signal
                s_buf(:,j) = sd%s_tot(:,j,0,1)
                call fill_all_masked(s_buf(:,j), sd%mask(:,j), sd%ntod, .false., real(tod%scans(i)%d(j)%N_psd%sigma0, sp), handle, tod%scans(i)%chunk_num)
                call tod%downsample_tod(s_buf(:,j), ext, s_invsqrtN(:,j))
             else
                s_buf(:,j) = sd%s_gain(:,j,0)
                call fill_all_masked(s_buf(:,j), sd%mask(:,j), sd%ntod, .false., real(tod%scans(i)%d(j)%N_psd%sigma0, sp), handle, tod%scans(i)%chunk_num)
                call tod%downsample_tod(s_buf(:,j), ext, s_invsqrtN(:,j))
             end if
          else if (trim(mode) == 'imbal' .and. tod%nhorn == 2) then
             ! Calibrator = common mode signal
             s_buf(:,j) = tod%scans(i)%d(j)%gain*(sd%s_tot(:,j,1,1) + sd%s_tot(:,j,2,1))
             call fill_all_masked(s_buf(:,j), sd%mask(:,j), sd%ntod, .false., &
               & real(tod%scans(i)%d(j)%N_psd%sigma0, sp), handle, tod%scans(i)%chunk_num)
             call tod%downsample_tod(s_buf(:,j), ext, s_invsqrtN(:,j))
          else
             ! Calibrator = total signal
             s_buf(:,j) = sd%s_tot(:,j,0,1)
             call fill_all_masked(s_buf(:,j), sd%mask(:,j), sd%ntod, .false., real(tod%scans(i)%d(j)%N_psd%sigma0, sp), handle, tod%scans(i)%chunk_num)
             call tod%downsample_tod(s_buf(:,j), ext, s_invsqrtN(:,j))
          end if
       end do

       !write(*,*) "sample 2: ext(1)= ", ext(1), " ext(2)= ", ext(2)

       ! [Debug] if (tod%myid == 0) write(*,*) '|    --> Passed the loop with downsampls tod'!(mode)
       call multiply_inv_N(tod, i, s_invsqrtN, sampfreq=tod%samprate_lowres, pow=0.5d0)

       if (trim(mode) == 'abscal' .or. trim(mode) == 'relcal' .or. trim(mode) == 'imbal') then
          ! Constant gain terms; accumulate contribution from this scan
          do j = 1, tod%ndet
             if (.not. tod%scans(i)%d(j)%accept) cycle
             if (trim(mode) == 'abscal') then
                if (trim(tod%abscal_comps) == 'orbital') then
                  s_buf(:,j) = real(tod%gain0(0),sp) * (sd%s_tot(:,j,0,1) - sd%s_orb(:,j,0)) + &
                       & real(tod%gain0(j) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j,0,1)
                else if (trim(tod%abscal_comps) == 'full') then
                  s_buf(:,j) = real(tod%gain0(j) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j,0,1)
                else
                  s_buf(:,j) = real(tod%gain0(0),sp) * (sd%s_tot(:,j,0,1) - sd%s_gain(:,j,0)) + &
                       & real(tod%gain0(j) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j,0,1)
                end if
             else if (trim(mode) == 'relcal') then
                s_buf(:,j) = real(tod%gain0(0) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j,0,1)
             else if (trim(mode) == 'imbal') then
                s_buf(:,j) = tod%scans(i)%d(j)%gain * (sd%s_tot(:,j,1,1) - sd%s_tot(:,j,2,1))
             end if
          end do
          call accumulate_abscal(tod, i, sd%mask, s_buf, s_invsqrtN, A, b, handle, &
               & out=.true., mask_lowres=mask_lowres, tod_arr=sd%tod)
       else
          ! Time-variable gain terms
            call calculate_gain_mean_std_per_scan(tod, i, s_invsqrtN, sd%mask, sd%s_tot(:,:,0,1), &
              & handle, mask_lowres=mask_lowres, tod_arr=sd%tod)
          do j = 1, tod%ndet
             if (.not. tod%scans(i)%d(j)%accept) cycle
             dipole_mod(tod%scanid(i),j) = masked_variance(sd%s_sky(:,j,0,1), sd%mask(:,j))
          end do
       end if
       call timer%stop(timer_id, tod%band)

       ![Debug] if (tod%myid == 0) write(*,*) '|    --> Passed another loop'!(mode)
       ! Clean up
       call wall_time(t2)
       tod%scans(i)%proctime   = tod%scans(i)%proctime   + t2-t1
       tod%scans(i)%n_proctime = tod%scans(i)%n_proctime + 1
       call dealloc_scan_data(sd)
       deallocate(s_invsqrtN, s_buf, mask_lowres)
    end do


    call timer%start(TOD_WAIT, tod%band)
    call mpi_barrier(tod%info%comm, ierr)
    call timer%stop(TOD_WAIT, tod%band)

    call timer%start(timer_id, tod%band)
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

    call timer%stop(timer_id, tod%band)

  end subroutine sample_calibration


!!$  ! Sample baseline
!!$  subroutine sample_baseline(tod, handle, map_sky, map_gain, procmask, procmask2, skip_nonlin)
!!$    implicit none
!!$    class(comm_tod),                              intent(inout) :: tod
!!$    type(planck_rng),                             intent(inout) :: handle
!!$    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
!!$    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_gain
!!$    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2
!!$    integer(i4b),                                 intent(in), optional :: skip_nonlin
!!$
!!$    integer(i4b) :: i, j, skip_nonlin_
!!$    real(dp)     :: t1, t2
!!$    type(comm_scandata) :: sd
!!$
!!$    if (tod%myid == 0) write(*,*) '   --> Sampling baseline'
!!$    skip_nonlin_ = 10; if (present(skip_nonlin)) skip_nonlin_=skip_nonlin
!!$
!!$    do i = 1, tod%nscan
!!$       if (.not. any(tod%scans(i)%d%accept)) cycle
!!$       call wall_time(t1)
!!$
!!$       ! Prepare data
!!$       if (tod%nhorn == 1) then
!!$          call init_scan_data_singlehorn(sd, tod, i, map_sky, map_gain, procmask, procmask2,handle_=handle,skip_nonlin=skip_nonlin_)
!!$       else
!!$          call init_scan_data_differential(sd, tod, i, map_sky, map_gain, procmask, procmask2)
!!$       end if
!!$
!!$       do j = 1, tod%ndet
!!$          ! Return the data to its raw state
!!$          sd%tod(:,j) = sd%tod(:,j) + tod%scans(i)%d(j)%baseline
!!$
!!$          ! Estimate the baseline and sample it if requested
!!$          tod%scans(i)%d(j)%baseline =sum((sd%tod(:,j) - tod%scans(i)%d(j)%gain*sd%s_tot(:,j)) &
!!$            & *sd%mask(:,j))/sum(sd%mask(:,j))
!!$          if (trim(tod%operation) == 'sample') then
!!$            tod%scans(i)%d(j)%baseline = tod%scans(i)%d(j)%baseline &
!!$             &  + rand_gauss(handle)/sqrt(sum(sd%mask(:,j)*tod%scans(i)%d(j)%N_psd%sigma0**2))
!!$          end if
!!$       end do
!!$
!!$       ! Clean up
!!$       call wall_time(t2)
!!$       tod%scans(i)%proctime   = tod%scans(i)%proctime   + t2-t1
!!$       tod%scans(i)%n_proctime = tod%scans(i)%n_proctime + 1
!!$       call dealloc_scan_data(sd)
!!$    end do
!!$
!!$  end subroutine sample_baseline
!!$
!!$||||||| 9ba1f051
!!$  ! Sample baseline
!!$  subroutine sample_baseline(tod, handle, map_sky, map_gain, procmask, procmask2)
!!$    implicit none
!!$    class(comm_tod),                              intent(inout) :: tod
!!$    type(planck_rng),                             intent(inout) :: handle
!!$    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
!!$    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_gain
!!$    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2
!!$
!!$    integer(i4b) :: i, j
!!$    real(dp)     :: t1, t2
!!$    type(comm_scandata) :: sd
!!$
!!$    if (tod%myid == 0) write(*,*) '   --> Sampling baseline'
!!$
!!$    do i = 1, tod%nscan
!!$       if (.not. any(tod%scans(i)%d%accept)) cycle
!!$       call wall_time(t1)
!!$
!!$       ! Prepare data
!!$       if (tod%nhorn == 1) then
!!$          call init_scan_data_singlehorn(sd, tod, i, map_sky, map_gain, procmask, procmask2)
!!$       else
!!$          call init_scan_data_differential(sd, tod, i, map_sky, map_gain, procmask, procmask2)
!!$       end if
!!$
!!$       do j = 1, tod%ndet
!!$          ! Return the data to its raw state
!!$          sd%tod(:,j) = sd%tod(:,j) + tod%scans(i)%d(j)%baseline
!!$
!!$          ! Estimate the baseline and sample it if requested
!!$          tod%scans(i)%d(j)%baseline =sum((sd%tod(:,j) - tod%scans(i)%d(j)%gain*sd%s_tot(:,j)) &
!!$            & *sd%mask(:,j))/sum(sd%mask(:,j))
!!$          if (trim(tod%operation) == 'sample') then
!!$            tod%scans(i)%d(j)%baseline = tod%scans(i)%d(j)%baseline &
!!$             &  + rand_gauss(handle)/sqrt(sum(sd%mask(:,j)*tod%scans(i)%d(j)%N_psd%sigma0**2))
!!$          end if
!!$       end do
!!$
!!$       ! Clean up
!!$       call wall_time(t2)
!!$       tod%scans(i)%proctime   = tod%scans(i)%proctime   + t2-t1
!!$       tod%scans(i)%n_proctime = tod%scans(i)%n_proctime + 1
!!$       call dealloc_scan_data(sd)
!!$    end do
!!$
!!$  end subroutine sample_baseline

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

    call timer%start(TOD_CHISQ, tod%band)
    ntod = size(flag,1)
    ndet = size(flag,2)
    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (count(iand(flag(:,j),tod%flag0) .ne. 0) > tod%accept_threshold*ntod) then    ! Discard scans with less than a given percentage of good data
          tod%scans(scan)%d(j)%accept = .false.
          write(*, fmt='(a, i8, a, i8, a, i8)') ' | Reject scan = ', &
            & tod%scanid(scan), ': ', count(iand(flag(:,j),tod%flag0) .ne. 0), &
            &  ' flagged data out of', ntod
       else if (abs(tod%scans(scan)%d(j)%chisq) > tod%chisq_threshold .or. &  ! Discard scans with high chisq or NaNs
            & isNaN(tod%scans(scan)%d(j)%chisq)) then
          write(*,fmt='(a,i8,i5,a,f12.1)') ' | Reject scan, det = ', &
               & tod%scanid(scan), j, ', chisq = ', tod%scans(scan)%d(j)%chisq
          tod%scans(scan)%d(j)%accept = .false.
       else if (abs(tod%scans(scan)%d(j)%N_psd%sigma0) > tod%sigma0_threshold .or. tod%scans(scan)%d(j)%N_psd%sigma0 .le.  0d0) then
          write(*,fmt='(a,i8,i5,a,f12.1)') ' | Reject scan, det = ', &
               & tod%scanid(scan), j, ', sigma0 = ', tod%scans(scan)%d(j)%N_psd%sigma0
          tod%scans(scan)%d(j)%accept = .false.
!!$       else if (abs(tod%scans(scan)%d(j)%N_psd%xi_n(2)) > 1.5d0) then
!!$          write(*,fmt='(a,i,i5,a,f12.1)') ' | Reject scan, det = ', &
!!$               & tod%scanid(scan), j, ', fknee = ', tod%scans(scan)%d(j)%N_psd%xi_n(2)
!!$          tod%scans(scan)%d(j)%accept = .false.
       end if
    end do
    if (tod%symm_flags) then
       do j = 1, ndet
           if (.not. tod%scans(scan)%d(j)%accept .and. tod%partner(j) >= 0) tod%scans(scan)%d(tod%partner(j))%accept = .false.
       end do
    end if
    call timer%stop(TOD_CHISQ, tod%band)
  end subroutine remove_bad_data

  subroutine remove_tod_outliers(tod, window, threshold)
    implicit none 
    class(comm_tod),               intent(inout) :: tod
    integer(i4b),                  intent(in)    :: window
    real(sp),        dimension(7), intent(in)    :: threshold

    integer(i4b) :: i, j, k, l, iter, scan, n, ierr
    integer(i4b), parameter :: nstat = 7
    real(dp)     :: mu, sigma
    real(sp),     allocatable, dimension(:,:,:) :: stat
    logical(lgt), allocatable, dimension(:,:)   :: accept
    character(len=6), dimension(nstat) :: label

    call timer%start(TOD_CHISQ, tod%band)
    label = ['chisq ', 'sigma0', 'fknee ', 'alpha ', 'base  ', 'base1 ', 'base2 ']
    
    ! Collect test statistics from all cores 
    allocate(stat(tod%last_scan,tod%ndet,-1:nstat), accept(tod%last_scan,tod%ndet))
    stat = 0.
    do i = 1, tod%nscan
       scan = tod%scanid(i)
       do j = 1, tod%ndet
          stat(scan,j,0)  = 0.; if (tod%scans(i)%d(j)%accept) stat(scan,j,0) = 1.
          stat(scan,j,-1) = tod%mod_phase(j,i)
          if (.not. tod%scans(i)%d(j)%accept) cycle
          if (threshold(1) > 0.) stat(scan,j,1) = tod%scans(i)%d(j)%chisq
          if (threshold(2) > 0.) stat(scan,j,2) = tod%scans(i)%d(j)%N_psd%xi_n(1)
          if (threshold(3) > 0.) stat(scan,j,3) = tod%scans(i)%d(j)%N_psd%xi_n(2)
          if (threshold(4) > 0.) stat(scan,j,4) = tod%scans(i)%d(j)%N_psd%xi_n(3)
          if (threshold(5) > 0.) stat(scan,j,5) = tod%scans(i)%d(j)%baseline(0)
          if (threshold(6) > 0.) then
             if (tod%mod_phase(j,i) == 1.) then
                stat(scan,j,6) = tod%scans(i)%d(j)%baseline1
             else
                stat(scan,j,6) = tod%scans(i)%d(j)%baseline2
             end if
          end if
          if (threshold(7) > 0.) then
             if (tod%mod_phase(j,i) == 1.) then
                stat(scan,j,7) = tod%scans(i)%d(j)%baseline2
             else
                stat(scan,j,7) = tod%scans(i)%d(j)%baseline1
             end if
          end if
       end do
    end do

    ! Find and exclude outliers
    if (tod%myid == 0) then
       call mpi_reduce(mpi_in_place, stat, size(stat), &
            & MPI_REAL, MPI_SUM, 0, tod%comm, ierr)

       ! Search for outliers; compare given scan to surrounding scans. Iterate,
       ! and let the window separation decrease in each step to eliminate extended
       ! regions of outliers
       accept = .true.
       do iter = 1, 5
          do j = 1, tod%ndet
             do l = 1, nstat
                if (threshold(l) <= 0.) cycle
                do i = 1, tod%nscan_tot
                   if (stat(i,j,0) <= 0. .or. .not. accept(i,j)) cycle
                   mu = 0.0; sigma = 0.0; n = 0
                   do k = max(i-iter*window,1), max(i-(iter-1)*window-1,1)
                      if (stat(k,j,0) == 0.) cycle
                      mu    = mu    + stat(k,j,l)
                      sigma = sigma + real(stat(k,j,l),dp)**2
                      n     = n     + 1
                   end do
                   do k = min(i+(iter-1)*window+1,tod%nscan_tot), min(i+iter*window,tod%nscan_tot)
                      if (stat(k,j,0) == 0.) cycle
                      mu    = mu    + stat(k,j,l)
                      sigma = sigma + real(stat(k,j,l),dp)**2
                      n     = n     + 1
                   end do
                   if (n > window/2 .and. n > 1) then
                      mu    = mu / n
                      sigma = sqrt((sigma/n-mu**2)*n/real(n-1,dp))
                      if (sigma > 0.) then
                         accept(i,j) =  (abs(stat(i,j,l)-mu) < threshold(l)*sigma)
                         if (.not. accept(i,j)) write(*,fmt='(a,i8,i4,a,a,a,f16.3,a,f8.3)') '  Rejecting scan = ', i, j, ', ', label(l), ' = ', stat(i,j,l), ', sigma = ', (stat(i,j,l)-mu)/sigma
                         !if (accept(i,j)) write(*,fmt='(a,i8,i4,a,a,a,f16.3,a,f8.3)') '  Accepting scan = ', i, j, ', ', label(l), ' = ', stat(i,j,l), ', sigma = ', (stat(i,j,l)-mu)/sigma
                      end if
                   end if
                end do
             end do
          end do
       end do

      !  open(58,file='stat.dat', recl=1024)
      !  do i = 1, tod%nscan_tot
      !     write(58,*) i, accept(i,1) .and. stat(i,1,0)==1., accept(i,1), stat(i,1,:)
      !  end do
      !  close(58)
      !  open(58,file='stat2.dat', recl=1024)
      !  do i = 1, tod%nscan_tot
      !     if (accept(i,1) .and. stat(i,1,0) == 1.) then
      !        write(58,*) i, accept(i,1), stat(i,1,:)
      !     end if
      !  end do
      !  close(58)

    else
       call mpi_reduce(stat, stat, size(stat), &
            & MPI_REAL, MPI_SUM, 0, tod%comm, ierr)
    end if

    ! Distribute final accept list, and update local data structure
    call mpi_bcast(accept,  size(accept) , MPI_LOGICAL, 0, tod%comm, ierr)
    do i = 1, tod%nscan
       scan = tod%scanid(i)
       do j = 1, tod%ndet
          tod%scans(i)%d(j)%accept = tod%scans(i)%d(j)%accept .and. accept(scan,j)
       end do
    end do

    ! Clean up
    deallocate(stat, accept)
    call timer%stop(TOD_CHISQ, tod%band)
    
  end subroutine remove_tod_outliers
  
  subroutine compute_chisq_abs_bp(tod, scan, sd, chisq)
    implicit none
    class(comm_tod),                       intent(inout) :: tod
    type(comm_scandata),                   intent(in)    :: sd
    real(dp),            dimension(:,:),   intent(inout) :: chisq

    integer(i4b) :: j, k, scan
    real(sp), allocatable, dimension(:,:) :: s_buf

    scan = sd%scan
    
    call timer%start(TOD_BP, tod%band)
    allocate(s_buf(sd%ntod,sd%ndet))
    do j = 1, tod%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       do k = 1, tod%n_bp_prop+1
          call tod%compute_tod_chisq(sd, j, absbp=.true.)
          chisq(j,k) = chisq(j,k) + tod%scans(scan)%d(j)%chisq_prop
       end do
    end do
    deallocate(s_buf)
    call timer%stop(TOD_BP, tod%band)

  end subroutine compute_chisq_abs_bp

  subroutine compute_calibrated_data(tod, scan, sd, d_calib)
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
    class(comm_tod),                       intent(inout)   :: tod
    integer(i4b),                          intent(in)      :: scan
    type(comm_scandata),                   intent(in)      :: sd
    real(sp),            dimension(:,:,:), intent(out)     :: d_calib
    integer(i4b) :: i, j, k, nout
    real(dp)     :: inv_gain
   !  write(*, *) "s_bp:", sd%s_sky(:,1)
    call timer%start(TOD_MAPBIN, tod%band)
    nout = size(d_calib,1)
    d_calib = 0.d0
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (abs(tod%scans(scan)%d(j)%gain) < tiny(1.0)) then
         inv_gain = 0.0
       else
         inv_gain = 1.0 / tod%scans(scan)%d(j)%gain
       end if
       ! sky signal
       d_calib(1,:,j) = sd%tod(:,j) * inv_gain
       if (allocated(sd%n_corr))   d_calib(1,:,j) = d_calib(1,:,j) - sd%n_corr(:,j) * inv_gain
       if (allocated(sd%s_sl))     d_calib(1,:,j) = d_calib(1,:,j) - sd%s_sl(:,j,0)
       if (allocated(sd%s_orb))    d_calib(1,:,j) = d_calib(1,:,j) - sd%s_orb(:,j,0)
       if (allocated(sd%s_zodi))   d_calib(1,:,j) = d_calib(1,:,j) - sd%s_zodi(:,j,0)
       if (allocated(sd%s_objctr)) d_calib(1,:,j) = d_calib(1,:,j) - sd%s_objctr(:,j,0)
       if (allocated(sd%s_bp))     d_calib(1,:,j) = d_calib(1,:,j) - sd%s_bp(:,j,0,1)
       ! residual
       if (tod%output_n_maps > 1) then
          d_calib(2,:,j) = d_calib(1,:,j)
          if (allocated(sd%s_sky)) d_calib(2,:,j) = d_calib(2,:,j) - sd%s_sky(:,j,0,1)
          if (allocated(sd%s_bp))  d_calib(2,:,j) = d_calib(2,:,j) + sd%s_bp(:,j,0,1)
       end if
       ! ncorr
       if (tod%output_n_maps > 2 .and. allocated(sd%n_corr)) d_calib(3,:,j) = sd%n_corr(:,j) * inv_gain
       ! bandpass
       if (tod%output_n_maps > 3 .and. allocated(sd%s_bp))   d_calib(4,:,j) = sd%s_bp(:,j,0,1)
       ! orbital dipole
       if (tod%output_n_maps > 4 .and. allocated(sd%s_orb))  d_calib(5,:,j) = sd%s_orb(:,j,0)
       ! sidelobes
       if (tod%output_n_maps > 5 .and. allocated(sd%s_sl))   d_calib(6,:,j) = sd%s_sl(:,j,0)          
       ! zodi
       if (tod%output_n_maps > 6 .and. allocated(sd%s_zodi)) d_calib(7,:,j) = sd%s_zodi(:,j,0)
       ! instrument specific
       if (tod%output_n_maps > 7 .and. allocated(sd%s_inst)) d_calib(8,:,j) = (sd%s_inst(:,j) - sum(real(sd%s_inst(:,j),dp)/sd%ntod)) * inv_gain
!!$       if (tod%output_n_maps > 8) then
!!$          do i = 1, zodi_model%n_comps
!!$             call get_s_tot_zodi(zodi_model, tod, j, scan, d_calib(8+i,:,j), pix_dynamic=sd%pix(:,j,:), exclude_static='all', comp=i)
!!$         end do
!!$       end if
!!$      !  Bandpass proposals
!!$       !if(tod%n_bp_prop > 1) then
!!$       do i = 1, nout-tod%output_n_maps
!!$          d_calib(tod%output_n_maps+i,:,j) = d_calib(1,:,j) + sd%s_bp(:,j) - sd%s_bp_prop(:,j,i+1)
!!$       end do
!!$      !end if

    end do

    call timer%stop(TOD_MAPBIN, tod%band)

  end subroutine compute_calibrated_data

end module comm_tod_driver_mod
