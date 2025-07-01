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

  ! Class for uncompressed data of a single detector over the full flight
  type :: comm_detdata

    integer(i4b) :: nscan
    integer(i4b), allocatable, dimension(:) :: ntod, scans
    real(sp), allocatable, dimension(:,:) :: tod
  contains
    procedure init_singlehorn => init_det_data_singlehorn
    procedure dealloc         => dealloc_det_data
  end type comm_detdata

contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Scan data routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_scan_data_singlehorn(sd, tod, scan, map_sky, map_gain, procmask, procmask2, procmask_zodi, &
       & init_s_bp, init_s_bp_prop, init_s_sky_prop, skip_nonlin, skip_zodi,  darkdata)
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
    logical(lgt),                              intent(in),   optional :: skip_zodi
    logical(lgt),                              intent(in),   optional :: darkdata

    integer(i4b) :: i, j, k, ndelta
    logical(lgt) :: init_s_bp_, init_s_bp_prop_, init_s_sky_prop_, skip_nonlin_, darkdata_, skip_zodi_

    call timer%start(TOD_ALLOC, tod%band)


    if (tod%nhorn /= 1) then
       write(*,*) 'Error: init_scan_data_singlehorn only applicable for 1-horn experiments'
       stop
    end if
        !if (.true. .or. tod%myid == 78) write(*,*) 'c', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    init_s_bp_ = .false.; if (present(init_s_bp)) init_s_bp_ = init_s_bp
    init_s_sky_prop_ = .false.; if (present(init_s_sky_prop)) init_s_sky_prop_ = init_s_sky_prop
    skip_nonlin_ = .false.; if (present(skip_nonlin)) skip_nonlin_ = skip_nonlin
    skip_zodi_ = .false.; if (present(skip_zodi)) skip_zodi_ = skip_zodi
    darkdata_ = .false.; if (present(darkdata)) darkdata_ = darkdata
 
    init_s_bp_prop_ = .false.
    if (present(init_s_bp_prop)) then
       init_s_bp_prop_  = init_s_bp_prop
       init_s_sky_prop_ = init_s_bp_prop
    end if
    !if (.true. .or. tod%myid == 78) write(*,*) 'c1', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    sd%ntod   = tod%scans(scan)%ntod
    sd%ndet   = tod%ndet
    sd%nhorn  = tod%nhorn
    sd%ndelta = 0; if (init_s_sky_prop_ .or. init_s_bp_prop_) sd%ndelta = size(map_sky,4)

    ! Allocate data structures
    allocate(sd%tod(sd%ntod, sd%ndet))
    allocate(sd%n_corr(sd%ntod, sd%ndet))
    allocate(sd%s_sl(sd%ntod, sd%ndet))
    allocate(sd%s_sky(sd%ntod, sd%ndet))
    allocate(sd%s_bp(sd%ntod, sd%ndet))
    allocate(sd%s_orb(sd%ntod, sd%ndet))
    allocate(sd%s_tot(sd%ntod, sd%ndet))
    allocate(sd%s_gain(sd%ntod, sd%ndet))
    allocate(sd%mask(sd%ntod, sd%ndet))
    allocate(sd%pix(sd%ntod, sd%ndet, sd%nhorn))
    allocate(sd%psi(sd%ntod, sd%ndet, sd%nhorn))
    allocate(sd%flag(sd%ntod, sd%ndet))
    if (init_s_sky_prop_)    allocate(sd%s_sky_prop(sd%ntod, sd%ndet, 2:sd%ndelta))
    if (init_s_bp_prop_)     allocate(sd%s_bp_prop(sd%ntod, sd%ndet, 2:sd%ndelta))
    if (init_s_sky_prop_)    allocate(sd%mask2(sd%ntod, sd%ndet))
    if (tod%sample_mono)     allocate(sd%s_mono(sd%ntod, sd%ndet))
    if (tod%apply_inst_corr) allocate(sd%s_inst(sd%ntod, sd%ndet))

    if (darkdata_) then
        allocate(sd%dark(sd%ntod, tod%ndark))
    end if

    if (tod%subtract_zodi) then
      allocate(sd%s_zodi(sd%ntod, sd%ndet))
      if (tod%sample_zodi) allocate(sd%mask_zodi(sd%ntod, sd%ndet))
    end if
    !call update_status(status, "todinit_alloc")
    call timer%stop(TOD_ALLOC, tod%band)

    !if (.true. .or. tod%myid == 78) write(*,*) 'c2', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    ! Decompress pointing, psi and flags for current scan
    call timer%start(TOD_DECOMP, tod%band)
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       call tod%decompress_pointing_and_flags(scan, j, sd%pix(:,j,:), &
            & sd%psi(:,j,:), sd%flag(:,j))
    end do

!!$    open(58,file='decomp.dat', recl=1024)
!!$    do j = 1, sd%ntod
!!$       write(58,*) j, sd%pix(j,1,1), sd%flag(j,1)
!!$    end do
!!$    close(58)


    if(darkdata_) then
      do j=1, tod%ndark
        call tod%decompress_dark_data(scan, j, sd%dark(:,j))
      end do
    end if
    
    call timer%stop(TOD_DECOMP, tod%band)
    !call update_status(status, "todinit_decomp")
    !if (tod%myid == 78) write(*,*) 'c3', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
    if (tod%symm_flags) call tod%symmetrize_flags(sd%flag)
    !call update_status(status, "todinit_symmflag")
    !if (.true. .or. tod%myid == 78) write(*,*) 'c4', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
    
    ! Prepare TOD
    if (tod%ndiode == 1) then
       call timer%start(TOD_DECOMP, tod%band)
       do j = 1, sd%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          if (tod%compressed_tod) then
             call tod%decompress_tod(scan, j, sd%tod(:,j))
          else
             sd%tod(:,j) = tod%scans(scan)%d(j)%tod
          end if
       end do
       call timer%stop(TOD_DECOMP, tod%band)
    else
       call tod%diode2tod_inst(scan, map_sky, procmask, sd%tod)
    end if
    !call update_status(status, "todinit_tod")
    !if (.true. .or. tod%myid == 78) write(*,*) 'c5', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    ! Construct sky signal template
    call timer%start(TOD_PROJECT, tod%band)
    if (init_s_bp_) then
       call project_sky(tod, map_sky(:,:,:,1), sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, &
            & procmask, scan, sd%s_sky, sd%mask, s_bp=sd%s_bp)
       call project_sky(tod, map_gain(:,:,:,1), sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, &
            & procmask, scan, sd%s_gain, sd%mask, s_bp=sd%s_bp)
    else
       call project_sky(tod, map_sky(:,:,:,1), sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, &
            & procmask, scan, sd%s_sky, sd%mask)
       call project_sky(tod, map_gain(:,:,:,1), sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, &
            & procmask, scan, sd%s_gain, sd%mask)
    end if
    !call update_status(status, "todinit_sky")
    !if (tod%myid == 78) write(*,*) 'c6', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    ! Set up (optional) bandpass sampling quantities (s_sky_prop, mask2 and bp_prop)
    if (init_s_bp_prop_) then
       do j = 2, sd%ndelta
          !if (.true. .or. tod%myid == 78) write(*,*) 'c61', j, tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires, size(map_sky,4)
          call project_sky(tod, map_sky(:,:,:,j), sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, &
               & procmask2, scan, sd%s_sky_prop(:,:,j), sd%mask2, s_bp=sd%s_bp_prop(:,:,j))
       end do
    else if (init_s_sky_prop_) then
       do j = 2, sd%ndelta
          !if (.true. .or. tod%myid == 78) write(*,*) 'c62', j, tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
          call project_sky(tod, map_sky(:,:,:,j), sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, &
               & procmask2, scan, sd%s_sky_prop(:,:,j), sd%mask2)
       end do
    end if
    call timer%stop(TOD_PROJECT, tod%band)

    ! Perform sanity tests
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (all(sd%mask(:,j) == 0)) tod%scans(scan)%d(j)%accept = .false.
       if (tod%scans(scan)%d(j)%N_psd%sigma0 <= 0.d0) tod%scans(scan)%d(j)%accept = .false.
    end do
    !call update_status(status, "todinit_sanity")
    !if (.true. .or. tod%myid == 78) write(*,*) 'c8', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
    
    ! Construct orbital dipole template
    if (tod%correct_orb) then
       call timer%start(TOD_ORBITAL, tod%band)
       call tod%construct_dipole_template(scan, sd%pix(:,:,1), sd%psi(:,:,1), sd%s_orb)
       call timer%stop(TOD_ORBITAL, tod%band)
    else
       sd%s_orb = 0.
    end if
    !if (.true. .or. tod%myid == 78) write(*,*) 'c9', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    ! Construct zodical light template
    if (tod%subtract_zodi) then
       if (skip_zodi_) then
          sd%s_zodi = 0.
       else
          call timer%start(TOD_ZODI, tod%band)
          if (tod%myid == 0) write(*, fmt='(a24, i3, a1)') '    --> Computing zodi: ', nint(real(scan-1, sp)/real(tod%nscan,sp) * 100, i4b), '%'
          do j = 1, sd%ndet
             call get_s_tot_zodi(zodi_model, tod, j, scan, sd%s_zodi(:, j), pix_dynamic=sd%pix(:,j,:))
!!$          if (tod%myid == 0) then
!!$             open(58,file='zodi.dat')
!!$             do k =  1, size(self%s_zodi(:,j))
!!$                write(58,*) k, self%s_zodi(k,j), self%mask(k,j)
!!$             end do
!!$             close(58)
!!$          end if
!!$          call mpi_finalize(k)
!!$          stop
          end do
          call timer%stop(TOD_ZODI, tod%band)
       end if
    end if

    ! Construct sidelobe template
    !if (.true. .or. tod%myid == 78) write(*,*) 'd', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
    if (tod%correct_sl) then
       call timer%start(TOD_SL_INT, tod%band)
       do j = 1, sd%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          !if (.true. .or. tod%myid == 78) write(*,*) 'e', tod%myid, j, tod%slconv(j)%p%psires, tod%slconv(j)%p%psisteps
          call tod%construct_sl_template(tod%slconv(j)%p, &
               & sd%pix(:,j,1), sd%psi(:,j,1), sd%s_sl(:,j), tod%mbang(j))
          sd%s_sl(:,j) = 2.d0 * sd%s_sl(:,j) ! Scaling by a factor of 2, by comparison with LevelS. Should be understood
       end do
       call timer%stop(TOD_SL_INT, tod%band)
    else
       do j = 1, sd%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          sd%s_sl(:,j) = 0.
       end do
    end if
!!$    if (tod%scanid(scan) == 3) then
!!$       open(58,file='sidelobe_BP10.dat')
!!$       do k = 1, size(sd%s_sl,1)
!!$          write(58,*) k, sd%s_sl(k,1)
!!$       end do
!!$       close(58)
!!$    end if


    !call update_status(status, "todinit_sl")

    ! Construct monopole correction template
    if (tod%sample_mono) then
       do j = 1, sd%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          !sd%s_mono(:,j) = tod%mono(j)
          sd%s_mono(:,j) = 0.d0 ! Disabled for now
       end do
    end if

    ! Generate and apply instrument-specific correction template
    if (tod%apply_inst_corr) then
       call tod%construct_corrtemp_inst(scan, sd%pix(:,:,1), sd%psi(:,:,1), sd%s_inst)
!!$       do j = 1, sd%ndet
!!$          if (.not. tod%scans(scan)%d(j)%accept) cycle
!!$          sd%tod(:,j) = sd%tod(:,j) - sd%s_inst(:,j)
!!$       end do
       call timer%stop(TOD_INSTCORR, tod%band)
    end if
    !call update_status(status, "todinit_instcorr")

    ! Construct total sky signal
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (tod%subtract_zodi) then 
           sd%s_tot(:,j) = sd%s_sky(:,j) + sd%s_sl(:,j) + sd%s_orb(:,j) + sd%s_zodi(:,j)
       else
           sd%s_tot(:,j) = sd%s_sky(:,j) + sd%s_sl(:,j) + sd%s_orb(:,j)
       end if
       if (tod%sample_mono) sd%s_tot(:,j) = sd%s_tot(:,j) + sd%s_mono(:,j)
       if (tod%apply_inst_corr) sd%s_tot(:,j) = sd%s_tot(:,j) + sd%s_inst(:,j)
    end do

    ! Apply non-linearity corrections
    if (.not. skip_nonlin_) call tod%apply_nonlin_corr_inst(scan, sd)
    
    !call update_status(status, "todinit_stot")

  end subroutine init_scan_data_singlehorn


  
  subroutine init_scan_data_singlehorn_singledet(sd, tod, det, scan, map_sky, procmask, &
         & skip_nonlin, skip_zodi, darkdata)
    implicit none
    class(comm_scandata),                      intent(inout)          :: sd    
    class(comm_tod),                           intent(inout)          :: tod
    integer(i4b),                              intent(in)             :: det
    integer(i4b),                              intent(in)             :: scan
    real(sp),          dimension(1:,1:),       intent(in)             :: map_sky
    real(sp),          dimension(0:),          intent(in)             :: procmask
    logical(lgt),                              intent(in),   optional :: skip_nonlin
    logical(lgt),                              intent(in),   optional :: skip_zodi
    logical(lgt),                              intent(in),   optional :: darkdata

    integer(i4b) :: i, j, k, ndelta
    logical(lgt) :: init_s_bp_, init_s_bp_prop_, init_s_sky_prop_, skip_nonlin_, darkdata_, skip_zodi_

    call timer%start(TOD_ALLOC, tod%band)


    if (tod%nhorn /= 1) then
       write(*,*) 'Error: init_scan_data_singlehorn_singlescan only applicable for 1-horn experiments'
       stop
    end if
        !if (.true. .or. tod%myid == 78) write(*,*) 'c', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    skip_nonlin_ = .false.; if (present(skip_nonlin)) skip_nonlin_ = skip_nonlin
    skip_zodi_ = .false.; if (present(skip_zodi)) skip_zodi_ = skip_zodi
    darkdata_ = .false.; if (present(darkdata)) darkdata_ = darkdata
 
    sd%ntod   = tod%scans(scan)%ntod
    sd%ndet   = 1
    sd%nhorn  = 1
    sd%ndelta = 0

    ! Allocate data structures
    allocate(sd%tod(sd%ntod, sd%ndet))
    allocate(sd%n_corr(sd%ntod, sd%ndet))
    allocate(sd%s_sl(sd%ntod, sd%ndet))
    allocate(sd%s_sky(sd%ntod, sd%ndet))
    allocate(sd%s_orb(sd%ntod, sd%ndet))
    allocate(sd%s_tot(sd%ntod, sd%ndet))
    allocate(sd%mask(sd%ntod, sd%ndet))
    allocate(sd%pix(sd%ntod, sd%ndet, sd%nhorn))
    allocate(sd%psi(sd%ntod, sd%ndet, sd%nhorn))
    allocate(sd%flag(sd%ntod, sd%ndet))
    if (tod%apply_inst_corr) allocate(sd%s_inst(sd%ntod, sd%ndet))
    if (darkdata_) allocate(sd%dark(sd%ntod, tod%ndark))

    if (tod%subtract_zodi) then
      allocate(sd%s_zodi(sd%ntod, sd%ndet))
      if (tod%sample_zodi) allocate(sd%mask_zodi(sd%ntod, sd%ndet))
    end if
    !call update_status(status, "todinit_alloc")
    call timer%stop(TOD_ALLOC, tod%band)

    !if (.true. .or. tod%myid == 78) write(*,*) 'c2', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    ! Decompress pointing, psi and flags for current scan
    call timer%start(TOD_DECOMP, tod%band)
    call tod%decompress_pointing_and_flags(scan, det, sd%pix(:,1,:), &
            & sd%psi(:,1,:), sd%flag(:,1))

    if(darkdata_) then
      do j=1, tod%ndark
        call tod%decompress_dark_data(scan, det, sd%dark(:,j))
      end do
    end if
    
    call timer%stop(TOD_DECOMP, tod%band)
    
    ! Prepare TOD
    call timer%start(TOD_DECOMP, tod%band)
    if (tod%compressed_tod) then
       call tod%decompress_tod(scan, det, sd%tod(:,1))
    else
       sd%tod(:,1) = tod%scans(scan)%d(det)%tod
    end if
    call timer%stop(TOD_DECOMP, tod%band)

    ! Construct sky signal template
    call timer%start(TOD_PROJECT, tod%band)
    call project_sky(tod, det, map_sky, sd%pix(:,1,1), &
         & sd%psi(:,1,1), sd%flag(:,1), procmask, scan, sd%s_sky(:,1), sd%mask(:,1))
    call timer%stop(TOD_PROJECT, tod%band)
    
    ! Construct orbital dipole template
    if (tod%correct_orb) then
       call timer%start(TOD_ORBITAL, tod%band)
       call tod%construct_dipole_template(scan, sd%pix(:,:,1), sd%psi(:,:,1), sd%s_orb, det=det)
       call timer%stop(TOD_ORBITAL, tod%band)
    else
       sd%s_orb = 0.
    end if

    ! Construct zodical light template
    if (tod%subtract_zodi) then
       if (skip_zodi_) then
          sd%s_zodi = 0.
       else
          call timer%start(TOD_ZODI, tod%band)
          call get_s_tot_zodi(zodi_model, tod, det, scan, sd%s_zodi(:, j), pix_dynamic=sd%pix(:,1,:))
          call timer%stop(TOD_ZODI, tod%band)
       end if
    end if

    ! Construct sidelobe template
    if (tod%correct_sl) then
       call timer%start(TOD_SL_INT, tod%band)
       call tod%construct_sl_template(tod%slconv(det)%p, &
               & sd%pix(:,1,1), sd%psi(:,1,1), sd%s_sl(:,1), tod%mbang(det))
       sd%s_sl(:,1) = 2.d0 * sd%s_sl(:,1) ! Scaling by a factor of 2, by comparison with LevelS. Should be understood
       call timer%stop(TOD_SL_INT, tod%band)
    else
       sd%s_sl = 0.
    end if

    ! Generate and apply instrument-specific correction template
    if (tod%apply_inst_corr) then
       call tod%construct_corrtemp_inst(scan, sd%pix(:,:,1), sd%psi(:,:,1), sd%s_inst)
       call timer%stop(TOD_INSTCORR, tod%band)
    end if

    ! Construct total sky signal
    if (tod%subtract_zodi) then 
       sd%s_tot(:,1) = sd%s_sky(:,1) + sd%s_sl(:,1) + sd%s_orb(:,1) + sd%s_zodi(:,1)
    else
       sd%s_tot(:,1) = sd%s_sky(:,1) + sd%s_sl(:,1) + sd%s_orb(:,1)
    end if
    if (tod%apply_inst_corr) sd%s_tot(:,1) = sd%s_tot(:,1) + sd%s_inst(:,1)

    ! Apply non-linearity corrections
    if (.not. skip_nonlin_) call tod%apply_nonlin_corr_inst(scan, sd, det)
    
    !call update_status(status, "todinit_stot")

  end subroutine init_scan_data_singlehorn_singledet


  subroutine init_scan_data_differential(sd, tod, scan, map_sky, map_gain, procmask, procmask2, &
       & init_s_bp, init_s_bp_prop, init_s_sky_prop, polang_)
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
    real(dp),                                  intent(in),   optional :: polang_

    integer(i4b) :: j, k, ndelta
    logical(lgt) :: init_s_bp_, init_s_bp_prop_, init_s_sky_prop_
    real(sp),     allocatable, dimension(:,:)     :: s_bufA, s_bufB, s_buf2A, s_buf2B      ! Buffer
    real(dp) :: polang
    !
    ! Note that procmask should be larger, cover the residuals in the galactic
    ! plane, so that it is not used for calibration and correlated noise. 
    ! procmask2 should be smaller and be used only for mapmaking.
    !

    if (present(polang_)) then
      polang = polang_
    else
      polang = 0d0
    endif

    call timer%start(TOD_ALLOC, tod%band)
    if (tod%nhorn /= 2) then
       write(*,*) 'Error: init_scan_data_differential only applicable for 2-horn experiments'
       stop
    end if

    init_s_bp_ = .false.; if (present(init_s_bp)) init_s_bp_ = init_s_bp
    init_s_sky_prop_ = .false.; if (present(init_s_sky_prop)) init_s_sky_prop_ = init_s_sky_prop
    init_s_bp_prop_ = .false.
    if (present(init_s_bp_prop)) then
       init_s_bp_prop_  = init_s_bp_prop
       init_s_sky_prop_ = init_s_bp_prop
    end if

    sd%ntod   = tod%scans(scan)%ntod
    sd%ndet   = tod%ndet
    sd%nhorn  = tod%nhorn
    sd%ndelta = 0; if (init_s_sky_prop_ .or. init_s_bp_prop_) sd%ndelta = size(map_sky,4)

    ! Allocate data structures
    allocate(sd%tod(sd%ntod, sd%ndet))
    allocate(sd%n_corr(sd%ntod, sd%ndet))
    allocate(sd%s_sl(sd%ntod, sd%ndet))
    allocate(sd%s_sky(sd%ntod, sd%ndet))
    allocate(sd%s_bp(sd%ntod, sd%ndet))
    allocate(sd%s_orb(sd%ntod, sd%ndet))
    allocate(sd%s_orbA(sd%ntod, sd%ndet))
    allocate(sd%s_orbB(sd%ntod, sd%ndet))
    allocate(sd%s_tot(sd%ntod, sd%ndet))
    allocate(sd%s_totA(sd%ntod, sd%ndet))
    allocate(sd%s_totB(sd%ntod, sd%ndet))
    allocate(sd%s_gain(sd%ntod, sd%ndet))
    allocate(sd%s_gainA(sd%ntod, sd%ndet))
    allocate(sd%s_gainB(sd%ntod, sd%ndet))
    allocate(sd%mask(sd%ntod, sd%ndet))
    allocate(sd%pix(sd%ntod, 1, sd%nhorn))
    allocate(sd%psi(sd%ntod, 1, sd%nhorn))
    allocate(sd%flag(sd%ntod, 1))
    allocate(sd%s_inst(sd%ntod, sd%ndet))
    if (init_s_sky_prop_)   allocate(sd%s_sky_prop(sd%ntod, sd%ndet, 2:sd%ndelta))
    if (init_s_bp_prop_)    allocate(sd%s_bp_prop(sd%ntod, sd%ndet, 2:sd%ndelta))
    if (init_s_sky_prop_)   allocate(sd%mask2(sd%ntod, sd%ndet))
    if (tod%sample_mono)    allocate(sd%s_mono(sd%ntod, sd%ndet))
    if (tod%subtract_zodi)  allocate(sd%s_zodi(sd%ntod, sd%ndet))
    sd%s_tot   = 0.
    sd%s_totA  = 0.
    sd%s_totB  = 0.
    sd%s_orb   = 0.
    sd%s_orbA  = 0.
    sd%s_orbB  = 0.
    sd%s_gain  = 0.
    sd%s_gainA = 0.
    sd%s_gainB = 0.

    allocate(s_bufA(sd%ntod, sd%ndet))
    allocate(s_bufB(sd%ntod, sd%ndet))
    allocate(s_buf2A(sd%ntod, sd%ndet))
    allocate(s_buf2B(sd%ntod, sd%ndet))
    ! Decompress pointing, psi and flags for current scan
    ! Only called for one detector, det=1, since the pointing and polarization
    ! angles are the same for all detectors
    call tod%decompress_pointing_and_flags(scan, 1, sd%pix(:,1,:), &
            & sd%psi(:,1,:), sd%flag(:,1))
    
    ! Prepare TOD
    if (tod%ndiode == 1 .or. trim(tod%level) == 'L2') then
       do j = 1, sd%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          if (tod%compressed_tod) then
             call tod%decompress_tod(scan, j, sd%tod(:,j))
          else
             sd%tod(:,j) = tod%scans(scan)%d(j)%tod
          end if
       end do
       call timer%stop(TOD_DECOMP, tod%band)
    else
       call tod%diode2tod_inst(scan, map_sky, procmask, sd%tod)
    end if

    ! Construct sky signal template
    call timer%start(TOD_PROJECT, tod%band)
    if (init_s_bp_) then
       call project_sky_differential_multi(tod, map_sky(:,:,:,1), sd%pix(:,1,:), sd%psi(:,1,:), sd%flag(:,1), &
            & procmask, scan, s_bufA, s_bufB, sd%mask, s_bpA=s_buf2A, s_bpB=s_buf2B)
    else
       call project_sky_differential_multi(tod, map_sky(:,:,:,1), sd%pix(:,1,:), sd%psi(:,1,:), sd%flag(:,1), &
            & procmask, scan, s_bufA, s_bufB, sd%mask)
    end if
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       sd%s_sky(:,j)  = (1.+tod%x_im(j))*s_bufA(:,j)  - (1.-tod%x_im(j))*s_bufB(:,j)
       sd%s_totA(:,j) = sd%s_totA(:,j) + s_bufA(:,j)
       sd%s_totB(:,j) = sd%s_totB(:,j) + s_bufB(:,j)
       sd%s_tot(:,j)  = sd%s_tot(:,j)  + sd%s_sky(:,j)
       if (init_s_bp_) sd%s_bp(:,j)  = (1.+tod%x_im(j))*s_buf2A(:,j) - (1.-tod%x_im(j))*s_buf2B(:,j)
    end do


    ! Set up (optional) bandpass sampling quantities (s_sky_prop, mask2 and bp_prop)
    if (init_s_bp_prop_) then
       do k = 2, sd%ndelta
          call project_sky_differential_multi(tod, map_sky(:,:,:,k), sd%pix(:,1,:), sd%psi(:,1,:), sd%flag(:,1), &
               & procmask2, scan, s_bufA, s_bufB, sd%mask2, s_bpA=s_buf2A, s_bpB=s_buf2B)
          do j = 1, sd%ndet
             if (.not. tod%scans(scan)%d(j)%accept) cycle
             sd%s_sky_prop(:,j,k) = (1.+tod%x_im(j))*s_bufA(:,j)  - (1.-tod%x_im(j))*s_bufB(:,j)
             if (init_s_bp_) sd%s_bp_prop(:,j,k)  = (1.+tod%x_im(j))*s_buf2A(:,j) - (1.-tod%x_im(j))*s_buf2B(:,j)
          end do
       end do
    else if (init_s_sky_prop_) then
       do k = 2, sd%ndelta
          call project_sky_differential_multi(tod, map_sky(:,:,:,k), sd%pix(:,1,:), sd%psi(:,1,:), sd%flag(:,1), &
               & procmask2, scan, s_bufA, s_bufB, sd%mask2)
          do j = 1, sd%ndet
             if (.not. tod%scans(scan)%d(j)%accept) cycle
             sd%s_sky_prop(:,j,k) = (1.+tod%x_im(j))*s_bufA(:,j) - (1.-tod%x_im(j))*s_bufB(:,j)
          end do
       end do
    end if
    call timer%stop(TOD_PROJECT, tod%band)

    ! Perform sanity tests
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (all(sd%mask(:,j) == 0)) tod%scans(scan)%d(j)%accept = .false.
       if (tod%scans(scan)%d(j)%N_psd%sigma0 <= 0.d0) tod%scans(scan)%d(j)%accept = .false.
    end do
    
    ! Construct orbital dipole template
    call tod%construct_dipole_template_diff(scan, sd%pix(:,:,1), sd%psi(:,:,1), .true., 1, s_bufA, 1d3)
    call tod%construct_dipole_template_diff(scan, sd%pix(:,:,2), sd%psi(:,:,2), .true., 2, s_bufB, 1d3)
    sd%s_orbA = s_bufA
    sd%s_orbB = s_bufB
    sd%s_totA = sd%s_totA + sd%s_orbA
    sd%s_totB = sd%s_totB + sd%s_orbB
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       sd%s_orb(:,j)  = (1.+tod%x_im(j))*sd%s_orbA(:,j)  - (1.-tod%x_im(j))*sd%s_orbB(:,j)
       sd%s_tot(:,j)  = sd%s_tot(:,j)  + sd%s_orb(:,j)
    end do

    ! Construct zodical light template
    if (tod%subtract_zodi) then
       do j = 1, sd%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          sd%s_tot(:,j)  = sd%s_tot(:,j) + sd%s_zodi(:,j)
          sd%s_totA(:,j) = sd%s_totA(:,j) + s_bufA(:,j)
          sd%s_totB(:,j) = sd%s_totB(:,j) + s_bufB(:,j)
       end do
    end if

    ! Construct sidelobe template
    if (tod%correct_sl) then
       do j = 1, sd%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          call tod%construct_sl_template(tod%slconvA(1)%p, sd%pix(:,1,1), sd%psi(:,1,1), s_bufA(:,j),  polang)
          call tod%construct_sl_template(tod%slconvB(3)%p, sd%pix(:,1,2), sd%psi(:,1,2), s_bufB(:,j), -polang)
          sd%s_sl(:,j)  = (1d0+tod%x_im(j))*s_bufA(:,j) - (1d0-tod%x_im(j))*s_bufB(:,j)
          sd%s_tot(:,j) = sd%s_tot(:,j) + sd%s_sl(:,j)
          sd%s_totA(:,j) = sd%s_totA(:,j) + s_bufA(:,j)
          sd%s_totB(:,j) = sd%s_totB(:,j) + s_bufB(:,j)
       end do
       call timer%stop(TOD_SL_INT, tod%band)
    else
       sd%s_sl = 0.
    end if

    ! Construct monopole correction template
    if (tod%sample_mono) then
       do j = 1, sd%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          !sd%s_mono(:,j) = tod%mono(j)
          sd%s_mono(:,j) = 0.d0 ! Disabled for now
          sd%s_tot(:,j)  = sd%s_tot(:,j) + sd%s_mono(:,j)
          sd%s_totA(:,j) = sd%s_totA(:,j) + sd%s_mono(:,j)
          sd%s_totB(:,j) = sd%s_totB(:,j) + sd%s_mono(:,j)
       end do
    else
       sd%s_mono = 0.d0 
    end if

    ! Generate and apply instrument-specific correction template
    if (tod%apply_inst_corr) then
       call timer%start(TOD_INSTCORR, tod%band)
       call tod%construct_corrtemp_inst(scan, sd%pix(:,:,1), sd%psi(:,:,1), sd%s_inst)
       do j = 1, sd%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          sd%tod(:,j) = sd%tod(:,j) - sd%s_inst(:,j)
       end do
       call timer%stop(TOD_INSTCORR, tod%band)
    end if

    ! Construct zodical light template
    if (tod%subtract_zodi) then
       call timer%start(TOD_ZODI, tod%band)
       do j = 1, sd%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
         !  call get_zodi_emission(tod, self%pix(:,1:1,1), scan, s_bufA)
         !  call get_zodi_emission(tod, self%pix(:,1:1,2), scan, s_bufB)
         !  self%s_zodi(:,j) = (1.+tod%x_im(j))*s_bufA(:,j) - (1.-tod%x_im(j))*s_bufB(:,j)
          sd%s_tot(:,j)  = sd%s_tot(:,j)! + sd%s_zodi(:,j)
          sd%s_totA(:,j) = sd%s_totA(:,j) + s_bufA(:,j)
          sd%s_totB(:,j) = sd%s_totB(:,j) + s_bufB(:,j)
       end do
       call timer%stop(TOD_ZODI, tod%band)
    end if

    ! Clean-up
    call timer%start(TOD_ALLOC, tod%band)
    deallocate(s_bufA, s_bufB, s_buf2A, s_buf2B)
    call timer%stop(TOD_ALLOC, tod%band)

  end subroutine init_scan_data_differential


  subroutine dealloc_scan_data(sd)
    implicit none
    class(comm_scandata), intent(inout) :: sd    

    sd%ntod = -1; sd%ndet = -1; sd%nhorn = -1

    ! Deallocate data structures
    if (allocated(sd%tod))           deallocate(sd%tod)
    if (allocated(sd%n_corr))        deallocate(sd%n_corr)
    if (allocated(sd%s_sl))          deallocate(sd%s_sl)
    if (allocated(sd%s_sky))         deallocate(sd%s_sky)
    if (allocated(sd%s_orb))         deallocate(sd%s_orb)
    if (allocated(sd%s_tot))         deallocate(sd%s_tot)
    if (allocated(sd%mask))          deallocate(sd%mask)
    if (allocated(sd%pix))           deallocate(sd%pix)
    if (allocated(sd%psi))           deallocate(sd%psi)
    if (allocated(sd%flag))          deallocate(sd%flag)
    if (allocated(sd%s_sky_prop))    deallocate(sd%s_sky_prop)
    if (allocated(sd%s_bp_prop))     deallocate(sd%s_bp_prop)
    if (allocated(sd%s_bp))          deallocate(sd%s_bp)
    if (allocated(sd%s_mono))        deallocate(sd%s_mono)
    if (allocated(sd%mask2))         deallocate(sd%mask2)
    if (allocated(sd%s_zodi))        deallocate(sd%s_zodi)
    if (allocated(sd%mask_zodi))     deallocate(sd%mask_zodi)
    if (allocated(sd%s_totA))        deallocate(sd%s_totA)
    if (allocated(sd%s_totB))        deallocate(sd%s_totB)
    if (allocated(sd%s_inst))        deallocate(sd%s_inst)
    if (allocated(sd%s_orbA))        deallocate(sd%s_orbA)
    if (allocated(sd%s_orbB))        deallocate(sd%s_orbB)
    if (allocated(sd%s_gain))        deallocate(sd%s_gain)
    if (allocated(sd%s_gainA))       deallocate(sd%s_gainA)
    if (allocated(sd%s_gainB))       deallocate(sd%s_gainB)
    if (allocated(sd%dark))          deallocate(sd%dark)
  end subroutine dealloc_scan_data

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! detdata routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! deallocates a det_data struture
  subroutine dealloc_det_data(dd)
    implicit none
    class(comm_detdata), intent(inout) :: dd

    dd%nscan = -1
    if (allocated(dd%scans)) deallocate(dd%scans)
    if (allocated(dd%ntod))  deallocate(dd%ntod) 
    if (allocated(dd%tod))   deallocate(dd%tod)
  end subroutine dealloc_det_data

  ! initializes a det_data structure for a single detector over the entire flight
  ! for a singlehorn experiment like planck
  subroutine init_det_data_singlehorn(dd, tod, det)
    implicit none
    class(comm_detdata),                       intent(inout)          :: dd
    class(comm_tod),                           intent(inout)          :: tod
    integer(i4b),                              intent(in)             :: det

    integer(i4b) :: i
    
    dd%nscan = tod%nscan
    allocate(dd%ntod(dd%nscan))

    do i = 1, tod%nscan
      ! decompress tod into dd%tod array
      dd%ntod(i) = tod%scans(i)%ntod
    end do

  end subroutine init_det_data_singlehorn


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Scandata 2 detdata
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine populate_sd_from_dd(sd, dd, scan, det)
    implicit none
    class(comm_scandata),                      intent(inout)          :: sd
    class(comm_detdata),                       intent(inout)          :: dd
    integer(i4b),                              intent(in)             :: scan
    integer(i4b),                              intent(in)             :: det

    allocate(sd%tod(det, dd%ntod(scan)))
    allocate(sd%s_inst(det, dd%ntod(scan)))

  end subroutine populate_sd_from_dd


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Sampling drivers etc.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sample_calibration(tod, mode, handle, map_sky, map_gain, procmask, procmask2, &
      & polang, smooth, mask_threshold)
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
       if (tod%nhorn == 1) then
          call init_scan_data_singlehorn(sd, tod, i, map_sky, map_gain, procmask, procmask2)
       else
         if (present(polang)) then
          call init_scan_data_differential(sd, tod, i, map_sky, map_gain, procmask, procmask2, polang_=polang)
        else
          call init_scan_data_differential(sd, tod, i, map_sky, map_gain, procmask, procmask2)
        end if
       end if

       ![Debug] if (tod%myid == 0) write(*,*) '|    --> Setup filtered calibration signal'! m(mode)
       ! Set up filtered calibration signal, conditional contribution and mask
       call timer%start(timer_id, tod%band)
       call tod%downsample_tod(sd%s_orb(:,1), ext)
       allocate(s_invsqrtN(ext(1):ext(2), tod%ndet))      ! s * invN
       allocate(s_buf(sd%ntod, sd%ndet))
       allocate(mask_lowres(ext(1):ext(2), tod%ndet))
       do j = 1, tod%ndet
          if (.not. tod%scans(i)%d(j)%accept) cycle
          call tod%downsample_tod(sd%mask(:,j), ext, mask_lowres(:,j), threshold=threshold)
          !if (size(sd%mask(:,j)) > 0) write(*,*) "fsky", sum(sd%mask(:,j))/size(sd%mask(:,j)), sum(mask_lowres(:,j))/size(mask_lowres(:,j))
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

       ! [Debug] if (tod%myid == 0) write(*,*) '|    --> Passed the loop with downsampls tod'!(mode)
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
             dipole_mod(tod%scanid(i),j) = masked_variance(sd%s_sky(:,j), sd%mask(:,j))
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
          call init_scan_data_singlehorn(sd, tod, i, map_sky, map_gain, procmask, procmask2)
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
       if (count(iand(flag(:,j),tod%flag0) .ne. 0) > tod%accept_threshold*ntod) then    ! Discard scans with less than a given percentage of good data
          tod%scans(scan)%d(j)%accept = .false.
          write(*, fmt='(a, i, a, i8, a, i8)') ' | Reject scan = ', &
            & tod%scanid(scan), ': ', count(iand(flag(:,j),tod%flag0) .ne. 0), &
            &  ' flagged data out of', ntod
!!$       else if (abs(tod%scans(scan)%d(j)%chisq) > tod%chisq_threshold .or. &  ! Discard scans with high chisq or NaNs
!!$            & isNaN(tod%scans(scan)%d(j)%chisq)) then
!!$          write(*,fmt='(a,i,i5,a,f12.1)') ' | Reject scan, det = ', &
!!$               & tod%scanid(scan), j, ', chisq = ', tod%scans(scan)%d(j)%chisq
!!$          tod%scans(scan)%d(j)%accept = .false.
!!$       else if (abs(tod%scans(scan)%d(j)%N_psd%sigma0) > tod%sigma0_threshold) then
!!$          write(*,fmt='(a,i,i5,a,f12.1)') ' | Reject scan, det = ', &
!!$               & tod%scanid(scan), j, ', sigma0 = ', tod%scans(scan)%d(j)%N_psd%sigma0
!!$          tod%scans(scan)%d(j)%accept = .false.
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

    label = ['chisq ', 'sigma0', 'fknee ', 'alpha ', 'base  ', 'base1 ', 'base2 ']
    
    ! Collect test statistics from all cores 
    allocate(stat(tod%nscan_tot,tod%ndet,-1:nstat), accept(tod%nscan_tot,tod%ndet))
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
                   mu = 0.0; n = 0
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
                   if (n > window/2) then
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

       open(58,file='stat.dat', recl=1024)
       do i = 1, tod%nscan_tot
          write(58,*) i, accept(i,1) .and. stat(i,1,0)==1., accept(i,1), stat(i,1,:)
       end do
       close(58)
       open(58,file='stat2.dat', recl=1024)
       do i = 1, tod%nscan_tot
          if (accept(i,1) .and. stat(i,1,0) == 1.) then
             write(58,*) i, accept(i,1), stat(i,1,:)
          end if
       end do
       close(58)

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
    
  end subroutine remove_tod_outliers
  
  subroutine compute_chisq_abs_bp(tod, scan, sd, chisq)
    implicit none
    class(comm_tod),                       intent(inout) :: tod
    integer(i4b),                          intent(in)    :: scan
    type(comm_scandata),                   intent(in)    :: sd
    real(dp),            dimension(:,:),   intent(inout) :: chisq

    integer(i4b) :: j, k
    real(sp), allocatable, dimension(:,:) :: s_buf

    call timer%start(TOD_BP, tod%band)
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
    call timer%stop(TOD_BP, tod%band)

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
    class(comm_tod),                       intent(inout)   :: tod
    integer(i4b),                          intent(in)      :: scan
    type(comm_scandata),                   intent(in)      :: sd
    real(sp),            dimension(:,:,:), intent(out)     :: d_calib
    real(sp), dimension(:,:), intent(in), optional         :: jump_template
    integer(i4b) :: i, j, k, nout
    real(dp)     :: inv_gain
   !  write(*, *) "s_bp:", sd%s_sky(:,1)
    call timer%start(TOD_MAPBIN, tod%band)
    nout = size(d_calib,1)
    d_calib = 0d0
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       inv_gain = 1.0 / tod%scans(scan)%d(j)%gain
       if (tod%compressed_tod) then
!!$          write(*,*) 'calib g', tod%scanid(scan), inv_gain
!!$          write(*,*) 'calib d', tod%scanid(scan), maxval(abs(sd%tod(:,j)))
!!$          write(*,*) 'calib n', tod%scanid(scan), maxval(abs(sd%n_corr(:,j)))
!!$          write(*,*) 'calib t', tod%scanid(scan), maxval(abs(sd%s_tot(:,j)))
!!$          write(*,*) 'calib s', tod%scanid(scan), maxval(abs(sd%s_sky(:,j)))
!!$          write(*,*) 'calib b', tod%scanid(scan), maxval(abs(sd%s_bp(:,j)))
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
       if (tod%output_n_maps > 8) then
          do i = 1, zodi_model%n_comps
             call get_s_tot_zodi(zodi_model, tod, j, scan, d_calib(8+i,:,j), pix_dynamic=sd%pix(:,j,:), exclude_static='all', comp=i)
         end do
       end if
      !  Bandpass proposals
      if(tod%n_bp_prop > 1) then
       do i = 1, nout-tod%output_n_maps
          d_calib(tod%output_n_maps+i,:,j) = d_calib(1,:,j) + sd%s_bp(:,j) - sd%s_bp_prop(:,j,i+1)
       end do
      end if

    end do

    call timer%stop(TOD_MAPBIN, tod%band)

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
      call timer%start(TOT_FFT)
      call sfftw_execute_dft_c2r(plan_back, dv, dt)
      call timer%stop(TOT_FFT)
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
