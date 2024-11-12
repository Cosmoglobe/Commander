module comm_tod_driver_mod
  use comm_tod_mod
  use comm_tod_gain_mod
  use comm_tod_noise_mod
  use comm_tod_pointing_mod
  use comm_tod_bandpass_mod
  use comm_tod_orbdipole_mod
  use comm_tod_simulations_mod
  use comm_tod_mapmaking_mod
  use comm_zodi_mod
  use comm_tod_cray_mod
  use comm_shared_arr_mod
  use omp_lib
  implicit none

contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Scan data routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_scan_data_singlehorn(sd, tod, scan, map_sky, procmask, procmask2, &
       & init_s_bp, init_s_bp_prop, init_s_sky_prop, skip_nonlin)
    implicit none
    class(comm_scandata),                      intent(inout)          :: sd    
    class(comm_tod),                           intent(inout)          :: tod
    integer(i4b),                              intent(in)             :: scan
    real(sp),          dimension(1:,1:,0:,1:), intent(in)             :: map_sky
    real(sp),          dimension(0:),          intent(in)             :: procmask
    real(sp),          dimension(0:),          intent(in)             :: procmask2
    logical(lgt),                              intent(in),   optional :: init_s_bp
    logical(lgt),                              intent(in),   optional :: init_s_bp_prop
    logical(lgt),                              intent(in),   optional :: init_s_sky_prop
    logical(lgt),                              intent(in),   optional :: skip_nonlin

    integer(i4b) :: j, k, ndelta
    logical(lgt) :: init_s_bp_, init_s_bp_prop_, init_s_sky_prop_, skip_nonlin_

    if (tod%nhorn /= 1) then
       write(*,*) 'Error: init_scan_data_singlehorn only applicable for 1-horn experiments'
       stop
    end if
    
    !if (.true. .or. tod%myid == 78) write(*,*) 'c', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    init_s_bp_ = .false.; if (present(init_s_bp)) init_s_bp_ = init_s_bp
    init_s_sky_prop_ = .false.; if (present(init_s_sky_prop)) init_s_sky_prop_ = init_s_sky_prop
    skip_nonlin_ = .false.; if (present(skip_nonlin)) skip_nonlin_ = skip_nonlin 
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
    allocate(sd%mask(sd%ntod, sd%ndet))
    allocate(sd%pix(sd%ntod, sd%ndet, sd%nhorn))
    allocate(sd%psi(sd%ntod, sd%ndet, sd%nhorn))
    allocate(sd%flag(sd%ntod, sd%ndet))
    if (init_s_sky_prop_)    allocate(sd%s_sky_prop(sd%ntod, sd%ndet, 2:sd%ndelta))
    if (init_s_bp_prop_)     allocate(sd%s_bp_prop(sd%ntod, sd%ndet, 2:sd%ndelta))
    if (init_s_sky_prop_)    allocate(sd%mask2(sd%ntod, sd%ndet))
    if (tod%sample_mono)     allocate(sd%s_mono(sd%ntod, sd%ndet))
    if (tod%subtract_zodi)   allocate(sd%s_zodi(sd%ntod, sd%ndet))
    if (tod%apply_inst_corr) allocate(sd%s_inst(sd%ntod, sd%ndet))
    !call update_status(status, "todinit_alloc")

    !if (.true. .or. tod%myid == 78) write(*,*) 'c2', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    ! Decompress pointing, psi and flags for current scan
    call timer%start(TOD_DECOMP, tod%band)
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       call tod%decompress_pointing_and_flags(scan, j, sd%pix(:,j,:), &
            & sd%psi(:,j,:), sd%flag(:,j))
    end do
    
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
    else
       call project_sky(tod, map_sky(:,:,:,1), sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, &
            & procmask, scan, sd%s_sky, sd%mask)
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
    !call update_status(status, "todinit_bp")
    !if (.true. .or. tod%myid == 78) write(*,*) 'c71', tod%myid, tod%correct_sl
    !if (.true. .or. tod%myid == 78) write(*,*) 'c72', tod%myid, tod%ndet
    !if (.true. .or. tod%myid == 78) write(*,*) 'c73', tod%myid, tod%slconv(1)%p%psires

    ! Perform sanity tests
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (all(sd%mask(:,j) == 0)) tod%scans(scan)%d(j)%accept = .false.
       if (tod%scans(scan)%d(j)%N_psd%sigma0 <= 0.d0) tod%scans(scan)%d(j)%accept = .false.
    end do
    !call update_status(status, "todinit_sanity")
    !if (.true. .or. tod%myid == 78) write(*,*) 'c8', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
    
    ! Construct orbital dipole template
    call timer%start(TOD_ORBITAL, tod%band)
    call tod%construct_dipole_template(scan, sd%pix(:,:,1), sd%psi(:,:,1), sd%s_orb)
    call timer%stop(TOD_ORBITAL, tod%band)
    !if (.true. .or. tod%myid == 78) write(*,*) 'c9', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    ! Construct zodical light template
    if (tod%subtract_zodi) then
       call timer%start(TOD_ZODI, tod%band)
       call compute_zodi_template(tod%nside, sd%pix(:,:,1), tod%scans(scan)%satpos, tod%nu_c, sd%s_zodi)
       call timer%stop(TOD_ZODI, tod%band)
    end if
    !if (.true. .or. tod%myid == 78) write(*,*) 'c10', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

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
    end if
    !call update_status(status, "todinit_instcorr")

    ! Construct total sky signal
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       sd%s_tot(:,j) = sd%s_sky(:,j) + sd%s_sl(:,j) + sd%s_orb(:,j)
       if (tod%sample_mono) sd%s_tot(:,j) = sd%s_tot(:,j) + sd%s_mono(:,j)
       if (tod%apply_inst_corr) sd%s_tot(:,j) = sd%s_tot(:,j) + sd%s_inst(:,j)
    end do

    ! Apply non-linearity corrections
    if (.not. skip_nonlin_) call tod%apply_nonlin_corr_inst(scan, sd)
    
    !call update_status(status, "todinit_stot")

  end subroutine init_scan_data_singlehorn


  subroutine init_scan_data_differential(sd, tod, scan, map_sky, procmask, procmask2, &
       & init_s_bp, init_s_bp_prop, init_s_sky_prop, polang)
    implicit none
    class(comm_scandata),                      intent(inout)          :: sd    
    class(comm_tod),                           intent(inout)          :: tod
    integer(i4b),                              intent(in)             :: scan
    real(sp),          dimension(1:,1:,0:,1:), intent(in)             :: map_sky
    real(sp),          dimension(0:),          intent(in)             :: procmask
    real(sp),          dimension(0:),          intent(in)             :: procmask2
    logical(lgt),                              intent(in),   optional :: init_s_bp
    logical(lgt),                              intent(in),   optional :: init_s_bp_prop
    logical(lgt),                              intent(in),   optional :: init_s_sky_prop
    real(dp),                                  intent(in),   optional :: polang

    integer(i4b) :: j, k, ndelta
    logical(lgt) :: init_s_bp_, init_s_bp_prop_, init_s_sky_prop_
    real(sp),     allocatable, dimension(:,:)     :: s_bufA, s_bufB, s_buf2A, s_buf2B      ! Buffer

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
    allocate(sd%s_tot(sd%ntod, sd%ndet))
    allocate(sd%s_totA(sd%ntod, sd%ndet))
    allocate(sd%s_totB(sd%ntod, sd%ndet))
    allocate(sd%s_orbA(sd%ntod, sd%ndet))
    allocate(sd%s_orbB(sd%ntod, sd%ndet))
    allocate(sd%mask(sd%ntod, sd%ndet))
    allocate(sd%pix(sd%ntod, 1, sd%nhorn))
    allocate(sd%psi(sd%ntod, 1, sd%nhorn))
    allocate(sd%flag(sd%ntod, 1))
    if (init_s_sky_prop_)   allocate(sd%s_sky_prop(sd%ntod, sd%ndet, 2:sd%ndelta))
    if (init_s_bp_prop_)    allocate(sd%s_bp_prop(sd%ntod, sd%ndet, 2:sd%ndelta))
    if (init_s_sky_prop_)   allocate(sd%mask2(sd%ntod, sd%ndet))
    if (tod%sample_mono)    allocate(sd%s_mono(sd%ntod, sd%ndet))
    if (tod%subtract_zodi)  allocate(sd%s_zodi(sd%ntod, sd%ndet))
    sd%s_tot  = 0.
    sd%s_totA = 0.
    sd%s_totB = 0.
    sd%s_orb  = 0.
    sd%s_orbA = 0.
    sd%s_orbB = 0.

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
    else
       call tod%diode2tod_inst(scan, map_sky, procmask, sd%tod)
    end if

    ! Construct sky signal template
    if (init_s_bp_) then
       call project_sky_differential(tod, map_sky(:,:,:,1), sd%pix(:,1,:), sd%psi(:,1,:), sd%flag(:,1), &
            & procmask, scan, s_bufA, s_bufB, sd%mask, s_bpA=s_buf2A, s_bpB=s_buf2B)
    else
       call project_sky_differential(tod, map_sky(:,:,:,1), sd%pix(:,1,:), sd%psi(:,1,:), sd%flag(:,1), &
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
          call project_sky_differential(tod, map_sky(:,:,:,k), sd%pix(:,1,:), sd%psi(:,1,:), sd%flag(:,1), &
               & procmask, scan, s_bufA, s_bufB, sd%mask, s_bpA=s_buf2A, s_bpB=s_buf2B)
          do j = 1, sd%ndet
             if (.not. tod%scans(scan)%d(j)%accept) cycle
             sd%s_sky_prop(:,j,k) = (1.+tod%x_im(j))*s_bufA(:,j)  - (1.-tod%x_im(j))*s_bufB(:,j)
             if (init_s_bp_) sd%s_bp_prop(:,j,k)  = (1.+tod%x_im(j))*s_buf2A(:,j) - (1.-tod%x_im(j))*s_buf2B(:,j)
          end do
       end do
    else if (init_s_sky_prop_) then
       do k = 2, sd%ndelta
          call project_sky_differential(tod, map_sky(:,:,:,k), sd%pix(:,1,:), sd%psi(:,1,:), sd%flag(:,1), &
               & procmask, scan, s_bufA, s_bufB, sd%mask)
          do j = 1, sd%ndet
             if (.not. tod%scans(scan)%d(j)%accept) cycle
             sd%s_sky_prop(:,j,k) = (1.+tod%x_im(j))*s_bufA(:,j) - (1.-tod%x_im(j))*s_bufB(:,j)
          end do
       end do
    end if

    ! Perform sanity tests
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (all(sd%mask(:,j) == 0)) tod%scans(scan)%d(j)%accept = .false.
       if (tod%scans(scan)%d(j)%N_psd%sigma0 <= 0.d0) tod%scans(scan)%d(j)%accept = .false.
    end do
    
    ! Construct orbital dipole template
    call tod%construct_dipole_template_diff(scan, sd%pix(:,:,1), sd%psi(:,:,1), s_bufA, 1d3)
    call tod%construct_dipole_template_diff(scan, sd%pix(:,:,2), sd%psi(:,:,2), s_bufB, 1d3)
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
          call compute_zodi_template(tod%nside, sd%pix(:,1:1,1), tod%scans(scan)%satpos, tod%nu_c(j:j), s_bufA)
          call compute_zodi_template(tod%nside, sd%pix(:,1:1,2), tod%scans(scan)%satpos, tod%nu_c(j:j), s_bufB)
          sd%s_zodi(:,j) = (1.+tod%x_im(j))*s_bufA(:,j) - (1.-tod%x_im(j))*s_bufB(:,j)
          sd%s_tot(:,j)  = sd%s_tot(:,j) + sd%s_zodi(:,j)
          sd%s_totA(:,j) = sd%s_totA(:,j) + s_bufA(:,j)
          sd%s_totB(:,j) = sd%s_totB(:,j) + s_bufB(:,j)
       end do
    else
       sd%s_zodi = 0.
    end if

    ! Construct sidelobe template
    if (tod%correct_sl) then
       do j = 1, sd%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          !call tod%construct_sl_template(tod%slconv(1)%p, sd%pix(:,1,1), sd%psi(:,1,1), s_bufA(:,j), 1.5707963267948966192d0)
          !call tod%construct_sl_template(tod%slconv(3)%p, sd%pix(:,1,2), sd%psi(:,1,2), s_bufB(:,j), -1.5707963267948966192d0)
          call tod%construct_sl_template(tod%slconv(1)%p, sd%pix(:,1,1), sd%psi(:,1,1), s_bufA(:,j),  polang)
          call tod%construct_sl_template(tod%slconv(3)%p, sd%pix(:,1,2), sd%psi(:,1,2), s_bufB(:,j), -polang)
          sd%s_sl(:,j)  = (1d0+tod%x_im(j))*s_bufA(:,j) - (1d0-tod%x_im(j))*s_bufB(:,j)
          sd%s_tot(:,j) = sd%s_tot(:,j) + sd%s_sl(:,j)
          sd%s_totA(:,j) = sd%s_totA(:,j) + s_bufA(:,j)
          sd%s_totB(:,j) = sd%s_totB(:,j) + s_bufB(:,j)
       end do
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

    ! Clean-up
    deallocate(s_bufA, s_bufB, s_buf2A, s_buf2B)

  end subroutine init_scan_data_differential


  subroutine dealloc_scandata(sd)
    implicit none
    class(comm_scandata), intent(inout) :: sd    

    sd%ntod = -1; sd%ndet = -1; sd%nhorn = -1

    ! Deallocate data structures
    deallocate(sd%tod, sd%n_corr, sd%s_sl, sd%s_sky)
    deallocate(sd%s_orb, sd%s_tot, sd%mask)
    deallocate(sd%pix, sd%psi, sd%flag)
    if (allocated(sd%s_sky_prop))  deallocate(sd%s_sky_prop)
    if (allocated(sd%s_bp_prop))   deallocate(sd%s_bp_prop)
    if (allocated(sd%s_bp))        deallocate(sd%s_bp)
    if (allocated(sd%s_mono))      deallocate(sd%s_mono)
    if (allocated(sd%mask2))       deallocate(sd%mask2)
    if (allocated(sd%s_zodi))      deallocate(sd%s_zodi)
    if (allocated(sd%s_totA))      deallocate(sd%s_totA)
    if (allocated(sd%s_totB))      deallocate(sd%s_totB)
    if (allocated(sd%s_inst))      deallocate(sd%s_inst)
    if (allocated(sd%s_orbA))      deallocate(sd%s_orbA)
    if (allocated(sd%s_orbB))      deallocate(sd%s_orbB)

  end subroutine dealloc_scandata


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Sampling drivers etc.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sample_calibration(tod, mode, handle, map_sky, procmask, procmask2, polang, smooth)
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

    if (trim(mode) == 'abscal') then
       timer_id = TOD_ABSCAL
    else if (trim(mode) == 'relcal') then
       timer_id = TOD_RELCAL
    else if (trim(mode) == 'imbal') then
       timer_id = TOD_IMBAL
    else if (trim(mode) == 'deltaG') then
       timer_id = TOD_DELTAG
    end if

    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       if (tod%nhorn == 1) then
          call init_scan_data_singlehorn(sd, tod, i, map_sky, procmask, procmask2)
       else
          call init_scan_data_differential(sd, tod, i, map_sky, procmask, procmask2)
       end if

       ! Set up filtered calibration signal, conditional contribution and mask
       call timer%start(timer_id, tod%band)
       call tod%downsample_tod(sd%s_orb(:,1), ext)
       allocate(s_invsqrtN(ext(1):ext(2), tod%ndet))      ! s * invN
       allocate(s_buf(sd%ntod, sd%ndet))
       allocate(mask_lowres(ext(1):ext(2), tod%ndet))
       do j = 1, tod%ndet
          if (.not. tod%scans(i)%d(j)%accept) cycle
          call tod%downsample_tod(sd%mask(:,j), ext, mask_lowres(:,j), threshold=0.9)
          if (trim(mode) == 'abscal' .and. tod%orb_abscal) then
             ! Calibrator = orbital dipole only
             call tod%downsample_tod(sd%s_orb(:,j), ext, s_invsqrtN(:,j))
          else if (trim(mode) == 'imbal' .and. tod%orb_abscal) then
             ! Calibrator = common mode signal
             ! Jarosik uses the orbital dipole for this.
             s_buf(:,j) = tod%scans(i)%d(j)%gain*(sd%s_orbA(:,j) + sd%s_orbB(:,j))
             call fill_all_masked(s_buf(:,j), sd%mask(:,j), sd%ntod, .false., &
               & real(tod%scans(i)%d(j)%N_psd%sigma0, sp), handle, tod%scans(i)%chunk_num)
             call tod%downsample_tod(s_buf(:,j), ext, s_invsqrtN(:,j))
          else if (trim(mode) == 'imbal' .and. .not. tod%orb_abscal) then
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
       call multiply_inv_N(tod, i, s_invsqrtN, sampfreq=tod%samprate_lowres, pow=0.5d0)

       if (trim(mode) == 'abscal' .or. trim(mode) == 'relcal' .or. trim(mode) == 'imbal') then
          ! Constant gain terms; accumulate contribution from this scan
          do j = 1, tod%ndet
             if (.not. tod%scans(i)%d(j)%accept) cycle
             if (trim(mode) == 'abscal' .and. tod%orb_abscal) then
                s_buf(:,j) = real(tod%gain0(0),sp) * (sd%s_tot(:,j) - sd%s_orb(:,j)) + &
                     & real(tod%gain0(j) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j)
             else if (trim(mode) == 'abscal' .and. .not. tod%orb_abscal) then
                s_buf(:,j) = real(tod%gain0(j) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j)
             else if (trim(mode) == 'relcal') then
                s_buf(:,j) = real(tod%gain0(0) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j)
             else if (trim(mode) == 'imbal' .and. tod%orb_abscal) then
                 s_buf(:,j) = real(tod%scans(i)%d(j)%gain,sp) * (  &
             &   sd%s_totA(:,j) - sd%s_totB(:,j) + &
             &   real(tod%x_im(j),sp)*(sd%s_totA(:,j) + sd%s_totB(:,j)   &
             &                       -(sd%s_orbA(:,j) + sd%s_orbB(:,j))) &
             &   )
             else if (trim(mode) == 'imbal' .and. .not. tod%orb_abscal) then
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

       ! Clean up
       call wall_time(t2)
       tod%scans(i)%proctime   = tod%scans(i)%proctime   + t2-t1
       tod%scans(i)%n_proctime = tod%scans(i)%n_proctime + 1
       call dealloc_scandata(sd)
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

    call timer%stop(timer_id, tod%band)

  end subroutine sample_calibration


  ! Sample baseline
  subroutine sample_baseline(tod, handle, map_sky, procmask, procmask2)
    implicit none
    class(comm_tod),                              intent(inout) :: tod
    type(planck_rng),                             intent(inout) :: handle
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
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
          call init_scan_data_singlehorn(sd, tod, i, map_sky, procmask, procmask2)
       else
          call init_scan_data_differential(sd, tod, i, map_sky, procmask, procmask2)
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
       call dealloc_scandata(sd)
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
       else if (abs(tod%scans(scan)%d(j)%chisq) > tod%chisq_threshold .or. &  ! Discard scans with high chisq or NaNs
            & isNaN(tod%scans(scan)%d(j)%chisq)) then
          write(*,fmt='(a,i8,i5,a,f12.1)') 'Reject scan, det = ', &
               & tod%scanid(scan), j, ', chisq = ', tod%scans(scan)%d(j)%chisq
          tod%scans(scan)%d(j)%accept = .false.
       end if
    end do
   !  if (any(.not. tod%scans(scan)%d%accept)) tod%scans(scan)%d%accept = .false. ! Do we actually want this..?
    do j = 1, ndet
        if (.not. tod%scans(scan)%d(j)%accept .and. tod%partner(j) >= 0) tod%scans(scan)%d(tod%partner(j))%accept = .false.
    end do

  end subroutine remove_bad_data

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
       call tod%compute_chisq(scan, j, sd%mask2(:,j), sd%s_sky(:,j), &
            & s_buf(:,j), sd%n_corr(:,j), sd%tod(:,j), absbp=.true.)
       chisq(j,1) = chisq(j,1) + tod%scans(scan)%d(j)%chisq_prop
       do k = 2, tod%n_bp_prop+1
          call tod%compute_chisq(scan, j, sd%mask2(:,j), sd%s_sky_prop(:,j,k), &
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
    !
    implicit none
    class(comm_tod),                       intent(in)   :: tod
    integer(i4b),                          intent(in)   :: scan
    type(comm_scandata),                   intent(in)   :: sd
    real(sp),            dimension(:,:,:), intent(out)  :: d_calib
    real(sp), dimension(:,:), intent(in), optional      :: jump_template

    integer(i4b) :: i, j, nout
    real(dp)     :: inv_gain

    call timer%start(TOD_MAPBIN, tod%band)
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
        if (present(jump_template)) d_calib(1,:,j) = d_calib(1,:,j) - jump_template(:,j) * inv_gain
       end if
       if (tod%output_n_maps > 1) d_calib(2,:,j) = d_calib(1,:,j) - sd%s_sky(:,j) + sd%s_bp(:,j)              ! residual
       if (tod%output_n_maps > 2) d_calib(3,:,j) = (sd%n_corr(:,j) - sum(sd%n_corr(:,j)/sd%ntod)) * inv_gain  ! ncorr
       if (tod%output_n_maps > 3) d_calib(4,:,j) = sd%s_bp(:,j)                                               ! bandpass
       if (tod%output_n_maps > 4) d_calib(5,:,j) = sd%s_orb(:,j)                                              ! orbital dipole
       if (tod%output_n_maps > 5) d_calib(6,:,j) = sd%s_sl(:,j)                                               ! sidelobes
       if (tod%output_n_maps > 6) then
          if (allocated(sd%s_zodi)) then
             d_calib(7,:,j) = sd%s_zodi(:,j)                                                     ! zodi
          else
             d_calib(7,:,j) = 0.
          end if
       end if
       if (tod%output_n_maps > 7) d_calib(8,:,j) = sd%s_inst(:,j)                                               ! instrument specific
       
       !Bandpass proposals
       do i = 1, nout-tod%output_n_maps
          d_calib(tod%output_n_maps+i,:,j) = d_calib(1,:,j) + sd%s_bp(:,j) - sd%s_bp_prop(:,j,i+1)
       end do

    end do
    call timer%stop(TOD_MAPBIN, tod%band)

  end subroutine compute_calibrated_data

  subroutine distribute_sky_maps(tod, map_in, scale, map_out, map_full)
    implicit none
    class(comm_tod),                       intent(in)     :: tod
    type(map_ptr), dimension(1:,1:),       intent(inout)  :: map_in       ! (ndet,ndelta)    
    real(sp),                              intent(in)     :: scale
    real(sp),      dimension(1:,1:,0:,1:), intent(out)    :: map_out
    real(dp),      dimension(0:), intent(out), optional   :: map_full

    integer(i4b) :: i, j, k, l, npix, nmaps
    real(dp),     allocatable, dimension(:,:) :: m_buf

    npix  = map_in(1,1)%p%info%npix
    nmaps = map_in(1,1)%p%info%nmaps
    allocate(m_buf(0:npix-1,nmaps))
    do j = 1, size(map_in,2)
       do i = 1, size(map_in,1)
          map_in(i,j)%p%map = scale * map_in(i,j)%p%map ! unit conversion
          call map_in(i,j)%p%bcast_fullsky_map(m_buf)
          do k = 1, tod%nobs
             map_out(:,k,i,j) = m_buf(tod%ind2pix(k),:)
          end do
          if (j == 1 .and. present(map_full)) map_full = map_full + m_buf(:,1)
       end do
       do k = 1, tod%nobs
          do l = 1, tod%nmaps
             map_out(l,k,0,j) = sum(map_out(l,k,1:tod%ndet,j))/tod%ndet
          end do
       end do
    end do
    deallocate(m_buf)

  end subroutine distribute_sky_maps


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
    ! HDF5 variables
    character(len=512) :: mystring, mysubstring !< dummy values for string manipulation
    integer(i4b)       :: myindex     !< dummy value for string manipulation
    character(len=512) :: currentHDFFile !< hdf5 file which stores simulation output
    character(len=6)   :: pidLabel
    character(len=3)   :: detectorLabel
    type(hdf_file)     :: hdf5_file   !< hdf5 file to work with
    integer(i4b)       :: hdf5_error  !< hdf5 error status
    integer(HID_T)     :: hdf5_file_id !< File identifier
    integer(HID_T)     :: dset_id     !< Dataset identifier
    integer(HSIZE_T), dimension(1) :: dims
    ! Other variables
    integer(i4b)                          :: i, j, k !< loop variables
    integer(i4b)       :: mpi_err !< MPI error status
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
        ! getting gain for each detector (units, V / K)
        ! (gain is assumed to be CONSTANT for EACH SCAN)
        gain   = self%scans(scan_id)%d(j)%gain
        !write(*,*) "gain ", gain
        sigma0 = self%scans(scan_id)%d(j)%N_psd%sigma0
        !write(*,*) "sigma0 ", sigma0
        ! Simulating tods
        !tod_per_detector(i,j) = n_corr(i, j) + sigma0 * rand_gauss(handle)
        tod_per_detector(i,j) = gain * s_tot(i,j) + n_corr(i, j) + sigma0 * rand_gauss(handle)
        !tod_per_detector(i,j) = gain * s_tot(i,j) + sigma0 * rand_gauss(handle)
        !tod_per_detector(i,j) = sigma0 * rand_gauss(handle)
        !tod_per_detector(i,j) = 0
      end do
    end do

    !write(*,*) 'a', self%scanid(scan_id), self%scans(scan_id)%d(1)%N_psd%sigma0, (sum((tod_per_detector(:,1)/self%scans(scan_id)%d(1)%N_psd%sigma0)**2)/ntod-1)/sqrt(2./ntod)

    !----------------------------------------------------------------------------------
    ! Saving stuff to hdf file
    ! Getting the full path and name of the current hdf file to overwrite
    !----------------------------------------------------------------------------------
    mystring = trim(self%hdfname(scan_id))
    mysubstring = 'LFI_0'
    myindex = index(trim(mystring), trim(mysubstring))
    call get_tokens(trim(mystring), "/", toks=toks, num=ntoks)
    currentHDFFile = trim(self%sims_output_dir)//'/'//trim(toks(ntoks))
    !write(*,*) "hdf5name "//trim(self%hdfname(scan_id))
    !write(*,*) "currentHDFFile "//trim(currentHDFFile)
    ! Converting PID number into string value
    call int2string(self%scanid(scan_id), pidLabel)
    call int2string(self%myid, processor_label)
    !write(*,*) "!  Process: "//trim(processor_label)//" started writing PID: "//trim(pidLabel)//", into:"
    write(*,*) "!  Process:", self%myid, "started writing PID: "//trim(pidLabel)//", into:"
    write(*,*) "!  "//trim(currentHDFFile)
    ! For debugging
    !call MPI_Finalize(mpi_err)
    !stop

    dims(1) = ntod
    ! Initialize FORTRAN interface.
    call h5open_f(hdf5_error)
    ! Open an existing file - returns hdf5_file_id
    call  h5fopen_f(currentHDFFile, H5F_ACC_RDWR_F, hdf5_file_id, hdf5_error)
    if (hdf5_error /= 0) call h5eprint_f(hdf5_error)
    do j = 1, ndet
      detectorLabel = self%label(j)
      ! Open an existing dataset.
      call h5dopen_f(hdf5_file_id, trim(pidLabel)//'/'//trim(detectorLabel)//'/'//'tod', dset_id, hdf5_error)
      ! Write tod data to a dataset
      call h5dwrite_f(dset_id, H5T_IEEE_F32LE, tod_per_detector(:,j), dims, hdf5_error)
      ! Close the dataset.
      call h5dclose_f(dset_id, hdf5_error)
    end do
    ! Close the file.
    call h5fclose_f(hdf5_file_id, hdf5_error)
    ! Close FORTRAN interface.
    call h5close_f(hdf5_error)


    !write(*,*) "hdf5_error",  hdf5_error
    ! freeing memory up
    deallocate(tod_per_detector)
    write(*,*) "!  Process:", self%myid, "finished writing PID: "//trim(pidLabel)//"."

    ! lastly, we need to copy an existing filelist.txt into simulation folder
    ! and change the pointers to new files
    !if (self%myid == 0) then
    !  call system("cp "//trim(filelist)//" "//trim(simsdir))
    !  !mystring = filelist
    !  !mysubstring = ".txt"
    !  !myindex = index(trim(mystring), trim(mysubstring))
    !end if

    ! For debugging
    !call MPI_Finalize(mpi_err)
    !stop
  end subroutine simulate_tod

end module comm_tod_driver_mod
