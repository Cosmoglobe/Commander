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
  use comm_shared_arr_mod
  implicit none


  ! Class for uncompressed data for a given scan
  type :: comm_scandata
     integer(i4b) :: ntod, ndet, nhorn, ndelta
     real(sp),     allocatable, dimension(:,:)     :: tod        ! Raw data
     real(sp),     allocatable, dimension(:,:)     :: n_corr     ! Correlated noise in V
     real(sp),     allocatable, dimension(:,:)     :: s_sl       ! Sidelobe correction
     real(sp),     allocatable, dimension(:,:)     :: s_sky      ! Stationary sky signal
     real(sp),     allocatable, dimension(:,:,:)   :: s_sky_prop ! Stationary sky signal proposal for bandpass sampling
     real(sp),     allocatable, dimension(:,:)     :: s_orb      ! Orbital dipole
     real(sp),     allocatable, dimension(:,:)     :: s_mono     ! Detector monopole correction 
     real(sp),     allocatable, dimension(:,:)     :: s_bp       ! Bandpass correction
     real(sp),     allocatable, dimension(:,:,:)   :: s_bp_prop  ! Bandpass correction proposal     
     real(sp),     allocatable, dimension(:,:)     :: s_zodi     ! Zodiacal light
     real(sp),     allocatable, dimension(:,:)     :: s_tot      ! Total signal
     real(sp),     allocatable, dimension(:,:)     :: mask       ! TOD mask (flags + main processing mask)
     real(sp),     allocatable, dimension(:,:)     :: mask2      ! Small TOD mask, for bandpass sampling
     integer(i4b), allocatable, dimension(:,:,:)   :: pix        ! Discretized pointing 
     integer(i4b), allocatable, dimension(:,:,:)   :: psi        ! Discretized polarization angle
     integer(i4b), allocatable, dimension(:,:)     :: flag       ! Quality flags

     real(sp),     allocatable, dimension(:,:)     :: s_totA     ! Total signal, horn A (differential only)
     real(sp),     allocatable, dimension(:,:)     :: s_totB     ! Total signal, horn B (differential only)
   contains
     procedure  :: init_singlehorn   => init_scan_data_singlehorn
     procedure  :: init_differential => init_scan_data_differential
     procedure  :: dealloc           => dealloc_scan_data
  end type comm_scandata


contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Scan data routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_scan_data_singlehorn(self, tod, scan, map_sky, procmask, procmask2, &
       & init_s_bp, init_s_bp_prop, init_s_sky_prop)
    implicit none
    class(comm_scandata),                      intent(inout)          :: self    
    class(comm_tod),                           intent(inout)          :: tod
    integer(i4b),                              intent(in)             :: scan
    real(sp),          dimension(1:,1:,0:,1:), intent(in)             :: map_sky
    real(sp),          dimension(0:),          intent(in)             :: procmask
    real(sp),          dimension(0:),          intent(in)             :: procmask2
    logical(lgt),                              intent(in),   optional :: init_s_bp
    logical(lgt),                              intent(in),   optional :: init_s_bp_prop
    logical(lgt),                              intent(in),   optional :: init_s_sky_prop

    integer(i4b) :: j, k, ndelta
    logical(lgt) :: init_s_bp_, init_s_bp_prop_, init_s_sky_prop_

    if (tod%nhorn /= 1) then
       write(*,*) 'Error: init_scan_data_singlehorn only applicable for 1-horn experiments'
       stop
    end if

    init_s_bp_ = .false.; if (present(init_s_bp)) init_s_bp_ = init_s_bp
    init_s_sky_prop_ = .false.; if (present(init_s_sky_prop)) init_s_sky_prop_ = init_s_sky_prop
    init_s_bp_prop_ = .false.
    if (present(init_s_bp_prop)) then
       init_s_bp_prop_ = init_s_bp_prop
       init_s_sky_prop_ = init_s_sky_prop
    end if

    self%ntod   = tod%scans(scan)%ntod
    self%ndet   = tod%ndet
    self%nhorn  = tod%nhorn
    self%ndelta = 0; if (present(init_s_sky_prop)) self%ndelta = size(map_sky,4)

    ! Allocate data structures
    allocate(self%tod(self%ntod, self%ndet))
    allocate(self%n_corr(self%ntod, self%ndet))
    allocate(self%s_sl(self%ntod, self%ndet))
    allocate(self%s_sky(self%ntod, self%ndet))
    allocate(self%s_bp(self%ntod, self%ndet))
    allocate(self%s_orb(self%ntod, self%ndet))
    allocate(self%s_tot(self%ntod, self%ndet))
    allocate(self%mask(self%ntod, self%ndet))
    allocate(self%pix(self%ntod, self%ndet, self%nhorn))
    allocate(self%psi(self%ntod, self%ndet, self%nhorn))
    allocate(self%flag(self%ntod, self%ndet))
    if (init_s_sky_prop_)   allocate(self%s_sky_prop(self%ntod, self%ndet, 2:self%ndelta))
    if (init_s_bp_prop_)    allocate(self%s_bp_prop(self%ntod, self%ndet, 2:self%ndelta))
    if (init_s_sky_prop_)   allocate(self%mask2(self%ntod, self%ndet))
    if (tod%sample_mono)    allocate(self%s_mono(self%ntod, self%ndet))
    if (tod%subtract_zodi)  allocate(self%s_zodi(self%ntod, self%ndet))

    ! Decompress pointing, psi and flags for current scan
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       call tod%decompress_pointing_and_flags(scan, j, self%pix(:,j,:), &
            & self%psi(:,j,:), self%flag(:,j))
    end do
    if (tod%symm_flags) call tod%symmetrize_flags(self%flag)
    
    ! Prepare TOD
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (tod%compressed_tod) then
          !call tod%decompress_tod(scan, j, self%tod(:,j))
       else
          self%tod(:,j) = tod%scans(scan)%d(j)%tod
       end if
    end do

    ! Construct sky signal template
    if (init_s_bp_) then
       call project_sky(tod, map_sky(:,:,:,1), self%pix(:,:,1), self%psi(:,:,1), self%flag, &
            & procmask, scan, self%s_sky, self%mask, s_bp=self%s_bp)
    else
       call project_sky(tod, map_sky(:,:,:,1), self%pix(:,:,1), self%psi(:,:,1), self%flag, &
            & procmask, scan, self%s_sky, self%mask)
    end if

    ! Set up (optional) bandpass sampling quantities (s_sky_prop, mask2 and bp_prop)
    if (init_s_bp_prop_) then
       do j = 2, size(map_sky,4)
          call project_sky(tod, map_sky(:,:,:,j), self%pix(:,:,1), self%psi(:,:,1), self%flag, &
               & procmask2, scan, self%s_sky_prop(:,:,j), self%mask2, s_bp=self%s_bp_prop(:,:,j))
       end do
    else if (init_s_sky_prop_) then
       do j = 2, size(map_sky,4)
          call project_sky(tod, map_sky(:,:,:,j), self%pix(:,:,1), self%psi(:,:,1), self%flag, &
               & procmask2, scan, self%s_sky_prop(:,:,j), self%mask2)
       end do
    end if

    ! Perform sanity tests
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (all(self%mask(:,j) == 0)) tod%scans(scan)%d(j)%accept = .false.
       if (tod%scans(scan)%d(j)%sigma0 <= 0.d0) tod%scans(scan)%d(j)%accept = .false.
    end do
    
    ! Construct orbital dipole template
    call tod%construct_dipole_template(scan, self%pix(:,:,1), self%psi(:,:,1), .true., self%s_orb)

    ! Construct zodical light template
    if (tod%subtract_zodi) then
       call compute_zodi_template(tod%nside, self%pix(:,:,1), tod%scans(scan)%satpos, tod%nu_c, self%s_zodi)
    end if

    ! Construct sidelobe template
    if (tod%correct_sl) then
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          call tod%construct_sl_template(tod%slconv(j)%p, &
               & self%pix(:,j,1), self%psi(:,j,1), self%s_sl(:,j), tod%polang(j))
          self%s_sl(:,j) = 2.d0 * self%s_sl(:,j) ! Scaling by a factor of 2, by comparison with LevelS. Should be understood
       end do
    else
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          self%s_sl(:,j) = 0.
       end do
    end if

    ! Construct monopole correction template
    if (tod%sample_mono) then
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          !self%s_mono(:,j) = tod%mono(j)
          self%s_mono(:,j) = 0.d0 ! Disabled for now
       end do
    end if

    ! Construct total sky signal
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       self%s_tot(:,j) = self%s_sky(:,j) + self%s_sl(:,j) + self%s_orb(:,j)
       if (tod%sample_mono) self%s_tot(:,j) = self%s_tot(:,j) + self%s_mono(:,j)
    end do

  end subroutine init_scan_data_singlehorn


  subroutine init_scan_data_differential(self, tod, scan, map_sky, procmask, procmask2, &
       & init_s_bp, init_s_bp_prop, init_s_sky_prop)
    implicit none
    class(comm_scandata),                      intent(inout)          :: self    
    class(comm_tod),                           intent(inout)          :: tod
    integer(i4b),                              intent(in)             :: scan
    real(sp),          dimension(1:,1:,0:,1:), intent(in)             :: map_sky
    real(sp),          dimension(0:),          intent(in)             :: procmask
    real(sp),          dimension(0:),          intent(in)             :: procmask2
    logical(lgt),                              intent(in),   optional :: init_s_bp
    logical(lgt),                              intent(in),   optional :: init_s_bp_prop
    logical(lgt),                              intent(in),   optional :: init_s_sky_prop

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
       init_s_bp_prop_ = init_s_bp_prop
       init_s_sky_prop_ = init_s_sky_prop
    end if

    self%ntod   = tod%scans(scan)%ntod
    self%ndet   = tod%ndet
    self%nhorn  = tod%nhorn
    self%ndelta = 0; if (present(init_s_sky_prop)) self%ndelta = size(map_sky,4)

    ! Allocate data structures
    allocate(self%tod(self%ntod, self%ndet))
    allocate(self%n_corr(self%ntod, self%ndet))
    allocate(self%s_sl(self%ntod, self%ndet))
    allocate(self%s_sky(self%ntod, self%ndet))
    allocate(self%s_bp(self%ntod, self%ndet))
    allocate(self%s_orb(self%ntod, self%ndet))
    allocate(self%s_tot(self%ntod, self%ndet))
    allocate(self%s_totA(self%ntod, self%ndet))
    allocate(self%s_totB(self%ntod, self%ndet))
    allocate(self%mask(self%ntod, self%ndet))
    allocate(self%pix(self%ntod, 1, self%nhorn))
    allocate(self%psi(self%ntod, 1, self%nhorn))
    allocate(self%flag(self%ntod, 1))
    if (init_s_sky_prop_)   allocate(self%s_sky_prop(self%ntod, self%ndet, 2:self%ndelta))
    if (init_s_bp_prop_)    allocate(self%s_bp_prop(self%ntod, self%ndet, 2:self%ndelta))
    if (init_s_sky_prop_)   allocate(self%mask2(self%ntod, self%ndet))
    if (tod%sample_mono)    allocate(self%s_mono(self%ntod, self%ndet))
    if (tod%subtract_zodi)  allocate(self%s_zodi(self%ntod, self%ndet))
    self%s_totA = 0.
    self%s_totB = 0.

    allocate(s_bufA(self%ntod, self%ndet))
    allocate(s_bufB(self%ntod, self%ndet))
    allocate(s_buf2A(self%ntod, self%ndet))
    allocate(s_buf2B(self%ntod, self%ndet))

    ! Decompress pointing, psi and flags for current scan
    call tod%decompress_pointing_and_flags(scan, j, self%pix(:,1,:), &
            & self%psi(:,1,:), self%flag(:,1))
    
    ! Prepare TOD
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (tod%compressed_tod) then
          !call tod%decompress_tod(scan, j, self%tod(:,j))
       else
          self%tod(:,j) = tod%scans(scan)%d(j)%tod
       end if
    end do

    ! Construct sky signal template
    if (init_s_bp_) then
       call project_sky_differential(tod, map_sky(:,:,:,1), self%pix(:,1,:), self%psi(:,1,:), self%flag(:,1), &
            & procmask, scan, s_bufA, s_bufB, self%mask, s_bpA=s_buf2A, s_bpB=s_buf2B)
    else
       call project_sky_differential(tod, map_sky(:,:,:,1), self%pix(:,1,:), self%psi(:,1,:), self%flag(:,1), &
            & procmask, scan, s_bufA, s_bufB, self%mask)
    end if
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       self%s_sky(:,j)  = (1.+tod%x_im(j))*s_bufA(:,j)  - (1.-tod%x_im(j))*s_bufB(:,j)
       self%s_totA(:,j) = self%s_totA(:,j) + s_bufA(:,j)
       self%s_totB(:,j) = self%s_totB(:,j) + s_bufB(:,j)
       self%s_tot(:,j)  = self%s_tot(:,j)  + self%s_sky(:,j)
       if (init_s_bp_) self%s_bp(:,j)  = (1.+tod%x_im(j))*s_buf2A(:,j) - (1.-tod%x_im(j))*s_buf2B(:,j)
    end do

    ! Set up (optional) bandpass sampling quantities (s_sky_prop, mask2 and bp_prop)
    if (init_s_bp_prop_) then
       do k = 2, size(map_sky,4)
          call project_sky_differential(tod, map_sky(:,:,:,k), self%pix(:,1,:), self%psi(:,1,:), self%flag(:,1), &
               & procmask, scan, s_bufA, s_bufB, self%mask, s_bpA=s_buf2A, s_bpB=s_buf2B)
          do j = 1, self%ndet
             if (.not. tod%scans(scan)%d(j)%accept) cycle
             self%s_sky_prop(:,j,k) = (1.+tod%x_im(j))*s_bufA(:,j)  - (1.-tod%x_im(j))*s_bufB(:,j)
             if (init_s_bp_) self%s_bp_prop(:,j,k)  = (1.+tod%x_im(j))*s_buf2A(:,j) - (1.-tod%x_im(j))*s_buf2B(:,j)
          end do
       end do
    else if (init_s_sky_prop_) then
       do k = 2, size(map_sky,4)
          call project_sky_differential(tod, map_sky(:,:,:,k), self%pix(:,1,:), self%psi(:,1,:), self%flag(:,1), &
               & procmask, scan, s_bufA, s_bufB, self%mask)
          do j = 1, self%ndet
             if (.not. tod%scans(scan)%d(j)%accept) cycle
             self%s_sky_prop(:,j,k) = (1.+tod%x_im(j))*s_bufA(:,j) - (1.-tod%x_im(j))*s_bufB(:,j)
          end do
       end do
    end if

    ! Perform sanity tests
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (all(self%mask(:,j) == 0)) tod%scans(scan)%d(j)%accept = .false.
       if (tod%scans(scan)%d(j)%sigma0 <= 0.d0) tod%scans(scan)%d(j)%accept = .false.
    end do
    
    ! Construct orbital dipole template
    call tod%construct_dipole_template(scan, self%pix(:,:,1), self%psi(:,:,1), .true., s_bufA)
    call tod%construct_dipole_template(scan, self%pix(:,:,2), self%psi(:,:,2), .true., s_bufB)
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       self%s_orb(:,j)  = (1.+tod%x_im(j))*s_bufA(:,j)  - (1.-tod%x_im(j))*s_bufB(:,j)
       self%s_totA(:,j) = self%s_totA(:,j) + s_bufA(:,j)
       self%s_totB(:,j) = self%s_totB(:,j) + s_bufB(:,j)
       self%s_tot(:,j)  = self%s_tot(:,j)  + self%s_orb(:,j)
    end do


    ! Construct zodical light template
    if (tod%subtract_zodi) then
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          call compute_zodi_template(tod%nside, self%pix(:,1:1,1), tod%scans(scan)%satpos, tod%nu_c(j:j), s_bufA)
          call compute_zodi_template(tod%nside, self%pix(:,1:1,2), tod%scans(scan)%satpos, tod%nu_c(j:j), s_bufB)
          self%s_zodi(:,j) = (1.+tod%x_im(j))*s_bufA(:,j) - (1.-tod%x_im(j))*s_bufB(:,j)
          self%s_tot(:,j)  = self%s_tot(:,j) + self%s_zodi(:,j)
          self%s_totA(:,j) = self%s_totA(:,j) + s_bufA(:,j)
          self%s_totB(:,j) = self%s_totB(:,j) + s_bufB(:,j)
       end do
    else
       self%s_zodi = 0.
    end if

    ! Construct sidelobe template
    if (tod%correct_sl) then
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          call tod%construct_sl_template(tod%slconv(1)%p, self%pix(:,1,1), self%psi(:,1,1), s_bufA(:,j), 0d0)
          call tod%construct_sl_template(tod%slconv(3)%p, self%pix(:,1,2), self%psi(:,1,2), s_bufB(:,j), 0d0)
          self%s_sl(:,j)  = 2.*((1d0+tod%x_im(j))*s_bufA(:,j) - (1d0-tod%x_im(j))*s_bufB(:,j))
          self%s_tot(:,j) = self%s_tot(:,j) + self%s_sl(:,j)
          self%s_totA(:,j) = self%s_totA(:,j) + 2.*s_bufA(:,j)
          self%s_totB(:,j) = self%s_totB(:,j) + 2.*s_bufB(:,j)
       end do
    else
       self%s_sl = 0.
    end if

    ! Construct monopole correction template
    if (tod%sample_mono) then
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          !self%s_mono(:,j) = tod%mono(j)
          self%s_mono(:,j) = 0.d0 ! Disabled for now
          self%s_tot(:,j)  = self%s_tot(:,j) + self%s_mono(:,j)
          self%s_totA(:,j) = self%s_totA(:,j) + self%s_mono(:,j)
          self%s_totB(:,j) = self%s_totB(:,j) + self%s_mono(:,j)
       end do
    else
       self%s_mono = 0.d0 
    end if

    ! Clean-up
    deallocate(s_bufA, s_bufB, s_buf2A, s_buf2B)

  end subroutine init_scan_data_differential


  subroutine dealloc_scan_data(self)
    implicit none
    class(comm_scandata), intent(inout) :: self    

    self%ntod = -1; self%ndet = -1; self%nhorn = -1

    ! Deallocate data structures
    deallocate(self%tod, self%n_corr, self%s_sl, self%s_sky)
    deallocate(self%s_orb, self%s_tot, self%mask)
    deallocate(self%pix, self%psi, self%flag)
    if (allocated(self%s_sky_prop))  deallocate(self%s_sky_prop)
    if (allocated(self%s_bp_prop))   deallocate(self%s_bp_prop)
    if (allocated(self%s_bp))        deallocate(self%s_bp)
    if (allocated(self%s_mono))      deallocate(self%s_mono)
    if (allocated(self%mask2))       deallocate(self%mask2)
    if (allocated(self%s_zodi))      deallocate(self%s_zodi)
    if (allocated(self%s_totA))      deallocate(self%s_totA)
    if (allocated(self%s_totB))      deallocate(self%s_totB)

  end subroutine dealloc_scan_data


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Sampling drivers etc.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Sample gain
  ! Supported modes = {abscal, relcal, deltaG}
  subroutine sample_calibration(tod, mode, handle, map_sky, procmask, procmask2)
    implicit none
    class(comm_tod),                              intent(inout) :: tod
    character(len=*),                             intent(in)    :: mode
    type(planck_rng),                             intent(inout) :: handle
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2

    integer(i4b) :: i, j, ext(2), ierr
    real(dp)     :: t1, t2
    real(dp), allocatable, dimension(:)   :: A, b
    real(sp), allocatable, dimension(:,:) :: s_invN, mask_lowres, s_buf
    real(dp), allocatable, dimension(:,:) :: dipole_mod
    type(comm_scandata) :: sd



    if (tod%myid == 0) then
       write(*,*) '   --> Sampling absolute calibration'
    end if

    if (trim(mode) == 'abscal' .or. trim(mode) == 'relcal') then
       allocate(A(tod%ndet), b(tod%ndet))
       A = 0.d0; b = 0.d0
    else if (trim(mode) == 'deltaG') then
       allocate(dipole_mod(tod%nscan_tot, tod%ndet))
       dipole_mod = 0.d0
    end if

    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       if (tod%nhorn == 1) then
          call sd%init_singlehorn(tod, i, map_sky, procmask, procmask2)
       else
          call sd%init_differential(tod, i, map_sky, procmask, procmask2)
       end if

       ! Set up filtered calibration signal, conditional contribution and mask
       call tod%downsample_tod(sd%s_orb(:,1), ext)
       allocate(s_invN(ext(1):ext(2), tod%ndet))      ! s * invN
       allocate(s_buf(sd%ntod, sd%ndet))
       allocate(mask_lowres(ext(1):ext(2), tod%ndet))
       do j = 1, tod%ndet
          if (.not. tod%scans(i)%d(j)%accept) cycle
          call tod%downsample_tod(sd%mask(:,j), ext, mask_lowres(:,j), threshold=0.9)
          if (trim(mode) == 'abscal' .and. tod%orb_abscal) then
             ! Calibrator = orbital dipole only
             call tod%downsample_tod(sd%s_orb(:,j), ext, s_invN(:,j))
          else
             ! Calibratior = total signal
             s_buf(:,j) = sd%s_tot(:,j)
             call fill_all_masked(s_buf(:,j), sd%mask(:,j), sd%ntod, .false., real(tod%scans(i)%d(j)%sigma0, sp), handle, tod%scans(i)%chunk_num)
             call tod%downsample_tod(s_buf(:,j), ext, s_invN(:,j))
          end if
       end do
       call multiply_inv_N(tod, i, s_invN, sampfreq=tod%samprate_lowres, pow=0.5d0)

       if (trim(mode) == 'abscal' .or. trim(mode) == 'relcal') then
          ! Constant gain terms; accumulate contribution from this scan
          do j = 1, tod%ndet
             if (.not. tod%scans(i)%d(j)%accept) cycle
             if (trim(mode) == 'abscal' .and. tod%orb_abscal) then
                s_buf(:,j) = real(tod%gain0(0),sp) * (sd%s_tot(:,j) - sd%s_orb(:,j)) + &
                     & real(tod%gain0(j) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j)
             else
                s_buf(:,j) = real(tod%gain0(j) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j)
             end if
             call accumulate_abscal(tod, i, sd%mask, s_buf, s_invN, s_invN, A, b, handle, .false., mask_lowres=mask_lowres)
          end do
       else
          ! Time-variable gain terms
          call calculate_gain_mean_std_per_scan(tod, i, s_invN, sd%mask, s_invN, sd%s_tot, handle, mask_lowres=mask_lowres)

          do j = 1, tod%ndet
             if (.not. tod%scans(i)%d(j)%accept) cycle
             dipole_mod(tod%scanid(i),j) = masked_variance(sd%s_sky(:,j), sd%mask(:,j))
          end do
       end if

       ! Clean up
       call wall_time(t2)
       tod%scans(i)%proctime   = tod%scans(i)%proctime   + t2-t1
       tod%scans(i)%n_proctime = tod%scans(i)%n_proctime + 1
       call sd%dealloc
       deallocate(s_invN, s_buf, mask_lowres)
    end do

    ! Perform sampling operations
    if (trim(mode) == 'abscal') then
       call sample_abscal_from_orbital(tod, handle, A, b)
    else if (trim(mode) == 'relcal') then
       call sample_relcal(tod, handle, A, b)
    else if (trim(mode) == 'deltaG') then
       call mpi_allreduce(mpi_in_place, dipole_mod, size(dipole_mod), MPI_DOUBLE_PRECISION, MPI_SUM, tod%info%comm, ierr)
       call sample_smooth_gain(tod, handle, dipole_mod)
    end if

    ! Clean up
    if (allocated(A))          deallocate(A)
    if (allocated(b))          deallocate(b)
    if (allocated(dipole_mod)) deallocate(dipole_mod)

  end subroutine sample_calibration

  subroutine remove_bad_data(tod, scan, flag)
    implicit none
    class(comm_tod),                   intent(inout) :: tod
    integer(i4b),    dimension(1:,1:), intent(in)    :: flag
    integer(i4b),                      intent(in)    :: scan

    integer(i4b) :: j, ntod, ndet
    
    ntod = size(flag,1)
    ndet = size(flag,2)
    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (count(iand(flag(:,j),tod%flag0) .ne. 0) > 0.1*ntod) then    ! Discard scans with less than 10% good data
          tod%scans(scan)%d(j)%accept = .false.
       else if (abs(tod%scans(scan)%d(j)%chisq) > tod%chisq_threshold .or. &  ! Discard scans with high chisq or NaNs
            & isNaN(tod%scans(scan)%d(j)%chisq)) then
          write(*,fmt='(a,i8,i5,a,f12.1)') 'Reject scan, det = ', &
               & tod%scanid(scan), j, ', chisq = ', tod%scans(scan)%d(j)%chisq
          tod%scans(scan)%d(j)%accept = .false.
       end if
    end do
    if (any(.not. tod%scans(scan)%d%accept)) tod%scans(scan)%d%accept = .false. ! Do we actually want this..?
    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) tod%scans(scan)%d(tod%partner(j))%accept = .false.
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

    allocate(s_buf(sd%ntod,sd%ndet))
    do j = 1, tod%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       s_buf(:,j) =  sd%s_sl(:,j) + sd%s_orb(:,j)
       call tod%compute_chisq(scan, j, sd%mask2(:,j), sd%s_sky(:,j), &
            & s_buf(:,j), sd%n_corr(:,j), absbp=.true.)
       chisq(j,1) = chisq(j,1) + tod%scans(scan)%d(j)%chisq_prop
       do k = 2, size(sd%s_sky_prop)
          call tod%compute_chisq(scan, j, sd%mask2(:,j), sd%s_sky_prop(:,j,k), &
               & s_buf(:,j), sd%n_corr(:,j), absbp=.true.)
          chisq(j,k) = chisq(j,k) + tod%scans(scan)%d(j)%chisq_prop
       end do
    end do
    deallocate(s_buf)

  end subroutine compute_chisq_abs_bp

  subroutine compute_calibrated_data(tod, scan, sd, d_calib)
    implicit none
    class(comm_tod),                       intent(in)   :: tod
    integer(i4b),                          intent(in)   :: scan
    type(comm_scandata),                   intent(in)   :: sd
    real(sp),            dimension(:,:,:), intent(out)  :: d_calib

    integer(i4b) :: j, nout
    real(sp)     :: inv_gain

    nout = size(d_calib,1)
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       inv_gain = 1.0 / real(tod%scans(scan)%d(j)%gain,sp)
       d_calib(1,:,j) = (tod%scans(scan)%d(j)%tod - sd%n_corr(:,j)) * inv_gain - sd%s_tot(:,j) + sd%s_sky(:,j) - sd%s_bp(:,j)
       if (nout > 1) d_calib(2,:,j) = d_calib(1,:,j) - sd%s_sky(:,j) + sd%s_bp(:,j)              ! residual
       if (nout > 2) d_calib(3,:,j) = (sd%n_corr(:,j) - sum(sd%n_corr(:,j)/sd%ntod)) * inv_gain  ! ncorr
       if (nout > 3) d_calib(4,:,j) = sd%s_bp(:,j)                                               ! bandpass
       if (nout > 4) d_calib(5,:,j) = sd%s_orb(:,j)                                              ! orbital dipole
       if (nout > 5) d_calib(6,:,j) = sd%s_sl(:,j)                                               ! sidelobes
       if (nout > 6) d_calib(7,:,j) = sd%s_zodi(:,j)                                             ! zodi
    end do

  end subroutine compute_calibrated_data

  subroutine distribute_sky_maps(tod, map_in, scale, map_out)
    implicit none
    class(comm_tod),                       intent(in)     :: tod
    type(map_ptr), dimension(1:,1:),       intent(inout)  :: map_in       ! (ndet,ndelta)    
    real(sp),                              intent(in)     :: scale
    real(sp),      dimension(1:,1:,0:,1:), intent(out)    :: map_out

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
       end do
       do k = 1, tod%nobs
          do l = 1, tod%nmaps
             map_out(l,k,0,j) = sum(map_out(l,k,1:tod%ndet,j))/tod%ndet
          end do
       end do
    end do
    deallocate(m_buf)

  end subroutine distribute_sky_maps

end module comm_tod_driver_mod
