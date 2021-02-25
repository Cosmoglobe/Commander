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
  type :: scan_data
     integer(i4b) :: ntod, ndet, nhorn, ndelta
     real(sp),     allocatable, dimension(:,:)     :: tod        ! Raw data
     real(sp),     allocatable, dimension(:,:)     :: n_corr     ! Correlated noise in V
     real(sp),     allocatable, dimension(:,:,:)   :: s_sl       ! Sidelobe correction
     real(sp),     allocatable, dimension(:,:,:)   :: s_sky      ! Stationary sky signal
     real(sp),     allocatable, dimension(:,:,:,:) :: s_sky_prop ! Stationary sky signal proposal for bandpass sampling
     real(sp),     allocatable, dimension(:,:,:)   :: s_orb      ! Orbital dipole
     real(sp),     allocatable, dimension(:,:)     :: s_mono     ! Detector monopole correction 
     real(sp),     allocatable, dimension(:,:,:)   :: s_bp       ! Bandpass correction
     real(sp),     allocatable, dimension(:,:,:,:) :: s_bp_prop  ! Bandpass correction proposal     
     real(sp),     allocatable, dimension(:,:,:)   :: s_zodi     ! Zodiacal light
     real(sp),     allocatable, dimension(:,:,:)   :: s_tot      ! Total signal
     real(sp),     allocatable, dimension(:,:)     :: mask       ! TOD mask (flags + main processing mask)
     real(sp),     allocatable, dimension(:,:)     :: mask2      ! Small TOD mask, for bandpass sampling
     integer(i4b), allocatable, dimension(:,:,:)   :: pix        ! Discretized pointing 
     integer(i4b), allocatable, dimension(:,:,:)   :: psi        ! Discretized polarization angle
     integer(i4b), allocatable, dimension(:,:)     :: flag       ! Quality flags
     real(sp),     allocatable, dimension(:,:)     :: s_buf      ! Buffer
   contains
     procedure  :: initialize => init_scan_data
     procedure  :: dealloc    => dealloc_scan_data
  end type scan_data


contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Scan data routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_scan_data(self, tod, scan, map_sky, sprocmask, sprocmask2, &
       & init_s_bp, init_s_bp_prop, init_s_sky_prop)
    implicit none
    class(scan_data),                      intent(inout)          :: self    
    class(comm_tod),                       intent(inout)          :: tod
    integer(i4b),                          intent(in)             :: scan
    real(sp),          dimension(:,:,:,:), intent(in)             :: map_sky
    type(shared_1d_int),                   intent(in)             :: sprocmask
    type(shared_1d_int),                   intent(in)             :: sprocmask2
    logical(lgt),                          intent(in),   optional :: init_s_bp
    logical(lgt),                          intent(in),   optional :: init_s_bp_prop
    logical(lgt),                          intent(in),   optional :: init_s_sky_prop

    integer(i4b) :: j, k, ndelta
    logical(lgt) :: init_s_bp_, init_s_bp_prop_, init_s_sky_prop_

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
    allocate(self%s_sl(self%ntod, self%ndet, self%nhorn))
    allocate(self%s_sky(self%ntod, self%ndet, self%nhorn))
    allocate(self%s_bp(self%ntod, self%ndet, self%nhorn))
    allocate(self%s_orb(self%ntod, self%ndet, self%nhorn))
    allocate(self%s_tot(self%ntod, self%ndet, self%nhorn))
    allocate(self%s_buf(self%ntod, self%ndet))
    allocate(self%mask(self%ntod, self%ndet))
    allocate(self%pix(self%ntod, self%ndet, self%nhorn))
    allocate(self%psi(self%ntod, self%ndet, self%nhorn))
    allocate(self%flag(self%ntod, self%ndet))
    if (init_s_sky_prop_)   allocate(self%s_sky_prop(self%ntod, self%ndet, self%nhorn, 2:self%ndelta))
    if (init_s_bp_prop_)    allocate(self%s_bp_prop(self%ntod, self%ndet, self%nhorn, 2:self%ndelta))
    if (init_s_sky_prop_)   allocate(self%mask2(self%ntod, self%ndet))
    if (tod%sample_mono)    allocate(self%s_mono(self%ntod, self%ndet))
    if (tod%subtract_zodi)  allocate(self%s_zodi(self%ntod, self%ndet, self%nhorn))

    ! Decompress pointing, psi and flags for current scan
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       call tod%decompress_pointing_and_flags(scan, j, self%pix(:,j,:), &
            & self%psi(:,j,:), self%flag(:,j))
    end do
    call tod%symmetrize_flags(self%flag)
    
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
    do k = 1, self%nhorn
       if (init_s_bp_) then
          call project_sky(tod, map_sky(:,:,:,1), self%pix(:,:,k:k), self%psi(:,:,k:k), self%flag, &
               & sprocmask%a, scan, self%s_sky(:,:,k), self%mask, s_bp=self%s_bp(:,:,k))
       else
          call project_sky(tod, map_sky(:,:,:,1), self%pix(:,:,k:k), self%psi(:,:,k:k), self%flag, &
               & sprocmask%a, scan, self%s_sky(:,:,k), self%mask)
       end if
    end do

    ! Set up (optional) bandpass sampling quantities (s_sky_prop, mask2 and bp_prop)
    do k = 1, self%nhorn
       if (init_s_bp_prop_) then
          do j = 2, size(map_sky,4)
             call project_sky(tod, map_sky(:,:,:,j), self%pix(:,:,k:k), self%psi(:,:,k:k), self%flag, &
                  & sprocmask2%a, scan, self%s_sky_prop(:,:,k,j), self%mask2, s_bp=self%s_bp_prop(:,:,k,j))
          end do
       else if (init_s_sky_prop_) then
          do j = 2, size(map_sky,4)
             call project_sky(tod, map_sky(:,:,:,j), self%pix(:,:,k:k), self%psi(:,:,k:k), self%flag, &
                  & sprocmask2%a, scan, self%s_sky_prop(:,:,k,j), self%mask2)
          end do
       end if
    end do

    ! Perform sanity tests
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (all(self%mask(:,j) == 0)) tod%scans(scan)%d(j)%accept = .false.
       if (tod%scans(scan)%d(j)%sigma0 <= 0.d0) tod%scans(scan)%d(j)%accept = .false.
    end do
    
    ! Construct orbital dipole template
    do k = 1, tod%nhorn
       if (tod%orb_4pi_beam) then
          !call tod%orb_dp%p%compute_orbital_dipole_4pi(scan, self%pix(:,:,k), self%psi(:,:,k), self%s_orb(:,:,k))
       else
          ! Duncan: Fix units on your side to get rid of the 1d3 scaling factor
          !call tod%orb_dp%p%compute_orbital_dipole_pencil(scan, self%pix(:,:,k), self%psi(:,:,k), self%s_orb(:,:,k))
       end if
    end do

    ! Construct zodical light template
    if (tod%subtract_zodi) then
       do k = 1, self%nhorn
          call compute_zodi_template(tod%nside, self%pix(:,:,k), tod%scans(scan)%satpos, tod%nu_c, self%s_zodi(:,:,k))
       end do
    end if

    ! Construct sidelobe template
    if (tod%correct_sl) then
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          do k = 1, self%nhorn
             call tod%construct_sl_template(tod%slconv(j)%p, &
                  & self%pix(:,j,1), self%psi(:,j,1), self%s_sl(:,j,k), tod%polang(j))
             self%s_sl(:,j,k) = 2.d0 * self%s_sl(:,j,k) ! Scaling by a factor of 2, by comparison with LevelS. Should be understood
          end do
       end do
    else
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          self%s_sl(:,j,:) = 0.
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
       self%s_tot(:,j,:) = self%s_sky(:,j,:) + self%s_sl(:,j,:) + self%s_orb(:,j,:)
       do k = 1, self%nhorn
          if (tod%sample_mono) self%s_tot(:,j,k) = self%s_tot(:,j,k) + self%s_mono(:,j)
       end do
    end do

  end subroutine init_scan_data

  subroutine dealloc_scan_data(self)
    implicit none
    class(scan_data), intent(inout) :: self    

    self%ntod = -1; self%ndet = -1; self%nhorn = -1

    ! Deallocate data structures
    deallocate(self%tod, self%n_corr, self%s_sl, self%s_sky)
    deallocate(self%s_orb, self%s_tot, self%s_buf, self%mask)
    deallocate(self%pix, self%psi, self%flag)
    if (allocated(self%s_sky_prop))  deallocate(self%s_sky_prop)
    if (allocated(self%s_bp_prop))   deallocate(self%s_bp_prop)
    if (allocated(self%s_bp))        deallocate(self%s_bp)
    if (allocated(self%s_mono))      deallocate(self%s_mono)
    if (allocated(self%mask2))       deallocate(self%mask2)
    if (allocated(self%s_zodi))      deallocate(self%s_zodi)

  end subroutine dealloc_scan_data


end module comm_tod_driver_mod
