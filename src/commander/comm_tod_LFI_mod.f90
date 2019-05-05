module comm_tod_LFI_mod
  use comm_tod_mod
  use comm_param_mod
  use comm_map_mod
  implicit none

  private
  public comm_LFI_tod

  type, extends(comm_tod) :: comm_LFI_tod 

   contains
     procedure     :: process_tod        => process_LFI_tod
     procedure     :: compute_binned_map
     procedure     :: fit_gain
     procedure     :: fit_n_corr
  end type comm_LFI_tod

  interface comm_LFI_tod
     procedure constructor
  end interface comm_LFI_tod

contains

  !**************************************************
  !             Constructor
  !**************************************************
  function constructor(cpar)
    implicit none
    type(comm_params),       intent(in) :: cpar
    class(comm_LFI_tod),     pointer    :: constructor

    ! Set up LFI specific parameters
    allocate(constructor)
    constructor%myid     = cpar%myid
    constructor%numprocs = cpar%numprocs

    ! Test code, just to be able to read a single file; need to settle on parameter structure
    call constructor%get_scan_ids("data/filelist_5files.txt")

    constructor%nmaps    = 3
    constructor%ndet     = 4
    constructor%nhorn    = 1
    constructor%samprate = 50.d0
    allocate(constructor%stokes(constructor%nmaps))
    allocate(constructor%w(constructor%nmaps,constructor%nhorn,constructor%ndet))
    constructor%stokes = [1,2,3]
    constructor%w      = 1.d0

    ! Read the actual TOD
    call constructor%read_tod

  end function constructor

  !**************************************************
  !             Driver routine
  !**************************************************
  subroutine process_LFI_tod(self, map_in, map_out, rms_out)
    implicit none
    class(comm_LFI_tod),               intent(in)  :: self
    class(comm_map),     dimension(:), intent(in)  :: map_in            ! One map per detector
    class(comm_map),                   intent(out) :: map_out, rms_out  ! Combined output map and rms

    integer(i4b) :: i, j, ntod, ndet, nside, npix, nmaps
    real(sp), allocatable, dimension(:,:) :: n_corr, s_sl, d_calib, s_sky, s_orb
    real(sp), allocatable, dimension(:,:) :: map_tot, rms_tot
    real(dp), allocatable, dimension(:)   :: A_abscal, b_abscal

    ! Set up full-sky map structures
    nside = self%scans(1)%d(1)%nside
    npix  = 12*nside**2
    nmaps = self%nmaps
    allocate(map_tot(0:npix-1,nmaps), rms_tot(0:npix-1,nmaps))
    allocate(A_abscal(self%ndet), b_abscal(self%ndet))

    ! Compute output map and rms
    map_tot  = 0.d0
    rms_tot  = 0.d0
    A_abscal = 0.d0
    b_abscal = 0.d0
    do i = 1, self%nscan

       ! Short-cuts to local variables
       ntod = self%scans(i)%ntod
       ndet = self%ndet

       ! Set up local data structure for current scan
       allocate(d_calib(ntod, ndet))            ! Calibrated and cleaned TOD in uKcmb
       allocate(n_corr(ntod, ndet))             ! Correlated noise in V
       allocate(s_sl(ntod, ndet))               ! Sidelobe in uKcmb
       allocate(s_sky(ntod, ndet))              ! Stationary sky signal in uKcmb
       allocate(s_orb(ntod, ndet))              ! Orbital dipole in uKcmb

       ! --------------------
       ! Analyze current scan
       ! --------------------

       ! Construct sky signal template

       ! Construct orbital dipole template

       ! Construct sidelobe template

       ! Fit correlated noise

       ! Fit gain for current scan

       ! Compute contribution to absolute calibration from current scan

       ! Compute clean and calibrated TOD

       ! Compute binned map from cleaned TOD
       
       ! Clean up
       deallocate(n_corr, s_sl, s_sky, s_orb, d_calib)

    end do

    ! Compute absolute calibration, summed over all scans, and rescale output maps

    ! Copy rescaled maps to final output structure
    map_out%map = map_tot(map_out%info%pix,:)
    rms_out%map = rms_tot(rms_out%info%pix,:)

    deallocate(map_tot, rms_tot)

  end subroutine process_LFI_tod


  !**************************************************
  !             Sub-process routines
  !**************************************************

  ! Compute map with white noise assumption from correlated noise 
  ! corrected and calibrated data, d' = (d-n_corr-n_temp)/gain 
  subroutine compute_binned_map(self, map, rms)
    implicit none
    class(comm_LFI_tod), intent(in)  :: self
    class(comm_map),     intent(out) :: map, rms

    map%map = 0.d0
    rms%map = 0.d0

  end subroutine compute_binned_map
  
  ! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  subroutine fit_gain(self, map)
    implicit none
    class(comm_LFI_tod), intent(inout)  :: self
    class(comm_map),     intent(in)     :: map

  end subroutine fit_gain

  ! Compute correlated noise term, n_corr
  subroutine fit_n_corr(self, map)
    implicit none
    class(comm_LFI_tod), intent(inout)  :: self
    class(comm_map),     intent(in)     :: map

  end subroutine fit_n_corr

end module comm_tod_LFI_mod
