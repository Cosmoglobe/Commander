module comm_tod_mod
  use comm_utils
  use comm_param_mod
  use comm_hdf_mod
  use comm_map_mod
  implicit none

  private
  public comm_tod

  type :: comm_detscan
     character(len=10) :: label                           ! Detector label
     real(dp)          :: gain                            ! Gain; assumed constant over scan
     real(dp)          :: sigma0, alpha, fknee            ! Noise parameters
     real(dp)          :: samprate                        ! Sample rate in Hz
     real(dp),     allocatable, dimension(:)   :: tod     ! Detector values in time domain, (ntod)     
     real(dp),     allocatable, dimension(:)   :: n_corr  ! Correlated noise estimate in time domain
     real(dp),     allocatable, dimension(:)   :: n_temp  ! Template estimate in time domain
     integer(i4b), allocatable, dimension(:)   :: flag    ! Detector flag; 0 is accepted, /= 0 is rejected
     integer(i4b), allocatable, dimension(:,:) :: pix     ! Pixelized pointing, (ntod,nhorn)
     real(dp),     allocatable, dimension(:,:) :: psi     ! Polarization angle, (ntod,nhorn)
     real(dp),     allocatable, dimension(:)   :: N_psd   ! Noise power spectrum density; in uncalibrated units
     real(dp),     allocatable, dimension(:)   :: nu_psd  ! Noise power spectrum bins; in Hz
  end type comm_detscan

  type :: comm_scan
     integer(i4b) :: ntod                                        ! Number of time samples
     real(dp)     :: v_sun(3)                                    ! Observatory velocity relative to Sun in km/s
     real(dp),            allocatable, dimension(:)     :: time  ! MJD for current scan
     class(comm_detscan), allocatable, dimension(:)     :: d     ! Array of all detectors
   contains
     procedure :: initHDF => read_hdf_scan
  end type comm_scan

  type :: comm_tod 
     integer(i4b) :: nmaps                                      ! Number of Stokes parameters
     integer(i4b) :: ndet                                       ! Number of active detectors
     integer(i4b) :: nhorn                                      ! Number of horns
     integer(i4b) :: nscan                                      ! Number of scans
     real(dp)     :: samprate                                   ! Sample rate in Hz
     integer(i4b),     allocatable, dimension(:)     :: stokes  ! List of Stokes parameters
     real(dp),         allocatable, dimension(:,:,:) :: w       ! Stokes weights per detector per horn, (nmaps,nhorn,ndet)
     class(comm_scan), allocatable, dimension(:)     :: scans   ! Array of all scans
   contains
     procedure     :: compute_binned_map
     procedure     :: fit_gain
     procedure     :: fit_n_corr
  end type comm_tod

  interface comm_tod
     procedure constructor
  end interface comm_tod

contains

  !**************************************************
  !             Constructor
  !**************************************************
  function constructor(cpar)
    implicit none
    type(comm_params),       intent(in) :: cpar
    class(comm_tod),         pointer    :: constructor

    allocate(constructor)

  end function constructor


  !**************************************************
  !             Main routines
  !**************************************************

  ! Compute map with white noise assumption from correlated noise 
  ! corrected and calibrated data, d' = (d-n_corr-n_temp)/gain 
  subroutine compute_binned_map(self, map, rms)
    implicit none
    class(comm_tod), intent(in)  :: self
    class(comm_map), intent(out) :: map, rms

    map%map = 0.d0
    rms%map = 0.d0

  end subroutine compute_binned_map

  ! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  subroutine fit_gain(self, map)
    implicit none
    class(comm_tod), intent(inout)  :: self
    class(comm_map), intent(in)     :: map

  end subroutine fit_gain

  ! Compute correlated noise term, n_corr
  subroutine fit_n_corr(self, map)
    implicit none
    class(comm_tod), intent(inout)  :: self
    class(comm_map), intent(in)     :: map

  end subroutine fit_n_corr

  !**************************************************
  !             Utility routines
  !**************************************************

  subroutine read_hdf_scan(self, filename, scan, ndet)
    implicit none
    class(comm_scan)                   :: self
    character(len=*),       intent(in) :: filename
    integer(i4b),           intent(in) :: scan, ndet

    integer(i4b)     :: i, j, n, nhorn, ext(2)
    character(len=6) :: slabel, dlabel
    type(hdf_file)   :: file

    call int2string(scan, slabel)

    call open_hdf_file(filename, file, "r")

    ! Find array sizes
    call get_size_hdf(file, slabel // "/000001/pix", ext)
    nhorn     = ext(1)
    n         = ext(2)
    self%ntod = n

    ! Read common scan data
    allocate(self%time(n))
    call read_hdf(file, slabel // "/vsun", self%v_sun)
    call read_hdf(file, slabel // "/time", self%time)

    ! Read detector scans
    allocate(self%d(ndet))
    do i = 1, ndet
       allocate(self%d(i)%tod(n), self%d(i)%flag(n), self%d(i)%pix(nhorn,n), self%d(i)%psi(nhorn,n))
       call int2string(i, dlabel)
       call read_hdf(file, slabel // "/" // dlabel // "/label",  self%d(i)%label)
       call read_hdf(file, slabel // "/" // dlabel // "/gain",   self%d(i)%gain)
       call read_hdf(file, slabel // "/" // dlabel // "/sigma0", self%d(i)%sigma0)
       call read_hdf(file, slabel // "/" // dlabel // "/alpha",  self%d(i)%alpha)
       call read_hdf(file, slabel // "/" // dlabel // "/fknee",  self%d(i)%fknee)
       call read_hdf(file, slabel // "/" // dlabel // "/tod",    self%d(i)%tod)
       call read_hdf(file, slabel // "/" // dlabel // "/flag",   self%d(i)%flag)
       call read_hdf(file, slabel // "/" // dlabel // "/pix",    self%d(i)%pix)
       call read_hdf(file, slabel // "/" // dlabel // "/psi",    self%d(i)%psi)
    end do

  end subroutine read_hdf_scan

end module comm_tod_mod
