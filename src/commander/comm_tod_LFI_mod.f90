module comm_tod_LFI_mod
  use comm_tod_mod
  use comm_param_mod
  use comm_map_mod
  implicit none

  private
  public comm_tod_LFI

  type, extends(comm_tod) :: comm_tod_LFI 

   contains
     procedure     :: compute_binned_map => compute_binned_LFI_map
     procedure     :: fit_gain           => fit_LFI_gain
     procedure     :: fit_n_corr         => fit_LFI_n_corr
  end type comm_tod_LFI

  interface comm_tod_LFI
     procedure constructor
  end interface comm_tod_LFI

contains

  !**************************************************
  !             Constructor
  !**************************************************
  function constructor(cpar)
    implicit none
    type(comm_params),       intent(in) :: cpar
    class(comm_tod_LFI),     pointer    :: constructor

    integer(i4b) :: i, j

    write(*,*) 'init tod_LFI'

    ! Test code, just to be able to read a single file; need to settle on parameter structure
    allocate(constructor)
    constructor%nscan    = 1
    constructor%nmaps    = 3
    constructor%ndet     = 1
    constructor%nhorn    = 1
    constructor%samprate = 50.d0
    allocate(constructor%stokes(constructor%nmaps))
    allocate(constructor%w(constructor%nmaps,constructor%nhorn,constructor%ndet))
    allocate(constructor%scans(constructor%nscan))
    constructor%stokes = [1,2,3]
    constructor%w      = 1.d0
    
    do i = 1, constructor%nscan
       call constructor%scans(i)%initHDF("filename.hdf", i, constructor%ndet)
    end do

  end function constructor


  ! Compute map with white noise assumption from correlated noise 
  ! corrected and calibrated data, d' = (d-n_corr-n_temp)/gain 
  subroutine compute_binned_LFI_map(self, map, rms)
    implicit none
    class(comm_tod_LFI), intent(in)  :: self
    class(comm_map),     intent(out) :: map, rms

    map%map = 0.d0
    rms%map = 0.d0

  end subroutine compute_binned_LFI_map
  
  ! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  subroutine fit_LFI_gain(self, map)
    implicit none
    class(comm_tod_LFI), intent(inout)  :: self
    class(comm_map),     intent(in)     :: map

  end subroutine fit_LFI_gain

  ! Compute correlated noise term, n_corr
  subroutine fit_LFI_n_corr(self, map)
    implicit none
    class(comm_tod_LFI), intent(inout)  :: self
    class(comm_map),     intent(in)     :: map

  end subroutine fit_LFI_n_corr

end module comm_tod_LFI_mod
