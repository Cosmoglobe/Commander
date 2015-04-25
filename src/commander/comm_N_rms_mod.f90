module comm_N_rms_mod
  use comm_N_mod
  use comm_param_mod
  implicit none

  private
  public comm_N_rms
  
  type, extends (comm_N) :: comm_N_rms
     real(dp), allocatable, dimension(:,:) :: iN, siN
   contains
     ! Data procedures
     procedure :: invN     => matmulInvN
     procedure :: sqrtInvN => matmulSqrtInvN
  end type comm_N_rms

  interface comm_N_rms
     procedure constructor
  end interface comm_N_rms

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, mask)
    implicit none
    type(comm_params),                  intent(in) :: cpar
    integer(i4b),                       intent(in) :: id
    real(dp),           dimension(:,:), intent(in) :: mask
    class(comm_N_rms),                  pointer    :: constructor

    character(len=512) :: dir
    
    ! General parameters
    allocate(constructor)
    dir = trim(cpar%datadir) // '/'

    ! Component specific parameters
    constructor%type  = cpar%ds_noise_format(id)
    constructor%nside = cpar%ds_nside(id)
    constructor%pol   = cpar%ds_polarization(id)
    constructor%nmaps = 1; if (constructor%pol) constructor%nmaps = 3

    call allocate_map(cpar%comm_chain, constructor%nside, constructor%nmaps, 0, &
         & constructor%rms, filename=trim(dir)//cpar%ds_noise_rms(id), np=constructor%np)
    allocate(constructor%iN(constructor%np,constructor%nmaps))
    constructor%iN = 1.d0 / constructor%rms**2 * mask
    constructor%iN = 1.d0 / constructor%rms    * mask
  end function constructor

  ! Return map_out = invN * map
  function matmulInvN(self, map)
    class(comm_N_rms),                                intent(in) :: self
    real(dp),          dimension(self%np,self%nmaps), intent(in) :: map
    real(dp),          dimension(self%np,self%nmaps)             :: matmulInvN
    matmulInvN = self%iN * map
  end function matmulInvN
  
  ! Return map_out = sqrtInvN * map
  function matmulSqrtInvN(self, map)
    class(comm_N_rms),                                intent(in) :: self
    real(dp),          dimension(self%np,self%nmaps), intent(in) :: map
    real(dp),          dimension(self%np,self%nmaps)             :: matmulSqrtInvN
    matmulSqrtInvN = self%siN * map
  end function matmulSqrtInvN

end module comm_N_rms_mod
