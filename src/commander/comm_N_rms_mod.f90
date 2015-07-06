module comm_N_rms_mod
  use comm_N_mod
  use comm_param_mod
  use comm_map_mod
  implicit none

  private
  public comm_N_rms
  
  type, extends (comm_N) :: comm_N_rms
     class(comm_map), pointer :: siN
   contains
     ! Data procedures
     procedure :: invN     => matmulInvN_1map
     procedure :: sqrtInvN => matmulSqrtInvN_1map
     procedure :: rms      => returnRMS
  end type comm_N_rms

  interface comm_N_rms
     procedure constructor
  end interface comm_N_rms

!!$  interface matmulInvN
!!$     module procedure matmulInvN_1map, matmulInvN_2map
!!$  end interface matmulInvN
!!$  
!!$  interface matmulSqrtInvN
!!$     module procedure matmulSqrtInvN_1map, matmulSqrtInvN_2map
!!$  end interface matmulSqrtInvN
  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, mask)
    implicit none
    type(comm_params),                  intent(in) :: cpar
    type(comm_mapinfo), target,         intent(in) :: info
    integer(i4b),                       intent(in) :: id
    class(comm_map),                    intent(in) :: mask
    class(comm_N_rms),                  pointer    :: constructor

    character(len=512) :: dir
    
    ! General parameters
    allocate(constructor)
    dir = trim(cpar%datadir) // '/'

    ! Component specific parameters
    constructor%type    = cpar%ds_noise_format(id)
    constructor%nside   = info%nside
    constructor%nmaps   = info%nmaps
    constructor%np      = info%np
    constructor%pol     = info%nmaps == 3
    constructor%siN     => comm_map(info, trim(dir)//trim(cpar%ds_noise_rms(id)))
    constructor%siN%map = 1.d0 / constructor%siN%map

    ! Apply mask
    constructor%siN%map = constructor%siN%map * mask%map

    ! Set up diagonal covariance matrix in both pixel and harmonic space
    constructor%invN_diag     => comm_map(info)
    constructor%invN_diag%map = constructor%siN%map**2
    call compute_invN_lm(constructor%invN_diag)
    
  end function constructor

  ! Return map_out = invN * map
  subroutine matmulInvN_1map(self, map)
    implicit none
    class(comm_N_rms), intent(in)     :: self
    class(comm_map),   intent(inout)  :: map
    map%map = (self%siN%map)**2 * map%map
  end subroutine matmulInvN_1map
  
  ! Return map_out = sqrtInvN * map
  subroutine matmulSqrtInvN_1map(self, map)
    implicit none
    class(comm_N_rms), intent(in)    :: self
    class(comm_map),   intent(inout) :: map
    map%map = self%siN%map * map%map
  end subroutine matmulSqrtInvN_1map

  ! Return map_out = invN * map
  subroutine matmulInvN_2map(self, map, res)
    implicit none
    class(comm_N_rms), intent(in)     :: self
    class(comm_map),   intent(in)     :: map
    class(comm_map),   intent(inout)  :: res
    res%map = (self%siN%map)**2 * map%map
  end subroutine matmulInvN_2map
  
  ! Return map_out = sqrtInvN * map
  subroutine matmulSqrtInvN_2map(self, map, res)
    implicit none
    class(comm_N_rms), intent(in)    :: self
    class(comm_map),   intent(in)    :: map
    class(comm_map),   intent(inout) :: res
    res%map = self%siN%map * map%map
  end subroutine matmulSqrtInvN_2map

  ! Return RMS map
  subroutine returnRMS(self, res)
    implicit none
    class(comm_N_rms), intent(in)    :: self
    class(comm_map),   intent(inout) :: res
    where (self%siN%map > 0.d0)
       res%map = 1.d0/self%siN%map
    elsewhere
       res%map = infinity
    end where
  end subroutine returnRMS

end module comm_N_rms_mod
