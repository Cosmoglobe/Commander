module comm_B_bl_mod
  use comm_param_mod
  use comm_map_mod
  use comm_B_mod
  use comm_utils
  implicit none

  private
  public comm_B_bl

  type, extends (comm_B) :: comm_B_bl
   contains
     ! Data procedures
     procedure :: B   => matmulB
     procedure :: B_t => matmulBt
  end type comm_B_bl

  interface comm_B_bl
     procedure constructor
  end interface comm_B_bl

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id)
    implicit none
    type(comm_params),                  intent(in) :: cpar
    type(comm_mapinfo), target,         intent(in) :: info
    integer(i4b),                       intent(in) :: id
    class(comm_B_bl),                  pointer    :: constructor

    character(len=512) :: dir
    
    ! General parameters
    allocate(constructor)
    dir = trim(cpar%datadir) // '/'

    ! Component specific parameters
    constructor%type  =  cpar%ds_beamtype(id)
    constructor%info  => info
    call read_beam(constructor%info%lmax, constructor%info%nmaps, constructor%b_l, &
         & trim(dir)//trim(cpar%ds_blfile(id)), trim(dir)//trim(cpar%ds_pixwin(id)))

  end function constructor
  
  subroutine matmulB(self, map, res)
    implicit none
    class(comm_B_bl), intent(in)    :: self
    class(comm_map),  intent(in)    :: map
    class(comm_map),  intent(inout) :: res
  end subroutine matmulB
  
  subroutine matmulBt(self, map, res)
    implicit none
    class(comm_B_bl), intent(in)    :: self
    class(comm_map),  intent(in)    :: map
    class(comm_map),  intent(inout) :: res
  end subroutine matmulBt

end module comm_B_bl_mod
