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
     procedure :: conv => matmulB
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
    constructor%lmax   = size(constructor%b_l,1)-1

  end function constructor
  
  function matmulB(self, alm_in, alm_out, trans, map)
    implicit none
    class(comm_B_bl), intent(in) :: self
    logical(lgt),     intent(in) :: alm_in, alm_out
    logical(lgt),     intent(in) :: trans
    class(comm_map),  intent(in) :: map
    class(comm_map),  pointer    :: matmulB

    integer(i4b) :: i, l

    matmulB => comm_map(map%info)
    if (.not. alm_in) then
       matmulB%map = map%map
       call matmulB%YtW
    else
       matmulB%alm = map%alm
    end if
    
    do i = 0, map%info%nalm-1
       l = map%info%lm(1,i)
       if (l <= self%lmax) then
          matmulB%alm(i,:) = matmulB%alm(i,:) * self%b_l(map%info%lm(1,i),:)
       end if
    end do

    if (.not. alm_out) call matmulB%Y

  end function matmulB
  
end module comm_B_bl_mod
