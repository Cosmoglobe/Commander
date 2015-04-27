module comm_F_mod
  use comm_comp_mod
  use comm_param_mod
  use comm_map_mod
  implicit none

  private
  public comm_F

  type :: comm_F
     ! Linked list variables
     class(comm_F), pointer :: nextLink => null()
     class(comm_F), pointer :: prevLink => null()

     ! Data variables
     class(comm_comp), pointer :: c
     class(comm_map),  pointer :: F_diag
   contains
     ! Linked list procedures
     procedure :: next    ! get the link after this link
     procedure :: prev    ! get the link before this link
     procedure :: setNext ! set the link after this link
     procedure :: add     ! add new link at the end

     ! Data procedures
     procedure :: F => matmulF
  end type comm_F

  interface comm_F
     procedure constructor
  end interface comm_F

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, c)
    implicit none
    type(comm_params),   intent(in)         :: cpar
    integer(i4b),        intent(in)         :: id
    class(comm_comp),    intent(in), target :: c
    class(comm_F),       pointer            :: constructor

    ! General parameters
    allocate(constructor)
    
  end function constructor

  ! Return map_out = F * map
  subroutine matmulF(self, map, res)
    implicit none
    class(comm_F),   intent(in)     :: self
    class(comm_map), intent(in)     :: map
    class(comm_map), intent(inout)  :: res
    res%map = self%F_diag%map * map%map
  end subroutine matmulF

  function next(self)
    class(comm_F) :: self
    class(comm_F), pointer :: next
    next => self%nextLink
  end function next

  function prev(self)
    class(comm_F) :: self
    class(comm_F), pointer :: prev
    prev => self%prevLink
  end function prev
  
  subroutine setNext(self,next)
    class(comm_F) :: self
    class(comm_F), pointer :: next
    self%nextLink => next
  end subroutine setNext

  subroutine add(self,link)
    class(comm_F)         :: self
    class(comm_F), target :: link

    class(comm_F), pointer :: c
    
    c => self%nextLink
    do while (associated(c%nextLink))
       c => c%nextLink
    end do
    link%prevLink => c
    c%nextLink    => link
  end subroutine add  
  
end module comm_F_mod
