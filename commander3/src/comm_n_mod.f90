module comm_N_mod
  use comm_map_mod
  implicit none

  type :: comm_N
     ! Data variables
     class(comm_map),     pointer :: invN_diag => null()
     class(comm_mapinfo), pointer :: info      => null()
     class(comm_map),     pointer :: rms0      => null()
   contains
     procedure :: update_N
  end type comm_N

  type comm_N_ptr
     class(comm_N), pointer :: p => null()
  end type comm_N_ptr

  interface comm_N
     procedure constructor
  end interface comm_N

contains

  function constructor(info)
    implicit none
    class(comm_N),                      pointer       :: constructor
    type(comm_mapinfo), target,         intent(in)    :: info

    allocate(constructor)
    call constructor%update_N(info)

  end function constructor

  subroutine update_N(self, info)
    implicit none
    class(comm_N),                       intent(inout)          :: self
    class(comm_mapinfo),                 intent(in)             :: info

    self%rms0     => comm_map(info)
    call compute_invN_lm(self%rms0)

  end subroutine update_N

  subroutine compute_invN_lm(invN_diag)
    implicit none

    class(comm_map),  intent(inout) :: invN_diag

    call invN_diag%YtW_scalar
    
  end subroutine compute_invN_lm
  
end module comm_N_mod
