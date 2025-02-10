module comm_N_rms_mod
  use comm_N_mod
  implicit none

  private
  public comm_N_rms, comm_N_rms_ptr
  
  type, extends (comm_N) :: comm_N_rms
     class(comm_map), pointer :: rms0       => null()
   contains
     procedure :: update_N    => update_N_rms
  end type comm_N_rms

  interface comm_N_rms
     procedure constructor
  end interface comm_N_rms

  type comm_N_rms_ptr
     type(comm_N_rms), pointer :: p => null()
  end type comm_N_rms_ptr
  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(info)
    implicit none
    class(comm_N_rms),                  pointer       :: constructor
    type(comm_mapinfo), target,         intent(in)    :: info

    allocate(constructor)
    call constructor%update_N(info)

  end function constructor

  subroutine update_N_rms(self, info)
    implicit none
    class(comm_N_rms),                   intent(inout)          :: self
    class(comm_mapinfo),                 intent(in)             :: info

    self%rms0     => comm_map(info)
    call compute_invN_lm(self%rms0)

  end subroutine update_N_rms

end module comm_N_rms_mod
