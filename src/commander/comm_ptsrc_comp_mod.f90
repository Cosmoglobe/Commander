module comm_ptsrc_comp_mod
  use comm_param_mod
  use comm_template_comp_mod
  implicit none

  private
  public comm_ptsrc_comp

  type Tnu
     real(dp), allocatable, dimension(:,:) :: map
  end type Tnu
  
  !**************************************************
  !            Template class
  !**************************************************
  type, extends (comm_template_comp) :: comm_ptsrc_comp
   contains
     procedure :: initPtsrc
  end type comm_ptsrc_comp

contains

  subroutine initPtsrc(self, cpar, id)
    implicit none
    class(comm_ptsrc_comp)                :: self
    type(comm_params),         intent(in) :: cpar
    integer(i4b),              intent(in) :: id

    call self%initTemplate(cpar, id)

    ! Initialize variables specific to template source type
    
  end subroutine initPtsrc

end module comm_ptsrc_comp_mod
