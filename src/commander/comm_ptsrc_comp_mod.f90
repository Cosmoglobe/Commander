module comm_ptsrc_comp_mod
  use comm_param_mod
  use comm_comp_mod
  implicit none

  private
  public comm_ptsrc_comp

  !**************************************************
  !            Point-source class
  !**************************************************
  type, abstract, extends (comm_comp) :: comm_ptsrc_comp
     integer(i4b) :: nsrc
     real(dp), allocatable, dimension(:,:,:) :: T
   contains
     procedure :: initPtsrc
  end type comm_ptsrc_comp

contains

  subroutine initPtsrc(self, cpar, id)
    implicit none
    class(comm_ptsrc_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id

    call self%initComp(cpar, id)

    ! Initialize variables specific to diffuse source type

  end subroutine initPtsrc
  
end module comm_ptsrc_comp_mod
