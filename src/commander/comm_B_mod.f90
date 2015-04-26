module comm_B_mod
  use comm_param_mod
  use comm_map_mod
  implicit none

  private
  public comm_B

  type, abstract :: comm_B
     ! Data variables
     character(len=512)           :: type
     class(comm_mapinfo), pointer :: info
     real(dp), allocatable, dimension(:,:) :: b_l
   contains
     ! Data procedures
     procedure(matmulB),  deferred :: B
     procedure(matmulBt), deferred :: B_t
  end type comm_B

  abstract interface

     subroutine matmulB(self, map, res)
       import comm_map, comm_B, dp, i4b
       implicit none
       class(comm_B),   intent(in)    :: self
       class(comm_map), intent(in)    :: map
       class(comm_map), intent(inout) :: res
     end subroutine matmulB

     subroutine matmulBt(self, map, res)
       import comm_map, comm_B, dp, i4b
       implicit none
       class(comm_B),   intent(in)    :: self
       class(comm_map), intent(in)    :: map
       class(comm_map), intent(inout) :: res
     end subroutine matmulBt
  end interface

end module comm_B_mod
