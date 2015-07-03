module comm_B_mod
  use comm_param_mod
  use comm_map_mod
  implicit none

  private
  public comm_B

  type, abstract :: comm_B
     ! Data variables
     character(len=512)           :: type
     integer(i4b)                 :: lmax
     class(comm_mapinfo), pointer :: info
     real(dp), allocatable, dimension(:,:) :: b_l
   contains
     ! Data procedures
     procedure(matmulB),  deferred :: conv
  end type comm_B

  abstract interface
     subroutine matmulB(self, alm_in, alm_out, trans, map)
       import comm_map, comm_B, dp, i4b, lgt
       implicit none
       class(comm_B),   intent(in)    :: self
       logical(lgt),    intent(in)    :: alm_in, alm_out
       logical(lgt),    intent(in)    :: trans
       class(comm_map), intent(inout) :: map
     end subroutine matmulB
  end interface

end module comm_B_mod
