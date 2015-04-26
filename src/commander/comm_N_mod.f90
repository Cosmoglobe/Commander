module comm_N_mod
  use comm_param_mod
  use comm_map_mod
  implicit none

  private
  public comm_N
  
  type, abstract :: comm_N
     ! Data variables
     character(len=512) :: type
     integer(i4b)       :: nside, nmaps, np
     logical(lgt)       :: pol
   contains
     ! Data procedures
     procedure(matmulInvN),     deferred :: invN
     procedure(matmulSqrtInvN), deferred :: sqrtInvN
     procedure(returnRMS),      deferred :: rms
  end type comm_N

  abstract interface
     ! Return map_out = invN * map
     subroutine matmulInvN(self, map, res)
       import comm_map, comm_N, dp, i4b
       implicit none
       class(comm_N),   intent(in)    :: self
       class(comm_map), intent(in)    :: map
       class(comm_map), intent(inout) :: res
     end subroutine matmulInvN

     ! Return map_out = sqrtInvN * map
     subroutine matmulSqrtInvN(self, map, res)
       import comm_map, comm_N, dp, i4b
       implicit none
       class(comm_N),   intent(in)    :: self
       class(comm_map), intent(in)    :: map
       class(comm_map), intent(inout) :: res
     end subroutine matmulSqrtInvN

     ! Return rms map
     subroutine returnRMS(self, res)
       import comm_map, comm_N, dp
       implicit none
       class(comm_N),   intent(in)    :: self
       class(comm_map), intent(inout) :: res
     end subroutine returnRMS
  end interface

end module comm_N_mod
