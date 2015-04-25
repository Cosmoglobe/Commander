module comm_N_mod
  use comm_param_mod
  implicit none

  private
  public comm_N
  
  type, abstract :: comm_N
     ! Data variables
     character(len=512) :: type
     integer(i4b)       :: nside, nmaps, np
     logical(lgt)       :: pol
     real(dp), allocatable, dimension(:,:) :: rms
   contains
     ! Data procedures
     procedure(matmulInvN),     deferred :: invN
     procedure(matmulSqrtInvN), deferred :: sqrtInvN
     procedure                           :: getRMS
  end type comm_N

  abstract interface
     ! Return map_out = invN * map
     function matmulInvN(self, map)
       import comm_N, dp, i4b
       class(comm_N),                                intent(in) :: self
       real(dp),      dimension(self%np,self%nmaps), intent(in) :: map
       real(dp),      dimension(self%np,self%nmaps)             :: matmulInvN
     end function matmulInvN

     ! Return map_out = sqrtInvN * map
     function matmulSqrtInvN(self, map)
       import comm_N, dp, i4b
       class(comm_N),                                intent(in) :: self
       real(dp),      dimension(self%np,self%nmaps), intent(in) :: map
       real(dp),      dimension(self%np,self%nmaps)             :: matmulSqrtInvN
     end function matmulSqrtInvN
  end interface

contains

  ! Return rms map
  function getRMS(self)
    class(comm_N),                               intent(in) :: self
    real(dp),      dimension(self%np,self%nmaps)            :: getRMS
    getRMS = self%rms
  end function getRMS
  
end module comm_N_mod
