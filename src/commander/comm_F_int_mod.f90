module comm_F_int_mod
  use comm_utils
  use comm_bp_mod
  implicit none

  private
  public comm_F_int, F_int_ptr
  
  type, abstract :: comm_F_int
     class(comm_bp), pointer :: bp
   contains
     ! Data procedures
     procedure(evalIntF),     deferred :: eval
     procedure(updateIntF),   deferred :: update
  end type comm_F_int

  abstract interface
     ! Evaluate SED in brightness temperature normalized to reference frequency
     function evalIntF(self, theta)
       import comm_F_int, dp
       class(comm_F_int),                intent(in) :: self
       real(dp),          dimension(1:), intent(in) :: theta
       real(dp)                                     :: evalIntF
     end function evalIntF

     ! Compute/update integration look-up tables
     subroutine updateIntF(self, f_precomp, pol)
       import comm_F_int, dp, i4b
       class(comm_F_int),    intent(inout)           :: self
       real(dp),             intent(in),    optional :: f_precomp
       integer(i4b),         intent(in),    optional :: pol
     end subroutine updateIntF
  end interface

  type F_int_ptr
     class(comm_F_int), pointer :: p
  end type F_int_ptr
  
end module comm_F_int_mod
