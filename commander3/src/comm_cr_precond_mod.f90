module comm_cr_precond_mod
  use comm_utils
  implicit none

  type invM
     integer(i4b)                              :: n
     integer(i4b), allocatable, dimension(:)   :: ind, comp2ind
     real(dp),     allocatable, dimension(:,:) :: M0, M
  end type invM

  type precond
     type(invM), allocatable, dimension(:,:)      :: invM_diff ! (0:nalm-1,nmaps)
     type(invM), allocatable, dimension(:,:)      :: invM_src  ! (1,nmaps)
     type(invM), allocatable, dimension(:,:)      :: invM_temp ! (1,1)
  end type precond

  type(precond) :: P_cr
  
end module comm_cr_precond_mod
