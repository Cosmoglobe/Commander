module comm_zodi_mod
  use comm_utils
  use comm_param_mod
  implicit none

  private
  public :: initialize_zodi_mod, compute_zodi_template

  real(dp) :: myvariable


contains

  subroutine initialize_zodi_mod(cpar)
    implicit none
    type(comm_params), intent(in) :: cpar

    myvariable = 16.4d0
    
  end subroutine initialize_zodi_mod

  subroutine compute_zodi_template(nside, pix, nu, s_zodi)
    implicit none
    integer(i4b),                   intent(in)  :: nside
    integer(i4b), dimension(1:,1:), intent(in)  :: pix
    real(dp),     dimension(1:),    intent(in)  :: nu
    real(sp),     dimension(1:,1:), intent(out) :: s_zodi

    integer(i4b) :: ndet, ntod

    ntod = size(pix,1)
    ndet = size(pix,2)

    s_zodi = 0.

  end subroutine compute_zodi_template


end module comm_zodi_mod
