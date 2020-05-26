module comm_camb_mod
  use healpix_types
  implicit none

  integer(i4b), parameter :: CAMB_NPAR      = 7
  integer(i4b), parameter :: CAMB_H         = 1
  integer(i4b), parameter :: CAMB_TAU       = 2
  integer(i4b), parameter :: CAMB_OMEGA_B   = 3
  integer(i4b), parameter :: CAMB_OMEGA_C   = 4
  integer(i4b), parameter :: CAMB_AS        = 5
  integer(i4b), parameter :: CAMB_NS        = 6
  integer(i4b), parameter :: CAMB_R         = 7

  integer(i4b), parameter :: CAMB_NSPEC     = 6
  integer(i4b), parameter :: CAMB_TT        = 1
  integer(i4b), parameter :: CAMB_TE        = 2
  integer(i4b), parameter :: CAMB_TB        = 3
  integer(i4b), parameter :: CAMB_EE        = 4
  integer(i4b), parameter :: CAMB_EB        = 5
  integer(i4b), parameter :: CAMB_BB        = 6
  
  type :: camb
     logical(lgt) :: do_lensing
     integer(i4b) :: lmax
     real(dp), dimension(CAMB_NPAR) :: theta
   contains
     procedure    :: getCls
     procedure    :: setParam
  end type camb

  interface camb
     procedure constructor
  end interface camb
  
contains

  !**************************************************
  !             Constructor routine
  !**************************************************
  function constructor(lmax, do_lensing)
    implicit none
    integer(i4b),         intent(in) :: lmax
    logical(lgt),         intent(in) :: do_lensing
    class(camb), pointer             :: constructor

    ! General parameters
    allocate(constructor)
    constructor%do_lensing          = do_lensing
    constructor%lmax                = lmax

    ! Initialize cosmological parameters to Planck 2015
    constructor%theta(CAMB_H)       = 0.d0
    constructor%theta(CAMB_TAU)     = 0.d0
    constructor%theta(CAMB_OMEGA_B) = 0.d0
    constructor%theta(CAMB_OMEGA_C) = 0.d0
    constructor%theta(CAMB_AS)      = 0.d0
    constructor%theta(CAMB_NS)      = 0.d0
    constructor%theta(CAMB_R)       = 0.d0
    
  end function constructor

  !**************************************************
  !             Main CAMB interface
  !**************************************************
  subroutine getCls(self, cls)
    implicit none
    class(camb),                                intent(in)  :: self
    real(dp),      allocatable, dimension(:,:), intent(out) :: cls

    ! Return power spectrum; filling in this array is Jeff's task
    if (.not. allocated(cls)) allocate(cls(0:self%lmax,CAMB_NSPEC))
    cls = 0.d0
    
  end subroutine getCls

  subroutine setParam(self, h, tau, Omega_b, Omega_c, A_s, n_s, r)
    implicit none
    class(camb),  intent(inout)        :: self
    real(dp),     intent(in), optional :: h, tau, Omega_b, Omega_c, A_s, n_s, r

    if (present(h))       self%theta(CAMB_H)       = h
    if (present(tau))     self%theta(CAMB_TAU)     = tau
    if (present(Omega_b)) self%theta(CAMB_OMEGA_B) = Omega_b
    if (present(Omega_c)) self%theta(CAMB_OMEGA_C) = Omega_c
    if (present(A_s))     self%theta(CAMB_AS)      = A_s
    if (present(n_s))     self%theta(CAMB_NS)      = n_s
    if (present(r))       self%theta(CAMB_R)       = r
    
  end subroutine setParam
  
end module comm_camb_mod
