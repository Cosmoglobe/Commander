module comm_cmb_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_F_int_0D_mod
  use comm_data_mod
  use comm_bp_utils
  implicit none

  private
  public comm_cmb_comp

  !**************************************************
  !           CMB component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_cmb_comp
   contains
     procedure :: S            => evalSED
     procedure :: update_F_int => updateIntF
  end type comm_cmb_comp

  interface comm_cmb_comp
     procedure constructor
  end interface comm_cmb_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_cmb_comp), pointer   :: constructor

    integer(i4b) :: i, j, k
    real(dp)     :: f
    
    ! General parameters
    allocate(constructor)
    call constructor%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    constructor%npar         = 0

    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(3,numband,0:constructor%ndet))
    call constructor%update_F_int
    
    ! Initialize mixing matrix
    call constructor%updateMixmat

  end function constructor

  ! Definition:
  !    SED  = conversion between thermodynamic and brightness temperature = 1/a2t
  function evalSED(self, nu, band, pol, theta)
    class(comm_cmb_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    real(dp) :: x
    x       = h*nu / (k_B*T_CMB)
    evalSED = (x**2 * exp(x)) / (exp(x)-1.d0)**2

  end function evalSED

  ! Update band integration lookup tables
  subroutine updateIntF(self, band)
    class(comm_cmb_comp),                    intent(inout)        :: self
    integer(i4b),                            intent(in), optional :: band

    integer(i4b) :: i, j, k
    real(dp)     :: f

    do k = 1, 3
       do i = 1, numband
          do j = 0, data(i)%ndet
             if (k > 1) then
                if (self%nu_ref(k) == self%nu_ref(k-1)) then
                   self%F_int(k,i,j)%p => self%F_int(k-1,i,j)%p
                   cycle
                end if
             end if
             f = comp_a2t(self%nu_ref(k)) / data(i)%bp(j)%p%a2t * data(i)%RJ2data(j)
             self%F_int(k,i,j)%p => comm_F_int_0D(self, data(i)%bp(j)%p, k, f_precomp=f)
          end do
       end do
    end do

  end subroutine updateIntF

end module comm_cmb_comp_mod
