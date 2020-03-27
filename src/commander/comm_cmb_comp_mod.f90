module comm_cmb_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_F_int_0D_mod
  use comm_data_mod
  use comm_bp_utils
  implicit none

  private
  public comm_cmb_comp, apply_cmb_dipole_prior

  !**************************************************
  !           CMB component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_cmb_comp
     real(dp) :: cmb_dipole_prior(3)
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

    ! Prepare CMB dipole prior
    if (trim(cpar%cmb_dipole_prior_mask) /= 'none') then
       constructor%priormask        => comm_map(constructor%x%info, trim(cpar%datadir)//'/'//trim(cpar%cmb_dipole_prior_mask))
       constructor%cmb_dipole_prior =  constructor%cmb_dipole_prior / constructor%RJ2unit_(1)
    end if

  end function constructor

  ! Definition:
  !    SED  = conversion between thermodynamic and brightness temperature = 1/a2t
  function evalSED(self, nu, band, pol, theta)
    implicit none
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
    implicit none
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
             !if (.not. associated(self%F_int(k,i,j)%p)) then
                self%F_int(k,i,j)%p => comm_F_int_0D(self, data(i)%bp(j)%p, k, f_precomp=f)
             !else
             !   call self%F_int(k,i,j)%p%update(f_precomp=f)
             !end if
          end do
       end do
    end do

  end subroutine updateIntF

  subroutine apply_cmb_dipole_prior(cpar, handle)
    implicit none
    type(comm_params),   intent(in)    :: cpar
    type(planck_rng),    intent(inout) :: handle

    integer(i4b) :: i, l, m
    real(dp)     :: md(0:3)
    class(comm_map),  pointer :: map
    class(comm_comp), pointer :: c => null()

    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_cmb_comp)
          call c%x%Y
          md = c%x%fit_MDpoles(c%priormask)
          
          do i = 0, c%x%info%nalm-1
             call c%x%info%i2lm(i,l,m)
             if (l == 1 .and. m == -1) then   ! Y dipole
                c%x%alm(i,1) = c%x%alm(i,1) - sqrt(4.d0*pi/3.d0) * md(2) + cpar%cmb_dipole_prior(2)
             end if
             if (l == 1 .and. m ==  0) then   ! Z dipole
                c%x%alm(i,1) = c%x%alm(i,1) - sqrt(4.d0*pi/3.d0) * md(3) + cpar%cmb_dipole_prior(3)
             end if
             if (l == 1 .and. m ==  1) then   ! X dipole
                c%x%alm(i,1) = c%x%alm(i,1) + sqrt(4.d0*pi/3.d0) * md(1) + cpar%cmb_dipole_prior(1)
             end if
          end do

       end select
       c => c%next()
    end do

  end subroutine apply_cmb_dipole_prior


end module comm_cmb_comp_mod
