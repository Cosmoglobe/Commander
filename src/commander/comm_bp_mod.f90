module comm_bp_mod
  use comm_param_mod
  use comm_bp_utils
  implicit none 

  private
  public comm_bp, comm_bp_ptr, initialize_bp_mod

  type :: comm_bp
     ! Data variables
     character(len=512) :: type, model
     integer(i4b)       :: n, npar
     real(dp)           :: threshold
     real(dp)           :: nu_c, a2t, f2t, a2sz
     real(dp), allocatable, dimension(:) :: nu0, nu, tau0, tau, delta
   contains
     ! Data procedures
     procedure     :: update_tau
     procedure     :: SED2F
     procedure     :: lineAmp_RJ
  end type comm_bp

  type comm_bp_ptr
     class(comm_bp), pointer :: p
  end type comm_bp_ptr


  interface comm_bp
     procedure constructor
  end interface comm_bp

  real(dp) :: ind_iras = 1.d0
  
contains

  subroutine initialize_bp_mod(cpar)
    implicit none

    type(comm_params), intent(in) :: cpar

    ! Set up global parameters
    T_CMB = cpar%T_CMB
    if (trim(cpar%MJysr_convention) == 'PSM') then
       ind_iras = 0.d0
    else if (trim(cpar%MJysr_convention) == 'IRAS') then
       ind_iras = 1.d0
    else
       write(*,*) 'Unsupported MJy/sr convention = ', trim(cpar%MJysr_convention)
       stop
    end if
    
  end subroutine initialize_bp_mod
  
  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs, detlabel)
    implicit none
    type(comm_params),           intent(in) :: cpar
    integer(i4b),                intent(in) :: id, id_abs
    character(len=*),            intent(in) :: detlabel
    class(comm_bp),     pointer             :: constructor

    integer(i4b)       :: i
    character(len=512) :: dir, label

    label = cpar%ds_label(id_abs)
    
    ! General parameters
    allocate(constructor)
    dir = trim(cpar%datadir) // '/'

    constructor%nu_c = cpar%ds_nu_c(id_abs)
    
    ! Define special case parameters
    constructor%type = cpar%ds_bptype(id_abs)
    select case (trim(constructor%type))
    case ('delta')
       constructor%threshold = 0.d0
    case ('LFI') 
       constructor%threshold = 0.d0
    case ('WMAP') 
       constructor%threshold = 0.d0
    case ('DIRBE') 
       constructor%threshold = 0.d0
    case ('HFI_cmb') 
       constructor%threshold = 1.d-7
    case ('PSM_LFI') 
       constructor%threshold = 1.d-7
    case ('HFI_submm') 
       constructor%threshold = 1.d-5
    case ('dame') 
       constructor%threshold = 0.d0
       
    case default
       call report_error('Error -- unsupported bandpass type = '//trim(constructor%type))
    end select

    ! Initialize raw bandpass
    if (trim(constructor%type) == 'delta') then
       allocate(constructor%nu0(1),constructor%tau0(1), constructor%nu(1), constructor%tau(1))
       constructor%n       = 1
       constructor%nu0(1)  = constructor%nu_c
       constructor%tau0(1) = 1.d0
    else
       call read_bandpass(trim(dir)//cpar%ds_bpfile(id_abs), detlabel, &
            & constructor%threshold, &
            & constructor%n, constructor%nu0, constructor%tau0)
       allocate(constructor%nu(constructor%n), constructor%tau(constructor%n))
    end if

    ! Initialize fitting model
    constructor%model = cpar%ds_bpmodel(id_abs)
    if (trim(constructor%model) == 'additive_shift') then
       constructor%npar = 1
       allocate(constructor%delta(constructor%npar))
       constructor%delta = 0.d0
    else if (trim(constructor%model) == 'powlaw_tilt') then
       constructor%npar = 1
       allocate(constructor%delta(constructor%npar))
       constructor%delta = 0.d0
    else
       call report_error('Error -- unsupported bandpass model = ' // trim(constructor%model))
    end if

    ! Read default delta from instrument parameter file
    call read_instrument_file(trim(cpar%datadir)//'/'//trim(cpar%cs_inst_parfile), &
         & 'delta', cpar%ds_label(id_abs), 0.d0, constructor%delta(1))

    ! Initialize active bandpass 
    call constructor%update_tau(constructor%delta)

  end function constructor
  

  
  subroutine update_tau(self, delta)
    implicit none

    class(comm_bp),                       intent(inout) :: self
    real(dp),       dimension(self%npar), intent(in)    :: delta

    integer(i4b) :: i, n
    real(dp), allocatable, dimension(:)  :: a, bnu_prime, bnu_prime_RJ, sz

    self%delta = delta
    
    n = self%n
    select case (trim(self%model))
    case ('powlaw_tilt')
       
       ! Power-law model, centered on nu_c
       self%nu = self%nu0
       do i = 1, n
          self%tau(i) = self%tau0(i) * (self%nu(i)/self%nu_c)**delta(1)
       end do
       
    case ('additive_shift') 
       
       ! Additive frequency shift
       self%tau = self%tau0
       do i = 1, n
          self%nu(i) = self%nu0(i) + 1d9*delta(1)
          if (self%nu(i) <= 0.d0) self%tau(i) = 0.d0
       end do
       
    end select

    ! Compute unit conversion factors
    allocate(a(n), bnu_prime(n), bnu_prime_RJ(n), sz(n))
    do i = 1, n
       a(i)            = comp_a2t(self%nu(i))          
       bnu_prime(i)    = comp_bnu_prime(self%nu(i))
       bnu_prime_RJ(i) = comp_bnu_prime_RJ(self%nu(i))
       sz(i)           = comp_sz_thermo(self%nu(i))
    end do

    select case (trim(self%type))
    case ('delta')
       
       self%a2t  = a(1)
       self%a2sz = 2.d0*self%nu_c**2*k_b/c**2 / &
            & (bnu_prime(1) * sz(1)) * 1.d-6
       self%f2t  = 1.d0 / bnu_prime(1) * 1.d-14
       
    case ('WMAP')
          
       ! See Appendix E of Bennett et al. (2013) for details
       self%a2t     = sum(self%tau) / sum(self%tau/a)
       self%a2sz    = sum(self%tau) / sum(self%tau/a * sz) * 1.d-6
       self%f2t     = sum(self%tau/self%nu**2 * (self%nu_c/self%nu)**ind_iras) * &
                          & 1.d-14 / sum(self%tau/self%nu**2 * bnu_prime)
       self%tau     = self%tau * a

    case ('LFI') 

       ! Multiply with normalization factor
       self%a2t     = tsum(self%nu, self%tau/self%nu**2 * bnu_prime_RJ) / &
                       & tsum(self%nu, self%tau/self%nu**2 * bnu_prime)
       self%a2sz    = tsum(self%nu, self%tau/self%nu**2 * bnu_prime_RJ) / &
                       & tsum(self%nu, self%tau/self%nu**2 * bnu_prime * sz) * 1.d-6
       self%f2t     = tsum(self%nu, self%tau/self%nu**2 * (self%nu_c/self%nu)**ind_iras) &
                       & * 1.d-14 / tsum(self%nu, self%tau/self%nu**2 * bnu_prime)
       self%tau     = self%tau / tsum(self%nu, self%tau/a)

    case ('HFI_cmb', 'PSM_LFI') 

       self%a2t     = tsum(self%nu, self%tau * bnu_prime_RJ) / tsum(self%nu, self%tau*bnu_prime)
       self%a2sz    = tsum(self%nu, self%tau * bnu_prime_RJ) / &
                       & tsum(self%nu, self%tau*bnu_prime*sz) * 1.d-6
       self%f2t     = tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * &
                       & 1.d-14 / tsum(self%nu, self%tau*bnu_prime)
       self%tau     = self%tau / tsum(self%nu, self%tau*bnu_prime)
       
    case ('HFI_submm') 

       self%a2t     = tsum(self%nu, self%tau * bnu_prime_RJ) / tsum(self%nu, self%tau*bnu_prime)
       self%a2sz    = tsum(self%nu, self%tau * bnu_prime_RJ) / &
                       & tsum(self%nu, self%tau*bnu_prime*sz) * 1.d-6
       self%f2t     = tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * &
                       & 1.d-14 / tsum(self%nu, self%tau*bnu_prime)
       self%tau     = self%tau / tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * 1.d14

    case ('DIRBE') 

       self%a2t     = tsum(self%nu, self%tau * bnu_prime_RJ) / tsum(self%nu, self%tau*bnu_prime)
       self%a2sz    = tsum(self%nu, self%tau * bnu_prime_RJ) / &
                       & tsum(self%nu, self%tau*bnu_prime*sz) * 1.d-6
       self%f2t     = tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * &
                       & 1.d-14 / tsum(self%nu, self%tau*bnu_prime)
       self%tau     = self%tau / tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * 1.d14

    ! NEW !
    case ('dame')
       self%a2t     = -1.d30 !tsum(self%nu, self%tau * bnu_prime_RJ) / tsum(self%nu, self%tau*bnu_prime)
       self%a2sz    = 1.0 !tsum(self%nu, self%tau * bnu_prime_RJ) / &
                       !& tsum(self%nu, self%tau*bnu_prime*sz) * 1.d-6
       self%f2t     = 1.0 !tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * &
                       !& 1.d-14 / tsum(self%nu, self%tau*bnu_prime)
       self%tau     = 1.0 !self%tau / tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * 1.d14

    end select
    deallocate(a, bnu_prime, bnu_prime_RJ, sz)

  end subroutine update_tau


  function SED2F(self, f)
    implicit none

    class(comm_bp),               intent(in) :: self
    real(dp),       dimension(:), intent(in) :: f
    real(dp)                                 :: SED2F

    select case (trim(self%type))
    case ('delta')
       SED2F = f(1) * self%a2t
    case ('LFI')
       SED2F = tsum(self%nu, self%tau * f)
    case ('HFI_cmb')
       SED2F = tsum(self%nu, self%tau * 2.d0*k_B*self%nu**2/c**2 * f)
    case ('PSM_LFI')
       SED2F = tsum(self%nu, self%tau * 2.d0*k_B*self%nu**2/c**2 * f)
    case ('HFI_submm') 
       SED2F = tsum(self%nu, self%tau * 2.d0*k_B*self%nu**2/c**2 * f)
    case ('DIRBE') 
       SED2F = tsum(self%nu, self%tau * 2.d0*k_B*self%nu**2/c**2 * f)
    case ('WMAP')
       SED2F = sum(self%tau * f)
    case ('dame') ! NEW
       SED2F = f(1) * self%a2t
    end select

  end function SED2F

  function lineAmp_RJ(self, nu)
    implicit none

    class(comm_bp), intent(in) :: self
    real(dp),       intent(in) :: nu
    real(dp)                   :: lineAmp_RJ

    integer(i4b) :: i
    real(dp)     :: x, tau
    real(dp), allocatable, dimension(:) :: bnu_prime_RJ, a2t

    if (nu < self%nu(1) .or. nu > self%nu(self%n)) then
       lineAmp_RJ = 0.d0
       return
    end if

    ! Read off correct bandpass coefficient; use linear interpolation
    i = 1
    do while (self%nu(i) < nu .and. i < self%n)
       i = i+1
    end do
    x   = (self%nu(i)-nu) / (self%nu(i)-self%nu(i-1))
    tau = x * self%tau(i-1) + (1.d0-x)*self%tau(i)

    select case (trim(self%type))
    case ('delta')

       if (nu /= self%nu_c) then
          lineAmp_RJ = 0.d0
       else
          lineAmp_RJ = nu/c 
       end if

    case ('WMAP') 

       ! See Appendix E of Bennett et al. (2013) for details
       lineAmp_RJ = tau/(self%nu(2)-self%nu(1)) * nu/c / sum(self%tau)

    case ('LFI') 

       ! Multiply with normalization factor
       allocate(bnu_prime_RJ(self%n))       
       bnu_prime_RJ = comp_bnu_prime_RJ(self%nu)
       lineAmp_RJ = tau/nu**2 * nu/c * compute_bnu_prime_RJ_single(nu) / &
            & tsum(self%nu, self%tau/self%nu**2 * bnu_prime_RJ)
       deallocate(bnu_prime_RJ)

    case ('HFI_cmb', 'PSM_LFI') 
          
       allocate(bnu_prime_RJ(self%n))
       bnu_prime_RJ = comp_bnu_prime_RJ(self%nu)
       lineAmp_RJ = tau * nu/c * compute_bnu_prime_RJ_single(nu) / &
            & tsum(self%nu, self%tau*bnu_prime_RJ)
       deallocate(bnu_prime_RJ)

    case ('HFI_submm') 
       
       lineAmp_RJ = tau * nu/c * compute_bnu_prime_RJ_single(nu) / &
            & (tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * 1.d-14)

    ! NEW
    case ('dame') 

       if (nu /= self%nu_c) then
          lineAmp_RJ = 0.d0
       else
          lineAmp_RJ = nu/c 
       end if

    end select

    lineAmp_RJ = lineAmp_RJ * 1.d9 ! Convert to uK_ant / (K_ant km/s)

  end function lineAmp_RJ

end module comm_bp_mod
