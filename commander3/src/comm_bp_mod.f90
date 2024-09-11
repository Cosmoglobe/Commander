!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
module comm_bp_mod
  use comm_param_mod
  use comm_bp_utils
  implicit none 

  type :: comm_bp
     ! Data variables
     character(len=512) :: type, model
     logical(lgt)       :: sample_bandpass
     integer(i4b)       :: n, npar
     real(dp)           :: threshold
     real(dp)           :: nu_c, a2t, f2t, a2sz, unit_scale, nu_eff, a2f
     real(dp), allocatable, dimension(:) :: nu0, nu, tau0, tau, delta, a2f_arr
   contains
     ! Data procedures
     procedure     :: update_tau
     procedure     :: SED2F
     procedure     :: lineAmp_RJ
  end type comm_bp

  type comm_bp_ptr
     class(comm_bp), pointer :: p => null()
  end type comm_bp_ptr


  interface comm_bp
     procedure constructor_bp
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
  function constructor_bp(cpar, id, id_abs, detlabel, subdets) result(c)
    !
    ! Initialization routine (constructor) for bandpass objects. Reads in bandpass 
    ! data, and precomputes default unit conversions etc.
    !
    ! Arguments:
    ! ----------
    ! cpar: type(comm_params)
    !    Commander parameter structure
    ! id: int (scalar)
    !    Frequency channel ID/counter, only counting active bands
    ! id_abs: int (scalar)
    !    Frequency channel ID/counter, counting all bands (as defined in parameter file)
    ! detlabel: string (scalar; optional)
    !    Detector label; typically used for full-frequency bands 
    ! subdets: string (scalar; optional)
    !    Comma-separated string with sub-detector labels. Used for TOD-type bands
    !
    ! Returns:
    ! --------
    ! constructor: class(comm_bp)
    !    Pointer to new comm_bp object
    ! 
    implicit none
    type(comm_params),                intent(in)           :: cpar
    integer(i4b),                     intent(in)           :: id, id_abs
    character(len=*),                 intent(in), optional :: detlabel, subdets
    class(comm_bp),     pointer                            :: c

    integer(i4b)       :: i, j, ndet
    character(len=512) :: label
    character(len=25)  :: dets(1500)
    real(dp), allocatable, dimension(:) :: nu0, tau0
    
    label = cpar%ds_label(id_abs)
    
    ! General parameters
    allocate(c)
    
    c%nu_c = cpar%ds_nu_c(id_abs)
    ! Define special case parameters
    c%type = cpar%ds_bptype(id_abs)
    select case (trim(c%type))
    case ('delta')
       c%threshold = 0.d0
    case ('LFI') 
       c%threshold = 0.d0
    case ('WMAP') 
       c%threshold = 0.d0
    case ('DIRBE') 
       c%threshold = 0.d0
    case ('AKARI') 
       c%threshold = 0.d0
    case ('HFI_cmb') 
       c%threshold = 1.d-7
    case ('PSM_LFI') 
       c%threshold = 1.d-7
    case ('HFI_submm') 
       c%threshold = 1.d-5
    case ('dame') 
       c%threshold = 0.d0
    case ('LB')
       c%threshold = 0.d0
    case ('SPIDER')
       c%threshold = 0.d0
    case default
       call report_error('Error -- unsupported bandpass type = '//trim(c%type))
    end select

    ! Initialize unit scale
    if (trim(cpar%ds_unit(id_abs)) == 'mK_cmb') then
       c%unit_scale = 1.d-3
    else if (trim(cpar%ds_unit(id_abs)) == 'K_cmb') then
       c%unit_scale = 1.d-6
    else
       c%unit_scale = 1.d0
    end if


    ! Initialize raw bandpass
    if (trim(c%type) == 'delta') then
       allocate(c%nu0(1),c%tau0(1), c%nu(1), c%tau(1))
       c%n       = 1
       c%nu0(1)  = c%nu_c
       c%tau0(1) = 1.d0
    else
       if (present(detlabel)) then
          call read_bandpass(trim(cpar%ds_bpfile(id_abs)), detlabel, &
               & c%threshold, &
               & c%n, c%nu0, c%tau0)
       else 
          if (index(subdets, '.txt') /=0) then
               ndet = count_detectors(trim(subdets))
               call get_detectors(trim(subdets), dets, ndet)
          else
               ndet = num_tokens(subdets, ",")
               call get_tokens(subdets, ",", dets, ndet)
          end if
          if (c%threshold == 0.d0) then
               call read_bandpass(trim(cpar%ds_bpfile(id_abs)), dets(1), &
                    & c%threshold, &
                    & c%n, c%nu0, c%tau0)
               do i = 2, ndet
                    call read_bandpass(trim(cpar%ds_bpfile(id_abs)), dets(i), &
                        & c%threshold, c%n, nu0, tau0)
                    c%tau0 = c%tau0 + tau0
                    deallocate(nu0, tau0)
               end do
               c%tau0 = c%tau0 / ndet
          else
               print *, "got to nonzero threshold, aborting"
               stop
               call read_bandpass_nonzero_threshold(cpar%ds_bpfile(id_abs), dets, ndet, &
                    & c%threshold, &
                    & c%n, c%nu0, c%tau0)
          end if
       end if
       allocate(c%nu(c%n), c%tau(c%n))
       if (trim(c%type) == 'DIRBE' .or. trim(c%type) == 'AKARI') then
          allocate(c%a2f_arr(c%n))
       end if
    end if
    ! Initialize fitting model
    c%model = cpar%ds_bpmodel(id_abs)
    if (trim(c%model) == 'additive_shift') then
       c%npar = 1
       allocate(c%delta(c%npar))
       c%delta = 0.d0
    else if (trim(c%model) == 'powlaw_tilt') then
       c%npar = 1
       allocate(c%delta(c%npar))
       c%delta = 0.d0
    else
       call report_error('Error -- unsupported bandpass model = ' // trim(c%model))
    end if

    ! Read default delta from instrument parameter file
    call read_instrument_file(trim(cpar%cs_inst_parfile), &
         & 'delta', cpar%ds_label(id_abs), 0.d0, c%delta(1))

    ! Initialize active bandpass 
    call c%update_tau(c%delta)

    ! WARNING! Should be replaced with proper integral. See planck2013 HFI spectral response eq. 2
    c%nu_eff = sum(c%tau*c%nu)/sum(c%tau)
    
  end function constructor_bp
  

  
  subroutine update_tau(self, delta)
    implicit none

    class(comm_bp),                       intent(inout) :: self
    real(dp),       dimension(self%npar), intent(in)    :: delta
    real(dp), allocatable, dimension(:)  :: a, bnu_prime, bnu_prime_RJ, sz

    integer(i4b) :: i, n

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
          !if (abs(self%nu(i))>1e15) write(*,*) "i, nu, nu0, delta: ", i, self%nu(i), self%nu0(i), 1d9*delta(1)
       end do
    case default
       call report_error('Error -- unsupported bandpass model = ' // trim(self%model))
    end select

    ! Compute unit conversion factors
    allocate(a(n), bnu_prime(n), bnu_prime_RJ(n), sz(n))
    do i = 1, n
       if (trim(self%type) == 'DIRBE' .or. trim(self%type) == 'AKARI') then
          bnu_prime_RJ(i) = comp_bnu_prime_RJ(self%nu(i))
          ! These overflow in exp(x) due to large x
          bnu_prime(i)    = 1.d0 !comp_bnu_prime(self%nu(i))
          sz(i)           = 1.d0 !comp_sz_thermo(self%nu(i))
       else if (trim(self%type) == 'HFI_submm') then
          bnu_prime(i)    = comp_bnu_prime(self%nu(i))
          bnu_prime_RJ(i) = comp_bnu_prime_RJ(self%nu(i))
          sz(i)           = comp_sz_thermo(self%nu(i))
       else
          a(i)            = comp_a2t(self%nu(i))          
          bnu_prime(i)    = comp_bnu_prime(self%nu(i))
          bnu_prime_RJ(i) = comp_bnu_prime_RJ(self%nu(i))
          sz(i)           = comp_sz_thermo(self%nu(i))
       end if
    end do

    select case (trim(self%type))
    case ('delta')
       
       self%a2t  = a(1)
       self%a2sz = 2.d0*self%nu_c**2*k_b/c**2 / &
            & (bnu_prime(1) * sz(1)) * 1.d-6
       self%f2t  = 1.d0 / bnu_prime(1) * 1.d-14
       self%a2f  = 2.d0*self%nu_c**2*k_b/c**2 * 1d14
       
    case ('WMAP')

       !write(*,*) self%nu
          
       ! See Appendix E of Bennett et al. (2013) for details
       self%tau     = self%tau / sum(self%tau)
       self%a2t     = sum(self%tau) / sum(self%tau/a)
       self%a2sz    = sum(self%tau) / sum(self%tau/a * sz) * 1.d-6
       self%f2t     = sum(self%tau/self%nu**2 * (self%nu_c/self%nu)**ind_iras) * &
                          & 1.d-14 / sum(self%tau/self%nu**2 * bnu_prime)
       self%a2f     = sum(self%tau) / sum(self%tau/self%nu**2 * (self%nu_c / self%nu)**ind_iras) * 1d14
       self%tau     = self%tau * a

    case ('LFI') 

       ! Multiply with normalization factor
       self%a2t     = tsum(self%nu, self%tau/self%nu**2 * bnu_prime_RJ) / &
                       & tsum(self%nu, self%tau/self%nu**2 * bnu_prime)
       self%a2sz    = tsum(self%nu, self%tau/self%nu**2 * bnu_prime_RJ) / &
                       & tsum(self%nu, self%tau/self%nu**2 * bnu_prime * sz) * 1.d-6
       self%f2t     = tsum(self%nu, self%tau/self%nu**2 * (self%nu_c/self%nu)**ind_iras) &
                       & * 1.d-14 / tsum(self%nu, self%tau/self%nu**2 * bnu_prime)
       self%a2f     = tsum(self%nu, self%tau/self%nu**2 * bnu_prime_RJ) / tsum(self%nu, self%tau/self%nu**2 * (self%nu_c / self%nu)**ind_iras) * 1d14
       self%tau     = self%tau / tsum(self%nu, self%tau/a)

    case ('HFI_cmb', 'PSM_LFI', 'SPIDER') 

       self%a2t     = tsum(self%nu, self%tau * bnu_prime_RJ) / tsum(self%nu, self%tau*bnu_prime)
       self%a2sz    = tsum(self%nu, self%tau * bnu_prime_RJ) / &
                       & tsum(self%nu, self%tau*bnu_prime*sz) * 1.d-6
       self%f2t     = tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * &
                       & 1.d-14 / tsum(self%nu, self%tau*bnu_prime)
       self%a2f     = tsum(self%nu, self%tau * bnu_prime_RJ) / tsum(self%nu, self%tau * (self%nu_c / self%nu)**ind_iras) * 1d14
       self%tau     = self%tau / tsum(self%nu, self%tau*bnu_prime)
       
    case ('HFI_submm') 

       self%a2t     = tsum(self%nu, self%tau * bnu_prime_RJ) / tsum(self%nu, self%tau*bnu_prime)
       self%a2sz    = tsum(self%nu, self%tau * bnu_prime_RJ) / &
                       & tsum(self%nu, self%tau*bnu_prime*sz) * 1.d-6
       self%f2t     = tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * &
                       & 1.d-14 / tsum(self%nu, self%tau*bnu_prime)
       self%a2f     = tsum(self%nu, self%tau * bnu_prime_RJ) / tsum(self%nu, self%tau * (self%nu_c / self%nu)**ind_iras) * 1d14
       self%tau     = self%tau / tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * 1.d14
 
    case ('DIRBE', 'AKARI') 
      ! a = brightness temperature (antenenna temperature) [K_RJ]
      ! t = thermodynamic temperature [K_CMB]
      ! f = flux intensity [MJy/sr]
      ! sz = ?

       self%a2t     = tsum(self%nu, self%tau * bnu_prime_RJ) / tsum(self%nu, self%tau*bnu_prime)
       self%a2sz    = tsum(self%nu, self%tau * bnu_prime_RJ) / &
                       & tsum(self%nu, self%tau*bnu_prime*sz) * 1.d-6
       self%f2t     = tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * &
                       & 1.d-14 / tsum(self%nu, self%tau*bnu_prime)
       self%a2f     = tsum(self%nu, self%tau * bnu_prime_RJ) / tsum(self%nu, self%tau * (self%nu_c / self%nu)) * 1d14
       self%a2f_arr = bnu_prime_RJ / (self%nu_c / self%nu)**ind_iras * 1d14
       self%tau     = self%tau / tsum(self%nu, self%tau)
       !self%tau     = self%tau / tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * 1.d14


    ! NEW !
    case ('dame')
       self%a2t     = -1.d30 !tsum(self%nu, self%tau * bnu_prime_RJ) / tsum(self%nu, self%tau*bnu_prime)
       self%a2sz    = 1.0 !tsum(self%nu, self%tau * bnu_prime_RJ) / &
                       !& tsum(self%nu, self%tau*bnu_prime*sz) * 1.d-6
       self%f2t     = 1.0 !tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * &
                       !& 1.d-14 / tsum(self%nu, self%tau*bnu_prime)
       self%tau     = 1.0 !self%tau / tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * 1.d14
       self%a2f     = 1.0
    case ('LB')
       
       self%a2t     = tsum(self%nu, self%tau/self%nu**2 * bnu_prime_RJ) / &
                       & tsum(self%nu, self%tau/self%nu**2 * bnu_prime)
       self%a2sz    = tsum(self%nu, self%tau/self%nu**2 * bnu_prime_RJ) / &
                       & tsum(self%nu, self%tau/self%nu**2 * bnu_prime * sz) * 1.d-6
       self%f2t     = tsum(self%nu, self%tau/self%nu**2 * (self%nu_c/self%nu)**ind_iras) &
                   & * 1.d-14 / tsum(self%nu, self%tau/self%nu**2 * bnu_prime)
       self%a2f     = tsum(self%nu, self%tau/self%nu**2 * bnu_prime_RJ) / tsum(self%nu, self%tau/self%nu**2 * (self%nu_c/self%nu)**ind_iras) * 1.d-14 
       self%tau     = self%tau / tsum(self%nu, self%tau/a)

   !  case ('SPIDER') 

   !     self%a2t     = tsum(self%nu, self%tau/self%nu**2 * bnu_prime_RJ) / &
   !                       & tsum(self%nu, self%tau/self%nu**2 * bnu_prime)
   !     self%a2sz    = tsum(self%nu, self%tau/self%nu**2 * bnu_prime_RJ) / &
   !                       & tsum(self%nu, self%tau/self%nu**2 * bnu_prime * sz) * 1.d-6
   !     self%f2t     = tsum(self%nu, self%tau/self%nu**2 * (self%nu_c/self%nu)**ind_iras) &
   !                       & * 1.d-14 / tsum(self%nu, self%tau/self%nu**2 * bnu_prime)
   !     self%tau     = self%tau / tsum(self%nu, self%tau/a)

    end select
    deallocate(a, bnu_prime, bnu_prime_RJ, sz)

  end subroutine update_tau

  function SED2F(self, f)
    implicit none

    ! Implementation of the mixing matrix (SED2F = M, f = frequency scaling of a component).
    ! Depending on the units of the bandpass a different function which converts from K_RJ sed to DATA units bandpass integrated.
    ! See BP1 footnote 7 for K_RJ -> MJy/sr (additional factor of 2.d0*k_B*self%nu**2/c**2).

    class(comm_bp),               intent(in) :: self
    real(dp),       dimension(:), intent(in) :: f
    real(dp)                                 :: SED2F

    integer(i4b) :: i, j
    real(dp)     :: a2f1, a2f2, a2fc, K, Inu0

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
    case ('DIRBE', 'AKARI') 
       !SED2F = tsum(self%nu, self%tau * 2.d0*k_B*self%nu**2/c**2 * f) !* 1d14
       i     = locate(self%nu, self%nu_c)
       Inu0  = (self%a2f_arr(i)*f(i) + (self%a2f_arr(i+1)*f(i+1)-self%a2f_arr(i)*f(i))*(self%nu_c-self%nu(i))/(self%nu(i+1)-self%nu(i)))
       !  if (Inu0 < TINY(1.0_dp)) then
       if (Inu0 == 0.d0) then
         SED2F = 0.d0
       else
         K     = tsum(self%nu, self%tau * self%a2f_arr*f/Inu0) / tsum(self%nu, self%tau * self%nu_c/self%nu)
         SED2F = K * Inu0
         !SED2F = K
         !  print *, SED2F, tsum(self%nu, self%tau * self%a2f_arr*f) / tsum(self%nu, self%tau * self%nu_c/self%nu)
      end if
    case ('WMAP')
       SED2F = sum(self%tau * f)
    case ('dame') ! NEW
       SED2F = f(1) * self%a2t
    case ('LB')
       SED2F = tsum(self%nu, self%tau * f)
    case ('SPIDER')
       SED2F = tsum(self%nu, self%tau * 2.d0*k_B*self%nu**2/c**2 * f)
    case default
       write(*,*) 'Unsupported bandpass type'
       stop
    end select
    SED2F = SED2F * self%unit_scale

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

    case ('HFI_cmb', 'HFI_submm', 'PSM_LFI', 'SPIDER', 'DIRBE', 'AKARI', 'FIRAS') 
          
       allocate(bnu_prime_RJ(self%n))
       bnu_prime_RJ = comp_bnu_prime_RJ(self%nu)
       lineAmp_RJ = tau * nu/c * compute_bnu_prime_RJ_single(nu) / &
            & tsum(self%nu, self%tau*bnu_prime_RJ)
       deallocate(bnu_prime_RJ)

!!$    case ('HFI_submm') 
!!$       
!!$       lineAmp_RJ = tau * nu/c * compute_bnu_prime_RJ_single(nu) / &
!!$            & (tsum(self%nu, self%tau * (self%nu_c/self%nu)**ind_iras) * 1.d-14)

    ! NEW
    case ('dame') 

       if (nu /= self%nu_c) then
          lineAmp_RJ = 0.d0
       else
          lineAmp_RJ = nu/c 
       end if

    case default
       write(*,*) 'Unsupported bandpass type'
       stop
    end select

    lineAmp_RJ = lineAmp_RJ * 1.d9 ! Convert to uK_ant / (K_ant km/s)

  end function lineAmp_RJ


end module comm_bp_mod
