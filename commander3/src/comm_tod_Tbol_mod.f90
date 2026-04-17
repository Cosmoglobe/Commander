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

module comm_tod_Tbol_mod
  use comm_utils
  use comm_status_mod
  implicit none
  private

  public comm_Tbol, Tbol_ptr

  integer(i4b), parameter :: NSPLINE = 1000 ! Number of nodes between 0 and samprate/2
  
  type :: comm_Tbol
     integer(i4b) :: npole
     real(dp)     :: samprate
     real(dp)     :: tau_stray, S_phase
     real(dp), allocatable, dimension(:) :: a
     real(dp), allocatable, dimension(:) :: tau
     type(spline_type) :: TF_real, TF_imag
  contains
    procedure :: estimate_params
    procedure :: get_HFI_electronic_transfunc
    procedure :: get_HFI_Homega
    procedure :: update           => update_Tbol
    procedure :: convolve         => convolve_Tbol
    procedure :: dealloc          => dealloc_Tbol
  end type comm_Tbol   

  interface comm_Tbol
    procedure Tbol_constructor
  end interface comm_Tbol

  type Tbol_ptr
    class(comm_Tbol), pointer :: p => null()
  end type Tbol_ptr

contains
 
  !initializes a comm_Tbol class
  function Tbol_constructor(samprate, par_init) result(c)
    implicit none
    real(dp),                       intent(in) :: samprate
    real(dp),         dimension(:), intent(in) :: par_init
    class(comm_Tbol), pointer                  :: c

    integer(i4b) :: i
    
    allocate(c)
    c%npole    = count(par_init(3::2) > 0.)
    c%samprate = samprate
    
    allocate(c%a(c%npole), c%tau(c%npole))
    call c%update(par_init)
    
  end function Tbol_constructor

  subroutine dealloc_Tbol(self)
    implicit none
    class(comm_Tbol), intent(inout)      :: self
    deallocate(self%a)
    deallocate(self%tau)
  end subroutine dealloc_Tbol
  
  ! estimates Tbol parameters for each detector
  subroutine estimate_params(self)
    implicit none
    class(comm_Tbol), intent(inout)      :: self
  end subroutine estimate_params

  ! Routine for convolving with bolometer transfer function, tod = T * tod
  subroutine convolve_Tbol(self, tod, plan_fwd, plan_back)
    implicit none
    class(comm_Tbol),               intent(in)    :: self
    real(sp),         dimension(:), intent(inout) :: tod
    integer*8,                      intent(inout) :: plan_fwd, plan_back
    
    integer(i4b) :: i, j, k, l, n, ntod, nomp, nfft, err

    real(dp)     :: f
    complex(dpc) :: TF
    real(sp),     allocatable, dimension(:)   :: dt
    complex(spc), allocatable, dimension(:)   :: dv

    ntod     = size(tod)
    nomp     = 1
    nfft     = 2 * ntod
    n        = nfft / 2 + 1

    !call sfftw_init_threads(err)
    !call sfftw_plan_with_nthreads(nomp)

    allocate(dt(nfft), dv(0:n-1))
    !call sfftw_plan_dft_r2c_1d(plan_fwd,  nfft, dt, dv, fftw_estimate + fftw_unaligned)
    !call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)

    ! FFT
    dt(1:ntod)           = tod
    dt(2*ntod:ntod+1:-1) = tod(1:ntod)
    call timer%start(TOT_FFT)
    call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
    call timer%stop(TOT_FFT)

    ! Apply transfer function in Fourier space
    do i = 0, n-1
       f     = i*(self%samprate/2)/(n-1) 
       TF    = cmplx(splint(self%TF_real, f), splint(self%TF_imag, f))
       dv(i) = TF * dv(i)
    end do

    ! Inverse FFT
    call timer%start(TOT_FFT)
    call sfftw_execute_dft_c2r(plan_back, dv, dt)
    call timer%stop(TOT_FFT)
    tod = dt(1:ntod) / nfft

    deallocate(dt, dv)
    !call dfftw_destroy_plan(plan_fwd)
    !call dfftw_destroy_plan(plan_back)
    
  end subroutine convolve_Tbol

  ! Routine for convolving with bolometer transfer function, tod = T * tod
  subroutine update_Tbol(self, par)
    implicit none
    class(comm_Tbol),               intent(inout) :: self
    real(dp),         dimension(:), intent(in)    :: par

    integer(i4b) :: i, j
    real(dp)     :: omega
    complex(dpc) :: F, Hp
    real(dp),     allocatable, dimension(:)   :: nu
    complex(dpc), allocatable, dimension(:)   :: TF

    ! Update parameters
    self%npole     = count(par(3::2) > 0.d0)
    self%tau_stray = par(1)
    self%S_phase   = par(2)
    do i = 1, self%npole
       self%a(i)   = par(2+2*i-1)
       self%tau(i) = par(2+2*i)
    end do

    ! Generate transfer function
    allocate(nu(0:NSPLINE), TF(0:NSPLINE)) ! position 0 is spin frequency
    do i = 0, NSPLINE
       if (i == 0) then
          nu(i) = 1.d0/60.d0 ! One rotation per minute
       else
          nu(i) = (i-1)*(self%samprate/2)/(NSPLINE-1)
       end if
       omega = 2.d0*pi * nu(i)

       ! Compute bolometer transfer function
       F = 0.d0
       do j = 1, self%npole
          F = F + self%a(j) / cmplx(1.d0, omega*self%tau(j))
       end do

       ! Compute electronic transfer function
       Hp = self%get_HFI_electronic_transfunc(omega)

       ! Compute total transfer function
       TF(i) = F*Hp
       !write(*,fmt='(i6,5f10.6)') i, nu(i), omega, abs(F), abs(Hp), abs(TF(i))
    end do

    ! Normalize to unity amplitude at spin (dipole) frequency; require real component to be positive
    TF = TF / abs(TF(0))
    if (real(TF(0)) < 0.d0) TF = -TF
    
    ! Spline real and imaginary components separately
    call free_spline(self%TF_real)
    call free_spline(self%TF_imag)
    call spline(self%TF_real, nu(1:NSPLINE),  real(TF(1:NSPLINE)))
    call spline(self%TF_imag, nu(1:NSPLINE), aimag(TF(1:NSPLINE)))

!!$    open(58, file='tf.dat')
!!$    do i = 1, NSPLINE
!!$       write(58,*) nu(i), real(TF(i)), aimag(TF(i))
!!$    end do
!!$    write(58,*)
!!$    do i = 1, 100000
!!$       omega = (i-1)*(self%samprate/2)/(100000-1)
!!$       write(58,*) omega, splint(self%TF_real, omega), splint(self%TF_imag, omega)
!!$    end do
!!$    close(58)

    deallocate(nu, TF)
    
  end subroutine update_Tbol

  
  function get_HFI_electronic_transfunc(self, omega) result(Hp)
    implicit none
    class(comm_Tbol), intent(in)    :: self
    real(dp),         intent(in)    :: omega
    complex(dpc)                    :: Hp

    integer(i4b) :: k, Nsamp
    real(dp)     :: omega_r, dt, om_pos, om_neg, fmod
    complex(dpc) :: i, H_pos, H_neg

    Nsamp   = 40.d0                   ! Number of samples per modulation semi-cycle
    fmod    = 90.18759d0              ! Modulation frequency
    omega_r = 2.d0*pi * fmod          ! Angular modulation frequency
    dt      = self%S_phase/Nsamp/fmod ! Time shift
    i       = cmplx(0.d0, 1.d0)
    
    Hp = 0.d0
    do k = 0, 4
       om_pos = omega + (2*k+1) * omega_r
       om_neg = omega - (2*k+1) * omega_r
       H_pos  = self%get_HFI_Homega(om_pos)
       H_neg  = self%get_HFI_Homega(om_neg)
       
       Hp = Hp + 0.5d0 * exp(-i*(pi*omega/(2*omega_r)+omega*dt)) * &
            & (H_pos*exp(i*om_pos*dt)/((2*k+1)*om_pos) * (1-exp(i*om_pos*pi/omega_r)) - &
            &  H_neg*exp(i*om_neg*dt)/((2*k+1)*om_neg) * (1-exp(i*om_neg*pi/omega_r)))
    end do
       
  end function get_HFI_electronic_transfunc

  function get_HFI_Homega(self, omega) result(H)
    implicit none
    class(comm_Tbol), intent(in)    :: self
    real(dp),         intent(in)    :: omega
    complex(dpc)                    :: H

    real(dp)     :: R1, C1, R2, C2, R4, C4, R9, R12, C, R78, C18, K1, K2, K3
    complex(dpc) :: s, h0, h1, h2, h3, h4, h5

    R1  = 1.d3   ! Ohm;   Low pass filter
    C1  = 1.d-7  ! Farad; Lowpass filter
    R2  = 51.d3  ! Ohm;   Sallen Key highpass filter
    C2  = 1.d-6  ! Farad; Sallen Key highpass filter
    R4  = 10.d3  ! Ohm;   Single pole lowpass filter with gain
    C4  = 1.d-8  ! Farad; Single pole lowpass filter with gain
    R9  = 18.7d3 ! Ohm;   Single pole high pass filter + Sallen Key low pass filter
    R12 = 37.4d3 ! Ohm;   Single pole high pass filter + Sallen Key low pass filter
    R78 = 5.1d5  ! Ohm;   Single pole high pass filter + Sallen Key low pass filter
    C   = 1.d-8  ! Farad; Single pole high pass filter + Sallen Key low pass filter
    C18 = 1.d-6  ! Farad; Single pole high pass filter + Sallen Key low pass filter

    K3 = R9**2 * R78 * R12**2 * C**2 * C18
    K2 = R9 * R12**2 * R78 * C**2 + R9**2 * R12**2 * C**2 + R9 * R12**2 * R78 * C18 * C
    K1 = R9 * R12**2 * C + R12 * R78 * R9 * C18
    
    s  = cmplx(0.d0,omega)
    h0 = 1.d0/(1.d0 + self%tau_stray*s)           ! Stray capacitance lowpass filter
    h1 = (2.d0 + R1*C1*s) / (2.d0*(1.d0+R1*C1*s)) ! Lowpass filter
    h2 = (R2*C2*s)**2 / (1.d0 + R2*C2*s)**2       ! Sallen Key highpass
    h3 = -5.1d0                                   ! Sign reversal
    h4 = 1.5d0 / (1.d0 + R4*C4*s)                 ! Single pole LP + gain
    h5 = 2.d0*R12*R9*R78*C18*s / (s**3*K3 + s**2*K2 + s*K1 + R12*R9) ! SP HPF+Sallen Key

    H = h0*h1*h2*h3*h4*h5
 
  end function get_HFI_Homega

end module comm_tod_Tbol_mod
