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
  use comm_fft_mod
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
     type(spline_type) :: TF_full_real,    TF_full_imag
     type(spline_type) :: TF_preproc_real, TF_preproc_imag
     type(spline_type) :: TF_cgmap_real,   TF_cgmap_imag
  contains
    procedure :: estimate_params
    procedure :: get_HFI_electronic_transfunc
    procedure :: get_HFI_Homega
    procedure :: get_K_regularization
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
  subroutine convolve_Tbol(self, tbol_type, operation, tod, plan_fwd, plan_back)
    implicit none
    class(comm_Tbol),               intent(in)    :: self
    character(len=*),               intent(in)    :: tbol_type
    character(len=*),               intent(in)    :: operation
    real(sp),         dimension(:), intent(inout) :: tod
    type(C_PTR),                    intent(inout) :: plan_fwd, plan_back
    
    integer(i4b) :: i, j, k, l, n, ntod, nomp, nfft, err

    real(dp)     :: f
    complex(dpc) :: TF
    real(sp),     allocatable, dimension(:)   :: dt
    complex(spc), allocatable, dimension(:)   :: dv

    ntod     = size(tod)
    nfft     = get_closest_fft_pow2(ntod)
    n        = nfft / 2 + 1

    allocate(dt(nfft), dv(0:n-1))

    ! FFT
    dt(1:ntod)           = tod
    dt(ntod+1:nfft)      = 0.    ! Zero pad
    call timer%start(TOT_FFT)
    call fftwf_execute_dft_r2c(plan_fwd, dt, dv)
    call timer%stop(TOT_FFT)

    ! Apply transfer function in Fourier space
    do i = 0, n-1
       f     = i*(self%samprate/2)/(n-1)
       if (trim(tbol_type) == 'full') then
          TF    = cmplx(splint(self%TF_full_real, f), splint(self%TF_full_imag, f))
       else if (trim(tbol_type) == 'preproc') then
          TF    = cmplx(splint(self%TF_preproc_real, f), splint(self%TF_preproc_imag, f))
       else if (trim(tbol_type) == 'cgmap') then
          TF    = cmplx(splint(self%TF_cgmap_real, f), splint(self%TF_cgmap_imag, f))
       end if

       if (trim(operation) == 'transpose') then
          TF = conjg(TF)
       else if (trim(operation) == 'inverse') then
          if (abs(TF) > 0.d0) then 
             TF = 1.d0 / TF
          else
             TF = 0.d0
          end if
       end if
       dv(i) = TF * dv(i)
    end do

    ! Inverse FFT
    call timer%start(TOT_FFT)
    call fftwf_execute_dft_c2r(plan_back, dv, dt)
    call timer%stop(TOT_FFT)
    tod = dt(1:ntod) / nfft

    deallocate(dt, dv)
    
  end subroutine convolve_Tbol

  ! Routine for convolving with bolometer transfer function, tod = T * tod
  subroutine update_Tbol(self, par)
    implicit none
    class(comm_Tbol),               intent(inout) :: self
    real(dp),         dimension(:), intent(in)    :: par

    integer(i4b) :: i, j
    real(dp)     :: omega, K
    complex(dpc) :: F, Hp
    real(dp),     allocatable, dimension(:)   :: nu
    complex(dpc), allocatable, dimension(:)   :: TF_full, TF_preproc, TF_cgmap

    ! Update parameters
    self%npole     = count(par(3::2) > 0.d0)
    self%tau_stray = par(1)
    self%S_phase   = par(2)
    do i = 1, self%npole
       self%a(i)   = par(2+2*i-1)
       self%tau(i) = par(2+2*i)
    end do

    ! Generate transfer function; position 0 is spin frequency
    allocate(nu(0:NSPLINE), TF_full(0:NSPLINE), TF_preproc(0:NSPLINE), TF_cgmap(0:NSPLINE))
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

       ! Compute regularization kernel
       K = self%get_K_regularization(nu(i))
       
       ! Compute total transfer function
       if (K > 0.d0) then
          TF_full(i)    = K
          TF_preproc(i) = K/(F*Hp)
       else
          TF_full(i)    = 0.d0
          TF_preproc(i) = 0.d0
       end if
       !write(*,fmt='(i6,5f10.6)') i, nu(i), omega, abs(F), abs(Hp), abs(TF(i))
    end do

    ! Normalize to unity amplitude at spin (dipole) frequency; require real component to be positive
    TF_full    = TF_full    / abs(TF_full(0))
    TF_preproc = TF_preproc / abs(TF_preproc(0))
    TF_cgmap   = 1.d0 !TF_full/TF_preproc 
    if (real(TF_preproc(0)) < 0.d0) TF_preproc = -TF_preproc

    ! Spline real and imaginary components separately
    call free_spline(self%TF_full_real)
    call free_spline(self%TF_full_imag)
    call spline(self%TF_full_real, nu(1:NSPLINE),  real(TF_full(1:NSPLINE)))
    call spline(self%TF_full_imag, nu(1:NSPLINE), aimag(TF_full(1:NSPLINE)))
    
    ! Regularize high-frequency real component for preprocessing, avoid noise boost
!!$    TF_preproc = TF_full
!!$    i = 1
!!$    do while (i < NSPLINE .and. (real(TF_preproc(i)) > real(TF_preproc(i+1)) .or. nu(i) < 80.d0))
!!$       i = i+1
!!$    end do
!!$    ! Set real component at hifher frequencies to local minimum
!!$    TF_preproc(i:) = cmplx(real(TF_preproc(i)), aimag(TF_preproc(i:)))

    ! Spline real and imaginary components separately
    call free_spline(self%TF_preproc_real)
    call free_spline(self%TF_preproc_imag)
    call spline(self%TF_preproc_real, nu(1:NSPLINE),  real(TF_preproc(1:NSPLINE)))
    call spline(self%TF_preproc_imag, nu(1:NSPLINE), aimag(TF_preproc(1:NSPLINE)))

    ! Compute residual ratio, to be left for CG mapmaking
    call free_spline(self%TF_cgmap_real)
    call free_spline(self%TF_cgmap_imag)
    call spline(self%TF_cgmap_real, nu(1:NSPLINE),  real(TF_cgmap(1:NSPLINE)))
    call spline(self%TF_cgmap_imag, nu(1:NSPLINE), aimag(TF_cgmap(1:NSPLINE)))     

!!$    open(58, file='tf.dat', recl=1024)
!!$    do i = 1, NSPLINE
!!$       write(58,*) nu(i), real(TF_full(i)), aimag(TF_full(i)), real(TF_preproc(i)), aimag(TF_preproc(i)), real(TF_cgmap(i)), aimag(TF_cgmap(i))
!!$    end do
!!$    close(58)
!!$    write(*,*) 'done'

    deallocate(nu, TF_full, TF_preproc, TF_cgmap)
    
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
    !H = h0*h1*h3*h4
 
  end function get_HFI_Homega

  function get_K_regularization(self, f) result(K)
    implicit none
    class(comm_Tbol), intent(in)    :: self
    real(dp),         intent(in)    :: f
    real(dp)                        :: K

    real(dp) :: f_gauss, f_c, k0, fmax, fsamp

    fsamp   = 2*90.18759d0
    f_gauss = 65.d0 ! Hz
    f_c     = 80.d0 ! Hz
    k0      = 0.9
    fmax    = f_c + k0*(fsamp/2-f_c)

    if (f > fmax) then
       K = 0.d0
    else
       K = exp(-0.5d0 * (f/f_gauss)**2)
       if (f > f_c) then
          K = K * cos(0.5d0*pi * (f-f_c)/(fmax-f_c))**2
       end if
    end if

    !write(*,*) f, K
    
  end function get_K_regularization
  
end module comm_tod_Tbol_mod
