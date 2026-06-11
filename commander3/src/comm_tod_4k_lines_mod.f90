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

module comm_tod_4k_lines_mod
  use comm_tod_mod
  use comm_status_mod
  use comm_tod_driver_mod
  implicit none
  private

  public comm_4k_lines

  type :: comm_4k_lines
    integer(i4b) :: window
    real(sp), dimension(:),   allocatable :: A_fit, f0_fit, sigma_fit, baseline
    real(sp), dimension(:,:), allocatable :: spike_profile
  contains
    procedure :: estimate_4k_baseline
    procedure :: A_lin_fit
  end type comm_4k_lines

  interface comm_4k_lines
    procedure cooler_4k_lines_constructor
  end interface comm_4k_lines

contains
 
  !initializes a comm_4K_lines class
  function cooler_4k_lines_constructor() result(c)
    implicit none
    class(comm_4K_lines), pointer :: c
    allocate(c)
  end function cooler_4k_lines_constructor
 
  subroutine estimate_4k_baseline(self,sub_ps)
     ! estimates the logaritmic baseline of the power spectrum in a 
     ! window around the spike
     implicit none
     class(comm_4k_lines),        intent(inout) :: self
     real(sp), dimension(1:,1:),  intent(inout) :: sub_ps

     integer(i4b) :: i, j, nsub, ex0, ex1
     real(sp) :: a, b, frac_excl, X, Y, XX, XY, meanX, meanY

     nsub = size(sub_ps(:,1))
     allocate(self%baseline(nsub))
     self%baseline = 0.d0
     
     frac_excl = 0.2 ! exclude +- 20% of the window around peak
     ex0 = nsub * (0.5 - frac_excl)
     ex1 = nsub * (0.5 + frac_excl)
     j = 0
     X = 0.d0; Y = 0.d0; XX = 0.d0; XY = 0.d0
     do i = 1, nsub
        ! BUGFIX: this exclusion test used .or., which is true for every i, so ALL samples were
        ! excluded and the log-linear fit below never ran (j stayed 0, always falling back to the
        ! constant baseline). Changed to .and. so only the central window [ex0,ex1] around the
        ! peak is excluded, as intended.
        if (i>=ex0 .and. i<=ex1) cycle
        X = X + log(sub_ps(i,1))
        Y = Y + log(sub_ps(i,2))
        XY = XY + log(sub_ps(i,1)) * log(sub_ps(i,2))
        XX = XX + log(sub_ps(i,1)) * log(sub_ps(i,1))
        j = j+1
     end do

     if (j<3) then ! baseline is constant
        a = log(max(1.d-30,(sub_ps(1,2) + sub_ps(nsub,2))/2))
        b = 0.d0
     else
        meanX = X/j
        meanY = Y/j
        ! BUGFIX: the least-squares slope denominator was XX - j*meanX*meanY; the correct
        ! expression is XX - j*meanX*meanX (= n*Var(x)). The old form gave a wrong baseline slope.
        b = (XY - j*meanX*meanY) / (XX - j*meanX*meanX)
        a = meanY - b*meanX
     end if

     do i = 1, nsub
        self%baseline(i) = exp(a + b*log(sub_ps(i,1)))
     end do
  end subroutine estimate_4k_baseline

  subroutine A_lin_fit(self,sub_ps,iter,n_f0,n_sig)
  ! fits amplitude of gaussian profile from fixing f0 and sigma
     implicit none
     class(comm_4k_lines),        intent(inout) :: self
     real(sp), dimension(1:,1:),  intent(inout) :: sub_ps
     integer(i4b),                intent(in)    :: iter, n_f0, n_sig

     integer(i4b) :: i, j, k, nsub, idx_peak
     real(sp) :: peak_val, f0_low, f0_high, sigma_low, sigma_high, chisq
     real(sp) :: A_try, f0_try, sigma_try, denom, numer, chisq_try
     real(sp), dimension(:), allocatable :: phi

     nsub = size(sub_ps(:,1))
     peak_val = 1.d-30; idx_peak = 1
     do i = 1, nsub
        if (sub_ps(i,2) > peak_val) then
             peak_val = sub_ps(i,2)
             idx_peak = i
        end if
     end do
    
     f0_low = sub_ps(idx_peak,1) - (sub_ps(nsub,1) - sub_ps(1,1))/2
     f0_low = max(sub_ps(1,1),f0_low)
     f0_high = sub_ps(idx_peak,1) + (sub_ps(nsub,1) - sub_ps(1,1))/2
     ! BUGFIX: this clamp previously assigned to f0_low instead of f0_high, overwriting the lower
     ! bound with a value near the upper bound and collapsing the f0 grid search to a single point.
     f0_high = min(sub_ps(nsub,1),f0_high)
     sigma_low = (sub_ps(2,1) - sub_ps(1,1))/5
     sigma_high = (sub_ps(2,1) - sub_ps(1,1))*2 

     ! Init chisq
     allocate(phi(nsub))
     self%A_fit(iter) = peak_val
     self%f0_fit(iter) = sub_ps(idx_peak,1)
     self%sigma_fit(iter) = 2.d-4 
     chisq = 0.d0
     do k = 1, nsub
        phi(k) = exp(-0.5 * ((sub_ps(k,1) - self%f0_fit(iter))/self%sigma_fit(iter))**2)
        ! BUGFIX: chisq was assigned ("chisq = ...") instead of accumulated, so the reference chisq
        ! only contained the last sample, and the grid search below almost never improved on it.
        chisq = chisq + (sub_ps(k,2) - self%A_fit(iter)*phi(k))**2
     end do


     do i = 0, n_f0
        f0_try = f0_low + (f0_high - f0_low)*i/n_f0
        do j = 0, n_sig
           sigma_try = sigma_low + (sigma_high - sigma_low)*j/n_sig
           
           ! Build phi
           denom = 0.d0; numer = 0.d0
           do k = 1, nsub
              phi(k) = exp(-0.5 * ((sub_ps(k,1) - f0_try)/sigma_try)**2)
              numer = numer + sub_ps(k,2) * phi(k)
              denom = denom + phi(k)**2
           end do

           if (denom <=0 .or. numer <=0) cycle
           A_try = numer / denom
           A_try = min(1.d5, A_try) ! avoid too strong correction
           chisq_try = 0.d0
           do k = 1, nsub
              chisq_try = chisq_try + (sub_ps(k,2) - A_try*phi(k))**2
           end do

           if (chisq_try < chisq) then
              chisq = chisq_try
              self%A_fit(iter) = A_try
              self%f0_fit(iter) = f0_try
              self%sigma_fit(iter) = sigma_try
           end if
        end do
     end do     

     deallocate(phi)

  end subroutine A_lin_fit

end module comm_tod_4k_lines_mod
