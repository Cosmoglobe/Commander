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
module comm_tod_gain_mod
  use comm_tod_mod
  use comm_tod_noise_mod
  use comm_utils
  use comm_fft_mod
  use math_tools
  implicit none

interface

  module subroutine calculate_gain_mean_std_per_scan(tod, scan_id, s_invsqrtN, mask, s_tot, handle, mask_lowres, tod_arr)
      ! 
      ! Calculate the scan-specific contributions to the gain estimation. 
      !
      ! This subroutine is called during the processing of each scan, and the
      ! results are stored in the TOD object. Then, after each scan has been
      ! processed in the main loop, these data are used to calculate the total
      ! gain estimate for each detector, since these gain estimates depend on
      ! knowing all the scan-specific gain contributions. The data will be
      ! downsampled before calculating the estimate.
      !
      ! Arguments:
      ! ----------
      ! tod:           derived class (comm_tod)
      !                TOD object containing all scan-relevant data. Will be
      !                modified.
      ! scan_id:       integer(i4b)
      !                The ID of the current scan.
      ! s_invsqrtN:    real(sp) array
      !                The product of the reference signal we want to calibrate
      !                on multiplied by the square root of the inverse noise
      !                matrix for this scan.
      ! handle:        derived class (planck_rng)
      !                Random number generator handle. Will be modified.
      ! mask_lowres:   real(sp) array
      !                The mask for the low-resolution downsampled data.
      !
      ! tod_arr:       To be deprecated soon.
      !
      ! Returns:
      ! --------
      ! tod:           derived class (comm_tod)
      !                Will update the fields tod%scans(scan_id)%d(:)%dgain and
      !                tod%scans(scan_id)%d(:)%gain_invsigma, which contain
      !                incremental estimates to be used for global gain
      !                estimation. 
      !                dgain = s_ref * n_invN * residual where residual is
      !                     d - (g_0 + g_i) * s_tot (g_0 being the absolute gain, and
      !                     g_i the absolute gain per detector)
      !                gain_invsigma =  s_ref * n_invN * s_ref

      
    implicit none
    class(comm_tod),                      intent(inout) :: tod
    real(sp),             dimension(:,:), intent(in)    :: tod_arr
    real(sp),             dimension(:,:), intent(in)    :: s_invsqrtN, mask, s_tot
    integer(i4b),                         intent(in)    :: scan_id
    type(planck_rng),                     intent(inout) :: handle
    real(sp),             dimension(:,:), intent(in), optional :: mask_lowres
  end subroutine calculate_gain_mean_std_per_scan


! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  ! Eirik: Update this routine to sample time-dependent gains properly; results should be stored in self%scans(i)%d(j)%gain, with gain0(0) and gain0(i) included
  module subroutine sample_smooth_gain(tod, handle, dipole_mods, smooth)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),   dimension(:, :),       intent(in)     :: dipole_mods
    logical(lgt), optional,            intent(in)     :: smooth
  end subroutine sample_smooth_gain

  module subroutine normalize_gain_variance(g1, g2, sigma0)
    implicit none
    real(dp), dimension(:), intent(inout) :: g1
    real(dp), dimension(:), intent(inout) :: g2
    real(dp),               intent(out)   :: sigma0
  end subroutine normalize_gain_variance


   ! This is implementing equation 16, adding up all the terms over all the sums
   ! the sum i is over the detector.
  module subroutine accumulate_abscal(tod, scan, mask, s_sub, s_invsqrtN, A_abs, b_abs, handle, out, s_highres, mask_lowres, tod_arr)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    integer(i4b),                      intent(in)     :: scan
    real(sp),          dimension(:,:), intent(in)     :: mask, s_sub
    real(sp),          dimension(:,:), intent(in)     :: s_invsqrtN
    real(dp),          dimension(:),   intent(inout)  :: A_abs, b_abs
    type(planck_rng),                  intent(inout)  :: handle
    logical(lgt), intent(in) :: out
    real(sp),          dimension(:,:), intent(in), optional :: s_highres
    real(sp),          dimension(:,:), intent(in), optional :: mask_lowres
    real(sp),          dimension(:,:), intent(in), optional :: tod_arr
  end subroutine accumulate_abscal

  ! Sample absolute gain from orbital dipole alone 
  module subroutine sample_abscal_from_orbital(tod, handle, A_abs, b_abs)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs
  end subroutine sample_abscal_from_orbital

  ! Sample absolute gain from orbital dipole alone 
  ! Eirik: Change this routine to sample the constant offset per radiometer; put the result in tod%gain0(i)
  module subroutine sample_relcal(tod, handle, A_abs, b_abs)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs
  end subroutine sample_relcal

  module subroutine sample_imbal_cal(tod, handle, A_abs, b_abs)
    !  Subroutine to sample the transmission imbalance parameters, defined in
    !  the WMAP data model as the terms x_im; given the definition
    !  d_{A/B} = T_{A/B} \pm Q_{A/B} cos(2 gamma_{A/B}) \pm U_{A/B} sin(2 gamma_{A/B})
    !  we have
    !  d = g[(1+x_im)*d_A - (1-x_im)*d_B]
    !    = g(d_A - d_B) + g*x_im*(d_A + d_B)
    !  Returns x_{im,1} for detectors 13/14, and x_{im,2} for detectors 23/24.
    !
    !
    !  Arguments (fixed):
    !  ------------------
    !  A_abs: real(dp)
    !     Accumulated A_abs = s_ref^T N^-1 s_ref for all scans
    !  b_abs: real(dp)
    !     Accumulated b_abs = s_ref^T N^-1 s_sub for all scans
    !
    !  
    !  Arguments (modified):
    !  ---------------------
    !  tod: comm_WMAP_tod
    !     The entire tod object. tod%x_im estimated and optionally sampled
    !  handle: planck_rng derived type 
    !     Healpix definition for random number generation
    !     so that the same sequence can be resumed later on from that same point
    !
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs
  end subroutine sample_imbal_cal

  module subroutine get_smoothing_windows(tod, windows, dipole_mods)
     implicit none

     integer(i4b), dimension(:, :), intent(out)         :: windows
     class(comm_tod),               intent(in)          :: tod
     real(dp), dimension(:, :),     intent(in)          :: dipole_mods

  end subroutine get_smoothing_windows


  module subroutine wiener_filtered_gain(b, inv_N_wn, sigma_0, alpha, fknee, sample, &
     & handle)
     !
     ! Given a spectral model for the gain, samples a wiener filtered (smoothed)
     ! estimate of the gain given data. The basic equation to be solved for is
     !  Ax = b where
     !  A = inv_G + s_ref * inv_N_wn * s_ref, where inv_G is the covariance
     !                                        matrix of the gain, given by the
     !                                        spectral model.
     !  and 
     !  b = s_ref * inv_N_wn * residual.
     !
     ! If a sample is to be drawn, a fluctuation term will also be added.
     !
     ! The spectral model is given by P(f) = sigma_0**2 * (f/f_knee) ** alpha
     !
     ! Arguments:
     ! ----------
     ! b:           real(dp) array
     !              The "right hand side" of the equation Ax = b, where x is the
     !              wiener filtered solution. b = dgain from
     !              calculate_gain_mean_std_per_scan. Will contain the solution.
     ! inv_N_wn:    real(dp) array
     !              The inverse white noise diagonal covariance matrix
     ! sigma_0:     real(dp)
     !              The current estimate of sigma_0 in the spectral model
     ! alpha:       real(dp)
     !              The current estimate of alpha in the spectral model
     ! fknee:       real(dp)
     !              The current estimate of fknee in the spectral model
     ! sample:      logical(lgt)
     !              Whether to draw a sample from the distribution or just
     !              return the best-fit solution.
     ! handle:      derived class (planck_rng)
     !              Random number generator handle. Will be modified.
     !
     ! Returns:
     ! --------
     ! b:           real(dp) array
     !              At exit, will contain the solution to Ax = b, including
     !              fluctuations if sampling is turned on.

     implicit none

     real(dp), dimension(:), intent(inout)   :: b
     real(dp), dimension(:), intent(in)      :: inv_N_wn
     real(dp), intent(in)                    :: sigma_0, alpha, fknee
     logical(lgt), intent(in)                :: sample
     type(planck_rng)                        :: handle
  end subroutine wiener_filtered_gain

  module subroutine fft_fwd(dt, dv, plan_fwd)
     !
     ! Given a fft plan, transforms a vector from time domain to Fourier domain.
     !
     ! Arguments:
     ! ----------
     ! vector:          real(dp) array
     !                  The original time-domain vector.
     ! fourier_vector:  complex(dpc) array
     !                  The array that will contain the Fourier vector.
     ! plan_fwd:        integer*8
     !                  The fft plan to carry out the transformation
     !
     ! Returns:
     ! --------
     ! fourier_vector:  complex(dpc) array
     !                  At exit will contain the Fourier domain vector.

     implicit none
     real(dp),     dimension(:), intent(in)     :: dt
     complex(dpc), dimension(:), intent(out)    :: dv
     integer*8,                  intent(in)     :: plan_fwd
   end subroutine fft_fwd

  module subroutine fft_back(dv, dt, plan_back)
     !
     ! Given a fft plan, transforms a vector from Fourier domain to time domain.
     !
     ! Arguments:
     ! ----------
     ! dv:              complex(dpc) array
     !                  The original Fourier domain vector.
     ! dt:              real(dp) array
     !                  The vector that will contain the transformed vector.
     ! plan_back:       integer*8
     !                  The fft plan to carry out the transformation.
     !
     ! Returns:
     ! --------
     ! vector:          real(dp) array
     !                  At exit will contain the time-domain vector.
     implicit none
     complex(dpc), dimension(:), intent(in)     :: dv
     real(dp),     dimension(:), intent(out)    :: dt
     integer*8,                  intent(in)     :: plan_back
   end subroutine fft_back

  module function solve_cg_gain(inv_N_wn, inv_N_corr, b, precond, plan_fwd, plan_back)!, &
!     & with_precond)
     !
     ! Specialized function for solving the gain Wiener filter equation using
     ! the CG algorithm.
     !
     ! Arguments:
     ! ----------
     ! inv_N_wn:    real(dp) array
     !              The inverse white noise covariance matrix in time domain.
     ! inv_N_corr:  real(dp) array
     !              The inverse gain covariance matrix in Fourier domain.
     ! b:           real(dp) array
     !              The right hand side of the Wiener filter equation.
     ! precond:     real(dp) array
     !              The preconditioner matrix in time domain.
     ! plan_fwd:    integer*8
     !              The FFT forward plan.
     ! plan_back:   integer*8
     !              The FFT backward plan.
     !
     ! Returns:
     ! --------
     ! solve_cg_gain real(dp) array
     !               The solution to the Wiener filter equation.
     implicit none
     real(dp), dimension(:), intent(in) :: inv_N_wn, inv_N_corr, b, precond
     integer*8             , intent(in) :: plan_fwd, plan_back
     real(dp), dimension(size(b))       :: solve_cg_gain
  end function solve_cg_gain

  module subroutine apply_cg_precond(vec_in, vec_out, precond, plan_fwd, plan_back)
    implicit none
    real(dp), dimension(:), intent(in)  :: vec_in, precond
    real(dp), dimension(:), intent(out) :: vec_out
    integer*8             , intent(in)  :: plan_fwd, plan_back    
  end subroutine apply_cg_precond

  module function tot_mat_mul_by_vector(time_mat, fourier_mat, vector, plan_fwd, &
     & plan_back, filewrite)
     !
     ! Multiplies a sum of a time-domain diagonal matrix and a Fourier-domain
     ! diagonal matrix by a time-domain vector.
     !
     ! Arguments:
     ! ----------
     ! time_mat:    real(dp) array
     !              The matrix that is diagonal in time domain.
     ! fourier_mat: real(dp) array
     !              The matrix that is diagonal in Fourier domain.
     ! vector:      real(dp) array
     !              The vector to multiply by.
     ! plan_fwd:    integer*8
     !              The FFT forward plan.
     ! plan_back:   integer*8
     !              The FFT backward plan.
     ! filewrite:   logical(lgt), optional
     !              If present, writes various debug info into files.
     !
     ! Returns:
     ! --------
     ! tot_mat_mul_by_vector:   real(dp) array
     !                          The time-domain result of multiplying the
     !                          matrices by the input vector.
     implicit none

     real(dp), dimension(:)                            :: vector
     real(dp), dimension(size(vector))                 :: tot_mat_mul_by_vector
     real(dp), dimension(size(vector)), intent(in)     :: time_mat
     real(dp), dimension(:) , intent(in)     :: fourier_mat
     integer*8              , intent(in)     :: plan_fwd, plan_back
     logical(lgt)           , optional       :: filewrite
  end function tot_mat_mul_by_vector

  module subroutine calculate_invcov(sigma_0, alpha, fknee, freqs, invcov) 
     !
     ! Calculate the inverse covariance matrix of the gain PSD given the PSD
     ! parameters as a diagonal matrix in Fourier domain.
     !
     ! The spectral model is given by P(f) = sigma_0**2 * (f/f_knee) ** alpha
     !
     ! Arguments:
     ! ----------
     ! sigma_0:     real(dp)
     !              The sigma_0 parameter of the spectral model. 
     ! alpha:       real(dp)
     !              The alpha parameter of the spectral model. 
     ! fknee:       real(dp)
     !              The fknee parameter of the spectral model. 
     ! freqs:       real(dp) array
     !              The frequencies at which to calculate the covariance matrix.
     !
     ! Returns:
     ! --------
     ! calculate_invcov:    real(dp) array
     !                      The covariance matrix in Fourier space, evaluated at
     !                      the frequencies given by 'freqs'.

     implicit none
     real(dp), dimension(:), intent(in)     :: freqs
     real(dp), intent(in)                   :: sigma_0, fknee, alpha
     real(dp), dimension(:), intent(out)    :: invcov
  end subroutine calculate_invcov

  module function calc_sigma_0(gain)
     !
     ! Function to estimate sigma_0 in the gain PSD model. Sigma_0 is estimated
     ! as the standard deviation of
     !      res(i) = gain(i) - gain(i-1)
     ! over all scans.
     !
     ! Arguments:
     ! ----------
     ! gain:        real(dp) array
     !              The gains per scan.
     !
     ! Returns:
     ! --------
     ! calc_sigma_0: real(dp)
     !               The estimated sigma_0.
     implicit none
     real(dp)       :: calc_sigma_0
     real(dp), dimension(:) :: gain
  end function calc_sigma_0

  module subroutine sample_gain_psd(tod, handle)
     !
     ! "Master" subroutine for gain PSD sampling. Uses the current gain solution
     ! and current values for the PSD parameters to sample new values for these
     ! parameters.
     !
     ! Arguments:
     ! ----------
     ! tod:         derived class (comm_tod)
     !              The TOD object that contains all relevant data. The fields
     !              containing the gain PSD parameters will be updated.
     ! handle:      derived class (planck_rng)
     !              Random number generator handle. Will be modified.
     !
     ! Returns:
     ! --------
     ! tod will be updated with a new set of values for the gain PSD model,
     ! drawn from their conditional distribution given the current gain
     ! solution.
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
  end subroutine sample_gain_psd

  module subroutine sample_psd_params_by_mh(gain_ps, freqs, sigma_0, fknee, alpha, &
     & sigma0_std, fknee_std, alpha_std, handle)
    ! 
    ! Uses the Metropolis-Hastings algorithm to draw a sample of the gain PSD
    ! parameters given the current gain solution.
    !
    ! The subroutine will do an initial run to estimate a proposal covariance
    ! matrix, then use that covariance matrix in the final run. The last sample
    ! of the chain is used as the drawn sample.
    !
    ! Arguments:
    ! ----------
    ! gain_ps:      real(dp) array
    !               The current gain power spectrum in fourier space.
    ! freqs:        real(dp) array
    !               The frequencies at which the power spectrum is calculated.
    ! sigma_0:      real(dp)
    !               The sigma_0 parameter of the spectral model. 
    ! fknee:        real(dp)
    !               The fknee parameter of the spectral model. Will be modified.
    ! alpha:        real(dp)
    !               The alpha parameter of the spectral model. Will be modified.
    ! fknee_std:    real(dp)
    !               The initial proposal standard deviation for fknee - will be
    !               re-estimated during the sampling.
    ! alpha_std:    real(dp)
    !               The initial proposal standard deviation for alpha - will be
    !               re-estimated during the sampling.
    ! handle:       derived class (planck_rng)
    !               Random number generator handle. Will be modified.
    !
    ! Returns:
    !---------
    ! fknee and alpha will contain new values sampled from the posterior given
    ! by the current gain power spectrum.
    implicit none

    real(dp), dimension(:), intent(in)      :: gain_ps
    real(dp), dimension(:), intent(in)      :: freqs
    real(dp), intent(inout)                 :: sigma_0, fknee, alpha
    real(dp), intent(in)                    :: sigma0_std, fknee_std, alpha_std
    type(planck_rng),                  intent(inout)  :: handle
  end subroutine sample_psd_params_by_mh

  module subroutine run_mh(sigma_0, fknee, alpha, propcov, num_samples, samples, &
     & gain_ps, freqs, adjust_scaling_full, handle)
  !
  ! The core Metropolis-Hastings routine used for gain PSD sampling.
  !
  ! Arguments:
  ! ----------
  ! fknee:                  real(dp)
  !                         The current value of fknee.
  ! alpha:                  real(dp)
  !                         The current value of alpha.
  ! propcov:                real(dp) 2x2 matrix
  !                         The covariance matrix of the proposal density.
  ! sigma_0:                real(dp)
  !                         The sigma_0 parameter.
  ! num_samples:            integer(i4b)
  !                         The number of samples to draw.
  ! samples:                real(dp) array
  !                         Output array for the samples.
  ! gain_ps:                real(dp) array
  !                         The current gain power spectrum.
  ! freqs:                  real(dp) array
  !                         The frequencies at which the gain power spectrum is
  !                         evaluated.
  ! adjust_scaling_full:    logical(lgt)
  !                         If true, will continue to adjust the scaling
  !                         parameter (the step length) throughout. If false,
  !                         will only adjust the step length during the first
  !                         half of the sampling.
  ! handle:                 derived class (planck_rng)
  !                         Random number generator handle. Will be modified.
  !
  ! Returns:
  ! --------
  ! the 'samples' array will contain num_samples samples drawn from the
  ! posterior distribution of the gain PSD given the current gain power
  ! spectrum. Only the last half of the samples are strictly speaking
  ! statistically proper samples, and only if adjust_scaling_full is false.
    implicit none

    real(dp), intent(in)       :: fknee, alpha, sigma_0
    real(dp), dimension(3, 3), intent(in)  :: propcov
    integer(i4b), intent(in)       :: num_samples
    real(dp), dimension(:, :), intent(out) :: samples
    real(dp), dimension(:),    intent(in)  :: freqs  
    real(dp), dimension(:), intent(in)  :: gain_ps
    logical(lgt)          ,    intent(in)  :: adjust_scaling_full
    type(planck_rng),                  intent(inout)  :: handle
  end subroutine run_mh

  module function psd_loglike(sigma_0, fknee, alpha, gain_ps, freqs)
     !
     ! Calculates the log-likelihood given the gain and psd parameters.
     !
     ! The spectral model is given by P(f) = sigma_0**2 * (f/f_knee) ** alpha,
     ! and the log-likelihood is given by taking the sum over all frequencies of
     ! the gain power spectrum times the inverse of the covariance matrix given
     ! by the spectral model, minus the logarithm of this covariance matrix.
     !
     ! Arguments:
     ! ----------
     ! fknee:       real(dp)
     !              The value of the fknee parameter.
     ! alpha:       real(dp)
     !              The value of the alpha parameter.
     ! sigma_0:     real(dp)
     !              The value of sigma_0
     ! gain_ps:     real(dp) array
     !              The power spectrum of the current gain solution.
     ! freqs:       real(dp) array
     !              The frequencies at which the gain power spectrum is
     !              evaluated.
     !
     ! Returns:
     ! --------
     ! psd_loglike      real(dp)
     !                  The log-likelihood value of the parameter combination
     !                  given the current gain power spectrum.
     implicit none

     real(dp)        :: psd_loglike
     real(dp)        :: fknee, alpha, sigma_0
     real(dp), dimension(:)  :: gain_ps, freqs
  end function psd_loglike

end interface    

end module comm_tod_gain_mod
