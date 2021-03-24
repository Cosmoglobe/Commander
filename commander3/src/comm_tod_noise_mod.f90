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
module comm_tod_noise_mod
  use comm_tod_mod
  use comm_utils
  use comm_fft_mod
  use InvSamp_mod
  use comm_tod_noise_psd_mod
  implicit none


contains

  subroutine sample_n_corr(self, tod, handle, scan, mask, s_sub, n_corr, pix, freqmask, dospike)
    ! 
    ! Routine for sample TOD-domain correlated noise given a pre-computed noise PSD, as defined by
    !    ((N_c^-1 + N_wn^-1) n_corr = d_prime + w1 * sqrt(N_wn) + w2 * sqrt(N_c) 
    ! 
    ! Arguments
    ! ---------
    ! self:    derived type (comm_tod)
    !          TOD object used for meta-data
    ! tod:     sp (ntod x ndet array)
    !          Decompressed TOD data
    ! handle:  type(planck_rng)
    !          Healpix random number type
    ! scan:    int (scalar)
    !          Scan number ID for local core
    ! mask:    sp (ntod x ndet array)
    !          TOD mask (0 = masked, 1 = included)
    ! s_sub:   sp (ntod x ndet array)
    !          TOD-domain signal template to be subtracted from tod to generate a residual
    ! pix:     int (ntod x ndet array)
    !          Pointing array in terms of pixel number per sample
    ! freqmask: sp (nfreq x ndet array)
    !          TOD frequency mask used to remove bad frequencies; these are re-filled with random noise
    ! dospike: lgt (scalar)
    !          Flag to identify spikes and output debug information (HKE: This is currently very non-intuitive...)
    !
    ! Returns
    ! -------
    ! n_corr:  sp (ntod x ndet array)
    !          TOD-domain correlated noise realization
    ! 
    implicit none
    class(comm_tod),                    intent(in)     :: self
    real(sp),         dimension(1:,1:), intent(in)     :: tod
    type(planck_rng),                   intent(inout)  :: handle
    integer(i4b),                       intent(in)     :: scan
    integer(i4b),     dimension(1:,1:), intent(in)     :: pix
    real(sp),         dimension(1:,1:), intent(in)     :: mask, s_sub
    real(sp),         dimension(1:,1:), intent(out)    :: n_corr
    real(sp),         dimension(0:,1:), intent(in), optional :: freqmask
    logical(lgt),                       intent(in), optional :: dospike

    integer(i4b) :: i, j, l, k, n, m, nomp, ntod, ndet, err, omp_get_max_threads
    integer(i4b) :: nfft, nbuff, j_end, j_start
    integer*8    :: plan_fwd, plan_back
    logical(lgt) :: init_masked_region, end_masked_region, pcg_converged
    real(sp)     :: sigma_0, alpha, nu_knee,  samprate, gain, mean, N_wn, N_c, nu
    real(dp)     :: power, fft_norm
    character(len=1024) :: filename
    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv
    real(sp),     allocatable, dimension(:) :: d_prime, ncorr2

    ntod     = self%scans(scan)%ntod
    ndet     = self%ndet
    nomp     = 1 !omp_get_max_threads()
    samprate = self%samprate
    nfft = get_closest_fft_magic_number(ceiling(ntod * 1.15d0))
    !nfft     = 2 * ntod
    fft_norm = sqrt(1.d0 * nfft)  ! used when adding fluctuation terms to Fourier coeffs (depends on Fourier convention)
    n        = nfft / 2 + 1

    call sfftw_init_threads(err)
    call sfftw_plan_with_nthreads(nomp)

    allocate(dt(nfft), dv(0:n-1), d_prime(ntod), ncorr2(ntod))
    call sfftw_plan_dft_r2c_1d(plan_fwd,  nfft, dt, dv, fftw_estimate + fftw_unaligned)
    call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)

    do i = 1, ndet
       if (.not. self%scans(scan)%d(i)%accept) cycle
       gain     = self%scans(scan)%d(i)%gain  ! Gain in V / K
       sigma_0  = abs(self%scans(scan)%d(i)%N_psd%sigma0)
       N_wn     = sigma_0**2  ! white noise power spectrum

       ! Prepare TOD residual
       d_prime = tod(:,i) - self%scans(scan)%d(i)%baseline - gain * S_sub(:,i)

       ! Fill gaps in data 
       init_masked_region = .true.
       end_masked_region  = .false.
       do j = 1, ntod
          if (mask(j,i) == 1.) then
             if (end_masked_region) then
                j_end = j - 1
                call fill_masked_region(d_prime, mask(:,i), j_start, j_end, ntod, self%scans(scan)%chunk_num)
                ! Add noise to masked region
                if (trim(self%operation) == "sample") then
                   do k = j_start, j_end
                      d_prime(k) = d_prime(k) + sigma_0 * rand_gauss(handle)
                   end do
                end if
                end_masked_region = .false.
                init_masked_region = .true.
             end if
          else
             if (init_masked_region) then
                init_masked_region = .false.
                end_masked_region = .true.
                j_start = j
             end if
          end if
       end do
       ! if the data ends with a masked region
       if (end_masked_region) then
          j_end = ntod
          call fill_masked_region(d_prime, mask(:,i), j_start, j_end, ntod, self%scans(scan)%chunk_num)
          if (trim(self%operation) == "sample") then
             do k = j_start, j_end
                d_prime(k) = d_prime(k) + sigma_0 * rand_gauss(handle)
             end do
          end if      
       end if

       ! Identify spikes
       !if (self%first_call .and. dospike) call find_d_prime_spikes(self, scan, i, d_prime, pix)

       !alpha    = self%scans(scan)%d(i)%N_psd%alpha
       !nu_knee  = self%scans(scan)%d(i)%N_psd%fknee

       
       pcg_converged = .false.
       call get_ncorr_sm_cg(handle, d_prime, ncorr2, mask(:,i), self%scans(scan)%d(i)%N_psd, samprate, nfft, plan_fwd, plan_back, pcg_converged, self%scanid(scan), i, trim(self%freq))
       n_corr(:,i) = ncorr2(:)

       if (.not. pcg_converged) then
          ! Preparing for fft
          dt(1:ntod)           = d_prime(:)
          !dt(2*ntod:ntod+1:-1) = dt(1:ntod)

          nbuff = nfft - ntod
          do j=1, nbuff
             dt(ntod+j) = sum(d_prime(ntod-20:ntod)) / 20.0 + (sum(d_prime(1:20)) - sum(d_prime(ntod-20:ntod))) / 20.0 * (j-1) / (nbuff - 1)
          end do

          call sfftw_execute_dft_r2c(plan_fwd, dt, dv)

          if (trim(self%operation) == "sample") then
             dv(0)    = dv(0) + fft_norm * sqrt(N_wn) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0)
          end if

          do l = 1, n-1                                                      
             nu  = l*(samprate/2)/(n-1)
             N_c = self%scans(scan)%d(i)%N_psd%eval_corr(nu)  ! correlated noise power spectrum

             ! Add random noise to masked frequencies
             if (present(freqmask)) then
                if (freqmask(l,i) == 0.) then
                   dv(l) = fft_norm * sqrt(N_c) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0) 
                end if
             end if

             ! Apply Wiener filter/constrained realization solution
             if (trim(self%operation) == "sample") then
                dv(l) = (dv(l) + fft_norm * (sqrt(N_wn) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0) &
                               + N_wn * sqrt(1.0 / N_c) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0))) &
                        & / (1.0 + N_wn/N_c)
             else
                dv(l) = dv(l) * 1.0/(1.0 + N_wn/N_c)
             end if
          end do
          call sfftw_execute_dft_c2r(plan_back, dv, dt)
          dt          = dt / nfft
          n_corr(:,i) = dt(1:ntod) 
       end if


       if (.true. .and. .true.) then !mod(self%scanid(scan),100) == 1) then
         write(filename, "(A, I0.3, A, I0.3, 3A)") 'ncorr_tods_new/ncorr_times', self%scanid(scan), '_', i, '_',trim(self%freq),'_test.dat' 
         open(65,file=trim(filename),status='REPLACE')
         do j = 1, ntod
            write(65, '(14(E15.6E3))') n_corr(j,i), s_sub(j,i), mask(j,i), d_prime(j), real(tod(j,i),sp), self%scans(scan)%d(i)%gain, self%scans(scan)%d(i)%N_psd%xi_n(3), self%scans(scan)%d(i)%N_psd%xi_n(2), self%scans(scan)%d(i)%N_psd%sigma0, self%scans(scan)%d(i)%N_psd%P_active(3,1), self%scans(scan)%d(i)%N_psd%P_active(2,1), self%samprate, ncorr2(j), nfft 
            ! if (present(tod_arr)) then
            !   write(65, '(14(E15.6E3))') n_corr(j,i), s_sub(j,i), mask(j,i), d_prime(j), real(tod_arr(j,i),sp), self%scans(scan)%d(i)%gain, self%scans(scan)%d(i)%N_psd%alpha, self%scans(scan)%d(i)%N_psd%fknee, self%scans(scan)%d(i)%N_psd%sigma0, self%scans(scan)%d(i)%N_psd%alpha_def, self%scans(scan)%d(i)%N_psd%fknee_def, self%scans(scan)%d(i)%N_psd%sigma0_def, self%samprate, ncorr2(j)
            ! else
            !   write(65, '(14(E15.6E3))') n_corr(j,i), s_sub(j,i), mask(j,i), d_prime(j), self%scans(scan)%d(i)%tod(j), self%scans(scan)%d(i)%gain, self%scans(scan)%d(i)%N_psd%alpha, self%scans(scan)%d(i)%N_psd%fknee, self%scans(scan)%d(i)%N_psd%sigma0, self%scans(scan)%d(i)%N_psd%alpha_def, self%scans(scan)%d(i)%N_psd%fknee_def, self%scans(scan)%d(i)%N_psd%sigma0_def, self%samprate, ncorr2(j)
            ! end if
         end do
         close(65)
         !stop
       end if

    end do
    deallocate(dt, dv)
    deallocate(d_prime)
    deallocate(ncorr2)

    call dfftw_destroy_plan(plan_fwd)                                           
    call dfftw_destroy_plan(plan_back)                                          
  
  end subroutine sample_n_corr


  subroutine get_ncorr_sm_cg(handle, d_prime, ncorr, mask, N_psd, samprate, nfft, plan_fwd, plan_back, converged, scan, det, band)
    ! 
    ! Routine implementing the Sherman-Morrison based CG method (from Keihänen et al. 2021 arXiv:2011.06024) to solve the sampling equation for correlated noise in the presence of gaps in the data
    !
    ! Arguments
    ! ---------
    ! handle:  type(planck_rng)
    !          Healpix random number type
    ! d_prime: sp (ntod array) 
    !          Signal subtracted (and gap-filled) raw timestream
    ! mask:    sp (ntodt array)
    !          TOD mask (0 = masked, 1 = included)
    ! N_psd:   class(comm_noise_psd)
    !          Class containing the relevant information about the noise model for this detector and scan
    ! samprate: sp (scalar)
    !          Sample rate of the current detector
    ! nfft:    int (scalar) 
    !          Number of tod points used for fft's. Larger than ntod for padding etc. 
    ! plan_fwd int*8 (scalar)
    !          fftw plan for the forward Fourier transform
    ! plan_back int*8 (scalar)
    !          fftw plan for the backward Fourier transform
    ! scan:    int (scalar)
    !          Scan number
    ! det:     int (scalar)
    !          Detector number
    ! band     string
    !          String denoting the current band (e.g. '030' for LFI)
    !
    ! Returns
    ! -------
    ! n_corr:  sp (ntod array)
    !          TOD-domain correlated noise realization
    ! converged: lgt (scalar)
    !          Flag to signify if the CG method converged or not
    ! 
    implicit none
    type(planck_rng),                intent(inout)  :: handle
    integer*8,                       intent(in)     :: plan_fwd, plan_back
    class(comm_noise_psd),           intent(in)     :: N_psd
    real(sp),                        intent(in)     :: samprate
    integer(i4b),                    intent(in)     :: nfft, scan, det
    logical(lgt),                    intent(out)    :: converged 
    character(len=*),                intent(in)     :: band
    real(sp),          dimension(:), intent(out)    :: ncorr
    real(sp),          dimension(:), intent(in)     :: d_prime, mask

    real(dp)            :: r2, r2new, alp, bet, eps, d2, sigma_bp
    real(sp)            :: freq
    integer(i4b)        :: i, j, k, l, ntod, n, n_iter, nmask
    character(len=1024) :: filename
    real(dp),     allocatable, dimension(:) :: x, b, r, d, Mr, Ad, bp, xp, p, rp
    real(dp),     allocatable, dimension(:) :: invNcorr, invM
    integer(i4b), allocatable, dimension(:) :: u

    n_iter    = 15
    n         = nfft / 2 + 1
    ntod      = size(d_prime, 1)
    eps       = 1.d-5
    converged = .false.
    nmask     = ntod - sum(mask)
    if (nmask == 0) then
       ncorr = d_prime
       return
    end if

    allocate(invM(0:n-1),invNcorr(0:n-1))
    allocate(x(ntod), b(ntod), r(ntod), d(ntod), Mr(ntod), Ad(ntod))
    allocate(u(nmask), bp(nmask), xp(nmask), rp(nmask), p(nmask))

    invNcorr(0) = 0.d0
    invM(0)     = 1.d0
    do l = 1, n-1
       freq        = l*(samprate/2)/(n-1)
       invNcorr(l) = N_psd%sigma0**2 / N_psd%eval_corr(freq)
       invM(l)     = 1.d0 / (1.d0 + invNcorr(l))
    end do

    j = 1
    do i = 1, ntod
       d(i) = rand_gauss(handle)
       r(i) = rand_gauss(handle)
       if (mask(i) == 0.d0) then
          u(j) = i
          j = j + 1
       end if
    end do

    x = d_prime / N_psd%sigma0
    b = (x+d)*mask
    
    call apply_fourier_mat(r, sqrt(invNcorr), Ad, nfft, plan_fwd, plan_back)
    b = b + Ad
    
    x = (x - mean(x)) / 10.d0 + mean(x) ! initial vector for PCG
    
    call apply_fourier_mat(b, invM, Mr, nfft, plan_fwd, plan_back)
    bp = Mr(u)

    if (nmask == 1) then
       sigma_bp = sqrt(variance(real(d_prime,dp)))
    else
       sigma_bp = sqrt(variance(bp))
    end if

    xp   = x(u)  !! not sure exactly how to get a good initial vector
    x    = 0.d0
    x(u) = xp
    call apply_fourier_mat(x, invM, Ad, nfft, plan_fwd, plan_back)
    rp   = bp + Ad(u) - xp
    
    p    = rp
    r2   = sum(rp**2)

    do k = 1, n_iter
       x(u) = p
       call apply_fourier_mat(x, invM, Ad, nfft, plan_fwd, plan_back)
       Ad(u) = p - Ad(u) 
       alp   = r2 / sum(p*Ad(u))
       xp    = xp + alp*p
       rp    = rp - alp * Ad(u)
       r2new = sum(rp**2)
       if (sqrt(abs(r2new)) < eps * sigma_bp * nmask) then  ! average error in each datapoint < eps * sigma_bp
          converged = .true.
          exit
       end if
       bet = r2new / r2
       p  = rp + bet*p
       r2 = r2new
    end do
    x(u) = xp
    call apply_fourier_mat(x + b, invM, Ad, nfft, plan_fwd, plan_back)

    x     = Ad
    ncorr = x * N_psd%sigma0
    
    deallocate(invNcorr, invM)
    deallocate(x, b, r, d, Mr, Ad)
    deallocate(u, bp, xp, p, rp)

  contains

    subroutine apply_fourier_mat(vec, mat, res, nfft, plan_fwd, plan_back)
      real(dp),                  dimension(:),  intent(in)    :: vec
      real(dp),                  dimension(0:), intent(in)    :: mat
      real(dp),                  dimension(:),  intent(out)   :: res
      integer(i4b),                             intent(in)    :: nfft 
      integer*8,                                intent(in)    :: plan_fwd, plan_back

      integer(i4b) :: i, j, k, l, ntod, n, nbuff
      real(sp),     allocatable, dimension(:) :: dt
      complex(spc), allocatable, dimension(:) :: dv

      n    = nfft / 2 + 1
      ntod = size(vec, 1)

      allocate(dt(nfft), dv(0:n-1))
      dt(1:ntod)           = vec(:)
      !dt(2*ntod:ntod+1:-1) = dt(1:ntod)
      
      nbuff = nfft - ntod
      do j=1, nbuff
         dt(ntod+j) = sum(d_prime(ntod-20:ntod)) / 20.0 + (sum(d_prime(1:20)) - sum(d_prime(ntod-20:ntod))) / 20.0 * (j-1) / (nbuff - 1)
      end do


      call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
      dv = dv * mat
      call sfftw_execute_dft_c2r(plan_back, dv, dt)
      dt = dt / nfft
      res(:) = dt(1:ntod)
      deallocate(dt, dv)

    end subroutine apply_fourier_mat
  end subroutine get_ncorr_sm_cg


  subroutine find_d_prime_spikes(self, scan, det, d_prime, pix)
    implicit none
    class(comm_tod),                   intent(in) :: self
    integer(i4b),                      intent(in) :: scan, det
    real(sp),        dimension(:),     intent(in) :: d_prime
    integer(i4b),    dimension(1:,1:), intent(in) :: pix

    integer(i4b) :: i, j, l, k, n, m, nomp, ntod, ndet, err, n_downsamp, n_short
    integer(i4b) :: nfft, nbuff, j_end, j_start, sampnum
    logical(lgt) :: found_spike
    real(dp)     :: nu, power, n_sigma, rms, avg
    character(len=1024) :: filename
    real(dp), allocatable, dimension(:) :: d_downsamp, backup
    
    ntod       = self%scans(scan)%ntod
    n_downsamp = floor(self%samprate)
    n          = ntod - mod(ntod, n_downsamp)
    n_short    = n / n_downsamp

    allocate(backup(n_short), d_downsamp(n_short))

    avg = mean(d_prime * 1.d0)
    do i = 1, n_short
       backup(i) = sum(d_prime((i-1)*n_downsamp+1:i*n_downsamp) - avg) / n_downsamp  
    end do
    
    do i = 1, n_short
       l = max(i-5, 1)
       k = min(i+5, n_short)
       d_downsamp(i) = backup(i) - median(backup(l:k))
    end do
    found_spike = .false.
    rms = sqrt(variance(d_downsamp(:)))
    n_sigma = 10
    do i = 1, n_short
       if (.false. .and. d_downsamp(i) > n_sigma * rms) then
          if (.not. found_spike) then
             write(filename, "(A, I0.3, A, I0.3, 3A)") 'spike_pix_', self%scanid(scan), '_', det, '_',trim(self%freq),'.dat' 
             open(62,file=filename, status='REPLACE')
             found_spike = .true.
          end if
          sampnum = i*n_downsamp - floor(n_downsamp / 2.d0)
          write(62, '(4I7, A)') pix(sampnum,det), sampnum, det, self%scanid(scan), trim(self%freq)
       end if
    end do
    
    if (.false. .and. found_spike) then
       close(62)
       write(filename, "(A, I0.3, A, I0.3, 3A)") 'spike_tod_', self%scanid(scan), '_', det, '_',trim(self%freq),'.dat' 
       open(63,file=filename, status='REPLACE')
       do i = 1, n_short
          sampnum = i*n_downsamp - floor(n_downsamp / 2.d0)
          write(63, '(I7, 2(E15.6E3))') sampnum, d_downsamp(i), rms
       end do
       close(63)
    end if

    deallocate(backup, d_downsamp)
  end subroutine find_d_prime_spikes


  ! Sample noise psd
  subroutine sample_noise_psd(self, tod, handle, scan, mask, s_tot, n_corr, freqmask)
    implicit none
    class(comm_tod),                    intent(inout)  :: self
    real(sp),         dimension(1:,1:), intent(in)     :: tod
    type(planck_rng),                   intent(inout)  :: handle
    integer(i4b),                       intent(in)     :: scan
    real(sp),         dimension(:,:),   intent(in)     :: mask, s_tot, n_corr
    real(sp),         dimension(0:),    intent(in), optional :: freqmask

    integer*8    :: plan_fwd
    integer(i4b) :: i, j, n, nval, n_bins, l, nomp, omp_get_max_threads, err, ntod, n_low, n_high, currdet, currpar
    integer(i4b) :: ndet
    real(sp)     :: f
    real(dp)     :: s, res, log_nu, samprate, gain, dlog_nu, nu, xi_n
    real(dp)     :: alpha, sigma0, fknee, x_in(3), prior_fknee(2), prior_alpha(2), alpha_dpc, fknee_dpc, P_uni(2)
    character(len=1024) :: filename
    real(sp),     allocatable, dimension(:) :: dt, ps
    complex(spc), allocatable, dimension(:) :: dv
    real(sp),     allocatable, dimension(:) :: d_prime
    
    ntod     = self%scans(scan)%ntod
    ndet     = self%ndet
    nomp     = 1 !omp_get_max_threads()
    n        = ntod/2 + 1
    samprate = self%samprate

    ! Sample sigma_0 from pairwise differenced TOD
    do i = 1, ndet
       if (.not. self%scans(scan)%d(i)%accept) cycle
       s    = 0.d0
       nval = 0

       do j = 1, self%scans(scan)%ntod-1
          if (any(mask(j:j+1,i) < 0.5)) cycle
          res = ((tod(j,i)   - self%scans(scan)%d(i)%gain * s_tot(j,i)   - n_corr(j,i))   - &
               & (tod(j+1,i) - self%scans(scan)%d(i)%gain * s_tot(j+1,i) - n_corr(j+1,i)))/sqrt(2.)
          s    = s    + res**2
          nval = nval + 1
       end do
       if (nval > 100) self%scans(scan)%d(i)%N_psd%xi_n(1) = sqrt(s/(nval-1))
    end do

    ! Initialize FFTW
    allocate(dt(ntod), dv(0:n-1), ps(0:n-1))
    call sfftw_init_threads(err)
    call sfftw_plan_with_nthreads(nomp)
    call sfftw_plan_dft_r2c_1d(plan_fwd,  ntod, dt, dv, fftw_estimate + fftw_unaligned)

    ! Sample non-linear spectral parameters
    do i = 1, ndet
       if (.not. self%scans(scan)%d(i)%accept) cycle
       currdet = i

       ! Commpute power spectrum
       n_low  = max(ceiling(self%scans(scan)%d(i)%N_psd%nu_fit(1) * (n-1) / (samprate/2)), 2) ! Never include offset
       n_high =     ceiling(self%scans(scan)%d(i)%N_psd%nu_fit(2) * (n-1) / (samprate/2)) 
       dt     = n_corr(:,i)
       call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
       do l = n_low, n_high
          ps(l) = abs(dv(l)) ** 2 / ntod          
       end do

       ! Perform sampling over all non-linear parameters
       do j = 2, self%scans(scan)%d(i)%N_psd%npar
          P_uni   = self%scans(scan)%d(i)%N_psd%P_uni(j,:)
          if (self%scans(scan)%d(i)%N_psd%P_active(j,2) <= 0.d0 .or. P_uni(2) == P_uni(1)) cycle

          currpar = j
          xi_n    = self%scans(scan)%d(i)%N_psd%xi_n(j)
          x_in(1) = max(xi_n - 0.5 * abs(xi_n), P_uni(1))
          x_in(3) = min(xi_n + 0.5 * abs(xi_n), P_uni(2))
          x_in(2) = 0.5 * (x_in(1) + x_in(3))

          xi_n = sample_InvSamp(handle, x_in, lnL_xi_n, P_uni)
          xi_n = min(max(xi_n,self%scans(scan)%d(i)%N_psd%P_uni(j,1)), self%scans(scan)%d(i)%N_psd%P_uni(j,2))
          self%scans(scan)%d(i)%N_psd%xi_n(j) = xi_n
       end do
    end do
    deallocate(dt, dv)
    deallocate(ps)
    
    call sfftw_destroy_plan(plan_fwd)
    
  contains

    function lnL_xi_n(x) 
      use healpix_types
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: lnL_xi_n

      real(dp) :: sconst, s, N_corr, mu, sigma
      real(sp) :: f, tmp

      if (x < P_uni(1) .or. x > P_uni(2)) then
         lnL_xi_n = -1.d30
         return
      end if

      ! Update xi_n with new proposal
      tmp = self%scans(scan)%d(i)%N_psd%xi_n(currpar) 
      self%scans(scan)%d(i)%N_psd%xi_n(currpar) = x

      ! Add likelihood term
      lnL_xi_n = 0.d0
      do l = n_low, n_high
         if (present(freqmask)) then
            if (freqmask(l) == 0.) cycle
         end if
         f         = l*(samprate/2)/(n-1)
         N_corr    = self%scans(scan)%d(currdet)%N_psd%eval_corr(f)
         lnL_xi_n  = lnL_xi_n - (ps(l) / N_corr + log(N_corr))
      end do

      ! Add prior
      mu    = self%scans(scan)%d(currdet)%N_psd%P_active(currpar,1)
      sigma = self%scans(scan)%d(currdet)%N_psd%P_active(currpar,2)
      if (self%scans(scan)%d(currdet)%N_psd%P_lognorm(currpar)) then 
         ! Log-normal prior
         lnL_xi_n = lnL_xi_n - 0.5d0 * (log(x) - log(mu))**2 / (sigma * log(10.d0))**2 - log(x)
      else
         ! Gaussian prior
         lnL_xi_n = lnL_xi_n - 0.5d0 * (x - mu)**2 / sigma**2
      end if

      ! Revert xi_n with old value
      self%scans(scan)%d(i)%N_psd%xi_n(currpar) = tmp
      
    end function lnL_xi_n
       
  end subroutine sample_noise_psd


  ! Routine for multiplying a set of timestreams with inverse noise covariance 
  ! matrix. inp and res have dimensions (ntime,ndet,ninput), where ninput is 
  ! the number of input timestreams (e.g. 2 if you input s_tot and s_orb). 
  ! Here inp and res are assumed to be already allocated. 
  subroutine multiply_inv_N(tod, scan, buffer, sampfreq, pow)
    implicit none
    class(comm_tod),                     intent(in)     :: tod
    integer(i4b),                        intent(in)     :: scan
    real(sp),          dimension(:,:),   intent(inout)  :: buffer !input/output
    real(dp),                            intent(in), optional :: sampfreq, pow
    integer(i4b) :: i, j, l, n, m, nomp, ntod, ndet, err, omp_get_max_threads
    integer*8    :: plan_fwd, plan_back
    real(sp)     :: sigma_0, alpha, nu_knee,  samprate, noise, signal
    real(sp)     :: nu, pow_
    real(sp),     allocatable, dimension(:,:) :: dt
    complex(spc), allocatable, dimension(:,:) :: dv
    
    ntod     = size(buffer, 1)
    ndet     = size(buffer, 2)
    nomp     = omp_get_max_threads()
    m        = count(tod%scans(scan)%d%accept)
    n        = ntod + 1
    pow_     = 1.d0; if (present(pow)) pow_ = pow
    samprate = real(tod%samprate,sp); if (present(sampfreq)) samprate = real(sampfreq,sp)
    
    allocate(dt(2*ntod,m), dv(0:n-1,m))
    call sfftw_plan_many_dft_r2c(plan_fwd, 1, 2*ntod, m, dt, &
         & 2*ntod, 1, 2*ntod, dv, n, 1, n, fftw_patient)
    call sfftw_plan_many_dft_c2r(plan_back, 1, 2*ntod, m, dv, &
         & n, 1, n, dt, 2*ntod, 1, 2*ntod, fftw_patient)
    
    j = 0
    do i = 1, ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle
       j = j+1
       dt(1:ntod,j)           = buffer(:,i)
       dt(2*ntod:ntod+1:-1,j) = dt(1:ntod,j)
    end do

    call sfftw_execute_dft_r2c(plan_fwd, dt, dv)

    j = 0
    do i = 1, ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle
       j = j+1
       noise    = tod%scans(scan)%d(i)%N_psd%sigma0**2 * samprate / tod%samprate
       
       if (pow_ >= 0.d0) dv(0,j) = 0.d0   ! If pow < 0, leave offset as is
       do l = 1, n-1                                                      
          nu      = l*(samprate/2)/(n-1)
          signal  = tod%scans(scan)%d(i)%N_psd%eval_corr(nu) * samprate / tod%samprate
          dv(l,j) = dv(l,j) * 1.0/(noise + signal)**pow_
       end do
    end do
       
    call sfftw_execute_dft_c2r(plan_back, dv, dt)
    dt = dt / (2*ntod)

    j = 0
    do i = 1, ndet
       if (.not. tod%scans(scan)%d(i)%accept) then
          buffer(:,i)  = 0.d0
       else
          j           = j+1
          buffer(:,i) = dt(1:ntod,j) 
       end if
    end do
    call sfftw_destroy_plan(plan_fwd)                                           
    call sfftw_destroy_plan(plan_back)                                          
    deallocate(dt, dv)

  end subroutine multiply_inv_N

end module comm_tod_noise_mod
