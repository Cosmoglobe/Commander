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
  use InvSamp_mod
  implicit none



contains

 ! Compute correlated noise term, n_corr from eq:
  ! ((N_c^-1 + N_wn^-1) n_corr = d_prime + w1 * sqrt(N_wn) + w2 * sqrt(N_c) 
  subroutine sample_n_corr(self, handle, scan, mask, s_sub, n_corr, pix, dospike)
    implicit none
    class(comm_tod),               intent(in)     :: self
    type(planck_rng),                  intent(inout)  :: handle
    integer(i4b),                      intent(in)     :: scan
    integer(i4b),   dimension(1:,1:),  intent(in)     :: pix
    real(sp),          dimension(:,:), intent(in)     :: mask, s_sub
    logical(lgt),            intent(in), optional     :: dospike
    real(sp),          dimension(:,:), intent(out)    :: n_corr
    integer(i4b) :: i, j, l, k, n, m, nomp, ntod, ndet, err, omp_get_max_threads
    integer(i4b) :: nfft, nbuff, j_end, j_start
    integer*8    :: plan_fwd, plan_back
    logical(lgt) :: init_masked_region, end_masked_region, pcg_converged
    real(sp)     :: sigma_0, alpha, nu_knee,  samprate, gain, mean, N_wn, N_c
    real(dp)     :: nu, power, fft_norm
    ! real(sp),     allocatable, dimension(:) :: dt
    ! complex(spc), allocatable, dimension(:) :: dv
    real(dp),     allocatable, dimension(:) :: dt
    complex(dpc), allocatable, dimension(:) :: dv
    real(sp),     allocatable, dimension(:) :: d_prime, ncorr2
    character(len=1024) :: filename

    ntod = self%scans(scan)%ntod
    ndet = self%ndet
    nomp = omp_get_max_threads()
    
    nfft = 2 * ntod
    !nfft = get_closest_fft_magic_number(ceiling(ntod * 1.05d0))
    
    n = nfft / 2 + 1

    call dfftw_init_threads(err)
    call dfftw_plan_with_nthreads(nomp)

    allocate(dt(nfft), dv(0:n-1))
    call dfftw_plan_dft_r2c_1d(plan_fwd,  nfft, dt, dv, fftw_estimate + fftw_unaligned)
    call dfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)

    !$OMP PARALLEL PRIVATE(i,j,l,k,dt,dv,nu,sigma_0,alpha,nu_knee,d_prime,init_masked_region,end_masked_region)
    allocate(d_prime(ntod))
    allocate(ncorr2(ntod))
    allocate(dt(nfft), dv(0:n-1))

    !$OMP DO SCHEDULE(guided)
    do i = 1, ndet
       if (.not. self%scans(scan)%d(i)%accept) cycle
       gain = self%scans(scan)%d(i)%gain  ! Gain in V / K
       d_prime(:) = self%scans(scan)%d(i)%tod(:) - self%scans(scan)%d(i)%baseline &
                &- S_sub(:,i) * gain

       sigma_0 = abs(self%scans(scan)%d(i)%sigma0)
       !write(*,*) "rms:", scan, sigma_0, sqrt(sum(d_prime**2)/size(d_prime))
       ! Fill gaps in data 
       init_masked_region = .true.
       end_masked_region = .false.
       do j = 1,ntod
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

       if (self%first_call .and. not(present(dospike))) then
          call find_d_prime_spikes(self, scan, i, d_prime, pix)
       end if

       samprate = self%samprate
       alpha    = self%scans(scan)%d(i)%alpha
       nu_knee  = self%scans(scan)%d(i)%fknee
       N_wn = sigma_0 ** 2  ! white noise power spectrum

       !call get_ncorr_pcg(handle, d_prime, ncorr2, mask(:,i), alpha, nu_knee, N_wn, samprate, nfft, plan_fwd, plan_back)
       !n_corr(:, i) = ncorr2(:)
       
       pcg_converged = .false.
!  subroutine get_ncorr_pcg(handle, d_prime, ncorr, mask, alpha, fknee, wn, samprate, nfft, plan_fwd, plan_back, converged, scan, det, freq)
       !!!! add choice between PCG and regular ncorr here
       if (.true.) then !(self%scanid(scan) == 2112) .and. (i == 1)) then
          call get_ncorr_sm_cg(handle, d_prime, ncorr2, mask(:,i), alpha, nu_knee, N_wn, samprate, nfft, plan_fwd, plan_back, pcg_converged, self%scanid(scan), i, trim(self%freq))
          n_corr(:, i) = ncorr2(:)
       end if

       if (.not. pcg_converged) then
          ! Preparing for fft
          dt(1:ntod)           = d_prime(:)
          dt(2*ntod:ntod+1:-1) = dt(1:ntod)
          ! nbuff = nfft - ntod
          ! do j=1, nbuff
          !    dt(ntod+j) = d_prime(ntod) + (d_prime(1) - d_prime(ntod)) * (j-1) / (nbuff - 1)
          ! end do


          call dfftw_execute_dft_r2c(plan_fwd, dt, dv)


          fft_norm = sqrt(1.d0 * nfft)  ! used when adding fluctuation terms to Fourier coeffs (depends on Fourier convention)

          if (trim(self%operation) == "sample") then
             dv(0)    = dv(0) + fft_norm * sqrt(N_wn) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0)
          end if


          do l = 1, n-1                                                      
             nu = l*(samprate/2)/(n-1)

             N_c = N_wn * (nu/(nu_knee))**(alpha)  ! correlated noise power spectrum

             if (N_c == 0) N_c = 1
             if (abs(nu-1.d0/60.d0)*60.d0 < 0.05d0) then
                dv(l) = fft_norm * sqrt(N_c) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0) ! Dont include scan frequency; replace with better solution
             end if


             if (trim(self%operation) == "sample") then
                dv(l) = (dv(l) + fft_norm * ( &
                     sqrt(N_wn) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0) &
                     + N_wn * sqrt(1.0 / N_c) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0) &
                     )) * 1.d0/(1.d0 + N_wn / N_c)
             else
                dv(l) = dv(l) * 1.0/(1.0 + N_wn/N_c)
             end if
          end do
          call dfftw_execute_dft_c2r(plan_back, dv, dt)
          dt          = dt / nfft
          n_corr(:,i) = dt(1:ntod) 
       end if

       if (mod(self%scanid(scan),100) == 1) then
          write(filename, "(A, I0.3, A, I0.3, 3A)") 'ncorr_tods/ncorr_times', self%scanid(scan), '_', i, '_',trim(self%freq),'_final_hundred.dat' 
          open(65,file=trim(filename),status='REPLACE')
          do j = 1, ntod
             write(65, '(14(E15.6E3))') n_corr(j,i), s_sub(j,i), mask(j,i), d_prime(j), self%scans(scan)%d(i)%tod(j), self%scans(scan)%d(i)%gain, self%scans(scan)%d(i)%alpha, self%scans(scan)%d(i)%fknee, self%scans(scan)%d(i)%sigma0, self%scans(scan)%d(i)%alpha_def, self%scans(scan)%d(i)%fknee_def, self%scans(scan)%d(i)%sigma0_def, self%samprate, ncorr2(j)
          end do
          close(65)
          !stop
       end if

    end do
    !$OMP END DO                                                          
    deallocate(dt, dv)
    deallocate(d_prime)
    deallocate(ncorr2)
    !deallocate(diff)
    !$OMP END PARALLEL
    
    

    call dfftw_destroy_plan(plan_fwd)                                           
    call dfftw_destroy_plan(plan_back)                                          
  
  end subroutine sample_n_corr



  ! subroutine test_fft(handle, d_prime, ncorr, mask, alpha, fknee, wn, samprate, nfft, plan_fwd, plan_back, converged, scan, det, band)
  !   implicit none
  !   type(planck_rng),    intent(inout)  :: handle
  !   integer*8,  intent(in)   :: plan_fwd, plan_back
  !   real(sp),   intent(in)   :: alpha, fknee, wn, samprate
  !   integer(i4b), intent(in) :: nfft, scan, det
  !   logical(lgt), intent(out) :: converged 
  !   character(len=*), intent(in) :: band
  !   real(sp),          dimension(:), intent(out)   :: ncorr
  !   real(sp),          dimension(:), intent(in)    :: d_prime, mask
  !   real(dp),    allocatable,      dimension(:) :: x, b, r, d, Mr, Ad, u, bp, xp, p, rp
  !   real(dp),    allocatable,      dimension(:) :: invNcorr, invM
  !   real(dp)           :: r2, r2new, alp, bet, eps, freq, d2
  !   real(dp),     allocatable, dimension(:) :: dt
  !   complex(dpc), allocatable, dimension(:) :: dv
  !   ! real(sp),     allocatable, dimension(:) :: dt
  !   ! complex(spc), allocatable, dimension(:) :: dv
  !   integer(i4b) :: i, j, k, l, ntod, n, n_iter, nmask
  !   character(len=1024) :: filename

  !   n_iter = 10
  !   n = nfft / 2 + 1
  !   ntod = size(d_prime, 1)
  !   nmask = ntod - sum(mask)
  !   write(*,*) ntod, nfft, n
  !   eps = 0.1d0

  !   converged = .false.

  !   allocate(dt(nfft), dv(0:n-1))
  !   allocate(invM(0:n-1),invNcorr(0:n-1))
  !   allocate(x(ntod), b(ntod), r(ntod), d(ntod), Mr(ntod), Ad(ntod))
  !   allocate(u(nmask), bp(nmask), xp(nmask), rp(nmask), p(nmask))

  !   do l=0, n-1
  !      freq = l*(samprate/2)/(n-1)
  !      invNcorr(l) = 1.d0 * (freq / fknee) ** (-alpha)
  !      invM(l) = 1.d0 / ( 1.d0 * (1.d0 + (freq / fknee) ** (-alpha)))
  !   end do
  !   ! b = d_prime * mask / wn + w1 * mask / np.sqrt(wn) + apply_sqrt_S_inv(w2, f, mask, alpha, fknee, wn)
  !   j = 1
  !   do i = 1, ntod
  !      r(i) = rand_gauss(handle)
  !   end do
  !   write(*,*) invM(0)
  !   call apply_fourier_mat(r, invM, Mr, dt, dv, nfft, plan_fwd, plan_back)
    
  !   call apply_fourier_mat(Mr, 1.d0 / invM, d, dt, dv, nfft, plan_fwd, plan_back)
    
  !   x(:) = r(:) - d(:)
  !   if (.true.) then !.not. converged) then
  !      write(filename, "(A)") 'fft_test.dat' 
  !      open(63,file=filename, status='REPLACE')
  !      do i = 1, ntod
  !         write(63, '(14(E21.11E3))') d_prime(i), b(i), x(i), r(i), d(i), Ad(i), Mr(i), r2, alp, bet, alpha, fknee, wn, mask(i)
  !      end do
  !      close(63)
  !   end if
    
  !   deallocate(dt, dv)
  !   deallocate(invNcorr, invM)
  !   deallocate(x, b, r, d, Mr, Ad)
  !   deallocate(u, bp, xp, p, rp)
  ! contains

  !   subroutine apply_fourier_mat(vec, mat, res, dt, dv, nfft, plan_fwd, plan_back)
  !     integer*8,  intent(in)   :: plan_fwd, plan_back
  !     !real(dp),   intent(in)   :: alpha, fknee, wn, samprate
  !     real(dp),  dimension(0:), intent(in)    :: mat
  !     real(dp),  dimension(:), intent(in)    :: vec
  !     real(dp),  dimension(:), intent(out)   :: res
  !     !real(sp),  dimension(:)  :: apply_fourier_mat
  !     integer(i4b), intent(in) :: nfft 
  !     real(dp),     dimension(:), intent(inout) :: dt
  !     complex(dpc), dimension(0:), intent(inout) :: dv
  !     ! real(sp),     allocatable, dimension(:), intent(inout) :: dt
  !     ! complex(spc), allocatable, dimension(:), intent(inout) :: dv
  !     integer(i4b) :: i, j, k, l, ntod, n
  !     n = nfft / 2 + 1
  !     ntod = size(vec, 1)
  !     dt(1:ntod)           = vec(:)
  !     dt(2*ntod:ntod+1:-1) = dt(1:ntod)
      
  !     call dfftw_execute_dft_r2c(plan_fwd, dt, dv)
  !     do l = 0, n-1
  !        dv(l) = dv(l) * mat(l)
  !     end do
  !     call dfftw_execute_dft_c2r(plan_back, dv, dt)
  !     dt = dt / nfft
  !     res(:) = dt(1:ntod)
      
  !   end subroutine apply_fourier_mat
  ! end subroutine test_fft




  subroutine get_ncorr_sm_cg(handle, d_prime, ncorr, mask, alpha, fknee, wn, samprate, nfft, plan_fwd, plan_back, converged, scan, det, band)
    implicit none
    type(planck_rng),    intent(inout)  :: handle
    integer*8,  intent(in)   :: plan_fwd, plan_back
    real(sp),   intent(in)   :: alpha, fknee, wn, samprate
    integer(i4b), intent(in) :: nfft, scan, det
    logical(lgt), intent(out) :: converged 
    character(len=*), intent(in) :: band
    real(sp),          dimension(:), intent(out)   :: ncorr
    real(sp),          dimension(:), intent(in)    :: d_prime, mask
    real(dp),    allocatable,      dimension(:) :: x, b, r, d, Mr, Ad, bp, xp, p, rp
    real(dp),    allocatable,      dimension(:) :: invNcorr, invM
    integer(i4b),    allocatable,  dimension(:) :: u
    real(dp)           :: r2, r2new, alp, bet, eps, freq, d2, sigma_bp
    real(dp),     allocatable, dimension(:) :: dt
    complex(dpc), allocatable, dimension(:) :: dv
    ! real(sp),     allocatable, dimension(:) :: dt
    ! complex(spc), allocatable, dimension(:) :: dv
    integer(i4b) :: i, j, k, l, ntod, n, n_iter, nmask
    character(len=1024) :: filename

    n_iter = 15
    n = nfft / 2 + 1
    ntod = size(d_prime, 1)
    nmask = ntod - sum(mask)
    
    eps = 1.d-5

    converged = .false.

    allocate(dt(nfft), dv(0:n-1))
    allocate(invM(0:n-1),invNcorr(0:n-1))
    allocate(x(ntod), b(ntod), r(ntod), d(ntod), Mr(ntod), Ad(ntod))
    allocate(u(nmask), bp(nmask), xp(nmask), rp(nmask), p(nmask))

    do l=0, n-1
       freq = l*(samprate/2)/(n-1)
       invNcorr(l) = 1.d0 * (freq / fknee) ** (-alpha)
       invM(l) = 1.d0 / ( 1.d0 * (1.d0 + (freq / fknee) ** (-alpha)))
    end do
    ! b = d_prime * mask / wn + w1 * mask / np.sqrt(wn) + apply_sqrt_S_inv(w2, f, mask, alpha, fknee, wn)
    j = 1
    do i = 1, ntod
       d(i) = rand_gauss(handle)
       r(i) = rand_gauss(handle)
       if (mask(i) == 0.d0) then
          u(j) = i
          j = j + 1
       end if
    end do

    x(:) = d_prime(:) / sqrt(wn) 
    b(:) = x * mask(:) + d(:) * mask(:) 
    
    
    call apply_fourier_mat(r, sqrt(invNcorr), Ad, dt, dv, nfft, plan_fwd, plan_back)
    b(:) = b(:) + Ad(:)
    
    x(:) = (x(:) - mean(x)) / 10.d0 + mean(x) ! initial vector for PCG
    
    call apply_fourier_mat(b, invM, Mr, dt, dv, nfft, plan_fwd, plan_back)
    bp(:) = Mr(u)

    sigma_bp = sqrt(variance(bp))

    xp(:) = x(u)  !! not sure exactly how to get a good initial vector
    x(:) = 0.d0
    x(u) = xp(:)
    call apply_fourier_mat(x, invM, Ad, dt, dv, nfft, plan_fwd, plan_back)
    rp(:) = bp(:) + Ad(u) - xp(:)
    
    p(:) = rp(:)
    r2 = sum(rp(:) ** 2)

    do k = 1, n_iter
       x(u) = p(:)
       call apply_fourier_mat(x, invM, Ad, dt, dv, nfft, plan_fwd, plan_back)
       Ad(u) = p(:) - Ad(u) 
       !d2 = sum(d(:) ** 2 * mask(:) + Ad(:) ** 2)

       alp = r2 / sum(p(:) * Ad(u))
       xp(:) = xp(:) + alp * p(:)
    
       rp(:) = rp(:) - alp * Ad(u)
    
       r2new = sum(rp(:) ** 2)
       if (sqrt(abs(r2new)) < eps * sigma_bp * nmask) then  ! average error in each datapoint < eps * sigma_bp
          converged = .true.
          !write(*,*) 'converged after ', k, ' iterations ', scan, det, trim(band)
          exit
       end if
       bet = r2new / r2
       !write(*,*) 'iteration ', k, r2, r2new, scan, det, trim(band)
       p(:) = rp(:) + bet * p(:)
       r2 = r2new
    end do
    x(u) = xp(:)
    call apply_fourier_mat(x + b, invM, Ad, dt, dv, nfft, plan_fwd, plan_back)

    ! !! extra
    ! call apply_fourier_mat(Ad, invNcorr, Mr, dt, dv, nfft, plan_fwd, plan_back)
    ! Mr(:) =  Mr(:) + Ad(:) * mask(:)  !!! this is Ax 
    ! d(:) = Mr(:) - b(:)

    ! call apply_fourier_mat(Ad, 1.d0 / invM, r, dt, dv, nfft, plan_fwd, plan_back)
    ! r(:) = r(:) - x(:) - b(:) 

    !if (.false. .and. .not. converged) then
    !   write(filename, "(A, I0.3, A, I0.3, 3A)") 'ms_cg_tod_', scan, '_', det, '_',trim(band),'.dat' 
    !   open(63,file=filename, status='REPLACE')
    !   do i = 1, ntod
    !      write(63, '(14(E21.11E3))') d_prime(i), b(i), x(i), r(i), d(i), Ad(i), Mr(i), r2, alp, bet, alpha, fknee, wn, mask(i)
    !   end do
    !   close(63)
    !end if
    x(:) = Ad(:)
    ncorr(:) = x(:) * sqrt(wn)
    
    deallocate(dt, dv)
    deallocate(invNcorr, invM)
    deallocate(x, b, r, d, Mr, Ad)
    deallocate(u, bp, xp, p, rp)
  contains

    subroutine apply_fourier_mat(vec, mat, res, dt, dv, nfft, plan_fwd, plan_back)
      integer*8,  intent(in)   :: plan_fwd, plan_back
      !real(dp),   intent(in)   :: alpha, fknee, wn, samprate
      real(dp),  dimension(:), intent(in)    :: vec
      real(dp),  dimension(0:), intent(in)    :: mat
      real(dp),  dimension(:), intent(out)   :: res
      !real(sp),  dimension(:)  :: apply_fourier_mat
      integer(i4b), intent(in) :: nfft 
      real(dp),      dimension(:), intent(inout) :: dt
      complex(dpc), dimension(0:), intent(inout) :: dv
      ! real(sp),     allocatable, dimension(:), intent(inout) :: dt
      ! complex(spc), allocatable, dimension(:), intent(inout) :: dv
      integer(i4b) :: i, j, k, l, ntod, n
      n = nfft / 2 + 1
      ntod = size(vec, 1)
      dt(1:ntod)           = vec(:)
      dt(2*ntod:ntod+1:-1) = dt(1:ntod)

      call dfftw_execute_dft_r2c(plan_fwd, dt, dv)
      do l = 0, n-1
         dv(l) = dv(l) * mat(l)
      end do
      call dfftw_execute_dft_c2r(plan_back, dv, dt)
      dt = dt / nfft
      res(:) = dt(1:ntod)
      
    end subroutine apply_fourier_mat
  end subroutine get_ncorr_sm_cg

!!!!!!!!!!!!!!!!!!!!! OLD PCG method below !!!!!!!!!!!!!!

  subroutine get_ncorr_pcg(handle, d_prime, ncorr, mask, alpha, fknee, wn, samprate, nfft, plan_fwd, plan_back, converged, scan, det, band)
    implicit none
    type(planck_rng),    intent(inout)  :: handle
    integer*8,  intent(in)   :: plan_fwd, plan_back
    real(sp),   intent(in)   :: alpha, fknee, wn, samprate
    integer(i4b), intent(in) :: nfft, scan, det
    logical(lgt), intent(out) :: converged 
    character(len=*), intent(in) :: band
    real(sp),          dimension(:), intent(out)   :: ncorr
    real(sp),          dimension(:), intent(in)    :: d_prime, mask
    real(dp),    allocatable,      dimension(:) :: x, b, r, d, Mr, Ad
    real(dp),    allocatable,      dimension(:) :: invNcorr, invM
    real(dp)           :: r2, r2new, alp, bet, eps, freq, d2
    real(dp),     allocatable, dimension(:) :: dt
    complex(dpc), allocatable, dimension(:) :: dv
    ! real(sp),     allocatable, dimension(:) :: dt
    ! complex(spc), allocatable, dimension(:) :: dv
    integer(i4b) :: i, j, k, l, ntod, n, n_iter
    character(len=1024) :: filename

    n_iter = 10
    n = nfft / 2 + 1
    ntod = size(d_prime, 1)
      
    eps = 0.1d0

    converged = .false.

    allocate(dt(nfft), dv(0:n-1))
    allocate(invNcorr(0:n-1),invM(0:n-1))
    allocate(x(ntod), b(ntod), r(ntod), d(ntod), Mr(ntod), Ad(ntod))

    do l=0, n-1
       freq = l*(samprate/2)/(n-1)
       invNcorr(l) = 1.d0 * (freq / fknee) ** (-alpha)
       invM(l) = 1.d0 / ( 1.d0 * (1.d0 + (freq / fknee) ** (-alpha)))
    end do
    ! b = d_prime * mask / wn + w1 * mask / np.sqrt(wn) + apply_sqrt_S_inv(w2, f, mask, alpha, fknee, wn)
    do i = 1, ntod
       d(i) = rand_gauss(handle)
       r(i) = rand_gauss(handle)
    end do
    x(:) = d_prime(:) / sqrt(wn) 
    b(:) = x * mask(:) + d(:) * mask(:) 
    
    
    call apply_fourier_mat(r, sqrt(invNcorr), Ad, dt, dv, nfft, plan_fwd, plan_back)
    b(:) = b(:) + Ad(:)
    
    x(:) = (x(:) - mean(x)) / 10.d0 + mean(x) ! initial vector for PCG
    call apply_fourier_mat(x, invNcorr, Ad, dt, dv, nfft, plan_fwd, plan_back)
    r(:) = b(:) - Ad(:) - x(:) * mask(:) 
    
    call apply_fourier_mat(r, invM, Mr, dt, dv, nfft, plan_fwd, plan_back)
    d(:) = Mr(:)
    r2 = sum(r(:) * Mr(:))
    
    do k = 1, n_iter
       call apply_fourier_mat(d, invNcorr, Ad, dt, dv, nfft, plan_fwd, plan_back)
       Ad(:) = Ad(:) + d(:) * mask(:) 
       !d2 = sum(d(:) ** 2 * mask(:) + Ad(:) ** 2)

       alp = r2 / sum(d(:) * Ad(:))
       x(:) = x(:) + alp * d(:)
    
       r(:) = r(:) - alp * Ad(:)
    
       call apply_fourier_mat(r, invM, Mr, dt, dv, nfft, plan_fwd, plan_back)
       r2new = sum(Mr(:) * r(:))
       if (sqrt(abs(r2new)) < eps) then
          converged = .true.
!          write(*,*) 'converged after ', k, ' iterations ', scan, det, trim(band)
          exit
       end if
       bet = r2new / r2
!       write(*,*) 'iteration ', k, r2, r2new, scan, det, trim(band)
       d(:) = Mr(:) + bet * d(:)
       r2 = r2new
    end do
    
    if (.not. converged) then
       write(filename, "(A, I0.3, A, I0.3, 3A)") 'pcg_tod_', scan, '_', det, '_',trim(band),'.dat' 
       open(63,file=filename, status='REPLACE')
       do i = 1, ntod
          write(63, '(13(E15.6E3))') d_prime(i), b(i), x(i), r(i), d(i), Ad(i), Mr(i), r2, alp, bet, alpha, fknee, wn 
       end do
       close(63)
    end if
    
    ncorr(:) = x(:) * sqrt(wn)
    
    deallocate(dt, dv)
    deallocate(invNcorr, invM)
    deallocate(x, b, r, d, Mr, Ad)
  contains

    subroutine apply_fourier_mat(vec, mat, res, dt, dv, nfft, plan_fwd, plan_back)
      integer*8,  intent(in)   :: plan_fwd, plan_back
      !real(dp),   intent(in)   :: alpha, fknee, wn, samprate
      real(dp),  dimension(:), intent(in)    :: mat, vec
      real(dp),  dimension(:), intent(out)   :: res
      !real(sp),  dimension(:)  :: apply_fourier_mat
      integer(i4b), intent(in) :: nfft 
      real(dp),     allocatable, dimension(:), intent(inout) :: dt
      complex(dpc), allocatable, dimension(:), intent(inout) :: dv
      ! real(sp),     allocatable, dimension(:), intent(inout) :: dt
      ! complex(spc), allocatable, dimension(:), intent(inout) :: dv
      integer(i4b) :: i, j, k, l, ntod, n
      n = nfft / 2 + 1
      ntod = size(vec, 1)
      dt(1:ntod)           = vec(:)
      dt(2*ntod:ntod+1:-1) = dt(1:ntod)

      call dfftw_execute_dft_r2c(plan_fwd, dt, dv)
      do l = 0, n-1
         dv(l) = dv(l) * mat(l)
      end do
      call dfftw_execute_dft_c2r(plan_back, dv, dt)
      dt = dt / nfft
      res(:) = dt(1:ntod)
      
    end subroutine apply_fourier_mat
  end subroutine get_ncorr_pcg


  ! subroutine get_ncorr_pcg(handle, d_prime, ncorr, mask, alpha, fknee, wn, samprate, nfft, plan_fwd, plan_back)
  !   implicit none
  !   type(planck_rng),    intent(inout)  :: handle
  !   integer*8,  intent(in)   :: plan_fwd, plan_back
  !   real(sp),   intent(in)   :: alpha, fknee, wn, samprate
  !   integer(i4b), intent(in) :: nfft 
  !   real(sp),          dimension(:), intent(out)   :: ncorr
  !   real(sp),          dimension(:), intent(in)    :: d_prime, mask
  !   real(dp),    allocatable,      dimension(:) :: x, b, r, d, Mr, Ad
  !   real(dp),    allocatable,      dimension(:) :: invNcorr, invM
  !   real(dp)           :: r2, r2new, alp, bet, eps, freq
  !   ! real(dp),     allocatable, dimension(:) :: dt
  !   ! complex(dpc), allocatable, dimension(:) :: dv
  !   real(sp),     allocatable, dimension(:) :: dt
  !   complex(spc), allocatable, dimension(:) :: dv
  !   integer(i4b) :: i, j, k, l, ntod, n, n_iter
  !   n_iter = 5
  !   n = nfft / 2 + 1
  !   ntod = size(d_prime, 1)
      
  !   eps = 0.01d0
  !   allocate(dt(nfft), dv(0:n-1))
  !   allocate(invNcorr(0:n-1),invM(0:n-1))
  !   allocate(x(ntod), b(ntod), r(ntod), d(ntod), Mr(ntod), Ad(ntod))


  !   do l=0, n-1
  !      freq = l*(samprate/2)/(n-1)
  !      invNcorr(l) = 1.d0 / wn * (freq / fknee) ** (-alpha)
  !      invM(l) = 1.d0 / ( 1.d0 / wn * (1.d0 + (freq / fknee) ** (-alpha)))
  !   end do
  !   ! b = d_prime * mask / wn + w1 * mask / np.sqrt(wn) + apply_sqrt_S_inv(w2, f, mask, alpha, fknee, wn)
  !   do i = 1, ntod
  !      d(i) = rand_gauss(handle)
  !      r(i) = rand_gauss(handle)
  !   end do
  !   b(:) = d_prime(:) * mask(:) / wn + d(:) * mask(:) / sqrt(wn)
    
  !   open(60,file='test_pcg.dat',status='REPLACE')
    
  !   write(60, '(*(E15.6E3))') d(:)
  !   write(60, '(*(E15.6E3))') r(:)

    
    
  !   call apply_fourier_mat(r, sqrt(invNcorr), Ad, dt, dv, nfft, plan_fwd, plan_back)
    
  !   write(60, '(*(E15.6E3))') Ad(:)
  !   b(:) = b(:) + Ad(:)
  !   write(60, '(*(E15.6E3))') b(:)
      
  !   !x(:) = d_prime(:)  ! initial vector for PCG

  !   call apply_fourier_mat(x, invNcorr, Ad, dt, dv, nfft, plan_fwd, plan_back)
  !   r(:) = b(:) - Ad(:) - x(:) * mask(:) / wn
  !   write(60, '(*(E15.6E3))') r(:)
  !   write(60, '(*(E15.6E3))') Ad(:)
    
  !   call apply_fourier_mat(r, invM, Mr, dt, dv, nfft, plan_fwd, plan_back)
  !   write(60, '(*(E15.6E3))') Mr(:)
  !   d(:) = Mr(:)
  !   r2 = sum(r(:) * Mr(:))
    
  !   do k = 1, n_iter
  !      call apply_fourier_mat(d, invNcorr, Ad, dt, dv, nfft, plan_fwd, plan_back)
  !      Ad(:) = Ad(:) + d(:) * mask(:) / wn
  !      write(60, '(*(E15.6E3))') Ad(:)
    
  !      if (.true.) then !mod(self%scanid(scan),100) == 1) then
  !         !write(filename, "(A, I0.3, A, I0.3, 3A)") 'ncorr_pcg.dat' 
  !         open(65,file='ncorr_pcg.dat',status='REPLACE')
  !         do j = 1, ntod
  !            write(65, '(9(E15.6E3))') d_prime(j), mask(j), b(j), x(j), r(j), d(j), Mr(j), Ad(j), r2
  !         end do
  !         close(65)
  !         open(65,file='ncorr_invmats.dat',status='REPLACE')
  !         do j = 0, n-1
  !            freq = j*(samprate/2)/(n-1)
  !            write(65, '(6(E15.6E3))') freq, invNcorr(j), invM(j), alpha, fknee, wn
  !         end do
  !         close(65)
  !         !stop
  !      end if


  !      alp = r2 / sum(d(:) * Ad(:))
  !      x(:) = x(:) + alp * d(:)
  !      write(60, '(*(E15.6E3))') x(:)
    
  !      ! can add if-test here to skip the final fourier transform
  !      r(:) = r(:) - alp * Ad(:)
  !      write(60, '(*(E15.6E3))') r(:)
    
  !      call apply_fourier_mat(r, invM, Mr, dt, dv, nfft, plan_fwd, plan_back)
  !      write(60, '(*(E15.6E3))') Mr(:)
  !      r2new = sum(r(:) * Mr(:))
  !      !! add convergence test here
  !      bet = r2new / r2
  !      d(:) = Mr(:) + bet * d(:)
  !      write(60, '(*(E15.6E3))') d(:)
  !      write(*,*) alp, bet, r2, r2new
  !      r2 = r2new
  !   end do
    
  !   ncorr(:) = x(:)
  !   write(60, '(*(E15.6E3))') ncorr(:)
  !   close(60)
    
  !   deallocate(dt, dv)
  !   deallocate(invNcorr, invM)
  !   deallocate(x, b, r, d, Mr, Ad)
  ! contains

  !   subroutine apply_fourier_mat(vec, mat, res, dt, dv, nfft, plan_fwd, plan_back)
  !     integer*8,  intent(in)   :: plan_fwd, plan_back
  !     !real(dp),   intent(in)   :: alpha, fknee, wn, samprate
  !     real(dp),  dimension(:), intent(in)    :: mat, vec
  !     real(dp),  dimension(:), intent(out)   :: res
  !     !real(sp),  dimension(:)  :: apply_fourier_mat
  !     integer(i4b), intent(in) :: nfft 
  !     ! real(dp),     allocatable, dimension(:), intent(inout) :: dt
  !     ! complex(dpc), allocatable, dimension(:), intent(inout) :: dv
  !     real(sp),     allocatable, dimension(:), intent(inout) :: dt
  !     complex(spc), allocatable, dimension(:), intent(inout) :: dv
  !     integer(i4b) :: i, j, k, l, ntod, n
  !     n = nfft / 2 + 1
  !     ntod = size(vec, 1)
  !     dt(1:ntod)           = vec(:)
  !     dt(2*ntod:ntod+1:-1) = dt(1:ntod)

  !     call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
  !     do l = 0, n-1
  !        dv(l) = dv(l) * mat(l)
  !     end do
  !     call sfftw_execute_dft_c2r(plan_back, dv, dt)
  !     dt = dt / nfft
  !     res(:) = dt(1:ntod)
      
  !   end subroutine apply_fourier_mat
  ! end subroutine get_ncorr_pcg




  subroutine find_d_prime_spikes(self, scan, det, d_prime, pix)
    implicit none
    class(comm_tod),                 intent(in)    :: self
    integer(i4b),                    intent(in)    :: scan, det
    real(sp),          dimension(:), intent(in)    :: d_prime
    integer(i4b),  dimension(1:,1:), intent(in)    :: pix
    integer(i4b) :: i, j, l, k, n, m, nomp, ntod, ndet, err, n_downsamp, n_short
    integer(i4b) :: nfft, nbuff, j_end, j_start, sampnum
    real(dp)     :: nu, power, n_sigma, rms, avg
    real(dp),     allocatable, dimension(:) :: d_downsamp, backup
    logical(lgt) :: found_spike
    character(len=1024) :: filename
    
    ntod = self%scans(scan)%ntod
    n_downsamp = floor(self%samprate)
    n = ntod - mod(ntod, n_downsamp)
    n_short = n / n_downsamp

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
  subroutine sample_noise_psd(self, handle, scan, mask, s_tot, n_corr)
    implicit none
    class(comm_tod),             intent(inout)  :: self
    type(planck_rng),                intent(inout)  :: handle
    integer(i4b),                    intent(in)     :: scan
    real(sp),        dimension(:,:), intent(in)     :: mask, s_tot, n_corr
    
    integer*8    :: plan_fwd
    integer(i4b) :: i, j, n, n_bins, l, nomp, omp_get_max_threads, err, ntod, n_f 
    integer(i4b) :: ndet
    real(dp)     :: s, res, log_nu, samprate, gain, dlog_nu, nu, f
    real(dp)     :: alpha, sigma0, fknee, x_in(3), prior_fknee(2), prior_alpha(2), alpha_dpc, fknee_dpc
    real(sp),     allocatable, dimension(:) :: dt, ps
    complex(spc), allocatable, dimension(:) :: dv
    real(sp),     allocatable, dimension(:) :: d_prime
    character(len=1024) :: filename

    
    ntod = self%scans(scan)%ntod
    ndet = self%ndet
    nomp = omp_get_max_threads()

    if (trim(self%freq) == '030') then
       prior_fknee = [0.010d0, 0.45d0]
       prior_alpha = [-2.5d0, -0.4d0]
    else if (trim(self%freq) == '044') then
       prior_fknee = [0.002d0, 0.40d0]
       prior_alpha = [-2.5d0, -0.4d0]
    else if (trim(self%freq) == '070') then
       prior_fknee = [0.001d0, 0.25d0]
       prior_alpha = [-3.0d0, -0.4d0]
    else if (trim(self%freq) == '023-WMAP_K') then
       prior_fknee = [0.01d0,1.0d0]
       prior_alpha = [-3d0, -0.4d0]
    else if (trim(self%freq) == '061-WMAP_V2') then
       prior_fknee = [0.0001d0,1.0d0]
       prior_alpha = [-3d0, -0.4d0]
    else 
       prior_fknee = [0.02d0,1.0d0]
       prior_alpha = [-3d0, -0.4d0]
    end if

    ! compute sigma_0 the old way
    do i = 1, ndet
       if (.not. self%scans(scan)%d(i)%accept) cycle
       s = 0.d0
       n = 0

       do j = 1, self%scans(scan)%ntod-1
          if (any(mask(j:j+1,i) < 0.5)) cycle
          res = (self%scans(scan)%d(i)%tod(j) - &
               & (self%scans(scan)%d(i)%gain * s_tot(j,i) + &
               & n_corr(j,i)) - &
               & (self%scans(scan)%d(i)%tod(j+1) - &
               & (self%scans(scan)%d(i)%gain * s_tot(j+1,i) + &
               & n_corr(j+1,i))))/sqrt(2.)
          s = s + res**2
          n = n + 1
       end do
       ! if ((i == 1) .and. (scan == 1)) then
       !    write(*,*) "sigma0: ", sqrt(s/(n-1))
       ! end if
       if (n > 100) self%scans(scan)%d(i)%sigma0 = sqrt(s/(n-1))
    end do

    ! return
    ! if (trim(self%freq) == '044') return
    ! if (trim(self%freq) == '070') return

!    n = ntod + 1
    n = ntod/2 + 1

    call sfftw_init_threads(err)
    call sfftw_plan_with_nthreads(nomp)

!    allocate(dt(2*ntod), dv(0:n-1))
!    call sfftw_plan_dft_r2c_1d(plan_fwd,  2*ntod, dt, dv, fftw_estimate + fftw_unaligned)
    allocate(dt(ntod), dv(0:n-1))
    call sfftw_plan_dft_r2c_1d(plan_fwd,  ntod, dt, dv, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)

    ! Commented out OMP since we had problems with global parameters 
    ! in the likelihood functions
!    !$OMP PARALLEL PRIVATE(i,l,j,dt,dv,f,d_prime,gain,ps,sigma0,alpha,fknee,samprate)
!    allocate(dt(2*ntod), dv(0:n-1))
    allocate(dt(ntod), dv(0:n-1))
    
    allocate(ps(0:n-1))
!    !$OMP DO SCHEDULE(guided)
    do i = 1, ndet
       if (.not. self%scans(scan)%d(i)%accept) cycle
       dt(1:ntod) = n_corr(:,i)
       !dt(2*ntod:ntod+1:-1) = dt(1:ntod)

       ps(:) = 0
       
       samprate = self%samprate
       alpha    = min(max(self%scans(scan)%d(i)%alpha,prior_alpha(1)), prior_alpha(2))
       sigma0   = abs(self%scans(scan)%d(i)%sigma0)
       fknee    = min(max(self%scans(scan)%d(i)%fknee,prior_fknee(1)), prior_fknee(2))
       
       call sfftw_execute_dft_r2c(plan_fwd, dt, dv)

       fknee_dpc = self%scans(scan)%d(i)%fknee_def

       ! n_f should be the index representing fknee
       ! we want to only use smaller frequencies than this in the likelihood
       n_f = ceiling(2 * fknee_dpc * (n-1) / (samprate/2)) !ceiling(0.01d0 * (n-1) / (samprate/2)) ! n-1 
       !n_f = 1000
       do l = 1, n_f !n-1
          ps(l) = abs(dv(l)) ** 2 / ntod          
       end do

       ! Sampling noise parameters given n_corr
       ! TODO: get prior parameters from parameter file
       ! Sampling fknee, in Hz
       ! FOr WMAP, they report the "optimal time-domain filters", i.e., noise
       ! autocorrelation functions, rather than PSDs. But Table 2 of Jarosik et
       ! al. (2003) (On-orbit radiometer characterization) they report in Table
       ! 1 for each of the 20 radiometers fknee in mHz;
       ! K11    K12     Ka11    Ka12
       ! 0.40   0.51    0.71    0.32    
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Q11    Q12     Q21     Q22
       ! 1.09   0.35    5.76    8.62
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! V11    V12     V21     V22
       ! 0.09   1.41    0.88    8.35    
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! W111   W112    W211    W212
       ! 7.88   0.66    9.02    7.47
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! W311   W312    W411    W412
       ! 0.93   0.28   46.5    26.0
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! As far as I can tell, this is the only place where the fknees are
       ! reported...

       x_in(1) = max(fknee - 0.5 * fknee, prior_fknee(1))
       x_in(3) = min(fknee + 0.5 * fknee, prior_fknee(2))
       x_in(2) = 0.5 * (x_in(1) + x_in(3))
       fknee = sample_InvSamp(handle, x_in, lnL_fknee, prior_fknee)
       if ((fknee < prior_fknee(1)) .or. (fknee > prior_fknee(2))) then
          fknee = min(max(self%scans(scan)%d(i)%fknee,prior_fknee(1)), prior_fknee(2))
       end if


       ! Sampling alpha
       alpha_dpc = self%scans(scan)%d(i)%alpha_def
       x_in(1)   = max(alpha - 0.2 * abs(alpha), prior_alpha(1))
       x_in(3)   = min(alpha + 0.2 * abs(alpha), prior_alpha(2))
       x_in(2)   = 0.5 * (x_in(1) + x_in(3))
       alpha = sample_InvSamp(handle, x_in, lnL_alpha, prior_alpha)
       if ((alpha < prior_alpha(1)) .or. (alpha > prior_alpha(2))) then
          alpha = min(max(self%scans(scan)%d(i)%alpha,prior_alpha(1)), prior_alpha(2))
       end if

       self%scans(scan)%d(i)%alpha = alpha
       self%scans(scan)%d(i)%fknee = fknee
    end do
!    !$OMP END DO
    deallocate(dt, dv)
    deallocate(ps)
    
!    !$OMP END PARALLEL
    
    call sfftw_destroy_plan(plan_fwd)
    
  contains

    function lnL_fknee(x) 
      use healpix_types
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: lnL_fknee, sconst, f, s

      if (x <= 1e-6 .or. x > 10.d0) then
         lnL_fknee = -1.d30
         return
      end if
      lnL_fknee = 0.d0
      sconst = sigma0 ** 2 * x ** (-alpha) 

      lnL_fknee = lnL_fknee - 0.5d0 * (log(x) - log(fknee_dpc)) ** 2 / (0.10d0 * log(10.d0)) ** 2 - log(x)

      do l = 1, n_f  ! n-1
         f = l*(samprate/2)/(n-1)
         if (abs(f-1.d0/60.d0)*60.d0 < 0.05d0) then
             continue
         end if

         s = sconst * f ** (alpha)
         lnL_fknee = lnL_fknee - (ps(l) / s + log(s))
      end do
      
    end function lnL_fknee

    function lnL_alpha(x) 
      use healpix_types
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: lnL_alpha, sconst, f, s
      
      if (abs(x) > 5.d0) then
         lnL_alpha = -1.d30
         return
      end if
      lnL_alpha = 0.d0

      lnL_alpha = lnL_alpha - 0.5d0 * (x - alpha_dpc) ** 2 / 0.2d0 ** 2
      
      ! if (trim(self%freq) == '070') then
      !    lnL_alpha = lnL_alpha - 0.5d0 * (x - alpha_dpc) ** 2 / 0.2d0 ** 2
      ! end if
      
      sconst = sigma0 ** 2 * fknee ** (-x) 
      
      do l = 1, n_f  ! n-1
         f = l*(samprate/2)/(n-1)
         if (abs(f-1.d0/60.d0)*60.d0 < 0.05d0) then
             continue
         end if

         s = sconst * f ** (x)
         lnL_alpha = lnL_alpha - (ps(l) / s + log(s))
      end do
      
    end function lnL_alpha

       
  end subroutine sample_noise_psd



  ! Sample noise psd
  ! TODO: Add fluctuation term if operation == sample
  subroutine sample_noise_psd2(tod, handle, scan, mask, s_tot, n_corr)
    implicit none
    class(comm_tod),                 intent(inout)  :: tod
    type(planck_rng),                intent(inout)  :: handle
    integer(i4b),                    intent(in)     :: scan
    real(sp),        dimension(:,:), intent(in)     :: mask, s_tot, n_corr
    
    integer*8    :: plan_fwd
    integer(i4b) :: i, j, n, n_bins, l, nomp, omp_get_max_threads, err, ntod 
    integer(i4b) :: ndet
    real(dp)     :: s, res, log_nu, samprate, gain, dlog_nu, nu
    
    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv
    real(sp),     allocatable, dimension(:) :: d_prime
    real(dp),     allocatable, dimension(:) :: log_nu_bin_edges, psd, nu_sum
    integer(i4b), allocatable, dimension(:) :: n_modes
    
    ntod = tod%scans(scan)%ntod
    ndet = tod%ndet
    nomp = omp_get_max_threads()
    

    ! compute sigma_0 the old way
    do i = 1, ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle
    
       s = 0.d0
       n = 0

       do j = 1, tod%scans(scan)%ntod-1
          if (any(mask(j:j+1,i) < 0.5)) cycle
          res = (tod%scans(scan)%d(i)%tod(j) - &
               & (tod%scans(scan)%d(i)%gain * s_tot(j,i) + &
               & n_corr(j,i)) - &
               & (tod%scans(scan)%d(i)%tod(j+1) - &
               & (tod%scans(scan)%d(i)%gain * s_tot(j+1,i) + &
               & n_corr(j+1,i))))/sqrt(2.)
          s = s + res**2
          n = n + 1
       end do
       ! if ((i == 1) .and. (scan == 1)) then
       !    write(*,*) "sigma0: ", sqrt(s/(n-1))
       ! end if
       if (n > 100) tod%scans(scan)%d(i)%sigma0 = sqrt(s/(n-1))
    end do

    return
    
    n = ntod + 1
    n_bins = 20

    call sfftw_init_threads(err)
    call sfftw_plan_with_nthreads(nomp)

    allocate(dt(2*ntod), dv(0:n-1))
    call sfftw_plan_dft_r2c_1d(plan_fwd,  2*ntod, dt, dv, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)
    
    do i = 1, ndet
       if (.not. allocated(tod%scans(scan)%d(i)%log_n_psd)) allocate(tod%scans(scan)%d(i)%log_n_psd(n_bins))
       if (.not. allocated(tod%scans(scan)%d(i)%log_nu)) allocate(tod%scans(scan)%d(i)%log_nu(n_bins))
       if (.not. allocated(tod%scans(scan)%d(i)%log_n_psd2)) allocate(tod%scans(scan)%d(i)%log_n_psd2(n_bins))
    end do
    
    !$OMP PARALLEL PRIVATE(i,l,j,dt,dv,nu,log_nu,d_prime,log_nu_bin_edges,n_modes,psd,nu_sum,gain,dlog_nu)
    allocate(dt(2*ntod), dv(0:n-1))
    allocate(d_prime(ntod))
    
    allocate(log_nu_bin_edges(n_bins + 1))
    allocate(n_modes(n_bins))
    allocate(psd(n_bins))
    allocate(nu_sum(n_bins))
    !$OMP DO SCHEDULE(guided)
    do i = 1, ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle
    
       !gain = tod%scans(scan)%d(i)%gain
       !d_prime(:) = tod%scans(scan)%d(i)%tod(:) - s_tot(:, i) * gain
       dt(1:ntod) = n_corr(:,i)
       ! do j = 1, ntod
       !    if (mask(j,i) == 0.) then
       !       dt(j) = n_corr(j,i) + rand_gauss(handle) * tod%scans(scan)%d(i)%sigma0

       !    else
       !       dt(j) = d_prime(j)
       !    end if
       ! end do
       !dt(1:ntod)           = d_prime(:)
       !dt(1:ntod)           = d_prime(:)
       dt(2*ntod:ntod+1:-1) = dt(1:ntod)

       ! if ((i == 1) .and. (scan == 1)) then
       !    open(65,file='ps_dt.dat')
       !    do j = 1, ntod
       !       write(65, *) dt(j)
       !    end do
       !    close(65)
       ! end if
              
       n_modes(:) = 0
       psd(:) = 0.d0
       nu_sum(:) = 0.d0
       samprate = tod%samprate
       call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
       call linspace(log(1*(samprate/2)/(n-1)),log(samprate/2), log_nu_bin_edges)

       dlog_nu = log_nu_bin_edges(2) - log_nu_bin_edges(1)
       
       !if ((i == 1) .and. (scan == 1)) open(65,file='ps.dat') 

       do l = 1, n-1
          nu = l*(samprate/2)/(n-1)
          log_nu = log(nu)
          j = min(ceiling((log_nu - log_nu_bin_edges(1) + 1d-10) / dlog_nu), n_bins)
          n_modes(j) = n_modes(j) + 1
          psd(j) = psd(j) + abs(dv(l)) ** 2
          nu_sum(j) = nu_sum(j) + nu
          !if ((i == 1) .and. (scan == 1)) write(65, *) abs(dv(l)) ** 2, log_nu
       end do
       !if ((i == 1) .and. (scan == 1)) close(65)
       
       if (trim(tod%operation) == "sample") then
          ! use abs to prevent rare cases of negative power spectra
          ! (should have been samples from inverse gamma distribution)
!          write(*,*) "sampling!!!!!!!"
          tod%scans(scan)%d(i)%log_n_psd(:) =  log(psd(:) / n_modes(:)) ! &
              ! & * abs(1.d0 + sqrt(2.d0 / n_modes(:)) * rand_gauss(handle))) 
       else
          tod%scans(scan)%d(i)%log_n_psd(:) = log(psd(:) / n_modes(:))
       end if
       tod%scans(scan)%d(i)%log_nu(:) = log(nu_sum(:) / n_modes(:)) !log_nu_bin_edges(1:n_bins) + 0.5d0 * dlog_nu
       !tod%scans(scan)%d(i)%log_nu(:) = log_nu_bin_edges(1:n_bins) + 0.5d0 * dlog_nu
       ! write(*,*) tod%scans(scan)%d(i)%log_n_psd
       ! write(*,*) tod%scans(scan)%d(i)%log_nu
       ! write(*,*) "After"
       ! if ((i == 1) .and. (scan == 1)) then
       !    open(65,file='ps_binned.dat')
       !    do j = 1, n_bins
       !       write(65, *) tod%scans(scan)%d(i)%log_n_psd(j), tod%scans(scan)%d(i)%log_nu(j)
       !    end do
       !    close(65)
       ! end if
       ! tod%scans(scan)%d(i)%sigma0 = sqrt(exp(tod%scans(scan)%d(i)%log_n_psd(n_bins)))
       ! if ((i == 1) .and. (scan == 1)) then
       !    write(*,*) "sigma0: ", tod%scans(scan)%d(i)%sigma0
       
       ! end if
       
       ! if ((i == 1) .and. (scan == 1)) then
       !    open(65,file='ps_tod.dat')
       !    do j = 1, ntod
       !       write(65, *) n_corr(j,i), tod%scans(scan)%d(i)%tod(j), &
       !            & s_tot(j,i), mask(j,i)
       !    end do
       !    close(65)
       ! end if
       ! if ((i == 1) .and. (scan == 1)) then
       !    open(65,file='ps_params.dat')
       !    write(65, '(4(E15.6E3))') tod%scans(scan)%d(i)%gain, tod%scans(scan)%d(i)%alpha, &
       !            & tod%scans(scan)%d(i)%sigma0, tod%scans(scan)%d(i)%fknee
       !    close(65)
       ! end if
    end do
    !$OMP END DO
    deallocate(dt, dv)
    deallocate(d_prime, psd, n_modes, nu_sum)
    deallocate(log_nu_bin_edges)
    
    !$OMP END PARALLEL
    
    call sfftw_destroy_plan(plan_fwd)
    
    do i = 1, ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle
       ! call spline(tod%scans(scan)%d(i)%log_nu,&
       !      tod%scans(scan)%d(i)%log_n_psd,&
       !      0.d0,0.d0,tod%scans(scan)%d(i)%log_n_psd2)
       call spline(tod%scans(scan)%d(i)%log_nu,&
            tod%scans(scan)%d(i)%log_n_psd,&
            1.d30,1.d30,tod%scans(scan)%d(i)%log_n_psd2)
    end do
       
  end subroutine sample_noise_psd2

  ! Compute correlated noise term, n_corr from eq:
  ! ((N_c^-1 + N_wn^-1) n_corr = d_prime + w1 * sqrt(N_wn) + w2 * sqrt(N_c) 
  subroutine sample_n_corr2(tod, handle, scan, mask, s_sub, n_corr)
    implicit none
    class(comm_tod),               intent(in)     :: tod
    type(planck_rng),                  intent(inout)  :: handle
    integer(i4b),                      intent(in)     :: scan
    real(sp),          dimension(:,:), intent(in)     :: mask, s_sub
    real(sp),          dimension(:,:), intent(out)    :: n_corr
    integer(i4b) :: i, j, l, k, n, m, nlive, nomp, ntod, ndet, err, omp_get_max_threads
    integer(i4b) :: nfft, nbuff, j_end, j_start
    integer*8    :: plan_fwd, plan_back
    logical(lgt) :: init_masked_region, end_masked_region
    real(sp)     :: sigma_0, alpha, nu_knee,  samprate, gain, mean, N_wn, N_c
    real(dp)     :: nu, power, fft_norm, t1, t2
    real(sp),     allocatable, dimension(:,:) :: dt
    complex(spc), allocatable, dimension(:,:) :: dv
    real(sp),     allocatable, dimension(:) :: d_prime
    
    ntod = tod%scans(scan)%ntod
    ndet = tod%ndet
!    nomp = omp_get_max_threads()
    m    = count(tod%scans(scan)%d%accept)
 !   nlive = count(tod%scans(scan)%d%accept)
    
    nfft = 2 * ntod
    !nfft = get_closest_fft_magic_number(ceiling(ntod * 1.05d0))
    
    n = nfft / 2 + 1

    call wall_time(t1)
!    call sfftw_init_threads(err)
!    call sfftw_plan_with_nthreads(nomp)

    allocate(dt(nfft,m), dv(0:n-1,m))
!!$    call sfftw_plan_dft_r2c_1d(plan_fwd,  nfft, dt, dv, fftw_patient) 
!!$    call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_patient)
    call sfftw_plan_many_dft_r2c(plan_fwd, 1, nfft, m, dt, &
         & nfft, 1, nfft, dv, n, 1, n, fftw_patient)
    call sfftw_plan_many_dft_c2r(plan_back, 1, nfft, m, dv, &
         & n, 1, n, dt, nfft, 1, nfft, fftw_patient)


!    deallocate(dt, dv)
    call wall_time(t2)
    !if (tod%myid == 0) write(*,*) ' fft1 =', t2-t1 

!    !$OMP PARALLEL PRIVATE(i,j,l,k,dt,dv,nu,sigma_0,alpha,nu_knee,d_prime,init_masked_region,end_masked_region)
!    !$OMP PARALLEL PRIVATE(i,j,l,k,dt,dv,nu,sigma_0,alpha,nu_knee,d_prime)
!    allocate(dt(nfft), dv(0:n-1))
    allocate(d_prime(ntod))

    !!$OMP DO SCHEDULE(guided)
    j = 0
    do i = 1, ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle
       j       = j+1
       gain    = tod%scans(scan)%d(i)%gain  ! Gain in V / K
       d_prime = tod%scans(scan)%d(i)%tod - tod%scans(scan)%d(i)%baseline - S_sub(:,i) * gain
       sigma_0 = tod%scans(scan)%d(i)%sigma0

       call wall_time(t1)       
       call fill_all_masked(d_prime, mask(:,i), ntod, (trim(tod%operation) == "sample"), sigma_0, handle, tod%scans(scan)%chunk_num)
       call wall_time(t2)
!    if (tod%myid == 0) write(*,*) ' fft2 =', t2-t1 

       ! Preparing for fft
       dt(1:ntod,j)           = d_prime
       dt(2*ntod:ntod+1:-1,j) = dt(1:ntod,j)
    end do
  
    call wall_time(t1)
    call sfftw_execute_dft_r2c(plan_fwd, dt, dv)

    j = 0
    do i = 1, ndet
       if (.not. tod%scans(scan)%d(i)%accept) cycle
       j       = j+1
    
       samprate = tod%samprate
       alpha    = tod%scans(scan)%d(i)%alpha
       nu_knee  = tod%scans(scan)%d(i)%fknee
       N_wn     = sigma_0 ** 2           ! white noise power spectrum
       fft_norm = sqrt(1.d0 * nfft)  ! used when adding fluctuation terms to Fourier coeffs (depends on Fourier convention)
       call wall_time(t2)
!    if (tod%myid == 0) write(*,*) ' fft3 =', t2-t1 

       if (trim(tod%operation) == "sample") then
          dv(0,j)    = dv(0,j) + fft_norm * sqrt(N_wn) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0)
       end if
       
       do l = 1, n-1                                                      
          nu = l*(samprate/2)/(n-1)
!!$          if (abs(nu-1.d0/60.d0)*60.d0 < 0.001d0) then
!!$             dv(l) = 0.d0 ! Dont include scan frequency; replace with better solution
!!$          end if
          
          N_c = N_wn * (nu/(nu_knee))**(alpha)  ! correlated noise power spectrum
          if (trim(tod%operation) == "sample") then
             dv(l,j) = (dv(l,j) + fft_norm * ( &
                  sqrt(N_wn) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0) &
                  + N_wn * sqrt(1.0 / N_c) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0) &
                  )) * 1.d0/(1.d0 + N_wn / N_c)
          else
             dv(l,j) = dv(l,j) * 1.0/(1.0 + N_wn/N_c)
          end if
          !if (abs(nu-1.d0/60.d0) < 0.01d0 .and. trim(tod%freq)== '070') dv(l) = 0.d0
       end do
    end do
    call wall_time(t1)
    call sfftw_execute_dft_c2r(plan_back, dv, dt)
    dt = dt / nfft
    call wall_time(t2)
!    if (tod%myid == 0) write(*,*) ' fft4 =', t2-t1 

    j = 0
    do i = 1, ndet
       if (.not. tod%scans(scan)%d(i)%accept) then
          n_corr(:,i) = 0.
          cycle
       end if
       j       = j+1
       n_corr(:,i) = dt(1:ntod,j)
       if (sum(mask(:,i)) > 0 .and. sum(1-mask(:,i)) > 0 .and. i == 1) then
!          write(*,*) tod%scanid(scan), sum(n_corr(:,i)*(1-mask(:,i)))/sum((1-mask(:,i)))-sum(n_corr(:,i)*mask(:,i))/sum(mask(:,i)), 'q4'
       end if
!!$        if (i == 1) then
!!$           open(65,file='ncorr_times.dat', recl=1024)
!!$           do j = 1, ntod
!!$              write(65,*) j, n_corr(j,i), s_sub(j,i), mask(j,i), d_prime(j), tod%scans(scan)%d(i)%tod(j), tod%scans(scan)%d(i)%gain
!!$           end do
!!$           close(65)
!!$           stop
!!$        end if

    end do
!    !$OMP END DO                                                          
    call sfftw_destroy_plan(plan_fwd)                                           
    call sfftw_destroy_plan(plan_back)                                          

    deallocate(dt, dv)
    deallocate(d_prime)
    !deallocate(diff)
 !   !$OMP END PARALLEL
    


  
  end subroutine sample_n_corr2


  ! Routine for multiplying a set of timestreams with inverse noise covariance 
  ! matrix. inp and res have dimensions (ntime,ndet,ninput), where ninput is 
  ! the number of input timestreams (e.g. 2 if you input s_tot and s_orb). 
  ! Here inp and res are assumed to be already allocated. 

  subroutine multiply_inv_N(tod, scan, buffer, sampfreq, pow)
    implicit none
    class(comm_tod),                     intent(in)     :: tod
    integer(i4b),                        intent(in)     :: scan
    real(sp),          dimension(:,:),   intent(inout)     :: buffer !input/output
    real(dp),                            intent(in), optional :: sampfreq, pow
    integer(i4b) :: i, j, l, n, m, nomp, ntod, ndet, err, omp_get_max_threads
    integer*8    :: plan_fwd, plan_back
    real(sp)     :: sigma_0, alpha, nu_knee,  samprate, noise, signal
    real(sp)     :: nu, pow_
    real(sp),     allocatable, dimension(:,:) :: dt
    complex(spc), allocatable, dimension(:,:) :: dv
    
    ntod = size(buffer, 1)
    ndet = size(buffer, 2)
    nomp = omp_get_max_threads()
    m    = count(tod%scans(scan)%d%accept)
    pow_ = 1.d0; if (present(pow)) pow_ = pow

!!$ !   if (tod%myid == 0) open(58,file='invN1.dat')
!!$    allocate(dt(ntod))
!!$    do i = 1, ndet
!!$       if (.not. tod%scans(scan)%d(i)%accept) cycle
!!$       sigma_0  = tod%scans(scan)%d(i)%sigma0
!!$       do j = 1, ninput
!!$          dt = buffer(:,i,j) / sigma_0**2
!!$          dt = dt - mean(1.d0*dt)
!!$          buffer(:,i,j) = dt
!!$       end do
!!$    end do
!!$    deallocate(dt)
!!$!    if (tod%myid == 0) close(58)
!!$    return

    
    n = ntod + 1
    
!    call sfftw_init_threads(err)
!    call sfftw_plan_with_nthreads(nomp)

!!$    allocate(dt(2*ntod), dv(0:n-1))
!!$    call sfftw_plan_dft_r2c_1d(plan_fwd,  2*ntod, dt, dv, fftw_patient)
!!$    call sfftw_plan_dft_c2r_1d(plan_back, 2*ntod, dv, dt, fftw_patient)
!    deallocate(dt, dv)

    allocate(dt(2*ntod,m), dv(0:n-1,m))
    call sfftw_plan_many_dft_r2c(plan_fwd, 1, 2*ntod, m, dt, &
         & 2*ntod, 1, 2*ntod, dv, n, 1, n, fftw_patient)
    call sfftw_plan_many_dft_c2r(plan_back, 1, 2*ntod, m, dv, &
         & n, 1, n, dt, 2*ntod, 1, 2*ntod, fftw_patient)
    
!    if (tod%myid == 0) open(58,file='invN2.dat')
!    !$OMP PARALLEL PRIVATE(i,j,l,dt,dv,nu,sigma_0,alpha,nu_knee,noise,signal)
!    allocate(dt(2*ntod), dv(0:n-1))
    
    !!$OMP DO SCHEDULE(guided)
    j = 0
    do i = 1, ndet
       sigma_0  = abs(real(tod%scans(scan)%d(i)%sigma0,sp))
       if (.not. tod%scans(scan)%d(i)%accept .or. sigma_0 <= 0.d0) cycle
       j = j+1
       dt(1:ntod,j)           = buffer(:,i)
       dt(2*ntod:ntod+1:-1,j) = dt(1:ntod,j)
    end do

    call sfftw_execute_dft_r2c(plan_fwd, dt, dv)

    j = 0
    do i = 1, ndet
       sigma_0  = abs(real(tod%scans(scan)%d(i)%sigma0,sp))
       if (.not. tod%scans(scan)%d(i)%accept .or. sigma_0 <= 0.d0) cycle
       j = j+1
       samprate = real(tod%samprate,sp); if (present(sampfreq)) samprate = real(sampfreq,sp)
       alpha    = real(tod%scans(scan)%d(i)%alpha,sp)
       nu_knee  = real(tod%scans(scan)%d(i)%fknee,sp)
       noise = sigma_0 ** 2 * samprate / tod%samprate
       
       dv(0,j) = 0.d0
       do l = 1, n-1                                                      
          nu      = l*(samprate/2)/(n-1)
          signal  = noise * (nu/(nu_knee))**(alpha)
          dv(l,j) = dv(l,j) * 1.0/(noise + signal)**pow_
       end do
    end do
       
    call sfftw_execute_dft_c2r(plan_back, dv, dt)
    dt = dt / (2*ntod)

!!$          if (tod%myid == 0 .and. i==1 .and. j==1) then
!!$             do k = 1, ntod
!!$                write(58,*) k, dt(k)/dt(30000), buffer(k,i,j)/buffer(30000,i,j)
!!$             end do
!!$          end if

    j = 0
    do i = 1, ndet
       sigma_0  = real(tod%scans(scan)%d(i)%sigma0,sp)
       if (.not. tod%scans(scan)%d(i)%accept .or. sigma_0 <= 0.d0) then
          buffer(:,i)  = 0.d0
       else
          j           = j+1
          buffer(:,i) = dt(1:ntod,j) 
       end if
    end do
    !!$OMP END DO                                                          
    call sfftw_destroy_plan(plan_fwd)                                           
    call sfftw_destroy_plan(plan_back)                                          
    deallocate(dt, dv)
!    !$OMP END PARALLEL
!    if (tod%myid == 0) close(58)


!!$    call mpi_finalize(i)
!!$    stop
  end subroutine multiply_inv_N



end module comm_tod_noise_mod
