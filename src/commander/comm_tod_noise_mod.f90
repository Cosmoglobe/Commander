module comm_tod_noise_mod
  use comm_tod_mod
  use comm_utils
  use InvSamp_mod
  implicit none



contains

 ! Compute correlated noise term, n_corr from eq:
  ! ((N_c^-1 + N_wn^-1) n_corr = d_prime + w1 * sqrt(N_wn) + w2 * sqrt(N_c) 
  subroutine sample_n_corr(self, handle, scan, mask, s_sub, n_corr)
    implicit none
    class(comm_tod),               intent(in)     :: self
    type(planck_rng),                  intent(inout)  :: handle
    integer(i4b),                      intent(in)     :: scan
    real(sp),          dimension(:,:), intent(in)     :: mask, s_sub
    real(sp),          dimension(:,:), intent(out)    :: n_corr
    integer(i4b) :: i, j, l, k, n, m, nomp, ntod, ndet, err, omp_get_max_threads
    integer(i4b) :: nfft, nbuff, j_end, j_start
    integer*8    :: plan_fwd, plan_back
    logical(lgt) :: init_masked_region, end_masked_region
    real(sp)     :: sigma_0, alpha, nu_knee,  samprate, gain, mean, N_wn, N_c
    real(dp)     :: nu, power, fft_norm
    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv
    real(sp),     allocatable, dimension(:) :: d_prime
    
    ntod = self%scans(scan)%ntod
    ndet = self%ndet
    nomp = omp_get_max_threads()
    
    nfft = 2 * ntod
    !nfft = get_closest_fft_magic_number(ceiling(ntod * 1.05d0))
    
    n = nfft / 2 + 1

    call sfftw_init_threads(err)
    call sfftw_plan_with_nthreads(nomp)

    allocate(dt(nfft), dv(0:n-1))
    call sfftw_plan_dft_r2c_1d(plan_fwd,  nfft, dt, dv, fftw_estimate + fftw_unaligned)
    call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)

    !$OMP PARALLEL PRIVATE(i,j,l,k,dt,dv,nu,sigma_0,alpha,nu_knee,d_prime,init_masked_region,end_masked_region)
    allocate(dt(nfft), dv(0:n-1))
    allocate(d_prime(ntod))

    !$OMP DO SCHEDULE(guided)
    do i = 1, ndet
       if (.not. self%scans(scan)%d(i)%accept) cycle
       gain = self%scans(scan)%d(i)%gain  ! Gain in V / K
       d_prime(:) = self%scans(scan)%d(i)%tod(:) - S_sub(:,i) * gain

       sigma_0 = self%scans(scan)%d(i)%sigma0
       
       ! Fill gaps in data 
       init_masked_region = .true.
       end_masked_region = .false.
       do j = 1,ntod
          if (mask(j,i) == 1.) then
             if (end_masked_region) then
                j_end = j - 1
                call fill_masked_region(d_prime, mask(:,i), j_start, j_end, ntod)
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
          call fill_masked_region(d_prime, mask(:,i), j_start, j_end, ntod)
          if (trim(self%operation) == "sample") then
             do k = j_start, j_end
                d_prime(k) = d_prime(k) + sigma_0 * rand_gauss(handle)
             end do
          end if      
       end if
      
       ! Preparing for fft
       dt(1:ntod)           = d_prime(:)
       dt(2*ntod:ntod+1:-1) = dt(1:ntod)
       ! nbuff = nfft - ntod
       ! do j=1, nbuff
       !    dt(ntod+j) = d_prime(ntod) + (d_prime(1) - d_prime(ntod)) * (j-1) / (nbuff - 1)
       ! end do
  
       call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
       samprate = self%samprate
       alpha    = self%scans(scan)%d(i)%alpha
       nu_knee  = self%scans(scan)%d(i)%fknee
       N_wn = sigma_0 ** 2  ! white noise power spectrum
       fft_norm = sqrt(1.d0 * nfft)  ! used when adding fluctuation terms to Fourier coeffs (depends on Fourier convention)
       
       if (trim(self%operation) == "sample") then
          dv(0)    = dv(0) + fft_norm * sqrt(N_wn) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0)
       end if
       
       
       do l = 1, n-1                                                      
          nu = l*(samprate/2)/(n-1)
!!$          if (abs(nu-1.d0/60.d0)*60.d0 < 0.001d0) then
!!$             dv(l) = 0.d0 ! Dont include scan frequency; replace with better solution
!!$          end if
          
          N_c = N_wn * (nu/(nu_knee))**(alpha)  ! correlated noise power spectrum

          if (trim(self%operation) == "sample") then
             dv(l) = (dv(l) + fft_norm * ( &
                  sqrt(N_wn) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0) &
                  + N_wn * sqrt(1.0 / N_c) * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0) &
                  )) * 1.d0/(1.d0 + N_wn / N_c)
          else
             dv(l) = dv(l) * 1.0/(1.0 + N_wn/N_c)
          end if
       end do
       call sfftw_execute_dft_c2r(plan_back, dv, dt)
       dt          = dt / nfft
       n_corr(:,i) = dt(1:ntod) 
       ! if (i == 1) then
       !    open(65,file='ncorr_times.dat')
       !    do j = i, ntod
       !       write(65, '(6(E15.6E3))') n_corr(j,i), s_sub(j,i), mask(j,i), d_prime(j), self%scans(scan)%d(i)%tod(j), self%scans(scan)%d(i)%gain
       !    end do
       !    close(65)
       !    !stop
       ! end if

    end do
    !$OMP END DO                                                          
    deallocate(dt, dv)
    deallocate(d_prime)
    !deallocate(diff)
    !$OMP END PARALLEL
    

    call sfftw_destroy_plan(plan_fwd)                                           
    call sfftw_destroy_plan(plan_back)                                          
  
  end subroutine sample_n_corr


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
    real(dp)     :: alpha, sigma0, fknee, x_in(3), prior(2)
    real(sp),     allocatable, dimension(:) :: dt, ps
    complex(spc), allocatable, dimension(:) :: dv
    real(sp),     allocatable, dimension(:) :: d_prime
    
    ntod = self%scans(scan)%ntod
    ndet = self%ndet
    nomp = omp_get_max_threads()

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

    return
    
    n = ntod + 1

    call sfftw_init_threads(err)
    call sfftw_plan_with_nthreads(nomp)

    allocate(dt(2*ntod), dv(0:n-1))
    call sfftw_plan_dft_r2c_1d(plan_fwd,  2*ntod, dt, dv, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)

    ! Commented out OMP since we had problems with global parameters 
    ! in the likelihood functions
!    !$OMP PARALLEL PRIVATE(i,l,j,dt,dv,f,d_prime,gain,ps,sigma0,alpha,fknee,samprate)
    allocate(dt(2*ntod), dv(0:n-1))
    allocate(d_prime(ntod))
    
    allocate(ps(0:n-1))
!    !$OMP DO SCHEDULE(guided)
    do i = 1, ndet
       if (.not. self%scans(scan)%d(i)%accept) cycle
       dt(1:ntod) = n_corr(:,i)
       dt(2*ntod:ntod+1:-1) = dt(1:ntod)

       ps(:) = 0
       
       samprate = self%samprate
       alpha = self%scans(scan)%d(i)%alpha
       sigma0 = self%scans(scan)%d(i)%sigma0
       fknee = self%scans(scan)%d(i)%fknee
       
       call sfftw_execute_dft_r2c(plan_fwd, dt, dv)

       ! n_f should be the index representing fknee
       ! we want to only use smaller frequencies than this in the likelihood
       n_f = ceiling(fknee * (n-1) / (samprate/2))  

       do l = 1, n_f !n-1
          ps(l) = abs(dv(l)) ** 2 / (2 * ntod)
       end do
 
       ! Sampling noise parameters given n_corr
       
       ! TODO: get prior parameters from parameter file
       ! Sampling fknee
       x_in(1) = fknee - 0.02
       x_in(2) = fknee
       x_in(3) = fknee + 0.02
       prior(1) = 0.05
       prior(2) = 0.3
       fknee = sample_InvSamp(handle, x_in, lnL_fknee, prior)
       
       ! Sampling alpha
       x_in(1) = alpha - 0.1
       x_in(2) = alpha
       x_in(3) = alpha + 0.1
       prior(1) = -0.5
       prior(2) = -2.0
       alpha = sample_InvSamp(handle, x_in, lnL_alpha, prior)
       
       self%scans(scan)%d(i)%alpha = alpha
       self%scans(scan)%d(i)%fknee = fknee
    end do
!    !$OMP END DO
    deallocate(dt, dv)
    deallocate(d_prime)
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
         lnL_fknee = 1.d30
         return
      end if
      lnL_fknee = 0.d0
      sconst = sigma0 ** 2 * x ** alpha 
      do l = 1, n_f  ! n-1
         f = l*(samprate/2)/(n-1)
         s = sconst * f ** (-alpha)
         lnL_fknee = lnL_fknee - 0.5 * (ps(l) / s + log(s))
      end do
    end function lnL_fknee

    function lnL_alpha(x) 
      use healpix_types
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: lnL_alpha, sconst, f, s
      if (abs(x) > 5.d0) then
         lnL_alpha = 1.d30
         return
      end if
      lnL_alpha = 0.d0
      sconst = sigma0 ** 2 * fknee ** x 
      do l = 1, n_f  ! n-1
         f = l*(samprate/2)/(n-1)
         s = sconst * f ** (-x)
         lnL_alpha = lnL_alpha - 0.5 * (ps(l) / s + log(s))
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
       d_prime = tod%scans(scan)%d(i)%tod - S_sub(:,i) * gain
       sigma_0 = tod%scans(scan)%d(i)%sigma0

       call wall_time(t1)       
       call fill_all_masked(d_prime, mask(:,i), ntod, (trim(tod%operation) == "sample"), sigma_0, handle)
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
       sigma_0  = real(tod%scans(scan)%d(i)%sigma0,sp)
       if (.not. tod%scans(scan)%d(i)%accept .or. sigma_0 <= 0.d0) cycle
       j = j+1
       dt(1:ntod,j)           = buffer(:,i)
       dt(2*ntod:ntod+1:-1,j) = dt(1:ntod,j)
    end do

    call sfftw_execute_dft_r2c(plan_fwd, dt, dv)

    j = 0
    do i = 1, ndet
       sigma_0  = real(tod%scans(scan)%d(i)%sigma0,sp)
       if (.not. tod%scans(scan)%d(i)%accept .or. sigma_0 <= 0.d0) cycle
       j = j+1
       samprate = real(tod%samprate,sp); if (present(sampfreq)) samprate = real(sampfreq,sp)
       alpha    = real(tod%scans(scan)%d(i)%alpha,sp)
       nu_knee  = real(tod%scans(scan)%d(i)%fknee,sp)
       noise    = sigma_0 ** 2
       
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
