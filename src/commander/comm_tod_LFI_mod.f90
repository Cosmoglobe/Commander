module comm_tod_LFI_mod
  use comm_tod_mod
  use comm_param_mod
  use comm_map_mod
  implicit none
  include 'fftw3.f'

  private
  public comm_LFI_tod

  type, extends(comm_tod) :: comm_LFI_tod 
   contains
     procedure     :: process_tod        => process_LFI_tod
     procedure     :: compute_binned_map
     procedure     :: project_sky
     procedure     :: compute_orbital_dipole
     procedure     :: sample_gain
     procedure     :: sample_n_corr
     procedure     :: compute_cleaned_tod
     procedure     :: construct_sl_template
     procedure     :: compute_chisq
     procedure     :: sample_noise_psd
  end type comm_LFI_tod

  interface comm_LFI_tod
     procedure constructor
  end interface comm_LFI_tod

contains

  !**************************************************
  !             Constructor
  !**************************************************
  function constructor(cpar, info)
    implicit none
    type(comm_params),       intent(in) :: cpar
    class(comm_mapinfo),     target     :: info
    class(comm_LFI_tod),     pointer    :: constructor

    ! Set up LFI specific parameters
    allocate(constructor)
    constructor%myid     = cpar%myid
    constructor%numprocs = cpar%numprocs
    constructor%info     => info

    ! Test code, just to be able to read a single file; need to settle on parameter structure
!    call constructor%get_scan_ids("data/filelist_5files.txt")
     call constructor%get_scan_ids("data/filelist_1year.txt")
!    call constructor%get_scan_ids("data/filelist.txt")
!    call constructor%get_scan_ids("data/filelist_2half.txt")

    constructor%procmask => comm_map(constructor%info, "data/test_procmask_expand20arcmin_n512.fits")

    constructor%nmaps    = info%nmaps
    constructor%ndet     = 4
    constructor%nhorn    = 1
    constructor%samprate = 50.d0
    allocate(constructor%stokes(constructor%nmaps))
    allocate(constructor%w(constructor%nmaps,constructor%nhorn,constructor%ndet))
    constructor%stokes = [1,2,3]
    constructor%w      = 1.d0

    ! Read the actual TOD
    call constructor%read_tod

  end function constructor

  !**************************************************
  !             Driver routine
  !**************************************************
  subroutine process_LFI_tod(self, map_in, map_out, rms_out)
    implicit none
    class(comm_LFI_tod),               intent(inout)    :: self
    type(map_ptr),       dimension(:), intent(in)    :: map_in            ! One map per detector
    class(comm_map),                   intent(inout) :: map_out, rms_out ! Combined output map and rms

    integer(i4b) :: i, j, k, ntod, ndet, nside, npix, nmaps, naccept, ntot, ierr
    real(dp)     :: t1, t2, t3, t4, t5, t6, chisq_threshold
    real(dp)     :: t_tot(10) 
    real(sp), allocatable, dimension(:,:)   :: n_corr, s_sl, d_calib, s_sky, s_orb, mask
    real(dp), allocatable, dimension(:,:,:) :: map_sky
    real(dp), allocatable, dimension(:)     :: A_abscal, b_abscal
    real(dp), allocatable, dimension(:,:,:) :: A_map
    real(dp), allocatable, dimension(:,:)   :: b_map
    real(dp), allocatable, dimension(:,:)   :: procmask

    t_tot   = 0.d0
    call wall_time(t5)

    ! Set up full-sky map structures
    chisq_threshold = 7.d0
    ndet  = self%ndet
    nside = map_out%info%nside
    nmaps = map_out%info%nmaps
    npix  = 12*nside**2
    allocate(A_map(nmaps,nmaps,0:npix-1), b_map(nmaps,0:npix-1))
    allocate(A_abscal(self%ndet), b_abscal(self%ndet))
    allocate(map_sky(0:npix-1,nmaps,ndet))
    allocate(procmask(0:npix-1,nmaps))
    ! This step should be optimized -- currently takes 6 seconds..
    call wall_time(t1)
    do i = 1, self%ndet
       call map_in(i)%p%bcast_fullsky_map(map_sky(:,:,i))
    end do
    call self%procmask%bcast_fullsky_map(procmask)
    map_sky = map_sky * 1.d-6 ! Gain in V/K
    call wall_time(t2)
    t_tot(9) = t2-t1

    ! Compute output map and rms
    A_map = 0.d0
    b_map = 0.d0
    !A_abscal = 0.d0
    !b_abscal = 0.d0
    naccept = 0
    ntot    = 0
    do i = 1, self%nscan

       call wall_time(t3)

       ! Short-cuts to local variables
       ntod = self%scans(i)%ntod
       ndet = self%ndet

       ! Set up local data structure for current scan
       allocate(d_calib(ntod, ndet))            ! Calibrated and cleaned TOD in uKcmb
       allocate(n_corr(ntod, ndet))             ! Correlated noise in V
       allocate(s_sl(ntod, ndet))               ! Sidelobe in uKcmb
       allocate(s_sky(ntod, ndet))              ! Stationary sky signal in uKcmb
       allocate(s_orb(ntod, ndet))              ! Orbital dipole in uKcmb
       allocate(mask(ntod, ndet))              ! Processing mask in time

       ! Initializing arrays to zero
       n_corr = 0.d0
       s_sl   = 0.d0
       s_sky  = 0.d0
       s_orb  = 0.d0

       ! --------------------
       ! Analyze current scan
       ! --------------------

       ! Construct sky signal template -- Maksym -- this week
       call wall_time(t1)
       do j = 1, ndet
          call self%project_sky(map_sky(:,:,j), procmask, i, j, s_sky(:,j), mask(:,j))  ! scan_id, det,  s_sky(:, j))
       end do
       call wall_time(t2)
       t_tot(1) = t_tot(1) + t2-t1
       !if (self%myid == 0) write(*,*) 'Project    = ', t2-t1

       ! Construct orbital dipole template -- Kristian -- this week-ish
       call wall_time(t1)
       call self%compute_orbital_dipole(i, s_orb)
       call wall_time(t2)
       t_tot(2) = t_tot(2) + t2-t1
       !if (self%myid == 0) write(*,*) 'Orb dipole = ', t2-t1

       ! Construct sidelobe template -- Mathew -- long term
       !call self%construct_sl_template()

       ! Fit correlated noise -- Haavard -- this week-ish
       call wall_time(t1)
       call self%sample_n_corr(i, mask, s_sky, s_sl, s_orb, n_corr)
       call wall_time(t2)
       t_tot(3) = t_tot(3) + t2-t1
       !if (self%myid == 0) write(*,*) 'ncorr      = ', t2-t1

       ! Fit gain for current scan -- Eirik -- this week
       call wall_time(t1)
       do j = 1, ndet
          call sample_gain(self, j, i, n_corr(:, j), mask(:,j), s_sky(:, j), s_sl(:, j), s_orb(:, j))
       end do
       call wall_time(t2)
       t_tot(4) = t_tot(4) + t2-t1
       !if (self%myid == 0) write(*,*) 'gain       = ', t2-t1

       ! .. Compute contribution to absolute calibration from current scan .. -- let's see

       ! Compute bandpass corrections, as in s_sky(i) - <s_sky> -- Trygve, after deadline

       ! Compute clean and calibrated TOD -- Mathew -- this week
       call wall_time(t1)
       call self%compute_cleaned_tod(ntod, i, s_orb, s_sl, n_corr, d_calib)
       call wall_time(t2)
       t_tot(5) = t_tot(5) + t2-t1
       !if (self%myid == 0) write(*,*) 'clean      = ', t2-t1

       ! Compute noise spectrum
       call wall_time(t1)
       do j = 1, ndet
          call self%sample_noise_psd(i, j, mask(:,j), s_sky(:,j), s_sl(:,j), s_orb(:,j), n_corr(:,j))
       end do
       call wall_time(t2)
       t_tot(6) = t_tot(6) + t2-t1
       !if (self%myid == 0) write(*,*) 'noise      = ', t2-t1

       ! Compute chisquare
       call wall_time(t1)
       do j = 1, ndet
          call self%compute_chisq(i, j, mask(:,j), s_sky(:,j), s_sl(:,j), s_orb(:,j), n_corr(:,j))
       end do
       call wall_time(t2)
       t_tot(7) = t_tot(7) + t2-t1
       !if (self%myid == 0) write(*,*) 'chisq      = ', t2-t1

       ! Compute binned map from cleaned TOD -- Marie -- this week
       call wall_time(t1)
       do j = 1, ndet
          ntot= ntot + 1
          if (abs(self%scans(i)%d(j)%chisq) > chisq_threshold .or. &
               & self%scans(i)%d(j)%chisq /= self%scans(i)%d(j)%chisq) &
               & cycle
          naccept = naccept + 1
          call self%compute_binned_map(d_calib, A_map, b_map, i, j)
       end do
       call wall_time(t2)
       t_tot(8) = t_tot(8) + t2-t1
       !if (self%myid == 0) write(*,*) 'bin        = ', t2-t1

       ! Clean up
       deallocate(n_corr, s_sl, s_sky, s_orb, d_calib, mask)

       call wall_time(t4)
       do j = 1, ndet
          if (abs(self%scans(i)%d(j)%chisq) > chisq_threshold .or. &
               & self%scans(i)%d(j)%chisq /= self%scans(i)%d(j)%chisq) &
               write(*,fmt='(a,i8,i5,a,f12.1)') 'Reject scan, det = ', &
               & self%scanid(i), j, ', chisq = ', self%scans(i)%d(j)%chisq
       end do
    end do

    ! Compute absolute calibration, summed over all scans, and rescale output maps

    ! Solve combined map, summed over all pixels -- Marie -- weeks
    !call finalize_binned_map(A_map, b_map, map_tot, rms_tot)
    !map_out%map = map_tot(map_out%info%pix,:)
    !rms_out%map = rms_tot(rms_out%info%pix,:)
    call wall_time(t1)
    call finalize_binned_map(A_map, b_map, map_out, rms_out)
    call wall_time(t2)
    t_tot(10) = t2-t1

    call mpi_reduce(ntot, MPI_IN_PLACE, 1, MPI_INTEGER, MPI_SUM, 0, self%info%comm, ierr)
    call mpi_reduce(naccept, MPI_IN_PLACE, 1, MPI_INTEGER, MPI_SUM, 0, self%info%comm, ierr)

    call wall_time(t6)
    if (self%myid == 0) then
       write(*,*) '  Time dist sky   = ', t_tot(9)
       write(*,*) '  Time project    = ', t_tot(1)
       write(*,*) '  Time orbital    = ', t_tot(2)
       write(*,*) '  Time ncorr      = ', t_tot(3)
       write(*,*) '  Time gain       = ', t_tot(4)
       write(*,*) '  Time clean      = ', t_tot(5)
       write(*,*) '  Time noise      = ', t_tot(6)
       write(*,*) '  Time chisq      = ', t_tot(7)
       write(*,*) '  Time bin        = ', t_tot(8)
       write(*,*) '  Time final      = ', t_tot(10)
       write(*,*) '  Time total      = ', t6-t5, ', accept rate = ', real(naccept,sp) / ntot
    end if

    call map_out%writeFITS('map.fits')
    call rms_out%writeFITS('rms.fits')

    ! Clean up temporary arrays
    deallocate(map_sky, A_map, b_map, procmask)
    deallocate(A_abscal, b_abscal)

  end subroutine process_LFI_tod


  !**************************************************
  !             Sub-process routines
  !**************************************************
  ! Sky signal template
  subroutine project_sky(self, map, pmask, scan_id, det, s_sky, tmask)
    implicit none
    class(comm_LFI_tod),                  intent(in)  :: self
    real(dp),            dimension(:,:),  intent(in)  :: map, pmask
    integer(i4b),                         intent(in)  :: scan_id, det
    real(sp),            dimension(:),    intent(out) :: s_sky, tmask

    integer(i4b)                                      :: i, pix
    real(dp)                                          :: psi

    ! s = T + Q * cos(2 * psi) + U * sin(2 * psi)
    ! T - temperature; Q, U - Stoke's parameters
    do i = 1, self%scans(scan_id)%ntod
       pix = self%scans(scan_id)%d(det)%pix(i)
       psi = self%scans(scan_id)%d(det)%psi(i) + self%scans(scan_id)%d(det)%polang * DEG2RAD
       s_sky(i) = map(pix, 1) + map(pix, 2) * cos(2.d0 * psi) + map(pix, 3) * sin(2.d0 * psi)
       if (any(pmask(pix,:) < 0.5d0)) then
          tmask(i) = 0.
       else
          tmask(i) = 1.
       end if
    end do

  end subroutine project_sky


  ! Compute map with white noise assumption from correlated noise 
  ! corrected and calibrated data, d' = (d-n_corr-n_temp)/gain 
  subroutine compute_orbital_dipole(self, ind, s_orb)
    implicit none
    class(comm_LFI_tod),                 intent(in)  :: self
    integer(i4b),                        intent(in)  :: ind !scan nr/index
    real(sp),            dimension(:,:), intent(out) :: s_orb
    real(dp)             :: x, T_0, q, pix_dir(3), b, b_dot
    real(dp), parameter  :: h = 6.62607015d-34   ! Planck's constant [Js]
    integer(i4b)         :: i,j

    !T_0 = T_CMB*k_b/h !T_0 = T_CMB frequency
    !x = freq * 1.d9 / (2.d0*T_0) !freq is the center bandpass frequancy of the detector
    !q = x * (exp(2.d0*x)+1) / (exp(2.d0*x)-1) !frequency dependency of the quadrupole
    !b = sqrt(sum(self%scans(ind)%v_sun**2))/c !beta for the given scan

    do i = 1,self%ndet
       do j=1,self%scans(ind)%ntod !length of the tod
          call pix2vec_ring(self%scans(ind)%d(i)%nside, self%scans(ind)%d(i)%pix(j), &
               & pix_dir)
          b_dot = dot_product(self%scans(ind)%v_sun, pix_dir)/c
          s_orb(j,i) = T_CMB  * b_dot !only dipole, 1.d6 to make it uK, as [T_CMB] = K
          !s_orb(j,i) = T_CMB  * 1.d6 * b_dot !only dipole, 1.d6 to make it uK, as [T_CMB] = K
          !s_orb(j,i) = T_CMB * 1.d6 * (b_dot + q*b_dot**2) ! with quadrupole
          !s_orb(j,i) = T_CMB * 1.d6 * (b_dot + q*((b_dot**2) - (1.d0/3.d0)*(b**2))) ! net zero monopole
       end do
   end do

  end subroutine compute_orbital_dipole

  ! Compute correlated noise term, n_corr
  subroutine sample_n_corr(self, scan, mask, s_sky, s_sl, s_orb, n_corr)
    implicit none
    class(comm_LFI_tod),               intent(in)     :: self
    integer(i4b),                      intent(in)     :: scan
    real(sp),          dimension(:,:), intent(in)     :: mask, s_sky, s_sl, s_orb
    real(sp),          dimension(:,:), intent(out)    :: n_corr
    integer(i4b) :: i, j, k, l, n, nomp, ntod, ndet, err, omp_get_max_threads
    integer*8    :: plan_fwd, plan_back
    real(sp)     :: sigma_0, alpha, nu_knee, nu, samprate, gain
    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv
    real(sp),     allocatable, dimension(:) :: d_prime
    
    ntod = self%scans(scan)%ntod
    ndet = self%ndet
    nomp = omp_get_max_threads()
    samprate = self%samprate

!    do i = 1, ndet
!       gain = 1.d-6 * self%scans(scan)%d(i)%gain  ! Gain in V / muK
!       n_corr(:,i) = self%scans(scan)%d(i)%tod(:) - S_sky(:,i) * gain 
!    end do
!    return

    
    n = ntod + 1

    call sfftw_init_threads(err)
    call sfftw_plan_with_nthreads(nomp)

    allocate(dt(2*ntod), dv(0:n-1))
    call sfftw_plan_dft_r2c_1d(plan_fwd,  2*ntod, dt, dv, fftw_estimate + fftw_unaligned)
    call sfftw_plan_dft_c2r_1d(plan_back, 2*ntod, dv, dt, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)
    !$OMP PARALLEL PRIVATE(i,l,dt,dv,nu,sigma_0,alpha,nu_knee)
    allocate(dt(2*ntod), dv(0:n-1))
    allocate(d_prime(ntod))
    !$OMP DO SCHEDULE(guided)
    do i = 1, ndet
       gain = self%scans(scan)%d(i)%gain  ! Gain in V / K
!       where (mask(:,i) == 1.)
          d_prime(:) = self%scans(scan)%d(i)%tod(:) - S_sl(:,i) - (S_sky(:,i) + S_orb(:,i)) * gain
!       elsewhere
          ! Do gap filling, placeholder
!          d_prime(:) = self%scans(scan)%d(i)%tod(:) - S_sl(:,i) - (S_sky(:,i) + S_orb(:,i)) * gain 
!       end where
       ! if (i == 1 .and. scan == 1) then
       !    open(22, file="tod.unf", form="unformatted") ! Adjusted open statement
       !    write(22) self%scans(scan)%d(i)%tod(:)
       !    close(22)
       ! end if
       
       ! if (i == 1 .and. scan == 1) then
       !    open(22, file="d_prime.unf", form="unformatted") ! Adjusted open statement
       !    write(22) d_prime
       !    close(22)
       ! end if

       ! if (i == 1 .and. scan == 1) then
       !    open(22, file="sky.unf", form="unformatted") ! Adjusted open statement
       !    write(22) S_sky(:,i) * gain
       !    close(22)
       ! end if

       dt(1:ntod)           = d_prime(:)
       dt(2*ntod:ntod+1:-1) = dt(1:ntod)
       call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
       sigma_0 = self%scans(scan)%d(i)%sigma0
       alpha = self%scans(scan)%d(i)%alpha
       nu_knee = self%scans(scan)%d(i)%fknee
       do l = 0, n-1                                                      
          nu = l*(samprate/2)/(n-1)
          dv(l) = dv(l) * 1.d0/(1.d0 + (nu/(nu_knee))**(-alpha))          
       end do
       call sfftw_execute_dft_c2r(plan_back, dv, dt)
       dt          = dt / (2*ntod)
       n_corr(:,i) = dt(1:ntod)
       ! if (i == 1 .and. scan == 1) then
       !    open(22, file="n_corr.unf", form="unformatted") ! Adjusted open statement
       !    write(22) n_corr(:,i)
       !    close(22)
       ! end if

!!$       if (self%myid == 0 .and. i == 2 .and. scan == 2 .and. .false.) then
!!$          open(58,file='tod.dat')
!!$          do j = 1, ntod
!!$             write(58,*) j, d_prime(j), n_corr(j,i)
!!$          end do
!!$          close(58)
!!$       end if

    end do
    !$OMP END DO                                                          
    deallocate(dt, dv)                                      
    deallocate(d_prime)
    !$OMP END PARALLEL

    call sfftw_destroy_plan(plan_fwd)                                           
    call sfftw_destroy_plan(plan_back)                                          
  
  end subroutine sample_n_corr


  ! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  subroutine sample_gain(self, det, scan_id, n_corr, mask, s_sky, s_sl, s_orb)
    implicit none
    class(comm_LFI_tod),               intent(inout)  :: self
    integer(i4b),                      intent(in)     :: det, scan_id
    real(sp),            dimension(:), intent(in)     :: n_corr, mask, s_sky, s_sl, s_orb
    real(dp),            allocatable,  dimension(:)   :: d_only_wn
    real(dp),            allocatable,  dimension(:)   :: gain_template
    real(dp)                                          :: curr_gain, ata
    real(dp)                                          :: curr_sigma

    allocate(d_only_wn(size(s_sl)))
    allocate(gain_template(size(s_sl)))

    d_only_wn     = self%scans(scan_id)%d(det)%tod - s_sl - n_corr
    gain_template = s_sky + s_orb
    ata           = sum(mask*gain_template**2)
    curr_gain     = sum(mask * d_only_wn * gain_template) / ata            
    curr_sigma    = self%scans(scan_id)%d(det)%sigma0 / sqrt(ata)  
    self%scans(scan_id)%d(det)%gain       = curr_gain
    self%scans(scan_id)%d(det)%gain_sigma = curr_sigma

    deallocate(d_only_wn,gain_template)

  end subroutine sample_gain


  !compute the cleaned TOD from the computed TOD components
  subroutine compute_cleaned_tod(self, ntod, scan_num, s_orb, s_sl, n_corr, d_calib)
    implicit none
    class(comm_LFI_tod),               intent(in)    :: self
    integer(i4b),                      intent(in)    :: ntod, scan_num
    real(sp),          dimension(:,:), intent(in)    :: n_corr, s_sl, s_orb
    real(sp),          dimension(:,:), intent(out)   :: d_calib

    integer(i4b) :: i

    !cleaned calibrated data = (rawTOD - corrNoise)/gain - orbitalDipole - sideLobes
    do i=1, self%ndet
      d_calib(:,i) = (self%scans(scan_num)%d(i)%tod - n_corr(:,i))/ self%scans(scan_num)%d(i)%gain - s_orb(:,i) - s_sl(:,i)
    end do
  end subroutine compute_cleaned_tod


  ! Compute chisquare
  subroutine compute_chisq(self, scan, det, mask, s_sky, s_sl, s_orb, n_corr)
    implicit none
    class(comm_LFI_tod),             intent(inout)  :: self
    integer(i4b),                    intent(in)     :: scan, det
    real(sp),          dimension(:), intent(in)     :: mask, s_sky, s_sl, s_orb, n_corr

    real(dp) :: chisq    
    integer(i4b) :: i, n

    chisq = 0.d0
    n     = 0
    do i = 1, self%scans(scan)%ntod
       if (mask(i) < 0.5) cycle
       chisq = chisq + (self%scans(scan)%d(det)%tod(i) - &
         & (self%scans(scan)%d(det)%gain * (s_sky(i) + s_sl(i) + s_orb(i)) + &
         & n_corr(i)))**2/self%scans(scan)%d(det)%sigma0**2
       n = n+1
       if (.false. .and. self%myid == 0) write(*,*) i, chisq/n, (self%scans(scan)%d(det)%tod(i) - &
         & (self%scans(scan)%d(det)%gain * (s_sky(i) + s_sl(i) + s_orb(i)) + &
         & n_corr(i)))**2/self%scans(scan)%d(det)%sigma0**2
    end do
    self%scans(scan)%d(det)%chisq = (chisq - n) / sqrt(2.d0*n)

  end subroutine compute_chisq

  ! Sample noise psd
  subroutine sample_noise_psd(self, scan, det, mask, s_sky, s_sl, s_orb, n_corr)
    implicit none
    class(comm_LFI_tod),             intent(inout)  :: self
    integer(i4b),                    intent(in)     :: scan, det
    real(sp),          dimension(:), intent(in)     :: mask, s_sky, s_sl, s_orb, n_corr
    
    integer(i4b) :: i, n
    real(dp)     :: s, res

    s = 0.d0
    n = 0
    do i = 1, self%scans(scan)%ntod-1
       if (any(mask(i:i+1) < 0.5)) cycle
       res = (self%scans(scan)%d(det)%tod(i) - &
         & (self%scans(scan)%d(det)%gain * (s_sky(i) + s_sl(i) + s_orb(i)) + &
         & n_corr(i)) - &
         & (self%scans(scan)%d(det)%tod(i+1) - &
         & (self%scans(scan)%d(det)%gain * (s_sky(i+1) + s_sl(i+1) + s_orb(i+1)) + &
         & n_corr(i+1))))/sqrt(2.)
       s = s + res**2
       n = n + 1
    end do

    self%scans(scan)%d(det)%sigma0 = sqrt(s/(n-1))

  end subroutine sample_noise_psd

  ! Compute map with white noise assumption from correlated noise 
  ! corrected and calibrated data, d' = (d-n_corr-n_temp)/gain 
  subroutine compute_binned_map(self, data, A, b, scan, det)
    implicit none
    class(comm_LFI_tod),                      intent(in)    :: self
    integer(i4b),                             intent(in)    :: scan, det
    real(sp),            dimension(:,:),      intent(in)    :: data
    real(dp),            dimension(0:,1:,1:), intent(inout) :: A
    real(dp),            dimension(0:,1:),    intent(inout) :: b

    integer(i4b) :: i, j, t, pix
    real(dp)     :: psi, cos_psi, sin_psi, inv_sigmasq


    inv_sigmasq = 1.d0 / self%scans(scan)%d(det)%sigma0**2 *1d12
    do t = 1, self%scans(scan)%ntod

       if (iand(self%scans(scan)%d(det)%flag(t),6111248) .ne. 0) cycle

       pix     = self%scans(scan)%d(det)%pix(t)
       psi     = self%scans(scan)%d(det)%psi(t) + self%scans(scan)%d(det)%polang * DEG2RAD
       cos_psi = cos(2.d0*psi)
       sin_psi = sin(2.d0*psi)
       
       A(1,1,pix) = A(1,1,pix) + 1.d0            * inv_sigmasq
       A(1,2,pix) = A(1,2,pix) + cos_psi         * inv_sigmasq
       A(1,3,pix) = A(1,3,pix) + sin_psi         * inv_sigmasq
       A(2,2,pix) = A(2,2,pix) + cos_psi**2      * inv_sigmasq
       A(2,3,pix) = A(2,3,pix) + cos_psi*sin_psi * inv_sigmasq
       A(3,3,pix) = A(3,3,pix) + sin_psi**2      * inv_sigmasq

       b(1,pix) = b(1,pix) + data(t,det)           * inv_sigmasq
       b(2,pix) = b(2,pix) + data(t,det) * cos_psi * inv_sigmasq
       b(3,pix) = b(3,pix) + data(t,det) * sin_psi * inv_sigmasq

    end do

  end subroutine compute_binned_map

  subroutine finalize_binned_map(A, b, map, rms)
    implicit none
    real(dp),        dimension(1:,1:,0:), intent(in)    :: A
    real(dp),        dimension(1:,0:),    intent(in)    :: b
    class(comm_map),                      intent(inout) :: map, rms

    integer(i4b) :: i, j, k, npix, nmaps, np, ierr
    real(dp), allocatable, dimension(:,:)   :: A_inv, b_tot
    real(dp), allocatable, dimension(:,:,:) :: A_tot
    integer(i4b), allocatable, dimension(:) :: pix

    npix  = size(map%map, dim=1)
    nmaps = size(map%map, dim=2)

    ! Collect contributions from all cores
    allocate(A_tot(nmaps,nmaps,0:map%info%np-1), b_tot(nmaps,0:map%info%np-1))
    do i = 0, map%info%nprocs-1
       if (map%info%myid == i) np = map%info%np
       call mpi_bcast(np, 1,  MPI_INTEGER, i, map%info%comm, ierr)
       allocate(pix(np))
       if (map%info%myid == i) pix = map%info%pix
       call mpi_bcast(pix, np,  MPI_INTEGER, i, map%info%comm, ierr)
       do j =1, nmaps
          do k = j, nmaps
             call mpi_reduce(A(j,k,pix), A_tot(j,k,:), np, MPI_DOUBLE_PRECISION, MPI_SUM, i, map%info%comm, ierr)
          end do
       end do
       call mpi_reduce(b(:,pix),   b_tot, np*nmaps,    MPI_DOUBLE_PRECISION, MPI_SUM, i, map%info%comm, ierr)
       deallocate(pix)
    end do

    ! Solve for local map and rms
    allocate(A_inv(nmaps,nmaps))
    map%map = 0.d0
    rms%map = 0.d0
    do i = 0, map%info%np-1
       if (all(b_tot(:,i) == 0.d0)) cycle

       ! rms
       do j = 1, nmaps
          do k = j, nmaps
             A_inv(j,k) = A_tot(j,k,i)
             A_inv(k,j) = A_tot(j,k,i)
          end do
       end do
       call invert_singular_matrix(A_inv, 1d-12)
       do j = 1, nmaps
          rms%map(i,j) = sqrt(A_inv(j,j))
       end do

       ! map
       map%map(i,:) = matmul(A_inv,b_tot(:,i)) * 1.d6 ! uK

    end do

    deallocate(A_inv,A_tot,b_tot)

  end subroutine finalize_binned_map

  



  !construct a sidelobe template in the time domain
  subroutine construct_sl_template(self)
    implicit none
    class(comm_LFI_tod),               intent(in)    :: self

  end subroutine construct_sl_template


end module comm_tod_LFI_mod
