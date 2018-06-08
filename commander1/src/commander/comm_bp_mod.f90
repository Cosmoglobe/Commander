module comm_bp_mod
  use healpix_types
  use comm_utils
  use sort_utils
  implicit none 

  real(dp), parameter :: k_B      = 1.3806503d-23
  real(dp), parameter :: h        = 1.0545726691251021d-34 * 2.d0*pi !6.626068d-34
  real(dp), parameter :: c        = 2.99792458d8
  real(dp)            :: T_CMB    = 2.7255d0
  character(len=128)  :: bp_model = 'additive_shift'

  interface sz_thermo
     module procedure sz_thermo_single, sz_thermo_array
  end interface

  type bandinfo
     character(len=14) :: label, id, unit
     integer(i4b)      :: n
     integer(i4b)      :: comp_calib
     real(dp)          :: nu_c, gain, gain_rms, delta, delta_rms
     real(dp)          :: a2t, f2t, co2t, a2sz
     real(dp), allocatable, dimension(:) :: nu0, nu, tau0, tau
  end type bandinfo

  real(dp)                                      :: ind_iras                                     
  integer(i4b)                                  :: numband
  type(bandinfo),     allocatable, dimension(:) :: bp
  integer(i4b),       allocatable, dimension(:) :: i2f

  character(len=256), private  :: MJySr_convention

contains


  subroutine initialize_bp_mod(myid, chain, comm_chain, paramfile, handle)
    implicit none
    
    integer(i4b),     intent(in)    :: myid, chain, comm_chain
    character(len=*), intent(in)    :: paramfile
    type(planck_rng), intent(inout) :: handle

    integer(i4b)        :: i, j, q, unit, myid_chain
    logical(lgt)        :: exist, apply_bp_corr, apply_gain_corr
    real(dp)            :: threshold, gain_init_rms, bp_init_rms, ierr
    character(len=2)    :: i_text
    character(len=256)  :: filename, chaindir, gain_init_file, bp_init_file, MJySr_convention
    character(len=5000) :: line
    real(dp), allocatable, dimension(:) :: freq, a2t

    call mpi_comm_rank(comm_chain, myid_chain, ierr)

    unit = getlun()
    call get_parameter(paramfile, 'NUMBAND',                par_int=numband)    
    call get_parameter(paramfile, 'T_CMB',                  par_dp=T_CMB)
    call get_parameter(paramfile, 'CHAIN_DIRECTORY',        par_string=chaindir)
    call get_parameter(paramfile, 'GAIN_INIT',              par_string=gain_init_file)
    call get_parameter(paramfile, 'BANDPASS_INIT',          par_string=bp_init_file)
    call get_parameter(paramfile, 'MJYSR_CONVENTION',       par_string=MJysr_convention)
    call get_parameter(paramfile, 'APPLY_BP_CORRECTIONS',   par_lgt=apply_bp_corr)
    call get_parameter(paramfile, 'APPLY_GAIN_CORRECTIONS', par_lgt=apply_gain_corr)
    call get_parameter(paramfile, 'GAIN_INIT_RMS',          par_dp=gain_init_rms)
    call get_parameter(paramfile, 'BP_INIT_RMS',            par_dp=bp_init_rms)
    if (trim(MJysr_convention) == 'PSM') then
       ind_iras = 0.d0
    else if (trim(MJysr_convention) == 'IRAS') then
       ind_iras = 1.d0
    else
       write(*,*) 'Unsupported MJy/sr convention = ', trim(MJysr_convention)
       stop
    end if

    allocate(bp(numband))
    do i = 1, numband
       call int2string(i, i_text)
       call get_parameter(paramfile, 'BANDPASS_TYPE'        // i_text, par_string=bp(i)%id)
       call get_parameter(paramfile, 'FREQ_C'               // i_text, par_dp=bp(i)%nu_c)
       bp(i)%nu_c = bp(i)%nu_c * 1.d9
       call get_parameter(paramfile, 'FREQ_UNIT'            // i_text, par_string=bp(i)%unit)       
       call get_parameter(paramfile, 'FREQ_LABEL'           // i_text, par_string=bp(i)%label)       
       call get_parameter(paramfile, 'GAIN_CALIB_COMPONENT' // i_text, par_int=bp(i)%comp_calib)
       call get_parameter(paramfile, 'GAIN_PRIOR_RMS'       // i_text, par_dp=bp(i)%gain_rms)
       call get_parameter(paramfile, 'A2T'                  // i_text, par_dp=bp(i)%a2t)
       call get_parameter(paramfile, 'F2T'                  // i_text, par_dp=bp(i)%f2t)
       call get_parameter(paramfile, 'CO2T'                 // i_text, par_dp=bp(i)%co2t)
       call get_parameter(paramfile, 'A2SZ'                 // i_text, par_dp=bp(i)%a2sz)
       bp(i)%delta     = 0.d0
       if (apply_bp_corr) then
          call get_parameter(paramfile, 'BP_RMS'         // i_text, par_dp=bp(i)%delta_rms)
       else
          bp(i)%delta_rms = 0.d0
       end if
       if (trim(bp(i)%id) == 'delta') then
          filename = ''
       else if (trim(bp(i)%id) == 'LFI') then
          call get_parameter(paramfile, 'BANDPASS'           // i_text, par_string=filename)
          threshold = 0.d0
       else if (trim(bp(i)%id) == 'HFI_cmb' .or. trim(bp(i)%id) == 'PSM_LFI') then
          call get_parameter(paramfile, 'BANDPASS'           // i_text, par_string=filename)
          threshold = 1.d-7
       else if (trim(bp(i)%id) == 'HFI_submm') then
          call get_parameter(paramfile, 'BANDPASS'           // i_text, par_string=filename)
          threshold = 1.d-5
       else if (trim(bp(i)%id) == 'WMAP') then
          call get_parameter(paramfile, 'BANDPASS'           // i_text, par_string=filename)
          threshold = 0.d0
       else if (trim(bp(i)%id) == 'DIRBE') then
          call get_parameter(paramfile, 'BANDPASS'           // i_text, par_string=filename)
          threshold = 0.d0
       else 
          write(*,*) 'Error -- bandpass type not supported: ', trim(bp(i)%id)
          write(*,*) '         Supported types are {delta, LFI, HFI, WMAP, DIRBE}'
          stop
       end if

       if (trim(filename) /= '') then
          call read_bandpass(filename, threshold, bp(i)%n, bp(i)%nu0, bp(i)%tau0)
          allocate(bp(i)%nu(bp(i)%n), bp(i)%tau(bp(i)%n))
       else
          allocate(bp(i)%nu0(1),bp(i)%tau0(1), bp(i)%nu(1), bp(i)%tau(1))
          bp(i)%n       = 1
          bp(i)%nu0(1)  = bp(i)%nu_c
          bp(i)%tau0(1) = 1.d0
       end if

!!$       allocate(a2t(bp(i)%n))
!!$       call compute_ant2thermo(bp(i)%nu, a2t)
!!$       if (myid == 0) write(*,*) bp(i)%nu
!!$       if (myid == 0) write(*,*) a2t
!!$       call mpi_finalize(ierr)
!!$       stop

    end do

    ! Sort bands according to frequencies
    allocate(i2f(numband), freq(numband))
    do i = 1, numband
       i2f(i)  = i
       freq(i) = bp(i)%nu_c
    end do
    call QuickSort(i2f,freq)
    deallocate(freq)

    if (apply_gain_corr) then
       ! Initialize gains
       inquire(file=trim(gain_init_file), exist=exist)
       if (exist) then
          open(unit,file=trim(gain_init_file), recl=5000)
          do while (.true.)
             read(unit,'(a)',end=87) line
             if (line(1:1) == '#') cycle
             read(line,*) j, bp%gain
          end do
87        close(unit)
       else
          bp%gain = 1.d0
       end if
    else
       bp%gain = 1.d0
    end if

    if (apply_bp_corr) then
       ! Initialize bandpasses
       inquire(file=trim(bp_init_file), exist=exist)
       if (exist) then
          open(unit,file=trim(bp_init_file), recl=5000)
          do while (.true.)
             read(unit,'(a)',end=88) line
             if (line(1:1) == '#') cycle
             read(line,*) j, bp%delta
          end do
88        close(unit)
       end if
    end if
    
    if (myid_chain == 0) then
       do i = 1, numband
          if (bp(i)%gain_rms > 0.d0) &
               & bp(i)%gain  = bp(i)%gain  + gain_init_rms  * rand_gauss(handle)
          if (bp(i)%delta_rms > 0.d0) &
               & bp(i)%delta = bp(i)%delta + bp_init_rms    * rand_gauss(handle)
       end do
    end if
    call mpi_bcast(bp%gain,  numband, MPI_DOUBLE_PRECISION, 0, comm_chain, ierr)
    call mpi_bcast(bp%delta, numband, MPI_DOUBLE_PRECISION, 0, comm_chain, ierr)
    call update_tau(bp%delta)

    if (myid == 0) then
       open(unit,file=trim(chaindir)//'/unit_conversions.dat',recl=1024)
       write(unit,*) '# Band        Convention           Nu_c (GHz)      a2t [K_cmb/K_RJ]' // &
            & '      t2f [MJy/K_cmb]   co2t [uK_cmb / (K km/s)]   a2sz [y_sz/K_RJ]'
       do i = 1, numband
          q = i2f(i)
          write(unit,fmt='(a,a,a,a,f16.3,4e20.8)') bp(q)%label, ' ', bp(q)%id, ' ', &
               & bp(q)%nu_c/1.d9, bp(q)%a2t, 1.d0/bp(q)%f2t*1e6, bp(q)%co2t, bp(q)%a2sz * 1.d6
       end do
       close(unit)
    end if

  end subroutine initialize_bp_mod


  subroutine update_tau(delta, overwrite)
    implicit none

    real(dp),     dimension(:), intent(in)           :: delta
    logical(lgt),               intent(in), optional :: overwrite

    integer(i4b) :: i, j, n
    logical(lgt) :: overwrite_
    real(dp), allocatable, dimension(:) :: a2t, bnu_prime, bnu_prime_RJ, sz
    real(dp), allocatable, dimension(:) :: buffer

    overwrite_ = .true.; if (present(overwrite)) overwrite_ = overwrite
    allocate(buffer(numband))
    buffer   = bp%delta
    bp%delta = delta

    do j = 1, numband

       n = bp(j)%n

       if (trim(bp_model) == 'linear_tilt') then
          
          ! Linear tilt model
          bp(j)%nu = bp(j)%nu0
          do i = 1, n
             bp(j)%tau(i) = bp(j)%tau0(i) * (1.d0 + 1d9*delta(j) * (bp(j)%nu(i)-bp(j)%nu_c))
          end do
          
       else if (trim(bp_model) == 'powlaw_tilt') then
          
          ! Linear tilt model
          bp(j)%nu = bp(j)%nu0
          do i = 1, n
             bp(j)%tau(i) = bp(j)%tau0(i) * (bp(j)%nu(i)/bp(j)%nu_c)**delta(j)
          end do

       else if (trim(bp_model) == 'additive_shift') then

          ! Additive frequency shift
          bp(j)%tau = bp(j)%tau0
          do i = 1, n
             bp(j)%nu(i) = bp(j)%nu0(i) + 1d9*delta(j)
             if (bp(j)%nu(i) <= 0.d0) bp(j)%tau(i) = 0.d0
          end do
!!$          write(*,*) bp(j)%nu
!!$          write(*,*) bp(j)%tau
!!$          call mpi_finalize(i)
!!$          stop
          
       else if (trim(bp_model) == 'multiplicative_shift') then
          
          ! Multiplicative frequency shift
          bp(j)%tau = bp(j)%tau0
          do i = 1, n
             bp(j)%nu(i) = bp(j)%nu0(i) * (1.d0 + 1d9*delta(j))
          end do
          
       else
          write(*,*) 'Error: Unknown bandpass model = ', trim(bp_model)
          stop
       end if

       
       ! Compute unit conversion factors
       if (trim(bp(j)%id) == 'delta') then

          bp(j)%nu      = bp(j)%nu0
          if (bp(j)%a2t  <= 0.d0) bp(j)%a2t  = compute_ant2thermo_single(bp(j)%nu_c)
          if (bp(j)%a2sz <= 0.d0) bp(j)%a2sz = 2.d0*bp(j)%nu_c**2*k_b/c**2 / &
               & (compute_bnu_prime_single(bp(j)%nu_c) * sz_thermo(bp(j)%nu_c)) * 1.d-6
          if (bp(j)%f2t  <= 0.d0) bp(j)%f2t  = 1.d0 / compute_bnu_prime_single(bp(j)%nu_c) * 1.d-14
          bp(j)%tau     = bp(j)%tau0 * bp(j)%a2t
          !write(*,*) 'nu', bp(j)%nu, bp(j)%nu0, bp(j)%nu_c
          !write(*,*) 'tau', bp(j)%tau, bp(j)%tau0, bp(j)%a2t

       else if (trim(bp(j)%id) == 'WMAP') then
          
          ! See Appendix E of Bennett et al. (2013) for details
          allocate(a2t(bp(j)%n), bnu_prime(bp(j)%n), bnu_prime_RJ(bp(j)%n), sz(bp(j)%n))
          do i = 1, bp(j)%n
             a2t(i)          = compute_ant2thermo_single(bp(j)%nu(i))
             bnu_prime(i)    = compute_bnu_prime_single(bp(j)%nu(i))
             bnu_prime_RJ(i) = compute_bnu_prime_RJ_single(bp(j)%nu(i))
             sz(i)           = sz_thermo_single(bp(j)%nu(i))
          end do
          !call compute_ant2thermo(bp(j)%nu, a2t)          
          !call compute_bnu_prime(bp(j)%nu, bnu_prime)
          !call compute_bnu_prime_RJ(bp(j)%nu, bnu_prime_RJ)
          !sz = sz_thermo(bp(j)%nu)
          if (bp(j)%a2t  <= 0.d0) bp(j)%a2t  = sum(bp(j)%tau) / sum(bp(j)%tau/a2t)
          if (bp(j)%a2sz <= 0.d0) bp(j)%a2sz = sum(bp(j)%tau) / sum(bp(j)%tau/a2t * sz) * 1.d-6
          if (bp(j)%f2t  <= 0.d0) bp(j)%f2t  = sum(bp(j)%tau/bp(j)%nu**2 * (bp(j)%nu_c/bp(j)%nu)**ind_iras) &
               & * 1.d-14 / sum(bp(j)%tau/bp(j)%nu**2 * bnu_prime)
          bp(j)%tau = bp(j)%tau * a2t 
          deallocate(a2t, bnu_prime, bnu_prime_RJ, sz)

       else if (trim(bp(j)%id) == 'LFI') then

          ! Multiply with normalization factor
          allocate(bnu_prime(bp(j)%n), bnu_prime_RJ(bp(j)%n), a2t(bp(j)%n), sz(bp(j)%n))
          do i = 1, bp(j)%n
             a2t(i)          = compute_ant2thermo_single(bp(j)%nu(i))
             bnu_prime(i)    = compute_bnu_prime_single(bp(j)%nu(i))
             bnu_prime_RJ(i) = compute_bnu_prime_RJ_single(bp(j)%nu(i))
             sz(i)           = sz_thermo_single(bp(j)%nu(i))
          end do
          !call compute_ant2thermo(bp(j)%nu, a2t)          
          !call compute_bnu_prime(bp(j)%nu, bnu_prime)
          !call compute_bnu_prime_RJ(bp(j)%nu, bnu_prime_RJ)
          !sz = sz_thermo(bp(j)%nu)
          if (bp(j)%a2t <= 0.d0) bp(j)%a2t = tsum(bp(j)%nu, bp(j)%tau/bp(j)%nu**2 * bnu_prime_RJ) / &
               & tsum(bp(j)%nu, bp(j)%tau/bp(j)%nu**2 * bnu_prime)
          if (bp(j)%a2sz <= 0.d0) bp(j)%a2sz = tsum(bp(j)%nu, bp(j)%tau/bp(j)%nu**2 * bnu_prime_RJ) / &
               & tsum(bp(j)%nu, bp(j)%tau/bp(j)%nu**2 * bnu_prime * sz) * 1.d-6
          if (bp(j)%f2t <= 0.d0) bp(j)%f2t = tsum(bp(j)%nu, bp(j)%tau/bp(j)%nu**2 * &
               & (bp(j)%nu_c/bp(j)%nu)**ind_iras) * 1.d-14 / tsum(bp(j)%nu, bp(j)%tau/bp(j)%nu**2 * bnu_prime)
          bp(j)%tau = bp(j)%tau / tsum(bp(j)%nu, bp(j)%tau/a2t)
          deallocate(a2t, bnu_prime, bnu_prime_RJ, sz)

       else if (trim(bp(j)%id) == 'HFI_cmb' .or. trim(bp(j)%id) == 'PSM_LFI') then
          
          allocate(bnu_prime(bp(j)%n), bnu_prime_RJ(bp(j)%n), a2t(bp(j)%n), sz(bp(j)%n))
          do i = 1, bp(j)%n
             a2t(i)          = compute_ant2thermo_single(bp(j)%nu(i))
             bnu_prime(i)    = compute_bnu_prime_single(bp(j)%nu(i))
             bnu_prime_RJ(i) = compute_bnu_prime_RJ_single(bp(j)%nu(i))
             sz(i)           = sz_thermo_single(bp(j)%nu(i))
          end do
          !call compute_ant2thermo(bp(j)%nu, a2t)          
          !call compute_bnu_prime(bp(j)%nu, bnu_prime)
          !call compute_bnu_prime_RJ(bp(j)%nu, bnu_prime_RJ)
          !sz = sz_thermo(bp(j)%nu)
          if (bp(j)%a2t <= 0.d0) bp(j)%a2t = tsum(bp(j)%nu, bp(j)%tau * bnu_prime_RJ) / &
               & tsum(bp(j)%nu, bp(j)%tau*bnu_prime)
          if (bp(j)%a2sz <= 0.d0) bp(j)%a2sz = tsum(bp(j)%nu, bp(j)%tau * bnu_prime_RJ) / &
               & tsum(bp(j)%nu, bp(j)%tau*bnu_prime*sz) * 1.d-6
          if (bp(j)%f2t <= 0.d0) bp(j)%f2t = tsum(bp(j)%nu, bp(j)%tau * (bp(j)%nu_c/bp(j)%nu)**ind_iras) &
               & * 1.d-14 / tsum(bp(j)%nu, bp(j)%tau*bnu_prime)
          bp(j)%tau = bp(j)%tau / tsum(bp(j)%nu, bp(j)%tau*bnu_prime)
          deallocate(bnu_prime, bnu_prime_RJ, a2t, sz)

       else if (trim(bp(j)%id) == 'HFI_submm') then

          allocate(bnu_prime(bp(j)%n), bnu_prime_RJ(bp(j)%n), a2t(bp(j)%n), sz(bp(j)%n))
          do i = 1, bp(j)%n
             a2t(i)          = compute_ant2thermo_single(bp(j)%nu(i))
             bnu_prime(i)    = compute_bnu_prime_single(bp(j)%nu(i))
             bnu_prime_RJ(i) = compute_bnu_prime_RJ_single(bp(j)%nu(i))
             sz(i)           = sz_thermo_single(bp(j)%nu(i))
          end do
          !call compute_ant2thermo(bp(j)%nu, a2t)          
          !call compute_bnu_prime(bp(j)%nu, bnu_prime)
          !call compute_bnu_prime_RJ(bp(j)%nu, bnu_prime_RJ)
          !sz = sz_thermo(bp(j)%nu)
          if (bp(j)%a2t <= 0.d0) bp(j)%a2t = tsum(bp(j)%nu, bp(j)%tau * bnu_prime_RJ) / &
               & tsum(bp(j)%nu, bp(j)%tau*bnu_prime)
          if (bp(j)%a2sz <= 0.d0) bp(j)%a2sz = tsum(bp(j)%nu, bp(j)%tau * bnu_prime_RJ) / &
               & tsum(bp(j)%nu, bp(j)%tau*bnu_prime*sz) * 1.d-6
          if (bp(j)%f2t <= 0.d0) bp(j)%f2t = tsum(bp(j)%nu, bp(j)%tau * (bp(j)%nu_c/bp(j)%nu)**ind_iras) * &
               & 1.d-14 / tsum(bp(j)%nu, bp(j)%tau*bnu_prime)
          bp(j)%tau = bp(j)%tau / tsum(bp(j)%nu, bp(j)%tau * (bp(j)%nu_c/bp(j)%nu)**ind_iras) * 1.d14
          deallocate(bnu_prime, bnu_prime_RJ, a2t, sz)

       else if (trim(bp(j)%id) == 'DIRBE') then

          allocate(bnu_prime(bp(j)%n), bnu_prime_RJ(bp(j)%n), a2t(bp(j)%n), sz(bp(j)%n))
          do i = 1, bp(j)%n
             a2t(i)          = compute_ant2thermo_single(bp(j)%nu(i))
             bnu_prime(i)    = compute_bnu_prime_single(bp(j)%nu(i))
             bnu_prime_RJ(i) = compute_bnu_prime_RJ_single(bp(j)%nu(i))
             sz(i)           = sz_thermo_single(bp(j)%nu(i))
          end do
          !call compute_ant2thermo(bp(j)%nu, a2t)          
          !call compute_bnu_prime(bp(j)%nu, bnu_prime)
          !call compute_bnu_prime_RJ(bp(j)%nu, bnu_prime_RJ)
          !sz = sz_thermo(bp(j)%nu)
          if (bp(j)%a2t <= 0.d0) bp(j)%a2t = tsum(bp(j)%nu, bp(j)%tau * bnu_prime_RJ) / &
               & tsum(bp(j)%nu, bp(j)%tau*bnu_prime)
          if (bp(j)%a2sz <= 0.d0) bp(j)%a2sz = tsum(bp(j)%nu, bp(j)%tau * bnu_prime_RJ) / &
               & tsum(bp(j)%nu, bp(j)%tau*bnu_prime*sz) * 1.d-6
          if (bp(j)%f2t <= 0.d0) bp(j)%f2t = tsum(bp(j)%nu, bp(j)%tau * (bp(j)%nu_c/bp(j)%nu)**ind_iras) * &
               & 1.d-14 / tsum(bp(j)%nu, bp(j)%tau*bnu_prime)
          bp(j)%tau = bp(j)%tau / tsum(bp(j)%nu, bp(j)%tau * (bp(j)%nu_c/bp(j)%nu)**ind_iras) * 1.d14
          deallocate(bnu_prime, bnu_prime_RJ, a2t, sz)

       end if

    end do

    ! Revert to old value if not overwrite; use with caution
    if (.not. overwrite_) bp%delta = buffer

  end subroutine update_tau


  function get_bp_avg_spectrum(band, f)
    implicit none

    integer(i4b),               intent(in) :: band
    real(dp),     dimension(:), intent(in) :: f
    real(dp)                               :: get_bp_avg_spectrum
    real(dp),     dimension(:)             :: hm(1)

    if (trim(bp(band)%id) == 'delta') then
       ! Return single value
       hm = f * bp(band)%tau 
       get_bp_avg_spectrum = hm(1) 
    else if (trim(bp(band)%id) == 'LFI') then
       get_bp_avg_spectrum = tsum(bp(band)%nu, bp(band)%tau * f)
    else if (trim(bp(band)%id) == 'HFI_cmb' .or. trim(bp(band)%id) == 'PSM_LFI') then
       get_bp_avg_spectrum = tsum(bp(band)%nu, bp(band)%tau * 2.d0*k_B*bp(band)%nu**2/c**2 * f)
    else if (trim(bp(band)%id) == 'HFI_submm') then
       get_bp_avg_spectrum = tsum(bp(band)%nu, bp(band)%tau * 2.d0*k_B*bp(band)%nu**2/c**2 * f)
    else if (trim(bp(band)%id) == 'DIRBE') then
       ! Assume same as HFI_submm 
       get_bp_avg_spectrum = tsum(bp(band)%nu, bp(band)%tau * 2.d0*k_B*bp(band)%nu**2/c**2 * f)
    else if (trim(bp(band)%id) == 'WMAP') then
       get_bp_avg_spectrum = sum(bp(band)%tau * f)
    end if

  end function get_bp_avg_spectrum

  function get_bp_line_ant(band, nu)
    implicit none

    integer(i4b), intent(in) :: band
    real(dp),     intent(in) :: nu
    real(dp)                 :: get_bp_line_ant

    integer(i4b) :: i
    real(dp)     :: norm, x, tau
    real(dp), allocatable, dimension(:) :: b, bnu_prime_RJ, a2t

    if (nu < bp(band)%nu(1) .or. nu > bp(band)%nu(bp(band)%n)) then
       get_bp_line_ant = 0.d0
       return
    end if

    ! Read off correct bandpass coefficient; use linear interpolation
    i = 1
    do while (bp(band)%nu(i) < nu .and. i < bp(band)%n)
       i = i+1
    end do
    x   = (bp(band)%nu(i)-nu) / (bp(band)%nu(i)-bp(band)%nu(i-1))
    tau = x * bp(band)%tau(i-1) + (1.d0-x)*bp(band)%tau(i)


    if (trim(bp(band)%id) == 'delta') then
       if (nu /= bp(band)%nu_c) then
          get_bp_line_ant = 0.d0
       else
          get_bp_line_ant = nu/c 
       end if

    else if (trim(bp(band)%id) == 'WMAP') then

       ! See Appendix E of Bennett et al. (2013) for details
       allocate(a2t(bp(band)%n))
       call compute_ant2thermo(bp(band)%nu, a2t)          
       get_bp_line_ant = tau/(bp(band)%nu(2)-bp(band)%nu(1)) * nu/c / sum(bp(band)%tau)
       deallocate(a2t)

    else if (trim(bp(band)%id) == 'LFI') then

       ! Multiply with normalization factor
       allocate(bnu_prime_RJ(bp(band)%n))       
       call compute_bnu_prime_RJ(bp(band)%nu, bnu_prime_RJ)
       get_bp_line_ant = tau/nu**2 * nu/c * compute_bnu_prime_RJ_single(nu) / &
            & tsum(bp(band)%nu, bp(band)%tau/bp(band)%nu**2 * bnu_prime_RJ)
       deallocate(bnu_prime_RJ)

    else if (trim(bp(band)%id) == 'HFI_cmb' .or. trim(bp(band)%id) == 'PSM_LFI') then
          
       allocate(bnu_prime_RJ(bp(band)%n))
       call compute_bnu_prime_RJ(bp(band)%nu, bnu_prime_RJ)
       get_bp_line_ant = tau * nu/c * compute_bnu_prime_RJ_single(nu) / &
            & tsum(bp(band)%nu, bp(band)%tau*bnu_prime_RJ)
       deallocate(bnu_prime_RJ)

    else if (trim(bp(band)%id) == 'HFI_submm') then
       
       get_bp_line_ant = tau * nu/c * compute_bnu_prime_RJ_single(nu) / &
            & (tsum(bp(band)%nu, bp(band)%tau * (bp(band)%nu_c/bp(band)%nu)**ind_iras) * 1.d-14)
       
    else 
       write(*,*) 'Bandpass type not yet supported for spectral lines: ', trim(bp(band)%id)
       stop
    end if

    get_bp_line_ant = get_bp_line_ant * 1.d9 ! Convert to uK_ant / (K_ant km/s)

  end function get_bp_line_ant

  function ant2unit(unit_out, nu_ref, band_ref)
    implicit none

    character(len=*), intent(in)           :: unit_out
    real(dp),         intent(in), optional :: nu_ref    
    integer(i4b),     intent(in), optional :: band_ref
    real(dp)                               :: ant2unit

    if (trim(unit_out) == 'uK_cmb') then
       ant2unit = compute_ant2thermo_single(nu_ref)
    else if (trim(unit_out) == 'MJy/sr') then
       ant2unit = compute_bnu_prime_RJ_single(nu_ref) * 1e14
    else if (trim(unit_out) == 'K km/s') then
       !ant2unit = bp(band_ref)%a2t / bp(band_ref)%co2t
       !ant2unit = 1.d0 / (nu_ref/c)
       ant2unit = 1.d0 / get_bp_line_ant(band_ref, nu_ref)
    else if (trim(unit_out) == 'y_SZ') then
       ant2unit = 2.d0*nu_ref**2*k_b/c**2 / &
               & (compute_bnu_prime_single(nu_ref) * sz_thermo(nu_ref))
    else if (trim(unit_out) == 'uK_ant') then
       ant2unit = 1.d0
    else
       write(*,*) 'Unsupported unit: ', trim(unit_out)
       stop
    end if

  end function ant2unit


  function ant2data(band)
    implicit none

    integer(i4b), intent(in) :: band
    real(dp)                 :: ant2data

    if (trim(bp(band)%unit) == 'uK_cmb') then
       ant2data = bp(band)%a2t
    else if (trim(bp(band)%unit) == 'MJy/sr') then
       ant2data = bp(band)%a2t / bp(band)%f2t
    else if (trim(bp(band)%unit) == 'K km/s') then
       ant2data = bp(band)%a2t / bp(band)%co2t
    else if (trim(bp(band)%unit) == 'y_SZ') then
       ant2data = bp(band)%a2sz
    else if (trim(bp(band)%unit) == 'uK_ant') then
       ant2data = 1.d0
    end if
    
  end function ant2data

  function sz_thermo_single(nu)
    implicit none

    real(dp), intent(in) :: nu
    real(dp)             :: sz_thermo_single

    real(dp) :: x

    x = h*nu/(k_b*T_cmb)

    sz_thermo_single = T_cmb * (x*(exp(x)+1.d0)/(exp(x)-1.d0)-4.d0)
    
  end function sz_thermo_single

  function sz_thermo_array(nu) result (y)
    implicit none

    real(dp),              dimension(:), intent(in)  :: nu
    real(dp), allocatable, dimension(:)              :: y

    real(dp), allocatable, dimension(:) :: x

    allocate(x(size(nu)), y(size(nu)))
    x = h*nu/(k_b*T_cmb)
    y = T_cmb * (x*(exp(x)+1.d0)/(exp(x)-1.d0)-4.d0)
    deallocate(x)

  end function sz_thermo_array
  

  subroutine compute_ant2thermo(frequencies, ant2thermo)
    implicit none

    real(dp), dimension(:), intent(in)  :: frequencies
    real(dp), dimension(:), intent(out) :: ant2thermo

    integer(i4b) :: i, numband
    real(dp)     :: x

    numband = size(frequencies)

    do i = 1, numband
       x = h*frequencies(i) / (k_B*T_CMB)
       ant2thermo(i) = (exp(x)-1.d0)**2 / (x**2 * exp(x))
    end do
    
  end subroutine compute_ant2thermo

  function compute_ant2thermo_single(frequency)
    implicit none

    real(dp), intent(in)  :: frequency
    real(dp)              :: compute_ant2thermo_single

    real(dp)     :: x

    x = h*frequency / (k_B*T_CMB)
    compute_ant2thermo_single = (exp(x)-1.d0)**2 / (x**2 * exp(x))
    
  end function compute_ant2thermo_single

  subroutine compute_bnu_prime(nu, bnu_prime)
    implicit none

    real(dp), dimension(:), intent(in)  :: nu
    real(dp), dimension(:), intent(out) :: bnu_prime

    integer(i4b) :: i
    real(dp)     :: x

    do i = 1, size(nu)
       x = h*nu(i) / (k_B*T_CMB)
       bnu_prime(i) = (2.d0*h*nu(i)**3/(c**2*(exp(x)-1.d0))) * &
            & (exp(x) / (exp(x)-1.d0)) * h*nu(i)/(k_B*T_CMB**2)
    end do
    
  end subroutine compute_bnu_prime

  function compute_bnu_prime_single(nu)
    implicit none

    real(dp), intent(in)  :: nu
    real(dp)              :: compute_bnu_prime_single

    integer(i4b) :: i
    real(dp)     :: x

    x = h*nu / (k_B*T_CMB)
    compute_bnu_prime_single = (2.d0*h*nu**3/(c**2*(exp(x)-1.d0))) * &
         & (exp(x) / (exp(x)-1.d0)) * h*nu/(k_B*T_CMB**2)
    
  end function compute_bnu_prime_single

  subroutine compute_bnu_prime_RJ(nu, bnu_prime_RJ)
    implicit none

    real(dp), dimension(:), intent(in)  :: nu
    real(dp), dimension(:), intent(out) :: bnu_prime_RJ

    bnu_prime_RJ = 2.d0*k_B*nu**2/c**2
    
  end subroutine compute_bnu_prime_RJ

  function compute_bnu_prime_RJ_single(nu)
    implicit none

    real(dp), intent(in)  :: nu
    real(dp)              :: compute_bnu_prime_RJ_single

    compute_bnu_prime_RJ_single = 2.d0*k_B*nu**2/c**2
    
  end function compute_bnu_prime_RJ_single

  function dBdnu(nu)
    implicit none

    real(dp), intent(in)  :: nu
    real(dp)              :: dBdnu

    real(dp)     :: x

    x = h*nu / (k_B*T_CMB)
    dBdnu = (2*h*nu**3/(c**2 * (exp(x)-1.d0))) *  (exp(x) / (exp(x)-1.d0)) * (h*nu/(k_b*T_CMB**2))
    
  end function dBdnu

  function dB_rj_dnu(nu)
    implicit none

    real(dp), intent(in)  :: nu
    real(dp)              :: dB_rj_dnu

    dB_rj_dnu = 2*nu**2*k_b / c**2
    
  end function dB_rj_dnu

  ! Routine for reading bandpass files
  subroutine read_bandpass(filename, threshold, n, nu, tau)
    implicit none

    character(len=*),                            intent(in)  :: filename
    real(dp),                                    intent(in)  :: threshold
    integer(i4b),                                intent(out) :: n
    real(dp),         allocatable, dimension(:), intent(out) :: nu, tau

    integer(i4b)        :: i, j, unit, first, last, m
    logical(lgt)        :: exist
    character(len=128)  :: string
    real(dp), allocatable, dimension(:) :: x, y

    unit = getlun()
    
    inquire(file=trim(filename), exist=exist)
    if (.not. exist) then
       write(*,*) 'Bandpass file does not exist = ', trim(filename)
       stop
    end if

    ! Find the number of entries
    m = 0
    open(unit, file=trim(filename))
    do while (.true.)
       read(unit,*,end=1) string
       if (string(1:1)=='#') cycle
       m = m + 1
    end do
1   close(unit)

    if (m == 0) then
       write(*,*) 'No valid data entries in spectrum file ', trim(filename)
       stop
    end if

    allocate(x(m), y(m))
    m = 0
    open(unit, file=trim(filename))
    do while (.true.)
       read(unit,fmt='(a)',end=2) string
       if (string(1:1)=='#') cycle
       m = m+1
       read(string,*) x(m), y(m)

       ! Drop double entries
       if (m > 1) then
          if (x(m) == x(m-1)) m = m-1
       end if
    end do
2   close(unit)

    x = x * 1.d9 ! Convert from GHz to Hz

    first = 1
    last  = m
    if (threshold > 0.d0) then
       do while (y(first) < threshold*maxval(y(1:m)))
          first = first+1
       end do
       do while (y(last) < threshold*maxval(y(1:m)))
          last = last-1
       end do
    end if
    
    n = last-first+1
    allocate(nu(n), tau(n))
    nu  = x(first:last)
    tau = y(first:last)

    deallocate(x, y)

  end subroutine read_bandpass



end module comm_bp_mod
