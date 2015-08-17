program comm_like_sampler
  use comm_lowl_mod
  use comm_like_model_mod
  use comm_proc_utils
  use sort_utils
  use pix_tools
  implicit none

  include 'mpif.h'

  integer(i4b)       :: i, j, k, l, m, ind, nsamp, seed, seed_i, unit, status, ierr
  integer(i4b)       :: nside_dir, marg_freq, burnin_mcmc
  integer(i4b)       :: err, myid, n_proc, root
  character(len=3)   :: i_s
  logical(lgt)       :: analyze_simulation
  real(dp)           :: dx, t1, t2
  character(len=256) :: parfile, likelihood_file, prefix, filename, outpropfile, operation
  type(planck_rng)   :: rng_handle, sim_handle
  real(dp), allocatable, dimension(:,:) :: samples

  ! Sampling parameters
  real(dp)                              :: x_init(3), lnL_prev, lnL
  integer(i4b)                          :: p_current, prop_output_freq
  real(dp), allocatable, dimension(:)   :: x_current, x_n, S_n, eta
  real(dp), allocatable, dimension(:,:) :: L_prop, sqrt_S

  call mpi_init(err)
  call mpi_comm_rank(mpi_comm_world, myid, err)
  call mpi_comm_size(mpi_comm_world, n_proc, err)
  root = 0


  ! Initialize model module and likelihood file
  call getarg(1,parfile)
  if(myid==root) write(*,*) 'Initializing with parameter file = ', trim(parfile)

  call initialize_model_mod(parfile)
  call comm_get_parameter(parfile, 'LIKELIHOOD_PARAMETER_FILE', par_string=likelihood_file)
  call comm_lowl_initialize_object(likelihood_file)

  ! Read general parameters
  call comm_get_parameter(parfile, 'OPERATION',                    par_string=operation)
  call comm_get_parameter(parfile, 'NUM_SAMPLES',                  par_int=nsamp)
  call comm_get_parameter(parfile, 'SEED',                         par_int=seed)
  call comm_get_parameter(parfile, 'OUT_PREFIX',                   par_string=prefix)
  call comm_get_parameter(parfile, 'OUTPUT_PROPOSAL_MATRIX',       par_string=outpropfile)
  call comm_get_parameter(parfile, 'PROP_MATRIX_OUTPUT_FREQUENCY', par_int=prop_output_freq)
  call comm_get_parameter(parfile, 'ANALYZE_INTERNAL_SIMULATION',  par_lgt=analyze_simulation)
  call rand_init(sim_handle, seed+1)
  if (myid == root) then
     call rand_init(rng_handle, seed)
     do i = 1, n_proc-1
        seed_i = rand_uni(rng_handle)*1000000
        call mpi_send(seed_i, 1, mpi_integer, i, 0, mpi_comm_world, err)
     end do
  else
     call mpi_recv(seed, 1, mpi_integer, root, 0, mpi_comm_world, status, err)
     call rand_init(rng_handle, seed)
  end if

  unit = comm_getlun()

  ! Initialize chain and proposal density
  allocate(samples(0:nsamp,npar), x_current(npar), L_prop(npar,npar), eta(npar))
  call initialize_model_parameters(rng_handle, samples(0,:), L_prop)

  if (analyze_simulation) then
     ! Set up simulation -> phi, theta, psi, q, n, alpha
     samples(0,1) = 220.d0*pi/180.d0 
     samples(0,2) = 105.d0*pi/180.d0 
     samples(0,3) = 0.d0
     samples(0,4) = 1.d0
     samples(0,5) = 0.d0
     samples(0,6) = 0.2d0
  
     allocate(sqrt_S(n_h,n_h))
     call get_sqrt_S(samples(0,:), sqrt_S, ierr)
     call generate_single_simulation(sim_handle, sqrt_S)
     deallocate(sqrt_S)
  end if

  ! Initialize chain on random position 
  samples(0,1) = 2.d0*pi*rand_uni(rng_handle) 
  samples(0,2) = acos(2.d0*rand_uni(rng_handle)-1.d0) 
  samples(0,3) = 0.d0
  samples(0,4) = 1.d0
  samples(0,5) = 0.d0
  samples(0,6) = 0.2d0 * rand_uni(rng_handle) 
  
  if (trim(operation) == 'MCMC') then
     ! Run sampler

     call comm_get_parameter(parfile, 'NSIDE_DIRECTION_MARGINAL',  par_int=nside_dir)
     call comm_get_parameter(parfile, 'MCMC_BURNIN',               par_int=burnin_mcmc)
     call comm_get_parameter(parfile, 'MARGINAL_OUTPUT_FREQUENCY', par_int=marg_freq)

     if (myid == root) write(*,*) 'Sampling model = ', trim(model)

     ! Open new chain file
     call comm_int2string(myid, i_s)
     open(unit,file=trim(prefix)//'_chain'//trim(adjustl(i_s))//'.dat',recl=1024)
     write(unit,fmt='(a)',advance='no') ' #     Sample     '
     do i = 1, npar
        write(unit,fmt='(a)',advance='no') par_label(i)
     end do
     write(unit,*) 'lnL'

     lnL_prev = lnL_lowl_full(samples(0,:))     
     if (verbosity > 0) call wall_time(t1)
     do i = 1, nsamp

        if (verbosity > 0) call wall_time(t2)
        if (mod(i,10) == 0) write(*,fmt='(a,i4,a,i8,a,f8.3)') 'myid = ', myid, ' -- generating sample no. ', i, &
             & ', wall time/sample = ', real(t2-t1)/(i-1)

        ! Propose changes to all parameters except direction
        do j = 1, npar
           eta(j) = rand_gauss(rng_handle)
        end do
        samples(i,:) = samples(i-1,:) + matmul(L_prop,eta)

        ! Propose changes to direction
        if (dir_prop_radius > 0.d0) then
           call propose_direction(rng_handle, dir_prop_radius, samples(i,1), samples(i,2))
        end if

        ! Check priors
        if (any(samples(i,:) < prior_uni(:,1)) .or. any(samples(i,:) > prior_uni(:,2))) then
           lnL = -1.d30
        else
           lnL          = lnL_lowl_full(samples(i,:))
        end if

        if (exp(lnL-lnL_prev) > rand_uni(rng_handle)) then
           ! Accept sample
           lnL_prev = lnL
        else
           ! Reject sample
           samples(i,:) = samples(i-1,:)
        end if
        write(unit,*) i, real(samples(i,:),sp), real(lnL_prev,sp)

        if (mod(i,prop_output_freq) == 0) call output_prop_matrix(outpropfile, 0.3d0, samples(i/2:i,:))
        if (mod(i,marg_freq) == 0)        call output_marginals(prefix, samples(1:i,:), burnin_mcmc)
     end do
     close(unit)

     ! Output final marginals
     call output_marginals(prefix, samples(1:i,:), burnin_mcmc)

  else if (trim(operation) == 'print_conditionals') then

     nsamp = 100
     allocate(x_n(0:nsamp), S_n(0:nsamp))
     do i = 4, npar ! Don't sample directional parameters
        if (par_rms(i) <= 0.d0) cycle
        dx           = (prior_uni(i,2)-prior_uni(i,1))/nsamp
        samples(1,:) = samples(0,:)
        do j = 0, nsamp
           x_n(j)       = prior_uni(i,1) + dx * j
           samples(1,i) = x_n(j)
           S_n(j)       = lnL_lowl_full(samples(1,:))
           write(*,*) i, j, real(x_n(j),sp), real(S_n(j),sp)
        end do
        S_n = exp(S_n - maxval(S_n))
        S_n = S_n / (sum(S_n)*dx)

        open(unit,file='cond_' // trim(par_label(i)) // '.dat')
        do j = 0, nsamp
           if (S_n(j) > 1.d-6 * maxval(S_n)) then
              write(unit,*) x_n(j), S_n(j)
           end if
        end do
        close(unit)
     end do
     deallocate(x_n, S_n)

  end if

  ! Output summary statistics


  ! Clean up
  call mpi_finalize(err)

contains

  function lnL_lowl_1D(x)
    use healpix_types
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: lnL_lowl_1D

    real(dp), allocatable, dimension(:) :: x_full

    allocate(x_full(npar))
    x_full            = x_current
    x_full(p_current) = x
    lnL_lowl_1D       = lnL_lowl_full(x_full)
    if (prior_gauss(p_current,2) > 0.d0) then
       lnL_lowl_1D  = lnL_lowl_1D - &
            & 0.5d0 * ((x-prior_gauss(p_current,1))/prior_gauss(p_current,2))**2
    end if
    deallocate(x_full)

  end function lnL_lowl_1D

  function lnL_lowl_full(x)
    use healpix_types
    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: lnL_lowl_full

    integer(i4b) :: l, ierr, nside
    real(dp)     :: chisq, t1, t2
    real(dp), allocatable, dimension(:,:) :: sqrt_S

    allocate(sqrt_S(n_h,n_h))
    if (verbosity > 0) call wall_time(t1)
    call get_sqrt_S(x, sqrt_S, ierr)
    if (verbosity > 0) call wall_time(t2)
!    if (verbosity > 0) write(*,*) 'CPU setup = ', real(t2-t1,sp), ' sec'

    !call generate_simulation_from_sqrtS(rng_handle, 32, 600.d0, sqrt_S, 'test.fits')

    if (verbosity > 0) call wall_time(t1)
    if (ierr == 0) then
       lnL_lowl_full = comm_lowl_compute_lnL(sqrt_S=sqrt_S(indmap,indmap), &
            & ierr=ierr, red_chisq=chisq)
    end if
    if (verbosity > 0) call wall_time(t2)
!    if (verbosity > 0) write(*,*) 'CPU lnL   = ', real(t2-t1,sp), ' sec'
    if (ierr /= 0) lnL_lowl_full = -1.d30
    deallocate(sqrt_S)

!    write(*,*) real(x(4:5),sp), real(lnL_lowl_full,sp), real(chisq,sp), ierr
    
  end function lnL_lowl_full

  subroutine propose_direction(rng_handle, radius, phi, theta)  
    implicit none

    type(planck_rng), intent(inout) :: rng_handle
    real(dp),         intent(in)    :: radius
    real(dp),         intent(inout) :: phi, theta

    real(dp) :: vec(3), euler(3,3), eta(3,3)

    if (radius == 0.d0) return

    call compute_euler_matrix_zyz(2.d0*pi*rand_uni(rng_handle), &
         & radius*rand_uni(rng_handle), 0.d0, eta)
    call compute_euler_matrix_zyz(phi, theta, 0.d0, euler)
    vec    = 0.d0
    vec(3) = 1.d0
    vec    = matmul(euler, matmul(eta, vec))
    call vec2ang(vec, theta, phi)
    
  end subroutine propose_direction

  subroutine generate_simulation_from_sqrtS(rng_handle, nside, fwhm, sqrt_S, filename)
    implicit none

    type(planck_rng),                 intent(inout) :: rng_handle
    integer(i4b),                     intent(in)    :: nside
    real(dp),                         intent(in)    :: fwhm
    real(dp),         dimension(:,:), intent(in)    :: sqrt_S
    character(len=*),                 intent(in)    :: filename

    integer(i4b) :: i, j, l, n, npix, nsim

    real(dp),     allocatable, dimension(:)     :: eta
    real(dp),     allocatable, dimension(:,:)   :: alms, map, variance
    complex(dpc), allocatable, dimension(:,:,:) :: alms_cmplx

     npix = 12*nside**2
     nsim = 100000

    allocate(eta(n_h), alms(numcomp,nmaps), alms_cmplx(nmaps,0:lmax,0:lmax))
    allocate(map(0:npix-1,nmaps), variance(0:npix-1,nmaps))

    variance = 0.d0
    do j = 1, nsim
       write(*,*) j, nsim

       ! Draw vector of standard Gaussian random variates, and multiply with sqrt_S
       do i = 1, n_h
          eta(i) = rand_gauss(rng_handle)
       end do
       eta = matmul(sqrt_S, eta)
       
       ! Perform spherical harmonics transform
       do i = 1, nmaps
          alms(:,i) = eta((i-1)*numcomp+1:i*numcomp)
       end do
       call convert_real_to_complex_alms_dp(alms, alms_cmplx)
       do l = 0, lmax
          alms_cmplx(:,l,:) = alms_cmplx(:,l,:) * &
               & exp(-0.5*l*(l+1.d0)*(fwhm*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
       end do
       if (nmaps == 1) then
          call alm2map(nside, lmax, lmax, alms_cmplx, map(:,1))
       else
          call alm2map(nside, lmax, lmax, alms_cmplx, map)
       end if
       call write_map3(filename, map)

       variance = variance + map**2

       if (mod(j,1000) == 0) then
          call write_map3('stddev.fits', sqrt(variance/j))
       end if

    end do

    variance = sqrt(variance/nsim)
    call write_map3('stddev.fits', variance)

    deallocate(map, eta, alms, alms_cmplx, variance)
    stop

  end subroutine generate_simulation_from_sqrtS


  subroutine generate_single_simulation(rng_handle, sqrt_S)
    implicit none

    type(planck_rng),                 intent(inout) :: rng_handle
    real(dp),         dimension(:,:), intent(in)    :: sqrt_S

    integer(i4b) :: i, j, l, m, n, ind, n_h, nmaps, lmax

    real(dp),     allocatable, dimension(:)     :: eta
    real(dp),     allocatable, dimension(:,:)   :: sqrt_N, map

    n     = comm_lowl(1)%n
    n_h   = comm_lowl(1)%n_h
    lmax  = comm_lowl(1)%lmax
    nmaps = comm_lowl(1)%nmaps

    ! Draw vector of standard Gaussian random variates, and multiply with sqrt_S
    allocate(eta(n_h))
    ind = 1
    eta = 0.d0
    do i = 1, nmaps
       do l = 0, lmax
          if (l < lmin_fit .or. l > lmax_fit) then
             ind = ind + 2*l+1
             cycle
          end if
          do m = -l, l
             eta(ind) = rand_gauss(rng_handle)
             ind      = ind+1
          end do
       end do
    end do
    eta = matmul(sqrt_S, eta)

    ! Multiply with beam
    ind = 1
    do j = 1, nmaps
       do l = 0, lmax
          do m = -l, l
             eta(ind) = eta(ind) * comm_lowl(1)%beam(l,j)
             ind = ind+1
          end do
       end do
    end do

    ! Project onto correct basis
    comm_lowl(1)%d(:,1) = matmul(transpose(comm_lowl(1)%P_harm), eta)
    deallocate(eta)

!!$    allocate(map(0:12*16**2-1,nmaps))
!!$    ind = 1
!!$    do j = 1, nmaps
!!$       do i = 0, 12*16**2-1
!!$          map(i,j) = comm_lowl(1)%d(ind,1)
!!$          ind = ind+1
!!$       end do
!!$    end do
!!$    call write_map3('sig.fits', map)
!!$    deallocate(map)

    ! Add noise
    !comm_lowl(1)%N_cov = comm_lowl(1)%N_cov
    allocate(sqrt_N(comm_lowl(1)%n,comm_lowl(1)%n), eta(comm_lowl(1)%n))
    sqrt_N = comm_lowl(1)%N_cov

    !call cholesky_decompose_with_mask_dp(sqrt_N)
    call compute_hermitian_root(sqrt_N, 0.5d0)
    do i = 1, n
       eta(i) = rand_gauss(rng_handle)
    end do
    comm_lowl(1)%d(:,1) = comm_lowl(1)%d(:,1) + matmul(sqrt_N, eta)
    !comm_lowl(1)%d(:,1) = comm_lowl(1)%d(:,1) 
    !comm_lowl(1)%d(:,1) = matmul(sqrt_N, eta)

!!$    allocate(map(0:12*16**2-1,nmaps))
!!$    ind = 1
!!$    do j = 1, nmaps
!!$       do i = 0, 12*16**2-1
!!$          map(i,j) = comm_lowl(1)%d(ind,1)
!!$          ind = ind+1
!!$       end do
!!$    end do
!!$    call write_map3('sim.fits', map)
!!$    deallocate(map)
!!$    stop
    
    deallocate(eta, sqrt_N)

  end subroutine generate_single_simulation

  subroutine output_marginals(prefix, samples_in, burnin)
    implicit none

    character(len=*),                   intent(in) :: prefix
    real(dp),         dimension(1:,1:), intent(in) :: samples_in
    integer(i4b),                       intent(in) :: burnin

    integer(i4b)       :: c, i, j, p, nsamp, npar, unit, npix, nchain, myid, root
    character(len=512) :: filename
    real(dp), allocatable, dimension(:,:)       :: samples
    real(dp), allocatable, dimension(:,:), save :: map
    integer(i4b), save :: currsamp
    
    root = 0
    call mpi_comm_rank(mpi_comm_world, myid, err)
    call mpi_comm_size(mpi_comm_world, nchain, err)

    nsamp = size(samples_in,1)
    npar  = size(samples_in,2)
    unit  = comm_getlun()
    if (nsamp <= burnin) return

    if (.not. allocated(map)) then 
       npix = 12*nside_dir**2
       allocate(map(0:npix-1,1))
       map = 0.d0
       currsamp = burnin
    end if

    allocate(samples(currsamp+1:nsamp,npar))
    samples = samples_in(currsamp+1:nsamp,:)
    if (myid == root) then

       write(*,fmt='(a,i7,a,i7,a)') 'Writing samples ', currsamp, ' to ', nsamp, ' to marginals'

       do c = 0, nchain-1

          if (c > 0) then
             call mpi_recv(samples, size(samples), mpi_double_precision, c, 0, mpi_comm_world, status, err)
          end if

          ! First convert directions into sky maps
          if (dir_prop_radius > 0.d0) then
             do j = currsamp+1, nsamp
                call ang2pix_ring(nside_dir, samples(j,2), samples(j,1), p)
                map(p,1) = map(p,1) + 1.d0
             end do
             filename = trim(prefix) // '_marg_dir.fits'
             call write_map3(filename, map)
          end if

          ! Then print other parameters to ASCII table
          do i = 3, npar
             if (par_rms(i) <= 0.d0) cycle
             filename = trim(prefix) // '_marg_' // trim(par_label(i)) // '.dat'
             if (currsamp == burnin .and. c == 0) then
                open(unit,file=trim(filename), status='REPLACE')
             else
                open(unit,file=trim(filename), status='OLD', position='APPEND')
             end if
             do j = currsamp+1, nsamp
                write(unit,*) real(samples(j,i),sp)
             end do
             close(unit)
          end do

       end do

    else
       call mpi_send(samples, size(samples), mpi_double_precision, 0, 0, mpi_comm_world, err)
    end if

    deallocate(samples)

    currsamp = nsamp

  end subroutine output_marginals

  
!!$  subroutine generate_dipole(A, theta, phi)
!!$    implicit none
!!$    real(dp), intent(in)                  :: A, theta, phi
!!$    real(dp), dimension(:,:), allocatable :: map
!!$    real(dp), dimension(:),   allocatable :: dip_vector, pix_vector
!!$    integer(i4b)                          :: i
!!$
!!$    print*, "generate dipole map"
!!$    
!!$    allocate(map(0:12*16**2-1,1), dip_vector(3), pix_vector(3))
!!$    
!!$    dip_vector(1) = sin(theta) * cos(phi)
!!$    dip_vector(2) = sin(theta) * sin(phi)
!!$    dip_vector(3) = cos(theta)
!!$
!!$    map = 0.d0
!!$    do i = 0, 12*16**2-1
!!$       call pix2vec_ring(16, i, pix_vector(:))
!!$       map(i,1) = A * sum(pix_vector*dip_vector)
!!$    end do
!!$    call write_map3('dip.fits', map)
!!$
!!$  end subroutine generate_dipole

end program comm_like_sampler
