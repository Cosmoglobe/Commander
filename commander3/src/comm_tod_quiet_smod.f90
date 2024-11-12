submodule(comm_tod_QUIET_mod) comm_tod_QUIET_smod 
  implicit none
  !
  !   Submodule which contains implementation of the procedures defined in its parent
  !   module (comm_tod_QUIET_mod)
  !
  !   Main Methods
  !   ------------
  !   constructor(cpar, id_abs, info, tod_type)
  !       Initialization routine that reads in, allocates and associates 
  !       all data needed for TOD processing
  !   process_QUIET_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
  !       Routine which processes the time ordered data
  !

contains

  !**************************************************
  !             Constructor
  !**************************************************
  !function constructor(cpar, id_abs, info, tod_type)
  module procedure constructor
    implicit none
    integer(i4b) :: i, nside_beam, lmax_beam, nmaps_beam
    logical(lgt) :: pol_beam

    ! Initialize common parameters
    allocate (constructor)

    ! Set up noise PSD type and priors
    constructor%freq            = cpar%ds_label(id_abs)
    constructor%n_xi            = 3
    ! oof -- one over f noise
    constructor%noise_psd_model = 'oof'
    allocate(constructor%xi_n_P_uni(constructor%n_xi,2))
    allocate(constructor%xi_n_P_rms(constructor%n_xi))

    !
    constructor%xi_n_P_rms      = [-1.0, 0.1, 0.2]   ! [sigma0, fknee, alpha]; sigma0 is not used
    if (trim(constructor%freq) == 'Q') then
       constructor%xi_n_nu_fit = [0.0, 0.200]    
       constructor%xi_n_P_uni(2,:)  = [0.00001, 0.005]  ! fknee
       constructor%xi_n_P_uni(3,:)  = [-3.0, -0.01]     ! alpha
    else if (trim(constructor%freq) == 'W') then
       constructor%xi_n_nu_fit = [0.0, 0.200]    
       constructor%xi_n_P_uni(2,:)  = [0.0001, 0.01]    ! fknee
       constructor%xi_n_P_uni(3,:)  = [-3.0, -0.01]     ! alpha
    else
       write(*,*) 'Invalid QUIET frequency label = ', trim(constructor%freq)
       stop
    end if

    !
    call constructor%tod_constructor(cpar, id_abs, info, tod_type)

    ! Set up WMAP specific parameters
    constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
    constructor%nhorn           = 1 ! number of horns 
    !constructor%ndiode          = 4 ! number of diodes in an amplifier (from adc branch)
    constructor%n_xi            = 3 ! number of power spectrum density parameters (sigma0, fknee, alpha ?)
    constructor%compressed_tod  = .false.
    constructor%correct_sl      = .false.    ! sidelobes to be removed
    constructor%orb_4pi_beam    = .false.    ! full beam convolution
    constructor%symm_flags      = .false.
    constructor%chisq_threshold = 400.d0     ! throwing away the tod above thhat value
    constructor%nmaps           = info%nmaps ! map labels (?) 
    constructor%ndet            = num_tokens(cpar%ds_tod_dets(id_abs), ",") ! number of detectors (retrieved from HDF5 file)
    constructor%verbosity       = cpar%verbosity ! how extensive your output should be
    !constructor%x_im            = 0d0 ! WMAP only parameter

    ! Iniitialize TOD labels
    allocate (constructor%labels(7))
    constructor%labels(1) = 'map'    ! actual solution to the map 
    constructor%labels(2) = 'res'    ! residual to the map
    constructor%labels(3) = 'ncorr'  ! correlated noise
    constructor%labels(4) = 'bpcorr' ! bandpass correction
    constructor%labels(5) = 'orb'    ! orbital dipole
    constructor%labels(6) = 'sl'     ! sidelobe
    constructor%labels(7) = 'zodi'   ! zodiacal light

    ! Initialize beams
    nside_beam                  = 512
    nmaps_beam                  = 3
    pol_beam                    = .true.
    constructor%nside_beam      = nside_beam


    ! Get detector labels
    call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label)

    ! Define detector partners
    do i = 1, constructor%ndet
       if (mod(i,2) == 1) then
          constructor%partner(i) = i+1
       else
          constructor%partner(i) = i-1
       end if
       constructor%horn_id(i) = (i+1)/2
    end do

    ! Read the actual TOD
    call constructor%read_tod(constructor%label)

    ! Initialize bandpass mean and proposal matrix
    call constructor%initialize_bp_covar(trim(cpar%datadir)//cpar%ds_tod_bp_init(id_abs))

    ! Construct lookup tables -- reads everything up from pointing information
    call constructor%precompute_lookups()

    ! Load the instrument file
    call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

    ! Need precompute the main beam precomputation for both the A-horn and
    ! B-horn.
    ! Allocate sidelobe convolution data structures
    allocate(constructor%slconv(constructor%ndet), constructor%orb_dp)
    constructor%orb_dp => comm_orbdipole(constructor%mbeam)

    ! Initialize all baseline corrections to zero
    !do i = 1, constructor%nscan
    !   constructor%scans(i)%d%baseline = 0.d0
    !end do

  end procedure constructor

  !**************************************************
  !             Driver routine
  !**************************************************
  module procedure process_QUIET_tod
    implicit none
    real(dp)     :: t1, t2
    integer(i4b) :: i, j, k, l, n
    integer(i4b) :: nside, npix, nmaps 
    integer(i4b) :: ierr, ndelta
    real(sp), allocatable, dimension(:, :)    :: s_buf
    real(sp), allocatable, dimension(:, :, :) :: d_calib
    real(dp), allocatable, dimension(:, :)    :: chisq_S, m_buf
    real(dp), allocatable, dimension(:, :)    :: M_diag
    real(dp), allocatable, dimension(:, :, :) :: b_map, b_mono, sys_mono
    character(len=512)                        :: prefix, postfix
    character(len=2048)                       :: Sfilename

    logical(lgt)        :: select_data, sample_abs_bandpass
    logical(lgt)        :: sample_rel_bandpass, output_scanlist
    type(comm_binmap)   :: binmap
    type(comm_scandata) :: sd

    character(len=4)   :: ctext, myid_text
    character(len=6)   :: samptext, scantext
    character(len=512), allocatable, dimension(:) :: slist
    real(sp),       allocatable, dimension(:)     :: procmask, procmask2, sigma0
    real(sp),  allocatable, dimension(:, :, :, :) :: map_sky
    class(map_ptr),     allocatable, dimension(:) :: outmaps

    ! biconjugate gradient-stab parameters
    integer(i4b) :: num_cg_iters
    real(dp) ::  epsil(6)
    real(dp), allocatable, dimension(:, :, :) :: bicg_sol
    real(dp), allocatable, dimension(:)       :: map_full

    ! Toggle optional operations
    sample_rel_bandpass   = size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
    sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses
    select_data           = .false.                ! only perform data selection the first time
    ! The first couple of iterations have quite bad chi-square for certain
    ! detectors, but it quickly settles down.
    output_scanlist       = mod(iter-1,10) == 0    ! only output scanlist every 10th iteration

    ! Initialize local variables
    ndelta          = size(delta,3)
    self%n_bp_prop  = ndelta-1
    nside           = map_out%info%nside
    nmaps           = map_out%info%nmaps
    npix            = 12*nside**2
    self%output_n_maps = 1
    if (self%output_aux_maps > 0) then
       if (mod(iter-1,self%output_aux_maps) == 0) self%output_n_maps = 7
    end if

    call int2string(chain, ctext)
    call int2string(iter, samptext)
    call int2string(self%myid, myid_text)
    prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
    postfix = '_c' // ctext // '_k' // samptext // '.fits'

    ! Distribute maps
    allocate(map_sky(nmaps,self%nobs,0:self%ndet,ndelta))
    allocate(map_full(0:npix-1))
    map_full = 0.d0
    !call distribute_sky_maps(self, map_in, 1.e-3, map_sky) ! uK to mK
    call distribute_sky_maps(self, map_in, 1., map_sky, map_full) ! K to K?

    ! Distribute processing masks
    allocate(m_buf(0:npix-1,nmaps), procmask(0:npix-1), procmask2(0:npix-1))
    call self%procmask%bcast_fullsky_map(m_buf);  procmask  = m_buf(:,1)
    call self%procmask2%bcast_fullsky_map(m_buf); procmask2 = m_buf(:,1)
    deallocate(m_buf)

    ! Allocate total map (for monopole sampling)

    ! Precompute far sidelobe Conviqt structures
    if (self%correct_sl) then
       do i = 1, self%ndet
          !TODO: figure out why this is rotated
          call map_in(i,1)%p%YtW()  ! Compute sky a_lms
          self%slconv(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
               & 100, 3, 100, self%slbeam(i)%p, map_in(i,1)%p, 2)
       end do
    end if

    call update_status(status, "tod_init")

    !------------------------------------
    ! Perform main sampling steps
    !------------------------------------
    ! For QUIET we need to only sample for absolute calibration
    !call sample_baseline(self, handle, map_sky, procmask, procmask2)
    call sample_calibration(self, 'abscal', handle, map_sky, procmask, procmask2)
    !call sample_calibration(self, 'relcal', handle, map_sky, procmask, procmask2)
    !call sample_calibration(self, 'deltaG', handle, map_sky, procmask, procmask2)
    !call sample_calibration(self, 'imbal',  handle, map_sky, procmask, procmask2)
    ! Write out the way that WMAP calculated the imbalance parameters.

    ! Prepare intermediate data structures
    if (sample_abs_bandpass .or. sample_rel_bandpass) then
       allocate(chisq_S(self%ndet,size(delta,3)))
       chisq_S = 0.d0
    end if
    if (output_scanlist) then
       allocate(slist(self%nscan))
       slist   = ''
    end if

    allocate (M_diag(0:npix-1, nmaps+1))
    allocate ( b_map(0:npix-1, nmaps, self%output_n_maps))
    M_diag = 0d0
    b_map = 0d0

    ! Perform loop over scans
    if (self%myid == 0) write(*,*) '   --> Sampling ncorr, xi_n, maps'
    do i = 1, self%nscan

       ! Skip scan if no accepted data
       if (.not. any(self%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       if (sample_rel_bandpass) then
          !call sd%init_differential(self, i, map_sky, procmask, procmask2, &
          !  & init_s_bp=.true., init_s_bp_prop=.true.)
         call init_scan_data_singlehorn(sd, self, i, map_sky, procmask, procmask2, init_s_bp=.true., init_s_bp_prop=.true.)
       else if (sample_abs_bandpass) then
         ! call sd%init_differential(self, i, map_sky, procmask, procmask2, &
         !   & init_s_bp=.true., init_s_sky_prop=.true.)
         call init_scan_data_singlehorn(sd, self, i, map_sky, procmask, procmask2, init_s_bp=.true., init_s_bp_prop=.true.)
       else
         ! call sd%init_differential(self, i, map_sky, procmask, procmask2, &
         !   & init_s_bp=.true.)
         call init_scan_data_singlehorn(sd, self, i, map_sky, procmask, procmask2, init_s_bp=.true.)
       end if
       allocate(s_buf(sd%ntod,sd%ndet))

       ! Calling Simulation Routine
       ! Not implemented for differential
       !if (self%enable_tod_simulations) then
       !   call simulate_tod(self, i, sd%s_tot, handle)
       !   call sd%dealloc
       !   cycle
       !end if

       ! Sample correlated noise
       call sample_n_corr(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr, &
         & sd%pix(:,1,:), dospike=.false.)

       ! Compute noise spectrum parameters
       call sample_noise_psd(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr)

       ! Compute chisquare
       do j = 1, sd%ndet
          if (.not. self%scans(i)%d(j)%accept) cycle
          call self%compute_chisq(i, j, sd%mask(:,j), sd%s_sky(:,j), &
            & sd%s_sl(:,j) + sd%s_orb(:,j), sd%n_corr(:,j), sd%tod(:,j))
       end do
       ! Should we be computing the chisq of the sums and differences?
       ! I also think I just need to double check the chisq calculation,
       ! because to me, the fit looks quite good.

       ! Select data
       if (select_data) call remove_bad_data(self, i, sd%flag)

       ! Compute chisquare for bandpass fit
       if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)

       ! Compute binned map
       allocate(d_calib(self%output_n_maps,sd%ntod, sd%ndet))
       call compute_calibrated_data(self, i, sd, d_calib)
       if (.false. .and. i==1 .and. mod(iter,10) == 0) then
          call int2string(self%scanid(i), scantext)
          if (self%myid == 0 .and. self%verbosity > 0) write(*,*) 'Writing tod to txt'
          do k = 1, self%ndet
             open(78,file=trim(chaindir)//'/tod_'//trim(self%label(k))//'_pid'//scantext//'_samp'//samptext//'.dat', recl=1024)
             write(78,*) "# Sample   uncal_TOD (mK)  n_corr (mK) cal_TOD (mK)  sky (mK)  "// &
                  & " s_orb (mK),  mask, baseline, sl, bp, gain, sigma0"
             do j = 1, sd%ntod
                write(78,*) j, sd%tod(j, k), sd%n_corr(j, k), d_calib(1,j,k), &
                 &  sd%s_totA(j,k), sd%s_orbA(j,k), &
                 &  sd%s_totB(j,k), sd%s_orbB(j,k), &
                 &  sd%mask(j, k), self%scans(i)%d(k)%baseline, &
                 &  sd%s_sl(j,k),  sd%s_bp(j,k), real(self%scans(i)%d(k)%gain, sp), &
                 &  real(self%scans(i)%d(k)%N_psd%sigma0, sp)
             end do
             close(78)
          end do
       end if
       
       ! Output 4D map; note that psi is zero-base in 4D maps, and one-base in Commander
       ! if (self%output_4D_map > 0) then
       !    if (mod(iter-1,self%output_4D_map) == 0) then
       !       allocate(sigma0(sd%ndet))
       !       do j = 1, sd%ndet
       !          sigma0(j) = self%scans(i)%d(j)%N_psd%sigma0/self%scans(i)%d(j)%gain
       !       end do
       !       call output_4D_maps_hdf(trim(chaindir) // '/tod_4D_chain'//ctext//'_proc' // myid_text // '.h5', &
       !            & samptext, self%scanid(i), self%nside, self%npsi, &
       !            & self%label, self%horn_id, real(self%polang*180/pi,sp), sigma0, &
       !            & sd%pix(:,:,1), sd%psi(:,:,1)-1, d_calib(1,:,:), iand(sd%flag,self%flag0), &
       !            & self%scans(i)%d(:)%accept)
       !       deallocate(sigma0)
       !    end if
       ! end if

       ! Bin TOD
       !call bin_differential_TOD(self, d_calib, sd%pix(:,1,:),  &
       !  & sd%psi(:,1,:), sd%flag(:,1), self%x_im, procmask, b_map, M_diag, i, &
       !  & .false.)
       call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib, binmap)

       ! Update scan list
       call wall_time(t2)
       self%scans(i)%proctime   = self%scans(i)%proctime   + t2-t1
       self%scans(i)%n_proctime = self%scans(i)%n_proctime + 1
       if (output_scanlist) then
          write(slist(i),*) self%scanid(i), '"',trim(self%hdfname(i)), &
               & '"', real(self%scans(i)%proctime/self%scans(i)%n_proctime,sp),&
               & real(self%spinaxis(i,:),sp)
       end if

       ! Clean up
       call dealloc_scandata(sd)
       deallocate(s_buf, d_calib)

    end do

    !------------------------------------
    ! TODO: From line 507 to the 568 is what you will need to rewrite for QUIET
    ! Look into `run_bicgstab` inside comm_tod_mapmaking_mod.f90 <= written for WMAP
    ! Look into `finalize_binned_map` inside comm_tod_mapmaking_mod.f90 <= written for LFI
    if (self%myid == 0) write(*,*) '   --> Finalizing maps, bp'

    ! Output latest scan list with new timing information
    if (output_scanlist) call self%output_scan_list(slist)


    call update_status(status, "Running allreduce on M_diag")
    call mpi_allreduce(mpi_in_place, M_diag, size(M_diag), &
         & MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)
    call update_status(status, "Running allreduce on b")
    call mpi_allreduce(mpi_in_place, b_map, size(b_map), &
         & MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)

    where (M_diag == 0d0)
       M_diag = 1d0
    end where
    ! If we want to not do the "better preconditioning"
    M_diag(:,4) = 0d0

    ! Conjugate Gradient solution to (P^T Ninv P) m = P^T Ninv d, or Ax = b
    call update_status(status, "Allocating cg arrays")
    ! bicg_sol is the actual map
    allocate (bicg_sol(0:npix-1, nmaps, self%output_n_maps))
    allocate(m_buf(0:npix-1,nmaps))
    call map_in(1,1)%p%bcast_fullsky_map(m_buf)
    bicg_sol = 0.0d0
    !bicg_sol(:,:,1) = m_buf
    deallocate(m_buf)

    epsil(1)   = 1d-10
    !epsil(1)   = 1d-8
    epsil(2:6) = 1d-6
    num_cg_iters = 0
    if (self%myid == 0) then 
       if (self%verbosity > 0) write(*,*) '  Running BiCG'
    end if

    ! Doing this now because it's still burning in...
    !if (mod(iter-1,10*self%output_aux_maps) == 0) then
      ! Solve for maps
      call update_status(status, "Starting bicg-stab")
      do l=1, self%output_n_maps
         if (self%verbosity > 0 .and. self%myid == 0) then
           write(*,*) '    Solving for ', trim(adjustl(self%labels(l)))
         end if
         call run_bicgstab(self, handle, bicg_sol, npix, nmaps, num_cg_iters, &
                        & epsil(l), procmask, map_full, M_diag, b_map, l, &
                        & prefix, postfix)
      end do
      if (self%verbosity > 0 .and. self%myid == 0) write(*,*) '  Finished BiCG'
    !end if

    call MPI_Bcast(bicg_sol, size(bicg_sol),  MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)
    call MPI_Bcast(num_cg_iters, 1,  MPI_INTEGER, 0, self%info%comm, ierr)
    allocate(outmaps(self%output_n_maps))
    do i = 1, self%output_n_maps
       outmaps(i)%p => comm_map(self%info)
    end do
    do k = 1, self%output_n_maps
       do j = 1, nmaps
          outmaps(k)%p%map(:, j) = bicg_sol(self%info%pix, j, k)
       end do
    end do
    !------------------------------------

    ! Output maps to the files
    map_out%map = outmaps(1)%p%map
    rms_out%map = M_diag(self%info%pix, 1:nmaps)**-0.5
    call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
    call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
    do n = 2, self%output_n_maps
      call outmaps(n)%p%writeFITS(trim(prefix)//trim(adjustl(self%labels(n)))//trim(postfix))
    end do

    ! Sample bandpass parameters
    if (sample_rel_bandpass .or. sample_abs_bandpass) then
       call sample_bp(self, iter, delta, map_sky, handle, chisq_S)
       self%bp_delta = delta(:,:,1)
    end if

    ! Clean up temporary arrays
    deallocate(procmask, procmask2)
    deallocate(b_map, M_diag)
    deallocate(map_full)
    if (allocated(chisq_S)) deallocate (chisq_S)
    if (allocated(b_mono)) deallocate (b_mono)
    if (allocated(sys_mono)) deallocate (sys_mono)
    if (allocated(slist)) deallocate (slist)

    if (allocated(outmaps)) then
       do i = 1, self%output_n_maps
          call outmaps(i)%p%dealloc
       end do
       deallocate (outmaps)
    end if

    deallocate (map_sky, bicg_sol)

    if (self%correct_sl) then
       do i = 1, self%ndet
          call self%slconv(i)%p%dealloc(); deallocate(self%slconv(i)%p)
       end do
    end if

    call int2string(iter, ctext)
    call update_status(status, "tod_end"//ctext)

    ! Parameter to check if this is first time routine has been
    self%first_call = .false.

  end procedure process_QUIET_tod

end submodule comm_tod_QUIET_smod
