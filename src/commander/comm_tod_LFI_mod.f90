module comm_tod_LFI_mod
  use comm_tod_mod
  use comm_param_mod
  use comm_map_mod
  use comm_conviqt_mod
  use pix_tools
  use healpix_types  
  use comm_huffman_mod
  use comm_hdf_mod
  use comm_shared_arr_mod
  implicit none

  private
  public comm_LFI_tod

  integer(i4b), parameter :: N_test      = 15
  integer(i4b), parameter :: samp_N      = 1
  integer(i4b), parameter :: prep_G      = 15
  integer(i4b), parameter :: samp_G      = 2
  integer(i4b), parameter :: prep_acal   = 3
  integer(i4b), parameter :: samp_acal   = 4
  integer(i4b), parameter :: prep_bp     = 5
  integer(i4b), parameter :: samp_bp     = 11
  integer(i4b), parameter :: samp_sl     = 6
  integer(i4b), parameter :: samp_N_par  = 7
  integer(i4b), parameter :: sel_data    = 8
  integer(i4b), parameter :: bin_map     = 9
  integer(i4b), parameter :: calc_chisq  = 10
  integer(i4b), parameter :: output_slist = 12 
  integer(i4b), parameter :: samp_mono   = 13
  integer(i4b), parameter :: sub_sl      = 14
  logical(lgt), dimension(N_test) :: do_oper


  type, extends(comm_tod) :: comm_LFI_tod 
   contains
     procedure     :: process_tod        => process_LFI_tod
     procedure     :: compute_binned_map
     procedure     :: project_sky
     procedure     :: compute_orbital_dipole
     procedure     :: sample_gain
     procedure     :: sample_smooth_gain
     procedure     :: sample_n_corr
     procedure     :: sample_bp
     procedure     :: compute_cleaned_tod
     procedure     :: construct_sl_template
     procedure     :: compute_chisq
     procedure     :: sample_noise_psd
     procedure     :: decompress_pointing_and_flags
     procedure     :: accumulate_absgain_from_orbital
     procedure     :: sample_absgain_from_orbital
     procedure     :: sample_mono
     procedure     :: get_total_chisq
     procedure     :: output_scan_list
     procedure     :: finalize_binned_map
  end type comm_LFI_tod

  interface comm_LFI_tod
     procedure constructor
  end interface comm_LFI_tod

contains

  !**************************************************
  !             Constructor
  !**************************************************
  function constructor(cpar, id_abs, info)
    implicit none
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id_abs
    class(comm_mapinfo),     target     :: info
    class(comm_LFI_tod),     pointer    :: constructor

    integer(i4b) :: i, j, l, m, nside_beam, lmax_beam, nmaps_beam, ndelta, np_vec
    character(len=512) :: datadir
    logical(lgt) :: pol_beam
    type(hdf_file) :: h5_file

    ! Set up LFI specific parameters
    allocate(constructor)
    constructor%myid          = cpar%myid_chain
    constructor%comm          = cpar%comm_chain
    constructor%numprocs      = cpar%numprocs
    constructor%myid_shared   = cpar%myid_shared
    constructor%comm_shared   = cpar%comm_shared
    constructor%myid_inter    = cpar%myid_inter
    constructor%comm_inter    = cpar%comm_inter
    constructor%info          => info
    constructor%output_all    = .true.
    constructor%init_from_HDF = cpar%ds_tod_initHDF(id_abs)
    constructor%freq          = cpar%ds_label(id_abs)
    constructor%operation     = cpar%operation

    datadir = trim(cpar%datadir)//'/' 
    constructor%filelist    = trim(datadir)//trim(cpar%ds_tod_filelist(id_abs))
    constructor%procmaskf1  = trim(datadir)//trim(cpar%ds_tod_procmask1(id_abs))
    constructor%procmaskf2  = trim(datadir)//trim(cpar%ds_tod_procmask2(id_abs))
    constructor%instfile    = trim(datadir)//trim(cpar%ds_tod_instfile(id_abs))

    call constructor%get_scan_ids(constructor%filelist)

    constructor%procmask => comm_map(constructor%info, constructor%procmaskf1)
    constructor%procmask2 => comm_map(constructor%info, constructor%procmaskf2)
    do i = 0, constructor%info%np-1
       if (any(constructor%procmask%map(i,:) < 0.5d0)) then
          constructor%procmask%map(i,:) = 0.d0
       else
          constructor%procmask%map(i,:) = 1.d0
       end if
       if (any(constructor%procmask2%map(i,:) < 0.5d0)) then
          constructor%procmask2%map(i,:) = 0.d0
       else
          constructor%procmask2%map(i,:) = 1.d0
       end if
    end do


    constructor%nmaps    = info%nmaps
    constructor%ndet     = num_tokens(cpar%ds_tod_dets(id_abs), ",")
    constructor%nhorn    = 1
    allocate(constructor%stokes(constructor%nmaps))
    allocate(constructor%w(constructor%nmaps,constructor%nhorn,constructor%ndet))
    allocate(constructor%label(constructor%ndet))
    constructor%stokes = [1,2,3]
    constructor%w      = 1.d0
    call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label)

    if (trim(cpar%ds_bpmodel(id_abs)) == 'additive_shift') then
       ndelta = 1
    else if (trim(cpar%ds_bpmodel(id_abs)) == 'powlaw_tilt') then
       ndelta = 1
    else
       write(*,*) 'Unknown bandpass model'
       stop
    end if
    allocate(constructor%bp_delta(constructor%ndet,ndelta))

    ! Read the actual TOD
    call constructor%read_tod(constructor%label)

    if (.true.) then
       ! Initialize far sidelobe beams
       call open_hdf_file(constructor%instfile, h5_file, 'r')
       nside_beam = 128
       nmaps_beam = 3
       pol_beam   = .true.
       call read_hdf(h5_file, trim(adjustl(constructor%label(1)))//'/'//'sllmax', lmax_beam)
       constructor%slinfo => comm_mapinfo(cpar%comm_chain, nside_beam, lmax_beam, &
            & nmaps_beam, pol_beam)
       allocate(constructor%slbeam(constructor%ndet), constructor%slconv(constructor%ndet))
       do i = 1, constructor%ndet
          constructor%slbeam(i)%p => comm_map(constructor%slinfo, h5_file, .true., .true., trim(constructor%label(i)))

!!$          call constructor%slbeam(i)%p%Y
!!$          call constructor%slbeam(i)%p%writeFITS("beam.fits")
!!$          call mpi_finalize(m)
!!$          stop
          

!!$          do j = 0, constructor%slbeam(i)%p%info%nalm-1
!!$             l = constructor%slbeam(i)%p%info%lm(1,j)
!!$             m = constructor%slbeam(i)%p%info%lm(2,j)
!!$             if (m == 0) then
!!$                constructor%slbeam(i)%p%alm(j,1) = exp(-0.5d0*l*(l+1)*(420.d0*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
!!$             else
!!$                constructor%slbeam(i)%p%alm(j,1) = 0.d0
!!$             end if
!!$          end do

       end do
       call close_hdf_file(h5_file)
    end if

    ! Lastly, create a vector pointing table for fast look-up for orbital dipole
    np_vec = 12*constructor%nside**2 !npix
    allocate(constructor%pix2vec(3,0:np_vec-1))
    do i = 0,np_vec -1
       call pix2vec_ring(constructor%nside, i, constructor%pix2vec(:,i))
    end do
           
  end function constructor

  !**************************************************
  !             Driver routine
  !**************************************************
  subroutine process_LFI_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
    implicit none
    class(comm_LFI_tod),                   intent(inout) :: self
    character(len=*),                      intent(in)    :: chaindir
    integer(i4b),                          intent(in)    :: chain, iter
    type(planck_rng),                      intent(inout) :: handle
    type(map_ptr),       dimension(:,:),   intent(inout) :: map_in            ! (ndet,ndelta)
    real(dp),            dimension(:,:,:), intent(inout) :: delta             ! (ndet,npar,ndelta) BP corrections
    class(comm_map),                       intent(inout) :: map_out, rms_out  ! Combined output map and rms

    integer(i4b) :: i, j, k, l, ntod, ndet, nside, npix, nmaps, naccept, ntot, ns
    integer(i4b) :: ierr, main_iter, n_main_iter, ndelta, scanfile, ncol, n_A, nout
    real(dp)     :: t1, t2, t3, t4, t5, t6, t7, t8, chisq_threshold, delta_temp, chisq_tot
    real(dp)     :: t_tot(22), inv_gain
    real(sp),     allocatable, dimension(:,:)     :: n_corr, s_sl, s_sky, s_orb, mask,mask2, s_sb, s_sky_prop, s_sb_prop, s_mono, s_buf, s_tot
    real(sp),     allocatable, dimension(:,:,:)   :: d_calib
    real(dp),     allocatable, dimension(:,:,:,:) :: map_sky
    real(dp),     allocatable, dimension(:)       :: A_abscal, b_abscal
    real(dp),     allocatable, dimension(:,:)     :: chisq_S
    real(dp),     allocatable, dimension(:,:)     :: A_map
    integer(i4b), allocatable, dimension(:,:)     :: nhit
    real(dp),     allocatable, dimension(:,:,:)   :: b_map, b_mono, sys_mono
    integer(i4b), allocatable, dimension(:,:)     :: pix, psi, flag
    logical(lgt), save :: first = .true.
    logical(lgt) :: correct_sl
    character(len=512) :: prefix, postfix, Sfilename
    character(len=4)   :: ctext
    character(len=6)   :: samptext
    character(len=512), allocatable, dimension(:) :: slist
    character(len=8) :: id, stext
    character(len=1) :: did
    character(5)                                  :: istr
    type(shared_1d_int) :: sprocmask, sprocmask2
    type(shared_2d_sp), allocatable, dimension(:,:) :: smap_sky
    class(comm_map), pointer :: map_sl, map_sl2, condmap
    class(map_ptr), allocatable, dimension(:) :: outmaps
    class(comm_mapinfo), pointer :: info_sl2, info_sl_map

    call update_status(status, "tod_start")

    t_tot   = 0.d0
    call wall_time(t5)

    if (self%myid == 0) then
       write(*,*) 'delta in 1=', real(delta(:,1,1),sp)
       write(*,*) 'delta in 2=', real(delta(:,1,2),sp)
    end if

    ! Set up full-sky map structures
    call wall_time(t1)
    correct_sl      = .false.
    chisq_threshold = 7.d0
    n_main_iter     = 3
    !chisq_threshold = 20000.d0 
    !this ^ should be 7.d0, is currently 2000 to debug sidelobes
    ndet            = self%ndet
    ndelta          = size(delta,3)
    nside           = map_out%info%nside
    nmaps           = map_out%info%nmaps
    npix            = 12*nside**2
    ncol            = nmaps+ndet-1
    n_A             = nmaps*(nmaps+1)/2 + 4*(ndet-1)
    allocate(A_map(n_A,0:npix-1))
    allocate(A_abscal(self%ndet), b_abscal(self%ndet))
    allocate(smap_sky(0:ndet, ndelta)) 
    allocate(chisq_S(ndet,2))
    allocate(slist(self%nscan))
    nout = 1; if (self%output_all) nout = 7
    allocate(outmaps(nout), b_map(nout,ncol,0:npix-1))
    do i = 1, nout
       outmaps(i)%p => comm_map(map_out%info)
    end do

    call int2string(chain, ctext)
    call int2string(iter, samptext)
    prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
    postfix = '_c' // ctext // '_k' // samptext // '.fits'

!!$    do j = 1, map_in(1,1)%p%info%nmaps
!!$       do i = 0, map_in(1,1)%p%info%np-1
!!$          map_in(1,1)%p%map(i,j) = sum(map_in(:,1)%p%map(i,j))/4 
!!$       end do
!!$    end do
!!$    call map_in(1,1)%p%writeFITS("refmap.fits")
!!$    call mpi_finalize(ierr)
!!$    stop

    ! Distribute fullsky maps
    do j = 1, ndelta
       do i = 1, self%ndet
          map_in(i,j)%p%map = 1d-6 * map_in(i,j)%p%map ! uK to K
          call init_shared_2d_sp(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, [nmaps,npix], smap_sky(i,j))
          call sync_shared_2d_sp_map(smap_sky(i,j), map_in(i,j)%p%info%pix, &
               & transpose(map_in(i,j)%p%map))
          call mpi_win_fence(0, smap_sky(i,j)%win, ierr)
       end do
       call init_shared_2d_sp(self%myid_shared, self%comm_shared, &
            & self%myid_inter, self%comm_inter, [nmaps,npix], smap_sky(0,j))
       if (self%myid_shared == 0) then
          smap_sky(0,j)%a = smap_sky(1,j)%a
          do i = 2, self%ndet
             smap_sky(0,j)%a = smap_sky(0,j)%a + smap_sky(i,j)%a 
          end do
          smap_sky(0,j)%a = smap_sky(0,j)%a/self%ndet
       end if
       call mpi_win_fence(0, smap_sky(0,j)%win, ierr)
    end do

    ! Set up shared processing mask
    call init_shared_1d_int(self%myid_shared, self%comm_shared, &
         & self%myid_inter, self%comm_inter, [npix], sprocmask)
    call sync_shared_1d_int_map(sprocmask, self%procmask%info%pix, &
         & nint(self%procmask%map(:,1)))
    call init_shared_1d_int(self%myid_shared, self%comm_shared, &
         & self%myid_inter, self%comm_inter, [npix], sprocmask2)
    call sync_shared_1d_int_map(sprocmask2, self%procmask2%info%pix, &
         & nint(self%procmask2%map(:,1)))
    call wall_time(t2); t_tot(9) = t2-t1


    ! Compute far sidelobe Conviqt structures
    call wall_time(t1)
!!$map_in(1,1)%p%map = 0
!!$if (self%myid == 32) map_in(1,1)%p%map(5200,1) = 1
!!$call map_in(1,1)%p%writeFITS("in.fits")
    if (self%myid == 0) write(*,*) 'Precomputing sidelobe convolved sky'
    do i = 1, self%ndet
       if (.not. correct_sl) exit

       !call map_in(i,1)%p%remove_MDpoles() !remove mono and dipole components
       !call map_in(i,1)%p%writeFITS('nodp.fits')
       call map_in(i,1)%p%YtW()  ! Compute sky a_lms
       self%slconv(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
            & self%myid_inter, self%comm_inter, 128, 256, 3, 256, &
            & self%slbeam(i)%p, map_in(i,1)%p, 2)
    end do
    call wall_time(t2); t_tot(13) = t2-t1

    call update_status(status, "tod_init")

    ! Compute output map and rms
    call wall_time(t3)
    do_oper             = .true.
    do main_iter = 1, n_main_iter

       call wall_time(t7)

       if (self%myid == 0) write(*,*) '  Performing main iteration = ', main_iter
          
       ! Select operations for current iteration
       do_oper(bin_map)      = (main_iter == n_main_iter  )
       do_oper(sel_data)     = (main_iter == n_main_iter  ) .and.       first
       do_oper(calc_chisq)   = (main_iter == n_main_iter  ) 
       do_oper(prep_acal)    = (main_iter == n_main_iter-1) .and. .not. first
       do_oper(samp_acal)    = (main_iter == n_main_iter  ) .and. .not. first
       do_oper(prep_bp)      = (main_iter == n_main_iter-2) .and. .not. first
       do_oper(samp_bp)      = (main_iter == n_main_iter-1) .and. .not. first
       do_oper(prep_G)       = .false. !(main_iter == n_main_iter-3) .and. .not. first
       do_oper(samp_G)       = .false. !(main_iter == n_main_iter-2) .and. .not. first
       do_oper(samp_N)       = .true.
       do_oper(samp_mono)    = do_oper(bin_map)             .and. .not. first
       do_oper(sub_sl)       = correct_sl
       do_oper(output_slist) = mod(iter, 10) == 0

       ! Perform pre-loop operations
       if (do_oper(bin_map) .or. do_oper(prep_bp)) then
          A_map = 0.d0; b_map = 0.d0         
       end if

       if (do_oper(samp_mono)) then
          allocate(b_mono(nmaps,0:npix-1,ndet), sys_mono(nmaps,nmaps+ndet,0:self%info%np-1))
          b_mono = 0.d0
       end if

       if (do_oper(sel_data)) then
          naccept = 0; ntot    = 0
       end if

       if (do_oper(prep_acal)) then
          A_abscal = 0.d0; b_abscal = 0.d0
       end if

       if (do_oper(samp_acal)) then
          call mpi_barrier(self%comm, ierr)
          call wall_time(t1)
          call self%sample_absgain_from_orbital(A_abscal, b_abscal)
          call wall_time(t2); t_tot(16) = t_tot(16) + t2-t1
       end if

       if (do_oper(samp_bp)) then
          call wall_time(t1)
          call self%sample_bp(delta, smap_sky, handle, chisq_S)
          call wall_time(t2); t_tot(17) = t_tot(17) + t2-t1
       end if

       if (do_oper(samp_G)) then
          call wall_time(t1)
          call self%sample_smooth_gain()
          call wall_time(t2); t_tot(4) = t_tot(4) + t2-t1
       end if


       call wall_time(t8); t_tot(19) = t_tot(19) + t8-t7

       ! Perform main analysis loop
       do i = 1, self%nscan

          call update_status(status, "tod_loop1")
          call wall_time(t7)

          if (.not. any(self%scans(i)%d%accept)) cycle

          ! Short-cuts to local variables
          call wall_time(t1)
          ntod = self%scans(i)%ntod
          ndet = self%ndet
          
          ! Set up local data structure for current scan
          allocate(d_calib(nout,ntod, ndet))      ! Calibrated TOD in uKcmb
          allocate(n_corr(ntod, ndet))       ! Correlated noise in V
          allocate(s_sl(ntod, ndet))         ! Sidelobe in uKcm
          allocate(s_sky(ntod, ndet))        ! Sky signal in uKcmb
          allocate(s_sky_prop(ntod, ndet))   ! Sky signal in uKcmb
          allocate(s_sb(ntod, ndet))         ! Signal minus mean
          allocate(s_sb_prop(ntod, ndet))    ! Signal minus mean
          allocate(s_orb(ntod, ndet))        ! Orbital dipole in uKcmb
          allocate(s_mono(ntod, ndet))       ! Monopole correction in uKcmb
          allocate(s_buf(ntod, ndet))       ! Buffer
          allocate(s_tot(ntod, ndet))       ! Sum of all sky compnents
          allocate(mask(ntod, ndet))         ! Processing mask in time
          allocate(mask2(ntod, ndet))        ! Processing mask in time
          allocate(pix(ntod, ndet))          ! Decompressed pointing
          allocate(psi(ntod, ndet))          ! Decompressed pol angle
          allocate(flag(ntod, ndet))         ! Decompressed flags
          
          ! Initializing arrays to zero
          n_corr      = 0.d0
          s_sl        = 0.d0
          s_sky       = 0.d0
          s_sky_prop  = 0.d0
          s_orb       = 0.d0
          s_mono      = 0.d0
          call wall_time(t2); t_tot(18) = t_tot(18) + t2-t1
          
          ! --------------------
          ! Analyze current scan
          ! --------------------
         
          ! Decompress pointing, psi and flags for current scan
          call wall_time(t1)
          do j = 1, ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             call self%decompress_pointing_and_flags(i, j, pix(:,j), &
                  & psi(:,j), flag(:,j))
          end do
          !call validate_psi(self%scanid(i), psi)
          call wall_time(t2); t_tot(11) = t_tot(11) + t2-t1
          call update_status(status, "tod_decomp")
          
          ! Construct sky signal template
          call wall_time(t1)
          if (do_oper(bin_map) .or. do_oper(prep_bp)) then 
             call self%project_sky(smap_sky(:,1), pix, psi, flag, &
                  & sprocmask%a, i, s_sky, mask, s_sb=s_sb)  
          else 
             call self%project_sky(smap_sky(:,1), pix, psi, flag, &
                  & sprocmask%a, i, s_sky, mask)
          end if
          if (do_oper(prep_bp)) then
             call self%project_sky(smap_sky(:,2), pix, psi, flag, &
                  & sprocmask%a, i, s_sky_prop, mask, s_sb=s_sb_prop)  
          end if
          call wall_time(t2); t_tot(1) = t_tot(1) + t2-t1
          call update_status(status, "tod_project")
          
          ! Construct orbital dipole template
          call wall_time(t1)
          call self%compute_orbital_dipole(i, pix, s_orb)
          call wall_time(t2); t_tot(2) = t_tot(2) + t2-t1
          call update_status(status, "tod_orb")
           
          ! Construct sidelobe template 
          call wall_time(t1)
          if (do_oper(sub_sl)) then
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                call self%construct_sl_template(self%slconv(j)%p, i, &
                     & nside, pix(:,j), psi(:,j), s_sl(:,j), &
                     & self%mbang(j)+self%polang(j))
             end do
          end if
          call wall_time(t2); t_tot(12) = t_tot(12) + t2-t1

          ! Construct monopole correction template 
          call wall_time(t1)
          do j = 1, ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             s_mono(:,j) = self%mono(j)
          end do
          s_tot = s_sky + s_sl + s_orb + s_mono
          call wall_time(t2); t_tot(1) = t_tot(1) + t2-t1

          ! Fit correlated noise
          if (do_oper(samp_N)) then
             call wall_time(t1)
             call self%sample_n_corr(handle, i, mask, s_tot, n_corr)
             call wall_time(t2); t_tot(3) = t_tot(3) + t2-t1
          end if
             
          ! Compute noise spectrum
          if (do_oper(samp_N_par)) then
             call wall_time(t1)
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                call self%sample_noise_psd(handle, i, j, mask(:,j), s_tot(:,j),&
                     & n_corr(:,j))
             end do
             call wall_time(t2); t_tot(6) = t_tot(6) + t2-t1
          end if

          ! Fit gain 
          if (do_oper(prep_G)) then
             call wall_time(t1)
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                call self%sample_gain(handle, j, i, n_corr(:, j), mask(:,j), &
                     & s_tot(:, j))
             end do
             call wall_time(t2); t_tot(4) = t_tot(4) + t2-t1
          end if

          ! Compute chisquare
          if (do_oper(calc_chisq)) then
             call wall_time(t1)
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                s_buf(:,j) =  s_sl(:,j) + s_orb(:,j) + s_mono(:,j)
                call self%compute_chisq(i, j, mask(:,j), s_sky(:,j), &
                     & s_buf(:,j), n_corr(:,j))
             end do
             call wall_time(t2); t_tot(7) = t_tot(7) + t2-t1
          end if

          ! Select data
          if (do_oper(sel_data)) then
             call wall_time(t1)
             do j = 1, ndet
                ntot= ntot + 1
                if (.not. self%scans(i)%d(j)%accept) cycle
                if (count(iand(flag(:,j),6111248) .ne. 0) > 0.1*ntod) then
                   !write(*,fmt='(a,i8,i5,a,f12.1)') 'Reject scan, det = ', &
                    !    & self%scanid(i), j, ', more than 10% flagged samples'
                   self%scans(i)%d(j)%accept = .false.
                else if (abs(self%scans(i)%d(j)%chisq) > chisq_threshold .or. &
                & isNaN(self%scans(i)%d(j)%chisq)) then
                   write(*,fmt='(a,i8,i5,a,f12.1)') 'Reject scan, det = ', &
                        & self%scanid(i), j, ', chisq = ', &
                        & self%scans(i)%d(j)%chisq
                   self%scans(i)%d(j)%accept = .false.
                   cycle
                else
                   naccept = naccept + 1
                end if
             end do
             call wall_time(t2); t_tot(15) = t_tot(15) + t2-t1
          end if

          ! Compute chisquare for bandpass fit
          if (do_oper(prep_bp)) then
             call wall_time(t1)
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                s_buf(:,j) =  s_tot(:,j) - s_sky(:,j) !s_sl(:,j) + s_orb(:,j) + s_mono(:,j)
                call self%compute_chisq(i, j, mask(:,j), s_sky(:,j), &
                     & s_buf(:,j), n_corr(:,j),  s_sky_prop(:,j))
             end do
             call wall_time(t2); t_tot(7) = t_tot(7) + t2-t1
          end if

          ! Prepare for absolute calibration
          if (do_oper(prep_acal)) then
             call wall_time(t1)
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                s_buf(:,j) =  s_tot(:,j) - s_orb(:,j) !s_sky(:,j) + s_sl(:,j) + s_mono(:,j)
                call self%accumulate_absgain_from_orbital(i, j, mask(:,j),&
                     & s_buf(:,j), s_orb(:,j), n_corr(:,j), &
                     & A_abscal(j), b_abscal(j))
             end do
             call wall_time(t2); t_tot(14) = t_tot(14) + t2-t1
          end if

          ! Compute binned map 
          if (do_oper(bin_map) .or. do_oper(prep_bp)) then
             call wall_time(t1)
             if (self%output_all .and. do_oper(bin_map)) then
                do j = 1, ndet
                   inv_gain       = 1.d0 / self%scans(i)%d(j)%gain
                   d_calib(1,:,j) = (self%scans(i)%d(j)%tod - n_corr(:,j)) * &
                        & inv_gain - s_tot(:,j) + s_sky(:,j) - s_sb(:,j)
                   d_calib(2,:,j) = d_calib(1,:,j) - s_sky(:,j) + s_sb(:,j) ! Residual
                   d_calib(3,:,j) = n_corr(:,j) * inv_gain
                   d_calib(4,:,j) = s_sb(:,j)                         
                   d_calib(5,:,j) = s_mono(:,j)
                   d_calib(6,:,j) = s_orb(:,j)
                   d_calib(7,:,j) = s_sl(:,j)
                end do
             else if (do_oper(prep_bp)) then
                do j = 1, ndet
                   inv_gain       = 1.d0 / self%scans(i)%d(j)%gain
                   d_calib(1,:,j) = (self%scans(i)%d(j)%tod - n_corr(:,j)) * &
                        & inv_gain - s_tot(:,j) + s_sky(:,j) - s_sb(:,j)
                   d_calib(2,:,j) = d_calib(1,:,j) + s_sb(:,j) - s_sb_prop(:,j)
                end do
             else
                do j = 1, ndet
                   inv_gain       = 1.d0 / self%scans(i)%d(j)%gain
                   d_calib(1,:,j) = (self%scans(i)%d(j)%tod - n_corr(:,j)) * &
                        & inv_gain - s_tot(:,j) + s_sky(:,j) - s_sb(:,j)
                end do
             end if
             call wall_time(t2); t_tot(5) = t_tot(5) + t2-t1

             call wall_time(t1)
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                if (do_oper(prep_bp)) then
                   call self%compute_binned_map(d_calib(1:2,:,j), pix(:,j), &
                        & psi(:,j), flag(:,j), A_map, b_map(1:2,:,:), i, j, .true.)
                else if (do_oper(samp_mono)) then
                   call self%compute_binned_map(d_calib(:,:,j), pix(:,j), &
                        & psi(:,j), flag(:,j), A_map, b_map, i, j, .false., b_mono(:,:,j))
                else
                   call self%compute_binned_map(d_calib(:,:,j), pix(:,j), &
                        & psi(:,j), flag(:,j), A_map, b_map, i, j, .false.)
                end if
             end do
             call wall_time(t2); t_tot(8) = t_tot(8) + t2-t1
             !if (self%myid == 0) write(*,*) 'bin = ', t2-t1
          end if

          ! Clean up
          call wall_time(t1)
          deallocate(n_corr, s_sl, s_sky, s_orb, s_tot, d_calib)
          deallocate(mask, mask2, pix, psi, flag, s_sb, s_sky_prop, s_sb_prop, s_buf, s_mono)
          call wall_time(t2); t_tot(18) = t_tot(18) + t2-t1

          call wall_time(t8); t_tot(19) = t_tot(19) + t8-t7

          if (do_oper(output_slist)) then
             self%scans(i)%proctime   = self%scans(i)%proctime   + t8-t7
             self%scans(i)%n_proctime = self%scans(i)%n_proctime + 1
             if (main_iter == n_main_iter) then
                write(slist(i),*) self%scanid(i), '"',trim(self%hdfname(i)), &
                     & '"', real(self%scans(i)%proctime/self%scans(i)%n_proctime,sp)
             end if
          end if
          call update_status(status, "tod_loop2")

       end do

       ! Output total chisquare
       call wall_time(t7)
       if (.false.) then
          do j = 1, ndet
             chisq_tot = self%get_total_chisq(j)
             if (self%myid == 0) then
                write(*,fmt='(a,i4,a,i4,a,f12.1)') 'iter = ', main_iter, &
                     & ', det = ', j, ' -- chisq = ', chisq_tot  
             end if
          end do
       end if
 
       if (do_oper(prep_bp)) then
          call mpi_barrier(self%comm, ierr)
          call wall_time(t1)
          Sfilename = trim(prefix) // 'Smap'// trim(postfix) 
          call self%finalize_binned_map(handle, A_map, b_map(1:2,:,:), outmaps, &
               & rms_out, chisq_S, Sfile=Sfilename, mask=sprocmask2%a)
          call wall_time(t2); t_tot(10) = t_tot(10) + t2-t1
       end if
       call update_status(status, "tod_prepbp")

    end do
    call wall_time(t4)


    ! Output latest scan list with new timing information
    if (do_oper(output_slist)) then
       call wall_time(t1)
       call self%output_scan_list(slist)
       call wall_time(t2); t_tot(20) = t_tot(20) + t2-t1
    end if

    ! Solve combined map, summed over all pixels 
    call mpi_barrier(self%comm, ierr)
    call wall_time(t1)
    if (do_oper(samp_mono)) then
       condmap => comm_map(self%info)
       call self%finalize_binned_map(handle, A_map, b_map, outmaps, rms_out, b_mono=b_mono, sys_mono=sys_mono, condmap=condmap)
       call condmap%writeFITS("cond.fits")
       call condmap%dealloc()
    else
       call self%finalize_binned_map(handle, A_map, b_map, outmaps, rms_out)
    end if
    map_out%map = outmaps(1)%p%map

    ! Sample monopole coefficients and update output map
    if (do_oper(samp_mono)) then
       call self%sample_mono(handle, sys_mono, outmaps(2)%p, rms_out, &
            & self%procmask)
    end if

    ! Update bandpass parameters
    self%bp_delta = delta(:,:,1)

    ! Output maps to disk
    call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
    call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))

    if (self%output_all) then
       call outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
       call outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
       call outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix))
       call outmaps(5)%p%writeFITS(trim(prefix)//'mono'//trim(postfix))
       call outmaps(6)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
       call outmaps(7)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
    end if

    if (first) then
       call mpi_reduce(ntot,    i, 1, MPI_INTEGER, MPI_SUM, &
            & self%numprocs/2, self%info%comm, ierr)
       ntot = i
       call mpi_reduce(naccept, i, 1, MPI_INTEGER, MPI_SUM, &
            & self%numprocs/2, self%info%comm, ierr)
       naccept = i
    end if
    call wall_time(t2); t_tot(10) = t_tot(10) + t2-t1

    call wall_time(t6)
    if (self%myid == self%numprocs/2) then
       write(*,*) '  Time dist sky   = ', t_tot(9)
       write(*,*) '  Time sl precomp = ', t_tot(13)
       write(*,*) '  Time decompress = ', t_tot(11)
       write(*,*) '  Time alloc      = ', t_tot(18)
       write(*,*) '  Time project    = ', t_tot(1)
       write(*,*) '  Time orbital    = ', t_tot(2)
       write(*,*) '  Time sl interp  = ', t_tot(12)
       write(*,*) '  Time ncorr      = ', t_tot(3)
       write(*,*) '  Time gain       = ', t_tot(4)
       write(*,*) '  Time absgain    = ', t_tot(14)
       write(*,*) '  Time sel data   = ', t_tot(15)
       write(*,*) '  Time clean      = ', t_tot(5)
       write(*,*) '  Time noise      = ', t_tot(6)
       write(*,*) '  Time samp abs   = ', t_tot(16)
       write(*,*) '  Time samp bp    = ', t_tot(17)
       write(*,*) '  Time chisq      = ', t_tot(7)
       write(*,*) '  Time bin        = ', t_tot(8)
       write(*,*) '  Time scanlist   = ', t_tot(20)
       write(*,*) '  Time final      = ', t_tot(10)
       if (first) then
          write(*,*) '  Time total      = ', t6-t5, &
               & ', accept rate = ', real(naccept,sp) / ntot
       else
          write(*,*) '  Time total      = ', t6-t5, sum(t_tot(1:18))
       end if
    end if

    !stop

    ! Clean up temporary arrays
    deallocate(A_map, b_map)
    deallocate(A_abscal, b_abscal, chisq_S)
    if (allocated(b_mono)) deallocate(b_mono)
    if (allocated(sys_mono)) deallocate(sys_mono)
    do i = 1, nout
       call outmaps(i)%p%dealloc
    end do
    deallocate(outmaps)

    call dealloc_shared_1d_int(sprocmask)
    call dealloc_shared_1d_int(sprocmask2)
    do j = 1, ndelta
       do i = 0, self%ndet
          call dealloc_shared_2d_sp(smap_sky(i,j))
       end do
    end do
    deallocate(smap_sky)

    do i = 1, self%ndet
       if (.not. correct_sl) exit
       call self%slconv(i)%p%dealloc()
    end do

    call update_status(status, "tod_end")

    ! Parameter to check if this is first time routine has been 
    first = .false.

  end subroutine process_LFI_tod


  !**************************************************
  !             Sub-process routines
  !**************************************************
  subroutine decompress_pointing_and_flags(self, scan, det, pix, psi, flag)
    implicit none
    class(comm_LFI_tod),                intent(in)  :: self
    integer(i4b),                       intent(in)  :: scan, det
    integer(i4b),        dimension(:),  intent(out) :: pix, psi, flag

    integer(i4b) :: i, j
    character(len=6) :: stext

    call huffman_decode(self%scans(scan)%hkey, self%scans(scan)%d(det)%pix,  pix)
    call huffman_decode(self%scans(scan)%hkey, self%scans(scan)%d(det)%psi,  psi)
    call huffman_decode(self%scans(scan)%hkey, self%scans(scan)%d(det)%flag, flag)

    do j = 2, self%scans(scan)%ntod
       pix(j)  = pix(j-1)  + pix(j)
       psi(j)  = psi(j-1)  + psi(j)
       flag(j) = flag(j-1) + flag(j)
    end do
    psi = modulo(psi,4096)

  end subroutine decompress_pointing_and_flags

  ! Sky signal template
  subroutine project_sky(self, map, pix, psi, flag, pmask, scan_id, &
       & s_sky, tmask, s_sb)
    implicit none
    class(comm_LFI_tod),                    intent(in)  :: self
    integer(i4b),        dimension(0:),     intent(in)  :: pmask
    type(shared_2d_sp),  dimension(0:),     intent(in)  :: map
    integer(i4b),        dimension(:,:),    intent(in)  :: pix, psi
    integer(i4b),        dimension(:,:),    intent(in)  :: flag
    integer(i4b),                           intent(in)  :: scan_id
    real(sp),            dimension(:,:),    intent(out) :: s_sky, tmask
    real(sp),            dimension(:,:),    intent(out), optional :: s_sb

    integer(i4b) :: i, ierr,det
    real(sp)     :: s

    ! s = T + Q * cos(2 * psi) + U * sin(2 * psi)
    ! T - temperature; Q, U - Stoke's parameters
    do det = 1, self%ndet
       if (.not. self%scans(scan_id)%d(det)%accept) cycle
       do i = 1, self%scans(scan_id)%ntod
          s_sky(i,det) = map(det)%a(1,pix(i,det)+1) + &
                       & map(det)%a(2,pix(i,det)+1) * self%cos2psi(psi(i,det)) + &
                       & map(det)%a(3,pix(i,det)+1) * self%sin2psi(psi(i,det))
          tmask(i,det) = pmask(pix(i,det)) 
          if (iand(flag(i,det),6111248) .ne. 0) tmask(i,det) = 0.
       end do
    end do

    if (present(s_sb)) then
       do det = 1, self%ndet
          if (.not. self%scans(scan_id)%d(det)%accept) cycle
          do i = 1, self%scans(scan_id)%ntod
             s =    map(0)%a(1,pix(i,det)+1) + &
                  & map(0)%a(2,pix(i,det)+1) * self%cos2psi(psi(i,det)) + &
                  & map(0)%a(3,pix(i,det)+1) * self%sin2psi(psi(i,det))
             s_sb(i,det)  = s_sky(i,det) - s
          end do
       end do
    end if

  end subroutine project_sky


  ! Compute map with white noise assumption from correlated noise 
  ! corrected and calibrated data, d' = (d-n_corr-n_temp)/gain 
  subroutine compute_orbital_dipole(self, ind, pix, s_orb)
    implicit none
    class(comm_LFI_tod),                 intent(in)  :: self
    integer(i4b),                        intent(in)  :: ind !scan nr/index
    integer(i4b),        dimension(:,:), intent(in)  :: pix
    real(sp),            dimension(:,:), intent(out) :: s_orb
    real(dp)             :: x, T_0, q, pix_dir(3), b, b_dot
    real(dp), parameter  :: h = 6.62607015d-34   ! Planck's constant [Js]
    integer(i4b)         :: i,j

    !T_0 = T_CMB*k_b/h !T_0 = T_CMB frequency
    !x = freq * 1.d9 / (2.d0*T_0) !freq is the center bandpass frequancy of the detector
    !q = x * (exp(2.d0*x)+1) / (exp(2.d0*x)-1) !frequency dependency of the quadrupole
    !b = sqrt(sum(self%scans(ind)%v_sun**2))/c !beta for the given scan

    do i = 1,self%ndet
       if (.not. self%scans(ind)%d(i)%accept) cycle
       do j=1,self%scans(ind)%ntod !length of the tod
          if(pix(j, i) < 0 .or. pix(j,i) > 12*(self%nside**2)) then
             !write(*,*) pix(j, i), self%nside
             cycle
          end if
          
          b_dot = dot_product(self%scans(ind)%v_sun, self%pix2vec(:,pix(j,i)))/c
          s_orb(j,i) = T_CMB  * b_dot !only dipole, 1.d6 to make it uK, as [T_CMB] = K
          !s_orb(j,i) = T_CMB  * 1.d6 * b_dot !only dipole, 1.d6 to make it uK, as [T_CMB] = K
          !s_orb(j,i) = T_CMB * 1.d6 * (b_dot + q*b_dot**2) ! with quadrupole
          !s_orb(j,i) = T_CMB * 1.d6 * (b_dot + q*((b_dot**2) - (1.d0/3.d0)*(b**2))) ! net zero monopole
       end do
   end do

  end subroutine compute_orbital_dipole

  ! Compute correlated noise term, n_corr
  ! TODO: Add fluctuation if operation == sample
  subroutine sample_n_corr(self, handle, scan, mask, s_sub, n_corr)
    implicit none
    class(comm_LFI_tod),               intent(in)     :: self
    type(planck_rng),                  intent(inout)  :: handle
    integer(i4b),                      intent(in)     :: scan
    real(sp),          dimension(:,:), intent(in)     :: mask, s_sub
    real(sp),          dimension(:,:), intent(out)    :: n_corr
    integer(i4b) :: i, j, k, l, n, nomp, ntod, ndet, err, omp_get_max_threads, meanrange, start, last
    integer*8    :: plan_fwd, plan_back
    real(sp)     :: sigma_0, alpha, nu_knee, nu, samprate, gain, mean
    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv
    real(sp),     allocatable, dimension(:) :: d_prime, diff
    
    ntod = self%scans(scan)%ntod
    ndet = self%ndet
    nomp = omp_get_max_threads()
    
    meanrange = 20
    
    n = ntod + 1

    call sfftw_init_threads(err)
    call sfftw_plan_with_nthreads(nomp)

    allocate(dt(2*ntod), dv(0:n-1))
    call sfftw_plan_dft_r2c_1d(plan_fwd,  2*ntod, dt, dv, fftw_estimate + fftw_unaligned)
    call sfftw_plan_dft_c2r_1d(plan_back, 2*ntod, dv, dt, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)
    !$OMP PARALLEL PRIVATE(i,l,dt,dv,nu,sigma_0,alpha,nu_knee,d_prime,diff,start,last)
    allocate(dt(2*ntod), dv(0:n-1))
    allocate(d_prime(ntod))
    allocate(diff(ntod+1))
    diff = 0.d0
    !$OMP DO SCHEDULE(guided)
    do i = 1, ndet
       if (.not. self%scans(scan)%d(i)%accept) cycle
       gain = self%scans(scan)%d(i)%gain  ! Gain in V / K
       d_prime(:) = self%scans(scan)%d(i)%tod(:) - S_sub(:,i) * gain
      
       if(isNaN(sum(d_prime))) then
         !write(*,*) 'dprime', sum(self%scans(scan)%d(i)%tod), sum(S_sl(:,i)), sum(S_sky(:,i)), sum(S_orb(:,i)), gain
       end if
 
       
       ! if (i == 4 .and. scan == 1) then
       !    open(23, file="d_prime1.unf", form="unformatted")
       !    write(23) d_prime
       !    close(23)
       ! end if
       ! if (i == 4 .and. scan == 1) then
       !    open(230, file="mask.unf", form="unformatted")
       !    write(230) mask(:,i)
       !    close(230)
       ! end if


       sigma_0 = self%scans(scan)%d(i)%sigma0 * gain
       do j = 1,ntod
          if (mask(j,i) == 0.d0) then
             start = max(j - meanrange, 1)
             last = min(j + meanrange, ntod)
             if (sum(mask(start:last, i)) > 2) then
                mean = sum(d_prime(start:last) * mask(start:last, i)) / sum(mask(start:last, i))
             else
                start = max(j - meanrange * 4, 1)
                last = min(j + meanrange * 4, ntod)
!                write(*,*) "wider mean area needed", sum(mask(start:last, i))
                if(sum(mask(start:last, i)) == 0.d0) then 
                !  write(*,*) sum(d_prime(start:last)), sum(mask(start:last, i)), start, last, i, sum(S_sl(:,i))
                  mean = 0.d0
                else
                  mean = sum(d_prime(start:last) * mask(start:last, i)) / sum(mask(start:last, i))
                end if
                if(isNaN(mean)) then
                  mean = 0.d0
                end if
             end if
             
             !if(isNaN(d_prime(j)) .or. isNaN(mean)) then
             !  d_prime(j) = 0.d0
             !  mean = 0.d0
             !  write(*,*) 'setting d_prime, mean to', d_prime(j), mean, sigma_0
             !end if
             if (abs(d_prime(j) - mean) > 0.d0 * sigma_0) then 
                d_prime(j) = mean
             end if
          end if
       end do
       
       ! if (i == 4 .and. scan == 1) then
       !    open(23, file="d_prime2.unf", form="unformatted")
       !    write(23) d_prime
       !    close(23)
       ! end if

       ! where (mask(:,i) == 1.)
       !    d_prime(:) = self%scans(scan)%d(i)%tod(:) - S_sl(:,i) - (S_sky(:,i) + S_orb(:,i)) * gain
       ! elsewhere
       !    ! Do gap filling, placeholder
       !    sigma_0 = self%scans(scan)%d(i)%sigma0
       !    d_prime(:) = self%scans(scan)%d(i)%tod(:) - S_sl(:,i) - (S_sky(:,i) + S_orb(:,i)) * gain 
       ! end where

       dt(1:ntod)           = d_prime(:)
       dt(2*ntod:ntod+1:-1) = dt(1:ntod)
       if(isNaN(sum(dt))) then
         !write(*,*) sum(dt), sum(d_prime(:))
       end if

       call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
       samprate = self%samprate
       sigma_0  = self%scans(scan)%d(i)%sigma0
       alpha    = self%scans(scan)%d(i)%alpha
       nu_knee  = self%scans(scan)%d(i)%fknee
       do l = 1, n-1                                                      
          nu = l*(samprate/2)/(n-1)
          dv(l) = dv(l) * 1.d0/(1.d0 + (nu/(nu_knee))**(-alpha))
       end do
       call sfftw_execute_dft_c2r(plan_back, dv, dt)
       dt          = dt / (2*ntod)
       n_corr(:,i) = dt(1:ntod)
       if(isNaN(sum(n_corr(:,i)))) then
         !write(*,*) sum(dt), sum(dv)
       end if

       ! if (i == 1 .and. scan == 1) then
       !    open(22, file="tod.unf", form="unformatted")
       !    write(22) self%scans(scan)%d(i)%tod(:)
       !    close(22)
       ! end if

       ! if (i == 1 .and. scan == 1) then
       !    open(24, file="sky.unf", form="unformatted")
       !    write(24) S_sky(:,i) * gain
       !    close(24)
       ! end if

       
       ! if (i == 1 .and. scan == 1) then
       !    open(23, file="d_prime.unf", form="unformatted")
       !    write(23) d_prime
       !    close(23)
       ! end if


       ! if (i == 1 .and. scan == 1) then
       !    open(25, file="n_corr.unf", form="unformatted")
       !    write(25) n_corr(:,i)
       !    close(25)
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
    deallocate(diff)
    !$OMP END PARALLEL

    call sfftw_destroy_plan(plan_fwd)                                           
    call sfftw_destroy_plan(plan_back)                                          
  
  end subroutine sample_n_corr

  ! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  subroutine sample_gain(self, handle, det, scan_id, n_corr, mask, s_ref)
    implicit none
    class(comm_LFI_tod),               intent(inout)  :: self
    type(planck_rng),                  intent(inout)  :: handle
    integer(i4b),                      intent(in)     :: det, scan_id
    real(sp),            dimension(:), intent(in)     :: n_corr, mask, s_ref
    real(dp),            allocatable,  dimension(:)   :: d_only_wn
    real(dp)                                          :: curr_gain, ata
    real(dp)                                          :: curr_sigma

    allocate(d_only_wn(size(s_ref)))

    d_only_wn     = self%scans(scan_id)%d(det)%tod - n_corr
    ata           = sum(mask*s_ref**2)
    curr_gain     = sum(mask * d_only_wn * s_ref) / ata            
    curr_sigma    = self%scans(scan_id)%d(det)%sigma0 / sqrt(ata)  
    if (trim(self%operation) == 'optimize') then
       self%scans(scan_id)%d(det)%gain       = curr_gain
       self%scans(scan_id)%d(det)%gain_sigma = curr_sigma
    else
       ! TODO: Add fluctuation
       self%scans(scan_id)%d(det)%gain       = curr_gain
       self%scans(scan_id)%d(det)%gain_sigma = curr_sigma
    end if

    if(isNaN(self%scans(scan_id)%d(det)%gain)) then
      write(*,*) "Gain estimate", sum(s_ref), sum(d_only_wn), sum(self%scans(scan_id)%d(det)%tod), sum(n_corr), sum(mask), ata, curr_gain, curr_sigma
    end if

    deallocate(d_only_wn)

  end subroutine sample_gain


  ! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  subroutine sample_smooth_gain(self)
    implicit none
    class(comm_LFI_tod),               intent(inout)  :: self

    integer(i4b) :: i, j, k, ndet, nscan_tot, ierr, binsize, nbin, b1, b2, n
    real(dp)     :: g0
    real(dp), allocatable, dimension(:)     :: vals
    real(dp), allocatable, dimension(:,:,:) :: g
    
    ndet       = self%ndet
    nscan_tot  = self%nscan_tot
    binsize    = 100

    ! Collect all gain estimates on the root processor
    allocate(g(nscan_tot,ndet,2), vals(binsize))
    g = 0.d0
    do j = 1, ndet
       do i = 1, self%nscan
          k        = self%scanid(i)
          if (.not. self%scans(i)%d(j)%accept) cycle
          g(k,j,1) = self%scans(i)%d(j)%gain
          g(k,j,2) = self%scans(i)%d(j)%gain_sigma
       end do
    end do
    if (self%myid == 0) then
       call mpi_reduce(mpi_in_place, g, size(g), MPI_DOUBLE_PRECISION, MPI_SUM, &
            & 0, self%comm, ierr)
    else
       call mpi_reduce(g,            g, size(g), MPI_DOUBLE_PRECISION, MPI_SUM, &
            & 0, self%comm, ierr)
    end if


    ! Perform joint sampling/smoothing
    if (self%myid == 0) then
       nbin = nscan_tot/binsize+1
       do j = 1, ndet
          do i = 1, nbin
             b1 = (i-1)*binsize+1
             b2 = min(i*binsize,nscan_tot)
             if (count(g(b1:b2,j,1) > 0) <= 1) cycle
             
             !g0 = 0.d0
             n  = 0
             do k = b1, b2
                if (g(k,j,1) == 0) cycle
                !g0 = g0 + g(k,j,1)
                n  = n+1
                vals(n) = g(k,j,1)
             end do
             g(b1:b2,j,1) = median(vals(1:n))  !g0/n
          end do
       end do
    end if

    ! Distribute and update results
    call mpi_bcast(g, size(g),  MPI_DOUBLE_PRECISION, 0, self%comm, ierr)    
    do j = 1, ndet
       do i = 1, self%nscan
          k        = self%scanid(i)
          self%scans(i)%d(j)%gain = g(k,j,1)
       end do
    end do

    deallocate(g, vals)

  end subroutine sample_smooth_gain

  subroutine accumulate_absgain_from_orbital(self, scan, det, mask, s_sub, &
       & s_orb, n_corr, A_abs, b_abs)
    implicit none
    class(comm_LFI_tod),             intent(in)     :: self
    integer(i4b),                    intent(in)     :: scan, det
    real(sp),          dimension(:), intent(in)     :: mask, s_sub, s_orb
    real(sp),          dimension(:), intent(in)     :: n_corr
    real(dp),                        intent(inout)  :: A_abs, b_abs

    integer(i4b) :: i
    real(dp)     :: data          ! Cleaned data in units of K
    real(dp)     :: inv_sigmasq  ! Inverse variance in units of K**-2

    inv_sigmasq = (self%scans(scan)%d(det)%gain / &
         & self%scans(scan)%d(det)%sigma0)**2
    do i = 1, self%scans(scan)%ntod
       data = (self%scans(scan)%d(det)%tod(i)-n_corr(i))/ &
            & self%scans(scan)%d(det)%gain - s_sub(i) 
       A_abs = A_abs + S_orb(i) * inv_sigmasq * mask(i) * S_orb(i)
       b_abs = b_abs + S_orb(i) * inv_sigmasq * mask(i) * data
       if(abs(A_abs) < 1.d0) then
         !write(*,*) "accumulate gain", A_abs, b_abs, mask(i), inv_sigmasq, data, S_orb(i)
       end if
    end do

  end subroutine accumulate_absgain_from_orbital

  ! Sample absolute gain from orbital dipole alone 
  subroutine sample_absgain_from_orbital(self, A_abs, b_abs)
    implicit none
    class(comm_LFI_tod),               intent(inout)  :: self
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs

    integer(i4b) :: i, j, ierr
    real(dp), allocatable, dimension(:) :: A, b, dgain, sigma

    ! Collect contributions from all cores
    allocate(A(self%ndet), b(self%ndet), dgain(self%ndet), sigma(self%ndet))
    call mpi_reduce(A_abs, A, self%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & self%info%comm, ierr)
    call mpi_reduce(b_abs, b, self%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & self%info%comm, ierr)

    ! Compute gain update and distribute to all cores
    if (self%myid == 0) then
       sigma = 1.d0/sqrt(A)
       dgain = b/A
       if (.false.) then
          ! Add fluctuation term if requested
       end if
       do j = 1, self%ndet
          write(*,fmt='(a,i5,a,3f8.3)') 'Orb gain -- d = ', j, ', dgain = ', &
               & 100*(dgain(j)-1), 100*sigma(j), (dgain(j)-1)/sigma(j)
       end do
    end if
    call mpi_bcast(dgain, self%ndet,  MPI_DOUBLE_PRECISION, 0, &
         & self%info%comm, ierr)
    call mpi_bcast(sigma, self%ndet,  MPI_DOUBLE_PRECISION, 0, &
         & self%info%comm, ierr)

    do j = 1, self%ndet
       do i = 1, self%nscan
          if(isNaN(dgain(j))) then
            !write(*,*) "sample absgain", dgain(j), self%scans(i)%d(j)%gain
            dgain(j) = 1.d0
          end if
          self%scans(i)%d(j)%gain = dgain(j) * self%scans(i)%d(j)%gain
       end do
    end do

    deallocate(A, b, dgain, sigma)

  end subroutine sample_absgain_from_orbital
  
  !construct a sidelobe template in the time domain
  subroutine construct_sl_template(self, slconv, scan_id, nside, pix, psi, s_sl, polangle)
    implicit none
    class(comm_LFI_tod),                 intent(in)    :: self
    class(comm_conviqt),                 intent(in)    :: slconv
    integer(i4b),                        intent(in)    :: scan_id, nside
    integer(i4b),        dimension(:),   intent(in)    :: pix, psi
    real(dp),                            intent(in)   :: polangle
    real(sp),            dimension(:),   intent(out)   :: s_sl
    
    integer(i4b) :: j
    real(dp)     :: psi_

    do j=1, size(pix)
       psi_    = self%psi(psi(j))-polangle 
       s_sl(j) = slconv%interp(pix(j), psi_) 
    end do

  end subroutine construct_sl_template

  !compute the cleaned TOD from the computed TOD components
  subroutine compute_cleaned_tod(self, ntod, scan_num, s_sub, &
       & n_corr, d_calib)
    implicit none
    class(comm_LFI_tod),               intent(in)    :: self
    integer(i4b),                      intent(in)    :: ntod, scan_num
    real(sp),          dimension(:,:), intent(in)    :: s_sub
    real(sp),          dimension(:,:), intent(in)    :: n_corr
    real(sp),          dimension(:,:), intent(out)   :: d_calib

    integer(i4b) :: i

    !cleaned calibrated data = (rawTOD - corrNoise)/gain - orbitalDipole - sideLobes
    do i=1, self%ndet
       if (.not. self%scans(scan_num)%d(i)%accept) cycle
      d_calib(:,i) = (self%scans(scan_num)%d(i)%tod - n_corr(:,i)) / &
           & self%scans(scan_num)%d(i)%gain - s_sub(:,i)
    end do
  end subroutine compute_cleaned_tod


  ! Compute chisquare
  subroutine compute_chisq(self, scan, det, mask, s_sky, s_spur, &
       & n_corr, s_prop)
    implicit none
    class(comm_LFI_tod),             intent(inout)  :: self
    integer(i4b),                    intent(in)     :: scan, det
    real(sp),          dimension(:), intent(in)     :: mask,s_sky, s_spur
    real(sp),          dimension(:), intent(in)     :: n_corr
    real(sp),          dimension(:), intent(in), optional :: s_prop

    real(dp)     :: chisq, chisq_prop, d0, g
    integer(i4b) :: i, n

    chisq       = 0.d0
    chisq_prop  = 0.d0
    n           = 0
    g           = self%scans(scan)%d(det)%gain 
    do i = 1, self%scans(scan)%ntod
       if (mask(i) < 0.5) cycle 
       n     = n+1
       d0    = self%scans(scan)%d(det)%tod(i) - &
            & (g * s_spur(i) + n_corr(i))
       chisq = chisq + (d0 - g * s_sky(i))**2 
       if (present(s_prop)) then
          chisq_prop = chisq_prop + (d0 - g * s_prop(i))**2 
       end if
    end do
    chisq      = chisq      / self%scans(scan)%d(det)%sigma0**2

    if (present(s_prop)) then
       chisq_prop = chisq_prop / self%scans(scan)%d(det)%sigma0**2
       self%scans(scan)%d(det)%chisq_masked = chisq
       self%scans(scan)%d(det)%chisq_prop   = chisq_prop
    else
       self%scans(scan)%d(det)%chisq        = (chisq - n) / sqrt(2.d0*n)
    end if

    if(self%scans(scan)%d(det)%chisq > 2000.d0 .or. isNaN(self%scans(scan)%d(det)%chisq)) then
      !write(*,*) "chisq", scan, det, sum(mask), sum(s_sky), sum(s_sl), sum(s_orb), sum(n_corr)
    end if

  end subroutine compute_chisq

  ! Sample noise psd
  ! TODO: Add fluctuation term if operation == sample
  subroutine sample_noise_psd(self, handle, scan, det, mask, s_tot, n_corr)
    implicit none
    class(comm_LFI_tod),             intent(inout)  :: self
    type(planck_rng),                intent(in)     :: handle
    integer(i4b),                    intent(in)     :: scan, det
    real(sp),          dimension(:), intent(in)     :: mask, s_tot, n_corr
    
    integer(i4b) :: i, n
    real(dp)     :: s, res

    s = 0.d0
    n = 0
    do i = 1, self%scans(scan)%ntod-1
       if (any(mask(i:i+1) < 0.5)) cycle
       res = (self%scans(scan)%d(det)%tod(i) - &
         & (self%scans(scan)%d(det)%gain * s_tot(i) + &
         & n_corr(i)) - &
         & (self%scans(scan)%d(det)%tod(i+1) - &
         & (self%scans(scan)%d(det)%gain * s_tot(i+1) + &
         & n_corr(i+1))))/sqrt(2.)
       s = s + res**2
       n = n + 1
    end do

    self%scans(scan)%d(det)%sigma0 = sqrt(s/(n-1))

  end subroutine sample_noise_psd

  ! Compute map with white noise assumption from correlated noise 
  ! corrected and calibrated data, d' = (d-n_corr-n_temp)/gain 
  subroutine compute_binned_map(self, data, pix, psi, flag, A, b, scan, det, comp_S, b_mono)
    implicit none
    class(comm_LFI_tod),                      intent(in)    :: self
    integer(i4b),                             intent(in)    :: scan, det
    real(sp),            dimension(1:,1:),    intent(in)    :: data
    integer(i4b),        dimension(1:),       intent(in)    :: pix, psi, flag
    real(dp),            dimension(1:,0:),    intent(inout) :: A
    real(dp),            dimension(1:,1:,0:), intent(inout) :: b
    real(dp),            dimension(1:,0:),    intent(inout), optional :: b_mono
    logical(lgt),                             intent(in)    :: comp_S

    integer(i4b) :: i, j, t, pix_, off, nout
    real(dp)     :: psi_, inv_sigmasq

    nout        = size(b,dim=1)
    inv_sigmasq = (self%scans(scan)%d(det)%gain/self%scans(scan)%d(det)%sigma0)**2
    do t = 1, self%scans(scan)%ntod

       if (iand(flag(t),6111248) .ne. 0) cycle

       pix_    = pix(t)
       psi_    = psi(t)
       
       A(1,pix_) = A(1,pix_) + 1.d0                                  * inv_sigmasq
       A(2,pix_) = A(2,pix_) + self%cos2psi(psi_)                    * inv_sigmasq
       A(3,pix_) = A(3,pix_) + self%cos2psi(psi_)**2                 * inv_sigmasq
       A(4,pix_) = A(4,pix_) + self%sin2psi(psi_)                    * inv_sigmasq
       A(5,pix_) = A(5,pix_) + self%cos2psi(psi_)*self%sin2psi(psi_) * inv_sigmasq
       A(6,pix_) = A(6,pix_) + self%sin2psi(psi_)**2                 * inv_sigmasq

       do i = 1, nout
          b(i,1,pix_) = b(i,1,pix_) + data(i,t)                      * inv_sigmasq
          b(i,2,pix_) = b(i,2,pix_) + data(i,t) * self%cos2psi(psi_) * inv_sigmasq
          b(i,3,pix_) = b(i,3,pix_) + data(i,t) * self%sin2psi(psi_) * inv_sigmasq
       end do

       if (present(b_mono)) then
          b_mono(1,pix_) = b_mono(1,pix_) +                      inv_sigmasq
          b_mono(2,pix_) = b_mono(2,pix_) + self%cos2psi(psi_) * inv_sigmasq
          b_mono(3,pix_) = b_mono(3,pix_) + self%sin2psi(psi_) * inv_sigmasq
       end if

       if (comp_S .and. det < self%ndet) then
          off           = 6 + 4*(det-1)
          A(off+1,pix_) = A(off+1,pix_) + 1.d0               * inv_sigmasq 
          A(off+2,pix_) = A(off+2,pix_) + self%cos2psi(psi_) * inv_sigmasq
          A(off+3,pix_) = A(off+3,pix_) + self%sin2psi(psi_) * inv_sigmasq
          A(off+4,pix_) = A(off+4,pix_) + 1.d0               * inv_sigmasq
          do i = 1, nout
             b(i,det+3,pix_) = b(i,det+3,pix_) + data(i,t) * inv_sigmasq 
          end do
       end if

    end do

  end subroutine compute_binned_map


  subroutine finalize_binned_map(self, handle, A, b, outmaps, rms, chisq_S, Sfile, mask, b_mono, sys_mono, condmap)
    implicit none
    class(comm_LFI_tod),                  intent(in)    :: self
    type(planck_rng),                     intent(inout) :: handle
    real(dp),        dimension(1:,0:),    intent(in)    :: A
    real(dp),        dimension(1:,1:,0:), intent(in)    :: b
    class(map_ptr),  dimension(1:),       intent(inout) :: outmaps
    class(comm_map),                      intent(inout) :: rms
    real(dp),        dimension(1:,1:),    intent(out),   optional :: chisq_S
    character(len=*),                     intent(in),    optional :: Sfile
    integer(i4b),    dimension(0:),       intent(in),    optional :: mask
    real(dp),        dimension(1:,0:,1:), intent(inout), optional :: b_mono
    real(dp),        dimension(1:,1:,0:), intent(out),   optional :: sys_mono
    class(comm_map),                      intent(inout), optional :: condmap

    integer(i4b) :: i, j, k, l, nmaps, np, ierr, ndet, ncol, n_A, off
    integer(i4b) :: det, nout, np0, comm, myid, nprocs
    real(dp), allocatable, dimension(:,:)   :: A_inv
    real(dp), allocatable, dimension(:,:,:) :: b_tot, buffer_b
    real(dp), allocatable, dimension(:)     :: buffer, W, eta
    real(dp), allocatable, dimension(:,:)   :: A_tot, buffer_A
    integer(i4b), allocatable, dimension(:) :: pix
    class(comm_mapinfo), pointer :: info
    class(comm_map), pointer :: smap

    myid  = outmaps(1)%p%info%myid
    nprocs= outmaps(1)%p%info%nprocs
    comm  = outmaps(1)%p%info%comm
    np0   = outmaps(1)%p%info%np
    nout  = size(b,dim=1)
    nmaps = size(outmaps(1)%p%map, dim=2)
    ndet  = size(b,dim=2)-nmaps+1
    if (present(chisq_S)) then
       ncol  = nmaps+ndet-1
       n_A   = nmaps*(nmaps+1)/2 + 4*(ndet-1)
    else
       ncol  = nmaps
       n_A   = nmaps*(nmaps+1)/2
    end if

    ! Collect contributions from all cores
    allocate(A_tot(n_A,0:np0-1), b_tot(nout,ncol,0:np0-1), W(ncol), eta(nmaps))
    
    call update_status(status, "tod_final1")
    do i = 0, nprocs-1
       if (myid == i) np = np0
       call mpi_bcast(np, 1,  MPI_INTEGER, i, comm, ierr)
       allocate(pix(np))
       if (myid == i) pix = rms%info%pix
       call mpi_bcast(pix, np,  MPI_INTEGER, i, comm, ierr)
       call mpi_reduce(A(1:n_A,pix), A_tot, n_A*np, MPI_DOUBLE_PRECISION, MPI_SUM, i, comm, ierr)
       call mpi_reduce(b(1:nout,1:ncol,pix), b_tot, ncol*np*nout, MPI_DOUBLE_PRECISION, MPI_SUM, i, comm, ierr)
       if (present(b_mono)) then
          if (myid == i) then
             allocate(buffer_b(nmaps,0:np-1,ndet))
             call mpi_reduce(b_mono(1:,pix,1:), buffer_b, size(buffer_b), MPI_DOUBLE_PRECISION, MPI_SUM, i, comm, ierr)
             do j = 1, ndet
                sys_mono(:,nmaps+j,:) = buffer_b(:,:,j)
             end do
             deallocate(buffer_b)
          else
             call mpi_reduce(b_mono(1:,pix,1:), b_mono(1:,pix,1:), nmaps*np*ndet, MPI_DOUBLE_PRECISION, MPI_SUM, i, comm, ierr)
          end if
       end if
       deallocate(pix)
    end do
    call update_status(status, "tod_final2")


    ! Solve for local map and rms
    if (present(Sfile)) then
       info => comm_mapinfo(self%comm, self%info%nside, 0, ndet-1, .false.)
       smap => comm_map(info)
    end if
    allocate(A_inv(ncol,ncol))
    A_inv   = 0.d0
    !if (present(chisq_S)) smap%map = 0.d0
    if (present(chisq_S)) chisq_S = 0.d0
    do i = 0, np0-1
       if (all(b_tot(1,:,i) == 0.d0)) cycle

       A_inv(1,1) = A_tot(1,i)
       A_inv(2,1) = A_tot(2,i)
       A_inv(1,2) = A_inv(2,1)
       A_inv(2,2) = A_tot(3,i)
       A_inv(3,1) = A_tot(4,i)
       A_inv(1,3) = A_inv(3,1)
       A_inv(3,2) = A_tot(5,i)
       A_inv(2,3) = A_inv(3,2)
       A_inv(3,3) = A_tot(6,i)
       if (present(chisq_S)) then
          do det = 1, ndet-1
             off           = 6 + 4*(det-1)
             A_inv(1,    3+det) = A_tot(off+1,i)
             A_inv(3+det,1    ) = A_inv(1,3+det)
             A_inv(2,    3+det) = A_tot(off+2,i)
             A_inv(3+det,2    ) = A_inv(2,3+det)
             A_inv(3,    3+det) = A_tot(off+3,i)
             A_inv(3+det,3    ) = A_inv(3,3+det)
             A_inv(3+det,3+det) = A_tot(off+4,i)
          end do
       end if

       if (present(condmap)) then
          call get_eigenvalues(A_inv, W)
          condmap%map(i,1) = log10(max(abs(maxval(W)/minval(W)),1.d0))
       end if
       
       call invert_singular_matrix(A_inv, 1d-12)
       do k = 1, nout
          b_tot(k,1:ncol,i) = matmul(A_inv,b_tot(k,1:ncol,i))
       end do
       if (present(b_mono)) sys_mono(1:nmaps,1:nmaps,i) = A_inv(1:nmaps,1:nmaps)
       if (present(Sfile)) then
          do j = 1, ndet-1
             smap%map(i,j) = b_tot(1,nmaps+j,i) / sqrt(A_inv(nmaps+j,nmaps+j))
          end do
       end if
       if (present(chisq_S)) then
          if (mask(rms%info%pix(i+1)) == 0) cycle
          do j = 1, ndet-1
             if (A_inv(nmaps+j,nmaps+j) == 0.d0) cycle
             chisq_S(j,1) = chisq_S(j,1) + b_tot(1,nmaps+j,i)**2 / A_inv(nmaps+j,nmaps+j)
             chisq_S(j,2) = chisq_S(j,2) + b_tot(2,nmaps+j,i)**2 / A_inv(nmaps+j,nmaps+j)
          end do
       else
          do j = 1, nmaps
             rms%map(i,j) = sqrt(A_inv(j,j))  * 1.d6 ! uK
             do k = 1, nout
                outmaps(k)%p%map(i,j) = b_tot(k,j,i) * 1.d6 ! uK
             end do
             if (trim(self%operation) == 'sample') then
                ! Add random fluctuation
                call compute_hermitian_root(A_inv, 0.5d0)
                do l = 1, nmaps
                   eta(l) = rand_gauss(handle)
                end do
                outmaps(k)%p%map(i,:) = outmaps(k)%p%map(i,:) + matmul(A_inv,eta)
             end if
          end do
       end if
    end do

    if (present(chisq_S)) then
       if (myid == 0) then
          call mpi_reduce(mpi_in_place, chisq_S, size(chisq_S), &
               & MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
       else
          call mpi_reduce(chisq_S,      chisq_S, size(chisq_S), &
               & MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
       end if

       if (present(Sfile)) call smap%writeFITS(trim(Sfile))
    end if

    if (present(Sfile)) call smap%dealloc

    deallocate(A_inv,A_tot,b_tot, W, eta)

  end subroutine finalize_binned_map


!!$  subroutine sample_bp2(self, delta, map_sky, handle, chisq_S)
!!$    implicit none
!!$    class(comm_LFI_tod),               intent(inout)  :: self
!!$    real(dp), dimension(1:,1:),        intent(inout)  :: delta
!!$    real(dp), dimension(0:,1:,0:,1:),  intent(inout)  :: map_sky
!!$    type(planck_rng),                  intent(inout)  :: handle
!!$    real(dp), dimension(1:,1:),        intent(in)     :: chisq_S
!!$
!!$    integer(i4b) :: i, j, ierr
!!$    logical(lgt) :: accept
!!$    real(dp)     :: chisq_prop, chisq_curr, cp, cc, accept_rate
!!$
!!$    do j = 1, self%ndet
!!$       ! Summing chisq for current delta
!!$       chisq_curr = 0.d0
!!$       chisq_prop = 0.d0
!!$       do i = 1, self%nscan ! Sum chisq for all scans
!!$          if (.not. self%scans(i)%d(j)%accept) cycle
!!$          chisq_curr = chisq_curr + self%scans(i)%d(j)%chisq_masked 
!!$          chisq_prop = chisq_prop + self%scans(i)%d(j)%chisq_prop 
!!$       end do
!!$
!!$       call mpi_reduce(chisq_prop, cp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
!!$            & 0, self%info%comm, ierr)
!!$       call mpi_reduce(chisq_curr, cc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
!!$            & 0, self%info%comm, ierr)
!!$             
!!$       if (self%myid == 0) then
!!$          cp          = cp + chisq_S(j,2) 
!!$          cc          = cc + chisq_S(j,1)
!!$          accept_rate = exp(-0.5d0*(cp-cc))  
!!$          accept = (rand_uni(handle) < accept_rate)
!!$          write(*,*) 'S=',chisq_S(j,:)
!!$          write(*,fmt='(a,i3,a,f16.1,a,f10.1,a,f10.6)') "det = ", j, &
!!$               & ', c0 = ', cc, ', diff = ', cp-cc, ', delta = ', delta(j,1)
!!$       end if
!!$
!!$       ! Broadcast new saved data
!!$       call mpi_bcast(accept, 1,  MPI_LOGICAL, 0, self%info%comm, ierr)
!!$       if (accept) then
!!$          ! Set current to proposal
!!$          map_sky(:,:,j,1) = map_sky(:,:,j,2)
!!$          delta(j,1)       = delta(j,2)
!!$       end if
!!$    end do
!!$
!!$    ! Update mean
!!$    map_sky(:,:,0,:) = map_sky(:,:,1,:)
!!$    do i = 2, self%ndet
!!$       map_sky(:,:,0,:) = map_sky(:,:,0,:) + map_sky(:,:,i,:) 
!!$    end do
!!$    map_sky(:,:,0,:) = map_sky(:,:,0,:)/self%ndet
!!$    
!!$  end subroutine sample_bp2

  subroutine sample_bp(self, delta, smap_sky, handle, chisq_S)
    implicit none
    class(comm_LFI_tod),               intent(inout)  :: self
    real(dp), dimension(1:,1:,1:),     intent(inout)  :: delta
    type(shared_2d_sp), dimension(0:,1:),  intent(inout)  :: smap_sky
    type(planck_rng),                  intent(inout)  :: handle
    real(dp), dimension(1:,1:),        intent(in)     :: chisq_S

    integer(i4b) :: i, j, ierr
    logical(lgt) :: accept
    real(dp)     :: chisq_prop, chisq_curr, cp, cc, accept_rate

    if (self%myid == 0) then
       cp          = sum(chisq_S(:,2))
       cc          = sum(chisq_S(:,1))
       accept_rate = exp(-0.5d0*(cp-cc))  
       if (trim(self%operation) == 'optimize') then
          accept = cp <= cc
       else
          accept = (rand_uni(handle) < accept_rate)
       end if
       if (accept) then
          write(*,fmt='(a,f16.1,a,f10.1)') 'bp c0 = ', cp, &
          & ', diff = ', cp-cc
       else
          write(*,fmt='(a,f16.1,a,f10.1)') 'bp c0 = ', cc, &
               & ', diff = ', cp-cc
       end if
    end if

    ! Broadcast new saved data
    call mpi_bcast(accept, 1,  MPI_LOGICAL, 0, self%info%comm, ierr)
    if (accept) then
       ! Set current to proposal
       do i = 0, self%ndet
          if (self%myid_shared == 0) smap_sky(i,1)%a = smap_sky(i,2)%a
       end do
       delta(:,:,1) =  delta(:,:,2)
    end if
    
  end subroutine sample_bp

  function get_total_chisq(self, det)
    implicit none
    class(comm_LFI_tod), intent(in)  :: self
    integer(i4b),        intent(in)  :: det
    real(dp)                         :: get_total_chisq

    integer(i4b) :: i, ierr
    real(dp)     :: chisq, buffer

    chisq = 0.d0
    do i = 1, self%nscan ! Sum chisq for all scans
       if (.not. self%scans(i)%d(det)%accept) cycle
       chisq = chisq + self%scans(i)%d(det)%chisq 
    end do
    call mpi_reduce(chisq, buffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         & 0, self%info%comm, ierr)

    get_total_chisq = buffer

  end function get_total_chisq

  subroutine output_scan_list(self, slist)
    implicit none
    class(comm_LFI_tod),                           intent(in)    :: self
    character(len=512), allocatable, dimension(:), intent(inout) :: slist

    integer(i4b)     :: i, j, mpistat(MPI_STATUS_SIZE), unit, ns, ierr
    character(len=4) :: pid

    if (self%myid == 0) then
       call int2string(self%myid, pid)
       unit = getlun()
       open(unit,file='newlist_debug.txt', recl=512)
       write(unit,*) sum(self%nscanprproc)
       do i = 1, self%nscan
          write(unit,*) trim(slist(i))
       end do
       deallocate(slist)
       do j = 1, self%numprocs-1
          ns = self%nscanprproc(j)
          allocate(slist(ns))
          call mpi_recv(slist, 512*ns, MPI_CHARACTER, j, 98, self%comm, mpistat, ierr)
          do i = 1, ns
             write(unit,*) trim(slist(i))
          end do
          deallocate(slist)
       end do
       close(unit)
    else
       call mpi_send(slist, 512*self%nscan, MPI_CHARACTER, 0, 98, self%comm, ierr)
       deallocate(slist)
    end if
  end subroutine output_scan_list

  subroutine sample_mono(self, handle, sys_mono, res, rms, mask)
    implicit none
    class(comm_LFI_tod),                   intent(inout) :: self
    type(planck_rng),                      intent(inout) :: handle
    real(dp),         dimension(1:,1:,0:), intent(in)    :: sys_mono
    class(comm_map),                       intent(in)    :: res, rms, mask

    integer(i4b) :: i, j, k, nstep, nmaps, ndet, ierr
    real(dp)     :: t1, t2, s(3), b(3), sigma, chisq, chisq0, chisq_prop
    logical(lgt) :: accept
    real(dp), allocatable, dimension(:) :: mono0, mono_prop

    call wall_time(t1)

    nstep = 500
    sigma = 0.03d-6 ! RMS proposal in K
    nmaps = size(sys_mono, dim=1)
    ndet  = size(sys_mono, dim=2)-nmaps
    allocate(mono0(ndet), mono_prop(ndet))
 
    ! Initialize monopoles on existing values
    if (self%myid == 0) then
       do j = 1, ndet
          mono0(j) = self%mono(j)
       end do
    end if
    call mpi_bcast(mono0, ndet,  MPI_DOUBLE_PRECISION, 0, res%info%comm, ierr)

    ! Compute chisquare for old values; only include Q and U in evaluation; 
    ! add old correction to output map
    chisq = 0.d0
    do k = 0, res%info%np-1
       b = 0.d0
       do j = 1, ndet
          b = b + mono0(j) * sys_mono(:,nmaps+j,k)
       end do
       s = matmul(sys_mono(1:3,1:3,k),b) * 1.d6 ! uK
       if (mask%map(k,1) == 1.d0) then
          chisq = chisq + sum((res%map(k,2:3)-s(2:3))**2 / rms%map(k,2:3)**2)
       end if
    end do
    call mpi_reduce(chisq, chisq0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, res%info%comm, ierr)

    do i = 1, nstep

       ! Propose new monopoles; force zero average
       if (res%info%myid == 0) then
          do j = 1, ndet
             mono_prop(j) = mono0(j) + sigma * rand_gauss(handle)
          end do
          mono_prop = mono_prop - mean(mono_prop)
       end if
       call mpi_bcast(mono_prop, ndet,  MPI_DOUBLE_PRECISION, 0, res%info%comm, ierr)

       ! Compute chisquare for new values; only include Q and U in evaluation
       chisq = 0.d0
       do k = 0, res%info%np-1
          b = 0.d0
          do j = 1, ndet
             b = b + mono_prop(j) * sys_mono(:,nmaps+j,k)
          end do
          s = matmul(sys_mono(1:3,1:3,k),b) * 1.d6 ! uK
          if (mask%map(k,1) == 1) then
             chisq = chisq + sum((res%map(k,2:3)-s(2:3))**2 / rms%map(k,2:3)**2)
          end if
       end do
       call mpi_reduce(chisq, chisq_prop, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, res%info%comm, ierr)

       ! Apply Metropolis rule or maximize
       if (res%info%myid == 0) then
          if (trim(self%operation) == 'optimize') then
             accept = chisq_prop <= chisq0
          else
             accept = (rand_uni(handle) < exp(-0.5d0 * (chisq_prop-chisq0)))
          end if
!!$          write(*,*) 'Mono chisq0 = ', chisq0, ', delta = ', chisq_prop-chisq0
!!$          write(*,*) 'm0 = ', real(mono0,sp)
!!$          write(*,*) 'm1 = ', real(mono_prop,sp)
       end if
       call mpi_bcast(accept, 1,  MPI_LOGICAL, 0, res%info%comm, ierr)

       if (accept) then
          mono0  = mono_prop
          chisq0 = chisq_prop
       end if

    end do

    ! Update parameters
    self%mono = mono0

    call wall_time(t2)
    if (res%info%myid == 0) then
       write(*,*) 'Time for monopole sampling = ', t2-t1
    end if

    deallocate(mono0, mono_prop)

  end subroutine sample_mono

  subroutine validate_psi(scan, psi)
    implicit none
    integer(i4b),                 intent(in)    :: scan
    integer(i4b), dimension(:,:), intent(inout) :: psi

    integer(i4b) :: i, j, ntod, diff0, diff
    logical(lgt) :: ok
    character(len=6) :: stext


    ntod = size(psi, dim=1)

!!$    diff0 = psi(1,1)-psi(1,2)
!!$    diff  = psi(1,3)-psi(1,4)
!!$    do i = 1, ntod
!!$       if (psi(i,1) == 0) psi(i,1) = modulo(psi(i,2) + diff0,4096) 
!!$       if (psi(i,2) == 0) psi(i,2) = modulo(psi(i,1) - diff0,4096)
!!$       if (psi(i,3) == 0) psi(i,3) = modulo(psi(i,4) + diff, 4096) 
!!$       if (psi(i,4) == 0) psi(i,4) = modulo(psi(i,3) - diff, 4096) 
!!$    end do
!!$
!!$    return

    diff0 = smallest_angle(psi(1,1), psi(1,2))
    do i = 2, ntod
       diff = smallest_angle(psi(i,1), psi(i,2))
       ok = abs(diff-diff0) < 50 
       if (.not. ok) then
          write(*,*) 'Potential psi problem in scan = ', scan, ', det = 1 or 2', diff0, diff, i
          call int2string(scan, stext)
          open(58,file='error12_'//stext//'.txt')
          do j = 1, ntod
             diff = smallest_angle(psi(j,1), psi(j,2))
             write(58,*) j, psi(j,1), psi(j,2), diff
          end do
          close(58)
          exit
       end if
    end do

    diff0 = smallest_angle(psi(1,3), psi(1,4))
    do i = 2, ntod
       diff = smallest_angle(psi(i,3), psi(i,4))
       ok = abs(diff-diff0) < 10 
       if (.not. ok) then
          write(*,*) 'Potential psi problem in scan = ', scan, ', det = 3 or 4', diff0, diff, i
          call int2string(scan, stext)
          open(58,file='error34_'//stext//'.txt')
          do j = 1, ntod
             diff = smallest_angle(psi(j,3), psi(j,4))
             write(58,*) j, psi(j,3), psi(j,4), diff
          end do
          close(58)

          exit
       end if
    end do

  end subroutine validate_psi

  function smallest_angle(x, y)
    integer(i4b), intent(in) :: x, y
    integer(i4b)             :: smallest_angle

    integer(i4b) :: a, b

    a = modulo(x-y, 4096)
    b = modulo(y-x, 4096)
    if (a < b) then
       smallest_angle = -a
    else
       smallest_angle =  b
    end if

  end function smallest_angle

end module comm_tod_LFI_mod
