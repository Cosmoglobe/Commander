module comm_tod_LFI_mod
  use comm_tod_mod
  use comm_param_mod
  use comm_map_mod
  use comm_conviqt_mod
  use pix_tools
  use healpix_types  
  use comm_huffman_mod
  use comm_hdf_mod
  implicit none

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
     procedure     :: decompress_pointing_and_flags
     procedure     :: accumulate_absgain_from_orbital
     procedure     :: sample_absgain_from_orbital
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

    integer(i4b) :: i, nside_beam, lmax_beam, nmaps_beam
    logical(lgt) :: pol_beam
    type(hdf_file) :: h5_file

    ! Set up LFI specific parameters
    allocate(constructor)
    constructor%myid     = cpar%myid_chain
    constructor%comm     = cpar%comm_chain
    constructor%numprocs = cpar%numprocs
    constructor%info     => info

    ! Test code, just to be able to read a single file; need to settle on parameter structure
    !call constructor%get_scan_ids("data/filelist_1file.txt")
    call constructor%get_scan_ids("data/filelist_1year.txt")
    !call constructor%get_scan_ids("data/filelist.txt")
!    call constructor%get_scan_ids("data/filelist_2half.txt")

    constructor%procmask => comm_map(constructor%info, "data/test_procmask_expand20arcmin_n512.fits")
    constructor%procmask2 => comm_map(constructor%info, "data/mask_fullsky_n512_tqu.fits")

    constructor%nmaps    = info%nmaps
    constructor%ndet     = 4
    constructor%nhorn    = 1
    allocate(constructor%stokes(constructor%nmaps))
    allocate(constructor%w(constructor%nmaps,constructor%nhorn,constructor%ndet))
    allocate(constructor%label(constructor%ndet))
    constructor%stokes = [1,2,3]
    constructor%w      = 1.d0
    constructor%label  = ['27M','27S','28M','28S']

    ! Read the actual TOD
    call constructor%read_tod

    ! Initialize far sidelobe beams
    call open_hdf_file('data/LFI_instrument.h5', h5_file, 'r')
    nside_beam = 128
    call read_hdf(h5_file, trim(adjustl(constructor%label(1)))//'/'//'sllmax', lmax_beam)

    nmaps_beam = 1
    pol_beam   = .false.
    constructor%slinfo => comm_mapinfo(cpar%comm_chain, nside_beam, lmax_beam, &
            & nmaps_beam, pol_beam, .true.)
    allocate(constructor%slbeam(constructor%ndet), constructor%slconv(constructor%ndet))
    do i = 1, constructor%ndet
        constructor%slbeam(i)%p => comm_map(constructor%slinfo, h5_file, .true., .true., trim(constructor%label(i)))
    end do
    call close_hdf_file(h5_file)

  end function constructor

  !**************************************************
  !             Driver routine
  !**************************************************
  subroutine process_LFI_tod(self, map_in, delta, map_out, rms_out)
    implicit none
    class(comm_LFI_tod),                 intent(inout) :: self
    type(map_ptr),       dimension(:,:), intent(inout) :: map_in            ! (ndet,ndelta)
    real(dp),            dimension(:,:), intent(inout)    :: delta             ! (ndet,ndelta) BP corrections
    class(comm_map),                     intent(inout) :: map_out, rms_out  ! Combined output map and rms

    integer(i4b) :: i, j, k, ntod, ndet, nside, npix, nmaps, naccept, ntot
    integer(i4b) :: ierr, iter, main_iter, n_main_iter, ndelta
    real(dp)     :: t1, t2, t5, t6, chisq_threshold, delta_temp, cc, cp
    real(dp)     :: t_tot(14), accept_rate
    real(sp),     allocatable, dimension(:,:)     :: n_corr, s_sl, d_calib, s_sky, s_orb, mask,mask2, s_sb, s_sky_prop
    real(dp),     allocatable, dimension(:,:,:,:) :: map_sky
    real(dp),     allocatable, dimension(:)       :: A_abscal, b_abscal, chisq_prop, chisq_curr
    real(dp),     allocatable, dimension(:,:,:)   :: A_map
    real(dp),     allocatable, dimension(:,:)     :: b_map, map_temp
    real(dp),     allocatable, dimension(:,:)     :: procmask, procmask2
    integer(i4b), allocatable, dimension(:,:)     :: pix, psi, flag
    logical(lgt) :: dobp, accept
    logical(lgt), save :: first_call = .true.
    type(planck_rng) :: handle 
    character(5)                                  :: istr
    
    t_tot   = 0.d0
    call wall_time(t5)

!!$    write(*,*) delta(:,1)
!!$    write(*,*) delta(:,2)
!!$    call map_in(1,1)%p%writefits("map1_1.fits")
!!$    call map_in(2,1)%p%writefits("map2_1.fits")
!!$    call map_in(3,1)%p%writefits("map3_1.fits")
!!$    call map_in(4,1)%p%writefits("map4_1.fits")
!!$    call map_in(1,2)%p%writefits("map1_2.fits")
!!$    call map_in(2,2)%p%writefits("map2_2.fits")
!!$    call map_in(3,2)%p%writefits("map3_2.fits")
!!$    call map_in(4,2)%p%writefits("map4_2.fits")
!!$    call mpi_finalize(i)
!!$    stop

    call rand_init(handle, 243789)
    ! Set up full-sky map structures
    n_main_iter     = 1
    chisq_threshold = 1000000.d0 
    !this ^ should be 7.d0, is currently 2000 to debug sidelobes
    ndet            = self%ndet
    ndelta          = size(map_in,2)
    nside           = map_out%info%nside
    nmaps           = map_out%info%nmaps
    npix            = 12*nside**2
    allocate(A_map(nmaps,nmaps,0:npix-1), b_map(nmaps,0:npix-1))
    allocate(A_abscal(self%ndet), b_abscal(self%ndet))
    allocate(chisq_prop(self%ndet), chisq_curr(self%ndet))
    allocate(map_sky(0:npix-1,nmaps,0:ndet, ndelta),map_temp(0:npix-1,nmaps)) 
    allocate(procmask(0:npix-1,nmaps))
    allocate(procmask2(0:npix-1,nmaps))
    ! This step should be optimized -- currently takes 6 seconds..
    call wall_time(t1)
    do j = 1, ndelta
       do i = 1, self%ndet
          call map_in(i,j)%p%bcast_fullsky_map(map_sky(:,:,i,j))
       end do
    end do
    call self%procmask%bcast_fullsky_map(procmask)
    call self%procmask2%bcast_fullsky_map(procmask2)
    
    ! Trygve - Calculate mean map at zeroth index 
    map_sky(:,:,0,:) = map_sky(:,:,1,:)
    do i = 2, self%ndet
       map_sky(:,:,0,:) = map_sky(:,:,0,:) + map_sky(:,:,i,:) 
    end do
    map_sky(:,:,0,:) = map_sky(:,:,0,:)/self%ndet
    
!!$    do i = 0, self%ndet
!!$       write(*,*) i, sum(abs(map_sky(:,:,i)))
!!$    end do
!!$    stop

    map_sky = map_sky * 1.d-6 ! Gain in V/K
    call wall_time(t2)
    t_tot(9) = t2-t1


    ! Compute far sidelobe Conviqt structures
    call wall_time(t1)
    do i = 1, self%ndet
       call map_in(i,1)%p%remove_MDpoles() !remove mono and dipole components
       !call map_in(i,1)%p%writeFITS('nodp.fits')
       call map_in(i,1)%p%YtW()  ! Compute sky a_lms
       !TODO: make this work with shared arrays instead because that makes more
       !sense
       call self%slbeam(i)%p%bcast_fullsky_alms()
       if (self%myid == 0) write(*,*) 'precomputing sky', i
       self%slconv(i)%p => comm_conviqt(self%slbeam(i)%p, map_in(i,1)%p, 2)
       call self%slbeam(i)%p%distribute_alms()
    end do
    call wall_time(t2)
    t_tot(13) = t2-t1

    ! Toggle bp sampling on or off
    dobp = .true.

    A_abscal = 0.d0
    b_abscal = 0.d0

    ! Compute output map and rms
    do main_iter = 1, n_main_iter

       write(*,*) "main_iter", main_iter, A_abscal, b_abscal

       if (main_iter == n_main_iter) then
          ! Initialize main map arrays
          A_map = 0.d0
          b_map = 0.d0
          naccept = 0
          ntot    = 0

          ! Compute contribution to absolute calibration from current scan
          ! If you are running only one main iter I think this will be 0
          ! and thus will return NaN - Mathew
          call self%sample_absgain_from_orbital(A_abscal, b_abscal)

       else if (main_iter == n_main_iter-1) then
          ! Prepare for final absolute calibration from orbital
          A_abscal = 0.d0
          b_abscal = 0.d0
       end if
       
       if (main_iter == n_main_iter-1 .and. dobp) then
          ! trygve2
          ! Sample bp adjustments with mcmc sampler
          ! new variables: dobp, chisq_prop,_current

          do j = 1, self%ndet
             ! Summing chisq for current delta
             chisq_curr(j) = 0.d0
             chisq_prop(j) = 0.d0
             do i = 1, self%nscan ! Sum chisq for all scans
                chisq_curr(j) = chisq_curr(j) + self%scans(i)%d(j)%chisq_masked 
                chisq_prop(j) = chisq_prop(j) + self%scans(i)%d(j)%chisq_prop 
             end do
             
             cp = 0.d0
             cc = 0.d0
             call mpi_reduce(chisq_prop(j),cp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%info%comm, ierr)
             call mpi_reduce(chisq_curr(j),cc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%info%comm, ierr)
             chisq_prop(j) = cp
             chisq_curr(j) = cc
             
             ! Gather all proposed data
!!$             map_temp = 0.d0
!!$             delta_temp = 0.d0
!!$             call mpi_reduce(map_sky(:,:,j,2), map_temp, size(map_sky(:,:,j,2)), MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%info%comm, ierr)
!!$             call mpi_reduce(delta(j,2),delta_temp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%info%comm, ierr)

             if (self%myid == 0) then
                accept_rate = exp(-0.5d0*(cp-cc))  
                accept = (rand_uni(handle) < accept_rate)
                write(*,fmt='(a,i3,a,f12.1,a,f12.1,a,f12.1)') "det = ", j, ', c0 = ', cc, ', cp = ', cp, ', diff = ', cp-cc
                !write(*,*) map_sky(:,:,j,1) - map_temp ! accept proposal map 
                !map_sky(:,:,j,1) = map_temp ! accept proposal map 
                !delta(j,1) = delta_temp             ! accept proposal bp shift
             end if

             ! Broadcast new saved data
             call mpi_bcast(accept, 1,  MPI_LOGICAL, 0, self%info%comm, ierr)
             if (accept) then
                ! Set current to proposal
                map_sky(:,:,j,1) = map_sky(:,:,j,2)
                delta(j,1)       = delta(j,2)
             end if

!             call mpi_bcast(map_sky(:,:,j,1), size(map_sky(:,:,j,1)),  MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)

!             call mpi_bcast(delta(j,1), 1,  MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)
          end do
          
          !if (self%myid == 0 ) write(*,*) map_sky(:,:,1,1)
          !stop
          !call mpi_finalize()

       end if

       do i = 1, self%nscan
          
          ! Short-cuts to local variables
          ntod = self%scans(i)%ntod
          ndet = self%ndet
          
          ! Set up local data structure for current scan
          allocate(d_calib(ntod, ndet))            ! Calibrated and cleaned TOD in uKcmb
          allocate(n_corr(ntod, ndet))             ! Correlated noise in V
          allocate(s_sl(ntod, ndet))               ! Sidelobe in uKcm
          allocate(s_sky(ntod, ndet))              ! Stationary sky signal in uKcmb
          allocate(s_sky_prop(ntod, ndet))              ! Stationary sky signal in uKcmb
          allocate(s_sb(ntod, ndet))               ! Signal minus mean
          allocate(s_orb(ntod, ndet))              ! Orbital dipole in uKcmb
          allocate(mask(ntod, ndet))               ! Processing mask in time
          allocate(mask2(ntod, ndet))               ! Processing mask in time
          allocate(pix(ntod, ndet))                ! Decompressed pointing
          allocate(psi(ntod, ndet))                ! Decompressed discretized polarization angle
          allocate(flag(ntod, ndet))               ! Decompressed flags
          
          ! Initializing arrays to zero
          n_corr = 0.d0
          s_sl   = 0.d0
          s_sky  = 0.d0
          s_sky_prop  = 0.d0
          s_orb  = 0.d0
          
          ! --------------------
          ! Analyze current scan
          ! --------------------
         
          ! Decompress pointing, psi and flags for current scan
          call wall_time(t1)
          do j = 1, ndet
             call self%decompress_pointing_and_flags(i, j, pix(:,j), psi(:,j), flag(:,j))
          end do
          call wall_time(t2)       
          t_tot(11) = t_tot(11) + t2-t1
          
          ! Construct sky signal template -- Maksym -- this week - Trygve added mean sky
          call wall_time(t1)
          do j = 1, ndet            
             call self%project_sky(map_sky(:,:,j,:), pix(:,j), psi(:,j), procmask,procmask2, i, j, s_sky(:,j), mask(:,j),mask2(:,j), map_sky(:,:,0,:), s_sb(:,j), s_sky_prop(:,j), dobp)  ! scan_id, det,  s_sky(:, j))
          end do
          call wall_time(t2)
          t_tot(1) = t_tot(1) + t2-t1
          !if (self%myid == 0) write(*,*) 'Project    = ', t2-t1
          
          ! Construct orbital dipole template -- Kristian -- this week-ish
          call wall_time(t1)
          call self%compute_orbital_dipole(i, pix, s_orb)
          call wall_time(t2)
          t_tot(2) = t_tot(2) + t2-t1
          !if (self%myid == 0) write(*,*) 'Orb dipole = ', t2-t1
          
          ! Construct sidelobe template -- Mathew -- long term
          call wall_time(t1)
          do j = 1, ndet
             call self%construct_sl_template(self%slconv(j)%p, i, nside, pix(:,j), psi(:,j), s_sl(:,j))
          end do
          call wall_time(t2)
          t_tot(12) = t_tot(12) + t2-t1
          
          do iter = 1, 1
             
             ! Fit correlated noise -- Haavard -- this week-ish
             call wall_time(t1)
             call self%sample_n_corr(i, mask, s_sky, s_sl, s_orb, n_corr)
             call wall_time(t2)
             t_tot(3) = t_tot(3) + t2-t1
             !if (self%myid == 0) write(*,*) 'ncorr      = ', t2-t1
             
             ! Fit gain 
             call wall_time(t1)
             if (main_iter < n_main_iter) then
                ! Fit both relative and absolute gain
                do j = 1, ndet
                   call sample_gain(self, j, i, n_corr(:, j), mask(:,j), s_sky(:, j), s_sl(:, j), s_orb(:, j))
                end do
             end if
             call wall_time(t2)
             t_tot(4) = t_tot(4) + t2-t1
             !if (self%myid == 0) write(*,*) 'gain       = ', t2-t1
             
             ! Compute noise spectrum
             call wall_time(t1)
             do j = 1, ndet
                call self%sample_noise_psd(i, j, mask(:,j), s_sky(:,j), s_sl(:,j), s_orb(:,j), n_corr(:,j))
             end do
             call wall_time(t2)
             t_tot(6) = t_tot(6) + t2-t1
             !if (self%myid == 0) write(*,*) 'noise      = ', t2-t1
             
          end do

          if (main_iter == n_main_iter-1 .or. dobp) then
             ! Compute contribution to absolute calibration from current scan
             call wall_time(t1)
             do j = 1, ndet
                call self%compute_chisq(i, j, mask(:,j), mask2(:,j), s_sky(:,j), s_sl(:,j), s_orb(:,j), n_corr(:,j), s_sky_prop(:,j), dobp)
                
               if (abs(self%scans(i)%d(j)%chisq) > chisq_threshold .or. &
                     & isNaN(self%scans(i)%d(j)%chisq)) then
                 !write(*,*) "testing chisq", self%scans(i)%d(j)%chisq
                 cycle
                end if
                call self%accumulate_absgain_from_orbital(i, j, mask(:,j), s_sky(:,j), &
                     & s_sl(:,j), s_orb(:,j), n_corr(:,j), A_abscal(j), b_abscal(j))
             end do
             call wall_time(t2)
             t_tot(14) = t_tot(14) + t2-t1
          end if


          if (main_iter == n_main_iter) then
             ! Compute chisquare
             call wall_time(t1)
             do j = 1, ndet
                call self%compute_chisq(i, j, mask(:,j), mask2(:,j), s_sky(:,j), s_sl(:,j), s_orb(:,j), n_corr(:,j), s_sky_prop(:,j), dobp)
             end do
             call wall_time(t2)
             t_tot(7) = t_tot(7) + t2-t1
             !if (self%myid == 0) write(*,*) 'chisq      = ', t2-t1
             
             ! Compute clean and calibrated TOD -- Mathew -- this week
             call wall_time(t1)
             call self%compute_cleaned_tod(ntod, i, s_orb, s_sl, n_corr, d_calib,s_sb)
             call wall_time(t2)
             t_tot(5) = t_tot(5) + t2-t1
             !if (self%myid == 0) write(*,*) 'clean      = ', t2-t1

!!$             open(58,file='poster_tod.dat',recl=1024)
!!$             do j = 1, ntod
!!$                write(58,*) j, self%scans(i)%d(1)%tod(j), n_corr(j,1), s_orb(j,1), s_sky(j,1), s_sb(j,i), d_calib(j,1)
!!$             end do
!!$             close(58)
!!$             call mpi_finalize(ierr)
!!$             stop

             ! Compute binned map from cleaned TOD -- Marie -- this week
             call wall_time(t1)
             do j = 1, ndet
                ntot= ntot + 1
                if (first_call) then
                   if (abs(self%scans(i)%d(j)%chisq) > chisq_threshold .or. &
                        & self%scans(i)%d(j)%chisq /= self%scans(i)%d(j)%chisq) then
                      self%scans(i)%d(j)%accept = .false.
                   else
                      self%scans(i)%d(j)%accept = .true.                   
                   end if
                end if
                if (self%scans(i)%d(j)%accept) then
                   naccept = naccept + 1
                   call self%compute_binned_map(d_calib, pix(:,j), psi(:,j), flag(:,j), A_map, b_map, i, j)
                end if
             end do
             
             do j = 1, ndet
                if (abs(self%scans(i)%d(j)%chisq) > chisq_threshold .or. &
                     & self%scans(i)%d(j)%chisq /= self%scans(i)%d(j)%chisq) &
                     write(*,fmt='(a,i8,i5,a,f12.1)') 'Reject scan, det = ', &
               & self%scanid(i), j, ', chisq = ', self%scans(i)%d(j)%chisq
             end do

          end if
          !Write out sl template somehow
          !write(istr, '(I5.5)') i
          !open(134, file= 'sl'//trim(istr)//'.txt', status='replace',recl=10000)
          !do j = 1, size(s_sl(:,0))
          !  write(134,*) s_sl(j,1), s_sky(j, 1), d_calib(j, 1), s_orb(j, 1), n_corr(j, 1), s_sb(j, 1)
          !end do
          !close(134)

          ! Clean up
          deallocate(n_corr, s_sl, s_sky, s_orb, d_calib, mask,mask2, pix, psi, flag, s_sb, s_sky_prop)!, s_sb_prop)


       end do


       call wall_time(t2)
       t_tot(8) = t_tot(8) + t2-t1
       !if (self%myid == 0) write(*,*) 'bin        = ', t2-t1
       
    end do

    ! Parameter to check if this is first time routine has been 
    first_call = .false.

    ! Solve combined map, summed over all pixels -- Marie 
    call wall_time(t1)
    call finalize_binned_map(A_map, b_map, map_out, rms_out)
    call wall_time(t2)
    t_tot(10) = t2-t1

    call mpi_reduce(ntot, i, 1, MPI_INTEGER, MPI_SUM, 0, self%info%comm, ierr)
    ntot = i
    call mpi_reduce(naccept, i, 1, MPI_INTEGER, MPI_SUM, 0, self%info%comm, ierr)
    naccept = i

    call wall_time(t6)
    if (self%myid == 0) then
       write(*,*) '  Time dist sky   = ', t_tot(9)
       write(*,*) '  Time sl precomp = ', t_tot(13)
       write(*,*) '  Time decompress = ', t_tot(11)
       write(*,*) '  Time project    = ', t_tot(1)
       write(*,*) '  Time orbital    = ', t_tot(2)
       write(*,*) '  Time sl interp  = ', t_tot(12)
       write(*,*) '  Time ncorr      = ', t_tot(3)
       write(*,*) '  Time gain       = ', t_tot(4)
       write(*,*) '  Time absgain    = ', t_tot(14)
       write(*,*) '  Time clean      = ', t_tot(5)
       write(*,*) '  Time noise      = ', t_tot(6)
       write(*,*) '  Time chisq      = ', t_tot(7)
       write(*,*) '  Time bin        = ', t_tot(8)
       write(*,*) '  Time final      = ', t_tot(10)
       write(*,*) '  Time total      = ', t6-t5, ', accept rate = ', real(naccept,sp) / ntot
       write(*,*) naccept, ntot
    end if

    call map_out%writeFITS('map.fits')
    call rms_out%writeFITS('rms.fits')

    !stop
    ! Clean up temporary arrays
    deallocate(map_sky, A_map, b_map, procmask, procmask2)
    deallocate(A_abscal, b_abscal)

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

!    write(*,*) self%scans(scan)%d(det)%pix(1:5)
    call huffman_decode(self%scans(scan)%hkey, self%scans(scan)%d(det)%pix,  pix)
    call huffman_decode(self%scans(scan)%hkey, self%scans(scan)%d(det)%psi,  psi)
    call huffman_decode(self%scans(scan)%hkey, self%scans(scan)%d(det)%flag, flag)
!    write(*,*) self%scans(scan)%ntod
!    write(*,*) pix(1:5)
!    write(*,*) pix(self%scans(scan)%ntod-4:self%scans(scan)%ntod)
    do j = 2, self%scans(scan)%ntod
       pix(j)  = pix(j-1)  + pix(j)
       psi(j)  = psi(j-1)  + psi(j)
       flag(j) = flag(j-1) + flag(j)
    end do
    psi = min(psi,4095)
!    write(*,*) pix(1:5)
    !write(*,*) psi(1:5)
    !write(*,*) flag(1:5)


  end subroutine decompress_pointing_and_flags

  ! Sky signal template
  subroutine project_sky(self, map, pix, psi, pmask,pmask2, scan_id, det, s_sky, tmask,tmask2, map_mean, s_sb, s_sky_prop, dobp)
    implicit none
    class(comm_LFI_tod),                  intent(in)  :: self
    real(dp),            dimension(0:,1:),  intent(in)  :: pmask,pmask2
    real(dp),            dimension(0:,1:,:),  intent(in)  :: map, map_mean !Added dimension trygve
    integer(i4b),        dimension(:),    intent(in)  :: pix, psi
    integer(i4b),                         intent(in)  :: scan_id, det
    real(sp),            dimension(:),    intent(out) :: s_sky, tmask,tmask2, s_sb, s_sky_prop !, s_sb_prop

    integer(i4b)                                      :: i!, pix
    logical(lgt), intent(in) :: dobp
    !real(dp)                                          :: psi

    ! s = T + Q * cos(2 * psi) + U * sin(2 * psi)
    ! T - temperature; Q, U - Stoke's parameters
    do i = 1, self%scans(scan_id)%ntod
       ! note that psi(i) is an index now; must be converted to a real number based on lookup table
!       s_sky(i) = map(pix(i), 1) + map(pix(i), 2) * cos(2.d0 * psi(i)) + map(pix(i), 3) * sin(2.d0 * psi(i))
       s_sky(i) = map(pix(i), 1,1) + map(pix(i), 2,1) * self%cos2psi(psi(i)) + map(pix(i), 3,1) * self%sin2psi(psi(i))
       s_sb(i)  = s_sky(i) - (map_mean(pix(i), 1,1) + map_mean(pix(i), 2,1) * self%cos2psi(psi(i)) + map_mean(pix(i), 3,1)*self%sin2psi(psi(i)))

       ! Calculate proposed sky if asked to
       if ( dobp ) then
          s_sky_prop(i) = map(pix(i), 1, 2) + map(pix(i), 2, 2) * self%cos2psi(psi(i)) + map(pix(i), 3, 2) * self%sin2psi(psi(i))
       end if

       if (s_sky(i) /= s_sky(i)) then
          write(*,*) 'Nan in s_sky', scan_id, i, map(pix(i),1,:), psi(i), self%cos2psi(psi(i)), self%sin2psi(psi(i))
          stop
       end if
       if (any(pmask(pix(i),:) < 0.5d0)) then
          tmask(i) = 0.
       else
          tmask(i) = 1.
       end if
       if (any(pmask2(pix(i),:) < 0.5d0)) then
          tmask2(i) = 0.
       else
          tmask2(i) = 1.
       end if
    end do

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
       do j=1,self%scans(ind)%ntod !length of the tod
          if(pix(j, i) < 0 .or. pix(j,i) > 12*(self%scans(ind)%d(i)%nside**2)) then
            !write(*,*) pix(j, i), self%scans(ind)%d(i)%nside
            cycle
        end if
          call pix2vec_ring(self%scans(ind)%d(i)%nside, pix(j,i), &
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
    integer(i4b) :: i, j, k, l, n, nomp, ntod, ndet, err, omp_get_max_threads, meanrange, start, last
    integer*8    :: plan_fwd, plan_back
    real(sp)     :: sigma_0, alpha, nu_knee, nu, samprate, gain, mean
    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv
    real(sp),     allocatable, dimension(:) :: d_prime, diff
    
    ntod = self%scans(scan)%ntod
    ndet = self%ndet
    nomp = omp_get_max_threads()
    
!    do i = 1, ndet
!       gain = 1.d-6 * self%scans(scan)%d(i)%gain  ! Gain in V / muK
!       n_corr(:,i) = self%scans(scan)%d(i)%tod(:) - S_sky(:,i) * gain 
!    end do
!    return
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
       gain = self%scans(scan)%d(i)%gain  ! Gain in V / K
       d_prime(:) = self%scans(scan)%d(i)%tod(:) - S_sl(:,i) - (S_sky(:,i) + S_orb(:,i)) * gain
      
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
       samprate = self%scans(i)%d(i)%samprate
       sigma_0 = self%scans(scan)%d(i)%sigma0
       alpha = self%scans(scan)%d(i)%alpha
       nu_knee = self%scans(scan)%d(i)%fknee
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

    d_only_wn     = self%scans(scan_id)%d(det)%tod - n_corr
    gain_template = s_sky + s_orb + s_sl
    ata           = sum(mask*gain_template**2)
    curr_gain     = sum(mask * d_only_wn * gain_template) / ata            
    curr_sigma    = self%scans(scan_id)%d(det)%sigma0 / sqrt(ata)  
    self%scans(scan_id)%d(det)%gain       = curr_gain
    self%scans(scan_id)%d(det)%gain_sigma = curr_sigma

    if(isNaN(self%scans(scan_id)%d(det)%gain)) then
      write(*,*) "Gain estimate", sum(s_sky), sum(s_orb), sum(s_sl), ata, curr_gain, curr_sigma
    end if

    deallocate(d_only_wn,gain_template)

  end subroutine sample_gain

  subroutine accumulate_absgain_from_orbital(self, scan, det, mask, s_sky, s_sl, s_orb, n_corr, A_abs, b_abs)
    implicit none
    class(comm_LFI_tod),             intent(in)     :: self
    integer(i4b),                    intent(in)     :: scan, det
    real(sp),          dimension(:), intent(in)     :: mask, s_sky, s_sl, s_orb, n_corr
    real(dp),                        intent(inout)  :: A_abs, b_abs

    integer(i4b) :: i
    real(dp)     :: data          ! Cleaned data in units of K
    real(dp)     :: inv_sigmasq  ! Inverse variance in units of K**-2

    inv_sigmasq = (self%scans(scan)%d(det)%gain/self%scans(scan)%d(det)%sigma0)**2
    do i = 1, self%scans(scan)%ntod
       data = (self%scans(scan)%d(det)%tod(i)-n_corr(i))/self%scans(scan)%d(det)%gain - s_sky(i) - s_sl(i)
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
    call mpi_reduce(A_abs, A, self%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%info%comm, ierr)
    call mpi_reduce(b_abs, b, self%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%info%comm, ierr)

    ! Compute gain update and distribute to all cores
    if (self%myid == 0) then
       sigma = 1.d0/sqrt(A)
       dgain = b/A
       write(*,*) "dgain", b, A
       if (.false.) then
          ! Add fluctuation term if requested
       end if
       do j = 1, self%ndet
          write(*,fmt='(a,i5,a,3f8.3)') 'Orb gain -- d = ', j, ', dgain = ', &
               & 100*(dgain(j)-1), 100*sigma(j), (dgain(j)-1)/sigma(j)
       end do
    end if
    call mpi_bcast(dgain, self%ndet,  MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)
    call mpi_bcast(sigma, self%ndet,  MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)

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
  subroutine construct_sl_template(self, slconv, scan_id, nside, pix, psi, s_sl)
    implicit none
    class(comm_LFI_tod),                 intent(in)    :: self
    class(comm_conviqt),                 intent(in)    :: slconv
    integer(i4b),                        intent(in)    :: scan_id, nside
    integer(i4b),        dimension(:),   intent(in)    :: pix, psi
    real(sp),            dimension(:),   intent(out)   :: s_sl
    
    integer(i4b) :: i, j, p
    real(dp)     :: psi_, theta, phi

    do j=1, size(pix)
       !if(pix(j) > 12*(nside**2) .or. pix(j) < 0) then
       !  write(*,*) pix(j), nside
       !end if
       call pix2ang_ring(nside, pix(j), theta, phi)
       psi_    = psi(j) ! HKE: Should this be defined without the detector angle
       s_sl(j) = slconv%interp(theta, phi, psi_)  ! Commented out until validated
       if(abs(s_sl(j)) > 10.d0) then
         !write(*,*) s_sl(j), pix, theta, phi, psi_
       end if
       !s_sl(j) = 0.d0
    end do

  end subroutine construct_sl_template

  !compute the cleaned TOD from the computed TOD components
  subroutine compute_cleaned_tod(self, ntod, scan_num, s_orb, s_sl, n_corr, d_calib,s_sb)
    implicit none
    class(comm_LFI_tod),               intent(in)    :: self
    integer(i4b),                      intent(in)    :: ntod, scan_num
    real(sp),          dimension(:,:), intent(in)    :: n_corr, s_sl, s_orb,s_sb
    real(sp),          dimension(:,:), intent(out)   :: d_calib

    integer(i4b) :: i

    !cleaned calibrated data = (rawTOD - corrNoise)/gain - orbitalDipole - sideLobes
    do i=1, self%ndet
      d_calib(:,i) = (self%scans(scan_num)%d(i)%tod - n_corr(:,i))/ self%scans(scan_num)%d(i)%gain - s_orb(:,i) - s_sb(:,i) -s_sl(:,i)
      !d_calib(:,i) = s_sl(:,i)
    end do
  end subroutine compute_cleaned_tod


  ! Compute chisquare
  subroutine compute_chisq(self, scan, det, mask,mask2, s_sky, s_sl, s_orb, n_corr, s_sky_prop, dobp)
    implicit none
    class(comm_LFI_tod),             intent(inout)  :: self
    integer(i4b),                    intent(in)     :: scan, det
    real(sp),          dimension(:), intent(in)     :: mask,mask2, s_sky, s_sl, s_orb, n_corr, s_sky_prop

    real(dp) :: chisq, chisq_prop, chisq_masked
    integer(i4b) :: i, n
    logical(lgt), intent(in) :: dobp

    chisq = 0.d0
    chisq_prop = 0.d0
    chisq_masked = 0.d0
    n     = 0
    do i = 1, self%scans(scan)%ntod
       if (mask(i) > 0.5) chisq = chisq + (self%scans(scan)%d(det)%tod(i) - &
            & (self%scans(scan)%d(det)%gain * (s_sky(i) + s_sl(i) + s_orb(i)) + &
            & n_corr(i)))**2/self%scans(scan)%d(det)%sigma0**2
       
       
       if ( dobp ) then
       
          if (mask2(i) > 0.5) chisq_masked = chisq_masked + (self%scans(scan)%d(det)%tod(i) - &
               & (self%scans(scan)%d(det)%gain * (s_sky(i) + s_sl(i) + s_orb(i)) + &
               & n_corr(i)))**2/self%scans(scan)%d(det)%sigma0**2
          
          if (mask2(i) > 0.5) chisq_prop = chisq_prop + (self%scans(scan)%d(det)%tod(i) - &
               & (self%scans(scan)%d(det)%gain * (s_sky_prop(i) + s_sl(i) + s_orb(i)) + &
               & n_corr(i)))**2/self%scans(scan)%d(det)%sigma0**2
       end if
       n = n+1
       if (.false. .and. self%myid == 0) write(*,*) i, chisq/n, (self%scans(scan)%d(det)%tod(i) - &
         & (self%scans(scan)%d(det)%gain * (s_sky(i) + s_sl(i) + s_orb(i)) + &
         & n_corr(i)))**2/self%scans(scan)%d(det)%sigma0**2
    end do
    !close(58)
    self%scans(scan)%d(det)%chisq = (chisq - n) / sqrt(2.d0*n)
    self%scans(scan)%d(det)%chisq_prop = chisq_prop
    self%scans(scan)%d(det)%chisq_masked = chisq_masked

    if(self%scans(scan)%d(det)%chisq > 2000.d0 .or. isNaN(self%scans(scan)%d(det)%chisq)) then
      !write(*,*) "chisq", scan, det, sum(mask), sum(s_sky), sum(s_sl), sum(s_orb), sum(n_corr)
    end if

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
  subroutine compute_binned_map(self, data, pix, psi, flag, A, b, scan, det)
    implicit none
    class(comm_LFI_tod),                      intent(in)    :: self
    integer(i4b),                             intent(in)    :: scan, det
    real(sp),            dimension(:,:),      intent(in)    :: data
    integer(i4b),        dimension(:),        intent(in)    :: pix, psi, flag
    real(dp),            dimension(1:,1:,0:), intent(inout) :: A
    real(dp),            dimension(1:,0:),    intent(inout) :: b

    integer(i4b) :: i, j, t, pix_
    real(dp)     :: psi_, inv_sigmasq

    inv_sigmasq = (self%scans(scan)%d(det)%gain/self%scans(scan)%d(det)%sigma0)**2
    do t = 1, self%scans(scan)%ntod

       if (iand(flag(t),6111248) .ne. 0) cycle

       pix_    = pix(t)
       psi_    = psi(t)
       
       A(1,1,pix_) = A(1,1,pix_) + 1.d0                                          * inv_sigmasq
       A(1,2,pix_) = A(1,2,pix_) + self%cos2psi(psi_)                        * inv_sigmasq
       A(1,3,pix_) = A(1,3,pix_) + self%sin2psi(psi_)                        * inv_sigmasq
       A(2,2,pix_) = A(2,2,pix_) + self%cos2psi(psi_)**2                     * inv_sigmasq
       A(2,3,pix_) = A(2,3,pix_) + self%cos2psi(psi_)*self%sin2psi(psi_) * inv_sigmasq
       A(3,3,pix_) = A(3,3,pix_) + self%sin2psi(psi_)**2                     * inv_sigmasq

       b(1,pix_) = b(1,pix_) + data(t,det)                      * inv_sigmasq
       b(2,pix_) = b(2,pix_) + data(t,det) * self%cos2psi(psi_) * inv_sigmasq
       b(3,pix_) = b(3,pix_) + data(t,det) * self%sin2psi(psi_) * inv_sigmasq

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
          rms%map(i,j) = sqrt(A_inv(j,j)) * 1.d6 ! uK
       end do

       ! map
       map%map(i,:) = matmul(A_inv,b_tot(:,i)) * 1.d6 ! uK

    end do

    deallocate(A_inv,A_tot,b_tot)

  end subroutine finalize_binned_map

end module comm_tod_LFI_mod
