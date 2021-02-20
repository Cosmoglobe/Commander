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
module comm_tod_LFI_mod
  use comm_tod_mod
  use comm_param_mod
  use comm_map_mod
  use comm_conviqt_mod
  use pix_tools
  use healpix_types
  use comm_huffman_mod
  use comm_hdf_mod
  use comm_fft_mod
  use comm_shared_arr_mod
  use spline_1D_mod
  use comm_4D_map_mod
  use comm_zodi_mod
  use comm_tod_mapmaking_mod
  use comm_tod_pointing_mod
  use comm_tod_gain_mod
  use comm_tod_bandpass_mod
  use comm_tod_orbdipole_mod
  use comm_utils
  implicit none

  private
  public comm_LFI_tod

  integer(i4b), parameter :: N_test      = 19
  integer(i4b), parameter :: samp_N      = 1
  integer(i4b), parameter :: prep_G      = 15
  integer(i4b), parameter :: samp_G      = 2
  integer(i4b), parameter :: prep_acal   = 3
  integer(i4b), parameter :: samp_acal   = 4
  integer(i4b), parameter :: prep_rcal   = 18
  integer(i4b), parameter :: samp_rcal   = 19
  integer(i4b), parameter :: prep_relbp  = 5
  integer(i4b), parameter :: prep_absbp  = 16
  integer(i4b), parameter :: samp_bp     = 11
  integer(i4b), parameter :: samp_sl     = 6
  integer(i4b), parameter :: samp_N_par  = 7
  integer(i4b), parameter :: sel_data    = 8
  integer(i4b), parameter :: bin_map     = 9
  integer(i4b), parameter :: calc_chisq  = 10
  integer(i4b), parameter :: output_slist = 12
  integer(i4b), parameter :: samp_mono   = 13
  integer(i4b), parameter :: sub_sl      = 14
  integer(i4b), parameter :: sub_zodi    = 17
  logical(lgt), dimension(N_test) :: do_oper


  type, extends(comm_tod) :: comm_LFI_tod
    class(orbdipole_pointer), allocatable :: orb_dp !orbital dipole calculator
   contains
     procedure     :: process_tod        => process_LFI_tod
     procedure     :: simulate_LFI_tod
  end type comm_LFI_tod

  interface comm_LFI_tod
     procedure constructor
  end interface comm_LFI_tod

contains

  !**************************************************
  !             Constructor
  !**************************************************
  function constructor(cpar, id_abs, info, tod_type)
    implicit none
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id_abs
    class(comm_mapinfo),     target     :: info
    character(len=128),      intent(in) :: tod_type
    class(comm_LFI_tod),     pointer    :: constructor

    integer(i4b) :: i, j, k, nside_beam, lmax_beam, nmaps_beam, ndelta, np_vec, ierr, offset
    real(dp)     :: f_fill, f_fill_lim(3), theta, phi, pixVal
    character(len=512) :: datadir
    logical(lgt) :: pol_beam
    type(hdf_file) :: h5_file
    integer(i4b), allocatable, dimension(:) :: pix
    real(dp), dimension(3) :: v
    character(len=10) :: fieldname
    class(tod_pointer), pointer :: selfptr

    ! Set up LFI specific parameters
    allocate(constructor)
    constructor%tod_type      = tod_type
    constructor%myid          = cpar%myid_chain
    constructor%comm          = cpar%comm_chain
    constructor%numprocs      = cpar%numprocs_chain
    constructor%myid_shared   = cpar%myid_shared
    constructor%comm_shared   = cpar%comm_shared
    constructor%myid_inter    = cpar%myid_inter
    constructor%comm_inter    = cpar%comm_inter
    constructor%info          => info
    constructor%init_from_HDF = cpar%ds_tod_initHDF(id_abs)
    constructor%freq          = cpar%ds_label(id_abs)
    constructor%operation     = cpar%operation
    constructor%outdir        = cpar%outdir
    constructor%first_call    = .true.
    constructor%first_scan    = cpar%ds_tod_scanrange(id_abs,1)
    constructor%last_scan     = cpar%ds_tod_scanrange(id_abs,2)
    constructor%flag0         = cpar%ds_tod_flag(id_abs)
    constructor%orb_abscal    = cpar%ds_tod_orb_abscal(id_abs)
    constructor%nscan_tot     = cpar%ds_tod_tot_numscan(id_abs)
    constructor%output_4D_map = cpar%output_4D_map_nth_iter
    constructor%output_aux_maps = cpar%output_aux_maps
    constructor%subtract_zodi = cpar%include_TOD_zodi
    constructor%central_freq  = cpar%ds_nu_c(id_abs)
    constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
    constructor%halfring_split = cpar%ds_tod_halfring(id_abs)
    constructor%nside_param   = cpar%ds_nside(id_abs)

    !----------------------------------------------------------------------------------
    ! Simulation Routine
    !----------------------------------------------------------------------------------
    constructor%sims_output_dir = cpar%sims_output_dir
    constructor%enable_tod_simulations = cpar%enable_tod_simulations
    !----------------------------------------------------------------------------------

    call mpi_comm_size(cpar%comm_shared, constructor%numprocs_shared, ierr)

    if (constructor%first_scan > constructor%last_scan) then
       write(*,*) 'Error: First scan larger than last scan for ', trim(constructor%freq)
       call mpi_finalize(ierr)
       stop
    end if

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
    allocate(constructor%partner(constructor%ndet))
    allocate(constructor%horn_id(constructor%ndet))
    constructor%stokes = [1,2,3]
    constructor%w      = 1.d0
    call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label)
    do i = 1, constructor%ndet
       if (mod(i,2) == 1) then
          constructor%partner(i) = i+1
       else
          constructor%partner(i) = i-1
       end if
       constructor%horn_id(i) = (i+1)/2
    end do


    if (trim(cpar%ds_bpmodel(id_abs)) == 'additive_shift') then
       ndelta = 1
    else if (trim(cpar%ds_bpmodel(id_abs)) == 'powlaw_tilt') then
       ndelta = 1
    else
       write(*,*) 'Unknown bandpass model'
       stop
    end if
    allocate(constructor%bp_delta(0:constructor%ndet,ndelta))
    constructor%bp_delta = 0.d0

    ! Read the actual TOD
    call constructor%read_tod(constructor%label)

    ! Initialize bandpass mean and proposal matrix
    call constructor%initialize_bp_covar(trim(datadir)//cpar%ds_tod_bp_init(id_abs))

    ! Initialize beams
    allocate(constructor%fwhm(constructor%ndet))
    allocate(constructor%elip(constructor%ndet))
    allocate(constructor%psi_ell(constructor%ndet))
    allocate(constructor%mb_eff(constructor%ndet))
    allocate(constructor%nu_c(constructor%ndet))

    call open_hdf_file(constructor%instfile, h5_file, 'r')
    nside_beam = 512
    nmaps_beam = 3
    pol_beam   = .true.
    call read_hdf(h5_file, trim(adjustl(constructor%label(1)))//'/'//'sllmax', lmax_beam)
    constructor%slinfo => comm_mapinfo(cpar%comm_chain, nside_beam, lmax_beam, &
         & nmaps_beam, pol_beam)

    allocate(constructor%slbeam(constructor%ndet), constructor%slconv(constructor%ndet))
    allocate(constructor%mbeam(constructor%ndet))

    do i = 1, constructor%ndet
       constructor%slbeam(i)%p => comm_map(constructor%slinfo, h5_file, .true., "sl", &
            & trim(constructor%label(i)))
       constructor%mbeam(i)%p => comm_map(constructor%slinfo, h5_file, .true., "beam", trim(constructor%label(i)))
       call constructor%mbeam(i)%p%Y()
    end do

    allocate(constructor%orb_dp)
    constructor%orb_dp%p => comm_orbdipole(constructor, constructor%mbeam)

    do i = 1, constructor%ndet
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'fwhm', constructor%fwhm(i))
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'elip', constructor%elip(i))
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'psi_ell', constructor%psi_ell(i))
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'mbeam_eff', constructor%mb_eff(i))
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'centFreq', constructor%nu_c(i))
    end do

    constructor%mb_eff = 1.d0
    constructor%mb_eff = constructor%mb_eff / mean(constructor%mb_eff)
    constructor%nu_c   = constructor%nu_c * 1d9

    call close_hdf_file(h5_file)

    ! Lastly, create a vector pointing table for fast look-up for orbital dipole
    np_vec = 12*constructor%nside**2 !npix
    allocate(constructor%pix2vec(3,0:np_vec-1))
    do i = 0, np_vec-1
       call pix2vec_ring(constructor%nside, i, constructor%pix2vec(:,i))
    end do

    ! Construct observed pixel array
    allocate(constructor%pix2ind(0:12*constructor%nside**2-1))
    constructor%pix2ind = -1
    do i = 1, constructor%nscan
       allocate(pix(constructor%scans(i)%ntod))
       offset = 0
       if(constructor%halfring_split == 2) then
         offset = constructor%scans(i)%ntod
       end if
       do j = 1, constructor%ndet
          call huffman_decode(constructor%scans(i)%hkey, &
               constructor%scans(i)%d(j)%pix(1)%p, pix)
          constructor%pix2ind(pix(1)) = 1
          do k = 2, constructor%scans(i)%ntod
             pix(k)  = pix(k-1)  + pix(k)
             constructor%pix2ind(pix(k)) = 1
          end do
       end do
       deallocate(pix)
    end do
    constructor%nobs = count(constructor%pix2ind == 1)
    allocate(constructor%ind2pix(constructor%nobs))
    allocate(constructor%ind2sl(constructor%nobs))
    allocate(constructor%ind2ang(2,constructor%nobs))
    j = 1
    do i = 0, 12*constructor%nside**2-1
       if (constructor%pix2ind(i) == 1) then
          constructor%ind2pix(j) = i
          constructor%pix2ind(i) = j
          call pix2ang_ring(constructor%nside, i, theta, phi)
          call ang2pix_ring(nside_beam, theta, phi, constructor%ind2sl(j))
          constructor%ind2ang(:,j) = [theta,phi]
          j = j+1
       end if
    end do
    f_fill = constructor%nobs/(12.*constructor%nside**2)
    call mpi_reduce(f_fill, f_fill_lim(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, constructor%info%comm, ierr)
    call mpi_reduce(f_fill, f_fill_lim(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, constructor%info%comm, ierr)
    call mpi_reduce(f_fill, f_fill_lim(3), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, constructor%info%comm, ierr)
    if (constructor%myid == 0) then
       write(*,*) '  Min/mean/max TOD-map f_sky = ', real(100*f_fill_lim(1),sp), real(100*f_fill_lim(3)/constructor%info%nprocs,sp), real(100*f_fill_lim(2),sp)
    end if

  end function constructor

  !**************************************************
  !             Driver routine
  !**************************************************
  subroutine process_LFI_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
    implicit none
    class(comm_LFI_tod),                      intent(inout) :: self
    character(len=*),                         intent(in)    :: chaindir
    integer(i4b),                             intent(in)    :: chain, iter
    type(planck_rng),                         intent(inout) :: handle
    type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
    real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
    class(comm_map),                          intent(inout) :: map_out      ! Combined output map
    class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms

    integer(i4b) :: i, j, k, l, start_chunk, end_chunk, chunk_size, ntod, ndet
    integer(i4b) :: nside, npix, nmaps, naccept, ntot, ext(2), nscan_tot, nhorn
    integer(i4b) :: ierr, main_iter, n_main_iter, ndelta, ncol, n_A, nout=1
    real(dp)     :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, chisq_threshold
    real(dp)     :: t_tot(23)
    real(sp)     :: inv_gain
    real(sp),     allocatable, dimension(:,:)     :: n_corr, s_sl, s_sky, s_orb, mask,mask2, s_bp
    real(sp),     allocatable, dimension(:,:)     :: s_mono, s_buf, s_tot, s_zodi
    real(sp),     allocatable, dimension(:,:)     :: s_invN, s_lowres
    real(sp),     allocatable, dimension(:,:,:)   :: s_sky_prop, s_bp_prop
    real(sp),     allocatable, dimension(:,:,:)   :: d_calib
    real(dp),     allocatable, dimension(:)       :: A_abscal, b_abscal
    real(dp),     allocatable, dimension(:,:)     :: chisq_S, m_buf
    real(dp),     allocatable, dimension(:,:)     :: A_map, dipole_mod
    real(dp),     allocatable, dimension(:,:,:)   :: b_map, b_mono, sys_mono
    integer(i4b), allocatable, dimension(:,:,:)   :: pix, psi 
    integer(i4b), allocatable, dimension(:,:)     :: flag
    logical(lgt)       :: correct_sl
    character(len=512) :: prefix, postfix, prefix4D, filename
    character(len=2048) :: Sfilename
    character(len=4)   :: ctext, myid_text
    character(len=6)   :: samptext, scantext
    character(len=512), allocatable, dimension(:) :: slist
    type(shared_1d_int) :: sprocmask, sprocmask2
    real(sp),           allocatable, dimension(:,:,:,:) :: map_sky
    type(shared_2d_dp) :: sA_map
    type(shared_3d_dp) :: sb_map, sb_mono
    class(comm_map), pointer :: condmap
    class(map_ptr), allocatable, dimension(:) :: outmaps

    !if (iter > 1) self%first_call = .false.
    call int2string(iter, ctext)
    call update_status(status, "tod_start"//ctext)

    t_tot   = 0.d0
    call wall_time(t5)


    ! Set up full-sky map structures
    call wall_time(t1)
    if (self%output_aux_maps <= 0) then
       self%output_n_maps = 3
    else
       if (mod(iter-1,self%output_aux_maps) == 0) then
          self%output_n_maps = 7
       else
          self%output_n_maps = 3
       end if
    end if

    correct_sl      = .true.
    chisq_threshold = 7.d0
    n_main_iter     = 4
    chisq_threshold = 9.d0  !3000.d0
    !this ^ should be 7.d0, is currently 2000 to debug sidelobes
    ndet            = self%ndet
    nhorn           = self%nhorn
    ndelta          = size(delta,3)
    nside           = map_out%info%nside
    nmaps           = map_out%info%nmaps
    npix            = 12*nside**2
    nscan_tot       = self%nscan_tot
    chunk_size      = npix/self%numprocs_shared
    if (chunk_size*self%numprocs_shared /= npix) chunk_size = chunk_size+1
    allocate(A_abscal(self%ndet), b_abscal(self%ndet))
    allocate(map_sky(nmaps,self%nobs,0:ndet,ndelta))
    allocate(chisq_S(ndet,ndelta))
    allocate(slist(self%nscan))
    allocate(dipole_mod(nscan_tot, ndet))
    slist   = ''

    call int2string(chain, ctext)
    call int2string(iter, samptext)
    prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
    postfix = '_c' // ctext // '_k' // samptext // '.fits'

    if (self%myid == 0) then
       do k = 1, ndelta
          write(*,*) 'delta in =', real(delta(:,1,k),sp)
       end do
    end if

    if (.false.) then
       do i = 1, self%ndet
          filename = trim(chaindir) // '/BP_fg_' // trim(self%label(i)) // '_v1.fits'
          call map_in(i,1)%p%writeFITS(filename)
          filename = trim(chaindir) // '/BP_fg_' // trim(self%label(i)) // '_v2.fits'
          call map_in(i,2)%p%writeFITS(filename)
       end do
       call mpi_finalize(ierr)
       stop
!!$       deallocate(A_abscal, chisq_S, slist)
!!$       return
    end if

    ! Distribute fullsky maps
    allocate(m_buf(0:npix-1,nmaps))
    do j = 1, ndelta
       do i = 1, self%ndet
          map_in(i,j)%p%map = 1d-6 * map_in(i,j)%p%map ! uK to K
          call map_in(i,j)%p%bcast_fullsky_map(m_buf)
          do k = 1, self%nobs
             map_sky(:,k,i,j) = m_buf(self%ind2pix(k),:) ! uK to K
          end do
       end do
       do k = 1, self%nobs
          do l = 1, nmaps
             map_sky(l,k,0,j) = sum(map_sky(l,k,1:ndet,j))/ndet
          end do
       end do
    end do
    deallocate(m_buf)

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
    if (self%myid == 0) write(*,*) 'Precomputing sidelobe convolved sky'
    do i = 1, self%ndet
       if (.not. correct_sl) exit

       !TODO: figure out why this is rotated
       call map_in(i,1)%p%YtW()  ! Compute sky a_lms
       self%slconv(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
            & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
            & 100, 3, 100, self%slbeam(i)%p, map_in(i,1)%p, 2)
!       write(*,*) i, 'b', sum(abs(self%slconv(i)%p%c%a))
    end do
    call wall_time(t2); t_tot(13) = t2-t1

    call update_status(status, "tod_init")

    ! Compute output map and rms
    call wall_time(t3)
    do_oper             = .true.
    do main_iter = 1, n_main_iter

       call wall_time(t7)
       call update_status(status, "tod_istart")

       if (self%myid == 0) write(*,*) '  Performing main iteration = ', main_iter

       ! Select operations for current iteration
       do_oper(samp_acal)    = (main_iter == n_main_iter-3) .and. .not. self%first_call
       do_oper(samp_rcal)    = (main_iter == n_main_iter-2) .and. .not. self%first_call
       do_oper(samp_G)       = (main_iter == n_main_iter-1) .and. .not. self%first_call
       do_oper(samp_N)       = (main_iter >= n_main_iter-0)
       do_oper(samp_N_par)   = do_oper(samp_N)
       do_oper(prep_relbp)   = ndelta > 1 .and. (main_iter == n_main_iter-0) .and. .not. self%first_call !.and. mod(iter,2) == 0
       do_oper(prep_absbp)   = .false. !ndelta > 1 .and. (main_iter == n_main_iter-0) .and. .not. self%first_call .and. mod(iter,2) == 1
       do_oper(samp_bp)      = ndelta > 1 .and. (main_iter == n_main_iter-0) .and. .not. self%first_call
       do_oper(samp_mono)    = .false. !do_oper(bin_map)             !.and. .not. self%first_call
       do_oper(bin_map)      = (main_iter == n_main_iter  )
       do_oper(sel_data)     = .false. !(main_iter == n_main_iter  ) .and.       self%first_call
       do_oper(calc_chisq)   = (main_iter == n_main_iter  )
       do_oper(sub_sl)       = correct_sl
       do_oper(sub_zodi)     = self%subtract_zodi
       do_oper(output_slist) = mod(iter, 1) == 0

       dipole_mod = 0

       ! Perform pre-loop operations
       if (do_oper(bin_map) .or. do_oper(prep_relbp)) then
          if (do_oper(prep_relbp)) then
             ncol = nmaps + ndet - 1
             n_A  = nmaps*(nmaps+1)/2 + 4*(ndet-1)
             nout = self%output_n_maps + ndelta - 1
          else
             ncol = nmaps
             n_A  = nmaps*(nmaps+1)/2
             nout = self%output_n_maps
          end if
          allocate(outmaps(nout))
          do i = 1, nout
             outmaps(i)%p => comm_map(map_out%info)
          end do
          if (allocated(A_map)) deallocate(A_map, b_map)
          allocate(A_map(n_A,self%nobs), b_map(nout,ncol,self%nobs))
          A_map = 0.d0; b_map = 0.d0
          if (sA_map%init) call dealloc_shared_2d_dp(sA_map)
          call init_shared_2d_dp(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, [n_A,npix], sA_map)
          call mpi_win_fence(0, sA_map%win, ierr)
          if (sA_map%myid_shared == 0) sA_map%a = 0.d0
          call mpi_win_fence(0, sA_map%win, ierr)
          if (sb_map%init) call dealloc_shared_3d_dp(sb_map)
          call init_shared_3d_dp(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, [nout,ncol,npix], sb_map)
          call mpi_win_fence(0, sb_map%win, ierr)
          if (sb_map%myid_shared == 0) sb_map%a = 0.d0
          call mpi_win_fence(0, sb_map%win, ierr)
       end if

       if (do_oper(samp_mono)) then
          allocate(b_mono(nmaps,self%nobs,ndet), sys_mono(nmaps,nmaps+ndet,0:self%info%np-1))
          b_mono = 0.d0
          call init_shared_3d_dp(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, [nmaps,npix,ndet], sb_mono)
          call mpi_win_fence(0, sb_mono%win, ierr)
          if (sb_map%myid_shared == 0) sb_mono%a = 0.d0
          call mpi_win_fence(0, sb_mono%win, ierr)
       end if

       if (do_oper(sel_data)) then
          naccept = 0; ntot    = 0
       end if

       if (do_oper(samp_acal) .or. do_oper(samp_rcal)) then
          A_abscal = 0.d0; b_abscal = 0.d0
       end if

       if (do_oper(prep_absbp)) then
          chisq_S = 0.d0
       end if

       call wall_time(t8); t_tot(19) = t_tot(19) + t8-t7

      !  print *, "Total number of scans:", self%nscan
      !  print *, "Total number of tot scans:", self%nscan_tot
       ! Perform main analysis loop
       do i = 1, self%nscan

          !call update_status(status, "tod_loop1")
          call wall_time(t7)

          if (.not. any(self%scans(i)%d%accept)) cycle

          ! Short-cuts to local variables
          call wall_time(t1)
          ntod = self%scans(i)%ntod
          ndet = self%ndet

          ! Set up local data structure for current scan
          allocate(n_corr(ntod, ndet))                 ! Correlated noise in V
          allocate(s_sl(ntod, ndet))                   ! Sidelobe in uKcmb
          allocate(s_sky(ntod, ndet))                  ! Sky signal in uKcmb
          allocate(s_sky_prop(ntod, ndet,2:ndelta))    ! Sky signal in uKcmb
          allocate(s_bp(ntod, ndet))                   ! Signal minus mean
          allocate(s_bp_prop(ntod, ndet, 2:ndelta))    ! Signal minus mean
          allocate(s_orb(ntod, ndet))                  ! Orbital dipole in uKcmb
          if (do_oper(samp_mono)) allocate(s_mono(ntod, ndet))                 ! Monopole correction in uKcmb
          allocate(s_buf(ntod, ndet))                  ! Buffer
          allocate(s_tot(ntod, ndet))                  ! Sum of all sky compnents
          if (do_oper(sub_zodi)) allocate(s_zodi(ntod, ndet))                 ! Zodical light
          allocate(mask(ntod, ndet))                   ! Processing mask in time
          allocate(mask2(ntod, ndet))                  ! Processing mask in time
          allocate(pix(ntod, ndet, nhorn))                    ! Decompressed pointing
          allocate(psi(ntod, ndet, nhorn))                    ! Decompressed pol angle
          allocate(flag(ntod, ndet))                   ! Decompressed flags

          ! Initializing arrays to zero
          call wall_time(t2); t_tot(18) = t_tot(18) + t2-t1

          ! --------------------
          ! Analyze current scan
          ! --------------------

          ! Decompress pointing, psi and flags for current scan
          call wall_time(t1)
          do j = 1, ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             call self%decompress_pointing_and_flags(i, j, pix(:,j,:), &
                  & psi(:,j,:), flag(:,j))
          end do
          call self%symmetrize_flags(flag)
          !call validate_psi(self%scanid(i), psi)
          call wall_time(t2); t_tot(11) = t_tot(11) + t2-t1
          !call update_status(status, "tod_decomp")

          ! Construct sky signal template
          call wall_time(t1)
          if (do_oper(bin_map) .or. do_oper(prep_relbp)) then
             call project_sky(self, map_sky(:,:,:,1), pix, psi, flag, &
                  & sprocmask%a, i, s_sky, mask, s_bp=s_bp)
          else
             call project_sky(self, map_sky(:,:,:,1), pix, psi, flag, &
                  & sprocmask%a, i, s_sky, mask)
          end if
          if (do_oper(prep_relbp)) then
             do j = 2, ndelta
                call project_sky(self, map_sky(:,:,:,j), pix, psi, flag, &
                     & sprocmask2%a, i, s_sky_prop(:,:,j), mask2, s_bp=s_bp_prop(:,:,j))
             end do
          else if (do_oper(prep_absbp)) then
             do j = 2, ndelta
                call project_sky(self, map_sky(:,:,:,j), pix, psi, flag, &
                     & sprocmask2%a, i, s_sky_prop(:,:,j), mask2)
             end do
          end if
          if (main_iter == 1 .and. self%first_call) then
             do j = 1, ndet
                if (all(mask(:,j) == 0)) self%scans(i)%d(j)%accept = .false.
                if (self%scans(i)%d(j)%sigma0 <= 0.d0) self%scans(i)%d(j)%accept = .false.
             end do
          end if
          do j = 1, ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             if (self%scans(i)%d(j)%sigma0 <= 0) write(*,*) main_iter, self%scanid(i), j, self%scans(i)%d(j)%sigma0
          end do
          call wall_time(t2); t_tot(1) = t_tot(1) + t2-t1
          !call update_status(status, "tod_project")

          ! Construct orbital dipole template
          call wall_time(t1)
          call self%orb_dp%p%compute_orbital_dipole_4pi(i, pix(:,:,1), psi(:,:,1), s_orb)
          call wall_time(t2); t_tot(2) = t_tot(2) + t2-t1
          !call update_status(status, "tod_orb")

         !  call wall_time(t9)
          ! Construct zodical light template
          if (do_oper(sub_zodi)) then
             call compute_zodi_template(self%nside, pix(:,:,1), self%scans(i)%satpos, [30.d9, 30.d9, 30.d9, 30.d9], s_zodi)
          end if
         !  call wall_time(t10)
         !  print *, "Zodi template took :", t10-t9, "sec"

          ! Construct sidelobe template
          call wall_time(t1)
          if (do_oper(sub_sl)) then
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                call self%construct_sl_template(self%slconv(j)%p, &
                     & pix(:,j,1), psi(:,j,1), s_sl(:,j), self%polang(j))
                s_sl(:,j) = 2 * s_sl(:,j) ! Scaling by a factor of 2, by comparison with LevelS. Should be understood
             end do
          else
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                s_sl(:,j) = 0.
             end do
          end if
          call wall_time(t2); t_tot(12) = t_tot(12) + t2-t1

          ! Construct monopole correction template
          call wall_time(t1)
          if (do_oper(samp_mono)) then
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                !s_mono(:,j) = self%mono(j)
                s_mono(:,j) = 0.d0 ! Disabled for now
             end do
          end if
          do j = 1, ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             s_tot(:,j) = s_sky(:,j) + s_sl(:,j) + s_orb(:,j)
             if (do_oper(samp_mono)) s_tot(:,j) = s_tot(:,j) + s_mono(:,j)
          end do
          call wall_time(t2); t_tot(1) = t_tot(1) + t2-t1

          ! Precompute filtered signal for calibration
          if (do_oper(samp_G) .or. do_oper(samp_rcal) .or. do_oper(samp_acal)) then
             call self%downsample_tod(s_orb(:,1), ext)
             allocate(s_invN(ext(1):ext(2), ndet))      ! s * invN
             allocate(s_lowres(ext(1):ext(2), ndet))      ! s * invN
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                if (do_oper(samp_G) .or. do_oper(samp_rcal) .or. .not. self%orb_abscal) then
                   s_buf(:,j) = s_tot(:,j)
                   call fill_all_masked(s_buf(:,j), mask(:,j), ntod, trim(self%operation)=='sample', real(self%scans(i)%d(j)%sigma0, sp), handle, self%scans(i)%chunk_num)
                   call self%downsample_tod(s_buf(:,j), ext, &
                        & s_lowres(:,j))!, mask(:,j))
                else
                   call self%downsample_tod(s_orb(:,j), ext, &
                        & s_lowres(:,j))!, mask(:,j))
                end if
             end do
             s_invN = s_lowres
             call multiply_inv_N(self, i, s_invN,   sampfreq=self%samprate_lowres, pow=0.5d0)
             call multiply_inv_N(self, i, s_lowres, sampfreq=self%samprate_lowres, pow=0.5d0)
          end if

          ! Prepare for absolute calibration
          if (do_oper(samp_acal) .or. do_oper(samp_rcal)) then
             call wall_time(t1)
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle

                if (do_oper(samp_acal)) then
                   if (self%orb_abscal) then
                      s_buf(:, j) = real(self%gain0(0),sp) * (s_tot(:, j) - s_orb(:, j)) + &
                           & real(self%gain0(j) + self%scans(i)%d(j)%dgain,sp) * s_tot(:, j)
                   else
                      s_buf(:, j) = real(self%gain0(j) + self%scans(i)%d(j)%dgain,sp) * s_tot(:, j)
                   end if
                else
!                   s_buf(:,j) =  s_tot(:,j) - s_orb(:,j) !s_sky(:,j) + s_sl(:,j) + s_mono(:,j)
                   s_buf(:,j) = real(self%gain0(0) + self%scans(i)%d(j)%dgain,sp) * s_tot(:, j)
                end if
             end do
             call accumulate_abscal(self, i, mask, s_buf, s_lowres, s_invN, A_abscal, b_abscal, handle, .false.)

             if (.false.) then
                call int2string(self%scanid(i), scantext)
                open(78,file='tod_'//trim(self%label(j))//'_pid'//scantext//'_k'//samptext//'.dat', recl=1024)
                write(78,*) "# Sample     Data (V)    Res (V)    s_sub (K)   s_orb (K)   mask"
                do k = 1, ntod, 60
                   write(78,*) k, mean(1.d0*self%scans(i)%d(j)%tod(k:k+59)), mean(1.d0*self%scans(i)%d(j)%tod(k:k+59) - &
                        & self%scans(i)%d(j)%gain*s_buf(k:k+59,j)), mean(1.d0*s_orb(k:k+59,j)),  mean(1.d0*s_buf(k:k+59,j)),  minval(mask(k:k+59,j))
                end do
                close(78)
             end if
             call wall_time(t2); t_tot(14) = t_tot(14) + t2-t1
          end if

          ! Fit gain
          if (do_oper(samp_G)) then
             call wall_time(t1)
             call calculate_gain_mean_std_per_scan(self, i, s_invN, mask, s_lowres, s_tot, handle)
             call wall_time(t2); t_tot(4) = t_tot(4) + t2-t1
          end if

          ! Fit correlated noise
          if (do_oper(samp_N)) then
             call wall_time(t1)
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                if (do_oper(samp_mono)) then
                   s_buf(:,j) = s_tot(:,j)-s_mono(:,j)
                else
                   s_buf(:,j) = s_tot(:,j)
                end if
             end do
             call sample_n_corr(self, handle, i, mask, s_buf, n_corr, pix(:,:,1), .true.)
!!$             do j = 1, ndet
!!$                n_corr(:,j) = sum(n_corr(:,j))/ size(n_corr,1)
!!$             end do
             call wall_time(t2); t_tot(3) = t_tot(3) + t2-t1
          else
             n_corr = 0.
          end if

          ! Compute noise spectrum
          if (do_oper(samp_N_par)) then
             call wall_time(t1)
             call sample_noise_psd(self, handle, i, mask, s_tot, n_corr)
             call wall_time(t2); t_tot(6) = t_tot(6) + t2-t1
          end if
          ! Compute chisquare
          if (do_oper(calc_chisq)) then
             call wall_time(t1)
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                s_buf(:,j) =  s_sl(:,j) + s_orb(:,j)
                if (do_oper(samp_mono)) s_buf(:,j) =  s_buf(:,j) + s_mono(:,j)
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
                if (count(iand(flag(:,j),self%flag0) .ne. 0) > 0.1*ntod) then
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
             if (any(.not. self%scans(i)%d%accept)) self%scans(i)%d%accept = .false.
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) then
                   self%scans(i)%d(self%partner(j))%accept = .false.
                end if
             end do
             call wall_time(t2); t_tot(15) = t_tot(15) + t2-t1
          end if

          ! Compute chisquare for bandpass fit
          if (do_oper(prep_absbp)) then
             call wall_time(t1)
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                s_buf(:,j) =  s_sl(:,j) + s_orb(:,j)
                if (do_oper(samp_mono)) s_buf(:,j) =  s_buf(:,j) + s_mono(:,j)
                call self%compute_chisq(i, j, mask2(:,j), s_sky(:,j), &
                     & s_buf(:,j), n_corr(:,j), absbp=.true.)
                chisq_S(j,1) = chisq_S(j,1) + self%scans(i)%d(j)%chisq_prop
                do k = 2, ndelta
                   call self%compute_chisq(i, j, mask2(:,j), s_sky_prop(:,j,k), &
                        & s_buf(:,j), n_corr(:,j), absbp=.true.)
                   chisq_S(j,k) = chisq_S(j,k) + self%scans(i)%d(j)%chisq_prop
                end do
             end do
             call wall_time(t2); t_tot(7) = t_tot(7) + t2-t1
          end if


          ! Compute binned map
          if (do_oper(bin_map)) then
             call wall_time(t1)
             allocate(d_calib(nout,ntod, ndet))
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                inv_gain = 1.0 / real(self%scans(i)%d(j)%gain,sp)
                d_calib(1,:,j) = (self%scans(i)%d(j)%tod - n_corr(:,j)) * &
                     & inv_gain - s_tot(:,j) + s_sky(:,j) - s_bp(:,j)
                if (do_oper(bin_map) .and. nout > 1) d_calib(2,:,j) = d_calib(1,:,j) - s_sky(:,j) + s_bp(:,j) ! Residual
                if (do_oper(bin_map) .and. nout > 2) d_calib(3,:,j) = (n_corr(:,j) - sum(n_corr(:,j)/ntod)) * inv_gain
                if (do_oper(bin_map) .and. nout > 3) d_calib(4,:,j) = s_bp(:,j)
                if (do_oper(bin_map) .and. nout > 4) then
                   if (do_oper(samp_mono)) then
                      d_calib(5,:,j) = s_mono(:,j)
                   else
                      d_calib(5,:,j) = 0.d0
                   end if
                end if
                if (do_oper(bin_map) .and. nout > 5) d_calib(6,:,j) = s_orb(:,j)
                if (do_oper(bin_map) .and. nout > 6) d_calib(7,:,j) = s_sl(:,j)
                if (do_oper(bin_map) .and. nout > 7 .and. do_oper(sub_zodi)) d_calib(8,:,j) = s_zodi(:,j)

                if (do_oper(prep_relbp)) then
                   do k = 2, ndelta
                      d_calib(self%output_n_maps+k-1,:,j) = d_calib(1,:,j) + s_bp(:,j) - s_bp_prop(:,j,k)
                   end do
                end if
             end do

             if (do_oper(bin_map) .and. self%output_4D_map > 0 .and. mod(iter-1,self%output_4D_map) == 0) then

                ! Output 4D map; note that psi is zero-base in 4D maps, and one-base in Commander
                call int2string(self%myid, myid_text)
                call output_4D_maps_hdf(trim(chaindir) // '/tod_4D_chain'//ctext//'_proc' // myid_text // '.h5', &
                     & samptext, self%scanid(i), self%nside, self%npsi, &
                     & self%label, self%horn_id, real(self%polang*180/pi,sp), &
                     & real(self%scans(i)%d%sigma0/self%scans(i)%d%gain,sp), &
                     & pix(:,:,1), psi(:,:,1)-1, d_calib(1,:,:), iand(flag,self%flag0), &
                     & self%scans(i)%d(:)%accept)
!                prefix4D = "!"// trim(chaindir) // '/.tod_' // trim(self%freq) // '_' // '4D_pid' // scantext
!                call output_4D_maps(prefix4D, postfix, self%scanid(i), self%nside, self%npsi, &
!                     & self%label, self%horn_id, real(self%polang*180/pi,sp), &
!                     & real(self%scans(i)%d%sigma0/self%scans(i)%d%gain,sp), &
!                     & pix, psi-1, d_calib(1,:,:), iand(flag,self%flag0), &
!                     & self%scans(i)%d(:)%accept)
             end if


             call wall_time(t2); t_tot(5) = t_tot(5) + t2-t1

             if (.false. .and. do_oper(bin_map) ) then
                call int2string(self%scanid(i), scantext)
                do k = 1, self%ndet
                   open(78,file='tod_'//trim(self%label(k))//'_pid'//scantext//'.dat', recl=1024)
                   write(78,*) "# Sample     Data (V)     Mask    cal_TOD (K)   res (K)"// &
                        & "   n_corr (K)   s_corr (K)   s_mono (K)   s_orb  (K)   s_sl (K)"
                   do j = 1, ntod
                      write(78,*) j, self%scans(i)%d(k)%tod(j), mask(j,1), d_calib(:,j,k)
                   end do
                   close(78)
                end do
             end if


             call wall_time(t1)

             if (do_oper(samp_mono)) then
                call bin_TOD(self, d_calib, pix(:,:,1), &
                     & psi(:,:,1), flag, A_map, b_map, i, do_oper(prep_relbp), b_mono=b_mono)
             else
                call bin_TOD(self, d_calib, pix(:,:,1), &
                     & psi(:,:,1), flag, A_map, b_map, i, do_oper(prep_relbp))
             end if
             deallocate(d_calib)
             call wall_time(t2); t_tot(8) = t_tot(8) + t2-t1
          end if

          do j = 1, ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             dipole_mod(self%scanid(i), j) = masked_variance(s_sky(:, j), mask(:, j))
          end do

          !----------------------------------------------------------------------------------
          ! Calling Simulation Routine
          !write(*,*) "Debug Message before simulation routine."
          if (self%enable_tod_simulations) then !.and. (main_iter == 1)) then
            call wall_time(t1)
            call self%simulate_LFI_tod(i, s_tot, handle)
            call wall_time(t2); t_tot(23) = t_tot(23) + t2-t1
          end if
          !----------------------------------------------------------------------------------
          ! Clean up
          call wall_time(t1)
          deallocate(n_corr, s_sl, s_sky, s_orb, s_tot, s_buf)
          deallocate(s_bp, s_sky_prop, s_bp_prop)
          deallocate(mask, mask2, pix, psi, flag)
          if (allocated(s_lowres)) deallocate(s_lowres)
          if (allocated(s_invN)) deallocate(s_invN)
          if (do_oper(sub_zodi)) deallocate(s_zodi)
          if (do_oper(samp_mono)) deallocate(s_mono)
          call wall_time(t2); t_tot(18) = t_tot(18) + t2-t1

          call wall_time(t8); t_tot(19) = t_tot(19) + t8-t7

          if (do_oper(output_slist)) then
             self%scans(i)%proctime   = self%scans(i)%proctime   + t8-t7
             self%scans(i)%n_proctime = self%scans(i)%n_proctime + 1
             if (main_iter == n_main_iter) then
                write(slist(i),*) self%scanid(i), '"',trim(self%hdfname(i)), &
                     & '"', real(self%scans(i)%proctime/self%scans(i)%n_proctime,sp),real(self%spinaxis(i,:),sp)
             end if
          end if
          !call update_status(status, "tod_loop2")

       end do

       call mpi_allreduce(mpi_in_place, dipole_mod, size(dipole_mod), MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)

!       if (self%myid == 0) then
!          write(*, *) "CHECKPOINT"
!          write(*, *) dipole_mod
!       end if

       if (do_oper(samp_acal)) then
          call wall_time(t1)
          call sample_abscal_from_orbital(self, handle, A_abscal, b_abscal)
          call wall_time(t2); t_tot(16) = t_tot(16) + t2-t1
       end if

       if (do_oper(samp_rcal)) then
          call wall_time(t1)
          call sample_relcal(self, handle, A_abscal, b_abscal)
          call wall_time(t2); t_tot(16) = t_tot(16) + t2-t1
       end if

       if (do_oper(samp_G)) then
          call wall_time(t1)
          call sample_smooth_gain(self, handle, dipole_mod)
          call wall_time(t2); t_tot(4) = t_tot(4) + t2-t1
       end if

       call wall_time(t7)

    end do
    call wall_time(t4)

    ! Output latest scan list with new timing information
    if (do_oper(output_slist)) then
       call update_status(status, "scanlist1")
       call wall_time(t1)
       call self%output_scan_list(slist)
       call wall_time(t2); t_tot(20) = t_tot(20) + t2-t1
       call update_status(status, "scanlist2")
    end if

    ! Solve combined map, summed over all pixels
    call wall_time(t1)
    call update_status(status, "shared1")
    if (sA_map%init) then
       do i = 0, self%numprocs_shared-1
          !call update_status(status, "tod_share_loop")
          start_chunk = mod(self%myid_shared+i,self%numprocs_shared)*chunk_size
          end_chunk   = min(start_chunk+chunk_size-1,npix-1)
          do while (start_chunk < npix)
             if (self%pix2ind(start_chunk) /= -1) exit
             start_chunk = start_chunk+1
          end do
          do while (end_chunk >= start_chunk)
             if (self%pix2ind(end_chunk) /= -1) exit
             end_chunk = end_chunk-1
          end do
          if (start_chunk < npix)       start_chunk = self%pix2ind(start_chunk)
          if (end_chunk >= start_chunk) end_chunk   = self%pix2ind(end_chunk)

          call mpi_win_fence(0, sA_map%win, ierr)
          call mpi_win_fence(0, sb_map%win, ierr)
          do j = start_chunk, end_chunk
             sA_map%a(:,self%ind2pix(j)+1) = sA_map%a(:,self%ind2pix(j)+1) + &
                  & A_map(:,j)
             sb_map%a(:,:,self%ind2pix(j)+1) = sb_map%a(:,:,self%ind2pix(j)+1) + &
                  & b_map(:,:,j)
          end do
          if (do_oper(samp_mono)) then
             call mpi_win_fence(0, sb_mono%win, ierr)
             do j = start_chunk, end_chunk
                sb_mono%a(:,self%ind2pix(j)+1,:) = sb_mono%a(:,self%ind2pix(j)+1,:) + &
                     & b_mono(:,j,:)
             end do
          end if
       end do
       call mpi_win_fence(0, sA_map%win, ierr)
       call mpi_win_fence(0, sb_map%win, ierr)
       if (do_oper(samp_mono)) call mpi_win_fence(0, sb_mono%win, ierr)

       call update_status(status, "shared2")

       call update_status(status, "finalize1")
       Sfilename = trim(prefix) // 'Smap'// trim(postfix)
       if (do_oper(samp_mono)) then
          if (do_oper(prep_relbp)) then
             call finalize_binned_map(self, handle, sA_map, sb_map, rms_out, outmaps=outmaps, sb_mono=sb_mono, sys_mono=sys_mono, chisq_S=chisq_S, Sfile=Sfilename, mask=sprocmask2%a)
          else
             call finalize_binned_map(self, handle, sA_map, sb_map, rms_out, outmaps=outmaps, sb_mono=sb_mono, sys_mono=sys_mono)
          end if
!!$          condmap => comm_map(self%info)
!!$          call self%finalize_binned_map(handle, sA_map, sb_map, outmaps, rms_out, sb_mono=sb_mono, sys_mono=sys_mono, condmap=condmap)
!!$          call condmap%writeFITS("cond.fits")
!!$          call condmap%dealloc()
       else
          !condmap => comm_map(self%info)
          if (do_oper(prep_relbp)) then
             call finalize_binned_map(self, handle, sA_map, sb_map, rms_out, outmaps=outmaps, chisq_S=chisq_S, Sfile=Sfilename, mask=sprocmask2%a)
          else
             call finalize_binned_map(self, handle, sA_map, sb_map, rms_out, outmaps=outmaps)
          end if
          !call condmap%writeFITS("cond.fits")
          !call condmap%dealloc()
       end if

       if (do_oper(samp_bp)) then
          call wall_time(t1)
          call sample_bp(self, iter, delta, map_sky, handle, chisq_S)
          call wall_time(t2); t_tot(17) = t_tot(17) + t2-t1
       end if

       call update_status(status, "finalize2")
       map_out%map = outmaps(1)%p%map

       ! Sample monopole coefficients
       if (do_oper(samp_mono)) then
          call sample_mono(self, handle, sys_mono, outmaps(2)%p, rms_out, &
               & self%procmask)
       end if

       ! Update bandpass parameters
       self%bp_delta = delta(:,:,1)

       ! Output maps to disk
       if (.false. .and. trim(self%freq) == '030') then
          if (self%myid == 0) write(*,*) 'Boosting rms 5x'
          rms_out%map = 5*rms_out%map 
       end if
       call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
       call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))

       if (self%output_n_maps > 1) call outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
       if (self%output_n_maps > 2) call outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
       if (self%output_n_maps > 3) call outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix))
       if (self%output_n_maps > 4) call outmaps(5)%p%writeFITS(trim(prefix)//'mono'//trim(postfix))
       if (self%output_n_maps > 5) call outmaps(6)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
       if (self%output_n_maps > 6) call outmaps(7)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
       if (self%output_n_maps > 7) call outmaps(8)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))

       call update_status(status, "finalize3")

    end if

    if (self%first_call) then
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
       if (self%enable_tod_simulations) write(*,*) '  Time tod sims   = ', t_tot(23)
       write(*,*) '  Time scanlist   = ', t_tot(20)
       write(*,*) '  Time final      = ', t_tot(10)
       if (self%first_call) then
!!$          write(*,*) '  Time total      = ', t6-t5, &
!!$               & ', accept rate = ', real(naccept,sp) / ntot
       else
          write(*,*) '  Time total      = ', t6-t5, sum(t_tot(1:18))
       end if
    end if


    ! Clean up temporary arrays
    deallocate(A_abscal, b_abscal, chisq_S)
    if (allocated(A_map)) deallocate(A_map, b_map)
    if (sA_map%init)  call dealloc_shared_2d_dp(sA_map)
    if (sb_map%init)  call dealloc_shared_3d_dp(sb_map)
    if (sb_mono%init) call dealloc_shared_3d_dp(sb_mono)
    if (allocated(b_mono)) deallocate(b_mono)
    if (allocated(sys_mono)) deallocate(sys_mono)
    if (allocated(slist)) deallocate(slist)
    if (allocated(dipole_mod)) deallocate(dipole_mod)

    if (allocated(outmaps)) then
       do i = 1, nout
          call outmaps(i)%p%dealloc
       end do
       deallocate(outmaps)
    end if

    if (sprocmask%init)  call dealloc_shared_1d_int(sprocmask)
    if (sprocmask2%init) call dealloc_shared_1d_int(sprocmask2)
    deallocate(map_sky)

    if (correct_sl) then
       do i = 1, self%ndet
          call self%slconv(i)%p%dealloc(); deallocate(self%slconv(i)%p)
       end do
    end if

    call int2string(iter, ctext)
    call update_status(status, "tod_end"//ctext)

    ! Parameter to check if this is first time routine has been
    self%first_call = .false.

  end subroutine process_LFI_tod

  ! ************************************************
  !
  !> @brief Commander3 native simulation module. It
  !! simulates correlated noise and then rewrites
  !! the original timestreams inside the files.
  !
  !> @author Maksym Brilenkov
  !
  !> @param[in]
  !> @param[out]
  !
  ! ************************************************
  subroutine simulate_LFI_tod(self, scan_id, s_tot, handle)
    use omp_lib
    implicit none
    class(comm_LFI_tod), intent(inout) :: self
    ! Parameter file variables
    !type(comm_params),                     intent(in)    :: cpar
    ! Other input/output variables
    real(sp), allocatable, dimension(:,:), intent(in)    :: s_tot   !< total sky signal
    integer(i4b),                          intent(in)    :: scan_id !< current PID
    type(planck_rng),                      intent(inout) :: handle
    ! Simulation variables
    real(sp), allocatable, dimension(:,:) :: tod_per_detector !< simulated tods per detector
    real(sp)                              :: gain   !< detector's gain value
    real(sp)                              :: sigma0
    real(sp) :: nu_knee
    real(sp) :: alpha
    real(sp) :: samprate
    real(sp) :: fft_norm
    integer(i4b)                          :: ntod !< total amount of ODs
    integer(i4b)                          :: ndet !< total amount of detectors
    ! HDF5 variables
    character(len=512) :: mystring, mysubstring !< dummy values for string manipulation
    integer(i4b)       :: myindex     !< dummy value for string manipulation
    character(len=512) :: currentHDFFile !< hdf5 file which stores simulation output
    character(len=6)   :: pidLabel
    character(len=3)   :: detectorLabel
    type(hdf_file)     :: hdf5_file   !< hdf5 file to work with
    integer(i4b)       :: hdf5_error  !< hdf5 error status
    integer(HID_T)     :: hdf5_file_id !< File identifier
    integer(HID_T)     :: dset_id     !< Dataset identifier
    integer(HSIZE_T), dimension(1) :: dims
    ! Other variables
    integer(i4b)                          :: i, j, k !< loop variables
    integer(i4b)       :: mpi_err !< MPI error status
    integer(i4b)       :: nomp !< Number of threads available
    integer(i4b)       :: omp_err !< OpenMP error status
    integer(i4b) :: omp_get_max_threads
    integer(i4b) :: n, nfft
    integer*8    :: plan_back
    real(sp) :: nu
    real(sp), allocatable, dimension(:,:) :: n_corr
    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv
    character(len=10) :: processor_label   !< to have a nice output to screen

    ! shortcuts
    ntod = self%scans(scan_id)%ntod
    ndet = self%ndet

    ! Simulating 1/f noise
    !write(*,*) "Simulating correlated noise"
    nfft = 2 * ntod
    n = nfft / 2 + 1
    nomp = omp_get_max_threads()
    call sfftw_init_threads(omp_err)
    call sfftw_plan_with_nthreads(nomp)
    ! planning FFTW - in principle we should do both forward and backward FFTW,
    ! but in this case we can omit forward one and go directly with backward to
    ! save some time on a whole operation.
    allocate(dt(nfft), dv(0:n-1))
    call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)

    !$OMP PARALLEL PRIVATE(i, j, k, dt, dv, sigma0, nu, nu_knee, alpha)
    allocate(dt(nfft), dv(0:n-1), n_corr(ntod, ndet))
    !$OMP DO SCHEDULE(guided)
    do j = 1, ndet
      ! skipping iteration if scan was not accepted
      if (.not. self%scans(scan_id)%d(j)%accept) cycle
      ! getting gain for each detector (units, V / K)
      ! (gain is assumed to be CONSTANT for EACH SCAN)
      gain   = self%scans(scan_id)%d(j)%gain
      sigma0 = self%scans(scan_id)%d(j)%sigma0
      samprate = self%samprate
      alpha    = self%scans(scan_id)%d(j)%alpha
      ! knee frequency
      nu_knee  = self%scans(scan_id)%d(j)%fknee
      ! used when adding fluctuation terms to Fourier coeffs (depends on Fourier convention)
      fft_norm = sqrt(1.d0 * nfft)
      !
      !dv(0) = dv(0) + fft_norm * sigma0 * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0)
      dv(0) = fft_norm * sigma0 * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0)
      do k = 1, (n - 1)
        nu = k * (samprate / 2) / (n - 1)
        !dv(k) = sigma0 * cmplx(rand_gauss(handle), rand_gauss(handle)) * sqrt(1 + (nu / nu_knee)**alpha) /sqrt(2          .0)
        dv(k) = sigma0 * cmplx(rand_gauss(handle), rand_gauss(handle)) * sqrt((nu / nu_knee)**alpha) /sqrt(2.0)
      end do
      ! Executing Backward FFT
      call sfftw_execute_dft_c2r(plan_back, dv, dt)
      dt = dt / nfft
      n_corr(:, j) = dt(1:ntod)
      !write(*,*) "n_corr ", n_corr(:, j)
    end do
    !$OMP END DO
    deallocate(dt, dv)
    !$OMP END PARALLEL

    call sfftw_destroy_plan(plan_back)

    ! Allocating main simulations' array
    allocate(tod_per_detector(ntod, ndet))       ! Simulated tod

    ! Main simulation loop
    do i = 1, ntod
      do j = 1, ndet
        ! skipping iteration if scan was not accepted
        if (.not. self%scans(scan_id)%d(j)%accept) cycle
        ! getting gain for each detector (units, V / K)
        ! (gain is assumed to be CONSTANT for EACH SCAN)
        gain   = self%scans(scan_id)%d(j)%gain
        !write(*,*) "gain ", gain
        sigma0 = self%scans(scan_id)%d(j)%sigma0
        !write(*,*) "sigma0 ", sigma0
        ! Simulating tods
        tod_per_detector(i,j) = gain * s_tot(i,j) + n_corr(i, j) + sigma0 * rand_gauss(handle)
        !tod_per_detector(i,j) = 0
      end do
    end do

    !----------------------------------------------------------------------------------
    ! Saving stuff to hdf file
    ! Getting the full path and name of the current hdf file to overwrite
    !----------------------------------------------------------------------------------
    mystring = trim(self%hdfname(scan_id))
    mysubstring = 'LFI_0'
    myindex = index(trim(mystring), trim(mysubstring))
    currentHDFFile = trim(self%sims_output_dir)//'/'//trim(mystring(myindex:))
    !write(*,*) "hdf5name "//trim(self%hdfname(scan_id))
    !write(*,*) "currentHDFFile "//trim(currentHDFFile)
    ! Converting PID number into string value
    call int2string(self%scanid(scan_id), pidLabel)
    call int2string(self%myid, processor_label)
    write(*,*) "Process: "//trim(processor_label)//" started writing PID: "//trim(pidLabel)//", into:"
    write(*,*) trim(currentHDFFile)
    ! For debugging
    !call MPI_Finalize(mpi_err)
    !stop

    dims(1) = ntod
    ! Initialize FORTRAN interface.
    call h5open_f(hdf5_error)
    ! Open an existing file - returns hdf5_file_id
    call  h5fopen_f(currentHDFFile, H5F_ACC_RDWR_F, hdf5_file_id, hdf5_error)
    do j = 1, ndet
      detectorLabel = self%label(j)
      ! Open an existing dataset.
      call h5dopen_f(hdf5_file_id, trim(pidLabel)//'/'//trim(detectorLabel)//'/'//'tod', dset_id, hdf5_error)
      ! Write tod data to a dataset
      call h5dwrite_f(dset_id, H5T_IEEE_F32LE, tod_per_detector(:,j), dims, hdf5_error)
      ! Close the dataset.
      call h5dclose_f(dset_id, hdf5_error)
    end do
    ! Close the file.
    call h5fclose_f(hdf5_file_id, hdf5_error)
    ! Close FORTRAN interface.
    call h5close_f(hdf5_error)

    !write(*,*) "hdf5_error",  hdf5_error
    ! freeing memory up
    deallocate(n_corr, tod_per_detector)
    write(*,*) "Process:", self%myid, "finished writing PID: "//trim(pidLabel)//"."

    ! lastly, we need to copy an existing filelist.txt into simulation folder
    ! and change the pointers to new files
    !if (self%myid == 0) then
    !  call system("cp "//trim(filelist)//" "//trim(simsdir))
    !  !mystring = filelist
    !  !mysubstring = ".txt"
    !  !myindex = index(trim(mystring), trim(mysubstring))
    !end if

    ! For debugging
    !call MPI_Finalize(mpi_err)
    !stop
  end subroutine simulate_LFI_tod

end module comm_tod_LFI_mod
