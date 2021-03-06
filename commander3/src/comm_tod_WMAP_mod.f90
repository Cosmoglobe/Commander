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
module comm_tod_WMAP_mod
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
  public comm_WMAP_tod

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
  integer(i4b), parameter :: sub_zodi    = 17
  logical(lgt), dimension(N_test) :: do_oper


  type, extends(comm_tod) :: comm_WMAP_tod
    class(orbdipole_pointer), allocatable :: orb_dp !orbital dipole calculator

   contains
     procedure     :: process_tod        => process_WMAP_tod
  end type comm_WMAP_tod

  interface comm_WMAP_tod
     procedure constructor
  end interface comm_WMAP_tod

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
    class(comm_WMAP_tod),    pointer    :: constructor

    integer(i4b) :: i,nside_beam, lmax_beam, nmaps_beam, ndelta
    character(len=512) :: datadir
    logical(lgt) :: pol_beam

    ! Set up WMAP specific parameters
    allocate(constructor)
    constructor%output_n_maps = 3
    constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
    constructor%nhorn    = 1

    ! Initialize beams
    nside_beam = 512
    nmaps_beam = 3
    pol_beam   = .true.
    constructor%nside_beam = nside_beam


    !initialize the common tod stuff
    call constructor%tod_constructor(cpar, id_abs, info, tod_type)
    
    !TODO: this is LFI specific, write something here for wmap
    call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label)
    do i = 1, constructor%ndet
       if (mod(i,2) == 1) then
          constructor%partner(i) = i+1
       else
          constructor%partner(i) = i-1
       end if
       constructor%horn_id(i) = (i+1)/2
    end do

    ! Read the actual TOD
    ! TODO: this probabl needs a seperate fucntion/ modifications for wmap
    call constructor%read_tod(constructor%label)

    datadir = trim(cpar%datadir)//'/'

    ! Initialize bandpass mean and proposal matrix
    call constructor%initialize_bp_covar(trim(datadir)//cpar%ds_tod_bp_init(id_abs))

    !load the instrument file 
    call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

    allocate(constructor%orb_dp)
    constructor%orb_dp%p => comm_orbdipole(constructor, constructor%mbeam)
  end function constructor

  !**************************************************
  !             Driver routine
  !**************************************************
  subroutine process_WMAP_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
    implicit none
    class(comm_WMAP_tod),                      intent(inout) :: self
    character(len=*),                         intent(in)    :: chaindir
    integer(i4b),                             intent(in)    :: chain, iter
    type(planck_rng),                         intent(inout) :: handle
    type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
    real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
    class(comm_map),                          intent(inout) :: map_out      ! Combined output map
    class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms

    integer(i4b) :: i, j, k, l, start_chunk, end_chunk, chunk_size, ntod, ndet
    integer(i4b) :: nside, npix, nmaps, naccept, ntot, ext(2), nscan_tot
    integer(i4b) :: ierr, main_iter, n_main_iter, ndelta, ncol, n_A, nout=1
    real(dp)     :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, chisq_threshold
    real(dp)     :: t_tot(22)
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
    integer(i4b), allocatable, dimension(:,:)     :: pix, psi, flag
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

    if (iter > 1) self%first_call = .false.
    call int2string(iter, ctext)
    call update_status(status, "tod_start"//ctext)

    t_tot   = 0.d0
    call wall_time(t5)


    ! Set up full-sky map structures
    call wall_time(t1)
    chisq_threshold = 7.d0
    ndet            = self%ndet
    ndelta          = size(delta,3)
    nside           = map_out%info%nside
    nmaps           = map_out%info%nmaps
    npix            = 12*nside**2
    nscan_tot       = self%nscan_tot
    chunk_size      = npix/self%numprocs_shared
    if (chunk_size*self%numprocs_shared /= npix) chunk_size = chunk_size+1
    allocate(map_sky(nmaps,self%nobs,0:ndet,ndelta))
    allocate(chisq_S(ndet,ndelta))

    call int2string(chain, ctext)
    call int2string(iter, samptext)
    prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
    postfix = '_c' // ctext // '_k' // samptext // '.fits'

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

    call update_status(status, "tod_init")

    ! Perform main analysis loop
    do i = 1, self%nscan

      !call update_status(status, "tod_loop1")
      if (.not. any(self%scans(i)%d%accept)) cycle

      ! Short-cuts to local variables
      call wall_time(t1)
      ntod = self%scans(i)%ntod
      ndet = self%ndet

      ! Set up local data structure for current scan
      allocate(n_corr(ntod, ndet))                 ! Correlated noise in V
      allocate(s_sky(ntod, ndet))                  ! Sky signal in uKcmb
      allocate(s_sky_prop(ntod, ndet,2:ndelta))    ! Sky signal in uKcmb
      allocate(s_orb(ntod, ndet))                  ! Orbital dipole in uKcmb
      allocate(s_buf(ntod, ndet))                  ! Buffer
      allocate(s_tot(ntod, ndet))                  ! Sum of all sky compnents
      allocate(mask(ntod, ndet))                   ! Processing mask in time
      allocate(mask2(ntod, ndet))                  ! Processing mask in time
      allocate(pix(ntod, ndet))                    ! Decompressed pointing
      allocate(psi(ntod, ndet))                    ! Decompressed pol angle
      allocate(flag(ntod, ndet))                   ! Decompressed flags

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
      call self%symmetrize_flags(flag)
      call wall_time(t2); t_tot(11) = t_tot(11) + t2-t1
      !call update_status(status, "tod_decomp")

      ! Construct sky signal template
      !TODO: we will want this but I don't know if it will be the same code or
      !not
      !call wall_time(t1)
      !if (do_oper(bin_map) .or. do_oper(prep_relbp)) then
      !   call project_sky(self, map_sky(:,:,:,1), pix, psi, flag, &
      !        & sprocmask%a, i, s_sky, mask, s_bp=s_bp)
      !else
      !   call project_sky(self, map_sky(:,:,:,1), pix, psi, flag, &
      !        & sprocmask%a, i, s_sky, mask)
      !end if
      !if (do_oper(prep_relbp)) then
      !   do j = 2, ndelta
      !      call project_sky(self, map_sky(:,:,:,j), pix, psi, flag, &
      !           & sprocmask2%a, i, s_sky_prop(:,:,j), mask2, s_bp=s_bp_prop(:,:,j))
      !   end do
      !else if (do_oper(prep_absbp)) then
      !   do j = 2, ndelta
      !      call project_sky(self, map_sky(:,:,:,j), pix, psi, flag, &
      !           & sprocmask2%a, i, s_sky_prop(:,:,j), mask2)
      !   end do
      !end if
      if (main_iter == 1 .and. self%first_call) then
         do j = 1, ndet
            if (all(mask(:,j) == 0)) self%scans(i)%d(j)%accept = .false.
            if (self%scans(i)%d(j)%sigma0 <= 0.d0) self%scans(i)%d(j)%accept = .false.
         end do
      end if
      
      ! Compute chisquare
      !TODO: this but probably different code?
      !if (do_oper(calc_chisq)) then
      !   call wall_time(t1)
      !   do j = 1, ndet
      !      if (.not. self%scans(i)%d(j)%accept) cycle
      !      s_buf(:,j) =  s_sl(:,j) + s_orb(:,j)
      !      if (do_oper(samp_mono)) s_buf(:,j) =  s_buf(:,j) + s_mono(:,j)
      !      call self%compute_chisq(i, j, mask(:,j), s_sky(:,j), &
      !           & s_buf(:,j), n_corr(:,j))
      !   end do
      !   call wall_time(t2); t_tot(7) = t_tot(7) + t2-t1
      !end if

      ! Select data
     call wall_time(t1)
     !TODO: cut scans here
     !if (self%first_scan) then 
         !do j = 1, ndet
         !   ntot= ntot + 1
         !   if (.not. self%scans(i)%d(j)%accept) cycle
         !   if (count(iand(flag(:,j),self%flag0) .ne. 0) > 0.1*ntod) then
         !      self%scans(i)%d(j)%accept = .false.
         !   else if (abs(self%scans(i)%d(j)%chisq) > chisq_threshold .or. &
         !   & isNaN(self%scans(i)%d(j)%chisq)) then
         !      write(*,fmt='(a,i8,i5,a,f12.1)') 'Reject scan, det = ', &
         !           & self%scanid(i), j, ', chisq = ', &
         !           & self%scans(i)%d(j)%chisq
         !      self%scans(i)%d(j)%accept = .false.
         !      cycle
         !   else
         !      naccept = naccept + 1
         !   end if
         !end do
     !end if
     
     if (any(.not. self%scans(i)%d%accept)) self%scans(i)%d%accept = .false.
     do j = 1, ndet
        if (.not. self%scans(i)%d(j)%accept) then
           self%scans(i)%d(self%partner(j))%accept = .false.
        end if
     end do
     call wall_time(t2); t_tot(15) = t_tot(15) + t2-t1

      ! Compute binned map
     allocate(d_calib(nout,ntod, ndet))
     do j = 1, ndet
       d_calib(1,:,j) = self%scans(i)%d(j)%tod
     !   if (.not. self%scans(i)%d(j)%accept) cycle
     !   inv_gain = 1.0 / real(self%scans(i)%d(j)%gain,sp)
     !   d_calib(1,:,j) = (self%scans(i)%d(j)%tod - n_corr(:,j)) * &
     !        & inv_gain - s_tot(:,j) + s_sky(:,j) - s_bp(:,j)
     !   if (do_oper(bin_map) .and. nout > 1) d_calib(2,:,j) = d_calib(1,:,j) - s_sky(:,j) + s_bp(:,j) ! Residual
     !   if (do_oper(bin_map) .and. nout > 2) d_calib(3,:,j) = (n_corr(:,j) - sum(n_corr(:,j)/ntod)) * inv_gain
     !   if (do_oper(bin_map) .and. nout > 3) d_calib(4,:,j) = s_bp(:,j)
     !   if (do_oper(bin_map) .and. do_oper(samp_mono) .and. nout > 4) d_calib(5,:,j) = s_mono(:,j)
     !   if (do_oper(bin_map) .and. nout > 5) d_calib(6,:,j) = s_orb(:,j)
     !   if (do_oper(bin_map) .and. nout > 6) d_calib(7,:,j) = s_sl(:,j)
     !   if (do_oper(bin_map) .and. nout > 7 .and. do_oper(sub_zodi)) d_calib(8,:,j) = s_zodi(:,j)

     !   end if
     end do
    
     !TODO: this routine almost certainly needs a wmap equivalent 
     call bin_TOD(self, d_calib, pix, psi, flag, A_map, b_map, i, .false.)
 

      ! Clean up
      call wall_time(t1)
      deallocate(n_corr, s_sl, s_sky, s_orb, s_tot, s_buf)
      deallocate(s_bp, s_sky_prop, s_bp_prop)
      deallocate(mask, mask2, pix, psi, flag)
      if (allocated(s_lowres)) deallocate(s_lowres)
      if (allocated(s_invN)) deallocate(s_invN)
      call wall_time(t2); t_tot(18) = t_tot(18) + t2-t1

      call wall_time(t8); t_tot(19) = t_tot(19) + t8-t7

   end do


   ! Output latest scan list with new timing information
   if (self%first_call) then
     call update_status(status, "scanlist1")
     call wall_time(t1)
     call self%output_scan_list(slist)
     !TODO: all this timing stuff is a disaster
     call wall_time(t2); t_tot(20) = t_tot(20) + t2-t1
     call update_status(status, "scanlist2")
   end if

   ! Solve combined map, summed over all pixels
   !TODO: this probably also needs changing, also why is this so long and not a
   !functions
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
     end do
     call mpi_win_fence(0, sA_map%win, ierr)
     call mpi_win_fence(0, sb_map%win, ierr)
     Sfilename = trim(prefix) // 'Smap'// trim(postfix)
     call finalize_binned_map(self, handle, sA_map, sb_map, rms_out, outmaps=outmaps)
    end if

   ! Output maps to disk
   call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
   call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))

   if (self%output_n_maps > 1) call outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
   !if (self%output_n_maps > 2) call outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
   !if (self%output_n_maps > 3) call outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix))
   !if (self%output_n_maps > 4) call outmaps(5)%p%writeFITS(trim(prefix)//'mono'//trim(postfix))
   !if (self%output_n_maps > 5) call outmaps(6)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
   !if (self%output_n_maps > 6) call outmaps(7)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
   !if (self%output_n_maps > 7) call outmaps(8)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))


  if (self%first_call) then
   call mpi_reduce(ntot,    i, 1, MPI_INTEGER, MPI_SUM, &
        & self%numprocs/2, self%info%comm, ierr)
   ntot = i
   call mpi_reduce(naccept, i, 1, MPI_INTEGER, MPI_SUM, &
        & self%numprocs/2, self%info%comm, ierr)
   naccept = i
  end if
  call wall_time(t2); t_tot(10) = t_tot(10) + t2-t1

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

    call int2string(iter, ctext)
    call update_status(status, "tod_end"//ctext)

    ! Parameter to check if this is first time routine has been
    self%first_call = .false.

  end subroutine process_WMAP_tod

end module comm_tod_WMAP_mod
