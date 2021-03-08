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

   integer(i4b), parameter :: N_test = 22
   integer(i4b), parameter :: samp_N = 1
   integer(i4b), parameter :: prep_G = 15
   integer(i4b), parameter :: samp_G = 2
   integer(i4b), parameter :: prep_acal = 3
   integer(i4b), parameter :: samp_acal = 4
   integer(i4b), parameter :: prep_rcal = 18
   integer(i4b), parameter :: samp_rcal = 19
   integer(i4b), parameter :: prep_relbp = 5
   integer(i4b), parameter :: prep_absbp = 16
   integer(i4b), parameter :: samp_bp = 11
   integer(i4b), parameter :: samp_sl = 6
   integer(i4b), parameter :: samp_N_par = 7
   integer(i4b), parameter :: sel_data = 8
   integer(i4b), parameter :: bin_map = 9
   integer(i4b), parameter :: calc_chisq = 10
   integer(i4b), parameter :: output_slist = 12
   integer(i4b), parameter :: samp_mono = 13
   integer(i4b), parameter :: sub_sl      = 14
   integer(i4b), parameter :: sub_zodi = 17
   integer(i4b), parameter :: sim_map = 20
   integer(i4b), parameter :: samp_imbal = 21
   integer(i4b), parameter :: samp_bline = 22
   logical(lgt), dimension(N_test) :: do_oper

   type, extends(comm_tod) :: comm_WMAP_tod
      class(orbdipole_pointer), allocatable :: orb_dp ! orbital dipole calculator
      character(len=20), allocatable, dimension(:) :: labels ! names of fields
   contains
      procedure     :: process_tod => process_WMAP_tod
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
      type(comm_params),      intent(in) :: cpar
      integer(i4b),           intent(in) :: id_abs
      class(comm_mapinfo),    target     :: info
      character(len=128),     intent(in) :: tod_type
      class(comm_WMAP_tod),   pointer    :: constructor

      integer(i4b) :: i, nside_beam, lmax_beam, nmaps_beam, ndelta
      character(len=512) :: datadir
      logical(lgt) :: pol_beam

      ! Set up WMAP specific parameters
      allocate (constructor)
      constructor%output_n_maps = 3
      constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
      constructor%nhorn = 2
      constructor%first_call = .true.
      constructor%verbosity = cpar%verbosity
      constructor%compressed_tod = .true.

      ! Iniitialize TOD labels
      allocate (constructor%labels(6))
      constructor%labels(1) = 'map'
      constructor%labels(2) = 'res'
      constructor%labels(3) = 'ncorr'
      constructor%labels(4) = 'orb_dp'
      constructor%labels(5) = 'sl'
      constructor%labels(6) = 'bpcorr'

      ! Initialize beams
      nside_beam = 512
      nmaps_beam = 3
      pol_beam = .true.
      constructor%nside_beam = nside_beam



      !initialize the common tod stuff
      call constructor%tod_constructor(cpar, id_abs, info, tod_type)


      !TODO: this is LFI specific, write something here for wmap
      call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label)



      ! Read the actual TOD
      call constructor%read_tod(constructor%label)

      call constructor%precompute_lookups()

      datadir = trim(cpar%datadir)//'/'

      ! Initialize bandpass mean and proposal matrix
      call constructor%initialize_bp_covar(trim(datadir)//cpar%ds_tod_bp_init(id_abs))

      !load the instrument file
      call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

      allocate(constructor%slconv(constructor%ndet))
      allocate (constructor%orb_dp)
      ! Need precompute the main beam precomputation for both the A-horn and
      ! B-horn.
      constructor%orb_dp%p => comm_orbdipole(constructor, constructor%mbeam)
   end function constructor

   !**************************************************
   !             Driver routine
   !**************************************************
   subroutine process_WMAP_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
      implicit none
      class(comm_WMAP_tod), intent(inout) :: self
      character(len=*), intent(in)    :: chaindir
      integer(i4b), intent(in)    :: chain, iter
      type(planck_rng), intent(inout) :: handle
      type(map_ptr), dimension(1:, 1:), intent(inout) :: map_in       ! (ndet,ndelta)
      real(dp), dimension(0:, 1:, 1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
      class(comm_map), intent(inout) :: map_out      ! Combined output map
      class(comm_map), intent(inout) :: rms_out      ! Combined output rms

      integer(i4b) :: i, j, k, l, m, n, t, ntod, ndet
      integer(i4b) :: nside, npix, nmaps, naccept, ntot, ext(2), nscan_tot, nhorn
      integer(i4b) :: ierr, main_iter, n_main_iter, ndelta, ncol, n_A, np0, nout
      character(len=40), allocatable, dimension(:) :: iter_labels ! names of iters 
      real(dp)     :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, chisq_threshold
      real(dp)     :: t_tot(22)
      real(sp)     :: inv_gain
      real(dp)     :: x_im(4)
      real(sp), allocatable, dimension(:, :)          :: n_corr, s_sky, s_skyA, s_skyB
      real(sp), allocatable, dimension(:, :)          :: s_sl, s_slA, s_slB
      real(sp), allocatable, dimension(:, :)          :: s_orbA, s_orbB, s_orb_tot
      real(sp), allocatable, dimension(:, :)          :: mask, mask2, mask_lowres, s_bp, s_bpA, s_bpB
      real(sp), allocatable, dimension(:, :)          :: s_mono, s_buf, s_tot, s_zodi
      real(sp), allocatable, dimension(:, :)          :: s_totA, s_totB
      real(sp), allocatable, dimension(:, :)          :: s_invN, s_lowres
      real(sp), allocatable, dimension(:, :, :)       :: s_sky_prop, s_bp_prop
      real(sp), allocatable, dimension(:, :, :)       :: d_calib
      real(dp), allocatable, dimension(:)             :: A_abscal, b_abscal
      real(dp), allocatable, dimension(:, :)          :: chisq_S, m_buf
      real(dp), allocatable, dimension(:, :)          :: A_map, dipole_mod, M_diag
      real(dp), allocatable, dimension(:, :, :)       :: b_map, b_mono, sys_mono
      integer(i4b), allocatable, dimension(:, :)      :: pix, psi
      integer(i4b), allocatable, dimension(:)         :: flag
      integer(i4b), allocatable, dimension(:,:)           :: tod
      real(dp), allocatable, dimension(:, :, :)       :: b_tot, M_diag_tot
      real(dp), allocatable, dimension(:, :)          :: cg_tot
      logical(lgt)       :: correct_sl, verbose, finished
      character(len=512) :: prefix, postfix, prefix4D, filename
      character(len=2048) :: Sfilename
      character(len=4)   :: ctext, myid_text
      character(len=6)   :: samptext, scantext
      character(len=512), allocatable, dimension(:) :: slist
      integer(i4b), allocatable, dimension(:)     :: procmask
      real(sp), allocatable, dimension(:, :, :, :) :: map_sky
      class(comm_map), pointer :: condmap
      class(map_ptr), allocatable, dimension(:) :: outmaps

      ! conjugate gradient parameters
      integer(i4b) :: i_max, i_min, num_cg_iters
      real(dp) :: delta_0, delta_old, delta_new, epsil(6)
      real(dp) :: alpha, beta, g, f_quad, sigma_mono
      real(dp), allocatable, dimension(:, :, :) :: cg_sol
      real(dp), allocatable, dimension(:, :)    :: r, s, d, q
      real(dp), allocatable, dimension(:)       :: map_full
      real(dp) :: monopole
      logical(lgt) :: write_cg_iter=.false.


      real(dp) :: phi, theta
      real(dp), dimension(3) :: vnorm
      integer(i4b) :: pixind


      real(dp) :: masked_var




      ! biconjugate gradient parameters
      real(dp) :: rho_old, rho_new
      real(dp) :: omega, delta_r, delta_s
      real(dp), allocatable, dimension(:, :) :: rhat, r0, shat, p, phat, v
      real(dp), allocatable, dimension(:)    :: determ

      call int2string(iter, ctext)
      call update_status(status, "tod_start"//ctext)

      t_tot = 0.d0
      call wall_time(t5)

      ! Set up full-sky map structures
      call wall_time(t1)
      correct_sl = .false.
      chisq_threshold = 6d0
      n_main_iter     = 6
      ndet = self%ndet
      nhorn = self%nhorn
      ndelta = size(delta, 3)
      nside = map_out%info%nside
      nmaps = map_out%info%nmaps
      npix = 12*nside**2
      nscan_tot = self%nscan_tot
      x_im(1) = self%x_im(1)
      x_im(2) = self%x_im(1)
      x_im(3) = self%x_im(2)
      x_im(4) = self%x_im(2)
      allocate(A_abscal(self%ndet), b_abscal(self%ndet))
      allocate (map_sky(nmaps, self%nobs, 0:ndet, ndelta))
      allocate (chisq_S(ndet, ndelta))
      allocate(dipole_mod(nscan_tot, ndet))
      allocate (slist(self%nscan))
      slist = ''
      allocate (outmaps(self%output_n_maps))
      do i = 1, self%output_n_maps
         outmaps(i)%p => comm_map(map_out%info)
      end do
      call int2string(chain, ctext)
      call int2string(iter, samptext)
      prefix = trim(chaindir)//'/tod_'//trim(self%freq)//'_'
      postfix = '_c'//ctext//'_k'//samptext//'.fits'

      ! Distribute fullsky maps
      allocate (m_buf(0:npix - 1, nmaps))
      !if (self%myid == 0) then
         allocate(map_full(0:npix-1))
         map_full = 0.d0
      !end if
      do j = 1, ndelta
         do i = 1, self%ndet
            map_in(i, j)%p%map = map_in(i, j)%p%map
            call map_in(i, j)%p%bcast_fullsky_map(m_buf)
            if (.true. .and. self%myid == 0) then
               filename = trim(chaindir) // '/fg_' // trim(self%label(i)) // '_v1.fits'
               call write_map2(filename, m_buf)
            end if
            do k = 1, self%nobs
               map_sky(:, k, i, j) = m_buf(self%ind2pix(k), :)
            end do
            if (self%myid == 0 .and. j == 1) map_full = map_full + m_buf(:,1)
         end do
         do k = 1, self%nobs
            do l = 1, nmaps
               map_sky(l, k, 0, j) = sum(map_sky(l, k, 1:ndet, j))/ndet
            end do
         end do
      end do
      if (self%myid == 0) map_full = map_full/self%ndet
      deallocate (m_buf)
!      call mpi_finalize(ierr)
!      stop

      allocate(procmask(0:npix-1))
      procmask = 0
      do i = 1, size(self%procmask%map(:,1))
         procmask(self%procmask%info%pix(i)) = nint(self%procmask%map(i-1,1))
      end do
      call mpi_allreduce(mpi_in_place, procmask, size(procmask), MPI_INTEGER, MPI_SUM, self%info%comm, ierr)
      !where (procmask .ge. 1)
      !    procmask = 1
      !end where

      call wall_time(t2); t_tot(9) = t2 - t1

      ! Compute far sidelobe Conviqt structures
      call wall_time(t1)
      do i = 1, self%ndet
         if (.not. correct_sl) exit
         !TODO: figure out why this is rotated
         call map_in(i,1)%p%YtW()  ! Compute sky a_lms
         self%slconv(i)%p => comm_conviqt(self%myid, self%comm_shared, &
              & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
              & 100, 3, 100, self%slbeam(i)%p, map_in(i,1)%p, 2)
      end do
      call wall_time(t2); t_tot(13) = t2-t1

      call update_status(status, "tod_init")
      call wall_time(t3)
      do_oper             = .true.
      ! There are five main iterations, for imbalance, absolute calibration, 
      ! relative calibration, time-variable calibration, and 
      ! correlated noise estimation.
      allocate(iter_labels(n_main_iter))
      iter_labels(1) = 'baseline'
      iter_labels(2) = 'absolute calibration'
      iter_labels(3) = 'relative calibration'
      iter_labels(4) = 'time-varying gain'
      iter_labels(5) = 'imbalance'
      iter_labels(6) = 'ncorr, npar, binmap, and chisq'

      main_it: do main_iter = 1, n_main_iter
         call update_status(status, "tod_istart")

         if (self%myid == 0 .and. self%verbosity > 0) write(*,'(A, I, X, A)') '  Performing main iteration = ', main_iter, iter_labels(main_iter)
         ! Select operations for current iteration
         do_oper(samp_bline)   = (main_iter == n_main_iter-5) ! .false. !  
         do_oper(samp_acal)    = (main_iter == n_main_iter-4) ! .false. !      
         do_oper(samp_rcal)    = (main_iter == n_main_iter-3) ! .false. !      
         do_oper(samp_G)       = (main_iter == n_main_iter-2) ! .false. !      
         do_oper(samp_imbal)   = (main_iter == n_main_iter-1) ! .false. !  
         do_oper(samp_N)       = (main_iter >= n_main_iter-0) ! .false. ! 
         do_oper(samp_N_par)   = do_oper(samp_N)
         do_oper(prep_relbp)   = ndelta > 1 .and. (main_iter == n_main_iter-0)
         do_oper(prep_absbp)   = .false. ! ndelta > 1 .and. (main_iter == n_main_iter-0) .and. .not. self%first_call .and. mod(iter,2) == 1
         do_oper(samp_bp)      = ndelta > 1 .and. (main_iter == n_main_iter-0)
         do_oper(samp_mono)    = .false.
         do_oper(bin_map)      = (main_iter == n_main_iter  )
         do_oper(sel_data)     = .false.
         do_oper(calc_chisq)   = (main_iter == n_main_iter  )
         do_oper(sub_sl)       = correct_sl
         do_oper(sub_zodi)     = self%subtract_zodi
         do_oper(output_slist) = mod(iter, 3) == 0
         do_oper(sim_map)      = .false. ! (main_iter == 1) !   

         dipole_mod = 0
         if (do_oper(bin_map) .or. do_oper(prep_relbp)) then
            if (do_oper(prep_relbp)) then
               nout = self%output_n_maps + ndelta - 1
            else
               nout = self%output_n_maps
            end if
         end if
         if (do_oper(bin_map)) then
           allocate (M_diag(0:npix-1, nmaps+1))
           allocate ( b_map(0:npix-1, nmaps, nout))
           M_diag = 0d0
           b_map = 0d0
         end if

         if (do_oper(samp_acal) .or. do_oper(samp_rcal) .or. do_oper(samp_imbal)) then
            A_abscal = 0.d0; b_abscal = 0.d0
         end if

         ! Perform main analysis loop
         naccept = 0; ntot = 0
         do i = 1, self%nscan
            call wall_time(t7)


            if (.not. any(self%scans(i)%d%accept)) cycle

            ! Short-cuts to local variables
            ndet = self%ndet
            ntod = self%scans(i)%ntod

            ! Set up local data structure for current scan
            allocate (n_corr(ntod, ndet))                 ! Correlated noise in V
            allocate (s_sl(ntod, ndet))                   ! Sidelobe in uKcmb 
            allocate (s_slA(ntod, ndet))                  ! Sidelobe in uKcmb (beam A)
            allocate (s_slB(ntod, ndet))                  ! Sidelobe in uKcmb (beam B)
            allocate (s_sky(ntod, ndet))                  ! Sky signal in uKcmb
            allocate (s_skyA(ntod, ndet))                  ! Sky signal in uKcmb
            allocate (s_skyB(ntod, ndet))                  ! Sky signal in uKcmb
            allocate (s_sky_prop(ntod, ndet, 2:ndelta))   ! Sky signal in uKcmb
            allocate (s_bp(ntod, ndet))                   ! Signal minus mean
            allocate (s_bpA(ntod, ndet))                   ! Signal minus mean
            allocate (s_bpB(ntod, ndet))                   ! Signal minus mean
            allocate (s_bp_prop(ntod, ndet, 2:ndelta))    ! Signal minus mean
            allocate (s_orbA(ntod, ndet))                 ! Orbital dipole (beam A)
            allocate (s_orbB(ntod, ndet))                 ! Orbital dipole (beam B)
            allocate (s_orb_tot(ntod, ndet))              ! Orbital dipole (both)
            allocate (s_buf(ntod, ndet))                  ! Buffer
            allocate (s_tot(ntod, ndet))                  ! Sum of all sky components
            allocate (s_totA(ntod, ndet))                  ! Sum of all sky components
            allocate (s_totB(ntod, ndet))                  ! Sum of all sky components
            allocate (mask(ntod, ndet))                   ! Processing mask in time
            allocate (mask2(ntod, ndet))                  ! Processing mask in time
            allocate (pix(ntod, nhorn))             ! Decompressed pointing
            allocate (psi(ntod, nhorn))             ! Decompressed pol angle
            allocate (tod(ntod, ndet))              ! Decompressed tod
            allocate (flag(ntod))                   ! Decompressed flags

            ! --------------------
            ! Analyze current scan
            ! --------------------

            ! Decompress pointing, psi and flags for current scan
            call wall_time(t1)
            call self%decompress_pointing_and_flags(i, 1, pix, psi, flag)
            do j = 1, ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               call self%decompress_tod(i, j, tod(:,j))
            end do
            ! Implement the tod 
            ! Every module that requires the tod, make the tod an argument
            call wall_time(t2); t_tot(11) = t_tot(11) + t2 - t1

            ! Construct sky signal template
            call wall_time(t1)
            if (do_oper(bin_map) .or. do_oper(prep_relbp)) then
               call project_sky_differential(self, map_sky(:,:,:,1), pix, psi, flag, &
                 & procmask, i, s_skyA, s_skyB, mask, s_bpA=s_bpA, s_bpB=s_bpB)
            else
               call project_sky_differential(self, map_sky(:,:,:,1), pix, psi, flag, &
                    & procmask, i, s_skyA, s_skyB, mask)
               s_bpA = 0.
               s_bpB = 0.
               s_bp  = 0.
            end if
            do j = 1, ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               s_sky(:, j) = (1d0+x_im(j))*s_skyA(:,j) - &
                           & (1d0-x_im(j))*s_skyB(:,j)
               if (do_oper(bin_map) .or. do_oper(prep_relbp)) then
                  s_bp(:, j)  = (1d0+x_im(j))*s_bpA(:,j) - &
                              & (1d0-x_im(j))*s_bpB(:,j)
               end if
            end do

            if (main_iter == 1 .and. self%first_call) then
               do j = 1, ndet
                  if (all(mask(:,j) == 0)) self%scans(i)%d(j)%accept = .false.
                  if (self%scans(i)%d(j)%sigma0 <= 0.d0) self%scans(i)%d(:)%accept = .false.
               end do
            end if
            call wall_time(t2); t_tot(1) = t_tot(1) + t2-t1

            ! Construct orbital dipole template
            call wall_time(t1)
            !call self%orb_dp%p%compute_orbital_dipole_pencil(i, pix(:,1), psi(:,1), s_orbA, 1d3)
            !call self%orb_dp%p%compute_orbital_dipole_pencil(i, pix(:,2), psi(:,2), s_orbB, 1d3)
            ! The ordering is a bit messed up, where 13/14 are beam A, 23/24 are beam B
            call self%orb_dp%p%compute_orbital_dipole_4pi_diff(i, pix(:,1), psi(:,1),s_orbA,1,1d3)
            call self%orb_dp%p%compute_orbital_dipole_4pi_diff(i, pix(:,2), psi(:,2),s_orbB,3,1d3)
            do j = 1, ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               s_orb_tot(:, j) = (1d0+x_im(j))*s_orbA(:,j) - &
                               & (1d0-x_im(j))*s_orbB(:,j)
            end do
            call wall_time(t2); t_tot(2) = t_tot(2) + t2-t1

            ! Construct sidelobe template
            call wall_time(t1)
            if (do_oper(sub_sl)) then
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) then
                    s_sl(:,j) = 0.
                    s_slA(:,j) = 0.
                    s_slB(:,j) = 0.
                    cycle
                  end if
                  ! K113/114 are horn A
                  ! K123/124 are horn B
                  ! Only for the sidelobe
                  call self%construct_sl_template(self%slconv(1)%p, &
                       & pix(:,1), psi(:,1), s_slA(:,j), 0d0)
                  call self%construct_sl_template(self%slconv(3)%p, &
                       & pix(:,2), psi(:,2), s_slB(:,j), 0d0)
                  s_slA(:,j) = 2 * s_slA(:,j) 
                  s_slB(:,j) = 2 * s_slB(:,j) 
                  s_sl(:, j) = (1d0+x_im(j))*s_slA(:,j) - &
                               (1d0-x_im(j))*s_slB(:,j)
               end do
            else
               s_sl = 0.
               s_slA = 0.
               s_slB = 0.
            end if
            call wall_time(t2); t_tot(12) = t_tot(12) + t2-t1

            ! Add orbital dipole and sidelobes to total signal
            s_totA = 0.
            s_totB = 0.
            s_tot = 0.
            do j = 1, ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               s_totA(:, j) = s_skyA(:, j) + s_slA(:, j) + s_orbA(:,j) 
               s_totB(:, j) = s_skyB(:, j) + s_slB(:, j) + s_orbB(:,j) 
               s_tot(:, j) = (1d0+x_im(j))*s_totA(:,j) - &
                           & (1d0-x_im(j))*s_totB(:,j)
            end do


            !!!!!!!!!!!!!!!!!!!!!
            ! Baseline estimation
            !!!!!!!!!!!!!!!!!!!!!

            ! n.b. this should be taken care of in the ncorr sampling. Still
            ! working on why this is the case. For now, just estimating it using
            ! Hinshaw et al. 2003 Eqn (21). If we are stuck with this, should
            ! replace it with a sampling step.

            do j = 1, ndet
              if (.not. self%scans(i)%d(j)%accept) cycle
              if (do_oper(samp_bline)) then
                self%scans(i)%d(j)%baseline =sum((tod(:,j) &
              &- self%scans(i)%d(j)%gain*s_tot(:,j))*mask(:,j))/sum(mask(:,j))
                if (trim(self%operation) == 'sample') then
                  self%scans(i)%d(j)%baseline = self%scans(i)%d(j)%baseline &
                   &  + rand_gauss(handle)/sqrt(sum(mask(:,j)*self%scans(i)%d(j)%sigma0**2))
                end if
              end if
            end do

            !!!!!!!!!!!!!!!!!!!
            ! Gain calculations
            !!!!!!!!!!!!!!!!!!!

            s_buf = 0.
            ! Precompute filtered signal for calibration
            if (do_oper(samp_G) .or. do_oper(samp_rcal) .or. do_oper(samp_acal)) then
               call self%downsample_tod(s_orb_tot(:,1), ext)
               allocate(     s_invN(ext(1):ext(2), ndet))      ! s * invN
               allocate(mask_lowres(ext(1):ext(2), ndet))
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  call self%downsample_tod(mask(:,j), ext, &
                       & mask_lowres(:,j), threshold=0.9)!, mask(:,j))
                  if (do_oper(samp_G) .or. do_oper(samp_rcal) .or. .not. self%orb_abscal) then
                     s_buf(:,j) = s_tot(:,j)
                     call fill_all_masked(s_buf(:,j), mask(:,j), ntod, &
                     &  .false., &   !trim(self%operation)=='sample', &
                     &  real(self%scans(i)%d(j)%sigma0, sp), &
                     &  handle, self%scans(i)%chunk_num)
                     call self%downsample_tod(s_buf(:,j), ext, &
                          & s_invN(:,j))!, mask(:,j))
                  else
                     call self%downsample_tod(s_orb_tot(:,j), ext, &
                          & s_invN(:,j))!, mask(:,j))
                  end if
               end do
               call multiply_inv_N(self, i, s_invN,   sampfreq=self%samprate_lowres, pow=0.5d0)

!!$                  write(*,*) 'hei',  self%myid
!!$                  if (self%myid == 62) then
!!$                     open(58,file='s_invN.dat')
!!$                     do k = 1, size(s_invN,1)
!!$                        write(58,*) k, s_invN(k,1)
!!$                     end do
!!$                     close(58)
!!$                  end if
!!$                  call mpi_finalize(ierr)
!!$                  stop

            else if (do_oper(samp_imbal)) then
               call self%downsample_tod(s_orb_tot(:,1), ext)
               allocate(  s_invN(ext(1):ext(2), ndet))      ! s * invN
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                   s_buf(:,j) = self%scans(i)%d(j)%gain*(s_totA(:,j)+s_totB(:,j))
                   call fill_all_masked(s_buf(:,j), mask(:,j), ntod, &
                   &  .false., &   !trim(self%operation)=='sample', &
                   &  real(self%scans(i)%d(j)%sigma0, sp), &
                   &  handle, self%scans(i)%chunk_num)
                   call self%downsample_tod(s_buf(:,j), ext, &
                        & s_invN(:,j))!, mask(:,j))
               end do
               call multiply_inv_N(self, i, s_invN,   sampfreq=self%samprate_lowres, pow=0.5d0)
            end if

            ! Prepare for absolute calibration
            call wall_time(t1)
            if (do_oper(samp_acal) .or. do_oper(samp_rcal)) then
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  if (do_oper(samp_acal)) then
                     if (self%orb_abscal) then
                        s_buf(:, j) = real(self%gain0(0),sp) * (s_tot(:, j) - s_orb_tot(:, j)) + &
                             & real(self%gain0(j) + self%scans(i)%d(j)%dgain,sp) * s_tot(:, j)
                     else
                        !if (self%scanid(i)==2114) write(*,*) j, self%gain0(j), self%scans(i)%d(j)%dgain, mean(abs(1.d0*s_tot(:, j)))
                        s_buf(:, j) = real(self%gain0(j) + self%scans(i)%d(j)%dgain,sp) * s_tot(:, j)
                     end if
                  else
                     s_buf(:,j) = real(self%gain0(0) + self%scans(i)%d(j)%dgain,sp) * s_tot(:, j)
                  end if
               end do

               call accumulate_abscal(self, i, mask, s_buf, s_invN, s_invN, A_abscal, b_abscal, handle, do_oper(samp_acal), mask_lowres=mask_lowres, tod_arr=tod)
            else if (do_oper(samp_imbal)) then
               do j = 1, ndet
                 ! this is subtracted from the data, equivalent to the first
                 ! equality in Eqn 17 of Gjerlow et al. 2021
                 s_buf(:,j) = self%scans(i)%d(j)%gain * (s_totA(:,j) - s_totB(:,j))
               end do

               call accumulate_imbal_cal(self, i, mask, s_buf, s_invN, s_invN, A_abscal, b_abscal, handle, tod_arr=tod)
            end if
            call wall_time(t2); t_tot(14) = t_tot(14) + t2-t1

            ! Fit gain
            if (do_oper(samp_G)) then
               call wall_time(t1)
               call calculate_gain_mean_std_per_scan(self, i, s_invN, mask, s_invN, s_tot, handle, mask_lowres=mask_lowres, tod_arr=tod)
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
               call sample_n_corr(self, handle, i, mask, s_buf, n_corr, pix, dospike=.false., tod_arr=tod)
               do j = 1, ndet
                  n_corr(:,j) = n_corr(:, j) - sum(n_corr(:,j))/ size(n_corr,1)
               end do
               call wall_time(t2); t_tot(3) = t_tot(3) + t2-t1
            else
               n_corr = 0.
            end if

            ! Compute noise spectrum
            if (do_oper(samp_N_par)) then
               call wall_time(t1)
               call sample_noise_psd(self, handle, i, mask, s_tot, n_corr, tod_arr=tod)
               call wall_time(t2); t_tot(6) = t_tot(6) + t2-t1
            end if

            !! Compute chisquare
            verbose = (self%verbosity > 2)
            if (do_oper(calc_chisq)) then
               call wall_time(t1)
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  s_buf(:,j) =  s_sl(:,j) + s_orb_tot(:,j)
                  if (do_oper(samp_mono)) s_buf(:,j) =  s_buf(:,j) + s_mono(:,j)
                  call self%compute_chisq(i, j, mask(:,j), s_sky(:,j), &
                       & s_buf(:,j), n_corr(:,j), verbose=.false., tod_arr=tod)
               end do
               call wall_time(t2); t_tot(7) = t_tot(7) + t2-t1
            end if

            !*******************
            ! Compute binned map
            !*******************

            ! Get calibrated map
            if (do_oper(bin_map)) then
               call wall_time(t1)
               allocate (d_calib(nout, ntod, ndet))
               d_calib = 0.
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  inv_gain = 1.0/real(self%scans(i)%d(j)%gain, sp)
!!$                  write(*,*) 'a', j, inv_gain
!!$                  write(*,*) 'b', j, sum(abs(n_corr(:,j)))
!!$                  write(*,*) 'c', j, sum(abs(s_tot(:,j)))
!!$                  write(*,*) 'd', j, sum(abs(s_sky(:,j)))
!!$                  write(*,*) 'e', j, sum(abs(s_bp(:,j)))
                  d_calib(1, :, j) = (tod(:,j) - &
                     & self%scans(i)%d(j)%baseline - n_corr(:, j))* &
                     & inv_gain - s_tot(:, j) + s_sky(:, j) - s_bp(:, j)
                  if (nout > 1) d_calib(2, :, j) = d_calib(1, :, j) - &
                    & s_sky(:, j) + s_bp(:, j) ! Residual

                  if (nout > 2) d_calib(3, :, j) = (n_corr(:, j) - sum(n_corr(:, j)/ntod))*inv_gain
                  if (do_oper(bin_map) .and. nout > 3) d_calib(4,:,j) = s_orb_tot(:,j)
                  if (do_oper(bin_map) .and. nout > 4) d_calib(5,:,j) = s_sl(:,j)
                  if (do_oper(bin_map) .and. nout > 5) d_calib(6,:,j) = s_bp(:,j)

                  if (do_oper(prep_relbp)) then
                     do k = 2, ndelta
                        d_calib(self%output_n_maps+k-1,:,j) = d_calib(1,:,j) + s_bp(:,j) - s_bp_prop(:,j,k)
                     end do
                  end if

               end do

               call wall_time(t2); t_tot(5) = t_tot(5) + t2-t1


               call wall_time(t1)
               if (.false. .and. i==1 .and. do_oper(bin_map) .and. self%first_call) then
                  call int2string(self%scanid(i), scantext)
                  if (self%myid == 0 .and. self%verbosity > 0) write(*,*) 'Writing tod to txt'
                  do k = 1, self%ndet
                     open(78,file=trim(chaindir)//'/tod_'//trim(self%label(k))//'_pid'//scantext//'.dat', recl=1024)
                     write(78,*) "# Sample   uncal_TOD (mK)  n_corr (mK) cal_TOD (mK)  skyA (mK)  skyB (mK)"// &
                          & " s_orbA (mK)  s_orbB (mK)  mask, baseline, invgain, flag"
                     do j = 1, ntod
                        inv_gain = 1.0/real(self%scans(i)%d(k)%gain, sp)
                        write(78,*) j, tod(j, k), n_corr(j, k), d_calib(1,j,k), &
                         &  s_skyA(j,k), s_skyB(j,k), s_orbA(j,k), s_orbB(j,k), mask(j, k), self%scans(i)%d(k)%baseline, inv_gain, flag(j)
                     end do
                     close(78)
                  end do
                  !call mpi_finalize(ierr)
                  !stop
               end if
               call wall_time(t2); t_tot(22) = t2 - t1

               call wall_time(t1)
               ! Bin the calibrated map
               call bin_differential_TOD(self, d_calib, pix,  &
                 & psi, flag, self%x_im, procmask, b_map, M_diag, i, &
                 & do_oper(prep_relbp))
               deallocate(d_calib)
               call wall_time(t2); t_tot(8) = t_tot(8) + t2-t1
            end if

            do j = 1, ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               masked_var = masked_variance(s_sky(:, j), mask(:, j))
               if (masked_var == 9999999999999) then
                 dipole_mod(self%scanid(i), j) = 0
               else
                 dipole_mod(self%scanid(i), j) = masked_var
               end if
            end do

            ! Clean up
            deallocate (n_corr, s_sky, s_orbA, s_orbB, s_orb_tot, s_tot, s_buf)
            deallocate (s_totA, s_totB)
            deallocate (s_sl, s_slA, s_slB, s_sky_prop)
            deallocate (s_skyA, s_skyB)
            deallocate (mask, mask2, pix, psi, flag, tod)
            if (allocated(s_invN)) deallocate (s_invN)
            if (allocated(mask_lowres)) deallocate (mask_lowres)
            deallocate(s_bp, s_bp_prop, s_bpA, s_bpB)

            call wall_time(t8); t_tot(19) = t_tot(19) + t8-t7
            if (do_oper(output_slist)) then
               self%scans(i)%proctime   = self%scans(i)%proctime   + t8-t7
               self%scans(i)%n_proctime = self%scans(i)%n_proctime + 1
               if (main_iter == n_main_iter) then
                  write(slist(i),*) self%scanid(i), '"',trim(self%hdfname(i)), &
                       & '"', real(self%scans(i)%proctime/self%scans(i)%n_proctime,sp),real(self%spinaxis(i,:),sp)
               end if
            end if
         end do

         call mpi_allreduce(mpi_in_place, dipole_mod, size(dipole_mod), MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)

         if (do_oper(samp_imbal)) then
            call wall_time(t1)
            call sample_imbal_cal(self, handle, A_abscal, b_abscal)
            x_im(1) = self%x_im(1)
            x_im(2) = self%x_im(1)
            x_im(3) = self%x_im(2)
            x_im(4) = self%x_im(2)
            call wall_time(t2); t_tot(16) = t_tot(16) + t2-t1
         end if

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

      end do main_it
      deallocate(iter_labels)
      call wall_time(t4)


      ! Output latest scan list with new timing information
      if (do_oper(output_slist)) then
         call update_status(status, "scanlist1")
         call wall_time(t1)
         call self%output_scan_list(slist)
         call wall_time(t2); t_tot(20) = t_tot(20) + t2-t1
         call update_status(status, "scanlist2")
      end if

      ! ************************
      ! Finalize maps
      ! ************************
      call wall_time(t1)
      call update_status(status, "Running allreduce on M_diag")
      call mpi_allreduce(mpi_in_place, M_diag, size(M_diag), &
           & MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)
      call update_status(status, "Running allreduce on b")
      call mpi_allreduce(mpi_in_place, b_map, size(b_map), &
           & MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)


      np0 = self%info%np
      allocate (cg_tot(0:np0 - 1, nmaps))

      ! write out M_diag, b_map to fits.
      !cg_tot = b_map(self%info%pix, 1:nmaps, 1)
      !call write_map2(trim(prefix)//'b2'//trim(postfix), b_map(:,:,1))
      !call self%write_fits_file_iqu(trim(prefix)//'b'//trim(postfix), cg_tot, outmaps)
      !cg_tot = M_diag(self%info%pix, 1:nmaps)
      !call write_fits_file_iqu(trim(prefix)//'M'//trim(postfix), cg_tot, outmaps)
      !cg_tot = b_map(self%info%pix, 1:nmaps, 5)
      !call write_fits_file_iqu(trim(prefix)//'b_sl'//trim(postfix), cg_tot, outmaps)


      where (M_diag == 0d0)
         M_diag = 1d0
      end where

      ! Conjugate Gradient solution to (P^T Ninv P) m = P^T Ninv d, or Ax = b
      call update_status(status, "Allocating cg arrays")
      allocate (cg_sol(0:npix-1, nmaps, nout))
      if (self%myid == 0) then 
         cg_sol = 0.0d0
         epsil(1)   = 1d-10
         epsil(2:6) = 1d-6
         num_cg_iters = 0
         if (self%verbosity > 0) write(*,*) '  Running BiCG'
      end if

      call wall_time(t9)
      call update_status(status, "Starting bicg-stab")
      do l=1, self%output_n_maps
         if (self%verbosity > 0 .and. self%myid == 0) then
           write(*,*) '    Solving for ', trim(adjustl(self%labels(l)))
         end if
         call run_bicgstab(self, handle, cg_sol, npix, nmaps, num_cg_iters, &
                          & epsil(l), procmask, map_full, M_diag, b_map, l)
      end do

      call wall_time(t10); t_tot(21) = (t10 - t9)

      call mpi_bcast(cg_sol, size(cg_sol),  MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)
      call mpi_bcast(num_cg_iters, 1,  MPI_INTEGER, 0, self%info%comm, ierr)
      do k = 1, self%output_n_maps
         do j = 1, nmaps
            outmaps(k)%p%map(:, j) = cg_sol(self%info%pix, j, k)
         end do
      end do

      map_out%map = outmaps(1)%p%map
      rms_out%map = M_diag(self%info%pix, 1:nmaps)**-0.5
      call outmaps(1)%p%writeFITS(trim(prefix)//'map'//trim(postfix))
      call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
      do n = 2, self%output_n_maps
        call outmaps(n)%p%writeFITS(trim(prefix)//trim(adjustl(self%labels(n)))//trim(postfix))
      end do

      if (self%first_call) then
         call mpi_reduce(ntot, i, 1, MPI_INTEGER, MPI_SUM, &
              & self%numprocs/2, self%info%comm, ierr)
         ntot = i
         call mpi_reduce(naccept, i, 1, MPI_INTEGER, MPI_SUM, &
              & self%numprocs/2, self%info%comm, ierr)
         naccept = i
      end if
      call wall_time(t2); t_tot(10) = t_tot(10) + t2 - t1
      call wall_time(t6)
      if (self%myid == self%numprocs/2 .and. self%verbosity > 0) then
         write(*,*) '  Time dist sky        = ', nint(t_tot(9))
         write(*,*) '  Time sl precomp      = ', nint(t_tot(13))
         write(*,*) '  Time decompress      = ', nint(t_tot(11))
         write(*,*) '  Time project         = ', nint(t_tot(1))
         write(*,*) '  Time orbital         = ', nint(t_tot(2))
         write(*,*) '  Time sl interp       = ', nint(t_tot(12))
         write(*,*) '  Time ncorr           = ', nint(t_tot(3))
         write(*,*) '  Time gain            = ', nint(t_tot(4))
         write(*,*) '  Time absgain         = ', nint(t_tot(14))
         write(*,*) '  Time sel data        = ', nint(t_tot(15))
         write(*,*) '  Time clean           = ', nint(t_tot(5))
         write(*,*) '  Time noise           = ', nint(t_tot(6))
         write(*,*) '  Time samp abs        = ', nint(t_tot(16))
         write(*,*) '  Time samp bp         = ', nint(t_tot(17))
         write(*,*) '  Time chisq           = ', nint(t_tot(7))
         write(*,*) '  Time bin             = ', nint(t_tot(8))
         write(*,*) '  Time scanlist        = ', nint(t_tot(20))
         write(*,*) '  Time final           = ', nint(t_tot(10))
         write(*,*) '  Time solving cg      = ', nint(t_tot(21))
         write(*,*) '  Time per cg iter     = ', nint(t_tot(21)/num_cg_iters)
         write(*,*) '  Number of cg iters   = ', num_cg_iters
         write(*,*) '  Time writing to txt  = ', nint(t_tot(22))
         write(*,*) '  Time total           = ', int(t6-t5), int(sum(t_tot(1:18)))
      end if

      ! Clean up temporary arrays
      deallocate(A_abscal, b_abscal, chisq_S, procmask)
      deallocate(b_map, M_diag, cg_tot)
      !if (self%myid ==0) deallocate(map_full)
      deallocate(map_full)
      if (allocated(b_mono)) deallocate (b_mono)
      if (allocated(sys_mono)) deallocate (sys_mono)
      if (allocated(slist)) deallocate (slist)
      if (allocated(dipole_mod)) deallocate (dipole_mod)

      if (allocated(outmaps)) then
         do i = 1, self%output_n_maps
            call outmaps(i)%p%dealloc
         end do
         deallocate (outmaps)
      end if

      deallocate (map_sky, cg_sol)
      if (self%myid == 0) deallocate (r, rhat, s, r0, q, shat, p, phat, v, m_buf, determ)

      if (correct_sl) then
         do i = 1, self%ndet
            call self%slconv(i)%p%dealloc(); deallocate(self%slconv(i)%p)
         end do
      end if

      call int2string(iter, ctext)
      call update_status(status, "tod_end"//ctext)

      ! Parameter to check if this is first time routine has been
      self%first_call = .false.

   end subroutine process_WMAP_tod





  subroutine sample_imbal_cal(tod, handle, A_abs, b_abs)
    implicit none
    class(comm_WMAP_tod),              intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs

    integer(i4b) :: i, j, ierr
    real(dp), allocatable, dimension(:) :: A, b

    ! Collect contributions from all cores
    allocate(A(tod%ndet), b(tod%ndet))
    call mpi_reduce(A_abs, A, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)
    call mpi_reduce(b_abs, b, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)

    ! Compute gain update and distribute to all cores
    if (tod%myid == 0) then
       tod%x_im(1) = sum(b(1:2))/sum(A(1:2))
       tod%x_im(2) = sum(b(3:4))/sum(A(3:4))
       if (trim(tod%operation) == 'sample') then
          ! Add fluctuation term if requested
          tod%x_im(1) = tod%x_im(1) + 1.d0/sqrt(sum(A(1:2))) * rand_gauss(handle)
          tod%x_im(2) = tod%x_im(2) + 1.d0/sqrt(sum(A(3:4))) * rand_gauss(handle)
       end if
       if (tod%verbosity > 1) then
         write(*,*) 'imbal sample =', tod%x_im
         write(*,*) 'b', sum(b(1:2)), sum(b(3:4))
         write(*,*) 'A', sum(A(1:2)), sum(A(3:4))
       end if
    end if
    call mpi_bcast(tod%x_im, 2,  MPI_DOUBLE_PRECISION, 0, &
         & tod%info%comm, ierr)


    deallocate(A, b)

  end subroutine sample_imbal_cal

  subroutine accumulate_imbal_cal(tod, scan, mask, s_sub, s_ref, s_invN, A_abs, b_abs, handle, tod_arr)
    implicit none
    class(comm_tod),                   intent(in)     :: tod
    integer(i4b),                      intent(in)     :: scan
    real(sp),          dimension(:,:), intent(in)     :: mask, s_sub, s_ref
    real(sp),          dimension(:,:), intent(in)     :: s_invN
    real(dp),          dimension(:),   intent(inout)  :: A_abs, b_abs
    type(planck_rng),                  intent(inout)  :: handle
    integer(i4b),      dimension(:,:), intent(in), optional :: tod_arr

    real(sp), allocatable, dimension(:,:)     :: residual
    real(sp), allocatable, dimension(:)       :: r_fill
    real(dp)     :: A, b
    integer(i4b) :: i, j, ext(2), ndet, ntod
    character(len=5) :: itext

    ndet = tod%ndet
    ntod = size(s_sub,1)
    call tod%downsample_tod(s_sub(:,1), ext)    
    allocate(residual(ext(1):ext(2),ndet), r_fill(size(s_sub,1)))
    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) then
          residual(:,j) = 0.
          cycle
       end if
       r_fill = tod_arr(:,j) - s_sub(:,j) - tod%scans(scan)%d(j)%baseline
       call fill_all_masked(r_fill, mask(:,j), ntod, trim(tod%operation) == 'sample', abs(real(tod%scans(scan)%d(j)%sigma0, sp)), handle, tod%scans(scan)%chunk_num)
       call tod%downsample_tod(r_fill, ext, residual(:,j))
    end do

    call multiply_inv_N(tod, scan, residual, sampfreq=tod%samprate_lowres, pow=0.5d0)

    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       A_abs(j) = A_abs(j) + sum(s_invN(:,j) * s_ref(:,j))
       b_abs(j) = b_abs(j) + sum(s_invN(:,j) * residual(:,j))
    end do

    deallocate(residual, r_fill)

  end subroutine accumulate_imbal_cal

end module comm_tod_WMAP_mod
