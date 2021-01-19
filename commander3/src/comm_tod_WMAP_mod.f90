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

   integer(i4b), parameter :: N_test = 20
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
   logical(lgt), dimension(N_test) :: do_oper

   type, extends(comm_tod) :: comm_WMAP_tod
      class(orbdipole_pointer), allocatable :: orb_dp ! orbital dipole calculator
      real(dp), allocatable, dimension(:)  :: x_im    ! feedhorn imbalance parameters
      character(len=20), allocatable, dimension(:) :: labels ! names of fields
   contains
      procedure     :: process_tod => process_WMAP_tod
   end type comm_WMAP_tod

   interface comm_WMAP_tod
      procedure constructor
   end interface comm_WMAP_tod

contains

   !*************************************************
   !    Convert integer to string
   !*************************************************
   character(len=20) function str(k)
       integer, intent(in) :: k
       write (str, *) k
       str = adjustl(str)
   end function str


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
      constructor%output_n_maps = 5
      constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
      constructor%nhorn = 2


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
      allocate (constructor%x_im(constructor%ndet/2))
      constructor%x_im(:) = 0.0d0
      ! For K-band
      constructor%x_im = [-0.00067, 0.00536]
      ! constructor%x_im = [-0.05, 0.05]


      !TODO: this is LFI specific, write something here for wmap
      call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label)



      ! Read the actual TOD
      call constructor%read_tod_WMAP(constructor%label)

      call constructor%precompute_lookups()

      datadir = trim(cpar%datadir)//'/'

      ! Initialize bandpass mean and proposal matrix
      call constructor%initialize_bp_covar(trim(datadir)//cpar%ds_tod_bp_init(id_abs))

      !load the instrument file
      call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

      allocate(constructor%slconv(constructor%ndet))
      allocate (constructor%orb_dp)
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
      real(dp)     :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, chisq_threshold
      real(dp)     :: t_tot(22)
      real(sp)     :: inv_gain
      real(sp), allocatable, dimension(:, :)          :: n_corr, s_sky
      real(sp), allocatable, dimension(:, :)          :: s_sl, s_slA, s_slB
      real(sp), allocatable, dimension(:, :)          :: s_orbA, s_orbB, s_orb_tot
      real(sp), allocatable, dimension(:, :)          :: s_solA, s_solB, s_sol_tot
      real(sp), allocatable, dimension(:, :)          :: mask, mask2, s_bp
      real(sp), allocatable, dimension(:, :)          :: s_mono, s_buf, s_tot, s_zodi
      real(sp), allocatable, dimension(:, :)          :: s_invN, s_lowres
      real(sp), allocatable, dimension(:, :, :)       :: s_sky_prop, s_bp_prop
      real(sp), allocatable, dimension(:, :, :)       :: d_calib
      real(dp), allocatable, dimension(:)             :: A_abscal, b_abscal
      real(dp), allocatable, dimension(:, :)          :: chisq_S, m_buf
      real(dp), allocatable, dimension(:, :)          :: A_map, dipole_mod, M_diag
      real(dp), allocatable, dimension(:, :, :)       :: b_map, b_mono, sys_mono
      integer(i4b), allocatable, dimension(:, :, :)   :: pix, psi
      integer(i4b), allocatable, dimension(:, :)      :: flag
      real(dp), allocatable, dimension(:, :, :)       :: b_tot, M_diag_tot
      real(dp), allocatable, dimension(:, :)          :: cg_tot
      logical(lgt)       :: correct_sl
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
      integer(i4b) :: i_max=100
      real(dp) :: delta_0, delta_old, delta_new, epsil
      real(dp) :: alpha, beta, g, f_quad
      real(dp), allocatable, dimension(:, :, :) :: cg_sol
      real(dp), allocatable, dimension(:, :)    :: r, s, d, q
      logical(lgt) :: write_cg_iter=.false., verbose=.false.


      logical(lgt) :: remove_solar_dipole=.false.
      real(dp) :: phi, theta
      real(dp), dimension(3) :: vnorm
      integer(i4b) :: pixind


      real(dp) :: masked_var




      ! biconjugate gradient parameters
      real(dp) :: rho_old, rho_new
      real(dp) :: omega, delta_r, delta_s
      real(dp), allocatable, dimension(:, :) :: r0, shat, p, phat, v

      if (iter > 1) self%first_call = .false.
      call int2string(iter, ctext)
      call update_status(status, "tod_start"//ctext)

      t_tot = 0.d0
      call wall_time(t5)

      ! Set up full-sky map structures
      call wall_time(t1)
      correct_sl = .true.
      chisq_threshold = 6d0
      n_main_iter     = 4
      ndet = self%ndet
      nhorn = self%nhorn
      ndelta = size(delta, 3)
      nside = map_out%info%nside
      nmaps = map_out%info%nmaps
      npix = 12*nside**2
      nout = self%output_n_maps
      nscan_tot = self%nscan_tot
      allocate(A_abscal(self%ndet), b_abscal(self%ndet))
      allocate (map_sky(nmaps, self%nobs, 0:ndet, ndelta))
      allocate (chisq_S(ndet, ndelta))
      allocate(dipole_mod(nscan_tot, ndet))
      allocate (slist(self%nscan))
      slist = ''
      allocate (outmaps(nout))
      do i = 1, nout
         outmaps(i)%p => comm_map(map_out%info)
      end do
      call int2string(chain, ctext)
      call int2string(iter, samptext)
      prefix = trim(chaindir)//'/tod_'//trim(self%freq)//'_'
      postfix = '_c'//ctext//'_k'//samptext//'.fits'

      ! Distribute fullsky maps
      allocate (m_buf(0:npix - 1, nmaps))
      do j = 1, ndelta
         do i = 1, self%ndet
            map_in(i, j)%p%map = map_in(i, j)%p%map
            call map_in(i, j)%p%bcast_fullsky_map(m_buf)
            do k = 1, self%nobs
               map_sky(:, k, i, j) = m_buf(self%ind2pix(k), :)
            end do
         end do
         do k = 1, self%nobs
            do l = 1, nmaps
               map_sky(l, k, 0, j) = sum(map_sky(l, k, 1:ndet, j))/ndet
            end do
         end do
      end do
      deallocate (m_buf)

      allocate(procmask(0:npix-1))
      do i = 1, size(self%procmask%map(:,1))
         procmask(self%procmask%info%pix(i)) = nint(self%procmask%map(i-1,1))
      end do
      call mpi_allreduce(mpi_in_place, procmask, size(procmask), MPI_INTEGER, MPI_SUM, self%info%comm, ierr)
      where (procmask .ge. 1)
          procmask = 1
      end where

      call wall_time(t2); t_tot(9) = t2 - t1


      ! Compute far sidelobe Conviqt structures
      call wall_time(t1)
      do i = 1, self%ndet
         if (.not. correct_sl) exit

         !TODO: figure out why this is rotated
         call map_in(i,1)%p%YtW()  ! Compute sky a_lms
         self%slconv(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
              & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
              & 100, 3, 100, self%slbeam(i)%p, map_in(i,1)%p, 2)
!         write(*,*) i, 'b', sum(abs(self%slconv(i)%p%c%a))
      end do
      call wall_time(t2); t_tot(13) = t2-t1




      call update_status(status, "tod_init")
      do_oper             = .true.
      allocate (M_diag(0:npix-1, nmaps+1))
      allocate ( b_map(0:npix-1, nmaps+1, nout))
      M_diag = 0d0
      b_map = 0d0
      ! There are four main iterations, for absolute calibration, relative
      ! calibration, time-variable calibration, and correlated noise estimation.
      main_it: do main_iter = 1, n_main_iter
         call update_status(status, "tod_istart")

         if (self%myid == 0) write(*,*) '  Performing main iteration = ', main_iter
         ! Select operations for current iteration
         do_oper(samp_acal)    = (main_iter == n_main_iter-3) ! .false. !      
         do_oper(samp_rcal)    = (main_iter == n_main_iter-2) ! .false. !      
         do_oper(samp_G)       = (main_iter == n_main_iter-1) ! .false. !      
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
         do_oper(output_slist) = mod(iter, 1) == 0
         do_oper(sim_map)      = .false. ! (main_iter == 1) !   

         dipole_mod = 0

         if (do_oper(samp_acal) .or. do_oper(samp_rcal)) then
            A_abscal = 0.d0; b_abscal = 0.d0
         end if

         ! Perform main analysis loop
         naccept = 0; ntot = 0
         do i = 1, self%nscan


            if (.not. any(self%scans(i)%d%accept)) cycle

            ! Short-cuts to local variables
            call wall_time(t1)
            ndet = self%ndet
            ntod = self%scans(i)%ntod

            ! Set up local data structure for current scan
            allocate (n_corr(ntod, ndet))                 ! Correlated noise in V
            allocate (s_sl(ntod, ndet))                   ! Sidelobe in uKcmb 
            allocate (s_slA(ntod, ndet))                  ! Sidelobe in uKcmb (beam A)
            allocate (s_slB(ntod, ndet))                  ! Sidelobe in uKcmb (beam B)
            allocate (s_sky(ntod, ndet))                  ! Sky signal in uKcmb
            allocate (s_sky_prop(ntod, ndet, 2:ndelta))   ! Sky signal in uKcmb
            allocate (s_bp(ntod, ndet))                   ! Signal minus mean
            allocate (s_bp_prop(ntod, ndet, 2:ndelta))    ! Signal minus mean
            allocate (s_orbA(ntod, ndet))                 ! Orbital dipole (beam A)
            allocate (s_orbB(ntod, ndet))                 ! Orbital dipole (beam B)
            allocate (s_orb_tot(ntod, ndet))              ! Orbital dipole (both)
            allocate (s_solA(ntod, ndet))                 ! Solar dipole (beam A)
            allocate (s_solB(ntod, ndet))                 ! Solar dipole (beam B)
            allocate (s_sol_tot(ntod, ndet))              ! Solar dipole (both)
            allocate (s_buf(ntod, ndet))                  ! Buffer
            allocate (s_tot(ntod, ndet))                  ! Sum of all sky components
            allocate (mask(ntod, ndet))                   ! Processing mask in time
            allocate (mask2(ntod, ndet))                  ! Processing mask in time
            allocate (pix(ntod, ndet, nhorn))             ! Decompressed pointing
            allocate (psi(ntod, ndet, nhorn))             ! Decompressed pol angle
            allocate (flag(ntod, ndet))                   ! Decompressed flags

            ! --------------------
            ! Analyze current scan
            ! --------------------

            ! Decompress pointing, psi and flags for current scan
            call wall_time(t1)
            do j = 1, ndet
               call self%decompress_pointing_and_flags(i, j, pix(:, j, :), &
                    & psi(:, j, :), flag(:, j))
            end do
            call wall_time(t2); t_tot(11) = t_tot(11) + t2 - t1

            ! Construct sky signal template
            call wall_time(t1)
            if (do_oper(bin_map) .or. do_oper(prep_relbp)) then
               call project_sky_differential(self, map_sky(:,:,:,1), pix, psi, flag, &
                 & self%x_im, procmask, i, s_sky, mask, do_oper(sim_map), s_bp=s_bp)
            else
               call project_sky_differential(self, map_sky(:,:,:,1), pix, psi, flag, &
                    & self%x_im, procmask, i, s_sky, mask, do_oper(sim_map))
            end if
            if (do_oper(prep_relbp)) then
               do j = 2, ndelta
                  call project_sky_differential(self, map_sky(:,:,:,j), pix, psi, flag, &
                       & self%x_im, procmask, i, s_sky_prop(:,:,j), mask, do_oper(sim_map), s_bp=s_bp_prop(:,:,j))
               end do
            else if (do_oper(prep_absbp)) then
               do j = 2, ndelta
                  call project_sky_differential(self, map_sky(:,:,:,j), pix, psi, flag, &
                       & self%x_im, procmask, i, s_sky_prop(:,:,j), mask, do_oper(sim_map))
               end do
            end if


            if (main_iter == 1 .and. self%first_call) then
               do j = 1, ndet
                  self%scans(i)%d(j)%accept = .true.
                  !if (all(mask(:,j) == 0)) self%scans(i)%d(j)%accept = .false.
                  !if (self%scans(i)%d(j)%sigma0 <= 0.d0) self%scans(i)%d(j)%accept = .false.
               end do
            end if

            ! Construct orbital dipole template
            call wall_time(t1)
            s_orbA = 0d0
            s_orbB = 0d0
            call self%orb_dp%p%compute_orbital_dipole_pencil(i, pix(:,:,1), psi(:,:,1), s_orbA)
            call self%orb_dp%p%compute_orbital_dipole_pencil(i, pix(:,:,2), psi(:,:,2), s_orbB)
            !call self%orb_dp%p%compute_orbital_dipole_4pi(i, pix(:,:,1), psi(:,:,1), s_orbA)
            !call self%orb_dp%p%compute_orbital_dipole_4pi(i, pix(:,:,2), psi(:,:,2), s_orbB)
            s_orbA = s_orbA * 1d6 ! MK -> mK
            s_orbB = s_orbB * 1d6 ! MK -> mK
            s_solA = 0d0
            s_solB = 0d0
            call self%orb_dp%p%compute_solar_dipole_pencil(i, pix(:,:,1), psi(:,:,1), s_solA)
            call self%orb_dp%p%compute_solar_dipole_pencil(i, pix(:,:,2), psi(:,:,2), s_solB)
            !call self%orb_dp%p%compute_solar_dipole_4pi(i, pix(:,:,1), psi(:,:,1), s_solA)
            !call self%orb_dp%p%compute_solar_dipole_4pi(i, pix(:,:,2), psi(:,:,2), s_solB)
            s_solA = s_solA * 1d3 ! K -> mK
            s_solB = s_solB * 1d3 ! K -> mK
            do j = 1, ndet
               s_orb_tot(:, j) = (1+self%x_im((j+1)/2))*s_orbA(:,j) - &
                               & (1-self%x_im((j+1)/2))*s_orbB(:,j)
               s_sol_tot(:, j) = (1+self%x_im((j+1)/2))*s_solA(:,j) - &
                               & (1-self%x_im((j+1)/2))*s_solB(:,j)
            end do

            if (do_oper(sim_map)) then
                do j = 1, ndet
                   inv_gain = 1.0/real(self%scans(i)%d(j)%gain, sp)
                   self%scans(i)%d(j)%tod = floor(self%scans(i)%d(j)%tod + s_orb_tot(:,j)/inv_gain)
                   if (inv_gain < 0) then
                      self%scans(i)%d(j)%tod = -self%scans(i)%d(j)%tod
                      self%scans(i)%d(j)%gain = - self%scans(i)%d(j)%gain
                   end if
                end do
            end if
            call wall_time(t2); t_tot(2) = t_tot(2) + t2-t1


            ! Construct sidelobe template
            call wall_time(t1)
            if (do_oper(sub_sl)) then
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  call self%construct_sl_template(self%slconv(j)%p, &
                       & pix(:,j,1), psi(:,j,1), s_slA(:,j), 0d0)
                  call self%construct_sl_template(self%slconv(j)%p, &
                       & pix(:,j,2), psi(:,j,2), s_slB(:,j), 0d0)
                  s_sl(:,j) = (1+self%x_im((j+1)/2))*s_slA(:,j) - &
                            & (1-self%x_im((j+1)/2))*s_slB(:,j)
                  !s_sl(:,j) = 0.5 * s_sl(:,j) ! Scaling by a factor of 1/2, following Barnes notation
                  s_sl(:,j) = 2 * s_sl(:,j)   ! Scaling by a factor of 2, need to understand why
               end do
            else
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  s_sl(:,j) = 0.
               end do
            end if
            call wall_time(t2); t_tot(12) = t_tot(12) + t2-t1


            ! Add orbital dipole and sidelobes to total signal
            s_buf = 0.d0
            do j = 1, ndet
               s_tot(:, j) = s_sky(:, j) + s_sl(:, j) + s_orb_tot(:,j)
               s_buf(:, j) = s_tot(:, j)
            end do


            !!!!!!!!!!!!!!!!!!!
            ! Gain calculations
            !!!!!!!!!!!!!!!!!!!

            ! Precompute filtered signal for calibration
            if (do_oper(samp_G) .or. do_oper(samp_rcal) .or. do_oper(samp_acal)) then
               call self%downsample_tod(s_orb_tot(:,1), ext)
               allocate(  s_invN(ext(1):ext(2), ndet))      ! s * invN
               allocate(s_lowres(ext(1):ext(2), ndet))      ! s * invN
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  if (do_oper(samp_G) .or. do_oper(samp_rcal) .or. .not. self%orb_abscal) then
                     s_buf(:,j) = s_tot(:,j)
                     call fill_all_masked(s_buf(:,j), mask(:,j), ntod, &
                     &  trim(self%operation)=='sample', &
                     &  real(self%scans(i)%d(j)%sigma0, sp), &
                     &  handle, self%scans(i)%chunk_num)
                     call self%downsample_tod(s_buf(:,j), ext, &
                          & s_lowres(:,j))!, mask(:,j))
                  else
                     call self%downsample_tod(s_orb_tot(:,j), ext, &
                          & s_lowres(:,j))!, mask(:,j))
                  end if
               end do
               s_invN = s_lowres
               call multiply_inv_N(self, i, s_invN,   sampfreq=self%samprate_lowres, pow=0.5d0)
               call multiply_inv_N(self, i, s_lowres, sampfreq=self%samprate_lowres, pow=0.5d0)
            end if

            ! Prepare for absolute calibration
            if (do_oper(samp_acal) .or. do_oper(samp_rcal)) then
               call update_status(status, "Prepping for absolute calibration")
               call wall_time(t1)
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  if (do_oper(samp_acal)) then
                     if (self%orb_abscal) then
                        s_buf(:, j) = real(self%gain0(0),sp) * (s_tot(:, j) - s_orb_tot(:, j)) + &
                             & real(self%gain0(j) + self%scans(i)%d(j)%dgain,sp) * s_tot(:, j)
                     else
                        s_buf(:, j) = real(self%gain0(j) + self%scans(i)%d(j)%dgain,sp) * s_tot(:, j)
                     end if
                  else
                     s_buf(:,j) = real(self%gain0(0) + self%scans(i)%d(j)%dgain,sp) * s_tot(:, j)
                  end if
               end do
               call accumulate_abscal(self, i, mask, s_buf, s_lowres, s_invN, A_abscal, b_abscal, handle)

            end if

            ! Fit gain
            if (do_oper(samp_G)) then
               call update_status(status, "Fitting gain")
               call wall_time(t1)
               call calculate_gain_mean_std_per_scan(self, i, s_invN, mask, s_lowres, s_tot, handle)
               call wall_time(t2); t_tot(4) = t_tot(4) + t2-t1
            end if

            ! Fit correlated noise
            if (do_oper(samp_N)) then
               call update_status(status, "Fitting correlated noise")
               call wall_time(t1)
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  if (do_oper(samp_mono)) then
                     s_buf(:,j) = s_tot(:,j)-s_mono(:,j)
                  else
                     s_buf(:,j) = s_tot(:,j)
                  end if
               end do
               call sample_n_corr(self, handle, i, mask, s_buf, n_corr, pix(:,:,1), .false.)
!!               do j = 1, ndet
!!                  n_corr(:,j) = sum(n_corr(:,j))/ size(n_corr,1)
!!               end do
               call wall_time(t2); t_tot(3) = t_tot(3) + t2-t1
            else
               n_corr = 0.
            end if

            ! Compute noise spectrum
            if (do_oper(samp_N_par)) then
               call update_status(status, "Computing noise power")
               call wall_time(t1)
               call sample_noise_psd(self, handle, i, mask, s_tot, n_corr)
               call wall_time(t2); t_tot(6) = t_tot(6) + t2-t1
            end if

            !! Compute chisquare
            if (do_oper(calc_chisq)) then
               call wall_time(t1)
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  s_buf(:,j) =  s_sl(:,j) + s_orb_tot(:,j)
                  if (do_oper(samp_mono)) s_buf(:,j) =  s_buf(:,j) + s_mono(:,j)
                  call self%compute_chisq(i, j, mask(:,j), s_sky(:,j), &
                       & s_buf(:,j), n_corr(:,j), verbose=verbose)
               end do
               call wall_time(t2); t_tot(7) = t_tot(7) + t2-t1
            end if

            !*******************
            ! Compute binned map
            !*******************

            ! Get calibrated map
            if (do_oper(bin_map)) then
               call update_status(status, "Computing binned map")
               allocate (d_calib(nout, ntod, ndet))
               d_calib = 0
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  inv_gain = 1.0/real(self%scans(i)%d(j)%gain, sp)
                  if (remove_solar_dipole) then
                     d_calib(1, :, j) = (self%scans(i)%d(j)%tod - n_corr(:, j))* &
                        & inv_gain - s_tot(:, j) + s_sky(:, j) - s_sol_tot(:, j)
                     if (nout > 1) d_calib(2, :, j) = d_calib(1, :, j) - (s_sky(:, j) - s_sol_tot(:, j)) - s_bp(:, j) ! Residual
                  else
                     d_calib(1, :, j) = (self%scans(i)%d(j)%tod - n_corr(:, j))* &
                        & inv_gain - s_tot(:, j) + s_sky(:, j) - s_bp(:, j)
                     if (nout > 1) d_calib(2, :, j) = d_calib(1, :, j) - s_sky(:, j) + s_bp(:, j) ! Residual
                  end if

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

               if (.true. .and. do_oper(bin_map) .and. self%scans(i)%chunk_num == 6878) then
                  call int2string(self%scanid(i), scantext)
                  do k = 1, self%ndet
                     open(78,file='chains_WMAP/tod_'//trim(self%label(k))//'_pid'//scantext//'.dat', recl=1024)
                     write(78,*) "# Sample   uncal_TOD (mK)  n_corr (mK) cal_TOD (mK)   res (mK)"// &
                          & "    sig (mK)    s_sol_dip  (mK)   s_orb_dip (mK)"
                     do j = 1, ntod
                        inv_gain = 1.0/real(self%scans(i)%d(k)%gain, sp)
                        !write(78,*) j, self%scans(i)%d(k)%tod*inv_gain, n_corr(j, k)*inv_gain, d_calib(1,j,k), d_calib(2,j,k),  s_sky(j,k), s_sol_tot(j,k), s_orb_tot(j,k)
                        write(78,*) j, self%scans(i)%d(k)%tod(j), n_corr(j, k), d_calib(1,j,k), s_sky(j,k), s_orb_tot(j,k), s_sol_tot(j,k), mask(j, k), inv_gain
                     end do
                     close(78)
                  end do
               end if

               ! Bin the calibrated map
               call bin_differential_TOD(self, d_calib, pix,  &
                 & psi, flag, self%x_im, procmask, b_map, M_diag, i, &
                 & do_oper(prep_relbp))
               deallocate(d_calib)
               call wall_time(t8); t_tot(19) = t_tot(19) + t8
            end if

            do j = 1, ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               masked_var = masked_variance(s_sky(:, j), mask(:, j))
               if (masked_var == 9999999999999) then
                 self%scans(i)%d(j)%accept = .false.
                 dipole_mod(self%scanid(i), j) = 0
                 !cycle
               else
                 dipole_mod(self%scanid(i), j) = masked_var
               end if
            end do

            ! Clean up
            deallocate (n_corr, s_sky, s_orbA, s_orbB, s_orb_tot, s_tot, s_buf)
            deallocate ( s_sl, s_slA, s_slB, s_sky_prop)
            deallocate (mask, mask2, pix, psi, flag, s_solA, s_solB, s_sol_tot)
            if (allocated(s_lowres)) deallocate (s_lowres)
            if (allocated(s_invN)) deallocate (s_invN)
            deallocate(s_bp, s_bp_prop)
         end do

         call mpi_allreduce(mpi_in_place, dipole_mod, size(dipole_mod), MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)

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

         call update_status(status, "Finished main loop iteration") 
      end do main_it

      call update_status(status, "Running allreduce on M_diag")
      call mpi_allreduce(mpi_in_place, M_diag, size(M_diag), &
           & MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)
      call update_status(status, "Running allreduce on b")
      call mpi_allreduce(mpi_in_place, b_map, size(b_map), &
           & MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)


      np0 = self%info%np
      allocate (cg_tot(0:np0 - 1, nmaps+1))

      ! write out M_diag, b_map to fits.
      cg_tot = b_map(self%info%pix, 1:nmaps, 1)
      call write_fits_file_iqu(trim(prefix)//'b'//trim(postfix), cg_tot, outmaps)
      cg_tot = M_diag(self%info%pix, 1:nmaps)
      call write_fits_file_iqu(trim(prefix)//'M'//trim(postfix), cg_tot, outmaps)

      where (M_diag == 0d0)
         M_diag = 1d0
      end where

      ! Conjugate Gradient solution to (P^T Ninv P) m = P^T Ninv d, or Ax = b
      call update_status(status, "Allocating cg arrays")
      allocate (r     (0:npix-1, nmaps+1))
      allocate (r0    (0:npix-1, nmaps+1))
      allocate (q     (0:npix-1, nmaps+1))
      allocate (v     (0:npix-1, nmaps+1))
      allocate (p     (0:npix-1, nmaps+1))
      allocate (s     (0:npix-1, nmaps+1))
      allocate (phat  (0:npix-1, nmaps+1))
      allocate (shat  (0:npix-1, nmaps+1))
      allocate (cg_sol(0:npix-1, nmaps+1, nout))

      cg_sol = 0.0d0
      epsil = 1.0d-2

      if (self%myid_shared ==0) write(*,*) '  Running BiCG'

      do l=1, nout
         if (self%myid_shared==0) write(*,*) '    Solving for ', trim(adjustl(self%labels(l)))
         call update_status(status, "Starting bicg-stab")
         r  = b_map(:, :, l)
         r0 = b_map(:, :, l)
         if (maxval(r) == 0) cycle
         delta_r = sum(r**2/M_diag)
         ! WMAP's metric was |Ax-b|/|b| < 10^-8, which essentially is 
         delta_0 = delta_r
         delta_s = delta_s

         omega = 1
         alpha = 1

         rho_new = sum(r0*r)
         bicg: do i = 1, i_max
            rho_old = rho_new
            rho_new = sum(r0*r)
            if (i==1) then
                p = r
            else
                beta = (rho_new/rho_old)/(alpha/omega)
                p = r + beta*(p - omega*v)
            end if
            phat = p/M_diag
            v = 0d0
            call update_status(status, "Calling p=Av")
            call compute_Ax(self, phat, v, self%x_im, procmask, i)
            call mpi_allreduce(MPI_IN_PLACE, v, size(v), &
                              & MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)
            alpha = rho_new/sum(r0*v)
            cg_sol(:,:,l) = cg_sol(:,:,l) + alpha*phat
            if (write_cg_iter) then
               cg_tot = cg_sol(self%info%pix, 1:nmaps+1, l)
               call write_fits_file_iqu(trim(prefix)//'cg'//trim(str(l))//'_iter'//trim(str(2*(i-1)))//trim(postfix), cg_tot, outmaps)
            end if
            s = r - alpha*v

            shat = s/M_diag

            delta_s = sum(s*shat)
            if (self%myid_shared==0 .and. verbose) then 
                write(*,101) 2*i-1, delta_s/delta_0
                101 format (6X, I4, ':   delta_s/delta_0:',  2X, ES9.2)
            end if
            if (delta_s .le. (delta_0*epsil)) exit bicg
            q = 0d0
            call update_status(status, "Calling  q= A shat")
            call compute_Ax(self, shat, q, self%x_im, procmask, i)
            call mpi_allreduce(MPI_IN_PLACE, q, size(q), &
                              & MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)
            omega = sum(q*s)/sum(q**2)
            cg_sol(:,:,l) = cg_sol(:,:,l) + omega*shat
            if (mod(i, 10) == 1) then
               call update_status(status, 'r = b - Ax')
               r = 0d0
               call compute_Ax(self, cg_sol(:,:,l), r, self%x_im, procmask, i)
               call mpi_allreduce(MPI_IN_PLACE, r, size(r), &
                                 & MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)
               r = b_map(:, :, l) - r
            else
               call update_status(status, 'r = s - omega*t')
               r = s - omega*q
            end if

            if (write_cg_iter) then
               cg_tot = cg_sol(self%info%pix, 1:nmaps+1, l)
               call write_fits_file_iqu(trim(prefix)//'cg'//trim(str(l))//'_iter'//trim(str(2*(i-1)+1))//trim(postfix), cg_tot, outmaps)
            end if

            delta_r = sum(r**2/M_diag)
            if (self%myid_shared==0 .and. verbose) then 
                write(*,102) 2*i, delta_r/delta_0
                102 format (6X, I4, ':   delta_r/delta_0:',  2X, ES9.2)
            end if
            if ((delta_r .le. (delta_0*epsil))) exit bicg
         end do bicg
      end do

      if (remove_solar_dipole) then
        ! Need to add solar dipole that was removed in timestream back in
        write(*,*) 'You have not yet added the signal you removed back in.'
        ! 
        ! In the tod_orb module, we were looping over detector i and tod index
        ! j. Here, I want to loop over pixel indices.
        ! b_dot = dot_product(vnorm, self%tod%pix2vec(:,pix(j,i)))
        ! s_orb(j,i) = T_CMB_DIP * b_dot
        !phi=4.607145626489432
        !theta=0.7278022980816355
        !vnorm = (/ sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) /)
        !do i=0, self%info%np
        !  cg_sol(self%info%pix(i), 1, 1) = cg_sol(self%info%pix(i), 1, 1) +  dot_product(vnorm, self%pix2vec(:,self%info%pix(i)))
        !end do
      end if



      do k = 1, self%output_n_maps
         do j = 1, nmaps
            outmaps(k)%p%map(:, j) = cg_sol(self%info%pix, j, k)
         end do
      end do

      cg_tot = cg_sol(self%info%pix, 1:nmaps+1, 1)
      !call write_fits_file(trim(prefix)//'S'//trim(postfix), cg_tot(:,nmaps+1), outmaps)

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

      ! Clean up temporary arrays
      deallocate(A_abscal, b_abscal, chisq_S)
      deallocate(b_map, M_diag, cg_tot)
      if (allocated(b_mono)) deallocate (b_mono)
      if (allocated(sys_mono)) deallocate (sys_mono)
      if (allocated(slist)) deallocate (slist)
      if (allocated(dipole_mod)) deallocate (dipole_mod)

      if (allocated(outmaps)) then
         do i = 1, nout
            call outmaps(i)%p%dealloc
         end do
         deallocate (outmaps)
      end if

      deallocate (map_sky)
      deallocate (cg_sol, r, s, q, r0, shat, p, phat, v)

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



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to save time-ordered-data chunk
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_tod_chunk(filename, tod)
     implicit none
     character(len=*),                   intent(in) :: filename
     real(sp),         dimension(:),     intent(in) :: tod
     ! Expects one-dimensional TOD chunk

     integer(i4b) :: unit, n_tod, t

     n_tod = size(tod)

     unit = getlun()
     open(unit,file=trim(filename), recl=1024)
     write(unit,*) '# TOD value in mK'
     do t = 1, n_tod
        write(unit,fmt='(e16.8)') tod(t)
     end do
     close(unit)
   end subroutine write_tod_chunk


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to save map array to fits file 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_fits_file(filename, array, outmaps)
     implicit none
     character(len=*),                   intent(in) :: filename
     real(dp),         dimension(0:),    intent(in) :: array
     class(map_ptr),   dimension(:),     intent(in) :: outmaps

     integer(i4b) :: np0, m

     do m = 0, size(array) - 1
        outmaps(1)%p%map(m, 1) = array(m)
     end do

     call outmaps(1)%p%writeFITS(filename)

   end subroutine write_fits_file

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to save map array to fits file 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_fits_file_iqu(filename, array, outmaps)
     implicit none
     character(len=*),                    intent(in) :: filename
     real(dp),         dimension(0:, 1:), intent(in) :: array
     class(map_ptr),   dimension(:),      intent(in) :: outmaps

     outmaps(1)%p%map = array

     call outmaps(1)%p%writeFITS(filename)

   end subroutine write_fits_file_iqu

end module comm_tod_WMAP_mod
