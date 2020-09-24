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

   integer(i4b), parameter :: N_test = 19
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
   integer(i4b), parameter :: sub_zodi = 17
   logical(lgt), dimension(N_test) :: do_oper

   type, extends(comm_tod) :: comm_WMAP_tod
      class(orbdipole_pointer), allocatable :: orb_dp ! orbital dipole calculator
      real(dp), allocatable, dimension(:)  :: x_im    ! feedhorn imbalance parameters

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
      type(comm_params), intent(in) :: cpar
      integer(i4b), intent(in) :: id_abs
      class(comm_mapinfo), target     :: info
      character(len=128), intent(in) :: tod_type
      class(comm_WMAP_tod), pointer    :: constructor

      integer(i4b) :: i, nside_beam, lmax_beam, nmaps_beam, ndelta
      character(len=512) :: datadir
      logical(lgt) :: pol_beam

      ! Set up WMAP specific parameters
      allocate (constructor)
      constructor%output_n_maps = 3
      constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
      constructor%nhorn = 2

      ! Initialize beams
      nside_beam = 512
      nmaps_beam = 3
      pol_beam = .true.
      constructor%nside_beam = nside_beam

      !initialize the common tod stuff
      call constructor%tod_constructor(cpar, id_abs, info, tod_type)
      allocate (constructor%x_im(constructor%ndet/2))
      constructor%x_im(:) = 0.0d0

      !TODO: this is LFI specific, write something here for wmap
      call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label)
      do i = 1, constructor%ndet
         if (mod(i, 2) == 1) then
            constructor%partner(i) = i + 1
         else
            constructor%partner(i) = i - 1
         end if
         constructor%horn_id(i) = (i + 1)/2
      end do

      ! Read the actual TOD
      ! TODO: this probabl needs a seperate fucntion/ modifications for wmap
      call constructor%read_tod(constructor%label)

      call constructor%precompute_lookups()

      datadir = trim(cpar%datadir)//'/'

      ! Initialize bandpass mean and proposal matrix
      call constructor%initialize_bp_covar(trim(datadir)//cpar%ds_tod_bp_init(id_abs))

      !load the instrument file
      call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

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

      integer(i4b) :: i, j, k, l, t, start_chunk, end_chunk, chunk_size, ntod, ndet
      integer(i4b) :: nside, npix, nmaps, naccept, ntot, ext(2), nscan_tot, nhorn
      integer(i4b) :: ierr, main_iter, n_main_iter, ndelta, ncol, n_A, np0, nout = 1
      real(dp)     :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, chisq_threshold
      real(dp)     :: t_tot(22)
      real(sp)     :: inv_gain
      real(sp), allocatable, dimension(:, :)     :: n_corr, s_sl, s_sky, s_orb, mask, mask2, s_bp
      real(sp), allocatable, dimension(:, :)     :: s_mono, s_buf, s_tot, s_zodi
      real(sp), allocatable, dimension(:, :)     :: s_invN, s_lowres
      real(sp), allocatable, dimension(:, :, :)   :: s_sky_prop, s_bp_prop
      real(sp), allocatable, dimension(:, :, :)   :: d_calib
      real(dp), allocatable, dimension(:)       :: A_abscal, b_abscal
      real(dp), allocatable, dimension(:, :)     :: chisq_S, m_buf
      real(dp), allocatable, dimension(:, :)     :: A_map, dipole_mod
      real(dp), allocatable, dimension(:, :, :)   :: b_map, b_mono, sys_mono, M_diag
      integer(i4b), allocatable, dimension(:, :, :)   :: pix, psi
      integer(i4b), allocatable, dimension(:, :)     :: flag
      real(dp), allocatable, dimension(:, :, :)   :: b_tot
      character(len=512) :: prefix, postfix, prefix4D, filename
      character(len=2048) :: Sfilename
      character(len=4)   :: ctext, myid_text
      character(len=6)   :: samptext, scantext
      character(len=512), allocatable, dimension(:) :: slist
      type(shared_1d_int) :: sprocmask, sprocmask2
      real(sp), allocatable, dimension(:, :, :, :) :: map_sky
      type(shared_2d_dp) :: sA_map
      type(shared_3d_dp) :: sb_map, sb_mono
      class(comm_map), pointer :: condmap
      class(map_ptr), allocatable, dimension(:) :: outmaps

      ! conjugate gradient parameters
      integer(i4b) :: i_max
      real(dp) :: delta_0, delta_old, delta_new, epsil
      real(dp) :: alpha, beta, g
      real(dp), allocatable, dimension(:, :, :) :: cg_sol, r, s, d, q
      real(dp), allocatable, dimension(:, :) :: A
      integer(i4b) :: lpoint, rpoint, lpsi, rpsi, sgn
      real(dp) ::  inv_sigmasq, x_im
      real(dp) :: dA, dB, d1

      if (iter > 1) self%first_call = .false.
      call int2string(iter, ctext)
      call update_status(status, "tod_start"//ctext)

      t_tot = 0.d0
      call wall_time(t5)

      ! Set up full-sky map structures
      call wall_time(t1)
      chisq_threshold = 100.d0
      ndet = self%ndet
      nhorn = self%nhorn
      ndelta = size(delta, 3)
      nside = map_out%info%nside
      nmaps = map_out%info%nmaps
      npix = 12*nside**2
      nout = self%output_n_maps
      nscan_tot = self%nscan_tot
      chunk_size = npix/self%numprocs_shared
      if (chunk_size*self%numprocs_shared /= npix) chunk_size = chunk_size + 1
      allocate (map_sky(nmaps, self%nobs, 0:ndet, ndelta))
      allocate (chisq_S(ndet, ndelta))
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
            map_in(i, j)%p%map = 1d-6*map_in(i, j)%p%map ! uK to K
            call map_in(i, j)%p%bcast_fullsky_map(m_buf)
            do k = 1, self%nobs
               map_sky(:, k, i, j) = m_buf(self%ind2pix(k), :) ! uK to K
            end do
         end do
         do k = 1, self%nobs
            do l = 1, nmaps
               map_sky(l, k, 0, j) = sum(map_sky(l, k, 1:ndet, j))/ndet
            end do
         end do
      end do
      deallocate (m_buf)

      ! Set up shared processing mask
      call init_shared_1d_int(self%myid_shared, self%comm_shared, &
           & self%myid_inter, self%comm_inter, [npix], sprocmask)
      call sync_shared_1d_int_map(sprocmask, self%procmask%info%pix, &
           & nint(self%procmask%map(:, 1)))
      call init_shared_1d_int(self%myid_shared, self%comm_shared, &
           & self%myid_inter, self%comm_inter, [npix], sprocmask2)
      call sync_shared_1d_int_map(sprocmask2, self%procmask2%info%pix, &
           & nint(self%procmask2%map(:, 1)))
      call wall_time(t2); t_tot(9) = t2 - t1

      call update_status(status, "tod_init")

      ! Perform main analysis loop
      ! In LFI, this was in the loop
      ! if (do_oper(sel_data)) then
      !    naccept = 0; ntot    = 0
      ! end if
      naccept = 0; ntot = 0
      allocate (M_diag(nout, nmaps, npix))
      allocate (b_map(nout, nmaps, npix))
      M_diag(:,:,:) = 1d-3
      do i = 1, self%nscan

         !write(*,*) "Processing scan: ", i, self%scans(i)%d%accept
         call update_status(status, "tod_loop1")

         ! Short-cuts to local variables
         call wall_time(t1)
         ntod = self%scans(i)%ntod
         ndet = self%ndet

         ! Set up local data structure for current scan
         allocate (n_corr(ntod, ndet))                 ! Correlated noise in V
         allocate (s_sky(ntod, ndet))                  ! Sky signal in uKcmb
         allocate (s_sky_prop(ntod, ndet, 2:ndelta))    ! Sky signal in uKcmb
         allocate (s_orb(ntod, ndet))                  ! Orbital dipole in uKcmb
         allocate (s_buf(ntod, ndet))                  ! Buffer
         allocate (s_tot(ntod, ndet))                  ! Sum of all sky compnents
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
         !call self%symmetrize_flags(flag)
         call wall_time(t2); t_tot(11) = t_tot(11) + t2 - t1
         call update_status(status, "tod_decomp")

         ! Construct sky signal template
         call project_sky_differential(self, map_sky(:, :, :, 1), pix, psi, flag, &
                  & self%x_im, sprocmask%a, i, s_sky, mask)
         call update_status(status, "Finished projecting sky")
         !if (self%first_call) then
         !   do j = 1, ndet
         !      if (all(mask(:, j) == 0)) self%scans(i)%d(j)%accept = .false.
         !      if (self%scans(i)%d(j)%sigma0 <= 0.d0) self%scans(i)%d(j)%accept = .false.
         !   end do
         !end if

         !estimate the correlated noise
         s_buf = 0.d0
         do j = 1, ndet
            s_tot(:, j) = s_sky(:, j)
            s_buf(:, j) = s_tot(:, j)
         end do
         n_corr(:, :) = 0d0
         !call sample_n_corr(self, handle, i, mask, s_buf, n_corr)

         ! Compute chisquare
         !TODO: this but probably different code?
         !if (do_oper(calc_chisq)) then
         !   call wall_time(t1)
         !do j = 1, ndet
         !   s_buf(:,j) =  s_sl(:,j) + s_orb(:,j)
         !         if (do_oper(samp_mono)) s_buf(:,j) =  s_buf(:,j) + s_mono(:,j)
         !   call self%compute_chisq(i, j, mask(:, j), s_sky(:, j), &
         !        & s_buf(:, j), n_corr(:, j))
         !end do
         call wall_time(t2); t_tot(7) = t_tot(7) + t2 - t1
         !end if

         ! Select data
         call wall_time(t1)

         ! Compute binned map
         ! By the end of this loop, you won't have access to the raw TOD loop
         ! have a pixA array, pixB, psiA, psiB, flags_A, flags_B
         ! Also, compute P^T Ninv d (over raw tod)
         ! Potentially decompress these arrays during the CG solving?
         allocate (d_calib(nout, ntod, ndet))
         do j = 1, ndet
            inv_gain = 1.0/real(self%scans(i)%d(j)%gain, sp)
            d_calib(1, :, j) = (self%scans(i)%d(j)%tod - n_corr(:, j))* &
                & inv_gain - s_tot(:, j) + s_sky(:, j)
            if (nout > 1) d_calib(2, :, j) = d_calib(1, :, j) - s_sky(:, j) ! Residual
            if (nout > 2) d_calib(3, :, j) = (n_corr(:, j) - sum(n_corr(:, j)/ntod))*inv_gain
         end do

         call bin_differential_TOD(self, d_calib, pix,  &
                & psi, flag, self%x_im, b_map, M_diag, i)
         call update_status(status, "Finished binning TOD")
         deallocate (n_corr, s_sky, s_orb, s_tot, s_buf, s_sky_prop, d_calib)
         deallocate (mask, mask2, pix, psi, flag)
         if (allocated(s_lowres)) deallocate (s_lowres)
         if (allocated(s_invN)) deallocate (s_invN)
         call wall_time(t8); t_tot(19) = t_tot(19) + t8

         !update scanlist with new timing info
         self%scans(i)%proctime = t8 - t1
         self%scans(i)%n_proctime = self%scans(i)%n_proctime + 1
         write (slist(i), *) self%scanid(i), '"', trim(self%hdfname(i)), '"', &
             & real(self%scans(i)%proctime/self%scans(i)%n_proctime, sp), &
             & real(self%spinaxis(i, :), sp)

      end do

      ! sum up all of the individual threads so that everybody has access to the
      ! hits map and the binned map P^T Ninv d.
      call mpi_allreduce(MPI_IN_PLACE, b_map, size(b_map), &
                        & MPI_DOUBLE_PRECISION, MPI_SUM, self%comm_shared, ierr)
      call mpi_allreduce(MPI_IN_PLACE, M_diag, size(M_diag), &
                        & MPI_DOUBLE_PRECISION, MPI_SUM, self%comm_shared, ierr)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Summary of the conjugate gradient algorithm as I am implementing it.
      ! The mapmaking problem is in sum, P^T Ninv P m = P^T Ninv d
      ! where m is a map, d is the time-ordered data, Ninv is the diagonal noise
      ! covariance matrix, and P is the pointing matrix. The right-hand side is
      ! b = P^T Ninv d
      ! and is constructed iteratively
      ! for t in time_index:
      !     b[pixA[t]] += d[t]/sigma**2
      !     b[pixB[t]] -= d[t]/sigma**2
      ! (for temperature only data)
      !
      ! For a given map m, you can construct b = P^T Ninv P m as follows:
      ! for t in time_index:
      !     b[pixA[t]] += (m[pixA[t]] - m[pixB[t]])/sigma**2
      !     b[pixB[t]] -= (m[pixA[t]] - m[pixB[t]])/sigma**2
      !
      ! In linear algebra terms, we can write this whole matrix as
      ! Am = b
      ! where b is the binned matrix P^T Ninv d. This is computed as we load the
      ! data, and is a map.
      ! 
      ! The CG algorithm iteratively solves for the m in Am = b.
      ! If you start with a guess m_0, there is a residual r_0 = b - Am_0
      ! r_0 also gives us a guess for which direction will reduce the error in
      ! the solution; we can call this direction d_0.
      !
      ! For every iteration, we start with a residual map, r_i, a direction map
      ! d_i, and an estimate of the solution x_i.
      ! 
      ! To get the i+1st estimate, compute the following;
      !
      ! q_i     = A.dot(d_i)
      ! You need to distribute q_i, because it includes contributions from
      ! pixels you don't already have
      !
      ! delta_i = r_i.dot(r_i)  can be done per core
      ! gamma_i = d_i.dot(q_i)
      ! mpi_allreduce(gamma_i) as well.
      ! alpha_i = delta_i/gamma_i
      !
      ! x_{i+1} = x_i + alpha_i*d_i
      ! ! 1 out of 50 times, compute r_{i+1} = b - A.dot(x_{i+1})
      ! r_{i+1} = r_i - alpha_i*q_i
      !
      ! beta_{i+1} = r_{i+1}.dot(r_{i+1})/r_i.dot(r_i)
      ! beta_{i+1} = delta_{i+1}/delta_i
      ! allreduce(delta_{i+1})
      !
      ! d_{i+1} = r_{i+1} + beta_{i+1}*d_i
      !
      ! 
      ! Note that this can be simplified by computing q_i = A.dot(d_i) once per
      ! loop, so that r_{i+1} = r_i - alpha_i q_i
      ! 
      ! As long as you distribute q_i, the rest should be done per core.
      ! A is the only operator that requires distributing the result.
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! start a new loop that is the CG solver loop
      allocate (r(nout, nmaps, npix))
      allocate (s(nout, nmaps, npix))
      allocate (q(nout, nmaps, npix))
      allocate (d(nout, nmaps, npix))
      allocate (cg_sol(nout, nmaps, npix))
      ! allocate as a shared memory array, so that everything can access it.
      cg_sol(:, :, :) = 0.0d0
      epsil = 1.0d-16
      ! Properties to test to ensure that the CG solver is basically doing the
      ! correct thing
      ! Each direction is conjugate (A-orthogonal) to each other
      ! d_{i}.dot(A.dot(d_j)) = delta_{ij}
      ! Each residual is orthogonal,
      ! r_{i}.dot(r_{j}) = delta_{ij}

      i_max = 100

      !do l = 1, nout
      l = 1
         ! All threads have access to the same b_map, hence the same residual
         ! map r and the direction map d.
         r(l, :, :) = b_map(l, :, :)
         d(l,:,:) = r(l,:,:)/M_diag(l,:,:)
         ! delta_new is the same for all threads.
         delta_new = sum(r(l, :, :)*d(l, :, :))
         if (self%myid_shared==0) then 
            write (*, *), delta_new, 'delta_0'
         end if
         delta_0 = delta_new
         i = 0
         do while ((i .lt. i_max) .and. (delta_new .ge. (epsil**2)*delta_0))
            call update_status(status, 'While loop iteration')
            ! This is evaluating the matrix product q = (P^T Ninv P) d
            q(l,:,:) = 0d0
            do j = 1, self%nscan
               ntod = self%scans(j)%ntod
               ndet = self%ndet
               allocate (pix(ntod, ndet, nhorn))             ! Decompressed pointing
               allocate (psi(ntod, ndet, nhorn))             ! Decompressed pol angle
               allocate (flag(ntod, ndet))                   ! Decompressed flags
               do k = 1, self%ndet
                  ! decompress the data so we have one chunk of TOD in memory
                  ! In the working code above, i is looping over nscan, j is ndet...
                  call self%decompress_pointing_and_flags(j, k, pix(:, k, :), &
                      & psi(:, k, :), flag(:, k))
                  do t = 1, ntod
                     inv_sigmasq = 1/self%scans(j)%d(k)%sigma0**2
                     lpoint = self%pix2ind(pix(t, k, 1))
                     rpoint = self%pix2ind(pix(t, k, 2))
                     lpsi = psi(t, k, 1)
                     rpsi = psi(t, k, 2)
                     x_im = self%x_im((k + 1)/2)
                     sgn = (-1)**((k + 1)/2 + 1)
                     ! This is the model for each timestream
                     ! The sgn parameter is +1 for timestreams 13 and 14, -1
                     ! for timestreams 23 and 24, and also is used to switch
                     ! the sign of the polarization sensitive parts of the
                     ! model
                     dA = d(l, 1, lpoint) + sgn*(d(l, 2, lpoint)*self%cos2psi(lpsi) + d(l, 3, lpoint)*self%sin2psi(lpsi))
                     dB = d(l, 1, rpoint) + sgn*(d(l, 2, rpoint)*self%cos2psi(rpsi) + d(l, 3, rpoint)*self%sin2psi(rpsi))
                     d1 = (1 + x_im)*dA - (1 - x_im)*dB
                     ! Temperature
                     q(l, 1, lpoint) = q(l, 1, lpoint) + (1 + x_im)*d1*inv_sigmasq
                     q(l, 1, rpoint) = q(l, 1, rpoint) - (1 - x_im)*d1*inv_sigmasq
                     ! Q
                     q(l, 2, lpoint) = q(l, 2, lpoint) + (1 + x_im)*d1*self%cos2psi(lpsi)*sgn*inv_sigmasq
                     q(l, 2, rpoint) = q(l, 2, rpoint) - (1 - x_im)*d1*self%cos2psi(rpsi)*sgn*inv_sigmasq
                     ! U
                     q(l, 2, lpoint) = q(l, 2, lpoint) + (1 + x_im)*d1*self%sin2psi(lpsi)*sgn*inv_sigmasq
                     q(l, 2, rpoint) = q(l, 2, rpoint) - (1 - x_im)*d1*self%sin2psi(rpsi)*sgn*inv_sigmasq
                  end do
               end do
               deallocate (pix, psi, flag)
            end do
            ! q is populated using the decompressed pointing and flags from
            ! self, which is distributed among the cores. So it is necessary to
            ! add up the different cores' results.
            write(*, *), minval(q), maxval(q), 'minmax(q)'
            call mpi_allreduce(MPI_IN_PLACE, q(l,:,:), size(q(l,:,:)), &
                              & MPI_DOUBLE_PRECISION, MPI_SUM, self%comm_shared, ierr)
            write(*, *) minval(q), maxval(q), 'minmax(q) (after reduce)'
            ! In the pre-loop operations, d was available to every core. It
            ! makes sense that this should be the case within the loop.
            ! Similarly, every core has the same access to q. Therefore, g
            ! should be the same in every core.
            g = sum(d(l,:,:)*q(l,:,:))
            write(*, *) g, 'g'
            !call mpi_allreduce(MPI_IN_PLACE, g, 1, &
            !                  & MPI_DOUBLE_PRECISION, MPI_SUM, self%comm_shared, ierr)
            alpha = delta_new/g
            ! Every core has access to the same d, so for this to make sense, we
            ! need to have the same cg_sol available to every core.
            cg_sol(l,:,:) = cg_sol(l,:,:) + alpha*d(l,:,:)
            write (*, *) minval(cg_sol(l,:,:)), maxval(cg_sol(l,:,:)), 'minmax(cg_sol)'
            ! evaluating r = b - Ax
            if (mod(i, 50) == 0) then
               call update_status(status, 'iter % 50 == 0, recomputing residual')
               r(l,:,:) = 0d0
               do j = 1, self%nscan
                  ntod = self%scans(j)%ntod
                  ndet = self%ndet
                  allocate (pix(ntod, ndet, nhorn))             ! Decompressed pointing
                  allocate (psi(ntod, ndet, nhorn))             ! Decompressed pol angle
                  allocate (flag(ntod, ndet))                   ! Decompressed flags
                  do k = 1, self%ndet
                     ! evaluating the matrix operation r = b_map - Ax
                     call self%decompress_pointing_and_flags(j, k, pix(:, k, :), &
                         & psi(:, k, :), flag(:, k))
                     do t = 1, ntod
                        inv_sigmasq = 1/self%scans(j)%d(k)%sigma0**2
                        lpoint = self%pix2ind(pix(t, k, 1))
                        rpoint = self%pix2ind(pix(t, k, 2))
                        lpsi = psi(t, k, 1)
                        rpsi = psi(t, k, 2)
                        x_im = self%x_im((k + 1)/2)
                        sgn = (-1)**((k + 1)/2 + 1)
                        ! This is the model for each timestream
                        ! The sgn parameter is +1 for timestreams 13 and 14, -1
                        ! for timestreams 23 and 24, and also is used to switch
                        ! the sign of the polarization sensitive parts of the
                        ! model
                        dA = cg_sol(l, 1, lpoint) + sgn*(cg_sol(l, 2, lpoint)*self%cos2psi(lpsi) &
                                                &      + cg_sol(l, 3, lpoint)*self%sin2psi(lpsi))
                        dB = cg_sol(l, 1, rpoint) + sgn*(cg_sol(l, 2, rpoint)*self%cos2psi(rpsi) &
                                                &      + cg_sol(l, 3, rpoint)*self%sin2psi(rpsi))
                        d1 = (1 + x_im)*dA - (1 - x_im)*dB
                        ! Temperature
                        r(l, 1, lpoint) = r(l, 1, lpoint) - (1 + x_im)*d1*inv_sigmasq
                        r(l, 1, rpoint) = r(l, 1, rpoint) + (1 - x_im)*d1*inv_sigmasq
                        ! Q
                        r(l, 2, lpoint) = r(l, 2, lpoint) - (1 + x_im)*d1*self%cos2psi(lpsi)*sgn*inv_sigmasq
                        r(l, 2, rpoint) = r(l, 2, rpoint) + (1 - x_im)*d1*self%cos2psi(rpsi)*sgn*inv_sigmasq
                        ! U
                        r(l, 2, lpoint) = r(l, 2, lpoint) - (1 + x_im)*d1*self%sin2psi(lpsi)*sgn*inv_sigmasq
                        r(l, 2, rpoint) = r(l, 2, rpoint) + (1 - x_im)*d1*self%sin2psi(rpsi)*sgn*inv_sigmasq
                     end do
                  end do
                  deallocate (pix, psi, flag)
               end do
               call mpi_allreduce(MPI_IN_PLACE, r(l,:,:), size(r(l,:,:)), &
                                 & MPI_DOUBLE_PRECISION, MPI_SUM, self%comm_shared, ierr)
               r(l, :, :) = b_map(l, :, :) + r(l, :, :)
            else
               r(l, :, :) = r(l, :, :) - alpha*q(l, :, :)
            end if
            s(l,:,:) = r(l,:,:)/M_diag(l,:,:)
            !s(l, :, :) = r(l, :, :)
            delta_old = delta_new
            delta_new = sum(r(l, :, :)*s(l, :, :))
            if (self%myid_shared==0) then 
                write (*, *) delta_new, 'delta'
                write (*, *) g, 'd.T.dot(q)'
                write (*, *), minval(q(l,:,:)), maxval(q(l,:,:)), 'minmax(q)'
                write (*, *), minval(r(l,:,:)), maxval(r(l,:,:)), 'minmax(r)'
                write (*, *), minval(s(l,:,:)), maxval(s(l,:,:)), 'minmax(s)'
                write (*, *), minval(cg_sol(l,:,:)), maxval(cg_sol(l,:,:)), 'minmax(cg)'
            end if
            beta = delta_new/delta_old
            d(l, :, :) = s(l, :, :) + beta*d(l, :, :)
            i = i + 1

         end do
      !end do

      ! cg_sol should be the total array at this point
      write (*, *), minval(cg_sol), maxval(cg_sol), 'minmax(cg_sol)'

      ! cg_sol = wmap_estimate
      ! distribute x to map object
      ! in comm_map_mod
      ! each core will have a subset of the pixels and alms at all times.
      ! figure out for each core, which cores that pixel needs, then populate the
      ! array
      ! np is the number of pixels that each core has
      ! Basically, it's ordered by rings.
      ! the pixel value should be included in the array pix, ordered from the first
      ! pixel that that core has and the last one that that one has.

      ! each core has a pix array as part of the comm_mapinfo
      ! it's an array of pixel indices.
      ! In the map array, it has the values of those particular values.
      ! Take the common array of the map solution, on each core go through the list
      ! of cores, look up the value, and put it in the map array.

      ! get a do loop; for each core, it'll go through its own pix array, and put
      ! it in the map array for that core.

      ! a sort of equivalent thing is in line 333 of comm_tod_mapmaking_mod, except
      ! the pixels have already been distributed here, so I need to do that for
      ! cg_sol. It's looping over np0, which is the number of pixels that each core
      ! has.

      ! In the LFI one, pixel i is the pixel that comm knows about, but I need it
      ! to be referring to that pixel. Create a look-up table that will turn a
      ! pixel number in the core's array to the right indexing convention.

      ! I also think that I need to use the "shared_2d_dp, shared_3d_dp" types.

      ! I will be "done" when I have put cg_sol into outmaps(1).

      ! there is the "shared bmap" sb_map that I think needs to be distributed.

      ! initialize shared map object, sb_map
      if (sb_map%init) call dealloc_shared_3d_dp(sb_map)
      call init_shared_3d_dp(self%myid_shared, self%comm_shared, &
           & self%myid_inter, self%comm_inter, [nout, nmaps, npix], sb_map)

      ! distributing the cg_sol to sb_map
      if (sb_map%init) then
         do i = 0, self%numprocs_shared - 1
            ! Define chunk indices
            start_chunk = mod(self%myid_shared + i, self%numprocs_shared)*chunk_size
            end_chunk = min(start_chunk + chunk_size - 1, npix - 1)
            write (*, *) start_chunk, end_chunk, 'start_chunk, end_chunk'
            do while (start_chunk < npix)
               if (self%pix2ind(start_chunk) /= -1) exit
               start_chunk = start_chunk + 1
            end do
            do while (end_chunk >= start_chunk)
               if (self%pix2ind(end_chunk) /= -1) exit
               end_chunk = end_chunk - 1
            end do
            if (start_chunk < npix) start_chunk = self%pix2ind(start_chunk)
            if (end_chunk >= start_chunk) end_chunk = self%pix2ind(end_chunk)

            !call mpi_win_fence(0, sA_map%win, ierr)
            call mpi_win_fence(0, sb_map%win, ierr)

            ! Assign cg_sol to each chunk of the data
            do j = start_chunk, end_chunk
               sb_map%a(:, :, self%ind2pix(j) + 1) = sb_map%a(:, :, self%ind2pix(j) + 1) + &
                    & cg_sol(:, :, j)
            end do
         end do
      end if
      np0 = self%info%np
      allocate (b_tot(nout, nmaps, 0:np0 - 1))
      b_tot = sb_map%a(:, 1:nmaps, self%info%pix + 1)
      do i = 0, np0 - 1
         do j = 1, nmaps
            do k = 1, self%output_n_maps
               outmaps(k)%p%map(i, j) = b_tot(k, j, i)*1.d3 ! convert from mK to uK
            end do
         end do
      end do

      map_out%map = outmaps(1)%p%map

      ! Output maps to disk
      ! call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))

      ! This is the thing that I can check to see if it makes sense.
      call outmaps(1)%p%writeFITS(trim(prefix)//'map'//trim(postfix))
      !if (self%output_n_maps > 1) call outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
      !if (self%output_n_maps > 2) call outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
      !if (self%output_n_maps > 3) call outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix))
      !if (self%output_n_maps > 4) call outmaps(5)%p%writeFITS(trim(prefix)//'mono'//trim(postfix))
      !if (self%output_n_maps > 5) call outmaps(6)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
      !if (self%output_n_maps > 6) call outmaps(7)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
      !if (self%output_n_maps > 7) call outmaps(8)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))

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
      !deallocate(A_abscal, b_abscal, chisq_S)
      if (allocated(b_map)) deallocate (b_map)
      if (sA_map%init) call dealloc_shared_2d_dp(sA_map)
      if (sb_map%init) call dealloc_shared_3d_dp(sb_map)
      if (sb_mono%init) call dealloc_shared_3d_dp(sb_mono)
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

      if (sprocmask%init) call dealloc_shared_1d_int(sprocmask)
      if (sprocmask2%init) call dealloc_shared_1d_int(sprocmask2)
      deallocate (map_sky)

      call int2string(iter, ctext)
      call update_status(status, "tod_end"//ctext)

      ! Parameter to check if this is first time routine has been
      self%first_call = .false.

   end subroutine process_WMAP_tod

end module comm_tod_WMAP_mod
