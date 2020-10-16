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

   !*************************************************
   !    Convert integer to string
   !*************************************************
   character(len=20) function str(k)
       integer, intent(in) :: k
       write (str, *) k
       str = adjustl(str)
   end function str

   !*************************************************
   !            Random normal generator
   !*************************************************
   function rand_normal(mean,stdev) result(c)
            double precision :: mean,stdev,c,temp(2),theta,r
            if (stdev <= 0.0d0) then
               write(*,*) "Standard Deviation must be positive."
            else
               call RANDOM_NUMBER(temp)
               r=(-2.0d0*log(temp(1)))**0.5
               theta = 2.0d0*PI*temp(2)
           c= mean+stdev*r*sin(theta)
         end if
       end function



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
      constructor%output_n_maps = 2
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
      ! For K-band
      ! constructor%x_im = [-0.00067, 0.00536]

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

      integer(i4b) :: i, j, k, l, m, n, t, start_chunk, end_chunk, chunk_size, ntod, ndet
      integer(i4b) :: nside, npix, nmaps, naccept, ntot, ext(2), nscan_tot, nhorn
      integer(i4b) :: ierr, main_iter, n_main_iter, ndelta, ncol, n_A, np0, nout = 1
      real(dp)     :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, chisq_threshold
      real(dp)     :: t_tot(22)
      real(sp)     :: inv_gain
      real(sp), allocatable, dimension(:, :)     :: n_corr, s_sl, s_sky, s_orbA, s_orbB, mask, mask2, s_bp
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
      real(dp), allocatable, dimension(:, :, :)   :: b_tot, M_diag_tot
      real(dp), allocatable, dimension(:, :)   :: cg_tot, r_tot
      character(len=512) :: prefix, postfix, prefix4D, filename
      character(len=2048) :: Sfilename
      character(len=4)   :: ctext, myid_text
      character(len=6)   :: samptext, scantext
      character(len=512), allocatable, dimension(:) :: slist
      type(shared_1d_int) :: sprocmask, sprocmask2
      real(sp), allocatable, dimension(:, :, :, :) :: map_sky
      type(shared_2d_dp) :: sA_map, scg_sol
      type(shared_3d_dp) :: sb_map, sb_mono, sM_diag
      class(comm_map), pointer :: condmap
      class(map_ptr), allocatable, dimension(:) :: outmaps

      ! conjugate gradient parameters
      integer(i4b) :: i_max
      real(dp) :: delta_0, delta_old, delta_new, epsil
      real(dp) :: alpha, beta, g, f_quad
      real(dp), allocatable, dimension(:, :, :) :: cg_sol, r, s, d, q
      logical(lgt) :: write_cg_iter=.false.

      ! biconjugate gradient parameters
      real(dp) :: rho_old, rho_new
      real(dp) :: omega, delta_r, delta_s
      real(dp), allocatable, dimension(:, :, :) :: r0, shat, p, phat, v

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
            map_in(i, j)%p%map = map_in(i, j)%p%map ! uK to K
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
      naccept = 0; ntot = 0
      ! Using nmaps + 1 to include the spurious component
      allocate (M_diag(nout, nmaps, 0:npix-1))
      allocate (b_map(nout, nmaps, 0:npix-1))
      M_diag(:,:,:) = 0d0
      do i = 1, self%nscan

         ! Short-cuts to local variables
         call wall_time(t1)
         ntod = self%scans(i)%ntod
         ndet = self%ndet

         ! Set up local data structure for current scan
         allocate (n_corr(ntod, ndet))                 ! Correlated noise in V
         allocate (s_sky(ntod, ndet))                  ! Sky signal in uKcmb
         allocate (s_sky_prop(ntod, ndet, 2:ndelta))   ! Sky signal in uKcmb
         allocate (s_orbA(ntod, ndet))                 ! Orbital dipole (beam A) in uKcmb
         allocate (s_orbB(ntod, ndet))                 ! Orbital dipole (beam B) in uKcmb
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
         call project_sky_differential(self, map_sky(:, :, :, 1), pix, psi, flag, &
                  & self%x_im, sprocmask%a, i, s_sky, mask)


         ! Construct orbital dipole template
         call wall_time(t1)
         call self%orb_dp%p%compute_orbital_dipole_4pi(i, pix(:,:,1), psi(:,:,1), s_orbA)
         call self%orb_dp%p%compute_orbital_dipole_4pi(i, pix(:,:,2), psi(:,:,2), s_orbB)
         call wall_time(t2); t_tot(2) = t_tot(2) + t2-t1

         !estimate the correlated noise
         ! Add orbital dipole to total signal
         s_buf = 0.d0
         do j = 1, ndet
            s_tot(:, j) = s_sky(:, j)! + (1+self%x_im((j+1)/2))*s_orbA(:,j) - &
                                     ! & (1-self%x_im((j+1)/2))*s_orbB(:,j)
            s_buf(:, j) = s_tot(:, j)
         end do
         n_corr(:, :) = 0d0
         call wall_time(t2); t_tot(7) = t_tot(7) + t2 - t1

         !*******************
         ! Compute binned map
         !*******************

         ! Get calibrated map
         allocate (d_calib(nout, ntod, ndet))
         do j = 1, ndet
            inv_gain = 1.0/real(self%scans(i)%d(j)%gain, sp)
            d_calib(1, :, j) = (self%scans(i)%d(j)%tod - n_corr(:, j))* &
               & inv_gain - s_tot(:, j) + s_sky(:, j)

            if (nout > 1) d_calib(2, :, j) = d_calib(1, :, j) - s_sky(:, j) ! Residual
            if (nout > 2) d_calib(3, :, j) = (n_corr(:, j) - sum(n_corr(:, j)/ntod))*inv_gain
            if (i == 1) then
               if (self%myid_shared == 0) then
                  filename = trim(prefix)//'calib_'//trim(str(j))// '.dat'
                  call write_tod_chunk(filename, d_calib(1,:,j))
                  filename = trim(prefix)//'res'//trim(str(j))//'.dat'
                  call write_tod_chunk(filename, d_calib(2,:,j))
                  filename = trim(prefix)//'sky_sig'//trim(str(j))//'.dat'
                  call write_tod_chunk(filename, s_sky(:,j))
               end if
            end if

         end do

         ! Bin the calibrated map
         call bin_differential_TOD(self, d_calib, pix,  &
                & psi, flag, self%x_im, sprocmask%a, b_map, M_diag, i)
         deallocate (n_corr, s_sky, s_orbA, s_orbB, s_tot, s_buf, s_sky_prop, d_calib)
         deallocate (mask, mask2, pix, psi, flag)
         if (allocated(s_lowres)) deallocate (s_lowres)
         if (allocated(s_invN)) deallocate (s_invN)
         call wall_time(t8); t_tot(19) = t_tot(19) + t8

      end do

      call update_status(status, "Running allreduce on M_diag")
      call mpi_allreduce(mpi_in_place, M_diag, size(M_diag), &
           & MPI_DOUBLE_PRECISION, MPI_SUM, self%comm_shared, ierr)
      call update_status(status, "Running allreduce on b")
      call mpi_allreduce(mpi_in_place, b_map, size(b_map), &
           & MPI_DOUBLE_PRECISION, MPI_SUM, self%comm_shared, ierr)
      ! Mean-subtracting b_map to remove fitting of spurious monopole
      !do j=1, nmaps
      !   b_map(l,j,:) = b_map(l,j,:) - sum(b_map(l,j,:))/size(b_map(l,j,:))
      !end do

      where (M_diag .eq. 0d0)
         M_diag = 1d0
      end where

      ! Conjugate Gradient solution to (P^T Ninv P) m = P^T Ninv d, or Ax = b
      call update_status(status, "Allocating cg arrays")
      np0 = self%info%np
      allocate (r(nout, nmaps, 0:npix-1))
      allocate (s(nout, nmaps, 0:npix-1))
      allocate (q(nout, nmaps, 0:npix-1))
      allocate (d(nout, nmaps, 0:npix-1))
      allocate (cg_sol(nout, nmaps, 0:npix-1))
      allocate (cg_tot(nmaps, 0:np0 - 1))

      cg_sol = 0.0d0
      epsil = 1.0d-3

      i_max = int(npix**0.5)

      allocate (r0(nout, nmaps, 0:npix-1))
      allocate (v(nout, nmaps, 0:npix-1))
      allocate (p(nout, nmaps, 0:npix-1))
      allocate (phat(nout, nmaps, 0:npix-1))
      allocate (shat(nout, nmaps, 0:npix-1))
      do l=1, nout
         call update_status(status, "Starting bicg-stab")
         r(l, :, :)  = b_map(l, :, :)
         r0(l, :, :) = b_map(l, :, :)
         delta_0 = sum(r(l,:,:)**2/M_diag(l,:,:))
         delta_r = delta_0
         delta_s = delta_0

         omega = 1
         alpha = 1

         rho_new = sum(r0(l,:,:)*r(l,:,:))
         bicg: do i = 1, i_max
            rho_old = rho_new
            rho_new = sum(r0(l,:,:)*r(l,:,:))
            if (i==1) then
                p(l,:,:) = r(l,:,:)
            else
                beta = (rho_new/rho_old)/(alpha/omega)
                p(l,:,:) = r(l,:,:) + beta*(p(l,:,:) - omega*v(l,:,:))
            end if
            phat(l,:,:) = p(l,:,:)/M_diag(l,:,:)
            v(l,:,:) = 0d0
            call update_status(status, "Calling p=Av")
            call compute_Ax(self, phat, v, self%x_im, sprocmask%a, i, l)
            call mpi_allreduce(MPI_IN_PLACE, v(l,:,:), size(v(l,:,:)), &
                              & MPI_DOUBLE_PRECISION, MPI_SUM, self%comm_shared, ierr)
            alpha = rho_new/sum(r0(l,:,:)*v(l,:,:))
            s(l,:,:) = r(l,:,:) - alpha*v(l,:,:)

            shat(l,:,:) = s(l,:,:)/M_diag(l,:,:)

            delta_s = sum(s(l,:,:)*shat(l,:,:))
            if (delta_s .le. (delta_0*epsil**2)) then
                cg_sol(l,:,:) = cg_sol(l,:,:) + alpha*phat(l,:,:)
                exit bicg
            end if
            if (self%myid_shared==0) then 
                write(*,101) i, delta_s/delta_0
                101 format (I3, ':   delta_s/delta_0:',  2X, ES9.2)
            end if
            q(l,:,:) = 0d0
            call update_status(status, "Calling  t= A shat")
            call compute_Ax(self, shat, q, self%x_im, sprocmask%a, i, l)
            call mpi_allreduce(MPI_IN_PLACE, q(l,:,:), size(q(l,:,:)), &
                              & MPI_DOUBLE_PRECISION, MPI_SUM, self%comm_shared, ierr)
            omega = sum(q(l,:,:)*s(l,:,:))/sum(q(l,:,:)**2)
            cg_sol(l,:,:) = cg_sol(l,:,:) + alpha*phat(l,:,:) + omega*shat(l,:,:)
            if (mod(i, 50) == 1) then
               call update_status(status, 'r = b - Ax')
               r(l,:,:) = 0d0
               call compute_Ax(self, cg_sol, r, self%x_im, sprocmask%a, i, l)
               call mpi_allreduce(MPI_IN_PLACE, r(l,:,:), size(r(l,:,:)), &
                                 & MPI_DOUBLE_PRECISION, MPI_SUM, self%comm_shared, ierr)
               r(l, :, :) = b_map(l, :, :) - r(l, :, :)
            else
               call update_status(status, 'r = s - omega*t')
               r(l,:,:) = s(l,:,:) - omega*q(l,:,:)
            end if

            if (write_cg_iter) then
               cg_tot = cg_sol(l, 1:nmaps, self%info%pix + 1)
               do m = 0, np0 - 1
                  do n = 1, nmaps
                     outmaps(1)%p%map(m, n) = cg_tot(n, m)*1.d3 ! convert from mK to uK
                  end do
               end do
               call outmaps(1)%p%writeFITS(trim(prefix)//'cg_iter_'//trim(str(i))//trim(postfix))
            end if

            delta_r = sum(r(l,:,:)**2/M_diag(l,:,:))
            if (self%myid_shared==0) then 
                write(*,102) i, delta_r/delta_0
                102 format (I3, ':   delta_r/delta_0:',  2X, ES9.2)
            end if
            if ((delta_r .le. (delta_0*epsil**2))) exit bicg

         end do bicg
      end do

      allocate (r_tot(nmaps, 0:np0 - 1))
      cg_tot = cg_sol(1, 1:nmaps, self%info%pix + 1)
      r_tot = cg_sol(2, 1:nmaps, self%info%pix + 1)
      call update_status(status, "Got total map arrays")
      ! I think there is a factor of ncpus that needs to be divided out, but
      ! want to be clear about that.
      do i = 0, np0 - 1
         do j = 1, nmaps
            outmaps(1)%p%map(i, j) = cg_tot(j, i)*1.d6 ! convert from K to uK
            outmaps(2)%p%map(i, j) = r_tot(j, i)*1.d6 ! convert from K to uK
         end do
      end do

      map_out%map = outmaps(1)%p%map
      call outmaps(1)%p%writeFITS(trim(prefix)//'cg_sol'//trim(postfix))
      call outmaps(2)%p%writeFITS(trim(prefix)//'resid'//trim(postfix))

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

end module comm_tod_WMAP_mod
