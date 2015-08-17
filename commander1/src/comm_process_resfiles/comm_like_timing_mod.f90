module comm_like_timing_mod
   use comm_like_utils
   use comm_proc_utils
   use comm_lowl_mod
   implicit none

  ! **************************************************************************
  ! *  comm_like_optimization_mod -- module for searching for best basis set *
  ! *                                                                        *
  ! *          H. K. Eriksen, E. Gjerløw (University of Oslo)                *
  ! *                                                                        * 
  ! *                               and                                      *
  ! *                                                                        *
  ! *                       I. K. Wehus  (JPL)                               *
  ! *                                                                        *
  ! *                                                                        *
  ! **************************************************************************

contains
   subroutine run_timing
      implicit none

      character(len=512) :: mapfile, maskfile, covfile, beamfile, clfile
      character(len=512) :: temp, prefix, binfile, fname
      integer(i4b)      :: lmax, lhigh_min, lhigh_max, lmax_comp, seed, n_h
      integer(i4b)      :: ordering, nside, nmaps, nval, n_t, numtemp, n_p
      integer(i4b)      :: curr_nmodes_T, curr_nmodes_P, npix
      integer(i4b)      :: nmodes_min, nmodes_max, num_eval
      integer(i4b)      :: npoints, currnmodes, numcomp, ind, ind2, ind1
      integer(i4b), allocatable, dimension(:)     ::    map2mask, bincomp_map
      integer(i4b), allocatable, dimension(:)     ::    lmax2lhigh, points_arr
      integer(i4b), allocatable, dimension(:)     ::    buffer_1d_int
      real(dp),     allocatable, dimension(:,:)   ::    cls, map, mask, T
      real(dp),     allocatable, dimension(:,:)   ::    weights, cov, Y, Yp
      real(dp),     allocatable, dimension(:,:)   ::    alms, beam 
      real(dp),     allocatable, dimension(:,:)   ::    PY, C, P_highl
      real(dp),     allocatable, dimension(:,:)   ::    buffer, S_mat, P
      real(dp),     allocatable, dimension(:,:)   ::    buffer2, C_transf
      real(dp),     allocatable, dimension(:,:)   ::    invF, V, map_1d
      real(dp),     allocatable, dimension(:,:,:)   ::  window
      real(dp),     allocatable, dimension(:)       ::  mask_1d, W
      real(dp),     allocatable, dimension(:)       ::  res_array
      real(dp),     allocatable, dimension(:)       ::  buffer_1d
      real(dp),     pointer,     dimension(:,:)     :: pixwin
      complex(dpc), allocatable, dimension(:,:,:)   :: alms_cmplx
      real(dp)  :: BBamp, scale, dummy
      type(planck_rng)  :: handle

      integer(i4b)      :: i, j, k, l, m

!      if (iargc() /= 19) then
!         write(*,*) '    Calculate the smallest basis that still gives consistent Fisher matrix error bars'
!         write(*,*) '    Options:  [mapfile] [maskfile] [covfile] [beamfile] [clfile_basis] [clfile_fisher] [binfile] [lmax_data] [lmax_basis_min] [lmax_basis_max] [l_maxcomp]'
!         write(*,*) '                 [basis] [BB amp] [T reg] [P reg] [seed] [outprefix] [relative tolerance]'
!         write(*,*) '        where [basis]     = {pseudoalm, eigen_N, eigen_StoN, eigen_S+N}'
!         stop
!      end if

      if (iargc() /= 11) then
         write(*,*) '    Do an iteratively more fine-grained timing test. Uses the pseudoalm basis'
         write(*,*) '     Options: [mapfile] [maskfile] [covfile] [beamfile] [clfile] [lmax] '
         write(*,*) '                  [min_nmodes] [max_nmodes] [num_eval] [outprefix]'
         stop
      end if

      !Read parameters
      call getarg(2,mapfile)
      call getarg(3,maskfile)
      call getarg(4,covfile)
      call getarg(5,beamfile)
      call getarg(6,clfile)
      call getarg(7,temp)
      read(temp,*) lmax
      call getarg(8,temp)
      read(temp,*) nmodes_min
      call getarg(9,temp)
      read(temp,*) nmodes_max
      call getarg(10,temp)
      read(temp,*) num_eval
      call getarg(11,prefix)

      call rand_init(handle, seed)

      write(*,*) 'Reading data:'
      write(*,*) '     Map        = ', trim(mapfile)
      write(*,*) '     Mask       = ', trim(maskfile)
      write(*,*) '     Covariance = ', trim(covfile)
      write(*,*) '     Cl file    = ', trim(clfile)

      call read_fiducial_spectrum(clfile, cls)

      ! Read map and mask
      i    = getsize_fits(trim(mapfile), ordering=ordering, nside=nside, nmaps=nmaps)
      numcomp = (lmax+1)**2
      npix = 12*nside**2
      n_p  = npix    * nmaps
      n_h  = numcomp * nmaps

      allocate(map(0:npix-1,nmaps), mask(0:npix-1,nmaps), mask_1d(n_p), map_1d(n_p, 1))
      call read_map(mapfile,  map)
      call read_map(maskfile, mask)
      where (mask < 0.5d0)
         mask = 0.d0
      elsewhere
         mask = 1.d0
      end where
      do i = 1, nmaps
         mask_1d((i-1)*npix+1:i*npix) = mask(:,i)
         !Add if regularization noise is needed
!         do j = 0, npix-1
!            map_1d((i-1)*npix+j+1, 1)  = map(j,i) + reg(i) * rand_gauss(handle)
!         end do
      end do

      ! Set up map2mask conversion array
      nval = count(mask > 0.5d0)
      n_t  = count(mask(:,1) > 0.5d0)
      allocate(map2mask(nval))
      k = 1
      do j = 1, nmaps
         do i = 0, npix-1
            if (mask(i,j) > 0.5d0) then
               map2mask(k) = (j-1)*npix+i+1
               k           = k+1
            end if
         end do
      end do
      print *, n_t

      ! Read beam and pixel window
      allocate(beam(0:lmax,nmaps), window(0:lmax,nmaps,nmaps), weights(2*nside,nmaps))
      call read_pixwin(nside, nmaps, pixwin)
      call read_ringweights(nside, weights)
      call read_beam(beamfile, beam)
      beam = beam * pixwin(0:lmax,:)
      do l = 0, lmax
         do i = 1, nmaps
            do j = 1, nmaps
               window(l,i,j) = beam(l,i)*beam(l,j)
            end do
         end do
      end do
   
      ! Read covariance matrix
      allocate(cov(n_p,n_p))
      call comm_read_covmat(covfile, nmaps, cov, scale)
      map_1d = map_1d * scale ! Assume map is in same units as covariance matrix
   
      ! Add regularization noise to covariance matrix
!      ind = 1
!      do j = 1, nmaps
!         do i = 1, npix
!            cov(ind,ind) = cov(ind,ind) + reg(j)**2
!            ind          = ind+1
!         end do
!      end do

      ! Compute spherical harmonics
      write(*,*) 'Computing spherical harmonics'
      allocate(Y(n_p, n_h), Yp(n_h,n_p), alms(numcomp,nmaps), alms_cmplx(nmaps,0:lmax,0:lmax))
      ind1 = 1
      do j = 1, nmaps
         ind2 = 1
         do l = 0, lmax
            do m = -l, l
               alms        = 0.d0
               alms(ind2,j) = 1.d0
               call convert_real_to_complex_alms_dp(alms, alms_cmplx)
               if (nmaps == 1) then
                  call alm2map(nside, lmax, lmax, alms_cmplx, map(:,1))
               else
                  call alm2map(nside, lmax, lmax, alms_cmplx, map)
               end if
               do i = 1, nmaps
                  Y((i-1)*npix+1:i*npix,ind1) = map(:,i)
               end do
               ind1 = ind1+1
               ind2 = ind2+1
            end do
         end do
      end do
   
      !Compute inverse spherical harmonic
      ind1 = 1
      do j = 1, nmaps
         do i = 0, npix-1
            map      = 0.d0
            map(i,j) = 1.d0
            if (nmaps == 1) then
               call map2alm(nside, lmax, lmax, map(:,1), alms_cmplx, [0.d0,0.d0], weights)
            else
               call map2alm(nside, lmax, lmax, map,      alms_cmplx, [0.d0,0.d0], weights)
            end if
            call convert_complex_to_real_alms(alms_cmplx, alms)
            do k = 1, nmaps
               Yp((k-1)*numcomp+1:k*numcomp,ind1) = alms(:,k)
            end do
            ind1 = ind1+1
         end do
      end do
      deallocate(alms, alms_cmplx)

      allocate(C(nval, nval))
      C = cov(map2mask, map2mask)

      allocate(V(nval, nval))
      allocate(W(nval))
      allocate(P(nval,nval))

      !Do pixel based reference first
      V = 0
      do i = 1, nval
         V(i, i) = 1.d0
      end do

      dummy = perform_test(nval, num_eval, W, V, n_t, beam, window, &
         & map_1d(map2mask, :), cls, 2, lmax, lmax, n_h, nmaps, &
         & Yp(:, map2mask), C)
      write(*, *) 'Pixel based:', dummy

      !Then form pseudoalm
      call dgemm('N','N',nval,nval,n_h,1.d0,Y(map2mask,:),nval,&
            & Yp(:,map2mask),n_h,0.d0,P,nval)
      P(1:n_t,n_t+1:nval) = 0.d0
      P(n_t+1:nval,1:n_t) = 0.d0
      write(*,*) 'Computing eigen vectors'
      allocate(buffer(n_t, n_t), buffer2(n_t, n_t), buffer_1d(n_t))
      buffer = P(1:n_t, 1:n_t)
      V = 0 
      call get_eigen_decomposition(buffer, buffer_1d, buffer2)
      W(1:n_t) = buffer_1d
      V(1:n_t, 1:n_t) = buffer2
      deallocate(buffer, buffer2, buffer_1d)
      allocate(buffer(nval-n_t, nval-n_t), buffer2(nval-n_t, nval-n_t), &
         & buffer_1d(nval-n_t))
      buffer = P(n_t+1:nval, n_t+1:nval)
      call get_eigen_decomposition(buffer, buffer_1d, buffer2)
      W(n_t+1:) = buffer_1d
      V(n_t+1:, n_t+1:) = buffer2
      deallocate(buffer, buffer2, buffer_1d, P)


      !Do the actual test
      write(*, *) 'Running test'
      allocate(res_array(2))
      allocate(points_arr(2))

      points_arr(1) = nmodes_min
      points_arr(2) = nmodes_max
      res_array(1) = perform_test(nmodes_min, num_eval, W, V, n_t, beam, window, &
         & map_1d(map2mask, :), cls, 2, lmax, lmax, n_h, nmaps, &
         & Yp(:,map2mask), C)
      res_array(2) = perform_test(nmodes_max, num_eval, W, V, n_t, beam, window, &
         & map_1d(map2mask, :), cls, 2, lmax, lmax, n_h, nmaps, &
         & Yp(:,map2mask), C)
      call output_results(prefix, points_arr, res_array)
      npoints = 2
      do while (npoints < nmodes_max-nmodes_min + 1)
         allocate(buffer_1d(2*npoints - 1))
         allocate(buffer_1d_int(2*npoints-1))
         j = 1
         do i = 1, npoints - 1
            buffer_1d(j) = res_array(i)
            buffer_1d_int(j) = points_arr(i)
            j = j + 1
            if (points_arr(i+1)+1 == points_arr(i)) cycle
            currnmodes = (points_arr(i+1) + points_arr(i)) / 2
            buffer_1d(j) = perform_test(currnmodes, num_eval, W, V, n_t, beam, window, &
               & map_1d(map2mask, :), cls, 2, lmax, lmax, n_h, nmaps, &
               & Yp(:,map2mask), C)
            buffer_1d_int(j) = currnmodes
            j = j + 1
         end do
         buffer_1d(j) = res_array(npoints)
         buffer_1d_int(j) = points_arr(npoints)
         npoints = j
         deallocate(res_array)
         deallocate(points_arr)
         allocate(res_array(npoints))
         allocate(points_arr(npoints))
         res_array = buffer_1d(:npoints)
         points_arr = buffer_1d_int(:npoints)
         call output_results(prefix, points_arr, res_array)
         deallocate(buffer_1d, buffer_1d_int)
      end do
   end subroutine run_timing

   function perform_test(nmodes, num_eval, eigvals, eigvecs, n_t, beam, &
        &  w, d_in, cls, lmin_like, lmax_like, lmax, n_h, nmaps, Yp, N_cov_in) &
        & result(comptime)
      implicit none

      integer(i4b), intent(in)  :: nmodes, num_eval, n_t, lmin_like
      integer(i4b), intent(in)  :: lmax_like, lmax, nmaps, n_h
      real(dp), dimension(:), intent(in)        :: eigvals
      real(dp), dimension(:, :), intent(in)     :: eigvecs, N_cov_in, Yp, d_in
      real(dp), dimension(0:, :), intent(in)    :: beam, cls
      real(dp), dimension(0:, :, :), intent(in)    :: w

      real(dp), dimension(:, :), allocatable    :: N_cov, P_harm, d
      real(dp)  :: time1, time2
      real(dp)  :: comptime, dummy
      integer(i4b)      :: handle, nmodes_t, nmodes_p, i

      write(*, *) 'Timing ', nmodes, ' modes'

      if (mod(nmodes, 3) == 0) then
         nmodes_t = nmodes / 3
         nmodes_p = 2 * nmodes / 3
      else if (mod(nmodes - 1, 3) == 0) then
         nmodes_t = (nmodes - 1) / 3 +1
         nmodes_p = (nmodes - 1) / 3 * 2
      else if (mod(nmodes - 2, 3) == 0) then
         nmodes_t = (nmodes - 2) / 3
         nmodes_p = (nmodes - 2) / 3 * 2 + 2
      else
         stop "This should never happen: Nothing divisible by 3!"
      end if

      !HACK
      nmodes_t = 2411
      nmodes_p = nmodes - nmodes_t

      allocate(N_cov(nmodes, nmodes))
      allocate(P_harm(n_h, nmodes))
      allocate(d(nmodes,1))
      call get_projection_operators(eigvecs, nmodes_t, nmodes_p, &
         & n_t, N_cov_in, N_cov, Yp, P_harm, d_in, d)

      call comm_lowl_initialize_object_direct(lmin_like, lmax_like, lmax, 1, &
         & n_h, nmodes, nmaps, d, N_cov, w, beam, P_harm, cls, 1.d0, 1.d30, &
         & handle)
      deallocate(N_cov, P_harm, d)
      call cpu_time(time1)
      do i = 1, num_eval
         dummy = comm_lowl_compute_lnL(cls, handle=handle)
         print *, 'lnl', dummy
      end do
      call cpu_time(time2)
      comptime = (time2-time1) / num_eval
      print *, 'time', comptime
      call comm_lowl_deallocate_object(handle)

   end function perform_test

   subroutine output_results(prefix, x, y)
      implicit none
      character(len=*), intent(in)      :: prefix
      integer(i4b), dimension(:),intent(in)     :: x
      real(dp), dimension(:), intent(in)        :: y

      integer(i4b)      :: unit, i

      unit = comm_getlun()

      open(unit, file=trim(prefix) // '_timing.dat')
      do i = 1, size(x)
         write(unit, fmt='(I10, F15.8)') x(i), y(i)
      end do
      close(unit)
   end subroutine output_results

   subroutine get_projection_operators(eigvecs, tmodes, pmodes, n_t, noise, & 
         & noise_transf, ylmp, pylm, d_in, d_out)
      implicit none

      real(dp), dimension(:, :), intent(in)     :: eigvecs, noise, ylmp, d_in
      real(dp), dimension(:, :), intent(out)    :: noise_transf, pylm, d_out
      integer(i4b), intent(in)      :: tmodes, pmodes, n_t

      real(dp), allocatable, dimension(:, :) :: projmat, buffer
      integer(i4b)      :: nval, nmodes, n_h

      nval = size(ylmp, 2)
      n_h = size(ylmp, 1)
      nmodes = tmodes+pmodes

      allocate(projmat(nmodes, nval))

      projmat = 0

      projmat(1:tmodes, :) = transpose(eigvecs(:, n_t-tmodes+1:n_t))
      projmat(tmodes+1:, :) = transpose(eigvecs(:, nval-pmodes+1:nval))

      call dgemm('N','N',nmodes,1,nval,1.d0,projmat,nmodes,d_in,nval,0.d0, &
         & d_out,nmodes)

      allocate(buffer(nval, nmodes))
      !Form PCP^T
      call dgemm('N', 'T', nval, nmodes, nval, 1.d0, noise, nval, projmat, & 
         & nmodes, 0.d0, buffer, nval)
      call dgemm('N', 'N', nmodes, nmodes, nval, 1.d0, projmat, nmodes, &
         & buffer, nval, 0.d0, noise_transf, nmodes)
      deallocate(buffer)

      !Form Y^TP^T
      call dgemm('N', 'T', n_h, nmodes, nval, 1.d0, ylmp, n_h, projmat, &
         & nmodes, 0.d0, pylm, n_h)

      deallocate(projmat)

   end subroutine get_projection_operators

end module comm_like_timing_mod
