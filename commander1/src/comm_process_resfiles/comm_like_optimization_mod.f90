module comm_like_optimization_mod
   use comm_like_utils
   use comm_proc_utils
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
   subroutine run_optimization_search
      implicit none

      character(len=512) :: mapfile, maskfile, covfile, beamfile, clfile_basis
      character(len=512) :: temp, prefix, clfile_fisher, binfile, fname
      integer(i4b)      :: lmax, lhigh_min, lhigh_max, lmax_comp, seed, n_h
      integer(i4b)      :: ordering, nside, nmaps, nval, n_t, numtemp, n_p
      integer(i4b)      :: ind, ind1, ind2, n, numcomp, unit, lhigh, npix, lmin
      integer(i4b)      :: npar_powspec, numcomp_cutoff_T, iter_count, l_maxcomp
      integer(i4b)      :: lmin_temp, lmax_temp, n_bincomp, numcomp_cutoff_P
      integer(i4b)      :: nhigh, curr_nmodes_T, curr_nmodes_P, curr_tot_nmodes
      integer(i4b)      :: high_nmodes_T, high_nmodes_P, low_nmodes_T
      integer(i4b)      :: low_nmodes_P
      integer(i4b), allocatable, dimension(:)     ::    map2mask, bincomp_map
      integer(i4b), allocatable, dimension(:)     ::    lmax2lhigh
      real(dp),     allocatable, dimension(:,:)   ::    cls_basis, map, mask, T
      real(dp),     allocatable, dimension(:,:)   ::    weights, cov, Y, Yp
      real(dp),     allocatable, dimension(:,:)   ::    alms, beam, cls_fisher
      real(dp),     allocatable, dimension(:,:)   ::    PY, C, P_highl
      real(dp),     allocatable, dimension(:,:)   ::    buffer, S_mat, P
      real(dp),     allocatable, dimension(:,:)   ::    buffer2, C_transf
      real(dp),     allocatable, dimension(:,:)   ::    invF, V
      real(dp),     allocatable, dimension(:,:,:)   ::  window
      real(dp),     allocatable, dimension(:)       ::  map_1d, mask_1d, W
      real(dp),     allocatable, dimension(:)       ::  buffer_1d, sigma_ref
      real(dp),     pointer,     dimension(:,:)     :: pixwin
      integer(i4b),              dimension(1000,3) :: bins
      complex(dpc), allocatable, dimension(:,:,:)   :: alms_cmplx
      logical(lgt), allocatable, dimension(:)   :: compres
      real(dp)  :: BBamp, scale, tolerance
      real(dp), dimension(3)    :: reg
      character(len=24) :: basis
      character(len=256)        :: line
      character(len=1), dimension(6) :: status
      character(len=5)  :: lhighstring
      type(planck_rng)  :: handle
      logical(lgt)      :: first, failed_on_first, tflag
      logical(lgt)      :: found_T, found_P, i_exist, inv_error

      integer(i4b)      :: i, j, k, l, m, dummy

      if (iargc() /= 19) then
         write(*,*) '    Calculate the smallest basis that still gives consistent Fisher matrix error bars'
         write(*,*) '    Options:  [mapfile] [maskfile] [covfile] [beamfile] [clfile_basis] [clfile_fisher]'
         write(*,*) '                 [binfile] [lmax_data] [lmax_basis_min] [lmax_basis_max] [l_maxcomp]'
         write(*,*) '                 [basis] [BB amp] [T reg] [P reg] [seed] [outprefix] [relative tolerance]'
         write(*,*) '        where [basis]     = {pseudoalm, eigen_N, eigen_StoN, eigen_S+N}'
         stop
      end if

      !Read parameters
      call getarg(2,mapfile)
      call getarg(3,maskfile)
      call getarg(4,covfile)
      call getarg(5,beamfile)
      call getarg(6,clfile_basis)
      call getarg(7,clfile_fisher)
      call getarg(8,binfile)
      call getarg(9,temp)
      read(temp,*) lmax
      call getarg(10,temp)
      read(temp,*) lhigh_min
      call getarg(11,temp)
      read(temp,*) lhigh_max
      call getarg(12,temp)
      read(temp,*) l_maxcomp
      call getarg(13,basis)
      call getarg(14,temp)
      read(temp,*) BBamp
      call getarg(15,temp)
      read(temp,*) reg(1)
      call getarg(16,temp)
      read(temp,*) reg(2)
      reg(3) = reg(2)
      call getarg(17,temp)
      read(temp,*) seed
      call getarg(18,prefix)
      call getarg(19, temp)
      read(temp,*) tolerance

      call rand_init(handle, seed)

      if (trim(basis) /= 'pseudoalm' .and. trim(basis) /= 'eigen_S+N' .and. & 
            & trim(basis) /= 'eigen_N' .and. trim(basis) /= 'eigen_StoN') then
         write(*,*) 'Invalid basis ', trim(basis)
         stop
      end if
         

      write(*,*) 'Reading data:'
      write(*,*) '     Map        = ', trim(mapfile)
      write(*,*) '     Mask       = ', trim(maskfile)
      write(*,*) '     Covariance = ', trim(covfile)
      write(*,*) '     Cl file (basis calculations)    = ', trim(clfile_basis)
      write(*,*) '     Cl file (Fisher calculations)    = ', trim(clfile_fisher)

      ! Read fiducial spectrum; only actually used for S/N bases
      call read_fiducial_spectrum(clfile_basis, cls_basis)
      !Read fiducial spectrum for Fisher matrix calculations
      call read_fiducial_spectrum(clfile_basis, cls_fisher)

      ! Read map and mask
      i    = getsize_fits(trim(mapfile), ordering=ordering, nside=nside, nmaps=nmaps)
      numcomp = (lmax+1)**2
      npix = 12*nside**2
      n_p  = npix    * nmaps
      n_h  = numcomp * nmaps

      allocate(map(0:npix-1,nmaps), mask(0:npix-1,nmaps), mask_1d(n_p), map_1d(n_p))
      call read_map(mapfile,  map)
      call read_map(maskfile, mask)
      where (mask < 0.5d0)
         mask = 0.d0
      elsewhere
         mask = 1.d0
      end where
      do i = 1, nmaps
         mask_1d((i-1)*npix+1:i*npix) = mask(:,i)
         do j = 0, npix-1
            map_1d((i-1)*npix+j+1)  = map(j,i) + reg(i) * rand_gauss(handle)
         end do
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
      ind = 1
      do j = 1, nmaps
         do i = 1, npix
            cov(ind,ind) = cov(ind,ind) + reg(j)**2
            ind          = ind+1
         end do
      end do

      ! Read templates for marginalization
      numtemp = iargc()-19
      if (numtemp > 0) then
         allocate(T(n_p,numtemp))
         T = 0.d0
         do k = 1, iargc()-19
            call getarg(19+k,mapfile)
            call read_map(mapfile, map)
            do i = 1, nmaps
               T((i-1)*npix+1:i*npix,k) = map(:,i)
            end do
         end do
         T   = T * scale ! Assume templates are in same units as covariance
         cov = cov + matmul(T,transpose(T))
         deallocate(T)
      end if

      if (nmaps == 1) then
         write(*,fmt='(a,f8.3)') '      I sky fractions = ', sum(mask(:,1))/npix
      else
         write(*,fmt='(a,3f8.3)') '      I/Q/U sky fractions = ', sum(mask(:,1))/npix, &
               & sum(mask(:,2))/npix, sum(mask(:,3))/npix
      end if

      ! Read binfile
      npar_powspec = 0
      unit = comm_getlun()
      open(unit,file=trim(binfile))
      do while (.true.)
         read(unit,'(A)',end=511) line
         if (line(1:1) == '#' .or. trim(line)=='') cycle
         read(line,*) lmin_temp, lmax_temp, status
         do i = 1, 6
            if (status(i) /= 'S') cycle
            npar_powspec = npar_powspec+1
            bins(npar_powspec,1) = lmin_temp
            bins(npar_powspec,2) = lmax_temp
            bins(npar_powspec,3) = i
         end do
      end do
511   close(unit)
      
      ! Compute spherical harmonics
      write(*,*) 'Computing spherical harmonics'
      n = n_h
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

      allocate(sigma_ref(npar_powspec))
      allocate(invF(npar_powspec, npar_powspec))
      fname = 'pixel_lhigh_00000_iter_00000_invF.dat'
      inquire(file=trim(fname), exist=i_exist)
      if (i_exist) then
         write(*, *) 'Reading pixel based information from file'
         unit = comm_getlun()
         open(unit, file=trim(fname))
         do i = 1, npar_powspec
            read(unit, *) dummy, dummy, dummy, dummy, sigma_ref(i)
         end do
         close(unit)
      else
         write(*, *) 'Calculating pixel based inverse Fisher'
         !Do pixel based first to provide reference
         allocate(PY(nval, n_h))
         PY = Y(map2mask, :)
   
         invF = 0
         call get_invfisher(cls_fisher, window, PY, C, bins(1:npar_powspec, :), & 
            & invF, inv_error)
         call dump_matrix(invF, 'pixbased_invF.dat')
         if (inv_error) stop "Cannot invert pixel based Fisher matrix. Stopping."
         do i = 1, npar_powspec
            sigma_ref(i) = sqrt(invF(i, i))
         end do
         deallocate(PY)
         call output_invfisher_diags('pixel', 0, bins, invF, 0)
      end if

      !Figure out which bins should be included in the comparison
      n_bincomp = 0
      do i = 1, npar_powspec
         if (bins(i, 2) <= l_maxcomp) then
            n_bincomp = n_bincomp+1
         end if
      end do
      allocate(bincomp_map(n_bincomp))
      allocate(compres(n_bincomp))
      n = 1
      do i = 1, npar_powspec
         if (bins(i, 2) <= l_maxcomp) then
            bincomp_map(n) = i
            n = n + 1
         end if
      end do

      allocate(V(nval, nval))
      allocate(W(nval))

      !Do the actual search - loop over cutoff ls
      do lhigh = lhigh_min, lhigh_max
         write(*, *) 'Current lhigh: ', lhigh
         allocate(P_highl(nval,nval))

         !This is the maximum number of modes we can admit (polarization will
         !have 4 modes less due to no monopole and dipole)
         numcomp_cutoff_T = (lhigh+1)**2
         numcomp_cutoff_P = 2*numcomp_cutoff_T - 8

         ! Set up lmax to lhigh conversion array
         nhigh = nmaps * (lhigh+1)**2
         if (nmaps == 3) nhigh = nhigh - 8 ! No polarization monopoles or dipoles
         allocate(lmax2lhigh(nhigh))
         k = 1
         j = 1
         do i = 1, nmaps
            do l = 0, lmax
               do m = -l, l
                  if (l <= lhigh .and. (i == 1 .or. l>=2)) then
                     lmax2lhigh(k) = j
                     k             = k+1
                  end if
                  j = j+1
               end do
            end do
         end do

         call dgemm('N','N',nval,nval,nhigh,1.d0,Y(map2mask,lmax2lhigh),nval,&
            & Yp(lmax2lhigh,map2mask),nhigh,0.d0,P_highl,nval)

         write(*, *) 'Setting up basis'
         allocate(P(nval, nval))
         if (trim(basis) == 'pseudoalm') then
            P = P_highl
         else if (trim(basis) == 'eigen_N') then
            P = C
            call invert_matrix(P)
         else if (trim(basis) == 'eigen_StoN') then
            P = C
            call invert_matrix(P)
            allocate(buffer(nval,nval), S_mat(nval,nval))
            ! Compute S^1/2 * invN * S^1/2
            call get_basis_S(BBamp, 0.5d0, cls_basis, beam, Y(map2mask,:), S_mat)
            call dgemm('N','N',nval,nval,nval,1.d0,S_mat,nval,P,nval,0.d0,buffer,nval)
            call dgemm('N','N',nval,nval,nval,1.d0,buffer,nval,S_mat,nval,0.d0,P,nval)
            deallocate(buffer, S_mat)
         else if (trim(basis) == 'eigen_S+N') then
            allocate(S_mat(nval,nval))
            call get_basis_S(BBamp, 1.d0, cls_basis, beam, Y(map2mask,:), S_mat)
            P = C + S_mat
            deallocate(S_mat)
         end if
         ! Project out high-l modes by P_highl * invN * P_highl^t
         if (lhigh < lmax .and. (trim(basis) == 'eigen_N' .or. & 
               & trim(basis) == 'eigen_StoN' .or. trim(basis) == 'eigen_S+N')) &
               & then
            write(*,*) 'Projecting out high-l modes'
            allocate(buffer(nval,nval))
            call dgemm('N','N',nval,nval,nval,1.d0,P_highl,nval,P,nval,0.d0,buffer,nval)
            call dgemm('N','T',nval,nval,nval,1.d0,buffer,nval,P_highl,nval,0.d0,P,nval)
            deallocate(buffer)
         end if
         deallocate(P_highl)
         ! Nullify temperature-polarization cross-elements
         P(1:n_t,n_t+1:nval) = 0.d0
         P(n_t+1:nval,1:n_t) = 0.d0

         !Get the set of eigenvectors and eigenvalues corresponding to this
         !choice of lhigh
         write(*,*) 'Computing eigen vectors for ', trim(basis)
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

         call int2string(lhigh, lhighstring)
         fname = trim(prefix) // '_lhigh_' // lhighstring  // '_W.dat'
         unit = comm_getlun()
         open(unit, file=trim(fname))
         do i = 1, nval
            if (i <= n_t) then
               write(unit, *) i, abs(W(i) / maxval(W(1:n_t))), W(i)
            else
               write(unit, *) i, abs(W(i) / maxval(W(n_t+1:nval))), W(i)
            end if
         end do
         close(unit)

         high_nmodes_T = min(numcomp_cutoff_T, n_t)
         high_nmodes_P = min(numcomp_cutoff_P, nval - n_t)
         curr_nmodes_T = high_nmodes_T
         curr_nmodes_P = high_nmodes_P
         curr_tot_nmodes = curr_nmodes_T+curr_nmodes_P
         low_nmodes_T = 0
         low_nmodes_P = 0
         found_T = .False.
         found_P = .False.
         failed_on_first = .False.
         first = .True.
         call open_subresfile(unit, prefix, lhigh)
         tflag = .True.
         iter_count = 1
         write(*, *) 'Doing bisection search between ', low_nmodes_T, ' and ', high_nmodes_T, ' for temperature, and ', low_nmodes_P, ' and ', high_nmodes_P, ' for polarization'
         do while ((.not. found_T .or. .not. found_P)  .and. .not. failed_on_first)
            if (.not. first) then
               if (found_T) then
                  tflag = .False.
               else if (found_P) then
                  tflag = .True.
               end if
               if (tflag) then
                  curr_nmodes_T = int((high_nmodes_T + low_nmodes_T) / 2)
               else
                  curr_nmodes_P = int((high_nmodes_P + low_nmodes_P) / 2)
               end if
            end if
            curr_tot_nmodes = curr_nmodes_T + curr_nmodes_P
            print *, 'iter', iter_count
            print *, 'curr_nmodes_T:', curr_nmodes_T
            print *, 'curr_nmodes_P:', curr_nmodes_P
            allocate(C_transf(curr_tot_nmodes, curr_tot_nmodes))
            allocate(PY(curr_tot_nmodes, n_h))
            PY = 0
            call get_projection_operators(V, curr_nmodes_T, curr_nmodes_P, &
               & n_t, C, C_transf, Y(map2mask, :), PY)
            invF = 0
            call get_invfisher(cls_fisher, window, PY, C_transf, bins(1:npar_powspec, :), &
              &  invF, inv_error)
            if (.not. inv_error) then
               call output_invfisher_diags(prefix, lhigh, bins, invF, iter_count)
            end if
            deallocate(C_transf, PY)

            if (inv_error) then
               compres = .False.
            else
               do i = 1, n_bincomp
                  compres(i) = sqrt(invF(bincomp_map(i), bincomp_map(i))) < &
                     & (1.d0 + tolerance) * sigma_ref(bincomp_map(i))
               end do
            end if

            if (.not. all(compres)) then
               write(*, *) 'Exceeded tolerance limit'
               if (tflag) then
                  low_nmodes_T = curr_nmodes_T
                  curr_nmodes_T = high_nmodes_T
               else 
                  low_nmodes_P = curr_nmodes_P
                  curr_nmodes_P = high_nmodes_P
               end if
            else
               write(*, *) 'Success!'
               if (tflag) then
                  high_nmodes_T = curr_nmodes_T
               else
                  high_nmodes_P = curr_nmodes_P
               end if
            end if
            if (high_nmodes_T-low_nmodes_T == 1) then
               found_T = .True.
               curr_nmodes_T = high_nmodes_T
            end if
            if (high_nmodes_P - low_nmodes_P == 1) then
               found_P = .True.
               curr_nmodes_P = high_nmodes_P
            end if
            if (.not. first) tflag = .not. tflag

            if ((high_nmodes_T == low_nmodes_T) .or. &
                  & (high_nmodes_P == low_nmodes_P)) then
               if (first) then 
                  failed_on_first = .True.
                  write(*, *) 'Failed on first try.'
               else
                  stop "This should never happen"
               end if
            end if
            if (first) first = .False.
            call output_subres(unit, high_nmodes_T, high_nmodes_P, inv_error, &
               & iter_count)
            iter_count = iter_count + 1
         end do
         close(unit)
         call output_results(prefix, lhigh, high_nmodes_T, high_nmodes_P, & 
            & failed_on_first)
         deallocate(lmax2lhigh)
      end do

   end subroutine run_optimization_search

   subroutine get_projection_operators(eigvecs, tmodes, pmodes, n_t, noise, & 
         & noise_transf, ylm, pylm)
      implicit none

      real(dp), dimension(:, :), intent(in)     :: eigvecs, noise, ylm
      real(dp), dimension(:, :), intent(out)    :: noise_transf, pylm
      integer(i4b), intent(in)      :: tmodes, pmodes, n_t

      real(dp), allocatable, dimension(:, :) :: projmat, buffer
      integer(i4b)      :: nval, nmodes, n_h

      nval = size(ylm, 1)
      n_h = size(ylm, 2)
      nmodes = tmodes+pmodes

      allocate(projmat(nmodes, nval))

      projmat = 0

      projmat(1:tmodes, :) = transpose(eigvecs(:, n_t-tmodes+1:n_t))
      projmat(tmodes+1:, :) = transpose(eigvecs(:, nval-pmodes+1:nval))

      allocate(buffer(nval, nmodes))
      !Form PCP^T
      call dgemm('N', 'T', nval, nmodes, nval, 1.d0, noise, nval, projmat, & 
         & nmodes, 0.d0, buffer, nval)
      call dgemm('N', 'N', nmodes, nmodes, nval, 1.d0, projmat, nmodes, &
         & buffer, nval, 0.d0, noise_transf, nmodes)
      deallocate(buffer)

      !Form PY
      call dgemm('N', 'N', nmodes, n_h, nval, 1.d0, projmat, nmodes, ylm, &
         & nval, 0.d0, pylm, nmodes)

      deallocate(projmat)

   end subroutine get_projection_operators

   subroutine open_subresfile(unit, prefix, lhigh)
      implicit none

      integer(i4b), intent(out)      :: unit
      integer(i4b), intent(in)       :: lhigh   
      character(len=*)  ::      prefix

      character(len=512)        :: fname
      character(len=4)  :: lhighstring

      call comm_int2string(lhigh, lhighstring)
      fname = trim(prefix) // '_' // lhighstring // '_subres.dat'
      unit = comm_getlun()
      open(unit, file=trim(fname))

   end subroutine open_subresfile

   subroutine output_subres(unit, tmodes, pmodes, inv_error, iter_count)
      implicit none
      integer(i4b), intent(in)      :: unit, tmodes, pmodes, iter_count
      logical(lgt), intent(in)      :: inv_error

      write(unit, "(I6, I6, I6, L10)") iter_count, tmodes, pmodes, inv_error

   end subroutine output_subres

   subroutine output_results(prefix, lhigh, tmodes, pmodes, failed_on_first)
      implicit none

      character(len=*), intent(in)      :: prefix
      integer(i4b), intent(in)  :: lhigh, tmodes, pmodes
      logical(lgt), intent(in)  :: failed_on_first

      integer(i4b)      :: unit 
      character(len=512)        :: fname
      character(len=4)        :: lhighstring

      call comm_int2string(lhigh, lhighstring)

      fname = trim(prefix) // '_lhigh_' // lhighstring // '.dat'
      unit = comm_getlun()
      open(unit, file=trim(fname))
      if (failed_on_first) then
         write(unit, "(I6, I6, L10)") 0, 0, failed_on_first
      else
         write(unit, "(I6, I6, L10)") tmodes, pmodes, failed_on_first
      end if
      close(unit)

   end subroutine output_results

   subroutine output_invfisher_diags(prefix, lhigh, bins, invF, iter_count)
      implicit none
      character(len=*), intent(in)      :: prefix
      integer(i4b), intent(in)  :: lhigh, iter_count
      integer(i4b), dimension(:, :), intent(in)     :: bins
      real(dp), dimension(:, :), intent(in)     :: invF

      integer(i4b)      :: nbins, i, unit
      character(len=5)  :: lhighstring
      character(len=5)  ::      iterstring
      character(len=512)        :: fname

      nbins = size(invF, 1)
      unit = comm_getlun()

      call comm_int2string(lhigh, lhighstring)
      call comm_int2string(iter_count, iterstring)
      fname = trim(prefix) // '_lhigh_' // lhighstring // '_iter_' // iterstring //  '_invF.dat'

      open(unit, file=trim(fname))
      do i = 1, nbins
         write(unit, "(I6, I6, I6, I6, E15.7)"), i, bins(i, 1), bins(i, 2), bins(i, 3), sqrt(invF(i, i))
      end do
      close(unit)

   end subroutine output_invfisher_diags
   
end module comm_like_optimization_mod
