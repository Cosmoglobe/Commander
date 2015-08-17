program comm_like_tools
  use comm_proc_utils
  use comm_br_mod
  use comm_lowl_mod
  use quasi_newton_mod
  use comm_like_optimization_mod
  use powell_mod
  use spline_1D_mod
  use sort_utils
  use comm_like_timing_mod
#ifdef WMAP
  use wmap_likelihood_9yr
  use wmap_options
#endif
  implicit none

  include "mpif.h"

  character(len=128) :: operation
  integer(i4b)       :: unit

  ! Global parameters for power spectrum search
  real(dp)                                     :: chisq_powspec
  integer(i4b)                                 :: npar_powspec
  integer(i4b),              dimension(1000,3) :: bins_powspec
  real(dp),     allocatable, dimension(:,:)    :: cls_fid

  integer(i4b) :: ierr, myid, numprocs, root

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
  root         = 0

  if (iargc() == 0) then
     write(*,*) 'Usage: comm_like_tools [operation] [args] ... '
     write(*,*) '  Valid operations are:'
     write(*,*) '     mapcov2gausslike     -- '
     write(*,*) '     sigma2cl_BR          -- '
     write(*,*) '     sigma2slice_BR       -- '
     write(*,*) '     sigma2par_BR         -- '
     write(*,*) '     sigma2tau_BR         -- '
     write(*,*) '     sigma2qn_BR          -- '
     write(*,*) '     slice_gauss          -- '
     write(*,*) '     powspec_gauss        -- '
     write(*,*) '     optimization_search  -- '
     write(*,*) '     print_fisher_mat     -- '
     write(*,*) '     timing               -- '
     write(*,*) '  For further usage information, type "comm_like_tools [operation]"'
     stop
  end if

  unit = 19
  call getarg(1,operation)
  if (trim(operation) == 'mapcov2gausslike') then
     if (myid == root) call mapcov2gausslike
!!$  else if (trim(operation) == 'sigma2cl_BR') then
!!$     call sigma2cl_BR
  else if (trim(operation) == 'sigma2slice_BR') then
     if (myid == root) call sigma2slice_BR
  else if (trim(operation) == 'slice_gauss') then
     if (myid == root) call slice_gauss
!!$  else if (trim(operation) == 'sigma2par_BR') then
!!$     call sigma2par_BR('par')
!!$  else if (trim(operation) == 'sigma2tau_BR') then
!!$     call sigma2par_BR('tau')
  else if (trim(operation) == 'sigma2qn_BR') then
     if (myid == root) call sigma2qn_BR
  else if (trim(operation) == 'powspec_gauss') then
     call powspec_gauss
  else if (trim(operation) == 'optimization_search') then
     if (myid == root) call run_optimization_search
  else if (trim(operation) == 'print_fisher_mat') then
     if (myid == root) call print_fisher_mat
  else if (trim(operation) == 'timing') then
     if (myid == root) call run_timing
  end if

  call mpi_finalize(ierr)

contains

  subroutine mapcov2gausslike
    implicit none

    integer(i4b)       :: i, j, k, l, m, q, nside, nmaps, npix, lmax, lhigh_T, lhigh_P, ordering, numcomp
    integer(i4b)       :: l1, l2, m1, m2, i1, i2, n_d, n_p, n_h, n
    integer(i4b)       :: ind, ind1, ind2, cov_order, seed, numtemp, nmask, n_t, n_pol
    integer(i4b)       :: lmax_basis, pol, nval, nhigh, num_fg_temp, nside_temp, npix_temp
    real(dp)           :: t1, t2, Tthreshold, Pthreshold, reg(3), scale, W_max_T, W_max_P, BBamp, cvec(3)
    character(len=5)   :: itext
    character(len=512) :: mapfile, maskfile, covfile, beamfile, outfile, filename, temp
    character(len=512) :: basis, prefix, clfile
    real(dp),     allocatable, dimension(:,:)     :: cov, P_harm, C, cl, B, invC
    real(dp),     allocatable, dimension(:)       :: W, mask_1d, test3, A_T
    real(dp),     allocatable, dimension(:,:)     :: d, mask, map, alms, V, P_l, V_red, S_mat, P, map_1d
    real(dp),     allocatable, dimension(:,:)     :: Y, Yp, beam, cls, T, buffer, buffer2, P_highl, weights
    real(dp),     allocatable, dimension(:,:)     :: PT, PT_invN_TP
    complex(dpc), allocatable, dimension(:,:,:)   :: alms_cmplx
    real(dp),     pointer,     dimension(:,:)     :: pixwin
    real(dp),     allocatable, dimension(:,:,:)   :: window
    integer(i4b), allocatable, dimension(:)       :: indmap, map2mask, lmax2lhigh
    logical(lgt), allocatable, dimension(:)       :: tempmode
    type(planck_rng) :: handle


    if (iargc() < 17) then
       write(*,*) '    Pre-process low-l likelihood inputs from map and covariance matrix'
       write(*,*) '    Options:  [mapfile] [maskfile] [covfile] [beamfile] [clfile] [lmax_data] [lmax_basis_T] [lmax_basis_P]'
       write(*,*) '                 [basis] [T threshold] [P threshold] [BB amp] [T reg] [P reg] [seed] [outprefix] [template1...]'
       write(*,*) '        where [basis]     = {pixel, pseudoalm, eigen_N, eigen_StoN, eigen_S+N}'
       write(*,*) '              [threshold] = basis dependent accuracy threshold'
       stop
    end if

    unit = 58
    call getarg(2,mapfile)
    call getarg(3,maskfile)
    call getarg(4,covfile)
    call getarg(5,beamfile)
    call getarg(6,clfile)
    call getarg(7,temp)
    read(temp,*) lmax
    call getarg(8,temp)
    read(temp,*) lhigh_T
    call getarg(9,temp)
    read(temp,*) lhigh_P
    call getarg(10,basis)
    call getarg(11,temp)
    read(temp,*) Tthreshold
    call getarg(12,temp)
    read(temp,*) Pthreshold
    call getarg(13,temp)
    read(temp,*) BBamp
    call getarg(14,temp)
    read(temp,*) reg(1)
    call getarg(15,temp)
    read(temp,*) reg(2)
    reg(3) = reg(2)
    call getarg(16,temp)
    read(temp,*) seed
    call getarg(17,prefix)
    numcomp = (lmax+1)**2

    call rand_init(handle, seed)

    write(*,*) 'Reading data:'
    write(*,*) '     Map file   = ', trim(mapfile)
    write(*,*) '     Mask       = ', trim(maskfile)
    write(*,*) '     Covariance = ', trim(covfile)
    write(*,*) '     Cl file    = ', trim(clfile)

    ! Read fiducial spectrum; only actually used for S/N bases
    call read_fiducial_spectrum(clfile, cls)

    ! Read map and mask
    j = len(trim(mapfile))
    if (mapfile(j-2:j) == 'txt') then
       open(unit,file=trim(mapfile))
       read(unit,*) n_d
       read(unit,*) mapfile
       backspace(unit)
       i    = getsize_fits(trim(mapfile), ordering=ordering, nside=nside, nmaps=nmaps)
       npix = 12*nside**2
       n_p  = npix    * nmaps
       n_h  = numcomp * nmaps
       allocate(map(0:npix-1,nmaps), map_1d(n_p,n_d))
       do i = 1, n_d
          read(unit,*) mapfile
          call read_map(mapfile,  map)
          do j = 1, nmaps
             map_1d((j-1)*npix+1:j*npix,i) = map(:,j)
          end do
       end do
    else
       i    = getsize_fits(trim(mapfile), ordering=ordering, nside=nside, nmaps=nmaps)
       npix = 12*nside**2
       n_d  = 1
       n_p  = npix    * nmaps
       n_h  = numcomp * nmaps
       allocate(map(0:npix-1,nmaps), map_1d(n_p,n_d))
       call read_map(mapfile,  map)
       do j = 1, nmaps
          map_1d((j-1)*npix+1:j*npix,1) = map(:,j)
       end do
    end if
    allocate(mask(0:npix-1,nmaps), mask_1d(n_p))
    call read_map(maskfile, mask)
    where (mask < 0.5d0)
       mask = 0.d0
    elsewhere
       mask = 1.d0
    end where

!!$    allocate(indmap(0:npix-1))
!!$    call ang2vec(0.5d0*pi, pi, cvec)
!!$    call query_disc(nside, cvec, 90.d0*pi/180.d0, indmap, i)
!!$    mask(indmap(0:i-1),2) = 0.d0
!!$    deallocate(indmap)

    do i = 1, nmaps
       mask_1d((i-1)*npix+1:i*npix) = mask(:,i)
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

    ! Set up lmax to lhigh conversion array
    if (nmaps == 3) then
      nhigh = (lhigh_T+1)**2 + 2*(lhigh_P+1) ** 2 - 8
   else
      nhigh = (lhigh_T+1) ** 2
   end if
!    if (nmaps == 3) nhigh = nhigh - 8 ! No polarization monopoles or dipoles
    allocate(lmax2lhigh(nhigh))
    k = 1
    j = 1
    do i = 1, nmaps
       do l = 0, lmax
          do m = -l, l
!             if (l <= lhigh .and. (i == 1 .or. l>=2)) then
             if ((l <= lhigh_T .and. i == 1) .or. (l>=2 .and. l <= lhigh_P .and. i > 1)) then
                lmax2lhigh(k) = j
                k             = k+1
             end if
             j = j+1
          end do
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
    cov(1:npix,:) = 0.d0
    cov(:,1:npix) = 0.d0

    ! Add regularization noise to covariance matrix
    ind = 1
    do j = 1, nmaps
       do i = 1, npix
          cov(ind,ind) = cov(ind,ind) + reg(j)**2
          do k = 1, n_d
             map_1d(ind,k)  = map_1d(ind,k)  + reg(j) * rand_gauss(handle)
          end do
          ind          = ind+1
       end do
    end do

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

    ind1 = 1
    allocate(buffer(numcomp,nmaps))
    do j = 1, nmaps
       do i = 0, npix-1
          map      = 0.d0
          map(i,j) = 1.d0
          if (nmaps == 1) then
             call map2alm(nside, lmax, lmax, map(:,1), alms_cmplx, [0.d0,0.d0], weights)
          else
             call map2alm(nside, lmax, lmax, map,      alms_cmplx, [0.d0,0.d0], weights)
          end if
          call convert_complex_to_real_alms(alms_cmplx, buffer)
          do k = 1, nmaps
             Yp((k-1)*numcomp+1:k*numcomp,ind1) = buffer(:,k)
          end do
          ind1 = ind1+1
       end do
    end do
    deallocate(buffer)

    ! Set up monopole and dipole templates if n_t > 0, and read foreground templates
    nside_temp  = 0
    if (nside_temp > 0) then
       npix_temp   = 12*nside_temp**2
    else
       npix_temp   = 1
    end if
    q           = npix/npix_temp
    num_fg_temp = iargc()-17
    numtemp     = num_fg_temp * npix_temp
    if (n_t > 0) numtemp = numtemp+4 
    allocate(T(n_p,numtemp))
    T = 0.d0
    if (num_fg_temp > 0) then
       do k = 1, num_fg_temp
          call getarg(17+k,filename)
          call read_map(filename, map)
          do i = 1, nmaps
             call convert_ring2nest(nside, map(:,i))
!             do j = (k-1)*npix_temp+1, k*npix_temp
             do j = 1, npix_temp
                !T((i-1)*npix+1:i*npix,j) = map(:,i) 
                T((i-1)*npix+(j-1)*q+1:(i-1)*npix+j*q,(k-1)*npix_temp+j) = map((j-1)*q:j*q-1,i) 
                call convert_nest2ring(nside, T((i-1)*npix+1:i*npix,(k-1)*npix_temp+j))
             end do
          end do
       end do
       T   = T * scale ! Assume templates are in same units as covariance
    end if
    if (n_t > 0) then
       T(1:npix,num_fg_temp+1) = 1.d0
       do i = 0, npix-1
          call pix2vec_ring(nside,i,T(i+1,num_fg_temp+2:num_fg_temp+4))
       end do
       !T(1:npix,num_fg_temp+1:num_fg_temp+4) = scale*T(1:npix,num_fg_temp+1:num_fg_temp+4) 
    end if

    ! Subtract best-fit templates from data, and add uncertainty to covariance matrix
    if (numtemp > 0) then
       allocate(invC(nval,nval), PT(nval,numtemp), PT_invN_TP(numtemp,numtemp), A_T(numtemp))
       PT = T(map2mask,:)
!       PT_invN_TP = matmul(transpose(PT), PT)
!       call invert_matrix(PT_invN_TP)

       allocate(S_mat(nval,nval))
       call get_basis_S(0.d0, 1.d0, cls, beam, Y(map2mask,:), S_mat)
       S_mat = S_mat + cov(map2mask,map2mask)
       call invert_matrix(S_mat)
       S_mat(1:n_t,n_t+1:nval) = 0.d0
       S_mat(n_t+1:nval,1:n_t) = 0.d0
       PT = T(map2mask,:)
       PT_invN_TP = matmul(transpose(PT), matmul(S_mat,PT))
       call invert_singular_matrix(PT_invN_TP, 1.d-12)

       do k = 1, n_d
          A_T = matmul(PT_invN_TP, matmul(transpose(PT), matmul(S_mat,map_1d(map2mask,k))))
       
          write(*,*) 'Subtracted template amplitudes for map ', k
          do i = 1, numtemp
             write(*,*) '   Template ', i, ' = ', A_T(i)
          end do
       
!!$       map_1d = -1.6375d30
!!$       map_1d(map2mask) = d
!!$       do j = 1, nmaps
!!$          map(:,j) = map_1D((j-1)*npix+1:j*npix) * mask_1d((j-1)*npix+1:j*npix)
!!$       end do
!!$       call write_map3('data_foer.fits', map)
!!$       write(*,*) 'foer = ', sum(map_1D(map2mask(1:n_t)))/n_t

          if (k == 1) then
             do j = 1, nmaps
                map(:,j) = map_1d((j-1)*npix+1:j*npix,k) * mask_1d((j-1)*npix+1:j*npix)
             end do
             call write_map3(trim(prefix)//'_rawmap.fits', map)
          end if

          map_1d(map2mask,k) = map_1d(map2mask,k) - matmul(PT, A_T)

          if (k == 1) then
             do j = 1, nmaps
                map(:,j) = map_1d((j-1)*npix+1:j*npix,k) * mask_1d((j-1)*npix+1:j*npix)
             end do
             call write_map3(trim(prefix)//'_cleanmap.fits', map)
          end if
       end do
       deallocate(S_mat)

       cov(map2mask,map2mask) = cov(map2mask,map2mask) + 1.d3 * matmul(PT, transpose(PT))
    end if
       


    if (nmaps == 1) then
       write(*,fmt='(a,f8.3)') '      I sky fractions = ', sum(mask(:,1))/npix
    else
       write(*,fmt='(a,3f8.3)') '      I/Q/U sky fractions = ', sum(mask(:,1))/npix, &
            & sum(mask(:,2))/npix, sum(mask(:,3))/npix
    end if


    ! Compute high-l projection operator
    allocate(P_highl(nval,nval))
!    call dgemm('N','T',nval,nval,nhigh,1.d0,Y(map2mask,lmax2lhigh),nval,&
!         & Y(map2mask,lmax2lhigh),nval,0.d0,P_highl,nval)
    call dgemm('N','N',nval,nval,nhigh,1.d0,Y(map2mask,lmax2lhigh),nval,&
         & Yp(lmax2lhigh,map2mask),nhigh,0.d0,P_highl,nval)

    write(*,*) 'Setting up basis set'
    allocate(P(nval,nval))
    if (trim(basis) == 'pixel') then

       P = 0.d0
       do i = 1, nval
          P(i,i) = 1.d0
       end do
       
    else if (trim(basis) == 'pseudoalm') then

       P = P_highl

    else if (trim(basis) == 'eigen_N') then

       P = cov(map2mask,map2mask)
       call invert_matrix(P)

    else if (trim(basis) == 'eigen_StoN') then

       P = cov(map2mask,map2mask)
       call invert_matrix(P)
       allocate(buffer(nval,nval), S_mat(nval,nval))
       ! Compute S^1/2 * invN * S^1/2
       call get_basis_S(BBamp, 0.5d0, cls, beam, Y(map2mask,:), S_mat)
       call dgemm('N','N',nval,nval,nval,1.d0,S_mat,nval,P,nval,0.d0,buffer,nval)
       call dgemm('N','N',nval,nval,nval,1.d0,buffer,nval,S_mat,nval,0.d0,P,nval)
       deallocate(buffer, S_mat)
       
    else if (trim(basis)=='eigen_S+N') then

       allocate(S_mat(nval,nval))
       call get_basis_S(BBamp, 1.d0, cls, beam, Y(map2mask,:), S_mat)
       P = cov(map2mask,map2mask) + S_mat
       deallocate(S_mat)

    else
       write(*,*) 'Unsupported basis type = ', trim(basis)
       stop
    end if

    ! Project out high-l modes by P_highl * invN * P_highl^t
    if (min(lhigh_T, lhigh_P) < lmax .and. (trim(basis) == 'eigen_N' .or. trim(basis) == 'eigen_StoN' .or. &
         & trim(basis) == 'eigen_S+N')) then
       write(*,*) 'Projecting out high-l modes'
       allocate(buffer(nval,nval))
       call dgemm('N','N',nval,nval,nval,1.d0,P_highl,nval,P,nval,0.d0,buffer,nval)
       call dgemm('N','T',nval,nval,nval,1.d0,buffer,nval,P_highl,nval,0.d0,P,nval)
       deallocate(buffer)
    end if
    
    ! Nullify temperature-polarization cross-elements
    P(1:n_t,n_t+1:nval) = 0.d0
    P(n_t+1:nval,1:n_t) = 0.d0

    ! Compute projection eigen-decomposition
    allocate(B(nval,nval), W(nval), V(nval,nval))
    if (trim(basis) == 'pixel') then
       V = P
       W = 1.d0
    else
       write(*,*) 'Computing eigen vectors for ', trim(basis)
       call get_eigen_decomposition(P, W, V)
    end if
    
    ! Find max temperature and polarization eigenvalues
    allocate(tempmode(nval))
    W_max_T = -1.d30; W_max_P = -1.d30
    do i = 1, nval
       if (n_t > 0 .and. sum(V(1:n_t,i)**2) > sum(V(n_t+1:nval,i)**2)) then
          W_max_T = max(W_max_T, W(i))
          tempmode(i) = .true.
       else
          W_max_P = max(W_max_P, W(i))
          tempmode(i) = .false.
       end if
    end do
       
    ! Print eigen spectrum, and find number of acceptable modes
    outfile = trim(prefix) // '_W.dat'
    write(*,*) '   Writing basis eigen spectrum to ', trim(outfile)
    open(unit,file=trim(outfile))
    ! First temperature
    n = 0
    j = 1
    do i = 1, nval
       if (tempmode(i)) then
          j = j+1
          write(unit,*) j, abs(W(i)/maxval(W)), W(i)
          if (W(i) > Tthreshold * W_max_T) n = n+1
       end if
    end do
    k = n
    write(*,*) '     Number of accepted temperature modes = ', n
    write(unit,*)
    ! Then polarization
    j = 0
    do i = 1, nval
       if (.not. tempmode(i)) then
          j = j+1
          write(unit,*) j, abs(W(i)/maxval(W)), W(i)
          if (W(i) > Pthreshold * W_max_P) n = n+1
       end if
    end do
    write(*,*) '     Number of accepted polarization modes = ', n-k
    close(unit)
          
    ! Select acceptable eigenvectors
    allocate(V_red(nval,n))
    j = 1
    do i = 1, nval
       if ((tempmode(i) .and. W(i) > Tthreshold * W_max_T) .or. &
            & (.not. tempmode(i) .and. W(i) > Pthreshold * W_max_P)) then
          V_red(:,j) = V(:,i)
          j          = j+1
       end if
    end do

    ! Transform data and covariance to new basis, and set up harmonic space projection operator
    allocate(d(n,n_d), C(n,n), P_harm(n_h,n), buffer(nval,n))
    !d = matmul(transpose(V_red), map_1d(map2mask))
    call dgemm('T','N',n,n_d,nval,1.d0,V_red,nval,map_1d(map2mask,:),nval,0.d0,d,n)
    ! Compute C = V^t * C * V
    call dgemm('N','N',nval,n,nval,1.d0,cov(map2mask,map2mask),nval,V_red,nval,0.d0,buffer,nval)
    call dgemm('T','N',n,n,nval,1.d0,V_red,nval,buffer,nval,0.d0,C,n)
    ! Compute projection vectors, P_harm = Y^t * V
    call dgemm('T','N',n_h,n,nval,1.d0,Y(map2mask,:),nval,V_red,nval,0.d0,P_harm,n_h)
    deallocate(buffer)

    ! Output data objects
    write(*,*) '   Number of sky maps      = ', n_d
    write(*,*) '   Number of pixels        = ', n_p
    write(*,*) '   Number of harmonics     = ', n_h
    write(*,*) '   Number of basis vectors = ', n

    ! Output number of modes to disk
    open(68,file='nmodes.dat')
    write(68,*) n
    close(68)

    ! Output datafile
    call write_lowl_datafile(trim(prefix)//'.fits', d, C, beam, P_harm, mapfile, maskfile, &
         & covfile, beamfile, clfile, basis, lhigh_T, lhigh_P, Tthreshold, Pthreshold)

  end subroutine mapcov2gausslike


!!$  subroutine sigma2cl_BR
!!$    implicit none
!!$
!!$    character(len=256) :: sigmafile, clfile, binfile, outfile, temp, line
!!$    real(dp)           :: cl_min, cl_max, sigma, top, left, right, int, tot_int, cl_in(4), dx
!!$    real(dp)           :: peak, lower, upper
!!$    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep
!!$    integer(i4b)       :: numbin, pos(3), ell, b_min, b_max, bincount
!!$    integer(i4b)       :: nspec, lmax, numchain, numiter, i, j, k, l, m, n, p, bin(2), spec, ierr
!!$    logical(lgt)       :: exist, marginal
!!$    character(len=1), dimension(6)   :: status
!!$    integer(i4b), allocatable, dimension(:) :: include_chain
!!$    real(sp), pointer, dimension(:,:,:,:) :: cls_pointer
!!$    real(dp), allocatable, dimension(:)   :: vals, cl, vals2, cl_hires, vals_hires, sigmas
!!$    real(dp), allocatable, dimension(:,:) :: lnL, cls, cls0
!!$
!!$    if (iargc() < 12) then
!!$       write(*,*) '    Compute spectrum with Blackwell-Rao'
!!$       write(*,*) '      Options: [sigmafile] [clfile] [binfile] [firstchain] [lastchain] '
!!$       write(*,*) '               [firstsample] [lastsample] [thinstep] [spec] [marginal] [outfile] [chains to remove]'
!!$       stop
!!$    end if
!!$
!!$    call getarg(2, sigmafile)
!!$    call getarg(3, clfile)
!!$    call getarg(4, binfile)
!!$    call getarg(5, temp)
!!$    read(temp,*) firstchain
!!$    call getarg(6, temp)
!!$    read(temp,*) lastchain
!!$    call getarg(7, temp)
!!$    read(temp,*) firstsample
!!$    call getarg(8, temp)
!!$    read(temp,*) lastsample
!!$    call getarg(9, temp)
!!$    read(temp,*) thinstep
!!$    call getarg(10, temp)
!!$    read(temp,*) spec
!!$    call getarg(11, temp)
!!$    read(temp,*) marginal
!!$    call getarg(12, outfile)
!!$
!!$    allocate(include_chain(firstchain:lastchain))
!!$    include_chain = .true.
!!$    do i = 1, iargc()-12
!!$       call getarg(12+i, temp)
!!$       read(temp,*) j
!!$       include_chain(j) = .false.
!!$    end do
!!$
!!$    numbin = 10000
!!$    if (spec == 1) then
!!$       cl_min = 1.d1
!!$       cl_max = 1.d4
!!$    else if (spec == 2) then
!!$       cl_min = -1.d3
!!$       cl_max =  1.d3
!!$    else if (spec == 4) then
!!$       cl_min = 1.d-5
!!$       cl_max = 1.d2
!!$    else if (spec == 6) then
!!$       cl_min = 1.d-5
!!$       cl_max = 1.d2
!!$    else
!!$       write(*,*) 'Spectrum not yet supported = ', spec
!!$       stop
!!$    end if
!!$
!!$    inquire(file=trim(binfile), exist=exist)
!!$    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(binfile)
!!$    inquire(file=trim(sigmafile), exist=exist)
!!$    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(sigmafile)
!!$    inquire(file=trim(clfile), exist=exist)
!!$    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(clfile)
!!$
!!$    ! Read fiducial spectrum file
!!$    unit = comm_br_getlun()
!!$    open(unit,file=trim(clfile))
!!$    lmax = -1
!!$    do while (.true.)
!!$       read(unit,'(a)',end=79) line
!!$       line = trim(line)
!!$       if (line(1:1) == '#') cycle
!!$       read(line,*) l
!!$       lmax = max(lmax, l)
!!$    end do
!!$79  close(unit)
!!$    if (lmax > -1) then
!!$       allocate(cls0(0:lmax,6), cls(0:lmax,6))
!!$       cls0 = 0.d0
!!$       open(unit,file=trim(clfile))
!!$       do while (.true.)
!!$          read(unit,'(a)',end=80) line
!!$          line = trim(line)
!!$          if (line(1:1) == '#') cycle
!!$          read(line,*) l, cl_in
!!$          cls0(l,1) = cl_in(1) ! Assume (TT, EE, BB, TE) ordering for now
!!$          cls0(l,2) = cl_in(4)          
!!$          cls0(l,4) = cl_in(2)
!!$          cls0(l,6) = cl_in(3)
!!$       end do
!!$80     close(unit)
!!$    else
!!$       write(*,*) 'Error -- comm_br_mod: No valid entries in clfile = ', trim(clfile)
!!$       stop
!!$    end if
!!$
!!$    allocate(cl(numbin), vals(numbin), vals2(numbin))
!!$    allocate(cl_hires(numbin), vals_hires(numbin))
!!$
!!$    open(59,file=trim(outfile), recl=1024)
!!$    open(60,file='slices_'//trim(outfile), recl=1024)
!!$    write(59,*) '#  l_center  l_min     l_max        D_l         upper error     lower error     mean error      no-CV error'    
!!$
!!$    open(58,file=trim(binfile))
!!$    bincount = 0
!!$    do while (.true.)
!!$       read(58,'(A)',end=341) line
!!$       if (line(1:1) == '#' .or. trim(line)=='') cycle
!!$       read(line,*) bin, status
!!$       if (status(spec) == '0') cycle
!!$       do j = 1, size(status)
!!$          if (j /= spec) status(j) = '0'
!!$       end do
!!$       bincount = bincount + 1
!!$       !if (bin(1) /= 234) cycle
!!$       if (marginal) then
!!$          call comm_br_initialize_object(sigmafile, clfile, firstchain, lastchain, &
!!$               & firstsample, lastsample, thinstep, binfile, binnum=bincount, spec_marg=spec, &
!!$               & include_chain=include_chain)
!!$       else
!!$          call comm_br_initialize_object(sigmafile, clfile, firstchain, lastchain, &
!!$               & firstsample, lastsample, thinstep, binfile, binnum=bincount, include_chain=include_chain)
!!$       end if
!!$
!!$       cls = cls0
!!$       do i = 1, numbin
!!$          cl(i)                   = cl_min + (i-0.5) * (cl_max - cl_min) / (numbin-1)
!!$          cls(bin(1):bin(2),spec) = cl(i)
!!$          vals(i)                 = comm_br_compute_lnL(cls, ierr)
!!$       end do
!!$       vals = exp(vals - maxval(vals))
!!$       vals = vals / (sum(vals) * (cl_max-cl_min)/(numbin-1))
!!$       call spline(cl, vals, 1.d30, 1.d30, vals2)
!!$
!!$       if (.false.) then
!!$          open(99,file='slice.dat')
!!$          do i = 1, numbin
!!$             write(99,*) cl(i), vals(i)
!!$          end do
!!$          close(99)
!!$          stop
!!$       end if
!!$
!!$       open(99,file='myslices.dat')
!!$       do i = 1, size(lnL,1)
!!$          do p = 1, numbin
!!$             if (lnL(i,p) > 0.d0) then
!!$                write(99,*) cl(p), lnL(i,p)
!!$             end if
!!$          end do
!!$          write(99,*)
!!$       end do
!!$       close(99)
!!$       stop
!!$
!!$       open(99,file='sigmas.dat')
!!$       do i = 1, size(sigmas)
!!$          write(99,*) sigmas(i)
!!$       end do
!!$       close(99)
!!$       stop
!!$       
!!$       ! Find maximum
!!$       pos(1:1) = maxloc(vals)
!!$       pos(2) = 1
!!$       do while (vals(pos(2))/vals(pos(1)) < 1.d-6 .and. pos(2) < pos(1))
!!$          pos(2) = pos(2) + 1
!!$       end do
!!$       pos(3) = numbin
!!$       do while (vals(pos(3))/vals(pos(1)) < 1.d-6 .and. pos(3) > pos(1))
!!$          pos(3) = pos(3) - 1
!!$       end do
!!$       pos(2) = max(pos(2)-3,1)
!!$       pos(3) = min(pos(3)+3,numbin)
!!$
!!$       ! Make a high-resolution grid
!!$       cls = cls0
!!$       do i = 1, numbin
!!$          cl_hires(i)             = cl(pos(2)) + (i-1) * (cl(pos(3))-cl(pos(2))) / (numbin-1.d0)
!!$          cls(bin(1):bin(2),spec) = cl_hires(i)
!!$          vals_hires(i)           = comm_br_compute_lnL(cls, ierr)
!!$!          vals_hires(i) = splint(cl, vals, vals2, cl_hires(i))
!!$       end do
!!$       vals_hires = exp(vals_hires - maxval(vals_hires))
!!$       vals_hires = vals_hires / (sum(vals_hires) * (cl(pos(3))-cl(pos(2)))/(numbin-1))
!!$
!!$       if (.true.) then
!!$          open(99,file='slice.dat')
!!$          do i = 1, numbin
!!$             write(99,*) cl_hires(i), vals_hires(i)
!!$          end do
!!$          close(99)
!!$          stop
!!$       end if
!!$       !do i = 1, numbin
!!$       !   write(60,*) cl_hires(i), vals_hires(i)
!!$       !end do
!!$       !write(60,*)
!!$       pos(1:1) = maxloc(vals_hires)
!!$       pos(2:3) = pos(1)
!!$       tot_int  = sum(vals_hires)
!!$       int      = 0.d0
!!$       do while (int < 0.68*tot_int)
!!$          if ((vals_hires(pos(2)) > vals_hires(pos(3)) .and. pos(2) > 1) .or. pos(3) == numbin) then
!!$             pos(2) = pos(2)-1
!!$             int    = int + 0.5d0 * (vals_hires(pos(2))+vals_hires(pos(2)+1))
!!$          else
!!$             pos(3) = pos(3)+1
!!$             int    = int + 0.5d0 * (vals_hires(pos(3))+vals_hires(pos(3)-1))
!!$          end if
!!$       end do
!!$       pos(3) = max(pos(3), pos(2)+1)
!!$
!!$       write(*,*) pos
!!$       dx    = cl(2)-cl(1)
!!$       peak  = cl_hires(pos(1))
!!$       lower = cl_hires(pos(1))-cl_hires(pos(2))+0.5d0*dx
!!$       upper = cl_hires(pos(3))-cl_hires(pos(1))+0.5d0*dx
!!$       if (pos(1) == 1) then
!!$          peak = 0.d0
!!$          lower = 0.d0
!!$       end if
!!$
!!$
!!$       b_min = bin(1); b_max = bin(2)
!!$       write(*,fmt='(f8.1,2i10,3f16.8)') (b_min+b_max)/2., b_min, b_max, peak, upper, lower
!!$       write(59,fmt='(f8.1,2i10,5f16.8)') (b_min+b_max)/2., b_min, b_max, peak, upper, lower, 0.5d0 * (upper+lower), sqrt(variance(sigmas))
!!$
!!$       call comm_br_deallocate_object
!!$
!!$    end do
!!$341 close(58)
!!$    close(59)
!!$    close(60)
!!$
!!$    deallocate(vals, vals2)
!!$
!!$  end subroutine sigma2cl_BR

  subroutine sigma2slice_BR
    implicit none

    character(len=256) :: infile, binfile, outfile, temp, line, clfile, outprefix
    character(len=2)   :: spectext, ptext
    character(len=4)   :: lmin_text, lmax_text
    real(dp)           :: cl_min, cl_max, sigma, top, left, right, int, tot_int, cl_in(4), lnL_max
    integer(i4b)       :: nspec, lmin, lmax, numchain, numiter, i, j, k, l, m, n, p, bins(4), id, numbin, ind
    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep, nmaps, ierr
    logical(lgt)       :: exist
    real(dp), allocatable, dimension(:,:) :: cls0, cls
    real(dp), allocatable, dimension(:)   :: cl, lnL
    character(len=2), dimension(6) :: stext = ['TT','TE','TB','EE','EB','BB']
    real(dp), dimension(6) :: clmin = [0.d0, -2.d1, -1.d0, 0.d0, -1.d0, 0.d0]
    real(dp), dimension(6) :: clmax = [1.d4,  2.d1,  1.d0, 1.d0,  1.d0, 0.01d0]

    if (iargc() /= 4) then
       write(*,*) '    Compute Blackwell-Rao slices'
       write(*,*) '      Usage: comm_like_tools sigma2slice_BR [BR file] [clfile] [outprefix]'
       stop
    end if

    call getarg(2, infile)
    call getarg(3, clfile)
    call getarg(4, outprefix)

    inquire(file=trim(infile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(infile)
    inquire(file=trim(clfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(clfile)

    call read_fiducial_spectrum(clfile, cls0)
    lmax = size(cls0,1)-1
    allocate(cls(0:lmax,6))

    ! Initialize BR module
    call comm_br_initialize_object(infile, handle=id)

    ! Loop over all free parameters, and output slices
    numbin = 10000
    allocate(cl(numbin), lnL(numbin))
    do p = 1, comm_br(id)%n_p
       nmaps = comm_br(id)%P(p)%nmaps
       do i = 1, comm_br(id)%P(p)%nfactor
          do j = 1, comm_br(id)%P(p)%factor(i)%p
             do k = j, comm_br(id)%P(p)%factor(i)%p
                if (comm_br(id)%P(p)%factor(i)%status(j,k) == 'S') then
                   ind = comm_br(id)%P(p)%factor(i)%ind(j)*(1-comm_br(id)%P(p)%factor(i)%ind(j))/2 + &
                        & (comm_br(id)%P(p)%factor(i)%ind(j)-1)*nmaps + &
                        &  comm_br(id)%P(p)%factor(i)%ind(k)
                   cl_min   = clmin(ind)
                   cl_max   = clmax(ind)
                   spectext = stext(ind)
                   lmin     = comm_br(id)%P(p)%factor(i)%lmin
                   lmax     = comm_br(id)%P(p)%factor(i)%lmax
                   
                   cls = cls0
                   lnL_max = -1.d30
                   lnL     = log(1.d-30)
                   do l = 1, numbin
                      cl(l)  = cl_min + (cl_max-cl_min)/(numbin-1) * (l-1)
                      cls(lmin:lmax,ind) = cl(l)
                      lnL(l) = comm_br_compute_lnL(cls, ierr, handle=id, P_id=p)
                      lnL_max = max(lnL_max, lnL(l))
                      !write(*,*) cl(l), lnL(l)
                      if (lnL_max-lnL(l) > 10 .and. l>=10) exit
                   end do
                   lnL = exp(lnL-maxval(lnL))
                   lnL = lnL / sum(lnL) / ((cl_max-cl_min)/(numbin-1))
                   
                   call int2string(p,    ptext)
                   call int2string(lmin, lmin_text)
                   call int2string(lmax, lmax_text)
                   outfile = trim(outprefix) // '_P' // ptext // '_' // spectext &
                        & // '_ell' // lmin_text // '_' // lmax_text // '.dat'
                   write(*,*) 'Writing ', trim(outfile)
                   open(unit,file=trim(outfile))
                   do l = 1, numbin
                      if (lnL(l) > 1.d-6 * maxval(lnL)) write(unit,*) cl(l), lnL(l)
                   end do
                   close(unit)
                end if
             end do
          end do
       end do
    end do
    
  end subroutine sigma2slice_BR


  subroutine slice_gauss
    implicit none

    character(len=256) :: infile, binfile, outfile, temp, line, clfile, outprefix
    character(len=2)   :: spectext, ptext
    character(len=4)   :: lmin_text, lmax_text
    real(dp)           :: cl_min, cl_max, sigma, top, left, right, int, tot_int, cl_in(4), lnL_max
    real(dp)           :: chisq, t1, t2, cl_maxlike
    integer(i4b)       :: nspec, lmin, lmax, numchain, numiter, g, i, j, k, l, m, n, p, b_min, b_max, id, numbin, ind
    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep, nmaps, ierr
    integer(i4b)       :: unit_summary
    logical(lgt)       :: exist
    real(dp),     allocatable, dimension(:,:) :: cls0, cls
    real(dp),     allocatable, dimension(:)   :: cl, lnL
    character(len=1), dimension(6) :: status
    character(len=2), dimension(6) :: stext = ['TT','TE','TB','EE','EB','BB']
    real(dp), dimension(6) :: clmin = [0.d0,   -10.d0, -1.d0, -0.1d0,  -1.d0,-0.02d0]
    real(dp), dimension(6) :: clmax = [5000.d0, 10.d0,  1.d0, 1.d0,   1.d0, 0.02d0]

    if (iargc() /= 5) then
       write(*,*) '    Compute Gaussian likelihood slices'
       write(*,*) '      Usage: comm_like_tools slice_gauss [BR file] [clfile] [binfile] [outprefix]'
       stop
    end if

    call getarg(2, infile)
    call getarg(3, clfile)
    call getarg(4, binfile)
    call getarg(5, outprefix)
    nspec        = 6
    numbin       = 41
    unit_summary = comm_getlun()

    inquire(file=trim(infile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(infile)
    inquire(file=trim(clfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(clfile)
    inquire(file=trim(binfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(binfile)

    call read_fiducial_spectrum(clfile, cls0)
    lmax = size(cls0,1)-1
    allocate(cls(0:lmax,6))

    ! Initialize BR module
    call comm_lowl_initialize_object(infile)

    ! Loop over all bins, and output slices
    allocate(cl(numbin), lnL(numbin))
    open(58,file=trim(binfile))
    open(unit_summary,file=trim(outprefix)//'_summary.dat', recl=1024)
    do while (.true.)
       read(58,'(A)',end=518) line
       if (line(1:1) == '#' .or. trim(line)=='') cycle
       read(line,*) lmin, lmax, status
       do ind = 1, nspec
!          if (ind /= 1) cycle
          if (status(ind) /= '0') then
             cl_min   = clmin(ind)
             cl_max   = clmax(ind)
             spectext = stext(ind)
             
             cls = cls0
             lnL_max = -1.d30
             lnL     = -1.d30
             do l = 1, numbin
                cl(l)  = cl_min + (cl_max-cl_min)/(numbin-1) * (l-1)
                cls(lmin:lmax,ind) = cl(l)
                call wall_time(t1)
                lnL(l) = comm_lowl_compute_lnL(cls=cls, ierr=ierr, red_chisq=chisq,&
                     & enforce_pos_def=.false.)
                call wall_time(t2)
                write(*,fmt='(a,f8.3,a,f12.3,a,f8.3,a,f8.3)') &
                     & 'Cl = ', cl(l), ', lnL = ', lnL(l), ', red_chisq = ', chisq, &
                     & ', time = ', t2-t1
                lnL_max = max(lnL_max, lnL(l))
                !if (lnL_max-lnL(l) > 10 .and. l>=10) exit
             end do
             lnL = exp(lnL-maxval(lnL))
             lnL = lnL / sum(lnL) / ((cl_max-cl_min)/(numbin-1))
             
             call int2string(lmin, lmin_text)
             call int2string(lmax, lmax_text)
             outfile = trim(outprefix) // '_' // spectext &
                  & // '_ell' // lmin_text // '_' // lmax_text // '.dat'
             write(*,*) 'Writing ', trim(outfile)
             open(unit,file=trim(outfile))
             do l = 1, numbin
                if (lnL(l) > 1.d-6 * maxval(lnL)) write(unit,*) cl(l), lnL(l)
             end do
             close(unit)

             call compute_conf_limits(cl, lnL, top, left, right)
             write(unit_summary,fmt='(f6.1,2i6,3f8.3,a,a)') 0.5d0*(lmin+lmax), lmin, lmax, top, left, right, '   ', spectext

          end if
          
       end do
    end do
518 close(58)
    close(unit_summary)

  end subroutine slice_gauss

  subroutine powspec_gauss
    implicit none

    character(len=256) :: infile, binfile, outfile, temp, line, clfile, outprefix, optimizer
    character(len=256) :: error_type
    character(len=2)   :: spectext, ptext, itertext
    character(len=4)   :: lmin_text, lmax_text, chain_text
    real(dp)           :: cl_min, cl_max, top, left, right, int, tot_int, cl_in(4), lnL_max, delta
    real(dp)           :: chisq, t1, t2, cl_maxlike, dC(1), alpha, alpha_prop, lnL_prop, lnL_old, xmin
    real(dp)           :: cl_prev, lnL_new(2), cl_new, eps, eps0, accept, reject
    integer(i4b)       :: nspec, lmin, lmax, numchain, numiter, g, i, j, k, l, m, n, p, p_max(1), n_lnL, nchain
    integer(i4b)       :: b_min, b_max, id, numbin, ind, ind1, ind2, n_spline, nsamp, ndraw
    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep, nmaps, ierr
    integer(i4b)       :: unit, unit_out, iter, n_QML, maxiter, seed, burnin, blocksize
    logical(lgt)       :: exist, asymmetric_errors
    real(dp),     allocatable, dimension(:,:)    :: cls, L_prop
    real(dp),     allocatable, dimension(:,:,:)  :: sigmas
    real(dp),     allocatable, dimension(:)      :: cl, dCl, dCl0, lnL, alphas, binmask, W, cl0, lnL0, cl1, lnL1, lnL2
    real(dp),     allocatable, dimension(:,:)    :: invF, invF2, sigma, lnL_cond
    real(dp),     allocatable, dimension(:,:,:)  :: samples, tot_samples
    real(dp),     allocatable, dimension(:)      :: cl_cond, p_cond, eta, currsamp
    real(dp),     allocatable, dimension(:)      :: condnum
    real(dp),     allocatable, dimension(:)      :: cl_spline, lnL_spline
    integer(i4b), allocatable, dimension(:)      :: indmap
    logical(lgt), allocatable, dimension(:)      :: converged
    character(len=1), dimension(6) :: status
    character(len=2), dimension(6) :: stext = ['TT','TE','TB','EE','EB','BB']
    type(planck_rng) :: handle

    if (iargc() < 7) then
       write(*,*) '    Compute power spectrum from Gaussian likelihood'
       write(*,*) '      Usage: comm_like_tools powspec_gauss [optimizer] [error type] '
       write(*,*) '                  [infofile] [binfile] [fid cl file] [outprefix]'
       write(*,*) '          where optimizer = {QML_single_iteration,QML_line_search}'
       stop
    end if

    call getarg(2, optimizer)
    call getarg(3, error_type)
    call getarg(4, infile)
    call getarg(5, binfile)
    call getarg(6, clfile)
    call getarg(7, outprefix)
    ierr     = 0
    nspec    = 6
    unit     = comm_getlun()
    unit_out = unit+1
    maxiter  = 6
    nchain   = numprocs
    seed     = 37841

    ! Initialize random number generator                                                     
    if (myid == root) then
       call rand_init(handle, seed)
       do i = 1, numprocs-1
          seed = nint(rand_uni(handle)*1000000.d0)
          call mpi_send(seed, 1, MPI_INTEGER, i, 98, MPI_COMM_WORLD, ierr)
       end do
    else
       call mpi_recv(seed, 1, MPI_INTEGER, root, 98, MPI_COMM_WORLD, status, ierr)
       call rand_init(handle, seed)
    end if


    if (trim(error_type) /= 'marginal' .and. myid /= root) return

    inquire(file=trim(infile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(infile)
    inquire(file=trim(binfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(binfile)
    inquire(file=trim(clfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(clfile)

    call read_fiducial_spectrum(clfile, cls_fid)

    ! Initialize low-l module
    call comm_lowl_initialize_object(infile)
#ifdef WMAP
    use_TT       = .false.
    use_TE       = .false.
    use_lowl_TT  = .true.
    use_lowl_pol = .true.
    call wmap_likelihood_init
#endif

    ! Read binfile
    npar_powspec = 0
    open(unit,file=trim(binfile))
    do while (.true.)
       read(unit,'(A)',end=511) line
       if (line(1:1) == '#' .or. trim(line)=='') cycle
       read(line,*) lmin, lmax, status
       do i = 1, 6
          if (status(i) /= 'S') cycle
          npar_powspec = npar_powspec+1
          bins_powspec(npar_powspec,1) = lmin
          bins_powspec(npar_powspec,2) = lmax
          bins_powspec(npar_powspec,3) = i
       end do
    end do
511 close(unit)

    ! Compute C_l with non-linear search; initialize on fiducial
    allocate(cl(npar_powspec), sigma(npar_powspec,2))
    do i = 1, npar_powspec
       cl(i) = mean(cls_fid(bins_powspec(i,1):bins_powspec(i,2),bins_powspec(i,3)))
    end do
    lnL_new(1) = lnL_powspec(cl)          
    if (lnL_new(1) == -1.d30) then
       write(*,*) '   Error: Initialization spectrum not acceptable'
       stop
    end if
    
    if (trim(optimizer) == 'quasinewton') then
!       call dfpmin(cl, 1.d-6, iter, chisq, chisq_powspec, dchisq_powspec, ierr)
    else if (trim(optimizer) == 'powell_WMAP') then
#ifdef WMAP
       allocate(invF(npar_powspec,npar_powspec))

       inquire(file='cl_wmap.unf',exist=exist)
       if (exist) then
          open(18,file='cl_wmap.unf',form='unformatted')
          read(18) cl
          close(18)
       else
          call powell(cl, neg_lnL_powspec_WMAP, ierr)
       end if

       open(18,file='cl_wmap.unf',form='unformatted')
       write(18) cl
       close(18)
#endif
    else if (trim(optimizer) == 'QML_single_iteration') then

       ! Cl_QML = Cl_old + invF * dC/dCl
       allocate(dCl(npar_powspec), dCl0(npar_powspec), invF(npar_powspec,npar_powspec))
       dCl  = dlnL_powspec(cl)
       invF = get_invfisher_powspec(cl)
       cl   = cl + matmul(invF, dCl)
       do j = 1, npar_powspec
          sigma(j,:) = sqrt(invF(j,j))
       end do

    else if (trim(optimizer) == 'QML_line_search') then

       n_QML = 50
       if (trim(error_type) == 'marginal') n_QML = 1
       allocate(dCl(npar_powspec), dCl0(npar_powspec), invF(npar_powspec,npar_powspec))
       allocate(cls(npar_powspec,n_QML), sigmas(npar_powspec,2,n_QML))
       allocate(converged(npar_powspec), lnL(0:n_QML), alphas(n_QML))
       allocate(binmask(npar_powspec))
       converged = .false.
       alpha     = 1.d0
       cls       = -1.d30
       lnL(0)    = lnL_powspec(cl)

       inquire(file='cl.unf',exist=exist)
       if (exist) then
!          open(18,file='cl.unf',form='unformatted')
!          read(18) cl
!          close(18)
       end if

       do i = 1, n_QML
          write(*,*) 'Calculating NR iteration no. ', i
          call wall_time(t1)
          ! Compute proposed step size
          dCl0  = dlnL_powspec(cl)
          call wall_time(t2)
          write(*,*) '     Computing derivatives -- wall time = ', real(t2-t1,sp), ' sec'
          invF = get_invfisher_powspec(cl)

          call invert_matrix(invF)
          do j = 1, npar_powspec
             if (converged(j)) then
                invF(j,:) = 0.d0
                invF(:,j) = 0.d0
                dCl0(j)   = 0.d0
             end if
          end do
          call invert_matrix_with_mask(invF)

          do j = 1, npar_powspec
             if (converged(j)) then
                sigma(j,:) = sigmas(j,:,i-1)
             else
                sigma(j,:) = sqrt(invF(j,j))
             end if
          end do

          ! Update power spectrum, but only allow steps that yield a positive-definite
          ! S+N covariance matrix
          dCl = matmul(invF, dCl0)

          call linmin(cl, dCl, lnL_prop, neg_lnL_powspec, ierr, eps=1.d-3, xmin_out=alpha)
          lnL_new(1) = lnL_powspec(cl)          
          write(*,*) '   Step length alpha = ', alpha, lnL_new(1)

          if (abs(alpha) < 1.d-3) then
             ! Check for sharp boundaries
             allocate(p_cond(npar_powspec))
             do j = 1, npar_powspec
                ! Lower edge
                p_cond    = cl
                p_cond(j) = cl(j) - max(1.d-3*abs(cl(j)),1.d-3)
                lnL_new(1) = lnL_powspec(p_cond)          
                if (lnL_prop-lnL_new(1) > 10.d0) then
                   write(*,*) '    Fixing ', bins_powspec(j,:), ' at sharp lower edge'
                   converged(j) = .true.
                   cycle
                end if
                ! Upper edge
                p_cond    = cl
                p_cond(j) = cl(j) + max(1.d-3*abs(cl(j)),1.d-3)
                lnL_new(1) = lnL_powspec(p_cond)          
                if (lnL_prop-lnL_new(1) > 10.d0) then
                   write(*,*) '    Fixing ', bins_powspec(j,:), ' at sharp lower edge'
                   converged(j) = .true.
                   cycle
                end if
             end do
	     deallocate(p_cond)
          end if

          alphas(i)     = alpha
          cls(:,i)      = cl
          sigmas(:,:,i) = sigma
          lnL(i)        = lnL_powspec(cl)

          write(*,*) '   Log-likelihood     = ', lnL(i)
          write(*,*) '   Reduced chisq      = ', chisq_powspec, sqrt(0.5d0/comm_lowl(1)%n)
          call comm_int2string(i,itertext)
          call comm_output_powspec(trim(outprefix)//'_iter' // &
               & itertext//'_cls.dat', bins_powspec(1:npar_powspec,:), cl, sigma)

          open(18,file='cl.unf',form='unformatted')
          write(18) cl
          close(18)

          ! Check for convergence
          if (i < 4 .or. alpha < 0.3d0) cycle
          if (lnL(i)-lnL(i-3) < 0.1d0) converged = .true.

          if (all(converged)) exit
       end do

    else
       write(*,*) 'Error -- Unsupported optimizer = ', trim(optimizer)
       write(*,*) '         Supported options     = {QML_single_iteration, QML_line_search}'
       stop
    end if

!!$    open(58,file='wmap_l2_EE.dat')
!!$    n = 1000
!!$    allocate(p_cond(npar_powspec))
!!$    p_cond = cl
!!$    do i = 1, n
!!$       cl_new    = -0.04d0 + (i-1.d0)*(2.d0+0.04d0)/(n-1.d0)
!!$       p_cond(3) = cl_new
!!$       lnL_max   = lnL_powspec(p_cond)
!!$       write(58,*) cl_new, lnL_max
!!$       write(*,*) i, cl_new, lnL_max
!!$    end do
!!$    close(58)
!!$    stop



    ! Compute error bars
    invF     = get_invfisher_powspec(cl)

    if (trim(error_type) == 'marginal') then
       burnin = 1000
       nsamp  = 10000
       allocate(eta(npar_powspec), samples(npar_powspec,0:nsamp,0:nchain-1), L_prop(npar_powspec,npar_powspec))
       ! Marginal errors by MCMC
       inquire(file=trim(outprefix)//'_propmat.unf', exist=exist)
       if (exist) then
          open(unit,file=trim(outprefix)//'_propmat.unf', form='unformatted')
          read(unit) L_prop
          close(unit) 
       else	
          call cholesky_decompose(invF, L_prop)
       end if
       call mpi_bcast(L_prop, size(L_prop), MPI_DOUBLE_PRECISION, root, mpi_comm_world, ierr)
       samples           = 0
       samples(:,0,myid) = cl
       accept            = 0.d0
       blocksize         = 10000
       lnL_old = lnL_powspec(samples(:,0,myid))
       call int2string(myid+1,chain_text)
       open(58,file=trim(outprefix)//'_c'//chain_text//'.dat', recl=10000)
       ndraw = 0
       do i = 1, nsamp
	  samples(:,i,myid) = samples(:,i-1,myid)
	  do k = 1, npar_powspec/blocksize+1
             call wall_time(t1)	
             eta = 0.d0
             do j = min((k-1)*blocksize+1,npar_powspec), min(k*blocksize,npar_powspec)
                eta(j) = rand_gauss(handle)	
             end do   
             cl = samples(:,i,myid) + 0.3d0*matmul(L_prop,eta)
             lnL_prop = lnL_powspec(cl)
             call wall_time(t2)
	     ndraw = ndraw+1
             if (exp(lnL_prop-lnL_old) > rand_uni(handle)) then
                accept       = accept+1.d0
                write(*,fmt='(i5,a,i8,2f16.3,2f8.3)') myid, ' accept', i, lnL_prop, lnL_old, t2-t1, accept/ndraw
                samples(:,i,myid) = cl
                lnL_old           = lnL_prop
             else
                write(*,fmt='(i5,a,i8,2f16.3,2f8.3)') myid, ' reject', i, lnL_prop, lnL_old, t2-t1, accept/ndraw
             end if
	  end do
          write(58,*) i, real(lnL_old,sp), real(samples(:,i,myid),sp)
          if (mod(i,1000) == 0) call output_prop_mat(unit, samples(:,1:i,myid), trim(outprefix)//'_propmat.unf')
       end do
       close(58)

       allocate(tot_samples(npar_powspec,0:nsamp,0:nchain-1))
       call mpi_reduce(samples, tot_samples, size(samples), MPI_DOUBLE_PRECISION, &
           & MPI_SUM, root, MPI_COMM_WORLD, ierr)
       samples = tot_samples

       if (myid == root) then
         do j = 1, npar_powspec
            call int2string(bins_powspec(j,1), lmin_text)
            call int2string(bins_powspec(j,2), lmax_text)
            open(58,file=trim(outprefix)//'_'//stext(bins_powspec(j,3))//'_'//&
                 & lmin_text // '_' // lmax_text // '.dat')
            write(58,*) '# Marginal distribution'
	    do l = 0, nchain-1
              do k = 1, nsamp
                 write(58,*) k, real(samples(j,k,l),sp)
              end do
	      write(58,*)
	    end do
            close(58)

!            cl(j)      = mean(samples(j,1:nsamp,0))
!            sigma(j,:) = sqrt(variance(samples(j,1:nsamp)))
!            sigma(j,:) = sqrt(variance(samples(j,1:nsamp)))

	    call compute_histogram(reshape(samples(j,burnin+1:nsamp,0:nchain-1),[(nsamp-burnin)*nchain]), &
	        n_spline, cl_spline, lnL_spline)
            call compute_asymmetric_errors(cl_spline, lnL_spline, cl(j), sigma(j,1), sigma(j,2))
!            write(*,fmt='(a,i5,a,3f10.3)') 'bin = ', j, ' -- ', cl(j), sigma(j,:)

            call int2string(bins_powspec(j,1), lmin_text)
            call int2string(bins_powspec(j,2), lmax_text)
            open(58,file=trim(outprefix)//'_'//stext(bins_powspec(j,3))//'_'//&
                 & lmin_text // '_' // lmax_text // '.dat')
            write(58,*) '# Histogram'
            do k = 1, n_spline
               write(58,*) cl_spline(k), lnL_spline(k)
            end do
            write(58,*) 
            write(58,*) cl(j), minval(lnL_spline)  ! Maximum point as determined by marginal
            write(58,*) cl(j), maxval(lnL_spline)
            write(58,*) 
            write(58,*) cl(j)+sigma(j,1), minval(lnL_spline)  ! Upper 68% limit
            write(58,*) cl(j)+sigma(j,1), maxval(lnL_spline)
            write(58,*) 
            write(58,*) cl(j)-sigma(j,2), minval(lnL_spline)  ! Lower 68% limit
            write(58,*) cl(j)-sigma(j,2), maxval(lnL_spline)
            close(58)
          
            write(*,fmt='(a,i5,a,3f10.3)') 'bin = ', j, ' -- ', cl(j), sigma(j,:)
         end do

       else	
          return
       end if

    else if (trim(error_type) == 'conditional') then
       
       n        = 10000
       n_spline = 10000
       allocate(cl_cond(n), lnL_cond(n,2), p_cond(npar_powspec))
       allocate(cl_spline(n_spline), lnL_spline(n_spline))
       open(87,file=trim(outprefix)//'_maxcond_vs_cls.dat',recl=1024)
       do j = 1, npar_powspec
          !       do j = 5, 5
          n_lnL         = 0
          p_cond        = cl
          cl_new        = cl(j)
          lnL_new(1)    = lnL_powspec(p_cond)
          call comm_lowl_get_recent_value(cond_number=lnL_new(2))
          call insert_entry(cl_new, lnL_new, cl_cond, lnL_cond, n_lnL, 1.d-4)
          
          ! Check that we're not on a sharp peak
          cl_new     = cl(j) - max(1.d-3*abs(cl(j)),1.d-3)
          p_cond(j)  = cl_new
          lnL_new(1) = lnL_powspec(p_cond)          
          call comm_lowl_get_recent_value(cond_number=lnL_new(2))
          write(*,*) n_lnL, real(cl_new,sp), real(lnL_new(1),sp), real(lnL_new(1)-maxval(lnL_cond(1:n_lnL,1)),sp)
          if (maxval(lnL_cond(1:n_lnL,1))-lnL_new(1) < 10.d0) then
             call insert_entry(cl_new, lnL_new, cl_cond, lnL_cond, n_lnL, 1.d-4)
             delta       = 0.5d0 * sqrt(invF(j,j))
             cl_prev     = cl(j)
             cl_min      = -1.d30
             do while (maxval(lnL_cond(1:n_lnL,1))-lnL_cond(1,1) < 10.d0 .and. n_lnL < n .and. &
                  & delta > max(1.d-3*cl(j),1.d-4))
                cl_new     = cl_cond(1) - delta
                p_cond(j)  = cl_new
                lnL_new(1) = lnL_powspec(p_cond)
                call comm_lowl_get_recent_value(cond_number=lnL_new(2))
                write(*,*) n_lnL, real(cl_new,sp), real(lnL_new(1),sp), real(lnL_new(1)-maxval(lnL_cond(1:n_lnL,1)),sp)
                if (lnL_new(1) == -1.d30) then
                   delta  = 0.5d0 * (cl_cond(1)-cl_new) 
                   cl_min = cl_new
                else
                   delta   = 2.d0*delta
                   call insert_entry(cl_new, lnL_new, cl_cond, lnL_cond, n_lnL, 1.d-8)
                   if (cl_new-delta <= cl_min) delta = 0.5d0 * (cl_new-cl_min)
                end if
             end do
          end if
          
          cl_new     = cl(j) + max(1.d-3*abs(cl(j)),1.d-3)
          p_cond(j)  = cl_new
          lnL_new(1) = lnL_powspec(p_cond)          
          call comm_lowl_get_recent_value(cond_number=lnL_new(2))
          write(*,*) n_lnL, real(cl_new,sp), real(lnL_new(1),sp), real(lnL_new(1)-maxval(lnL_cond(1:n_lnL,1)),sp)
          if (maxval(lnL_cond(1:n_lnL,1))-lnL_new(1) < 10.d0) then
             call insert_entry(cl_new, lnL_new, cl_cond, lnL_cond, n_lnL, 1.d-8)
             delta       = 0.5d0 * sqrt(invF(j,j))
             cl_prev     = cl(j)
             cl_max      = 1.d30
             do while (maxval(lnL_cond(1:n_lnL,1))-lnL_cond(n_lnL,1) < 10.d0 .and. n_lnL < n .and. &
                  & delta > max(1.d-3*cl(j),1.d-4))
                cl_new     = cl_cond(n_lnL) + delta
                p_cond(j)  = cl_new
                lnL_new(1) = lnL_powspec(p_cond)
                call comm_lowl_get_recent_value(cond_number=lnL_new(2))
                write(*,*) n_lnL, real(cl_new,sp), real(lnL_new(1),sp), real(lnL_new(1)-maxval(lnL_cond(1:n_lnL,1)),sp)
                if (lnL_new(1) == -1.d30) then
                   delta  = 0.5d0 * (cl_new-cl_cond(n_lnL))
                   cl_max = cl_new
                else
                   delta   = 2.d0*delta
                   call insert_entry(cl_new, lnL_new, cl_cond, lnL_cond, n_lnL, 1.d-4)
                   if (cl_new+delta >= cl_max) delta = 0.5d0 * (cl_max-cl_new)
                end if
             end do
          end if
          
          ! Refine until good spline stability
          eps0 = 1.d30
          do while (.true.) 
             allocate(cl0(n_lnL-1), lnL0(n_lnL-1), lnL2(n_lnL-1), cl1(n_lnL), lnL1(n_lnL))
             eps0 = 0.d0
             m    = n_lnL
             cl1  = cl_cond(1:m)
             lnL1 = lnL_cond(1:m,1)
             do k = 2, m-1
                cl0(1:k-1)  = cl1(1:k-1)
                cl0(k:m-1)  = cl1(k+1:m)
                lnL0(1:k-1) = lnL1(1:k-1)
                lnL0(k:m-1) = lnL1(k+1:m)
                call spline(cl0, lnL0, 1.d30, 1.d30, lnL2)             
                eps = abs(lnL1(k)-splint(cl0, lnL0, lnL2, cl1(k)))
                eps0 = max(eps, eps0)
                if (eps > 0.2d0) then
                   cl_new     = 0.3333d0*cl1(k-1) + 0.6667d0*cl1(k)
                   p_cond(j)  = cl_new
                   lnL_new(1) = lnL_powspec(p_cond)
                   call comm_lowl_get_recent_value(cond_number=lnL_new(2))
                   call insert_entry(cl_new, lnL_new, cl_cond, lnL_cond, n_lnL, 1.d-8)          
                   write(*,*) 'a', n_lnL, real(cl_new,sp), real(lnL_new(1),sp), eps
                   cl_new     = 0.3333d0*cl1(k+1) + 0.6667d0*cl1(k)
                   p_cond(j)  = cl_new
                   lnL_new(1) = lnL_powspec(p_cond)
                   call comm_lowl_get_recent_value(cond_number=lnL_new(2))
                   call insert_entry(cl_new, lnL_new, cl_cond, lnL_cond, n_lnL, 1.d-8)                   
                   write(*,*) 'b', n_lnL, real(cl_new,sp), real(lnL_new(1),sp), eps
                end if
             end do
             write(*,*) 
             write(*,*) n_lnL
             write(*,*) real(cl_cond(1:n_lnL),sp)
             write(*,*) real(lnL_cond(1:n_lnL,1),sp)
             write(*,*) 
             deallocate(cl0, lnL0, lnL2, cl1, lnL1)
             if (eps0 < 0.2) exit
          end do
          
          allocate(lnL2(n_lnL))
          call spline(cl_cond(1:n_lnL), lnL_cond(1:n_lnL,1), 1.d30, 1.d30, lnL2(1:n_lnL))
          do k = 1, n_spline
             cl_spline(k)  = cl_cond(1) + (k-1.d0) * (cl_cond(n_lnL)-cl_cond(1))/(n_spline-1.d0)
             lnL_spline(k) = splint(cl_cond(1:n_lnL), lnL_cond(1:n_lnL,1), lnL2(1:n_lnL), cl_spline(k))
          end do
          deallocate(lnL2)
          lnL_spline = exp(lnL_spline-maxval(lnL_spline))
          lnL_spline = lnL_spline / tsum(cl_spline, lnL_spline)
          
          call compute_asymmetric_errors(cl_spline, lnL_spline, cl_max, sigma(j,1), sigma(j,2))

          call int2string(bins_powspec(j,1), lmin_text)
          call int2string(bins_powspec(j,2), lmax_text)
          open(58,file=trim(outprefix)//'_'//stext(bins_powspec(j,3))//'_'//&
               & lmin_text // '_' // lmax_text // '.dat')
          write(58,*) '# Unsplined log-likelihood'
          write(58,*) '# Columns are Cl, lnL, condition number'
          do k = 1, n_lnL
             write(58,*) '#  ', cl_cond(k), lnL_cond(k,1), lnL_cond(k,2)
          end do
          write(58,*) 
          write(58,*) '# Splined likelihood'
          do k = 1, n_spline
             write(58,*) cl_spline(k), lnL_spline(k)
          end do
          write(58,*) 
          write(58,*) cl_max, minval(lnL_spline)  ! Maximum point as determined by slice
          write(58,*) cl_max, maxval(lnL_spline)
          write(58,*) 
          write(58,*) cl_max+sigma(j,1), minval(lnL_spline)  ! Upper 68% limit
          write(58,*) cl_max+sigma(j,1), maxval(lnL_spline)
          write(58,*) 
          write(58,*) cl_max-sigma(j,2), minval(lnL_spline)  ! Lower 68% limit
          write(58,*) cl_max-sigma(j,2), maxval(lnL_spline)
          write(58,*) 
          write(58,*) cl(j), minval(lnL_spline)  ! Maximum point as determined by multivariate search
          write(58,*) cl(j), maxval(lnL_spline)
          close(58)
          
!          if (trim(optimizer) == 'QML_line_search') cl(j) = cl_max
          write(87,fmt='(i6,4f10.3)') j, cl(j), cl_max, sigma(j,:)
          write(*,fmt='(a,i5,a,3f10.3)') 'bin = ', j, ' -- ', cl(j), sigma(j,:)
       end do
       close(87)
       
    else if (trim(error_type) == 'fisher') then
       
       do j = 1, npar_powspec
          sigma(j,:) = sqrt(invF(j,j))
       end do
       
    else
       write(*,*) 'Unknown error type. Supported types are {marginal,conditional,fisher}'
       stop
    end if


    if (ierr /= 0) then
       write(*,*) 'Error -- powell search failed, ierr = ', ierr
       stop
    end if

    ! Output spectrum 
    call comm_output_powspec(trim(outprefix)//'_cls.dat', bins_powspec(1:npar_powspec,:), cl, sigma)

    ! Output error correlation matrix
    do i = 1, npar_powspec
       do j = i+1, npar_powspec
          invF(i,j) = invF(i,j) / sqrt(invF(i,i)*invF(j,j))
          invF(j,i) = invF(i,j)
       end do
    end do

    open(unit_out,file=trim(outprefix) // '_corrmat.dat', recl=1024)
    do i = 1, npar_powspec
       write(unit_out,fmt='(4i6)',advance='no') i, bins_powspec(i,:)
       do j = 1, i-1
          write(unit_out,fmt='(f8.3)',advance='no') invF(i,j)
       end do
       write(unit_out,fmt='(f8.3)') 1.d0
    end do
    close(unit_out)

    if (allocated(cls)) then
       open(unit_out,file=trim(outprefix) // '_convergence.dat',recl=1024)
       k = count(cls(1,:) /= -1.d30)
       do i = 1, npar_powspec
          write(unit_out,*) '#  lmin = ', bins_powspec(i,1), ', lmax = ', bins_powspec(i,2), &
               & ', spec = ', bins_powspec(i,3)
          do j = 1, k
             write(unit_out,*) j, cls(i,j), sigmas(i,:,j)
          end do
          write(unit_out,*) 
       end do
       close(unit_out)

       open(unit_out,file=trim(outprefix) // '_lnL.dat')
       k = count(cls(1,:) /= -1.d30)
       write(unit_out,*) '# Iteration      lnL      Step length'
       do i = 1, k
          write(unit_out,*) i, lnL(i), alpha
       end do
       close(unit_out)

    end if

    deallocate(cl, sigma)

  end subroutine powspec_gauss

  subroutine output_prop_mat(unit, samples, filename)
    implicit none

    integer(i4b),                   intent(in) :: unit
    real(dp),     dimension(1:,1:), intent(in) :: samples
    character(len=*),               intent(in) :: filename

    integer(i4b) :: i, j, n, m, m_tot, ierr
    real(dp), allocatable, dimension(:)   :: mu, mu_tot
    real(dp), allocatable, dimension(:,:) :: sigma, sigma_tot

    n = size(samples,1)
    m = size(samples,2)

    allocate(mu(n), mu_tot(n))
    do i = 1, n
       mu(i) = sum(samples(i,:))
    end do
    call mpi_reduce(mu, mu_tot, size(mu), MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    mu = mu_tot / (m*numprocs)
    call mpi_bcast(mu, size(mu), MPI_DOUBLE_PRECISION, root, mpi_comm_world, ierr)

    allocate(sigma(n,n), sigma_tot(n,n))
    do i = 1, n
       do j = 1, n
          sigma(i,j) = sum((samples(i,:)-mu(i))*(samples(j,:)-mu(j)))
       end do
    end do
    call mpi_reduce(sigma, sigma_tot, size(sigma), MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    sigma = sigma_tot / (m*numprocs-1)

    if (myid == 0) then
       call cholesky_decompose(sigma, sigma_tot, ierr)
       if (ierr == 0) then
          open(unit,file=trim(filename),form='unformatted')
          write(unit) sigma_tot	
          close(unit) 
       else
          open(unit,file=trim(filename))
          call get_eigenvalues(sigma, mu)
          do i = 1, n
	     write(unit,*) i, mu(i)
  	  end do
	  close(unit)
       end if
    end if

    deallocate(mu, mu_tot, sigma, sigma_tot)

  end subroutine


  subroutine print_fisher_mat
    implicit none

    character(len=256) :: infile, binfile, outfile, clfile, outprefix, line
    character(len=2)   :: spectext, ptext
    character(len=4)   :: lmin_text, lmax_text
    integer(i4b)       :: nspec, lmin, lmax, numchain, numiter, g, i, j, k, l, m, n, p, n_h
    integer(i4b)       :: b_min, b_max, id, numbin, ind, npar
    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep, nmaps, ierr
    integer(i4b)       :: unit, unit_out
    logical(lgt)       :: exist, inv_error
    real(dp),     allocatable, dimension(:)      :: d
    real(dp),     allocatable, dimension(:,:)    :: invF, invF2, PY, cls_fid
    real(dp),     allocatable, dimension(:,:)    :: beam, N_cov
    real(dp),     allocatable, dimension(:,:,:)  :: w
    integer(i4b), allocatable, dimension(:,:)    :: bins
    character(len=1), dimension(6) :: status

    if (iargc() /= 5) then
       write(*,*) '    Print Fisher matrix and error bars'
       write(*,*) '      Usage: comm_like_tools print_fisher_mat [datafile] [binfile] [fid cl file] [outprefix]'
       stop
    end if

    call getarg(2, infile)
    call getarg(3, binfile)
    call getarg(4, clfile)
    call getarg(5, outprefix)
    ierr     = 0
    nspec    = 6
    unit     = comm_getlun()
    unit_out = unit+1

    inquire(file=trim(infile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(infile)
    inquire(file=trim(binfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(binfile)
    inquire(file=trim(clfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(clfile)

    call read_fiducial_spectrum(clfile, cls_fid)

    ! Read data file
    write(*,*) 'Reading input data'
    open(unit,file=trim(infile),form='unformatted')
    read(unit) lmax, nmaps
    read(unit) n_h, n
    allocate(w(0:lmax,nmaps,nmaps))
    allocate(beam(0:lmax,nmaps))
    allocate(d(n), N_cov(n,n))
    allocate(PY(n_h,n))
    read(unit) d
    read(unit) N_cov
    read(unit) w
    read(unit) beam
    read(unit) PY
    close(unit)

    ! Read binfile
    allocate(bins(1000,3))
    npar = 0
    open(unit,file=trim(binfile))
    do while (.true.)
       read(unit,'(A)',end=513) line
       if (line(1:1) == '#' .or. trim(line)=='') cycle
       read(line,*) lmin, lmax, status
       do i = 1, 6
          if (status(i) /= 'S') cycle
          npar = npar+1
          bins(npar,1) = lmin
          bins(npar,2) = lmax
          bins(npar,3) = i
       end do
    end do
513 close(unit)

    ! Compute Fisher matrix
    allocate(invF(npar,npar), invF2(npar,npar))
    write(*,*) 'Computing Fisher matrix'
    call get_invfisher(cls_fid, w, transpose(PY), N_cov, bins(1:npar,:), invF, inv_error)
    if (inv_error) stop
    !call get_invfisher_safe(cls_fid, w, transpose(PY), N_cov, bins(1:npar,:), invF2)

!!$    open(58,file='fisher_comparsion.dat', recl=1024)
!!$    do i = 1, npar
!!$       do j = 1, npar
!!$          write(58,*) i, j, invF(i,j), invF2(i,j), abs(invF(i,j)-invF2(i,j))/abs(invF(i,j)+invF2(i,j))
!!$       end do
!!$    end do
!!$    close(58)
!!$    stop

    ! Output error bars
    write(*,*) 'Outputting Fisher matrix and uncertainties'
    open(unit,file=trim(outprefix)//'_sigma.dat')
    do i = 1, nspec
       n = 0
       do j = 1, npar
          if (bins(j,3) == i) then
             write(unit,*) 0.5*sum(bins(j,1:2)), sqrt(invF(j,j))
             n = n+1
          end if
       end do
       if (n > 0) write(unit,*)
    end do

    ! Output error correlation matrix
    do i = 1, npar
       do j = i+1, npar
          invF(i,j) = invF(i,j) / sqrt(invF(i,i)*invF(j,j))
          invF(j,i) = invF(i,j)
       end do
    end do

    open(unit_out,file=trim(outprefix) // '_corrmat.dat', recl=1024)
    do i = 1, npar
       write(unit_out,fmt='(4i6)',advance='no') i, bins_powspec(i,:)
       do j = 1, i-1
          write(unit_out,fmt='(f8.3)',advance='no') invF(i,j)
       end do
       write(unit_out,fmt='(f8.3)') 1.d0
    end do
    close(unit_out)

    deallocate(cls_fid)

  end subroutine print_fisher_mat

  function neg_lnL_powspec(p)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp)                                     :: neg_lnL_powspec

    integer(i4b) :: i, lmax, ierr
    real(dp)     :: lnL, chisq
    real(dp), allocatable, dimension(:,:) :: cls

    ierr = 0
    lmax = comm_lowl(1)%lmax
    allocate(cls(0:lmax,6))
    cls = cls_fid(0:lmax,:)
    do i = 1, npar_powspec
       cls(bins_powspec(i,1):bins_powspec(i,2),bins_powspec(i,3)) = p(i)
    end do

    neg_lnL_powspec = -comm_lowl_compute_lnL(cls=cls, ierr=ierr, enforce_pos_def=.false., &
         & red_chisq=chisq_powspec)
    if (ierr /= 0) neg_lnL_powspec = 1.d30

!    write(*,fmt='(e16.8,4e10.3)') neg_lnL_powspec, cls(2,[1,2,4,6])

    deallocate(cls)

  end function neg_lnL_powspec

#ifdef WMAP
  function neg_lnL_powspec_WMAP(p)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp)                                     :: neg_lnL_powspec_WMAP

    integer(i4b) :: i, lmax, ierr
    real(dp)     :: lnL, chisq
    real(dp), allocatable, dimension(:,:) :: cls, cls_WMAP
    real(dp)     :: like(num_WMAP)

    ierr = 0
    lmax = comm_lowl(1)%lmax
    allocate(cls(0:lmax,6))
    cls = cls_fid(0:lmax,:)
    do i = 1, npar_powspec
       cls(bins_powspec(i,1):bins_powspec(i,2),bins_powspec(i,3)) = p(i)
    end do

    like = 0.d0
    allocate(cls_WMAP(ttmin:ttmax,6))
    cls_WMAP           = 0.d0
    cls_WMAP(2:lmax,:) = cls(2:lmax,:)
    ! WMAP likelihood returns -lnL
    call wmap_likelihood_compute(cls_WMAP(:,1), cls_WMAP(:,2), cls_WMAP(:,4), cls_WMAP(:,6), like)
    neg_lnL_powspec_WMAP = sum(like)
!    write(*,*) 'like = ', real(like,sp)
    if (like(ttlowllike) > 60.d0 .or. like(lowllike) == 0.d0) ierr = 1
    deallocate(cls_WMAP)
    if (ierr /= 0) neg_lnL_powspec_WMAP = 1.d30

    write(*,fmt='(e16.8,4e10.3)') neg_lnL_powspec_WMAP, cls(2,[1,2,4,6])

    deallocate(cls)

  end function neg_lnL_powspec_WMAP
#endif


  function lnL_powspec(p)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp)                                     :: lnL_powspec

    lnL_powspec = -neg_lnL_powspec(p)

  end function lnL_powspec

!!$  function lnL_powspec(p)
!!$    use healpix_types
!!$    implicit none
!!$    real(dp), dimension(:), intent(in), optional :: p
!!$    real(dp)                                     :: lnL_powspec
!!$
!!$    integer(i4b) :: i, lmax, ierr
!!$    real(dp)     :: lnL, chisq
!!$    real(dp), allocatable, dimension(:,:) :: cls
!!$
!!$    lmax = comm_lowl(1)%lmax
!!$    allocate(cls(0:lmax,6))
!!$    cls = cls_fid(0:lmax,:)
!!$    do i = 1, npar_powspec
!!$       cls(bins_powspec(i,1):bins_powspec(i,2),bins_powspec(i,3)) = p(i)
!!$    end do
!!$    lnL_powspec = comm_lowl_compute_lnL(cls=cls, ierr=ierr, enforce_pos_def=.false., &
!!$         & red_chisq=chisq_powspec)
!!$    if (ierr /= 0) lnL_powspec = -1.d30
!!$    deallocate(cls)
!!$
!!$  end function lnL_powspec

  function dlnL_powspec(p)
    USE healpix_types
    IMPLICIT NONE
    REAL(dp), DIMENSION(:), INTENT(IN) :: p
    REAL(dp), DIMENSION(size(p)) :: dlnL_powspec

    integer(i4b) :: i, lmax, ierr, n, n_h, l, m, ind, col, row, numcomp, s
    real(dp)     :: t1, t2
    real(dp), allocatable, dimension(:,:) :: cls, C, invC, d, P_invC, P_invC_D, A, invC_Pt, D_invC_Pt
    integer(i4b), dimension(6) :: ind1 = [1, 1, 1, 2, 2, 3]
    integer(i4b), dimension(6) :: ind2 = [1, 2, 3, 2, 3, 3]
    real(dp),     dimension(6) :: mult = [1.d0, 2.d0, 2.d0, 1.d0, 2.d0, 1.d0]

    ! Compute covariance matrix
    lmax = comm_lowl(1)%lmax
    allocate(cls(0:lmax,6))
    cls = cls_fid(0:lmax,:)
    do i = 1, npar_powspec
       cls(bins_powspec(i,1):bins_powspec(i,2),bins_powspec(i,3)) = p(i)
    end do

    ! Precompute A = P*invC*D*invC*P^t
    n       = comm_lowl(1)%n
    n_h     = comm_lowl(1)%n_h
    numcomp = n_h / comm_lowl(1)%nmaps
    allocate(A(n_h,n_h), invC(n,n), d(n,1), D_invC_Pt(n,n_h), invC_Pt(n,n_h))
    call comm_lowl_getC(.false., C, ierr, cls=cls)
    invC    = C
    d       = comm_lowl(1)%d(:,1:1)
    C = matmul(d, transpose(d)) - C
    call dpotrf('L', n, invC, n, ierr)
    if (ierr /= 0) then
       write(*,*) 'dchisq -- Cholesky decomposition failed, ierr = ', ierr
       stop
    end if
    invC_Pt = transpose(comm_lowl(1)%P_harm)
    call dpotrs('L', n, n_h, invC, n, invC_Pt, n, ierr)    
    call dgemm('N', 'N', n, n_h, n, 1.d0, C, n, invC_Pt, n, 0.d0, D_invC_Pt, n)
    call dgemm('T', 'N', n_h, n_h, n, 1.d0, invC_Pt, n, D_invC_Pt, n, 0.d0, A, n_h)

    ! Compute derivatives
    dlnL_powspec = 0.d0
    do i = 1, npar_powspec
       s   = bins_powspec(i,3)
       row = (ind1(s)-1)*numcomp + bins_powspec(i,1)**2 + 1
       col = (ind2(s)-1)*numcomp + bins_powspec(i,1)**2 + 1
       do l = bins_powspec(i,1), bins_powspec(i,2)
          do m = -l, l
             dlnL_powspec(i) = dlnL_powspec(i) + 0.5d0 * mult(s) * A(row,col) * &
                  & comm_lowl(1)%beam(l,ind1(s)) * comm_lowl(1)%beam(l,ind2(s)) * &
                  & 2.d0*pi/real(l*(l+1),dp)
             row = row + 1
             col = col + 1
          end do
       end do
    end do

    deallocate(cls, A, invC, d, D_invC_Pt, invC_Pt)

  end function dlnL_powspec

  function get_sigma_powspec(p)
    use healpix_types
    IMPLICIT NONE
    REAL(dp), DIMENSION(:), INTENT(IN) :: p
    REAL(dp), DIMENSION(size(p)) :: get_sigma_powspec

    integer(i4b) :: i
    real(dp), allocatable, dimension(:,:) :: invF

    ! Get inverse Fisher matrix
    allocate(invF(npar_powspec,npar_powspec))
    invF = get_invfisher_powspec(p)

    ! Return diagonal elements
    do i = 1, npar_powspec
       get_sigma_powspec(i) = sqrt(invF(i,i))
    end do

    deallocate(invF)

  end function get_sigma_powspec


  function get_invfisher_powspec(p)
    use healpix_types
    implicit none

    real(dp), DIMENSION(:), INTENT(IN) :: p
    real(dp), DIMENSION(size(p),size(p)) :: get_invfisher_powspec

    integer(i4b) :: i, j, lmax, ierr
    logical(lgt) :: inv_error
    real(dp), allocatable, dimension(:,:) :: cls, invF

    ! Compute inverse covariance matrix
    lmax    = comm_lowl(1)%lmax
    allocate(cls(0:lmax,6), invF(npar_powspec,npar_powspec))
    cls = cls_fid(0:lmax,:)
    do i = 1, npar_powspec
       cls(bins_powspec(i,1):bins_powspec(i,2),bins_powspec(i,3)) = p(i)
    end do

    call get_invfisher(cls, comm_lowl(1)%w, transpose(comm_lowl(1)%P_harm), comm_lowl(1)%N_cov, &
         & bins_powspec(1:npar_powspec,:), invF, inv_error)    

    get_invfisher_powspec = invF

    deallocate(cls, invF)

!!$    ! Precompute P * invC * P^t
!!$    allocate(P_invC(n_h,n), P_invC_Pt(n_h,n_h), A(n_h,n_h), P_invC_Pt_dC(n_h,n_h))
!!$    allocate(indmap(n_h), invC_Pt(n,n_h))
!!$
!!$    call wall_time(t1)
!!$    call comm_lowl_getC(.false., invC, ierr, cls=cls)
!!$    call dpotrf('L', n, invC, n, ierr)
!!$    invC_Pt = transpose(comm_lowl(1)%P_harm)
!!$    call dpotrs('L', n, n_h, invC, n, invC_Pt, n, ierr)    
!!$    call dgemm('N', 'N', n_h, n_h, n, 1.d0, comm_lowl(1)%P_harm, n_h, &
!!$         & invC_Pt, n, 0.d0, P_invC_Pt, n_h)
!!$    call wall_time(t2)
!!$    write(*,*) '     Computing P*invC*P^t -- wall time = ', real(t2-t1,sp), ' sec'
!!$
!!$    ! Compute Fisher matrix
!!$    allocate(F(npar_powspec,npar_powspec))
!!$    F = 0.d0
!!$    call wall_time(t1)
!!$    do i = 1, npar_powspec
!!$
!!$       ! Precompute P*invC*Pt*dC
!!$       s   = bins_powspec(i,3)
!!$       row  = (ind1(s)-1)*numcomp + bins_powspec(i,1)**2 + 1
!!$       col  = (ind2(s)-1)*numcomp + bins_powspec(i,1)**2 + 1
!!$       P_invC_Pt_dC = 0.d0
!!$       do l = bins_powspec(i,1), bins_powspec(i,2)
!!$          do m = -l, l
!!$             do j = 1, n_h
!!$                P_invC_Pt_dC(j,col) = P_invC_Pt(j,row) * 2.d0*pi/real(l*(l+1),dp) * &
!!$                     & comm_lowl(1)%beam(l,ind1(s)) * comm_lowl(1)%beam(l,ind2(s))
!!$                if (row /= col) then
!!$                   P_invC_Pt_dC(j,row) = P_invC_Pt(j,col) * 2.d0*pi/real(l*(l+1),dp) * &
!!$                        & comm_lowl(1)%beam(l,ind1(s)) * comm_lowl(1)%beam(l,ind2(s))
!!$                end if
!!$             end do
!!$             row = row + 1
!!$             col = col + 1
!!$          end do
!!$       end do
!!$
!!$       nmode = 0
!!$       do j = 1, n_h
!!$          if (any(P_invC_Pt_dC(:,j) /= 0.d0)) then
!!$             nmode         = nmode+1
!!$             indmap(nmode) = j
!!$          end if
!!$       end do
!!$
!!$       call dgemm('N', 'N', n_h, n_h, nmode, 1.d0, P_invC_Pt_dC(:,indmap(1:nmode)), n_h, &
!!$            & P_invC_Pt(indmap(1:nmode),:), nmode, 0.d0, A, n_h)
!!$
!!$       do j = i, npar_powspec
!!$          s   = bins_powspec(j,3)
!!$          row = (ind1(s)-1)*numcomp + bins_powspec(j,1)**2 + 1
!!$          col = (ind2(s)-1)*numcomp + bins_powspec(j,1)**2 + 1
!!$          do l = bins_powspec(j,1), bins_powspec(j,2)
!!$             do m = -l, l
!!$                F(i,j) = F(i,j) + 0.5d0 * mult(s) * A(row,col) * &
!!$                     & comm_lowl(1)%beam(l,ind1(s)) * comm_lowl(1)%beam(l,ind2(s)) * &
!!$                     & 2.d0*pi/real(l*(l+1),dp)
!!$                row = row + 1
!!$                col = col + 1
!!$             end do
!!$          end do
!!$          F(j,i) = F(i,j)
!!$       end do
!!$
!!$    end do
!!$    call wall_time(t2)
!!$    write(*,fmt='(a,f8.3,a)') &
!!$         & '      Computing Fisher matrix -- wall time = ', real(t2-t1,sp), ' sec'
!!$
!!$    allocate(W(npar_powspec))
!!$    call get_eigenvalues(F, W)
!!$    if (any(W <= 0.d0) .or. W(npar_powspec)/W(1) > 1.d12) then
!!$       write(*,*) 'Error -- Fisher matrix not invertible'
!!$       write(*,*) '         Condition number = ', W(npar_powspec)/W(1)
!!$       write(*,*) '         Either remove spectrum bins, or add more basis vectors'
!!$       stop
!!$    end if
!!$    call invert_matrix(F)
!!$
!!$    get_invfisher_powspec = F
!!$
!!$    deallocate(cls, invC, P_invC, P_invC_Pt, A, F, indmap, W)

  end function get_invfisher_powspec


!!$  subroutine sigma2par_BR(partype)
!!$    implicit none
!!$
!!$    character(len=*), intent(in)   :: partype
!!$
!!$    character(len=256) :: infile, binfile, outfile, temp, line, clfile, outprefix, templatefile
!!$    character(len=2)   :: spectext
!!$    character(len=4)   :: lmin_text, lmax_text
!!$    real(dp)           :: cl_min, cl_max, sigma, top, left, right, int, tot_int, cl_in(4)
!!$    real(dp)           :: A_min, A_max, A0, dA, mu, norm
!!$    integer(i4b)       :: nspec, lmin, lmax, numchain, numiter, i, j, k, l, m, n, p, bins(4), id, numbin, ind
!!$    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep, nmaps, ierr, spec_fit
!!$    logical(lgt)       :: exist, include_BB
!!$    real(dp), allocatable, dimension(:,:) :: cls0, cls
!!$    real(dp), allocatable, dimension(:)   :: amp, lnL, F
!!$
!!$    if (iargc() /= 14) then
!!$       write(*,*) '    Compute spectrum with Blackwell-Rao'
!!$       write(*,*) '    Options: [sigmafile] [clfile] [binfile] [firstchain] [lastchain] '
!!$       write(*,*) '                 [firstsample] [lastsample] [thinstep] [templatefile] '
!!$       write(*,*) '                 [default_amplitude] [spec_fit] [include_BB] [outprefix]'
!!$       stop
!!$    end if
!!$
!!$    call getarg(2, infile)
!!$    call getarg(3, clfile)
!!$    call getarg(4, binfile)
!!$    call getarg(5, temp)
!!$    read(temp,*) firstchain
!!$    call getarg(6, temp)
!!$    read(temp,*) lastchain
!!$    call getarg(7, temp)
!!$    read(temp,*) firstsample
!!$    call getarg(8, temp)
!!$    read(temp,*) lastsample
!!$    call getarg(9, temp)
!!$    read(temp,*) thinstep
!!$    call getarg(10, templatefile)
!!$    call getarg(11, temp)
!!$    read(temp,*) A0
!!$    call getarg(12, temp)
!!$    read(temp,*) spec_fit
!!$    call getarg(13, temp)
!!$    read(temp,*) include_BB
!!$    call getarg(14, outprefix)
!!$
!!$    inquire(file=trim(binfile), exist=exist)
!!$    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(binfile)
!!$    inquire(file=trim(infile), exist=exist)
!!$    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(infile)
!!$    inquire(file=trim(clfile), exist=exist)
!!$    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(clfile)
!!$    inquire(file=trim(clfile), exist=exist)
!!$    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(clfile)
!!$
!!$    ! Read fiducial spectrum file
!!$    unit = comm_br_getlun()
!!$    open(unit,file=trim(templatefile))
!!$    lmax = -1
!!$    do while (.true.)
!!$       read(unit,'(a)',end=79) line
!!$       line = trim(line)
!!$       if (line(1:1) == '#') cycle
!!$       read(line,*) l
!!$       lmax = max(lmax, l)
!!$    end do
!!$79  close(unit)
!!$    if (lmax > -1) then
!!$       allocate(cls0(0:lmax,6), cls(0:lmax,6))
!!$       cls0 = 0.d0
!!$       open(unit,file=trim(templatefile))
!!$       do while (.true.)
!!$          read(unit,'(a)',end=80) line
!!$          line = trim(line)
!!$          if (line(1:1) == '#') cycle
!!$          read(line,*) l, cl_in
!!$          cls0(l,1) = cl_in(1) ! Assume (TT, EE, BB, TE) ordering for now
!!$          cls0(l,2) = cl_in(4)          
!!$          cls0(l,4) = cl_in(2)
!!$          cls0(l,6) = cl_in(3)
!!$       end do
!!$80     close(unit)
!!$    else
!!$       write(*,*) 'Error -- comm_br_mod: No valid entries in clfile = ', trim(clfile)
!!$       stop
!!$    end if
!!$    if (.not. include_BB) cls0(:,6) = 0.d0
!!$
!!$    ! Initialize BR module
!!$    nmaps = 3
!!$    call comm_br_initialize_object(infile, clfile, firstchain, lastchain, &
!!$       & firstsample, lastsample, thinstep, binfile=binfile, handle=id, prefix=outprefix)
!!$
!!$    ! Loop over all free parameters, and output slices
!!$    A_min =  0.d0
!!$    A_max =  10.d0
!!$    numbin = 10000
!!$    allocate(amp(numbin), lnL(numbin), F(numbin))
!!$    cls = cls0
!!$    do l = 1, numbin
!!$       amp(l) = A_min + (A_max-A_min)/(numbin-1.d0) * (l-1)
!!$       cls(:,spec_fit) = amp(l) * cls0(:,spec_fit)
!!$       lnL(l) = comm_br_compute_lnL(cls, ierr, id)
!!$       if (trim(partype) == 'tau') then
!!$          amp(l) = sqrt(amp(l)) * A0
!!$       else
!!$          amp(l) = amp(l) * A0
!!$       end if
!!$    end do
!!$    lnL = exp(lnL-maxval(lnL))
!!$    norm = 0.d0
!!$    do i = 1, numbin-1
!!$       norm = norm + 0.5d0*(lnL(i)+lnL(i+1)) * (amp(i+1)-amp(i))
!!$    end do
!!$    lnL = lnL / norm
!!$
!!$    outfile = trim(outprefix) // '.dat'
!!$    write(*,*) 'Writing ', trim(outfile)
!!$    open(unit,file=trim(outfile))
!!$    do l = 1, numbin
!!$       if (lnL(l) > 1.d-6 * maxval(lnL)) write(unit,*) amp(l), lnL(l)
!!$    end do
!!$    close(unit)
!!$
!!$    mu    = 0.d0
!!$    do i = 1, numbin-1
!!$       mu = mu + 0.5d0 * (amp(i)*lnL(i)+amp(i+1)*lnL(i+1)) * (amp(i+1)-amp(i))
!!$    end do
!!$    sigma = 0.d0
!!$    do i = 1, numbin-1
!!$       sigma = sigma + 0.5d0*((amp(i)-mu)**2*lnL(i) + (amp(i+1)-mu)**2*lnL(i+1)) * (amp(i+1)-amp(i))
!!$    end do
!!$    sigma = sqrt(sigma)
!!$    write(*,fmt='(a,f8.4,a,f8.4)') ' Estimated Gaussian parameter = ', real(mu,sp), ' +/-', real(sigma,sp)
!!$
!!$    ! Compute cumulative distribution
!!$    F(1) = 0.d0
!!$    do i = 2, numbin
!!$       F(i) = F(i-1) + 0.5d0 * (lnL(i)+lnL(i-1)) * (amp(i)-amp(i-1))
!!$    end do
!!$    F = F / F(numbin)
!!$    
!!$    do i = 1, numbin
!!$       if (F(i) > 0.95d0) exit
!!$    end do
!!$
!!$    write(*,fmt='(a,f8.4)') ' Upper 95% confidence limit   = ', real(amp(i),sp)
!!$    
!!$  end subroutine sigma2par_BR

  subroutine sigma2qn_BR
    implicit none

    character(len=256) :: BRfile, outfile, temp, line, clfile, outprefix, templatefile
    character(len=2)   :: spectext
    character(len=4)   :: lmin_text, lmax_text
    real(dp)           :: cl_min, cl_max, sigma, top, left, right, int, tot_int, cl_in(4)
    real(dp)           :: q_min, q_max, n_min, n_max, mu, l_pivot
    integer(i4b)       :: nspec, lmin, lmax, numchain, numiter, i, j, k, l, m, p, bins(4), id, numbin, ind
    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep, nmaps, ierr, spec_fit
    logical(lgt)       :: exist
    real(dp), allocatable, dimension(:,:) :: cls0, cls, lnL
    real(dp), allocatable, dimension(:)   :: q, n, lnL_marg

    if (iargc() /= 11) then
       write(*,*) '    Compute 2D q-n contours with Blackwell-Rao'
       write(*,*) '    Options: [BR file] [templatefile] [q_min] [q_max] [n_min] [n_max] [numbin]'
       write(*,*) '            [l_pivot] [spec_fit] [outprefix]'
       stop
    end if

    call getarg(2, BRfile)
    call getarg(3, templatefile)
    call getarg(4, temp)
    read(temp,*) q_min
    call getarg(5, temp)
    read(temp,*) q_max
    call getarg(6, temp)
    read(temp,*) n_min
    call getarg(7, temp)
    read(temp,*) n_max
    call getarg(8, temp)
    read(temp,*) numbin
    call getarg(9, temp)
    read(temp,*) l_pivot
    call getarg(10, temp)
    read(temp,*) spec_fit
    call getarg(11, outprefix)

    inquire(file=trim(BRfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(BRfile)
    inquire(file=trim(templatefile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(templatefile)

    call read_fiducial_spectrum(templatefile, cls0)
    lmax = size(cls0,1)-1
    allocate(cls(0:lmax,6))

    ! Initialize BR module
    call comm_br_initialize_object(BRfile, handle=id)

    ! Loop over all free parameters, and output slices
    allocate(q(numbin), n(numbin), lnL(numbin,numbin))
    cls = cls0
    do i = 1, numbin
       q(i) = q_min + (q_max-q_min)/(numbin-1.d0) * (i-1)
       do j = 1, numbin
          n(j) = n_min + (n_max-n_min)/(numbin-1.d0) * (j-1)
          do l = 2, lmax
             cls(l,spec_fit) = q(i) * (l/l_pivot)**n(j) * cls0(l,spec_fit)
          end do
          lnL(i,j) = comm_br_compute_lnL(cls, ierr, handle=id)
       end do
    end do
    lnL = exp(lnL-maxval(lnL))
    lnL = lnL / sum(lnL) / ((q_max-q_min)/(numbin-1)) / ((n_max-n_min)/(numbin-1))

    outfile = trim(outprefix) // '_qn.dat'
    write(*,*) 'Writing ', trim(outfile)
    open(unit,file=trim(outfile))
    write(unit,*) numbin, numbin
    do i = 1, numbin
       do j = 1, numbin
          write(unit,*) q(i), n(j), lnL(i,j)
       end do
    end do
    close(unit)

    ! Compute marginals
    allocate(lnL_marg(numbin))
    do i = 1, numbin
       lnL_marg(i) = sum(lnL(i,:))
    end do
    lnL_marg = lnL_marg / sum(lnL_marg) / ((q_max-q_min)/(numbin-1))
    mu    = sum(q * lnL_marg) * (q_max-q_min)/(numbin-1)
    sigma = sqrt(sum((q-mu)**2 * lnL_marg) * (q_max-q_min)/(numbin-1))
    write(*,fmt='(a,f8.4,a,f8.4)') 'Estimated q = ', real(mu,sp), ' +/-', real(sigma,sp)

    outfile = trim(outprefix) // '_q.dat'
    open(unit,file=trim(outfile))
    do i = 1, numbin
       write(unit,*) q(i), lnL_marg(i)
    end do
    close(unit)

    do i = 1, numbin
       lnL_marg(i) = sum(lnL(:,i))
    end do
    lnL_marg = lnL_marg / sum(lnL_marg) / ((n_max-n_min)/(numbin-1))
    mu    = sum(n * lnL_marg) * (n_max-n_min)/(numbin-1)
    sigma = sqrt(sum((n-mu)**2 * lnL_marg) * (n_max-n_min)/(numbin-1))
    write(*,fmt='(a,f8.4,a,f8.4)') 'Estimated n = ', real(mu,sp), ' +/-', real(sigma,sp)
    
    outfile = trim(outprefix) // '_n.dat'
    open(unit,file=trim(outfile))
    do i = 1, numbin
       write(unit,*) n(i), lnL_marg(i)
    end do
    close(unit)

  end subroutine sigma2qn_BR

  subroutine compute_conf_limits(x, P, xmax, error_low, error_high)
    implicit none
    
    real(dp), dimension(:), intent(in)  :: x, P
    real(dp),               intent(out) :: xmax, error_low, error_high

    integer(i4b) :: b_max(1), b_low, b_high, n
    real(dp)     :: integral, tot_integral

    b_max   = maxloc(P)
    b_low   = b_max(1)
    b_high  = b_max(1)
    n       = size(P)

    tot_integral = sum(P)
    integral     = P(b_max(1))
    do while (integral < 0.68*tot_integral)
       if (b_low == 1) then
          b_high   = b_high+1
          integral = integral + P(b_high)
       else if (b_high == n) then
          b_low    = b_low-1
          integral = integral + P(b_low)
       else if (P(b_low-1) > P(b_high+1)) then
          b_low    = b_low-1
          integral = integral + P(b_low)
       else
          b_high   = b_high+1
          integral = integral + P(b_high)
       end if
    end do

    xmax       = x(b_max(1))
    error_low  = xmax-x(b_low)
    error_high = x(b_high)-xmax

  end subroutine compute_conf_limits

  subroutine print_coupling_kernel(lmax, nmaps, lhigh, P, prefix)
    implicit none
    
    integer(i4b),                     intent(in) :: lmax, nmaps, lhigh
    real(dp),         dimension(:,:), intent(in) :: P
    character(len=*),                 intent(in) :: prefix
    
    integer(i4b) :: l, k, j, ind, i1, i2, l1, l2, m1, m2, ind1, ind2
    real(dp), allocatable, dimension(:,:,:,:) :: K_llp
    
    ! Compute (l,l') coupling matrix
    allocate(K_llp(0:lmax,nmaps,0:lmax,nmaps))
    K_llp = 0.d0
    ind1  = 0
    do i1 = 1, nmaps
       do l1 = 0, lmax
          do m1 = -l1, l1
             ind1 = ind1+1
             
             ind2 = 0
             do i2 = 1, nmaps
                do l2 = 0, lmax
                   do m2 = -l2, l2
                      ind2 = ind2+1
                      K_llp(l1,i1,l2,i2) = K_llp(l1,i1,l2,i2) + P(ind1,ind2)
                   end do
                end do
             end do
             
          end do
       end do
    end do
    do i1 = 1, nmaps
       do l1 = 0, lmax
          do i2 = 1, nmaps
             do l2 = 0, lmax
                K_llp(l1,i1,l2,i2) = K_llp(l1,i1,l2,i2) / real((2*l1+1)*(2*l2+1),dp)
             end do
          end do
       end do
    end do
    do i1 = 1, nmaps
       do l1 = 0, lmax
          do i2 = 1, nmaps
             do l2 = 0, lmax
                if (i1 /= i2 .or. l1 /= l2) then
                   if (K_llp(l1,i1,l1,i1)*K_llp(l2,i2,l2,i2) > 0.d0) then
                      K_llp(l1,i1,l2,i2) = K_llp(l1,i1,l2,i2) / &
                           & sqrt(K_llp(l1,i1,l1,i1)*K_llp(l2,i2,l2,i2))
                   end if
                end if
             end do
          end do
       end do
    end do
    do i1 = 1, nmaps
       do l1 = 0, lmax
          if (K_llp(l1,i1,l1,i1) /= 0.d0) K_llp(l1,i1,l1,i1) = 1.d0
       end do
    end do
    
    open(unit,file=trim(prefix) // '_lhigh_kernel.dat')
    do j = 1, nmaps
       ind = 0
       do k = 1, nmaps
          do l = 0, lmax
             write(unit,*) ind, abs(K_llp(l,k,lhigh,j))
             ind = ind+1
          end do
       end do
       write(unit,*)
    end do
    close(unit)
    
    deallocate(K_llp)
    
  end subroutine print_coupling_kernel

end program comm_like_tools
