program comm_like_tools
  use comm_proc_utils
  use comm_br_mod
  use comm_br_old_mod
  use comm_gauss_br_mod
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

  character(len=128) :: operation
  integer(i4b)       :: unit

  ! Global parameters for power spectrum search
  real(dp)                                     :: chisq_powspec
  integer(i4b)                                 :: npar_powspec
  integer(i4b),              dimension(1000,3) :: bins_powspec
  real(dp),     allocatable, dimension(:,:)    :: cls_fid

  if (iargc() == 0) then
     write(*,*) 'Usage: comm_like_tools [operation] [args] ... '
     write(*,*) '  Valid operations are:'
     write(*,*) '     mapcov2gausslike     -- '
     write(*,*) '     coadd_dataset        -- '
     write(*,*) '     sigma2cl_BR          -- '
     write(*,*) '     sigma2cl_gauss_BR    -- '
     write(*,*) '     sigma2slice_BR       -- '
     write(*,*) '     sigma2par_BR         -- '
     write(*,*) '     sigma2tau_BR         -- '
     write(*,*) '     sigma2qn_BR          -- '
     write(*,*) '     sigma2gauss          -- '
     write(*,*) '     slice_gauss          -- '
     write(*,*) '     powspec_gauss        -- '
     write(*,*) '     optimization_search  -- '
     write(*,*) '     print_fisher_mat     -- '
     write(*,*) '     compute_delta_lnL    -- '
     write(*,*) '     print_lnL            -- '
     write(*,*) '     print_lnL_BR         -- '
     write(*,*) '     print_lnL_old_BR     -- '
     write(*,*) '     timing               -- '
     write(*,*) '     fast_par_estimation  -- '
     write(*,*) '     output_covmat_slice  -- '
     write(*,*) '  For further usage information, type "comm_like_tools [operation]"'
     stop
  end if

  unit = 19
  call getarg(1,operation)
  if (trim(operation) == 'mapcov2gausslike') then
     call mapcov2gausslike
  else if (trim(operation) == 'coadd_dataset') then
     call coadd_dataset
  else if (trim(operation) == 'sigma2cl_BR') then
     call sigma2cl_BR
  else if (trim(operation) == 'sigma2cl_gauss_BR') then
     call sigma2cl_gauss_BR
  else if (trim(operation) == 'sigma2slice_BR') then
     call sigma2slice_BR
  else if (trim(operation) == 'slice_gauss') then
     call slice_gauss
  else if (trim(operation) == 'sigma2par_BR') then
     call sigma2par_BR('par')
  else if (trim(operation) == 'sigma2tau_BR') then
     call sigma2par_BR('tau')
  else if (trim(operation) == 'sigma2qn_BR') then
     call sigma2qn_BR
  else if (trim(operation) == 'powspec_gauss') then
     call powspec_gauss
  else if (trim(operation) == 'optimization_search') then
     call run_optimization_search
  else if (trim(operation) == 'print_fisher_mat') then
     call print_fisher_mat
  else if (trim(operation) == 'timing') then
     call run_timing
  else if (trim(operation) == 'compute_delta_lnL') then
     call compute_delta_lnL
  else if (trim(operation) == 'fast_par_estimation') then
     call fast_par_estimation
  else if (trim(operation) == 'print_lnL') then
     call print_lnL
  else if (trim(operation) == 'print_lnL_BR') then
     call print_lnL_BR
  else if (trim(operation) == 'print_lnL_old_BR') then
     call print_lnL_old_BR
  else if (trim(operation) == 'sigma2gauss') then
     call sigma2gauss
  else if (trim(operation) == 'output_covmat_slice') then
     call output_covmat_slice
  end if

contains

  subroutine mapcov2gausslike
    implicit none

    integer(i4b)       :: i, j, k, l, m, q, nside, nmaps, npix, lmax, lhigh_T, lhigh_P, ordering, numcomp
    integer(i4b)       :: l1, l2, m1, m2, i1, i2, n_d, n_p, n_h, n
    integer(i4b)       :: ind, ind1, ind2, cov_order, seed, numtemp, nmask, n_t, n_pol
    integer(i4b)       :: lmax_basis, pol, nval, nhigh, num_fg_temp, nside_temp, npix_temp
    real(dp)           :: t1, t2, Tthreshold, Pthreshold, reg(3), scale, W_max_T, W_max_P
    real(dp)           :: BBamp, cvec(3), weight_T, weight_P
    logical(lgt)       :: comp_cov
    character(len=5)   :: itext
    character(len=512) :: mapfile, maskfile, covfile, beamfile, outfile, filename, temp
    character(len=512) :: basis, prefix, clfile, partext, tempcovmat(10)
    real(dp),     allocatable, dimension(:,:)     :: cov, P_harm, C, cl, B, invC, tempcov
    real(dp),     allocatable, dimension(:)       :: W, mask_1d, test3, A_T
    real(dp),     allocatable, dimension(:,:)     :: d, mask, map, alms, V, P_l, V_red, S_mat, P, map_1d
    real(dp),     allocatable, dimension(:,:)     :: Y, Yp, beam, cls, T, buffer, buffer2, P_highl, weights
    real(dp),     allocatable, dimension(:,:)     :: PT, PT_invN_TP
    complex(dpc), allocatable, dimension(:,:,:)   :: alms_cmplx
    real(dp),     pointer,     dimension(:,:)     :: pixwin
    real(dp),     allocatable, dimension(:,:,:)   :: window
    integer(i4b), allocatable, dimension(:)       :: indmap, map2mask, lmax2lhigh
    logical(lgt), allocatable, dimension(:)       :: tempmode, rescale_CMB
    type(planck_rng) :: handle


    if (iargc() < 17) then
       write(*,*) '    Pre-process low-l likelihood inputs from map and covariance matrix'
       write(*,*) '    Options:  [mapfile] [maskfile] [covfile] [compressed cov T/F] [beamfile] [clfile] [lmax_data]'
       write(*,*) '              [lmax_basis_T] [lmax_basis_P] [basis] [T threshold] [P threshold]'
       write(*,*) '              [BB amp] [T reg] [P reg] [seed] [outprefix] [template1] [rescale CMB]'
       write(*,*) '              [template2] [tempcovmat/none] [rescale CMB]...'
       write(*,*) '        where [basis]     = {pixel, pseudoalm, eigen_N, eigen_StoN, eigen_S+N}'
       write(*,*) '              [threshold] = basis dependent accuracy threshold'
       write(*,*) '        If T_reg or P_reg are negative, the corresponding covar blocks are nulled'
       write(*,*) '        before adding the regularization noise in quadrature'
       stop
    end if

    unit = 58
    call getarg(2,mapfile)
    call getarg(3,maskfile)
    call getarg(4,covfile)
    call getarg(5,temp)	
    read(temp,*) comp_cov
    call getarg(6,beamfile)
    call getarg(7,clfile)
    call getarg(8,temp)
    read(temp,*) lmax
    call getarg(9,temp)
    read(temp,*) lhigh_T
    call getarg(10,temp)
    read(temp,*) lhigh_P
    call getarg(11,basis)
    call getarg(12,temp)
    read(temp,*) Tthreshold
    call getarg(13,temp)
    read(temp,*) Pthreshold
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
    beam(:,1)   = beam(:,1) * pixwin(0:lmax,1) ! Only pixwin in T
    beam(:,2)   = beam(:,2) * pixwin(0:lmax,2) ! Only pixwin in T
    beam(:,3)   = beam(:,3) * pixwin(0:lmax,2) ! Only pixwin in T
    do l = 0, lmax
       do i = 1, nmaps
          do j = 1, nmaps
             window(l,i,j) = beam(l,i)*beam(l,j)
          end do
       end do
    end do

    ! Read covariance matrix
    if (comp_cov) then
       allocate(cov(nval,nval))
       call comm_read_covmat(covfile, nmaps, cov, scale)
    else	  
       allocate(cov(n_p,n_p))
       call comm_read_covmat(covfile, nmaps, cov, scale)
       do j = 1, nmaps
          if (reg(j) < 0) then
             cov((j-1)*npix+1:j*npix,:) = 0.d0
             cov(:,(j-1)*npix+1:j*npix) = 0.d0
          end if
       end do
    end if
    map_1d = map_1d * scale ! Assume map is in same units as covariance matrix

    ! Add regularization noise to covariance matrix; DOES NOT YET SUPPORT COMPRESSED COVMATS
    allocate(W(nval))
    call get_eigenvalues(cov,W)
    write(*,*) maxval(W)/minval(W)
    write(*,*) reg
    if (any(reg > 0.d0)) then
       ind = 1
       do j = 1, nmaps
          do i = 1, npix
             cov(ind,ind) = cov(ind,ind) + reg(j)**2
             do k = 1, n_d
                map_1d(ind,k)  = map_1d(ind,k)  + abs(reg(j)) * rand_gauss(handle)
             end do
             ind          = ind+1
          end do
       end do
   end if
    call get_eigenvalues(cov,W)
    write(*,*) maxval(W)/minval(W)
   deallocate(W)
   !stop

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
    num_fg_temp = (iargc()-18)/3
    numtemp     = num_fg_temp * npix_temp
    if (n_t > 0) numtemp = numtemp+4 
    allocate(T(n_p,numtemp), rescale_CMB(numtemp))
    T           = 0.d0
    rescale_CMB = .false.
    if (num_fg_temp > 0) then
       do k = 1, num_fg_temp
          call getarg(18+3*k-2,filename)
          call getarg(18+3*k-1,partext)
          call getarg(18+3*k,  tempcovmat(k))
          read(partext,*) rescale_CMB(k)
          call read_map(filename, map)
          do i = 1, nmaps
             call convert_ring2nest(nside, map(:,i))
             do j = 1, npix_temp
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
    end if

    ! Subtract best-fit templates from data, and add uncertainty to covariance matrix
    if (numtemp > 0) then
       allocate(invC(nval,nval), PT(nval,numtemp), PT_invN_TP(numtemp,numtemp), A_T(numtemp))
!!$       PT = T(map2mask,:)
!!$       PT_invN_TP = matmul(transpose(PT), PT)
!!$       call invert_matrix(PT_invN_TP)
       allocate(S_mat(nval,nval))
       call get_basis_S(0.d0, 1.d0, cls, beam, Y(map2mask,:), S_mat)
       if (comp_cov) then
          S_mat = S_mat + cov
       else
          S_mat = S_mat + cov(map2mask,map2mask)
       end if
       call invert_matrix(S_mat)
       S_mat(1:n_t,n_t+1:nval) = 0.d0
       S_mat(n_t+1:nval,1:n_t) = 0.d0
       PT = T(map2mask,:)
       PT_invN_TP = matmul(transpose(PT), matmul(S_mat,PT))
       call invert_singular_matrix(PT_invN_TP, 1.d-12)

       do k = 1, n_d
           A_T = matmul(PT_invN_TP, matmul(transpose(PT), matmul(S_mat,map_1d(map2mask,k))))
!          A_T = matmul(PT_invN_TP, matmul(transpose(PT), map_1d(map2mask,k)))
	   
	   ! 44 GHz, BolPol
           !A_T(1) = 0.204
           !A_T(2) = 0.0048

	   ! 44 GHz; beta_s = -3.1, beta_d = 1.6, T_d = 18K
	   !A_T(1) = 0.256d0
	   !A_T(2) = 0.035d0

	   ! 70 GHz; beta_s = -3.1, beta_d = 1.6, T_d = 18K
	   !A_T(1) = 0.0656d0
	   !A_T(2) = 0.0763d0

          write(*,*) 'Subtracted template amplitudes for map ', k
          weight_T = 0.d0
          weight_P = 0.d0
          do i = 1, numtemp
             write(*,*) '   Template ', i, ' = ', real(A_T(i),sp), '+/-', real(sqrt(Pt_invN_TP(i,i)),sp)
             if (rescale_CMB(i) .and. any(PT(1:n_t,i) /= 0.d0))      weight_T = weight_T + A_t(i)
             if (rescale_CMB(i) .and. any(PT(n_t+1:nval,i) /= 0.d0)) weight_P = weight_P + A_t(i)
          end do
          write(*,*) '   CMB weights = ', real(1.d0-weight_T,sp), real(1.d0-weight_P,sp)
       
          if (k == 1) then
             do j = 1, nmaps
                map(:,j) = map_1d((j-1)*npix+1:j*npix,k) * mask_1d((j-1)*npix+1:j*npix)
             end do
             call write_map3(trim(prefix)//'_rawmap.fits', map)
          end if

          map_1d(map2mask,k)    = map_1d(map2mask,k) - matmul(PT, A_T)
          map_1d(1:npix,k)      = map_1d(1:npix,k)     / (1.d0-weight_T)
          map_1d(npix+1:+n_p,k) = map_1d(npix+1:n_p,k) / (1.d0-weight_P)

          if (k == 1) then
             do j = 1, nmaps
                map(:,j) = map_1d((j-1)*npix+1:j*npix,k) * mask_1d((j-1)*npix+1:j*npix)
             end do
             call write_map3(trim(prefix)//'_cleanmap.fits', map)
          end if
       end do
       deallocate(S_mat)

       do i = 1, numtemp
          if (comp_cov) then
             cov = cov + &
                  & Pt_invN_TP(i,i) * matmul(PT(:,i:i), transpose(PT(:,i:i)))
	  else
             cov(map2mask,map2mask) = cov(map2mask,map2mask) + &
                  & Pt_invN_TP(i,i) * matmul(PT(:,i:i), transpose(PT(:,i:i)))
          end if
          if (trim(tempcovmat(i)) /= 'none' .and. i <= num_fg_temp) then
             allocate(tempcov(n_p,n_p))
             call comm_read_covmat(tempcovmat(i), nmaps, tempcov, scale)
             tempcov(1:npix,:) = 0.d0 ! Don't include temperature correlations in covariance
             tempcov(:,1:npix) = 0.d0
	     if (comp_cov) then
	        cov = cov + A_T(i)**2 * tempcov(map2mask,map2mask)
	     else
	        cov(map2mask,map2mask) = cov(map2mask,map2mask) + A_T(i)**2 * tempcov(map2mask,map2mask)
             end if
             deallocate(tempcov)
          end if
       end do
       if (weight_T /= 0.d0) then
          cov(1:npix,:)     = cov(1:npix,:)     / (1.d0-weight_T)
          cov(:,1:npix)     = cov(:,1:npix)     / (1.d0-weight_T)
       end if
       if (weight_P /= 0.d0) then
          cov(npix+1:n_p,:) = cov(npix+1:n_p,:) / (1.d0-weight_P)
          cov(:,npix+1:n_p) = cov(:,npix+1:n_p) / (1.d0-weight_P)
       end if
    end if

    if (nmaps == 1) then
       write(*,fmt='(a,f8.3)') '      I sky fractions = ', sum(mask(:,1))/npix
    else
       write(*,fmt='(a,3f8.3)') '      I/Q/U sky fractions = ', sum(mask(:,1))/npix, &
            & sum(mask(:,2))/npix, sum(mask(:,3))/npix
    end if


    ! Compute high-l projection operator
    write(*,*) 'nval = ', nval
    allocate(P_highl(nval,nval))
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

       if (comp_cov) then
          P = cov(map2mask,map2mask)
       else
          P = cov
       end if
       call invert_matrix(P)

    else if (trim(basis) == 'eigen_StoN') then

       if (comp_cov) then
          P = cov
       else
          P = cov(map2mask,map2mask)
       end if
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
       if (comp_cov) then
          P = cov + S_mat
       else
          P = cov(map2mask,map2mask) + S_mat
       end if
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
    if (comp_cov) then
       call dgemm('N','N',nval,n,nval,1.d0,cov,nval,V_red,nval,0.d0,buffer,nval)
    else 
       call dgemm('N','N',nval,n,nval,1.d0,cov(map2mask,map2mask),nval,V_red,nval,0.d0,buffer,nval)
    end if
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

    write(*,*) n
    deallocate(W)
    allocate(W(n))
    call get_eigenvalues(C,W)

    open(68,file='W.dat')
    do i = 1, n
       write(68,*) W(i)/maxval(W)
    end do
    close(68)
    
    write(*,*) maxval(W)/minval(W)
    deallocate(W)

    ! Output datafile
    call write_lowl_datafile(trim(prefix)//'.fits', d, C, beam, P_harm, mapfile, maskfile, &
         & covfile, beamfile, clfile, basis, lhigh_T, lhigh_P, Tthreshold, Pthreshold)

  end subroutine mapcov2gausslike

  subroutine coadd_dataset
    implicit none

    integer(i4b)       :: i, j, k, l, n, unit, nside, npix, nmaps, m
    character(len=512) :: filelist, prefix, mapname, covname, maskfile
    real(dp)           :: scale, t1, t2
    real(dp), allocatable, dimension(:)   :: map_1d, map_tot
    real(dp), allocatable, dimension(:,:) :: map_in, map, cov_in, cov_tot, mask

    call getarg(2,filelist)
    call getarg(3,prefix)

    unit = comm_getlun()
    open(unit,file=trim(filelist))
    read(unit,*) nside
    read(unit,*) nmaps
    read(unit,*) m
    npix = 12*nside**2
    n    = nmaps*npix

    allocate(map(0:npix-1,nmaps), map_1d(n), map_tot(n), cov_tot(n,n), cov_in(n,n), mask(0:npix-1,1))
    map_tot = 0.d0
    cov_tot = 0.d0
    do i = 1, m
       call wall_time(t1)
       read(unit,*) mapname, covname, maskfile
       write(*,*) 'Processing ', trim(mapname)
       call comm_read_covmat(covname, nmaps, cov_in, scale)

       call read_map(mapname, map)
       call read_map(maskfile, mask)

       where (map > -1.637d30) 
          map = map * scale
       end where
       do j = 1, nmaps
          where (mask(:,1) < 0.5) 
             map(:,j) = -1.6375d30
          end where       
       end do

       do j = 1, nmaps
          map_1d((j-1)*npix+1:j*npix) = map(:,j)
       end do
       do j = 1, n
          if (map_1d(j) < -1.637d30) then
	     cov_in(:,j) = 0.d0
	     cov_in(j,:) = 0.d0
             cov_in(j,j) = 1d12
             map_1d(j)    = 0.d0
          end if
       end do
       call invert_matrix(cov_in)
       cov_tot = cov_tot + cov_in
       map_tot = map_tot + matmul(cov_in, map_1d)
       call wall_time(t2)
       write(*,*) '     Wall time = ', real(t2-t1,sp)
    end do
    call invert_matrix(cov_tot)
    map_tot = matmul(cov_tot, map_tot)

    do j = 1, nmaps
       map(:,j) = map_tot((j-1)*npix+1:j*npix)
    end do

    call write_map3(trim(prefix)//'_map.fits', map)

    ! Write to unformatted F90 file with appropriate header entries
    open(unit,file=trim(prefix)//'_cov.unf',form='unformatted')
    write(unit) n
    write(unit) 1
    write(unit) 1
    do i = 1, n
       write(unit) cov_tot(i,:)
    end do
    write(unit) .false.
    close(unit)

    k = 1
    do j = 1, nmaps
       do i = 0, npix-1
          map(i,j) = sqrt(cov_tot(k,k))
          k        = k+1
       end do
    end do

    call write_map3(trim(prefix)//'_rms.fits', map)
    
  end subroutine coadd_dataset


  subroutine fast_par_estimation
    implicit none

    character(len=256) :: infile, maxlnLfile, temp, clfile,filename, filelist, outfile
    real(dp)           :: lnL_max, lnL_fid, lnL0, val, vals(1000), lnLs(1000), lnL2(1000), delta_lnL
    real(dp)           :: peak, upper, lower, limit, red_chisq, chisq
    integer(i4b)       :: nspec, lmin_eval, lmax_eval, l, i, n
    integer(i4b)       :: ierr, ndof
    integer(i4b)       :: unit, unit_out, nspline
    logical(lgt)       :: exist, flag(6), output_detrate
    real(dp),     allocatable, dimension(:,:)    :: cls_maxlnL, cls_in, cls0
    real(dp),     allocatable, dimension(:)      :: cl_cond, p_cond, eta
    real(dp),     allocatable, dimension(:)      :: P_spline, x_spline
    real(dp),     allocatable, dimension(:)      :: condnum
    real(dp),     allocatable, dimension(:)      :: cl_spline, lnL_spline
    integer(i4b), allocatable, dimension(:)      :: indmap
    logical(lgt), allocatable, dimension(:)      :: converged
    character(len=1), dimension(6) :: status
    character(len=2), dimension(6) :: stext = ['TT','TE','TB','EE','EB','BB']
    type(planck_rng) :: handle

    if (iargc() /= 14) then
       write(*,*) '    Estimate parameters given tabulated spectra'
       write(*,*) '      Usage: comm_like_tools fast_par_estimation [infofile] [maxLn clfile]'
       write(*,*) '                 [spectrum list] [lmin] [lmax] [6 spectrum flags] '
       write(*,*) '                 [output detection level] [outfile]'

       stop
    end if

    call getarg(2, infile)
    call getarg(3, maxlnlfile)
    call getarg(4, filelist)
    call getarg(5, temp)
    read(temp,*) lmin_eval
    call getarg(6, temp)
    read(temp,*) lmax_eval
    do i = 1, 6
       call getarg(6+i, temp)
       read(temp,*) flag(i)
    end do
    call getarg(13, temp)
    read(temp,*) output_detrate
    call getarg(14, outfile)
    ierr     = 0
    nspec    = 6
    unit     = comm_getlun()
    unit_out = unit+1

    inquire(file=trim(infile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(infile)
    inquire(file=trim(filelist), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(filelist)

    ! Initialize low-l module
    call comm_lowl_initialize_object(infile)

    ! Read fiducial spectrum
    call read_fiducial_spectrum(maxlnlfile, cls0)             
    call read_fiducial_spectrum(maxlnlfile, cls_fid)             
    !cls0 = 0.d0

    open(unit,file=trim(filelist))
    n = 0
    do while (.true.)
       read(unit,*,end=91) val, filename
       n       = n+1
       vals(n) = val
       call read_fiducial_spectrum(filename, cls_maxlnL)
       cls_fid = cls0
       do i = 1, 6
          if (.not. flag(i)) cycle
          write(*,*) n, i
          do l = lmin_eval, lmax_eval
             cls_fid(l,i) = cls_maxlnL(l,i)
          end do
       end do
       lnLs(n) = comm_lowl_compute_lnL(cls=cls_fid, ierr=ierr, enforce_pos_def=.false., red_chisq=red_chisq, chisq=chisq)
       write(*,*) n, trim(filename), lnLs(n), chisq, red_chisq
    end do
91  close(unit)
    
    nspline = 10000
    allocate(x_spline(nspline), P_spline(nspline))
    call spline(vals(1:n), lnLs(1:n), 1.d30, 1.d30, lnL2(1:n))
    do i = 1, nspline
       x_spline(i) = vals(1) + (vals(n)-vals(1)) / (nspline-1.d0) * (i-1)
       P_spline(i) = splint(vals(1:n), lnLs(1:n), lnL2(1:n), x_spline(i))
    end do
    
    P_spline = exp(P_spline-maxval(P_spline))
    P_spline = P_spline / sum(P_spline) / ((vals(n)-vals(1))/(nspline-1))
    open(unit,file=trim(outfile))
    do i = 1, nspline
       write(unit,*) x_spline(i), P_spline(i)
    end do
    close(unit)

    call compute_asymmetric_errors(x_spline, P_spline, peak, upper, lower)
    write(*,fmt='(a,3f12.6)') 'Estimate = ', peak, upper, lower

    ! Output upper 2 and 3 sigma limits
    limit = 0.d0
    i     = 0
    do while (limit < 0.954d0)
        i = i+1
        limit = limit + P_spline(i) * ((vals(n)-vals(1))/(nspline-1))
    end do
    write(*,fmt='(a,f12.6)') 'Upper (one-sided) 2-sigma limit = ', x_spline(i)

    do while (limit < 0.997d0)
        i = i+1
        limit = limit + P_spline(i) * ((vals(n)-vals(1))/(nspline-1))
    end do
    write(*,fmt='(a,f12.6)') 'Upper (one-sided) 3-sigma limit = ', x_spline(i)

    if (output_detrate) then
       cls_fid = cls0
       do i = 1, 6
          if (.not. flag(i)) cycle
          do l = lmin_eval, lmax_eval
             cls_fid(l,i) = 0.d0
          end do
       end do
       delta_lnL = comm_lowl_compute_lnL(cls=cls_fid, ierr=ierr, enforce_pos_def=.false., red_chisq=red_chisq) - maxval(lnLs(1:n))
       write(*,fmt='(a,f8.2,a)') 'Detection level = ', sqrt(-2d0*delta_lnL), ' sigma (Gaussianized)'
       write(*,fmt='(a,f8.2)') 'lnL(max) = ', maxval(lnLs(1:n))
       write(*,fmt='(a,2f8.2)') 'lnL(0)   = ', maxval(lnLs(1:n)) + delta_lnL, red_chisq
    end if

  end subroutine fast_par_estimation
  

  subroutine compute_delta_lnL
    implicit none

    character(len=256) :: infile, maxlnLfile, temp, clfile,filename
    real(dp)           :: lnL_max, lnL_fid, lnL0, val, vals(1000), lnLs(1000), lnL2(1000)
    integer(i4b)       :: nspec, lmin_eval, lmax_eval, l, i, n
    integer(i4b)       :: ierr, ndof
    integer(i4b)       :: unit, unit_out, nspline
    logical(lgt)       :: exist, flag(6)
    real(dp),     allocatable, dimension(:,:)    :: cls_maxlnL, cls_in, cls0
    real(dp),     allocatable, dimension(:)      :: cl_cond, p_cond, eta
    real(dp),     allocatable, dimension(:)      :: P_spline, x_spline
    real(dp),     allocatable, dimension(:)      :: condnum
    real(dp),     allocatable, dimension(:)      :: cl_spline, lnL_spline
    integer(i4b), allocatable, dimension(:)      :: indmap
    logical(lgt), allocatable, dimension(:)      :: converged
    character(len=1), dimension(6) :: status
    character(len=2), dimension(6) :: stext = ['TT','TE','TB','EE','EB','BB']
    type(planck_rng) :: handle

    if (iargc() /= 12) then
       write(*,*) '    Compute delta lnL (and chisq) from low-l likelihood'
       write(*,*) '      Usage: comm_like_tools compute_delta_lnL [infofile] '
       write(*,*) '                 [max lnL file] [theory C_l file] [lmin] [lmax] [6 spectrum flags]'

       stop
    end if

    call getarg(2, infile)
    call getarg(3, maxlnLfile)
    call getarg(4, clfile)
    call getarg(5, temp)
    read(temp,*) lmin_eval
    call getarg(6, temp)
    read(temp,*) lmax_eval
    do i = 1, 6
       call getarg(6+i, temp)
       read(temp,*) flag(i)
    end do
    ierr     = 0
    nspec    = 6
    unit     = comm_getlun()
    unit_out = unit+1

    inquire(file=trim(infile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(infile)
    inquire(file=trim(maxlnLfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(maxlnLfile)
    inquire(file=trim(clfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(clfile)

    ! Initialize low-l module
    call comm_lowl_initialize_object(infile)

    if (trim(clfile) == 'camblist.txt') then
       open(unit,file=trim(clfile))
       n = 0
       do while (.true.)
          read(unit,*,end=91) val, filename
          if (.not. allocated(cls0)) then
             call read_fiducial_spectrum(filename, cls0)             
             call read_fiducial_spectrum(filename, cls_fid)             
          end if
          n       = n+1
          vals(n) = val
          call read_fiducial_spectrum(filename, cls_maxlnL)
          cls_fid = cls0
          do i = 1, 6
             if (.not. flag(i)) cycle
             do l = lmin_eval, lmax_eval
                cls_fid(l,i) = cls_maxlnL(l,i)
             end do
          end do
          lnLs(n) = comm_lowl_compute_lnL(cls=cls_fid, ierr=ierr, enforce_pos_def=.false.)
          write(*,*) vals(n), lnLs(n)
       end do
91     close(unit)

       nspline = 1000
       allocate(x_spline(nspline), P_spline(nspline))
       call spline(vals(1:n), lnLs(1:n), 1.d30, 1.d30, lnL2(1:n))
       do i = 1, nspline
          x_spline(i) = vals(n) / (nspline-1.d0) * (i-1)
          P_spline(i) = splint(vals(1:n), lnLs(1:n), lnL2(1:n), x_spline(i))
       end do

       P_spline = exp(P_spline-maxval(P_spline))
       open(unit,file='lnL.dat')
       do i = 1, n
          write(unit,*) vals(i), lnLs(i)
       end do
       close(unit)

    else

       call read_fiducial_spectrum(clfile, cls_fid)
       call read_fiducial_spectrum(clfile, cls_maxlnL)
       call read_fiducial_spectrum(clfile, cls0)
       call read_fiducial_spectrum(maxlnLfile, cls_in)
       do i = 1, 6
          if (.not. flag(i)) cycle
          do l = lmin_eval, lmax_eval
             cls_maxlnL(l,i) = cls_in(l,i)
             cls0(l,i)       = 0.d0
          end do
       end do
       ndof    = count(cls_maxlnL /= cls_fid)
       
       lnL_max = comm_lowl_compute_lnL(cls=cls_maxlnL, ierr=ierr, enforce_pos_def=.false.)
       lnL_fid = comm_lowl_compute_lnL(cls=cls_fid,    ierr=ierr, enforce_pos_def=.false.)
       lnL0    = comm_lowl_compute_lnL(cls=cls0,       ierr=ierr, enforce_pos_def=.false.)
       
       write(*,*)
       write(*,*) 'lnL_max           = ', lnL_max
       write(*,*) 'lnL_fid           = ', lnL_fid
       write(*,*) 'lnL0              = ', lnL0
       write(*,*) 'Delta lnL vs fid  = ', lnL_fid-lnL_max
       write(*,*) 'chisq vs fid      = ', -2.d0*(lnL_fid-lnL_max)
       write(*,*) 'Delta lnL vs null = ', lnL0-lnL_max
       write(*,*) 'chisq vs null     = ', -2.d0*(lnL0-lnL_max)
       write(*,*) 'ndof              = ', ndof
    end if


  end subroutine compute_delta_lnL



  subroutine print_lnL
    implicit none

    character(len=256) :: infile, maxlnLfile, temp, clfile,filename
    real(dp)           :: lnL_max, lnL_fid, lnL0, val, vals(1000), lnLs(1000), lnL2(1000), lnL
    real(dp)           :: chisq, red_chisq
    integer(i4b)       :: nspec, lmin_eval, lmax_eval, l, i, n
    integer(i4b)       :: ierr, ndof
    integer(i4b)       :: unit, unit_out, nspline
    logical(lgt)       :: exist, flag(6)
    real(dp),     allocatable, dimension(:,:)    :: cls_maxlnL, cls_in, cls0
    real(dp),     allocatable, dimension(:)      :: cl_cond, p_cond, eta
    real(dp),     allocatable, dimension(:)      :: P_spline, x_spline
    real(dp),     allocatable, dimension(:)      :: condnum
    real(dp),     allocatable, dimension(:)      :: cl_spline, lnL_spline
    integer(i4b), allocatable, dimension(:)      :: indmap
    logical(lgt), allocatable, dimension(:)      :: converged
    character(len=1), dimension(6) :: status
    character(len=2), dimension(6) :: stext = ['TT','TE','TB','EE','EB','BB']
    type(planck_rng) :: handle

    if (iargc() /= 3) then
       write(*,*) '    Print lnL from low-l likelihood'
       write(*,*) '      Usage: comm_like_tools print_lnL [infofile] [C_l file]'
       stop
    end if

    call getarg(2, infile)
    call getarg(3, clfile)
    unit     = comm_getlun()

    inquire(file=trim(infile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(infile)
    inquire(file=trim(clfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(clfile)

    ! Initialize low-l module
    call comm_lowl_initialize_object(infile)

    call read_fiducial_spectrum(clfile, cls_fid)
    lnL = comm_lowl_compute_lnL(cls=cls_fid, chisq=chisq, red_chisq=red_chisq, ierr=ierr)	

    write(*,*) 'filename  = ', trim(clfile)
    write(*,*) 'lnL       = ', lnL
    write(*,*) 'chisq     = ', chisq
    write(*,*) 'chisq_red = ', red_chisq
    write(*,*) 'ierr      = ', ierr

  end subroutine print_lnL


  subroutine print_lnL_BR
    implicit none

    character(len=256) :: infile, maxlnLfile, temp, clfile,filename, brfile, BRtype
    real(dp)           :: lnL_max, lnL_fid, lnL0, val, vals(1000), lnLs(1000), lnL2(1000), lnL
    integer(i4b)       :: nspec, lmin_eval, lmax_eval, l, i, n, lmin, lmax, dl
    integer(i4b)       :: ierr, ndof
    integer(i4b)       :: unit, unit_out, nspline
    logical(lgt)       :: exist, flag(6)
    real(dp),     allocatable, dimension(:,:)    :: cls_maxlnL, cls_in, cls0
    real(dp),     allocatable, dimension(:)      :: cl_cond, p_cond, eta
    real(dp),     allocatable, dimension(:)      :: P_spline, x_spline
    real(dp),     allocatable, dimension(:)      :: condnum
    real(dp),     allocatable, dimension(:)      :: cl_spline, lnL_spline
    integer(i4b), allocatable, dimension(:)      :: indmap
    logical(lgt), allocatable, dimension(:)      :: converged
    character(len=1), dimension(6) :: status
    character(len=2), dimension(6) :: stext = ['TT','TE','TB','EE','EB','BB']
    type(planck_rng) :: handle

    if (iargc() < 4) then
       write(*,*) '    Print lnL from Blackwell-Rao likelihood'
       write(*,*) '      Usage: comm_like_tools print_lnL_BR [BRtype] [infofile] [C_l file]'
       write(*,*) '                 [lmin (for gauss)] [lmax (for gauss)] [dl_band (for gauss)]'
       stop
    end if

    call getarg(2, BRtype)
    call getarg(3, brfile)
    call getarg(4, clfile)
    unit     = comm_getlun()
    if (trim(BRtype) == 'gauss') then
       call getarg(5, temp)
       read(temp,*) lmin
       call getarg(6, temp)
       read(temp,*) lmax
       call getarg(7, temp)
       read(temp,*) dl
    end if

    inquire(file=trim(brfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: BR file does not exist = ', trim(brfile)
    inquire(file=trim(clfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: Cl file does not exist = ', trim(clfile)

    call read_fiducial_spectrum(clfile, cls_fid)

    ierr = 0
    if (trim(BRtype) == 'exact') then
       call comm_br_initialize_object(brfile)
       lnL = comm_br_compute_lnL(cls_fid, ierr)
    else if (trim(BRtype) == 'gauss') then
       call comm_gauss_br_initialize_object(brfile, lmin, lmax, dl)
       lnL = comm_gauss_br_compute_lnL(cls_fid(2:,1))
    end if


    write(*,*) 'sigmafile = ', trim(brfile)
    write(*,*) 'clfile    = ', trim(clfile)
    if (trim(BRtype) == 'gauss') then
       write(*,*) 'lmin      = ', lmin
       write(*,*) 'lmax      = ', lmax
       write(*,*) 'delta_l   = ', dl
    end if
    write(*,*) 'lnL       = ', lnL
    write(*,*) 'ierr      = ', ierr

  end subroutine print_lnL_BR


  subroutine print_lnL_old_BR
    implicit none

    character(len=256) :: infile, maxlnLfile, temp, clfile,filename, brfile, sigmafile
    real(dp)           :: lnL_max, lnL_fid, lnL0, val, vals(1000), lnLs(1000), lnL2(1000), lnL
    integer(i4b)       :: nspec, lmin_eval, lmax_eval, l, i, n
    integer(i4b)       :: ierr, ndof
    integer(i4b)       :: unit, unit_out, nspline
    integer(i4b)       :: lmin, lmax, firstsample, lastsample, firstchain, lastchain, thinstep
    logical(lgt)       :: exist, flag(6)
    real(dp),     allocatable, dimension(:,:)    :: cls_maxlnL, cls_in, cls0
    real(dp),     allocatable, dimension(:)      :: cl_cond, p_cond, eta
    real(dp),     allocatable, dimension(:)      :: P_spline, x_spline
    real(dp),     allocatable, dimension(:)      :: condnum
    real(dp),     allocatable, dimension(:)      :: cl_spline, lnL_spline
    integer(i4b), allocatable, dimension(:)      :: indmap
    logical(lgt), allocatable, dimension(:)      :: converged
    character(len=1), dimension(6) :: status
    character(len=2), dimension(6) :: stext = ['TT','TE','TB','EE','EB','BB']
    type(planck_rng) :: handle

    if (iargc() /= 10) then
       write(*,*) '    Print lnL from old Blackwell-Rao likelihood'
       write(*,*) '      Usage: comm_like_tools print_lnL_old_BR [sigmafile] [clfile] '
       write(*,*) '                    [lmin] [lmax] [firstchain] [lastchain] '
       write(*,*) '                    [firstsample] [lastsample] [thinstep] '
       stop
    end if

    call getarg(2, sigmafile)
    call getarg(3, clfile)
    call getarg(4, temp)
    read(temp,*) lmin
    call getarg(5, temp)
    read(temp,*) lmax
    call getarg(6, temp)
    read(temp,*) firstchain
    call getarg(7, temp)
    read(temp,*) lastchain
    call getarg(8, temp)
    read(temp,*) firstsample
    call getarg(9, temp)
    read(temp,*) lastsample
    call getarg(10, temp)
    read(temp,*) thinstep

    unit     = comm_getlun()

    inquire(file=trim(sigmafile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(sigmafile)
    inquire(file=trim(clfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(clfile)

    ! Initialize BR module
    call comm_br_old_initialize_object(sigmafile, clfile, lmin, lmax, firstchain, lastchain, &
       & firstsample, lastsample, thinstep)

    call read_fiducial_spectrum(clfile, cls_fid)
    lnL = comm_br_old_compute_lnL(cls_fid(2:lmax,1))

    write(*,*) 'sigmafile   = ', trim(sigmafile)
    write(*,*) 'clfile      = ', trim(clfile)
    write(*,*) 'lmin        = ', lmin
    write(*,*) 'lmax        = ', lmax
    write(*,*) 'firstchain  = ', firstchain
    write(*,*) 'lastchain   = ', lastchain
    write(*,*) 'firstsample = ', firstsample
    write(*,*) 'lastsample  = ', lastsample
    write(*,*) 'thinstep    = ', thinstep
    write(*,*) 'lnL         = ', lnL
    write(*,*) 'ierr        = ', ierr

  end subroutine print_lnL_old_BR
  

  subroutine sigma2cl_BR
    implicit none

    character(len=256) :: brfile, clfile, outfile, outprefix
    character(len=2)   :: spectext
    real(dp)           :: cl_min, cl_max, sigma, top, left, right, int, tot_int, cl_in(4), dx
    real(dp)           :: peak, lower, upper
    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep, nmaps, ind
    integer(i4b)       :: numbin, pos(3), ell, b_min, b_max, bincount, id_bin, nspline
    integer(i4b)       :: nspec, lmin, lmax, numchain, numiter, i, j, k, l, m, n, p, spec, ierr, id
    integer(i4b)       :: ind1, ind2
    logical(lgt)       :: exist
    character(len=1), dimension(6)   :: status
    real(dp), allocatable, dimension(:,:) :: cls, cls0
    real(dp), allocatable, dimension(:)   :: cl, lnL, cl_spline, lnL_spline, lnL2
    character(len=2), dimension(6) :: stext = ['TT','TE','TB','EE','EB','BB']
    real(dp), dimension(6) :: clmin = [0.d0, -2.d1, -1.d0, 0.d0, -1.d0, 0.d0]
    real(dp), dimension(6) :: clmax = [1.d4,  2.d1,  1.d0, 1.d0,  1.d0, 0.01d0]

    if (iargc() < 3) then
       write(*,*) '    Compute spectrum with Blackwell-Rao'
       write(*,*) '      Usage: comm_like_tools sigma2cl_BR [BR file] [fid cl file] [outprefix]'
       stop
    end if

    call getarg(2, brfile)
    call getarg(3, clfile)
    call getarg(4, outprefix)

    inquire(file=trim(brfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(brfile)
    inquire(file=trim(clfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(clfile)

    call read_fiducial_spectrum(clfile, cls0)
    lmax = size(cls0,1)-1
    allocate(cls(0:lmax,6))

    ! Initialize BR module
    call comm_br_initialize_object(brfile, handle=id)

!!$    open(58,file='cls.dat')
!!$    cls = 0.d0
!!$    do l = 2, 153
!!$       read(58,*) i, cls(l,1)
!!$       cls(l,1) = cls(l,1) * l*(l+1)/2.d0/pi
!!$    end do
!!$    peak = comm_br_compute_lnL(cls, ierr)
!!$    write(*,*) peak, ierr
!!$    stop

    ! Loop over all free parameters, and output slices
    numbin  = 1000
    nspline = 100000
    allocate(cl(numbin), lnL(numbin), lnL2(numbin), cl_spline(nspline), lnL_spline(nspline))
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

                   call comm_br_initialize_object(BRfile, handle=id_bin, &
                        & ell_min=lmin, ell_max=lmax)
                   
                   cls = cls0
                   lnL     = log(1.d-30)
                   do l = 1, numbin
                      cl(l)  = cl_min + (cl_max-cl_min)/(numbin-1) * (l-1)
                      cls(lmin:lmax,ind) = cl(l)
                      lnL(l) = comm_br_compute_lnL(cls, ierr, handle=id_bin, check_posdef=.false.)
                      write(*,*) cl(l), lnL(l)
                   end do
                   ind1 = 1
                   do while ((lnL(ind1) < -1.d29 .or. lnL(ind1) /= lnL(ind1)) .and. ind1 < numbin)
                      ind1 = ind1+1
                   end do
                   ind2 = numbin
                   do while ((lnL(ind2) < -1.d29 .or. (lnL(ind2) /= lnL(ind2))) .and. ind2 > ind1+1)
                      ind2 = ind2-1
                   end do
                   write(*,*) ind1, ind2
                   write(*,*) cl(ind1), cl(ind2)
                   write(*,*) cl(ind1:ind2)
                   write(*,*) lnL(ind1:ind2)
                   !lnL = exp(lnL-maxval(lnL))
                   !lnL = lnL / sum(lnL) / ((cl_max-cl_min)/(numbin-1))
                   call spline(cl(ind1:ind2), lnL(ind1:ind2), 1.d30, 1.d30, lnL2(ind1:ind2))
                   do l = 1, nspline
                      cl_spline(l)  = cl_min + (cl_max-cl_min)/(nspline-1) * (l-1)
                      if (cl_spline(l) < cl(ind1) .or. cl_spline(l) > cl(ind2)) then
                         lnL_spline(l) = -1.d30
                      else
                         lnL_spline(l) = splint(cl(ind1:ind2), lnL(ind1:ind2), lnL2(ind1:ind2), &
                              & cl_spline(l))
                      end if
                   end do
                   lnL_spline = exp(lnL_spline-maxval(lnL_spline))
                   !write(*,*) lnL_spline
                   lnL_spline = lnL_spline / sum(lnL_spline) / ((cl_max-cl_min)/(nspline-1))

                   !open(58,file='slice.dat')
                   !do l = 1, nspline
                   !   write(58,*) cl_spline(l), lnL_spline(l)
                   !end do
                   !close(58)

                   call compute_asymmetric_errors(cl_spline, lnL_spline, peak, upper, lower)
                   
                   outfile = trim(outprefix) // '_' // spectext // '.dat'
                   if (i == 1) then
                      open(unit,file=trim(outfile))
                      write(unit,*) '# l_c   lmin  lmax   D_l      +dD_l     -dD_l'
                   else
                      open(unit,file=trim(outfile),position='append')
                   end if
                   write(unit,fmt='(f6.1,2i6,3f10.3)') 0.5*(lmin+lmax), lmin, lmax, &
                        & peak, upper, lower
                   write(*,fmt='(f6.1,2i6,3f10.3)') 0.5*(lmin+lmax), lmin, lmax, &
                        & peak, upper, lower
                   close(unit)

                   call comm_br_deallocate_object(handle=id_bin)
                end if
             end do
          end do
       end do
    end do

  end subroutine sigma2cl_BR

  subroutine sigma2cl_gauss_BR
    implicit none

    character(len=256) :: brfile, clfile, outfile, outprefix, tmp
    character(len=2)   :: spectext
    real(dp)           :: cl_min, cl_max, sigma, top, left, right, int, tot_int, cl_in(4), dx
    real(dp)           :: peak, lower, upper
    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep, nmaps, ind
    integer(i4b)       :: numbin, pos(3), ell, b_min, b_max, bincount, id_bin, nspline
    integer(i4b)       :: nspec, lmin, lmax, lmax_file, numchain, numiter, i, j, k, l, m, n, p, spec, ierr, id
    integer(i4b)       :: ind1, ind2
    logical(lgt)       :: exist
    character(len=1), dimension(6)   :: status
    real(dp), allocatable, dimension(:)   :: cl, lnL, cl_spline, lnL_spline, lnL2, cls
    character(len=2), dimension(6) :: stext = ['TT','TE','TB','EE','EB','BB']
    real(dp), dimension(6) :: clmin = [0.d0, -2.d1, -1.d0, 0.d0, -1.d0, 0.d0]
    real(dp), dimension(6) :: clmax = [1.d4,  2.d1,  1.d0, 1.d0,  1.d0, 0.01d0]

    if (iargc() < 3) then
       write(*,*) '    Compute spectrum with Gaussianized Blackwell-Rao'
       write(*,*) '      Usage: comm_like_tools sigma2cl_BR [BR file] [lmax] [outprefix]'
       stop
    end if

    call getarg(2, brfile)
    call getarg(3, tmp)
    read(tmp,*) lmax_file
    call getarg(4, outprefix)

    inquire(file=trim(brfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(brfile)

    ! Loop over all free parameters, and output slices
    spectext = stext(1)
    numbin   = 1000
    nspline  = 100000
    allocate(cl(numbin), lnL(numbin), lnL2(numbin), cl_spline(nspline), lnL_spline(nspline))
    allocate(cls(2:lmax_file))
    do l = 2, lmax_file
       lmin = l
       lmax = l
       call comm_gauss_br_initialize_object(BRfile, lmin, lmax, 1, handle=id_bin)
       do i = 1, numbin
          cl(i)  = clmin(1) + (clmax(1)-clmin(1))/(numbin-1) * (i-1)
          cls(l) = cl(i)
          lnL(i) = comm_gauss_br_compute_lnL(cls, handle=id_bin)
       end do
       ind1 = 1
       do while ((abs(lnL(ind1)) > 1.d29 .or. lnL(ind1) /= lnL(ind1)) .and. ind1 < numbin)
          ind1 = ind1+1
       end do
       ind2 = numbin
       do while ((abs(lnL(ind2)) > 1.d29 .or. (lnL(ind2) /= lnL(ind2))) .and. ind2 > ind1+1)
          ind2 = ind2-1
       end do
       call spline(cl(ind1:ind2), lnL(ind1:ind2), 1.d30, 1.d30, lnL2(ind1:ind2))
       do i = 1, nspline
          cl_spline(i)  = cl(ind1) + (cl(ind2)-cl(ind1))/(nspline-1) * (i-1)
          lnL_spline(i) = splint(cl(ind1:ind2), lnL(ind1:ind2), lnL2(ind1:ind2), &
               & cl_spline(i))
       end do
       lnL_spline = exp(lnL_spline-maxval(lnL_spline))
       lnL_spline = lnL_spline / sum(lnL_spline) / ((cl(ind2)-cl(ind1))/(nspline-1))

!!$       open(58,file='slice.dat')
!!$       do i = 1, nspline
!!$          write(58,*) cl_spline(i), lnL_spline(i)
!!$       end do
!!$       close(58)
!!$       stop

       call compute_asymmetric_errors(cl_spline, lnL_spline, peak, upper, lower)

       outfile = trim(outprefix) // '_' // spectext // '.dat'
       if (l == 2) then
          open(unit,file=trim(outfile))
          write(unit,*) '# l_c   lmin  lmax   D_l      +dD_l     -dD_l'
       else
          open(unit,file=trim(outfile),position='append')
       end if
       write(unit,fmt='(f6.1,2i6,3f10.3)') 0.5*(lmin+lmax), lmin, lmax, &
            & peak, upper, lower
       write(*,fmt='(f6.1,2i6,3f10.3)') 0.5*(lmin+lmax), lmin, lmax, &
            & peak, upper, lower
       close(unit)

       call comm_gauss_br_deallocate_object(handle=id_bin)
    end do

  end subroutine sigma2cl_gauss_BR

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
!       do i = 1, 1
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
                      lnL(l) = comm_br_compute_lnL(cls(:,1:1), ierr, handle=id, P_id=p)
                      !write(*,*) cl(l), lnL(l)
                      lnL_max = max(lnL_max, lnL(l))
!                      if (lnL_max-lnL(l) > 10 .and. l>=10) exit
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
    integer(i4b)       :: nspec, lmin, lmax, numchain, numiter, g, i, j, k, l, m, n, p, b_min, b_max, id, numbin, ind, nspline, imin, imax
    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep, nmaps, ierr
    integer(i4b)       :: unit_summary
    logical(lgt)       :: exist
    real(dp),     allocatable, dimension(:,:) :: cls0, cls
    real(dp),     allocatable, dimension(:)   :: cl, lnL, lnL2
    real(dp),     allocatable, dimension(:)      :: cl_spline, lnL_spline
    character(len=1), dimension(6) :: status
    character(len=2), dimension(6) :: stext = ['TT','TE','TB','EE','EB','BB']
    real(dp), dimension(6) :: clmin = [0.d0,   -10.d0, -10.d0 , -1d0,  -1.d0,-1d0]
    real(dp), dimension(6) :: clmax = [10000.d0, 10.d0,  10.d0, 1.5d0,   1.d0, 1.5d0]

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
    numbin       = 10001
    nspline      = 100*numbin
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
    allocate(cl(numbin), lnL(numbin), lnL2(numbin))
    allocate(cl_spline(nspline), lnL_spline(nspline))
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

             imin = 1
             do while (lnL(imin) < 1d-6*maxval(lnL))
                write(*,*) imin, cl(imin), lnL(imin)
                imin = imin+1
             end do

             imax = numbin
             do while (lnL(imax) < 1d-6*maxval(lnL))
                imax = imax-1
             end do
     
             call spline(cl(imin:imax), lnL(imin:imax), 1.d30, 1.d30, lnL2(imin:imax))             
             do l = 1, nspline
                cl_spline(l) = cl(imin) + (cl(imax)-cl(imin))*(l-1.d0)/(nspline-1.d0)
                lnL_spline(l) = splint(cl(imin:imax), lnL(imin:imax), lnL2(imin:imax), cl_spline(l))
             end do
             

             call compute_conf_limits(cl_spline, lnL_spline, top, left, right)
             write(unit_summary,fmt='(f6.1,2i6,3f10.3,a,a)') 0.5d0*(lmin+lmax), lmin, lmax, top, left, right, '   ', spectext

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
    character(len=4)   :: lmin_text, lmax_text
    real(dp)           :: cl_min, cl_max, top, left, right, int, tot_int, cl_in(4), lnL_max, delta
    real(dp)           :: chisq, t1, t2, cl_maxlike, dC(1), alpha, alpha_prop, lnL_prop, lnL_old, xmin
    real(dp)           :: cl_prev, lnL_new(2), cl_new, eps, eps0, accept, reject, rms_prop(10000), clinit(10000)
    integer(i4b)       :: nspec, lmin, lmax, numchain, numiter, g, i, j, k, l, m, n, p, p_max(1), n_lnL, main_iter
    integer(i4b)       :: b_min, b_max, id, numbin, ind, ind1, ind2, n_spline, nsamp
    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep, nmaps, ierr
    integer(i4b)       :: unit, unit_out, iter, n_QML, maxiter
    logical(lgt)       :: exist, asymmetric_errors
    real(dp),     allocatable, dimension(:,:)    :: cls, L_prop
    real(dp),     allocatable, dimension(:,:,:)  :: sigmas
    real(dp),     allocatable, dimension(:)      :: cl, dCl, dCl0, lnL, alphas, binmask, W, cl0, lnL0, cl1, lnL1, lnL2
    real(dp),     allocatable, dimension(:,:)    :: invF0, invF, invF2, sigma, lnL_cond, samples
    real(dp),     allocatable, dimension(:)      :: cl_cond, p_cond, eta, mu
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
    call rand_init(handle, 34171)

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
          clinit(npar_powspec) = 0.
          select case (i)
          case (1)
             rms_prop(npar_powspec) = 50
             clinit(npar_powspec) = 0!1000
          case (2)
             rms_prop(npar_powspec) = 0.5
          case (3)
             rms_prop(npar_powspec) = 0.5
          case (4)
             rms_prop(npar_powspec) = 0.005
          case (5)
             rms_prop(npar_powspec) = 0.005
          case (6)
             rms_prop(npar_powspec) = 0.005
          end select
       end do
    end do
511 close(unit)

    ! Compute C_l with non-linear search; initialize on fiducial
    allocate(cl(npar_powspec), sigma(npar_powspec,2))
    do i = 1, npar_powspec
       cl(i) = mean(cls_fid(bins_powspec(i,1):bins_powspec(i,2),bins_powspec(i,3)))
    end do
    write(*,*) cl
    lnL_new(1) = lnL_powspec(cl)          
    write(*,*) 'lnL in =', lnL_new(1)
    if (lnL_new(1) == -1.d30) then
       write(*,*) '   Error: Initialization spectrum not acceptable'
       stop
    end if

    if (trim(optimizer) == 'quasinewton') then
!       call dfpmin(cl, 1.d-6, iter, chisq, chisq_powspec, dchisq_powspec, ierr)
    else if (trim(optimizer) == 'powell') then

       ! Cl_QML = Cl_old + invF * dC/dCl
       allocate(invF0(npar_powspec,npar_powspec))
       allocate(invF(npar_powspec,npar_powspec))
       invF0  = get_invfisher_powspec(cl)

       call powell(cl, neg_lnL_powspec, ierr, tolerance=1.d-3)
       !invF  = get_invfisher_powspec(cl)

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
       allocate(invF0(npar_powspec,npar_powspec))
       dCl   = dlnL_powspec(cl)
       invF  = get_invfisher_powspec(cl)
!       do j = 1, npar_powspec
!          write(*,*) j, cl(j), sqrt(invF(j,j))
!       end do
       invF0 = invF
       cl    = cl + matmul(invF, dCl)
       do j = 1, npar_powspec
          sigma(j,:) = sqrt(invF(j,j))
       end do

       lnL_prop = lnL_powspec(cl)
       write(*,*) '   Log-likelihood     = ', lnL_prop
       write(*,*) '   Reduced chisq      = ', chisq_powspec, sqrt(0.5d0/comm_lowl(1)%n)

    else if (trim(optimizer) == 'QML_line_search') then

       allocate(invF0(npar_powspec,npar_powspec))
       invF0  = get_invfisher_powspec(cl)

       n_QML = 50
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
          if (i < 3) cycle
          if (lnL(i)-lnL(i-2) < 0.1d0) converged = .true.

          if (all(converged)) exit
       end do

    else if (trim(optimizer) == 'conditional') then
    else if (trim(error_type) == 'marginal') then
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
    if (trim(error_type) == 'marginal') then
       !nsamp = 1000000
       !allocate(eta(npar_powspec), mu(npar_powspec), samples(npar_powspec,0:nsamp), L_prop(npar_powspec,npar_powspec))
       ! Marginal errors by MCMC
!       invF     = get_invfisher_powspec(cl)
       !call cholesky_decompose(invF, L_prop)
          allocate(eta(npar_powspec), mu(npar_powspec), L_prop(npar_powspec,npar_powspec), cl1(npar_powspec))
       do main_iter = 1, 3


          if (main_iter == 1) then
             nsamp = 10000
             allocate(samples(npar_powspec,0:nsamp))

             L_prop = 0.d0
             do i = 1, npar_powspec
                L_prop(i,i)  = rms_prop(i)
             end do

          else
             do i = 1, npar_powspec
                mu(i)  = mean(samples(i,nsamp/2+1:nsamp))
             end do
             L_prop = 0.d0
             do i = nsamp/2+1, nsamp
                eta = samples(:,i) -mu
                do j = 1, npar_powspec
                   do k = 1, npar_powspec
                      L_prop(j,k)  = L_prop(j,k) + eta(j)*eta(k)
                   end do
                end do
             end do
             L_prop = L_prop/(nsamp/2)
             call compute_hermitian_root(L_prop, 0.5d0)
             L_prop = L_prop * 0.2d0
             do j = 1, npar_powspec
                write(*,*) j, L_prop(j,j)
             end do
             
             if (main_iter == 2) then
                nsamp = 10000
             else
                nsamp = 10000
             end if
             deallocate(samples)
             allocate(samples(npar_powspec,0:nsamp))

             !stop
          end if
          samples(:,0) = clinit(1:npar_powspec) !0 !cl

       accept       = 0.d0
       lnL_old = lnL_powspec(samples(:,0))
       open(58,file='samples.dat', recl=10000)
       do i = 1, nsamp
          call wall_time(t1)
          do j = 1, npar_powspec
             eta(j) = rand_gauss(handle)
          end do
          cl1 = samples(:,i-1) + 1d0*matmul(L_prop,eta)
          lnL_prop = lnL_powspec(cl1)
          call wall_time(t2)
          if (exp(lnL_prop-lnL_old) > rand_uni(handle)) then
             accept       = accept+1.d0
             if (mod(i,100)== 0) write(*,fmt='(i7,a,2f16.3,2f8.3)') i, ' accept', lnL_prop, lnL_old, t2-t1, accept/i
             samples(:,i) = cl1
             lnL_old      = lnL_prop
          else
             if (mod(i,100) == 0) write(*,fmt='(i7,a,2f16.8,2f8.3)') i, ' reject', lnL_prop, lnL_old, t2-t1, accept/i
             samples(:,i) = samples(:,i-1)
          end if
          if (mod(i,10) == 0) write(58,*) i, real(lnL_old,sp), real(samples(:,i),sp)
       end do
       close(58)
       clinit(1:npar_powspec) = samples(:,nsamp)
       end do

       do j = 1, npar_powspec
          call int2string(bins_powspec(j,1), lmin_text)
          call int2string(bins_powspec(j,2), lmax_text)
          open(58,file=trim(outprefix)//'_'//stext(bins_powspec(j,3))//'_'//&
               & lmin_text // '_' // lmax_text // '.dat')
          write(58,*) '# Marginal distribution'
          do k = 1, nsamp
             write(58,*) k, real(samples(j,k),sp)
          end do
          close(58)

          if (trim(optimizer) == 'marginal') cl(j)      = mean(samples(j,nsamp/2+1:nsamp))
          sigma(j,:) = sqrt(variance(samples(j,nsamp/2+1:nsamp)))
          write(*,fmt='(a,i5,a,3f10.3)') 'bin = ', j, ' -- ', cl(j), sigma(j,:)
       end do

    else if (trim(error_type) == 'conditional') then

       invF     = get_invfisher_powspec(cl)       
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
          
          if (trim(optimizer) == 'QML_line_search') cl(j) = cl_max
          if (trim(optimizer) == 'conditional')     cl(j) = cl_max
          write(87,fmt='(i6,4f10.3)') j, cl(j), cl_max, sigma(j,:)
          write(*,fmt='(a,i5,a,3f10.3)') 'bin = ', j, ' -- ', cl(j), sigma(j,:)
       end do
       close(87)
       
    else if (trim(error_type) == 'fisher_fid') then
       
       do j = 1, npar_powspec
          sigma(j,:) = sqrt(invF0(j,j))
       end do

    else if (trim(error_type) == 'fisher') then
       
       invF = get_invfisher_powspec(cl)       
       do j = 1, npar_powspec
          sigma(j,:) = sqrt(invF(j,j))
       end do

    else if (trim(error_type) == 'noise') then

       invF0 = get_invfisher_powspec(0.d0*cl)       
       do j = 1, npar_powspec
          sigma(j,:) = sqrt(invF0(j,j))
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
!!$    do i = 1, npar_powspec
!!$       do j = i+1, npar_powspec
!!$          invF(i,j) = invF(i,j) / sqrt(invF(i,i)*invF(j,j))
!!$          invF(j,i) = invF(i,j)
!!$       end do
!!$    end do

!!$    open(unit_out,file=trim(outprefix) // '_corrmat.dat', recl=1024)
!!$    do i = 1, npar_powspec
!!$       write(unit_out,fmt='(4i6)',advance='no') i, bins_powspec(i,:)
!!$       do j = 1, i-1
!!$          write(unit_out,fmt='(f8.3)',advance='no') invF(i,j)
!!$       end do
!!$       write(unit_out,fmt='(f8.3)') 1.d0
!!$    end do
!!$    close(unit_out)

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

    open(unit_out,file=trim(outprefix) // '_corrmat.dat', recl=2048)
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
    !cls(:,[2,3,5]) = 0.d0

    neg_lnL_powspec = -comm_lowl_compute_lnL(cls=cls, ierr=ierr, enforce_pos_def=.false., &
         & red_chisq=chisq_powspec)
    if (ierr /= 0) neg_lnL_powspec = 1.d30

    !write(*,fmt='(e16.8,4e10.3)') neg_lnL_powspec, cls(2,[1,2,4,6])

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
    !call get_invfisher_safe(cls, comm_lowl(1)%w, transpose(comm_lowl(1)%P_harm), comm_lowl(1)%N_cov, &
    !     & bins_powspec(1:npar_powspec,:), invF)    
    !call get_invfisher_analytic(cls, comm_lowl(1)%w, 1.d0**2 * 4*pi/3072.d0, &
    !     & bins_powspec(1:npar_powspec,:), invF)    

    get_invfisher_powspec = invF

    ! 

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


  subroutine sigma2par_BR(partype)
    implicit none

    character(len=*), intent(in)   :: partype

    character(len=256) :: infile, binfile, outfile, temp, line, clfile, outprefix, templatefile
    character(len=2)   :: spectext
    character(len=4)   :: lmin_text, lmax_text
    real(dp)           :: cl_min, cl_max, sigma, top, left, right, int, tot_int, cl_in(4)
    real(dp)           :: A_min, A_max, A0, dA, mu, norm
    integer(i4b)       :: nspec, lmin, lmax, numchain, numiter, i, j, k, l, m, n, p, bins(4), id, numbin, ind
    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep, nmaps, ierr, spec_fit
    logical(lgt)       :: exist, include_BB
    real(dp), allocatable, dimension(:,:) :: cls0, cls
    real(dp), allocatable, dimension(:)   :: amp, lnL, F

    if (iargc() /= 7) then
       write(*,*) '    Compute spectrum with Blackwell-Rao'
       write(*,*) '    Options: [BR infofile] [templatefile] '
       write(*,*) '                 [default_amplitude] [spec_fit] [include_BB] [outprefix]'
       stop
    end if

    call getarg(2, infile)
    call getarg(3, templatefile)
    call getarg(4, temp)
    read(temp,*) A0
    call getarg(5, temp)
    read(temp,*) spec_fit
    call getarg(6, temp)
    read(temp,*) include_BB
    call getarg(7, outprefix)

    inquire(file=trim(infile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(infile)
    inquire(file=trim(templatefile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(templatefile)

    ! Read fiducial spectrum file
    unit = 58
    open(unit,file=trim(templatefile))
    lmax = -1
    do while (.true.)
       read(unit,'(a)',end=79) line
       line = trim(line)
       if (line(1:1) == '#') cycle
       read(line,*) l
       lmax = max(lmax, l)
    end do
79  close(unit)
    if (lmax > -1) then
       allocate(cls0(0:lmax,6), cls(0:lmax,6))
       cls0 = 0.d0
       open(unit,file=trim(templatefile))
       do while (.true.)
          read(unit,'(a)',end=80) line
          line = trim(line)
          if (line(1:1) == '#') cycle
          read(line,*) l, cl_in
          cls0(l,1) = cl_in(1) ! Assume (TT, EE, BB, TE) ordering for now
          cls0(l,2) = cl_in(4)          
          cls0(l,4) = cl_in(2)
          cls0(l,6) = cl_in(3)
       end do
80     close(unit)
    else
       write(*,*) 'Error -- comm_br_mod: No valid entries in clfile = ', trim(templatefile)
       stop
    end if
    if (.not. include_BB) cls0(:,6) = 0.d0

    ! Initialize BR module
    nmaps = 3
    call comm_br_initialize_object(infile)

    ! Loop over all free parameters, and output slices
    id    =  1
    A_min =  0.d0
    A_max =  10.d0
    numbin = 10000
    allocate(amp(numbin), lnL(numbin), F(numbin))
    cls = cls0
    do l = 1, numbin
       amp(l) = A_min + (A_max-A_min)/(numbin-1.d0) * (l-1)
       cls(:,spec_fit) = amp(l) * cls0(:,spec_fit)
       lnL(l) = comm_br_compute_lnL(cls, ierr, id)
       if (trim(partype) == 'tau') then
          amp(l) = sqrt(amp(l)) * A0
       else
          amp(l) = amp(l) * A0
       end if
    end do
    lnL = exp(lnL-maxval(lnL))
    norm = 0.d0
    do i = 1, numbin-1
       norm = norm + 0.5d0*(lnL(i)+lnL(i+1)) * (amp(i+1)-amp(i))
    end do
    lnL = lnL / norm

    outfile = trim(outprefix) // '.dat'
    write(*,*) 'Writing ', trim(outfile)
    open(unit,file=trim(outfile))
    do l = 1, numbin
       if (lnL(l) > 1.d-6 * maxval(lnL)) write(unit,*) amp(l), lnL(l)
    end do
    close(unit)

    mu    = 0.d0
    do i = 1, numbin-1
       mu = mu + 0.5d0 * (amp(i)*lnL(i)+amp(i+1)*lnL(i+1)) * (amp(i+1)-amp(i))
    end do
    sigma = 0.d0
    do i = 1, numbin-1
       sigma = sigma + 0.5d0*((amp(i)-mu)**2*lnL(i) + (amp(i+1)-mu)**2*lnL(i+1)) * (amp(i+1)-amp(i))
    end do
    sigma = sqrt(sigma)
    write(*,fmt='(a,f8.4,a,f8.4)') ' Estimated Gaussian parameter = ', real(mu,sp), ' +/-', real(sigma,sp)

    ! Compute cumulative distribution
    F(1) = 0.d0
    do i = 2, numbin
       F(i) = F(i-1) + 0.5d0 * (lnL(i)+lnL(i-1)) * (amp(i)-amp(i-1))
    end do
    F = F / F(numbin)
    
    do i = 1, numbin
       if (F(i) > 0.95d0) exit
    end do

    write(*,fmt='(a,f8.4)') ' Upper 95% confidence limit   = ', real(amp(i),sp)
    
  end subroutine sigma2par_BR

  subroutine sigma2qn_BR
    implicit none

    character(len=256) :: BRfile, outfile, temp, line, clfile, outprefix, templatefile, BRtype
    character(len=2)   :: spectext
    character(len=4)   :: lmin_text, lmax_text
    real(dp)           :: cl_min, cl_max, sigma, top, left, right, int, tot_int, cl_in(4)
    real(dp)           :: q_min, q_max, n_min, n_max, mu, l_pivot
    integer(i4b)       :: nspec, lmin, lmax, numchain, numiter, i, j, k, l, m, p, bins(4), id, numbin, ind
    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep, nmaps, ierr, spec_fit
    integer(i4b)       :: lmin_gauss, lmax_gauss, delta_l
    logical(lgt)       :: exist
    real(dp), allocatable, dimension(:,:) :: cls0, cls, lnL
    real(dp), allocatable, dimension(:)   :: q, n, lnL_marg

    if (iargc() < 12) then
       write(*,*) '    Compute 2D q-n contours with Blackwell-Rao'
       write(*,*) '    Options: [BRtype] [BR file] [templatefile] [q_min] [q_max] '
       write(*,*) '                [n_min] [n_max] [numbin] [l_pivot] [spec_fit] [outprefix]'
       write(*,*) '                [lmin (for gauss)] [lmax (for gauss)] [delta_l (for gauss)]'
       write(*,*) '       where BRtype = {exact, gauss}'
       stop
    end if

    call getarg(2, BRtype)
    call getarg(3, BRfile)
    call getarg(4, templatefile)
    call getarg(5, temp)
    read(temp,*) q_min
    call getarg(6, temp)
    read(temp,*) q_max
    call getarg(7, temp)
    read(temp,*) n_min
    call getarg(8, temp)
    read(temp,*) n_max
    call getarg(9, temp)
    read(temp,*) numbin
    call getarg(10, temp)
    read(temp,*) l_pivot
    call getarg(11, temp)
    read(temp,*) spec_fit
    call getarg(12, outprefix)
    if (trim(BRtype) == 'gauss') then
       call getarg(13, temp)
       read(temp,*) lmin_gauss
       call getarg(14, temp)
       read(temp,*) lmax_gauss
       call getarg(15, temp)
       read(temp,*) delta_l
    end if

    inquire(file=trim(BRfile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(BRfile)
    inquire(file=trim(templatefile), exist=exist)
    if (.not. exist) write(*,*) 'Error: file does not exist = ', trim(templatefile)

    call read_fiducial_spectrum(templatefile, cls0)
    lmax = size(cls0,1)-1
    allocate(cls(0:lmax,6))

    ! Initialize BR module
    if (trim(BRtype) == 'exact') then
       call comm_br_initialize_object(BRfile, handle=id)
    else if (trim(BRtype) == 'gauss') then
       call comm_gauss_br_initialize_object(BRfile, lmin_gauss, lmax_gauss, delta_l, handle=id)
    else if (trim(BRtype) == 'lowl') then
       call comm_lowl_initialize_object(BRfile)
    end if

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
          if (trim(BRtype) == 'exact') then
             lnL(i,j) = comm_br_compute_lnL(cls, ierr, handle=id)
          else if (trim(BRtype) == 'gauss') then
             lnL(i,j) = comm_gauss_br_compute_lnL(cls(2:lmax,1), handle=id)
          else if (trim(BRtype) == 'lowl') then
             lnL(i,j) = comm_lowl_compute_lnL(cls=cls)
             write(*,*) i, j, lnL(i,j)
          end if
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

  subroutine sigma2gauss
    implicit none

    character(len=256) :: sigmafile, gaussfile, temp
    real(dp)           :: cl
    integer(i4b)       :: nspec, lmin, lmax, numchain, numiter, i, j, k, l, m, n, lp, numbin
    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, lmax_chain, numchains
    integer(i4b)       :: numsamples, nbin, n_MC, seed, nsamp, delta_l
    real(dp)           :: cl_min, cl_max, rho
    logical(lgt)       :: exist, debug
    type(planck_rng)   :: handle
    real(dp), allocatable, dimension(:)       :: cls, F, mu_x, P_tot, mu_sigma
    real(dp), allocatable, dimension(:,:,:,:) :: sigma
    real(dp), allocatable, dimension(:,:)     :: sigma_1D, cov_x, P
    real(dp), allocatable, dimension(:,:)     :: samples
    real(dp), allocatable, dimension(:,:,:)   :: cl2x

    if (iargc() /= 10) then
       write(*,*) '    Compute Gaussianized Blackwell-Rao file from sigma file'
       write(*,*) '    Options: [sigmafile] [gaussfile] [lmin] [lmax]'
       write(*,*) '              [firstchain] [lastchain] [firstsample] [lastsample] '
       write(*,*) '              [n_MC for covar]'
       stop
    end if

    unit = comm_getlun()
    debug = .false.

    call getarg(2, sigmafile)
    call getarg(3, gaussfile)
    call getarg(4, temp)
    read(temp,*) lmin
    call getarg(5, temp)
    read(temp,*) lmax
    call getarg(6, temp)
    read(temp,*) firstchain
    call getarg(7, temp)
    read(temp,*) lastchain
    call getarg(8, temp)
    read(temp,*) firstsample
    call getarg(9, temp)
    read(temp,*) lastsample
    call getarg(10, temp)
    read(temp,*) n_MC

    inquire(file=trim(sigmafile), exist=exist)
    if (.not. exist) then
       write(*,*) 'Error: file does not exist = ', trim(sigmafile)
    end if
    inquire(file=trim(gaussfile), exist=exist)
    if (exist) then
       write(*,*) 'Error: Output file already exist = ', trim(gaussfile)
       stop
    end if

    call comm_read_chain(sigmafile, unit, lmax_chain, numchains, numsamples, sigma)

    ! Find number of valid samples
    nsamp = 0
    do j = firstchain, lastchain
       do k = firstsample, min(lastsample,numsamples)
          if (sigma(1,1,j,k) /= 0.d0) nsamp = nsamp+1
       end do
    end do
    allocate(sigma_1D(0:lmax_chain,nsamp))
    n = 0
    do j = firstchain, lastchain
       do k = firstsample, min(lastsample,numsamples)
          if (sigma(1,1,j,k) /= 0.d0) then
             n     = n+1
             sigma_1D(:,n) = sigma(:,1,j,k)
          end if
       end do
    end do

    ! Compute mean of sigma samples
    allocate(mu_sigma(lmin:lmax))
    do l = lmin, lmax
       mu_sigma(l) = mean(sigma_1D(l,1:n))
       write(*,*) l, real(sigma_1D(l,1:4),sp)
    end do

    ! Draw Monte Carlo samples
    seed = 12481
    call rand_init(handle, seed)
    allocate(samples(lmin:lmax,n_MC))
    do i = 1, n_MC
       if (mod(i,100000) == 0) write(*,*) 'Generating sample ', i, ' of ', n_MC
       do l = lmin, lmax
          rho = 0.d0
          do m = 1, 2*l-1
             rho = rho + rand_gauss(handle)**2
          end do
          samples(l,i) = (2*l+1)*sigma_1D(l,mod(i,nsamp)+1) / rho
       end do
    end do

    ! Compute cumulative distributions
    nbin = 1000
    allocate(cl2x(nbin,lmin:lmax,3), P(nbin,nsamp), P_tot(nbin), F(nbin))
    do l = lmin, lmax

       !write(*,*) 'Processing multipole ', l, ' of ', lmax

       if (debug) then
          open(58,file='samples.dat')
          do i = 1, n_MC
             write(58,*) i, samples(l,i)
          end do
          close(58)
       end if

       ! Set up C_l grid
       cl_min = minval(samples(l,:))
       cl_max = maxval(samples(l,:))
       cl2x(1,l,1) = cl_min
       do i = 2, nbin
          cl2x(i,l,1) = cl2x(i-1,l,1) * (cl_max/cl_min)**(1.d0/real(nbin-1,dp))
       end do

       ! Compute marginal distribution
       P = 0.d0
       do i = 1, nbin
          cl   = cl2x(i,l,1)
          do j = 1, nsamp
             P(i,j) = P(i,j) + 0.5d0 * ((2*l-1)*log(sigma_1D(l,j)) - (2*l+1)*log(cl) - &
                  (2*l+1) * sigma_1D(l,j)/cl)
          end do
       end do
       P = exp(P - maxval(P))
       do i = 1, nbin
          P_tot(i) = sum(P(i,:))
       end do

       if (debug) then
          open(58,file='P.dat')
          do i = 1, nbin
             write(58,*) cl2x(i,l,1), P_tot(i)
          end do
       end if

       ! Compute cumulative distribution
       F(1) = 0.d0
       do i = 2, nbin
          F(i) = F(i-1) + 0.5d0 * (P_tot(i)+P_tot(i-1)) * (cl2x(i,l,1)-cl2x(i-1,l,1))
       end do
       F = F / maxval(F)

       if (debug) then
          open(58,file='F.dat')
          do i = 1, nbin
             write(58,*) cl2x(i,l,1), F(i)
          end do
       end if

       ! Compute cl2x
       do i = 2, nbin-1
          if (F(i) < 0.5d0) then
             call convert_fract2sigma(cl2x(i,l,2), 1.d0-2.d0*F(i))
             cl2x(i,l,2) = -cl2x(i,l,2)
          else
             call convert_fract2sigma(cl2x(i,l,2), 2.d0*F(i)-1.d0)
          end if
       end do
       cl2x(1,l,2)    = cl2x(2,l,2)
       cl2x(nbin,l,2) = cl2x(nbin-1,l,2)

       if (debug) then
          open(58,file='cl2x.dat')
          do i = 1, nbin
             write(58,*) cl2x(i,l,1), cl2x(i,l,2)
          end do
       end if

       ! Spline conversion rule
       call spline(cl2x(:,l,1), cl2x(:,l,2), 0.d0, 0.d0, cl2x(:,l,3))

       ! Gaussianize samples
       do i = 1, n_MC
          samples(l,i) = splint(cl2x(:,l,1), cl2x(:,l,2), cl2x(:,l,3), samples(l,i))
       end do

       if (debug) then
          open(58,file='samples_x.dat')
          do i = 1, n_MC
             write(58,*) samples(l,i)
          end do
          close(58)
       end if
    end do

    ! Compute mean and covariance
    allocate(mu_x(lmin:lmax), cov_x(lmin:lmax,lmin:lmax))
    do l = lmin, lmax
       mu_x(l)    = sum(samples(l,:)) / n_MC
    end do

    do l = lmin, lmax
       do lp = lmin, lmax
          cov_x(l,lp) = sum((samples(l,:)-mu_x(l))*(samples(lp,:)-mu_x(lp)))/(n_MC-1)
       end do
    end do

    write(*,*) 'mu  = ', real(mu_x,sp)
    write(*,*) 'cov = ', real(cov_x,sp)

    ! Output compressed likelihood objects to file
    call write_gauss_BR_datafile(gaussfile, lmin, lmax, cl2x, mu_x, cov_x, mu_sigma)

  end subroutine sigma2gauss

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

  subroutine output_covmat_slice
    implicit none

    integer(i4b)       :: i, j, k, l, m, q, nside, nmaps, npix, lmax, lhigh_T, lhigh_P, ordering, numcomp, col
    integer(i4b)       :: l1, l2, m1, m2, i1, i2, n_d, n_p, n_h, n
    integer(i4b)       :: ind, ind1, ind2, cov_order, seed, numtemp, nmask, n_t, n_pol
    integer(i4b)       :: lmax_basis, pol, nval, nhigh, num_fg_temp, nside_temp, npix_temp
    real(dp)           :: t1, t2, Tthreshold, Pthreshold, reg(3), scale, W_max_T, W_max_P
    real(dp)           :: BBamp, cvec(3), weight_T, weight_P
    logical(lgt)       :: comp_cov
    character(len=5)   :: itext
    character(len=512) :: mapfile, maskfile, covfile, beamfile, outfile, filename, temp
    character(len=512) :: basis, prefix, clfile, partext, tempcovmat(10)
    real(dp),     allocatable, dimension(:,:)     :: cov, P_harm, C, cl, B, invC, tempcov
    real(dp),     allocatable, dimension(:)       :: W, mask_1d, test3, A_T, map_1d
    real(dp),     allocatable, dimension(:,:)     :: d, mask, map, alms, V, P_l, V_red, S_mat, P
    real(dp),     allocatable, dimension(:,:)     :: Y, Yp, beam, cls, T, buffer, buffer2, P_highl, weights
    real(dp),     allocatable, dimension(:,:)     :: PT, PT_invN_TP
    complex(dpc), allocatable, dimension(:,:,:)   :: alms_cmplx
    real(dp),     pointer,     dimension(:,:)     :: pixwin
    real(dp),     allocatable, dimension(:,:,:)   :: window
    integer(i4b), allocatable, dimension(:)       :: indmap, map2mask, lmax2lhigh
    logical(lgt), allocatable, dimension(:)       :: tempmode, rescale_CMB
    type(planck_rng) :: handle


    if (iargc() < 5) then
       write(*,*) '    Pre-process low-l likelihood inputs from map and covariance matrix'
       write(*,*) '    Options:  [mapfile] [maskfile] [covfile] [col]'
       stop
    end if

    unit = 58
    call getarg(2,mapfile)
    call getarg(3,maskfile)
    call getarg(4,covfile)
    call getarg(5,temp)	
    read(temp,*) col

    write(*,*) 'Reading data:'
    write(*,*) '     Map file   = ', trim(mapfile)
    write(*,*) '     Mask       = ', trim(maskfile)
    write(*,*) '     Covariance = ', trim(covfile)

    i    = getsize_fits(trim(mapfile), ordering=ordering, nside=nside, nmaps=nmaps)
    npix = 12*nside**2
    lmax = 2*nside
    n_d  = 1
    n_p  = npix    * nmaps
    allocate(map(0:npix-1,nmaps), d(0:npix-1,nmaps))
    call read_map(mapfile,  map)

    allocate(mask(0:npix-1,nmaps), mask_1d(n_p))
    call read_map(maskfile, mask)
    where (mask < 0.5d0)
       mask = 0.d0
    elsewhere
       mask = 1.d0
    end where

    do i = 1, nmaps
       mask_1d((i-1)*npix+1:i*npix) = mask(:,i)
    end do

    ! Set up map2mask conversion array
    nval = count(mask > 0.5d0)
    allocate(map2mask(nval), map_1d(nval))
    k = 1
    do j = 1, nmaps
       do i = 0, npix-1
          if (mask(i,j) > 0.5d0) then
             map2mask(k) = (j-1)*npix+i+1
             map_1d(k)   = map(i,j)
             k           = k+1
          end if
       end do
    end do

    ! Read covariance matrix
    allocate(cov(nval,nval))
    call comm_read_covmat(covfile, nmaps, cov, scale)
    !cov = cov/2

    map_1d = map_1d * scale ! Assume map is in same units as covariance matrix

    allocate(W(nval))
    call get_eigenvalues(cov,W)
    write(*,*) maxval(W)/minval(W)
    deallocate(W)

    d = 0.d0
    k = 1
    do j = 1, nmaps
       do i = 0, npix-1
          if (mask(i,j) > 0.5d0) then
             d(i,j) = cov(k,col)
             if (k == col) d(i,j) = -1.6375d30 
             k           = k+1
          end if
       end do
    end do
    
    call write_map3('slice.fits', d)

    d = 0.d0
    k = 1
    do j = 1, nmaps
       do i = 0, npix-1
          if (mask(i,j) > 0.5d0) then
             d(i,j) = sqrt(cov(k,k))
             k           = k+1
          end if
       end do
    end do
    
    call write_map3('rms.fits', d)

    k = count(mask(:,1) > 0.5d0)
   do i = 1, k
      cov(i,k+1:) = 0.d0
      cov(k+1:,i) = 0.d0
   end do
    


    write(*,*) minval(cov), maxval(cov)
    call compute_hermitian_root(cov, -0.5d0)
    write(*,*) minval(cov), maxval(cov)
    map_1d = matmul(cov, map_1d)
    write(*,*) 'chisq =', (sum(map_1d**2)- nval)/sqrt(2.d0*nval)

    d = 0.d0
    k = 1
    do j = 1, nmaps
       do i = 0, npix-1
          if (mask(i,j) > 0.5d0) then
             d(i,j) = map_1d(k)
             k           = k+1
          end if
       end do
    end do
    
    call write_map3('whitened.fits', d)


    if (nmaps == 1) then
       write(*,fmt='(a,f8.3)') '      I sky fractions = ', sum(mask(:,1))/npix
    else
       write(*,fmt='(a,3f8.3)') '      I/Q/U sky fractions = ', sum(mask(:,1))/npix, &
            & sum(mask(:,2))/npix, sum(mask(:,3))/npix
    end if


    allocate(alms_cmplx(3,0:lmax,0:lmax))

    call compute_hermitian_root(cov, -1d0)
    



  end subroutine output_covmat_slice


end program comm_like_tools
