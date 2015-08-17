program comm_like_tools
  use comm_proc_utils
  use comm_br_mod
  use comm_lowl_mod
  implicit none

  character(len=128) :: operation
  integer(i4b)       :: unit

  if (iargc() == 0) then
     write(*,*) 'Usage: comm_like_tools [operation] [args] ... '
     write(*,*) '  Valid operations are:'
     write(*,*) '     sigma2cl_BR        -- '
     write(*,*) '     sigma2slice_BR     -- '
     write(*,*) '     sigma2slice_gauss  -- '
     write(*,*) '     sigma2par_BR       -- '
     write(*,*) '     sigma2tau_BR       -- '
     write(*,*) '     sigma2qn_BR        -- '
     write(*,*) '  For further usage information, type "comm_like_tools [operation]"'
     stop
  end if

  unit = 19
  call getarg(1,operation)
  if (trim(operation) == 'mapcov2gausslike') then
     call mapcov2gausslike
!!$  else if (trim(operation) == 'sigma2cl_BR') then
!!$     call sigma2cl_BR
  else if (trim(operation) == 'sigma2slice_BR') then
     call sigma2slice_BR
  else if (trim(operation) == 'sigma2slice_gauss') then
     call sigma2slice_gauss
!!$  else if (trim(operation) == 'sigma2par_BR') then
!!$     call sigma2par_BR('par')
!!$  else if (trim(operation) == 'sigma2tau_BR') then
!!$     call sigma2par_BR('tau')
  else if (trim(operation) == 'sigma2qn_BR') then
     call sigma2qn_BR
  end if

contains

  subroutine mapcov2gausslike
    implicit none

    integer(i4b)       :: i, j, k, l, m, nside, nmaps, npix, lmax, ordering, numcomp, n_p, n_h, n
    integer(i4b)       :: ind, ind1, ind2, cov_order, seed, numtemp
    logical(lgt)       :: rem_monodipole
    real(dp)           :: t1, t2, threshold, reg(3)
    character(len=5)   :: itext
    character(len=512) :: mapfile, maskfile, covfile, beamfile, outfile, filename, temp, eig
    character(len=512) :: basis, prefix, clfile
    real(dp),     allocatable, dimension(:,:)     :: cov, P_pix, P_harm, P_full, C, C_full, cl, B, CB
    real(dp),     allocatable, dimension(:)       :: W, d, d_full, mask_1d, map_1d, integral, buffer
    real(dp),     allocatable, dimension(:,:)     :: mask, map, alms, V
    real(dp),     allocatable, dimension(:,:)     :: Y, beam, cls, T
    complex(dpc), allocatable, dimension(:,:,:)   :: alms_cmplx
    real(dp),     pointer,     dimension(:,:)     :: pixwin
    real(dp),     allocatable, dimension(:,:,:)   :: window
    integer(i4b), allocatable, dimension(:)       :: indmap
    type(planck_rng) :: handle


    if (iargc() /= 14) then
       write(*,*) '    Pre-process low-l likelihood inputs from map and covariance matrix'
       write(*,*) '    Options:  [mapfile] [maskfile] [covfile] [beamfile] [clfile] [lmax] '
       write(*,*) '                 [basis] [threshold] [project mono/dipole] '
       write(*,*) '                 [T reg] [P reg] [seed] [outprefix]'
       write(*,*) '        where [basis]     = {pixel, pseudoalm, eigen_N, eigen_S+N}'
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
    call getarg(8,basis)
    call getarg(9,temp)
    read(temp,*) threshold
    call getarg(10,temp)
    read(temp,*) rem_monodipole
    call getarg(11,temp)
    read(temp,*) reg(1)
    call getarg(12,temp)
    read(temp,*) reg(2)
    reg(3) = reg(2)
    call getarg(13,temp)
    read(temp,*) seed
    call getarg(14,prefix)
    numcomp = (lmax+1)**2

    call rand_init(handle, seed)

    write(*,*) 'Reading data:'
    write(*,*) '     Map        = ', trim(mapfile)
    write(*,*) '     Mask       = ', trim(maskfile)
    write(*,*) '     Covariance = ', trim(covfile)
    write(*,*) '     Cl file    = ', trim(clfile)

    ! Read fiducial spectrum; only actually used for S+N or S/N bases
    call read_fiducial_spectrum(clfile, cls)

    ! Read map and mask
    i    = getsize_fits(trim(mapfile), ordering=ordering, nside=nside, nmaps=nmaps)
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

    ! Read templates to remove
    numtemp = iargc()-14
    if (rem_monodipole) numtemp = numtemp+4
    if (numtemp > 0) then
       allocate(T(n_p,numtemp))
       T = 0.d0
    end if
    do k = 1, iargc()-14
       call getarg(14+k,filename)
       call read_map(mapfile, map)
       do i = 1, nmaps
          T((i-1)*npix+1:i*npix,k) = map(:,i)
       end do
    end do
    if (rem_monodipole) then
       j = iargc()-14
       T(1:npix,j+1) = 1.d0 ! Monopole
       do i = 0, npix-1
          call pix2vec_ring(nside, i, T(i+1,j+2:j+4)) ! Dipoles
       end do
    end if

    ! Read beam and pixel window
    allocate(beam(0:lmax,nmaps), window(0:lmax,nmaps,nmaps))
    call read_pixwin(nside, nmaps, pixwin)
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
    if (.true.) then
       open(unit,file=trim(covfile),form='unformatted')
       read(unit) i
       read(unit) cov_order
       read(unit) i
       do i = 1, n_p
          read(unit) cov(:,i)
       end do
       close(unit)
       if (cov_order == 2) then
          allocate(buffer(0:npix-1))
          ! Convert columns from nest to ring
          do i = 1, n_p
             do j = 1, nmaps
                buffer = cov((j-1)*npix+1:j*npix,i)
                call convert_nest2ring(nside, buffer)
                cov((j-1)*npix+1:j*npix,i) = buffer
             end do
          end do
          ! Convert rows from nest to ring
          do i = 1, n_p
             do j = 1, nmaps
                buffer = cov(i,(j-1)*npix+1:j*npix)
                call convert_nest2ring(nside, buffer)
                cov(i,(j-1)*npix+1:j*npix) = buffer
             end do
          end do
          deallocate(buffer)
       end if
    else
       cov = 0.d0
       do i = 1, n_p
          cov(i,i) = 2.d0**2
       end do
    end if

    ! Add regularization noise to covariance matrix
    ind = 1
    do j = 1, nmaps
       do i = 1, npix
          cov(ind,ind) = cov(ind,ind) + reg(j)**2
          ind          = ind+1
       end do
    end do

    if (nmaps == 1) then
       write(*,fmt='(a,f8.3)') '      I sky fractions = ', sum(mask(:,1))/npix
    else
       write(*,fmt='(a,3f8.3)') '      I/Q/U sky fractions = ', sum(mask(:,1))/npix, &
            & sum(mask(:,2))/npix, sum(mask(:,3))/npix
    end if

    ! Apply mask 
    map_1d = map_1d * mask_1d
    do i = 1, n_p
       if (mask_1d(i) == 0.d0) then
          cov(i,:) = 0.d0
          cov(:,i) = 0.d0
       end if
    end do
    do i = 1, numtemp
       T(:,i) = T(:,i) * mask_1d
    end do

    ! Fit and subtract templates from data, and add to covariance matrix
    if (numtemp > 0) then

    end if

    ! Compute spherical harmonics
    write(*,*) 'Computing spherical harmonics'
    n = n_h
    allocate(Y(n_p, n_h), alms(numcomp,nmaps), alms_cmplx(nmaps,0:lmax,0:lmax))
    allocate(P_pix(n_p,n))
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
                Y((i-1)*npix+1:i*npix,ind1) = mask(:,i)*map(:,i)
             end do
             ind1 = ind1+1
             ind2 = ind2+1
          end do
       end do
    end do

    if (trim(basis) == 'pixel') then

       n = count(mask > 0.5d0)
       allocate(indmap(n), d(n), C(n,n), P_harm(n_h,n))
       k   = 0
       ind = 1
       do j = 1, nmaps
          do i = 0, npix-1
             if (mask(i,j) > 0.5d0) then
                k = k+1
                indmap(k) = ind
             end if
             ind = ind+1
          end do
       end do
       d = map_1d(indmap)
       C = cov(indmap,indmap)
       P_harm = transpose(Y(indmap,:))
       deallocate(indmap)

       !d = matmul(transpose(P_pix), map_1d) !/ (npix/(4.d0*pi))
       !C = matmul(transpose(P_pix), matmul(cov, P_pix)) !/ (npix/(4.d0*pi))**2
       !call wall_time(t2)
       !write(*,*) sum(abs(d)), sum(abs(C)), t2-t1

    else 

       n = n_h
       allocate(d_full(n), C_full(n,n), P_full(n_h,n_h))
       d_full = matmul(transpose(Y), map_1d) * (4.d0*pi/npix)
       C_full = matmul(transpose(Y), matmul(cov,Y)) * (4.d0*pi/npix)**2
       !call dgemm('N','N',n_p,n,n_p,1.d0,cov,n_p,B,n_p,0.d0,CB,n_p)
       !call dgemm('T','N',n,n,n_p,1.d0,B,n_p,CB,n_p,0.d0,C,n)
       call dgemm('T','N',n_h,n,n_p,1.d0,Y,n_p,Y,n_p,0.d0,P_full,n_h)
       P_full = P_full * 4.d0*pi/npix

       ! Set up basis vectors
       write(*,*) 'Computing basis vectors for ', trim(basis)
       allocate(B(n_h,n_h), W(n_h), V(n_h,n_h))
       if (trim(basis) == 'pseudoalm') then
          B = P_full
       else if (trim(basis) == 'eigen_N') then
          B = C_full
       else if (trim(basis) == 'eigen_S+N') then
          write(*,*) 'Not yet supported = ', trim(basis)
       else
          write(*,*) 'Unsupported basis vector type = ', trim(basis)
          stop
       end if
       call get_eigen_decomposition(B, W, V)

       ! Print eigen spectrum
       outfile = trim(prefix) // '_W.dat'
       write(*,*) '   Writing basis eigen spectrum to ', trim(outfile)
       open(unit,file=trim(outfile))
       do i = 1, n_h
          write(unit,*) i, abs(W(i)/maxval(W)), W(i)
       end do
       close(unit)

       ! Transform data, covariance and projection operator to new basis
       n = count(W/maxval(W) > threshold)
       allocate(d(n), C(n,n), P_harm(n_h,n))
       d = matmul(transpose(V(:,n_h-n+1:n_h)), d_full)
       C = matmul(transpose(V(:,n_h-n+1:n_h)), matmul(C_full, V(:,n_h-n+1:n_h)))
       P_harm = matmul(P_full, V(:,n_h-n+1:n_h))
!!$      ! call dgemm('N','N',n_h,n,n_h,1.d0,B,n_h,V(:,n_h-n+1:n_h),n_h,0.d0,P_harm,n_h)
!!$       call dgemm('N','N',n_h,n,n_h,1.d0,B,n_h,V(:,n_h-n+1:n_h),n_h,0.d0,P_harm,n_h)
!!$
!!$       ! Transform data and covariance to new basis
!!$       call wall_time(t1)
!!$       deallocate(B); allocate(B(n_p,n), CB(n_p,n))
!!$       !B = matmul(Y, P_harm)
!!$       call dgemm('N','N',n_p,n,n_h,1.d0,Y,n_p,P_harm,n_h,0.d0,B,n_p)
!!$       d = matmul(transpose(B), map_1d)
!!$       call dgemm('N','N',n_p,n,n_p,1.d0,cov,n_p,B,n_p,0.d0,CB,n_p)
!!$       call dgemm('T','N',n,n,n_p,1.d0,B,n_p,CB,n_p,0.d0,C,n)
!!$!       d = d * (4.d0*pi/npix)
!!$!       C = C * (4.d0*pi/npix)**2
!!$       !C = matmul(transpose(B), matmul(cov, B)) !/ (npix/(4.d0*pi))**2
!!$       deallocate(B, CB)
!!$       call wall_time(t2)
!!$!       write(*,*) sum(abs(d)), sum(abs(C)), t2-t1
!       stop

    end if

    write(*,*) '   Number of pixels        = ', n_p
    write(*,*) '   Number of harmonics     = ', n_h
    write(*,*) '   Number of basis vectors = ', n

    ! Output data objects
    outfile = trim(prefix) // '.unf'
    write(*,*) 'Writing data objects to ', trim(outfile)
    open(unit,file=trim(outfile),form='unformatted')
    write(unit) lmax, nmaps
    write(unit) n_h, n
    write(unit) d
    write(unit) C
    write(unit) window
    write(unit) P_harm
    close(unit)

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


  subroutine sigma2slice_gauss
    implicit none

    character(len=256) :: infile, binfile, outfile, temp, line, clfile, outprefix
    character(len=2)   :: spectext, ptext
    character(len=4)   :: lmin_text, lmax_text
    real(dp)           :: cl_min, cl_max, sigma, top, left, right, int, tot_int, cl_in(4), lnL_max
    integer(i4b)       :: nspec, lmin, lmax, numchain, numiter, g, i, j, k, l, m, n, p, b_min, b_max, id, numbin, ind
    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep, nmaps, ierr
    logical(lgt)       :: exist
    real(dp), allocatable, dimension(:,:) :: cls0, cls
    real(dp), allocatable, dimension(:)   :: cl, lnL
    character(len=1), dimension(6) :: status
    character(len=2), dimension(6) :: stext = ['TT','TE','TB','EE','EB','BB']
    real(dp), dimension(6) :: clmin = [0.d0,   -10.d0, -1.d0, 0.d0,  -1.d0, 0.d0]
    real(dp), dimension(6) :: clmax = [5000.d0, 10.d0,  1.d0, 0.2d0,  1.d0, 0.01d0]

    if (iargc() /= 5) then
       write(*,*) '    Compute Gaussian likelihood slices'
       write(*,*) '      Usage: comm_like_tools sigma2slice_gauss [BR file] [clfile] [binfile] [outprefix]'
       stop
    end if

    call getarg(2, infile)
    call getarg(3, clfile)
    call getarg(4, binfile)
    call getarg(5, outprefix)
    nspec  = 6
    numbin = 100

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
    do while (.true.)
       read(58,'(A)',end=518) line
       if (line(1:1) == '#' .or. trim(line)=='') cycle
       read(line,*) lmin, lmax, status
       do ind = 1, nspec
          if (ind /= 4) cycle
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
                lnL(l) = comm_lowl_compute_lnL(cls, ierr)
                write(*,*) cl(l), lnL(l)
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
          end if
          
       end do
    end do
518 close(58)

  end subroutine sigma2slice_gauss


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


end program comm_like_tools
