program comm_process_resfiles
  use comm_proc_utils
  use comm_br_mod
  use comm_lowl_mod
  implicit none

  character(len=128) :: operation
  integer(i4b)       :: unit

  if (iargc() == 0) then
     write(*,*) 'Usage: comm_process_resfiles [operation] [args] ... '
     write(*,*) '  Valid operations are:'
     write(*,*) '     cl2fits            -- '
     write(*,*) '     sigma2fits         -- '
     write(*,*) '     map2mean_cov       -- '
     write(*,*) '     pix2mean_cov       -- '
     write(*,*) '     alm2mean_cov       -- '
     write(*,*) '     free_temp2fits     -- '
     write(*,*) '     ind_temp2fits      -- '
     write(*,*) '     fg_amp2fits        -- '
     write(*,*) '     fg_ind2fits        -- '
     write(*,*) '     cl2single_chain    -- '
     write(*,*) '     sigma2single_chain -- '
     write(*,*) '     mcmc2single_chain  -- '
     write(*,*) '     mean_rms2fits      -- '
     write(*,*) '     fg_amp2fits        -- '
     write(*,*) '     fg_ind2fits        -- '
     write(*,*) '     cmb2fits           -- '
     write(*,*) '     IQU_mean2fits      -- '
     write(*,*) '     TEB_mean2fits      -- '
     write(*,*) '     mean_stddev2fits   -- '
     write(*,*) '     trim               -- '
     write(*,*) '     flatten            -- '
     write(*,*) '     cl2fits_filelist   -- '
     write(*,*) '     cl2ascii_binned    -- '
     write(*,*) '     sigma2gelman_rubin -- '
     write(*,*) '  For further usage information, type "comm_process_resfiles [operation]"'
     stop
  end if

  unit = 19
  call getarg(1,operation)
  if (trim(operation) == 'cl2fits') then
     call cl2fits('cls')
  else if (trim(operation) == 'sigma2fits') then
     call cl2fits('sigma')
  else if (trim(operation) == 'free_temp2fits') then
     call fg_amp2fits('free')
  else if (trim(operation) == 'ind_temp2fits') then
     call fg_amp2fits('ind')
  else if (trim(operation) == 'fg_amp2fits') then
     call map2fits('amp')
  else if (trim(operation) == 'fg_ind2fits') then
     call map2fits('ind')
  else if (trim(operation) == 'cmb2fits') then
     call cmb2fits
  else if (trim(operation) == 'mean_rms2fits') then
     call mean_rms2fits
  else if (trim(operation) == 'IQU_mean2fits') then
     call cmb_mean2fits('IQU')
  else if (trim(operation) == 'TEB_mean2fits') then
     call cmb_mean2fits('TEB')
  else if (trim(operation) == 'cl2single_chain') then
     call cl2fits_singlechain('cls')
  else if (trim(operation) == 'sigma2single_chain') then
     call cl2fits_singlechain('sigma')
  else if (trim(operation) == 'mcmc2single_chain') then
     call cl2fits_singlechain('mcmc')
  else if (trim(operation) == 'trim') then
     call trim_file
  else if (trim(operation) == 'map2mean_cov') then
     call maps2mean_cov
  else if (trim(operation) == 'pix2mean_cov') then
     call chain2mean_cov('pix')
  else if (trim(operation) == 'alm2mean_cov') then
     call chain2mean_cov('alm')
  else if (trim(operation) == 'flatten') then
     call flatten_chain
  else if (trim(operation) == 'cl2fits_filelist') then
     call cl2fits_from_files
  else if (trim(operation) == 'cl2ascii_binned') then
     call cl2ascii_binned
  else if (trim(operation) == 'ASCII2fits') then
     call ASCII2fits
  else if (trim(operation) == 'sigma2gelman_rubin') then
     call gelman_rubin
  else if (trim(operation) == 'coadd_mapcov') then
     call coadd_mapcov
  end if

contains


  ! ************************************************************************
  !                        Operational routines
  ! ************************************************************************

  subroutine gelman_rubin
    implicit none

    character(len=256) :: infile, outfile, temp
    integer(i4b)       :: i, j, k, l, m, n
    integer(i4b)       :: burnin, numchain, numiter, nspec, lmax, thinstep
    real(dp)           :: W, mu, B, V, R
    real(dp), allocatable, dimension(:)       :: mu_chain
    real(sp), pointer,     dimension(:,:,:,:) :: cls

    if (iargc() /= 4) then
       write(*,*) '    Compute Gelman-Rubin statistic from sigma/cls FITS file'
       write(*,*) '    Options: [infile] [burnin] [outfile] '
       stop
    end if

    call getarg(2,infile)
    call getarg(3,temp)
    read(temp,*) burnin
    call getarg(4,outfile)

    unit = 58
    call read_resfile(infile, unit, nspec, lmax, numchain, numiter, cls)

    allocate(mu_chain(numchain))
    open(unit,file=trim(outfile))
    do j = 1, nspec

       n = 100000000
       do k = 1, numchain
          n = min(n, count(cls(2,j,k,burnin+1:0+numiter) /= 0.d0))
       end do
       write(*,*) 'Spectrum = ', j, ', n = ', n
       if (n == 0) cycle
       write(unit,*) '#  Spectrum = ', j, ', nchain = ', n, ', ntot = ', n*numchain

       do l = 2, lmax
          
          ! Compute mean for each chain
          mu_chain = 0.d0
          W  = 0.d0
          do k = 1, numchain
             mu_chain(k) = sum(cls(l,j,k,burnin+1:burnin+0+n)) / n
             W = W + sum((cls(l,j,k,burnin+1:burnin+0+n)-mu_chain(k))**2)
          end do
          mu = sum(mu_chain)/numchain
          B  = n / (numchain-1.d0) * sum((mu_chain-mu)**2)
          W  = W / numchain / (n-1.d0)
          V  = (n-1.d0)/n * W + 1.d0/n * B
          R  = sqrt(V/W * n/(n-1.d0))
          write(unit,*) l, R
       end do
       write(unit,*) 
    end do
    close(unit)

  end subroutine gelman_rubin

  subroutine ASCII2fits
    implicit none

    character(len=256) :: outfile, temp, infile
    character(len=4)   :: chain_text
    integer(i4b)       :: i, j, k, l, stat, unit
    integer(i4b)       :: numiter, numchain, lmax, nspec, numburnin, max_num_samples, order(6)
    logical(lgt)       :: reached_limit, exist
    real(dp)           :: val

    real(sp),     allocatable, dimension(:,:,:,:)  :: pre_cls, cls

    if (iargc() /= 5) then
       write(*,*) '    Create power spectrum FITS file'
       write(*,*) '    Options:  [infile] [lmax] [nspec] [outfile]'
       stop
    end if

    unit = 58

    call getarg(2,infile)
    call getarg(3,temp)
    read(temp,*) lmax
    call getarg(4,temp)
    read(temp,*) nspec
    call getarg(5,outfile)

    numchain = 1
    numiter  = 1
    allocate(cls(0:lmax, nspec, numchain, 0:numiter))

    order        = [1, 4, 6, 2, 3, 5]
    cls          = 0.d0
    cls(:,:,:,0) = 1.d0
    open(unit,file=trim(infile))
    do j = 1, nspec
       do i = 0, lmax
          read(unit,*) l, val
          cls(l,order(j),1,1) = val*1d12
       end do
    end do
    close(unit)

    call write_resfile(outfile, lmax, nspec, 1, numchain, "Gibbs sampled power spectra", cls)

    deallocate(cls)

  end subroutine ASCII2fits

  subroutine cl2fits(filetype)
    implicit none

    character(len=*), intent(in) :: filetype

    character(len=256) :: outfile, temp, infile
    character(len=4)   :: chain_text
    integer(i4b)       :: i, j, k, l, stat
    integer(i4b)       :: numiter, numchain, lmax, nspec, numburnin, max_num_samples
    logical(lgt)       :: reached_limit, exist

    real(sp),     allocatable, dimension(:,:,:,:)  :: pre_cls, cls

    if (iargc() /= 7) then
       write(*,*) '    Create power spectrum FITS file'
       write(*,*) '    Options:  [outfile] [max number of iterations per chain] '
       write(*,*) '                  [number of chains] [lmax] [nspec] [burnin]'
       stop
    end if

    call getarg(2,outfile)
    call getarg(3,temp)
    read(temp,*) numiter
    call getarg(4,temp)
    read(temp,*) numchain
    call getarg(5,temp)
    read(temp,*) lmax
    call getarg(6,temp)
    read(temp,*) nspec
    call getarg(7,temp)
    read(temp,*) numburnin

    allocate(pre_cls(0:lmax, nspec, numchain, 0:numiter))

    pre_cls = 0.d0

    max_num_samples = 0
    do j = 1, numchain

       call int2string(j, chain_text)

       if (trim(filetype) == 'cls') then
          infile = 'chain_cls_no' // chain_text // '.unf'
       else if (trim(filetype) == 'sigma') then
          infile = 'chain_sigma_no' // chain_text // '.unf'
       else if (trim(filetype) == 'mcmc') then
          infile = 'chain_mcmc_no' // chain_text // '.unf'
       else
          write(*,*) 'Unknown file type.'
          stop
       end if

       inquire(file=trim(infile), exist=exist)
       open(unit, file=trim(infile), form='unformatted')
       read(unit) k
       l = 1
       reached_limit = .false.
       do i = 1, numiter
          read(unit, iostat=stat) pre_cls(:,:,j,l)
          if (stat < 0) then
             pre_cls(:, :, j, 0) = l - 1
             reached_limit = .true.
             exit
          end if
          if (i > numburnin) l = l + 1
       end do

       if (.not. reached_limit) pre_cls(:, :, j, 0) = numiter - numburnin

       close(unit)

    end do

    max_num_samples = numiter
    do while (all(pre_cls(:,:,:,max_num_samples) == 0.d0))
       max_num_samples = max_num_samples-1
    end do

    allocate(cls(0:lmax, nspec, numchain, 0:max_num_samples))
    cls = pre_cls(:, :, :, 0:max_num_samples)

    call write_resfile(outfile, lmax, nspec, max_num_samples, numchain, &
         & "Gibbs sampled power spectra", cls)

    deallocate(cls, pre_cls)

  end subroutine cl2fits

  subroutine fg_amp2fits(filetype)
    implicit none

    character(len=*), intent(in) :: filetype

    character(len=256) :: outfile, temp, infile
    character(len=4)   :: chain_text
    integer(i4b)       :: numtemp, numband, numiter, numchain, i, j, k, l
    real(sp),     allocatable, dimension(:,:,:,:)  :: cls

    if (iargc() /= 6) then
       write(*,*) '    Create foreground template coefficient FITS file'
       write(*,*) '    Options:  [outfile] [number of templates] [number of channels]'
       write(*,*) '              [number of iterations] [number of chains]'
       stop
    end if

    call getarg(2,outfile)
    call getarg(3,temp)
    read(temp,*) numtemp
    call getarg(4,temp)
    read(temp,*) numband
    call getarg(5,temp)
    read(temp,*) numiter
    call getarg(6,temp)
    read(temp,*) numchain

    allocate(cls(numtemp, numband, numchain, numiter))
    do j = 1, numchain

       call int2string(j, chain_text)

       if (trim(filetype) == 'free') then
          infile = 'chain_free_fgtemp_no' // chain_text // '.unf'
       else if (trim(filetype) == 'ind') then
          infile = 'chain_ind_fgtemp_no' // chain_text // '.unf'
       end if
       open(unit, file=trim(infile), form='unformatted')
       read(unit) k
       do i = 1, numiter
          read(unit) cls(:,:,j,i)
       end do
       close(unit)

    end do

    call write_fgt_resfile(outfile, cls)

    deallocate(cls)

  end subroutine fg_amp2fits

  subroutine chain2mean_cov(filetype)
    implicit none

    character(len=*), intent(in) :: filetype

    character(len=256) :: outfile, temp, infile, prefix, beamfile, maskfile
    character(len=4)   :: chain_text
    character(len=5)   :: pix_text
    character(len=1)   :: s_text(3) = ['T', 'Q', 'U']
    integer(i4b)       :: i, j, k, l, m, n, p, unit, nreal, ordering, chain, lmax, counter, nval
    integer(i4b)       :: nside, npix, nmaps, comp, ncomp, burnin, numiter, numchain, seed, numcomp, nmax
    real(dp)           :: sigma_reg(3)
    logical(lgt)       :: compute_covar
    type(planck_rng)   :: handle
    integer(i4b), allocatable, dimension(:)        :: ind, n_per_chain
    integer(i4b), allocatable, dimension(:,:)      :: mask2map
    real(dp),     allocatable, dimension(:)        :: map_1D, mean_1D, W
    real(dp),     allocatable, dimension(:,:)      :: mean, cov, map, beam, alms_mean, map_1d_, mask
    real(dp),     pointer,     dimension(:,:)      :: pixwin
    real(sp),     allocatable, dimension(:,:)      :: alms
    complex(dpc), allocatable, dimension(:,:,:)    :: alms_cmplx
    real(sp),     allocatable, dimension(:,:,:)    :: fg_amp

    if (iargc() < 5) then
       if (trim(filetype) == 'pix') then
          write(*,*) '    Compute mean and covariance from pixel amplitude chain files'
          write(*,*) '    Options:  [outprefix] [nside] [nmaps] [ordering] [thiscomp] [ncomp] [burnin]'
          write(*,*) '                 [T regnoise] [QU regnoise] [seed] [compute_covar] [maskfile] [chain1] [chain2] ...'
       else if (trim(filetype) == 'alm') then
          write(*,*) '    Compute mean and covariance from CMB alms chain files'
          write(*,*) '    Options:  [outprefix] [nside] [nmaps] [ordering] [lmax] [beamfile] [burnin]'
          write(*,*) '                 [T regnoise] [QU regnoise] [seed] [compute_covar] [chain1] [chain2] ...'
       else
          write(*,*) 'Unknown file type = ', trim(filetype)
       end if
       stop
    end if

    call getarg(2,prefix)
    call getarg(3,temp)
    read(temp,*) nside
    call getarg(4,temp)
    read(temp,*) nmaps
    call getarg(5,temp)
    read(temp,*) ordering
    if (trim(filetype) == 'pix') then
       call getarg(6,temp)
       read(temp,*) comp
       call getarg(7,temp)
       read(temp,*) ncomp
    else 
       call getarg(6,temp)
       read(temp,*) lmax
       call getarg(7,temp)
       read(temp,*) beamfile
       numcomp = (lmax+1)**2
    end if
    call getarg(8,temp)
    read(temp,*) burnin
    call getarg(9,temp)
    read(temp,*) sigma_reg(1)
    call getarg(10,temp)
    read(temp,*) sigma_reg(2)
    sigma_reg(3) = sigma_reg(2)
    call getarg(11,temp)
    read(temp,*) seed
    call getarg(12,temp)
    read(temp,*) compute_covar
    call getarg(13,maskfile)
    numchain = iargc()-13
    npix     = 12*nside**2
    unit     = 58
    nmax     = 10005000

    if (trim(filetype) == 'pix') then
       allocate(fg_amp(0:npix-1,nmaps,ncomp))
    else
       allocate(alms(numcomp,nmaps), alms_cmplx(nmaps,0:lmax,0:lmax))
       allocate(alms_mean(numcomp,nmaps))
       allocate(beam(0:lmax,nmaps))
       call read_pixwin(nside, nmaps, pixwin)
       if (trim(beamfile) == 'none') then
          beam = 1.d0
       else
          call read_beam(beamfile, beam)
       end if
       beam = beam * pixwin(0:lmax,1:nmaps)
       deallocate(pixwin)
       alms_mean = 0.d0
    end if

    ! Read maskfile
    allocate(mask(0:12*nside**2-1,nmaps))
    if (trim(maskfile) /= 'none') then
       call read_map(maskfile, mask)
    else
       mask = 1.d0
    end if
    n = count(mask > 0.5d0)
    allocate(mask2map(n,2))
    k = 1
    do j = 1, nmaps
       do i = 0, 12*nside**2-1
          if (mask(i,j) > 0.5d0) then
             mask2map(k,1) = i
             mask2map(k,2) = j
             k             = k+1
          end if
       end do
    end do

    ! Compute mean 
    allocate(mean(0:npix-1,nmaps), mean_1D(n), cov(n,n), map(0:npix-1,nmaps), map_1D(n), map_1D_(n,1))
    allocate(n_per_chain(numchain))
    mean  = 0.d0
    nreal = 0
    n_per_chain = 0
    do chain = 1, numchain
       
       call getarg(13+chain,infile)
       open(unit, file=trim(infile), form='unformatted')
       read(unit) k
       counter = 0
       do while (.true.)
          if (trim(filetype) == 'pix') then
             read(unit,end=48) fg_amp
             counter = counter+1
             n_per_chain(chain) = n_per_chain(chain) + 1
             if (counter > burnin .and. counter <= nmax) then
                mean  = mean + fg_amp(:,:,comp)
                nreal = nreal + 1
                n_per_chain(chain) = n_per_chain(chain) + 1
             end if
          else
             read(unit,end=48) alms
             counter = counter+1
             if (counter > burnin .and. counter <= nmax) then
                alms_mean = alms_mean + alms
                nreal = nreal+1
                n_per_chain(chain) = n_per_chain(chain) + 1
             end if
          end if
          write(*,*) 'Computing mean -- chain = ', chain, ', n = ', nreal
       end do
48     close(unit)
    end do

    if (trim(filetype) == 'alm') then
       call convert_real_to_complex_alms_dp(alms_mean, alms_cmplx)
       do k = 1, nmaps
          do l = 0, lmax
             alms_cmplx(k,l,0:l) = alms_cmplx(k,l,0:l) * beam(l,k) 
          end do
       end do

       if (nmaps == 1) then
          call alm2map(nside, lmax, lmax, alms_cmplx, mean(:,1))
       else
          call alm2map(nside, lmax, lmax, alms_cmplx, mean)
       end if
    end if
    mean = mean / nreal
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_ring2nest(nside, mean(:,i))
       end do
    end if

    if (compute_covar) then
       ! Compute covariance matrix
       cov   = 0.d0
       m     = 0
       do chain = 1, numchain
          
          call getarg(13+chain,infile)       
          open(unit, file=trim(infile), form='unformatted')
          read(unit) k
          do counter = 1, n_per_chain(chain)
             if (trim(filetype) == 'pix') then
                read(unit,end=49) fg_amp
                if (counter <= burnin .or. counter > nmax) cycle
                map = fg_amp(:,:,comp)
             else
                read(unit,end=49) alms
                if (counter <= burnin .or. counter > nmax) cycle
                alms_mean = alms
                call convert_real_to_complex_alms_dp(alms_mean, alms_cmplx)
                do k = 1, nmaps
                   do l = 0, lmax
                      alms_cmplx(k,l,0:l) = alms_cmplx(k,l,0:l) * beam(l,k) 
                   end do
                end do
                if (nmaps == 1) then
                   call alm2map(nside, lmax, lmax, alms_cmplx, map(:,1))
                else
                   call alm2map(nside, lmax, lmax, alms_cmplx, map)
                end if
             end if
             m = m+1
             write(*,*) 'Computing covmat -- iter = ', m, ' of ', nreal
             
             if (ordering == 2) then
                do j = 1, nmaps
                   call convert_ring2nest(nside, map(:,j))
                end do
             end if
             
             ! Linearize map
             do j = 1, n
                map_1D(j) = map(mask2map(j,1),mask2map(j,2)) - mean(mask2map(j,1),mask2map(j,2))
             end do
             
             ! Compute outer product 
             !$OMP PARALLEL PRIVATE(i,j) 
             !$OMP DO SCHEDULE(guided)                   
             do i = 1, n
                do j = i, n
                   cov(j,i) = cov(j,i) + map_1D(i)*map_1D(j)
                end do
             end do
             !$OMP END DO                                                                                      
             !$OMP END PARALLEL             
          end do
49        close(unit)
       end do
       cov = cov / (nreal-1)
       do i = 1, n
          do j = 1, i
             cov(j,i) = cov(i,j)
          end do
       end do
    end if

    ! Apply mask
    where (mask < 0.5d0)
       mean = -1.6375d30
    end where
    
    ! Add regularization noise if requested
    call rand_init(handle, seed)
    if (any(sigma_reg > 0.d0)) then
       do i = 1, n
          j = mask2map(i,1)
          k = mask2map(i,2)
          if (sigma_reg(k) > 0.d0) then
             mean(j,k) = mean(j,k) + sigma_reg(k)*rand_gauss(handle)
             cov(i,i)  = cov(i,i) + sigma_reg(k)**2
          end if
       end do
    end if

    ! Output to standard formats
    outfile = trim(prefix) // '_mean.fits'
    call write_map3(outfile, mean, ordering=ordering)    

    if (compute_covar) then
       outfile = trim(prefix) // '_N.unf'
       open(unit,file=trim(outfile),form='unformatted')
       write(unit) n
       write(unit) ordering ! Ordering
       write(unit) nmaps ! Polarization status
       do i = 1, n
          write(unit) cov(:,i)
       end do
       write(unit) .false. ! Not inverse
       close(unit)
       
!!$       open(54,file='test.dat')
!!$       do i = 1, n
!!$          write(58,*) i, cov(i,7418)
!!$       end do
!!$       close(54)
       
!!$       do k = 1, nmaps
!!$          do p = nint(0.3*npix), npix-1, nint(0.3*npix)
!!$             do i = 1, nmaps
!!$                mean(:,i) = cov((i-1)*npix+1:i*npix,(k-1)*npix+1+p)
!!$             end do
!!$             if (any(mean /= 0.d0)) then
!!$                call int2string(p,pix_text)
!!$                outfile = trim(prefix) // '_' // s_text(k) // &
!!$                     & '_pix' // pix_text // '.fits'
!!$                call write_map3(outfile, mean, ordering=ordering)    
!!$             end if
!!$          end do
!!$       end do
       
       k = 0
       mean = -1.6375d30
       do i = 1, n
          j = mask2map(i,1)
          k = mask2map(i,2)
          mean(j,k) = sqrt(cov(i,i))
       end do
       outfile = trim(prefix) // '_rms.fits'
       call write_map3(outfile, mean, ordering=ordering)    
       
       write(*,*) 'Computing eigenspectrum'
       allocate(W(n))
       call get_eigenvalues(cov, W)
       open(unit,file=trim(prefix)//'_W.dat')
       do i = 1, n
          if (W(i) == 0.d0) then
             write(unit,*) '# ', i, 0., 0.
          else
             write(unit,*) i, abs(real(W(i),sp)), real(W(i),sp)
          end if
       end do
       close(unit)
       j = 1
       do while (W(j) == 0.d0)
          j = j+1
       end do
       write(*,*) '   Condition number = ', W(n)/W(j)
    end if

  end subroutine chain2mean_cov

  subroutine maps2mean_cov
    implicit none

    character(len=256) :: outfile, temp, infile, prefix, beamfile, maskfile
    character(len=4)   :: chain_text
    character(len=5)   :: pix_text
    character(len=1)   :: s_text(3) = ['T', 'Q', 'U']
    integer(i4b)       :: i, j, k, l, m, n, p, unit, nreal, ordering, chain, lmax, counter, nval
    integer(i4b)       :: nside, npix, nmaps, comp, ncomp, burnin, numiter, numchain, seed, numcomp, nmax
    real(dp)           :: sigma_reg(3)
    logical(lgt)       :: compute_covar
    type(planck_rng)   :: handle
    integer(i4b), allocatable, dimension(:)        :: ind, n_per_chain
    integer(i4b), allocatable, dimension(:,:)      :: mask2map
    real(dp),     allocatable, dimension(:)        :: map_1D, mean_1D, W
    real(dp),     allocatable, dimension(:,:)      :: mean, cov, map, beam, alms_mean, map_1d_, mask
    real(dp),     pointer,     dimension(:,:)      :: pixwin
    real(sp),     allocatable, dimension(:,:)      :: alms
    complex(dpc), allocatable, dimension(:,:,:)    :: alms_cmplx
    real(sp),     allocatable, dimension(:,:,:)    :: fg_amp

    if (iargc() < 8) then
       write(*,*) '    Compute mean and covariance from pixel amplitude chain files'
       write(*,*) '    Options:  [outprefix] [T regnoise] [QU regnoise] [seed] [compute_covar] [maskfile] [map1] [map22] ...'
       stop
    end if

    call getarg(2,prefix)
    call getarg(3,temp)
    read(temp,*) sigma_reg(1)
    call getarg(4,temp)
    read(temp,*) sigma_reg(2)
    sigma_reg(3) = sigma_reg(2)
    call getarg(5,temp)
    read(temp,*) seed
    call getarg(6,temp)
    read(temp,*) compute_covar
    call getarg(7,maskfile)
    nreal    = iargc()-7

    i = getsize_fits(trim(maskfile), nside=nside, nmaps=nmaps, ordering=ordering)
    npix = 12*nside**2

    ! Read maskfile
    allocate(mask(0:12*nside**2-1,nmaps))
    call read_map(maskfile, mask)
    n = count(mask > 0.5d0)
    allocate(mask2map(n,2))
    k = 1
    do j = 1, nmaps
       do i = 0, 12*nside**2-1
          if (mask(i,j) > 0.5d0) then
             mask2map(k,1) = i
             mask2map(k,2) = j
             k             = k+1
          end if
       end do
    end do

    ! Compute mean 
    allocate(mean(0:npix-1,nmaps), mean_1D(n), cov(n,n), map(0:npix-1,nmaps), map_1D(n), map_1D_(n,1))
    mean  = 0.d0
    do i = 1, nreal       
       call getarg(7+i,infile)
       call read_map(infile, map)
       mean  = mean + map
       write(*,*) 'Computing mean -- i = ', i, ', n = ', nreal
    end do
    mean = mean / nreal

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_ring2nest(nside, mean(:,i))
       end do
    end if

    if (compute_covar) then
       ! Compute covariance matrix
       cov   = 0.d0
       do k = 1, nreal
          call getarg(7+k,infile)
          call read_map(infile, map)
          write(*,*) 'Computing covmat -- i = ', k, ' of ', nreal
             
          if (ordering == 2) then
             do j = 1, nmaps
                call convert_ring2nest(nside, map(:,j))
             end do
          end if
             
          ! Linearize map
          do j = 1, n
             map_1D(j) = map(mask2map(j,1),mask2map(j,2)) - mean(mask2map(j,1),mask2map(j,2))
          end do
             
          ! Compute outer product 
          !$OMP PARALLEL PRIVATE(i,j) 
          !$OMP DO SCHEDULE(guided)                   
          do i = 1, n
             do j = i, n
                cov(j,i) = cov(j,i) + map_1D(i)*map_1D(j)
             end do
          end do
          !$OMP END DO                             
          !$OMP END PARALLEL           

       end do
       cov = cov / (nreal-1)
       do i = 1, n
          do j = 1, i
             cov(j,i) = cov(i,j)
          end do
       end do
    end if

    ! Apply mask
    where (mask < 0.5d0)
       mean = -1.6375d30
    end where
    
    ! Add regularization noise if requested
    call rand_init(handle, seed)
    if (any(sigma_reg > 0.d0)) then
       do i = 1, n
          j = mask2map(i,1)
          k = mask2map(i,2)
          if (sigma_reg(k) > 0.d0) then
             mean(j,k) = mean(j,k) + sigma_reg(k)*rand_gauss(handle)
             cov(i,i)  = cov(i,i) + sigma_reg(k)**2
          end if
       end do
    end if

    ! Output to standard formats
    outfile = trim(prefix) // '_mean.fits'
    call write_map3(outfile, mean, ordering=ordering)    

    if (compute_covar) then
       outfile = trim(prefix) // '_N.unf'
       open(unit,file=trim(outfile),form='unformatted')
       write(unit) n
       write(unit) ordering ! Ordering
       write(unit) nmaps ! Polarization status
       do i = 1, n
          write(unit) cov(:,i)
       end do
       write(unit) .false. ! Not inverse
       close(unit)
       
!!$       open(54,file='test.dat')
!!$       do i = 1, n
!!$          write(58,*) i, cov(i,7418)
!!$       end do
!!$       close(54)
       
!!$       do k = 1, nmaps
!!$          do p = nint(0.3*npix), npix-1, nint(0.3*npix)
!!$             do i = 1, nmaps
!!$                mean(:,i) = cov((i-1)*npix+1:i*npix,(k-1)*npix+1+p)
!!$             end do
!!$             if (any(mean /= 0.d0)) then
!!$                call int2string(p,pix_text)
!!$                outfile = trim(prefix) // '_' // s_text(k) // &
!!$                     & '_pix' // pix_text // '.fits'
!!$                call write_map3(outfile, mean, ordering=ordering)    
!!$             end if
!!$          end do
!!$       end do
       
       k = 0
       mean = -1.6375d30
       do i = 1, n
          j = mask2map(i,1)
          k = mask2map(i,2)
          mean(j,k) = sqrt(cov(i,i))
       end do
       outfile = trim(prefix) // '_rms.fits'
       call write_map3(outfile, mean, ordering=ordering)    
       
       write(*,*) 'Computing eigenspectrum'
       allocate(W(n))
       call get_eigenvalues(cov, W)
       open(unit,file=trim(prefix)//'_W.dat')
       do i = 1, n
          if (W(i) == 0.d0) then
             write(unit,*) '# ', i, 0., 0.
          else
             write(unit,*) i, abs(real(W(i),sp)), real(W(i),sp)
          end if
       end do
       close(unit)
       j = 1
       do while (W(j) == 0.d0)
          j = j+1
       end do
       write(*,*) '   Condition number = ', W(n)/W(j)
    end if

  end subroutine maps2mean_cov


  subroutine mean_rms2fits
    implicit none

    character(len=256) :: outfile, temp, infile, prefix, meanfile, rmsfile, outprefix
    character(len=4)   :: chain_text
    character(len=2)   :: comp_text
    integer(i4b)       :: numtemp, numband, numiter, numchain, i, j, k, l, n, nside, nmaps
    integer(i4b)       :: ncomp, comp, npix, burnin, firstchain, lastchain
    real(sp),     allocatable, dimension(:,:,:)  :: vals
    real(sp),     allocatable, dimension(:,:,:)  :: mean, rms

    if (iargc() /= 9) then
       write(*,*) '    Compute mean and stddev maps from chain files'
       write(*,*) '    Options:  [chain_prefix] [nside] [nmaps] [ncomp] '
       write(*,*) '                  [firstchain] [lastchain] [burnin] [outprefix]'
       stop
    end if

    call getarg(2,prefix)
    call getarg(3,temp)
    read(temp,*) nside
    call getarg(4,temp)
    read(temp,*) nmaps
    call getarg(5,temp)
    read(temp,*) ncomp
    call getarg(6,temp)
    read(temp,*) firstchain
    call getarg(7,temp)
    read(temp,*) lastchain
    call getarg(8,temp)
    read(temp,*) burnin
    call getarg(9, outprefix)
    npix = 12*nside**2

    allocate(vals(0:npix-1,nmaps,ncomp))
    allocate(mean(0:npix-1,nmaps,ncomp), rms(0:npix-1,nmaps,ncomp))
    mean = 0.d0
    n    = 0
    do j = firstchain, lastchain
       call int2string(j, chain_text)
       infile = trim(prefix) // '_no' // chain_text // '.unf'
       open(unit, file=trim(infile), form='unformatted')
       read(unit) k
       l = 0
       do while (.true.)
          read(unit,end=110) vals
          l = l+1
          if (l <= burnin) cycle
          write(*,*) 'Computing mean -- ', j, l
          mean = mean + vals
          n    = n    + 1
       end do
110    close(unit)

    end do
    mean = mean / n

    rms  = 0.d0
    do j = firstchain, lastchain
       call int2string(j, chain_text)
       infile = trim(prefix) // '_no' // chain_text // '.unf'
       open(unit, file=trim(infile), form='unformatted')
       read(unit) k
       l = 0
       do while (.true.)
          read(unit,end=111) vals
          l = l+1
          if (l <= burnin) cycle
          write(*,*) 'Computing rms  -- ', j, l
          rms = rms + (vals-mean)**2
       end do
111    close(unit)

    end do
    rms = sqrt(rms/(n-1))

    do i = 1, ncomp
       call int2string(i, comp_text)
       outfile = trim(outprefix) // '_c' // comp_text // '_mean.fits'
       call write_map3(outfile, real(mean(:,:,i),dp))
       outfile = trim(outprefix) // '_c' // comp_text // '_rms.fits'
       call write_map3(outfile,  real(rms(:,:,i),dp))
    end do

    deallocate(vals, mean, rms)

  end subroutine mean_rms2fits



  subroutine map2fits(filetype)
    implicit none

    character(len=*), intent(in) :: filetype

    character(len=256) :: outfile, temp, infile
    character(len=4)   :: chain_text
    integer(i4b)       :: component, num_component, nside, nmaps, numiter, numburnin
    integer(i4b)       :: first_chain, last_chain
    integer(i4b)       :: i, j, k
    real(sp),     allocatable, dimension(:,:,:,:)  :: cls
    real(sp),     allocatable, dimension(:,:,:)    :: map_in

    if (iargc() == 10) then
       write(*,*) '  Create sky map FITS file (warning -- can become very big!)'
       write(*,*) '    Options:  [outfile] [component] [num_component] [nside] [nmaps] '
       write(*,*) '              [number of iterations] [first_chain] [last_chain]'
       stop
    end if

    call getarg(2,outfile)
    call getarg(3,temp)
    read(temp,*) component
    call getarg(4,temp)
    read(temp,*) num_component
    call getarg(5,temp)
    read(temp,*) nside
    call getarg(6,temp)
    read(temp,*) nmaps
    call getarg(7,temp)
    read(temp,*) numiter
    call getarg(8,temp)
    read(temp,*) numburnin
    call getarg(9,temp)
    read(temp,*) first_chain
    call getarg(10,temp)
    read(temp,*) last_chain

    allocate(cls(0:12*nside**2-1, nmaps, last_chain-first_chain+1, numiter-numburnin))
    allocate(map_in(0:12*nside**2-1, nmaps, num_component))

    do j = first_chain, last_chain

       write(*,*) 'Chain = ', j
       call int2string(j, chain_text)

       if (trim(filetype) == 'amp') then
          infile = 'chain_fg_amps_no' // chain_text // '.unf'
       else if (trim(filetype) == 'ind') then
          infile = 'chain_fg_ind_no' // chain_text // '.unf'
       end if
       open(unit, file=trim(infile), form='unformatted')
       read(unit) k
       do i = 1, numiter
          if (mod(i,100) == 0) write(*,*) 'Map = ', i
          read(unit) map_in
          if (i > numburnin) cls(:,:,j-first_chain+1,i-numburnin) = map_in(:,:,component)
       end do
       close(unit)
    end do

    call write_fgt_resfile(outfile, cls)

    deallocate(cls, map_in)

  end subroutine map2fits

  subroutine cmb2fits
    implicit none

    character(len=256) :: outfile, temp, infile
    character(len=4)   :: chain_text
    integer(i4b)       :: lmax, nside, nmaps, numiter, numchain, i, j, k, l, m, numcomp
    real(dp)           :: fwhm
    real(sp),     allocatable, dimension(:,:,:,:)  :: cls    
    real(sp),     allocatable, dimension(:,:)      :: map, alms_in
    real(dp),     pointer,     dimension(:,:)      :: ringweights, pixwin
    real(dp),     allocatable, dimension(:,:)      :: gb
    complex(spc), allocatable, dimension(:,:,:)    :: alms, alms_sum

    if (iargc() /= 8) then
       write(*,*) '    Create CMB map data base'
       write(*,*) '    Options: [outfile] [lmax] [nside] [nmaps] [fwhm] '
       write(*,*) '             [number of iterations] [number of chains]'
       stop
    end if

    call getarg(2,outfile)
    call getarg(3,temp)
    read(temp,*) lmax
    call getarg(4,temp)
    read(temp,*) nside
    call getarg(5,temp)
    read(temp,*) nmaps
    call getarg(6,temp)
    read(temp,*) fwhm
    call getarg(7,temp)
    read(temp,*) numiter
    call getarg(8,temp)
    read(temp,*) numchain


    numcomp = (lmax+1)**2
    allocate(cls(0:12*nside**2-1, nmaps, numchain, numiter))
    allocate(alms_in(numcomp, nmaps))
    allocate(alms(nmaps, 0:lmax, 0:lmax))
    allocate(map(0:12*nside**2-1, nmaps))
    allocate(gb(0:lmax,nmaps))
    call read_pixwin(nside, nmaps, pixwin)
    call gaussbeam(fwhm, lmax, gb)

    do j = 1, numchain

       call int2string(j, chain_text)

       infile = 'chain_alms_no' // chain_text // '.unf'
       open(unit, file=trim(infile), form='unformatted')
       read(unit) k
       do i = 1, numiter
          write(*,*) j, i 
          read(unit) alms_in
          call convert_real_to_complex_alms(alms_in, alms)

          do k = 1, nmaps
             do l = 0, lmax
                do m = 0, l
                   alms(k,l,m) = alms(k,l,m) * gb(l,k) * pixwin(l,k)
                end do
             end do
          end do

          if (nmaps == 1) then
             call alm2map(nside, lmax, lmax, alms, map(:,1))
          else
             call alm2map(nside, lmax, lmax, alms, map)
          end if
          cls(:,:,j,i) = map
       end do
       close(unit)

    end do

    call write_fgt_resfile(outfile, cls)

    deallocate(cls, alms, alms_in, pixwin, gb)

  end subroutine cmb2fits

  subroutine chisq2fits
    implicit none

    character(len=256) :: outfile, temp, infile
    character(len=4)   :: chain_text
    integer(i4b)       :: numiter, numchain, i, j, k
    logical(lgt)       :: exist
    real(sp),     allocatable, dimension(:,:,:,:)  :: cls    

    if (iargc() /= 4) then
       write(*,*) '    Create chi-square FITS database'
       write(*,*) '    Options:  [outfile] [number of iterations] [number of chains]'
       stop
    end if

    call getarg(2,outfile)
    call getarg(3,temp)
    read(temp,*) numiter
    call getarg(4,temp)
    read(temp,*) numchain

    allocate(cls(1, 1, numchain, numiter))
    cls = 0.
    do j = 1, numchain
       call int2string(j, chain_text)
       infile = 'chain_chisq_no' // chain_text // '.unf'
       inquire(file=trim(infile), exist=exist)
       open(unit, file=trim(infile), form='unformatted')
       read(unit) k
       do i = 1, numiter
          read(unit) cls(1,1,j,i)
       end do
       close(unit)
    end do

    call write_resfile(outfile, 1, 1, numiter, numchain, &
         & "Chi-squares from Gibbs samples", cls)

    deallocate(cls)

  end subroutine chisq2fits

  subroutine cmb_mean2fits(filetype)
    implicit none

    character(len=*), intent(in) :: filetype

    character(len=256) :: outfile, temp, maskfile, infile
    character(len=4)   :: chain_text
    integer(i4b)       :: lmax_file, lmin, lmax, nside, nmaps, numiter, numchain
    integer(i4b)       :: numburnin, numcomp, i, j, k, l, m
    real(dp)           :: fwhm
    real(sp),     allocatable, dimension(:,:)      :: map, alms_in
    real(dp),     allocatable, dimension(:,:)      :: mask, map_dp
    real(sp),     allocatable, dimension(:,:,:,:)  :: cls
    real(dp),     pointer,     dimension(:,:)      :: ringweights, pixwin
    real(dp),     allocatable, dimension(:,:)      :: gb
    complex(spc), allocatable, dimension(:,:,:)    :: alms, alms_sum

    if (iargc() /= 12) then
       write(*,*) '    Compute CMB map average file'
       write(*,*) '    Options: [outfile] [lmax_file] [lmin] [lmax] [nside] [nmaps] [fwhm] '
       write(*,*) '             [number of iterations] [number of chains] [burn_in] '
       write(*,*) '             [maskfile (or "none")]'
       stop
    end if

    call getarg(2,outfile)
    call getarg(3,temp)
    read(temp,*) lmax_file
    call getarg(4,temp)
    read(temp,*) lmin
    call getarg(5,temp)
    read(temp,*) lmax
    call getarg(6,temp)
    read(temp,*) nside
    call getarg(7,temp)
    read(temp,*) nmaps
    call getarg(8,temp)
    read(temp,*) fwhm
    call getarg(9,temp)
    read(temp,*) numiter
    call getarg(10,temp)
    read(temp,*) numchain
    call getarg(11,temp)
    read(temp,*) numburnin
    call getarg(12,maskfile)

    ! Read mask
    if (trim(maskfile) /= 'none') then
       allocate(mask(0:12*nside**2-1,nmaps))
       call read_map(maskfile, mask)
    end if

    numcomp = (lmax_file+1)**2
    allocate(alms_in(numcomp, nmaps))
    allocate(alms(nmaps, 0:lmax_file, 0:lmax_file))
    allocate(alms_sum(nmaps, 0:lmax, 0:lmax))
    allocate(map(0:12*nside**2-1, nmaps))
    allocate(map_dp(0:12*nside**2-1, nmaps))
    allocate(gb(0:lmax,nmaps))
    call read_pixwin(nside, nmaps, pixwin)
    call gaussbeam(fwhm, lmax, gb)

    alms_sum = cmplx(0., 0.)
    m        = 0
    do j = 1, numchain
       call int2string(j, chain_text)
       infile = 'chain_alms_no' // chain_text // '.unf'
       open(unit, file=trim(infile), form='unformatted')
       read(unit) k
       do i = 1, numiter
          read(unit) alms_in
          call convert_real_to_complex_alms(alms_in, alms)
          alms_sum = alms_sum + alms(:,0:lmax,0:lmax)
          m        = m + 1
       end do
       close(unit)
    end do

    alms_sum = alms_sum / real(m,dp)
    do k = 1, nmaps
       do l = 0, lmax
          do m = 0, l
             alms_sum(k,l,m) = alms_sum(k,l,m) * gb(l,k) * pixwin(l,k)
          end do
       end do
    end do

    ! Remove multipoles up to lmin
    alms_sum(:,0:lmin-1,:) = cmplx(0.d0, 0.d0)

    if (trim(filetype) == 'IQU') then
       if (nmaps == 1) then
          call alm2map(nside, lmax, lmax, alms_sum, map(:,1))
       else
          call alm2map(nside, lmax, lmax, alms_sum, map)
       end if
    else if (trim(filetype) == 'TEB') then
       ! Make TEB maps
       do i = 1, nmaps
          call alm2map(nside, lmax, lmax, alms_sum(i:i,:,:), map(:,i))
       end do
    else
       write(*,*) 'Unknown map type.'
       stop
    end if

    map_dp = real(map,dp)

    if (trim(maskfile) /= 'none') then
       where (mask == 0.d0) map_dp = -1.6375d30
    end if

    call write_map3(outfile, map_dp)

  end subroutine cmb_mean2fits

  subroutine cl2fits_singlechain(filetype)
    implicit none

    character(len=*), intent(in) :: filetype

    character(len=256) :: outfile, temp, infile
    character(len=4)   :: chain_text
    integer(i4b)       :: max_num_samples, numchain, lmax, nspec, numburnin, stat
    integer(i4b)       :: i, j, k, l, counter, counter_local
    logical(lgt)       :: exist
    real(sp),     allocatable, dimension(:,:,:,:)  :: cls    
    real(sp),     allocatable, dimension(:,:)      :: in_cls

    if (iargc() /= 7) then
       write(*,*) '    Create power spectrum FITS files; collapse all samples into single chain'
       write(*,*) '    Options: [outfile] [maximum number of samples]'
       write(*,*) '             [number of chains] [lmax] [nspec] [burnin]'
       stop
    end if

    call getarg(2,outfile)
    call getarg(3,temp)
    read(temp,*) max_num_samples
    call getarg(4,temp)
    read(temp,*) numchain
    call getarg(5,temp)
    read(temp,*) lmax
    call getarg(6,temp)
    read(temp,*) nspec
    call getarg(7,temp)
    read(temp,*) numburnin

    allocate(cls(0:lmax, nspec, 1, 0:max_num_samples))
    allocate(in_cls(0:lmax, nspec))

    cls = 0.
    do l = 0, lmax
       cls(l,:,:,0) = real(l,sp)
    end do

    counter = 1
    do j = 1, numchain
       call int2string(j, chain_text)
       if (trim(filetype) == 'cls') then
          infile = 'chain_cls_no' // chain_text // '.unf'
       else if (trim(filetype) == 'sigma') then
          infile = 'chain_sigma_no' // chain_text // '.unf'
       else if (trim(filetype) == 'mcmc') then
          infile = 'chain_mcmc_no' // chain_text // '.unf'
       else
          write(*,*) 'Unknown file type.'
          stop
       end if

       inquire(file=trim(infile), exist=exist)
       open(unit, file=trim(infile), form='unformatted')
       read(unit) k
       counter_local = 0
       do while (.true.)
          read(unit,iostat=stat) in_cls
          if (stat == 0) then
             counter_local = counter_local + 1
             if (counter_local > numburnin) then
                cls(:,:,1,counter) = in_cls
                counter            = counter + 1
             end if
          else
             close(unit)
             exit
          end if
       end do

    end do

    write(*,*) 'Total number of samples = ', counter-1
    call write_resfile(outfile, lmax, nspec, counter-1, 1, &
         & "Gibbs sampled power spectra", cls(:,:,:,0:counter-1))

    deallocate(cls)

  end subroutine cl2fits_singlechain

  subroutine trim_file
    implicit none

    character(len=256) :: infile, outfile, temp
    integer(i4b)       :: first_chain, last_chain, first_sample, last_sample
    integer(i4b)       :: numchain, numiter, nspec, lmax, thinstep
    real(sp), pointer,     dimension(:,:,:,:) :: cls_pointer
    real(sp), allocatable, dimension(:,:,:,:) :: buffer

    if (iargc() /= 8) then
       write(*,*) '    Trim result file'
       write(*,*) '    Options: [infile] [outfile] [first chain] [last chain] '
       write(*,*) '             [first sample] [last sample] [thin step]'
       stop
    end if

    call getarg(2,infile)
    call getarg(3,outfile)
    call getarg(4,temp)
    read(temp,*) first_chain
    call getarg(5,temp)
    read(temp,*) last_chain
    call getarg(6,temp)
    read(temp,*) first_sample
    call getarg(7,temp)
    read(temp,*) last_sample
    call getarg(8,temp)
    read(temp,*) thinstep

    unit = 58
    call read_resfile(infile, unit, nspec, lmax, numchain, numiter, cls_pointer)

    numchain = last_chain-first_chain+1
    numiter  = (last_sample-first_sample+1)/thinstep
    allocate(buffer(0:lmax,nspec,numchain,numiter+1))
    buffer(:,:,:,1) = cls_pointer(:,:,first_chain:last_chain,1)
    buffer(:,:,1:numchain,2:numiter+1:thinstep) = &
         & cls_pointer(:,:,first_chain:last_chain,first_sample+1:last_sample+1:thinstep)
    call write_resfile(outfile, lmax, nspec, numiter, numchain, 'Commander', buffer)
    deallocate(cls_pointer, buffer)

  end subroutine trim_file

  subroutine flatten_chain
    implicit none

    character(len=256) :: infile, outfile, temp
    integer(i4b)       :: first_chain, last_chain, first_sample, last_sample
    integer(i4b)       :: numchain, numiter, numiter_new, numburnin, nspec, lmax, i, j, k
    real(sp), pointer,     dimension(:,:,:,:) :: cls_pointer
    real(sp), allocatable, dimension(:,:,:,:) :: buffer

    if (iargc() /= 4) then
       write(*,*) '    Flatten multi-chain FITS file'
       write(*,*) '    Options: [infile] [outfile] [burn-in per chain]'
       stop
    end if

    call getarg(2, infile)
    call getarg(3, outfile)
    call getarg(4, temp)
    read(temp, *) numburnin

    unit = 58
    call read_resfile(infile, unit, nspec, lmax, numchain, numiter, cls_pointer)
    numiter_new = 0
    do i = 1, numchain 
       if (cls_pointer(0, 1, i, 1) > numburnin) then
          numiter_new = numiter_new + cls_pointer(0, 1, i, 1) - numburnin
       end if
    end do

    allocate(buffer(0:lmax, nspec, 1, 0:numiter_new))
    buffer = 0
    !Saves number of iterations in first sample column
    buffer(:, :, 1, 0) = numiter_new
    k = 0
    do i = 1, numchain
       if (cls_pointer(0, 1, i, 1) > numburnin) then
          do j = numburnin + 1, nint(cls_pointer(0, 1, i, 1))
             k = k + 1
             buffer(0:lmax, :, 1, k) = cls_pointer(0:lmax, :, i, j)
          end do
       end if
    end do
    call write_resfile(outfile, lmax, nspec, numiter_new, 1, 'Commander', buffer)

    deallocate(buffer, cls_pointer)

  end subroutine flatten_chain

  subroutine cl2fits_from_files
    implicit none

    character(len=256) :: outfile, temp, infile
    character(len=4)   :: chain_text
    integer(i4b)       :: numiter, lmax, nspec, numburnin, numchain, stat
    integer(i4b)       :: i, j, k, l, max_num_samples
    logical(lgt)       :: reached_limit, exist
    real(sp),     allocatable, dimension(:,:,:,:)  :: cls, pre_cls

    if (iargc() <= 6) then
       write(*,*) '    Create power spectrum FITS file from multiple chain files'
       write(*,*) '    Options:  [outfile] [max number of iterations per chain] '
       write(*,*) '              [lmax] [nspec] [burnin] [chain1] [chain2] ...'
       stop
    end if

    call getarg(2,outfile)
    call getarg(3,temp)
    read(temp,*) numiter
    call getarg(4,temp)
    read(temp,*) lmax
    call getarg(5,temp)
    read(temp,*) nspec
    call getarg(6,temp)
    read(temp,*) numburnin
    numchain = iargc()-6

    allocate(pre_cls(0:lmax, nspec, numchain, 0:numiter))
    pre_cls = 0.d0
    do j = 1, numchain

       call int2string(j, chain_text)
       call getarg(6+j, infile)

       inquire(file=trim(infile), exist=exist)
       open(unit, file=trim(infile), form='unformatted')
       read(unit) k
       l = 1
       reached_limit = .false.
       do i = 1, numiter
          read(unit, iostat=stat) pre_cls(:,:,j,l)
          if (stat < 0) then
             pre_cls(:, :, j, 0) = l - 1
             reached_limit = .true.
             exit
          end if
          if (i > numburnin) l = l + 1
       end do

       if (.not. reached_limit) pre_cls(:, :, j, 0) = numiter - numburnin

       close(unit)

    end do

    max_num_samples = numiter
    do while (all(pre_cls(:,:,:,max_num_samples) == 0.d0))
       max_num_samples = max_num_samples-1
    end do

    allocate(cls(0:lmax, nspec, numchain, 0:max_num_samples))
    cls = pre_cls(:, :, :, 0:max_num_samples)

    call write_resfile(outfile, lmax, nspec, max_num_samples, numchain, &
         & "Gibbs sampled power spectra", cls)

    deallocate(pre_cls)

    l = sum(cls(0,1,:,0))
    allocate(pre_cls(0:lmax,nspec,1,0:l))
    l = 1
    do i = 1, numchain
       do j = 1, nint(cls(0,1,i,0))
          pre_cls(:,:,1,l) = cls(:,:,i,j)
          l                = l+1
       end do
    end do

    call write_resfile(trim(outfile)//'_flat', lmax, nspec, l-1, 1, &
         & "Gibbs sampled power spectra", pre_cls)

    deallocate(cls)

  end subroutine cl2fits_from_files

  subroutine cl2ascii_binned
    implicit none

    character(len=256) :: infile, binfile, outfile, temp, line
    integer(i4b)       :: lmin_col, lmax_col, numburnin
    integer(i4b)       :: nspec, lmax, numchain, numiter, i, j, k, s, bins(2)
    real(sp), pointer, dimension(:,:,:,:) :: cls_pointer
    real(dp), allocatable, dimension(:,:)   :: vals

    if (iargc() /= 7) then
       write(*,*) '    Compute binned spectrum'
       write(*,*) '    Options: [infile] [binfile] [lmin_col] [lmax_col] [burnin] [outfile]'
       stop
    end if

    call getarg(2, infile)
    call getarg(3, binfile)
    call getarg(4, temp)
    read(temp,*) lmin_col
    call getarg(5, temp)
    read(temp,*) lmax_col
    call getarg(6, temp)
    read(temp, *) numburnin
    call getarg(7, outfile)


    call read_resfile(infile, unit, nspec, lmax, numchain, numiter, cls_pointer)
    open(58,file=trim(binfile))
    open(59,file=trim(outfile), recl=1024)
    do while (.true.)
       read(58,'(A)',end=97) line
       if (line(1:1) == '#' .or. trim(line)=='') cycle
       read(line,*) bins
       if (bins(lmax_col) > lmax) exit
       if (count(cls_pointer(bins(lmax_col),1,:,numburnin+2:numiter+1) /= 0.d0) == 0) cycle
       if (all(cls_pointer(bins(lmax_col),1,:,numburnin+2:numiter+1) == 0.d0)) cycle
       allocate(vals(count(cls_pointer(bins(lmax_col),1,:,numburnin+2:numiter+1) /= 0.d0),nspec))
       do s = 1, nspec
          k = 1
          do i = 1, numchain
             do j = 1+numburnin, numiter
                if (cls_pointer(bins(lmin_col),1,i,j+1) /= 0.d0) then
                   vals(k,s) = sum(cls_pointer(bins(lmin_col):bins(lmax_col),&
                        & 1,i,j+1)) / size(cls_pointer(bins(lmin_col):bins(lmax_col),s,i,j+1))
                   k = k+1
                end if
             end do
          end do
       end do
       write(59,fmt='(f6.1)',advance='no') (bins(lmin_col)+bins(lmax_col))/2.d0
       do s = 1, nspec
          write(59,fmt='(2f12.4)',advance='no') mean(vals(:,s)), sqrt(variance(vals(:,s)))
       end do
       write(59,*)
       deallocate(vals)
    end do
97  close(58)
    close(59)

    deallocate(cls_pointer)

  end subroutine cl2ascii_binned

  subroutine coadd_mapcov
    implicit none

    character(len=256) :: outfile, temp, infile, prefix, beamfile, maskfile, operation, mapfile1, mapfile2, Ncovfile1, Ncovfile2
    character(len=4)   :: chain_text
    character(len=5)   :: pix_text
    character(len=1)   :: s_text(3) = ['T', 'Q', 'U']
    integer(i4b)       :: i, j, k, l, m, n, p, unit, nreal, ordering, chain, lmax, counter, nval
    integer(i4b)       :: nside, npix, nmaps, comp, ncomp, burnin, numiter, numchain, seed, numcomp, nmax
    real(dp)           :: sigma_reg(3)
    logical(lgt)       :: compute_covar
    type(planck_rng)   :: handle
    integer(i4b), allocatable, dimension(:)        :: ind, n_per_chain
    integer(i4b), allocatable, dimension(:,:)      :: mask2map
    real(dp),     allocatable, dimension(:)        :: map_1D, mean_1D, W
    real(dp),     allocatable, dimension(:,:)      :: mean, cov, map, beam, alms_mean, map_1d_, mask, map1, map2, cov1, cov2
    real(dp),     pointer,     dimension(:,:)      :: pixwin
    real(sp),     allocatable, dimension(:,:)      :: alms
    complex(dpc), allocatable, dimension(:,:,:)    :: alms_cmplx
    real(sp),     allocatable, dimension(:,:,:)    :: fg_amp

    if (iargc() < 7) then
       write(*,*) '    Coadd two maps+covmat'
       write(*,*) '    Options:  [halfsum/halfdiff] [map1] [Ncov1] [map2] [Ncov2] [out prefix]'
       stop
    end if

    call getarg(2,operation)
    call getarg(3,mapfile1)
    call getarg(4,Ncovfile1)
    call getarg(5,mapfile2)
    call getarg(6,Ncovfile2)
    call getarg(7,prefix)

    i = getsize_fits(trim(mapfile1), nside=nside, nmaps=nmaps, ordering=ordering)
    npix = 12*nside**2

    ! Read maskfile
    allocate(map1(0:12*nside**2-1,nmaps), map2(0:12*nside**2-1,nmaps))
    call read_map(mapfile1, map1)
    call read_map(mapfile2, map2)

    open(unit,file=trim(Ncovfile1),form='unformatted')
    read(unit) n
    allocate(cov1(n,n), cov2(n,n))
    read(unit) ordering ! Ordering
    read(unit) nmaps ! Polarization status
    do i = 1, n
       read(unit) cov1(:,i)
    end do
    close(unit)

    open(unit,file=trim(Ncovfile2),form='unformatted')
    read(unit) n
    read(unit) ordering ! Ordering
    read(unit) nmaps ! Polarization status
    do i = 1, n
       read(unit) cov2(:,i)
    end do
    close(unit)

    if (trim(operation) == 'halfsum') then
       map1 = 0.5*(map1+map2)
       outfile = trim(prefix) // '_halfsum.fits'
    else 
       map1 = 0.5*(map1-map2)
       outfile = trim(prefix) // '_halfdiff.fits'
    end if
    cov1 = 0.25d0*(cov1+cov2)

    ! Output to standard formats
    call write_map3(outfile, map1, ordering=ordering)    

    outfile = trim(prefix) // '_Nsum.unf'
    open(unit,file=trim(outfile),form='unformatted')
    write(unit) n
    write(unit) ordering ! Ordering
    write(unit) nmaps ! Polarization status
    do i = 1, n
       write(unit) cov1(:,i)
    end do
    write(unit) .false. ! Not inverse
    close(unit)
       

  end subroutine coadd_mapcov



end program comm_process_resfiles
