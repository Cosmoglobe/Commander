module map_editor
  use healpix_types
  use pix_tools
  use alm_tools
  use fitstools
  use rngmod
  use map_editor_utils
  use math_tools
  use spline_1D_mod
  use spline_2D_mod
  use comm_utils
  use map_editor_simple_ops_mod
  implicit none

  real(dp), allocatable, dimension(:), private :: d, nu
  ! real(dp), parameter, private :: k_B      = 1.3806503d-23
  ! real(dp), parameter, private :: h        = 1.0545726691251021d-34 * 2.d0*pi !6.626068d-34
  ! real(dp), parameter, private :: c        = 2.99792458d8
  ! real(dp),            private :: T_CMB    = 2.7255d0

contains

  subroutine fix_monopole(mapname_in1, mapname_in2, maskfile, mapname_out)
    implicit none
    character(len=*) :: mapname_in1, mapname_in2, maskfile, mapname_out
    integer(i4b)     :: i, nside, order1, nmaps, order2, order3, npix
    real(dp)         :: mu, nullval
    logical(lgt)     :: anynull
    real(dp), allocatable, dimension(:,:) :: map1, map2, mask
    character(len=80), dimension(180)         :: header

    ! Read input map
    i = getsize_fits(mapname_in1, nside=nside, ordering=order1, nmaps=nmaps)
    i = getsize_fits(mapname_in2, ordering=order2)
    i = getsize_fits(maskfile,    ordering=order3)
    npix = nside2npix(nside)
    allocate(map1(0:npix-1,nmaps), map2(0:npix-1,nmaps), mask(0:npix-1,1))
    call read_bintab(mapname_in1, map1, npix, nmaps, nullval, anynull, header=header)
    call read_bintab(mapname_in2, map2, npix, nmaps, nullval, anynull)
    call read_bintab(maskfile,    mask, npix, 1, nullval, anynull)

    if (order2 /= order1) then
       do i = 1, nmaps
          if (order2 == 1) then
             call convert_ring2nest(nside, map2(:,i))
          else
             call convert_nest2ring(nside, map2(:,i))
          end if
       end do
    end if

    if (order3 /= order1) then
       do i = 1, 1
          if (order3 == 1) then
             call convert_ring2nest(nside, mask(:,i))
          else
             call convert_nest2ring(nside, mask(:,i))
          end if
       end do
    end if

    mu        = sum(mask(:,1)*(map1(:,1)-map2(:,1))) / sum(mask(:,1))
    write(*,*) 'mu = ', mu
    map1(:,1) = map1(:,1) - mu

    call write_result_map(trim(mapname_out), nside, order1, header, map1)

  end subroutine fix_monopole

  subroutine fit_line
    implicit none

    character(len=512) :: filename, string
    integer(i4b) :: i, j, k, n, setsize, nset
    real(dp) :: sigma_x, sigma_y
    real(dp), allocatable, dimension(:,:) :: data, data_red
    real(dp), allocatable, dimension(:)   :: a

    real(dp) :: C, b, D, V1, V2, sigma_a
    type(planck_rng) :: handle
    
    if (iargc() /= 5) then
       write(*,*) 'Usage: map_editor fit_line_to_ASCII_data [filename] [numpt] [sample size] [num sets]'
       stop
    end if
    
    call getarg(2, filename)
    call getarg(3, string)
    read(string,*) n
    call getarg(4, string)
    read(string,*) setsize
    call getarg(5, string)
    read(string,*) nset

    call rand_init(handle, 182941)

    allocate(data(n,2), data_red(setsize,2), a(nset))
    open(58,file=trim(filename))
    do i = 1, n
       read(58,*) data(i,1), data(i,2)
    end do
    close(58)

    do i = 1, nset
       do j = 1, setsize
          k = max(min(int(rand_uni(handle)*n)+1,n),1)
          data_red(j,:) = data(k,:)
       end do
       C  = mean(data_red(:,1)*data_red(:,2)) - mean(data_red(:,1))*mean(data_red(:,2))
       V2 = mean(data_red(:,2)**2)-mean(data_red(:,2))**2
       V1 = mean(data_red(:,1)**2)-mean(data_red(:,1))**2
       D  = (V1-V2)/(2.d0*C)
       a(i)  = D + C/abs(C) * sqrt(1-D**2)
    end do

    write(*,*) 'Slope = ', mean(a), '+/-', sqrt(variance(a))

    
  end subroutine fit_line
  
  subroutine print_scaled_gal_avg
    implicit none

    logical(lgt) :: anynull
    real(dp)     :: nullval, T0, beta0, x, alpha, alpha_min, alpha_max, dalpha, dl, db, theta, phi
    integer(i4b) :: numband, nside, npix, ordering, nmaps, i, j, k, ierr, n, comp, m, b
    integer(i4b) :: listpix(0:1000000), nlist
    real(dp), allocatable, dimension(:,:) :: maps, par
    real(dp), allocatable, dimension(:,:) :: map1, map2, diff
    real(dp), allocatable, dimension(:)   :: tot, numpt
    character(len=512) :: infile1, infile2, outprefix, temp
    character(len=10)  :: alpha_text
    character(len=80), dimension(180)         :: header

    call getarg(2,infile1)
    call getarg(3,infile2)
    call getarg(4,outprefix)

    ! Read input map
    i = getsize_fits(infile1, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map1(0:npix-1,nmaps), map2(0:npix-1,nmaps), diff(0:npix-1,nmaps))
    call read_bintab(infile1, map1, npix, nmaps, nullval, anynull)
    call read_bintab(infile2, map2, npix, nmaps, nullval, anynull)

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map1(:,i))
          call convert_nest2ring(nside, map2(:,i))
       end do
    end if
    
    comp      = 1
    alpha_min = 0.88d0
    alpha_max = 1.12d0
    n         = 9
    dalpha    = (alpha_max-alpha_min)/(n-1)
    dl        = 1.d0*pi/180.d0 ! In degrees
    m         = 360.d0 / dl
    db        = 1.5d0*pi/180.d0
    
    allocate(tot(-180:179), numpt(-180:179))
    open(58,file=trim(outprefix)//'_mean.dat')
    do i = 1, n
       alpha = alpha_min + (i-1) * dalpha
       diff = 0.d0
       where (map1 /= -1.6375d30 .and. map2 /= -1.6375d30)
          diff  = map1 - alpha*map2
       end where

       write(alpha_text,fmt='(f4.2)') alpha
       call write_minimal_header(header, 'MAP', nside=nside, order=1, polar=.true.)
       call write_result_map(trim(outprefix)//'_'//trim(alpha_text)//'.fits',   &
            & nside, 1, header, diff)

       call query_strip(nside, 0.5d0*pi-db, 0.5d0*pi+db, listpix, nlist)

       tot   = 0.d0
       numpt = 0.d0
       do j = 0, nlist-1
          if (diff(listpix(j),comp) /= 0.d0) then
             call pix2ang_ring(nside, listpix(j), theta, phi)
             b = int((phi+0.5d0*dl)/dl)
             if (b > 179) b = b-360
             tot(b)   = tot(b)   + diff(listpix(j),comp)
             numpt(b) = numpt(b) + 1
          end if
       end do

       open(59,file=trim(outprefix)//'_'//trim(alpha_text)//'.dat')
       do b = -180, 179
          if (numpt(b) > 0) then
             write(59,*) b, tot(b)/numpt(b)
          end if
       end do
       close(59)
       
    end do
    close(58)
    

  end subroutine print_scaled_gal_avg
  


  subroutine fit_gain_offset
    implicit none

    logical(lgt) :: anynull
    real(dp)     :: nullval, T0, beta0, x, a, b
    real(dp)     :: missval = -1.6375d30
    integer(i4b) :: numband, nside, npix, ordering, nmaps, i, j, k, ierr, n, comp, m, order_mask
    integer(i4b) :: listpix(0:1000000), nlist, maxiter, sum1, sum2
    real(dp), allocatable, dimension(:,:) :: maps, par
    real(dp), allocatable, dimension(:,:) :: map1, map2, mask, gainmask
    real(dp), allocatable, dimension(:)   :: tot, numpt, gain, offset
    character(len=512) :: infile1, infile2, outprefix, temp, maskfile, gainmaskfile
    character(len=10)  :: alpha_text
    character(len=80), dimension(180)         :: header

    if (iargc() /= 6) then
       write(*,*) 'Usage: map_editor fit_gain_offset [map1] [map2] [gain mask] [offset mask] [output residual map]'
       stop
    end if

    call getarg(2,infile1)
    call getarg(3,infile2)
    call getarg(4,gainmaskfile)
    call getarg(5,maskfile)
    call getarg(6,outprefix)

    ! Read input map
    i = getsize_fits(infile1, nside=nside, ordering=ordering, nmaps=nmaps)
    npix  = nside2npix(nside)
    nmaps = 1
    allocate(map1(0:npix-1,nmaps), map2(0:npix-1,nmaps), mask(0:npix-1,nmaps), gainmask(0:npix-1,nmaps))
    call read_bintab(infile1, map1, npix, nmaps, nullval, anynull)
    i = getsize_fits(infile2, nside=nside, ordering=order_mask)
    call read_bintab(infile2, map2, npix, nmaps, nullval, anynull)
    if (ordering /= order_mask) then
       if (order_mask == 1) then
          call convert_ring2nest(nside, map2(:,1))
       else
          call convert_nest2ring(nside, map2(:,1))
       end if
    end if
    i = getsize_fits(maskfile, ordering=order_mask)
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    if (ordering /= order_mask) then
       if (order_mask == 1) then
          call convert_ring2nest(nside, mask(:,1))
       else
          call convert_nest2ring(nside, mask(:,1))
       end if
    end if
    i = getsize_fits(gainmaskfile, ordering=order_mask)
    call read_bintab(gainmaskfile, gainmask, npix, nmaps, nullval, anynull)
    if (ordering /= order_mask) then
       if (order_mask == 1) then
          call convert_ring2nest(nside, gainmask(:,1))
       else
          call convert_nest2ring(nside, gainmask(:,1))
       end if
    end if

    where(abs((map1-missval)/missval) < 1.d-8 .or. abs((map2-missval)/missval) < 1.d-8 .or. &
          abs((mask-missval)/missval) < 1.d-8 .or. abs((gainmask-missval)/missval) < 1.d-8)
       map1 = 0.d0
       map2 = 0.d0
       mask = 0.d0
       gainmask = 0.d0
    end where

    where(abs(map1+1.6375d30) < 1d25 .or. map1 == 0.d0 .or. map1 > 1.d30)
       map1 = 0.d0
    end where

    where(abs(map2+1.6375d30) < 1.d25 .or. abs(map1+1.6375d30) < 1.d25)
       mask = 0.d0
    end where

    maxiter = 10000
    allocate(gain(0:maxiter), offset(0:maxiter))

    gain(0)   = 1.d0
    offset(0) = 0.d0

    do i = 1, maxiter

       ! Estimate offset of masked sky
       offset(i) = sum(mask(:,1)*(map2(:,1) - gain(i-1)*map1(:,1)))/sum(mask(:,1))

       ! Estimate gain over full sky
       !gain(i) = sum(map1(:,1)*(map2(:,1)-offset(i)))/sum(map1(:,1)**2)
       gain(i)   = sum(gainmask(:,1)*(map1(:,1)*(map2(:,1)-offset(i))))/sum(gainmask(:,1)*(map1(:,1)**2))

       if (i > 100) then
          exit
       endif

       if (i > 10) then
          if (abs(gain(i)-gain(i-10))/abs(gain(i)) < 1d-6 .and. abs(offset(i)-offset(i-10)) < 1d-6) exit
       end if
    end do

    write(*,fmt='(a,f8.4,a,f10.4)') 'gain = ', gain(i), ', offset = ', offset(i)
    write(*,*) gain(i), offset(i)
    map2 = (gain(i)*map1+offset(i)) - map2
    !map2 = map2 - map1

    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=.false.)
    call write_result_map(trim(outprefix), nside, ordering, header, map2(:,1:1))

  end subroutine fit_gain_offset

  subroutine fit_gain_offset_dipole
    implicit none

    logical(lgt) :: anynull
    real(dp)     :: nullval, T0, beta0, x, A(0:3,0:3), b(0:3), vector(3)
    integer(i4b) :: numband, nside, npix, ordering, nmaps, i, j, k, ierr, n, comp, m, order_mask
    integer(i4b) :: listpix(0:1000000), nlist, maxiter
    real(dp), allocatable, dimension(:,:) :: maps, par, offsetmap
    real(dp), allocatable, dimension(:,:) :: map1, map2, mask, gainmask, offset, harmonics, harmonics2
    real(dp), allocatable, dimension(:)   :: tot, numpt, gain
    character(len=512) :: infile1, infile2, outprefix, temp, maskfile, gainmaskfile
    character(len=10)  :: alpha_text
    character(len=80), dimension(180)         :: header

    if (iargc() /= 6) then
       write(*,*) 'Usage: map_editor fit_gain_offset_dipole [xmap] [ymap] [gain mask] [offset/dipole mask] [output residualfile]'
       stop
    end if

    call getarg(2,infile1)
    call getarg(3,infile2)
    call getarg(4,gainmaskfile)
    call getarg(5,maskfile)
    call getarg(6,outprefix)

    ! Read input map
    i = getsize_fits(infile1, nside=nside, ordering=ordering, nmaps=nmaps)
    npix  = nside2npix(nside)
    nmaps = 1
    allocate(map1(0:npix-1,nmaps), map2(0:npix-1,nmaps), mask(0:npix-1,nmaps), gainmask(0:npix-1,nmaps))
    call read_bintab(infile1, map1, npix, nmaps, nullval, anynull)
    i = getsize_fits(infile2, nside=nside, ordering=order_mask)
    call read_bintab(infile2, map2, npix, nmaps, nullval, anynull)
    if (ordering /= order_mask) then
       if (order_mask == 1) then
          call convert_ring2nest(nside, map2(:,1))
       else
          call convert_nest2ring(nside, map2(:,1))
       end if
    end if
    i = getsize_fits(maskfile, ordering=order_mask)
    call read_bintab(maskfile, mask, npix, 1, nullval, anynull)
    if (ordering /= order_mask) then
       if (order_mask == 1) then
          call convert_ring2nest(nside, mask(:,1))
       else
          call convert_nest2ring(nside, mask(:,1))
       end if
    end if
    i = getsize_fits(gainmaskfile, ordering=order_mask)
    call read_bintab(gainmaskfile, gainmask, npix, nmaps, nullval, anynull)
    if (ordering /= order_mask) then
       if (order_mask == 1) then
          call convert_ring2nest(nside, gainmask(:,1))
       else
          call convert_nest2ring(nside, gainmask(:,1))
       end if
    end if

    where(abs(map2+1.6375d30) < 1d25 .or. map2 == 0.d0 .or. map2 > 1.d30)
       map2 = 0.d0
    end where

    allocate(harmonics(0:npix-1,0:3))
    allocate(harmonics2(0:npix-1,0:3))
    do i = 0, npix-1
       if (ordering == 1) then
          call pix2vec_ring(nside, i, vector)
       else
          call pix2vec_nest(nside, i, vector)          
       end if
          
       if (mask(i,1) == 1.) then
          harmonics(i,0) = 1.d0
          harmonics(i,1) = vector(1)
          harmonics(i,2) = vector(2)
          harmonics(i,3) = vector(3)
       else
          harmonics(i,:) = 0.d0
       end if

       harmonics2(i,0) = 1.d0
       harmonics2(i,1) = vector(1)
       harmonics2(i,2) = vector(2)
       harmonics2(i,3) = vector(3)

    end do

    A = 0.d0
    do j = 0, 3
       do k = 0, 3
          A(j,k) = sum(harmonics(:,j) * harmonics(:,k))
       end do
    end do


    maxiter = 10000
    allocate(gain(0:maxiter), offset(0:3,0:maxiter), offsetmap(0:npix-1,1))

    gain(0)     = 1.d0
    offset(:,0) = 0.d0
    
    do i = 1, maxiter

       ! Estimate offset of masked sky
       b = 0.d0
       do j = 0, 3
          b(j) = sum((map2(:,1)-gain(i-1)*map1(:,1)) * harmonics(:,j))
       end do
       call solve_linear_system(A, offset(:,i), b)
       !offset(i) = sum(mask(:,1)*(map2(:,1) - gain(i-1)*map1(:,1)))/sum(mask(:,1))
       offsetmap = offset(0,i)
       do j = 1, 3
          offsetmap(:,1) = offsetmap(:,1) + harmonics2(:,j) * offset(j,i)
       end do

       ! Estimate gain over full sky
       !gain(i) = sum(map1(:,1)*(map2(:,1)-offset(i)))/sum(map1(:,1)**2)
       gain(i) = sum(gainmask(:,1)*(map1(:,1)*(map2(:,1)-offsetmap(:,1))))/sum(gainmask(:,1)*map1(:,1)**2)
       
       write(*,fmt='(a,i6,a,f40.12,a,4f40.12)') 'iter = ', i, ', gain = ', gain(i), ', offset = ', offset(:,i)
       !write(*,*,advance='no') 'iter = ', i, ', gain = ', gain(i), ', offset = ', offset(:,i)
       !write(*,*)
       if (i > 10) then
          if (abs(gain(i)-gain(i-10))/abs(gain(i)) < 1d-6 .and. abs(offset(0,i)-offset(0,i-10)) < 1d-6) exit
       end if
    end do

    write(*,fmt='(a,f40.12,a,4f40.12)') 'gain = ', gain(i), ', offset = ', real(offset(:,i),sp)
    !write(*,*,advance='no') 'gain = ', gain(i), ', offset = ', real(offset(:,i),sp)
    !write(*,*)
    !write(*,*) gain(i), offset(:,i)
    map2 = (gain(i)*map1+offsetmap) - map2
    !map2 = map2 - map1

    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=.false.)
    call write_result_map(trim(outprefix), nside, ordering, header, map2(:,1:1))

  end subroutine fit_gain_offset_dipole


  subroutine fit_gain
    implicit none

    logical(lgt) :: anynull
    real(dp)     :: nullval, T0, beta0, x, a, b, gain
    integer(i4b) :: numband, nside, npix, ordering, nmaps, i, j, k, ierr, n, comp, m, order_mask
    integer(i4b) :: listpix(0:1000000), nlist, maxiter
    real(dp), allocatable, dimension(:,:) :: maps, par
    real(dp), allocatable, dimension(:,:) :: map1, map2, mask, gainmask
    real(dp), allocatable, dimension(:)   :: tot, numpt
    character(len=512) :: infile1, infile2, outprefix, temp, maskfile, gainmaskfile
    character(len=10)  :: alpha_text
    character(len=80), dimension(180)         :: header

    if (iargc() /= 5) then
       write(*,*) 'Usage: map_editor fit_gain_offset [xmap] [ymap] [mask] [output residualfile]'
       stop
    end if

    call getarg(2,infile1)
    call getarg(3,infile2)
    call getarg(4,maskfile)
    call getarg(5,outprefix)

    ! Read input map
    i = getsize_fits(infile1, nside=nside, ordering=ordering, nmaps=nmaps)
    npix  = nside2npix(nside)
    nmaps = 3
    allocate(map1(0:npix-1,nmaps), map2(0:npix-1,nmaps), mask(0:npix-1,nmaps), gainmask(0:npix-1,nmaps))
    call read_bintab(infile1, map1, npix, nmaps, nullval, anynull)
    i = getsize_fits(infile2, nside=nside, ordering=order_mask)
    call read_bintab(infile2, map2, npix, nmaps, nullval, anynull)
    if (ordering /= order_mask) then
       do j = 1, nmaps
          if (order_mask == 1) then
             call convert_ring2nest(nside, map2(:,j))
          else
             call convert_nest2ring(nside, map2(:,j))
          end if
       end do
    end if
    i = getsize_fits(maskfile, ordering=order_mask)
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    if (ordering /= order_mask) then
       do j = 1, nmaps
          if (order_mask == 1) then
             call convert_ring2nest(nside, mask(:,j))
          else
             call convert_nest2ring(nside, mask(:,j))
          end if
       end do
    end if
    mask(:,1) = 0.d0

    where(abs(map2+1.6375d30) < 1d25 .or. map2 == 0.d0 .or. map2 > 1.d30)
       map2 = 0.d0
    end where

    ! Estimate gain over full sky
    gain = sum(mask*map1*map2)/sum(mask*map1**2)

    write(*,fmt='(a,f8.4)') 'gain = ', gain
    map2 = gain*map1 - map2

    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=.true.)
    call write_result_map(trim(outprefix), nside, ordering, header, map2)

  end subroutine fit_gain


  subroutine fit_ideal_dust
    implicit none

    logical(lgt) :: anynull
    real(dp)     :: nullval, T0, beta0, x
    integer(i4b) :: numband, nside, npix, ordering, nmaps, i, j, k, ierr
    real(dp), allocatable, dimension(:,:) :: maps, par
    real(dp), allocatable, dimension(:,:) :: cmb, fg
    character(len=512) :: infile1, infile2, outprefix, temp
    character(len=80), dimension(180)         :: header

    call getarg(2,infile1)
    call getarg(3,infile2)
    call getarg(4,outprefix)
    numband = iargc()-4
    T0      = 18.d0
    beta0   = 1.5d0

    if (iargc() < 3) then
       write(*,*) 'Usage: map_editor fit_ideal_dust [band1] [band2] [output prefix]'
       write(*,*) '  Support for additional bands can be added. This function fits a dust spectrum to input maps.     '
       write(*,*) '  Outputs amplitude, beta, T for the model, as well as CMB and chisq.' 
       stop
    end if

    ! Read input map
    i = getsize_fits(infile1, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(maps(0:npix-1,numband), par(0:npix-1,5), nu(numband), d(numband), cmb(0:npix-1,1), fg(0:npix-1,1))
    call read_bintab(infile1, fg, npix, nmaps, nullval, anynull)
    call read_bintab(infile2, cmb, npix, nmaps, nullval, anynull)

    x = h / (k_B*T0)
    do i = 1, numband
       call getarg(4+i,temp)
       read(temp,*) nu(i)
       nu(i) = nu(i) * 1d9
       maps(:,i) = fg(:,1) * &
            & (exp(x*nu(1))-1.d0) / (exp(x*nu(i))-1.d0) * (nu(i)/nu(1))**(beta0+1.d0) + & 
            & cmb(:,1) / ant2thermo(nu(i)/1d9)
    end do

    maps(:,1) = 1.00d0 * (maps(:,1) + 0.)
!    maps(:,4) = 1.01d0 * (maps(:,4) + 0.)

    ! Perform the fit
    do i = 0, npix-1
       if (mod(i,1000) == 0) write(*,*) i, npix
       d = maps(i,:)
       par(i,1) = fg(i,1)
       par(i,2) = beta0
       par(i,3) = T0
       par(i,4) = cmb(i,1)
       !call powell(par(i,1:4), chisq_ideal_dust, ierr)       
       par(i,5) = chisq_ideal_dust(par(i,1:3))
    end do

    ! Output result files
    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=.false.)
    call write_result_map(trim(outprefix)//'_amp.fits',   nside, ordering, header, par(:,1:1))
    call write_result_map(trim(outprefix)//'_beta.fits',  nside, ordering, header, par(:,2:2))
    call write_result_map(trim(outprefix)//'_T.fits',     nside, ordering, header, par(:,3:3))
    call write_result_map(trim(outprefix)//'_cmb.fits',     nside, ordering, header, par(:,4:4))
    call write_result_map(trim(outprefix)//'_chisq.fits', nside, ordering, header, par(:,5:5))

  end subroutine fit_ideal_dust

  function chisq_ideal_dust(p)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp)                           :: chisq_ideal_dust

    integer(i4b) :: i
    real(dp)     :: x, scale

    chisq_ideal_dust = 0.d0
    x                = h / (k_B*p(3))
    do i = 1, size(d)
       scale = 1.d0; if (i == 1) scale = 1.01d0
       chisq_ideal_dust = chisq_ideal_dust + &
            & (d(i) - &
            &   (scale*p(1)*(exp(x*nu(1))-1.d0) / (exp(x*nu(i))-1.d0) * (nu(i)/nu(1))**(p(2)+1.d0) + &
            &    p(4)/ ant2thermo(nu(i)/1d9)))**2
    end do

  end function chisq_ideal_dust



  subroutine compute_index_map
    implicit none

    character(len=512) :: mapname1, mapname2, outfile, tmp
    integer(i4b)       :: i, nside, nside2, ordering, ordering2, nmaps, nmaps2, npix, npix2
    real(dp)           :: nu1, nu2, nullval
    logical(lgt)       :: anynull
    real(dp), allocatable, dimension(:,:) :: map1, map2, ind
    character(len=80), dimension(180)         :: header

    if (iargc() < 6) then
       write(*,*) 'Usage: map_editor compute_spectral_index_map [map1] [nu1] [map2] [nu2] [index map] '
       stop
    end if

    call getarg(2,mapname1)
    call getarg(3,tmp)
    read(tmp,*) nu1
    call getarg(4,mapname2)
    call getarg(5,tmp)
    read(tmp,*) nu2
    call getarg(6,outfile)

    ! Read input map
    i = getsize_fits(mapname1, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map1(0:npix-1,nmaps), map2(0:npix-1,nmaps))
    call read_bintab(mapname1, map1, npix, nmaps, nullval, anynull)

    i = getsize_fits(mapname2, nside=nside2, ordering=ordering2, nmaps=nmaps2)
    npix2 = nside2npix(nside2)
    call read_bintab(mapname2, map2, npix2, nmaps2, nullval, anynull)
    if (nside /= nside2) then
       write(*,*) 'Error -- incompatible Nside. Exiting'
       stop
    end if
    write(*,*) ordering, ordering2
    if (ordering2 /= ordering) then
       write(*,*) 'a'
       do i = 1, nmaps2
          if (ordering2 == 1) then
             call convert_ring2nest(nside, map2(:,i))
          else
             call convert_nest2ring(nside, map2(:,i))
          end if
       end do
    end if

    nmaps = min(nmaps, nmaps2)
    allocate(ind(0:npix-1,nmaps))
    ind = -1.6375d30
    do i = 1, nmaps
       if (nu1 > 0) then
          where (map1(:,i)*map2(:,i) > 0.d0)
             ind(:,i) = log((map1(:,i)/ant2thermo(nu1)) / (map2(:,i)/ant2thermo(nu2))) / log(nu1/nu2)
          end where
       else
          where (map1(:,i)*map2(:,i) > 0.d0)
             ind(:,i) = log(map1(:,i) / map2(:,i)) / log(abs(nu1/nu2))
          end where
       end if
    end do

    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=(nmaps==3))
    call write_result_map(outfile, nside, ordering, header, ind)

  end subroutine compute_index_map

  subroutine make_procmask
    implicit none

    logical(lgt)       :: anynull
    integer(i4b)       :: i, j, k, n_neigh, n_valid, listpix(8), list_valid(8)
    integer(i4b)       :: nside, npix, nmaps, ordering
    character(len=512) :: partext, infile, outfile
    real(dp)           :: vec(3), cvec(3), threshold, nullval, missval = -1.6375e30
    real(dp)           :: sum_n, avg_n
    real(dp),     allocatable, dimension(:,:) :: map, mask, mask2
    character(len=80), dimension(180)         :: header


    if (iargc() < 4) then
       write(*,*) 'Usage: map_editor make_procmap [map] [threshold] [mask out] '
       write(*,*) ''       
       stop
    end if

    call getarg(2,infile)
    call getarg(3,partext)
    read(partext,*) threshold    
    call getarg(4,outfile)

    write(*,*) '   Creating procmask from: '//trim(infile)
    write(*,*) '   Threshold: '//trim(partext)
    write(*,*) '   Proc. mask name: '//trim(outfile)

    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps), mask(0:npix-1,nmaps))
    call read_bintab(infile, map, npix, nmaps, nullval, anynull)
    
    !make sure we have nested ordering
    if (ordering == 1) then
       do j = 1, nmaps
          call convert_ring2nest(nside, map(:,j))
       end do
    end if

    mask = 1.d0

    !Loop over maps
    do j = 1, nmaps
       !loop over pixels, take average of neighbouring pixels, mask if diff > threshold
       do i = 0,npix-1
          if (abs((map(i,j)-missval)/missval) < 1.d-6) then
             mask(i,j)=0.d0 !mask out missing pixels in inmap
          else
             ! get pixel numbers of neighbours
             call neighbours_nest(nside,i,listpix,n_neigh)

             !check validity of neighbours
             n_valid=0
             do k = 1,n_neigh
                if (abs((map(listpix(k),j)-missval)/missval) > 1.d-6) then
                   n_valid=n_valid+1
                   list_valid(n_valid)=listpix(k)
                end if
             end do

             !check if pixel should be masked out
             if (n_valid==0) then
                mask(i,j)=0.d0 !mask out lone pixel inside missing pixels
             else
                avg_n = sum(map(list_valid(:n_valid),j))/n_valid
                if (abs(map(i,j)-avg_n)>threshold) mask(i,j)=0.d0
             end if
          end if
       end do
       
    end do

    !if input map was ring ordered, transform mask back to ring
    if (ordering == 1) then 
       do j = 1, nmaps
          call convert_nest2ring(nside, mask(:,j))
       end do
    end if
        
    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=(nmaps==3))
    call write_result_map(trim(outfile), nside, ordering, header, mask)


  end subroutine make_procmask

  subroutine expand_mask_neighbours
    implicit none

    logical(lgt)       :: anynull
    integer(i4b)       :: i, j, k, n_neigh, n_valid, listpix(8), list_valid(8)
    integer(i4b)       :: nside, npix, nmaps, ordering, n_threshold
    character(len=512) :: partext, infile, outfile
    real(dp)           :: vec(3), cvec(3), nullval, missval = -1.6375e30
    real(dp)           :: sum_n, avg_n
    real(dp),     allocatable, dimension(:,:) :: map, mask, mask2
    character(len=80), dimension(180)         :: header


    if ((iargc() < 3) .or. (iargc() > 4) ) then
       write(*,*) 'Usage: map_editor expand_mask_neighbours [mask in] [mask out] [n_neigh (optional)] '
       write(*,*) ''       
       stop
    end if

    call getarg(2,infile)
    call getarg(3,outfile)
    if (iargc()==4) then
       call getarg(4,partext)
       read(partext,*) n_threshold    
    else
       n_threshold = 3
    end if

    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(mask(0:npix-1,nmaps),mask2(0:npix-1,nmaps))
    call read_bintab(infile, mask, npix, nmaps, nullval, anynull)
    
    !make sure we have nested ordering
    if (ordering == 1) then
       do j = 1, nmaps
          call convert_ring2nest(nside, mask(:,j))
       end do
    end if

    mask2(:,:) = mask(:,:)

    !Loop over maps
    do j = 1, nmaps
       do i = 0,npix-1
          if (mask2(i,j) > 0.5d0) then
             ! get pixel numbers of neighbours
             call neighbours_nest(nside,i,listpix,n_neigh)

             !check if >= n_threshold neighbors are 0
             n_valid=0
             do k = 1,n_neigh
                if (mask2(listpix(k),j) < 0.5d0) n_valid=n_valid+1
             end do
             if (n_valid >= n_threshold) then
                mask(i,j)=0.d0
             end if
          end if
       end do

    end do

    !if input map was ring ordered, transform mask back to ring
    if (ordering == 1) then 
       do j = 1, nmaps
          call convert_nest2ring(nside, mask(:,j))
       end do
    end if
        
    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=(nmaps==3))
    call write_result_map(trim(outfile), nside, ordering, header, mask)

  end subroutine expand_mask_neighbours

  subroutine convolve_masks
    implicit none

    integer(i4b)     :: i, j, k, n, npix, nside, nmaps, nmaps2, order, order2
    character(len=512) :: infile1, infile2, outfile
    real(dp)     :: nullval, missval = -1.6375e30
    logical(lgt) :: anynull
    real(dp),     allocatable, dimension(:,:) :: mask, mask2, mask_out
    character(len=80), dimension(180)         :: header

    if ((iargc() < 4) .or. (iargc() > 5) ) then
       write(*,*) 'Usage: map_editor convolve_masks [mask1] [mask2] [mask out] '
       write(*,*) ''       
       stop
    end if

    call getarg(2,infile1)
    call getarg(3,infile2)
    call getarg(4,outfile)

    i = getsize_fits(infile1, nside=nside, ordering=order, nmaps=nmaps)
    npix = nside2npix(nside)
    i = getsize_fits(infile2,ordering=order2)
    
    allocate(mask(0:npix-1,nmaps),mask2(0:npix-1,nmaps),mask_out(0:npix-1,nmaps))
    call read_bintab(infile1, mask, npix, nmaps, nullval, anynull)
    call read_bintab(infile1, mask2, npix, nmaps, nullval, anynull)

    if (order2 /= order) then
       do i = 1, nmaps
          if (order2 == 1) then
             call convert_ring2nest(nside, mask2(:,i))
          else
             call convert_nest2ring(nside, mask2(:,i))
          end if
       end do
    end if

    do j = 1,nmaps
       do i = 0, npix-1
          if (mask(i,j) == missval) mask(i,j) = 0.d0
          if (mask2(i,j) == missval) mask2(i,j) = 0.d0
          mask_out(i,j) = mask(i,j)*mask2(i,j)
       end do
    end do
        
    call write_minimal_header(header, 'MAP', nside=nside, order=order, polar=(nmaps==3))
    call write_result_map(trim(outfile), nside, order, header, mask_out)

  end subroutine convolve_masks

  subroutine make_ptsrc_map
    implicit none

    integer(i4b)       :: i, n, nside, npix, nmaps, unit, col, nlist, nsrc
    character(len=512) :: partext, convention, infile, outfile, beamfile
    character(len=1024) :: line
    real(dp)           :: vec(3), cvec(3), r, sigma, dOmega, Omega, fwhm, norm, unit_conv
    real(dp)           :: u(3), v(3), len_u, len_v, theta, flux, z, sgn
    integer(i4b), allocatable, dimension(:)   :: listpix
    real(dp),     allocatable, dimension(:,:) :: map
    real(dp),     allocatable, dimension(:)   :: data, src
    real(dp),     allocatable, dimension(:)   :: b_r, b_prof, b_prof2
    character(len=80), dimension(180)         :: header

    if (iargc() < 8) then
       write(*,*) 'Usage: map_editor make_ptsrc_map [ptsrc catalog] {WMAP,Planck} [flux column] '
       write(*,*) '            [nside] [fwhm] [Kcmb to MJy/sr] [outfile] (radial beam profile for WMAP)'
       stop
    end if
     
    call getarg(2,infile)
    call getarg(3,convention)
    call getarg(4,partext)
    read(partext,*) col
    call getarg(5,partext)
    read(partext,*) nside
    call getarg(6,partext)
    read(partext,*) fwhm
    call getarg(7,partext)
    read(partext,*) unit_conv
    call getarg(8,outfile)
    if (trim(convention) == 'WMAP') then
       call getarg(9,beamfile)
    end if
    npix   = 12*nside**2
    nmaps  = 3
    unit   = 58
    fwhm   = fwhm * pi/180.d0/60.
    sigma  = fwhm / sqrt(8.d0*log(2.d0))
    dOmega = 4.d0*pi / npix

    allocate(map(0:npix-1,nmaps), listpix(0:npix-1), src(0:npix-1))
    if (trim(convention) == 'WMAP') then
       allocate(data(col+2))
       open(unit,file=trim(beamfile))
       read(unit,*) n
       allocate(b_r(n), b_prof(n), b_prof2(n))
       do i = 1, n
          read(unit,*) b_r(i), b_prof(i)
       end do
       close(unit)
       b_r = b_r * pi/180.d0
!!$       open(58,file='test.dat')
!!$       do i = 1, n
!!$          write(58,*) b_r(i)*180.d0/pi, b_prof(i)
!!$       end do
       !call spline(b_r, b_prof, 0.d0, 0.d0, b_prof2)
       call smooth_spline('uniform', 1.d-12, b_r, b_prof, 0.d0, 0.d0, b_prof2)
!!$       write(58,*)
!!$       do i = 1, n
!!$          write(58,*) b_r(i)*180.d0/pi, b_prof(i)
!!$       end do
!!$       write(58,*)
!!$       do i = 0, 50000
!!$          write(58,*) i*1d-4, splint(b_r, b_prof, b_prof2, i*1.d-4*pi/180.d0)
!!$       end do
!!$       close(58)
!!$       stop
    else
       allocate(data(7))
    end if
    map = 0.d0
    open(unit,file=trim(infile))
    nsrc = 0
    do while (.true.)
       if (trim(convention) == 'WMAP') then
          read(unit,*,end=91) data
          flux = data(2+col)
          call ang2vec(0.5d0*pi-data(2)*pi/180.d0, data(1)*pi/180.d0, cvec)
          call query_disc(nside, cvec, 5*fwhm, listpix, nlist)       
          src = 0.d0
          do i = 0, nlist-1
             call pix2vec_ring(nside, listpix(i), vec)
             r      = acos(sum(vec*cvec))
             src(i) = max(splint(b_r, b_prof, b_prof2, r),0.d0)
          end do
          ! Listed values are peak flux density
          norm = max(data(2+col),0.d0) / (sum(src(0:nlist-1))*dOmega)
       else if (trim(convention) == 'Planck') then
          read(unit,*,end=91) data
          flux      = data(4) * 1.d-3 * 1.d-6
          data(1:2) = data(1:2) * pi/180.d0
          data(5:6) = data(5:6) * pi/180.d0/60.d0 / sqrt(8.d0*log(2.d0))
          data(7)   = data(7)   * pi/180.d0
          nsrc = nsrc+1
          if (mod(nsrc,100) == 0) write(*,*) 'Processing source no. ', nsrc
          call ang2vec(0.5d0*pi-data(2), data(1), cvec)
          call query_disc(nside, cvec, 10*fwhm, listpix, nlist)       
          src = 0.d0
          do i = 0, nlist-1
             call pix2vec_ring(nside, listpix(i), vec)
             r         = acos(sum(vec*cvec))
             if (data(5) /= data(6) .and. .false.) then
                z         = sum(cvec*vec)
                sgn    = 1.d0
                if (cvec(1)*vec(2)-cvec(2)*vec(1) < 0.d0) sgn = -1.d0
                u         = cvec(3) * cvec
                u(3)      = u(3) - 1.d0
                v         = vec - z * cvec
                len_u     = sqrt(sum(u*u))
                len_v     = sqrt(sum(v*v))
                theta     = sgn*acos(max(min((z * cvec(3) - vec(3)) / (len_u*len_v), 1.d0), -1.d0))
                theta     = theta + data(7)
             else
                theta = 0.d0
             end if
             src(i)    = exp(-0.5d0 * ((r*cos(theta)/data(5))**2 + (r*sin(theta)/data(6))**2))
          end do
       else
          write(*,*) 'Unsupported flux convention: ', trim(convention)
          stop
       end if
       ! Listed values are peak flux density
       norm           = max(flux,0.d0) / (sum(src(0:nlist-1))*dOmega)
       src(0:nlist-1) = norm / unit_conv * src(0:nlist-1)
       map(listpix(0:nlist-1),1) = map(listpix(0:nlist-1),1) + src(0:nlist-1)
    end do
91  close(unit)

    call write_minimal_header(header, 'MAP', nside=nside, order=1, polar=.true.)
    call write_result_map(outfile, nside, 1, header, map)
    deallocate(map)

  end subroutine make_ptsrc_map

  subroutine convert_beam(fwhm, lmax, nmaps, outfile, infile)
    implicit none

    integer(i4b),     intent(in) :: lmax, nmaps
    real(dp),         intent(in) :: fwhm
    character(len=*), intent(in) :: outfile
    character(len=*), intent(in), optional :: infile

    integer(i4b) :: nlheader, i, l
    real(dp),     allocatable, dimension(:)     :: b_in
    real(dp),     allocatable, dimension(:,:)   :: beam
    character(len=80), dimension(80) :: header

    allocate(beam(0:lmax, nmaps), b_in(nmaps))
    if (present(infile)) then
       if (infile(len(trim(infile))-3:len(trim(infile))) == 'fits') then
          call fits2cl(infile, beam, lmax, nmaps, header)
       else
          open(58,file=trim(outfile))
          beam = 0.d0
          do while (.true.)
             read(58,*,end=94) i, b_in
             beam(i,:) = b_in
          end do
94        close(58)          
       end if
    else
       call generate_beam(fwhm, lmax, beam)
    end if

    if (outfile(len(trim(outfile))-3:len(trim(outfile))) == 'fits') then
       header   = ''
       nlheader = 1
       call write_asctab(beam, lmax, nmaps, header, nlheader, outfile)
    else
       open(58,file=trim(outfile))
       do l = 0, lmax
          write(58,*) l, real(beam(l,:),sp)
       end do
       close(58)
    end if
    deallocate(beam)

  end subroutine convert_beam

  subroutine output_pointsource
    implicit none

    character(len=100) :: string_int
    integer(i4b)       :: nside, nmaps, npix, i, j, nlist
    real(dp)           :: lat, lon, fwhm_Q, ell_Q, ell_dir_Q, amp_Q, psi, beta, nu, nu0, aq, bq, cq
    real(dp)           :: fwhm_U, ell_U, ell_dir_U, au, bu, cu, amp_U, sigma_U
    real(dp)           :: vec(3), cvec(3), x, y, f, sigma_Q, sx_Q, sy_Q, sx_U, sy_U, cos2psi, sin2psi, ratio
    character(len=512) :: filename
    character(len=80), dimension(180) :: header

    real(dp),     allocatable, dimension(:,:) :: map
    integer(i4b), allocatable, dimension(:)   :: listpix

    call getarg(2,string_int)
    read(string_int,*) lon
    call getarg(3,string_int)
    read(string_int,*) lat
    call getarg(4,string_int)
    read(string_int,*) fwhm_Q
    call getarg(5,string_int)
    read(string_int,*) ell_Q     
    call getarg(6,string_int)
    read(string_int,*) ell_dir_Q
    call getarg(7,string_int)
    read(string_int,*) amp_Q
    call getarg(8,string_int)
    read(string_int,*) fwhm_U
    call getarg(9,string_int)
    read(string_int,*) ell_U     
    call getarg(10,string_int)
    read(string_int,*) ell_dir_U
    call getarg(11,string_int)
    read(string_int,*) amp_U
!    call getarg(12,string_int)
!    read(string_int,*) psi
!    call getarg(13,string_int)
!    read(string_int,*) beta
!    call getarg(14,string_int)
!    read(string_int,*) nu
!    call getarg(15,string_int)
!    read(string_int,*) nu0
    call getarg(12,string_int)
    read(string_int,*) nside
    call getarg(13,filename)

    lat       = lat * DEG2RAD
    lon       = lon * DEG2RAD
    ell_dir_Q = ell_dir_Q * DEG2RAD
    ell_dir_U = ell_dir_U * DEG2RAD
    psi       = psi * DEG2RAD
    sigma_Q   = fwhm_Q/sqrt(8.d0*log(2.d0))
    sigma_U   = fwhm_U/sqrt(8.d0*log(2.d0))
    sx_Q      = 2*sigma_Q / (2.d0-ell_Q)
    sx_U      = 2*sigma_U / (2.d0-ell_U)
    sy_Q      = 2*sigma_Q*(1.d0-ell_Q) / (2.d0-ell_Q)
    sy_U      = 2*sigma_U*(1.d0-ell_U) / (2.d0-ell_U)

    ! Allocate data structures
    npix  = 12*nside**2
    nmaps = 3
    allocate(map(0:npix-1,nmaps), listpix(0:npix-1))
    call ang2vec(0.5d0*pi-lat, lon, cvec)
    call query_disc(nside, cvec, 5*max(sigma_Q,sigma_U)*DEG2RAD, listpix, nlist, nest=1)

!    f     = (nu/nu0)**beta
!    ratio = (2.d0*pi*(0.88/sqrt(8*log(2.d0))*DEG2RAD)**2) / (2.d0*pi*(sigma*DEG2RAD)**2)
    
    map = 0.d0
    open(58,file='dist.dat')
    do i = 0, nlist-1
       call pix2vec_nest(nside, listpix(i), vec)
       call project2xy(vec, cvec, x, y)

       aq  = cos(ell_dir_Q)**2/(2*sx_Q**2) + sin(ell_dir_Q)**2/(2*sy_Q**2)
       bq  = - sin(2*ell_dir_Q)/(4*sx_Q**2) + sin(2*ell_dir_Q)/(4*sy_Q**2)
       cq  = sin(ell_dir_Q)**2/(2*sx_Q**2) + cos(ell_dir_Q)**2/(2*sy_Q**2)
       au  = cos(ell_dir_U)**2/(2*sx_U**2) + sin(ell_dir_U)**2/(2*sy_U**2)
       bu  = - sin(2*ell_dir_U)/(4*sx_U**2) + sin(2*ell_dir_U)/(4*sy_U**2)
       cu  = sin(ell_dir_U)**2/(2*sx_U**2) + cos(ell_dir_U)**2/(2*sy_U**2)

       map(listpix(i),2) = amp_Q * exp(-(aq*x**2 + 2*bq*x*y + cq*y**2))
       map(listpix(i),3) = amp_U * exp(-(au*x**2 + 2*bu*x*y + cu*y**2))

!       map(listpix(i),1) = amp_T * ratio * f * exp(-(a*x**2 + 2*b*x*y + c*y**2))
!       map(listpix(i),2) = amp_P * ratio * f * exp(-(a*x**2 + 2*b*x*y + c*y**2)) * cos(2.d0*psi)
!       map(listpix(i),3) = amp_P * ratio * f * exp(-(a*x**2 + 2*b*x*y + c*y**2)) * sin(2.d0*psi)

       !map(listpix(i),1) = amp_T * f * exp(-0.5d0 * (x**2 + y**2)/sigma**2)
       !map(listpix(i),2) = amp_P * f * exp(-0.5d0 * (x**2 + y**2)/sigma**2) * cos(2.d0*psi)
       !map(listpix(i),3) = amp_P * f * exp(-0.5d0 * (x**2 + y**2)/sigma**2) * sin(2.d0*psi)
       write(58,*) sqrt(x**2+y**2), map(listpix(i),2)
    end do
    close(58)
!    write(*,*) 'amp_T = ', amp_T
!    write(*,*) 'amp_P = ', amp_P
!    write(*,*) 'f     = ', f
!    write(*,*) maxval(map(:,1)), maxval(map(:,2:3))

    call write_minimal_header(header, 'MAP', nside=nside, ordering='NESTED', polar=nmaps==3)
    call write_result_map(filename, nside, 2, header, map)
    deallocate(map, listpix)

  end subroutine output_pointsource



  subroutine project2xy(mpixvec, centervec, x, y)
    implicit none

    real(dp)               :: sigma, radius, pixdist, northpixang, psi_true, x, y
    real(dp), dimension(3) :: centervec,  northvec, centerpixvec, centernorthvec, somevec, mpixvec, testvec
   
    psi_true = 0.d0
    northvec=[0,0,1]
    somevec=[centervec(2),-centervec(1),0.d0] 
    call crossproduct(somevec, centervec, centernorthvec)
    call crossproduct(centervec, mpixvec, somevec)
    call crossproduct(somevec, centervec, centerpixvec)
    call angdist(centernorthvec, centerpixvec, northpixang)
    call angdist(centervec, mpixvec, pixdist) 
    call crossproduct(centernorthvec, centerpixvec, testvec)
    if (sum(testvec*centervec) > 0) northpixang=-northpixang
    x = 2.d0*sin(pixdist/2.d0)*cos(psi_true-northpixang)*180.d0/pi
    y = 2.d0*sin(pixdist/2.d0)*sin(psi_true-northpixang)*180.d0/pi
     
  end subroutine project2xy

  !---------------------------------------------------------------------
  ! Crossproduct of two 3-vectors  
  !---------------------------------------------------------------------

  subroutine crossproduct(vector1, vector2, crossvector)
    implicit none

    real(dp), dimension(3), intent(in)  :: vector1, vector2
    real(dp), dimension(3), intent(out) :: crossvector

    crossvector=[vector1(2)*vector2(3)-vector1(3)*vector2(2), vector1(3)*vector2(1)-vector1(1)*vector2(3), &
                                                            & vector1(1)*vector2(2)-vector1(2)*vector2(1) ]

  end subroutine crossproduct

  ! Returns cos(2*psi) and sin(2*psi) of the rotation induced by
  ! parallel transport from vec1 to vec2.
  subroutine qu_transport_rot_v1(vec1, vec2, cos2psi, sin2psi)
    implicit none

    real(dp), dimension(3), intent(in)  :: vec1, vec2
    real(dp),               intent(out) :: cos2psi, sin2psi

    integer(i4b) :: i, j, nfield
    real(dp) :: len_u, len_v, z, cos_theta, sgn, c1, s1, c2, s2
    real(dp), dimension(3) :: u, v

    z = sum(vec1*vec2)
    if (abs(z) >= 1.d0-1.d-8) then
       cos2psi = 1; sin2psi = 0
       return
    end if

    sgn    = 1.d0
    if (vec1(1)*vec2(2)-vec1(2)*vec2(1) < 0.d0) sgn = -1.d0

    ! Rotation from vec1 to vec 2
    u         = vec1(3) * vec1 
    u(3)      = u(3) - 1.d0
    v         = vec2 - z * vec1
    len_u     = sqrt(sum(u*u))
    len_v     = sqrt(sum(v*v))
    ! Local angle from vertical
    cos_theta = max(min((z * vec1(3) - vec2(3)) / (len_u*len_v), 1.d0), -1.d0)
    ! Double it and calculate cos and sin
    c1 = 2*cos_theta**2-1
    s1 = 2*sgn*sqrt(1-cos_theta**2)*cos_theta

    ! Rotation from vec2 to vec 1; sgn is opposite from 1->2
    u          = vec2(3) * vec2
    u(3)       = u(3) - 1.d0
    v          = vec1 - z * vec2
    len_u      = sqrt(sum(u*u))
    len_v      = sqrt(sum(v*v))
    cos_theta  = max(min((z * vec2(3) - vec1(3)) / (len_u*len_v),1.d0),-1.d0)
    c2 =  2*cos_theta**2-1
    s2 = -2*sgn*sqrt(1-cos_theta**2)*cos_theta

    cos2psi = c1*c2+s1*s2
    sin2psi = c1*s2-s1*c2
  end subroutine


  subroutine shift_columns(filename, ncol, nside, ordering, header, map)
    implicit none

    character(len=*),                         intent(in)  :: filename
    integer(i4b),                             intent(in)  :: ncol
    integer(i4b),                             intent(out) :: nside, ordering
    real(dp),         pointer, dimension(:,:)             :: map
    character(len=80),         dimension(180)             :: header

    integer(i4b) :: j, nmaps, npix
    real(dp)     :: nullval
    logical(lgt) :: anynull

    j    = getsize_fits(filename, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = 12*nside**2

    allocate(map(0:npix-1,nmaps+ncol))
    map = 0.d0
    call read_bintab(filename, map(:,ncol+1:ncol+nmaps), npix, nmaps, nullval, anynull, header=header)

  end subroutine shift_columns

  subroutine compute_weighted_sum(infoname, nside, ordering, nmaps, outmap, header)
    implicit none

    character(len=*),                         intent(in)  :: infoname
    integer(i4b),                             intent(out) :: nside, ordering, nmaps
    real(dp),         pointer, dimension(:,:)             :: outmap
    character(len=80),         dimension(180), intent(out) :: header

    integer(i4b)       :: i, j, unit, npix, nummaps, nside_in, ordering_in, nmaps_in
    real(dp)           :: temp_r, weight, nullval, missval = -1.6375e30
    logical(lgt)       :: anynull
    character(len=512) :: filename

    real(dp), allocatable, dimension(:,:)   :: map

    unit = 24

    open(unit, file=trim(infoname))
    read(unit,*) nside
    read(unit,*) nummaps

    if (nummaps < 1) then
       write(*,*) 'Error: The number of input maps must be larger than 0.'
       stop
    end if

    ! Get parameters from the first map file, and check consistency
    read(unit,*) filename
    j = getsize_fits(filename, nside=nside, ordering=ordering, nmaps=nmaps)
    do i = 2, nummaps
       read(unit,*) filename
       j = getsize_fits(filename, nside=nside_in, ordering=ordering_in, nmaps=nmaps_in)
       if (nside_in /= nside) then
          write(*,*) 'Error: Nsides of input maps are not consistent'
          stop
       end if
       if (nmaps_in /= nmaps) then
          write(*,*) 'Error: Nmaps of input maps are not consistent'
          nmaps = 1
!          stop
       end if
    end do

    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    allocate(outmap(0:npix-1,nmaps))

    close(unit)
    open(unit, file=trim(infoname))
    read(unit,*) nside
    read(unit,*) nummaps

    outmap = 0.
    do i = 1, nummaps
       read(unit,*) filename, weight

       j = getsize_fits(filename, ordering=ordering_in)
       if (i == 1) then
          call read_bintab(filename, map, npix, nmaps, nullval, anynull, header=header)
       else
          call read_bintab(filename, map, npix, nmaps, nullval, anynull)
       end if

       if (ordering_in /= ordering) then
          if (ordering == 1) then
             do j = 1, nmaps
                call convert_nest2ring(nside, map(:,j))
             end do
          else
             do j = 1, nmaps
                call convert_ring2nest(nside, map(:,j))
             end do
          end if
       end if

       where (outmap /= missval .and. map /= missval) 
          outmap = outmap + weight * map
       elsewhere
          outmap = missval
       end where

    end do
    close(unit)

    deallocate(map)

  end subroutine compute_weighted_sum

  subroutine maskcount(maskfile)
    implicit none

    character(len=*),                   intent(in)    :: maskfile

    integer(i4b) :: npix, nside, ordering, nmaps, i
    real(dp)     :: nullval
    logical(lgt) :: anynull
    real(dp), allocatable, dimension(:,:) :: mask

    npix = getsize_fits(maskfile, nside=nside, ordering=ordering, nmaps=nmaps)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    do i=1,nmaps
       write(*,*) int(i,i2b), '% Unmasked =', count(mask(:,i) > 0.5d0)/real(npix,dp)*100.d0, &
            & '% Masked =', count(mask(:,i)<=0.5d0)/real(npix,dp)*100.d0
    end do
    write(*,*) ' Total % unmasked =', count(mask>0.5d0)/real(size(mask),dp)*100.d0, &
         & '% Masked =', count(mask<=0.5d0)/real(size(mask),dp)*100.d0
    write(*,*)
    do i=1,nmaps
       write(*,*) int(i,i2b), 'Number unmasked =', count(mask(:,i)>0.5d0), &
            & ' Number masked =', count(mask(:,i)<=0.5d0)
    end do
    write(*,*) int(i,i2b), 'Total number unmasked =', count(mask>0.5d0), &
         & ', Total number masked =', count(mask<=0.5d0)

  end subroutine maskcount

  subroutine zero_count
    implicit none

    integer(i4b)       :: i,j,count
    real(dp)           :: nullval,missval
    real(dp)           :: nullval_dp
    logical(lgt)       :: anynull
    character(len=512) :: infile
    integer(i4b)       :: ordering, nmaps

    integer(i4b) :: nside, npix

    real(dp),     allocatable, dimension(:,:)   :: inmap

    if (iargc() < 2) then
       write(*,*) 'Usage: map_editor zero_count [map_in]'
       stop
    end if
    missval=-1.6375e30
    call getarg(2,infile)

    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(inmap(0:npix-1,nmaps))
    call read_bintab(infile, inmap, npix, nmaps, nullval, anynull)
    write(*,*)
    do i =1,nmaps
       count=0
       do j=0,npix-1
          if (inmap(j,i)==0.0) then
             count = count + 1
          end if
       end do
       write(*,*)"  Map:", i,",   #pixels equal to zero:",count
    end do
    write(*,*)
    deallocate(inmap)
  end subroutine zero_count

  subroutine value_count
    implicit none

    integer(i4b)       :: i,j,count
    real(dp)           :: nullval,missval
    real(dp)           :: nullval_dp,val_comp,zero_val
    logical(lgt)       :: anynull
    character(len=512) :: infile,temptex
    integer(i4b)       :: ordering, nmaps

    integer(i4b) :: nside, npix

    real(dp),     allocatable, dimension(:,:)   :: inmap

    if (iargc() /= 3) then
       write(*,*) 'Usage: map_editor valuecount [map_in] [value]'
       stop
    end if
    missval=-1.6375e30
    call getarg(2,infile)
    call getarg(3,temptex)
    read(temptex,*) val_comp
    
    !covering the basics!
    temptex="0.0"
    read(temptex,*) zero_val
    if (zero_val==val_comp) then
       call zero_count !we call zero_count instead
       stop
    end if
    temptex="0"
    read(temptex,*) zero_val
    if (zero_val==val_comp) then
       call zero_count !we call zero_count instead
       stop
    end if
    temptex="0."
    read(temptex,*) zero_val
    if (zero_val==val_comp) then
       call zero_count !we call zero_count instead
       stop
    end if
    
    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(inmap(0:npix-1,nmaps))
    call read_bintab(infile, inmap, npix, nmaps, nullval, anynull)
    write(*,*)
    do i =1,nmaps
       count=0
       do j=0,npix-1
          if (abs((inmap(j,i)-val_comp)/val_comp) < 1d-6 ) then
             count = count + 1
          end if
       end do
       write(*,*)"  Map:", i,",   #pixels equal to value:",count
    end do
    write(*,*)
    deallocate(inmap)
  end subroutine value_count

  subroutine misspix_count
    implicit none

    integer(i4b)       :: i,j,count
    real(dp)           :: nullval,missval
    real(dp)           :: nullval_dp,val_comp
    logical(lgt)       :: anynull
    character(len=512) :: infile,temptex
    integer(i4b)       :: ordering, nmaps

    integer(i4b) :: nside, npix

    real(dp),     allocatable, dimension(:,:)   :: inmap

    if (iargc() /= 2) then
       write(*,*) 'Usage: map_editor misspixcount [map_in]'
       stop
    end if
    missval=-1.6375e30
    call getarg(2,infile)

    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(inmap(0:npix-1,nmaps))
    call read_bintab(infile, inmap, npix, nmaps, nullval, anynull)
    write(*,*)
    do i =1,nmaps
       count=0
       do j=0,npix-1
          if (abs((inmap(j,i)-missval)/missval) < 1d-5) then
             count = count + 1
          end if
       end do
       write(*,*)"  Map:", i,",   missing pixels:",count
    end do
    write(*,*)
    deallocate(inmap)
  end subroutine misspix_count

  subroutine badcount(mapfile)
    implicit none

    character(len=*),                   intent(in)    :: mapfile

    integer(i4b) :: npix, nside, ordering, nmaps, i
    real(dp)     :: nullval, missval = -1.6375e30
    logical(lgt) :: anynull
    real(dp), allocatable, dimension(:,:) :: mask

    npix = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(mapfile, mask, npix, nmaps, nullval, anynull)
    write(*,*)
    write(*,*) 'Map nr.     Bad pixels       Good pixels'
    do i=1,nmaps

       write(*,'(I5,A,I12,A,I12)') int(i,i2b) ,'      ', count(mask(:,i)<=missval), '      ', count(mask(:,i)/=missval)
    end do
    write(*,*)
    write(*,*) 'Total number bad =', count(mask<=missval), ', Total number good =', count(mask/=missval)
    write(*,*)
  end subroutine badcount


  subroutine polarization_fraction
    implicit none

    integer(i4b)       :: i,j,count
    real(dp)           :: nullval,missval
    real(dp)           :: nullval_dp,val_comp,zero_val
    logical(lgt)       :: anynull
    character(len=512) :: infile, temptex, outfile
    integer(i4b)       :: ordering, ordering2, nmaps, nmaps2

    integer(i4b) :: nside, nside2, npix

    real(dp),     allocatable, dimension(:,:)   :: polA, polB, mapT, fracmap, P_map
    character(len=80),         dimension(1:180) :: header

    if (iargc() < 5 .or. iargc() > 6 ) then
       write(*,*) 'Usage: map_editor polarization_fraction [T-map] [polA] [polB] [outmap] [P-map out (optional)]'
       stop
    end if
    missval=-1.6375e30
    
    call getarg(2,infile)
    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(mapT(0:npix-1,nmaps))
    call read_bintab(infile, mapT, npix, nmaps, nullval, anynull)
    
    call getarg(3,infile)
    i = getsize_fits(infile, nside=nside2, ordering=ordering2, nmaps=nmaps2)
    if (nside /= nside2) then
       write(*,*) 'Nside of T-map and polA are not equal. Exiting!'
       stop
    else if (nmaps2 < 3) then
       write(*,*) 'Nmaps of polA is less than 3. Exiting!'
       stop
    end if
    allocate(polA(0:npix-1,nmaps2))
    call read_bintab(infile, polA, npix, nmaps2, nullval, anynull)
    if (ordering /= ordering2) then
       do i = 1, nmaps2
          if (ordering2 == 1) then
             call convert_ring2nest(nside2, polA(:,i))
          else
             call convert_nest2ring(nside2, polA(:,i))
          end if
       end do
    end if
    
    call getarg(4,infile)
    i = getsize_fits(infile, nside=nside2, ordering=ordering2, nmaps=nmaps2)
    if (nside /= nside2) then
       write(*,*) 'Nside of T-map and polB are not equal. Exiting!'
       stop
    else if (nmaps2 < 3) then
       write(*,*) 'Nmaps of polB is less than 3. Exiting!'
       stop
    end if
    allocate(polB(0:npix-1,nmaps2))
    call read_bintab(infile, polB, npix, nmaps2, nullval, anynull)
    if (ordering /= ordering2) then
       do i = 1, nmaps2
          if (ordering2 == 1) then
             call convert_ring2nest(nside2, polB(:,i))
          else
             call convert_nest2ring(nside2, polB(:,i))
          end if
       end do
    end if

    call getarg(5,outfile)
    allocate(fracmap(0:npix-1,nmaps),P_map(0:npix-1,nmaps))

    P_map(:,1) = dsqrt(max(polA(:,2)*polB(:,2) + polA(:,3)*polB(:,3),0.0))
    fracmap(:,1) = P_map(:,1) / mapT(:,1)

    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=nmaps==3)
    call write_result_map(outfile, nside, ordering, header, fracmap(:,:nmaps))
    if (iargc() == 6) then
       call getarg(6,outfile)
       call write_result_map(outfile, nside, ordering, header, P_map(:,:nmaps))
    end if
    deallocate(mapT,polA,polB,fracmap,P_map)
  end subroutine polarization_fraction

  subroutine apply_mask1(maskfile, nside, ordering, map, fact)
    implicit none

    integer(i4b),                       intent(in)    :: nside, ordering
    character(len=*),                   intent(in)    :: maskfile
    real(dp),         dimension(0:,1:), intent(inout) :: map
    real(dp),                           intent(in), optional :: fact

    integer(i4b) :: i, j, nside_mask, ordering_mask, npix, n, nmaps, nmaps_mask
    real(dp)     :: nullval
    logical(lgt) :: anynull
    real(dp)     :: vec(3), radius
    integer(i4b), allocatable, dimension(:) :: listpix
    real(dp), allocatable, dimension(:,:) :: mask

    nmaps = size(map,2)

    ! Read rms file and check consistency
    npix = getsize_fits(maskfile, nside=nside_mask, ordering=ordering_mask, nmaps=nmaps_mask)

    if (nside_mask /= nside) then
       write(*,*) 'Error: Different Nsides of the two files.'
       stop
    end if
    if (nmaps_mask /= nmaps) then
       write(*,*) 'Error: Different Nmaps of the two files. Putting nmaps=1'
       write(*,*) 'Mask will be applied to all layers of map.'
       nmaps=1
    end if

    allocate(mask(0:npix-1,nmaps))
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    
    do i = 1, nmaps
       if (ordering_mask /= ordering) then
          if (ordering_mask == 1) call convert_ring2nest(nside, mask(:,i))
          if (ordering_mask == 2) call convert_nest2ring(nside, mask(:,i))
       end if
    end do

    if (present(fact)) then
       allocate(listpix(0:npix-1))
       radius = fact * pi/180.d0/60.
       do i = 0, npix-1
          if (mod(i,1000) == 0) write(*,*) i, npix
          if (ordering == 1) then
             call pix2vec_ring(nside, i, vec)
          else 
             call pix2vec_nest(nside, i, vec)            
          end if
          call query_disc(nside, vec, radius, listpix, n, nest=ordering-1) !query_disc uses ring by default
          do j = 1,nmaps
             if (any(mask(listpix(0:n-1),j) == 1.d0) .and. any(mask(listpix(0:n-1),j) == 0.d0)) then
                ! Border point
                if (nmaps > 1) then
                   map(listpix(0:n-1),j) = -1.6375d30
                else
                   map(listpix(0:n-1),:) = -1.6375d30
                end if
             end if
          end do
       end do
       deallocate(listpix)
    else
       do j = 1, nmaps
          do i = 0, npix-1
             if (mask(i,j) < 0.5d0) then
                if (nmaps > 1) then
                   map(i,j) = -1.6375d30
                else
                   map(i,:) = -1.6375d30
                end if
             end if
          end do
       end do
    end if

    ! TMR hack to mask gc
!    call ang2vec(pi/2.,0.d0,vec)
!    write(*,*) vec
!    allocate(listpix(0:npix-1))
!    call query_disc(nside, vec, 0.8d0*DEG2RAD, listpix, n, nest=1)
!    write(*,*) listpix(0:n-1)
!    write(*,*) 'Masking out gc pixels'
!    map(listpix(0:n-1),:) = -1.6375e30

  end subroutine apply_mask1


  subroutine output_mean_and_stddev(nside, ordering, nmaps, header, mu, rms)
    implicit none

    integer(i4b),                               intent(out)          :: nside, nmaps, ordering
    real(dp),           pointer, dimension(:,:)                      :: mu, rms
    character(len=80),           dimension(180)                      :: header

    integer(i4b)       :: npix
    integer(i4b)       :: i, j, l, m, numfiles
    real(dp)           :: nullval
    real(dp)           :: sigma_sq, nullval_dp
    logical(lgt)       :: anynull
    character(len=512) :: winfile_in, winfile_out
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    real(dp),     allocatable, dimension(:,:)   :: map
    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     pointer,     dimension(:,:)   :: weights
    real(dp),     allocatable, dimension(:,:)   :: beam
    character(len=512), allocatable, dimension(:) :: filenames
    real(dp),                  dimension(2)     :: zbounds = 0.d0

    numfiles = iargc()-3
    allocate(filenames(numfiles))
    do i = 1, numfiles
       call getarg(i+3, filenames(i))
    end do

    ! Read input map
    i = getsize_fits(filenames(1), nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps), mu(0:npix-1,nmaps), rms(0:npix-1,nmaps))
    
    ! Compute mu
    mu = 0.d0
    do i = 1, numfiles
       if (i == 1) then
          call read_bintab(filenames(i), map, npix, nmaps, nullval, anynull, header=header)
       else
          call read_bintab(filenames(i), map, npix, nmaps, nullval, anynull)          
       end if

       !write(*,*) i, mean(map(:,1)), sqrt(variance(map(:,1)))
       write(*,*) 'a', i, numfiles
       !call input_map(filenames(i), map, npix, nmaps)
       mu = mu + map
    end do
    mu = mu / numfiles

    ! Compute standard deviation
    rms = 0.d0
    do i = 1, numfiles
       call read_bintab(filenames(i), map, npix, nmaps, nullval, anynull)
       write(*,*) 'b', i, numfiles
!       call input_map(filenames(i), map, npix, nmaps)
       rms = rms + (map-mu)**2
    end do
    rms = sqrt(rms / (numfiles-1))

    deallocate(map)

  end subroutine output_mean_and_stddev

  subroutine smooth_map(infile, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
       & output_EB, beamfile_in, beamfile_out, fwhm_in, fwhm_out)
    implicit none

    character(len=512),                         intent(in)           :: infile
    real(dp),                                   intent(inout)        :: r_fill
    integer(i4b),                               intent(in)           :: nside, lmin, lmax
    real(dp),                                   intent(in), optional :: fwhm_in, fwhm_out
    character(len=512),                         intent(in), optional :: beamfile_in, beamfile_out
    integer(i4b),                               intent(out)          :: nmaps, ordering
    logical(lgt),                               intent(in)           :: output_EB
    real(dp),          pointer, dimension(:,:)                       :: map
    character(len=80),          dimension(180)                       :: header

    integer(i4b)       :: nside_in, npix, npix_in
    integer(i4b)       :: i, j, k, l, m, r, nlist
    real(dp)           :: nullval, tot
    real(dp)           :: sigma_sq, nullval_dp, vec(3), theta, phi, rf
    logical(lgt)       :: anynull
    character(len=512) :: winfile_in, winfile_out
    character(len=4)   :: ntext

    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     allocatable, dimension(:,:)   :: map_in, map_buffer
    real(dp),     pointer,     dimension(:,:)   :: pixwin_in, pixwin_out
    real(dp),     pointer,     dimension(:,:)   :: weights
    real(dp),     allocatable, dimension(:,:)   :: beam_in, beam_out
    real(dp),                  dimension(2)     :: zbounds = 0.d0
    integer(i4b), allocatable, dimension(:)     :: listpix

    ! Read input map
    i = getsize_fits(infile, nside=nside_in, ordering=ordering, nmaps=nmaps)
    if (nmaps /= 1 .and. nmaps /= 3) then
       if (nmaps < 3) then
          nmaps = 1
       else
          nmaps = 3
       end if
    end if

    npix    = nside2npix(nside)
    npix_in = nside2npix(nside_in)
    allocate(map_in(0:npix_in-1,nmaps), map_buffer(0:npix_in-1,nmaps))
    allocate(map(0:npix-1,nmaps))

    call read_bintab(infile, map_in, npix_in, nmaps, nullval, anynull, header=header)

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside_in, map_in(0:npix_in-1,i))
       end do
    end if

    ! Read pixel windows
    call read_pixwin(nside_in, nmaps, pixwin_in)
    call read_pixwin(nside, nmaps, pixwin_out)
    call read_ringweights(nside_in, weights)

    ! Create or read beams
    allocate(beam_in(0:lmax, nmaps))
    allocate(beam_out(0:lmax, nmaps))

    if (present(fwhm_in)) then
       call generate_beam(abs(fwhm_in), lmax, beam_in)
       if (fwhm_in < 0.) beam_in = 1.d0 / beam_in
    else
       call generate_beam(0.d0, lmax, beam_in, beamfile_in)
       !if (nmaps == 3) then
       !   beam_in(:,2) = beam_in(:,1)
       !   beam_in(:,3) = beam_in(:,1)
       !   write(*,*) 'Warning: Setting polarized beams equal to temperature beam'
       !end if
    end if

    if (present(fwhm_out)) then
       call generate_beam(abs(fwhm_out), lmax, beam_out)
    else
       call generate_beam(0.d0, lmax, beam_out, beamfile_out)
       !if (nmaps == 3) then
       !   beam_out(:,2) = beam_out(:,1)
       !   beam_out(:,3) = beam_out(:,1)
       !   write(*,*) 'Warning: Setting polarized beams equal to temperature beam'
       !end if
    end if

    allocate(listpix(0:npix_in-1))
    map_buffer = map_in
    do j = 1, nmaps
       do i = 0, npix_in-1
          if (abs(map_in(i,j)) > 1e30) then
             if (r_fill >= 0.d0) then
                call pix2vec_ring(nside_in, i, vec)
                do r = 1, r_fill
                   rf = (r+0.2d0) * sqrt(4.d0*pi/npix_in)
                   call query_disc(nside_in, vec, rf, listpix, nlist)
                   tot = 0.d0
                   m   = 0.d0
                   do k = 0, nlist-1
                      if (abs(map_buffer(listpix(k),j)) < 1.d30) then
                         tot = tot + map_buffer(listpix(k),j)
                         m   = m   + 1
                      end if
                   end do
                   if (m > 0) then
                      map_in(i,j) = tot / m
                      exit
                   end if
                end do

                if (r > r_fill) then
                   call pix2ang_ring(nside_in, i, theta, phi)
                   write(*,*) 'Error: Filter did not return valid pixel; increase r_fill'
                   write(*,*) '     p         = ', i
                   write(*,*) '     theta     = ', 90.d0-theta*180/pi
                   write(*,*) '     phi       = ', phi*180/pi
                   write(*,*) '     radius    = ', r_fill
                   write(*,*) '     neighbors = ', listpix(0:nlist-1)
                   write(*,*) '     vals      = ', map_buffer(listpix(0:nlist-1),j)
                   stop
                end if
             else
                map_in(i,j) = 0.d0
             end if
          end if
       end do
    end do
    deallocate(listpix, map_buffer)

    ! Compute the spherical harmonics transform
    allocate(alms(nmaps, 0:lmax, 0:lmax))

    alms = cmplx(0.,0.)
    if (nmaps == 1) then
       call map2alm(nside_in, lmax, lmax, map_in(:,1), alms, zbounds, weights)
    else
       call map2alm(nside_in, lmax, lmax, map_in, alms, zbounds, weights)
    end if

    ! Deconvolve old beam and pixel window, and convolve new beam and pixel window
    do i = 1, nmaps
       alms(i,0:lmin-1,:) = cmplx(0.,0.)
       do l = lmin, lmax
          if (abs(pixwin_in(l,i)*beam_in(l,i)) < 0.00000001d0) then
             if (l > 1 .or. i == 1) write(*,*) 'Multipole l = ', l, ' set to zero'
             alms(i,l,0:l) = cmplx(0.,0.)
          else
             alms(i,l,0:l) = alms(i,l,0:l) * &
                  & (pixwin_out(l,i)*beam_out(l,i)) / (pixwin_in(l,i)*beam_in(l,i))
          end if
       end do
    end do

    ! Compute the inverse spherical harmonics
    if (nmaps == 1) then
       call alm2map(nside, lmax, lmax, alms, map(:,1))
    else
       if (output_EB) then
          do i = 1, nmaps
             call alm2map(nside, lmax, lmax, alms(i:i,:,:), map(:,i))
          end do
       else
          call alm2map(nside, lmax, lmax, alms, map)       
       end if
    end if

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_ring2nest(nside, map(:,i))
       end do
    end if

    deallocate(weights)
    deallocate(pixwin_in)
    deallocate(pixwin_out)
    deallocate(beam_in)
    deallocate(beam_out)
    deallocate(map_in)
    deallocate(alms)

  end subroutine smooth_map


  subroutine smooth_map_zerospin(infile, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
       & output_EB, beamfile_in, beamfile_out, fwhm_in, fwhm_out)
    implicit none

    character(len=512),                         intent(in)           :: infile
    real(dp),                                   intent(inout)        :: r_fill
    integer(i4b),                               intent(in)           :: nside, lmin, lmax
    real(dp),                                   intent(in), optional :: fwhm_in, fwhm_out
    character(len=512),                         intent(in), optional :: beamfile_in, beamfile_out
    integer(i4b),                               intent(out)          :: nmaps, ordering
    logical(lgt),                               intent(in)           :: output_EB
    real(dp),           pointer, dimension(:,:)                      :: map
    character(len=80),           dimension(180)                      :: header

    integer(i4b)       :: nside_in, npix, npix_in
    integer(i4b)       :: i, j, k, l, m, r, nlist, map_count, nmaps2
    real(dp)           :: nullval, tot
    real(dp)           :: sigma_sq, nullval_dp, vec(3), theta, phi, rf
    logical(lgt)       :: anynull
    character(len=512) :: winfile_in, winfile_out
    character(len=4)   :: ntext

    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     allocatable, dimension(:,:)   :: map_in, map_buffer, map_in2
    real(dp),     pointer,     dimension(:,:)   :: pixwin_in, pixwin_out
    real(dp),     pointer,     dimension(:,:)   :: weights
    real(dp),     allocatable, dimension(:,:)   :: beam_in, beam_out, beam_in2, beam_out2
    real(dp),                  dimension(2)     :: zbounds = 0.d0
    integer(i4b), allocatable, dimension(:)     :: listpix

    ! Read input map
    i = getsize_fits(infile, nside=nside_in, ordering=ordering, nmaps=nmaps2)
    if (nmaps2 /= 1 .and. nmaps2 /= 3) then
       if (nmaps2 < 3) then
          nmaps2 = 1
       else
          nmaps2 = 3
       end if
    end if

    nmaps = 1

    npix    = nside2npix(nside)
    npix_in = nside2npix(nside_in)
    allocate(map_in2(0:npix_in-1,nmaps2))
    allocate(map(0:npix-1,nmaps2))
    call read_bintab(infile, map_in2, npix_in, nmaps2, nullval, anynull, header=header)

    ! Read pixel windows
    call read_pixwin(nside_in, nmaps, pixwin_in)
    call read_pixwin(nside, nmaps, pixwin_out)
    call read_ringweights(nside_in, weights)


    do map_count = 1, nmaps2 
       write(*,*) ""
       write(*,*) "map", map_count, "of", nmaps2
       write(*,*) ""

       allocate(map_buffer(0:npix_in-1,nmaps))
       allocate(map_in(0:npix_in-1,nmaps))
       map_in(:,1)=map_in2(:,map_count)
       if (ordering == 2) then
          do i = 1, nmaps
             call convert_nest2ring(nside_in, map_in(0:npix_in-1,i))
          end do
       end if

       ! Create or read beams
       allocate(beam_in2(0:lmax, nmaps2))
       allocate(beam_out2(0:lmax, nmaps2))

       ! Create or read beams
       allocate(beam_in(0:lmax, nmaps))
       allocate(beam_out(0:lmax, nmaps))

       if (present(fwhm_in)) then
          call generate_beam(abs(fwhm_in), lmax, beam_in)
          if (fwhm_in < 0.) beam_in = 1.d0 / beam_in
       else
          call generate_beam(0.d0, lmax, beam_in2, beamfile_in)
          beam_in(:,1)=beam_in2(:,map_count)
          !if (nmaps == 3) then
          !   beam_in(:,2) = beam_in(:,1)
          !   beam_in(:,3) = beam_in(:,1)
          !   write(*,*) 'Warning: Setting polarized beams equal to temperature beam'
          !end if
       end if

       if (present(fwhm_out)) then
          call generate_beam(abs(fwhm_out), lmax, beam_out)
       else
          call generate_beam(0.d0, lmax, beam_out2, beamfile_out)
          beam_out(:,1)=beam_out2(:,map_count)
          !if (nmaps == 3) then
          !   beam_out(:,2) = beam_out(:,1)
          !   beam_out(:,3) = beam_out(:,1)
          !   write(*,*) 'Warning: Setting polarized beams equal to temperature beam'
          !end if
       end if



       allocate(listpix(0:npix_in-1))
       map_buffer = map_in
       do j = 1, nmaps
          do i = 0, npix_in-1
             if (abs(map_in(i,j)) > 1e30) then
                if (r_fill >= 0.d0) then
                   call pix2vec_ring(nside_in, i, vec)
                   do r = 1, r_fill
                      rf = (r+0.2d0) * sqrt(4.d0*pi/npix_in)
                      call query_disc(nside_in, vec, rf, listpix, nlist)
                      tot = 0.d0
                      m   = 0.d0
                      do k = 0, nlist-1
                         if (abs(map_buffer(listpix(k),j)) < 1.d30) then
                            tot = tot + map_buffer(listpix(k),j)
                            m   = m   + 1
                         end if
                      end do
                      if (m > 0) then
                         map_in(i,j) = tot / m
                         exit
                      end if
                   end do

                   if (r > r_fill) then
                      call pix2ang_ring(nside_in, i, theta, phi)
                      write(*,*) 'Error: Filter did not return valid pixel; increase r_fill'
                      write(*,*) '     p         = ', i
                      write(*,*) '     theta     = ', 90.d0-theta*180/pi
                      write(*,*) '     phi       = ', phi*180/pi
                      write(*,*) '     radius    = ', r_fill
                      write(*,*) '     neighbors = ', listpix(0:nlist-1)
                      write(*,*) '     vals      = ', map_buffer(listpix(0:nlist-1),j)
                      stop
                   end if
                else
                   map_in(i,j) = 0.d0
                end if
             end if
          end do
       end do
       deallocate(listpix, map_buffer)

       ! Compute the spherical harmonics transform
       allocate(alms(nmaps, 0:lmax, 0:lmax))

       alms = cmplx(0.,0.)
       if (nmaps == 1) then
          call map2alm(nside_in, lmax, lmax, map_in(:,1), alms, zbounds, weights)
       end if

       ! Deconvolve old beam and pixel window, and convolve new beam and pixel window
       do i = 1, nmaps
          alms(i,0:lmin-1,:) = cmplx(0.,0.)
          do l = lmin, lmax
             if (abs(pixwin_in(l,i)*beam_in(l,i)) < 0.00000001d0) then
                if (l > 1 .or. i == 1) write(*,*) 'Multipole l = ', l, ' set to zero'
                alms(i,l,0:l) = cmplx(0.,0.)
             else
                alms(i,l,0:l) = alms(i,l,0:l) * &
                     & (pixwin_out(l,i)*beam_out(l,i)) / (pixwin_in(l,i)*beam_in(l,i))
             end if
          end do
       end do

       ! Compute the inverse spherical harmonics
       if (nmaps == 1) then
          call alm2map(nside, lmax, lmax, alms, map(:,map_count))
       end if

  
       deallocate(beam_in)
       deallocate(beam_out)
       deallocate(beam_in2)
       deallocate(beam_out2)
       deallocate(map_in)
       deallocate(alms)

    end do
    ! after each map of map_in2 has been smoothed as if nmaps = 1, then we return and deallocate the rest
    if (ordering == 2) then
       do i = 1, nmaps2
          call convert_ring2nest(nside, map(:,i))
       end do
    end if

    nmaps=nmaps2 !for output nmaps

    deallocate(map_in2)
    deallocate(weights)
    deallocate(pixwin_in)
    deallocate(pixwin_out)

  end subroutine smooth_map_zerospin

  subroutine smooth_rms_true(infile, r_fill, lmin, lmax, nside, ordering, & 
       nmaps_in, map_out, header, verbose, beamfile_in, beamfile_out, fwhm_in, fwhm_out)
    implicit none

    character(len=512),                         intent(in)           :: infile
    real(dp),                                   intent(inout)        :: r_fill
    integer(i4b),                               intent(in)           :: nside, lmin, lmax
    integer(i4b),                               intent(in)           :: verbose
    real(dp),                                   intent(in), optional :: fwhm_in, fwhm_out
    character(len=512),                         intent(in), optional :: beamfile_in, beamfile_out
    integer(i4b),                               intent(out)          :: nmaps_in, ordering
    real(dp),          pointer, dimension(:,:)                       :: map_out
    character(len=80),          dimension(180)                       :: header

    integer(i4b)       :: nside_in, npix, npix_in, nmaps, lmax2
    integer(i4b)       :: i, j, k, l, m, p, r, s, nlist, sum_pix
    real(dp)           :: nullval, tot, missval, a_scale
    real(dp)           :: sigma_sq, nullval_dp, vec(3), theta, phi, rf
    logical(lgt)       :: anynull
    character(len=512) :: winfile_in, winfile_out
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     allocatable, dimension(:,:)   :: map_in, map_buffer
    real(dp),     pointer,     dimension(:,:)   :: pixwin_in, pixwin_out
    real(dp),     pointer,     dimension(:,:)   :: weights
    real(dp),     allocatable, dimension(:,:)   :: beam_in, beam_out
    real(dp),                  dimension(2)     :: zbounds = 0.d0
    integer(i4b), allocatable, dimension(:)     :: listpix

    missval = -1.6375d30
    
    ! Read input map
    i = getsize_fits(infile, nside=nside_in, ordering=ordering, nmaps=nmaps_in)
    if (nmaps_in /= 1 .and. nmaps_in /= 3) then
       if (nmaps_in < 3) then
          nmaps_in = 1
       else
          nmaps_in = 3
       end if
    end if

    nmaps   = 1  ! we are smoothong all RMS maps like temperature maps (spin-zero-maps)
    npix    = nside2npix(nside)
    npix_in = nside2npix(nside_in)

    if (present(fwhm_in) .and. present(fwhm_out)) then
       if (fwhm_out < fwhm_in) then
          write(*,*) 'Beam of output is smaller than that of input. Exiting.'
          stop
       else if (fwhm_out == fwhm_in) then
          write(*,*) "Beams of input and output are equal. Use 'ud_grade_rms' instead. Exiting."
          stop
       end if
    end if

    allocate(map_in(0:npix_in-1,nmaps_in))
    allocate(map_buffer(0:npix_in-1,nmaps_in))
    allocate(map_out(0:npix-1,nmaps_in))
    allocate(listpix(0:npix_in-1))
    allocate(alms(nmaps, 0:lmax, 0:lmax))

    call read_bintab(infile, map_in, npix_in, nmaps_in, nullval, anynull, header=header)

    if (ordering == 2) then
       do i = 1, nmaps_in
          call convert_nest2ring(nside_in, map_in(0:npix_in-1,i))
       end do
    end if

    ! Read pixel windows
    call read_pixwin(nside_in, nmaps, pixwin_in)
    call read_pixwin(nside, nmaps, pixwin_out)
    call read_ringweights(nside_in, weights)

    ! Create or read beams
    allocate(beam_in(0:lmax, nmaps))
    allocate(beam_out(0:lmax, nmaps))

    if (present(fwhm_in)) then
       call generate_beam(abs(fwhm_in), lmax, beam_in)
       if (fwhm_in < 0.) beam_in = 1.d0 / beam_in
    else
       call generate_beam(0.d0, lmax, beam_in, beamfile_in)
       !if (nmaps == 3) then
       !   beam_in(:,2) = beam_in(:,1)
       !   beam_in(:,3) = beam_in(:,1)
       !   write(*,*) 'Warning: Setting polarized beams equal to temperature beam'
       !end if
    end if

    if (present(fwhm_out)) then
       call generate_beam(abs(fwhm_out), lmax, beam_out)
    else
       call generate_beam(0.d0, lmax, beam_out, beamfile_out)
       !if (nmaps == 3) then
       !   beam_out(:,2) = beam_out(:,1)
       !   beam_out(:,3) = beam_out(:,1)
       !   write(*,*) 'Warning: Setting polarized beams equal to temperature beam'
       !end if
    end if

    !removing missing pixels and invalid values
    ! fill in missing pixels before smoothing
    map_buffer = map_in
    do j = 1, nmaps
       do i = 0, npix_in-1
          if (abs(map_in(i,j)) > 1e30) then
             if (r_fill >= 0.d0) then
                call pix2vec_ring(nside_in, i, vec)
                do r = 1, r_fill
                   rf = (r+0.2d0) * sqrt(4.d0*pi/npix_in)
                   call query_disc(nside_in, vec, rf, listpix, nlist)
                   tot = 0.d0
                   m   = 0.d0
                   do k = 0, nlist-1
                      if (abs(map_buffer(listpix(k),j)) < 1.d30) then
                         tot = tot + map_buffer(listpix(k),j)
                         m   = m   + 1
                      end if
                   end do
                   if (m > 0) then
                      map_in(i,j) = tot / m
                      exit
                   end if
                end do

                if (r > r_fill) then
                   call pix2ang_ring(nside_in, i, theta, phi)
                   write(*,*) 'Error: Filter did not return valid pixel; increase r_fill'
                   write(*,*) '     p         = ', i
                   write(*,*) '     theta     = ', 90.d0-theta*180/pi
                   write(*,*) '     phi       = ', phi*180/pi
                   write(*,*) '     radius    = ', r_fill
                   write(*,*) '     neighbors = ', listpix(0:nlist-1)
                   write(*,*) '     vals      = ', map_buffer(listpix(0:nlist-1),j)
                   stop
                end if
             else
                map_in(i,j) = 0.d0
             end if
          end if
       end do
    end do


    if (verbose > 0) then
       write(*,*) ""
       write(*,*) "Smoothing input (in Variance domain)"
       write(*,*) ""
    end if

    !smooth input RMS map to map_out
    ! Do one map at a time; has to be smoothed as spin-zero-maps (i.e. Temp maps)
    ! The smoothing has to be done on the variance map with the square beam,
    ! i.e., a beam that is FWHM/sqrt(2),
    ! where FWHM is the effective beam given by
    ! FWHM = sqrt(1 - (FWHM_in/FWHM_out)**2) * FWHM_out
    ! The effective beam has the same beam weights (per l) as beam(FWHM_out)/beam(FWHM_in)
    ! so we do not need to calculate FWHM_effective, and it is not possible if we have beam files
    ! Further we use the fact that (in Stokes I, i.e. temperature)
    ! beam(fwhm') = beam(fwhm)**[(fwhm' /fwhm)**2] for all l > 0
    ! so that beam(fwhm/sqrt(2)) = beam(fwhm)**1/2 = sqrt(beam(fwhm))
    
    map_buffer = map_in*map_in !getting the variance of the input

    do p = 1, nmaps_in
       if (verbose > 1) write(*,*) "      Polarization: ",p
          
       alms = cmplx(0.,0.)
       call map2alm(nside_in, lmax, lmax, map_buffer(:,p), alms, zbounds, weights)
       
       ! Deconvolve old beam and pixel window, and convolve new beam and pixel window
       do i = 1, 1
          alms(i,0:lmin-1,:) = cmplx(0.,0.)
          do l = lmin, lmax
             if (abs(pixwin_in(l,i)*beam_in(l,i)) < 0.00000001d0) then
                if (l > 1 .or. i == 1) write(*,*) 'Multipole l = ', l, ' set to zero'
                alms(i,l,0:l) = cmplx(0.,0.)
             else
                !smooth with 1/sqrt(2) of effective beam
                alms(i,l,0:l) = alms(i,l,0:l) * &
                     & (pixwin_out(l,i)/pixwin_in(l,i)) * dsqrt(beam_out(l,i)/beam_in(l,i)) 
             end if
          end do
       end do

       ! Compute the inverse spherical harmonics to get the smoother variance map
       call alm2map(nside, lmax, lmax, alms, map_out(:,p))
    end do


    !take square-root of variance map to get RMS
    map_out=dsqrt(map_out)

    ! if nside in /= nside out, rescale the smoothed RMS map to ensure equal
    ! chi-squared values 
    if ( nside /= nside_in ) then
       map_out=map_out*(nside*1.d0/nside_in)

    end if


    ! go back to nested ordering if input was nested
    if (ordering == 2) then
       do i = 1, nmaps_in
          call convert_ring2nest(nside, map_out(:,i))
       end do
    end if

    deallocate(listpix)
    deallocate(weights)
    deallocate(pixwin_in)
    deallocate(pixwin_out)
    deallocate(beam_in)
    deallocate(beam_out)
    deallocate(map_in)
    deallocate(alms)
    deallocate(map_buffer)

  end subroutine smooth_rms_true

  function smooth_rms_map(rms_in, r_fill, lmin, lmax, nside, ordering, & 
       nmaps_in, map_out, header, verbose, beamfile_in, beamfile_out, fwhm_in, fwhm_out)
    implicit none

    real(dp),                   dimension(:,:), intent(in)           :: rms_in
    real(dp),                                   intent(inout)        :: r_fill
    integer(i4b),                               intent(in)           :: nside, lmin, lmax
    ! integer(i4b),                               intent(in)           :: verbose
    real(dp),                                   intent(in), optional :: fwhm_in, fwhm_out
    character(len=512),                         intent(in), optional :: beamfile_in, beamfile_out
    integer(i4b),                               intent(out)          :: nmaps_in, ordering
    real(dp),          pointer, dimension(:,:)                       :: map_out
    character(len=80),          dimension(180)                       :: header

    integer(i4b)       :: nside_in, npix, npix_in, nmaps, lmax2
    integer(i4b)       :: i, j, k, l, m, p, r, s, nlist, sum_pix
    real(dp)           :: nullval, tot, missval, a_scale
    real(dp)           :: sigma_sq, nullval_dp, vec(3), theta, phi, rf
    logical(lgt)       :: anynull
    character(len=512) :: winfile_in, winfile_out
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     allocatable, dimension(:,:)   :: map_in, map_buffer
    real(dp),     pointer,     dimension(:,:)   :: pixwin_in, pixwin_out
    real(dp),     pointer,     dimension(:,:)   :: weights
    real(dp),     allocatable, dimension(:,:)   :: beam_in, beam_out
    real(dp),                  dimension(2)     :: zbounds = 0.d0
    integer(i4b), allocatable, dimension(:)     :: listpix

    missval = -1.6375d30
    
    ! Read input map
    i = getsize_fits(infile, nside=nside_in, ordering=ordering, nmaps=nmaps_in)
    if (nmaps_in /= 1 .and. nmaps_in /= 3) then
       if (nmaps_in < 3) then
          nmaps_in = 1
       else
          nmaps_in = 3
       end if
    end if

    nmaps   = 1  ! we are smoothong all RMS maps like temperature maps (spin-zero-maps)
    npix    = nside2npix(nside)
    npix_in = nside2npix(nside_in)

    if (present(fwhm_in) .and. present(fwhm_out)) then
       if (fwhm_out < fwhm_in) then
          write(*,*) 'Beam of output is smaller than that of input. Exiting.'
          stop
       else if (fwhm_out == fwhm_in) then
          write(*,*) "Beams of input and output are equal. Use 'ud_grade_rms' instead. Exiting."
          stop
       end if
    end if

    allocate(map_in(0:npix_in-1,nmaps_in))
    allocate(map_buffer(0:npix_in-1,nmaps_in))
    allocate(map_out(0:npix-1,nmaps_in))
    allocate(listpix(0:npix_in-1))
    allocate(alms(nmaps, 0:lmax, 0:lmax))

    call read_bintab(infile, map_in, npix_in, nmaps_in, nullval, anynull, header=header)

    if (ordering == 2) then
       do i = 1, nmaps_in
          call convert_nest2ring(nside_in, map_in(0:npix_in-1,i))
       end do
    end if

    ! Read pixel windows
    call read_pixwin(nside_in, nmaps, pixwin_in)
    call read_pixwin(nside, nmaps, pixwin_out)
    call read_ringweights(nside_in, weights)

    ! Create or read beams
    allocate(beam_in(0:lmax, nmaps))
    allocate(beam_out(0:lmax, nmaps))

    if (present(fwhm_in)) then
       call generate_beam(abs(fwhm_in), lmax, beam_in)
       if (fwhm_in < 0.) beam_in = 1.d0 / beam_in
    else
       call generate_beam(0.d0, lmax, beam_in, beamfile_in)
       !if (nmaps == 3) then
       !   beam_in(:,2) = beam_in(:,1)
       !   beam_in(:,3) = beam_in(:,1)
       !   write(*,*) 'Warning: Setting polarized beams equal to temperature beam'
       !end if
    end if

    if (present(fwhm_out)) then
       call generate_beam(abs(fwhm_out), lmax, beam_out)
    else
       call generate_beam(0.d0, lmax, beam_out, beamfile_out)
       !if (nmaps == 3) then
       !   beam_out(:,2) = beam_out(:,1)
       !   beam_out(:,3) = beam_out(:,1)
       !   write(*,*) 'Warning: Setting polarized beams equal to temperature beam'
       !end if
    end if

    !removing missing pixels and invalid values
    ! fill in missing pixels before smoothing
    map_buffer = map_in
    do j = 1, nmaps
       do i = 0, npix_in-1
          if (abs(map_in(i,j)) > 1e30) then
             if (r_fill >= 0.d0) then
                call pix2vec_ring(nside_in, i, vec)
                do r = 1, r_fill
                   rf = (r+0.2d0) * sqrt(4.d0*pi/npix_in)
                   call query_disc(nside_in, vec, rf, listpix, nlist)
                   tot = 0.d0
                   m   = 0.d0
                   do k = 0, nlist-1
                      if (abs(map_buffer(listpix(k),j)) < 1.d30) then
                         tot = tot + map_buffer(listpix(k),j)
                         m   = m   + 1
                      end if
                   end do
                   if (m > 0) then
                      map_in(i,j) = tot / m
                      exit
                   end if
                end do

                if (r > r_fill) then
                   call pix2ang_ring(nside_in, i, theta, phi)
                   write(*,*) 'Error: Filter did not return valid pixel; increase r_fill'
                   write(*,*) '     p         = ', i
                   write(*,*) '     theta     = ', 90.d0-theta*180/pi
                   write(*,*) '     phi       = ', phi*180/pi
                   write(*,*) '     radius    = ', r_fill
                   write(*,*) '     neighbors = ', listpix(0:nlist-1)
                   write(*,*) '     vals      = ', map_buffer(listpix(0:nlist-1),j)
                   stop
                end if
             else
                map_in(i,j) = 0.d0
             end if
          end if
       end do
    end do


    ! if (verbose > 0) then
    !    write(*,*) ""
    !    write(*,*) "Smoothing input (in Variance domain)"
    !    write(*,*) ""
    ! end if

    !smooth input RMS map to map_out
    ! Do one map at a time; has to be smoothed as spin-zero-maps (i.e. Temp maps)
    ! The smoothing has to be done on the variance map with the square beam,
    ! i.e., a beam that is FWHM/sqrt(2),
    ! where FWHM is the effective beam given by
    ! FWHM = sqrt(1 - (FWHM_in/FWHM_out)**2) * FWHM_out
    ! The effective beam has the same beam weights (per l) as beam(FWHM_out)/beam(FWHM_in)
    ! so we do not need to calculate FWHM_effective, and it is not possible if we have beam files
    ! Further we use the fact that (in Stokes I, i.e. temperature)
    ! beam(fwhm') = beam(fwhm)**[(fwhm' /fwhm)**2] for all l > 0
    ! so that beam(fwhm/sqrt(2)) = beam(fwhm)**1/2 = sqrt(beam(fwhm))
    
    map_buffer = map_in*map_in !getting the variance of the input

    do p = 1, nmaps_in
       ! if (verbose > 1) write(*,*) "      Polarization: ",p
          
       alms = cmplx(0.,0.)
       call map2alm(nside_in, lmax, lmax, map_buffer(:,p), alms, zbounds, weights)
       
       ! Deconvolve old beam and pixel window, and convolve new beam and pixel window
       do i = 1, 1
          alms(i,0:lmin-1,:) = cmplx(0.,0.)
          do l = lmin, lmax
             if (abs(pixwin_in(l,i)*beam_in(l,i)) < 0.00000001d0) then
                if (l > 1 .or. i == 1) write(*,*) 'Multipole l = ', l, ' set to zero'
                alms(i,l,0:l) = cmplx(0.,0.)
             else
                !smooth with 1/sqrt(2) of effective beam
                alms(i,l,0:l) = alms(i,l,0:l) * &
                     & (pixwin_out(l,i)/pixwin_in(l,i)) * dsqrt(beam_out(l,i)/beam_in(l,i)) 
             end if
          end do
       end do

       ! Compute the inverse spherical harmonics to get the smoother variance map
       call alm2map(nside, lmax, lmax, alms, map_out(:,p))
    end do


    !take square-root of variance map to get RMS
    map_out=dsqrt(map_out)

    ! if nside in /= nside out, rescale the smoothed RMS map to ensure equal
    ! chi-squared values 
    if ( nside /= nside_in ) then
       map_out=map_out*(nside*1.d0/nside_in)

    end if


    ! go back to nested ordering if input was nested
    if (ordering == 2) then
       do i = 1, nmaps_in
          call convert_ring2nest(nside, map_out(:,i))
       end do
    end if

    deallocate(listpix)
    deallocate(weights)
    deallocate(pixwin_in)
    deallocate(pixwin_out)
    deallocate(beam_in)
    deallocate(beam_out)
    deallocate(map_in)
    deallocate(alms)
    deallocate(map_buffer)

  end function smooth_rms_map


  subroutine print_maximum(mapfile, component, fwhm)
    implicit none

    character(len=*),                          intent(in)           :: mapfile
    integer(i4b),                              intent(in)           :: component
    real(dp),                                  intent(in)           :: fwhm

    integer(i4b) :: npix, pix, nside, nmaps, ordering, i, nlist, nest
    real(dp)     :: nullval, max_val
    real(dp)     :: theta, phi, t1, t2, integral, radius
    logical(lgt) :: anynull
    type(planck_rng) :: rng_handle
    real(dp),    allocatable, dimension(:,:)   :: map
    real(dp),                 dimension(3)     :: vector0
    integer(i4b), allocatable, dimension(:) :: listpix

    call cpu_time(t1)

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull)

    if (nmaps < component) then
       write(*,*) 'Requested component is not present in map file'
       stop
    end if

    pix    = 0
    max_val = -1.d30
    do i = 0, npix-1
       if (map(i,component) /= -1.6375e30) then
          if (map(i,component) > max_val) then
             pix    = i
             max_val = map(i,component)
          end if
       end if
    end do

    goto 1010

    ! Search for best-fit pixel
    radius = fwhm * pi/180.d0/60 * 4.d0
    allocate(listpix(0:npix-1))
    if (ordering == 1) then
       call pix2vec_ring(nside, pix, vector0)
       call query_disc(nside, vector0, radius, listpix, nlist)
    else
       call pix2vec_nest(nside, pix, vector0)
       call query_disc(nside, vector0, radius, listpix, nlist, nest=1)
    end if

    integral = convolve_with_gaussian(nside, ordering, map(:,component), pix, fwhm)

    pix     = 0
    max_val = -1.d30
    do i = 0, nlist-1
       integral = convolve_with_gaussian(nside, ordering, map(:,component), listpix(i), fwhm)
       if (integral > max_val) then
          pix     = listpix(i)
          max_val = integral
!          write(*,*) pix, real(map(pix,component),sp), max_val
       end if
    end do

    call cpu_time(t2)
!    write(*,*) t2-t1

1010 if (ordering == 1) then
       call pix2ang_ring(nside, pix, theta, phi)
    else
       call pix2ang_nest(nside, pix, theta, phi)
    end if

    write(*,*) pix, phi*180./pi, 90.-180./pi*theta

    deallocate(map)

  end subroutine print_maximum

  subroutine print_stats(mapfile, maskfile)
    implicit none

    character(len=*),                          intent(in)           :: mapfile, maskfile

    integer(i4b) :: npix, pix, nside, nmaps, ordering, i, j, k 
    integer(i4b) :: nlist, nest, order_mask, misspix(3)
    real(dp)     :: nullval, max_val, med
    real(dp)     :: theta, phi, t1, t2, integral, radius, mu, sigma, my_min, my_max, scale
    real(dp), parameter :: missval=-1.6375d30
    logical(lgt) :: anynull
    type(planck_rng) :: rng_handle
    character(len=512) :: par
    real(dp),    allocatable, dimension(:,:)   :: map, mask
    real(dp),    allocatable, dimension(:)     :: vals

    call cpu_time(t1)

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull)

    ! Read input mask
    if (maskfile == 'dummy') then
       allocate(mask(0:npix-1,nmaps))
       mask = 1.0
    else
       i = getsize_fits(maskfile, nside=nside, ordering=order_mask)
       allocate(mask(0:npix-1,nmaps))
       call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
       if (order_mask /= ordering) then
          if (order_mask == 1) then
             do i = 1, nmaps
                call convert_ring2nest(nside, mask(:,i))
             end do
          else
             do i = 1, nmaps
                call convert_nest2ring(nside, mask(:,i))
             end do
          end if
       end if
    endif
    misspix=0
    do i = 1,nmaps
       do j = 0, npix-1
          if (mask(j,i) > 0.5) then
             if (abs((map(j,i)-missval)/missval) < 1d-5) misspix(i)=misspix(i)+1
          end if
       end do
    end do
    where (abs((map-missval)/missval) < 1d-5)
       mask = 0.d0
    end where
    
    allocate(vals(npix))

    write(*,*) '   -------------------------'
    write(*,*) '   nside    = ', nside
    write(*,*) '   nmaps    = ', nmaps
    write(*,*) '   ordering = ', ordering
    write(*,*) '   -------------------------'

    do i = 1, nmaps
    !do i = 2, 2
       
       k = 0
       do j = 0, npix-1
          if (mask(j,i) > 0.5d0) then
             k = k+1
             vals(k) = map(j,i)
          end if
       end do
       mu    = sum(vals(1:k))/k
       sigma = sqrt(sum((vals(1:k)-mu)**2)/(k-1))
       if (iargc()==4) then
          med   = median(abs(vals(1:k)))
          call getarg(4, par)
          read(par,*) scale
          if (scale > 0.d0) then
             write(*,fmt='(f6.2)') sigma*scale
          else if (scale < 0) then
             write(*,fmt='(f6.3)') med*abs(scale)
          else if (scale == 0.d0) then
             write(*,fmt='(a,f8.3,a,f8.3,a)') '$', mu, "\pm", sigma, '$'
          end if
       else
          write(*,*) ''
          write(*,*) '   Pol         = ', i
          write(*,*) '   Mean        = ', mu
          write(*,*) '   RMS         = ', sigma
          write(*,*) '   Min         = ', minval(map(:,i), mask(:,i) > 0.5d0)
          write(*,*) '   Max         = ', maxval(map(:,i), mask(:,i) > 0.5d0)
          write(*,*) '   Missing pix = ', misspix(i)
          write(*,*) ''
       end if
    end do
    
    

    deallocate(map, mask)

  end subroutine print_stats

  subroutine print_stats_col(mapfile, maskfile)
    implicit none

    character(len=*),                          intent(in)           :: mapfile, maskfile

    integer(i4b) :: npix, pix, nside, nmaps, ordering, i, j, k, nlist, nest, order_mask, col, extno
    real(dp)     :: nullval, max_val, med
    real(dp)     :: theta, phi, t1, t2, integral, radius, mu, sigma, my_min, my_max, scale
    logical(lgt) :: anynull
    type(planck_rng) :: rng_handle
    character(len=512) :: par
    real(dp),    allocatable, dimension(:,:)   :: map, mask
    real(dp),    allocatable, dimension(:)     :: vals

    call cpu_time(t1)

    call getarg(5,par)
    read(par,*) extno

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps, extno=extno)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull, extno=extno)

    ! Read input mask
    i = getsize_fits(maskfile, nside=nside, ordering=order_mask)
    allocate(mask(0:npix-1,1))
    call read_bintab(maskfile, mask, npix, 1, nullval, anynull)
    if (order_mask /= ordering) then
       if (order_mask == 1) then
          call convert_ring2nest(nside, mask(:,1))
       else
          call convert_nest2ring(nside, mask(:,1))
       end if
    end if

    where (map < -1.637d30)
       mask = 0.d0
    end where

    call getarg(4,par)
    read(par,*) col
    call getarg(6,par)
    read(par,*) scale

    allocate(vals(npix))
    k = 0
    do j = 0, npix-1
       if (mask(j,1) > 0.5d0) then
          k = k+1
          if (trim(mapfile) == 'data/products/COM_CompMap_DustPol-commander_1024_R2.00.fits' .or. &
               & trim(mapfile) == 'data/products/COM_CompMap_SynchrotronPol-commander_0256_R2.00.fits') then
             vals(k) = sqrt(map(j,col)**2 + map(j,col+1)**2)
          else
             vals(k) = map(j,col)
          end if
       end if
    end do
    
    mu    = sum(vals(1:k))/k
    sigma = sqrt(sum((vals(1:k)-mu)**2)/(k-1))
    write(*,fmt='(a,f8.2,a,f8.2,a)') '$', mu*scale, "\pm", sigma*scale, '$'

    deallocate(map, mask)

  end subroutine print_stats_col

  function convolve_with_gaussian(nside, ordering, map, pix, fwhm)
    implicit none

    real(dp),     dimension(0:), intent(in) :: map
    integer(i4b),                intent(in) :: pix, nside, ordering
    real(dp),                    intent(in) :: fwhm
    real(dp)                                :: convolve_with_gaussian
    
    real(dp) :: s, radius, sigma
    integer(i4b) :: i, nlist, npix, nest, n
    integer(i4b), allocatable, dimension(:) :: listpix
    real(dp), dimension(3) :: vector0, vector

    npix   = size(map)
    sigma  = fwhm * pi/180.d0/60.d0 / sqrt(8.d0*log(2.d0))
    radius = 2.d0 * fwhm * pi/180.d0/60.d0
    nest   = ordering-1

    allocate(listpix(0:npix-1))
    
    if (ordering == 1) then
       call pix2vec_ring(nside, pix, vector0)
       call query_disc(nside, vector0, radius, listpix, nlist)
    else
       call pix2vec_nest(nside, pix, vector0)
       call query_disc(nside, vector0, radius, listpix, nlist,nest=1)
    end if
    

    s = 0.d0
    n = 0
    do i = 0, nlist-1
       if (map(listpix(i)) /= -1.6375e30) then
          if (ordering == 1) then
             call pix2vec_ring(nside, listpix(i), vector)
          else
             call pix2vec_nest(nside, listpix(i), vector)
          end if
          if (listpix(i) /= pix) then
             radius = acos(sum(vector * vector0))
          else
             radius = 0.d0
          end if
          s = s + map(listpix(i)) * exp(-0.5d0*(radius/sigma)**2)
          n = n + 1
       end if
    end do

    if (n == nlist) then
       convolve_with_gaussian = s / real(n,dp)
    else
       convolve_with_gaussian = -1.d30
    end if

    deallocate(listpix)

  end function convolve_with_gaussian

  subroutine add_gaussian_noise(mapfile, rmsfile, seed, nside, ordering, nmaps, map, header, sigma_0)
    implicit none

    character(len=*),                          intent(in)           :: mapfile, rmsfile
    integer(i4b),                              intent(inout)        :: seed
    integer(i4b),                              intent(out)          :: nside, ordering, nmaps
    real(dp),         pointer, dimension(:,:)                       :: map   
    character(len=80),         dimension(180)                       :: header
    real(dp),                                  intent(in), optional :: sigma_0

    integer(i4b) :: npix, nside_rms, ordering_rms, nmaps_rms, i, j
    real(dp)     :: nullval
    logical(lgt) :: anynull
    type(planck_rng) :: rng_handle
    real(dp),    allocatable, dimension(:,:)   :: rmsmap

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull, header=header)

    ! Read rms file and check consistency
    i = getsize_fits(rmsfile, nside=nside_rms, ordering=ordering_rms, nmaps=nmaps_rms)

    if (nside_rms /= nside) then
       write(*,*) 'Error: Different Nsides of the two files.'
       stop
    end if

    if (nmaps_rms /= nmaps) then
       write(*,*) 'Error: Different Nmaps of the two files.'
       stop
    end if

    allocate(rmsmap(0:npix-1,nmaps))
    call read_bintab(rmsfile, rmsmap, npix, nmaps, nullval, anynull)

    if (ordering_rms /= ordering) then
       if (ordering_rms == 1) then
          do i = 1, nmaps
             call convert_ring2nest(nside, rmsmap(:,i))
          end do
       else
          do i = 1, nmaps
             call convert_nest2ring(nside, rmsmap(:,i))
          end do
       end if
    end if

    call rand_init(rng_handle, seed)

    ! Add Gaussian noise to the input map
    if (present(sigma_0)) then
       ! Interpret RMS map as Nobs, not RMS
       do j = 1, nmaps
          do i = 0, npix-1
             if (abs((map(i,j) -( -1.6375e30))/-1.6375e30) > 1.d-5 .and. rmsmap(i,j) >= 0.) then
                map(i,j) = map(i,j) + sigma_0 / sqrt(rmsmap(i,j)) * rand_gauss(rng_handle)
             else 
                map(i,j) = -1.6375e30
             end if
          end do
       end do
    else
       do j = 1, nmaps
          do i = 0, npix-1
             if (map(i,j) /= -1.6375e30 .and. rmsmap(i,j) >= 0.) then
                map(i,j) = map(i,j) + rmsmap(i,j) * rand_gauss(rng_handle)
             else 
                map(i,j) = -1.6375e30
             end if
          end do
       end do
    end if

    deallocate(rmsmap)

  end subroutine add_gaussian_noise


  subroutine add_gaussian_noise_sqrt_N(mapfile, rmsfile, map2mask_file, seed, nside, ordering, nmaps, map, header, sigma_0)
    implicit none

    character(len=*),                          intent(in)           :: mapfile, rmsfile, map2mask_file
    integer(i4b),                              intent(inout)        :: seed
    integer(i4b),                              intent(out)          :: nside, ordering, nmaps
    real(dp),         pointer, dimension(:,:)                       :: map   
    real(dp),                                  intent(in), optional :: sigma_0
    character(len=80),         dimension(180)                       :: header

    integer(i4b) :: npix, nside_rms, ordering_rms, nmaps_rms, i, j, k, polarization
    real(dp)     :: nullval
    integer(i8b) :: n
    logical(lgt) :: anynull, inv
    type(planck_rng) :: rng_handle
    real(dp),    allocatable, dimension(:,:)   :: rmsmap, map2mask

    real(dp), allocatable, dimension(:,:) :: sqrt_N
    real(dp), allocatable, dimension(:)   :: eta

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull, header=header)

    ! Read map2mask file
    i = getsize_fits(map2mask_file, nside=nside_rms, ordering=ordering_rms, nmaps=nmaps_rms)

    if (nside_rms /= nside) then
       write(*,*) 'Error: Different Nsides of the two files.'
       stop
    end if

    if (nmaps_rms /= nmaps) then
       write(*,*) 'Error: Different Nmaps of the two files.'
       stop
    end if

    allocate(map2mask(0:npix-1,nmaps))
    call read_bintab(map2mask_file, map2mask, npix, nmaps, nullval, anynull)

    ! Read covariance file
    call read_covmatrix(58, rmsfile, ordering_rms, polarization, sqrt_N, inv, n)

    if (ordering_rms /= ordering) then
       if (ordering_rms == 1) then
          do i = 1, nmaps
             call convert_ring2nest(nside, rmsmap(:,i))
          end do
       else
          do i = 1, nmaps
             call convert_nest2ring(nside, rmsmap(:,i))
          end do
       end if
    end if

    call rand_init(rng_handle, seed)

    ! Draw a random gaussian vector
    allocate(eta(n))
    do i = 1, n
       eta(i) = rand_gauss(rng_handle)
    end do
    eta = matmul(sqrt_N, eta)

    ! Add noise to the input map
    k = 1
    do j = 1, nmaps
       do i = 0, npix-1
          if (nint(map2mask(i,j)) > -1) then
             map(i,j) = map(i,j) + eta(k)
             k        = k+1
          else
             map(i,j) = -1.6375d30
          end if
       end do
    end do

    deallocate(map2mask)
    deallocate(eta)
    deallocate(sqrt_N)

  end subroutine add_gaussian_noise_sqrt_N


  subroutine subtract_mono_dipole(mapfile, maskfile, nside, ordering, map, header, md)
    implicit none

    character(len=*),                          intent(in)           :: mapfile, maskfile
    integer(i4b),                              intent(out)          :: nside, ordering
    real(dp),         pointer, dimension(:,:)                       :: map   
    character(len=80),         dimension(180)                       :: header
    real(dp),                  dimension(4),   intent(in), optional :: md

    integer(i4b) :: npix, nside_mask, ordering_mask, nmaps, i, j, k, degree
    real(dp)     :: nullval, tot
    real(dp)     :: fmissval 
    logical(lgt) :: anynull
    type(planck_rng) :: rng_handle
    real(dp),    allocatable, dimension(:,:)       :: mask, harmonics, harmonics2
    real(dp),                 dimension(0:3)       :: multipoles, b
    real(dp),                 dimension(3)         :: alb         ! amplitude, longitude and lattitude
    real(dp),                 dimension(3)         :: vector
    real(dp),                 dimension(0:3,0:3)   :: A
    real(dp),                 dimension(2)         :: zbounds = 0.d0

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(:,i))
       end do
    end if
    ordering = 1

    if (present(md)) then
       write(*,*) 'md input = ', real(md,sp)
       do i = 0, npix-1
          if (ordering == 1) then
             call pix2vec_ring(nside, i, vector)
          else
             call pix2vec_nest(nside, i, vector)          
          end if
          map(i,1) = map(i,1) - md(1)
          do j = 1, 3
             map(i,1) = map(i,1) - md(j+1) * vector(j)
          end do

!          map(i,1) = md(1)
!          do j = 1, 3
!             map(i,1) = map(i,1) + md(j+1) * vector(j)
!          end do
       end do
       return
    end if

    ! Read rms file and check consistency
    i = getsize_fits(maskfile, nside=nside_mask, ordering=ordering_mask)

    if (nside_mask /= nside) then
       write(*,*) 'Error: Different Nsides of the two files.'
       stop
    end if

    allocate(mask(0:npix-1,1))
    call read_bintab(maskfile, mask, npix, 1, nullval, anynull)

    if (ordering_mask == 2) then
       call convert_nest2ring(nside, mask(:,1))
    end if

    allocate(harmonics(0:npix-1,0:3))
    allocate(harmonics2(0:npix-1,0:3))
    do i = 0, npix-1
       if (ordering == 1) then
          call pix2vec_ring(nside, i, vector)
       else
          call pix2vec_nest(nside, i, vector)          
       end if
          
       if (mask(i,1) == 1.) then
          harmonics(i,0) = 1.d0
          harmonics(i,1) = vector(1)
          harmonics(i,2) = vector(2)
          harmonics(i,3) = vector(3)
       else
          harmonics(i,:) = 0.d0
       end if

       harmonics2(i,0) = 1.d0
       harmonics2(i,1) = vector(1)
       harmonics2(i,2) = vector(2)
       harmonics2(i,3) = vector(3)

    end do

    A = 0.d0
    b = 0.d0
    do j = 0, 3
       do k = 0, 3
          A(j,k) = sum(harmonics(:,j) * harmonics(:,k))
       end do
       b(j) = sum(map(:,1) * harmonics(:,j))
    end do


    call solve_system_real(A, multipoles, b)

    do i = 0, npix-1
!       if (mask(i,1) == 1) then
          do j = 0, 3
             map(i,1) = map(i,1) - multipoles(j) * harmonics2(i,j)
          end do
!       else
!          map(i,:) = -1.6375e30
!       end if
    end do

    write(*,*) 'Coefficients = ', real(multipoles,sp)

    alb(1) = dsqrt(sum(multipoles(1:3)*multipoles(1:3)))
    do i = 1,3 
       if (abs(multipoles(i)/alb(1))<1.d-10) multipoles(i)=0.d0
    end do
    alb(3) = asind(multipoles(3)/alb(1))
    if (multipoles(1)==0.d0) then
       if (multipoles(2) < 0.d0) then
          alb(2) = 270.d0
       elseif(multipoles(2) > 0.d0) then
          alb(2) = 90.d0
       else
          alb(2) = 0.d0
       end if
    else
       alb(2) = atand(multipoles(2)/multipoles(1))
       if ( multipoles(1) < 0.d0) then
          alb(2) = alb(2) + 180.d0
       else if ( multipoles(2) < 0.d0) then
          alb(2) = alb(2) + 360.d0
       end if
    end if
    write(*,*) 'Dipole parameters: Amplitude    Longitude    Latitude '
    write(*,fmt='(a,f13.6,f12.6,f12.6)') '                  ', alb

!    degree = 2
!    call remove_dipole(nside, map(:,1), ordering, degree, multipoles, zbounds, fmissval, mask(:,1))

!    write(*,*) 'multipoles = ', multipoles

    deallocate(mask)

  end subroutine subtract_mono_dipole


  subroutine subtract_mono_dipole_highl(mapfile, maskfile, lmax, lcut, nside, ordering, map, header, md)
    implicit none

    character(len=*),                          intent(in)           :: mapfile, maskfile
    integer(i4b),                              intent(in)           :: lmax, lcut
    integer(i4b),                              intent(out)          :: nside, ordering
    real(dp),         pointer, dimension(:,:)                       :: map   
    real(dp),                  dimension(4),   intent(in), optional :: md
    character(len=80),         dimension(180)                       :: header

    integer(i4b) :: npix, nside_mask, ordering_mask, nmaps, i, j, k, l, m, degree, numcomp
    real(dp)     :: nullval, tot
    real(dp)     :: nullval_sp
    real(dp)     :: fmissval 
    logical(lgt) :: anynull
    type(planck_rng) :: rng_handle
    real(dp),    allocatable, dimension(:,:)       :: mask, harmonics
    complex(dpc), allocatable, dimension(:,:,:)    :: alms
    real(dp),                 dimension(3)         :: vector
    real(dp), allocatable,    dimension(:)         :: b, multipoles
    real(dp), allocatable,    dimension(:,:)       :: A
    real(dp),                 dimension(2)         :: zbounds = 0.d0

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(:,i))
       end do
    end if
    ordering = 1

    ! Read rms file and check consistency
    i = getsize_fits(maskfile, nside=nside_mask, ordering=ordering_mask)

    if (nside_mask /= nside) then
       write(*,*) 'Error: Different Nsides of the two files.'
       stop
    end if

    allocate(mask(0:npix-1,1))
    call read_bintab(maskfile, mask, npix, 1, nullval, anynull)

    if (ordering_mask == 2) then
       call convert_nest2ring(nside, mask(:,1))
    end if

    numcomp = (lmax+1)**2
    allocate(harmonics(0:npix-1,numcomp), A(numcomp,numcomp), b(numcomp), multipoles(numcomp), alms(1,0:lmax,0:lmax))
    do i = 0, npix-1
       if (ordering == 1) then
          call pix2vec_ring(nside, i, vector)
       else
          call pix2vec_nest(nside, i, vector)          
       end if
       harmonics(i,1) = 1.d0
       harmonics(i,2) = vector(1)
       harmonics(i,3) = vector(2)
       harmonics(i,4) = vector(3)
    end do

    i = 5
    do l = 2, lmax
       do m = 0, l
          write(*,*) l, m
          alms = cmplx(0.d0, 0.d0)
          alms(1,l,m) = 1.d0
          call alm2map(nside, lmax, lmax, alms, harmonics(:,i))
          i = i+1
          if (m > 0) then
             alms(1,l,m) = cmplx(0.d0, 1.d0)
             call alm2map(nside, lmax, lmax, alms, harmonics(:,i))
             i = i+1
          end if
       end do
    end do

    A = 0.d0
    b = 0.d0
    do j = 1, numcomp
       write(*,*) numcomp, j
       do k = j, numcomp
          A(j,k) = sum(harmonics(:,j) * mask(:,1) * harmonics(:,k))
          A(k,j) = A(j,k)
       end do
       b(j) = sum(map(:,1) * mask(:,1) * harmonics(:,j))
    end do


    call solve_linear_system(A, multipoles, b)

    do i = 0, npix-1
       if (map(i,1) > -1.637d30) then
          do j = 1, (lcut+1)**2
             map(i,1) = map(i,1) - multipoles(j) * harmonics(i,j)
          end do
       end if
    end do

    write(*,*) 'Coefficients = ', real(multipoles(1:4),sp)

!    degree = 2
!    call remove_dipole(nside, map(:,1), ordering, degree, multipoles, zbounds, fmissval, mask(:,1))

!    write(*,*) 'multipoles = ', multipoles

    deallocate(mask)

  end subroutine subtract_mono_dipole_highl

  subroutine extract_multipole_range(mapfile, lmin, lmax, nside, ordering, map, header)
    implicit none

    character(len=*),                          intent(in)           :: mapfile
    integer(i4b),                              intent(in)           :: lmin, lmax
    integer(i4b),                              intent(out)          :: nside, ordering
    real(dp),         pointer, dimension(:,:)                       :: map   
    character(len=80),         dimension(180)                       :: header

    integer(i4b) :: npix, nmaps, i, j, k
    real(dp)     :: nullval, tot
    real(dp)     :: fmissval 
    logical(lgt) :: anynull
    logical(lgt), allocatable, dimension(:,:)       :: mask
    real(dp),     allocatable, dimension(:,:)       :: map_lowl
    real(dp),                  dimension(2)         :: zbounds = 0.d0
    complex(dpc), allocatable, dimension(:,:,:)     :: alms
    real(dp),     pointer,     dimension(:,:)      :: weights

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(:,i))
       end do
    end if
    ordering = 1

    ! Set undefined values to zero
    allocate(mask(0:npix,nmaps))
    mask = .true.
    do j = 1, nmaps
       do i = 0, npix-1
          if (map(i,j) == -1.6375e30) then
             mask(i,j) = .false.
             map(i,j)  = 0.d0
          end if
       end do
    end do

    ! Read Healpix ring weights
    call read_ringweights(nside, weights)

    ! Compute spherical harmonic expansion
    allocate(alms(nmaps,0:lmax,0:lmax))
    if (nmaps == 1) then
       call map2alm(nside, lmax, lmax, map(:,1), alms, zbounds, weights)
    else
       call map2alm(nside, lmax, lmax, map, alms, zbounds, weights)
    end if

    alms(:,0:lmin-1,0:lmin-1) = cmplx(0.d0, 0.d0)

    ! Compute low-l map
    allocate(map_lowl(0:npix-1,nmaps))
    if (nmaps == 1) then
       call alm2map(nside, lmax, lmax, alms, map(:,1))
    else
       call alm2map(nside, lmax, lmax, alms, map)       
    end if

    ! Apply mask
    do j = 1, nmaps
       do i = 0, npix-1
          if (.not. mask(i,j)) then
             map(i,j)  = -1.6375d30
          end if
       end do
    end do
    

  end subroutine extract_multipole_range



  subroutine print_map_to_ascii(mapfile_in, maskfile, mapfile_out)
    implicit none

    character(len=512), intent(in) :: mapfile_in, maskfile, mapfile_out

    integer(i4b) :: npix, nside, ordering, ordering2, nmaps, i, j, c
    real(dp)     :: nullval, mu, sigma
    real(dp)     :: missval = -1.6375d30
    real(dp)     :: fmissval = 0.
    logical(lgt) :: anynull

    real(dp), allocatable, dimension(:,:) :: map, mask

    ! Read input map
    i = getsize_fits(mapfile_in, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile_in, map, npix, nmaps, nullval, anynull)

    ! Read input mask
    i = getsize_fits(maskfile, nside=nside, ordering=ordering2, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    if (ordering2 /= ordering) then
       do i = 1, nmaps
          if (ordering2 == 1) then
             call convert_ring2nest(nside, mask(:,i))
          else
             call convert_nest2ring(nside, mask(:,i))             
          end if
       end do
    end if

    open(48,file=trim(mapfile_out))
    mu = 0.d0
    c  = 1
    do i = 0, npix-1
       if (mask(i,c) /= 0. .and. mask(i,c) /= missval) then
          write(48,*) map(i,:)
          mu = mu + map(i,1)
       end if
    end do
    close(48)

    mu = mu / sum(mask(:,c))
    sigma = 0.d0
    do i = 0, npix-1
       if (mask(i,c) /= 0. .and. mask(i,c) /= missval) then
          sigma = sigma + (map(i,1)-mu)**2
       end if
    end do
    sigma = sqrt(sigma/(sum(mask(:,c))-1))

    write(*,*) 'Average outside mask = ', mu, ' +/- ', sigma

    deallocate(map)

  end subroutine print_map_to_ascii

  subroutine print_two_maps_to_ascii(mapfile1, mapfile2, maskfile, mapfile_out)
    implicit none

    character(len=512), intent(in) :: mapfile1, mapfile2, maskfile, mapfile_out

    integer(i4b) :: npix, nside, ordering, nmaps, i, j, ordering2
    real(dp)     :: nullval
    real(dp)     :: fmissval = 0.
    logical(lgt) :: anynull

    real(dp), allocatable, dimension(:,:) :: map1, map2, mask

    ! Read input map
    i = getsize_fits(mapfile1, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map1(0:npix-1,nmaps))
    call read_bintab(mapfile1, map1, npix, nmaps, nullval, anynull)

    ! Read input map
    i = getsize_fits(mapfile2, nside=nside, ordering=ordering2, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map2(0:npix-1,nmaps))
    call read_bintab(mapfile2, map2, npix, nmaps, nullval, anynull)
    if (ordering2 /= ordering) then
       do j = 1, nmaps
          if (ordering2 == 1) then
             call convert_ring2nest(nside, map2(:,j))
          else
             call convert_nest2ring(nside, map2(:,j))
          end if
       end do
    end if

    ! Read input mask
    i = getsize_fits(maskfile, nside=nside, ordering=ordering2, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    if (ordering2 /= ordering) then
       do j = 1, nmaps
          if (ordering2 == 1) then
             call convert_ring2nest(nside, mask(:,j))
          else
             call convert_nest2ring(nside, mask(:,j))
          end if
       end do
    end if

    open(48,file=trim(mapfile_out))
    do j = 1, 1 !nmaps
       do i = 0, npix-1
          if (mask(i,j) /= 0. .and. mask(i,j) /= 0.d0 .and. map1(i,j) /= -1.6375e30 .and. &
               & map2(i,j) /= -1.6375e30 .and. map1(i,j) /= -1.6375d30 .and. &
               & map2(i,j) /= -1.6375d30) write(48,*) map1(i,j), map2(i,j)
       end do
       write(48,*)
    end do
    close(48)

    deallocate(map1, map2)

  end subroutine print_two_maps_to_ascii

  subroutine print_isolatitude(mapfile_in, maskfile, outfile)
    implicit none

    character(len=512), intent(in) :: mapfile_in, maskfile, outfile

    integer(i4b) :: npix, nside, ordering, ordering2, nmaps, i, j, n, nlist, nring
    real(dp)     :: nullval, mu, sigma, dtheta
    real(dp)     :: fmissval = 0.
    logical(lgt) :: anynull
    integer(i4b), allocatable, dimension(:) :: pixlist

    real(dp), allocatable, dimension(:,:) :: map, mask

    ! Read input map
    i = getsize_fits(mapfile_in, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile_in, map, npix, nmaps, nullval, anynull)

    ! Read input mask
    i = getsize_fits(maskfile, nside=nside, ordering=ordering2, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    if (ordering2 /= ordering) then
       do i = 1, nmaps
          if (ordering2 == 1) then
             call convert_ring2nest(nside, mask(:,i))
          else
             call convert_nest2ring(nside, mask(:,i))             
          end if
       end do
    end if

    nring  = 35
    dtheta = pi/nring
    nring = nint(pi / dtheta) - 1
    allocate(pixlist(0:npix-1))
    open(58,file=trim(outfile))
    do i = 1, nring
       call query_strip(nside, (i-1)*dtheta, i*dtheta, pixlist, nlist)
       
       mu = 0.d0
       n  = 0
       do j = 0, nlist-1
          if (mask(pixlist(j),1) > 0.5d0) then
             mu = mu + map(pixlist(j),1)
             n  = n+1
          end if
       end do
       if (n > 0) then
          mu = mu / n
       else
          mu = 0.d0
       end if

       sigma = 0.d0
       do j = 0, nlist-1
          if (mask(pixlist(j),1) > 0.5d0) then
             sigma = sigma + (map(pixlist(j),1)-mu)**2
          end if
       end do
       sigma = sqrt(sigma/(n-1))

       write(58,*) (i-0.5d0)*dtheta*180.d0/pi, mu, sigma

    end do
    close(58)

  end subroutine print_isolatitude


  subroutine print_isolatitude_var(clfile, maskfile, outfile)
    implicit none

    character(len=512), intent(in) :: clfile, maskfile, outfile

    integer(i4b) :: npix, nside, ordering, ordering2, nmaps, i, j, n, nlist, nring
    real(dp)     :: nullval, mu, sigma, dtheta
    real(dp)     :: fmissval = 0.
    logical(lgt) :: anynull
    integer(i4b), allocatable, dimension(:) :: pixlist

    real(dp), allocatable, dimension(:,:) :: map, mask

!!$    ! Read input map
!!$    i = getsize_fits(mapfile_in, nside=nside, ordering=ordering, nmaps=nmaps)
!!$    npix = nside2npix(nside)
!!$    allocate(map(0:npix-1,nmaps))
!!$    call read_bintab(mapfile_in, map, npix, nmaps, nullval, anynull)

    ! Read input mask
    i = getsize_fits(maskfile, nside=nside, ordering=ordering2, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    if (ordering2 /= ordering) then
       do i = 1, nmaps
          if (ordering2 == 1) then
             call convert_ring2nest(nside, mask(:,i))
          else
             call convert_nest2ring(nside, mask(:,i))             
          end if
       end do
    end if

    nring  = 35
    dtheta = pi/nring
    nring = nint(pi / dtheta) - 1
    allocate(pixlist(0:npix-1))
    open(58,file=trim(outfile))
    do i = 1, nring
       call query_strip(nside, (i-1)*dtheta, i*dtheta, pixlist, nlist)
       
       mu = 0.d0
       n  = 0
       do j = 0, nlist-1
          if (mask(pixlist(j),1) > 0.5d0) then
             mu = mu + map(pixlist(j),1)
             n  = n+1
          end if
       end do
       if (n > 0) then
          mu = mu / n
       else
          mu = 0.d0
       end if

       sigma = 0.d0
       do j = 0, nlist-1
          if (mask(pixlist(j),1) > 0.5d0) then
             sigma = sigma + (map(pixlist(j),1)-mu)**2
          end if
       end do
       sigma = sqrt(sigma/(n-1))

       write(58,*) (i-0.5d0)*dtheta*180.d0/pi, mu, sigma

    end do
    close(58)

  end subroutine print_isolatitude_var


  subroutine summarize_detector_angles(filelist, s_max, nside, prefix)
    implicit none

    integer(i4b),     intent(in) :: s_max, nside
    character(len=*), intent(in) :: filelist, prefix

    integer(i4b)       :: unit1, unit2, npix, pix, s
    character(len=1)   :: s_string
    character(len=512) :: filename
    real(dp), allocatable, dimension(:,:,:) :: map
    real(dp), allocatable, dimension(:,:)   :: nhits
    real(dp),              dimension(6)     :: ang
    character(len=80), dimension(180)       :: header

    unit1 = 21
    unit2 = 22
    npix  = 12*nside**2

    allocate(map(0:npix-1,2,0:s_max))
    allocate(nhits(0:npix-1,1))

    map   = 0.d0
    nhits = 0.d0
    open(unit1,file=trim(filelist))
    do while (.true.)
       read(unit1,'(a)',end=990) filename
       open(unit2,file=trim(filename))
       write(*,*) 'Processing ', trim(filename)
       do while (.true.)
          read(unit2,*,end=991) ang
          ang = ang * pi / 180.d0

          ! Add first set
          call ang2pix_ring(nside, ang(2), ang(1), pix)
          nhits(pix,1) = nhits(pix,1) + 1.d0
          do s = 0, s_max
             map(pix,1,s) = map(pix,1,s) + cos(s*ang(3))
             map(pix,2,s) = map(pix,2,s) + sin(s*ang(3))
          end do

          ! Add second set
          call ang2pix_ring(nside, ang(5), ang(4), pix)
          nhits(pix,1) = nhits(pix,1) + 1.d0
          do s = 0, s_max
             map(pix,1,s) = map(pix,1,s) + cos(s*ang(6))
             map(pix,2,s) = map(pix,2,s) + sin(s*ang(6))
          end do
       end do
991    close(unit2)
    end do
990 close(unit1)

    filename = trim(prefix) // '_nhits.fits'
    call write_minimal_header(header, 'MAP', nside=nside, ordering='RING', polar=.false.)
    call write_result_map(filename, nside, 1, header, nhits)

    do s = 0, s_max
       write(*,*) 'c', s
       map(:,1,s) = map(:,1,s) / nhits(:,1)
       map(:,2,s) = map(:,2,s) / nhits(:,1)

       call int2string(s, s_string)
       filename = trim(prefix) // '_s' // s_string // '.fits'
       call write_minimal_header(header, 'MAP', nside=nside, ordering='RING')
       call write_result_map(filename, nside, 1, header, map(:,:,s))
    end do

    

  end subroutine summarize_detector_angles

  subroutine ud_grade_map_editor(filename_in, nside_out, filename_out, double_precision)
    implicit none

    character(len=*), intent(in) :: filename_in, filename_out
    integer(i4b),     intent(in) :: nside_out
    logical(lgt),     intent(in) :: double_precision

    real(dp)     :: nullval
    logical(lgt) :: anynull
    integer(i4b) :: i, j, nside_in, ordering, nmaps, npix_in, npix_out
    real(dp), allocatable, dimension(:,:) :: map_in, map_out
    character(len=80), dimension(180)     :: header

    if (iargc() /= 4) then
       write(*,*) 'Usage: map_editor ud_grade [input map] [nside_out] [output map]'
       stop
    end if

    i = getsize_fits(filename_in, nside=nside_in, ordering=ordering, nmaps=nmaps)
    npix_in  = nside2npix(nside_in)
    npix_out = nside2npix(nside_out)

    write(*,*) 'Input map nside = ', nside_in

    allocate(map_in(0:npix_in-1,nmaps))
    allocate(map_out(0:npix_out-1,nmaps))
    call read_bintab(filename_in, map_in, npix_in, nmaps, nullval, anynull, header=header)
    where (map_in < -1.6d30) 
       map_in = -1.6375d30
    end where

    if (ordering == 1) then
       call udgrade_ring(map_in, nside_in, map_out, nside_out)
    else
       call udgrade_nest(map_in, nside_in, map_out, nside_out)
    end if

    call write_result_map(filename_out, nside_out, ordering, header, map_out, double_precision)

    if (count(map_in==-1.6375d30) > 0) then
       do i=1,nmaps
          write(*,*) count(map_in(:,i)==-1.6375d30), ' missing pixels out of ', npix_in
       end do
    end if
    if (count(map_out==-1.6375d30) > 0) then
       do i=1,nmaps
          write(*,*) count(map_out(:,i)==-1.6375d30), ' missing pixels out of ', npix_out
       end do
    end if
    deallocate(map_in)
    deallocate(map_out)

  end subroutine ud_grade_map_editor

  subroutine ud_grade_rms_map_editor(filename_in, nside_out, filename_out, double_precision)
    implicit none

    character(len=*), intent(in) :: filename_in, filename_out
    integer(i4b),     intent(in) :: nside_out
    logical(lgt),     intent(in) :: double_precision

    real(dp)     :: nullval, missval
    logical(lgt) :: anynull
    integer(i4b) :: i, j, nside_in, ordering, nmaps, npix_in, npix_out
    real(dp), allocatable, dimension(:,:) :: map_in, map_out
    character(len=80), dimension(180)     :: header

    missval = -1.6375d30

    if (iargc() /= 4) then
       write(*,*) 'Usage: map_editor ud_grade_rms [input map] [nside_out] [output map]'
       stop
    end if

    i = getsize_fits(filename_in, nside=nside_in, ordering=ordering, nmaps=nmaps)
    npix_in  = nside2npix(nside_in)
    npix_out = nside2npix(nside_out)

    write(*,*) 'Input map nside = ', nside_in

    allocate(map_in(0:npix_in-1,nmaps))
    allocate(map_out(0:npix_out-1,nmaps))
    call read_bintab(filename_in, map_in, npix_in, nmaps, nullval, anynull, header=header)

    if (nside_in == nside_out) then
       map_out = map_in
    else
       where (map_in < -1.6d30) 
          map_in = -1.6375d30
       elsewhere
          map_in = map_in*map_in !get variance
       end where

       if (ordering == 1) then
          call udgrade_ring(map_in, nside_in, map_out, nside_out)
       else
          call udgrade_nest(map_in, nside_in, map_out, nside_out)
       end if

       where (abs((map_out-missval)/missval) > 1.d-5)
          map_out = sqrt(map_out)*(nside_out*1.d0/nside_in) !get RMS from average variance and scale with nside ratio
       end where
    end if

    call write_result_map(filename_out, nside_out, ordering, header, map_out, double_precision)

    if (count(map_in==-1.6375d30) > 0) then
       do i=1,nmaps
          write(*,*) count(map_in(:,i)==-1.6375d30), ' missing pixels out of ', npix_in
       end do
    end if
    if (count(map_out==-1.6375d30) > 0) then
       do i=1,nmaps
          write(*,*) count(map_out(:,i)==-1.6375d30), ' missing pixels out of ', npix_out
       end do
    end if
    deallocate(map_in)
    deallocate(map_out)

  end subroutine ud_grade_rms_map_editor


  subroutine read_covmatrix(unit, filename, ordering, polarization, matrix, inv, n)

    integer(i4b),                          intent(in)  :: unit
    integer(i8b),                          intent(out) :: n
    integer(i4b),                          intent(out) :: ordering, polarization
    character(len=*),                      intent(in)  :: filename
    logical(lgt),                          intent(out) :: inv
    real(dp), allocatable, dimension(:,:), intent(out) :: matrix

    integer(i4b) :: n_in, i

    ! polarization = 1  => I only
    ! polarization = 2  => QU only
    ! polarization = 3  => IQU

    open(unit, file=trim(filename), form='unformatted')
    read(unit) n_in
    n = n_in
    write(*,*) n, '= n'
    read(unit) ordering
    write(*,*) ordering, '= ordering'
    read(unit) polarization
    write(*,*) polarization, '= polarisation'
    allocate(matrix(n,n))
    do i = 1, n
       read(unit) matrix(:,i)
    end do
    read(unit) inv
    close(unit)

  end subroutine read_covmatrix


  subroutine make_co_region_map(mapname, nside, ordering, reg)
    implicit none

    character(len=*),                         intent(in)  :: mapname
    integer(i4b),                             intent(out) :: nside, ordering
    real(dp),         pointer, dimension(:,:)             :: reg

    integer(i4b) :: npix, nmaps, i, j, k, lmax, nreg, p(1), listpix(0:100000), nlist, r
    real(dp)     :: nullval, tot, fmissval, threshold, radius, vec(3), threshold0, dist, dist0
    logical(lgt) :: anynull
    real(dp),     allocatable, dimension(:,:)       :: map, map_seed, seeds

    threshold0 = 1.d0 ! K km/s
    radius     = 10.d0 * pi/180.d0 ! Distance between seeds in radians

    ! Read input map
    i = getsize_fits(mapname, nside=nside, ordering=ordering)
    npix = nside2npix(nside)
    nmaps = 1
    allocate(map(0:npix-1,nmaps), reg(0:npix-1,nmaps), map_seed(0:npix-1,nmaps), seeds(10000,3))
    call read_bintab(mapname, map, npix, nmaps, nullval, anynull)
    if (ordering == 2) call convert_nest2ring(nside, map(:,1))
    ordering = 1

    nreg = 0
    reg     = -1.6375d30

    ! Start by setting all low-amplitude pixels to region 1
    nreg = nreg+1
    where (map < threshold0) 
       reg = nreg
       map = -1.d30
    end where

    ! Set up a grid of seed points
    map_seed = map
    do while (any (map_seed >= 0.d0))
       p = maxloc(map_seed(:,1))
       call pix2vec_ring(nside, p(1), vec)
       call query_disc(nside, vec, 2.d0*radius, listpix, nlist)       

       nreg = nreg+1
       seeds(nreg,:) = vec
       reg(p,1)  = nreg
       map_seed(listpix(0:nlist-1),1) = -1.d30
       write(*,*) 'p = ', p, ', nreg = ', nreg
    end do

    ! Put remaining pixels in closest seed
    do i = 0, npix-1
       if (reg(i,1) > 0.d0) cycle
       call pix2vec_ring(nside, i, vec)

       dist0 = 1.d30
       do j = 2, nreg
          call angdist(vec, seeds(j,:), dist)
          if (dist < dist0) then
             reg(i,1) = j
             dist0    = dist
          end if
       end do
    end do

    deallocate(map, map_seed)

  end subroutine make_co_region_map



  ! Parallel transport polarization angles from a given map to centers
  ! of a lower resolution map at nside_out. 
  ! Note: Maps must be in nested format
  subroutine qu_transport_map(nside_out, map)
    implicit none

    integer(i4b),                   intent(in)    :: nside_out
    real(dp),     dimension(0:,1:), intent(inout) :: map

    integer(i4b) :: i, j, q
    real(dp)     :: cos2psi, sin2psi, m(2)
    integer(i4b), save :: nside, npix, nmaps
    real(dp), allocatable, dimension(:,:,:), save :: vec

    ! Precompute pixel vectors
    if (.not. allocated(vec)) then
       npix  = size(map,1)
       nmaps = size(map,2)
       nside = npix2nside(npix)
       q     = (nside/nside_out)**2
       if (nmaps /= 3) then
          write(*,*) 'partrans -- nmaps must be equal to 3'
          stop
       end if
       allocate(vec(0:npix-1,3,2))
       do i = 0, npix-1
          call pix2vec_nest(nside_out, i/q, vec(i,:,1)) ! Low-res pixel centers
          call pix2vec_nest(nside,     i,   vec(i,:,2)) ! High-res pixel centers
       end do
    end if

    ! Perform the parallel transport
    do i = 0, npix-1
       call qu_transport_rot(vec(i,:,1), vec(i,:,2), cos2psi, sin2psi)
       m        = map(i,2:3)
       map(i,2) =  cos2psi*m(1) - sin2psi*m(2)
       map(i,3) =  sin2psi*m(1) + cos2psi*m(2)
    end do

  end subroutine qu_transport_map

  subroutine merge_maps(double_precision)
    implicit none

    logical(lgt), intent(in) :: double_precision

    integer(i4b)       :: npix
    integer(i4b)       :: i, j, l, m
    real(dp)           :: nullval
    real(dp)           :: sigma_sq, nullval_dp, fwhm_trans, threshold, cl
    logical(lgt)       :: anynull
    character(len=512) :: beamfile1, beamfile2, infile1, infile2, outfile, partext, prefix, maskfile
    integer(i4b)       :: ltrans_start, ltrans_stop, ordering, nmaps
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    integer(i4b) :: nside1, nside2, npix1, npix2

    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     pointer,     dimension(:,:)   :: weights1, weights2
    real(dp),     allocatable, dimension(:,:)   :: beam1, beam2, map1, map2, map, map_out, w_l
    character(len=80),         dimension(1:180) :: header
    real(dp),                  dimension(2)     :: zbounds = 0.d0
    real(dp),     pointer,     dimension(:,:)   :: pixwin1, pixwin2

    if (iargc() /= 10) then
       write(*,*) 'Usage: map_editor merge_maps [map_lowl] [beam_lowl] [map_highl] [beam_highl]'
       write(*,*) '                    [l_trans_start] [l_trans_stop] [threshold] [outmap] [prefix]'
       stop
    end if

    call getarg(2,infile1)
    call getarg(3,beamfile1)
    call getarg(4,infile2)
    call getarg(5,beamfile2)
    call getarg(6,partext)
    read(partext,*) ltrans_start
    call getarg(7,partext)
    read(partext,*) ltrans_stop
    call getarg(8,partext)
    read(partext,*) threshold
    call getarg(9,outfile)
    call getarg(10,prefix)

    ! Read input data set 1
    i = getsize_fits(infile1, nside=nside1, ordering=ordering, nmaps=nmaps)
    npix1 = nside2npix(nside1)
    allocate(map1(0:npix1-1,nmaps))
    call read_bintab(infile1, map1, npix1, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside1, map1(0:npix-1,i))
       end do
    end if
    call read_ringweights(nside1, weights1)
    allocate(beam1(0:ltrans_stop, nmaps))
    call generate_beam(0.d0, ltrans_stop, beam1, beamfile1)
    call read_pixwin(nside1, nmaps, pixwin1)

    ! Read input data set 2
    i = getsize_fits(infile2, nside=nside2, ordering=ordering, nmaps=nmaps)
    npix2 = nside2npix(nside2)
    allocate(map2(0:npix2-1,nmaps), map_out(0:npix2-1,nmaps))
    call read_bintab(infile2, map2, npix2, nmaps, nullval, anynull)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside2, map2(0:npix2-1,i))
       end do
    end if
    call read_ringweights(nside2, weights2)
    allocate(beam2(0:ltrans_stop, nmaps))
    call generate_beam(0.d0, ltrans_stop, beam2, beamfile2)
    call read_pixwin(nside2, nmaps, pixwin2)

    ! Threshold maps
    if (.false. .and. threshold > 0.d0) then
       map1 = max(min(map1,threshold),-threshold)
       map2 = max(min(map2,threshold),-threshold)
    end if

    ! Compute transition beam, and effective weight function
    allocate(w_l(0:ltrans_stop,nmaps))
    do i = 1, nmaps
       do l = 0, ltrans_stop
          ! Gaussian smoothing
          w_l(l,i) = 1.d0 
          if (l > ltrans_start) then
             ! Cosine apodization
             w_l(l,i) = w_l(l,i) * 0.5d0*(1.d0-cos(pi*real(ltrans_stop-l,dp)/real(ltrans_stop-ltrans_start,dp)))          
          end if
       end do
    end do

    ! Compute the spherical harmonics transform
    allocate(alms(nmaps, 0:ltrans_stop, 0:ltrans_stop))
    allocate(map(0:npix1-1,nmaps))
    map = map1
    !if (threshold > 0) map = max(min(map,threshold),-threshold)
    if (nmaps == 1) then
       call map2alm(nside1, ltrans_stop, ltrans_stop, map(:,1), alms, zbounds, weights1)
    else
       call map2alm(nside1, ltrans_stop, ltrans_stop, map, alms, zbounds, weights1)
    end if
    deallocate(map)

    ! Deconvolve old beam, convolve with new
    do i = 1, nmaps
       do l = 0, ltrans_stop
          if (beam1(l,i)*pixwin1(l,i) > 0.d0) then
             alms(i,l,0:l) = alms(i,l,0:l) / (beam1(l,i)*pixwin1(l,i)) * (beam2(l,i)*pixwin2(l,i))
          end if
       end do
    end do

    open(58,file=trim(prefix)//'_cls1.dat')
    do l = 0, ltrans_stop
       cl = abs(alms(1,l,0))**2
       do m = 1, l
          cl = cl + 2.d0*abs(alms(1,l,m))**2
       end do
       write(58,*) l, cl/real(2*l+1,dp) * l*(l+1)/2/pi
    end do
    close(58)

    ! Cosine apodization
    open(58,file=trim(prefix)//'_apod.dat')
    do i = 1, nmaps
       do l = 0, ltrans_stop
          alms(i,l,0:l) = alms(i,l,0:l) * w_l(l,i)
          if (i == 1) write(58,*) l, w_l(l,i), 1.d0-w_l(l,i)
       end do
    end do
    close(58)

    ! Compute the inverse spherical harmonics
    if (nmaps == 1) then
       call alm2map(nside2, ltrans_stop, ltrans_stop, alms, map_out(:,1))
    else
       call alm2map(nside2, ltrans_stop, ltrans_stop, alms, map_out)
    end if
    deallocate(map1)
    allocate(map1(0:npix2-1,nmaps))
    map1 = map_out

    call write_result_map(trim(prefix)//'_comp1.fits', nside2, ordering, header, map_out, double_precision)


    ! Add second map
    map_out = map_out + map2
    !map_out = map2


    ! Compute the spherical harmonics transform
    allocate(map(0:npix2-1,nmaps))
    map = map2
    !if (threshold > 0) map = max(min(map,threshold),-threshold)
    if (nmaps == 1) then
       call map2alm(nside2, ltrans_stop, ltrans_stop, map(:,1), alms, zbounds, weights2)
    else
       call map2alm(nside2, ltrans_stop, ltrans_stop, map, alms, zbounds, weights2)
    end if
    deallocate(map)

    open(58,file=trim(prefix)//'_cls2.dat')
    do l = 0, ltrans_stop
       cl = abs(alms(1,l,0))**2
       do m = 1, l
          cl = cl + 2.d0*abs(alms(1,l,m))**2
       end do
       write(58,*) l, cl/real(2*l+1,dp) * l*(l+1)/2/pi
    end do
    close(58)

    ! Cosine apodization
    do i = 1, nmaps
       do l = 0, ltrans_stop
          alms(i,l,0:l) = alms(i,l,0:l) * w_l(l,i)
       end do
    end do

    ! Compute the inverse spherical harmonics
    if (nmaps == 1) then
       call alm2map(nside2, ltrans_stop, ltrans_stop, alms, map2(:,1))
    else
       call alm2map(nside2, ltrans_stop, ltrans_stop, alms, map2)
    end if

    ! Subtract large scales to avoid double counting
    map_out = map_out - map2

    call write_result_map(trim(prefix)//'_comp2.fits', nside2, ordering, header, map_out-map1, double_precision)

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_ring2nest(nside2, map_out(:,i))
       end do
    end if

    call write_result_map(outfile, nside2, ordering, header, map_out, double_precision)

    if (nmaps == 1) then
       call map2alm(nside2, ltrans_stop, ltrans_stop, map_out(:,1), alms, zbounds, weights2)
    else
       call map2alm(nside2, ltrans_stop, ltrans_stop, map_out, alms, zbounds, weights2)
    end if

    open(58,file=trim(prefix)//'_cls.dat')
    do l = 0, ltrans_stop
       cl = abs(alms(1,l,0))**2
       do m = 1, l
          cl = cl + 2.d0*abs(alms(1,l,m))**2
       end do
       write(58,*) l, cl/real(2*l+1,dp) * l*(l+1)/2/pi
    end do
    close(58)

  end subroutine merge_maps

  subroutine rescale_dust_amplitude_map
    implicit none

    integer(i4b)       :: i, j, k, l, m, n1, n2, lmax
    real(dp)           :: nullval, old_nu, new_nu,h_div_k
    real(dp)           :: nullval_dp, missval
    logical(lgt)       :: anynull
    character(len=512) :: infile, infile2, outfile, partext
    integer(i4b)       :: ordering, nmaps, nmaps2, nmaps3, ordering2
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    integer(i4b) :: nside, npix, nside2

    real(dp),     allocatable, dimension(:,:)   :: map, outmap, T_map, b_map, x_map, mask
    character(len=80),         dimension(1:180) :: header, header2

    missval = -1.6375d30
    h_div_k = 6.626069d-34/1.38065d-23*1.d9  !h/k_B * 1e9 due to GHz


    if (iargc() /= 7) then
       write(*,*) 'Usage: map_editor rescale_dust_amplitude_map [ampl_inmap] [Td_map] [beta_map]'
       write(*,*) '             [old nu_ref (GHz)] [new nu_ref (GHz)] [ampl_outmap]'
       stop
    end if

    call getarg(2,infile)
    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)

    allocate(map(0:npix-1,nmaps), outmap(0:npix-1,nmaps), mask(0:npix-1,nmaps), x_map(0:npix-1,nmaps))
    call read_bintab(infile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(0:npix-1,i))
       end do
    end if

    outmap=0.d0
    ! creating mask:
    mask=1.0
    do j=1,nmaps
       do i=0,npix-1
          if (abs((map(i,j)-missval)/missval) < 1.d-6) then
             mask(i,j)=0.d0
          end if
       end do
    end do

    call getarg(3,infile2)
    i = getsize_fits(infile2, nside=nside2, ordering=ordering2, nmaps=nmaps2)
    if (nside /= nside2) then
       write(*,*) "Temperature map must have the same nside as amplitude map"
       write(*,*) "nside ampl:", nside, "   nside T:", nside2
       stop
    end if
    allocate(T_map(0:npix-1,nmaps2))
    call read_bintab(infile2, T_map, npix, nmaps2, nullval, anynull, header=header)
    
    if (ordering2 == 2) then
       do i = 1, nmaps2
          call convert_nest2ring(nside, T_map(0:npix-1,i))
       end do
    end if

    if (nmaps2 > nmaps) nmaps2=nmaps
    do j=1,nmaps2
       do i=0,npix-1
          if (abs((T_map(i,j)-missval)/missval) < 1.d-6) then
             mask(i,j)=0.d0
             x_map(i,j)=0.d0
          else
             x_map(i,j)=h_div_k/t_map(i,j)
          end if
       end do
    end do

    call getarg(4,infile2)
    i = getsize_fits(infile2, nside=nside2, ordering=ordering2, nmaps=nmaps3)
    if (nside /= nside2) then
       write(*,*) "Beta map must have the same nside as amplitude map"
       write(*,*) "nside ampl:", nside, "   nside beta:", nside2
       stop
    end if
    allocate(b_map(0:npix-1,nmaps3))
    call read_bintab(infile2, b_map, npix, nmaps3, nullval, anynull, header=header)
    if (ordering2 == 2) then
       do i = 1, nmaps3
          call convert_nest2ring(nside, b_map(0:npix-1,i))
       end do
    end if

    if (nmaps3 > nmaps) nmaps3=nmaps
    do j=1,nmaps3
       do i=0,npix-1
          if (abs((b_map(i,j)-missval)/missval) < 1.d-6) then
             mask(i,j)=0.d0
          end if
       end do
    end do

    call getarg(5,partext)
    read(partext,*) old_nu ! In GHz
    call getarg(6,partext)
    read(partext,*) new_nu ! In GHz
    call getarg(7,outfile)

    if (nmaps2 > nmaps3) nmaps2=nmaps3 !we find the smallest nmaps
    do j=1,nmaps2
       do i=0,npix-1
          if (mask(i,j) < 0.5) then
             outmap(i,j) = missval
          else
             outmap(i,j) = map(i,j)*(exp(x_map(i,j)*old_nu)-1)/(exp(x_map(i,j)*new_nu)-1)*(new_nu/old_nu)**(b_map(i,j)+1.d0)
          end if
       end do
    end do

    if (ordering == 2) then
       do i = 1, nmaps2
          call convert_ring2nest(nside, outmap(0:npix-1,i))
       end do
    end if
    
    ! Output result map
    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=nmaps2==3)
    call write_result_map(outfile, nside, ordering, header, outmap(:,:nmaps2))
    deallocate(map,outmap,mask,T_map,b_map,x_map)
  end subroutine rescale_dust_amplitude_map

  subroutine median_filter_holes
    implicit none

    integer(i4b)       :: i, j, k, l, m, n1, n2, lmax, j_m
    real(dp)           :: nullval, mean1, mean2
    real(dp)           :: nullval_dp, radius, vec(3), threshold
    logical(lgt)       :: anynull, output_mask, mask_present
    character(len=512) :: infile, outfile, partext, maskfile
    integer(i4b)       :: ordering, nmaps, ordering2, nmaps2
    integer(i4b)       :: nside, npix, nside2
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle


    real(dp),     pointer,     dimension(:,:)   :: weights
    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     allocatable, dimension(:,:)   :: map, outmap, mask, mask2, mask3, inmask
    character(len=80),         dimension(1:180) :: header, header2
    integer(i4b), allocatable, dimension(:)     :: listpix1, listpix2

    if (iargc() < 5) then
       write(*,*) ''
       write(*,*) 'Usage: map_editor median_filter_source_holes [inmap] [radius arcmin] [threshold] [outmap] (options)'
       write(*,*) ''
       write(*,*) 'Options:'
       write(*,*) '  -mask [mask]   (a mask to filter inside, specific filtering)'
       write(*,*) '  -output_mask   (also output the two-radii mask and the apodization mask)'
       write(*,*) ''
       stop
    end if

    output_mask=.false.
    mask_present=.false.

    call getarg(2,infile)
    call getarg(3,partext)
    read(partext,*) radius ! In arcmin
    call getarg(4,partext)
    read(partext,*) threshold ! In arcmin
    call getarg(5,outfile)
    radius = radius/60.d0 * pi/180.d0

    if (iargc() > 5) then
       do i = 6,iargc()
          call getarg(i,partext)
          if (trim(partext) == '-mask') then
             call getarg(i+1,maskfile)
             mask_present=.true.
          elseif (trim(partext) == '-output_mask') then
             output_mask = .true.
          else
             call getarg(i-1,partext)
             if (trim(partext) /= '-mask') then
                call getarg(i,partext)
                write(*,*) 'Unknown option: '//trim(partext)
                stop
             end if
          end if
       end do
    end if

    ! Read input data set 1
    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)

    if (mask_present) then
       i = getsize_fits(maskfile, nside=nside2, ordering=ordering2, nmaps=nmaps2)
       if (nside /= nside2) then
          write(*,*) 'Error: nside of input map and mask is not equal. Exiting.'
          stop
       end if
    end if

    allocate(map(0:npix-1,nmaps), outmap(0:npix-1,nmaps), mask(0:npix-1,nmaps), mask2(0:npix-1,nmaps), mask3(0:npix-1,nmaps))
    call read_bintab(infile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(0:npix-1,i))
       end do
    end if
    if (mask_present) then
       allocate(inmask(0:npix-1,nmaps2))
       call read_bintab(maskfile, inmask, npix, nmaps2, nullval, anynull, header=header2)
       if (ordering2 == 2) then
          do i = 1, nmaps2
             call convert_nest2ring(nside, inmask(0:npix-1,i))
          end do
       end if
    else
       allocate(inmask(0:npix-1,nmaps))
       inmask=1.d0
    end if

    ! Loop through all pixels (inside input mask if given), searching for holes
    ! If input mask is given and nmaps of mask is less than that of the input map, 
    ! use map #1 of the mask for the overflowing maps.
    allocate(listpix1(0:1000000), listpix2(0:1000000))
    do j = 1, nmaps
       if (j > nmaps2) then
          j_m = 1
       else
          j_m = j
       end if
       do i = 0, npix-1
          if (mod(i,1000000) == 0) write(*,*) j, i, npix
          if (inmask(i,j_m) > 0.5) then
             call pix2vec_ring(nside, i, vec)
             call query_disc(nside, vec, radius, listpix1, n1)
             mean1 = mean(map(listpix1(0:n1-1),j))
             if (map(i,j) < 1.d-4 * mean1 .and. mean1 > threshold) then
                mask(i,j)   = 1.d0 
             else
                !outmap(i,j) = map(i,j)
                mask(i,j)   = 0.d0 
             end if
          end if
       end do
    end do

    ! Expand relevant region by 4 radii 
    ! First map 2*radii, then median filter those over 2 radii.
    mask2 = 0.d0
    mask3 = 0.d0
    do j = 1, nmaps
       do i = 0, npix-1
          if (mask(i,j) == 1.d0) then
             call pix2vec_ring(nside, i, vec)
             call query_disc(nside, vec, 2.d0*radius, listpix1, n1)             
             mask2(listpix1(0:n1-1),j) = 1.d0 ! mask with all pixels inside 2 radii of holes
             call query_disc(nside, vec, radius, listpix1, n1)             
             mask3(listpix1(0:n1-1),j) = 1.d0 ! mask with all pixels inside 1 radii of holes
          end if
       end do
    end do

    ! Median filter selected pixels (inside 2 radii of holes) over 2 radii
    ! 'outmap' is equal to input outside 2 radii of the holes, 
    ! inside it is the median of 2 radii of the pixel.
    do j = 1, nmaps
       do i = 0, npix-1
          if (mask2(i,j) == 1.d0) then
             call pix2vec_ring(nside, i, vec)
             call query_disc(nside, vec, 2.d0*radius, listpix1, n1)
             outmap(i,j) = median(map(listpix1(0:n1-1),j))
          else
             outmap(i,j) = map(i,j)
          end if
       end do
    end do

    ! Apodize mask, i.e. create a smoothed mask around source holes (the single radius mask)
    lmax = 3*nside
    mask = mask3
    call read_ringweights(nside, weights)    
    allocate(alms(1,0:lmax,0:lmax))
    do j = 1, nmaps
       call map2alm(nside, lmax, lmax, mask(:,j), alms, [0.d0,0.d0], weights)
       do l = 0, lmax
          alms(1,l,0:l) = alms(1,l,0:l) * &
               & exp(-0.5d0*l*(l+1.d0)*(radius*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
       end do
       call alm2map(nside, lmax, lmax, alms, mask(:,j))
    end do

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_ring2nest(nside, outmap(0:npix-1,i))
          call convert_ring2nest(nside, mask(0:npix-1,i))
          call convert_ring2nest(nside, mask2(0:npix-1,i))
       end do
    end if

    ! Weight pixels (around source holes)
    do j = 1, nmaps
       do i = 0, npix-1
          if (mask(i,1) > 1.d-2) then
             outmap(i,j) = (1.d0-mask(i,1))*map(i,j) + mask(i,1) * outmap(i,j)
          else
             outmap(i,j) = map(i,j)
          end if
       end do
    end do

    ! Output result map
    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=nmaps==3)
    call write_result_map('apomask_'//trim(outfile), nside, ordering, header, mask)
    call write_result_map('mask_'//trim(outfile), nside, ordering, header, mask2)
    call write_result_map(outfile, nside, ordering, header, outmap)

  end subroutine median_filter_holes

  subroutine median_filter_specific_value
    implicit none

    integer(i4b)       :: i, j, k, l, m, n1, n2, lmax, pixcount
    real(dp)           :: nullval, mean1
    real(dp)           :: nullval_dp, radius, vec(3), threshold
    logical(lgt)       :: anynull
    character(len=512) :: infile, outfile, partext
    integer(i4b)       :: ordering, nmaps
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    integer(i4b) :: nside, npix

    real(dp),     allocatable, dimension(:,:)   :: map, outmap
    character(len=80),         dimension(1:180) :: header
    integer(i4b), allocatable, dimension(:)     :: listpix1
    real(dp),     allocatable, dimension(:)     :: listmap1

    if (iargc() /= 5) then
       write(*,*) 'Usage: map_editor median_filter_specific_value [inmap] [radius (arcmin)]'
       write(*,*) '                                               [value] [outmap]'
       stop
    end if

    call getarg(2,infile)
    call getarg(3,partext)
    read(partext,*) radius ! In arcmin
    call getarg(4,partext)
    read(partext,*) threshold ! In arcmin
    call getarg(5,outfile)
    radius = radius/60.d0 * pi/180.d0

    ! Read input data set 1
    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps), outmap(0:npix-1,nmaps))
    call read_bintab(infile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(0:npix-1,i))
       end do
    end if

    ! Loop through all pixels, searching for holes
    allocate(listpix1(0:1000000),listmap1(0:1000000))
    write(*,*)
    do j = 1, nmaps
       pixcount=0
       do i = 0, npix-1
          if (threshold==0.0) then
    
             if (map(i,j) == threshold) then
                call pix2vec_ring(nside, i, vec)
                call query_disc(nside, vec, radius, listpix1, n1)
                m=-1
                do k = 0,n1-1
                   if ((map(listpix1(k),j)==threshold) > 1d-6) then
                      m = m + 1
                      listmap1(m)=map(listpix1(k),j)
                   end if
                end do
                if (m==-1) then
                   write(*,*) "No values outside the input value found, increase radius"
                   stop
                end if
                mean1 = median(listmap1(0:m))
                outmap(i,j)=mean1

                pixcount = pixcount+1
             else
                outmap(i,j)=map(i,j)
             end if

          else
    
             if (abs((map(i,j)-threshold)/threshold) < 1d-6) then
                call pix2vec_ring(nside, i, vec)
                call query_disc(nside, vec, radius, listpix1, n1)
                m=-1
                do k = 0,n1-1
                   if (abs((map(listpix1(k),j)-threshold)/threshold) > 1d-6) then
                      m = m + 1
                      listmap1(m)=map(listpix1(k),j)
                   end if
                end do
                if (m==-1) then
                   write(*,*) "No values outside the input value found, increase radius"
                   stop
                end if
                mean1 = median(listmap1(0:m))
                outmap(i,j)=mean1

                pixcount = pixcount+1
             else
                outmap(i,j)=map(i,j)
             end if

          end if
                 
       end do
       write(*,*) "    Map:",j,"Pixels with specific value:",pixcount
    end do
    write(*,*)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_ring2nest(nside, outmap(0:npix-1,i))
          call convert_ring2nest(nside, map(0:npix-1,i))
       end do
    end if


    ! Output result map
    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=nmaps==3)
    call write_result_map(outfile, nside, ordering, header, outmap)

    deallocate(map, outmap)
  end subroutine median_filter_specific_value


  subroutine median_filter_misspix
    implicit none

    integer(i4b)       :: i, j, k, l, m, n1, n2, lmax, nummiss, r, nlist
    real(dp)           :: nullval, mean1, mean2, theta, phi, missval
    real(dp)           :: nullval_dp, radius, vec(3), threshold
    logical(lgt)       :: anynull, anymiss, dpcfix
    character(len=512) :: infile, outfile, partext
    integer(i4b)       :: ordering, nmaps, nref
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    integer(i4b) :: nside, npix, r_fill
    real(dp)     :: rf, tot

    real(dp),     allocatable, dimension(:,:)   :: map, outmap
    character(len=80),         dimension(1:180) :: header
    real(dp),     allocatable, dimension(:)     :: buffer
    integer(i4b), allocatable, dimension(:)     :: listpix

    if (iargc() < 4) then
       write(*,*) '  Usage: map_editor median_filter_misspix [map_in] [radius in pixels]'
       write(*,*) '                                          [map_out] [[dpc]]'
       stop
    end if

    call getarg(2,infile)
    call getarg(3,partext)
    read(partext,*) r_fill ! In pixels
    call getarg(4,outfile)
    dpcfix = .false.
    if (iargc() > 4) then
       call getarg(5,partext)
       read(partext,*) dpcfix
    end if

    missval = -1.6375d30 !missing pixel value
    ! Read input data set 1
    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps), outmap(0:npix-1,nmaps))
    call read_bintab(infile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(0:npix-1,i))
       end do
    end if
    write(*,*) "LFI dpc =", dpcfix, "nmaps =", nmaps

    ! Loop through all pixels, searching for holes
    allocate(listpix(0:1000000), buffer(1000000))
    anymiss = .false.
    do j = 1, nmaps
       nref = j
       if (dpcfix) then
          if (nmaps==3 .and. j==3) then
             nref = 1
          else if (nmaps==10 .and. j==5) then
             nref = 1
          else if (nmaps==10 .and. j==8) then
             nref = 2
          else if (nmaps==10 .and. j==10) then
             nref = 3
          end if
       end if
       do i = 0, npix-1
          if (abs((map(i,nref)-missval)/missval) < 1d-5) then  !pixel has missing pixel value
             anymiss = .true.
             call pix2vec_ring(nside, i, vec)
             do r = 1, r_fill
                rf = (r+0.2d0) * sqrt(4.d0*pi/npix)
                call query_disc(nside, vec, rf, listpix, nlist)
                tot = 0.d0
                m   = 0.d0
                do k = 0, nlist-1
                   if (abs((map(listpix(k),nref)-missval)/missval) > 1d-5) then !pixel does not have missing pixel value
                      tot = tot + map(listpix(k),j)
                      m   = m   + 1
                      buffer(m) = map(listpix(k),j)
                   end if
                end do
                if (m > 4) then    ! more than 4 pixels that is not missing pixels 
                   !outmap(i,j) = tot / m
                   outmap(i,j) = median(buffer(1:m))  !set missing pix to median.
                   exit
                end if
             end do

             if (r > r_fill) then
                call pix2ang_ring(nside, i, theta, phi)
                write(*,*) 'Error: Filter did not return valid pixel; increase radius'
!                write(*,*) '     p         = ', i
!                write(*,*) '     theta     = ', 90.d0-theta*180/pi
!                write(*,*) '     phi       = ', phi*180/pi
!                write(*,*) '     radius    = ', r_fill
!                write(*,*) '     neighbors = ', listpix(0:nlist-1)
!                write(*,*) '     vals      = ', map(listpix(0:nlist-1),j)
                stop
             end if
          else
             outmap(i,j) = map(i,j)
          end if
       end do
    end do

!    if (anymiss) then
       do i = 1, nmaps
          write(*,*) i, count(abs(map(:,i)) > 1e30), ' missing pixels out of ', npix
          if (ordering==2) call convert_ring2nest(nside, outmap(0:npix-1,i))
       end do
       ! Output result map
       call write_result_map(outfile, nside, ordering, header, outmap)
!    end if

  end subroutine median_filter_misspix


  subroutine median_filter_misspix_neighbor
    !get the (n'th pixel radius) neighboring pixels of missing pixels and give 
    !the pixel value equal to the median of the none missing neighbors.
    !do so recursively until all is filled up. 
    implicit none

    integer(i4b)       :: i, j, k, l, m, n1,n2, lmax, nummiss, r, nlist
    real(dp)           :: nullval, mean1, mean2, theta, phi, missval
    real(dp)           :: nullval_dp, radius, vec(3), threshold, miss_sum
    logical(lgt)       :: anynull, anymiss, dpcfix
    character(len=512) :: infile, outfile, partext
    integer(i4b)       :: ordering, nmaps, nref
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    integer(i4b) :: nside, npix, r_fill
    real(dp)     :: rf, tot

    real(dp),     allocatable, dimension(:,:)   :: map, outmap, tempmap, mask
    character(len=80),         dimension(1:180) :: header
    real(dp),     allocatable, dimension(:)     :: buffer
    integer(i4b), allocatable, dimension(:)     :: listpix

    if (iargc() < 3) then
       write(*,*) '  Usage: map_editor median_filter_misspix_neigbor [map_in] [map_out]'
       stop
    end if

    call getarg(2,infile)
    call getarg(3,outfile)

    missval = -1.6375d30 !missing pixel value
    ! Read input data set 1
    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps), outmap(0:npix-1,nmaps))
    allocate(tempmap(0:npix-1,nmaps),mask(0:npix-1,nmaps))
    call read_bintab(infile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 1) then
       do i = 1, nmaps
          call convert_ring2nest(nside, map(0:npix-1,i))
       end do
    end if

    m=0
    allocate(listpix(8),buffer(8))
    anymiss=.true.
    do while (anymiss)
       !create a mask with misspix, using where
       where (abs((missval-map)/missval) < 1.0d-5)
          mask=1.0
       elsewhere
          mask=0.0
       end where
       
       !sum mask, if sum = 0, anymiss=.false. and coutinue
       miss_sum=sum(mask)
       if (miss_sum < 0.5) then
          anymiss=.false.
          cycle
       end if

       m=m+1
       write(*,*) 'Missing pixels, cycle nr.',m
       tempmap=map
       !pixel by pixel
       do j = 1,nmaps
          do i = 0,npix-1
             if (abs((missval-map(i,j))/missval) < 1.0d-5) then
                
                call neighbours_nest(nside, i, listpix,n1)
                n2=0
                buffer=0
                do k = 1,n1
                   if (abs((missval-map(listpix(k),j))/missval) > 1.0d-5) then
                      n2=n2+1
                      buffer(n2)=map(listpix(k),j)
                   end if
                end do
                if (n2==0) then
                   tempmap(i,j)=missval !neighbors are all missval
                else
                   tempmap(i,j)=median(buffer(1:n2))
                end if
             end if
          end do
       end do
       map=tempmap
    end do

    outmap=map
    do i = 1, nmaps
       if (ordering==1) call convert_nest2ring(nside, outmap(0:npix-1,i))
    end do
    ! Output result map
    call write_result_map(outfile, nside, ordering, header, outmap)

  end subroutine median_filter_misspix_neighbor

  subroutine mask2misspix
    implicit none

    integer(i4b)       :: i, j 
    real(dp)           :: nullval, missval
    real(dp)           :: nullval_dp
    logical(lgt)       :: anynull
    character(len=512) :: infile, maskfile, outfile, partext
    integer(i4b)       :: ordering, nmaps, ordering2, nmaps2
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    integer(i4b) :: nside, npix, nside2

    real(dp),     allocatable, dimension(:,:)   :: map, outmap, mask
    character(len=80),         dimension(1:180) :: header, header2

    missval = -1.6375d30

    if (iargc() /= 4) then
       write(*,*) 'Usage: map_editor mask2misspix [inmap] [mask] [outmap]'
       stop
    end if

    call getarg(2,infile)
    call getarg(3,maskfile)
    call getarg(4,outfile)

    ! Read input data set 1
    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps), outmap(0:npix-1,nmaps), mask(0:npix-1,nmaps))
    call read_bintab(infile, map, npix, nmaps, nullval, anynull, header=header)
    
    i = getsize_fits(maskfile, nside=nside2, ordering=ordering2, nmaps=nmaps2)

    if (nside/=nside2) then
       write(*,*) ""
       write(*,*) "Nsides of input map and mask do not match"
       write(*,*) ""
       stop
    elseif (nmaps/=nmaps2) then
       write(*,*) ""
       write(*,*) "Nmaps of input map and mask do not match"
       write(*,*) ""
       stop
    endif

    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull, header=header2)

    if (ordering /= ordering2) then
       if (ordering2 == 2) then
          do i = 1, nmaps
             call convert_nest2ring(nside, mask(0:npix-1,i))
          end do
       else
          do i = 1, nmaps
             call convert_ring2nest(nside, mask(0:npix-1,i))
          end do
       end if
    end if

    where (mask < 0.5)
       outmap=missval
    elsewhere
       outmap=map
    end where

    ! Output result map
    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=nmaps==3)
    call write_result_map(outfile, nside, ordering, header, outmap)

  end subroutine mask2misspix

  subroutine expand_mask
    implicit none

    integer(i4b)       :: i, j, k, l, m, n1, n2, lmax
    real(dp)           :: nullval, mean1, mean2
    real(dp)           :: nullval_dp, radius, vec(3), threshold
    logical(lgt)       :: anynull, expand_pos
    character(len=512) :: infile, outfile, partext
    integer(i4b)       :: ordering, nmaps
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    integer(i4b) :: nside, npix

    real(dp),     pointer,     dimension(:,:)   :: weights
    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     allocatable, dimension(:,:)   :: map, outmap, mask, mask2, mask3
    character(len=80),         dimension(1:180) :: header
    integer(i4b), allocatable, dimension(:)     :: listpix1, listpix2

    if (iargc() < 4) then
       write(*,*) ''
       write(*,*) 'Usage: map_editor expand_mask [mask_in] [radius arcmin] [mask_out] (options)'
       write(*,*) ''
       write(*,*) 'Options:'
       write(*,*) '  -pos         Expands the positive part of the mask'
       write(*,*) '               (the negative, < 0.5, by default)'
       write(*,*) ''
       stop
    end if

    expand_pos=.false.

    call getarg(2,infile)
    call getarg(3,partext)
    read(partext,*) radius ! In arcmin
    call getarg(4,outfile)
    radius = radius/60.d0 * pi/180.d0
    if (iargc() > 4) then
       call getarg(5,partext)
       if (trim(partext) == '-pos') expand_pos=.true.
    end if

    ! Read input data set 1
    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps), outmap(0:npix-1,nmaps), mask(0:npix-1,nmaps))
    call read_bintab(infile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(0:npix-1,i))
       end do
    end if

    ! Loop through all pixels, search for holes (zeroes) or tops (ones, expand_pos=.true.)
    allocate(listpix1(0:1000000), listpix2(0:1000000))
    mask = map
    do j = 1, nmaps
       do i = 0, npix-1
          if (mod(i,1000000) == 0) write(*,*) j, i, npix
          if (expand_pos) then
             if (mask(i,j) < 0.5d0) cycle
             call pix2vec_ring(nside, i, vec)
             call query_disc(nside, vec, radius, listpix1, n1) 
             !expand ones in map/ set pixels within radius of 'tops' to one  
             map(listpix1(0:n1-1),j) = 1.d0
          else
             if (mask(i,j) > 0.5d0) cycle
             call pix2vec_ring(nside, i, vec)
             call query_disc(nside, vec, radius, listpix1, n1) 
             !expand zeroes in map/ set pixels within radius of 'holes' to zero  
             map(listpix1(0:n1-1),j) = 0.d0
          end if
       end do
    end do

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_ring2nest(nside, map(0:npix-1,i))
       end do
    end if

    ! Output result map
    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=nmaps==3)
    call write_result_map(outfile, nside, ordering, header, map)

  end subroutine expand_mask

  subroutine ud_grade_mask
    implicit none

    character(len=512)  :: filename_in, filename_out
    character(len=512)  :: partext
    integer(i4b)        :: nside_out
    logical(lgt)        :: double_precision
    real(dp)            :: threshold_value

    real(dp)     :: nullval
    logical(lgt) :: anynull
    integer(i4b) :: i, j, nside_in, ordering, nmaps, npix_in, npix_out
    real(dp), allocatable, dimension(:,:) :: map_in, map_out
    character(len=80), dimension(180)     :: header

    if (iargc() == 2) then
       call getarg(2,partext)
       if (trim(partext) == 'help') then
          write(*,*) 'Usage: map_editor ud_grade_mask [input map] [nside_out] [output map] [threshold value]'
          write(*,*) 'In the case of downgrading all downgraded pixels with value less than [threshold] is set to zero'
          stop
       end if
    elseif (iargc() < 4 .or. iargc() > 5) then
       write(*,*) 'Usage: map_editor ud_grade_mask [input map] [nside_out] [output map] [threshold value (downgrading)]'
       stop
    end if
    

    call getarg(2,filename_in)
    call getarg(3,partext)
    read(partext,*) nside_out
    call getarg(4,filename_out)
    if (iargc() == 5) then
       call getarg(5,partext)
       read(partext,*) threshold_value
    end if

    i = getsize_fits(filename_in, nside=nside_in, ordering=ordering, nmaps=nmaps)

    if (nside_in > nside_out) then
       if (iargc() /= 5) then
          write(*,*) 'Error: Threshold value must be specified in case of downgrading'
          stop
       elseif (threshold_value < 0 .or. threshold_value > 1) then
          write(*,*) 'Error: Threshold value must be in range (0,1)'
          stop
       end if
    end if
    

    npix_in  = nside2npix(nside_in)
    npix_out = nside2npix(nside_out)

    allocate(map_in(0:npix_in-1,nmaps))
    allocate(map_out(0:npix_out-1,nmaps))
    call read_bintab(filename_in, map_in, npix_in, nmaps, nullval, anynull, header=header)
    where (map_in < -1.6d30) 
       map_in = -1.6375d30
    end where

    if (ordering == 1) then
       call udgrade_ring(map_in, nside_in, map_out, nside_out)
    else
       call udgrade_nest(map_in, nside_in, map_out, nside_out)
    end if

    !thresholding out_map
    if (nside_in > nside_out) then    
       where (map_out < threshold_value)
          map_out = 0.d0
       elsewhere
          map_out = 1.d0
       end where
    end if

    call write_result_map(filename_out, nside_out, ordering, header, map_out, double_precision)

    if (count(map_in==-1.6375d30) > 0) then
       do i=1,nmaps
          write(*,*) count(map_in(:,i)==-1.6375d30), ' missing pixels out of ', npix_in
       end do
    end if
    if (count(map_out==-1.6375d30) > 0) then
       do i=1,nmaps
          write(*,*) count(map_out(:,i)==-1.6375d30), ' missing pixels out of ', npix_out
       end do
    end if
    deallocate(map_in)
    deallocate(map_out)

  end subroutine ud_grade_mask

  subroutine ud_grade_rms
    implicit none

    character(len=512)  :: filename_in, filename_out
    character(len=512)  :: partext
    integer(i4b)        :: nside_out
    logical(lgt)        :: double_precision

    real(dp)     :: nullval
    logical(lgt) :: anynull
    integer(i4b) :: i, j, nside_in, ordering, nmaps, npix_in, npix_out, rel_diff_npix
    real(dp), allocatable, dimension(:,:) :: map_in, map_out
    character(len=80), dimension(180)     :: header

    if (iargc() == 2) then
       call getarg(2,partext)
       if (trim(partext) == 'help') then
          write(*,*) 'Usage: map_editor ud_grade_rms [input map] [nside_out] [output map]'
          write(*,*) 'In the case of downgrading all downgraded pixels with value less than [threshold] is set to zero'
          stop
       end if
    elseif (iargc() /= 4 ) then
       write(*,*) 'Usage: map_editor ud_grade_rms [input map] [nside_out] [output map]'
       stop
    end if
    

    call getarg(2,filename_in)
    call getarg(3,partext)
    read(partext,*) nside_out
    call getarg(4,filename_out)

    i = getsize_fits(filename_in, nside=nside_in, ordering=ordering, nmaps=nmaps)
    
    npix_in  = nside2npix(nside_in)
    npix_out = nside2npix(nside_out)

    rel_diff_npix = npix_in/npix_out
    allocate(map_in(0:npix_in-1,nmaps))
    allocate(map_out(0:npix_out-1,nmaps))
    call read_bintab(filename_in, map_in, npix_in, nmaps, nullval, anynull, header=header)
    
    map_in=map_in(:,:)*map_in(:,:) !when ud_grading rms maps, one must use variances



    if (ordering == 1) then
       call udgrade_ring(map_in, nside_in, map_out, nside_out)
    else
       call udgrade_nest(map_in, nside_in, map_out, nside_out)
    end if

    map_out(:,:)=sqrt(map_out(:,:)/rel_diff_npix)

    call write_result_map(filename_out, nside_out, ordering, header, map_out, double_precision)

    deallocate(map_in)
    deallocate(map_out)

  end subroutine ud_grade_rms


  subroutine fill_zero_mask_udgrade
    implicit none

    character(len=512)  :: filename_in, maskname, filename_out
    character(len=512)  :: partext
    integer(i4b)        :: nside_upper, nside_temp
    logical(lgt)        :: double_precision, mainzeropix, nsidenotzero

    real(dp)     :: nullval, missval
    logical(lgt) :: anynull
    integer(i4b) :: i, j, k, nside_in, ordering, ordering2, nmaps2, nmaps, npix_in, npix_temp, npix_lower
    real(dp), allocatable, dimension(:,:)   :: map_in, map_out, mask, mask_lower, map_ud, temp
    integer(i4b), allocatable, dimension(:) :: maskzerocount, nside_lower
    character(len=80), dimension(180)       :: header

    missval=-1.6375d30

    if (iargc() == 2) then
       call getarg(2,partext)
       if (trim(partext) == 'help') then
          write(*,*) 'Usage: map_editor fill_zero_mask_udgrade [input map] [mask] [output map]'
          write(*,*) 'Subroutine fills in new values where the mask < 0.5 based on ud-grading.'
          write(*,*) 'Mask is downgraded until there is no more "zero" pixels or nside = 1.'
          write(*,*) 'The map is ud-graded down to given nside, excluding pixels where mask = zero,'
          write(*,*) 'and back up again to full resolution. Output map have the values where'
          write(*,*) 'mask < 0.5 replaced by the ud-graded map.'
          stop
       end if
    elseif (iargc() /= 4 ) then
       write(*,*) 'Usage: map_editor ud_grade_rms [input map] [mask] [output map]'
       stop
    end if
    

    call getarg(2,filename_in)
    call getarg(3,maskname)
    call getarg(4,filename_out)

    i = getsize_fits(filename_in, nside=nside_upper, ordering=ordering, nmaps=nmaps)
    
    npix_in  = nside2npix(nside_upper)
    
    allocate(map_in(0:npix_in-1,nmaps), map_ud(0:npix_in-1,nmaps))
    call read_bintab(filename_in, map_in, npix_in, nmaps, nullval, anynull, header=header)
    
    i = getsize_fits(maskname, nside=nside_in, ordering=ordering2, nmaps=nmaps2)

    if (nside_in /= nside_upper) then
       write(*,*) "Nside of mask and input map doesn't match. Terminaling."
       stop
    else if (nmaps2 < nmaps) then
       write(*,*) "Nmaps of mask is less than nmaps of input map. Using map1 of mask for all nmaps."
    end if

    allocate(temp(0:npix_in-1,nmaps2))
    allocate(mask(0:npix_in-1,nmaps))

    call read_bintab(maskname, temp, npix_in, nmaps2, nullval, anynull, header=header)

    if (ordering /= ordering2) then
       do i=1,nmaps2
          if (ordering == 1) then
             call convert_nest2ring(nside_in, temp(:,i))
          else
             call convert_ring2nest(nside_in, temp(:,i))
          end if
       end do
    end if
       

    do i = 1,nmaps
       if (nmaps2 < nmaps) then
          mask(:,i)=temp(:,1)
          j=1
       else
          mask(:,i)=temp(:,i)
          j=nmaps
       end if
    end do
    deallocate(temp)
    allocate(maskzerocount(j), nside_lower(j))
    maskzerocount(:)=0
    map_ud(:,:)=map_in(:,:)
    where (mask < 0.5)
       mask = 0.d0
       map_ud = missval
    elsewhere
       mask = 1.d0
    end where
        

    do i=1,j
       mainzeropix=.false.
       maskzerocount(i)=count(mask(:,i) < 0.5)
       nside_lower(i)=nside_upper
       nside_temp=nside_upper
       if (maskzerocount(i)==0) cycle
       nsidenotzero=.false.
       if (nside_temp >= 1) nsidenotzero=.true.
       allocate(mask_lower(0:npix_in-1,1))
       mask_lower(:,1)=mask(:,i)
       do while (maskzerocount(i) > 0 .and. nsidenotzero)
          if (nside_temp == 1) then
             nsidenotzero=.false.
             mainzeropix=.true.
             cycle
          end if
          npix_temp=12*nside_temp**2
          nside_lower(i)=nside_temp/2
          npix_lower=12*nside_lower(i)**2
          allocate(temp(0:npix_temp-1,1))
          temp(:,1)=mask_lower(:,1)
          deallocate(mask_lower)
          allocate(mask_lower(0:npix_lower-1,1))
          
          if (ordering == 1) then
             call udgrade_ring(temp(:,1), nside_temp, mask_lower(:,1), nside_lower(i))
          else
             call udgrade_nest(temp(:,1), nside_temp, mask_lower(:,1), nside_lower(i))
          end if

          where (mask_lower > 0.1)
             mask_lower = 1.d0
          elsewhere
             mask_lower = 0.d0
          end where
          
          maskzerocount(i) = count(mask_lower(:,1) < 0.5)
          nside_temp=nside_lower(i)
          deallocate(temp)
       end do

       if (mainzeropix) then
          write(*,*) "mask has a main pixel equal to zero. Terminating."
          stop
       end if
    end do


    do i = 1,nmaps
       if (nmaps2 < nmaps) then
          j=1
       else
          j=i
       end if

       if (nside_lower(j) == nside_upper) cycle !no missing pixels/mask pixels equal to zero
       npix_lower=12*nside_lower(j)**2
       allocate(temp(0:npix_lower-1,1))
       !ud_grade the input map
       if (ordering == 1) then
          call udgrade_ring(map_ud(:,i), nside_upper, temp(:,1), nside_lower(j))
       else
          call udgrade_nest(map_ud(:,i), nside_upper, temp(:,1), nside_lower(j))
       end if

       if (ordering == 1) then
          call udgrade_ring(temp(:,1), nside_lower(j), map_ud(:,i), nside_upper)
       else
          call udgrade_nest(temp(:,1), nside_lower(j), map_ud(:,i), nside_upper)
       end if
       deallocate(temp)
    end do
    
    where (mask < 0.5)
       map_out=map_ud
    elsewhere
       map_out=map_in
    end where
    
    call write_result_map(filename_out, nside_upper, ordering, header, map_out, double_precision)

    if (allocated(map_in)) deallocate(map_in)
    if (allocated(map_out)) deallocate(map_out)
    if (allocated(mask)) deallocate(mask)
    if (allocated(mask_lower)) deallocate(mask_lower)
    if (allocated(maskzerocount)) deallocate(maskzerocount)
    if (allocated(nside_lower)) deallocate(nside_lower)
    if (allocated(map_ud)) deallocate(map_ud)


  end subroutine fill_zero_mask_udgrade



  subroutine map_read_header
    implicit none

    character(len=512)  :: filename_in
    character(len=512)  :: partext
    logical(lgt)        :: double_precision

    real(dp)     :: nullval
    logical(lgt) :: anynull
    integer(i4b) :: i,j, nside_in, ordering, nmaps, npix_in, rel_diff_npix
    real(dp), allocatable, dimension(:,:) :: map_in
    character(len=80), dimension(180)     :: header

    if (iargc() == 2) then
       call getarg(2,partext)
       if (trim(partext) == 'help') then
          write(*,*) 'Usage: map_editor read_header [map]'
          stop
       end if
    else
       write(*,*) 'Usage: map_editor read_header [map]'
       stop
    end if
    

    call getarg(2,filename_in)

    i = getsize_fits(filename_in, nside=nside_in, ordering=ordering, nmaps=nmaps)
    
    npix_in  = nside2npix(nside_in)
 
    allocate(map_in(0:npix_in-1,nmaps))
    call read_bintab(filename_in, map_in, npix_in, nmaps, nullval, anynull, header=header)

    write(*,*) ""
    write(*,*) "Header length =", len(header)
    write(*,*) ""
    write(*,*) "---------------------------------- HEADER ---------------------------------"
    write(*,*) ""

    do i =1,len(header)
       if (header(i)==" ") then
          j=1
       else
          write(*,*) header(i)
       end if

    end do
    deallocate(map_in)

  end subroutine map_read_header


  subroutine map_size_info
    implicit none
    
    character(len=512)          :: string, infile
    integer(i4b)                :: npix, ordering, nside, nmaps
    real(dp)                    :: nullval
    logical(lgt)                :: anynull

    ! Get parameters
    if (iargc() /= 2) then
       write(*,*) 'Usage: map_editor size_info [input map]'
       return
    else 
       call getarg(2, string)
    end if
    read(string,*) infile
    
    ! Read map
    npix = getsize_fits(infile, nmaps=nmaps, ordering=ordering, nside=nside)
    write(*,*) 'Nside = ', nside
    write(*,*) 'Nmaps = ', nmaps
    write(*,*) 'Ordering = ', ordering
    
  end subroutine map_size_info

  subroutine create_source_mask_from_file
    implicit none

    integer(i4b)       :: i, j, k, l, m, n1, n2, stat
    real(dp)           :: nullval, mean1, mean2
    real(dp)           :: nullval_dp, radius, vec(3), gal_cor(2)
    logical(lgt)       :: anynull
    character(len=512) :: infile, sourcefile, outfile, partext
    character(len=512) :: line
    integer(i4b)       :: ordering, nmaps
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    integer(i4b) :: nside, npix

    real(dp),     allocatable, dimension(:,:)   :: outmap, mask
    character(len=80),         dimension(1:180) :: header
    integer(i4b), allocatable, dimension(:)     :: listpix1

    if (iargc() == 2) then
       call getarg(2,infile)
       if (infile .eq. 'help') then
          write(*,*) 'Usage: map_editor create_source_mask [fullsky mask] [source file in (txt)] [map_out]'
          write(*,*) ''
          write(*,*) 'Outputs a mask file of nside and nmaps given by the fullsky mask.' 
          write(*,*) 'Sources are masked out with corresponding radius given in the source file.'
          write(*,*) ''
          write(*,*) 'Source file setup:'
          write(*,*) '#  Galactic coordinates | radius of source | optional comment '
          write(*,*) '#  Lat (deg)  Lon (deg)   radius (arcmin)'
          write(*,*) '   -6.341      185.32       80.0              source1 '
          write(*,*) '  -17.723      121.64       60.0              source2 '
          write(*,*) '   57.02        83.1       120.0              source3 '
          write(*,*) ''
          write(*,*) ''
          write(*,*) '!!!IMPORTANT!!!'
          write(*,*) "Commented lines MUST start with '#' as the FIRST character"
          stop
       end if
    end if
    
    if (iargc() /= 4) then
       write(*,*) 'Usage: map_editor create_source_mask [fullsky mask] [source file in (txt)] [map_out]'
       write(*,*) 'For extended help, use: map_editor create_source_mask help'
       stop
    end if

    call getarg(2,infile)
    call getarg(3,sourcefile)
    call getarg(4,outfile)

    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(mask(0:npix-1,nmaps), outmap(0:npix-1,nmaps))
    call read_bintab(infile, mask, npix, nmaps, nullval, anynull, header=header)    
    
    mask(:,:) = 1.d0 !removing any zero values
    
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, mask(0:npix-1,i))
       end do
    end if

    allocate(listpix1(0:npix-1))
    
    ! get source vectors
    ! read from file
    open(1,file=sourcefile)
    write(*,*) 'Source file read successfull'
    j=1
    do 
       read(1, '(a)', IOSTAT=stat) line !There is an error here!!! check it out!!!
       if (stat > 0) then
          write(*,*) "Something went wrong while trying to read the file"
          return
       end if
       if (stat < 0) exit !end of file reached

       !else we extract the source data
       line=trim(line)
       if (line(1:1)=='#') cycle ! a commented line
       write(*,*) 'Reading source #', j
       j=j+1
       read(line,*) gal_cor(1),gal_cor(2),radius,partext
       !galactic coord -> spherical coord
       gal_cor(1) = 90.d0 - gal_cor(1)
       gal_cor(1) = gal_cor(1)*pi/180.d0
       gal_cor(2) = gal_cor(2)*pi/180.d0
       call ang2vec(gal_cor(1),gal_cor(2),vec(1:3))

       radius = radius/60.d0 * pi/180.d0
       call query_disc(nside, vec, radius, listpix1, n1)

       do i =1,nmaps
          mask(listpix1(0:n1-1),i) = 0.d0
       end do
    end do


    if (ordering == 2) then
       do i = 1, nmaps
          call convert_ring2nest(nside, mask(0:npix-1,i))
       end do
    end if
    
    outmap(:,:)=mask(:,:)

    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=nmaps==3)
    call write_result_map(outfile, nside, ordering, header, outmap)

  end subroutine create_source_mask_from_file
  
  subroutine copy_T_to_QU
        implicit none

    character(len=512)  :: filename_in, filename_out
    character(len=512)  :: partext
    logical(lgt)        :: double_precision
    real(dp)            :: threshold_value

    real(dp)     :: nullval
    logical(lgt) :: anynull
    integer(i4b) :: i, j, nside_in, ordering, nmaps, npix_in, npix_out
    real(dp), allocatable, dimension(:,:) :: map_in, map_out
    character(len=80), dimension(180)     :: header

    if (iargc() > 1) then
       call getarg(2,partext)
       if (trim(partext) == 'help') then
       write(*,*) '' 
       write(*,*) '   The copy_T_to_QU operation reads an input map and outputs a map of'
       write(*,*) '   nside = 3, where map1 is a copy of map1 (T polarization) of the input'
       write(*,*) '   map, and map 2 and 3 is copies of map1'
       write(*,*) ''
       write(*,*) '   Usage: map_editor copy_T_to_QU [input map] [output map]'
       write(*,*) ''
       stop
       end if
    end if
    if (iargc() /= 3) then
       write(*,*) 'Usage: map_editor copy_T_to_QU [input map] [output map]'
       stop
    end if
    

    call getarg(2,filename_in)
    call getarg(3,filename_out)

    i = getsize_fits(filename_in, nside=nside_in, ordering=ordering, nmaps=nmaps)

    npix_in  = nside2npix(nside_in)

    allocate(map_in(0:npix_in-1,nmaps))
    allocate(map_out(0:npix_in-1,3))
    call read_bintab(filename_in, map_in, npix_in, nmaps, nullval, anynull, header=header)
    where (map_in < -1.6d30) 
       map_in = -1.6375d30
    end where
    nmaps=3
    map_out(:,1)=map_in(:,1)
    map_out(:,2)=map_in(:,1)
    map_out(:,3)=map_in(:,1)
    call write_minimal_header(header, 'MAP', nside=nside_in, order=ordering, polar=nmaps==3)
    call write_result_map(filename_out, nside_in, ordering, header, map_out, double_precision)
    deallocate(map_in)
    deallocate(map_out)

  end subroutine copy_T_to_QU

  subroutine copy_pol_to_T
        implicit none

    character(len=512)  :: filename_in, filename_out
    character(len=512)  :: partext
    logical(lgt)        :: double_precision
    real(dp)            :: threshold_value

    real(dp)     :: nullval
    logical(lgt) :: anynull
    integer(i4b) :: i, map_num, nside_in, ordering, nmaps, npix_in, npix_out
    real(dp), allocatable, dimension(:,:) :: map_in, map_out
    character(len=80), dimension(180)     :: header

    if (iargc() > 1) then
       call getarg(2,partext)
       if (trim(partext) == 'help') then
          write(*,*) ''
          write(*,*) '   The copy_pol_to_T operation reads an input map and outputs a map of'
          write(*,*) '   nside = 1, where map1 is a copy one polarization of the input map.'
          write(*,*) ''  
          write(*,*) '   Usage: map_editor copy_pol_to_T [input map] [pol] [output map]'
          write(*,*) ''
          stop
       end if
    end if
    if (iargc() /= 4) then
       write(*,*) 'Usage: map_editor copy_pol_to_T [input map] [pol/map#] [output map]'
       stop
    end if
    

    call getarg(2,filename_in)
    call getarg(3,partext)
    read(partext,*) map_num
    call getarg(4,filename_out)

    i = getsize_fits(filename_in, nside=nside_in, ordering=ordering, nmaps=nmaps)

    if (map_num > nmaps) then
       write(*,*) 'Error: The specified map# (signal) does not exist'
       stop
    endif
    npix_in  = nside2npix(nside_in)

    allocate(map_in(0:npix_in-1,nmaps))
    allocate(map_out(0:npix_in-1,1))
    call read_bintab(filename_in, map_in, npix_in, nmaps, nullval, anynull, header=header)
    where (map_in < -1.6d30) 
       map_in = -1.6375d30
    end where
    nmaps=1
    map_out(:,1)=map_in(:,map_num)

    call write_minimal_header(header, 'MAP', nside=nside_in, order=ordering, polar=nmaps==3)
    call write_result_map(filename_out, nside_in, ordering, header, map_out, double_precision)
    deallocate(map_in)
    deallocate(map_out)

  end subroutine copy_pol_to_T

  subroutine firas_calib(bandpass_in,weight_file,output_map) !,units)
    implicit none
    
    integer(i4b)       :: count,start,length,i,j,k,bands
    integer(i4b)       :: nside, ordering, nmaps, npix
    integer(i4b)       :: nside_mask, ordering_mask, nmaps_mask,npix_mask
    real(dp)           :: norm, nullval
    character(len=3)   :: suffix
    character(len=512) :: bandpass_in,string,firas_file,firas_maps,weight_file,output_map !,units
    character(len=512) :: firmap1, firmap2,firas_path,firas_mask !,firas_conv
    logical(lgt)       :: exist
    logical(lgt)       :: double_precision
    logical(lgt)       :: anynull
    real(dp)           :: missval = -1.6375d30
    
    real(dp), allocatable, dimension(:)   :: firas,freq,weights !,convert
    real(dp), allocatable, dimension(:)   :: firas_weights,test1,test2
    real(dp), allocatable, dimension(:,:) :: map_1,map_2,new_map,mask
    character(len=80), dimension(180)     :: header
    character(len=512), allocatable       :: str(:) 

    firas_file='/mn/stornext/d14/Planck1/daniher/data/firas/firas_frequencies.dat'
    firas_maps='/mn/stornext/d14/Planck1/daniher/data/firas/firas_maps.dat'
    firas_path='/mn/stornext/u3/hke/xsan/firas/'
    firas_mask='/mn/stornext/u3/hke/xsan/firas/FIRAS_mask.fits'
    ! firas_conv='/mn/stornext/d14/Planck1/daniher/data/firas/firas_conversions.dat'

    inquire(file=firas_file,exist=exist)
    if (.not. exist) then
       write(*,*) ''
       write(*,*) '   No FIRAS frequency file found!'
       write(*,*) ''
       write(*,*) '   The FIRAS frequency file contains a list of the FIRAS'
       write(*,*) '   frequencies. A copy can be found at:'
       write(*,*) '   /mn/stornext/d14/Planck1/daniher/data/firas/firas_frequencies.dat'
       write(*,*) ''
       stop
    endif

    inquire(file=firas_maps,exist=exist)
    if (.not. exist) then
       write(*,*) ''
       write(*,*) '   No FIRAS maps file found!'
       write(*,*) ''
       write(*,*) '   The FIRAS maps file contains a list of the FIRAS map names'
       write(*,*) '   A copy can be found at:'
       write(*,*) '   /mn/stornext/d14/Planck1/daniher/data/firas/firas_maps.dat'
       write(*,*) ''
       stop
    endif

    ! Reads the number of FIRAS frequency entries (should be 210)
    length = 0
    open(29,file=firas_file)
    do while (.true.)
       read(29,*,END=1) string
       if (string(1:1)=='#') cycle
       length = length + 1
    end do

1   close(29)

    ! Reads the number of map frequency bands
    bands = 0
    open(31,file=bandpass_in)
    do while (.true.)
       read(31,*,END=2) string
       if (string(1:1)=='#') cycle
       bands = bands + 1
    end do

2   close(31)

    ! Allocating arrays
    allocate(firas(length),firas_weights(length))
    ! allocate(convert(length))
    allocate(freq(bands),weights(bands))
    allocate(test1(length)) !,test2(length))
    allocate(str(length))
    test1 = 0.d0
    test2 = 0.d0
        
    ! Creating the FIRAS array
    k = 0
    open(29,file=firas_file)
    do while (.true.)
       read(29,fmt='(a)',END=3) string
       if (string(1:1)=='#') cycle
       k = k + 1
       read(string,*) firas(k)
    end do

3   close(29)

    ! Calculating the f2t conversion factors for each frequency (MJy/sr -> uK)
    
    ! open(45,file=firas_conv)
    ! do i=1,length
       ! nu = firas(i)*1.d9
       ! write(*,*) nu
       ! x = h*nu/(k_B*T_CMB)
       ! convert(i) = 1.d0/(1.d14*((2.d0*h*nu**3/(c**2*(exp(x)-1.d0))) * &
       !      & (exp(x) / (exp(x)-1.d0)) * h*nu/(k_B*T_CMB**2)))
       ! write(45,'(2(E17.8))') firas(i), convert(i)
    ! end do
    ! close(45)

    ! Creating array of FIRAS map names
    open(30,file=firas_maps)
    do i=1,length
       read(30,fmt='(a)') str(i)
       str(i) = trim(firas_path) // trim(str(i))
    end do

    ! Creating the bandpsas arrays
    k = 0
    open(31,file=bandpass_in)
    do while (.true.)
       read(31,fmt='(a)',END=4) string
       if (string(1:1)=='#') cycle
       k = k + 1
       read(string,*) freq(k), weights(k)
    end do

4   close(31)

    ! Normalizing the bandpass area
    norm = sum(weights)
    weights = weights/norm
    
    start = 1 
    do while(firas(start) .LT. freq(1))
       start = start + 1
    end do
    
    count = 1 
    do while (firas(start+count) .GT. freq(1) .and. firas(start+count) .LT. freq(bands))
       count = count + 1
    end do

    do i=1,count
       k = 0
       do j=1,bands-1
          if ((firas(i+start) + 6.8d0) .GT. freq(j) .and. (firas(i+start) - 6.8d0) .LT. freq(j)) then
             test1(i+start) = test1(i+start) + weights(j)
             k = k + 1
          endif
       end do
       if (k .GT. 1) then
          test1(i+start) = test1(i+start)/k
       endif
    end do

    firas_weights = test1

    norm = sum(firas_weights)
    firas_weights = firas_weights/norm

    do i=1,count
       if (firas_weights(i+start) .LT. 1.0d-4) then
          firas_weights(i+start) = 0.d0
       end if
    end do

    norm = sum(firas_weights)
    firas_weights = firas_weights/norm

    open(32,file=weight_file)
       do i=1,length
          write(32,'(I4,F8.4)') int(firas(i)), firas_weights(i)
       end do
    close(32)

    deallocate(test1)
    deallocate(freq)
    deallocate(weights)

    ! Adding together weighted FIRAS maps

    firmap1 = str(start)
    
    i = getsize_fits(firmap1, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    npix_mask = getsize_fits(firas_mask, nside=nside_mask, ordering=ordering_mask, nmaps=nmaps_mask)
    if (nmaps_mask /= nmaps) then
       write(*,*) 'Error: Different nmaps of the FIRAS_mask and FIRAS_map. Setting nmaps=1.'
       nmaps=1
    end if

    allocate(map_1(0:npix-1,nmaps))
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(firmap1, map_1, npix, nmaps, nullval, anynull, header=header)
    call read_bintab(firas_mask, mask, npix, nmaps, nullval, anynull)

    ! if (unit == 2) then
    !    map_1 = map_1*convert(start)
    ! end if

    do i = 1, nmaps
       if (ordering_mask /= ordering) then
          if (ordering_mask == 1) call convert_ring2nest(nside, mask(:,i))
          if (ordering_mask == 2) call convert_nest2ring(nside, mask(:,i))
       end if
    end do

    allocate(new_map(0:npix-1,nmaps))

    ! Initializing the mask on the first firas map
    do k = 1, nmaps
       do j = 0, npix-1
          if (mask(j,k) < 0.5d0) then
             if (nmaps > 1) then
                map_1(j,k) = -1.6375e30
             else
                map_1(j,:) = -1.6375e30
             end if
          end if
       end do
    end do

    ! Begin adding the weighted maps together
    do k = 1, nmaps
       do j = 0, npix-1
          if (abs((map_1(j,k)-missval)/missval) > 1d-5) then
             new_map(j,k) = firas_weights(start)*map_1(j,k)
          end if
       end do
    end do

    do i=1,count-1
       firmap2 = str(start+i)
       call read_bintab(firmap2, map_1, npix, nmaps, nullval, anynull, header=header)
       ! if (unit == 2) then
       !    map_1 = map_1*convert(start+i)
       ! end if
       do k = 1, nmaps
          do j = 0, npix-1
             if (mask(j,k) < 0.5d0) then
                if (nmaps > 1) then
                   map_1(j,k) = -1.6375e30
                else
                   map_1(j,:) = -1.6375e30
                end if
             end if
          end do
       end do
       do k = 1, nmaps
          do j = 0, npix-1
             new_map(j,k) = new_map(j,k) + firas_weights(start+i)*map_1(j,k)
          end do
       end do
    end do

    call write_result_map(output_map, nside, ordering, header, new_map, double_precision)

    deallocate(new_map)
    deallocate(firas)
    deallocate(firas_weights)
    deallocate(str)
    ! deallocate(convert)

  end subroutine firas_calib


  subroutine calculate_rms_half_ring_diff
    implicit none

    integer(i4b) :: npix, nside, nside2, ordering, ordering2, nmaps, nmaps2
    integer(i4b) :: i, j, n, nmask, npol
    real(dp)     :: nullval, mu_half, mu_exp, sigma, dtheta, var_half, var_exp
    real(dp)     :: fmissval = 0.
    logical(lgt) :: anynull, double_precision
    character(len=512) :: partext, outfile
    real(dp), allocatable, dimension(:,:) :: half1, half2, rms1, rms2, map, mask, outmap, inmap
    real(dp), allocatable, dimension(:)   :: Ascale
    character(len=80), dimension(180)     :: header
    if (iargc() < 8 .or. iargc() > 9) then
       write(*,*) "Usage: map_editor rms_half_ring_diff [half_ring1] [half_ring2] [rms1] [rms2]"
       write(*,*) "      [mask] [fullfreq rms map] [outmap] [npol (1/3), optional]"
       write(*,*) ""
       stop
    end if

    
    ! Read input maps
    
    call getarg(2, partext)
    i = getsize_fits(trim(partext), nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    if (iargc() == 9) then
       call getarg(9,partext)
       read(partext,*) npol
    else
       npol = nmaps
    end if
    if (nmaps < npol) then
       write(*,*) "nmaps of half_ring1 too small for the number of polarizations (npol)"
       stop
    end if
    call getarg(2, partext)
    allocate(half1(0:npix-1,nmaps), half2(0:npix-1,nmaps), rms1(0:npix-1,nmaps), rms2(0:npix-1,nmaps))
    call read_bintab(trim(partext), half1, npix, nmaps, nullval, anynull)
    call getarg(3, partext)
    i = getsize_fits(trim(partext), nside=nside2, ordering=ordering2, nmaps=nmaps2)
    if (nmaps /= nmaps2) then
       write(*,*) "nmaps between half_ring1 and half_ring2 doesn't match"
       stop
    else if (nside /= nside2) then
       write(*,*) "nside between half_ring1 and half_ring2 doesn't match"
       stop
    end if
    call read_bintab(trim(partext), half2, npix, nmaps, nullval, anynull)
    if (ordering /= ordering2) then
       do j = 1,nmaps
          if (ordering == 1) then
             call convert_nest2ring(nside, half2(:,j))
          else
             call convert_ring2nest(nside, half2(:,j))
          end if
       end do
    end if

    call getarg(4, partext)
    i = getsize_fits(trim(partext), nside=nside2, ordering=ordering2, nmaps=nmaps2)
    if (nmaps /= nmaps2) then
       write(*,*) "nmaps between half_ring1 and rms1 doesn't match"
       stop
    else if (nside /= nside2) then
       write(*,*) "nside between half_ring1 and rms1 doesn't match"
       stop
    end if
    call read_bintab(trim(partext), rms1, npix, nmaps, nullval, anynull)
    if (ordering /= ordering2) then
       do j = 1,nmaps
          if (ordering == 1) then
             call convert_nest2ring(nside, rms1(:,j))
          else
             call convert_ring2nest(nside, rms1(:,j))
          end if
       end do
    end if

    call getarg(5, partext)
    i = getsize_fits(trim(partext), nside=nside2, ordering=ordering2, nmaps=nmaps2)
    if (nmaps /= nmaps2) then
       write(*,*) "nmaps between half_ring2 and rms2 doesn't match"
       stop
    else if (nside /= nside2) then
       write(*,*) "nside between half_ring2 and rms2 doesn't match"
       stop
    end if
    call read_bintab(trim(partext), rms2, npix, nmaps, nullval, anynull)
    if (ordering /= ordering2) then
       do j = 1,nmaps
          if (ordering == 1) then
             call convert_nest2ring(nside, rms2(:,j))
          else
             call convert_ring2nest(nside, rms2(:,j))
          end if
       end do
    end if

    call getarg(6, partext)
    i = getsize_fits(trim(partext), nside=nside2, ordering=ordering2, nmaps=nmaps2)
    if (nside /= nside2) then
       write(*,*) "nside between half_ring1 and mask doesn't match"
       stop
    end if
    allocate(inmap(0:npix-1,nmaps2))
    call read_bintab(trim(partext), inmap, npix, nmaps2, nullval, anynull)
    allocate(mask(0:npix-1,nmaps))
    if (nmaps2 < nmaps) then
       write(*,*) "nmaps of mask is smaller than that of half_ring1." 
       write(*,*) "Copying pol=1 in mask to all npol"
       do j = 1,nmaps
          mask(:,j)=inmap(:,1)
       end do
    else
       do j = 1,nmaps
          mask(:,j)=inmap(:,j)
       end do
    end if
    deallocate(inmap)
    if (ordering /= ordering2) then
       do j = 1,nmaps
          if (ordering == 1) then
             call convert_nest2ring(nside, mask(:,j))
          else
             call convert_ring2nest(nside, mask(:,j))
          end if
       end do
    end if

    call getarg(7, partext)
    i = getsize_fits(trim(partext), nside=nside2, ordering=ordering2, nmaps=nmaps2)
    if (nside /= nside2) then
       write(*,*) "nside between half_ring1 and full freq rms map doesn't match"
       stop
    end if
    allocate(map(0:npix-1,nmaps2),outmap(0:npix-1,nmaps2))
    call read_bintab(trim(partext), map, npix, nmaps2, nullval, anynull)    
    if (npol==1) then
       write(*,*) "Calculating new rms using only scaling factor from pol 1"
    else if (nmaps2 < nmaps) then
       write(*,*) "nmaps of full freq rms map is smaller than that of half_ring1"
       write(*,*) "Calculating new rms only for the nmaps of the full freq rms map"
    else if (nmaps2 > nmaps) then
       write(*,*) "nmaps of full freq rms map is larger than that of half_ring1"
       write(*,*) "Calculating new rms for the extra nmaps of the full freq rms map"
       write(*,*) "using the scaling factor from pol 1."
    end if
    if (ordering /= ordering2) then
       do j = 1,nmaps2
          if (ordering == 1) then
             call convert_nest2ring(nside, map(:,j))
          else
             call convert_ring2nest(nside, map(:,j))
          end if
       end do
    end if

    call getarg(8,outfile)

    allocate(Ascale(npol))
    do j=1,min(npol,nmaps2)
       rms1(:,j)=rms1(:,j)*rms1(:,j) + rms2(:,j)*rms2(:,j) ! setting it to sum of variances
       half1(:,j)=half1(:,j)-half2(:,j) !half ring diff
       rms1(:,j)=sqrt(rms1(:,j))
       half1(:,j)=half1(:,j)/rms1(:,j)
       ! we now expect half1(:,j) to have STD (or variance) = 1.

       ! calculate var of half1(:,j)
       nmask=0
       mu_half=0
       ! need the mean of half1 and number of pixels in the mask with value > 1
       do i = 0,npix-1
          if (mask(i,j) < 0.5) cycle
          nmask=nmask+1
          mu_half=mu_half+half1(i,j)
       end do
       mu_half=mu_half/real(nmask) !mean value
       var_half=0
       do i = 0,npix-1
          if (mask(i,j) < 0.5) cycle
          var_half = var_half + (half1(i,j)-mu_half)**2
       end do
       var_half=var_half/(nmask-1) ! getting the variance of half_diff
       
       Ascale(j)= sqrt(var_half) ! as A*rms --> A**2 * var --> 
                              ! var( var_half/(A**2 * var_exp) ) = 1
    end do

    do j =1,nmaps2
       if (npol==1) then
          write(*,*) "pol =",j,",   A_scale =",Ascale(1)
          outmap(:,j)=map(:,j)*Ascale(1)
       else if (j > npol) then
          write(*,*) "pol =",j,",   A_scale =",Ascale(1)
          outmap(:,j)=map(:,j)*Ascale(1)
       else
          write(*,*) "pol =",j,",   A_scale =",Ascale(j)
          outmap(:,j)=map(:,j)*Ascale(j)
       end if
    end do

    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=nmaps2==3)
    call write_result_map(outfile, nside, ordering, header, outmap, double_precision)
  end subroutine calculate_rms_half_ring_diff


  subroutine fit_CO_lineratio
    implicit none

    logical(lgt) :: anynull, converged,double_precision, debug
    real(dp)     :: nullval, T0, beta0, x, A, B, temp_LR
    real(dp)     :: chisq_curr, chisq_in, chisq_out, chisq_new, chisq_diff
    integer(i4b) :: numband, nside, npix, ordering, nmaps, i, j, k, ierr, n, comp, m
    integer(i4b) :: order_mask, temp_int
    integer(i4b) :: n_co, n_sky, ordering2, nmaps2, nside2, nlist, counter
    real(dp),     allocatable, dimension(:,:) :: co_maps, sky_maps, rms_maps, res_maps
    real(dp),     allocatable, dimension(:,:) :: co_maps2, diff_maps, var_maps, inv_N    
    real(dp),     allocatable, dimension(:,:) :: mask, temp_map, LR, diff_maps_full
    real(dp),     allocatable, dimension(:)   :: gain, offset, avg_co, avg_sky, data
    integer(i4b), allocatable, dimension(:)   :: listpix
    character(len=512) :: infile1, filename, outprefix, temp, maskfile, gainmaskfile
    character(len=10)  :: alpha_text
    character(len=2)   :: i_text
    character(len=80), dimension(180)         :: header
    !For mono/dipole correction
    real(dp),    allocatable, dimension(:,:)       :: md_mask, harmonics, harmonics2
    real(dp),                 dimension(0:3)       :: multipoles, b_md
    real(dp),                 dimension(3)         :: vector
    real(dp),                 dimension(0:3,0:3)   :: A_md


    if (iargc() < 8) then
       write(*,*) 'Usage: map_editor fit_CO_lineratio [n_CO_maps] [n_sky_maps] [mask] [MD_mask]'

       write(*,*) '  [CO_map1] [CO_map2] ... [sky_map1] [sky_map2] ... [rms_map1] [rms_map2] ...'
       stop
    end if

    double_precision=.true.
    write(*,*) 'Reading files'
    call getarg(2,temp)
    read(temp,*) n_co
    call getarg(3,temp)
    read(temp,*) n_sky

    debug=.false.
    if (iargc()==5+n_co+2*n_sky+1) then
       call getarg(5+n_co+2*n_sky+1,temp) 
       if (trim(temp)=='debug') debug=.true.
    end if

    call getarg(4,maskfile)
    if (debug) write(*,*) trim(maskfile)
    i = getsize_fits(trim(maskfile), nside=nside, ordering=ordering, nmaps=nmaps)
    npix  = nside2npix(nside)
    if (nmaps /= 1) then
       write(*,*) '  all maps need nmaps = 1, mask has nmaps =', nmaps
       stop
    end if
    allocate(co_maps(0:npix-1,n_co), sky_maps(0:npix-1,n_sky), rms_maps(0:npix-1,n_sky), &
         & mask(0:npix-1,1), md_mask(0:npix,1), temp_map(0:npix-1,1))
    call read_bintab(trim(maskfile), mask, npix, nmaps, nullval, anynull)
    if (ordering == 2) then
       call convert_nest2ring(nside, mask(:,1))
       ordering=1
    end if

    call getarg(5,infile1)
    if (debug) write(*,*) trim(infile1)
    i = getsize_fits(trim(infile1), nside=nside2, ordering=ordering2, nmaps=nmaps2)
    if (nmaps2 /= 1) then
       write(*,*) '  all maps need nmaps = 1, md_mask',infile1,' has nmaps =', nmaps2
       stop
    elseif (nside2 /= nside) then
       write(*,*) '  all maps need same nside, md_mask',infile1,' has nside =', nmaps2,&
            & 'mask has nside =',nside
       stop
    end if

    call read_bintab(trim(infile1), temp_map, npix, nmaps, nullval, anynull)
    md_mask(:,1)=temp_map(:,1)
    if (ordering2 == 2) then
       call convert_nest2ring(nside, md_mask(:,1))
       ordering2=1
    end if

    do j = 1,n_co
       call getarg(j+5,infile1)
       if (debug) write(*,*) trim(infile1)
       i = getsize_fits(trim(infile1), nside=nside2, ordering=ordering2, nmaps=nmaps2)
       if (nmaps2 /= 1) then
          write(*,*) '  all maps need nmaps = 1, co_map',infile1,' has nmaps =', nmaps2
          stop
       elseif (nside2 /= nside) then
          write(*,*) '  all maps need same nside, co_map',infile1,' has nside =', nmaps2,&
               & 'mask has nside =',nside
          stop
       end if

       call read_bintab(trim(infile1), temp_map, npix, nmaps, nullval, anynull)
       co_maps(:,j)=temp_map(:,1)
       if (ordering2 == 2) then
          call convert_nest2ring(nside, co_maps(:,j))
          ordering2=1
       end if
    end do

    k=5+n_co
    do j = 1,n_sky
       call getarg(j+k,infile1)
       if (debug) write(*,*) trim(infile1)
       i = getsize_fits(trim(infile1), nside=nside2, ordering=ordering2, nmaps=nmaps2)
       if (nmaps2 /= 1) then
          write(*,*) '  all maps need nmaps = 1, sky_map',infile1,' has nmaps =', nmaps2
          stop
       elseif (nside2 /= nside) then
          write(*,*) '  all maps need same nside, sky_map',infile1,' has nside =', nmaps2,&
               & 'mask has nside =',nside
          stop
       end if

       call read_bintab(trim(infile1), temp_map, npix, nmaps, nullval, anynull)
       sky_maps(:,j)=temp_map(:,1)
       if (ordering2 == 2) then
          call convert_nest2ring(nside, sky_maps(:,j))
          ordering2=1
       end if
    end do

    k=k+n_sky
    do j = 1,n_sky
       call getarg(j+k,infile1)
       if (debug) write(*,*) trim(infile1)
       i = getsize_fits(trim(infile1), nside=nside2, ordering=ordering2, nmaps=nmaps2)
       if (nmaps2 /= 1) then
          write(*,*) '  all maps need nmaps = 1, rms_map',infile1,' has nmaps =', nmaps2
          stop
       elseif (nside2 /= nside) then
          write(*,*) '  all maps need same nside, rms_map',infile1,' has nside =', nmaps2,&
               & 'mask has nside =',nside
          stop
       end if

       call read_bintab(trim(infile1), temp_map, npix, nmaps, nullval, anynull)
       rms_maps(:,j)=temp_map(:,1)
       if (ordering2 == 2) then
          call convert_nest2ring(nside, rms_maps(:,j))
          ordering2=1
       end if
    end do


    write(*,*) 'setting up diffmaps, variance maps, and mask'

    allocate(listpix(npix))
    j=0
    do i=0,npix-1
       if (mask(i,1) > 0.5) then
          j=j+1
          listpix(j)=i
       end if
    end do
    nlist=j !number of pixels inside mask
    !NOTE: we only care about these pixels from now on

    !creating diff_maps and variance_maps
    allocate(diff_maps_full(0:npix-1,n_sky-1))
    if (debug) write(*,*) 'writing diff maps and skymaps'
    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=.false.)
    do j = 1,n_sky-1
       diff_maps_full(:,j) = sky_maps(:,1)-sky_maps(:,j+1)
       if (debug) then
          temp_map(:,1)=diff_maps_full(:,j)

          temp_int = j
          do i = 1, len(i_text)
             n = temp_int / 10**(len(i_text)-i)
             write(i_text(i:i),'(I1)') n
             temp_int = temp_int - n * 10**(len(i_text)-i)
          end do

          filename='diffmap_' // trim(i_text) // '_uK.fits'
          call write_result_map(trim(filename), nside, ordering, header,temp_map(:,:1),double_precision)
       end if
    end do
    if (debug) then
       do j = 1,n_sky-1
          temp_map(:,1)=sky_maps(:,j)

          temp_int = j
          do i = 1, len(i_text)
             n = temp_int / 10**(len(i_text)-i)
             write(i_text(i:i),'(I1)') n
             temp_int = temp_int - n * 10**(len(i_text)-i)
          end do

          filename='skymap_' // trim(i_text) // '_uK.fits'
          call write_result_map(trim(filename), nside, ordering, header,temp_map(:,:1),double_precision)
       end do
    end if
    if (debug) then
       do j = 1,n_co
          temp_map(:,1)=co_maps(:,j)

          temp_int = j
          do i = 1, len(i_text)
             n = temp_int / 10**(len(i_text)-i)
             write(i_text(i:i),'(I1)') n
             temp_int = temp_int - n * 10**(len(i_text)-i)
          end do

          filename='comap_' // trim(i_text) // '_uK.fits'
          call write_result_map(trim(filename), nside, ordering, header,temp_map(:,:1),double_precision)
       end do
    end if


    !Now we need to adjust the monopoles (and dipoles) on the diffmaps so they dissappear
    allocate(harmonics(0:npix-1,0:3))
    allocate(harmonics2(0:npix-1,0:3))
    do i = 0, npix-1
       if (ordering == 1) then
          call pix2vec_ring(nside, i, vector)
       else
          call pix2vec_nest(nside, i, vector)          
       end if
          
       if (md_mask(i,1) > 0.5d0) then
          harmonics(i,0) = 1.d0
          harmonics(i,1) = vector(1)
          harmonics(i,2) = vector(2)
          harmonics(i,3) = vector(3)
       else
          harmonics(i,:) = 0.d0
       end if

       harmonics2(i,0) = 1.d0
       harmonics2(i,1) = vector(1)
       harmonics2(i,2) = vector(2)
       harmonics2(i,3) = vector(3)

    end do

    !Loop over all diffmaps
    do n = 1,n_sky-1
       write(*,*) '  Mono-/dipole correcting diffmap nr.', n
       temp_map(:,1) = diff_maps_full(:,n)
       A_md = 0.d0
       b_md = 0.d0
       do j = 0, 3
          do k = 0, 3
             A_md(j,k) = sum(harmonics(:,j) * harmonics(:,k))
          end do
          b_md(j) = sum(temp_map(:,1) * harmonics(:,j))
       end do


       call solve_system_real(A_md, multipoles, b_md)

       do i = 0, npix-1
          !       if (mask(i,1) == 1) then
          do j = 0, 3
             diff_maps_full(i,n) = temp_map(i,1) - multipoles(j) * harmonics2(i,j)
          end do
          !       else
          !          map(i,:) = -1.6375e30
          !       end if
       end do
       write(*,*) '  Coefficients = ', real(multipoles,sp)
    end do

    allocate(avg_co(n_co-1), avg_sky(n_sky-1),co_maps2(nlist,n_co), &
         & diff_maps(nlist,n_sky-1),var_maps(nlist,n_sky-1), inv_N(nlist,n_sky-1))
    do j = 1,n_sky-1
       do i = 1,nlist
          diff_maps(i,j) = diff_maps_full(listpix(i),j) !mono-/dipole corrected
          var_maps(i,j) = rms_maps(listpix(i),1)**2 + rms_maps(listpix(i),j+1)**2
          inv_N(i,:j) = 1.d0/var_maps(i,j)
       end do
       avg_sky(j) = sum(diff_maps(:,j))/nlist
    end do


    do j = 1,n_co
       do i = 1,nlist
          co_maps2(i,j) = co_maps(listpix(i),j)
       end do
       avg_co(j) = sum(co_maps2(:,j))/nlist
    end do
    if (debug) then
       do j = 1,n_co
          temp_map(:,1)=0.d0
          temp_map(listpix(:nlist),1)=co_maps2(:,j)

          temp_int = j
          do i = 1, len(i_text)
             n = temp_int / 10**(len(i_text)-i)
             write(i_text(i:i),'(I1)') n
             temp_int = temp_int - n * 10**(len(i_text)-i)
          end do

          filename='comap_cut_' // trim(i_text) // '_uK.fits'
          call write_result_map(trim(filename), nside, ordering, header,temp_map(:,:1),double_precision)
       end do
    end if

    !Now we have our reduced maps. Now we start the setup for finding the line ratios

    !Model 1:  diff = LR_12 * A_12CO + LR_13 * A_CO13 
    ! as an initial point, use the mean value of the diff and A_12CO
    ! LR_init = <diff> / <A_co>
    ! Set LR_13 so that the the effective starting diff in LR12 and LR13 is approx. 20:1
    ! Then we loop over all bands and LRs, estimating new LRs, and computing the chisq.
    ! loop until the change in chisq drops below some value (i.e. 1.d-5).
    allocate(LR(n_co,n_sky-1), data(nlist), res_maps(nlist,n_sky-1))
    
    LR(:,:)=0.d0 !assume no LR to begin with
    ! assume LR of first co component is <map>/<co_comp1> inside mask
    !LR(1,:)=avg_sky(:)/avg_co(1) 


    ! compute initial residual maps
    res_maps(:,:)=diff_maps(:,:)
    do j = 1,n_sky-1
       do i = 1,n_co
          res_maps(:,j)=res_maps(:,j)-LR(i,j)*co_maps2(:,i)
       end do
    end do

    !calc initial chisq
    chisq_curr = sum(sum(res_maps(:,:)*inv_N(:,:)*res_maps(:,:),dim=1),dim=1)
    converged=.false.

    write(*,*) 'Starting iterating over line ratios'
    chisq_diff = 1.d0
    counter=0
    do while (converged==.false. .and. counter<=500)
       counter = counter + 1
       do j=1,n_co
          do k = 1, n_sky-1
             res_maps(:,k)=diff_maps(:,k)
             do i = 1,n_co
                res_maps(:,k)=res_maps(:,k)-LR(i,k)*co_maps2(:,i)
             end do
             
             chisq_in = sum(res_maps(:,k)*inv_N(:,k)*res_maps(:,k))
             if (isnan(chisq_in)) then
                write(*,*) 'Chisq_in is NaN'
                write(*,*) 'iter',counter,'  CO comp',j,'  diff map',k
                stop
             end if
             data(:)=res_maps(:,k)+LR(j,k)*co_maps2(:,j)
             A = 0.d0
             B = 0.d0
             do i = 1,nlist
                A = A + co_maps2(i,j)*inv_N(i,k)*co_maps2(i,j)
                B = B + co_maps2(i,j)*inv_N(i,k)*data(i)
                !if (isnan(B)) then
                !   write(*,*) 'B is NaN. CO_map(pix) =',co_maps2(i,j),'  inv_N(pix) =', &
                !        & inv_N(i,k),'  data(pix) =',data(i)
                !   write(*,*) 'iter',counter,'  CO comp',j,'  diff map',k
                !   write(*,*) 'listpix nr',i,'  pix nr',listpix(i)
                !   write(*,*) 'data(:)=res_maps(:,k)+LR(j,k)*co_maps2(:,j)'
                !   write(*,*) 'res_maps(i,k)=',res_maps(i,k),'LR(j,k)=',LR(j,k),&
                !        & 'co_maps2(i,j)=',co_maps2(i,j)                 
                !   stop
                !end if
             end do

             if (A > 0.d0) then
                temp_LR = B/A
             else
                temp_LR = 0.d0
             end if
             data(:)=data(:)-temp_LR*co_maps2(:,j) !new res_map
             chisq_out = sum(data(:)*inv_N(:,k)*data(:),dim=1)
             if (isnan(chisq_out)) then
                write(*,*) 'Chisq_out is NaN'
                write(*,*) 'iter',counter,'  CO comp',j,'  diff map',k
                write(*,*) 'A =',A,'  B =',B,'  LR out =',temp_LR
                stop
             elseif (debug) then
                write(*,*) 'iter',counter,'  CO comp',j,'  diff map',k
                write(*,*) 'A =',A,'  B =',B,'  LR out =',temp_LR
             end if
             if (chisq_in < chisq_out) then
                write(*,*) 'chisq increases in LR fit iter',counter,'  CO comp',j,'  diff map',k
                write(*,*) 'chisq_in',chisq_in,'   chisq_out',chisq_out
                write(*,*) 'LR in',LR(j,k),'    LR out',temp_LR
                if (debug) then
                   call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=.false.)
                   temp_map=0.d0
                   do i = 1,nlist
                      temp_map(listpix(i),1) = res_maps(i,k)
                   end do
                   call write_result_map('residual_cut_pre_fit.fits', nside, ordering, header,temp_map(:,:1),double_precision)
                   do i = 1,nlist
                      temp_map(listpix(i),1) = data(i)
                   end do
                   call write_result_map('residual_cut_post_fit.fits', nside, ordering, header,temp_map(:,:1),double_precision)
                end if

             !   stop
             end if
             LR(j,k)=temp_LR
          end do
       end do
       
       res_maps(:,:)=diff_maps(:,:)
       do k = 1,n_sky-1
          do j = 1,n_co
             res_maps(:,k)=res_maps(:,k)-LR(j,k)*co_maps2(:,j)
          end do
       end do

       chisq_new = sum(sum(res_maps(:,:)*inv_N(:,:)*res_maps(:,:),dim=1),dim=1)

    
       if (mod(counter,5)==0) then
          chisq_diff = abs(chisq_curr-chisq_new)/min(chisq_curr,chisq_new)
          if ( chisq_diff < 1.d-7) then
             converged = .true.
          else
             chisq_curr = chisq_new
          end if
       end if
       write(*,*) 'iteration =', counter,' chisq old =',chisq_curr,' chisq =',chisq_new,' delta_chisq =',chisq_diff 

    end do

    write(*,*) 'converged line ratios, writing to file CO_LR_diffmaps.dat'

    open(1,file='CO_LR_diffmaps.dat')
    write(1,*) '# CO_comp     diff_map     LR'
    write(*,*) '# CO_comp     diff_map     LR'
    do j = 1,n_co
       do k = 1,n_sky-1
          write(1,*) j, k, LR(j,k)
          write(*,*) j, k, LR(j,k)
       end do
    end do
    close(1)

    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=.false.)
    write(*,*) 'writing residual maps'
    do k = 1,n_sky-1
       temp_map(:,1)=diff_maps_full(:,k)
       do j = 1,n_co
          temp_map(:,1)=temp_map(:,1)-LR(j,k)*co_maps(:,j)
       end do

       temp_int = k
       do i = 1, len(i_text)
          n = temp_int / 10**(len(i_text)-i)
          write(i_text(i:i),'(I1)') n
          temp_int = temp_int - n * 10**(len(i_text)-i)
       end do

       filename='residual_diff' // trim(i_text) // '.fits'
       call write_result_map(trim(filename), nside, ordering, header,temp_map(:,:1),double_precision)
    end do
  end subroutine fit_CO_lineratio



end module map_editor
