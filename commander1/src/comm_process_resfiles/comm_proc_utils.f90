module comm_proc_utils
  use pix_tools
  use fitstools
  use head_fits
  use healpix_types
  use math_tools
  use alm_tools
  use rngmod
  use comm_like_utils
  use sort_utils
  implicit none

  ! *********************************************************************
  ! *  comm_proc_mod -- Support module for Commander processing utils   *
  ! *                                                                   *
  ! *          H. K. Eriksen, E. Gjerløw (University of Oslo)           *
  ! *                                                                   *   
  ! *                               and                                 *
  ! *                                                                   *   
  ! *                       I. K. Wehus  (JPL)                          *   
  ! *                                                                   *
  ! *  History:                                                         *
  ! *      January 10th, 2014 -- First fully functional version         *
  ! *                                                                   *
  ! *********************************************************************

!!$  integer, parameter, private :: i4b = selected_int_kind(9)
!!$  integer, parameter, private :: sp  = selected_real_kind(5,30)
!!$  integer, parameter, private :: dp  = selected_real_kind(12,200)
!!$  integer, parameter, private :: lgt = kind(.true.)
!!$
!!$  real(dp), parameter, private :: pi = 3.141592653589793238462643383279502d0

contains


  ! ************************************************************************
  !                        Utility routines
  ! ************************************************************************

  ! Routine for writing the results into a FITS-file in HKE's standard result file format
  subroutine write_resfile(filename, lmax, nspec, numsamp, numchain, funcname, data)
    implicit none

    character(len=*),             intent(in) :: filename, funcname
    integer(i4b),                 intent(in) :: lmax, nspec, numsamp, numchain
    real(sp), dimension(:,:,:,:), intent(in) :: data

    integer(i4b)         :: status, blocksize, readwrite, unit
    integer(i4b)         :: nelements, fpixel, group, bitpix, naxis
    logical(lgt)         :: simple, extend, anyf, exist
    real(sp)             :: nullval
    character(len=80)    :: comment, errorline
    integer(i4b), dimension(4) :: naxes

    unit = 21

    ! Output to a fits-file
    status = 0
    bitpix = -32
    simple = .true.
    naxis = 4
    extend = .true.

    naxes(1) = size(data(:,1,1,1))
    naxes(2) = size(data(1,:,1,1))
    naxes(3) = size(data(1,1,:,1))
    naxes(4) = size(data(1,1,1,:))

    inquire(file=trim(filename),exist=exist)
    if (exist) call system('mv ' // trim(filename) // ' ' // trim(filename) // '_old')

    ! Open a new fitsfile
    blocksize = 1
    call ftinit(unit,trim(filename),blocksize,status)

    ! Write the required headers
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    ! Write a few keywords
    call ftpkys(unit,'FUNCNAME',trim(funcname),'Full function name',status)
    call ftpkyj(unit,'LMAX',lmax,'Maximum multipole moment',status)
    call ftpkyj(unit,'NUMSAMP',numsamp,'Number of samples',status)
    call ftpkyj(unit,'NUMCHAIN',numchain,'Number of independent chains',status)
    call ftpkyj(unit,'NUMSPEC',nspec,'Number of power spectra',status)


    ! Output the array
    group=1
    fpixel=1
    nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
    call ftppre(unit,group,fpixel,nelements,data,status)


    ! Write the current date
    call ftpdat(unit,status)

    ! close the file and free the unit number
    call ftclos(unit, status)
    call ftfiou(unit, status)

    if (status /= 0) then
       write(*,*) 'Error in FITS operations -- status = ', status

       call ftgmsg(errorline)
       do while (trim(errorline) /= '')
          write(*,*) errorline
          call ftgmsg(errorline)
       end do

    end if

  end subroutine write_resfile


  ! Routine for reading the Gibbs sigma samples 
  subroutine read_resfile(filename, unit, numspec, lmax, numchains, numsamples, data)
    implicit none

    character(len=128),                           intent(in)  :: filename
    integer(i4b),                                 intent(in)  :: unit
    integer(i4b),                                 intent(out) :: numspec, lmax, numchains, numsamples
    real(sp),         pointer, dimension(:,:,:,:)             :: data

    integer(i4b)         :: l, status, blocksize, readwrite
    integer(i4b)         :: fpixel, group, numargs
    logical(lgt)         :: anyf
    real(sp)             :: nullval
    character(len=80)    :: comment

    status = 0
    readwrite = 0
    nullval = 0.

    numargs = 1

    ! Open the result file
    call ftopen(unit,trim(filename),readwrite,blocksize,status)

    ! Read keywords
    call ftgkyj(unit,'LMAX',     lmax,       comment,status)
    call ftgkyj(unit,'NUMSAMP',  numsamples, comment,status)
    call ftgkyj(unit,'NUMCHAIN', numchains,  comment,status)
    call ftgkyj(unit,'NUMSPEC',  numspec,    comment,status)

    allocate(data(0:lmax,numspec,numchains,numargs+numsamples))

    ! Read the binned power spectrum array
    group  = 1
    fpixel = 1
    call ftgpve(unit,group,fpixel,size(data),nullval,data,anyf,status)

    call ftclos(unit,status)

  end subroutine read_resfile


  ! Routine for writing the results into a FITS-file in HKE's standard result file format
  subroutine write_fgt_resfile(filename, data)
    implicit none

    character(len=*),             intent(in) :: filename
    real(sp), dimension(:,:,:,:), intent(in) :: data

    integer(i4b)         :: status, blocksize, readwrite, unit
    integer(i4b)         :: nelements, fpixel, group, bitpix, naxis
    logical(lgt)         :: simple, extend, anyf, exist
    real(sp)             :: nullval
    character(len=80)    :: comment, errorline

    integer(i4b), dimension(4) :: naxes

    unit = 21

    ! Output to a fits-file
    status = 0
    bitpix = -32
    simple = .true.
    naxis = 4
    extend = .true.

    naxes(1) = size(data(:,1,1,1))
    naxes(2) = size(data(1,:,1,1))
    naxes(3) = size(data(1,1,:,1))
    naxes(4) = size(data(1,1,1,:))

    inquire(file=trim(filename),exist=exist)
    if (exist) call system('mv ' // trim(filename) // ' ' // trim(filename) // '_old')

    ! Open a new fitsfile
    blocksize = 1
    call ftinit(unit,trim(filename),blocksize,status)

    ! Write the required headers
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    ! Write a few keywords
    call ftpkys(unit,'FUNCNAME','Gibbs sampled foreground template weights','Full function name',status)
    call ftpkyj(unit,'NUMTEMP',naxes(1),'Number of templates fitted',status)
    call ftpkyj(unit,'NUMBAND',naxes(2),'Number of bands',status)
    call ftpkyj(unit,'NUMCHAIN',naxes(3),'Number of chains',status)
    call ftpkyj(unit,'NUMITER',naxes(4),'Number of iterations',status)

    ! Output the array
    group=1
    fpixel=1
    nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
    call ftppre(unit,group,fpixel,nelements,data,status)

    ! Write the current date
    call ftpdat(unit,status)

    ! close the file and free the unit number
    call ftclos(unit, status)
    call ftfiou(unit, status)

    if (status /= 0) then
       write(*,*) 'Error in FITS operations -- status = ', status

       call ftgmsg(errorline)
       do while (trim(errorline) /= '')
          write(*,*) errorline
          call ftgmsg(errorline)
       end do

    end if

  end subroutine write_fgt_resfile



  subroutine write_map3(filename, map, ordering)
    implicit none

    character(len=*),                   intent(in)  :: filename
    real(dp),         dimension(0:,1:), intent(in)  :: map
    integer(i4b),                       intent(in), optional :: ordering

    integer(i4b)   :: npix, nlheader, nmaps, i, nside, order
    logical(lgt)   :: exist, polarization

    character(len=80), dimension(1:120)    :: header

    npix         = size(map(:,1))
    nside        = nint(sqrt(real(npix,sp)/12.))
    nmaps        = size(map(0,:))
    polarization = (nmaps == 3)
    order        = 1; if (present(ordering)) order = ordering


    inquire(file=trim(filename),exist=exist)
    if (exist) call system('rm ' // trim(filename))


    !-----------------------------------------------------------------------
    !                      write the map to FITS file
    !  This is copied from the synfast.f90 file in the Healpix package
    !-----------------------------------------------------------------------

    nlheader = SIZE(header)
    do i=1,nlheader
       header(i) = ""
    enddo

    ! start putting information relative to this code and run
    call add_card(header)
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Sky Map Pixelisation Specific Keywords    ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"PIXTYPE","HEALPIX","HEALPIX Pixelisation")
    if (order == 1) then
       call add_card(header,"ORDERING","RING",  "Pixel ordering scheme, either RING or NESTED")
    else
       call add_card(header,"ORDERING","NEST",  "Pixel ordering scheme, either RING or NESTED")
    end if
    call add_card(header,"NSIDE"   ,nside,   "Resolution parameter for HEALPIX")
    call add_card(header,"FIRSTPIX",0,"First pixel # (0 based)")
    call add_card(header,"LASTPIX",npix-1,"Last pixel # (0 based)")
    call add_card(header,"BAD_DATA",  HPX_DBADVAL ,"Sentinel value given to bad pixels")
    call add_card(header) ! blank line

    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Planck Simulation Specific Keywords      ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"EXTNAME","Commander Gibbs sample")
    call add_card(header,"POLCCONV","COSMO"," Coord. convention for polarisation (COSMO/IAU)")

    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Data Description Specific Keywords       ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"INDXSCHM","IMPLICIT"," Indexing : IMPLICIT or EXPLICIT")
    call add_card(header,"GRAIN", 0, " Grain of pixel indexing")
    call add_card(header,"COMMENT","GRAIN=0 : no indexing of pixel data                         (IMPLICIT)")
    call add_card(header,"COMMENT","GRAIN=1 : 1 pixel index -> 1 pixel data                     (EXPLICIT)")
    call add_card(header,"COMMENT","GRAIN>1 : 1 pixel index -> data of GRAIN consecutive pixels (EXPLICIT)")
    call add_card(header) ! blank line
    call add_card(header,"POLAR",polarization," Polarisation included (True/False)")
    call add_card(header,"DERIV",0," Derivative included (0, 1 or 2)")

    call add_card(header) ! blank line
    call add_card(header,"TTYPE1", "TEMPERATURE","Temperature map")
    call add_card(header,"TUNIT1", 'muK',"map unit")
    call add_card(header)

    if (polarization) then
       call add_card(header,"TTYPE2", "Q-POLARISATION","Q Polarisation map")
       call add_card(header,"TUNIT2", 'muK',"map unit")
       call add_card(header)

       call add_card(header,"TTYPE3", "U-POLARISATION","U Polarisation map")
       call add_card(header,"TUNIT3", 'muK',"map unit")
       call add_card(header)
    endif
    call add_card(header,"COMMENT","*************************************")

    call write_bintab(map, npix, nmaps, header, nlheader, filename)

  end subroutine write_map3

  ! Small utility for converting an integer to a string
  subroutine int2string(integer, string)
    implicit none

    integer(i4b),     intent(in)  :: integer
    character(len=*), intent(out) :: string

    integer(i4b)               :: temp_int, i, k

    temp_int = integer
    do i = 1, len(string)
       k = temp_int / 10**(len(string)-i)
       write(string(i:i),'(I1)') k
       temp_int = temp_int - k * 10**(len(string)-i)
    end do

  end subroutine int2string


  subroutine read_pixwin(nside, nmaps, pixwin)
    implicit none

    integer(i4b),                          intent(in)  :: nside, nmaps
    real(dp),     pointer, dimension(:,:)              :: pixwin

    integer(i4b)        :: nc
    character(len=128)  :: pixwin_file
    character(len=4)    :: nside_text
    logical(lgt)        :: exist, anynull
    real(dp)            :: nullval

    allocate(pixwin(0:4*nside,nmaps))

    if (nmaps == 3) then
       nc = 2
    else
       nc = 1
    end if

    call int2string(nside, nside_text)
    pixwin_file = 'pixel_window_n' // nside_text // '.fits'
    inquire(file=pixwin_file, exist=exist)
    if (exist) then
       nullval = 0.d0
       call read_dbintab(pixwin_file, pixwin(0:4*nside,1:nc), 4*nside+1, nc, nullval, anynull)
       if (nmaps == 3) pixwin(:,3) = pixwin(:,2)
    else
       pixwin = 1.d0
       write(*,*) ''
       write(*,*) 'Warning! Pixel window file ', trim(pixwin_file), ' not found. '
       write(*,*) 'Using unity weights.'
       write(*,*) ''
    end if

  end subroutine read_pixwin

  subroutine read_ringweights(nside, weights)
    implicit none

    integer(i4b),             intent(in)  :: nside
    real(dp), dimension(:,:), intent(out) :: weights

    character(len=128)  :: weight_file
    character(len=5)    :: nside_text
    logical(lgt)        :: exist, anynull
    integer(i4b)        :: nmaps
    real(dp)            :: nullval

    nmaps = size(weights(1,:))

    call int2string(nside, nside_text)
    weight_file = 'weight_ring_n' // nside_text // '.fits'
    inquire(file=weight_file, exist=exist)
    if (exist) then
       nullval = 0.d0
       call read_dbintab(weight_file, weights, 2*nside, nmaps, nullval, anynull)
       weights = 1.d0 + weights
    else
       weights = 1.d0
       write(*,*) ''
       write(*,*) 'Warning! Weight file ', trim(weight_file), ' not found. '
       write(*,*) 'Using unity weights in the spherical harmonic transforms.'
       write(*,*) ''
    end if

  end subroutine read_ringweights

  subroutine read_map(filename, map)
    implicit none

    character(len=128),                 intent(in)  :: filename
    real(dp),         dimension(0:,1:), intent(out) :: map

    integer(i4b)   :: nside, nmaps, ordering, int_npix, temp_i, i, npix, nmaps_in, nside_in
    real(dp)       :: nullval
    logical(lgt)   :: anynull

    npix  = size(map(:,1))
    nmaps = size(map(0,:))
    nside = nint(sqrt(real(npix,sp)/12.))

    temp_i = getsize_fits(trim(filename), ordering=ordering, nside=nside_in, nmaps=nmaps_in)
    if ((nmaps_in < nmaps) .or. (nside_in /= nside)) then
       write(*,*) 'Incorrect nside or nmaps for file called ', trim(filename)
    end if

    call read_bintab(trim(filename), map, npix, nmaps, nullval, anynull)

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(:,i))
       end do
    end if

  end subroutine read_map


  subroutine read_beam(beamfile, beam)
    implicit none

    character(len=128),                   intent(in)  :: beamfile
    real(dp),           dimension(0:,1:), intent(out) :: beam

    integer(i4b) :: l, lmax, nmaps
    real(dp)     :: sigma_sq
    real(dp),          allocatable, dimension(:,:)   :: inbeam
    character(len=80),              dimension(1:180) :: header
    
    lmax  = size(beam(:,1))-1
    nmaps = size(beam(0,:))

    ! Seem to remember there is something weird going on with the WMAP beams when reading only 
    ! one component at a time. Remove this wrapper if you feel more comfortable with that...
    allocate(inbeam(0:lmax,4))

!    call read_asctab(beamfile, inbeam, lmax, 4, header)
    call fits2cl(beamfile, inbeam, lmax, 4, header)

    beam(0:lmax,1:nmaps) = inbeam(0:lmax,1:nmaps)

    if (nmaps > 1) then
       if (sum(beam(:,2)) < 1.d0) beam(:,2) = beam(:,1)
       if (sum(beam(:,3)) < 1.d0) beam(:,3) = beam(:,2)
    end if

    deallocate(inbeam)

  end subroutine read_beam

   subroutine comm_read_covmat(covfile, nmaps, cov, scale)
      implicit none
   
      character(len=*),                 intent(in)  :: covfile
      integer(i4b),                     intent(in)  :: nmaps
      real(dp),         dimension(:,:), intent(out) :: cov
      real(dp),                         intent(out) :: scale
   
      integer(i4b) :: i, j, k, cov_order, unit, nside, npix, n, n_p
      character(len=4) :: suffix
      real(dp), allocatable, dimension(:)   :: buffer
      real(dp), allocatable, dimension(:,:) :: map
   
      unit   = comm_getlun()
      n      = len(trim(covfile))
      n_p    = size(cov,1)
      suffix = covfile(n-3:n)
   
      if (suffix == '.unf') then
         ! Assume unformatted F90 file with appropriate header entries
         scale = 1.d0 ! Assume data are in uK
         open(unit,file=trim(covfile),form='unformatted')
            read(unit) i
         read(unit) cov_order
         read(unit) i
         do i = 1, n_p
            read(unit) cov(i,:)
         end do
         close(unit)
      else if (suffix == '.dat') then
         ! Assume binary C format
         npix      = size(cov,1) / nmaps
         nside     = nint(sqrt(real(npix,dp)/12.d0))
         scale     = 1.d6 ! Assume data are in K
         cov_order = 2    ! Assume covariance matrix is in NESTED format
         open(unit,file=trim(covfile),form='unformatted',access='direct', recl=8*n_p)
         do i = 1, n_p
            read(unit,rec=i) cov(i,:)
         end do
         close(unit)
         cov = cov * scale**2
      else if (suffix == 'fits') then
         ! Assume binary C format
         npix      = size(cov,1) / nmaps
         nside     = nint(sqrt(real(npix,dp)/12.d0))
         scale     = 1.d0 ! Assume data are in uK
         cov_order = 1    ! Map is always converted to RING format
         allocate(map(0:npix-1,nmaps))
         call read_map(covfile, map)
         cov = 0.d0
         k = 1
         do j = 1, nmaps
            do i = 0, npix-1
               cov(k,k) = map(i,j)**2
               k        = k+1
            end do
         end do
         cov = cov * scale**2
      else
         write(*,*) 'Unsupported covariance matrix type = ', suffix
         stop
      end if
      if (cov_order == 2) then
         npix   = size(cov,1) / nmaps
         nside  = nint(sqrt(real(npix,dp)/12.d0))
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
   
   end subroutine comm_read_covmat

  subroutine comm_output_powspec(filename, bins, cls, sigma)
    implicit none

    character(len=*),                   intent(in) :: filename
    integer(i4b),     dimension(1:,1:), intent(in) :: bins
    real(dp),         dimension(1:),    intent(in) :: cls
    real(dp),         dimension(1:,1:), intent(in) :: sigma

    integer(i4b) :: i, j, unit, npar

    unit = comm_getlun()
    npar = size(bins,1)

    ! Output results
    open(unit,file=trim(filename), recl=1024)
    write(unit,fmt='(a)',advance='no') '#  l_c   lmin  lmax      TT         +dTT        -dTT      '
    write(unit,fmt='(a)',advance='no') 'TE      +dTE     -dTE       TB     +dTB     -dTB       '
    write(unit,fmt='(a)',advance='no') 'EE     +dEE     -dEE       EB     +dEB     -dEB       '
    write(unit,fmt='(a)')              'BB     +dBB     -dBB '
    j = 1
    do while (j <= npar)
       write(unit,fmt='(f6.1,2i6)',advance='no') 0.5*(bins(j,1)+bins(j,2)),bins(j,1), bins(j,2)
       do i = 1, 6
          if (i == 1) then
             if (bins(j,3) /= i) then
                write(unit,fmt='(3f12.3)',advance='no') 0.d0, 0.d0, 0.d0
             else
                write(unit,fmt='(3f12.3)',advance='no') cls(j), sigma(j,1), sigma(j,2)
                j = j+1
             end if
          else
             if (bins(j,3) /= i) then
                write(unit,fmt='(3f12.6)',advance='no') 0.d0, 0.d0, 0.d0
             else
                write(unit,fmt='(3f12.6)',advance='no') cls(j), sigma(j,1), sigma(j,2)
                j = j+1
             end if
          end if
       end do
       write(unit,*)
    end do
    close(unit)

  end subroutine comm_output_powspec

  subroutine dump_matrix(mat, fname, unit, fmt, idx)
    implicit none
    real(dp),         intent(in)           :: mat(:,:)
    character(len=*), intent(in), optional :: fname, fmt
    integer(i4b),     intent(in), optional :: unit
    logical(lgt),     intent(in), optional :: idx
    character(len=256)                     :: fmt_
    integer(i4b)                           :: i, j, unit_
    logical(lgt)                           :: doclose, idx_
    fmt_ = '(e12.4)'; if(present(fmt)) fmt_ = fmt
    idx_ = .false.;   if(present(idx)) idx_ = idx
    doclose = .false.
    unit_ = 6
    if(present(unit)) then
       unit_ = unit
    elseif(present(fname)) then
       unit_ = comm_getlun()
       open(unit_,file=fname)
       doclose = .true.
    end if
    do i = 1, size(mat,1)
       if(idx_) write(unit_,fmt="(i9)",advance="no") i
       do j = 1, size(mat,2)
          write(unit_,fmt=fmt_,advance="no") mat(i,j)
       end do
       write(unit_,*)
    end do
    if(doclose) then
       close(unit_)
    end if
  end subroutine

  subroutine compute_histogram(samples, n, x, P)
    implicit none

    real(dp),                  dimension(:), intent(in)  :: samples
    integer(i4b),                            intent(out) :: n
    real(dp),     allocatable, dimension(:), intent(out) :: x, P

    integer(i4b) :: i, j, nsamp
    real(dp)     :: xmin, xmax, dx
    real(dp), allocatable, dimension(:) :: sorted

    nsamp = size(samples)
    n     = min(nint(sqrt(1.d0*nsamp)), 100)

    allocate(sorted(nsamp))
    sorted = samples
    call QuickSort_real(sorted)
    xmin = sorted(max(nint(0.005*nsamp),1))
    xmax = sorted(min(nint(0.995*nsamp),nsamp))
    dx   = (xmax-xmin) / n

    allocate(x(n), P(n))
    x = 0.d0
    P = 0.d0
    do i = 1, nsamp
       j = max(min(floor((samples(i)-xmin)/dx)+1,n),1)
       P(j) = P(j)+1.d0
    end do
    P = P / (n*dx)
    do i = 1, n
       x(i) = xmin + (i-0.5d0)*dx
    end do
    
    deallocate(sorted)

  end subroutine compute_histogram
    

end module comm_proc_utils
