module mask_tools
  use healpix_types
  use pix_tools
  use fitstools
  implicit none

  ! *******************************************************
  ! *                                                     *
  ! *             Module for handling masks               *
  ! *                                                     *
  ! *   written by Hans Kristian Eriksen, April 30, 2003  *
  ! *                                                     *
  ! *         Copyright 2003. All rights reserved.        *
  ! *                                                     *
  ! *******************************************************

  !
  ! ChangeLog:
  ! 
  ! April 30, 2003 -- Implemented write_mask, read_mask and 
  !                   get_maskparam
  ! 
  !


contains

  ! Routine for writing a mask to a FITS-file
  ! 
  ! fileformat = 1   ===>   straigth, run-length coding 
  !
  subroutine write_mask(unit, maskfile, ordering, mask, fileformat)
    implicit none

    character(len=*),                intent(in)   :: maskfile
    integer(i4b),                    intent(in)   :: unit, ordering, fileformat
    logical(lgt),     dimension(0:), intent(in)   :: mask

    integer(i4b)  :: status, blocksize, group, fpixel, included
    integer(i4b)  :: naxis, bitpix, i, j, num_true, npix, nside, numinc, numval
    logical(lgt)  :: simple, extend, current
    integer(i4b),              dimension(1) :: naxes
    integer(i4b), allocatable, dimension(:) :: pix_num

    npix  = size(mask)
    nside = nint(sqrt(real(npix,sp)/12.))

    if (fileformat == 1) then

       allocate(pix_num(npix))
       pix_num = 0

       numinc = 1
       current = .true.

       i = 0
       do while (i < npix)
          
          if (((mask(i)) .and. (.not. current)) .or. ((.not. mask(i)) .and. (current))) then
!          if (mask(i) /= current) then
             numinc = numinc + 1
             current = .not. current
          end if

          pix_num(numinc) = pix_num(numinc) + 1

          i = i+1

       end do

    end if

    numval = 0
    do i = 0, npix-1
       if (mask(i)) numval = numval+1
    end do

    !Output the reduced mask to a FITS-file
    status = 0
!    unit = 150
    bitpix = 32
    simple = .true.
    naxis = 1
    extend = .true.
    naxes(1) = numinc
    group = 1
    fpixel = 1
    blocksize = 1

    ! Open a new file with given filename
    call ftinit(unit,trim(maskfile),blocksize,status)
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    ! Write the current date
    call ftpdat(unit,status)

    ! Output the data array
    call ftpprj(unit,group,fpixel,numinc,pix_num(1:numinc),status)

    ! And finally som useful keywords
    call ftpkyj(unit,'NSIDE',nside,'',status)
    call ftpkyj(unit,'NUMPIX',numinc,'Number of listed values',status)
    call ftpkyj(unit,'ORDERING',ordering,'RING=1, NESTED=2',status)
    call ftpkyj(unit,'NUMVAL',numval,'Number of .true. pixels',status)
    call ftpkyj(unit,'FFORMAT',fileformat,&
         & 'Compression status',status)

    ! Close the file
    call ftclos(unit,status)
    call ftfiou(unit,status)

    if (status .gt. 0) then
       call printerror(status)
    end if

    deallocate(pix_num)

  end subroutine write_mask

  ! Routine for reading a mask from a FITS-file
  subroutine read_mask(unit, maskfile, mask)
    implicit none

    integer(i4b),                intent(in)   :: unit
    character(len=*),            intent(in)   :: maskfile
    logical(lgt), dimension(0:), intent(out)  :: mask

    integer(i4b)      :: status, blocksize, group, fpixel, nullval
    integer(i4b)      :: bitpix, i, j, k, npix, readwrite, included, fileformat
    logical(lgt)      :: anyf, current
    character(len=80) :: comment

    integer(i4b), allocatable, dimension(:)  :: pix_num

    status = 0
!    unit = 150
    readwrite = 0
    nullval = -1

    ! Open the mask file
    call ftopen(unit,trim(maskfile),readwrite,blocksize,status)

    ! Read the necessary keywords
    call ftgkyj(unit,'NUMPIX',npix,comment,status)
    call ftgkyj(unit,'FFORMAT',fileformat,comment,status)

    ! Read the actual data
    group = 1
    fpixel = 1
    allocate(pix_num(npix))
    call ftgpvj(unit,group,fpixel,npix,nullval,pix_num,anyf,status)

    call ftclos(unit,status)
    
    if (fileformat == 1) then

       current = .true.
       j = 0

       do i = 1, npix

          do k = 1, pix_num(i)
             mask(j) = current
             j = j+1
          end do
          
          current = .not. current

       end do

    end if
    deallocate(pix_num)

    if (status .gt. 0) then
       call printerror(status)
    end if

  end subroutine read_mask


  ! Routine for getting parameters from a mask FITS-filexs
  subroutine get_maskparam(unit, maskfile, nside, ordering, numval)
    implicit none

    integer(i4b),       intent(in)            :: unit
    character(len=128), intent(in)            :: maskfile
    integer(i4b),       intent(out)           :: nside, ordering
    integer(i4b),       intent(out), optional :: numval

    integer(i4b)      :: status, blocksize
    integer(i4b)      :: readwrite
    character(len=80) :: comment

    status = 0
!    unit = 150
    readwrite = 0

    ! Open the mask file
    call ftopen(unit,trim(maskfile),readwrite,blocksize,status)

    ! Read the necessary keywords
    call ftgkyj(unit,'NSIDE',   nside,   comment,status)
    call ftgkyj(unit,'ORDERING',ordering,comment,status)

    if (present(numval)) then
       call ftgkyj(unit,'NUMVAL',numval,comment,status)
    end if

    call ftclos(unit,status)

    if (status .gt. 0) then
       call printerror(status)
    end if

  end subroutine get_maskparam

  subroutine convert_nest2ring_mask(mask)
    implicit none

    logical(lgt), dimension(0:), intent(inout) :: mask

    integer(i4b) :: npix, nside, i, j
    logical(lgt), allocatable, dimension(:) :: mask_conv

    npix  = size(mask)
    nside = nint(sqrt(real(npix,sp)/12.))

    allocate(mask_conv(0:npix-1))
    do i = 0, npix-1
       call nest2ring(nside,i,j)
       mask_conv(j) = mask(i)
    end do

    mask = mask_conv

    deallocate(mask_conv)

  end subroutine convert_nest2ring_mask

  subroutine convert_ring2nest_mask(mask)
    implicit none

    logical(lgt), dimension(0:), intent(inout) :: mask

    integer(i4b) :: npix, nside, i, j
    logical(lgt), allocatable, dimension(:) :: mask_conv

    npix  = size(mask)
    nside = nint(sqrt(real(npix,sp)/12.))

    allocate(mask_conv(0:npix-1))
    do i = 0, npix-1
       call ring2nest(nside,i,j)
       mask_conv(j) = mask(i)
    end do

    mask = mask_conv

    deallocate(mask_conv)

  end subroutine convert_ring2nest_mask

  subroutine udgrade_mask(inmask, outmask)
    implicit none

    logical(lgt), dimension(0:), intent(in)  :: inmask
    logical(lgt), dimension(0:), intent(out) :: outmask

    integer(i4b) :: i, j, npix_in, npix_out, fac

    npix_in  = size(inmask)
    npix_out = size(outmask)

    if (npix_out > npix_in) then

       fac     = npix_out / npix_in
       outmask = .false.

       do i = 0, npix_in-1
          if (inmask(i)) then
             do j = 0, fac-1
                outmask(fac*i+j) = .true.
             end do
          end if
       end do

    else

       fac     = npix_in / npix_out
       outmask = .true.

       do i = 0, npix_in-1
          if (.not. (inmask(i))) then
             outmask(int(real(i,sp)/real(fac,sp))) = .false.
          end if
       end do

    end if

  end subroutine udgrade_mask


end module mask_tools
