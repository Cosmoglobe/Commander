!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon, 
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------
!--------------------------------------------------------------
!
! generic body of the subroutines convert_inplace_*_nd 
! to be inserted as is in pix_tools.f90
!
!--------------------------------------------------------------
    integer(kind=I1B), parameter :: TODO = 0_I1B, DONE = 1_I1B
    !
    external  subcall ! required by some compilers (DEC, Portland, ...)
    integer(kind=i4b)                       :: nside, nd, j
    integer(kind=i4b)                       :: npix, ilast, i1, i2
    integer(kind=I1B), dimension(:), allocatable::check
    !------------------------------------------------------------------
    npix=size(map,1)
    nd  =size(map,2)
    nside=npix2nside(npix)
    call assert (nside<=NS_MAX, code//": map too big")
    call assert (nside>0,       code//": invalid Nside")
    call assert (nd<=ND_MAX,    code//": map 2nd dim too large")

    print*, "Convert: Converting map pixelisation scheme"
    allocate(check(0:npix-1))
    check = TODO

    ilast=0                   !start at first pixel
    do
       pixbuf2(1:nd) = map(ilast,1:nd)      !initialise
       i1 = ilast
       call subcall(nside, i1, i2)
       do
          if (check(i2) == DONE) exit
          do j = 1, nd ! should be faster than 2 separate (implicit) loops
             pixbuf1(j) = map(i2,j)
             map(i2,j)  = pixbuf2(j)
          enddo
          check(i2) = DONE
          pixbuf2(1:nd) = pixbuf1(1:nd)
          i1 = i2
          call subcall(nside, i1, i2)
       enddo
       do
          if (.not. (check(ilast)==DONE .and. ilast<npix-1)) exit ! npix-1 or npix
          ilast = ilast + 1
       enddo
       if(ilast == npix-1) exit ! npix-1 or npix
    enddo
    deallocate(check)
    return

