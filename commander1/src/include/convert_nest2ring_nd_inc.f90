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
! generic body of the subroutines convert_nest2ring_*_nd 
! to be inserted as is in pix_tools.f90
!
!--------------------------------------------------------------
    !=======================================================================
    !     makes the conversion NEST to RING
    !=======================================================================
    integer(kind=I4B), intent(in) :: nside
    integer(kind=I4B) :: nd, j
    integer(kind=I4B) :: npix, ipn, ipr
    integer(kind=I4B),  dimension(:), allocatable :: mapping
    !=======================================================================

    npix = nside2npix(nside)
    call assert (npix>0,       code//": invalid Nside")

    nd = size(map,2)

    ! case where 2nd dimension=1
    if (nd == 1) then
       map1 => map(0:npix-1,1)
       call convert_nest2ring(nside, map1)
       return
    endif

    allocate(map_tmp(0:npix-1))
    allocate(mapping(0:npix-1))

    ! do N->R mapping only once
!$OMP parallel default(none) &
!$OMP   shared(mapping, npix, nside) &
!$OMP   private(ipr, ipn)
!$OMP do schedule(dynamic,64)
    do ipn = 0, npix-1
       call nest2ring(nside, ipn, ipr)
       mapping(ipn) = ipr
    enddo
!$OMP end do
!$OMP end parallel

    ! convert maps one by one
    do j = 1, nd
       do ipn = 0, npix-1
          map_tmp(mapping(ipn)) = map(ipn,j)
       enddo
       map(0:npix-1, j) = map_tmp(0:npix-1)
    enddo

    deallocate(map_tmp)
    deallocate(mapping)

    return
