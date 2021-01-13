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
! generic body of the subroutines remove_dipole_*
! to be inserted as is in pix_tools.f90
!
!--------------------------------------------------------------
    ! dummy
    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree
    real   (kind=DP),   dimension(0:),  intent(out)   :: multipoles
!    real   (kind=DP),                   intent(in)    :: cos_theta_cut
    real   (kind=DP),   dimension(1:),  intent(in)    :: zbounds
    real   (kind=KMAP), dimension(0:),  intent(inout) :: map
    real   (kind=KMAP),                 intent(in), optional :: fmissval
    real   (kind=KMAP), dimension(0:),  intent(in), optional :: mask
    ! local
    real   (kind=KMAP)                :: fmiss_effct
    integer(kind=i4b)                 :: ipix, npix, n_mult
    integer(kind=i4b)                 :: i
    logical(lgt)                      :: dodipole, do_mask
    real(kind=dp)                     :: flag, temp, wmin
    real(kind=dp), dimension(1:3)     :: vec
    real(kind=dp), dimension(0:3)     :: b, wdiag
    real(kind=dp), dimension(0:3,0:3) :: mat, umat, vmat
    real(kind=dp), dimension(0:3)     :: dmultipoles
    real(kind=dp)                     :: theta1, theta2
    integer(kind=i4b)                 :: ncpix, ncp, nbad
    integer(kind=i4b), dimension(:), allocatable :: cut_pixel
    !============================================================

    npix = nside2npix(nside)
    multipoles = 0.0_dp
    n_mult = (degree)**2
    fmiss_effct = fbad_value
    if (present(fmissval)) fmiss_effct = fmissval
    do_mask = .false.
    if (present(mask)) do_mask = (size(mask) == size(map))

    if (degree == 0) then
       print*," No monopole nor dipole removal"
       return
    elseif (degree == 1) then
       dodipole = .false.
    else if (degree == 2) then
       dodipole = .true.
    else
       print*,code//"> degree can only be "
       print*,"      1: monopole (l=0) removal or "
       print*,"      2: monopole and dipole (l=0,1) removal"
       print*,code//"> ABORT ! "
       call fatal_error
    endif

    !----------------------------------------------
    ! flag out pixels excluded by cut
    !----------------------------------------------
    if (zbounds(2) <= zbounds(1)) then ! remove one strip
       ncpix = npix * (zbounds(1)-zbounds(2))/2. * 1.1 + 10* nside
    else ! remove 2 strips
       ncpix = npix * (2.0 + zbounds(1)-zbounds(2))/2 * 1.1 + 10*nside
    endif
    allocate(cut_pixel(0:ncpix))
    theta1 = acos(zbounds(1)) ; theta2 = acos(zbounds(2))
    call query_strip(nside, theta1, theta2, cut_pixel, ncp, nest=ordering-1, inclusive=0)
    map(cut_pixel(0:ncp-1)) = fmiss_effct
    deallocate(cut_pixel)

    !----------------------------------------------
    ! generate least square linear system
    !----------------------------------------------
    mat = 0.0_dp
    b   = 0.0_dp
    nbad = 0
    do ipix = 0, npix-1

       ! flag = 1 for good values
       !      = 0 for bad values
       flag = 1.0_dp
       if ( abs(map(ipix) - fmiss_effct) <= abs(1.e-5*fmiss_effct) ) then
          nbad = nbad + 1
          goto 20
       endif
       if (do_mask) then
          if (abs(mask(ipix)) <= 1.e-10) then
             nbad = nbad + 1
             goto 20
          endif
       endif

       if (dodipole) then
          ! computes dipole basis functions
          ! pixel -> vector
          if (ordering == 1) call pix2vec_ring( nside, ipix, vec)
          if (ordering == 2) call pix2vec_nest( nside, ipix, vec)
       endif

       ! construct vector T*(1,x,y,z)
       temp = map(ipix) * flag
       b(0) = b(0) + temp
       if (dodipole) then
          ! computes dipole basis functions
          b(1:3) = b(1:3) + temp * vec(1:3)
       endif

       ! construct matrix (1,x,y,z)#(1,x,y,z)
       mat(0,0) = mat(0,0) + flag
       if (dodipole) then
          do i = 1, 3
             mat(i,0)   = mat(i,0)   + vec(i) * flag
             mat(i,1:3) = mat(i,1:3) + vec(i) * vec(1:3) * flag
          enddo
       endif

20     continue
    enddo
    print*,'Excluding ',nbad,' pixels when computing monopole and/or dipole'

    ! first row = first column (and vice versa)
    mat(0,1:3) = mat(1:3,0)


    !----------------------------------------------
    ! solve system    mat . (mono, dip_x, dip_y, dip_z) = b
    !----------------------------------------------

    if (dodipole) then
       ! SVD decomposition
       umat = mat
       call dsvdcmp(umat, 4, 4, 4, 4, wdiag, vmat)
       ! thresholding
       wmin = maxval(wdiag)* 1.e-6_dp
       where (wdiag < wmin)
          wdiag = 0.0_dp
       end where
       ! back substitution
       call dsvbksb(umat, wdiag, vmat, 4, 4, 4, 4, b, dmultipoles)

    else
       dmultipoles(0) = b(0) / mat(0,0) ! average temperature
    endif

    !----------------------------------------------
    ! remove monopole and dipole
    !----------------------------------------------
    do ipix = 0, npix-1

!        if ( abs(map(ipix) - fmiss_effct) < 1.e-5*fmiss_effct ) goto 10
       if ( abs(map(ipix) - fmiss_effct) < abs(1.e-5*fmiss_effct) ) goto 10

       map(ipix) = map(ipix) - dmultipoles(0)
       if (dodipole) then
          ! computes dipole basis functions
          ! pixel -> vector
          if (ordering == 1) call pix2vec_ring( nside, ipix, vec)
          if (ordering == 2) call pix2vec_nest( nside, ipix, vec)
          map(ipix) = map(ipix) - SUM( dmultipoles(1:3) * vec(1:3))
       endif

10     continue

    enddo

    multipoles(0:n_mult-1) = dmultipoles(0:n_mult-1)

    return
