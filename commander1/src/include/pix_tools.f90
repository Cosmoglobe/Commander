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
module pix_tools
!==================================================================
!    Various subroutines for converting angles and unit vectors
!    to pixel numbers in both the NESTED and RING schemes,
!    as well as looking for neighbours and making pixel queries
!
!    last update : Oct 4, 2002, EH, computation of pixel vertex in
!        pix2vec_ring and pix2vec_nest
!    2004-05-28, EH: edition of next_in_line_nest for nside=1
!                    addition of optional 'mask' in remove_dipole
!    2004-06-01      correction of bug at low phi in in_ring
!    2004-08-25      added multi-D convert_inplace_* and convert_nest2ring, convert_ring2nest
!    2004-10-07      added template_pixel_*, same_shape_pixels_*
!    2004-12-15      generic forms for convert_inplace (I,R,D with 1 or N dimensions), 
!                     remove_dipole (R,D)
!    2005-01-25      ensure backward compatibility of remove_dipole
!    2005-08-03      edit same_shape_pixels_nest to comply with GFORTRAN
!==================================================================

  USE healpix_types
  USE misc_utils
  IMPLICIT none

  INTEGER(KIND=i4b), private, PARAMETER :: ns_max=8192 ! 2^13 : largest nside available

  !initialise array x2pix, y2pix and pix2x, pix2y used in several routines
  integer(KIND=i4b), private, save, dimension(128) :: x2pix=0,y2pix=0

  integer(KIND=i4b), private, save, dimension(0:1023) :: pix2x=0, pix2y=0

  ! obsolete
  interface convert_inplace_real
     module procedure convert_inplace_real_1d, convert_inplace_real_nd
  end interface 
  interface convert_inplace_int
     module procedure convert_inplace_int_1d, convert_inplace_int_nd
  end interface 

  ! generic form
  interface convert_inplace
     module procedure convert_inplace_real_1d, convert_inplace_real_nd, &
          &           convert_inplace_int_1d, convert_inplace_int_nd, &
          &           convert_inplace_double_1d, convert_inplace_double_nd
  end interface 

  interface remove_dipole
     module procedure remove_dipole_real, remove_dipole_double, &
          &           remove_dipole_real_old, remove_dipole_double_old, &
          &           remove_dipole_real_v12, remove_dipole_double_v12
  end interface

  interface add_dipole
     module procedure add_dipole_real, add_dipole_double
  end interface

  interface convert_nest2ring
     module procedure convert_nest2ring_int_1d, &
          &           convert_nest2ring_real_1d, &
          &           convert_nest2ring_double_1d, &
          &           convert_nest2ring_int_nd, &
          &           convert_nest2ring_real_nd, &
          &           convert_nest2ring_double_nd
  end interface 

  interface convert_ring2nest
     module procedure convert_ring2nest_int_1d, &
          &           convert_ring2nest_real_1d, &
          &           convert_ring2nest_double_1d, &
          &           convert_ring2nest_int_nd, &
          &           convert_ring2nest_real_nd, &
          &           convert_ring2nest_double_nd
  end interface 

  interface medfiltmap
     module procedure medfiltmap_s, medfiltmap_d
  end interface 

  private

  public :: remove_dipole, add_dipole, &
       & query_strip, &
       & query_polygon, &
       & query_triangle, &
       & query_disc, &
       & pix2ang_ring, pix2vec_ring, ang2pix_ring, vec2pix_ring, &
       & pix2ang_nest, pix2vec_nest, ang2pix_nest, vec2pix_nest, &
       & convert_nest2ring, convert_ring2nest, &
       & convert_inplace, convert_inplace_real, convert_inplace_int, &
       & nest2ring, ring2nest, xy2pix_nest, pix2xy_nest, &
       & mk_pix2xy, mk_xy2pix, &
       & neighbours_nest, &
       & next_in_line_nest, &
       & ang2vec, vec2ang, &
       & npix2nside, nside2npix, &
       & surface_triangle, angdist, vect_prod

  public :: nside2ntemplates, &
       & template_pixel_ring, same_shape_pixels_ring, &
       & template_pixel_nest, same_shape_pixels_nest

  public :: medfiltmap

  public :: intrs_intrv, in_ring, ring_num, ring2z ! arcane usage

  public :: getdisc_ring  ! obsolete


contains

  !=======================================================================
  subroutine query_strip ( nside, theta1, theta2, listpix, nlist, nest, inclusive)
  !=======================================================================
    ! query_strip ( nside, theta1, theta2, listpix, nlist, nest, inclusive)
    !
    ! finds pixels having a colatitude (measured from North Pole):
    !  theta1 < colatitude < theta2
    !  with 0 <= theta1 < theta2 <= Pi
    ! if theta2 < theta1 then pixels with
    ! 0 <= colatitude < theta2   or   theta1 < colatitude < Pi are returned
    !
    ! nside :          I4,       input,  resolution parameter
    ! theta1, theta2 : DP,       input,  bounds
    ! listpix :        I4 array, input,  list of pixels
    ! nlist :          I4,       output, number of pixels in list
    ! nest  :          I4 opt.,  input, =0 : RING scheme, =1 : NESTED scheme
    ! inclusive:       I4 opt.,  input, if=1 include pixel touched by strip, 
    !                                   if 0 keep pixels with center in strip
    !
    ! v1.0, EH, Caltech, Jan-2002
    ! v2.0: 2004-12 : added inclusive keyword
    !=======================================================================
    integer(kind=I4B), intent(in)                    :: nside
    real(kind=DP),     intent(in)                    :: theta1, theta2
    integer(kind=I4B), intent(out), dimension(0:)    :: listpix
    integer(kind=I4B), intent(out)                   :: nlist
    integer(kind=I4B), intent(in),  optional         :: nest
    integer(kind=i4b), intent(in),  optional         :: inclusive

    integer(kind=I4B)                  :: npix, nstrip
    integer(kind=I4B)                  :: iz, ip, is, irmin, irmax
    integer(kind=I4B)                  :: ilist, nir, nlost, list_size
    real(kind=DP)                      :: phi0, dphi, zu, zd, zring
    real(kind=DP), dimension(1:4)      :: colrange
    character(len=*), parameter        :: code = "query_strip"
    integer(kind=I4B), dimension(:), allocatable :: listir
    logical(kind=LGT) :: be_inclusive

    !=======================================================================
    list_size = size(listpix)
    !     ---------- check inputs ----------------
    npix = nside2npix(nside)
    if (npix < 0) then
       print*,code//"> Nside should be a power of 2"
       print*,code//"> current value = ",nside
       call fatal_error("> program abort ")
    endif

    if    (theta1 < 0.0_dp .or. theta1 > PI .or. &
         & theta2 < 0.0_dp .or. theta2 > PI) then
       write(unit=*,fmt="(a)") code//"> the colatitudes are in RADIAN "
       write(unit=*,fmt="(a)") code//"> and should lie in [0,Pi] "
       print*,code//"> current value = ", theta1, theta2
       call fatal_error("> program abort ")
    endif

    if (theta1 <= theta2) then
       nstrip = 1
       colrange(1:2*nstrip) = (/ theta1, theta2 /)
    else
       nstrip = 2
       colrange(1:2*nstrip) = (/ 0.0_dp, theta2, theta1, PI/)
    endif

    be_inclusive = .false.
    if (present(inclusive)) be_inclusive = (inclusive==1)

    !     ------------- loop on strips ---------------------
    ilist = -1
    allocate(listir(0:4*nside-1))
    do is =0, nstrip-1
       zu = cos(colrange(2*is+1)) ! upper bound in z
       zd = cos(colrange(2*is+2)) ! lower bound in z
       irmin = ring_num(nside, zu)
       irmax = ring_num(nside, zd)


       !     ------------- loop on ring number ---------------------
       do iz = irmin, irmax

          zring = ring2z(nside, iz) ! z of ring being considered

          if ((zring >= zd .and. zring <= zu) .or. be_inclusive) then
             !        ------- finds pixels in the ring ---------
             phi0 = 0
             dphi = PI
             call in_ring(nside, iz, phi0, dphi, listir, nir, nest)
             
             nlost = ilist + nir + 1 - list_size
             if ( nlost > 0 ) then
                print*,code//"> listpix is too short, it will be truncated at ",nir
                print*,"                         pixels lost : ", nlost
                nir = nir - nlost
             endif
             do ip = 0, nir-1
                ilist = ilist + 1
                listpix(ilist) = listir(ip)
             enddo
          endif

       enddo

    enddo

    !     ------ total number of pixel in the strip --------
    nlist = ilist + 1

    !     ------- deallocate memory and exit ------
    DEALLOCATE(listir)

    return
  end subroutine query_strip
  !=======================================================================
  subroutine query_polygon ( nside, vlist, nv, listpix, nlist, nest, inclusive)
    !=======================================================================
    ! finds pixels that lie within a CONVEX polygon defined by its vertex on the sphere
    !
    ! nside             : IN
    ! vlist(1:3, 0:n-1) : IN, list of vertices
    ! nv                : IN, number of vertices to be used (nv <= n)
    ! listpix           : OUT
    ! nlist             : OUT
    ! nest              : IN, OPTIONAL
    ! inclusive         : IN, OPTIONAL
    !
    ! algorithm:
    !   the polygon is divided into triangles
    !   vertex 0 belongs to all the triangles
    !
    ! v1.0, EH, Caltech, Dec-2001
    !=======================================================================
    USE num_rec, ONLY : isort
    integer(kind=I4B), intent(in)                    :: nside
    real(kind=DP),     intent(in),  dimension(1:,1:) :: vlist
    integer(kind=I4B), intent(in)                    :: nv
    integer(kind=I4B), intent(out), dimension(0:)    :: listpix
    integer(kind=I4B), intent(out)                   :: nlist
    integer(kind=I4B), intent(in),  optional         :: nest
    integer(kind=i4b), intent(in),  optional         :: inclusive

    real(kind=DP),     dimension(:), pointer         :: vp0, vp1, vp2
    real(kind=DP),     dimension(1:3)                :: vo
    real(kind=DP),     dimension(:,:), allocatable, target   :: vvlist
    real(kind=DP),     dimension(:),   allocatable   :: ss
    real(kind=DP)                                    :: surface, fsky, hand
    integer(kind=I4B)                                :: npix, n_in_trg, ilist, ntl
    integer(kind=I4B)                                :: n_remain
    integer(kind=I4B)                                :: i0, i1, i2, i, k
    integer(kind=I4B)                                :: np, nm, nlow
    integer(kind=I4B)                                :: list_size, nlost
    integer(kind=I4B), dimension(1:1)                :: ix
    integer(kind=I4B), dimension(:), allocatable     :: templist
    character(len=*), parameter :: code = "QUERY_POLYGON"

    !=======================================================================

    list_size = size(listpix)
    n_remain = nv
    ilist = -1
    listpix(:) = -1

    call assert(n_remain>=3, &
      code//"> The polygon should have at least 3 vertices")

    if (size(vlist) < n_remain*3) then
       print*,code//"> ",size(vlist)/3," vertices are given"
       print*,code//"> expected : ",n_remain
       call fatal_error
    endif

    allocate(vvlist(1:3, 1:n_remain))
    vvlist = vlist

    !-----------------------------------------------------------------
    ! check that the polygon is convex or has only one concave vertex
    !-----------------------------------------------------------------
    if (n_remain == 3) goto 1000 ! a triangle is always convex

    allocate(ss(1:n_remain))
    do i1 = 0, n_remain-1
       i0 = modulo(i1-1, n_remain)  ! in [0, n_remain-1]
       i2 = modulo(i1+1, n_remain)  ! in [0, n_remain-1]

       vp0 => vvlist(:,i0+1)
       vp1 => vvlist(:,i1+1)
       vp2 => vvlist(:,i2+1)

       ! computes the handedness   (v0 x v2) . v1  for each vertex v1
       call vect_prod(vp0, vp2, vo)
       hand = dot_product(vo, vp1)
       ss(i1+1) = sign(1.0_dp, hand) ! either +1 or -1
    enddo
    np = count( ss > 0.0_dp) ! number of vertex with positive handedness
    nm = n_remain - np

    nlow = min(np,nm)
    if (nlow == 0) goto 1000
    if (nlow == 1) then
       ! only one concave vertex
       if (np == 1) then
          ix = maxloc(ss)
       else
          ix = minloc(ss)
       endif
       ! rotate pixel list to put that vertex in #0
       vvlist = cshift(vvlist, ix(1)-1, dim = 2)
    endif
    if (nlow > 1) then
       ! more than 1 concave vertex
       print*,"***************************************"
       print*,"The polygon is not convex"
       print*," and has more than one concave vertex"
       print*,"The result is unpredictable"
       print*,"***************************************"
    endif
    deallocate(ss)

1000 continue

    !--------------------------------------------
    ! fill the polygon, one triangle at a time
    !--------------------------------------------
    npix = nside2npix(nside)
    do
       ! build triangle from vertices #0, n-2, and n-1
       vp0 => vvlist(1:3,1)
       vp1 => vvlist(1:3,n_remain-1)
       vp2 => vvlist(1:3,n_remain  )

       ! computes its surface
       call surface_triangle( vp0, vp1, vp2, surface )
       fsky = surface/FOURPI
       n_in_trg = npix * (fsky * 1.4) + 12*nside
       n_in_trg = min(n_in_trg, npix)

       ! find pixels within triangle
       allocate(templist(0:n_in_trg-1))
       call query_triangle( nside, vp0, vp1, vp2, templist, ntl, nest=nest,inclusive=inclusive )

       ! merge new list with existing one
       nlost = ilist + ntl + 1 - list_size
       if ( nlost > 0 ) then
          print*,code//"> listpix is too short, it will be truncated at ",list_size
          print*,"                         pixels lost : ", nlost
          print*, list_size
          ntl = ntl - nlost
       endif
       do i = 0, ntl - 1
          ilist = ilist + 1
          listpix(ilist) = templist(i)
       enddo
       deallocate(templist)

       if (n_remain == 3) exit

       ! prune vertex list
       n_remain = n_remain - 1

    enddo
    deallocate(vvlist)

    !-------------------------
    ! make final pixel list
    !-------------------------
    nlist = ilist + 1

    ! sort final list
    call isort(nlist, listpix)

    ! remove redondant pixels
    ! (we keep 0th element of the list)
    k = 1
    do i = 1, nlist-1
       if (listpix(i) > listpix(i-1)) then
          listpix(k) = listpix(i)
          k = k + 1
       endif
    enddo
    nlist = k

    ! pad end of list with -1
    listpix(nlist:list_size-1) = -1

    return
  end subroutine query_polygon
  !=======================================================================
  subroutine query_triangle ( nside, v1, v2, v3, listpix, nlist, nest, inclusive)
    !=======================================================================
    !
    !    query_triangle ( nside, v1, v2, v3, listpix, nlist, nest, inclusive)
    !    --------------
    !     nside       = resolution parameter (a power of 2)
    !     v1, v2, v3  = 3D vector location of the 3 triangle vertices
    !     list_pix    = list of pixel lying in the triangle
    !     nlist       = number of pixels in the list
    !     nest  (OPT), :0 by default, the output list is in RING scheme
    !                  if set to 1, the output list is in NESTED scheme
    !     inclusive (OPT) , :0 by default, only the pixels whose center
    !                       lie in the triangle are listed on output
    !                  if set to 1, all pixels overlapping the triangle are output
    !
    !     NB : the dimension of the listpix array is fixed in the calling
    !     routine and should be large enough for the specific configuration
    !
    !
    ! v1.0, EH, Caltech, Nov-Dec-2001
    ! v1.1,  Aug-Sep-2002 : added nest and inclusive
    !=======================================================================
    integer(kind=I4B), intent(in) :: nside
    real(kind=DP),     intent(in),  dimension(1:)  :: v1, v2, v3
    integer(kind=I4B), intent(out), dimension(0:)  :: listpix
    integer(kind=I4B), intent(out)                 :: nlist
    integer(kind=I4B), intent(in), optional        :: nest
    integer(kind=I4B), intent(in), optional        :: inclusive

    integer(kind=I4B) :: npix, ilist !, ip1, ip2, ip3
    integer(kind=I4B) :: iz, irmin, irmax
    integer(kind=I4B) :: n12, n123a, n123b, ndom
    integer(kind=I4B), dimension(:),   allocatable  :: listir
    logical(kind=LGT) :: test1, test2, test3
    logical(kind=LGT) :: test1a, test1b, test2a, test2b, test3a, test3b
    real(kind=DP) :: dth1, dth2, determ, sdet
    real(kind=DP) :: zmax, zmin, z1max, z1min, z2max, z2min, z3max, z3min
    real(kind=DP) :: z, zz
    real(kind=DP) :: tgth, st
    real(kind=DP) :: offset, sin_off
    real(kind=DP), dimension(1:3,1:3) :: vv, vo
    real(kind=DP), dimension(1:3) :: sprod, sto, phi0i, tgthi
    real(kind=DP), dimension(1:3) :: dc
    real(kind=DP), dimension(1:2,1:3) :: dom
    real(kind=DP), dimension(1:4) :: dom12, dom123a, dom123b
    real(kind=DP), dimension(1:6) :: alldom
    real(kind=DP) :: a_i, b_i, phi0, dphiring
    integer(kind=I4B) :: idom, nir, ip
    integer(kind=I4B) :: status
    integer(kind=I4B) :: j
    character(len=*), parameter :: code = "QUERY_TRIANGLE"
    integer(kind=I4B) :: list_size, nlost
    logical(LGT)      :: do_inclusive

    !=======================================================================

    list_size = size(listpix)

    npix = nside2npix(nside)
    if (npix < 0) then
       print*,"Invalid Nside = ",nside
       call fatal_error
    endif

    do_inclusive = .false.
    if (present(inclusive)) then
       if (inclusive == 1) do_inclusive = .true.
    endif

    !   ! sort pixels by number
    !   ip1 = MIN(ipix1, ipix2, ipix3)
    !   ip3 = MAX(ipix1, ipix2, ipix3)
    !   ip2 = ipix1 + ipix2 + ipix3 - ip1 - ip3

    !   !     ---------- check inputs ----------------
    !   if (ip1 < 0 .or. ip3 > npix-1) then
    !      write(unit=*,fmt="(a)") " > Non valid choice for pixel number :"
    !      write(unit=*,fmt="(a)") " > ",ipix1,ipix2,ipix3
    !      write(unit=*,fmt="(a,i2,i10)") " > valid range : ",0,npix-1
    !      nlist = 0
    !      listpix(0) = -1
    !   endif

    !   call pix2vec_ring( nside, ip1, vv(1:3,1))
    !   call pix2vec_ring( nside, ip2, vv(1:3,2))
    !   call pix2vec_ring( nside, ip3, vv(1:3,3))

    vv(1:3,1) = v1(1:3) / sqrt(dot_product(v1,v1))
    vv(1:3,2) = v2(1:3) / sqrt(dot_product(v2,v2))
    vv(1:3,3) = v3(1:3) / sqrt(dot_product(v3,v3))

    !     --------- allocate memory -------------
    ALLOCATE( listir(0: 4*nside-1), STAT = status)
    if (status /= 0) then
       write(unit=*,fmt="(a)") code//" > can not allocate memory for listir :"
       call fatal_error(" > program abort ")
    endif

    dth1 = 1.0_dp / (3.0_dp*REAL(nside,kind=dp)**2)
    dth2 = 2.0_dp / (3.0_dp*REAL(nside,kind=dp))


    ! determ = (vect1 X vect2) . vect3
    ! determines the left(<0)/right(>0) handedness of the triangle
    determ = vv(1,1)*vv(2,2)*vv(3,3) + vv(1,2)*vv(2,3)*vv(3,1) + vv(1,3)*vv(2,1)*vv(3,2) &
         & - vv(3,1)*vv(2,2)*vv(1,3) - vv(3,2)*vv(2,3)*vv(1,1) - vv(3,3)*vv(2,1)*vv(1,2)

    if (abs(determ) < 1.e-20_dp) then
       print*,' ************************************************************'
       print*,' The triangle is degenerated (2 of the vertices are antipodal)'
       print*,' The query can not be performed '
       print*,' ************************************************************'
       call fatal_error
    endif

    !   print*,determ
    sdet = SIGN(1.0_dp,determ) ! = +1 or -1, the sign of determ

    ! scalar product of vertices vectors
    sprod(1) = dot_product(vv(1:3,2),vv(1:3,3))
    sprod(2) = dot_product(vv(1:3,3),vv(1:3,1))
    sprod(3) = dot_product(vv(1:3,1),vv(1:3,2))

    ! vector orthogonal to the great circle containing the vertex doublet
    call vect_prod(vv(1:3,2), vv(1:3,3), vo(1:3,1))
    call vect_prod(vv(1:3,3), vv(1:3,1), vo(1:3,2))
    call vect_prod(vv(1:3,1), vv(1:3,2), vo(1:3,3))

    ! normalize the orthogonal vector
    vo(1:3,1) = vo(1:3,1) /  SQRT(SUM(vo(1:3,1)**2))
    vo(1:3,2) = vo(1:3,2) /  SQRT(SUM(vo(1:3,2)**2))
    vo(1:3,3) = vo(1:3,3) /  SQRT(SUM(vo(1:3,3)**2))

    ! test presence of poles in the triangle
    zmax = -1.0_dp
    zmin =  1.0_dp
    test1 = (vo(3,1) * sdet >= 0.0_dp) ! north pole in hemisphere defined by 2-3
    test2 = (vo(3,2) * sdet >= 0.0_dp) ! north pole in hemisphere defined by 1-2
    test3 = (vo(3,3) * sdet >= 0.0_dp) ! north pole in hemisphere defined by 1-3
    if (test1 .and. test2 .and. test3) then
       zmax = 1.0_dp ! north pole in the triangle
    endif
    if ((.not.test1) .and. (.not.test2) .and. (.not.test3)) then
       zmin = -1.0_dp ! south pole in the triangle
    endif

    ! look for northernest and southernest points in the triangle
    !   node(1,2) = vector of norm=1, in the plane defined by (1,2) and with z=0
    test1a = ((vv(3,3) - sprod(1) * vv(3,2)) >= 0.0_dp) ! segment 2-3 : -vector(3) . node(2,3)
    test1b = ((vv(3,2) - sprod(1) * vv(3,3)) >= 0.0_dp) !                vector(2) . node(2,3)
    test2a = ((vv(3,3) - sprod(2) * vv(3,1)) >= 0.0_dp) ! segment 1-3 : -vector(3) . node(1,3)
    test2b = ((vv(3,1) - sprod(2) * vv(3,3)) >= 0.0_dp) !                vector(1) . node(1,3)
    test3a = ((vv(3,2) - sprod(3) * vv(3,1)) >= 0.0_dp) ! segment 1-2 : -vector(2) . node(1,2)
    test3b = ((vv(3,1) - sprod(3) * vv(3,2)) >= 0.0_dp) !                vector(1) . node(1,2)

    ! sin of theta for orthogonal vector
    sto(1:3) = SQRT( (1.0_dp-vo(3,1:3))*(1.0_dp+vo(3,1:3)) )

    ! for each segment (=side of the triangle) the extrema are either
    ! - the 2 vertices
    ! - one of the vertices and a point within the segment

    ! segment 2-3
    z1max = vv(3,2)
    z1min = vv(3,3)
    if ( test1a .EQV. test1b ) then
       zz = sto(1)
       if ((vv(3,2)+vv(3,3)) >= 0.0_dp) then
          z1max =  zz
       else
          z1min = -zz
       endif
    endif

    ! segment 1-3
!     z2max = vv(3,1)
!     z2min = vv(3,3)
    z2max = vv(3,3)
    z2min = vv(3,1)
    if ( test2a .EQV. test2b ) then
       zz = sto(2)
       if ((vv(3,1)+vv(3,3)) >= 0.0_dp) then
          z2max =  zz
       else
          z2min = -zz
       endif
    endif

    ! segment 1-2
    z3max = vv(3,1)
    z3min = vv(3,2)
    if ( test3a .EQV. test3b ) then
       zz = sto(3)
       if ((vv(3,1)+vv(3,2)) >= 0.0_dp) then
          z3max =  zz
       else
          z3min = -zz
       endif
    endif

    zmax = MAX(z1max, z2max, z3max, zmax)
    zmin = MIN(z1min, z2min, z3min, zmin)

    ! if we are inclusive, move the upper point up, and the lower point down, by a half pixel size
    offset = 0.0_dp
    sin_off = 0.0_dp
    if (do_inclusive) then
       offset = PI / (4.0_dp*nside) ! half pixel size
       sin_off = sin(offset)
       zmax = min( 1.0_dp, cos( acos(zmax) - offset) )
       zmin = max(-1.0_dp, cos( acos(zmin) + offset) )
    endif

    !   print*,"zmin, zmax ",zmin,zmax

    ! northernest and sourthernest ring number
    irmin = ring_num(nside, zmax)
!!!  irmin = MAX(1, irmin - 1) ! start from a higher point, to be safe
    irmax = ring_num(nside, zmin)
!!!  irmax = MIN(4*nside-1, irmax + 1) ! go down to a lower point

    ilist = -1
    !    print*,"irmin, irmax ",irmin,irmax

    ! -------- loop on the rings -------------------------

    tgthi(1:3) = -1.0e30_dp * vo(3,1:3)
    phi0i(1:3) =  0.0_dp
    do j=1,3
       if (sto(j) > 1.0e-10_dp) then
          tgthi(j) = -vo(3,j) / sto(j)  ! -cotan(theta_orth)
          phi0i(j) = ATAN2(vo(2,j),vo(1,j))
       endif
    enddo
    !   print*,tgthi,phi0i

    ! the triangle boundaries are geodesics : intersection of the sphere with plans going thru (0,0,0)
    ! if we are inclusive, the boundaries are the intersecion of the sphere with plans pushed outward
    ! by sin(offset)
    do iz = irmin, irmax
       if (iz <= nside-1) then      ! north polar cap
          z = 1.0_dp  - REAL(iz,kind=dp)**2 * dth1
       else if (iz <= 3*nside) then    ! tropical band + equat.
          z = REAL(2*nside-iz,kind=dp) * dth2
       else
          z = - 1.0_dp + REAL(4*nside-iz,kind=dp)**2 * dth1
       endif
       ! computes the 3 intervals described by the 3 great circles
       st = SQRT((1.0_dp - z)*(1.0_dp + z))
       tgth = z / st ! cotan(theta_ring)
!        dc(1:3)  = tgthi(1:3) * tgth - sdet * sin_off / (sto(1:3) * st)
       dc(1:3)  = tgthi(1:3) * tgth - sdet * sin_off / ((sto(1:3)+1.e-30_dp) * st) ! sto is slightly offset to avoid division by 0

       do j=1,3
          if (dc(j)*sdet <= -1.0_dp) then  ! the whole iso-latitude ring is on the right side of the great circle
             dom(1:2, j) = (/ 0.0_dp, twopi /)
          else if (dc(j)*sdet >= 1.0_dp) then ! all on the wrong side
             dom(1:2, j) = (/ -1.000001_dp, -1.0_dp /) * j
          else ! some is good, some is bad
             dom(1:2, j) = MODULO( phi0i(j) + (ACOS(dc(j)) * sdet) * (/-1.0_dp, 1.0_dp /), twopi)
          endif
       enddo

       ! identify the intersections (0,1,2 or 3) of the 3 intervals
       call intrs_intrv( dom(1:2,1), dom(1:2,2), dom12, n12)
       if (n12 == 0) goto 20
       if (n12 == 1) then
          call intrs_intrv( dom(1:2,3), dom12, dom123a, n123a)
          if (n123a == 0) goto 20
          alldom(1:2*n123a) = dom123a(1:2*n123a)
          ndom = n123a ! 1 or 2
       endif
       if (n12 == 2) then
          call intrs_intrv( dom(1:2,3), dom12(1:2), dom123a, n123a)
          call intrs_intrv( dom(1:2,3), dom12(3:4), dom123b, n123b)
          ndom = n123a + n123b ! 0, 1, 2 or 3
          if (ndom == 0) goto 20
          if (n123a /= 0) alldom(1:2*n123a) = dom123a(1:2*n123a)
          if (n123b /= 0) alldom(2*n123a+1:2*ndom)  = dom123b(1:2*n123b)
          if (ndom > 3) then
             print*,code//"> too many intervals found"
          endif
       endif
       do idom=0,ndom-1
          a_i = alldom(2*idom+1)
          b_i = alldom(2*idom+2)
          phi0 = (a_i + b_i) * 0.5_dp
          dphiring = (b_i - a_i) * 0.5_dp
          if (dphiring < 0.0_dp) then
             phi0 = phi0 + pi
             dphiring = dphiring + pi
          endif

          !        ------- finds pixels in the triangle on that ring ---------
          call in_ring(nside, iz, phi0, dphiring, listir, nir, nest=nest)

          !        ----------- merge pixel lists -----------
          nlost = ilist + nir + 1 - list_size
          if ( nlost > 0 ) then
             print*,code//"> listpix is too short, it will be truncated at ",nir
             print*,"                         pixels lost : ", nlost
             print*, list_size
             nir = nir - nlost
          endif
          do ip = 0, nir-1
             ilist = ilist + 1
             listpix(ilist) = listir(ip)
          enddo
       enddo
20     continue
    enddo !-----------------------------------------

    !     ------ total number of pixel in the disc --------
    nlist = ilist + 1

    !     ------- deallocate memory and exit ------
    DEALLOCATE(listir)

    return
  end subroutine query_triangle

  !=======================================================================
  subroutine getdisc_ring ( nside, vector0, radius, listpix, nlist)
    integer(kind=I4B), intent(in)                 :: nside
    real(kind=DP),     intent(in), dimension(1:)  :: vector0
    real(kind=DP),     intent(in)                 :: radius
    integer(kind=I4B), intent(out), dimension(0:) :: listpix
    integer(kind=I4B), intent(out)                :: nlist
  !=======================================================================

    print*,"-------------------------------------------------------------"
    print*,"WARNING : the routine getdisc_ring is now obsolete"
    print*,"  Use "
    print*," call query_disc(nside, vector0, radius_radian, listpix, nlist [, nest]) "
    print*,"  instead (same module)"
    print*," It lets you choose the indexing scheme (RING or NESTED) "
    print*," used for the outuput"
    print*,"-------------------------------------------------------------"

    call query_disc(nside, vector0, radius, listpix, nlist)

    return
  end subroutine getdisc_ring

  !=======================================================================
  subroutine query_disc ( nside, vector0, radius, listpix, nlist, nest, inclusive)
    !=======================================================================
    !
    !      query_disc (Nside, Vector0, Radius, Listpix, Nlist[, Nest, Inclusive])
    !      ----------
    !      routine for pixel query in the RING or NESTED scheme
    !      all pixels within an angular distance Radius of the center
    !
    !     Nside    = resolution parameter (a power of 2)
    !     Vector0  = central point vector position (x,y,z in double precision)
    !     Radius   = angular radius in RADIAN (in double precision)
    !     Listpix  = list of pixel closer to the center (angular distance) than Radius
    !     Nlist    = number of pixels in the list
    !     nest  (OPT), :0 by default, the output list is in RING scheme
    !                  if set to 1, the output list is in NESTED scheme
    !     inclusive (OPT) , :0 by default, only the pixels whose center
    !                       lie in the triangle are listed on output
    !                  if set to 1, all pixels overlapping the triangle are output
    !
    !      * all pixel numbers are in {0, 12*Nside*Nside - 1}
    !     NB : the dimension of the listpix array is fixed in the calling
    !     routine and should be large enough for the specific configuration
    !
    !      lower level subroutines called by getdisc_ring :
    !       (you don't need to know them)
    !      ring_num (nside, ir)
    !      --------
    !      in_ring(nside, iz, phi0, dphi, listir, nir, nest=nest)
    !      -------
    !
    ! v1.0, EH, TAC, ??
    ! v1.1, EH, Caltech, Dec-2001
    !=======================================================================
    integer(kind=I4B), intent(in)                 :: nside
    real(kind=DP),     intent(in), dimension(1:)  :: vector0
    real(kind=DP),     intent(in)                 :: radius
    integer(kind=I4B), intent(out), dimension(0:) :: listpix
    integer(kind=I4B), intent(out)                :: nlist
    integer(kind=I4B), intent(in), optional       :: nest
    integer(kind=I4B), intent(in), optional       :: inclusive

    INTEGER(KIND=I4B) :: irmin, irmax, ilist, iz, ip, nir, npix
    REAL(KIND=DP) :: norm_vect0
    REAL(KIND=DP) :: x0, y0, z0, radius_eff, fudge
    REAL(KIND=DP) :: a, b, c, cosang
    REAL(KIND=DP) :: dth1, dth2
    REAL(KIND=DP) :: phi0, cosphi0, cosdphi, dphi
    REAL(KIND=DP) :: rlat0, rlat1, rlat2, zmin, zmax, z
    INTEGER(KIND=I4B), DIMENSION(:),   ALLOCATABLE  :: listir
    INTEGER(KIND=I4B) :: status
    character(len=*), parameter :: code = "QUERY_DISC"
    integer(kind=I4B) :: list_size, nlost
    logical(LGT) :: do_inclusive

    !=======================================================================

    list_size = size(listpix)
    !     ---------- check inputs ----------------
    npix = 12 * nside * nside

    if (radius < 0.0_dp .or. radius > PI) then
       write(unit=*,fmt="(a)") code//"> the angular radius is in RADIAN "
       write(unit=*,fmt="(a)") code//"> and should lie in [0,Pi] "
       call fatal_error("> program abort ")
    endif

    do_inclusive = .false.
    if (present(inclusive)) then
       if (inclusive == 1) do_inclusive = .true.
    endif

    !     --------- allocate memory -------------
    ALLOCATE( listir(0: 4*nside-1), STAT = status)
    if (status /= 0) then
       write(unit=*,fmt="(a)") code//"> can not allocate memory for listir :"
       call fatal_error("> program abort ")
    endif

    dth1 = 1.0_dp / (3.0_dp*real(nside,kind=dp)**2)
    dth2 = 2.0_dp / (3.0_dp*real(nside,kind=dp))

    radius_eff = radius
    if (do_inclusive) then
!        fudge = PI / (4.0_dp*nside) ! increase radius by half pixel size
       fudge = acos(TWOTHIRD) / real(nside,kind=dp) ! 1.071* half pixel size
       radius_eff = radius + fudge
    endif
    cosang = COS(radius_eff)

    !     ---------- circle center -------------
    norm_vect0 =  SQRT(DOT_PRODUCT(vector0,vector0))
    x0 = vector0(1) / norm_vect0
    y0 = vector0(2) / norm_vect0
    z0 = vector0(3) / norm_vect0

    phi0=0.0_dp
    if ((x0/=0.0_dp).or.(y0/=0.0_dp)) phi0 = ATAN2 (y0, x0)  ! in ]-Pi, Pi]
    cosphi0 = COS(phi0)
    a = x0*x0 + y0*y0

    !     --- coordinate z of highest and lowest points in the disc ---
    rlat0  = ASIN(z0)    ! latitude in RAD of the center
    rlat1  = rlat0 + radius_eff
    rlat2  = rlat0 - radius_eff
    if (rlat1 >=  halfpi) then
       zmax =  1.0_dp
    else
       zmax = SIN(rlat1)
    endif
    irmin = ring_num(nside, zmax)
    irmin = MAX(1, irmin - 1) ! start from a higher point, to be safe

    if (rlat2 <= -halfpi) then
       zmin = -1.0_dp
    else
       zmin = SIN(rlat2)
    endif
    irmax = ring_num(nside, zmin)
    irmax = MIN(4*nside-1, irmax + 1) ! go down to a lower point

    ilist = -1

    !     ------------- loop on ring number ---------------------
    do iz = irmin, irmax

       if (iz <= nside-1) then      ! north polar cap
          z = 1.0_dp  - real(iz,kind=dp)**2 * dth1
       else if (iz <= 3*nside) then    ! tropical band + equat.
          z = real(2*nside-iz,kind=dp) * dth2
       else
          z = - 1.0_dp + real(4*nside-iz,kind=dp)**2 * dth1
       endif

       !        --------- phi range in the disc for each z ---------
       b = cosang - z*z0
       c = 1.0_dp - z*z
       if ((x0==0.0_dp).and.(y0==0.0_dp)) then
          cosdphi=-1.0_dp
          dphi=PI
          goto 500
       endif
       cosdphi = b / SQRT(a*c)
       if (ABS(cosdphi) <= 1.0_dp) then
          dphi = ACOS (cosdphi) ! in [0,Pi]
       else
          if (cosphi0 < cosdphi) goto 1000 ! out of the disc
          dphi = PI ! all the pixels at this elevation are in the disc
       endif
500    continue

       !        ------- finds pixels in the disc ---------
       call in_ring(nside, iz, phi0, dphi, listir, nir, nest)

       !        ----------- merge pixel lists -----------
       nlost = ilist + nir + 1 - list_size
       if ( nlost > 0 ) then
          print*,code//"> listpix is too short, it will be truncated at ",nir
          print*,"                         pixels lost : ", nlost
          nir = nir - nlost
       endif
       do ip = 0, nir-1
          ilist = ilist + 1
          listpix(ilist) = listir(ip)
       enddo

1000   continue
    enddo

    !     ------ total number of pixel in the disc --------
    nlist = ilist + 1


    !     ------- deallocate memory and exit ------
    DEALLOCATE(listir)

    return
  end subroutine query_disc
  !=======================================================================
  function ring_num (nside, z) result(ring_num_result)
    !=======================================================================
    !     returns the ring number in {1, 4*nside-1}
    !     from the z coordinate
    !=======================================================================
    INTEGER(KIND=I4B) :: ring_num_result
    REAL(KIND=DP), INTENT(IN) :: z
    INTEGER(KIND=I4B), INTENT(IN) :: nside

    INTEGER(KIND=I4B) :: iring
    !=======================================================================

    !     ----- equatorial regime ---------
    iring = NINT( nside*(2.0_dp-1.500_dp*z))

    !     ----- north cap ------
    if (z > twothird) then
       iring = NINT( nside* SQRT(3.0_dp*(1.0_dp-z)))
       if (iring == 0) iring = 1
    endif

    !     ----- south cap -----
    if (z < -twothird   ) then
       iring = NINT( nside* SQRT(3.0_dp*(1.0_dp+z)))
       if (iring == 0) iring = 1
       iring = 4*nside - iring
    endif

    ring_num_result = iring

    return
  end function ring_num
  !=======================================================================
  function ring2z (nside, ir) result(z)
    !=======================================================================
    !     returns the z coordinate of ring ir for Nside
    !=======================================================================
    integer(kind=I4B), intent(in) :: ir
    integer(kind=I4B), intent(in) :: nside
    real(kind=DP)                 :: z

    real(DP)     :: dth1, dth2, fn
    !=======================================================================

    fn = real(nside,kind=DP)
    dth1 = 1.0_dp / (3.0_dp * fn * fn)
    dth2 = 2.0_dp / (3.0_dp * fn)

    if (ir < nside) then  ! polar cap (north)
       z = 1.0_dp - ir**2 * dth1
    else if (ir < 3*nside) then ! tropical band
       z = real( 2*nside-ir, kind=DP) * dth2
    else                  ! polar cap (south)
       z = -1.0_dp + (4*nside-ir)**2 * dth1
    endif
    
    return
  end function ring2z
  !=======================================================================
  subroutine in_ring (nside, iz, phi0, dphi, listir, nir, nest)
    !=======================================================================
    !     returns the list of pixels in RING or NESTED scheme (listir)
    !     and their number (nir)
    !     with latitude in [phi0-dphi, phi0+dphi] on the ring ir
    !     (in {1,4*nside-1})
    !     the pixel id-numbers are in {0,12*nside^2-1}
    !     the indexing is RING, unless NEST is set to 1
    !=======================================================================
    integer(kind=i4b), intent(in)                 :: nside, iz
    integer(kind=i4b), intent(out)                :: nir
    real(kind=dp),     intent(in)                 :: phi0, dphi
    integer(kind=i4b), intent(out), dimension(0:) :: listir
    integer(kind=i4b), intent(in), optional       :: nest

!     logical(kind=lgt) :: conservative = .true.
    logical(kind=lgt) :: conservative = .false.
    logical(kind=lgt) :: take_all, to_top, do_ring

    integer(kind=i4b) :: ip_low, ip_hi, i, in, inext, diff
    integer(kind=i4b) :: npix, nr, nir1, nir2, ir, ipix1, ipix2, kshift, ncap
    real(kind=dp)     :: phi_low, phi_hi, shift
    !=======================================================================

    take_all = .false.
    to_top   = .false.
    do_ring  = .true.
    if (present(nest)) then
       do_ring = (nest == 0)
    endif
    npix = 12 * nside * nside
    ncap  = 2*nside*(nside-1) ! number of pixels in the north polar cap
    listir = -1
    nir = 0

    phi_low = MODULO(phi0 - dphi, twopi)
    phi_hi  = MODULO(phi0 + dphi, twopi)
    if (ABS(dphi-PI) < 1.0e-6_dp) take_all = .true.

    !     ------------ identifies ring number --------------
    if (iz >= nside .and. iz <= 3*nside) then ! equatorial region
       ir = iz - nside + 1  ! in {1, 2*nside + 1}
       ipix1 = ncap + 4*nside*(ir-1) !  lowest pixel number in the ring
       ipix2 = ipix1 + 4*nside - 1   ! highest pixel number in the ring
       kshift = MODULO(ir,2)
       nr = nside*4
    else
       if (iz < nside) then       !    north pole
          ir = iz
          ipix1 = 2*ir*(ir-1)        !  lowest pixel number in the ring
          ipix2 = ipix1 + 4*ir - 1   ! highest pixel number in the ring
       else                          !    south pole
          ir = 4*nside - iz
          ipix1 = npix - 2*ir*(ir+1) !  lowest pixel number in the ring
          ipix2 = ipix1 + 4*ir - 1   ! highest pixel number in the ring
       endif
       nr = ir*4
       kshift = 1
    endif

    !     ----------- constructs the pixel list --------------
    if (take_all) then
       nir    = ipix2 - ipix1 + 1
       if (do_ring) then
          listir(0:nir-1) = (/ (i, i=ipix1,ipix2) /)
       else
          call ring2nest(nside, ipix1, in)
          listir(0) = in
          do i=1,nir-1
             call next_in_line_nest(nside, in, inext)
             in = inext
             listir(i) = in
          enddo
       endif
       return
    endif

    shift = kshift * 0.5_dp
    if (conservative) then
       ! conservative : include every intersected pixels,
       ! even if pixel CENTER is not in the range [phi_low, phi_hi]
       ip_low = nint (nr * phi_low / TWOPI - shift)
       ip_hi  = nint (nr * phi_hi  / TWOPI - shift)
       ip_low = modulo (ip_low, nr) ! in {0,nr-1}
       ip_hi  = modulo (ip_hi , nr) ! in {0,nr-1}
    else
       ! strict : include only pixels whose CENTER is in [phi_low, phi_hi]
       ip_low = ceiling (nr * phi_low / TWOPI - shift)
       ip_hi  = floor   (nr * phi_hi  / TWOPI - shift)
!        if ((ip_low - ip_hi == 1) .and. (dphi*nr < PI)) then ! EH, 2004-06-01
       diff = modulo(ip_low - ip_hi, nr) ! in {-nr+1, nr-1} or {0,nr-1} ???
       if (diff < 0) diff = diff + nr    ! in {0,nr-1}
       if ((diff == 1) .and. (dphi*nr < PI)) then
          ! the interval is so small (and away from pixel center)
          ! that no pixel is included in it
          nir = 0
          return
       endif
!        ip_low = min(ip_low, nr-1) !  EH, 2004-05-28
!        ip_hi  = max(ip_hi , 0   )
       if (ip_low >= nr) ip_low = ip_low - nr
       if (ip_hi  <  0 ) ip_hi  = ip_hi  + nr
    endif
    !
    if (ip_low > ip_hi) to_top = .true.
    ip_low = ip_low + ipix1
    ip_hi  = ip_hi  + ipix1

    if (to_top) then
       nir1 = ipix2 - ip_low + 1
       nir2 = ip_hi - ipix1  + 1
       nir  = nir1 + nir2
       if (do_ring) then
          listir(0:nir1-1)   = (/ (i, i=ip_low, ipix2) /)
          listir(nir1:nir-1) = (/ (i, i=ipix1, ip_hi) /)
       else
          call ring2nest(nside, ip_low, in)
          listir(0) = in
          do i=1,nir-1
             call next_in_line_nest(nside, in, inext)
             in = inext
             listir(i) = in
          enddo
       endif
    else
       nir = ip_hi - ip_low + 1
       if (do_ring) then
          listir(0:nir-1) = (/ (i, i=ip_low, ip_hi) /)
       else
          call ring2nest(nside, ip_low, in)
          listir(0) = in
          do i=1,nir-1
             call next_in_line_nest(nside, in, inext)
             in = inext
             listir(i) = in
          enddo
       endif
    endif

    return
  end subroutine in_ring
  !=======================================================================
  subroutine intrs_intrv( d1, d2, di, ni)
    !=======================================================================
    ! computes the intersection di
    ! of 2 intervals d1 (= [a1,b1]) and d2 (= [a2,b2])
    ! on the periodic domain ( = [A,B], where A and B are arbitrary)
    ! ni is the resulting number of intervals (0,1, or 2)
    !
    ! if a1<b1 then d1 = {x | a1 <= x <= b1}
    ! if a1>b1 then d1 = {x | a1 <= x <= B  U  A <= x <= b1}
    !=======================================================================
    real(kind=DP), dimension(1:), INTENT(IN)  :: d1, d2
    real(kind=DP), dimension(1:), INTENT(OUT) :: di
    integer(kind=I4B), INTENT(OUT) :: ni

    real(kind=DP), dimension(1:4) :: dk
    integer(kind=I4B) :: ik
    logical(kind=LGT) :: tr12, tr21, tr34, tr43, tr13, tr31, tr24, tr42, tr14, tr32
    !=======================================================================

    tr12 = (d1(1) < d1(2))
    tr21 = .NOT. tr12
    tr34 = (d2(1) < d2(2))
    tr43 = .NOT. tr34
    tr13 = (d1(1) < d2(1))
    tr31 = .NOT. tr13
    tr24 = (d1(2) < d2(2))
    tr42 = .NOT. tr24
    tr14 = (d1(1) < d2(2))
    tr32 = (d2(1) < d1(2))

    ik = 0
    dk(1:4) = -1.0e9_dp


    if ((tr34.AND.tr31.AND.tr14) .OR. (tr43.AND.(tr31.OR.tr14))) then
       ik = ik + 1
       dk(ik) = d1(1)  ! a1
    endif
    if ((tr12.AND.tr13.AND.tr32) .OR. (tr21.AND.(tr13.OR.tr32))) then
       ik = ik + 1
       dk(ik) = d2(1)  ! a2
    endif
    if ((tr34.AND.tr32.AND.tr24) .OR. (tr43.AND.(tr32.OR.tr24))) then
       ik = ik + 1
       dk(ik) = d1(2)  ! b1
    endif
    if ((tr12.AND.tr14.AND.tr42) .OR. (tr21.AND.(tr14.OR.tr42))) then
       ik = ik + 1
       dk(ik) =  d2(2)  ! b2
    endif

    di(1:4) = 0.0_dp
    select case (ik)
    case (0)
       ni = 0
    case (2)
       ni = 1
       di(1:2) = (/ dk(1), dk(2) /) ! [a1,b1] or [a1,b2] or [a2,b1] or [a2,b2]
    case (4)
       ni = 2
       di(1:4) = (/ dk(1), dk(4), dk(2), dk(3) /) ! [a1,b2] U [a2,b1]
    case default
       print*,"error in intrs_intrv", ik
       print*,dk
       print*,d1,d2
    end select

    return
  end subroutine intrs_intrv

  !=======================================================================
  subroutine pix2ang_ring(nside, ipix, theta, phi)
    !=======================================================================
    !     renders theta and phi coordinates of the nominal pixel center
    !     for the pixel number ipix (RING scheme)
    !     given the map resolution parameter nside
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: ipix, nside
    REAL(KIND=DP), INTENT(OUT) ::  theta, phi

    INTEGER(KIND=I4B) ::  nl2, nl4, npix, ncap, iring, iphi, ip, ipix1
    REAL(KIND=DP) ::  fodd, hip, fihip
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error ("nside out of range")
    npix = 12*nside**2       ! total number of points
    if (ipix <0 .or. ipix>npix-1) call fatal_error ("ipix out of range")

    ipix1 = ipix + 1 ! in {1, npix}
    nl2 = 2*nside
    ncap = 2*nside*(nside-1) ! points in each polar cap, =0 for nside =1

    if (ipix1 <= ncap) then ! North Polar cap -------------

       hip   = ipix1*0.5_dp
       fihip = AINT ( hip , kind=DP)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1 ! counted from North pole
       iphi  = ipix1 - 2*iring*(iring - 1)

       theta = ACOS( 1.0_dp - iring**2 / (3.0_dp*nside**2) )
       phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

    elseif (ipix1 <= nl2*(5*nside+1)) then ! Equatorial region ------

       ip    = ipix1 - ncap - 1
       nl4   = 4*nside
       iring = INT( ip / nl4 ) + nside ! counted from North pole
       iphi  = MODULO(ip,nl4) + 1

       fodd  = 0.5_dp * (1 + MODULO(iring+nside,2))  ! 1 if iring+nside is odd, 1/2 otherwise
       theta = ACOS( (nl2 - iring) / (1.5_dp*nside) )
       phi   = (real(iphi,kind=dp) - fodd) * PI /(2.0_dp*nside)

    else ! South Polar cap -----------------------------------

       ip    = npix - ipix1 + 1
       hip   = ip*0.5_dp
       fihip = AINT ( hip , kind=DP)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1     ! counted from South pole
       iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1))

       theta = ACOS( -1.0_dp + iring**2 / (3.0_dp*nside**2) )
       phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

    endif

    return
  end subroutine pix2ang_ring ! pix2ang_ring
!   !=======================================================================
!   subroutine pix2vec0_ring(nside, ipix, vector)
!     !=======================================================================
!     !     renders vector (x,y,z) coordinates of the nominal pixel center
!     !     for the pixel number ipix (RING scheme)
!     !     given the map resolution parameter nside
!     !=======================================================================
!     INTEGER(KIND=I4B), INTENT(IN) :: ipix, nside
!     REAL(KIND=DP), INTENT(OUT),dimension(1:) :: vector

!     INTEGER(KIND=I4B) :: nl2, nl4, npix, ncap, iring, iphi, ip, ipix1
!     REAL(KIND=DP) ::  fact1, fact2, fodd, hip, fihip, z, sth, phi
!     !-----------------------------------------------------------------------
!     if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
!     npix = 12*nside**2       ! total number of points
!     if (ipix <0 .or. ipix>npix-1) call fatal_error("ipix out of range")

!     ipix1 = ipix + 1 ! in {1, npix}
!     nl2 = 2*nside
!     nl4 = 4*nside
!     ncap = 2*nside*(nside-1) ! points in each polar cap, =0 for nside =1
!     fact1 = 1.5_dp*nside
!     fact2 = 3.0_dp*nside**2

!     if (ipix1 <= ncap) then ! North Polar cap -------------

!        hip   = ipix1/2.0_dp
!        fihip = AINT ( hip ,kind=DP)
!        iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1 ! counted from North pole
!        iphi  = ipix1 - 2*iring*(iring - 1)

!        z =  1.0_dp - iring**2 / fact2
!        phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

!     elseif (ipix1 <= nl2*(5*nside+1)) then ! Equatorial region ------

!        ip    = ipix1 - ncap - 1
!        iring = INT( ip / nl4 ) + nside ! counted from North pole
!        iphi  = MODULO(ip,nl4) + 1

!        fodd  = 0.5_dp * (1 + MODULO(iring+nside,2))  ! 1 if iring+nside is odd, 1/2 otherwise
!        z = (nl2 - iring) / fact1
!        phi   = (real(iphi,kind=dp) - fodd) * PI /(2.0_dp*nside)

!     else ! South Polar cap -----------------------------------

!        ip    = npix - ipix1 + 1
!        hip   = ip/2.0_dp
!        fihip = AINT ( hip ,kind=DP)
!        iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1     ! counted from South pole
!        iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1))

!        z = -1.0_dp + iring**2 / fact2
!        phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

!     endif

!     sth = SQRT((1.0_dp-z)*(1.0_dp+z))
!     vector(1) = sth * COS(phi)
!     vector(2) = sth * SIN(phi)
!     vector(3) = z

!     return
!   end subroutine pix2vec0_ring
  !=======================================================================
  subroutine pix2vec_ring(nside, ipix, vector, vertex)
    !=======================================================================
    !     renders vector (x,y,z) coordinates of the nominal pixel center
    !     for the pixel number ipix (RING scheme)
    !     given the map resolution parameter nside
    !     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
    !     in the order N,W,S,E
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN)                             :: ipix, nside
    REAL(KIND=DP),     INTENT(OUT),dimension(1:)              :: vector
    REAL(KIND=DP),     INTENT(OUT),dimension(1:,1:), optional :: vertex

    INTEGER(KIND=I4B) :: nl2, nl4, npix, ncap, iring, iphi, ip, ipix1
    REAL(KIND=DP) ::  fact1, fact2, fodd, hip, fihip, z, sth, phi

    real(kind=DP) :: phi_nv, phi_wv, phi_sv, phi_ev
    real(kind=DP) :: z_nv, z_sv, sth_nv, sth_sv
    real(kind=DP) :: hdelta_phi
    integer(kind=I4B) :: iphi_mod, iphi_rat
    logical(kind=LGT) :: do_vertex
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    npix = 12*nside**2       ! total number of points
    if (ipix <0 .or. ipix>npix-1) call fatal_error("ipix out of range")

    ipix1 = ipix + 1 ! in {1, npix}
    nl2 = 2*nside
    nl4 = 4*nside
    ncap = 2*nside*(nside-1) ! points in each polar cap, =0 for nside =1
    fact1 = 1.5_dp*nside
    fact2 = 3.0_dp*nside**2

    do_vertex = .false.
    if (present(vertex)) then
       if (size(vertex,dim=1) >= 3 .and. size(vertex,dim=2) >= 4) then
          do_vertex = .true.
       else
          call fatal_error(" pix2vec_ring : vertex array has wrong size ")
       endif
    endif

    phi_nv = 0.0_dp
    phi_sv = 0.0_dp
    if (ipix1 <= ncap) then ! North Polar cap -------------

       hip   = ipix1/2.0_dp
       fihip = AINT ( hip ,kind=DP)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1 ! counted from North pole
       iphi  = ipix1 - 2*iring*(iring - 1)

       z =  1.0_dp - iring**2 / fact2
       phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

       if (do_vertex) then
          hdelta_phi = PI/(4.0_dp*iring)   ! half pixel width
          z_nv = 1.0_dp - (iring-1)**2 / fact2
          z_sv = 1.0_dp - (iring+1)**2 / fact2
          iphi_mod = MODULO(iphi-1, iring) ! in {0,1,... iring-1}
          iphi_rat = (iphi-1) / iring      ! in {0,1,2,3}
          if (iring > 1) phi_nv = HALFPI * (iphi_rat +  iphi_mod   /real(iring-1,kind=dp))
          phi_sv                = HALFPI * (iphi_rat + (iphi_mod+1)/real(iring+1,kind=dp))
       endif


    elseif (ipix1 <= nl2*(5*nside+1)) then ! Equatorial region ------

       ip    = ipix1 - ncap - 1
       iring = INT( ip / nl4 ) + nside ! counted from North pole
       iphi  = MODULO(ip,nl4) + 1

       fodd  = 0.5_dp * (1 + MODULO(iring+nside,2))  ! 1 if iring+nside is odd, 1/2 otherwise
       z = (nl2 - iring) / fact1
       phi   = (real(iphi,kind=dp) - fodd) * PI /(2.0_dp*nside)

       if (do_vertex) then
          hdelta_phi = PI/(4.0_dp*nside)   ! half pixel width
          phi_nv = phi
          phi_sv = phi
          z_nv = (nl2 - iring +1) / fact1
          z_sv = (nl2 - iring -1) / fact1
          if (iring == nside) then ! northern transition
             z_nv = 1.0_dp - (nside-1)**2 / fact2
             iphi_mod = MODULO(iphi-1, nside) ! in {0,1,... nside-1}
             iphi_rat = (iphi-1) / nside      ! in {0,1,2,3}
             if (nside > 1) phi_nv = HALFPI * (iphi_rat +  iphi_mod   /real(nside-1,kind=dp))
          elseif (iring == 3*nside) then ! southern transition
             z_sv = -1.0_dp + (nside-1)**2 / fact2
             iphi_mod = MODULO(iphi-1, nside) ! in {0,1,... iring-1}
             iphi_rat = (iphi-1) / nside      ! in {0,1,2,3}
             if (nside > 1) phi_sv = HALFPI * (iphi_rat +  iphi_mod   /real(nside-1,kind=dp))
          endif
       endif

    else ! South Polar cap -----------------------------------

       ip    = npix - ipix1 + 1
       hip   = ip/2.0_dp
       fihip = AINT ( hip ,kind=DP)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1     ! counted from South pole
       iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1))

       z = -1.0_dp + iring**2 / fact2
       phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

       if (do_vertex) then
          hdelta_phi = PI/(4.0_dp*iring)   ! half pixel width
          z_nv = -1.0_dp + (iring+1)**2 / fact2
          z_sv = -1.0_dp + (iring-1)**2 / fact2
          iphi_mod = MODULO(iphi-1, iring) ! in {0,1,... iring-1}
          iphi_rat = (iphi-1) / iring      ! in {0,1,2,3}
          phi_nv                = HALFPI * (iphi_rat + (iphi_mod+1)/real(iring+1,kind=dp))
          if (iring > 1) phi_sv = HALFPI * (iphi_rat +  iphi_mod   /real(iring-1,kind=dp))
       endif

    endif

    ! pixel center
    sth = SQRT((1.0_dp-z)*(1.0_dp+z))
    vector(1) = sth * COS(phi)
    vector(2) = sth * SIN(phi)
    vector(3) = z

    if (do_vertex) then
       ! west vertex
       phi_wv      = phi - hdelta_phi
       vertex(1,2) = sth * COS(phi_wv)
       vertex(2,2) = sth * SIN(phi_wv)
       vertex(3,2) = z

       ! east vertex
       phi_ev      = phi + hdelta_phi
       vertex(1,4) = sth * COS(phi_ev)
       vertex(2,4) = sth * SIN(phi_ev)
       vertex(3,4) = z

       ! north vertex
       sth_nv = SQRT((1.0_dp-z_nv)*(1.0_dp+z_nv))
       vertex(1,1) = sth_nv * COS(phi_nv)
       vertex(2,1) = sth_nv * SIN(phi_nv)
       vertex(3,1) = z_nv

       ! south vertex
       sth_sv = SQRT((1.0_dp-z_sv)*(1.0_dp+z_sv))
       vertex(1,3) = sth_sv * COS(phi_sv)
       vertex(2,3) = sth_sv * SIN(phi_sv)
       vertex(3,3) = z_sv
    endif


    return
  end subroutine pix2vec_ring
  !=======================================================================
  subroutine ang2pix_ring(nside, theta, phi, ipix)
    !=======================================================================
    !     renders the pixel number ipix (RING scheme) for a pixel which contains
    !     a point on a sphere at coordinates theta and phi, given the map
    !     resolution parameter nside
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=I4B), INTENT(OUT) :: ipix
    REAL(KIND=DP), INTENT(IN) ::  theta, phi

    INTEGER(KIND=I4B) ::  nl4, jp, jm
    REAL(KIND=DP) ::  z, za, tt, tp, tmp, temp1, temp2
    INTEGER(KIND=I4B) ::  ir, ip, kshift

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error ("nside out of range")
    if (theta<0.0_dp .or. theta>pi)  then
       print*,"ANG2PIX_RING: theta : ",theta," is out of range [0, Pi]"
       call fatal_error
    endif

    z = COS(theta)
    za = ABS(z)
    tt = MODULO( phi, twopi) / halfpi  ! in [0,4)


    if ( za <= twothird ) then ! Equatorial region ------------------
       temp1 = nside*(.5_dp+tt)
       temp2 = nside*.75_dp*z
       jp = int(temp1-temp2) ! index of  ascending edge line
       jm = int(temp1+temp2) ! index of descending edge line

       ir = nside + 1 + jp - jm ! in {1,2n+1} (ring number counted from z=2/3)
       kshift = 1 - modulo(ir,2) ! kshift=1 if ir even, 0 otherwise

       nl4 = 4*nside
       ip = INT( ( jp+jm - nside + kshift + 1 ) / 2 ) ! in {0,4n-1}
       if (ip >= nl4) ip = ip - nl4

       ipix = 2*nside*(nside-1) + nl4*(ir-1) + ip

    else ! North & South polar caps -----------------------------

       tp = tt - INT(tt)      !MODULO(tt,1.0_dp)
       tmp = nside * SQRT( 3.0_dp*(1.0_dp - za) )

       jp = INT(tp          * tmp ) ! increasing edge line index
       jm = INT((1.0_dp - tp) * tmp ) ! decreasing edge line index

       ir = jp + jm + 1        ! ring number counted from the closest pole
       ip = INT( tt * ir )     ! in {0,4*ir-1}
       if (ip >= 4*ir) ip = ip - 4*ir

       if (z>0._dp) then
          ipix = 2*ir*(ir-1) + ip
       else
          ipix = 12*nside**2 - 2*ir*(ir+1) + ip
       endif

    endif

    return
  end subroutine ang2pix_ring

  !=======================================================================
  subroutine vec2pix_ring(nside, vector, ipix)
    !=======================================================================
    !     renders the pixel number ipix (RING scheme) for a pixel which contains
    !     a point on a sphere at coordinate vector (=x,y,z), given the map
    !     resolution parameter nside
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=I4B), INTENT(OUT) :: ipix
    REAL(KIND=DP), INTENT(IN), dimension(1:) :: vector

    INTEGER(KIND=I4B) :: nl2, nl4, ncap, npix, jp, jm, ipix1
    REAL(KIND=DP) ::  z, za, tt, tp, tmp, dnorm, phi
    INTEGER(KIND=I4B) :: ir, ip, kshift

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")

    dnorm = SQRT(vector(1)**2+vector(2)**2+vector(3)**2)
    z = vector(3) / dnorm
    phi = 0.0_dp
    if (vector(1) /= 0.0_dp .or. vector(2) /= 0.0_dp) &
         &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]

    za = ABS(z)
    if (phi < 0.0)     phi = phi + twopi ! phi in [0,2pi[
    tt = phi / halfpi   ! in [0,4)

    nl2 = 2*nside
    nl4 = 4*nside
    ncap  = nl2*(nside-1) ! number of pixels in the north polar cap
    npix  = 12*nside**2

    if ( za <= twothird ) then ! Equatorial region ------------------

       jp = INT(nside*(0.5_dp + tt - z*0.75_dp)) ! index of  ascending edge line
       jm = INT(nside*(0.5_dp + tt + z*0.75_dp)) ! index of descending edge line

       ir = nside + 1 + jp - jm ! in {1,2n+1} (ring number counted from z=2/3)
       kshift = 0
       if (MODULO(ir,2) == 0) kshift = 1 ! kshift=1 if ir even, 0 otherwise

       ip = INT( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1 ! in {1,4n}
       if (ip > nl4) ip = ip - nl4

       ipix1 = ncap + nl4*(ir-1) + ip

    else ! North & South polar caps -----------------------------

       tp = tt - INT(tt)      !MODULO(tt,1.0_dp)
       tmp = SQRT( 3.0_dp*(1.0_dp - za) )

       jp = INT( nside * tp          * tmp ) ! increasing edge line index
       jm = INT( nside * (1.0_dp - tp) * tmp ) ! decreasing edge line index

       ir = jp + jm + 1        ! ring number counted from the closest pole
       ip = INT( tt * ir ) + 1 ! in {1,4*ir}
       if (ip > 4*ir) ip = ip - 4*ir

       ipix1 = 2*ir*(ir-1) + ip
       if (z <= 0.0_dp) then
          ipix1 = npix - 2*ir*(ir+1) + ip
       endif

    endif

    ipix = ipix1 - 1 ! in {0, npix-1}

    return
  end subroutine vec2pix_ring
  !=======================================================================
  subroutine pix2ang_nest(nside, ipix, theta, phi)
    !=======================================================================
    !     renders theta and phi coordinates of the nominal pixel center
    !     for the pixel number ipix (NESTED scheme)
    !     given the map resolution parameter nside
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside, ipix
    REAL(KIND=DP), INTENT(OUT) :: theta, phi

    INTEGER(KIND=I4B) :: npix, npface, &
         &     ipf, ip_low, ip_trunc, ip_med, ip_hi, &
         &     jrt, jr, nr, jpt, jp, kshift, nl4
    REAL(KIND=DP) :: z, fn, fact1, fact2

    INTEGER(KIND=I4B) :: ix, iy, face_num
!     common /xy_nest/ ix, iy, face_num ! can be useful to calling routine

    ! coordinate of the lowest corner of each face
    INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    npix = 12 * nside**2
    if (ipix <0 .or. ipix>npix-1) call fatal_error("ipix out of range")

    !     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    fn = real(nside,kind=dp)
    fact1 = 1.0_dp/(3.0_dp*fn*fn)
    fact2 = 2.0_dp/(3.0_dp*fn)
    nl4   = 4*nside

    !     finds the face, and the number in the face
    npface = nside**2

    face_num = ipix/npface  ! face number in {0,11}
    ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    ip_low = MODULO(ipf,1024)       ! content of the last 10 bits
    ip_trunc =   ipf/1024        ! truncation of the last 10 bits
    ip_med = MODULO(ip_trunc,1024)  ! content of the next 10 bits
    ip_hi  =     ip_trunc/1024   ! content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}

    nr = nside                  ! equatorial region (the most frequent)
    z  = (2*nside-jr)*fact2
    kshift = MODULO(jr - nside, 2)
    if (jr < nside) then     ! north pole region
       nr = jr
       z = 1.0_dp - nr*nr*fact1
       kshift = 0
    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       z = - 1.0_dp + nr*nr*fact1
       kshift = 0
    endif
    theta = ACOS(z)

    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4

    phi = (jp - (kshift+1)*0.5_dp) * (halfpi / nr)

    return
  end subroutine pix2ang_nest
!   !=======================================================================
!   subroutine pix2vec_nest(nside, ipix, vector)
!     !=======================================================================
!     !     renders vector (x,y,z) coordinates of the nominal pixel center
!     !     for the pixel number ipix (NESTED scheme)
!     !     given the map resolution parameter nside
!     !=======================================================================
!     INTEGER(KIND=I4B), INTENT(IN) :: nside, ipix
!     REAL(KIND=DP), INTENT(OUT), dimension(1:) :: vector

!     INTEGER(KIND=I4B) :: npix, npface, &
!          &     ipf, ip_low, ip_trunc, ip_med, ip_hi, &
!          &     jrt, jr, nr, jpt, jp, kshift, nl4
!     REAL(KIND=DP) :: z, fn, fact1, fact2, sth, phi

!     INTEGER(KIND=I4B) ::  ix, iy, face_num
! !     common /xy_nest/ ix, iy, face_num ! can be useful to calling routine

!     ! coordinate of the lowest corner of each face
!     INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
!     INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2
!     !-----------------------------------------------------------------------
!     if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
!     npix = 12 * nside**2
!     if (ipix <0 .or. ipix>npix-1) call fatal_error("ipix out of range")

!     !     initiates the array for the pixel number -> (x,y) mapping
!     if (pix2x(1023) <= 0) call mk_pix2xy()

!     fn = real(nside,kind=dp)
!     fact1 = 1.0_dp/(3.0_dp*fn*fn)
!     fact2 = 2.0_dp/(3.0_dp*fn)
!     nl4   = 4*nside

!     !     finds the face, and the number in the face
!     npface = nside**2

!     face_num = ipix/npface  ! face number in {0,11}
!     ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}

!     !     finds the x,y on the face (starting from the lowest corner)
!     !     from the pixel number
!     ip_low = MODULO(ipf,1024)       ! content of the last 10 bits
!     ip_trunc =   ipf/1024        ! truncation of the last 10 bits
!     ip_med = MODULO(ip_trunc,1024)  ! content of the next 10 bits
!     ip_hi  =     ip_trunc/1024   ! content of the high weight 10 bits

!     ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
!     iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)

!     !     transforms this in (horizontal, vertical) coordinates
!     jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
!     jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

!     !     computes the z coordinate on the sphere
!     jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}

!     nr = nside                  ! equatorial region (the most frequent)
!     z  = (2*nside-jr)*fact2
!     kshift = MODULO(jr - nside, 2)
!     if (jr < nside) then     ! north pole region
!        nr = jr
!        z = 1.0_dp - nr*nr*fact1
!        kshift = 0
!     else if (jr > 3*nside) then ! south pole region
!        nr = nl4 - jr
!        z = - 1.0_dp + nr*nr*fact1
!        kshift = 0
!     endif
!     !ccc      theta = ACOS(z)

!     !     computes the phi coordinate on the sphere, in [0,2Pi]
!     jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
!     if (jp > nl4) jp = jp - nl4
!     if (jp < 1)   jp = jp + nl4

!     phi = (jp - (kshift+1)*0.5_dp) * (halfpi / nr)

!     sth = SQRT((1.0_dp-z)*(1.0_dp+z))
!     vector(1) = sth * COS(phi)
!     vector(2) = sth * SIN(phi)
!     vector(3) = z

!     return
!   end subroutine pix2vec_nest
  !=======================================================================
  subroutine pix2vec_nest(nside, ipix, vector, vertex)
    !=======================================================================
    !     renders vector (x,y,z) coordinates of the nominal pixel center
    !     for the pixel number ipix (NESTED scheme)
    !     given the map resolution parameter nside
    !     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
    !     in the order N,W,S,E
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside, ipix
    REAL(KIND=DP), INTENT(OUT), dimension(1:) :: vector
    REAL(KIND=DP),     INTENT(OUT),dimension(1:,1:), optional :: vertex

    INTEGER(KIND=I4B) :: npix, npface, &
         &     ipf, ip_low, ip_trunc, ip_med, ip_hi, &
         &     jrt, jr, nr, jpt, jp, kshift, nl4
    REAL(KIND=DP) :: z, fn, fact1, fact2, sth, phi

    INTEGER(KIND=I4B) ::  ix, iy, face_num
!     common /xy_nest/ ix, iy, face_num ! can be useful to calling routine

    ! coordinate of the lowest corner of each face
    INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2

    real(kind=DP) :: phi_nv, phi_wv, phi_sv, phi_ev, phi_up, phi_dn
    real(kind=DP) :: z_nv, z_sv, sth_nv, sth_sv
    real(kind=DP) :: hdelta_phi
    integer(kind=I4B) :: iphi_mod, iphi_rat
    logical(kind=LGT) :: do_vertex
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    npix = 12 * nside**2
    if (ipix <0 .or. ipix>npix-1) call fatal_error("ipix out of range")

    !     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    fn = real(nside,kind=dp)
    fact1 = 1.0_dp/(3.0_dp*fn*fn)
    fact2 = 2.0_dp/(3.0_dp*fn)
    nl4   = 4*nside

    do_vertex = .false.
    if (present(vertex)) then
       if (size(vertex,dim=1) >= 3 .and. size(vertex,dim=2) >= 4) then
          do_vertex = .true.
       else
          call fatal_error(" pix2vec_ring : vertex array has wrong size ")
       endif
    endif

    !     finds the face, and the number in the face
    npface = nside**2

    face_num = ipix/npface  ! face number in {0,11}
    ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    ip_low = MODULO(ipf,1024)       ! content of the last 10 bits
    ip_trunc =   ipf/1024        ! truncation of the last 10 bits
    ip_med = MODULO(ip_trunc,1024)  ! content of the next 10 bits
    ip_hi  =     ip_trunc/1024   ! content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}

    nr = nside                  ! equatorial region (the most frequent)
    z  = (2*nside-jr)*fact2
    kshift = MODULO(jr - nside, 2)
    if (do_vertex) then
       z_nv = (2*nside-jr+1)*fact2
       z_sv = (2*nside-jr-1)*fact2
       if (jr == nside) then ! northern transition
          z_nv =  1.0_dp - (nside-1)**2 * fact1
       elseif (jr == 3*nside) then  ! southern transition
          z_sv = -1.0_dp + (nside-1)**2 * fact1
       endif
    endif
    if (jr < nside) then     ! north pole region
       nr = jr
       z = 1.0_dp - nr*nr*fact1
       kshift = 0
       if (do_vertex) then
          z_nv = 1.0_dp - (nr-1)**2*fact1
          z_sv = 1.0_dp - (nr+1)**2*fact1
       endif
    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       z = - 1.0_dp + nr*nr*fact1
       kshift = 0
       if (do_vertex) then
          z_nv = - 1.0_dp + (nr+1)**2*fact1
          z_sv = - 1.0_dp + (nr-1)**2*fact1
       endif
    endif

    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4

    phi = (jp - (kshift+1)*0.5_dp) * (halfpi / nr)

    sth = SQRT((1.0_dp-z)*(1.0_dp+z))
    vector(1) = sth * COS(phi)
    vector(2) = sth * SIN(phi)
    vector(3) = z

    if (do_vertex) then
       phi_nv = phi
       phi_sv = phi

       phi_up = 0.0_dp
       iphi_mod = MODULO(jp-1, nr) ! in {0,1,... nr-1}
       iphi_rat = (jp-1) / nr      ! in {0,1,2,3}
       if (nr > 1) phi_up = HALFPI * (iphi_rat +  iphi_mod   /real(nr-1,kind=dp))
       phi_dn             = HALFPI * (iphi_rat + (iphi_mod+1)/real(nr+1,kind=dp))
       if (jr < nside) then            ! North polar cap
          phi_nv = phi_up
          phi_sv = phi_dn
       else if (jr > 3*nside) then     ! South polar cap
          phi_nv = phi_dn
          phi_sv = phi_up
       else if (jr == nside) then      ! North transition
          phi_nv = phi_up
       else if (jr == 3*nside) then    ! South transition
          phi_sv = phi_up
       endif

       hdelta_phi = PI / (4.0_dp*nr)

       ! west vertex
       phi_wv      = phi - hdelta_phi
       vertex(1,2) = sth * COS(phi_wv)
       vertex(2,2) = sth * SIN(phi_wv)
       vertex(3,2) = z

       ! east vertex
       phi_ev      = phi + hdelta_phi
       vertex(1,4) = sth * COS(phi_ev)
       vertex(2,4) = sth * SIN(phi_ev)
       vertex(3,4) = z

       ! north vertex
       sth_nv = SQRT((1.0_dp-z_nv)*(1.0_dp+z_nv))
       vertex(1,1) = sth_nv * COS(phi_nv)
       vertex(2,1) = sth_nv * SIN(phi_nv)
       vertex(3,1) = z_nv

       ! south vertex
       sth_sv = SQRT((1.0_dp-z_sv)*(1.0_dp+z_sv))
       vertex(1,3) = sth_sv * COS(phi_sv)
       vertex(2,3) = sth_sv * SIN(phi_sv)
       vertex(3,3) = z_sv
    endif

    return
  end subroutine pix2vec_nest
  !=======================================================================
  subroutine ang2pix_nest(nside, theta, phi, ipix)
    !=======================================================================
    !     renders the pixel number ipix (NESTED scheme) for a pixel which contains
    !     a point on a sphere at coordinates theta and phi, given the map
    !     resolution parametr nside
    !
    !     the computation is made to the highest resolution available (nside=8192)
    !     and then degraded to that required (by integer division)
    !     this doesn't cost more, and it makes sure
    !     that the treatement of round-off will be consistent
    !     for every resolution
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=I4B), INTENT(OUT) :: ipix
    REAL(KIND=DP), INTENT(IN) ::  theta, phi

    REAL(KIND=DP) ::  z, za, tt, tp, tmp
    INTEGER(KIND=I4B) :: jp, jm, ifp, ifm, face_num, &
         &     ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    if (theta<0.0_dp .or. theta>pi)  then
       print*,"ANG2PIX_NEST: theta : ",theta," is out of range [0,Pi]"
       call fatal_error
    endif
    if (x2pix(128) <= 0) call mk_xy2pix()

    z  = COS(theta)
    za = ABS(z)
    tt = MODULO(phi, twopi) / halfpi  ! in [0,4[

    if (za <= twothird) then ! equatorial region

       !        (the index of edge lines increase when the longitude=phi goes up)
       jp = INT(ns_max*(0.5_dp + tt - z*0.75_dp)) !  ascending edge line index
       jm = INT(ns_max*(0.5_dp + tt + z*0.75_dp)) ! descending edge line index

       !        finds the face
       ifp = jp / ns_max  ! in {0,4}
       ifm = jm / ns_max
       if (ifp == ifm) then          ! faces 4 to 7
          face_num = MODULO(ifp,4) + 4
       else if (ifp < ifm) then     ! (half-)faces 0 to 3
          face_num = MODULO(ifp,4)
       else                            ! (half-)faces 8 to 11
          face_num = MODULO(ifm,4) + 8
       endif

       ix = MODULO(jm, ns_max)
       iy = ns_max - MODULO(jp, ns_max) - 1

    else ! polar region, za > 2/3

       ntt = INT(tt)
       if (ntt >= 4) ntt = 3
       tp = tt - ntt
       tmp = SQRT( 3.0_dp*(1.0_dp - za) )  ! in ]0,1]

       !        (the index of edge lines increase when distance from the closest pole goes up)
       jp = INT( ns_max * tp          * tmp ) ! line going toward the pole as phi increases
       jm = INT( ns_max * (1.0_dp - tp) * tmp ) ! that one goes away of the closest pole
       jp = MIN(ns_max-1, jp) ! for points too close to the boundary
       jm = MIN(ns_max-1, jm)

       !        finds the face and pixel's (x,y)
       if (z >= 0) then
          face_num = ntt  ! in {0,3}
          ix = ns_max - jm - 1
          iy = ns_max - jp - 1
       else
          face_num = ntt + 8 ! in {8,11}
          ix =  jp
          iy =  jm
       endif

       !         print*,z,face_num,ix,iy
    endif

    ix_low = MODULO(ix,128)
    ix_hi  =     ix/128
    iy_low = MODULO(iy,128)
    iy_hi  =     iy/128

    ipf =  (x2pix(ix_hi +1)+y2pix(iy_hi +1)) * (128 * 128) &
         &     + (x2pix(ix_low+1)+y2pix(iy_low+1))

    ipf = ipf / ( ns_max/nside ) **2  ! in {0, nside**2 - 1}

    ipix = ipf + face_num* nside **2    ! in {0, 12*nside**2 - 1}

    return
  end subroutine ang2pix_nest
  !=======================================================================
  subroutine vec2pix_nest(nside, vector, ipix)
    !=======================================================================
    !     renders the pixel number ipix (NESTED scheme) for a pixel which contains
    !     a point on a sphere at coordinate vector (=x,y,z), given the map
    !     resolution parameter nside
    !
    !     the computation is made to the highest resolution available (nside=8192)
    !     and then degraded to that required (by integer division)
    !     this doesn't cost more, and it makes sure
    !     that the treatement of round-off will be consistent
    !     for every resolution
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=I4B), INTENT(OUT) :: ipix
    REAL(KIND=DP), INTENT(IN), dimension(1:) ::  vector

    REAL(KIND=DP) ::  z, za, tt, tp, tmp, dnorm, phi
    INTEGER(KIND=I4B) ::  jp, jm, ifp, ifm, face_num, &
         &     ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    if (x2pix(128) <= 0) call mk_xy2pix()

    dnorm = SQRT(vector(1)**2+vector(2)**2+vector(3)**2)
    z = vector(3) / dnorm
    phi = 0.0_dp
    if (vector(1) /= 0.0_dp .or. vector(2) /= 0.0_dp) &
         &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]

    za = ABS(z)
    if (phi < 0.0)    phi = phi + twopi ! phi in [0,2pi[
    tt = phi / halfpi ! in [0,4[

    if (za <= twothird) then ! equatorial region

       !        (the index of edge lines increase when the longitude=phi goes up)
       jp = INT(ns_max*(0.5_dp + tt - z*0.75_dp)) !  ascending edge line index
       jm = INT(ns_max*(0.5_dp + tt + z*0.75_dp)) ! descending edge line index

       !        finds the face
       ifp = jp / ns_max  ! in {0,4}
       ifm = jm / ns_max
       if (ifp == ifm) then          ! faces 4 to 7
          face_num = MODULO(ifp,4) + 4
       else if (ifp < ifm) then     ! (half-)faces 0 to 3
          face_num = MODULO(ifp,4)
       else                            ! (half-)faces 8 to 11
          face_num = MODULO(ifm,4) + 8
       endif

       ix = MODULO(jm, ns_max)
       iy = ns_max - MODULO(jp, ns_max) - 1

    else ! polar region, za > 2/3

       ntt = INT(tt)
       if (ntt >= 4) ntt = 3
       tp = tt - ntt
       tmp = SQRT( 3.0_dp*(1.0_dp - za) )  ! in ]0,1]

       !        (the index of edge lines increase when distance from the closest pole goes up)
       jp = INT( ns_max * tp          * tmp ) ! line going toward the pole as phi increases
       jm = INT( ns_max * (1.0_dp - tp) * tmp ) ! that one goes away of the closest pole
       jp = MIN(ns_max-1, jp) ! for points too close to the boundary
       jm = MIN(ns_max-1, jm)

       !        finds the face and pixel's (x,y)
       if (z >= 0) then
          face_num = ntt  ! in {0,3}
          ix = ns_max - jm - 1
          iy = ns_max - jp - 1
       else
          face_num = ntt + 8 ! in {8,11}
          ix =  jp
          iy =  jm
       endif

       !         print*,z,face_num,ix,iy
    endif

    ix_low = MODULO(ix,128)
    ix_hi  =     ix/128
    iy_low = MODULO(iy,128)
    iy_hi  =     iy/128

    ipf =  (x2pix(ix_hi +1)+y2pix(iy_hi +1)) * (128 * 128) &
         &     + (x2pix(ix_low+1)+y2pix(iy_low+1))

    ipf = ipf / ( ns_max/nside ) **2  ! in {0, nside**2 - 1}

    ipix = ipf + face_num* nside **2    ! in {0, 12*nside**2 - 1}

    return
  end subroutine vec2pix_nest

  !*************************************************************
  !
  !                    MAP Manipulations
  !
  !*************************************************************
  subroutine warning_oldbounds(code, cos_theta_cut, zbounds)
    character(len=*), intent(in)          :: code
    real(DP), intent(in)                  :: cos_theta_cut
    real(DP), intent(out), dimension(1:2) :: zbounds

    if (cos_theta_cut <= 0.0_dp) then ! no cut
       zbounds(1) =  -1.0_dp
       zbounds(2) =   1.0_dp
    else
       zbounds(1) =   cos_theta_cut
       zbounds(2) = - cos_theta_cut
    endif
    print*,' -------------------------------------'
    print*,'WARNING: obsolete interface to '//code
    print*,'    cos_theta_cut currently a DP scalar with value'
    write(*,9000) '    ',cos_theta_cut
    print*,'    shoud now be replaced with a 2-element vector with values:'
    write(*,9001) '    ',zbounds(1),zbounds(2)
    print*,'    See documentation for details.'
    print*,' -------------------------------------'
9000 format (a,g12.6)
9001 format (a,g12.6,g12.6)

    return
  end subroutine warning_oldbounds
  !==============================================================
  ! REMOVE_DIPOLE( nside, map, ordering, degree, multipoles, zbounds, fmissval, mask)
  !
  ! removes monopole (and dipole) from a map
  !
  ! Nside:     I4,       IN   : Healpix resolution parameter
  ! map:       KMAP, array,INOUT: Heapix map (see Notes below)
  ! ordering:  I4,       IN:   Healpix scheme 1:RING, 2: NESTED
  ! degree:    I4,       IN:   multipole to remove, 1: monopole, 2: monopole and dipole
  ! multipoles:R8, array,OUT:  value of monopole and dipole
  ! zbounds   :R8,      IN:  2-el vector
  ! range of z in [-1,1] on which to estimate monopole and dipole
  ! fmissval:  KMAP, Option, IN: value used to flag bad pixel on input, default=-1.6375e30
  !                            Pixels with map = fmissval are not used for fit
  ! mask    :  KMAP, Option, IN: Pixels with |mask|<1.e-10 are not used for fit
  !                              others are kept as is
  !                               note : the mask in NOT applied to the map
  !
  ! KMAP: either R4 or R8
  !
  ! note : if degree= 1, or 2, the map is modified on output
  !     * the monopole (and dipole) is/are removed
  !     * pixels within the symmetric cut parameterized
  !       by cos_theta_cut are set to fmissval (or its default value)
  !  if degree = 0, nothing is done
  !  all other values of degree are invalid
  !
  ! v1.0, EH, Caltech, Jan-2002, based on homonyme IDL routine
  !==============================================================
  subroutine remove_dipole_real( nside, map, ordering, degree, multipoles, zbounds, fmissval, mask)
    !============================================================
    use num_rec, only : dsvdcmp, dsvbksb
    ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_REAL"
    integer(kind=I4B),     parameter :: KMAP = SP
    real   (kind=KMAP),    parameter :: fbad_value = -1.6375e30_KMAP
    !
    include 'remove_dipole_inc.f90'
  end subroutine remove_dipole_real
  !==============================================================
  subroutine remove_dipole_double( nside, map, ordering, degree, multipoles, zbounds, fmissval, mask)
    !============================================================
    use num_rec, only : dsvdcmp, dsvbksb
    ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_DOUBLE"
    integer(kind=I4B),     parameter :: KMAP = DP
    real   (kind=KMAP),    parameter :: fbad_value = -1.6375e30_KMAP
    !
    include 'remove_dipole_inc.f90'
  end subroutine remove_dipole_double

  !==============================================================
  subroutine remove_dipole_real_old( nside, map, ordering, degree, multipoles, cos_theta_cut, fmissval, mask)
    !============================================================
    ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_REAL"
    integer(kind=I4B),     parameter :: KMAP = SP
    ! dummy
    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree
    real   (kind=DP),   dimension(0:),  intent(out)   :: multipoles
    real   (kind=KMAP), dimension(0:),  intent(inout) :: map
    real   (kind=DP),                   intent(in)    :: cos_theta_cut
    real   (kind=KMAP),                 intent(in), optional :: fmissval
    real   (kind=KMAP), dimension(0:),  intent(in), optional :: mask
    ! local
    real(DP),                          dimension(1:2) :: zbounds

    call warning_oldbounds(code, cos_theta_cut, zbounds)
    call remove_dipole(nside, map, ordering, degree, multipoles, zbounds, fmissval, mask)
  end subroutine remove_dipole_real_old
  !==============================================================
  subroutine remove_dipole_double_old( nside, map, ordering, degree, multipoles, cos_theta_cut, fmissval, mask)
    !============================================================
     ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_DOUBLE"
    integer(kind=I4B),     parameter :: KMAP = DP
    ! dummy
    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree
    real   (kind=DP),   dimension(0:),  intent(out)   :: multipoles
    real   (kind=KMAP), dimension(0:),  intent(inout) :: map
    real   (kind=DP),                   intent(in)    :: cos_theta_cut
    real   (kind=KMAP),                 intent(in), optional :: fmissval
    real   (kind=KMAP), dimension(0:),  intent(in), optional :: mask
    ! local
    real(DP),                          dimension(1:2) :: zbounds

    call warning_oldbounds(code, cos_theta_cut, zbounds)
    call remove_dipole(nside, map, ordering, degree, multipoles, zbounds, fmissval, mask) 
  end subroutine remove_dipole_double_old

  !==============================================================
  subroutine remove_dipole_real_v12( nside, map, nmaps, ordering, degree, multipoles, cos_theta_cut, fmissval, mask)
    !============================================================
    ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_REAL"
    integer(kind=I4B),     parameter :: KMAP = SP
    ! dummy
    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree, nmaps
    real   (kind=DP),   dimension(0:),  intent(out)   :: multipoles
    real   (kind=KMAP), dimension(0:),  intent(inout) :: map
    real   (kind=DP),                   intent(in)    :: cos_theta_cut
    real   (kind=KMAP),                 intent(in), optional :: fmissval
    real   (kind=KMAP), dimension(0:),  intent(in), optional :: mask
    ! local
    real(DP),                          dimension(1:2) :: zbounds

    print*,'=========================================================='
    print*,'WARNING: Interface to remove_dipole has changed'
    print*,' from'
    print*,'call remove_dipole(nside, map, NMAPS, ordering, degree, multipoles, COS_THETA_CUT, fmissval, mask)'
    print*,' to'
    print*,'call remove_dipole(nside, map,        ordering, degree, multipoles, ZBOUNDS,       fmissval, mask)'
    print*,'=========================================================='
    
    call warning_oldbounds(code, cos_theta_cut, zbounds)
    call remove_dipole(nside, map, ordering, degree, multipoles, zbounds, fmissval, mask)
  end subroutine remove_dipole_real_v12
  !==============================================================
  subroutine remove_dipole_double_v12( nside, map, nmaps, ordering, degree, multipoles, cos_theta_cut, fmissval, mask)
    !============================================================
     ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_DOUBLE"
    integer(kind=I4B),     parameter :: KMAP = DP
    ! dummy
    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree, nmaps
    real   (kind=DP),   dimension(0:),  intent(out)   :: multipoles
    real   (kind=KMAP), dimension(0:),  intent(inout) :: map
    real   (kind=DP),                   intent(in)    :: cos_theta_cut
    real   (kind=KMAP),                 intent(in), optional :: fmissval
    real   (kind=KMAP), dimension(0:),  intent(in), optional :: mask
    ! local
    real(DP),                          dimension(1:2) :: zbounds

    print*,'=========================================================='
    print*,'WARNING: Interface to remove_dipole has changed'
    print*,' from'
    print*,'call remove_dipole(nside, map, NMAPS, ordering, degree, multipoles, COS_THETA_CUT, fmissval, mask)'
    print*,' to'
    print*,'call remove_dipole(nside, map,        ordering, degree, multipoles, ZBOUNDS,       fmissval, mask)'
    print*,'=========================================================='

    call warning_oldbounds(code, cos_theta_cut, zbounds)
    call remove_dipole(nside, map, ordering, degree, multipoles, zbounds, fmissval, mask) 
  end subroutine remove_dipole_double_v12

  !=======================================================================
  ! ADD_DIPOLE( nside, map, ordering, degree, multipoles, fmissval)
  !
  ! removes monopole (and dipole) from a map
  !
  ! Nside:     I4,       IN   : Healpix resolution parameter
  ! map:       KMAP, array,INOUT: Heapix map (see Notes below)
  ! ordering:  I4,       IN:   Healpix scheme 1:RING, 2: NESTED
  ! degree:    I4,       IN:   multipole to remove, 1: monopole, 2: monopole and dipole
  ! multipoles:R8, array,IN:  value of monopole and dipole
  ! fmissval:  KMAP, Option, IN: value used to flag bad pixel on input, default=-1.6375e30
  !                            Pixels with map = fmissval are left unchanged
  !=======================================================================
  subroutine add_dipole_real(nside, map, ordering, degree, multipoles, fmissval)
    !=======================================================================
    ! single precision
    !=======================================================================
    character(len=*),      parameter :: code = "ADD_DIPOLE_REAL"
    integer(kind=I4B),     parameter :: KMAP = SP
    real   (kind=KMAP),    parameter :: fbad_value = -1.6375e30_KMAP
    ! dummy
    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree
    real   (kind=DP),   dimension(0:),  intent(in)    :: multipoles
    real   (kind=KMAP), dimension(0:),  intent(inout) :: map
    real   (kind=KMAP),                 intent(in), optional :: fmissval

    real   (kind=KMAP)                :: fmiss_effct
    integer(kind=i4b)                 :: ipix, npix
    logical(lgt)                      :: dodipole !, do_mask
    real(kind=dp), dimension(1:3)     :: vec
    !=======================================================================
    
    npix = nside2npix(nside)
    fmiss_effct = fbad_value
    if (present(fmissval)) fmiss_effct = fmissval

    if (degree == 0) then
       print*," No monopole nor dipole to add"
       return
    elseif (degree == 1) then
       dodipole = .false.
    else if (degree == 2) then
       dodipole = .true.
    else
       print*,code//"> degree can only be "
       print*,"      1: monopole (l=0) addition or "
       print*,"      2: monopole and dipole (l=0,1) addition"
       print*,code//"> ABORT ! "
       call fatal_error
    endif

    do ipix = 0, npix-1
       if ( abs(map(ipix) - fmiss_effct) <= abs(1.e-5*fmiss_effct) ) goto 20
!        if (do_mask) then
!           if (abs(mask(ipix)) <= 1.e-10) goto 20
!        endif
       map(ipix) = map(ipix) + multipoles(0)
       if (dodipole) then
          ! computes dipole basis functions
          ! pixel -> vector
          if (ordering == 1) call pix2vec_ring( nside, ipix, vec)
          if (ordering == 2) call pix2vec_nest( nside, ipix, vec)
          map(ipix) = map(ipix) + sum(multipoles(1:3) * vec(1:3))
       endif

20     continue
    enddo

    return
  end subroutine add_dipole_real
  !=======================================================================
  subroutine add_dipole_double(nside, map, ordering, degree, multipoles, fmissval)
    !=======================================================================
    ! single precision
    !=======================================================================
    character(len=*),      parameter :: code = "ADD_DIPOLE_DOUBLE"
    integer(kind=I4B),     parameter :: KMAP = DP
    real   (kind=KMAP),    parameter :: fbad_value = -1.6375e30_KMAP
    ! dummy
    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree
    real   (kind=DP),   dimension(0:),  intent(in)    :: multipoles
    real   (kind=KMAP), dimension(0:),  intent(inout) :: map
    real   (kind=KMAP),                 intent(in), optional :: fmissval

    real   (kind=KMAP)                :: fmiss_effct
    integer(kind=i4b)                 :: ipix, npix
    logical(lgt)                      :: dodipole !, do_mask
    real(kind=dp), dimension(1:3)     :: vec
    !=======================================================================
    
    npix = nside2npix(nside)
    fmiss_effct = fbad_value
    if (present(fmissval)) fmiss_effct = fmissval

    if (degree == 0) then
       print*," No monopole nor dipole to add"
       return
    elseif (degree == 1) then
       dodipole = .false.
    else if (degree == 2) then
       dodipole = .true.
    else
       print*,code//"> degree can only be "
       print*,"      1: monopole (l=0) addition or "
       print*,"      2: monopole and dipole (l=0,1) addition"
       print*,code//"> ABORT ! "
       call fatal_error
    endif

    do ipix = 0, npix-1
       if ( abs(map(ipix) - fmiss_effct) <= abs(1.e-5*fmiss_effct) ) goto 20
!        if (do_mask) then
!           if (abs(mask(ipix)) <= 1.e-10) goto 20
!        endif
       map(ipix) = map(ipix) + multipoles(0)
       if (dodipole) then
          ! computes dipole basis functions
          ! pixel -> vector
          if (ordering == 1) call pix2vec_ring( nside, ipix, vec)
          if (ordering == 2) call pix2vec_nest( nside, ipix, vec)
          map(ipix) = map(ipix) + sum(multipoles(1:3) * vec(1:3))
       endif

20     continue
    enddo

    return
  end subroutine add_dipole_double
  !====================================================
  ! medfiltmap
  !   compute the median filtered map of a given Healpix map
  !   in_map: SP/DP  input Healpix full sky map
  !   radius: DP     radius in Radians
  !   med_map: SP/DP output Healpix full sky map
  !   nest:    I4B, OPT   either 0 (ring scheme) or 1 (nested scheme)
  !   fmissval:SP/DP, OPT sentinel value given to missing pixels
  !   fill_holes: LGT, OPT 
  !=================================================================
  subroutine medfiltmap_S( in_map, radius, med_map, nest, fmissval, fill_holes)
  !=================================================================
    use statistics, only: median
    integer(I4B), parameter :: KMAP = SP
    !
    real(KMAP), dimension(0:),intent(in), target       :: in_map
    real(DP),                 intent(in)               :: radius
    real(KMAP), dimension(0:),intent(out)              :: med_map
    integer(I4B),             intent(in),  optional    :: nest
    real(KMAP),               intent(in),  optional    :: fmissval
    logical(LGT),             intent(in),  optional    :: fill_holes
    !
    integer(I4B) :: nside, npix, np, p, nlist, status
    logical(LGT) :: do_nest, do_fill
    real(DP)     :: fraction
    real(KMAP)   :: fmissval_in
    integer(I4B), dimension(:),  allocatable :: listpix
    real(DP),     dimension(1:3)             :: vector
    character(len=*), parameter :: code = 'medfiltmap'
    !-----------------------------------------------
    npix = size(in_map)
    nside = npix2nside(npix)
    call assert(nside > 0, code//": invalid map size")

    fraction = 0.5_DP * (1.0_dp - cos(radius))
    np = npix * fraction * 1.2 + 50
    allocate(listpix(0:np-1),stat=status)
    call assert_alloc(status,code,'listpix')

    do_nest = .false.
    if (present(nest)) then
       call assert(nest>=0 .and. nest <=1,code//': invalid NEST flag')
       do_nest = (nest == 1)
    endif

    do_fill = .false.
    if (present(fill_holes)) do_fill = fill_holes

    fmissval_in = hpx_Sbadval
    if (present(fmissval)) fmissval_in = fmissval

    do p = 0, npix-1
       ! find pixel location
       if (do_nest) then
          call pix2vec_nest( nside, p, vector)
       else
          call pix2vec_ring( nside, p, vector)
       endif

       ! find disc centered on pixel
       call query_disc(nside, vector, radius, listpix, nlist, nest=nest)

       if (do_fill .or. abs(in_map(p)-fmissval_in) > abs(fmissval_in*1.e-7)) then
          med_map(p) = median(in_map(listpix(0:nlist-1)), badval = fmissval_in, even= .true.)
       else
          med_map(p) = in_map(p)
       endif
    enddo

    deallocate(listpix)
    return
  end subroutine medfiltmap_S
  !=================================================================
  subroutine medfiltmap_D( in_map, radius, med_map, nest, fmissval, fill_holes)
  !=================================================================
    use statistics, only: median
    integer(I4B), parameter :: KMAP = DP
    !
    real(KMAP), dimension(0:),intent(in), target       :: in_map
    real(DP),                 intent(in)               :: radius
    real(KMAP), dimension(0:),intent(out)              :: med_map
    integer(I4B),             intent(in),  optional    :: nest
    real(KMAP),               intent(in),  optional    :: fmissval
    logical(LGT),             intent(in),  optional    :: fill_holes
    !
    integer(I4B) :: nside, npix, np, p, nlist, status
    logical(LGT) :: do_nest, do_fill
    real(DP)     :: fraction
    real(KMAP)   :: fmissval_in
    integer(I4B), dimension(:),  allocatable :: listpix
    real(DP),     dimension(1:3)             :: vector
    character(len=*), parameter :: code = 'medfiltmap'
    !-----------------------------------------------
    npix = size(in_map)
    nside = npix2nside(npix)
    call assert(nside > 0, code//": invalid map size")

    fraction = 0.5_DP * (1.0_dp - cos(radius))
    np = npix * fraction * 1.2 + 50
    allocate(listpix(0:np-1),stat=status)
    call assert_alloc(status,code,'listpix')

    do_nest = .false.
    if (present(nest)) then
       call assert(nest>=0 .and. nest <=1,code//': invalid NEST flag')
       do_nest = (nest == 1)
    endif

    do_fill = .false.
    if (present(fill_holes)) do_fill = fill_holes

    fmissval_in = hpx_Dbadval
    if (present(fmissval)) fmissval_in = fmissval

    do p = 0, npix-1
       ! find pixel location
       if (do_nest) then
          call pix2vec_nest( nside, p, vector)
       else
          call pix2vec_ring( nside, p, vector)
       endif

       ! find disc centered on pixel
       call query_disc(nside, vector, radius, listpix, nlist, nest=nest)

       if (do_fill .or. abs(in_map(p)-fmissval_in) > abs(fmissval_in*1.e-7)) then
          med_map(p) = median(in_map(listpix(0:nlist-1)), badval = fmissval_in, even= .true.)
       else
          med_map(p) = in_map(p)
       endif
    enddo

    deallocate(listpix)
    return
  end subroutine medfiltmap_D
  !**************************************************************
  ! following 2x3 routines: out of place conversions for 1D maps
  ! if size(map) = Npix, peak memory = 2*Npix
  ! These routines are parallelized for shared memory architecture 
  !   (using OpenMP directives)
  !**************************************************************
  !=======================================================================
  !     makes the conversion NEST to RING
  !=======================================================================
  subroutine convert_nest2ring_int_1d(nside, map)
    !=======================================================================
    integer(kind=I4B),   parameter :: KMAP = I4B
    integer(kind=I4B),                  intent(IN) :: nside
    integer(kind=KMAP),  dimension(0:), intent(inout) ::  map

    integer(kind=I4B) :: npix, ipn, ipr
    integer(kind=KMAP),  dimension(:), allocatable :: map_tmp
    !=======================================================================
    npix = 12*nside*nside
    allocate(map_tmp(0:npix-1))

!$OMP parallel default(none) &
!$OMP   shared(map, map_tmp, npix, nside) &
!$OMP   private(ipr, ipn)
!$OMP do schedule(dynamic,64)
    do ipn = 0, npix-1
       call nest2ring(nside, ipn, ipr)
       map_tmp(ipr) = map(ipn)
    enddo
!$OMP end do
!$OMP end parallel
    map = map_tmp

    deallocate(map_tmp)
    return
  end subroutine convert_nest2ring_int_1d
  !=======================================================================
  subroutine convert_nest2ring_real_1d(nside, map)
    !=======================================================================
    integer(kind=I4B),   parameter :: KMAP = SP
    integer(kind=I4B),                  intent(IN) :: nside
    real   (kind=KMAP),  dimension(0:), intent(inout) ::  map

    integer(kind=I4B) :: npix, ipn, ipr
    real   (kind=KMAP),  dimension(:), allocatable :: map_tmp
    !=======================================================================
    npix = 12*nside*nside
    allocate(map_tmp(0:npix-1))

!$OMP parallel default(none) &
!$OMP   shared(map, map_tmp, npix, nside) &
!$OMP   private(ipr, ipn)
!$OMP do schedule(dynamic,64)
    do ipn = 0, npix-1
       call nest2ring(nside, ipn, ipr)
       map_tmp(ipr) = map(ipn)
    enddo
!$OMP end do
!$OMP end parallel
    map = map_tmp

    deallocate(map_tmp)
    return
  end subroutine convert_nest2ring_real_1d
  !=======================================================================
  subroutine convert_nest2ring_double_1d(nside, map)
    !=======================================================================
    integer(kind=I4B),   parameter :: KMAP = DP
    integer(kind=I4B),                  intent(IN) :: nside
    real   (kind=KMAP),  dimension(0:), intent(inout) ::  map

    integer(kind=I4B) :: npix, ipn, ipr
    real   (kind=KMAP),  dimension(:), allocatable :: map_tmp
    !=======================================================================
    npix = 12*nside*nside
    allocate(map_tmp(0:npix-1))

!$OMP parallel default(none) &
!$OMP   shared(map, map_tmp, npix, nside) &
!$OMP   private(ipr, ipn)
!$OMP do schedule(dynamic,64)
    do ipn = 0, npix-1
       call nest2ring(nside, ipn, ipr)
       map_tmp(ipr) = map(ipn)
    enddo
!$OMP end do
!$OMP end parallel
    map = map_tmp

    deallocate(map_tmp)
    return
  end subroutine convert_nest2ring_double_1d
  !=======================================================================
  !     makes the conversion RING to NEST
  !=======================================================================
  subroutine convert_ring2nest_int_1d(nside, map)
    !=======================================================================
    integer(kind=I4B),   parameter :: KMAP = I4B
    integer(kind=I4B),                  intent(IN) :: nside
    integer(kind=KMAP),  dimension(0:), intent(inout) ::  map

    integer(kind=I4B) :: npix, ipn, ipr
    integer(kind=KMAP),  dimension(:), allocatable :: map_tmp
    !=======================================================================
    npix = 12*nside*nside
    allocate(map_tmp(0:npix-1))

!$OMP parallel default(none) &
!$OMP   shared(map, map_tmp, npix, nside) &
!$OMP   private(ipr, ipn)
!$OMP do schedule(dynamic,64)
    do ipr = 0, npix-1
       call ring2nest(nside, ipr, ipn)
       map_tmp(ipn) = map(ipr)
    enddo
!$OMP end do
!$OMP end parallel
    map = map_tmp

    deallocate(map_tmp)
    return
  end subroutine convert_ring2nest_int_1d
  !=======================================================================
  subroutine convert_ring2nest_real_1d(nside, map)
    !=======================================================================
    integer(kind=I4B),   parameter :: KMAP = SP
    integer(kind=I4B),                  intent(IN) :: nside
    real   (kind=KMAP),  dimension(0:), intent(inout) ::  map

    integer(kind=I4B) :: npix, ipn, ipr
    real   (kind=KMAP),  dimension(:), allocatable :: map_tmp
    !=======================================================================
    npix = 12*nside*nside
    allocate(map_tmp(0:npix-1))

!$OMP parallel default(none) &
!$OMP   shared(map, map_tmp, npix, nside) &
!$OMP   private(ipr, ipn)
!$OMP do schedule(dynamic,64)
    do ipr = 0, npix-1
       call ring2nest(nside, ipr, ipn)
       map_tmp(ipn) = map(ipr)
    enddo
!$OMP end do
!$OMP end parallel
    map = map_tmp

    deallocate(map_tmp)
    return
  end subroutine convert_ring2nest_real_1d
  !=======================================================================
  subroutine convert_ring2nest_double_1d(nside, map)
    !=======================================================================
    integer(kind=I4B),   parameter :: KMAP = DP
    integer(kind=I4B),                  intent(IN) :: nside
    real   (kind=KMAP),  dimension(0:), intent(inout) ::  map

    integer(kind=I4B) :: npix, ipn, ipr
    real   (kind=KMAP),  dimension(:), allocatable :: map_tmp
    !=======================================================================
    npix = 12*nside*nside
    allocate(map_tmp(0:npix-1))

!$OMP parallel default(none) &
!$OMP   shared(map, map_tmp, npix, nside) &
!$OMP   private(ipr, ipn)
!$OMP do schedule(dynamic,64)
    do ipr = 0, npix-1
       call ring2nest(nside, ipr, ipn)
       map_tmp(ipn) = map(ipr)
    enddo
!$OMP end do
!$OMP end parallel
    map = map_tmp

    deallocate(map_tmp)
    return
  end subroutine convert_ring2nest_double_1d

  !**************************************************************
  ! following 6 routines: out of place conversions for N-Dim maps
  ! if size(map) = Npix*Nd, peak memory = (Nd+2)*Npix
  ! 2004-08-25, EH
  !**************************************************************
  !=======================================================================
  subroutine convert_nest2ring_int_nd(nside, map)
    !=======================================================================
    !   NEST to RING conversion: 1D, integer
    !=======================================================================
    character(len=*),  parameter :: code = "convert_nest2ring_int_nd"
    integer(kind=I4B), parameter :: KMAP = I4B
    integer(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    integer(kind=KMAP), dimension(:),               allocatable :: map_tmp
    integer(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_nest2ring_nd_inc.f90'
  end subroutine convert_nest2ring_int_nd
  !=======================================================================
  subroutine convert_nest2ring_real_nd(nside, map)
    !=======================================================================
    !   NEST to RING conversion: 1D, real
    !=======================================================================
    character(len=*),  parameter :: code = "convert_nest2ring_real_nd"
    integer(kind=I4B), parameter :: KMAP = SP
    real(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    real(kind=KMAP), dimension(:),               allocatable :: map_tmp
    real(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_nest2ring_nd_inc.f90'
  end subroutine convert_nest2ring_real_nd
  !=======================================================================
  subroutine convert_nest2ring_double_nd(nside, map)
    !=======================================================================
    !   NEST to RING conversion: 1D, double
    !=======================================================================
    character(len=*),  parameter :: code = "convert_nest2ring_double_nd"
    integer(kind=I4B), parameter :: KMAP = DP
    real(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    real(kind=KMAP), dimension(:),               allocatable :: map_tmp
    real(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_nest2ring_nd_inc.f90'
  end subroutine convert_nest2ring_double_nd
  !=======================================================================
  subroutine convert_ring2nest_int_nd(nside, map)
    !=======================================================================
    !   RING to NEST conversion: 1D, integer
    !=======================================================================
    character(len=*),  parameter :: code = "convert_ring2nest_int_nd"
    integer(kind=I4B), parameter :: KMAP = I4B
    integer(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    integer(kind=KMAP), dimension(:),               allocatable :: map_tmp
    integer(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_ring2nest_nd_inc.f90'
  end subroutine convert_ring2nest_int_nd
  !=======================================================================
  subroutine convert_ring2nest_real_nd(nside, map)
    !=======================================================================
    !   RING to NEST conversion: 1D, real
    !=======================================================================
    character(len=*),  parameter :: code = "convert_ring2nest_real_nd"
    integer(kind=I4B), parameter :: KMAP = SP
    real(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    real(kind=KMAP), dimension(:),               allocatable :: map_tmp
    real(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_ring2nest_nd_inc.f90'
  end subroutine convert_ring2nest_real_nd
  !=======================================================================
  subroutine convert_ring2nest_double_nd(nside, map)
    !=======================================================================
    !   RING to NEST conversion: 1D, double
    !=======================================================================
    character(len=*),  parameter :: code = "convert_ring2nest_double_nd"
    integer(kind=I4B), parameter :: KMAP = DP
    real(kind=KMAP), dimension(0:,1:), intent(inout), target ::  map
    real(kind=KMAP), dimension(:),               allocatable :: map_tmp
    real(kind=KMAP), dimension(:),                   pointer :: map1
    include 'convert_ring2nest_nd_inc.f90'
  end subroutine convert_ring2nest_double_nd
  !====================================================================
  ! The following 6 routines convert in place integer, real, and double
  ! arrays between the NESTED and RING schemes.
  !
  ! in place: without allocating a temporary map. This routine is more general,
  ! but slower than convert_nest2ring.
  !
  !     This is a wrapper for the toolbox functions "ring2nest" and
  !     "nest2ring". Their names are supplied in the "subcall"
  !     argument.
  !
  ! Author: Benjamin D. Wandelt October 1997
  ! Added to pix_tools for version 1.00 in March 1999.
  ! 2004-08-25: EH, added N-Dim facility and double precision IO
  !====================================================================
  !====================================================================
  subroutine convert_inplace_int_1d(subcall,map)
    !==================================================================
    ! 1D, integer implementation
    !==================================================================
    character(len=*),  parameter :: code = "convert_inplace_int_1d"
    integer(kind=I4B), parameter :: KMAP = I4B
    integer(kind=KMAP), dimension(0:)    :: map
    integer(kind=KMAP)                   :: pixbuf1,pixbuf2
    include 'convert_inplace_1d_inc.f90'
  end subroutine convert_inplace_int_1d
  !====================================================================
  subroutine convert_inplace_real_1d(subcall,map)
    !====================================================================
    ! 1D, real implementation
    !==================================================================
    character(len=*),  parameter :: code = "convert_inplace_real_1d"
    integer(kind=I4B), parameter :: KMAP = SP
    real   (kind=KMAP), dimension(0:)    :: map
    real   (kind=KMAP)                   :: pixbuf1,pixbuf2
    include 'convert_inplace_1d_inc.f90'
  end subroutine convert_inplace_real_1d
  !====================================================================
  subroutine convert_inplace_double_1d(subcall,map)
    !====================================================================
    ! 1D, double precision implementation
    !==================================================================
    character(len=*),  parameter :: code = "convert_inplace_double_1d"
    integer(kind=I4B), parameter :: KMAP = DP
    real   (kind=KMAP), dimension(0:)    :: map
    real   (kind=KMAP)                   :: pixbuf1,pixbuf2
    include 'convert_inplace_1d_inc.f90'
  end subroutine convert_inplace_double_1d
  !====================================================================
  subroutine convert_inplace_int_nd(subcall,map)
    !==================================================================
    ! ND, integer implementation
    !==================================================================
    character(len=*),  parameter :: code = "convert_inplace_int_nd"
    integer(kind=i4b), parameter :: ND_MAX = 10
    integer(kind=I4B), parameter :: KMAP = I4B
    integer(kind=KMAP), dimension(0:,1:)    :: map
    integer(kind=KMAP), dimension(1:ND_MAX) :: pixbuf1,pixbuf2
    include 'convert_inplace_nd_inc.f90'
  end subroutine convert_inplace_int_nd
  !====================================================================
  subroutine convert_inplace_real_nd(subcall,map)
    !====================================================================
    ! ND, real implementation
    !==================================================================
    character(len=*),  parameter :: code = "convert_inplace_real_nd"
    integer(kind=i4b), parameter :: ND_MAX = 10
    integer(kind=I4B), parameter :: KMAP = SP
    real   (kind=KMAP), dimension(0:,1:)    :: map
    real   (kind=KMAP), dimension(1:ND_MAX) :: pixbuf1,pixbuf2
    include 'convert_inplace_nd_inc.f90'
  end subroutine convert_inplace_real_nd
  !====================================================================
  subroutine convert_inplace_double_nd(subcall,map)
    !====================================================================
    ! ND, double precision implementation
    !==================================================================
    character(len=*),  parameter :: code = "convert_inplace_double_nd"
    integer(kind=i4b), parameter :: ND_MAX = 10
    integer(kind=I4B), parameter :: KMAP = DP
    real   (kind=KMAP), dimension(0:,1:)    :: map
    real   (kind=KMAP), dimension(1:ND_MAX) :: pixbuf1,pixbuf2
    include 'convert_inplace_nd_inc.f90'
  end subroutine convert_inplace_double_nd
  !=======================================================================
  subroutine nest2ring(nside, ipnest, ipring)
    !=======================================================================
    !     performs conversion from NESTED to RING pixel number
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) ::  nside, ipnest
    INTEGER(KIND=I4B), INTENT(OUT) :: ipring
    INTEGER(KIND=I4B) ::  npix, npface, face_num, ncap, n_before, &
         &     ipf, ip_low, ip_trunc, ip_med, ip_hi, &
         &     ix, iy, jrt, jr, nr, jpt, jp, kshift, nl4

    ! coordinate of the lowest corner of each face
    INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    npix = 12 * nside**2
    if (ipnest<0 .or. ipnest>npix-1) call fatal_error("ipnest out of range")

    !     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    ncap  = 2*nside*(nside-1) ! number of points in the North Polar cap
    nl4   = 4*nside

    !     finds the face, and the number in the face
    npface = nside**2
    !ccccc      ip = ipnest - 1         ! in {0,npix-1}

    face_num = ipnest/npface  ! face number in {0,11}
    ipf = MODULO(ipnest,npface)  ! pixel number in the face {0,npface-1}

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    ip_low = MODULO(ipf,1024)       ! content of the last 10 bits
    ip_trunc =   ipf/1024        ! truncation of the last 10 bits
    ip_med = MODULO(ip_trunc,1024)  ! content of the next 10 bits
    ip_hi  =     ip_trunc/1024   ! content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}

    nr = nside                  ! equatorial region (the most frequent)
    n_before = ncap + nl4 * (jr - nside)
    kshift = MODULO(jr - nside, 2)
    if (jr < nside) then     ! north pole region
       nr = jr
       n_before = 2 * nr * (nr - 1)
       kshift = 0
    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       n_before = npix - 2 * (nr + 1) * nr
       kshift = 0
    endif

    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}

    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4

    ipring = n_before + jp - 1 ! in {0, npix-1}
    return
  end subroutine nest2ring
  !=======================================================================
  subroutine ring2nest(nside, ipring, ipnest)
    !=======================================================================
    !     performs conversion from RING to NESTED pixel number
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside, ipring
    INTEGER(KIND=I4B), INTENT(OUT) :: ipnest

    REAL(KIND=DP) :: fihip, hip
    INTEGER(KIND=I4B) :: npix, nl2, nl4, ncap, ip, iphi, ipt, ipring1, &
         &     kshift, face_num, nr, &
         &     irn, ire, irm, irs, irt, ifm , ifp, &
         &     ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf

    ! coordinate of the lowest corner of each face
    INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    npix = 12*nside**2      ! total number of points
    if (ipring <0 .or. ipring>npix-1) call fatal_error("ipring out of range")
    if (x2pix(128) <= 0) call mk_xy2pix()

    nl2 = 2*nside
    nl4 = 4*nside
    ncap = nl2*(nside-1) ! points in each polar cap, =0 for nside =1
    ipring1 = ipring + 1

    !     finds the ring number, the position of the ring and the face number
    if (ipring1 <= ncap) then ! north polar cap

       hip   = ipring1/2.0_dp
       fihip = AINT ( hip ,kind=DP)
       irn   = INT( SQRT( hip - SQRT(fihip) ) ) + 1 ! counted from North pole
       iphi  = ipring1 - 2*irn*(irn - 1)

       kshift = 0
       nr = irn                  ! 1/4 of the number of points on the current ring
       face_num = (iphi-1) / irn ! in {0,3}

    elseif (ipring1 <= nl2*(5*nside+1)) then ! equatorial region

       ip    = ipring1 - ncap - 1
       irn   = INT( ip / nl4 ) + nside               ! counted from North pole
       iphi  = MODULO(ip,nl4) + 1

       kshift  = MODULO(irn+nside,2)  ! 1 if irn+nside is odd, 0 otherwise
       nr = nside
       ire =  irn - nside + 1 ! in {1, 2*nside +1}
       irm =  nl2 + 2 - ire
       ifm = (iphi - ire/2 + nside -1) / nside ! face boundary
       ifp = (iphi - irm/2 + nside -1) / nside
       if (ifp == ifm) then          ! faces 4 to 7
          face_num = MODULO(ifp,4) + 4
       else if (ifp + 1 == ifm) then ! (half-)faces 0 to 3
          face_num = ifp
       else if (ifp - 1 == ifm) then ! (half-)faces 8 to 11
          face_num = ifp + 7
       endif

    else ! south polar cap

       ip    = npix - ipring1 + 1
       hip   = ip/2.0_dp
       fihip = AINT ( hip ,kind=DP)
       irs   = INT( SQRT( hip - SQRT(fihip) ) ) + 1  ! counted from South pole
       iphi  = 4*irs + 1 - (ip - 2*irs*(irs-1))

       kshift = 0
       nr = irs
       irn   = nl4 - irs
       face_num = (iphi-1) / irs + 8 ! in {8,11}

    endif

    !     finds the (x,y) on the face
    irt =   irn  - jrll(face_num+1)*nside + 1       ! in {-nside+1,0}
    ipt = 2*iphi - jpll(face_num+1)*nr - kshift - 1 ! in {-nside+1,nside-1}
    if (ipt >= nl2) ipt = ipt - 8*nside ! for the face #4

    ix =  (ipt - irt ) / 2
    iy = -(ipt + irt ) / 2

    ix_low = MODULO(ix,128)
    ix_hi  = ix/128
    iy_low = MODULO(iy,128)
    iy_hi  = iy/128

    ipf =  (x2pix(ix_hi +1)+y2pix(iy_hi +1)) * (128 * 128) &
         &     + (x2pix(ix_low+1)+y2pix(iy_low+1))        ! in {0, nside**2 - 1}


    ipnest = ipf + face_num* nside **2    ! in {0, 12*nside**2 - 1}

    return
  end subroutine ring2nest
  !=======================================================================
  subroutine xy2pix_nest(nside, ix, iy, face_num, ipix)
    !=======================================================================
    !     gives the pixel number ipix (NESTED)
    !     corresponding to ix, iy and face_num
    !
    !     Benjamin D. Wandelt 13/10/97
    !     using code from HEALPIX toolkit by K.Gorski and E. Hivon
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) ::  nside, ix, iy, face_num
    INTEGER(KIND=I4B), INTENT(OUT) :: ipix
    INTEGER(KIND=I4B) ::  ix_low, ix_hi, iy_low, iy_hi, ipf

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    if (ix<0 .or. ix>(nside-1)) call fatal_error("ix out of range")
    if (iy<0 .or. iy>(nside-1)) call fatal_error("iy out of range")
    if (x2pix(128) <= 0) call mk_xy2pix()

    ix_low = MODULO(ix,128)
    ix_hi  =     ix/128
    iy_low = MODULO(iy,128)
    iy_hi  =     iy/128

    ipf =  (x2pix(ix_hi +1)+y2pix(iy_hi +1)) * (128 * 128) &
         &     + (x2pix(ix_low+1)+y2pix(iy_low+1))

    ipix = ipf + face_num* nside **2    ! in {0, 12*nside**2 - 1}
    return
  end subroutine xy2pix_nest
  !=======================================================================
  subroutine pix2xy_nest(nside, ipf, ix, iy)
    !=======================================================================
    !     gives the x, y coords in a face from pixel number within the face (NESTED)
    !
    !     Benjamin D. Wandelt 13/10/97
    !
    !     using code from HEALPIX toolkit by K.Gorski and E. Hivon
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: nside, ipf
    INTEGER(KIND=I4B), INTENT(OUT) :: ix, iy

    INTEGER(KIND=I4B) ::  ip_low, ip_trunc, ip_med, ip_hi

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    if (ipf <0 .or. ipf>nside*nside-1) &
         &     call fatal_error("ipix out of range")
    if (pix2x(1023) <= 0) call mk_pix2xy()

    ip_low = MODULO(ipf,1024)       ! content of the last 10 bits
    ip_trunc =   ipf/1024        ! truncation of the last 10 bits
    ip_med = MODULO(ip_trunc,1024)  ! content of the next 10 bits
    ip_hi  =     ip_trunc/1024   ! content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
    return
  end subroutine pix2xy_nest
  !=======================================================================
  subroutine mk_pix2xy()
    !=======================================================================
    !     constructs the array giving x and y in the face from pixel number
    !     for the nested (quad-cube like) ordering of pixels
    !
    !     the bits corresponding to x and y are interleaved in the pixel number
    !     one breaks up the pixel number by even and odd bits
    !=======================================================================
    INTEGER(KIND=I4B) ::  kpix, jpix, ix, iy, ip, id

    !cc cf block data      data      pix2x(1023) /0/
    !-----------------------------------------------------------------------
    !      print *, 'initiate pix2xy'
    do kpix=0,1023          ! pixel number
       jpix = kpix
       IX = 0
       IY = 0
       IP = 1               ! bit position (in x and y)
!        do while (jpix/=0) ! go through all the bits
       do
          if (jpix == 0) exit ! go through all the bits
          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in ix
          jpix = jpix/2
          IX = ID*IP+IX

          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in iy
          jpix = jpix/2
          IY = ID*IP+IY

          IP = 2*IP         ! next bit (in x and y)
       enddo
       pix2x(kpix) = IX     ! in 0,31
       pix2y(kpix) = IY     ! in 0,31
    enddo

    return
  end subroutine mk_pix2xy
  !=======================================================================
  subroutine mk_xy2pix()
    !=======================================================================
    !     sets the array giving the number of the pixel lying in (x,y)
    !     x and y are in {1,128}
    !     the pixel number is in {0,128**2-1}
    !
    !     if  i-1 = sum_p=0  b_p * 2^p
    !     then ix = sum_p=0  b_p * 4^p
    !          iy = 2*ix
    !     ix + iy in {0, 128**2 -1}
    !=======================================================================
    INTEGER(KIND=I4B):: k,ip,i,j,id
    !=======================================================================

    do i = 1,128           !for converting x,y into
       j  = i-1            !pixel numbers
       k  = 0
       ip = 1

       do
          if (j==0) then
             x2pix(i) = k
             y2pix(i) = 2*k
             exit
          else
             id = MODULO(J,2)
             j  = j/2
             k  = ip*id+k
             ip = ip*4
          endif
       enddo

    enddo

    RETURN
  END subroutine mk_xy2pix
  !=======================================================================

  ! The following is a routine which finds the 7 or 8 neighbours of
  ! any pixel in the nested scheme of the HEALPIX pixelisation.
  !====================================================================
  subroutine neighbours_nest(nside,ipix,n,nneigh)
    !====================================================================
    !   Returns list n(8) of neighbours of pixel ipix (in NESTED scheme)
    !   the neighbours are ordered in the following way:
    !   First pixel is the one to the south (the one west of the south
    ! direction is taken
    ! for the pixels which don't have a southern neighbour). From
    ! then on the neighbours are ordered in the clockwise direction
    ! about the pixel with number ipix.
    !
    !   nneigh is the number of neighbours (mostly 8, 8 pixels have 7 neighbours)
    !
    !   Benjamin D. Wandelt October 1997
    !   Added to pix_tools in March 1999
    !   added 'return' for case nside=1, EH, Oct 2005
    !====================================================================
    use bit_manipulation
    integer(kind=i4b), intent(in)::nside, ipix
    integer(kind=i4b), intent(out), dimension(1:):: n
    integer(kind=i4b), intent(out):: nneigh

    integer(kind=i4b) :: npix,ipf,ipo,ix,ixm,ixp,iy,iym,iyp,ixo,iyo
    integer(kind=i4b) :: face_num,other_face
    integer(kind=i4b) :: ia,ib,ibp,ibm,ib2,icase,nsidesq
    integer(kind=i4b) :: local_magic1,local_magic2

!     integer(kind=i4b), intrinsic :: IAND

    !--------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    nsidesq=nside*nside
    npix = 12*nsidesq       ! total number of points
    if (ipix <0 .or. ipix>npix-1) call fatal_error("ipix out of range")

    if (nside == 1) then
       nneigh = 6
       if (ipix==0 ) n(1:6) = (/ 8, 4, 3, 2, 1, 5 /)
       if (ipix==1 ) n(1:6) = (/ 9, 5, 0, 3, 2, 6 /)
       if (ipix==2 ) n(1:6) = (/10, 6, 1, 0, 3, 7 /)
       if (ipix==3 ) n(1:6) = (/11, 7, 2, 1, 0, 4 /)
       if (ipix==4 ) n(1:6) = (/11, 7, 3, 0, 5, 8 /)
       if (ipix==5 ) n(1:6) = (/ 8, 4, 0, 1, 6, 9 /)
       if (ipix==6 ) n(1:6) = (/ 9, 5, 1, 2, 7,10 /)
       if (ipix==7 ) n(1:6) = (/10, 6, 2, 3, 8,11 /)
       if (ipix==8 ) n(1:6) = (/10,11, 4, 0, 5, 9 /)
       if (ipix==9 ) n(1:6) = (/11, 8, 5, 1, 5,10 /)
       if (ipix==10) n(1:6) = (/ 8, 9, 6, 2, 7,11 /)
       if (ipix==11) n(1:6) = (/ 9,10, 6, 3, 4, 8 /)
       return
    endif

    !     initiates array for (x,y)-> pixel number -> (x,y) mapping
    if (x2pix(128) <= 0) call mk_xy2pix()

    local_magic1=(nsidesq-1)/3
    local_magic2=2*local_magic1
    face_num=ipix/nsidesq

    ipf=modulo(ipix,nsidesq)   !Pixel number in face

    call pix2xy_nest(nside,ipf,ix,iy)
    ixm=ix-1
    ixp=ix+1
    iym=iy-1
    iyp=iy+1

    nneigh=8                  !Except in special cases below

    !     Exclude corners
    if(ipf==local_magic2)     then !WestCorner
       icase=5
       goto 100
    endif
    if(ipf==(nsidesq-1)) then !NorthCorner
       icase=6
       goto 100
    endif
    if(ipf==0)           then !SouthCorner
       icase=7
       goto 100
    endif
    if(ipf==local_magic1)     then !EastCorner
       icase=8
       goto 100
    endif

    !     Detect edges
    if(IAND(ipf,local_magic1)==local_magic1) then !NorthEast
       icase=1
       goto 100
    endif
    if(IAND(ipf,local_magic1)==0)      then !SouthWest
       icase=2
       goto 100
    endif
    if(IAND(ipf,local_magic2)==local_magic2) then !NorthWest
       icase=3
       goto 100
    endif
    if(IAND(ipf,local_magic2)==0)      then !SouthEast
       icase=4
       goto 100
    endif

    !     Inside a face
    call xy2pix_nest(nside, ixm, iym, face_num, n(1))
    call xy2pix_nest(nside, ixm, iy , face_num, n(2))
    call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
    call xy2pix_nest(nside, ix , iyp, face_num, n(4))
    call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
    call xy2pix_nest(nside, ixp, iy , face_num, n(6))
    call xy2pix_nest(nside, ixp, iym, face_num, n(7))
    call xy2pix_nest(nside, ix , iym, face_num, n(8))
    return

100 continue

    ia= face_num/4            !in {0,2}
    ib= modulo(face_num,4)       !in {0,3}
    ibp=modulo(ib+1,4)
    ibm=modulo(ib+4-1,4)
    ib2=modulo(ib+2,4)

    if(ia==0) then          !North Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ibp
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          ipo=modulo(swapLSBMSB(ipf),nsidesq)    !East-West flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1 , iyo, other_face, n(5))
          n(6)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(7))
       case(2)              !SouthWest edge
          other_face=4+ib
          ipo=modulo(invLSB(ipf),nsidesq)        !SW-NE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
          n(2)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(3))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       case(3)              !NorthWest edge
          other_face=0+ibm
          ipo=modulo(swapLSBMSB(ipf),nsidesq)    !East-West flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(3))
          n(4)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(5))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       case(4)              !SouthEast edge
          other_face=4+ibp
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(7))
          n(8)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
       case(5)              !West corner
          nneigh=7
          other_face=4+ib
          n(2)=other_face*nsidesq+nsidesq-1
          n(1)=n(2)-2
          other_face=0+ibm
          n(3)=other_face*nsidesq+local_magic1
          n(4)=n(3)+2
          n(5)=ipix+1
          n(6)=ipix-1
          n(7)=ipix-2
       case(6)              !North corner
          n(1)=ipix-3
          n(2)=ipix-1
          n(8)=ipix-2
          other_face=0+ibm
          n(4)=other_face*nsidesq+nsidesq-1
          n(3)=n(4)-2
          other_face=0+ib2
          n(5)=other_face*nsidesq+nsidesq-1
          other_face=0+ibp
          n(6)=other_face*nsidesq+nsidesq-1
          n(7)=n(6)-1
       case(7)              !South corner
          other_face=8+ib
          n(1)=other_face*nsidesq+nsidesq-1
          other_face=4+ib
          n(2)=other_face*nsidesq+local_magic1
          n(3)=n(2)+2
          n(4)=ipix+2
          n(5)=ipix+3
          n(6)=ipix+1
          other_face=4+ibp
          n(8)=other_face*nsidesq+local_magic2
          n(7)=n(8)+1
       case(8)              !East corner
          nneigh=7
          n(2)=ipix-1
          n(3)=ipix+1
          n(4)=ipix+2
          other_face=0+ibp
          n(6)=other_face*nsidesq+local_magic2
          n(5)=n(6)+1
          other_face=4+ibp
          n(7)=other_face*nsidesq+nsidesq-1
          n(1)=n(7)-1
       end select ! north

    elseif(ia==1) then      !Equatorial region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ib
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo , iyo+1, other_face, n(5))
          n(6)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(7))
       case(2)              !SouthWest edge
          other_face=8+ibm
          ipo=modulo(invLSB(ipf),nsidesq)        !SW-NE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
          n(2)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(3))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       case(3)              !NorthWest edge
          other_face=0+ibm
          ipo=modulo(invMSB(ipf),nsidesq)    !NW-SE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(3))
          n(4)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(5))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       case(4)              !SouthEast edge
          other_face=8+ib
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(7))
          n(8)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
       case(5)              !West corner
          other_face=8+ibm
          n(2)=other_face*nsidesq+nsidesq-1
          n(1)=n(2)-2
          other_face=4+ibm
          n(3)=other_face*nsidesq+local_magic1
          other_face=0+ibm
          n(4)=other_face*nsidesq
          n(5)=n(4)+1
          n(6)=ipix+1
          n(7)=ipix-1
          n(8)=ipix-2
       case(6)              !North corner
          nneigh=7
          n(1)=ipix-3
          n(2)=ipix-1
          other_face=0+ibm
          n(4)=other_face*nsidesq+local_magic1
          n(3)=n(4)-1
          other_face=0+ib
          n(5)=other_face*nsidesq+local_magic2
          n(6)=n(5)-2
          n(7)=ipix-2
       case(7)              !South corner
          nneigh=7
          other_face=8+ibm
          n(1)=other_face*nsidesq+local_magic1
          n(2)=n(1)+2
          n(3)=ipix+2
          n(4)=ipix+3
          n(5)=ipix+1
          other_face=8+ib
          n(7)=other_face*nsidesq+local_magic2
          n(6)=n(7)+1
       case(8)              !East corner
          other_face=8+ib
          n(8)=other_face*nsidesq+nsidesq-1
          n(1)=n(8)-1
          n(2)=ipix-1
          n(3)=ipix+1
          n(4)=ipix+2
          other_face=0+ib
          n(6)=other_face*nsidesq
          n(5)=n(6)+2
          other_face=4+ibp
          n(7)=other_face*nsidesq+local_magic2
       end select ! equator
    else                    !South Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=4+ibp
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo , iyo+1, other_face, n(5))
          n(6)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(7))
       case(2)              !SouthWest edge
          other_face=8+ibm
          ipo=modulo(swapLSBMSB(ipf),nsidesq)        !W-E flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
          n(2)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(3))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       case(3)              !NorthWest edge
          other_face=4+ib
          ipo=modulo(invMSB(ipf),nsidesq)    !NW-SE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(3))
          n(4)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(5))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       case(4)              !SouthEast edge
          other_face=8+ibp
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          ipo=modulo(swapLSBMSB(ipf),nsidesq) !E-W flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(7))
          n(8)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
       case(5)              !West corner
          nneigh=7
          other_face=8+ibm
          n(2)=other_face*nsidesq+local_magic1
          n(1)=n(2)-1
          other_face=4+ib
          n(3)=other_face*nsidesq
          n(4)=n(3)+1
          n(5)=ipix+1
          n(6)=ipix-1
          n(7)=ipix-2
       case(6)              !North corner
          n(1)=ipix-3
          n(2)=ipix-1
          other_face=4+ib
          n(4)=other_face*nsidesq+local_magic1
          n(3)=n(4)-1
          other_face=0+ib
          n(5)=other_face*nsidesq
          other_face=4+ibp
          n(6)=other_face*nsidesq+local_magic2
          n(7)=n(6)-2
          n(8)=ipix-2
       case(7)              !South corner
          other_face=8+ib2
          n(1)=other_face*nsidesq
          other_face=8+ibm
          n(2)=other_face*nsidesq
          n(3)=n(2)+1
          n(4)=ipix+2
          n(5)=ipix+3
          n(6)=ipix+1
          other_face=8+ibp
          n(8)=other_face*nsidesq
          n(7)=n(8)+2
       case(8)              !East corner
          nneigh=7
          other_face=8+ibp
          n(7)=other_face*nsidesq+local_magic2
          n(1)=n(7)-2
          n(2)=ipix-1
          n(3)=ipix+1
          n(4)=ipix+2
          other_face=4+ibp
          n(6)=other_face*nsidesq
          n(5)=n(6)+2
       end select ! south
    endif

    return
  end subroutine neighbours_nest
  !====================================================================
  subroutine next_in_line_nest(nside, ipix, inext)
    !====================================================================
    !   given nside and a NESTED pixel number ipix, returns in inext
    !  the pixel that lies on the East side (and the same latitude) as ipix
    !
    !   Hacked by EH from BDW's neighbours_nest, 2001-12-18
    !   Hacked for Nside=1 by EH, 2004-05-28
    !====================================================================
    use bit_manipulation
    integer(kind=i4b), intent(in)::nside, ipix
    integer(kind=i4b), intent(out):: inext

    integer(kind=i4b) :: npix,ipf,ipo,ix,ixp,iy,iym,ixo,iyo
    integer(kind=i4b) :: face_num,other_face
    integer(kind=i4b) :: ia,ib,ibp,ibm,ib2,icase,nsidesq
    integer(kind=i4b) :: local_magic1,local_magic2

    !--------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    nsidesq=nside*nside
    npix = 12*nsidesq       ! total number of points
    if (ipix <0 .or. ipix>npix-1) call fatal_error("ipix out of range")

    ! quick and dirty hack for Nside=1
    if (nside == 1) then
       inext = ipix + 1
       if (ipix == 3)  inext = 0
       if (ipix == 7)  inext = 4
       if (ipix == 11) inext = 8
       return
    endif
    !     initiates array for (x,y)-> pixel number -> (x,y) mapping
    if (x2pix(128) <= 0) call mk_xy2pix()

    local_magic1=(nsidesq-1)/3
    local_magic2=2*local_magic1
    face_num=ipix/nsidesq

    ipf=modulo(ipix,nsidesq)   !Pixel number in face

    call pix2xy_nest(nside,ipf,ix,iy)
    ixp=ix+1
    iym=iy-1

    !     Exclude corners
    if(ipf==local_magic2)     then !WestCorner
       inext = ipix - 1
       return
    endif
    if(ipf==(nsidesq-1)) then !NorthCorner
       icase=6
       goto 100
    endif
    if(ipf==0)           then !SouthCorner
       icase=7
       goto 100
    endif
    if(ipf==local_magic1)     then !EastCorner
       icase=8
       goto 100
    endif

    !     Detect edges
    if(IAND(ipf,local_magic1)==local_magic1) then !NorthEast
       icase=1
       goto 100
    endif
    if(IAND(ipf,local_magic2)==0)      then !SouthEast
       icase=4
       goto 100
    endif

    !     Inside a face
    call xy2pix_nest(nside, ixp, iym, face_num, inext)
    return

100 continue

    ia= face_num/4            !in {0,2}
    ib= modulo(face_num,4)       !in {0,3}
    ibp=modulo(ib+1,4)
    ibm=modulo(ib+4-1,4)
    ib2=modulo(ib+2,4)

    if(ia==0) then          !North Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ibp
          ipo=modulo(swapLSBMSB(ipf),nsidesq)    !East-West flip
          inext = other_face*nsidesq+ipo         ! (6)
       case(4)              !SouthEast edge
          other_face=4+ibp
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, inext)
       case(6)              !North corner
          other_face=0+ibp
          inext=other_face*nsidesq+nsidesq-1
       case(7)              !South corner
          other_face=4+ibp
          inext=other_face*nsidesq+local_magic2+1
       case(8)              !East corner
          other_face=0+ibp
          inext=other_face*nsidesq+local_magic2
       end select ! north

    elseif(ia==1) then      !Equatorial region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ib
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, inext)
       case(4)              !SouthEast edge
          other_face=8+ib
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, inext)
       case(6)              !North corner
          other_face=0+ib
          inext=other_face*nsidesq+local_magic2-2
       case(7)              !South corner
          other_face=8+ib
          inext=other_face*nsidesq+local_magic2+1
       case(8)              !East corner
          other_face=4+ibp
          inext=other_face*nsidesq+local_magic2
       end select ! equator
    else                    !South Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=4+ibp
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, inext)
       case(4)              !SouthEast edge
          other_face=8+ibp
          ipo=modulo(swapLSBMSB(ipf),nsidesq) !E-W flip
          inext = other_face*nsidesq+ipo   ! (8)
       case(6)              !North corner
          other_face=4+ibp
          inext=other_face*nsidesq+local_magic2 -2
       case(7)              !South corner
          other_face=8+ibp
          inext=other_face*nsidesq
       case(8)              !East corner
          other_face=8+ibp
          inext=other_face*nsidesq+local_magic2
       end select ! south
    endif

    return
  end subroutine next_in_line_nest
  !=======================================================================
  subroutine ang2vec(theta, phi, vector)
    !=======================================================================
    !     renders the vector (x,y,z) corresponding to angles
    !     theta (co-latitude measured from North pole, in [0,Pi] radians)
    !     and phi (longitude measured eastward, in radians)
    !     North pole is (x,y,z)=(0,0,1)
    !     added by EH, Feb 2000
    !=======================================================================
    REAL(KIND=DP), INTENT(IN) :: theta, phi
    REAL(KIND=DP), INTENT(OUT), dimension(1:) :: vector

    REAL(KIND=DP) :: sintheta
    !=======================================================================

    if (theta<0.0_dp .or. theta>pi)  then
       print*,"ANG2VEC: theta : ",theta," is out of range [0, Pi]"
       call fatal_error
    endif
    sintheta = SIN(theta)

    vector(1) = sintheta * COS(phi)
    vector(2) = sintheta * SIN(phi)
    vector(3) = COS(theta)

    return
  end subroutine ang2vec
  !=======================================================================
  subroutine vec2ang(vector, theta, phi)
    !=======================================================================
    !     renders the angles theta, phi corresponding to vector (x,y,z)
    !     theta (co-latitude measured from North pole, in [0,Pi] radians)
    !     and phi (longitude measured eastward, in [0,2Pi[ radians)
    !     North pole is (x,y,z)=(0,0,1)
    !     added by EH, Feb 2000
    !=======================================================================
    REAL(KIND=DP), INTENT(IN), dimension(1:) :: vector
    REAL(KIND=DP), INTENT(OUT) :: theta, phi

    REAL(KIND=DP) :: dnorm, z
    !=======================================================================

    dnorm = SQRT(vector(1)**2+vector(2)**2+vector(3)**2)

    z = vector(3) / dnorm
    theta = ACOS(z)

    phi = 0.0_dp
    if (vector(1) /= 0.0_dp .or. vector(2) /= 0.0_dp) &
         &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]
    if (phi < 0.0)     phi = phi + twopi ! phi in [0,2pi[

    return
  end subroutine vec2ang
  !=======================================================================
  function npix2nside(npix) result(nside_result)
    !=======================================================================
    ! given npix, returns nside such that npix = 12*nside^2
    !  nside should be a power of 2 smaller than 8192
    !  if not, -1 is returned
    ! EH, Feb-2000
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: npix
!     INTEGER(KIND=I4B) :: npix2nside
    INTEGER(KIND=I4B) :: nside_result
    INTEGER(KIND=I4B) :: nside, ilog
    REAL(KIND=DP) :: fnside, fnpix, flog
    CHARACTER(LEN=*), PARAMETER :: code = "npix2nside"
    INTEGER(KIND=I4B), PARAMETER :: npix_max = 12*ns_max*ns_max
    !=======================================================================

    fnside = sqrt(npix/12.0_dp)
    nside = NINT(fnside)

    if (npix < 12) then
       print*, code//": Npix is too small :",npix
       print*, "                       < 12 "
       nside = -1
       goto 1
    endif

    if (npix > npix_max) then
       print*, code//": Npix is too large :",npix
       print*, "                       > ",npix_max
       nside = -1
       goto 1
    endif

    fnpix = 12.0_dp*nside*nside
    if (abs(fnpix-npix) > 1.0e-2) then
       print*, code//": wrong Npix ",npix
       print*, "    not 12*nside*nside "
       nside = -1
       goto 1
    endif

    flog = log(real(nside,kind=dp))/log(2.0_dp)
    ilog = NINT(flog)
    if (abs(ilog-flog) > 1.0e-6_dp) then
       print*, code//": wrong Nside ",nside
       print*, "    not a power of 2 "
       nside = -1
       goto 1
    endif

1   continue
!     npix2nside = nside
    nside_result = nside
    return

  end function npix2nside
  !=======================================================================
  function nside2npix(nside) result(npix_result)
    !=======================================================================
    ! given nside, returns npix such that npix = 12*nside^2
    !  nside should be a power of 2 smaller than 8192
    !  if not, -1 is returned
    ! EH, Feb-2000
    !=======================================================================
    INTEGER(KIND=I4B) :: npix_result
    INTEGER(KIND=I4B), INTENT(IN) :: nside

    INTEGER(KIND=I4B), dimension(1:14) :: listnside
    INTEGER(KIND=I4B) :: npix
    CHARACTER(LEN=*), PARAMETER :: code = "nside2npix"
    !=======================================================================

    listnside = (/1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192/)
    npix = 12*nside*nside

    if (MINVAL(abs(listnside-nside)) > 0) then
       print*,code//": invalid nside ",nside
!        write(unit=*,fmt="(a,4(i2),3(i3),3(i4),4(i5))") " Nside should be among ",listnside
       print "(a,4(i2),3(i3),3(i4),4(i5))", " Nside should be among ",listnside
       npix = -1
    endif

    npix_result = npix

    return
  end function nside2npix
  !=======================================================================
  subroutine surface_triangle(vec1, vec2, vec3, surface)
    !=======================================================================
    ! returns the surface in steradians
    !  of the spherical triangle with vertices vec1, vec2, vec3
    !
    ! algorithm : finds triangle sides and uses l'Huilier formula to compute
    ! "spherical excess" = surface area of triangle on a sphere of radius one
    ! see, eg Bronshtein, Semendyayev Eq 2.86
    !=======================================================================
    real(kind=dp), dimension(1:), intent(in) :: vec1, vec2, vec3
    real(kind=dp), intent(out) :: surface

    !   real(kind=dp), dimension(1:3) :: v1, v2, v3
    real(kind=dp), dimension(1:3) :: side
    !  real(kind=dp) :: hp
    real(kind=dp) :: x0, x1, x2, x3
    !=======================================================================

    ! half perimeter
!   hp = 0.5_dp * (side1 + side2 + side3)

!   ! l'Huilier formula
!   x0 = tan( hp          * 0.5_dp)
!   x1 = tan((hp - side1) * 0.5_dp)
!   x2 = tan((hp - side2) * 0.5_dp)
!   x3 = tan((hp - side3) * 0.5_dp)

    ! find triangle sides
    call angdist(vec2, vec3, side(1))
    call angdist(vec3, vec1, side(2))
    call angdist(vec1, vec2, side(3))
    ! divide by 4
    side(1:3) = side(1:3) * 0.25_dp

    ! l'Huilier formula
    x0 = tan( side(1) + side(2) + side(3) )
    x1 = tan(-side(1) + side(2) + side(3) )
    x2 = tan( side(1) - side(2) + side(3) )
    x3 = tan( side(1) + side(2) - side(3) )
    surface = 4.0_dp * atan( sqrt(x0 * x1 * x2 * x3) )

    return
  end subroutine surface_triangle
  !=======================================================================
  subroutine angdist(v1, v2, dist)
    !=======================================================================
    ! call angdist(v1, v2, dist)
    ! computes the angular distance dist (in rad) between 2 vectors v1 and v2
    ! in general dist = acos ( v1 . v2 )
    ! except if the 2 vectors are almost aligned.
    !=======================================================================
    real(kind=DP), intent(IN), dimension(1:) :: v1, v2
    real(kind=DP), intent(OUT) :: dist

    real(kind=DP), dimension(1:3) :: r1, r2, vdiff
    real(kind=DP) :: diff, sprod
    !=======================================================================

    ! normalize both vectors
    r1(1:3) = v1(1:3) / sqrt(dot_product(v1,v1))
    r2(1:3) = v2(1:3) / sqrt(dot_product(v2,v2))

    sprod = DOT_PRODUCT(r1, r2)

    if (sprod > 0.999_dp) then
       ! almost colinear vectors
       vdiff(1:3) = r1(1:3) - r2(1:3)
       diff = sqrt(dot_product(vdiff,vdiff)) ! norm of difference
       dist = 2.0_dp * asin(diff * 0.5_dp)

    else if (sprod < -0.999_dp) then
       ! almost anti-colinear vectors
       vdiff(1:3) = r1(1:3) + r2(1:3)
       diff = sqrt(dot_product(vdiff,vdiff)) ! norm of sum
       dist = PI - 2.0_dp * asin(diff * 0.5_dp)

    else
       ! other cases
       dist = acos( sprod )
    endif


    return
  end subroutine angdist
  !=======================================================================
  subroutine vect_prod(v1, v2, v3) !
    !=======================================================================
    !     returns in v3 the vectorial product of the 2 3-dimensional
    !     vectors v1 and v2
    !=======================================================================
    real(kind=DP), dimension(1:), INTENT(IN)  :: v1, v2
    real(kind=DP), dimension(1:), INTENT(OUT) :: v3
    !=======================================================================

    v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
    v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
    v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

    return
  end subroutine vect_prod

  !****************************************************************************

  !=======================================================================
  function nside2ntemplates(nside) result(ntemplates)
    !=======================================================================

    integer(kind=i4b) :: ntemplates
    integer(kind=i4b), intent(IN) :: nside

!    ntemplates = max( (nside*(nside+6))/4, 2)
    ntemplates = 1 + nside * (nside + 6)
    ntemplates = ntemplates / 4

    return
  end function nside2ntemplates

    !=======================================================================
  subroutine template_pixel_ring(nside, pixel, templ8, reflexion)
    !=======================================================================
    !    Returns template-pixel index corresponding to each (RING ordered) pixel provided
    !    The reflexion (E-W, or N-W or both) to apply to the pixel to match the
    !    template is optionally returned
    !-------------------------------------------------------------
    !   For a given resolution Nside, there are 12*Nside*Nside pixels on the
    !   sphere, whose shape belong to one of the (Nside+6)*Nside/4 templates.
    !   (2 in the case Nside=1).
    !   Any pixel can be matched to a single of these templates by a combination of
    !    - a rotation around the polar axis with 
    !    - reflexion(s) around a meridian and/or the equator.
    !   The reflexion number returned by this routine for each pixel is
    !   Reflexion         Operation to do to match template
    !     0                 rotation only
    !     1                 rotation + East-West swap
    !     2                 rotation + North-South swap
    !     3                 rotation + East-West + North-South swaps
    !   
    !   The template pixels are all located in the Northern Hemisphere, or on the
    !   Equator.
    !   They are chosen to have
    !       z >= 2/3,      0< phi <= Pi/2               [Nside*(Nside+2)/4]
    !       2/3 > z >= 0,  phi=0, or phi=Pi/(4*Nside)   [Nside]
    !   They are numbered from 0, starting at the North Pole, with the index
    !   increasing in phi, and then increasing for decreasing z.
    !-------------------------------------------------------------
    character(len=*),           parameter :: code = "Template_pixel_ring"
    integer(I4B),    intent(IN)           :: nside
    integer(I4B),    intent(IN)           :: pixel
    integer(I4B),    intent(OUT)          :: templ8
    integer(I4B),    intent(OUT), optional:: reflexion
    ! local variable
    integer(I4B) :: npix, ncap, ns4, n0
    integer(I4B) :: ip, iring, ifi, ifd, iff
    integer(I4B) :: refl
    real(DP)     :: fip
    !-------------------------------------------------------------

    npix = nside2npix(nside)
    if (npix < 0) then
       print*,'Invalid Nside = ',nside
       call fatal_error(code//' Abort')
    endif

    if (pixel < 0 .or. pixel >= npix) then
       print*,'Invalid Pixel = ',pixel
       call fatal_error(code//' Abort')
    endif

    ncap = 2*nside*(nside+1)
    ns4  = 4*nside
!     n0   = max( (nside*(nside+2))/4, 1)
    n0   = (1 + nside * (nside + 2)) / 4


    refl = 0 ! no swap

    ! North polar cap
    if (pixel < ncap) then
       ip    = pixel + 1
       fip   = real(int(ip/2) , kind=dp)
       iring = int(sqrt(ip/2.0_dp - sqrt(fip) )) + 1 ! counted from NORTH pole, starting at 1
       ifi   = MOD( (ip - 1 - 2*iring*(iring-1)), iring) ! in [0,iring-1], upwards
       ifd   = iring - 1 - ifi                           ! in [0,iring-1], downwards
       iff   = min(ifi, ifd)                             ! in [0,(iring-1)/2]
       templ8 = (iring*iring)/4 + iff
       if (ifd < ifi) refl = 1 ! East-West swap


    ! North-Equatorial region
    elseif (pixel < (npix+ns4)/2) then
       templ8 = n0 + (pixel - ncap)/ns4

    ! South-Equatorial region
    elseif (pixel < (npix-ncap)) then
       templ8 = n0 + (npix - 1 - pixel - ncap)/ns4
       refl = 2  ! North-South swap

    ! South polar cap
    else
       ip =  npix - pixel
       fip   = real(int(ip/2) , kind=dp)
       iring = int(sqrt( ip/2.0_dp - sqrt(fip) )) + 1 ! counted from SOUTH pole, starting at 1
       ifi   = MOD( (2*iring*(iring+1) - ip), iring)    ! in [0,iring-1], upwards
       ifd  = iring - 1 - ifi                           ! in [0,iring-1], downwards
       iff  = min(ifi, ifd)                             ! in [0,(iring-1)/2]
       templ8 = (iring*iring)/4 + iff
       if (ifd < ifi) refl = 1 ! East-West swap
       refl = refl + 2    ! North-South swap
    endif

    if (present(reflexion)) reflexion = refl

    return
  end subroutine template_pixel_ring


  !=======================================================================
  subroutine same_shape_pixels_ring(nside, templ8, list, reflexion, nrep)
    !   Returns the list of RING ordered pixels having the same shape as the
    !   template provided
    !-------------------------------------------------------------

    character(len=*), parameter :: code = "same_shape_pixels_ring"
    integer(I4B),    intent(IN)                        :: nside
    integer(I4B),    intent(IN)                        :: templ8
    integer(I4B),    pointer,   dimension(:), optional :: list
    integer(I4B),    pointer,   dimension(:), optional :: reflexion
    integer(I4B),    intent(OUT)            , optional :: nrep

    integer(i4b) :: npix, ntplt, ncap, ns4, n0
    integer(i4b) :: i, in0, is0, iring, ifi
    integer(i4b), dimension(0:3) :: ramp4 = (/ 0, 1, 2, 3 /)
    integer(i4b), dimension(0:7) :: alt8  = (/ 0, 1, 0, 1, 0, 1, 0, 1 /)
    integer(i4b), dimension(0:7) :: rr8
    logical(lgt) :: do_list, do_rfx, do_nrep
    integer(i4b) :: nrep_in
    integer(i4b) :: ll=1, lr=1
    !-------------------------------------------------------------

    do_list = (present(list))
    do_rfx  = (present(reflexion))
    do_nrep = (present(nrep))

    if (do_rfx .and. .not. do_list) then
       print*,'Error in '//code//'. Can not have Reflexion without pixel list'
       call fatal_error()
    endif

    npix = nside2npix(nside)
    if (npix < 0) then
       print*,'Invalid Nside = ',nside
       call fatal_error(code//' Abort')
    endif

    ntplt= nside2ntemplates(nside)

    if (templ8 < 0 .or. templ8 >= ntplt) then
       print*,'Error on template argument'
       print*,'Nside = ',nside,', Template = ',templ8
       print*,'Template should be in [0, (1+Nside*(Nside+6))/4-1=',ntplt-1,']'
       call fatal_error(code//' Abort')
    endif

    ncap = 2*nside*(nside+1)
    ns4  = 4*nside
!    n0   = max( (nside*(nside+2))/4, 1)
    n0 = (1 + nside * (nside + 2))/ 4

    ! find number of repetition of template around the sphere
    if (templ8 >= n0) then
       nrep_in = ns4 * 2 ! Non- Equator: nrep = 8*Nside
       if (templ8 == ntplt - 1) nrep_in = ns4 ! Equator : nrep = 4*Nside
    else
       nrep_in  = 16
       iring = int(sqrt(4*templ8+1.0_dp)) ! in [1,Nside]
       ifi   = templ8 - (iring*iring)/4
       if (2*ifi+1 == iring) nrep_in = 8
    endif
    if (do_nrep) nrep = nrep_in
          

    ! create or expand (or leave alone) list array
    if (do_list) then
       if (associated(list)) then
          ll = lbound(list,dim=1)
          if (size(list) < nrep_in) then
             deallocate(list)
             allocate(list(ll:ll+nrep_in-1))
          endif
       else
          allocate(list(ll:ll+nrep_in-1))
       endif
       list = -1 ! default value for unused array elements
    endif

    ! create or expand (or leave alone) reflexion array
    if (do_rfx) then
       if (associated(reflexion)) then
          lr = lbound(reflexion,dim=1)
          if (size(reflexion) < nrep_in) then
             deallocate(reflexion)
             allocate(reflexion(lr:lr+nrep_in-1))
          endif
       else
          allocate(reflexion(lr:lr+nrep_in-1))
       endif
       reflexion = -1 ! default value for unused array elements
    endif


    ! fill up list and reflexion arrays
    if (templ8 >= n0) then 
       ! Non Polar pixels (most frequent)
    
       in0 = ncap + (templ8 - n0)*ns4

       if (templ8 == ntplt-1) then 
          ! Equator : nrep = 4*Nside
          if (do_list) list      = in0 + (/ (i, i=0, ns4-1) /)
          if (do_rfx)  reflexion = 0
       else 
          ! Non- Equator: nrep = 8*Nside
          is0 = npix - in0 - ns4
          if (do_list) then
             list(ll    :ll+  ns4-1) = in0 + (/ (i, i=0, ns4-1) /)
             list(ll+ns4:ll+2*ns4-1) = is0 + (/ (i, i=0, ns4-1) /)
          endif
          if (do_rfx) then
             reflexion(lr    :lr+  ns4-1) = 0
             reflexion(lr+ns4:lr+2*ns4-1) = 2
          endif
       endif

    else 
       ! Polar pixels
       in0   = 2*iring*(iring-1)
       is0   = npix - in0 - 4*iring

       if (2*ifi+1 == iring) then 
          ! pixel on phi = Pi/4, Nrep = 8
          if (do_list) then
             list(ll  :ll+3) = in0 + ifi + ramp4*iring
             list(ll+4:ll+7) = is0 + ifi + ramp4*iring
          endif
          if (do_rfx) then
             reflexion = (/ 0, 0, 0, 0, 2, 2, 2, 2 /)
          endif
       else 
          ! Nrep = 16
          rr8(0:6:2) = (/ 0, 1, 2, 3 /) * iring + ifi
          rr8(1:7:2) = (/ 1, 2, 3, 4 /) * iring - ifi - 1
          if (do_list) then
             list(ll  :ll+7 ) = in0 + rr8
             list(ll+8:ll+15) = is0 + rr8
          endif
          if (do_rfx) then
             reflexion(lr  :lr+7 ) = alt8
             reflexion(lr+8:lr+15) = alt8 + 2
          endif
       endif
    endif

    return

  end subroutine same_shape_pixels_ring

  !=======================================================================
  subroutine template_pixel_nest(nside, pixel, templ8, reflexion)
    !=======================================================================
    !    Returns template-pixel index corresponding to each (NESTED ordered) pixel provided
    !    The reflexion (E-W, or N-W or both) to apply to the pixel to match the
    !    template is optionally returned
    !-------------------------------------------------------------
    character(len=*), parameter :: code = "Template_pixel_nest"
    integer(I4B),    intent(IN)           :: nside
    integer(I4B),    intent(IN)           :: pixel
    integer(I4B),    intent(OUT)          :: templ8
    integer(I4B),    intent(OUT), optional:: reflexion
    !
    integer(I4B) :: pix_ring, npix
    
    npix = nside2npix(nside)
    if (npix < 0) then
       print*,'Invalid Nside = ',nside
       call fatal_error(code//' Abort')
    endif

    if (pixel < 0 .or. pixel >= npix) then
       print*,'Invalid Pixel = ',pixel
       call fatal_error(code//' Abort')
    endif

    call nest2ring(nside, pixel, pix_ring)
    call template_pixel_ring(nside, pix_ring, templ8, reflexion)

    return
  end subroutine template_pixel_nest

  !=======================================================================
  subroutine same_shape_pixels_nest(nside, templ8, list, reflexion, nrep)
    !=======================================================================
    !   Returns the list of NESTED ordered pixels having the same shape as the
    !   template provided
    !-------------------------------------------------------------
    use num_rec, only : iindexx

    character(len=*), parameter :: code = "same_shape_pixels_nest"
    integer(I4B),    intent(IN)                        :: nside
    integer(I4B),    intent(IN)                        :: templ8
    integer(I4B),    pointer,   dimension(:), optional :: list
    integer(I4B),    pointer,   dimension(:), optional :: reflexion
    integer(I4B),    intent(OUT)            , optional :: nrep

    integer(i4b) :: ntplt
    integer(i4b) :: i
    integer(I4B),   allocatable,   dimension(:)   :: idx, tmp
    integer(i4b) :: ll, lr
    logical(lgt) :: do_list, do_rfx, do_nrep
    integer(i4b) :: nr
!     integer(I4B),    pointer,   dimension(:) :: list_in
!     integer(I4B),    pointer,   dimension(:) :: reflexion_in
!     integer(I4B),    intent(OUT)             :: nrep_in
    !-------------------------------------------------------------

    do_list = (present(list))
    do_rfx  = (present(reflexion))
    do_nrep = (present(nrep))
  
    if (do_rfx .and. .not. do_list) then
       print*,'Error in '//code//'. Can not have Reflexion without pixel list'
       call fatal_error
    endif

    ntplt= nside2ntemplates(nside)

    if (templ8 < 0 .or. templ8 >= ntplt) then
       print*,'Error on template argument'
       print*,'Nside = ',nside,', Template = ',templ8
       print*,'Template should be in [0, (1+Nside*(Nside+6))/4-1=',ntplt-1,']'
       call fatal_error(code//' Abort')
    endif

    ! makes pixel list in RING indexing
    call same_shape_pixels_ring( nside, templ8, list, reflexion, nr)
    if (do_nrep) nrep = nr

    if (do_list .or. do_rfx) then
       allocate(idx  (0:nr-1))
       allocate(tmp  (0:nr-1))

       ! convert to NESTED
       if (do_list) then
          ll = lbound(list,      dim=1)
          do i = 0, nr-1
             call ring2nest(nside, list(ll+i), list(ll+i))
          enddo

          ! sort arrays in increasing order of pixel NESTED index
          call iindexx(nr, list, idx)
          !!!!!  list(ll:nr-1+ll) = list(idx-1+ll) ! replaced by loops for Gfortran
          do i=0, nr-1
             tmp(i) = list(idx(i)-1+ll)
          enddo
          do i=0, nr-1
             list(ll+i) = tmp(i)
          enddo
       endif
  
       if (do_rfx) then
          lr = lbound(reflexion,      dim=1)
          !!!!! reflexion(lr:nr-1+lr) = reflexion(idx-1+lr) ! replaced by loops for Gfortran
          do i=0, nr-1
             tmp(i) = reflexion(idx(i)-1+lr)
          enddo
          do i=0, nr-1
             reflexion(lr+i) = tmp(i)
          enddo
       endif

       deallocate(tmp)
       deallocate(idx)
    endif

    return
  end subroutine same_shape_pixels_nest


end module pix_tools
