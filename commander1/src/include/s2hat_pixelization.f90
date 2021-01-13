
module s2hat_pixelization

use healpix_types
use s2hat_types
use alm_tools

! contains pixelization specific routines which create structures needed
! for the s2hat transforms to work - RS@APC February, 2007
! pixelization structure updated: kphi -> dp; qwght added - rs@apc Oct, 2007
! mask2scan fixed - rs@apc Oct, 2007

implicit none

! include "mpif.h"

public :: set_pixelization, zbounds2scan, zbounds2mask, mask2scan, destroy_pixelization, destroy_scan
public :: mpi_scanbcast, mpi_pixelizationbcast

contains

! generic

subroutine set_pixelization( pixchoice, pixpars, pixelization)

  implicit none

  !input 
  integer(i4b), intent(in) :: pixchoice
  type( pixparameters), intent(in) :: pixpars
  ! output
  type( pixeltype), intent( out) :: pixelization

  if( pixchoice == PIXCHOICE_HEALPIX) then
     call set_healpix_pixelization( pixpars%par1, pixelization)
  elseif (pixchoice == PIXCHOICE_GLESP) then
     call set_glesp_pixelization( pixpars%par1, pixpars%par2, pixelization)
  elseif (pixchoice == PIXCHOICE_ECP) then
     call set_ecp_pixelization( pixpars%par1, pixpars%par2, pixelization)
  elseif (pixchoice == PIXCHOICE_GLCP) then
     call set_glcp_pixelization( pixpars%par1, pixpars%par2, pixelization)
  else
     write(*,*)'wrong pixelization choice'
     stop 
  endif

  return

end subroutine set_pixelization

subroutine destroy_pixelization( pixelization)

  implicit none

  ! input
  type( pixeltype) :: pixelization

  deallocate( pixelization%fpix)
  deallocate( pixelization%nph)
  deallocate( pixelization%kphi)
  deallocate( pixelization%qwght)
  deallocate( pixelization%pixphi)
  deallocate( pixelization%parea)
  deallocate( pixelization%cth)
  deallocate( pixelization%sth)

end subroutine destroy_pixelization

subroutine mpi_scanbcast( pixelization, scan, root, my_rank, current_comm)

  implicit none

  ! input
  type( pixeltype), intent( in) :: pixelization
  integer(i4b) :: root, my_rank, current_comm

  ! output (and input on the root)
  type( scandef) :: scan

  ! internal 
  integer(i4b) :: ntmp, mpierr, nringsall
  integer(i4b), dimension(:), allocatable :: ivect

  if( my_rank == root) ntmp = int( scan%npixsobs, i4b)
  call mpi_bcast( ntmp, 1, mpi_integer, root, current_comm, mpierr)
  scan%npixsobs = int( ntmp, i8b)

  if( my_rank == root) ntmp = int( scan%nringsobs, i4b)
  call mpi_bcast( ntmp, 1, mpi_integer, root, current_comm, mpierr)
  scan%nringsobs = int( ntmp, i8b)  

  nringsall = int( pixelization%nringsall, i4b)

  if( my_rank /= root) then

     allocate( scan%nfl(1:nringsall))
     allocate( scan%sfl(1:nringsall))
     allocate( scan%fl(1:nringsall))

  end if

  allocate( ivect(1:nringsall))
  
  if( my_rank == root) then
     ivect(1:nringsall) = scan%nfl(1:nringsall)
  endif

  call mpi_bcast( ivect, nringsall, mpi_integer, root, current_comm, mpierr);

  if( my_rank /= root) then
     scan%nfl(1:nringsall) = int( ivect(1:nringsall), i8b)
  endif

  if( my_rank == root) then
     ivect(1:nringsall) = scan%sfl(1:nringsall)
  endif

  call mpi_bcast( ivect, nringsall, mpi_integer, root, current_comm, mpierr);

  if( my_rank /= root) then
     scan%sfl(1:nringsall) = int( ivect(1:nringsall), i8b)
  endif

  if( my_rank == root) then
     ivect(1:nringsall) = scan%fl(1:nringsall)
  endif

  call mpi_bcast( ivect, nringsall, mpi_integer, root, current_comm, mpierr);

  if( my_rank /= root) then
     scan%fl(1:nringsall) = int( ivect(1:nringsall), i8b)
  endif

  deallocate( ivect)

  return

end subroutine mpi_scanbcast

subroutine mpi_pixelizationbcast( pixelization, root, my_rank, current_comm)

  implicit none

  ! input
  integer(i4b) :: root, my_rank, current_comm

  ! output (and input on the root)
  type( pixeltype) :: pixelization

  ! internal 
  integer(i4b) :: ntmp, mpierr
  integer( i4b), dimension(:), allocatable :: ivect
  real( dp), dimension(:), allocatable :: fvect

  if( my_rank == root) ntmp = int( pixelization%npixsall, i4b)
  call mpi_bcast( ntmp, 1, mpi_integer, root, current_comm, mpierr)
  pixelization%npixsall = int( ntmp, i8b)

  if( my_rank == root) ntmp = int( pixelization%nringsall, i4b)
  call mpi_bcast( ntmp, 1, mpi_integer, root, current_comm, mpierr)
  pixelization%nringsall = int( ntmp, i8b)  

  if( my_rank /= root) then

     allocate( pixelization%fpix(1:ntmp))
     allocate( pixelization%nph(1:ntmp))
     allocate( pixelization%kphi(1:ntmp))
     allocate( pixelization%qwght(1:ntmp))
     allocate( pixelization%pixphi(1:ntmp))
     allocate( pixelization%parea(1:ntmp))
     allocate( pixelization%cth(1:ntmp))
     allocate( pixelization%sth(1:ntmp))

  end if

  allocate( ivect(1:ntmp))
  
  if( my_rank == root) then
     ivect(1:ntmp) = pixelization%fpix(1:ntmp)
  endif

  call mpi_bcast( ivect, ntmp, mpi_integer, root, current_comm, mpierr);

  if( my_rank /= root) then
     pixelization%fpix(1:ntmp) = int( ivect(1:ntmp), i8b)
  endif

  if( my_rank == root) then
     ivect(1:ntmp) = pixelization%nph(1:ntmp)
  endif

  call mpi_bcast( ivect, ntmp, mpi_integer, root, current_comm, mpierr);

  if( my_rank /= root) then
     pixelization%nph(1:ntmp) = int( ivect(1:ntmp), i8b)
  endif

  if( my_rank == root) then
     ivect(1:ntmp) = pixelization%kphi(1:ntmp)
  endif

  call mpi_bcast( ivect, ntmp, mpi_integer, root, current_comm, mpierr);

  if( my_rank /= root) then
     pixelization%kphi(1:ntmp) = int( ivect(1:ntmp), i8b)
  endif

  deallocate( fvect)

  ! double precision objects now

  allocate( fvect(1:ntmp))

  if( my_rank == root) then
      fvect(1:ntmp) = pixelization%qwght(1:ntmp)
  endif

  call mpi_bcast( fvect, ntmp, mpi_double_precision, root, current_comm, mpierr);

  if( my_rank /= root) then
      pixelization%qwght(1:ntmp) = fvect(1:ntmp)
  endif

  if( my_rank == root) then
      fvect(1:ntmp) = pixelization%pixphi(1:ntmp)
  endif

  call mpi_bcast( fvect, ntmp, mpi_double_precision, root, current_comm, mpierr);

  if( my_rank /= root) then
      pixelization%pixphi(1:ntmp) = fvect(1:ntmp)
  endif

  if( my_rank == root) then
      fvect(1:ntmp) = pixelization%parea(1:ntmp)
  endif

  call mpi_bcast( fvect, ntmp, mpi_double_precision, root, current_comm, mpierr);

  if( my_rank /= root) then
      pixelization%parea(1:ntmp) = fvect(1:ntmp)
  endif

  if( my_rank == root) then
      fvect(1:ntmp) = pixelization%cth(1:ntmp)
  endif

  call mpi_bcast( fvect, ntmp, mpi_double_precision, root, current_comm, mpierr);

  if( my_rank /= root) then
      pixelization%cth(1:ntmp) = fvect(1:ntmp)
  endif

  if( my_rank == root) then
      fvect(1:ntmp) = pixelization%sth(1:ntmp)
  endif

  call mpi_bcast( fvect, ntmp, mpi_double_precision, root, current_comm, mpierr);

  if( my_rank /= root) then
      pixelization%sth(1:ntmp) = fvect(1:ntmp)
  endif

  deallocate( fvect)

  return

end subroutine mpi_pixelizationbcast

subroutine zbounds2scan( zbounds, pixelization, scan)

  implicit none

  !input
  real(dp), dimension(2) :: zbounds
  type( pixeltype) :: pixelization
  ! output
  type( scandef) :: scan

  ! internal
  integer(i4b) :: iring, npixsobs, nringsobs, nringsall, equatorIn
  logical(lgt) :: keep_north, keep_south, keep_it
  real(dp) :: cth

  nringsall = pixelization%nringsall

  allocate( scan%nfl(1:nringsall))
  allocate( scan%sfl(1:nringsall))
  allocate( scan%fl(1:nringsall))

  ! write(*,*)pixelization%fpix(3), pixelization%nph(3), pixelization%cth(1), pixelization%sth(1)

  ! check if there is an equatorial ring of pixels
  equatorIn = equatorialRing( pixelization)

  npixsobs = 0
  nringsobs = 0

  do iring = 1, nringsall

     cth = pixelization%cth(iring)
     call select_rings( cth, zbounds, keep_north, keep_south, keep_it)

     if( keep_north) then
        scan%nfl(iring) = 1
        npixsobs = npixsobs + pixelization%nph(iring)
     else
        scan%nfl(iring) = 0
     endif

     if( keep_south .and. (iring <= nringsall-equatorIn)) then  ! make sure to include the equator (if present) only in the north hemisphere
        scan%sfl(iring) = 1
        npixsobs = npixsobs + pixelization%nph(iring)
     else
        scan%sfl(iring) = 0
     endif

     if( keep_it) then
        scan%fl(iring) = 1
        nringsobs = nringsobs + 1
     else
        scan%fl(iring) = 0
     endif

  end do
  
  scan%npixsobs = npixsobs
  scan%nringsobs = nringsobs

  return

end subroutine zbounds2scan

subroutine zbounds2mask( zbounds, pixelization, mask)

  implicit none

  !input
  real(dp), dimension(2), intent(in) :: zbounds
  type( pixeltype), intent(in) :: pixelization

  ! output
  integer(i4b), dimension(:), pointer :: mask

  ! internal
  integer(i4b) :: iring, iringFirst, iringLast
  real(dp) :: costh

  mask(:) = 0

  ! consistent with healpix's select_rings 

  if( abs( zbounds(1)-zbounds(2)) < 1.d-6) then  ! include all

      mask(:) = 1

  else

      if( zbounds(1) < zbounds(2)) then    ! only a strip retained

          do iring = 1, pixelization%nringsall

             ! north hemisphere

             costh = pixelization%cth(iring)

             iringFirst = int( pixelization%fpix(iring), i4b)  
             iringLast = int( pixelization%fpix(iring)+pixelization%nph(iring), i4b)-1

             if( (costh <= zbounds(2)) .and. (costh >= zbounds(1))) then
                 mask( iringFirst:iringLast) = 1
             endif

             ! south hemisphere

             if( costh /= 0.0) then

                 costh = -costh
                 iringFirst = int( pixelization%npixsall, i4b) - iringFirst-1  ! first <-> last is realy swaped here
                 iringLast = int( pixelization%npixsall, i4b) - iringLast-1

                 if( (costh <= zbounds(2)) .and. (costh >= zbounds(1))) then
                     mask( iringLast:iringFirst) = 1
                 endif

             endif

         end do

      else       ! strip to be removed from the map

         do iring = 1, pixelization%nringsall

            ! north hemisphere

            costh = pixelization%cth(iring)

            iringFirst = int( pixelization%fpix(iring), i4b)  
            iringLast = int( pixelization%fpix(iring)+pixelization%nph(iring), i4b)-1

            if( (costh < zbounds(2)) .or. (costh > zbounds(1))) then
                mask( iringFirst:iringLast) = 1
            endif

            ! south hemisphere

            if( costh /= 0.0) then

                costh = -costh
                iringFirst = int( pixelization%npixsall, i4b) - iringFirst-1  ! first <-> last is realy swaped here
                iringLast = int( pixelization%npixsall, i4b) - iringLast-1

                if( (costh < zbounds(2)) .and. (costh > zbounds(1))) then
                    mask( iringLast:iringFirst) = 1
                endif

            endif

         end do

    end if

  end if

  return

end subroutine zbounds2mask

subroutine mask2scan( mask, pixelization, scan)

  implicit none

  ! input
  integer(i4b), dimension(:), pointer :: mask
  type(pixeltype), intent(in) :: pixelization
  ! output
  type( scandef), intent(out) :: scan

  ! internal

  integer(i4b) :: ip, iring, iring_last, lastpix, lastring, nrings, npix, equatorIn
  integer(i8b), dimension(:), allocatable :: rings_north, rings_south
 
  ! northern hemisphere
 
  allocate( rings_north(1:pixelization%nringsall))

  rings_north(:) = 0

  iring_last = 1
  lastring = int( pixelization%nringsall, i4b)
  lastpix = pixelization%fpix(lastring)+pixelization%nph(lastring)

  ip = 0

  do while (ip < lastpix)

     if( mask(ip) == 1) then

        do iring = iring_last, pixelization%nringsall

           if( (pixelization%fpix(iring) <= ip) .and. (pixelization%fpix(iring)+pixelization%nph(iring)) > ip) exit

        end do

        rings_north( iring) = 1

        ! go right away to the next ring and a pixel
        ip = pixelization%fpix(iring)+pixelization%nph(iring)

        iring_last = iring+1

     end if

     ip = ip + 1

  end do

  ! southern hemisphere

  allocate( rings_south(1:pixelization%nringsall))

  rings_south(:) = 0

  ! check if there is an equatorial ring of pixels
  equatorIn = equatorialRing( pixelization)

  iring_last = 1

  lastpix = ip

  ip = pixelization%npixsall-1

  do while (ip >= lastpix)

     if( mask(ip) == 1) then

        do iring = iring_last, pixelization%nringsall-int( equatorIn, i8b)  ! do not double equator ring if it is there

           if( (pixelization%fpix(iring) <= pixelization%npixsall-ip-1) .and. (pixelization%fpix(iring)+pixelization%nph(iring)) > pixelization%npixsall-ip-1) exit

        end do

        rings_south( iring) = 1

        ! go right away to the next ring and a pixel
        ip = pixelization%npixsall-pixelization%fpix(iring)-pixelization%nph(iring)-1

        iring_last = iring+1

     end if

     ip = ip - 1

  end do

  allocate( scan%nfl(1:pixelization%nringsall))
  allocate( scan%sfl(1:pixelization%nringsall))
  allocate( scan%fl(1:pixelization%nringsall))

  nrings = 0
  do iring = 1, pixelization%nringsall
 
     scan%nfl( iring) = rings_north( iring)
     scan%sfl( iring) = rings_south( iring)

     if( (rings_north( iring) + rings_south( iring)) /= 0) then 
         scan%fl( iring) = 1
         nrings = nrings+1
     else
         scan%fl( iring) = 0
     endif
  end do

  scan%nringsobs = nrings

  npix = 0
  do ip = 0, pixelization%npixsall-1
     if( mask(ip) == 1) npix = npix + 1
  end do

  scan%npixsobs = npix
  
  deallocate( rings_north, rings_south)

  return

end subroutine mask2scan

function equatorialRing( pixelization)

   implicit none

   !input 
   type( pixeltype), intent(in) :: pixelization

   ! internal
   integer(i4b) :: nringlast

   ! output
   integer(i4b) :: equatorialRing

   nringlast = pixelization%nringsall

   if( pixelization%npixsall < (2*(pixelization%fpix( nringlast) + pixelization%nph( nringlast)-1))) then
       equatorialRing = 1
   else
       equatorialRing = 0
   end if

   return

end function equatorialRing

subroutine destroy_scan( scan)

  implicit none
  
  ! input
  type( scandef) :: scan

  deallocate( scan%nfl)
  deallocate( scan%sfl)
  deallocate( scan%fl)

end subroutine destroy_scan

! Healpix specific

subroutine set_healpix_pixelization( nsmax, pixelization)

  implicit none

  !input 
  integer(i4b) :: nsmax
  ! output
  type( pixeltype) :: pixelization

  ! internal
  integer(i4b) :: iring, nringsall, kphi, nph
  integer(i8b) :: fpix
  real(dp) :: cth, sth

  pixelization%type = PIXCHOICE_HEALPIX

  nringsall = 2*nsmax
  pixelization%nringsall = nringsall
  pixelization%npixsall = 12*nsmax*nsmax

  pixelization%nphmx = 4*nsmax

  allocate( pixelization%fpix(1:nringsall))
  allocate( pixelization%nph(1:nringsall))
  allocate( pixelization%kphi(1:nringsall))
  allocate( pixelization%qwght(1:nringsall))
  allocate( pixelization%pixphi(1:nringsall))
  allocate( pixelization%parea(1:nringsall))
  allocate( pixelization%cth(1:nringsall))
  allocate( pixelization%sth(1:nringsall))

  do iring = 1, nringsall

     call get_pixel_layout( nsmax, iring, cth, sth, nph, fpix, kphi)

     pixelization%cth(iring) = cth
     pixelization%sth(iring) = sth
     pixelization%nph(iring) = nph
     pixelization%kphi(iring) = dble( kphi)*PI/dble( nph) ! first pixel offset in radians (half pix size)
     pixelization%fpix(iring) = fpix
     pixelization%qwght(iring) = 1.0    ! equal weights

     pixelization%parea(iring) = 4.0*PI/(1.0*pixelization%npixsall)
     pixelization%pixphi(iring) = 2.0*PI/(1.0*nph)

  enddo
  
  return

end subroutine set_healpix_pixelization

! GLESP pixelization

subroutine set_glesp_pixelization( ntheta, nphi, pixelization)

  implicit none

  !input
  integer(i4b)     :: ntheta, nphi
  ! output
  type( pixeltype) :: pixelization

  ! internal
  integer(i4b) :: iring, nringsall
  integer(i8b) :: fpix
  real(dp) :: sin_curnt, sin_next, sin_prev, thetaFact, pixArea
  real(dp) :: cosdth1, cosdth2, sindth1, sindth2, totAng, totArea
  real(dp), dimension(:), allocatable :: zeros, wghts

  pixelization%type = PIXCHOICE_GLESP

  nringsall = (ntheta+1)/2                        ! only North Hemisphere and the equator if ntheta odd
  pixelization%nringsall = nringsall

  pixelization%nphmx = nphi                       ! max number of pixs per ring = number of pixs at the equator

  allocate( pixelization%fpix(1:nringsall))
  allocate( pixelization%nph(1:nringsall))
  allocate( pixelization%kphi(1:nringsall))
  allocate( pixelization%qwght(1:nringsall))
  allocate( pixelization%pixphi(1:nringsall))
  allocate( pixelization%parea(1:nringsall))
  allocate( pixelization%cth(1:nringsall))
  allocate( pixelization%sth(1:nringsall))

  allocate( zeros(1:ntheta))
  allocate( wghts(1:ntheta))
  call setGaussLegQuad( ntheta, zeros, wghts)

  ! approximate pixel area = equatorial ring (or the one next to it if equator missing) pixel area

  sin_next = dsqrt ( 1.d0 - zeros( nringsall+1)*zeros( nringsall+1))
  sin_curnt = dsqrt ( 1.d0 - zeros( nringsall)*zeros( nringsall))
  sin_prev = dsqrt ( 1.d0 - zeros( nringsall-1)*zeros( nringsall-1))

  cosdth1 = dsqrt( 0.5d0*(1.d0+zeros(nringsall)*zeros(nringsall-1)+sin_curnt*sin_prev))
  cosdth2 = dsqrt( 0.5d0*(1.d0+zeros(nringsall+1)*zeros(nringsall)+sin_next*sin_curnt))

  sindth1 = dsqrt( 0.5d0*(1.d0-zeros(nringsall)*zeros(nringsall-1)-sin_curnt*sin_prev))
  sindth2 = dsqrt( 0.5d0*(1.d0-zeros(nringsall+1)*zeros(nringsall)-sin_next*sin_curnt))

  thetaFact = zeros(nringsall)*(cosdth1-cosdth2)+sin_curnt*(sindth1+sindth2)
  pixArea = thetaFact*2.0*PI/pixelization%nphmx

  ! and now ring-by-ring 

  sin_next = dsqrt( 1.d0-zeros(1)*zeros(1))

  fpix = 0

  totArea = 0.d0
  totAng  = 0.d0

  do iring = 1, nringsall

     pixelization%cth(iring) = zeros( iring)

     pixelization%sth(iring) = sin_next
     sin_next = dsqrt( 1.d0-zeros(iring+1)*zeros(iring+1))

     if( iring > 1) then
         cosdth1 = sqrt( 0.5d0*(1.d0+zeros(iring)*zeros(iring-1)+pixelization%sth(iring)*sin_prev))
         cosdth2 = sqrt( 0.5d0*(1.d0+zeros(iring+1)*zeros(iring)+sin_next*pixelization%sth(iring)))

         sindth1 = sqrt( 0.5d0*(1.d0-zeros(iring)*zeros(iring-1)-pixelization%sth(iring)*sin_prev))
         sindth2 = sqrt( 0.5d0*(1.d0-zeros(iring+1)*zeros(iring)-sin_next*pixelization%sth(iring)))

         thetaFact = pixelization%cth(iring)*(cosdth1-cosdth2)+pixelization%sth(iring)*(sindth1+sindth2)

     else

         cosdth2 = sqrt( 0.5d0*(1.d0+zeros(iring+1)*zeros(iring)+sin_next*pixelization%sth(iring)))
         sindth2 = sqrt( 0.5d0*(1.d0-zeros(iring+1)*zeros(iring)-sin_next*pixelization%sth(iring)))

         thetaFact = 1.d0 - pixelization%cth( iring)*cosdth2 + pixelization%sth(iring)*sindth2             ! pole ring of pixels

     end if 

     pixelization%nph(iring) = int( 2*PI*thetaFact/pixArea+0.5, i8b)

     if( pixelization%nph(iring) <= 0) then
        write(*,*)'oops', pixelization%nph(iring)
     end if

     pixelization%parea(iring) = 2.0*PI/pixelization%nph(iring)*thetaFact                                  ! true area for this row of pixs

     pixelization%qwght(iring) = dabs( wghts( iring))/thetaFact

     pixelization%kphi(iring) = 0.d0                                               ! first pix of each ring is centered on meridian zero !
     pixelization%fpix(iring) = fpix

     pixelization%pixphi(iring) = 2.0*PI/pixelization%nph(iring)

     sin_prev = pixelization%sth(iring)

     fpix = fpix + pixelization%nph(iring)

     totAng  = totAng + 2.0*thetaFact
     totArea = totArea + 2.0*pixelization%nph(iring)*pixelization%parea(iring)

  enddo

  pixelization%npixsall = 2*fpix                           ! total number of pixs

  ! do not double count the equatorial ring
  if( ntheta /= 2*(ntheta/2)) then
      pixelization%npixsall = pixelization%npixsall - pixelization%nph(nringsall)
      totArea = totArea - pixelization%parea(nringsall)*pixelization%nph(nringsall)
      totAng = totAng - thetaFact
  end if

  deallocate( zeros, wghts)
  
  return

end subroutine set_glesp_pixelization


! ECP definitions

subroutine set_ecp_pixelization( ntheta, nphi, pixelization)

  implicit none

  !input
  integer(i4b) :: ntheta, nphi
  ! output
  type( pixeltype) :: pixelization

  ! internal
  integer(i4b) :: iring, nringsall, npixperring
  real(dp) :: theta, dtheta

  pixelization%type = PIXCHOICE_ECP

  nringsall = (ntheta+1)/2                        ! only North Hemisphere and the equator, if ntheta odd
  pixelization%nringsall = nringsall

  npixperring = nphi
  pixelization%nphmx = npixperring

  pixelization%npixsall = npixperring*ntheta      ! total number of pixs

  allocate( pixelization%fpix(1:nringsall))
  allocate( pixelization%nph(1:nringsall))
  allocate( pixelization%kphi(1:nringsall))
  allocate( pixelization%qwght(1:nringsall))
  allocate( pixelization%pixphi(1:nringsall))
  allocate( pixelization%parea(1:nringsall))
  allocate( pixelization%cth(1:nringsall))
  allocate( pixelization%sth(1:nringsall))

  dtheta = PI/ntheta
  theta = 0.5*dtheta

  do iring = 1, nringsall

     pixelization%cth(iring) = dcos( theta)
     pixelization%sth(iring) = dsin( theta)
     pixelization%nph(iring) = npixperring
     pixelization%kphi(iring) = 0.d0            ! first pix of each ring is centered on meridian zero !
     pixelization%fpix(iring) = (iring-1)*npixperring

     pixelization%parea(iring) = 2.0*PI/npixperring*(dcos( theta-0.5*dtheta) - dcos( theta+0.5*dtheta))
     pixelization%qwght(iring) = 1.0

     pixelization%pixphi(iring) = 2.0*PI/npixperring

     theta = theta + dtheta

  enddo
  
  return

end subroutine set_ecp_pixelization

! GLCP definitions

subroutine set_glcp_pixelization( ntheta, nphi, pixelization)

  implicit none

  !input
  integer(i4b) :: ntheta, nphi
  ! output
  type( pixeltype) :: pixelization

  ! internal
  integer(i4b) :: iring, nringsall, npixperring
  real(dp) :: sin_next, sin_prev, tmp, totArea, totAng
  real(dp) :: cosdth1, cosdth2, sindth1, sindth2
  real(dp), dimension(:), allocatable :: zeros, wghts

  pixelization%type = PIXCHOICE_GLCP

  nringsall = (ntheta+1)/2                        ! only North Hemisphere and the equator if ntheta odd
  pixelization%nringsall = nringsall

  npixperring = nphi
  pixelization%nphmx = npixperring

  pixelization%npixsall = npixperring*ntheta      ! total number of pixs

  allocate( pixelization%fpix(1:nringsall))
  allocate( pixelization%nph(1:nringsall))
  allocate( pixelization%kphi(1:nringsall))
  allocate( pixelization%qwght(1:nringsall))
  allocate( pixelization%pixphi(1:nringsall))
  allocate( pixelization%parea(1:nringsall))
  allocate( pixelization%cth(1:nringsall))
  allocate( pixelization%sth(1:nringsall))

  allocate( zeros(1:ntheta))
  allocate( wghts(1:ntheta))
  call setGaussLegQuad( ntheta, zeros, wghts)

  sin_next = dsqrt( 1.d0-zeros(1)*zeros(1))

  totArea = 0.d0
  totAng = 0.d0

  do iring = 1, nringsall

     pixelization%cth(iring) = zeros( iring)

     pixelization%sth(iring) = sin_next
     sin_next = dsqrt( 1.d0-zeros(iring+1)*zeros(iring+1))

     pixelization%nph(iring) = npixperring
     pixelization%kphi(iring) = 0.d0            ! first pix of each ring is centered on meridian zero !
     pixelization%fpix(iring) = (iring-1)*npixperring

     if( iring > 1) then
         cosdth1 = sqrt( 0.5d0*(1.d0+zeros(iring)*zeros(iring-1)+pixelization%sth(iring)*sin_prev))
         cosdth2 = sqrt( 0.5d0*(1.d0+zeros(iring+1)*zeros(iring)+sin_next*pixelization%sth(iring)))

         sindth1 = sqrt( 0.5d0*(1.d0-zeros(iring)*zeros(iring-1)-pixelization%sth(iring)*sin_prev))
         sindth2 = sqrt( 0.5d0*(1.d0-zeros(iring+1)*zeros(iring)-sin_next*pixelization%sth(iring)))

         tmp = pixelization%cth(iring)*(cosdth1-cosdth2)+pixelization%sth(iring)*(sindth1+sindth2)
         pixelization%parea(iring) = tmp*2.0*PI/npixperring
     else

         cosdth2 = sqrt( 0.5d0*(1.d0+zeros(iring+1)*zeros(iring)+sin_next*pixelization%sth(iring)))
         sindth2 = sqrt( 0.5d0*(1.d0-zeros(iring+1)*zeros(iring)-sin_next*pixelization%sth(iring)))

         tmp = 1.d0 - pixelization%cth( iring)*cosdth2 + pixelization%sth(iring)*sindth2             ! pole ring of pixels
         pixelization%parea(iring) = tmp*2.0*PI/npixperring
     end if 

     pixelization%qwght(iring) = dabs( wghts( iring))/tmp

     pixelization%pixphi(iring) = 2.0*PI/npixperring

     sin_prev = pixelization%sth(iring)

     totArea = totArea + 2.d0*pixelization%parea(iring)*pixelization%nph(iring)
     totAng = totAng + 2.0*tmp

  enddo

  if( ntheta /= 2*(ntheta/2)) totArea = totArea - pixelization%parea(nringsall)*pixelization%nph(nringsall)   ! do not double count the equatorial ring
  if( ntheta /= 2*(ntheta/2)) totAng = totAng - tmp   ! do not double count the equatorial ring

  deallocate( zeros, wghts)
  
  return

end subroutine set_glcp_pixelization

subroutine setGaussLegQuad( n, costh, wghts)

   ! adopted from the NR C-book ...

   use healpix_types

   implicit none

   ! input
   integer(i4b) :: n

   ! output
   real(dp), dimension(:) :: costh, wghts

   ! internal
   real(dp), parameter :: EPS = 1.d-16
   real(dp) :: med, hwdth, z, pp, p3, p2, p1
   integer(i4b) :: i, j, mpts

   mpts = (n+1)/2
   med =  0.d0  ! (upperlim(=-1)+bottomlim(=1))/2
   hwdth = -1.d0  ! (upperlim-botomlim)/2
  
   do i = 1, mpts

     ! starting approx to the roots
     z = dcos( PI*(dble(i)-0.25d0)/(dble(n)+0.5d0))
   
     ! and its refinement
     do
   	! recurrence of the Legendre polynomials
        p1 = 1.d0
	p2 = 0.d0
		
        do j = 1, n
           p3 = p2
	   p2 = p1
	   p1 = ((2.d0*j-1.d0)*z*p2-(dble(j)-1.d0)*p3)/dble(j)
        end do
		
        pp = n*(z*p1-p2)/(z*z-1.d0)

	z = z - p1/pp

	if( dabs(p1/pp) < EPS) exit

    end do
	
    ! prepare the output
     costh(i) = med - hwdth*z
     costh(n+1-i) = med + hwdth*z
     wghts(i) = 2.d0*hwdth/(1.d0-z*z)/pp/pp
     wghts(n+1-i) = wghts(i)
  end do

end subroutine setGaussLegQuad

end module s2hat_pixelization  
