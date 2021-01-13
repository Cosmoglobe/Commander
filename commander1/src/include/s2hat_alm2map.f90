
!--------------------------------------------------------------------------
!
!  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon, 
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
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
!--------------------------------------------------------------------------
! Written by Hans Kristian Eriksen and Snorre Boasson, 
! but copying large parts from the serial HEALPix code.
!--------------------------------------------------------------------------
!
! December 2006
!
! Modified by Radek Stompor (APC, Paris)
!
! The current version allows to distribute both pixel and harmonic domain 
! object over all processors. The spherical harmonic transform is split 
! into two steps separated by the data redistribution. The latter includes 
! an mpi process (all2all) over all the procs (routines redistInterProducts*).
!
! This file contains alm-to-map transforms. The complementary file
! s2hat_map2alm.f90 contains the alm-to-map transforms.
! Routines performing the load balanced data distribution and some other 
! auxiliary procedures are included in the file: s2hat_toolbox.f90.
!
!--------------------------------------------------------------------------
!
! February, 2007 RS@APC
!
! - load balancing improved for partial sky experiments
! - generalized to arbitrary isolatitudal pixelization scheme with same area pixels 
!   equidistant in azimuth for each isolatitude ring, and symmetric wrt to
!   the equator.
!
!--------------------------------------------------------------------------
!
! Oct, 2007 - rs@apc
!
! quadrature wghts and arbitrary pixels ring offsest added (i.e., adapted
! to the new pixelization struct format)
!
!--------------------------------------------------------------------------

module s2hat_alm2map_mod

  use healpix_types
  use s2hat_types
  ! use s2hat_pixelization
  use healpix_fft
  use misc_utils
  use alm_tools
  use s2hat_toolbox_mod
  
  implicit none

  ! include "mpif.h"

  !--------------------------------------------------------------------------

  public :: s2hat_alm2map, s2hat_alm2map_spin

contains
 
 
  ! ***************************************************************
  !            Wrapper for the computational routines
  ! ***************************************************************
 
  subroutine s2hat_alm2map( precompute_plms, pixelization, scan, nlmax, nmmax, nmvals, mvals, nmaps, nstokes, first_ring, last_ring, &
                           & map_size, local_map, lda, local_alm, nplm, local_plm, numprocs, myid, comm)

    implicit none
  
    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan
    integer(i4b), intent( in) :: nlmax, nmmax, nmvals, nstokes, nmaps, first_ring, last_ring, lda, map_size, numprocs, myid, comm, precompute_plms
    integer(i8b), intent(in) :: nplm
    integer(i4b), dimension(0:nmvals-1), intent(in) :: mvals
    complex(kind=dp), dimension(:,:,:,:), pointer :: local_alm
    real(kind=dp), dimension(:,:), pointer :: local_plm

    ! output
    real(kind=dp), dimension(:,:,:), pointer :: local_map

    ! internal
    integer(i4b) :: imap

    if (nstokes == 3) then
       ! Polarization routines
       if (precompute_plms == 1) then

          do imap = 1, nmaps
             call s2hat_alm2map_pol_pre1( pixelization, scan, nlmax, nmmax, nmvals, mvals, imap, first_ring, last_ring, map_size, local_map, lda, nplm, local_plm, local_alm, numprocs, myid, comm)
          enddo

       else if (precompute_plms == 2) then

          do imap = 1, nmaps
             call s2hat_alm2map_pol_pre2( pixelization, scan, nlmax, nmmax, nmvals, mvals, imap, first_ring, last_ring, map_size, local_map, lda, nplm, local_plm, local_alm, numprocs, myid, comm)
          enddo

       else
          call s2hat_alm2map_pol( pixelization, scan, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, map_size, local_map, lda, local_alm, numprocs, myid, comm)
       end if
    else
       ! Temperature routines
       if (precompute_plms == 1 .or. precompute_plms == 2) then

          do imap = 1, nmaps
             call s2hat_alm2map_sc_pre( pixelization, scan, nlmax, nmmax, nmvals, mvals, imap, first_ring, last_ring, map_size, local_map, lda, nplm, local_plm, local_alm, numprocs, myid, comm)
          enddo

       else
          call s2hat_alm2map_sc( pixelization, scan, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, map_size, local_map, lda, local_alm, numprocs, myid, comm)
       end if
    end if

 end subroutine s2hat_alm2map
 
 ! rearranges and redistributes over procs the data vector "data". On the input on each proc this
 ! vecotr stores results of the inverse theta transform for a set of m values assigned to the proc
 ! and for **all** pixel ring values.
 ! On the output the data vector contains the info for **all** m values and a subset of rings assigned
 ! to the proc. The m values are ordered in such a way that first go values of the proc 0, then 1, etc.
 
 subroutine redistInterProductsAlm2Map( nrings, nmvals, nmaps, ndata, data, myid, numprocs, comm)

   implicit none

   ! input/output
   integer(i4b), intent( in) :: nrings, nmvals, myid, nmaps, ndata, numprocs, comm
   complex(kind=dp), dimension(:,:), pointer :: data

   ! internal
   integer(i4b) :: i, ii, imap, j, jj, k, nringsall, ierr
   integer(i4b), dimension(:), allocatable :: send_counts, recv_counts, send_disps, recv_disps, ringsperproc, mperproc
   complex(kind=dp), dimension(:), pointer :: tmp

   call mpi_allreduce( nrings, nringsall, 1, mpi_integer, mpi_sum, comm, ierr)

   allocate( mperproc( 0:numprocs-1))
   call mpi_allgather( nmvals, 1, mpi_integer, mperproc, 1, mpi_integer, comm, ierr)

   allocate( ringsperproc( 0:numprocs-1))
   call mpi_allgather( nrings, 1, mpi_integer, ringsperproc, 1, mpi_integer, comm, ierr)

   allocate( send_counts( 0:numprocs-1))
   allocate( recv_counts( 0:numprocs-1))

   do i = 0, numprocs-1
      ringsperproc(i) = ringsperproc(i)
      recv_counts(i) = nrings*mperproc(i)
      send_counts(i) = nmvals*ringsperproc(i)
   enddo

   deallocate( mperproc)

   allocate( send_disps( 0:numprocs-1))
   allocate( recv_disps( 0:numprocs-1))

   recv_disps(0) = 0
   send_disps(0) = 0
      
   do i = 1, numprocs-1
      send_disps(i) = send_disps(i-1)+send_counts(i-1)
      recv_disps(i) = recv_disps(i-1)+recv_counts(i-1)
   enddo

   do imap = 1, nmaps

      allocate( tmp( 0:nmvals*nringsall-1))

      ! reshuffle to a ring-wise ordering here

      ii = 0
      do i = 0, nmvals-1
         jj = 0
         do j = 0, numprocs-1
            jj = jj+i*ringsperproc(j)
            do k = 0, ringsperproc(j)-1
               tmp( jj+k) = data( ii, imap)
               ii = ii+1
            end do
            jj = jj+(nmvals-i)*ringsperproc(j)
         end do
      end do

      call mpi_alltoallv( tmp, send_counts, send_disps, mpi_double_complex, data(0,imap), recv_counts, recv_disps, mpi_double_complex, comm, ierr)

     deallocate( tmp)

   end do ! over maps

   deallocate( send_counts, send_disps, recv_counts, recv_disps)
     
   deallocate( ringsperproc)

 end subroutine redistInterProductsAlm2Map

 ! actual alm2map transform 

 subroutine s2hat_alm2map_sc( pixelization, scan, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, map_size, local_map, lda, local_alms, numprocs, myid, comm)

    !=================================================================================
    !     computes the a(l,m) from a Temperature map for the HEALPIX pixelisation
    !     calculating spherical harmonics "on the fly"
    !=================================================================================

    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan
    integer (i4b) :: nlmax, nmmax, nmvals, nmaps, first_ring, last_ring, lda, map_size, numprocs, myid, comm
    integer (i4b), dimension(0:nmvals-1), intent(in) :: mvals
    complex (dp), dimension(:,:,:,:), pointer :: local_alms

    ! output
    real(dp), dimension(:,:,:), pointer :: local_map

    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !      for the Temperature field
    !     map(theta,phi) = sum_l_m a_lm Y_lm(theta,phi)
    !                    = sum_m {e^(i*m*phi) sum_l a_lm*lambda_lm(theta)}
    !
    !     where Y_lm(theta,phi) = lambda(theta) * e^(i*m*phi)
    !
    !     * the recurrence of Ylm is the standard one (cf Num Rec)
    !     * the sum over m is done by FFT
    !
    !=======================================================================

    integer(I4B) :: l, m, ith, scalel, scalem ! alm related
    integer(I8B) :: istart_south, istart_north, npix   ! map related
    integer(I4B) :: imap, iring, jring, kring, nrings, nphmx, nringsall, nringsobs

    integer(I4B) :: par_lm
    real(DP)     :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)     :: ovflow, unflow
    complex(DPC), dimension(-1:1,1:nmaps) :: b_ns
    real(DP),     dimension(:,:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: mfac
    real(DP),     dimension(:,:), allocatable :: recfac

    integer(i4b) :: i, ij, j, ll, ll1, l_min, mindx, mindx1, nlen
    complex(DPC), dimension(:,:), pointer :: b_north_ptr, b_south_ptr
    complex(DPC), dimension(:,:), allocatable, target :: b_north, b_south
    complex(DPC), dimension(:), allocatable :: phring
    real(DP),     dimension(:),   allocatable :: ring, cth, sth, kphi, kphi0
    integer(i4b) :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    integer(i4b), dimension(:), allocatable :: nph, nph0, indx
    integer(i4b), dimension(:), allocatable :: keep_north, keep_north_all, keep_south, keep_south_all, keep_it

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status, ierr

    complex(DPC), dimension(:,:,:), pointer :: local_alms_comp_a

    !=======================================================================

    if( lda == 1) then       ! i.e., healpix convention for the alm array

       local_alms_comp_a => local_alms(1,0:nlmax,0:nmvals-1,1:nmaps)

    else                     ! i.e., s2hat convention for the alm array

       local_alms_comp_a => local_alms(0:nlmax,0:nmvals-1,1,1:nmaps)

    end if

    nrings = (last_ring-first_ring+1)

    npix = pixelization%npixsall
    nphmx = pixelization%nphmx
    nringsall = pixelization%nringsall
    nringsobs = scan%nringsobs

    !     --- allocates space for arrays ---

    nchunks   = nringsall/SMAXCHK + 1  ! number of chunks
    chunksize = (nringsall+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    nlen = max( nringsobs*nmvals, nrings*(nmmax+1))
    allocate(b_north(0:nlen-1,1:nmaps),stat = status)
    call assert_alloc(status,code,'b_north')

    allocate(b_south(0:nlen-1,1:nmaps),stat = status)
    call assert_alloc(status,code,'b_south')

    allocate(nph(1:nringsall),stat = status)
    call assert_alloc(status,code,'nph')

    allocate(kphi(1:nringsall),stat = status)
    call assert_alloc(status,code,'kphi')

    allocate(cth(0:SMAXCHK),stat = status)
    call assert_alloc(status,code,'cth')
    allocate(sth(0:SMAXCHK),stat = status)
    call assert_alloc(status,code,'sth')

    allocate(keep_it(0:SMAXCHK),stat = status)
    call assert_alloc(status,code,'keep_it')

    allocate(keep_north_all(1:nringsall),stat = status)
    call assert_alloc(status,code,'keep_north_all')
    allocate(keep_south_all(1:nringsall),stat = status)
    call assert_alloc(status,code,'keep_south_all')

    allocate(recfac(0:1,0:nlmax),stat = status)
    call assert_alloc(status,code,'recfac')
	
    allocate(dalm(0:1,0:nlmax,1:nmaps), stat = status)
    call assert_alloc(status,code,'dalm')

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax,mfac)

    call init_rescale_local()   ! innitialize rescale table

    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    local_map = 0.d0     ! set the whole map to zero
    b_north(:,:) = 0_dpc ! pad with zeros
    b_south(:,:) = 0_dpc

    ! loop on chunks
    kring = 0
    do ichunk = 0, nchunks-1

       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nringsall)

       do ith = lchk, uchk
          ithl = ith - lchk !local index

          ! extract pixel location information
          kphi(ith) = pixelization%kphi(ith)
          nph(ith) = pixelization%nph(ith)
          cth(ithl) = pixelization%cth(ith)
          sth(ithl) = pixelization%sth(ith)

          keep_north_all( ith) = scan%nfl(ith)
          keep_south_all( ith) = scan%sfl(ith)
          keep_it( ithl) = scan%fl(ith)		  
       enddo

       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)

       do mindx = 0, nmvals-1
	   
          mindx1 = mindx+1
	  m = mvals( mindx)
		  
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)

          ! extract needed alm under memory and CPU efficient form
          do imap = 1, nmaps
             do ll = m, nlmax
                ll1 = ll+1
                dalm(0,ll,imap) = real( local_alms_comp_a(ll1,mindx1,imap),kind=dp)
                dalm(1,ll,imap) = aimag( local_alms_comp_a(ll1,mindx1,imap))
             enddo
          enddo

          jring = kring
          do ithl = 0, uchk-lchk

             if( keep_it( ithl) == 1) then            ! do not calculate rings which are not observed 
                l_min = l_min_ylm(m, sth(ithl))
                if (nlmax >= l_min) then ! skip calculations when Ylm too small

                   ! determine lam_mm
                   call mpi_compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)
               
                   !           ---------- l = m ----------

                   par_lm = 1  ! = (-1)^(l+m)

                   if (m >= l_min) then
                      lam_lm = corfac * lam_mm

                      b_ns( 1, 1:nmaps) = cmplx(lam_lm * dalm(0,m, 1:nmaps), lam_lm * dalm(1,m, 1:nmaps), kind=DP)
                      b_ns(-1, 1:nmaps) = 0.0_dpc
                   else
                      b_ns = 0.0_dpc
                   endif

                   !           ---------- l > m ----------

                   lam_0 = 0.0_dp
                   lam_1 = 1.0_dp
                   scalel=0
                   cth_ring = cth(ithl)
                   lam_2 = cth_ring * lam_1 * recfac(0,m)

                   do l = m+1, nlmax
                      par_lm = - par_lm  ! = (-1)^(l+m)
                      if (l >= l_min) then
                         lam_lm = lam_2 * corfac * lam_mm                   
                         b_ns(par_lm, 1:nmaps)  = b_ns(par_lm, 1:nmaps) &
                                    &     + cmplx(lam_lm * dalm(0,l, 1:nmaps), lam_lm * dalm(1,l, 1:nmaps), kind=DP)  ! increment even or odd
                      
                      endif

                      lam_0 = lam_1 * recfac(1,l-1)
                      lam_1 = lam_2
                      lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)

                      if (abs(lam_2) > OVFLOW) then
                         lam_1 = lam_1*UNFLOW
                         lam_2 = lam_2*UNFLOW
                         scalel= scalel + 1
                         corfac = rescale_tab(max(scalel+scalem,RSMIN))
                      elseif (abs(lam_2) < UNFLOW) then
                         lam_1 = lam_1*OVFLOW
                         lam_2 = lam_2*OVFLOW
                         scalel= scalel - 1
                         corfac = rescale_tab(max(scalel+scalem,RSMIN))
                      endif
                   
                   enddo ! loop on l

                   ith = mindx*nringsobs+jring   ! global index 

                   b_north(ith,1:nmaps) = b_ns(1,1:nmaps) + b_ns(-1,1:nmaps)
                   b_south(ith,1:nmaps) = b_ns(1,1:nmaps) - b_ns(-1,1:nmaps)

                endif ! test on nlmax
			 
                jring = jring +1

	     endif ! if ring observed
			 
          enddo ! loop over chunk rings (ithl)

       enddo ! loop over m

       kring = jring

    enddo    ! loop over chunks

    deallocate( cth, sth)
    deallocate( recfac, dalm)
    deallocate( mfac)
    deallocate( keep_it)

    ! and now the phi integrals ...

    allocate( nph0(0:nrings-1))	
    nph0 = nph(first_ring:last_ring)
    deallocate( nph)
    allocate( kphi0(0:nrings-1))	
    kphi0 = kphi(first_ring:last_ring)
    deallocate( kphi)

    allocate( keep_north(0:nrings-1))	
    keep_north = keep_north_all(first_ring:last_ring)
    deallocate( keep_north_all)

    allocate( keep_south(0:nrings-1))	
    keep_south = keep_south_all(first_ring:last_ring)
    deallocate( keep_south_all)
		
    !! the redistribution

    b_north_ptr => b_north
    b_south_ptr => b_south

    call redistInterProductsAlm2Map( nrings, nmvals, nmaps, nlen, b_north_ptr, myid, numprocs, comm)
    call redistInterProductsAlm2Map( nrings, nmvals, nmaps, nlen, b_south_ptr, myid, numprocs, comm)

    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')

    allocate(phring(0:nmmax),stat = status)
    call assert_alloc(status,code,'phring')

    allocate(indx(0:nmmax), stat=status)
    call assert_alloc(status,code,'indx')

    call fix_mmap( 1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

    do imap = 1, nmaps

       istart_north = 0
       istart_south = map_size

       do iring = 0, nrings-1

          nphl = nph0(iring)

          if( keep_north( iring) == 1) then
	   
             !        ---------------------------------------------------------------
             !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
	     !        ---------------------------------------------------------------

	     do i = 0, nmmax
	        ij = iring+indx(i)*nrings
                phring(i) = b_north( ij, imap)
	     enddo

	     call s2hat_ring_synthesis(nlmax,nmmax,phring,nphl,ring,kphi0(iring))   ! north hemisph. + equator
	     local_map(istart_north:istart_north+nphl-1,1,imap) = ring(0:nphl-1)
		  
          endif

          istart_north = istart_north+nphl

          istart_south = istart_south-nphl
	   
          if( keep_south(iring) == 1) then

	     do i = 0, nmmax
                ij = iring+indx(i)*nrings
                phring(i) = b_south(ij,imap)
   	     enddo

             call s2hat_ring_synthesis(nlmax,nmmax,phring,nphl,ring,kphi0(iring)) ! south hemisph. w/o equat
             local_map(istart_south:istart_south+nphl-1,1,imap) = ring(0:nphl-1)

          endif

       enddo    ! loop over rings
    enddo       ! loop over maps

    !     --------------------
    !     free memory and exit
    !     --------------------
     
    deallocate( b_north, b_south)
    deallocate( phring)
    deallocate( ring)
    deallocate( kphi0)
    deallocate( nph0)
    deallocate( indx)
    deallocate( keep_north, keep_south)

    return

  end subroutine s2hat_alm2map_sc
 
  subroutine s2hat_alm2map_sc_pre( pixelization, scan, nlmax, nmmax, nmvals, mvals, imap, first_ring, last_ring, map_size, local_map, &
                                 & lda, nplm, local_plm, local_alms, numprocs, myid, comm)

    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !     for the Temperature field, with precomputed and stored Ylm
	!
    !     map(theta,phi) = sum_l_m a_lm Y_lm(theta,phi)
    !                    = sum_m {e^(i*m*phi) sum_l a_lm*lambda_lm(theta)}
    !
    !     where Y_lm(theta,phi) = lambda(theta) * e^(i*m*phi)
    !
    !     * the recurrence of Ylm is the standard one (cf Num Rec)
    !     * the sum over m is done by FFT
    !
    !=======================================================================

    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan
    integer (i4b) :: nlmax, nmmax, nmvals, imap, first_ring, last_ring, lda, map_size, numprocs, myid, comm
    integer (i4b), dimension(0:nmvals-1) :: mvals
    integer (i8b), intent(in) :: nplm
    real(kind=dp), dimension(:,:), pointer :: local_plm
    complex (dp), dimension(:,:,:,:), pointer :: local_alms
	
    ! output
    real(dp), dimension(:,:,:), pointer :: local_map

    integer(I4B) :: l, m, ith, scalel, scalem          ! alm related
    integer(I8B) :: istart_south, istart_north, npix   ! map related
    integer(I4B) :: iring, jring, kring, nrings, nphmx, nringsall, nringsobs

    integer(I4B)  :: par_lm
    real(DP)      :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)      :: ovflow, unflow
    complex(DPC), dimension(-1:1) :: b_ns
    real(DP),     dimension(:,:), allocatable :: dalm

    integer(i4b)  :: i, ij, j, ll, ll1, l_min, l_start, mindx, mindx1, nlen, mshift
    complex(DPC), dimension(:,:), pointer :: b_north_ptr, b_south_ptr
    complex(DPC), dimension(:,:), allocatable, target :: b_north, b_south
    complex(DPC), dimension(:), allocatable :: phring
    real(DP),     dimension(:), allocatable :: ring, cth, sth, kphi0, kphi
    integer(i4b) :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    integer(I8B) :: n_lm, i_mm, i_mm1
    integer(i4b), dimension(:), allocatable :: nph, nph0, indx
    integer(i4b), dimension(:), allocatable :: keep_north, keep_north_all, keep_south, keep_south_all, keep_it

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status, ierr

    complex(DPC), dimension(:,:), pointer :: local_alms_comp_a
    real(DP), dimension(:), pointer :: local_plm_comp_a

    !=======================================================================

    if( lda == 1) then       ! i.e., healpix convention for the alm array

       local_alms_comp_a => local_alms(1,0:nlmax,0:nmvals-1,imap)
       local_plm_comp_a => local_plm(1,0:nplm-1)

    else                     ! i.e., s2hat convention for the alm array

       local_alms_comp_a => local_alms(0:nlmax,0:nmvals-1,1,imap)
       local_plm_comp_a => local_plm(0:nplm-1,1)

    end if

    npix = pixelization%npixsall
    nphmx = pixelization%nphmx
    nringsall = pixelization%nringsall
    nringsobs = scan%nringsobs

    nrings = (last_ring-first_ring+1)

    n_lm = nummmodes( nlmax, nmvals, mvals) ! an actual size of the local_plm array 	

    !     --- allocates space for arrays ---

    nchunks   = nringsall/SMAXCHK + 1         ! number of chunks
    chunksize = (nringsall+nchunks-1)/nchunks ! <= SMAXCHK

    nlen = max( nringsobs*nmvals, nrings*(nmmax+1))
    allocate(b_north(0:nlen-1,1:1),stat = status)
    call assert_alloc(status,code,'b_north')

    allocate(b_south(0:nlen-1,1:1),stat = status)
    call assert_alloc(status,code,'b_south')

    allocate(nph(1:nringsall),stat = status)
    call assert_alloc(status,code,'nph')

    allocate(kphi(1:nringsall),stat = status)
    call assert_alloc(status,code,'kphi')

    allocate(cth(0:SMAXCHK),stat = status)
    call assert_alloc(status,code,'cth')
    allocate(sth(0:SMAXCHK),stat = status)
    call assert_alloc(status,code,'sth')

    allocate(keep_it(0:SMAXCHK),stat = status)
    call assert_alloc(status,code,'keep_it')

    allocate(keep_north_all(1:nringsall),stat = status)
    call assert_alloc(status,code,'keep_north_all')

    allocate(keep_south_all(1:nringsall),stat = status)
    call assert_alloc(status,code,'keep_south_all')

    allocate(dalm(0:1,0:nlmax), stat = status)
    call assert_alloc(status,code,'dalm')

    !     ------------ initiate variables and arrays ----------------

    b_north(:,1:1) = 0_dpc          ! pad with zeros
    b_south(:,1:1) = 0_dpc

    ! loop over chunks
    kring = 0
    do ichunk = 0, nchunks-1

       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nringsall)

       do ith = lchk, uchk
          ithl = ith - lchk !local index

          ! extract pixel location information
          kphi(ith) = pixelization%kphi(ith)
          nph(ith) = pixelization%nph(ith)
          cth(ithl) = pixelization%cth(ith)
          sth(ithl) = pixelization%sth(ith)

          keep_north_all( ith) = scan%nfl(ith)
          keep_south_all( ith) = scan%sfl(ith)
          keep_it( ithl) = scan%fl(ith)		  
       enddo

       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)

       mshift = 0

       do mindx = 0, nmvals-1

          mindx1 = mindx+1
          m = mvals( mindx)
		  
          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             ll1 = ll+1
             dalm(0,ll) = real(local_alms_comp_a(ll1,mindx1),kind=dp)
             dalm(1,ll) = aimag(local_alms_comp_a(ll1,mindx1))
          enddo

          jring = kring		  
          do ithl = 0, uchk-lchk
		  
             if( keep_it( ithl) == 1) then           ! do not calculate rings which are not observed 

                l_min = l_min_ylm(m, sth(ithl))

                if (nlmax >= l_min) then        ! skip calculations when Ylm too small

                   i_mm = n_lm * jring + mshift ! location of Ym,m for ring ith
                   i_mm1 = i_mm+1

                   !           ---------- l = m ----------

		   par_lm = 1

                   if (m >= l_min) then
                      b_ns( 1) = cmplx(local_plm_comp_a(i_mm1) * dalm(0,m), local_plm_comp_a(i_mm1) * dalm(1,m), kind=DP)
                      b_ns(-1) = 0.0_dpc
                   else
                      b_ns = 0.0_dpc
                   endif

                   !           ---------- l > m ----------

                   l_start = max(m+1, l_min)
                   if (mod(m+l_start,2) == 0) l_start = l_start+1
                   do l = l_start, nlmax, 2
                      b_ns(-1) = b_ns(-1) + &
                               & cmplx(local_plm_comp_a(i_mm1+l-m) * dalm(0,l), local_plm_comp_a(i_mm1+l-m)*dalm(1,l), kind=DP)
                   enddo ! loop over l

                   l_start = max(m+2, l_min)
                   if (mod(m+l_start,2) == 1) l_start = l_start+1
                   do l = l_start, nlmax, 2
                      b_ns(1)  = b_ns(1) + &
                               & cmplx(local_plm_comp_a(i_mm1+l-m) * dalm(0,l), local_plm_comp_a(i_mm1+l-m)*dalm(1,l),kind=DP) 
                   enddo ! loop over l

                   iring = mindx*nringsobs+jring

                   b_north( iring,1) = b_ns(1) + b_ns(-1)
                   b_south( iring,1) = b_ns(1) - b_ns(-1)
				
                endif ! test on nlmax

                jring = jring + 1

             endif ! if the ring observed

          enddo ! loop over chunk rings (ithl)

	  mshift = mshift + nlmax-m+1		  

       enddo ! loop over m

       kring = jring

    enddo ! loop over chunks

    deallocate( cth, sth)
    deallocate( dalm)
    deallocate( keep_it)

    allocate( kphi0(0:nrings-1), stat=status)
    call assert_alloc(status,code,'kphi0')
    kphi0 = kphi( first_ring:last_ring)
    deallocate( kphi)
    
    allocate( nph0(0:nrings-1), stat=status)
    nph0 = nph(first_ring:last_ring)
    deallocate( nph)

    allocate( keep_north(0:nrings-1), stat=status)
    keep_north = keep_north_all(first_ring:last_ring)
    deallocate( keep_north_all)
    allocate( keep_south(0:nrings-1), stat=status)
    keep_south = keep_south_all(first_ring:last_ring)
    deallocate( keep_south_all)

    !! the redistribution
	
    b_north_ptr => b_north
    b_south_ptr => b_south

    call redistInterProductsAlm2Map( nrings, nmvals, 1, nlen, b_north_ptr, myid, numprocs, comm)
    call redistInterProductsAlm2Map( nrings, nmvals, 1, nlen, b_south_ptr, myid, numprocs, comm)	
       
    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')

    allocate(phring(0:nmmax),stat = status)
    call assert_alloc(status,code,'phring')

    allocate(indx(0:nmmax), stat=status)
    call assert_alloc(status,code,'indx')

    call fix_mmap( -1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

    local_map(:,:,imap) = 0.d0 ! set the whole map to zero
    istart_north = 0
    istart_south = map_size

    do iring = 0, nrings-1
       
       nphl = nph0(iring)

       if( keep_north(iring) == 1) then

          !        ---------------------------------------------------------------
          !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
          !        ---------------------------------------------------------------
             
          do j = 0, nmmax
             ij = indx(j)
             phring(ij) = b_north( iring+j*nrings,1)
          enddo
          call s2hat_ring_synthesis(nlmax,nmmax,phring,nphl,ring,kphi0(iring))   ! north hemisph. + equator
          local_map( istart_north:istart_north+nphl-1, 1, imap) = ring(0:nphl-1)
          
       endif

       istart_north = istart_north+nphl

       istart_south = istart_south-nphl
											                                  
       if( keep_south(iring) == 1) then

	  do j = 0, nmmax
             ij = indx(j)
             phring(ij) = b_south(iring+j*nrings,1)
	  enddo

	  call s2hat_ring_synthesis(nlmax,nmmax,phring,nphl,ring,kphi0(iring)) ! south hemisph. w/o equat
          local_map(istart_south:istart_south+nphl-1,1,imap) = ring(0:nphl-1)

       endif

    enddo    ! loop over rings

    !     --------------------
    !     free memory and exit
    !     --------------------
    
    deallocate( ring)
    deallocate( phring)
    deallocate( b_north, b_south)
    deallocate( kphi0)
    deallocate( nph0)
    deallocate( indx)
    deallocate( keep_south, keep_north)

    return

  end subroutine s2hat_alm2map_sc_pre

  subroutine s2hat_alm2map_spin( pixelization, scan, spin_in, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, map_size, local_map, lda, local_alms, numprocs, myid, comm)

    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !     for the Temperature+Polarization field and spherical harmonics
    !     calculated "on the fly"
    !=======================================================================

    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan
    integer (i4b), intent(in) :: spin_in, nlmax, nmmax, nmvals, nmaps, first_ring, last_ring, lda, map_size, numprocs, myid, comm
    integer (i4b), dimension(0:nmvals-1), intent(in) :: mvals
    complex (dp), dimension(:,:,:,:), pointer :: local_alms

    ! output
    real(dp), dimension(:,:,:), pointer :: local_map

    ! internal
    integer(I4B) :: imap, l, ll1, m, ith, scalel, scalem    ! alm related
    integer(I4B) :: istart_south, istart_north, l_start     ! map related
    integer(I4B) :: nrings, nringsall, nringsobs, npix, nphmx, mm, ss
    integer(i4b), dimension(1:2) :: compSign

    real(DP),     dimension(0:15,1:nmaps) :: b_ns
    real(DP),     dimension(0:7) :: tmp
    real(DP),     dimension(:,:), allocatable :: recfac, slam_lm
    real(DP),     dimension(:,:,:), allocatable :: dalm
    real(DP),     dimension(:), allocatable :: mfac, sfac

    integer(i4b)  :: ll, l_min, mindx, mindx1, nlen
    complex(DPC), dimension(:,:,:), allocatable, target :: b_TQU
    complex(DPC), dimension(:,:), pointer :: b_TQU_ptr
    complex(DPC), dimension(:,:), allocatable, target :: column
    complex(DPC), dimension(:), allocatable :: phring
    real(DP),     dimension(:), allocatable :: ring
    integer(i4b)  :: spin, ithl, jth, nphl, i, j, k, ij, iring, ith1
    integer(I4B), dimension(:), allocatable:: nph0, nph, indx
    real(DP),     dimension(:), allocatable :: cth, sth, kphi, kphi0
    real(DP) :: spinPow, signFix, mfac0
    integer(i4b), dimension(:), allocatable :: keep_north, keep_north_all, keep_south, keep_south_all, keep_it

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status, ierr

    complex(DPC), dimension(:,:,:), pointer :: local_alms_comp_a, local_alms_comp_b

    !=======================================================================

    spin = abs( spin_in)

    spinPow = (-1.d0) ** spin

   ! define how the input real maps are to be combined into a complex field of |s| (>0) spin

    compSign(1) = 1
    compSign(2) = -SPIN_CONV_SIGN

    if( lda == 2) then      ! i.e., healpix convention for the alm array

       local_alms_comp_a => local_alms(1,0:nlmax,0:nmvals-1,1:nmaps)
       local_alms_comp_b => local_alms(2,0:nlmax,0:nmvals-1,1:nmaps)

    else                    ! i.e., s2hat convention for the alm array

      local_alms_comp_a => local_alms(0:nlmax,0:nmvals-1,1,1:nmaps)
      local_alms_comp_b => local_alms(0:nlmax,0:nmvals-1,2,1:nmaps)

    end if

    nrings = (last_ring-first_ring+1)

    npix = pixelization%npixsall
    nphmx = pixelization%nphmx
    nringsall = pixelization%nringsall
    nringsobs = scan%nringsobs

    !     --- allocates space for arrays ---

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(sfac(0:spin),stat = status)
    call assert_alloc(status,code,'sfac')

    nlen = max( (nmmax+1)*nrings, nmvals*nringsobs)
    allocate(b_TQU(0:nlen-1, 1:4, 1:nmaps),stat = status)
    call assert_alloc(status,code,'b_TQU')

    allocate(recfac(0:1,0:nlmax), dalm(0:3,0:nlmax,1:nmaps),stat = status)
    call assert_alloc(status,code,'recfac, dalm')

    allocate( slam_lm(1:2,0:nlmax), stat=status)
    call assert_alloc(status,code,'lam_lm')  
    
    allocate(nph(1:nringsall), stat=status)
    call assert_alloc(status,code,'nph')

    allocate(kphi(1:nringsall), stat=status)
    call assert_alloc(status,code,'kphi')

    allocate(cth(0:nringsall-1), stat=status)
    call assert_alloc(status,code,'cth')

    allocate(sth(0:nringsall-1), stat=status)
    call assert_alloc(status,code,'sth')

    allocate(keep_it(0:nringsall-1), stat=status)
    call assert_alloc(status,code,'keep_it')
    allocate(keep_north_all(1:nringsall), stat=status)
    call assert_alloc(status,code,'keep_north_all')
    allocate(keep_south_all(1:nringsall), stat=status)
    call assert_alloc(status,code,'keep_south_all')

    !     ------------ initiate variables and arrays ----------------

    ! generate Polarization normalisation, etc
    call s2hat_gen_mfac(spin,nmmax,mfac)

    ! generate Polarization normalisation, etc
    call s2hat_gen_sfac( spin, sfac)

    call init_rescale_local()   ! innitialize rescale table

    b_TQU(:,:,:) = 0_dpc ! pad with zeros

    do ith = 0, nringsall-1
       ith1 = ith + 1
       ! extract pixel location information
       kphi(ith1) = pixelization%kphi(ith1)
       nph(ith1) = pixelization%nph(ith1)
       cth(ith) = pixelization%cth(ith1)
       sth(ith) = pixelization%sth(ith1)

       keep_north_all( ith1) = scan%nfl(ith1)
       keep_south_all( ith1) = scan%sfl(ith1)
       keep_it( ith) = scan%fl(ith1)
    enddo

       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------

    iring = 0
    do mindx = 0, nmvals-1
             
       m = mvals( mindx)

       if( m >= spin) then

          signFix = 1.d0
          mm = m
          ss = spin

          mfac0 = mfac(m)

       else

          signFix = (-1.d0)**(spin+m)
          mm = spin ! swap
          ss = m

          mfac0 = sfac(m)

       end if

       l_start = mm

       ! extract needed alm in a memory and CPU efficient form

       mindx1 = mindx+1

       do imap = 1, nmaps
          do ll = mm, nlmax
             ll1 = ll+1

             dalm(0,ll,imap) =  real(local_alms_comp_a(ll1,mindx1,imap),kind=dp) ! 1st spin component
             dalm(1,ll,imap) = aimag(local_alms_comp_a(ll1,mindx1,imap))
             dalm(2,ll,imap) =  real(local_alms_comp_b(ll1,mindx1,imap),kind=dp) ! 2nd spin component
             dalm(3,ll,imap) = aimag(local_alms_comp_b(ll1,mindx1,imap))
           enddo
        enddo

       ! generate recursion factors (recfac) for Ylm of degree mm and spin ss
       call s2hat_gen_recfac(ss, nlmax, mm, recfac)

       jth = 0
       do ith = 0, nringsall-1

          if( keep_it( ith) == 1) then 

              !! l_min =  l_min_ylm(m, sth(ith)) ! lowest l of non-negligible Ylm

              l_min = mm ! no attempt to save flops here yet ...

              if (nlmax >= l_min) then

              ! compute ss_lam_lmm(p,theta) for all l>=mm

              call mpi_do_slam_lm(spin, nlmax, m, cth(ith), sth(ith), mfac0, recfac, slam_lm)

              ! synthetize real spin component maps 
    
              b_ns = 0.0_dp
 
              do imap = 1, nmaps
   
                 do l = l_start, nlmax-1, 2

                    ! l+m of parity as that of l_start + m

                    tmp(0:3) = slam_lm(1,l) * dalm(0:3,l,imap)
                    tmp(4:7) = slam_lm(2,l) * dalm(0:3,l,imap)

                    b_ns(0:7,imap) = b_ns(0:7,imap) + tmp(0:7) ! (0:1) -> F^+ * s_E_lm; (2:3) ->  F^+ * s_B_lm
                                                               ! (4:5) -> F^- * s_E_lm; (6:7) ->  F^- * s_B_lm

                    ! l+m of parity different than that of l_start + m

                    tmp(0:3) = slam_lm(1,l+1) * dalm(0:3,l+1,imap)
                    tmp(4:7) = slam_lm(2,l+1) * dalm(0:3,l+1,imap)

                    b_ns(8:15,imap) = b_ns(8:15,imap) + tmp(0:7)

                 enddo

                 if( mod( nlmax-l_start, 2) == 0) then

                     tmp(0:3) = slam_lm(1, nlmax) * dalm(0:3, nlmax, imap)
                     tmp(4:7) = slam_lm(2, nlmax) * dalm(0:3, nlmax, imap)

                     b_ns(0:7,imap) = b_ns(0:7,imap) + tmp(0:7) ! (0:1) -> F^+ * s_E_lm; (2:3) ->  F^+ * s_B_lm
                                                                ! (4:5) -> F^- * s_E_lm; (6:7) ->  F^- * s_B_lm

                  end if

               enddo

               iring = mindx*nringsobs+jth
				               
               b_TQU(iring,1,1:nmaps) = 0.5*SPIN_CONV_SIGN * cmplx(b_ns(0,1:nmaps) + b_ns( 8, 1:nmaps) - b_ns(7,1:nmaps) - b_ns(15,1:nmaps), &
            &                                                     b_ns(1,1:nmaps) + b_ns( 9, 1:nmaps) + b_ns(6,1:nmaps) + b_ns(14,1:nmaps), kind = DP) ! 1st component, north
               b_TQU(iring,2,1:nmaps) = 0.5*SPIN_CONV_SIGN * cmplx(b_ns(2,1:nmaps) + b_ns(10,1:nmaps) + b_ns(5,1:nmaps) + b_ns(13,1:nmaps), &
            &                                                     b_ns(3,1:nmaps) + b_ns(11,1:nmaps) - b_ns(4,1:nmaps) - b_ns(12,1:nmaps), kind = DP) ! 2nd component, north
               b_TQU(iring,3,1:nmaps) = 0.5*SPIN_CONV_SIGN * signFix * spinPow * cmplx(b_ns(0,1:nmaps) - b_ns(8,1:nmaps) + b_ns(7,1:nmaps) - b_ns(15,1:nmaps), &
            &                                                     b_ns(1,1:nmaps) - b_ns(9,1:nmaps) - b_ns(6,1:nmaps) + b_ns(14,1:nmaps), kind = DP) ! 1st component, south
               b_TQU(iring,4,1:nmaps) = 0.5*SPIN_CONV_SIGN * signFix * spinPow * cmplx(b_ns(2,1:nmaps) - b_ns(10,1:nmaps) - b_ns( 5, 1:nmaps) + b_ns( 13, 1:nmaps),  &
            &                                                     b_ns(3,1:nmaps) - b_ns(11,1:nmaps) + b_ns(4,1:nmaps) - b_ns(12,1:nmaps), kind = DP) ! 2nd component, south

            endif ! test on nlmax

            jth = jth + 1

         endif ! test if ring is observed

      enddo ! loop over rings (ith)

    enddo ! loop over m

    deallocate( cth, sth)
    deallocate( recfac, dalm, slam_lm)
    deallocate( mfac, sfac)
    deallocate( keep_it)

    allocate( kphi0(0:nrings-1), stat=status)
    call assert_alloc(status,code,'kphi0')
    kphi0 = kphi(first_ring:last_ring)
    deallocate( kphi)
    
    allocate( nph0(0:nrings-1), stat=status)
    nph0 = nph(first_ring:last_ring)
    deallocate( nph)
	
    allocate( keep_north(0:nrings-1), stat=status)
    keep_north = keep_north_all(first_ring:last_ring)
    deallocate( keep_north_all)
    allocate( keep_south(0:nrings-1), stat=status)
    keep_south = keep_south_all(first_ring:last_ring)
    deallocate( keep_south_all)

    !! redistribute

    allocate(indx(0:nmmax), stat=status)
    call assert_alloc(status,code,'indx')

    call fix_mmap( -1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

    allocate(column(0:nlen-1,1:nmaps),stat = status)
    call assert_alloc(status,code,'column')
    
    do i = 1, 4
       column = b_TQU(0:nlen-1,i,1:nmaps)
       b_TQU_ptr => column 
       call redistInterProductsAlm2Map( nrings, nmvals, nmaps, nlen, b_TQU_ptr, myid, numprocs, comm)

       do imap = 1, nmaps
          do j = 0, nmmax
             do k = 0, nrings-1
                b_TQU( (nmmax+1)*k+indx(j),i,imap) = b_TQU_ptr(k+j*nrings,imap)
             enddo
          enddo
       enddo
    enddo

    deallocate( column)

    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')

    allocate(phring(0:nmmax),stat = status)
    call assert_alloc(status,code,'phring')

    local_map = 0.d0 ! set the whole map to zero


    do imap = 1, nmaps

       do i=1,2

          istart_north = 0

          do iring = 0, nrings-1

             nphl = nph0(iring)

             if( keep_north( iring) == 1) then
		         
                !        ---------------------------------------------------------------
                !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
                !        ---------------------------------------------------------------
 
                phring(0:nmmax) = b_TQU(iring*(nmmax+1):(iring+1)*(nmmax+1)-1,i,imap)

                call s2hat_ring_synthesis(nlmax,nmmax,phring,nphl,ring,kphi0(iring))   ! north hemisph. + equator
                local_map(istart_north:istart_north+nphl-1, i, imap) = ring(0:nphl-1)*compSign(i)
			 
	     endif

             istart_north = istart_north+nphl

	   enddo
  
	enddo                                       

	do i = 1, 2

           istart_south = map_size

           do iring = 0, nrings-1

              nphl = nph0(iring)
              istart_south = istart_south-nphl

              if( keep_south( iring) == 1) then

                  phring(0:nmmax) = b_TQU(iring*(nmmax+1):(iring+1)*(nmmax+1)-1,2+i,imap)

                  call s2hat_ring_synthesis(nlmax,nmmax,phring,nphl,ring,kphi0(iring)) ! south hemisph. w/o equat
                  local_map(istart_south:istart_south+nphl-1, i, imap) = ring(0:nphl-1)*compSign(i)
 
              endif
           enddo   ! loop over rings

        enddo  

    enddo  ! over maps

    !     --------------------
    !     free memory and exit
    !     --------------------
    
    deallocate( ring)
    deallocate( phring)
    deallocate( b_TQU)
    deallocate( kphi0)
    deallocate( nph0)
    deallocate( indx)
    deallocate( keep_south, keep_north)

    return

  end subroutine s2hat_alm2map_spin

  subroutine s2hat_alm2map_pol( pixelization, scan, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, map_size, local_map, lda, local_alms, numprocs, myid, comm)

    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !     for the Temperature+Polarization field and spherical harmonics
	!     calculated "on the fly"
    !=======================================================================

    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan
    integer (i4b), intent(in) :: nlmax, nmmax, nmvals, nmaps, first_ring, last_ring, lda, map_size, numprocs, myid, comm
    integer (i4b), dimension(0:nmvals-1), intent(in) :: mvals
    complex (dp), dimension(:,:,:,:), pointer :: local_alms

    ! output
    real(dp), dimension(:,:,:), pointer :: local_map

    ! internal
    integer(I4B) :: imap, l, m, ith, scalel, scalem    ! alm related
    integer(I4B) :: istart_south, istart_north   ! map related
    integer(I4B) :: nrings, nringsall, nringsobs, npix, nphmx

    integer(I4B)  :: par_lm, l_start
    real(DP)      :: lam_mm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)      :: lambda_w, lambda_x, normal_m, lam_lm1m
    real(DP)      :: fm2, fl, flm1
    real(DP)      :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP)      :: ovflow, unflow, c_on_s2, one_on_s2
    real(DP),     dimension(-3:8,1:nmaps) :: b_ns
    real(DP),     dimension(:,:), allocatable :: recfac, lam_lm
    real(DP),     dimension(:,:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: lam_fact
    real(DP),     dimension(:), allocatable :: mfac

    integer(i4b)  :: ll, ll1, l_min, mindx, mindx1, nlen
    real(DP),     dimension(:), allocatable :: normal_l
    complex(DPC), dimension(:,:,:), allocatable, target :: b_TQU
    complex(DPC), dimension(:,:), pointer :: b_TQU_ptr
    complex(DPC), dimension(:,:), allocatable, target :: column
    complex(DPC), dimension(:), allocatable :: phring
    real(DP),     dimension(:), allocatable :: ring
    integer(i4b)  :: ithl, jth, nphl, i, j, k, ij, iring, ith1
    integer(I4B), dimension(:), allocatable:: nph0, nph, indx
    real(DP),     dimension(:), allocatable :: cth, sth, kphi, kphi0
    integer(i4b), dimension(:), allocatable :: keep_north, keep_north_all, keep_south, keep_south_all, keep_it

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status, ierr

    complex(DPC), dimension(:,:,:), pointer :: local_alms_comp_a, local_alms_comp_b, local_alms_comp_c

    !=======================================================================


    if( lda == 3) then       ! i.e., healpix convention for the alm array

       local_alms_comp_a => local_alms(1,0:nlmax,0:nmvals-1,1:nmaps)
       local_alms_comp_b => local_alms(2,0:nlmax,0:nmvals-1,1:nmaps)
       local_alms_comp_c => local_alms(3,0:nlmax,0:nmvals-1,1:nmaps)

    else                     ! i.e., s2hat convention for the alm array

       local_alms_comp_a => local_alms(0:nlmax,0:nmvals-1,1,1:nmaps)
       local_alms_comp_b => local_alms(0:nlmax,0:nmvals-1,2,1:nmaps)
       local_alms_comp_c => local_alms(0:nlmax,0:nmvals-1,3,1:nmaps)

    end if

    nrings = (last_ring-first_ring+1)

    npix = pixelization%npixsall
    nphmx = pixelization%nphmx
    nringsall = pixelization%nringsall
    nringsobs = scan%nringsobs

    !     --- allocates space for arrays ---

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(normal_l(0:nlmax),stat = status)
    call assert_alloc(status,code,'normal_l')

    nlen = max( (nmmax+1)*nrings, nmvals*nringsobs)
    allocate(b_TQU(0:nlen-1, 1:6, 1:nmaps),stat = status)
    call assert_alloc(status,code,'b_TQU')

    allocate(recfac(0:1,0:nlmax), dalm(0:5,0:nlmax,1:nmaps), lam_fact(0:nlmax),stat = status)
    call assert_alloc(status,code,'recfac, dalm & lam_fact')

    allocate( lam_lm(1:3,0:nlmax), stat=status)
    call assert_alloc(status,code,'lam_lm')  
    
    allocate(nph(1:nringsall), stat=status)
    call assert_alloc(status,code,'nph')

    allocate(kphi(1:nringsall), stat=status)
    call assert_alloc(status,code,'kphi')

    allocate(cth(0:nringsall-1), stat=status)
    call assert_alloc(status,code,'cth')

    allocate(sth(0:nringsall-1), stat=status)
    call assert_alloc(status,code,'sth')

    allocate(keep_it(0:nringsall-1), stat=status)
    call assert_alloc(status,code,'keep_it')
    allocate(keep_north_all(1:nringsall), stat=status)
    call assert_alloc(status,code,'keep_north_all')
    allocate(keep_south_all(1:nringsall), stat=status)
    call assert_alloc(status,code,'keep_south_all')

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax,mfac)
    
    call init_rescale_local()   ! innitialize rescale table

    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)

    ! generate Polarization normalisation
    call gen_normpol(nlmax, normal_l)

    b_TQU(:,:,:) = 0_dpc ! pad with zeros

    do ith = 0, nringsall-1
       ith1 = ith + 1
       ! extract pixel location information
       kphi(ith1) = pixelization%kphi(ith1)
       nph(ith1) = pixelization%nph(ith1)
       cth(ith) = pixelization%cth(ith1)
       sth(ith) = pixelization%sth(ith1)

       keep_north_all( ith1) = scan%nfl(ith1)
       keep_south_all( ith1) = scan%sfl(ith1)
       keep_it( ith) = scan%fl(ith1)
    enddo

       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)

    iring = 0
    do mindx = 0, nmvals-1
             
       mindx1 = mindx+1
       m = mvals( mindx)

       ! generate recursion factors (recfac) for Ylm of degree m
       call gen_recfac(nlmax, m, recfac)
       ! generate Ylm relation factor for degree m
       call gen_lamfac(nlmax, m, lam_fact)
       fm2 = real(m * m, kind = DP)
       normal_m = (2.0_dp * m) * (1 - m)

       ! extract needed alm under memory and CPU efficient form

       do imap = 1, nmaps
          do ll = m, nlmax
             ll1 = ll+1
             dalm(0,ll,imap) =  real(local_alms_comp_a(ll1,mindx1,imap),kind=dp) ! T, real
             dalm(1,ll,imap) = aimag(local_alms_comp_a(ll1,mindx1,imap))         ! T, imag
             dalm(2,ll,imap) =  real(local_alms_comp_b(ll1,mindx1,imap),kind=dp)*normal_l(ll) ! G, real
             dalm(3,ll,imap) = aimag(local_alms_comp_b(ll1,mindx1,imap))        *normal_l(ll)
             dalm(4,ll,imap) =  real(local_alms_comp_c(ll1,mindx1,imap),kind=dp)*normal_l(ll) ! C, real
             dalm(5,ll,imap) = aimag(local_alms_comp_c(ll1,mindx1,imap))        *normal_l(ll)
          enddo
       enddo

       jth = 0
       do ith = 0, nringsall-1

          if( keep_it( ith) == 1) then 
             l_min =  l_min_ylm(m, sth(ith)) ! lowest l of non-negligible Ylm
             if (nlmax >= l_min) then

                ! compute lam_lm(p,theta) for all l>=m

                call mpi_do_lam_lm_pol(nlmax, m, cth(ith), sth(ith), mfac(m), recfac, lam_fact, lam_lm)

                !      T =   alm_T * Ylm
                !      Q = - alm_E * Wlm - i alm_B * Xlm
                !      U = i alm_E * Xlm -   alm_B * Wlm

                b_ns = 0.0_dp

                l_start = max(m, l_min) ! even values of (l+m)
                if (mod(m+l_start,2) == 1) l_start = l_start + 1

                do imap = 1, nmaps
                   do l = l_start, nlmax, 2
                      b_ns(3: 4,imap) = b_ns(3: 4,imap) + lam_lm(1,l) * dalm(0:1,l,imap) ! T even
                      b_ns(5: 8,imap) = b_ns(5: 8,imap) - lam_lm(2,l) * dalm(2:5,l,imap) ! Q, U  even
                      b_ns(-1,imap)   = b_ns(-1,imap)   + lam_lm(3,l) * dalm(5,l,imap)   ! Q odd
                      b_ns( 0,imap)   = b_ns( 0,imap)   - lam_lm(3,l) * dalm(4,l,imap)
                      b_ns( 1,imap)   = b_ns( 1,imap)   - lam_lm(3,l) * dalm(3,l,imap)   ! U odd
                      b_ns( 2,imap)   = b_ns( 2,imap)   + lam_lm(3,l) * dalm(2,l,imap)
                   enddo
                enddo

                l_start = max(m+1, l_min) ! odd values of (l+m)
                if (mod(m+l_start,2) == 0) l_start = l_start + 1
 
                do imap = 1, nmaps
                   do l = l_start, nlmax, 2
                      b_ns(-3:-2,imap) = b_ns(-3:-2,imap) + lam_lm(1,l) * dalm(0:1,l,imap) ! T odd
                      b_ns(-1: 2,imap) = b_ns(-1: 2,imap) - lam_lm(2,l) * dalm(2:5,l,imap) ! Q, U  odd
                      b_ns( 5,imap)    = b_ns( 5,imap)    + lam_lm(3,l) * dalm(5,l,imap)   ! Q even
                      b_ns( 6,imap)    = b_ns( 6,imap)    - lam_lm(3,l) * dalm(4,l,imap)
                      b_ns( 7,imap)    = b_ns( 7,imap)    - lam_lm(3,l) * dalm(3,l,imap)   ! U even
                      b_ns( 8,imap)    = b_ns( 8,imap)    + lam_lm(3,l) * dalm(2,l,imap)
                   end do
                enddo

                iring = mindx*nringsobs+jth
				               
                b_TQU(iring,1,1:nmaps) = cmplx(b_ns(3,1:nmaps) + b_ns(-3,1:nmaps), b_ns(4,1:nmaps) + b_ns(-2,1:nmaps), kind = DP) ! T north
                b_TQU(iring,4,1:nmaps) = cmplx(b_ns(3,1:nmaps) - b_ns(-3,1:nmaps), b_ns(4,1:nmaps) - b_ns(-2,1:nmaps), kind = DP) ! T south
                b_TQU(iring,2,1:nmaps) = cmplx(b_ns(5,1:nmaps) + b_ns(-1,1:nmaps), b_ns(6,1:nmaps) + b_ns( 0,1:nmaps), kind = DP) ! Q n
                b_TQU(iring,5,1:nmaps) = cmplx(b_ns(5,1:nmaps) - b_ns(-1,1:nmaps), b_ns(6,1:nmaps) - b_ns( 0,1:nmaps), kind = DP) ! Q s
                b_TQU(iring,3,1:nmaps) = cmplx(b_ns(7,1:nmaps) + b_ns( 1,1:nmaps), b_ns(8,1:nmaps) + b_ns( 2,1:nmaps), kind = DP) ! U n
                b_TQU(iring,6,1:nmaps) = cmplx(b_ns(7,1:nmaps) - b_ns( 1,1:nmaps), b_ns(8,1:nmaps) - b_ns( 2,1:nmaps), kind = DP) ! U s

             endif ! test on nlmax

             jth = jth + 1

          endif ! test if ring is observed

       enddo ! loop over rings (ith)

    enddo ! loop over m

    deallocate( cth, sth)
    deallocate( recfac, dalm, lam_fact, lam_lm)
    deallocate( mfac, normal_l)
    deallocate( keep_it)

    allocate( kphi0(0:nrings-1), stat=status)
    call assert_alloc(status,code,'kphi0')
    kphi0 = kphi(first_ring:last_ring)
    deallocate( kphi)
    
    allocate( nph0(0:nrings-1), stat=status)
    nph0 = nph(first_ring:last_ring)
    deallocate( nph)
	
    allocate( keep_north(0:nrings-1), stat=status)
    keep_north = keep_north_all(first_ring:last_ring)
    deallocate( keep_north_all)
    allocate( keep_south(0:nrings-1), stat=status)
    keep_south = keep_south_all(first_ring:last_ring)
    deallocate( keep_south_all)

    !! redistribute

    allocate(indx(0:nmmax), stat=status)
    call assert_alloc(status,code,'indx')

    call fix_mmap( -1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

    allocate(column(0:nlen-1,1:nmaps),stat = status)
    call assert_alloc(status,code,'column')
    
    do i = 1, 6
       column = b_TQU(0:nlen-1,i,1:nmaps)
       b_TQU_ptr => column 
       call redistInterProductsAlm2Map( nrings, nmvals, nmaps, nlen, b_TQU_ptr, myid, numprocs, comm)

       do imap = 1, nmaps
          do j = 0, nmmax
             do k = 0, nrings-1
                b_TQU( (nmmax+1)*k+indx(j),i,imap) = b_TQU_ptr(k+j*nrings,imap)
             enddo
          enddo
       enddo
    enddo

    deallocate( column)

    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')

    allocate(phring(0:nmmax),stat = status)
    call assert_alloc(status,code,'phring')

    local_map = 0.d0 ! set the whole map to zero


    do imap = 1, nmaps

       do i=1,3

          istart_north = 0

          do iring = 0, nrings-1

             nphl = nph0(iring)

             if( keep_north( iring) == 1) then
		         
                !        ---------------------------------------------------------------
                !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
                !        ---------------------------------------------------------------
 
                phring(0:nmmax) = b_TQU(iring*(nmmax+1):(iring+1)*(nmmax+1)-1,i,imap)

                call s2hat_ring_synthesis(nlmax,nmmax,phring,nphl,ring,kphi0(iring))   ! north hemisph. + equator
                local_map(istart_north:istart_north+nphl-1, i, imap) = ring(0:nphl-1)
			 
	     endif

             istart_north = istart_north+nphl

	   enddo
  
	enddo                                       

	do i = 1, 3

           istart_south = map_size

           do iring = 0, nrings-1

              nphl = nph0(iring)
              istart_south = istart_south-nphl

              if( keep_south( iring) == 1) then

                  phring(0:nmmax) = b_TQU(iring*(nmmax+1):(iring+1)*(nmmax+1)-1,3+i,imap)

                  call s2hat_ring_synthesis(nlmax,nmmax,phring,nphl,ring,kphi0(iring)) ! south hemisph. w/o equat
                  local_map(istart_south:istart_south+nphl-1, i, imap) = ring(0:nphl-1)
 
              endif
           enddo   ! loop over rings

        enddo  

    enddo  ! over maps

    !     --------------------
    !     free memory and exit
    !     --------------------
    
    deallocate( ring)
    deallocate( phring)
    deallocate( b_TQU)
    deallocate( kphi0)
    deallocate( nph0)
    deallocate( indx)
    deallocate( keep_south, keep_north)

    return

  end subroutine s2hat_alm2map_pol

  subroutine s2hat_alm2map_pol_pre1( pixelization, scan, nlmax, nmmax, nmvals, mvals, imap, first_ring, last_ring, map_size, local_map, &
                                   & lda, nplm, local_plm, local_alms, numprocs, myid, comm)

    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !     for the Temperature+Polarization field using precomputed
	!     and stored scalar spherical harmonics (Ylm_T(theta))
    !=======================================================================

    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan	
    integer (i4b), intent(in) :: nlmax, nmmax, nmvals, imap, first_ring, last_ring, lda, map_size, numprocs, myid, comm
    integer(i8b), intent(in) :: nplm
    integer (i4b), dimension(0:nmvals-1), intent(in) :: mvals
    real(kind=dp), dimension(:,:), pointer :: local_plm
    complex (dp), dimension(:,:,:,:), pointer :: local_alms

    ! output
    real(dp), dimension(:,:,:), pointer :: local_map

    ! internal
    integer(I4B) :: l, m, ith, scalel, scalem ! alm related
    integer(I4B) :: istart_south, istart_north   ! map related
    integer(I4B) :: nrings, nringsall, nringsobs, npix, nphmx

    integer(I4B)  :: par_lm
    real(DP)      :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)      :: lambda_w, lambda_x, normal_m, lam_lm1m
    real(DP)      :: fm2, fl, flm1
    real(DP)      :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP)      :: ovflow, unflow, one_on_s2, c_on_s2
    real(DP),     dimension(-3:8)             :: b_ns
    real(DP),     dimension(:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: lam_fact

    integer(i4b)  :: ll, ll1, l_min, l_start, mindx, mindx1, nlen, mshift
    real(DP),     dimension(:),   allocatable :: normal_l
    complex(DPC), dimension(:,:), allocatable, target :: b_TQU
    complex(DPC), dimension(:,:), pointer :: b_TQU_ptr
    complex(DPC), dimension(:,:), allocatable, target :: column
    complex(DPC), dimension(:),   allocatable :: phring
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)  :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i, j, ij, iring, jring, kring
    integer(I8B)  :: n_lm, i_mm, i_mm1
    integer(I4B), dimension(:), allocatable:: nph0, nph, indx
    real(DP),     dimension(:), allocatable :: cth, sth, kphi, kphi0
    integer(i4b), dimension(:), allocatable :: keep_north, keep_north_all, keep_south, keep_south_all, keep_it

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status

        complex(DPC), dimension(:,:), pointer :: local_alms_comp_a, local_alms_comp_b, local_alms_comp_c
    real(DP), dimension(:), pointer :: local_plm_comp_a

    !=======================================================================

    if( lda == 3) then       ! i.e., healpix convention for the alm array

       local_alms_comp_a => local_alms(1,0:nlmax,0:nmvals-1,imap)
       local_alms_comp_b => local_alms(2,0:nlmax,0:nmvals-1,imap)
       local_alms_comp_c => local_alms(3,0:nlmax,0:nmvals-1,imap)

       local_plm_comp_a => local_plm(1,0:nplm-1)

    else                     ! i.e., s2hat convention for the alm array

       local_alms_comp_a => local_alms(0:nlmax,0:nmvals-1,1,imap)
       local_alms_comp_b => local_alms(0:nlmax,0:nmvals-1,2,imap)
       local_alms_comp_c => local_alms(0:nlmax,0:nmvals-1,3,imap)

       local_plm_comp_a => local_plm(0:nplm-1,1)

    end if

    nrings = (last_ring-first_ring+1)
	
    npix = pixelization%npixsall
    nphmx = pixelization%nphmx
    nringsall = pixelization%nringsall
    nringsobs = scan%nringsobs	

    n_lm = nummmodes( nlmax, nmvals, mvals) ! an actual size of the local_plm array	

    !     --- allocates space for arrays ---
    nchunks   = nringsall/SMAXCHK + 1  ! number of chunks
    chunksize = (nringsall+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(normal_l(0:nlmax),stat = status)
    call assert_alloc(status,code,'normal_l')

    nlen = max( (nmmax+1)*nrings, nmvals*nringsobs)
    allocate(b_TQU(0:nlen-1, 1:6),stat = status)
    call assert_alloc(status,code,'b_TQU')

    allocate(dalm(0:5,0:nlmax), lam_fact(0:nlmax),stat = status)
    call assert_alloc(status,code,'dalm & lam_fact')
    
    allocate(nph(1:nringsall), stat=status)
    call assert_alloc(status,code,'nph')

    allocate(kphi(1:nringsall), stat=status)
    call assert_alloc(status,code,'kphi')

    allocate(cth(0:SMAXCHK), stat=status)
    call assert_alloc(status,code,'cth')

    allocate(sth(0:SMAXCHK), stat=status)
    call assert_alloc(status,code,'sth')

    allocate(keep_it(0:SMAXCHK), stat=status)
    call assert_alloc(status,code,'keep_it')
    allocate(keep_north_all(1:nringsall), stat=status)
    call assert_alloc(status,code,'keep_north_all')
    allocate(keep_south_all(1:nringsall), stat=status)
    call assert_alloc(status,code,'keep_south_all')

    !     ------------ initiate variables and arrays ----------------
    
    ! generate Polarization normalisation
    call gen_normpol(nlmax, normal_l)

    b_TQU(:,:) = 0_dpc ! pad with zeros

    kring = 0
    do ichunk = 0, nchunks-1

       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nringsall)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
		  
          ! extract pixel location information
          kphi(ith) = pixelization%kphi(ith)
          nph(ith) = pixelization%nph(ith)
          cth(ithl) = pixelization%cth(ith)
          sth(ithl) = pixelization%sth(ith)

          keep_north_all( ith) = scan%nfl(ith)
          keep_south_all( ith) = scan%sfl(ith)
          keep_it( ithl) = scan%fl(ith)				  
       enddo

       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)

       mshift = 0
       do mindx = 0, nmvals-1
             
          mindx1 = mindx+1
          m = mvals( mindx)

          ! generate Ylm relation factor for degree m
          call gen_lamfac(nlmax, m, lam_fact)
          fm2 = real(m * m, kind = DP)
          normal_m = (2.0_dp * m) * (1 - m)

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             ll1 = ll+1
             dalm(0,ll) =  real(local_alms_comp_a(ll1,mindx1),kind=dp) ! T, real
             dalm(1,ll) = aimag(local_alms_comp_a(ll1,mindx1))         ! T, imag
             dalm(2,ll) =  real(local_alms_comp_b(ll1,mindx1),kind=dp)*normal_l(ll) ! G, real
             dalm(3,ll) = aimag(local_alms_comp_b(ll1,mindx1))        *normal_l(ll)
             dalm(4,ll) =  real(local_alms_comp_c(ll1,mindx1),kind=dp)*normal_l(ll) ! C, real
             dalm(5,ll) = aimag(local_alms_comp_c(ll1,mindx1))        *normal_l(ll)
          enddo

          jring = kring
          do ithl = 0, uchk-lchk

             if( keep_it( ithl) == 1) then
			 
                 l_min =  l_min_ylm(m, sth(ithl)) ! lowest l of non-negligible Ylm
                 if (nlmax >= l_min) then
                    ith = ithl + lchk
                    i_mm = n_lm * jring + mshift ! location of Ym,m for ring ith
                    i_mm1 = i_mm+1

                    !           ---------- l = m ----------
                    par_lm = 3  ! = 3 * (-1)^(l+m)
                    lam_lm = local_plm_comp_a(i_mm1)

                    one_on_s2 =    1.0_dp / sth(ithl)**2
                    c_on_s2 = cth(ithl) *one_on_s2

                    if (m >= l_min) then ! skip Ymm if too small
                       lambda_w =  (normal_m * lam_lm)  * (0.5_dp - one_on_s2)
                       lambda_x =  (normal_m * lam_lm)  * c_on_s2

                       !      T =   alm_T * Ylm
                       !      Q = - alm_E * Wlm - i alm_B * Xlm
					   !      U = i alm_E * Xlm -   alm_B * Wlm
                       b_ns(-3:-2) = 0.0_dp               ! T odd
                       b_ns(-1: 0) =   lambda_x * (/   dalm(5,m), - dalm(4,m) /) ! Q odd
                       b_ns( 1: 2) =   lambda_x * (/ - dalm(3,m),   dalm(2,m) /) ! U odd
                       b_ns( 3: 4) =   lam_lm   * dalm(0:1,m) ! T even
                       b_ns( 5: 6) = - lambda_w * dalm(2:3,m) ! Q even
                       b_ns( 7: 8) = - lambda_w * dalm(4:5,m) ! U even
                    else
                       b_ns = 0.0_dp
                    endif

                    !           ---------- l > m ----------
                    cth_ring = cth(ithl)
                    fm_on_s2     =      m * one_on_s2
                    two_on_s2    = 2.0_dp * one_on_s2
                    two_cth_ring = 2.0_dp * cth_ring
                    b_w          =  c_on_s2

                    l_start = max(m+1, l_min)
                    if (mod(m+l_start,2) == 0) l_start = l_start+1
                    do l = l_start, nlmax, 2

                       lam_lm1m = local_plm_comp_a(i_mm1+l-m-1) * lam_fact(l)
                       lam_lm   = local_plm_comp_a(i_mm1+l-m)

                       fl = real(l, kind = DP)
                       flm1 = fl - 1.0_dp
                       a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                       a_x =  two_cth_ring * flm1
                       lambda_w =                b_w * lam_lm1m - a_w * lam_lm 
                       lambda_x = fm_on_s2 * (         lam_lm1m - a_x * lam_lm)

                       b_ns(-3) = b_ns(-3) + lam_lm * dalm(0,l) 
                       b_ns(-2) = b_ns(-2) + lam_lm * dalm(1,l) 

                       b_ns(-1) = b_ns(-1) - lambda_w * dalm(2,l) 
                       b_ns(0)  = b_ns(0)  - lambda_w * dalm(3,l) 
                       b_ns(1)  = b_ns(1)  - lambda_w * dalm(4,l) 
                       b_ns(2)  = b_ns(2)  - lambda_w * dalm(5,l) 

                       b_ns(5) = b_ns(5) + lambda_x * dalm(5,l) ! Q odd (or even)
                       b_ns(6) = b_ns(6) - lambda_x * dalm(4,l)
                       b_ns(7) = b_ns(7) - lambda_x * dalm(3,l) ! U odd (or even)
                       b_ns(8) = b_ns(8) + lambda_x * dalm(2,l)

                    enddo ! loop on l

                    l_start = max(m+2, l_min)
                    if (mod(m+l_start,2) == 1) l_start = l_start+1
                    do l = l_start, nlmax, 2

                       lam_lm1m = local_plm_comp_a(i_mm1+l-m-1) * lam_fact(l)
                       lam_lm   = local_plm_comp_a(i_mm1+l-m)
                   
                       fl = real(l, kind = DP)
                       flm1 = fl - 1.0_dp
                       a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                       a_x =  two_cth_ring * flm1
                       lambda_w =                b_w * lam_lm1m - a_w * lam_lm 
                       lambda_x = fm_on_s2 * (         lam_lm1m - a_x * lam_lm)
                   
                       b_ns(3) = b_ns(3) + lam_lm * dalm(0,l) 
                       b_ns(4) = b_ns(4) + lam_lm * dalm(1,l) 
                   
                       b_ns(5) = b_ns(5) - lambda_w * dalm(2,l) 
                       b_ns(6) = b_ns(6) - lambda_w * dalm(3,l) 
                       b_ns(7) = b_ns(7) - lambda_w * dalm(4,l) 
                       b_ns(8) = b_ns(8) - lambda_w * dalm(5,l) 
                   
                       b_ns(-1) = b_ns(-1) + lambda_x * dalm(5,l) 
                       b_ns(0)  = b_ns(0)  - lambda_x * dalm(4,l)
                       b_ns(1)  = b_ns(1)  - lambda_x * dalm(3,l) 
                       b_ns(2)  = b_ns(2)  + lambda_x * dalm(2,l)
                   
                    enddo ! loop on l

                    iring = mindx*nringsobs+jring
                    b_TQU(iring, 1) = cmplx(b_ns(3) + b_ns(-3), b_ns(4) + b_ns(-2), kind = DP) ! T north
                    b_TQU(iring, 4) = cmplx(b_ns(3) - b_ns(-3), b_ns(4) - b_ns(-2), kind = DP) ! T south
                    b_TQU(iring, 2) = cmplx(b_ns(5) + b_ns(-1), b_ns(6) + b_ns( 0), kind = DP) ! Q n
                    b_TQU(iring, 5) = cmplx(b_ns(5) - b_ns(-1), b_ns(6) - b_ns( 0), kind = DP) ! Q s
                    b_TQU(iring, 3) = cmplx(b_ns(7) + b_ns( 1), b_ns(8) + b_ns( 2), kind = DP) ! U n
                    b_TQU(iring, 6) = cmplx(b_ns(7) - b_ns( 1), b_ns(8) - b_ns( 2), kind = DP) ! U s

                 endif ! test on nlmax

                 jring = jring + 1

	      endif ! test if the rings is observed

           enddo ! loop over chunk rings (ithl)

           mshift = mshift + nlmax-m+1

       enddo ! loop over m

       kring = jring

    enddo ! loop over chunks

    deallocate( cth, sth)
    deallocate( dalm, lam_fact)
    deallocate( normal_l)
    deallocate( keep_it)

    allocate( kphi0(0:nrings-1), stat=status)
    call assert_alloc(status,code,'kphi0')
    kphi0 = kphi(first_ring:last_ring)
    deallocate( kphi)
    
    allocate( nph0(0:nrings-1), stat=status)
    nph0 = nph(first_ring:last_ring)
    deallocate( nph)
	
    allocate( keep_north(0:nrings-1), stat=status)
    keep_north = keep_north_all(first_ring:last_ring)
    deallocate( keep_north_all)
    allocate( keep_south(0:nrings-1), stat=status)
    keep_south = keep_south_all(first_ring:last_ring)
    deallocate( keep_south_all)
	
    !! redistribute data here

    allocate(column(0:nlen-1,1:1),stat = status)
    call assert_alloc(status,code,'column')
    
    do i = 1, 6
       column = b_TQU(0:nlen-1,i:i)
       b_TQU_ptr => column 
       call redistInterProductsAlm2Map( nrings, nmvals, 1, nlen, b_TQU_ptr, myid, numprocs, comm)
       b_TQU(0:nlen-1,i) = b_TQU_ptr(0:nlen-1,1)
    enddo

    deallocate( column)
       
    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')

    allocate(phring(0:nmmax),stat = status)
    call assert_alloc(status,code,'phring')

    allocate(indx(0:nmmax), stat=status)
    call assert_alloc(status,code,'indx')

    call fix_mmap( 1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

    local_map(:,:,imap) = 0.d0 ! set the whole map to zero
    istart_north = 0
    istart_south = map_size

    do iring = 0, nrings-1
       
       nphl = nph0(iring)

       if( keep_north( iring) == 1) then

          !        ---------------------------------------------------------------
          !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
          !        ---------------------------------------------------------------
             
          do i=1,3
             do j = 0, nmmax
                ij = iring+indx(j)*nrings
                phring(j) = b_TQU(ij,i)
             enddo
             call s2hat_ring_synthesis(nlmax,nmmax,phring,nphl,ring,kphi0(iring))   ! north hemisph. + equator
             local_map(istart_north:istart_north+nphl-1, i, imap) = ring(0:nphl-1)
          enddo
                     
       endif

       istart_north = istart_north+nphl

       istart_south = istart_south-nphl
								                                     
       if( keep_south( iring) == 1) then

           do i=1,3
              do j = 0, nmmax
                 ij = iring+indx(j)*nrings
                 phring(j) = b_TQU(ij,3+i)
              enddo
              call s2hat_ring_synthesis(nlmax,nmmax,phring,nphl,ring,kphi0(iring)) ! south hemisph. w/o equat
              local_map(istart_south:istart_south+nphl-1, i, imap) = ring(0:nphl-1)
           enddo

        endif

    enddo    ! loop over rings

    !     --------------------
    !     free memory and exit
    !     --------------------
    
    deallocate( ring)
    deallocate( phring)
    deallocate( b_TQU)
    deallocate( kphi0)
    deallocate( nph0)
    deallocate( indx)
    deallocate( keep_north, keep_south)

    return

  end subroutine s2hat_alm2map_pol_pre1

  subroutine s2hat_alm2map_pol_pre2( pixelization, scan, nlmax, nmmax, nmvals, mvals, imap, first_ring, last_ring, map_size, local_map, &
                                   & lda, nplm, local_plm, local_alms, numprocs, myid, comm)

    !=======================================================================
    !     computes a map from its alm for the HEALPIX pixelisation
    !     for the Temperature+Polarization field using the precomputed
    !     and stored scalar and spin spherical harmonics 
    !=======================================================================

    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan	
    integer (i4b), intent(in) :: nlmax, nmmax, nmvals, imap, first_ring, last_ring, lda, map_size, numprocs, myid, comm
    integer(i8b), intent(in) :: nplm
    integer (i4b), dimension(0:nmvals-1), intent(in) :: mvals
    real(kind=dp), dimension(:,:), pointer :: local_plm
    complex (dp), dimension(:,:,:,:), pointer :: local_alms

    ! output
    real(dp), dimension(:,:,:), pointer :: local_map

    ! internal
    integer(I4B) :: l, m, ith, scalel, scalem    ! alm related
    integer(I4B) :: istart_south, istart_north   ! map related
    integer(I4B) :: nrings, nringsall, nringsobs, npix, nphmx

    integer(I4B) :: par_lm
    real(DP),  dimension(-3:8)  :: b_ns
    real(DP),  dimension(:,:), allocatable ::  dalm, plm_sub
    real(DP) :: cth0

    integer(i4b)  :: ll, ll1, l_min, mindx, mindx1, nlen, mshift
    complex(DPC), dimension(:,:), allocatable, target :: b_TQU
    complex(DPC), dimension(:,:), pointer :: b_TQU_ptr
    complex(DPC), dimension(:,:), allocatable, target :: column
    complex(DPC), dimension(:),   allocatable :: phring
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)  :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i, j, ij, iring, jring, kring
    integer(I8B)  :: n_lm, i_mm, i_mm1, i_up, i_up1
    integer(I4B), dimension(:), allocatable:: nph0, nph, indx
    integer(i4b), dimension(:), allocatable :: keep_north, keep_north_all, keep_south, keep_south_all, keep_it

    real(kind = dp), dimension(:), allocatable :: sth, kphi, kphi0

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status, ierr

    complex(DPC), dimension(:,:), pointer :: local_alms_comp_a, local_alms_comp_b, local_alms_comp_c
    real(DP), dimension(:), pointer :: local_plm_comp_a, local_plm_comp_b, local_plm_comp_c

    !=======================================================================

    if( lda == 3) then       ! i.e., healpix convention for the alm array

       local_alms_comp_a => local_alms(1,0:nlmax,0:nmvals-1, imap)
       local_alms_comp_b => local_alms(2,0:nlmax,0:nmvals-1, imap)
       local_alms_comp_c => local_alms(3,0:nlmax,0:nmvals-1, imap)

       local_plm_comp_a => local_plm(1,0:nplm-1)
       local_plm_comp_b => local_plm(2,0:nplm-1)
       local_plm_comp_b => local_plm(3,0:nplm-1)

    else                     ! i.e., s2hat convention for the alm array

       local_alms_comp_a => local_alms(0:nlmax,0:nmvals-1,1,imap)
       local_alms_comp_b => local_alms(0:nlmax,0:nmvals-1,2,imap)
       local_alms_comp_c => local_alms(0:nlmax,0:nmvals-1,3,imap)

       local_plm_comp_a => local_plm(0:nplm-1,1)
       local_plm_comp_b => local_plm(0:nplm-1,2)
       local_plm_comp_c => local_plm(0:nplm-1,3)

    end if

    nrings = (last_ring-first_ring+1)

    npix = pixelization%npixsall
    nphmx = pixelization%nphmx
    nringsall = pixelization%nringsall
    nringsobs = scan%nringsobs

    n_lm = nummmodes( nlmax, nmvals, mvals) ! an actual size of the local_plm array	

    !     --- allocates space for arrays ---
	
    nchunks   = nringsall/SMAXCHK + 1         ! number of chunks
    chunksize = (nringsall+nchunks-1)/nchunks ! <= SMAXCHK

    nlen = max( (nmmax+1)*nrings, nmvals*nringsobs)
    allocate(b_TQU(0:nlen-1, 1:6),stat = status)
    call assert_alloc(status,code,'b_TQU')

    allocate(dalm(0:5,0:nlmax),stat = status)
    call assert_alloc(status,code,'dalm')

    allocate(plm_sub(1:3,0:nlmax),stat = status)
    call assert_alloc(status,code,'plm_sub')
    
    allocate(nph(1:nringsall), stat=status)
    call assert_alloc(status,code,'nph')

    allocate(kphi(1:nringsall), stat=status)
    call assert_alloc(status,code,'kphi')

    allocate(sth(0:SMAXCHK), stat=status)
    call assert_alloc(status,code,'sth')

    allocate(keep_it(0:SMAXCHK), stat=status)
    call assert_alloc(status,code,'keep_it')
    allocate(keep_north_all(1:nringsall), stat=status)
    call assert_alloc(status,code,'keep_north_all')
    allocate(keep_south_all(1:nringsall), stat=status)
    call assert_alloc(status,code,'keep_south_all')

    !     ------------ initiate variables and arrays ----------------
    
    b_TQU(:,:) = 0_dpc                    ! pad with zeros

    kring = 0
    do ichunk = 0, nchunks-1

       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nringsall)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
		  
          ! extract pixel location information
          kphi(ith) = pixelization%kphi(ith)
          nph(ith) = pixelization%nph(ith)
          sth(ithl) = pixelization%sth(ith)

          keep_north_all( ith) = scan%nfl(ith)
          keep_south_all( ith) = scan%sfl(ith)
          keep_it( ithl) = scan%fl(ith)
		  		  
       enddo

       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)

       mshift = 0

       do mindx = 0, nmvals-1
             
          mindx1 = mindx+1
          m = mvals( mindx)

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             ll1 = ll+1
             dalm(0,ll) =  real(local_alms_comp_a(ll1,mindx1),kind=dp) ! T, real
             dalm(1,ll) = aimag(local_alms_comp_a(ll1,mindx1))         ! T, imag
             dalm(2,ll) =  real(local_alms_comp_b(ll1,mindx1),kind=dp) ! G, real
             dalm(3,ll) = aimag(local_alms_comp_b(ll1,mindx1))
             dalm(4,ll) =  real(local_alms_comp_c(ll1,mindx1),kind=dp) ! C, real
             dalm(5,ll) = aimag(local_alms_comp_c(ll1,mindx1))
          enddo

          jring = kring		  
          do ithl = 0, uchk - lchk
		  
             if( keep_it( ithl) == 1) then
			 
                l_min =  l_min_ylm(m, sth(ithl)) ! lowest l of non-negligible Ylm
                if (nlmax >= l_min) then
                   ith = ithl + lchk
                   i_mm = n_lm * jring + mshift ! location of Ym,m for ring ith
                   i_up = i_mm + nlmax - m

                   i_mm1 = i_mm+1
                   i_up1 = i_up+1

                   plm_sub(1,m:nlmax) = local_plm_comp_a(i_mm1:i_up1)
                   plm_sub(2,m:nlmax) = local_plm_comp_b(i_mm1:i_up1)
                   plm_sub(3,m:nlmax) = local_plm_comp_c(i_mm1:i_up1)

                   !           ---------- l = m ----------

                   par_lm = 3  ! = 3 * (-1)^(l+m)

                   if (m >= l_min) then ! skip Ymm if too small
                      !      T =   alm_T * Ylm
                      !      Q = - alm_E * Wlm - i alm_B * Xlm
                      !      U = i alm_E * Xlm -   alm_B * Wlm
                      b_ns(-3:-2) = 0.0_dp               ! T odd
                      b_ns(-1: 0) =   plm_sub(3,m) * (/   dalm(5,m), - dalm(4,m) /) ! Q odd
                      b_ns( 1: 2) =   plm_sub(3,m) * (/ - dalm(3,m),   dalm(2,m) /) ! U odd
                      b_ns( 3: 4) =   plm_sub(1,m)   * dalm(0:1,m) ! T even
                      b_ns( 5: 6) = - plm_sub(2,m) * dalm(2:3,m) ! Q even
                      b_ns( 7: 8) = - plm_sub(2,m) * dalm(4:5,m) ! U even
                   else
                      b_ns = 0.0_dp
                   endif

                   !           ---------- l > m ----------

                   do l = m+1, nlmax
                      par_lm = - par_lm  ! = 3 * (-1)^(l+m)
                      if (l >= l_min) then

                         b_ns(par_lm:  par_lm+1) = b_ns(par_lm:  par_lm+1) + plm_sub(1,l) * dalm(0:1,l) ! T even or odd
                         b_ns(par_lm+2:par_lm+5) = b_ns(par_lm+2:par_lm+5) - plm_sub(2,l) * dalm(2:5,l) ! Q, U  even or odd
                         b_ns(2-par_lm) = b_ns(2-par_lm) + plm_sub(3,l) * dalm(5,l) ! Q odd (or even)
                         b_ns(3-par_lm) = b_ns(3-par_lm) - plm_sub(3,l) * dalm(4,l)
                         b_ns(4-par_lm) = b_ns(4-par_lm) - plm_sub(3,l) * dalm(3,l) ! U odd (or even)
                         b_ns(5-par_lm) = b_ns(5-par_lm) + plm_sub(3,l) * dalm(2,l)
                      endif

                   enddo ! loop on l

                   iring = mindx*nringsobs+jring
                   b_TQU(iring,1) = cmplx(b_ns(3) + b_ns(-3), b_ns(4) + b_ns(-2), kind = DP) ! T north
                   b_TQU(iring,4) = cmplx(b_ns(3) - b_ns(-3), b_ns(4) - b_ns(-2), kind = DP) ! T south
                   b_TQU(iring,2) = cmplx(b_ns(5) + b_ns(-1), b_ns(6) + b_ns( 0), kind = DP) ! Q n
                   b_TQU(iring,5) = cmplx(b_ns(5) - b_ns(-1), b_ns(6) - b_ns( 0), kind = DP) ! Q s
                   b_TQU(iring,3) = cmplx(b_ns(7) + b_ns( 1), b_ns(8) + b_ns( 2), kind = DP) ! U n
                   b_TQU(iring,6) = cmplx(b_ns(7) - b_ns( 1), b_ns(8) - b_ns( 2), kind = DP) ! U s

                endif ! test on nlmax
      
                jring = jring + 1

             endif ! test if the ring is observed

          enddo ! loop over chunk rings (ithl)

	  mshift = mshift + nlmax-m+1

       enddo ! loop over m
 
      kring = jring

    enddo ! loop over chunks

    deallocate( dalm)
    deallocate( plm_sub)
    deallocate( keep_it)
    deallocate( sth)

    allocate( kphi0(0:nrings-1), stat=status)
    call assert_alloc(status,code,'kphi0')
    kphi0 = kphi(first_ring:last_ring)
    deallocate( kphi)

    allocate( nph0(0:nrings-1), stat=status)
    nph0 = nph(first_ring:last_ring)
    deallocate( nph)

    allocate( keep_north(0:nrings-1), stat=status)
    keep_north = keep_north_all(first_ring:last_ring)
    deallocate( keep_north_all)
    allocate( keep_south(0:nrings-1), stat=status)
    keep_south = keep_south_all(first_ring:last_ring)
    deallocate( keep_south_all)

    !! redistribute data here

    allocate(column(0:nlen-1,1:1),stat = status)
    call assert_alloc(status,code,'column')
    
    do i = 1, 6
       column = b_TQU(0:nlen-1,i:i)
       b_TQU_ptr => column 
       call redistInterProductsAlm2Map( nrings, nmvals, 1, nlen, b_TQU_ptr, myid, numprocs, comm)
       b_TQU(0:nlen-1,i) = b_TQU_ptr(0:nlen-1,1)
    enddo
 
    deallocate( column)
      
    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')

    allocate(phring(0:nmmax),stat = status)
    call assert_alloc(status,code,'phring')

    allocate(indx(0:nmmax), stat=status)
    call assert_alloc(status,code,'indx')

    call fix_mmap( 1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

    local_map(:,:,imap) = 0.d0 ! set the whole map to zero

    do i = 1,3

       istart_north = 0

       do iring = 0, nrings-1
       
          nphl = nph0(iring)

          if( keep_north( iring) == 1) then

             !        ---------------------------------------------------------------
             !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
             !        ---------------------------------------------------------------

             do j = 0, nmmax
                ij = iring+indx(j)*nrings
                phring(j) = b_TQU(ij,i)
             enddo
             call s2hat_ring_synthesis(nlmax,nmmax,phring,nphl,ring,kphi0(iring))   ! north hemisph. + equator
             local_map(istart_north:istart_north+nphl-1, i, imap) = ring(0:nphl-1)
          endif

          istart_north = istart_north+nphl

       enddo

    enddo

    do i = 1,3

       istart_south = map_size

       do iring = 0, nrings-1
       
          nphl = nph0(iring)
          istart_south = istart_south-nphl

          if( keep_south( iring) == 1) then

              do j = 0, nmmax
                  ij = iring+indx(j)*nrings
                 phring(j) = b_TQU(ij,3+i)
              enddo
              call s2hat_ring_synthesis(nlmax,nmmax,phring,nphl,ring,kphi0(iring)) ! south hemisph. w/o equat
              local_map(istart_south:istart_south+nphl-1, i, imap) = ring(0:nphl-1)

          endif

       enddo    ! loop over rings

    enddo ! over stokes

    !     --------------------
    !     free memory and exit
    !     --------------------
    
    deallocate( ring)
    deallocate( phring)
    deallocate( b_TQU)
    deallocate( kphi0)
    deallocate( nph0)
    deallocate( indx)
    deallocate( keep_north, keep_south)

    return

  end subroutine s2hat_alm2map_pol_pre2

end module s2hat_alm2map_mod
