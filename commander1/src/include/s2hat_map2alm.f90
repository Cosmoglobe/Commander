
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
! objects over all processors. To facilitate that the spherical harmonic 
! transforms are split into two disjoint steps which are separated by 
! a data redistribution (a "transpose operation"), involving a global 
! (all2all) mpi communication over all the procs as it is performed in the 
! routines redistInterProducts*.
!
! This file contains map-to-alm transforms. The complementary file
! s2hat_alm2map.f90 contains the alm-to-map transforms.
! Routines performing the load balanced data distribution and some other
! auxiliary procedures are included in the file: s2hat_toolbox.f90.
!
!--------------------------------------------------------------------------
! February, 2007 RS@APC
!
! - load balancing improved for partial sky experiments
! - generalized to arbitrary isolatitudal pixelization scheme with same area pixels 
!   equidistant in azimuth for each isolatitude ring, and symmetric wrt to
!   the equator.
!
!--------------------------------------------------------------------------
!
! Jul/Aug, 2007 RS@APC
!
! spherical harmonic transforms for arbitrary spin values added.
!   
!--------------------------------------------------------------------------
!
! Oct, 2007 - rs@apc
!
! quadrature wghts and arbitrary pixels ring offsest added (i.e., adapted
! to the new pixelization struct format)
!
!--------------------------------------------------------------------------

module s2hat_map2alm_mod

  use healpix_types
  use s2hat_types
  ! use s2hat_pixelization
  use healpix_fft
  use misc_utils
  use alm_tools
  use s2hat_toolbox_mod
  use s2hat_alm2map_mod

  implicit none

  !--------------------------------------------------------------------------

  public :: s2hat_map2alm, s2hat_map2alm_spin

contains

  ! ***************************************************************
  !            Wrapper for the computational routines
  ! ***************************************************************
  
  subroutine s2hat_map2alm( precompute_plms, pixelization, scan, nlmax, nmmax, nmvals, mvals, nmaps, nstokes, first_ring, last_ring, local_w8ring, &
                          & map_size, local_map, lda, local_alm, nplm, local_plm, numprocs, myid, comm)

    implicit none
  
    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan
    integer(i4b), intent( in) :: nlmax, nmmax, nmvals, nstokes, nmaps, first_ring, last_ring, lda, map_size, numprocs, myid, comm, precompute_plms
    integer(i8b), intent(in) :: nplm
    integer(i4b), dimension(0:nmvals-1), intent(in) :: mvals
    real(kind=dp), dimension(:,:,:), pointer :: local_map
    real(kind=dp), dimension(:,:), pointer :: local_plm
    real(kind=dp), dimension(:,:), pointer :: local_w8ring

    ! output
    complex(kind=dp), dimension(:,:,:,:), pointer :: local_alm

    ! internal
    integer(i4b) :: imap
  
    if (nstokes == 3) then
       ! Polarization routines
       if (precompute_plms == 1) then
          do imap = 1, nmaps
             call s2hat_map2alm_pol_pre1( pixelization, scan, nlmax, nmmax, nmvals, mvals, imap, first_ring, last_ring, local_w8ring, map_size, local_map, lda, nplm, local_plm, local_alm, numprocs, myid, comm)
          enddo
       else if (precompute_plms == 2) then
          do imap = 1, nmaps
             call s2hat_map2alm_pol_pre2( pixelization, scan, nlmax, nmmax, nmvals, mvals, imap, first_ring, last_ring, local_w8ring, map_size, local_map, lda, nplm, local_plm, local_alm, numprocs, myid, comm)
          enddo
       else
          call s2hat_map2alm_pol( pixelization, scan, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, local_w8ring, map_size, local_map, lda, local_alm, numprocs, myid, comm)
       end if
    else
       ! Temperature routines
       if (precompute_plms == 1 .or. precompute_plms == 2) then
          do imap = 1, nmaps
             call s2hat_map2alm_sc_pre( pixelization, scan, nlmax, nmmax, nmvals, mvals, imap, first_ring, last_ring, local_w8ring, map_size, local_map, lda, nplm, local_plm, local_alm, numprocs, myid, comm)
          enddo 
       else
          call s2hat_map2alm_sc( pixelization, scan, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, local_w8ring, map_size, local_map, lda, local_alm, numprocs, myid, comm)
       end if
    end if

  end subroutine s2hat_map2alm

  ! rearranges and redistributes over procs the data vector "data". On the input on each proc this
  ! vecotr stores results of the \phi transform for each pixel ring assigned to the proc and for 
  ! **all** m values. The latter are aranged so that first go m values stored on proc 0, then 1 etc.
  !
  ! On the output the data vector contains the info for every m mode to be processed by the procs
  ! and **all** pixel rings.

  subroutine redistInterProductsMap2Alm( nrings, nmvals, nmaps, ndata, data, myid, numprocs, comm)

   implicit none

   ! input
   integer(i4b), intent( in) :: nrings, nmvals, myid, nmaps, ndata, numprocs, comm
   ! input/output
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
      send_counts(i) = nrings*mperproc(i)
      recv_counts(i) = nmvals*ringsperproc(i)
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

      call mpi_alltoallv( data(0,imap), send_counts, send_disps, mpi_double_complex, tmp, recv_counts, recv_disps, mpi_double_complex, comm, ierr)

      ! reshuffle to a ring-wise ordering here

      ii = 0
      do i = 0, nmvals-1
         jj = 0
         do j = 0, numprocs-1
            jj = jj+i*ringsperproc(j)
            do k = 0, ringsperproc(j)-1
               data( ii, imap) = tmp( jj+k)
               ii = ii+1
            end do
            jj = jj+(nmvals-i)*ringsperproc(j)
          end do
       end do
     
      deallocate( tmp)

   enddo

   deallocate( send_counts, send_disps, recv_counts, recv_disps)
   deallocate( ringsperproc)

  end subroutine redistInterProductsMap2Alm
  
  ! ***************************************************************
  !                Computational routines
  ! ***************************************************************

  subroutine s2hat_map2alm_sc( pixelization, scan, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, local_w8ring, map_size, local_map, lda, local_alms, numprocs, myid, comm)

    !=================================================================================
    !     computes the a(l,m) from a Temperature map for the HEALPIX pixelisation
    !     calculating spherical harmonics "on the fly"
    !=================================================================================

    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan

    integer (i4b) :: nlmax, nmmax, nmvals, nmaps, first_ring, last_ring, lda, map_size, numprocs, myid, comm
    integer (i4b), dimension(0:nmvals-1) :: mvals
    real(dp), dimension(:,:,:), pointer :: local_map
    real(dp), dimension(:,:), pointer :: local_w8ring

    ! output
    complex (dp), dimension(:,:,:,:), pointer :: local_alms

    ! internal
    integer(I4B) :: l, l1, m, ith, ithl, scalem, scalel   ! alm related
    integer(I4B) :: istart_south, istart_north, istart_south_save, istart_north_save  ! map related
    integer(I4B) :: nrings, nphmx

    integer(I4B) :: par_lm
    real(DP)     :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)     :: ovflow, unflow, cth0, sth0
    real(DP),     dimension(-1:2,1:nmaps)     :: phas_sd
    real(DP),     dimension(:,:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: mfac
    real(DP),     dimension(:,:), allocatable :: recfac
    integer(I4B)  :: l_min, kphi1, nphas, mindx, mindx1, chunksize, nchunks, lchk, uchk, ichunk
    complex(DPC), dimension(1:nmaps)  :: php, phm
    complex(DPC), dimension(:,:), allocatable, target :: phas_n, phas_s
    complex(DPC), dimension(:,:), pointer :: phas_n_ptr, phas_s_ptr
    complex(DPC), dimension(:), allocatable :: phring
    real(DP),     dimension(:), allocatable :: ring
    integer(I4B)  :: i, imap, iring, jring, kring, nphl, ij, nringsall, nringsobs
    integer(I8B)  :: npix
    real(DP),     dimension(:), allocatable :: cth, sth, kphi0, totWght
    integer(i4b), dimension(:), allocatable :: keep_it, keep_it_all, keep_north, keep_south
    integer(i4b), dimension(:), allocatable :: indx, nph

    character(LEN=*), PARAMETER :: code = 'MAP2ALM'
    integer(I4B) :: status, ierr

    complex(DPC), dimension(:,:,:), pointer :: local_alms_comp_a

    !=======================================================================

    if( lda == 1) then       ! i.e., healpix convention for the alm array

       local_alms_comp_a => local_alms(1,0:nlmax,0:nmvals-1,1:nmaps)

    else                     ! i.e., s2hat convention for the alm array

       local_alms_comp_a => local_alms(0:nlmax,0:nmvals-1,1,1:nmaps)

    end if

    ! Healpix definitions

    nrings = (last_ring-first_ring+1)
    npix   = pixelization%npixsall
    nphmx  = pixelization%nphmx
    nringsall = pixelization%nringsall
    nringsobs = scan%nringsobs

    !     --- allocates space for arrays ---

    nphas = max( (nmmax+1)*nrings, nmvals*nringsobs)

    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(phas_n(0:nphas-1,1:nmaps),stat = status)  ! extra memory as needed for the theta inegral
    call assert_alloc(status,code,'phas_n')

    allocate(phas_s(0:nphas-1,1:nmaps),stat = status)   ! extra memory as needed for the theta inegral
    call assert_alloc(status,code,'phas_s')

    allocate( keep_it(0:nrings-1), stat = status)
    call assert_alloc(status,code,'keep_it')

    allocate( kphi0(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'kphi0')

    allocate( totWght(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'totWght')

    allocate( keep_north(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'keep_north')

    allocate( keep_south(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'keep_south')

    allocate( nph(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'nph')

    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')

    allocate(phring(0:nmmax),stat = status)
    call assert_alloc(status,code,'phring')

    allocate(indx(0:nmmax), stat=status)
    call assert_alloc(status,code,'indx')

    call fix_mmap( 1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

    !     ------------ initiate variables and arrays ----------------

    istart_north = 0
    istart_south = map_size

    phas_n(:,:) = 0_dpc
    phas_s(:,:) = 0_dpc

    ! loop on chunks

    do ichunk = 0, nchunks-1

       lchk  = first_ring + ichunk * chunksize 
       uchk  = min(lchk+chunksize - 1, last_ring)

       do ith = lchk, uchk
          ithl = ith - lchk !local index

          ! extract pixel location information
          nph(ithl) = pixelization%nph(ith)
          cth0 = pixelization%cth(ith)
          sth0 = pixelization%sth(ith)
          kphi0(ithl) = pixelization%kphi(ith)
          totWght(ithl) = pixelization%parea(ith)*pixelization%qwght(ith)

          ! find out which rings are to be analysed
          keep_north( ithl) = scan%nfl( ith)
          keep_south( ithl) = scan%sfl( ith)
          keep_it( ith-first_ring) = scan%fl( ith)
       enddo

       !-----------------------------------------------------------------------
       !  computes the integral in phi : phas_m(theta)
       !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
       !-----------------------------------------------------------------------

       istart_north_save = istart_north
       istart_south_save = istart_south

       do imap = 1, nmaps

          istart_north = istart_north_save
          istart_south = istart_south_save

          do ith = lchk, uchk

             ithl = ith - lchk       ! local ring index in the chunk
             iring = ith-first_ring  ! in the procs rings

             nphl = nph(ithl)

             ! do Fourier Transform on rings

             if (keep_north(ithl) == 1) then
                 ring(0:nphl-1) = local_map(istart_north:istart_north+nphl-1,1,imap) * &
                                & local_w8ring(ith-first_ring+1,1)*totWght( ithl)
                 call s2hat_ring_analysis(nlmax, nmmax, ring, nphl, phring, kphi0(ithl))
                 istart_north = istart_north+nphl
              
                 ! store in a transposed order for the future alltoall
                 do i = 0, nmmax
                    ij = iring+indx(i)*nrings
                    phas_n(ij,imap) = phring(i)
                 enddo 
              endif

              if (keep_south(ithl) == 1) then
                 istart_south = istart_south-nphl
                 ring(0:nphl-1) = local_map(istart_south:istart_south+nphl-1,1,imap) * &
                                & local_w8ring(ith-first_ring+1,1)*totWght(ithl)
                 call s2hat_ring_analysis(nlmax, nmmax, ring, nphl, phring, kphi0(ithl))

                 ! store in a transposed order for the future alltoall
                 do i = 0, nmmax
                    ij = iring+indx(i)*nrings
                    phas_s(ij,imap) = phring(i)
                 enddo 
              endif

          enddo ! loop over local ring

      enddo  

   enddo ! loop over chunks

   deallocate( indx)
   deallocate ( nph)
   deallocate( ring)
   deallocate( phring)
   deallocate( kphi0, totWght)
   deallocate( keep_north, keep_south)

   ! redistribute intermediate products over procs

   phas_n_ptr => phas_n
   phas_s_ptr => phas_s

   call redistInterProductsMap2Alm( nrings, nmvals, nmaps, nphas, phas_n_ptr, myid, numprocs, comm)
   call redistInterProductsMap2Alm( nrings, nmvals, nmaps, nphas, phas_s_ptr, myid, numprocs, comm)

   ! and get global info about all the rings

   allocate( keep_it_all(0:nringsall-1))
   keep_it_all(:) = 0
   call gatherRingInfo( first_ring, nrings, keep_it, nringsall, keep_it_all, myid, numprocs, comm)
   deallocate( keep_it)

   !-----------------------------------------------------------------------
   !              computes the a_lm by integrating over theta
   !                  lambda_lm(theta) * phas_m(theta)
   !                         for each m and l
   !-----------------------------------------------------------------------

   allocate(mfac(0:nmmax),stat = status)
   call assert_alloc(status,code,'mfac')

   call gen_mfac(nmmax,mfac)

   call init_rescale_local()   ! innitialize rescale table

   ovflow = rescale_tab(1)
   unflow = rescale_tab(-1)

   allocate(recfac(0:1,0:nlmax),stat = status)
   call assert_alloc(status,code,'recfac')

   allocate(dalm(1:2,1:nmaps,0:nlmax),stat = status)
   call assert_alloc(status,code,'dalm')

   nchunks   = nringsall/SMAXCHK + 1  ! number of chunks
   chunksize = (nringsall+nchunks-1)/nchunks ! <= SMAXCHK

   local_alms = 0.0 ! set the whole alm array to zero


   kring = 0
   do ichunk = 0, nchunks-1

      lchk  = ichunk * chunksize + 1
      uchk  = min(lchk+chunksize - 1, nringsall)

      do mindx = 0, nmvals-1

         mindx1 = mindx+1

         m = mvals( mindx)

         ! generate recursion factors (recfac) for Ylm of degree m

         call gen_recfac(nlmax, m, recfac)

         ! introduce double precision vector to perform summation over ith for each l

         dalm(1:2, 1:nmaps, m:nlmax) = 0.0_dp

         jring = kring
         do ithl = 0, uchk-lchk

            iring = ithl+lchk-1

            l_min = l_min_ylm(m, pixelization%sth(iring+1))

            if (keep_it_all(iring) == 1) then

                if (nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001) split in two (RS, 02-2007)

                   ! determine lam_mm
                   call mpi_compute_lam_mm(mfac(m), m, pixelization%sth(iring+1), lam_mm, corfac, scalem)

                   !           ---------- l = m ----------

                   par_lm = 1
                   php(1:nmaps) = phas_n(mindx*nringsobs+jring,1:nmaps) + phas_s(mindx*nringsobs+jring,1:nmaps) ! sum  (if (l+m) even)
                   phm(1:nmaps) = phas_n(mindx*nringsobs+jring,1:nmaps) - phas_s(mindx*nringsobs+jring,1:nmaps) ! diff (if (l+m) odd)
                   do imap = 1, nmaps
                      phas_sd(-1:0,imap) =  (/ real(phm(imap), kind=dp), aimag(phm(imap)) /)
                      phas_sd( 1:2,imap) =  (/ real(php(imap), kind=dp), aimag(php(imap)) /)
                   enddo

                   if (m >= l_min) then
                      lam_lm = corfac * lam_mm !Actual lam_mm 
                      dalm(1:2, 1:nmaps, m) = dalm(1:2, 1:nmaps, m) + lam_lm * phas_sd(par_lm:par_lm+1,1:nmaps)
                   endif

                   !           ---------- l > m ----------

                   lam_0 = 0.0_dp
                   lam_1 = 1.0_dp
                   scalel=0
                   cth_ring = pixelization%cth(iring+1)
                   lam_2 = cth_ring * lam_1 * recfac(0,m)

                   do l = m+1, nlmax
                      par_lm = - par_lm  ! = (-1)^(l+m)
                      if (l >= l_min) then
                         lam_lm = lam_2 * corfac * lam_mm
                         dalm(1:2, 1:nmaps, l) = dalm(1:2, 1:nmaps, l) &
                              &       + lam_lm * phas_sd(par_lm:par_lm+1, 1:nmaps)
                      endif

                      lam_0 = lam_1 * recfac(1,l-1)
                      lam_1 = lam_2
                      lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)

                      if (abs(lam_2) > OVFLOW) then
                         lam_1 = lam_1*UNFLOW
                         lam_2 = lam_2*UNFLOW
                         scalel= scalel + 1
                         corfac= rescale_tab(max(scalem+scalel,RSMIN))
                      elseif (abs(lam_2) < UNFLOW) then
                         lam_1 = lam_1*OVFLOW
                         lam_2 = lam_2*OVFLOW
                         scalel= scalel - 1
                         corfac= rescale_tab(max(scalem+scalel,RSMIN))
                      endif

                   enddo ! loop over l

                endif ! test for l values

                jring = jring + 1
    
             endif ! test on cut sky and nlmax

          enddo ! loop over chunks rings (ithl)

          do imap = 1, nmaps
             do l = m, nlmax
                l1 = l+1
                local_alms_comp_a(l1,mindx1,imap) = local_alms_comp_a(l1,mindx1,imap) + &
                                                  & cmplx(dalm(1, imap, l), dalm(2, imap, l), kind=DP)
             enddo
          enddo

       enddo ! loop over m

       kring = jring ! store the number of obs rings so far       

   enddo ! loop over chunks

   !     --------------------
   !     free memory and exit
   !     --------------------

   deallocate( recfac,dalm)
   deallocate( keep_it_all)
   deallocate( phas_n,phas_s)
   deallocate( mfac)

   return

  end subroutine s2hat_map2alm_sc


  subroutine s2hat_map2alm_sc_pre( pixelization, scan, nlmax, nmmax, nmvals, mvals, imap, first_ring, last_ring, local_w8ring, map_size, local_map, &
                                 & lda, nplm, local_plm, local_alms, numprocs, myid, comm)

    !=================================================================================
    !     computes the a(l,m) from a Temperature map for the HEALPIX pixelisation
    !     using precomputed and stored Ylm
    !=================================================================================

    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan
    integer (i4b) ::  nlmax, nmmax, nmvals, imap, first_ring, last_ring, lda, map_size, numprocs, myid, comm
    integer (i4b), dimension(0:nmvals-1) :: mvals
    real(dp), dimension(:,:,:), pointer :: local_map
    real(dp), dimension(:,:), pointer :: local_w8ring
    integer(i8b), intent(in) :: nplm
    real(dp), dimension(:,:), pointer :: local_plm	

    ! output
    complex (dp), dimension(:,:,:,:), pointer :: local_alms

    ! internal
    integer(I4B) :: l, l1, m, ith, ithl, scalem, scalel   ! alm related
    integer(I4B) :: istart_south, istart_north  ! map related
    integer(I4B) :: nrings, nphmx

    integer(I4B)  :: par_lm, l_start
    real(DP)      :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)      :: cth0, sth0
    real(DP),     dimension(-1:2)             :: phas_sd
    real(DP),     dimension(:,:), allocatable :: dalm
    integer(I4B)  :: l_min, kphi1, nphas, mindx, mindx1, chunksize, nchunks, lchk, uchk, ichunk
    complex(DPC)  :: php, phm
    complex(DPC), dimension(:,:), pointer :: phas_n_ptr, phas_s_ptr
    complex(DPC), dimension(:,:), allocatable, target :: phas_n, phas_s
    complex(DPC), dimension(:), allocatable :: phring
    real(DP),     dimension(:), allocatable :: ring
    integer(I4B)                            :: i, iring, jring, kring, nphl, ij, nringsall, nringsobs
    integer(I8B)                            :: npix, n_lm, mshift, i_mm, i_mm1
    real(DP),     dimension(:), allocatable :: cth, sth, kphi0, totWght
    integer(i4b), dimension(:), allocatable :: keep_it, keep_it_all, keep_north, keep_south
    integer(i4b), dimension(:), allocatable :: indx, nph

    character(LEN=*), PARAMETER :: code = 'MAP2ALM'
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

    ! Healpix definitions

    nrings = (last_ring-first_ring+1)
    npix   = pixelization%npixsall
    nphmx  = pixelization%nphmx

    n_lm = nummmodes( nlmax, nmvals, mvals)

    !     --- allocates space for arrays ---

    nringsall = pixelization%nringsall
    nringsobs = scan%nringsobs

    nphas = max( (nmmax+1)*nrings, nmvals*nringsobs)

    nchunks   = nrings/SMAXCHK + 1         ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(phas_n(0:nphas-1,1:1),stat = status)  ! extra memory as needed for the theta inegral
    call assert_alloc(status,code,'phas_n')

    allocate(phas_s(0:nphas-1,1:1),stat = status)   ! extra memory as needed for the theta inegral
    call assert_alloc(status,code,'phas_s')

    allocate( keep_it(0:nrings-1), stat = status)
    call assert_alloc(status,code,'keep_it')

    allocate( kphi0(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'kphi0')

    allocate( totWght(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'totWght')

    allocate( keep_north(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'keep_north')

    allocate( keep_south(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'keep_south')

    allocate( nph(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'nph')

    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')

    allocate(phring(0:nmmax),stat = status)
    call assert_alloc(status,code,'phring')

    allocate(indx(0:nmmax), stat=status)
    call assert_alloc(status,code,'indx')

    call fix_mmap( 1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

    !     ------------ initiate variables and arrays ----------------

    istart_north = 0
    istart_south = map_size

    phas_n(:,:) = 0_dpc
    phas_s(:,:) = 0_dpc

    ! loop on chunks

    do ichunk = 0, nchunks-1

       lchk  = first_ring + ichunk * chunksize 
       uchk  = min(lchk+chunksize - 1, last_ring)

       do ith = lchk, uchk
          ithl = ith - lchk !local index

          ! extract pixel location information
          nph(ithl) = pixelization%nph(ith)
          cth0 = pixelization%cth(ith)
          sth0 = pixelization%sth(ith)
          kphi0(ithl) = pixelization%kphi(ith)
          totWght(ithl) = pixelization%parea(ith)*pixelization%qwght(ith)

          ! find out which rings are to be analysed
          keep_north( ithl) = scan%nfl( ith)
          keep_south( ithl) = scan%sfl( ith)
          keep_it( ith-first_ring) = scan%fl( ith)
       enddo

       !-----------------------------------------------------------------------
       !  computes the integral in phi : phas_m(theta)
       !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
       !-----------------------------------------------------------------------

       do ith = lchk, uchk

          ithl = ith - lchk       ! local ring index in the chunk
          iring = ith-first_ring  ! in the procs rings

          nphl = nph(ithl)

          ! do Fourier Transform on rings

          if (keep_north(ithl) == 1) then
             ring(0:nphl-1) = local_map(istart_north:istart_north+nphl-1,1,imap) * &
                            & local_w8ring(ith-first_ring+1,1)*totWght( ithl)
             call s2hat_ring_analysis(nlmax, nmmax, ring, nphl, phring, kphi0(ithl))
             istart_north = istart_north+nphl
              
             ! store in a transposed order for the future alltoall
             do i = 0, nmmax
                ij = iring+indx(i)*nrings
                phas_n(ij,1) = phring(i)
             enddo 
          endif

          if (keep_south(ithl) == 1) then
             istart_south = istart_south-nphl
             ring(0:nphl-1) = local_map(istart_south:istart_south+nphl-1,1,imap) * &
                            & local_w8ring(ith-first_ring+1,1)*totWght( ithl)

             call s2hat_ring_analysis(nlmax, nmmax, ring, nphl, phring, kphi0(ithl))

              ! store in a transposed order for the future alltoall
              do i = 0, nmmax
                 ij = iring+indx(i)*nrings
                 phas_s(ij,1) = phring(i)
              enddo 
           endif

       enddo ! loop over local ring
  
   enddo ! loop over chunks

   deallocate( indx)
   deallocate( nph)
   deallocate( ring)
   deallocate( phring)
   deallocate( kphi0, totWght)
   deallocate( keep_north, keep_south)

   ! redistribute intermediate products over procs

   phas_n_ptr => phas_n
   phas_s_ptr => phas_s

   call redistInterProductsMap2Alm( nrings, nmvals, 1, nphas, phas_n_ptr, myid, numprocs, comm)
   call redistInterProductsMap2Alm( nrings, nmvals, 1, nphas, phas_s_ptr, myid, numprocs, comm)

   ! and get global info about all the rings

   allocate( keep_it_all(0:nringsall-1))
   keep_it_all(:) = 0
   call gatherRingInfo( first_ring, nrings, keep_it, nringsall, keep_it_all, myid, numprocs, comm)
   deallocate( keep_it)

   !-----------------------------------------------------------------------
   !              computes the a_lm by integrating over theta
   !                  lambda_lm(theta) * phas_m(theta)
   !                         for each m and l
   !-----------------------------------------------------------------------

   allocate(dalm(1:2,0:nlmax),stat = status)
   call assert_alloc(status,code,'dalm')

   nchunks   = nringsall/SMAXCHK + 1  ! number of chunks
   chunksize = (nringsall+nchunks-1)/nchunks ! <= SMAXCHK

   local_alms(:,:,:,imap) = 0.0 ! set the whole alm array to zero

   kring = 0
   do ichunk = 0, nchunks-1

      lchk  = ichunk * chunksize + 1
      uchk  = min(lchk+chunksize - 1, nringsall)

      mshift = 0

      do mindx = 0, nmvals-1

          mindx1 = mindx+1
	  
	  m = mvals( mindx)
 
          ! introduce double precision vector to perform summation over ith for each l
          dalm(1:2, m:nlmax ) = 0.0_dp

          jring = kring
          do ithl = 0, uchk-lchk 
             
             iring = lchk+ithl-1

             l_min = l_min_ylm(m, pixelization%sth(iring+1))

             if (keep_it_all(iring) == 1) then

                if (nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001) split in two (RS, 02-2007)

                   i_mm = n_lm * jring + mshift        ! location of Ym,m for ring ith
                   i_mm1 = i_mm+1

                   !           ---------- l = m ----------
                   par_lm = 1

                   php = phas_n(mindx*nringsobs+jring,1) + phas_s(mindx*nringsobs+jring,1) ! sum  (if (l+m) even)
                   phm = phas_n(mindx*nringsobs+jring,1) - phas_s(mindx*nringsobs+jring,1) ! diff (if (l+m) odd)
                   phas_sd(-1:0) =  (/ real(phm, kind=dp), aimag(phm) /)
                   phas_sd( 1:2) =  (/ real(php, kind=dp), aimag(php) /)
                
                   if (m >= l_min) then
                      dalm(1:2, m) = dalm(1:2, m) + local_plm_comp_a(i_mm1) * phas_sd(par_lm:par_lm+1)
                   endif

                   !           ---------- l > m ----------

                   l_start = max(m+1, l_min)
                   if (mod(m+l_start,2) == 0) l_start = l_start+1
                   do l = l_start, nlmax, 2
                      dalm(1,l) = dalm(1,l) + local_plm_comp_a(i_mm1+l-m) * phas_sd(-1)
                      dalm(2,l) = dalm(2,l) + local_plm_comp_a(i_mm1+l-m) * phas_sd(0)
                   enddo ! loop over l

                   l_start = max(m+2, l_min)
                   if (mod(m+l_start,2) == 1) l_start = l_start+1
                   do l = l_start, nlmax, 2
                      dalm(1,l) = dalm(1,l) + local_plm_comp_a(i_mm1+l-m) * phas_sd(1)
                      dalm(2,l) = dalm(2,l) + local_plm_comp_a(i_mm1+l-m) * phas_sd(2)
                   enddo ! loop over l

                endif ! test for l values

                jring = jring + 1

             endif ! test for the cut sky

          enddo ! loop over chunk rings (ithl)
		  
          do l = m, nlmax
             l1 = l+1
             local_alms_comp_a(l1,mindx1) = local_alms_comp_a(l1,mindx1) + &
                                          & cmplx(dalm(1, l), dalm(2, l), kind=DP)
          enddo

          mshift = mshift + nlmax-m+1

       enddo ! loop over m

       kring = jring ! store the number of obs rings so far       

    enddo ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate( keep_it_all)
    deallocate( dalm)
    deallocate( phas_n)
    deallocate( phas_s)

    return

  end subroutine s2hat_map2alm_sc_pre


  subroutine s2hat_map2alm_pol( pixelization, scan, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, local_w8ring, map_size, local_map, lda, local_alms, numprocs, myid, comm)

    !=======================================================================
    !     computes the a(l,m) from a polarized map for the HEALPIX pixelisation
    !     calculating spherical harmonics "on the fly"
    !=======================================================================

    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan

    integer (i4b) ::  nlmax, nmmax, nmvals, nmaps, first_ring, last_ring, lda, map_size, numprocs, myid, comm
    integer(i4b), dimension(0:nmvals-1) :: mvals
    real(dp), dimension(:,:), pointer :: local_w8ring
    real(dp), dimension(:,:,:), pointer :: local_map	

    ! output
    complex (dp), dimension(:,:,:,:), pointer :: local_alms

    ! internal
    integer(I4B) :: l, l1, m, ith, ithl, scalel, scalem, nrings, nphmx, mindx, mindx1
    integer(I8B) :: istart_south, istart_north, istart_south_save, istart_north_save, npix

    integer(I4B) :: par_lm, imap, ichunk, lchk, uchk, chunksize, nchunks
    real(DP)     :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)     :: lambda_w, lambda_x, normal_m, lam_lm1m
    real(DP)     :: fm2, fl, flm1, cth0, sth0, one_on_s2, c_on_s2
    real(DP)     :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP)     :: ovflow, unflow
    real(DP),     dimension(-3:8,1:nmaps)  :: phas_sd
    real(DP),     dimension(-3:6,1:nmaps)  :: phas_sdx
    real(DP),     dimension(:,:), allocatable :: recfac
    real(DP),     dimension(:,:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: lam_fact, mfac

    integer(i4b)  :: l_min
    real(DP),     dimension(:),   allocatable :: normal_l
    complex(DPC)  :: php, phm
    complex(DPC), dimension(:,:), pointer :: tmphas_n_ptr, tmphas_s_ptr
    complex(DPC), dimension(:,:), allocatable, target :: tmphas_n, tmphas_s
    complex(DPC), dimension(:), allocatable :: phring
    complex(DPC), dimension(:,:,:), allocatable :: phas_ns
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)    :: nphl, i, ii, ij, j, iring, jring, kring, nringsall, nringsobs, nphas
    integer(i4b)    :: l_start
    real(DP),     dimension(:), allocatable :: cth, sth, kphi0, totWght
    integer(i4b), dimension(:), allocatable :: keep_it, keep_north, keep_south, keep_it_all
    integer(i4b), dimension(:), allocatable :: indx, nph

    character(LEN=*), parameter :: code = 'MAP2ALM'
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

    nringsall = pixelization%nringsall
    nringsobs = scan%nringsobs
    npix   = pixelization%npixsall
    nphmx  = pixelization%nphmx

    !     --- allocates space for arrays ---
	
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    nphas = max( (nmmax+1)*nrings, nmvals*nringsobs)

    allocate(phas_ns(0:5,0:nphas-1,1:nmaps),stat = status)
    call assert_alloc(status,code,'phas_ns')

    allocate(tmphas_n(0:nphas-1,1:nmaps),stat = status)
    call assert_alloc(status,code,'tmphas_n')

    allocate(tmphas_s(0:nphas-1,1:nmaps),stat = status)
    call assert_alloc(status,code,'tmphas_s')

    allocate(keep_it(0:nrings-1),stat=status) 
    call assert_alloc(status,code,'keep_it')

    allocate( kphi0(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'kphi0')

    allocate( totWght(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'totWght')

    allocate( keep_north(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'keep_north')

    allocate( keep_south(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'keep_south')

    allocate( nph(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'nph')

    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')
	
    allocate(phring(0:nmmax),stat = status)	   
    call assert_alloc(status,code,'phring')		
	
    allocate( indx( 0:nmmax), stat = status)
    call assert_alloc(status,code,'indx')

    call fix_mmap( 1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

    !     ------------ initiate variables and arrays ----------------

    ! loop over Stokes params	
 
   do i = 1, 3
	
       istart_north = 0
       istart_south = map_size
       tmphas_n(:,:) = 0_dpc
       tmphas_s(:,:) = 0_dpc	   

       do ichunk = 0, nchunks-1

          lchk = first_ring + ichunk * chunksize 
          uchk = min(lchk+chunksize - 1, last_ring)

          do ith = lchk, uchk
             ithl = ith - lchk !local index

             ! extract pixel location information
             nph(ithl) = pixelization%nph(ith)
             cth0 = pixelization%cth(ith)
             sth0 = pixelization%sth(ith)
             kphi0(ithl) = pixelization%kphi(ith)
             totWght(ithl) = pixelization%parea(ith)*pixelization%qwght(ith)

             ! find out which rings are to be analysed
             keep_north( ithl) = scan%nfl( ith)
             keep_south( ithl) = scan%sfl( ith)
             keep_it( ith-first_ring) = scan%fl( ith)             

          enddo

          !-----------------------------------------------------------------------
          !  computes the integral in phi : phas_m(theta)
          !  for each parallel from Pole to Equator (use symmetry for S. Hemisphere)
          !-----------------------------------------------------------------------

          istart_north_save = istart_north
          istart_south_save = istart_south

          do imap = 1, nmaps

             istart_north = istart_north_save
             istart_south = istart_south_save

             ! do Fourier Transform on rings
             do ith = lchk, uchk

                ithl = ith - lchk !local index
	        iring = ith - first_ring

                nphl = nph(ithl)

                if (keep_north(ithl) == 1) then
                   ring(0:nphl-1) = local_map(istart_north:istart_north+nphl-1,i,imap) * &
                                  & local_w8ring(iring+1,i)*totWght(ithl)
                   call s2hat_ring_analysis(nlmax, nmmax, ring, nphl, phring, kphi0(ithl))
			  
                   do j = 0, nmmax
                      ij = iring+indx(j)*nrings
                      tmphas_n(ij,imap) = phring(j)
	           enddo

                   istart_north = istart_north+nphl
                endif

                if (keep_south(ithl) == 1) then
                   istart_south = istart_south-nphl

                   ring(0:nphl-1) = local_map(istart_south:istart_south+nphl-1,i,imap) * &
                                  & local_w8ring(iring+1,i)*totWght(ithl)
                   call s2hat_ring_analysis(nlmax, nmmax, ring, nphl, phring, kphi0(ithl))

                   do j = 0, nmmax
                      ij = iring+indx(j)*nrings 
                      tmphas_s(ij,imap) = phring(j)
	           enddo
                endif
             enddo ! over rings of chunks

          enddo ! over maps

       enddo ! over chunks

       ! redistribute intermediate products over procs

       tmphas_n_ptr => tmphas_n
       tmphas_s_ptr => tmphas_s

       call redistInterProductsMap2Alm( nrings, nmvals, nmaps, nphas, tmphas_n_ptr, myid, numprocs, comm)
       call redistInterProductsMap2Alm( nrings, nmvals, nmaps, nphas, tmphas_s_ptr, myid, numprocs, comm)

       do imap = 1, nmaps
          do j = 0, nphas-1
	     phas_ns(2*(i-1), j, imap) = tmphas_n(j, imap)   ! north cap
	     phas_ns(2*i-1, j, imap) = tmphas_s(j, imap)     ! south cap
          enddo
       enddo

    enddo ! loop over Stokes 

    deallocate( tmphas_n, tmphas_s)
    deallocate( keep_north, keep_south)
    deallocate( kphi0, totWght)
    deallocate( nph)
    deallocate( ring)
    deallocate( phring)
    deallocate( indx)

    ! and get global info about all the rings
  
    allocate( keep_it_all(0:nringsall-1))
    keep_it_all(:) = 0
    call gatherRingInfo( first_ring, nrings, keep_it, nringsall, keep_it_all, myid, numprocs, comm)
    deallocate( keep_it)

    !-----------------------------------------------------------------------
    !              computes the a_lm by integrating over theta
    !                  lambda_lm(theta) * phas_m(theta)
    !                         for each m and l
    !-----------------------------------------------------------------------

    allocate(normal_l(0:nlmax),stat = status)
    call assert_alloc(status,code,'normal_l')

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    ! generate Polarization normalisation, etc
    call gen_mfac(nmmax,mfac)
    call gen_normpol(nlmax, normal_l)

    call init_rescale_local()   ! initialize rescale table

    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)

    allocate(recfac(0:1,0:nlmax),stat = status)
    call assert_alloc(status,code,'recfac')

    allocate(dalm(0:5,1:nmaps,0:nlmax), stat = status)
    call assert_alloc(status,code,'dalm')

    allocate(lam_fact(0:nlmax),stat = status)
    call assert_alloc(status,code,'lam_fact')

    nchunks   = nringsall/SMAXCHK + 1  ! number of chunks
    chunksize = (nringsall+nchunks-1)/nchunks ! <= SMAXCHK

    local_alms = 0.0            ! set the whole alm array to zero

    kring = 0
    do ichunk = 0, nchunks-1

       lchk  = ichunk * chunksize + 1
       uchk  = min(lchk+chunksize - 1, nringsall)

       do mindx = 0, nmvals-1

          mindx1 = mindx+1 
          m = mvals( mindx)  ! actual m value

          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)
          ! generate Ylm relation factor for degree m
          call gen_lamfac(nlmax, m, lam_fact)
          fm2 = real(m * m, kind = DP)
          normal_m = (2.0_dp * m) * (1 - m)

          ! introduce double precision vector to perform summation over ith for each l
          dalm(0:5, 1:nmaps, m:nlmax) = 0.0_dp

          jring = kring
          do ithl = 0, uchk-lchk

             iring = ithl+lchk-1

             cth_ring = pixelization%cth(iring+1)
             one_on_s2 =  1.0_dp / pixelization%sth(iring+1)**2
             c_on_s2 = cth_ring*one_on_s2
             two_on_s2 = 2.0_dp * one_on_s2
             two_cth_ring = 2.0_dp * cth_ring
             b_w =  c_on_s2

             l_min = l_min_ylm(m, pixelization%sth(iring+1)) ! lower l bound of non-negligible Ylm

             if (keep_it_all(iring) == 1) then

                if (nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001) split in two (RS, 02-2007)

                   ! determine lam_mm
                   call mpi_compute_lam_mm(mfac(m), m, pixelization%sth(iring+1), lam_mm, corfac, scalem)

                   do imap = 1, nmaps
                      do i=0,2 ! loop on T, Q, U
                         phm = phas_ns(2*i,mindx*nringsobs+jring, imap) - phas_ns(2*i+1,mindx*nringsobs+jring, imap) ! N-S: (l+m) odd
                         phas_sd(-3+2*i, imap)   = real(phm, kind=dp)
                         phas_sd(-3+2*i+1, imap) = aimag(phm)
                         php = phas_ns(2*i,mindx*nringsobs+jring,imap) + phas_ns(2*i+1,mindx*nringsobs+jring,imap) ! N+S: (l+m) even
                         phas_sd( 3+2*i, imap)   = real(php, kind=dp) 
                         phas_sd( 3+2*i+1, imap) = aimag(php)
                      enddo
                      phas_sdx(-3: 0, imap) = (/ phas_sd(8, imap), - phas_sd(7,imap), - phas_sd(6,imap),   phas_sd(5,imap) /)
                      phas_sdx( 3: 6, imap) = (/ phas_sd(2, imap), - phas_sd(1,imap), - phas_sd(0,imap),   phas_sd(-1,imap) /)
                   enddo

                   !           ---------- l = m ----------

                   par_lm = 3  ! = 3 * (-1)^(l+m)
                   lam_lm = corfac * lam_mm            ! Actual lam_mm

                   if (m >= l_min) then
                      lambda_w =  (normal_m * lam_lm)  * ( 0.5_dp - one_on_s2 )
                      lambda_x =  (normal_m * lam_lm)  * c_on_s2

                      ! alm_T = \int T Ylm
                      ! alm_G = \int ( - Q Wlm - i U Xlm )
                      ! alm_C = \int ( - U Wlm + i Q Xlm )
                      dalm(0:1, 1:nmaps, m) = dalm(0:1, 1:nmaps, m) + lam_lm * phas_sd(par_lm:par_lm+1,1:nmaps)
                      dalm(2:5, 1:nmaps, m) = dalm(2:5, 1:nmaps, m)  &
                          &       - lambda_w * phas_sd(par_lm+2:par_lm+5, 1:nmaps) &
                          &       + lambda_x * phas_sdx(par_lm:par_lm+3, 1:nmaps)
                   endif

	           !           ---------- l > m ----------

                   lam_0 = 0.0_dp
                   lam_1 = 1.0_dp
                   scalel = scalem
                   lam_2 = cth_ring * lam_1 * recfac(0,m)
                   fm_on_s2 = m * one_on_s2

                   do l = m+1, l_min-1

                      lam_0 = lam_1 * recfac(1,l-1)
                      lam_1 = lam_2
                      lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)
                      if (abs(lam_2) > OVFLOW) then
                         lam_1 = lam_1*UNFLOW
                         lam_2 = lam_2*UNFLOW
                         scalel= scalel + 1
            !             corfac = rescale_tab(max(scalel+scalem,RSMIN))

                         if( (scalel) > RSMIN) then
                            corfac = rescale_tab( scalel)
                         else
                            corfac = rescale_tab( RSMIN)
                         endif

                     elseif (abs(lam_2) < UNFLOW) then
                         lam_1 = lam_1*OVFLOW
	                 lam_2 = lam_2*OVFLOW
                         scalel= scalel - 1
            !             corfac = rescale_tab(max(scalel+scalem,RSMIN))

                         if( (scalel) > RSMIN) then
                            corfac = rescale_tab( scalel)
                         else
                            corfac = rescale_tab( RSMIN)
                         endif

                      endif   

      	           enddo ! loop over l up to lmin ...

                   ! l = max( m+1, l_min)
                   if( (m+1) >= l_min) then
                       l = m+1
                   else
                       l = l_min
                   endif

                   ! par_lm = -3 * (-1)^(l+m)
                   if( 2*((l+m)/2) == (l+m)) then 
                       par_lm = -3
                   else
                       par_lm = 3
                   endif

                   do l = l, nlmax

                      par_lm = - par_lm  ! = 3 * (-1)^(l+m)

                      lam_lm1m = lam_lm * lam_fact(l)
                      lam_lm = lam_2 * corfac * lam_mm

                      fl = real(l, kind = DP)
                      flm1 = fl - 1.0_dp
                      a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                      a_x =  two_cth_ring * flm1
                      lambda_w = b_w * lam_lm1m - a_w * lam_lm 
                      lambda_x = fm_on_s2 * (lam_lm1m - a_x * lam_lm)

                      dalm(0:1, 1:nmaps, l) = dalm(0:1, 1:nmaps, l) + lam_lm * phas_sd(par_lm:par_lm+1, 1:nmaps)
                      dalm(2:5, 1:nmaps, l) = dalm(2:5, 1:nmaps, l)  &
                            &       - lambda_w * phas_sd(par_lm+2:par_lm+5, 1:nmaps) &
                            &       + lambda_x * phas_sdx(par_lm:par_lm+3, 1:nmaps)

                      lam_0 = lam_1 * recfac(1,l-1)
                      lam_1 = lam_2
                      lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)
                      if (abs(lam_2) > OVFLOW) then
                         lam_1 = lam_1*UNFLOW
                         lam_2 = lam_2*UNFLOW
                         scalel= scalel + 1
     !                    corfac = rescale_tab(max(scalel+scalem,RSMIN))

                         if( scalel > RSMIN) then
                            corfac = rescale_tab( scalel)
                         else
                            corfac = rescale_tab( RSMIN)
                         endif

                     elseif (abs(lam_2) < UNFLOW) then
                         lam_1 = lam_1*OVFLOW
	                 lam_2 = lam_2*OVFLOW
                         scalel= scalel - 1
    !                     corfac = rescale_tab(max(scalel+scalem,RSMIN))

                         if( scalel > RSMIN) then
                            corfac = rescale_tab( scalel)
                         else
                            corfac = rescale_tab( RSMIN)
                         endif

                      endif   
      	           enddo ! loop over l

                 endif ! test for l values

                jring = jring + 1

             endif ! test for cut sky

          enddo ! over the rings of the chunk

          do imap = 1, nmaps
             do l = m, nlmax
                l1 = l+1

                local_alms_comp_a(l1,mindx1,imap) = local_alms_comp_a(l1,mindx1,imap) + cmplx(dalm(0, imap, l), dalm(1, imap, l))
                local_alms_comp_b(l1,mindx1,imap) = local_alms_comp_b(l1,mindx1,imap) + &
                                                  & cmplx(dalm(2, imap, l), dalm(3, imap, l)) * normal_l(l)
                local_alms_comp_c(l1,mindx1,imap) = local_alms_comp_c(l1,mindx1,imap) + &
                                                  & cmplx(dalm(4, imap, l), dalm(5, imap, l)) * normal_l(l)
             enddo
          enddo

	enddo ! loop over m

        kring = jring ! store the number of obs rings so far

    enddo    ! loop over chunks

    !     --------------------
    !     free memory and exit
    !     --------------------

    deallocate (recfac, dalm, lam_fact)
    deallocate( mfac, normal_l)
    deallocate( phas_ns)
    deallocate( keep_it_all)

    return

  end subroutine s2hat_map2alm_pol

  subroutine s2hat_map2alm_spin( pixelization, scan, spin_in, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, local_w8ring, map_size, local_map, lda, local_alms, numprocs, myid, comm)

    !=================================================================================================
    !
    !     computes the a(l,m) for a complex field of spin s.
    !
    !     the input field(s) is in a form of a real valued array [*, 1:2, nmaps] with real part 
    !     corresponding to the index 1 and the imaginary part to 2.
    !
    !     The complex field is then assumed to be :    re + i * (-conv_sign) * im.
    !
    !     by default the convention is chosen so SPIN_CONV_SIGN = -1 (Healpix convention)
    !     that can be changed only by hand (by changing the variable SPIN_CONV_SIGN) prior to compilation.
    !
    !     Note that whatever the input spin sign the routine actually deals always with
    !     fields of the positive spin. Given the convention above for negative spin s (<0) 
    !     the corresponding (-s) (>0) field is just:  re + i * SPIN_CONV_SIGN * im.
    !
    !     The dimensions of the output array are either: [ nlmax, nmvals, 2, nmaps] -> lda == nlmax;
    !     or                                             [ 2, nlmax, nmvals, nmaps] -> lda == 2
    !
    !     The spin spherical harmonics are calculated "on the fly".
    !
    !     rs@apc - sept 2007
    !
    !=================================================================================================

    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan

    integer (i4b), intent(in) ::  spin_in, nlmax, nmmax, nmvals, nmaps, first_ring, last_ring, lda, map_size, numprocs, myid, comm
    integer(i4b), dimension(0:nmvals-1) :: mvals
    real(dp), dimension(:,:), pointer :: local_w8ring
    real(dp), dimension(:,:,:), pointer :: local_map	

    ! output
    complex (dp), dimension(:,:,:,:), pointer :: local_alms

    ! internal
    integer(I4B) :: l, l1, m, ith, ithl, nrings, nphmx, mindx, mindx1
    integer(i4b), dimension( 0:1) :: scalel, scalem
    integer(i4b), dimension( 1:2) :: compSign
    integer(I8B) :: istart_south, istart_north, istart_south_save, istart_north_save, npix

    integer(I4B) :: par_lm, imap, ispin, ichunk, lchk, uchk, chunksize, nchunks, ss, mm
    real(DP)     :: cth_ring, cth0, sth0, spinPow, signFix, mfac0
    real(DP)     :: ovflow, unflow
    real(DP),     dimension(0:1) :: corfac, lam_0, lam_1, lam_2, lam_mm, lam_lm, lam_comb, fact1
    real(DP),     dimension(0:7)  :: phas_sd
    real(DP),     dimension(0:7,1:nmaps,-1:1)  :: phas_sdx
    real(DP),     dimension(:,:), allocatable :: recfac
    real(DP),     dimension(:,:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: mfac, sfac

    integer(i4b)  :: l_min
    ! real(DP),     dimension(:),   allocatable :: normal_l
    complex(DPC) :: php, phm
    complex(DPC), dimension(:,:), pointer :: tmphas_n_ptr, tmphas_s_ptr
    complex(DPC), dimension(:,:), allocatable, target :: tmphas_n, tmphas_s
    complex(DPC), dimension(:), allocatable :: phring
    complex(DPC), dimension(:,:,:), allocatable :: phas_ns
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)    :: spin, nphl, i, ii, ij, j, iring, jring, kring, nringsall, nringsobs, nphas
    integer(i4b)    :: l_start
    real(DP),     dimension(:), allocatable :: cth, sth, kphi0, totWght
    integer(i4b), dimension(:), allocatable :: keep_it, keep_north, keep_south, keep_it_all
    integer(i4b), dimension(:), allocatable :: indx, nph

    character(LEN=*), parameter :: code = 'MAP2ALM'
    integer(I4B) :: status, ierr

    complex(DPC), dimension(:,:,:), pointer :: local_alms_comp_a, local_alms_comp_b

    !=======================================================================

    spin = abs( spin_in)
    spinPow = (-1.d0) ** spin

    ! define how the input real maps are to be combined into a complex field of |s| (>0) spin

    compSign(1) = 1
    compSign(2) = -SPIN_CONV_SIGN

    if( lda == 2) then       ! i.e., healpix convention for the alm array

       local_alms_comp_a => local_alms(1,0:nlmax,0:nmvals-1,1:nmaps)
       local_alms_comp_b => local_alms(2,0:nlmax,0:nmvals-1,1:nmaps)

    else                     ! i.e., s2hat convention for the alm array

       local_alms_comp_a => local_alms(0:nlmax,0:nmvals-1,1,1:nmaps)
       local_alms_comp_b => local_alms(0:nlmax,0:nmvals-1,2,1:nmaps)

    end if

    nrings = (last_ring-first_ring+1)

    nringsall = pixelization%nringsall
    nringsobs = scan%nringsobs
    npix   = pixelization%npixsall
    nphmx  = pixelization%nphmx

    !     --- allocates space for arrays ---
	
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    nphas = max( (nmmax+1)*nrings, nmvals*nringsobs)

    allocate(phas_ns(0:5,0:nphas-1,1:nmaps),stat = status)
    call assert_alloc(status,code,'phas_ns')

    allocate(tmphas_n(0:nphas-1,1:nmaps),stat = status)
    call assert_alloc(status,code,'tmphas_n')

    allocate(tmphas_s(0:nphas-1,1:nmaps),stat = status)
    call assert_alloc(status,code,'tmphas_s')

    allocate(keep_it(0:nrings-1),stat=status) 
    call assert_alloc(status,code,'keep_it')

    allocate( kphi0(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'kphi0')

    allocate( totWght(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'totWght')

    allocate( keep_north(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'keep_north')

    allocate( keep_south(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'keep_south')

    allocate( nph(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'nph')

    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')
	
    allocate(phring(0:nmmax),stat = status)	   
    call assert_alloc(status,code,'phring')		
	
    allocate( indx( 0:nmmax), stat = status)
    call assert_alloc(status,code,'indx')

    call fix_mmap( 1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

    !     ------------ initiate variables and arrays ----------------

    ! loop over spin components
 
   do i = 1, 2
	
       istart_north = 0
       istart_south = map_size
       tmphas_n(:,:) = 0_dpc
       tmphas_s(:,:) = 0_dpc	   

       do ichunk = 0, nchunks-1

          lchk = first_ring + ichunk * chunksize 
          uchk = min(lchk+chunksize - 1, last_ring)

          do ith = lchk, uchk
             ithl = ith - lchk !local index

             ! extract pixel location information
             nph(ithl) = pixelization%nph(ith)
             cth0 = pixelization%cth(ith)
             sth0 = pixelization%sth(ith)
             kphi0(ithl) = pixelization%kphi(ith)
             totWght(ithl) = pixelization%parea(ith)*pixelization%qwght(ith)

             ! find out which rings are to be analysed
             keep_north( ithl) = scan%nfl( ith)
             keep_south( ithl) = scan%sfl( ith)
             keep_it( ith-first_ring) = scan%fl( ith)             

          enddo

          !-----------------------------------------------------------------------
          !  computes the integral in phi : phas_m(theta)
          !  for each parallel from Pole to Equator (use symmetry for S. Hemisphere)
          !-----------------------------------------------------------------------

          istart_north_save = istart_north
          istart_south_save = istart_south

          do imap = 1, nmaps

             istart_north = istart_north_save
             istart_south = istart_south_save

             ! do Fourier Transform on rings
             do ith = lchk, uchk

                ithl = ith - lchk !local index
	        iring = ith - first_ring

                nphl = nph(ithl)

                if (keep_north(ithl) == 1) then
                   ring(0:nphl-1) = compSign(i) * local_map(istart_north:istart_north+nphl-1,i,imap) * &
                                  & local_w8ring(iring+1,i)*totWght(ithl)
                   call s2hat_ring_analysis(nlmax, nmmax, ring, nphl, phring, kphi0(ithl))
			  
                   do j = 0, nmmax
                      ij = iring+indx(j)*nrings
                      tmphas_n(ij,imap) = phring(j)
	           enddo

                   istart_north = istart_north+nphl
                endif

                if (keep_south(ithl) == 1) then
                   istart_south = istart_south-nphl

                   ring(0:nphl-1) = compSign(i) * local_map(istart_south:istart_south+nphl-1,i,imap) * &
                                  & local_w8ring(iring+1,i)*totWght(ithl)
                   call s2hat_ring_analysis(nlmax, nmmax, ring, nphl, phring, kphi0(ithl))

                   do j = 0, nmmax
                      ij = iring+indx(j)*nrings 
                      tmphas_s(ij,imap) = phring(j)
	           enddo
                endif
             enddo ! over rings of chunks

          enddo ! over maps

       enddo ! over chunks

       ! redistribute intermediate products over procs

       tmphas_n_ptr => tmphas_n
       tmphas_s_ptr => tmphas_s

       call redistInterProductsMap2Alm( nrings, nmvals, nmaps, nphas, tmphas_n_ptr, myid, numprocs, comm)
       call redistInterProductsMap2Alm( nrings, nmvals, nmaps, nphas, tmphas_s_ptr, myid, numprocs, comm)

       do imap = 1, nmaps
          do j = 0, nphas-1
	     phas_ns(2*(i-1), j, imap) = tmphas_n(j, imap)   ! north cap
	     phas_ns(2*i-1, j, imap) = tmphas_s(j, imap)     ! south cap
          enddo
       enddo

    enddo ! loop over Stokes 

    deallocate( tmphas_n, tmphas_s)
    deallocate( keep_north, keep_south)
    deallocate( kphi0, totWght)
    deallocate( nph)
    deallocate( ring)
    deallocate( phring)
    deallocate( indx)

    ! and get global info about all the rings
  
    allocate( keep_it_all(0:nringsall-1))
    keep_it_all(:) = 0
    call gatherRingInfo( first_ring, nrings, keep_it, nringsall, keep_it_all, myid, numprocs, comm)
    deallocate( keep_it)

    !-----------------------------------------------------------------------
    !              computes the a_lm by integrating over theta
    !                  lambda_lm(theta) * phas_m(theta)
    !                         for each m and l
    !-----------------------------------------------------------------------

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    ! generate Polarization normalisation, etc
    call s2hat_gen_mfac(spin,nmmax,mfac)

    allocate(sfac(0:spin),stat = status)
    call assert_alloc(status,code,'sfac')

    ! generate Polarization normalisation, etc
    call s2hat_gen_sfac( spin, sfac)

    call init_rescale_local()   ! initialize rescale table

    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)

    allocate(recfac(0:1,0:nlmax),stat = status)
    call assert_alloc(status,code,'recfac')

    allocate(dalm(0:3, 1:nmaps, 0:nlmax), stat = status)
    call assert_alloc(status,code,'dalm')

    nchunks   = nringsall/SMAXCHK + 1  ! number of chunks
    chunksize = (nringsall+nchunks-1)/nchunks ! <= SMAXCHK

    local_alms = 0.0            ! set the whole alm array to zero

    kring = 0
    do ichunk = 0, nchunks-1

       lchk  = ichunk * chunksize + 1
       uchk  = min(lchk+chunksize - 1, nringsall)

       do mindx = 0, nmvals-1
 
          mindx1 = mindx+1
          m = mvals( mindx)  ! actual m value

          if( m >= spin) then      ! the cases with m < spin are treated below

              signFix = 1.0
              mm = m               ! values for the recurrence
              ss = spin

              mfac0 = mfac( m)

          else         

              signFix = (-1.0)**(spin+m)
              mm = spin  ! m <-> s interchange
              ss = m

              mfac0 = sfac( m)

          endif

          ! generate recursion factors (recfac) for sYlm of degree m
          call s2hat_gen_recfac( ss, nlmax, mm, recfac)

          ! introduce double precision vector to perform summation over ith for each l
          dalm(0:3, 1:nmaps, mm:nlmax) = 0.0_dp

          jring = kring
          do ithl = 0, uchk-lchk

             iring = ithl+lchk-1

             cth_ring = pixelization%cth(iring+1)

             !!! l_min = l_min_ylm( m, pixelization%sth(iring+1)) ! lower l bound of non-negligible Ylm

             l_min = mm ! either m or spin - no optimization for time being here ...

             if (keep_it_all(iring) == 1) then

                if (nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001) split in two (RS, 02-2007)

                   ! determine lam_mm
                   call mpi_compute_slam_mm( mfac0, ss, mm, pixelization%cth(iring+1), pixelization%sth(iring+1), lam_mm(0), corfac(0), scalem(0))
                   call mpi_compute_slam_mm( mfac0, -ss, mm, pixelization%cth(iring+1), pixelization%sth(iring+1), lam_mm(1), corfac(1), scalem(1))

                   do imap = 1, nmaps
                      do i=0,1 ! loop over spin components
                         phm = phas_ns(2*i,mindx*nringsobs+jring, imap) - spinPow * phas_ns(2*i+1,mindx*nringsobs+jring, imap) ! N-S:
                         phas_sd( 4*i)   = real(phm, kind=dp)
                         phas_sd( 4*i+1) = aimag(phm)

                         php = phas_ns(2*i,mindx*nringsobs+jring,imap) + spinPow * phas_ns(2*i+1,mindx*nringsobs+jring,imap)   ! N+S:
                         phas_sd( 4*i+2) = real(php, kind=dp) 
                         phas_sd( 4*i+3) = aimag(php)
                      enddo

                      ! even l+m data
                      phas_sdx(0:7, imap, -1) = (/ phas_sd(2), phas_sd(3), -SPIN_CONV_SIGN*phas_sd(6), -SPIN_CONV_SIGN*phas_sd(7), SPIN_CONV_SIGN*phas_sd(5), -SPIN_CONV_SIGN*phas_sd(4), phas_sd(1), -phas_sd(0) /)
                      ! odd l+m data
                      phas_sdx(0:7, imap,  1) = (/ phas_sd(0), phas_sd(1), -SPIN_CONV_SIGN*phas_sd(4), -SPIN_CONV_SIGN*phas_sd(5), SPIN_CONV_SIGN*phas_sd(7), -SPIN_CONV_SIGN*phas_sd(6), phas_sd(3), -phas_sd(2) /)

                   enddo

                   !           ---------- l = m ----------
                 
                   par_lm = -signFix
                   lam_lm(0) = signFix * corfac(0) * lam_mm(0)            ! Actual slam_mm with s > 0
                   lam_lm(1) = corfac(1) * lam_mm(1)                      ! Actual slam_mm with s < 0

                   lam_comb(0) = lam_lm(0) + spinPow * lam_lm(1)
                   lam_comb(1) = lam_lm(0) - spinPow * lam_lm(1)

                   if (mm >= l_min) then
                       dalm(0:3, 1:nmaps, mm) = dalm(0:3, 1:nmaps, mm)  &
                               &     + lam_comb(0) * phas_sdx( 0:3, 1:nmaps, par_lm) &
                               &     + lam_comb(1) * phas_sdx( 4:7, 1:nmaps, par_lm)
                   endif

	           !           ---------- l > m ----------

                   lam_0(0:1) = 0.0_dp
                   lam_1(0:1) = 1.0_dp

                   scalel(0:1) = 0

                   lam_2(0) = (cth_ring + ss/dble(mm+1)) * lam_1(0) * recfac(0,mm)   ! spin > 0, l = m
                   lam_2(1) = (cth_ring - ss/dble(mm+1)) * lam_1(1) * recfac(0,mm)   ! spin < 0, l = m

                   do l = mm+1, nlmax

                      par_lm = - par_lm

                      fact1(0) = cth_ring + ss*mm/dble(l+1)/dble(l)
                      fact1(1) = cth_ring - ss*mm/dble(l+1)/dble(l)

                      if (l >= l_min) then

                         lam_lm(0) = signFix * lam_2(0) * corfac(0) * lam_mm(0)
                         lam_lm(1) = lam_2(1) * corfac(1) * lam_mm(1)

                         lam_comb(0) = lam_lm(0) + spinPow * lam_lm(1)
                         lam_comb(1) = lam_lm(0) - spinPow * lam_lm(1)

                         dalm(0:3, 1:nmaps, l) = dalm(0:3, 1:nmaps, l)  &
                                 &       + lam_comb(0) * phas_sdx(0:3, 1:nmaps, par_lm) &
                                 &       + lam_comb(1) * phas_sdx(4:7, 1:nmaps, par_lm)
                      endif

                      ! actual sPlm reccurence

                      do ispin = 0, 1

                         lam_0(ispin) = lam_1(ispin) * recfac(1,l-1)
                         lam_1(ispin) = lam_2(ispin)
                         lam_2(ispin) = (fact1(ispin) * lam_1(ispin) - lam_0(ispin)) * recfac(0,l)

                         ! and rescaling if needed
                        
                         if (abs(lam_2(ispin)) > OVFLOW) then
                            lam_1(ispin) = lam_1(ispin)*UNFLOW
                            lam_2(ispin) = lam_2(ispin)*UNFLOW
                            scalel(ispin) = scalel(ispin) + 1
                            corfac(ispin) = rescale_tab(max(scalel(ispin)+scalem(ispin),RSMIN))
                         elseif (abs(lam_2(ispin)) < UNFLOW) then
                            lam_1(ispin) = lam_1(ispin)*OVFLOW
	                    lam_2(ispin) = lam_2(ispin)*OVFLOW
                            scalel(ispin) = scalel(ispin) - 1
                            corfac(ispin) = rescale_tab(max(scalel(ispin)+scalem(ispin),RSMIN))
                         endif   

                      enddo ! loop over ispin

      	           enddo ! loop over l

                 endif ! test for l values

                 jring = jring + 1

                 endif ! test for cut sky

              enddo ! over the rings of the chunk

              do imap = 1, nmaps
                 do l = mm, nlmax
                    l1 = l+1

                    local_alms_comp_a(l1,mindx1,imap) = local_alms_comp_a(l1,mindx1,imap) + &
                                                      & SPIN_CONV_SIGN * cmplx(dalm(0, imap, l), dalm(1, imap, l))/2.0
                    local_alms_comp_b(l1,mindx1,imap) = local_alms_comp_b(l1,mindx1,imap) + &
                                                      & SPIN_CONV_SIGN * cmplx(dalm(2, imap, l), dalm(3, imap, l))/2.0
                 enddo
              enddo

	enddo ! loop over m

        kring = jring ! store the number of obs rings so far

    enddo    ! loop over chunks

    !     --------------------
    !     free memory and exit
    !     --------------------

    deallocate (recfac, dalm)
    deallocate( mfac, sfac)
    deallocate( phas_ns)
    deallocate( keep_it_all)

    return

  end subroutine s2hat_map2alm_spin

  subroutine s2hat_map2alm_pol_pre1( pixelization, scan, nlmax, nmmax, nmvals, mvals, imap, first_ring, last_ring, local_w8ring, map_size, local_map, &
                                   & lda, nplm, local_plm, local_alms, numprocs, myid, comm)

    !=======================================================================
    !     computes the a(l,m) from a polarized map for the HEALPIX pixelisation
    !     using precomputed and stored scalar Ylm_T(theta)
    !=======================================================================

    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan

    integer (i4b) :: nlmax, nmmax, nmvals, imap, first_ring, last_ring, lda, map_size, numprocs, myid, comm
    integer(i4b), dimension(0:nmvals-1) :: mvals
    real(dp), dimension(:,:), pointer :: local_w8ring
    real(dp), dimension(:,:,:), pointer :: local_map
    integer(i8b) :: nplm
    real(dp), dimension(:,:), pointer :: local_plm	

    ! output
    complex (dp), dimension(:,:,:,:), pointer :: local_alms

    ! internal
    integer(I4B) :: l, l1, m, ith, ithl, scalel, scalem, nrings, nphmx, mindx, mindx1
    integer(I8B) :: istart_south, istart_north, npix, i_mm, i_mm1, mshift

    integer(I4B)  :: par_lm, ichunk, lchk, uchk, chunksize, nchunks
    real(DP)      :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)      :: lambda_w, lambda_x, normal_m, lam_lm1m
    real(DP)      :: fm2, fl, flm1, cth0, sth0, one_on_s2, c_on_s2
    real(DP)      :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP),     dimension(-3:8)             :: phas_sd
    real(DP),     dimension(-3:6)             :: phas_sdx
    real(DP),     dimension(:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: lam_fact

    integer(i4b)                              :: l_min
    real(DP),     dimension(:),   allocatable :: normal_l
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:), pointer :: tmphas_n_ptr, tmphas_s_ptr
    complex(DPC), dimension(:,:), allocatable, target :: tmphas_n, tmphas_s
    complex(DPC), dimension(:), allocatable :: phring
    complex(DPC), dimension(:,:), allocatable :: phas_ns
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)    :: nphl, i, ii, ij, j, iring, jring, kring, nringsall, nringsobs, nphas
    integer(i4b)    :: l_start
    integer(I8B)    :: startpix, n_lm
    real(DP),     dimension(:), allocatable :: cth, sth, kphi0, totWght
    integer(i4b), dimension(:), allocatable :: keep_it, keep_north, keep_south, keep_it_all
    integer(i4b), dimension(:), allocatable :: indx, nph

    character(LEN=*), parameter :: code = 'MAP2ALM'
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

    nringsall = pixelization%nringsall
    nringsobs = scan%nringsobs
    npix   = pixelization%npixsall
    nphmx  = pixelization%nphmx

    n_lm = nummmodes( nlmax, nmvals, mvals)

    !     --- allocates space for arrays ---
	
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    nphas = max( (nmmax+1)*nrings, nmvals*nringsobs)

    allocate(phas_ns(0:5,0:nphas-1),stat = status)
    call assert_alloc(status,code,'phas_ns')

    allocate(tmphas_n(0:nphas-1,1:1),stat = status)
    call assert_alloc(status,code,'tmphas_n')

    allocate(tmphas_s(0:nphas-1,1:1),stat = status)
    call assert_alloc(status,code,'tmphas_s')

    allocate(keep_it(0:nrings-1),stat=status) 
    call assert_alloc(status,code,'keep_it')

    allocate( kphi0(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'kphi0')

    allocate( totWght(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'totWght')

    allocate( keep_north(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'keep_north')

    allocate( keep_south(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'keep_south')

    allocate( nph(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'nph')

    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')
	
    allocate(phring(0:nmmax),stat = status)	   
    call assert_alloc(status,code,'phring')		
	
    allocate( indx( 0:nmmax), stat = status)
    call assert_alloc(status,code,'indx')

    call fix_mmap( 1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

    !     ------------ initiate variables and arrays ----------------

    ! loop over Stokes params	
 
   do i=1,3
	
       istart_north = 0
       istart_south = map_size
       tmphas_n(:,:) = 0_dpc
       tmphas_s(:,:) = 0_dpc	   

       do ichunk = 0, nchunks-1

          lchk = first_ring + ichunk * chunksize 
          uchk = min(lchk+chunksize - 1, last_ring)

          do ith = lchk, uchk
             ithl = ith - lchk !local index

             ! extract pixel location information
             nph(ithl) = pixelization%nph(ith)
             cth0 = pixelization%cth(ith)
             sth0 = pixelization%sth(ith)
             kphi0(ithl) = pixelization%kphi(ith)
             totWght(ithl) = pixelization%parea(ith)*pixelization%qwght(ith)   ! combine all the ring weights to be applied

             ! find out which rings are to be analysed
             keep_north( ithl) = scan%nfl( ith)
             keep_south( ithl) = scan%sfl( ith)
             keep_it( ith-first_ring) = scan%fl( ith)

          enddo

          !-----------------------------------------------------------------------
          !  computes the integral in phi : phas_m(theta)
          !  for each parallel from Pole to Equator (use symmetry for S. Hemisphere)
          !-----------------------------------------------------------------------

          ! do Fourier Transform on rings
          do ith = lchk, uchk
             ithl = ith - lchk !local index
	     iring = ith - first_ring

             nphl = nph(ithl)

             if (keep_north(ithl) == 1) then
                ring(0:nphl-1) = local_map(istart_north:istart_north+nphl-1,i,imap) * &
                               & local_w8ring(iring+1,i)*totWght( ithl)
                call s2hat_ring_analysis(nlmax, nmmax, ring, nphl, phring, kphi0(ithl))
			  
                do j = 0, nmmax
                   ij = iring+indx(j)*nrings
                   tmphas_n(ij,1) = phring(j)
	        enddo

                istart_north = istart_north+nphl
             endif

             if (keep_south(ithl) == 1) then
                istart_south = istart_south-nphl

                ring(0:nphl-1) = local_map(istart_south:istart_south+nphl-1,i,imap) * &
                               & local_w8ring(iring+1,i)*totWght( ithl)
                call s2hat_ring_analysis(nlmax, nmmax, ring, nphl, phring, kphi0(ithl))

                do j = 0, nmmax
                   ij = iring+indx(j)*nrings 
                   tmphas_s(ij,1) = phring(j)
	        enddo
             endif
          enddo ! over rings of chunks
       enddo ! over chunks

       ! redistribute intermediate products over procs

       tmphas_n_ptr => tmphas_n
       tmphas_s_ptr => tmphas_s

       call redistInterProductsMap2Alm( nrings, nmvals, 1, nphas, tmphas_n_ptr, myid, numprocs, comm)
       call redistInterProductsMap2Alm( nrings, nmvals, 1, nphas, tmphas_s_ptr, myid, numprocs, comm)

       do j = 0, nphas-1
	  phas_ns(2*(i-1), j) = tmphas_n(j,1)   ! north cap
	  phas_ns(2*i-1, j) = tmphas_s(j,1)     ! south cap
       enddo

    enddo ! loop over Stokes 

    deallocate( tmphas_n, tmphas_s)
    deallocate( keep_north, keep_south)
    deallocate( kphi0, totWght)
    deallocate( nph)
    deallocate( ring)
    deallocate( phring)
    deallocate( indx)

    ! and get global info about all the rings
  
    allocate( keep_it_all(0:nringsall-1))
    keep_it_all(:) = 0
    call gatherRingInfo( first_ring, nrings, keep_it, nringsall, keep_it_all, myid, numprocs, comm)
    deallocate( keep_it)

    !-----------------------------------------------------------------------
    !              computes the a_lm by integrating over theta
    !                  lambda_lm(theta) * phas_m(theta)
    !                         for each m and l
    !-----------------------------------------------------------------------

    allocate(normal_l(0:nlmax),stat = status)
    call assert_alloc(status,code,'normal_l')

    ! generate Polarization normalisation, etc
    call gen_normpol(nlmax, normal_l)

    allocate(dalm(0:5,0:nlmax), stat = status)
    call assert_alloc(status,code,'dalm')

    allocate(lam_fact(0:nlmax),stat = status)
    call assert_alloc(status,code,'lam_fact')

    nchunks   = nringsall/SMAXCHK + 1  ! number of chunks
    chunksize = (nringsall+nchunks-1)/nchunks ! <= SMAXCHK

    local_alms(:,:,:,imap) = 0.0            ! set the whole alm array to zero

    kring = 0
    do ichunk = 0, nchunks-1

       lchk  = ichunk * chunksize + 1
       uchk  = min(lchk+chunksize - 1, nringsall)
       
       mshift = 0

       do mindx = 0, nmvals-1
	  
          mindx1 = mindx+1 
	  m = mvals( mindx)

          ! generate Ylm relation factor for degree m
          call gen_lamfac(nlmax, m, lam_fact)
          fm2 = real(m * m, kind = DP)
          normal_m = (2.0_dp * m) * (1 - m)

          ! introduce double precision vector to perform summation over ith for each l
          dalm(0:5, m:nlmax) = 0.0_dp

          jring = kring
          do ithl = 0, uchk-lchk

             iring = lchk+ithl-1
			 
  	     cth_ring = pixelization%cth(iring+1)
             one_on_s2 =  1.0_dp / pixelization%sth(iring+1)**2
             c_on_s2 = cth_ring*one_on_s2
	     two_on_s2 = 2.0_dp * one_on_s2
	     two_cth_ring = 2.0_dp * cth_ring
	     b_w  =  c_on_s2   			 

             l_min = l_min_ylm(m, pixelization%sth(iring+1)) ! lower l bound of non-negligible Ylm

             if (keep_it_all(iring) == 1) then

                if (nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001) split in two (RS, 02-2007)

                   i_mm = n_lm * jring + mshift ! location of Ym,m for ring ith
                   i_mm1 = i_mm+1

                   do i=0,2 ! loop on T, Q, U
                      phm = phas_ns(2*i, mindx*nringsobs+jring) - phas_ns(2*i+1, mindx*nringsobs+jring) ! N-S: (l+m) odd
                      phas_sd(-3+2*i)   =  real(phm, kind=dp) 
                      phas_sd(-3+2*i+1) = aimag(phm)
                      php = phas_ns(2*i, mindx*nringsobs+jring) + phas_ns(2*i+1, mindx*nringsobs+jring) ! N+S: (l+m) even
                      phas_sd( 3+2*i)   =  real(php, kind=dp) 
                      phas_sd( 3+2*i+1) = aimag(php)
                   enddo
                   phas_sdx(-3: 0) = (/ phas_sd(8), - phas_sd(7), - phas_sd(6),   phas_sd(5) /)
                   phas_sdx( 3: 6) = (/ phas_sd(2), - phas_sd(1), - phas_sd(0),   phas_sd(-1) /)

                   !           ---------- l = m ----------
                   par_lm = 3  ! = 3 * (-1)^(l+m)
                   lam_lm = local_plm_comp_a(i_mm1)

                   if (m >= l_min) then
                      lambda_w =  (normal_m * lam_lm)  * ( 0.5_dp - one_on_s2)
                      lambda_x =  (normal_m * lam_lm)  * c_on_s2

                      ! alm_T = \int T Ylm
                      ! alm_G = \int ( - Q Wlm - i U Xlm )
                      ! alm_C = \int ( - U Wlm + i Q Xlm )
                      dalm(0:1, m) = dalm(0:1, m) + lam_lm * phas_sd(par_lm:par_lm+1)
                      dalm(2:5, m) = dalm(2:5, m)  &
                           &       - lambda_w * phas_sd(par_lm+2:par_lm+5) &
                           &       + lambda_x * phas_sdx(par_lm:par_lm+3)
                   endif

                   !           ---------- l > m ----------

                   fm_on_s2 = m * one_on_s2

                   l_start = max(m+1, l_min)
                   if (mod(m+l_start,2) == 0) l_start = l_start+1
                   do l = l_start, nlmax, 2

                      lam_lm1m = local_plm_comp_a(i_mm1+l-m-1) * lam_fact(l)
                      lam_lm   = local_plm_comp_a(i_mm1+l-m)

                      fl = real(l, kind = DP)
                      flm1 = fl - 1.0_dp
                      a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                      a_x =  two_cth_ring * flm1
                      lambda_w = b_w * lam_lm1m - a_w * lam_lm 
                      lambda_x = fm_on_s2 * ( lam_lm1m - a_x * lam_lm)
                   
                      dalm(0,l) = dalm(0,l) + lam_lm * phas_sd(-3)
                      dalm(1,l) = dalm(1,l) + lam_lm * phas_sd(-2)

                      dalm(2,l) = dalm(2,l) - lambda_w * phas_sd(-1) + lambda_x * phas_sdx(-3)
                      dalm(3,l) = dalm(3,l) - lambda_w * phas_sd(0) + lambda_x * phas_sdx(-2)
                      dalm(4,l) = dalm(4,l) - lambda_w * phas_sd(1) + lambda_x * phas_sdx(-1)
                      dalm(5,l) = dalm(5,l) - lambda_w * phas_sd(2) + lambda_x * phas_sdx(0)
                   enddo ! loop over l

                   l_start = max(m+2, l_min)
                   if (mod(m+l_start,2) == 1) l_start = l_start+1
                   do l = l_start, nlmax, 2

                      lam_lm1m = local_plm_comp_a(i_mm1+l-m-1) * lam_fact(l)
                      lam_lm   = local_plm_comp_a(i_mm1+l-m)

                      fl = real(l, kind = DP)
                      flm1 = fl - 1.0_dp
                      a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                      a_x =  two_cth_ring * flm1
                      lambda_w = b_w * lam_lm1m - a_w * lam_lm
                      lambda_x = fm_on_s2 * (lam_lm1m - a_x * lam_lm)

                      dalm(0,l) = dalm(0,l) + lam_lm * phas_sd(3)
                      dalm(1,l) = dalm(1,l) + lam_lm * phas_sd(4)
                      dalm(2,l) = dalm(2,l) - lambda_w * phas_sd(5) + lambda_x * phas_sdx(3)
                      dalm(3,l) = dalm(3,l) - lambda_w * phas_sd(6) + lambda_x * phas_sdx(4)
                      dalm(4,l) = dalm(4,l) - lambda_w * phas_sd(7) + lambda_x * phas_sdx(5)
                      dalm(5,l) = dalm(5,l) - lambda_w * phas_sd(8) + lambda_x * phas_sdx(6)
                   enddo ! loop over l

                endif ! test for l values

                jring = jring + 1

             endif ! test on cut sky

          enddo ! loop over chunk rings (ithl)

          do l = m, nlmax
             l1 = l+1

             local_alms_comp_a(l1,mindx1) = local_alms_comp_a(l1,mindx1) + &
                                          & cmplx(dalm(0,l), dalm(1,l))
             local_alms_comp_b(l1,mindx1) = local_alms_comp_b(l1,mindx1) + &
                                          & cmplx(dalm(2,l), dalm(3,l)) * normal_l(l)
             local_alms_comp_c(l1,mindx1) = local_alms_comp_c(l1,mindx1) + &
                                          & cmplx(dalm(4,l), dalm(5,l)) * normal_l(l)
          enddo

	  mshift = mshift +nlmax-m+1

       enddo ! loop over m

       kring = jring ! store the number of obs rings so far

    enddo    ! loop over chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate (dalm, lam_fact)
    deallocate(normal_l, phas_ns)
    deallocate( keep_it_all)

    return

  end subroutine s2hat_map2alm_pol_pre1


  subroutine s2hat_map2alm_pol_pre2( pixelization, scan, nlmax, nmmax, nmvals, mvals, imap, first_ring, last_ring, local_w8ring, map_size, local_map, &
                                   & lda, nplm, local_plm, local_alms, numprocs, myid, comm)

    !====================================================================================
    !     computes the a(l,m) from a polarized map for the HEALPIX pixelisation
    !     using precomputed and stored scalar and spin-2 Ylm_T(theta) and Ylm_P(theta)
    !====================================================================================

    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan
   
    integer (i4b) :: nlmax, nmmax, nmvals, imap, first_ring, last_ring, lda, map_size, numprocs, myid, comm
    integer(i4b), dimension(0:nmvals-1) :: mvals
    real(dp), dimension(:,:), pointer :: local_w8ring
    real(dp), dimension(:,:,:), pointer :: local_map
    integer(i8b) :: nplm
    real(dp), dimension(:,:), pointer :: local_plm	

    ! output
    complex (dp), dimension(:,:,:,:), pointer :: local_alms

    ! internal
    integer(I4B) :: l, l1, m, ith, ithl, scalel, scalem, nrings, nphmx, mindx, mindx1
    integer(I8B) :: istart_south, istart_north, npix, i_mm, i_mm1, mshift, i_up, i_up1

    integer(I4B)   :: par_lm, ichunk, lchk, uchk, chunksize, nchunks
    real(DP)       :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)       :: lambda_w, lambda_x, normal_m, lam_lm1m
    real(DP)       :: fm2, fl, flm1, cth0, sth0, one_on_s2, c_on_s2
    real(DP)       :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP),     dimension(-3:8)             :: phas_sd
    real(DP),     dimension(-3:6)             :: phas_sdx
    real(DP),     dimension(:,:), allocatable :: dalm, plm_sub
    real(DP),     dimension(:),   allocatable :: lam_fact

    integer(i4b)                              :: l_min
    real(DP),     dimension(:),   allocatable :: normal_l
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:), pointer :: tmphas_n_ptr, tmphas_s_ptr
    complex(DPC), dimension(:,:), allocatable, target :: tmphas_n, tmphas_s
    complex(DPC), dimension(:), allocatable :: phring
    complex(DPC), dimension(:,:), allocatable :: phas_ns
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)    :: nphl, i, ii, ij, j, iring, jring, kring, nringsall, nringsobs, nphas
    integer(i4b)    :: l_start
    integer(I8B)    :: startpix, n_lm
    real(DP),     dimension(:), allocatable :: cth, sth, kphi0, totWght
    integer(i4b), dimension(:), allocatable :: keep_it, keep_north, keep_south, keep_it_all
    integer(i4b), dimension(:), allocatable :: indx, nph

    character(LEN=*), parameter :: code = 'MAP2ALM'
    integer(I4B) :: status

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

    ! Healpix definitions
    nrings = (last_ring-first_ring+1)

    nringsall = pixelization%nringsall
    nringsobs = scan%nringsobs
    npix   = pixelization%npixsall
    nphmx  = pixelization%nphmx

    n_lm = nummmodes( nlmax, nmvals, mvals)

    !     --- allocates space for arrays ---
	
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    nphas = max( (nmmax+1)*nrings, nmvals*nringsobs)

    allocate(phas_ns(0:5,0:nphas-1),stat = status)
    call assert_alloc(status,code,'phas_ns')

    allocate(tmphas_n(0:nphas-1,1:1),stat = status)
    call assert_alloc(status,code,'tmphas_n')

    allocate(tmphas_s(0:nphas-1,1:1),stat = status)
    call assert_alloc(status,code,'tmphas_s')

    allocate(keep_it(0:nrings-1),stat=status) 
    call assert_alloc(status,code,'keep_it')

    allocate( kphi0(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'kphi0')

    allocate( totWght(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'totWght')

    allocate( keep_north(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'keep_north')

    allocate( keep_south(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'keep_south')

    allocate( nph(0:SMAXCHK-1), stat = status)
    call assert_alloc(status,code,'nph')

    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')
	
    allocate(phring(0:nmmax),stat = status)	   
    call assert_alloc(status,code,'phring')		
	
    allocate( indx( 0:nmmax), stat = status)
    call assert_alloc(status,code,'indx')

    call fix_mmap( 1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

    !     ------------ initiate variables and arrays ----------------

    ! loop over Stokes params	
 
   do i=1,3
	
       istart_north = 0
       istart_south = map_size
       tmphas_n(:,:) = 0_dpc
       tmphas_s(:,:) = 0_dpc	   

       do ichunk = 0, nchunks-1

          lchk = first_ring + ichunk * chunksize 
          uchk = min(lchk+chunksize - 1, last_ring)

          do ith = lchk, uchk
             ithl = ith - lchk !local index

            ! extract pixel location information
             nph(ithl) = pixelization%nph(ith)
             cth0 = pixelization%cth(ith)
             sth0 = pixelization%sth(ith)
             kphi0(ithl) = pixelization%kphi(ith)
             totWght(ithl) = pixelization%parea(ith)*pixelization%qwght(ith)

             ! find out which rings are to be analysed
             keep_north( ithl) = scan%nfl( ith)
             keep_south( ithl) = scan%sfl( ith)
             keep_it( ith-first_ring) = scan%fl( ith)

          enddo

          !-----------------------------------------------------------------------
          !  computes the integral in phi : phas_m(theta)
          !  for each parallel from Pole to Equator (use symmetry for S. Hemisphere)
          !-----------------------------------------------------------------------

          ! do Fourier Transform on rings
          do ith = lchk, uchk
             ithl = ith - lchk !local index
	     iring = ith - first_ring

             nphl = nph(ithl)

             if (keep_north(ithl) == 1) then
                ring(0:nphl-1) = local_map(istart_north:istart_north+nphl-1,i,imap) * &
                               & local_w8ring(iring+1,i)*totWght( ithl)
                call s2hat_ring_analysis(nlmax, nmmax, ring, nphl, phring, kphi0(ithl))
			  
                do j = 0, nmmax
                   ij = iring+indx(j)*nrings
                   tmphas_n(ij, 1) = phring(j)
	        enddo

                istart_north = istart_north+nphl
             endif

             if (keep_south(ithl) == 1) then
                istart_south = istart_south-nphl

                ring(0:nphl-1) = local_map(istart_south:istart_south+nphl-1,i,imap) * &
                               & local_w8ring(iring+1,i)*totWght( ithl)
                call s2hat_ring_analysis(nlmax, nmmax, ring, nphl, phring, kphi0(ithl))

                do j = 0, nmmax
                   ij = iring+indx(j)*nrings 
                   tmphas_s(ij, 1) = phring(j)
	        enddo
             endif
          enddo ! over rings of chunks
       enddo ! over chunks

       ! redistribute intermediate products over procs

       tmphas_n_ptr => tmphas_n
       tmphas_s_ptr => tmphas_s

       call redistInterProductsMap2Alm( nrings, nmvals, 1, nphas, tmphas_n_ptr, myid, numprocs, comm)
       call redistInterProductsMap2Alm( nrings, nmvals, 1, nphas, tmphas_s_ptr, myid, numprocs, comm)

       do j = 0, nphas-1
	  phas_ns(2*(i-1), j) = tmphas_n(j,1)   ! north cap
	  phas_ns(2*i-1, j) = tmphas_s(j,1)     ! south cap
       enddo

    enddo ! loop over Stokes 

    deallocate( tmphas_n, tmphas_s)
    deallocate( keep_north, keep_south)
    deallocate( kphi0, totWght)
    deallocate( nph)
    deallocate( ring)
    deallocate( phring)
    deallocate( indx)

    ! and get global info about all the rings
  
    allocate( keep_it_all(0:nringsall-1))
    keep_it_all(:) = 0
    call gatherRingInfo( first_ring, nrings, keep_it, nringsall, keep_it_all, myid, numprocs, comm)
    deallocate( keep_it)

    !-----------------------------------------------------------------------
    !              computes the a_lm by integrating over theta
    !                  lambda_lm(theta) * phas_m(theta)
    !                         for each m and l
    !-----------------------------------------------------------------------

    allocate(dalm(0:5,0:nlmax), stat = status)
    call assert_alloc(status,code,'dalm')

    allocate( plm_sub(1:3,0:nlmax), stat = status)
    call assert_alloc(status,code,'plm_sub')

    nchunks   = nringsall/SMAXCHK + 1  ! number of chunks
    chunksize = (nringsall+nchunks-1)/nchunks ! <= SMAXCHK

    local_alms(:,:,:,imap) = 0.0            ! set the whole alm array to zero

    kring = 0    
    do ichunk = 0, nchunks-1

       lchk  = ichunk * chunksize + 1
       uchk  = min(lchk+chunksize - 1, nringsall)
       
       mshift = 0

       do mindx = 0, nmvals-1
	  
          mindx1 = mindx+1 
	  m = mvals( mindx)
		  
          ! introduce double precision vector to perform summation over ith for each l
          dalm(0:5, m:nlmax) = 0.0_dp

          jring = kring
          do ithl = 0, uchk - lchk

             iring = lchk+ithl-1

             l_min = l_min_ylm(m, pixelization%sth(iring+1)) ! lower l bound of non-negligible Ylm

              if (keep_it_all(iring) == 1) then

                if (nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001) split in two (RS, 02-2007)

                   i_mm = n_lm * jring + mshift ! location of Ym,m for ring ith
                   i_mm1 = i_mm+1
                
                   i_up = i_mm + nlmax - m
                   i_up1 = i_up+1

                   plm_sub(1,m:nlmax) = local_plm_comp_a(i_mm1:i_up1)
                   plm_sub(2,m:nlmax) = local_plm_comp_b(i_mm1:i_up1)
                   plm_sub(3,m:nlmax) = local_plm_comp_c(i_mm1:i_up1)

                   do i=0,2 ! loop over T, Q, U
                      phm = phas_ns(2*i, mindx*nringsobs+jring) - phas_ns(2*i+1, mindx*nringsobs+jring) ! N-S: (l+m) odd
                      phas_sd(-3+2*i)   =  real(phm, kind=dp) 
                      phas_sd(-3+2*i+1) = aimag(phm)
                      php = phas_ns(2*i, mindx*nringsobs+jring) + phas_ns(2*i+1, mindx*nringsobs+jring) ! N+S: (l+m) even
                      phas_sd( 3+2*i)   =  real(php, kind=dp) 
                      phas_sd( 3+2*i+1) = aimag(php)
                   enddo
                   phas_sdx(-3: 0) = (/ phas_sd(8), - phas_sd(7), - phas_sd(6),   phas_sd(5) /)
                   phas_sdx( 3: 6) = (/ phas_sd(2), - phas_sd(1), - phas_sd(0),   phas_sd(-1) /)

                   !           ---------- l = m ----------
                   par_lm = 3  ! = 3 * (-1)^(l+m)

                   if (m >= l_min) then
                   
                      ! alm_T = \int T Ylm
                      ! alm_G = \int ( - Q Wlm - i U Xlm )
                      ! alm_C = \int ( - U Wlm + i Q Xlm )
                      dalm(0:1, m) = dalm(0:1, m) + plm_sub(1,m) * phas_sd(par_lm:par_lm+1)
                      dalm(2:5, m) = dalm(2:5, m)  &
                           &       - plm_sub(2,m) * phas_sd(par_lm+2:par_lm+5) &
                           &       + plm_sub(3,m) * phas_sdx(par_lm:par_lm+3)
                   endif

                   !           ---------- l > m ----------

                   do l = m+1, nlmax
                      par_lm = - par_lm  ! = 3 * (-1)^(l+m)
                      if (l >= l_min) then
                         dalm(0:1, l) = dalm(0:1, l) + plm_sub(1,l) * phas_sd(par_lm:par_lm+1)
                         dalm(2:5, l) = dalm(2:5, l)  &
                              &       - plm_sub(2,l) * phas_sd(par_lm+2:par_lm+5) &
                              &       + plm_sub(3,l) * phas_sdx(par_lm:par_lm+3)
                      endif
                   
                   enddo ! loop over l

                endif ! test for l values

                jring = jring + 1

             endif ! test on cut sky

          enddo ! loop over chunk rings (ithl)

          do l = m, nlmax
             l1 = l+1

             local_alms_comp_a(l1,mindx1) = local_alms_comp_a(l1,mindx1) + cmplx(dalm(0,l), dalm(1,l))
             local_alms_comp_b(l1,mindx1) = local_alms_comp_b(l1,mindx1) + cmplx(dalm(2,l), dalm(3,l))
             local_alms_comp_c(l1,mindx1) = local_alms_comp_c(l1,mindx1) + cmplx(dalm(4,l), dalm(5,l))
          enddo

          mshift = mshift +nlmax-m+1

       enddo ! loop over m

       kring = jring ! store the number of obs rings so far

    enddo    ! loop over chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate(phas_ns)
    deallocate (dalm, plm_sub)
    deallocate( keep_it_all)

    return		  
		  
 end subroutine s2hat_map2alm_pol_pre2
  
end module s2hat_map2alm_mod
