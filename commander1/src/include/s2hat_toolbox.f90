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
! January, 08, 2007, Radek Stompor (APC)
!
! The following routines enable computation of the direct and inverse
! spherical harmonic transforms of polarized (nearly) full sky maps.
! They were developed on the basis of the HEALPix software (see legal matters
! above) and make use of the MPI library to perform message passing
! communication.
! The routines provided here permit distributing both pixel and harmonic 
! domain object over all processors. The spherical harmonic transform is 
! splitted into two steps separated by the data redistribution. The latter 
! involves a global mpi communication (mpi_alltoallv) (routines redistInterProducts*).
!
! This file contains outines performing the load balanced data distribution 
! and other auxiliary procedures. The complementary files: s2hat_map2alm.f90
! and s2hat_alm2map.f90 contains the map-to-alm and alm-to-map transforms.
!
!--------------------------------------------------------------------------
!
! February, 2007 RS@APC
!
! - load balancing improved for partial sky experiments
! - generalized to aribtrary isolatitudal pixelization scheme with pixels 
!   equidistant in azimuth for each isolatidude ring, and symmetric wrt to
!   the equator.
!
!--------------------------------------------------------------------------
!
! July/Aug, 2007 RS@APC
!
! - arbitrary spin functions added.
!
!--------------------------------------------------------------------------
!
! October 2007 - rs@apc
!
! ring_analysis and ring_synthesis adapted from Healpix to allow for the
! more general pixelization schemes.
!
!--------------------------------------------------------------------------

module s2hat_toolbox_mod

  use healpix_types
  use s2hat_types
  ! use s2hat_pixelization
  use healpix_fft
  use misc_utils
  use alm_tools
  
  implicit none

  ! include 'mpif.h'

  private
  
  ! define large and small numbers used to renormalise the recursion 
  ! on the Legendre Polynomials
  integer(I4B),      public, parameter :: LOG2LG   = 100
  real(KIND=DP),     public, parameter :: FL_LARGE = 2.0_dp **   LOG2LG
  real(KIND=DP),     public, parameter :: FL_SMALL = 2.0_dp ** (-LOG2LG)
  ! declare array for dynamic rescaling of the Ylm
  integer(kind=i4b), public, parameter :: RSMAX = 20, RSMIN = -20
  real(dp),          public, dimension(RSMIN:RSMAX) :: rescale_tab
  real(DP),          public, parameter :: ALN2_INV = 1.4426950408889634073599246810_dp ! 1/log(2)
  integer(i4b),      public, parameter :: SMAXCHK = 50 ! maximum size of chunk (in ring pairs)
  
  !--------------------------------------------------------------------------

  public :: get_local_data_sizes
  public :: distribute_alms, distribute_w8ring, distribute_map, distribute_local_data_objects_map2alm, distribute_local_data_objects_alm2map
  public :: plm_mvalues_gen
  public :: collect_alms, collect_cls, collect_xls, collect_map
  public :: local_map_size, gatherRingInfo
  public :: find_mvalues, find_ring_range, fix_mmap, find_pixel_numbers, nummvalues, nummmodes
  public :: s2hat_ring_analysis, s2hat_ring_synthesis
  public :: mpi_do_lam_lm_pol, mpi_compute_lam_mm, init_rescale_local
  public :: mpi_do_slam_lm, mpi_compute_slam_mm, s2hat_gen_recfac, s2hat_gen_mfac, s2hat_gen_sfac

contains
  
  ! *****************************************************************
  ! *           Distribution routines
  ! *****************************************************************
  
  subroutine get_local_data_sizes( precompute_plms, pixelization, scan, nlmax, nmmax, myid, numprocs, nmvals, &
                                 & first_ring, last_ring, map_size, nplm, root, comm)

    ! computes the local object sizes needed for their allocation

    implicit none
    
    ! input
    type( pixeltype), intent(in) :: pixelization
    type( scandef), intent(in)  :: scan
    integer(i4b), intent(in)  :: precompute_plms, nlmax, nmmax, myid, numprocs, root, comm
    ! output
    integer(i4b), intent(out) :: first_ring, last_ring, map_size, nmvals
    integer(i8b), intent(out) :: nplm
	
    ! internal
    integer(i4b) ::  errcode, mpierr, outerror

    ! Find which rings to compute
    call find_ring_range( pixelization, scan, nmmax, myid, numprocs, first_ring, last_ring, outerror)

    if( outerror == 1) then
	  if(myid == root) write(*,*)'ring distribution failed -- too many procs ?!'
	  call mpi_abort( mpi_comm_world, errcode, mpierr )
    endif

    ! Find local map size
    map_size = local_map_size( pixelization, first_ring, last_ring)

    ! Find the m interval for each proc
    nmvals =  nummvalues( myid, numprocs, nmmax) 

    ! Estimate size of the plm array
    if( precompute_plms == 0) then 
	  nplm = 1
	else 
	  nplm = (nmvals+1)/2*(nlmax+2)*scan%nringsobs               ! that may overshoot for even nmmax values	
    endif

  end subroutine get_local_data_sizes  

  subroutine distribute_local_data_objects_map2alm( precompute_plms, pixelization, scan, nlmax, nmmax, nmaps, nstokes, map, map_size, local_map, &
                                                  & nmvals, mvals, lda, nplm, local_plm, first_ring, last_ring, local_w8ring, w8ring, myid, numprocs, root, comm)

    ! distributes and/or precomputes the data on all procs - a simple global "do-it-all" driver - (modified HEALPix routine)

    ! input
    type( pixeltype), intent( in) :: pixelization
    type( scandef), intent( in) :: scan
    integer(i4b), intent( in) :: precompute_plms, nmaps, nstokes, map_size, nlmax, nmmax, nmvals, first_ring, last_ring, lda
    integer(i4b), intent( in) :: myid, numprocs, root, comm
    integer(i8b), intent( in) :: nplm
    real(dp), dimension(:,:,:), pointer :: map
    real(dp), dimension(:,:), pointer :: w8ring

    ! output
    real(dp), dimension(:,:,:), pointer :: local_map
    integer(i4b), dimension(0:nmvals-1), intent( out) :: mvals
    real(dp), dimension(:,:), pointer :: local_w8ring
    real(dp), dimension(:,:), pointer :: local_plm

    ! internal

    integer( i4b) :: npixsall, nringsall, imap
    real(dp), dimension(:,:), allocatable, target :: tmpmap
    real(dp), dimension(:,:), pointer :: tmpmap_ptr

    nringsall = pixelization%nringsall
    npixsall = pixelization%npixsall

    ! Distribute ring weights
    call distribute_w8ring( nstokes, first_ring, last_ring, local_w8ring, w8ring, myid, numprocs, root, comm)

    ! and map
    allocate( tmpmap(0:npixsall-1,1:nstokes))
    do imap = 1, nmaps
       tmpmap(:,:) = map(:,:,imap)
       tmpmap_ptr => tmpmap
       call distribute_map( pixelization, nmaps, imap, nstokes, first_ring, last_ring, map_size, local_map, tmpmap_ptr, myid, numprocs, root, comm)	
    enddo
    deallocate( tmpmap)

    ! define the mvalues assignment
    call find_mvalues(myid, numprocs, nmmax, nmvals, mvals)

    ! and pre-compute plms if requested
    if (precompute_plms == 1) then
       call plm_mvalues_gen( pixelization, scan, 1, nlmax, nmmax, nmvals, mvals, lda, nplm, local_plm)
    else if (precompute_plms == 2) then
       call plm_mvalues_gen( pixelization, scan, 3, nlmax, nmmax, nmvals, mvals, lda, nplm, local_plm)
    end if

  end subroutine distribute_local_data_objects_map2alm
  
  subroutine distribute_local_data_objects_alm2map( precompute_plms, pixelization, scan, nlmax, nmmax, nmaps, nstokes, lda, alms, nmvals, mvals, local_alm, &
                                                  &  nplm, local_plm, myid, numprocs, root, comm)

    ! distributes and/or precomputes the data on all procs - a simple global "do-it-all" driver - (modified HEALPix routine)

    ! input
    type( pixeltype), intent( in) :: pixelization
    type( scandef), intent( in) :: scan
    integer(i4b), intent( in) :: precompute_plms, nmaps, nstokes, nlmax, nmmax, nmvals, lda
    integer(i4b), intent( in) :: myid, numprocs, root, comm
    integer(i8b), intent( in) :: nplm
    complex(dp), dimension(:,:,:,:), pointer :: alms

    ! output
    complex(dp), dimension(:,:,:,:), pointer :: local_alm
    integer(i4b), dimension(0:nmvals-1), intent( out) :: mvals
    real(dp), dimension(:,:), pointer :: local_plm

    ! internal
    integer(i4b) :: imap
    complex(dp), dimension(:,:,:), allocatable, target :: tmpalms
    complex(dp), dimension(:,:,:), pointer :: tmpalms_ptr

    ! define the mvalues assignment
    call find_mvalues(myid, numprocs, nmmax, nmvals, mvals)

    ! distribute alms
    allocate( tmpalms(0:nlmax,0:nmmax,1:nstokes))
    do imap = 1, nmaps
        tmpalms = alms(:,:,:,imap)
        tmpalms_ptr => tmpalms
	call distribute_alms( nlmax, nmmax, nmaps, imap, nstokes, nmvals, mvals, lda, local_alm, tmpalms_ptr, myid, numprocs, root, comm)
    enddo 
    deallocate( tmpalms)

    ! and pre-compute plms if requested
    if (precompute_plms == 1) then
       call plm_mvalues_gen( pixelization, scan, 1, nlmax, nmmax, nmvals, mvals, lda, nplm, local_plm)
    else if (precompute_plms == 2) then
       call plm_mvalues_gen( pixelization, scan, 3, nlmax, nmmax, nmvals, mvals, lda, nplm, local_plm)
    end if

  end subroutine distribute_local_data_objects_alm2map
  
  subroutine collect_map( pixelization, nmaps, mapnum, nstokes, map, first_ring, last_ring, map_size, local_map, myid, numprocs, root, comm)

    ! collect the map data which are distributed over all procs
    ! Note that if the map covers only part of the sky, the unobserved reminder will be left
    ! as on the input.

    implicit none
    
    ! input
    type( pixeltype) :: pixelization
    integer(i4b), intent(in) :: nmaps, mapnum, nstokes, numprocs, myid, root, comm
    integer(i4b), intent(in) :: first_ring, last_ring, map_size
    real(DP), dimension(:,:,:), pointer :: local_map
    ! output
    real(DP), dimension(:,:), pointer  :: map

    ! internal
    integer(i4b) :: i, istokes, j, seg_size, ierr, nsegment, send_count0, send_count1
    integer(i4b), dimension(MPI_STATUS_SIZE) :: status 
    integer(i4b), dimension(0:1,0:1) :: segment
    integer(i4b), dimension(:), allocatable :: num_segments, tmp, recv_count, recv_disps
    integer(i4b), dimension(:,:,:), allocatable :: segments

    ! get pixel ranges
    call find_pixel_numbers( pixelization, first_ring, last_ring, nsegment, segment)

    if( myid == root) then
       allocate( num_segments(0:numprocs-1))
       allocate( tmp(0:numprocs-1))
       allocate( segments(0:1,0:1,0:numprocs-1))
    end if

    call mpi_gather(nsegment, 1, mpi_integer, num_segments, 1, mpi_integer, root, comm, ierr)

    do i = 0, 1
       do j = 0, 1
          call mpi_gather(segment(i,j), 1, mpi_integer, tmp, 1, mpi_integer, root, comm, ierr)
          if( myid == root) segments(i,j,:) = tmp(:)
       end do
    end do

    if( myid == root) deallocate( tmp)

    ! send first set of rows now ...

    if( myid == root) then

       allocate( recv_count(0:numprocs-1))
       allocate( recv_disps(0:numprocs-1))

       do i = 0, numprocs-1
          recv_count(i) = segments(0,1,i)-segments(0,0,i)+1
          recv_disps(i) = segments(0,0,i)
       enddo

    endif

    send_count0 = segment(0,1)-segment(0,0)+1

    do istokes = 1, nstokes
       call mpi_gatherv( local_map(0,istokes,mapnum), send_count0, MPI_DOUBLE_PRECISION, map(0, istokes), recv_count, recv_disps, MPI_DOUBLE_PRECISION, root, comm, ierr)
    enddo

    ! and the second row now ...

    if( myid == root) then

       do i = 0, numprocs-1

          if( num_segments(i) == 1) then
             recv_count(i) = 0
             if( i /= 0) then 
                recv_disps(i) = recv_disps(i-1)+recv_count(i-1) 
             else 
                recv_disps(0) = 0
             endif
          else  
             recv_count(i) = segments(1,1,i)-segments(1,0,i)+1
             recv_disps(i) = segments(1,0,i)
          endif
       enddo

    endif

    send_count1 = 0
    if( nsegment == 2) send_count1 = segment(1,1)-segment(1,0)+1

    do istokes = 1, nstokes
       call mpi_gatherv( local_map( send_count0, istokes, mapnum), send_count1, MPI_DOUBLE_PRECISION, map(0, istokes), recv_count, recv_disps, MPI_DOUBLE_PRECISION, root, comm, ierr)
    enddo

    if( myid == root) deallocate( num_segments, segments, recv_count, recv_disps)

  end subroutine collect_map
  
  subroutine collect_alms( nlmax, nmmax, nmaps, mapnum, nstokes, nmvals, mvals, lda, local_alm, alms, myid, numprocs, root, comm)

    ! collects all the alms distributed over proc on a proc root

    implicit none

    !input
    integer(i4b), intent(in) :: nlmax, nmmax, nmaps, mapnum, nstokes, nmvals, lda, myid, numprocs, root, comm
    integer(i4b), dimension(0:nmvals-1), intent(in) :: mvals
    complex(DP),  dimension(:,:,:,:), pointer :: local_alm

    !output
    complex(DP),  dimension(:,:,:), pointer :: alms

    !internal
    integer(I4B) :: i, ip, im, ierr, status, send_count
    integer(I4B), dimension(:), allocatable :: mperproc, recv_count, recv_disps, indx

    if( myid == root) then
       allocate( mperproc(0:numprocs-1), stat=status)
    endif
    
    call mpi_gather( nmvals, 1, mpi_integer, mperproc, 1, mpi_integer, root, comm, ierr)

    if( myid == root) then

       alms(:,:,:) = 0.d0

       allocate( recv_count(0:numprocs-1))
       allocate( recv_disps(0:numprocs-1))

       recv_count(0) = (nlmax+1)*mperproc(0)
       recv_disps(0) = 0

       do ip = 1, numprocs-1
          recv_count(ip) = (nlmax+1)*mperproc(ip)
          recv_disps(ip) = recv_disps(ip-1)+recv_count(ip-1)
       enddo

     endif

     send_count = (nlmax+1)*nmvals

     if( (lda == 1) .or. (lda == 2) .or. (lda == 3)) then

        send_count = send_count*nstokes

        if( myid == root) then
            recv_count(0:numprocs-1) = nstokes *  recv_count(0:numprocs-1)
            recv_disps(0:numprocs-1) = nstokes *  recv_disps(0:numprocs-1)
        end if

        call mpi_gatherv( local_alm( 1, 0, 0, mapnum), send_count, MPI_DOUBLE_COMPLEX, alms(1, 0, 0), recv_count, recv_disps, MPI_DOUBLE_COMPLEX, root, comm, ierr)
     else
        do im = 1, nstokes
           call mpi_gatherv( local_alm( 0, 0, im, mapnum), send_count, MPI_DOUBLE_COMPLEX, alms(0, 0, im), recv_count, recv_disps, MPI_DOUBLE_COMPLEX, root, comm, ierr)
        enddo
     endif

     if( myid == root)  deallocate( recv_count, recv_disps, mperproc)

     allocate( indx(0:nmmax))
     call fix_mmap( -1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

     if( myid == root) then
        call sortAlms( nlmax, nmmax, nstokes, indx, lda, alms)
     endif

     deallocate( indx)

     return

  end subroutine collect_alms

  subroutine collect_cls( nmaps, mapnum, nstokes, nlmax, nmvals, mvals, lda, local_alm, nspec, cl, myid, numprocs, root, comm)

  ! computes power spectra cl from the alms m-distributed over all procs and stores the output on a
  ! proc root

   !input
   integer(i4b), intent(in) :: nmaps, mapnum, nstokes, nlmax, nmvals, nspec, lda, myid, numprocs, root, comm
   integer(i4b), dimension(0:nmvals-1), intent(in) :: mvals
   complex(DP), dimension(:,:,:,:), POINTER  :: local_alm
   
   !output
   real(DP), dimension(0:nlmax,1:nspec), intent( out) :: cl

   !internal
   integer(I4B) :: i, l, l1, ierr, mstart, mstart1
   integer(i4b) :: is, js, ks, ncomp1, ncomp2
   real(DP), dimension(0:nlmax,1:nspec)  :: cltmp

   TYPE ROW
       complex(DP), dimension(:,:), POINTER :: comp
   END TYPE ROW

   type(ROW), dimension(1:nstokes) :: current_alm

   !=======================================================================

   if( (lda == 1) .or. (lda == 2) .or. (lda == 3)) then       ! i.e., healpix convention for the alm array

      do i = 1, lda
         current_alm(i)%comp => local_alm(i,0:nlmax,0:nmvals-1,mapnum)
      end do

   else                                       ! i.e., s2hat convention for the alm array

      do i = 1, nstokes
         current_alm(i)%comp => local_alm(0:nlmax,0:nmvals-1,i,mapnum)
      end do

   end if

   cltmp = 0.d0

   mstart = 0
   if( mvals(0) == 0) then
     mstart = 1
   endif

   ncomp1 = min( nstokes, nspec)

   mstart1 = mstart+1

   do l = mvals(mstart), nlmax
         
       l1 = l+1
    
       do is = 1, ncomp1
          cltmp(l,is) = 2.*SUM( current_alm(is)%comp(l1,mstart1:nmvals)*CONJG(current_alm(is)%comp(l1,mstart1:nmvals)))
       end do
   enddo

   if( (nstokes >= 2) .and. (nspec >= nstokes+1)) then                                 ! and now do cross-spectra if requested

       ncomp1 = min( nstokes, nspec-nstokes)
       ncomp2 = min( nstokes, nspec-nstokes+1)

       mstart1 = mstart+1
       do l = mvals(mstart), nlmax
          l1 = l+1

          ks = nstokes+1

          do is = 1, ncomp1                                                           ! do all x-spectra
              do js = is+1, ncomp2
                  cltmp(l,ks) = 2.*SUM( current_alm(is)%comp(l1,mstart1:nmvals)*CONJG(current_alm(js)%comp(l1,mstart1:nmvals)))
                  ks = ks+1
              end do
          end do
      enddo

   end if

   if( mvals(0) == 0) then

       ncomp1 = min( nstokes, nspec)

       do l = 0, nlmax
          l1 = l+1
          do is = 1, ncomp1
             cltmp(l,is) = cltmp( l, is) + REAL( current_alm(is)%comp(l1,1))**2
          end do
       enddo

       if( (nstokes >= 2) .and. (nspec >= nstokes+1)) then                                 ! and now do cross-spectra if requested

          ncomp1 = min( nstokes, nspec-nstokes)
          ncomp2 = min( nstokes, nspec-nstokes+1)

          do l = 0, nlmax
             l1 = l+1
             ks = nstokes+1

             do is = 1, ncomp1                                                           ! do all x-spectra
                do js = is+1, ncomp2
                   cltmp(l,ks) = cltmp(l,ks) + REAL( current_alm(is)%comp(l1,1)*current_alm(js)%comp(l1,1))
                   ks = ks+1
                end do
             end do
          end do
       endif

   end if

   call mpi_reduce( cltmp, cl, nspec*(nlmax+1), mpi_double_precision, mpi_sum, root, comm, ierr)
	
   if( myid == root) then
      do i = 1, nspec
         do l = 0, nlmax
            cl(l,i) = cl( l, i)/(2.d0*l+1.d0)
         end do
      end do
   endif

end subroutine collect_cls

subroutine collect_xls( nmaps1, mapnum1, nmaps2, mapnum2, nstokes, nlmax, nmvals, mvals, lda, local_alm1, local_alm2, nxspec, xcl, myid, numprocs, root, comm)

  ! computes power cross-spectra cl from the two sets of alms m-distributed over all procs and stores the output on a
  ! proc root

   !input
   integer(i4b), intent(in) :: nmaps1, mapnum1, nmaps2, mapnum2, nstokes, nlmax, nmvals, nxspec, lda, myid, numprocs, root, comm
   integer(i4b), dimension(0:nmvals-1), intent(in) :: mvals
   complex(DP), dimension(:,:,:,:), POINTER :: local_alm1, local_alm2
   
   !output
   real(DP), dimension(0:nlmax,1:nxspec), intent( out) :: xcl

   !internal
   integer(I4B)  :: i, l, l1, ierr, mstart, mstart1, is, js, ks, ncomp1, ncomp2
   real(DP), dimension(0:nlmax,1:nxspec)  :: xcltmp

   TYPE ROW
       complex(DP), dimension(:,:), POINTER :: comp
   END TYPE ROW

   type(ROW), dimension(1:nstokes) :: current_alm1, current_alm2

   !=======================================================================

   if( (lda == 1) .or. (lda == 2) .or. (lda == 3)) then       ! i.e., healpix convention for the alm array

      do i = 1, lda
         current_alm1(i)%comp => local_alm1(i,0:nlmax,0:nmvals-1,mapnum1)
         current_alm2(i)%comp => local_alm2(i,0:nlmax,0:nmvals-1,mapnum2)
      end do

   else                                       ! i.e., s2hat convention for the alm array

      do i = 1, nstokes
         current_alm1(i)%comp => local_alm1(0:nlmax,0:nmvals-1,i,mapnum1)
         current_alm2(i)%comp => local_alm2(0:nlmax,0:nmvals-1,i,mapnum2)
      end do

   end if

   xcltmp = 0.d0

   mstart = 0
   if( mvals(0) == 0) then
     mstart = 1
   endif

   ncomp1 = min( nstokes, nxspec)

   mstart1 = mstart+1
   do l = mvals(mstart), nlmax
      l1 = l+1

      do is = 1, ncomp1
          xcltmp(l,is) = 2.*SUM( current_alm1(is)%comp(l1,mstart1:nmvals)*CONJG(current_alm2(is)%comp(l1,mstart1:nmvals)))
      end do

   end do

   if( (nstokes >= 2) .and. (nxspec >= nstokes+1)) then                                 ! and now do cross-spectra if requested

       ncomp1 = min( nstokes, nxspec-nstokes-1)
       ncomp2 = min( nstokes, nxspec-nstokes)

       mstart1 = mstart+1
       do l = mvals(mstart), nlmax
          l1 = l+1

          ks = nstokes+1

          do is = 1, ncomp1                                                           ! do all x-spectra
              do js = is+1, ncomp2
                  xcltmp(l,ks) = 2.*SUM( current_alm1(is)%comp(l1,mstart1:nmvals)*CONJG(current_alm2(js)%comp(l1,mstart1:nmvals)))
                  ks = ks+1
                  xcltmp(l,ks) = 2.*SUM( current_alm1(js)%comp(l1,mstart1:nmvals)*CONJG(current_alm2(is)%comp(l1,mstart1:nmvals)))
                  ks = ks+1
              end do
          end do

      end do

   end if

   if( mvals(0) == 0) then

      ncomp1 = min( nstokes, nxspec)

      do l = 0, nlmax
         l1 = l+1

         do is = 1, ncomp1
            xcltmp(l,is) = xcltmp(l,is) + REAL( current_alm1(is)%comp(l1,1)*current_alm2(is)%comp(l1,1))
         end do

      end do

      if( (nstokes >= 2) .and. (nxspec >= nstokes+1)) then                                 ! and now do cross-spectra if requested

         ncomp1 = min( nstokes, nxspec-nstokes-1)
         ncomp2 = min( nstokes, nxspec-nstokes)

         do l = 0, nlmax
            l1 = l+1

            ks = nstokes+1

            do is = 1, ncomp1                                                              ! do all x-spectra
               do js = is+1, ncomp2
                  xcltmp(l,ks) = xcltmp(l,ks) + REAL( current_alm1(is)%comp(l1,1)*current_alm2(js)%comp(l1,1))
                  ks = ks+1
                  xcltmp(l,ks) = xcltmp(l,ks) + REAL( current_alm1(js)%comp(l1,1)*current_alm2(is)%comp(l1,1))
                  ks = ks+1
               end do
            end do
         end do

      end if

   end if

   call mpi_reduce( xcltmp, xcl, nxspec*(nlmax+1), mpi_double_precision, mpi_sum, root, comm, ierr)
	
   if( myid == root) then
       do i = 1, nxspec
          do l = 0, nlmax
             xcl(l,i) = xcl( l, i)/(2.d0*l+1.d0)
          end do
       end do
   endif

end subroutine collect_xls

subroutine distribute_alms( nlmax, nmmax, nmaps, mapnum, nstokes, nmvals, mvals, lda, local_alm, alms, myid, numprocs, root, comm)

    ! distributes alms from proc root to all the others. It is m value-based distribution with the
    ! m values for each proc defined by mvals

    implicit none

    !input
    integer(i4b), intent(in) :: nlmax, nmmax, nmaps, mapnum, nstokes, nmvals, lda, myid, numprocs, root, comm
    integer(i4b), dimension(0:nmvals-1), intent(in) :: mvals
    complex(DP),  dimension(:,:,:), pointer :: alms

    !output
    complex(DP),  dimension(:,:,:,:), pointer :: local_alm

    !internal
    integer(I4B) :: i, ip, im, ierr, status, recv_count
    integer(I4B), dimension(:), allocatable :: mperproc, send_count, send_disps, indx

    ! sort the input alms to match their distribution over procs

    allocate( indx(0:nmmax))

    call fix_mmap( 1, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

    if( myid == root) then
       call sortAlms( nlmax, nmmax, nstokes, indx, lda, alms)
    endif

    deallocate( indx)

    ! dsitribute over procs

    if( myid == root) then
       allocate( mperproc(0:numprocs-1), stat=status)
    endif
    
    call mpi_gather( nmvals, 1, mpi_integer, mperproc, 1, mpi_integer, root, comm, ierr)

    if( myid == root) then

       allocate( send_count(0:numprocs-1))
       allocate( send_disps(0:numprocs-1))

       send_count(0) = (nlmax+1)*mperproc(0)
       send_disps(0) = 0

       do ip = 1, numprocs-1
          send_count(ip) = (nlmax+1)*mperproc(ip)
          send_disps(ip) = send_disps(ip-1)+send_count(ip-1)
       enddo

    endif

    recv_count = (nlmax+1)*nmvals

    if( (lda == 1) .or. (lda == 2) .or. (lda == 3)) then   ! healpix convention

        recv_count = recv_count*nstokes

        if( myid == root) then
            send_count(0:numprocs-1) = nstokes * send_count(0:numprocs-1)
            send_disps(0:numprocs-1) = nstokes * send_disps(0:numprocs-1)
        end if

       call mpi_scatterv( alms( 1, 0, 0), send_count, send_disps, MPI_DOUBLE_COMPLEX, local_alm(1, 0, 0, mapnum), recv_count, MPI_DOUBLE_COMPLEX, root, comm, ierr)

    else                                   ! s2hat convention

       do im = 1, nstokes
          call mpi_scatterv( alms( 0, 0, im), send_count, send_disps, MPI_DOUBLE_COMPLEX, local_alm(0, 0, im, mapnum), recv_count, MPI_DOUBLE_COMPLEX, root, comm, ierr)
       enddo

    end if

    if( myid == root)  deallocate( send_count, send_disps, mperproc)

    return

  end subroutine distribute_alms

  subroutine distribute_map( pixelization, nmaps, mapnum, nstokes, first_ring, last_ring, map_size, local_map, map, myid, numprocs, root, comm)

    ! the routine distributes the map 'map' from proc 'root' to all the other procs (modified MPI HEALpix routine)

    implicit none
    
    ! input
    type( pixeltype), intent(in) :: pixelization
    integer(i4b), intent(in) :: nmaps, mapnum, nstokes, numprocs, myid, root, comm
    integer(i4b), intent(in) :: first_ring, last_ring, map_size
    real(DP), dimension(:,:), pointer  :: map
    ! output
    real(DP), dimension(:,:,:), pointer :: local_map

    ! internal
    integer(i4b) :: i, istokes, j, seg_size, ierr, nsegment, recv_count0, recv_count1
    integer(i4b), dimension(MPI_STATUS_SIZE) :: status 
    integer(i4b), dimension(0:1,0:1) :: segment
    integer(i4b), dimension(:), allocatable :: num_segments, tmp, send_count, send_disps
    integer(i4b), dimension(:,:,:), allocatable :: segments

    ! get pixel ranges
    call find_pixel_numbers( pixelization, first_ring, last_ring, nsegment, segment)

    !Process root distributes the map data to the other proc's

    if( myid == root) then
       allocate( num_segments(0:numprocs-1))
       allocate( tmp(0:numprocs-1))
       allocate( segments(0:1,0:1,0:numprocs-1))
    endif

    call mpi_gather(nsegment, 1, mpi_integer, num_segments, 1, mpi_integer, root, comm, ierr)

    do i = 0, 1
       do j = 0, 1
          call mpi_gather(segment(i,j), 1, mpi_integer, tmp, 1, mpi_integer, root, comm, ierr)
          if( myid == root) segments(i,j,:) = tmp(:)
       enddo
    enddo

    if( myid == root) deallocate( tmp)

    ! send first set of rows now ...

    if( myid == root) then

       allocate( send_count(0:numprocs-1))
       allocate( send_disps(0:numprocs-1))

       do i = 0, numprocs-1
          send_count(i) = segments(0,1,i)-segments(0,0,i)+1
          send_disps(i) = segments(0,0,i)
       enddo

    endif

    recv_count0 = segment(0,1)-segment(0,0)+1

    do istokes = 1, nstokes
       call mpi_scatterv( map(0,istokes), send_count, send_disps, MPI_DOUBLE_PRECISION, local_map(0, istokes, mapnum), recv_count0, MPI_DOUBLE_PRECISION, root, comm, ierr)
    enddo

    ! and the second row now ...

    if( myid == root) then

       do i = 0, numprocs-1

          if( num_segments(i) == 1) then
             send_count(i) = 0
             if( i /= 0) then 
                send_disps(i) = send_disps(i-1)+send_count(i-1) 
             else 
                send_disps(0) = 0
             endif
          else  
             send_count(i) = segments(1,1,i)-segments(1,0,i)+1
             send_disps(i) = segments(1,0,i)
          endif
       enddo

    endif

    recv_count1 = 0
    if( nsegment == 2) recv_count1 = segment(1,1)-segment(1,0)+1

    do istokes = 1, nstokes
       call mpi_scatterv( map(0, istokes), send_count, send_disps, MPI_DOUBLE_PRECISION, local_map(recv_count0, istokes, mapnum), recv_count1, MPI_DOUBLE_PRECISION, root, comm, ierr)
    enddo   
 
    if( myid == root) deallocate( num_segments, segments, send_count, send_disps)

  end subroutine distribute_map

  subroutine distribute_w8ring(npol, first_ring, last_ring, local_w8ring, w8ring, myid, numprocs, root, comm)

    ! distributes the ring weights over the procs (modified MPI HEALPix routine)

    implicit none

    ! input
    integer(i4b), intent( in) :: first_ring, last_ring, npol
    integer(i4b), intent( in) :: myid, numprocs, root, comm
    real(dp), dimension(:,:), pointer :: w8ring
    ! output
    real(dp), dimension(:,:), pointer :: local_w8ring
    
    ! internal
    integer(i4b) :: i, ierr, imap, length, recv_count
    integer(i4b), dimension(:), allocatable           :: leftend, rightend, send_count, send_disps

    if( myid == root) then
      allocate( leftend(0:numprocs-1))
      allocate( rightend(0:numprocs-1))
    endif

    call mpi_gather( first_ring, 1, mpi_integer, leftend, 1, mpi_integer, root, comm, ierr) 
    call mpi_gather( last_ring, 1, mpi_integer, rightend, 1, mpi_integer, root, comm, ierr) 

    if (myid == root) then
      
       allocate( send_count(0:numprocs-1))
       allocate( send_disps(0:numprocs-1))

       do i = 0, numprocs-1
          send_count(i) = rightend(i)-leftend(i)+1
          send_disps(i) = leftend(i)-1
       end do

    endif

    recv_count = last_ring-first_ring+1

    do imap = 1, npol
       call mpi_scatterv( w8ring(1,imap), send_count, send_disps, MPI_DOUBLE_PRECISION, local_w8ring(1, imap), recv_count, MPI_DOUBLE_PRECISION, root, comm, ierr)
    enddo

    if( myid == root) deallocate( leftend, rightend, send_count, send_disps)

  end subroutine distribute_w8ring

  subroutine  gatherRingInfo( first_ring, nrings, keep_it, nringsall, keep_it_all, myid, numprocs, comm)

  ! gathers info determining which pixel ring is to be included

  implicit none

  ! input
   integer(i4b), intent( in) :: first_ring, nrings, nringsall, myid, numprocs, comm
   integer( i4b), dimension(0:nrings-1), intent(in) :: keep_it
  ! output
   integer( i4b), dimension(0:nringsall-1), intent(out) :: keep_it_all

  ! internal
   integer(i4b), dimension(:), allocatable :: frings, ringsperproc
   integer(i4b) :: i, ierr

   allocate( ringsperproc(0:numprocs-1))
   call mpi_allgather( nrings, 1, mpi_integer, ringsperproc, 1, mpi_integer, comm, ierr)
   allocate( frings(0:numprocs-1))
   i = first_ring-1
   call mpi_allgather( i, 1, mpi_integer, frings, 1, mpi_integer, comm, ierr)

   call mpi_allgatherv(keep_it, nrings, mpi_integer, keep_it_all, ringsperproc, frings, mpi_integer, comm, ierr)

   deallocate( ringsperproc, frings)

   return

  end subroutine gatherRingInfo

  function nummvalues( myid, numprocs, nmmax)

    ! determines a number of m modes assigned to a proc myid

    implicit none
  
    ! input
    integer(i4b), intent(in) :: myid, numprocs, nmmax 
    ! output
    integer(i4b) :: nummvalues

    ! internal
    integer(i4b) :: mperproc, npairs

    npairs = (nmmax+1)/2

    mperproc = 2*(npairs/numprocs)

    if( myid < mod( npairs, numprocs)) mperproc = mperproc + 2

    if( myid == mod( npairs, numprocs)) mperproc = mperproc+mod( nmmax+1, 2) ! will add a single m value if nmmax is even

    nummvalues = mperproc

    return

  end function nummvalues

  function nummmodes( nlmax, nmvals, mvals)

    ! computes a number of l,m modes corresponding to a mvals set of nmvals values (added routine)

    ! input
    integer(i4b), intent(in) :: nlmax, nmvals
    integer(i4b), dimension(0:nmvals-1), intent(in) :: mvals
    ! output
    integer(i4b) :: nummmodes

    ! internal 
    integer(i4b) :: i

    nummmodes = 0
    do i = 0, nmvals-1
      nummmodes = nummmodes + nlmax+1-mvals(i)
    end do

  return

  end function nummmodes

  subroutine find_mvalues(myid, numprocs, nmmax, nmvals, mvals)

    ! Assigns the m values to a proc myid to load balance the transforms

    implicit none
  
    ! input
    integer(i4b), intent(in) :: myid, numprocs, nmmax, nmvals
    ! output
    integer(i4b), dimension(0:nmvals-1), intent( out) :: mvals

    ! internal
    integer(i4b) :: i

    do i = 0, min( nmmax/2/numprocs, (nmvals-1)/2)
      mvals( i) = numprocs*i+myid
      mvals( nmvals-1-i) = nmmax-(numprocs*i+myid)
    enddo

    return

  end subroutine find_mvalues
  
  subroutine fix_mmap( dir, nmmax, indx, nmvals, mvals, myid, numprocs, comm)

     ! Determines the global mapping of m values onto procs, i.e., for each m value determines 
     ! its global consecutive number

     implicit none

     ! input
     integer(i4b), intent(in) :: dir, nmmax, nmvals, myid, numprocs, comm
     integer(i4b), dimension(0:nmvals-1), intent(in) :: mvals
     ! output
     integer(i4b), dimension(0:nmmax), intent(out) :: indx

     ! internal
     integer(i4b) :: i, myShift, ierr
     integer(i4b), dimension(:), allocatable :: mperproc, local_indx

     allocate( mperproc(0:numprocs-1))
     call mpi_allgather( nmvals, 1, mpi_integer, mperproc, 1, mpi_integer, comm, ierr)

     myShift = 0
     do i = 0, myid-1
        myShift = myShift + mperproc(i)
     end do

     deallocate( mperproc)

     allocate( local_indx( 0:nmmax))

     local_indx = 0
     if( dir == 1) then         ! for every m value assign the index defining the data storgae  
       do i = 0, nmvals-1
          local_indx( mvals(i)) = myShift+i
       end do
     else                       ! for every index value assign the corresponding m value
       do i = 0, nmvals-1      
          local_indx( myShift+i) = mvals(i)
       end do	 
     endif
	 
     call mpi_allreduce( local_indx, indx, nmmax+1, mpi_integer, mpi_sum, comm, ierr)

     deallocate( local_indx)

     return

  end subroutine fix_mmap
  
  subroutine find_ring_range( pixelization, scan, nmmax, myid, numprocs, first_ring, last_ring, outerror)

  ! Determines the range of rings to be assigned to a proc myid
  ! Takes into account scan pattern - 08/02/2007 RS@APC

    implicit none

    ! input
    type( pixeltype) :: pixelization
    type( scandef) :: scan
    integer(i4b), intent(in) :: nmmax, myid, numprocs
    ! output
    integer(i4b), intent(out) :: first_ring, last_ring, outerror

    ! internal
    integer(i4b) :: iring, jring, iproc, nph, nringsall, nringsobs
    integer(i4b) :: keep_north, keep_south
    real(kind = sp) :: totalcost, costperproc, extracost, coeff
    real(kind = sp), dimension(:), allocatable :: costfun
    integer(i4b), dimension(:), allocatable :: fring, lring, ringnums

    coeff = 4.0   ! coefficient of the load balance cost function

    nringsall = pixelization%nringsall

    ! compute a cost function

    allocate( costfun(1:nringsall))
    allocate( ringnums(1:nringsall))

    ! get pixel location information and find out which rings are to be analysed

    jring = 1
    do iring = 1, nringsall

       keep_north = scan%nfl(iring)
       keep_south = scan%sfl(iring)   
       nph = pixelization%nph(iring)

       extracost = 0.0
       if( keep_north == 1) extracost = 0.5*(coeff*nmmax + nph*log( 1.0*nph)/ALN2_INV)    ! half of the total cost
       if( keep_south == 1) extracost = extracost + 0.5*(coeff*nmmax + nph*log( 1.0*nph)/ALN2_INV)
	   
       if( (keep_north+keep_south) /= 0) then
           if( jring > 0) then
               costfun( jring) = costfun( jring-1) + extracost
           else 
               costfun( jring) = extracost
           endif
	   ringnums( jring) = iring
	   jring = jring+1
        endif
    end do

    nringsobs = jring-1

    totalcost = costfun( nringsobs)

    costperproc = totalcost/numprocs   ! include the equatorial ring

    allocate( fring(0:numprocs-1))

    do iproc = 0, numprocs-1
	
       do jring = 1, nringsobs
          if( iproc*costperproc < costfun(jring)) exit
       enddo

       if( iproc == 0) then
	  fring(iproc) = ringnums( 1)
       else
          if( (iproc*costperproc-costfun(jring-1)) > (costfun(jring)-iproc*costperproc)) then
	     fring( iproc) = ringnums( jring) 
          else
             fring( iproc) = ringnums( jring-1)
          endif
       endif

    end do

    do iproc = 1, numprocs-1
       if( fring( iproc) == fring( iproc-1)) fring( iproc) = fring( iproc)+1 
    end do
	    
    if( fring( numprocs-1) > ringnums( nringsobs)) then
	   
	if( fring( numprocs-2) >= nringsobs) then
           outerror = 1
           return 
	else
           fring( numprocs-1) = nringsobs
	endif

    endif
 
    allocate( lring( 0:numprocs-1))

    do iproc = 0, numprocs -1
	    
       do jring = 1, nringsobs
          if( (iproc+1)*costperproc < costfun(jring)) exit
       end do

      if( iproc == numprocs-1) then
         lring(iproc) = ringnums( nringsobs)
      else
         if( ((iproc+1)*costperproc-costfun(jring-1)) > (costfun(jring)-(iproc+1)*costperproc)) then
            lring(iproc) = ringnums( jring-1)
	 else
            lring(iproc) = ringnums( jring-2)
         endif
      endif

   end do
	    
   do iproc = 0, numprocs-1
      if( lring( iproc) < fring( iproc)) lring(iproc) = fring(iproc)
   end do

   deallocate( costfun)
   deallocate( ringnums)
	
   first_ring = fring( myid)
   last_ring = lring( myid)	
	
   deallocate( lring)
   deallocate( fring)
	
   outerror = 0  ! all fine
	
   return

  end subroutine find_ring_range
  
  subroutine find_pixel_numbers( pixelization, first_ring, last_ring, num_segment, segment)

  ! Determines the map pixel segment(s) corresponding to given ring numbers

    implicit none

    ! input
    type( pixeltype) :: pixelization
    integer(i4b), intent(in) :: first_ring, last_ring
    ! output
    integer(i4b), intent(out) :: num_segment
    integer(i4b), dimension(0:1,0:1), intent(out) :: segment

    num_segment = 2

    segment(0,0) = pixelization%fpix( first_ring) 
    segment(0,1) = pixelization%fpix( last_ring)+pixelization%nph( last_ring)-1

    segment(1,0) = pixelization%npixsall - segment(0,1)-1
    segment(1,1) = pixelization%npixsall - segment(0,0)-1

    if (segment(1,0) < segment(0,1)) then
       ! End of segment 1 is after the start of segment 2 => overlapping
       ! at equator.
       num_segment  = 1
       segment(0,1) = segment(1,1)
       segment(1,0) = 0
       segment(1,1) = 0
    end if

    return

  end subroutine find_pixel_numbers


  function local_map_size( pixelization, first_ring, last_ring)

    ! get a number of pixs in a local map

    ! input
    type( pixeltype) :: pixelization
    integer(i4b), intent(in) :: first_ring, last_ring
    ! output
    integer(i4b) :: local_map_size

    ! internal
    integer(i4b) :: num_segment
    integer(i4b), dimension(0:1,0:1) :: segment

    call find_pixel_numbers( pixelization, first_ring, last_ring, num_segment, segment)

    if (num_segment == 2) then 
       local_map_size = segment(0,1) - segment(0,0) + segment(1,1) - segment(1,0) + 2
    else
       local_map_size = segment(0,1) - segment(0,0) + 1
    end if
     
    return
	
  end function local_map_size

  !**************************************************************************
  !   FFTs-like
  !**************************************************************************  

  subroutine s2hat_ring_synthesis(nlmax,nmmax,datain,nph,dataout,kphi0)

    !=======================================================================
    ! hacked from Healpix and adapted to conform with the more general pixelization schemes
    ! kphi0 (real(dp)) is now an azimuthal offset of 1st pixel in each row
    ! wrt the meridian 0. It is assumed to be in radians (!).
    !
    !             - rs@apc, Oct 2007
    !
    !       RING_SYNTHESIS
    !       called by alm2map
    !       calls     real_fft
    !
    !       NB nph is not necessarily a power of 2
    !
    !=======================================================================


    INTEGER(I4B), INTENT(IN) :: nlmax, nmmax
    INTEGER(I4B), INTENT(IN) :: nph
    REAL(DP),     INTENT(IN) :: kphi0
    COMPLEX(DPC), DIMENSION(0:nmmax), INTENT(IN)  :: datain
    REAL(DP),     DIMENSION(0:nph-1), INTENT(OUT) :: dataout

    INTEGER(I4B) :: im, jm, jmax, iw, ksign, m, k, status
    COMPLEX(DPC), DIMENSION(:), allocatable :: bw
    type(planck_fft2_plan) :: plan
    COMPLEX(DPC) :: dat
    real(DP)     :: arg, dcosarg0, dcosarg1, dsinarg0, dsinarg1, tmp

    character(len=*), parameter :: code = 'S2HAT_RING_SYN'

    !=======================================================================

    ksign = + 1
    k = nph/2
    allocate( bw(0:k), stat = status)
    call assert_alloc(status,code,'bw')
    bw(0:k) = CMPLX(0.0_dp, 0.0_dp, KIND=DP)

    ! if (nmmax > nph) all frequencies [-m,m] are wrapped in [0,nph-1]
    ! otherwise the missing high frequncies are set to zero

    bw(0) = datain(0)

    dcosarg1 = 1.d0
    dsinarg1 = 0.d0

    arg = nph*kphi0
    dcosarg0 = dcos( arg)
    dsinarg0 = dsin( arg)

    k = nmmax/nph

    ! computes only positive freqs as needed for FFT

    do im  = 0, k

       jm = 0
       m = im*nph + jm

       if( (m <= nmmax) .and. (m > 0)) then         ! use only relevant m vals

           bw( jm) = bw( jm) + 2.0*datain(m)*cmplx( dcosarg1, dsinarg1, kind = dp)

       end if

       jmax = min( nph/2, nmmax-im*nph)

       do jm  = 1, jmax                             ! positive m contribution

          m = im*nph + jm

          bw( jm) = bw( jm) + datain(m)*cmplx( dcosarg1, dsinarg1, kind = dp)

       end do

       tmp = dcosarg1*dcosarg0 + dsinarg1*dsinarg0
       dsinarg1 = dcosarg1*dsinarg0 - dsinarg1*dcosarg0
       dcosarg1 = tmp

       jmax = min( nph-1, nmmax-im*nph)

       do jm = (nph+1)/2, jmax                      ! and now negative one (-m)

          m = im*nph + jm

          bw( nph-jm) = bw( nph-jm) + CONJG(datain(m))*cmplx( dcosarg1, -dsinarg1, kind = dp)

       end do

    enddo

    ! and now account for the shift of the first pixel

    dataout(0) = REAL(bw(0), kind=DP)

    arg = ksign*kphi0
    dcosarg0 = dcos( arg)
    dsinarg0 = dsin( arg)

    dcosarg1 = dcosarg0
    dsinarg1 = dsinarg0

    do iw = 1, (nph-1)/2           ! to allow for odd nph values

       dat = bw(iw) * CMPLX( dcosarg1, dsinarg1, kind=DP)

       tmp  = dcosarg1*dcosarg0-dsinarg1*dsinarg0
       dsinarg1 = dcosarg1*dsinarg0+dsinarg1*dcosarg0
       dcosarg1 = tmp

       dataout(iw*2-1) = REAL(dat, kind=DP)
       dataout(iw*2  ) = AIMAG(dat)
    enddo

    ! if nph is even then add the Nyquist mode now i.e., iw = nph/2

    if( 2*(nph/2) == nph) then

        iw = nph/2
        dat = bw(iw) * CMPLX( dcosarg1, dsinarg1, kind=DP)
        dataout(iw*2-1) = REAL(dat, kind=DP)

    endif

    deallocate( bw)

    call make_fft2_plan(plan,length=nph,direction=fft2_backward)
    call real_fft2 (plan, dataout)
    !     ^^^^^^^^^^^^
    call destroy_fft2_plan (plan)

    RETURN

  END subroutine s2hat_ring_synthesis

  subroutine s2hat_ring_analysis(nlmax,nmmax,datain,nph,dataout,kphi0)
    !==========================================================================
    ! hacked from Healpix to conform with the more general pixelization scheme
    ! kphi0 (real(dp)) is now an azimuthal offset of 1st pixel in each row
    ! wrt the meridian 0. It is assumed to be in radians (!).
    !
    !        - rs@apc Oct, 2007
    !
    !     ring_analysis
    !       called by map2alm
    !       calls     real_fft
    !
    !     integrates (data * phi-dependence-of-Ylm) over phi
    !     --> function of m can be computed by FFT
    !     with  0<= m <= npoints/2 (: Nyquist)
    !     because the data is real the negative m are the conjugate of the
    !     positive ones
    !=======================================================================

    INTEGER(I4B), INTENT(IN) :: nlmax, nmmax
    INTEGER(I4B), INTENT(IN) :: nph
    REAL(DP),     INTENT(IN) :: kphi0
    REAL(DP),     DIMENSION(0:nph-1), INTENT(IN)  :: datain
    COMPLEX(DPC), DIMENSION(0:nmmax), INTENT(OUT) :: dataout

    INTEGER(I4B) :: i,m,im_max,ksign,k
    real(dp) :: dcosarg, dcosarg0, dsinarg, dsinarg0, tmp
    REAL(DP), DIMENSION(0:nph-1) :: data
    type(planck_fft2_plan) :: plan
    real(DP)  :: arg

    !=======================================================================

    ksign = - 1
    data= 0.d0
    data( 0 : nph-1) = datain( 0 : nph-1)

    call make_fft2_plan(plan,length=nph,direction=fft2_forward)
    call real_fft2(plan,data)
    call destroy_fft2_plan(plan)

    arg = ksign * kphi0

    dcosarg0 = dcos( arg)
    dsinarg0 = dsin( arg)

    dcosarg = dcosarg0
    dsinarg = dsinarg0

    ! frequency reassignment and phase shift due to the first pixel offset wrt to the meridian 0

    dataout(0) = cmplx(  data(0), 0.d0, kind = dp)

    do m = 1, nmmax

       k = m - nph*(m/nph)

       if( k > nph/2) then
           k = nph - k
           dataout( m) = cmplx( data( 2*k-1), -data( 2*k), kind = dp)*cmplx( dcosarg, dsinarg, kind = dp)   ! complex conjugate
       else if ( 2*k < nph) then
           if( k > 0) then
              dataout( m) = cmplx( data( 2*k-1), data( 2*k), kind = dp)*cmplx( dcosarg, dsinarg, kind = dp)
           else
              dataout( m) = cmplx( data( 0), 0.d0, kind = dp)*cmplx( dcosarg, dsinarg, kind = dp)
           end if
       else
           ! Nyquist
           dataout( m) = cmplx( data( 2*k-1), 0.d0 , kind = dp)*cmplx( dcosarg, dsinarg, kind = dp)
       end if

       tmp = dcosarg*dcosarg0 - dsinarg*dsinarg0
       dsinarg = dsinarg*dcosarg0 + dcosarg*dsinarg0
       dcosarg = tmp

    end do

    RETURN

  END subroutine s2hat_ring_analysis

  !**************************************************************************
  !             PLM GENERATION
  !**************************************************************************

  subroutine plm_mvalues_gen( pixelization,  scan,  npols, nlmax, nmmax, nmvals, mvals, lda, nplm, plm)

  !========================================================
  ! Precomputes plms for each of Healpix constant theta rings [0,2nsmax] 
  ! and for each of nmvals values of m (as defined in mvals)
  ! and for all \ell [m,nlmax].
  !
  ! (adapted HEALPix routine)
  !========================================================

    type( pixeltype) :: pixelization  
    type( scandef) :: scan 
    integer(i4b)   :: lda, npols, nlmax, nmmax, nmvals
    integer(i4b), dimension( 0:nmvals-1) :: mvals
    integer(i8b) :: nplm
    real(dp), dimension(:,:), pointer :: plm

    INTEGER(I4B) :: l, m, ith, scalem, scalel, nd2, nringsall, nringsobs
    integer(I8B) :: nd1, n_lm, n_plm, i_mm, i_mm1, i_up, i_up1, mshift
    real(DP)     :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)     :: normal_m, lam_lm1m
    real(DP)     :: fm2, fl, flm1, fpol
    real(DP)     :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP)     :: ovflow, unflow, one_on_s2, c_on_s2
    real(DP),     dimension(:,:,:), allocatable :: plm_sub
    real(DP),     dimension(:,:), allocatable :: recfac
    real(DP),     dimension(:),   allocatable :: lam_fact
    real(DP),     dimension(:), allocatable :: mfac

    real(DP),     dimension(:),     allocatable :: normal_l
    integer(i4b)  :: nchunks, chunksize, ichunk, lchk, uchk, ithl, jring, kring
    integer(i4b)  :: nph, i, mindx, ndummy
    real(DP),     dimension(0:SMAXCHK) :: cth, sth
    integer(i4b), dimension(:), allocatable ::  keep_it

    INTEGER(I4B) :: status, ierr
    LOGICAL(LGT) :: polarisation
    character(len=*), parameter :: code = 'PLM_MVALS_GEN'

    TYPE ROW
       REAL(DP), dimension(:), POINTER :: comp
    END TYPE ROW

    type(ROW), dimension(npols) :: local_plm

    !=======================================================================

    if( (lda == 3) .or. (lda == 1)) then       ! i.e., healpix convention for the alm array

       do i = 1, lda
          local_plm(i)%comp => plm(i,0:nplm-1)
       end do

    else                                       ! i.e., s2hat convention for the alm array

       do i = 1, npols
          local_plm(i)%comp => plm(0:nplm-1,i)
       end do

    end if

    ! number of (l,m) with m in mvals and l in [m,L]

    nringsobs = scan%nringsobs
    nringsall = pixelization%nringsall

    n_lm  = nummmodes( nlmax, nmvals, mvals)
    n_plm = n_lm * nringsobs  ! output size of plm 
    nd1 = size(plm, 1)
    nd2 = size(plm, 2)

    if (nd1 < n_plm) then
       print*,code//' > Plm array too small:', nd1, n_plm
       stop
    endif
    if (nd2 /= 1 .and. nd2 /= 3) then
       print*,code//' > Plm array should have dimension 1 or 3:',nd2
       stop
    endif
    polarisation = (nd2 == 3)

    !     --- allocates space for arrays ---
    nchunks   = nringsall/SMAXCHK + 1  ! number of chunks
    chunksize = (nringsall+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')
    if (polarisation) then
       allocate(normal_l(0:nlmax),stat = status)
       call assert_alloc(status,code,'normal_l')
    endif
 
    allocate(recfac(0:1,0:nlmax), plm_sub(1:nd2,0:nlmax,0:chunksize-1), stat = status)
    call assert_alloc(status,code,'recfac & plm_sub')
    if (polarisation) then
       allocate(lam_fact(0:nlmax),stat = status)
       call assert_alloc(status,code,'lam_fact')
    endif

    allocate(keep_it(0:SMAXCHK),stat = status)
    call assert_alloc(status,code,'keep_it')

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax, mfac)
    ! generate Polarization normalisation
    if (polarisation) call gen_normpol(nlmax, normal_l)
    call init_rescale_local()
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    plm = 0.0_dp

    kring = 0
    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nringsall)

       do ith = lchk, uchk
          ithl = ith - lchk !local index

          ! extract pixel location information

          cth( ithl) = pixelization%cth( ith)
          sth( ithl) = pixelization%sth( ith)
          nph = pixelization%nph( ith)
       
          ! extract scan info
          keep_it( ithl) = scan%fl(ith)

       enddo

       mshift = 0

       do mindx = 0, nmvals-1

          m = mvals( mindx)  ! actual m value
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)
          if (polarisation) then
             ! generate Ylm relation factor for degree m
             call gen_lamfac(nlmax, m, lam_fact)
             fm2 = real(m * m, kind = DP)
             normal_m = (2.0_dp * m) * (1 - m)
          endif

          jring = kring
          do ithl = 0, uchk - lchk
       
             if( keep_it( ithl) == 1) then
                ! determine lam_mm
                call mpi_compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)
                ! ---------- l = m ----------
                !           temperature 
                lam_lm = corfac * lam_mm !Actual lam_mm 
                plm_sub(1, m, ithl) = lam_lm

                lam_0 = 0.0_dp
                lam_1 = 1.0_dp
                scalel=0
			 
                cth_ring = cth(ithl)
                one_on_s2 =  1.0_dp / sth(ithl)**2
                c_on_s2 = cth(ithl)*one_on_s2
		  
                lam_2 = cth_ring * lam_1 * recfac(0,m)

                if (polarisation) then
                  fpol = normal_m * normal_l(m) * lam_lm
                  plm_sub(2, m, ithl) =  fpol * ( 0.5_dp - one_on_s2)
                  plm_sub(3, m, ithl) =  fpol * c_on_s2
                  !
                  fm_on_s2     =      m * one_on_s2
                  two_on_s2    = 2.0_dp * one_on_s2
                  two_cth_ring = 2.0_dp * cth_ring
                  b_w          =  c_on_s2 
               endif

               ! ---------- l > m ----------

               do l = m+1, nlmax
                  if (polarisation) lam_lm1m = lam_lm * lam_fact(l)
                  lam_lm = lam_2 * corfac * lam_mm
                  plm_sub(1, l, ithl) = lam_lm
                
                  if (polarisation) then
                     fl = real(l, kind = DP)
                     flm1 = fl - 1.0_dp
                     a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                     a_x =  two_cth_ring * flm1
                     plm_sub(2, l, ithl) =            (   b_w * lam_lm1m - a_w * lam_lm) * normal_l(l)
                     plm_sub(3, l, ithl) = fm_on_s2 * (         lam_lm1m - a_x * lam_lm) * normal_l(l)
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
               
               jring = jring + 1

             endif ! check if the ring observed

          enddo ! loop on rings (ithl)

          ! do memory skipping operations outside inner loops
          jring = kring
          do ith = lchk, uchk
             ithl = ith - lchk
             if( keep_it( ithl) == 1) then
                ! location of Ym,m for ring ith
                i_mm = n_lm * jring + mshift
                i_up = i_mm + nlmax - m                     ! location of Ynlmax,m for ring ith

                i_mm1 = i_mm+1
                i_up1 = i_up+1

                do i=1, nd2
                  local_plm(i)%comp(i_mm1:i_up1) = plm_sub(i, m:nlmax, ithl)
                enddo

                jring = jring + 1

             endif  ! check if the ring observed

          enddo

          mshift  = mshift + nlmax-m+1                   ! offset in the output matrix

       enddo ! loop on m

       kring = jring

    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------

    deallocate (recfac, plm_sub)
    if (polarisation) deallocate(lam_fact)
    deallocate(mfac)
    if (polarisation) deallocate(normal_l)

    return

  end subroutine plm_mvalues_gen

  subroutine mpi_compute_lam_mm(mfac, m, sth, lam_mm, corfac, scalem)

    !=======================================================================
    ! computes lam_mm
    ! the true lam_mm is     lam_mm * corfac
    !
    ! (duplicated HEALpix routine compute_lam_mm)
    !=======================================================================

    integer(I4B),            intent(in)  :: m
    real(DP),                intent(in)  :: sth, mfac
    real(DP),                intent(out) :: lam_mm, corfac
    integer(I4B),            intent(out) :: scalem
    !
    real(DP) :: log2val, dlog2lg

    dlog2lg = real(LOG2LG, kind=DP)

    log2val = mfac + m*log(sth) * ALN2_INV ! log_2(lam_mm)
    scalem = int (log2val / dlog2lg)
    corfac = rescale_tab(max(scalem,RSMIN))
    lam_mm = 2.0_dp **(log2val - scalem * dlog2lg) ! rescaled lam_mm
    if (IAND(m,1)>0) lam_mm = -lam_mm ! negative for odd m

    return
  end subroutine mpi_compute_lam_mm

  subroutine s2hat_gen_mfac( s, m_max, m_fact)

    !=======================================================================
    ! generates factor used in slam_mm calculation
    ! for all m in |s| <= m <= m_max
    ! from HEALPIX now spin number included - rs@apc 19/07/2007
    !=======================================================================
    !
    integer(I4B), intent(IN)                       :: s, m_max
    real(DP),     intent(OUT), dimension(0:m_max)  :: m_fact
    !
    integer(I4B) :: m, spos

    ! fact(m) = fact(m-1) * sqrt( (2m+1) * (2m)/ 2/ (m+s) / (m-s) )

    spos = abs( s)

    m_fact(0:spos-1) = 0.0_dp
    
    m_fact(spos) = sqrt( 2.0*spos+1)/( 2.d0**dble( spos))   ! m = |s|

    do m = spos+1, m_max
       m_fact(m) = m_fact(m-1)*sqrt(dble(2*m+1)/dble(m+s)*dble( m)/dble(m-s)/2.0)
    end do

    ! Log_2 ( fact(m) / sqrt(4 Pi) )

    do m=spos,m_max
       m_fact(m) = log(SQ4PI_INV * m_fact(m)) * ALN2_INV
    enddo

    return

  end subroutine s2hat_gen_mfac

  subroutine s2hat_gen_sfac( spin, s_fact)

    !=======================================================================
    ! generates factors used in slam_mm calculation
    ! for all s such as 0 <= s <= spin and m = spin
    ! from HEALPIX now spin number included - rs@apc 19/07/2007
    !=======================================================================
    !
    integer(I4B), intent(IN)                      :: spin
    real(DP),     intent(OUT), dimension(0:spin)  :: s_fact
    !
    integer(I4B) :: m, s, spos

    ! fact(m) = fact(m-1) * sqrt( (2m+1) * (2m)/ 2/ (m+s) / (m-s) )

    spos = abs( spin)

    s_fact(0) = 1.d0

    do m = 1, spos
       s_fact(0) = s_fact(0)*sqrt(dble(2*m+1)/dble(m)/2.0)   ! s = 0 and l = m = spin
    end do

    do s = 1, spos
       s_fact(s) = s_fact(s-1)*sqrt(dble(spos-s+1)/dble(spos+s))  ! s (<= |spin|), and l = m = spin
    end do

    ! Log_2 ( fact(s) / sqrt(4 Pi) )

    do s=0,spos
       s_fact(s) = log(SQ4PI_INV * s_fact(s)) * ALN2_INV
    enddo

    return

  end subroutine s2hat_gen_sfac

  !=======================================================================

  subroutine s2hat_gen_recfac( s, l_max, m, recfac)

    !=======================================================================
    ! generates recursion factors used to computes the sYlm of degree m
    ! for all l in m <= l <= l_max
    !
    ! added spin number - it computes Lewis, Challinor, Turok (2001) C coefficient
    ! formula (B8) - 13/08/2007 rs@apc
    !
    ! m >= |s|
    !
    !=======================================================================

    integer(I4B), intent(IN)                            :: s, l_max, m
    real(DP),     intent(OUT), dimension(0:1, 0:l_max)  :: recfac

    !

    real(DP) :: fm2, fl2, fs2
    integer(I4B) :: l

    recfac(0:1,0:m) = 0.0_dp
    fm2 = DBLE(m) **2
    fs2 = DBLE(s) **2
    do l = m, l_max
       fl2 = DBLE(l+1) **2
       recfac(0,l) = DSQRT( fl2/(fl2-fs2) * (4.0_dp * fl2 - 1.0_dp) / (fl2-fm2) )
    enddo
    ! put outside the loop because of problem on some compilers
    recfac(1,m:l_max) = 1.0_dp / recfac(0,m:l_max)

    return

  end subroutine s2hat_gen_recfac

  subroutine mpi_compute_slam_mm(mfac, s, m, cth, sth, lam_mm, corfac, scalem)

    !=======================================================================
    ! computes s_lam_mm
    ! the true s_lam_mm is     lam_mm * corfac
    !
    ! (duplicated HEALpix routine compute_lam_mm
    !  with the spin number added (formula B7 of Lewis et al (2001) 
    !  - rs@apc, 18/07/2007)
    !
    ! (m >= |s|)
    !
    !=======================================================================

    integer(I4B),            intent(in)  :: s, m
    real(DP),                intent(in)  :: cth, sth, mfac
    real(DP),                intent(out) :: lam_mm, corfac
    integer(I4B),            intent(out) :: scalem
    !
    real(DP) :: log2val, dlog2lg

    dlog2lg = real(LOG2LG, kind=DP)

    log2val = mfac + (m*log(sth) + 0.5d0*s*(log(1.d0-cth)-log(1.d0+cth))) * ALN2_INV ! log_2(lam_mm)
    scalem = int (log2val / dlog2lg)
    corfac = rescale_tab(max(scalem,RSMIN))
    lam_mm = 2.0_dp **(log2val - scalem * dlog2lg) ! rescaled lam_mm
    if (IAND(m,1)>0) lam_mm = -lam_mm ! negative for odd m

    return

  end subroutine mpi_compute_slam_mm

  subroutine mpi_do_lam_lm_pol(lmax, m, cth, sth, mfac, recfac, lam_fact, lam_lm)
    !=======================================================================
    ! computes temperature&polar lambda_lm(p,theta) for all l in [m,lmax] for a given m, and given theta
    ! input: lmax, m, cos(theta), sin(theta)
    !        mfac: precomputed (by gen_mfac) quantity useful for lambda_mm calculation
    !        recfac: precomputed (by gen_recfac) quantities useful for 
    !            lambda_lm recursion for a given m
    !        lam_fact: precomputed (by gen_lamfac) factor useful for polarised lambda recursion
    ! output: lam_lm for T and P
    ! the routine also needs the array rescale_tac initialized by init_rescale
    !
    ! (duplicated HEALPix routine do_lam_lm_pol)
    !=======================================================================

    integer(I4B),                    intent(in)  :: lmax,  m
    real(DP),                        intent(in)  :: cth, sth, mfac
    real(DP), dimension(0:1,0:lmax), intent(in)  :: recfac
    real(DP), dimension(    0:lmax), intent(in)  :: lam_fact
    real(DP), dimension(1:3,0:lmax), intent(out) :: lam_lm
    !
    real(DP) :: log2val, dlog2lg
    real(DP) :: ovflow, unflow, corfac
    real(DP) :: lam_mm, lam_0, lam_1, lam_2, lam_lm1m
    integer(I4B) :: scalel, l, l_min
    real(DP) :: normal_m, fm2, fl, flm1
    real(DP) :: two_cth, one_on_s2, fm_on_s2, two_on_s2, c_on_s2
    real(DP) :: a_w, a_x, b_w
    !---------------------------------------------------------------

    ! define constants
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    l_min = l_min_ylm(m, sth)
    dlog2lg = real(LOG2LG, kind=DP)
    
    fm2       = real(m * m, kind = DP)
    normal_m  = (2.0_dp * m) * (1 - m)
    two_cth   = 2.0_dp * cth
    one_on_s2 = 1.0_dp / (sth * sth)
    fm_on_s2  =      m * one_on_s2
    two_on_s2 = 2.0_dp * one_on_s2
    c_on_s2   = cth    * one_on_s2
    b_w       =  c_on_s2 
    

    ! computes lamba_mm
    log2val = mfac + m*log(sth) * ALN2_INV ! log_2(lam_mm)
    scalel = int (log2val / dlog2lg)
    corfac = rescale_tab(max(scalel,RSMIN))
    lam_mm = 2.0_dp **(log2val - scalel * dlog2lg) ! rescaled lam_mm
    if (IAND(m,1)>0) lam_mm = -lam_mm ! negative for odd m
    
    lam_lm(1:3,m:lmax) = 0.0_dp
    ! --- l = m ---
    lam_lm(1,m) = corfac * lam_mm !Actual lam_mm 

    if (m >= l_min) then ! skip Ymm if too small
       lam_lm(2,m) =  (normal_m * lam_lm(1,m))  * ( 0.5_dp - one_on_s2 )
       lam_lm(3,m) =  (normal_m * lam_lm(1,m))  *            c_on_s2
    endif

    ! --- l > m ---
    lam_0 = 0.0_dp
    lam_1 = 1.0_dp
    lam_2 = cth * lam_1 * recfac(0,m)

    do l = m+1, lmax
       ! do recursion
       lam_lm1m = lam_lm(1,l-1) * lam_fact(l) ! must be incremented even if not used
       lam_lm(1,l) = lam_2 * corfac * lam_mm
       if (l >= l_min) then
          fl = real(l, kind = DP)
          flm1 = fl - 1.0_dp
          a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
          a_x =  two_cth * flm1
          lam_lm(2,l) =                b_w * lam_lm1m - a_w * lam_lm(1,l)
          lam_lm(3,l) = fm_on_s2 * (         lam_lm1m - a_x * lam_lm(1,l))
       endif
       lam_0 = lam_1 * recfac(1,l-1)
       lam_1 = lam_2
       lam_2 = (cth * lam_1 - lam_0) * recfac(0,l)

       ! do dynamic rescaling
       if (abs(lam_2) > ovflow) then
          lam_1 = lam_1*unflow
          lam_2 = lam_2*unflow
          scalel= scalel + 1
          corfac = rescale_tab(max(scalel,RSMIN))
       elseif (abs(lam_2) < unflow) then
          lam_1 = lam_1*ovflow
          lam_2 = lam_2*ovflow
          scalel= scalel - 1
          corfac = rescale_tab(max(scalel,RSMIN))
       endif
                   
    enddo ! loop on l
  end subroutine mpi_do_lam_lm_pol

  subroutine mpi_do_slam_lm_old( s, lmax, m, cth, sth, mfac, recfac, slam_comb)

    !=======================================================================
    !
    ! computes the spin associated Legendre functions s_lambda_lm(theta) for all l in [m,lmax] 
    !          for a given s, m, and theta
    !
    ! input: s, lmax, m, cos(theta), sin(theta)
    !        mfac: precomputed (by s2hat_gen_mfac) quantity used for s_lambda_mm calculation
    !        recfac: precomputed (by s2hat_gen_recfac) quantities useful for 
    !                s_lambda_lm recursion for a given m 
    !        (NB: these two are different for s < m and s > m cases !);
    !
    ! output: for each l>= m and +/-s two linear combinations of s_lam_lm given by:
    !         (0) : slam_lm + (-1)^s (-s)lam_lm
    !         (1) : slam_lm - (-1)^s (-s)lam_lm
    !
    ! the routine also needs the array rescale_tab initialized by init_rescale in the
    ! calling routine.
    !
    !
    ! (follows the logic of the healpix do_lam_lm routine)
    !=======================================================================

    integer(I4B),                    intent(in)  :: s, lmax,  m
    real(DP),                        intent(in)  :: cth, sth, mfac
    real(DP), dimension(0:1,0:lmax), intent(in)  :: recfac
    real(DP), dimension(1:2,0:lmax), intent(out) :: slam_comb
    !
    real(DP) :: ovflow, unflow, spinPow, signFix
    real(DP), dimension(0:1) :: corfac, slam_mm, slam_lm, slam_0, slam_1, slam_2
    real(DP), dimension(0:1) :: fact1, fact2
    integer(I4B) :: l, is
    integer(i4b), dimension( 0:1) :: scalel, scalem

    !---------------------------------------------------------------

    ! define constants
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    
    spinPow = (-1.d0) ** s

    if( s <= m) then

       slam_comb(1:2,m:lmax) = 0.0_dp
    
       call mpi_compute_slam_mm( mfac, s, m, cth, sth, slam_mm(0), corfac(0), scalem(0))
       call mpi_compute_slam_mm( mfac, -s, m, cth, sth, slam_mm(1), corfac(1), scalem(1))

       !           ---------- l = m ----------

       slam_lm(0) = corfac(0) * slam_mm(0)            ! Actual slam_mm with s > 0
       slam_lm(1) = corfac(1) * slam_mm(1)            ! Actual slam_mm with s < 0

       slam_comb(1,m) = slam_lm(0) + spinPow * slam_lm(1)
       slam_comb(2,m) = slam_lm(0) - spinPow * slam_lm(1)

       !           ---------- l > m ----------

       slam_0(0:1) = 0.0_dp
       slam_1(0:1) = 1.0_dp

       scalel(0:1) = 0

       slam_2(0) = (cth + s/dble(m+1)) * slam_1(0) * recfac(0,m)   ! spin > 0, l = m
       slam_2(1) = (cth - s/dble(m+1)) * slam_1(1) * recfac(0,m)   ! spin < 0, l = m

       do l = m+1, lmax

          fact1(0) = cth + s*m/dble(l+1)/dble(l)
          fact1(1) = cth - s*m/dble(l+1)/dble(l)

          slam_lm(0) = slam_2(0) * corfac(0) * slam_mm(0)
          slam_lm(1) = slam_2(1) * corfac(1) * slam_mm(1)

          slam_comb(1, l) = slam_lm(0) + spinPow * slam_lm(1)
          slam_comb(2, l) = slam_lm(0) - spinPow * slam_lm(1)

          ! actual sPlm reccurence

          do is = 0, 1

             slam_0(is) = slam_1(is) * recfac(1,l-1)
             slam_1(is) = slam_2(is)
             slam_2(is) = (fact1(is) * slam_1(is) - slam_0(is)) * recfac(0,l)

             ! and rescaling if needed

             if (abs(slam_2(is)) > OVFLOW) then
                slam_1(is) = slam_1(is)*UNFLOW
                slam_2(is) = slam_2(is)*UNFLOW
                scalel(is) = scalel(is) + 1
                corfac(is) = rescale_tab(max(scalel(is)+scalem(is),RSMIN))
             elseif (abs(slam_2(is)) < UNFLOW) then
                slam_1(is) = slam_1(is)*OVFLOW
                slam_2(is) = slam_2(is)*OVFLOW
                scalel(is) = scalel(is) - 1
                corfac(is) = rescale_tab(max(scalel(is)+scalem(is),RSMIN))
             endif

          enddo ! loop over is

       end do ! over l

    else  ! now the case when m < s i.e. recurrence wrt l for fixed m = s and for s = m ...

       signFix = (-1.0)**(s+m)

       call mpi_compute_slam_mm( mfac, m, s, cth, sth, slam_mm(0), corfac(0), scalem(0))  !  -> mPss
       call mpi_compute_slam_mm( mfac, -m, s, cth, sth, slam_mm(1), corfac(1), scalem(1)) !  -> -mPss

       !           ---------- l = spin ----------

       slam_lm(0) = signFix * corfac(0) * slam_mm(0)                   ! Actual lam_mm with s > 0 and m < s
       slam_lm(1) = corfac(1) * slam_mm(1)                             ! Actual lam_mm with s < 0 and m < -s

       slam_comb(1,s) = slam_lm(0) + spinPow * slam_lm(1)
       slam_comb(2,s) = slam_lm(0) - spinPow * slam_lm(1)

       !           ---------- l > spin ----------

       slam_0(0:1) = 0.0_dp
       slam_1(0:1) = 1.0_dp

       scalel(0:1) = 0

       slam_2(0) = (cth + m/dble(s+1)) * slam_1(0) * recfac(0,s)   ! s > 0, l = m
       slam_2(1) = (cth - m/dble(s+1)) * slam_1(1) * recfac(0,s)   ! s < 0, l = m

       do l = s+1, lmax

          fact1(0) = cth + s*m/dble(l+1)/dble(l)
          fact1(1) = cth - s*m/dble(l+1)/dble(l)

          slam_lm(0) = signFix * slam_2(0) * corfac(0) * slam_mm(0)
          slam_lm(1) = slam_2(1) * corfac(1) * slam_mm(1)

          slam_comb(1, l) = slam_lm(0) + spinPow * slam_lm(1)
          slam_comb(2, l) = slam_lm(0) - spinPow * slam_lm(1)

          ! actual sPlm reccurence

          do is = 0, 1

             slam_0(is) = slam_1(is) * recfac(1,l-1)
             slam_1(is) = slam_2(is)
             slam_2(is) = (fact1(is) * slam_1(is) - slam_0(is)) * recfac(0,l)

             ! and rescaling if needed

             if (abs(slam_2(is)) > OVFLOW) then
                slam_1(is) = slam_1(is)*UNFLOW
                slam_2(is) = slam_2(is)*UNFLOW
                scalel(is) = scalel(is) + 1
                corfac(is) = rescale_tab(max(scalel(is)+scalem(is),RSMIN))
             elseif (abs(slam_2(is)) < UNFLOW) then
                slam_1(is) = slam_1(is)*OVFLOW
                slam_2(is) = slam_2(is)*OVFLOW
                scalel(is) = scalel(is) - 1
                corfac(is) = rescale_tab(max(scalel(is)+scalem(is),RSMIN))
             endif

          enddo ! over is

       enddo ! over l

    endif ! s > m case

  end subroutine mpi_do_slam_lm_old

  subroutine mpi_do_slam_lm( s_in, lmax, m_in, cth, sth, mfac, recfac, slam_comb)

    !=======================================================================
    !
    ! computes the spin associated Legendre functions s_lambda_lm(theta) for all l in [m,lmax] 
    !          for a given s, m, and theta
    !
    ! input: s, lmax, m, cos(theta), sin(theta)
    !        mfac: precomputed (by s2hat_gen_mfac) quantity used for s_lambda_mm calculation
    !        recfac: precomputed (by s2hat_gen_recfac) quantities useful for 
    !                s_lambda_lm recursion for a given m 
    !        (NB: these two are different for s < m and s > m cases !);
    !
    ! output: for each l>= m and +/-s two linear combinations of s_lam_lm given by:
    !         (0) : slam_lm + (-1)^s (-s)lam_lm
    !         (1) : slam_lm - (-1)^s (-s)lam_lm
    !
    ! the routine also needs the array rescale_tab initialized by init_rescale in the
    ! calling routine.
    !
    !
    ! (follows the logic of the healpix do_lam_lm routine)
    !=======================================================================

    ! input
    integer(I4B),                    intent(in)  :: s_in, lmax,  m_in
    real(DP),                        intent(in)  :: cth, sth, mfac
    real(DP), dimension(0:1,0:lmax), intent(in)  :: recfac

    ! output
    real(DP), dimension(1:2,0:lmax), intent(out) :: slam_comb

    ! internal
    real(DP) :: ovflow, unflow, spinPow, signFix
    real(DP), dimension(0:1) :: corfac, slam_mm, slam_lm, slam_0, slam_1, slam_2
    real(DP), dimension(0:1) :: fact1, fact2
    integer(I4B) :: l, is, s, m
    integer(i4b), dimension( 0:1) :: scalel, scalem

    !---------------------------------------------------------------

    ! define constants
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    
    spinPow = (-1.d0) ** s_in

    if( s_in <= m_in) then

       signFix = 1.d0

       m = m_in
       s = s_in

    else

       m = s_in ! m <-> s interchanged
       s = m_in

       signFix = (-1.d0)**(s+m)

    end if

    slam_comb(1:2,m:lmax) = 0.0_dp
    
    call mpi_compute_slam_mm( mfac, s, m, cth, sth, slam_mm(0), corfac(0), scalem(0))
    call mpi_compute_slam_mm( mfac, -s, m, cth, sth, slam_mm(1), corfac(1), scalem(1))

    !           ---------- l = m ----------

    slam_lm(0) = signFix * corfac(0) * slam_mm(0)            ! Actual slam_mm with s > 0
    slam_lm(1) = corfac(1) * slam_mm(1)                      ! Actual slam_mm with s < 0

    slam_comb(1,m) = slam_lm(0) + spinPow * slam_lm(1)
    slam_comb(2,m) = slam_lm(0) - spinPow * slam_lm(1)

    !           ---------- l > m ----------

    slam_0(0:1) = 0.0_dp
    slam_1(0:1) = 1.0_dp

    scalel(0:1) = 0

    slam_2(0) = (cth + s/dble(m+1)) * slam_1(0) * recfac(0,m)   ! spin > 0, l = m
    slam_2(1) = (cth - s/dble(m+1)) * slam_1(1) * recfac(0,m)   ! spin < 0, l = m

    do l = m+1, lmax

       fact1(0) = cth + s*m/dble(l+1)/dble(l)
       fact1(1) = cth - s*m/dble(l+1)/dble(l)

       slam_lm(0) = signFix * slam_2(0) * corfac(0) * slam_mm(0)
       slam_lm(1) = slam_2(1) * corfac(1) * slam_mm(1)

       slam_comb(1, l) = slam_lm(0) + spinPow * slam_lm(1)
       slam_comb(2, l) = slam_lm(0) - spinPow * slam_lm(1)

       ! actual sPlm reccurence

       do is = 0, 1

          slam_0(is) = slam_1(is) * recfac(1,l-1)
          slam_1(is) = slam_2(is)
          slam_2(is) = (fact1(is) * slam_1(is) - slam_0(is)) * recfac(0,l)

          ! and rescaling if needed

          if (abs(slam_2(is)) > OVFLOW) then
             slam_1(is) = slam_1(is)*UNFLOW
             slam_2(is) = slam_2(is)*UNFLOW
             scalel(is) = scalel(is) + 1
             corfac(is) = rescale_tab(max(scalel(is)+scalem(is),RSMIN))
          elseif (abs(slam_2(is)) < UNFLOW) then
             slam_1(is) = slam_1(is)*OVFLOW
             slam_2(is) = slam_2(is)*OVFLOW
             scalel(is) = scalel(is) - 1
             corfac(is) = rescale_tab(max(scalel(is)+scalem(is),RSMIN))
          endif

       enddo ! loop over is

    end do ! over l

  end subroutine mpi_do_slam_lm

  subroutine init_rescale_local()

    ! (duplicated HEALPix routine init_rescale)
    !
    !================================================
    ! local variables
    integer(i4b) :: s, smax
    real(dp) :: logOVFLOW
    character(len=*), parameter :: code = 'gen_rescale'
    !------------------------------------------------
    logOVFLOW=log(FL_LARGE)
    smax = INT( log(MAX_DP) / logOVFLOW )

    if (smax > (RSMAX-1)) then
       print*,'Array rescale_tab too small in '//code
       print*,smax ,'>', RSMAX
       stop
    endif

    rescale_tab(RSMIN:RSMAX) = 0.0_dp
    do s = -smax, smax
       rescale_tab(s) = FL_LARGE ** s
    enddo
    rescale_tab(0) = 1.0_dp

    return
  end subroutine init_rescale_local

  subroutine sortAlms( nlmax, nmmax, nstokes, mvect, lda, alms)

     ! sorts the output alm array with respect to m - based on the quicksort routine
     ! from Numerical Recipies

     ! input
     implicit none
     integer(i4b), intent(in) :: nlmax, nmmax, nstokes, lda
 
     ! input/output
     integer(i4b), dimension(0:nmmax), intent(inout) :: mvect
     complex(DP),  dimension(:,:,:), pointer :: alms

     ! internal
     integer(i4b), parameter :: M = 7, NSTACK = 50
     integer(i4b) :: i, ir, imap, j, k, l, jstack, ma, itmp
     integer(i4b), dimension(0:NSTACK) :: istack
     complex(dp), dimension( 0:nlmax, 1:nstokes) :: bcolumn
     complex(dp), dimension(0:nlmax) :: coltmp

    TYPE ROW
       complex(DP), dimension(:,:), POINTER :: comp
    END TYPE ROW

    type(ROW), dimension(1:nstokes) :: local_alm

    !=======================================================================

    if( (lda == 3) .or. (lda == 2) .or. (lda == 1)) then       ! i.e., healpix convention for the alm array

       do i = 1, lda
          local_alm(i)%comp => alms(i,0:nlmax,0:nmmax)
       end do

     else                                       ! i.e., s2hat convention for the alm array

       do i = 1, nstokes
          local_alm(i)%comp => alms(0:nlmax,0:nmmax,i)
       end do

     end if

     l = 0
     ir = nmmax
     jstack = 0

     do

        if (ir-l < M) then
           do j = l+1, ir
 	      ma = mvect(j)
              do imap = 1, nstokes
                 bcolumn( 0:nlmax, imap) = local_alm(imap)%comp( 1:nlmax+1, j+1)
              enddo

	      do i = j-1, l, -1
	         if (mvect(i) <= ma) exit
                 mvect(i+1) = mvect(i)
                 do imap = 1, nstokes
	            local_alm(imap)%comp( 1:nlmax+1, i+2) = local_alm(imap)%comp( 1:nlmax+1, i+1)
                 enddo
              enddo

  	      mvect(i+1) = ma
              do imap = 1, nstokes
                 local_alm(imap)%comp( 1:nlmax+1,i+2) = bcolumn( 0:nlmax,imap)
              enddo
           enddo

           if( jstack == 0) return

           ir = istack( jstack)
           l = istack( jstack-1)
           jstack = jstack-2
        else
           k=(l+ir)/2

           ! swap k <-> l+1

           itmp = mvect(k)
           mvect(k) =  mvect(l+1)
           mvect(l+1) = itmp

           do imap = 1, nstokes
              coltmp(0:nlmax) = local_alm(imap)%comp(1:nlmax+1, k+1)
              local_alm(imap)%comp( 1:nlmax+1, k+1) = local_alm(imap)%comp( 1:nlmax+1, l+2)
              local_alm(imap)%comp( 1:nlmax+1, l+2) = coltmp(0:nlmax)
           enddo

           if (mvect(l) > mvect(ir)) then

              ! swap l <-> ir

              itmp = mvect(l)
              mvect(l) = mvect( ir)
              mvect( ir) = itmp

              do imap = 1, nstokes
	         coltmp(0:nlmax) = local_alm(imap)%comp( 1:nlmax+1, l+1)
                 local_alm(imap)%comp( 1:nlmax+1, l+1) = local_alm(imap)%comp( 1:nlmax+1, ir+1)
                 local_alm(imap)%comp( 1:nlmax+1, ir+1) = coltmp(0:nlmax)
              enddo
           endif

           if (mvect(l+1) > mvect( ir)) then

              ! swap: ir <-> l+1

              itmp = mvect(l+1)
              mvect(l+1) = mvect( ir)
              mvect( ir) = itmp

              do imap = 1, nstokes
                 coltmp( 0:nlmax) = local_alm(imap)%comp( 1:nlmax+1, l+2)
                 local_alm(imap)%comp( 1:nlmax+1, l+2) = local_alm(imap)%comp(1:nlmax+1, ir+1)
                 local_alm(imap)%comp( 1:nlmax+1, ir+1) = coltmp(0:nlmax)
              enddo
           endif

           if (mvect(l) > mvect(l+1)) then

              ! swap: l+1 <-> l

              itmp = mvect(l)
              mvect(l) = mvect(l+1)
              mvect(l+1) = itmp

              do imap = 1, nstokes
                 coltmp(0:nlmax) = local_alm(imap)%comp( 1:nlmax+1, l+1)
                 local_alm(imap)%comp( 1:nlmax+1, l+1) =  local_alm(imap)%comp( 1:nlmax+1, l+2)
                 local_alm(imap)%comp( 1:nlmax+1, l+2) = coltmp(0:nlmax)
              enddo
           endif

           i = l+1
           j = ir

           ma = mvect(l+1)
           do imap = 1, nstokes
              bcolumn(0:nlmax,imap) = local_alm(imap)%comp(1:nlmax+1, l+2)
           enddo

           do 
              do 
                 i = i+1 
                 if( mvect(i) >= ma) exit
              enddo

              do
                 j = j-1 
                 if( mvect(j) <= ma) exit
              enddo

	      if (j < i) exit

              ! swap: i <-> j

              itmp = mvect(i)
              mvect(i) = mvect(j)
              mvect(j) = itmp

              do imap = 1, nstokes
                 coltmp(0:nlmax) = local_alm(imap)%comp(1:nlmax+1, i+1)
                 local_alm(imap)%comp(1:nlmax+1, i+1) = local_alm(imap)%comp(1:nlmax+1, j+1)
                 local_alm(imap)%comp(1:nlmax+1, j+1) = coltmp( 0:nlmax)
              enddo
           enddo

           mvect(l+1) = mvect(j)
           mvect(j) = ma
           do imap = 1, nstokes
              local_alm(imap)%comp(1:nlmax+1, l+2) = local_alm(imap)%comp( 1:nlmax+1, j+1)
              local_alm(imap)%comp(1:nlmax+1, j+1) = bcolumn( 0:nlmax, imap)
           enddo

           jstack = jstack + 2

           if (jstack > NSTACK) then
              write(*,*)'NSTACK too small in sortAlms.'
              stop
           endif

           if (ir-i+1 >= j-l) then
	          istack( jstack) = ir
	          istack( jstack-1) = i
	          ir = j-1
           else 
	          istack( jstack) = j-1
	          istack( jstack-1) = l
	          l = i
           endif
        endif
     enddo

     return

  end subroutine sortAlms

end module s2hat_toolbox_mod
