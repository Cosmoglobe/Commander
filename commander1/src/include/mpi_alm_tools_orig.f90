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
! Written by Hans Kristian Eriksen and Snorre Boasson, 
! but copying large parts from the serial HEALPix code.
module mpi_alm_tools

  use healpix_types
  use healpix_fft
  use misc_utils
  use alm_tools
  use pix_tools
  
  implicit none

  include "mpif.h"

  ! define large and small numbers used to renormalise the recursion 
  ! on the Legendre Polynomials
  integer(I4B),      private, parameter :: LOG2LG   = 100
  real(KIND=DP),     private, parameter :: FL_LARGE = 2.0_dp **   LOG2LG
  real(KIND=DP),     private, parameter :: FL_SMALL = 2.0_dp ** (-LOG2LG)
  ! declare array for dynamic rescaling of the Ylm
  integer(kind=i4b), private, parameter :: RSMAX = 20, RSMIN = -20
  real(dp),          private, dimension(RSMIN:RSMAX) :: rescale_tab
  real(DP),          private, parameter :: ALN2_INV = 1.4426950408889634073599246810_dp ! 1/log(2)

  ! Do not change SMAXCHK -- this will affect the load balancing, and the
  ! corresponding constants in find_rings must be re-adjusted
  integer(i4b),      private, parameter :: SMAXCHK = 50 ! maximum size of chunk (in ring pairs)
  

  !--------------------------------------------------------------------------

  ! Internal parameters for the mpi_alm_tools module
  integer(i4b), private :: nsmax, nlmax, nmmax, precompute_plms
  logical(lgt), private :: polarization
  real(dp),     dimension(2), private :: zbounds

  ! MPI and parallelization variables
  integer(i4b), private :: myid, numprocs, comm, ierr, root = 0
  integer(i4b), private :: first_ring, last_ring, num_segment, map_size, nmaps

  real(dp),     allocatable, dimension(:,:),     private :: local_map, buffer
  complex(dpc), allocatable, dimension(:,:,:),   private :: local_alms, alms_work
  real(dp),     allocatable, dimension(:,:),     private :: local_w8ring, plm
  integer(i4b),              dimension(0:1,0:1), private :: segment

  integer(i4b), allocatable, dimension(:),       private :: num_segments, map_sizes
  integer(i4b), allocatable, dimension(:,:),     private :: ring_procs
  integer(i4b), allocatable, dimension(:,:,:),   private :: segments

  interface mpi_map2alm_simple
     module procedure mpi_map2alm_simple_s, mpi_map2alm_simple_d
  end interface

  interface mpi_alm2map_simple
     module procedure mpi_alm2map_simple_s, mpi_alm2map_simple_d
  end interface

  interface mpi_map2alm
     module procedure mpi_map2alm_s, mpi_map2alm_d
  end interface

  interface mpi_alm2map
     module procedure mpi_alm2map_s, mpi_alm2map_d
  end interface

  interface mpi_map2alm_dist
     module procedure mpi_map2alm_dist_s, mpi_map2alm_dist_d
  end interface

  interface mpi_alm2map_dist
     module procedure mpi_alm2map_dist_s, mpi_alm2map_dist_d
  end interface

  interface mpi_map2alm_dist_slave
     module procedure mpi_map2alm_dist_slave_s, mpi_map2alm_dist_slave_d
  end interface

  interface mpi_alm2map_dist_slave
     module procedure mpi_alm2map_dist_slave_s, mpi_alm2map_dist_slave_d
  end interface

  interface distribute_map
     module procedure distribute_map_s, distribute_map_d
  end interface

  interface gather_map
     module procedure gather_map_s, gather_map_d
  end interface


contains

  ! ***************************************************************
  !            Initialization and destruction routines
  ! ***************************************************************

  subroutine mpi_initialize_alm_tools(comm_in, nsmax_in, nlmax_in, &
       & nmmax_in, zbounds_in, polarization_in, &
       & precompute_plms_in, w8ring, mask)
    implicit none
    
    integer(i4b),                   intent(in)           :: comm_in
    integer(i4b),                   intent(in), optional :: nsmax_in, nlmax_in, nmmax_in
    integer(i4b),                   intent(in), optional :: precompute_plms_in
    logical(lgt),                   intent(in), optional :: polarization_in
    real(dp),     dimension(2),     intent(in), optional :: zbounds_in
    real(dp),     dimension(1:,1:), intent(in), optional :: w8ring
    real(dp),     dimension(0:),    intent(in), optional :: mask

    integer(i4b) :: i, s, smax, num_inc_rings
    real(dp)     :: logOVFLOW
    real(dp), dimension(3) :: vector
    
    comm            = comm_in

    ! Initialize MPI variables
    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, numprocs, ierr)

    if (myid == root) then
       nsmax            = nsmax_in
       nlmax            = nlmax_in
       nmmax            = nmmax_in
       zbounds          = zbounds_in
       polarization     = polarization_in
       precompute_plms  = precompute_plms_in

       if (present(mask)) then
          num_inc_rings = 0
          do i = 0, 6*nsmax**2-1
             if (mask(i)==1 .or. mask(12*nsmax**2-i-1)==1) then
                call pix2vec_ring(nsmax, i, vector)
                num_inc_rings = max(num_inc_rings, ring_num(nsmax, vector(3)))
             end if
          end do
       else
          num_inc_rings = 2*nsmax
       end if
       num_inc_rings = 2*nsmax
    end if

    call mpi_bcast(num_inc_rings,    1, MPI_INTEGER,          root, comm, ierr)
    call mpi_bcast(nsmax,            1, MPI_INTEGER,          root, comm, ierr)
    call mpi_bcast(nlmax,            1, MPI_INTEGER,          root, comm, ierr)
    call mpi_bcast(nmmax,            1, MPI_INTEGER,          root, comm, ierr)
    call mpi_bcast(precompute_plms,  1, MPI_INTEGER,          root, comm, ierr)
    call mpi_bcast(zbounds,          2, MPI_DOUBLE_PRECISION, root, comm, ierr)
    call mpi_bcast(polarization,     1, MPI_LOGICAL,          root, comm, ierr)

    ! Find which rings to compute
    call find_rings(myid, num_inc_rings, first_ring, last_ring)

!    write(*,*) myid, first_ring, last_ring
!    call mpi_finalize(ierr)
!    stop

    ! Find corresponding pixel numbers
    call find_pixel_numbers(first_ring, last_ring, num_segment, map_size, segment)

    if (polarization) then
       nmaps = 3
    else
       nmaps = 1
    end if

    ! Allocate work space for both map and ring weights
    allocate(local_map(0:map_size-1,1:nmaps))
    allocate(local_alms(1:nmaps,0:nlmax,0:nmmax))
    allocate(alms_work(1:nmaps,0:nlmax,0:nmmax))
    allocate(local_w8ring(1:last_ring-first_ring+1,1:nmaps)) 
    
    if (myid == root) then
       allocate(ring_procs(2,0:numprocs-1))
       allocate(num_segments(0:numprocs-1))
       allocate(map_sizes(0:numprocs-1))
       allocate(segments(0:1,0:1,0:numprocs-1))

       do i = 0, numprocs-1
          call find_rings(i, num_inc_rings, ring_procs(1,i), ring_procs(2,i))
          call find_pixel_numbers(ring_procs(1,i), ring_procs(2,i), &
               & num_segments(i), map_sizes(i), segments(0:1,0:1,i))
       end do

       allocate(buffer(0:2*maxval(segments(0,1,:)-segments(0,0,:))+1,1:nmaps))

    end if

    ! Distribute ring weights
    if (myid == root) then
       call distribute_w8ring(local_w8ring, w8ring)
    else
       call distribute_w8ring(local_w8ring)
    end if

    
    ! Initialize rescale_tab
    logOVFLOW=log(FL_LARGE)
    smax = INT( log(MAX_DP) / logOVFLOW )

    rescale_tab(RSMIN:RSMAX) = 0.0_dp
    do s = -smax, smax
       rescale_tab(s) = FL_LARGE ** s
    enddo
    rescale_tab(0) = 1.0_dp

    ! Pre-compute if requested
    if (precompute_plms == 1) then
       call mpi_plm_gen(.false.)
    else if (precompute_plms == 2) then
       call mpi_plm_gen(.true.)
    end if

  end subroutine mpi_initialize_alm_tools


  subroutine mpi_cleanup_alm_tools
    implicit none

    if (allocated(alms_work))       deallocate(alms_work)
    if (allocated(local_alms))      deallocate(local_alms)
    if (allocated(local_map))       deallocate(local_map)
    if (allocated(plm))             deallocate(plm)
    if (allocated(local_w8ring))    deallocate(local_w8ring)

    if (myid == root) then
       if (allocated(ring_procs))   deallocate(ring_procs)
       if (allocated(num_segments)) deallocate(num_segments)
       if (allocated(map_sizes))    deallocate(map_sizes)
       if (allocated(segments))     deallocate(segments)
       if (allocated(buffer))       deallocate(buffer)
    end if

  end subroutine mpi_cleanup_alm_tools


  subroutine get_pixels(pixels)
    implicit none

    integer(i4b), pointer, dimension(:) :: pixels

    integer(i4b) :: i, j, k

    allocate(pixels(0:map_size-1))

    j = 0
    do k = 0, num_segment-1
       do i = segment(k,0), segment(k,1)
          pixels(j) = i
          j         = j+1
       end do
    end do

  end subroutine get_pixels


  ! ***************************************************************
  !            Include the SP/DP specific routines
  ! ***************************************************************

  include 'mpi_alm_tools_ss.f90'

  include 'mpi_alm_tools_dd.f90'


  ! ***************************************************************
  !            Wrappers for the computational routines
  ! ***************************************************************

  
  subroutine mpi_map2alm_slave
    implicit none

    call distribute_map_slave

    if (polarization) then
       ! Polarization routines
       if (precompute_plms == 1) then
          call mpi_map2alm_pol_pre1
       else if (precompute_plms == 2) then
          call mpi_map2alm_pol_pre2
       else
          call mpi_map2alm_pol
       end if
    else
       ! Temperature routines
       if (precompute_plms == 1 .or. precompute_plms == 2) then
          call mpi_map2alm_sc_pre
       else
          call mpi_map2alm_sc
       end if
    end if

    call mpi_reduce(local_alms, alms_work, size(local_alms), &
         & MPI_DOUBLE_COMPLEX, MPI_SUM, root, comm, ierr) 

  end subroutine mpi_map2alm_slave


  subroutine mpi_alm2map_slave
    implicit none

    call mpi_bcast(local_alms, size(local_alms), MPI_DOUBLE_COMPLEX, root, comm, ierr)

    if (polarization) then
       ! Polarization routines
       if (precompute_plms == 1) then
          call mpi_alm2map_pol_pre1
       else if (precompute_plms == 2) then
          call mpi_alm2map_pol_pre2
       else
          call mpi_alm2map_pol
       end if
    else
       ! Temperature routines
       if (precompute_plms == 1 .or. precompute_plms == 2) then
          call mpi_alm2map_sc_pre
       else
          call mpi_alm2map_sc
       end if
    end if

    call gather_map_slave

  end subroutine mpi_alm2map_slave


  ! *****************************************************************
  ! *           Distribution routines
  ! *****************************************************************

  subroutine distribute_map_slave
    ! Distributes the map data
    implicit none

    integer(i4b) :: i, seg_size, ierr

    integer(i4b), dimension(MPI_STATUS_SIZE)          :: status

    !Processes 1 and up recieves map data
    seg_size = segment(0,1)-segment(0,0)+1
    if (num_segment == 2) then
       
       call mpi_recv(local_map(0:2*seg_size-1,1:nmaps), 2*seg_size*nmaps, &
            & MPI_DOUBLE_PRECISION, root, 90, comm, status, ierr)
       
    else
       
       call mpi_recv(local_map(0:seg_size-1,1:nmaps), seg_size*nmaps, &
            & MPI_DOUBLE_PRECISION, root, 90, comm, status, ierr)
       
    end if
    
  end subroutine distribute_map_slave


  subroutine gather_map_slave
    ! Gathers the map data
    implicit none

    integer(i4b) :: i, seg_size, ierr

    integer(i4b), dimension(MPI_STATUS_SIZE)          :: status

    !Processes 1 and up sends map data
    if(num_segment == 2) then
       
       seg_size = segment(0,1)-segment(0,0)+1

       call mpi_send(local_map(0:2*seg_size-1,1:nmaps), 2*seg_size*nmaps, &
            & MPI_DOUBLE_PRECISION, root, 91, comm, ierr)
       
    else
       
       seg_size = segment(0,1)-segment(0,0)+1
       call mpi_send(local_map(0:seg_size-1,1:nmaps), seg_size*nmaps, &
            & MPI_DOUBLE_PRECISION, root, 91, comm, ierr)
       
    end if

  end subroutine gather_map_slave


  subroutine distribute_w8ring(local_w8ring, w8ring)
    implicit none

    real(dp), dimension(1:,1:), intent(out)           :: local_w8ring
    real(dp), dimension(1:,1:), intent(in), optional  :: w8ring
    
    integer(i4b) :: i, ierr

    integer(i4b), dimension(MPI_STATUS_SIZE)          :: status

    if (myid == root) then

       do i = 1, numprocs-1
          call mpi_send(w8ring(ring_procs(1,i):ring_procs(2,i),1:nmaps), &
               & (ring_procs(2,i)-ring_procs(1,i)+1)*nmaps, MPI_DOUBLE_PRECISION, &
               & i, 92, comm, ierr)
       end do

       local_w8ring = w8ring(ring_procs(1,0):ring_procs(2,0),1:nmaps)

    else

       call mpi_recv(local_w8ring, (last_ring-first_ring+1)*nmaps, &
            & MPI_DOUBLE_PRECISION, root, 92, comm, status, ierr)

    end if
    
  end subroutine distribute_w8ring

  
  ! Finds the rings to compute for the processor with ID "myid"
  subroutine find_rings(myid, num_inc_rings, first_ring, last_ring)
    implicit none

    integer(i4b), intent(in)  :: myid, num_inc_rings
    integer(i4b), intent(out) :: first_ring, last_ring

    integer(i4b) :: i, base_size, remainder, n_iter, ring
    real(dp)     :: integral, x, dx, int_per_proc

    real(dp)     :: alpha = 1.8
    real(dp)     :: beta  = 0.45
    real(dp)     :: gamma = 1.1

    if (numprocs == 1) then
       first_ring = 1
!       last_ring  = num_inc_rings
       last_ring  = 2*nsmax
    else

       ring         = 1
       integral     = 0.d0
!       dx           = 1.d0 / real(num_inc_rings,dp)
       dx           = 1.d0 / real(2*nsmax,dp)
       x            = 0.d0
       int_per_proc = 1.d0 / real(numprocs,dp)

       do i = 0, myid-1
          do while (integral < real(i+1,dp)*int_per_proc)
             ring = ring+1
             x    = x+dx
             if (ring <= nsmax) then
                integral = integral + (alpha*x+beta)*dx
             else
                integral = integral + gamma*dx
             end if
          end do
       end do

       first_ring = ring
       do while (integral < real(myid+1,dp)*int_per_proc)
          ring = ring+1
          x = x+dx
          if (ring <= nsmax) then
             integral = integral + (alpha*x+beta)*dx
          else
             integral = integral + gamma*dx
          end if
       end do
       last_ring = ring-1

    end if

!    if (myid == numprocs-1) last_ring = num_inc_rings
    if (myid == numprocs-1) last_ring = 2*nsmax

!    n_iter    = 2*nsmax
!    base_size = n_iter/numprocs
!    remainder = n_iter - base_size*numprocs 

!    if (myid < remainder) then 
!       first_ring = base_size*myid + myid + 1
!       last_ring  = first_ring + base_size 
!    else 
!       first_ring = base_size*myid + remainder + 1
!       last_ring  = first_ring + base_size - 1
!    end if

  end subroutine find_rings


  !Returns the map segment(s) needed given ring numbers
  subroutine find_pixel_numbers(first_ring, last_ring, num_segment, map_size, segment)
    implicit none

    integer(i4b), intent(in)  :: first_ring, last_ring
    integer(i4b), intent(out) :: num_segment, map_size
    integer(i4b), dimension(0:1,0:1), intent(out) :: segment

    num_segment = 2

    segment(0,0) = ring_start(first_ring) 
    segment(0,1) = ring_start(last_ring+1)-1
    segment(1,0) = (12*nsmax**2) - segment(0,1)-1
    segment(1,1) = (12*nsmax**2) - segment(0,0)-1

    if (segment(1,0) < segment(0,1)) then
       ! End of segment 1 is after the start of segment 2 => overlapping
       ! at equator.
       num_segment  = 1
       segment(0,1) = segment(1,1)
       segment(1,0) = 0
       segment(1,1) = 0
    end if

    if (num_segment == 2) then 
       map_size = segment(0,1) - segment(0,0) + segment(1,1) - segment(1,0) + 2
    else
       map_size = segment(0,1) - segment(0,0) + 1
    end if
    
  end subroutine find_pixel_numbers

  ! Finds the first pixel in ring number i
  function ring_start(i)
    implicit none

    integer(i4b), intent(in) :: i

    integer(i4b) :: ring_start

    if (i <= nsmax) then
       ring_start = 2 * (i**2-i)
    else
       ring_start = 2 * (nsmax**2-nsmax) + 4*nsmax * (i-nsmax)
    endif

  end function ring_start



  ! ***************************************************************
  !                Computational routines
  ! ***************************************************************

  subroutine mpi_map2alm_sc
    !=======================================================================
    !     computes the a(l,m) from a Temperature map for the HEALPIX pixelisation
    !        all from scratch
    !=======================================================================

    integer(I4B) :: l, m, ith, scalem, scalel   ! alm related
    integer(I4B) :: istart_south, istart_north  ! map related
    integer(I4B) :: nrings, npix, nphmx
    real(DP)     :: omega_pix

    integer(I4B)                              :: par_lm
    real(DP)              :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)                      :: ovflow, unflow
    real(DP),     dimension(-1:2)             :: phas_sd
    real(DP),     dimension(:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: mfac
    real(DP),     dimension(:,:), allocatable :: recfac

    integer(I4B)                              :: l_min, kphi1
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:), allocatable :: phas_n, phas_s
    real(DP),     dimension(:),   allocatable :: ring
    integer(I4B)                   :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK-1) :: cth, sth
    logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it

    character(LEN=*), PARAMETER :: code = 'MAP2ALM'
    integer(I4B) :: status

    real(dp) :: t1, t2

    !=======================================================================

    ! Healpix definitions
    nrings = (last_ring-first_ring+1)
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / real(npix, DP)  ! pixel area (identical for all pixels)

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(phas_n(0:nmmax,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_n')

    allocate(phas_s(0:nmmax,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_s')

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax,mfac)

    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    local_alms = 0.0 ! set the whole alm array to zero

    istart_north = 0
    istart_south = map_size

    ! loop on chunks
    do ichunk = 0, nchunks-1

       lchk  = first_ring + ichunk * chunksize 
       uchk  = min(lchk+chunksize - 1, last_ring)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          ! find out which rings are to be analysed
          call select_rings(cth(ithl), zbounds, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       enddo

       !-----------------------------------------------------------------------
       !  computes the integral in phi : phas_m(theta)
       !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
       !-----------------------------------------------------------------------
       phas_n = 0_dpc
       phas_s = 0_dpc

       allocate(ring(0:nphmx-1),stat = status)
       call assert_alloc(status,code,'ring')

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          nphl = nph(ithl)
          ! do Fourier Transform on rings
          if (keep_north(ithl)) then
             ring(0:nphl-1) = local_map(istart_north:istart_north+nphl-1,1) * &
                  & local_w8ring(ith-first_ring+1,1)
             call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_n(0,ithl), kphi0(ithl))
             istart_north = istart_north+nphl
          endif

          if (ith < 2*nsmax .and. keep_south(ithl)) then
             istart_south = istart_south-nphl
             ring(0:nphl-1) = local_map(istart_south:istart_south+nphl-1,1) * &
                  & local_w8ring(ith-first_ring+1,1)
             call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_s(0,ithl), kphi0(ithl))
          endif
       enddo ! loop on ring

       deallocate(ring)

       !-----------------------------------------------------------------------
       !              computes the a_lm by integrating over theta
       !                  lambda_lm(theta) * phas_m(theta)
       !                         for each m and l
       !-----------------------------------------------------------------------

       allocate(recfac(0:1,0:nlmax),stat = status)
       call assert_alloc(status,code,'recfac')
       allocate(dalm(1:2,0:nlmax),stat = status)
       call assert_alloc(status,code,'dalm')

       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)

          ! introduce double precision vector to perform summation over ith for each l
          dalm(1:2, m:nlmax ) = 0.0_dp

          do ithl = 0, uchk-lchk
             l_min = l_min_ylm(m, sth(ithl))
             if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)
                ! determine lam_mm
                call mpi_compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)

                !           ---------- l = m ----------
                par_lm = 1

                php = phas_n(m,ithl) + phas_s(m,ithl) ! sum  (if (l+m) even)
                phm = phas_n(m,ithl) - phas_s(m,ithl) ! diff (if (l+m) odd)
                phas_sd(-1:0) =  (/ real(phm, kind=dp), aimag(phm) /)
                phas_sd( 1:2) =  (/ real(php, kind=dp), aimag(php) /)
                
                if (m >= l_min) then
                   lam_lm = corfac * lam_mm !Actual lam_mm 
                   dalm(1:2, m) = dalm(1:2, m) + lam_lm * phas_sd(par_lm:par_lm+1)
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
                      dalm(1:2, l) = dalm(1:2, l) &
                           &       + lam_lm * phas_sd(par_lm:par_lm+1)
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

                enddo ! loop on l
             endif ! test on cut sky and nlmax
          enddo ! loop on ithl
          do l = m, nlmax
             local_alms(1,l,m) = local_alms(1,l,m) + &
                  & cmplx(dalm(1, l), dalm(2, l), kind=DP) * omega_pix
          enddo
       enddo ! loop on m

       deallocate (recfac,dalm)

    enddo ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate(phas_n,phas_s)
    deallocate(mfac)

    return

  end subroutine mpi_map2alm_sc


  subroutine mpi_map2alm_sc_pre
    !=======================================================================
    !     computes the a(l,m) from a Temperature map for the HEALPIX pixelisation
    !       with precomputed Ylm(theta)
    !=======================================================================

    integer(I8B) :: istart_south, istart_north, npix
    integer(I4B) :: nrings, nphmx
    real(DP)     :: omega_pix

    integer(I8B)                              :: n_lm, i_mm
    integer(I4B)                              :: l, m, ith
    integer(I4B)                              :: par_lm
    real(DP),     dimension(-1:2)             :: phas_sd
    real(DP),     dimension(:,:), allocatable :: dalm
    integer(I4B)                              :: l_min
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:), allocatable :: phas_n, phas_s
    real(DP),     dimension(:),   allocatable :: ring
    real(DP)                                  :: cth
    integer(I4B)                   :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    integer(i4b)                   :: plm_ind
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    real(DP),     dimension(0:SMAXCHK-1) :: sth
    logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it

    character(LEN=*), PARAMETER :: code = 'MAP2ALM'
    integer(I4B) :: i, l_start, status, seed
    real(dp) :: t1, t2

    !=======================================================================

    ! Healpix definitions
    nrings = (last_ring-first_ring+1)
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / real(npix, DP)  ! pixel area (identical for all pixels)
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(phas_n(0:nmmax,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_n')

    allocate(phas_s(0:nmmax,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_s')

    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')
    allocate(dalm(1:2,0:nlmax),stat = status)
    call assert_alloc(status,code,'dalm')

    !     ------------ initiate variables and arrays ----------------
    local_alms = 0.0 ! set the whole alm array to zero
    istart_north = 0
    istart_south = map_size

    ! loop on chunks
    do ichunk = 0, nchunks-1

       call cpu_time(t1)

       lchk = first_ring + ichunk * chunksize 
       uchk = min(lchk+chunksize - 1, last_ring)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth, sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          ! find out which rings are to be analysed
          call select_rings(cth, zbounds, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       enddo

       !-----------------------------------------------------------------------
       !  computes the integral in phi : phas_m(theta)
       !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
       !-----------------------------------------------------------------------
       phas_n = 0_dpc
       phas_s = 0_dpc

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          nphl = nph(ithl)
          ! do Fourier Transform on rings
          if (keep_north(ithl)) then
             ring(0:nphl-1) = local_map(istart_north:istart_north+nphl-1,1) * &
                  & local_w8ring(ith-first_ring+1,1)
             call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_n(0,ithl), kphi0(ithl))
             istart_north = istart_north+nphl
          endif

          if (ith < 2*nsmax .and. keep_south(ithl)) then
             istart_south = istart_south-nphl
             ring(0:nphl-1) = local_map(istart_south:istart_south+nphl-1,1) * &
                  & local_w8ring(ith-first_ring+1,1)
             call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_s(0,ithl), kphi0(ithl))
          endif
       enddo ! loop on ring

       !-----------------------------------------------------------------------
       !              computes the a_lm by integrating over theta
       !                  lambda_lm(theta) * phas_m(theta)
       !                         for each m and l
       !-----------------------------------------------------------------------

       do m = 0, nmmax
 
          ! introduce double precision vector to perform summation over ith for each l
          dalm(1:2, m:nlmax ) = 0.0_dp

          do ithl = 0, uchk-lchk 

             l_min = l_min_ylm(m, sth(ithl))
             if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)
                ith = ithl + lchk
                i_mm = n_lm * (ith-first_ring) + &
                     & ((2_I8B*nlmax+3-m)*m)/2 ! location of Ym,m for ring ith

                !           ---------- l = m ----------
                par_lm = 1

                php = phas_n(m,ithl) + phas_s(m,ithl) ! sum  (if (l+m) even)
                phm = phas_n(m,ithl) - phas_s(m,ithl) ! diff (if (l+m) odd)
                phas_sd(-1:0) =  (/ real(phm, kind=dp), aimag(phm) /)
                phas_sd( 1:2) =  (/ real(php, kind=dp), aimag(php) /)
                
                if (m >= l_min) then
                   dalm(1:2, m) = dalm(1:2, m) + plm(i_mm,1) * phas_sd(par_lm:par_lm+1)
                endif

                !           ---------- l > m ----------
                l_start = max(m+1, l_min)
                if (mod(m+l_start,2) == 0) l_start = l_start+1
                do l = l_start, nlmax, 2
                   dalm(1,l) = dalm(1,l) + plm(i_mm+l-m,1) * phas_sd(-1)
                   dalm(2,l) = dalm(2,l) + plm(i_mm+l-m,1) * phas_sd(0)
                enddo ! loop on l

                l_start = max(m+2, l_min)
                if (mod(m+l_start,2) == 1) l_start = l_start+1
                do l = l_start, nlmax, 2
                   dalm(1,l) = dalm(1,l) + plm(i_mm+l-m,1) * phas_sd(1)
                   dalm(2,l) = dalm(2,l) + plm(i_mm+l-m,1) * phas_sd(2)
                enddo ! loop on l

             endif ! test on cut sky and nlmax
          enddo ! loop on ithl
          do l = m, nlmax
             local_alms(1,l,m) = local_alms(1,l,m) + &
                  & cmplx(dalm(1, l), dalm(2, l), kind=DP) * omega_pix
          enddo
       enddo ! loop on m

       call cpu_time(t2)
!       write(*,*) (lchk+uchk)/2, t2-t1

    enddo ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate(dalm, ring)
    deallocate(phas_n,phas_s)

    return

  end subroutine mpi_map2alm_sc_pre


  subroutine mpi_map2alm_pol
    !=======================================================================
    !     computes the a(l,m) from a polarized map for the HEALPIX pixelisation
    !       all from scratch
    !=======================================================================

    integer(I4B) :: l, m, ith, scalel, scalem, nrings, nphmx
    integer(I8B) :: istart_south, istart_north, npix
    real(DP)     :: omega_pix

    integer(I4B)                              :: par_lm
    real(DP)            :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)                 :: lambda_w, lambda_x, normal_m, lam_lm1m
    real(DP)                 :: fm2, fl, flm1
    real(DP)                 :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP)                 :: ovflow, unflow
    real(DP),     dimension(-3:8)             :: phas_sd
    real(DP),     dimension(-3:6)             :: phas_sdx
    real(DP),     dimension(:,:), allocatable :: recfac, dalm
    real(DP),     dimension(:),   allocatable :: lam_fact, mfac

    integer(i4b)                              :: l_min
    real(DP),     dimension(:),   allocatable :: normal_l
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:,:), allocatable :: phas_ns
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)    :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i, ii
    integer(i4b)    :: l_start
    integer(I8B), dimension(0:SMAXCHK) :: startpix
    integer(I4B), dimension(0:SMAXCHK) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK) :: cth, sth
    real(DP),     dimension(0:SMAXCHK) :: one_on_s2, c_on_s2
    logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it

    character(LEN=*), parameter :: code = 'MAP2ALM'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = (last_ring-first_ring+1)
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / real(npix, DP)

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(normal_l(0:nlmax),stat = status)
    call assert_alloc(status,code,'normal_l')

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(phas_ns(0:nlmax,0:5,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_ns')

    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')
    allocate(recfac(0:1,0:nlmax),stat = status)
    call assert_alloc(status,code,'recfac')
    allocate(dalm(0:5,0:nlmax), stat = status)
    call assert_alloc(status,code,'dalm')
    allocate(lam_fact(0:nlmax),stat = status)
    call assert_alloc(status,code,'lam_fact')

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax,mfac)

    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    local_alms = 0.0 ! set the whole alm array to zero
    istart_north = 0
    istart_south = map_size

    ! generate Polarization normalisation
    call gen_normpol(nlmax, normal_l)

    ! loop on chunks
    do ichunk = 0, nchunks-1

       lchk = first_ring + ichunk * chunksize 
       uchk = min(lchk+chunksize - 1, last_ring)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          one_on_s2(ithl) =    1.0_dp / sth(ithl)**2
            c_on_s2(ithl) = cth(ithl) / sth(ithl)**2
          ! find out which rings are to be analysed
          call select_rings(cth(ithl), zbounds, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       enddo

       !-----------------------------------------------------------------------
       !  computes the integral in phi : phas_m(theta)
       !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
       !-----------------------------------------------------------------------
       phas_ns = 0_dpc

       ! do Fourier Transform on rings
       do ith = lchk, uchk
          ithl = ith - lchk !local index
          nphl = nph(ithl)
          if (keep_north(ithl)) then
             do i=1,3
                ii = 2*(i-1) ! 0,2,4
                ring(0:nphl-1) = local_map(istart_north:istart_north+nphl-1,i) * &
                     & local_w8ring(ith-first_ring+1,i)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
             enddo
             istart_north = istart_north+nphl
          endif

          if (ith < 2*nsmax .and. keep_south(ithl)) then
             istart_south = istart_south-nphl
             do i=1,3
                ii = 2*i - 1 ! 1,3,5
                ring(0:nphl-1) = local_map(istart_south:istart_south+nphl-1,i) * &
                     & local_w8ring(ith-first_ring+1,i)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
             enddo
          endif
       enddo

       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)
          ! generate Ylm relation factor for degree m
          call gen_lamfac(nlmax, m, lam_fact)
          fm2 = real(m * m, kind = DP)
          normal_m = (2.0_dp * m) * (1 - m)

          ! introduce double precision vector to perform summation over ith for each l
          dalm(0:5, m:nlmax) = 0.0_dp

          do ithl = 0, uchk-lchk
             l_min = l_min_ylm(m, sth(ithl)) ! lower l bound of non-negligible Ylm
             if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)
                ! determine lam_mm
                call mpi_compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)

                do i=0,2 ! loop on T, Q, U
                   phm = phas_ns(m, 2*i, ithl) - phas_ns(m, 2*i+1, ithl) ! N-S: (l+m) odd
                   phas_sd(-3+2*i)   =  real(phm, kind=dp) 
                   phas_sd(-3+2*i+1) = aimag(phm)
                   php = phas_ns(m, 2*i, ithl) + phas_ns(m, 2*i+1, ithl) ! N+S: (l+m) even
                   phas_sd( 3+2*i)   =  real(php, kind=dp) 
                   phas_sd( 3+2*i+1) = aimag(php)
                enddo
                phas_sdx(-3: 0) = (/ phas_sd(8), - phas_sd(7), - phas_sd(6),   phas_sd(5) /)
                phas_sdx( 3: 6) = (/ phas_sd(2), - phas_sd(1), - phas_sd(0),   phas_sd(-1) /)

                !           ---------- l = m ----------
                par_lm = 3  ! = 3 * (-1)^(l+m)
                lam_lm = corfac * lam_mm            ! Actual lam_mm

                if (m >= l_min) then
                   lambda_w =  (normal_m * lam_lm)  * ( 0.5_dp - one_on_s2(ithl) )
                   lambda_x =  (normal_m * lam_lm)  *            c_on_s2(ithl)

                   ! alm_T = \int T Ylm
                   ! alm_G = \int ( - Q Wlm - i U Xlm )
                   ! alm_C = \int ( - U Wlm + i Q Xlm )
                   dalm(0:1, m) = dalm(0:1, m) + lam_lm * phas_sd(par_lm:par_lm+1)
                   dalm(2:5, m) = dalm(2:5, m)  &
                        &       - lambda_w * phas_sd(par_lm+2:par_lm+5) &
!                         &       + lambda_x * (/ phas_sd(-par_lm+5), - phas_sd(-par_lm+4), &
!                         &                     - phas_sd(-par_lm+3),   phas_sd(-par_lm+2) /)
                        &       + lambda_x * phas_sdx(par_lm:par_lm+3)
                endif

                !           ---------- l > m ----------
                lam_0 = 0.0_dp
                lam_1 = 1.0_dp
                scalel=0
                cth_ring = cth(ithl)
                lam_2 = cth_ring * lam_1 * recfac(0,m)
                fm_on_s2     =      m * one_on_s2(ithl)
                two_on_s2    = 2.0_dp * one_on_s2(ithl)
                two_cth_ring = 2.0_dp * cth_ring
                b_w          =  c_on_s2(ithl) 
                do l = m+1, nlmax
                   par_lm = - par_lm  ! = 3 * (-1)^(l+m)
                   if (l >= l_min) then
                      lam_lm1m = lam_lm * lam_fact(l)
                      lam_lm = lam_2 * corfac * lam_mm

                      fl = real(l, kind = DP)
                      flm1 = fl - 1.0_dp
                      a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                      a_x =  two_cth_ring * flm1
                      lambda_w =                b_w * lam_lm1m - a_w * lam_lm 
                      lambda_x = fm_on_s2 * (         lam_lm1m - a_x * lam_lm)

                      dalm(0:1, l) = dalm(0:1, l) + lam_lm * phas_sd(par_lm:par_lm+1)
                      dalm(2:5, l) = dalm(2:5, l)  &
                           &       - lambda_w * phas_sd(par_lm+2:par_lm+5) &
!                            &       + lambda_x * (/ phas_sd(-par_lm+5), - phas_sd(-par_lm+4), &
!                            &                     - phas_sd(-par_lm+3),   phas_sd(-par_lm+2) /)
                           &       + lambda_x * phas_sdx(par_lm:par_lm+3)
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

             endif ! test on cut sky
          enddo ! loop on ithl
          do l = m, nlmax
             local_alms(1,l,m) = local_alms(1,l,m) + cmplx(dalm(0,l), dalm(1,l)) * omega_pix
             local_alms(2,l,m) = local_alms(2,l,m) + &
                  & cmplx(dalm(2,l), dalm(3,l)) * normal_l(l) * omega_pix
             local_alms(3,l,m) = local_alms(3,l,m) + &
                  & cmplx(dalm(4,l), dalm(5,l)) * normal_l(l) * omega_pix
          enddo
       enddo ! loop on m

    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate (recfac,dalm,lam_fact,ring)
    deallocate(mfac, normal_l)
    deallocate(phas_ns)

    return

  end subroutine mpi_map2alm_pol


  subroutine mpi_map2alm_pol_pre1
    !=======================================================================
    !     computes the a(l,m) from a polarized map for the HEALPIX pixelisation
    !       with precomputed Ylm_T(theta)
    !=======================================================================

    integer(I4B) :: l, m, ith
    integer(I8B) :: istart_south, istart_north, npix, n_lm, i_mm
    integer(I4B) :: nrings, nphmx
    real(DP)     :: omega_pix

    integer(I4B)                              :: par_lm
    real(DP)            :: lam_lm, cth_ring
    real(DP)                 :: lambda_w, lambda_x, normal_m, lam_lm1m
    real(DP)                 :: fm2, fl, flm1
    real(DP)                 :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP),     dimension(-3:8)             :: phas_sd
    real(DP),     dimension(-3:6)             :: phas_sdx
    real(DP),     dimension(:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: lam_fact

    integer(i4b)                              :: l_min
    real(DP),     dimension(:),   allocatable :: normal_l
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:,:), allocatable :: phas_ns
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl, i, ii
    integer(i4b)        :: l_start
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK-1) :: cth, sth
    real(DP),     dimension(0:SMAXCHK-1) :: one_on_s2, c_on_s2
    logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it
    integer(I4B)                       :: nphl

    character(LEN=*), parameter :: code = 'MAP2ALM'
    integer(I4B) :: status
    logical(lgt) :: start_iter_even, start_iter_odd

    real(dp) :: t1, t2, t3

    !=======================================================================

    ! Healpix definitions
    nrings = (last_ring-first_ring+1)
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / real(npix, DP)
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(normal_l(0:nlmax),stat = status)
    call assert_alloc(status,code,'normal_l')

    allocate(phas_ns(0:nlmax,0:5,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_ns')

    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')
    allocate(dalm(0:5,0:nlmax), stat = status)
    call assert_alloc(status,code,'dalm')
    allocate(lam_fact(0:nlmax),stat = status)
    call assert_alloc(status,code,'lam_fact')

    !     ------------ initiate variables and arrays ----------------

    local_alms = 0.0 ! set the whole alm array to zero
    istart_north = 0
    istart_south = map_size

    ! generate Polarization normalisation
    call gen_normpol(nlmax, normal_l)

    ! loop on chunks
    do ichunk = 0, nchunks-1

       lchk = first_ring + ichunk * chunksize 
       uchk = min(lchk+chunksize - 1, last_ring)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          one_on_s2(ithl) =    1.0_dp / sth(ithl)**2
            c_on_s2(ithl) = cth(ithl) / sth(ithl)**2
          ! find out which rings are to be analysed
          call select_rings(cth(ithl), zbounds, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       enddo

       !-----------------------------------------------------------------------
       !  computes the integral in phi : phas_m(theta)
       !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
       !-----------------------------------------------------------------------
       phas_ns = 0_dpc

       ! do Fourier Transform on rings
       do ith = lchk, uchk
          ithl = ith - lchk !local index
          nphl = nph(ithl)
          if (keep_north(ithl)) then
             do i=1,3
                ii = 2*(i-1) ! 0,2,4
                ring(0:nphl-1) = local_map(istart_north:istart_north+nphl-1,i) * &
                     & local_w8ring(ith-first_ring+1,i)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
             enddo
             istart_north = istart_north+nphl
          endif

          if (ith < 2*nsmax .and. keep_south(ithl)) then
             istart_south = istart_south-nphl
             do i=1,3
                ii = 2*i - 1 ! 1,3,5
                ring(0:nphl-1) = local_map(istart_south:istart_south+nphl-1,i) * &
                     & local_w8ring(ith-first_ring+1,i)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
             enddo
          endif
       enddo

       do m = 0, nmmax
          ! generate Ylm relation factor for degree m
          call gen_lamfac(nlmax, m, lam_fact)
          fm2 = real(m * m, kind = DP)
          normal_m = (2.0_dp * m) * (1 - m)

          ! introduce double precision vector to perform summation over ith for each l
          dalm(0:5, m:nlmax) = 0.0_dp

          do ithl = 0, uchk-lchk
             l_min = l_min_ylm(m, sth(ithl)) ! lower l bound of non-negligible Ylm
             if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)
                ith = ithl + lchk
                i_mm = n_lm * (ith-first_ring) + ((2_I8B*nlmax + 3 - m)*m)/2 ! location of Ym,m for ring ith

                do i=0,2 ! loop on T, Q, U
                   phm = phas_ns(m, 2*i, ithl) - phas_ns(m, 2*i+1, ithl) ! N-S: (l+m) odd
                   phas_sd(-3+2*i)   =  real(phm, kind=dp) 
                   phas_sd(-3+2*i+1) = aimag(phm)
                   php = phas_ns(m, 2*i, ithl) + phas_ns(m, 2*i+1, ithl) ! N+S: (l+m) even
                   phas_sd( 3+2*i)   =  real(php, kind=dp) 
                   phas_sd( 3+2*i+1) = aimag(php)
                enddo
                phas_sdx(-3: 0) = (/ phas_sd(8), - phas_sd(7), - phas_sd(6),   phas_sd(5) /)
                phas_sdx( 3: 6) = (/ phas_sd(2), - phas_sd(1), - phas_sd(0),   phas_sd(-1) /)

                !           ---------- l = m ----------
                par_lm = 3  ! = 3 * (-1)^(l+m)
                lam_lm = plm(i_mm,1)

                if (m >= l_min) then
                   lambda_w =  (normal_m * lam_lm)  * ( 0.5_dp - one_on_s2(ithl) )
                   lambda_x =  (normal_m * lam_lm)  *            c_on_s2(ithl)

                   ! alm_T = \int T Ylm
                   ! alm_G = \int ( - Q Wlm - i U Xlm )
                   ! alm_C = \int ( - U Wlm + i Q Xlm )
                   dalm(0:1, m) = dalm(0:1, m) + lam_lm * phas_sd(par_lm:par_lm+1)
                   dalm(2:5, m) = dalm(2:5, m)  &
                        &       - lambda_w * phas_sd(par_lm+2:par_lm+5) &
!                         &       + lambda_x * (/ phas_sd(-par_lm+5), - phas_sd(-par_lm+4), &
!                         &                     - phas_sd(-par_lm+3),   phas_sd(-par_lm+2) /)
                        &       + lambda_x * phas_sdx(par_lm:par_lm+3)
                endif

                !           ---------- l > m ----------
                cth_ring = cth(ithl)
                fm_on_s2     =      m * one_on_s2(ithl)
                two_on_s2    = 2.0_dp * one_on_s2(ithl)
                two_cth_ring = 2.0_dp * cth_ring
                b_w          =  c_on_s2(ithl) 

                l_start = max(m+1, l_min)
                if (mod(m+l_start,2) == 0) l_start = l_start+1
                do l = l_start, nlmax, 2

                   lam_lm1m = plm(i_mm+l-m-1,1) * lam_fact(l)
                   lam_lm   = plm(i_mm+l-m,1)

                   fl = real(l, kind = DP)
                   flm1 = fl - 1.0_dp
                   a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                   a_x =  two_cth_ring * flm1
                   lambda_w =                b_w * lam_lm1m - a_w * lam_lm 
                   lambda_x = fm_on_s2 * (         lam_lm1m - a_x * lam_lm)
                   
                   dalm(0,l) = dalm(0,l) + lam_lm * phas_sd(-3)
                   dalm(1,l) = dalm(1,l) + lam_lm * phas_sd(-2)

                   dalm(2,l) = dalm(2,l) - lambda_w * phas_sd(-1) + lambda_x * phas_sdx(-3)
                   dalm(3,l) = dalm(3,l) - lambda_w * phas_sd(0) + lambda_x * phas_sdx(-2)
                   dalm(4,l) = dalm(4,l) - lambda_w * phas_sd(1) + lambda_x * phas_sdx(-1)
                   dalm(5,l) = dalm(5,l) - lambda_w * phas_sd(2) + lambda_x * phas_sdx(0)
                enddo ! loop on l

                l_start = max(m+2, l_min)
                if (mod(m+l_start,2) == 1) l_start = l_start+1
                do l = l_start, nlmax, 2

                   lam_lm1m = plm(i_mm+l-m-1,1) * lam_fact(l)
                   lam_lm   = plm(i_mm+l-m,1)

                   fl = real(l, kind = DP)
                   flm1 = fl - 1.0_dp
                   a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                   a_x =  two_cth_ring * flm1
                   lambda_w =                b_w * lam_lm1m - a_w * lam_lm 
                   lambda_x = fm_on_s2 * (         lam_lm1m - a_x * lam_lm)

                   dalm(0,l) = dalm(0,l) + lam_lm * phas_sd(3)
                   dalm(1,l) = dalm(1,l) + lam_lm * phas_sd(4)
                   dalm(2,l) = dalm(2,l) - lambda_w * phas_sd(5) + lambda_x * phas_sdx(3)
                   dalm(3,l) = dalm(3,l) - lambda_w * phas_sd(6) + lambda_x * phas_sdx(4)
                   dalm(4,l) = dalm(4,l) - lambda_w * phas_sd(7) + lambda_x * phas_sdx(5)
                   dalm(5,l) = dalm(5,l) - lambda_w * phas_sd(8) + lambda_x * phas_sdx(6)
                enddo ! loop on l

             endif ! test on cut sky
          enddo ! loop on ithl

          do l = m, nlmax
             local_alms(1,l,m) = local_alms(1,l,m) + &
                  & cmplx(dalm(0,l), dalm(1,l)) * omega_pix
             local_alms(2,l,m) = local_alms(2,l,m) + &
                  & cmplx(dalm(2,l), dalm(3,l)) * normal_l(l) * omega_pix
             local_alms(3,l,m) = local_alms(3,l,m) + &
                  & cmplx(dalm(4,l), dalm(5,l)) * normal_l(l) * omega_pix
          enddo

       enddo ! loop on m

    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate (dalm,lam_fact,ring)
    deallocate(normal_l, phas_ns)

    return

  end subroutine mpi_map2alm_pol_pre1


  subroutine mpi_map2alm_pol_pre2

    !=======================================================================
    !     computes the a(l,m) from a polarized map for the HEALPIX pixelisation
    !       with precomputed Ylm_T(theta) an dYlm_P(theta)
    !=======================================================================

    integer(I4B) :: l, m, ith
    integer(I8B) :: istart_south, istart_north, npix, n_lm, i_mm, i_up
    integer(I4B) :: nrings, nphmx
    real(DP)     :: omega_pix

    integer(I4B)                              :: par_lm
    real(DP),     dimension(-3:8)             :: phas_sd
    real(DP),     dimension(-3:6)             :: phas_sdx
    real(DP),     dimension(:,:), allocatable :: dalm, plm_sub

    integer(i4b)                              :: l_min
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:,:), allocatable :: phas_ns
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)    :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i, ii
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    real(DP)                           :: cth
    real(DP),     dimension(0:SMAXCHK-1) :: sth
    logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it

    character(LEN=*), parameter :: code = 'MAP2ALM'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = (last_ring-first_ring+1)
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / real(npix, DP)
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(phas_ns(0:nlmax,0:5,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_ns')

    allocate(dalm(0:5,0:nlmax), plm_sub(1:3,0:nlmax), stat = status)
    call assert_alloc(status,code,'dalm & plm_sub')
    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')

    !     ------------ initiate variables and arrays ----------------

    local_alms = 0.0 ! set the whole alm array to zero
    istart_north = 0
    istart_south = map_size

    ! loop on chunks
    do ichunk = 0, nchunks-1

       lchk = first_ring + ichunk * chunksize 
       uchk = min(lchk+chunksize - 1, last_ring)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth, sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          ! find out which rings are to be analysed
          call select_rings(cth, zbounds, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       enddo

       !-----------------------------------------------------------------------
       !  computes the integral in phi : phas_m(theta)
       !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
       !-----------------------------------------------------------------------
       phas_ns = 0_dpc

       ! do Fourier Transform on rings
       do ith = lchk, uchk
          ithl = ith - lchk !local index
          nphl = nph(ithl)
          if (keep_north(ithl)) then
             do i=1,3
                ii = 2*(i-1) ! 0,2,4
                ring(0:nphl-1) = local_map(istart_north:istart_north+nphl-1,i) * &
                     & local_w8ring(ith-first_ring+1,i)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
             enddo
             istart_north = istart_north + nphl
          endif

          if (ith < 2*nsmax .and. keep_south(ithl)) then
             istart_south = istart_south - nphl
             do i=1,3
                ii = 2*i - 1 ! 1,3,5
                ring(0:nphl-1) = local_map(istart_south:istart_south+nphl-1,i) * &
                     & local_w8ring(ith-first_ring+1,i)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
             enddo
          endif
       enddo

       do m = 0, nmmax 

          ! introduce double precision vector to perform summation over ith for each l
          dalm(0:5, m:nlmax) = 0.0_dp

          do ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl)) ! lower l bound of non-negligible Ylm
             if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)
                ith = ithl + lchk
                i_mm = n_lm * (ith-first_ring) + ((2_I8B*nlmax + 3 - m)*m)/2 ! location of Ym,m for ring ith
                i_up = i_mm + nlmax - m
                plm_sub(1,m:nlmax) = plm(i_mm:i_up,1)
                plm_sub(2,m:nlmax) = plm(i_mm:i_up,2)
                plm_sub(3,m:nlmax) = plm(i_mm:i_up,3)

                do i=0,2 ! loop on T, Q, U
                   phm = phas_ns(m, 2*i, ithl) - phas_ns(m, 2*i+1, ithl) ! N-S: (l+m) odd
                   phas_sd(-3+2*i)   =  real(phm, kind=dp) 
                   phas_sd(-3+2*i+1) = aimag(phm)
                   php = phas_ns(m, 2*i, ithl) + phas_ns(m, 2*i+1, ithl) ! N+S: (l+m) even
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
!                         &       + plm_sub(3,m) * (/ phas_sd(-par_lm+5), - phas_sd(-par_lm+4), &
!                         &                     - phas_sd(-par_lm+3),   phas_sd(-par_lm+2) /)
                        &       + plm_sub(3,m) * phas_sdx(par_lm:par_lm+3)
                endif

                !           ---------- l > m ----------
                do l = m+1, nlmax
                   par_lm = - par_lm  ! = 3 * (-1)^(l+m)
                   if (l >= l_min) then
                      dalm(0:1, l) = dalm(0:1, l) + plm_sub(1,l) * phas_sd(par_lm:par_lm+1)
                      dalm(2:5, l) = dalm(2:5, l)  &
                           &       - plm_sub(2,l) * phas_sd(par_lm+2:par_lm+5) &
!                            &       + plm_sub(3,l) * (/ phas_sd(-par_lm+5), - phas_sd(-par_lm+4), &
!                            &                     - phas_sd(-par_lm+3),   phas_sd(-par_lm+2) /)
                        &       + plm_sub(3,l) * phas_sdx(par_lm:par_lm+3)
                   endif
                   
                enddo ! loop on l
             endif ! test on cut sky
          enddo ! loop on ithl
          do l = m, nlmax
             local_alms(1,l,m) = local_alms(1,l,m) + cmplx(dalm(0,l), dalm(1,l)) * omega_pix
             local_alms(2,l,m) = local_alms(2,l,m) + cmplx(dalm(2,l), dalm(3,l)) * omega_pix
             local_alms(3,l,m) = local_alms(3,l,m) + cmplx(dalm(4,l), dalm(5,l)) * omega_pix
          enddo
       enddo ! loop on m
    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate(phas_ns)
    deallocate (ring, dalm, plm_sub)

    return

  end subroutine mpi_map2alm_pol_pre2





  subroutine mpi_alm2map_sc
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
    integer(I4B) :: nrings, nphmx

    integer(I4B)                              :: par_lm
    real(DP)              :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)                                  :: ovflow, unflow
    complex(DPC), dimension(-1:1)             :: b_ns
    real(DP),     dimension(:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: mfac
    real(DP),     dimension(:,:), allocatable :: recfac

    integer(i4b)                                :: ll, l_min
    complex(DPC), dimension(:,:), allocatable :: b_north, b_south
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)             :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK-1) :: cth, sth

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status

    real(dp) :: t1, t2
    !=======================================================================

    ! Healpix definitions
    nrings = (last_ring-first_ring+1)
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(b_north(0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_north')

    allocate(b_south(0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_south')

    allocate(recfac(0:1,0:nlmax),stat = status)
    call assert_alloc(status,code,'recfac')
    allocate(dalm(0:1,0:nlmax), stat = status)
    call assert_alloc(status,code,'dalm')
    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax,mfac)

    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    local_map = 0.d0 ! set the whole map to zero

    istart_north = 0
    istart_south = map_size

    ! loop on chunks
    do ichunk = 0, nchunks-1

       lchk = first_ring + ichunk * chunksize 
       uchk = min(lchk+chunksize - 1, last_ring)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
       enddo
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)


       b_north(:,:) = 0_dpc ! pad with zeros
       b_south(:,:) = 0_dpc

       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             dalm(0,ll) =  real(local_alms(1,ll,m),kind=dp)
             dalm(1,ll) = aimag(local_alms(1,ll,m))
          enddo
          do ithl = 0, uchk-lchk
             l_min = l_min_ylm(m, sth(ithl))
             if (nlmax >= l_min) then ! skip calculations when Ylm too small
                ! determine lam_mm
                call mpi_compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)
               
                !           ---------- l = m ----------
                par_lm = 1  ! = (-1)^(l+m)

                if (m >= l_min) then
                   lam_lm = corfac * lam_mm
                   b_ns( 1) = cmplx(lam_lm * dalm(0,m), lam_lm * dalm(1,m), kind=DP)
                   b_ns(-1) = 0.0_dpc
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
                      b_ns(par_lm)  = b_ns(par_lm) &
                           &        + cmplx(lam_lm * dalm(0,l), lam_lm * dalm(1,l), kind=DP)  ! increment even or odd
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

                b_north(m,ithl) = b_ns(1) + b_ns(-1)
                b_south(m,ithl) = b_ns(1) - b_ns(-1)
             endif ! test on nlmax
          enddo ! loop on rings (ithl)
       enddo ! loop on m

       do ithl = 0, uchk-lchk
          nphl = nph(ithl)
          ith  = ithl + lchk
          !        ---------------------------------------------------------------
          !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
          !        ---------------------------------------------------------------
          call ring_synthesis(nsmax,nlmax,nmmax,b_north(0,ithl),nphl,ring,kphi0(ithl))   ! north hemisph. + equator
          local_map(istart_north:istart_north+nphl-1,1) = ring(0:nphl-1)
          istart_north = istart_north+nphl

          if (ith < 2*nsmax) then
             istart_south = istart_south-nphl
             call ring_synthesis(nsmax,nlmax,nmmax,b_south(0,ithl),nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
             local_map(istart_south:istart_south+nphl-1,1) = ring(0:nphl-1)
          endif
       enddo ! loop on ithl

    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate(ring,recfac,dalm)
    deallocate(mfac)
    deallocate(b_north, b_south)

    return

  end subroutine mpi_alm2map_sc


  subroutine mpi_alm2map_sc_pre
    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !      for the Temperature field, with precomputed Ylm
    !     map(theta,phi) = sum_l_m a_lm Y_lm(theta,phi)
    !                    = sum_m {e^(i*m*phi) sum_l a_lm*lambda_lm(theta)}
    !
    !     where Y_lm(theta,phi) = lambda(theta) * e^(i*m*phi)
    !
    !     * the recurrence of Ylm is the standard one (cf Num Rec)
    !     * the sum over m is done by FFT
    !
    !=======================================================================

    integer(I4B) :: l, m, ith
    integer(I8B) :: istart_south, istart_north, npix
    integer(I4B) :: nrings, nphmx

    integer(I8B) :: n_lm, i_mm
    integer(I4B)                              :: par_lm
    complex(DPC), dimension(-1:1)             :: b_ns
    real(DP),     dimension(:,:), allocatable :: dalm
    integer(i4b)                              :: ll, l_min
    real(DP)                                  :: cth
    complex(DPC), dimension(:,:), allocatable :: b_north, b_south
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)             :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK-1) :: sth

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status, l_start

    real(dp) :: t1, t2

    !=======================================================================

    ! Healpix definitions
    nrings = (last_ring-first_ring+1)
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(b_north(0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_north')

    allocate(b_south(0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_south')

    allocate(dalm(0:1,0:nlmax), stat = status)
    call assert_alloc(status,code,'dalm')
    allocate(ring(0:nphmx-1),stat = status)
    call assert_alloc(status,code,'ring')

    !     ------------ initiate variables and arrays ----------------

    local_map = 0.d0 ! set the whole map to zero
    istart_north = 0
    istart_south = map_size

    ! loop on chunks
    do ichunk = 0, nchunks-1

       call cpu_time(t1)

       lchk = first_ring + ichunk * chunksize 
       uchk = min(lchk+chunksize - 1, last_ring)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth, sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
       enddo

       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------

       b_north(:,:) = 0_dpc ! pad with zeros
       b_south(:,:) = 0_dpc

       do m = 0, nmmax

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             dalm(0,ll) =  real(local_alms(1,ll,m),kind=dp)
             dalm(1,ll) = aimag(local_alms(1,ll,m))
          enddo
          do ithl = 0, uchk-lchk
             l_min = l_min_ylm(m, sth(ithl))
             if (nlmax >= l_min) then ! skip calculations when Ylm too small
                ith = ithl + lchk
                i_mm = n_lm * (ith-first_ring) + ((2_I8B*nlmax + 3 - m)*m)/2 ! location of Ym,m for ring ith

                !           ---------- l = m ----------
                par_lm = 1  ! = (-1)^(l+m)

                if (m >= l_min) then
                   b_ns( 1) = cmplx(plm(i_mm,1) * dalm(0,m), plm(i_mm,1) * dalm(1,m), kind=DP)
                   b_ns(-1) = 0.0_dpc
                else
                   b_ns = 0.0_dpc
                endif

                !           ---------- l > m ----------
                l_start = max(m+1, l_min)
                if (mod(m+l_start,2) == 0) l_start = l_start+1
                do l = l_start, nlmax, 2
                   b_ns(-1) = b_ns(-1) + &
                        & cmplx(plm(i_mm+l-m,1) * dalm(0,l), plm(i_mm+l-m,1)*dalm(1,l), kind=DP)
                enddo ! loop on l

                l_start = max(m+2, l_min)
                if (mod(m+l_start,2) == 1) l_start = l_start+1
                do l = l_start, nlmax, 2
                   b_ns(1)  = b_ns(1) + &
                        & cmplx(plm(i_mm+l-m,1) * dalm(0,l), plm(i_mm+l-m,1)*dalm(1,l),kind=DP) 
                enddo ! loop on l

                b_north(m,ithl) = b_ns(1) + b_ns(-1)
                b_south(m,ithl) = b_ns(1) - b_ns(-1)
             endif ! test on nlmax
          enddo ! loop on rings (ithl)
       enddo ! loop on m

       do ithl = 0, uchk-lchk
          nphl = nph(ithl)
          ith  = ithl + lchk
          !        ---------------------------------------------------------------
          !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
          !        ---------------------------------------------------------------
          call ring_synthesis(nsmax,nlmax,nmmax,b_north(0,ithl),nphl,ring,kphi0(ithl))   ! north hemisph. + equator
          local_map(istart_north:istart_north+nphl-1,1) = ring(0:nphl-1)
          istart_north = istart_north+nphl

          if (ith < 2*nsmax) then
             call ring_synthesis(nsmax,nlmax,nmmax,b_south(0,ithl),nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
             istart_south = istart_south-nphl
             local_map(istart_south:istart_south+nphl-1,1) = ring(0:nphl-1)
          endif
       enddo ! loop on ithl

       call cpu_time(t2)
!       write(*,*) (lchk+uchk)/2, t2-t1

    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate (ring,dalm)
    deallocate(b_north, b_south)

    return

  end subroutine mpi_alm2map_sc_pre



  subroutine mpi_alm2map_pol
    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !      for the Temperature+Polarization field
    !=======================================================================

    integer(I4B) :: l, m, ith, scalel, scalem ! alm related
    integer(I4B) :: istart_south, istart_north   ! map related
    integer(I4B) :: nrings, npix, nphmx

    integer(I4B)                              :: par_lm
    real(DP)            :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)                 :: lambda_w, lambda_x, normal_m, lam_lm1m
    real(DP)                 :: fm2, fl, flm1
    real(DP)                 :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP)                 :: ovflow, unflow
    real(DP),     dimension(-3:8)             :: b_ns
    real(DP),     dimension(:,:), allocatable :: recfac, dalm
    real(DP),     dimension(:),   allocatable :: lam_fact
    real(DP),     dimension(:), allocatable :: mfac

    integer(i4b)                                :: ll, l_min
    real(DP),     dimension(:),     allocatable :: normal_l
    complex(DPC), dimension(:,:,:), allocatable :: b_TQU
    complex(DPC), dimension(:),     allocatable :: bsub
    real(DP),     dimension(:),     allocatable :: ring
    integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i
    integer(I8B), dimension(0:SMAXCHK) :: startpix
    integer(I4B), dimension(0:SMAXCHK) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK) :: cth, sth
    real(DP),     dimension(0:SMAXCHK) :: one_on_s2, c_on_s2

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = (last_ring-first_ring+1)
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(normal_l(0:nlmax),stat = status)
    call assert_alloc(status,code,'normal_l')

    allocate(b_TQU(1:6, 0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_TQU')

    allocate(recfac(0:1,0:nlmax), dalm(0:5,0:nlmax), lam_fact(0:nlmax),stat = status)
    call assert_alloc(status,code,'recfac, dalm & lam_fact')
    allocate(ring(0:nphmx-1),bsub(0:nlmax),stat = status)
    call assert_alloc(status,code,'ring & bsub')

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax,mfac)

    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    local_map = 0.d0 ! set the whole map to zero
    istart_north = 0
    istart_south = map_size

    ! generate Polarization normalisation
    call gen_normpol(nlmax, normal_l)

    do ichunk = 0, nchunks-1

       lchk = first_ring + ichunk * chunksize 
       uchk = min(lchk+chunksize - 1, last_ring)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          one_on_s2(ithl) =    1.0_dp / sth(ithl)**2
            c_on_s2(ithl) = cth(ithl) / sth(ithl)**2
       enddo
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)


       b_TQU = 0_dpc ! pad with zeros

       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)
          ! generate Ylm relation factor for degree m
          call gen_lamfac(nlmax, m, lam_fact)
          fm2 = real(m * m, kind = DP)
          normal_m = (2.0_dp * m) * (1 - m)

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             dalm(0,ll) =  real(local_alms(1,ll,m),kind=dp) ! T, real
             dalm(1,ll) = aimag(local_alms(1,ll,m))         ! T, imag
             dalm(2,ll) =  real(local_alms(2,ll,m),kind=dp)*normal_l(ll) ! G, real
             dalm(3,ll) = aimag(local_alms(2,ll,m))        *normal_l(ll)
             dalm(4,ll) =  real(local_alms(3,ll,m),kind=dp)*normal_l(ll) ! C, real
             dalm(5,ll) = aimag(local_alms(3,ll,m))        *normal_l(ll)
          enddo
          do ithl = 0, uchk-lchk
             l_min =  l_min_ylm(m, sth(ithl)) ! lowest l of non-negligible Ylm
             if (nlmax >= l_min) then
                ! determine lam_mm
                call mpi_compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)

                !           ---------- l = m ----------
                par_lm = 3  ! = 3 * (-1)^(l+m)
                lam_lm = corfac * lam_mm !Actual lam_mm 

                if (m >= l_min) then ! skip Ymm if too small
                   lambda_w =  (normal_m * lam_lm)  * ( 0.5_dp - one_on_s2(ithl) )
                   lambda_x =  (normal_m * lam_lm)  *            c_on_s2(ithl)

                   !      T =   alm_T * Ylm
                   !      Q = - alm_E * Wlm - i alm_B * Xlm
                   !      U = i alm_E * Xlm -   alm_B * Wlm
                   b_ns(-3:-2) = 0.0_dp               ! T odd
                   b_ns(-1: 0) =   lambda_x * (/   dalm(5,m), - dalm(4,m) /) ! Q odd
                   b_ns( 1: 2) =   lambda_x * (/ - dalm(3,m),   dalm(2,m) /) ! U odd
!                    b_ns(-1) =   lambda_x * dalm(5,m) ! Q odd
!                    b_ns( 0) =  -lambda_x * dalm(4,m) ! Q odd
!                    b_ns( 1) =  -lambda_x * dalm(3,m) ! U odd
!                    b_ns( 2) =   lambda_x * dalm(2,m) ! U odd

                   b_ns( 3: 4) =   lam_lm   * dalm(0:1,m) ! T even
                   b_ns( 5: 6) = - lambda_w * dalm(2:3,m) ! Q even
                   b_ns( 7: 8) = - lambda_w * dalm(4:5,m) ! U even
                else
                   b_ns = 0.0_dp
                endif

                !           ---------- l > m ----------
                lam_0 = 0.0_dp
                lam_1 = 1.0_dp
                scalel=0
                cth_ring = cth(ithl)
                lam_2 = cth_ring * lam_1 * recfac(0,m)
                fm_on_s2     =      m * one_on_s2(ithl)
                two_on_s2    = 2.0_dp * one_on_s2(ithl)
                two_cth_ring = 2.0_dp * cth_ring
                b_w          =  c_on_s2(ithl) 
                do l = m+1, nlmax
                   par_lm = - par_lm  ! = 3 * (-1)^(l+m)
                   if (l >= l_min) then
                      lam_lm1m = lam_lm * lam_fact(l)
                      lam_lm = lam_2 * corfac * lam_mm

                      fl = real(l, kind = DP)
                      flm1 = fl - 1.0_dp
                      a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                      a_x =  two_cth_ring * flm1
                      lambda_w =                b_w * lam_lm1m - a_w * lam_lm 
                      lambda_x = fm_on_s2 * (         lam_lm1m - a_x * lam_lm)

                      b_ns(par_lm:  par_lm+1) = b_ns(par_lm:  par_lm+1) + lam_lm   * dalm(0:1,l) ! T even or odd
                      b_ns(par_lm+2:par_lm+5) = b_ns(par_lm+2:par_lm+5) - lambda_w * dalm(2:5,l) ! Q, U  even or odd
                      b_ns(2-par_lm) = b_ns(2-par_lm) + lambda_x * dalm(5,l) ! Q odd (or even)
                      b_ns(3-par_lm) = b_ns(3-par_lm) - lambda_x * dalm(4,l)
                      b_ns(4-par_lm) = b_ns(4-par_lm) - lambda_x * dalm(3,l) ! U odd (or even)
                      b_ns(5-par_lm) = b_ns(5-par_lm) + lambda_x * dalm(2,l)
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

                b_TQU(1,m,ithl) = cmplx(b_ns(3) + b_ns(-3), b_ns(4) + b_ns(-2), kind = DP) ! T north
                b_TQU(4,m,ithl) = cmplx(b_ns(3) - b_ns(-3), b_ns(4) - b_ns(-2), kind = DP) ! T south
                b_TQU(2,m,ithl) = cmplx(b_ns(5) + b_ns(-1), b_ns(6) + b_ns( 0), kind = DP) ! Q n
                b_TQU(5,m,ithl) = cmplx(b_ns(5) - b_ns(-1), b_ns(6) - b_ns( 0), kind = DP) ! Q s
                b_TQU(3,m,ithl) = cmplx(b_ns(7) + b_ns( 1), b_ns(8) + b_ns( 2), kind = DP) ! U n
                b_TQU(6,m,ithl) = cmplx(b_ns(7) - b_ns( 1), b_ns(8) - b_ns( 2), kind = DP) ! U s
             endif ! test on nlmax
          enddo ! loop on rings (ithl)
       enddo ! loop on m


       do ithl = 0, uchk-lchk
          nphl = nph(ithl)
          ith  = ithl + lchk
          !        ---------------------------------------------------------------
          !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
          !        ---------------------------------------------------------------
          do i=1,3
             bsub = b_TQU(i,:,ithl)
             call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl))   ! north hemisph. + equator
             local_map(istart_north:istart_north+nphl-1, i) = ring(0:nphl-1)
          enddo
          istart_north = istart_north+nphl

          if (ith < 2*nsmax) then
             istart_south = istart_south-nphl
             do i=1,3
                bsub = b_TQU(3+i,:,ithl)
                call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
                local_map(istart_south:istart_south+nphl-1, i) = ring(0:nphl-1)
             enddo
          endif
       enddo ! loop on ithl
    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate(recfac, dalm, lam_fact, ring, bsub)
    deallocate(mfac, normal_l)
    deallocate(b_TQU)

    return

  end subroutine mpi_alm2map_pol


  subroutine mpi_alm2map_pol_pre1
    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !      for the Temperature+Polarization field
    !=======================================================================

    integer(I4B) :: l, m, ith
    integer(I8B) :: istart_south, istart_north, npix, n_lm, i_mm
    integer(I4B) :: nrings, nphmx

    integer(I4B)                              :: par_lm
    real(DP)            :: lam_lm, cth_ring, lambda_w, lambda_x, normal_m, lam_lm1m
    real(DP)                 :: fm2, fl, flm1
    real(DP)                 :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP),     dimension(-3:8)             :: b_ns
    real(DP),     dimension(:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: lam_fact

    integer(i4b)                                :: ll, l_min
    real(DP),     dimension(:),     allocatable :: normal_l
    complex(DPC), dimension(:,:,:), allocatable :: b_TQU
    complex(DPC), dimension(:),     allocatable :: bsub
    real(DP),     dimension(:),     allocatable :: ring
    integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i
    integer(i4b)        :: l_start
    integer(I8B), dimension(0:SMAXCHK) :: startpix
    integer(I4B), dimension(0:SMAXCHK) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK) :: cth, sth
    real(DP),     dimension(0:SMAXCHK) :: one_on_s2, c_on_s2

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = (last_ring-first_ring+1)
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(normal_l(0:nlmax),stat = status)
    call assert_alloc(status,code,'normal_l')

    allocate(b_TQU(1:6, 0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_TQU')

    allocate(dalm(0:5,0:nlmax), lam_fact(0:nlmax),stat = status)
    call assert_alloc(status,code,'dalm & lam_fact')
    allocate(ring(0:nphmx-1),bsub(0:nlmax),stat = status)
    call assert_alloc(status,code,'ring & bsub')

    !     ------------ initiate variables and arrays ----------------

    local_map = 0.d0 ! set the whole map to zero
    istart_north = 0
    istart_south = map_size

    ! generate Polarization normalisation
    call gen_normpol(nlmax, normal_l)

    ! loop on chunks
    do ichunk = 0, nchunks-1

       lchk = first_ring + ichunk * chunksize 
       uchk = min(lchk+chunksize - 1, last_ring)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          one_on_s2(ithl) =    1.0_dp / sth(ithl)**2
            c_on_s2(ithl) = cth(ithl) / sth(ithl)**2
       enddo
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)


       b_TQU = 0_dpc ! pad with zeros

       do m = 0, nmmax
          ! generate Ylm relation factor for degree m
          call gen_lamfac(nlmax, m, lam_fact)
          fm2 = real(m * m, kind = DP)
          normal_m = (2.0_dp * m) * (1 - m)

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             dalm(0,ll) =  real(local_alms(1,ll,m),kind=dp) ! T, real
             dalm(1,ll) = aimag(local_alms(1,ll,m))         ! T, imag
             dalm(2,ll) =  real(local_alms(2,ll,m),kind=dp)*normal_l(ll) ! G, real
             dalm(3,ll) = aimag(local_alms(2,ll,m))        *normal_l(ll)
             dalm(4,ll) =  real(local_alms(3,ll,m),kind=dp)*normal_l(ll) ! C, real
             dalm(5,ll) = aimag(local_alms(3,ll,m))        *normal_l(ll)
          enddo
          do ithl = 0, uchk-lchk
             l_min =  l_min_ylm(m, sth(ithl)) ! lowest l of non-negligible Ylm
             if (nlmax >= l_min) then
                ith = ithl + lchk
                i_mm = n_lm * (ith-first_ring) + ((2_I8B*nlmax + 3 - m)*m)/2 ! location of Ym,m for ring ith

                !           ---------- l = m ----------
                par_lm = 3  ! = 3 * (-1)^(l+m)
                lam_lm = plm(i_mm,1)

                if (m >= l_min) then ! skip Ymm if too small
                   lambda_w =  (normal_m * lam_lm)  * ( 0.5_dp - one_on_s2(ithl) )
                   lambda_x =  (normal_m * lam_lm)  *            c_on_s2(ithl)

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
                fm_on_s2     =      m * one_on_s2(ithl)
                two_on_s2    = 2.0_dp * one_on_s2(ithl)
                two_cth_ring = 2.0_dp * cth_ring
                b_w          =  c_on_s2(ithl) 

                l_start = max(m+1, l_min)
                if (mod(m+l_start,2) == 0) l_start = l_start+1
                do l = l_start, nlmax, 2

                      lam_lm1m = plm(i_mm+l-m-1,1) * lam_fact(l)
                      lam_lm   = plm(i_mm+l-m,1)

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

                   lam_lm1m = plm(i_mm+l-m-1,1) * lam_fact(l)
                   lam_lm   = plm(i_mm+l-m,1)
                   
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

                b_TQU(1,m,ithl) = cmplx(b_ns(3) + b_ns(-3), b_ns(4) + b_ns(-2), kind = DP) ! T north
                b_TQU(4,m,ithl) = cmplx(b_ns(3) - b_ns(-3), b_ns(4) - b_ns(-2), kind = DP) ! T south
                b_TQU(2,m,ithl) = cmplx(b_ns(5) + b_ns(-1), b_ns(6) + b_ns( 0), kind = DP) ! Q n
                b_TQU(5,m,ithl) = cmplx(b_ns(5) - b_ns(-1), b_ns(6) - b_ns( 0), kind = DP) ! Q s
                b_TQU(3,m,ithl) = cmplx(b_ns(7) + b_ns( 1), b_ns(8) + b_ns( 2), kind = DP) ! U n
                b_TQU(6,m,ithl) = cmplx(b_ns(7) - b_ns( 1), b_ns(8) - b_ns( 2), kind = DP) ! U s
             endif ! test on nlmax
          enddo ! loop on rings (ithl)
       enddo ! loop on m

       do ithl = 0, uchk-lchk
          nphl = nph(ithl)
          ith  = ithl + lchk
          !        ---------------------------------------------------------------
          !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
          !        ---------------------------------------------------------------
          do i=1,3
             bsub = b_TQU(i,:,ithl)
             call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl))   ! north hemisph. + equator
             local_map(istart_north:istart_north+nphl-1, i) = ring(0:nphl-1)
          enddo
          istart_north = istart_north+nphl


          if (ith < 2*nsmax) then
             istart_south = istart_south-nphl
             do i=1,3
                bsub = b_TQU(3+i,:,ithl)
                call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
                local_map(istart_south:istart_south+nphl-1, i) = ring(0:nphl-1)
             enddo
          endif
       enddo ! loop on ithl
    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate(dalm, lam_fact, ring, bsub)
    deallocate(normal_l)
    deallocate(b_TQU)

    return

  end subroutine mpi_alm2map_pol_pre1


  subroutine mpi_alm2map_pol_pre2

    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !      for the Temperature+Polarization field
    !=======================================================================

    integer(I4B) :: l, m, ith
    integer(I8B) :: istart_south, istart_north, npix
    integer(I4B) :: nrings, nphmx

    integer(I8B) :: n_lm, i_mm, i_up
    integer(I4B)                              :: par_lm
    real(DP),     dimension(-3:8)             :: b_ns
    real(DP),     dimension(:,:), allocatable :: dalm, plm_sub

    integer(i4b)                                :: ll, l_min
    complex(DPC), dimension(:,:,:), allocatable :: b_TQU
    complex(DPC), dimension(:),     allocatable :: bsub
    real(DP),     dimension(:),     allocatable :: ring
    integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i
    real(DP)                           :: cth
    integer(I4B), dimension(0:SMAXCHK) :: nph, kphi0
    integer(I8B), dimension(0:SMAXCHK) :: startpix
    real(DP),     dimension(0:SMAXCHK) :: sth

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = last_ring - first_ring + 1
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(b_TQU(1:6, 0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_TQU')

    allocate(dalm(0:5,0:nlmax),plm_sub(1:3,0:nlmax), stat = status)
    call assert_alloc(status,code,'dalm & plm_sub')
    allocate(ring(0:nphmx-1),bsub(0:nlmax),stat = status)
    call assert_alloc(status,code,'ring & bsub')

    !     ------------ initiate variables and arrays ----------------
    local_map = 0.d0 ! set the whole map to zero
    istart_north = 0
    istart_south = map_size

    ! loop on chunks
    do ichunk = 0, nchunks-1

       lchk = first_ring + ichunk * chunksize 
       uchk = min(lchk+chunksize - 1, last_ring)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth, sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
       enddo
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)


       b_TQU = 0_dpc ! pad with zeros

       do m = 0, nmmax

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             dalm(0,ll) =  real(local_alms(1,ll,m),kind=dp) ! T, real
             dalm(1,ll) = aimag(local_alms(1,ll,m))         ! T, imag
             dalm(2,ll) =  real(local_alms(2,ll,m),kind=dp) ! G, real
             dalm(3,ll) = aimag(local_alms(2,ll,m))        
             dalm(4,ll) =  real(local_alms(3,ll,m),kind=dp) ! C, real
             dalm(5,ll) = aimag(local_alms(3,ll,m))        
          enddo
          do ithl = 0, uchk - lchk
             l_min =  l_min_ylm(m, sth(ithl)) ! lowest l of non-negligible Ylm
             if (nlmax >= l_min) then
                ith = ithl + lchk
                i_mm = n_lm * (ith-first_ring) + ((2_I8B*nlmax + 3 - m)*m)/2 ! location of Ym,m for ring ith
                i_up = i_mm + nlmax - m
                plm_sub(1,m:nlmax) = plm(i_mm:i_up,1)
                plm_sub(2,m:nlmax) = plm(i_mm:i_up,2)
                plm_sub(3,m:nlmax) = plm(i_mm:i_up,3)
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

                b_TQU(1,m,ithl) = cmplx(b_ns(3) + b_ns(-3), b_ns(4) + b_ns(-2), kind = DP) ! T north
                b_TQU(4,m,ithl) = cmplx(b_ns(3) - b_ns(-3), b_ns(4) - b_ns(-2), kind = DP) ! T south
                b_TQU(2,m,ithl) = cmplx(b_ns(5) + b_ns(-1), b_ns(6) + b_ns( 0), kind = DP) ! Q n
                b_TQU(5,m,ithl) = cmplx(b_ns(5) - b_ns(-1), b_ns(6) - b_ns( 0), kind = DP) ! Q s
                b_TQU(3,m,ithl) = cmplx(b_ns(7) + b_ns( 1), b_ns(8) + b_ns( 2), kind = DP) ! U n
                b_TQU(6,m,ithl) = cmplx(b_ns(7) - b_ns( 1), b_ns(8) - b_ns( 2), kind = DP) ! U s
             endif ! test on nlmax
          enddo ! loop on rings (ithl)
       enddo ! loop on m

       do ithl = 0, uchk - lchk
          nphl = nph(ithl)
          ith  = ithl + lchk
          !        ---------------------------------------------------------------
          !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
          !        ---------------------------------------------------------------
          do i=1,3
             bsub = b_TQU(i,:,ithl)
             call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl))   ! north hemisph. + equator
             local_map(istart_north:istart_north+nphl-1, i) = ring(0:nphl-1)
          enddo
          istart_north = istart_north + nphl

          if (ith < 2*nsmax) then
             istart_south = istart_south - nphl
             do i=1,3
                bsub = b_TQU(3+i,:,ithl)
                call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
                local_map(istart_south:istart_south+nphl-1, i) = ring(0:nphl-1)
             enddo
          endif
       enddo ! loop on ithl
    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate(dalm, ring, bsub)
    deallocate(b_TQU)

    return

  end subroutine mpi_alm2map_pol_pre2







  !**************************************************************************
  !
  !             PLM GENERATION
  !
  !**************************************************************************
  !========================================================
  subroutine mpi_plm_gen(polarization)

    logical(lgt), intent(in) :: polarization

    INTEGER(I4B) :: l, m, ith, scalem, scalel, nd2, nrings
    integer(I8B) :: nd1, n_lm, n_plm, i_mm, i_up
    real(DP)            :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)                 :: normal_m, lam_lm1m
    real(DP)                 :: fm2, fl, flm1, fpol
    real(DP)                 :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP)                 :: ovflow, unflow
    real(DP),     dimension(:,:,:), allocatable :: plm_sub
    real(DP),     dimension(:,:), allocatable :: recfac
    real(DP),     dimension(:),   allocatable :: lam_fact
    real(DP),     dimension(:), allocatable :: mfac

    real(DP),     dimension(:),     allocatable :: normal_l
    integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl
    integer(i4b)        :: nph, kphi0, i, j, plm_ind
    integer(i8b)        :: startpix
    real(DP),     dimension(0:SMAXCHK) :: cth, sth
    real(DP),     dimension(0:SMAXCHK) :: one_on_s2, c_on_s2

    INTEGER(I4B) :: status
    character(len=*), parameter :: code = 'PLM_GEN'

    !=================================================================

    ! Healpix definitions
    nrings = (last_ring-first_ring+1)
    n_lm   = ((nmmax+1)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]
    n_plm  = n_lm * (last_ring-first_ring+1)

    if (polarization) then
       allocate(plm(0:n_plm-1,3))
       nd2 = 3
    else
       allocate(plm(0:n_plm-1,1))
       nd2 = 1
    end if

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')
    if (polarization) then
       allocate(normal_l(0:nlmax),stat = status)
       call assert_alloc(status,code,'normal_l')
    endif

    allocate(recfac(0:1,0:nlmax), plm_sub(1:nd2,0:nlmax,0:chunksize-1), stat = status)
    call assert_alloc(status,code,'recfac & plm_sub')
    if (polarization) then
       allocate(lam_fact(0:nlmax),stat = status)
       call assert_alloc(status,code,'lam_fact')
    endif

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax, mfac)
    ! generate Polarization normalisation
    if (polarization) call gen_normpol(nlmax, normal_l)

    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    plm = 0.0_dp

    do ichunk = 0, nchunks-1

       lchk = first_ring + ichunk * chunksize 
       uchk = min(lchk+chunksize - 1, last_ring)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph, startpix, kphi0)
          one_on_s2(ithl) =    1.0_dp / sth(ithl)**2
            c_on_s2(ithl) = cth(ithl) / sth(ithl)**2
       enddo

       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)
          if (polarization) then
             ! generate Ylm relation factor for degree m
             call gen_lamfac(nlmax, m, lam_fact)
             fm2 = real(m * m, kind = DP)
             normal_m = (2.0_dp * m) * (1 - m)
          endif

          do ithl = 0, uchk-lchk
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
             lam_2 = cth_ring * lam_1 * recfac(0,m)

             if (polarization) then
                fpol = normal_m * normal_l(m) * lam_lm
                plm_sub(2, m, ithl) =  fpol * ( 0.5_dp - one_on_s2(ithl) )
                plm_sub(3, m, ithl) =  fpol *              c_on_s2(ithl)
                !
                fm_on_s2     =      m * one_on_s2(ithl)
                two_on_s2    = 2.0_dp * one_on_s2(ithl)
                two_cth_ring = 2.0_dp * cth_ring
                b_w          =  c_on_s2(ithl) 
             endif
             ! ---------- l > m ----------
             do l = m+1, nlmax
                if (polarization) lam_lm1m = lam_lm * lam_fact(l)
                lam_lm = lam_2 * corfac * lam_mm
                plm_sub(1, l, ithl) = lam_lm
                
                if (polarization) then
                   fl = real(l, kind = DP)
                   flm1 = fl - 1.0_dp
                   a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                   a_x =  two_cth_ring * flm1
                   plm_sub(2, l, ithl) =            ( b_w * lam_lm1m - a_w * lam_lm) * normal_l(l)
                   plm_sub(3, l, ithl) = fm_on_s2 * (       lam_lm1m - a_x * lam_lm) * normal_l(l)
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
          enddo ! loop on rings (ithl)

          ! do memory skipping operations outside inner loops
          do ith = lchk, uchk
             i_mm = n_lm*(ith-first_ring) + ((2*nlmax+3-m)*m)/2 ! location of Ym,m for ring ith
             i_up = i_mm + nlmax - m ! location of Ynlmax,m for ring ith
             ithl = ith - lchk
             do i = 1, nd2
                plm(i_mm:i_up, i) = plm_sub(i, m:nlmax, ithl)
             enddo
          enddo

       enddo ! loop on m
    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate (recfac, plm_sub)
    if (polarization) deallocate(lam_fact)
    deallocate(mfac)
    if (polarization) deallocate(normal_l)

    return

  end subroutine mpi_plm_gen




  subroutine mpi_compute_lam_mm(mfac, m, sth, lam_mm, corfac, scalem)
    !=======================================================================
    ! computes lam_mm
    ! the true lam_mm is     lam_mm * corfac
    !=======================================================================
    integer(I4B),            intent(in)  :: m
    real(DP),                intent(in)  :: sth, mfac
!!!    real(DP), dimension(0:), intent(in)  :: mfac
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





end module mpi_alm_tools
