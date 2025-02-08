!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
module comm_map_mod
  use sharp
  use healpix_types
  use fitstools
  use pix_tools
  use udgrade_nr
  use iso_c_binding, only : c_ptr, c_double
  use head_fits
  use comm_hdf_mod
  use extension
  use comm_param_mod
  implicit none

!  include "mpif.h"
      
  public comm_map, comm_mapinfo, map_ptr


  type :: comm_mapinfo
     ! Linked list variables
     class(comm_mapinfo), pointer :: nextLink => null()
     class(comm_mapinfo), pointer :: prevLink => null()

     ! Data variables
     type(sharp_alm_info)  :: alm_info
     type(sharp_geom_info) :: geom_info_T, geom_info_P
     logical(lgt) :: pol, dist
     integer(i4b) :: comm, myid, nprocs
     integer(i4b) :: nside, npix, nmaps, nspec, nring, np, lmax, nm, nalm, mmax
     integer(c_int), allocatable, dimension(:)   :: rings
     integer(c_int), allocatable, dimension(:)   :: ms
     integer(c_int), allocatable, dimension(:)   :: mind
     integer(c_int), allocatable, dimension(:,:) :: lm
     integer(c_int), allocatable, dimension(:)   :: pix
     real(c_double), allocatable, dimension(:,:) :: W
     integer(i4b),   allocatable, dimension(:)   :: nalms
  end type comm_mapinfo

  type :: comm_map
     ! Linked list variables
     class(comm_map), pointer :: nextLink => null()
     class(comm_map), pointer :: prevLink => null()

     ! Data variables
     class(comm_mapinfo), pointer :: info => null()
     real(c_double), allocatable, dimension(:,:) :: map
     real(c_double), allocatable, dimension(:,:) :: alm
     real(c_double), allocatable, dimension(:,:) :: alm_buff
   contains
     ! Data routines
     procedure     :: YtW_scalar  => exec_sharp_YtW_scalar
     procedure     :: readFITS

  end type comm_map

  type map_ptr
     class(comm_map), pointer :: p => null()
  end type map_ptr
  
  interface comm_mapinfo
     procedure constructor_mapinfo
  end interface comm_mapinfo

  interface comm_map
     procedure constructor_map, constructor_clone
  end interface comm_map


  ! Library of mapinfos; resident in memory during the analysis
  class(comm_mapinfo), pointer, private :: mapinfos => null()

contains

  !**************************************************
  !             Constructors
  !**************************************************
  function constructor_mapinfo(comm, nside, lmax, nmaps, pol)
    implicit none
    integer(i4b),                 intent(in) :: comm, nside, lmax, nmaps
    logical(lgt),                 intent(in) :: pol
    class(comm_mapinfo), pointer             :: constructor_mapinfo

    integer(i4b) :: myid, nprocs, ierr
    integer(i4b) :: l, m, i, j, k, np, ind
    real(dp)     :: nullval
    logical(lgt) :: anynull, distval
    integer(i4b), allocatable, dimension(:) :: pixlist
    !real(dp),     allocatable, dimension(:,:) :: weights
    character(len=5)   :: nside_text
    character(len=512) :: weight_file, healpixdir
    class(comm_mapinfo), pointer :: p => null(), p_prev => null(), p_new => null()

    distval = .true.

    ! Check if requested mapinfo already exists in library; if so return a link to that object
    p => mapinfos
    do while (associated(p)) 
       if ((p%nside == nside) .and. (p%lmax == lmax) .and. &
            & (p%nmaps == nmaps) .and. (p%pol .eqv. pol) .and. (p%dist .eqv. distval)) then
          write(*,*) 'Reusing old', nmaps, p%nmaps, nmaps == p%nmaps
          constructor_mapinfo => p
          return
       else
          write(*,*) "Map hasnt existed yet"
          p_prev => p
          p      => p%nextLink
       end if
    end do

    ! Set up new mapinfo object
    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, nprocs, ierr)

    allocate(p_new)
    p_new%comm   = comm
    p_new%myid   = myid
    p_new%nprocs = nprocs
    p_new%nside  = nside
    p_new%nmaps  = nmaps
    p_new%nspec  = nmaps*(nmaps+1)/2
    p_new%lmax   = lmax
    p_new%pol    = pol
    p_new%npix   = 12*nside**2
    p_new%dist   = distval


    ! Select rings and pixels
    allocate(pixlist(0:4*nside-1))
    p_new%nring = 0
    p_new%np    = 0
    do i = 1+myid, 2*nside, nprocs
       call in_ring(nside, i, 0.d0, pi, pixlist, np)
       p_new%nring = p_new%nring + 1
       p_new%np    = p_new%np    + np
       if (i < 2*nside) then ! Symmetric about equator
          p_new%nring = p_new%nring + 1
          p_new%np    = p_new%np    + np
       end if
    end do
    allocate(p_new%rings(p_new%nring))
    allocate(p_new%pix(p_new%np))
    j = 1
    k = 1
    do i = 1+myid, 2*nside, nprocs
       call in_ring(nside, i, 0.d0, pi, pixlist, np)
       p_new%rings(k) = i
       p_new%pix(j:j+np-1) = pixlist(0:np-1)
       k = k + 1
       j = j + np
       if (i < 2*nside) then ! Symmetric about equator
          call in_ring(nside, 4*nside-i, 0.d0, pi, pixlist, np)
          p_new%rings(k) = 4*nside-i
          p_new%pix(j:j+np-1) = pixlist(0:np-1)
          k = k + 1
          j = j + np
       end if
    end do
    deallocate(pixlist)
    call QuickSort_int(p_new%rings)
    call QuickSort_int(p_new%pix)

    ! Select m's
    p_new%nm   = 0
    p_new%nalm = 0
    do m = myid, lmax, nprocs
       p_new%nm   = p_new%nm   + 1
       if (m == 0) then
          p_new%nalm = p_new%nalm + lmax+1
       else
          p_new%nalm = p_new%nalm + 2*(lmax-m+1)
       end if
    end do
    allocate(p_new%ms(p_new%nm))
    allocate(p_new%mind(0:p_new%lmax))
    allocate(p_new%lm(2,0:p_new%nalm-1))
    ind = 0
    p_new%mind = -1
    do i = 1, p_new%nm
       m                           = myid + (i-1)*nprocs
       p_new%ms(i)   = m
       p_new%mind(m) = ind
       if (m == 0) then
          do l = m, lmax
             p_new%lm(:,ind) = [l,m]
             ind                           = ind+1
          end do
       else
          do l = m, lmax
             p_new%lm(:,ind) = [l,+m]
             ind                           = ind+1
             p_new%lm(:,ind) = [l,-m]
             ind                           = ind+1
          end do
       end if
    end do
    
    ! Read ring weights, and create SHARP info structures
    call sharp_make_mmajor_real_packed_alm_info(lmax, ms=p_new%ms, &
         & alm_info=p_new%alm_info)
    call getEnvironment('HEALPIX', healpixdir) 
    call int2string(nside, nside_text)
    weight_file = trim(healpixdir)//'/data/weight_ring_n' // nside_text // '.fits'
    if (nmaps == 1) then
       allocate(p_new%W(2*nside,1))
       call read_dbintab(weight_file, p_new%W, 2*nside, 1, nullval, anynull)
       p_new%W = p_new%W + 1.d0
       call sharp_make_healpix_geom_info(nside, rings=p_new%rings, &
            & weight=p_new%W(:,1), geom_info=p_new%geom_info_T)
    else
       allocate(p_new%W(2*nside,2))
       call read_dbintab(weight_file, p_new%W, 2*nside, 2, nullval, anynull)
       p_new%W(:,:) = p_new%W + 1.d0
       call sharp_make_healpix_geom_info(nside, rings=p_new%rings, &
            & weight=p_new%W(:,1), geom_info=p_new%geom_info_T)
       call sharp_make_healpix_geom_info(nside, rings=p_new%rings, &
            & weight=p_new%W(:,2), geom_info=p_new%geom_info_P)
    end if

    ! Set up new object, and add to list
    p => mapinfos
    if (.not. associated(p)) then
       mapinfos => p_new
    else
       do while (associated(p%nextLink)) 
          p => p%nextLink
       end do
       p%nextLink => p_new
    end if

    !stop lying to mpi_comm about rank and id
    if (.not. distval) then
      call mpi_comm_rank(comm, myid, ierr)
      call mpi_comm_size(comm, nprocs, ierr)
      p_new%myid = myid
      p_new%nprocs = nprocs
    end if
    constructor_mapinfo => p_new

  end function constructor_mapinfo

  function constructor_map(info, filename, mask_misspix, udgrade)
    implicit none
    class(comm_mapinfo),                              intent(in),  target   :: info
    character(len=*),                                 intent(in),  optional :: filename
    logical(lgt),                                     intent(in),  optional :: udgrade
    real(dp),            allocatable, dimension(:,:), intent(out), optional :: mask_misspix
    class(comm_map),     pointer                                            :: constructor_map

    allocate(constructor_map)
    constructor_map%info => info
    allocate(constructor_map%map(0:info%np-1,info%nmaps))
    allocate(constructor_map%alm(0:info%nalm-1,info%nmaps))
    allocate(constructor_map%alm_buff(0:info%nalm-1,info%nmaps))

    constructor_map%alm_buff = 0d0

    if (present(filename)) then
       if (present(mask_misspix)) then
          allocate(mask_misspix(0:info%np-1,info%nmaps))
          call constructor_map%readFITS(filename, mask=mask_misspix, udgrade=udgrade)
       else
          call constructor_map%readFITS(filename, udgrade=udgrade)
       end if
    else
       if (info%np > 0) constructor_map%map = 0.d0
    end if
    if (info%nalm > 0) constructor_map%alm = 0.d0
    
  end function constructor_map

  
  function constructor_clone(map)
    implicit none
    class(comm_map),         intent(in)  :: map
    class(comm_map), pointer             :: constructor_clone

    allocate(constructor_clone)
    constructor_clone%info => map%info
    allocate(constructor_clone%map(0:map%info%np-1,map%info%nmaps))
    allocate(constructor_clone%alm(0:map%info%nalm-1,map%info%nmaps))
    constructor_clone%map = map%map
    constructor_clone%alm = map%alm
    
  end function constructor_clone


  
  !**************************************************
  !             Spherical harmonic transforms
  !**************************************************



  subroutine exec_sharp_YtW_scalar(self)
    implicit none

    class(comm_map), intent(inout) :: self
    integer(i4b) :: i

    do i = 1, self%info%nmaps
       call sharp_execute(SHARP_YtW, 0, 1, self%alm(:,i:i), self%info%alm_info, &
            & self%map(:,i:i), self%info%geom_info_T, comm=self%info%comm)
    end do
    
  end subroutine exec_sharp_YtW_scalar
  
  !**************************************************
  !                   IO routines
  !**************************************************




  
  subroutine readFITS(self, filename, mask, udgrade)
    implicit none

    class(comm_map),                    intent(inout) :: self
    real(dp),         dimension(0:,1:), intent(out), optional :: mask
    logical(lgt),                       intent(in),  optional :: udgrade
    character(len=*),                   intent(in)    :: filename

    integer(i4b) :: i, j, np, npix, ordering, nside, nmaps, ierr, badpix
    real(dp),     allocatable, dimension(:,:) :: map, buffer, map_in
    integer(i4b), allocatable, dimension(:)   :: p
    integer(i4b), dimension(MPI_STATUS_SIZE)  :: mpistat

    ! Check file consistency 
    npix = int(getsize_fits(trim(filename), ordering=ordering, nside=nside, nmaps=nmaps),i4b)
    if (nmaps < self%info%nmaps) then
       if (self%info%myid == 0) write(*,*) 'Incorrect nmaps in ' // trim(filename), nmaps, self%info%nmaps
       call mpi_finalize(ierr)
       stop
    end if
    if (nside /= self%info%nside .and. .not. (present(udgrade))) then
       if (self%info%myid == 0) write(*,*) 'Incorrect nside in ' // trim(filename), 'Expected ', self%info%nside
       call mpi_finalize(ierr)
       stop
    end if

    ! Only the root actually reads from disk; data are distributed via MPI
    if (self%info%myid == 0) then

       ! Read map and convert to RING format if necessary
       allocate(map(0:self%info%npix-1,self%info%nmaps))
       if (present(udgrade)) then
          allocate(map_in(0:npix-1,self%info%nmaps))
          call input_map(filename, map_in, npix, self%info%nmaps, ignore_polcconv=.true.)
          if (ordering == 1) then
             call udgrade_ring(map_in, nside, map, nside_out=self%info%nside)
          else
             call udgrade_nest(map_in, nside, map, nside_out=self%info%nside)
          end if
          deallocate(map_in)
       else
          call input_map(filename, map, self%info%npix, self%info%nmaps, ignore_polcconv=.true.)
       end if

       if (present(mask)) then
          badpix = count(abs(map+1.6375d30) < 1d25 .or. map > 1.d30)
          if (badpix > 0) then
             write(*,fmt='(a,i8)') '   Number of bad pixels = ', badpix
          end if
       end if
       
       if (ordering == 2) then
          do i = 1, self%info%nmaps
             call convert_nest2ring(self%info%nside, map(:,i))
          end do
       end if

       ! Distribute to other nodes
       allocate(p(self%info%npix))
       self%map = map(self%info%pix,:)
       do i = 1, self%info%nprocs-1
          call mpi_recv(np,       1, MPI_INTEGER, i, 98, self%info%comm, mpistat, ierr)
          call mpi_recv(p(1:np), np, MPI_INTEGER, i, 98, self%info%comm, mpistat, ierr)
          allocate(buffer(np,self%info%nmaps))
          buffer = map(p(1:np),:)
          call mpi_send(buffer, np*self%info%nmaps, MPI_DOUBLE_PRECISION, i, 98, &
               & self%info%comm, ierr)
          map(p(1:np),:) = buffer
          deallocate(buffer)
       end do
       deallocate(p, map)
    else
       call mpi_send(self%info%np,               1, MPI_INTEGER, 0, 98, self%info%comm, ierr)
       call mpi_send(self%info%pix, self%info%np,   MPI_INTEGER, 0, 98, self%info%comm, ierr)
       call mpi_recv(self%map,      size(self%map), MPI_DOUBLE_PRECISION, 0, 98, &
            &  self%info%comm, mpistat, ierr)
    end if

    if (present(mask)) then
       do j = 1, self%info%nmaps
          do i = 0, self%info%np-1
             if (abs(self%map(i,j)+1.6375d30) < 1d25 .or. self%map(i,j) == 0.d0 .or. self%map(i,j) > 1.d30) then
                mask(i,j) = 0.d0
             else
                mask(i,j) = 1.d0
             end if
          end do
       end do
    end if
    
  end subroutine readFITS








end module comm_map_mod
