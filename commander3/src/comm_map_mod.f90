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
  use extension
  use comm_param_mod
  implicit none

  public comm_map, comm_mapinfo, map_ptr


  type :: comm_mapinfo
     ! Data variables
     type(sharp_alm_info)  :: alm_info
     type(sharp_geom_info) :: geom_info_T, geom_info_P
     logical(lgt) :: pol, dist
     integer(i4b) :: comm, myid, nprocs
     integer(i4b) :: nside, npix, nmaps, nring, np, lmax, nm, nalm, mmax
     integer(c_int), allocatable, dimension(:)   :: rings
     integer(c_int), allocatable, dimension(:)   :: ms
     integer(c_int), allocatable, dimension(:)   :: mind
     integer(c_int), allocatable, dimension(:,:) :: lm
     integer(c_int), allocatable, dimension(:)   :: pix
     real(c_double), allocatable, dimension(:,:) :: W
     integer(i4b),   allocatable, dimension(:)   :: nalms
  end type comm_mapinfo

  type :: comm_map
     ! Data variables
     class(comm_mapinfo), pointer :: info => null()
     real(c_double), allocatable, dimension(:,:) :: map
     real(c_double), allocatable, dimension(:,:) :: alm
     real(c_double), allocatable, dimension(:,:) :: alm_buff
   contains
     ! Data routines
     procedure     :: YtW_scalar  => exec_sharp_YtW_scalar

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



contains

  !**************************************************
  !             Constructors
  !**************************************************
  function constructor_mapinfo(comm, nside, lmax, nmaps)
    implicit none
    integer(i4b),                 intent(in) :: comm, nside, lmax, nmaps
    class(comm_mapinfo), pointer             :: constructor_mapinfo

    integer(i4b) :: myid, nprocs, ierr
    integer(i4b) :: l, m, i, j, k, np, ind
    integer(i4b), allocatable, dimension(:) :: pixlist
    class(comm_mapinfo), pointer :: p_new => null()

    ! Set up new mapinfo object
    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, nprocs, ierr)

    allocate(p_new)
    p_new%comm   = comm
    p_new%myid   = myid
    p_new%nprocs = nprocs
    p_new%nside  = nside
    p_new%nmaps  = nmaps
    p_new%lmax   = lmax


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
             p_new%lm(1,ind) = l
             p_new%lm(2,ind) = m
             ind                           = ind+1
          end do
       else
          do l = m, lmax
             p_new%lm(1,ind) = l
             p_new%lm(2,ind) = m
             ind                           = ind+1
             p_new%lm(1,ind) = l
             p_new%lm(2,ind) = -m
             ind                           = ind+1
          end do
       end if
    end do
    
    ! Read ring weights, and create SHARP info structures
    call sharp_make_mmajor_real_packed_alm_info(lmax, ms=p_new%ms, &
         & alm_info=p_new%alm_info)

    allocate(p_new%W(2*nside,1))
    p_new%W = 0d0 
    call sharp_make_healpix_geom_info(nside, rings=p_new%rings, &
         & weight=p_new%W(:,1), geom_info=p_new%geom_info_T)


    constructor_mapinfo => p_new

  end function constructor_mapinfo

  function constructor_map(info)
    implicit none
    class(comm_mapinfo),                intent(in),  target   :: info
    class(comm_map),     pointer                              :: constructor_map

    allocate(constructor_map)
    constructor_map%info => info
    allocate(constructor_map%map(0:info%np-1,info%nmaps))
    allocate(constructor_map%alm(0:info%nalm-1,info%nmaps))
    allocate(constructor_map%alm_buff(0:info%nalm-1,info%nmaps))

    constructor_map%alm_buff = 0d0
    if (info%np > 0) constructor_map%map = 0.d0
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
    constructor_clone%map = 0d0
    constructor_clone%alm = 0d0
    
  end function constructor_clone
  
  !**************************************************
  !             Spherical harmonic transforms
  !**************************************************

  subroutine exec_sharp_YtW_scalar(self)
    implicit none

    class(comm_map), intent(inout) :: self

    call sharp_execute(SHARP_YtW, 0, 1, self%alm(:,1:1), self%info%alm_info, &
         & self%map(:,1:1), self%info%geom_info_T, comm=self%info%comm)
    
  end subroutine exec_sharp_YtW_scalar
  

end module comm_map_mod
