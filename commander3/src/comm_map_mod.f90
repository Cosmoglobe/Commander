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
  use iso_c_binding, only : c_ptr, c_double
  implicit none

  public comm_map, comm_mapinfo

  type :: comm_mapinfo
     ! Data variables
     type(sharp_alm_info)  :: alm_info
     type(sharp_geom_info) :: geom_info_T
     integer(c_int), allocatable, dimension(:)   :: rings
     integer(c_int), allocatable, dimension(:)   :: ms
  end type comm_mapinfo

  type :: comm_map
     ! Data variables
     class(comm_mapinfo), pointer :: info => null()
     real(c_double), allocatable, dimension(:,:) :: map
     real(c_double), allocatable, dimension(:,:) :: alm
   contains
     ! Data routines
     procedure     :: YtW_scalar  => exec_sharp_YtW_scalar

  end type comm_map

  interface comm_mapinfo
     procedure constructor_mapinfo
  end interface comm_mapinfo

  interface comm_map
     procedure constructor_map
  end interface comm_map

contains

  !**************************************************
  !             Constructors
  !**************************************************
  function constructor_mapinfo(nside, lmax)
    implicit none
    integer(i4b),                 intent(in) :: nside, lmax
    class(comm_mapinfo), pointer             :: constructor_mapinfo
    class(comm_mapinfo), pointer :: p_new => null()

    allocate(p_new)
    allocate(p_new%rings(3))
    allocate(p_new%ms(4))

    p_new%ms    = [0,1,2,3]
    p_new%rings = [1,3,2]
    
    ! Read ring weights, and create SHARP info structures
    call sharp_make_mmajor_real_packed_alm_info(lmax, ms=p_new%ms, &
         & alm_info=p_new%alm_info)

    call sharp_make_healpix_geom_info(nside, rings=p_new%rings, &
         & geom_info=p_new%geom_info_T)

    constructor_mapinfo => p_new

  end function constructor_mapinfo

  function constructor_map(info)
    implicit none
    class(comm_mapinfo),                intent(in),  target   :: info
    class(comm_map),     pointer                              :: constructor_map

    allocate(constructor_map)
    constructor_map%info => info
    allocate(constructor_map%map(0:11,1), source=0d0)
    allocate(constructor_map%alm(0:15,1), source=0d0)

  end function constructor_map
  
  !**************************************************
  !             Spherical harmonic transforms
  !**************************************************

  subroutine exec_sharp_YtW_scalar(self)
    implicit none

    class(comm_map), intent(inout) :: self

    call sharp_execute(self%alm(:,1:1), self%info%alm_info, &
         & self%map(:,1:1), self%info%geom_info_T)
    
  end subroutine exec_sharp_YtW_scalar
  

end module comm_map_mod
