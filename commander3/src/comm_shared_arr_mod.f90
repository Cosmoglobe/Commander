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
module comm_shared_arr_mod
  use comm_utils
  implicit none

  type shared_2d_dp
     logical(lgt) :: init = .false.
     integer(i4b) :: myid_shared, comm_shared, myid_inter, comm_inter
     !integer(i4b) :: win, wsize, disp_unit
     integer(KIND=MPI_ADDRESS_KIND) :: win, wsize, disp_unit
     type(C_PTR)  :: baseptr
     integer(i4b), allocatable, dimension(:)   :: arrshape
     real(dp),     pointer,     dimension(:,:) :: a => null()
  end type shared_2d_dp


  type shared_3d_dp
     logical(lgt) :: init = .false.
     integer(i4b) :: myid_shared, comm_shared, myid_inter, comm_inter
     !integer(i4b) :: win, wsize, disp_unit
     integer(KIND=MPI_ADDRESS_KIND) :: win, wsize, disp_unit
     type(C_PTR)  :: baseptr
     integer(i4b), allocatable, dimension(:)   :: arrshape
     real(dp),     pointer,     dimension(:,:,:) :: a => null()
  end type shared_3d_dp

  type shared_2d_sp
     logical(lgt) :: init = .false.
     integer(i4b) :: myid_shared, comm_shared, myid_inter, comm_inter
     !integer(i4b) :: win, wsize, disp_unit
     integer(KIND=MPI_ADDRESS_KIND) :: win, wsize, disp_unit
     type(C_PTR)  :: baseptr
     integer(i4b), allocatable, dimension(:)   :: arrshape
     real(sp),     pointer,     dimension(:,:) :: a => null()
  end type shared_2d_sp

  type shared_2d_spc
     logical(lgt) :: init = .false.
     integer(i4b) :: myid_shared, comm_shared, myid_inter, comm_inter
     !integer(i4b) :: win, wsize, disp_unit
     integer(KIND=MPI_ADDRESS_KIND) :: win, wsize, disp_unit
     type(C_PTR)  :: baseptr
     integer(i4b), allocatable, dimension(:)   :: arrshape
     complex(spc), pointer,     dimension(:,:) :: a => null()
  end type shared_2d_spc


  type shared_1d_int
     logical(lgt) :: init = .false.
     integer(i4b) :: myid_shared, comm_shared, myid_inter, comm_inter
     !integer(i4b) :: win, wsize, disp_unit
     integer(KIND=MPI_ADDRESS_KIND) :: win, wsize, disp_unit
     type(C_PTR)  :: baseptr
     integer(i4b), allocatable, dimension(:)   :: arrshape
     integer(i4b), pointer,     dimension(:)   :: a => null()
  end type shared_1d_int

  type shared_2d_int
     logical(lgt) :: init = .false.
     integer(i4b) :: myid_shared, comm_shared, myid_inter, comm_inter
     !integer(i4b) :: win, wsize, disp_unit
     integer(KIND=MPI_ADDRESS_KIND) :: win, wsize, disp_unit
     type(C_PTR)  :: baseptr
     integer(i4b), allocatable, dimension(:)   :: arrshape
     integer(i4b), pointer,     dimension(:,:)  :: a => null()
  end type shared_2d_int


contains

  
  subroutine init_shared_2d_dp(myid_shared, comm_shared, myid_inter, comm_inter, &
       & n, arr)
    implicit none
    integer(i4b),       intent(in)  :: myid_shared, comm_shared, myid_inter, comm_inter
    integer(i4b),       intent(in)  :: n(:)
    type(shared_2d_dp), intent(out) :: arr

    integer(i4b) :: ierr, i


  end subroutine init_shared_2d_dp

  subroutine dealloc_shared_2d_dp(arr)
    implicit none
    type(shared_2d_dp), intent(inout) :: arr

    integer(i4b) :: ierr
  

  end subroutine dealloc_shared_2d_dp

  subroutine sync_shared_2d_dp_map(arr, ind, val)
    implicit none
    type(shared_2d_dp),                 intent(inout) :: arr
    integer(i4b),       dimension(:),   intent(in)    :: ind 
    real(dp),           dimension(:,:), intent(in)    :: val
    
    integer(i4b) :: ierr


  end subroutine sync_shared_2d_dp_map

  subroutine init_shared_3d_dp(myid_shared, comm_shared, myid_inter, comm_inter, &
       & n, arr)
    implicit none
    integer(i4b),       intent(in)  :: myid_shared, comm_shared, myid_inter, comm_inter
    integer(i4b),       intent(in)  :: n(:)
    type(shared_3d_dp), intent(out) :: arr

    integer(i4b) :: ierr, i


  end subroutine init_shared_3d_dp

  subroutine dealloc_shared_3d_dp(arr)
    implicit none
    type(shared_3d_dp), intent(inout) :: arr

    integer(i4b) :: ierr
  

  end subroutine dealloc_shared_3d_dp


  subroutine init_shared_2d_sp(myid_shared, comm_shared, myid_inter, comm_inter, &
       & n, arr)
    implicit none
    integer(i4b),       intent(in)  :: myid_shared, comm_shared, myid_inter, comm_inter
    integer(i4b),       intent(in)  :: n(:)
    type(shared_2d_sp), intent(out) :: arr

    integer(i4b) :: ierr, i


  end subroutine init_shared_2d_sp

  subroutine dealloc_shared_2d_sp(arr)
    implicit none
    type(shared_2d_sp), intent(inout) :: arr

    integer(i4b) :: ierr
  

  end subroutine dealloc_shared_2d_sp

  subroutine sync_shared_2d_sp_map(arr, ind, val)
    implicit none
    type(shared_2d_sp),                 intent(inout) :: arr
    integer(i4b),       dimension(:),   intent(in)    :: ind 
    real(dp),           dimension(:,:), intent(in)    :: val
    
    integer(i4b) :: ierr

  end subroutine sync_shared_2d_sp_map



  subroutine init_shared_2d_spc(myid_shared, comm_shared, myid_inter, comm_inter, &
       & n, arr)
    implicit none
    integer(i4b),       intent(in)  :: myid_shared, comm_shared, myid_inter, comm_inter
    integer(i4b),       intent(in)  :: n(:)
    type(shared_2d_spc), intent(out) :: arr

    integer(i4b) :: ierr, i


  end subroutine init_shared_2d_spc

  subroutine dealloc_shared_2d_spc(arr)
    implicit none
    type(shared_2d_spc), intent(inout) :: arr

    integer(i4b) :: ierr


  end subroutine dealloc_shared_2d_spc

  subroutine sync_shared_2d_spc_alm(arr, ind, val)
    implicit none
    type(shared_2d_spc),                 intent(inout) :: arr
    integer(i4b),        dimension(:),   intent(in)    :: ind 
    complex(spc),        dimension(:,:), intent(in)    :: val
    
    integer(i4b) :: ierr


  end subroutine sync_shared_2d_spc_alm



  subroutine init_shared_1d_int(myid_shared, comm_shared, myid_inter, comm_inter, &
       & n, arr)
    implicit none
    integer(i4b),        intent(in)  :: myid_shared, comm_shared, myid_inter, comm_inter
    integer(i4b),        intent(in)  :: n(:)
    type(shared_1d_int), intent(out) :: arr

    integer(i4b) :: ierr


  end subroutine init_shared_1d_int

  subroutine dealloc_shared_1d_int(arr)
    implicit none
    type(shared_1d_int), intent(inout) :: arr

    integer(i4b) :: ierr
  

  end subroutine dealloc_shared_1d_int

  subroutine sync_shared_1d_int_map(arr, ind, val)
    implicit none
    type(shared_1d_int),                intent(inout) :: arr
    integer(i4b),       dimension(:),   intent(in)    :: ind 
    integer(i4b),       dimension(:),   intent(in)    :: val
    
    integer(i4b) :: ierr

  end subroutine sync_shared_1d_int_map



  subroutine init_shared_2d_int(myid_shared, comm_shared, myid_inter, comm_inter, &
       & n, arr)
    implicit none
    integer(i4b),       intent(in)  :: myid_shared, comm_shared, myid_inter, comm_inter
    integer(i4b),       intent(in)  :: n(:)
    type(shared_2d_int), intent(out) :: arr

    integer(i4b) :: ierr, i



  end subroutine init_shared_2d_int

  subroutine dealloc_shared_2d_int(arr)
    implicit none
    type(shared_2d_int), intent(inout) :: arr

    integer(i4b) :: ierr

  end subroutine dealloc_shared_2d_int




end module comm_shared_arr_mod
