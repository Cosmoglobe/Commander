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
     real(dp),     pointer,     dimension(:,:) :: a
  end type shared_2d_dp


  type shared_3d_dp
     logical(lgt) :: init = .false.
     integer(i4b) :: myid_shared, comm_shared, myid_inter, comm_inter
     !integer(i4b) :: win, wsize, disp_unit
     integer(KIND=MPI_ADDRESS_KIND) :: win, wsize, disp_unit
     type(C_PTR)  :: baseptr
     integer(i4b), allocatable, dimension(:)   :: arrshape
     real(dp),     pointer,     dimension(:,:,:) :: a
  end type shared_3d_dp

  type shared_2d_sp
     logical(lgt) :: init = .false.
     integer(i4b) :: myid_shared, comm_shared, myid_inter, comm_inter
     !integer(i4b) :: win, wsize, disp_unit
     integer(KIND=MPI_ADDRESS_KIND) :: win, wsize, disp_unit
     type(C_PTR)  :: baseptr
     integer(i4b), allocatable, dimension(:)   :: arrshape
     real(sp),     pointer,     dimension(:,:) :: a
  end type shared_2d_sp

  type shared_2d_spc
     logical(lgt) :: init = .false.
     integer(i4b) :: myid_shared, comm_shared, myid_inter, comm_inter
     !integer(i4b) :: win, wsize, disp_unit
     integer(KIND=MPI_ADDRESS_KIND) :: win, wsize, disp_unit
     type(C_PTR)  :: baseptr
     integer(i4b), allocatable, dimension(:)   :: arrshape
     complex(spc), pointer,     dimension(:,:) :: a
  end type shared_2d_spc


  type shared_1d_int
     logical(lgt) :: init = .false.
     integer(i4b) :: myid_shared, comm_shared, myid_inter, comm_inter
     !integer(i4b) :: win, wsize, disp_unit
     integer(KIND=MPI_ADDRESS_KIND) :: win, wsize, disp_unit
     type(C_PTR)  :: baseptr
     integer(i4b), allocatable, dimension(:)   :: arrshape
     integer(i4b), pointer,     dimension(:)  :: a
  end type shared_1d_int

  type shared_2d_int
     logical(lgt) :: init = .false.
     integer(i4b) :: myid_shared, comm_shared, myid_inter, comm_inter
     !integer(i4b) :: win, wsize, disp_unit
     integer(KIND=MPI_ADDRESS_KIND) :: win, wsize, disp_unit
     type(C_PTR)  :: baseptr
     integer(i4b), allocatable, dimension(:)   :: arrshape
     integer(i4b), pointer,     dimension(:,:)  :: a
  end type shared_2d_int


contains

  
  subroutine init_shared_2d_dp(myid_shared, comm_shared, myid_inter, comm_inter, &
       & n, arr)
    implicit none
    integer(i4b),       intent(in)  :: myid_shared, comm_shared, myid_inter, comm_inter
    integer(i4b),       intent(in)  :: n(:)
    type(shared_2d_dp), intent(out) :: arr

    integer(i4b) :: ierr

    arr%myid_shared = myid_shared
    arr%comm_shared = comm_shared
    arr%myid_inter  = myid_inter
    arr%comm_inter  = comm_inter
    arr%init        = .true.
    
    if (arr%myid_shared == 0) then
       arr%wsize = 8*product(n)
    else
       arr%wsize = 0
    end if
    allocate(arr%arrshape(2))
    arr%arrshape  = n
    arr%disp_unit = 1
    call mpi_win_allocate_shared(arr%wsize, arr%disp_unit, MPI_INFO_NULL, &
         & arr%comm_shared, arr%baseptr, arr%win, ierr)
    if (arr%myid_shared /= 0) then
       call mpi_win_shared_query(arr%win, 0, arr%wsize, arr%disp_unit, &
            & arr%baseptr, ierr)
    end if
    call c_f_pointer(arr%baseptr, arr%a, arr%arrshape)

  end subroutine init_shared_2d_dp

  subroutine dealloc_shared_2d_dp(arr)
    implicit none
    type(shared_2d_dp), intent(inout) :: arr

    integer(i4b) :: ierr
  
    call mpi_win_fence(0, arr%win, ierr)
    call mpi_win_free(arr%win, ierr)
!    call mpi_free_mem(arr%baseptr,ierr)
    nullify(arr%a)
    deallocate(arr%arrshape)

  end subroutine dealloc_shared_2d_dp

  subroutine sync_shared_2d_dp_map(arr, ind, val)
    implicit none
    type(shared_2d_dp),                 intent(inout) :: arr
    integer(i4b),       dimension(:),   intent(in)    :: ind 
    real(dp),           dimension(:,:), intent(in)    :: val
    
    integer(i4b) :: ierr

    if (arr%myid_shared == 0) arr%a = 0.d0
    call mpi_win_fence(0, arr%win, ierr)
    arr%a(ind+1,:) = val
    call mpi_win_fence(0, arr%win, ierr)
    if (arr%myid_shared == 0) then
       call mpi_allreduce(MPI_IN_PLACE, arr%a, size(arr%a), &
            & MPI_DOUBLE_PRECISION, MPI_SUM, arr%comm_inter, ierr)
    end if
    call mpi_win_fence(0, arr%win, ierr)

  end subroutine sync_shared_2d_dp_map

  subroutine init_shared_3d_dp(myid_shared, comm_shared, myid_inter, comm_inter, &
       & n, arr)
    implicit none
    integer(i4b),       intent(in)  :: myid_shared, comm_shared, myid_inter, comm_inter
    integer(i4b),       intent(in)  :: n(:)
    type(shared_3d_dp), intent(out) :: arr

    integer(i4b) :: ierr

    arr%myid_shared = myid_shared
    arr%comm_shared = comm_shared
    arr%myid_inter  = myid_inter
    arr%comm_inter  = comm_inter
    arr%init        = .true.

    if (arr%myid_shared == 0) then
       arr%wsize = 8*product(n)
    else
       arr%wsize = 0
    end if
    allocate(arr%arrshape(3))
    arr%arrshape  = n
    arr%disp_unit = 1
    call mpi_win_allocate_shared(arr%wsize, arr%disp_unit, MPI_INFO_NULL, &
         & arr%comm_shared, arr%baseptr, arr%win, ierr)
    if (arr%myid_shared /= 0) then
       call mpi_win_shared_query(arr%win, 0, arr%wsize, arr%disp_unit, &
            & arr%baseptr, ierr)
    end if
    call c_f_pointer(arr%baseptr, arr%a, arr%arrshape)

  end subroutine init_shared_3d_dp

  subroutine dealloc_shared_3d_dp(arr)
    implicit none
    type(shared_3d_dp), intent(inout) :: arr

    integer(i4b) :: ierr
  
    call mpi_win_fence(0, arr%win, ierr)
    call mpi_win_free(arr%win, ierr)
!    call mpi_free_mem(arr%baseptr,ierr)
    nullify(arr%a)
    deallocate(arr%arrshape)

  end subroutine dealloc_shared_3d_dp


  subroutine init_shared_2d_sp(myid_shared, comm_shared, myid_inter, comm_inter, &
       & n, arr)
    implicit none
    integer(i4b),       intent(in)  :: myid_shared, comm_shared, myid_inter, comm_inter
    integer(i4b),       intent(in)  :: n(:)
    type(shared_2d_sp), intent(out) :: arr

    integer(i4b) :: ierr

    arr%myid_shared = myid_shared
    arr%comm_shared = comm_shared
    arr%myid_inter  = myid_inter
    arr%comm_inter  = comm_inter
    arr%init        = .true.
    
    if (arr%myid_shared == 0) then
       arr%wsize = 4*product(n)
    else
       arr%wsize = 0
    end if
    allocate(arr%arrshape(2))
    arr%arrshape  = n
    arr%disp_unit = 1
    call mpi_win_allocate_shared(arr%wsize, arr%disp_unit, MPI_INFO_NULL, &
         & arr%comm_shared, arr%baseptr, arr%win, ierr)
    if (arr%myid_shared /= 0) then
       call mpi_win_shared_query(arr%win, 0, arr%wsize, arr%disp_unit, &
            & arr%baseptr, ierr)
    end if
    call c_f_pointer(arr%baseptr, arr%a, arr%arrshape)

  end subroutine init_shared_2d_sp

  subroutine dealloc_shared_2d_sp(arr)
    implicit none
    type(shared_2d_sp), intent(inout) :: arr

    integer(i4b) :: ierr
  
    call mpi_win_fence(0, arr%win, ierr)
    call mpi_win_free(arr%win, ierr)
!    call mpi_free_mem(arr%baseptr,ierr)
    nullify(arr%a)
    deallocate(arr%arrshape)

  end subroutine dealloc_shared_2d_sp

  subroutine sync_shared_2d_sp_map(arr, ind, val)
    implicit none
    type(shared_2d_sp),                 intent(inout) :: arr
    integer(i4b),       dimension(:),   intent(in)    :: ind 
    real(dp),           dimension(:,:), intent(in)    :: val
    
    integer(i4b) :: ierr

    if (arr%myid_shared == 0) arr%a = 0.
    call mpi_win_fence(0, arr%win, ierr)
    arr%a(:,ind+1) = val
    call mpi_win_fence(0, arr%win, ierr)
    if (arr%myid_shared == 0) then
       call mpi_allreduce(MPI_IN_PLACE, arr%a, size(arr%a), &
            & MPI_REAL, MPI_SUM, arr%comm_inter, ierr)
    end if
    call mpi_win_fence(0, arr%win, ierr)

  end subroutine sync_shared_2d_sp_map



  subroutine init_shared_2d_spc(myid_shared, comm_shared, myid_inter, comm_inter, &
       & n, arr)
    implicit none
    integer(i4b),       intent(in)  :: myid_shared, comm_shared, myid_inter, comm_inter
    integer(i4b),       intent(in)  :: n(:)
    type(shared_2d_spc), intent(out) :: arr

    integer(i4b) :: ierr

    arr%myid_shared = myid_shared
    arr%comm_shared = comm_shared
    arr%myid_inter  = myid_inter
    arr%comm_inter  = comm_inter
    arr%init        = .true.
    
    if (arr%myid_shared == 0) then
       arr%wsize = 8*product(n)
    else
       arr%wsize = 0
    end if
    allocate(arr%arrshape(2))
    arr%arrshape  = n
    arr%disp_unit = 1
    call mpi_win_allocate_shared(arr%wsize, arr%disp_unit, MPI_INFO_NULL, &
         & arr%comm_shared, arr%baseptr, arr%win, ierr)
    if (arr%myid_shared /= 0) then
       call mpi_win_shared_query(arr%win, 0, arr%wsize, arr%disp_unit, &
            & arr%baseptr, ierr)
    end if
    call c_f_pointer(arr%baseptr, arr%a, arr%arrshape)

  end subroutine init_shared_2d_spc

  subroutine dealloc_shared_2d_spc(arr)
    implicit none
    type(shared_2d_spc), intent(inout) :: arr

    integer(i4b) :: ierr
  
    call mpi_win_fence(0, arr%win, ierr)
    call mpi_win_free(arr%win, ierr)
!    call mpi_free_mem(arr%baseptr,ierr)
    nullify(arr%a)
    deallocate(arr%arrshape)

  end subroutine dealloc_shared_2d_spc

  subroutine sync_shared_2d_spc_alm(arr, ind, val)
    implicit none
    type(shared_2d_spc),                 intent(inout) :: arr
    integer(i4b),        dimension(:),   intent(in)    :: ind 
    complex(spc),        dimension(:,:), intent(in)    :: val
    
    integer(i4b) :: ierr

    if (arr%myid_shared == 0) arr%a = 0.
    call mpi_win_fence(0, arr%win, ierr)
    arr%a(:,ind) = val
    call mpi_win_fence(0, arr%win, ierr)
    if (arr%myid_shared == 0) then
       call mpi_allreduce(MPI_IN_PLACE, arr%a, size(arr%a), &
            & MPI_COMPLEX, MPI_SUM, arr%comm_inter, ierr)
    end if
    call mpi_win_fence(0, arr%win, ierr)

  end subroutine sync_shared_2d_spc_alm



  subroutine init_shared_1d_int(myid_shared, comm_shared, myid_inter, comm_inter, &
       & n, arr)
    implicit none
    integer(i4b),       intent(in)  :: myid_shared, comm_shared, myid_inter, comm_inter
    integer(i4b),       intent(in)  :: n(:)
    type(shared_1d_int), intent(out) :: arr

    integer(i4b) :: ierr

    arr%myid_shared = myid_shared
    arr%comm_shared = comm_shared
    arr%myid_inter  = myid_inter
    arr%comm_inter  = comm_inter
    arr%init        = .true.
    
    if (arr%myid_shared == 0) then
       arr%wsize = 4*product(n)
    else
       arr%wsize = 0
    end if
    allocate(arr%arrshape(1))
    arr%arrshape  = n
    arr%disp_unit = 1
    call mpi_win_allocate_shared(arr%wsize, arr%disp_unit, MPI_INFO_NULL, &
         & arr%comm_shared, arr%baseptr, arr%win, ierr)
    if (arr%myid_shared /= 0) then
       call mpi_win_shared_query(arr%win, 0, arr%wsize, arr%disp_unit, &
            & arr%baseptr, ierr)
    end if
    call c_f_pointer(arr%baseptr, arr%a, arr%arrshape)

  end subroutine init_shared_1d_int

  subroutine dealloc_shared_1d_int(arr)
    implicit none
    type(shared_1d_int), intent(inout) :: arr

    integer(i4b) :: ierr
  
    call mpi_win_fence(0, arr%win, ierr)
    call mpi_win_free(arr%win, ierr)
!    call mpi_free_mem(arr%baseptr,ierr)
    nullify(arr%a)
    deallocate(arr%arrshape)

  end subroutine dealloc_shared_1d_int

  subroutine sync_shared_1d_int_map(arr, ind, val)
    implicit none
    type(shared_1d_int),                intent(inout) :: arr
    integer(i4b),       dimension(:),   intent(in)    :: ind 
    integer(i4b),       dimension(:),   intent(in)    :: val
    
    integer(i4b) :: ierr

    if (arr%myid_shared == 0) arr%a = 0
    call mpi_win_fence(0, arr%win, ierr)
    arr%a(ind+1) = val
    call mpi_win_fence(0, arr%win, ierr)
    if (arr%myid_shared == 0) then
       call mpi_allreduce(MPI_IN_PLACE, arr%a, size(arr%a), &
            & MPI_INTEGER, MPI_SUM, arr%comm_inter, ierr)
    end if
    call mpi_win_fence(0, arr%win, ierr)

  end subroutine sync_shared_1d_int_map



  subroutine init_shared_2d_int(myid_shared, comm_shared, myid_inter, comm_inter, &
       & n, arr)
    implicit none
    integer(i4b),       intent(in)  :: myid_shared, comm_shared, myid_inter, comm_inter
    integer(i4b),       intent(in)  :: n(:)
    type(shared_2d_int), intent(out) :: arr

    integer(i4b) :: ierr

    arr%myid_shared = myid_shared
    arr%comm_shared = comm_shared
    arr%myid_inter  = myid_inter
    arr%comm_inter  = comm_inter
    arr%init        = .true.
    
    if (arr%myid_shared == 0) then
       arr%wsize = 4*product(n)
    else
       arr%wsize = 0
    end if
    allocate(arr%arrshape(2))
    arr%arrshape  = n
    arr%disp_unit = 1
    call mpi_win_allocate_shared(arr%wsize, arr%disp_unit, MPI_INFO_NULL, &
         & arr%comm_shared, arr%baseptr, arr%win, ierr)
    if (arr%myid_shared /= 0) then
       call mpi_win_shared_query(arr%win, 0, arr%wsize, arr%disp_unit, &
            & arr%baseptr, ierr)
    end if
    call c_f_pointer(arr%baseptr, arr%a, arr%arrshape)

    if (arr%myid_shared == 0) arr%a = 0
    call mpi_win_fence(0, arr%win, ierr)


  end subroutine init_shared_2d_int

  subroutine dealloc_shared_2d_int(arr)
    implicit none
    type(shared_2d_int), intent(inout) :: arr

    integer(i4b) :: ierr
  
    call mpi_win_fence(0, arr%win, ierr)
    call mpi_win_free(arr%win, ierr)
    nullify(arr%a)
    deallocate(arr%arrshape)

  end subroutine dealloc_shared_2d_int




end module comm_shared_arr_mod
