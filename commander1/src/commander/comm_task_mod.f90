module comm_task_mod
!  use comm_system_mod
  use comm_utils
  use mpi_alm_tools
  implicit none

  type task_list
     integer(i4b) :: fd, n, id, comm
     logical(lgt), dimension(:), allocatable :: done ! Local, not global
  end type

contains

  subroutine init_task_list(tasks, filename, n, comm)
    implicit none
    type(task_list)  :: tasks
    character(len=*) :: filename
    integer(i4b)     :: n, id, ierr, comm

    ! Don't overwrite the old file until everybody is done with it
    call mpi_barrier(comm, ierr)

    call free_task_list(tasks)
    tasks%n    = n
    tasks%comm = comm
    call mpi_comm_rank(comm, tasks%id, ierr)
    tasks%id   = tasks%id+1
    if(tasks%id == 1) then
       tasks%fd = open_atomic_file(trim(filename), 1)
    end if
    call mpi_barrier(comm, ierr)
    if(tasks%id /= 1) then
       tasks%fd = open_atomic_file(trim(filename), 0)
    end if
    allocate(tasks%done(n))
    tasks%done = .false.
    call assert(tasks%fd > 0, "Error opening task list file!")
  end subroutine

  subroutine free_task_list(tasks)
    implicit none
    type(task_list) :: tasks
    tasks%n = 0
    if(allocated(tasks%done)) deallocate(tasks%done)
    if(tasks%fd <= 0) return
    call close_atomic_file(tasks%fd)
    tasks%fd = 0
  end subroutine

  function get_next_task(tasks, res) result(ok)
    implicit none
    type(task_list)                         :: tasks
    integer(i4b), dimension(:), allocatable :: stati
    integer(i4b) :: n, m, i, res
    logical(lgt) :: ok
    n = tasks%n; m = 4*tasks%n
    allocate(stati(n))
    stati = 0
    call lock_atomic_file(tasks%fd, 0, m)
    call read_atomic_file(tasks%fd, 0, m, stati)
    do i = 1, n; if(stati(i) == 0) exit; end do
    res = i
    if(res <= n) then
       stati(res) = tasks%id
       call write_atomic_file(tasks%fd, 0, m, stati)
       ok = .true.
    else
       res = 0
       ok  = .false.
    end if
    call unlock_atomic_file(tasks%fd, 0, m)
    deallocate(stati)
  end function

  function get_task_owner(tasks, taski) result(owner)
    implicit none
    type(task_list) :: tasks
    integer(i4b)    :: taski, owner
    call read_atomic_file(tasks%fd, (taski-1)*4, 4, owner)
    owner = owner-1
  end function

  function claim_task(tasks, taski) result(ok)
    implicit none
    type(task_list) :: tasks
    integer(i4b)    :: taski, owner, m
    logical(lgt)    :: ok
    ok = .false.
    call lock_atomic_file(tasks%fd, 4*(taski-1), 4)
    call read_atomic_file(tasks%fd, 4*(taski-1), 4, owner)
    if(owner == 0) call write_atomic_file(tasks%fd, 4*(taski-1), 4, tasks%id)
    call unlock_atomic_file(tasks%fd, 4*taski, 4)
  end function

  subroutine claim_tasks(tasks, mask)
    implicit none
    type(task_list) :: tasks
    logical(lgt)    :: mask(:)
    integer(i4b), dimension(:), allocatable :: stati
    integer(i4b)    :: n, m
    n = tasks%n; m = 4*tasks%n
    allocate(stati(n))
    call lock_atomic_file(tasks%fd, 0, m)
    call read_atomic_file(tasks%fd, 0, m, stati)
    where(mask) stati = tasks%id
    call write_atomic_file(tasks%fd, 0, m, stati)
    call unlock_atomic_file(tasks%fd, 0, m)
    deallocate(stati)
  end subroutine

end module
