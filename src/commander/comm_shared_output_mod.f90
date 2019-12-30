! Provide shared output to a single file from multiple processes.
! Basically the same as mpi shared writes, but we were asked not to use
! those, so.. Anyway, we could make this a wrapper for mpi write if we want.
!
! The current implementation uses posix locks and maintains an internal
! buffer to avoid excessive lock operations.
!
! Imported from QUIET repository; written by Sigurd Næss
!
module comm_shared_output_mod
  use comm_utils
  use comm_system_mod
  implicit none

  type shared_ofile
     character(len=512)     :: filename
     character(len=1)       :: buffer(1048576)
     integer(i4b)           :: bufind
  end type shared_ofile

contains

  subroutine open_shared_ofile(ofile, fname, communicator)
    implicit none
    type(shared_ofile)     :: ofile
    character(len=*)       :: fname
    integer(i4b), optional :: communicator
    integer(i4b)           :: comm, id, ierr
    comm = mpi_comm_world; if(present(communicator)) comm = communicator
    ofile%filename = fname
    ofile%bufind   = 0
    call mpi_comm_rank(comm, id, ierr)
    if(id == 0) call truncate(trim(fname),0)
    call mpi_barrier(comm, ierr)
  end subroutine

  subroutine write_shared_ofile(ofile, line, advance)
    implicit none
    type(shared_ofile)         :: ofile
    character(len=*)           :: line
    character(len=len(line)+1) :: lcopy
    character(len=*), optional :: advance
    character(len=64)          :: adv
    integer(i4b)               :: i, m, n, free, j
    adv = "yes"; if(present(advance)) adv = advance
    n   = len_trim(line)
    lcopy = line
    if(adv == "yes") then
       n = n+1
       lcopy(n:n) = achar(10)
    end if
    i = 0
    do
       ! Eat as much data from the line as the buffer can take.
       ! Flush and repeat if buffer gets full.
       free = size(ofile%buffer)-ofile%bufind
       m    = min(n-i, free)
       do j = 1, m
          ofile%buffer(ofile%bufind+j) = lcopy(i+j:i+j)
       end do
       ofile%bufind = ofile%bufind + m
       i    = i + m
       if(i == n) then
          exit
       elseif(i < n) then
          call flush_shared_ofile(ofile)
       else
          write(*,*) "Error in write shared ofile"
          write(*,*) n, m, free, i, ofile%bufind
          stop
       end if
    end do
  end subroutine

  subroutine close_shared_ofile(ofile)
    implicit none
    type(shared_ofile) :: ofile
    call flush_shared_ofile(ofile)
  end subroutine

  subroutine flush_shared_ofile(ofile)
    implicit none
    type(shared_ofile) :: ofile
    integer(i4b)       :: cfile
    cfile = open_atomic_file(trim(ofile%filename), 0)
    call lock_atomic_file(cfile, 0, 0)
    call append_atomic_file(cfile, ofile%bufind, ofile%buffer)
    ofile%bufind = 0
    call unlock_atomic_file(cfile, 0, 0)
    call close_atomic_file(cfile)
  end subroutine

end module
