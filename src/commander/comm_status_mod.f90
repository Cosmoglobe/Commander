! This module implements a handy type of shared output file called a status
! file, which prints time, id, memory and lun information at given checkpoints in
! the file. Very useful for debugging resource leaks and timing the program.
! If, for some reason, you want a status file that does not actually do anything,
! pass an empty string for the file name.
!
! Imported from QUIET repository; written by Sigurd Næss
!
module comm_status_mod
  use comm_utils
  use comm_system_mod
  use comm_shared_output_mod
  implicit none
  type status_file
     type(shared_ofile) :: file
     integer(i4b)       :: id
     logical(lgt)       :: active
     real(dp)           :: start_time
  end type
contains
  subroutine init_status(status, fname, communicator)
    implicit none
    type(status_file)      :: status
    character(len=*)       :: fname
    integer(i4b), optional :: communicator
    integer(i4b)           :: comm, ierr
    status%active = .false.
    if(fname == "") return
    comm = mpi_comm_world; if(present(communicator)) comm = communicator
    call open_shared_ofile(status%file, fname, comm)
    call mpi_comm_rank(comm, status%id, ierr)
    call wall_time(status%start_time)
    status%active = .true.
  end subroutine

  subroutine update_status(status, tag)
    implicit none
    type(status_file) :: status
    character(len=*)  :: tag
    character(len=512):: str
    real(dp)          :: mem, max_mem, t
    if(.not. status%active) return
    call wall_time(t)
    mem     = get_mem_use2()    /1024d0**3
    !max_mem = get_max_mem_use()/1024d0**3
    write(str,fmt="(f10.2,i5,f10.5,i6,a)") t-status%start_time, status%id, mem, getlun(), " " // trim(tag)
    call write_shared_ofile(status%file, str)
    call flush_shared_ofile(status%file)
  end subroutine

  subroutine free_status(status)
    implicit none
    type(status_file) :: status
    if(.not. status%active) return
    call close_shared_ofile(status%file)
    status%active = .false.
  end subroutine
end module
