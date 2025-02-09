module comm_param_mod
  use healpix_types
  use mpi
  implicit none

  type comm_params

     ! MPI info
     integer(i4b) :: myid
     integer(i4b) :: myid_chain, comm_chain, mychain

  end type comm_params


contains

  ! ********************************************************
  !                     Driver routines
  ! ********************************************************

  subroutine initialize_mpi_struct(cpar)
    implicit none
    type(comm_params), intent(inout) :: cpar

    integer(i4b) :: ierr

    cpar%mychain    = 1
    cpar%myid_chain = 0

    call mpi_comm_split(MPI_COMM_WORLD, cpar%mychain, cpar%myid_chain, cpar%comm_chain,  ierr) 

  end subroutine initialize_mpi_struct


end module comm_param_mod
