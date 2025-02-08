program commander
  use comm_data_mod
  use comm_param_mod
  implicit none

  integer(i4b)        :: ierr, length
  type(comm_params)   :: cpar
   
  !----------------------------------------------------------------------------------
  ! Command line arguments
  character(len=1000)           :: arg
  integer                     :: arg_indx

  ! Giving the simple command line arguments for user to chose from.
  comm3_args: do arg_indx = 1, command_argument_count()
    call get_command_argument(arg_indx, arg, length, ierr)
  end do comm3_args
  !----------------------------------------------------------------------------------

  ! **************************************************************
  ! *          Get parameters and set up working groups          *
  ! **************************************************************
  call MPI_Init(ierr)

  call read_comm_params(cpar)
  
  call initialize_mpi_struct(cpar)

  ! ************************************************
  ! *               Initialize modules             *
  ! ************************************************

  call initialize_data_mod(cpar)

  ! Wait for everybody to exit
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! And exit
  call mpi_finalize(ierr)



end program commander
