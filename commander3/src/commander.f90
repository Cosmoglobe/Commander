program commander
  use comm_data_mod
  use comm_param_mod
  implicit none

  integer(i4b)        :: ierr
  type(comm_params)   :: cpar
   
  call MPI_Init(ierr)
  if (ierr .ne. 0) write(*,*) 'mpi init error'

  call initialize_mpi_struct(cpar)

  call initialize_data_mod(cpar)

  call mpi_barrier(MPI_COMM_WORLD, ierr)
  if (ierr .ne. 0) write(*,*) 'mpi barrier error'

  ! And exit
  call mpi_finalize(ierr)
  if (ierr .ne. 0) write(*,*) 'mpi finalize error'

  write(*,*) 'Completed successfully'

end program commander
