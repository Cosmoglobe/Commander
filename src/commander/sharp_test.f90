program sharp_test
  use comm_map_mod
  use comm_param_mod
  implicit none

  integer(i4b)         :: i, j, n, nside, lmax, nmaps, myid, comm, numprocs, ierr
  real(dp)             :: t1, t2
  class(comm_map), pointer      :: map
  class(comm_mapinfo), pointer :: info
  real(dp), allocatable, dimension(:) :: arr1, arr2

  n = 1
  nside = 2048
  lmax  = 3*nside
  nmaps = 1

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)

  nside = 1000000
  allocate(arr1(nside), arr2(nside))

  arr1 = myid
  call wall_time(t1)
  call mpi_allreduce(arr1, arr2, nside, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call wall_time(t2)
  if (myid == 0) write(*,*) 'wall time = ', t2-t1, ', sum = ', sum(arr2)


!!$  info => comm_mapinfo(MPI_COMM_WORLD, nside, lmax, nmaps, .false.)
!!$  map  => comm_map(info)
!!$  do i = 0, map%info%np-1
!!$     map%map(i,1) = i
!!$  end do
!!$  call map%writeFITS('test.fits')
!!$
!!$  n = 5
!!$  call wall_time(t1)
!!$  do i = 1, n
!!$     call map%Yt()
!!$     call map%Y()
!!$  end do
!!$  call wall_time(t2)
!!$  if (myid == 0) write(*,*) 'wall time per round-trip sht = ', (t2-t1)/n
!!$  call map%writeFITS('test2.fits')


  call mpi_finalize(ierr)

end program sharp_test
