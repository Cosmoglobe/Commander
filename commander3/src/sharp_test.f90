program sharp_test
  use comm_map_mod 
  use comm_param_mod
  use ARS_mod
  implicit none

  integer(i4b)         :: i, j, n, nside, lmax, nmaps, myid, comm, numprocs, ierr
  real(dp)             :: t1, t2, x0, x1, tot
  class(comm_map), pointer      :: map
  class(comm_mapinfo), pointer :: info
  real(dp), allocatable, dimension(:) :: arr1, arr2
  type(planck_rng) :: handle


  n = 1
  nside = 2048
  lmax  = 3*nside
  nmaps = 1

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)


  info => comm_mapinfo(MPI_COMM_WORLD, nside, lmax, nmaps, .false.)
  map  => comm_map(info)
  do i = 0, map%info%np-1
     map%map(i,1) = i
  end do
  call map%writeFITS('test.fits')

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  n = 100
  do i = 1, n
     if (myid == 0) write(*,*) i, n
     do j = 0, map%info%np-1
        map%map(j,1) = j
     end do
     call wall_time(t1)
     call map%Yt()
     call map%Y()
     call wall_time(t2)
     tot = tot + t2-t1
  end do

  if (myid == 0) write(*,*) 'wall time per round-trip sht = ', (tot)/n
  call map%writeFITS('test2.fits')


  call mpi_finalize(ierr)

contains

       function lnL(x)
         use healpix_types
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: lnL

         lnL = -2.5d0 * (1000.d0/x + log(x)) 
         write(*,*) x, lnL

       end function lnL


end program sharp_test
