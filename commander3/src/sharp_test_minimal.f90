program sharp_test
  use iso_c_binding
  use sharp
  implicit none
  integer :: i, j, k, nside, lmax, nmap, npix, nalm
  real(c_double), allocatable, target :: map(:,:), alm(:,:)
  type(sharp_alm_info)    :: ainfo
  type(sharp_geom_info)   :: minfo
  type(c_ptr), allocatable:: alm_ptrs(:), map_ptrs(:)


  nside = 8
  lmax  = 3*nside
  nmap  = 1
  npix  = 12*nside**2
  call sharp_make_healpix_geom_info(nside, geom_info=minfo)
  call sharp_make_mmajor_real_packed_alm_info(lmax, alm_info=ainfo)
  !call sharp_make_triangular_alm_info(lmax, lmax, 1, ainfo%handle)
  write(*,*) ainfo%n_local
  nalm = c_sharp_alm_count(ainfo%handle)
  allocate(map(npix,nmap))
  map = 0
  allocate(alm(nalm,nmap))
  alm = 0
  alm(1,1) = 1

  allocate(alm_ptrs(nmap))
  allocate(map_ptrs(nmap))
  do i = 1, nmap
    alm_ptrs(i) = c_loc(alm(1,i))
    map_ptrs(i) = c_loc(map(1,i))
  end do

  call c_sharp_execute(SHARP_Y, 0, alm_ptrs, map_ptrs, &
          geom_info=minfo%handle, alm_info=ainfo%handle, flags=SHARP_DP)

  write(*,*) map(1,1)

end program

!
!
!  allocate(map(npix,nmap))
!
!
!
!  n = 1
!  nside = 64
!  lmax  = 3*nside
!  nmaps = 1
!
!  call mpi_init(ierr)
!  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
!  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
!
!
!  info => comm_mapinfo(MPI_COMM_WORLD, nside, lmax, nmaps, .false.)
!  map  => comm_map(info)
!  do i = 0, map%info%np-1
!     map%map(i,1) = i
!  end do
!  call map%writeFITS('test.fits')
!
!  call mpi_barrier(MPI_COMM_WORLD, ierr)
!
!  n = 100
!  do i = 1, n
!     if (myid == 0) write(*,*) i, n
!     do j = 0, map%info%np-1
!        map%map(j,1) = j
!     end do
!     call wall_time(t1)
!     call map%Yt()
!     call map%Y()
!     call wall_time(t2)
!     tot = tot + t2-t1
!  end do
!
!  if (myid == 0) write(*,*) 'wall time per round-trip sht = ', (tot)/n
!  call map%writeFITS('test2.fits')
!
!
!  call mpi_finalize(ierr)
!
!contains
!
!       function lnL(x)
!         use healpix_types
!         implicit none
!         real(dp), intent(in) :: x
!         real(dp)             :: lnL
!
!         lnL = -2.5d0 * (1000.d0/x + log(x)) 
!         write(*,*) x, lnL
!
!       end function lnL
!
!
!end program sharp_test
