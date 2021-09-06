program test
  use iso_c_binding
  use comm_mpi_mod
  use sharp
  implicit none

  interface
     subroutine c_sharp_make_weighted_healpix_geom_info (nside, stride, weight, geom_info) bind(c, name='sharp_make_weighted_healpix_geom_info')
       use iso_c_binding
       integer(c_int), value, intent(in)    :: nside, stride
       real(c_double), intent(in)           :: weight(:)
       type(c_ptr), intent(out)             :: geom_info
     end subroutine
  end interface

  integer :: ierr, i
  integer(c_int) :: type, spin, nside, lmax, flags
  type(sharp_geom_info) :: minfo
  type(sharp_alm_info)  :: ainfo
  real(c_double), target, allocatable :: map(:), alm(:), weight(:)
  integer(c_int), allocatable :: ms(:)
  type(c_ptr) :: pmap(1), palm(1)

  call mpi_init(ierr)

  nside= 512
  type = 2 !!YtW map2alm
  spin = 0
  lmax = nside*2
  flags= 16

  allocate(weight(2*nside))
  weight = 100000000000

  call c_sharp_make_weighted_healpix_geom_info(nside, 1, weight, minfo%handle)
  minfo%n_local = c_sharp_map_size(minfo%handle)
  allocate(ms(lmax+1))
  do i = 0, lmax
    ms(i+1) = i
  end do
  call c_sharp_make_mmajor_real_packed_alm_info(lmax, 1, lmax+1, ms, ainfo%handle)
  ainfo%n_local = c_sharp_alm_count(ainfo%handle)
  allocate(map(minfo%n_local))
  allocate(alm(ainfo%n_local))
  map = 1
  alm = 7
  pmap(1) = c_loc(map)
  palm(1) = c_loc(alm)

  call c_sharp_execute_mpi(MPI_COMM_WORLD, type, spin, palm, pmap, minfo%handle, ainfo%handle, flags)
  write(*,*) map(1)
  write(*,*) alm(1)

end program
