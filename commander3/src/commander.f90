program commander
  use sharp
  use healpix_types
  use iso_c_binding
  implicit none

  real(c_double), allocatable, dimension(:,:) :: map
  real(c_double), allocatable, dimension(:,:) :: alm
  integer(c_int), allocatable, dimension(:)   :: rings
  integer(c_int)                       :: nm
  integer(c_int)                       :: nrings
  integer(c_int)                       :: lmax
  integer(c_int)                       :: nside
  real(c_double), allocatable          :: weight(:)
  type(sharp_alm_info)  :: alm_info
  type(sharp_geom_info) :: geom_info

  allocate(rings(3))
  rings = [1,3,2]
  
  nm = 4
  lmax = 3
  call c_sharp_make_mmajor_real_packed_alm_info(lmax, 1, nm, alm_info=alm_info%handle)
  alm_info%n_local = c_sharp_alm_count(alm_info%handle)

  nrings = 3
  nside = 1
  allocate(weight(2*nside), source=0d0)

  call sharp_make_subset_healpix_geom_info(nside, 1, nrings, rings, &
                                           weight, geom_info%handle)
  geom_info%n_local = c_sharp_map_size(geom_info%handle)

  allocate(map(0:11,1), source=0d0)
  allocate(alm(0:15,1), source=0d0)

  call sharp_execute(alm(:,1:1), alm_info, map(:,1:1), geom_info)

  write(*,*) 'Completed successfully'


end program commander
