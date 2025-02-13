module sharp
  use iso_c_binding
  implicit none
  ! alm_info flags

  ! sharp job types
  integer(c_int) ::  SHARP_YtW = 0

  ! sharp job flags
  integer(c_int), parameter :: SHARP_DP             = ISHFT(1, 4)

  type sharp_geom_info
     type(c_ptr) :: handle
     integer(c_intptr_t) :: n_local
  end type sharp_geom_info

  type sharp_alm_info
     type(c_ptr) :: handle
     integer(c_intptr_t) :: n_local
  end type sharp_alm_info

  interface

     subroutine c_sharp_make_mmajor_real_packed_alm_info( &
         lmax, stride, nm, ms, alm_info) bind(c, name='sharp_make_mmajor_real_packed_alm_info')
       use iso_c_binding
       integer(c_int), value, intent(in)    :: lmax, nm, stride
       integer(c_int), intent(in), optional :: ms(nm)
       type(c_ptr), intent(out)             :: alm_info
     end subroutine c_sharp_make_mmajor_real_packed_alm_info

     function c_sharp_alm_count(alm_info) bind(c, name='sharp_alm_count')
       use iso_c_binding
       integer(c_intptr_t)           :: c_sharp_alm_count
       type(c_ptr), value, intent(in) :: alm_info
     end function c_sharp_alm_count

     ! geom_info
     subroutine sharp_make_subset_healpix_geom_info ( &
          nside, stride, nrings, rings, weight, geom_info) bind(c)
       use iso_c_binding
       integer(c_int), value, intent(in)    :: nside, stride, nrings
       integer(c_int), intent(in), optional :: rings(nrings)
       real(c_double), intent(in), optional :: weight(2 * nside)
       type(c_ptr), intent(out)             :: geom_info
     end subroutine sharp_make_subset_healpix_geom_info

     function c_sharp_map_size(info) bind(c, name='sharp_map_size')
       use iso_c_binding
       integer(c_intptr_t) :: c_sharp_map_size
       type(c_ptr), value   :: info
     end function c_sharp_map_size

     subroutine c_sharp_execute(type, spin, alm, map, geom_info, alm_info, &
                                flags, time, opcnt) bind(c, name='sharp_execute')
       use iso_c_binding
       integer(c_int), value                        :: type, spin, flags
       type(c_ptr), value                           :: alm_info, geom_info
       real(c_double), intent(out), optional        :: time
       integer(c_long_long), intent(out), optional  :: opcnt
       type(c_ptr), intent(in)                      :: alm(*), map(*)
     end subroutine c_sharp_execute

  end interface

contains
  subroutine sharp_make_mmajor_real_packed_alm_info(lmax, ms, alm_info)
    use iso_c_binding
    integer(c_int), value, intent(in)    :: lmax
    integer(c_int), intent(in), optional :: ms(:)
    type(sharp_alm_info), intent(out)    :: alm_info
    !--
    integer(c_int), allocatable          :: ms_copy(:)
    integer(c_int)                       :: nm

    if (present(ms)) then
       nm = size(ms)
       allocate(ms_copy(nm))
       ms_copy = ms
       call c_sharp_make_mmajor_real_packed_alm_info(lmax, 1, nm, ms_copy, alm_info=alm_info%handle)
       deallocate(ms_copy)
    else
       call c_sharp_make_mmajor_real_packed_alm_info(lmax, 1, lmax + 1, alm_info=alm_info%handle)
    end if
    alm_info%n_local = c_sharp_alm_count(alm_info%handle)
  end subroutine sharp_make_mmajor_real_packed_alm_info

  ! geom info
  subroutine sharp_make_healpix_geom_info(nside, rings, weight, geom_info)
    integer(c_int), value                :: nside
    integer(c_int), optional             :: rings(:)
    real(c_double), intent(in), optional :: weight(2 * nside)
    type(sharp_geom_info), intent(out)   :: geom_info
    !--
    integer(c_int) :: nrings
    integer(c_int), allocatable :: rings_copy(:)

    nrings = size(rings)
    allocate(rings_copy(nrings))
    rings_copy = rings
    call sharp_make_subset_healpix_geom_info(nside, 1, nrings, rings_copy, &
                                             weight, geom_info%handle)
    deallocate(rings_copy)
    geom_info%n_local = c_sharp_map_size(geom_info%handle)

  end subroutine sharp_make_healpix_geom_info



  subroutine sharp_execute(alm, alm_info, map, geom_info, &
                             time, opcnt)
    use iso_c_binding
    implicit none

    type(sharp_alm_info)                         :: alm_info
    type(sharp_geom_info)                        :: geom_info
    real(c_double), intent(out), optional        :: time
    integer(c_long_long), intent(out), optional  :: opcnt
    real(c_double), target, intent(inout)        :: alm(0:alm_info%n_local - 1, 1:1)
    real(c_double), target, intent(inout)        :: map(0:geom_info%n_local - 1, 1:1)


    integer(c_int)         :: mod_flags
    type(c_ptr), target    :: alm_ptr(1)
    type(c_ptr), target    :: map_ptr(1)

    real(c_double)        :: time_
    integer(c_long_long)  :: opcnt_

    mod_flags = SHARP_DP

    ! Set up pointer table to access maps
    alm_ptr(:) = c_null_ptr
    map_ptr(:) = c_null_ptr
    alm_ptr(1) = c_loc(alm(0, 1))
    map_ptr(1) = c_loc(map(0, 1))

    call c_sharp_execute(SHARP_YtW, 0, alm_ptr, map_ptr, &
        geom_info=geom_info%handle, &
        alm_info=alm_info%handle, &
        flags=mod_flags, &
        time=time_, &
        opcnt=opcnt_)
    write(*,*) "Passed the bug"

   if (present(time)) time_ = time
   if (present(opcnt)) opcnt_ = opcnt
  end subroutine sharp_execute
end module
