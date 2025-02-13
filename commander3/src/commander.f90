program commander
  use comm_map_mod
  implicit none

  type comm_data_set
     class(comm_mapinfo), pointer :: rmsinfo   => null()
     class(comm_map),     pointer :: map       => null()
     class(comm_N),       pointer :: N         => null()
  end type comm_data_set

  type(comm_data_set), allocatable, dimension(:) :: data

  call initialize_data_mod()
  write(*,*) 'Completed successfully'

  deallocate(data(1)%rmsinfo%ms)
  deallocate(data(1)%rmsinfo%rings)
  deallocate(data(1)%rmsinfo)

  deallocate(data(1)%N%rms0%map)
  deallocate(data(1)%N%rms0%alm)
  deallocate(data(1)%N%rms0)
  deallocate(data(1)%N)

  deallocate(data)

contains

  subroutine initialize_data_mod()
    !
    ! Routine to initialise Commander3 data
    !
    implicit none

    allocate(data(1))
    !                           nside, lmax
    data(1)%rmsinfo  => comm_mapinfo()
    data(1)%N        => comm_N(data(1)%rmsinfo)

  end subroutine initialize_data_mod

end program commander
