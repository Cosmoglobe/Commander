program commander
  use comm_N_mod
  use comm_N_rms_mod
  implicit none

  type comm_data_set
     class(comm_mapinfo), pointer :: rmsinfo   => null()
     class(comm_map),     pointer :: map       => null()
     class(comm_N),       pointer :: N         => null()
  end type comm_data_set

  type(comm_data_set), allocatable, dimension(:) :: data

  call initialize_data_mod()
  write(*,*) 'Completed successfully'

contains

  subroutine initialize_data_mod()
    !
    ! Routine to initialise Commander3 data
    !
    implicit none

    allocate(data(1))
    !                           nside, lmax
    data(1)%rmsinfo  => comm_mapinfo(1, 3)
    data(1)%N        => comm_N_rms(data(1)%rmsinfo)

  end subroutine initialize_data_mod

end program commander
