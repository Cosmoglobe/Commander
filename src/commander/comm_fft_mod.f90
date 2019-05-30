module comm_fft_mod
  use comm_utils
  use comm_param_mod
  use locate_mod
  implicit none

  private
  public initialize_fft_mod, get_closest_fft_magic_number

  integer(i4b) :: min_fft_magic_number, max_fft_magic_number
  integer(i4b), allocatable, dimension(:) :: fft_magic_numbers

contains

  subroutine initialize_fft_mod(cpar)
    implicit none
    type(comm_params),       intent(in) :: cpar

    integer(i4b)       :: i, n, unit

    if (trim(cpar%fft_magic_number_file) == 'none') then
       min_fft_magic_number = 0
       max_fft_magic_number = 0
       return
    end if

    ! Read magic numbers list
    unit = getlun()
    open(unit, file=trim(cpar%fft_magic_number_file))
    read(unit,*) n
    allocate(fft_magic_numbers(n))
    do i = 1, n
       read(unit,*) fft_magic_numbers(i)
    end do
    close(unit)

    min_fft_magic_number = minval(fft_magic_numbers)
    max_fft_magic_number = maxval(fft_magic_numbers)

  end subroutine initialize_fft_mod

  function get_closest_fft_magic_number(i)
    implicit none
    integer(i4b), intent(in) :: i
    integer(i4b)             :: get_closest_fft_magic_number

    integer(i4b) :: ind

    if (i > max_fft_magic_number .or. i < min_fft_magic_number) then
       get_closest_fft_magic_number = i
    else
       ind                          = locate(fft_magic_numbers, i)
       get_closest_fft_magic_number = fft_magic_numbers(ind)
    end if

  end function get_closest_fft_magic_number


end module comm_fft_mod
