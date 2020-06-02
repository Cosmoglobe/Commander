module locate_mod
  use healpix_types
  implicit none

  interface locate
     module procedure locate_int, locate_dp
  end interface

contains

  ! Find the index of the last element in xx which is < x
  ! The result will always be a valid index in x. Except
  ! for result == n, x will be between(xx(result) and xx(result+1)).
  function locate_int(xx, x)
    implicit none

    integer(i4b),               intent(in) :: x
    integer(i4b), dimension(:), intent(in) :: xx
    integer(i4b)                           :: locate_int

    integer(i4b) :: n, jl, jm, ju
    logical(lgt) :: ascnd

    n     = size(xx)
    ascnd = (xx(n) >= xx(1))
    jl    = 0
    ju    = n+1

    do 
       if (ju-jl <= 1) exit
       jm = (ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       end if
    end do

    if (x == xx(1)) then
       locate_int = 1
    else if (x == xx(n)) then
       locate_int = n-1
    else
       locate_int = jl
    end if

  end function locate_int

  function locate_dp(xx, x)
    implicit none

    real(dp),               intent(in) :: x
    real(dp), dimension(:), intent(in) :: xx
    integer(i4b)                       :: locate_dp

    integer(i4b) :: n, jl, jm, ju
    logical(lgt) :: ascnd

    n     = size(xx)
    ascnd = (xx(n) >= xx(1))
    jl    = 0
    ju    = n+1

    do 
       if (ju-jl <= 1) exit
       jm = (ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       end if
    end do

    if (x == xx(1)) then
       locate_dp = 1
    else if (x == xx(n)) then
       locate_dp = n-1
    else
       locate_dp = jl
    end if

  end function locate_dp


end module
