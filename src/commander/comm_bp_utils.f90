module comm_bp_utils
  use sort_utils
  use comm_utils
  implicit none 

  interface comp_sz_thermo
     module procedure compute_sz_thermo_single, compute_sz_thermo_array
  end interface comp_sz_thermo

  interface comp_a2t
     module procedure compute_ant2thermo_single, compute_ant2thermo_array
  end interface comp_a2t

  interface comp_bnu_prime
     module procedure compute_bnu_prime_single, compute_bnu_prime_array
  end interface comp_bnu_prime

  interface comp_bnu_prime_RJ
     module procedure compute_bnu_prime_RJ_single, compute_bnu_prime_RJ_array
  end interface comp_bnu_prime_RJ
  
contains

  function compute_sz_thermo_single(nu)
    implicit none
    real(dp), intent(in) :: nu
    real(dp)             :: compute_sz_thermo_single

    real(dp) :: x

    x = h*nu/(k_b*T_cmb)

    compute_sz_thermo_single = T_cmb * (x*(exp(x)+1.d0)/(exp(x)-1.d0)-4.d0)
    
  end function compute_sz_thermo_single

  function compute_sz_thermo_array(nu) result (y)
    implicit none
    real(dp),              dimension(:), intent(in)  :: nu
    real(dp), allocatable, dimension(:)              :: y

    real(dp), allocatable, dimension(:) :: x

    allocate(x(size(nu)), y(size(nu)))
    x = h*nu/(k_b*T_cmb)
    y = T_cmb * (x*(exp(x)+1.d0)/(exp(x)-1.d0)-4.d0)
    deallocate(x)

  end function compute_sz_thermo_array
  

  function compute_ant2thermo_array(nu)
    implicit none
    real(dp), dimension(:),        intent(in)  :: nu
    real(dp), dimension(size(nu))              :: compute_ant2thermo_array

    integer(i4b) :: i
    real(dp)     :: x

    do i = 1, size(nu)
       x = h*nu(i) / (k_B*T_CMB)
       compute_ant2thermo_array(i) = (exp(x)-1.d0)**2 / (x**2 * exp(x))
    end do
    
  end function compute_ant2thermo_array

  function compute_ant2thermo_single(nu)
    implicit none
    real(dp), intent(in)  :: nu
    real(dp)              :: compute_ant2thermo_single

    real(dp)     :: x

    x = h*nu / (k_B*T_CMB)
    compute_ant2thermo_single = (exp(x)-1.d0)**2 / (x**2 * exp(x))
    
  end function compute_ant2thermo_single

  function compute_bnu_prime_array(nu)
    implicit none

    real(dp), dimension(:),        intent(in)  :: nu
    real(dp), dimension(size(nu))              :: compute_bnu_prime_array

    integer(i4b) :: i
    real(dp)     :: x

    do i = 1, size(nu)
       x = h*nu(i) / (k_B*T_CMB)
       compute_bnu_prime_array(i) = (2.d0*h*nu(i)**3/(c**2*(exp(x)-1.d0))) * &
            & (exp(x) / (exp(x)-1.d0)) * h*nu(i)/(k_B*T_CMB**2)
    end do
    
  end function compute_bnu_prime_array

  function compute_bnu_prime_single(nu)
    implicit none

    real(dp), intent(in)  :: nu
    real(dp)              :: compute_bnu_prime_single

    real(dp)     :: x

    x = h*nu / (k_B*T_CMB)
    compute_bnu_prime_single = (2.d0*h*nu**3/(c**2*(exp(x)-1.d0))) * &
         & (exp(x) / (exp(x)-1.d0)) * h*nu/(k_B*T_CMB**2)
    
  end function compute_bnu_prime_single

  function compute_bnu_prime_RJ_array(nu)
    implicit none

    real(dp), dimension(:),        intent(in)  :: nu
    real(dp), dimension(size(nu))              :: compute_bnu_prime_RJ_array

    compute_bnu_prime_RJ_array = 2.d0*k_B*nu**2/c**2
    
  end function compute_bnu_prime_RJ_array

  function compute_bnu_prime_RJ_single(nu)
    implicit none

    real(dp), intent(in)  :: nu
    real(dp)              :: compute_bnu_prime_RJ_single

    compute_bnu_prime_RJ_single = 2.d0*k_B*nu**2/c**2
    
  end function compute_bnu_prime_RJ_single

  function dBdnu(nu)
    implicit none

    real(dp), intent(in)  :: nu
    real(dp)              :: dBdnu

    real(dp)     :: x

    x = h*nu / (k_B*T_CMB)
    dBdnu = (2*h*nu**3/(c**2 * (exp(x)-1.d0))) *  (exp(x) / (exp(x)-1.d0)) * (h*nu/(k_b*T_CMB**2))
    
  end function dBdnu

  function dB_rj_dnu(nu)
    implicit none

    real(dp), intent(in)  :: nu
    real(dp)              :: dB_rj_dnu

    dB_rj_dnu = 2*nu**2*k_b / c**2
    
  end function dB_rj_dnu

  ! Routine for reading bandpass files
  subroutine read_bandpass(filename, threshold, n, nu, tau)
    implicit none

    character(len=*),                            intent(in)  :: filename
    real(dp),                                    intent(in)  :: threshold
    integer(i4b),                                intent(out) :: n
    real(dp),         allocatable, dimension(:), intent(out) :: nu, tau

    integer(i4b)        :: unit, first, last, m, ierr
    logical(lgt)        :: exist
    character(len=128)  :: string
    real(dp), allocatable, dimension(:) :: x, y

    unit = getlun()
    
    inquire(file=trim(filename), exist=exist)
    if (.not. exist) call report_error('Bandpass file does not exist = ' // trim(filename))

    ! Find the number of entries
    m = 0
    open(unit, file=trim(filename))
    do while (.true.)
       read(unit,*,end=1) string
       if (string(1:1)=='#') cycle
       m = m + 1
    end do
1   close(unit)

    if (m == 0) call report_error('No valid data entries in bandpass file ' // trim(filename))

    allocate(x(m), y(m))
    m = 0
    open(unit, file=trim(filename))
    do while (.true.)
       read(unit,fmt='(a)',end=2) string
       if (string(1:1)=='#') cycle
       m = m+1
       read(string,*) x(m), y(m)

       ! Drop double entries
       if (m > 1) then
          if (x(m) == x(m-1)) m = m-1
       end if
    end do
2   close(unit)

    x = x * 1.d9 ! Convert from GHz to Hz

    first = 1
    last  = m
    if (threshold > 0.d0) then
       do while (y(first) < threshold*maxval(y))
          first = first+1
       end do
       do while (y(last) < threshold*maxval(y))
          last = last-1
       end do
    end if
    
    n = last-first+1
    allocate(nu(n), tau(n))
    nu  = x(first:last)
    tau = y(first:last)

    deallocate(x, y)

  end subroutine read_bandpass


end module comm_bp_utils
