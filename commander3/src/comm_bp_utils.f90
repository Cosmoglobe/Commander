!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
module comm_bp_utils
  use sort_utils
  use comm_utils
  use comm_hdf_mod
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
  

  function compute_ant2thermo_array(nu) result(a2t)
    implicit none
    real(dp), dimension(:),        intent(in)  :: nu
    real(dp), dimension(size(nu))              :: a2t

    integer(i4b) :: i
    real(dp)     :: x

    write(*,*) 'Do NOT use this function -- results in NaNs on some compilers/systems!'
    stop
    do i = 1, size(nu)
       x = h*nu(i) / (k_B*T_CMB)
       a2t(i) = (exp(x)-1.d0)**2 / (x**2 * exp(x))
    end do
    
  end function compute_ant2thermo_array

  function compute_ant2thermo_single(nu)
    implicit none
    real(dp), intent(in)  :: nu
    real(dp)              :: compute_ant2thermo_single

    real(dp)     :: x

    if (k_b <= 0.d0 .or. T_CMB <= 0.d0) then
      write(*,*) 'h, nu, k_b, T_CMB'
      write(*,*) h, nu, k_b, T_CMB
    end if
    x = h*nu / (k_B*T_CMB)
    if (x > 200) then
      write(*,*) 'h, nu, k_B, T_CMB, x = ', h, nu, k_b, T_CMB, x
    end if
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

  ! Routine for reading bandpass files for one detecor with threshold
  subroutine read_bandpass(filename, label, threshold, n, nu, tau)
    implicit none

    character(len=*),                            intent(in)  :: filename
    character(len=*),                            intent(in)  :: label
    real(dp),                                    intent(in)  :: threshold
    integer(i4b),                                intent(out) :: n
    real(dp),         allocatable, dimension(:), intent(out) :: nu, tau

    integer(i4b)        :: unit, first, last, m, ierr, l, ext(1)
    logical(lgt)        :: exist
    character(len=128)  :: string
    type(hdf_file)     :: file
    real(dp), allocatable, dimension(:) :: x, y

    unit = getlun()
    
    inquire(file=trim(filename), exist=exist)
    if (.not. exist) call report_error('Bandpass file does not exist = ' // trim(filename))

    l = len(trim(filename))

    if (filename(l-2:l) == '.h5' .or. filename(l-3:l) == '.hd5') then

       call open_hdf_file(filename, file, "r")
       call get_size_hdf(file, trim(label) // "/bandpass", ext)
       m = ext(1)

       allocate(x(m), y(m))
       !write(*,*) "About to read bandpass"
       call read_hdf(file, trim(label) // "/bandpassx",x)
       call read_hdf(file, trim(label) // "/bandpass", y)
       call close_hdf_file(file)

       ! Drop double entries
       l = 1
       do while (l < m)
          if (x(l) == x(l+1)) then
             x(l:m-1) = x(l+1:m)
             y(l:m-1) = y(l+1:m)
             m        = m-1
          else
             l = l+1
          end if
       end do
       
    else 
       ! Assume ASCII

       ! Find the number of entries
       m = 0
       open(unit, file=trim(filename))
       do while (.true.)
          read(unit,*,end=1) string
          if (string(1:1)=='#') cycle
          m = m + 1
       end do
1      close(unit)
       
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
2      close(unit)
    end if

    x(1:m) = x(1:m) * 1.d9 ! Convert from GHz to Hz

    first = 1
    last  = m
    if (threshold > 0.d0) then
       do while (y(first) < threshold*maxval(y(1:m)))
          first = first+1
       end do
       do while (y(last) < threshold*maxval(y(1:m)))
          last = last-1
       end do
    end if
    
    n = last-first+1
    allocate(nu(n), tau(n))
    nu  = x(first:last)
    tau = y(first:last)

    deallocate(x, y)

  end subroutine read_bandpass

  ! Routine for reading bandpass files
  subroutine read_bandpass_dirbe(filename, label, threshold, n, nu, tau)
    implicit none

    character(len=*),                            intent(in)  :: filename
    character(len=*),                            intent(in)  :: label
    real(dp),                                    intent(in)  :: threshold
    integer(i4b),                                intent(out) :: n
    real(dp),         allocatable, dimension(:), intent(out) :: nu, tau
    real(dp),         allocatable, dimension(:)              :: um, tau_Um

    integer(i4b)        :: unit, first, last, m, ierr, l, ext(1), i
    logical(lgt)        :: exist
    character(len=128)  :: string
    type(hdf_file)     :: file
    real(dp), allocatable, dimension(:) :: x, y

    unit = getlun()
    
    inquire(file=trim(filename), exist=exist)
    if (.not. exist) call report_error('Bandpass file does not exist = ' // trim(filename))

    l = len(trim(filename))

    if (filename(l-2:l) == '.h5' .or. filename(l-3:l) == '.hd5') then

       call open_hdf_file(filename, file, "r")
       call get_size_hdf(file, trim(label) // "/bandpass", ext)
       m = ext(1)

       allocate(x(m), y(m))
       !write(*,*) "About to read bandpass"
       call read_hdf(file, trim(label) // "/bandpassx",x)
       call read_hdf(file, trim(label) // "/bandpass", y)
       call close_hdf_file(file)

       ! Drop double entries
       l = 1
       do while (l < m)
          if (x(l) == x(l+1)) then
             x(l:m-1) = x(l+1:m)
             y(l:m-1) = y(l+1:m)
             m        = m-1
          else
             l = l+1
          end if
       end do
       
    else 
       ! Assume ASCII

       ! Find the number of entries
       m = 0
       open(unit, file=trim(filename))
       do while (.true.)
          read(unit,*,end=1) string
          if (string(1:1)=='#') cycle
          m = m + 1
       end do
1      close(unit)
       
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
2      close(unit)
    end if

    ! x(1:m) = x(1:m) * 1.d9 ! Convert from GHz to Hz

    first = 1
    last  = m
    if (threshold > 0.d0) then
       do while (y(first) < threshold*maxval(y(1:m)))
          first = first+1
       end do
       do while (y(last) < threshold*maxval(y(1:m)))
          last = last-1
       end do
    end if
    
    n = last-first+1
    allocate(nu(n), tau(n))
    allocate(um(n), tau_um(n))
    um     = x(first:last)
    tau_um = y(first:last)

    deallocate(x, y)

    ! Seems that according to the DIRBE team, the quoted spectral response
    ! does not get converted when doing color corrections, so we will keep
    ! tau = tau_um
    !
    ! Ref: https://lambda.gsfc.nasa.gov/data/cobe/dirbe/color_correction/dirbe_colorcorr.pro

    do i = 1, n
       nu(i)  = 2.99792458d14/um(n-i+1) ! Convert from micron to Hz
       tau(i) = tau_um(n-i+1)
    end do

  end subroutine read_bandpass_dirbe


  function micron_to_GHz(um)
    implicit none

    real(dp), intent(in) :: um
    real(dp) :: c_um_GHz
    real(dp) :: micron_to_GHz

    c_um_GHz = 2.99792458d5

    micron_to_GHz = um**2/c_um_GHz

  end function micron_to_GHz


  ! Routine for reading bandpass files all detectors at the same time
  subroutine read_bandpass_nonzero_threshold(filename, label, ndet, threshold, n, nu, tau)
    implicit none

    character(len=*),                            intent(in)  :: filename
    character(len=*),              dimension(:), intent(in)  :: label
    integer(i4b),                                intent(in)  :: ndet
    real(dp),                                    intent(in)  :: threshold
    integer(i4b),                                intent(out) :: n
    real(dp),         allocatable, dimension(:), intent(out) :: nu, tau

    integer(i4b)        :: unit, first, last, m, ierr, l, ext(1), i
    logical(lgt)        :: exist
    character(len=128)  :: string
    type(hdf_file)     :: file
    real(dp), allocatable, dimension(:) :: x, y
    integer(i4b), allocatable, dimension(:) :: first_array, last_array, m_array

    unit = getlun()
    
    inquire(file=trim(filename), exist=exist)
    if (.not. exist) call report_error('Bandpass file does not exist = ' // trim(filename))

    l = len(trim(filename))
    allocate(m_array(ndet), first_array(ndet), last_array(ndet))

    if (filename(l-2:l) == '.h5' .or. filename(l-3:l) == '.hd5') then

       call open_hdf_file(filename, file, "r")

       do i = 1, ndet      
           call get_size_hdf(file, trim(label(i)) // "/bandpass", ext)
           m = ext(1)

           allocate(x(m), y(m))
           !write(*,*) "About to read bandpass"
           call read_hdf(file, trim(label(i)) // "/bandpassx",x)
           call read_hdf(file, trim(label(i)) // "/bandpass", y)
           !call close_hdf_file(file)

           ! Drop double entries
           l = 1
           do while (l < m)
               if (x(l) == x(l+1)) then
                   x(l:m-1) = x(l+1:m)
                   y(l:m-1) = y(l+1:m)
                   m        = m-1
               else
                   l = l+1
               end if
           end do
           m_array(i) = m ! it only works because m is the same for all det in HFI
           deallocate(x, y)
       end do
       
    end if

    do i = 1, ndet
        allocate(x(m_array(i)), y(m_array(i)))
        call read_hdf(file, trim(label(i)) // "/bandpassx",x)
        call read_hdf(file, trim(label(i)) // "/bandpass", y)
        x = x * 1.d9 ! Convert from GHz to Hz

        first = 1
        last  = m_array(i)
        if (threshold > 0.d0) then
            do while (y(first) < threshold*maxval(y(1:m)))
                first = first+1
            end do
            do while (y(last) < threshold*maxval(y(1:m)))
                last = last-1
            end do
        end if
        first_array(i) = first
        last_array(i) = last
        deallocate(x, y)
    end do

    ! choosing the range first:last suited for all detectors 
    first = MAXVAL(first_array)
    last = MINVAL(last_array)

    
    n = last-first+1
    allocate(nu(n), tau(n))

    do i = 1, ndet
        allocate(x(m_array(i)), y(m_array(i)))
        call read_hdf(file, trim(label(i)) // "/bandpassx",x)
        call read_hdf(file, trim(label(i)) // "/bandpass", y)
        x = x * 1.d9 ! Convert from GHz to Hz

        if (i == 1) then
            nu  = x(first:last)
            tau = y(first:last)
        else
            nu  = nu  + x(first:last)
            tau = tau + y(first:last)
        end if
        deallocate(x, y)
    end do

    nu  = nu  / ndet
    tau = tau / ndet

    deallocate(m_array, first_array, last_array)

  end subroutine read_bandpass_nonzero_threshold


end module comm_bp_utils
