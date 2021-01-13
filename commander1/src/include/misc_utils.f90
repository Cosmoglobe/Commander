!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon, 
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------
module misc_utils
!   subroutine fatal_error
!   function file_present
!   subroutine assert_present
!   subroutine assert_directory_present
!   subroutine assert_not_present
!   subroutine assert_alloc
!   subroutine assert
!   function strupcase
!   function strlowcase
!   function string
!   subroutine wall_clock_time
!   subroutine brag_openmp

  use healpix_types
  use extension, only : exit_with_status
  implicit none
  private

  integer, parameter, private :: LCH=48
  interface string
     module procedure string_i, string_s, string_d
  end interface

  public :: fatal_error, assert, assert_present, assert_not_present, &
       & assert_alloc, file_present, assert_directory_present
  public :: upcase, lowcase
  public :: wall_clock_time
  public :: brag_openmp
  public :: strupcase, strlowcase, string

contains
  !-----------------------------------------------------
  subroutine fatal_error (msg)
    character(len=*), intent(in), optional :: msg

    if (present(msg)) then
       print *,'Fatal error: ', trim(msg)
    else
       print *,'Fatal error'
    endif
    call exit_with_status(1)
  end subroutine fatal_error

  !-----------------------------------------------------
  function file_present (filename)
    character(len=*), intent(in) :: filename
    logical :: file_present

    inquire(file=trim(filename),exist=file_present)
  end function file_present

  !-----------------------------------------------------
  subroutine assert_present (filename)
    character(len=*), intent(in) :: filename

    if (.not. file_present(trim(filename))) then
       print *, 'Error:  file ' // trim(filename) // ' does not exist!'
       call exit_with_status(1)
    end if
  end subroutine assert_present

  !-----------------------------------------------------
  subroutine assert_directory_present (filename)
    character(len=*), intent(in) :: filename
    integer pos

    pos = scan(filename,'/',.true.)
    if (pos<=0) return

    if (.not. file_present(filename(:pos-1))) then
       print *, 'Error:  directory ' // filename(:pos-1) // ' does not exist!'
       call exit_with_status(1)
    end if
  end subroutine assert_directory_present

  !-----------------------------------------------------
  subroutine assert_not_present (filename)
    character(len=*), intent(in) :: filename

    if (file_present(trim(filename))) then
       print *, 'Error:  file ' // trim(filename) // ' already exists!'
       call exit_with_status(1)
    end if
  end subroutine assert_not_present

  !-----------------------------------------------------
  subroutine assert_alloc (stat,code,arr)
    integer, intent(in) :: stat
    character(len=*), intent(in) :: code, arr

    if (stat==0) return

    print *, trim(code)//'> cannot allocate memory for array: '//trim(arr)
    call exit_with_status(1)
  end subroutine assert_alloc

  !-----------------------------------------------------
  subroutine assert (testval,msg,errcode)
    logical, intent(in) :: testval
    character(len=*), intent(in), optional :: msg
    integer, intent(in), optional :: errcode

    if (testval) return

    print *,"Assertion failed: "
    if (present(msg)) print *, trim(msg)
    if (present(errcode)) call exit_with_status (errcode)
    call exit_with_status(1)
  end subroutine assert

  !-----------------------------------------------------
  subroutine upcase(instr, outstr)
    ! turns a string to upper case
    ! instr and outstr can be the same variable
    character(len=*), intent(in)  :: instr
    character(len=*), intent(out) :: outstr

    integer(i4b) :: i, j, ascii, ll, la, ua

    la = iachar('a')
    ua = iachar('A')

    ll = len_trim(outstr)
    do i = 1, min(len_trim(instr),ll)
       ascii = iachar( instr(i:i) )
       if (ascii >= la .and. ascii < la+26) then ! in [a,z]
          outstr(i:i) = achar( ascii - la + ua )
       else
          outstr(i:i) = instr(i:i)
       endif
    enddo
    do j = i, ll ! pad with blanks
       outstr(j:j) = ' '
    enddo
    return
  end subroutine upcase

  !-----------------------------------------------------
  subroutine lowcase(instr, outstr)
    ! turns a string to lower case
    ! instr and outstr can be the same variable
    character(len=*), intent(in)  :: instr
    character(len=*), intent(out) :: outstr

    integer(i4b) :: i, j, ascii, ll, la, ua

    la = iachar('a')
    ua = iachar('A')

    ll = len_trim(outstr)
    do i = 1, min(len_trim(instr),ll)
       ascii = iachar( instr(i:i) )
       if (ascii >= ua .and. ascii < ua+26) then ! in [A,Z]
          outstr(i:i) = achar( ascii - ua + la )
       else
          outstr(i:i) = instr(i:i)
       endif
    enddo
    do j = i, ll ! pad with blanks
       outstr(j:j) = ' '
    enddo
    return
  end subroutine lowcase
  !-----------------------------------------------------
  function strupcase(instr) result(outstr)
    ! turns a character string to upper case
    character(len=*), intent(in)  :: instr
    character(len=FILENAMELEN)    :: outstr

    integer(i4b) :: i, ascii, la, ua

    la = iachar('a')
    ua = iachar('A')

    outstr = instr
    do i = 1, min(len_trim(instr),len_trim(outstr))
       ascii = iachar( instr(i:i) )
       if (ascii >= la .and. ascii < la+26) then ! in [a,z]
          outstr(i:i) = achar( ascii - la + ua )
       endif
    enddo

    return
  end function strupcase

  !-----------------------------------------------------
  function strlowcase(instr) result(outstr)
    ! turns a string to lower case
    character(len=*), intent(in)  :: instr
    character(len=FILENAMELEN)    :: outstr

    integer(i4b) :: i, ascii, la, ua

    la = iachar('a')
    ua = iachar('A')

    outstr = instr
    do i = 1, min(len_trim(instr),len_trim(outstr))
       ascii = iachar( instr(i:i) )
       if (ascii >= ua .and. ascii < ua+26) then ! in [A,Z]
          outstr(i:i) = achar( ascii - ua + la )
       endif
    enddo

    return
  end function strlowcase

  !========================================
  ! function string(arg, format)
  !=======================================
  function string_i(arg, format) result(str)
    integer(i4b) :: arg
    character(len=*),   optional :: format
    character(len=LCH)           :: str
    if (present(format)) then
       write(str,format) arg
    else
       write(str,*) arg
    endif
    return
  end function string_i
  !--------------------------------
  function string_s(arg, format) result(str)
    real(sp) :: arg
    character(len=*),   optional :: format
    character(len=LCH)           :: str
    if (present(format)) then
       write(str,format) arg
    else
       write(str,*) arg
    endif
    return
  end function string_s
  !--------------------------------
  function string_d(arg, format) result(str)
    real(dp) :: arg
    character(len=*),   optional :: format
    character(len=LCH)           :: str
    if (present(format)) then
       write(str,format) arg
    else
       write(str,*) arg
    endif
    return
  end function string_d
  !--------------------------------

  !-----------------------------------------------------
  subroutine wall_clock_time(time_sec)
    real(sp), intent(out) :: time_sec

    integer :: clock, clock_rate, clock_max
    integer, dimension(8) :: values_time

    time_sec = 0.

    call system_clock(count=clock, count_rate=clock_rate, count_max=clock_max)

    if (clock < 0 .or. clock_rate <= 0 .or. clock_max <= 0) then
       call date_and_time(values = values_time) ! y, m, d, x, h, m, s, ms
       if (minval(values_time) >= 0) then
          time_sec = ((values_time(3)*24. &
               &     + values_time(5)    )*60. &
               &     + values_time(6)        )*60. &
               &     + values_time(7) + values_time(8)/1000.
       endif
    else
       time_sec = clock/real(clock_rate)
    endif

    return
  end subroutine wall_clock_time
  !================================================
  subroutine brag_openmp()
  !================================================
    ! OpenMP bragging
    !================================================
    ! OpenMP variables
    !$     integer :: omp_get_thread_num, omp_get_num_threads, omp_get_num_procs
    !IBMP  integer :: omp_get_thread_num, omp_get_num_threads, omp_get_num_procs

!$OMP parallel
!
!$   if (omp_get_thread_num() == 0) then
!$       write(*,9000) ' --------------------------------------'
!$       write(*,9010) ' Number of OpenMP threads in use: ', omp_get_num_threads()
!$       write(*,9010) ' Number of CPUs available:        ', omp_get_num_procs()
!$       write(*,9000) ' --------------------------------------'
!$   end if
!
!IBMP   if (omp_get_thread_num() == 0) then
!IBMP       write(*,9000) ' --------------------------------------'
!IBMP       write(*,9010) ' Number of OpenMP threads in use: ', omp_get_num_threads()
!IBMP       write(*,9010) ' Number of CPUs available:        ', omp_get_num_procs()
!IBMP       write(*,9000) ' --------------------------------------'
!IBMP   end if
!
!$OMP end parallel
9000 format(a)
9010 format(a,i4)

    return
  end subroutine brag_openmp

end module misc_utils
