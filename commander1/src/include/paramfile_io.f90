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
! -*- f90 -*-

!
! v1.0: M. Reinecke
! v1.1: 2002-09, E. Hivon, added concatnl, scan_directories, numeric_string
!                     made interactive mode more user friendly
!
module paramfile_io
  use healpix_types
  use extension
  use misc_utils
  implicit none
  private

  public paramfile_handle, parse_init, parse_real, parse_double, parse_int, &
         parse_long, parse_lgt, parse_string, parse_summarize, parse_finish
  public concatnl, scan_directories

  type paramfile_handle
    character(len=filenamelen) filename
    character(len=filenamelen), pointer, dimension(:) :: keylist=>NULL()
    character(len=filenamelen), pointer, dimension(:) :: valuelist=>NULL()
    logical interactive
  end type paramfile_handle

  character(len=*), parameter, public :: ret = achar(10)//' '
  character(len=*), parameter, private :: swdef = ' <default>'

contains

!=====================================================================
subroutine notify_user (keyname, rdef, rmin, rmax, ddef, dmin, dmax, &
  idef, imin, imax, ldef, lmin, lmax, logdef, chdef, descr)
  !=====================================================================
  ! prompts user for next parameter when in interactive mode
  !=====================================================================
  character(len=*), intent(in) :: keyname
  real(sp), intent(in), optional :: rdef, rmin, rmax
  real(dp), intent(in), optional :: ddef, dmin, dmax
  integer(i4b), intent(in), optional :: idef, imin, imax
  integer(i8b), intent(in), optional :: ldef, lmin, lmax
  logical, intent(in), optional :: logdef
  character(len=*), intent(in), optional :: chdef, descr

  if (present(descr)) then
     write(*,'(a)') trim(descr)
  else
     print *, 'Please enter a value for the key ', keyname
  endif
  if (present(rmin) .and. present(rmax)) then
     print *, "allowed range: ", rmin, rmax
  else
     if (present(rmin)) print *, "min value: ", rmin
     if (present(rmax)) print *, "max value: ", rmax
  endif
  if (present(dmin) .and. present(dmax)) then
     print *, "allowed range: ", dmin, dmax
  else
     if (present(dmin)) print *, "min value: ", dmin
     if (present(dmax)) print *, "max value: ", dmax
  endif
  if (present(imin) .and. present(imax)) then
     print *, "allowed range: ", imin, imax
  else
     if (present(imin)) print *, "min value: ", imin
     if (present(imax)) print *, "max value: ", imax
  endif
  if (present(lmin) .and. present(lmax)) then
     print *, "allowed range: ", lmin, lmax
  else
     if (present(lmin)) print *, "min value: ", lmin
     if (present(lmax)) print *, "max value: ", lmax
  endif
  if (present(rdef)) print *, "default value: ", rdef
  if (present(ddef)) print *, "default value: ", ddef
  if (present(idef)) print *, "default value: ", idef
  if (present(ldef)) print *, "default value: ", ldef
  if (present(logdef)) print *, "default value: ", logdef
  if (present(chdef)) print *, "default value: ", trim(chdef)
end subroutine notify_user

!===================================================================
function parse_init (fname)
  !===================================================================
  character(len=*), intent(in) :: fname
  type(paramfile_handle) parse_init
  integer  :: i,cnt
  character(len=filenamelen) :: line, name, value

  if (len(trim(fname))==0) then
    parse_init%filename = ''
    parse_init%interactive=.true.
    parse_init%keylist     => NULL()
    parse_init%valuelist   => NULL()
    cnt = 30
    allocate(parse_init%keylist(cnt),parse_init%valuelist(cnt))
    parse_init%keylist = ''
    parse_init%valuelist = ''
  else
    call assert_present (fname)
    call assert(len(fname)<=filenamelen, 'Parser: error: file name too long')
    parse_init%filename    = fname
    parse_init%interactive = .false.
    parse_init%keylist     => NULL()
    parse_init%valuelist   => NULL()
    ! count valid lines
    open  (1, file=trim(fname))
    cnt=0
    do
      read (1,'(a)',end=2) line
      line = adjustl(line)
      i=scan(line,'=')
      if (i/=0 .and. line(1:1)/='#' .and. line(1:1)/='!') cnt=cnt+1
    end do
2   close (1)
    ! read and parse valid lines
    allocate(parse_init%keylist(cnt),parse_init%valuelist(cnt))
    open  (1, file=trim(fname))
    cnt=0
    do
      read (1,'(a)',end=3) line
      line = adjustl(line)
      i=scan(line,'=')
      if (i/=0 .and. line(1:1)/='#' .and. line(1:1)/='!') then
        cnt=cnt+1
        name = trim(adjustl(line(:i-1)))
        value = trim(adjustl(line(i+1:)))
        if (trim(value)=="") then
           write(*,'(a)') ' '
           write(*,'(a)') 'ERROR: Inputs of the form '
           write(*,'(a)') trim(name)//' = '
           write(*,'(a)') ' (ie, defined as a blank value) are not valid'
           write(*,'(a)') 'To get the default value, comment out the keyword in '&
                &         //trim(parse_init%filename)
           write(*,'(a)') '# '//trim(name)//' = '
           write(*,'(a)') "If you mean 'No file', use"
           write(*,'(a)') trim(name)//" = ''   "
           write(*,'(a)') ' '
           call fatal_error
        endif
        parse_init%keylist(cnt) = name 
        parse_init%valuelist(cnt) = value
      endif
    end do
3   close (1)
  endif

  ! be verbose
  if (parse_init%interactive) then
     write(*,'(a)') 'Interactive mode. Answer the following questions.'
     write(*,'(a)') 'If no answer is entered, the default value will be taken'
  else
     write(*,'(a)') 'Reading run parameters from '//trim(parse_init%filename)
     write(*,'(a)') ' parameters not defined in that file will be set to their default value'
  endif
end function parse_init
!===================================================================
subroutine parse_summarize (handle, code, prec)
  !===================================================================
  type(paramfile_handle),     intent(in) :: handle
  character(len=*), optional, intent(in) :: code
  integer(i4b),     optional, intent(in) :: prec
  !
  integer(i4b) :: i
  character(len=filenamelen) :: name, value, next_name, command

  if (handle%interactive) then
     command = ''
     if (present(code)) then
        command = trim(code)
        if (present(prec)) then
           if (prec == SP) command = trim(command)//' --single'
           if (prec == DP) command = trim(command)//' --double'
        endif
     endif
     if (trim(command) /= '') then
        print*,'  This run can be reproduced in non-interactive mode, with the command'
        print*,trim(command)//' paramfile'
        print*,'where paramfile contains'
     else
        print*,'  This run can be reproduced in non-interactive mode'
        print*,'if a parameter file with the following content is provided:'
     endif
     do i=1,size(handle%keylist)
        name = handle%keylist(i)
        next_name = handle%keylist(i+1)
        value = handle%valuelist(i)
        if (trim(name) /= '' .and. trim(name) /= trim(next_name)) then 
           if (trim(value) == '') then
              write(*,'(a)') '# '//trim(name)
           else
              write(*,'(a)') trim(name)//' = '//trim(value)
           endif
        endif
     enddo
     print*,' '
  endif
end subroutine parse_summarize
!===================================================================
subroutine parse_finish (handle)
  !===================================================================
  type(paramfile_handle), intent(inout) :: handle

  if (associated(handle%keylist)) &
    deallocate(handle%keylist, handle%valuelist)
end subroutine parse_finish

!===================================================================
subroutine find_param (handle,keyname,result,found,rdef,rmin,rmax, &
    ddef,dmin,dmax,idef,imin,imax,ldef,lmin,lmax,logdef,chdef,descr)
  !===================================================================
  ! extract parameter from file or read from standard input
  !===================================================================
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  character(len=*), intent(out) :: result
  logical, intent(out) :: found
  real(sp), intent(in), optional :: rdef, rmin, rmax
  real(dp), intent(in), optional :: ddef, dmin, dmax
  integer(i4b), intent(in), optional :: idef, imin, imax
  integer(i8b), intent(in), optional :: ldef, lmin, lmax
  logical, intent(in), optional :: logdef
  character(len=*), intent(in), optional :: chdef, descr

  character(len=filenamelen) :: line, name, value
  integer i
  !===================================================================

  found=.false.

  if (handle%interactive) then
     call notify_user (keyname,rdef,rmin,rmax,ddef,dmin,dmax, &
          &            idef,imin,imax,ldef,lmin,lmax,logdef,chdef,descr)
     read (*,'(a)',err=5) result
     found = (trim(result)/='')
     do i=1,size(handle%keylist)
        if (trim(handle%keylist(i))=='') then
           handle%keylist(i)   = trim(keyname)
           if (found) then
              handle%valuelist(i) = trim(result)
           else
              if (present(rdef)) write(handle%valuelist(i),*) rdef
              if (present(ddef)) write(handle%valuelist(i),*) ddef
              if (present(idef)) write(handle%valuelist(i),*) idef
              if (present(ldef)) write(handle%valuelist(i),*) ldef
              if (present(logdef)) write(handle%valuelist(i),*) logdef
              if (present(chdef))  handle%valuelist(i) = chdef
           endif
           exit
        end if
     end do
  else
     do i=1,size(handle%keylist)
        if (trim(handle%keylist(i))==keyname) then
           result=trim(handle%valuelist(i))
           found=.true.
        end if
     end do
2    close (1)
  endif
  return

5 print*,'Parser: find_param: error reading value'
  call fatal_error
end subroutine find_param
!===================================================================

!===================================================================
function parse_real (handle, keyname, default, vmin, vmax, descr)
  !===================================================================
  !===================================================================
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  real(sp), intent(in), optional :: default, vmin, vmax
  character(len=*), intent(in), optional :: descr
  real(sp) :: parse_real

  character(len=filenamelen) :: result
  character(len=30)          :: about_def
  logical found
  !===================================================================

10 continue
  about_def = ''
  call find_param (handle, trim(keyname), result, found, rdef=default, &
       &           rmin=vmin, rmax=vmax, descr=descr)
  if (found) then
     read (result,*,err=5) parse_real
  else
     if (present(default)) then
!      print *,'Parser: warning: using default value for ',trim(keyname)
        about_def = swdef
        parse_real = default
     else
        print *,'Parser: error: ',trim(keyname),' not found.'
        goto 2
     endif
  endif
  print *,'Parser: ',trim(keyname),' = ',parse_real, trim(about_def)
  if (present(vmin)) then
     if (parse_real<vmin) then
        print *,'Parser: error: value for ', trim(keyname),' too small.'
        goto 2
     endif
  endif
  if (present(vmax)) then
     if (parse_real>vmax) then
        print *,'Parser: error: value for ', trim(keyname),' too large.'
        goto 2
     endif
  endif
  
  return ! normal exit
  
5 print*,'Parser: parse_real: error reading value'
2 if (handle%interactive) goto 10 ! try again
  call fatal_error
  
end function parse_real

!===================================================================
function parse_double (handle, keyname, default, vmin, vmax, descr)
  !===================================================================
  !===================================================================
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  real(dp), intent(in), optional :: default, vmin, vmax
  character(len=*), intent(in), optional :: descr
  real(dp) :: parse_double

  character(len=filenamelen) :: result
  character(len=30)          :: about_def
  logical found
  !===================================================================

10 continue
  about_def = ''
  call find_param (handle, trim(keyname), result, found, ddef=default, &
       &           dmin=vmin, dmax=vmax, descr=descr)
  if (found) then
    read (result,*,err=5) parse_double
  else
    if (present(default)) then
!      print *,'Parser: warning: using default value for ',trim(keyname)
        about_def = swdef
      parse_double = default
    else
      print *,'Parser: error: ',trim(keyname),' not found.'
      goto 2
    endif
  endif
  print *,'Parser: ',trim(keyname),' = ',parse_double, about_def
  if (present(vmin)) then
    if (parse_double<vmin) then
      print *,'Parser: error: value for ', trim(keyname),' too small.'
      goto 2
    endif
  endif
  if (present(vmax)) then
    if (parse_double>vmax) then
      print *,'Parser: error: value for ', trim(keyname),' too large.'
      goto 2
    endif
  endif

  return ! normal exit

5 print*,'Parser: parse_double: error reading value'
2 if (handle%interactive) goto 10 ! try again
  call fatal_error

end function parse_double

!==================================================================
function parse_int (handle, keyname, default, vmin, vmax, descr)
  !==================================================================
  ! parse 4 byte integer parameter
  !==================================================================
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  integer(i4b), intent(in), optional :: default, vmin, vmax
  character(len=*), intent(in), optional :: descr
  integer(i4b) :: parse_int

  character(len=filenamelen) :: result
  character(len=30)          :: about_def
  logical :: found
  !==================================================================

10 continue
  about_def = ''
  call find_param (handle, trim(keyname), result, found, idef=default, &
       &           imin=vmin, imax=vmax, descr=descr)
  if (found) then
    read (result,*,err=5) parse_int
  else
    if (present(default)) then
!      print *,'Parser: warning: using default value for ',trim(keyname)
        about_def = swdef
      parse_int = default
    else
      print *,'Parser: error: ',trim(keyname),' not found.'
      goto 2
    endif
  endif
  print *,'Parser: ',trim(keyname),' = ',parse_int, trim(about_def)
  if (present(vmin)) then
    if (parse_int<vmin) then
      print *,'Parser: error: value for ', trim(keyname),' too small.'
      goto 2
    endif
  endif
  if (present(vmax)) then
    if (parse_int>vmax) then
      print *,'Parser: error: value for ', trim(keyname),' too large.'
      goto 2
    endif
  endif

  return ! normal exit

5 print*,'Parser: parse_int: error reading value'
2 if (handle%interactive) goto 10 ! try again
  call fatal_error

end function parse_int
!==================================================================

!==================================================================
function parse_long (handle, keyname, default, vmin, vmax, descr)
  !==================================================================
  ! parse 8 byte integer parameter
  !==================================================================
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  integer(i8b), intent(in), optional :: default, vmin, vmax
  character(len=*), intent(in), optional :: descr
  integer(i8b) :: parse_long

  character(len=filenamelen) :: result
  character(len=30)          :: about_def
  logical found
  !==================================================================

10 continue
  about_def = ''
  call find_param (handle, trim(keyname), result, found, ldef=default, &
                   lmin=vmin, lmax=vmax, descr=descr)
  if (found) then
    read (result,*,err=5) parse_long
  else
    if (present(default)) then
!      print *,'Parser: warning: using default value for ',trim(keyname)
        about_def = swdef
      parse_long = default
    else
      print *,'Parser: error: ',trim(keyname),' not found.'
      goto 2
    endif
  endif
  print *,'Parser: ',trim(keyname),' = ',parse_long, trim(about_def)
  if (present(vmin)) then
    if (parse_long<vmin) then
      print *,'Parser: error: value for ', trim(keyname),' too small.'
      goto 2
    endif
  endif
  if (present(vmax)) then
    if (parse_long>vmax) then
      print *,'Parser: error: value for ', trim(keyname),' too large.'
      goto 2
    endif
  endif

  return ! normal exit

5 print*,'Parser: parse_long: error reading value'
2 if (handle%interactive) goto 10 ! try again
  call fatal_error

end function parse_long

!===================================================================
function parse_lgt (handle, keyname, default, descr)
  !===================================================================
  ! parse (1 byte) logical parameter
  !===================================================================
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  logical, intent(in), optional :: default
  character(len=*), intent(in), optional :: descr
  logical :: parse_lgt

  character(len=filenamelen) :: result
  character(len=30)          :: about_def
  logical found
  !===================================================================

10 continue
  about_def = ''
  call find_param (handle, trim(keyname), result, found, logdef=default, &
       &           descr=descr)
  if (found) then
     select case (strupcase(result))
     case ('Y','YES','T','TRUE', '.TRUE.','1')
        parse_lgt = .true.
     case ('N','NO', 'F','FALSE','.FALSE.','0')
        parse_lgt= .false.
     case default
        goto 5
     end select
  else
    if (present(default)) then
!       print *,'Parser: warning: using default value for ',trim(keyname)
      parse_lgt = default
    else
      print *,'Parser: error: ',trim(keyname),' not found.'
      goto 2
    endif
  endif
  print *,'Parser: ',trim(keyname),' = ',parse_lgt, trim(about_def)

  return ! normal exit

5 print*,'Parser: parse_lgt: error reading value'
2 if (handle%interactive) goto 10 ! try again
  call fatal_error

end function parse_lgt

!===================================================================
function parse_string (handle, keyname, default, descr, filestatus, options)
  !===================================================================
  ! parse a character string parameter
  !
  ! if filestatus is 'old', look for an existing file having the name of the string
  !
  ! if filestatus is 'new', no file with the exact same name as the string should exist
  !
  ! options is the list of valid options
  !
  !===================================================================
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  character(len=*), intent(in), optional :: default
  character(len=*), intent(in), optional :: descr
  character(len=*), intent(in), optional :: filestatus
  character(len=*), intent(in), optional, dimension(1:) :: options

  character(len=filenamelen) :: parse_string

  character(len=filenamelen) :: result
  character(len=30)          :: about_def
  logical :: found, there
  integer :: i
  !===================================================================

10 continue
  about_def = ''
  call find_param (handle, trim(keyname), result, found, chdef=default, &
                   descr=descr)
  if (found) then
    parse_string = trim(result)
  else
    if (present(default)) then
!       write(*,'(1x,a)') 'Parser: warning: using default value for '//trim(keyname)
        about_def = swdef
      parse_string = trim(default)
    else
      write(*,'(1x,a)') 'Parser: error: '//trim(keyname)//' not found.'
      goto 2
    endif
  endif
  write(*,'(1x,a)') 'Parser: '//trim(keyname)//' = '//trim(parse_string)//trim(about_def)

  ! 0 (zero), '' and ' ' (2 single quotes with nothing or one space in between)
  !         are interpreted as "No File"
  if (trim(adjustl(parse_string)) == "0" )  parse_string = ''
  if (trim(adjustl(parse_string)) == "''")  parse_string = ''
  if (trim(adjustl(parse_string)) == "' '") parse_string = ''

  if (present(filestatus) .and. trim(parse_string) /= '') then
     if (trim(filestatus)=='new' .or. trim(filestatus)=='NEW') then
        inquire(file=trim(parse_string),exist=there)
        if (there) then
           print *, 'Parser: error: output file ' // trim(parse_string) // &
                ' already exists!'
           goto 2
        end if
     else if (trim(filestatus)=='old' .or. trim(filestatus)=='OLD') then
        inquire(file=trim(parse_string),exist=there)
        if (.not. there) then
           print *, 'Parser: error: input file ' // trim(parse_string) // &
                ' does not exist!'
           goto 2
        end if
     else
        print *, 'Parser: error: wrong value for filestatus :',filestatus
        call fatal_error
     endif
  endif

  if (present(options)) then
     do i=1, size(options) 
        if (trim(adjustl(parse_string)) == trim(adjustl(options(i)))) goto 5
     enddo
     print*,'Invalid choice'
     goto 2
5   continue
  endif

  return ! normal exit

2 if (handle%interactive) goto 10 ! try again
  call fatal_error

end function parse_string



!========================================================================
function concatnl(line1,line2,line3,line4,line5,line6,line7,line8,line9,line10)
  !========================================================================
  ! concatenate line1, line2, line3,... into one string,
  ! while putting a char(10) Line Feed in between
  !========================================================================

  character(len=*), intent(in)           :: line1
  character(len=*), intent(in), optional ::       line2,line3,line4,line5
  character(len=*), intent(in), optional :: line6,line7,line8,line9,line10

  character(len=filenamelen) :: concatnl

  concatnl = trim(line1)
  if (present(line2)) concatnl = trim(concatnl)//ret//trim(line2)
  if (present(line3)) concatnl = trim(concatnl)//ret//trim(line3)
  if (present(line4)) concatnl = trim(concatnl)//ret//trim(line4)
  if (present(line5)) concatnl = trim(concatnl)//ret//trim(line5)
  if (present(line6)) concatnl = trim(concatnl)//ret//trim(line6)
  if (present(line7)) concatnl = trim(concatnl)//ret//trim(line7)
  if (present(line8)) concatnl = trim(concatnl)//ret//trim(line8)
  if (present(line9)) concatnl = trim(concatnl)//ret//trim(line9)
  if (present(line10)) concatnl = trim(concatnl)//ret//trim(line10)


end function concatnl
!========================================================================

!========================================================================
function scan_directories(directories, filename, fullpath)
  !========================================================================
  ! scan directories in search of filename,
  ! if found, returns .true. and the full path is in fullpath.
  ! The search is *NOT* reccursive
  !
  ! it assumes that the given directory and filename are separated by either
  !  nothing, a / (slash) or a \ (backslash)
  !
  ! if several directories are to be searched (up to 20),
  ! concatenate them into 'directories',
  ! putting a special character (ASCII < 32) between them.
  !  see concatnl
  ! NB: a space is not a special character
  !========================================================================
  logical(LGT) :: scan_directories
  character(len=*), intent(in)  :: filename, directories
  character(len=*), intent(out) :: fullpath

  logical :: found
  integer(I4B), dimension(1:20) :: index
  integer(I4B)                  :: i, k, nc, nspecial
  character(len=1)              :: ch
  character(len=filenamelen)    :: directory
  character(len=3000)           :: string
  character(LEN=1), DIMENSION(1:3) :: separator
  character(len=*), parameter   :: code = 'scan_directories'
  !========================================================================

  ! define separators (this is the only way that works for all compilers)
  separator(1) = char(32) ! ' '
  separator(2) = char(47) ! '/'
  separator(3) = char(92) ! '\'

  ! find location of special characters
  nc = len_trim(directories)
  index(1) = 0
  nspecial = 2
  do i=1,nc
     ch = directories(i:i)
     if (iachar(ch) < 32) then
        index(nspecial) = i
        nspecial        = nspecial + 1
     endif
  enddo
  index(nspecial) = nc + 1

  ! test string between special character as potential directory
  fullpath = ''
  found = .false.
  do i = 1, nspecial-1
     directory=trim(adjustl(directories(index(i)+1:index(i+1)-1)))
        do k = 1, size(separator)
           string = trim(directory)//trim(separator(k))//trim(filename)
           inquire(&
                &  file=string, &
                &  exist=found)
           if (found) goto 10
        enddo
  enddo

10 continue
  if (found) then
     if (len(fullpath) >= len_trim(string)) then
        fullpath = trim(string)
     else
        print*,code
        print*,'variable fullpath is not large enough'
        print*,'requires ',len_trim(string),' characters'
        print*,'has only ',trim(fullpath)
        call fatal_error
     endif
  endif

  scan_directories = found

end function scan_directories

end module paramfile_io
