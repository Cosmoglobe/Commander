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

module obsolete

  ! ==============================================================
  ! contains routines that were initially in
  ! module fitstools (routines ask_*)
  !     and
  ! module utilities (*_parser, *_setpar, *_getpar ... and related)
  !
  ! and have become obsolete in version 1.2
  !
  ! They won't be supported any more
  ! See module paramfile_io for substitutes
  !===============================================================

  ! Authors: Banday and Bartelmann
  ! CHANGES:
  ! Changed 9/MAR/1999 by B.D.Wandelt
  !     Incorporated non-interactive functionality for anafast.
  !     Added: type anapar, subroutines anafast_parser, anafast_setpar and anafast_getpar.
  !     Incorporated non-interactive functionality for HotSpots.
  !     Added: type hotpar, subroutines hotspot_parser, hotspots_setpar and hotspots_getpar.
  !     Incorporated non-interactive functionality for ud_grade.
  !     Added: type udgpar, subroutines udgrade_parser, udgrade_setpar and udgrade_getpar.
  ! Changed 9/MAY/1999 by B.D.Wandelt
  ! Changed all routines to new, more user-friendly parser,
  ! provided by Bartelmann and Banday.
  ! Added parser for smoothing routine, E. Hivon, 2001-11
  !     Incorporated non-interactive functionality for smoothing.
  !     Added: type udgpar, subroutines smoothing_parser, smoothing_setpar and smoothing_getpar.
  !

  use healpix_types
  implicit none

  type synpar
     integer(i4b) :: simtp
     integer(i4b) :: nsmax
     integer(i4b) :: nlmax
     integer(i4b) :: iseed
     real(sp) :: fwhma
     character(len=filenamelen) :: inpfl
     character(len=filenamelen) :: outfl
     character(len=filenamelen) :: almsfile
     character(len=filenamelen) :: plmfile
     character(len=filenamelen) :: windowfile
     character(len=filenamelen) :: winfiledir
  end type synpar

  type(synpar) :: par

  type anapar
     integer(i4b) :: simtp
     integer(i4b) :: nlmax
     real(dp) :: theta_cut_deg
     integer(i4b) :: iter_order
     integer(i4b) :: regression
     character(len=filenamelen) :: w8file
     character(len=filenamelen) :: inpfl
     character(len=filenamelen) :: outfl
     character(len=filenamelen) :: plmfile
     character(len=filenamelen) :: outfl_alms
     character(len=filenamelen) :: w8filedir
     integer(i4b) :: won
  end type anapar

  type(anapar) :: anafast_par

  type hotpar
     character(len=filenamelen) :: inpfl
     character(len=filenamelen) :: outfl_extrema
     character(len=filenamelen) :: outfl_max
     character(len=filenamelen) :: outfl_min
  end type hotpar

  type(hotpar) :: hotspots_par

  type udgpar
     integer(i4b) :: nside_out
     character(len=filenamelen) :: inpfl
     character(len=filenamelen) :: outfl
  end type udgpar

  type(udgpar) :: udgrade_par

  type smopar
     integer(i4b) :: simtp
     integer(i4b) :: iter_order
     integer(i4b) :: nlmax
     real(sp)     :: fwhma
     integer(i4b) :: won
     character(len=filenamelen) :: plmfile
     character(len=filenamelen) :: w8file
     character(len=filenamelen) :: w8filedir
     character(len=filenamelen) :: inpfl
     character(len=filenamelen) :: outfl
  end type smopar

  type(smopar) :: smoothing_par

  private

  public :: ask_inputmap, ask_outputmap, ask_lrange

  public :: synpar, par
  public :: parser, setpar, getpar
  public :: anapar, anafast_par
  public :: anafast_parser, anafast_setpar, anafast_getpar
  public :: hotpar, hotspots_par
  public :: hotspots_parser, hotspots_setpar, hotspots_getpar
  public :: udgpar, udgrade_par
  public :: udgrade_parser, udgrade_setpar, udgrade_getpar
  public :: smopar, smoothing_par
  public :: smoothing_parser, smoothing_setpar, smoothing_getpar


contains
  !========================================================================
  ! subroutine ask_inputmap
  ! subroutine ask_outputmap
  ! subroutine ask_lrange
  !========================================================================
  subroutine ask_inputmap(code, inputfile)
    !========================================================================
    use healpix_types
    CHARACTER(LEN=*), intent(in) :: code
    CHARACTER(LEN=*), intent(out) :: inputfile
    LOGICAL(LGT) :: ok
    !========================================================================

3   PRINT *,' Enter input file name (map, eg: map.fits): '
    READ 9000, inputfile

    inquire(file=inputfile, exist=ok)
    if (.not.ok) then
       PRINT *,' '//code//'> '//TRIM(inputfile)//' not found'
       goto 3
    endif

9000 FORMAT(A)

    return
  end subroutine ask_inputmap
  !========================================================================
  subroutine ask_outputmap(code, outputfile)
    !========================================================================
    use healpix_types
    CHARACTER(LEN=*), intent(in) :: code
    CHARACTER(LEN=*), intent(out) :: outputfile
    LOGICAL(LGT) :: bad
    !========================================================================

2   PRINT *,' Enter Output map file name (eg, map_smoothed.fits) :'
    READ 9000, outputfile

    inquire(file=outputfile, exist=bad)
    if ( bad ) then
       PRINT *,' '//code//'> '//TRIM(outputfile)//' already exists'
       PRINT *,' '//code//'> choose a new output file name.'
       goto 2
    endif
    PRINT *,''

9000 FORMAT(A)

    return
  end subroutine ask_outputmap
  !========================================================================
  subroutine ask_lrange(code, nsmax, nlmax)
    !========================================================================
    use healpix_types
    CHARACTER(LEN=*), INTENT(IN) :: code
    INTEGER(I4B), INTENT(IN) :: nsmax
    INTEGER(I4B), INTENT(OUT) :: nlmax
    !========================================================================

    PRINT *,' Enter the maximum l range for the analysis: '
    WRITE(*,'(a,i5)') '  The map has Nside = ',nsmax
    WRITE(*,'(a,i5,a)') '  (0 <= l <= l_max <= ',2*nsmax,') l_max = '
    READ *, nlmax

    return
  end subroutine ask_lrange


!=======================================================================
!=======================================================================
!=======================================================================
!   subroutine setpar
!   subroutine getpar
!   subroutine anafast_parser
!   subroutine anafast_setpar
!   subroutine anafast_getpar
!   subroutine hotspots_parser
!   subroutine hotspots_setpar
!   subroutine hotspots_getpar
!   subroutine udgrade_parser
!   subroutine udgrade_setpar
!   subroutine udgrade_getpar
!   subroutine smoothing_parser
!   subroutine smoothing_setpar
!   subroutine smoothing_getpar
!=======================================================================
!=======================================================================
!=======================================================================


  subroutine parser(filename)
    implicit none
    character(len=*), intent(in) :: filename
    logical(lgt) :: exist
    character(len=filenamelen) :: line,name,value
    integer(i4b) :: i


    inquire(file=filename, exist=exist)
    if (.not. exist) then
       print 1, trim(filename)
1      format(/,&
            & " Error in parser:",/,&
            & " File ",a," does not exist.")
       stop 1
    end if

    ! --- set defaults:
    call setpar()

    open  (1,file=filename)
    do

       read (1,'(a)',end=11) line
       i=scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))

       select case (trim(name))
       case ('simul_type')
          read (value,*) par%simtp
       case ('nsmax')
          read (value,*) par%nsmax
       case ('nlmax')
          read (value,*) par%nlmax
       case ('iseed')
          read (value,*) par%iseed
       case ('fwhm_arcmin')
          read (value,*) par%fwhma
       case ('infile')
          par%inpfl = trim(value)
       case ('outfile')
          par%outfl = trim(value)
       case ('almsfile')
          par%almsfile = trim(value)
       case ('plmfile')
          par%plmfile = trim(value)
       case ('windowfile')
          par%windowfile = trim(value)
       case ('winfiledir')
          par%winfiledir = trim(value)
       end select
    end do
11  close (1)

  end subroutine parser


  subroutine setpar()
    implicit none

    par%simtp =   1
    par%nsmax =  32
    par%nlmax =  64
    par%iseed =  -1
    par%fwhma = 420.0
    par%inpfl = "cl.fits"
    par%outfl = "map.fits"
    par%almsfile = ""
    par%plmfile = ""
    par%windowfile = ""
    par%winfiledir = ""

  end subroutine setpar

  subroutine getpar(all)
    implicit none
    logical(lgt), optional :: all

    print *,'Non-interactive operation. The input file and defaults imply:'
    if (present(all)) then
       print 1, par%simtp,  &
            &   par%nsmax,  &
            &   par%nlmax,  &
            &   par%iseed,  &
            &   par%fwhma, &
            &   trim(par%inpfl),  &
            &   trim(par%outfl),  &
            &   trim(par%almsfile), &
            &   trim(par%plmfile),  &
            &   trim(par%windowfile),  &
            &   trim(par%winfiledir)
    else
       print 2, par%simtp, &
            &   par%nsmax, &
            &   par%nlmax, &
            &   par%iseed, &
            &   par%fwhma
    end if

1   format(/,&
         & " simulation type    = ",i5,/,&
         & " N_side             = ",i5,/,&
         & " l_max              = ",i5,/,&
         & " random number seed = ",i10,/,&
         & " FWHM (arc min.)    = ",f8.2,/,&
         & " input file         = ",a,/,&
         & " output file        = ",a,/,&
         & " alm file           = ",a,/,&
         & " precomputed P_lm   = ",a,/,&
         & " pixel window       = ",a,/,&
         & " pixel window dir.  = ",a,/)

2   format(/,&
         & " simulation type    = ",i5,/,&
         & " N_side             = ",i5,/,&
         & " l_max              = ",i5,/,&
         & " random number seed = ",i10,/,&
         & " FWHM (arc min.)    = ",f8.2,/)

  end subroutine getpar
!======================================================================================
  subroutine anafast_parser(filename)
    implicit none
    character(len=*) :: filename
    logical(lgt) :: exist
    character(len=filenamelen) :: line,name,value
    integer(i4b) :: i

    inquire(file=filename, exist=exist)
    if (.not. exist) then
       print 1, trim(filename)
1      format(/,&
            & " Error in anafast_parser:",/,&
            & " File ",a," does not exist.")
       stop 1
    end if

    ! --- set defaults:
    call anafast_setpar()

    open  (1,file=filename)
    do
       read (1,'(a)',end=11) line
       i=scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))

       select case (trim(name))
       case ('simul_type')
          read (value,*) anafast_par%simtp
       case ('nlmax')
          read (value,*) anafast_par%nlmax
       case ('theta_cut_deg')
          read (value,*) anafast_par%theta_cut_deg
       case ('regression')
          read (value,*) anafast_par%regression
       case ('iter_order')
          read (value,*) anafast_par%iter_order
       case ('w8file')
          anafast_par%w8file = trim(value)
       case ('infile')
          anafast_par%inpfl = trim(value)
       case ('outfile')
          anafast_par%outfl = trim(value)
       case ('plmfile')
          anafast_par%plmfile = trim(value)
       case ('outfile_alms')
          anafast_par%outfl_alms = trim(value)
       case ('w8filedir')
          anafast_par%w8filedir=trim(value)
       case ('won')
          read (value,*) anafast_par%won
       end select
    end do
11  close (1)

  end subroutine anafast_parser

  subroutine anafast_setpar()
    implicit none

    anafast_par%simtp =   1
    anafast_par%nlmax =  64
    anafast_par%theta_cut_deg = 0.0
    anafast_par%regression = 0
    anafast_par%iter_order = 0
    anafast_par%w8file = ''
    anafast_par%inpfl = "map.fits"
    anafast_par%outfl = "cl_out.fits"
    anafast_par%plmfile = ""
    anafast_par%outfl_alms = ""
    anafast_par%w8filedir = ''
    anafast_par%won = 2

  end subroutine anafast_setpar

  subroutine anafast_getpar(all)
    implicit none
    logical(lgt), optional :: all

    print *,'Non-interactive operation. The input file and defaults imply:'

    if (present(all)) then
       print 1, anafast_par%simtp,  &
            &   anafast_par%nlmax,  &
            &   anafast_par%regression,  &
            &   anafast_par%theta_cut_deg,  &
            &   anafast_par%iter_order, &
            &   trim(anafast_par%w8file),  &
            &   trim(anafast_par%inpfl),  &
            &   trim(anafast_par%outfl), &
            &   trim(anafast_par%plmfile),  &
            &   trim(anafast_par%outfl_alms),  &
            &   trim(anafast_par%w8filedir), &
            &   anafast_par%won
    else
       print 2, anafast_par%simtp, &
            &   anafast_par%nlmax, &
            &   anafast_par%theta_cut_deg, &
            &   anafast_par%iter_order
    end if

1   format(/,&
         & " simulation type     = ",i5,/,&
         & " l_max               = ",i5,/,&
         & " regres. low multipol= ",i2,/,&
         & " symmetric cut (deg) = ",f8.2,/,&
         & " iteration order     = ",i5,/,&
         & " weight file         = ",a,/,&
         & " input file          = ",a,/,&
         & " output file         = ",a,/,&
         & " precomputed Plm     = ",a,/,&
         & " alm output file     = ",a,/,&
         & " weight file dir.    = ",a,/,&
         & " weight parameter    = ",i1,/)

2   format(/,&
         & " simulation type     = ",i5,/,&
         & " l_max               = ",i5,/,&
         & " symmetric cut (deg) = ",f8.2,/,&
         & " iteration order     = ",i5,/)

  end subroutine anafast_getpar
!======================================================================================
  subroutine hotspots_parser(filename)
    implicit none
    character(len=*) :: filename
    logical(lgt) :: exist
    character(len=filenamelen) :: line,name,value
    integer(i4b) :: i

    inquire(file=filename, exist=exist)
    if (.not. exist) then
       print 1, trim(filename)
1      format(/,&
            & " Error in hotspots_parser:",/,&
            & " File ",a," does not exist.")
       stop 1
    end if

    ! --- set defaults:
    call hotspots_setpar()

    open  (1,file=filename)
    do
       read (1,'(a)',end=11) line
       i=scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))

       select case (trim(name))
       case ('infile')
          hotspots_par%inpfl = trim(value)
       case ('extrema_outfile')
          hotspots_par%outfl_extrema = trim(value)
       case ('maxima_outfile')
          hotspots_par%outfl_max = trim(value)
       case ('minima_outfile')
          hotspots_par%outfl_min = trim(value)
       end select
    end do
11  close (1)

  end subroutine hotspots_parser

  subroutine hotspots_setpar()
    implicit none

    hotspots_par%inpfl = "map.fits"
    hotspots_par%outfl_extrema = "pixlminmax.fits"
    hotspots_par%outfl_max = "maxima.dat"
    hotspots_par%outfl_min = "minima.dat"

  end subroutine hotspots_setpar

  subroutine hotspots_getpar(all)
    implicit none
    logical(lgt), optional :: all

    print *,'Non-interactive operation. The input file and defaults imply:'
    print 1, trim(hotspots_par%inpfl), &
         &   trim(hotspots_par%outfl_extrema), &
         &   trim(hotspots_par%outfl_max), &
         &   trim(hotspots_par%outfl_min)

1   format(/,&
         & " input file          = ",a,/,&
         & " extrema output file = ",a,/,&
         & " maxima output file  = ",a,/,&
         & " minima output file  = ",a,/)

  end subroutine hotspots_getpar
!======================================================================================
  subroutine udgrade_parser(filename)
    implicit none
    character(len=*) :: filename
    logical(lgt) :: exist
    character(len=filenamelen) :: line,name,value
    integer(i4b) :: i

    inquire(file=filename, exist=exist)
    if (.not. exist) then
       print 1, trim(filename)
1      format(/,&
            & " Error in udgrade_parser:",/,&
            & " File ",a," does not exist.")
       stop 1
    end if

    ! --- set defaults:
    call udgrade_setpar()

    open  (1,file=filename)
    do
       read (1,'(a)',end=11) line
       i=scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))

       select case (trim(name))
       case ('nside_out')
          read (value,*) udgrade_par%nside_out
       case ('infile')
          udgrade_par%inpfl = trim(value)
       case ('outfile')
          udgrade_par%outfl = trim(value)
       end select
    end do
11  close (1)

  end subroutine udgrade_parser

  subroutine udgrade_setpar()
    implicit none

    udgrade_par%nside_out = 64
    udgrade_par%inpfl = "map.fits"
    udgrade_par%outfl = "outmap.fits"

  end subroutine udgrade_setpar

  subroutine udgrade_getpar(all)
    implicit none
    logical(lgt), optional :: all

    print *,'Non-interactive operation. The input file and defaults imply:'

    print 1, udgrade_par%nside_out, trim(udgrade_par%inpfl), trim(udgrade_par%outfl)

1   format(/,&
         & " final nsmax         = ",i5,/,&
         & " input file          = ",a,/,&
         & " output file         = ",a,/)

  end subroutine udgrade_getpar

!======================================================================================
  subroutine smoothing_parser(filename)
    implicit none
    character(len=*) :: filename
    logical(lgt) :: exist
    character(len=filenamelen) :: line,name,value
    integer(i4b) :: i

    inquire(file=filename, exist=exist)
    if (.not. exist) then
       print 1, trim(filename)
1      format(/,&
            & " Error in smoothing_parser:",/,&
            & " File ",a," does not exist.")
       stop 1
    end if

    ! --- set defaults:
    call smoothing_setpar()

    open  (1,file=filename)
    do
       read (1,'(a)',end=11) line
       i=scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))

       select case (trim(name))
       case ('simul_type')
          read (value,*) smoothing_par%simtp
       case ('iter_order')
          read (value,*) smoothing_par%iter_order
       case ('nlmax')
          read (value,*) smoothing_par%nlmax
       case ('fwhm_arcmin')
          read (value,*) smoothing_par%fwhma
       case ('won')
          read (value,*) smoothing_par%won
       case ('plmfile')
          smoothing_par%plmfile = trim(value)
       case ('w8file')
          smoothing_par%w8file = trim(value)
       case ('w8filedir')
          smoothing_par%w8filedir=trim(value)
       case ('infile')
          smoothing_par%inpfl = trim(value)
       case ('outfile')
          smoothing_par%outfl = trim(value)
       end select
    end do
11  close (1)

  end subroutine smoothing_parser

  subroutine smoothing_setpar()
    implicit none

    smoothing_par%simtp = 1
    smoothing_par%nlmax = 64
    smoothing_par%fwhma = 420.0
    smoothing_par%iter_order = 0
    smoothing_par%inpfl = "map.fits"
    smoothing_par%outfl = "map_smoothed.fits"
    smoothing_par%w8file = ''
    smoothing_par%w8filedir = ''
    smoothing_par%plmfile = ""
    smoothing_par%won = 2

  end subroutine smoothing_setpar

  subroutine smoothing_getpar(all)
    implicit none
    logical(lgt), optional :: all

    print *,'Non-interactive operation. The input file and defaults imply:'

    print 1, &
         &   smoothing_par%simtp, &
         &   trim(smoothing_par%inpfl), &
         &   smoothing_par%fwhma, &
         &   smoothing_par%nlmax, &
         &   smoothing_par%iter_order, &
         &   trim(smoothing_par%outfl), &
         &   trim(smoothing_par%plmfile), &
         &   trim(smoothing_par%w8file), &
         &   trim(smoothing_par%w8filedir), &
         &   smoothing_par%won

1   format(/,&
         & " analysis type       = ",i5,/,&
         & " input file          = ",a,/,&
         & " FWHM (arc. min.)    = ",f8.2,/,&
         & " analysis l_max      = ",i5,/,&
         & " iteration order     = ",i5,/,&
         & " output file         = ",a,/,&
         & " precomputed Plm     = ",a,/,&
         & " weight file         = ",a,/,&
         & " weight file dir.    = ",a,/,&
         & " weight parameter    = ",i1,/)

  end subroutine smoothing_getpar

end module obsolete
