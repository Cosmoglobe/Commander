	integer*4 function FUT_VecPlt(xbuf,ybuf,npts,
     &		  plotmode,if2,if3,if4,if5,if6,if7,if8,if9,
     &		  plt_com, plt_com_file )

C-----------------------------------------------------------------------
C
C	Module Name: FUT_Vecplt
C
C	Purpose: To invoke the PLT routine for input vectors.
C
C	Author: Jody Caldwell, STX, September 1988.
C
C	Change Log:
C		 
C		Shirley M. Read, STX, September 1988.
C			Converted the subroutine to a function, reported
C			errors via Lib$Signal, and returned function status.
C			(SPRs 1566 and 1763)
C		        Determined min and max values for plot.
C		        Added command for automatic hardcopy.
C			Additional information for outer titles. (SPR 2430)
C
C               Shirley M. Read, STX, October, 1988 : Corrected the inclusion of
C		        the FUT_Labels text file in the subroutine. (SPR 2558)
C
C		Steven Alexander, STX, March 5, 1990 : Add capability to
C			pass in a PLT command file.
C
C----------------------------------------------------------------------

	implicit none

C	Include Files

	include '($ssdef)'
	include '(FUT_VecPlt_Labels)'
	include '(fut_params)'

C	External Messages

	external FUT_Normal
	external FUT_Aberr
	external FUT_TooManyPts

C	Functions

	integer*4 FUT_Erase	!Regis graphics screen erase
	
C	Local declarations

	integer*4 npts		!number of points in buffers
	real*4 xbuf(npts)	!buffer containing x points
	real*4 ybuf(npts)	!buffer containing y points
	integer*4 plotmode	! plotmode: screen, auto hardcopy
	integer*4 plt_com	!flag to indicate presence of PLT
				!command file
	character*64 plt_com_file	!PLT command file name
	integer*4 if2		!unknown flag?
	integer*4 if3		!unknown flag?
	integer*4 if4		!unknown flag?
	integer*4 if5		!unknown flag?
	integer*4 if6		!unknown flag?
	integer*4 if7		!unknown flag?
	integer*4 if8		!unknown flag?
	integer*4 if9		!unknown flag?
	real*4 xybuf(2048,2)	!input buffer to plt
	integer*4 mx1 /2048/	!Actual first dimension of y array
	integer*4 nveci /2/	!Actual number of vectors to be plotted
	integer*4 iery(5) /0,0,0,0,0/	! Flags for plot errors
	character*86 cmd(55)
	integer*4 ncmd		!number of commands to plt
	integer*4 ier		!error flag
	integer*4 status        !return status
	integer*4 i		!do loop index
	real*4 delta		!fraction beyond min and max range
	real*4 ymin, ymax       !min and max data values
	real*4 ypltmin, ypltmax !min and max for plot
	integer*4 screen   /1/  !write to CRT screen
	integer*4 autohard /2/  !automatic hardcopy
	logical*1 proceed       !flag to continue processing

C	Set function status and initialize.

	FUT_Vecplt = %loc(FUT_Normal)	
	proceed = .true.

	if(npts.gt.2048)then
	  call lib$signal(FUT_TooManyPts)
	  FUT_Vecplt = %loc(FUT_Aberr)
	  proceed = .false.
	endif

	if ( proceed ) then
	  cmd(1) = 'LA OT '//title
	  cmd(2) = 'LA X '//xlabl
	  cmd(3) = 'LA Y '//ylabl
	  cmd(4) = 'LA T '//title2
	  cmd(5) = 'LA OX '//oxlabl
	  cmd(6) = 'LA OY '//oylabl
	  cmd(7) = 'VP .1 .1 .9 .8'
	  do i=1,npts
	    xybuf(i,1)=xbuf(i)
	    xybuf(i,2)=ybuf(i)
	  enddo
	  ncmd = 8

C	Determine y min and max.

	  ymin = ybuf(1)
	  ymax = ybuf(1)
	  do i=2,npts
	    if ( ybuf(i) .lt. ymin ) ymin = ybuf(i)
	    if ( ybuf(i) .gt. ymax ) ymax = ybuf(i)
	  enddo
	  delta = 1.05 * (ymax - ymin)
	  if (delta .gt. 0) then
	    ypltmax = ymin + delta
	    ypltmin = ymax - delta
	  else
	    ypltmax = ymax + 0.1
	    ypltmin = ymin - 0.1
	  endif
	  write (cmd(ncmd),'(a, 2(1x,e14.4))') 'R Y ', ypltmin, ypltmax

C	If mode is auto hardcopy, set device for automatic writing to a file
C	of filename and device defined as logical Plt_Hardcopy.

	  if (plotmode .eq. autohard ) then
	    ncmd = ncmd + 1
	    cmd(ncmd) = 'DEVICE PLT_HARDCOPY'
	  endif

          if (plt_com .eq. fac_present) then
              ncmd = ncmd + 1
              cmd(ncmd) = '@' // plt_com_file
          endif

	  if (plotmode .eq. autohard ) then
	    ncmd = ncmd + 1
	    cmd(ncmd) = 'plot'
	    ncmd = ncmd + 1
	    cmd(ncmd) = 'quit'
	  endif

C	PLT interface.

	  call Plt(xybuf,iery,mx1,npts,nveci,cmd,ncmd,ier)

C	REGIS graphics erase screen.

	  status = FUT_Erase()
	  if (status .ne. %loc(FUT_Normal)) FUT_Vecplt = %loc(FUT_Aberr)
	endif

	return
	end
