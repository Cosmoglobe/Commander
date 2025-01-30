
	PROGRAM FTB_JDATE

c**************************************************************************
c
c  This program accepts either a Julian date (e.g. 87201) or a Gregorian
c  date (e.g. 23-APR-87) from the terminal and returns the same date
c  in the opposite form. It handles leap years and screwy inputs (like
c  nonexistent months or months with too many days).
c
c  Written 2 April 1987 by Rich Isaacman, Applied Research Corp.
c
c  Input:		from terminal only
c
c  Include files:	
c	   $ssdef
c
c  Subroutines called:	g_to_j, j_to_g
c	
c  Modified by:
c  Shirley M. Read
c  August 26, 1988
c  STX Inc.
c  Reason: SPRs 1566 and 1763 requested that all FIRAS facilities 
c          should tell of successful completion and signal errors
c	   via $status for batch processing. Also added interface 
c	   with Fut_Error, the condition handler.
c
c************************************************************************
	implicit none

c	FTB messages

	external   	ftb_normal
	external	ftb_aberr

c	Condition Handler

	external        fut_error

c	Include Files

	include         '($ssdef)'

c       Local variables

	integer*4 	status		      ! Status of processing
	integer*4 	success/1/, error/2/  ! Local status values
	character*40	input
	character*8	jstring
	character*9	gstring
	integer*4	idash
	integer*4	i
	integer*4	ipos
	integer*4	ndigit
	integer*4	iasc_val


c	Set status for FTB processing to success.

	status = success

c     Establish condition handler.            

	call lib$establish ( fut_error )

	type 80
80	format (/' Enter Julian (e.g. 87201) or Gregorian (e.g. ',
     .		 '23-APR-87) date: ',$)
	accept 180, input
180	format (a40)

c
c  Look for a dash to determine whether input is Gregorian or Julian date
c
	idash = index (input, '-')
	if (idash .ne. 0) then
c
c      Date is Gregorian. Fill gstring array with left-justified character 
c      string. Stick in extra zero at left side if day number < 10.
c
	   if (idash .eq. 2) then
	      do i=9,2,-1
	         input(i:i) = input(i-1:i-1)
	      enddo
	      input(1:1) = '0'
	      idash = 3
	   endif

	   do i=1,9
	      ipos = idash + i - 3
	      gstring(i:i) = input(ipos:ipos)
	   enddo
c
c      Now send the 9-character string for conversion
c      
	   call g_to_j (gstring)

	else
c
c      Date is Julian. Pick out the first five valid digits (ascii value
c      48-57 for characters 0-9) to fill left-justified jstring.
c
	   ndigit = 0
	   do i=1,40
	      iasc_val = ichar (input(i:i))
	      if (iasc_val.ge.48 .and. iasc_val.le.57
     .				 .and. ndigit.le.4) then
	         ndigit = ndigit + 1
		 jstring(ndigit:ndigit) = input(i:i)
	      endif
	   enddo
c
c      Convert the string after finding five digits. Otherwise give up.
c
	   if (ndigit .eq. 5) then
	      call j_to_g (jstring)
	   else
	      type *, 'That has too few digits for a Julian date!'
	      type *, ' '
	   endif

	endif

c       Exit the program.

	if ( status .eq. success ) then
	  call lib$signal(ftb_normal)
	  call exit(ss$_normal)
	else
	  call lib$signal(ftb_aberr)
	  call exit(ss$_abort)
	endif
	end
             
