	subroutine g_to_j (gstring)
c***********************************************************************
c
c  This routine takes takes the nine-character Gregorian date and converts
c  it to Julian format, for output to the terminal.
c
c  Written 2 April 1987 by Rich Isaacman, Applied Research Corp.
c
c  Input:		character*9 gstring (Gregorian date DD-MMM-YY)
c
c  Output:		Julian date YYDDD at terminal
c
c  Subroutines called:	None
c 
c  Include files:	None
c
c**************************************************************************
	implicit none
	character*3 	months(12)
	character*9 	gstring
	integer*4 	days_in_mo(12)
	integer*4 	j
	integer*4 	iasc_val
	integer*4 	iyear
	integer*4 	leapyr
	integer*4 	iday
	integer*4 	daytot
	integer*4 	jmo
	logical*1	foundit

	data months	/'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
     .			 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/

	data days_in_mo /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
c
c  Convert year and day numbers to integer, then check for leap year
c
	iyear = 10*(ichar(gstring(8:8))-48) + ichar(gstring(9:9))-48
	iday = 10*(ichar(gstring(1:1))-48) + ichar(gstring(2:2))-48
	leapyr = mod (iyear, 4)
	if (leapyr .eq. 0) days_in_mo(2) = 29
c
c  Convert month part of string into uppercase letters by adjusting 
c  ascii values as required.
c
	do j=4,6
	   iasc_val = ichar(gstring(j:j))
	   if (iasc_val .ge. 97) gstring(j:j) = char (iasc_val-32)
	enddo
c
c  Step thru list of months to find a match for the input string. Each 
c  time a month doesn't match, accumulate its number of days in the total.
c  "Foundit" flag is updated every time a comparison is made; loop 
c  terminates when month is found or all 12 months have been checked.
c
	daytot = 0
	jmo = 0
	do while (jmo.lt.12 .and. .not.foundit)
	   jmo = jmo + 1
	   if (gstring(4:6) .ne. months(jmo)) then
	      daytot = daytot + days_in_mo(jmo)
	      foundit = .false.
	   else
	      foundit = .true.
	   endif
	enddo
c
c  Bail out if no month match was found.
c
	if (.not.foundit) then
	   type 200, gstring(4:6)
200	   format (' My calendar seems to be missing a page. I can''t ',
     .		   'find a month called "',a3,'".'/)
	   return
	endif
c
c  Month has been found. Bail out if day number is too large for that
c  month (e.g. February 31....)
c
	if (iday .gt. days_in_mo(jmo)) then
	   type 300
300	   format (' There seem to have been more than the usual number',
     .		   ' of days in that month....'/)
	   return
	endif
c
c  Output integer year and total day number with appropriate format so 
c  there are no blanks in output record.
c
	iday = iday + daytot
	if (iday .ge. 100) then
	   type 400, iyear, iday
400	   format (' Julian day: ',i2,i3/)
	elseif (iday .ge. 10) then
	   type 500, iyear, iday
500	   format (' Julian day: ',i2,1h0,i2/)
	else
	   type 600, iyear, iday
600	   format (' Julian day: ',i2,2h00,i1/)
	endif
	return
	end
