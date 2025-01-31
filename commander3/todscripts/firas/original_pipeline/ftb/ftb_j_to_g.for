	subroutine j_to_g (jstring)
c***********************************************************************
c
c  This routine takes takes the five-character Julian date and converts
c  it to Gregorian format, for output to the terminal.
c
c  Written 2 April 1987 by Rich Isaacman, Applied Research Corp.
c
c  Input:		character*5 jstring (Julian date YYDDD)
c
c  Output:		Gregorian date DD-MMM-YY at terminal
c
c  Subroutines called:	None
c 
c  Include files:	None
c
c**************************************************************************
	implicit none
	character*3 	months(12)
	character*5	jstring
	integer*4	days_in_mo(12)
	integer*4 	iyear
	integer*4 	leapyr
	integer*4 	maxdays
	integer*4 	whichday
	integer*4 	daytot
	integer*4 	iday
	integer*4 	jmo

	data months	/'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
     .			 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/

	data days_in_mo /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
c
c  First convert character year to integer, and check for leap year
c
	iyear = 10*(ichar(jstring(1:1))-48) + ichar(jstring(2:2))-48
	leapyr = mod (iyear, 4)
	if (leapyr .eq. 0) then
	   days_in_mo(2) = 29
	   maxdays = 366
	else
	   maxdays = 365
	endif
c
c  Convert character day number to an integer and make sure it's less
c  than 365 or 366.
c
	iday = 100*(ichar(jstring(3:3))-48) + 10*(ichar(jstring(4:4))-48) + 
     .		    ichar(jstring(5:5))-48

	if (iday .gt. maxdays) then
	   type *, 'It''s been a long year, but not THAT long...!'
	   type *, ' '
	   return
	endif
c
c  Add up number of days in successive months until total exceeds the 
c  day number. At that point we know which month we're in.
c
	jmo = 0
	daytot = 0
	do while (jmo.le.12 .and. daytot.lt.iday)
	   jmo = jmo + 1
	   daytot = daytot + days_in_mo(jmo)
	enddo
c
c  Calculate how many days into the month we are, then show result
c
	daytot = daytot - days_in_mo(jmo)
	whichday = iday - daytot
	type 300, whichday, months(jmo), iyear
300	format (' Gregorian date:  ',i2,'-',a3,'-',i2/)
	return
	end
