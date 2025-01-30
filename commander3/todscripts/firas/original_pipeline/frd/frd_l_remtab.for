
	Integer*4 Function FRD_L_Remtab ( Buffer )

!	Program Name : FRD_L_Remtab
!
!	Programmer: Shirley M. Read, STX, January 1988
!
!	Program Description:
!	    This subroutine removes tabs or other unprintable characters 
!	    from the input buffer string and replaces them with blanks.
!
!----------------------------------------------------------------------------

	Implicit None

!	Passed Parameters.

	Character*80	Buffer	  ! String of 80 characters.

!	Include files.

	Include '($SSDef)'

!	Functions.

	Integer*4 Str$Translate    ! String translate characters.

!	Local variables.

	Integer*4 retstat	   ! Program Processing status
	Integer*4 success / 1 /, error / 2 /
	Integer*4 status	   ! System returned status
	Integer*4 ix               ! Index
	Character*5 unprintable(2)
	Data unprintable(1)(1:1) / '08'x /	! Backspace
	Data unprintable(1)(2:2) / '09'x /	! Horizontal tab
	Data unprintable(1)(3:3) / '0A'x /	! Linefeed
	Data unprintable(1)(4:4) / '0B'x /	! Vertical tab
	Data unprintable(1)(5:5) / '0C'x /	! Formfeed
	Data unprintable(2) / '     ' /         ! Translation table

!       Check the characters in the input buffer and replace unprintable
!	characters with blanks.

D	write(6,100) buffer(1:40), buffer(41:80)

	status = str$translate ( buffer, buffer, unprintable(2),
	1	 unprintable(1) )

	if ( status .eq. SS$_Normal ) then
	  retstat = success
D	  write (6,200) buffer(1:40), buffer(41:80)
	else
	  retstat = error
	  write(6,300) status
	endif
	
!	Set function return status.

	FRD_L_Remtab = retstat

	return

!	Debug Formats:
 100    format(1x,'String before:'/1x,a/1x,a)
 200    format(1x,'String after: '/1x,a/1x,a)
!	Processing Formats:
 300    format(1x,'Error: Bad return from Str$Translate. Status=',z8.8)

	end

