	Subroutine FUT_Format_Time ( binary_time, formatted_time )

C------------------------------------------------------------------------
C    PURPOSE: Formats binary time into an easily readable time string.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX
C            June 20, 1988
C
C    INVOCATION: CALL FUT_FORMAT_TIME ( BINARY_TIME, FORMATTED_TIME )
C
C    INPUT PARAMETERS:
C	BINARY_TIME(2)		I*4	Binary time to format.
C
C    OUTPUT PARAMETERS: 
C	FORMATTED_TIME		C*19	Formatted time.
C
C    SUBROUTINES CALLED: 
C	CT_BINARY_TO_GMT
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: 
C	CTUSER.INC
C
C----------------------------------------------------------------------

	Implicit None

	Include		'CT$Library:CTUser.Inc'

	Integer		*4	binary_time(2)
	Character	*19	formatted_time

	Character	*14	ascii_time

c
c Convert the input binary time to ASCII format.
c
	Call CT_Binary_To_GMT ( binary_time, ascii_time )

c
c Insert hyphens into the time string.
c
	formatted_time = ascii_time(1:2)   // '-' //
     .			 ascii_time(3:5)   // '-' //
     .			 ascii_time(6:7)   // '-' //
     .			 ascii_time(8:9)   // '-' //
     .			 ascii_time(10:11) // '-' //
     .			 ascii_time(12:14)

	Return
	End
