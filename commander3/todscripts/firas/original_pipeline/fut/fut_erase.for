C-------------------------------------------------------------------------------

	Integer*4 Function FUT_Erase

C-------------------------------------------------------------------------------
C
C	Purpose:
C
C		To erase the screen after plots are made on the terminal.
C
C	Author: Shirley M. Read
C		STX, September 1988
C
C	Invocation: status = FUT_Erase()
C
C	Modificaton History:
C
C	  Author	    Date	  Modification
C	  ----------------------------------------------------------------------
C
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	
C	Subroutines Called:
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C
C	Processing Method:
C
C		Write the Regis graphics erase sequence to the screen.
C
C------------------------------------------------------------------------------
C	
        implicit none

C	FUT messages

	external    FUT_normal
	external    FUT_aberr

	logical*1   normal 		! Process status flag
	character*6 inreg
	character*4 outreg
	byte inregis(6),outregis(4)
	equivalence (inregis(1),inreg), (outregis(1),outreg)
	data inregis/0,0,27,'P','1','p'/, outregis/0,0,27,'\'/

C	Initialize

	normal = .true.

C	REGIS graphics erase screen.

	write (6,50) inreg,outreg
50	format (a6,'S(E)W(I3)W(P1)',a4)

C       Exit the function.

	if (normal) then
	  FUT_erase = %loc(FUT_normal)
	else
	  FUT_erase = %loc(FUT_aberr)
	end if
	return
	end
