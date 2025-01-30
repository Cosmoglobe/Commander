	INTEGER*4 FUNCTION FUT_XTALK_CHECK (DUMMY3,DUMMY4)
C/
C/	PROGRAM NAME:
C/	  FUT_XTALK_CHECK
C/
C/	PROGRAM DESCRIPTION:
C/	  The routine will be provided for COBE subsystem interaction data 
C/        quality check
C/
C/	AUTHOR:
C/	  H. WANG
C/	  STX
C/	  FEB 6, 1991
C/
C/	CALLING SEQUENCE:
C/	  STATUS =  FUT_XTALK_CHECK (DUMMY3,DUMMY4)
C/
C/	INPUT PARAMETERS:
C/
C/	OUTPUT PARAMETERS:
C/
C/	METHOD USED:
C/	  The following is the PDL --
C/          This is STUB routine
C/	  
C/	    Return;
C/	    End.
C/

	IMPLICIT	NONE

!	Passed Parameters
        BYTE DUMMY3(*), DUMMY4(*)


        EXTERNAL FDQ_NORMAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	Set return status to success.

	FUT_XTALK_CHECK = %loc(FDQ_NORMAL)


	return
	end
