	INTEGER*4 FUNCTION FUT_TEMP_CORRECT (INPUT_TEMP,CORRECT_TEMP)
C/
C/	PROGRAM NAME:
C/	  FUT_TEMP_CORRECT
C/
C/	PROGRAM DESCRIPTION:
C/	  The routine will allow for correction of converted GRT temperatures
C/        due to errors in their calibration.
C/
C/	AUTHOR:
C/	  H. WANG
C/	  STX
C/	  FEB 6, 1991
C/
C/	CALLING SEQUENCE:
C/	  STATUS = FUT_TEMP_CORRECT (INPUT_TEMP, CORRECT_TEMP)
C/
C/	INPUT PARAMETERS:
C/	  INPUT_TEMP    record  GRT temperature
C/
C/	OUTPUT PARAMETERS:
C/	  CORRECT_TEMP  record	Corrected temperatures
C/	METHOD USED:
C/	  The following is the PDL --
C/          This is STUB routine
C/	  
C/	    Return;
C/	    End.
C/

	IMPLICIT	NONE

!	Passed Parameters

	DICTIONARY 'FUT_enganlg'
	RECORD /FUT_ENGANLG/ INPUT_TEMP,CORRECT_TEMP


        EXTERNAL FDQ_NORMAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	Set return status to success.

	FUT_TEMP_CORRECT = %loc(FDQ_NORMAL)


	return
	end
