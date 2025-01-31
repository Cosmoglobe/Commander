	INTEGER*4   FUNCTION FUT_WRITE_ERROR (LINE)
c
c  Purpose:
c  This routine is called by the FUT_ERROR subroutine
c  to write error message(s) into a facility report file 
c  using the SYS$PUTMSG system command.
c  
c  Written by: Nilo G. Gonzales/STX, 12-19-90
c
c  Input:
C    LINE    C*        error message
C  Output:
C    NONE 
C  PDL:
C     IF fortran unit number is greater than 0
C     Then
C        write error messages into a facility report file
C     endif
C     return
C     end
CH---------------------------------------------------------------
C  Changes:
C         SPR 2805, Add second line if error message is greater 
C                 than 80 characters in length.
C                 Nilo G. Gonzales/STX, Feb. 27, 1991.
CH---------------------------------------------------------------
C
C     Code begin
C    
	IMPLICIT NONE

	INCLUDE    '(FUT_ERROR)'
	CHARACTER*(*) LINE
	INTEGER*4  LENL
	
        LENL = LEN(LINE)
        IF (FUT_REPORT_LUN .GT. 0) THEN
	    IF (LENL .GT. 79) THEN
	       WRITE (FUT_REPORT_LUN, 15) LINE(1:79), LINE(80:LENL)
15	       FORMAT(/1X,A/1X,A)
	    ELSE
               WRITE (FUT_REPORT_LUN,20) LINE
20	       FORMAT(/1X,A)
	    ENDIF
	ENDIF

        RETURN
        END
