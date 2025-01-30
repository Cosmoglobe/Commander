	INTEGER*4 FUNCTION FDQ_WITHIN (TIME1, TIME2, THRES, DIFF, 
	1	BIGGER )
C/
C/	PROGRAM NAME:
C/	  FDQ_WITHIN
C/
C/	PROGRAM DESCRIPTION:
C/	  This is a logical function that compares 2 input VAX binary times
C/	  and see if they are within a specified threshold.  The positive
C/	  difference is also returned through an output parameter.
C/
C/	AUTHOR:
C/	  EDWIN H. FUNG
C/	  GSFC
C/	  APRIL 24, 1987
C/
C/	MODIFIED BY:
C/	  N/A
C/
C/	CALLING SEQUENCE:
C/	  Return status = FDQ_WITHIN(TIME1, TIME2, THRES, DIFF, BIGGER)
C/
C/	INPUT PARAMETERS:
C/	  TIME1(2)	I*4	First input time in VAX binary format.
C/	  TIME2(2)	I*4	Second input time in VAX binary format.
C/	  THRES(2)	I*4	The threshold or comparison criterion between
C/				the 2 times (in VAX binary format).
C/
C/	OUTPUT PARAMETERS:
C/	  DIFF(2)	I*4	The positive difference between the 2 times
C/				(in VAX binary format).
C/	  BIGGER	I*2	Indicates which time is "bigger".  If TIME1
C/				is bigger (later), BIGGER = 1; else BIGGER = 2.
C/	  STATUS         I*4    Status returned if system function error occurs
C/
C/	INCLUDE FILES USED:
C/	  NONE
C/
C/	SUBROUTINES CALLED:
C/	  LIB$SUBX
C/
C/	METHOD USED:
C/	  Trivial
C/
C/      CHANGES:
C/
C/        J. Durachta
C/           ARC
C/         10/23/87
C/        Reason:        THR and DIF added in order to prevent alteration of
C/                       DIFF and THRES in the process of performing the check
C/                       as per the discovery of Edwin Fung.
C/
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 12, 1988
C/	  REASON: 	Added error checking for system quadword subtraction.
C/			Returned a system error status from this function 
C/			for interface with Fut_Error condition handler. Added 
C/			calls to Lib$Signal in case of error.
C/

	IMPLICIT	NONE

	INTEGER*2	BIGGER
	INTEGER*4	DIFF(2), THRES(2), TIME1(2), TIME2(2)
	integer*4	STATUS		! Return status variable

	include		'($SSDEF)'

	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_ERSUBX
	EXTERNAL	FDQ_WITHINTRUE
	EXTERNAL	FDQ_WITHINFALS

	integer*4 	LIB$SUBX

!!!!!!!!!!!!!!!!!!!!!!!!!
!			!
!     Local storage     !
!			!
!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	INTEGER*4	EARLIER(2), LATER(2)
        integer*4       dif, thr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	IF (TIME1(2) .GT. TIME2(2)) THEN
	  CALL LIB$MOVC3 (8, TIME1, LATER)
	  CALL LIB$MOVC3 (8, TIME2, EARLIER)
	  BIGGER = 1
	ELSE IF (TIME1(2) .LT. TIME2(2)) THEN
	  CALL LIB$MOVC3 (8, TIME1, EARLIER)
	  CALL LIB$MOVC3 (8, TIME2, LATER)
	  BIGGER = 2
	ELSE IF (TIME1(2) .EQ. TIME2(2)) THEN

!     Go to lower order I*4 for tie-breaker

	  IF (TIME1(1) .LT. 0) THEN
	    IF (TIME2(1) .LT. 0) THEN

!     Both have high order bit set; clear both so both numbers are then positive

	      TIME1(1) = IBCLR (TIME1(1), 15)
	      TIME2(1) = IBCLR (TIME2(1), 15)
	      IF (TIME1(1) .GT. TIME2(1)) THEN
	        CALL LIB$MOVC3 (8, TIME1, LATER)
	        CALL LIB$MOVC3 (8, TIME2, EARLIER)
	        BIGGER = 1
	      ELSE
	        CALL LIB$MOVC3 (8, TIME1, EARLIER)
	        CALL LIB$MOVC3 (8, TIME2, LATER)
	        BIGGER = 2
	      ENDIF
	    ELSE

!     TIME1(1) has high order bit set but TIME2 does not; TIME1 is later

	      CALL LIB$MOVC3 (8, TIME1, LATER)
	      CALL LIB$MOVC3 (8, TIME2, EARLIER)
	      BIGGER = 1
	    ENDIF
	  ELSE
	    IF (TIME2(1) .LT. 0) THEN

!     TIME1(1) has high order bit clear(positive) but TIME2(1) has high order bit set

	      CALL LIB$MOVC3 (8, TIME1, EARLIER)
	      CALL LIB$MOVC3 (8, TIME2, LATER)
	      BIGGER = 2
	    ELSE

!     Both has high order bit clear
	      IF (TIME1(1) .GT. TIME2(1)) THEN
	        CALL LIB$MOVC3 (8, TIME1, LATER)
	        CALL LIB$MOVC3 (8, TIME2, EARLIER)
	        BIGGER = 1
	      ELSE
	        CALL LIB$MOVC3 (8, TIME1, EARLIER)
	        CALL LIB$MOVC3 (8, TIME2, LATER)
	        BIGGER = 2
	      ENDIF
	    ENDIF
	  ENDIF
	ENDIF

	STATUS = LIB$SUBX (LATER, EARLIER, DIFF)

	IF ( STATUS .NE. SS$_NORMAL ) THEN
	  CALL LIB$SIGNAL(FDQ_ERSUBX, %VAL(1), %VAL(STATUS))
	  FDQ_WITHIN = %loc(FDQ_ABERR)
	ENDIF

	IF ( STATUS .EQ. SS$_NORMAL ) THEN
	  IF (DIFF(2) .EQ. THRES(2)) THEN

!     Tie, need to go to lower order longword for "tie-breaker"

	    if (DIFF(1) .LT. 0) then
	      if (THRES(1) .LT. 0) then

!     Both have high order bit set; clear both bits to make the lower order
!     longwords positive

	        dif = ibclr (diff(1), 15)
	        thr = ibclr (thres(1), 15)

	        if (dif .lt. thr) then
	          FDQ_WITHIN = %loc(FDQ_WITHINTRUE)
	        else
	          FDQ_WITHIN = %LOC(FDQ_WITHINFALS) 
	        endif
	      elseif (THRES(1) .GE. 0) then
	        FDQ_WITHIN = %loc(FDQ_WITHINFALS)
	      endif
	    elseif (DIFF(1) .GE. 0) then
	      if (THRES(1) .LT. 0) then
	        FDQ_WITHIN = %loc(FDQ_WITHINTRUE)
	      elseif (THRES(1) .GE. 0) then
!     Both lower order longwords positive
	        if (DIFF(1) .LT. THRES(1)) then
	          FDQ_WITHIN = %loc(FDQ_WITHINTRUE)
	        else
	          FDQ_WITHIN = %loc(FDQ_WITHINFALS)
	        endif
	      endif
	    endif
	  ELSEIF (DIFF(2) .GT. THRES(2)) THEN
	    FDQ_WITHIN = %loc(FDQ_WITHINFALS)
	  ELSEIF (DIFF(2) .LT. THRES(2)) THEN
	    FDQ_WITHIN = %loc(FDQ_WITHINTRUE)
	  ENDIF
	ENDIF		! Status = SS$_Normal
	RETURN
	END


