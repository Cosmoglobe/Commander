	INTEGER*4 FUNCTION FDQ_STAT_INIT (ENG_REC, period, start_time, 
     1				stop_time)

C/	PURPOSE:  Calculate Avg, Sigma, Min, Max, for TRENDPLOTS for 
C/		  each orbit and write into FDQ_ETR (Eng Stat Rec).
C/
C/	PROGRAM NAME:
C/	  FDQ_STAT_INIT
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine initalizes the engineering statistics arrays.
C/
C/	AUTHOR:
C/	  D.T.WARD
C/	  GSFC/STX
C/	  OCT. 28, 1987
C/
C/	MODIFIED BY:
C/
C/	  D.WARD
C/	  GSFC/STX
C/	  NOV. 13, 1987
C/	  REASON:	To change the summation of element to a summation of
C/			differences from element1.
C/
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 15, 1988
C/	  REASON: 	Converted for interface with the Fut_Error 
C/			condition handler. Added error checking
C/		        and calls to Lib$Signal.
CH	CHANGE LOG:
CH
CH	Version 4.2.1 02/09/89, SPR 3139, Shirley M. Read, STX
CH		FDQ does not weed out the bad data flag value of -9999.0 when
CH		it forms the engineering trends records. Modules affected may be
CH		FDQ_Init_Stats, FDQ_Gen_Stats and FDQ_Form_Stats. In computing 
CH		the statistics, the flag value should not be included in the
CH		number, sum or sum of squares.
C/
C/	CALLING SEQUENCE FROM FDQ_GEN_STAT:
C/	  return status = FDQ_STAT_INIT (ENG_REC, period, start_time, stop_time)
C/
C/	INPUT PARAMETERS:
C/	  ENG_REC(1024)		BYTE	1024-byte record buffer to be
C/					output to FIRAS Eng. Archive
C/	  PERIOD(2)		I*4	I*8 form of the period
C/	  START_TIME(2)		I*4	I*8 form of the start time
C/	  STOP_TIME(2)		I*4	I*8 form of the stop time
C/
C/	OUTPUT PARAMETERS:
C/	  NONE
C/
C/	INPUT FILES:
C/	  NONE
C/
C/	OUTPUT FILES:
C/	  NONE
C/
C/	INCLUDE FILES:
C/	  CT$LIBRARY:CTUSER.INC
C/	  FDQ_STATBUFF.TXT
C/	  ($SSDEF)
C/
C/	SUBROUTINES CALLED:
C/	  CT_BINARY_TO_GMT
C/	  
C/	FUNCTIONS CALLED:
C/	  NONE
C/
C/	ERROR HANDLING:
C/	  Use DEBUG and type statements
C/
C/	METHOD USED:
C/
C/	PDL by Doug Ward  October 15, 1987
C/
C/	  Set the trends record count to 1.
C/	  Set start time of the trends record to FDQ_ENG start time 
C/	      (of the first segment/ first orbit?)
C/	  CALL LIB$ADDX to set stop time to start time plus orbital period
C/	  Initialize minimum, maximum, sum, sum_sq and first field value for
C/	      each engineering field.
C/	  RETURN  


	IMPLICIT	NONE

	DICTIONARY 'FDQ_ENG'
	RECORD /FDQ_ENG/ ENG_REC

	INCLUDE		'CT$LIBRARY:CTUSER.INC/NOLIST'

	INCLUDE		'(FDQ_STATBUFF)'

	INCLUDE		'($SSDEF)'

	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_ERADDX

	INTEGER*4	start_time(2), stop_time(2), 
	1		period(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			  !
!     Local Variables     !
!			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

	LOGICAL*1	FIRST_RECORD /.TRUE./

	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status

	CHARACTER*14	GMT_stop_time

	INTEGER*4	I, ISTAT, LIB$ADDX

	REAL*4		 n_orbit
                 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SAVE FIRST_RECORD	!=TRUE 1st time through only

C	Set return status to success.

	RETSTAT = SUCCESS

!!!!!!!!!!!!!!!!!!!!!!!!!
!			!
!	INITIALIZE	!
!			!
!!!!!!!!!!!!!!!!!!!!!!!!!

	  n_rec = 1

	  rec_time (1,n_rec, 1) = eng_rec.ct_head.TIME(1)
	  rec_time (2,n_rec, 1) = eng_rec.ct_head.TIME(2)

!
!	RESET start_time and stop_time
!

	  call LIB$MOVC3 (8, eng_rec.ct_head.TIME, start_time)

	  istat = LIB$ADDX ( start_time, period, stop_time)

	  if ( ISTAT .ne. SS$_NORMAL ) then
	    call lib$signal(FDQ_ERADDX,%val(1),%val(istat))
	    RETSTAT = ERROR
!	    write (6, 6150)  FIRST_RECORD, istat
 6150	FORMAT (' **** IN FDQ_GEN_STAT, FIRST_RECORD =',8A,/
     1	' UNSUCCESSFUL RETURN FROM LIB$ADDX ****.  istat=',I4)
	  endif

	  if ( retstat .eq. success ) then
	    call CT_BINARY_TO_GMT (stop_time, GMT_stop_time)
!
! Initialize min, max, sum, sum_sq, first engineering field value and rec count
!
	    DO	I = 1, 128
	      max(i) = 0.0
	      min(i) = 0.0
	      sum(i) = 0.0
	      sum_sq(i) = 0.0
	      element1(i) = 0.0
	      nf_rec(i) = 0
	    ENDDO

	    FIRST_RECORD = .FALSE.
	  endif

!	Set function to return status

	  IF (RETSTAT.EQ.SUCCESS) THEN
	    FDQ_STAT_INIT = %loc(FDQ_NORMAL)
	  ELSE
	    FDQ_STAT_INIT = %loc(FDQ_ABERR)
	  ENDIF
	  
	  RETURN
  	  END
