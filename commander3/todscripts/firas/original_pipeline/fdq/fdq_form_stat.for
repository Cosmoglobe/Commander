	INTEGER*4 FUNCTION FDQ_FORM_STAT (ETR_LUN,ETR_WRIT)
C-----------------------------------------------------------------------------
C/
C/	PROGRAM NAME:
C/	  FDQ_FORM_STAT
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine will form the engineering statistics as FDQ
C/	  has finished processing the required interval of Science records.
C/	  It is used at the end of every orbit and at the end of the
C/	  segment set read by FDQ.
C/
C/	AUTHOR:
C/	  E. FUNG
C/	  GSFC
C/	  DEC. 16, 1985
C/
C/      Modified by: 
C/              H. Wang, STX, 1/29/91
C/              New requirements for FDQ
C/              Output the number of FDQ_ETR records that have been written 
C/              to the FDQ_ETR file.
C/  
C/	CALLING SEQUENCE FROM FDQ_GEN_STAT AND FDQ:
C/	  FDQ_FORM_STAT (ETR_LUN)
C/
C/	INPUT PARAMETERS
C/	  ETR_LUN		I*2	logical unit no. of Eng. Trends
C/
C/	OUTPUT PARAMETERS:
C/	  ETR_WRIT              I*4     Number of NEW ETR records 
C/
C/	INPUT FILES:
C/	  NONE
C/
C/	OUTPUT FILES:
C/	  FIRAS ENG. TRENDS ARCHIVE (ETR_REC(4))
C/
C/	INCLUDE FILES:
C/	  CT$LIBRARY:CTUSER.INC
C/	  FDQ_STATBUFF.TXT
C/
C/	SUBROUTINES CALLED:
C/	  LIB$MOVC5
C/	  LIB$MOVC3
C/	  CT_BINARY_TO_GMT
C/
C/	FUNCTIONS CALLED:
C/	  FUT_AVE_BIN_TIMES
C/
C/	ERROR HANDLING:
C/	  Use DEBUG compiler and examine variables computed in code below.
C/
C/	METHOD USED:
C/	  The following is the PDL --
C/
C/	  IF (record count > 0) THEN
C/	    DOFOR each parameter (128)
C/	      Calculate Avg
C/	      Calculate Sigma
C/	      WRITE avg, sigma, min, max into data structure
C/            Increment the number of FDQ_ETR records  
C/	    ENDDOFOR
C/	    CALL FUT_AVE_BIN_TIMES to calculate avg time from time array
C/	    WRITE time into data structure
C/	  ENDIF
C/	  RETURN
C/	  END
C-----------------------------------------------------------------------------
CH	CHANGE LOG:
CH
CH	  E.FUNG
CH	  GSFC
CH	  DEC. 28, 1985
CH	  REASON:	(1) Call SYS$ASCTIM and SYS$BINTIM to get system time at
CH			the time of calling CT_WRITE_ARCV to the Eng. Statistics
CH			File.
CH			(2) Correct bug in reading nonexistent mission statistics
CH			record.
CH			(3) Correct bug in calculating GRANDSUMSQ(I).
CH
CH	  E.FUNG
CH	  GSFC
CH	  JAN. 2, 1985
CH	  REASON:	To avoid "taking the square root of a neg. number" due
CH			to truncation errors, add check before taking the square
CH			root.  If the number inside the square root sign is
CH			sufficiently small, the RMS will be set to zero.
CH
CH	  E.FUNG
CH	  GSFC
CH	  JAN. 6, 1986
CH	  REASON:	To go from 1 eng. statistics file to 2 (hourly and daily
CH			being separate files); this necessitates one additional
CH			input parameter (HES_LUN and DES_LUN replacing EST_LUN).
CH
CH	  D.WARD
CH	  GSFC/STX
CH	  OCT. 27, 1987
CH	  REASON:	To replace HES and DES with ETR, orbital only.
CH
CH	  D.WARD
CH	  GSFC/STX
CH	  OCT. 28, 1987
CH	  REASON:	To upgrade the function to use ETR_REC instead of
CH			FDQ_STATBUFF.TXT, to form orbital statistics only,
CH			and to be used for every orbit as well as at the
CH			end of the segment set read by FDQ.
CH
CH	  Shirley M. Read
CH	  STX
CH	  January 25, 1988
CH	  REASON:	Converted from subroutine to function for interface
CH			with Fut_Error condition handler. Added error checking
CH		        and calls to Lib$Signal.
CH
CH
CH	Version 4.2.1 02/09/89, SPR 3139, Shirley M. Read, STX
CH		FDQ does not weed out the bad data flag value of -9999.0 when
CH		it forms the engineering trends records. Modules affected may be
CH		FDQ_Init_Stats, FDQ_Gen_Stats and FDQ_Form_Stats. In computing
CH		the statistics, the flag value should not be included in the
CH		number, sum or sum of squares.
CH
CH	Version 4.5,  1989 Oct 10,  SPR 4670,  Fred Shuman,  STX.
CH	    FIT, which reads and plots FDQ_ETR records, one for each of five
CH	    (previously four) types of statistics, was making garbled plots.
CH	    This was caused by code in FIT which did not recognize the fifth
CH	    ETR record.  The fix for this included making the number of
CH	    statistics types a parameter, and this has been carried over
CH	    to FDQ, as well.  FDQ_FORM_STAT and FDQ_GEN_STAT.
C-----------------------------------------------------------------------------

	IMPLICIT	NONE

	DICTIONARY 'FDQ_ENG'
	RECORD /FDQ_ENG/ ENG_REC

	INTEGER*2	ETR_LUN

	INCLUDE		'CT$LIBRARY:CTUSER.INC/NOLIST'

	INCLUDE		'(FDQ_STATBUFF)'

	INCLUDE		'($SSDEF)'

	EXTERNAL	FUT_NORMAL
	EXTERNAL	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_CTWRITERR
	EXTERNAL        FDQ_STATLOOPER
	EXTERNAL        FDQ_STATMODER
	EXTERNAL        FDQ_AVBINTIMER
	EXTERNAL        FDQ_NEGRADER


!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			  !
!     Local variables     !
!			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

	INTEGER*4	num_stats_types
	PARAMETER	( num_stats_types = 5 )
	DICTIONARY 'FDQ_ETR'
	RECORD /FDQ_ETR/ ETR_REC(num_stats_types)
	RECORD /FDQ_ETR/ DUM_REC

	CHARACTER*14	GMT_stop_time, GMT_avg

	INTEGER*2       K, CT_STAT(20)

	INTEGER*4	n_loops, n_modulo,ETR_WRIT

	INTEGER*4	avg_bin_time(2), aweights(3)

	INTEGER*4	interm_avg_bin_time(2,3)

	INTEGER*4	FUT_AVE_BIN_TIMES

	INTEGER*4	I, ISTAT, BAD, ZERO / 0 /

	REAL*4		avg(128), sigma(128), var, rec,
	1		avg_diff, sumXavg_diff
	REAL*4	        FV / -9999.0 /  ! Engineering bad value flag

	PARAMETER	( BAD	= 2 )
	integer*4	RETSTAT		! Return status
	integer*4	SUCCESS / 1 /, ERROR / 2 /  ! Values for status
	integer*4	STATUS		! Dummy status variable


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set return status to success.

	RETSTAT = SUCCESS

	CALL LIB$MOVC5 (0, ,0, 4*1024, ETR_REC)		! Initialize ETR_REC

	IF (n_rec .GT. 0) THEN

!
!	Compute the average time for the orbital record set
!	 (FUT_AVE_BIN_TIMES takes only 100 times at once.)
!

	  n_loops = (n_rec - 1)/100

	  IF (n_loops .LT. 0) then
	     Call Lib$Signal (FDQ_STATLOOPER, %val(1), %val(n_loops))
	  ENDIF

	  IF (n_loops .GT. 0) THEN
	    do k = 1, n_loops
	     If ( retstat .eq. success ) then
	      istat = FUT_AVE_BIN_TIMES (rec_time (1, 1, k),
     1		        100, weights (1, k), interm_avg_bin_time (1, k) )
	      IF (istat .ne. %loc(FUT_NORMAL)) THEN
	        call Lib$Signal ( FDQ_AVBINTIMER, %val(1), %val(istat))
	        RETSTAT = ERROR
	      ENDIF
	      aweights (k) = 100
	     Endif	! Retstat is success
	    ENDDO
	  ENDIF

	  if ( retstat .eq.success ) then
	    n_modulo = n_rec - n_loops*100

	    If (n_modulo .LT. 0 ) Then
	      Call Lib$Signal (FDQ_STATMODER, %val(1), %val(n_modulo))
	    End If
	  endif  ! Retstat is success

	  if (Retstat .eq. success ) then
	    istat = FUT_AVE_BIN_TIMES (rec_time (1, 1, n_loops + 1),
     1	    n_modulo, weights (1, n_loops+1),
     1	    interm_avg_bin_time (1, n_loops+1) )
	    if (istat .ne. %loc(fut_normal)) then
	       call Lib$Signal (FDQ_AVBINTIMER, %val(1), %val(istat))
	       RETSTAT = ERROR
	    endif
	    aweights (n_loops + 1) = n_modulo
	  endif		! Retstat is success

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								     !
!	AVERAGE THE (N_LOOPS + 1) INTERMEDIATE_AVERAGES WEIGHTED     !
!	ACCORDING TO THE NUMBER OF AVERAGES USED IN MAKING EACH.     !
!								     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  if ( retstat .eq. success ) then

	    istat = FUT_AVE_BIN_TIMES (interm_avg_bin_time (1, 1),
     1	    n_loops + 1, aweights (1), avg_bin_time )
	    IF (istat .ne. %loc(fut_normal)) then
	      RETSTAT = ERROR
	      Call Lib$Signal (FDQ_AVBINTIMER, %val(1), %val(istat))
	    ENDIF
	  endif	! Retstat is success

	  if ( Retstat .eq. success ) then

	    call CT_BINARY_TO_GMT (avg_bin_time, GMT_avg)

	    DO	I = 1, 128
	     if ((Retstat .eq. success) .and. (nf_rec(i) .ne. zero)) then
	      rec = FLOAT (nf_rec(i))
	      avg_diff = sum(i)/rec
	      avg(i) = avg_diff + element1(i)
	      sumXavg_diff = sum(i)*avg_diff

	      IF ( nf_rec(i) .EQ. 1 ) THEN
	        sigma(i) = 0.0		! NO unbiased standard deviation
	      ELSE
	        IF (sum_sq(i) .LT. sumXavg_diff) then
	          Call Lib$Signal(FDQ_NEGRADER, %val(1), %val(i))
		  RETSTAT = ERROR
	        ELSE
	          sigma(i) = SQRT((sum_sq(i) - sumXavg_diff)/(rec - 1))
	        ENDIF
	      ENDIF
	     else	! All eng values were flagged
	      avg(i) = FV
	      sigma(i) = 0.0
	      min(i) = FV
	      max(i) = FV
	     endif	! Retstat is success and nf_rec(i) ne zero
	    ENDDO

	  endif		! Retstat is success

!	Another trends record containing the good value counts for each
!	engineering field was added to the FDQ_ETR file.

	  DO K = 1,num_stats_types
	   if ( retstat .eq. success ) then
	    call LIB$MOVC3 (8, avg_bin_time, ETR_REC(K).CT_HEAD.TIME)
	    ETR_REC(K).CT_HEAD.GMT = GMT_avg
	    ETR_REC(K).EN_HEAD.KEY_ID = K           ! I*2
	    ETR_REC(K).EN_HEAD.NUMBER_OF_RECORDS = n_rec
	    call LIB$MOVC3 (8, rec_time(1, 1, 1),
     1		 ETR_REC(K).EN_HEAD.FIRST_ENG_TIME)
	    call LIB$MOVC3 (8, rec_time(1, n_modulo, n_loops + 1),
     1		 ETR_REC(K).EN_HEAD.LAST_ENG_TIME)
             
            Do i=1, 64
              IF (K .EQ. 1) ETR_REC(K).EN_ANALOG.GRT(i) = avg(i)
              IF (K .EQ. 2) ETR_REC(K).EN_ANALOG.GRT(i) = sigma(i)
	      IF (K .EQ. 3) ETR_REC(K).EN_ANALOG.GRT(i) = min(i)
	      IF (K .EQ. 4) ETR_REC(K).EN_ANALOG.GRT(i) = max(i)
	      IF (K .EQ. 5) ETR_REC(K).EN_ANALOG.GRT(i) = nf_rec(i)
	    ENDDO 
	    DO i = 1,62
		  IF (K .EQ. 1) ETR_REC(K).EN_ANALOG.GROUP1(i) = avg(i+64)
		  IF (K .EQ. 2) ETR_REC(K).EN_ANALOG.GROUP1(i) = sigma(i+64)
		  IF (K .EQ. 3) ETR_REC(K).EN_ANALOG.GROUP1(i) = min(i+64)
		  IF (K .EQ. 4) ETR_REC(K).EN_ANALOG.GROUP1(i) = max(i+64)
		  IF (K .EQ. 5) ETR_REC(K).EN_ANALOG.GROUP1(i) = nf_rec(i+64)
            Enddo
	    IF (K .EQ. 1) ETR_REC(K).EN_TAIL.LMAC_ANALOG_TEMP = avg(127)
            IF (K .EQ. 2) ETR_REC(K).EN_TAIL.LMAC_ANALOG_TEMP = sigma(127)
            IF (K .EQ. 3) ETR_REC(K).EN_TAIL.LMAC_ANALOG_TEMP = min(127)
            IF (K .EQ. 4) ETR_REC(K).EN_TAIL.LMAC_ANALOG_TEMP = max(127)
            IF (K .EQ. 5) ETR_REC(K).EN_TAIL.LMAC_ANALOG_TEMP = nf_rec(127)
	    IF (K .EQ. 1) ETR_REC(K).EN_TAIL.LMAC_DIGITAL_TEMP = avg(128)
	    IF (K .EQ. 2) ETR_REC(K).EN_TAIL.LMAC_DIGITAL_TEMP = sigma(128)
	    IF (K .EQ. 3) ETR_REC(K).EN_TAIL.LMAC_DIGITAL_TEMP = min(128)
            IF (K .EQ. 4) ETR_REC(K).EN_TAIL.LMAC_DIGITAL_TEMP = max(128)
	    IF (K .EQ. 5) ETR_REC(K).EN_TAIL.LMAC_DIGITAL_TEMP = nf_rec(128)

	    DUM_REC = ETR_REC(K)
	    call CT_WRITE_ARCV (, ETR_LUN, DUM_REC, CT_STAT)

	    IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	      status = ct_stat(1)
	      call Lib$Signal(FDQ_CTWRITERR, %val(2), %val(status),
	1	'FDQ_ETR')
	      RETSTAT = ERROR
            ELSE
              ETR_WRIT = ETR_WRIT + 1 
	    ENDIF

	   endif	! Retstat is success
	  ENDDO		! K = 1, num_stats_types

	ENDIF		! n_rec .GT. 0

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_FORM_STAT = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_FORM_STAT = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END
