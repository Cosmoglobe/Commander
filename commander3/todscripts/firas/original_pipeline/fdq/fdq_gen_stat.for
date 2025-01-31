	INTEGER*4 FUNCTION FDQ_GEN_STAT (ETR_LUN, ENG_REC, ORBITAL_PD,ETR_WRIT)

C-----------------------------------------------------------------------------
C/	PURPOSE:  Calculate Avg, Sigma, Min, Max, for TRENDPLOTS for
C/		  each orbit and write into FDQ_ETR (Eng Stat Rec).
C/
C/	PROGRAM NAME:
C/	  FDQ_GEN_STAT
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine handles the engineering statistics.
C/
C/	AUTHOR:
C/	  E.FUNG
C/	  GSFC
C/	  DEC. 16, 1985
C/
C/
C/	CALLING SEQUENCE FROM FDQ:
C/	  Return status = FDQ_GEN_STAT (ETR_LUN, ENG_REC, ORBITAL_PD)
C/
C/	INPUT PARAMETERS:
C/	  ETR_LUN		I*2	logical unit no. for Eng. Trends
C/	  ENG_REC(1024)		BYTE	1024-byte record buffer to be
C/					output to FIRAS Eng. Archive
C/	  ORBITAL_PD            R*4     Time for one orbit in minutes
C/
C/	OUTPUT PARAMETERS:
C/        ETR_WRIT              I*4     Number of new ETR records written
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
C/	  CT_GMT_TO_BINARY
C/	  CT_BINARY_TO_GMT
C/
C/	FUNCTIONS CALLED:
C/	  FDQ_STAT_INIT
C/	  FDQ_FORM_STAT
C/
C/	ERROR HANDLING:
C/	  Call Lib$Signal
C/
C/	METHOD USED:
C/	  PDL by Doug Ward  October 28, 1987
C/
C/	  SAVE current start time, stop time, and eng stat counters
C/	  Convert ASCII period to binary
C/	  IF (FIRST_RECORD) THEN
C/	   *Set FIRST_RECORD FALSE
C/	   *CALL FDQ_STAT_INIT to set stop_time and start_time, and to
C/	    initialize min, max, sum, sum_sq, element1, nf_rec and n_rec
C/	  ELSEIF (time > stop_time)
C/	 *CALL FDQ_FORM_STAT to compute avg, sigma, and avg_time and
C/	  write these to ETR_LUN
C/	 *See if an orbit or more has been skipped and report skips
C/	 *CALL FDQ_STAT_INIT to reset stop_time and start_time, and to
C/	  re-initialize min, max, sum, sum_sq, element1, nf_rec and n_rec
C/	  ENDIF
C/	  IF (not first record in orbital period) then
C/	 *Increment n_rec
C/	 *Write time to rec_time
C/	  ENDIF
C/	 *DOFOR each parameter (128)
C/	   If the first field value is not set and field is not flagged, then
C/	     Set the first field value and increment field record count.
C/	   Endif
C/	   If the first field value is set, then
C/	     Accumulate min, max, sum, sum_sq
C/	   Endif
C/	 *ENDDOFOR
C/	RETURN
C-----------------------------------------------------------------------------
CH	CHANGE LOG:
CH	  E.FUNG
CH	  GSFC
CH	  DEC. 27, 1985
CH	  REASON:	(1) Left out a "%REF" in 2 LIB$MOVC3 calls for
CH			    character variables, this version to correct;
CH			(2) Program bomb with floating point overflow when
CH			    trying to calculate the square.  Put in check
CH			    to see if number is greater than 1.3E+19 (the
CH			    square root of 1.7E+38, the biggest floating
CH			    point possible in VAX floating form numbers)
CH			    before squaring it.
CH
CH	  E.FUNG
CH	  GSFC
CH	  DEC. 28, 1985
CH	  REASON:	Caught a typo for the first executable statement:
CH			should be ENG_REC(195) instead of ENG_REC(105)!
CH
CH			Also make sure that the system binary time fields
CH			in both the hourly and daily buffer are there; use
CH			system time when writing to the archive rather than
CH			the start of the interval.
CH
CH	  E.Fung
CH	  GSFC
CH	  Jan. 6, 1986
CH	  Reason:	To correct the bug about not resetting the sums
CH			and the sums of squares after a statistics record
CH			has been written.
CH
CH			Also, to avoid "taking the square root of a neg.
CH			number" due to truncation errors by checking
CH			the number inside the sq. root sign.  If the
CH			number is sufficiently small whether positive
CH			or negative, the RMS will be set to zero.
CH
CH			Also, to correct bug in calling LIB$MOVC5 (formerly
CH			mistakenly called with only 4 arguments instead of 5).
CH
CH	  E. FUNG
CH	  GSFC
CH	  JAN. 6, 1986
CH	  REASON:	To go to 2 engineering statistics files (hourly and
CH			daily being separated); requires one more input parameter
CH			(HES_LUN and DES_LUN replacing EST_LUN).
CH
CH	  E. FUNG
CH	  GSFC
CH	  FEB. 25, 1986
CH	  REASON:	To adjust to the new ENG archive format (analog ENG
CH			values start at offset 252 now instead of 194).
CH
CH	  D. WARD
CH	  GSFC/STX
CH	  OCT. 6, 1987
CH	  REASON:	To change computations to orbital increments in place
CH			of hourly and daily increments.  ETR_ replaces HES_ and
CH			DES_.  Data structure ENG_REC is also used now.
CH
CH	  D. WARD
CH	  GSFC/STX
CH	  OCT. 27, 1987
CH	  REASON:	To update FDQ_GEN_STAT to produce only one statistical
CH			set of quantities every orbit, roughly.  Function
CH			FDQ_STAT_INIT was added to trim the module as well
CH			as the maintenance.
CH
CH	  D.WARD
CH	  GSFC/STX
CH	  NOV. 13, 1987
CH	  REASON:	To add variable diff in order to avoid negative
CH			SQRT radicands made possible by round off error
CH			when all elements are equal.
CH
CH	  Shirley M. Read
CH	  STX
CH	  January 15, 1988
CH	  REASON:	Converted for interface with the Fut_Error
CH			condition handler. Added error checking
CH		        and calls to Lib$Signal.
CH
CH	Version 4.2.1 02/09/89, SPR 3139, Shirley M. Read, STX
CH		FDQ does not weed out the bad data flag value of -9999.0 when
CH		it forms the engineering trends records. Modules affected may be
CH		FDQ_Init_Stats, FDQ_Gen_Stats and FDQ_Form_Stats. In computing
CH		the statistics, the flag value should not be included in the
CH		number, sum or sum of squares.
CH
CH	Version 4.2.1 02/09/89, SPR 2955, Shirley M. Read, STX
CH		FDQ should be given the correct COBE orbital period which is
CH		used to form engineering trends. Make the orbital period a
CH		command line option.
CH
CH	Version 4.5,  1989 Oct 10,  SPR 4670,  Fred Shuman,  STX.
CH	    FIT, which reads and plots FDQ_ETR records, one for each of five
CH	    (previously four) types of statistics, was making garbled plots.
CH	    Most of this turned out to be in FIT and was caused by code there
CH	    which did not recognize the fifth ETR record, but when this was
CH	    cleared up, there remained some bad data in the ETR records written
CH	    by FDQ.  FDQ_FORM_STAT and FDQ_GEN_STAT.
C-----------------------------------------------------------------------------

	IMPLICIT	NONE

	DICTIONARY 'FDQ_ENG'
	RECORD /FDQ_ENG/ ENG_REC

	INCLUDE		'CT$LIBRARY:CTUSER.INC/NOLIST'

	INCLUDE		'(FDQ_STATBUFF)'

	INCLUDE		'($SSDEF)'

	LOGICAL*1	FIRST_RECORD /.TRUE./

	EXTERNAL	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL        FDQ_FRSTSTATER
	EXTERNAL	FDQ_STATINITER
	EXTERNAL        FDQ_FORMSTATER
	EXTERNAL	FDQ_ERBINTIM
	EXTERNAL	FDQ_ENGANLGDIF

	INTEGER*2	ETR_LUN
	REAL*4	        ORBITAL_PD	! Time for one orbit in minutes
        INTEGER*4       ETR_WRIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			  !
!     Local Variables     !
!			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

	LOGICAL*1	TIME_GT

	INTEGER*4	FDQ_STAT_INIT, FDQ_FORM_STAT

	INTEGER*4	I, istat, BAD, SYS$BINTIM

	INTEGER*4	start_time(2), stop_time(2), period(2), n_orbits

	INTEGER*4	k, n_loops, n_modulo

	REAL*4		diff, element
	REAL*8          R8_Orb_Pd	! R*8 orbital period
	Logical*1	first_orb_rec / .true. /
	Logical*1       flag		! Eng field flagged or not
	REAL*4		FV / -9999.0 /  ! Engineering bad value flag

	integer*4	RETSTAT		! Return status
	integer*4	SUCCESS / 1 /, ERROR / 2 /  ! Values for status
	integer*4	ZERO   / 0 /
	integer*2       one / 1 /
	PARAMETER	( BAD	= 2 )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!SAVE for next loop thru
!

	SAVE FIRST_RECORD, period, start_time, stop_time

C	Set return status to success.

	RETSTAT = SUCCESS

	IF (FIRST_RECORD) THEN		! only for the first record sent

	  FIRST_RECORD = .FALSE.

	  R8_Orb_Pd = Orbital_Pd * 60.0 * 1.D7		! Cnvt to 100-nsec units
	  Call Aut_Dfloat2Adt ( R8_Orb_Pd, Period )     ! Convert to quadword

	  if ( retstat .eq. success ) then
	    DO k = 1, 3
	      DO i = 1, 100
	        weights(i, k) = 1
	      ENDDO
	    ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!				!
!    1st CALL FDQ_STAT_INIT	!
!				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    istat = FDQ_STAT_INIT (ENG_REC, period, start_time, stop_time)
	    if ( istat .ne. %LOC(FDQ_NORMAL) ) then
	      if (istat .eq. %loc(FDQ_ABERR)) istat = zero
	      call lib$signal(FDQ_FRSTSTATER,%val(1),%val(istat))
	      retstat = error
	    else
	      first_orb_rec = .true.
	    endif

	  endif	! Status is success

	ELSEIF ( TIME_GT(eng_rec.ct_head.TIME, stop_time) ) THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						     !
!		CALL FDQ_FORM_STAT		     !
!						     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  istat = FDQ_FORM_STAT (ETR_LUN,ETR_WRIT)
	  if ( istat .ne. %LOC(FDQ_NORMAL) ) then
	      if (istat .eq. %loc(FDQ_ABERR)) istat = zero
	      call lib$signal(FDQ_FORMSTATER,%val(1),%val(istat))
	      retstat = error
	  endif

!
! reset stop_time, start_time and process this current file
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!				!
!    2nd CALL FDQ_STAT_INIT	!
!				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  if ( retstat .eq. success ) then
	    istat = FDQ_STAT_INIT (ENG_REC, period, start_time, stop_time)
	    if ( istat .ne. %LOC(FDQ_NORMAL) ) then
	      if (istat .eq. %loc(FDQ_ABERR)) istat = zero
	      call lib$signal(FDQ_STATINITER,%val(1),%val(istat))
	      retstat = error
	    endif
	    first_orb_rec = .true.
	  endif	! Retstat is success

	ENDIF		! (FIRST_RECORD)
!			!  ELSEIF (eng_rec.ct_head.TIME .GT. stop_time)

!!!!!!!!!!!!!!!!!!!!!!!!!
!			!
!	ACCUMULATE	!
!			!
!!!!!!!!!!!!!!!!!!!!!!!!!

	IF ( retstat .eq. success ) THEN
	  If ( .not. first_orb_rec ) then
	    n_rec = n_rec + 1
	    n_loops = (n_rec - 1)/100
	    n_modulo = n_rec - n_loops*100
	    rec_time(1, n_modulo, n_loops +1) = eng_rec.ct_head.TIME(1)
	    rec_time(2, n_modulo, n_loops +1) = eng_rec.ct_head.TIME(2)
	  Endif
          Do i = 1, 64
	    flag = .false.
	    If (eng_rec.en_analog.GRT(i) .ne. FV ) then
	      element = eng_rec.en_analog.GRT(i)
	      If (nf_rec(i) .eq. zero) then
	        element1(i) = element
	      Endif
	      nf_rec(i) = nf_rec(i) + 1
	    Else
	      flag = .true.
	    Endif
	    If ( .not. flag ) then
	      If ( nf_rec(i) .eq. one ) then
		min(i) = element
		max(i) = element
	      Else
	        IF (element .GT. max(i)) max(i) = element
	        IF (element .LT. min(i)) min(i) = element
	        diff = element - element1(i)
	        sum(i) = sum(i) + diff
	        IF (abs(diff) .gt. 1.3e19) then
	          call lib$signal(FDQ_ENGANLGDIF)
	          sum_sq(i) = sum_sq(i) + diff
	        ELSE
	          sum_sq(i) = sum_sq(i) + diff * diff
	        ENDIF
	      Endif      ! nf_rec(i) .eq. one
	    Endif	! not flag
          Enddo

	  DO I = 1, 62
	    flag = .false.
	    If (eng_rec.en_analog.GROUP1(i) .ne. FV ) then
	      element = eng_rec.en_analog.GROUP1(i)
	      If (nf_rec(i+64) .eq. zero) then
	        element1(i+64) = element
	      Endif
	      nf_rec(i+64) = nf_rec(i+64) + 1
	    Else
	      flag = .true.
	    Endif
	    If ( .not. flag ) then
	      If ( nf_rec(i+64) .eq. one ) then
		min(i+64) = element
		max(i+64) = element
	      Else
	        IF (element .GT. max(i+64)) max(i+64) = element
	        IF (element .LT. min(i+64)) min(i+64) = element
	        diff = element - element1(i+64)
	        sum(i+64) = sum(i+64) + diff
	        IF (abs(diff) .gt. 1.3e19) then
	          call lib$signal(FDQ_ENGANLGDIF)
	          sum_sq(i+64) = sum_sq(i+64) + diff
	        ELSE
	          sum_sq(i+64) = sum_sq(i+64) + diff * diff
	        ENDIF
	      Endif      ! nf_rec(i+64) .eq. one
	    Endif	! not flag
          Enddo
	  flag = .false.
	  If (eng_rec.en_tail.lmac_analog_temp .ne. FV ) then
	     element = eng_rec.en_tail.lmac_analog_temp
	     If (nf_rec(127) .eq. zero) then
	            element1(127) = element
	     Endif
	     nf_rec(127) = nf_rec(127) + 1
	  Else
	    flag = .true.
	  Endif
	  If ( .not. flag ) then
	    If ( nf_rec(127) .eq. one ) then
	      min(127) = element
	      max(127) = element
	    Else
	      IF (element .GT. max(127)) max(127) = element
	      IF (element .LT. min(127)) min(127) = element
	      diff = element - element1(127)
	      sum(127) = sum(127) + diff
	      IF (abs(diff) .gt. 1.3e19) then
	        call lib$signal(FDQ_ENGANLGDIF)
	        sum_sq(127) = sum_sq(127) + diff
	      ELSE
	        sum_sq(127) = sum_sq(127) + diff * diff
	      ENDIF
	    Endif      ! nf_rec(127) .eq. one
	  Endif	! not flag
          Flag = .false.
	  If (eng_rec.en_tail.lmac_digital_temp .ne. FV ) then
	     element = eng_rec.en_tail.lmac_digital_temp
	     If (nf_rec(128) .eq. zero) then
	          element1(128) = element
	     Endif
	     nf_rec(128) = nf_rec(128) + 1
	  Else
	     flag = .true.
	  Endif
	  If ( .not. flag ) then
	    If ( nf_rec(128) .eq. one ) then
	      min(128) = element
	      max(128) = element
	    Else
	      IF (element .GT. max(128)) max(128) = element
	      IF (element .LT. min(128)) min(128) = element
	      diff = element - element1(128)
	      sum(128) = sum(128) + diff
	      IF (abs(diff) .gt. 1.3e19) then
	        call lib$signal(FDQ_ENGANLGDIF)
	        sum_sq(128) = sum_sq(128) + diff
	      ELSE
	        sum_sq(128) = sum_sq(128) + diff * diff
	      ENDIF
	    Endif      ! nf_rec(128) .eq. one
	  Endif	! not flag

	  first_orb_rec = .false.
	ENDIF		! Retstat is success and accumulate

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_GEN_STAT = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_GEN_STAT = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END
