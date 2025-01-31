  	INTEGER*4 FUNCTION FEP_GRT_LOOKUP (GIVEN, IDIR, OHMS, TEMPS,
	2                                  first,LEN, RESULT, DR)

C--------------------------------------------------------------------------
C
C              GRT TABLE LOOKUP SUBROUTINE
C
C     Written 9 Oct 84 by Rich Isaacman (Applied Research Corp.)
C     Modified 16 Oct 84 to go in both directions
C
C     Modified 21-Jun-1989 by Don Stevens-Rayburn (Applied Research Corp.)
C	in an attempt to speed it up.  This version only handles conversions
C	from ohms to temperature because the search technique relys on the
C	fact that resistance is a monotonically increasing function of table
C	index.  The calling parameter IDIR is retained for compatibility but
C	this version ignores it.  In addition, this routine does not recompute
C	polynomial fitting coefficients if the GRT resistance is in the same
C	interval as during the previous invocation, because the old ones are
C	still correct.
C
C
c  Receives GRT resistance or temperature and searches appropriate lookup
c  table using a binary search. Having found the table entries that
c  bracket the resistance (temperature), it interpolates the temperature
c  (resistance) lookup table with a 2ND-degree polynomial.
c  The polynomial solution is done twice, once with table entries XOXX and
c  once with entries XXOX, where X represents the GRT lookup table entries 
c  and O represents the measured GRT value. The difference in the
c  derived solution in the two cases is used to estimate the accuracy
c  of the result. The input parameter IDIR tells the routine whether it is
c  go from resistance to temperature or vice versa.
C
C  Input Parameters:  GIVEN = Input value, either temperature or resistance
C
C		      IDIR = 0 if GIVEN = resistance;
C			   = 1 if GIVEN = temperature
C
C		      OHMS = Lookup table of resistances
C
C		      TEMPS = Lookup table of temperatures
C
C		      LEN = Length of lookup tables
C		     
C  Output Parameters: RESULT = Temperature if IDIR = 0;
C			     = Resistance if IDIR = 1
C
C		      DR = Estimated uncertainty in result in Kelvins or Ohms
C
C  Average cpu usage for this subroutine is 0.9 millisec per solution; average
C  accuracy is 0.15 millikelvins for a 200-entry lookup table.
C
C--------------------------------------------------------------------------
C
C Changes:
C
C	Treat FEP_POLYSLNERR as an informational error instead of fatal.
C	November 4, 1987, R. Kummerer.
C
C	SPR 1507, Signal TBLEDGE message once per plot.  R. Kummerer,
C	December 11, 1987.
C
C	SPR 2910, Plot all fields dwelled in the timerange rather than just
C	the one with the plurality of hits on each side (A/B).  Fred Shuman,
C	STX,  1989 Feb 17.
C
C	SPR 5130, Spurious temperature jumps in plots 
C	Harte Wang, STX, 1989 dec. 8
C--------------------------------------------------------------------------

	IMPLICIT	NONE

	REAL		*4	OHMS(400)
	REAL		*4	TEMPS(400)
	REAL		*4	POLYL(3)	! Polynomial coeffs for lower fit.
	REAL		*4	POLYH(3)	! Polynomial coeffs for upper fit.
	REAL		*4	SOLN(2)		! Two solutions.
        INTEGER		*4	IL  / 1 /	! Table index below GIVEN.
	INTEGER		*4	IU  / 1 /	! Table index above GIVEN.
	INTEGER		*4	LEN
	INTEGER		*4	ITEST
	INTEGER		*4	IDIR
        logical         *1      first
	REAL		*4	DR
	REAL		*4	GIVEN
	REAL		*4	RESULT
	INTEGER		*4	STATUS

	INTEGER		*4	FEP_INVERT

	EXTERNAL	FEP_NORMAL
	EXTERNAL	FEP_TBLEDGE
	EXTERNAL	FEP_POLYSLNERR

	FEP_GRT_LOOKUP = %LOC( FEP_NORMAL )

C
C  Return error if independent variable is at or beyond edge of table.
C
        if (first) then 
         IU = 1
         IL = 1
        Endif

	IF ( GIVEN .LT. OHMS ( 2 ) ) THEN
	   RESULT = TEMPS ( 1 )
	   DR = 999.
	   FEP_GRT_LOOKUP = %LOC( FEP_TBLEDGE )
	   RETURN
	ELSE IF ( GIVEN .GT.  OHMS ( LEN - 1 ) ) THEN
	   RESULT = TEMPS ( LEN )
	   DR = 999.
	   FEP_GRT_LOOKUP = %LOC( FEP_TBLEDGE )
	   RETURN
	ENDIF
        IF ( ( GIVEN .LT. OHMS ( IL ) ) .OR.		! Have we moved out of
	2    ( GIVEN .GT. OHMS ( IU ) )     ) THEN	! an interval?
							! If so, find the new 
	   IL = 1					! interval and compute
	   IU = LEN					! new coefficients.
           ITEST = ( IU + IL ) / 2

	   DO WHILE ( ITEST .NE. IL )			! Standard Binary search.
              IF ( GIVEN .GT. OHMS ( ITEST ) ) THEN
		 IL = ITEST
	      ELSE
	         IU = ITEST
	      END IF
              ITEST = ( IU + IL ) / 2
	   END DO

C  We can no longer quit early because we may need the coeffiecients later.
C  This is no big loss because the probability that the resistance is exactly
C  equal to a table entry is negligible.
C
C  Quit early if input value exactly equals a table entry
C
C	   IF ( GIVEN .EQ. OHMS ( IL ) ) THEN
C	      RESULT = TEMPS ( IL )
C	      DR = 0.
C	      RETURN
C	   ENDIF
C
C
C  Get the 3 points that bracket the value of the independent variable;
C  send these three xy pairs to the polynomial solver to get interpolation
C  coefficients.
C
C  First fit two points below and one above.

	   STATUS = FEP_INVERT (POLYL, OHMS( IL-1 ), TEMPS(IL-1) )
	   IF (STATUS .NE. %LOC(FEP_NORMAL)) THEN
	      FEP_GRT_LOOKUP = STATUS
	      RETURN
	   ENDIF

C  Now fit one below and two above.

	   STATUS = FEP_INVERT (POLYH, OHMS(IL), TEMPS(IL) )
	   IF (STATUS .NE. %LOC(FEP_NORMAL)) THEN
	      FEP_GRT_LOOKUP = STATUS
	      RETURN
	   ENDIF

	END IF

C  Whether we have recomputed the coefficients or not the arrays POLYL and POLYH
C  now contain the "correct" coefficients.

	SOLN(1) = POLYL(1) + POLYL(2)*GIVEN + POLYL(3)*GIVEN*GIVEN
	SOLN(2) = POLYH(1) + POLYH(2)*GIVEN + POLYH(3)*GIVEN*GIVEN

C  Average the two solutions and estimate the possible error from their
C  difference.

	RESULT = (SOLN(1) + SOLN(2))/2.	               !Average 2 sol'ns
	DR = ABS(SOLN(1) - SOLN(2))/2.	               !Estimate error from
	IF (DR/RESULT .GT. 6.E-4) THEN
	   CALL LIB$SIGNAL(FEP_POLYSLNERR)
	END IF

	RETURN

	END
