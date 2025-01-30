	INTEGER*4 FUNCTION FDQ_OL_POLY_CVT ( NBYTES, VALUE, NC, COEF, 
	1	CVALUE, MJF, CURPOLY) 
C/
C/	PROGRAM NAME:
C/	  FDQ_OL_POLY_CVT
C/
C/	PROGRAM DESCRIPTION:
C/	  This subroutine take the raw count of a data point as input
C/	  and return the value in engineering units as output assuming
C/	  the curve type is a polynomial.
C/
C/	AUTHOR:
C/	  EDWIN FUNG
C/	  GSFC
C/	  October 25, 1985
C/	  Converted from POLYCVT.
C/
C/	MODIFIED BY:
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 7, 1988
C/	  REASON: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal.
C/
C/	  Shirley M. Read
C/	  STX
C/	  June 3, 1988
C/	  REASON: 	Several FDQ text files were moved to the FUT Facility.
C/			The include files were renamed to FUT.
C/
CH      CHANGE LOG: 	New Format for Build 4.2  STX, 10/15/88
CH
CH      Version 4.1.1 10/15/88, SPR 2622, Shirley M. Read, STX
CH              Flag values for missing conversions need changes in FDQ.
CH	        The flag values for converted fields which are out of range
CH		need the same changes. GRTs in dwell mode need a flag value
CH              also. The SWG has selected a flag value of -9999 for all cases.
CH              The limit checking algorithms need to bypass setting the quality
CH              flags if the GRT or engineering analog has the flag value.
CH		Interpolation and averaging algorithms need to be modified to
CH              include the flag value. 
CH
C-----------------------------------------------------------------------------
C/
C/	CALLING SEQUENCE:
C/	  STATUS = FDQ_OL_POLY_CVT ( NBYTES, VALUE, NC, COEF, CVALUE, 
C/		   MJF, CURPOLY) 
C/
C/	INPUT PARAMETERS:
C/
C/	  NBYTES	I*2	The no. of bytes the data point occupies.
C/	  VALUE(*)	L*1	The value of the data point in raw counts.
C/	  NC		I*2	Number of polynomial coefficients.
C/	  COEF(6)	R*4	Polynomial coefficients in descending order;
C/				up to 6 (or 5th degree poly).
C/				COEF(1) is the coefficent of the highest
C/				degree term, NOT the 5th degree term.
C/	  MJF		I*2     Major frame of HKP record
C/	  CURPOLY       I*2     POLY_COEF array index of current HKP 
C/				engineering parameter being converted to E.U. 
C/
C/	OUTPUT PARAMETERS:
C/	  CVALUE	R*4	The value of the data point in engineering units.
C/
C/	ERROR HANDLING:
C/	  NONE
C/
C/	SUBROUTINES CALLED:
C/	  LIB$MOVC3
C/	  LIB$POLYF
C/
C/	MACRO/INCLUDE FILES:
C/	  FUT_CONVBUFF.TXT
C/	  FUT_FIRNAMES.TXT
C/
C/	METHOD USED:
C/	  The following is the PDL --
C/
C/	  Adjust number of bytes according to NC
C/	  Call LIB$POLYF to get polynomial conversion
C/

	IMPLICIT	NONE

	include 	'($SSDEF)'

	LOGICAL*1	VALUE(*)

	INTEGER*2	NBYTES, NC

	REAL*4		COEF(6), CVALUE

	INTEGER*2	MJF, CURPOLY

	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_POLYCVTER

	INCLUDE	        '(FUT_CONVBUFF)'
	INCLUDE		'(FUT_FIRNAMES)'

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			  !
!     Local variables     !
!			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

	integer*4       LIB$POLYF	! System lib function for poly conv
	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status
	integer*4       STATUS

	INTEGER*2	I2, NDEG

	INTEGER*4	I4

	REAL*4		X

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
! 			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set return status to success.

	RETSTAT = SUCCESS

	I4 = 0				! Zero out argument for polynomial

        IF (NBYTES .EQ. 1) THEN	! Unsigned 1-byte integer
	  I4 = ZEXT(VALUE(1))
	ELSEIF (NBYTES .EQ. 2) THEN	! Signed 2-byte integer
	  CALL LIB$MOVC3(2,VALUE,I2)
	  I4 = I2
	ELSE 			! Signed 4-byte integer
	  CALL LIB$MOVC3(4,VALUE,I4)
	ENDIF
	X = FLOAT(I4)
	NDEG = NC - 1
	IF ( NDEG .LT. 0 ) THEN
	  CVALUE = -9999.0
	ELSE
	  STATUS = LIB$POLYF (X,		! Argument for the polynomial
	1		NDEG,			! Degree of polynomial
	2		COEF,	! Array of coefficients (high order term first)
	3		CVALUE	)		! The result of the evaluation of the polynomial

	  IF ( STATUS .NE. SS$_Normal ) THEN
	      CALL LIB$SIGNAL(FDQ_POLYCVTER,%VAL(3),%VAL(STATUS),
	1	%VAL(MJF), POLY_NAMES(CURPOLY))
	      CVALUE = -9999.0
!	      RETSTAT = ERROR	! Remove error status and set cvalue to 0
	  ENDIF
	ENDIF	! Ndeg LT 0

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_OL_POLY_CVT = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_OL_POLY_CVT = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END

