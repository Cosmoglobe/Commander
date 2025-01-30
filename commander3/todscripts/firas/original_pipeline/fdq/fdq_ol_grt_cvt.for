	INTEGER*4 FUNCTION FDQ_OL_GRT_CVT (CNTS, RES_CNTS, CALSET,
	1	 NPTS, RESIS_TBL, TEMP_TBL, QUAD_COEF, GRT_TEMP, 
	1	 STATUS, MJF, telm_flag,GRTIDX)
C/
C/	PROGRAM NAME:
C/	  FDQ_OL_GRT_CVT
C/	
C/	PROGRAM DESCRIPTION:
C/	  This program will convert a GRT value from counts to temperature
C/	  given the count, the GRT lookup-table, the quadratic coefficients
C/	  to fit the table, and the counts of the 4 calibrator resistors
C/	  correponding to the GRT.
C/	
C/	AUTHOR:
C/	  E. Fung
C/	  GSFC
C/	  October 25, 1985
C/	
C/	MODIFIED BY:
C/	  E.FUNG
C/	  GSFC
C/	  November 6, 1985
C/	  REASON:	Add new return status for situation when none of the
C/			cal. resistors is in range.  (Flagged by a negative
C/			resistance value from COUNTS_TO_OHMS.)
C/
C/	  E.FUNG
C/	  GSFC
C/	  December 4, 1985
C/	  REASON:	To change the calling sequence to COUNTS_TO_OHMS for
C/			the new version of COUNTS_TO_OHMS.
C/	
C/	  E.FUNG
C/	  GSFC
C/	  DEC. 7, 1985
C/	  REASON:	To change the writing to the log file to D-lines.
C/			Also, a new input parameter LOGGING to indicate whether
C/			any logging is desired at all.
C/
C/	  E.FUNG
C/	  GSFC
C/	  DEC. 11, 1985
C/	  REASON:	To change to the "database driven" scheme in
C/			getting the cal resistors resistance values
C/			instead of hard-coding them in COUNTS_TO_OHMS.
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
C/	  February 19, 1988
C/	  REASON:       Data being used for testing has some GRTs inactive.
C/			To avoid a long series of error messages, each message
C/			related to GRTs will only be output once. This may
C/			be changed later when GRTs are fixed. 
C/
C/	  Shirley M. Read
C/	  STX
C/	  May 17, 1988
C/	  REASON:       The logging option was never fully implemented.
C/			Other reporting functions replaced logging. 
C/			Several FDQ text files were moved to the FUT Facility.
C/			The include files were renamed to FUT.
C/
C/	  Shirley M. Read
C/	  STX
C/	  August 10, 1988
C/	  REASON:       Removed the open and read of the calibrator resisitor
C/			file since the file has been put under configuration
C/			in a COBETRIEVE time tagged dataset which is now
C/			opened for use in FDQ via CCT_Get_Config. The calibrator
C/			resistor information is now passed in the FUT_Convbuff
C/			include file.
C/
CH	CHANGE LOG:		New Format for Build 4.2 STX, 10/15/88
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
CH      Version 4.2.whatever 03/28/89, SPR 3482, R. Kummerer, STX
CH		Initialize variable STATUS to SUCCESS such that its status
CH		is not remembered from the last call.
CH
CH	6 September 1991, SPR 8955, L. Rosen, STX
CH		Removal of grt_correction routine.
CH
C----------------------------------------------------------------------------------
C/
C/	CALLING SEQUENCE:
C/	  RETSTATUS = FDQ_OL_GRT_CVT (CNTS, RES_CNTS, CALSET, NPTS, RESIS_TBL,
C/		TEMP_TBL, QUAD_COEF, GRT_TEMP, STATUS, MJF, GRTIDX)
C/	
C/	INPUT PARAMETERS:
C/	  CNTS			I*2	Value of GRT in counts
C/	  RES_CNTS(4)		R*4	Counts of the 4 cal. corresponding resistors
C/	  CALSET		I*2	Calibrator resistor group,
C/					1 = A low, 2 = A high, 3 = B low, 4 = B high
C/	  NPTS			I*2	No. of points in lookup table
C/	  RESIS_TBL(NPTS)	R*4	Resistance entries in lookup table
C/	  TEMP_TBL(NPTS)	R*4	Temperature entries in lookup table
C/	  QUAD_COEF(3, NPTS)	R*4	Quadratic coefficients corresponding to
C/					table entries
C/	  MJF			I*2     Major frame of HKP record
C/	  GRTIDX  		I*4     GRT array index for GRT_Names
C/
C/	OUTPUT PARAMETERS:
C/	  GRT_TEMP		R*4	Converted and corrected GRT temperature
C/	  RETSTAT               I*4     Function value for return from lookup
C/					1 = GOOD
C/					17 = exceeds table bounds
C/					18 = accuracy tolerance in fit not satisfied
C/					19 = no cal. resistors in range to calibrate;
C/					     this is manifested by a negative resistance
C/					     returned from COUNTS_TO_OHMS.
C/	  STATUS			I*2	Return status
C/					1 = GOOD
C/					17 = exceeds table bounds
C/					18 = accuracy tolerance in fit not satisfied
C/					19 = no cal. resistors in range to calibrate;
C/					     this is manifested by a negative resistance
C/					     returned from COUNTS_TO_OHMS.
C/
C/	INPUT/OUTPUT FILES:
C/	  NONE
C/
C/	INCLUDE FILES:
C/	  FUT_CONVBUFF.TXT
C/	  FUT_FIRNAMES.TXT
C/
C/	SUBROUTINES CALLED:
C/	  FUT_COUNTS_TO_OHMS
C/	  GRT_LOOKUP
C/	
C/	METHOD USED:
C/
C/	  The following is the PDL --
C/	  
C/	  Call FUT_COUNTS_TO_OHMS using GRT count, Cal resistors' counts, 
C/		Cal resistors' resistances, etc.;
C/	  Call GRT_LOOKUP with input parameters	GRT resistance from (COUNTS_TO_OHMS)
C/		& lookup table;
C/	  Return;
C/	  End.

	IMPLICIT	NONE

	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_NOCALRESEX
	EXTERNAL	FDQ_GRTEXCTBL
	EXTERNAL	FDQ_GRTOLNOSAT

	INCLUDE		'(FUT_CONVBUFF)'
	INCLUDE		'(FUT_FIRNAMES)'
	INCLUDE		'($SSDEF)'

	INTEGER*2	CALSET, CNTS, NPTS,
	1	        STATUS
        Real*4      	RES_CNTS(4)
	integer*2 	MJF
        Logical*1       telm_flag
	integer*4	GRTIDX

	REAL*4		GRT_TEMP, QUAD_COEF(3, *),
	1		RESIS_TBL(*), TEMP_TBL(*)

	INTEGER*4	OTS$CNVOUT	! System function to convert fp to Ascii
	INTEGER*4	GRT_LOOKUP

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			  !
!     Local variables     !
!			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

	CHARACTER*8	CALRES_NAME(4, 4)
	1		/'FACR1TLO', 'FACR2TLO',
	1		'FACR3TLO', 'FACR4TLO',
	1		'FACR1THI', 'FACR2THI',
	1		'FACR3THI', 'FACR4THI',
	1		'FBCR1TLO', 'FBCR2TLO',
	1		'FBCR3TLO', 'FBCR4TLO',
	1		'FBCR1THI', 'FBCR2THI',
	1		'FBCR3THI', 'FBCR4THI' /

	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status
	integer*4	STATUS4		! System  return status
	integer*4       ZERO  / 0 /

	INTEGER*2	A_HIGH, A_LOW,
	1		B_HIGH, B_LOW,
	1		EXC_TBL_BDS, GOOD,
	1		HIGH, LOW, I,
	1		mode /2/,		! 2 for off-line
	1		NO_CAL_RESIS,
	1		TOL_NOT_SAT

	INTEGER*4	GRT_CNTS_I4,
	1		ICUR, IERR, IST 

	REAL*4		DR, GRT_RESIS, TOL,RES_CNTS_I4(4)

	CHARACTER*16	ASCGRTVAL	! Floating point value conv. to Ascii 
	CHARACTER*16    ASCDR		! Floating point value conv. to Ascii

	PARAMETER	(A_HIGH	=   2,
	1		 A_LOW	=   1,
	1		 B_HIGH	=   4,
	1		 B_LOW	=   3,
	1		 EXC_TBL_BDS = 17,
	1		 GOOD	=  1,
	1		 HIGH	=   2,
	1		 LOW	=   1,
	1		 NO_CAL_RESIS = 19,
	1		 TOL_NOT_SAT = 18 )

	DATA		TOL /3.E-4/		! Tolerance for the
						! uncertainty in quadratic
						! fit solutions
C	Temporary flags to signal messages only once.

	logical*1       first_nocalres /.true./,
	1		first_outable /.true./,
	1		first_notok /.true./

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set return status to success.

	RETSTAT = SUCCESS
	STATUS = SUCCESS

	GRT_CNTS_I4 = ZEXT(CNTS)		! Convert GRT_CNTS to
						! I*4 because Rich's routine
						! COUNTS_TO_OHMS expect an I*4
	DO I = 1, 4
	  RES_CNTS_I4(I) = RES_CNTS(I)		! Ditto for resistors' counts
	ENDDO

	IF (CALSET .EQ. A_HIGH .OR. CALSET .EQ. B_HIGH) THEN
	  ICUR = HIGH
	ELSE
	  ICUR = LOW
	ENDIF

!
!     Now call COUNTS_TO_OHMS using the new calling sequence:
!
	CALL FUT_COUNTS_TO_OHMS (RES_CNTS_I4, GRT_CNTS_I4, 1,
	1	CALLO(1, CALSET), CALHI(1, CALSET),
	1	CALRES(1, CALSET), GRT_RESIS)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								!
!     Now that we have the GRT resistance in GRT_RESIS		!
!     call GRT_LOOKUP to find the corresponding temperature     !
!								!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IF (GRT_RESIS .NE. -1.) THEN
	  CALL GRT_LOOKUP (mode,	! Mode = 2 for off-line (1 for on-line)
	1		GRT_RESIS,		! Input resistance
	1		RESIS_TBL,		! Resistance table entries
	1		TEMP_TBL,		! Temperature table entries
	1		QUAD_COEF,		! Quadratic fits coefficients
	1		NPTS,			! No. of table entries
	1		TOL,			! Tolerance for uncertainty; = 3.E-4
	1		GRT_TEMP,		! Temperature of GRT
	1		DR,			! Estimated uncertainty
	1		IERR )			! Status of call
	  STATUS = SUCCESS	! This is reset to override the status from the
				! last call of FDQ_OL_CVT_GRT.  This does not
				! mean the status of GRT_LOOKUP is SUCCESS as
				! IERR reflects its status and is properly
				! interpreted later.
	ELSE
	  STATUS = NO_CAL_RESIS
	  GRT_TEMP = -9999.0		! Give flag value for temperature
	  if ( first_nocalres ) then
	    CALL LIB$SIGNAL(FDQ_NOCALRESEX, %VAL(2), GRT_NAMES(GRTIDX),
	1	%VAL(MJF))
	    first_nocalres = .false.
	  endif
	ENDIF


	IF ( RETSTAT .EQ. SUCCESS .AND. STATUS .NE. NO_CAL_RESIS ) THEN
	  IF (IERR .EQ. 1) THEN
C	    if ( first_outable ) then 
C ---------------- check telm quality
	      STATUS4 = OTS$CNVOUT(GRT_RESIS,ASCGRTVAL,%VAL(5))
	      IF ( .NOT.STATUS4) THEN
	        ASCGRTVAL = 'OTS Conv. Error'
	      ENDIF
              If ( .not.telm_flag) then
       	        CALL LIB$SIGNAL(FDQ_GRTEXCTBL, %VAL(3), GRT_NAMES(GRTIDX),
	1	%VAL(MJF), ASCGRTVAL)
              Endif
	      first_outable = .false.
C	    endif
	    grt_temp = -9999.0
	    STATUS = EXC_TBL_BDS
	  ELSEIF (IERR .EQ. 2) THEN
C	    if ( first_notok) then 
C ---------------- check telm quality
	      STATUS4 = OTS$CNVOUT(GRT_TEMP,ASCGRTVAL,%VAL(5))
	      IF ( .NOT.STATUS4) THEN
	        ASCGRTVAL = 'OTS Conv. Error'
	      ENDIF 
	      STATUS4 = OTS$CNVOUT(DR,ASCDR,%VAL(5))
	      IF ( .NOT.STATUS4) THEN
	        ASCGRTVAL = 'OTS Conv. Error'
	      ENDIF 
              If (.not.telm_flag) then
	        CALL LIB$SIGNAL(FDQ_GRTOLNOSAT, %VAL(4), GRT_NAMES(GRTIDX),
	1	  %VAL(MJF), ASCGRTVAL, ASCDR)
              Endif
	      first_notok = .false.
C	    endif
	    STATUS = TOL_NOT_SAT
	  ELSEIF (IERR .EQ. 0) THEN
	    STATUS = GOOD
	  ENDIF
	ENDIF		! Return status is success and calres is available

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_OL_GRT_CVT = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_OL_GRT_CVT = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END

