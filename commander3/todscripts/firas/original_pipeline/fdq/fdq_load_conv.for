
 	INTEGER*4 FUNCTION FDQ_LOAD_CONV( Config_Lun, New_Segment )
C/
C/	PROGRAM NAME:
C/	  FDQ_LOAD_CONV
C/
C/	PROGRAM DESCRIPTION:
C/	  This subroutine reads the database disk files and loads the
C/	  polynomial coefficients, the GRT lookup tables and calibrator 
C/	  resistorinformation into the common area CONVBUFF, for 
C/	  conversions of Housekeeping data to engineering units later
C/	  (in subroutine FDQ_CONVERT).
C/
C/	AUTHOR:
C/	  EDWIN H. FUNG
C/	  GSFC
C/	  OCTOBER 23, 1985
C/
C/	MODIFIED BY:
C/	  E.FUNG
C/	  GSFC
C/	  NOV. 3, 1985
C/	  REASON:	Add write statements to unit 16 (log file).
C/			Also, move freeing and closing of LIMCUVKEY.DB
C/			until the GRT coefficients have been read in.
C/
C/	  E. FUNG
C/	  GSFC
C/	  NOV. 12, 1985
C/	  REASON:	To add the option of writing a brief ASCII conversions
C/			file.
C/
C/	  E. FUNG
C/	  NOV. 14, 1985
C/	  REASON:	To correct the bug on reading GRTRICH.DB (using GRT_NAMES(1)
C/			instead of GRT_NAMES(I)).
C/
C/	  E. FUNG
C/	  DEC. 7, 1985
C/	  REASON:	To change the writing to the log file to D-lines.
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 15, 1988
C/	  REASON: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal. Modified calling sequence,
C/		        added the missing parameters. These are user options.
C/
C/	  Shirley M. Read
C/	  STX
C/	  June 3, 1988
C/	  REASON: 	Several FDQ text files were moved to the FUT Facility.
C/		        The include files were renamed to FUT.
C/
C/	  Shirley M. Read
C/	  STX
C/	  August 10, 1988
C/	  REASON: 	The conversion coefficients, GRT tables and calibrator
C/			resistor information which used to be read from the
C/			DB$ RMS files are now under COBETRIEVE configuration
C/			as time tagged datasets. The time tagged datasets which
C/			match the data time will now be opened and made 
C/			available to FDQ via CCT_Get_Config. Since these are
C/			indexed files, FDQ_Load_Conv must still read each file
C/			whenever the segment is changed. The reading of the
C/			calibrator resistor information has been moved to this
C/			routine from FDQ_OL_GRT_Cvt.
C/
C/	  H. Wang
C/	  STX
C/	  Jan. 29, 1991
C/	  REASON: 	New requirements for FDQ
C/                      The routine FDQ_PRINT_XYB will no longer
C/                      be called.
C/    
C/	CALLING SEQUENCE:
C/	  Return_Status = FDQ_LOAD_CONV(Config_Lun, New_Segment)
C/
C/	INPUT PARAMETERS:
C/	  Config_Lun(3)   I*4  Unit numbers for reading conversion information
C/			       Unit 1: Conversion coefficients
C/			       Unit 2: GRT tables
C/			       Unit 3: Calibrator resistors
C/	  New_Segment(3)  L*1  Flag indicating if a new file of time tagged
C/			       data should be read
C/
C/	OUTPUT PARAMETERS:
C/	  The common area CONVBUFF is filled
C/
C/	INPUT FILES:
C/
C/	INCLUDE FILES USED:
C/	  FUT_CONVBUFF.TXT
C/	  FUT_FIRNAMES.TXT
C/	  FUT_GRTPARS.TXT
C/	  FUT_LCKEYPARS.TXT
C/	  SYS$LIBRARY:FORIOSDEF
C/
C/	SUBROUTINES USED:
C/	  BINARY_TO_GMT (COBETRIEVE ROUTINE)
C/	  LIB$GET_LUN (from system library)
C/	  LIB$MOVC3 (from system library)
C/	  MVBITS (from system library)
C/	  SYS$GETTIM (from system library)
C/
C/	ERROR HANDLING:
C/	  Passed back in output parameter ISTAT
C/
C/	METHOD USED:
C/	 The following is the PDL --
C/
C/
C/	 If a new segment of time tagged conversions is to be read;
C/	  Do until all polynomial coefficients are loaded
C/	    Read LIMCUVKEY.DB with database name as key;
C/	    Move coefficients into CONVBUFF;
C/	    From 3-word-index move NBYTES into VALSIZE;
C/	  Enddo;
C/	 Endif;
C/
C/	 If a new segment of time tagged GRT tables are to be read:
C/	  Do until all GRT lookup tables are loaded
C/	    Read GRTCOEFF.DB with database name as key;
C/	    Move lookup table and quadratic coefficients to CONVBUFF;
C/	 Endif;
C/
C/	 If a new segment of time tagged calibrator resistors is to be read;
C/        Do for all 4 cal sets
C/	    Do for all 4 cal resistors
C/	      Read configured file FDB_GRTCAL using cal resistors database name
C/	      Remember cal resistance values as well as valid lower and upper
C/	      counts;
C/	    Enddo;
C/	  Enddo;
C/	 Endif;
C/
C/	  Return;
C/	  End.
C/

	IMPLICIT	NONE

	INCLUDE		'(FUT_CONVBUFF)'
	INCLUDE		'(FUT_FIRNAMES)'
	INCLUDE		'(FUT_GRTPARS)'
	INCLUDE		'(FUT_LCKEYPARS)'
	INCLUDE		'SYS$LIBRARY:FORIOSDEF/NOLIST'
	INCLUDE		'($SSDEF)'

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         !
!	Passed Parameters !
! 			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

	INTEGER*4	CONFIG_LUN(3)		!Units for configuration files
	LOGICAL*1	NEW_SEGMENT(3)          !New time tagged file flags

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			  !
!     Local Variables     !
!			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!


	BYTE		GRT_BUFF(GRT_RECSZ),
	1		LCK_BUFF(LCK_RECSZ)

	INTEGER*2	I, II, J, NOC_EQV, TWI_EQV(3)

	INTEGER*4	 CALGRP_EQV,
	1		 IST, NPTS_EQV

	REAL*4		CFF_EQV(6)

	EQUIVALENCE	(CFF_EQV,	LCK_BUFF(LCK_CFF))	! Poly. Coef.
	EQUIVALENCE	(NOC_EQV,	LCK_BUFF(LCK_NOC))	! Number of Coef.
	EQUIVALENCE	(TWI_EQV,	LCK_BUFF(LCK_TWI))	! 3-word index

	EQUIVALENCE	(CALGRP_EQV,	GRT_BUFF(GRT_DB_SET))	! Cal group in DB file
	EQUIVALENCE	(NPTS_EQV,	GRT_BUFF(GRT_DB_KNCT))	! # points in lookup table

	LOGICAL*1	BUFF20(20)  	! Read buffer

	CHARACTER*8	CALRES_NAME(4, 4)
	1		/'FACR1TLO', 'FACR2TLO',
	1		'FACR3TLO', 'FACR4TLO',
	1		'FACR1THI', 'FACR2THI',
	1		'FACR3THI', 'FACR4THI',
	1		'FBCR1TLO', 'FBCR2TLO',
	1		'FBCR3TLO', 'FBCR4TLO',
	1		'FBCR1THI', 'FBCR2THI',
	1		'FBCR3THI', 'FBCR4THI' /

	INTEGER*2	CALSET


	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status
	integer*4	STATUS		! Dummy status variable
	integer*4       ZERO / 0 /

	logical*1 	FIRST / .TRUE. /  ! First time in subroutine

	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_ERGETTIM
	EXTERNAL	FDQ_LUNGETERR
	EXTERNAL	FDQ_OPENERR
	EXTERNAL	FDQ_NOPOLYDEF
	EXTERNAL	FDQ_READERR
	EXTERNAL	FDQ_CLOSERR
	EXTERNAL	FDQ_NOGRTLOOK
	EXTERNAL	FDQ_PRINTXYB

	Integer*4 sys$gettim
	Integer*4 lib$get_lun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set return status to success.

	RETSTAT = SUCCESS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									 !
!    The following DO loop will load conversions coefficients		 !
!    for all HSKP quantities' which have polynomial-type conversions     !
!									 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       If ( NEW_SEGMENT(1) ) Then

	DO	I = 1, NPOLY

	 If ( retstat .eq. success ) then

	  READ (CONFIG_LUN(1), KEY=POLY_NAMES(I), IOSTAT=IST) LCK_BUFF
	  IF (IST .NE. 0) THEN
	    IF (IST .EQ. FOR$IOS_ATTACCNON) THEN
	      Call Lib$Signal(FDQ_NOPOLYDEF,%val(1),poly_names(i))

	      POLY_NC(I) = 0
	    ELSE
	        call lib$signal(FDQ_READERR,%val(2),%val(status),
	1	'CSDR$FIRAS_REF:FDB_LIMCUVKY')
	        RETSTAT = ERROR
	    ENDIF
	  ELSE						! GOOD READ
	    POLY_NC(I) = NOC_EQV
	    DO	II=1, NOC_EQV				! Move coef.
	      POLY_COEF(II, I) = CFF_EQV(NOC_EQV-II+1)	! *** NOTE ***
	    ENDDO					! The coef. in DB file
							! are in ascending order;
							! but in POLY_COEF it is
							! in descending order because
							! later on LIB$POLY is expecting
							! coef. in descending order
!
!     Get # bytes in the database word in the 3-word index
!
	    CALL MVBITS (TWI_EQV(2), 6, 7, POLY_VALSIZE(I), 0) 
	  ENDIF		! IST .NE. 0

	 Endif    ! retstat is success
	ENDDO		! I = 1, NPOLY
       Endif	  ! New segment (1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							   !
!     The next DO loop will load the GRT lookup tables     !
!							   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       If ( NEW_SEGMENT(2) ) Then	
	DO	I = 1, NGRT
	 If ( Retstat .eq. success ) then

	  READ (CONFIG_LUN(2), KEY=GRT_NAMES(I), IOSTAT=IST) GRT_BUFF
	  IF (IST .NE. 0) THEN
	    IF (IST .EQ. FOR$IOS_ATTACCNON) THEN
	      call lib$signal(FDQ_NOGRTLOOK,%val(1),grt_names(i))
	      
	      GRT_CALGRP(I) = 0
	    ELSE
	        call lib$signal(FDQ_READERR,%val(2),%val(status),
	1	'CSDR$FIRAS_REF:FDB_GRTRICH')
	        RETSTAT = ERROR

	    ENDIF
	  ELSE						! GOOD READ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								  !
!     Fill in GRT_CALGRP, GRT_TBL with info read from DB file     !
!								  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    GRT_CALGRP(I) = CALGRP_EQV		! Even though CALGRP_EQV
	    GRT_NPTS(I) = NPTS_EQV		! and NPTS_EQV are I*4 and
						! GRT_CALGRP and GRT_NPTS are
						! I*2, it is not anticipated
						! that the Cal group number
						! or the #points in lookup tables
						! will exceed 32767 (or something like that)
	    CALL LIB$MOVC3 (4*NPTS_EQV, 	! Move all resistance values
	1		GRT_BUFF(GRT_DB_KNTS),
	1		GRT_XTAB(1, I))

	    CALL LIB$MOVC3 (4*NPTS_EQV, 	! Move all temperatures values
	1		GRT_BUFF(GRT_DB_KNTS + 4 * NPTS_EQV),
	1		GRT_YTAB(1, I))

	    CALL LIB$MOVC3 (4*3*(NPTS_EQV-2),	! Move quad. coef.; since the 1st
	1					! and last point in the table have no coefficients
	1		GRT_BUFF(GRT_DB_KNTS + 8 * NPTS_EQV),	! Start of coef. in record buffer
	1		GRT_COEF(1, 2, I))	! The middle index starts with
						! 2 again because of the 1st point
						! has no quadratic fit coefficients


	  ENDIF		! IST .NE. 0

	 endif   ! Retstat is success
	ENDDO			! I = 1, NGRT
       Endif	 ! New Segment (2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!	Read the calibrator resistor values from the configuration dataset.  !
!	  CALLO(4,4)            I*2     Low count                            !
!	  CALHI(4,4)            I*2     High count                           !
!	  CALRES(4,4)           R*4     Calibrator resistor values           ! 
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       If ( NEW_SEGMENT(3) ) Then
	DO CALSET = 1, 4	 
	  DO	I = 1, 4
	   IF ( RETSTAT .EQ. SUCCESS ) THEN
	    READ (CONFIG_LUN(3), KEY=CALRES_NAME(I, CALSET), IOSTAT=IST)
	1	BUFF20
	    IF (IST .NE. 0) THEN
	      CALL LIB$SIGNAL(FDQ_READERR,%VAL(2),%VAL(IST),
	1		'CSDR$FIRAS_REF:FDB_GRTCAL')
	      RETSTAT = ERROR
	    ELSE
	      CALL LIB$MOVC3(2, BUFF20(11), CALLO(I, CALSET))
	      CALL LIB$MOVC3(2, BUFF20(13), CALHI(I, CALSET))
	      CALL LIB$MOVC3(4, BUFF20(15), CALRES(I, CALSET))
	    ENDIF
	   ENDIF	! If return status is success
	  ENDDO
	ENDDO		! CALSET = 1, 4


       Endif		! New Segment (3)


C	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_LOAD_CONV = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_LOAD_CONV = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END

