	INTEGER*4 FUNCTION FDQ_CLOSE_ARCV (SCI_OPEN, CT_LUN, 
	1					ENGLIM)
C/
C/	PROGRAM NAME:
C/	  FDQ_CLOSE_ARCV
C/	
C/	PROGRAM DESCRIPTION:
C/	  This routine closes all archive files used in DATA QUALIFY.
C/
C/	AUTHOR:
C/	  EDWIN H. FUNG
C/	  GSFC
C/	  OCT. 25, 1985
C/
C/	MODIFIED BY:
C/	  E. FUNG
C/	  GSFC
C/	  JAN. 6, 1986
C/	  REASON:	To add an additional parameter as the eng. statistics
C/			went from 1 to 2 CT files (hourly and daily).
C/
C/	  E. FUNG
C/	  GSFC
C/	  FEB. 26, 1986
C/	  REASON:	To close 4 Science files instead of one (due to new
C/			Science formats).
C/
C/	  E. FUNG
C/	  GSFC
C/	  JULY 23, 1987
C/	  REASON:	To use 40 bytes for CT status return (formerly CT return
C/			was only to an I*4 variable; will have access violation
C/			in new CT).
C/
C/	  D. WARD
C/	  GSFC/STX
C/	  OCTOBER 20, 1987
C/	  REASON:	Replace HES AND DES with ETR.
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 6, 1988
C/	  REASON: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal.
C/
C/	  R. Kummerer
C/	  STX
C/	  June 27, 1988
C/	  REASON: 	Science Catalog / Matrix archive closes.
C/
C/        H. Wang
C/        STX
C/        Jan. 29/91
C/        Reason:     New requirements for FDQ
C/                    Tracking file is no longer be required.
C/                    Close The FEX_AV_CALRS Archive
C/         
CH	Version 4.4.1 07/22/89, SER 4168, R. Kummerer, STX
CH		There have been problems during test operations for FIRAS
CH		processing due to the required clean-up of the archives after
CH		an FPP or FDQ abort. The FPR tracking system compounds the
CH		problems. Files with non-matching version numbers seem often
CH		to result from improper clean-up. Bad record times cause
CH		SEGCTL to abort and mess up the tracking system. It was
CH		decided to change the modify of the science records in FPP
CH		and FDQ to a simple COBETRIEVE read of the existing records
CH		from a dataset and write a modifed dataset with the same
CH		information which was entered on the modify. Two new science
CH		data sets will be required: a science dataset of raw science
CH		data plus FPP input and a science dataset with FPP science
CH		data plus FDQ input. These datasets will be FPP_SDF_xx, where
CH		xx is the channel id (RH, RL, LH or LL) and FDQ_SDF_xx, where
CH		xx is the channel id. The new datasets must be opened and
CH		processed in FPP and FDQ. 
CH
CH	Version 4.4.3 11/14/89, SPR 5032, R. Kummerer STX
CH		Skip processing raw science segments from missing channels.
CH	Version 5.0 12/11/89, SPR 5313, R. Kummerer STX
CH		Perform appropriate raw science and IFG tracker archive closes.
CH
C/	CALLING SEQUENCE:
C/	  STATUS = FDQ_CLOSE_ARCV (SCI_OPEN, CT_LUN,  ENGLIM)
C/
C/	INPUT PARAMETERS:
C/	  SCI_OPEN(4)	L*1	Flags indicating raw science archive opens
C/	  CT_LUN(22)	I*2	Logical units for open archives
C/	  ENGLIM	I*4	Engineering limits violation counts.
C/
C/	OUTPUT PARAMETERS:
C/	  NONE
C/
C/	INPUT/OUTPUT FILES:
C/	  NONE
C/
C/	INCLUDE FILES USED:
C/	  CT$LIBRARY:CTUSER.INC
C/
C/	SUBROUTINES CALLED:
C/	  CT_CLOSE_ARCV
C/
C/	ERROR HANDLING:
C/	  Call Lib$Signal for errors and set return status for function.
C/
C/	METHOD USED:
C/	  Do for all channels
C/	    Close Science Archive;
C/	  Enddo;
C/	  Close Index Archive;
C/	  Close Housekeeping Archive;
C/	  Close Ave. Cal. Resistor Archive;
C/	  Close Engineering Archive;
C/	  Close Eng. Stat. Archive;
C/	  Return;
C/	  End.
C/
	IMPLICIT	NONE

	INCLUDE		'(FUT_PARAMS)'


	INTEGER*2	CT_LUN(22)

	INCLUDE		'CT$LIBRARY:CTUSER.INC'

	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_CTCLOSERR
	EXTERNAL	FDQ_RMSWRITECTL
	EXTERNAL	FDQ_RMSCLOSECTL

	INTEGER*4 	RETSTAT		! Return status
	INTEGER*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status
	INTEGER*4	ZERO / 0 /
	INTEGER*4	STATUS		! Dummy status variable
	INTEGER*4	IOSTATUS	! Dummy status variable
	INTEGER*2	CHAN, CT_ISTAT(20)
	INTEGER*4	ENGLIM
	LOGICAL*1	SCI_OPEN(4)

	INTEGER*2	ETR, ENG, IDX,CAL,
	1		HSK, RS1, RS2, RS3, RS4, RPT,
	2		SM1, SM2, SM3, SM4,
	3		SC1, SC2, SC3, SC4,
	4		ORS1, ORS2, ORS3, ORS4

	PARAMETER	(RS1 = 1)
	PARAMETER	(RS2 = 2)
	PARAMETER	(RS3 = 3)
	PARAMETER	(RS4 = 4)
	PARAMETER	(HSK = 5)
	PARAMETER	(IDX = 6)
	PARAMETER	(ENG = 7)
	PARAMETER	(ETR = 8)
	PARAMETER	(CAL = 9)
	PARAMETER       (RPT = 10)  ! Ct_Lun 10 is reserved for RMS report file
	PARAMETER	(SM1 = 11)
	PARAMETER	(SM2 = 12)
	PARAMETER	(SM3 = 13)
	PARAMETER	(SM4 = 14)
	PARAMETER	(SC1 = 15)
	PARAMETER	(SC2 = 16)
	PARAMETER	(SC3 = 17)
	PARAMETER	(SC4 = 18)
	PARAMETER	(ORS1 = 19)
	PARAMETER	(ORS2 = 20)
	PARAMETER	(ORS3 = 21)
	PARAMETER	(ORS4 = 22)

	CHARACTER*12	DATASETS(22)
	DATA DATASETS / 'FPP_SDF_RH',
	1		'FPP_SDF_RL',
	1		'FPP_SDF_LH',
	1		'FPP_SDF_LL',
	1		'NFS_HKP',
	1		'FDQ_IDX',
	1		'FDQ_ENG',
	1		'FDQ_ETR',
	1		'FEX_AV_CALRS',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		'FDQ_SDF_RH',
	1		'FDQ_SDF_RL',
	1		'FDQ_SDF_LH',
	1		'FDQ_SDF_LL' /



!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	Set return status to success.

	RETSTAT = SUCCESS

	DO	CHAN = 1, 4

	  IF ( SCI_OPEN(CHAN) ) THEN

	    CALL CT_CLOSE_ARCV (, CT_LUN(CHAN), CT_ISTAT)
	    IF (CT_ISTAT(1) .NE. CTP_NORMAL) THEN
	      RETSTAT = ERROR
	      STATUS = CT_ISTAT(1)
	      CALL LIB$SIGNAL(FDQ_CTCLOSERR,%VAL(2),%VAL(STATUS),
	1	DATASETS(CHAN))
	    ENDIF

	    CALL CT_CLOSE_ARCV (, CT_LUN(CHAN+18), CT_ISTAT)
	    IF (CT_ISTAT(1) .NE. CTP_NORMAL) THEN
	      RETSTAT = ERROR
	      STATUS = CT_ISTAT(1)
	      CALL LIB$SIGNAL(FDQ_CTCLOSERR,%VAL(2),%VAL(STATUS),
	1	DATASETS(CHAN+18))
	    ENDIF


	  ENDIF

	ENDDO

	IF ( CT_LUN(IDX) .NE. 0 ) THEN
	  CALL CT_CLOSE_ARCV (, CT_LUN(IDX), CT_ISTAT)
	  IF (CT_ISTAT(1) .NE. CTP_NORMAL) THEN
	    RETSTAT = ERROR
	    STATUS = CT_ISTAT(1)
	    CALL LIB$SIGNAL(FDQ_CTCLOSERR,%VAL(2),%VAL(STATUS),
	1	DATASETS(IDX))
	  ENDIF
	ENDIF

	IF ( CT_LUN(HSK) .NE. 0 ) THEN
	  CALL CT_CLOSE_ARCV (, CT_LUN(HSK), CT_ISTAT)
	  IF (CT_ISTAT(1) .NE. CTP_NORMAL) THEN
	     RETSTAT = ERROR
	     STATUS = CT_ISTAT(1)
	     CALL LIB$SIGNAL(FDQ_CTCLOSERR,%VAL(2),%VAL(STATUS),
	1	DATASETS(HSK))
	  ENDIF
	ENDIF

	IF ( CT_LUN(ENG) .NE. 0 ) THEN
	  CALL CT_CLOSE_ARCV (, CT_LUN(ENG), CT_ISTAT)
	  IF (CT_ISTAT(1) .NE. CTP_NORMAL) THEN
	     RETSTAT = ERROR
	     STATUS = CT_ISTAT(1)
	     CALL LIB$SIGNAL(FDQ_CTCLOSERR,%VAL(2),%VAL(STATUS),
	1	DATASETS(ENG))
	  ENDIF
	ENDIF

	IF ( CT_LUN(ETR) .NE. 0 ) THEN
	  CALL CT_CLOSE_ARCV (, CT_LUN(ETR), CT_ISTAT)
	  IF (CT_ISTAT(1) .NE. CTP_NORMAL) THEN
	     RETSTAT = ERROR
	     STATUS = CT_ISTAT(1)
	     CALL LIB$SIGNAL(FDQ_CTCLOSERR,%VAL(2),%VAL(STATUS),
	1	DATASETS(ETR))
	  ENDIF
	ENDIF
	IF ( CT_LUN(CAL) .NE. 0 ) THEN
	  CALL CT_CLOSE_ARCV (, CT_LUN(CAL), CT_ISTAT)
	  IF (CT_ISTAT(1) .NE. CTP_NORMAL) THEN
	     RETSTAT = ERROR
	     STATUS = CT_ISTAT(1)
	     CALL LIB$SIGNAL(FDQ_CTCLOSERR,%VAL(2),%VAL(STATUS),
	1	DATASETS(CAL))
	  ENDIF
	ENDIF

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_CLOSE_ARCV = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_CLOSE_ARCV = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END
