	SUBROUTINE FGA_LIST_SCI (WRITE_LUN, SCI_REC, NREC, LINE,
	1	CAT_ENTRY, WRITE_IFG)

C-----------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FGA_LIST_SCI
C
C	PROGRAM DESCRIPTION:
C	  The program will get a formatted listing of all the header
C	  info in a given FIRAS Science record.
C
C	AUTHOR:
C	  R. Kummerer
C	  STX
C	  November 13, 1987
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_SCI (WRITE_LUN, SCI_REC, NREC, LINE,
c		CAT_ENTRY, WRITE_IFG)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file.
C	  SCI_REC		RECORD	FIRAS Science record
C	  NREC			I*2	Count of the number of records in this run
C	  LINE			C*80	Comment line to identify the list file
C	  CAT_ENTRY		I*4	CT catalog entry #
C	  WRITE_IFG		L*1	Flag to WRITE out IFG as well as header info
C
C	OUTPUT PARAMETERS:
C	  NONE
C
C	INPUT FILES:
C	  NONE
C
C	OUTPUT FILES:
C	  Listing file or FOR006
C
C	INCLUDE FILES USED:
C	  NONE
C
C	SUBROUTINES CALLED:
C	  BINARY_TO_GMT
C
C
C-----------------------------------------------------------------------
C
C Changes:
C
C	SPR 2789, November 21, 1988, R. Kummerer.  Channel indictor
C	not displayed properly.
C
C	SPR 2700, February 22, 1988, Shirley M. Read. Getarchive is 
C	needed in the testing of the time tag fix for the FIRAS science
C	archive. The NFS_SDF Rdl has been updated to include new fields 
C	for the midpoint of collect time and badtime flag. These fields
C	must be dumped in order to verify that the FIRAS stripper and 
C	preprocessor are performing correctly.
C
C	SPR 3442, March 22, 1989, R. Kummerer. Dump data quality flags.
C
C       Add fields SC_HEAD23, HSKP1_TLM_FMT, HSKP2_TLM_FMT (SPR 4234)
C                  N. Gonzales, STX 7-17-90
C       SPR 8413, May 8, 1991, N. Gonzales. Corrected MIN_FRM_CNT equation
C                  by swapping fields SCI_HEAD12 and SCI_HEAD13. 
C                  DATA_TRANS_TIME is also corrected by swapping SC_HEAD4 and
C                  SC_HEAD5 within the equation.
C-----------------------------------------------------------------------

	IMPLICIT	NONE

	CHARACTER	*80	LINE

	INTEGER		*2	NREC
	INTEGER		*4	WRITE_LUN
	INTEGER		*4	CAT_ENTRY
	INTEGER		*4	STATUS

	CHARACTER	*14	GMT
	BYTE			WRITE_IFG

	INTEGER		*4	SW_VERSION_EQV
	BYTE			SW_VERSION(2)
	INTEGER		*4	CHAN_IND_EQV
	BYTE			CHAN_IND(2)

	INTEGER		*4	SC5
	INTEGER		*4	SC13

	INTEGER		*2	I
	INTEGER		*2	J
	INTEGER		*2	K
	INTEGER		*2	L

	INTEGER		*4	SUM
	INTEGER		*4	UNSIGNED_CHECKSUM
	INTEGER		*4	MIN_FRM_CNT
	INTEGER		*4	DATA_TRANS_TIME

	CHARACTER	*512	GLITCH
	CHARACTER       *14     COLLECT_GMT

	DICTIONARY	'NFS_SDF'
	RECORD /NFS_SDF/ SCI_REC

	EXTERNAL		FGA_WRITE_ERR

	EQUIVALENCE	(SW_VERSION_EQV, SW_VERSION)
	EQUIVALENCE	(CHAN_IND_EQV, CHAN_IND)

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6

	WRITE (WRITE_LUN, 6000, IOSTAT=STATUS) LINE, NREC, CAT_ENTRY

6000	FORMAT ('1', / 5X, A // 
	1	' Records number in this run: ', I5,
	1	10x, 'CT Catalog Entry: ', i7 //
	1	15X, 'FIRAS Science Record Formatted Dump' /)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Cobetrieve standard header.
C
	IF ((SCI_REC.COLLECT_TIME.MIDPOINT_TIME(1) .EQ. 0) .AND.
	1   (SCI_REC.COLLECT_TIME.MIDPOINT_TIME(2) .EQ. 0)) THEN
	   COLLECT_GMT = '00000000000000'
	ELSE
	   CALL CT_BINARY_TO_GMT (SCI_REC.COLLECT_TIME.MIDPOINT_TIME,
	1	COLLECT_GMT)
	ENDIF

	WRITE (WRITE_LUN, 6001, IOSTAT=STATUS)
	1			(SCI_REC.CT_HEAD.GMT(I:I), I=1,14),
	1			SCI_REC.CT_HEAD.TIME(2),
	1			SCI_REC.CT_HEAD.TIME(1),
	1			SCI_REC.CT_HEAD.SPACE_TIME,
	1			SCI_REC.CT_HEAD.MJR_FRM_NO,
	1			(COLLECT_GMT(I:I), I=1,14),
	1			SCI_REC.CT_HEAD.ORBIT,
	1                       SCI_REC.CT_HEAD.HSKP1_TLM_FMT,
	1                       SCI_REC.CT_HEAD.HSKP2_TLM_FMT,
	1			SCI_REC.COLLECT_TIME.MIDPOINT_TIME(2),
	1			SCI_REC.COLLECT_TIME.MIDPOINT_TIME(1),
	1			SCI_REC.COLLECT_TIME.BADTIME_FLAG

6001	FORMAT (/' Processed time (in 3 formats): ' /
	1	' <GMT> : ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	' <Binary> : ', Z8.8, 1X, Z8.8, 7X,
	1	' <PB5> : ',  6(Z2.2, 1X) /
	1	' Major Frame Number : ', I6, 8x, 
	1	'<Midpoint of Collect> : ',2A1,'-',3A1,'-',3(2A1,'-'),3A1/
	1	' Orbit number : ', I6, /
	1       ' TLM FMT Major Frame 1: ',10X,I4 /
	1       ' TLM FMT Major Frame 2: ',10X,I4 / 
	1	'<Binary Collect> : ', Z8.8, 1x, Z8.8, 6x,' Time Flag: ',I4)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Science record header.
C
	SW_VERSION_EQV = SCI_REC.SCI_HEAD.SC_HEAD2

	UNSIGNED_CHECKSUM = ZEXT(SCI_REC.SCI_HEAD.SC_HEAD6)

	SUM = 0
	DO J=1,512
	  SUM = SUM + SCI_REC.IFG_DATA.IFG(J)
	ENDDO

	SC13 = SCI_REC.SCI_HEAD.SC_HEAD12
	IF (SC13 .LT. 0) THEN
	   SC13 = SC13 + 65536
	END IF

	MIN_FRM_CNT = 65536 * SCI_REC.SCI_HEAD.SC_HEAD13 + SC13

	SC5 = SCI_REC.SCI_HEAD.SC_HEAD4
	IF (SC5 .LT. 0) THEN
	   SC5 = SC5 + 65536
	END IF

	DATA_TRANS_TIME = 65536 * SCI_REC.SCI_HEAD.SC_HEAD5 + SC5

	CHAN_IND_EQV = SCI_REC.SCI_HEAD.SC_HEAD25


	WRITE (WRITE_LUN, 6002, IOSTAT=STATUS)
	1			SCI_REC.SCI_HEAD.CHAN_ID,
	1			SCI_REC.SCI_HEAD.GAIN,
	1			SCI_REC.SCI_HEAD.MTM_SPEED,
	1			SCI_REC.SCI_HEAD.MTM_LENGTH,
	1			SCI_REC.SCI_HEAD.DATA_QUAL,
	1			SCI_REC.SCI_HEAD.DATA_READY,
	1			SCI_REC.SCI_HEAD.SC_HEAD0,
	1			SCI_REC.SCI_HEAD.SC_HEAD1A,
	1			SCI_REC.SCI_HEAD.SC_HEAD1B,
	1			(SW_VERSION(I),I=2,1,-1),
	1			SCI_REC.SCI_HEAD.SC_HEAD3,
	1			SCI_REC.SCI_HEAD.SC_HEAD4,
	1			SCI_REC.SCI_HEAD.SC_HEAD5,
	1			DATA_TRANS_TIME,
	1			SCI_REC.SCI_HEAD.SC_HEAD6,
	1			UNSIGNED_CHECKSUM,
	1			SUM,
	1			SCI_REC.SCI_HEAD.SC_HEAD7,
	1			SCI_REC.SCI_HEAD.SC_HEAD8,
	1			SCI_REC.SCI_HEAD.SC_HEAD9,
	1			SCI_REC.SCI_HEAD.SC_HEAD10,
	1			SCI_REC.SCI_HEAD.SC_HEAD11,
	1			SCI_REC.SCI_HEAD.SC_HEAD12,
	1			SCI_REC.SCI_HEAD.SC_HEAD13,
	1			MIN_FRM_CNT,
	1			SCI_REC.SCI_HEAD.SC_HEAD14,
	1			SCI_REC.SCI_HEAD.SC_HEAD15,
	1			SCI_REC.SCI_HEAD.SC_HEAD16,
	1			SCI_REC.SCI_HEAD.SC_HEAD17,
	1			SCI_REC.SCI_HEAD.SC_HEAD18,
	1			SCI_REC.SCI_HEAD.SC_HEAD19,
	1			SCI_REC.SCI_HEAD.SC_HEAD20,
	1			SCI_REC.SCI_HEAD.SC_HEAD21,
	1			SCI_REC.SCI_HEAD.SC_HEAD22,
	1	                SCI_REC.SCI_HEAD.SC_HEAD23,
	1			(CHAN_IND(I),I=2,1,-1)

6002	FORMAT (//'***** Science Header Information *****'/
	1	' Channel Number: ', T30, I12 /
	1	' Gain: ', T30, I12 /
	1	' MTM Speed: ', T30, I12 /
	1	' MTM Length: ', T30, I12 //
	1	' Telemetry Quality Flags: '  /
	1	'   1 -  10: ', 5X, 10(I4, 5X) /
	1	'  11 -  20: ', 5X, 10(I4, 5X) /
	1	'  21 -  30: ', 5X, 10(I4, 5X) /
	1	'  31 -  40: ', 5X, 10(I4, 5X) /
	1	'  41 -  50: ', 5X, 10(I4, 5X) /
	1	'  51 -  60: ', 5X, 10(I4, 5X) //
	1	' Data Ready Flags (HEX): ' /
	1	'   1 -   8: ', 5X, 8(2X, Z4.4, 5X) //
	1	' Data Block Sync (HEX): ', T38, Z4.4 /
	1	' Block ID: ', T30, I12/
	1	' Block Type: ', T41, A1/
	1	' Software Version Number: ', T40, 2A1 /
	1	' Status Bits (HEX): ', T38, Z4.4 /
	1	' Transmit Time Counter LSW: ', T30, I12 /
	1	' Transmit Time Counter MSW: ', T30, I12 /
	1	' Data Transmit Time: ', T30, I12 /
	1	' Data Checksum (HEX): ', T38, Z4.4 /
	1	' Data Checksum (Unsigned): ', T30, I12 /
	1	' Calculated Checksum: ', T30, I12 /
	1	' Number of A/D Samples: ', T30, I12 /
	1	' Data Points Processed: ', T30, I12 /
	1	' Adds per Group: ', T30, I12 /
	1	' Data Points per Mirror Sweep: ', T30, I12 /
	1	' Number of Mirror Sweeps: ', T30, I12 /
	1	' Frame Counter LSW: ', T30, I12 /
	1	' Frame Counter MSW: ', T30, I12 /
	1	' Minor Frame Count: ', T30, I12 /
	1	' Deglitcher Threshold Factor: ', T30, I12 /
	1	' Deglitcher Noise Level: ', T30, I12 /
	1	' Sampling Counter Constant: ', T30, I12 /
	1	' Deglitcher Speed Sample: ', T30, I12 /
	1	' Command Counter: ', T30, I12 /
	1	' Program Checksum (HEX): ', T38, Z4.4 /
	1	' Saturated Sample Count: ', T30, I12 /
	1	' Total Glitch Count: ', T30, I12 /
	1	' Deglitcher Overflow Address: ', T30, I12 /
	1       ' Division of IFG Performed: ',T30, I12 /
	1	' Channel Indicator: ', T40, 2A1)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	CALL CT_BINARY_TO_GMT (SCI_REC.DQ_DATA.ENG_TIME, GMT)

	WRITE (WRITE_LUN, 6010, IOSTAT=STATUS)
	1			SCI_REC.DQ_DATA.FAKE,
	1			SCI_REC.DQ_DATA.XCAL_POS,
	1			SCI_REC.DQ_DATA.IREF_TEMP,
	1			(GMT(I:I),I=1,14),
	1			SCI_REC.DQ_DATA.ENG_TIME(2),
	1			SCI_REC.DQ_DATA.ENG_TIME(1),
	1			SCI_REC.DQ_DATA.ENG_REC,
	1			SCI_REC.DQ_DATA.DATA_QUALITY(110),
	1			SCI_REC.DQ_DATA.DATA_QUALITY(109),
	1			SCI_REC.DQ_DATA.DATA_QUALITY

6010	FORMAT (//'***** Data Qualify Information *****' /
	1	' Fake-it Mode: ', T30, I6 /
	1	' XCal Position: ', T30, I6 /
	1	' ICal Temperature: ', T26, G14.5 /
	1	' Engineering Time: ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	' <Binary>: ', Z8.8, 1X, Z8.8 /
	1	' Engineering Record Number: ', T30, I6 /
	1	' Instrument Data Quality: ', T30, I6 /
	1	' Attitude Data Quality: ', T30, I6 /
	1	' Data Quality Flags: ' //
	1	'   1 -  10: ', 5X, 10(I4, 5X) /
	1	'  11 -  20: ', 5X, 10(I4, 5X) /
	1	'  21 -  30: ', 5X, 10(I4, 5X) /
	1	'  31 -  40: ', 5X, 10(I4, 5X) /
	1	'  41 -  50: ', 5X, 10(I4, 5X) /
	1	'  51 -  60: ', 5X, 10(I4, 5X) /
	1	'  61 -  70: ', 5X, 10(I4, 5X) /
	1	'  71 -  80: ', 5X, 10(I4, 5X) /
	1	'  81 -  90: ', 5X, 10(I4, 5X) /
	1	'  91 - 100: ', 5X, 10(I4, 5X) /
	1	' 101 - 110: ', 5X, 10(I4, 5X))

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Dump the attitude information.
C
	CALL FGA_LIST_ATT (WRITE_LUN, SCI_REC.ATTITUDE)

C
C IFG and glitch data.
C
	GLITCH=' '
	K=0
	DO I=1,32
	  DO J=15,0,-1
	    K=K+1
	    IF (BTEST(SCI_REC.IFG_DATA.GLTCH(I),J)) THEN
	      GLITCH(K:K)='*'
	    ENDIF
	  ENDDO
	ENDDO


	WRITE (WRITE_LUN, 6020, IOSTAT=STATUS)

6020	FORMAT (//'***** Onboard Glitch Map *****'/)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6025, IOSTAT=STATUS) SCI_REC.IFG_DATA.GLTCH

6025	FORMAT (1X,8Z6.4)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	IF (WRITE_IFG) THEN

	  WRITE (WRITE_LUN, 6200, IOSTAT=STATUS) (J,J=1,10)

6200	  FORMAT ('1 ***** Science Interferogram *****'/
	1	  '             (* = Glitch)'
	1	   ///18X, 10(I2, 10X))

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF


	  DO J=0,50
	    K = J*10
	    WRITE (WRITE_LUN, 6101, IOSTAT=STATUS) K,
	1	(SCI_REC.IFG_DATA.IFG(J*10+L),GLITCH(K+L:K+L),L=1,10)
	    IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	    END IF
	  ENDDO

	  K = 510

	  WRITE (WRITE_LUN, 6101, IOSTAT=STATUS) K,
	1		SCI_REC.IFG_DATA.IFG(511), GLITCH(511:511),
	1	 	SCI_REC.IFG_DATA.IFG(512), GLITCH(512:512)

6101	  FORMAT (1X, I7, 10(2X, I8, 1X, A1))

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

	ENDIF

	RETURN
	END
