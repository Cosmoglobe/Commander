	SUBROUTINE FGA_LIST_ENG (WRITE_LUN, ENG_REC, NREC, LINE,
	1	CAT_ENTRY)

C------------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FGA_LIST_ENG
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given FIRAS
C	  Engineering record.
C
C	AUTHOR:
C	  E.FUNG
C	  GSFC
C	  November 12, 1985
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_ENG (WRITE_LUN, ENG_REC, NREC, LINE, CAT_ENTRY)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file
C
C	  ENG_REC		RECORD	FIRAS Engineering record.
C
C	  NREC			I*2	The number of records processed in this run.
C
C	  LINE			C*80	Comment line for listing
C
C	  CAT_ENTRY		I*4	CT catalog entry
C
C	OUTPUT PARAMETERS:
C	  NONE
C
C	INPUT FILES:
C	  NONE
C
C	OUTPUT FILES:
C	  Listing file or FOR006.
C
C	INCLUDE FILES USED:
C	  NONE
C
C	SUBROUTINES CALLED:
C	  LIB$MOVC3
C
C------------------------------------------------------------------------
C   Changes:
C	Add LMAC fields--subroutine FGA_LIST_ENGANLG  (SPR 2634)
C	                 F. Shuman, STX, 1988 Oct 20.
C	SPR 3442, Dump data quality flags. R. Kummerer,  March 22, 1989.
C       Add fields IFG_NO, TLM_QUAL_MAJ_FRM, HSKP1_TLM_FMT, HSKP2_TLM_FMT,
C           DWELL_STAT --subroutine FGA_ENGSTAT (SPR 4234)
C                        N. Gonzales, STX, 1990 July 10.
C       Add Fields structure EN_Tempdiff (SPR 8372) 
C                        N. Gonzales/STX September 19, 1991
C------------------------------------------------------------------------

	IMPLICIT	NONE

	CHARACTER	*80	LINE
	INTEGER		*2	NREC
	INTEGER		*4	WRITE_LUN
	INTEGER		*4	CAT_ENTRY
	INTEGER		*4	STATUS
	CHARACTER	*7	REC_STAT
	CHARACTER	*14	GMT(4)
	INTEGER		*2	I
	INTEGER		*2	J

	DICTIONARY	'FDQ_ENG'
	RECORD /FDQ_ENG/ ENG_REC

	EXTERNAL		FGA_WRITE_ERR

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6

	IF (ENG_REC.EN_TAIL.HSKP_FLAG .EQ. 0) THEN
	  REC_STAT = 'Valid'
	ELSE
	  REC_STAT = 'Invalid'
	ENDIF


	WRITE (WRITE_LUN, 6000, IOSTAT=STATUS)	LINE, NREC, CAT_ENTRY, REC_STAT

6000	FORMAT ( '1' / 7X, 'FIRAS ENGINEERING FORMATTED LISTING' /
	1	'***  ', A, '  ***'
	1	/ '   Record Count for this Run: ***', I6, '  ***'
	1	10x, 'CT Catalog Entry: ', I7, 2x, 'Record is: ',
	1	'*** ', a7, ' ***')

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C	COBETRIEVE standard header
C
	WRITE (WRITE_LUN, 6001, IOSTAT=STATUS)
	1			ENG_REC.CT_HEAD.GMT(1:2),
	1			ENG_REC.CT_HEAD.GMT(3:5),
	1			ENG_REC.CT_HEAD.GMT(6:7),
	1			ENG_REC.CT_HEAD.GMT(8:9),
	1			ENG_REC.CT_HEAD.GMT(10:11),
	1			ENG_REC.CT_HEAD.GMT(12:14),
	1			ENG_REC.CT_HEAD.TIME(2),
	1			ENG_REC.CT_HEAD.TIME(1),
	1			ENG_REC.CT_HEAD.SPACE_TIME,
	1			ENG_REC.CT_HEAD.ORBIT,
	1                       ENG_REC.CT_HEAD.HSKP1_TLM_FMT,
	1                       ENG_REC.CT_HEAD.HSKP2_TLM_FMT

6001	FORMAT (/' ASCII GMT:                             ', 4X,
	1	A2, '-', A3, '-', 3(A2, '-'), A3,
	1	5X, ' GMT in binary (high order 4 bytes first): ', 2X,
	1	2(Z8.8,2X)
	1	/' S/C time in PB4 (HEX):                 ', 6X,
	1	6(Z2.2,1X), 5X,
	1	' Orbit Number:                          ',12X,
	1	I11, /
	1       ' TLM FMT Major Frame 1: ',10X,I4 /
	1       ' TLM FMT Major Frame 2: ',10X,I4 ) 

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Limits and Conversions File Names and SCI times
C
	  WRITE (WRITE_LUN, 6002, IOSTAT=STATUS) 
	1           ENG_REC.EN_HEAD.LIMITS,
	1           ENG_REC.EN_HEAD.CONVERT

6002	  FORMAT (/' Limits File Name:      ', A23,
	1 	  /' Conversions File Name: ', A23 //)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	IF (ENG_REC.EN_HEAD.KEY_ID .EQ. 0) THEN

	   DO	I = 1, 4
	      CALL CT_BINARY_TO_GMT (ENG_REC.EN_HEAD.SCI_TIME(I).BIN_TIME,
	1				GMT(I))
	   ENDDO

	   WRITE (WRITE_LUN, 6010, IOSTAT=STATUS)
	1	ENG_REC.EN_HEAD.SCI_TIME(1).BIN_TIME(2),
	1	ENG_REC.EN_HEAD.SCI_TIME(1).BIN_TIME(1),
	1	(GMT(1)(I:I),I=1,14),
	1	ENG_REC.EN_HEAD.SCI_TIME(2).BIN_TIME(2),
	1	ENG_REC.EN_HEAD.SCI_TIME(2).BIN_TIME(1),
	1	(GMT(2)(I:I),I=1,14),
	1	ENG_REC.EN_HEAD.SCI_TIME(3).BIN_TIME(2),
	1	ENG_REC.EN_HEAD.SCI_TIME(3).BIN_TIME(1),
	1	(GMT(3)(I:I),I=1,14),
	1	ENG_REC.EN_HEAD.SCI_TIME(4).BIN_TIME(2),
	1	ENG_REC.EN_HEAD.SCI_TIME(4).BIN_TIME(1),
	1	(GMT(4)(I:I),I=1,14)

6010	   FORMAT (' RH binary time: ', Z8.8, 2X, Z8.8,
	1	'  (', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, ')',
	1	T65, ' RL binary time: ', Z8.8, 2X, Z8.8,
	1	'  (', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, ')'/
	1	  ' LH binary time: ', Z8.8, 2X, Z8.8,
	1	'  (', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, ')',
	1	T65, ' LL binary time: ', Z8.8, 2X, Z8.8,
	1	'  (', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, ')')

	   WRITE (WRITE_LUN, 6020, IOSTAT=STATUS)
	1	ENG_REC.EN_HEAD.DQ_SUMMARY_FLAG(1),
	1       ENG_REC.EN_HEAD.IFG_NO(1),
	1	ENG_REC.EN_HEAD.ATT_SUMMARY_FLAG(1),
	1	ENG_REC.EN_HEAD.DQ_SUMMARY_FLAG(2),
	1       ENG_REC.EN_HEAD.IFG_NO(2),
	1	ENG_REC.EN_HEAD.ATT_SUMMARY_FLAG(2),
	1	ENG_REC.EN_HEAD.DQ_SUMMARY_FLAG(3),
	1       ENG_REC.EN_HEAD.IFG_NO(3),
	1	ENG_REC.EN_HEAD.ATT_SUMMARY_FLAG(3),
	1	ENG_REC.EN_HEAD.DQ_SUMMARY_FLAG(4),
	1       ENG_REC.EN_HEAD.IFG_NO(4),
	1	ENG_REC.EN_HEAD.ATT_SUMMARY_FLAG(4)

6020	   FORMAT (/' RH instrument quality: ', I7,10x
	1           '    IFG_Number segment: ', I7,10x
	1	    '    attitude quality: ', I7, /
	1	    ' RL instrument quality: ', I7,10x
	1           '    IFG_Number segment: ', I7,10x
	1	    '    attitude quality: ', I7, /
	1	    ' LH instrument quality: ', I7,10x
	1           '    IFG_Number segment: ', I7,10x
	1	    '    attitude quality: ', I7, /
	1	    ' LL instrument quality: ', I7,10x
	1           '    IFG_Number segment: ', I7,10x
	1	    '    attitude quality: ', I7 /)

	ELSE

	   CALL CT_BINARY_TO_GMT (ENG_REC.EN_HEAD.FIRST_ENG_TIME,
	1				GMT(1))
	   CALL CT_BINARY_TO_GMT (ENG_REC.EN_HEAD.LAST_ENG_TIME,
	1				GMT(2))

	   WRITE (WRITE_LUN, 6030, IOSTAT=STATUS)
	1	ENG_REC.EN_HEAD.KEY_ID,
	1	ENG_REC.EN_HEAD.FIRST_ENG_TIME(2),
	1	ENG_REC.EN_HEAD.FIRST_ENG_TIME(1),
	1	(GMT(1)(I:I),I=1,14),
	1	ENG_REC.EN_HEAD.LAST_ENG_TIME(2),
	1	ENG_REC.EN_HEAD.LAST_ENG_TIME(1),
	1	(GMT(2)(I:I),I=1,14),
	1	ENG_REC.EN_HEAD.NUMBER_OF_RECORDS

6030	   FORMAT (' Type of Engineering Statistic: ', I4,/
	1	' Binary time of first engineering record: ', Z8.8, 2X, Z8.8,
	1	'  (', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, ')'/
	1	' Binary time of last engineering record:  ', Z8.8, 2X, Z8.8,
	1	'  (', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, ')'/
	1	' Number of engineering records used to form statistic: ', I6)

	END IF

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C     Status Monitor Command Words
C
	WRITE (WRITE_LUN, 6040, IOSTAT=STATUS)

6040	FORMAT (/' ******** Engineering Status Fields ********')

	CALL FGA_LIST_ENGSTAT (WRITE_LUN, ENG_REC.EN_STAT)

C
C GRT's
C
	WRITE (WRITE_LUN, 6050, IOSTAT=STATUS)

6050	FORMAT (/' ******** Engineering Analog Quantities ********')

	CALL FGA_LIST_ENGANLG (WRITE_LUN, ENG_REC.EN_ANALOG)

C
C The following are new fields in the new ENG format
C
	WRITE (WRITE_LUN, 6060, IOSTAT=STATUS)	ENG_REC.EN_XCAL.POS

6060	FORMAT (/' **** Channel Specific Information *****'/
	1	' XCAL position, Side A: ', I7, T65,
	1	' XCAL position, Side B: ', I7)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6070, IOSTAT=STATUS)
	1			ENG_REC.CHAN(1).UP_SCI_MODE,
	1			ENG_REC.CHAN(1).FAKEIT,
	1			ENG_REC.CHAN(1).UP_ADDS_PER_GROUP,
	1			ENG_REC.CHAN(1).UP_SWPS_PER_IFG,
	1			ENG_REC.CHAN(2).UP_SCI_MODE,
	1			ENG_REC.CHAN(2).FAKEIT,
	1			ENG_REC.CHAN(2).UP_ADDS_PER_GROUP,
	1			ENG_REC.CHAN(2).UP_SWPS_PER_IFG,
	1			ENG_REC.CHAN(3).UP_SCI_MODE,
	1			ENG_REC.CHAN(3).FAKEIT,
	1			ENG_REC.CHAN(3).UP_ADDS_PER_GROUP,
	1			ENG_REC.CHAN(3).UP_SWPS_PER_IFG,
	1			ENG_REC.CHAN(4).UP_SCI_MODE,
	1			ENG_REC.CHAN(4).FAKEIT,
	1			ENG_REC.CHAN(4).UP_ADDS_PER_GROUP,
	1			ENG_REC.CHAN(4).UP_SWPS_PER_IFG

6070	FORMAT (' Chan. ', '< Science Mode >', X,'  < Fake-it >   ',
	1	X, ' < Adds/Group > ', X, ' < Sweeps/IFG > ' /
	1	'  RH:  ', 3(I10, 8X), I10 /
	1	'  RL:  ', 3(I10, 8X), I10 /
	1	'  LH:  ', 3(I10, 8X), I10 /
	1	'  LL:  ', 3(I10, 8X), I10)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6080, IOSTAT=STATUS)
	1			ENG_REC.CHAN(1).XMIT_MTM_SPEED,
	1			ENG_REC.CHAN(1).XMIT_MTM_LEN,
	1			ENG_REC.CHAN(1).SCI_GAIN,
	1			ENG_REC.CHAN(1).DITHER,
	1			ENG_REC.CHAN(2).XMIT_MTM_SPEED,
	1			ENG_REC.CHAN(2).XMIT_MTM_LEN,
	1			ENG_REC.CHAN(2).SCI_GAIN,
	1			ENG_REC.CHAN(2).DITHER,
	1			ENG_REC.CHAN(3).XMIT_MTM_SPEED,
	1			ENG_REC.CHAN(3).XMIT_MTM_LEN,
	1			ENG_REC.CHAN(3).SCI_GAIN,
	1			ENG_REC.CHAN(3).DITHER,
	1			ENG_REC.CHAN(4).XMIT_MTM_SPEED,
	1			ENG_REC.CHAN(4).XMIT_MTM_LEN,
	1			ENG_REC.CHAN(4).SCI_GAIN,
	1			ENG_REC.CHAN(4).DITHER,
	1			ENG_REC.EN_TAIL.HSKP_FLAG

6080	FORMAT (/' Chan. ', ' < MTM Speed >  ', 2X, ' < MTM Length > '
	1	2X, '  < Sci Gain >  ', 2X, '   < Dither >   ' /
	1	'  RH:  ', 3(I10, 8X), I10 /
	1	'  RL:  ', 3(I10, 8X), I10 /
	1	'  LH:  ', 3(I10, 8X), I10 /
	1	'  LL:  ', 3(I10, 8X), I10 /
	1	' No HSK flag (0 = OK, 1 = No HSK): ', i5)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C	Telemetry Quality
	
	WRITE (WRITE_LUN, 6090, IOSTAT=STATUS)
	1                       ENG_REC.EN_TAIL.TLM_QUAL_MAJ_FRM(1),
	1                       ENG_REC.EN_TAIL.TLM_QUAL_MAJ_FRM(2)

6090	FORMAT (/'***** Overall Telemetry Quality *****' /
	1        ' Telemetry Quality Major Frame 1: ',8x, I7 /
	1        ' Telemetry Quality Major Frame 2: ',8x, I7)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	WRITE (WRITE_LUN, 7000, IOSTAT=STATUS)

7000	FORMAT (/' ****  TEMP_DIFFERENCE: Side A  ****', T65,
	1	' ****  TEMP_DIFFERENCE: Side B  **** '  )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	WRITE (WRITE_LUN, 7010, IOSTAT=STATUS)
	1		       ( ENG_REC.EN_TEMPDIFF(1).XCAL,
	1                        ENG_REC.EN_TEMPDIFF(2).XCAL,
	1			 ENG_REC.EN_TEMPDIFF(1).ICAL,
	1                        ENG_REC.EN_TEMPDIFF(2).ICAL,
	1			 ENG_REC.EN_TEMPDIFF(1).SKYHORN,
	1                        ENG_REC.EN_TEMPDIFF(2).SKYHORN,
	1			 ENG_REC.EN_TEMPDIFF(1).REFHORN,
	1                        ENG_REC.EN_TEMPDIFF(2).REFHORN,
	1			 ENG_REC.EN_TEMPDIFF(1).DIHEDRAL,
	1                        ENG_REC.EN_TEMPDIFF(2).DIHEDRAL,
	1			 ENG_REC.EN_TEMPDIFF(1).COLLIMATOR_MIRROR,
	1                        ENG_REC.EN_TEMPDIFF(2).COLLIMATOR_MIRROR,
	1			 ENG_REC.EN_TEMPDIFF(1).BOL_ASSEM(1),
	1			 ENG_REC.EN_TEMPDIFF(2).BOL_ASSEM(1),
	1			 ENG_REC.EN_TEMPDIFF(1).BOL_ASSEM(2),
	1			 ENG_REC.EN_TEMPDIFF(2).BOL_ASSEM(2),
	1			 ENG_REC.EN_TEMPDIFF(1).BOL_ASSEM(3),
	1			 ENG_REC.EN_TEMPDIFF(2).BOL_ASSEM(3),
	1			 ENG_REC.EN_TEMPDIFF(1).BOL_ASSEM(4),
	1			 ENG_REC.EN_TEMPDIFF(2).BOL_ASSEM(4))

7010	FORMAT (' External Calibrator:        ', I5, t65,
	1	' External Calibrator:          ', I5
	1	/ ' Internal Reference:         ', I5, t65,
	1	' Internal Reference:           ', I5
	1	/ ' Sky Horn:                   ', I5, t65,
	1	' Sky Horn:                     ', I5
	1	/ ' Reference Horn:             ', I5, t65,
	1	' Reference Horn:               ', I5
	1	/ ' Left Dihedral:              ', I5, t65,
	1	' Right Dihedral:               ', I5
	1	/ ' Collimator Mirror:          ', I5, t65,
	1	' Collimator Mirror:            ', I5
	1	/ ' Bolometer Assembly, RH:       ',I7, t65,
	1	' Bolometer Assembly, RH:       ', I7
	1	/ ' Bolometer Assembly, RL:       ',I7, t65,
	1	' Bolometer Assembly, RL:       ', I7
	1	/ ' Bolometer Assembly, LH:       ',I7, t65,
	1	' Bolometer Assembly, LH:       ', I7
	1	/ ' Bolometer Assembly, LL:       ',I7, t65,
	1	' Bolometer Assembly, LL:       ', I7)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF
	
	RETURN
	END
