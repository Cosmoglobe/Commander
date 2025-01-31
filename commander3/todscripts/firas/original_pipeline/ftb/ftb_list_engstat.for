	SUBROUTINE FTB_LIST_ENGSTAT (WRITE_LUN,ENG_REC, STAT_REC)

C------------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FTB_LIST_ENGSTAT
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given FIRAS
C	  Engineering status fields.
C
C	AUTHOR:
C	  R. Kummerer
C	  STX
C	  April 7, 1988
C      
C         H. Wang,  Modified it for FTB_CHECK_DQ to use
C         STX
C         OCT 10, 1989
C
C		Adapted from work done by Ed Fung.
C
C	CALLING SEQUENCE:
C	  CALL FTB_LIST_ENGSTAT (WRITE_LUN,ENG_REC, STAT_REC)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file
C         ENG_REC
C	  STAT_REC		RECORD	FIRAS Engineering record.
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

	IMPLICIT	NONE

	INTEGER		*4	WRITE_LUN
	INTEGER		*4	STATUS
	CHARACTER	*7	REC_STAT
	CHARACTER	*14	GMT(4)
	INTEGER		*2	I
	INTEGER		*2	J

	DICTIONARY	'FUT_ENGSTAT'
	DICTIONARY	'FDQ_ENG'
	RECORD /FUT_ENGSTAT/ STAT_REC
        record/FDQ_ENG/ ENG_REC

	EXTERNAL		FTB_WRITE

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6

C
	WRITE (WRITE_LUN, 6001, IOSTAT=STATUS)
	1			ENG_REC.CT_HEAD.GMT(1:2),
	1			ENG_REC.CT_HEAD.GMT(3:5),
	1			ENG_REC.CT_HEAD.GMT(6:7),
	1			ENG_REC.CT_HEAD.GMT(8:9),
	1			ENG_REC.CT_HEAD.GMT(10:11),
	1			ENG_REC.CT_HEAD.GMT(12:14)

6001	FORMAT (/' ASCII GMT:                             ', 4X,
	1	A2, '-', A3, '-', 3(A2, '-'), A3)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FTB_WRITE, %VAL(1), %VAL(STATUS) )
	END IF

C
C      write channel 4 science time 
C



	IF (ENG_REC.EN_HEAD.KEY_ID .EQ. 0) THEN

	   DO	I = 1, 4
	      CALL CT_BINARY_TO_GMT (ENG_REC.EN_HEAD.SCI_TIME(I).BIN_TIME,
	1				GMT(I))
	   ENDDO

	   WRITE (WRITE_LUN, 6010, IOSTAT=STATUS)
	1	(GMT(1)(I:I),I=1,14),
	1	(GMT(2)(I:I),I=1,14),
	1	(GMT(3)(I:I),I=1,14),
	1	(GMT(4)(I:I),I=1,14)
        ENDIF
6010	   FORMAT (' RH GMT time: ', 
	1	'  (', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, ')',
	1	T65, ' RL GMT time: ', 
	1	'  (', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, ')'/
	1	  ' LH GMT time: ', 
	1	'  (', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, ')',
	1	T65, ' LL GMT time: ',
	1	'  (', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, ')')

C
C     Status Monitor Command Words
C
	WRITE (WRITE_LUN, 6003, IOSTAT=STATUS)
	1			 (STAT_REC.GROUP1(I),
	1			 STAT_REC.GROUP1(I+8), I = 1, 8)

6003	FORMAT (/' <<<<< STATUS MONITOR COMMAND WORDS >>>>> '
	1	/'   -----   SIDE A   -----', T65,
	1	'   -----   SIDE B   -----' /
	1	' Status Word 1 (HEX):                  ', Z4.4, T65,
	1	' Status Word 9 (HEX):                  ', Z4.4/
	1	' Status Word 2, Int. Ref. TC 1:      ', I6, T65,
	1	' Status Word 10, Int. Ref. TC 5:     ', I6/
	1	' Status Word 3, Ref. Horn TC 2:      ', I6, T65,
	1	' Status Word 11, Ref. Horn TC 6:     ', I6/
	1	' Status Word 4 (HEX):                  ', Z4.4, T65,
	1	' Status Word 12 (HEX):                 ', Z4.4 /
	1	' Status Word 5 (HEX):                  ', Z4.4, T65,
	1	' Status Word 13 (HEX):                 ', Z4.4 /
	1	' Status Word 6, Sky Horn TC 3:       ', I6, T65,
	1	' Status Word 14, Sky Horn TC 7:      ', I6/
	1	' Status Word 7, Ext. Cal. TC 4:      ', I6, T65,
	1	' Status Word 15, Ext. Cal. TC 8:     ', I6/
	1	' Status Word 8 (HEX):                  ', Z4.4, T65,
	1	' Status Word 16 (HEX):                 ', Z4.4 )

	IF (STATUS .NE. 0) THEN
          type *, ' ftb_list_engstat write error '
	  CALL LIB$SIGNAL ( Ftb_WRITE, %VAL(1), %VAL(STATUS) )
	END IF

C
C Other Status Words
C
	WRITE (WRITE_LUN, 6004, IOSTAT=STATUS)
	1			STAT_REC.GRT_ADDR(1),
	1			STAT_REC.GRT_ADDR(2),
	1			STAT_REC.BOL_CMD_BIAS(1),
	1			STAT_REC.BOL_CMD_BIAS(2),
	1			STAT_REC.BOL_CMD_BIAS(3),
	1			STAT_REC.BOL_CMD_BIAS(4),
	1			STAT_REC.MICRO_STAT_BUS(1),
	1			STAT_REC.MICRO_STAT_BUS(2),
	1			STAT_REC.MICRO_STAT_BUS(3),
	1			STAT_REC.MICRO_STAT_BUS(4),
	1			STAT_REC.POWER_A_STATUS(1),
	1			STAT_REC.POWER_A_STATUS(2),
	1			STAT_REC.POWER_B_STATUS(1),
	1			STAT_REC.POWER_B_STATUS(2),
	1			STAT_REC.LVDT_STAT(1),
	1			STAT_REC.LVDT_STAT(2)

6004	FORMAT (' ***** OTHER STATUSES ***** '
	1	/' Dwell GRT Address, Side A (HEX):      ',2x, Z2.2,
	1	T65,' Dwell GRT Address, Side B (HEX):      ',2x, Z2.2,
	1	/' Bolometer bias status, RH (HEX):    ', 2X, Z2.2
	1	T65,' Bolometer bias status, RL (HEX):    ', 2X, Z2.2
	1	/' Bolometer bias status, LH (HEX):    ', 2X, Z2.2
	1	T65,' Bolometer bias status, LL (HEX):    ', 2X, Z2.2
	1	/' Microprocessor bus readout, RH (HEX): ', Z2.2,
	1	T65,' Microprocessor bus readout, RL (HEX): ', Z2.2,
	1	/' Microprocessor bus readout, LH (HEX): ', Z2.2,
	1	T65,' Microprocessor bus readout, LL (HEX): ', Z2.2,
	1	/' Power Status, Side A, Frm 1 (HEX): ', Z2.2,
	1	T65,' Power Status, Side A, Frm 2 (HEX):    ', 2X, Z2.2,
	1	/' Power Status, Side B, Frm 1 (HEX): ', Z2.2,
	1	T65,' Power Status, Side B, Frm 2 (HEX):    ', 2X, Z2.2,
	1	/' LVDT Status, Side A (HEX): ', Z2.2,
	1	T65,' LVDT Status, Side B (HEX):    ', 2X, Z2.2)

	IF (STATUS .NE. 0) THEN
          type *, ' ftb_list_Engsata write error '
	  CALL LIB$SIGNAL ( FTB_WRITE, %VAL(1), %VAL(STATUS) )
	END IF


	RETURN
	END
