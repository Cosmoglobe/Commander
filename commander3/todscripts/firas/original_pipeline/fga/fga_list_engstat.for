	SUBROUTINE FGA_LIST_ENGSTAT (WRITE_LUN, STAT_REC)

C------------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FGA_LIST_ENGSTAT
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
C		Adapted from work done by Ed Fung.
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_ENGSTAT (WRITE_LUN, STAT_REC)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file
C
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
C      Changes:
C       Add DWELL_STAT fields (SPR 4234), N. Gonzales, STX, 1990 July 10.
C       ADD HOT_SPOT_CMD fields (SPR 8372), N. Gonzales,STX, 1991 May 8.
C------------------------------------------------------------------------

	IMPLICIT	NONE

	INTEGER		*4	WRITE_LUN
	INTEGER		*4	STATUS
	CHARACTER	*7	REC_STAT
	CHARACTER	*14	GMT(4)
	INTEGER		*2	I
	INTEGER		*2	J

	DICTIONARY	'FUT_ENGSTAT'
	RECORD /FUT_ENGSTAT/ STAT_REC

	EXTERNAL		FGA_WRITE_ERR

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6

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
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Other Status Words
C
	WRITE (WRITE_LUN, 6004, IOSTAT=STATUS)
	1			STAT_REC.GRT_ADDR(1),
	1			STAT_REC.GRT_ADDR(2),
	1			STAT_REC.MICRO_STAT_BUS(1),
	1			STAT_REC.MICRO_STAT_BUS(2),
	1			STAT_REC.MICRO_STAT_BUS(3),
	1			STAT_REC.MICRO_STAT_BUS(4),
	1			STAT_REC.BOL_CMD_BIAS(1),
	1			STAT_REC.BOL_CMD_BIAS(2),
	1			STAT_REC.BOL_CMD_BIAS(3),
	1			STAT_REC.BOL_CMD_BIAS(4),
	1                       STAT_REC.DWELL_STAT(1),
	1                       STAT_REC.DWELL_STAT(2),
	1			STAT_REC.LVDT_STAT(1),
	1			STAT_REC.LVDT_STAT(2),
	1			STAT_REC.POWER_A_STATUS(1),
	1			STAT_REC.POWER_A_STATUS(2),
	1			STAT_REC.POWER_B_STATUS(1),
	1			STAT_REC.POWER_B_STATUS(2),
	1                       STAT_REC.HOT_SPOT_CMD(1),
	1                       STAT_REC.HOT_SPOT_CMD(2)	
6004	FORMAT (' ***** OTHER STATUSES ***** '
	1	/' Dwell GRT Address, Side A (HEX):      ',2x, Z2.2,
	1	T65,' Dwell GRT Address, Side B (HEX):      ',2x, Z2.2,
	1	/' Microprocessor bus readout, RH (HEX): ', Z2.2,
	1	T65,' Microprocessor bus readout, RL (HEX): ', Z2.2,
	1	/' Microprocessor bus readout, LH (HEX): ', Z2.2,
	1	T65,' Microprocessor bus readout, LL (HEX): ', Z2.2,
	1	/' Bolometer bias status, RH (HEX):    ', 2X, Z2.2
	1	T65,' Bolometer bias status, RL (HEX):    ', 2X, Z2.2
	1	/' Bolometer bias status, LH (HEX):    ', 2X, Z2.2
	1	T65,' Bolometer bias status, LL (HEX):    ', 2X, Z2.2
	1	/' Dwell Status, Side A (HEX):         ',2x, Z2.2,
	1	T65,' Dwell Status, Side B (HEX):         ',2x, Z2.2,
	1	/' LVDT Status, Side A (HEX): ', Z2.2,
	1	T65,' LVDT Status, Side B (HEX):    ', 2X, Z2.2,
	1	/' Power Status, Side A, Frm 1 (HEX): ', Z2.2,
	1	T65,' Power Status, Side A, Frm 2 (HEX):    ', 2X, Z2.2,
	1	/' Power Status, Side B, Frm 1 (HEX): ', Z2.2,
	1	T65,' Power Status, Side B, Frm 2 (HEX):    ', 2X, Z2.2,
	1	/' Hot Spot Command, Side A, (HEX):   ', Z2.2,
	1	T65,' Hot Spot Command, Side B, (HEX):      ', 2X, Z2.2)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	RETURN
	END
