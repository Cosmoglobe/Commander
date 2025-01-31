	Integer*4 Function FGA_LIST_MTMSWEEP (WRITE_LUN, MTMSWEEP_REC,
	1                  MTMSWEEP_FILE, CLINE, NREC)

C---------------------------------------------------------------------------
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given FIRAS 
C	  reference data of MTMSWEEP.
C
C	AUTHOR:
C	  Nilo G. Gonzales
C	  STX
C	  September 23, 1991   Reference: SPR 8372
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_MTMSWEEP FGA_LIST_MTMSWEEP (WRITE_LUN, MTMSWEEP_REC,
C                                MTMSWEEP_FILE, CLINE, NREC)
C	INPUT PARAMETERS:
C	  WRITE_LUN	 I*4    Logical unit number for a list file;
C				if not supplied assume 6.  Assume the
C				calling routine has opened this list file.
C	  MTMSWEEP_REC         	MTMSWEEP record.
C	  CLINE		 CH*80  A user supplied comment line.
C	  NREC		 I*2    The number of records processed in this run.
C                                      
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
C	  CT_BINARY_TO_GMT
C
C	ERROR HANDLING:
C	  NONE
C
C---------------------------------------------------------------------------
C Changes:
C---------------------------------------------------------------------------

	IMPLICIT	NONE

	INTEGER		*2	NREC
        INTEGER         *4      SUCCESS/1/, ERROR/2/
        LOGICAL         *1	FIRST/.TRUE./
	INTEGER		*4	WRITE_LUN
	INTEGER         *2      I
	INTEGER		*4	STATUS, IOSTAT
	REAL            *4      SS, LS, SF, LF
	REAL            *4      FSS, FLS, FSF, FLF
	CHARACTER       *64     MTMSWEEP_FILE
	CHARACTER	*14	GMT
	CHARACTER       *80	CLINE
	CHARACTER	*1	LINE(130)/ 130 * '-'/
	EXTERNAL		FGA_WRITE_ERR

	DICTIONARY	'FEX_MTMSWEEP'
	RECORD /FEX_MTMSWEEP/ MTMSWEEP_REC

        FGA_LIST_MTMSWEEP = SUCCESS

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6
	IF ( FIRST) THEN
	  WRITE (WRITE_LUN, 10, IOSTAT=STATUS) MTMSWEEP_FILE,CLINE,LINE
          FIRST = .FALSE.
        ENDIF
 10   FORMAT(	// 1X, 'FEX_MtmSweep FILE : ',A/,X,A80/,130A/)

	IF (IOSTAT .NE. 0) THEN
	  PRINT *,' WRITE ERROR!!! UNIT = ', WRITE_LUN 
          FGA_LIST_MTMSWEEP = ERROR          
	END IF

	WRITE ( WRITE_LUN,20) NREC
  20   FORMAT(/,T30,' REC # : ',I4)
C
C Cobetrieve standard header.
C
	WRITE (WRITE_LUN, 70, IOSTAT=STATUS)
	1     (MTMSWEEP_REC.CT_HEAD.GMT(I:I), I=1,14),
	1     MTMSWEEP_REC.CT_HEAD.TIME(2),
	1     MTMSWEEP_REC.CT_HEAD.TIME(1),
	1     MTMSWEEP_REC.CT_HEAD.SPACE_TIME,
	1     MTMSWEEP_REC.CT_HEAD.MJR_FRM_NO,
	1     MTMSWEEP_REC.CT_HEAD.ORBIT,
	1     MTMSWEEP_REC.CT_HEAD.HSKP1_TLM_FMT,
	1     MTMSWEEP_REC.CT_HEAD.HSKP2_TLM_FMT

70	FORMAT (/' Processed time (in 3 formats): ' /
	1	' <GMT> : ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	' <Binary> : ', Z8.8, 1X, Z8.8, 7X,
	1	' <PB5> : ',  6(Z2.2, 1X) /
	1	' Major Frame Number : ', T30, I6 /
	1	' Orbit number : ', T30, I6 /
	1       ' TLM FMT Major Frame 1: ',10X,I4 /
	1       ' TLM FMT Major Frame 2: ',10X,I4) 

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	WRITE (WRITE_LUN, 75, IOSTAT=STATUS)
75	  FORMAT (//' *****  MTM Sweeps in unit seconds  *****'/,
	1           ' Total time = Turn Around + Sweep + Flyback'/)    

C       Convert Total time to unit seconds

	SS = FLOAT (MTMSWEEP_REC.TOTAL_SWEEP_FLYBACK(1)) / 1.E7
	LS = FLOAT (MTMSWEEP_REC.TOTAL_SWEEP_FLYBACK(2)) / 1.E7
	SF = FLOAT (MTMSWEEP_REC.TOTAL_SWEEP_FLYBACK(3)) / 1.E7
	LF = FLOAT (MTMSWEEP_REC.TOTAL_SWEEP_FLYBACK(4)) / 1.E7
C
C Dump MTM sweeps record.
C
	WRITE (WRITE_LUN, 80, IOSTAT=STATUS) SS, LS, SF, LF

80	FORMAT (10X,'Short Slow: ',F7.3/
	1      10X,'Long Slow: ',F7.3/
	1      10X,'Short Fast: ',F7.3/
	1      10X,'Long Fast: ',F7.3)

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	WRITE (WRITE_LUN, 85, IOSTAT=STATUS)
85	FORMAT (//' *****  MTM FLYBACK in unit seconds  *****'/)

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C       Convert Flyback time to unit seconds

	FSS = FLOAT (MTMSWEEP_REC.FLYBACK(1)) / 1.E7
	FLS = FLOAT (MTMSWEEP_REC.FLYBACK(2)) / 1.E7
	FSF = FLOAT (MTMSWEEP_REC.FLYBACK(3)) / 1.E7
	FLF = FLOAT (MTMSWEEP_REC.FLYBACK(4)) / 1.E7

	WRITE (WRITE_LUN, 90, IOSTAT=STATUS) FSS, FLS, FSF, FLF 

90      FORMAT (10X,'Short Slow: ',F7.3/
	1     10X,'Long Slow: ',F7.3/
	1     10X,'Short Fast: 'F7.3/
	1     10X,'Long Fast: ',F7.3)

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	RETURN
	END
