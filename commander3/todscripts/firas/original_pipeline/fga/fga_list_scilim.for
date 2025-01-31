	Integer*4 Function FGA_LIST_SCILIM (WRITE_LUN, SCILIM, SCILIM_FILE,
	1				    CLINE, NREC)

C---------------------------------------------------------------------------
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given raw science
C	  limits record.
C
C	AUTHOR:
C	  QUOC CHUNG
C	  STX
C	  OCT 24, 1989
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_SCILIM (WRITE_LUN, SCILIM, SCILIM_FILE, CLINE, NREC)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file
C	  SCILIM         		SCILIM record.
C	  CLINE			CH*80	A user supplied comment line.
C	  NREC			I*2	The number of records processed in this run.
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
C
C	SPR 6159.  May 25, 1990, R. Kummerer, STX@GSFC, Version 6.1.
C		   Clean up reference data set lists.
C       SPR 7422.  Sept. 21, 1990, N. Gonzales,STX
C                  Limits are switched from Yellow to Red and Red to Yellow. 
C       SPR 8372.  Sept. 30, 1991, N. Gonzales/STX
C                  Add glitch rate field.
C	SPR 9162.  Oct. 15, 1991, S. Alexander/STX
C		   Remove misleading total glitch count output line.
C---------------------------------------------------------------------------

	IMPLICIT	NONE
        
	INCLUDE    '(FUT_PARAMS)'
    
	CHARACTER	*14	GMT

	INTEGER		*4	SW_VERSION_EQV
	BYTE			SW_VERSION(2)
	INTEGER		*4	CHAN_IND_EQV
	BYTE			CHAN_IND(2)

	INTEGER		*2	SC5
	INTEGER		*2	SC13
	INTEGER		*2	I
	INTEGER		*2	J
	INTEGER		*2	K
	INTEGER		*2	L

	INTEGER		*4	SUM
	INTEGER		*4	UNSIGNED_CHECKSUM
	INTEGER		*4	MIN_FRM_CNT
	INTEGER		*4	DATA_TRANS_TIME

	CHARACTER       *14     COLLECT_GMT

	INTEGER		*2	NREC
        INTEGER         *2      HEAD3
	INTEGER		*4	WRITE_LUN
	INTEGER		*4	IOSTATUS
        INTEGER         *4      STATUS
        INTEGER         *4      MJR_FRAME
        INTEGER         *4      SUCCESS/1/,ERROR/2/


	CHARACTER       *64     SCILIM_FILE		
	CHARACTER	*1	LINE(130)/ 130 * '-'/
	CHARACTER	*80	CLINE
        CHARACTER       *1      CBITS16(16)
        CHARACTER       *1      CBITS32(32)
	CHARACTER       *14     LTEXT(2)/'RED','YELLOW'/
	INTEGER		*4	LLEN(2)/3,6/

	REAL		*4	ELIMB
	REAL		*4	SUNANG
	REAL		*4	MOONANG
	CHARACTER	*29	PDEF
	CHARACTER	*25	INTRO
	INTEGER		*4	INDEX
	REAL		*4	ORIGIN
	
        LOGICAL         *1      FIRST/.TRUE./
	EQUIVALENCE	(SW_VERSION_EQV, SW_VERSION)
	EQUIVALENCE	(CHAN_IND_EQV, CHAN_IND)
	EXTERNAL		FGA_WRITE_ERR

	DICTIONARY	'FEX_SCILIM'
	RECORD /FEX_SCILIM/      SCILIM

C
C Begin the dump.
C
	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6
	IF ( FIRST) THEN
          L = 0
	  WRITE (WRITE_LUN, 600, IOSTAT=IOSTATUS) SCILIM_FILE,CLINE,LINE
          FIRST = .FALSE.
        ENDIF
 600    FORMAT(	// 1X, 'FEX_SCILIM FILE : ',A/,X,A80/,130A/)

	IF (IOSTATUS .NE. 0) THEN
	  PRINT *,' WRITE ERROR!!! UNIT = ', WRITE_LUN 
	END IF

	WRITE ( WRITE_LUN,601) NREC
 601    FORMAT(/,T50,' REC # : ',I4/)
C
C Cobetrieve standard header.
C
	WRITE (WRITE_LUN, 6000, IOSTAT=STATUS)
	1			(SCILIM.CT_HEAD.GMT(I:I), I=1,14),
	1			SCILIM.CT_HEAD.TIME(2),
	1			SCILIM.CT_HEAD.TIME(1),
	1			SCILIM.CT_HEAD.SPACE_TIME,
	1			SCILIM.CT_HEAD.MJR_FRM_NO,
	1			SCILIM.CT_HEAD.ORBIT,
	1                       SCILIM.CT_HEAD.HSKP1_TLM_FMT,
	1                       SCILIM.CT_HEAD.HSKP2_TLM_FMT

6000	FORMAT (/' Processed time (in 3 formats): ' /
	1	' <GMT> : ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	' <Binary> : ', Z8.8, 1X, Z8.8, 7X,
	1	' <PB5> : ',  6(Z2.2, 1X) /
	1	' Major Frame Number : ', I6, 8x, 
	1	' Orbit number : ', I6 /
	1       ' TLM FMT Major Frame 1: ',10X,I4 /
	1       ' TLM FMT Major Frame 2: ',10X,I4 ) 

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	DO I = 1 , 2
C
C Science LIMIT record header.
C
	  WRITE (WRITE_LUN, 6001, IOSTAT=STATUS) LTEXT(I)(1:LLEN(I)),
	1			SCILIM.LIM(I).SCI_HEAD.DATA_READY(1),
	1			SCILIM.LIM(I).SCI_HEAD.SC_HEAD20

6001	  FORMAT (//' ***** Science Header Information For ',A,' Limits *****'/,
	1	' Data Ready: ', T30, I12 /
	1	' Saturated Sample Count: ', T30, I12)

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF
C
C Dump the attitude limits.
C Convert the attitude to recognizable units and display.
C
	  ELIMB = FLOAT(SCILIM.LIM(I).ATTITUDE.EARTH_LIMB) * FAC_ATT_CONV
	  SUNANG = FLOAT(SCILIM.LIM(I).ATTITUDE.SUN_ANGLE) * FAC_ATT_CONV
	  MOONANG = FLOAT(SCILIM.LIM(I).ATTITUDE.MOON_ANGLE) * FAC_ATT_CONV

	  WRITE (WRITE_LUN, 6020, IOSTAT=STATUS) LTEXT(I)(1:LLEN(I)),
	1		ELIMB,
	1		SUNANG,
	1		MOONANG

6020	  FORMAT (//' ***** Attitude Information For ',A,' Limits *****'/,
	1	' Earth Limb Angle (degrees): ', T50, G13.4 /
	1	' Sun Angle (degrees): ', T50, G13.4 /
	1	' Moon Angle (degrees): ', T50, G13.4)

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

	  WRITE (WRITE_LUN, 6023, IOSTAT=STATUS) LTEXT(I)(1:LLEN(I)),
	1			SCILIM.LIM(I).GLITCH_RATE(1),
	1			SCILIM.LIM(I).GLITCH_RATE(2),
	1			SCILIM.LIM(I).GLITCH_RATE(3),
	1			SCILIM.LIM(I).GLITCH_RATE(4)

6023	  FORMAT (//' ***** Glitch Rate Information For ',A,' Limits *****'/,
	1	' Channel RH: ', T30, G13.4 /
	1       '         RL: ', T30, G13.4 /
	1       '         LH: ', T30, G13.4 /
	1       '         LL: ', T30, G13.4)

	  IF (STATUS .NE. 0) THEN
	     CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF
	END DO

	WRITE (WRITE_LUN, 6025, IOSTAT=STATUS) LINE
6025    FORMAT (/,130A/)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	RETURN
	END
