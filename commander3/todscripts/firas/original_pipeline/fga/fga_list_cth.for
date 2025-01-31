	INTEGER*4 FUNCTION FGA_LIST_CTH (WRITE_LUN, CTH_REC, CTH_FILE, CLINE,
	1				 NREC)

C---------------------------------------------------------------------------
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given FEX_CTH
C	  record.
C
C	AUTHOR:
C	  QUOC CHUNG
C	  STX
C	  Nov 15, 89 SPR 4750
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_CTH (WRITE_LUN, CTH_REC, CTH_FILE, NREC)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file
C	  CTH_REC         	RECORD	FEC_CTH record.
C	  CTH_FILE		CH*64	FEX_CTH file specification.
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
C
C---------------------------------------------------------------------------

	IMPLICIT	NONE

	INTEGER		*2	NREC
        INTEGER         *2      HEAD3
	INTEGER		*4	WRITE_LUN
	INTEGER		*4	I
	INTEGER		*4	J
	INTEGER		*4	K
	INTEGER		*4	L
	INTEGER		*4	M
	INTEGER		*4	IOSTATUS
        INTEGER         *4      STATUS
        INTEGER         *4      SUCCESS/1/, ERROR/2/
        LOGICAL         *1	FIRST/.TRUE./

	CHARACTER       *64     CTH_FILE		
	CHARACTER       *80	CLINE
	CHARACTER	*1	LINE(130)/ 130 * '-'/
	CHARACTER       *30 	TEXT(12)/'MIN_IFG_COADD','MIN_IFG_SIGMA',
	1                   		 'MAX_IFG_SIGMA','MAX_POINT_DEVIATION',
	2                   		 'MAX_BAD_POINTS','MAX_SATURATED_POINTS',
	3                   		 'MAX_GLITCHES','MAX_DATA_DROPOUT',
	4                   		 'BOLOMETER_VOLTAGE_TOLERANCES',
	5                   		 'GRT_TOLERANCES','RESERVE',
	6				 'TIME_TOLERANCE'/

	INTEGER		*2	TIME_LEN
	CHARACTER	*32	TIME
	INTEGER		*4	SYS$ASCTIM

	DICTIONARY		'FEX_CTH'
	RECORD /FEX_CTH/	CTH_REC

C
C Begin the dump.
C
        FGA_LIST_CTH = SUCCESS

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6
	IF ( FIRST) THEN
	  WRITE (WRITE_LUN, 6000, IOSTAT=IOSTATUS) CTH_FILE,CLINE,LINE
          FIRST = .FALSE.
        ENDIF
 6000   FORMAT(	// 1X, 'FEX_CTH FILE : ',A/,X,A80/,130A/)

	IF (IOSTATUS .NE. 0) THEN
	  PRINT *,' WRITE ERROR!!! UNIT = ', WRiTE_LUN 
          FGA_LIST_CTH = error          
	END IF

	WRITE ( WRITE_LUN,600) NREC
  600   FORMAT(/,T50,' REC # : ',I4/)
C	
C Header information.
C
	WRITE (WRITE_LUN, 6001, IOSTAT=IOSTATUS)
	1              (CTH_REC.CT_HEAD.GMT(I:I),I=1,14),
	1		CTH_REC.CT_HEAD.TIME(2),
	1		CTH_REC.CT_HEAD.TIME(1),
	1		CTH_REC.CT_HEAD.SPACE_TIME,
	1		CTH_REC.CT_HEAD.MJR_FRM_NO,
	1		CTH_REC.CT_HEAD.ORBIT,
	1		CTH_REC.CT_HEAD.DATASET_ID,
	1               LINE

6001	FORMAT (/' <Start GMT> : ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	 ' <Binary TIME> : ', Z8.8, 1X, Z8.8, 7X, /,
	1        ' <PB5> : ', 6(Z2.2, 1X) /
	1        ' Major Frame Number : ', I6, 8x,
	1        ' Orbit Number : ', I12, 12x,
	1	 ' Data Set Id  : ', I10/,130A/)

	STATUS = SYS$ASCTIM ( TIME_LEN, TIME, CTH_REC.TIME_TOLERANCE, 0 )

	WRITE (WRITE_LUN, 6002, IOSTAT=IOSTATUS) 
	1        TEXT(1),(CTH_REC.MIN_IFG_COADD(I),I=1,4),
	1        TEXT(2),(CTH_REC.MIN_IFG_SIGMA(I), I=1,4),
	1        TEXT(3),(CTH_REC.MAX_IFG_SIGMA(I), I=1,4),
	1        TEXT(4),(CTH_REC.MAX_POINT_DEVIATION(I),I=1,4),
	1        TEXT(5),(CTH_REC.MAX_BAD_POINTS(I),I=1,4),
	1        TEXT(6),(CTH_REC.MAX_SATURATED_POINTS(I),I=1,4),
	1        TEXT(7),(CTH_REC.MAX_GLITCHES(I),I=1,4),
	1        TEXT(8),(CTH_REC.MAX_DATA_DROPOUT(I),I=1,4),
	1        TEXT(9),(CTH_REC.BOLOMETER_VOLTAGE_TOLERANCES(I),I=1,4),
	1        TEXT(10),(CTH_REC.GRT_TOLERANCES(I),I=1,16),
	1        TEXT(11),(CTH_REC.RESERVE(I),I=1,4),
	1        TEXT(12),CTH_REC.TIME_TOLERANCE(2),
	1        CTH_REC.TIME_TOLERANCE(1),TIME,
	1               LINE

6002    FORMAT  (1X,'FEX_CTH Data Record : '/,
	1         1X,A,' : ',4I12/,
	1         1X,A,' : ',4F12.2/,
	1         1X,A,' : ',4F12.2/,
	1         1X,A,' : ',4F12.2/,
	1         1X,A,' : ',4I12/,
	1         1X,A,' : ',4I12/,
	1         1X,A,' : ',4I12/,
	1         1X,A,' : ',4I12/,
	1         1X,A,' : ',4F12.2/,
	1         1X,A,' : ',8F10.4/,34X,8F10.4/,
	1         1X,A,' : ',4I12/,
	1         1X,A,' : ',Z8.8,1X,Z8.8,5X,A32//,130A/)

        RETURN
	END
