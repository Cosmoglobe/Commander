	Integer*4 Function FGA_LIST_IDXTOLS (WRITE_LUN, IDXTOLS, FILE_NAME,
	1				     CLINE, NREC)

C---------------------------------------------------------------------------
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given engineering tolerances
C	  record.
C
C	AUTHOR:
C	  QUOC CHUNG
C	  STX
C	  NOV 08, 1989
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_IDXTOLS (WRITE_LUN, IDXTOLS, FILE_NAME, CLINE, NREC)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file
C	  IDXTOLS         		FEX_IDX_tols record.
C	  NREC			I*2	The number of records processed in this run.
C	  CLINE			CH*80	A user supplied comment line.
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
    
	CHARACTER	*14	GMT

	INTEGER		*2	I
	INTEGER		*2	J
	INTEGER		*2	K
	INTEGER		*2	L

	INTEGER		*2	NREC
	INTEGER		*4	WRITE_LUN
	INTEGER		*4	IOSTATUS
        INTEGER         *4      STATUS

	CHARACTER       *64     FILE_NAME		
	CHARACTER       *80     CLINE
	CHARACTER	*1	LINE(130)/ 130 * '-'/

	LOGICAL         *1      NORMAL/.TRUE./
        LOGICAL         *1      FIRST/.TRUE./

	EXTERNAL		FGA_WRITE_ERR

	DICTIONARY	'FEX_IDX_TOLS'
	RECORD /FEX_IDX_TOLS/ IDXTOLS

C
C Begin the dump.
C
	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6
	IF ( FIRST) THEN
          L = 0
	  WRITE (WRITE_LUN, 600, IOSTAT=IOSTATUS) FILE_NAME,CLINE,LINE
          FIRST = .FALSE.
        ENDIF
 600    FORMAT(	// 1X, 'FEX_IDX_TOLS File : ',A/,X,A80/,130A/)

	IF (IOSTATUS .NE. 0) THEN
	  PRINT *,' WRITE ERROR!!! UNIT = ', WRiTE_LUN 
	END IF

	WRITE ( WRITE_LUN,601) NREC
 601    FORMAT(/,T50,' REC # : ',I4/)
C
C Cobetrieve standard header.
C
	WRITE (WRITE_LUN, 6000, IOSTAT=STATUS)
	1			(IDXTOLS.CT_HEAD.GMT(I:I), I=1,14),
	1			IDXTOLS.CT_HEAD.TIME(2),
	1			IDXTOLS.CT_HEAD.TIME(1),
	1			IDXTOLS.CT_HEAD.SPACE_TIME,
	1			IDXTOLS.CT_HEAD.MJR_FRM_NO,
	1			IDXTOLS.CT_HEAD.ORBIT

6000	FORMAT (/' Processed time (in 3 formats): ' /
	1	' <GMT> : ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	' <Binary> : ', Z8.8, 1X, Z8.8, 7X,
	1	' <PB5> : ',  6(Z2.2, 1X) /
	1	' Major Frame Number : ', I6, 8x, 
	1	' Orbit number : ', I6)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF
        
	WRITE (WRITE_LUN, 6001, IOSTAT=IOSTATUS) 
	1      IDXTOLS.IDX_TOLS.TOLS,LINE
 6001   FORMAT(/,' Index Tolerances For Engineering Fields : '/,40x,
	1	'  01 -- 10  : ',10i6/,40x,
	1	'  11 -- 20  : ',10i6/,40x,
	1	'  21 -- 30  : ',10i6/,40x,
	1	'  31 -- 40  : ',10i6/,40x,
	1	'  41 -- 50  : ',10i6/,40x,
	1	'  51 -- 60  : ',10i6/,40x,
	1	'  61 -- 70  : ',10i6/,40x,
	1	'  71 -- 80  : ',10i6/,40x,
	1	'  81 -- 90  : ',10i6/,40x,
	1	'  91 -- 100 : ',10i6/,40x,
	1	' 101 -- 110 : ',10i6/,40x,
	1	' 111 -- 120 : ',10i6/,40x,
	1	' 121 -- 130 : ',10i6/,40x,
	1	' 131 -- 140 : ',10i6/,40x,
	1	' 141 -- 150 : ',10i6/,40x,
	1	' 151 -- 160 : ',10i6/,40x,
	1	' 161 -- 170 : ',10i6/,40x,
	1	' 171 -- 180 : ',10i6/,40x,
	1	' 181 -- 190 : ',10i6/,40x,
	1	' 191 -- 200 : ',10i6/,40x,
	1	' 201 -- 210 : ',10i6/,40x,
	1	' 211 -- 220 : ',10i6/,40x,
	1	' 221 -- 230 : ',10i6/,40x,
	1	' 231 -- 240 : ',10i6/,40x,
	1	' 241 -- 250 : ',10i6/,40x,
	1	' 251 -- 256 : ',6i6//,130a)

        IF (IOSTATUS .NE. 0 ) THEN
          NORMAL = .FALSE.
          PRINT *, ' Formatted write error ', WRITE_LUN
        END IF

        RETURN
        END
