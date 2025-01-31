	Integer*4 Function FGA_LIST_IDXFLAGS (WRITE_LUN, IDXFLAGS, FILE_NAME, CLINE, NREC)

C---------------------------------------------------------------------------
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given engineering index
C	  flags record.
C
C	AUTHOR:
C	  QUOC CHUNG
C	  STX
C	  NOV 08, 1989
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_IDXFLAGS (WRITE_LUN, IFXFLAGS, FILE_NAME, CLINE, NREC)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file
C	  IDXFLAGS         		FEX_IDX_FLAG record.
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
C       SPR 8372.  May 14, l991, N. Gonzales, STX
C                  Increase ITEM_TEXT array from 274 to 284.       
C---------------------------------------------------------------------------

	IMPLICIT	NONE
        
        INCLUDE    '(FUT_TEXT_ITEM)'
    
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
	CHARACTER       *30     ITEM_TEXT(284)
	CHARACTER       *80     CLINE
	CHARACTER	*1	LINE(130)/ 130 * '-'/

	LOGICAL         *1	NORMAL/.TRUE./
        LOGICAL         *1	FIRST/.TRUE./

	EXTERNAL		FGA_WRITE_ERR

	DICTIONARY	'FEX_IDX_FLAG'
	RECORD /FEX_IDX_FLAG/ IDXFLAGS

C
C Begin the dump.
C
	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6
	IF ( FIRST) THEN
          L = 0
	  WRITE (WRITE_LUN, 600, IOSTAT=IOSTATUS) FILE_NAME,CLINE,LINE
          FIRST = .FALSE.
        ENDIF
 600    FORMAT(	// 1X, 'FEX_IDX_FLAG FILE : ',A/,X,A80/,130A/)

	IF (IOSTATUS .NE. 0) THEN
	  PRINT *,' WRITE ERROR!!! UNIT = ', WRiTE_LUN 
	END IF

	WRITE ( WRITE_LUN,601) NREC
  601   FORMAT(/,T50,' REC # : ',I4/)
C
C Cobetrieve standard header.
C
	WRITE (WRITE_LUN, 6000, IOSTAT=STATUS)
	1			(IDXFLAGS.CT_HEAD.GMT(I:I), I=1,14),
	1			IDXFLAGS.CT_HEAD.TIME(2),
	1			IDXFLAGS.CT_HEAD.TIME(1),
	1			IDXFLAGS.CT_HEAD.SPACE_TIME,
	1			IDXFLAGS.CT_HEAD.MJR_FRM_NO,
	1			IDXFLAGS.CT_HEAD.ORBIT

6000	FORMAT (/' Processed time (in 3 formats): ' /
	1	' <GMT> : ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	' <Binary> : ', Z8.8, 1X, Z8.8, 7X,
	1	' <PB5> : ',  6(Z2.2, 1X) /
	1	' Major Frame Number : ', I6, 8x, 
	1	' Orbit number : ', I6)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	WRITE ( WRITE_LUN,1600)
 1600   FORMAT(//,10x,' Item Name ', 36x, 'Flag'/130('-') )  
	DO L = 1 , 284
	  WRITE (WRITE_LUN, 6001, IOSTAT=STATUS) ITEM_TEXT(L),
	1               			IDXFLAGS.IDX_FLAGS.FLAGS(L)
        END DO
 6001   FORMAT(4x,a,24x,l1)

	WRITE (WRITE_LUN,6002) LINE
 6002   FORMAT(/,130a)

        IF (STATUS .NE. 0 ) THEN
          NORMAL = .FALSE.
          PRINT *, ' Formatted write error ', WRITE_LUN
        END IF

        RETURN
        END
