	SUBROUTINE FGA_LIST_ANC (WRITE_LUN, ANC_REC, NREC, LINE,
	1	CAT_ENTRY, WRITE_ANC)

C-----------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FGA_LIST_ANC
C
C	PROGRAM DESCRIPTION:
C	  The program will get a formatted listing of all the header
C	  info in a given FIRAS ancillary housekeeping record.
C
C	AUTHOR:
C	  R. Kummerer
C	  STX
C	  April 27, 1989
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_ANC (WRITE_LUN, ANC_REC, NREC, LINE,
c		CAT_ENTRY, WRITE_ANC)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list
C					file.
C	  ANC_REC		RECORD	FIRAS ancillary housekeeping record
C	  NREC			I*2	Count of the number of records in this
C					run
C	  LINE			C*80	Comment line to identify the list file
C	  CAT_ENTRY		I*4	CT catalog entry #
C	  WRITE_ANC		L*1	Flag to write out housekeeping as well
C					as header info
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
C	  CT_BINARY_TO_GMT
C
C-----------------------------------------------------------------------

	IMPLICIT	NONE

	CHARACTER	*80	LINE

	INTEGER		*2	NREC
	INTEGER		*4	WRITE_LUN
	INTEGER		*4	CAT_ENTRY
	BYTE			WRITE_ANC
	INTEGER		*4	STATUS

	CHARACTER	*14	GMT

	INTEGER		*2	I
	INTEGER		*2	J
	INTEGER		*2	K
	INTEGER		*2	L

	DICTIONARY	'NFS_ANC'
	RECORD /NFS_ANC/ ANC_REC

	EXTERNAL		FGA_WRITE_ERR

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6

	WRITE (WRITE_LUN, 6000, IOSTAT=STATUS) LINE, NREC, CAT_ENTRY

6000	FORMAT ('1', / 5X, A // 
	1	' Records number in this run: ', I5,
	1	10x, 'CT Catalog Entry: ', i7 //
	1	15X, 'FIRAS Ancillary Housekeeping Record Formatted Dump' /)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Cobetrieve standard header.
C
	WRITE (WRITE_LUN, 6001, IOSTAT=STATUS)
	1			(ANC_REC.CT_HEAD.GMT(I:I), I=1,14),
	1			ANC_REC.CT_HEAD.TIME(2),
	1			ANC_REC.CT_HEAD.TIME(1),
	1			ANC_REC.CT_HEAD.SPACE_TIME,
	1			ANC_REC.CT_HEAD.MJR_FRM_NO,
	1			ANC_REC.CT_HEAD.ORBIT

6001	FORMAT (/' Processed time (in 3 formats): ' /
	1	' <GMT> : ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	' <Binary> : ', Z8.8, 1X, Z8.8, 7X,
	1	' <PB5> : ',  6(Z2.2, 1X) /
	1	' Major Frame Number : ', T30, I6 /
	1	' Orbit number : ', T30, I6)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Dump Science/Engineering Data.
C
	IF (WRITE_ANC) THEN

	  WRITE (WRITE_LUN, 6200, IOSTAT=STATUS) (J,J=1,10)

6200	  FORMAT (/,'***** First Major Frame Minor Frame Status Bits *****'
	1	   ///18X, 10(I2, 10X))

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF


	  DO J=0,11
	    K = J*10
	    WRITE (WRITE_LUN, 6400, IOSTAT=STATUS) K,
	1	(ANC_REC.FRAME(1).MINOR_FRAME_STATUS_BITS(J*10+L),L=1,10)
	    IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	    END IF
	  ENDDO

	  K = 120

	  WRITE (WRITE_LUN, 6400, IOSTAT=STATUS) K,
	1	(ANC_REC.FRAME(1).MINOR_FRAME_STATUS_BITS(J*10+L),L=1,8)

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

	  WRITE (WRITE_LUN, 6300, IOSTAT=STATUS) (J,J=1,10)

6300	  FORMAT (/,'***** Second Major Frame Minor Frame Status Bits *****'
	1	   ///18X, 10(I2, 10X))

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF


	  DO J=0,11
	    K = J*10
	    WRITE (WRITE_LUN, 6400, IOSTAT=STATUS) K,
	1	(ANC_REC.FRAME(2).MINOR_FRAME_STATUS_BITS(J*10+L),L=1,10)
	    IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	    END IF
	  ENDDO

	  K = 120

	  WRITE (WRITE_LUN, 6400, IOSTAT=STATUS) K,
	1	(ANC_REC.FRAME(2).MINOR_FRAME_STATUS_BITS(J*10+L),L=1,8)

6400	  FORMAT (1X, I7, 10(2X, Z10))

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

	ENDIF

C
C	Write GMT of 2nd major frame.
C
	CALL CT_BINARY_TO_GMT (ANC_REC.GMT_MJF2, GMT)

	WRITE (WRITE_LUN, 6500, IOSTAT=STATUS)
	1			GMT(1:2),
	1			GMT(3:5),
	1			GMT(6:7),
	1			GMT(8:9),
	1			GMT(10:11),
	1			GMT(12:14)

6500	FORMAT (//
	1	' GMT of 2nd Major frame: ', 
	1	A2, '-', A3, '-', 3(A2, '-'), A3)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	RETURN
	END
