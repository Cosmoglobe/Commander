	SUBROUTINE Ftb_LIST_SCI (WRITE_LUN, SCI_REC, NREC, LINE,
	1	CAT_ENTRY, WRITE_IFG)
C-----------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  Ftb_LIST_SCI
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
C         Modified by
C	  Q. CHUNG
C	  STX
C	  OCT. 10, 1989
C
C	CALLING SEQUENCE:
C	  CALL FTB_LIST_SCI (WRITE_LUN, SCI_REC, NREC, LINE,
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
C       SPR xxxx, OCT 5, 1989, Q. Chung  modify to dump Science data
c                                        quality flag in bits .
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

	CHARACTER	*512	GLITCH
	CHARACTER       *14     COLLECT_GMT
        Character       *1      c2bits(2,8),c18bits(18,8),c4bits(4,8)
        Integer         *2      B_quality(110)/110*0/
        CHARACTER       *2      WORK
        INTEGER        *2      IWORK/0/

	DICTIONARY	'NFS_SDF'
	RECORD /NFS_SDF/ SCI_REC

	EXTERNAL		FTB_WRITE

	EQUIVALENCE	(SW_VERSION_EQV, SW_VERSION)
	EQUIVALENCE	(CHAN_IND_EQV, CHAN_IND)
        EQUIVALENCE      (IWORK, WORK(1:1))
	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6

	WRITE (WRITE_LUN, 6000, IOSTAT=STATUS) LINE, NREC, CAT_ENTRY

6000	FORMAT ( / 5X, A // 
	1	' Records number in this run: ', I5,
	1	10x, 'CT Catalog Entry: ', i7 //
	1	15X, 'FIRAS Science Record Formatted Dump' /)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FTB_WRITE, %VAL(1), %VAL(STATUS) )
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
	1			SCI_REC.COLLECT_TIME.MIDPOINT_TIME(2),
	1			SCI_REC.COLLECT_TIME.MIDPOINT_TIME(1),
	1			SCI_REC.COLLECT_TIME.BADTIME_FLAG

6001	FORMAT (/' Processed time (in 3 formats): ' /
	1	' <GMT> : ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	' <Binary> : ', Z8.8, 1X, Z8.8, 7X,
	1	' <PB5> : ',  6(Z2.2, 1X) /
	1	' Major Frame Number : ', I6, 8x, 
	1	'<Midpoint of Collect> : ',2A1,'-',3A1,'-',3(2A1,'-'),3A1/
	1	' Orbit number : ', I6, 14x,
	1	'<Binary Collect> : ', Z8.8, 1x, Z8.8, 6x,' Time Flag: ',I4)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FTB_WRITE, %VAL(1), %VAL(STATUS) )
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

	SC13 = SCI_REC.SCI_HEAD.SC_HEAD13
	IF (SC13 .LT. 0) THEN
	   SC13 = SC13 + 65536
	END IF

	MIN_FRM_CNT = 65536 * SCI_REC.SCI_HEAD.SC_HEAD12 + SC13

	SC5 = SCI_REC.SCI_HEAD.SC_HEAD5
	IF (SC5 .LT. 0) THEN
	   SC5 = SC5 + 65536
	END IF

	DATA_TRANS_TIME = 65536 * SCI_REC.SCI_HEAD.SC_HEAD4 + SC5

	CHAN_IND_EQV = SCI_REC.SCI_HEAD.SC_HEAD25

	CALL CT_BINARY_TO_GMT (SCI_REC.DQ_DATA.ENG_TIME, GMT)

	WRITE (WRITE_LUN, 6010, IOSTAT=STATUS)
	1			SCI_REC.DQ_DATA.FAKE,
	1			SCI_REC.DQ_DATA.XCAL_POS,
	1			SCI_REC.DQ_DATA.IREF_TEMP,
	1			(GMT(I:I),I=1,14),
	1			SCI_REC.DQ_DATA.ENG_TIME(2),
	1			SCI_REC.DQ_DATA.ENG_TIME(1),
	1			SCI_REC.DQ_DATA.ENG_REC


6010	FORMAT (//'***** Data Qualify Information *****' /
	1	' Fake-it Mode: ', T30, I6 /
	1	' XCal Position: ', T30, I6 /
	1	' ICal Temperature: ', T26, G14.5 /
	1	' Engineering Time: ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	' <Binary>: ', Z8.8, 1X, Z8.8 /
	1	' Engineering Record Number: ', T30, I6 /
	1	' Instrument Data Quality: ', T30, I6 /
	1	' Attitude Data Quality: ', T30, I6 /)

	                DO l = 1 , 110
	                  iwork=0
	                 WRITE(WORK(1:1),'(a1)') SCI_REC.DQ_DATA.data_QUALITY(l)
	                  B_QUALITY(l)= IWORK
	                ENDDO
                

	                J = 0
	  		DO l = 4 , 5
	                   J = J + 1
	                  iwork=0
                        WRITE(WORK(1:1), '(a1)') SCI_REC.DQ_DATA.data_QUALITY(l)
                        b_quality(l)=IWORK
                       do k = 0 , 7
                         c2bits(j,8-k) = '0'
                         if ( btest(b_quality(l),k) ) c2bits(j,8-k) = '1'
                       enddo
	                enddo

	                J = 0
	  		Do l = 11 , 28
                          J = J + 1
                          iwork=0
                       WRITE(WORK(1:1), '(a1)') SCI_REC.DQ_DATA.data_QUALITY(l)
                       b_quality(L) = iwork
                       do k = 0,7
                         c18bits(j,8-k) = '0'
                         if ( btest(b_quality(l),k) ) c18bits(j,8-k) = '1'
                       enddo
	                enddo

	                J = 0
	  		Do l = 58 , 61
	                   J = J + 1
	                  iwork=0
                    WRITE(WORK(1:1), '(a1)') SCI_REC.DQ_DATA.data_QUALITY(l)
                     b_quality(l) = IWORK
                       do k = 0,7
                         c4bits(j,8-k) = '0'
                         if ( btest(b_quality(l),k) ) c4bits(j,8-k) = '1'
                       enddo
	                enddo

	write (write_lun,6020)
	1	(B_QUALITY(I),I = 1,3),
	1       ((c2bits(i,j),j=1,8),i=1,2),
	1	(B_QUALITY(I),I = 6,10),
	1       ((c18bits(i,j),j=1,8),i=1,18),
	1	(B_QUALITY(I),I = 29,30),
	1	(B_QUALITY(I), I = 31,50),
	1	(B_QUALITY(I), I = 51,57),
	1       ((c4bits(i,j),j=1,8),i=1,4),
	1	(B_QUALITY(I),I = 62,110)
 6020   format (1x,' Data Quality Flags: ' //
	1	'    1 -  3: ', 5X, 3(I4, 5X) /
	1       '    4 -  5: ', 5x, 2(8a,5X) /
	1       '    6 - 10: ', 5x, 5(i4, 5x) /
	1	'   11 - 16: ', 5X, 6(8a,1x) /
	1	'   17 - 22: ', 5X, 6(8a, 1X) /
	1	'   23 - 28: ', 5X, 6(8a, 1X) /
	1       '   29 - 30: ', 5x, 2(I4,5x) /
	1	'   31 - 40: ', 5X, 10(I4, 5X) /
	1	'   41 - 50: ', 5X, 10(I4, 5X) /
	1	'   51 - 57: ', 5X, 7(I4, 5X) /
	1       '   58 - 61: ', 5x, 4(8a,5x) /
	1	'   62 - 70: ', 5X, 9(I4, 5X) /
	1	'   71 - 80: ', 5X, 10(I4, 5X) /
	1	'   81 - 90: ', 5X, 10(I4, 5X) /
	1	'  91 - 100: ', 5X, 10(I4, 5X) /
	1	' 101 - 110: ', 5X, 10(I4, 5X) )


	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( Ftb_WRITE, %VAL(1), %VAL(STATUS) )
	END IF


	IF (WRITE_IFG) THEN

	  WRITE (WRITE_LUN, 6200, IOSTAT=STATUS) (J,J=1,10)

6200	  FORMAT ('1 ***** Science Interferogram *****'/
	1	  '             (* = Glitch)'
	1	   ///18X, 10(I2, 10X))

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FTB_WRITE, %VAL(1), %VAL(STATUS) )
	  END IF


	  DO J=0,50
	    K = J*10
	    WRITE (WRITE_LUN, 6101, IOSTAT=STATUS) K,
	1	(SCI_REC.IFG_DATA.IFG(J*10+L),GLITCH(K+L:K+L),L=1,10)
	    IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FTB_WRITE, %VAL(1), %VAL(STATUS) )
	    END IF
	  ENDDO

	  K = 510

	  WRITE (WRITE_LUN, 6101, IOSTAT=STATUS) K,
	1		SCI_REC.IFG_DATA.IFG(511), GLITCH(511:511),
	1	 	SCI_REC.IFG_DATA.IFG(512), GLITCH(512:512)

6101	  FORMAT (1X, I7, 10(2X, I8, 1X, A1))

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FTB_WRITE, %VAL(1), %VAL(STATUS) )
	  END IF

	ENDIF

	RETURN
	END
