	Integer*4 Function FGA_LIST_VABSAA (WRITE_LUN, VABSAA_REC,
	1                  VABSAA_FILE, CLINE, NREC)

C---------------------------------------------------------------------------
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given FIRAS 
C	  reference data of VABSAA.
C
C	AUTHOR:
C	  Nilo G. Gonzales
C	  STX
C	  September 23, 1991  Reference: SPR 8372
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_VABSAA FGA_LIST_VABSAA (WRITE_LUN, VABSAA_REC,
C                                VABSAA_FILE, CLINE, NREC)
C	INPUT PARAMETERS:
C	  WRITE_LUN	 I*4    Logical unit number for a list file;
C				if not supplied assume 6.  Assume the
C				calling routine has opened this list file.
C	  VABSAA_REC         	VABSAA record.
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
	INTEGER		*4	STATUS, IOSTAT
	CHARACTER	*14	GMT
	CHARACTER       *64     VABSAA_FILE
	INTEGER		*2	I, J, K, L, N
	CHARACTER       *80	CLINE
	CHARACTER	*1	LINE(130)/ 130 * '-'/
	EXTERNAL		FGA_WRITE_ERR

	DICTIONARY	'FEX_VABSAA'
	RECORD /FEX_VABSAA/ VABSAA_REC

        FGA_LIST_VABSAA = SUCCESS

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6
	IF ( FIRST) THEN
	  WRITE (WRITE_LUN, 10, IOSTAT=STATUS) VABSAA_FILE,CLINE,LINE
          FIRST = .FALSE.
        ENDIF
 10   FORMAT(	// 1X, 'FEX_VABSAA FILE : ',A/,X,A80/,130A/)

	IF (IOSTAT .NE. 0) THEN
	  PRINT *,' WRITE ERROR!!! UNIT = ', WRITE_LUN 
          FGA_LIST_VABSAA = ERROR          
	END IF

	WRITE ( WRITE_LUN,20) NREC
  20   FORMAT(/,T30,' REC # : ',I4/)
C
C Cobetrieve standard header.
C
	WRITE (WRITE_LUN, 70, IOSTAT=STATUS)
	1     (VABSAA_REC.CT_HEAD.GMT(I:I), I=1,14),
	1     VABSAA_REC.CT_HEAD.TIME(2),
	1     VABSAA_REC.CT_HEAD.TIME(1),
	1     VABSAA_REC.CT_HEAD.SPACE_TIME,
	1     VABSAA_REC.CT_HEAD.MJR_FRM_NO,
	1     VABSAA_REC.CT_HEAD.ORBIT,
	1     VABSAA_REC.CT_HEAD.HSKP1_TLM_FMT,
	1     VABSAA_REC.CT_HEAD.HSKP2_TLM_FMT

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
C
C Dump MTM VAB record.
C
	WRITE (WRITE_LUN, 80, IOSTAT=STATUS)
	1     VABSAA_REC.VAB(1).LONSTEP,
	1     VABSAA_REC.VAB(1).LATMIN,
	1     VABSAA_REC.VAB(1).LATMAX

80	FORMAT (//' ***** VAB Longitude/Latitude North boundary *****'/,
	1      5x,'Longitude step interval: ', T30, G14.5 /
	1      5x,'Minimum Latitude: ', T30, G14.5 /
	1      5x,'Maximum Latitude: ', T30, G14.5)

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	WRITE (WRITE_LUN, 90, IOSTAT=STATUS)
	1	VABSAA_REC.VAB(2).LONSTEP,
	1       VABSAA_REC.VAB(2).LATMIN,
	1       VABSAA_REC.VAB(2).LATMAX

90	FORMAT (//' ***** VAB Longitude/Latitude South boundary *****'/,
	1	5x,'Longitude step interval: ', T30, G14.5 /
	1       5x,'Minimum Latitude: ', T30, G14.5 /
	1       5x,'Maximum Latitude: ', T30, G14.5)

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	WRITE (WRITE_LUN, 100, IOSTAT=STATUS)
100	FORMAT (//' ***** VAB North Latitude values/North boundary *****'/)

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	DO J=0,4
	   K = J*8
	   WRITE (WRITE_LUN, 200, IOSTAT=STATUS) K,
	1        (VABSAA_REC.VAB(1).LATN(J*8+L),L=1,8)
200	   FORMAT (1X, I3, 8(2X, G14.6))

	   IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	   END IF
	ENDDO

	N = 40
	WRITE (WRITE_LUN, 250, IOSTAT=STATUS) N,
	1     VABSAA_REC.VAB(1).LATN(41)
250     FORMAT(1X, I3,2X, G14.6)	   

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	WRITE (WRITE_LUN, 500, IOSTAT=STATUS)
500	FORMAT (//' ***** VAB North Latitude values/South boundary *****'/)

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	DO J=0,4
	   K = J*8
	   WRITE (WRITE_LUN, 600, IOSTAT=STATUS) K,
	1        (VABSAA_REC.VAB(1).LATS(J*8+L),L=1,8)
600	   FORMAT (1X, I3, 8(2X, G14.6))

	   IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	   END IF
	ENDDO

	N = 40
	WRITE (WRITE_LUN, 650, IOSTAT=STATUS) N,
	1      VABSAA_REC.VAB(1).LATS(41)
650     FORMAT(1X, I3,2X, G14.6)	   

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	WRITE (WRITE_LUN, 300, IOSTAT=STATUS)
300	FORMAT (//' ***** VAB South Latitude values/North boundary *****'/)

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	DO J=0,4
	   K = J*8
	   WRITE (WRITE_LUN, 400, IOSTAT=STATUS) K,
	1        (VABSAA_REC.VAB(2).LATN(J*8+L),L=1,8)
400	   FORMAT (1X, I3, 8(2X, G14.6))

	   IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	   END IF
	ENDDO

	N = 40
	WRITE (WRITE_LUN, 450, IOSTAT=STATUS) N,
	1      VABSAA_REC.VAB(2).LATN(41)
450     FORMAT(1X, I3,2X, G14.6)	   

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	WRITE (WRITE_LUN, 700, IOSTAT=STATUS)
700	FORMAT (//' ***** VAB South Latitude values/South boundary *****'/)

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	DO J=0,4
	   K = J*8
	   WRITE (WRITE_LUN, 800, IOSTAT=STATUS) K,
	1        (VABSAA_REC.VAB(2).LATS(J*8+L),L=1,8)
800	   FORMAT (1X, I3, 8(2X, G14.6))

	   IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	   END IF
	ENDDO

	N = 40
	WRITE (WRITE_LUN, 850, IOSTAT=STATUS) N,
	1      VABSAA_REC.VAB(2).LATS(41)
850     FORMAT(1X, I3,2X, G14.6)	   

	WRITE (WRITE_LUN, 900, IOSTAT=STATUS)
	1      VABSAA_REC.SAA.LONSTEP,
	1      VABSAA_REC.SAA.LONMIN,
	1      VABSAA_REC.SAA.LONMAX,
	1      VABSAA_REC.SAA.LATMIN,
	1      VABSAA_REC.SAA.LATMAX

900	FORMAT (//' ***** Information of South Alantic Anomaly *****'/,
	1	5x,'Longitude step interval: ', T30, G14.5 /
	1       5x,'Minimum Longitude: ', T30, G14.5 /
	1       5x,'Maximum Longitude: ', T30, G14.5 /
	1       5x,'Minimum Latitude: ', T30, G14.5 /
	1       5x,'Maximum Latitude: ', T30, G14.5 )

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	WRITE (WRITE_LUN, 1000, IOSTAT=STATUS)
1000	FORMAT (//' ***** SAA Latitude values/North boundary *****'/)

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	DO J=0,4
	   K = J*6
	   WRITE (WRITE_LUN, 2000, IOSTAT=STATUS) K,
	1        (VABSAA_REC.SAA.LATN(J*6+L),L=1,6)
2000	   FORMAT (1X, I3, 6(3X, G14.6))

	   IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	   END IF
	ENDDO

	N = 30
	WRITE (WRITE_LUN, 2500, IOSTAT=STATUS) N,
	1      VABSAA_REC.SAA.LATN(31)
2500    FORMAT(1X, I3,3X, G14.6)	   

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	WRITE (WRITE_LUN, 3000, IOSTAT=STATUS)
3000	FORMAT (//' ***** SAA Latitude values/South boundary *****'/)

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	DO J=0,4
	   K = J*6
	   WRITE (WRITE_LUN, 4000, IOSTAT=STATUS) K,
	1        (VABSAA_REC.SAA.LATS(J*6+L),L=1,6)
4000	   FORMAT (1X, I3, 6(3X, G14.6))

	   IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	   END IF
	ENDDO

	N = 30
	WRITE (WRITE_LUN, 4500, IOSTAT=STATUS) N,
	1      VABSAA_REC.SAA.LATS(31)
4500    FORMAT(1X, I3,3X, G14.6)	   

	IF (STATUS .NE. 0) THEN
	   CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	RETURN
	END
