	PROGRAM FTB_GETTIM

C-------------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FTB_GETTIM
C
C	PROGRAM DESCRIPTION:
C	  This program will inquire a user for an interval of FIRAS
C	  archive records, and then reads and dumps the timetags stored
C	  in each record.
C
C	AUTHOR:
C	  R. Kummerer
C	  STX
C	  March 28, 1988
C
C	  Adapted from work done by Ed Fung and Fred Shuman.
C
C	CALLING SEQUENCE:
C	  This is a main program
C
C	INPUT/OUTPUT PARAMETERS:
C	  NONE
C
C	INPUT FILES:
C	  FIRAS ARCHIVE
C
C	OUTPUT FILES:
C	  List files.
C
C	INCLUDE FILES USED:
C	  CTUSER.INC
C	  FUT_PARAMS
C	  CCT_STATUS_RECORD
C	  $SSDEF
C
C	SUBROUTINES CALLED:
C	  FUT_QGET_ENG
C	  FUT_QGET_HKP
C	  FUT_QGET_IDX
C	  FUT_QGET_EMF
C	  FUT_QGET_NOS
C	  FUT_QGET_ANC
C	  LIB$GET_LUN
C	  LIB$FREE_LUN
C	  STR$UPCASE
C
C---------------------------------------------------------------------------
C
C CHANGES:
C
C	Shirley M. Read
C	August 1988
C	STX
C	Reason: SPRs 1566 and 1763 requested that all FIRAS facilities 
C	        should tell of successful completion and signal errors 
C	        via $status for batch processing. Also added interface 
C		with Fut_Error, the condition handler and abort message
C	 	for Lib$Signal.
C
C	R. Kummerer, February  4, 1989. SPR 2848. Dump ancillary housekeeping.
C	R. Kummerer, May       4, 1989. SPR 3622. Appropriate dump file names.
C	R. Kummerer, July     22, 1989. SER 4176. Dump new archives FPP_SDF_xx
C						  and FDQ_SDF_xx.
C       N. Gonzales, June     29, 1992. SPR 9798. Deleted obsolete FIRAS
C                                                 datasets from FCI, FES, FCS,
C                                                 FFC, FPR, FPS, and FSF
C---------------------------------------------------------------------------

	IMPLICIT	NONE

C	Include Files

	INCLUDE		'CT$LIBRARY:CTUSER.INC'
	INCLUDE		'(CCT_STATUS_RECORD)'
	INCLUDE		'(FUT_PARAMS)'
	INCLUDE		'($SSDEF)'

C	Condition Handler

	EXTERNAL        FUT_ERROR

	LOGICAL		*1	NORMAL
	CHARACTER	*80	COMMENT
	CHARACTER	*16	FNAME
	CHARACTER	*80	LINE

	CHARACTER	*14	START_TIME
	CHARACTER	*14	END_TIME
	INTEGER		*4	STARTING_TIME(2)
	INTEGER		*4	ENDING_TIME(2)
	CHARACTER	*30	TIME_RANGE
	INTEGER		*2	DATASET_ID
	CHARACTER	*128	RSE(16)
	LOGICAL		*1	OLDRSE/.FALSE./

	CHARACTER	*8	ANS
	INTEGER		*4	DATATYP
	INTEGER		*4	NL

	INTEGER		*2	CT_LUN
	INTEGER		*2	CT_STAT(20)
	INTEGER		*4	CAT_ENTRY
	INTEGER		*4	STATUS
	INTEGER		*4	CHAN
	INTEGER		*2	ISTAT
	INTEGER		*2	NMJFR
	INTEGER		*2	NREC/0/
	INTEGER		*4	LUN
	INTEGER		*4	WRITE_LUN
	CHARACTER	*14	GMT

	INTEGER		*4	CUT_READ_DAFS_RECORD_IGSE_KEY
	INTEGER		*4	FUT_QGET_ENG
	INTEGER		*4	FUT_QGET_HKP
	INTEGER		*4	FUT_QGET_IDX
	INTEGER		*4	FUT_QGET_SCI
	INTEGER		*4	FUT_QGET_EMF
	INTEGER		*4	FUT_QGET_UCS
	INTEGER		*4	FUT_QGET_NOS
	INTEGER		*4	FUT_QGET_ANC
	INTEGER		*4	FUT_QGET_SSC
	INTEGER		*4	FUT_QGET_EXT
	INTEGER		*4	LIB$GET_LUN
	INTEGER		*4	LIB$FREE_LUN
	INTEGER		*4	STR$UPCASE

	RECORD /CCT_STATUS/ STRUC

	DICTIONARY 'CDU_DAFS'
	RECORD /CDU_DAFS/ DAFS_REC

	DICTIONARY 'NFS_HKP'
	RECORD /NFS_HKP/ HKP_REC
	DICTIONARY 'NFS_ANC'
	RECORD /NFS_ANC/ ANC_REC
	DICTIONARY 'NFS_SDF'
	RECORD /NFS_SDF/ SCI_REC
	DICTIONARY 'NFS_EMF'
	RECORD /NFS_EMF/ EMF_REC

	DICTIONARY 'FDQ_IDX'
	RECORD /FDQ_IDX/ IDX_REC
	DICTIONARY 'FDQ_ENG'
	RECORD /FDQ_ENG/ ENG_REC
	DICTIONARY 'FDQ_ETR'
	RECORD /FDQ_ETR/ ETR_REC
	DICTIONARY 'FXT_ENG_XTRM'
	RECORD /FXT_ENG_XTRM/ EXT_REC

	DICTIONARY 'FEC_SSCAL'
	RECORD /FEC_SSCAL/ SS_REC

	DICTIONARY 'FNT_NOISE'
	RECORD /FNT_NOISE/ NOS_REC

	EXTERNAL	CUT_NORMAL
	EXTERNAL	FUT_NORMAL
	EXTERNAL	FUT_EOF
	EXTERNAL	FTB_NORMAL
	EXTERNAL        FTB_ABERR
	EXTERNAL	FTB_OPEN
	EXTERNAL	FTB_WRITE
	EXTERNAL	FTB_CLOSE
	EXTERNAL	FTB_CTINIT
	EXTERNAL	FTB_CTOPEN
	EXTERNAL	FTB_CTQBLD
	EXTERNAL	FTB_CTCLOS


C     Establish condition handler.            

	CALL LIB$ESTABLISH ( FUT_ERROR )


5	CALL LIB$ERASE_PAGE (1, 1)

	TYPE *
	TYPE *,' FTB_GETTIME for producing timetag dumps'
        TYPE *,'             of the FIRAS archive'
	TYPE *
	TYPE *,' CSDR Version 4.2   November 1, 1988'
	TYPE *
	TYPE *

	NORMAL = .TRUE.

	CALL CT_INIT (STRUC)

	IF (STRUC.CTERR .EQ. CTP_NORMAL) THEN

	 CALL CT_QUERY_BLD(CTU_$FIRAS, TIME_RANGE, RSE, STRUC, OLDRSE)
	 IF (STRUC.CTERR .EQ. CTP_NORMAL) THEN

	   DATASET_ID = STRUC.DATASET_ID

	   STATUS = CUT_READ_DAFS_RECORD_IGSE_KEY(CTU_$FIRAS,DATASET_ID,
	2					  DAFS_REC)
	   IF (STATUS .EQ. %LOC(CUT_NORMAL)) THEN

	     NL = INDEX(DAFS_REC.DATASET,' ')
	     IF (NL .EQ. 0) THEN
	       NL = LEN(DAFS_REC.DATASET)
	     ELSE
	       NL = NL - 1
	     END IF

	     FNAME = DAFS_REC.DATASET(1:NL) // '.LIS'
	     DATATYP = -1

	     IF (DATASET_ID .EQ. CTU_$FIR_HKP) THEN
	       DATATYP = CTU_$FIR_HKP

	     ELSE IF ((DATASET_ID .GE. CTU_$FIR_RS1 .AND.
	1              DATASET_ID .LE. CTU_$FIR_RS4) .OR.
	2	      (DATASET_ID .GE. FAC_FPP_RH .AND.
	3	       DATASET_ID .LE. FAC_FPP_LL) .OR.
	4	      (DATASET_ID .GE. FAC_FDQ_RH .AND.
	5              DATASET_ID .LE. FAC_FDQ_LL)) THEN
	       DATATYP = CTU_$FIR_RS1

	     ELSE IF (DATASET_ID .GE. CTU_$FIR_SS1 .AND.
	1             DATASET_ID .LE. CTU_$FIR_SS4) THEN
	       DATATYP = CTU_$FIR_SS1

	     ELSE IF (DATASET_ID .GE. CTU_$FIR_ED1 .AND.
	1             DATASET_ID .LE. CTU_$FIR_ED4) THEN
	       DATATYP = CTU_$FIR_ED1

	     ELSE IF (DATASET_ID .EQ. CTU_$FIR_EDF .OR.
	1 	      DATASET_ID .EQ. FAC_HKP_EDF) THEN
	       DATATYP = CTU_$FIR_EDF

	     ELSE IF (DATASET_ID .EQ. CTU_$FIR_ETR) THEN
	       DATATYP = CTU_$FIR_ETR

	     ELSE IF (DATASET_ID .EQ. CTU_$FIR_IDX .OR.
	1   	      DATASET_ID .EQ. FAC_HKP_IDX) THEN
	       DATATYP = CTU_$FIR_IDX

	     ELSE IF (DATASET_ID .GE. CTU_$FIR_NRH .AND.
	1             DATASET_ID .LE. CTU_$FIR_NLL) THEN
	       DATATYP = CTU_$FIR_NRH
 
	     ELSE IF (DATASET_ID .EQ. FAC_FIR_ANC) THEN
	       DATATYP = FAC_FIR_ANC

	     ELSE IF (DATASET_ID .EQ. FAC_ENG_XTRM) THEN
	       DATATYP = FAC_ENG_XTRM

	     ELSE
	       TYPE 6002
	       GOTO 5
	     END IF

6002	     FORMAT (//' Invalid response, please try again!')

C
C Dump the records in the selected timerange.
C
	     STATUS = LIB$GET_LUN (LUN)

	     IF (STATUS .EQ. SS$_NORMAL) THEN

	       WRITE (6, 100)
100   	       FORMAT (//,' Make a hard copy dump (Y/[N], Default=TT)? ', $)
	       ACCEPT 20,ANS
 20    	       FORMAT (A)
	       STATUS = STR$UPCASE(ANS,ANS)

	       IF (ANS .EQ. 'Y') THEN
	         WRITE_LUN = LUN
	         OPEN (UNIT=WRITE_LUN, FILE=FNAME, TYPE='NEW', IOSTAT=ISTAT)
	       ELSE
	         WRITE_LUN = 6
	         ISTAT = 0
	       END IF

	       TYPE *

	       IF (ISTAT .EQ. 0) THEN

	         CALL CT_OPEN_ARCV(, CT_LUN, CTU_$FIRAS, DATASET_ID, CTU_$QNR,
	1                          STRUC, %REF(TIME_RANGE), )

	         IF (STRUC.CTERR .EQ. CTP_NORMAL) THEN

	           CALL CT_QUERY_ARCV(, CT_LUN, RSE, STRUC)

	  	   IF (STRUC.CTERR .EQ. CTP_NORMAL) THEN

	             DO WHILE (NORMAL)
C
C Dump housekeeping.
C
	               IF (DATATYP .EQ. CTU_$FIR_HKP) THEN
	                 STATUS = FUT_QGET_HKP (CT_LUN, HKP_REC)
	                 IF (STATUS .EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC =  NREC + 1
	     	           GMT = HKP_REC.CT_HEAD.GMT
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
C
C Dump science.
C
	               ELSE IF (DATATYP .EQ. CTU_$FIR_RS1) THEN
	                 STATUS = FUT_QGET_SCI (CT_LUN, SCI_REC)
	                 IF (STATUS .EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
 	    	           GMT = SCI_REC.CT_HEAD.GMT
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
C
C Dump short science.
C
	               ELSE IF (DATATYP .EQ. CTU_$FIR_SS1) THEN
	                 STATUS = FUT_QGET_SSC (CT_LUN, SS_REC)
	                 IF (STATUS .EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	     	           CALL CT_BINARY_TO_GMT (SS_REC.TIME, GMT)
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
C
C Dump engineering.
C
	               ELSE IF (DATATYP .EQ. CTU_$FIR_EDF) THEN
	                 STATUS = FUT_QGET_ENG (CT_LUN, ENG_REC)
	                 IF (STATUS .EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	  	           GMT = ENG_REC.CT_HEAD.GMT
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF

C
C Dump engineering trends.
C
	               ELSE IF (DATATYP .EQ. CTU_$FIR_ETR) THEN
	                 STATUS = FUT_QGET_ENG (CT_LUN, ETR_REC)
	                 IF (STATUS .EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	     	           GMT = ETR_REC.CT_HEAD.GMT
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF

C
C Dump engineeering mode.
C
	               ELSE IF (DATATYP .EQ. CTU_$FIR_ED1) THEN
	                 STATUS = FUT_QGET_EMF (CT_LUN, EMF_REC)
	                 IF (STATUS .EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	   	           GMT = EMF_REC.CT_HEAD.GMT
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF

C
C Dump engineering index.
C
	               ELSE IF (DATATYP .EQ. CTU_$FIR_IDX) THEN
	                 STATUS = FUT_QGET_IDX (CT_LUN, IDX_REC)
	                 IF (STATUS .EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
		           GMT = IDX_REC.HEADER.GMT_START_TIME_ASCII
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
C
C Dump noise spectra.
C
	               ELSE IF (DATATYP .EQ. CTU_$FIR_NRH) THEN
	                 STATUS = FUT_QGET_NOS (CT_LUN, NOS_REC)
	                 IF (STATUS .EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	   	           GMT = NOS_REC.CT_HEAD.GMT
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
C
C Dump ancillary housekeeping
C
	               ELSE IF (DATATYP .EQ. FAC_FIR_ANC) THEN
	                 STATUS = FUT_QGET_ANC (CT_LUN, ANC_REC)
	                 IF (STATUS .EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
		           GMT = ANC_REC.CT_HEAD.GMT
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
C
C Dump engineering extrema
C
	               ELSE IF (DATATYP .EQ. FAC_ENG_XTRM) THEN
	                 STATUS = FUT_QGET_EXT (CT_LUN, EXT_REC)
	                 IF (STATUS .EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	  	           GMT = EXT_REC.CT_HEAD.GMT
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
	               END IF
C
C Display the timetags.
C

	 	       IF (NORMAL) THEN

	 	         WRITE (WRITE_LUN, 6001) NREC, GMT
6001		         FORMAT (' Record number: ', I7, '  GMT:  ', A)

		       END IF

	             END DO

		   END IF

		   NORMAL = .TRUE.
	           IF (STATUS .NE. %LOC(FUT_EOF)) THEN
	             CALL LIB$SIGNAL (STATUS)
	   	     NORMAL = .FALSE.
	           END IF

	           CALL CT_CLOSE_ARCV (, CT_LUN, STRUC)

	           IF (STRUC.CTERR .NE. CTP_NORMAL) THEN
	             CALL LIB$SIGNAL (FTB_CTCLOS, %VAL(1), %VAL(STRUC.CTERR))
	             NORMAL = .FALSE.
	           END IF

	         ELSE

	           CALL LIB$SIGNAL (FTB_CTOPEN, %VAL(1), %VAL(STRUC.CTERR))
	           NORMAL = .FALSE.

	         END IF

	         IF (ANS .EQ. 'Y') THEN

	           CLOSE (UNIT=WRITE_LUN, IOSTAT=ISTAT)

	           IF (ISTAT .NE. 0) THEN
	             CALL LIB$SIGNAL (FTB_CLOSE, %VAL(1), %VAL(ISTAT))
	           END IF

	         END IF

	       ELSE

	         CALL LIB$SIGNAL (FTB_OPEN, %VAL(1), %VAL(ISTAT))

	       END IF

	     ELSE

	       CALL LIB$SIGNAL (%VAL(STATUS))
	       NORMAL = .FALSE.

	     END IF

	   ELSE

	     CALL LIB$SIGNAL (FTB_CTQBLD, %VAL(1), %VAL(STRUC.CTERR))
	     NORMAL = .FALSE.

	   END IF

	 ELSE

	   CALL LIB$SIGNAL (%VAL(STATUS))
	   NORMAL = .FALSE.

	 END IF

	 IF (ANS .EQ. 'Y') THEN

	   STATUS = LIB$FREE_LUN(LUN)

	   IF (STATUS .NE. SS$_NORMAL) THEN
	     CALL LIB$SIGNAL (%VAL(STATUS))
	     NORMAL = .FALSE.
	   END IF

         END IF

	ELSE

	  CALL LIB$SIGNAL (FTB_CTINIT, %VAL(1), %VAL(STRUC.CTERR))
	  NORMAL = .FALSE.

	END IF


	IF (NREC .GT. 0 .AND. ANS .EQ. 'Y') THEN
	  WRITE(6, 111) FNAME
111	  FORMAT (/' List file is named ', A  //)
	END IF

C       Exit the program.

	IF (NORMAL) THEN
	  CALL LIB$SIGNAL(FTB_NORMAL)
	  CALL EXIT(SS$_NORMAL)
	ELSE
	  CALL LIB$SIGNAL(FTB_ABERR)
	  CALL EXIT(SS$_ABORT)
	END IF

	END
