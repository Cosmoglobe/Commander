	SUBROUTINE FGA_LIST_NOS (WRITE_LUN, NOS_REC, NREC, LINE,
	1	CAT_ENTRY, WRITE_SPC, WRITE_SIG)

C---------------------------------------------------------------------------
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given FIRAS
C	  noise spectrum record.
C
C	AUTHOR:
C	  R. Kummerer
C	  STX
C	  November 10, 1987
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_NOS (WRITE_LUN, NOS_REC, NREC, LINE, CAT_ENTRY)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file
C	  NOS_REC		RECORD	FIRAS noise spectra record.
C	  NREC			I*2	The number of records processed in this run.
C	  LINE			C*80	Comment line for listing
C	  CAT_ENTRY		I*4	CT catalog entry number
C	  WRITE_SPC		L*1	Dump the spectrum.
C	  WRITE_SIG		L*1	Dump the spectrum sigmas.
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

	IMPLICIT	NONE

	INCLUDE		'(FUT_PARAMS)'

	INTEGER		*2	NREC
	INTEGER		*4	WRITE_LUN
	INTEGER		*4	CAT_ENTRY
	INTEGER		*4	I
	INTEGER		*4	J
	INTEGER		*4	K
	INTEGER		*4	L
	CHARACTER	*80	LINE
	CHARACTER	*14	GMT
	INTEGER		*4	STATUS
	INTEGER		*4	CHAN
	LOGICAL		*1	WRITE_SPC
	LOGICAL		*1	WRITE_SIG

	EXTERNAL		FGA_WRITE_ERR

	DICTIONARY 'FNT_NOISE'
	RECORD /FNT_NOISE/ NOS_REC

C
C Begin the dump.
C
	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6

	WRITE (WRITE_LUN, 6000, IOSTAT=STATUS)	LINE, NREC, CAT_ENTRY
6000	FORMAT ( '1' /3x, a // ' Record Count for this Run: ', I6,
	1	10x, 'CT catalog entry: ', i7
	1	// 15X, 'FIRAS Noise Spectrum Formatted Listing' //)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C	
C Header information.
C
	WRITE (WRITE_LUN, 6001, IOSTAT=STATUS)
	1		(NOS_REC.CT_HEAD.GMT(I:I), I=1,14),
	1		NOS_REC.CT_HEAD.TIME(2),
	1		NOS_REC.CT_HEAD.TIME(1),
	1		NOS_REC.CT_HEAD.SPACE_TIME,
	1		(NOS_REC.COA_HEAD.FIRST_GMT(I:I), I=1,14),
	1		NOS_REC.COA_HEAD.FIRST_BIN_GMT(2),
	1		NOS_REC.COA_HEAD.FIRST_BIN_GMT(1),
	1		NOS_REC.COA_HEAD.FIRST_SPACE_TIME,
	1		NOS_REC.COA_HEAD.FIRST_MJR_FRM_NUM,
	1		(NOS_REC.COA_HEAD.LAST_GMT(I:I), I=1,14),
	1		NOS_REC.COA_HEAD.LAST_BIN_GMT(2),
	1		NOS_REC.COA_HEAD.LAST_BIN_GMT(1),
	1		NOS_REC.COA_HEAD.LAST_SPACE_TIME,
	1		NOS_REC.COA_HEAD.LAST_MJR_FRM_NUM,
	1		NOS_REC.COA_HEAD.NUM_FIRST_IFG,
	1		NOS_REC.COA_HEAD.NUM_LAST_IFG,
	1		NOS_REC.COA_HEAD.LABEL,
	1		NOS_REC.ATTIT.PIXEL_NO,
	1		NOS_REC.COA_HEAD.IREF_TEMP,
	1		NOS_REC.COA_HEAD.XCAL_POS,
	1		NOS_REC.COA_HEAD.NUM_IFGS

6001	FORMAT (/' Processed time (in 3 formats): ' /
	1	' <GMT> : ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	' <Binary> : ', Z8.8, 1X, Z8.8, 7X,
	1	' <PB5> : ',  6(Z2.2, 1X) //
	1	' Information for the first coadded IFG: ' /
	1	' <GMT> : ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	' <Binary> : ', Z8.8, 1X, Z8.8, 7X,
	1	' <PB5> : ',  6(Z2.2, 1X), 7X, ' MJR FR #: ', I11 //
	1	' Information for the last coadded IFG: ' /
	1	' <GMT> : ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	' <Binary> : ', Z8.8, 1X, Z8.8, 7X,
	1	' <PB5> : ',  6(Z2.2, 1X), 7X, ' MJR FR #: ', I11 //
	1	' First Coadded IFG #: ', I12, 15X,
	1	' Last Coadded IFG #: ', I12 //
	1	' Coadd label: ', A / 14X, 7('1234567890') //
	1	' Pixel #: ', T61, I6 /
	1	' Internal reference source temp. X 1000 (deg. K): ', 
	1	T59, G12.4 /
	1	' External calibrator position: ', T61, I6 /
	1	' Number of science IFG''s in coadd: ', T61, I6 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C IFG specific data.
C
	CHAN = NOS_REC.CHAN.CHAN_ID

	CALL CT_BINARY_TO_GMT (NOS_REC.CHAN.AV_COLLECT_TIME,
	1				GMT)

	WRITE (WRITE_LUN, 6010, IOSTAT=STATUS)
	1	FAC_CHANNEL_IDS(CHAN),
	1	NOS_REC.CHAN.AV_COLLECT_TIME(2),
	1	NOS_REC.CHAN.AV_COLLECT_TIME(1),
	1	(GMT(I:I),I=1,14)

6010	FORMAT (//' Coadd Time ', A2, ': ', Z8.8, 2X, Z8.8, 
	1	'  (', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, ')'/)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6995, IOSTAT=STATUS) FAC_CHANNEL_IDS(CHAN)

6995	FORMAT (/' ***** IFG Specific Data:  Channel ', A2, ' *****')

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 7000, IOSTAT=STATUS)
	1			NOS_REC.CHAN.CHAN_ID,
	1			NOS_REC.CHAN.NGROUP,
	1			NOS_REC.CHAN.MTM_SPEED,
	1			NOS_REC.CHAN.MTM_LENGTH,
	1			NOS_REC.CHAN.GAIN,
	1			NOS_REC.CHAN.SCI_MODE,
	1			NOS_REC.CHAN.FAKEIT,
	1			NOS_REC.CHAN.BOL_CMD_BIAS,
	1			NOS_REC.CHAN.BOL_VOLT,
	1			NOS_REC.CHAN.SWEEPS

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


7000	FORMAT (' Channel Number:            ', I6/
	1	' Adds Per Group:            ', I6/
	1	' MTM Speed:                 ', I6/
	1	' MTM Length:                ', I6/
	1	' Gain:                      ', I6/
	1	' Microprocessor Mode:       ', I6/
	1	' Fake-it Mode:              ', I6/
	1	' Commanded Bolometer Bias:  ', I6/
	1	' Bolometer Voltage:         ', G14.5/
	1	' Number of Sweeps:          ', I6)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Dump the attitude.
C
	CALL FGA_LIST_ATT (WRITE_LUN, NOS_REC.ATTIT)

C
C Dump the spectrum.
C
	IF (WRITE_SPC) THEN

	    WRITE (WRITE_LUN, 6400, IOSTAT=STATUS) (J,J=1,8)

6400	    FORMAT ('1 ***** Noise Spectrum *****'/
	1	    ///14X, 8(I2, 13X))

	    IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	    END IF


	    DO J=0,63
	      K = J*8
	      WRITE (WRITE_LUN, 6401, IOSTAT=STATUS) K,
	1		  (NOS_REC.CHAN.SPEC(J*8+L),L=1,8)
6401	      FORMAT (1X, I3, 8(2X, G13.6))
	      IF (STATUS .NE. 0) THEN
	        CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	      END IF
	    ENDDO

	    IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	    END IF

	ENDIF

C
C Dump the spectrum sigmas.
C
	IF (WRITE_SIG) THEN

	    WRITE (WRITE_LUN, 6500, IOSTAT=STATUS) (J,J=1,8)

6500	    FORMAT ('1 ***** Noise Spectrum Sigmas *****'/
	1	    ///14X, 8(I2, 13X))

	    IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	    END IF


	    DO J=0,63
	      K = J*8
	      WRITE (WRITE_LUN, 6501, IOSTAT=STATUS) K,
	1		  (NOS_REC.CHAN.SIGMAS(J*8+L),L=1,8)
6501	      FORMAT (1X, I3, 8(2X, G13.6))
	      IF (STATUS .NE. 0) THEN
	        CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	      END IF
	    ENDDO

	    IF (STATUS .NE. 0) THEN
	      CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	    END IF

	ENDIF


	RETURN
	END
