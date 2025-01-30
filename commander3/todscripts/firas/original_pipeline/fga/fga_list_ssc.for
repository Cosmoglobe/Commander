	SUBROUTINE FGA_LIST_SSC (WRITE_LUN, SSC_REC, NREC, LINE, CAT_ENTRY)

C---------------------------------------------------------------------------
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given FIRAS
C	  Short Science record.
C
C	AUTHOR:
C	  R. Kummerer
C	  STX
C	  March 28, 1988
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_SSC (WRITE_LUN, SSC_REC, NREC, LINE, CAT_ENTRY)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file
C	  SSC_REC		RECORD	FIRAS record.
C	  NREC			I*2	The number of records processed in this run.
C	  LINE			C*80	Comment line for listing
C	  CAT_ENTRY		I*4	CT catalog entry number
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
C	  FUT_PARMS.TXT
C
C	SUBROUTINES CALLED:
C	  CT_BINARY_TO_GMT
C
C	ERROR HANDLING:
C	  NONE
C---------------------------------------------------------------------------
Ch      N. Gonzales, Oct. 2, l990. SPR 7482.
C          Add fields SCI_MODE, PIXEL_DEFINITION, 
C                     SKY_MAP_INDEX, EXC_GALACTIC_LAT, IFG_NO
C
Ch      N. Gonzales, Dec. 5, 1990. SPR 7764.
C          Add field RAW_SCI_EXT
C
C       N. Gonzales, April 25, 1991. SPR 8404.
C          Deletted fields IFG_NO, GALACTIC_LAT and RAW_SCI_EXT.
C          Add fields Data_Quality, Attitude_quality, Xcal_temp,
C                     Skyhorn_temp, Ical_temp, Dihedral_temp, 
C                     Bolometer_temp and Bol_cmd_bias.
C---------------------------------------------------------------------------

	IMPLICIT	NONE

	INCLUDE         '(FUT_PARAMS)'

	INTEGER		*2	NREC
	INTEGER		*4	WRITE_LUN
	INTEGER		*4	CAT_ENTRY
	INTEGER		*4	I
	INTEGER		*4	J
	INTEGER		*4	K
	INTEGER		*4	L
	INTEGER		*4	M
	CHARACTER	*80	LINE
	CHARACTER	*14	GMT(3)
	INTEGER         *4      INDEX
	INTEGER         *4      GEO
	REAL            *4      ORBIT
	REAL            *4      SCAN
	REAL            *4      GLAT
	INTEGER		*4	STATUS

	EXTERNAL		FGA_WRITE_ERR


	DICTIONARY	'FEC_SSCAL'
	RECORD /FEC_SSCAL/ SSC_REC


C
C Begin the dump.
C
	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6

	WRITE (WRITE_LUN, 6000, IOSTAT=STATUS)	LINE, NREC, CAT_ENTRY
6000	FORMAT ( '1' /3x, a // ' Record Count for this Run: ', I6,
	1	10x, 'CT Catalog Entry: ', I7
	1	// 15X, 'FIRAS Short Science Formatted Listing' //)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C	
C Header information.
C
	CALL CT_BINARY_TO_GMT ( SSC_REC.TIME, GMT(1) )
	
	WRITE (WRITE_LUN, 6001, IOSTAT=STATUS)
	1		(GMT(1)(I:I), I=1,14),
	1		SSC_REC.TIME(2),
	1		SSC_REC.TIME(1),
	1		SSC_REC.CHANNEL_ID,
	1		SSC_REC.PIXEL_NO,
	1		SSC_REC.MTM_SCAN_SPEED,
	1		SSC_REC.MTM_SCAN_LENGTH,
	1               SSC_REC.SCI_MODE,
	1		SSC_REC.ADDS_PER_GROUP,
	1		SSC_REC.TRANSITION_FLAG

6001	FORMAT (/' <GMT> : ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	 ' <Binary> : ', Z8.8, 1X, Z8.8, 7X, /,
	1	 ' Channel : ', I6, /,
	1	 ' Pixel Number : ', I6, /,
	1	 ' MTM Scan Speed : ', I6, /,
	1	 ' MTM Scan Length : ', I6, /,
	1        ' Science Mode: ', I6, /,
	1	 ' Adds per Group : ', I6, /,
	1	 ' Transition Flag : ', I6 )

	If (SSC_REC.PIXEL_DEFINITION .Eq. 'q') Then
	  INDEX = SSC_REC.SKYMAP_INDEX

	Else If (SSC_REC.PIXEL_DEFINITION .Eq. 'O') Then
	  ORBIT = SSC_REC.SKYMAP_INDEX / 10.0

	Else If (SSC_REC.PIXEL_DEFINITION .Eq. 'S') Then
	  SCAN = SSC_REC.SKYMAP_INDEX / 10.0

	Else If (SSC_REC.PIXEL_DEFINITION .Eq. 'E') Then
	  GEO = SSC_REC.SKYMAP_INDEX
	End If

	If (SSC_REC.PIXEL_DEFINITION .Eq. 'q') Then
 	    WRITE (WRITE_LUN,6010,IOSTAT=STATUS) index

	Else If (SSC_REC.PIXEL_DEFINITION .Eq. 'O') Then
	    WRITE (WRITE_LUN,6020,IOSTAT=STATUS) orbit

	Else If (SSC_REC.PIXEL_DEFINITION .Eq. 'S') Then
	    WRITE (WRITE_LUN,6030,IOSTAT=STATUS) scan

	Else If (SSC_REC.PIXEL_DEFINITION .Eq. 'E') Then
	    WRITE (WRITE_LUN,6040,IOSTAT=STATUS) geo
	End If

6010	FORMAT (/1x,'Standard Pixels - Skymap_Index : ',I4)

6020	FORMAT (/1x,'Orbit Averaged Pixels - Orbit Origin : ',F9.1)

6030	FORMAT (/1x,'Scan Angle Pixels - Scan Angle : ',F9.1)

6040	FORMAT (/1x,'Geocentric Pixels; Skymap_Index : ',I4)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

 	WRITE (WRITE_LUN,6070,IOSTAT=STATUS)
	1         ( SSC_REC.DATA_QUALITY,
	1           SSC_REC.ATTITUDE_QUALITY,
	1           SSC_REC.XCAL_TEMP,
	1           SSC_REC.SKYHORN_TEMP,
	1           SSC_REC.REFHORN_TEMP,
	1           SSC_REC.ICAL_TEMP,
	1           SSC_REC.DIHEDRAL_TEMP,
	1           SSC_REC.BOLOMETER_TEMP,
	1           SSC_REC.BOL_CMD_BIAS )

6070	FORMAT (/10x,'Data Quality:        ',I6
	1       /10x,'Attitude Quality:    ',I6
	1       /10x,'External Calibrator Temperature: ',G14.5
	1       /10x,'Sky Horn Temperature:            ',G14.5
	1       /10x,'Reference Horn Temperature:      ',G14.5
	1       /10x,'Internal Calibrator Temperature: ',G14.5
	1       /10x,'Dihedral Temperature:            ',G14.5
	1       /10x,'Bolometer Temperature:           ',G14.5
	1       /10x,'Bolometer Bias Status:           ',G14.5 )
                                    
	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	RETURN
	END
