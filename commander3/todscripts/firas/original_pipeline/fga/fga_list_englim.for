	Integer*4 Function FGA_LIST_ENGLIM (WRITE_LUN, ENGLIM_REC, ENGLIM_FILE,
	1				    CLINE, NREC)

C---------------------------------------------------------------------------
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given engineering
C	  limits record.
C
C	AUTHOR:
C	  QUOC CHUNG
C	  STX
C	  AUG 22, 1989
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_ENGLIM (WRITE_LUN, ENGLIM_REC, ENGLIM_FILE, CLINE, NREC)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file
C	  ENGLIM_REC         		ENGLIM record.
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
C
C	SPR 9161.  October 16, 1991, S. Alexander, STX.
C		   Fix error on B_HI_GRT output; add LMAC output; remove
C		   temperature controller and engineering status writes,
C		   as these quantities are never put in the englims file.
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
        INTEGER         *4      MJR_FRAME
        INTEGER         *4      SUCCESS/1/,ERROR/2/

	REAL            *4      A_LO_TEMP_GRT(16)
	REAL            *4      A_HI_TEMP_GRT(16)
	REAL            *4      B_LO_TEMP_GRT(16)
	REAL            *4      B_HI_TEMP_GRT(16)
	REAL            *4      IPDU_TMP(2)
	REAL            *4      CNA_TMP(4)
	REAL            *4      DBX_TMP(2)
	REAL            *4      STAT_MON_TMP(2)
	REAL            *4      H_SPOT(2)
	REAL            *4      MTM_CAL_MTOR(2)
	REAL            *4      MTM_POSIT(2)
	REAL            *4      BMTR_VLT(4)
	REAL	        *4      IPDU_VLT(20)
	REAL            *4      IPDU_AMPS(12)
	REAL		*4	LMAC(2)

	CHARACTER       *64     ENGLIM_FILE		
	CHARACTER	*80	CLINE
	CHARACTER	*1	LINE(130)/ 130 * '-'/
        CHARACTER       *1      CBITS16(16)
        CHARACTER       *1      CBITS32(32)
        CHARACTER       *14     ASCIGMT,MJRFRAME2_GMT
        CHARACTER       *14     GMT_START,GMT_STOP
	CHARACTER       *11     TEXT(4)/'RED LOW','YELLOW LOW',
	1                               'YELLOW HIGH','RED HIGH'/

        LOGICAL         *1	FIRST/.TRUE./

	DICTIONARY	'FEX_ENGLIM'
	RECORD /FEX_ENGLIM/      ENGLIM_REC

C
C Begin the dump.
C
        FGA_LIST_ENGLIM = SUCCESS

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6
	IF ( FIRST) THEN
	  WRITE (WRITE_LUN, 6000, IOSTAT=IOSTATUS) ENGLIM_FILE,CLINE,LINE
          FIRST = .FALSE.
        ENDIF
 6000   FORMAT(	// 1X, 'FEX_ENGLIM FILE : ',A/,X,A80/,130A/)

	IF (IOSTATUS .NE. 0) THEN
	  PRINT *,' WRITE ERROR!!! UNIT = ', WRITE_LUN 
          FGA_LIST_ENGLIM = ERROR          
	END IF

	WRITE ( WRITE_LUN,600) NREC
  600   FORMAT(/,T50,' REC # : ',I4/)
C	
C Header information.
C
	WRITE (WRITE_LUN, 6001, IOSTAT=IOSTATUS)
	1              (ENGLIM_REC.CT_HEAD.GMT(I:I),I=1,14),
	1		ENGLIM_REC.CT_HEAD.TIME(2),
	1		ENGLIM_REC.CT_HEAD.TIME(1),
	1		ENGLIM_REC.CT_HEAD.SPACE_TIME,
	1		ENGLIM_REC.CT_HEAD.MJR_FRM_NO,
	1		ENGLIM_REC.CT_HEAD.ORBIT,
	1		ENGLIM_REC.CT_HEAD.DATASET_ID,
	1               LINE

6001	FORMAT (/' <Start GMT> : ', 2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, 7X,
	1	 ' <Binary TIME> : ', Z8.8, 1X, Z8.8, 7X, /,
	1        ' <PB5> : ', 6(Z2.2, 1X) /
	1        ' Major Frame Number : ', I6, 8x,
	1        ' Orbit Number : ', I12, 12x,
	1	 ' Data Set Id  : ', I10/,130A/)
C
C Dump the FEX_ENGLIM record.
C Loop through four limits: Red Low, Yellow Low, Yellow High, Red High.
C
        DO I=1,4

		WRITE (WRITE_LUN, 6002, IOSTAT=IOSTATUS) TEXT(I),
	1               LINE

6002        	FORMAT  (/1X,'FEX_ENGLIM Data Record for ',A,' : '/,
	1         130A/)

		DO J = 1 , 16
                  A_LO_TEMP_GRT(J) = ENGLIM_REC.LIM(I).EN_ANALOG.A_LO_GRT(J)
                  A_HI_TEMP_GRT(J) = ENGLIM_REC.LIM(I).EN_ANALOG.A_HI_GRT(J)
                  B_LO_TEMP_GRT(J) = ENGLIM_REC.LIM(I).EN_ANALOG.B_LO_GRT(J)
                  B_HI_TEMP_GRT(J) = ENGLIM_REC.LIM(I).EN_ANALOG.B_HI_GRT(J)
		ENDDO

		WRITE(WRITE_LUN,6004) 
	4    		A_LO_TEMP_GRT,
    	4    		A_HI_TEMP_GRT,
	4    		B_LO_TEMP_GRT,
	4    		B_HI_TEMP_GRT,
	4    		LINE
	      
 6004		FORMAT(' 16 A_LO_GRT : ',8F10.2/,15X,8F10.2/,
	1      	' 16 A_HI_GRT : ',8F10.2/,15X,8F10.2/,
	1      	' 16 B_LO_GRT : ',8F10.2/,15X,8F10.2/,
	1      	' 16 B_HI_GRT : ',8F10.2/,15X,8F10.2/,
	1      	130A )

		DO K = 1 , 20
	          IPDU_VLT(K) = ENGLIM_REC.LIM(I).EN_ANALOG.IPDU_VOLT(K)
                ENDDO

		DO M = 1 , 12
	          IPDU_AMPS(M) = ENGLIM_REC.LIM(I).EN_ANALOG.IPDU_AMP(M)
		ENDDO

		DO J = 1 , 4
	          CNA_TMP(J) = ENGLIM_REC.LIM(I).EN_ANALOG.CNA_TEMP(J)
                  BMTR_VLT(J) = ENGLIM_REC.LIM(I).EN_ANALOG.BOL_VOLT(J)
		ENDDO

		DO K = 1 , 2
		  IPDU_TMP(K) = ENGLIM_REC.LIM(I).EN_ANALOG.IPDU_TEMP(K)
                  DBX_TMP(K) = ENGLIM_REC.LIM(I).EN_ANALOG.DBX_TEMP(K)
                  STAT_MON_TMP(K) = ENGLIM_REC.LIM(I).EN_ANALOG.STAT_MON_TEMP(K)
                  H_SPOT(K) = ENGLIM_REC.LIM(I).EN_ANALOG.HOT_SPOT(K)
                  MTM_CAL_MTOR(K) = ENGLIM_REC.LIM(I).EN_ANALOG.MTM_CAL_MTR(K)
                  MTM_POSIT(K) = ENGLIM_REC.LIM(I).EN_ANALOG.MTM_POS(K)
	        ENDDO

	        LMAC(1) = ENGLIM_REC.LIM(I).EN_TAIL.LMAC_ANALOG_TEMP
	        LMAC(2) = ENGLIM_REC.LIM(I).EN_TAIL.LMAC_DIGITAL_TEMP

		WRITE(WRITE_LUN,6005)
	5    	IPDU_TMP,CNA_TMP,DBX_TMP,STAT_MON_TMP,
	5    	ENGLIM_REC.LIM(I).EN_ANALOG.PAMP_CHAN,
	5    	ENGLIM_REC.LIM(I).EN_ANALOG.PAMP_OP,
	5    	H_SPOT,MTM_CAL_MTOR,MTM_POSIT,BMTR_VLT,LMAC,
	5    	IPDU_VLT,IPDU_AMPS,
	5	LINE

 6005		FORMAT(/,' < ENGINEERING ANALOG > '/,
	5 	' IPDU TEMPERATURE A , B  : ',2F10.2/,
	5 	' CHANNEL TEMPERATURE     : ',4F10.2/,
	5 	' DRIVE BOX TEMPERATURE   : ',2F10.2/,
	5 	' STATUS MONITOR TEMP.    : ',2F10.2/,
	5 	' CHANNEL PRE-AMP         : ',F10.2/,
	5 	' OPTICAL PRE-AMP         : ',F10.2/,
	5 	' HOT SPOT HEATER A, B    : ',2F10.2/,
	5 	' MTM/CAL MOTORS          : ',2F10.2/,
	5 	' MTM POSITION            : ',2F10.2/,
	5 	' BOLOMETER VOLTAGES      : ',4F10.3/,
	5	' LMAC TEMPERATURE        : ',2F10.2/,
	5 	' IPDU VOLTAGES           : ',10F10.2/,27X,10F10.2/,
	5 	' IPDU CURRENTS           : ',6F10.4/,27X,6F10.4/,130A/)

	ENDDO

        IF (FGA_LIST_ENGLIM .NE. SUCCESS) THEN
          PRINT *, ' FGA_LIST_ENGLIM Fail !!!!'
        END IF

        RETURN
	END
