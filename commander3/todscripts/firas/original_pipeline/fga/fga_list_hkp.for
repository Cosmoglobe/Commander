	SUBROUTINE FGA_LIST_HKP (WRITE_LUN, HKP_REC, NREC, NMJFR,
	1		LINE, CAT_ENTRY)

C------------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FGA_LIST_HKP
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given FIRAS
C	  Housekeeping record.
C
C	AUTHOR:
C	  E.FUNG
C	  GSFC
C	  October 31, 1985
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_HKP (WRITE_LUN, HKP_REC, NREC, NMJFR, LINE, CAT_ENTRY)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file
C
C	  HKP_REC		RECORD	FIRAS Housekeeping record.
C
C	  NREC			I*2	The number of records processed in this run.
C
C	  NMJFR			I*2	1 = Just list the first major frame
C					2 = List both major frames
C
C         IPDU_AD_1             I*2     1 = power_a_status first major frame 
C         IPDU_AD_2             I*2     2 = power_b_status first major frame 
C         IPDU_AD_3             I*2     3 = power_a_status second major frame 
C         IPDU_AD_4             I*2     4 = power_b_status second major frame 
C
C	  LINE			C*80	Comment line for listing
C
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
C	  NONE
C
C	SUBROUTINES CALLED:
C	  LIB$MOVC3
C
C------------------------------------------------------------------------
C   Changes:
C	Add LMAC fields  (SPR 2634)  F. Shuman, STX, 1988 Oct 20.
C	Add telemetry quality field  (SPR 5594)  R. Kummerer, STX, 1990 Jan 9
C       Add fields HSKP1_TLM_FMT, HSKP2_TLM_HSKP2_TLM_FMT, POWER_A_STATUS,
C           POWER_B_STATUS (SPR 4234) N. Gonzales, STX, 7-12-90       
C       Corrected the wrong labels for channel temperature, Ref. SPR 9220.
C           N. Gonzales,Hughes/STX, November 1, 1991.
C------------------------------------------------------------------------

	IMPLICIT	NONE

	CHARACTER	*80	LINE

	INTEGER		*2	I
	INTEGER		*2	J
	INTEGER		*2	K
	INTEGER		*2	L
	INTEGER		*2	M
	INTEGER		*2	NREC
	INTEGER		*2	NMJFR
	INTEGER		*2      IPDU_AD_1
	INTEGER		*2      IPDU_AD_2
     	INTEGER		*2      IPDU_AD_3
      	INTEGER		*2      IPDU_AD_4
	INTEGER		*2      ANAL_DIG(4)
	INTEGER		*2      INT_CONV(4)
	INTEGER		*2      BIAS_PRE(4)
	INTEGER		*2      MTM_XCAL(4)
	INTEGER		*4	WRITE_LUN
	INTEGER		*4	CAT_ENTRY
	INTEGER		*4	STATUS

	INTEGER		*2	UHSK_LVDT_STAT_A
	INTEGER		*2	UHSK_LVDT_STAT_B
	INTEGER		*2	USBI(32)

	DICTIONARY 'NFS_HKP'
	RECORD /NFS_HKP/ HKP_REC

	EXTERNAL		FGA_WRITE_ERR

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6

	WRITE (WRITE_LUN, 6000, IOSTAT=STATUS)	NREC, CAT_ENTRY, LINE

6000	FORMAT ( '1','*****  FIRAS Housekeeping Record Formatted Listing  *****', 3x,
	1	'Record: ', i6, 10x, 'CT catalog entry: ', i7/3x, 
	1	'Comment line: ', a )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C	
C Cobetrieve standard header
C
	WRITE (WRITE_LUN, 6001, IOSTAT=STATUS)
	1			HKP_REC.CT_HEAD.GMT(1:2),
	1			HKP_REC.CT_HEAD.GMT(3:5),
	1			HKP_REC.CT_HEAD.GMT(6:7),
	1			HKP_REC.CT_HEAD.GMT(8:9),
	1			HKP_REC.CT_HEAD.GMT(10:11),
	1			HKP_REC.CT_HEAD.GMT(12:14),
	1			HKP_REC.CT_HEAD.TIME(2), 
	1			HKP_REC.CT_HEAD.TIME(1), 
	1			HKP_REC.CT_HEAD.SPACE_TIME,
	1			HKP_REC.CT_HEAD.ORBIT,
	1                       HKP_REC.CT_HEAD.HSKP1_TLM_FMT,
	1                       HKP_REC.CT_HEAD.HSKP2_TLM_FMT
6001	FORMAT (/' ASCII GMT:                        ', 2X,
	1	A2, '-', A3, '-', 3(A2, '-'), A3, t65,
	1	' GMT in binary (high order first): ', 3X,
	1	2(Z8.8,2X)
	1	/' S/C time in PB4 (HEX):            ', 4X,
	1	6(Z2.2,1X), t65,
	1	' Orbit Number:                     ',10X,
	1	I11, / 
	1       ' TLM FMT Major Frame 1: ',10X,I4 /
	1       ' TLM FMT Major Frame 2: ',10X,I4 ) 


	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C The 2 Major Frames
C
	DO I=1,NMJFR

	  IF (I .EQ. 1) THEN
	    WRITE (WRITE_LUN, 6100, IOSTAT=STATUS)
	  ELSE
	    WRITE (WRITE_LUN, 6200, IOSTAT=STATUS)
	  ENDIF

6100	  FORMAT (/' >>> First Major Frame in Record <<<' /)
6200	  FORMAT ('1'//' >>> Second Major Frame in Record <<<' /)

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

C
C Telemetry Quality.
C
	  WRITE (WRITE_LUN, 6008, IOSTAT=STATUS)
	1		HKP_REC.MJ_FRM(I).TLM_QUAL_MAJ_FRM

6008	  FORMAT (' Telemetry Quality:', I7, /)

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

C IPDU Power status bits for side A & B.
C
        IPDU_AD_1 = HKP_REC.MJ_FRM(1).POWER_A_STATUS
        IPDU_AD_2 = HKP_REC.MJ_FRM(1).POWER_B_STATUS
        IPDU_AD_3 = HKP_REC.MJ_FRM(2).POWER_A_STATUS
        IPDU_AD_4 = HKP_REC.MJ_FRM(2).POWER_B_STATUS

c       Bits 0:1  IPDU Analog & Digital Converters

        CALL mvbits (ipdu_ad_1,0,2,anal_dig(1),0)
        CALL mvbits (ipdu_ad_2,0,2,anal_dig(2),0)
        CALL mvbits (ipdu_ad_3,0,2,anal_dig(3),0)
        CALL mvbits (ipdu_ad_4,0,2,anal_dig(4),0)

c       Bits 2:3  IPDU Internal Converters

        CALL mvbits (ipdu_ad_1,2,2,int_conv(1),0)
        CALL mvbits (ipdu_ad_2,2,2,int_conv(2),0)
        CALL mvbits (ipdu_ad_3,2,2,int_conv(3),0)
        CALL mvbits (ipdu_ad_4,2,2,int_conv(4),0)

c       Bits 4:5  IPDU Bias Pre-Regulators

        CALL mvbits (ipdu_ad_1,4,2,bias_pre(1),0)
        CALL mvbits (ipdu_ad_2,4,2,bias_pre(2),0)
        CALL mvbits (ipdu_ad_3,4,2,bias_pre(3),0)
        CALL mvbits (ipdu_ad_4,4,2,bias_pre(4),0)

c       Bits 6:7  IPDU MTM & Xcal Motors

        CALL mvbits (ipdu_ad_1,6,2,mtm_xcal(1),0)
        CALL mvbits (ipdu_ad_2,6,2,mtm_xcal(2),0)
        CALL mvbits (ipdu_ad_3,6,2,mtm_xcal(3),0)
        CALL mvbits (ipdu_ad_4,6,2,mtm_xcal(4),0)

	  WRITE (WRITE_LUN, 6666, IOSTAT=STATUS)
	1        anal_dig(1),
	1        anal_dig(2),
	1        anal_dig(3),
	1        anal_dig(4),
	1        int_conv(1),
	1        int_conv(2),
	1        int_conv(3),
	1        int_conv(4),
	1        bias_pre(1),
	1        bias_pre(2),
	1        bias_pre(3),
	1        bias_pre(4),
	1        mtm_xcal(1),
	1        mtm_xcal(2),
	1        mtm_xcal(3),
	1        mtm_xcal(4)

6666	  FORMAT (/' ***** IPDU Power Status Bits ***** ' //
	1	' IPDU Anal & Dig_Conv MJFRM_1 side A: ', I4 / 
	1	' IPDU Anal & Dig_Conv MJFRM_1 side B: ', I4 /
	1	' IPDU Anal & Dig_Conv MJFRM_2 side A: ', I4 / 
	1	' IPDU Anal & Dig_Conv MJFRM_2 side B: ', I4 /
	1	' IPDU Internal   Conv MJFRM_1 side A: ', I4 / 
	1	' IPDU Internal   Conv MJFRM_1 side B: ', I4 /
	1	' IPDU Internal   Conv MJFRM_2 side A: ', I4 / 
	1	' IPDU Internal   Conv MJFRM_2 side B: ', I4 /
	1	' IPDU Bias Pre-Reg    MJFRM_1 side A: ', I4 /
 	1	' IPDU Bias Pre-Reg    MJFRM_1 side B: ', I4 /
	1	' IPDU Bias Pre-Reg    MJFRM_2 side A: ', I4 /
 	1	' IPDU Bias Pre-Reg    MJFRM_2 side B: ', I4 /
	1	' IPDU MTM & XCAL MTR  MJFRM_1 side A: ', I4 /
 	1	' IPDU MTM & XCAL MTR  MJFRM_1 side B: ', I4 /
	1	' IPDU MTM & XCAL MTR  MJFRM_2 side A: ', I4 /
 	1	' IPDU MTM & XCAL MTR  MJFRM_2 side B: ', I4 /)

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

C
C Status Monitor.
C
	  WRITE (WRITE_LUN, 6011, IOSTAT=STATUS)

6011	  FORMAT (' ***** Status Monitor Command Words ',
	1	'(all in HEX) ***** ' /)

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

	  IF (BTEST(HKP_REC.FRAME(I).HSKP_HEAD.STAT_MONITOR_CMD(1), 14)) THEN
	    WRITE (WRITE_LUN, 6003, IOSTAT=STATUS)
	1	(HKP_REC.FRAME(I).HSKP_HEAD.STAT_MONITOR_CMD(J), J = 1, 4)
	  ELSE
	    WRITE (WRITE_LUN, 6002, IOSTAT=STATUS)
	1	(HKP_REC.FRAME(I).HSKP_HEAD.STAT_MONITOR_CMD(J), J = 1, 4)
	  ENDIF

6002	  FORMAT ('    +++++  Side A in FIRAS Major Frame 1  +++++ '//
	1	' COMMAND 0, SIDE A: ', Z4.4, 5x,
	1	' COMMAND 1, SIDE A: ', Z4.4, 5x,
	1	' COMMAND 2, SIDE A: ', Z4.4, 5x,
	1	' COMMAND 3, SIDE A: ', Z4.4)

6003	  FORMAT ('    +++++  Side A in FIRAS Major Frame 2  +++++  '//
	1	' COMMAND 4, SIDE A: ', Z4.4, 5x,
	1	' COMMAND 5, SIDE A: ', Z4.4, 5x,
	1	' COMMAND 6, SIDE A: ', Z4.4, 5x,
	1	' COMMAND 7, SIDE A: ', Z4.4 )

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF


	  IF (BTEST(HKP_REC.FRAME(I).HSKP_HEAD.STAT_MONITOR_CMD(5), 14)) THEN
	    WRITE (WRITE_LUN, 6103, IOSTAT=STATUS)
	1	(HKP_REC.FRAME(I).HSKP_HEAD.STAT_MONITOR_CMD(J), J = 5, 8)
	  ELSE
	    WRITE (WRITE_LUN, 6102, IOSTAT=STATUS)
	1	(HKP_REC.FRAME(I).HSKP_HEAD.STAT_MONITOR_CMD(J), J = 5, 8)
	  ENDIF

6102	  FORMAT (/'    +++++  Side B in FIRAS Major Frame 1  +++++ '//
	1	' COMMAND 0, SIDE B: ', Z4.4, 5x,
	1	' COMMAND 1, SIDE B: ', Z4.4, 5x,
	1	' COMMAND 2, SIDE B: ', Z4.4, 5x,
	1	' COMMAND 3, SIDE B: ', Z4.4 )

6103	  FORMAT (/'    +++++  Side B in FIRAS Major Frame 2  +++++  '//
	1	' COMMAND 4, SIDE B: ', Z4.4, 5x,
	1	' COMMAND 5, SIDE B: ', Z4.4, 5x,
	1	' COMMAND 6, SIDE B: ', Z4.4, 5x,
	1	' COMMAND 7, SIDE B: ', Z4.4 )

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

C
C IPDU statuses
C
	  WRITE (WRITE_LUN, 6004, IOSTAT=STATUS)
	1		(HKP_REC.FRAME(I).HSKP_HEAD.IPDU_STAT(J),
	1		HKP_REC.FRAME(I).HSKP_HEAD.IPDU_STAT(J+4),  J=1, 4)

6004	FORMAT (/' ***** IPDU Statuses (all in HEX) ***** ' //
	1	' IPDU Power Relay 4 Pole A1:     ', Z2.2, t65,
	1	' IPDU Power Relay 2 Pole A3:     ', Z2.2 /
	1	' IPDU Power Relay 4 Pole B1:     ', Z2.2, t65,
	1	' IPDU Power Relay 2 Pole B3:     ', Z2.2 /
	1	' IPDU Power Relay 4 Pole A2:     ', Z2.2, t65,
	1	' IPDU Other Status A4:           ', Z2.2 /
	1	' IPDU Power Relay 4 Pole B2:     ', Z2.2 , t65
	1	' IPDU Other Status B4:           ', Z2.2 )

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

C
C Dwell, MTM status and position, microprocess bus readouts,
C bolometer bias statuses
C
	  UHSK_LVDT_STAT_A = ZEXT (HKP_REC.FRAME(I).HSKP_HEAD.LVDT_STAT(1))
	  UHSK_LVDT_STAT_B = ZEXT (HKP_REC.FRAME(I).HSKP_HEAD.LVDT_STAT(2))

	  DO L=1,5
	    USBI(L) = ZEXT (HKP_REC.MJ_FRM(I).NTCH_FLT_A(L))
	    USBI(L+5) = ZEXT (HKP_REC.MJ_FRM(I).NTCH_FLT_B(L))
	  END DO

	  WRITE (WRITE_LUN, 6005, IOSTAT=STATUS)
	1		HKP_REC.FRAME(I).HSKP_HEAD.DWELL_STAT(1),
	1		HKP_REC.FRAME(I).HSKP_HEAD.DWELL_STAT(2),
	1		UHSK_LVDT_STAT_A,
	1		UHSK_LVDT_STAT_B,
	1		(HKP_REC.FRAME(I).HSKP_HEAD.U_PROC_STAT(J), J=1,4),
	1		(HKP_REC.FRAME(I).HSKP_HEAD.BOL_CMD_BIAS(K), K=1,4),
	1		(USBI(L), USBI(L+5), L=1,5)

6005	FORMAT (/' ***** Other Statuses ***** ' /
	1	/' Dwell Status, Side A (HEX):         ',2x, Z2.2
	1	, t65, ' Dwell Status, Side B (HEX):         ',2x, Z2.2
	1	/' LVDT Status, Side A:   ', I5
	1	, t65, ' LVDT Status, Side B:   ', I5
	1	/' Microprocessor bus readout, RH (HEX): ', Z2.2
	1	, t65, ' Microprocessor bus readout, RL (HEX): ', Z2.2
	1	/' Microprocessor bus readout, LH (HEX): ', Z2.2
	1	, t65, ' Microprocessor bus readout, LL (HEX): ', Z2.2
	1	/' Bolometer bias status, RH (HEX):    ', 2X, Z2.2
	1	, t65, ' Bolometer bias status, RL (HEX):    ', 2X, Z2.2
	1	/' Bolometer bias status, LH (HEX):    ', 2X, Z2.2, t65,
	1	' Bolometer bias status, LL (HEX):    ', 2X, Z2.2
	1	/' Notch Filter 1, Side A (HEX):         ',2x, Z2.2
	1	, t65, ' Notch Filter 1, Side B (HEX):         ',2x, Z2.2
	1	/' Notch Filter 2, Side A (HEX):         ',2x, Z2.2
	1	, t65, ' Notch Filter 2, Side B (HEX):         ',2x, Z2.2
	1	/' Notch Filter 3, Side A (HEX):         ',2x, Z2.2
	1	, t65, ' Notch Filter 3, Side B (HEX):         ',2x, Z2.2
	1	/' Notch Filter 4, Side A (HEX):         ',2x, Z2.2
	1	, t65, ' Notch Filter 4, Side B (HEX):         ',2x, Z2.2
	1	/' Notch Filter 5, Side A (HEX):         ',2x, Z2.2
	1	, t65, ' Notch Filter 5, Side B (HEX):         ',2x, Z2.2)

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

C
C GRT's
C
	  WRITE (WRITE_LUN, 6006, IOSTAT=STATUS)

6006	  FORMAT (/' ***** GRTs:  Low Current, A-side *****' , T65,
	1	' ***** GRTs:  High Current, A-side *****' /)

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF


	  WRITE (WRITE_LUN, 6007, IOSTAT=STATUS)
	1	(HKP_REC.FRAME(I).TEMPS.SIDE_AMP(1).GRT(J),
	1	HKP_REC.FRAME(I).TEMPS.SIDE_AMP(2).GRT(J), J=1, 16)

6007	  FORMAT (' External Calibrator:          ', I7, T65,
	1	' External Calibrator:          ', I7
	1	/ ' Sky Horn:                     ', I7, T65,
	1	 ' Sky Horn:                     ', I7
	1	/ ' Reference Horn:               ', I7, T65,
	1	 ' Reference Horn:               ', I7
	1	/ ' Internal Reference Source:    ', I7, T65,
	1	 ' Internal Reference Source:    ', I7
	1	/ ' Right Dihedral:               ', I7, T65,
	1	 ' Right Dihedral:               ', I7
	1	/ ' Bolometer Assembly, RH:       ', I7, T65,
	1	 ' Bolometer Assembly, RH:       ', I7
	1	/ ' Bolometer Assembly, RL:       ', I7, T65,
	1	 ' Bolometer Assembly, RL:       ', I7
	1	/ ' Bolometer Assembly, LH:       ', I7, T65,
	1	 ' Bolometer Assembly, LH:       ', I7
	1	/ ' Bolometer Assembly, LL:       ', I7, T65,
	1	 ' Bolometer Assembly, LL:       ', I7
	1	/ ' Right Mirror:                 ', I7, T65,
	1	 ' Right Mirror:                 ', I7
	1	/ ' Calibrator Resistor 1:        ', I7, T65,
	1	 ' Calibrator Resistor 1:        ', I7
	1	/ ' Calibrator Resistor 2:        ', I7, T65,
	1	 ' Calibrator Resistor 2:        ', I7
	1	/ ' Calibrator Resistor 3:        ', I7, T65,
	1	 ' Calibrator Resistor 3:        ', I7
	1	/ ' Calibrator Resistor 4:        ', I7, T65,
	1	 ' Calibrator Resistor 4:        ', I7
	1	/ ' XCal Segment 5:               ', I7, T65,
	1	 ' XCal Segment 5:               ', I7
	1	/ ' Right Collimator:             ', I7, T65,
	1	 ' Right Collimator:             ', I7 )

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF


	  WRITE (WRITE_LUN, 6009, IOSTAT=STATUS)

6009	  FORMAT (/' ***** GRTs:  Low Current, B-side *****' , T65,
	1	' ***** GRTs:  High Current, B-side *****' /)

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF



	  WRITE (WRITE_LUN, 6010, IOSTAT=STATUS)
	1	(HKP_REC.FRAME(I).TEMPS.SIDE_AMP(3).GRT(J),
	1	HKP_REC.FRAME(I).TEMPS.SIDE_AMP(4).GRT(J), J=1, 16)

6010	  FORMAT (' External Calibrator:          ', I7, T65,
	1	' External Calibrator:          ', I7
	1	/ ' Sky Horn:                     ', I7, T65,
	1	 ' Sky Horn:                     ', I7
	1	/ ' Reference Horn:               ', I7, T65,
	1	 ' Reference Horn:               ', I7
	1	/ ' Internal Reference Source:    ', I7, T65,
	1	 ' Internal Reference Source:    ', I7
	1	/ ' Left Dihedral:                ', I7, T65,
	1	 ' Left Dihedral:                ', I7
	1	/ ' Bolometer Assembly, RH:       ', I7, T65,
	1	 ' Bolometer Assembly, RH:       ', I7
	1	/ ' Bolometer Assembly, RL:       ', I7, T65,
	1	 ' Bolometer Assembly, RL:       ', I7
	1	/ ' Bolometer Assembly, LH:       ', I7, T65,
	1	 ' Bolometer Assembly, LH:       ', I7
	1	/ ' Bolometer Assembly, LL:       ', I7, T65,
	1	 ' Bolometer Assembly, LL:       ', I7
	1	/ ' Left Mirror:                  ', I7, T65,
	1	 ' Left Mirror:                  ', I7
	1	/ ' Calibrator Resistor 1:        ', I7, T65,
	1	 ' Calibrator Resistor 1:        ', I7
	1	/ ' Calibrator Resistor 2:        ', I7, T65,
	1	 ' Calibrator Resistor 2:        ', I7
	1	/ ' Calibrator Resistor 3:        ', I7, T65,
	1	 ' Calibrator Resistor 3:        ', I7
	1	/ ' Calibrator Resistor 4:        ', I7, T65,
	1	 ' Calibrator Resistor 4:        ', I7
	1	/ ' XCal Segment 6:               ', I7, T65,
	1	 ' XCal Segment 6:               ', I7
	1	/ ' Left Collimator (Not used):   ', I7, T65,
	1	 ' Left Collimator (Not used):   ', I7 )

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

C
C Other temperatures/currents
C
	  WRITE (WRITE_LUN, 6012, IOSTAT=STATUS)

6012	  FORMAT (/' ***** Other Temperatures/Currents *****' /)

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF


	  USBI(1) = ZEXT(HKP_REC.FRAME(I).TEMPS.IPDU(1))
	  USBI(2) = ZEXT(HKP_REC.FRAME(I).TEMPS.IPDU(2))
	  USBI(3) = ZEXT(HKP_REC.FRAME(I).TEMPS.DRIVE_BOX_A)
	  USBI(4) = ZEXT(HKP_REC.FRAME(I).TEMPS.DRIVE_BOX_B)
	  USBI(5) = ZEXT(HKP_REC.FRAME(I).TEMPS.CHAN_TEMP(1))
	  USBI(6) = ZEXT(HKP_REC.FRAME(I).TEMPS.CHAN_TEMP(2))
	  USBI(7) = ZEXT(HKP_REC.FRAME(I).TEMPS.CHAN_TEMP(3))
	  USBI(8) = ZEXT(HKP_REC.FRAME(I).TEMPS.CHAN_TEMP(4))
	  USBI(9) = ZEXT(HKP_REC.FRAME(I).TEMPS.STAT_MON(1))
	  USBI(10) = ZEXT(HKP_REC.FRAME(I).TEMPS.STAT_MON(2))
	  USBI(11) = ZEXT(HKP_REC.FRAME(I).TEMPS.HOT_SPOT_CURRENT(1))
	  USBI(12) = ZEXT(HKP_REC.FRAME(I).TEMPS.HOT_SPOT_CURRENT(2))
	  USBI(13) = ZEXT(HKP_REC.FRAME(I).TEMPS.OPTICAL_PREAMP)
	  USBI(14) = ZEXT(HKP_REC.FRAME(I).TEMPS.CHAN_PRE_AMP)
	  USBI(15) = ZEXT(HKP_REC.MJ_FRM(I).LMAC_ANALOG_TEMP)
	  USBI(16) = ZEXT(HKP_REC.MJ_FRM(I).LMAC_DIGITAL_TEMP)


	  WRITE (WRITE_LUN, 6013, IOSTAT=STATUS)  (USBI(J),J=1,16)

6013	FORMAT (/' IPDU temperature - Side A:       ', I7, T65,
	1	' IPDU temperature - Side B:       ', I7
	1	/' Drive Box temperature - Side A:  ', I7, T65,
	1	' Drive Box temperature - Side B:  ', I7
	1	/' Channel temperature - RH:        ', I7, T65,
	1	' Channel temperature - RL:        ', I7
	1	/' Channel temperature - LH:        ', I7, T65,
	1	' Channel temperature - LL:        ', I7
	1	/' Status Monitor temp. - Side A:   ', I7, T65,
	1	' Status Monitor temp. - Side B:   ', I7
	1	/' Hot Spot Heater Current-Side A:  ', I7, T65,
	1	' Hot Spot Heater Current-Side B:  ', I7
	1	/' Optical Pre-Amp temperature:     ', I7, T65,
	1	' Channel Pre-Amp temperature:     ', I7
	1	/' LMAC Analog temperature:         ', I7, T65,
	1	' LMAC Digital temperature:        ', I7 /)

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF


C
C IPDU Voltages
C
	  DO	J = 1, 20
	    USBI(J) = ZEXT(HKP_REC.FRAME(I).V_AND_I.IPDU_VOLT(J))
	  ENDDO

	  WRITE (WRITE_LUN, 6014, IOSTAT=STATUS) (USBI(J),J=1,20)

6014	  FORMAT ( /' ***** IPDU Voltages *****' /
	1	/' Digital Converter, -15V, Side A: ', I7, T65,
	1	' Digital Converter, -15V, Side B: ', I7
	1	/' Digital Converter, +15V, Side A: ', I7, T65,
	1	' Digital Converter, +15V, Side B: ', I7
	1	/' Digital Converter, +5V, Side A:  ', I7, T65,
	1	' Digital Converter, +5V, Side B:  ', I7
	1	/' Analog Converter, +15V, Side A:  ', I7, T65,
	1	' Analog Converter, +15V, Side B:  ', I7
	1	/' Analog Converter, -15V, Side A:  ', I7, T65,
	1	' Analog Converter, -15V, Side B:  ', I7
	1	/' Bias Pre-Reg. +25V, Side A:      ', I7, T65,
	1	' Bias Pre-Reg. +25V, Side B:      ', I7
	1	/' Int. P.S., +28V, Side A:         ', I7, T65,
	1	' Int. P.S., +28V, Side B:         ', I7
	1	/' Int. P.S., +15V, Side A:         ', I7, T65,
	1	' Int. P.S., +15V, Side B:         ', I7
	1	/' Int. P.S., -15V, Side A:         ', I7, T65,
	1	' Int. P.S., -15V, Side B:         ', I7
	1	/' Int. P.S., +5V, Side A:          ', I7, T65,
	1	' Int. P.S., +5V, Side B:          ', I7 )

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

C
C IPDU Currents
C
	  DO	J = 1, 12
	    USBI(J) = ZEXT(HKP_REC.FRAME(I).V_AND_I.IPDU_AMP(J))
	  ENDDO

	  WRITE (WRITE_LUN, 6015, IOSTAT=STATUS) (USBI(J),J=1,12)

6015	  FORMAT ( /' ***** IPDU Currents *****' /
	1	/' Bias Pre-Reg. Current, Side A:   ', I7    ,T65,
	1	' Bias Pre-Reg. Current, Side B:   ', I7
	1	/' Analog Conv. Current, Side A:    ', I7, T65,
	1	' Analog Conv. Current, Side B:    ', I7
	1	/' Digital Conv. Current, Side A:   ', I7, T65,
	1	' Digital Conv. Current, Side B:   ', I7
	1	/' Constant Current RH:             ', I7, T65,
	1	' Constant Current LH:             ', I7
	1	/' Constant Current RL:             ', I7, T65,
	1	' Constant Current LL:             ', I7
	1	/' Int. Conv. Current, Side A:      ', I7, T65,
	1	' Int. Conv. Current, Side B:      ', I7 /)

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

C
C MTM/Cal Motor and Bolometer Bias Voltages
C
	  USBI(1) = ZEXT (HKP_REC.FRAME(I).V_AND_I.MTM_CAL_MOTOR(1))
	  USBI(2) = ZEXT (HKP_REC.FRAME(I).V_AND_I.MTM_CAL_MOTOR(2))

	  DO J=1,4
	    USBI(J+2) = ZEXT(HKP_REC.FRAME(I).V_AND_I.BOL_BIAS_VOLT(J))
	  ENDDO

	  WRITE (WRITE_LUN, 6016, IOSTAT=STATUS) (USBI(J),J=1,6)

6016	  FORMAT (/' ***** MTM/Cal Motor and Bolometer Bias Voltages *****'/
	1	/' MTM/Cal Motor, Side A:           ', I7, T65,
	1	' MTM/Cal Motor, Side B:           ', I7
	1	/' Bolometer Bias Voltage, RH:      ', I7, T65,
	1	' Bolometer Bias Voltage, RL:      ', I7
	1	/' Bolometer Bias Voltage, LH:      ', I7, T65,
	1	' Bolometer Bias Voltage, LL:      ', I7 )

	  IF (STATUS .NE. 0) THEN
	    CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	  END IF

	END DO

C
C	Now do the trailer: GMT of 2nd major frame, and SIS status
C
	WRITE (WRITE_LUN, 6017, IOSTAT=STATUS)
	1			HKP_REC.HSKP_TAIL.GMT_MJF2(1:2),
	1			HKP_REC.HSKP_TAIL.GMT_MJF2(3:5),
	1			HKP_REC.HSKP_TAIL.GMT_MJF2(6:7),
	1			HKP_REC.HSKP_TAIL.GMT_MJF2(8:9),
	1			HKP_REC.HSKP_TAIL.GMT_MJF2(10:11),
	1			HKP_REC.HSKP_TAIL.GMT_MJF2(12:14)

6017	FORMAT (/ ' ***** Trailer Information *****'//
	1	' GMT of 2nd Major frame: ', 
	1	A2, '-', A3, '-', 3(A2, '-'), A3)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	RETURN
	END
