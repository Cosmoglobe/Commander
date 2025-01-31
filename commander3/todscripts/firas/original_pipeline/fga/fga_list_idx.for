	SUBROUTINE FGA_LIST_IDX (WRITE_LUN, IDX_REC, NREC, LINE,
	1	CAT_ENTRY)

C------------------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FGA_LIST_IDX
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given FIRAS
C	  Index record.
C
C	AUTHOR:
C	  E.FUNG
C	  GSFC
C	  November 21, 1985
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_IDX (WRITE_LUN, IDX_REC, NREC, LINE, cat_entry)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file
C
C	  IDX_REC		RECORD	FIRAS Index record.
C
C	  NREC			I*2	The number of records processed in this run.
C
C	  LINE			C*80	Comment line for listing
C
C	  CAT_ENTRY		I*4	CT catalog entry
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
C------------------------------------------------------------------------------
C  Changes:
C   Add IDX_TAIL Structure for Tolerance values (SPR 4234)
C              N. Gonzales, STX, 10-July-90     
C------------------------------------------------------------------------------

	IMPLICIT	NONE

	CHARACTER	*80	LINE
	CHARACTER	*2	CHAN(4)/'RH', 'RL', 'LH', 'LL'/
	INTEGER		*2	NREC
	INTEGER		*4	WRITE_LUN
	INTEGER		*4	CAT_ENTRY
	INTEGER		*4	STATUS
	INTEGER		*2	I
	INTEGER		*2	J
	INTEGER		*2	ST
	INTEGER		*2	USBI(20)

	DICTIONARY 'FDQ_IDX'
	RECORD /FDQ_IDX/ IDX_REC

	EXTERNAL		FGA_WRITE_ERR

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6

	WRITE (WRITE_LUN, 6000, IOSTAT=STATUS)	LINE, NREC, CAT_ENTRY

6000	FORMAT ( '1' /15X, '***** FIRAS Engineering Index Formatted Listing *****'  
	1	/3x, a / ' Record Count for this Run: ', I6, 
	1	10x, 'CT Catalog Entry: ', I7  /)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C	
C Cobetrieve standard header
C
	WRITE (WRITE_LUN, 6001, IOSTAT=STATUS)
	1			IDX_REC.HEADER.GMT_START_TIME_ASCII(1:2),
	1			IDX_REC.HEADER.GMT_START_TIME_ASCII(3:5),
	1			IDX_REC.HEADER.GMT_START_TIME_ASCII(6:7),
	1			IDX_REC.HEADER.GMT_START_TIME_ASCII(8:9),
	1			IDX_REC.HEADER.GMT_START_TIME_ASCII(10:11),
	1			IDX_REC.HEADER.GMT_START_TIME_ASCII(12:14),
	1			IDX_REC.HEADER.BINARY_START_TIME(2), 
	1			IDX_REC.HEADER.BINARY_START_TIME(1), 
	1			IDX_REC.HEADER.SPACECRAFT_START_TIME_PB5,
	1			IDX_REC.HEADER.GMT_END_TIME_ASCII(1:2),
	1			IDX_REC.HEADER.GMT_END_TIME_ASCII(3:5),
	1			IDX_REC.HEADER.GMT_END_TIME_ASCII(6:7),
	1			IDX_REC.HEADER.GMT_END_TIME_ASCII(8:9),
	1			IDX_REC.HEADER.GMT_END_TIME_ASCII(10:11),
	1			IDX_REC.HEADER.GMT_END_TIME_ASCII(12:14),
	1			IDX_REC.HEADER.BINARY_END_TIME(2), 
	1			IDX_REC.HEADER.BINARY_END_TIME(1), 
	1			IDX_REC.HEADER.SPACECRAFT_END_TIME_PB5

6001	FORMAT (3X, ' START GMT (ASCII):    ', 4X,
	1	A2, '-', A3, '-', 3(A2, '-'), A3, 7X, '(BINARY): ',
	1	2(Z8.8, 2X), 7X, '(PB5): ', 6(Z2.2,1X) / 3X,
	1	' END GMT (ASCII):      ', 4X,
	1	A2, '-', A3, '-', 3(A2, '-'), A3, 7X, '(BINARY): ',
	1	2(Z8.8, 2X), 7X, '(PB5): ', 6(Z2.2,1X) )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C XCAL position
C
	WRITE (WRITE_LUN, 6011, IOSTAT=STATUS)
	1	IDX_REC.XCAL.POS(1), IDX_REC.XCAL.POS(2)

6011	FORMAT (/3x, ' ***** External Calibrator Position ***** (0 = ',
	1	'Travel; 1 = Stow; 2 = Cal; 3 = Error): ' / 3x,
	1	' Side A: ', I3, t65, ' Side B: ', I3)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Channel-specific Stuff
C
	WRITE (WRITE_LUN, 6300, IOSTAT=STATUS)

6300	FORMAT (/3X, ' ***** Channel Specific Fields *****' )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	DO I=1,4
	  USBI(I) = ZEXT(IDX_REC.CHAN(I).BOL_READOUT_VOLT)
	END DO

	WRITE (WRITE_LUN, 6301, IOSTAT=STATUS)
	1			IDX_REC.CHAN(1).UP_SCI_MODE,
	1			IDX_REC.CHAN(2).UP_SCI_MODE,
	1			IDX_REC.CHAN(3).UP_SCI_MODE,
	1			IDX_REC.CHAN(4).UP_SCI_MODE,
	1			IDX_REC.CHAN(1).FAKEIT,
	1			IDX_REC.CHAN(2).FAKEIT,
	1			IDX_REC.CHAN(3).FAKEIT,
	1			IDX_REC.CHAN(4).FAKEIT,
	1			IDX_REC.CHAN(1).UP_ADDS_PER_GRP,
	1			IDX_REC.CHAN(2).UP_ADDS_PER_GRP,
	1			IDX_REC.CHAN(3).UP_ADDS_PER_GRP,
	1			IDX_REC.CHAN(4).UP_ADDS_PER_GRP,
	1			IDX_REC.CHAN(1).UP_SWEEPS_PER_IFG,
	1			IDX_REC.CHAN(2).UP_SWEEPS_PER_IFG,
	1			IDX_REC.CHAN(3).UP_SWEEPS_PER_IFG,
	1			IDX_REC.CHAN(4).UP_SWEEPS_PER_IFG,
	1			IDX_REC.CHAN(1).MTM_SPEED_AT_XMIT,
	1			IDX_REC.CHAN(2).MTM_SPEED_AT_XMIT,
	1			IDX_REC.CHAN(3).MTM_SPEED_AT_XMIT,
	1			IDX_REC.CHAN(4).MTM_SPEED_AT_XMIT,
	1			IDX_REC.CHAN(1).MTM_LENGTH_AT_XMIT,
	1			IDX_REC.CHAN(2).MTM_LENGTH_AT_XMIT,
	1			IDX_REC.CHAN(3).MTM_LENGTH_AT_XMIT,
	1			IDX_REC.CHAN(4).MTM_LENGTH_AT_XMIT,
	1			IDX_REC.CHAN(1).SCIENCE_GAIN,
	1			IDX_REC.CHAN(2).SCIENCE_GAIN,
	1			IDX_REC.CHAN(3).SCIENCE_GAIN,
	1			IDX_REC.CHAN(4).SCIENCE_GAIN,
	1			IDX_REC.CHAN(1).BOL_CMD_BIAS,
	1			IDX_REC.CHAN(2).BOL_CMD_BIAS,
	1			IDX_REC.CHAN(3).BOL_CMD_BIAS,
	1			IDX_REC.CHAN(4).BOL_CMD_BIAS,
	1			(USBI(J),J=1,4)

6301	FORMAT	(32x, '<< Channel RH >> ', 5X,
	1	'<< Channel RL >>', 5X,
	1	'<< Channel LH >>', 5X,
	1	'<< Channel LL >>' /
	1	' Microprocessor Science Mode:   ', 4(I10, 11X) /
	1	' Fake-it:                       ', 4(I10,11X) /
	1	' Microprocessor Adds/Group:     ', 4(I10,11X) /
	1	' Microprocessor Sweeps/IFG:     ', 4(I10,11X) /
	1	' MTM Speed at Transmit:         ', 4(I10,11X) /
	1	' MTM Length at Transmit:        ', 4(I10,11X) /
	1	' Science Gain:                  ', 4(I10,11X) /
	1	' Bolometer Commanded Bias:      ', 4(I10,11X) /
	1	' Bolometer Voltage Readout:     ', 4(I10,11X))

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C GRTs
C
	WRITE (WRITE_LUN, 6006, IOSTAT=STATUS)

6006	FORMAT (/3x,' <<<<<  GRT Temperature in Milli-Kelvins  >>>>>'/
	1	3x, ' ***** GRTs:  Low Current, A-side *****', t65,
	1	'***** GRTs: High Current, A-side *****')

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6007, IOSTAT=STATUS)
	1	(IDX_REC.GRTS(1).GROUP1(J),
	1	 IDX_REC.GRTS(2).GROUP1(J),
	1	 J=1,16)

6007	FORMAT ( 3X, ' External Calibrator:          ', I7, T65,
	1	' External Calibrator:          ', I7
	1	/ 3X, ' Sky Horn:                     ', I7, T65,
	1	' Sky Horn:                     ', I7
	1	/ 3X, ' Reference Horn:               ', I7, T65,
	1	' Reference Horn:               ', I7
	1	/ 3X, ' Internal Reference Source:    ', I7, T65,
	1	' Internal Reference Source:    ', I7
	1	/ 3X, ' Right Dihedral:               ', I7, T65,
	1	' Right Dihedral:               ', I7
	1	/ 3X, ' Bolometer Assembly, RH:       ', I7, T65,
	1	' Bolometer Assembly, RH:       ', I7
	1	/ 3X, ' Bolometer Assembly, RL:       ', I7, T65,
	1	' Bolometer Assembly, RL:       ', I7
	1	/ 3X, ' Bolometer Assembly, LH:       ', I7, T65,
	1	' Bolometer Assembly, LH:       ', I7
	1	/ 3X, ' Bolometer Assembly, LL:       ', I7, T65,
	1	' Bolometer Assembly, LL:       ', I7
	1	/ 3X, ' Right Mirror:                 ', I7, T65,
	1	' Right Mirror:                 ', I7
	1	/ 3X, ' Calibrator Resistor 1:        ', I7, T65,
	1	' Calibrator Resistor 1:        ', I7
	1	/ 3X, ' Calibrator Resistor 2:        ', I7, T65,
	1	' Calibrator Resistor 2:        ', I7
	1	/ 3X, ' Calibrator Resistor 3:        ', I7, T65,
	1	' Calibrator Resistor 3:        ', I7
	1	/ 3X, ' Calibrator Resistor 4:        ', I7, T65,
	1	' Calibrator Resistor 4:        ', I7
	1	/ 3X, ' XCal Segment 5:               ', I7, T65,
	1	' XCal Segment 5:               ', I7
	1	/ 3X, ' Right Collimator:             ', I7, T65,
	1	' Right Collimator:             ', I7 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6009, IOSTAT=STATUS)

6009	FORMAT (/3x, '***** GRTs:  Low Current, B-side *****', t65,
	1	'***** GRTs: High Current, B-side *****')

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6010, IOSTAT=STATUS)
	1	(IDX_REC.GRTS(3).GROUP1(J),
	1	 IDX_REC.GRTS(4).GROUP1(J),
	1	 J=1,16)

6010	FORMAT ( 3X, ' External Calibrator:          ', I7, T65,
	1	' External Calibrator:          ', I7
	1	/ 3X, ' Sky Horn:                     ', I7, T65,
	1	' Sky Horn:                     ', I7
	1	/ 3X, ' Reference Horn:               ', I7, T65,
	1	' Reference Horn:               ', I7
	1	/ 3X, ' Internal Reference Source:    ', I7, T65,
	1	' Internal Reference Source:    ', I7
	1	/ 3X, ' Left Dihedral:                ', I7, T65,
	1	' Left Dihedral:                ', I7
	1	/ 3X, ' Bolometer Assembly, RH:       ', I7, T65,
	1	' Bolometer Assembly, RH:       ', I7
	1	/ 3X, ' Bolometer Assembly, RL:       ', I7, T65,
	1	' Bolometer Assembly, RL:       ', I7
	1	/ 3X, ' Bolometer Assembly, LH:       ', I7, T65,
	1	' Bolometer Assembly, LH:       ', I7
	1	/ 3X, ' Bolometer Assembly, LL:       ', I7, T65,
	1	' Bolometer Assembly, LL:       ', I7
	1	/ 3X, ' Left Mirror:                  ', I7, T65,
	1	' Left Mirror:                  ', I7
	1	/ 3X, ' Calibrator Resistor 1:        ', I7, T65,
	1	' Calibrator Resistor 1:        ', I7
	1	/ 3X, ' Calibrator Resistor 2:        ', I7, T65,
	1	' Calibrator Resistor 2:        ', I7
	1	/ 3X, ' Calibrator Resistor 3:        ', I7, T65,
	1	' Calibrator Resistor 3:        ', I7
	1	/ 3X, ' Calibrator Resistor 4:        ', I7, T65,
	1	' Calibrator Resistor 4:        ', I7
	1	/ 3X, ' XCal Segment 6:               ', I7, T65,
	1	' XCal Segment 6:               ', I7
	1	/ 3X, ' Left Collimator:              ', I7, T65,
	1	' Left Collimator:              ', I7 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Temperature Controllers
C
	WRITE (WRITE_LUN, 6021, IOSTAT=STATUS)
	1			(IDX_REC.TEMP_CTRL.GROUP2(J),
	1			 IDX_REC.TEMP_CTRL.GROUP2(J+4), J=1,4)

6021	FORMAT (/3x, ' ***** Temperature Controllers *****'/
	1	3x, ' << Side A >>', t65, ' << Side B >>' /
	1	' Internal Calibrator: ', I7, t65,
	1	' Internal Calibrator: ', i7 /
	1	' Reference Horn:      ', i7, t65,
	1	' Reference Horn:      ', i7 /
	1	' Sky Horn:            ', i7, t65,
	1	' Sky Horn:            ', i7 /
	1	' External Calibrator: ', i7, t65,
	1	' External Calibrator: ', i7 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Status Monitor Individual Statuses
C
	WRITE (WRITE_LUN, 6022, IOSTAT=STATUS)
	1			(IDX_REC.STAT_MON(1).GROUP3(J),
	1			 IDX_REC.STAT_MON(2).GROUP3(J),
	1			 J=1,29),
	1			 IDX_REC.IDX_SYNC.SIDE(1),
	1			 IDX_REC.IDX_SYNC.SIDE(2)

6022	FORMAT (/3x, ' ***** Status Monitor Individual Statuses *****'/
	1	3x, ' << Side A >>', t65, ' << Side B >>' /
	1	' MTM Operating Mode:         ', I2, t65,
	1	' MTM Operating Mode:         ', I2 /
	1	' ICAL TC Integrator Enable:  ', I2, t65,
	1	' ICAL TC Integrator Enable:  ', I2 /
	1	' REFH TC Integrator Enable:  ', I2, t65,
	1	' REFH TC Integrator Enable:  ', I2 /
	1	' SKYH TC Integrator Enable:  ', I2, t65,
	1	' SKYH TC Integrator Enable:  ', I2 /
	1	' XCAL TC Integrator Enable:  ', I2, t65,
	1	' XCAL TC Integrator Enable:  ', I2 /
	1	' Dither Channel RH:          ', I2, T65,
	1	' Dither Channel LH:          ', I2 /
	1	' Dither Channel RL:          ', I2, T65,
	1	' Dither Channel LL:          ', I2 /
	1	' ROM Power RH:               ', I2, T65,
	1	' ROM Power LH:               ', I2 /
	1	' ROM Power RL:               ', I2, T65,
	1	' ROM Power LL:               ', I2 /
	1	' MTM Latch:                  ', I2, t65,
	1	' MTM Latch:                  ', I2 /
	1	' XCAL Latch:                 ', i2, t65,
	1	' XCAL Latch:                 ', i2 /
	1	' ICAL TC Current Range:      ', I2, t65,
	1	' ICAL TC Current Range:      ', I2 /
	1	' REFH TC Current Range:      ', I2, t65,
	1	' REFH TC Current Range:      ', I2 /
	1	' SKYH TC Current Range:      ', I2, t65,
	1	' SKYH TC Current Range:      ', I2 /
	1	' XCAL TC Current Range:      ', I2, t65,
	1	' XCAL TC Current Range:      ', I2 /
	1	' ICAL TC Integrator Gain:    ', I2, t65,
	1	' ICAL TC Integrator Gain:    ', I2 /
	1	' REFH TC Integrator Gain:    ', I2, t65,
	1	' REFH TC Integrator Gain:    ', I2 /
	1	' SKYH TC Integrator Gain:    ', I2, t65,
	1	' SKYH TC Integrator Gain:    ', I2 /
	1	' XCAL TC Integrator Gain:    ', I2, t65,
	1	' XCAL TC Integrator Gain:    ', I2 /
	1	' ICAL TC Proportional Gain:  ', I2, t65,
	1	' ICAL TC Proportional Gain:  ', I2 /
	1	' REFH TC Proportional Gain:  ', I2, t65,
	1	' REFH TC Proportional Gain:  ', I2 /
	1	' SKYH TC Proportional Gain:  ', I2, t65,
	1	' SKYH TC Proportional Gain:  ', I2 /
	1	' XCAL TC Proportional Gain:  ', I2, t65,
	1	' XCAL TC Proportional Gain:  ', I2 /
	1	' MTM Stow:                   ', I2, T65,
	1	' MTM Stow:                   ', I2 /
	1	' MTM Speed (Commanded):      ', I2, t65,
	1	' MTM Speed (Commanded):      ', I2 /
	1	' MTM Length (Commanded):     ', I2, t65,
	1	' MTM Length (Commanded):     ', I2 /
	1	' XCAL Stow:                  ', I2, t65,
	1	' XCAL Stow:                  ', I2 /
	1	' XCAL Latch (Commanded):     ', I2, t65,
	1	' XCAL Latch (Commanded):     ', I2 /
	1	' MTM Latch (Commanded):      ', I2, t65,
	1	' MTM Latch (Commanded):      ', I2 /
	1	' Stat.Mon. Sync (1=in sync): ', i2, t65,
	1	' Stat.Mon. Sync (1=in sync): ', i2)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C IPDU Power Relay Status
C
	WRITE (WRITE_LUN, 6023, IOSTAT=STATUS)
	1		(IDX_REC.IPDU_RELAY.GROUP4(J), J=1,8)

6023	FORMAT (/3x, ' ***** IPDU Power Relay Status ***** ' /
	1	' IPDU Power Relay 4-Pole A1 (HEX): ', 2X, Z2.2, t65,
	1	' IPDU Power Relay 4-Pole B1 (HEX): ', 2X, Z2.2 /
	1	' IPDU Power Relay 4-Pole A2 (HEX): ', 2X, Z2.2, t65,
	1	' IPDU Power Relay 4-Pole B2 (HEX): ', 2X, Z2.2 /
	1	' IPDU Power Relay 2-Pole A3 (HEX): ', 2X, Z2.2, t65,
	1	' IPDU Power Relay 2-Pole B3 (HEX): ', 2X, Z2.2 /
	1	' IPDU Other Status A4 (HEX):       ', 2X, Z2.2, t65,
	1	' IPDU Other Status B4 (HEX):       ', 2X, Z2.2 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C IPDU Individual Statuses
C
	WRITE (WRITE_LUN, 6024, IOSTAT=STATUS)
	1			(IDX_REC.IPDU_STAT(1).GROUP5(J),
	1			 IDX_REC.IPDU_STAT(2).GROUP5(J),
	1			 J=1,15)

6024	FORMAT (/3X, ' ***** IPDU Individual Statuses *****'/
	1	' << Side A >>', T65, ' << Side B >>'/
	1	' MTM Motor Power Enable:     ', I4, t65,
	1	' MTM Motor Power Enable:     ', I4 /
	1	' Main Elec. Converters RH:   ', i4, t65,
	1	' Main Elec. Converters LH:   ', i4 /
	1	' Main Elec. Converters RL:   ', i4, t65,
	1	' Main Elec. Converters LL:   ', i4 /
	1	' Hot Spot Heater:            ', i4, t65,
	1	' Hot Spot Heater:            ', i4 /
	1	' MTM Drive Motor Elect.:     ', i4, t65,
	1	' MTM Drive Motor Elect.:     ', i4 /
	1	' XCAL TC Power:              ', i4, t65,
	1	' XCAL TC Power:              ', i4 /
	1	' SKYH TC Power:              ', i4, t65,
	1	' SKYH TC Power:              ', i4 /
	1	' REFH TC Power:              ', i4, t65,
	1	' REFH TC Power:              ', i4 /
	1	' ICAL TC Power:              ', i4, t65,
	1	' ICAL TC Power:              ', i4 /
	1	' MTM Latch Motor Elect.:     ', i4, t65,
	1	' MTM Latch Motor Elect.:     ', i4 /
	1	' XCAL Latch Motor Elect.:    ', i4, t65,
	1	' XCAL Latch Motor Elect.:    ', i4 /
	1	' XCAL Drive Motor Elect.:    ', i4, t65,
	1	' XCAL Drive Motor Elect.:    ', i4 /
	1	' XCAL Power Enable:          ', i4, t65,
	1	' XCAL Power Enable:          ', i4 /
	1	' LMAC for DC Convert.:       ', i4, t65,
	1	' LMAC for DC Convert.:       ', i4 /
	1	' LMAC for AC Convert.:       ', i4, t65,
	1	' LMAC for AC Convert.:       ', i4 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Dwell, Microprocessor Bus Readouts, SIS statuses
C
	WRITE (WRITE_LUN, 6025, IOSTAT=STATUS)
	1			IDX_REC.MISC_STAT.DWELL(1),
	1			IDX_REC.MISC_STAT.DWELL(3),
	1			IDX_REC.MISC_STAT.DWELL(2),
	1			IDX_REC.MISC_STAT.DWELL(4),
	1			IDX_REC.MISC_STAT.LVDT_STAT(1),
	1			IDX_REC.MISC_STAT.LVDT_STAT(2),
	1			IDX_REC.MISC_STAT.STAT_RD_BUS(1),
	1			IDX_REC.MISC_STAT.STAT_RD_BUS(3),
	1			IDX_REC.MISC_STAT.STAT_RD_BUS(2),
	1			IDX_REC.MISC_STAT.STAT_RD_BUS(4)

6025	FORMAT (/3x, ' ***** Dwell, Microprocessor Bus Readouts *****'/
	1	' << Side A >>', t65, ' << Side B >>' /
	1	' Dwell Mode (1=Dwell):       ', I4, t65,
	1	' Dwell Mode (1=Dwell):       ', I4 /
	1	' Dwell Address:              ', I4, t65,
	1	' Dwell Address:              ', I4 /
	1	' LVDT Status:                ', I4, t65,
	1	' LVDT Status:                ', I4,/
	1	' Microproc. Bus Readout RH:  ', I4, T65,
	1	' Microproc. Bus Readout LH:  ', I4 /
	1	' Microproc. Bus Readout RL:  ', I4, T65,
	1	' Microproc. Bus Readout LL:  ', I4 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6026, IOSTAT=STATUS)
	1		((IDX_REC.MISC_STAT.POWER_A_STATUS(J),
	2		  IDX_REC.MISC_STAT.POWER_B_STATUS(J)),J=1,4)

6026	FORMAT (/3x, ' ***** Special SIS Statuses *****'/
	1	' IPDU Int.Converter, Side A: ', I4, t65,
	1	' IPDU Int.Converter, Side B: ', I4 /
	1	' IPDU A-D Converter, Side A: ', I4, t65,
	1	' IPDU A-D Converter, Side B: ', I4 /
	1	' IPDU Bias Pre-Reg.Conv. A:  ', i4, t65,
	1	' IPDU Bias Pre-Reg.Conv. B:  ', i4 /
	1	' IPDU MTM and XCAL, Side A:  ', i4, t65,
	1	' IPDU MTM and XCAL, Side B:  ', i4 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6031, IOSTAT=STATUS)
	1		(IDX_REC.MISC_STAT.NTCH_FLTR_A(J),
	1		IDX_REC.MISC_STAT.NTCH_FLTR_B(J), J=1,5)

6031	FORMAT (/3x, ' ***** Notch Filter Statuses *****'/
	1	' Notch Filter 1, Side A: ', Z2.2, t65,
	1	' Notch Filter 1, Side B: ', Z2.2, /
	1	' Notch Filter 2, Side A: ', Z2.2, t65,
	1	' Notch Filter 2, Side B: ', Z2.2, /
	1	' Notch Filter 3, Side A: ', Z2.2, t65,
	1	' Notch Filter 3, Side B: ', Z2.2, /
	1	' Notch Filter 4, Side A: ', Z2.2, t65,
	1	' Notch Filter 4, Side B: ', Z2.2, /
	1	' Notch Filter 5, Side A: ', Z2.2, t65,
	1	' Notch Filter 5, Side B: ', Z2.2 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


C
C Other temperatures/currents
C
	DO	J = 1, 14
	  USBI(J) = ZEXT(IDX_REC.IDX_ANALOG.TEMP(J))
	ENDDO

	WRITE (WRITE_LUN, 6027, IOSTAT=STATUS) (USBI(J), J=1,14)

6027	FORMAT (/3x, ' ***** Other Temperatures/Currents *****' /
	1	' IPDU Temperature, Side A:  ', i4, t65,
	1	' IPDU Temperature, Side B:  ', i4 /
	1	' Drive Box Temp., Side A:   ', i4, t65,
	1	' Channel Pre-Amp:           ', i4 /
	1	' Channel Temperature, RH:   ', i4, t65,
	1	' Channel Temperature, RL:   ', i4 /
	1	' Channel Temperature, LH:   ', i4, t65,
	1	' Channel Temperature, LL:   ', i4 /
	1	' Stat.Mon. Temp. Side A:    ', i4, t65,
	1	' Stat.Mon. Temp. Side B:    ', i4 /
	1	' Hot Spot Heater Curr., A:  ', i4, t65,
	1	' Hot Spot Heater Curr., B:  ', i4 /
	1	' Optical Pre-Amp.:          ', i4, t65,
	1	' Drive Box Temp. Side B:    ', i4 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C IPDU Voltages
C
	DO	J = 1, 20
	  USBI(J) = ZEXT(IDX_REC.IDX_ANALOG.V_AND_I(J))
	ENDDO

	WRITE (WRITE_LUN, 6014, IOSTAT=STATUS)	(USBI(J),J=1,20)

6014	FORMAT ( /3x, ' ***** IPDU Voltages *****' /
	1	' Digital Converter, -15V, Side A: ', I7, t65,
	1	' Digital Converter, -15V, Side B: ', I7
	1	/' Digital Converter, +15V, Side A: ', I7, t65,
	1	' Digital Converter, +15V, Side B: ', I7
	1	/' Digital Converter, +5V, Side A:  ', I7, t65,
	1	' Digital Converter, +5V, Side B:  ', I7
	1	/' Analog Converter, +15V, Side A:  ', I7, t65,
	1	' Analog Converter, +15V, Side B:  ', I7
	1	/' Analog Converter, -15V, Side A:  ', I7, t65,
	1	' Analog Converter, -15V, Side B:  ', I7
	1	/' Bias Pre-Reg. +25V, Side A:      ', I7, t65,
	1	' Bias Pre-Reg. +25V, Side B:      ', I7
	1	/' Int. P.S., +28V, Side A:         ', I7, t65,
	1	' Int. P.S., +28V, Side B:         ', I7
	1	/' Int. P.S., +15V, Side A:         ', I7, t65,
	1	' Int. P.S., +15V, Side B:         ', I7
	1	/' Int. P.S., -15V, Side A:         ', I7, t65,
	1	' Int. P.S., -15V, Side B:         ', I7
	1	/' Int. P.S., +5V, Side A:          ', I7, t65
	1	' Int. P.S., +5V, Side B:          ', I7 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C IPDU Currents
C
	DO J=1,12
	  USBI(J) = ZEXT(IDX_REC.IDX_ANALOG.V_AND_I(J+20))
	ENDDO

	WRITE (WRITE_LUN, 6015, IOSTAT=STATUS)	(USBI(J),J=1,12)

6015	FORMAT ( /3x, ' ***** IPDU Currents *****' 
	1	/' Bias Pre-Reg. Current, Side A:   ', I7, t65,
	1	' Bias Pre-Reg. Current, Side B:   ', I7
	1	/' Analog Conv. Current, Side A:    ', I7, t65,
	1	' Analog Conv. Current, Side B:    ', I7
	1	/' Digital Conv. Current, Side A:   ', I7, t65,
	1	' Digital Conv. Current, Side B:   ', I7
	1	/' Constant Current RH:             ', I7, t65,
	1	' Constant Current  RL:             ', I7
	1	/' Constant Current LH:             ', I7, t65,
	1	' Constant Current LL:             ', I7
	1	/' Int. Conv. Current, Side A:      ', I7, t65,
	1	' Int. Conv. Current, Side B:      ', I7 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C MTM Cal Motor
C
	WRITE (WRITE_LUN, 6302, IOSTAT=STATUS)
	1			IDX_REC.IDX_ANALOG.MTM_CAL_MOTOR_SIDE_A,
	1			IDX_REC.IDX_ANALOG.MTM_CAL_MOTOR_SIDE_B

6302	FORMAT (/3X, ' MTM Cal Motor, Side A:      ', i7, t65,
	1	     ' MTM Cal Motor, Side B:      ', i7 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Tolerance Values...
C
	DO	J = 1, 22
	  USBI(J) = ZEXT(IDX_REC.IDX_TAIL.TOLS(J))
	ENDDO

	WRITE (WRITE_LUN, 6303, IOSTAT=STATUS) (USBI(J), J=1,22)

6303	FORMAT (/3x, ' ***** Tolerance Values *****' /
	1	' CONTLD_TEMP_LO_TOL:   ', i4, t65,
	1	' CONTLD_TEMP_HI_TOL:   ', i4 /
	1	' DIHED_TEMP_LO_TOL:    ', i4, t65,
	1	' DIHED_TEMP_HI_TOL:    ', i4 /
	1	' BOLOM_TEMP_LO_TOL:    ', i4, t65,
	1	' BOLOM_TEMP_HI_TOL:    ', i4 /
	1	' MIRR_TEMP_LO_TOL:     ', i4, t65,
	1	' MIRR_TEMP_HI_TOL:     ', i4 /
	1	' TEMP_CONTROLLER_TOL:  ', i4, t65,
	1	' OTHER_TEMP_TOL  :     ', i4 /
	1	' BOLOM_BIAS_RDOUT_TOL: ', i4, t65,
	1	' DIG_CONV_VOL_TOL:     ', i4 /
	1	' ANLG_CONV_VOL_TOL:    ', i4, t65,
	1       ' BIAS_PRE_REG_VOL_TOL: ', i4, /
	1	' INT_PWR_SUPP_28V_TOL: ', i4, t65,
	1	' INT_PWR_SUPP_15V_TOL: ', i4, /
	1	' INT_PWR_SUPP_5V_TOL:  ', i4, t65,
	1	' BIAS_PRE_REG_CUR_TOL: ', i4, /
	1	' ANLG_CONV_CURR_TOL:   ', i4, t65,
	1	' DIG_CONV_CURR_TOL:    ', i4, /
	1	' CONST_CURR_TOL:       ', i4, t65,
	1	' INTERNAL_CONV_TOL:    ', i4 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

	RETURN
	END
