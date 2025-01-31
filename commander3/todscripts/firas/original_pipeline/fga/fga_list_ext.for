	SUBROUTINE FGA_LIST_EXT (WRITE_LUN, EXT_REC, NREC, LINE,
	1	CAT_ENTRY)

C------------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FGA_LIST_EXT
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given FIRAS
C	  Engineering Extrema record.
C
C	AUTHOR:
C	  R. Kummerer
C	  STX
C	  April 26, 1989
C
C	  Adapted from work performed by E. Fung.
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_EXT (WRITE_LUN, EXT_REC, NREC, LINE, CAT_ENTRY)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list
C					file
C
C	  EXT_REC		RECORD	FIRAS Engineering Extrema record.
C
C	  NREC			I*2	The number of records processed in this
C					run.
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
C------------------------------------------------------------------------

	IMPLICIT	NONE

	CHARACTER	*80	LINE
	INTEGER		*2	NREC
	INTEGER		*4	WRITE_LUN
	INTEGER		*4	CAT_ENTRY
	INTEGER		*4	STATUS
	INTEGER		*2	I
	INTEGER		*2	J
	CHARACTER	*14	GMT(64)

	DICTIONARY	'FXT_ENG_XTRM'
	RECORD /FXT_ENG_XTRM/ EXT_REC

	EXTERNAL		FGA_WRITE_ERR

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6


	WRITE (WRITE_LUN, 6000, IOSTAT=STATUS)	LINE, NREC, CAT_ENTRY

6000	FORMAT ( '1' / 7X, 'FIRAS Engineering Extrema Formatted Listing' /
	1	'***  ', A, '  ***'
	1	/ '   Record Count for this Run: ***', I6, '  ***'
	1	10X, 'CT Catalog Entry: ', I7)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C	COBETRIEVE standard header
C
	WRITE (WRITE_LUN, 6001, IOSTAT=STATUS)
	1			EXT_REC.CT_HEAD.GMT(1:2),
	1			EXT_REC.CT_HEAD.GMT(3:5),
	1			EXT_REC.CT_HEAD.GMT(6:7),
	1			EXT_REC.CT_HEAD.GMT(8:9),
	1			EXT_REC.CT_HEAD.GMT(10:11),
	1			EXT_REC.CT_HEAD.GMT(12:14),
	1			EXT_REC.CT_HEAD.TIME(2),
	1			EXT_REC.CT_HEAD.TIME(1),
	1			EXT_REC.CT_HEAD.SPACE_TIME,
	1			EXT_REC.CT_HEAD.ORBIT

6001	FORMAT (/' ASCII GMT:                             ', 4X,
	1	A2, '-', A3, '-', 3(A2, '-'), A3,
	1	5X, ' GMT in binary (high order 4 bytes first): ', 2X,
	1	2(Z8.8,2X)
	1	/' S/C time in PB4 (HEX):                 ', 6X,
	1	6(Z2.2,1X), 5X,
	1	' Orbit Number:                          ',12X,
	1	I11)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Minimum engineering analog quantities.
C
	WRITE (WRITE_LUN, 6200, IOSTAT=STATUS)

6200	FORMAT (/' ******** Minimum Engineering Analog Quantities ********')

C
C GRT's
C
	WRITE (WRITE_LUN, 6005, IOSTAT=STATUS)

6005	FORMAT (/' ***** GRT''S:  Low Current, A-side *****' ,
	1	T65, ' ***** GRT''S:  High Current, A-side *****' )
     
	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	DO I = 1, 64
	  CALL CT_BINARY_TO_GMT(EXT_REC.MINMAX(1).GRT_TIME(1,I),GMT(I))
	END DO

	WRITE (WRITE_LUN, 6006, IOSTAT=STATUS)
	1			(EXT_REC.MINMAX(1).GRT(I), GMT(I),
	1			 EXT_REC.MINMAX(1).GRT(I+16), GMT(I+16),
	1			 I=1,16)

6006	FORMAT (' External Calibrator:          ', G14.5, X, A, T65,
	1	' External Calibrator:          ', G14.5, X, A
	1	/ ' Sky Horn:                     ', G14.5, X, A, T65,
	1	' Sky Horn:                     ', G14.5, X, A
	1	/ ' Reference Horn:               ', G14.5, X, A, T65,
	1	' Reference Horn:               ', G14.5, X, A
	1	/ ' Internal Reference Source:    ', G14.5, X, A, T65,
	1	' Internal Reference Source:    ', G14.5, X, A
	1	/ ' Right Dihedral:               ', G14.5, X, A, T65,
	1	' Right Dihedral:               ', G14.5, X, A
	1	/ ' Bolometer Assembly, RH:       ', G14.5, X, A, T65,
	1	' Bolometer Assembly, RH:       ', G14.5, X, A
	1	/ ' Bolometer Assembly, RL:       ', G14.5, X, A, T65,
	1	' Bolometer Assembly, RL:       ', G14.5, X, A
	1	/ ' Bolometer Assembly, LH:       ', G14.5, X, A, T65,
	1	' Bolometer Assembly, LH:       ', G14.5, X, A
	1	/ ' Bolometer Assembly, LL:       ', G14.5, X, A, T65,
	1	' Bolometer Assembly, LL:       ', G14.5, X, A
	1	/ ' Right Mirror:                 ', G14.5, X, A, T65,
	1	' Right Mirror:                 ', G14.5, X, A
	1	/ ' Calibrator Resistor 1:        ', G14.5, X, A, T65,
	1	' Calibrator Resistor 1:        ', G14.5, X, A
	1	/ ' Calibrator Resistor 2:        ', G14.5, X, A, T65,
	1	' Calibrator Resistor 2:        ', G14.5, X, A
	1	/ ' Calibrator Resistor 3:        ', G14.5, X, A, T65,
	1	' Calibrator Resistor 3:        ', G14.5, X, A
	1	/ ' Calibrator Resistor 4:        ', G14.5, X, A, T65,
	1	' Calibrator Resistor 4:        ', G14.5, X, A
	1	/ ' XCal Segment 5:               ', G14.5, X, A, T65,
	1	' XCal Segment 5:               ', G14.5, X, A
	1	/ ' Right Collimator:             ', G14.5, X, A, T65,
	1	' Right Collimator:             ', G14.5, X, A )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6008, IOSTAT=STATUS)

6008	FORMAT (/' ***** GRT''S:  Low Current, B-side *****', T65,
	1	' ***** GRT''S:  High Current, B-side *****'  )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6009, IOSTAT=STATUS)
	1			(EXT_REC.MINMAX(1).GRT(I+32), GMT(I+32),
	1			 EXT_REC.MINMAX(1).GRT(I+48), GMT(I+48),
	1			 I=1,16)

6009	FORMAT (' External Calibrator:          ', G14.5, X, A, T65,
	1	' External Calibrator:          ', G14.5, X, A
	1	/ ' Sky Horn:                     ', G14.5, X, A, T65,
	1	' Sky Horn:                     ', G14.5, X, A
	1	/ ' Reference Horn:               ', G14.5, X, A, T65,
	1	' Reference Horn:               ', G14.5, X, A
	1	/ ' Internal Reference Source:    ', G14.5, X, A, T65,
	1	' Internal Reference Source:    ', G14.5, X, A
	1	/ ' Left Dihedral:                ', G14.5, X, A, T65,
	1	' Left Dihedral:                ', G14.5, X, A
	1	/ ' Bolometer Assembly, RH:       ', G14.5, X, A, T65,
	1	' Bolometer Assembly, RH:       ', G14.5, X, A
	1	/ ' Bolometer Assembly, RL:       ', G14.5, X, A, T65,
	1	' Bolometer Assembly, RL:       ', G14.5, X, A
	1	/ ' Bolometer Assembly, LH:       ', G14.5, X, A, T65,
	1	' Bolometer Assembly, LH:       ', G14.5, X, A
	1	/ ' Bolometer Assembly, LL:       ', G14.5, X, A, T65,
	1	' Bolometer Assembly, LL:       ', G14.5, X, A
	1	/ ' Left Mirror:                  ', G14.5, X, A, T65,
	1	' Left Mirror:                  ', G14.5, X, A
	1	/ ' Calibrator Resistor 1:        ', G14.5, X, A, T65,
	1	' Calibrator Resistor 1:        ', G14.5, X, A
	1	/ ' Calibrator Resistor 2:        ', G14.5, X, A, T65,
	1	' Calibrator Resistor 2:        ', G14.5, X, A
	1	/ ' Calibrator Resistor 3:        ', G14.5, X, A, T65,
	1	' Calibrator Resistor 3:        ', G14.5, X, A
	1	/ ' Calibrator Resistor 4:        ', G14.5, X, A, T65,
	1	' Calibrator Resistor 4:        ', G14.5, X, A
	1	/ ' XCal Segment 6:               ', G14.5, X, A, T65,
	1	' XCal Segment 6:               ', G14.5, X, A
	1	/ ' Left Collimator (Not used):   ', G14.5, X, A, T65,
	1	' Left Collimator (Not used):   ', G14.5, X, A )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C IPDU temperatures, etc.
C
	WRITE (WRITE_LUN, 6012, IOSTAT=STATUS)

6012	FORMAT (/' ***** OTHER TEMPERATURES/CURRENTS *****' )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	DO I = 1, 16
	  CALL CT_BINARY_TO_GMT(EXT_REC.MINMAX(1).T_AND_I_TIME(1,I),GMT(I))
	END DO

	WRITE (WRITE_LUN, 6013, IOSTAT=STATUS)
	1			(EXT_REC.MINMAX(1).T_AND_I(I), GMT(I), I=1,16)

6013	FORMAT (' IPDU temperature - Side A:       ', G14.5, X, A, T65,
	1	' IPDU temperature - Side B:       ', G14.5, X, A
	1	/' Channel temperature - RH:        ', G14.5, X, A, T65,
	1	' Channel temperature - LH:        ', G14.5, X, A
	1	/' Channel temperature - RL:        ', G14.5, X, A, T65,
	1	' Channel temperature - LL:        ', G14.5, X, A
	1	/' Drive Box temperature - Side A:  ', G14.5, X, A, T65,
	1	' Drive Box temperature - Side B:  ', G14.5, X, A
	1	/' Status Monitor temp. - Side A:   ', G14.5, X, A, T65,
	1	' Status Monitor temp. - Side B:   ', G14.5, X, A
	1	/' Channel Pre-Amp temperature:     ', G14.5, X, A, T65,
	1	' Optical Pre-Amp temperature:     ', G14.5, X, A
	1	/' Hot Spot Heater Current-Side A:  ', G14.5, X, A, T65,
	1	' Hot Spot Heater Current-Side B:  ', G14.5, X, A
	1	/' MTM Cal Motor - Side A:          ', G14.5, X, A, T65,
	1	' MTM Cal Motor - Side B:          ', G14.5, X, A )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Voltages and Currents
C
	DO I = 1, 36
	  CALL CT_BINARY_TO_GMT(EXT_REC.MINMAX(1).V_AND_I_TIME(1,I),GMT(I))
	END DO

	WRITE (WRITE_LUN, 6014, IOSTAT=STATUS)	(EXT_REC.MINMAX(1).V_AND_I(I),
	1					 GMT(I), I=1,36)

6014	FORMAT ( /' ***** Various Voltages *****' 
	1	/' Bolometer Bias Voltage, RH:      ', G14.5, X, A, T65,
	1	' Bolometer Bias Voltage, RL:      ', G14.5, X, A
	1	/' Bolometer Bias Voltage, LH:      ', G14.5, X, A, T65,
	1	' Bolometer Bias Voltage, LL:      ', G14.5, X, A
	1	/' Digital Converter, -15V, Side A: ', G14.5, X, A, T65,
	1	' Digital Converter, -15V, Side B: ', G14.5, X, A
	1	/' Digital Converter, +15V, Side A: ', G14.5, X, A, T65,
	1	' Digital Converter, +15V, Side B: ', G14.5, X, A
	1	/' Digital Converter, +5V, Side A:  ', G14.5, X, A, T65,
	1	' Digital Converter, +5V, Side B:  ', G14.5, X, A
	1	/' Analog Converter, +15V, Side A:  ', G14.5, X, A, T65,
	1	' Analog Converter, +15V, Side B:  ', G14.5, X, A
	1	/' Analog Converter, -15V, Side A:  ', G14.5, X, A, T65,
	1	' Analog Converter, -15V, Side B:  ', G14.5, X, A
	1	/' Bias Pre-Reg. +25V, Side A:      ', G14.5, X, A, T65,
	1	' Bias Pre-Reg. +25V, Side B:      ', G14.5, X, A
	1	/' Int. P.S., +28V, Side A:         ', G14.5, X, A, T65,
	1	' Int. P.S., +28V, Side B:         ', G14.5, X, A
	1	/' Int. P.S., +15V, Side A:         ', G14.5, X, A, T65,
	1	' Int. P.S., +15V, Side B:         ', G14.5, X, A
	1	/' Int. P.S., -15V, Side A:         ', G14.5, X, A, T65,
	1	' Int. P.S., -15V, Side B:         ', G14.5, X, A
	1	/' Int. P.S., +5V, Side A:          ', G14.5, X, A, T65,
	1	' Int. P.S., +5V, Side B:          ', G14.5, X, A
	1	//' ***** Various Currents *****' 
	1	/' Bias Pre-Reg. Current, Side A:   ', G14.5, X, A, T65,
	1	' Bias Pre-Reg. Current, Side B:   ', G14.5, X, A
	1	/' Analog Conv. Current, Side A:    ', G14.5, X, A, T65,
	1	' Analog Conv. Current, Side B:    ', G14.5, X, A
	1	/' Digital Conv. Current, Side A:   ', G14.5, X, A, T65,
	1	' Digital Conv. Current, Side B:   ', G14.5, X, A
	1	/' Constant Current RH:             ', G14.5, X, A, T65,
	1	' Constant Current LH:             ', G14.5, X, A
	1	/' Constant Current RL:             ', G14.5, X, A, T65,
	1	' Constant Current LL:             ', G14.5, X, A
	1	/' Int. Conv. Current, Side A:      ', G14.5, X, A, T65,
	1	' Int. Conv. Current, Side B:      ', G14.5, X, A )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C LMACs
C
	DO I = 1, 2
	  CALL CT_BINARY_TO_GMT(EXT_REC.MINMAX(1).LMAC_TIME(1,I),GMT(I))
	END DO

	WRITE (WRITE_LUN, 6400, IOSTAT=STATUS)	(EXT_REC.MINMAX(1).LMACS(I),
	1					 GMT(I), I=1,2)

6400	FORMAT (/' LMAC Analog Temperature:    ', G14.5, X, A, T65,
	1	/' LMAC Digital Temperature:   ', G14.5, X, A)


C
C Maximum engineering analog quantities.
C
	WRITE (WRITE_LUN, 6300, IOSTAT=STATUS)

6300	FORMAT (/' ******** Maximum Engineering Analog Quantities ********')

C
C GRT's
C
	WRITE (WRITE_LUN, 6005, IOSTAT=STATUS)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	DO I = 1, 64
	  CALL CT_BINARY_TO_GMT(EXT_REC.MINMAX(2).GRT_TIME(1,I),GMT(I))
	END DO

	WRITE (WRITE_LUN, 6006, IOSTAT=STATUS)
	1			(EXT_REC.MINMAX(2).GRT(I), GMT(I),
	1			 EXT_REC.MINMAX(2).GRT(I+16), GMT(I+16),
	1			 I=1,16)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6008, IOSTAT=STATUS)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6009, IOSTAT=STATUS)
	1			(EXT_REC.MINMAX(2).GRT(I+32), GMT(I+32),
	1			 EXT_REC.MINMAX(2).GRT(I+48), GMT(I+48),
	1			 I=1,16)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C IPDUs, etc.
C
	WRITE (WRITE_LUN, 6012, IOSTAT=STATUS)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	DO I = 1, 16
	  CALL CT_BINARY_TO_GMT(EXT_REC.MINMAX(2).T_AND_I_TIME(1,I),GMT(I))
	END DO

	WRITE (WRITE_LUN, 6013, IOSTAT=STATUS)
	1			(EXT_REC.MINMAX(2).T_AND_I(I), GMT(I), I=1,16)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Voltages and Currents
C
	DO I = 1, 36
	  CALL CT_BINARY_TO_GMT(EXT_REC.MINMAX(2).V_AND_I_TIME(1,I),GMT(I))
	END DO

	WRITE (WRITE_LUN, 6014, IOSTAT=STATUS)	(EXT_REC.MINMAX(2).V_AND_I(I),
	1					 GMT(I), I=1,36)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C LMACs
C
	DO I = 1, 2
	  CALL CT_BINARY_TO_GMT(EXT_REC.MINMAX(2).LMAC_TIME(1,I),GMT(I))
	END DO

	WRITE (WRITE_LUN, 6400, IOSTAT=STATUS)	(EXT_REC.MINMAX(2).LMACS(I),
	1					 GMT(I), I=1,2)


	RETURN
	END
