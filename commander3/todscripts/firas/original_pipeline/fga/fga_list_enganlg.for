	SUBROUTINE FGA_LIST_ENGANLG (WRITE_LUN, ANLG_REC)

C------------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FGA_LIST_ENGANLG
C
C	PROGRAM DESCRIPTION:
C	  This program will get a formatted listing of a given FIRAS
C	  Engineering analog fields.
C
C	AUTHOR:
C	  R. Kummerer
C	  STX
C	  April 7, 1988
C
C		Adapted from work done by Ed Fung.
C
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_ENGANLG (WRITE_LUN, ANLG_REC)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN		I*4	Logical unit number for a list file;
C					if not supplied assume 6.  Assume the
C					calling routine has opened this list file
C
C	  ANLG_REC		RECORD	FIRAS Engineering record.
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
C------------------------------------------------------------------------

	IMPLICIT	NONE

	INTEGER		*4	WRITE_LUN
	INTEGER		*4	STATUS
	CHARACTER	*7	REC_STAT
	CHARACTER	*14	GMT(4)
	INTEGER		*2	I
	INTEGER		*2	J

	DICTIONARY	'FUT_ENGANLG'
	RECORD /FUT_ENGANLG/ ANLG_REC

	DICTIONARY	'FDQ_ENG'
	RECORD /FDQ_ENG/ ENG_REC

	EXTERNAL		FGA_WRITE_ERR

	IF (WRITE_LUN .EQ. 0) WRITE_LUN = 6

C
C GRT's
C
	WRITE (WRITE_LUN, 6005, IOSTAT=STATUS)

6005	FORMAT (/' ***** GRT''S:  Low Current, A-side *****' ,
	1	T65, ' ***** GRT''S:  High Current, A-side *****' )
     
	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6006, IOSTAT=STATUS)
	1			(ANLG_REC.A_LO_GRT(I),
	1			 ANLG_REC.A_HI_GRT(I), I=1,16)

6006	FORMAT (' External Calibrator:          ', G14.5, T65,
	1	' External Calibrator:          ', G14.5
	1	/ ' Sky Horn:                     ', G14.5, T65,
	1	' Sky Horn:                     ', G14.5
	1	/ ' Reference Horn:               ', G14.5, T65,
	1	' Reference Horn:               ', G14.5
	1	/ ' Internal Reference Source:    ', G14.5, T65,
	1	' Internal Reference Source:    ', G14.5
	1	/ ' Right Dihedral:               ', G14.5, T65,
	1	' Right Dihedral:               ', G14.5
	1	/ ' Bolometer Assembly, RH:       ', G14.5, T65,
	1	' Bolometer Assembly, RH:       ', G14.5
	1	/ ' Bolometer Assembly, RL:       ', G14.5, T65,
	1	' Bolometer Assembly, RL:       ', G14.5
	1	/ ' Bolometer Assembly, LH:       ', G14.5, T65,
	1	' Bolometer Assembly, LH:       ', G14.5
	1	/ ' Bolometer Assembly, LL:       ', G14.5, T65,
	1	' Bolometer Assembly, LL:       ', G14.5
	1	/ ' Right Mirror:                 ', G14.5, T65,
	1	' Right Mirror:                 ', G14.5
	1	/ ' Calibrator Resistor 1:        ', G14.5, T65,
	1	' Calibrator Resistor 1:        ', G14.5
	1	/ ' Calibrator Resistor 2:        ', G14.5, T65,
	1	' Calibrator Resistor 2:        ', G14.5
	1	/ ' Calibrator Resistor 3:        ', G14.5, T65,
	1	' Calibrator Resistor 3:        ', G14.5
	1	/ ' Calibrator Resistor 4:        ', G14.5, T65,
	1	' Calibrator Resistor 4:        ', G14.5
	1	/ ' XCal Segment 5:               ', G14.5, T65,
	1	' XCal Segment 5:               ', G14.5
	1	/ ' Right Collimator:             ', G14.5, T65,
	1	' Right Collimator:             ', G14.5 )

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
	1			(ANLG_REC.B_LO_GRT(I),
	1			 ANLG_REC.B_HI_GRT(I), I=1,16)

6009	FORMAT (' External Calibrator:          ', G14.5, t65,
	1	' External Calibrator:          ', G14.5
	1	/ ' Sky Horn:                     ', G14.5, t65,
	1	' Sky Horn:                     ', G14.5
	1	/ ' Reference Horn:               ', G14.5, t65,
	1	' Reference Horn:               ', G14.5
	1	/ ' Internal Reference Source:    ', G14.5, t65,
	1	' Internal Reference Source:    ', G14.5
	1	/ ' Left Dihedral:                ', G14.5, t65,
	1	' Left Dihedral:                ', G14.5
	1	/ ' Bolometer Assembly, RH:       ', G14.5, t65,
	1	' Bolometer Assembly, RH:       ', G14.5
	1	/ ' Bolometer Assembly, RL:       ', G14.5, t65,
	1	' Bolometer Assembly, RL:       ', G14.5
	1	/ ' Bolometer Assembly, LH:       ', G14.5, t65,
	1	' Bolometer Assembly, LH:       ', G14.5
	1	/ ' Bolometer Assembly, LL:       ', G14.5, t65,
	1	' Bolometer Assembly, LL:       ', G14.5
	1	/ ' Left Mirror:                  ', G14.5, t65,
	1	' Left Mirror:                  ', G14.5
	1	/ ' Calibrator Resistor 1:        ', G14.5, t65,
	1	' Calibrator Resistor 1:        ', G14.5
	1	/ ' Calibrator Resistor 2:        ', G14.5, t65,
	1	' Calibrator Resistor 2:        ', G14.5
	1	/ ' Calibrator Resistor 3:        ', G14.5, t65,
	1	' Calibrator Resistor 3:        ', G14.5
	1	/ ' Calibrator Resistor 4:        ', G14.5, t65,
	1	' Calibrator Resistor 4:        ', G14.5
	1	/ ' XCal Segment 6:               ', G14.5, T65,
	1	' XCal Segment 6:               ', G14.5
	1	/ ' Left Collimator (Not used):   ', G14.5, t65,
	1	' Left Collimator (Not used):   ', G14.5 ) 

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Temperature Controllers
C
	WRITE (WRITE_LUN, 6011, IOSTAT=STATUS)
	1			(ANLG_REC.TEMP_CTRL(I),
	1			 ANLG_REC.TEMP_CTRL(I+4), I=1,4)

6011	FORMAT (/'  *****  Temperature Controllers Temperatures  '
	1	,'*****' 
	1	/' Internal Reference Source, Side A:   ', G14.5, t65,
	1	' Internal Reference Source, Side B:   ', G14.5
	1	/' Reference Horn, Side A:              ', G14.5, t65,
	1	' Reference Horn, Side B:              ', G14.5
	1	/' Sky Horn, Side A:                    ', G14.5, t65,
	1	' Sky Horn, Side B:                    ', G14.5
	1	/' External Calibrator, Side A:         ', G14.5, t65,
	1	' External Calibrator, Side B:         ', G14.5)

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C Other temperatures/currents
C
	WRITE (WRITE_LUN, 6012, IOSTAT=STATUS)

6012	FORMAT (/' ***** OTHER TEMPERATURES/CURRENTS *****' )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	WRITE (WRITE_LUN, 6013, IOSTAT=STATUS)
	1			ANLG_REC.IPDU_TEMP,
	1			ANLG_REC.CNA_TEMP(1),
	1			ANLG_REC.CNA_TEMP(3),
	1			ANLG_REC.CNA_TEMP(2),
	1			ANLG_REC.CNA_TEMP(4),
	1			ANLG_REC.DBX_TEMP,
	1			ANLG_REC.STAT_MON_TEMP,
	1			ANLG_REC.PAMP_CHAN,
	1			ANLG_REC.PAMP_OP,
	1			ANLG_REC.HOT_SPOT,
	1			ANLG_REC.MTM_CAL_MTR,
	1			ANLG_REC.MTM_POS,
	1			ANLG_REC.BOL_VOLT,
	1			ENG_REC.EN_TAIL.LMAC_ANALOG_TEMP,
	1			ENG_REC.EN_TAIL.LMAC_DIGITAL_TEMP

6013	FORMAT (' IPDU temperature - Side A:       ', G14.5, t65,
	1	' IPDU temperature - Side B:       ', G14.5
	1	/' Channel temperature - RH:        ', G14.5, t65,
	1	' Channel temperature - LH:        ', G14.5
	1	/' Channel temperature - RL:        ', G14.5, t65,
	1	' Channel temperature - LL:        ', G14.5
	1	/' Drive Box temperature - Side A:  ', G14.5, t65,
	1	' Drive Box temperature - Side B:  ', G14.5
	1	/' Status Monitor temp. - Side A:   ', G14.5, t65,
	1	' Status Monitor temp. - Side B:   ', G14.5
	1	/' Channel Pre-Amp temperature:     ', G14.5, t65,
	1	' Optical Pre-Amp temperature:     ', G14.5
	1	/' Hot Spot Heater Current-Side A:  ', G14.5, t65,
	1	' Hot Spot Heater Current-Side B:  ', G14.5
	1	/' MTM Cal Motor - Side A:          ', G14.5, t65,
	1	' MTM Cal Motor - Side B:          ', G14.5
	1	/' MTM Position - Side A:           ', G14.5, t65,
	1	' MTM Position - Side B:           ', G14.5
	1	/' Bolometer Bias Voltage, RH:      ', G14.5, t65,
	1	' Bolometer Bias Voltage, RL:      ', G14.5
	1	/' Bolometer Bias Voltage, LH:      ', G14.5, t65,
	1	' Bolometer Bias Voltage, LL:      ', G14.5
	1	/' LMAC analog temperature:         ', G14.5, t65,
	1	' LMAC digital temperature:        ', G14.5 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C IPDU Voltages
C
	WRITE (WRITE_LUN, 6014, IOSTAT=STATUS)	ANLG_REC.IPDU_VOLT

6014	FORMAT ( /' ***** IPDU VOLTAGES *****' 
	1	/' Digital Converter, -15V, Side A: ', G14.5, t65,
	1	' Digital Converter, -15V, Side B: ', G14.5
	1	/' Digital Converter, +15V, Side A: ', G14.5, t65,
	1	' Digital Converter, +15V, Side B: ', G14.5
	1	/' Digital Converter, +5V, Side A:  ', G14.5, t65,
	1	' Digital Converter, +5V, Side B:  ', G14.5
	1	/' Analog Converter, +15V, Side A:  ', G14.5, t65,
	1	' Analog Converter, +15V, Side B:  ', G14.5
	1	/' Analog Converter, -15V, Side A:  ', G14.5, t65,
	1	' Analog Converter, -15V, Side B:  ', G14.5
	1	/' Bias Pre-Reg. +25V, Side A:      ', G14.5, t65,
	1	' Bias Pre-Reg. +25V, Side B:      ', G14.5
	1	/' Int. P.S., +28V, Side A:         ', G14.5, t65,
	1	' Int. P.S., +28V, Side B:         ', G14.5
	1	/' Int. P.S., +15V, Side A:         ', G14.5, t65,
	1	' Int. P.S., +15V, Side B:         ', G14.5
	1	/' Int. P.S., -15V, Side A:         ', G14.5, t65,
	1	' Int. P.S., -15V, Side B:         ', G14.5
	1	/' Int. P.S., +5V, Side A:          ', G14.5, t65,
	1	' Int. P.S., +5V, Side B:          ', G14.5 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF

C
C IPDU Currents
C
	WRITE (WRITE_LUN, 6015, IOSTAT=STATUS)	ANLG_REC.IPDU_AMP

6015	FORMAT ( /' ***** IPDU CURRENTS *****' 
	1	/' Bias Pre-Reg. Current, Side A:   ', G14.5, t65,
	1	' Bias Pre-Reg. Current, Side B:   ', G14.5
	1	/' Analog Conv. Current, Side A:    ', G14.5, t65,
	1	' Analog Conv. Current, Side B:    ', G14.5
	1	/' Digital Conv. Current, Side A:   ', G14.5, t65,
	1	' Digital Conv. Current, Side B:   ', G14.5
	1	/' Constant Current RH:             ', G14.5, t65,
	1	' Constant Current LH:             ', G14.5
	1	/' Constant Current RL:             ', G14.5, t65,
	1	' Constant Current LL:             ', G14.5
	1	/' Int. Conv. Current, Side A:      ', G14.5, t65,
	1	' Int. Conv. Current, Side B:      ', G14.5 )

	IF (STATUS .NE. 0) THEN
	  CALL LIB$SIGNAL ( FGA_WRITE_ERR, %VAL(1), %VAL(STATUS) )
	END IF


	RETURN
	END
