C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Get_Fldval (Eng_Rec, Index, Fieldval, Goodval)

C-------------------------------------------------------------------------------
C
C	Purpose: To get the field value of the selected engineering field.
C
C	Author: Shirley M. Read
C		STX, January, 1990
C
C	Invocation: Status = FTB_Get_Fldval ( Eng_Rec, Index, Filedval, Goodval)
C
CH	Change Log:
CH
C	  ----------------------------------------------------------------------
C
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Eng_Rec       FDQ_ENG         FIRAS engineering record
C	  Index         I*2             Index of field name in array
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C         Fieldval      R*4             Value of engineering field
C	  Goodval       L*1             Flag for good or bad value
C	
C	Subroutines Called:
C
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C
C	Processing Method:
C	  Set function return.
C	  Set Goodval to true.
C	  Using the index to the plot names array, get the field value.
C	  If the value is the flag value, set the goodval flag to false.
C	  Return.
C
C------------------------------------------------------------------------------
C	
	Implicit None

	Include		'(Fut_Plot_Names)'

C	Passed Parameters.

	Dictionary 'FDQ_ENG'
	Record /FDQ_ENG/ Eng_Rec     ! Current engineering record
	Integer*2    Index           ! Index in the Plot_Names array
	Real*4       Fieldval	     ! Value of engineering field
	Logical*1    Goodval	     ! Flag indicating good value
	
C	Local Declarations.

	Integer*2 Ix
	Real*4 FV / - 9999.0 /       ! Flag value for invalid data

C	External Parameters.

	External FTB_Normal
	External FTB_Aberr
	External FTB_Badfld

C 	Set status to normal  and goodval flag to true.

	FTB_Get_Fldval = %loc(FTB_Normal)
	Goodval = .True.

C	Get the value of the engineering field using the Index.

	If ( Index .LE. 64 ) Then
	  Fieldval = Eng_Rec.En_Analog.GRT(Index)

	Elseif ( Index .EQ. 65 ) Then
	  Fieldval = Eng_Rec.En_Analog.Ipdu_Temp(1)
	Elseif ( Index .EQ. 66 ) Then
	  Fieldval = Eng_Rec.En_Analog.Ipdu_Temp(2)

	Elseif ( Index .EQ. 67 ) Then
	  Fieldval = Eng_Rec.En_Analog.CNA_Temp(1)
	Elseif ( Index .EQ. 68 ) Then
	  Fieldval = Eng_Rec.En_Analog.CNA_Temp(2)
	Elseif ( Index .EQ. 69 ) Then
	  Fieldval = Eng_Rec.En_Analog.CNA_Temp(3)
	Elseif ( Index .EQ. 70 ) Then
	  Fieldval = Eng_Rec.En_Analog.CNA_Temp(4)

	Elseif ( Index .EQ. 71 ) Then
	  Fieldval = Eng_Rec.En_Analog.Dbx_Temp(1)
	Elseif ( Index .EQ. 72 ) Then
	  Fieldval = Eng_Rec.En_Analog.Dbx_Temp(2)

	Elseif ( Index .EQ. 73 ) Then
	  Fieldval = Eng_Rec.En_Analog.Stat_Mon_Temp(1)
	Elseif ( Index .EQ. 74 ) Then
	  Fieldval = Eng_Rec.En_Analog.Stat_Mon_Temp(2)

	Elseif ( Index .EQ. 75 ) Then
	  Fieldval = Eng_Rec.En_Analog.Pamp_Chan
	
	Elseif ( Index .EQ. 76 ) Then
	  Fieldval = Eng_Rec.En_Analog.Pamp_Op

	Elseif ( Index .EQ. 77 ) Then
	  Fieldval = Eng_Rec.En_Analog.Hot_Spot(1)
	Elseif ( Index .EQ. 78 ) Then
	  Fieldval = Eng_Rec.En_Analog.Hot_Spot(2)

	Elseif ( Index .EQ. 79 ) Then
	  Fieldval = Eng_Rec.En_Analog.MTM_Cal_Mtr(1)
	Elseif ( Index .EQ. 80 ) Then
	  Fieldval = Eng_Rec.En_Analog.MTM_Cal_Mtr(2)

	Elseif ( Index .EQ. 81 ) Then
	  Fieldval = Eng_Rec.En_Analog.bol_volt(1)
	Elseif ( Index .EQ. 82 ) Then
	  Fieldval = Eng_Rec.En_Analog.bol_volt(2)
	Elseif ( Index .EQ. 83 ) Then
	  Fieldval = Eng_Rec.En_Analog.bol_volt(3)
	Elseif ( Index .EQ. 84) Then
	  Fieldval = Eng_Rec.En_Analog.bol_volt(4)

	Elseif ( (Index .GE. 85) .AND. (Index .LE. 104) ) Then
	  Fieldval = Eng_Rec.En_Analog.IPDU_Volt(Index - 84)

	Elseif ( (Index .GE. 105) .AND. (Index .LE. 116) ) Then
	  Fieldval = Eng_Rec.En_Analog.IPDU_Amp(Index - 104)

	Else
	  Call Lib$Signal(FTB_BadFld)
	  FTB_Get_Fldval = %loc(FTB_Aberr)
	Endif

C	Check to see if the engineering field has been flagged as bad data.

	If ((FTB_Get_Fldval .EQ. %loc(FTB_Normal)) .AND. 
	1	( Fieldval .EQ. FV )) Goodval = .False.

	Return
	End
