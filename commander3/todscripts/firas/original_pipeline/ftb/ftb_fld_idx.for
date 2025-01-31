C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Fld_Idx ( Field, Index )

C-------------------------------------------------------------------------------
C
C	Purpose: To determine the index of the selected engineering field
C		 in the FUT_Plot_Names array.
C
C	Author: Shirley M. Read
C		STX, January, 1990
C
C	Invocation: Status = FTB_Fld_Idx ( Field, Index )
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
C	  Field          C*30           FIRAS engineering field
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Index          I*2            Index of field name in array
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
C	  Set function return             
C	  Search the Plot_Names array for a match to the selected engineering
C	    field name.
C	  If there is no match, signal an error and set function return status.
C	  Return.
C
C------------------------------------------------------------------------------
C	
	Implicit None

	Include		'(Fut_Plot_Names)'

C	Passed Parameters.

	Character*30 Field           ! Engineering field name
	Integer*2    Index           ! Index in the Plot_Names array

C       Functions

	Integer*4 Lib$Locc

C	Local Declarations.

	Integer*4 Pos, Pos1          ! Line position
	Logical*1 Match /.False./    ! Matched field name to plot name
	Integer*2 Ix

C	External Parameters.

	External FTB_Normal
	External FTB_Aberr
	External FTB_NoMatch

	FTB_Fld_Idx = %loc(FTB_Normal)

C	Search for a match of field names.

	Pos = Lib$Locc(' ',Field )
	Pos1 = Pos - 1
	Do Ix = 1,118
	   If ( Field(1:pos1) .EQ. Plot_Name(Ix) ) Then
	     Index = Ix
	     Match = .True.
	   Endif
	Enddo 

	If ( .NOT. Match ) Then
	   Call Lib$Signal(FTB_NoMatch)
	   FTB_Fld_Idx = %loc(FTB_Aberr)
	Endif

	Return
	End
