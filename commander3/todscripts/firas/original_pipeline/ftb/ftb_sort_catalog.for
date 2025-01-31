C------------------------------------------------------------------------

      Integer*4 Function FTB_Sort_Catalog ( Cats, NRecs )

C------------------------------------------------------------------------
C    PURPOSE: Sorts the catalog records into ascending binary initial
C	      time (primary key) and activity time (secondary key).
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            ST Systems Corporation
C            July 31, 1989
C
C    INVOCATION: 
C
C    INPUT PARAMETERS:
C
C    OUTPUT PARAMETERS: 
C
C    SUBROUTINES CALLED: 
C
C    COMMON VARIABLES USED:
C
C    INCLUDE FILES: 
C
C----------------------------------------------------------------------

	Implicit None

	Include 'CT$Library:CTUser.Inc'
	Include	'($SSDef)'
	Include	'(CUT_Params)'

	Integer		*2	NRecs

	Dictionary 'CCM_CME_Catalog_Entry'
	Record /CCM_CME_Catalog_Entry/ Cats(CUT$_Max_Segments)

	Integer		*2 	Recl / 146 /
	Character	*146	SortRec

	Integer		*4	i
	Integer		*4	status
	Integer		*2	Sort_Keys(9)

	Integer		*4	SOR$Begin_Sort
	Integer		*4	SOR$Release_Rec
	Integer		*4	SOR$Sort_Merge
	Integer		*4	SOR$Return_Rec
	Integer		*4	SOR$End_Sort

	External	DSC$K_DType_Q
	External	DSC$K_DType_B
	External	FTB_Normal

C
C Initialize the sorting process.
C
	Sort_Keys(1) = 2
	Sort_Keys(2) = %Loc(DSC$K_DType_Q)
	Sort_Keys(3) = 0
	Sort_Keys(4) = 52
	Sort_Keys(5) = 8
	Sort_Keys(6) = %Loc(DSC$K_DType_Q)
	Sort_Keys(7) = 0
	Sort_Keys(8) = 127
	Sort_Keys(9) = 8

	status = SOR$Begin_Sort ( Sort_Keys, Recl )

	If (status .Ne. SS$_Normal) Then
          FTB_Sort_Catalog = status
          Call LIB$Signal ( %Val(status) )
	  Return
	End If

C
C Introduce records to the sort.
C
	Do i=1,NRecs

	   Call LIB$Movc3 ( 146, Cats(i), %Ref(SortRec) )

	   status = SOR$Release_Rec ( SortRec )

	   If (status .Ne. SS$_Normal) Then
             FTB_Sort_Catalog = status
             Call LIB$Signal ( %Val(status) )
	     Return
	   End If

	End Do

C
C Perform the sort.
C
	status = SOR$Sort_Merge ( )

	If (status .Ne. SS$_Normal) Then
          FTB_Sort_Catalog = status
          Call LIB$Signal ( %Val(status) )
	  Return
	End If

C
C Fetch records from the sort.
C
	Do i=1,NRecs

	   status = SOR$Return_Rec ( SortRec )

	   Call LIB$Movc3 ( 146, %Ref(SortRec), Cats(i) )

	   If (status .Ne. SS$_Normal) Then
             FTB_Sort_Catalog = status
             Call LIB$Signal ( %Val(status) )
	     Return
	   End If

	End Do

C
C Finish the sort.
C
	status = SOR$End_Sort ( )

	If (status .Ne. SS$_Normal) Then
           FTB_Sort_Catalog = status
           Call LIB$Signal ( %Val(status) )
	   Return
	End If


	FTB_Sort_Catalog = %Loc(FTB_Normal)

	Return
	End
