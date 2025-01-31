	Integer*4 Function FUT_Get_LUN ( lun )

C------------------------------------------------------------------------
C    PURPOSE:
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            ST Systems Corporation
C            August 3, 1988
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

	Implicit	None

	Integer		*4	lun
	Logical		*1	free_lun(21:99)
	Logical		*1	first_time/.True./

	External	FUT_Normal
	External	FUT_NoLUN

	Common /FUT_LUN/ free_lun

	If (first_time) Then
	  Do lun=21,99
	    free_lun(lun) = .True.
	  End Do
	  first_time = .False.
	End If

	lun = 99

	Do While (.Not. free_lun(lun) .And. lun .Ge. 21)
	  lun = lun - 1
	End Do

	If (lun .Ge. 21 .And. lun .Le. 99) Then
	  free_lun(lun) = .False.
	Else
	  lun = 0
	  FUT_Get_LUN = %Loc(FUT_NoLUN)
	  Call LIB$Signal(FUT_NoLUN)
	  Return
	End If

	FUT_Get_LUN = %Loc(FUT_Normal)

	Return
	End
