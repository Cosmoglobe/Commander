	Integer*4 Function FUT_Free_LUN ( lun )

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

	External	FUT_Normal
	External	FUT_InvLUN
	External	FUT_LUNFree

	Common /FUT_LUN/ free_lun

	If (lun .Ge. 21 .And. lun .Le. 99) Then

	  If (free_lun(lun)) Then
	    Call LIB$Signal(FUT_LUNFree)
	    FUT_Free_LUN = %Loc(FUT_LUNFree)
	    Return
	  Else
	    free_lun(lun) = .True.
	  End If

	Else
	  FUT_Free_LUN = %Loc(FUT_InvLUN)
	  Call LIB$Signal(FUT_InvLUN)
	  Return
	End If

	FUT_Free_LUN = %Loc(FUT_Normal)

	Return
	End
