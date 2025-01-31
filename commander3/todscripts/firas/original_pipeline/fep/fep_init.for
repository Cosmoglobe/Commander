      Integer*4 Function FEP_Init ( fr_lun )

C------------------------------------------------------------------------
C    PURPOSE: Initialize program FEP_Engplots.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX
C            April 24, 1987
C
C    INVOCATION: status = FEP_Init ( )
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS: 
C	STATUS		I*4		Status flag.
C
C    SUBROUTINES CALLED: 
C	FEP_Get_Command
C	LIB$Set_Logical
C	LIB$Get_LUN
C	UStart
C
C    COMMON VARIABLES USED: 
C	FCC_INTERACTIVE		I*4		Interactive/Batch flag.
C	FCC_REPORT		I*4		Report flag.
C	FCC_REPORT_LUN		I*4		Report file lun.
C	FCC_REPORT_FILE		CH*64		Report file name.
C
C    INCLUDE FILES: 
C	FEP_Invoc
C	FUT_Error
C	FUT_Params
C	$SSDEF
C
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FEP_Invoc)'
	Include		'(FUT_Error)'
	Include		'(FUT_Params)'
	Include		'($SSDEF)'

	Integer         *4	status
	Integer         *4	fr_lun
	Integer         *2	istat

	Logical		*4	FR_Init
	Integer         *4	FEP_Get_Command
	Integer         *4	LIB$Set_Logical
	Integer         *4	LIB$Get_LUN

	External        FEP_Normal
	External        FEP_OpenRep
	External        FEP_CFRInit

C
C Parse the input command line.
C
	status = FEP_Get_Command ( )

	If (status .Ne. %Loc(FEP_Normal)) Then
	   FEP_Init = status
	   Return
	End If

C
C Open the report file.
C
	If (fcc_report .Eq. fac_present) Then

	   status = LIB$Get_LUN ( fcc_report_lun )

	   If (status .Ne. SS$_Normal) Then
	      FEP_Init = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   Open ( Unit=fcc_report_lun, File=fcc_report_file,
     .            Status='New', IOstat=istat )

	   If (istat .Ne. 0) Then
	      FEP_Init = %Loc(FEP_OpenRep)
	      Call LIB$Signal(FEP_OpenRep,%Val(2),%Val(istat),
     .				fcc_report_file)
	      Return
	   End If

	End If

C
C Initialize the Field Retriever.
C
	status = FR_Init ( fr_lun )

	If (.Not. status) Then
	   FEP_Init = %Loc(FEP_CFRInit)
	   Call LIB$Signal(FEP_CFRInit)
	   Return
	End If

	FEP_Init = %Loc(FEP_Normal)

	Return
	End
