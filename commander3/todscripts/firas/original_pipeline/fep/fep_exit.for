      Integer*4 Function FEP_Exit ( )

C------------------------------------------------------------------------
C    PURPOSE: Close up and exit FEP_Engplots.
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
C    INVOCATION: status = FEP_Exit ( )
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS: 
C	STATUS		I*4		Status flag.
C
C    SUBROUTINES CALLED: 
C	LIB$Delete_Logical
C	LIB$Free_LUN
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
C	$SSDef
C
C----------------------------------------------------------------------

	Implicit        None

	Include         '(FEP_Invoc)'
	Include         '(FUT_Error)'
	Include         '(FUT_Params)'
	Include		'($SSDef)'

	Integer		*4	status
	Integer		*2	istat

	Logical		*4	FR_Close_RDF
	Integer         *4      LIB$Delete_Logical
	Integer         *4      LIB$Free_LUN

	External        FEP_Normal
	External        FEP_CloseRep
	External	FEP_CFRClose

C
C Close the report file.
C
	If (fcc_report .Eq. fac_present) Then

	   Close ( Unit=fcc_report_lun, IOstat=istat )

	   If (istat .Ne. 0) Then
	      FEP_Exit = %Loc(FEP_CloseRep)
	      Call LIB$Signal(FEP_CloseRep,%Val(1),%Val(istat))
	   End If

	   status = LIB$Free_LUN ( fcc_report_lun )

	   If (status .Ne. SS$_Normal) Then
	      Call LIB$Signal(%Val(status))
	   End If

	End If

c
c End the Field Retriever session.
c
!	Call FR_End_Open

c
c Shut down Field Retriever.
c
	status = FR_Close_RDF (0,0)

	If (.Not. status) Then
	   Call LIB$Signal(FEP_CFRClose)
	End If

	FEP_Exit = %Loc(FEP_Normal)

	Return
	End
