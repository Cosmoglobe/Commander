      Integer*4 Function FIT_Exit ( )

C------------------------------------------------------------------------
C    PURPOSE: Close up and exit FIT_Instrument_Trends.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX			revision: 1987 Nov 20
C            April 24, 1987			to comment out UEND
C
C    INVOCATION: status = FIT_Exit ( )
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS: 
C	STATUS		I*4		Status flag.
C
C    SUBROUTINES CALLED: 
C	FUT_DeaDev
C	LIB$Delete_Logical
C	LIB$Free_LUN
C    C		UEnd
C
C    COMMON VARIABLES USED: 
C	FCC_INTERACTIVE		I*4		Interactive/Batch flag.
C	FCC_REPORT		I*4		Report flag.
C	FCC_REPORT_LUN		I*4		Report file lun.
C	FCC_REPORT_FILE		CH*64		Report file name.
C	FCC_PLOTS		I*4		Plots flag.
C
C    INCLUDE FILES: 
C	FIT_Invoc
C	FUT_Error
C	FUT_Params
C	$SSDef
C
C----------------------------------------------------------------------

	Implicit        None

	Include         '(FIT_Invoc)'
	Include         '(FUT_Error)'
	Include         '(FUT_Params)'
	Include		'($SSDef)'

	Integer		*4	status
	Integer		*2	istat

	Logical		*4      FR_Close_RDF
	Integer         *4      LIB$Free_LUN

	External        FIT_Normal
	External        FIT_CloseRep
	External        FIT_CFRClose

C
C Close the report file.
C
	If (fcc_report .Eq. fac_present) Then

	   Close ( Unit=fcc_report_lun, IOstat=istat )

	   If (istat .Ne. 0) Then
	      FIT_Exit = %Loc(FIT_CloseRep)
	      Call LIB$Signal(FIT_CloseRep,%Val(1),%Val(istat))
	   End If

	   status = LIB$Free_LUN ( fcc_report_lun )

	   If (status .Ne. SS$_Normal) Then
	      Call LIB$Signal(%Val(status))
	   End If

	End If

c
c Shut down Field Retriever.
c
	status = FR_Close_RDF (0,0)

	If (.Not. status) Then
	   Call LIB$Signal(FIT_CFRClose)
	End If


	FIT_Exit = %Loc(FIT_Normal)

	Return
	End
