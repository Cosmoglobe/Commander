      Integer*4 Function FIT_Init ( fr_lun )

C------------------------------------------------------------------------
C    PURPOSE: Initialize program FIT_Instrument_Trends.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Fred Shuman
C            STX			revision: 1987 Nov 20
C            September 11, 1987			to comment out USTART
C
C    INVOCATION: status = FIT_Init ( )
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS:
C	STATUS		I*4		Status flag.
C
C    SUBROUTINES CALLED:
C	FIT_Get_Command
C	LIB$Set_Logical
C	LIB$Get_LUN
C
C    COMMON VARIABLES USED:
C	FCC_INTERACTIVE		 I*4		Interactive/Batch flag.
C	FCC_REPORT		 I*4		Report flag.
C	FCC_REPORT_LUN		 I*4		Report file lun.
C	FCC_REPORT_FILE		Ch*64		Report file name.
C	FCC_PLOTS		 I*4		Plots flag.
C	FCC_PLOTS_FILE		Ch*64		Plots file name.
C
C    INCLUDE FILES:
C	FIT_Invoc
C	FUT_Error
C	FUT_Params
C	$SSDEF
C
C----------------------------------------------------------------------
C   Changes:
C
C	VERSION 4.4.1 SER 3306 Q. CHUNG STX 08/10/89
C                    PROVIDE VERSION NUMBER TO TRACK SOFTWARE UPDATE.
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FIT_Invoc)'
	Include		'(FUT_Error)'
	Include		'(FUT_Params)'
	Include		'($SSDEF)'

	Integer		*4	fr_lun
	Integer		*4	status
	Integer		*2	istat

	Logical		*4	FR_Init
	Integer		*4	FIT_Get_Command
	Integer		*4	LIB$Get_LUN

	Integer		*4	CUT_Register_Version
	Integer		*4	CUT_Display_Banner
	Integer		*4	num_vol/80/
	Character	*5	version
	Parameter		(version='4.4.1')

	External        FIT_Normal
	External        FIT_OpenRep
	External        FIT_CFRInit

C
C Parse the input command line.
C
	status = FIT_Get_Command ( )

	If (status .Ne. %Loc(FIT_Normal)) Then
	   FIT_Init = status
	   Return
	End If

C
C Open the report file.
C
	If (fcc_report .Eq. fac_present) Then

	   status = LIB$Get_LUN ( fcc_report_lun )

	   If (status .Ne. SS$_Normal) Then
	      FIT_Init = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   Open ( UNIT=fcc_report_lun, FILE=fcc_report_file,
	2         STATUS='New', IOSTAT=istat )

	   If (istat .Ne. 0) Then
	      FIT_Init = %Loc(FIT_OpenRep)
	      Call LIB$Signal(FIT_OpenRep, %Val(2),%Val(istat), fcc_report_file)
	      Return
	   End If

	   status = CUT_Register_Version(version)
	   status = CUT_Display_Banner(fcc_report_lun, num_vol,
	2                              ' FIRAS Facility FIT_Instrument_Trends')
	End If

C
C Initialize the Field Retriever.
C
	status = FR_Init ( fr_lun )

	If (.Not. status) Then
	   FIT_Init = %Loc(FIT_CFRInit)
	   Call LIB$Signal(FIT_CFRInit)
	   Return
	End If


	FIT_Init = %Loc(FIT_Normal)

	Return
	End
