C------------------------------------------------------------------------

	Program FTB_Overlap

C------------------------------------------------------------------------
C    PURPOSE: Produce a list of non-overlapping timeranges for the NFS_SDF
C	      raw science data on which FPP is to be run.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            ST Systems Corporation
C            October 13, 1989
C
C    CHANGE:
C       SPR 7020, FTB_OVERLAP CT read error during the period: 90175-90179
C                 H. Wang, STX, 7/2/90
C
C       SPR 7346, FTB_OVERLAP CT should take the overlap case:
C                 |------------| < A segment delivered after B
C                 |-----------------| < B segment
C                 |------------||---| < create new 2 segments
C                 H. Wang, STX, 8/22/90
C 
C    INVOCATION: 
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS: None
C
C    SUBROUTINES CALLED: 
C	FTB_Overlap_Search
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES:
C	$SSDef
C
C----------------------------------------------------------------------

	Implicit None

	Include		'($SSDef)'

	Integer		*4	jstart(2)
	Integer		*4	jstop(2)
	Logical		*1	report
	Integer		*4	report_lun

	Integer		*4	status
	Character	*14	start
	Character	*14	stop

	Integer		*4	num_vol/80/
        Integer		*4	lun_out/6/
	Character	*5	version
	Parameter	(version='6.8')

	Integer		*4	CUT_Register_Version
	Integer		*4	CUT_Display_Banner
	Integer		*4	FTB_Overlap_Init
	Integer		*4	FTB_Overlap_Search

	External		LIB$Establish
	External		FUT_Error
	External		FTB_Normal
	External		FTB_Aberr

C Initialize OVERLAP.

	Call LIB$Establish ( FUT_Error )

	status = CUT_Register_Version(version)
	status = CUT_Display_Banner(lun_out,num_vol,
	1			'FIRAS Facility FTB_Overlap')

	status = FTB_Overlap_Init ( jstart, jstop, report, report_lun )

	Call CT_Binary_To_GMT ( jstart, start )
	Call CT_Binary_To_GMT ( jstop, stop )

	Write (6,100) start, stop

	If (report) Then
	  Write (report_lun, 100) start, stop
	End If

C Produce a list of timeranges on which FPP can be run such that the
C resulting FPP_SDF raw science data will not overlap.

	If (status) Then
	   status = FTB_Overlap_Search ( jstart, jstop, report, report_lun )
	End If

	If (status) Then
	   If (report) Then
	      Write (report_lun, 200)
	   End If
	   Call LIB$Signal(FTB_Normal)
	   Call Exit(SS$_Normal)
	Else
	   If (report) Then
	      Write (report_lun, 300)
	   End If
	   Call LIB$Signal(FTB_Aberr)
	   Call Exit(SS$_Abort)
	End If


100	Format (//, x, 'FTB_Overlap processing timerange: ', a, ' to ', a)
200	Format (//, x, 'FTB_Overlap completes successfully.')
300	Format (//, x, 'FTB_Overlap terminates abnormally.')

	Stop
	End
