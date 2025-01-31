C-----------------------------------------------------------------------------
	Program FRD_GFG
C-----------------------------------------------------------------------------
C 	Program Description:
C	    The purpose of this facility is to create reference files,
C	  FEX_GAIN.DAT and FEX_FAKEIT.DAT.  These files contain records
C	  of gain settings and fake_it statuses respectively.  Each record
C         will have times indicating the range over which the value was set.
C         Each record ends due to value changes, gaps in the housekeeping
C         data, housekeeping frames which had bad telemetry quality, or the
C         end of data is reached.  The FPP facility will use these "time line"
C         reference files in determining the gains and fakeit statuses 
C         associated with each science record, and occurences of data gaps.
C
C 	Calling Sequence:
C 	  This is a main program.
C 
C       Programmer:  Nilo G. Gonzales, Larry Rosen/STX
C                    April 5, 1991    
C       Change Log:
CH
CH
C 	Input Files:
C 	  FIRAS Raw Housekeeping Data Archive Files : NFS_HKP
C       
C 	Output Files:
C 	  FIRAS Reference Archive Files : FEX_GAIN.DAT, FEX_FAKEIT.DAT 
C 	  Report file for run : FRD_GFG.REP_yydddhhmm 
C 
C 	Include Files Used:
C 	  CT$Library:CTUser.Inc (COBETRIEVE return status definitions)
C	  '($SSDef)'
C 
C 	Subroutines and Functions Called:
C 	  FRD_GFG_Get_Options
C	  FRD_Gfg_Open
C 	  FRD_Gfg_Check
C	  FRD_Gfg_Close
C	  Lib$Establish
C	  Lib$Signal
C 
C 	Method Used:  PDL for FRD_GFG
C
C    Begin  
C 	  Establish the condition handler.
C
C	  Call FRD_GFG_GET_OPTIONS to parse user options from the command line.
C
C         Call FRD_GFG_OPEN to open housekeeping record in FIRAS archive
C              for read, open FEX_GAIN and FEX_FAKEIT for write.
C              If a report is requested, Open the report file, then write
C              initial information to report file.
C
C         Call FRD_GFG_CHECK, This subroutine checks for telemetry quality,
C              gain change, fakeit change, time gaps.  It also calls two
C              subroutines to get the fakeit and gain values.
C       
C         Call FRD_GFG_CLOSE to close FIRAS archive
C      End
C----------------------------------------------------------------------------- 
	Implicit	None

C	Include Files

	Include		'CT$Library:CTUser.Inc'
	Include		'($SSDef)'
	Include		'(Fut_Error)'

C	Functions

 	Integer*4 FRD_GFG_Get_Options
 	Integer*4 FRD_Gfg_Open
 	Integer*4 FRD_Gfg_Check
 	Integer*4 FRD_Gfg_Close

C	Externals

	External  FUT_Error    ! Condition handler
	External  FRD_Normal
	External  FRD_CTHkpEof
	External  frd_aberr

C	Local Variables

	Integer*4 	Status		! Status of processing
	Integer*4 	Retstat		! Return status from function call
	Integer*4 	Success / 1 /, Err / 2 /  ! Values for status
	Integer*4       Max_Sec         ! Maximum number of seconds allowed
                                        ! between MJF before a gap is declared
	Integer*4	LUN_HKP  	! Unit numbers for NFS_HKP CT archive
	Integer*4	LUN_R_Gain	! Unit numbers for FEX_GAIN CT archive
	Integer*4	LUN_L_Gain	! Unit numbers for FEX_GAIN CT archive
	Integer*4       LUN_R_Fakeit    ! Unit numbers for FEX_FAKEIT archive 
	Integer*4       LUN_L_Fakeit    ! Unit numbers for FEX_FAKEIT archive 
	Integer*4	LUN_Rpt         ! Report unit number
	Character*21    Report_File     ! Report file name
	Logical*1       Report          ! Flag to enable writing a report file
	Character*79	CMD_Line(2)	! Command line string with defaults
	Character*14	GMT_Start	! Start time for running GFG
	Character*14	GMT_Stop	! Stop time for running GFG

C	Set status for FRD processing to Success.

	Status = Success

C       Establish condition handler.       

	Call Lib$Establish ( FUT_ERROR )

C       Get processing options from command line.

	Retstat = FRD_GFG_Get_Options ( Max_Sec, Report, Report_File, CMD_Line,
	1   GMT_Start, GMT_Stop)
	If ( Retstat .ne. %loc(FRD_Normal) ) then
	  Status = Err
	Endif

C	Open housekeeping record (NFS_SDF) in the FIRAS archive for read,
C       Open FEX_GAIN and FEX_FAKEIT for write.  If report is requested,
C       Open the report file, then write initial information to report file.

	If ( Status .eq. Success ) Then
	   Retstat = FRD_Gfg_Open ( LUN_HKP, LUN_R_Gain, LUN_L_Gain,
 	1       LUN_R_Fakeit, LUN_L_Fakeit, Report_File, LUN_Rpt, Report,
	2        CMD_Line, GMT_Start, GMT_Stop)
	      If ( Retstat .ne. %loc(FRD_Normal) ) Then
	        Status = Err
	      Endif
	Endif

	If ( Status .eq. Success ) Then
	   Retstat = FRD_Gfg_Check ( LUN_HKP, LUN_R_Gain, LUN_L_Gain,
	1            LUN_R_Fakeit, LUN_L_Fakeit, LUN_Rpt, Report, Max_Sec )
	   If ( Retstat .ne. %loc(FRD_CTHkpEof) ) Then
	      Status = Err
	   Endif
	Endif

	If ( Status .eq. Success ) Then
	   Retstat = FRD_Gfg_Close ( LUN_HKP, LUN_R_Gain, LUN_L_Gain,
	1            LUN_R_Fakeit, LUN_L_Fakeit, LUN_Rpt, Report )
	      If ( Retstat .Ne. %loc(FRD_Normal) ) Then
	        Status = Err
	      Endif
	Endif

C  Signal the processing status.

	IF (status .EQ. success) THEN
	  CALL lib$signal (frd_normal)
	ELSE
	  CALL lib$signal (frd_aberr)
	ENDIF


C  Turn FUT_Report_Lun to off.  Close report.
 
	If (Report) Then
	   FUT_Report_Lun = 0
	   Close (Lun_RPT)
	Endif

C  Exit with status.

	IF (status .EQ. success) THEN
	  CALL exit (ss$_normal)
	ELSE
	  CALL exit (ss$_abort)
	ENDIF

	End
