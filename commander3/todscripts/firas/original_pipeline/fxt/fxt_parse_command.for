C-------------------------------------------------------------------------------

	Integer*4 Function FXT_Parse_Command ( Gmt_Start, Gmt_Stop,
	1         File_Seg, Filter, Init_Xtrm, Plot, Min_Offtime, 
	2	  Max_Offtime, Flags, Report, Report_File, CurGMT)

C-------------------------------------------------------------------------------
C
C	Purpose: To parse the DCL command line for FXT_Log_Extrema.
C
C	Author: Shirley M. Read
C		STX, November 21, 1988
C
C	Invocation: Status = FXT_Parse_Command ( GMT_Start, Gmt_Stop, File_Seg,
C               Filter, Init_Xtrm, Plot, Min_Offset, Max_Offset, Flags, Report,
C		Report_File, CurGMT )
C
CH	Change Log:
CH        SPR 6727, FXT must support talaris 1590t printer
CH              H. Wang, STX, May 15, 1990
CH	  SPR 4171, Standardize report file name,
CH		Larry P. Rosen, STX, 29 August 1990
C	  ----------------------------------------------------------------------
C

C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Gmt_Start     C*14            Selected start time
C	  Gmt_Stop      C*14            Selected stop time
C	  File_Seg      C*39            Selected data segment (file)
C	  Filter        I*4             Filter value for # consecutive extrema
C	  Init_Xtrm     L*1             Flag for initialization of extrema 
C	 				record to start over with extrema checks
C	  Plot          I*4             Plot type for Engplots command file:
C				        0=none, 1=lineprinter, 2=laserprinter
C	  Min_Offset    I*4             Number of hours to offset the first 
C	                                extrema exceeded for start plot time 
C	  Max_Offset    I*4             Number of hours to offset the first 
C	                                extrema exceeded for stop plot time 
C	  Flags         L*1             Enable/disable flags to be used or not
C	  Report        L*1             Flag for printing a report file or not
C	  Report_File	C*33            File name for report
C	  CurGMT	C*14            Current run time GMT
C	
C	Subroutines Called:
C
C	  UPM_Present
C	  UPM_Get_Value
C	  UPM_Get_Longword
C	  Time_LT
C	  CT_Gmt_to_Binary
C	  Lib$Signal
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C	  FXT_Msg -- External message params needed for FXT
C	  UPM_Stat_Msg -- UPM facility status and message file
C	  $SSdef     -- System status values
C
C	Processing Method:
C
C	  Set the function status to normal.
C	  If Engfile is on the command line,
C	    Get the name of the FDQ_Eng file.
C	  Else if Jstart is on the command line,
C	    Get the start and stop times for the run.
C	  Else
C	    Signal an error.
C	    Set the function status to abort.
C	  Endif
C	  Get the value of the filter.
C	  If the Init qualifier is on the command line, set the Xtrm_Init flag.
C         If the Plots qualifier is on the command line, get the Plot type.
C	  If the Plot flag is greater than zero,
C	    Get the Minimum offset time.
C	    Get the Maximum offset time.
C	  Endif
C         If the Flags qualifier is on the command line, set the Flags.
C         Get the current run time.
C         Convert the run time into GMT.
C         Create the default report file name.
C	  If the Report qualifier is on the command line,
C            Set the report qualifier to true,
C            Get the report file name if entered,
C               If none entered, use the default report file name.
C         Else if a report is defaulted in the CLD file:
C            Set the report qualifier to true,
C            Use the default report file name.
C         Else set the report flag to false
C	  Endif
C	  Return with function status.
C
C------------------------------------------------------------------------------
C	
	Implicit None

!	Passed Parameters.

	Character*14 Gmt_Start		! Selected start time
	Character*14 Gmt_Stop		! Selected stop time
	Character*39 File_Seg           ! File segment name
	Integer*4    Filter             ! Filter value for checking extrema
	Logical*1    Init_Xtrm          ! Flag to initialize extrema record
	Integer*4    Plot               ! Flag for type of plot if an Engplots
                                        ! command file is requested
        Integer*4    Min_Offtime        ! Hours to subtract from first extrema
                                        ! exceeded time for plot start
        Integer*4    Max_Offtime        ! Hours to add to first extrema
                                        ! exceeded time for plot stop
	Logical*1    Flags              ! Flag to use enable/disable flags
	Logical*1    Report             ! Flag to print a report
	Character*33 Report_File	! File name for report
	Character*14 CurGMT		! Current run time GMT

!	Include Files.

	Include '(FXT_Msg)'
	Include '($SSDef)'
	Include '(Upm_Stat_Msg)'

!	Functions.

	Integer*4 UPM_Present
	Integer*4 UPM_Get_Longword
	Integer*4 UPM_Get_Value
	Logical*1 Time_LT
	Integer*4 Sys$Gettim

!	Local Declarations.

	Integer*4 Status	! Return status
	Integer*4 Ret_Status	! Return status
	Integer*2 Txtlen        ! Length of text string
	Character*14 Txtval     ! String to hold text
	Character*39 Blank      ! Blanks
	Character*1  Blanks(39) / 39 * ' ' /
	Equivalence  ( Blank, Blanks(1) )
	Integer*4 Start(2)      ! Binary start time
	Integer*4 Stop(2)       ! Binary stop time
	Character*33 Report_Default	! Default Report filename
	Integer*4    Curtime(2)         ! System binary time
	
!	Set function return to normal.

	FXT_Parse_Command = %loc(FXT_Normal)

!	Get the FDQ_ENG file segment or the time range.

	File_Seg(1:39) = Blank(1:39)

	If (UPM_Present('ENGFILE')) Then
	  Status = UPM_Get_Value('ENGFILE', File_Seg, Txtlen)
	  If (Status .Ne. SS$_Normal) Then
	    FXT_Parse_Command = %loc(FXT_Aberr)
	    Call Lib$Signal (%val(Status))
	  Endif
	Elseif (UPM_Present('JSTART')) Then
	  Gmt_Start = '85001000000000'
	  Status = UPM_Get_Value('JSTART', Txtval, Txtlen)
	  If (Status .Eq. SS$_Normal) Then
	    Gmt_Start(1:Txtlen) = Txtval(1:Txtlen)
	  Else
            FXT_Parse_Command = %loc(FXT_Aberr)
	    Call Lib$Signal (%val(Status))
	  Endif
	  Gmt_Stop = '99365235959900'
	  Status = UPM_Get_Value('JSTOP', Txtval, Txtlen)
	  If (Status .Eq. SS$_Normal) Then
	    Gmt_Stop(1:Txtlen) = Txtval(1:Txtlen)
	  Else
            FXT_Parse_Command = %loc(FXT_Aberr)
	    Call Lib$Signal (%val(Status))
	  Endif
	  Call CT_Gmt_to_Binary( Gmt_Start, Start)
	  Call CT_Gmt_to_Binary( Gmt_Stop, Stop)
	  If (Time_LT(Stop, Start)) Then
            FXT_Parse_Command = %loc(FXT_Aberr)
	    Call Lib$Signal (FXT_InvTim)
	  Endif	    
	Else
          FXT_Parse_Command = %loc(FXT_Aberr)
	  Call Lib$Signal (FXT_MissEngSel)
	Endif

!	Get the filter value for checking extrema.

	If ((UPM_Present('FILTER') .Eq. UPM_Pres) .Or. 
	1   (UPM_Present('FILTER') .Eq. UPM_Defaulted) .Or.
	2   (UPM_Present('FILTER') .Eq. SS$_Normal)) Then
	  Status = UPM_Get_Longword ('FILTER', Filter)
	  If (Status .Ne. SS$_Normal) Then
	    FXT_Parse_Command = %loc(FXT_Aberr)
	    Call Lib$Signal (%val(Status))
	  Endif
	Endif

!	Get the initiialization qualifier.

	If (UPM_Present('INIT')) Then
	  Init_Xtrm = .True.
	Else
	  Init_Xtrm = .False.
	Endif
C

!	Get the plots qualifier.

	Plot = 0		! Default = no plots

	If (UPM_Present('PLOTS')) Then
	  Status = UPM_Get_Value ('PLOTS', Txtval, Txtlen)
	  If (Status .Ne. SS$_Normal) Then
	    FXT_Parse_Command = %loc(FXT_Aberr)
	    Call Lib$Signal (%val(Status))
	  Else
	    If ( Txtval(1:3) .Eq. 'NON' ) Then
	      Plot = 0
	    Elseif ( Txtval(1:3) .Eq. 'LIN' ) Then
	      Plot = 1
	    Elseif ( Txtval(1:6) .Eq. 'LASQMS' ) Then
	      Plot = 2
	    Elseif ( Txtval(1:5) .Eq. 'LASTF' ) Then
	      Plot = 3
	    Else
	      Plot = - 1
	    Endif	 
	    If ((Plot .Lt. 0) .Or. (Plot .Gt. 3)) Then
	      FXT_Parse_Command = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_InvPlot)
	    Endif
	  Endif
	Endif

!	Get the minimum offset time for the plot start and maximum offset
!	time for the plot stop.

	If ((Plot .Gt. 0) .And. 
	1   (FXT_Parse_Command .Eq. %loc(FXT_Normal))) Then
	  If ((UPM_Present('MIN_OFFTIME') .Eq. UPM_Pres) .Or. 
	1   (UPM_Present('MIN_OFFTIME') .Eq. UPM_Defaulted) .Or.
	2   (UPM_Present('MIN_OFFTIME') .Eq. SS$_Normal)) Then
	    Status = UPM_Get_Longword ('MIN_OFFTIME', Min_Offtime)
	    If (Status .Ne. SS$_Normal) Then
	      FXT_Parse_Command = %loc(FXT_Aberr)
	      Call Lib$Signal (%val(Status))
	    Endif
	  Endif
	  If ((UPM_Present('MAX_OFFTIME') .Eq. UPM_Pres) .Or. 
	1   (UPM_Present('MAX_OFFTIME') .Eq. UPM_Defaulted) .Or.
	2   (UPM_Present('MAX_OFFTIME') .Eq. SS$_Normal)) Then

	    Status = UPM_Get_Longword ('MAX_OFFTIME', Max_Offtime)
	    If (Status .Ne. SS$_Normal) Then
	      FXT_Parse_Command = %loc(FXT_Aberr)
	      Call Lib$Signal (%val(Status))
	    Endif
	  Endif
	Endif

!	Get the flags qualifier.

	If (UPM_Present('FLAGS')) Then
	  Flags = .True.
	Else
	  Flags = .False.
	Endif

! Check status before processing and get the system time.

	  Status = Sys$Gettim( Curtime )
	  If ( Status .Eq. SS$_Normal ) Then
	    Call CT_Binary_To_Gmt( Curtime, Curgmt )
	  Else
	    FXT_Parse_Command = %loc(Fxt_Aberr)
	    Call Lib$Signal(FXT_GettimErr, %val(1), %val(Status))
	  Endif
	  Report_Default = 'FXT_' // GMT_Start(1:7) // '_' // 
	1     GMT_Stop(1:7) // '.REP_' // Curgmt(1:9)

!	Get the report qualifier.

	   status = upm_present('REPORT')
	   IF (status .EQ. upm_pres) THEN
              Report = .True.
              ret_status = upm_get_value('REPORT',Report_File,Txtlen)
	      IF (ret_status .EQ. upm_absent) THEN
		 Report_File = Report_Default
	      END IF
	   ELSE IF (status .EQ. upm_defaulted) THEN
	      Report = .True.
	      Report_File = Report_Default
	   ELSE
              Report = .False.
	   END IF

	Return
	End
