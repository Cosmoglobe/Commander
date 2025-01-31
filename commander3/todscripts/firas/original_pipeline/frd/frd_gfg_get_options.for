C--------------------------------------------------------------------------
	Integer*4 Function FRD_GFG_Get_Options ( Max_Sec, Report, Report_File,
	1				         CMD_Line, GMT_Start, GMT_Stop)
C--------------------------------------------------------------------------
C       Purpose: To parse the DCL command line for FRD_GFG.
C       
C       Programmer: Larry Rosen & Nilo Gonzales, STX, 14 March 1991
C
C       Invocation: Status = FRD_GFG_Get_Options (Max_Sec, Report, Report_File,
C				CMD_Line, GMT_Start, GMT_Stop)
C
CH      Change Log:
CH	  20 Sept. 1991, Larry P. Rosen, STX
CH	  Added JSTART, JSTOP qualifiers to command line.
CH
C--------------------------------------------------------------------------
C       Output Parameters:
C         Name             Type        Description
C         -----------------------------------------------------------------
C         Max_Sec          I*4         Maximum number of seconds between
C                                      major frames before a gap is declared.
C                                      The default time is 40 seconds.
C         Report           L*1         Flag for printing a report file or not
C         Report_File      C*21        Report file name
C	  CMD_Line(2)	   C*79	       Command lines string
C	  Gmt_Start	   C*14        Data Start time for running GFG
C	  GMT_Stop	   C*14        Data Stop time for GFG
C                                             
C      Subroutine Called:
C         
C         UPM_Present
C	  UPM_Get_Value
C 	  UPM_Get_Longword
C	  CT_Gmt_to_Binary
C	  sys$gettim
C
C       Externals
C
C         FRD_Normal
C         FRD_ABerr
C         FRD_GettimErr
C
C       Include Files:
C       $SSdef              -- System status values
C	upm_stat_msg	    -- Return values for UPM function
C	FUT_Params
C
C      Processing Method: PDL for FRD_GFG_Get_Options
C   
C      Get the Max_Sec qualifier, the maximum number of seconds allowed
C        between major frames before a gap in the telemetry is established,
C        from the command line. The qualifier is defaulted.
C
C	If Qualifier Report is present, then
C   	   Set the report flag to true.
C	   If report name entered
C	      Use entered name
C	   Else
C	      Use default name
C	   Endif
C	Elseif Qualifier Report is negated or not present, then
C 	   Set the report flag to false.
C	Endif
C
C	Return with normal or error status. 
C      
C--------------------------------------------------------------------------
	Implicit None

C   	Passed Parameters:
	Integer*4	Max_Sec		! Maximum number of seconds between
                               		! major frames before a gap is declared.
	Logical*1	Report		! Flag for printing a report file or not
	Character*21	Report_File	! Report file name
	Character*79	CMD_Line(2)	! Command line with defaults
	Character*14	GMT_Start	! Data start time
	Character*14	GMT_Stop	! Data stop time

C      Include Files:
	
	Include '(upm_stat_msg)'
	Include '($ssdef)'
	Include '(FUT_Params)'

C      Functions:

	Integer*4 UPM_Present
	Integer*4 UPM_Get_Value
  	Integer*4 UPM_Get_Longword
	Integer*4 CT_Gmt_to_Binary
	Integer*4 sys$gettim

	External  FRD_Normal
	External  FRD_ABerr
	External  FRD_GettimErr

C      Local Variables:

	Integer*4	Status		! Return status
	Integer*4	Ret_Status	! Return status
	Character*14	CurGMT		! Current run GMT
	Character*21	Report_Default	! Default Report filename
	Integer*2	Txtlen		! String to hold text
	Character*14	Txtval		! String read from command line
	Integer*4	Curtime(2)	! System binary time
	Character*6	Mstring		! Max_sec as a string
	Integer*2	lnum		! length of integer as string
	Integer*2	Nstart, Nend	! start and end of string of integer
 
C     Set function return to normal.

	FRD_GFG_Get_Options = %loc(FRD_Normal)

	CMD_Line(1)(1:13) = 'FGFG/MAX_SEC='

C     Get the Max_Sec qualifier, the number of seconds allowed for the major 
C	frame before a gap is declared.

	If ((UPM_Present('MAX_SEC') .Eq. UPM_Pres) .Or. 
	1   (UPM_Present('MAX_SEC') .Eq. UPM_Defaulted) .Or.
	2   (UPM_Present('MAX_SEC') .Eq. SS$_Normal)) Then
	    Status = UPM_Get_Longword ('MAX_SEC', Max_Sec)
	    If (Status .Ne. SS$_Normal) Then
	      FRD_GFG_Get_Options = %loc(FRD_Aberr)
	      Call Lib$Signal (%val(Status))
	    Else
	       lnum = int(log10(real(max_sec))) + 1
	       Write (Mstring,10) max_sec
  10	       Format (I6)
	       nstart = 7 - lnum
	       Nend = 13 + lnum
	       CMD_Line(1)(14:Nend) = Mstring (nstart:6)
	    EndIf
	Endif
	IF (upm_present('jstart')) THEN
	   gmt_start = fac_jstart_default
	   status = upm_get_value('jstart', txtval, txtlen)
	   IF (status .EQ. ss$_normal) THEN
	      Nstart = Nend + 1
	      Nend = Nstart + 21
	      gmt_start(1:txtlen) = txtval(1:txtlen)
	      CMD_Line(1)(Nstart:Nend) = '/JSTART=' // gmt_start
	   ELSE
	      FRD_GFG_Get_Options = %loc(FRD_Aberr)
	      Call Lib$Signal (%val(Status))
	   ENDIF
	ELSE
	   Write (6,*) ' Must use JSTART on command line.'
	   FRD_GFG_Get_Options = %loc(FRD_Aberr)
	   Call Lib$Signal (%val(Status))
	ENDIF
	IF (upm_present('jstop')) THEN
	   gmt_stop = fac_jstop_default
	   status = upm_get_value('jstop', txtval, txtlen)
	   IF (status .EQ. ss$_normal) THEN
	      Nstart = Nend + 1
	      Nend = Nstart + 20
	      gmt_stop(1:txtlen) = txtval(1:txtlen)
	      CMD_Line(1)(Nstart:Nend) = '/JSTOP=' // gmt_stop
	   ELSE
	      FRD_GFG_Get_Options = %loc(FRD_Aberr)
	      Call Lib$Signal (%val(Status))
	   ENDIF
	ELSE
	   Write (6,*) ' Must use JSTOP on command line.'
	   FRD_GFG_Get_Options = %loc(FRD_Aberr)
	   Call Lib$Signal (%val(Status))
	ENDIF

C      Check status before processing and get the system time.

	Status = Sys$Gettim( Curtime )
	If ( Status .Eq. SS$_Normal ) Then
	   Call CT_Binary_To_Gmt( Curtime, Curgmt )
	Else
	   FRD_GFG_Get_Options = %loc(FRD_Aberr)
	   Call Lib$Signal(FRD_GettimErr, %val(1), %val(Status))
	Endif
	Report_Default = 'FRD_GFG.REP_' // Curgmt(1:9)

C	Get the report qualifier.

	Nstart = 15			! indent continue line of command line
	status = upm_present('report')
	IF (status .EQ. upm_pres) THEN
           Report = .True.
           ret_status = upm_get_value('report',Report_File,Txtlen)
	   IF (ret_status .EQ. upm_absent) THEN
	      Report_File = Report_Default
	      Txtlen = 21
	   END IF
	   Nend = Nstart + txtlen + 7
	   CMD_Line(2)(Nstart:Nend) = '/REPORT=' // Report_File
	ELSE IF (status .EQ. upm_defaulted) THEN
	   Report = .True.
	   Report_File = Report_Default
	   Txtlen = 21
	   Nend = Nstart + txtlen + 7
	   CMD_Line(2)(Nstart:Nend) = '/REPORT=' // Report_File
	ELSE
           Report = .False.
	   Nend = Nstart + 8
	   CMD_Line(2)(Nstart:Nend) = '/NOREPORT'
	END IF

	Return
	End
