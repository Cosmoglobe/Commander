	Integer*4  Function  FAD_Parse_Command  ( chan, scan, filext, modelext,
	1                                         report, report_name,
	2                                         cmdline, current_adt)

c------------------------------------------------------------------------------
c   Purpose: Parse the DCL command line for qualifier values for FAD.
c
c   Input Parameters: none
c   Output Parameters:
c      character*2   chan          -- channel (RH, RL, LH, or LL)
c      character*2   scan          -- scan mode (SS, SF, LS, or LF)
c      character*22  filext        -- file extension of FCF_SKY file
c      character*22  modelext      -- FISH model name - FEX_EJ file extension
c      logical*1     report        -- true if report is to be written (default)
c      character*45  report_name   -- name of the report to write
c      character*79  cmdline(2)    -- command line as a character string
c      integer*4     current_adt(2)  -- current system time in ADT format
c
c   Subroutines and functions Called:
c
c      UPM_Present
c      UPM_Get_Value
c      Str$Upcase
c      Sys$GetTim
c      CT_Binary_To_GMT
c
c   Include Files:
c      FAD_msg       -- External message params needed for FAD
c      UPM_stat_msg  -- UPM facility status and message file
c      $SSdef        -- System status values
c
c   Author:  Larry Paris Rosen, Hughes STX, 19 April 1993
c------------------------------------------------------------------------------
	Implicit None

c Include files

	Include		'($ssdef)'
	Include		'(fad_msg)'
	Include		'(upm_stat_msg)'

c Functions

	Integer*4	Upm_Present
	Integer*4	Upm_Get_Longword
	Integer*4	Upm_Get_Value
	Logical*1	Time_LT

c Passed parameters

	Character*2	chan, scan
	Character*22	filext
	Character*22	modelext
	Logical*1	report
	Character*45	report_name
	Character*79	cmdline(2)
	Integer*4	current_adt(2)		! current time in adt

c Local

	Character*79	blank_line		! blank line
	Character*1	lin(79) / 79 * ' ' /	! blanks
	Equivalence	(blank_line, lin(1))
	Integer*4	rstat			! return status
	Integer*4	txtlen
	Integer*2	pos			! position in string
	Integer*2	cend			! end of cmd string
	Integer*2	cstart			! start of cmd string
	Logical*1	repdef		! flag to use default report name
	Character*14	current_time		! current time in gmt

c------------------------------------------------------------------------------

c Begin

c  Set function return to normal and initialize command line.

	FAD_Parse_Command = %Loc (FAD_Normal)
	cmdline (1) = blank_line
	cmdline (2) = blank_line
	cmdline (1)(1:3) = 'FAD'

c Get channel

	rstat = Upm_Present ('chan')
	If ((rstat .EQ. upm_pres) .OR. (rstat .EQ. upm_defaulted)) Then
	   rstat = Upm_Get_Value ('chan', chan, txtlen)
	   If (rstat .EQ. SS$_Normal) Then
	      Call Str$Upcase (chan, chan)
	      cmdline (1)(4:11) = '/CHAN=' // chan
	   Else
	      Call lib$signal(fad_parserr)
	      FAD_Parse_Command = %Loc (FAD_Abort)
	   Endif
	Else
      	      Call lib$signal(fad_parserr)
	   FAD_Parse_Command = %Loc (FAD_Abort)
	Endif

c Get scan mode

	If (FAD_Parse_Command .EQ. %Loc (FAD_Normal)) Then
	   rstat = Upm_Present ('scan')
	   If ((rstat .EQ. upm_pres) .OR. (rstat .EQ. upm_defaulted)) Then
	      rstat = Upm_Get_Value ('scan', scan, txtlen)
	      If (rstat .EQ. SS$_Normal) Then
	         Call Str$Upcase (scan, scan)
	         cmdline (1)(12:19) = '/SCAN=' // scan
	      Else
       	         Call lib$signal(fad_parserr)
	         FAD_Parse_Command = %Loc (FAD_Abort)
	      Endif
	   Else
	      Call lib$signal(fad_parserr)
	      FAD_Parse_Command = %Loc (FAD_Abort)
	   Endif
	Endif

c Get file extension of FCF_SKY

	If (FAD_Parse_Command .EQ. %Loc (FAD_Normal)) Then
	   rstat = Upm_Present ('filext')
	   If (rstat .EQ. upm_pres) Then
	      rstat = Upm_Get_Value ('filext', filext, txtlen)
	      If (rstat .EQ. SS$_Normal) Then
	         Call Str$Upcase (filext, filext)
	         pos = Index (filext,' ')
	         If (pos .NE. 0) Then
	            txtlen = pos
	         Endif
	         cend = 20 + 8 + txtlen - 1
	         cmdline (1)(20:cend) = '/FILEXT=' // filext (1:txtlen)
	      Else
	         Call lib$signal(fad_parserr)
	         FAD_Parse_Command = %Loc (FAD_Abort)
	      Endif
	   Else
	      Call lib$signal(fad_parserr)
	      FAD_Parse_Command = %Loc (FAD_Abort)
	   Endif
	Endif

c Get the file extension of FEX_EJ file

	If (FAD_Parse_Command .EQ. %Loc (FAD_Normal)) Then
	   rstat = Upm_Present ('model')
	   If (rstat .EQ. upm_pres) Then
	      rstat = Upm_Get_Value ('model', modelext, txtlen)
	      If (rstat .EQ. SS$_Normal) Then
	         Call Str$Upcase (modelext, modelext)
	         pos = Index (modelext,' ')
	         If (pos .NE. 0) Then
	            txtlen = pos
	         Endif
	         cstart = cend + 1
	         cend = cend + 7 + txtlen
	         cmdline (1)(cstart:cend) = '/MODEL=' // modelext (1:txtlen)
	      Else
	         Call lib$signal(fad_parserr)
	         FAD_Parse_Command = %Loc (FAD_Abort)
	      Endif
	   Else
	      Call lib$signal(fad_parserr)
	      FAD_Parse_Command = %Loc (FAD_Abort)
	   Endif
	Endif

c Get report name

	If (FAD_Parse_Command .EQ. %Loc (FAD_Normal)) Then

	   Call Sys$GetTim (current_adt)      ! get current system time

	   rstat = Upm_Present ('report')
	   If (rstat .EQ. upm_pres) Then
	      report = .TRUE.
	      rstat = Upm_Get_Value ('report', report_name, txtlen)
	      If (rstat .EQ. upm_absent) Then
	         repdef = .TRUE.
	      Else
	         repdef = .FALSE.
	         Call Str$Upcase (report_name, report_name)
	      Endif
	   Elseif (rstat .EQ. Upm_Defaulted) Then
	      report = .TRUE.
	      repdef = .TRUE.
	   Else
	      report = .FALSE.
	      cstart = 4
	      cend = cstart + 9 - 1
	      cmdline (2)(cstart:cend) = '/NOREPORT'
	   Endif
	   If (repdef) Then

c Make default report name from FCF_SKY file extension.
c Assumes file extentsion is like 'ED_8912345_9012345'    or
c                            like 'ED_8912345_9012345;78'

	      Call CT_Binary_To_GMT (current_adt, current_time)

	      report_name = 'FAD_' // filext (4:18) // '_' // chan
	1        // scan // '.REP_' // current_time

	   Endif
	   If (report) Then
	      pos = Index (report_name,' ')
	      If (pos .NE. 0) Then
	         txtlen = pos
	      Else
	         txtlen = len (report_name)
	      Endif
	      cstart = 4
	      cend = cstart + 8 + txtlen - 1
	      cmdline (2)(cstart:cend) = '/REPORT=' // report_name (1:txtlen)
	   Endif
	Endif
	Return
	End
