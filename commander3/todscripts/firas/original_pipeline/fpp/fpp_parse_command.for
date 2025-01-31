C-------------------------------------------------------------------------------
	INTEGER*4 FUNCTION FPP_Parse_Command ( time_present, gmt_start,
	1	gmt_stop, file_present, file_name, proc_chan, num_chan, max_sec,
	2	report, report_file, command_line, sci_gap, anc_gap,sweep_check)
C-------------------------------------------------------------------------------
C	Purpose: To parse the DCL command line for FPP_FIRAS_Pre_Processor.
C
C	Author: Shirley M. Read
C		STX, January, 1989
C
C	Invocation: Status = FPP_Parse_Command (Time_Present, Gmt_Start,
C		Gmt_Stop, File_Present, File_Name, Proc_Chan, Num_Chan, Max_Sec,
C		Report, Report_File, Command_line, sci_gap, anc_gap,sweep_check)
C
CH	Change Log:
CH
CH		Version 4.2.2 04/29/89, SPR 3415, Shirley M. Read, STX
CH			FPP tracking should not be optional. Tracking should
CH			be the default on the command line.
CH
CH		Version 5.7 2/23/90, SPR 5876, H. Wang STX
CH			IFG tracking default does not conform with
CH			rest of pipeline.
CH
CH		Version    7/24/90, SPR 4171, L. Rosen STX
CH			Get current run time and convert to GMT.
CH			Get or create default report file name.
CH
CH		Version 7.1 10/25/90, SPR 7294, H. Wang STX
CH			FPP needs to check for gaps of FPP_SDF records.
CH
CH		New Version  2/27/91, Larry P. Rosen, STX
CH			New FPP requirements and design. Fewer qualifiers here.
CH			Construct command_line with defaults.
CH
CH		New Version.1 17 October 1991, Larry P. Rosen, Hughes STX
CH		SPR 9167.  Add SCREEN_SWEEP qualifier. Pass logical SWEEP_CHECK.
CH		True = flag ifg's with sweeps not equal to 0, 1, 4, or 16.
C	  ----------------------------------------------------------------------
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Time_Present	I*4		/JSTART, /JSTOP presence on command line
C	  Gmt_Start	C*14		Selected start time
C	  Gmt_Stop	C*14		Selected stop time
C	  File_Present	I*4		/FILENAME presence on command line
C	  File_Name	C*39		Selected data segment (file)
C	  Proc_Chan(4)	I*2		Channels to be processed.
C	  Num_Chan	I*2		Number of channels to be processed
C	  Max_Sec	I*4		Maximum number of seconds between major
C					frames before a gap is declared
C	  Report	L*1		Flag for printing a report file or not
C	  Report_file	C*33		Report file name
C	  Command_line(3) C*79		Command line including defaults
C	  Sci_Gap	I*4		Max SDF data gap else report
C	  Anc_gap	I*4		Max ANC data gap else report
C	  Sweep_check	L*1		Flag ifg's with sweeps not = 0,1,4,16
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
C	Include Files:
C	  FPP_Msg      -- External message params needed for FPP
C	  UPM_Stat_Msg -- UPM facility status and message file
C	  $SSdef       -- System status values
C	  FUT_Params   -- FIRAS global parameters
C
C	Processing Method: PDL for FPP_Parse_Command
C
C	(construct command_line with defaults during routine)
C
C	If Qualifier Filename is present, then
C	   Get the complete NFS_SDF filename from the command line.
C	Endif
C	If Qualifier Jstart is present, then
C	   Get the start time and stop time from the command line.
C	Endif
C
C	Get the Channel from the command line. The default is 'ALL'.
C
C	Get the Max_Sec qualifier, the maximum number of seconds allowed
C          between major frames before a gap in the telemetry is established,
C          from the command line. This qualifier is defaulted at 40.
C
C	If Qualifier Sci_Gap is present, Then
C	   If a value is entered Then
C	      Set sci_gap to that value.
C	   Else
C	      Set sci_gap to the default, 300.
C	   Endif
C	Else
C	   Set sci_gap to the default, 300.
C	Endif
C
C	If Qualifier Anc_Gap is present, Then
C	   If a value is entered Then
C	      Set anc_gap to that value.
C	   Else
C	      Set anc_gap to the default, 90.
C	   Endif
C	Else
C	   Set anc_gap to the default, 90.
C	Endif
C
C	If Qualifier Screen_Sweeps is present, Then
C	   Set Sweep_Check to True.
C	Else If Qualifier is defaulted, Then
C	   Set Sweep_Check to True.
C	Else
C	   Set Sweep_Check to False.
C	Endif
C
C	Get system time and generate default report file name.
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
C------------------------------------------------------------------------------
	IMPLICIT NONE

C Passed Parameters.

	INTEGER*4	file_present	! /filename selected
	INTEGER*4	time_present	! /jstart, /jstop selected
	CHARACTER*14	gmt_start	! selected start time
	CHARACTER*14	gmt_stop	! selected stop time
	CHARACTER*39	file_name	! file segment name
	INTEGER*2	proc_chan(4)	! channels to be processed
	INTEGER*2	num_chan	! number of channels to process
	CHARACTER*4	channel		! user specified channel
	INTEGER*2	chan		! channel counter
	INTEGER*4	max_sec		! maximum number of seconds between
					! major frames before a gap is declared
	LOGICAL*1	report		! flag to print a report
	CHARACTER*33	report_file	! file name for report
	CHARACTER*79	command_line(3)	! command line including defaults
	INTEGER*4	sci_gap		! max sdf data gap else report
	INTEGER*4	anc_gap		! max anc data gap else report
	LOGICAL*1	sweep_check	! Whether to perform check of sweeps.

C  Include Files.

	INCLUDE '(fut_params)'
	INCLUDE '(fpp_msg)'
	INCLUDE '($ssdef)'
	INCLUDE '(upm_stat_msg)'

C  Functions.

	INTEGER*4	upm_present
	INTEGER*4	upm_get_longword
	INTEGER*4	upm_get_value
	LOGICAL*1	time_LT
	INTEGER*4	sys$gettim

C  Local Declarations.

	INTEGER*4	status			! return status
	INTEGER*4	ret_status		! return status
	INTEGER*2	txtlen			! length of text string
	CHARACTER*14	txtval			! string to hold text
	CHARACTER*39	blank			! blanks
	CHARACTER*1	blanks(39) / 39 * ' ' /
	EQUIVALENCE	( blank, blanks(1) )
	INTEGER*4	tstart(2)		! binary start time
	INTEGER*4	tstop(2)		! binary stop time
	CHARACTER*33	report_default		! default report filename
	INTEGER*4	curtime(2)		! system binary time
	CHARACTER*14	curgmt			! system time in gmt
	CHARACTER*79	line			! blank line
	CHARACTER*1	lin(79) / 79 * ' ' /	! blanks
	EQUIVALENCE	( line, lin(1) )
	INTEGER*2	clen			! length of current command line
	INTEGER*2	cnex			! next position in command line
	INTEGER*2	clin /1/		! which line of command line
	INTEGER*2	lnum			! length of number as string
	INTEGER*2	nstart			! position of 1st digit string
	CHARACTER*10	numstr			! number as string
	CHARACTER*41	repstr			! string report qualif and name
	
C  Set function return to normal and initialize command line.

	fpp_parse_command = %loc(fpp_normal)
	command_line(1)(1:79) = line
	command_line(2)(1:79) = line
	command_line(3)(1:79) = line
	command_line(1)(1:3) = 'FPP'
	clen = 3
	cnex = clen + 1

C  Get the science file segment or the time range.

	file_name(1:39) = blank(1:39)
	file_present = fac_not_present
	IF (upm_present('filename')) THEN
	  file_present = fac_present
	  status = upm_get_value('filename', file_name, txtlen)
	  IF (status .NE. ss$_normal .AND. status .NE. upm_absent) THEN
	    fpp_parse_command = %loc(fpp_aberr)
	    CALL lib$signal (%val(status))
	  ELSE
	    CALL str$upcase (file_name,file_name)
	    clen = clen + 10 + txtlen
	    command_line(1)(cnex:clen) = '/FILENAME=' // file_name
	    cnex = clen + 1
	  ENDIF
	ENDIF
	time_present = fac_not_present
	IF (upm_present('jstart')) THEN
	  time_present = fac_present
	  gmt_start = fac_jstart_default
	  status = upm_get_value('jstart', txtval, txtlen)
	  IF (status .EQ. ss$_normal) THEN
	    gmt_start(1:txtlen) = txtval(1:txtlen)
	    clen = clen + 22
	    command_line(1)(cnex:clen) = '/JSTART=' // gmt_start
	    cnex = clen + 1
	  ELSE
            fpp_parse_command = %loc(fpp_aberr)
	    CALL lib$signal (%val(status))
	  ENDIF
	  gmt_stop = fac_jstop_default
	  status = upm_get_value('jstop', txtval, txtlen)
	  IF (status .EQ. ss$_normal) THEN
	    gmt_stop(1:txtlen) = txtval(1:txtlen)
	    IF (clen .GT. 58) THEN
	      clin = clin + 1
	      clen = 0
	      cnex = 1
	    ENDIF
	    clen = clen + 21
	    command_line(clin)(cnex:clen) = '/JSTOP=' // gmt_stop
	    cnex = clen + 1
	  ELSE
            fpp_parse_command = %loc(fpp_aberr)
	    CALL lib$signal (%val(status))
	  ENDIF
	  CALL ct_gmt_to_binary( gmt_start, tstart)
	  CALL ct_gmt_to_binary( gmt_stop, tstop)
	  IF (time_lt(tstop, tstart)) THEN
            fpp_parse_command = %loc(fpp_aberr)
	    CALL lib$signal (fpp_invtim)
	  ENDIF	
	ENDIF
	IF (file_present.EQ.fac_not_present.AND.time_present.EQ.fac_not_present)
	1   Then
	    CALL lib$signal (fpp_MissSciSel)
            fpp_parse_command = %loc(fpp_aberr)
	EndIf

C  Get the channel selection.

	IF ((upm_present('channel') .EQ. upm_pres) .OR.
	1   (upm_present('channel') .EQ. upm_defaulted) .OR.
	2   (upm_present('channel') .EQ. ss$_normal)) THEN
	  status = upm_get_value ('channel', channel, txtlen)
	  IF (status .NE. ss$_normal) THEN
	    fpp_parse_command = %loc(fpp_aberr)
	    CALL lib$signal (%val(status))
	  ELSE
	    CALL str$upcase (channel,channel)
	    IF (clen .GT. 67) THEN
	      clin = clin + 1
	      clen = 0
	      cnex = 1
	    ENDIF
	    clen = clen + 9 + txtlen
	    command_line(clin)(cnex:clen) = '/CHANNEL=' // channel
	    cnex = clen + 1
	  ENDIF
	ENDIF

C  Parse channel spec

	IF (channel(1:3) .EQ. 'ALL') THEN
	  DO chan = 1 , 4
	    proc_chan(chan) = chan
	  ENDDO
	  num_chan = 4
	ELSE
	  num_chan = 1
	  IF ( channel(1:2) .EQ. 'RH' ) THEN
	    proc_chan(1) = 1
	  elseif ( channel(1:2).EQ. 'RL') THEN
	    proc_chan(1) = 2
	  elseif (channel(1:2) .EQ. 'LH' ) THEN
	    proc_chan(1) = 3
	  elseif (channel(1:2) .EQ. 'LL' ) THEN
	    proc_chan(1) = 4
	  ELSE
	    CALL lib$signal(fpp_bchannel,%val(1),channel)
	    fpp_parse_command = %loc(fpp_aberr)
	  ENDIF  ! channel(1:2)
	ENDIF    ! channel

C  Get the Max_Sec qualifier, the number of seconds allowed for the major
C  frame before a gap is declared.

	IF ((upm_present('max_sec') .EQ. upm_pres) .OR.
	1   (upm_present('max_sec') .EQ. upm_defaulted) .OR.
	2   (upm_present('max_sec') .EQ. ss$_normal)) THEN
	  status = upm_get_longword ('max_sec', max_sec)
	  IF (status .NE. ss$_normal) THEN
	    fpp_parse_command = %loc(fpp_aberr)
	    CALL lib$signal (%val(status))
	  ELSE
	    lnum = int(log10(real(max_sec))) + 1
	    IF (clen .GT. 70-lnum) THEN
	      clin = clin + 1
	      clen = 0
	      cnex = 1
	    ENDIF
	    clen = clen + 9 + lnum
	    WRITE (numstr,10) max_sec
  10	    FORMAT (I10)
	    nstart = 11 - lnum
	    command_line(clin)(cnex:clen) = '/MAX_SEC=' // numstr(nstart:10)
	    cnex = clen + 1
	  ENDIF
	ENDIF

C  Get SCI_GAP, the max time allowed between science records before gap rept.

	IF ((upm_present('sci_gap') .EQ. upm_pres) .OR.
	1   (upm_present('sci_gap') .EQ. upm_defaulted) .OR.
	2   (upm_present('sci_gap') .EQ. ss$_normal)) THEN
	  status = upm_get_longword ('sci_gap', sci_gap)
	  IF (status .NE. ss$_normal) THEN
	    fpp_parse_command = %loc(fpp_aberr)
	    CALL lib$signal (%val(status))
	  ELSE
	    lnum = int(log10(real(sci_gap))) + 1
	    IF (clen .GT. 70-lnum) THEN
	      clin = clin + 1
	      clen = 0
	      cnex = 1
	    ENDIF
	    clen = clen + 9 + lnum
	    WRITE (numstr,10) sci_gap
	    nstart = 11 - lnum
	    command_line(clin)(cnex:clen) = '/SCI_GAP=' // numstr(nstart:10)
	    cnex = clen + 1
	  ENDIF
	ENDIF

C  Get ANC_GAP, the max time allowed between ANC major frames before gap report.

	IF ((upm_present('anc_gap') .EQ. upm_pres) .OR.
	1   (upm_present('anc_gap') .EQ. upm_defaulted) .OR.
	2   (upm_present('anc_gap') .EQ. ss$_normal)) THEN
	    status = upm_get_longword ('anc_gap', anc_gap)
	  IF (status .NE. ss$_normal) THEN
	    fpp_parse_command = %loc(fpp_aberr)
	    CALL lib$signal (%val(status))
	  ELSE
	    lnum = int(log10(real(anc_gap))) + 1
	    IF (clen .GT. 70-lnum) THEN
	      clin = clin + 1
	      clen = 0
	      cnex = 1
	    ENDIF
	    clen = clen + 9 + lnum
	    WRITE (numstr,10) anc_gap
	    nstart = 11 - lnum
	    command_line(clin)(cnex:clen) = '/ANC_GAP=' // numstr(nstart:10)
	    cnex = clen + 1
	  ENDIF
	ENDIF

C  Get SCREEN_SWEEPS qualifier, set Sweep_check.

	IF ((upm_present('SCREEN_SWEEPS') .EQ. upm_pres) .OR.
	1   (upm_present('SCREEN_SWEEPS') .EQ. upm_defaulted) .OR.
	2   (upm_present('SCREEN_SWEEPS') .EQ. ss$_normal)) THEN
	   Sweep_check = .TRUE.
	   IF (clen .GT. 65) THEN
	      clin = clin + 1
	      clen = 0
	      cnex = 1
	   ENDIF
	   clen = clen + 14
	   command_line(clin)(cnex:clen) = '/SCREEN_SWEEPS'
	ELSE
	   Sweep_check = .FALSE.
	   IF (clen .GT. 63) THEN
	      clin = clin + 1
	      clen = 0
	      cnex = 1
	   ENDIF
	   clen = clen + 16
	   command_line(clin)(cnex:clen) = '/NOSCREEN_SWEEPS'
	ENDIF
	cnex = clen + 1

C  Get the system time for report name.

	status = sys$gettim( curtime )
	IF ( status .EQ. ss$_normal ) THEN
	  CALL ct_binary_to_gmt( curtime, curgmt )
	ELSE
	  fpp_parse_command = %loc(fpp_aberr)
	  CALL lib$signal(fpp_gettimerr, %val(1), %val(status))
	ENDIF
	report_default = 'FPP_' // gmt_start(1:7) // '_' //
	1     gmt_stop(1:7) // '.REP_' // curgmt(1:9)

C  Get the report qualifier.

	status = upm_present('report')
	IF (status .EQ. upm_pres) THEN
          report = .TRUE.
          ret_status = upm_get_value('report',report_file,txtlen)
	  IF (ret_status .EQ. upm_absent) THEN
	    report_file = report_default
	    lnum = 41
	  ELSE
	    lnum = 8 + txtlen
	  END IF
	  repstr = '/REPORT=' // report_file
	ELSE IF (status .EQ. upm_defaulted) THEN
	  report = .TRUE.
	  report_file = report_default
	  repstr = '/REPORT=' // report_file
	  lnum = 41
	ELSE
          report = .FALSE.
	  repstr = '/NOREPORT'
	  lnum = 9
	END IF
	IF (clen .GT. 79-lnum) THEN
	      clin = clin + 1
	      clen = 0
	      cnex = 1
	ENDIF
	clen = clen + lnum
	command_line(clin)(cnex:clen) = repstr(1:lnum)

	RETURN
	END
