C-------------------------------------------------------------------------------
	INTEGER*4 FUNCTION FRD_FACR_Get_Options ( gmt_start, gmt_stop, Avetime,
	1	Report, Report_File, Cmd_Line )
C-------------------------------------------------------------------------------
C       Purpose: To parse the DCL command line for FRD_FACR.
C       
C       Programmer: Nilo G. Gonzales/STX, April 5, 1991
C
C       Invocation: Status = FRD_FACR_Get_Options (gmt_start, gmt_stop, 
C                                      Avetime, Report, Report_File, Cmd_Line)
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Gmt_Start	C*14		Selected start time
C	  Gmt_Stop	C*14		Selected stop time
C         Avetime       R*4             Averaging time 
C	  Report_file	C*34		Report file name
C	  Report	L*1		Flag for printing a report file or not
C   	  Cmd_Line      C*79            Command line
C
C	Subroutines Called:
C
C	  UPM_Present
C	  UPM_Get_Value
C	  UPM_Get_Longword
C	  UPM_Get_Float
C	  Time_LT
C	  CT_Gmt_to_Binary
C	  Lib$Signal
C
C	Include Files:
C	  UPM_Stat_Msg -- UPM facility status and message file
C	  $SSdef       -- System status values
C	  FUT_Params   -- FIRAS global parameters
C
C	Processing Method: PDL for FRD_FACR_Get_Options
C  Begin
C	If Qualifier Jstart is present, then
C	   Get the start time and stop time from the command line.
C	Endif
C
C	If Qualifier Avetime is present, then
C	   Get the Avetime value from the command line.
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
C  End
C------------------------------------------------------------------------------
	IMPLICIT NONE

C Passed Parameters.

	CHARACTER*14	gmt_start	! selected start time
	CHARACTER*14	gmt_stop	! selected stop time
	REAL*4          AVETIME         ! Averaging time 
	LOGICAL*1	report		! flag to print a report
	CHARACTER*34	report_file	! file name for report
	CHARACTER*79	cmd_line(2)	! command line including defaults
	
C  Include Files.

	INCLUDE '(fut_params)'
	INCLUDE '($ssdef)'
	INCLUDE '(upm_stat_msg)'

C  Functions.

	INTEGER*4	upm_present
	INTEGER*4	upm_get_longword
	INTEGER*4	upm_get_value
	INTEGER*4	upm_get_float
	LOGICAL*1	time_LT
	INTEGER*4	sys$gettim
	EXTERNAL	frd_normal
	EXTERNAL	frd_aberr
	EXTERNAL	frd_gettimerr

C  Local Declarations.

	INTEGER*4	status			! return status
	INTEGER*4	ret_status		! return status
	INTEGER*2	txtlen			! length of text string
	CHARACTER*14	txtval			! string to hold text
	INTEGER*4	tstart(2)		! binary start time
	INTEGER*4	tstop(2)		! binary stop time
	CHARACTER*34	report_default		! default report filename
	INTEGER*4	curtime(2)		! system binary time
	CHARACTER*14	curgmt			! system time in gmt
	INTEGER*2	clen			! length of current command line
	INTEGER*2	lnum			! length of number as string
	CHARACTER*8	numstr			! number as string
	CHARACTER*42	repstr			! string report qualif and name
	
C  Set function return to normal and initialize command line.

	FRD_FACR_Get_Options = %loc(FRD_Normal)
	cmd_line(1)(1:4) = 'FACR'
	IF (upm_present('jstart')) THEN
	    gmt_start = fac_jstart_default
	    status = upm_get_value('jstart', txtval, txtlen)
	    IF (status .EQ. ss$_normal) THEN
	       gmt_start(1:txtlen) = txtval(1:txtlen)
	       cmd_line(1)(5:26) = '/JSTART=' // gmt_start
	    ELSE
               FRD_FACR_Get_Options = %loc(FRD_Aberr)
	       CALL lib$signal (%val(status))
	    ENDIF
	    gmt_stop = fac_jstop_default
	    status = upm_get_value('jstop', txtval, txtlen)
	    IF (status .EQ. ss$_normal) THEN
	       gmt_stop(1:txtlen) = txtval(1:txtlen)
	       cmd_line(1)(27:47) = '/JSTOP=' // gmt_stop
	    ELSE
               FRD_FACR_Get_Options = %loc(FRD_Aberr)
	       CALL lib$signal (%val(status))
	    ENDIF
	    CALL ct_gmt_to_binary( gmt_start, tstart)
	    CALL ct_gmt_to_binary( gmt_stop, tstop)
	    IF (time_lt(tstop, tstart)) THEN
               FRD_FACR_Get_Options = %loc(FRD_Aberr)
	       CALL lib$signal (FRD_Aberr)
	    ENDIF
	Else
           FRD_FACR_Get_Options = %loc(FRD_Aberr)
	   Call Lib$signal(FRD_Aberr)
	ENDIF

C  Get the Avetime qualifier, the number of seconds to be used for averaging
C  the calibrator resistors. The default value is 0.5 day.

	IF ((upm_present('avetime') .EQ. upm_pres) .OR.
	1   (upm_present('avetime') .EQ. upm_defaulted) .OR.
	2   (upm_present('avetime') .EQ. ss$_normal)) THEN
	  status = upm_get_float ('avetime', avetime)
	  IF (status .NE. ss$_normal) THEN
	    FRD_FACR_Get_Options = %loc(FRD_Aberr)
	    CALL lib$signal (%val(status))
	  ELSE
	    WRITE (numstr,10) avetime
  10	    FORMAT (F8.3)
	    cmd_line(1)(48:64) = '/AVETIME=' // numstr
	  ENDIF
	ENDIF

C  Get the system time for report name.

	status = sys$gettim( curtime )
	IF ( status .EQ. ss$_normal ) THEN
	  CALL ct_binary_to_gmt( curtime, curgmt )
	ELSE
	    FRD_FACR_Get_Options = %loc(FRD_Aberr)
	    CALL lib$signal(FRD_gettimerr, %val(1), %val(status))
	ENDIF
	report_default = 'FACR_' // gmt_start(1:7) // '_' // 
     .                    gmt_stop(1:7) // '.REP_' // curgmt(1:9)

C  Get the report qualifier.

	status = upm_present('report')
	IF (status .EQ. upm_pres) THEN
          report = .TRUE.
          ret_status = upm_get_value('report',report_file,txtlen)
	  IF (ret_status .EQ. upm_absent) THEN
	    report_file = report_default
	    lnum = 42
	  ELSE
	    lnum = 8 + txtlen
	  END IF
	  repstr = '/REPORT=' // report_file
	ELSE IF (status .EQ. upm_defaulted) THEN
	  report = .TRUE.
	  report_file = report_default
	  repstr = '/REPORT=' // report_file
	  lnum = 42
	ELSE
          report = .FALSE.
	  repstr = '/NOREPORT'
	  lnum = 9
	ENDIF
	clen = 4 + lnum
	cmd_line(2)(5:clen) = repstr(1:lnum)
	RETURN
	END
