	Integer * 4 Function  FFL_Initialize_Report (current_gmt, ncmd,
     &				 cmd_line, cmdlen, parse_status, ffli)
c-------------------------------------------------------------------------------
c
c	Function FFL_INITIALIZE_REPORT
c
c	This function initializes the FFL processing report.  The report is
c	opened, then the error handler is invoked.  Finally, the account
c	information, command line invocation, and logical pointer translations
c	are written to the report.
c
c	Author:	 Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 9 March 1993
c		 SER 10763
c
c-------------------------------------------------------------------------------
c
c	Input:
c	    current_gmt	    Ch*14	! GMT time of invocation
c	    ncmd	    I * 4	! number of command lines in invocation
c	    cmd_line(3)	    Ch*79	! command line invocation
c	    cmdlen(3)	    I * 4	! length of command lines in invocation
c	    parse_status    I * 4	! status from command line parse
c	    ffli	    Structure (declared in the FFL_Invoc include file)
c
c	Output:
c		none
c
c	Subroutines called:
c		cut_display_banner
c		cut_register_version
c		cut_translate_archive_id
c		fut_get_lun
c		LIB$Establish
c		LIB$Getjpi
c		LIB$Signal
c		STR$Trim
c
c	Include files:
c		ffl_invoc.txt
c		fut_error.txt
c		fut_params.txt
c		$jpidef
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Name of facility changed from FFI to FFL.  Fred Shuman, HSTX
c	1995 May 19.
c
c	In order to remove some of the obscurity arising from Include files that
c	   harbor hidden Common blocks, converted these Commons to Structures.
c	   This necessitates adding their Record names to the calling lists of
c	   functions that use their variables.
c	Fred Shuman, HSTX, 1995 June 14.
c
c       S. Brodd, HSTX, 12/20/95, SPR 12284.  Change banner length from 90
c                                             to 80 to improve report files.
c-------------------------------------------------------------------------------

	Implicit None

	Include '(fut_error)'
	Include '(fut_params)'		! defines FIRAS parameters, fac_*
	Include '($jpidef)'
c
c  Call arguments:
c
	Character*14	current_gmt	! GMT time of invocation
	Integer  * 4	ncmd		! number of command lines in invocation
	Character*79	cmd_line(3)	! command line invocation
	Integer  * 4	cmdlen(3)	! length of command lines in invocation
	Integer  * 4	parse_status	! command line parse return status
	Include '(ffl_invoc)'		! defines record ffli (structure invoc)
					!   and Ch*6 parameter "version"
c
c  All other variables and functions:
c
	Character*72	log_name(4)	! logical name to be translated
	Character*72	temp_log	! logical name buffer
	Character*72	trans_name(4)	! translated logical name
	Character*20	username	! user name of invoking account

	Integer  * 2	len		! string length
	Integer  * 2	temp_len	! logical name buffer length
	Integer  * 2	trans_len(4)	! translated logical name length

	Integer  * 4	io_stat		! I/O return status
	Integer  * 4	j, k		! counters
	Integer  * 4	n_names		! number of logical names to be
					!   translated
	Integer  * 4	status, rstatus, tstatus    ! return statuses

	Integer	 * 4	fut_get_lun
	Integer	 * 4	cut_display_banner
	Integer	 * 4	cut_register_version
	Integer  * 4	cut_translate_archive_id

	External	ffl_invaltime
	External	ffl_nooutfile
	External	ffl_normal
	External	ffl_repopen
	External	ffl_repwrite
	External	fut_normal
	External	fut_error
	External	cut_normal

C
C  Open the report file.
C
	rstatus = fut_get_lun(fut_report_lun)
	If (rstatus .Ne. %Loc(fut_normal)) Then
	   Call LIB$Signal (%Val(rstatus))
	EndIf

	Open (unit=fut_report_lun, file=ffli.report_file, status='new',
     &            iostat=io_stat)

	If (io_stat .Eq. 0) Then
	   status = %Loc(ffl_normal)
	Else
	   status = %Loc(ffl_repopen)
	   Call LIB$Signal (ffl_repopen, %Val(2),
     &			    ffli.report_file(1:ffli.replen), %Val(io_stat))
	EndIf


	If (status .Eq. %Loc(ffl_normal)) Then
C
C  Initialize the report.
C

c
c  Initialize the error handler.
c
	   Call LIB$Establish(fut_error)

c
c  Get the the username of the account invoking the program.
c
	   Call LIB$Getjpi (jpi$_username,,,,username,)

c
c  Write the banner and the account and invocation information.
c
	   tstatus = cut_register_version (version)
	   tstatus = cut_display_banner (fut_report_lun, 80,
     &					'FIRAS Facility FFL_Fishinput_Long')
	   Write (fut_report_lun,10,iostat=io_stat)
	   Write (fut_report_lun,20,iostat=io_stat) username
	   Write (fut_report_lun,30,iostat=io_stat) current_gmt(1:11)
	   Write (fut_report_lun,40,iostat=io_stat) ffli.scan_mode
	   Write (fut_report_lun,50,iostat=io_stat) ffli.jstart_time(1:11),
     &						    ffli.jstop_time(1:11)
	   Write (fut_report_lun,60,iostat=io_stat) ffli.npts
	   Write (fut_report_lun,70,iostat=io_stat)
     &						 ffli.report_file(1:ffli.replen)
  10	   Format (31x, 'Processing Report', //)
  20	   Format (x, 'Run by:                     ', a)
  30	   Format (x, 'Run time:                   ', a)
  40	   Format (x, 'Channel / Scan Mode:        ', a)
  50	   Format (x, 'Timerange:                  ', a, ' - ', a)
  60	   Format (x, 'Spectrum length:            ', I3, ' points')
  70	   Format (x, 'Processing Report File:     ', a)
	   If (io_stat .Ne. 0) Then
	      status = %Loc(ffl_repwrite)
	      Call LIB$Signal (ffl_repwrite, %Val(2),
     &			       ffli.report_file(1:ffli.replen), %Val(io_stat))
	   EndIf

	EndIf	!	(status from open


	If (status .Eq. %Loc(ffl_normal)) Then
C
C  Translate the software logical names.
C

c
c  Set up the logical names to be translated.
c
	   If (ffli.hybrid .Eq. fac_present) Then
	      n_names = 4
	      log_name(1) = 'CSDR$FIRAS_IN1'
	      log_name(2) = 'CSDR$FIRAS_IN2'
	      log_name(3) = 'CSDR$FIRAS_OUT'
	      log_name(4) = 'CSDR$FIRAS_REF'
	   Else
	      n_names = 3
	      log_name(1) = 'CSDR$FIRAS_IN1'
	      log_name(2) = 'CSDR$FIRAS_OUT'
	      log_name(3) = 'CSDR$FIRAS_REF'
	   EndIf

c
c  Get the translation.
c
	   Do j = 1,n_names
	      tstatus = cut_translate_archive_id (log_name(j), temp_log,
     &				      temp_len, trans_name(j), trans_len(j))
	      if (tstatus .Ne. %Loc(cut_normal)) Call LIB$Signal (%Val(tstatus))
	      Call STR$Trim (trans_name(j), trans_name(j), trans_len(j))
	   EndDo

c
c  Write out the translated logicals.
c
	   Write (fut_report_lun,80,iostat=io_stat)
	   Do j = 1,n_names
	      Write (fut_report_lun, 90,iostat=io_stat) log_name(j)(1:14),
     &						   trans_name(j)(1:trans_len(j))
	   EndDo
  80	   Format (//, x, 'Logical Name Translations:', /)
  90	   Format (4x, a, 11x, a)
	   If (io_stat .Ne. 0) Then
	      status = %Loc(ffl_repwrite)
	      Call LIB$Signal (ffl_repwrite, %Val(2),
     &			       ffli.report_file(1:ffli.replen), %Val(io_stat))
	   EndIf

	EndIf	!	(status from open and writes


	If (status .Eq. %Loc(ffl_normal)) Then
C
C  Write out the command line invocation
C

	   Write (fut_report_lun, 100,iostat=io_stat)
	   Do j = 1,ncmd
	      Write (fut_report_lun, 110,iostat=io_stat) cmd_line(j)
	   EndDo
 100	   Format (//, x, 'Command Line Invocation:', /)
 110	   Format (x, a)
	   If (io_stat .Ne. 0) Then
	      status = %Loc(ffl_repwrite)
	      Call LIB$Signal (ffl_repwrite, %Val(2),
     &			       ffli.report_file(1:ffli.replen), %Val(io_stat))
	   EndIf

	EndIf	!	(status from open and writes


	If (status .Eq. %Loc(ffl_normal)) Then
C
C  Report processing errors.
C

c
c  If return status from the command line parse is bad:
c
	   If (parse_status .Ne. %Loc(ffl_normal)) Then
	      Write (fut_report_lun,120,iostat=io_stat)
	      If (parse_status .Eq. %Loc(ffl_invaltime)) Then
	         Write (fut_report_lun,130,iostat=io_stat) ffli.jstop_time,
     &							   ffli.jstart_time
	      ElseIf (parse_status .Eq. %Loc(ffl_nooutfile)) Then
	         Write (fut_report_lun,140,iostat=io_stat)
	      EndIf
	   EndIf
 120	   Format (//, x, 'Error from FFL_PARSE:', /)
 130	   Format (4x, 'Jstop time ', a, ' is less than Jstart time ', a, '.')
 140	   Format (4x, 'No output file name extension specified.')
	   If (io_stat .Ne. 0) Then
	      status = %Loc(ffl_repwrite)
	      Call LIB$Signal (ffl_repwrite, %Val(2),
     &			       ffli.report_file(1:ffli.replen), %Val(io_stat))
	   EndIf

	EndIf	!	(status from open and writes


	FFL_Initialize_Report = status

	Return
	End
