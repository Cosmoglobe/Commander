	Integer * 4 Function  FFL_Parse (current_gmt, ncmd, cmd_line, cmdlen,
     &                                   ffli)

c-------------------------------------------------------------------------------
c
c	Function FFL_PARSE
c
c	This function parses the FFL command line, identifying qualifiers and
c	keywords.  Invocation flags are set the in the include file FFL_INVOC.
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
c		None
c
c	Output:
c		current_gmt	Character * 14	! GMT time of invocation
c		ncmd		Integer * 4	! number of command lines in
c						!   invocation
c		cmd_line(3)	Character * 79	! command line invocation
c		cmdlen(3)	Integer * 4	! length of command lines in
c						!   invocation
c		ffli		Record Structure defined in FFL_Invoc.txt
c
c	Subroutines called:
c		ct_binary_to_gmt
c		ct_gmt_to_binary
c		time_lt
c		upm_get_float
c		upm_get_longword
c		upm_get_value
c		upm_present
c		LIB$Signal
c		STR$UpCase
c		STR$Trim
c		SYS$GetTim
c
c	Include files:
c		ffl_invoc.txt
c		fut_params.txt
c		upm_stat_msg.txt
c		$ssdef
c-------------------------------------------------------------------------------
c   Changes:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11690
c
c	Added apodization qualifier      Alice Trenholme, GSC, 28 December 1994
c
c	Name of facility changed from FFI to FFL.  Fred Shuman, HSTX
c	1995 May 19.
c
c	Changed FL from scan mode 5 to 6; added scan mode 5: FS.
c	Fred Shuman, HSTX, 1995 June 9.
c
c	In order to remove some of the obscurity arising from Include files that
c	   harbor hidden Common blocks, converted these Commons to Structures.
c	   This necessitates adding their Record names to the calling lists of
c	   functions that use their variables.
c	Fred Shuman, HSTX, 1995 June 14.
c-------------------------------------------------------------------------------

	Implicit None

	Include '($ssdef)'
	Include '(fut_params)'		! defines FIRAS parameters, fac_*
	Include '(ffl_invoc)'		! defines record ffli (structure invoc)
					!   and parameter "version"
	Include '(upm_stat_msg)'
c
c  Call arguments:
c
	Character *14	current_gmt	!  GMT time of invocation
	Integer   * 4	ncmd		!  number of command lines in invocation
	Character *79	cmd_line(3)	!  command line invocation
	Integer   * 4	cmdlen(3)	!  length of command lines in invocation
cccccc	Record Structure ffli is defined in Include file ffl_invoc.txt
c
c  All other variables and functions:
c
	Character * 4	thresh_val	!  hybrid input threshold
	Character * 2	qual_val(2)	!  data quality thresholds (string)
	Integer   * 2	exlen		!  length of input file extension
	Integer   * 2	len		!  length of string from command line

	Integer   * 4	clen		!  length of command line
	Integer   * 4	current_time(2)	!  VAX ADT time of invocation
	Integer   * 4	j		!  a counter
	Integer   * 4	ljstart		!  length of start time string
	Integer   * 4	ljstop		!  length of stop time string
	Integer   * 4	parse_status	!  return status
	Integer   * 4	quality		!  data quality threshold (numeric)
	Integer   * 4	status		!  return status

	Real	  * 4	threshold	!  hybrid input threshold

	Integer   * 4	upm_get_float
	Integer   * 4	upm_get_longword
	Integer   * 4	upm_get_value
	Integer   * 4	upm_present
	Logical   * 1	time_lt

	External	ffl_invaltime
	External	ffl_invscnmode
	External	ffl_nooutfile
	External	ffl_normal

C
C  Initialize the parse.
C
	Call SYS$GetTim (current_time)		!  Get the time of invocation
	Call ct_binary_to_gmt (current_time, current_gmt)
	parse_status = %Loc(ffl_normal)
	ncmd = 1				!  Number of command lines
	cmd_line(ncmd)(1:3) = 'FFL'		!  Initial command line
	cmdlen(ncmd) = 3			!  Initial command length


C
C  Parse the command line.
C

c
c
c  Set the start timetag invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('JSTART')
	If (status .Eq. upm_pres) Then
	   status = upm_get_value ('JSTART', ffli.jstart_time, len)
	   If (status .Eq. upm_absent) Then
	      ffli.jstart_time = fac_jstart_default
	   EndIf
	   ljstart = index(ffli.jstart_time,' ') - 1
	   if (ljstart .Eq. -1) ljstart = 14
	   ffli.jstart_time = ffli.jstart_time(1:ljstart) //
     &			      fac_jstart_default(ljstart+1:)
	Else
	   ffli.jstart_time = fac_jstart_default
	EndIf
	cmd_line(ncmd)(clen+1:clen+19) = '/JSTART=' // ffli.jstart_time(1:11)
	cmdlen(ncmd) = cmdlen(ncmd) + 19

c
c  Set the stop timetag invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('JSTOP')
	If (status .Eq. upm_pres) Then
	   status = upm_get_value ('JSTOP', ffli.jstop_time, len)
	   If (status .Eq. upm_absent) Then
	      ffli.jstop_time = fac_jstop_default
	   EndIf
	   ljstop = index(ffli.jstop_time,' ') - 1
	   if (ljstop .Eq. -1) ljstop = 14
	   ffli.jstop_time  = ffli.jstop_time(1:ljstop) //
     &			      fac_jstop_default(ljstop+1:)
	Else
	   ffli.jstop_time = fac_jstop_default
	EndIf
	cmd_line(ncmd)(clen+1:clen+18) = '/JSTOP=' // ffli.jstop_time(1:11)
	cmdlen(ncmd) = cmdlen(ncmd) + 18

c
c  Set the data quality invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('QUALITY')
	If (status .Eq. upm_pres  .Or.  status .Eq. upm_defaulted) Then
	   ffli.instr_qual = 3
	   ffli.attit_qual = 32
	   status = upm_get_longword ('QUALITY.INSTRUMENT', quality)
	   if (status .Eq. ss$_normal) ffli.instr_qual = quality
	   status = upm_get_longword ('QUALITY.ATTITUDE', quality)
	   if (status .Eq. ss$_normal) ffli.attit_qual = quality
	EndIf
	Write(qual_val(1),20) ffli.instr_qual
	Write(qual_val(2),20) ffli.attit_qual
  20	Format (i2)
	cmd_line(ncmd)(clen+1:clen+36) = '/QUALITY=(INSTRUMENT=' // qual_val(1)
     &				      // ',ATTITUDE=' // qual_val(2) // ')'
	cmdlen(ncmd) = cmdlen(ncmd) + 36
	ncmd = ncmd + 1
	cmd_line(ncmd)(1:3) = '   '
	cmdlen(ncmd) = 3

c
c  Set the channel specifier invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('CHANNEL')
	If (status .Eq. upm_pres) Then
	   status = upm_present ('CHANNEL.RH')
	   if (status .Eq. upm_pres) ffli.chan = 1
	   status = upm_present ('CHANNEL.RL')
	   if (status .Eq. upm_pres) ffli.chan = 2
	   status = upm_present ('CHANNEL.LH')
	   if (status .Eq. upm_pres) ffli.chan = 3
	   status = upm_present ('CHANNEL.LL')
	   if (status .Eq. upm_pres) ffli.chan = 4
	Else
	   ffli.chan = 1
	EndIf
	cmd_line(ncmd)(clen+1:clen+11) = '/CHANNEL=' //
     &					   fac_channel_ids(ffli.chan)
	cmdlen(ncmd) = cmdlen(ncmd) + 11

c
c  Set the scan mode specifier invocation flags, and set the short flag where
c   appropriate.
c
	clen = cmdlen(ncmd)
	status = upm_present ('SCAN_MODE')
	If (status .Eq. upm_pres) Then
	   status = upm_present ('SCAN_MODE.SS')
	   If (status .Eq. upm_pres) Then
	      ffli.smode  = 1
	   EndIf
	   status = upm_present ('SCAN_MODE.SF')
	   If (status .Eq. upm_pres) Then
	      ffli.smode  = 2
	   EndIf
	   status = upm_present ('SCAN_MODE.LS')
	   If (status .Eq. upm_pres) Then
	      ffli.smode  = 3
	   EndIf
	   status = upm_present ('SCAN_MODE.LF')
	   If (status .Eq. upm_pres) Then
	      ffli.smode  = 4
	   EndIf
	   status = upm_present ('SCAN_MODE.FS')
	   If (status .Eq. upm_pres) Then
	      ffli.smode  = 5
	      If (ffli.chan .Eq. 1  .Or.  ffli.chan .Eq. 3) Then
	         parse_status = %Loc(ffl_invscnmode)
	         Call LIB$Signal(ffl_invscnmode)
	      EndIf
	   EndIf
	   status = upm_present ('SCAN_MODE.FL')
	   If (status .Eq. upm_pres) Then
	      ffli.smode  = 6
	      If (ffli.chan .Eq. 1  .Or.  ffli.chan .Eq. 3) Then
	         parse_status = %Loc(ffl_invscnmode)
	         Call LIB$Signal(ffl_invscnmode)
	      EndIf
	   EndIf
	Else
	   ffli.smode  = 1
	EndIf
	cmd_line(ncmd)(clen+1:clen+13) = '/SCAN_MODE=' //
     &					   fac_scan_mode_idsl(ffli.smode)
	cmdlen(ncmd) = cmdlen(ncmd) + 13
	ffli.scan_mode = fac_channel_ids(ffli.chan) //
     &			 fac_scan_mode_idsl(ffli.smode)
	ffli.npts = fac_spec_length(ffli.smode)

c
c  Set the hybrid input invocation flag.
c
	clen = cmdlen(ncmd)
	status = upm_present ('HYBRID')
	If (status .Eq. upm_pres) Then
	   ffli.hybrid = fac_present
	   cmd_line(ncmd)(clen+1:clen+7) = '/HYBRID'
	   cmdlen(ncmd) = cmdlen(ncmd) + 7
	Else
	   ffli.hybrid = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+9) = '/NOHYBRID'
	   cmdlen(ncmd) = cmdlen(ncmd) + 9
	EndIf

c
c  Set the hybrid input threshold value.
c
	If (ffli.hybrid .Eq. fac_present) Then
	   clen = cmdlen(ncmd)
	   status = upm_present ('THRESHOLD')
	   If (status .Eq. upm_pres  .Or.  status .Eq. upm_defaulted) Then
	      status = upm_get_float ('THRESHOLD', threshold)
	      if (status .Eq. ss$_normal) ffli.threshold = threshold
	      Write(thresh_val,30) ffli.threshold
	      cmd_line(ncmd)(clen+1:clen+15) = '/THRESHOLD=' // thresh_val
	      cmdlen(ncmd) = cmdlen(ncmd) + 15
	   EndIf
	EndIf
  30	Format (f4.1)

c
c  Set the output file extension invocation flag.
c
	clen = cmdlen(ncmd)
	status = upm_present ('FILE_EXT')
	If (status .Eq. upm_pres) Then
	   status = upm_get_value ('FILE_EXT', ffli.file_ext, exlen)
	   Call STR$UpCase (ffli.file_ext, ffli.file_ext)
	   cmd_line(ncmd)(clen+1:clen+10+exlen) = '/FILE_EXT=' //
     &						    ffli.file_ext(1:exlen)
	   cmdlen(ncmd) = cmdlen(ncmd) + 10 + exlen
	Else
	   parse_status = %Loc(ffl_nooutfile)
	   Call LIB$Signal (ffl_nooutfile)
	EndIf
	ncmd = ncmd + 1
	cmd_line(ncmd)(1:3) = '   '
	cmdlen(ncmd) = 3

c
c  Set the processing report invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('REPORT')
	If (status .Eq. upm_negated) Then
	   ffli.report = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+9) = '/NOREPORT'
	   cmdlen(ncmd) = cmdlen(ncmd) + 9
	Else
	   ffli.report = fac_present
	   status = upm_get_value ('REPORT', ffli.report_file, ffli.replen)
	   If (status .Eq. upm_absent) Then  !  Set default report file name
	      ffli.report_file = 'FFL_' // ffli.scan_mode //
     &				 ffli.file_ext(1:exlen) // '.' //
     &				'REP_' // current_gmt(1:9)
	      Call STR$Trim (ffli.report_file, ffli.report_file, ffli.replen)
	      cmd_line(ncmd)(clen+1:clen+8+ffli.replen) = '/REPORT=' //
     &						 ffli.report_file(1:ffli.replen)
	      cmdlen(ncmd) = cmdlen(ncmd) + 8 + ffli.replen
	   Else				     !  Report file name from invocation
	      Call STR$UpCase (ffli.report_file, ffli.report_file)
	      cmd_line(ncmd)(clen+1:clen+8+ffli.replen) = '/REPORT=' //
     &						 ffli.report_file(1:ffli.replen)
	      cmdlen(ncmd) = cmdlen(ncmd) + 8 + ffli.replen
	   EndIf
	EndIf


C
C  Define the timerange variables.
C
	ffli.time_range  = ffli.jstart_time // ';' //
     &			   ffli.jstop_time // ';'

	Call ct_gmt_to_binary (ffli.jstart_time, ffli.jstart)
	Call ct_gmt_to_binary (ffli.jstop_time, ffli.jstop)

	If (time_lt(ffli.jstop, ffli.jstart)) Then
	   parse_status = %Loc(ffl_invaltime)
	   Call LIB$Signal(ffl_invaltime, %Val(2), ffli.jstop_time,
     &						   ffli.jstart_time)
	EndIf


	FFL_Parse = parse_status

	Return
	End
