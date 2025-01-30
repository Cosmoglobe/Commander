	integer * 4 function  ffi_parse (current_gmt, ncmd, cmd_line, cmdlen)

c-------------------------------------------------------------------------------
c
c	Function FFI_PARSE
c
c	This function parses the FFI command line, identifying qualifiers and
c	keywords.  Invocation flags are set the in the include file FFI_INVOC.
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
c		none
c
c	Output:
c		current_gmt		character * 14		GMT time of
c								invocation
c		ncmd			integer * 4		number of 
c								command lines
c								in invocation
c		cmd_line(3)		character * 79		command line
c								invocation
c		cmdlen(3)		integer * 4		length of 
c								command lines
c								in invocation
c
c	Subroutines called:
c		ct_binary_to_gmt
c		ct_gmt_to_binary
c		lib$signal
c		str$upcase
c		str$trim
c		sys$gettim
c		time_lt
c		upm_get_float
c		upm_get_longword
c		upm_get_value
c		upm_present
c
c	Include files:
c		$ssdef
c		ffi_invoc.txt
c		fut_params.txt
c		upm_stat_msg.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11690
c
c-------------------------------------------------------------------------------

	implicit none

	include '($ssdef)'
	include '(fut_params)'
	include '(ffi_invoc)'
	include '(upm_stat_msg)'


	character * 79  cmd_line(3)	!  command line invocation
	character * 14	current_gmt	!  GMT time of invocation
	character *  2	qual_val(2)	!  data quality thresholds
	character *  4	thresh_val	!  hybrid input threshold

	integer * 2	exlen		!  length of input file extension
	integer * 2	len		!  length of string from command line

	integer * 4	clen		!  length of command line
	integer * 4	cmdlen(3)	!  length of command lines in invocation
	integer * 4	current_time(2)	!  VAX ADT time of invocation
	integer * 4	j		!  a counter
	integer * 4	ljstart		!  length of start time string
	integer * 4	ljstop		!  length of stop time string
	integer * 4	ncmd		!  number of command lines in invocation
	integer * 4	parse_status	!  return status
	integer * 4	quality		!  data quality threshold
	integer * 4	status		!  return status

	real	* 4	threshold	!  hybrid input threshold

	integer * 4	upm_get_float
	integer * 4	upm_get_longword
	integer * 4	upm_get_value
	integer * 4	upm_present
	logical * 1     time_lt

	external	ffi_invaltime
	external	ffi_nooutfile
	external	ffi_normal

C
C  Initialize the parse.
C
	call sys$gettim (current_time)		!  Get the time of invocation
	call ct_binary_to_gmt (current_time, current_gmt)
	parse_status = %loc(ffi_normal)
	fcc_xhsf = fac_not_present		!  Set the initial XHSF flag
	fcc_xlsf = fac_not_present		!  Set the initial XLSF flag
	fcc_xllf = fac_not_present		!  Set the initial XLLF flag
	ncmd = 1				!  Number of command lines  
	cmd_line(ncmd)(1:3) = 'FFI'		!  Initial command line
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
	if (status .eq. upm_pres) then
	   status = upm_get_value ('JSTART', fcc_jstart_time, len)
	   if (status .eq. upm_absent) then
	      fcc_jstart_time = fac_jstart_default
	   endif
	   ljstart = index(fcc_jstart_time,' ') - 1
	   if (ljstart .eq. -1) ljstart = 14
	   fcc_jstart_time = fcc_jstart_time(1:ljstart) //
     .			     fac_jstart_default(ljstart+1:)
	else
	   fcc_jstart_time = fac_jstart_default
	endif
	cmd_line(ncmd)(clen+1:clen+19) = '/JSTART=' // fcc_jstart_time(1:11)
	cmdlen(ncmd) = cmdlen(ncmd) + 19

c
c  Set the stop timetag invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('JSTOP')
	if (status .eq. upm_pres) then
	   status = upm_get_value ('JSTOP', fcc_jstop_time, len)
	   if (status .eq. upm_absent) then
	      fcc_jstop_time = fac_jstop_default
	   endif
	   ljstop = index(fcc_jstop_time,' ') - 1
	   if (ljstop .eq. -1) ljstop = 14
	   fcc_jstop_time  = fcc_jstop_time(1:ljstop) //
     .			     fac_jstop_default(ljstop+1:)
	else
	   fcc_jstop_time = fac_jstop_default
	endif
	cmd_line(ncmd)(clen+1:clen+18) = '/JSTOP=' // fcc_jstop_time(1:11)
	cmdlen(ncmd) = cmdlen(ncmd) + 18

c
c  Set the data quality invocation flags.
c
	clen = cmdlen(ncmd)
	fcc_instr_qual = fac_many_yellow
	fcc_attit_qual = fac_many_y_some_r
	status = upm_present ('QUALITY')
	if (status .eq. upm_pres) then
	   status = upm_get_longword ('QUALITY.INSTRUMENT', quality)
	   if (status .eq. ss$_normal) fcc_instr_qual = quality
	   status = upm_get_longword ('QUALITY.ATTITUDE', quality)
	   if (status .eq. ss$_normal) fcc_attit_qual = quality
	endif
	write(qual_val(1),20) fcc_instr_qual
	write(qual_val(2),20) fcc_attit_qual
	cmd_line(ncmd)(clen+1:clen+36) = '/QUALITY=(INSTRUMENT=' // qual_val(1) 
     .				      // ',ATTITUDE=' // qual_val(2) // ')'
	cmdlen(ncmd) = cmdlen(ncmd) + 36
	ncmd = ncmd + 1
	cmd_line(ncmd)(1:3) = '   '
	cmdlen(ncmd) = 3
  20	format (i2)

c
c  Set the channel specifier invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('CHANNEL')
	if (status .eq. upm_pres) then
	   status = upm_present ('CHANNEL.RH')
	   if (status .eq. upm_pres) fcc_chan = 1
	   status = upm_present ('CHANNEL.RL')
	   if (status .eq. upm_pres) fcc_chan = 2
	   status = upm_present ('CHANNEL.LH')
	   if (status .eq. upm_pres) fcc_chan = 3
	   status = upm_present ('CHANNEL.LL')
	   if (status .eq. upm_pres) fcc_chan = 4
	else
	   fcc_chan = 1
	endif
	cmd_line(ncmd)(clen+1:clen+11) = '/CHANNEL=' //
     .					   fac_channel_ids(fcc_chan)
	cmdlen(ncmd) = cmdlen(ncmd) + 11

c
c  Set the scan mode specifier invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('SCAN_MODE')
	if (status .eq. upm_pres) then
	   status = upm_present ('SCAN_MODE.SS')
	   if (status .eq. upm_pres) then
	      fcc_smode  = 1
	      fcc_length = 0
	      fcc_speed  = 0
	   endif
	   status = upm_present ('SCAN_MODE.SF')
	   if (status .eq. upm_pres) then
	      fcc_smode  = 2
	      fcc_length = 0
	      fcc_speed  = 1
	   endif
	   status = upm_present ('SCAN_MODE.LS')
	   if (status .eq. upm_pres) then
	      fcc_smode  = 3
	      fcc_length = 1
	      fcc_speed  = 0
	   endif
	   status = upm_present ('SCAN_MODE.LF')
	   if (status .eq. upm_pres) then
	      fcc_smode  = 4
	      fcc_length = 1
	      fcc_speed  = 1
	   endif
	   status = upm_present ('SCAN_MODE.FL')
	   if (status .eq. upm_pres) then
	      fcc_smode  = 5
	      fcc_length = 0
	      fcc_speed  = 1
	   endif
	else
	   fcc_smode  = 1
	   fcc_length = 0
	   fcc_speed  = 0
	endif
	cmd_line(ncmd)(clen+1:clen+13) = '/SCAN_MODE=' //
     .					   fac_scan_mode_ids(fcc_smode)
	cmdlen(ncmd) = cmdlen(ncmd) + 13
	fcc_scan_mode = fac_channel_ids(fcc_chan) //
     .			fac_scan_mode_ids(fcc_smode)

c
c  Set the hybrid input invocation flag.
c
	clen = cmdlen(ncmd)
	status = upm_present ('HYBRID')
	if (status .eq. upm_pres) then
	   fcc_hybrid = fac_present
	   cmd_line(ncmd)(clen+1:clen+7) = '/HYBRID'
	   cmdlen(ncmd) = cmdlen(ncmd) + 7
	else
	   fcc_hybrid = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+9) = '/NOHYBRID'
	   cmdlen(ncmd) = cmdlen(ncmd) + 9
	endif

c
c  Set the hybrid input threshold value.
c
	if (fcc_hybrid .eq. fac_present) then
	   clen = cmdlen(ncmd)
	   fcc_threshold = 10.0
	   status = upm_present ('THRESHOLD')
	   if (status .eq. upm_pres) then
	      status = upm_get_float ('THRESHOLD', threshold)
	      if (status .eq. ss$_normal) fcc_threshold = threshold
	      write(thresh_val,30) fcc_threshold
	      cmd_line(ncmd)(clen+1:clen+15) = '/THRESHOLD=' // thresh_val
	      cmdlen(ncmd) = cmdlen(ncmd) + 15
	   endif
	endif
  30	format (f4.1)

c
c  Set the output file extension invocation flag.
c
	clen = cmdlen(ncmd)
	status = upm_present ('FILE_EXT')
	if (status .eq. upm_pres) then
	   status = upm_get_value ('FILE_EXT', fcc_file_ext, exlen)
	   call str$upcase (fcc_file_ext, fcc_file_ext)
	   cmd_line(ncmd)(clen+1:clen+10+exlen) = '/FILE_EXT=' //
     .						    fcc_file_ext(1:exlen)
	   cmdlen(ncmd) = cmdlen(ncmd) + 10 + exlen
	else
	   parse_status = %loc(ffi_nooutfile)
	   call lib$signal (ffi_nooutfile)
	endif
	ncmd = ncmd + 1
	cmd_line(ncmd)(1:3) = '   '
	cmdlen(ncmd) = 3

c
c  Set the processing report invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('REPORT')
	if (status .eq. upm_negated) then
	   fcc_report = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+9) = '/NOREPORT'
	   cmdlen(ncmd) = cmdlen(ncmd) + 9
	else
	   fcc_report = fac_present
	   status = upm_get_value ('REPORT', fcc_report_file, fcc_replen)
	   if (status .eq. upm_absent) then  !  Set default report file name
	      fcc_report_file = 'FFI_' // fcc_scan_mode // '_' //
     .				 fcc_file_ext(1:exlen) // '.' //
     .				'REP_' // current_gmt(1:9)
	      call str$trim (fcc_report_file, fcc_report_file, fcc_replen)
	      cmd_line(ncmd)(clen+1:clen+8+fcc_replen) = '/REPORT=' //
     .						   fcc_report_file(1:fcc_replen)
	      cmdlen(ncmd) = cmdlen(ncmd) + 8 + fcc_replen
	   else				     !  Report file name from invocation
	      call str$upcase (fcc_report_file, fcc_report_file)
	      cmd_line(ncmd)(clen+1:clen+8+fcc_replen) = '/REPORT=' //
     .						   fcc_report_file(1:fcc_replen)
	      cmdlen(ncmd) = cmdlen(ncmd) + 8 + fcc_replen
	   endif
	endif


C
C  Define the timerange variables.
C
	fcc_time_range  = fcc_jstart_time // ';' //
     .			  fcc_jstop_time // ';'

	call ct_gmt_to_binary (fcc_jstart_time, fcc_jstart)
	call ct_gmt_to_binary (fcc_jstop_time, fcc_jstop)

	if (time_lt(fcc_jstop,fcc_jstart)) then
	   parse_status = %loc(ffi_invaltime)
	   call lib$signal(ffi_invaltime, %val(2), fcc_jstop_time,
     .						   fcc_jstart_time)
	end if


C
C  Determine the correct adds per group
C
	if ((fcc_chan .eq. 1)  .or.  (fcc_chan .eq. 3)) then
	   if (fcc_speed .eq. 0) then
	      fcc_ngroup = 3
	   elseif (fcc_speed .eq. 1) then
	      fcc_ngroup = 2
	      if (fcc_length .eq. 0) fcc_xhsf = fac_present  ! Set the XHSF flag
	   endif
	elseif ((fcc_chan .eq. 2)  .or.  (fcc_chan .eq. 4)) then
	   if (fcc_smode .eq. 1) then
	      fcc_ngroup = 3
	   elseif (fcc_smode .eq. 2) then
	      fcc_ngroup = 2
	      fcc_xlsf = fac_present  	! Set the XLSF flag
	   elseif (fcc_smode .eq. 3) then
	      fcc_ngroup = 12
	   elseif (fcc_smode .eq. 4) then
	      fcc_ngroup = 8
	   elseif (fcc_smode .eq. 5) then
	      fcc_ngroup = 2
	      fcc_xllf   = fac_present	!  Set the XLLF flag
	   endif
	endif


	ffi_parse = parse_status

	return
	end
