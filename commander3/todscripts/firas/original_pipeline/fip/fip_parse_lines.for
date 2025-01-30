	integer * 4 function  fip_parse_lines (current_gmt)

c-------------------------------------------------------------------------------
c
c	Function FIP_PARSE_LINES
c
c	This function parses the command line for FIP_LINES.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  7 June 1993
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		current_gmt	character * 14		invocation time
c
c	Subroutines called:
c		ct_binary_to_gmt
c		ct_gmt_to_binary
c		sys$gettim
c		str$upcase
c		upm_get_longword
c		upm_get_value
c		upm_present
c
c	Incude files:
c		$ssdef
c		fip_config_freq.txt
c		fip_frequency.txt
c		fip_invoc_lines.txt
c		fut_params.txt
c		upm_stat_msg.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Replace FIP_CONFIG_LINES.TXT by FIP_CONFIG_FREQ.TXT.
c	Added FCC_FCHAN and FCC_FSMODE to FIP_FREQUENCY.TXT
c	Gene Eplee, GSC, 10 February 1994.
c
c-------------------------------------------------------------------------------

	implicit none

	include '($ssdef)'
	include '(fut_params)'
	include '(fip_invoc_lines)'
	include '(fip_config_freq)'
	include '(fip_frequency)'
	include '(upm_stat_msg)'

	character * 14	current_gmt		!  GMT time of invocation
	character * 20  rep_mid			!  report file name middle part

	integer * 2	midlen			!  rep_mid string length

	integer * 4	bar			!  index of '_' in string
	integer * 4	current_time(2)		!  VAX ADT time of invocation
	integer * 4	freq			!  input frequency cutoff
	integer * 4	pstatus			!  parse return status
	integer * 4	status			!  return status

	integer * 4	upm_get_longword
	integer * 4	upm_get_value
	integer * 4	upm_present

	external	fip_nolines
	external	fip_noreftime
	external	fip_normal

C
C  Parse the command line.
C

c
c  Get the time of the invocation
c
	call sys$gettim (current_time)
	call ct_binary_to_gmt (current_time, current_gmt)

c
c  Set the line file extension invocation flag.
c
	status = upm_present ('LINES_EXT')
	if (status .eq. upm_pres) then
	   status = upm_get_value ('LINES_EXT', fcc_lines_ext, fcc_extlen)
	   call str$upcase (fcc_lines_ext, fcc_lines_ext)
	   pstatus = %loc(fip_normal)
	else
	   pstatus = %loc(fip_nolines)
	   call lib$signal (fip_nolines)
	endif

c
c  Set the channel specifier invocation flag.
c
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
	fcc_fchan = fcc_chan

c
c  Set the scan mode specifier invocation flag.
c
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
	else
	   fcc_smode  = 1
	   fcc_length = 0
	   fcc_speed  = 0
	endif
	fcc_fsmode = fcc_smode
	fcc_scan_mode = fac_channel_ids(fcc_chan) //
     .			fac_scan_mode_ids(fcc_smode)

c
c  Set the reference timetag flag.
c
	fcc_flight = fac_endoffile
	ref_gmt_time = '90001000000000'
	status = upm_present ('FLIGHT')
	if (status .eq. upm_pres) then
	   fcc_flight = fac_present
	   ref_gmt_time = '90001000000000'
	endif
	status = upm_present ('INT')
	if (status .eq. upm_pres) then
	   fcc_flight = fac_not_present
	   ref_gmt_time = '89251000000000'
	endif
	if ((fcc_flight .eq. fac_present)  .or.
     .	     (fcc_flight .eq. fac_not_present)) then
	   call ct_gmt_to_binary(ref_gmt_start,ref_start)
	   call ct_gmt_to_binary(ref_gmt_stop,ref_stop)
	   call ct_gmt_to_binary(ref_gmt_time,ref_time)
	else
	   pstatus = %loc(fip_noreftime)
	   call lib$signal (fip_noreftime)
	endif

c
c  Set the frequency range invocation flags.
c
	status = upm_present ('FREQ_RANGE')
	if (status .eq. upm_pres) then
	   status = upm_present ('FREQ_RANGE.ALL')
	   if (status .eq. upm_pres) then
	      fcc_freq = fac_not_present
	   else
	      fcc_freq = fac_present
	      status = upm_present ('FREQ_RANGE.LOW')
	      if (status .eq. upm_pres) then
	         status = upm_get_longword ('FREQ_RANGE.LOW', freq)
	         if (status .eq. ss$_normal) fcc_lofreq = freq
	      endif
	      status = upm_present ('FREQ_RANGE.HIGH')
	      if (status .eq. upm_pres) then
	         status = upm_get_longword ('FREQ_RANGE.HIGH', freq)
	         if (status .eq. ss$_normal) fcc_hifreq = freq
	      endif
	   endif
	else
	   fcc_freq = fac_not_present
	   fcc_lofreq = fcc_lofreq_default
	   fcc_hifreq = fcc_hifreq_default
	endif

c
c  Set the processing report invocation flags.
c
	status = upm_present ('REPORT')
	if (status .eq. upm_negated) then
	   fcc_report = fac_not_present
	else
	   fcc_report = fac_present
	   status = upm_get_value ('REPORT', fcc_report_file, fcc_replen)
	   if (status .eq. upm_absent) then  !  Set default report file name
	      bar = index(fcc_lines_ext, '_')
	      rep_mid = fcc_lines_ext(bar+1:fcc_extlen)
	      call str$trim (rep_mid, rep_mid, midlen)
	      fcc_report_file = 'FIP_LINES_' // fcc_scan_mode // '_' //
     .				 rep_mid(1:midlen) // '.' // 'REP_' //
     .				 current_gmt(1:9)
	   else				     !  Report file name from invocation
	      call str$upcase (fcc_report_file, fcc_report_file)
	   endif
	   call str$trim (fcc_report_file, fcc_report_file, fcc_replen)
	endif


	fip_parse_lines = pstatus

	return
	end
