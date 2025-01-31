	integer * 4 function  frd_parse_modell ()

c-------------------------------------------------------------------------------
c
c	Function FRD_PARSE_MODELL
c
c	This function parses the command line for FRD_EXTRACT_MODELL.
c	This parse determines the channel and scan mode for which 
c	FRD_EXTRACT_MODELL is to be run, the file extension of the model
c	solution to be extracted, and the reference timetag for the model
c	solution.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  2 October 1992
c		  SER 8178
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		none
c
c	Subroutines called:
c		ct_gmt_to_binary
c		str$trim
c		str$upcase
c		upm_get_value
c		upm_present
c
c	Include files:
c		frd_model_configl.txt
c		frd_model_invoc.txt
c		fut_params.txt
c		upm_stat_msg.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 22 October 1993
c	SER 11394
c
c       Changed adds per group determination for LSF/LLF data to 
c       reflect fact LSF data is decimated and "funny" LLF data is
c       truncated.  Changed parsing of scan mode to include FS and FL.
c       Alice Trenholme, GSC, 12 December 1995, SPR 12282.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(frd_model_configl)'
	include '(frd_model_invoc)'
	include '(upm_stat_msg)'

	integer * 2	len		!  input string length

	integer * 4	parse_status	!  return status
	integer * 4	status		!  return status

	integer * 4	upm_get_value
	integer * 4	upm_present

	external	frd_noinfile
	external	frd_normal


C
C  Parse the command line.
C

	parse_status = %loc(frd_normal)

c
c  Set the channel specifier invocation flags.
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

c
c  Set the scan mode specifier invocation flags.
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
	   status = upm_present ('SCAN_MODE.FS')
	   if (status .eq. upm_pres) then
	      fcc_smode  = 5
	      fcc_length = 0
	      fcc_speed  = 1
	   endif
	   status = upm_present ('SCAN_MODE.FL')
	   if (status .eq. upm_pres) then
	      fcc_smode  = 6
	      fcc_length = 0
	      fcc_speed  = 1
	   endif
	else
	   fcc_smode  = 1
	   fcc_length = 0
	   fcc_speed  = 0
	endif

c
c  Set the input file extension invocation flags.
c
	status = upm_present ('FILE_EXT')
	if (status .eq. upm_pres) then
	   status = upm_get_value ('FILE_EXT', fcc_infile_ext, len)
	   call str$upcase (fcc_infile_ext, fcc_infile_ext)
	else
	   parse_status = %loc(frd_noinfile)
	   call lib$signal (frd_noinfile)
	endif
	call str$trim (fcc_infile_ext, fcc_infile_ext, len)

c
c  Set the reference timetag flag.
c
	ref_gmt_time = '90001000000000'
	status = upm_present ('FLIGHT')
	if (status .eq. upm_pres) then
	   ref_gmt_time = '90001000000000'
	endif
	status = upm_present ('INT')
	if (status .eq. upm_pres) then
	   ref_gmt_time = '89251000000000'
	endif
	call ct_gmt_to_binary(ref_gmt_start,ref_start)
	call ct_gmt_to_binary(ref_gmt_stop,ref_stop)
	call ct_gmt_to_binary(ref_gmt_time,ref_time)


C
C  Determine the correct adds per group
C
	if (fcc_chan .eq. 1  .or.  fcc_chan .eq. 3) then
	   if (fcc_speed .eq. 0) then
	      fcc_ngroup = 3
	   elseif (fcc_speed .eq. 1) then
	      fcc_ngroup = 2
	   endif
	elseif (fcc_chan .eq. 2  .or.  fcc_chan .eq. 4) then
	   if (fcc_smode .eq. 1) then
	      fcc_ngroup = 3
	   elseif (fcc_smode .eq. 2) then
	      fcc_ngroup = 2
	   elseif (fcc_smode .eq. 3) then
	      fcc_ngroup = 12
	   elseif (fcc_smode .eq. 4) then
	      fcc_ngroup = 8
	   elseif (fcc_smode .eq. 5) then
	      fcc_ngroup = 8             !"FS" LSF data is decimated ....
	   elseif (fcc_smode .eq. 6) then
	      fcc_ngroup = 8
	   endif
	endif


	frd_parse_modell = parse_status

	return
	end
