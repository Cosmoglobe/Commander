	integer * 4 function  fsl_parse (current_gmt, ncmd, cmd_line, cmdlen)

c-------------------------------------------------------------------------------
c
c	Function FSL_PARSE
c
c	This function parses the FSL command line, identifying qualifiers and
c	keywords.  Invocation flags are returned to the FSL caller program.
c
c       Author:
c                FCF_Parse
c                Gene Eplee
c		 General Sciences Corp.
c		 10 March 1993
c
c                FSL_Parse
c                Shirley M. Read
c                Hughes STX Corporation
c                12 July 1995
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
c		cmd_line(11)		character * 79		command line
c								invocation
c		cmdlen(11)		integer * 4		length of
c								command lines
c								in invocation
c
c	Subroutines called:
c		ct_binary_to_gmt
c		ct_gmt_to_binary
c		lib$signal
c		lib$sys_trnlog
c		str$upcase
c		str$trim
c		sys$gettim
c		time_lt
c		upm_get_longword
c		upm_get_value
c		upm_present
c
c	Include files:
c		$ssdef
c		fsl_invoc.txt
c		fut_params.txt
c		upm_stat_msg.txt
c
c-------------------------------------------------------------------------------
c
c	Changes for FCF:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11395
c
c	Add DVECTOR command line switch.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11397
c
c	Add DIFFERENTIAL command line switch.
c	Gene Eplee, GSC, 11 July 1994
c	SPR 11826
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, July 12, 1995
c       Modified FCF_Parse to FSL_Parse for the new FIRAS pipeline which
c       will process long spectra to get improved frequency resolution.
c           1. Added error checking and new FUT error messages for command
c              line input. Input is required for data type, channel, and
c              scan mode; no default is allowed.
c           2. FSL_Scan 3 is no longer allowed on the command line since
c              there is no corresponding calibration model solution. Code
c              relating to FSL_Scan 3 has been removed.
c           3. Some of the code used to set fcc type parameters relating to
c              the 5 allowed FSL_Scan types of coadded IFGs was modified to
c              meet new FSL requirements.
c           4. Set fcc parameter FFT length and spectra length according to
c              the scan mode of data being processed.
c           5. Added qualifier /FIL_VAR to use the all-sky variances from FIL.
c           6. Added qualifier /TSIG_ZERO to set the temperature sigmas to zero.
c
c	Fred Shuman, Hughes STX Corporation, 1995 Sep 15
c	    7. Relocated /TSIG_ZERO and /FIL_VAR qualifiers in order to group
c	       them more logically and so that they would fit better into text
c	       lines in the report file.
c
c	Fred Shuman, Hughes STX Corporation, 1995 Sep 27
c	    8. Changed  fac_scan_mode_ids  to  fac_scan_mode_idsL.
c-------------------------------------------------------------------------------

	implicit none

	include '($ssdef)'
	include '(fut_params)'
	include '(fsl_invoc)'
	include '(upm_stat_msg)'


	character * 79  cmd_line(11)	!  command line invocation
	character * 14	current_gmt	!  GMT time of invocation
	character *  4	pixno		!  pixel list entry
	character * 32	plot_device	!  PLT plot device for report
	character *  2	qual_val(2)	!  data quality thresholds
	character * 15  rep_mid		!  report file name middle part

	integer * 2	exlen		!  length of input file extension
	integer * 2	len		!  length of string from command line

	integer * 4	clen		!  length of command line
	integer * 4	cmdlen(11)	!  length of command lines in invocation
	integer * 4	current_time(2)	!  VAX ADT time of invocation
	integer * 4	j		!  a counter
	integer * 4	ljstart		!  length of start time string
	integer * 4	ljstop		!  length of stop time string
	integer * 4	ncmd		!  number of command lines in invocation
	integer * 4	parse_status	!  return status
	integer * 4	quality		!  data quality threshold
	integer * 4	status		!  return status

	integer * 4	upm_get_longword
	integer * 4	upm_get_value
	integer * 4	upm_present
	logical * 1     time_lt

	external	fsl_invaltime
	external	fsl_nomodfile
	external	fsl_nopixlist
	external	fsl_normal
	external	fsl_noinput
	external	fsl_nochan
	external	fsl_noscan

C
C  Initialize the parse.
C

	call sys$gettim (current_time)		!  Get the time of invocation
	call ct_binary_to_gmt (current_time, current_gmt)
	parse_status = %loc(fsl_normal)
	fcc_num_pix = 0				!  Set pixel list at zero
	fcc_xhsf = fac_not_present		!  Set the initial XHSF flag
						!  to trigger phase_shift
	ncmd = 1				!  Number of command lines
	cmd_line(ncmd)(1:3) = 'FSL'		!  Initial command line
	cmdlen(ncmd) = 3			!  Initial command length

C
C  Parse the command line.
C

c
c  Set the input file extension invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('FILE_EXT')
	if (status .eq. upm_pres) then
	   status = upm_get_value ('FILE_EXT', fcc_file_ext, exlen)
	   fcc_file = fac_present
	   call str$upcase (fcc_file_ext, fcc_file_ext)
	   cmd_line(ncmd)(clen+1:clen+10+exlen) = '/FILE_EXT=' //
     .						    fcc_file_ext(1:exlen)
	   cmdlen(ncmd) = cmdlen(ncmd) + 10 + exlen
	else
	   fcc_file = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+11) = '/NOFILE_EXT'
	   cmdlen(ncmd) = cmdlen(ncmd) + 11
	endif

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
	elseif (fcc_file .eq. fac_not_present) then
	   fcc_jstart_time = fac_jstart_default
	endif
	if (fcc_file .eq. fac_not_present) then
	   cmd_line(ncmd)(clen+1:clen+19) = '/JSTART=' // fcc_jstart_time(1:11)
	   cmdlen(ncmd) = cmdlen(ncmd) + 19
	endif

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
	elseif (fcc_file .eq. fac_not_present) then
	   fcc_jstop_time = fac_jstop_default
	endif
	if (fcc_file .eq. fac_not_present) then
	   cmd_line(ncmd)(clen+1:clen+18) = '/JSTOP=' // fcc_jstop_time(1:11)
	   cmdlen(ncmd) = cmdlen(ncmd) + 18
	endif

c
c  Set the input data type invocation flags. Data type required on command line.
c
	clen = cmdlen(ncmd)
	status = upm_present ('INPUT')
	if (status .eq. upm_pres) then
	   status = upm_present ('INPUT.SKY')
	   if (status .eq. upm_pres) then
	      fcc_sky = fac_present
	      fcc_cal = fac_not_present
	      fcc_data_type = 'SKY'
	      cmd_line(ncmd)(clen+1:clen+10) = '/INPUT=SKY'
	   endif
	   status = upm_present ('INPUT.CAL')
	   if (status .eq. upm_pres) then
	      fcc_cal = fac_present
	      fcc_sky = fac_not_present
	      fcc_data_type = 'CAL'
	      cmd_line(ncmd)(clen+1:clen+10) = '/INPUT=CAL'
	   endif
	else
	   parse_status = %loc(fsl_noinput)
	   call lib$signal(fsl_noinput)
	endif
	cmdlen(ncmd) = cmdlen(ncmd) + 10

c
c  Set the input pixel list invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('PIXEL')
	if (status .eq. upm_pres) then
	   fcc_pixel = fac_present
	   fcc_num_pix = 1
	   do while (status  .and.  (fcc_num_pix .le. fcc_max_pix))
	      status = upm_get_longword ('PIXEL', fcc_plist(fcc_num_pix))
	      if (status  .and.
     .		  (fcc_plist(fcc_num_pix) .ge. 0)  .and.
     .		  (fcc_plist(fcc_num_pix) .le. 6143)) then
	         fcc_num_pix = fcc_num_pix + 1
	      endif
	   enddo
	   fcc_num_pix = fcc_num_pix - 1
	   if (fcc_num_pix .eq. 0) then
	      parse_status = %loc(fsl_nopixlist)
	      call lib$signal(fsl_nopixlist)
	   elseif (fcc_num_pix .eq. 1) then
	      write (pixno,10) fcc_plist(1)
	      cmd_line(ncmd)(clen+1:clen+11) = '/PIXEL=' // pixno
	      cmdlen(ncmd) = cmdlen(ncmd) + 11
	   elseif (((fcc_num_pix .le. 5) .and. (fcc_file .eq. fac_present)) .or.
     .	      ((fcc_num_pix .le. 2) .and. (fcc_file .eq. fac_not_present))) then
	      cmd_line(ncmd)(clen+1:clen+8) = '/PIXEL=('
	      clen = clen + 8
	      do j = 1,fcc_num_pix-1
	         write (pixno,10) fcc_plist(j)
	         cmd_line(ncmd)(clen+1:clen+5) = pixno // ','
	         clen = clen + 5
	      enddo
	      write (pixno,10) fcc_plist(fcc_num_pix)
	      cmd_line(ncmd)(clen+1:clen+5) = pixno // ')'
	      cmdlen(ncmd) = clen + 5
	   else

	      ncmd = ncmd + 1
	      cmd_line(ncmd)(1:3) = '   '

	      cmd_line(ncmd)(4:11) = '/PIXEL=('
	      clen = 11
	      do j = 1,fcc_num_pix-1
	         write (pixno,10) fcc_plist(j)
	         cmd_line(ncmd)(clen+1:clen+5) = pixno // ','
	         clen = clen + 5
	         if (amod(floatj(j),13.0) .eq. 0.0) then
	            cmdlen(ncmd) = clen

	            ncmd = ncmd + 1
	            cmd_line(ncmd)(1:11) = '           '

	            clen = 11
	         endif
	      enddo
	      write (pixno,10) fcc_plist(fcc_num_pix)
	      cmd_line(ncmd)(clen+1:clen+5) = pixno // ')'
	      cmdlen(ncmd) = clen + 5
	   endif
	else
	   fcc_pixel = fac_not_present
	   if (fcc_sky .eq. fac_present) then
	      cmd_line(ncmd)(clen+1:clen+10) = '/PIXEL=ALL'
	      cmdlen(ncmd) = cmdlen(ncmd) + 10
	   endif
	endif

	ncmd = ncmd + 1
	cmd_line(ncmd)(1:3) = '   '
	cmdlen(ncmd) = 3
  10	format (I4)

c
c  Set the data quality invocation flags.
c
	clen = cmdlen(ncmd)
	fcc_instr_qual = fac_many_yellow
	fcc_attit_qual = fac_many_yellow
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
  20	format (i2)

c
c  Set the channel specifier invocation flags. Channel required on command line.
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
	   parse_status = %loc(fsl_nochan)
	   call lib$signal(fsl_nochan)
	endif
	cmd_line(ncmd)(clen+1:clen+11) = '/CHANNEL=' //
     .					   fac_channel_ids(fcc_chan)
	cmdlen(ncmd) = cmdlen(ncmd) + 11

c
c  Set the scan mode specifier invocation flags. Scan mode is required on
c  command line. Set the fcc flags for scan mode, MTM length and speed,
c  and adds per group accordingly. Also set the high channel SF data
c  phase shift flag.
c
	clen = cmdlen(ncmd)
	status = upm_present ('SCAN_MODE')
	if (status .eq. upm_pres) then
	   status = upm_present ('SCAN_MODE.SS')
	   if (status .eq. upm_pres) then
	      fcc_smode  = 1
	      fcc_length = 0
	      fcc_speed  = 0
	      fcc_ngroup = 3
	   endif
	   status = upm_present ('SCAN_MODE.SF')
	   if (status .eq. upm_pres) then
	      fcc_smode  = 2
	      fcc_length = 0
	      fcc_speed  = 1
	      fcc_ngroup = 2
	   endif
c
c  Scan Mode LS is not allowed as input to FSL because there is no
c  corresponding calibration solution.
c
	   status = upm_present ('SCAN_MODE.LF')
	   if (status .eq. upm_pres) then
	      fcc_smode  = 4
	      fcc_length = 1
	      fcc_speed  = 1
	      if ((fcc_chan .eq. 1) .or. (fcc_chan .eq. 3)) fcc_ngroup = 2
	      if ((fcc_chan .eq. 2) .or. (fcc_chan .eq. 4)) fcc_ngroup = 8
	   endif
	   status = upm_present ('SCAN_MODE.FS')
	   if (status .eq. upm_pres) then
	      fcc_smode  = 5
	      fcc_length = 0
	      fcc_speed  = 1
	      fcc_ngroup = 8
	   endif
	   status = upm_present ('SCAN_MODE.FL')
	   if (status .eq. upm_pres) then
	      fcc_smode  = 6
	      fcc_length = 0
	      fcc_speed  = 1
	      fcc_ngroup = 8
	   endif
	else
	   parse_status = %loc(fsl_noscan)
	   call lib$signal(fsl_noscan)
	endif
c
c  Set XHSF flag for short fast high channels.
c
	if (((fcc_chan .eq. 1) .or. (fcc_chan .eq. 3)) .and.
     .	   (fcc_length .eq. 0) .and. (fcc_speed .eq. 1)) fcc_xhsf = fac_present
	cmd_line(ncmd)(clen+1:clen+13) = '/SCAN_MODE=' //
     .					   fac_scan_mode_idsL(fcc_smode)
	cmdlen(ncmd) = cmdlen(ncmd) + 13
	fcc_scan_mode = fac_channel_ids(fcc_chan) //
     .			fac_scan_mode_idsL(fcc_smode)

	ncmd = ncmd + 1
	cmd_line(ncmd)(1:3) = '   '
	cmdlen(ncmd) = 3

c
c  Set the single ifg invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('SINGLE_IFG')
	if (status .eq. upm_pres) then
	   fcc_single = fac_present
	   cmd_line(ncmd)(clen+1:clen+11) = '/SINGLE_IFG'
	   cmdlen(ncmd) = cmdlen(ncmd) + 11
	else
	   fcc_single = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+13) = '/NOSINGLE_IFG'
	   cmdlen(ncmd) = cmdlen(ncmd) + 13
	endif

c
c  Set the calibrate spectra invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('CALIBRATE')
	if (status .eq. upm_negated) then
	   fcc_calibrate = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+12) = '/NOCALIBRATE'
	   cmdlen(ncmd) = cmdlen(ncmd) + 12
	else
	   fcc_calibrate = fac_present
	   cmd_line(ncmd)(clen+1:clen+10) = '/CALIBRATE'
	   cmdlen(ncmd) = cmdlen(ncmd) + 10
	endif

c
c  Set the calibration model file extension invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('MODEL_EXT')
	if (status .eq. upm_pres) then
	   status = upm_get_value ('MODEL_EXT', fcc_model_ext, len)
	   call str$upcase (fcc_model_ext, fcc_model_ext)
	   cmd_line(ncmd)(clen+1:clen+11+len) = '/MODEL_EXT=' //
     .						    fcc_model_ext(1:len)
	   cmdlen(ncmd) = cmdlen(ncmd) + 11 + len
	elseif ((status .eq. upm_absent)  .and.
     .		(fcc_calibrate .eq. fac_not_present)) then
	   cmd_line(ncmd)(clen+1:clen+12) = '/NOMODEL_EXT'
	   cmdlen(ncmd) = cmdlen(ncmd) + 12
	else
	   parse_status = %loc(fsl_nomodfile)
	   call lib$signal (fsl_nomodfile)
	endif

	ncmd = ncmd + 1
	cmd_line(ncmd)(1:3) = '   '
	cmdlen(ncmd) = 3

c
c  Set the D-Vector weighting invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('DVECTOR')
	if (status .eq. upm_pres) then
	   fcc_dvec = fac_present
	   cmd_line(ncmd)(clen+1:clen+8) = '/DVECTOR'
	   cmdlen(ncmd) = cmdlen(ncmd) + 8
	else
	   fcc_dvec = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+10) = '/NODVECTOR'
	   cmdlen(ncmd) = cmdlen(ncmd) + 10
	endif

c
c  Set the differential spectra invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('DIFFERENTIAL')
	if (status .eq. upm_pres) then
	   fcc_diff = fac_present
	   if (fcc_cal .eq. fac_present) then
	      fcc_data_type = 'DCL'
	   else
	      fcc_data_type = 'DSK'
	   endif
	   cmd_line(ncmd)(clen+1:clen+13) = '/DIFFERENTIAL'
	   cmdlen(ncmd) = cmdlen(ncmd) + 13
	else
	   fcc_diff = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+15) = '/NODIFFERENTIAL'
	   cmdlen(ncmd) = cmdlen(ncmd) + 15
	endif

c
c  Set the all-sky FIL_VAR weighting invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('FIL_VAR')
	if (status .eq. upm_pres) then
	   fcc_flv = fac_present
	   cmd_line(ncmd)(clen+1:clen+8) = '/FIL_VAR'
	   cmdlen(ncmd) = cmdlen(ncmd) + 8
	else
	   fcc_flv = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+10) = '/NOFIL_VAR'
	   cmdlen(ncmd) = cmdlen(ncmd) + 10
	endif
c
c  Set the flag to reset temperature sigmas to zero.
c
	clen = cmdlen(ncmd)
	status = upm_present ('TSIG_ZERO')
	if (status .eq. upm_pres) then
	   fcc_tsig0 = fac_present
	   cmd_line(ncmd)(clen+1:clen+10) = '/TSIG_ZERO'
	   cmdlen(ncmd) = cmdlen(ncmd) + 10
	else
	   fcc_tsig0 = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+12) = '/NOTSIG_ZERO'
	   cmdlen(ncmd) = cmdlen(ncmd) + 12
	endif

	ncmd = ncmd + 1
	cmd_line(ncmd)(1:3) = '   '
	cmdlen(ncmd) = 3
c
c  Set the plot model, spectra invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('PLOT_DEVICE')
	if (status .eq. upm_pres) then
	   status = upm_get_value ('PLOT_DEVICE', fcc_plot_device, len)
	   if (status .eq. upm_absent) then
	      fcc_plot_device = 'PLT_DEVICE'
	      call lib$sys_trnlog ('PLT_DEVICE', len, plot_device)
	      plot_device = '"' // plot_device(1:len) // '"'
	      call str$upcase (plot_device, plot_device)
	      len = len + 2
	   else
	      call str$upcase (fcc_plot_device, fcc_plot_device)
	      plot_device = fcc_plot_device
	   endif
	   fcc_plot = fac_present
	   cmd_line(ncmd)(clen+1:clen+13+len) = '/PLOT_DEVICE=' //
     .						 plot_device(1:len)
	   cmdlen(ncmd) = clen + 13 + len
	else
	   fcc_plot = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+14) = '/NOPLOT_DEVICE'
	   cmdlen(ncmd) = clen + 14
	endif

c
c  Set the PLT command file invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('PLT_FILE')
	if (status .eq. upm_pres) then
	   status = upm_get_value ('PLT_FILE', fcc_plt_com_file, len)
	   fcc_plt_com = fac_present
	   call str$upcase (fcc_plt_com_file, fcc_plt_com_file)
	   cmd_line(ncmd)(clen+1:clen+10+len) = '/PLT_FILE=' //
     .						  fcc_plt_com_file(1:len)
	   cmdlen(ncmd) = cmdlen(ncmd) + 10 + len
	else
	   fcc_plt_com = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+11) = '/NOPLT_FILE'
	   cmdlen(ncmd) = cmdlen(ncmd) + 11
	endif

	ncmd = ncmd + 1
	cmd_line(ncmd)(1:3) = '   '
	cmdlen(ncmd) = 3

c
c  Set the write spectra invocation flags.
c
	clen = cmdlen(ncmd)
	status = upm_present ('WRITE')
	if (status .eq. upm_pres) then
	   status = upm_present ('WRITE.VOLT')
	   if (status .eq. upm_pres) then
	      fcc_write_vs = fac_present
	      fcc_write_cs = fac_not_present
	      cmd_line(ncmd)(clen+1:clen+11) = '/WRITE=VOLT'
	      cmdlen(ncmd) = cmdlen(ncmd) + 11
	   endif
	   status = upm_present ('WRITE.CAL')
	   if (status .eq. upm_pres) then
	      fcc_write_vs = fac_not_present
	      if (fcc_diff .eq. fac_present) then
	         fcc_write_ds = fac_present
	         fcc_write_cs = fac_not_present
	      else
	         fcc_write_ds = fac_not_present
	         fcc_write_cs = fac_present
	      endif
	      cmd_line(ncmd)(clen+1:clen+10) = '/WRITE=CAL'
	      cmdlen(ncmd) = cmdlen(ncmd) + 10
	   endif
	   status = upm_present ('WRITE.BOTH')
	   if (status .eq. upm_pres) then
	      fcc_write_vs = fac_present
	      if (fcc_diff .eq. fac_present) then
	         fcc_write_ds = fac_present
	         fcc_write_cs = fac_not_present
	      else
	         fcc_write_ds = fac_not_present
	         fcc_write_cs = fac_present
	      endif
	      cmd_line(ncmd)(clen+1:clen+11) = '/WRITE=BOTH'
	      cmdlen(ncmd) = cmdlen(ncmd) + 11
	   endif
	else
	   fcc_write_vs = fac_not_present
	   fcc_write_ds = fac_not_present
	   fcc_write_cs = fac_not_present
	   cmd_line(ncmd)(clen+1:clen+8) = '/NOWRITE'
	   cmdlen(ncmd) = cmdlen(ncmd) + 8
	endif

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
	   status = upm_get_value ('REPORT', fcc_report_file, len)
	   if (status .eq. upm_absent) then  !  Set default report file name
	      if (fcc_file .eq. fac_present) then
	         rep_mid = fcc_file_ext(4:exlen)
	      else
	         rep_mid = fcc_jstart_time(1:7) // '_' // fcc_jstop_time(1:7)
	      endif
	      fcc_report_file = 'FSL_' // fcc_data_type // '_' //
     .				 fcc_scan_mode // '_' // rep_mid // '.' //
     .				'REP_' // current_gmt(1:9)
	      cmd_line(ncmd)(clen+1:clen+50) = '/REPORT=' //
     .						 fcc_report_file(1:42)
	      cmdlen(ncmd) = cmdlen(ncmd) + 50
	   else				     !  Report file name from invocation
	      call str$upcase (fcc_report_file, fcc_report_file)
	      cmd_line(ncmd)(clen+1:clen+8+len) = '/REPORT=' //
     .						    fcc_report_file(1:len)
	      cmdlen(ncmd) = cmdlen(ncmd) + 8 + len
	   endif
	   call str$trim (fcc_report_file, fcc_report_file, fcc_replen)
	endif


	if (fcc_file .eq. fac_not_present) then
C
C  Define the timerange variables.
C
	   fcc_time_range  = fcc_jstart_time // ';' //
     .			     fcc_jstop_time // ';'

	   call ct_gmt_to_binary (fcc_jstart_time, fcc_jstart)
	   call ct_gmt_to_binary (fcc_jstop_time, fcc_jstop)

	   if (time_lt(fcc_jstop,fcc_jstart)) then
	      parse_status = %loc(fsl_invaltime)
	      call lib$signal(fsl_invaltime, %val(2), fcc_jstop_time,
     .						      fcc_jstart_time)
	   end if
	endif
C
C  Define the fcc parameters for FFT length and spectra length accordong to
C  the scan mode of data being processed.
C
	fcc_fft_length = fac_fft_length(fcc_smode)
	fcc_spec_length = fac_spec_length(fcc_smode)

	fsl_parse = parse_status

	return
	end
