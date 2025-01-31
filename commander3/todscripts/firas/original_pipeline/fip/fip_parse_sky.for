	integer * 4 function  fip_parse_sky (current_gmt)

c-------------------------------------------------------------------------------
c
c	Function FIP_PARSE_SKY
c
c	This function parses the command line for FIP_SKY.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  16 June 1993
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
c		index
c		sys$gettim
c		str$upcase
c		upm_get_float
c		upm_get_longword
c		upm_get_value
c		upm_present
c
c	Incude files:
c		$ssdef
c		fip_frequency.txt
c		fip_invoc_sky.txt
c		fut_params.txt
c		upm_stat_msg.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Fixed bug in report file name.
c	Gene Eplee, GSC, 16 June 1993
c
c	Change input datatype, channel, scan mode parse to input filename parse.
c	Gene Eplee, GSC, 10 February 1994
c
c	Set new FCC_DIFF flag to FAC_NOT_PRESENT for FCS and FMS files.
c	Gene Eplee, GSC, 7 July 1994
c
c	Set new FCC_NUM_TYPE flag to FAC_PRESENT for FCS and FMS files.
c	Gene Eplee, GSC, 27 July 1994
c
c       Set FCC_DATA_TYPE to 1 for FMS data merging channels and scan mode
c       (HIGH, LOWF, LRES, or HRES), to 2 for FMS data merging scan modes
c       (RHFA, RLFA, LHFA, LLFA), and to 1 for FCS data.
c-------------------------------------------------------------------------------

	implicit none

	include '($ssdef)'
	include '(fut_params)'
	include '(fip_invoc_sky)'
	include '(fip_frequency)'
	include '(upm_stat_msg)'

	character * 14	current_gmt		!  GMT time of invocation
	character * 20	file_ext		!  input file extension
	character *  3  input			!  input data type
	character * 15  rep_mid			!  report file name middle part

	integer * 2	len			!  input string length
	integer * 2	extlen			!  file extension length
	integer * 2	midlen			!  rep_mid string length

	integer * 4	bar			!  index of '_' in string
	integer * 4	current_time(2)		!  VAX ADT time of invocation
	integer * 4	freq			!  input frequency cutoff
	integer * 4	period			!  index of '.' in string
	integer * 4	pstatus			!  parse return status
	integer * 4	status			!  return status

	real	* 4	glat			!  input latitude cutoff

	integer * 4	upm_get_float
	integer * 4	upm_get_longword
	integer * 4	upm_get_value
	integer * 4	upm_present

	external	fip_galexc
	external	fip_invalfile
	external	fip_nofile
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
c  Get the input file name.
c
	status = upm_present ('FILE')
	if (status .eq. upm_pres) then
	   status = upm_get_value ('FILE', fcc_infile, fcc_inlen)
	   call str$upcase (fcc_infile, fcc_infile)
c
c  Verify the input data type.
c
	   input = fcc_infile(1:3)
	   if (input .eq. 'FCS') then
	      pstatus = %loc(fip_normal)
	      fcc_data_type = 3
	   elseif (input .eq. 'FMS') then
	      pstatus = %loc(fip_normal)
              fcc_data_type = 0
	   else
	      pstatus = %loc(fip_invalfile)
	      call lib$signal (fip_invalfile, %val(1), fcc_infile(1:fcc_inlen))
	   endif
	   fcc_num_type = fac_present
	   fcc_diff = fac_not_present
	   fcc_infile = 'CSDR$FIRAS_IN:' // fcc_infile
	   call str$trim (fcc_infile, fcc_infile, fcc_inlen)
c
c  Get the input scan mode and filename extension.
c
	   period = index(fcc_infile, '.')
	   fcc_scan_mode = fcc_infile(period-4:period-1)
	   file_ext = fcc_infile(period+1:fcc_inlen)
	   call str$trim (file_ext, file_ext, extlen)
c
c  Set the output filename
c
	   bar = index(file_ext, '_')
	   fcc_outfile = 'CSDR$FIRAS_OUT:FIP_SKY_' // fcc_scan_mode // '.'
     .		          // input // file_ext(bar:extlen)
	   call str$trim (fcc_outfile, fcc_outfile, fcc_outlen)
           if (fcc_scan_mode(3:4) .eq. 'FA') then
              fcc_data_type = 2
           else if (fcc_data_type .eq. 0) then
              fcc_data_type = 1
           end if
	else
	   pstatus = %loc(fip_nofile)
	   call lib$signal (fip_nofile)
	endif





c
c  Set the galactic latitude exclusion invocation flags
c
	status = upm_present ('EXC_GLAT')
	if (status .eq. upm_pres) then
	   fcc_galexc = fac_present
	   status = upm_get_float ('EXC_GLAT', glat, len)
	   if (status .eq. ss$_normal) then
	      if ((glat .ge. 0.0)  .or.  (glat .le. fcc_glat_default)) then
	         fcc_glat = glat
	      else
	         pstatus = %loc(fip_galexc)
	         call lib$signal (fip_galexc, %val(1), %val(status))
	      endif
	   else
	      pstatus = %loc(fip_galexc)
	      call lib$signal (fip_galexc, %val(1), %val(status))
	   endif
	else
	   fcc_galexc = fac_not_present
	   fcc_glat = fcc_glat_default
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
	      bar = index(file_ext, '_')
	      rep_mid = file_ext(bar+1:extlen)
	      call str$trim (rep_mid, rep_mid, midlen)
	      fcc_report_file = 'FIP_SKY_' // input // '_' // fcc_scan_mode //
     .				'_' // rep_mid(1:midlen) // '.' //
     .				'REP_' // current_gmt(1:9)
	   else				     !  Report file name from invocation
	      call str$upcase (fcc_report_file, fcc_report_file)
	   endif
	   call str$trim (fcc_report_file, fcc_report_file, fcc_replen)
	endif


	fip_parse_sky = pstatus

	return
	end
