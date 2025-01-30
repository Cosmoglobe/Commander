	integer * 4 function  frd_parse_variances (current_time, current_gmt)

c-------------------------------------------------------------------------------
c
c	Function FRD_PARSE_VARIANCES
c
c
c	This function parses the command line for FRD_VARIANCES.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  27 October 1993
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		current_time	integer * 4 array	ADT time of invocation
c		current_gmt	character * 14		GMT time of invocation
c
c	Subroutines called:
c		ct_binary_to_gmt
c		fut_free_lun
c		fut_get_lun
c		lib$signal
c		str$trim
c		str$upcase
c		sys$gettim
c		upm_get_float
c		upm_get_value
c		upm_present
c
c	Include files:
c		$ssdef
c		frd_invoc_variances.txt
c		fut_params.txt
c		upm_stat_msg.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '($ssdef)'
	include '(fut_params)'
	include '(frd_invoc_variances)'
	include '(upm_stat_msg)'

	character * 14	current_gmt		!  GMT time of invocation

	integer * 4	current_time(2)		!  ADT time of invocation
	integer * 4	io_stat			!  I/O return status
	integer * 4	j			!  a counter
	integer * 4	pstatus			!  parse return status
	integer * 4	rstatus			!  return status
	integer * 4	slun			!  script file lun
	integer * 4	status			!  return status

	real	* 4	glat			!  cutoff latitude

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun
	integer * 4	upm_get_float
	integer * 4	upm_get_value
	integer * 4	upm_present

c
c  Set up the namelist for the script file.
c
	namelist / smaps / skymap, filext
	character*20	skymap (fcc_max_skymaps)	! Input skymaps
	character*20	filext				! output filename ext

	external	frd_galexc
	external	frd_normal
	external	frd_nofilext
	external	frd_noscript
	external	frd_rmsclose
	external	frd_rmsopen
	external	frd_rmsread
	external	fut_normal

C
C  Parse the command line.
C

c
c  Get the time of the invocation
c
	call sys$gettim (current_time)
	call ct_binary_to_gmt (current_time, current_gmt)

c
c  Set the namelist script file invocation flag.
c
	status = upm_present ('SCRIPT')
	if (status .eq. upm_pres) then
	   status = upm_get_value ('SCRIPT', fcc_sfile, fcc_sflen)
	   call str$upcase (fcc_sfile, fcc_sfile)
	   pstatus = %loc(frd_normal)
	else
	   pstatus = %loc(frd_noscript)
	   call lib$signal (frd_noscript)
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
	fcc_scan_mode = fac_channel_ids(fcc_chan) //
     .			fac_scan_mode_ids(fcc_smode)

c
c  Set the galactic latitude exclusion invocation flags
c
	status = upm_present ('EXC_GLAT')
	if (status .eq. upm_pres) then
	   fcc_galexc = fac_present
	   status = upm_get_float ('EXC_GLAT', glat)
	   if (status .eq. ss$_normal) then
	      if ((glat .ge. fcc_glat_default)  .or.  (glat .le. 90.0)) then
	         fcc_glat = glat
	      else
	         pstatus = %loc(frd_galexc)
	         call lib$signal (frd_galexc, %val(1), %val(pstatus))
	      endif
	   else
	      pstatus = %loc(frd_galexc)
	      call lib$signal (frd_galexc, %val(1), %val(status))
	   endif
	else
	   fcc_galexc = fac_not_present
	   fcc_glat = fcc_glat_default
	endif


	if (pstatus .eq. %loc(frd_normal)) then
C
C  Read the namelist file.
C

c
c  Initialize the input arrays.
c
	   do j=1,fcc_max_skymaps
	      skymap(j) = ' '
	      filext    = ' '
	   enddo
	   fcc_mapnum = 0

c
c  Open the file.
c
	   rstatus = fut_get_lun(slun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	   open (unit=slun, file=fcc_sfile, access='sequential', status='old',
     .	         readonly, iostat=io_stat)

	   if (io_stat .eq. 0) then
c
c  Read the skymap filenames from the file.
c
	      read (slun, nml=smaps, iostat=io_stat)
	      if (io_stat .eq. 0) then
	         j = 1
	         do while (j .le. fcc_max_skymaps)
	            if (skymap(j)(1:1) .ne. ' ') then
	               fcc_mapnum = fcc_mapnum + 1
	               fcc_infile(j) = 'CSDR$FIRAS_IN:FCF_SKY_' // fcc_scan_mode
     .				       // '.' // skymap(j)
	               call str$upcase (fcc_infile(j), fcc_infile(j))
	               call str$trim (fcc_infile(j), fcc_infile(j),
     .						     fcc_inlen(j))
	            endif
	            j = j + 1
	         enddo
	         if (filext(1:1) .eq. ' ') then
	            pstatus = %loc(frd_nofilext)
	            call lib$signal (frd_nofilext, %val(1),
     .				     fcc_sfile(1:fcc_sflen))
	         endif
c
c  Determine the output file name.
c
	         fcc_outfile = 'CSDR$FIRAS_OUT:FEX_VAR_' // fcc_scan_mode //
     .			       '.' // filext
	         call str$upcase (fcc_outfile, fcc_outfile)
	         call str$trim (fcc_outfile, fcc_outfile, fcc_outlen)

c
c  Close the namelist file.
c
	         close (slun, iostat=io_stat)
	         if (io_stat .ne. 0) then
	            pstatus = %loc(frd_rmsread)
	            call lib$signal (frd_rmsclose, %val(2),
     .				     fcc_sfile(1:fcc_sflen), %val(io_stat))
	         endif

	      else
	         pstatus = %loc(frd_rmsread)
	         call lib$signal (frd_rmsread, %val(2), fcc_sfile(1:fcc_sflen),
     .					       %val(io_stat))
	      endif

	   else
              pstatus = %loc(frd_rmsopen)
	      call lib$signal (frd_rmsopen, %val(2), fcc_sfile(1:fcc_sflen),
     .					    %val(io_stat))
	   endif

	endif		!  pstatus from parse


	frd_parse_variances = pstatus

	return
	end
