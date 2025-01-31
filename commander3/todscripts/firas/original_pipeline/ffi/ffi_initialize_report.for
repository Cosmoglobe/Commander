	integer * 4 function  ffi_initialize_report (current_gmt, ncmd,
     .						 cmd_line, cmdlen, parse_status)

c-------------------------------------------------------------------------------
c
c	Function FFI_INITIALIZE_REPORT
c
c	This function initializes the FFI processing report.  The report is
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
c		parse_status		integer * 4		status from
c								command line
c								parse
c
c	Output:
c		none
c
c	Subroutines called:
c		cut_display_banner
c		cut_register_version
c		cut_translate_archive_id
c		fut_get_lun
c		lib$establish
c		lib$getjpi
c		lib$signal
c		str$trim
c
c	Include files:
c		ffi_invoc.txt
c		fut_error.txt
c		fut_params.txt
c		$jpidef
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_error)'
	include '(fut_params)'
	include '(ffi_invoc)'
	include '($jpidef)'

	character * 79	cmd_line(3)		!  command line invocation
	character * 14	current_gmt		!  GMT time of invocation
	character * 72	log_name(4)		!  logical name to be translated
	character * 72	temp_log		!  logical name buffer
	character * 72	trans_name(4)		!  translated logical name
	character * 20	username		!  user name of invoking account

	integer * 2	len			!  string length
	integer * 2	temp_len		!  logical name buffer length
	integer * 2	trans_len(4)		!  translated logical name len

	integer * 4	cmdlen(3)		!  length of command lines in
						!    invocation
	integer * 4	io_stat			!  I/O return status
	integer * 4	j			!  a counter
	integer * 4	k			!  a counter
	integer * 4	ncmd			!  number of command lines in
						!    invocation
	integer * 4	n_names			!  number of logical names
						!    to be translated
	integer * 4	parse_status		!  command line parse return
						!    status
	integer * 4	rstatus			!  return status
	integer * 4	status			!  return status
	integer * 4	tstatus			!  return status


        integer	* 4	cut_display_banner
        integer	* 4	cut_register_version
	integer * 4	cut_translate_archive_id
	integer	* 4	fut_get_lun

	external	fut_error

	external	cut_normal
	external	ffi_invaltime
	external	ffi_nooutfile
	external	ffi_normal
	external	ffi_repopen
	external	ffi_repwrite
	external	fut_normal

C
C  Open the report file.
C
	rstatus = fut_get_lun(fut_report_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=fut_report_lun, file=fcc_report_file, status='new',
     .            iostat=io_stat)

	if (io_stat .eq. 0) then
	   status = %loc(ffi_normal)
	else
	   status = %loc(ffi_repopen)
	   call lib$signal (ffi_repopen, %val(2), fcc_report_file(1:fcc_replen),
     .					 %val(io_stat))
	endif


	if (status .eq. %loc(ffi_normal)) then
C
C  Initialize the report.
C

c
c  Initialize the error handler.
c
	   call lib$establish(fut_error)

c
c  Get the the username of the account invoking the program.
c
	   call lib$getjpi (jpi$_username,,,,username,)

c
c  Write the banner and the account and invocation information.
c
	   tstatus = cut_register_version (fcc_version)
	   tstatus = cut_display_banner (fut_report_lun, 80,
     .					'FIRAS Facility FFI_Fish_Input')
	   write (fut_report_lun,10,iostat=io_stat)
	   write (fut_report_lun,20,iostat=io_stat) username
	   write (fut_report_lun,30,iostat=io_stat) current_gmt(1:11)
	   write (fut_report_lun,40,iostat=io_stat) fcc_scan_mode
	   write (fut_report_lun,50,iostat=io_stat) fcc_jstart_time(1:11),
     .						    fcc_jstop_time(1:11)
	   write (fut_report_lun,60,iostat=io_stat)
     .						   fcc_report_file(1:fcc_replen)
  10	   format (31x, 'Processing Report', //)
  20	   format (x, 'Run by:                     ', a)
  30	   format (x, 'Run time:                   ', a)
  40	   format (x, 'Channel / Scan Mode:        ', a)
  50	   format (x, 'Timerange:                  ', a, ' - ', a)
  60	   format (x, 'Processing Report File:     ', a)
	   if (io_stat .ne. 0) then
	      status = %loc(ffi_repwrite)
	      call lib$signal (ffi_repwrite, %val(2),
     .			       fcc_report_file(1:fcc_replen), %val(io_stat))
	   endif

	endif	!	(status from open


	if (status .eq. %loc(ffi_normal)) then
C
C  Translate the software logical names.
C

c
c  Set up the logical names to be translated.
c
	   if (fcc_hybrid .eq. fac_present) then
	      n_names = 4
	      log_name(1) = 'CSDR$FIRAS_IN1'
	      log_name(2) = 'CSDR$FIRAS_IN2'
	      log_name(3) = 'CSDR$FIRAS_OUT'
	      log_name(4) = 'CSDR$FIRAS_REF'
	   else
	      n_names = 3
	      log_name(1) = 'CSDR$FIRAS_IN1'
	      log_name(2) = 'CSDR$FIRAS_OUT'
	      log_name(3) = 'CSDR$FIRAS_REF'
	   endif

c
c  Get the translation.
c
	   do j = 1,n_names
	      tstatus = cut_translate_archive_id (log_name(j), temp_log,
     .				      temp_len, trans_name(j), trans_len(j))
	      if (tstatus .ne. %loc(cut_normal)) call lib$signal (%val(tstatus))
	      call str$trim (trans_name(j), trans_name(j), trans_len(j))
	   enddo

c
c  Write out the translated logicals.
c
	   write (fut_report_lun,70,iostat=io_stat)
	   do j = 1,n_names
	      write (fut_report_lun,80,iostat=io_stat) log_name(j)(1:14),
     .						   trans_name(j)(1:trans_len(j))
	   enddo
  70	   format (//, x, 'Logical Name Translations:', /)
  80	   format (4x, a, 11x, a)
	   if (io_stat .ne. 0) then
	      status = %loc(ffi_repwrite)
	      call lib$signal (ffi_repwrite, %val(2),
     .			       fcc_report_file(1:fcc_replen), %val(io_stat))
	   endif

	endif	!	(status from open and writes


	if (status .eq. %loc(ffi_normal)) then
C
C  Write out the command line invocation
C

	   write (fut_report_lun,90,iostat=io_stat)
	   do j = 1,ncmd
	      write (fut_report_lun, 100,iostat=io_stat) cmd_line(j)
	   enddo
  90	   format (//, x, 'Command Line Invocation:', /)
 100	   format (x, a)
	   if (io_stat .ne. 0) then
	      status = %loc(ffi_repwrite)
	      call lib$signal (ffi_repwrite, %val(2),
     .			       fcc_report_file(1:fcc_replen), %val(io_stat))
	   endif

	endif	!	(status from open and writes


	if (status .eq. %loc(ffi_normal)) then
C
C  Report processing errors.
C

c
c  If return status from the command line parse is bad:
c
	   if (parse_status .ne. %loc(ffi_normal)) then
	      write (fut_report_lun,110,iostat=io_stat)
	      if (parse_status .eq. %loc(ffi_invaltime)) then
	         write (fut_report_lun,120,iostat=io_stat) fcc_jstop_time,
     .							   fcc_jstart_time
	      elseif (parse_status .eq. %loc(ffi_nooutfile)) then
	         write (fut_report_lun,130,iostat=io_stat)
	      endif
	   endif
 110	   format (//, x, 'Error from FFI_PARSE:', /)
 120	   format (4x, 'Jstop time ', a, ' is less than Jstart time ', a, '.')
 130	   format (4x, 'No output file name extension specified.')
	   if (io_stat .ne. 0) then
	      status = %loc(ffi_repwrite)
	      call lib$signal (ffi_repwrite, %val(2),
     .			       fcc_report_file(1:fcc_replen), %val(io_stat))
	   endif

	endif	!	(status from open and writes


	ffi_initialize_report = status

	return
	end
