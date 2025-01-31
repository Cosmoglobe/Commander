	integer * 4 function  fip_init_ejv_report (current_gmt, parse_status)

c-------------------------------------------------------------------------------
c
c	Function FIP_INIT_EJV_REPORT
c
c	This function initializes the FIP_EJV processing report.  The report
c	is opened, then the error handler is invoked.  Finally, the account
c	information, command line invocation, and logical pointer translations
c	are written to the report.
c
c	Author:	 Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 8 July 1994
c
c-------------------------------------------------------------------------------
c
c	Input:
c		current_gmt		character * 14		invocation time
c
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
c		fip_invoc_ejv.txt
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
	include '(fip_invoc_ejv)'
	include '($jpidef)'

	character * 14	current_gmt		!  GMT time of invocation
	character * 72	log_name(3)		!  logical name to be translated
	character * 72	temp_log		!  logical name buffer
	character * 72	trans_name(3)		!  translated logical name
	character * 20	username		!  user name of invoking account

	integer * 2	len			!  string length
	integer * 2	temp_len		!  logical name buffer length
	integer * 2	trans_len(3)		!  translated logical name len

	integer * 4	io_stat			!  I/O return status
	integer * 4	j			!  a counter
	integer * 4	k			!  a counter
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
	external	fip_noejv
	external	fip_noreftime
	external	fip_normal
	external	fip_repopen
	external	fip_repwrite
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
	   status = %loc(fip_normal)
	else
	   status = %loc(fip_repopen)
	   call lib$signal (fip_repopen, %val(2), fcc_report_file(1:fcc_replen),
     .					 %val(io_stat))
	endif


	if (status .eq. %loc(fip_normal)) then
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
     .					'FIRAS Facility FIP_EJV')
	   write (fut_report_lun,10,iostat=io_stat)
	   write (fut_report_lun,20,iostat=io_stat) username
	   write (fut_report_lun,30,iostat=io_stat) current_gmt(1:11)
	   write (fut_report_lun,40,iostat=io_stat) fcc_scan_mode
	   if (fcc_flight .eq. fac_present) then
	      write (fut_report_lun,50,iostat=io_stat)
	   else
	      write (fut_report_lun,60,iostat=io_stat)
	   endif
	   write (fut_report_lun,70,iostat=io_stat) fcc_file_ext
	   write (fut_report_lun,80,iostat=io_stat)
     .						   fcc_report_file(1:fcc_replen)
  10	   format (31x, 'Processing Report', //)
  20	   format (x, 'Run by:                         ', a)
  30	   format (x, 'Run time:                       ', a)
  40	   format (x, 'Channel / Scan Mode:            ', a)
  50       format (x, 'Reference Dateset Time:         FLIGHT')
  60       format (x, 'Reference Dateset Time:         INT')
  70	   format (x, 'Offset / Gain File Extension:   ', a)
  80	   format (x, 'Processing Report File:         ', a)
	   if (io_stat .ne. 0) then
	      status = %loc(fip_repwrite)
	      call lib$signal (fip_repwrite, %val(2),
     .			       fcc_report_file(1:fcc_replen), %val(io_stat))
	   endif

	endif	!	(status from open


	if (status .eq. %loc(fip_normal)) then
C
C  Translate the software logical names.
C

c
c  Set up the logical names to be translated.
c
	   log_name(1) = 'CSDR$FIRAS_IN'
	   log_name(2) = 'CSDR$FIRAS_OUT'
	   log_name(3) = 'CSDR$FIRAS_REF'

c
c  Get the translation.
c
	   do j = 1,3
	      tstatus = cut_translate_archive_id (log_name(j), temp_log,
     .				      temp_len, trans_name(j), trans_len(j))
	      if (tstatus .ne. %loc(cut_normal)) call lib$signal (%val(tstatus))
	      call str$trim (trans_name(j), trans_name(j), trans_len(j))
	   enddo

c
c  Write out the translated logicals.
c
	   write (fut_report_lun,90,iostat=io_stat)
	   do j = 1,3
	      write (fut_report_lun,100,iostat=io_stat) log_name(j)(1:14),
     .						   trans_name(j)(1:trans_len(j))
	   enddo
  90	   format (//, x, 'Logical Name Translations:', /)
 100	   format (4x, a, 15x, a)
	   if (io_stat .ne. 0) then
	      status = %loc(fip_repwrite)
	      call lib$signal (fip_repwrite, %val(2),
     .			       fcc_report_file(1:fcc_replen), %val(io_stat))
	   endif

	endif	!	(status from open and writes


	if (status .eq. %loc(fip_normal)) then
C
C  Report processing errors.
C

c
c  If return status from the command line parse is bad:
c
	   if (parse_status .ne. %loc(fip_normal)) then
	      write (fut_report_lun,110,iostat=io_stat)
	      if (parse_status .eq. %loc(fip_noejv)) then
	         write (fut_report_lun,120,iostat=io_stat)
	      elseif (parse_status .eq. %loc(fip_noreftime)) then
	         write (fut_report_lun,130,iostat=io_stat)
	      endif
	   endif
 110	   format (//, x, 'Error from FIP_PARSE_EJV:', /)
 120	   format (4x, 'No ejv file specified for /FILE_EXT qualifier.')
 130	   format (4x, 'Neither /FLIGHT nor /INT reference time qualifier set.')
	   if (io_stat .ne. 0) then
	      status = %loc(fip_repwrite)
	      call lib$signal (fip_repwrite, %val(2),
     .			       fcc_report_file(1:fcc_replen), %val(io_stat))
	   endif

	endif	!	(status from open and writes


	fip_init_ejv_report = status

	return
	end
