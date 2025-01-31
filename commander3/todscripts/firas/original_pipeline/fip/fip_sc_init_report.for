	Integer*4  Function  FIP_SC_INIT_REPORT ( parse_status, current_gmt,
	1                                         input, sfile, slen,
	2                                         cmdline, filexts, fnum,
	3                                         archin, archout )

c------------------------------------------------------------------------------
c
c	Function FIP_SC_INIT_REPORT
c
c	This function initializes the FIP_Spectra_Coadd processing report.  The
c	report is opened, then the error handler is invoked.  Finally, the
c	account information, command line invocation, and logical pointer
c	translations are written to the report.
c
c	Author:	 Larry P. Rosen, Hughes STX, 13 July 1994, November 1994
c       Copied from FIP_INIT_SKY_REPORT by Gene Eplee with modifications
c------------------------------------------------------------------------------
c
c	Input:
c	integer*4	parse_status		! Command line parse status
c	character*14	current_gmt		! GMT time of invocation
c	character*12	input			! Input filename base
c	character*80	sfile			! Script file name
c	integer*2	slen			! Length of script file name
c	character*79	cmdline (3)		! Command line with defaults
c	character*20	filexts (fac_max_num)	! Input file extensions
c	integer*2	fnum			! Number of input files
c	character*13	archin			! Input archive name
c	character*14	archout			! Output archive name
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
c		fip_invoc_sky.txt
c		fut_error.txt
c		fut_params.txt
c		$jpidef
c
c------------------------------------------------------------------------------

	Implicit None

	Include '(fut_error)'
	Include '(fut_params)'
	Include '(FIP_invoc_sky)'
	Include '($jpidef)'

c Passed parameters:
	Integer*4	parse_status	!  command line parse return status
	Character*14	current_gmt		! GMT time of invocation
	Character*12	input			! Input filename base
	Character*80	sfile			! Script file name
	Integer*2	slen			! Length of script file name
	Character*79	cmdline (3)		! Command line with defaults
	Character*20	filexts (fac_max_num)	! Input file extensions
	Integer*2	fnum			! Number of input files
	Character*13	archin			! Input archive
	Character*14	archout			! Output archive

c Local:
	Integer*4	status			!  return status
	Integer*4	rstatus			!  return status
	Integer*4	io_stat			!  I/O return status
	Character*20	username		!  user name invoking account
	Integer*4	tstatus			!  return status
	Integer*4	j			!  counter
	Character*72	log_name(2)	!  logical name to be translated
	Character*72	temp_log		!  logical name buffer
	Character*72	trans_name(2)		!  translated logical name
	Integer*2	len			!  string length
	Integer*2	temp_len		!  logical name buffer length
	Integer*2	trans_len(2)		!  translated logical name len

c Functions:
	Integer*4	FUT_get_lun
        Integer*4	CUT_display_banner
        Integer*4	CUT_register_version
	Integer*4	CUT_translate_archive_id

C Externals:
	External	FIP_normal
	External	FUT_normal
	External	FIP_lunerr
	External	FUT_error
	External	FIP_repopen
	External	FIP_repWrite
	External	CUT_normal
	External	FIP_parserr
	External	FIP_galexc
	External	FIP_openerr
	External	FIP_readerr
	External	FIP_invalfile
	External	FIP_nofile
	External	FIP_nofilext
	External	FIP_noinput
	External	FIP_closerr

C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C Begin
	status = %loc (FIP_Normal)
C
C  Open the report file.
C
	rstatus = FUT_Get_Lun (FUT_Report_Lun)
	If (rstatus .NE. %loc(FUT_Normal)) Then
	   Call LIB$Signal (FIP_lunerr, %val(1), %val(rstatus))
	Else
	   Open ( unit=FUT_Report_Lun, file=fcc_report_file, status='new',
	1         iostat=io_stat )

	   If (io_stat .EQ. 0) Then
	      status = %loc(FIP_Normal)
	   Else
	      status = %loc(FIP_RepOpen)
	      call Lib$Signal ( FIP_RepOpen, %val(2),
	1                       fcc_report_file(1:fcc_replen), %val(io_stat) )
	   EndIf
	EndIf

	If (status .EQ. %loc(FIP_Normal)) Then
C
C  Initialize the report.
C
c
c  Initialize the error handler.
c
	   Call Lib$establish (FUT_Error)
c
c  Get the the username of the account invoking the program.
c
	   Call Lib$getjpi (jpi$_username,,,,username,)
c
c  Write the banner and the account and invocation information.
c
	   tstatus = CUT_register_version (fcc_version)
	   tstatus = CUT_display_banner ( FUT_Report_Lun, 80,
	1                                 'FIRAS Facility FIP_SPECTRA_COADD' )
	   Write (FUT_Report_Lun,10,iostat=io_stat)
	   Write (FUT_Report_Lun,20,iostat=io_stat) username
	   Write (FUT_Report_Lun,30,iostat=io_stat) current_gmt(1:11)
	   Write (FUT_Report_Lun,40,iostat=io_stat) fcc_scan_mode
	   Write (FUT_Report_Lun, 50, iostat=io_stat)
	1         fcc_report_file (1:fcc_replen)
  10	   Format (31x, 'Processing Report', //)
  20	   Format (x, 'Run by:                     ', a)
  30	   Format (x, 'Run time:                   ', a)
  40	   Format (x, 'Channel / Scan Mode:        ', a)
  50	   Format (x, 'Processing Report File:     ', a, // )
	   If (io_stat .NE. 0) Then
	      status = %loc (FIP_RepWrite)
	      Call Lib$Signal ( FIP_RepWrite, %val(2),
	1                       fcc_report_file(1:fcc_replen), %val(io_stat) )
	      FIP_SC_Init_Report = status
	      Return
	   EndIf

c Write out complete command line with defaults

	   Do j = 1, 3
	      Write (FUT_Report_Lun, 60) cmdline (j)
	   EndDo
  60	   Format (1x,a)
	EndIf	!	(status from open

	If (status .EQ. %loc(FIP_Normal)) Then
C
C  Translate the software logical names.
C
c
c  Set up the logical names to be translated.
c
	   log_name(1) = archin
	   log_name(2) = archout
c
c  Get the translation.
c
	   Do j = 1,2
	      tstatus = CUT_Translate_Archive_ID ( log_name(j), temp_log,
	1                 temp_len, trans_name(j), trans_len(j) )
	      If (tstatus .NE. %loc(CUT_Normal)) Call Lib$Signal(%val(tstatus))
	      Call str$trim (trans_name(j), trans_name(j), trans_len(j))
	   EndDo
c
c  Write out the translated logicals.
c
	   Write (FUT_Report_Lun, 70, iostat=io_stat)
	   Do j = 1,2
	      Write (FUT_Report_Lun,80,iostat=io_stat) log_name(j)(1:14),
	1        trans_name(j)(1:trans_len(j))
	   EndDo
  70	   Format (/, x, 'Logical Name Translations:')
  80	   Format (4x, a, 11x, a)
	   If (io_stat .ne. 0) Then
	      status = %loc(FIP_RepWrite)
	      call Lib$Signal (FIP_RepWrite, %val(2),
	1                      fcc_report_file(1:fcc_replen), %val(io_stat))
	   EndIf
	EndIf	!	(status from open and Writes

	If (status .EQ. %loc(FIP_normal)) Then
	
C  Report file names

	   Write (FUT_Report_Lun,90) sfile (1:slen)
 90	   Format (/,1x,'Script file of input files to reformat:',/,4x,a)
	   Write (FUT_Report_Lun,100)
 100	   Format (/,1x,'Input filenames to process:')
	   Do j = 1, fnum
	      Write (FUT_Report_Lun,110) input, filexts(j)
	   EndDo
 110	   Format (4x,a,'.',a)
C
C  Report processing errors.
C
c
c  If return status from the command line parse is bad:
c
	   If (parse_status .NE. %loc(FIP_normal)) Then
	      If (parse_status .EQ. %loc (FIP_parserr)) Then
	         call Lib$Signal (FIP_parserr, %val(1), '<see log file>')
	      ElseIf (parse_status .EQ. %loc (FIP_galexc)) Then
	         call Lib$Signal (FIP_galexc, %val(1), 0)
	      ElseIf (parse_status .EQ. %loc (FIP_lunerr)) Then
	         call Lib$Signal (FIP_lunerr, %val(1), 0)
	      ElseIf (parse_status .EQ. %loc (FIP_openerr)) Then
	         call Lib$Signal (FIP_openerr, %val(2), sfile(1:slen), 0)
	      ElseIf (parse_status .EQ. %loc (FIP_readerr)) Then
	         call Lib$Signal (FIP_readerr, %val(2), sfile(1:slen), 0)
	      ElseIf (parse_status .EQ. %loc (FIP_closerr)) Then
	         call Lib$Signal (FIP_closerr, %val(2), sfile(1:slen), 0)
	      ElseIf (parse_status .EQ. %loc (FIP_nofilext)) Then
	         call Lib$Signal (FIP_nofilext)
	      ElseIf (parse_status .EQ. %loc (FIP_noinput)) Then
	         call Lib$Signal (FIP_noinput)
	      Else
	         Write (FUT_Report_Lun,*) ' Unknown error during parse: ',
	1           parse_status
	      EndIf
	   EndIf
	EndIf	!	(status from open and Writes

	FIP_SC_Init_Report = status

	Return
	End
