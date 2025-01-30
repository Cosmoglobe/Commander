	Integer*4  Function  FIP_INIT_REPORT_COV  ( filext, chan, scan,
	1                                           freq_low, freq_high,
	1                                           repfile, cmdline,
	1                                           curtime, arch_in, arch_out,
	1                                           arch_ref, version, lun_rpt)

C  Open and initialize the report file; writing banner, user, translated
C  logicals, and fully defaulted command line.
C
C  Larry P. Rosen, Hughes STX, 19 July 1993.
C
C    Input:
C	Character*20  filext               ! Input file name extension
C	Character*2   chan                 ! Channel to process (RH,RL,LH,LL)
C	Character*2   scan                 ! Scan mode to process (SS,SF,LS,LF)
C	Real*4        freq_low             ! Low frequency cut-off (icm)
C	Real*4        freq_high            ! High frequency cut-off (icm)
C	Character*42  repfile              ! File name for report
C	Character*79  cmdline (3)          ! Command line with defaults
C	Integer*4     curtime(2)           ! System time in ADT
C	Character*(*) arch_in              ! input data archive
C	Character*(*) arch_out             ! output data archive
C	Character*(*) arch_ref             ! reference data archive
C	Character*6   version              ! Version of software.
C   Output:
C	Integer*4     lun_rpt              ! Logical unit number for report
C

	Implicit None

C  Include Files

	Include		'($ssdef)'
	Include		'(fut_error)'
	Include		'($jpidef)'

C Passed Parameters

	Character*20	filext             ! Input file name extension
	Character*2	chan               ! Channel to process (RH,RL,LH,LL)
	Character*2	scan               ! Scan mode to process (SS,SF,LS,LF)
	Real*4		freq_low           ! Low frequency cut-off (icm)
	Real*4		freq_high          ! High frequency cut-off (icm)
	Logical*1	report             ! Flag whether to write report
	Character*42	repfile            ! File name for report
	Character*79	cmdline (3)        ! Command line with defaults
	Integer*4	curtime (2)        ! System time in ADT
	Character*(*)	arch_in            ! input data archive
	Character*(*)	arch_out           ! output data archive
	Character*(*)	arch_ref           ! reference data archive
	Character*6	version            ! Version of software.
	Integer*4	lun_rpt            ! Logical unit number for report

C Function

	Integer*4	CUT_Register_Version, CUT_Display_Banner
	Integer*4	Lib$Get_Lun
	Integer*4	Time_LT
	Integer*4	Lib$GetJPI
	Integer*4	Sys$AscTim
	Integer*4 	CUT_Translate_Archive_ID
	Integer*4	Lib$Rename_File

C External

	External	fut_error
	External  	cut_translate_archive_id
	External	fip_normal
	External	fip_abort
	External	fip_lunerr
	External	fip_openerr

C Local

	Integer*4	rstat
	Character*8	owner			! invoking user name
	Integer*2	time_len		! length of time string
	Character*32	time			! current system time string
	Integer*2	i			! counter
	Character*72	logn, flogn, tlogn, blog ! logical name and translated
	Integer*4	flen, tlen		! length of translated logicals
	Integer*4	larc, larct		! length of archive names

C Begin

	FIP_INIT_REPORT_COV = %loc (fip_normal)

C Get logical unit number for report.

	rstat = Lib$Get_Lun (lun_rpt)
	If (rstat .NE. ss$_normal) Then
	   FIP_INIT_REPORT_COV = %loc (fip_abort)
	   Call Lib$Signal (fip_lunerr, %val(1), %val(rstat))
	Else

C Set up error handler.

	   fut_report_lun = lun_rpt
	   Call Lib$Establish (fut_error)

C Open the report file

	   Open ( Unit=lun_rpt, File=repfile, Status='new', Form='formatted',
	1         Access='sequential', Iostat=rstat)

	   If (rstat .NE. 0) Then
	      FIP_INIT_REPORT_COV = %loc (fip_abort)
	      Call Lib$Signal (fip_openerr, %val(2), repfile, %val(rstat))
	   Endif
	Endif

	If (FIP_INIT_REPORT_COV .EQ. %loc (fip_normal)) Then
	   rstat = CUT_Register_Version (version)
	   rstat  = CUT_Display_Banner ( lun_rpt, 80,
	1                                'FIRAS  Facility  FIP_COVAR')
	   Write (lun_rpt, 10) 'Processing Report'
  10	   Format (30X, A, /)

C Write command line with defaults to the report file.

	   Do i = 1, 3
	      Write (lun_rpt,*) cmdline (i)
	   Enddo
	   Write (lun_rpt,*)

C Write user and current time to the report file.

	   rstat = Lib$GetJPI (jpi$_username,,,,owner,)
	   rstat = Sys$AscTim (time_len, time, curtime, 0)
	   Write (lun_rpt, 20) owner, time (1:time_len)
  20	   Format (' Run by:', 21X, A, /, ' At time:', 20X, A)

C Write parsed info: chan, scan, filext, report name.

	   Write (lun_rpt, 25) chan, scan, filext, repfile
  25	   Format (' Channel / Scan Mode:', 8X, A2, A2, /,
	1          ' Input File Extension:', 7X, A, /,
	1          ' Processing Report File:', 5X, A, / )

c Write translation of logical names.

	   larc = Len (arch_in)
	   logn (1:larc) = arch_in
	   rstat = CUT_Translate_Archive_ID ( logn, flogn, flen, tlogn, tlen )
	   Write (lun_rpt, *) ' Logical Name Translations:'
	   Write (lun_rpt, 30) arch_in, tlogn (1:tlen)
  30	   Format ( 3X, A14, 9X, A )
	   larc = Len (arch_out)
	   logn = blog
	   logn (1:larc) = arch_out
	   rstat = CUT_Translate_Archive_ID (logn, flogn, flen, tlogn, tlen)
	   Write (lun_rpt, 35) arch_out, tlogn (1:tlen)
  35	   Format ( 3X, A15, 8X, A )
	   larc = Len (arch_ref)
	   logn = blog
	   logn (1:larc) = arch_ref
	   rstat = CUT_Translate_Archive_ID (logn, flogn, flen, tlogn, tlen)
	   Write (lun_rpt, 35) arch_ref, tlogn (1:tlen)
	   Write (lun_rpt, *)
	   Write (lun_rpt, 40) 'Input Data File: ', arch_in, 'FCC_COV_', chan,
	1                      scan, '.',  filext
  40	   Format (1X, A17, 3X, A, A, A2, A2, A1, A)
	   Write (lun_rpt, 40) 'Output Data File:', arch_out, 'FIP_COV_', chan,
	1                      scan, '.',  filext
	Endif
	Return
	End
