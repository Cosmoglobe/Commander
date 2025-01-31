	Integer*4  Function  FCS_Report_Check  ( report, reportf, repdef,
	1                                        version, cmdline, clin, clen,
	2                                        input, chan, scan, sfile,
	3                                        in_skymap, tstart, tstop,
	4                                        fext, snum, arch_in, arch_out,
	5                                        lun_rpt, goodmap )

c Function to check the script file data, open, and initialize the processing
c report file.
c Author: Larry P. Rosen, March 1992, Hughes STX
c Modified: Larry P. Rosen, October 1993.  Removed two scan modes possibility.

	Implicit None

c Include

	Include		'($ssdef)'
	Include		'(fut_error)'
	Include		'(fut_params)'
	Include		'(fcs_msg)'
	Include		'($jpidef)'

c Passed parameters

	Logical*1	report				! flag to write report
	Character*36	reportf				! Filename for report
	Logical*1	repdef				! Report flag
	Character*6	version
	Character*79	cmdline(3)		! Command line with defaults
	Integer*2	clin			! which line of command line
	Integer*2	clen			! length of current command line
	Character*3	input			! Input data type; FAD or FCF
	Character*2	chan, scan			! Channel, scan mode
	Character*80	sfile				! Script file name
	Character*56	in_skymap (fac_max_num)		! Input skymaps
	Character*14	tstart (fac_max_num)		! Start times
	Character*14	tstop (fac_max_num)		! Stop times
	Character*20	fext				! Output file extension
	Integer*2	snum				! # of input skymaps
	Character*(*)	arch_in				! input data archive
	Character*(*)	arch_out			! output data archive
	Integer*4	lun_rpt				! logic unit for report
	Logical*1	goodmap (fac_max_num)	! map has right chan & scan

c Function

	Integer*4	CUT_Register_Version, CUT_Display_Banner
	Integer*4	Lib$Get_Lun
	Integer*4	Time_LT
	Integer*4	Lib$GetJPI
	Integer*4	Sys$AscTim
	Integer*4 	CUT_Translate_Archive_ID
	Integer*4	Lib$Rename_File

c External

	External	fut_error
	External  	cut_translate_archive_id

c Local

	Integer*4	rstat
	Character*20	tempfil /'FCS_REPORT.TEMPORARY'/
	Character*14	gmtmin, gmtmax
	Integer*4	tmin(2), tmax(2), adtstart(2), adtstop(2)
	Integer*2	i, j
	Character*8	owner			! invoking user name
	Integer*4	current_time(2)		! current system adt time
	Integer*2	time_len		! length of time string
	Character*32	time			! current system time string
	Character*72	logn, flogn, tlogn, blog ! logical name and translated
	Integer*4	flen, tlen		! length of translated logicals
	Integer*4	larc, larct		! length of archive names
	Integer*2	numsky /0/, numchan /0/, numscan /0/, numgood /0/
	Character*3	icheck			! check input data type
	Character*2	ccheck, scheck		! check channel, scan mode
	Logical*1	changood, scangood	! for checking chan & scan
	Logical*1	inpgood			! for checking input data type
c-----------------------------------------------------------------------------
c Begin

	FCS_Report_Check = %loc (fcs_normal)
	If (report) Then
	   rstat = Lib$Get_Lun (lun_rpt)
	   If (rstat .NE. ss$_normal) Then
	      FCS_Report_Check = %loc (fcs_abort)
	      Call Lib$Signal (fcs_lunerr, %val(1), %val(rstat))
	   Else
	      fut_report_lun = lun_rpt
	      Call Lib$Establish (fut_error)

c If repdef is true, report name was defaulted and there was no file extension
c given in the script file.  In this case, open a temporary report file
c 'FCS_REPORT.TEMPORARY' to catch any error messages.  Later close, rename, and
c reopen the report file for appending.

	      If (repdef) Then
	         Open (Unit=lun_rpt, File=tempfil, Status='new',
	1              Form='formatted', Access='sequential',
	2              Organization='sequential', Iostat=rstat )

	         If (rstat .NE. 0) Then
	            FCS_Report_Check = %loc (fcs_abort)
	            Call Lib$Signal (fcs_openrep, %val(2), tempfil,%val(rstat))
	         Endif
	      Else
	         Open (Unit=lun_rpt, File=reportf, Status='new',
	1              Form='formatted', Access='sequential',
	2              Organization='sequential', Iostat=rstat )
	         If (rstat .NE. 0) Then
	            FCS_Report_Check = %loc (fcs_abort)
	            Call Lib$Signal (fcs_openrep, %val(2), reportf,%val(rstat))
	         Endif
	      Endif
	   Endif
	   If (FCS_Report_Check .EQ. %loc (fcs_normal)) Then
	      rstat = CUT_Register_Version (version)
	      rstat  = CUT_Display_Banner (lun_rpt, 80,
	1        'FIRAS  Facility  FCS_Combine_Spectra')

c Write user and current time to the report file.

	      rstat = Lib$GetJPI (jpi$_username,,,,owner,)
	      Call Sys$GetTim (current_time)
	      rstat = Sys$AscTim (time_len, time, current_time, 0)
	      Write (lun_rpt,10) owner, time (1:time_len)
  10	      Format (' Run by:   ', A, '   at  Time: ',A,/)

c Write command line with defaults to the report file.

	      Do i = 1, clin
	         Write (lun_rpt,20) cmdline(i)
	      Enddo
  20	      Format (1X,A)
	      Write (lun_rpt,*)

c Write translation of logical names.

	      Call Str$Upcase (arch_in, arch_in)
	      Larc = Len (arch_in)
	      logn (1:larc) = arch_in
	      rstat = CUT_Translate_Archive_ID (logn, flogn, flen, tlogn, tlen)
	      Write (lun_rpt,30) 'Input', arch_in, tlogn (1:tlen)
  30	      Format (2X, 'Logical Translation for ', A6, ' Archive:', /, 5X,
	1             A, ' = ', A)
	      Call Str$Upcase (arch_out, arch_out)
	      larc = len (arch_out)
	      logn = blog
	      logn (1:larc) = arch_out
	      rstat = CUT_Translate_Archive_ID (logn, flogn, flen, tlogn, tlen)
	      Write (lun_rpt,30) 'Output', arch_out, tlogn (1:tlen)

c Write data type to report.

	      Write (lun_rpt,40) input
  40	      Format (/,2X,'Input data type: ',A3,/)
	   Endif
	Endif			! if report

c Whether or not a report is requested, check the skymap names and times in the 
c script file.  Also find the earliest and latest, and use them for the file
c extension if no extension is given in the script file.

	If (FCS_Report_Check .EQ. %loc (fcs_normal)) Then
	   gmtmin = fac_jstop_default
	   gmtmax = fac_jstart_default
	   Call CT_GMT_To_Binary (gmtmin, tmin)
	   Call CT_GMT_To_Binary (gmtmax, tmax)
	   Do i = 1, snum
	      goodmap(i) = .FALSE.
	      j = Index (in_skymap(i),' ')

c Check that file has correct channel and scan mode.

	      If (j .NE. 0) Then
	         numsky = numsky + 1
	         icheck = in_skymap(i)(1:3)
	         Call Str$Upcase (icheck, icheck)
	         ccheck = in_skymap(i)(9:10)
	         Call Str$Upcase (ccheck, ccheck)
	         scheck = in_skymap(i)(11:12)
	         Call Str$Upcase (scheck, scheck)
	         If (icheck .EQ. input) Then
	            inpgood = .TRUE.
	            numgood = numgood + 1
	         Else
	            inpgood = .FALSE.
	         Endif
	         If (ccheck .EQ. chan) Then
	            changood = .TRUE.
	            numchan = numchan + 1
	         Else
	            changood = .FALSE.
	         Endif
	         If (scheck .EQ. scan) Then
	            scangood = .TRUE.
	            numscan = numscan + 1
	         Else
	            scangood = .FALSE.
	         Endif
	         If (changood .AND. scangood .AND. inpgood) Then
	            goodmap(i) = .TRUE.
	            Call CT_GMT_To_Binary (tstart(i), adtstart)
	            Call CT_GMT_To_Binary (tstop(i), adtstop)
	            If ( Time_LT (adtstop, adtstart) ) Then
	               FCS_Report_Check = %loc (fcs_abort)
	               Call Lib$Signal (fcs_invtime, %val(1), in_skymap(i))
	            Endif
	            If ( Time_LT (adtstart, tmin) ) Then
	               tmin(1) = adtstart(1)
	               tmin(2) = adtstart(2)
	            Endif	         
	            If ( Time_LT (tmax, adtstop) ) Then
	               tmax(1) = adtstop(1)
	               tmax(2) = adtstop(2)
	            Endif
	         Endif
	         If (.NOT. inpgood) Then
	            Call Lib$Signal (fcs_wrong_input, %val(1), in_skymap(i) )
	         Endif
	         If (.NOT. changood) Then
	            Call Lib$Signal (fcs_wrong_chan, %val(1), in_skymap(i) )
	         Endif
	         If (.NOT. scangood) Then
	            Call Lib$Signal (fcs_wrong_scan, %val(1), in_skymap(i) )
	         Endif
	      Endif
	   Enddo
	   If (numchan .EQ. 0) Then
	      FCS_Report_Check = %loc (fcs_abort)
	      Call Lib$Signal (fcs_zchannel)
	   Endif
	   If (numscan .EQ. 0) Then
	      FCS_Report_Check = %loc (fcs_abort)
	      Call Lib$Signal (fcs_zscan)
	   Endif
	   If (numgood .EQ. 0) Then
	      FCS_Report_Check = %loc (fcs_abort)
	      Call Lib$Signal (fcs_zdatatype, %val(1), input)
	   Endif
	   If (numsky .EQ. 0) Then
	      FCS_Report_Check = %loc (fcs_abort)
	      Call Lib$Signal (fcs_noinput)
	   Endif
	Endif
	If (FCS_Report_Check .EQ. %loc(fcs_normal) .AND. fext(1:1) .EQ. ' ')Then
	   Call CT_Binary_To_GMT (tmin, gmtmin)
	   Call CT_Binary_To_GMT (tmax, gmtmax)
	   fext = input // '_' // gmtmin(1:7) // '_' //
	1         gmtmax(1:7)
	   If (repdef) Then
	      reportf = reportf(1:16) // fext(1:19)
	      Close (lun_rpt, Iostat=rstat)
	      If (rstat .NE. 0) Then
	         FCS_Report_Check = %loc (fcs_abort)
	         Call Lib$Signal (fcs_closerep, %val(2), reportf, %val(rstat))
	      Else
	         rstat = Lib$Rename_File (tempfil, reportf)
	         If (rstat .NE. ss$_normal) Then
	            FCS_Report_Check = %loc (fcs_abort)
	            Call Lib$Signal (fcs_renerr, %val(1), %val(rstat))
	         Else
	            Open (Unit=lun_rpt, File=reportf, Status='old',
	1                 Form='formatted', Access='append',
	2                 Organization='sequential', Iostat=rstat)

	            If (rstat .NE. 0) Then
	               FCS_Report_Check = %loc (fcs_abort)
	               Call Lib$Signal (fcs_openrep,%val(2),reportf,%val(rstat))
	            Endif
	         Endif
	      Endif
	   Endif
	Endif
	If (FCS_Report_Check .EQ. %loc (fcs_normal) .AND. report) Then
	   Write (lun_rpt,50) 'Report', reportf
   50	   Format (2X, A6, 1X, 'file name is : ',A)
	   j = Index (sfile,' ')
	   Write (lun_rpt,50) 'Script', sfile(1:j-1)

c Check script file name for a logical.  If there might be one (indicated by
c presence of ':', then do a cut_translate.  If the translated logical, tlogn,
c is the same as the script logical, logn, then it wasn't a logical, but was a
c directory and so nothing is written.

	   j = Index (sfile, ':')
	   If (j .GE. 1) Then
	      logn = blog
	      logn (1:j-1) = sfile (1:j-1)
	      rstat = CUT_Translate_Archive_ID (logn, flogn, flen, tlogn, tlen)
	      If (tlogn(1:j-1) .NE. logn(1:j-1)) Then
	         Write (lun_rpt,30) 'Script', sfile(1:j-1), tlogn(1:tlen)
	      Endif
	   Endif
	   Write (lun_rpt,*)
	1     ' FCS interprets the script file information as:'
	   Write (lun_rpt,*) '    File extension = ',fext
	   Do i=1,snum
	      Write (lun_rpt,60) i, in_skymap(i), tstart(i), tstop(i)
   60	      Format (5X, 'Skymap ',I3,2X,':',2X,A,/,20X,'jstart = ',A14,5X,
	1             'jstop = ',A14)
	   Enddo
	   Write (lun_rpt,*)' '
	Endif
	Return
	End
