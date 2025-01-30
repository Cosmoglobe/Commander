	Integer*4  Function  FCS_Parse  ( sfile, input, chan, scan, report,
	1                                 reportf, cmdline, in_skymap, tstart,
	2                                 tstop, fext, snum, repdef, clin,
	3                                 clen ) 

c------------------------------------------------------------------------------
c Function to parse the fcs command line and read the script file
c Author: Larry P. Rosen, March 1992, Hughes STX
c Modifications: Larry P. Rosen, 6/11/93, HSTX, SPR 11057.  Script file must
c                be opened with the readonly qualifier.
c                Larry P. Rosen, October 1993.  No more two scan mode option.
c------------------------------------------------------------------------------

	Implicit None

c Include

	Include		'($ssdef)'
	Include		'(upm_stat_msg)'
	Include		'(fut_params)'
	Include		'(fcs_msg)'
	Include		'(cct_query_catalog_record)'

c Passed parameters

	Character*80	sfile			! Script file name
	Character*3	input			! Input data type FAD or FCF
	Character*2	chan, scan		! Channel, scan mode
	Logical*1	report			! flag whether to write report
	Character*36	reportf			! Filename for report
	Character*79	cmdline (3)		! Command line with defaults
	Character*56	in_skymap (fac_max_num)		! Input skymaps
	Character*14	tstart (fac_max_num)		! Start times
	Character*14	tstop (fac_max_num)		! Stop times
	Character*20	fext				! Output file extension
	Integer*2	snum				! # of input skymaps
	Logical*1	repdef			! True = use default report name
	Integer*2	clin			! which line of command line
	Integer*2	clen			! length of current command line

c Namelist

	Namelist / smaps / skymap, jstart, jstop, filext
	Character*56	skymap (fac_max_num)			! Input skymaps
	Character*14	jstart (fac_max_num)			! Start times
	Character*14	jstop (fac_max_num)			! Stop time
	Character*20	filext				! Output file extension

c Functions

	Integer*4	UPM_Present
	Integer*4	UPM_Get_Value
	Integer*4	UPM_Get_Float
	Integer*4	Lib$Get_Lun
	Integer*4	CCT_Query_Catalog
	Integer*4	CLI$Present
	Integer*4	CLI$Get_Value

c Local

	Integer*4	rstat, txtlen
	Integer*2	cnex			! next position in command line
	Integer*2	rlen			! length of report string
	Character*34	report_default		! default report filename
	Character*42	repstr			! string for report qualifier
	Character*6	numstr		! string of number for command line
	Integer*2	i, j			! counter
	Integer*4	lun			! logic unit # for script file
	Character*79	blankline
	Character*1	blanks(79) / 79 * ' '/
	Equivalence	(blankline, blanks(1) )
	Character*13	archive /'CSDR$FIRAS_IN'/
	Dictionary 'ccm_cme_catalog_entry'
	Record /query_catalog/		query_cat
	Record /ccm_cme_catalog_entry/	cats(50)

c Counters
	Integer*2	lenfext, lenrepf, rep_i, fex_i, lenrepstr, repstr_i

c Begin

	FCS_Parse = %loc (fcs_normal)
	cmdline(1) = blankline
	cmdline(2) = blankline
	cmdline(3) = blankline
	snum = 0
	repdef = .false.
	cmdline(1)(1:3) = 'FCS'
	clin = 1
	clen = 3
	rstat = UPM_Present ('skyfile')
	If (rstat .EQ. upm_pres) Then
	   rstat = UPM_Get_Value ('skyfile', sfile, txtlen)
	   cnex = clen + 1
	   clen = clen + 9 + txtlen
	   cmdline(1)(cnex:clen) = '/SKYFILE=' // sfile (1:txtlen)
	   If (rstat .NE. ss$_normal) Then
	      Call Lib$Signal (fcs_parserr, %val(1),
	1        ' User must enter file name containing skymap names')
	      FCS_Parse = %loc (fcs_abort)
	   Endif
	Else
	   Call Lib$Signal (fcs_parserr, %val(1),
	1     ' User must enter file name containing skymap names')
	   FCS_Parse = %loc (fcs_abort)
	Endif
	If (FCS_Parse .EQ. %loc (fcs_normal)) Then
	   rstat = UPM_Present ('input')
	   If (rstat .EQ. upm_pres .OR. rstat .EQ. upm_defaulted) Then
	      rstat = UPM_Get_Value ('input', input, txtlen)
	      If (rstat .EQ. ss$_normal) Then
	         Call STR$Upcase (input, input)
	         cnex = clen + 1
	         clen = clen + 10
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 1
	            clen = 10
	         Endif
	         cmdline(clin)(cnex:clen) = '/INPUT=' // input
	      Else
	         Call Lib$Signal (fcs_parserr, %val(1),
	1           ' User must choose FAD or FCF input')
	         FCS_Parse = %loc (fcs_abort)
	      Endif
	   Else
	      Call Lib$Signal (fcs_parserr, %val(1),
	1        ' User must choose FAD or FCF input')
	      FCS_Parse = %loc (fcs_abort)
	   Endif
	Endif
	If (FCS_Parse .EQ. %loc (fcs_normal)) Then
	   rstat = UPM_Present ('channel')
	   If (rstat .EQ. upm_pres .OR. rstat .EQ. upm_defaulted) Then
	      rstat = UPM_Get_Value ('channel', chan, txtlen)
	      If (rstat .EQ. ss$_normal) Then
	         Call STR$Upcase (chan, chan)
	         cnex = clen + 1
	         clen = clen + 11
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 1
	            clen = 11
	         Endif
	         cmdline(clin)(cnex:clen) = '/CHANNEL=' // chan
	      Else
	         Call Lib$Signal (fcs_parserr, %val(1),
	1           ' User must enter a valid channel')
	         FCS_Parse = %loc (fcs_abort)
	      Endif
	   Else
	      Call Lib$Signal (fcs_parserr, %val(1),
	1        ' User must enter a valid channel')
	      FCS_Parse = %loc (fcs_abort)
	   Endif
	Endif
	If (FCS_Parse .EQ. %loc (fcs_normal)) Then
	   rstat = UPM_Present ('scan_mode')
	   If (rstat .EQ. upm_pres .OR. rstat .EQ. upm_defaulted) Then
	      rstat = UPM_Get_Value ('scan_mode', scan, txtlen)
	      If (rstat .EQ. ss$_normal) Then
	         Call STR$Upcase (scan, scan)
	         cnex = clen + 1
	         clen = clen + 13
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 1
	            clen = 13
	         Endif
	         cmdline(clin)(cnex:clen) = '/SCAN_MODE=' // scan
	      Else
	         Call Lib$Signal (fcs_parserr, %val(1),
	1           ' User must enter a valid scan mode')
	         FCS_Parse = %loc (fcs_abort)
	      Endif
	   Else
	      Call Lib$Signal (fcs_parserr, %val(1),
	1        ' User must enter a valid scan mode')
	      FCS_Parse = %loc (fcs_abort)
	   Endif
	Endif
	If (FCS_Parse .EQ. %loc (fcs_normal)) Then
	   report_default = 'FCS_' // chan // scan // '_REPORT.'
	   rstat = UPM_Present ('report')
	   If (rstat .EQ. upm_pres) Then
	      report = .TRUE.
	      rstat = UPM_Get_Value ('report', reportf, txtlen)
	      If (rstat .EQ. upm_absent) Then
	         reportf = report_default
	         rlen = 24
	         repdef = .true.
	         repstr = '/REPORT=' // reportf
	      Else
	         rlen = 8 + txtlen
	         repdef = .false.
	         cnex = clen + 1
	         clen = clen + rlen
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 1
	            clen = rlen
	         Endif
	         Call STR$Upcase (reportf, reportf)
	         cmdline(clin)(cnex:clen) = '/REPORT=' // reportf (1:txtlen)
	      Endif
	   Elseif (rstat .EQ. upm_defaulted) Then
	      report = .TRUE.
	      reportf = report_default
	      repstr = '/REPORT=' // reportf
	      repdef = .true.
	      rlen = 24
	   Else
	      report = .FALSE.
	   Endif

c If the report name is defaulted then we need to read the script file for a
c file extension or time ranges to complete the file name.  The check of the
c time ranges is done in FCS_REPORT_CHECK so that error messages can be
c saved in a report file.
c
c Read namelist file.  Get skymap file names, time ranges, and file extension.

	   Do i=1,fac_max_num
	      tstart(i)  = fac_jstart_default
	      tstop(i)   = fac_jstop_default
	      jstart(i)  = ' '
	      jstop(i)   = ' '
	      skymap(i)  = ' '
	   Enddo
	   filext = ' '
	   rstat = Lib$Get_Lun (lun)
	   If ( rstat .NE. ss$_normal ) Then
	      FCS_Parse = %loc (fcs_abort)
	      Call Lib$Signal (fcs_lunerr, %val(1), %val(rstat))
	   Else
	      Open (Unit=lun, File=sfile, Status='old', Access='sequential',
	1        Iostat=rstat, READONLY)
	      If (rstat .NE. 0) Then
	         FCS_Parse = %loc (fcs_abort)
	         Call Lib$Signal (fcs_openerr, %val(1), %val(rstat))
	      Else
	         Read (lun, Nml=smaps, Iostat=rstat)
	         If (rstat .NE. 0) Then
	            FCS_Parse = %loc (fcs_abort)
	            Call Lib$Signal (fcs_readerr, %val(1), %val(rstat))
	         Else
	            i = 1
	            Do While (i .LE. fac_max_num .AND.
	1                     FCS_Parse .EQ. %loc (fcs_normal))
	               If (skymap(i)(1:1) .NE. ' ') Then
	                  snum = snum + 1
	                  in_skymap(i) = skymap(i)
	                  query_cat.archive_id = archive
	                  query_cat.filename = skymap(i)
	                  rstat = cct_query_catalog (query_cat, cats(1))
	                  If (.NOT. rstat) Then
	                     FCS_Parse = %loc (fcs_abort)
	                     Call Lib$Signal(fcs_query_cat,%val(1),%val(rstat))
	                  Else
	                     If (jstart(i)(1:1) .NE. ' ') Then
	                        j = Index (jstart(i),' ')
	                        If (j .EQ. 0) Then
	                           tstart(i) = jstart(i)
	                        Else
	                           tstart(i)(1:j-1) = jstart(i)(1:j-1)
	                        Endif
	                        j = Index (jstop(i),' ')
	                        If (j .EQ. 0) Then
	                           tstop(i) = jstop(i)
	                        Else
	                           tstop(i)(1:j-1) = jstop(i)(1:j-1)
	                        Endif
	                     Else
	                        Call CT_Binary_To_GMT (cats(1).initial_time,
	1                                                          tstart(i))
	                        Call CT_Binary_To_GMT (cats(1).final_time,
	1                                                         tstop(i))
	                     Endif
	                  Endif
	               Endif
	               i = i + 1
	            Enddo
	            fext = filext
	            Close (lun, Iostat=rstat)
	            If (rstat .NE. 0) Then
	               FCS_Parse = %loc (fcs_abort)
	               Call Lib$Signal (fcs_closerr, %val(1), %val(rstat))
	            Else
	               Call Lib$Free_Lun (lun)
	            Endif
	         Endif
	      Endif
	   Endif
	Endif

c If the script file contains an extension then it is used on the report name.

	If (FCS_Parse .EQ. %loc (fcs_normal)) Then
	   If (repdef) Then
	      If (fext(1:1) .NE. ' ') Then
	         Call STR$Upcase (fext, fext)
	         lenfext = Index (fext,' ') - 1
	         If (lenfext .EQ. 0) lenfext = Len (fext)
	         lenrepf = Len (reportf)
	         rep_i = Index (reportf,' ')
	         fex_i = 1
	         Do While ((rep_i .LE. lenrepf) .AND.
	1                  (fex_i .LE. lenfext))
	            reportf (rep_i:rep_i) = fext (fex_i:fex_i)
	            rep_i = rep_i + 1
	            fex_i = fex_i + 1
	            rlen = rlen + 1
	         Enddo
	         fex_i = 1
	         lenrepstr = Len (repstr)
	         repstr_i = Index (repstr,' ')
	         Do While ((repstr_i .LE. lenrepstr) .AND. (fex_i .LE. lenfext))
	            repstr (repstr_i:repstr_i) = fext (fex_i:fex_i)
	            fex_i = fex_i + 1
	            repstr_i = repstr_i + 1
	         Enddo
	         repdef = .FALSE.
	      Endif
	      cnex = clen + 1
	      clen = clen + rlen
	      If (clen .GE. 80) Then
	         clin = clin + 1
	         cnex = 1
	         clen = rlen
	      Endif
	      cmdline(clin)(cnex:clen) = repstr
	   Endif
	Endif
	Return
	End
