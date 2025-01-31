	Integer*4  Function  FIP_PARSE_COV  ( filext, chan, scan,
	1                                     freq_low, freq_high,
	1                                     report, repfile, cmdline,
	1                                     curtime )

C Function to parse the fip_covar command line.
C
C Author: Larry P. Rosen, July 1993, Hughes STX
C
C    Input: None
C    Output:
C       Character*20  filext               ! Input file name extension
C       Character*2   chan                 ! Channel to process (RH,RL,LH,LL)
C       Character*2   scan                 ! Scan mode to process (SS,SF,LS,LF)
C       Real*4        freq_low             ! Low frequency cut-off (icm)
C       Real*4        freq_high            ! High frequency cut-off (icm)
C       Logical*1     report               ! Flag whether to write report
C       Character*42  repfile              ! File name for report
C       Character*79  cmdline (3)          ! Command line with defaults
C       Integer*4     curtime(2)           ! System time in ADT

	Implicit None

C  Include Files

	Include       '($ssdef)'
	Include       '(upm_stat_msg)'

C Passed Parameters (all output)

	Character*20  filext               ! Input file name extension
	Character*2   chan                 ! Channel to process (RH,RL,LH,LL)
	Character*2   scan                 ! Scan mode to process (SS,SF,LS,LF)
	Real*4        freq_low             ! Low frequency cut-off (icm)
	Real*4        freq_high            ! High frequency cut-off (icm)
	Logical*1     report               ! Flag whether to write report
	Character*42  repfile              ! File name for report
	Character*79  cmdline (3)          ! Command line with defaults
	Integer*4     curtime (2)          ! System time in ADT

C Functions

	Integer*4     UPM_Present
	Integer*4     UPM_Get_Value
	Integer*4     UPM_Get_Float
	Integer*4     Lib$Get_Lun
	Integer*4     CCT_Query_Catalog
	Integer*4     CLI$Present
	Integer*4     CLI$Get_Value
	Integer*4     Sys$GetTim

C Local

	Character*79  blankline
	Character*1   blanks(79) / 79 * ' '/
	Equivalence   ( blankline, blanks(1) )
	Integer*4     rstat             ! Return status from function calls
	Integer*4     txtlen            ! Length of text read from command line
	Integer*2     clin              ! Which line of command line
	Integer*2     clen              ! Length of current command line
	Integer*2     cnex              ! Next position in command line
	Character*6   lowstr            ! Low frequency converted to string
	Character*6   highstr           ! High frequency converted to string
	Character*42  default_repfile   ! Default file name for report
	Integer*2     replen            ! Length of report name
	Character*14  curgmt            ! System time in gmt
	Character*50  repstr            ! String for report name

C External

	External	fip_normal
	External	fip_parserr
	External	fip_abort
	External	fip_freqerr
	External	fip_gettimerr

C Begin

	FIP_PARSE_COV = %loc (fip_normal)
	cmdline (1) = blankline
	cmdline (2) = blankline
	cmdline (3) = blankline
	cmdline(1)(1:4) = 'FIPC'
	clin = 1
	clen = 4

C Parse FILE_EXT qualifier.

	rstat = UPM_Present ('file_ext')
	If (rstat .EQ. upm_pres) Then
	   rstat = UPM_Get_Value ('file_ext', filext, txtlen)
	   cnex = clen + 1
	   clen = clen + 11 + txtlen
	   cmdline (clin)(cnex:clen) = ' /FILE_EXT=' // filext (1:txtlen)
	   If (rstat .NE. ss$_normal) Then
	      Call Lib$Signal (fip_parserr, %val(1), 'FILE_EXT')
	      FIP_PARSE_COV = %loc (fip_abort)
	   Endif
	Else
	   Call Lib$Signal (fip_parserr, %val(1), 'FILE_EXT')
	   FIP_PARSE_COV = %loc (fip_abort)
	Endif

C Parse CHANNEL qualifier.

	If (FIP_PARSE_COV .EQ. %loc (fip_normal)) Then
	   rstat = UPM_Present ('channel')
	   If (rstat .EQ. upm_pres .OR. rstat .EQ. upm_defaulted) Then
	      rstat = UPM_Get_Value ('channel', chan, txtlen)
	      If (rstat .EQ. ss$_normal) Then
	         Call STR$Upcase (chan, chan)
	         cnex = clen + 1
	         clen = clen + 12
	         cmdline (clin)(cnex:clen) = ' /CHANNEL=' // chan
	      Else
	         Call Lib$Signal (fip_parserr, %val(1), 'CHANNEL')
	         FIP_PARSE_COV = %loc (fip_abort)
	      Endif
	   Else
	      Call Lib$Signal (fip_parserr, %val(1), 'CHANNEL')
	      FIP_PARSE_COV = %loc (fip_abort)
	   Endif
	Endif

C Parse SCAN_MODE qualifier.

	If (FIP_PARSE_COV .EQ. %loc (fip_normal)) Then
	   rstat = UPM_Present ('scan_mode')
           If (rstat .EQ. upm_pres .OR. rstat .EQ. upm_defaulted) Then
	      rstat = UPM_Get_Value ('scan_mode', scan, txtlen)
	      If (rstat .EQ. ss$_normal) Then
	         Call STR$Upcase (scan, scan)
	         cnex = clen + 1
	         clen = clen + 14
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 5
	            clen = 18
	         Endif
	         cmdline (clin)(cnex:clen) = ' /SCAN_MODE=' // scan
	      Else
	         Call Lib$Signal (fip_parserr, %val(1), 'SCAN_MODE')
	         FIP_PARSE_COV = %loc (fip_abort)
	      Endif
	   Else
	      Call Lib$Signal (fip_parserr, %val(1), 'SCAN_MODE')
	      FIP_PARSE_COV = %loc (fip_abort)
	   Endif
	Endif

C Parse FREQ_RANGE qualifier.

	If (FIP_PARSE_COV .EQ. %loc (fip_normal)) Then
	   rstat = UPM_Present ('freq_range')
	   If (rstat .EQ. upm_pres) Then
	      rstat = upm_present ('freq_range.low')
	      If (rstat .eq. upm_pres) Then
	         rstat = upm_get_float ('freq_range.low', freq_low)
	         If (rstat .EQ. ss$_normal) Then
	            cnex = clen + 1
	            clen = clen + 25
	            If (clen .GE. 80) Then
	               clin = clin + 1
	               cnex = 5
	               clen = 29
	            Endif
	            Write (lowstr, 10) freq_low
  10	            Format (F6.2)
	            cmdline (clin)(cnex:clen) = ' /FREQ_RANGE=(LOW=' //
	1                                       lowstr // ','
	         Else
	            Call Lib$Signal (fip_parserr, %val(1),
	1                            'FREQ_RANGE.LOW')
	            FIP_PARSE_COV = %loc (fip_abort)
	         Endif
	      Else
	         Call Lib$Signal (fip_parserr, %val(1),
	1                         'FREQ_RANGE.LOW')
	         FIP_PARSE_COV = %loc (fip_abort)
	      Endif
	      rstat = upm_present ('freq_range.high')
	      If (rstat .eq. upm_pres) Then
	         rstat = upm_get_float ('freq_range.high', freq_high)
	         If (rstat .EQ. ss$_normal) Then
	            cnex = clen + 1
	            clen = clen + 13
	            If (clen .GE. 80) Then
	               clin = clin + 1
	               cnex = 5
	               clen = 17
	            Endif
	            Write (highstr, 10) freq_high
	            cmdline (clin)(cnex:clen) = ' HIGH=' //
	1                                       highstr // ')'
	         Else
	            Call Lib$Signal (fip_parserr, %val(1),
	1                            'FREQ_RANGE.HIGH')
	            FIP_PARSE_COV = %loc (fip_abort)
	         Endif
	      Else
	         Call Lib$Signal (fip_parserr, %val(1),
	1                         'FREQ_RANGE.HIGH')
	         FIP_PARSE_COV = %loc (fip_abort)
	      Endif
	      If (FIP_PARSE_COV .EQ. %loc (fip_normal)) Then
	         If (freq_high .LE. freq_low) Then
	            Call Lib$Signal (fip_freqerr, %val(2),
	1                            lowstr, highstr)
	            FIP_PARSE_COV = %loc (fip_abort)
	         Endif
	      Endif
	   Endif
	Endif

C Get the system time for report name.

	rstat = Sys$GetTim ( curtime )
	If ( rstat .EQ. ss$_normal ) Then
	   Call CT_Binary_To_GMT ( curtime, curgmt )
	Else
	   FIP_PARSE_COV = %loc (fip_abort)
	   Call Lib$Signal (fip_gettimerr, %val(1), %val(rstat))
	Endif

C Parse REPORT qualifier.

	default_repfile = 'FIP_COV_' // chan // scan // '_' //
	1                 filext (4:10) // '_' // filext (12:18) //
	1                 '.REP_' // curgmt (1:9)

	If (FIP_PARSE_COV .EQ. %loc (fip_normal)) Then
	   rstat = UPM_Present ('report')
	   If (rstat .EQ. upm_pres) Then
	      report = .TRUE.
	      rstat = upm_get_value ('report', repfile, txtlen)
	      If (rstat .EQ. upm_absent) Then
	         repfile = default_repfile
	         replen = 50
	      Else
	         replen = 8 + txtlen
	      Endif
	      repstr = ' /REPORT=' // repfile
	   Elseif (rstat .EQ. upm_defaulted) Then
	      report = .TRUE.
	      repfile = default_repfile
	      repstr = ' /REPORT=' // repfile
	      replen = 50
	   Else
	      report = .FALSE.
	      repstr = ' /NOREPORT'
	      replen = 9
	   Endif
	   If (clen .GT. 79-replen) Then
	      clin = clin + 1
	      cnex = 5
	      clen = 4
	   Endif
	   clen = clen + replen
	   cmdline (clin)(cnex:clen) = repstr (1:replen)
	Endif
	Return
	End
