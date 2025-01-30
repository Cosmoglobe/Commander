	Integer*4  Function  FFP_SC_Parse  ( current_gmt, data_type, input,
	1                                    sfile, slen, cmdline, filexts,
	2                                    fnum )
c------------------------------------------------------------------------------
c
c	Function FFP_SC_PARSE
c
c	This function parses the command line for FFP_SPECTRA_COADD and
c	opens and reads the script file containing files to process.
c
c	Author:  S. Brodd, HSTX, 3/21/96
c
c	Input:
c		none
c
c	Output:
c		character*14	current_gmt		! GMT of invocation
c		character*3	data_type		! Cal or sky type
c		character*12	input			! Input filename base
c		character*80	sfile			! Script file name
c		integer*2	slen			! Length of sfile name
c		character*79	cmdline (3)		! Command line 
c		character*20	filexts (fac_max_num)	! Input file extensions
c		integer*2	fnum			! Number of input files
c
c		Plus quite a few parameters set in the common block in the
c		ffp include files.
c
c	Subroutines called:
c		ct_binary_to_gmt
c		sys$gettim
c		str$upcase
c		upm_get_float
c		upm_get_longword
c		upm_get_value
c		upm_present
c		lib$get_lun
c
c	Incude files:
c		$ssdef
c		ffp_frequency.txt
c		ffp_invoc_sky.txt
c		fut_params.txt
c		upm_stat_msg.txt
c
c       Modifications:
c               S. Brodd, HSTX, 1/23/97, SPR 12339. Now delivering only
c               merged channel and scan mode destriped spectra.
c      
c------------------------------------------------------------------------------
	Implicit None

	Include '($ssdef)'
	Include '(fut_params)'
	Include '(ffp_invoc_sky)'
	Include '(ffp_frequency)'
	Include '(upm_stat_msg)'

c Passed parameters (plus infilext of namelist is passed too)

	Character*14	current_gmt		! GMT time of invocation
	Character*3	data_type		! Cal or sky type
	Character*12	input			! Input filename base
	Character*80	sfile			! Script file name
	Integer*2	slen			! Length of script file name
	Character*79	cmdline (3)		! Command line with defaults
	Character*20	filexts (fac_max_num)	! Input file extensions
	Integer*2	fnum			! Number of input files

c Local
	Integer*4	pstatus			!  parse return status
	Integer*4	current_time(2)		!  VAX ADT time of invocation
	Integer*4	status			!  return status
	Integer*2	clin			!  which line of command line
	Integer*2	clen		!  length of current command line
	Character*12	spectra_type		!  diff., cal., or coadd
	Integer*2	len			!  length of string
	Character*2	scan			!  scan mode
	Character*2	chan			!  channel
	Integer*4	lun			!  logical unit for script
	Integer*2	i			!  counter
	Real*4		glat			!  input latitude cutoff
	Integer*4	freq			!  input frequency cutoff
	Character*79	blankline
	Character*1	blanks(79) / 79 * ' '/
	Equivalence	(blankline, blanks(1) )
	Integer*2	cnex
	Character*7	num_str				!  number as a string
	Character*2	num_str1			!  number as a string
	Character*3	num_str2			!  number as a string

c Namelist

	Namelist / sfiles / infilext, outfilext
	Character*20	infilext (fac_max_num)	! Input files extensions
	Character*20	outfilext		! Output file extension

c Functions

	Integer*4	upm_get_float
	Integer*4	upm_get_longword
	Integer*4	upm_get_value
	Integer*4	upm_present
	Integer*4	lib$get_lun

c Externals

	External	FFP_Lunerr
	External	FFP_Normal
	External	FFP_Parserr
	External	FFP_Galexc
	External	FFP_Openerr
	External	FFP_Readerr
	External	FFP_Closerr
	External	FFP_Nofilext
	External	FFP_Noinput
C
C  Parse the command line.
C
	pstatus = %loc (FFP_Normal)
	cmdline(1) = blankline
	cmdline(2) = blankline
	cmdline(3) = blankline
	cmdline(1)(1:4) = 'FFPA'
	clin = 1
	clen = 4
	cnex = 5
c
c  Get the time of the invocation
c
	Call sys$gettim (current_time)
	Call ct_binary_to_gmt (current_time, current_gmt)
c
c  Get the script file name.
c
	status = upm_present ('SCRIPT')
	If (status .EQ. upm_pres) Then
	   status = upm_get_value ('SCRIPT', sfile, slen)
	   If (status .EQ. ss$_normal) Then
	      Call str$upcase (sfile, sfile)
	      clen = clen + 8 + slen
	      cmdline (1)(cnex:clen) = '/SCRIPT=' // sfile (1:slen)
	      cnex = clen + 1
	   Else
	      Call lib$signal (FFP_Parserr, %val(1), 'SCRIPT')
	      pstatus = %loc (FFP_Parserr)
	   EndIf
	Else
	   Call lib$signal (FFP_Parserr, %val(1), 'SCRIPT')
	   pstatus = %loc (FFP_Parserr)
	EndIf
c
c  Flags per input file type and qualifiers:
c
c                            Flags                            Qualifiers
c              ---------------------------------------  ---------------------
c  Input file  data_type    fcc_diff         fcc_coadd  /Spectra        /Data
c     FSL_SKY     SKY       fac_not_present   .FALSE.     CALIBRATED       SKY
c     FSL_CAL     CAL       fac_not_present   .FALSE.     CALIBRATED       CAL
c     FIL_SKY     SKY       fac_not_present   .TRUE.      COADD            SKY
c     FIL_CAL     CAL       fac_not_present   .TRUE.      COADD            CAL
c     FSL_DSK     SKY       fac_present       .FALSE.     DIFFERENTIAL     SKY
c     FSL_DCL     CAL       fac_present       .FALSE.     DIFFERENTIAL     CAL
c
	If (pstatus .EQ. %loc (FFP_Normal)) Then
	   status = upm_present ('SPECTRA')
	   If (status .EQ. upm_pres) Then
	      status = upm_get_value ('SPECTRA', spectra_type, len)
	      If (status .EQ. ss$_normal) Then
	         Call str$upcase (spectra_type, spectra_type)
	         clen = cnex + 9 + len
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 5
	            clen = 13 + len
	         EndIf
	         cmdline (clin)(cnex:clen) = '/SPECTRA=' // spectra_type(1:len)
	         If (spectra_type (1:len) .EQ. 'DIFFERENTIAL') Then
	            fcc_diff = fac_present
	            fcc_coadd = .FALSE.
	         ElseIf(spectra_type (1:len) .EQ. 'CALIBRATED') Then
	            fcc_diff = fac_not_present
	            fcc_coadd = .FALSE.
	         Else				! 'COADD'
	            fcc_diff = fac_not_present
	            fcc_coadd = .TRUE.
	         EndIf
	      Else
	         Call lib$signal (FFP_Parserr, %val(1), 'SPECTRA')
	         pstatus = %loc (FFP_Parserr)
	      EndIf
	   Else
	      Call lib$signal (FFP_Parserr, %val(1), 'SPECTRA')
	      pstatus = %loc (FFP_Parserr)
	   EndIf
	EndIf

c  data_type = 'SKY' or 'CAL'

	If (pstatus .EQ. %loc (FFP_Normal)) Then
	   status = upm_present ('DATA')
	   If (status .EQ. upm_pres) Then
	      status = upm_get_value ('DATA', data_type, len)
	      If (status .EQ. ss$_normal) Then
	         Call str$upcase (data_type, data_type)
	         cnex = clen + 1
	         clen = clen + 9
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 5
	            clen = 13
	         EndIf
	         cmdline (clin)(cnex:clen) = '/DATA=' // data_type
	      Else
	         Call lib$signal (FFP_Parserr, %val(1), 'DATA_TYPE')
	         pstatus = %loc (FFP_Parserr)
	      EndIf
	   Else
	      Call lib$signal (FFP_Parserr, %val(1), 'DATA_TYPE')
	      pstatus = %loc (FFP_Parserr)
	   EndIf
	EndIf

c  scan mode

	If (pstatus .EQ. %loc (FFP_Normal)) Then
	   status = upm_present ('SCAN_MODE')
	   If (status .EQ. upm_pres .or. status .EQ. upm_defaulted) Then
	      status = upm_get_value ('SCAN_MODE', scan, len)
	      If (status .EQ. ss$_normal) Then
	         Call str$upcase (scan, scan)
                 if (scan(1:2) .eq. 'SS') then
                    fcc_fsmode = 1
                 elseif (scan(1:2) .eq. 'SF') then
                    fcc_fsmode = 2
                 elseif (scan(1:2) .eq. 'LS') then
                    fcc_fsmode = 3
                 elseif (scan(1:2) .eq. 'LF') then
                    fcc_fsmode = 4
                 elseif (scan(1:2) .eq. 'FS') then
                    fcc_fsmode = 5
                 elseif (scan(1:2) .eq. 'FL') then
                    fcc_fsmode = 6
                 endif
	         cnex = clen + 1
	         clen = clen + 13
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 5
	            clen = 17
	         EndIf
	         cmdline (clin)(cnex:clen) = '/SCAN_MODE=' // scan
	      Else
	         Call lib$signal (FFP_Parserr, %val(1), 'SCAN_MODE')
	         pstatus = %loc (FFP_Parserr)
	      EndIf
	   Else
	      Call lib$signal (FFP_Parserr, %val(1), 'SCAN_MODE')
	      pstatus = %loc (FFP_Parserr)
	   EndIf
	EndIf

c  channel

	If (pstatus .EQ. %loc (FFP_Normal)) Then
	   status = upm_present ('CHANNEL')
	   If (status .EQ. upm_pres .or. status .EQ. upm_defaulted) Then
	      status = upm_get_value ('CHANNEL', chan, len)
	      If (status .EQ. ss$_normal) Then
	         Call str$upcase (chan, chan)
                 if (chan(1:2) .eq. 'RH') then
                    fcc_fchan = 1
                 elseif (chan(1:2) .eq. 'RL') then
                    fcc_fchan = 2
                 elseif (chan(1:2) .eq. 'LH') then
                    fcc_fchan = 3
                 elseif (chan(1:2) .eq. 'LL') then
                    fcc_fchan = 4
                 endif
	         cnex = clen + 1
	         clen = clen + 11
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 5
	            clen = 15
	         EndIf
	         cmdline (clin)(cnex:clen) = '/CHANNEL=' // chan
	      Else
	         Call lib$signal (FFP_Parserr, %val(1), 'CHANNEL')
	         pstatus = %loc (FFP_Parserr)
	      EndIf
	   Else
	      Call lib$signal (FFP_Parserr, %val(1), 'CHANNEL')
	      pstatus = %loc (FFP_Parserr)
	   EndIf
	EndIf
c
c  Construct the input file type
c
	If (pstatus .EQ. %loc (FFP_Normal)) Then
	   fcc_scan_mode = chan // scan
	   fcc_outfile (1:4) = 'FFP_'
	   If (.not. fcc_dIff) Then
	      If (fcc_coadd) Then
	         input (1:4) = 'FIL_'
	         fcc_outfile (5:5) = 'I'
   	         If (data_type .EQ. 'SKY') Then
	            input (5:8) = 'SKY_'
	            fcc_outfile (6:8) = 'SK_'
                    fcc_data_type = 6
	         Else
	            input (5:8) = 'CAL_'
	            fcc_outfile (6:8) = 'CL_'
                    fcc_data_type = 7
	         EndIf
	      Else
	         input (1:4) = 'FSL_'
	         fcc_outfile (5:5) = 'C'
	         If (data_type .EQ. 'SKY') Then
	            input (5:8) = 'SKY_'
	            fcc_outfile (6:8) = 'SK_'
                    fcc_data_type = 2
	         Else
	            input (5:8) = 'CAL_'
	            fcc_outfile (6:8) = 'CL_'
                    fcc_data_type = 4
	         EndIf
	      EndIf
	   Else
	      input (1:4) = 'FSL_'
	      If (data_type .EQ. 'SKY') Then
	         input (5:8) = 'DSK_'
                 fcc_data_type = 3
	      Else
	         input (5:8) = 'DCL_'
                 fcc_data_type = 5
	      EndIf
	      fcc_outfile (5:8) = input (5:8)
	   EndIf
	   input (9:12) = fcc_scan_mode			! full chanscan
	   fcc_outfile (9:12) = fcc_scan_mode
	EndIf
c
c  Set the galactic latitude exclusion invocation flags
c
	If (pstatus .EQ. %loc (FFP_Normal)) Then
	   status = upm_present ('EXC_GLAT')
	   If (status .EQ. upm_pres) Then
	      fcc_galexc = fac_present
	      status = upm_get_float ('EXC_GLAT', glat)
	      If (status .EQ. ss$_normal) Then
	         If ((glat .GE. 0.0)  .or.  (glat .LE. fcc_glat_default)) Then
	            fcc_glat = glat
	            write (num_str,10) glat		! convert # to string
  10	            format (f6.2)			! eg. '359.95'
	            cnex = clen + 1
	            clen = clen + 16
	            If (clen .GE. 80) Then
	               clin = clin + 1
	               cnex = 5
	               clen = 20
	            EndIf
	            cmdline (clin)(cnex:clen) = '/EXC_GLAT=' // num_str
	         Else
	            pstatus = %loc (FFP_Galexc)
	            Call lib$signal (FFP_Galexc, %val(1), %val(status))
	         EndIf
	      Else
	         pstatus = %loc(FFP_Galexc)
	         Call lib$signal (FFP_Galexc, %val(1), %val(status))
	      EndIf
	   Else
	      fcc_galexc = fac_not_present
	      fcc_glat = fcc_glat_default
	      write (num_str,10) fcc_glat		! convert # to string
	      cnex = clen + 1
	      clen = clen + 17
	      If (clen .GE. 80) Then
	         clin = clin + 1
	         cnex = 5
	         clen = 21
	      EndIf
	      cmdline (clin)(cnex:clen) = '/EXC_GLAT=' // num_str
	   EndIf
	EndIf
c
c  Set the frequency range invocation flags.
c
	If (pstatus .EQ. %loc (FFP_Normal)) Then
	   cnex = clen + 1
	   clen = clen + 12
	   If (clen .GE. 80) Then
	      clin = clin + 1
	      cnex = 5
	      clen = 16
	   EndIf
	   cmdline (clin)(cnex:clen) = '/FREQ_RANGE='
	   status = upm_present ('FREQ_RANGE')
	   If (status .EQ. upm_pres) Then
	      status = upm_present ('FREQ_RANGE.ALL')
	      If (status .EQ. upm_pres) Then
	         fcc_freq = fac_not_present
	         cnex = clen + 1
	         clen = clen + 3
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 5
	            clen = 7
	         EndIf
	         cmdline (clin)(cnex:clen) = 'ALL'
	      Else
	         fcc_freq = fac_present
	         status = upm_present ('FREQ_RANGE.LOW')
	         If (status .EQ. upm_pres) Then
	            status = upm_get_longword ('FREQ_RANGE.LOW', freq)
	            If (status .EQ. ss$_normal) fcc_jlo = freq
	         Else
	            fcc_jlo = fcc_jlo_default
	         EndIf
	         status = upm_present ('FREQ_RANGE.HIGH')
	         If (status .EQ. upm_pres) Then
	            status = upm_get_longword ('FREQ_RANGE.HIGH', freq)
	            If (status .EQ. ss$_normal) fcc_jhi = freq
	         Else
	            fcc_jhi = fcc_jhi_default
	         EndIf
	         write (num_str1,20) fcc_jlo ! convert # to string
  20	         format (i2)			! eg. ' 2', '22'
	         cnex = clen + 1
	         clen = clen + 8
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 5
	            clen = 12
	         EndIf
	         cmdline (clin)(cnex:clen) = '(LOW=' // num_str1 // ','
	         write (num_str2,30) fcc_jhi ! convert # to string
  30	         format (i3)			! eg. ' 80', '100'
	         cnex = clen + 1
	         clen = clen + 9
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 5
	            clen = 13
	         EndIf
	         cmdline (clin)(cnex:clen) = 'HIGH=' // num_str2 // ')'
	      EndIf
	   Else
	      fcc_freq = fac_not_present
	      fcc_jlo = fcc_jlo_default
	      fcc_jhi = fcc_jhi_default
	      write (num_str1,20) fcc_jlo ! convert # to string
	      cnex = clen + 1
	      clen = clen + 8
	      If (clen .GE. 80) Then
	         clin = clin + 1
	         cnex = 5
	         clen = 12
	      EndIf
	      cmdline (clin)(cnex:clen) = '(LOW=' // num_str1 // ','
	      write (num_str2,30) fcc_jhi ! convert # to string
	      cnex = clen + 1
	      clen = clen + 9
	      If (clen .GE. 80) Then
	         clin = clin + 1
	         cnex = 5
	         clen = 13
	      EndIf
	      cmdline (clin)(cnex:clen) = 'HIGH=' // num_str2 // ')'
	   EndIf
	EndIf
c
c Read namelist file.  Get input file names.
c
c	Namelist / sfiles / infilext, outfilext
c	Character*20	infilext (fac_max_num)	! Input files extensions
c	Character*20	outfilext		! Output file extension

	If (pstatus .EQ. %loc (FFP_Normal)) Then
	   Do i=1, fac_max_num
	      infilext (i)  = ' '
	      filexts (i) = ' '
	   EndDo
	   outfilext = ' '
	   status = lib$get_lun (lun)
	   If ( status .NE. ss$_normal ) Then
	      pstatus = %loc (FFP_Lunerr)
	      Call lib$signal (FFP_Lunerr, %val(1), %val(status))
	   Else
	      Open ( Unit=lun, File=sfile(1:slen), Status='old',
	1            Access='sequential', Iostat=status, READONLY )
	      If (status .NE. 0) Then
	         pstatus = %loc (FFP_Openerr)
	         Call lib$signal (FFP_Openerr, %val(2), sfile(1:slen),
	1                         %val(status))
	      Else
	         Read (lun, nml=sfiles, iostat=status)
	         If (status .NE. 0) Then
	            pstatus = %loc (FFP_Readerr)
	            Call lib$signal (FFP_Readerr, %val(2), sfile(1:slen),
	1                            %val(status))
	         Else
	            fnum = 0		! Count the number of input files
	            i = 1
	            Do while (i .LE. fac_max_num)
	               If (infilext(i)(1:1) .NE. ' ') Then
	                  fnum = fnum + 1
	                  filexts (fnum) = infilext (i)
	               EndIf
	               i = i + 1
	            EndDo
	            Close (lun, iostat=status)
	            If (status .NE. 0) Then
	               pstatus = %loc (FFP_Closerr)
	               Call lib$signal (FFP_Closerr, %val(2), sfile(1:slen),
	1                            %val(status))
	            Else
	               Call lib$free_lun (lun)
	               If (outfilext(1:1) .EQ. ' ') Then
	                  pstatus = %loc (FFP_Nofilext)
	                  Call lib$signal (FFP_Nofilext)
	               Else
	                  Call str$trim (outfilext, outfilext, len)
	                  fcc_outlen = len + 14
	                  fcc_outfile (13:fcc_outlen) = '.' // outfilext(1:len)
	               EndIf
	               If (fnum .EQ. 0) Then
	                  pstatus = %loc (FFP_Noinput)
	                  Call lib$signal (FFP_Noinput)
	               EndIf
	            EndIf
	         EndIf
	      EndIf
	   EndIf
	EndIf
c
c  Get report name and set the processing report flags.
c  If report name is defaulted, construct it.
c
	If (pstatus .EQ. %loc (FFP_Normal)) Then
	   status = upm_present ('REPORT')
	   If (status .EQ. upm_negated) Then
	      fcc_report = fac_not_present
	   Else
	      fcc_report = fac_present
	      status = upm_get_value ('REPORT', fcc_report_file, fcc_replen)
	      If (status .EQ. upm_absent) Then  !  Set default report file name
	         fcc_report_file = 'FFP_SPECTRA_COADD_' // input // '.REP_' //
	1                          current_gmt(1:9)
	      Else			   !  Report file name from invocation
	         Call str$upcase (fcc_report_file, fcc_report_file)
	      EndIf
	      Call str$trim (fcc_report_file, fcc_report_file, fcc_replen)
	   EndIf
	EndIf

	FFP_SC_Parse = pstatus

	Return
	End
