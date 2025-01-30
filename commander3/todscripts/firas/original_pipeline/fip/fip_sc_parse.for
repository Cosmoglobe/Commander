	Integer*4  Function  FIP_SC_Parse  ( current_gmt, data_type, input,
	1                                    sfile, slen, cmdline, filexts,
	2                                    fnum )

c------------------------------------------------------------------------------
c
c	Function FIP_SC_PARSE
c
c	This function parses the command line for FIP_SPECTRA_COADD and
c	opens and reads the script file containing files to process.
c
c	Author:   Larry P. Rosen, Hughes STX, July-November 1994
c       Based in part on FIP_PARSE_SKY by Gene Eplee, GSC
c
c------------------------------------------------------------------------------
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
c		fip include files.
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
c		fip_frequency.txt
c		fip_invoc_sky.txt
c		fut_params.txt
c		upm_stat_msg.txt
c
c       Modifications:
c      
c               Steve Brodd, HSTX, 6/26/95, SPR 12208.  Add fcc_data_type
c               definition (4 for FCF_SKY, 5 for FCF_DSK, 6 for FCF_CAL,
c               7 for FCF_DCL, 8 for FIC_SKY, and 9 for FIC_CAL).
c------------------------------------------------------------------------------
	Implicit None

	Include '($ssdef)'
	Include '(fut_params)'
	Include '(fip_invoc_sky)'
	Include '(fip_frequency)'
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

	External	FIP_Lunerr
	External	FIP_Normal
	External	FIP_Parserr
	External	FIP_Galexc
	External	FIP_Openerr
	External	FIP_Readerr
	External	FIP_Closerr
	External	FIP_Nofilext
	External	FIP_Noinput
C------------------------------------------------------------------------------
C Begin
C
C  Parse the command line.
C
	pstatus = %loc (FIP_Normal)
	cmdline(1) = blankline
	cmdline(2) = blankline
	cmdline(3) = blankline
	cmdline(1)(1:4) = 'FIPA'
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
	      Call lib$signal (FIP_Parserr, %val(1), 'SCRIPT')
	      pstatus = %loc (FIP_Parserr)
	   EndIf
	Else
	   Call lib$signal (FIP_Parserr, %val(1), 'SCRIPT')
	   pstatus = %loc (FIP_Parserr)
	EndIf
c
c  Flags per input file type and qualifiers:
c
c                            Flags                            Qualifiers
c              ---------------------------------------  ---------------------
c  Input file  data_type    fcc_diff         fcc_coadd  /Spectra        /Data
c     FCF_SKY     SKY       fac_not_present   .FALSE.     CALIBRATED       SKY
c     FCF_CAL     CAL       fac_not_present   .FALSE.     CALIBRATED       CAL
c     FIC_SKY     SKY       fac_not_present   .TRUE.      COADD            SKY
c     FIC_CAL     CAL       fac_not_present   .TRUE.      COADD            CAL
c     FCF_DSK     SKY       fac_present       .FALSE.     DIFFERENTIAL     SKY
c     FCF_DCL     CAL       fac_present       .FALSE.     DIFFERENTIAL     CAL
c
	If (pstatus .EQ. %loc (FIP_Normal)) Then
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
	         Call lib$signal (FIP_Parserr, %val(1), 'SPECTRA')
	         pstatus = %loc (FIP_Parserr)
	      EndIf
	   Else
	      Call lib$signal (FIP_Parserr, %val(1), 'SPECTRA')
	      pstatus = %loc (FIP_Parserr)
	   EndIf
	EndIf

c  data_type = 'SKY' or 'CAL'

	If (pstatus .EQ. %loc (FIP_Normal)) Then
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
	         Call lib$signal (FIP_Parserr, %val(1), 'DATA_TYPE')
	         pstatus = %loc (FIP_Parserr)
	      EndIf
	   Else
	      Call lib$signal (FIP_Parserr, %val(1), 'DATA_TYPE')
	      pstatus = %loc (FIP_Parserr)
	   EndIf
	EndIf

c  scan mode

	If (pstatus .EQ. %loc (FIP_Normal)) Then
	   status = upm_present ('SCAN_MODE')
	   If (status .EQ. upm_pres .or. status .EQ. upm_defaulted) Then
	      status = upm_get_value ('SCAN_MODE', scan, len)
	      If (status .EQ. ss$_normal) Then
	         Call str$upcase (scan, scan)
	         cnex = clen + 1
	         clen = clen + 13
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 5
	            clen = 17
	         EndIf
	         cmdline (clin)(cnex:clen) = '/SCAN_MODE=' // scan
	      Else
	         Call lib$signal (FIP_Parserr, %val(1), 'SCAN_MODE')
	         pstatus = %loc (FIP_Parserr)
	      EndIf
	   Else
	      Call lib$signal (FIP_Parserr, %val(1), 'SCAN_MODE')
	      pstatus = %loc (FIP_Parserr)
	   EndIf
	EndIf

c  channel

	If (pstatus .EQ. %loc (FIP_Normal)) Then
	   status = upm_present ('CHANNEL')
	   If (status .EQ. upm_pres .or. status .EQ. upm_defaulted) Then
	      status = upm_get_value ('CHANNEL', chan, len)
	      If (status .EQ. ss$_normal) Then
	         Call str$upcase (chan, chan)
	         cnex = clen + 1
	         clen = clen + 11
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 5
	            clen = 15
	         EndIf
	         cmdline (clin)(cnex:clen) = '/CHANNEL=' // chan
	      Else
	         Call lib$signal (FIP_Parserr, %val(1), 'CHANNEL')
	         pstatus = %loc (FIP_Parserr)
	      EndIf
	   Else
	      Call lib$signal (FIP_Parserr, %val(1), 'CHANNEL')
	      pstatus = %loc (FIP_Parserr)
	   EndIf
	EndIf
c
c  Construct the input file type
c
	If (pstatus .EQ. %loc (FIP_Normal)) Then
	   fcc_scan_mode = chan // scan
	   fcc_outfile (1:4) = 'FIP_'
	   If (.not. fcc_dIff) Then
	      If (fcc_coadd) Then
	         input (1:4) = 'FIC_'
	         fcc_outfile (5:5) = 'I'
   	         If (data_type .EQ. 'SKY') Then
	            input (5:8) = 'SKY_'
	            fcc_outfile (6:8) = 'SK_'
                    fcc_data_type = 8
	         Else
	            input (5:8) = 'CAL_'
	            fcc_outfile (6:8) = 'CL_'
                    fcc_data_type = 9
	         EndIf
	      Else
	         input (1:4) = 'FCF_'
	         fcc_outfile (5:5) = 'C'
	         If (data_type .EQ. 'SKY') Then
	            input (5:8) = 'SKY_'
	            fcc_outfile (6:8) = 'SK_'
                    fcc_data_type = 4
	         Else
	            input (5:8) = 'CAL_'
	            fcc_outfile (6:8) = 'CL_'
                    fcc_data_type = 6
	         EndIf
	      EndIf
	   Else
	      input (1:4) = 'FCF_'
	      If (data_type .EQ. 'SKY') Then
	         input (5:8) = 'DSK_'
                 fcc_data_type = 5
	      Else
	         input (5:8) = 'DCL_'
                 fcc_data_type = 7
	      EndIf
	      fcc_outfile (5:8) = input (5:8)
	   EndIf
	   input (9:12) = fcc_scan_mode			! full chanscan
	   fcc_outfile (9:12) = fcc_scan_mode
	EndIf
c
c  Set the galactic latitude exclusion invocation flags
c
	If (pstatus .EQ. %loc (FIP_Normal)) Then
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
	            pstatus = %loc (FIP_Galexc)
	            Call lib$signal (FIP_Galexc, %val(1), %val(status))
	         EndIf
	      Else
	         pstatus = %loc(FIP_Galexc)
	         Call lib$signal (FIP_Galexc, %val(1), %val(status))
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
	If (pstatus .EQ. %loc (FIP_Normal)) Then
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
	            If (status .EQ. ss$_normal) fcc_lofreq = freq
	         Else
	            fcc_lofreq = fcc_lofreq_default
	         EndIf
	         status = upm_present ('FREQ_RANGE.HIGH')
	         If (status .EQ. upm_pres) Then
	            status = upm_get_longword ('FREQ_RANGE.HIGH', freq)
	            If (status .EQ. ss$_normal) fcc_hIfreq = freq
	         Else
	            fcc_hifreq = fcc_hifreq_default
	         EndIf
	         write (num_str1,20) fcc_lofreq		! convert # to string
  20	         format (i2)			! eg. ' 2', '22'
	         cnex = clen + 1
	         clen = clen + 8
	         If (clen .GE. 80) Then
	            clin = clin + 1
	            cnex = 5
	            clen = 12
	         EndIf
	         cmdline (clin)(cnex:clen) = '(LOW=' // num_str1 // ','
	         write (num_str2,30) fcc_hifreq		! convert # to string
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
	      fcc_lofreq = fcc_lofreq_default
	      fcc_hifreq = fcc_hifreq_default
	      write (num_str1,20) fcc_lofreq		! convert # to string
	      cnex = clen + 1
	      clen = clen + 8
	      If (clen .GE. 80) Then
	         clin = clin + 1
	         cnex = 5
	         clen = 12
	      EndIf
	      cmdline (clin)(cnex:clen) = '(LOW=' // num_str1 // ','
	      write (num_str2,30) fcc_hifreq		! convert # to string
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

	If (pstatus .EQ. %loc (FIP_Normal)) Then
	   Do i=1, fac_max_num
	      infilext (i)  = ' '
	      filexts (i) = ' '
	   EndDo
	   outfilext = ' '
	   status = lib$get_lun (lun)
	   If ( status .NE. ss$_normal ) Then
	      pstatus = %loc (FIP_Lunerr)
	      Call lib$signal (FIP_Lunerr, %val(1), %val(status))
	   Else
	      Open ( Unit=lun, File=sfile(1:slen), Status='old',
	1            Access='sequential', Iostat=status, READONLY )
	      If (status .NE. 0) Then
	         pstatus = %loc (FIP_Openerr)
	         Call lib$signal (FIP_Openerr, %val(2), sfile(1:slen),
	1                         %val(status))
	      Else
	         Read (lun, nml=sfiles, iostat=status)
	         If (status .NE. 0) Then
	            pstatus = %loc (FIP_Readerr)
	            Call lib$signal (FIP_Readerr, %val(2), sfile(1:slen),
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
	               pstatus = %loc (FIP_Closerr)
	               Call lib$signal (FIP_Closerr, %val(2), sfile(1:slen),
	1                            %val(status))
	            Else
	               Call lib$free_lun (lun)
	               If (outfilext(1:1) .EQ. ' ') Then
	                  pstatus = %loc (FIP_Nofilext)
	                  Call lib$signal (FIP_Nofilext)
	               Else
	                  Call str$trim (outfilext, outfilext, len)
	                  fcc_outlen = len + 14
	                  fcc_outfile (13:fcc_outlen) = '.' // outfilext(1:len)
	               EndIf
	               If (fnum .EQ. 0) Then
	                  pstatus = %loc (FIP_Noinput)
	                  Call lib$signal (FIP_Noinput)
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
	If (pstatus .EQ. %loc (FIP_Normal)) Then
	   status = upm_present ('REPORT')
	   If (status .EQ. upm_negated) Then
	      fcc_report = fac_not_present
	   Else
	      fcc_report = fac_present
	      status = upm_get_value ('REPORT', fcc_report_file, fcc_replen)
	      If (status .EQ. upm_absent) Then  !  Set default report file name
	         fcc_report_file = 'FIP_SPECTRA_COADD_' // input // '.REP_' //
	1                          current_gmt(1:9)
	      Else			   !  Report file name from invocation
	         Call str$upcase (fcc_report_file, fcc_report_file)
	      EndIf
	      Call str$trim (fcc_report_file, fcc_report_file, fcc_replen)
	   EndIf
	EndIf

	FIP_SC_Parse = pstatus

	Return
	End
