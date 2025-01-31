	Program FFP_SPECTRA_COADD

c------------------------------------------------------------------------------
c
c	Program FFP_SPECTRA_COADD
c
c	This program reformats FSL and FIL sky and cal data into the Final
c	Release ADB format.  FFP reformats FSL calibrated spectra and
c       differential spectra skymaps into FFP skymaps, and calibration data
c       sets.  The FFP output is ready for conversion to FITS files.  The
c       input/output data sets are:
c          InputFile         Cal/Sky           Output
c          ----------        -------        --------------
c          FSL_SKY_xxxx        Sky           FFP_CSK_xxxx
c          FSL_CAL_xxxx        Cal           FFP_CCL_xxxx
c          FSL_DSK_xxxx        Sky           FFP_DSK_xxxx
c          FSL_DCL_xxxx        Cal           FFP_DCL_xxxx
c          FIL_SKY_xxxx        Sky           FFP_ISK_xxxx
c          FIL_CAL_xxxx        Cal           FFP_ICL_xxxx
c
c	Author:  S. Brodd, HSTX, 3/21/96, SPR 12316
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
c------------------------------------------------------------------------------
c
c	Subroutines called:
c		cut_display_banner
c		cut_register_version
c		ffp_sc_parse
c		ffp_sc_init_report
c		ffp_sc_process_sky
c		ffp_sc_process_cal
c		ffp_sc_update_report
c		ct_init
c		lib$signal
c
c	Include files:
c		ct$library:ctuser.inc
c		ffp_invoc_sky.txt
c		fut_error.txt
c		fut_params.txt
c
c------------------------------------------------------------------------------
	Implicit None

	Include 'ct$library:ctuser.inc'
	Include '(fut_error)'
	Include '(fut_params)'
	Include '(ffp_invoc_sky)'

c Functions

	Integer*4	CUT_Display_Banner
	Integer*4	CUT_Register_Version
	Integer*4	FFP_SC_Parse
	Integer*4	FFP_SC_Init_Report
	Integer*4	FFP_SC_Process_Sky
	Integer*4	FFP_SC_Process_Cal
	Integer*4	FFP_SC_Update_Report
	integer*4	fut_free_lun

c Local

	Integer*4	status		!  processing status
	Integer*4	rstatus		!  return status
	Integer*4	terminal /6/	!  terminal lun
	Character*14	current_gmt		! GMT time of invocation
	Character*3	data_type		! Cal or sky type
	Character*12	input			! Input filename base
	Character*80	sfile			! Script file name
	Integer*2	slen			! Length of script file name
	Character*79	cmdline (3)		! Command line with defaults
	Character*20	filexts (fac_max_num)	! Input file extensions
	Integer*2	fnum			! Number of input files
	Integer*2	ct_stat(20)		!  Cobetrieve return status
	Character*13	archin /'CSDR$FIRAS_IN'/	! Input archive
	Character*14	archout /'CSDR$FIRAS_OUT'/	! Output archive
	Integer*4	nrec			! Number of records processed
	Integer*4	npix			! Number of pixels processed
	Integer*4	io_stat			! IO status

c Externals

	External	FFP_Normal
	External	FUT_Error
	External	FUT_Normal
	External	FFP_Failure
	External	FFP_CTInit
	External	FFP_Repclose
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c Begin
	status = %loc (FFP_Normal)
c
c  Print the banner.
c
	rstatus = CUT_Register_Version (fcc_version)
	rstatus = CUT_Display_Banner ( terminal, 80,
	1                             'FIRAS Facility FFP_Spectra_Coadd')
	Write (terminal,10)
 10	Format (/)
c
c  Parse the command line.
c
	rstatus = FFP_SC_Parse ( current_gmt, data_type, input, sfile, slen,
	1                        cmdline, filexts, fnum )
c
c  Initialize the processing report.
c
	If ( fcc_report .EQ. fac_present ) Then

	   status = FFP_SC_Init_Report ( rstatus, current_gmt, input, sfile,
	1                                slen, cmdline, filexts, fnum,
	2                                archin, archout )

	   If (status .EQ. %loc(FFP_Normal)) Call Lib$Establish (FUT_Error)
	EndIf

	If ( (status .EQ. %loc(FFP_Normal)) .AND.
	1    (rstatus .EQ. %loc(FFP_Normal)) ) Then
C
C  Initialize Cobetrieve.
C
	   Call CT_Init (ct_stat)

	   If (ct_stat(1) .EQ. CTP_Normal) Then

C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
c  Processing for sky data is done with CSA functions for skymaps.
c  Processing for cal data is done with Cobetrieve reads and writes.
C
	      If ( data_type .EQ. 'SKY' ) Then
	         status = FFP_SC_Process_Sky ( input, filexts, fnum, archin,
	1                                      archout, nrec, npix )
	      Else
	         status = FFP_SC_Process_Cal ( input, filexts, fnum, archin,
	1                                      archout, nrec )
	         npix = 0
	      EndIf

	   Else
	      status = %loc(ffp_ctinit)
	      Call lib$signal (ffp_ctinit, %val(1), %val(ct_stat(1)))
	   EndIf		!  return status from CT_INIT
	EndIf			!  return status from parse and report init
c
c  Update the processing report with the number of spectra processed.
c
	If ( (status .EQ. %loc(ffp_normal))  .AND.
	1    (fcc_report .EQ. fac_present) ) Then
	   status = FFP_SC_Update_Report ( data_type, nrec, npix, input,
	1                                  filexts, fnum )
	EndIf
c
c  Print out the number of records and pixels processed.
c
	Type 20, nrec
 20	Format (/, 4x, 'Number of records reformated:     ', I5)
	If (data_type .EQ. 'SKY') Then
	   Type 30, npix
 30	   Format (4x, 'Number of pixels containing spectra:  ', I5)
	EndIf
c
c  Signal the program completion status.
c
	If ( (status .EQ. %loc(ffp_normal)) .AND.
	1    (rstatus .EQ. %loc(ffp_normal))) Then
	   Call lib$signal (ffp_normal)
	Else
	   Call lib$signal (ffp_failure)
	EndIf
c
c  Close the processing report file.
c
	If (fcc_report .EQ. fac_present) Then
	   Close (unit=fut_report_lun, iostat=io_stat)
	   If (io_stat .NE. 0) Then
	      Call lib$signal (ffp_repclose, %val(2),
	1                      fcc_report_file(1:fcc_replen), %val(io_stat))
	   EndIf
	   rstatus = fut_free_lun (fut_report_lun)
	   If (rstatus .NE. %loc(fut_normal)) Then
	      Call lib$signal (%val(rstatus))
	   EndIf
	EndIf
	End
