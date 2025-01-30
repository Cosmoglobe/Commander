C******************************************************************************
C
C                                  FIP_COVAR
C
C                     FIRAS Initial Product COVARiance matrix
C
C******************************************************************************
C
C  The purpose of FIP_COVAR is to convert the calibrated covariance matrix
C  output by FCC_Calibrate_Covariances into proper units and write only the
C  frequency bins that are specified for the FIRAS initial product data
C  release.
C
C  FIP_COVAR will retrieve the covariance matrix and associated information
C  from the FCC output, for the specified channel and scan mode, convert the
C  units from Eplees2 to MJy/sr2, and write the submatrix for the specified
C  frequency range for the deliverable initial product.
C
C  Author:  Larry Paris Rosen, Hughes STX, 13 July 1993, SER 11259.
C
C  This is the main module of the FIP_COVAR program.
C
C  Here is the PDL:
C
C    Register software version number and display banner.
C    Call FIP_PARSE_COV to parse command line.
C    Call FIP_INIT_REPORT_COV to initialize report file.
C    Call Lib$Establish to establish the error handler.
C    Call CT_INIT to initialize cobetrieve.
C    Call FIP_READ_COV to open, read, and close the FCC covariance matrix file.
C    Call FIP_CONVERT_COV to convert the units.  Uses FIP_FREQUENCY_CUT to
C       determine index range given frequency range.  
C    Call FIP_WRITE_COV to open, write the converted covariance matrix and
C       associated quantities to, and close the output file.
C    Close the report file.
C    Signal completion status.
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Program  FIP_COVAR

	Implicit None

C  Include Files

	Include		'($ssdef)'
	Include		'CSDR$LIBRARY:CTPARAMS.INC'
	Include		'(fut_error)'

C  External

	External	fut_error
	External	fip_normal
	External	fip_ctinit
	External	fip_closerr
	External	fip_abort

C  Functions

	Integer*4	CUT_REGISTER_VERSION, CUT_DISPLAY_BANNER
	Integer*4	FIP_PARSE_COV
	Integer*4	FIP_INIT_REPORT_COV
	Integer*4	FIP_READ_COV
	Integer*4	FIP_CONVERT_COV
	Integer*4	FIP_WRITE_COV

C  Local Variables

	Integer*4     rstatus              ! Return status for function calls.
	Integer*2     status /0/           ! Processing status
	Integer*2     err /2/, ok /0/      ! Processing status codes
	Character*6   version              ! Version of software.
	Parameter     (version = '11.4  ')
	Character*20  filext               ! Input file name extension
	Character*2   chan                 ! Channel to process (RH,RL,LH,LL)
	Character*2   scan                 ! Scan mode to process (SS,SF,LS,LF)
	Real*4        freq_low             ! Low frequency cut-off (icm)
	Real*4        freq_high            ! High frequency cut-off (icm)
	Logical*1     report               ! Flag whether to write report
	Character*42  repfile              ! File name for report
	Character*79  cmdline (3)          ! Command line with defaults
	Integer*4     curtime(2)           ! Current system time in ADT
	Character*14  arch_in  /'CSDR$FIRAS_IN:'/    ! input data archive
	Character*15  arch_out /'CSDR$FIRAS_OUT:'/   ! output data archive
	Character*15  arch_ref /'CSDR$FIRAS_REF:'/   ! reference data archive
	Integer*4     lun_rpt              ! Logical unit number - report file
	Integer*2     ct_stat (20)         ! Cobetrieve status block
	Dictionary    'fcc_cov'            ! FCC data structure.
	Record / fcc_cov / related_data    ! Data in the first FCC record
	Integer*4     bin_total (256)      ! Number of voltage spectra
	Real*4        wtd_bin_total (256)  ! Weighted total voltage spectra
	Real*8        eff_wt (256)         ! Effective weight of each volt spec
	Real*8        rdisp (257,256)      ! Total real spectra dispersion
	Real*8        idisp (257,256)      ! Total imag spectra dispersion
	Real*8        temp (8,256)         ! Temp and glitch rate disp. total
	Real*8        fcc_covar (522,522)  ! FCC Covariance Matrix
	Integer*2     index_low            ! Low frequency cut-off array index
	Integer*2     index_high           ! High frequency cut-off array index
	Real*4        delta_nu             ! Optical frequency interval, GHz.
	Real*4        nu_zero              ! Optical frequency of initial data
                                           !    point, in GHz.

C Begin Code

	rstatus = CUT_REGISTER_VERSION ( version )
	rstatus = CUT_DISPLAY_BANNER (6, 80, 'FIRAS  FIP_COVAR')

C Call FIP_PARSE_COV to parse the command line.

	rstatus = FIP_PARSE_COV ( filext, chan, scan, freq_low, freq_high,
	1                         report, repfile, cmdline, curtime )

	If ( rstatus .NE. %Loc (fip_normal) ) status = err

C Call FIP_INIT_REPORT_COV to open and initialize report file.

	If ( ( status .EQ. ok ) .AND. report ) Then

	   rstatus = FIP_INIT_REPORT_COV ( filext, chan, scan, freq_low,
	1                                  freq_high, repfile, cmdline,
	1                                  curtime, arch_in, arch_out,
	1                                  arch_ref, version, lun_rpt )

	   If ( rstatus .NE. %Loc (fip_normal) ) status = err
	Endif

C Call Lib$Establish to establish error handler.

	Call Lib$Establish ( fut_error )

C Call CT_INIT to initialize Cobetrieve.

	If ( status .EQ. ok ) Then

	   Call CT_INIT ( ct_stat )

	   If (ct_stat(1) .NE. ctp_normal) Then
	      Call Lib$Signal (fip_ctinit, %val(1), %val(ct_stat(1)))
	      status = err
	   Endif
	Endif

C Call FIP_READ_COV to open, read, unpack, and close the FCC output
C covariance matrix.

	If ( status .EQ. ok ) Then

	   rstatus = FIP_READ_COV ( arch_in, lun_rpt, report, filext, chan,
	1                           scan, related_data, bin_total,
	1                           wtd_bin_total, eff_wt, rdisp, idisp,
	1                           temp, fcc_covar )

	   If ( rstatus .NE. %Loc (fip_normal) ) status = err
	Endif

C Call FIP_CONVERT_COV to convert the data.

	If ( status .EQ. ok ) Then

	   rstatus = FIP_CONVERT_COV ( lun_rpt, report, chan, scan,
	1                              related_data, bin_total, wtd_bin_total,
	1                              eff_wt, rdisp, idisp, temp, fcc_covar,
	1                              freq_low, freq_high, index_low,
	1                              index_high, arch_ref, delta_nu, nu_zero)

	   If ( rstatus .NE. %Loc (fip_normal) ) status = err
	Endif

C Call FIP_WRITE_COV to open, write the converted covariance matrix and
C associated quantities to, and close the output file.

	If ( status .EQ. ok ) Then

	   rstatus = FIP_WRITE_COV ( lun_rpt, report, related_data,
	1                            bin_total, wtd_bin_total, eff_wt,
	1                            rdisp, idisp, temp, fcc_covar,
	1                            index_low, index_high, arch_out,
	1                            chan, scan, filext, delta_nu, nu_zero )

	   If ( rstatus .NE. %Loc (fip_normal) ) status = err
	Endif

C Signal the processing status.

	If (status .EQ. ok) Then
	   Call Lib$Signal (fip_normal)
	Else
	   Call Lib$Signal (fip_abort)
	Endif

C Close report

	If (report) Then
	   Close (Unit=lun_rpt, Iostat=rstatus)
	   fut_report_lun = 0
	   If (rstatus .NE. 0) Then
	      status = err
	      Call Lib$Signal (fip_closerr, %val(2), repfile,
	1                      %val(rstatus))
	   Endif
	Endif

C Exit

	If (status .EQ. ok) Then
	   Call Exit (ss$_normal)
	Else
	   Call Exit (ss$_abort)
	Endif
	Stop
	End
