	Program  FAD

C==============================================================================
C                       FAD - FIRAS Apply Destriper
C                       ---------------------------
C   FAD applies the "destriper" time dependent corrections, created by FDS,
C   to the FIRAS calibrated spectra, FCF_SKY_ccss, producing the corrected
C   data set FAD_SKY_ccss.
C
C   Input:   CSDR$FIRAS_CAL:  FEX_EJV_ccss         - destriper corrections.
C            CSDR$FIRAS_IN:   FCF_SKY_ccss         - calibrated spectra.
C
C   Output:  CSDR$FIRAS_OUT:  FAD_SKY_ccss         - corrected spectra.
C            local directory: FAD_jstart_jstop_ccss.REP_rungmt  - report file.
C
C   (cc = channel = (RH,RL,LH, or LL), ss = scan mode = (SS,SF,LS, or LF)
C    jstart_jstop = data times from FCF_SKY file extension,
C    rungmt = FAD run time in GMT format.)
C
C   Algorithm:
C
C      Correction (i,f) = Sum_ufo Ai exp (-t/tau(i)) +
C                         Sum_top Bj Topj (t) +
C                         V0 + V1 ty + V2 ty^2 + V3 ty^3 + V4 ty^4
C
C      Corrected Spectra (i,f) = Spectra (i,f) - Correction (i,f)
C
C      where: i = spectrum number, f = frequency index, t is time in days
C      since cover ejection (893251118), ty = time in years since cover
C      ejection.  The A's, B's, V's and tau's are all read from the FEX_EJV
C      reference file.
C
C   Author:  Larry Paris Rosen, Hughes STX, 21 April 1993, SER 10798
C------------------------------------------------------------------------------
C   Modifications:
C      L. Rosen, February 1994.  Add Gain correction.  New reference file
C      FEX_GN contains gain correction.  SER 11688
C      L. Rosen, May 1994.  Removed all Gain correction references since
C      the destriper with gain proved faulty.
C      L. Rosen, June 1994.  New correction formula above includes several
C      taus (ufo terms) and new vibration correction quartic.
C      L. Rosen, June 1994.  Change FEX_EJ to FEX_EJV.  SER 11796.
C      L. Rosen, July 29 1994.  At frequencies where the variance is zero, the
C      spectrum should be zero.  It isn't because of FCF's general code, but
C      this caused floating point underflow warning errors in FCS and FMS.
C      The simple change is to not copy spectrum where variance is zero.
C      SER 11847
C------------------------------------------------------------------------------
	Implicit None

c Include files

	Include		'($ssdef)'
	Include		'(FUT_error)'
	Include		'(FAD_msg)'

c External

	External	fut_error

c Functions

	Integer*4	Cut_Register_Version, Cut_Display_Banner
	Integer*4	FAD_Parse_Command, FAD_Init_Report, FAD_Open_Archives
	Integer*4	FAD_Read_Reference_Data, FAD_Process_Spectra
	Integer*4	FAD_Close_Archives, FAD_Close_Report

c Local

	Integer*4	rstat                                   ! Return status
	Character*6	version
	Parameter	(version = '12.3  ')
	Character*2	chan, scan     ! FIRAS channel and scan mode (2 letter)
	Character*22	filext                 ! File extension of FCF_SKY file
	Character*22	modelext     ! FISH model name - FEX_EJV file extension
	Logical*1	report                            ! Flag to make report
	Character*45	report_name                            ! Name of report
	Character*79	cmdline(2)                  ! Command line with default
	Integer*4	current_adt(2)                        ! run time in adt
	Integer*4	lun_rpt, lun_in, lun_out                 ! unit numbers
	Dictionary	'FEX_EJV'
	Record /FEX_EJV/ fex_ejv_rec
	Character*14	arch_in                            ! input data archive
	Parameter	(arch_in='CSDR$FIRAS_IN:')
	Character*15	arch_cal                            ! reference archive
	Parameter	(arch_cal='CSDR$FIRAS_CAL:')
	Character*15	arch_out                          ! output data archive
	Parameter	(arch_out='CSDR$FIRAS_OUT:')
	Character*50	filein, fileout                            ! file names
	Character*51	file_ejv                                   ! file name
	Integer*4	access_time (2)		! time for reference data
c------------------------------------------------------------------------------
c Begin

	rstat = Cut_Register_Version (version)
	rstat = Cut_Display_Banner (6, 80, 'FIRAS  FAD_Apply_Destriper')
	rstat = %Loc (FAD_Normal)

c Call FAD_Parse_Command to get the command line qualifiers.

	rstat = FAD_Parse_Command ( chan, scan, filext, modelext, report,
	1                           report_name, cmdline, current_adt )

c Call FAD_Init_Report to open and write initial info to report file.

	If (report .AND. (rstat .EQ. %Loc (FAD_Normal))) Then

	   rstat = FAD_Init_Report ( report_name, lun_rpt, cmdline,
	1                            current_adt, version, arch_in, arch_cal,
	2                            arch_out )
	Endif

c Establish error handler (writes errors to report file and screen).

	Call Lib$Establish ( fut_error )

c  Call FAD_Open_Archives.

	If (rstat .EQ. %Loc (FAD_Normal)) Then

	   rstat = FAD_Open_Archives ( lun_in, lun_out, lun_rpt, report, chan,
	1                              scan, filext, arch_in, arch_out,
	2                              filein, fileout, access_time )
	Endif

c Call FAD_Read_Reference_Data.

	If (rstat .EQ. %Loc (FAD_Normal)) Then

	   rstat = FAD_Read_Reference_Data ( arch_cal, modelext, chan, scan,
	1                                    lun_rpt, report, fex_ejv_rec,
	2                                    file_ejv, access_time )
	Endif

c Call FAD_Process_Spectra to apply correction to all spectra.

	If (rstat .EQ. %Loc (FAD_Normal)) Then

	   rstat = FAD_Process_Spectra ( lun_in, lun_out, lun_rpt, report,
	1                                fex_ejv_rec )
	Endif

c Call FAD_Close_Archives

	If (rstat .EQ. %Loc (FAD_Normal)) Then
	   rstat = FAD_Close_Archives ( lun_in, lun_out, lun_rpt, report,
	1                               filein, fileout )
	Endif

c Close report file, signal the processing status and exit.

	If (rstat .EQ. %Loc (FAD_Normal)) Then
	   Call Lib$Signal (FAD_Normal)
	Else
	   Call Lib$Signal (FAD_Abort)
	Endif
	If (report .AND. (rstat .EQ. %Loc (FAD_Normal))) Then
	   rstat = FAD_Close_Report ( lun_rpt, report_name )
	Endif

	If (rstat .EQ. %Loc (FAD_Normal)) Then
	   Call Exit (SS$_Normal)
	Else
	   Call Exit (SS$_Abort)
	Endif
	Stop
	End
