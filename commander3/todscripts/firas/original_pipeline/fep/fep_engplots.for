	Program FEP_EngPlots

C-------------------------------------------------------------------------------
C    PURPOSE: Collect housekeeping data and plot any set of
C             housekeeping data in engineering units.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer,  STX,  March 5, 1987
C
C    INVOCATION: ENGPLOTS
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS: None
C
C    SUBROUTINES CALLED:
C	STR$UpCase
C	FEP_Init- - - - - - -<FEP_Get_Command
C	FUT_Timerange
C	FEP_Fix_Counts_to_Ohms
C	                    _/FEP_Display_Main
C	FEP_Display_Menu- -' \FEP_Display_Plot_Menu
C	FEP_Match_Fields      /fep_get_grt_attr
C	                     /FEP_Get_GRT_Conv- -*            /FEP_Decode_DBname
C	                    //FEP_Counts_To_Ohms-<FEP_Get_Calres
C	FEP_Convert_To_Eng-<<-FEP_GRT_Lookup- - -<FEP_Invert
C	                    \\GRT_Correction
C	                     \*FUT_Field_Attributes
C                             \FEP_GET_CURVE
C
C	FEP_Plot- - - - - - -<FEP_Convolve
C	FEP_Report
C	FEP_Exit
C	Time_GT
C	CT_Init
C	CT_Open_Arcv
C	CT_Read_Arcv
C	CT_Close_Arcv
C	CT_Binary_To_GMT
C	LIB$Movc3
C	LIB$Signal
C	SYS$GetTim
C
C    COMMON VARIABLES USED:
C	NRECS			I*4		Number of HSKP records.
C	HSKP_DATA		BYTE(248,3000)	Hskp data (248 bytes = 1 mjf)
C
C    INCLUDE FILES:
C	FEP_Invoc
C	FEP_Menu
C	FEP_Data
C	FUT_Params
C	FUT_Error
C	$SSDef
C	$LNMDef
C	CTUser.Inc
C
C    Software Maintainers, CAUTION:  Hardcoded numbers used in tests of
C    menu_element() in modules FEP_PLOT and FEP_CONVERT_TO_ENG.  Look for
C    this comment block in those modules:
C	c
C	c  BEWARE--HARDCODED MENU ELEMENT NUMBERS!  REVISE WHEN MENU CHANGES!
C	c
C-------------------------------------------------------------------------------
C
C Changes:
C
C	Make room in the data holding buffer when the data limit is
C	exceeded in monitor mode.  R. Kummerer, August 20, 1987.
C
C	Move wait for next major frame (when monitoring) to FEP_PLOT.
C	November 4, 1987, R. Kummerer.
C
C	SPR 1692, Quit directly to timerange prompt. December 2, 1987,
C	R. Kummerer.
C
C	SPR 1820, Use FAC_MAX_MENU_NUM to maintain plot menu pages.
C	R. Kummerer, December 11, 1987.
C
C	Convert from plotting with Vecplt to PLT.  F. Shuman, 1988 Apr 8.
C
C	Changed status to FEP_Normal for loop in Build 4.1 problem.
C						   Shirley M. Read, 1988 Aug 22.
C
C	Label each curve just inside the left and right plot boundaries with
C	its numerical order in the field list.  This module and FEP_PLOT.
C	(SPR 2580)  F. Shuman, STX, 1988 Nov 1.
C
C	Split up command line invocation into 3 strings to write to report file
C	to prevent FORTRAN 'Output Statement Overflows Record' error.  Module
C	FEP_REPORT.  (SPR 2605)  F. Shuman, STX, 1988 Nov 1.
C
C	Reverse the order of arguments in LineStyle command sent to PLT, due to
C	syntax change in newly delivered PLT.  Old:  LS <vector#> <LStyle#>;
C	new:  LS <LStyle#> ON <vector#>.  Module FEP_PLOT.  (SPR 2650)
C	F. Shuman, STX, 1988 Nov 1.
C
C	Collect dwell_mode(*) returned by FEP_Convert_To_Eng in dwell_3(*,*),
C	the second index going from 1 to 3 to indicate separately the presence
C	of dwell data on A side, B side, and A or B side.  Using this array, the
C	presence of GRTA dwell will not blank out GRTB or non-GRT fields.  This
C	module and FEP_PLOT.  (SPR 2651) F. Shuman, STX, 1988 Nov 1.
C
C	Increase the index limits in declarations of plot_buff() to allow the
C	full 9 fields to be plotted.  This module and FEP_PLOT.  (SPR 2673)
C	F. Shuman, STX, 1988 Nov 1.
C
C	SPR 2377, Provide masking capability for digital status words.
C	Fred Shuman,  1989 Feb 20.
C
C	SPR 3662.  Engplots labels are unreadable when there are more fields
C	than will fit in the title of a QMS (e.g., TALARIS) plot.  Remedy:
C	reduce the number of selectable fields by:
C	o Reducing the dimension of fcc_fields() and fcc_fields_len() in
C	  FEP_Invoc.txt from 10 to 4
C	o Reducing the limit on fcc_fields_count in FEP_Get_Command from 10 to 4
C	o Reducing the value of the parameter 'maxfields' in FEP_Match_Fields
C	  from 9 to 4.
C	F. Shuman, STX, 1989 May 12.
C
C	SPR 3563.  Engplots/Monitor aborts when /Jstart & /Jstop are on the
C	command line.  Remedy:  revise the logic for getting the data selection
C	timerange and checking for the 'realtime' condition.
C	F. Shuman, STX, 1989 May 24.
C
C	Checks the Real_time flag.  If we are not running in real time, calls
C	the routine FEP_Fix_Counts_to_Ohms to "pre-compute" the counts to ohms
C	conversion coefficients.   Don Stevens-Rayburn, ARC, 13-Jul-1989.
C
C	Previous change modified so that it pre-computes the counts to ohms
C	conversion coefficients only if we are:
C	o  Not running in real time, and
C	o  The user did not specify "/NOAVERAGE_CALIBRATION" on the command
C	   line.
C	Don Stevens-Rayburn, ARC, 07-Aug-1989.
C
C	SER 3306, Call CUT_Register_Version and CUT_Display_Banner to generate
C	a 'banner'.  F. Shuman, STX, 1989 Aug 30.
C
C	SPR 4463.  Engplots aborts:  divide by zero in FEP_Fix_Counts_to_Ohms.
C	When all values are equal, their sample std dev will be 0.  Test for
C	non-outliers was "x .gt. avg-3s .and. x .lt. avg+3s", so when s=0,
C	nothing passed.  Changed .gt. and .lt. to .ge. and .le.  Also put in a
C	check for n=0 or 1, so as not to divide by 0 when computing avg and s.
C	Normally, with lots of data points, we needn't worry about s=0, but when
C	the chosen timerange has only two points, it becomes pretty likely.
C	Qfix 422. Fred Shuman,  STX,  1989 Aug 31.
C
C       SPR 4852, pass variable 'Fr_lun' into routine fep_FIX_COUNTS_TO_OHMS,
C       in order to handle the dewell data right in that routine.
C       Qfix #, Harte Wang, STX, 1989, Nov. 3
C
C       QFIX 636 / SPR 4958.  Engplots issues YMIN = YMAX error from PLT.
C       Rescale values sent to PLT were increased in precision from E14.4 to
C       E16.6.  Module FEP_Plot.  Fred Shuman,  STX,  1989 Nov 6.
C
C       SPR 5109, When combining GRTs, convert the low current GRT regardless
C	of the success of the high current GRT and then combine temperatures
C	according to the temperature regime AND the success of either
C	conversions.  Version 4.4.05, R. Kummerer, STX, Nov 20, 1989.
C
C       SPR 5068, Avoid plotting data points when telemetry quality is bad
C	                 Version 4.4.05, H. Wang, STX, DEC. 2, 1989.
C
C       SPR 5110, Cannot specifiy a mask for plotting statuses in batch
C	                 Version 4.4.05, H. Wang, STX, DEC. 2, 1989.
C
C       SPR 3959,5078, Hard code FEP menu pages.
C	                 Version 5.0, H. Wang, STX, DEC. 20, 1989.
C       SPR 5472, Unflagged bad data in houskeeping
C	                 Version 5.2, H. Wang, STX, Jan. 11, 1990.
C       SPR 5713, Since the NOTCH FILTER data are same for both
C                  first and second major frames, it is not necessary to
C                  copy notch filter data from first major frame to
C                  second major frame
C	                 Version 5.5, H. Wang, STX, Jan. 29, 1990.
C
C       SPR 3921,  add the capabilities of reading the reference data
C                  from reference archive.
C	           Version 5.7, H. Wang, STX, Feb. 22, 1990.
C
C       SPR 5713,  FEP Will plot bad telemetry data when it plots NOTCH
C                  FILTER FIELDS.
C	           Version 5.7, H. Wang, STX, Feb. 22, 1990.
C
C       SPR 6253,  Xcal_Temp_S6_b GRT plot is Scaled to 6 x E04 for
C                  Xcal motion.
C	           Version 5.7, H. Wang, STX, Feb. 22, 1990.
C
C       SER 5725,  Add capability to access a plt command from CLD
C	           Version 5.8, H. Wang, STX, Mar. 2, 1990.
C
C       SPR 6884, Correct the hot spot engplots label on the "Y" axis
C                 Version 6.3, H. Wang, STX, JUNE 14, 1990
C
C	SPR 4171,  Standardize report file name.  Changes in FEP.CLD and
C	           FEP_GET_COMMAND.FOR.  L.P. Rosen, STX, August 6, 1990.
C
C       SPR 7369, Incorrect labelling of the CAL. resistors.
C                 Version 6.8, H. Wang, STX, Aug. 30, 1990
C
C       SPR 7514,7673, Rename the calibration resistors.
C                 Version 7.1, H. Wang, STX, Nov. 13, 1990
C
C	SPR 2690, Change the time axis origin from the first HKP rec time to
C		the user-specified start time, so the plots will be easier to
C		compare.  F. Shuman,  STX,  1990 Nov 12.
C
C       SPR 7769, Rename the 'LVDT_STATUS' to 'MTM_POSITION'
C                 Version 7.3, H. Wang, STX, Dec. 5, 1990
C
C       SPR 9079, Invoke FUT_COMBINE_HI LO to combine Hi, Lo current readings
C		  for GRT temperatures.  Change made to FEP_Convert_to_Eng.
C		  STX, H.WANG, 1991 Sep. 27.
C
C       SPR 9290, FEP using new GRT ref. files in wrong order
C		  Hughes STX, H.WANG, 1991 Nov. 26.
C
C       SPR 9368, FEP need free logical unit number.
C		  Hughes STX, H.WANG, 1992 Jan. 7.
C-------------------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Params)'
	Include		'(FUT_Error)'
	Include		'(FEP_Invoc)'
	Include		'(FEP_Menu)'
	Include		'(FEP_Firnames)'
	Include		'(FEP_Data)'
	Include		'($SSDef)'
	Include		'(CCT_GET_CONFIG)'
	Include		'CT$Library:CTUser.Inc'
C
	Dictionary		'NFS_HKP'
	Record /NFS_HKP/	hskp_rec
	Record /Config_Status/ Stat(3)
C
	Integer		*4	CUT_Register_Version
	Integer		*4	CUT_Display_Banner
	Integer		*4	num_vol/80/
	Integer		*4	lun_out/6/
	Character	*6	version
	Parameter		(version='9.1')
C
	Character       *32     Diset(3)
	Character       *32     Dset(2)
	Data Diset(1)  /'Csdr$Firas_Ref:FDB_GRTCAL'/
	Data Diset(2)  /'Csdr$Firas_Ref:FDB_GRTRICH'/
	Data Diset(3)  /'Csdr$Firas_Ref:FDB_LIMCUVKY'/
	Data Dset(1)  /'Csdr$Firas_Ref:Fex_GRTrawwt'/
	Data Dset(2)  /'Csdr$Firas_Ref:Fex_GRTtrans'/
	Integer         *4      Nidset/3/
	Integer         *4      Ndset/2/
	Integer         *4      Config_lun(3)
	Integer         *4      Con_lun(2)
	Integer         *4      Index(3)
	Integer         *4      Indx(2)
	Integer         *4      Ncache/1/
	Integer         *4      Ref_Count
	Integer         *4      Size(3)/20,7120,128/
	Integer         *4      Size1(2)/256,384/
	Character       *5      Access_mode/'KEYED'/
	Character       *1      Access_mode1/' '/

C
	Integer         *4      CCT_OPEN_CONFIG
	Integer         *4      CCT_CLOSE_CONFIG
C
	Integer		*4	eng_lun
	Integer		*4	ios
	Character	*64	infile
	Character       *64     telm_field
	Integer         *2      telm_length
	Integer         *2      telm_offset
	Byte                    telm_mode(fac_max_hskp_recs)
	Integer		*4	fr_lun
	Integer		*2	ct_lun
	Integer		*2	ct_status(20)
	Logical		*1	first_time
	Character	*4	ans
	Logical		*1	ans_ok
	Integer		*4	i
	Integer		*4	j
	Integer		*4	k
	Integer		*4	status
	Integer		*4	exit_status

	Integer		*4	start_time(2)
	Integer		*4	end_time(2)
	Integer		*4	current_time(2)
	Character	*14	starting_time
	Character	*14	ending_time
	Logical		*1	another_time
	Logical		*1	real_time
	Character	*14	tag_first
	Character	*14	tag_last
	Character	*14	Gmt_time
	Logical		*1	display_menu
	Logical		*1	monitor_display		! Display the menu once
							!   when monitoring
	Real		*4	timeaxis(fac_max_hskp_recs)
	Real		*4	timeorigin
	Integer		*4	iy
	Integer		*4	refyr
	Integer		*4	res
	Integer		*4	iyr
	Integer		*4	iday
	Integer		*4	ihr
	Integer		*4	imin
	Integer		*4	msec

	Integer		*2	two_hskp_frames

	Integer		*2	menu_page(fac_max_menu_num)
	Integer		*2	menu_element(fac_max_menu_num)
	Logical		*1	convert_flag(fac_max_menu_num)
	Integer		*2	digital_mask(fac_max_menu_num)
	Integer		*4	nfields
	Integer		*2	length
	Integer		*4	lth
	Integer		*4	off
	Integer		*4	dbuf
	Real		*4	conv_buff(fac_max_hskp_recs)
	Real		*4	plot_buff(fac_max_hskp_recs, 9)
	Byte			dwell_mode(fac_max_hskp_recs)
	Byte			dwell_3(fac_max_hskp_recs, 3)
	Real		*4	tjump
	Integer		*4	ngaps
	Integer		*4	gaps(fac_max_hskp_recs)
	Byte			dl1
	Integer		*2	dl2
	Integer		*4	plot_group
	Integer		*4	plot_limit

	Integer		*4	ljstart
	Integer		*4	ljstop
	Character	*14	jstart_default
	Character	*14	jstop_default
	Character	*14	trans
	Integer         *4      Time_Lt
	Integer		*4	STR$UpCase
	Integer		*4	FEP_Init
	Integer		*4	FEP_Get_HSKP_Fields
	Integer		*4	FUT_field_attributes
	Integer		*4	FEP_Fix_Counts_to_Ohms
	Integer		*4	FEP_Display_Menu
	Integer		*4	FEP_Match_Fields
	Integer		*4	FEP_Convert_To_Eng
	Integer		*4	FEP_Plot
	Integer		*4	FEP_Report
	Integer		*4	FEP_Exit
	Integer		*4	FUT_Timerange
	Integer         *4      Npoints
	Logical		*2	Time_GT
	Integer		*4	Lib$Get_Lun
	Integer		*4	Lib$Free_Lun
	Real		*4	CO_Coeffs ( 3, 4 )

	Parameter	gap_thresh = 5. * fac_minute

	Integer		*4	CT_Connect_Read
	External		CT_Connect_Read

	External		FUT_Error
	External                FUT_Normal
	External		FEP_Normal
	External		FEP_CTInit
	External		FEP_CTOpenHskp
	External                FEP_Clsconfigerr
	External                FEP_Opnconfigerr
	External		FEP_CTReadHskp
	External		FEP_CTCloseHskp
	External		FEP_ExitMenu
	External		FEP_InvTim
	External		FEP_TruncData
	External		FEP_NoData
	External		FEP_InvDataLen
	External		FEP_MonRealTime

	Data jstart_default/'85001000000000'/
	Data jstop_default/'99365235959999'/

c
c Initialize.
c
	status = CUT_Register_Version(version)
	status = CUT_Display_Banner(lun_out, num_vol,
	1                           'FIRAS Facility FEP_Engplots')
	Write(lun_out,61)
61	Format(//)

	another_time = .True.
	first_time = .True.
	display_menu = .True.
	monitor_display = .True.
	two_hskp_frames = 2*(fac_hskp_frame_size + fac_notch_size)

	nfields = 0
	Do i=1,fac_max_menu_num
	   menu_page(i) = 0
	   menu_element(i) = 0
	   convert_flag(i) = .False.
	   digital_mask(i) = -1
	End Do

	Call LIB$Establish ( FUT_Error )

	status = FEP_Init ( fr_lun )

	If (status .Eq. %Loc(FEP_Normal)) Then

	   Call CT_Init ( ct_status )

	   If (ct_status(1) .Eq. CTP_Normal) Then
c
c
c
c Retrieve the housekeeping data.
c
	      Do While (another_time .And. status .Eq. %Loc(FEP_Normal))

	         If (first_time) Then

c
c	                Get the data selection timerange.
c
	            If (fcc_interactive .Eq. fac_present) Then

	               If (fcc_jstart_select .Eq. fac_not_present .Or.
	2                  fcc_jstop_select  .Eq. fac_not_present) Then

	                  status = FUT_Timerange ( starting_time, start_time,
	2                                          ending_time, end_time )
	               Else
	                  start_time(1) = fcc_jstart(1)
	                  start_time(2) = fcc_jstart(2)
	                  end_time(1) = fcc_jstop(1)
	                  end_time(2) = fcc_jstop(2)
	               End If

	               Call SYS$GetTim(current_time)

	               If (Time_GT(current_time,end_time)) Then
	                  real_time = .False.
	               Else
	                  real_time = .True.
	               End If

	            Else

	               start_time(1) = fcc_jstart(1)
	               start_time(2) = fcc_jstart(2)
	               end_time(1) = fcc_jstop(1)
	               end_time(2) = fcc_jstop(2)

	            End If

	            If (fcc_plot_dev_select .Eq. fac_not_present) Then
	               If (fcc_interactive .Eq. fac_present) Then
	                  Type 5
5	                  Format (' Select a PLT device type [/VT]: ', $)
	                  Accept 7, fcc_plot_device
7	                  Format (a)
	                  If (fcc_plot_device .Eq. ' ') Then
	                     fcc_plot_device = '/VT'
	                  End If
	               Else
	                  fcc_plot_device = '/QMS'
	               End If
	            End If

	            Call CT_Binary_To_GMT(start_time,starting_time)
	            If (real_time) Then
	               Call CT_Binary_To_GMT(current_time,ending_time)
	            Else
	               Call CT_Binary_To_GMT(end_time,ending_time)
	            End If

	            fcc_time_range = starting_time // ';' //
	2                             ending_time // ';'
	            first_time = .False.

c
c	                Check that ENGPLOTS is running real-time
c	                when monitoring data.
c
	            If (fcc_monitor .Eq. fac_present) Then
	               If (.Not. real_time) Then
	                  status = %Loc(FEP_MonRealTime)
	               End If
	            Else
	               real_time = .False.
	            End If

	         End If

c
c Open the Reference archive.
c
	         Status = CCT_OPEN_CONFIG(Start_Time,End_time,Nidset,
	1	          diset,size,access_mode,ncache,config_lun,
	2	          index,stat,ref_count)
	         If (.not. status) Then
	           Call Lib$Signal(Fep_Opnconfigerr,%val(1),%val(status))
	           Call Exit(SS$_Abort)
	         End If
	         Status = CCT_OPEN_CONFIG(Start_Time,End_time,Ndset,
	1	          dset,size1,access_mode1,ncache,con_lun,
	2	          indx,stat,ref_count)
	         If (.not. status) Then
	           Call Lib$Signal(Fep_Opnconfigerr,%val(1),%val(status))
	           Call Exit(SS$_Abort)
	         End If
c
c Open the housekeeping archive.
c
	         infile = 'CSDR$FIRAS_RAW:NFS_HKP/' // fcc_time_range
	         status = Lib$Get_Lun(eng_lun)
	         ct_lun = eng_lun
	         Open ( UNIT = eng_lun,
	1               FILE = infile,
	2               STATUS = 'old', IOSTAT = ios,
	3               USEROPEN = CT_Connect_Read )

	         If ((ct_status(1) .Eq. CTP_Normal .Or.
	2             (real_time .And.
	3                ct_status(1) .Eq. CTP_EndofFile)) .And.
	4	     status .Ne. %Loc(FEP_MonRealTime)) Then

c
c Read a record from the archive.
c
	             If (.Not. real_time) Then
	                Type 10
10	                Format (/, ' Fetching housekeeping data.')
	             Else
	                Type 15
15	                Format (' Fetching housekeeping data.')
	             End If

	             nrecs = 0

	             Call CT_Read_Arcv ( , ct_lun, hskp_rec, ct_status )

	             Do While (ct_status(1) .Eq. CTP_Normal .And.
	2                            nrecs .lt. fac_max_hskp_recs-1)

c
c Get first major frame from housekeeping record.
c
	                nrecs = nrecs + 1
	                Call LIB$Movc3 ( fac_hskp_frame_size,
	2                       hskp_rec.frame(1), hskp_data(1,nrecs) )
	                Call LIB$Movc3 ( fac_notch_size,
	2                       hskp_rec.mj_frm(1),
	3                       hskp_data(fac_hskp_frame_size+1,nrecs))
	                timetags(nrecs).time(1)=hskp_rec.Ct_Head.time(1)
	                timetags(nrecs).time(2)=hskp_rec.Ct_Head.time(2)
C	                timetags(nrecs) = hskp_rec.ct_head.gmt
c
c Convert the timetags to real numbers in seconds.  Note that the last 5 digits
c  can be taken together as milliseconds.  "refyr" is the yr field of the first
c  GMT, used as a reference year; "res" is the residue mod 4 of (refyr-1), so
c  chosen that, when added to the expression inside the INT, it will trigger a
c  leap day at the correct year-end, and 89001 minus 88366 will be 1 day.

	                Read (hskp_rec.ct_head.gmt,20) iyr, iday, ihr, imin,msec
20	                Format (i2,i3,i2,i2,i5)
	                If (nrecs .Eq. 1) Then
	                   refyr = iyr
	                   res = refyr-1 - 4*INT( (refyr-1)/4 )
	                End If
	                iy = iyr - refyr
	                timeaxis(nrecs) = (INT((iy*1461+res)/4) + iday)*fac_day
	2                        + ihr*fac_hour + imin*fac_minute + msec*1.e-3

c
c Get second major frame from housekeeping record.
c
	                nrecs = nrecs + 1
	                Call LIB$Movc3 ( fac_hskp_frame_size,
	2                       hskp_rec.frame(2), hskp_data(1,nrecs) )
	                Call LIB$Movc3 ( fac_notch_size,
	2                       hskp_rec.mj_frm(2),
	3                       hskp_data(fac_hskp_frame_size+1,nrecs))
C	                Call LIB$Movc3 ( 10,
C	2                       hskp_rec.mj_frm(1),
C	3                       hskp_data(fac_hskp_frame_size+1,nrecs))
	                Call Ct_GMT_To_Binary(Hskp_rec.hskp_tail.gmt_mjf2,
	1	                timetags(nrecs).time)
C	                timetags(nrecs) = hskp_rec.hskp_tail.gmt_mjf2
c
c Convert the timetags to real numbers in seconds.  Note that the last 5 digits
c  can be taken together as milliseconds.

	                Read (hskp_rec.hskp_tail.gmt_mjf2,20) iyr, iday, ihr,
	2                     imin, msec
	                iy = iyr - refyr
	                timeaxis(nrecs) = (INT((iy*1461+res)/4) + iday)*fac_day
	2                        + ihr*fac_hour + imin*fac_minute + msec*1.e-3

	                Call CT_Read_Arcv ( , ct_lun, hskp_rec,
	2                                       ct_status )

	             End Do

	             If (Time_lt(End_time, Timetags(nrecs).time)) Then
	                 nrecs = nrecs - 1
	             End If
c
c Reset the timeaxis reference (origin) to the user-specified start time.
c
	             Read (starting_time,20) iyr, iday, ihr, imin, msec
	             iy = iyr - refyr
	             timeorigin = (INT((iy*1461+res)/4) + iday)*fac_day
	2                        + ihr*fac_hour + imin*fac_minute + msec*1.e-3
	             Do i=1,nrecs
	                If (timeaxis(i) .Ne. fac_no_data) Then
	                   timeaxis(i) = timeaxis(i) - timeorigin
	                End If
	             End Do

	             tag_first = starting_time
c
c
c Display plots for the currently selected timerange.
c
	             If ((ct_status(1) .Eq. CTP_Normal .Or.
	2                   ct_status(1) .Eq. CTP_EndOfFile) .And.
	3                nrecs .Ne. 0) Then

	                If (nrecs .Eq. fac_max_hskp_recs) Then
	                   Call LIB$Signal ( FEP_TruncData )
	                   real_time = .False.
	                End If
c
c Look for gaps, expanding timeaxis() for each one found.
c     No need to expand timetags(), since we're done using it!
c
	                ngaps = 0

	                Do i = 2, nrecs
	                   tjump = timeaxis(i) - timeaxis(i-1)
	                   If (tjump .Gt. gap_thresh) Then
	                      ngaps = ngaps + 1
	                      gaps(ngaps) = i
	                   End If
	                End Do
	                gaps(ngaps+1) = nrecs + 1

	                Do i = ngaps, 1, -1
	                   Do j = gaps(i+1)-1, gaps(i), -1
	                      timeaxis(i+j) = timeaxis(j)
	                   End Do
	                   timeaxis(i+gaps(i)-1) = fac_no_data
	                End Do
	                npoints = nrecs + ngaps

	                Call ct_binary_to_Gmt(timetags(nrecs).time,gmt_time)
	                tag_last  = gmt_time
c
c If the user either specified or defaulted /AVERAGE_CAL_COUNTS on the command
c line and we are not running in real time, precompute the counts to ohms
c conversion coefficients.
c
	                If ( fcc_average_cal_counts .Eq. fac_present .And.
	2                    .Not. Real_time ) Then
	                   status = FEP_Fix_Counts_To_Ohms (fr_lun,
	1                         config_lun(1), co_coeffs )
	                   If ( .Not. status ) Then
	                      Call LIB$Signal ( %Val ( status ) )
	                   End If
	                End If

	                Do While (display_menu)
c
c Select fields from a menu or by matching /FIELDS
c with menu fields.
c
	                   If (fcc_interactive .eq. fac_present) Then
c
c	                        Always display a menu interactively; display
c	                        the first time only when monitoring.
c
	                      If (monitor_display) Then

	                         status = FEP_Display_Menu ( menu_page,
	2                                                    menu_element,
	3                                                    convert_flag,
	4                                                    digital_mask,
	5                                                    nfields )
c
c	                            Turn off menu display when quitting.
c
	                         If (status .Ne. %Loc(FEP_ExitMenu)) Then
	                            display_menu = .False.
	                         End If
c
c	                            Turn off menu display after selection
c	                            when monitoring.
c
	                         If (fcc_monitor .Eq. fac_present) Then
	                            display_menu = .False.
	                            monitor_display = .False.
	                         End If

	                      Else
	                         display_menu = .False.
	                      End If

	                   Else
c
c	                        Get fields from /FIELDS.
c
	                      status = FEP_Match_Fields ( menu_page,
	2                                                 menu_element,
	3                                                 convert_flag,
	4                                                 nfields )
	                      display_menu = .False.

	                   End If

c
c Convert to engineering units and make plots.
c
	                   If (status .Eq. %Loc(FEP_Normal) .Or.
	2                      status .Eq. %Loc(FEP_ExitMenu)) Then

c
c Collect the data to be plotted.
c
	                      Do j = 1, fac_max_hskp_recs
	                         Do k = 1, 3
	                            dwell_3(j,k) = 0
	                         End Do
	                      End Do

	                      Do i = 1, nfields

	                         status = FEP_Convert_To_Eng (
	1	                             fr_lun,Config_lun,con_lun,
	2                                    menu_page(i),
	3                                    menu_element(i),
	4                                    convert_flag(i),
	5                                    digital_mask(i),
	6                                    conv_buff,
	7                                    dwell_mode,
	8                                    co_coeffs  )

	                         Do j = 1, fac_max_hskp_recs
	                            plot_buff(j,i) = conv_buff(j)
	                         End Do

	                         If (menu_page(i) .Eq. GRTA_page .Or.
	2                            (menu_page(i) .Eq. COMB_page .And.
	3                             menu_element(i) .Le.
	4                               Comb_menu_num/2) ) Then
c
c	           ...for A-side GRTs, flag A-side dwell and overall dwell.
c
	                            Do j = 1, fac_max_hskp_recs
	                               If (dwell_mode(j) .Eq. 1) Then
	                                  dwell_3(j,1) = 1
	                                  dwell_3(j,3) = 1
	                               End If
	                            End Do
	                         Else If (menu_page(i) .Eq. GRTB_page .Or.
	2                            (menu_page(i) .Eq. COMB_page .And.
	3                             menu_element(i) .Gt.
	4                               Comb_menu_num/2)) Then
c
c	           ...for B-side GRTs, flag B-side dwell and overall dwell.
c
	                            Do j = 1, fac_max_hskp_recs
	                               If (dwell_mode(j) .Eq. 1) Then
	                                  dwell_3(j,2) = 1
	                                  dwell_3(j,3) = 1
	                               End If
	                            End Do
	                         End If

	                      End Do

c
c Use the gaps found in timeaxis() to expand dwell_3() and plot_buff().
c
	                      Do i = ngaps, 1, -1
	                         Do j = gaps(i+1)-1, gaps(i), -1
	                            dwell_3(i+j, 1) = dwell_3(j,1)
	                            dwell_3(i+j, 2) = dwell_3(j,2)
	                            Do k = 1,nfields
	                               plot_buff(i+j, k) = plot_buff(j,k)
	                            End Do
	                         End Do
	                         dwell_3(i+gaps(i)-1, 1) = fac_no_data
	                         dwell_3(i+gaps(i)-1, 2) = fac_no_data
	                         Do k = 1,nfields
	                            plot_buff(i+gaps(i)-1, k) = fac_no_data
	                         End Do
	                      End Do
c
c Find the telemetry quality bad record
c
	                     If (status .eq. %loc(Fep_Normal)) Then
	                      telm_field = 'FEP_HKP.CURR_MENU.TELEMETRY_QUALITY'
	                      Status=FUT_FIELD_Attributes(Fr_lun,Telm_field,
	2                                              telm_length,telm_offset)
	                     End If
	                      If (status .eq. %loc(fut_normal)) Then
	                       Do k=1, nrecs
	                         telm_mode(k) = 0
	                         If (hskp_data(telm_offset+1,k) .ne. 0) Then
	                          telm_mode(k) = 1
	                         End If
	                       End Do
	                      End If
	                      If (status .Eq. %Loc(Fut_Normal)) Then
	                         status = FEP_Plot (
	2                                    plot_buff,telm_mode,
	3                                    dwell_3, timeaxis,
	4                                    tag_first, tag_last,
	5                                    npoints, nfields,
	6                                    menu_page, menu_element,
	7                                    convert_flag )
	                      End If

	                      If (nfields .Ne. 0 .And.
	2                         status .Eq. %Loc(FEP_Normal)) Then

	                         status = FEP_Report ( menu_page,
	2                                           menu_element,
	3                                           convert_flag,
	4                                           nfields )

	                      End If

	                   End If

	                End Do

	             Else

	                If (nrecs .Eq. 0) Then
	                   Call LIB$Signal ( FEP_NoData )
	                Else
	                   status = %Loc(FEP_CTReadHskp)
	                   Call LIB$Signal(FEP_CTReadHskp,%Val(1),
	2                               %Val(ct_status(1)))
	                End If

	             End If

c
c Close the housekeeping & reference data archive.
c
	             Call CT_Close_Arcv ( , ct_lun, ct_status )

	             If (ct_status(1) .Ne. CTP_Normal) Then
	                Call LIB$Signal(FEP_CTCloseHskp,%Val(1),
	2                               %Val(ct_status(1)))
	             End If

	             Status = CCT_Close_Config(nidset,config_lun,index)
	             If (.not. status) Then
	               Call Lib$Signal(Fep_Clsconfigerr,%val(1),%val(status))
	               Call Exit(SS$_Abort)
	             End If

c
c Prepare for the next set of plots.
c
	             If (real_time) Then

	                Call SYS$GetTim(current_time)
	                If (Time_GT(current_time,end_time)) Then
	                   real_time = .False.
	                Else
	                   real_time = .True.
	                End If
	                Call CT_Binary_To_GMT(current_time,ending_time)
	                fcc_time_range = starting_time // ';' //
	2                                ending_time // ';'
	                display_menu = .True.

	             Else

	                If (fcc_interactive .Eq. fac_present .And.
	2                  (fcc_jstart_select .Eq. fac_not_present .And.
	3                   fcc_jstop_select .Eq. fac_not_present)) Then
	                   ans_ok = .False.
	                   Do While (.Not. ans_ok)
	                      Type 40
40	                      Format (' Do you want to select '
	2                              'another timerange? (Y/[N]): ', $)
	                      Accept 45, ans
45	                      Format (a)
	                      status = STR$UpCase(ans,ans)
	                      If (ans(1:1) .Eq. 'Y') Then
	                         ans_ok = .True.
        	                 status = Lib$Free_Lun(eng_lun)
        	                 status = Lib$Free_Lun(config_lun(1))
        	                 status = Lib$Free_Lun(config_lun(2))
        	                 status = Lib$Free_Lun(config_lun(3))
        	                 status = Lib$Free_Lun(con_lun(1))
        	                 status = Lib$Free_Lun(con_lun(2))
	                         another_time = .True.
	                         first_time = .True.
	                         display_menu = .True.
	                         monitor_display = .True.
	                      Else If (ans(1:1) .Eq. 'N' .Or.
	2                              ans(1:1) .Eq. ' ') Then
	                         ans_ok = .True.
	                         another_time = .False.
	                         first_time = .False.
	                      Else
	                         ans_ok = .False.
	                      End If
	                   End Do
	                Else
	                   another_time = .False.
	                End If
	                status = %Loc(FEP_Normal)

	             End If

	         Else

	            If (status .Eq. %Loc(FEP_MonRealTime)) Then
	               status = %Loc(FEP_MonRealTime)
	               Call LIB$Signal(FEP_MonRealTime)
	            Else
	               status = %Loc(FEP_CTOpenHskp)
	               Call LIB$Signal(FEP_CTOpenHskp,%Val(1),
	2                              %Val(ct_status(1)))
	            End If

	         End If

	      End Do

	   Else

	      status = %Loc(FEP_CTInit)
	      Call LIB$Signal(FEP_CTInit,%Val(1),%Val(ct_status(1)))

	   End If

	End If

c
c Terminate ENGPLOTS.
c
	exit_status = FEP_EXIT ( )

	If (status .Eq. %Loc(FEP_Normal) .And.
	2       exit_status .Eq. %Loc(FEP_Normal)) Then
	   Call LIB$Signal(FEP_Normal)
	Else
	   Call LIB$Signal(%Val(SS$_Abort))
	End If


	Stop
	End
