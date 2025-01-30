	Program FIT_Main

C-------------------------------------------------------------------------------
C    PURPOSE: Plot orbit-by-orbit engineering data in user-specified
C             time bins, showing min, max, ave+-sigma, and ave.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Fred Shuman, STX
C            October 23, 1987
C
C    INVOCATION: INPUT PARAMETERS: OUTPUT PARAMETERS: This is a main pgm
C
C    SUBROUTINES CALLED:
C	FIT_Init >- - - - - - - -<FIT_Get_Command
C	FIT_Display_Menu >- - _ _/FIT_Display_Main
C	FIT_Match_Fields	 \FIT_Display_Plot_Menu
C	FIT_Extract_ENG_Stats>- _/FUT_Field_Attributes
C	FIT_Merge_ENG_Stats      \FUT_Combine_HiLo
C	FIT_Plot >- - - - - - - -<PLT
C	FIT_Report
C	FIT_Exit
C	FUT_Timerange
C	FUT_Orbital_Period
C	AUT_DFloat2ADT
C	AUT_ADT2DFloat
C	CT_Init
C	CT_Connect_Read
C	CT_Read_Arcv
C	CT_Close_Arcv
C	CT_Binary_To_GMT
C	LIB$Establish
C	LIB$Get_Lun
C	LIB$Movc3
C	LIB$Signal
C	LIB$Erase_Page
C	STR$UpCase
C
C
C    VARIABLES USED FROM STRUCTURE ORB (declared in FIT_Data):
C	NRECS		I*4		    Number of ENG_STATS records.
C	STATS_DATA	BYTE(512,max,4)    Stats data (512 bytes = 1 orbit)
C	NUMRECS_ORB	R*4(max)\  ^
C	                        |  |
C    INCLUDE FILES:	        |  |
C	FIT_Invoc	        |   \___________________
C	FIT_Data		|			\
C	FUT_Params    (includes fac_engtemp_size=504, fac_max_stats_recs=1000)
C	FUT_Error		+fac_lmac_size=8	(as of 1987 Nov 19)
C	$SSDef
C	CTUser.Inc
C-----------------------------------------------------------------------------
C	Revised:
C	                                        F. Shuman, STX,  1987 Nov 21
C	   Replace Vecplt with PLT for plots.   F. Shuman, STX,  1988 Mar 31
C	   Reset status in loop. Build 4.1 Quickfix
C	                                        Shirley M. Read, 1988 Aug 21
C	   Batch mode failure (SPR 2870)--default Binsize was not being
C	         picked up.                     F. Shuman, STX,  1988 Dec 1
C
C	version 4.1.1 12/01/88, ser 2379, J.T.Bonnell, STX@GSFC
C		This program was modified to refer to the
C		new firas archive logical names in response
C		to SER 2379 (csdr$firas_in, _out, _raw,
C		_ref, _uref, and _cal).
C
C	   SPR 3139.  Reject records with FDQ's "bad record" flag.
C	              This module and FIT_MERGE_ENG_STATS.
C	              F. Shuman, STX,  1989 Feb 23.
C	   SPR 3329.  Get correct Combined Menu field names.  Module FIT_PLOT.
C	              F. Shuman, STX,  1989 Mar 03.
C
C          spr 3690   D. Bouler  May 16, 1989 STX - Infinite loop
C                     If an open failure for the engineering statistics
C                     archive is encountered, program reprompts for
C                     start and stop time if in interactive mode and
C                     times have not been specified; otherwise, exits.
C
C	   SPR 3285.  Get COBE orbital period from FUT_Orbital_Period.
C	              This module and FIT_Get_Command, FIT_Display_Menu,
C	              and FIT_Display_Main were revised.
C	              Fred Shuman, STX  1989 Jul 20.
C
C	version 4.4.1 SER 3306 Q. Chung STX 08/10/89
C		      Provide version number to track software update.
C
C       version 4.4.2 spr 4358 D. Bouler STX Aug 28, 1989
C                              F. Shuman STX Aug 29, 1989
C                     Get start and stop time from command line if /NOINT
C                     qualifier is specified.  Compute binsize after getting
C                     orbital period.  Modified FIT, FIT_Get_Command, FIT.CLD,
C                     FIT_INVOC.TXT.
C
C       version 4.4.2 spr 4356 D. Bouler STX Aug 28, 1989
C                     Add /orbital_pd qualifier to specify orbital period
C                     on command line. Modified fit_invoc.txt to add
C                     flag for qualifier to common block. Modified
C                     fit_get_command so value for fcc_orbit is obtained
C                     from command line if qualifier is specified.
C
C       version 4.5.1 QFix 530, SPRs 4669, 4670.  F. Shuman,  STX,  1989 Oct 2.
C                     Successive cropping of endtime; shuffling of data.
C                     Changes to FIT.for, FIT_Extract_Eng_Stats.for,
C                     FIT_Merge_Eng_Stats.for, and FIT_Data.txt.
C
C       version 4.5.1 QFix xxx, SPR 4762.  F. Shuman,  STX,  1989 Oct 24.
C           FIT was using wrong method to combine hi/lo GRTs.  Eliminate call
C           to FUT_Combine_Temps in favor of a new module FIT_Combine_Temps,
C           written to follow the method used in FUT_Temperature_List.
C           FIT, FIT_Extract_Eng_Stats, FIT_Combine_Temps, and FIT_MSG.
C
C       Spr 6309, version 5.8, Rename the FIT to FIT_MAIN to avoid to
C           use the same name as the plt using.
C           H. Wang, Mar. 8, 1990
C
C       Spr 5655, version 5.8, Add the hard coded menu text file.
C           H. Wang, Mar. 8, 1990
C
C       Spr 6590, version 5.9, correct the menu from 'CHAN_PREAMP_TEMP' to
C              'CHAN_PREAMP_TEMP_A'
C           H. Wang, STX,Mar. 29, 1990
C
C	SPR 6623, vsn 5.10.  When the user specifies bin in orbits, bin size is
C	    computed incorrectly.  Needed to convert orbital period returned
C	    by FUT_Orbital_Period from minutes to seconds.
C	  F. Shuman,  STX,  1990 Apr 03.
C
C       SPR 1860, 7251, 7319.
C              1860: Extend the limit on trendplots from 70 days to
C                    full period of operation.
C              7251, 7319: correct trendplots label
C              H. WANG, STX, Aug. 29, 1990
C
C       Spr 7316, version 6.9, correct the menu from 'CHAN_PREAMP_TEMP_A' to
C              'CHAN_PREAMP_TEMP'
C           H. Wang, STX,Sept. 20, 1990
C
C       Spr 7694,7706, version 7.1.
C            Rename the calibration resistors.
C           H. Wang, STX, Nov. 13, 1990
C
C	SPR 7758, Change the time axis origin from the first trend rec time
C	   to the user-specified start time, so the plots will be easier to
C	   compare.  Change to FIT_Plot.  F. Shuman,  STX,  1990 Nov 21.
C
C       Spr 7692,7720, version 7.3.
C            Correct MTM_POSITION & LMAC Digital Y label.
C           H. Wang, STX, Dec. 5, 1990
C
C       Spr 7892, version 7.3.
C            Correct MTM_CAL_POWER Y label.
C           H. Wang, STX, Dec. 17, 1990
C
C    1991 Sep 25, F. Shuman, STX, SPR 9078.
C       Rename FIT_Combine_Temps routine to FUT_Combine_HiLo and move into FUT
C       to be used by both FIT and FEP.  Change to FIT_Extract_Eng_Stats.
C
C    1991 Nov 26, F. Shuman,  Hughes STX,  SPR 9291.
C       The two reference files GRTwt & GRTtrans that GRTswt split into were
C       erroneously being switched in the routine FIT_Extract_Eng_Stats.
C
C    1991 Dec 3,  F. Shuman,  Hughes STX,  SPR 9327.
C       Make hi/lo-current-combination method in FUT_Combine_HiLo consistent
C       w. that in FUT_Temp_List; remove Ch*1 SIDE from FUT_Combine_HiLo call
C       in FIT_Extract_Eng_Stats.
C
C    1991 Dec 11,  F. Shuman,  Hughes STX,  SPR 9334.
C       Remove obsolete logical name 'CSDR$FIRAS_URef' in FIT_Extract_Eng_Stats
C       --replace with 'CSDR$FIRAS_Ref'.
C
C    1991 Dec 17,  F. Shuman,  Hughes STX,  SPR 9352.
C       Correct some hard-coded menu numbers in FIT_Extract_Eng_Stats which are
C       causing trouble with Bolo-B plots.
C-----------------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Params)'
	Include		'(FUT_Error)'
	Include		'(FIT_Invoc)'
	Include		'(FIT_Data)'
	Include		'($SSDef)'
	Include		'CT$Library:CTUser.Inc'

	Dictionary		'FDQ_ETR'
	Record /FDQ_ETR/	stats_rec
	Record /Orb_Data/	orb

	Integer		*4	eng_lun
	Integer		*4	ios
	Character	*64	infile

	Integer		*4	fr_lun
	Integer		*2	ct_lun
	Integer		*2	ct_status(20)
	Logical		*1	first_time
	Character	*4	ans
	Logical		*1	ans_ok
	Integer		*4	i
	Integer		*4	j
	Integer		*4	status
	Integer		*4	exit_status

	Integer		*4	start_time(2)
	Integer		*4	end_time(2)
	Integer         *4      diff_time(2)
	Character	*14	starting_time
	Character	*14	ending_time
	Character       *14     zero
	Real		*8	R8Start
	Real		*8	R8End
	Real		*8	R8Center
	Integer		*4	ADTCenter(2)
	Character	*14	GMTCenter
	Logical		*1	another_time
	Logical		*1	display_menu

	Integer		*2	menu_page(fac_max_menu_num)
	Integer		*2	menu_element(fac_max_menu_num)
	Integer		*2	nfields
	Integer		*2	length
	Integer         *2      len1
	Integer         *2      len2
	Real		*4	stats_buff(fac_max_stats_recs, 5)
	Character	*14	stats_gmts(fac_max_stats_recs, 3)
	Character	*14	gmtime
	Integer		*4	adtime(2)
	Real		*4	bin_buff(fac_max_stats_recs, num_stats_types)
	Character	*14	bin_gmts(fac_max_stats_recs, 3)
	Real		*4	bin_size
	Integer		*4	nbins
	Integer		*4	plot_group
	Integer		*4	plot_limit

	Character	*14	jstart_default
	Character	*14	jstop_default

	Integer		*4	STR$UpCase
	Integer		*4	FIT_Init
	Integer		*4	FIT_Get_ENG_Fields
	Integer		*4	FIT_Display_Menu
	Integer		*4	FIT_Match_Fields
	Integer		*4	FIT_Extract_ENG_Stats
	Integer		*4	FIT_Merge_ENG_Stats
	Integer		*4	FIT_Plot
	Integer		*4	FIT_Report
	Integer		*4	FIT_Exit
	Integer		*4	FUT_Timerange
	Integer		*4	FUT_Orbital_Period
	Real		*8	AUT_ADT2DFloat
	Integer		*4	LIB$Get_Lun
	Integer		*4	CT_Connect_Read
	Integer		*4	CLI$Get_Value

	Integer		*4	CUT_Register_Version
	Integer		*4	CUT_Display_Banner
	Integer		*4	num_vol/80/
	Integer		*4	lun_out/6/
	Character	*6	version
	Parameter		(version='9.0')

	Data jstart_default/'85001000000000'/
	Data jstop_default/'99365235959990'/
	Data zero         /'00000000000000'/

	External		FIT_Normal
	External		FIT_CTInit
	External		FIT_CTOpenStats
	External		FIT_CTReadStats
	External		FIT_CTCloseStats
	External		FIT_ExitMenu
	External		FIT_QuitMenu
	External		FIT_InvTim
	External		FIT_TruncData
	External		FIT_NoGooDat
	External		FIT_NoData
	External		FIT_InvDataLen
	External		FUT_Error
	External		CT_Connect_Read
	External                CLI$Get_Value

c
c Initialize.
c
	status = CUT_Register_Version(version)
	status = CUT_Display_Banner(lun_out, num_vol,
	1                           'FIRAS Facility  FIT_Instrument_Trends ')
	Write(lun_out,61)
61	Format(//)

	another_time  = .True.
	first_time    = .True.

	nfields = 0
	Do i=1,fac_max_menu_num
	   menu_page(i) = 0
	   menu_element(i) = 0
	End Do

	Call LIB$Establish ( FUT_Error )

	status = FIT_Init ( fr_lun )

	If (status .Eq. %Loc(FIT_Normal)) Then

	   Call CT_Init ( ct_status )

	   If (ct_status(1) .Eq. CTP_Normal) Then

c
c Get the engineering statistics field names to form a menu.
c
C	      status = FIT_Get_ENG_Fields ()

c
c Retrieve the engineering statistics data.
c
	      Do While ( another_time .And. (status .Eq. %Loc(FIT_Normal) .Or.
	2                              status .Eq. %Loc(FIT_CTOpenStats)) )

	         If (first_time) Then

c
c Get the data selection timerange.
c
	            If (fcc_interactive .Eq. fac_present .And.
	1               fcc_jstart_select .Eq. fac_not_present) Then

	               status = FUT_Timerange ( starting_time,
	1                                start_time, ending_time, end_time )
	            Else
	               start_time(1)     = fcc_jstart(1)
	               start_time(2)     = fcc_jstart(2)
	               end_time(1)       = fcc_jstop(1)
	               end_time(2)       = fcc_jstop(2)
	               starting_time     = fcc_jstart_time
	               ending_time       = fcc_jstop_time
	            End If
c
c Get orbital period if it has not been specified on command line.
c FUT_Orbital_Period returns the orb pd in minutes; we deal in sec, so convert:
c
	            If (fcc_orbit_select .Eq. fac_not_present) Then
	               R8Start  = AUT_ADT2DFloat (start_time)
	               R8End    = AUT_ADT2DFloat (end_time)
	               R8Center = (R8Start + R8End)/2.D0
	               Call AUT_DFloat2ADT  (R8Center, ADTCenter)
	               Call CT_Binary_To_GMT(ADTCenter, GMTCenter)
	               status = FUT_Orbital_Period (GMTCenter, fcc_orbit)
	            End If
	            fcc_orbit = 60. * fcc_orbit

	            If (fcc_binorb_select .Eq. fac_present) Then
	               fcc_binsize = fcc_binsize_orbit * fcc_orbit
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
	            Call CT_Binary_To_GMT(end_time,ending_time)
	            fcc_time_range = starting_time // ';' //
	1                             ending_time // ';'
	            first_time = .False.

c
c Open the engineering trends archive.
c
	            infile = 'CSDR$FIRAS_RAW:fdq_etr/' // fcc_time_range
	            status = LIB$Get_Lun(eng_lun)
	            ct_lun = eng_lun

	            Open ( UNIT = eng_lun,
	1                  FILE = infile,
	2                  STATUS = 'old', IOSTAT = ios,
	3                  USEROPEN = CT_Connect_Read )

	         End If   ! (first_time)

	         If ( ios .Eq. 0) Then

c
c Read a record from the archive.
c
	            Type 35
35	            Format (/, ' Fetching engineering trends data.')

c
c Reset the data counter to 0.
c
	            orb.nrecs = 0

	            Do While (ct_status(1) .Ne. CTP_EndOfFile .And.
	1                     orb.nrecs .Lt. fac_max_stats_recs-1)
c
c Get an orbit from the engineering trends record.
c
	               orb.nrecs = orb.nrecs + 1

	               Do i=1,num_stats_types
	                  Call CT_Read_Arcv ( ,
	1                                     ct_lun,
	2                                     stats_rec,
	3                                     ct_status )
	                  Call LIB$Movc3 ( fac_engtemp_size,
	1                                stats_rec.en_analog,
	2                                orb.stats_data(1,orb.nrecs,i) )
	                  Call LIB$Movc3 ( fac_lmaceng_size,
	1                    stats_rec.en_tail.lmac_analog_temp,
	2                    orb.stats_data(fac_engtemp_size+1,orb.nrecs,i) )
	               End Do

	               ADTime(1) = stats_rec.en_head.first_eng_time(1)
	               ADTime(2) = stats_rec.en_head.first_eng_time(2)
	               Call CT_Binary_To_GMT ( ADTime, GMTime )
	               stats_gmts(orb.nrecs, 1) = GMTime

	               stats_gmts(orb.nrecs, 2) = stats_rec.ct_head.gmt

	               ADTime(1) = stats_rec.en_head.last_eng_time(1)
	               ADTime(2) = stats_rec.en_head.last_eng_time(2)
	               Call CT_Binary_To_GMT ( ADTime, GMTime )
	               stats_gmts(orb.nrecs, 3) = GMTime

	            End Do

	            orb.nrecs = orb.nrecs - 1
c
c Display plots for the currently selected timerange.
c
	            If ((ct_status(1) .Eq. CTP_Normal .Or.
	1                ct_status(1) .Eq. CTP_EndOfFile) .And.
	2                orb.nrecs .Ne. 0) Then

	               If (orb.nrecs .Eq. fac_max_stats_recs) Then
	                  Call LIB$Signal ( FIT_TruncData )
	               End If

	               display_menu = .True.

	               Do While (display_menu)

c
c Select fields from a menu or by matching /FIELDS
c with menu fields.
c
	                  If (fcc_interactive .Eq. fac_present) Then
	                     status = FIT_Display_Menu ( menu_page,
	1                                                menu_element,
	2                                                bin_size,
	3                                                nfields )
	                     If (status .Ne. %Loc(FIT_ExitMenu)) Then
	                        display_menu = .False.
	                     End If
	                  Else
	                     status = FIT_Match_Fields ( menu_page,
	1                                                menu_element,
	2                                                nfields )
	                     bin_size = fcc_binsize
	                     display_menu = .False.
	                  End If

c
c Now we will make plots...
c
	                  If (status .Eq. %Loc(FIT_Normal) .Or.
	1                     status .Eq. %Loc(FIT_ExitMenu)) Then

	                     plot_limit = 0

	                     Do While (plot_limit .Lt. nfields)

	                        If (fcc_interactive .Eq. fac_present) Then
	                           plot_group = 1
	                           plot_limit = nfields
	                        Else
	                           plot_group = plot_limit + 1
	                           plot_limit = plot_group
	                        End If

c
c ...First, collect the raw orbit data;
c
	                        Do i=plot_group,plot_limit

	                           status = FIT_Extract_Eng_Stats (
	1                                       fr_lun,
	2                                       menu_page(i),
	3                                       menu_element(i),
	4                                       stats_gmts,
	5                                       orb,
	6                                       stats_buff )
c
c    then, 'compute' them into the form to be sent to PLT for plotting;
c
	                           If (status .Eq. %Loc(FIT_Normal)) Then
	                              status = FIT_Merge_Eng_Stats (
	1                                       stats_buff, stats_gmts,
	2                                       bin_size,
	3                                       starting_time, ending_time,
	4                                orb.nrecs, bin_buff, bin_gmts, nbins)
c
c    finally, plot 'em!
c
	                              If (status .Eq. %Loc(FIT_Normal)) Then
	                                 status = FIT_Plot (
	1                                          bin_buff,
	2                                          bin_gmts,
	3                                          bin_size,
	4                                          nbins,starting_time,
	5                                          menu_page(i),
	6                                          menu_element(i) )
	                              Else If (status .Eq. %Loc(FIT_NoGooDat))
	1                                     Then
	                                 Call LIB$Signal (FIT_NoGooDat)
	                                 If(fcc_interactive.Eq.fac_present)Then
	                                    Type *,'Press any key to continue'
	                                    Accept 45,ans
45	                                    Format (a4)
	                                 End If
	                              End If
	                           End If

	                        End Do

	                     End Do

	                     If (nfields .Ne. 0 .And.
	1                        status .Eq. %Loc(FIT_Normal)) Then
	                        status = FIT_Report ( menu_page,
	1                                             menu_element,
	2                                             nfields )
	                     End If

	                  End If

	               End Do

	            Else

	               If (orb.nrecs .Eq. 0) Then
	                  Call LIB$Signal ( FIT_NoData )
	               Else
	                  status = %Loc(FIT_CTReadStats)
	                  Call LIB$Signal(FIT_CTReadStats,%Val(1),
	1                                 %Val(ct_status(1)))
	               End If

	            End If

c
c Prepare for the next set of plots.
c
	            Call CT_Close_Arcv ( , ct_lun, ct_status )

	            If (ct_status(1) .Eq. CTP_Normal) Then

	               If (fcc_interactive .Eq. fac_present .And.
	1                    (fcc_jstart_select .Eq. fac_not_present .And.
	2                     fcc_jstop_select .Eq. fac_not_present)) Then
	                  ans_ok = .False.
	                  Do While (.Not. ans_ok)
	                     Call LIB$Erase_Page(1,1)
	                     Type 30
30	                     Format (' Do you want to select '
	1                            'another timerange? (Y/[N]): ', $)
	                     Accept 25, ans
25	                     Format (a)
	                     status = STR$UpCase(ans,ans)
	                     If (ans(1:1) .Eq. 'Y') Then
	                        ans_ok = .True.
	                        another_time = .True.
	                        first_time = .True.
	                     Else If (ans(1:1) .Eq. 'N' .Or.
	1                             ans(1:1) .Eq. ' ') Then
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

	               status = %Loc(FIT_Normal)

	            Else

	               status = %Loc(FIT_CTCloseStats)
	               Call LIB$Signal(FIT_CTCloseStats,%Val(1),
	1                              %Val(ct_status(1)))

	            End If

	         Else

		    status = %Loc(FIT_CTOpenStats)
	            Call LIB$Signal(FIT_CTOpenStats,%Val(1),
	1                           %Val(ios))

c
c If interactive, and start/stop not present,
c reset flag for first time
c so that start and stop times
c will again be prompted for.
c Otherwise, set another_time flag to exit loop.
c

	           If (fcc_interactive   .Eq. fac_present     .And.
	1              fcc_jstart_select .Eq. fac_not_present .And.
	2              fcc_jstop_select  .Eq. fac_not_present) Then

	                Write(6,*) 'No data found. Enter valid time range.'
	                first_time = .True.
	           Else
	                another_time = .False.
	           End If

	         End If  ! (ios .Eq. 0)

	      End Do     ! (another_time ...)

	   Else

	      status = %Loc(FIT_CTInit)
	      Call LIB$Signal(FIT_CTInit,%Val(1),%Val(ct_status(1)))

	   End If   ! (ct_status .Eq. ctp_normal)

	End If      ! (status .Eq. fit_normal)

c
c Terminate TRENDPLOTS.
c
	exit_status = FIT_EXIT ( )

	If ( status .Eq. %Loc(FIT_Normal)       .And.
	2    exit_status .Eq. %Loc(FIT_Normal) ) Then
	   Call LIB$Signal(FIT_Normal)
	Else
	   Call LIB$Signal(%Val(SS$_Abort))
	End If

	Stop
	End
