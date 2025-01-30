	Integer*4 Function FIT_Plot ( plot_buff, bin_gmts, bin_size, nrecs,
	2                             starting_time, menu_page, menu_element )

C------------------------------------------------------------------------------
C    PURPOSE: Generates the engineering data plots for the TRENDPLOTS facility.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Fred Shuman	(adapted from FEP_Plot
C            STX		 by Rob Kummerer, STX, April 13, 1987)
C            October 5, 1987
C
C    INVOCATION: STATUS = FIT_PLOT ( PLOT_BUFF, BIN_GMTS, BIN_SIZE, NRECS,
C				     STARTING_TIME, MENU_PAGE, MENU_ELEMENT )
C
C    INPUT PARAMETERS:
C	PLOT_BUFF(max,5)	 R*4		Plot data buffer (Y) for the 4
C						   stats types--min,mean,s,max.
C	BIN_GMTS(max, 3)	Ch*14		Times(X) corres to each data pt
C						   3 times/pt: min, mean, max
C	BIN_SIZE		 R*4		User-chosen time bin (sec).
C	NRECS			 I*4		Number of data points to plot.
C	STARTING_TIME		Ch*14		User-supplied start time.
C	MENU_PAGE		 I*2		Menu page selected.
C	MENU_ELEMENT		 I*2		Menu field selected.
C		Note: "max" is fac_max_stats_recs, declared in FUT_PARAMS.TXT
C				= 1000 as of 1988 Jun 6
C    OUTPUT PARAMETERS: None
C
C    SUBROUTINES CALLED:
C	STR$Trim
C	PLT
C	LIB$Erase_Page
C
C    COMMON VARIABLES USED:
C	FCC_INTERACTIVE		 I*4		Interactive/Batch flag.
C	FCC_PLOT_DEVICE		Ch*32		Initial PLT device type.
C	FCC_PLT_COM		 I*4		Select flag.
C	FCC_PLT_COM_FILE	Ch*64		PLT command file name.
C	GRTA_MENU(32)		Ch*32		GRT side A menu.
C	GRTB_MENU(32)		Ch*32		GRT side B menu.
C	COMB_MENU(32)		Ch*32		Combined GRT menu.
C	IPDU_MENU(14)		Ch*32		IPDU temps menu.
C	VOLT_MENU(20)		Ch*32		Voltages menu.
C	CURR_MENU(28)		Ch*32		Currents menu.
C
C    INCLUDE FILES:
C	FUT_Error
C	FUT_Params
C	FIT_Invoc
C	FIT_Menu
C
C-------------------------------------------------------------------------------
C	Revisions:                           F. Shuman, 1987 Nov 21
C	   Replace Vecplt with PLT.  F. Shuman, 1988 Mar 28
C	   SPR 3329.  Get correct Combined Menu field names.
C	              F. Shuman,  STX,  1989 Mar 03.
C
C	   SER 5727.  Add capability to pass in a command file to PLT.
C		      Steven Alexander, STX,  1990 Feb 28.
C
C	   SPR 7251,7319.  correct the label of FIT plots
C                      H. WANG, STX, Aug. 29, 1990
C
C	   SPR 7758, Change the time axis origin from the first trend rec time
C		to the user-specified start time, so the plots will be easier to
C		compare.  F. Shuman,  STX,  1990 Nov 21.
C
C	   SPR 7892.  correct the label of FIT plots
C                      H. WANG, STX, Dec. 17, 1990
C
C-------------------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FUT_Params)'
	Include		'(FIT_Invoc)'
	Include		'(FIT_Menu)'

	Integer		*4	i
	Integer		*4	j
	Integer		*4	k
	Integer		*4	status

	Real		*4	plot_buff(fac_max_stats_recs, 5)
	Character	*14	bin_gmts(fac_max_stats_recs, 3)
	Real		*4	bin_size
	Integer		*4	nrecs
	Character	*14	starting_time
	Integer		*2	menu_page
	Integer		*2	menu_element

	Real		*4	y(fac_max_stats_recs, 6)
	Integer		*4	iery(5)
	Integer		*4	nveci
	Character	*80	cmd(20)
	Integer		*4	ncmd
	Integer		*4	ier

	Logical		*1	moreplot
	Real		*4	maxtime
	Real		*4	time_unit
	Real		*4	timeaxis(fac_max_stats_recs)
	Real		*4	timeplotmin
	Real		*4	timeplotmax
	Real		*4	ymin
	Real		*4	yplotmin
	Real		*4	ymax
	Real		*4	yplotmax

	Real		*4	timeorigin
	Integer		*4	iy
	Integer		*4	refyr
	Integer		*4	res
	Integer		*4	iyr
	Integer		*4	iday
	Integer		*4	ihr
	Integer		*4	imin
	Integer		*4	msec

	Character	*100	title1
	Character	*100	title2
c	Character	*70	legend
	Integer		*4	len
	Character	*80	xlabel
	Character	*30	ylabel
	Character	*30	sylabel

	Integer		*4	STR$Trim

	External	FIT_Normal
	External	FIT_MenuElErr

c
c Set the plot title and axis labels.
c
	If (menu_page .Eq. GRTA_page) Then

	   title1 = GRTA_menu(menu_element)
	   status = STR$Trim ( title1, title1, len )
	   sylabel = 'Kelvin'
	   If (menu_element .ge.5 .and. menu_element .le. 8) Then
	      sylabel = 'Counts'
           Endif
	   If (menu_element .ge.21 .and. menu_element .le. 24) Then
	      sylabel = 'Counts'
           Endif
	Else If (menu_page .Eq. GRTB_page) Then

	   title1 = GRTB_menu(menu_element)
	   status = STR$Trim ( title1, title1, len )
	   sylabel = 'Kelvin'
	   If (menu_element .ge.5 .and. menu_element .le. 8) Then
	      sylabel = 'Counts'
           Endif
	   If (menu_element .ge.21 .and. menu_element .le. 24) Then
	      sylabel = 'Counts'
           Endif

	Else If (menu_page .Eq. COMB_page) Then

	   title1 = COMB_menu(menu_element)
	   status = STR$Trim ( title1, title1, len )
	   sylabel = 'Kelvin'
	   If (menu_element .ge.5 .and. menu_element .le. 8) Then
	      sylabel = 'Counts'
           Endif
	   If (menu_element .ge.21 .and. menu_element .le. 24) Then
	      sylabel = 'Counts'
           Endif

	Else If (menu_page .Eq. IPDU_page) Then

	   title1 = IPDU_menu(menu_element)
	   status = STR$Trim ( title1, title1, len )
	   sylabel = 'Degrees C'
	   If (menu_element .ge.1 .and. menu_element .le. 4) Then
	      sylabel = 'VOLTS'
	   Else If (menu_element .Eq.12 .Or. menu_element .Eq. 13) Then
	      sylabel = 'Micro-Amps\'
	   Else If (menu_element .Eq.24 .Or. menu_element .Eq. 25) Then
	      sylabel = 'Micro-Amps'
	   Else If (menu_element .Ge. 14 .And. menu_element .Le. 21) Then
	      sylabel = 'STATUS\'
	   Else If (menu_element .Eq. 26 .or. menu_element .Eq. 27) Then
	      sylabel = 'Eng Units'
	   End If

	Else If (menu_page .Eq. VOLT_page) Then

	   title1 = VOLT_menu(menu_element)
	   status = STR$Trim ( title1, title1, len )
	   sylabel = 'Volts'

	Else If (menu_page .Eq. CURR_page) Then

	   title1 = CURR_menu(menu_element)
	   status = STR$Trim ( title1, title1, len )
	   sylabel = 'Amps'
	   If (menu_element .ge. 13 .and. menu_element .le. 14) Then
	      sylabel = 'Degrees C'
	   End If

	End If

	title2 = 'Minimum, Mean - s, Mean, Mean + s, Maximum'

c
c Convert the timetags to real numbers in seconds.  Note that the last 5 digits
c  can be taken together as milliseconds.  "refyr" is the yr field of the first
c  GMT, used as a reference year; "res" is the residue mod 4 of (refyr-1), so
c  chosen that, when added to the expression inside the INT, it will trigger a
c  leap day at the correct year-end, and 89001 minus 88366 will be 1 day.
c
	Do i=1,nrecs
	   Read (bin_gmts(i,2),35) iyr, iday, ihr, imin, msec
35	   Format (i2,i3,i2,i2,i5)
	   If (i .Eq. 1) Then
	      refyr = iyr
	      res = refyr-1 - 4*INT( (refyr-1)/4 )
	   End If
	   iy = iyr - refyr
	   timeaxis(i) = (INT((iy*1461+res)/4) + iday)*fac_day
	2                + ihr*fac_hour + imin*fac_minute + msec*1.e-3
	End Do
c
c Reset the timeaxis reference (origin) to the user-supplied start time.
c
	Read (starting_time,35) iyr, iday, ihr, imin, msec
	iy = iyr - refyr
	timeorigin = (INT((iy*1461+res)/4) + iday)*fac_day
	2            + ihr*fac_hour + imin*fac_minute + msec*1.e-3

	Do i=1,nrecs
	   timeaxis(i) = timeaxis(i) - timeorigin
	End Do

	If (nrecs .Eq. 1) Then
	   maxtime = bin_size
	Else
	   maxtime = timeaxis(nrecs)
	End If

	If (maxtime .Ge. 3*fac_day) Then
	   time_unit = fac_day
	   xlabel = 'Time (days)'
	Else If (maxtime .Ge. 3*fac_hour) Then
	   time_unit = fac_hour
	   xlabel = 'Time (hours)'
	Else If (maxtime .Ge. 3*fac_minute) Then
	   time_unit = fac_minute
	   xlabel = 'Time (minutes)'
	Else
	   time_unit = 1.
	   xlabel = 'Time (seconds)'
	End If

	status = STR$Trim ( xlabel, xlabel, len )
	xlabel = xlabel(1:len) // ' from ' // starting_time
	2                      // ' to ' // bin_gmts(nrecs,2)

	Do i=1,nrecs
	   timeaxis(i) = timeaxis(i) / time_unit
	End Do
c
c Make the plot...
c
c     Find the left and right plot boundaries.
c
	timeplotmin = - (fac_margin * maxtime) / time_unit
	timeplotmax = (1. + fac_margin)* maxtime / time_unit

	moreplot = .True.

c
c Set up the PLT calling arguments, y(), iery(), nveci, cmd(), and ncmd
c     and find the top and bottom plot boundaries.
c
	ymin = plot_buff(1,1)
	ymax = plot_buff(1,1)
	Do j=1,nrecs
	   y(j,1) = timeaxis(j)
	   Do k=2,6
	      y(j,k) = plot_buff(j,k-1)
	   End Do
	   ymin = amin1(ymin, y(j,2))
	   ymin = amin1(ymin, y(j,3)-y(j,4))
	   ymax = amax1(ymax, y(j,5))
	   ymax = amax1(ymax, y(j,3)+y(j,4))
	End Do

	If (ymin .Eq. ymax) Then
	   If (ymin .Gt. 0.) Then
	      yplotmin = (1. - fac_margin)*ymin
	      yplotmax = (1. + fac_margin)*ymax
	   Else If (ymin .Eq. 0.) Then
	      yplotmin = -1.
	      yplotmax =  1.
	   Else
	      yplotmin = (1. + fac_margin)*ymin
	      yplotmax = (1. - fac_margin)*ymax
	   End If
	Else
	   yplotmin = ymin - fac_margin*(ymax - ymin)
	   yplotmax = ymax + fac_margin*(ymax - ymin)
	End If

	ylabel = sylabel

	if (fcc_interactive .eq. fac_present) Call LIB$Erase_Page ( 1, 1 )
	If (nrecs .Eq. 1) Then
	   nrecs = 2
	   y(2,1) = timeaxis(1) + bin_size/time_unit
	   y(2,2) = -1.
	   y(2,3) = -1.
	   y(2,4) =  0.
	   y(2,5) = -1.
	   y(2,6) = -1.
	End If

	iery(1) = 0
	iery(2) = 0
	iery(3) = 1
	iery(4) = 0
	iery(5) = 0

	nveci = 5

	cmd(1)  = 'D ' // fcc_plot_device
	cmd(2)  = 'LA OT ' // title1
	cmd(3)  = 'LA  T ' // title2
	cmd(4)  = 'LA  X ' // xlabel
	cmd(5)  = 'LA  Y ' // ylabel
	cmd(6)  = 'LA  F'

	If (fcc_interactive .Eq. fac_present) Then
	   cmd(7) = 'CS 1.5'
	Else
	   cmd(7) = 'CS 1'
	End If

	cmd(8)  = 'V .2 .1 .9 .8'
	Write (cmd(9), '(a, 4(1x, e14.4))')
	2      'R', timeplotmin, timeplotmax, yplotmin, yplotmax
	cmd(10) = 'LI ON'
	cmd(11) = 'CO 0 ON 5'
	cmd(12) = 'PLOT OVERLAY'

        ncmd = 12

        If (fcc_plt_com .Eq. fac_present) Then
            ncmd = ncmd + 1
            cmd(ncmd) = '@' // fcc_plt_com_file
        End If

        If (fcc_interactive .Ne. fac_present) Then
	    ncmd = ncmd + 1
            cmd(ncmd) = 'p'
            ncmd = ncmd + 1
            cmd(ncmd) = 'q'
        End If
c
c Display the plot until the user finishes looking at it.
c
	Call PLT (y, iery, fac_max_stats_recs, nrecs, nveci, cmd, ncmd, ier)
        If (fcc_interactive .Eq. fac_present) Call LIB$Erase_Page ( 1, 1 )

	FIT_Plot = %Loc(FIT_Normal)

	Return
	End
