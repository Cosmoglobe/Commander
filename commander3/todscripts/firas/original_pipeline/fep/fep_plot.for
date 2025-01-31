	Integer*4 Function FEP_Plot ( plot_buff, telm_mode,dwell_mode, timeaxis,
	2                             tag_first, tag_last, nrecs, nfields,
	3                             menu_page, menu_element, convert_flag )

C-------------------------------------------------------------------------------
C    PURPOSE: Generates the engineering data plots.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX
C            April 13, 1987
C
C    INVOCATION: STATUS = FEP_PLOT ( PLOT_BUFF, Telm_mode,DWELL_MODE, TIMEAXIS,
C	                             TAG_FIRST, TAG_LAST, NRECS, NFIELDS,
C	                             MENU_PAGE, MENU_ELEMENT, CONVERT_FLAG )
C
C    INPUT PARAMETERS:
C	PLOT_BUFF(max, 9)	 R*4		Plot data buffer (Y).
C	DWELL_MODE(max, 3)	Byte		Dwell flag for sel'd grp of flds
C				       (*,1) --> dwell on some chosen A field
C				       (*,2) --> dwell on some chosen B field
C				       (*,3) --> dwell on some chosen field
C	TIMEAXIS()		 R*4		Time(X) in sec since 86001.
C	TAG_FIRST		Ch*14		Start time in "Julian" format.
C	TAG_LAST		Ch*14		Stop time in "Julian" format.
C	NRECS			 I*4		Number of data points.
C	NFIELDS			 I*4		Number of fields to plot.
C	MENU_PAGE(*)		 I*2		Menu page(s) selected.
C	MENU_ELEMENT(*)		 I*2		Menu fields selected.
C	CONVERT_FLAG(*)		 L*1		[No]convert interactive flags.
C
C    OUTPUT PARAMETERS: None
C
C    SUBROUTINES CALLED:
C	STR$Trim
C	PLT
C	FEP_Convolve
C	LIB$Erase_Page
C
C    COMMON VARIABLES USED:
C	FCC_INTERACTIVE		 I*4		Interactive/Batch flag.
C	FCC_CONVERT		 I*4		Convert to EU flag.
C	FCC_SMOOTH		 I*4		Smoothing selection flag.
C	FCC_CROSSGAPS		 I*4		Gap-crossing selection flag.
C	FCC_MONITOR		 I*4		Monitor field flag.
C	FCC_PLOT_DEVICE		Ch*32		Initial PLT device type.
C	FCC_PLOT_COM_FILE	Ch*64	        PLT command file name
C	FCC_PLOT_COM		I*4		Select flag
C	HSKP_MENU(28)		Ch*28		Housekeeping header menu.
C	GRTA_MENU(32)		Ch*32		GRT side A menu.
C	GRTB_MENU(32)		Ch*32		GRT side B menu.
C	COMB_MENU(32)		Ch*32		Combined GRT menu.
C	IPDU_MENU(14)		Ch*32		IPDU temps menu.
C	VOLT_MENU(20)		Ch*32		Voltages menu.
C	CURR_MENU(18)		Ch*32		Currents menu.
C
C    INCLUDE FILES:
C	FUT_Error
C	FUT_Params
C	FEP_Invoc
C	FEP_Menu
C
C    Software Maintainers, CAUTION:  Hardcoded numbers used in tests of
C    menu_element().  Look for the comment block:
C	c
C	c  BEWARE--HARDCODED MENU ELEMENT NUMBERS!  REVISE WHEN MENU CHANGES!
C	c
C
C-------------------------------------------------------------------------------
C
C Changes:
C
C	Wait 32 seconds (one major frame) for next set of data to arrive.
C	November 4, 1987, R. Kummerer.
C
C	Display dwell mode regimes.  November 5, 1987, R. Kummerer.
C
C	Replace Vecplt with PLT for doing plots.  1988 April 8, F. Shuman.
C
C	Added the capability of plotting the fields in the overlay PLT mode
C	with up to 5 different line styles available for the plot package
C	for the interactive runs as well as batch jobs. The current limit
C	of fields selected by the query mode is 4. For batch jobs with more
C	than 5 fields in overlay mode, the line style sequence is repeated.
C	A command was added to label the line sequence key below the X axis.
C	In order to display this line, the character size was changed to 1.0
C	for both batch or interactive modes. This change also ensures that
C	the field names in the title will not be truncated as they were with
C	the previous version. Eventually the number of fields allowed on
C	overlay plots may have to be limited to 4 or 5 because of the limit
C	of 5 on line styles and the length of the field names. The new PLT
C	may allow more flexibility.  September 21, 1988, Shirley M. Read, STX.
C
C	Label each curve just inside the left and right plot boundaries with
C	its numerical order in the field list.  (SPR 2580)  F. Shuman, STX,
C	1988 Nov 1.
C
C	Reverse the order of arguments in LineStyle command sent to PLT, due to
C	syntax change in newly delivered PLT.  Old:  LS <vector#> <LStyle#>;
C	new:  LS <LStyle#> ON <vector#>.  Module FEP_PLOT.  (SPR 2650)
C	F. Shuman, STX, 1988 Nov 1.
C
C	Add second index to dwell_mode(*,*) going 1 to 3 indicating separately
C	the presence of dwell data on A (*,1) side, B (*,2) side, and A or B
C	side (*,3).  Using this information, the presence of GRTA dwell will not
C	blank out GRTB or non-GRT fields.  (SPR 2651)  F. Shuman, STX,
C	1988 Nov 1.
C
C	Increase the index limits in declarations of plot_buff() and cmd() to
C	allow the full 9 fields to be plotted.  (SPR 2673)  F. Shuman, STX,
C	1988 Nov 1.
C
C	Correct y-scale and "dwellval" determinations.  (SPRs 3101, 3102)
C	F. Shuman, STX, 1989 Jan 6.
C
C       SPR 4852, correct to plot the good data and don't plot the bad and
C       dwell data. Qfix 617, Harte Wang, STX, 1989 Nov 3.
C
C       QFIX 636 / SPR 4958.  Engplots issues YMIN = YMAX error from PLT.
C       Rescale values sent to PLT were increased in precision from E14.4 to
C       E16.6.  Module FEP_Plot.  Fred Shuman, STX, 1989 Nov 6.
C
C       SPR 5068, avoid plotting the data points when telm. qual. is bad
C       Qfix 766, Harte Wang, STX, 1989 Dec 2.
C
C       SPR 5122, bad label units, The plots of the LMAC temperatures plotted
C       with Engplots incorrectly label them as in Kelvin(rather than correct 
C       degrees C). 
C       Harte Wang, STX, 1989 Dec 8.
C
C       SPR 5266, FEP incorrect label for Y-axis
C       Harte Wang, STX 1989 Dec. 20 
C
C       SER 5725, Add the capability to access PLT command file from CLD.
C       Harte Wang, STX 1990 Mar. 2 
C
C       SPR 6884, Correct the hot spot engplots label on the "y" axis
C       Harte Wang, STX 1990 June. 14 
C
C       SPR 7369, Incorrect labelling of the CAL resistors
C       Harte Wang, STX 1990 Aug. 30 
C
C	SPR 2690, Change the time axis origin from the first HKP rec time to
C	the user-specified start time, so the plots will be easier to compare.
C	F. Shuman,  STX,  1990 Nov 12.
C-------------------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FUT_Params)'
	Include		'(FEP_Invoc)'
	Include		'(FEP_Menu)'

	Integer		*4	i
	Integer		*4	j
	Integer		*4	k
	Integer		*4	m
	Integer		*4	status

	Real		*4	plot_buff(fac_max_hskp_recs, 9)
	Byte			dwell_mode(fac_max_hskp_recs, 3)
	Byte			telm_mode(fac_max_hskp_recs)
	Real		*4	timeaxis(fac_max_hskp_recs)
	Character	*14	tag_first
	Character	*14	tag_last
	Integer		*4	nrecs
	Integer		*4	nfields
	Integer		*2	menu_page(*)
	Integer		*2	menu_element(*)
	Logical		*1	convert_flag(*)

	Character	*300	title
	Character	*100	title1
	Character	*100	title2
	Character	*100	title3
	Character	*100	Legend
	Integer		*4	len
	Integer		*4	num1
	Integer		*4	num2
	Character	*80	xlabel
	Character	*30	ylabel(9)
	Character	*30	sylabel(9)

	Byte			dwellbuff(fac_max_hskp_recs)
	Integer		*4	nplotpts
	Integer		*4	ismooth
	Integer		*4	iwindow
	Integer		*4	select_window
	Logical		*1	moreplot
	Character	*4	smooth_label
	Logical		*1	crossgaps
	Character	*2	ans
	Real		*4	time_unit
	Logical		*1	same_units

	Real		*4	maxtime
	Real		*4	timeplotmin
	Real		*4	timeplotmax
	Real		*4	xleftlabel
	Real		*4	xrightlabel
	Real		*4	ymin
	Real		*4	ymax
	Real		*4	yplotmin
	Real		*4	yplotmax
	Real		*4	ymn(9)
	Real		*4	ymx(9)
	Real		*4	yplotmn(9)
	Real		*4	yplotmx(9)
	Real		*4	y(fac_max_hskp_recs, 10)
	Integer		*4	iery(10)
	Integer		*4	nveci
	Character	*150	cmd(55)
	Integer		*4	ncmd
	Integer		*4	ier
	Integer		*4      ls
	Logical		*1	dwell_found
	Logical         *1      found
	Logical		*1	dwell_present
	Logical		*1	telm_present
	Real		*4      dwellval
	Real		*4      Telmval
        
	Integer		*4	STR$Trim
	Integer         *4      kk
	Integer         *4      kkl
	Integer         *4      kkr
	External	FEP_Normal
	External	FEP_MenuElErr

c
c Set the plot title and axis labels.  "m" is keeping track of the position of
c                                       the "\" marking the end of 'title'.
	title  = '\'
	m = 1
	num1 = 0
	num2 = 0

	Do k = 1, nfields
	   sylabel(k) = 'Counts\'
c
c  BEWARE--HARDCODED MENU ELEMENT NUMBERS!  REVISE WHEN MENU CHANGES!
c
	   If (menu_page(k) .Eq. HSKP_page) Then

	      sylabel(k) = 'Status\'
	      if(menu_element(k) .Eq. 15 .Or. menu_element(k) .Eq. 16)then
	        sylabel(k) = 'Eng Units\'
	      endif
	      title(m:) = HSKP_menu(menu_element(k))
	      status = STR$Trim ( title, title, len )
	      title(len+1:) = ', \'
	      If (((fcc_interactive .Eq. fac_present .And.
	2           .not. convert_flag(k)) .Or.
	3          (fcc_interactive .Eq. fac_not_present .And.
	4            fcc_convert .Eq. fac_present)) .And.
	5         (menu_element(k) .Eq. 15 .Or. menu_element(k) .Eq. 16)) Then
	         sylabel(k) = 'Counts\'
	      End If

	   Else If (menu_page(k) .Eq. GRTA_page) Then

	      title(m:) = GRTA_menu(menu_element(k))
	      status = STR$Trim ( title, title, len )
	      title(len+1:) = ', \'
	      If ((fcc_interactive .Eq. fac_present .And.
	2          convert_flag(k)) .Or.
	3         (fcc_interactive .Eq. fac_not_present .And.
	4          fcc_convert .Eq. fac_present)) Then
	         sylabel(k) = 'Kelvin\'
                if (menu_element(k) .ge. 5 .and. menu_element(k) .le. 8)
	1	   sylabel(k) = 'Counts\'
                if (menu_element(k) .ge. 21 .and. menu_element(k) .le. 24)
	1	   sylabel(k) = 'Counts\'
	      End If

	   Else If (menu_page(k) .Eq. GRTB_page) Then

	      title(m:) = GRTB_menu(menu_element(k))
	      status = STR$Trim ( title, title, len )
	      title(len+1:) = ', \'
	      If ((fcc_interactive .Eq. fac_present .And.
	2          convert_flag(k)) .Or.
	3         (fcc_interactive .Eq. fac_not_present .And.
	4          fcc_convert .Eq. fac_present)) Then
	         sylabel(k) = 'Kelvin\'
                if (menu_element(k) .ge. 5 .and. menu_element(k) .le. 8)
	1	   sylabel(k) = 'Counts\'
                if (menu_element(k) .ge. 21 .and. menu_element(k) .le. 24)
	1	   sylabel(k) = 'Counts\'
	      End If

	   Else If (menu_page(k) .Eq. COMB_page) Then

	      title(m:) = COMB_menu(menu_element(k))
	      status = STR$Trim ( title, title, len )
	      title(len+1:) = ', \'
	      If ((fcc_interactive .Eq. fac_present .And.
	2          convert_flag(k)) .Or.
	3         (fcc_interactive .Eq. fac_not_present .And.
	4          fcc_convert .Eq. fac_present)) Then
	         sylabel(k) = 'Kelvin\'
                if (menu_element(k) .ge. 5 .and. menu_element(k) .le. 8)
	1	   sylabel(k) = 'Counts\'
                if (menu_element(k) .ge. 21 .and. menu_element(k) .le. 24)
	1	   sylabel(k) = 'Counts\'
	      End If

	   Else If (menu_page(k) .Eq. IPDU_page) Then
c
c  BEWARE--HARDCODED MENU ELEMENT NUMBERS!  REVISE WHEN MENU CHANGES!
c
	      title(m:) = IPDU_menu(menu_element(k))
	      status = STR$Trim ( title, title, len )
	      title(len+1:) = ', \'
	      If ((fcc_interactive .Eq. fac_present .And.
	2          convert_flag(k)) .Or.
	3         (fcc_interactive .Eq. fac_not_present .And.
	4          fcc_convert .Eq. fac_present)) Then
	              sylabel(k) = 'Degrees C\'
	              if(menu_element(k) .Eq. 8 .Or. menu_element(k) .Eq. 9)then
	                sylabel(k) = 'Milli-Amps\'
	              endif
	      End If

	   Else If (menu_page(k) .Eq. VOLT_page) Then

	      title(m:) = VOLT_menu(menu_element(k))
	      status = STR$Trim ( title, title, len )
	      title(len+1:) = ', \'
	      If ((fcc_interactive .Eq. fac_present .And.
	2          convert_flag(k)) .Or.
	3         (fcc_interactive .Eq. fac_not_present .And.
	4          fcc_convert .Eq. fac_present)) Then
	         sylabel(k) = 'Volts\'
	      End If

	   Else If (menu_page(k) .Eq. CURR_page) Then
c
c  BEWARE--HARDCODED MENU ELEMENT NUMBERS!  REVISE WHEN MENU CHANGES!
c
	      if(menu_element(k) .ge. 21 .and. menu_element(k) .le. 32)then
	        sylabel(k) = 'Status\'
	      endif
	      title(m:) = CURR_menu(menu_element(k))
	      status = STR$Trim ( title, title, len )
	      title(len+1:) = ', \'
	      If ((fcc_interactive .Eq. fac_present .And.
	2          convert_flag(k)) .Or.
	3         (fcc_interactive .Eq. fac_not_present .And.
	4          fcc_convert .Eq. fac_present)) Then
	             sylabel(k) = 'Amps\'
	             if(menu_element(k) .ge. 5 .and. menu_element(k) .le. 8)then
	               sylabel(k) = 'Volts\'
	             elseif(menu_element(k) .ge. 21 .and. menu_element(k) .le. 32)then
	               sylabel(k) = 'Status\'
	             elseif(menu_element(k) .eq. 33)then
	               sylabel(k) = 'Quality\'
	             elseif(menu_element(k) .ge. 17 .and. menu_element(k) .le. 18)then
	               sylabel(k) = 'Degrees C\'
	             endif
	      End If

	   End If

	   m = Index(title, '\')
c
c Find the character positions, num1 and num2, in the string "title", of the end
c   of the last field names that fit into the 1st two 66-chr lines (68 includes
c   the " \")
c
	   If (m .Le. 68) Then
	      num1 = m - 2
	      num2 = m - 2
	   Else If (m-num1-1 .Le. 68) Then
	      num2 = m - 2
	   End If

	End Do  ! k = 1, nfields
c
c Trim the final ", \" off the title string
c
	len = m - 3
	title(len+1:) = ' '

	If (nrecs .Eq. 1) Then
	   maxtime = fac_minute
	Else
	   maxtime = timeaxis(nrecs)
	End If

	If (timeaxis(nrecs) .Le. 3.*fac_hour) Then
	   time_unit = fac_minute
	   xlabel = 'Time (minutes)  from ' // tag_first // ' to ' // tag_last
	Else
	   time_unit = fac_hour
	   xlabel = 'Time (hours)  from ' // tag_first // ' to ' // tag_last
	End If
c
c Load time into y(*,1)
c
	Do i=1,nrecs
	   If (timeaxis(i) .Eq. fac_no_data) Then
	      y(i,1) = fac_no_data
	   Else
	      y(i,1) = timeaxis(i) / time_unit
	   End If
	End Do

	nplotpts = nrecs
	moreplot = .True.
	ismooth = 0
	iwindow = 1
	select_window = 1
	crossgaps = .True.
c
c Force smoothing if selected.
c
	If (fcc_smooth .Eq. fac_present) Then
	   ismooth = 1
	   select_window = 0
	   If (fcc_crossgaps .Eq. fac_not_present) Then
	      crossgaps = .False.
	   End If
	End If
c
c Display the plot until the user finishes looking at it.
c
	Do While (moreplot)

	   same_units = .True.
	   dwell_present = .False.
	   Telm_present = .False.
	   dwell_found = .False.
          
	   Do k=1,nfields
               
	      Do j=1,nplotpts
                 
c
c Load eng. data into y(*,*); dwell flag into dwellbuff(*).
c
	         y(j,k+1) = plot_buff(j,k)
	         If (plot_buff(j,k) .eq. -99999.0) y(j,k+1) = fac_no_data
	         dwellbuff(j) = 0
 	         If (dwell_mode(j,3) .Eq. 1) Then
	            dwell_present = .True.
	            y(j, nfields+2) = 1.
	            If (.Not. dwell_found) Then
	               dwell_found = .True.
	            End If
c
c Mark dwell mode data.
c
c  BEWARE--HARDCODED MENU ELEMENT NUMBERS!  REVISE WHEN MENU CHANGES!
c
	            If ( dwell_mode(j,1) .Eq. 1 .And.
	2                (menu_page(k) .Eq. 2 .Or. (menu_page(k) .Eq. 4 .And.
	3                                      menu_element(k) .Le. 16)) ) Then
	               dwellbuff(j) = 1
	               y(j,k+1) = fac_no_data
	            Else If ( dwell_mode(j,2) .Eq. 1 .And.
	2                (menu_page(k) .Eq. 3 .Or. (menu_page(k) .Eq. 4 .And.
	3                                      menu_element(k) .Ge. 17)) ) Then
	               dwellbuff(j) = 1
	               y(j,k+1) = fac_no_data
	            End If

	         Else
	            y(j, nfields+2) = fac_no_data
	         End If

	      End Do
   
c Set up plot labels and smooth the data if it was requested.
c
	      If (ismooth .Eq. 1) Then
	         Call FEP_Convolve (y(1,k+1), dwellbuff, nplotpts,
	2                           select_window, k, crossgaps, iwindow)
	      End If

	      If (sylabel(k) .Ne. sylabel(1)) Then
	         same_units = .False.
	      End If

	      ylabel(k) = sylabel(k)
	      m = Index(ylabel(k),'\')
	      If (iwindow .Gt. 1) Then
	         ylabel(k) (m:) = ' (smooth=xxxx)'
	         Write (smooth_label,20) iwindow
20	         Format (i4)
	         m = Index(ylabel(k),'xxxx')
	         ylabel(k) (m:m+3) = smooth_label
	      Else
	         ylabel(k) (m:) = ' '
	      End If

	   End Do
        Do k = 1, nplotpts
            
           if (dwell_present) then
            if (telm_mode(k) .eq. 1) then 
              if (.not. telm_present) telm_present = .true.
              y(k,nfields +3) = 1.
            else
              y(k,nfields + 3) = fac_no_data
            endif 
           endif    
           if (.not. dwell_present) then
            if (telm_mode(k) .eq. 1) then 
              if (.not. telm_present) telm_present = .true.
              y(k,nfields +2) = 1.
            else
              y(k,nfields + 2) = fac_no_data
            endif 
           endif    
           Enddo
c
	   ismooth = 0
	   select_window = 1

	   If (fcc_monitor .Eq. fac_not_present) Then
	      Call LIB$Erase_Page ( 1, 1 )
	   End If

c
c Make the plot...
c
c     Find the left and right plot boundaries and curve-labelling boundaries.
c
	   timeplotmin = - (fac_margin * maxtime) / time_unit
	   timeplotmax = (1. + fac_margin)* maxtime / time_unit

	   xleftlabel = - (0.5*fac_margin * maxtime) / time_unit
	   xrightlabel = (1. + 0.5*fac_margin)* maxtime / time_unit

	   moreplot = .True.

c
c Set up the PLT calling arguments, y(), iery(), nveci, cmd(), and ncmd
c     and find the top and bottom plot boundaries.
c
	   j = 1
	   Do While (abs(y(j,2) - fac_no_data) .Le. 1.e-35)
	      j = j + 1
	   End Do
	   ymin = y(j,2)
	   ymax = y(j,2)

	   Do k=1,nfields
	      ymn(k) = y(1,k+1)
	      ymx(k) = y(1,k+1)
	      Do j=1,nplotpts
	         If (abs(y(j,k+1) - fac_no_data) .Gt. 1.e-35) Then
	            ymn(k) = amin1(ymn(k), y(j,k+1))
	            ymx(k) = amax1(ymx(k), y(j,k+1))
	         End If
	      End Do
	      ymin = amin1(ymin, ymn(k))
	      ymax = amax1(ymax, ymx(k))

	      If (ymn(k) .Eq. ymx(k)) Then
	         If (ymn(k) .Gt. 0.) Then
	            yplotmn(k) = (1. - fac_margin)*ymn(k)
	            yplotmx(k) = (1. + fac_margin)*ymx(k)
	         Else if (ymn(k) .Eq. 0.) Then
	            yplotmn(k) = -1.
	            yplotmx(k) =  1.
	         Else
	            yplotmn(k) = (1. + fac_margin)*ymn(k)
	            yplotmx(k) = (1. - fac_margin)*ymx(k)
	         End If
	      Else
	         yplotmn(k) = ymn(k) - fac_margin*(ymx(k) - ymn(k))
	         yplotmx(k) = ymx(k) + fac_margin*(ymx(k) - ymn(k))
	      End If
	   End Do

	   If (ymin .Eq. ymax) Then
	      If (ymin .Gt. 0.) Then
	         yplotmin = (1. - fac_margin)*ymin
	         yplotmax = (1. + fac_margin)*ymax
	      Else if (ymin .Eq. 0.) Then
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
c
c Scale dwell mode field.
c
           
	   Telmval = ymax + (fac_margin - fac_telm_pad)*(ymax - ymin)
	   dwellval = ymax + (fac_margin - fac_dwell_pad)*(ymax - ymin)
           if (ymax .eq. ymin) then
             telmval = yplotmx(1) -( (yplotmx(1) - ymx(1))*fac_margin)
           endif
	   Do i=1,nplotpts
	      If (dwell_present) then
               If (y(i, nfields+2) .Ne. fac_no_data) Then
	         y(i, nfields+2) = dwellval
                 if( (y(i-1,nfields+2) .eq. fac_no_data) .and.
	1           (y(i+1,nfields+2) .eq. fac_no_data)) then
                   if ( i .ne. 1) then 
                     y(i-1,nfields+2) = dwellval
                   else
                     y(i+1,nfields+2) = dwellval
                   endif 
                 endif 
	       Endif
              if (telm_present) then
	       If (y(i, nfields+3) .Ne. fac_no_data) Then
	         y(i, nfields+3) = telmval
                 if( (y(i-1,nfields+3) .eq. fac_no_data) .and.
	1           (y(i+1,nfields+3) .eq. fac_no_data)) then
                   if ( i .ne. 1) then 
                     y(i-1,nfields+3) = telmval
                   else
                     y(i+1, nfields+3)= telmval 
                   endif   
                 endif 
	       End If
              endif
             else
               if (telm_present) then
                 If (y(i, nfields+2) .Ne. fac_no_data) Then
	           y(i, nfields+2) = telmval
                 if( (y(i-1,nfields+2) .eq. fac_no_data) .and.
	1           (y(i+1,nfields+2) .eq. fac_no_data)) then
                   if (i .ne. 1) then
                     y(i-1,nfields+2) = telmval
                   else
                    y(i+1, nfields+2) = telmval
                   endif  
                 endif 
	         End If
               endif
            endif      
	   End Do
c
c
c
	   Call LIB$Erase_Page ( 1, 1 )
	   If (nplotpts .Eq. 1) Then
	      nplotpts = 2
	      y(2,1) =  maxtime / time_unit
	      Do k=1,nfields
	         y(j,k+1) = -1.
	      End Do
	   End If
c
c If Dwell data, append ', DWELL' (= 7 chars) to "title"; add 7 to "len".
c    If "title" still fits into the one- or two-line limit, augment "num1"
c    or "num2" by 7, too.
c
	   If (dwell_present) Then
C	      title(len+1:) = ', DWELL'
	      nveci = nfields + 2
	      len = len + 7
              if (telm_present) nveci = nveci + 1
           Endif
           if (.not. dwell_present .and. telm_present) then
              nveci = nfields + 2
           endif	   
          
             if (.not. telm_present .and. .not. dwell_present) then 
	      nveci = nfields + 1
	     End If

	   If (len .Le. 66) Then
	      num1 = len
	      num2 = len
	   Else If (len-num1-1 .Le. 66) Then
	      num2 = len
	   End If

	   title1 = title(1:num1)
	   title2 = title(num1+2:num2)
	   title3 = title(num2+2:)

	   Do k=1,nveci
	      iery(k) = 0.
	   End Do
	   cmd(1) = 'D ' // fcc_plot_device
	   cmd(2) = 'LA OT ' // title1
	   cmd(3) = 'LA T '  // title2
	   cmd(4) = 'LA F '  // title3
	   cmd(5) = 'LA X '  // xlabel
c
c When fields are not all in the same units, a stacked plot will be made, so
c   suppress the outer y-label sent to cmd(6).
c
	   If (same_units) Then
	      cmd(6) = 'LA OY ' // ylabel(1)
	   Else
	      cmd(6) = 'LA OY '
	   End If

c   Changed size to 1.0 for both modes. 09/21/88.
!	   If (fcc_interactive .Eq. fac_present) Then
!	      cmd(7) = 'CS 1.5'
!	   Else
	      cmd(7) = 'CS 1'
!	   End If

	   cmd(8) = 'V .2 .1 .9 .8'

	   write (cmd(9), '(a, 4(1x, e16.6))')
	2         'R', timeplotmin, timeplotmax, yplotmin, yplotmax
	   ncmd = 9

	   Do k = 1, nfields
	      write (cmd(ncmd+k), '(a4, i1, a31)')
	2      'LA Y', k+1, ' ' // ylabel(k)
	   End Do
	   ncmd = ncmd + nfields

	   Do k = 1, nfields
	      write ( cmd(ncmd+k), '(a3, i1, 2(1x, e16.6))' )
	2      'R Y', k+1, yplotmn(k), yplotmx(k)
	   End Do
	   ncmd = ncmd + nfields

	   If (dwell_present) Then
	      write (cmd(ncmd+1), '(a4, i1, a9)')
	2      'LA Y', nfields+2, ' ' // 'Dwell on'
	      write ( cmd(ncmd+2), '(a3, i1, 2(1x, e16.6))' )
	2      'R Y', nfields+2, 0., yplotmax
	      ncmd = ncmd + 2
	   End If
C	   If (telm_present) Then
C	      write (cmd(ncmd+1), '(a4, i1, a9)')
C	2      'LA Y', nfields+3, ' ' // 'Telm on'
C	      write ( cmd(ncmd+2), '(a3, i1, 2(1x, e16.6))' )
C	2      'R Y', nfields+3, 0., yplotmax
C	      ncmd = ncmd + 2
C	   End If
c
c On Overlay plots, distinguish the various traces by LineStyle, cycling 1 to 5.
c   Also place numeric LAbels at the left and right endpoints of each trace.
c
	   If (same_units) Then
	      Do k = 1, nfields
	        ncmd = ncmd + 1
	        write (cmd(ncmd), '(a3, i1, a4, i1)')
	2       'LS ', mod(k-1, 5) + 1, ' ON ', k+1
	      End Do
	      m = min(nfields,5)
	      If (nfields .Gt. 1) Then
	        Do k = 1, m
	          found = .false.
	          kk = 1
	          Do While (.Not. found)
	            If (y(kk, k+1).eq. fac_no_data) Then
	              kk = kk + 1
	            Else
	              found = .true.
	            End If
	          End Do
	          kkl = kk
	          found = .false.
	          kk = nplotpts
	          Do While (.Not. found)
	            If (y(kk, k+1).eq. fac_no_data) Then
	              kk = kk - 1
	            Else
	              found = .true.
	            End If
	          End Do
	          kkr = kk
	          ncmd = ncmd + 2
	          write ( cmd(ncmd-1), '(a3, i2, a3, 2(e16.6, 1x), a1, i1, a1)')
	2         'LA ', 2*k-1, ' P ', xleftlabel, y(kkl,k+1), '"', k, '"'
	          write ( cmd(ncmd), '(a3, i2, a3, 2(e16.6, 1x), a1, i1, a1)' )
	2         'LA ', 2*k, ' P ', xrightlabel, y(kkr,k+1), '"', k, '"'
	        End Do
	      End If

	      If (dwell_present .And. m .Lt. 5) Then
	        ncmd = ncmd + 2
	        write ( cmd(ncmd-1), '(a3, i2, a3, 2(e16.6, 1x), a4)' )
	2       'LA ', 2*k-1, ' P ', xleftlabel, dwellval, '"DW"'
	        write ( cmd(ncmd), '(a3, i2, a3, 2(e16.6, 1x), a4)' )
	2       'LA ', 2*k, ' P ', xrightlabel, dwellval, '"DW"'
	      End If
              if (dwell_present) k = k + 1
	      If (telm_present .And. m .Lt. 5) Then
	        ncmd = ncmd + 2
	        write ( cmd(ncmd-1), '(a3, i2, a3, 2(e16.6, 1x), a4)' )
	2       'LA ', 2*k-1, ' P ', xleftlabel, telmval, '"TQ"'
	        write ( cmd(ncmd), '(a3, i2, a3, 2(e16.6, 1x), a4)' )
	2       'LA ', 2*k, ' P ', xrightlabel, telmval, '"TQ"'
	      End If

	      If (nfields .Gt. 1) Then
	        Legend = 'Fields: ' //
	2              '1=solid; 2=dash; 3=dot_dash; 4=dot; 5=dot_dot_dot_dash'
	        If (nfields .Eq. 2) Then
	          len = 23
	        Else If (nfields .Eq. 3) Then
	          len = 35
	        Else If (nfields .Eq. 4) Then
	          len = 42
	        Else
	          len = 62
	        End If

	        ncmd = ncmd + 1
	        cmd(ncmd) = 'LA OX ' // Legend(1:len)
	      End If
	   End If

	   If (same_units) Then
              ncmd = ncmd + 1
	      cmd(ncmd) = 'PLOT OVERLAY'
	   Else
              k = nfields + 1

       	      If (telm_present ) Then
                nveci = nveci - 1
	        ncmd = ncmd + 2
	        write ( cmd(ncmd-1), '(a3, i2, a3, 2(e16.6, 1x), a4)' )
	2       'LA ', 2*k-1, ' P ', xleftlabel, telmval, '"TQ"'
	        write ( cmd(ncmd), '(a3, i2, a3, 2(e16.6, 1x), a4)' )
	2       'LA ', 2*k, ' P ', xrightlabel, telmval, '"TQ"'
	      End If
              
              ncmd = ncmd + 1   
              cmd(ncmd) = 'PLOT VERTICAL'
	 
	   End If

	   If (fcc_plt_com .Eq. fac_present) Then
	      ncmd = ncmd + 1
	      cmd(ncmd)   = '@' // Fcc_plt_com_file
	   End If

	   If (fcc_interactive .Ne. fac_present) Then
	      ncmd = ncmd + 2
	      cmd(ncmd-1) = 'p'
	      cmd(ncmd)   = 'q'
	   End If
c
c Display the plot until the user finishes looking at it.
c
	   Call PLT (y, iery, fac_max_hskp_recs, nplotpts, nveci,
	2            cmd, ncmd, ier)

c
c Interactively select smoothing if not in /MONITOR mode.
c In /MONITOR mode, update the plot after timing out 32 seconds.
c No smoothing prompt in batch mode.
c
	   If (fcc_interactive .Eq. fac_present) Then
	      If (fcc_monitor .Eq. fac_not_present) Then
	         Type 40
40	         Format (1x, 'Choose smooth (S); ',/,
	2                1x, 'Smooth, but Not across time gaps (SN); ',/,
	3                1x, 'or return to menu (<CR>): '$)
	         Accept 80, ans
80	         Format (a)
	      Else
	         Call LIB$Wait(32.)
	      End If
	   Else
	      ans = ' '
	   End If

c
c Select the smooth window (interactive only).
c
	   If (ans(1:1) .Eq. 'S' .Or. ans(1:1) .Eq. 's') Then
	      ismooth = 1
	      select_window = 1
	      moreplot = .True.
	      If (ans(2:2) .Eq. 'N' .Or. ans(2:2) .Eq. 'n') Then
	         crossgaps = .False.
	      Else
	         crossgaps = .True.
	      End If
	   Else
	      moreplot = .False.
	   End If

	End Do

	If (fcc_monitor .Eq. fac_not_present) Then
	   Call LIB$Erase_Page ( 1, 1 )
	End If


	FEP_Plot = %Loc(FEP_Normal)

	Return
	End
