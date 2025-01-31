	Integer*4 Function FIT_Display_Main ( inp, menu_page,
	2                           menu_element, bin_size, binstring, nfields )

C------------------------------------------------------------------------
C    PURPOSE: Displays the main menu from which the field menus are
C	      selected.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Fred Shuman, STX	(adapted from FEP_Display_Main
C            October 1, 1987	 by Rob Kummerer, STX, March 26, 1987)
C            version: Nov 2, 1987
C
C    INVOCATION: STATUS = FIT_DISPLAY_MAIN ( INP, MENU_PAGE, MENU_ELEMENT,
C					     BIN_SIZE, BINSTRING,
C					     NFIELDS )
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS:
C	INP			I*4		Fields menu selection.
C	MENU_PAGE(*)		I*2		Selected menu page.
C	MENU_ELEMENT(*)		I*2		Selected menu element.
C	BIN_SIZE		R*4		Selected bin size (seconds).
C	BINSTRING(*)		C*32		Selected bin display string.
C	NFIELDS			I*2		Number of fields.
C
C    SUBROUTINES CALLED:
C	STR$UpCase
C	OTS$Cvt_TI_L
C	OTS$Cvt_T_F
C	SMG$Create_Virtual_Display
C	SMG$Create_Virtual_Keyboard
C	SMG$Delete_Virtual_Display
C	SMG$Delete_Virtual_Keyboard
C	SMG$Paste_Virtual_Display
C	SMG$Draw_Rectangle
C	SMG$Draw_Line
C	SMG$Put_Line
C	SMG$Ring_Bell
C	SMG$Flush_Buffer
C	SMG$Read_String
C	SMG$Set_Cursor_Abs
C	SMG$Erase_Chars
C	SMG$Erase_PasteBoard
C
C    COMMON VARIABLES USED:
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
C	FIT_Menu
C	$SMGDEF
C	$SSDEF
C
C----------------------------------------------------------------------
C	Revised:
C
C	   SPR 3285.  Get COBE orbital period from FUT_Orbital_Period.
C	              This module and FIT, FIT_Get_Command, and
C	              FIT_Display_Menu were revised.
C	              Fred Shuman, STX  1989 Jul 20.
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FUT_Params)'
	Include		'(FIT_Menu)'
	Include		'(FIT_Invoc)'
	Include		'($SMGDEF)'
	Include		'($SSDEF)'

	Integer		*4	i
	Integer		*4	k
	Integer		*4	status
	Integer		*4	inp
	Character	*4	cp
	Integer		*4	irow
	Integer		*4	icol
	Integer		*4	pbid
	Integer		*4	vdid
	Integer		*4	kbid
	Character	*38	bin_hdr
	Character	*22	indent
	Character	*38	MAIN_items(14)		! 14 = Number of
	Integer		*4	MAIN_items_num/14/	!      main menu items
	Character	*80	MAIN_lines(10)		! 10 = Number of
	Integer		*4	MAIN_lines_num/10/	!      main menu lines
	Character	*12	prompt
	Integer		*2	m
	Integer		*2	n
	Character	*32	selprmt
	Logical		*1	display_menu
	Logical		*1	bin_prompt

	Integer		*2	nfields
	Character	*32	field
	Integer		*2	menu_page(*)
	Integer		*2	menu_element(*)
	Real		*4	bin_size
	Character	*32	binstring(*)
	Character	*7	duration
	Real		*4	dur

	Integer		*4	STR$UpCase
	Integer		*4	OTS$Cvt_TI_L
	Integer		*4	OTS$Cvt_T_F
	Integer		*4	SMG$Create_Virtual_Display
	Integer		*4	SMG$Create_Virtual_Keyboard
	Integer		*4	SMG$Delete_Virtual_Display
	Integer		*4	SMG$Delete_Virtual_Keyboard
	Integer		*4	SMG$Paste_Virtual_Display
	Integer		*4	SMG$Draw_Rectangle
	Integer		*4	SMG$Draw_Line
	Integer		*4	SMG$Put_Line
	Integer		*4	SMG$Ring_Bell
	Integer		*4	SMG$Flush_Buffer
	Integer		*4	SMG$Read_String
	Integer		*4	SMG$Erase_Chars
	Integer		*4	SMG$Set_Cursor_Abs
	Integer		*4	SMG$Erase_PasteBoard

	External	FIT_Normal

	Data MAIN_items /' 1 A Side GRT Menu       ',
	2                ' 2 B Side GRT Menu       ',
	3                ' 3 GRTs Combined Menu    ',
	4                ' 4 IPDU Temperatures Menu',
	5                ' 5 Voltages Menu         ',
	6                ' 6 Currents Menu         ',
	7                ' 7 Reset plots selection ',
	8                ' 8 ORBIT (how many--default = 1.0)',
	9                ' 9 DAILY',
	1                '10 MONTHLY (30 days)',
	1                '11 TOTAL MISSION',
	2                '12 Length of bin as DDDHHMM',
	3                '13 Plot the chosen fields',
	4                '14 Quit Menus; get a new timerange'/
	Data bin_hdr    /'Bin size (default = 1 orbit):'/
	Data indent     /' '/           !Actual depth of indent is controlled
C                                       ! not here, but in declared length
	Common /PB/	pbid, irow, icol
c
c Set up for drawing the main menu on the screen.
c
	status = SMG$Erase_PasteBoard ( pbid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Create_Virtual_Display ( irow, icol, vdid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Create_Virtual_Keyboard ( kbid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Set_Cursor_Abs ( vdid, 2, 35 )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If
c
c"                      Main Menu                       "
c"______________________________________________________"
c
	status = SMG$Put_Line ( vdid, 'Main Menu' )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Draw_Line ( vdid, 3, 1, 3, icol )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If
C
C Construct the Main menu selection lines
C
	MAIN_lines(1) = MAIN_items(1) // bin_hdr
	MAIN_lines(2) = MAIN_items(2)
	MAIN_lines(3) = MAIN_items(3) // '  ' // MAIN_items(8)
	MAIN_lines(4) = MAIN_items(4) // '  ' // MAIN_items(9)
	MAIN_lines(5) = MAIN_items(5) // '  ' // MAIN_items(10)
	MAIN_lines(6) = MAIN_items(6) // '  ' // MAIN_items(11)
	MAIN_lines(7) = MAIN_items(7) // '  ' // MAIN_items(12)
	MAIN_lines(8) = ' '
	MAIN_lines(9) = indent // MAIN_items(13)
	MAIN_lines(10)= indent // MAIN_items(14)
C
C Write them on the screen.
C
	Do i=1,MAIN_lines_num

	   status = SMG$Set_Cursor_Abs ( vdid, i+3, 3 )
	   If (status .Ne. SS$_Normal) Then
	      FIT_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   status = SMG$Put_Line ( vdid, MAIN_lines(i) )
	   If (status .Ne. SS$_Normal) Then
	      FIT_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	End Do
C
C Write the selected plot fields on the screen.
C
	status = SMG$Set_Cursor_Abs ( vdid, irow-5, 3 )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Put_Line ( vdid,
	2               'Selected Plot Fields:                 '//
	3               'Selected Bin Size:' )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Set_Cursor_Abs ( vdid, irow-4, 44 )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Put_Line ( vdid, binstring(1))
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	Do i = 1,nfields

	   status = SMG$Set_Cursor_Abs ( vdid, irow+i-5, 6 )
	   If (status .Ne. SS$_Normal) Then
	      FIT_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   If (menu_page(i) .Eq. GRTA_page) Then
	      field = GRTA_menu(menu_element(i))
	   Else If (menu_page(i) .Eq. GRTB_page) Then
	      field = GRTB_menu(menu_element(i))
	   Else If (menu_page(i) .Eq. COMB_page) Then
	      field = COMB_menu(menu_element(i))
	   Else If (menu_page(i) .Eq. IPDU_page) Then
	      field = IPDU_menu(menu_element(i))
	   Else If (menu_page(i) .Eq. VOLT_page) Then
	      field = VOLT_menu(menu_element(i))
	   Else If (menu_page(i) .Eq. CURR_page) Then
	      field = CURR_menu(menu_element(i))
	   End If

	   status = SMG$Put_Line ( vdid, field//'      '//binstring(i))
	   If (status .Ne. SS$_Normal) Then
	      FIT_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	End Do

	status = SMG$Draw_Line ( vdid, irow-6, 1, irow-6, icol )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If
c
c Draw the Main menu frame on the screen.
c
	status = SMG$Draw_Rectangle ( vdid, 1, 1, irow, icol )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Paste_Virtual_Display ( vdid, pbid, 1, 1 )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Flush_Buffer ( pbid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	display_menu = .True.

	Do While (display_menu)

	   status = SMG$Erase_Chars ( vdid, 76, irow-9, 3 )
	   If (status .Ne. SS$_Normal) Then
	      FIT_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   status = SMG$Set_Cursor_Abs ( vdid, irow-9, 3 )
	   If (status .Ne. SS$_Normal) Then
	      FIT_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   selprmt = 'Select a menu item: '
	   status = SMG$Read_String ( kbid, cp, selprmt,,,,,,, vdid )
	   If (status .Ne. SS$_Normal) Then
	      FIT_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   status = OTS$Cvt_TI_L(cp,inp,%val(4),%val(17))
c
c Check whether choice# or choice letter was input
c
	   If (inp .Eq. 0) Then         ! A non-numeric was input...

	      status = STR$UpCase ( cp, cp )
	      If (status .Ne. SS$_Normal) Then
	         FIT_Display_Main = status
	         Call LIB$Signal(%Val(status))
	         Return
	      End If
c
c ...so check its initial letter against the available choice letters
c
	      Do i=1,MAIN_items_num
	         If (MAIN_items(i)(4:4) .Eq. cp(1:1)) Then
	            inp = i
	         End If
	      End Do

	   End If
c
c  Now check for bin size selection.
c
	   If (inp .Eq. 8) Then
c
c  Bin size = User's choice (R*4 no. of orbits).
c
	      bin_prompt = .True.

	      Do While (bin_prompt)

	         status = SMG$Erase_Chars ( vdid, 76, irow-9, 3 )
	         If (status .Ne. SS$_Normal) Then
	            FIT_Display_Main = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If

	         selprmt = 'Select a bin size (orbits):'

	         status = SMG$Set_Cursor_Abs ( vdid, irow-9, 3 )
	         If (status .Ne. SS$_Normal) Then
	            FIT_Display_Main = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If

	         status = SMG$Read_String ( kbid, prompt, selprmt,
	2                                   ,,,,,, vdid )
	         If (status .Ne. SS$_Normal) Then
	            FIT_Display_Main = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If

	         status = OTS$Cvt_T_F(prompt, bin_size,,,%val(19))
	         bin_size = abs(bin_size)

	         m = 1
	         Do While (prompt(m:m) .Eq. ' ' .And. m .Le. 12)
	            m = m+1
	         End Do

	         n = 12
	         Do While (prompt(n:n) .Eq. ' ' .And. n .Ge. m)
	            n = n-1
	         End Do

	         If (abs(bin_size-1.) .Lt. 1.E-10) Then
	            binstring(1) = '1 ORBIT'
	         Else
	            binstring(1) = prompt(m:n) // ' ORBITS'
	         End If

	         bin_size = bin_size * fcc_orbit
	         If (bin_size .Lt. fcc_orbit .Or.
	2            bin_size .Gt. fac_mission) Then
	            status = SMG$Ring_Bell ( vdid, 3 )
	         Else
	            bin_prompt = .False.
	         End If

	      End Do
c
c  Bin size = One day.
c
	   Else If (inp .Eq. 9) Then
	      bin_size = fac_day
	      binstring(1) = '1 DAY'
c
c  Bin size = One month.
c
	   Else If (inp .Eq. 10) Then
	      bin_size = fac_month
	      binstring(1) = '30 DAYS'
c
c  Bin size = Mission-to-date.
c
	   Else If (inp .Eq. 11) Then
	      bin_size = fac_mission
	      binstring(1) = 'TOTAL MISSION'
c
c  Bin size = User's choice (DDDHHMM).
c
	   Else If (inp .Eq. 12) Then

	      bin_prompt = .True.

	      Do While (bin_prompt)

	         status = SMG$Erase_Chars ( vdid, 76, irow-9, 3 )
	         If (status .Ne. SS$_Normal) Then
	            FIT_Display_Main = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If

	         selprmt = 'Select a bin size (DDDHHMM): '

	         status = SMG$Set_Cursor_Abs ( vdid, irow-9, 3 )
	         If (status .Ne. SS$_Normal) Then
	            FIT_Display_Main = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If

	         status = SMG$Read_String ( kbid, duration, selprmt,
	2                                   ,,,,,, vdid )
	         If (status .Ne. SS$_Normal) Then
	            FIT_Display_Main = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If

	         status = OTS$Cvt_T_F(duration(1:3), dur,,,%val(19))
	         bin_size = dur*fac_day
	         status = OTS$Cvt_T_F(duration(4:5), dur,,,%val(19))
	         bin_size = bin_size + dur*fac_hour
	         status = OTS$Cvt_T_F(duration(6:7), dur,,,%val(19))
	         bin_size = bin_size + dur*fac_minute

	         If (bin_size .Lt. fcc_orbit .Or.
	2            bin_size .Gt. fac_mission) Then
	            status = SMG$Ring_Bell ( vdid, 3 )
	         Else
	            bin_prompt = .False.
	            binstring(1) = 'DURATION=' // duration
	         End If

	      End Do

	   Else If (inp .Lt. 1 .Or. inp .Gt. MAIN_items_num) Then

	      status = SMG$Ring_Bell ( vdid, 3 )	! input out-of-bounds

	   Else

	      display_menu = .False.	! Matched a submenu# or Reset Plots;
	                                !  set up to exit DOWHILE and Return
	   End If
c
c  Write the bin selection string to the selection portion of the Main Menu
c
	   status = SMG$Set_Cursor_Abs ( vdid, irow-4, 44 )
	   If (status .Ne. SS$_Normal) Then
	      FIT_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   status = SMG$Put_Line ( vdid, binstring(1))
	   If (status .Ne. SS$_Normal) Then
	      FIT_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	End Do

	status = SMG$Erase_PasteBoard ( pbid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Delete_Virtual_Keyboard ( kbid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Delete_Virtual_Display ( vdid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If


	FIT_Display_Main = %Loc(FIT_Normal)

	Return
	End
