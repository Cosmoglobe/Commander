	Integer*4 Function FIT_Display_Menu ( menu_page, menu_element,
	2                                     bin_size, nfields )

C------------------------------------------------------------------------
C    PURPOSE: Displays the engineering statistics record fields in menu
C	      form and allows the user to select fields to be plotted.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Fred Shuman, STX	(adapted from FEP_Display_Menu
C            October 1, 1987	 by Rob Kummerer, STX, March 26, 1987)
C            version: Nov 21, 1987
C
C    INVOCATION: STATUS = FIT_DISPLAY_MENU ( MENU_PAGE, MENU_ELEMENT,
C					     BIN_SIZE, NFIELDS )
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS:
C	MENU_PAGE(*)		I*2		Selected menu page.
C	MENU_ELEMENT(*)		I*2		Selected menu element.
C	BIN_SIZE		R*4		Selected bin size (seconds).
C	NFIELDS			I*2		Number of fields.
C
C    SUBROUTINES CALLED:
C	FIT_Display_Main
C	FIT_Display_Plot_Menu
C	SMG$Create_PasteBoard
C	SMG$Delete_PasteBoard
C	SMG$Erase_PasteBoard
C	LIB$Signal
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
C	              FIT_Display_Main were revised.
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
	Integer		*4	status
	Integer		*4	inp
	Integer		*4	pbid
	Integer		*4	vdid
	Integer		*4	kbid
	Integer		*4	irow
	Integer		*4	icol

	Integer		*2	nfields
	Integer		*2	menu_page(*)
	Integer		*2	menu_element(*)
	Real		*4	bin_size
	Character	*32	binstring(24)
	Character	*32	selprmt

	Logical		*1	exit_menu
	Logical		*1	quit_menu

	Integer		*4	FIT_Display_Main
	Integer		*4	FIT_Display_Plot_Menu
	Integer		*4	SMG$Create_PasteBoard
	Integer		*4	SMG$Delete_PasteBoard
	Integer		*4	SMG$Erase_PasteBoard

	External	FIT_Normal
	External	FIT_ExitMenu
	External	FIT_QuitMenu

	Common /PB/	pbid, irow, icol
c
c Initialize.
c
	status = SMG$Create_PasteBoard ( pbid, 'SYS$OUTPUT', irow,
     .						icol, )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = %Loc(FIT_Normal)
	exit_menu = .False.
	quit_menu = .False.

	bin_size = fcc_orbit
	binstring(1) = '1 ORBIT'
	Do i=2,24
	   binstring(i) = ' '
	End Do

	Do While ((.Not. exit_menu .And. .Not. quit_menu) .And.
     .		  status .Eq. %Loc(FIT_Normal))

c
c  Select from the main menu.
c
	   status = FIT_Display_Main ( inp, menu_page, menu_element,
     .					bin_size, binstring, nfields )
	   If (status .Ne. %Loc(FIT_Normal)) Then
	      FIT_Display_Menu = status
	      Return
	   End If

	   If (inp .Eq. GRTA_page) Then
c
c  Select from the A side GRTs menu.
c
	      status = FIT_Display_Plot_Menu ( GRTA_page,
     .				GRTA_menu, GRTA_menu_num,
     .				menu_page, menu_element,
     .				nfields )
	      If (status .Ne. %Loc(FIT_Normal)) Then
	         FIT_Display_Menu = status
	         Return
	      End If

	   Else If (inp .Eq. GRTB_page) Then

c
c  Select from the B side GRTs menu.
c
	      status = FIT_Display_Plot_Menu ( GRTB_page,
     .				GRTB_menu, GRTB_menu_num,
     .				menu_page, menu_element,
     .				nfields )
	      If (status .Ne. %Loc(FIT_Normal)) Then
	         FIT_Display_Menu = status
	         Return
	      End If

	   Else If (inp .Eq. COMB_page) Then

c
c  Select from the Combined GRTs menu.
c
	      status = FIT_Display_Plot_Menu ( COMB_page,
     .				COMB_menu, COMB_menu_num,
     .				menu_page, menu_element,
     .				nfields )
	      If (status .Ne. %Loc(FIT_Normal)) Then
	         FIT_Display_Menu = status
	         Return
	      End If

	   Else If (inp .Eq. IPDU_page) Then

c
c  Select from the IPDU temperatures menu.
c
	      status = FIT_Display_Plot_Menu ( IPDU_page,
     .				IPDU_menu, IPDU_menu_num,
     .				menu_page, menu_element,
     .				nfields )
	      If (status .Ne. %Loc(FIT_Normal)) Then
	         FIT_Display_Menu = status
	         Return
	      End If

	   Else If (inp .Eq. VOLT_page) Then

c
c  Select from the voltages menu.
c
	      status = FIT_Display_Plot_Menu ( VOLT_page,
     .				VOLT_menu, VOLT_menu_num,
     .				menu_page, menu_element,
     .				nfields )
	      If (status .Ne. %Loc(FIT_Normal)) Then
	         FIT_Display_Menu = status
	         Return
	      End If

	   Else If (inp .Eq. CURR_page) Then

c
c  Select from the currents menu.
c
	      status = FIT_Display_Plot_Menu ( CURR_page,
     .				CURR_menu, CURR_menu_num,
     .				menu_page, menu_element,
     .				nfields )
	      If (status .Ne. %Loc(FIT_Normal)) Then
	         FIT_Display_Menu = status
	         Return
	      End If

	   Else If (inp .Eq. 7) Then

c
c  Reset the field plot selection.
c
	      nfields = 0
	      Do i=1,total_menu_num
	         menu_page(i) = 0
	         menu_element(i) = 0
	      End Do
	      binstring(1) = '1 ORBIT'
	      Do i=2,24
	         binstring(i) = ' '
	      End Do

	   Else If (inp .Eq. 13) Then
c
c  Exit and plot.
c
	      exit_menu = .True.

	   Else If (inp .Eq. 14) Then
c
c  Quit (without plotting).
c
	      quit_menu = .True.

	   End If

	End Do

	status = SMG$Erase_PasteBoard ( pbid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Delete_PasteBoard ( pbid, )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If


	If (exit_menu) Then
	   FIT_Display_Menu = %Loc(FIT_ExitMenu)
	Else
	   FIT_Display_Menu = %Loc(FIT_QuitMenu)
	End If

	Return
	End
