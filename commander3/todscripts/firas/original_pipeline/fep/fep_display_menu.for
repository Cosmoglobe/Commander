	Integer*4 Function FEP_Display_Menu ( menu_page, menu_element,
	2                                     convert_flag, digital_mask,
	3                                     nfields )

C------------------------------------------------------------------------
C    PURPOSE: Displays the housekeeping record fields in menu form and
C	      allows the user to select fields to be plotted.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX
C            March 26, 1987
C
C    INVOCATION: STATUS = FEP_DISPLAY_MENU ( MENU_PAGE, MENU_ELEMENT,
C					     CONVERT_FLAG, DIGITAL_MASK,
C					     NFIELDS )
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS: 
C	MENU_PAGE(*)		I*2		Selected menu page.		
C	MENU_ELEMENT(*)		I*2		Selected menu element.
C	CONVERT_FLAG(*)		L*1		[No]Convert selected field.
C	DIGITAL_MASK(*)		I*2		Mask for digital status flds.
C	NFIELDS			I*2		Number of fields.
C
C    SUBROUTINES CALLED: 
C	FEP_Display_Main
C	FEP_Display_Plot_Menu
C	SMG$Create_PasteBoard
C	SMG$Erase_PasteBoard
C	SMG$Delete_PasteBoard
C
C    COMMON VARIABLES USED:
C	HSKP_MENU(28)		CH*28		Housekeeping header menu.
C	HSKP_DB(28)		CH*8		STOL database words.
C	GRTA_MENU(32)		CH*32		GRT side A menu.
C	GRTA_DB(32)		CH*8		STOL database words.
C	GRTB_MENU(32)		CH*32		GRT side B menu.
C	GRTB_DB(32)		CH*8		STOL database words.
C	COMB_MENU(32)		CH*32		Combined GRT menu.
C	IPDU_MENU(14)		CH*32		IPDU temps menu.
C	IPDU_DB(14)		CH*8		STOL database words.
C	VOLT_MENU(20)		CH*32		Voltages menu.
C	VOLT_DB(20)		CH*8		STOL database words.
C	CURR_MENU(18)		CH*32		Currents menu.
C	CURR_DB(18)		CH*8		STOL database words.
C
C    INCLUDE FILES: 
C	FUT_Error
C	FEP_Menu
C	$SMGDEF
C	$SSDEF
C
C----------------------------------------------------------------------
C
C Changes:
C
C	SPR 1820, Use parameter FAC_MAX_MENU_NUM to maintain menu pages.
C	R. Kummerer,  December 11, 1987.
C
C	SPR 2377, Provide masking capability for digital status words.
C	Fred Shuman,  1989 Feb 20.
C
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FUT_Params)'
	Include 	'(FEP_Menu)'
	Include 	'($SMGDEF)'
	Include 	'($SSDEF)'

	Integer 	*4	i
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
	Logical		*1	convert_flag(*)
	Integer		*2	digital_mask(*)

	Logical		*1	exit_menu
	Logical		*1	quit_menu

	Integer		*4	FEP_Display_Main
	Integer		*4	FEP_Display_Plot_Menu
	Integer		*4	SMG$Create_PasteBoard
	Integer		*4	SMG$Erase_PasteBoard
	Integer		*4	SMG$Delete_PasteBoard

	External	FEP_Normal
	External	FEP_ExitMenu
	External	FEP_QuitMenu

	Common /PB/ 	pbid, irow, icol
c
c Initialize.
c
	status = SMG$Create_PasteBoard ( pbid, 'SYS$OUTPUT', irow,
     .						icol, )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = %Loc(FEP_Normal)
	exit_menu = .False.
	quit_menu = .False.

	Do While ((.Not. exit_menu .And. .Not. quit_menu) .And.
     .		  status .Eq. %Loc(FEP_Normal))

c
c	  Select from the main menu.
c
	   status = FEP_Display_Main ( inp, menu_page, menu_element,
     .					convert_flag, digital_mask, nfields )
	   If (status .Ne. %Loc(FEP_Normal)) Then
	      FEP_Display_Menu = status
	      Return
	   End If

	   If (inp .Eq. HSKP_page) Then

c
c	  Select from the housekeeping header menu.
c
	      status = FEP_Display_Plot_Menu ( HSKP_page,
     .				HSKP_menu, HSKP_menu_num,
     .				menu_page, menu_element,
     .				convert_flag, digital_mask, nfields )
	      If (status .Ne. %Loc(FEP_Normal)) Then
	         FEP_Display_Menu = status
	         Return
	      End If

	   Else If (inp .Eq. GRTA_page) Then
c
c	  Select from the A side GRTs menu.
c
	      status = FEP_Display_Plot_Menu ( GRTA_page,
     .				GRTA_menu, GRTA_menu_num,
     .				menu_page, menu_element,
     .				convert_flag, digital_mask, nfields )
	      If (status .Ne. %Loc(FEP_Normal)) Then
	         FEP_Display_Menu = status
	         Return
	      End If

	   Else If (inp .Eq. GRTB_page) Then

c
c	  Select from the B side GRTs menu.
c
	      status = FEP_Display_Plot_Menu ( GRTB_page,
     .				GRTB_menu, GRTB_menu_num,
     .				menu_page, menu_element,
     .				convert_flag, digital_mask, nfields )
	      If (status .Ne. %Loc(FEP_Normal)) Then
	         FEP_Display_Menu = status
	         Return
	      End If

	   Else If (inp .Eq. COMB_page) Then

c
c	  Select from the Combined GRTs menu.
c
	      status = FEP_Display_Plot_Menu ( COMB_page,
     .				COMB_menu, COMB_menu_num,
     .				menu_page, menu_element,
     .				convert_flag, digital_mask, nfields )
	      If (status .Ne. %Loc(FEP_Normal)) Then
	         FEP_Display_Menu = status
	         Return
	      End If

	   Else If (inp .Eq. IPDU_page) Then

c
c	  Select from the IPDU temperatures menu.
c
	      status = FEP_Display_Plot_Menu ( IPDU_page,
     .				IPDU_menu, IPDU_menu_num,
     .				menu_page, menu_element,
     .				convert_flag, digital_mask, nfields )
	      If (status .Ne. %Loc(FEP_Normal)) Then
	         FEP_Display_Menu = status
	         Return
	      End If

	   Else If (inp .Eq. VOLT_page) Then

c
c	  Select from the voltages menu.
c
	      status = FEP_Display_Plot_Menu ( VOLT_page,
     .				VOLT_menu, VOLT_menu_num,
     .				menu_page, menu_element,
     .				convert_flag, digital_mask, nfields )
	      If (status .Ne. %Loc(FEP_Normal)) Then
	         FEP_Display_Menu = status
	         Return
	      End If

	   Else If (inp .Eq. CURR_page) Then

c
c	  Select from the currents menu.
c
	      status = FEP_Display_Plot_Menu ( CURR_page,
     .				CURR_menu, CURR_menu_num,
     .				menu_page, menu_element,
     .				convert_flag, digital_mask, nfields )
	      If (status .Ne. %Loc(FEP_Normal)) Then
	         FEP_Display_Menu = status
	         Return
	      End If

	   Else If (inp .Eq. 8) Then

c
c	  Reset the field plot selection.
c
	      nfields = 0
	      Do i=1,fac_max_menu_num
		 menu_page(i) = 0
		 menu_element(i) = 0
		 convert_flag(i) = .False.
	      End Do

	   Else If (inp .Eq. 9) Then

c
c	  Exit and plot.
c
	      exit_menu = .True.

	   Else

c
c	  Quit and plot.
c
	      quit_menu = .True.

	   End If

	End Do

	status = SMG$Erase_PasteBoard ( pbid )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Delete_PasteBoard ( pbid, )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If


	If (exit_menu) Then
	   FEP_Display_Menu = %Loc(FEP_ExitMenu)
	Else
	   FEP_Display_Menu = %Loc(FEP_QuitMenu)
	End If

	Return
	End
