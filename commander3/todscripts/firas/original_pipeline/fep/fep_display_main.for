	Integer*4 Function FEP_Display_Main ( inp, menu_page, menu_element,
     .                                 convert_flag, digital_mask, nfields )

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
C    AUTHOR: Rob Kummerer
C            STX
C            March 26, 1987
C
C    INVOCATION: STATUS = FEP_DISPLAY_MAIN ( INP, MENU_PAGE,
C					     MENU_ELEMENT, CONVERT_FLAG,
C					     DIGITAL_MASK, NFIELDS )
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS: 
C	INP			I*4		Fields menu selection.
C	MENU_PAGE(*)		I*2		Selected menu page.		
C	MENU_ELEMENT(*)		I*2		Selected menu element.
C	CONVERT_FLAG(*)		L*1		[No]Convert selected field.
C	DIGITAL_MASK(*)		I*2		Mask for digital status flds.
C	NFIELDS			I*2		Number of fields.
C
C    SUBROUTINES CALLED: 
C	STR$UpCase
C	OTS$Cvt_TI_L
C	SMG$Create_Virtual_Display
C	SMG$Create_Virtual_Keyboard
C	SMG$Delete_Virtual_Keyboard
C	SMG$Paste_Virtual_Display
C	SMG$Draw_Rectangle
C	SMG$Draw_Line
C	SMG$Put_Line
C	SMG$Ring_Bell
C	SMG$Flush_Buffer
C	SMG$Read_String
C	SMG$Set_Cursor_Abs
C	SMG$Set_Cursor_Rel
C	SMG$Create_PasteBoard
C	SMG$Erase_PasteBoard
C	SMG$Delete_PasteBoard
C	SMG$Delete_Virtual_Display
C
C    COMMON VARIABLES USED:
C	HSKP_MENU(28)		CH*28		Housekeeping header menu.
C	GRTA_MENU(32)		CH*32		GRT side A menu.
C	GRTB_MENU(32)		CH*32		GRT side B menu.
C	COMB_MENU(32)		CH*32		Combined GRT menu.
C	IPDU_MENU(14)		CH*32		IPDU temps menu.
C	VOLT_MENU(20)		CH*32		Voltages menu.
C	CURR_MENU(18)		CH*32		Currents menu.
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
C	SPR 2377, Provide masking capability for digital status words.
C	Fred Shuman,  1989 Feb 20.
C
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include 	'(FEP_Menu)'
	Include 	'($SMGDEF)'
	Include 	'($SSDEF)'

	Integer 	*4	i
	Integer 	*4	j
	Integer 	*4	k
	Integer 	*4	first
	Integer 	*4	second
	Integer		*4	status
	Integer		*4	inp
	Character	*4	cp
	Integer		*4	irow
	Integer		*4	icol
	Integer		*4	pbid
	Integer		*4	vdid
	Integer		*4	kbid
	Character	*32	MAIN_menu(10)
	Integer		*4	MAIN_menu_num/10/
	Logical		*1	display_menu

	Integer		*2	nfields
	Character	*64	field
	Character	*12	dum
	Integer		*2	menu_page(*)
	Integer		*2	menu_element(*)
	Logical		*1	convert_flag(*)
	Integer		*2	digital_mask(*)

	Integer		*4	STR$UpCase
	Integer		*4	OTS$Cvt_TI_L
	Integer		*4	SMG$Create_Virtual_Display
	Integer		*4	SMG$Create_Virtual_Keyboard
	Integer		*4	SMG$Delete_Virtual_Keyboard
	Integer		*4	SMG$Paste_Virtual_Display
	Integer		*4	SMG$Draw_Rectangle
	Integer		*4	SMG$Draw_Line
	Integer		*4	SMG$Put_Line
	Integer		*4	SMG$Ring_Bell
	Integer		*4	SMG$Flush_Buffer
	Integer		*4	SMG$Read_String
	Integer		*4	SMG$Set_Cursor_Abs
	Integer		*4	SMG$Set_Cursor_Rel
	Integer		*4	SMG$Erase_PasteBoard
	Integer		*4	SMG$Delete_Virtual_Display

	External	FEP_Normal

	Data MAIN_menu /'1  Housekeeping Header Menu',
     .		        '2  A Side GRT Menu',
     .		        '3  B Side GRT Menu',
     .		        '4  GRTs Combined Menu',
     .		        '5  IPDU Temperatures Menu',
     .		        '6  Voltages Menu',
     .		        '7  Currents Menu',
     .		        '8  Reset plots selection',
     .		        '9  Exit and make plots',
     .		        '10 Quit ENGPLOTS'/

	Common /PB/	pbid, irow, icol

c
c Select from the main menu.
c
	status = SMG$Erase_PasteBoard ( pbid )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Create_Virtual_Display ( irow, icol, vdid )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Create_Virtual_Keyboard ( kbid )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Set_Cursor_Abs ( vdid, 2, 35 )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Put_Line ( vdid, 'Main Menu' )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Draw_Line ( vdid, 3, 1, 3, icol )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	Do i=1,MAIN_menu_num

	   status = SMG$Set_Cursor_Abs ( vdid, i+3, 3 )
	   If (status .Ne. SS$_Normal) Then
	      FEP_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   status = SMG$Put_Line ( vdid, MAIN_menu(i) )
	   If (status .Ne. SS$_Normal) Then
	      FEP_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	End Do

	status = SMG$Set_Cursor_Abs ( vdid, irow-5, 3 )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Put_Line ( vdid, 'Selected Plot Fields:' )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	j = 0
	Do i=irow-4,irow-(4-nfields)-1

	   status = SMG$Set_Cursor_Abs ( vdid, i, 6 )
	   If (status .Ne. SS$_Normal) Then
	      FEP_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   j = j + 1

	   If (menu_page(j) .Eq. HSKP_page) Then
	      field = HSKP_menu(menu_element(j))
	   Else If (menu_page(j) .Eq. GRTA_page) Then
	      field = GRTA_menu(menu_element(j))
	   Else If (menu_page(j) .Eq. GRTB_page) Then
	      field = GRTB_menu(menu_element(j))
	   Else If (menu_page(j) .Eq. COMB_page) Then
	      field = COMB_menu(menu_element(j))
	   Else If (menu_page(j) .Eq. IPDU_page) Then
	      field = IPDU_menu(menu_element(j))
	   Else If (menu_page(j) .Eq. VOLT_page) Then
	      field = VOLT_menu(menu_element(j))
	   Else If (menu_page(j) .Eq. CURR_page) Then
	      field = CURR_menu(menu_element(j))
	   End If

	   If (convert_flag(j)) Then
	      field = 'Convert    :  ' // field
	   Else
	      field = 'No Convert :  ' // field
	   End If

	   If (digital_mask(j) .Ne. 0) Then
	      Write (dum, '(a6, i6)') 'Mask= ', digital_mask(j)
	      field(53:64) = dum
	   End If

	   status = SMG$Put_Line ( vdid, field )
	   If (status .Ne. SS$_Normal) Then
	      FEP_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	End Do

	status = SMG$Draw_Line ( vdid, irow-6, 1, irow-6, icol )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Draw_Rectangle ( vdid, 1, 1, irow, icol )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Paste_Virtual_Display ( vdid, pbid, 1, 1 )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Flush_Buffer ( pbid )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	display_menu = .True.

	Do While (display_menu)

	   status = SMG$Set_Cursor_Abs ( vdid, irow-9, 3 )
	   If (status .Ne. SS$_Normal) Then
	      FEP_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   status = SMG$Read_String ( kbid, cp, 'Select a menu item: ',
     .				      ,,,,,, vdid )
	   If (status .Ne. SS$_Normal) Then
	      FEP_Display_Main = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   status = OTS$Cvt_TI_L(cp,inp,%val(4),%val(17))

	   If (inp .Eq. 0) Then

	      status = STR$UpCase ( cp, cp )
	      If (status .Ne. SS$_Normal) Then
	         FEP_Display_Main = status
	         Call LIB$Signal(%Val(status))
	         Return
	      End If

	      Do i=1,MAIN_menu_num
		 If (MAIN_menu(i)(4:4) .Eq. cp(1:1)) Then
		    inp = i
		 End If
	      End Do

	      If (inp .Eq. 0) Then
		 status = SMG$Ring_Bell ( vdid, 3 )
	         If (status .Ne. SS$_Normal) Then
	            FEP_Display_Main = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If
	      Else
	         display_menu = .False.
	      End If

	   Else If (inp .Lt. 1 .Or. inp .Gt. MAIN_menu_num) Then

	      status = SMG$Ring_Bell ( vdid, 3 )

	   Else

	      display_menu = .False.

	   End If

	End Do

	status = SMG$Erase_PasteBoard ( pbid )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Delete_Virtual_Keyboard ( kbid )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Delete_Virtual_Display ( vdid )
	If (status .Ne. SS$_Normal) Then
	   FEP_Display_Main = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If


	FEP_Display_Main = %Loc(FEP_Normal)

	Return
	End
