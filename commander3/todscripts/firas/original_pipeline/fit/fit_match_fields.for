	Integer*4 Function FIT_Match_Fields ( menu_page, menu_element,
     .						nfields )

C------------------------------------------------------------------------
C    PURPOSE: Match menu fields with those selected by the /FIELDS batch
C	      option.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Fred Shuman	(adapted from FEP_Match_Fields
C            STX		 by Rob Kummerer, STX, March 25, 1987)
C            October 5, 1987
C
C    INVOCATION: STATUS = FIT_MATCH_FIELDS ( MENU_PAGE, MENU_ELEMENT,
C					     NFIELDS )
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS: 
C	MENU_PAGE(*)		I*2		Selected menu page number.
C	MENU_ELEMENT(*)		I*2		Selected menu field number.
C	NFIELDS			I*2		Number of selected fields.
C
C    SUBROUTINES CALLED: 
C	STR$Match_Wild
C
C    COMMON VARIABLES USED:
C	GRTA_MENU(32)		CH*32		GRT side A menu.
C	GRTB_MENU(32)		CH*32		GRT side B menu.
C	COMB_MENU(32)		CH*32		Combined GRT menu.
C	IPDU_MENU(14)		CH*32		IPDU temps menu.
C	VOLT_MENU(20)		CH*32		Voltages menu.
C	CURR_MENU(28)		CH*32		Currents menu.
C
C    INCLUDE FILES: 
C	FUT_Error
C	FUT_Params
C	FIT_Invoc
C	FIT_Menu
C	$STRDEF
C  
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FUT_Params)'
	Include 	'(FIT_Invoc)'
	Include 	'(FIT_Menu)'
	Include 	'($STRDEF)'

	Integer		*4	i
	Integer		*4	j
	Integer		*4	status
	Integer		*2	menu_page(*)
	Integer		*2	menu_element(*)
	Integer		*2	ofields
	Integer		*2	nfields
	Integer		*4	menu_len

	Integer		*4	STR$Match_Wild

	External	FIT_Normal
	External	FIT_MenuElErr

	nfields = 0
	ofields = 0
c
c Try to match GRTA menu fields.
c

	Do i=1,GRTA_menu_num

	   If (i .Eq. 1) Then
	      Type 5, ' Selected A Side GRT Menu fields:'
5	      Format (//, a, /)
	   End If

	   j = 0
	   status = STR$_NoMatch

	   Do While (j .Lt. fcc_fields_count .And.
	2               status .Ne. STR$_Match)

	      j = j + 1

	      menu_len = index ( GRTA_menu(i), ' ' ) - 1
	      If (menu_len .Lt. 0) Then
		 Call LIB$Signal ( FIT_MenuElErr )
		 menu_len = 32
	      End If

	      status = STR$Match_Wild ( GRTA_menu(i)(1:menu_len),
	2                       fcc_fields(j)(1:fcc_fields_len(j)) )

	      If (status .Eq. STR$_Match) Then
	         nfields = nfields + 1
	         menu_page(nfields) = GRTA_page
	         menu_element(nfields) = i
	         Type 10, GRTA_menu(i)(1:menu_len)
10	         Format (1x,a)
	      End If

	   End Do

	End Do

	If (nfields .Eq. ofields) Then
	   Type *, '     No fields selected.'
	End If

c
c Try to match GRTB menu fields.
c
	ofields = nfields

	Do i=1,GRTB_menu_num

	   If (i .Eq. 1) Then
	      Type 5, ' Selected B Side GRT Menu fields:'
	   End If

	   j = 0
	   status = STR$_NoMatch

	   Do While (j .Lt. fcc_fields_count .And.
	2               status .Ne. STR$_Match)

	      j = j + 1

	      menu_len = index ( GRTB_menu(i), ' ' ) - 1
	      If (menu_len .Lt. 0) Then
		 Call LIB$Signal ( FIT_MenuElErr )
		 menu_len = 32
	      End If

	      status = STR$Match_Wild ( GRTB_menu(i)(1:menu_len),
	2                       fcc_fields(j)(1:fcc_fields_len(j)) )

	      If (status .Eq. STR$_Match) Then
		 nfields = nfields + 1
		 menu_page(nfields) = GRTB_page
		 menu_element(nfields) = i
	         Type 10, GRTB_menu(i)(1:menu_len)
	      End If

	   End Do

	End Do

	If (nfields .Eq. ofields) Then
	   Type *, '     No fields selected.'
	End If

c
c Try to match Combined GRT menu fields.
c
	ofields = nfields

	Do i=1,COMB_menu_num

	   If (i .Eq. 1) Then
	      Type 5, ' Selected Combined GRT Menu fields:'
	   End If

	   j = 0
	   status = STR$_NoMatch

	   Do While (j .Lt. fcc_fields_count .And.
	2               status .Ne. STR$_Match)

	      j = j + 1

	      menu_len = index ( COMB_menu(i), ' ' ) - 1
	      If (menu_len .Lt. 0) Then
		 Call LIB$Signal ( FIT_MenuElErr )
		 menu_len = 32
	      End If

	      status = STR$Match_Wild ( COMB_menu(i)(1:menu_len),
	2                       fcc_fields(j)(1:fcc_fields_len(j)) )

	      If (status .Eq. STR$_Match) Then
		 nfields = nfields + 1
		 menu_page(nfields) = COMB_page
		 menu_element(nfields) = i
	         Type 10, COMB_menu(i)(1:menu_len)
	      End If

	   End Do

	End Do

	If (nfields .Eq. ofields) Then
	   Type *, '     No fields selected.'
	End If

c
c Try to match IPDU menu fields.
c
	ofields = nfields

	Do i=1,IPDU_menu_num

	   If (i .Eq. 1) Then
	      Type 5, ' Selected IPDU Temperature Menu fields:'
	   End If

	   j = 0
	   status = STR$_NoMatch

	   Do While (j .Lt. fcc_fields_count .And.
	2               status .Ne. STR$_Match)

	      j = j + 1

	      menu_len = index ( IPDU_menu(i), ' ' ) - 1
	      If (menu_len .Lt. 0) Then
		 Call LIB$Signal ( FIT_MenuElErr )
		 menu_len = 32
	      End If

	      status = STR$Match_Wild ( IPDU_menu(i)(1:menu_len),
	2                       fcc_fields(j)(1:fcc_fields_len(j)) )

	      If (status .Eq. STR$_Match) Then
		 nfields = nfields + 1
		 menu_page(nfields) = IPDU_page
		 menu_element(nfields) = i
	         Type 10, IPDU_menu(i)(1:menu_len)
	      End If

	   End Do

	End Do

	If (nfields .Eq. ofields) Then
	   Type *, '     No fields selected.'
	End If

c
c Try to match VOLT menu fields.
c
	ofields = nfields

	Do i=1,VOLT_menu_num

	   If (i .Eq. 1) Then
	      Type 5, ' Selected Voltages Menu fields:'
	   End If

	   j = 0
	   status = STR$_NoMatch

	   Do While (j .Lt. fcc_fields_count .And.
	2               status .Ne. STR$_Match)

	      j = j + 1

	      menu_len = index ( VOLT_menu(i), ' ' ) - 1
	      If (menu_len .Lt. 0) Then
		 Call LIB$Signal ( FIT_MenuElErr )
		 menu_len = 32
	      End If

	      status = STR$Match_Wild ( VOLT_menu(i)(1:menu_len),
	2                       fcc_fields(j)(1:fcc_fields_len(j)) )

	      If (status .Eq. STR$_Match) Then
		 nfields = nfields + 1
		 menu_page(nfields) = VOLT_page
		 menu_element(nfields) = i
	         Type 10, VOLT_menu(i)(1:menu_len)
	      End If

	   End Do

	End Do

	If (nfields .Eq. ofields) Then
	   Type *, '     No fields selected.'
	End If

c
c Try to match CURR menu fields.
c
	ofields = nfields

	Do i=1,CURR_menu_num

	   If (i .Eq. 1) Then
	      Type 5, ' Selected Currents Menu fields:'
	   End If

	   j = 0
	   status = STR$_NoMatch

	   Do While (j .Lt. fcc_fields_count .And.
	2               status .Ne. STR$_Match)

	      j = j + 1

	      menu_len = index ( CURR_menu(i), ' ' ) - 1
	      If (menu_len .Lt. 0) Then
		 Call LIB$Signal ( FIT_MenuElErr )
		 menu_len = 32
	      End If

	      status = STR$Match_Wild ( CURR_menu(i)(1:menu_len),
	2                       fcc_fields(j)(1:fcc_fields_len(j)) )

	      If (status .Eq. STR$_Match) Then
		 nfields = nfields + 1
		 menu_page(nfields) = CURR_page
		 menu_element(nfields) = i
	         Type 10, CURR_menu(i)(1:menu_len)
	      End If

	   End Do

	End Do

	If (nfields .Eq. ofields) Then
	   Type *, '     No fields selected.'
	End If


	FIT_Match_Fields = %Loc(FIT_Normal)

	Return
	End
