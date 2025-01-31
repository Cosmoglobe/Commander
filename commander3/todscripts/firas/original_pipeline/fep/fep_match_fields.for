	Integer*4 Function FEP_Match_Fields ( menu_page, menu_element,
     .						convert_flag, nfields )

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
C    AUTHOR: Rob Kummerer
C            STX
C            March 25, 1987
C
C    INVOCATION: STATUS = FEP_MATCH_FIELDS ( MENU_PAGE, MENU_ELEMENT,
C					     CONVERT_FLAG, NFIELDS )
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS:
C	MENU_PAGE(*)		I*2		Selected menu page number.
C	MENU_ELEMENT(*)		I*2		Selected menu field number.
C	CONVERT_FLAG(*)		L*1		[No]Convert flags.
C	NFIELDS			I*2		Number of selected fields.
C
C    SUBROUTINES CALLED:
C	STR$Match_Wild
C	Lib$Signal
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
C	FUT_Params
C	FEP_Invoc
C	FEP_Menu
C	$STRDEF
C
C----------------------------------------------------------------------
C    Changes:
C	SPR 3328. Choosing 12 fields (in /NOINTER) results in ERROR
C	     RETRIEVING POLYNOMIAL COEFFS.  Limit fields to a max of 9.
C	     Fred Shuman,  STX,  1989 Mar 03.
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
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FUT_Params)'
	Include 	'(FEP_Invoc)'
	Include 	'(FEP_Menu)'
	Include 	'($STRDEF)'

	Integer		*4	i
	Integer		*4	j
	Integer		*4	status
	Integer		*2	menu_page(*)
	Integer		*2	menu_element(*)
	Logical		*1	convert_flag(*)
	Integer		*2	ofields
	Integer		*2	nfields
	Integer		*2	maxfields
	Parameter	(maxfields = 4)
	Integer		*4	menu_len

	Integer		*4	STR$Match_Wild

	External	FEP_Normal
	External	FEP_MenuElErr
	External	FEP_MaxFldSel

	nfields = 0
	ofields = 0
c
c Try to match housekeeping header menu fields.
c
	Do i=1,HSKP_menu_num

           If (i .Eq. 1) Then
	      Type 5, ' Selected Housekeeping Header Menu fields:'
5	      Format (//, a, /)
	   End If

	   j = 0
	   status = STR$_NoMatch

	   Do While (j .Lt. fcc_fields_count .And. nfields .Lt. maxfields .And.
     .			status .Ne. STR$_Match)

	      j = j + 1

	      menu_len = index ( HSKP_menu(i), ' ' ) - 1
	      If (menu_len .Lt. 0) Then
		 Call LIB$Signal ( FEP_MenuElErr )
		 menu_len = 32
	      End If

	      status = STR$Match_Wild ( HSKP_menu(i)(1:menu_len),
     .				fcc_fields(j)(1:fcc_fields_len(j)) )

	      If (status .Eq. STR$_Match) Then
		 nfields = nfields + 1
		 menu_page(nfields) = HSKP_page
		 menu_element(nfields) = i
		 If (fcc_convert .Eq. fac_present) Then
		    convert_flag(nfields) = .True.
		    Type 10, HSKP_menu(i)(1:menu_len)
10		    Format ('      Convert   : ', a)
		 Else
		    convert_flag(nfields) = .False.
		    Type 20, HSKP_menu(i)(1:menu_len)
20		    Format ('      Noconvert : ', a)
		 End If
	      End If

	   End Do

	End Do

	If (nfields .Eq. ofields) Then
	   Type *, '     No Housekeeping Header Menu fields selected.'
	End If

c
c Try to match GRTA menu fields.
c
	ofields = nfields

	Do i=1,GRTA_menu_num

           If (i .Eq. 1) Then
	      Type 5, ' Selected A Side GRT Menu fields:'
	   End If

	   j = 0
	   status = STR$_NoMatch

	   Do While (j .Lt. fcc_fields_count .And. nfields .Lt. maxfields .And.
     .			status .Ne. STR$_Match)

	      j = j + 1

	      menu_len = index ( GRTA_menu(i), ' ' ) - 1
	      If (menu_len .Lt. 0) Then
		 Call LIB$Signal ( FEP_MenuElErr )
		 menu_len = 32
	      End If

	      status = STR$Match_Wild ( GRTA_menu(i)(1:menu_len),
     .				fcc_fields(j)(1:fcc_fields_len(j)) )

	      If (status .Eq. STR$_Match) Then
		 nfields = nfields + 1
		 menu_page(nfields) = GRTA_page
		 menu_element(nfields) = i
		 If (fcc_convert .Eq. fac_present) Then
		    convert_flag(nfields) = .True.
		    Type 10, GRTA_menu(i)(1:menu_len)
		 Else
		    convert_flag(nfields) = .False.
		    Type 20, GRTA_menu(i)(1:menu_len)
		 End If
	      End If

	   End Do

	End Do

	If (nfields .Eq. ofields) Then
	   Type *, '     No A Side GRT Menu fields selected.'
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

	   Do While (j .Lt. fcc_fields_count .And. nfields .Lt. maxfields .And.
     .			status .Ne. STR$_Match)

	      j = j + 1

	      menu_len = index ( GRTB_menu(i), ' ' ) - 1
	      If (menu_len .Lt. 0) Then
		 Call LIB$Signal ( FEP_MenuElErr )
		 menu_len = 32
	      End If

	      status = STR$Match_Wild ( GRTB_menu(i)(1:menu_len),
     .				fcc_fields(j)(1:fcc_fields_len(j)) )

	      If (status .Eq. STR$_Match) Then
		 nfields = nfields + 1
		 menu_page(nfields) = GRTB_page
		 menu_element(nfields) = i
		 If (fcc_convert .Eq. fac_present) Then
		    convert_flag(nfields) = .True.
		    Type 10, GRTB_menu(i)(1:menu_len)
		 Else
		    convert_flag(nfields) = .False.
		    Type 20, GRTB_menu(i)(1:menu_len)
		 End If
	      End If

	   End Do

	End Do

	If (nfields .Eq. ofields) Then
	   Type *, '     No B Side GRT Menu fields selected.'
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

	   Do While (j .Lt. fcc_fields_count .And. nfields .Lt. maxfields .And.
     .			status .Ne. STR$_Match)

	      j = j + 1

	      menu_len = index ( COMB_menu(i), ' ' ) - 1
	      If (menu_len .Lt. 0) Then
		 Call LIB$Signal ( FEP_MenuElErr )
		 menu_len = 32
	      End If

	      status = STR$Match_Wild ( COMB_menu(i)(1:menu_len),
     .				fcc_fields(j)(1:fcc_fields_len(j)) )

	      If (status .Eq. STR$_Match) Then
		 nfields = nfields + 1
		 menu_page(nfields) = COMB_page
		 menu_element(nfields) = i
		 If (fcc_convert .Eq. fac_present) Then
		    convert_flag(nfields) = .True.
		    Type 10, COMB_menu(i)(1:menu_len)
		 Else
		    convert_flag(nfields) = .False.
		    Type 20, COMB_menu(i)(1:menu_len)
		 End If
	      End If

	   End Do

	End Do

	If (nfields .Eq. ofields) Then
	   Type *, '     No Combined GRT Menu fields selected.'
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

	   Do While (j .Lt. fcc_fields_count .And. nfields .Lt. maxfields .And.
     .			status .Ne. STR$_Match)

	      j = j + 1

	      menu_len = index ( IPDU_menu(i), ' ' ) - 1
	      If (menu_len .Lt. 0) Then
		 Call LIB$Signal ( FEP_MenuElErr )
		 menu_len = 32
	      End If

	      status = STR$Match_Wild ( IPDU_menu(i)(1:menu_len),
     .				fcc_fields(j)(1:fcc_fields_len(j)) )

	      If (status .Eq. STR$_Match) Then
		 nfields = nfields + 1
		 menu_page(nfields) = IPDU_page
		 menu_element(nfields) = i
		 If (fcc_convert .Eq. fac_present) Then
		    convert_flag(nfields) = .True.
		    Type 10, IPDU_menu(i)(1:menu_len)
		 Else
		    convert_flag(nfields) = .False.
		    Type 20, IPDU_menu(i)(1:menu_len)
		 End If
	      End If

	   End Do

	End Do

	If (nfields .Eq. ofields) Then
	   Type *, '     No IPDU Temperature Menu fields selected.'
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

	   Do While (j .Lt. fcc_fields_count .And. nfields .Lt. maxfields .And.
     .			status .Ne. STR$_Match)

	      j = j + 1

	      menu_len = index ( VOLT_menu(i), ' ' ) - 1
	      If (menu_len .Lt. 0) Then
		 Call LIB$Signal ( FEP_MenuElErr )
		 menu_len = 32
	      End If

	      status = STR$Match_Wild ( VOLT_menu(i)(1:menu_len),
     .				fcc_fields(j)(1:fcc_fields_len(j)) )

	      If (status .Eq. STR$_Match) Then
		 nfields = nfields + 1
		 menu_page(nfields) = VOLT_page
		 menu_element(nfields) = i
		 If (fcc_convert .Eq. fac_present) Then
		    convert_flag(nfields) = .True.
		    Type 10, VOLT_menu(i)(1:menu_len)
		 Else
		    convert_flag(nfields) = .False.
		    Type 20, VOLT_menu(i)(1:menu_len)
		 End If
	      End If

	   End Do

	End Do

	If (nfields .Eq. ofields) Then
	   Type *, '     No Voltages Menu fields selected.'
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

	   Do While (j .Lt. fcc_fields_count .And. nfields .Lt. maxfields .And.
     .			status .Ne. STR$_Match)

	      j = j + 1

	      menu_len = index ( CURR_menu(i), ' ' ) - 1
	      If (menu_len .Lt. 0) Then
		 Call LIB$Signal ( FEP_MenuElErr )
		 menu_len = 32
	      End If

	      status = STR$Match_Wild ( CURR_menu(i)(1:menu_len),
     .				fcc_fields(j)(1:fcc_fields_len(j)) )

	      If (status .Eq. STR$_Match) Then
		 nfields = nfields + 1
		 menu_page(nfields) = CURR_page
		 menu_element(nfields) = i
		 If (fcc_convert .Eq. fac_present) Then
		    convert_flag(nfields) = .True.
		    Type 10, CURR_menu(i)(1:menu_len)
		 Else
		    convert_flag(nfields) = .False.
		    Type 20, CURR_menu(i)(1:menu_len)
		 End If
	      End If

	   End Do

	End Do

	If (nfields .Eq. ofields) Then
	   Type *, '     No Currents Menu fields selected.'
	End If

	If (nfields .Eq. maxfields) Then
	   Call Lib$Signal(FEP_MaxFldSel, %val(1), %val(maxfields))
	End If

	FEP_Match_Fields = %Loc(FEP_Normal)

	Return
	End
