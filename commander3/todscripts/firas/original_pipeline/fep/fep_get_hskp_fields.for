	Integer*4 Function FEP_Get_HSKP_Fields ( fr_lun )

C------------------------------------------------------------------------
C    PURPOSE: Retrieve field names from the housekeeping record and
C	      group according to menu.
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
C    INVOCATION: STATUS = FEP_GET_HSKP_FIELDS ()
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS: None
C
C    SUBROUTINES CALLED: 
C	FR_Get_Fields
C	LIB$Signal
C	STR$UpCase
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
C	FEP_Invoc
C	FEP_Menu
C	FUT_Params
C	FUT_Error
C	$SSDEF
C
C----------------------------------------------------------------------

	Implicit	None

	Include 	'(FEP_Invoc)'
	Include 	'(FEP_Menu)'
	Include		'(FUT_Params)'
	Include		'(FUT_Error)'
	Include		'($SSDEF)'

	Dictionary	'FEP_STOLDB'
	Record /FEP_STOLDB/ STOL_DB

	Integer		*4	i
	Integer		*4	j
	Integer		*4	k
	Integer		*4	status
	Logical		*4	fstatus
	Integer		*4	fr_lun
	Integer		*4	lun
	Character	*64	field_name
	Character	*64	xref_dataset_file
	Integer		*2	length
	Integer		*2	offset

	Character	*32	menu(32)

	Byte			db_hkp_xref(8*159)
	Byte			db_eqv(8)
	Character	*8	db_word

	Integer		*4	FUT_Field_Attributes
	Logical		*4	FR_Get_Fields
	Integer		*4	LIB$Get_LUN
	Integer		*4	LIB$Free_LUN
	Integer		*4	STR$UpCase
	Integer		*4	FR_Init

	External	FEP_Normal
	External	FEP_MenuElErr
	External	FEP_CFRGetFld
	External	FEP_RMSOpenXref
	External	FEP_RMSReadXref
	External	FEP_RMSCloseXref

	Equivalence	(db_eqv, db_word)

c
c Read the Housekeeping fields / STOL Database words cross-reference list.
c
        status = LIB$Get_LUN ( lun )

        If (status .Ne. SS$_Normal) Then
           FEP_Get_HSKP_Fields = status
	   Call LIB$Signal(%Val(status))
           Return
        End If

	xref_dataset_file = 'CSDR$FIRAS:HKP_STOLDB_XREF.DAT'

        Open ( Unit=lun, File=xref_dataset_file, Status='Old',
     .         Form='Unformatted', Shared, ReadOnly, IOstat=status )

        If (status .Ne. 0) Then
          FEP_Get_HSKP_Fields = %Loc(FEP_RMSOpenXref)
          Call LIB$Signal(FEP_RMSOpenXref,%Val(2),%Val(status),
     .				xref_dataset_file)
	  Return
        End If

	Read (lun, iostat=status) db_hkp_xref

        If (status .Ne. 0) Then
          FEP_Get_HSKP_Fields = %Loc(FEP_RMSReadXref)
          Call LIB$Signal(FEP_RMSReadXref,%Val(2),%Val(status),
     .				xref_dataset_file)
	  Return
        End If

	Close (lun, iostat=status)

        If (status .Ne. 0) Then
          FEP_Get_HSKP_Fields = %Loc(FEP_RMSCloseXref)
          Call LIB$Signal(FEP_RMSCloseXref,%Val(2),%Val(status),
     .				xref_dataset_file)
        End If

        status = LIB$Free_LUN ( lun )

        If (status .Ne. SS$_Normal) Then
           FEP_Get_HSKP_Fields = status
	   Call LIB$Signal(%Val(status))
        End If

c
c Get the housekeeping header menu.
c
	field_name = 'FEP_HKP.HSKP_Menu'

	fstatus = FR_Get_Fields ( field_name, HSKP_menu,
     .			          HSKP_menu_num )

	If (.Not. fstatus) Then
	   FEP_Get_HSKP_Fields = %Loc(FEP_CFRGetFld)
	   Call LIB$Signal(FEP_CFRGetFld)
	   Return
	End If

	Do i=1,HSKP_menu_num
	   field_name = 'FEP_STOLDB.HSKP_MENU.' // HSKP_menu(i)
	   status = FUT_Field_Attributes ( fr_lun, field_name, length, 
     .						offset )
	   Call LIB$Movc3 ( length, db_hkp_xref(offset+1), db_eqv )
	   HSKP_db(i) = db_word
	End Do

c
c Get the GRTA_MENU fields.
c
	field_name = 'FEP_HKP.GRTA_Menu'

	fstatus = FR_Get_Fields ( field_name, menu, GRTA_menu_num )

	If (.Not. fstatus) Then
	   FEP_Get_HSKP_Fields = %Loc(FEP_CFRGetFld)
	   Call LIB$Signal(FEP_CFRGetFld)
	   Return
	End If

	offset = GRTA_menu_num/2 + 1

	Do i=1,GRTA_menu_num,2
	   j = i/2
	   GRTA_menu(j+1) = menu(i)
	   GRTA_menu(j+offset) = menu(i+1)
	End Do

	Do i=1,GRTA_menu_num
	   field_name = 'FEP_STOLDB.GRTA_MENU.' // GRTA_menu(i)
	   status = FUT_Field_Attributes ( fr_lun, field_name, length, 
     .						offset )
	   Call LIB$Movc3 ( length, db_hkp_xref(offset+1), db_eqv )
	   GRTA_db(i) = db_word
	End Do

c
c Get the GRTB_MENU fields.
c
	field_name = 'FEP_HKP.GRTB_Menu'

	fstatus = FR_Get_Fields ( field_name, menu, GRTB_menu_num )

	If (.Not. fstatus) Then
	   FEP_Get_HSKP_Fields = %Loc(FEP_CFRGetFld)
	   Call LIB$Signal(FEP_CFRGetFld)
	   Return
	End If

	offset = GRTB_menu_num/2 + 1

	Do i=1,GRTB_menu_num,2
	   j = i/2
	   GRTB_menu(j+1) = menu(i)
	   GRTB_menu(j+offset) = menu(i+1)
	End Do

	Do i=1,GRTB_menu_num
	   field_name = 'FEP_STOLDB.GRTB_MENU.' // GRTB_menu(i)
	   status = FUT_Field_Attributes ( fr_lun, field_name, length, 
     .						offset )
	   Call LIB$Movc3 ( length, db_hkp_xref(offset+1), db_eqv )
	   GRTB_db(i) = db_word
	End Do

c
c Make the Combined GRT menu fields, using side A.
c
	COMB_menu_num = GRTA_menu_num

	Do i=1,COMB_menu_num/2
	   j = Index ( GRTA_menu(i), ' ' )
	   If (j .Lt. 1) Then
	      Call LIB$Signal ( FEP_MenuElErr )
	      j = 32
	   End If
	   COMB_menu(i) = GRTA_menu(i)(1:j-2)
	End Do

	k = 0
	Do i=COMB_menu_num/2+1,COMB_menu_num
	   k = k + 1
	   j = Index ( GRTB_menu(k), ' ' )
	   If (j .Lt. 1) Then
	      Call LIB$Signal ( FEP_MenuElErr )
	      j = 32
	   End If
	   COMB_menu(i) = GRTB_menu(k)(1:j-2)
	End Do

c
c Get the IPDU temperature menu.
c
	field_name = 'FEP_HKP.IPDU_Menu'

	fstatus = FR_Get_Fields ( field_name, IPDU_menu,
     .			          IPDU_menu_num )

	If (.Not. fstatus) Then
	   FEP_Get_HSKP_Fields = %Loc(FEP_CFRGetFld)
	   Call LIB$Signal(FEP_CFRGetFld)
	   Return
	End If

	Do i=1,IPDU_menu_num
	   field_name = 'FEP_STOLDB.IPDU_MENU.' // IPDU_menu(i)
	   status = FUT_Field_Attributes ( fr_lun, field_name, length, 
     .						offset )
	   Call LIB$Movc3 ( length, db_hkp_xref(offset+1), db_eqv )
	   IPDU_db(i) = db_word
	End Do

c
c Get the voltages menu.
c
	field_name = 'FEP_HKP.VOLT_Menu'

	fstatus = FR_Get_Fields ( field_name, VOLT_menu,
     .			          VOLT_menu_num )

	If (.Not. fstatus) Then
	   FEP_Get_HSKP_Fields = %Loc(FEP_CFRGetFld)
	   Call LIB$Signal(FEP_CFRGetFld)
	   Return
	End If

	Do i=1,VOLT_menu_num
	   field_name = 'FEP_STOLDB.VOLT_MENU.' // VOLT_menu(i)
	   status = FUT_Field_Attributes ( fr_lun, field_name, length, 
     .						offset )
	   Call LIB$Movc3 ( length, db_hkp_xref(offset+1), db_eqv )
	   VOLT_db(i) = db_word
	End Do

c
c Get the currents menu.
c
	field_name = 'FEP_HKP.CURR_Menu'

	fstatus = FR_Get_Fields ( field_name, CURR_menu,
     .			          CURR_menu_num )

	If (.Not. fstatus) Then
	   FEP_Get_HSKP_Fields = %Loc(FEP_CFRGetFld)
	   Call LIB$Signal(FEP_CFRGetFld)
	   Return
	End If

	Do i=1,CURR_menu_num
	   field_name = 'FEP_STOLDB.CURR_MENU.' // CURR_menu(i)
	   status = FUT_Field_Attributes ( fr_lun, field_name, length,
     .						offset )
	   Call LIB$Movc3 ( length, db_hkp_xref(offset+1), db_eqv )
	   CURR_db(i) = db_word
	End Do


	FEP_Get_HSKP_Fields = %Loc(FEP_Normal)

	Return
	End
