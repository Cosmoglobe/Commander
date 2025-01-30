	Integer*4 Function FIT_Extract_Eng_Stats ( fr_lun, menu_page,
	2                                          menu_element, stats_gmts,
	3                                          orb, stats_buff )

C-------------------------------------------------------------------------------
C    PURPOSE: Fetches the requested engineering statistics fields.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Fred Shuman
C            STX
C            October 5, 1987
C
C    INVOCATION: STATUS = FIT_Extract_Eng_Stats ( FR_LUN, MENU_PAGE,
C				  MENU_ELEMENT, STATS_GMTS, ORB, STATS_BUFF )
C
C    INPUT PARAMETERS:
C	FR_LUN			 I*4		Field retriever lgcl unit.
C	MENU_PAGE		 I*2		Selected menu page.
C	MENU_ELEMENT		 I*2		Selected menu element.
C	STATS_GMTS(max, 3)	Ch*14		GMTs; for each orbit:
C						       min, mean, max values.
C	ORB			Structure	(Defined in FIT_DATA.TXT)
C	  [ NRECS		    I*4		    No. of orbit records.
C	  [ STATS_DATA(504,max,
C	         num_stats_types)   B		    Entire orbit stats data.
C	  [ NUMRECS_ORB(max)	    R*4		    No. records in each orbit.
C				End Structure
C		Note: "max" is fac_max_stats_recs, declared in FUT_PARAMS.TXT
C				= 1000 as of 1988 Jun 6
C    OUTPUT PARAMETERS:
C	STATS_BUFF(*)		 R*4		ENG STATS data to be merged.
C
C    SUBROUTINES CALLED:
C	FUT_Field_Attributes
C	LIB$Movc3
C	LIB$Signal
C	CT_GMT_to_Binary
C	CCT_Open_Config
C	CCT_Get_Config_TOD
C	FUT_Combine_HiLo
C	CCT_Close_Config
C
C    COMMON VARIABLES USED:
C	GRTA_MENU(32)		Ch*32		GRT side A menu.
C	GRTB_MENU(32)		Ch*32		GRT side B menu.
C	COMB_MENU(32)		Ch*32		Combined GRT menu.
C	IPDU_MENU(14)		Ch*32		IPDU temps menu.
C	VOLT_MENU(20)		Ch*32		Voltages menu.
C	CURR_MENU(28)		Ch*32		Currents menu.
C	FCC_INTERACTIVE		 I*4		Interactive/Batch flag.
C
C    INCLUDE FILES:
C	FUT_Error
C	FUT_Params
C	FIT_Invoc
C	FIT_Menu
C	FIT_Data
C	$SSDEF
C
C-----------------------------------------------------------------------------
C	Revised:
C
C       version 4.5.1 QFix 530, SPRs 4669, 4670.  F. Shuman,  STX,  1989 Oct 2.
C           Successive cropping of endtime; shuffling of data.
C           Changes to FIT.for, FIT_Extract_Eng_Stats.for,
C           FIT_Merge_Eng_Stats.for, and FIT_Data.txt.
C
C       version 4.5.1 QFix xxx, SPR 4762.  F. Shuman,  STX,  1989 Oct 24.
C           FIT was using wrong method to combine hi/lo GRTs.  Eliminate call
C           to FUT_Combine_Temps in favor of a new module FIT_Combine_Temps,
C           written to follow the method used in FUT_Temperature_List.
C           FIT, FIT_Extract_Eng_Stats, FIT_Combine_Temps, and FIT_MSG.
C
C       Spr 5655, Add the hard coded menu text file.
C           Version 5.8, H. Wang, STX, Mar. 8, 1990.
C
C    1991 Sep 25, F. Shuman, STX, SPR 9078.
C       Rename FIT_Combine_Temps routine to FUT_Combine_HiLo and move into FUT
C       to be used by both FIT and FEP.
C
C    1991 Nov 26, F. Shuman,  Hughes STX,  SPR 9291.
C       The two reference files GRTwt & GRTtrans that GRTswt split into were
C       erroneously being switched in this routine.
C
C    1991 Dec 3,  F. Shuman,  Hughes STX,  SPR 9327.
C       Make hi/lo-current-combination method in FUT_Combine_HiLo consistent
C       w. that in FUT_Temp_List; remove Ch*1 SIDE from FUT_Combine_HiLo call.
C
C    1991 Dec 11,  F. Shuman,  Hughes STX,  SPR 9334.
C       Remove obsolete logical name 'CSDR$FIRAS_URef'--replace with
C       'CSDR$FIRAS_Ref'.
C
C    1991 Dec 17,  F. Shuman,  Hughes STX,  SPR 9352.
C       Correct some hard-coded menu numbers causing trouble with Bolo-B plots.
C-------------------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FUT_Params)'
	Include		'(FIT_Invoc)'
	Include		'(FIT_Menu)'
	Include		'(FIT_Data)'
	Include		'(CCT_Status_Record)'
	Include		'(CCT_Get_Config)'
	Include		'($SSDEF)'

C  Calling arguments:
	Record /Orb_Data/	orb	! Structure defined in FIT_Data

	Integer		*4	fr_lun
	Integer		*2	menu_page
	Integer		*2	menu_element
	Character	*14	stats_gmts(fac_max_stats_recs, 3)
	Real		*4	stats_buff(fac_max_stats_recs, num_stats_types)

C  Internal variables:
	Integer		*4	i
	Integer		*4	j
	Integer		*4	k
	Integer		*4	ios
	Integer		*4	status

	Character	*64	field
	Character	*64	lo_field
	Integer		*2	length
	Integer		*2	lo_length
	Integer		*2	offset
	Integer		*2	lo_offset
	Integer		*4	adtstart(2)	! starting ADT for the timerange
	Integer		*4	adtmean(2)	! mean ADT for current orbit
	Integer		*4	adtend(2)	! ending ADT for the timerange

C  Number and names of reference data sets:
	Integer		*4	number/2/
	Character	*32	name(2) /'CSDR$FIRAS_Ref:FEX_GRTRawWt',
	2                                'CSDR$FIRAS_Ref:FEX_GRTTrans'/

	Integer		*4	size(2)/256,384/    ! sizes of data sets
	Character	*1	access_mode/' '/    ! data set access mode
	Integer		*4	ncache/1/
	Integer		*4	lun(2)		! logical unit numbers
	Integer		*4	cindex(2)	! initial cache pointers
	Logical		*1	new_segment(2)	! flags for new segments
	Integer		*4	ref_count

	Integer		*4	dl
	Real		*4	dl_eqv
	Real		*4	grtwgt	   ! GRT weight from FEX_GRTRawWt
c                                            used as flag value
	Real		*4	trans_temp ! hi/lo current temperature
					   !   transition midpt
	Real		*4	trans_hwid ! hi/lo current temperature
					   !   transition half-width
	Real		*4	wt	   !weighting factor for combining temps
	Integer		*4	off

	Integer		*4	FUT_Combine_HiLo
	Integer		*4	FUT_Field_Attributes
	Integer		*4	CT_GMT_to_Binary
	Integer		*4	CCT_Open_Config
	Integer		*4      CCT_Get_Config_TOD
	Integer		*4	CCT_Close_Config

	External		FUT_Normal
	External		FUT_BadTemps
	External		FIT_Normal
	External		FIT_InvDataLen
	External		FIT_OpnConfigErr
	External		FIT_GetConfigErr
	External		FIT_ClsConfigErr

	Equivalence	(dl, dl_eqv)

	Record /config_status/	stat(1)
	Dictionary 'FEX_GrtRawWt'
	Dictionary 'FEX_GrtTrans'
	Structure /refdata/
	   Record /FEX_GrtRawWt/ GRTWt	    !Ref file containing GRT weights
	   Record /FEX_GrtTrans/ GRTTrans   !Ref file of GRT transition
	                                    !  temperature midpts & 1/2-widths
	End Structure
	Record /refdata/ ref

	FIT_Extract_Eng_Stats = %Loc(FIT_Normal)
	off = COMB_menu_num/2
C
C Retrieve the field name.
C
	If (menu_page .Eq. GRTA_page) Then
	   field = 'FIT_ETR.GRTA_MENU.' // GRTA_menu(menu_element)
	Else If (menu_page .Eq. GRTB_page) Then
	   field = 'FIT_ETR.GRTB_MENU.' // GRTB_menu(menu_element)
	Else If (menu_page .Eq. COMB_page) Then
	   If (menu_element .Le. COMB_menu_num/2) Then
	      If (menu_element .Gt. 4) Then
	         field = 'FIT_ETR.GRTA_MENU.' // GRTA_menu(menu_element)
	      End If
	      lo_field = 'FIT_ETR.GRTA_MENU.' // GRTA_menu(menu_element+off)
	   End If
	   If (menu_element .Gt. COMB_menu_num/2) Then
	      If (menu_element .Gt. 20) Then
	         field = 'FIT_ETR.GRTB_MENU.' // GRTB_menu(menu_element-off)
	      Endif
	      lo_field = 'FIT_ETR.GRTB_MENU.' // GRTB_menu(menu_element)
	   End If
	Else If (menu_page .Eq. IPDU_page) Then
	   field = 'FIT_ETR.IPDU_MENU.' // IPDU_menu(menu_element)
	Else If (menu_page .Eq. VOLT_page) Then
	   field = 'FIT_ETR.VOLT_MENU.' // VOLT_menu(menu_element)
	Else If (menu_page .Eq. CURR_page) Then
	   field = 'FIT_ETR.CURR_MENU.' // CURR_menu(menu_element)
	End If
C
C Get the attributes of the field being requested.
C
	If (menu_page .Eq. COMB_page) Then
	   If (menu_element .Le. COMB_menu_num/2) Then
	      If (menu_element .Gt. 4) Then
	         status = FUT_Field_Attributes ( fr_lun, field, length, offset )
	      End If
	      status=FUT_Field_Attributes( fr_lun,lo_field,lo_length,lo_offset )
	   End If
	   If (menu_element .Gt. COMB_menu_num/2) Then
	      If (menu_element .Gt. 20) Then
	         status = FUT_Field_Attributes ( fr_lun, field, length, offset )
	      Endif
	      status=FUT_Field_Attributes( fr_lun,lo_field,lo_length,lo_offset )
	   End If
	   If (status .Ne. %Loc(FUT_Normal)) Then
	      FIT_Extract_Eng_Stats = status
	      Return
	   End If
	Else
	   status = FUT_Field_Attributes ( fr_lun, field, length, offset )
	   If (status .Ne. %Loc(FUT_Normal)) Then
	      FIT_Extract_Eng_Stats = status
	      Return
	   End If
	End If
C
C Fetch the field contents.
C
	If (menu_page .Eq. COMB_page) Then
	   If (menu_element .Le. COMB_menu_num/2 .And. menu_element .Le. 4) Then
	      length = lo_length
	      offset = lo_offset
	   End If
	   If (menu_element .Gt. COMB_menu_num/2 .And. menu_element .Le. 20) Then
	      length = lo_length
	      offset = lo_offset
	   End If
	End If
	If (length .Eq. 4) Then
	   Do j=1,orb.nrecs
	      Do k=1,num_stats_types
	         Call LIB$Movc3 ( length,orb.stats_data(offset+1,j,k), dl )
	         stats_buff(j,k) = dl_eqv
	      End Do
	   End Do
	Else
	   Call LIB$Signal ( FIT_InvDataLen )
	End If
C
C If "Combined-GRTs", fetch the contents of the second field...
C
	If (menu_page .Eq. COMB_page .And. (menu_element .Gt. COMB_menu_num/2+4 .Or.
	2    menu_element .Gt. 4 .And. menu_element .Le. COMB_menu_num/2)) Then
	   If (lo_length .Eq. 4) Then
C
C   Open FEX_GRTTrans and FEX_GRTRawWt reference files, which will be used
C      to Get GRT weights from the module FUT_Combine_HiLo.
C
	      status = CT_GMT_to_Binary ( stats_gmts(1,2), adtstart )
	      status = CT_GMT_to_Binary ( stats_gmts(orb.nrecs,2), adtend )
	      ios = CCT_Open_Config ( adtstart, adtend,
	2                             number, name, size, access_mode,
	3                             ncache, lun, cindex, stat, ref_count )
	      If (.Not. ios) Then
	         Call LIB$Signal(FIT_OpnConfigErr, %Val(1), %Val(ios))
	         FIT_Extract_Eng_Stats = %Loc(FIT_OpnConfigErr)
	      End If

	      Do j=1,orb.nrecs
C
C   Get GRT reference data using GET_CONFIG.
C
	         status = CT_GMT_to_Binary ( stats_gmts(j, 2), adtmean )
	         ios = CCT_Get_Config_TOD ( adtmean,
	2                                   number, size, lun, cindex,
	3                                   ref, new_segment, stat )

	         If (ios) Then
C
C  Extract GRT weights and transition temps from the input FEX_GRTRawWt and
C    FEX_GRTTrans records.
C
	            If (offset .Lt. 128) Then
	               i = offset/4 - 15
	               grtwgt = ref.GRTWT.GRT_A_Weight(i)
	               trans_temp = ref.GRTTRANS.GRT_A_Trans_Temp(i)
	               trans_hwid = ref.GRTTRANS.GRT_A_Trans_HWid(i)
	            Else
	               i = offset/4 - 47
	               grtwgt = ref.GRTWT.GRT_B_Weight(i)
	               trans_temp = ref.GRTTRANS.GRT_B_Trans_Temp(i)
	               trans_hwid = ref.GRTTRANS.GRT_B_Trans_HWid(i)
	            End If
	            Do k=1,num_stats_types-1
	               Call LIB$Movc3 ( lo_length,
	2                               orb.stats_data(lo_offset+1,j,k), dl )
C
C ...and combine the lo & hi fields using the GRT weights (for std. dev. [k=2]
C     use the same weight applied for the mean [k=1] )...
C
	               If (k .Eq. 2 .And. grtwgt .Ne. 0.) Then
	                  stats_buff(j,k) = (1.0 - wt)*dl_eqv +
	2                                           wt*stats_buff(j,k)
	               Else If (grtwgt .Ne. 0.) Then
	                  status = FUT_Combine_HiLo ( stats_buff(j,k), dl,
	2                                 trans_temp, trans_hwid, wt)
	                  If (wt.Eq.-9999. .Or. status .Eq. %Loc(FUT_BadTemps))
	2                 Then
	                     stats_buff(j,k) = -9999.
	                  Else
	                     stats_buff(j,k) = (1.0 - wt)*dl_eqv +
	2                                              wt*stats_buff(j,k)
	                  End If
	               End If
	            End Do
	            Call LIB$Movc3 ( lo_length,
	2                       orb.stats_data(lo_offset+1,j,num_stats_types),
	3                       stats_buff(j,num_stats_types) )
	         Else
	            Call LIB$Signal(FIT_GetConfigErr, %Val(1), %Val(ios))
	            FIT_Extract_Eng_Stats = %Loc(FIT_GetConfigErr)
	         End If
	      End Do
C
C   Close configuration file.
C
	      ios = CCT_Close_Config ( number, lun, cindex )
	      If (.Not. ios) Then
	         Call LIB$Signal(FIT_ClsConfigErr, %Val(1), %Val(ios))
	         FIT_Extract_Eng_Stats = %Loc(FIT_ClsConfigErr)
	      End If

	   Else
	      Call LIB$Signal ( FIT_InvDataLen )
	   End If
	End If

	Return
	End
