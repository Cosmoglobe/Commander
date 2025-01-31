	Integer*4 Function FEP_Convert_To_Eng ( fr_lun,Config_Lun,con_lun,
	2                                       menu_page, menu_element,
	3                                       convert_flag, digital_mask,
	4                                       conv_buff, dwell_mode,
	5                                       co_coeffs )

C------------------------------------------------------------------------------
C    PURPOSE: Converts the housekeeping data to engineering data.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX
C            April 7, 1987
C
C    INVOCATION: STATUS = FEP_CONVERT_TO_ENG ( MENU_PAGE, Config_Lun,
C                                              Con_lun,MENU_ELEMENT,
C					       CONVERT_FLAG, DIGITAL_MASK,
C					       CONV_BUFF, DWELL_MODE )
C
C    INPUT PARAMETERS:
C       Config_Lun(3)           I*4             Fortran unit number of
C                                               reference files
C       Con_Lun(2   )           I*4             Fortran unit number of
C                                               reference files
C	MENU_PAGE		I*2		Selected menu page.
C	MENU_ELEMENT		I*2		Selected menu element.
C	CONVERT_FLAG		L*1		[No]Convert selected field.
C	DIGITAL_MASK		I*2		Mask for digital status flds.
C	CONV_BUFF(*)		R*4		HSKP data to be converted.
C	CO_COEFFS(3,4)		R*4		Counts to ohms conversion
C						coefficients.
C
C    OUTPUT PARAMETERS:
C	CONV_BUFF(*)		R*4		Converted HSKP data.
C	DWELL_MODE(*)		I*1		Dwell mode data?
C
C    SUBROUTINES CALLED:
C	FUT_Field_Attributes
C	FEP_Get_GRT_Conv
C	LIB$Polyf
C	FEP_Counts_To_Ohms
C	FEP_GRT_Lookup
C	LIB$Movc3
C	GETCURVES
C	LIB$Signal
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
C	FCC_INTERACTIVE		I*4		Interactive/Batch flag.
C	FCC_CONVERT		I*4		[No]Convert flag.
C
C    INCLUDE FILES:
C	FEP_Invoc
C	FEP_Menu
C	FEP_Data
C	FUT_Error
C	FUT_Params
C	$SSDEF
C
C    Software Maintainers, CAUTION:  Hardcoded numbers used in tests of
C    menu_element().  Look for the comment block:
C	c
C	c  BEWARE--HARDCODED MENU ELEMENT NUMBERS!  REVISE WHEN MENU CHANGES!
C	c
C
C----------------------------------------------------------------------
C PDL:
C
C	If ((menu_page equal  GRT A side page, or GRT B side page,
C        or GRT A side B side combine page
C       Then
C
C         Call Fut_Field_attributes to
C           get the dwell status mode offsets, side A and B.
C         If return status is bad
C         then
C            set return status to bad
C            return
C         Endif
C
C         Convert any GRT menu item except cal resistors.
C
C	   If (Option Conversion is true)
C          Then
C            Convert GRT temperatures.
C
C           Retrieve the field names, database words, and field
C           GRT conversion information.
C
C	      If (menu_page .Eq. GRT Side A) Then
C
C               Get the data_base word
C               Call Fep_Get_Grt_Attr to
C                  get the attributes of the field being converted.
C              If return status is bad
C              then
C                set return status to bad
C                return
C              Endif
C
C	    Else If (menu_page .Eq. GRT B side)
C           Then
C
C               Get the data_base word
C               Call Fep_Get_Grt_Attr to
C                  get the attributes of the field being converted.
C              If return status is bad
C              then
C                set return status to bad
C                return
C              Endif
C
C	      Else If (menu_page .Eq. GRT A and B side combined)
C             Then
C
C                Combine low/high currents of either A or B.
C
C             Endif
C
C             Do for each record
C
C               Check the dwell status.  Don't convert if in dwell mode.
C
C               Do the GRTA menu field, GRTB menu field or
C               the low current portion of GRT combination.
C
C              Get the calibrator resistor settings.
C
C              Get the GRT counts.
C
C              Call FEP_Counts_To_Ohms to
C                 Convert from counts to ohms.
C              Call  cct_get_config_idx_tod to get the
C                index and new_ref_data_flag
C              if new_ref_data_flag is true
C              Then
C                Call Fep_Get_Grt_Conv to get
C                 Counts to ohms conversion table, ohms to temps
C                 conversion table, number of points in conversion tables
C              Endif
C              Call Fep_Grt_lookup to
C                  Convert from ohms to temperature.
C
C           END DO FOR each record( for GRT CONVERSION)
C	   Else
C             do no grt conversion
C          endif
C       Else  ( Not for GRT menu,
C             Convert via polynomial conversion or do no conversion at all)
C
C           IF Do no conversion at all:
C           Then
C             Retrieve the field name.
C
C             Get the attributes of the field being converted.
C             Fetch the field contents.
C             For digital status fields, mask and return.
C           ENDIF
C           If (convert via polynomial conversion)
C           Then
C	      Do for each record
C              Call  cct_get_config_idx_tod to get the
C                index and new_ref_data_flag
C              if new_ref_data_flag is true
C              Then
C                Call Fep_Get_Curve to get poly. Coeffs. from the archive
C              Endif
C              Call LIB$Polyf to do polynomial conversion of data
C	      End Do
C           Endif
C         Endif
C	Return
C	End
C----------------------------------------------------------------------
C
C Changes:
C
C	SPR 1507, Signal TBLEDGE message once per plot.  R. Kummerer,
C	December 11, 1987.
C
C	SPR 1723, Bolometer GRT use low current only; Combine Menu
C	must avoid using high current for these.  R. Kummerer, December
C	15, 1987.
C
C	SPR 1985, Transition range for merging hi and lo current GRT
C	readings changed from 500-700 to 15450-15660 due to change in
C	status monitor electronics.  R. Isaacman, 7 Oct 88
C
C	SPR 3074, Low/high current GRTs temperature conversion correction
C	caused by GET_FIELDS now returning alphabetical order fields.
C	ENGPLOTS interpreted low current temperatures as high current.
C	R. Kummerer, January 4, 1988.
C
C	SPR 2377, Provide masking capability for digital status words.
C	Fred Shuman,  1989 Feb 20.
C
C	22-Jun-1989	Added the housekeeping record index, i, to the calling
C			sequence of FEP_Counts_to_ohms.  Don Stevens-Rayburn.
C
C	14-Jul-1898	Added the "pre-computed" counts to ohms coefficients
C			to the calling sequence.  This routine simply passes
C			those coefficients through to FEP_Counst_to_Ohms.
C
C       SPR 4852, corrections for handling the dwell data right.
C       Qfix #, Harte Wang, STX, 1989 Nov 3
C
C       SPR 4971, corrections for handling the dwell data right for the
C        cal resistor plots of menu pages of grta, grtb, comb
C        Qfix #, Harte Wang, STX, 1989 Nov 10
C
C       SPR 5109, When combining GRTs, convert the low current GRT regardless
C	of the success of the high current GRT and then combine temperatures
C	according to the temperature regime AND the success of either
C	conversions.  Version 4.4.05, R. Kummerer, STX, Nov 20, 1989.
C
C       SPR 5068, avoid plotting data points when telm. qual. is bad
C       Qfix #, Harte Wang, STX, 1989 dec. 2
C
C       SPR 5130, Spurious temperature jumps in plots
C       Qfix #, Harte Wang, STX, 1989 dec. 8
C
C       SPR 5078, failure to plot some unavailable items
C       Qfix #, Harte Wang, STX, 1989 dec. 20
C
C       SPR 5540, Confusion over DWELL MODE.
C                 Harte Wang, STX, 1990 Jan. 11
C
C       SPR 3921, Add the capabalities to read reference data from the
C                 reference archive.
C                 Harte Wang, STX, 1990 Feb. 10.
C
C       SPR 6253, Xcal_temp_s6_b plot is scaled to 6 x E04 for Xcal motion.
C                 Harte Wang, STX, 1990 Feb. 22.
C
C	SPR 8956, remove grt_correction routine.
C		  L. Rosen, STX, 6 sept. 1991
C
C	SPR 9079, Invoke FUT_COMBINE_HI LO to combine Hi, Lo current readings
C		  for GRT temperatures, STX, H.WANG, 1991 Sep. 27.
C
C	SPR 9290, FEP using new GRT ref. files in wrong order
C		  Hughes STX, H.WANG, 1991 Nov. 26.
C
C	SPR 9327, Make the hi/lo current combination method in FUT_Combine_HiLo
C		  consistent with that in FUT_Temp_List; remove tran_side from
C		  FUT_Combine_HiLo call.  Hughes STX,  F. Shuman,  1991 Dec 3.
C------------------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FUT_Params)'
	Include		'(FEP_Invoc)'
	Include		'(FEP_Menu)'
	Include		'(FEP_firnames)'
	Include		'(FEP_Data)'
	Include		'($SSDEF)'
	Include		'(CCT_GET_CONFIG)'
	Record/config_status/stat1(3)
	Record/config_status/stat2(2)

	Integer		*4	i
	Integer		*4	j
	Integer		*4	k
	Integer		*4	status
	Integer		*4	Config_lun(3)
	Integer		*4	Con_lun(2)
	Integer		*4	Size(2)/256,384/
	Integer		*4	index(3)
	Integer		*4	indx(2)
	Logical		*1	New_GRTWT
	Logical		*1	New_GRT
	Logical		*1	New_Poly
	Integer		*4	fr_lun
	Integer		*2	menu_page
	Integer		*2	menu_element_b
	Integer		*2	menu_element
	Logical		*1	convert_flag
	Logical		*1	first
	Logical		*1	lo_first
	Logical		*1	poly_first
	Integer		*2	digital_mask
	Real		*4	conv_buff(*)
	Byte			dwell_mode(*)
	Byte			dwell_byte
	Byte			Null(2)/0,0/
	Integer		*2	Equiv,dwecont
	Equivalence (null(1), equiv)
	Character	*64	field
	Character	*64	lo_field
	Character	*8	dbword
	Character	*8	lo_dbword
	Integer		*2	length
	Integer		*2	lo_length
	Integer		*2	offset
	Integer		*2	lo_offset
	Character	*64	dwell_field
	Character	*64	telm_field
	Integer		*2	dwell_length
	Integer		*2	telm_length
	Integer		*2	dwell_offset_a
	Integer		*2	dwell_offset_b
	Integer		*2	telm_offset
	Byte			dwell_stat_a
	Byte			telm_stat
	Byte			dwell_stat_b
	Byte			dl1
	Integer		*2	dl2
	Integer		*4	off

	Integer		*4	calres(4)
	Integer		*4	lo_calres(4)
	Integer		*4	caloff
	Integer		*4	lo_caloff
	Integer		*4	temp_cnts
	Integer		*4	lo_temp_cnts
	Real		*4	resistance
	Integer		*4	iside
	Real		*4	ohms(400)
	Real		*4	ctemps(400)
	Integer		*4	npts
	Integer		*4	lo_iside
	Real		*4	lo_ohms(400)
	Real		*4	lo_ctemps(400)
	Integer		*4	lo_npts
	Real		*4	grt_temp
	Real		*4	lo_grt_temp
	Real		*4	dt
	Integer		*4	icurrent
	Real		*4	grt_wt
	Real		*4	co_coeffs ( 3, 4 )

	Integer		*2	indices(2)
	Integer		*2	npoly
	Integer		*2	istat
	Real		*4	coeffs(6)
	Real		*4	poly_coeffs(6)
	Integer		*2	ndeg

	Integer		*2	itemp
	Logical		*1	signaled
	Logical		*1	not_avail
	Integer		*4	FUT_Field_Attributes
	Integer		*4	FEP_Get_GRT_Conv
	Integer		*4	FEP_Get_GRT_ATTR
	Integer		*4	FEP_Counts_To_Ohms
	Integer		*4	FEP_GRT_Lookup
	Integer		*4	CCT_GET_CONFIG_IDX_TOD
	Integer		*4	CCT_GET_CONFIG_TOD
	Integer		*4	FUT_COMBINE_HILO
	Integer		*4	tranidx(16) / 6, 7, 8, 9,11,12,13,14,
	2                                    16, 5, 4,10, 3, 2, 1,15/
	Real		*4	trans_temp,trans_hwid

C
	Dictionary    'fex_grtrawwt'
	Dictionary    'fex_grttrans'
C
	Structure /FEP_config/
	   Record /fex_grtrawwt/ grtrawwt
	   Record /fex_grttrans/ grttrans
	End Structure
	Record /fep_config/ config
C
	External		FEP_Normal
	External		FEP_InvDataLen
	External		FEP_GetPolyErr
	External		FEP_TblEdge
	External		FEP_GetConfigerr
	External		FUT_Normal

c
c Perform conversion to engineering units. Note that housekeeping
c fields do not require conversion.
c
c  BEWARE--HARDCODED MENU ELEMENT NUMBERS!  REVISE WHEN MENU CHANGES!
c
	Do j=1, nrecs
	  Dwell_mode(j) = 0
	End Do
	Dwell_offset_a = 0
	Dwell_offset_b = 0
	Dwell_stat_a = 0
	Dwell_stat_b = 0
	telm_offset=0
	telm_stat = 0
	telm_field = 'FEP_HKP.CURR_MENU.TELEMETRY_QUALITY'

	status = FUT_Field_Attributes ( fr_lun, telm_field,
	2                               telm_length, telm_offset )
	If (status .Ne. %Loc(FUT_Normal)) Then
	  FEP_Convert_To_Eng = status
	  Return
	End If

	If ((menu_page .Eq. GRTA_page .Or.
	2    menu_page .Eq. GRTB_page .Or.
	3    menu_page .Eq. COMB_page) .And.
	4   (menu_element .Lt. 5 .Or.
	5    (menu_element .Gt. 8 .And. menu_element .Lt. 21) .Or.
	6    menu_element .Gt. 24)) Then
c
c Get the dwell status mode offsets, side A and B.
c
	  dwell_field = 'FEP_HKP.HSKP_MENU.DWELL_STATUS_A'

	  status = FUT_Field_Attributes ( fr_lun, dwell_field,
	2                                 dwell_length, dwell_offset_a )
	  If (status .Ne. %Loc(FUT_Normal)) Then
	    FEP_Convert_To_Eng = status
	    Return
	  End If

	  dwell_field = 'FEP_HKP.HSKP_MENU.DWELL_STATUS_B'

	  status = FUT_Field_Attributes ( fr_lun, dwell_field,
	2                                 dwell_length, dwell_offset_b )
	  If (status .Ne. %Loc(FUT_Normal)) Then
	    FEP_Convert_To_Eng = status
	    Return
	  End If
c
c Convert any GRT menu item except cal resistors.
c
	  If (convert_flag) Then
c
c Convert GRT temperatures.
c
c Retrieve field names, database words, and field GRT conversion information.
c
	    If (menu_page .Eq. GRTA_page) Then
c
c Side A.
c
	      field = 'FEP_HKP.GRTA_MENU.' // GRTA_menu(menu_element)
	      dbword = GRTA_db(menu_element)
	      status = FEP_Get_GRT_Attr ( fr_lun, field,
	2                                   length, offset, caloff, iside)
	      If (status .Ne. %Loc(FUT_Normal)) Then
	        FEP_Convert_To_Eng = status
	        Return
	      End If
	      caloff = caloff - 2

	    Else If (menu_page .Eq. GRTB_page) Then
c
c Side B.
c
	      field = 'FEP_HKP.GRTB_MENU.' // GRTB_menu(menu_element)
	      dbword = GRTB_db(menu_element)
	      status = FEP_Get_GRT_Attr ( fr_lun, field,
	2                                 length, offset, caloff, iside)
	      If (status .Ne. %Loc(FUT_Normal)) Then
	        FEP_Convert_To_Eng = status
	        Return
	      End If
	      caloff = caloff - 2

	    Else If (menu_page .Eq. COMB_page) Then
c
c Combine low/high currents of either A or B.
c
	      off = COMB_menu_num/2
	      If (menu_element .Le. off ) Then
c
c Side A.
c
	        If (menu_element .Gt. 4) Then
	          field = 'FEP_HKP.GRTA_MENU.' // GRTA_menu(menu_element)
	          dbword = GRTA_db(menu_element)
	          status = FEP_Get_GRT_Attr ( fr_lun, field,
	2                                     length, offset, caloff, iside)
	          If (status .Ne. %Loc(FUT_Normal)) Then
	            FEP_Convert_To_Eng = status
	            Return
	          End If
	          caloff = caloff - 2
	        End If

	        lo_field = 'FEP_HKP.GRTA_MENU.' // GRTA_menu(menu_element+off)
	        lo_dbword = GRTA_db(menu_element+off)
	        status = FEP_Get_GRT_Attr ( fr_lun, lo_field, lo_length,
	2                                   lo_offset, lo_caloff, lo_iside)
	        If (status .Ne. %Loc(FUT_Normal)) Then
	          FEP_Convert_To_Eng = status
	          Return
	        End If
	        lo_caloff = lo_caloff - 2
	      Else
c
c Side B.
c
	        If ((menu_element - off) .Gt. 4) Then
	          field = 'FEP_HKP.GRTB_MENU.' // GRTB_menu(menu_element-off)
	          dbword = GRTB_db(menu_element-off)
	          status = FEP_Get_GRT_Attr (fr_lun, field,
	2                                    length, offset, caloff, iside)
	          If (status .Ne. %Loc(FUT_Normal)) Then
	            FEP_Convert_To_Eng = status
	            Return
	          End If
	          caloff = caloff - 2
	        End If
	        lo_field = 'FEP_HKP.GRTB_MENU.' // GRTB_menu(menu_element)
	        lo_dbword = GRTB_db(menu_element)
	        status = FEP_Get_GRT_Attr ( fr_lun, lo_field,
	2                            lo_length, lo_offset, lo_caloff, lo_iside)
	        If (status .Ne. %Loc(FUT_Normal)) Then
	          FEP_Convert_To_Eng = status
	          Return
	        End If
	        lo_caloff = lo_caloff - 2
	      End If        !(menu_element .Le. off)

	    End If          !(menu_page .Eq. GRTA_page)
c
c Do the conversion.
c
	    signaled = .False.              ! Signal error messages just once.
	    first = .True.
	    Lo_first = .True.
	    not_avail = .False.

	    Do i=1,nrecs

	      dwell_mode(i) = 0
c
c Check the dwell status.  Don't convert if in dwell mode.
c
	      dwell_stat_a = 0
	      dwell_stat_b = 0
	      Telm_stat = 0
	      Status = CCT_Get_Config_Idx_Tod(timetags(i).time,1,
	2                         config_lun(2),index(2),New_Grt,Stat1(2))
	      If (.Not. Status) Then
	        Call lib$signal(Fep_GetConfigErr,%val(1), %val(status))
	        Fep_Convert_to_eng = Status
	        Return
	      End If
	      status = CCT_Get_Config_Tod(timetags(i).time,2,size,
	2                           con_lun,indx,config,New_GrtWT,Stat2)
	      If (.Not. Status) Then
	        Call lib$signal(Fep_GetConfigErr,%val(1), %val(status))
	        Fep_Convert_to_eng = Status
	        Return
	      End If
	      If (hskp_data(telm_offset+1,i) .Ne. 0) Then
	        telm_stat = 1
	      End If

	      If (hskp_data(dwell_offset_a+1,i) .Lt. 0) Then
	        dwell_byte=hskp_data(dwell_offset_a+1,i)
	        null(1) = dwell_byte
	        Call mvbits(equiv,0,5,dwecont,0)
	        if (dwecont .Lt. 31) dwell_stat_a = 1
	      End If

	      If (hskp_data(dwell_offset_b+1,i) .Lt. 0) Then
	        dwell_byte=hskp_data(dwell_offset_b+1,i)
	        null(1) = dwell_byte
	        Call mvbits(equiv,0,5,dwecont,0)
	        if (dwecont .Lt. 31) dwell_stat_b = 1
	      End If

	      If (telm_stat .Eq. 0) Then
	        If ((dwell_stat_a .Eq. 0 .And.
	2            (menu_page .Eq. GRTA_page .Or.
	3             (menu_page .Eq. COMB_page .And. menu_element .Le. off)))
	4      .Or. (dwell_stat_b .Eq. 0 .And.
	5            (menu_page .Eq. GRTB_page .Or.
	6             (menu_page .Eq. COMB_page .And. menu_element .Gt. off))))
	7       Then
c
c Do the GRTA field, GRTB field or the low current portion of GRT combination.
c
c Get the calibrator resistor settings.
c
	          If (menu_page .Eq. comb_page .And.
	2             (menu_element .Le. off .And. menu_element .Lt. 5   .Or.
	3              menu_element .Gt. off .And. menu_element .Lt. 21)) Then
	            not_avail = .True.
	          End If
	          If (.Not. not_avail) Then
	            Do j=1,4
	              Call LIB$Movc3 ( 2, hskp_data(caloff+j*2+1,i), dl2 )
	              calres(j) = dl2
	              If (dl2 .Lt. 0) Then
	                calres(j) = 65536 + dl2
	              End If
	            End Do
c
c Get the GRT counts.
c
	            Call LIB$Movc3 ( 2, hskp_data(offset+1,i), dl2 )
	            temp_cnts = dl2
	            If (dl2 .Lt. 0) Then
	              temp_cnts = 65536 + dl2
	            End If
	            If (New_Grt .Or. first) Then
	              Status= Fep_Get_Grt_Conv(Config_Lun(2),dbword,
	2                                      Ohms,Ctemps,Npts)
	              If (Status .Ne. %Loc(FEP_Normal)) Then
	                Fep_Convert_To_Eng = status
	                Return
	              End If
	            End If

	            If (temp_cnts .Ne. 0) Then
c
c Convert from counts to ohms.
c
	              status = FEP_Counts_To_Ohms ( config_lun(1), first,
	2                                           calres, temp_cnts, 1, iside,
	3                                           resistance, i, co_coeffs )
c
c Convert from ohms to temperature.
c
	              status = FEP_GRT_LookUp ( resistance, 0, ohms, ctemps,
	2                                       first,npts, grt_temp, dt )
	              first = .False.

	              If (status .Eq. %Loc(FEP_TBLEDGE)) Then
	                If (.Not. Signaled) Then
	                  Call Lib$Signal(FEP_TblEdge)
	                  Signaled = .True.
	                End If
	                grt_temp = -99999.0
	              Else
	                If (iside .Le. 2) Then
	                  icurrent = 1
	                Else
	                  icurrent = 2
	                End If
	              End If
	            Else
	              grt_temp = -99999.0
	            End If     !(temp_cnts .Ne. 0)

	            conv_buff(i) = grt_temp
	            dwell_mode(i) = 0

	          End If         !(.Not. not_avail)
c
c  BEWARE--HARDCODED MENU ELEMENT NUMBERS!  REVISE WHEN MENU CHANGES!
c
	          If (menu_page .Eq. COMB_page) Then
c
c Repeat the above sequence on the high current GRT if combining GRT currents.
c  Avoid combining bolometer GRTs, as only low current is used.
c
	            If ( menu_element .Lt. 17) Then
	trans_temp = config.grttrans.grt_a_trans_temp(tranidx(menu_element))
	trans_hwid = config.grttrans.grt_a_trans_hwid(tranidx(menu_element))
	            Else
	trans_temp = config.grttrans.grt_b_trans_temp(tranidx(menu_element-16))
	trans_hwid = config.grttrans.grt_b_trans_hwid(tranidx(menu_element-16))
	            End If

	            Do j=1,4
	              Call LIB$Movc3( 2, hskp_data(lo_caloff+j*2+1,i), dl2 )
	              lo_calres(j) = dl2
	              If (dl2 .Lt. 0) Then
	                lo_calres(j) = 65536 + dl2
	              End If
	            End Do

	            Call LIB$Movc3 ( 2, hskp_data(lo_offset+1,i), dl2 )
	            lo_temp_cnts = dl2
	            If (dl2 .Lt. 0) Then
	              lo_temp_cnts = 65536 + dl2
	            End If
	            If (New_Grt .Or. lo_first) Then
	              status = Fep_Get_Grt_Conv(Config_Lun(2),lo_dbword,
	2                                       lo_Ohms,lo_Ctemps,lo_Npts)
	              If (Status .Ne. %Loc(FEP_Normal)) Then
	                Fep_Convert_To_Eng = status
	                Return
	              End If
	            End If

	            If (lo_temp_cnts .Ne. 0) Then
	              status = FEP_Counts_To_Ohms ( Config_lun(1),lo_first,
	2                                           lo_calres, lo_temp_cnts, 1,
	3                                           lo_iside, resistance, i,
	4                                           co_coeffs )
	              status = FEP_GRT_LookUp ( resistance, 0, lo_ohms,
	2                                       lo_ctemps, lo_first,lo_npts,
	3                                       lo_grt_temp, dt )
	              lo_first = .False.

	              If (status .Eq. %Loc(FEP_TblEdge)) Then
	                If (.Not. signaled) Then
	                  Call Lib$Signal(FEP_TblEdge)
	                  signaled = .True.
	                End If
	                lo_grt_temp = -99999.0
	              Else
	                If (lo_iside .Le. 2) Then
	                  icurrent = 1
	                Else
	                  icurrent = 2
	                End If
	              End If
	            Else
	              lo_grt_temp = -99999.0
	            End If            !(lo_temp_cnts .Ne. 0)
c
c Combine low/high current temperatures, except when doing bolometers.
c Always take the low current bolometer GRT.
c
	            If (menu_element .Gt. 20 .Or.
	2               (menu_element .Gt. 4 .And. menu_element .Lt. 17)) Then
	              status = FUT_Combine_HiLo(grt_temp,lo_grt_temp,
	2                                       trans_temp,trans_hwid,grt_wt)
	              If (grt_wt .Eq. -9999.) Then
	                Conv_buff(i) = -99999.0
	              Else
	                If (grt_temp .Ne. -99999.0 .And.
	2                   lo_grt_temp .Ne. -99999.0) Then
	                  conv_buff(i)=(1.-grt_wt)*lo_grt_temp + grt_wt*grt_temp
	                Else
	                  conv_buff(i) = -99999.0
	                End If
	              End If
C------------------------------------------------------------------------------
C	              If (lo_temp_cnts .Le. 15450) Then
C	                conv_buff(i) = lo_grt_temp
C	              Else If (lo_temp_cnts .Ge. 15660) Then
C	                conv_buff(i) = grt_temp
C	              Else
C                       If (grt_temp .Ne. -99999.0 .And.
C	2                   lo_grt_temp .Ne. -99999.0) Then
C	                  grt_wt = lo_temp_cnts/210 - 73.57143
C	                  conv_buff(i) = (1.0-grt_wt)*lo_grt_temp +
C	2                                      grt_wt*grt_temp
C                       Else
C                         conv_buff(i) = -99999.0
C                       End If
C	              End If
C------------------------------------------------------------------------------
	              if (conv_buff(i) .Le. 0) conv_buff(i) = -99999.0
	            Else
	              conv_buff(i) = lo_grt_temp
	            End If       !(menu_element .Gt. 20 .Or....
	          End If         !(menu_page .Eq. COMB_page)
	        Else
	          conv_buff(i) = -99999.0
	          dwell_mode(i) = 1
	        End If           !((dwell_stat_a .Eq. 0 .And....
	      Else
	        conv_buff(i) = -99999.0
	      End If             !(telm_stat .Eq. 0)
	    End Do               !i=1,nrecs

c
c No GRT conversion.
c
	  Else

	    If (menu_page .Eq. GRTA_page) Then
	      field = 'FEP_HKP.GRTA_MENU.' // GRTA_menu(menu_element)
	      status = FUT_Field_Attributes ( fr_lun, field, length, offset )
	      If (status .Ne. %Loc(FUT_Normal)) Then
	        FEP_Convert_To_Eng = status
	        Return
	      End If
	    Else If (menu_page .Eq. GRTB_page) Then
	      field = 'FEP_HKP.GRTB_MENU.' // GRTB_menu(menu_element)
	      status = FUT_Field_Attributes ( fr_lun, field, length, offset )
	      If (status .Ne. %Loc(FUT_Normal)) Then
	        FEP_Convert_To_Eng = status
	        Return
	      End If
	    Else If (menu_page .Eq. COMB_page) Then
	      off = COMB_menu_num/2
	      If (menu_element .Le. off) Then
	        If (menu_element .Gt. 4) Then
	          field = 'FEP_HKP.GRTA_MENU.' // GRTA_menu(menu_element)
	          status=FUT_Field_Attributes ( fr_lun, field, length, offset )
	          If (status .Ne. %Loc(FUT_Normal)) Then
	            FEP_Convert_To_Eng = status
	            Return
	          End If
	        End If
	        lo_field = 'FEP_HKP.GRTA_MENU.' // GRTA_menu(menu_element+off)
	        status=FUT_Field_Attributes( fr_lun,lo_field,length,lo_offset )
	        If (status .Ne. %Loc(FUT_Normal)) Then
	          FEP_Convert_To_Eng = status
	          Return
	        End If
	      Else
	        If ((menu_element - off) .Gt. 4) Then
	          field = 'FEP_HKP.GRTB_MENU.' // GRTB_menu(menu_element-off)
	          status=FUT_Field_Attributes ( fr_lun, field, length, offset )
	          If (status .Ne. %Loc(FUT_Normal)) Then
	            FEP_Convert_To_Eng = status
	            Return
	          End If
	        End If
	        lo_field = 'FEP_HKP.GRTB_MENU.' // GRTB_menu(menu_element)
	        status=FUT_Field_Attributes( fr_lun,lo_field,length,lo_offset )
	        If (status .Ne. %Loc(FUT_Normal)) Then
	          FEP_Convert_To_Eng = status
	          Return
	        End If
	      End If      !(menu_element .Le. off)

	    End If        !(menu_page .Eq. GRTA_page)
c
c Fetch the field contents.
c
	    not_avail= .False.
	    Do j=1,nrecs
c
c Check the dwell status.  Don't do anything if in dwell mode.
c
	      telm_stat = 0
	      dwell_stat_a = 0
	      dwell_stat_b = 0
	      If (hskp_data(telm_offset+1,j) .Ne. 0) Then
	        telm_stat = 1
	      End If

	      If (hskp_data(dwell_offset_a+1,j) .Lt. 0) Then
	        dwell_byte=hskp_data(dwell_offset_a+1,j)
	        null(1) = dwell_byte
	        Call mvbits(equiv,0,5,dwecont,0)
	        if (dwecont .Lt. 31) dwell_stat_a = 1
	      End If

	      If (hskp_data(dwell_offset_b+1,j) .Lt. 0) Then
	        dwell_byte=hskp_data(dwell_offset_b+1,j)
	        null(1) = dwell_byte
	        Call mvbits(equiv,0,5,dwecont,0)
	        if (dwecont .Lt. 31) dwell_stat_b = 1
	      End If

	      If (telm_stat .Eq. 0) Then
	        If ( (dwell_stat_a .Eq. 0 .And.
	2             (menu_page .Eq. GRTA_page .Or.
	3               menu_page .Eq. COMB_page .And. menu_element .Le. off))
	4           .Or.
	5            (dwell_stat_b .Eq. 0 .And.
	6             (menu_page .Eq. GRTB_page .Or.
	7               menu_page .Eq. COMB_page .And. menu_element .Gt. off)))
	8       Then
	          If (menu_page .Eq. comb_page .And.
	2             (menu_element .Le. off .And. menu_element .Lt. 5  .Or.
	3              menu_element .Gt. off .And. menu_element .Lt. 21)) Then
	            not_avail = .True.
	          End If
	          If (.Not. not_avail) Then
	            Call LIB$Movc3 ( 2, hskp_data(offset+1,j), dl2 )
	            temp_cnts = dl2
	            If (dl2 .Lt. 0) Then
	              temp_cnts = 65536 + dl2
	            End If
	            If (temp_cnts .Eq. 0) Then
	              temp_cnts = -99999.0
	            End If
	            conv_buff(j) = temp_cnts
	            dwell_mode(j) = 0
	          End If
	          If (menu_page .Eq. COMB_page) Then
	            Call LIB$Movc3 ( 2, hskp_data(lo_offset+1,j), dl2 )
	            lo_temp_cnts = dl2
	            If (dl2 .Lt. 0) Then
	              lo_temp_cnts = 65536 + dl2
	            End If
	            if (lo_temp_cnts .Eq. 0) lo_temp_cnts = -99999.0
	            If (menu_element .Gt. 20 .Or.
	2                menu_element .Gt. 4 .And. menu_element .Lt. 17) Then

	              If (lo_temp_cnts .Le. 15450) Then
	                conv_buff(j) = lo_temp_cnts
	              Else If (lo_temp_cnts .Ge. 15660) Then
	                conv_buff(j) = temp_cnts
	              Else
	                If (lo_temp_cnts .Ne. -99999.0 .And.
	2                   temp_cnts .Ne. -99999.0) Then
	                  grt_wt = lo_temp_cnts/210 - 73.57143
	                  conv_buff(j) = (1.0-grt_wt)*lo_temp_cnts +
	2                                      grt_wt*temp_cnts
	                Else
	                  conv_buff(j) = -99999.0
	                End If
	              End If
	              if (conv_buff(j) .Le. 0) conv_buff(j) = -99999.0
	            Else
	              If (lo_temp_cnts .Gt. 0) Then
	                conv_buff(j) = lo_temp_cnts
	              Else
	                conv_buff(j) = -99999.0
	              End If
	            End If
	          End If
	        Else
	          conv_buff(j) = -99999.0
	          dwell_mode(j) = 1
	        End If  !( (dwell_stat_a .Eq. 0 .And.
	      Else
	        conv_buff(j) = -99999.0
	      End If    !(telm_stat .Eq. 0)
	    End Do      !j=1,nrecs

	  End If        !(convert_flag)

	Else            !((menu_page .Eq. GRTA_page .Or.
c
c Convert via polynomial conversion or do no conversion at all.
c
c Retrieve the field name.
c
	  dwell_field = 'FEP_HKP.HSKP_MENU.DWELL_STATUS_A'

	  status = FUT_Field_Attributes ( fr_lun, dwell_field,
	2                       dwell_length, dwell_offset_a )
	  If (status .Ne. %Loc(FUT_Normal)) Then
	    FEP_Convert_To_Eng = status
	    Return
	  End If
	  dwell_field = 'FEP_HKP.HSKP_MENU.DWELL_STATUS_B'

	  status = FUT_Field_Attributes ( fr_lun, dwell_field,
	2                       dwell_length, dwell_offset_b )
	  If (status .Ne. %Loc(FUT_Normal)) Then
	    FEP_Convert_To_Eng = status
	    Return
	  End If
	  If (menu_page .Eq. HSKP_page) Then
	    field = 'FEP_HKP.HSKP_MENU.' // HSKP_menu(menu_element)
	    dbword = HSKP_db(menu_element)
	  Else If (menu_page .Eq. GRTA_page) Then
	    field = 'FEP_HKP.GRTA_MENU.' // GRTA_menu(menu_element)
	    dbword = GRTA_db(menu_element)
	  Else If (menu_page .Eq. GRTB_page) Then
	    field = 'FEP_HKP.GRTB_MENU.' // GRTB_menu(menu_element)
	    dbword = GRTB_db(menu_element)
	  Else If (menu_page .Eq. COMB_page) Then
	    If (menu_element .Le. COMB_menu_num/2) Then
	      field = 'FEP_HKP.GRTA_MENU.' // GRTA_menu(menu_element)
	      lo_field = 'FEP_HKP.GRTA_MENU.' // GRTA_menu(menu_element+off)
	      dbword = GRTA_db(menu_element)
	      lo_dbword = GRTA_db(menu_element+off)
	    Else
	      field = 'FEP_HKP.GRTB_MENU.' // GRTB_menu(menu_element-off)
	      lo_field = 'FEP_HKP.GRTB_MENU.' // GRTB_menu(menu_element)
	      dbword = GRTB_db(menu_element-off)
	      lo_dbword = GRTB_db(menu_element)
	    End If
	  Else If (menu_page .Eq. IPDU_page) Then
	    field = 'FEP_HKP.IPDU_MENU.' // IPDU_menu(menu_element)
	    dbword = IPDU_db(menu_element)
	  Else If (menu_page .Eq. VOLT_page) Then
	    field = 'FEP_HKP.VOLT_MENU.' // VOLT_menu(menu_element)
	    dbword = VOLT_db(menu_element)
	  Else If (menu_page .Eq. CURR_page) Then
	    field = 'FEP_HKP.CURR_MENU.' // CURR_menu(menu_element)
	    dbword = CURR_db(menu_element)
	  End If
c
c Get the attributes of the field being converted.
c
	  status = FUT_Field_Attributes ( fr_lun, field, length, offset )
	  If (status .Ne. %Loc(FUT_Normal)) Then
	    FEP_Convert_To_Eng = status
	    Return
	  End If
c
c Fetch the field contents.
c
	  If (length .Eq. 1) Then

	    Do j=1,nrecs
	      telm_stat = 0
	      if (hskp_data(telm_offset+1,j) .Ne. 0) telm_stat = 1
	      Call LIB$Movc3 ( length, hskp_data(offset+1,j), dl1 )
	      conv_buff(j) = dl1
	      If (dl1 .Lt. 0) Then
	        conv_buff(j) = 256 + dl1
	      End If
	      If (telm_stat .Eq. 1 ) Then
	        If (menu_page .Eq. CURR_page .And. menu_element .Eq. 33) Then
	          conv_buff(j) = conv_buff(j)
	        Else
	          conv_buff(j) = -99999.0
	        End If
	      End If
	    End Do

	  Else If (length .Eq. 2) Then

	    Do j=1,nrecs
	      dwell_mode(j) = 0
	      telm_stat = 0
	      if (hskp_data(telm_offset+1,j) .Ne. 0) telm_stat = 1
	      If (menu_page .Eq. grta_page .Or. menu_page .Eq. grtb_page
	2         .Or. menu_page .Eq. comb_page) Then
	        dwell_stat_a = 0
	        dwell_stat_b = 0
	        If (hskp_data(dwell_offset_a+1,j) .Lt. 0) Then
	          dwell_byte=hskp_data(dwell_offset_a+1,j)
	          null(1) = dwell_byte
	          Call mvbits(equiv,0,5,dwecont,0)
	          if (dwecont .Lt. 31) dwell_stat_a = 1
	        End If

	        If (hskp_data(dwell_offset_b+1,j) .Lt. 0) Then
	          dwell_byte=hskp_data(dwell_offset_b+1,j)
	          null(1) = dwell_byte
	          Call mvbits(equiv,0,5,dwecont,0)
	          if (dwecont .Lt. 31) dwell_stat_b = 1
	        End If
	        if (menu_page .Eq. comb_page) off=comb_menu_num/2
	        If (telm_stat .Eq. 0) Then
	          If ( (dwell_stat_a .Eq. 0 .And.
	2               (menu_page .Eq. GRTA_page .Or.
	3                (menu_page .Eq. COMB_page .And.
	4                        menu_element .Le. off))) .Or.
	5              (dwell_stat_b .Eq. 0 .And.
	6               (menu_page .Eq. GRTB_page .Or.
	7                (menu_page .Eq. COMB_page .And.
	8                        menu_element .Gt. off))) ) Then

	            Call LIB$Movc3 ( length, hskp_data(offset+1,j), dl2 )
	            conv_buff(j) = dl2
	            If (dl2 .Lt. 0) Then
	              conv_buff(j) = 65536 + dl2
	            End If
	            if (dl2 .Eq. 0) conv_buff(j) = -99999.0
	          Else
	            conv_buff(j) = -99999.0
	            dwell_mode(j) = 1
	          End If    !( (dwell_stat_a .Eq. 0 .And....
	        Else
	          conv_buff(j) = -99999.0
	        End If      !(telm_stat .Eq. 0)
	      Else
	        Call LIB$Movc3 ( length, hskp_data(offset+1,j), dl2 )
	        conv_buff(j) = dl2
	        If (dl2 .Lt. 0) Then
	          conv_buff(j) = 65536 + dl2
	        End If
	        if (telm_stat .Eq. 1) conv_buff(j) = -99999.0
	      End If        !(menu_page .Eq. grta_page .Or. menu_page .Eq....
	    End Do          !j=1,nrecs

	  Else
	    Call LIB$Signal ( FEP_InvDataLen )
	  End If            !(length .Eq. 1)
c
c For digital status fields, mask and return.
c
	  If ( (menu_page .Eq. HSKP_Page .And.
	2       (menu_element .Ge.  5 .And. menu_element .Le. 14 .Or.
	3        menu_element .Ge. 17 .And. menu_element .Le. 20)) .Or.
	4      (menu_page .Eq. CURR_Page .And.
	5       menu_element .Ge. 21 .And. menu_element .Le. 32) ) Then

	    Do j=1,nrecs
	      If (conv_buff(j) .Ne. -99999.0) Then
	        itemp = conv_buff(j)
	        conv_buff(j) = IIAND(itemp, digital_mask)
	      End If
C------------------------------------------------------------------------------
C                 telm_stat = 0
C                 if (hskp_data(telm_offset+1,j) .Ne. 0) telm_stat = 1
C                 if (telm_stat .Eq. 1) conv_buff(j) = -99999.0
C------------------------------------------------------------------------------
	    End Do

	  Else If (dbword .Ne. '        ' .And. convert_flag) Then
C
C Convert via polynomial conversion.
C
	    Poly_first = .True.
	    Do i=1,nrecs
	      If (conv_buff(i) .Ne. -99999.0) Then
	        status = CCT_Get_Config_Idx_Tod(timetags(i).time,1,
	2                     config_lun(3),index(3),New_poly,Stat1(3))
	        If (.Not. Status) Then
	          Call lib$signal(Fep_GetConfigErr, %val(1), %val(status))
	          Fep_Convert_to_eng = Status
	          Return
	        End If
	        If (New_poly .Or. poly_first) Then
	          Call Fep_Get_curve(Config_Lun(3),dbword,
	2                            Ndeg,poly_coeffs,Istat)
	          If (Istat .Ne. 1) Then
	            Fep_Convert_To_Eng=%loc(Fep_GetPolyErr)
	            Call lib$signal(Fep_GetpolyErr, %val(1),%val(istat))
	            Return
	          End If
	        End If
	        Poly_first = .False.
	        Call LIB$Polyf (conv_buff(i), ndeg, poly_coeffs, conv_buff(i))

	      End If        !(conv_buff(i) .Ne. -99999.0)
	    End Do          !i=1,nrecs

	  End If            !( (menu_page .Eq. HSKP_Page .And....

	End If


	FEP_Convert_To_Eng = %Loc(FEP_Normal)

	Return
	End
