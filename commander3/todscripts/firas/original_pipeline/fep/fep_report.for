	Integer*4 Function FEP_Report ( menu_page, menu_element,
	2                               convert_flag, nfields )

C------------------------------------------------------------------------
C    PURPOSE:
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX
C            April 24, 1987
C
C    INVOCATION: STATUS = FEP_REPORT ( MENU_PAGE, MENU_ELEMENT, CONVERT_FLAG,
C				       NFIELDS )
C
C    INPUT PARAMETERS:
C	MENU_PAGE(*)		I*2	Menu page number.
C	MENU_ELEMENT(*)		I*2	Menu element number.
C	CONVERT_FLAG(*)		L*1	Convertion flag.
C	NFIELDS			I*4	Number of fields.
C
C    OUTPUT PARAMETERS: None
C
C    SUBROUTINES CALLED: None
C
C    COMMON VARIABLES USED:
C	HSKP_MENU(28)		Ch*32		Housekeeping header menu.
C	GRTA_MENU(32)		Ch*32		GRT side A menu.
C	GRTB_MENU(32)		Ch*32		GRT side B menu.
C	COMB_MENU(32)		Ch*32		Combined GRT menu.
C	IPDU_MENU(14)		Ch*32		IPDU temps menu.
C	VOLT_MENU(20)		Ch*32		Voltages menu.
C	CURR_MENU(34)		Ch*32		Currents menu.
C
C    INCLUDE FILES:
C	FEP_Menu
C	FEP_Invoc
C	FUT_Params
C	FUT_Error
C	$SSDef
C	$JPIDef
C
C----------------------------------------------------------------------
C Changes:
C	1988 Nov 1  Break up the writing of fcc_command_line to avoid
C	            'record overflow error'.  This was occurring when
C	            the cmd string caused an overflow of the 132-char
C	            limit of the 'a' format.  (SPR 2605)  F. Shuman, STX
C----------------------------------------------------------------------

	Implicit	None

	Include         '(FEP_Menu)'
	Include         '(FEP_Invoc)'
	Include         '(FUT_Params)'
	Include         '(FUT_Error)'
	Include         '($SSDef)'
	Include         '($JPIDef)'

	Integer		*2	menu_page(*)
	Integer		*2	menu_element(*)
	Logical		*1	convert_flag(*)
	Integer		*2	nfields

	Integer         *4	status
	Integer         *4      i
	Integer         *4      j
	Integer         *4      k
	Logical		*1	first_time/.True./

	Integer		*4	current_time(2)
	Integer         *2      time_len
	Character       *32     time

	Character	*8      owner
	Character       *9      day
	Character	*64	field

	Integer		*4	LIB$GetJPI
	Integer		*4	SYS$GetTim
	Integer		*4	SYS$ASCTim

	External	FEP_Normal
	External	FEP_WriteRep

c
c Write a report.
c
	If (first_time) Then

	   If (fcc_report) Then

	      status = LIB$GetJPI (JPI$_UserName,,,,owner,)
	      Write (fcc_report_lun,5)
	      Write (fcc_report_lun,10) owner

	      Call SYS$GetTim ( current_time )

	      status = SYS$ASCTim ( time_len, time, current_time, 0 )
	      If (status .Ne. SS$_Normal) Then
	         FEP_Report = status
	         Call LIB$Signal(%Val(status))
	         Return
	      End If

	      Write (fcc_report_lun,15) time

	      Write (fcc_report_lun,20) fcc_command_line(1:64)
	      Do j=65,fcc_command_len,64
	         Write (fcc_report_lun,21) fcc_command_line(j:j+63)
	      End Do

	      first_time = .False.

	   End If

	End If

c
c Write the selected plot fields to the report log.
c
	If (fcc_report .Eq. fac_present) Then

	   Write (fcc_report_lun,25) fcc_time_range

	   Do j=1,nfields

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

	      Write (fcc_report_lun,30) field

	   End Do

	End If


5	Format (' ',21x,'FIRAS ENGPLOTS Report',//)
10	Format (' Run by:     ',a)
15	Format (' Run Time:   ',a)
20	Format (' Invocation: ',a)
21	Format ('             ',a)
25	Format (///,' Selected plot fields for time range: ',a)
30	Format ('    ',a)


	FEP_Report = %Loc(FEP_Normal)

	Return
	End
