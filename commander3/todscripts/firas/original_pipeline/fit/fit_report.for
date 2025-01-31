	Integer*4 Function FIT_Report ( menu_page, menu_element, nfields )

C------------------------------------------------------------------------
C    PURPOSE: 
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Fred Shuman	(adapted from FEP_Report
C            STX		 by Rob Kummerer, STX, April 24, 1987)
C            October 5, 1987
C
C    INVOCATION: STATUS = FIT_REPORT ( MENU_PAGE, MENU_ELEMENT, NFIELDS )
C
C    INPUT PARAMETERS:
C	MENU_PAGE(*)		I*2	Menu page number.
C	MENU_ELEMENT(*)		I*2	Menu element number.
C	NFIELDS			I*4	Number of fields.
C
C    OUTPUT PARAMETERS: None
C
C    SUBROUTINES CALLED: None
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
C	FIT_Menu
C	FIT_Invoc
C	FUT_Params
C	FUT_Error
C	$SSDef
C	$JPIDef
C

CH	VERSION 4.4.1 SER 3306 Q. CHUNG STX 08/10/89
CH                    PROVIDE VERSION NUMBER TO TRACK SOFTWARE UPDATE.
C----------------------------------------------------------------------

	Implicit	None

	Include         '(FIT_Menu)'
	Include         '(FIT_Invoc)'
	Include         '(FUT_Params)'
	Include         '(FUT_Error)'
	Include         '($SSDef)'
	Include         '($JPIDef)'

	Integer		*2	menu_page(*)
	Integer		*2	menu_element(*)
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

	External	FIT_Normal
	External	FIT_WriteRep

	Integer*4  Cut_Register_Version
	Integer*4  Cut_Display_Banner
	Integer*4  num_vol/80/
	Integer*4  rstatus
	Character*5 Version
	Parameter   (version='4.4.1')

c
c Write a report.
c
	If (first_time) Then

	   If (fcc_report) Then

	Rstatus = Cut_Register_Version(version)
	Rstatus = Cut_Display_Banner(fcc_report_lun,Num_vol,
	1 ' FIRAS Facility FIT_Instrument_Trends')

	      status = LIB$GetJPI (JPI$_UserName,,,,owner,)
	      Write (fcc_report_lun,5)
	      Write (fcc_report_lun,10) owner

	      Call SYS$GetTim ( current_time )

	      status = SYS$ASCTim ( time_len, time, current_time, 0 )
	      If (status .Ne. SS$_Normal) Then
	         FIT_Report = status
                 Call LIB$Signal(%Val(status))
	         Return
	      End If

	      Write (fcc_report_lun,15) time

	      Write (fcc_report_lun,20)
     .				fcc_command_line(1:fcc_command_len)

	      first_time = .False.

	   End If

	End If

c
c Write the selected plot fields to the report log.
c
	If (fcc_report .Eq. fac_present) Then

	   Write (fcc_report_lun,25) fcc_time_range

	   Do j=1,nfields

	      If (menu_page(j) .Eq. GRTA_page) Then
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

	      Write (fcc_report_lun,30) field

	   End Do

	End If


5	Format (' ',21x,'FIRAS TRENDPLOTS Report',//)
10	Format (' Run by:     ',a)
15	Format (' Run Time:   ', a)
20	Format (' Invocation: ',a,//)
25	Format (/,' Selected plot fields for time range: ',a)
30	Format ('    ',a)


	FIT_Report = %Loc(FIT_Normal)

	Return
	End
