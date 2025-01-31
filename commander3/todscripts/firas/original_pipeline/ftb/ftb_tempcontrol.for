	Program FTB_TempControl

C---------------------------------------------------------------------------
C
C      Purpose:
C
C          To produce a graphic display of the temperature controller status 
C          bits in the status monitor.
C
C      Written by:
C
C          J. W. Durachta
C              ARC
C          December, 1988
C
C---------------------------------------------------------------------------
C	Revisions:
C
C	   Version 4.4  1989 May 10, SPR 3045.  Fred Shuman, STX.
C		New facility BOZO to be brought into CSDR standards
C		and renamed FTB_TempC.
C
C	   Version 4.6  1989 Nov 21, SPR 5023.  Fred Shuman, STX.
C		Should not require the integrators to be on to plot.
C
C---------------------------------------------------------------------------


	Implicit None

C   FTB messages

	External	FTB_Normal
	External	FTB_Aberr

C   External Functions

	External	CLI$_Defaulted
	External	CLI$_Present

C   Include file

	Include		'($SSdef)'

C   Local variables and functions

	Integer*4	i
	Integer*4	status
	Integer*4	CLI$Present
	Integer*4	CLI$Get_Value
	Integer*4	fetch_type
	Integer*4	data_address
	Integer*4	num_recs
	Integer*4	flag_address
	Character*14	start_time
	Character*14	end_time
	Character*14	data_start
	Character*20	ylabel

	Integer*4	FTB_Fetch_Data
	Integer*4	FTB_Graph_TempC
	Integer*4	UPM_Get_Longword
	Integer*4	LIB$Free_VM

	Integer*4	on
	Integer*4	enable
	Integer*4	mjfm
	Integer*4	stat_mon
	Integer*4	start_bit

	Logical*1	proceed / .true. /

C   FTB Parameters

	Integer*4	Word31
	Integer*4	TempC
	Parameter	(Word31 = 1)
	Parameter	(TempC = 2)

	Common /tc_info/ on, enable, mjfm, stat_mon, start_bit

C
C   Retrieve the data from COBETRIEVE.
C
	fetch_type = TempC
	status = FTB_Fetch_Data(fetch_type, data_address, num_recs,
	2                       flag_address, start_time, end_time, data_start)

	If (status .Ne. 1) Then
	  proceed = .False.
	End If

C   Determine the fields for display

	If ( proceed ) Then

	  If (CLI$Present('ICAL') .Eq. %loc(CLI$_Present)) Then

	    ylabel = 'INTERNAL CALIBRATOR'
	    mjfm = 1
	    stat_mon = 2
	    start_bit = 0
	    on = 7
	    enable = 2

	  Else If (CLI$Present('XCAL') .Eq. %loc(CLI$_Present)) Then

	    ylabel = 'EXTERNAL CALIBRATOR'
	    mjfm = 2
	    stat_mon = 3
	    start_bit = 2
	    on = 4
	    enable = 3

	  Else If (CLI$Present('SKYHORN') .Eq. %loc(CLI$_Present)) Then

	    ylabel = 'SKY HORN'
	    mjfm = 2
	    stat_mon = 2
	    start_bit = 0
	    on = 5
	    enable = 2

	  Else If (CLI$Present('REFHORN') .Eq. %loc(CLI$_Present)) Then

	    ylabel = 'REFERENCE HORN'
	    mjfm = 1
	    stat_mon = 3
	    start_bit = 2
	    on = 6
	    enable = 3

	  Else

	    proceed = .False.

	  End If

	End If  !proceed

C   Graph the status bits

	If (proceed) Then
	  status = FTB_Graph_TempC(data_address, flag_address, num_recs,
	2                          start_time, end_time, data_start, ylabel)
	End If

C   Deallocate the HKP and time flag memory.

	status = LIB$Free_VM(576*num_recs, data_address)
	status = LIB$Free_VM(4*num_recs, flag_address)

C   Exit the program.

	If ( proceed ) Then
	  Call LIB$Signal(FTB_Normal)
	  Call Exit(SS$_Normal)
	Else
	  Call LIB$Signal(FTB_Aberr)
	  Call Exit(SS$_Abort)
	End If

	End
