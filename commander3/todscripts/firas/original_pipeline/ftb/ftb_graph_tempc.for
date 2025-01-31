	Integer*4 Function FTB_Graph_TempC(data_address, flag_address, num_recs,
	2                              start_time, end_time, data_start, ylabel)

C----------------------------------------------------------------------------
C
C      Purpose:
C
C          Fills the arrays for display using FTB_Fill_TempC and creates
C          the command lines for PLT.
C
C      Written by:
C
C          J. W. Durachta
C              ARC
C          August, 1988
C
C---------------------------------------------------------------------------
C	Revisions:
C
C	   Version 4.4  1989 May 10, SPR 3045.  Fred Shuman, STX.
C			New facility BOZO to be brought into CSDR standards
C			and renamed FTB_TempC.
C
C       version 4.4.2   SER 4180.  Add batch mode.  Aug 30,1989 D. Bouler, STX.
C                       Add retrieval of plot name from command line.
C
C	   Version 4.6  1989 Nov 21, SPR 5023.  Fred Shuman, STX.
C		Should not require the integrators to be on to plot.
C
C	   Version 4.7  1990 Mar 6, SER 5726.  Steven Alexander, STX.
C		Add capability to pass in a PLT command file.
C
C----------------------------------------------------------------------------

	Implicit None

	Integer*4	data_address
	Integer*4	flag_address
	Integer*4	num_recs
	Character*14	start_time
	Character*14	end_time
	Character*14	data_start
	Character*20	ylabel

	Integer*4	i
	Integer*4	num_pts
	Integer*4	size
	Integer*4	status
	Integer*4	FUT_Erase
	Integer*4	LIB$Get_VM
	Integer*4	LIB$Free_VM
	Integer*4	tc_address
	Integer*4	iery(9)
	Integer*4	ier, ncmd
	Integer*4	LIB$Erase_Page
	Integer*4	FTB_Fill_TempC

	Real*4		timeplotmin, timeplotmax
	Real*4		max

	Character*150	cmd(30)
	Character*38	title /' TEMPERATURE CONTROLLER STATUS - SIDE '/
	Character*1	side
	Character*66	xlabel /'Time (hours) beginning at '/
	Integer*4	xlen /26/

	Logical*4	data_present

	Integer*4       CLI$Defaulted
	Integer*4       CLI$Present
	Integer*4       UPM_Get_Value
	Integer*4       name_len
	Character*64    plot_name
	Character*64	plt_com_file

	External        CLI$Defaulted
	External        CLI$Present
	External        CLI$_Defaulted
	External        CLI$_Present
	External	FTB_NoData

C   Allocate the memory necessary for the temperature controller status words.
C   Fill the array to be plotted by PLT.

	status = LIB$Get_VM(4*5*num_recs*2, tc_address)

	status = FTB_Fill_TempC(%val(tc_address), num_recs, num_pts,
	2                       %val(flag_address), %val(data_address),
	3                       max, data_present, side)

	If (.Not. data_present) Then
	  Call LIB$Signal(FTB_NoData)
	  Return
	End If

	xlabel(xlen+1:) = data_start

C   Create the PLT commands

	Do i=1,5
	  iery(i) = 0.
	End Do

	cmd(1) = 'Plot Vertical'
	cmd(2) = 'D /VT'
	cmd(3) = 'LA OT ' // ylabel // title // side
	cmd(4) = 'LA T Temperature, Current Range, Integrator Gain,'//
	2        ' Proportional Gain'
	cmd(5) = 'LA F Timerange selected: '//start_time//' - '//end_time
	cmd(6) = 'LA X '//xlabel
	cmd(7) = 'CS 1.0'
	cmd(8) = 'V .2 .1 .9 .8'

	timeplotmin = 0.
	timeplotmax = max

	Write (cmd(9), '(a, 1x, e14.4, 1x, e14.4)')
	2                             'R', timeplotmin, timeplotmax

	cmd(10) = 'LS 1 ON 2'
	cmd(11) = 'LS 2 ON 3'
	cmd(12) = 'LS 3 ON 4'
	cmd(13) = 'LS 4 ON 5'

	cmd(14) = 'LA OX FIELDS: 1=solid; 2=dash; 3=dot_dash; 4=dot'
	ncmd = 14
C
C   Get plot name if not interactive
C
	If ( CLI$Present('INTERACTIVE') .Ne. %loc(CLI$_Present) .And.
	2    CLI$Present('INTERACTIVE') .Ne. %loc(CLI$_Defaulted) ) Then
	  status = UPM_Get_Value('PLOTDEVICE', plot_name, name_len)
	  cmd(2) = 'D ' // plot_name(1:name_len)
	End If
C
C   Get PLT command file if specified
C
	If ( CLI$Present('PLTFILE') .EQ. %loc(CLI$_Present) ) Then
	  status = UPM_Get_Value('PLTFILE', plt_com_file, name_len)
	  ncmd = ncmd + 1
	  cmd(ncmd) = '@' // plt_com_file
	End If
C
C   Make plots if not interactive
C
	If ( CLI$Present('INTERACTIVE') .Ne. %loc(CLI$_Present) .And.
	2    CLI$Present('INTERACTIVE') .Ne. %loc(CLI$_Defaulted) ) Then
	  ncmd = ncmd + 1
	  cmd(ncmd) = 'P'
	  ncmd = ncmd + 1
	  cmd(ncmd) = 'Q'
	End If

	status = LIB$Erase_Page(1,1)
	size = num_recs*2
	Call PLT (%val(tc_address), iery, size, num_pts, 5, cmd, ncmd, ier)

	status = FUT_Erase()

C   Free the display data memory.

	status = LIB$Free_VM(4*5*num_recs*2, tc_address)

	Return
	End
