	Integer*4 Function FIT_Get_Command ()

C------------------------------------------------------------------------
C    PURPOSE: Parse FIT_Instrument_Trends invocation line.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Fred Shuman
C            STX
C            September 11, 1987		revision:  1987 Nov 21
C
C    INVOCATION: status = FIT_GET_COMMAND ()
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS: None
C
C    SUBROUTINES CALLED:
C	LIB$Date_Time
C	Time_LT
C	SYS$BinTim
C	UPM_Present
C	UPM_Get_Value
C	UPM_Get_Float
C	LIB$Signal
C	CT_Binary_To_GMT
C
C    COMMON VARIABLES USED:
C	FCC_INTERACTIVE		 I*4		Interactive/Batch flag.
C	FCC_REPORT		 I*4		Report flag.
C	FCC_REPORT_FILE		Ch*64		Report file name.
C	FCC_PLOT_DEVICE		Ch*32		Initial PLT device type.
C	FCC_PLOT_DEV_SELECT	 I*4		Select flag.
C	FCC_PLT_COM_FILE	Ch*64		PLT command file name.
C	FCC_PLT_COM		 I*4		Select flag
C	FCC_FIELDS(10)		Ch*32		Batch STATS field selection.
C	FCC_FIELDS_LEN(10)	 I*2		Field lengths.
C	FCC_FIELDS_COUNT	 I*4		Number of selected fields.
C	FCC_JSTART_TIME		Ch*14		Batch starttime value.
C	FCC_JSTART(2)		 I*4		VAX internal representation.
C	FCC_JSTART_SELECT	 I*4		Select flag.
C	FCC_JSTOP_TIME		Ch*14		Batch stoptime value.
C	FCC_JSTOP(2)		 I*4		VAX internal representation.
C	FCC_JSTOP_SELECT	 I*4		Select flag.
C	FCC_ORBIT		 R*4		Orbital period in min.
C       FCC_ORBIT_SELECT  	 I*4            Orbital period select flag.
C	FCC_BINSIZE		 R*4		Bin size value.
C       FCC_BINORB_SELECT  	 I*4            Bin size is in orbital periods.
C	FCC_BINSIZE_ORBIT	 R*4		Bin size in orbits.
C	FCC_BINSIZE_LENGTH	Ch*7		Bin size in DDDHHMM.
C	FCC_BINSIZE_SELECT	 I*4		Select flag.
C
C    INCLUDE FILES:
C	FUT_Error
C	FIT_Invoc
C	FUT_Params
C	UPM_Stat_Msg
C	CTUser.Inc
C	$SSDEF
C
C----------------------------------------------------------------------
C	Revised:
C
C	   SPR 3285.  Get COBE orbital period from FUT_Orbital_Period.
C	              This module and FIT, FIT_Display_Menu, and
C	              FIT_Display_Main were revised.
C	              Fred Shuman, STX  1989 Jul 20.
C
C          SPRs 4356, 4358.  Add orbital period qualifier to command line.
C                     Get start and stop time from command line if /NOINT
C                     qualifier is specified.  Compute binsize after getting
C                     orbital period.  Modified FIT, FIT_Get_Command, FIT.CLD,
C                     FIT_INVOC.TXT.
C                     David Bouler, STX  1989 Aug 25.
C                     Fred Shuman,  STX  1989 Aug 29.
C	
C	   SER 5727.  Add capability to pass in a command file to PLT. 
C		      Added qualifier PLTFILE, revised FIT_PLOT, this module,
C		      FIT.CLD, and FIT_INVOC.TXT.
C		      Steven Alexander, STX  1990 Feb 28.
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FIT_Invoc)'
	Include		'(FUT_Params)'
	Include		'(UPM_Stat_Msg)'
	Include		'CT$LIBRARY:CTUser.Inc'
	Include		'($SSDEF)'

	Integer		*4	i
	Integer		*4	j
	Integer		*4	k
	Integer		*4	num_qual/9/	!number of qualifiers
	Character	*32	qualifier(9)	!qualifiers
	Character	*32	keyword(9,5)    ! number of keywords
	Integer		*4	num_key(7)	!number of keywords per qualifier
	Integer		*4	ljstart		!length of start time
	Integer		*4	ljstop		!length of stop time
	Integer		*4	status		!return status from function
	Integer		*2	len		!length of string from command line

	Real		*4	dur		!to convert DDDHHMM to binsize
	Integer		*4	time(2)
	Character	*14	jstart_default	!default /JSTART
	Character	*14	jstop_default	!default /JSTOP
	Character	*23	date

	Integer		*4	LIB$Date_Time
	Logical		*2      Time_LT
	Integer		*4	SYS$BinTim
	Integer		*4	UPM_Present
	Integer		*4	UPM_Get_Value
	Integer		*4	UPM_Get_Float
	Integer		*4	OTS$Cvt_T_F

	External	FIT_Normal
	External	FIT_InvTim
	External	FIT_InvalBin

	Data jstart_default/'85001000000000'/
	Data jstop_default/'99365235959990'/

c
c Fetch the invocation.
c
	FIT_Get_Command = %Loc(FIT_Normal)

	status = UPM_Get_Value ( '$LINE', fcc_command_line,
	2                         fcc_command_len )

	If (status .Ne. SS$_Normal) Then
	   FIT_Get_Command = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

c
c Assign qualifiers.
c
	qualifier(1) = 'Interactive'
	qualifier(2) = 'Report'
	qualifier(3) = 'Plotdevice'
	qualifier(4) = 'Fields'
	qualifier(5) = 'JStart'
	qualifier(6) = 'JStop'
	qualifier(7) = 'Binsize'
	qualifier(8) = 'Orbital_pd'
	qualifier(9) = 'Pltfile'

c
c Assign keywords.
c
	Do i=1,num_qual
	   num_key(i) = 0
	End Do

	num_key(7) = 5
	keyword(7,1) = 'Binsize.Orbit'
	keyword(7,2) = 'Binsize.Daily'
	keyword(7,3) = 'Binsize.Monthly'
	keyword(7,4) = 'Binsize.Mission'
	keyword(7,5) = 'Binsize.Length'

c
c Parse the command line.
c
	Do i=1,num_qual

	   status = UPM_Present ( qualifier(i) )

	   If (status .Eq. UPM_Pres) Then

c
c	If it's in the command line, set the proper flags.
c
	      If (i .Eq. 1) Then
	         fcc_interactive = fac_present

	      Else If (i .Eq. 2) Then
	         fcc_report = fac_present
	         status = UPM_Get_Value(qualifier(i),fcc_report_file,
	2                               len)
	         If (status .Ne. SS$_Normal) Then
	            FIT_Get_Command = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If

	      Else If (i .Eq. 3) Then
	         fcc_plot_dev_select = fac_present
	         status = UPM_Get_Value(qualifier(i),fcc_plot_device,
	2                               len)
	         If (status .Eq. UPM_Absent) Then

	            If (fcc_interactive .Eq. fac_present) Then
	               fcc_plot_dev_select = fac_not_present
	            Else
	               fcc_plot_device = '/QMS'
	            End If

	         Else

	            If (status .Ne. SS$_Normal) Then
	               FIT_Get_Command = status
	               Call LIB$Signal(%Val(status))
	               Return
	            End If

	         End If

	      Else If (i .Eq. 4) Then
	         fcc_fields_count = 0
	         status = UPM_Comma
	         Do While (status .Eq. UPM_Comma .And. k .lt. 10)
	           fcc_fields_count = fcc_fields_count + 1
	           status = UPM_Get_Value(qualifier(i),
	2                       fcc_fields(fcc_fields_count),
	3                       fcc_fields_len(fcc_fields_count))
	         End Do
	         If (status .Ne. SS$_Normal) Then
	            FIT_Get_Command = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If

	      Else If (i .Eq. 5) Then
	         fcc_jstart_select = fac_present
	         status = UPM_Get_Value(qualifier(i),fcc_jstart_time,
	2                               len)
	         If (status .Eq. UPM_Absent) Then
c
c		Assign the start time to be the start of the current day
c
	            status = LIB$Date_Time(date)
	            If (status .Ne. SS$_Normal) Then
	               FIT_Get_Command = status
	               Call LIB$Signal(%Val(status))
	               Return
	            End If
	            date = date(1:12) // '00:00:00.00'
	            status = SYS$BinTim(date,time)
	            If (status .Ne. SS$_Normal) Then
	               FIT_Get_Command = status
	               Call LIB$Signal(%Val(status))
	               Return
	            End If
	            Call CT_Binary_To_GMT(time,fcc_jstart_time)
	            fcc_jstart_select = fac_not_present

	         Else

	            If (status .Ne. SS$_Normal) Then
	               FIT_Get_Command = status
	               Call LIB$Signal(%Val(status))
	               Return
	            End If

	         End If

	      Else If (i .Eq. 6) Then
	         fcc_jstop_select = fac_present
	         status = UPM_Get_Value(qualifier(i),fcc_jstop_time,
	2                               len)
	         If (status .Eq. UPM_Absent) Then
c
c		Assign the stop time to be the end of the current day
c
	            status = LIB$Date_Time(date)
	            If (status .Ne. SS$_Normal) Then
	               FIT_Get_Command = status
	               Call LIB$Signal(%Val(status))
	               Return
	            End If
	            date = date(1:12) // '23:59:59.99'
	            status = SYS$BinTim(date,time)
	            If (status .Ne. SS$_Normal) Then
	               FIT_Get_Command = status
	               Call LIB$Signal(%Val(status))
	               Return
	            End If
	            Call CT_Binary_To_GMT(time,fcc_jstop_time)
	            fcc_jstop_select = fac_not_present

	         Else

	            If (status .Ne. SS$_Normal) Then
	               FIT_Get_Command = status
	               Call LIB$Signal(%Val(status))
	               Return
	            End If

	         End If

	      Else If (i .Eq. 7) Then
	         fcc_binsize_select = fac_present
	         Do j=1,num_key(i)
	            status = UPM_Present(keyword(i,j))
	            If (status .Eq. UPM_Pres) Then
	               If (j .Eq. 1) Then
	                  fcc_binorb_select = fac_present
	                  status = UPM_Get_Float(keyword(i,j),
	2                               fcc_binsize_orbit,len)
	                  If (status .Ne. SS$_Normal) Then
	                     FIT_Get_Command = status
	                     Call LIB$Signal(%Val(status))
	                     Return
	                  End If
	               Else If (j .Eq. 2) Then
	                  fcc_binsize = fac_day
	               Else If (j .Eq. 3) Then
	                  fcc_binsize = fac_month
	               Else If (j .Eq. 4) Then
	                  fcc_binsize = fac_mission
	               Else If (j .Eq. 5) Then
	                  status = UPM_Get_Value(keyword(i,j),
	2                               fcc_binsize_length,len)
	                  If (status .Ne. SS$_Normal) Then
	                     FIT_Get_Command = status
	                     Call LIB$Signal(%Val(status))
	                     Return
	                  End If
	                  status = OTS$Cvt_T_F(fcc_binsize_length(1:3),
	2                                      dur,,,%val(19))
	                  fcc_binsize = dur*fac_day
	                  status = OTS$Cvt_T_F(fcc_binsize_length(4:5),
	2                                      dur,,,%val(19))
	                  fcc_binsize = fcc_binsize + dur*fac_hour
	                  status = OTS$Cvt_T_F(fcc_binsize_length(6:7),
	2                                      dur,,,%val(19))
	                  fcc_binsize = fcc_binsize + dur*fac_minute

	               End If

	               If (fcc_binsize .Lt. fcc_orbit .Or.
	2                  fcc_binsize .Gt. fac_mission) Then
	                  FIT_Get_Command = %Loc(FIT_InvalBin)
	               End If

	            End If
	         End Do

	      Else If (i .Eq. 8) Then
	         fcc_orbit_select = fac_present
	         status = UPM_Get_float(qualifier(i),fcc_orbit)

	      Else If (i .Eq. 9) Then
	         fcc_plt_com = fac_present
	         status = UPM_Get_Value(qualifier(i),fcc_plt_com_file,
	2                               len)
	         If (status .Ne. SS$_Normal) Then
	            FIT_Get_Command = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If

	      End If

	   Else If (status .Eq. UPM_Defaulted) Then
c
c	If something's been defaulted, then set proper flag.
c
	      If (i .Eq. 1) Then
	         fcc_interactive = fac_present

	      Else If (i .Eq. 2) Then
	         fcc_report = fac_present
	         status = UPM_Get_Value(qualifier(i),fcc_report_file,
	2                               len)
	         If (status .Ne. SS$_Normal) Then
	            FIT_Get_Command = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If

	      Else If (i .Eq. 4) Then
	         fcc_fields_count = 0
	         status = UPM_Comma
	         Do While (status .Eq. UPM_Comma .And. k .lt. 10)
	           fcc_fields_count = fcc_fields_count + 1
	           status = UPM_Get_Value(qualifier(i),
	2                       fcc_fields(fcc_fields_count),
	3                       fcc_fields_len(fcc_fields_count))
	         End Do
	         If (status .Ne. SS$_Normal) Then
	            FIT_Get_Command = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If

	      Else If (i .Eq. 5) Then
c
c		Set default start time to 00:00 of current day
c
	         status = LIB$Date_Time(date)
	         If (status .Ne. SS$_Normal) Then
	            FIT_Get_Command = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If
	         date = date(1:12) // '00:00:00.00'
	         status = SYS$BinTim(date,time)
	         If (status .Ne. SS$_Normal) Then
	            FIT_Get_Command = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If
	         Call CT_Binary_To_GMT(time,fcc_jstart_time)
	         fcc_jstart_select = fac_not_present

	      Else If (i .Eq. 6) Then
c
c		Set default stop time to 23:59 of current day
c
	         status = LIB$Date_Time(date)
	         If (status .Ne. SS$_Normal) Then
	            FIT_Get_Command = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If
	         date = date(1:12) // '23:59:59.99'
	         status = SYS$BinTim(date,time)
	         If (status .Ne. SS$_Normal) Then
	            FIT_Get_Command = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If
	         Call CT_Binary_To_GMT(time,fcc_jstop_time)
	         fcc_jstop_select = fac_not_present

	      Else If (i .Eq. 7) Then
	         fcc_binorb_select = fac_present

	      End If

	   Else If (status .Eq. UPM_Negated) Then
c
c	If something's been negated, then set the proper flag.
c
	      If (i .Eq. 1) Then
	         fcc_interactive = fac_not_present

	      Else If (i .Eq. 2) Then
	         fcc_report = fac_not_present

	      Else If (i .Eq. 4) Then
	         fcc_fields_count = 0

	      Else If (i .Eq. 5) Then
	         fcc_jstart_select = fac_not_present

	      Else If (i .Eq. 6) Then
	         fcc_jstop_select = fac_not_present

	      Else If (i .Eq. 7) Then
	         fcc_binsize_select = fac_not_present

	      End If

	   Else
c
c	Assuming no defaults, set the processing flags.
c
	      If (i .Eq. 1) Then
	         fcc_interactive = fac_present

	      Else If (i .Eq. 3) Then
	         fcc_plot_dev_select = fac_not_present

	      Else If (i .Eq. 4) Then
	         fcc_fields_count = 0

	      Else If (i .Eq. 5) Then
c
c		Set default start time to 00:00 of current day
c
	         status = LIB$Date_Time(date)
	         If (status .Ne. SS$_Normal) Then
	            FIT_Get_Command = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If
	         date = date(1:12) // '00:00:00.00'
	         status = SYS$BinTim(date,time)
	         If (status .Ne. SS$_Normal) Then
	            FIT_Get_Command = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If
	         Call CT_Binary_To_GMT(time,fcc_jstart_time)
	         fcc_jstart_select = fac_not_present

	      Else If (i .Eq. 6) Then
c
c		Set default stop time to 23:59 of current day
c
	         status = LIB$Date_Time(date)
	         If (status .Ne. SS$_Normal) Then
	            FIT_Get_Command = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If
	         date = date(1:12) // '23:59:59.99'
	         status = SYS$BinTim(date,time)
	         If (status .Ne. SS$_Normal) Then
	            FIT_Get_Command = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If
	         Call CT_Binary_To_GMT(time,fcc_jstop_time)
	         fcc_jstop_select = fac_not_present

	      Else If (i .Eq. 7) Then
	         fcc_binsize_select = fac_not_present
                 fcc_binorb_select = fac_present

	      Else If (i .Eq. 8) Then
	         fcc_orbit_select = fac_not_present

	      Else If (i .Eq. 9) Then
	         fcc_plt_com = fac_not_present

	      End If

	   End If

	End Do

c
c Setup the timerange variables.
c
	ljstart = index(fcc_jstart_time,' ')-1
	If (ljstart .Eq. -1)ljstart=14
	ljstop = index(fcc_jstop_time,' ')-1
	If (ljstop .Eq. -1) ljstop=14

	fcc_jstart_time = fcc_jstart_time(1:ljstart) //
	2                    jstart_default(ljstart+1:)
	fcc_jstop_time = fcc_jstop_time(1:ljstop) //
	2                    jstop_default(ljstop+1:)

	Call CT_GMT_To_Binary(fcc_jstart_time,fcc_jstart)
	Call CT_GMT_To_Binary(fcc_jstop_time,fcc_jstop)

	If (Time_LT(fcc_jstop,fcc_jstart)) Then
	   FIT_Get_Command = %Loc(FIT_InvTim)
	   Call LIB$Signal(FIT_InvTim)
	   Return
	End If

	Return
	End
