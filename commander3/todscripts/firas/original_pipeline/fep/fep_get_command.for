	Integer*4 Function FEP_Get_Command ()

C------------------------------------------------------------------------------
C    PURPOSE: Parse FEP_ENGPLOTS invocation line.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX
C            March 9, 1987
C
C    INVOCATION: status = FEP_GET_COMMAND ()
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
C	LIB$Signal
C	Sys$Gettim
C	CT_Binary_To_GMT
C
C    COMMON VARIABLES USED:
C	FCC_COMMAND_LINE	Ch*256		Invocation of ENGplots.
C	FCC_COMMAND_LEN		 I*2		Length of invoc string.
C	FCC_INTERACTIVE		 I*4		Interactive/Batch flag.
C	FCC_REPORT		 I*4		Report flag.
C	FCC_REPORT_FILE		Ch*64		Report file name.
C	FCC_REPORT_LUN		 I*4		Report file lun.
C	FCC_CONVERT		 I*4		Convert-to-EU flag.
C	FCC_MONITOR		 I*4		Monitor field flag.
C	FCC_PLOT_DEVICE		Ch*32		Initial PLT device type.
C	FCC_PLOT_DEV_SELECT	 I*4		Select flag.
C	FCC_PLOT_COM_FILE	Ch*32		put the Plt command file
C	FCC_PLOT_COM    	 I*4		Select Command file
C	FCC_FIELDS(4)		Ch*32		Batch HSKP field selection.
C	FCC_FIELDS_LEN(4)	 I*2		Field lengths.
C	FCC_FIELDS_COUNT	 I*4		Number of selected fields.
C	FCC_SMOOTH		 I*4		Smoothing select flag.
C	FCC_WINDOW		 I*4		Width of smoothing window.
C	FCC_CROSSGAPS		 I*4		Flag(smooth to cross datagaps?).
C	FCC_JSTART_SELECT	 I*4		Select flag.
C	FCC_JSTOP_SELECT	 I*4		Select flag.
C	FCC_JSTART_TIME		Ch*14		Batch starttime select.
C	FCC_JSTART(2)		 I*4		VAX internal representation.
C	FCC_JSTOP_TIME		Ch*14		Batch stoptime select.
C	FCC_JSTOP(2)		 I*4		VAX internal representation.
C	FCC_TIME_RANGE		Ch*30		CT timerange select.
C	FCC_AVERAGE_CAL_COUNTS	 I*4		Average calibration counts.
C
C    INCLUDE FILES:
C	FUT_Error
C	FEP_Invoc
C	FUT_Params
C	UPM_Stat_Msg
C	CTUser.Inc
C	$SSDEF
C
C------------------------------------------------------------------------------
C    Revisions:
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
C    07-Aug-1989
C	Added the qualifier 'Average_calibration' to allow the user to choose
C	whether or not to average the calibration resistor counts prior to 
C	computing the counts to ohms conversion coefficients.
C	Don Stevens-Rayburn, ARC.
C
C       SER 5725, Add the capability to access the plt command file from
C                 CLD.
C                 H. Wang, STX, Mar. 2, 1990
C
C	SPR 4171, Standardize report file names.  Larry P. Rosen, STX,
C		July-August 1990
C------------------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include 	'(FEP_Invoc)'
	Include 	'(FUT_Params)'
	Include 	'(UPM_Stat_Msg)'
	Include 	'CT$LIBRARY:CTUser.Inc'
	Include 	'($SSDEF)'

	Integer 	*4	i
	Integer 	*4	j
	Integer 	*4	k
	Character	*32	qualifier(12)	!qualifiers
c	Character	*32	keyword(12,3)
	Integer 	*4	num_qual/12/	!number of qualifiers
c	Integer 	*4	num_key(12)	!number of keywords
	Integer 	*4	ljstart		!length of start time
	Integer 	*4	ljstop		!length of stop time
	Integer 	*4	status		!return status from function
	Integer		*4	ret_status	!return status from other fn.
	Integer 	*2	len		!length of string from cmd line
	Integer		*4	irem		!check window size

	Integer 	*4	time(2)
	Character	*14	jstart_default	!default /JSTART
	Character	*14	jstop_default	!default /JSTOP
	Character	*23	date
	Character	*14	CurGMT		! current run time in GMT

	Integer 	*4	LIB$Date_Time
	Logical 	*2      Time_LT
	Integer 	*4	SYS$BinTim
	Integer		*4	Sys$Gettim
	Integer		*2	CurTime(2)
	Integer 	*4	UPM_Present
	Integer 	*4	UPM_Get_Value
	Integer 	*4	UPM_Get_LongWord
	Character	*33	Report_Default

	External	FEP_Normal
	External	FEP_InvTim
	External	FEP_Window

	Data jstart_default/'85001000000000'/
	Data jstop_default/'99365235959999'/

c
c Fetch the invocation.
c
	FEP_Get_Command = %Loc(FEP_Normal)

	status = UPM_Get_Value ( '$LINE', fcc_command_line,
     .				  fcc_command_len )

	If (status .Ne. SS$_Normal) Then
	   FEP_Get_Command = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If
c
c Assign qualifiers.
c
	qualifier(1) = 'Interactive'
	qualifier(2) = 'Pltfile'
	qualifier(3) = 'Convert'
	qualifier(4) = 'Monitor'
	qualifier(5) = 'Plotdevice'
	qualifier(6) = 'Fields'
	qualifier(7) = 'JStart'
	qualifier(8) = 'JStop'
	qualifier(9) = 'Smooth'
	qualifier(10)= 'Crossgaps'
	qualifier(11)= 'Average_calibration'
	qualifier(12)= 'Report'
c
c Get current run time and set default file name.
c
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
                 fcc_plt_com = fac_present
                 status = UPM_Get_Value(qualifier(i),fcc_plt_com_file,
     .					len)
		 If (status .Ne. SS$_Normal) Then
		    FEP_Get_Command = status
		    Call LIB$Signal(%Val(status))
		    Return
		 End If

	      Else If (i .Eq. 3) Then
	         fcc_convert = fac_present

	      Else If (i .Eq. 4) Then
	         fcc_monitor = fac_present

	      Else If (i .Eq. 5) Then
		 fcc_plot_dev_select = fac_present
	         status = UPM_Get_Value(qualifier(i),fcc_plot_device,
     .					len)
	         If (status .Eq. UPM_Absent) Then

	            If (fcc_interactive .Eq. fac_present) Then
		       fcc_plot_dev_select = fac_not_present
		    Else
		       fcc_plot_device = '/QMS'
	            End If

		 Else

		    If (status .Ne. SS$_Normal) Then
		       FEP_Get_Command = status
		       Call LIB$Signal(%Val(status))
		       Return
		    End If

	         End If

	      Else If (i .Eq. 6) Then
                 fcc_fields_count = 0
		 status = UPM_Comma
		 Do While (status .Eq. UPM_Comma .And. fcc_fields_count .Lt. 4)
                   fcc_fields_count = fcc_fields_count + 1
                   status = UPM_Get_Value(qualifier(i),
     .				fcc_fields(fcc_fields_count),
     .				fcc_fields_len(fcc_fields_count))
                 End Do
		 If (status .Ne. SS$_Normal) Then
		    FEP_Get_Command = status
		    Call LIB$Signal(%Val(status))
		    Return
		 End If

	      Else If (i .Eq. 7) Then
		 fcc_jstart_select = fac_present
	         status = UPM_Get_Value(qualifier(i),fcc_jstart_time,
     .					len)
	         If (status .Eq. UPM_Absent) Then
c
c   Assign the start time to be the start of the current day
c
		      status = LIB$Date_Time(date)
		      If (status .Ne. SS$_Normal) Then
		         FEP_Get_Command = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      date = date(1:12) // '00:00:00.00'
		      status = SYS$BinTim(date,time)
		      If (status .Ne. SS$_Normal) Then
		         FEP_Get_Command = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      Call CT_Binary_To_GMT(time,fcc_jstart_time)
		      fcc_jstart_select = fac_not_present

		 Else

		    If (status .Ne. SS$_Normal) Then
		       FEP_Get_Command = status
		       Call LIB$Signal(%Val(status))
		       Return
		    End If

	         End If

	      Else If (i .Eq. 8) Then
		 fcc_jstop_select = fac_present
	         status = UPM_Get_Value(qualifier(i),fcc_jstop_time,
     .					len)
	         If (status .Eq. UPM_Absent) Then
c
c   Assign the stop time to be the end of the current day
c
		      status = LIB$Date_Time(date)
		      If (status .Ne. SS$_Normal) Then
		         FEP_Get_Command = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      date = date(1:12) // '23:59:59.99'
		      status = SYS$BinTim(date,time)
		      If (status .Ne. SS$_Normal) Then
		         FEP_Get_Command = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      Call CT_Binary_To_GMT(time,fcc_jstop_time)
		      fcc_jstop_select = fac_not_present

		 Else

		    If (status .Ne. SS$_Normal) Then
		       FEP_Get_Command = status
		       Call LIB$Signal(%Val(status))
		       Return
		    End If

	         End If

	      Else If (i .Eq. 9) Then

		 fcc_smooth = fac_present
	         status = UPM_Get_LongWord(qualifier(i),fcc_window)
	         If (status .Eq. UPM_Absent) Then
		    fcc_window = 1
		 End If

	         irem = Mod(fcc_window,2)

	         If (.Not. (irem .Eq. 1 .And. fcc_window .Ge. 0) ) Then
		    FEP_Get_Command = %Loc(FEP_Window)
	            Call LIB$Signal ( FEP_Window )
		    Return
	         End If

	      Else If (i .Eq. 10) Then
		 fcc_crossgaps = fac_present

	      Else If (i .Eq. 11) Then
	         fcc_average_cal_counts = fac_present

	      Else If (i .Eq. 12) Then
                 fcc_report = fac_present
                 Status = Sys$Gettim( Curtime )
                 Call CT_Binary_To_Gmt( Curtime, Curgmt )
                 Report_Default = 'FEP_' // FCC_Jstart_Time(1:7) // '_' // 
     1              FCC_Jstop_Time(1:7) // '.REP_' // Curgmt(1:9)
                 status = UPM_Get_Value(qualifier(i),fcc_report_file,
     .					len)
                 If (status .EQ. upm_absent) Then
                    fcc_report_file = Report_Default
                 End If
 		write(6,222)fcc_report_file
 222		format(1x,'The file name is generated ',a40)
	      End If

	   Else If (status .Eq. UPM_Defaulted) Then
c
c   If something's been defaulted, then set proper flag.
c
	      If (i .Eq. 1) Then
	         fcc_interactive = fac_present

              Else If (i .Eq. 2) Then
                 fcc_plt_com = fac_present

	      Else If (i .Eq. 3) Then
	         fcc_convert = fac_present

	      Else If (i .Eq. 4) Then
	         fcc_monitor = fac_not_present

	      Else If (i .Eq. 6) Then
                 fcc_fields_count = 0
		 status = UPM_Comma
		 Do While (status .Eq. UPM_Comma .And. fcc_fields_count .Lt. 4)
                   fcc_fields_count = fcc_fields_count + 1
                   status = UPM_Get_Value(qualifier(i),
     .				fcc_fields(fcc_fields_count),
     .				fcc_fields_len(fcc_fields_count))
                 End Do
		 If (status .Ne. SS$_Normal) Then
		    FEP_Get_Command = status
		    Call LIB$Signal(%Val(status))
		    Return
		 End If

	      Else If (i .Eq. 7) Then
c
c		Set default start time to 00:00 of current day
c
		      status = LIB$Date_Time(date)
		      If (status .Ne. SS$_Normal) Then
		         FEP_Get_Command = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      date = date(1:12) // '00:00:00.00'
		      status = SYS$BinTim(date,time)
		      If (status .Ne. SS$_Normal) Then
		         FEP_Get_Command = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      Call CT_Binary_To_GMT(time,fcc_jstart_time)
		      fcc_jstart_select = fac_not_present

	      Else If (i .Eq. 8) Then
c
c		Set default stop time to 23:59 of current day
c
		      status = LIB$Date_Time(date)
		      If (status .Ne. SS$_Normal) Then
		         FEP_Get_Command = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      date = date(1:12) // '23:59:59.99'
		      status = SYS$BinTim(date,time)
		      If (status .Ne. SS$_Normal) Then
		         FEP_Get_Command = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      Call CT_Binary_To_GMT(time,fcc_jstop_time)
		      fcc_jstop_select = fac_not_present

	      Else If (i .Eq. 9) Then
		 fcc_smooth = fac_not_present
		 fcc_window = 1

	      Else If (i .Eq. 10) Then
		 fcc_crossgaps = fac_present

	      Else If (i .Eq. 11) Then
	         fcc_average_cal_counts = fac_present

              Else If (i .Eq. 12) Then
                 fcc_report = fac_present
                 Status = Sys$Gettim( Curtime )
                 Call CT_Binary_To_Gmt( Curtime, Curgmt )
                 Report_Default = 'FEP_' // FCC_Jstart_Time(1:7) // '_' // 
     1              FCC_Jstop_Time(1:7) // '.REP_' // Curgmt(1:9)
                 fcc_report_file = Report_Default
	      End If

	   Else If (status .Eq. UPM_Negated) Then
c
c   If something's been negated, then set the proper flag.
c
	      If (i .Eq. 1) Then
	         fcc_interactive = fac_not_present

	      Else If (i .Eq. 2) Then
                 fcc_plt_com = fac_present

	      Else If (i .Eq. 3) Then
	         fcc_convert = fac_not_present

	      Else If (i .Eq. 4) Then
	         fcc_monitor = fac_not_present

	      Else If (i .Eq. 6) Then
                 fcc_fields_count = 0

	      Else If (i .Eq. 7) Then
                 fcc_jstart_select = fac_not_present

	      Else If (i .Eq. 8) Then
                 fcc_jstop_select = fac_not_present

	      Else If (i .Eq. 9) Then
		 fcc_smooth = fac_not_present
		 fcc_window = 1

	      Else If (i .Eq. 10) Then
		 fcc_crossgaps = fac_not_present

	      Else If (i .Eq. 11) Then
		 fcc_average_cal_counts = fac_not_present

              Else If (i .Eq. 12) Then
                 fcc_report = fac_not_present

	      End If

	   Else
c
c   Assuming no defaults, set the processing flags.
c
	      If (i .Eq. 1) Then
	         fcc_interactive = fac_present

	      Else If (i .Eq. 2) Then
		 fcc_plt_com = fac_not_present

	      Else If (i .Eq. 3) Then
	         fcc_convert = fac_present

	      Else If (i .Eq. 4) Then
	         fcc_monitor = fac_not_present

	      Else If (i .Eq. 5) Then
	         fcc_plot_dev_select = fac_not_present

	      Else If (i .Eq. 6) Then
                 fcc_fields_count = 0

	      Else If (i .Eq. 7) Then
c
c		Set default start time to 00:00 of current day
c
		      status = LIB$Date_Time(date)
		      If (status .Ne. SS$_Normal) Then
		         FEP_Get_Command = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      date = date(1:12) // '00:00:00.00'
		      status = SYS$BinTim(date,time)
		      If (status .Ne. SS$_Normal) Then
		         FEP_Get_Command = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      Call CT_Binary_To_GMT(time,fcc_jstart_time)
		      fcc_jstart_select = fac_not_present

	      Else If (i .Eq. 8) Then
c
c		Set default stop time to 23:59 of current day
c
		      status = LIB$Date_Time(date)
		      If (status .Ne. SS$_Normal) Then
		         FEP_Get_Command = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      date = date(1:12) // '23:59:59.99'
		      status = SYS$BinTim(date,time)
		      If (status .Ne. SS$_Normal) Then
		         FEP_Get_Command = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      Call CT_Binary_To_GMT(time,fcc_jstop_time)
		      fcc_jstop_select = fac_not_present

	      Else If (i .Eq. 9) Then
		 fcc_smooth = fac_not_present
		 fcc_window = 1

	      Else If (i .Eq. 10) Then
		 fcc_crossgaps = fac_present

	      Else If (i .Eq. 11) Then
		 fcc_average_cal_counts = fac_present

              Else If (i .Eq. 12) Then
                 fcc_report = fac_not_present

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
     .			     jstart_default(ljstart+1:)
	fcc_jstop_time = fcc_jstop_time(1:ljstop) //
     .			     jstop_default(ljstop+1:)

	Call CT_GMT_To_Binary(fcc_jstart_time,fcc_jstart)
	Call CT_GMT_To_Binary(fcc_jstop_time,fcc_jstop)

	If (Time_LT(fcc_jstop,fcc_jstart)) Then
	   FEP_Get_Command = %Loc(FEP_InvTim)
	   Call LIB$Signal(FEP_InvTim)
	   Return
	End If

	Return
	End
