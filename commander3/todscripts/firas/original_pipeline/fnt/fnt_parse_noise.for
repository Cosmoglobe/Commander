	integer * 4 function fnt_parse_noise(start_chan, stop_chan,
	1                                    skip_chan)

c--------------------------------------------------------------------------
c
c	Function fnt_parse_noise
c		This function parses the command line
c	Author: W. K. Young
c		STI Inc.
c		July 1986
c
c	Input:
c		None
c	Output:
c		Start and stop channels
c		channel increment
c
c	Include files:
c		fnt_invoc
c		fut_params
c		upm_stat_msg
c		$ssdef
c	Calling Sequence
c		status = fnt_parse_noise()
c	Gist of routine
c--------------------------------------------------------------------------
c
c Changes:
c
c	Add /PLOTS, /WRITE. Remove /FILE_SPEC. R. Kummerer, May 6, 1987.
c
c	Remove Write to Archive; change chans to RH, RL, etc.; "CHAN_ID" to
c		"CHANNEL".    F. Shuman  1988 Mar 29.
c
c	SER 3493, Use new PLT graphics to display IFGs and spectra.
c		R. Kummerer, August 18, 1989.
c
c	Version 4.4.2 SPR 5041, 5057, R. Kummerer, Nov 15, 1989.  Inappropriate
c		FIRAS logical for fetching raw data; change CSDR$FIRAS_ARCHIVE 
c		to CSDR$FIRAS_RAW.  Correct interpretation of /INPUT.
c
c	SER 5779, Add capability to pass in a PLT command file.
c		S. Alexander, March 5, 1990.
c
c	SPR 5129,6538, Correction of INPUT default of FNT command line 
c		H. WANG, March 29, 1990.
c--------------------------------------------------------------------------

	implicit none

	include '(fnt_invoc)'
	include '(fut_params)'
	include '(upm_stat_msg)'
	include '($ssdef)'

	integer * 4	i		!a counter
	integer * 4	j		!another counter
	integer * 4	num_qual/12/	!number of qualifiers
	integer * 4	num_key(12)	!number of keywords
	integer * 4	start_chan	!starting channel
	integer * 4	stop_chan	!stopping channel
	integer * 4	skip_chan	!channel increment
	integer * 4	status		!return status from function
	integer * 2	len		!length of string from command line

	character * 32 qualifier(12)		!qualifiers
	character * 32 keyword(12,9)

	Integer 	*4	time(2)
	Integer 	*4	ljstart		!length of start time
	Integer 	*4	ljstop		!length of stop time
	Character 	*23  	date

	integer * 4	upm_present	!find qualifier-keyword presences
	integer * 4	upm_get_value	!get value from command line
	Integer 	*4	LIB$Date_Time
	Logical 	*2      Time_LT
	Integer 	*4	SYS$BinTim

	External	FNT_InvTim
c
c =================================================================
	fnt_parse_noise = fac_normal

	start_chan = 1
	stop_chan = 4
	skip_chan = 1
	ftc_preamp = fac_present
	ftc_plots = fac_present
	ftc_interactive = fac_present
	ftc_plt_com = fac_not_present
c
c	Assign qualifiers
c
	qualifier(1) = 'LOGLOG'
	qualifier(2) = 'SEMILOG'
	qualifier(3) = 'LINEAR'
	qualifier(4) = 'PREAMP'
	qualifier(5) = 'CHANNEL'
	qualifier(6) = 'WRITE'
	qualifier(7) = 'INTERACTIVE'
	qualifier(8) = 'PLOTDEVICE'
	qualifier(9) = 'JSTART'
	qualifier(10) = 'JSTOP'
	qualifier(11) = 'INPUT'
	qualifier(12) = 'PLTFILE'
c
c	Assign keywords
c
	do i=1,num_qual
	   num_key(i) = 0
	end do

	keyword(5,1) = 'CHANNEL.RH'
	keyword(5,2) = 'CHANNEL.RL'
	keyword(5,3) = 'CHANNEL.LH'
	keyword(5,4) = 'CHANNEL.LL'
	keyword(5,5) = 'CHANNEL.RIGHT'
	keyword(5,6) = 'CHANNEL.LEFT'
	keyword(5,7) = 'CHANNEL.HIGH'
	keyword(5,8) = 'CHANNEL.LOW'
	keyword(5,9) = 'CHANNEL.ALL'
	num_key(5) = 9

	keyword(11,1) = 'INPUT.RAW'
	keyword(11,2) = 'INPUT.FPP'
	keyword(11,3) = 'INPUT.FDQ'
	num_key(11) = 3

	do i=1,num_qual
	   status = upm_present(qualifier(i))
	   if(status .eq. upm_pres)then
	      if(i .eq. 1)then
	         ftc_plot = fac_loglog
	      else if(i .eq. 2)then
	         ftc_plot = fac_semilog
	      else if(i .eq. 3)then
	         ftc_plot = fac_linear
	      else if(i .eq. 4)then
	         ftc_preamp = fac_present
	      else if(i .eq. 5)then
	         do j=1,num_key(i)
	            status = upm_present(keyword(i,j))
	            if(status .eq. upm_pres)then
	               if(j .eq. 1)then
	                  start_chan = 1
	                  stop_chan = 1
	                  skip_chan = 1
	               else if(j .eq. 2)then
	                  start_chan = 2
	                  stop_chan = 2
	                  skip_chan = 1
	               else if(j .eq. 3)then
	                  start_chan = 3
	                  stop_chan = 3
	                  skip_chan = 1
	               else if(j .eq. 4)then
	                  start_chan = 4
	                  stop_chan = 4
	                  skip_chan = 1
	               else if(j .eq. 5)then
	                  start_chan = 1
	                  stop_chan = 2
	                  skip_chan = 1
	               else if(j .eq. 6)then
	                  start_chan = 3
	                  stop_chan = 4
	                  skip_chan = 1
	               else if(j .eq. 7)then
	                  start_chan = 1
	                  stop_chan = 3
	                  skip_chan = 2
	               else if(j .eq. 8)then
	                  start_chan = 2
	                  stop_chan = 4
	                  skip_chan = 2
	               else if(j .eq. 9)then
	                  start_chan = 1
	                  stop_chan = 4
	                  skip_chan = 1
	               end if
	            end if
	         end do
	      else if(i .eq. 6)then
	         ftc_write = fac_present
	      else if(i .eq. 7)then
		 ftc_interactive = fac_present
		 ftc_plots = fac_present
		 ftc_plots_device = 'PLT_DEVICE'
	      else if(i .eq. 8)then
	         ftc_plots = fac_present
	         status = upm_get_value(qualifier(i),ftc_plots_device,len)
	      else if(i .eq. 9)then

		 ftc_jstart_select = fac_present
	         status = UPM_Get_Value(qualifier(i),ftc_jstart_time,
     .					len)
	         If (status .Eq. UPM_Absent) Then
c
c		Assign the start time to be the start of the current day
c
		      status = LIB$Date_Time(date)
		      If (status .Ne. SS$_Normal) Then
		         FNT_Parse_Noise = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      date = date(1:12) // '00:00:00.00'
		      status = SYS$BinTim(date,time)
		      If (status .Ne. SS$_Normal) Then
		         FNT_Parse_Noise = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      Call CT_Binary_To_GMT(time,ftc_jstart_time)
		      ftc_jstart_select = fac_not_present

		 Else

		    If (status .Ne. SS$_Normal) Then
	   	       FNT_Parse_Noise = status
	   	       Call LIB$Signal(%Val(status))
	   	       Return
		    End If

	         End If

	      Else If (i .Eq. 10) Then

		 ftc_jstop_select = fac_present
	         status = UPM_Get_Value(qualifier(i),ftc_jstop_time,
     .					len)
	         If (status .Eq. UPM_Absent) Then
c
c		Assign the stop time to be the end of the current day
c
		      status = LIB$Date_Time(date)
		      If (status .Ne. SS$_Normal) Then
		         FNT_Parse_Noise = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      date = date(1:12) // '23:59:59.99'
		      status = SYS$BinTim(date,time)
		      If (status .Ne. SS$_Normal) Then
		         FNT_Parse_Noise = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      Call CT_Binary_To_GMT(time,ftc_jstop_time)
		      ftc_jstop_select = fac_not_present

		 Else

		    If (status .Ne. SS$_Normal) Then
	   	       FNT_Parse_Noise = status
	   	       Call LIB$Signal(%Val(status))
	   	       Return
		    End If

	         End If

	      Else If (i .Eq. 11) Then
	         do j=1,num_key(i)
	            status = upm_present(keyword(i,j))
	            if(status .eq. upm_pres)then
	               if(j .eq. 1)then
			  ftc_input = 'RAW'
	               else if(j .eq. 2)then
			  ftc_input = 'FPP'
	               else if(j .eq. 3)then
			  ftc_input = 'FDQ'
	               end if
	            end if
	         end do
	      else if(i .eq. 12)then
	         ftc_plt_com = fac_present
	         status = upm_get_value(qualifier(i),ftc_plt_com_file,len)
	      end if

	   else if(status .eq. upm_defaulted)then

	         if(i .eq. 4)then
	            ftc_preamp = fac_present
	         else if(i .eq. 6)then
		    ftc_write = fac_not_present
	         else if(i .eq. 7)then
		    ftc_interactive = fac_present
		    ftc_plots = fac_present
		    ftc_plots_device = 'PLT_DEVICE'
	         else if(i .eq. 8)then
		    ftc_plots = fac_present
		    ftc_plots_device = 'PLT_DEVICE'
	         Else If (i .Eq. 9) Then
c
c		Set default start time to 00:00 of current day
c
		      status = LIB$Date_Time(date)
		      If (status .Ne. SS$_Normal) Then
		         FNT_Parse_Noise = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      date = date(1:12) // '00:00:00.00'
		      status = SYS$BinTim(date,time)
		      If (status .Ne. SS$_Normal) Then
		         FNT_Parse_Noise = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      Call CT_Binary_To_GMT(time,ftc_jstart_time)
		      ftc_jstart_select = fac_not_present

	          Else If (i .Eq. 10) Then
c
c		Set default stop time to 23:59 of current day
c
		      status = LIB$Date_Time(date)
		      If (status .Ne. SS$_Normal) Then
		         FNT_Parse_Noise = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      date = date(1:12) // '23:59:59.99'
		      status = SYS$BinTim(date,time)
		      If (status .Ne. SS$_Normal) Then
		         FNT_Parse_Noise = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      Call CT_Binary_To_GMT(time,ftc_jstop_time)
		      ftc_jstop_select = fac_not_present

	         Else If (i .Eq. 11) Then
		    ftc_input = 'FDQ'
	         end if

	   else if(status .eq. upm_negated)then

	         if(i .eq. 4)then
	            ftc_preamp = fac_not_present
	         else if(i .eq. 6)then
		    ftc_write = fac_not_present
	         else if(i .eq. 7)then
		    ftc_interactive = fac_not_present
		    ftc_plots = fac_present
		    ftc_plots_device = 'PLT_HARDCOPY'
	         else if(i .eq. 8)then
		    ftc_plots = fac_not_present
	         Else If (i .Eq. 9) Then
                    ftc_jstart_select = fac_not_present
	         Else If (i .Eq. 10) Then
                    ftc_jstop_select = fac_not_present
	         end if

           else

	           If (i .eq. 7) Then
		      ftc_interactive = fac_present
		      ftc_plots = fac_present
		      ftc_plots_device = 'PLT_DEVICE'
	           Else If (i .eq. 8) Then
		      ftc_plots = fac_present
		      ftc_plots_device = 'PLT_DEVICE'
	           Else If (i .Eq. 9) Then
c
c		Set default start time to 00:00 of current day
c
		      status = LIB$Date_Time(date)
		      If (status .Ne. SS$_Normal) Then
		         FNT_Parse_Noise = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      date = date(1:12) // '00:00:00.00'
		      status = SYS$BinTim(date,time)
		      If (status .Ne. SS$_Normal) Then
		         FNT_Parse_Noise = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      Call CT_Binary_To_GMT(time,ftc_jstart_time)
		      ftc_jstart_select = fac_not_present

	         Else If (i .Eq. 10) Then
c
c		Set default stop time to 23:59 of current day
c
		      status = LIB$Date_Time(date)
		      If (status .Ne. SS$_Normal) Then
		         FNT_Parse_Noise = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      date = date(1:12) // '23:59:59.99'
		      status = SYS$BinTim(date,time)
		      If (status .Ne. SS$_Normal) Then
		         FNT_Parse_Noise = status
		         Call LIB$Signal(%Val(status))
		         Return
	              End If
		      Call CT_Binary_To_GMT(time,ftc_jstop_time)
		      ftc_jstop_select = fac_not_present

	         Else If (i .Eq. 11) Then
		      ftc_input = 'RAW'
	         end if

	   end if
	end do

c
c Setup the timerange variables.
c
	ljstart = index(ftc_jstart_time,' ')-1
	If (ljstart .Eq. -1)ljstart=14
	ljstop = index(ftc_jstop_time,' ')-1
	If (ljstop .Eq. -1) ljstop=14

	ftc_jstart_time = ftc_jstart_time(1:ljstart) //
     .			     fac_jstart_default(ljstart+1:)
	ftc_jstop_time = ftc_jstop_time(1:ljstop) //
     .			     fac_jstop_default(ljstop+1:)

	Call CT_GMT_To_Binary(ftc_jstart_time,ftc_jstart)
	Call CT_GMT_To_Binary(ftc_jstop_time,ftc_jstop)

	If (Time_LT(ftc_jstop,ftc_jstart)) Then
	   FNT_Parse_Noise = %Loc(FNT_InvTim)
	   Call LIB$Signal(FNT_InvTim)
	   Return
	End If

	return
	end
