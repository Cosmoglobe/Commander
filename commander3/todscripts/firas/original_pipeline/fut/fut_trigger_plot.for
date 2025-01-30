	INTEGER*4 FUNCTION FUT_TRIGGER_PLOT ( SCI_REC, CHANNEL, PLOT, 
	1	  STATMON_FLAG, VOLTS_FLAG, CUR_FLAG, MIN_OFFTIME, 
	2	  MAX_OFFTIME, PLOT_DEVICE, PLOT_TABLE,
	3	  PLOT_START, PLOT_STOP, ENGLIM )
C/
C/ PROGRAM NAME: 
C/	FUT_TRIGGER_PLOT
C/
C/ PROGRAM DESCRIPTION:
C/	This routine accumulates logical flags for plotting engineering words
C/      which exceeded read or yellow limits. Whether yellow and red or just
C/	red limits depends on the value of the plot flag. When the average 
C/	science time is max hours beyond the first time a limit is exceeded,
C/	the text file to produce Engplots is written. The tables are then reset.
C/
C/ AUTHOR:
C/	Shirley M. Read ( STX )  May 1988.
C/
CH CHANGE LOG:
CH
CH	Version 4.2.1 02/09/89, SPR 2550, Shirley M. Read, STX
CH		Provide user option of setting the time interval for triggered
CH		plots.
CH	Version 4.2.1 02/09/89, SPR 2990, Shirley M. Read, STX
CH		FDQ needs a command line option for plot device. The batch mode
CH	        runs should be able to go to the lineprinter or laserprinter.
CH	Version 4.2.1 02/25/89, SPR 3323, Shirley M. Read, STX
CH	        Index for GRT bit flags goes out of bounds.
CH
C/
C/ CALLING SEQUENCE:
C/	Status = FUT_Trigger_Plot( Sci_Rec, Channel, Plot, Statmon_Flag,
C/		 Volts_Flag, Cur_Flag, Min_Offtime, Max_Offtime, Plot_Device,
C/		 Plot_Table, Plot_Start, Plot_Stop, EngLim )
C/
C/ INPUT PARAMETERS:
C/	Name             Type    Description
C/	-----------------------------------------------------------------
C/	Sci_Rec          Record  Science record with quality flags.
C/	Channel          I*2     Current channel being processed.
C/	Plot             I*2     Plot if red only or yellow and red exceeded .
C/	Statmon_Flag(2)  L*1     Flag for status monitor limits exceeded, A/B.
C/	Volts_Flag(2,10) L*1     Flag for IPDU voltage limits exceeded.
C/	Cur_Flag(2,6)    L*1     Flag for IPDU current limits exceeded.
C/	Min_Offtime      I*4     Hours to offset plot start
C/	Max_Offtime      I*4     Hours to offset plot stop
C/	Plot_Device      I*4     Plot device: line or laser printer
C/
C/ OUTPUT PARAMETERS:
C/	Name             Type    Description
C/	-----------------------------------------------------------------
C/	Plot_Table(118)  L*1     Table of flags for fields exceeding limits.
C/			         To be used to trigger Engplots.
C/	Plot_Start(2)    I*4     Plot start time.
C/	Plot_Stop(2)     I*4     Plot stop time.
C/	EngLim		 I*4	 Flag indicating engineering limits violations.
C/
C/ ENTRY/EXIT:
C/	Normal function entry and return.
C/
C/ ERROR HANDLING
C/	Calls to Lib$Signal and interface with Fut_Error condition handler.
C/	
C/ INPUT/OUTPUT DISK FILES: 
C/	None
C/
C/ PROCS/MACROS/INCLUDE FILES:
C/	FUT_QUALFLAGS
C/
C/ SUBROUTINES CALLED:
C/	FUT_Write_Engplot
C/	Btest
C/	Time_Gt
C/
C/ SHARED DATA AREAS:
C/
C/ METHOD USED:
C/	If ( First time in routine ) Initialize buffers, counters and time. 
C/	If ( New_Plot and at least one limit was exceeded ) 
C/	   Record the time of the first exceeded limit.
C/	   Set the start time for the plots  as min_offtime hours earlier. 
C/	   Set the final time as max_offtime hours later.
C/	Else if ( Trigger and the new record time exceeds the final time )
C/	   Write the text file to produce Engplots.
C/	   Reset all buffers and counters.
C/	Endif
C/	Do for all data quality flags on science record relating to 
C/	   engineering fields which have red and yellow limits
C/	   If plot type is red limits only, check all flags for exceeding red
C/	      limits and set the corresponding flag in the plot table
C/	      (Keep track of all channels exceeded for channel related flags
C/	      and both sides for FIRAS A/B side related flags.)
C/	   Elseif plot type is yellow and red, check flags for exceeding yellow
C/	      limits and set the corresponding flag in the plot table
C/	      (Keep track of all channels exceeded for channel related flags
C/	      and both sides for FIRAS A/B side related flags.)
C/	   Endif
C?	Enddo
C/
C/ SPECIAL HANDLING:
C/	None
C/
C**************************************************************************************

	IMPLICIT NONE

!	External Parameters and Include Files
	
	External   	FUT_NORMAL
	External	FUT_ABERR
	External	FUT_WRTPLOTERR
	External	FUT_INVALPLOT
	External        FUT_TRIGPLOTERR
	Include         '(FUT_QUALFLAGS)'

!	Input/Output Passed Parameters

	Dictionary 	'NFS_SDF'
	Record	/ NFS_SDF / Sci_Rec  ! Science record

	Integer*2	Channel	  ! Micro-processor channel
	Integer*2       Plot	  ! Plot if red only or yellow and red exceeded
	Logical*1       Statmon_Flag(2) !Status monitor limit exceeded, A/B side
	Logical*1       Volts_Flag(2,10) !IPDU voltage limit exceeded, A/B by 10
	LogicaL*1       Cur_Flag(2,6)   !IPDU current limit exceeded, A/B by 6
	Integer*4       Min_Offtime     ! Hours to offset plot start
	Integer*4       Max_Offtime     ! Hours to offset plot stop
	Integer*4       Plot_Device     ! Plot device: line or laser printer
	Integer*4	EngLim

!	Functions

	Integer*4	FUT_Write_Engplot    ! Write command line to file
	Logical*1       Time_Gt		     ! time greater than

!	Local Variables
	
	Integer*4 	Retstat		! Return status
	Integer*4 	Success / 1 /, ERROR / 2 /  ! Values for status
	Integer*4	Status		! Dummy status variable

	Integer*4       B_Min_Hr(2)    ! Binary quad min hour delta time
	Integer*4       B_Max_Hr(2)    ! Binary quad max hour delta time
	Real*8          R_Min_Hr       ! Real min hour delta time
	Real*8          R_Max_hr       ! Real max hour delta time

	Integer*4       Plot_Start(2) ! Binary quad start time for plot
	Integer*4       Plot_Stop(2)  ! Binary quad stop time for plot
	Integer*4       First_Badlim(2) ! Binary time for first limit exceeded

	Logical*1       Plot_Table(118)  ! Flags for plotting engineering 
				         ! fields
	Integer*2	Bit_Flags(8)     ! I*2 Storage for byte size bit masks
	Byte		Bit_Ext(2,8)     
	Equivalence     (Bit_Flags(1), Bit_Ext(1,1))
	Integer*2	Count		 ! Counter
	Byte		Pval(2) / 1, 0 / ! GT Values of qualflags for plot trig
	Integer*2       Chan_Tmp(4) / 67, 68, 69, 70 /  !Array pos. channel temp
	Integer*2	Bol_Vol(4)  / 81, 82, 83, 84 /  !Array pos. chan bol vol
	Integer*2	Flag		 ! Flag position in DQ array
	Logical*1 	First /.true./	 ! First time
	Logical*1       New_Plot         ! Start new plot flag
	Logical*1       Trigger          ! Trigger plot flag 
	Integer*2 	RED / 1 /, YELLOW / 2 /  ! Red and yellow flags for plot
	Byte            Zero	        ! Zero value
	Parameter       ( Zero = 0 )
	Byte		One
	Parameter       ( One = 1 )	! Good quality flag or valid Plot
	Byte		Two
	Parameter       ( Two = 2 )	! Valid Plot

	Integer*2	Ix, Jx

!	Set return status to success.

	Retstat = Success

!       First Time: Initialize
 
	If ( First ) then
	    Trigger = .false.
	    New_Plot = .true.
	    DO Jx = 1, 110
	      Plot_Table(Jx) = .false.
	    ENDDO
	    R_Min_Hr = Float(Min_Offtime) * 3600.0  * 10000000.0
	    R_Max_Hr = FLoat(Max_Offtime) * 3600.0  * 10000000.0
	    Call AUT_Dfloat2Adt ( R_Min_Hr, B_Min_Hr )
	    Call AUT_Dfloat2Adt ( R_Max_Hr, B_Max_Hr )
	    first = .false.
	Endif		! First time
!
!	If a new plot is starting and a limits are exceeded on the new record,
!       save the time. Set the plot start and stop. Else if a plot is triggered
!	and the plot stop time is exceeded, write the command line for Engplots.
!	Check the new record for limits exceeded. If so, start a new plot.

	If ((New_Plot) .and. (Sci_Rec.DQ_Data.Data_Quality(110) .gt. One)) then
	   Do ix = 1, 2
	     First_Badlim(ix) = Sci_Rec.Ct_Head.Time(ix)
	   Enddo
	   Call Lib$Subx(First_Badlim,B_Min_Hr,Plot_Start)
	   Call Lib$Addx(First_Badlim,B_Max_Hr,Plot_Stop)
	   New_Plot = .false.
  	   Trigger = .true.
	Elseif ((Trigger) .and. (Time_Gt(Sci_Rec.Ct_Head.Time,Plot_Stop))) Then
	   Status = FUT_Write_Engplot(Plot_Table,Plot_Start,Plot_Stop,
	1	    EngLim, Plot_Device)
	   If ( Status .ne. %loc(FUT_NORMAL)) Then
	     If ( Status .eq. %loc(FUT_ABERR)) Status = Zero
	     Retstat = error
	     Call Lib$Signal(FUT_TRIGPLOTERR, %val(1), %val(status))
	   Else
	     If (Sci_Rec.DQ_Data.Data_Quality(110) .gt. One) Then
	        Do ix = 1, 2
	          First_Badlim(ix) = Sci_Rec.Ct_Head.Time(ix)
	        Enddo
	        Call Lib$Subx(First_Badlim,B_Min_Hr,Plot_Start)
	        Call Lib$Addx(First_Badlim,B_Max_Hr,Plot_Stop)
	        New_Plot = .false.
  	        Trigger = .true.
	     Else
		New_Plot = .true.
	        Trigger = .false.
	     Endif
	   Endif 	! Status from FUT_Write_Engplot
	Endif		! New plot to start or existing plot to be written

!	Check the value of Plot and if not valid, signal an error.

	If (( Plot .ne. One ) .and. ( Plot .ne. Two )) Then
	  Retstat = Error
	  Call Lib$Signal(FUT_INVALPLOT, %val(1), %val(Plot))
	Endif

	If ( Retstat .eq. Success ) Then

!	Process the GRTs. Each data quality flag is a bit mask.



	  flag = Flg_Grtred_St - 1

	  Do Ix = 1, 8
	    flag = flag + 1
	    Call Lib$Movc3(1,Sci_Rec.Dq_Data.Data_Quality(flag),
	1	Bit_Ext(1,Ix))
	  Enddo
	  Count = 0
	  Do Ix = 1, 8
	    Do Jx = 0 , 7
	      Count = Count + 1
	      If (BTest(Bit_Flags(Ix),Jx)) Plot_Table(Count) = .true.
	    Enddo
	  Enddo

!	If plot = yellow, check yellow limits flags also.

	  If( Plot .eq. YELLOW ) Then

	    flag = Flg_Grtyel_St - 1

	    Do Ix = 1, 8
	      flag = flag + 1
	      Call Lib$Movc3(1,Sci_Rec.Dq_Data.Data_Quality(flag),
	1	Bit_Ext(1,Ix))
	    Enddo
	    Count = 0
	    Do Ix = 1, 8
	      Do Jx = 0 , 7
	        Count = Count + 1
	        If (BTest(Bit_Flags(Ix),Jx)) Plot_Table(Count) = .true.
	      Enddo
	    Enddo
	  Endif

!	Process IPDU temperatures.

	  Flag = Flg_IPDU_Tmp - 1
	  Do Ix = 65,66
	    Flag = Flag + 1
	    If (Sci_Rec.DQ_Data.Data_Quality(Flag) .gt. Pval(Plot))
	1	Plot_Table(Ix) = .true.
	  Enddo

!	Process Channel Temperature.

	  If (Sci_Rec.DQ_Data.Data_Quality(Flg_Chan_Tmp_St)
	1  .gt. Pval(Plot)) Plot_Table(Chan_Tmp(Channel)) = .true.

!	Process Drive Box Temperatures.

	  Flag = Flg_DBox_Tmp - 1
	  Do Ix = 71, 72
	    Flag = Flag + 1
	    If (Sci_Rec.DQ_Data.Data_Quality(Flag) .gt. Pval(Plot))
	1	Plot_Table(Ix) = .true.
	  Enddo

!	Process Status Monitor Temperatures. One flag can be set by either side.

	  If(Sci_Rec.DQ_Data.Data_Quality(Flg_Stmon_Tmp) .gt. Pval(Plot)) Then
	    Count = 72
	    Do Ix = 1, 2
	      Count = Count + 1
	      If (Statmon_Flag(Ix)) Plot_Table(Count) = .true.
	    Enddo
	  Endif

!	Process the remainder of assorted engineering analog temps and currents

	  Flag = Flg_Chpa_Tmp - 1
	  Do Ix = 75, 80
	    Flag = Flag + 1
	    If (Sci_Rec.DQ_Data.Data_Quality(Flag) .gt. Pval(Plot))
	1	Plot_Table(Ix) = .true.
	  Enddo

!	Process the Bolometer Bias Voltages.

	  If (Sci_Rec.DQ_Data.Data_Quality(Flg_Bol_Vol_St)
	1  .gt. Pval(Plot)) Plot_Table(Bol_Vol(Channel)) = .true.

!	Process the IPDU Voltages. Only 10 flags, set by either A or B side.

	  Flag = Flg_IPDU_Vol_St - 1
	  Count = 84
	  Do Ix = 1, 10
	    Flag = Flag + 1
	    Do Jx = 1, 2
	      Count = Count + 1
	      If (Sci_Rec.DQ_Data.Data_Quality(Flag) .gt. Pval(Plot)) Then
		 If (Volts_Flag(Jx,Ix)) Plot_Table(Count) = .true.
	      Endif
	    Enddo
	  Enddo

!	Process the IPDU currents. Only 6 flags, set by either A or B side.

	  Flag = Flg_IPDU_Cur_St - 1
	  Count = 104
	  Do Ix = 1, 6
	    Flag = Flag + 1
	    Do Jx = 1, 2
	      Count = Count + 1
	      If (Sci_Rec.DQ_Data.Data_Quality(Flag) .gt. Pval(Plot)) Then
		 If (Cur_Flag(Jx,Ix)) Plot_Table(Count) = .true.
	      Endif
	    Enddo
	  Enddo

!	Process the LMACs.

	  Flag = Flg_AcLmac_Tmp - 1
	  Do Ix = 117, 118
	    Flag = Flag + 1
	    If (Sci_Rec.DQ_Data.Data_Quality(Flag) .gt. Pval(Plot))
	1	Plot_Table(Ix) = .true.
	  Enddo

	Endif		! Retstat is success

!	Set function to return status

	If ( Retstat .eq. Success ) then
	  FUT_Trigger_Plot = %loc(FUT_NORMAL)
	Else
	  FUT_Trigger_Plot = %loc(FUT_ABERR)
	Endif

	Return
	End

