C-------------------------------------------------------------------------------

	Integer*4 Function FXT_Trigger_Plot ( Eng_Time, Cur_Table, 
	1	  Min_Offtime, Max_Offtime, Plot, Report, Lun_Rpt, 
	2	  Plot_Table, Plot_Start, Plot_Stop)

C-------------------------------------------------------------------------------
C
C	Purpose: This routine accumulates logical flags for plotting 
C       engineering words which exceeded the current minima or maxima. 
C       When the time of the engineering record is Max_Offtime hours 
C       beyond the first time an extrema is exceeded, the command text 
C       file to produce Engplots is written. The tables are then reset.
C
C
C	Author: Shirley M. Read
C		STX, November, 1988
C
C	Invocation: Status = FXT_Trigger_Plot( Eng_Time, Cur_Table, 
C	            Min_Offtime, Max_Offtime, Plot, Report, Lun_Rpt, 
C	            Plot_Table, Plot_Start, Plot_Stop )
C
CH	Change Log:
CH
C	  ----------------------------------------------------------------------
C
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name		   Type	   Description
C 	  ----------------------------------------------------------------------
C	  Eng_Time(2)      I*4     Time of current engineering record.   
C         Cur_Table(118)   L*1     Flags for extremea exceeded on current record
C         Min_Offtime      I*4     Minimum offset in hours from time when first 
C	                           extrema was exceeded.
C         Max_Offtime      I*4     Maximum offset in hours from time when first 
C	                           extrema was exceeded.
C	  Plot             I*4     Flag for plot type: 0=none,1=line,2=laser
C	  Report           L*1     Flag to write report.
C	  Lun_Rpt          I*4     logical unit for writing report.
C
C	Output Parameters:
C	  Name		   Type	   Description
C 	  ----------------------------------------------------------------------
C	  Plot_Table(118)  L*1     Table of flags for fields exceeding limits.
C			           To be used to trigger Engplots.
C	  Plot_Start(2)    I*4     Plot start time.
C	  Plot_Stop(2)     I*4     Plot stop time.
C	
C	Subroutines Called:
C         FXT_Write_Engplot
C	  AUT_Dfloat2Adt
C	  Time_Gt
C	  Lib$Signal
C         Lib$Subx
C         Lib$Addx
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C	  FXT_Msg.Txt
C
C	Processing Method:
C
C	Set the function status to FXT_Normal.
C	If ( First time in routine ) Initialize buffers, counters and time. 
C	If ( New_Plot ) Then
C	   Record the time of the first exceeded limit.
C	   Set the start time for the plots as Min_Offtime hours earlier. 
C	   Set the final time as Max_Offtime hours later.
C	Else if ( the new record time exceeds the final time )
C	   Write the text file to produce Engplots.
C	   Reset all buffers and counters.
C	   Reset the times for the next plot.
C	Endif
C	Do for all the engineering analogs
C	  If the flag is set in the current exceed table, 
C           set the corresponding flag in the plot table.
C	Enddo
C
C****************************************************************************

	Implicit None

!	External Parameters and Include Files
	
	Include '(FXT_Msg)'

!	Input/Output Passed Parameters

	Integer*4 Eng_Time(2)                ! Engineering record time
	Logical*1 Cur_Table(118)             ! Flags for extrema on current rec
	Integer*4 Min_Offtime                ! Minimum hours offset in plot time
	Integer*4 Max_Offtime                ! Maximum hours offset in plot time
	Integer*4 Plot                       ! Flag for plot type
	Logical*1 Report 		     ! Flag for writing report
	Integer*4 Lun_Rpt                    ! Logical unit for writing report 
	Logical*1 Plot_Table(118)            ! Flags for plots
	Integer*4 Plot_Start(2)		     ! Plot start time
	Integer*4 Plot_Stop(2)		     ! Plot stop time
!	Functions

	Integer*4	FXT_Write_Engplot    ! Write command line to file
	Logical*1       Time_Gt		     ! Time greater than

!	Local Variables
	
	Integer*4	Status		! Dummy status variable

	Integer*4       B_Min_Hr(2)       ! Binary quad min hour delta time
	Integer*4       B_Max_Hr(2)       ! Binary quad max hour delta time
	Real*8          R_Min_Hr          ! Real min hour delta time
	Real*8          R_Max_hr          ! Real max hour delta time

	Logical*1 	First /.true./	 ! First time
	Logical*1       New_Plot         ! Start new plot flag
	Integer*2	Ix

!	Set function return status to normal.

	FXT_Trigger_Plot = %loc(FXT_Normal)

!       First Time: Initialize
 
	If ( First ) Then
	    New_Plot = .True.
	    Do Ix = 1, 118
	      Plot_Table(Ix) = .False.
	    Enddo
	    R_Min_Hr = Float(Min_Offtime) * 3600.0  * 10000000.0
	    R_Max_Hr = Float(Max_Offtime) * 3600.0  * 10000000.0
	    Call AUT_Dfloat2Adt ( R_Min_Hr, B_Min_Hr )
	    Call AUT_Dfloat2Adt ( R_Max_Hr, B_max_Hr )
	    First = .False.
	Endif		! First time
!
!	If a new plot is starting, save the time. Set the plot start and stop.
!       Else if the plot stop time is exceeded, write the command line for 
!       Engplots. Reset flags and plot time.

	If ( New_Plot ) Then
	   Call Lib$Subx(Eng_Time,B_Min_Hr,Plot_Start)
	   Call Lib$Addx(Eng_Time,B_Max_Hr,Plot_Stop)
	   New_Plot = .False.
 	Elseif ( Time_Gt(Eng_Time,Plot_Stop) ) Then
	   Status = FXT_Write_Engplot(Plot_Table,Plot_Start,Plot_Stop,
	1	    Plot, Report, Lun_Rpt)
	   If ( Status .ne. %loc(FXT_Normal)) Then
	     FXT_Trigger_Plot = %loc(FXT_Aberr)
	   Else
	     Call Lib$Subx(Eng_Time,B_Min_Hr,Plot_Start)
	     Call Lib$Addx(Eng_Time,B_Max_Hr,Plot_Stop)
	     New_Plot = .False.
	   Endif 	! Status from FXT_Write_Engplot
	Endif		! New plot to start or existing plot to be written

!	Set the flags in the plot table according to the current exceed table.

	If ( FXT_Trigger_Plot .Eq. %loc(FXT_Normal) ) Then

	  Do Ix = 1, 118

	    If ( Cur_Table(Ix) ) Plot_Table(Ix) = .True.

	  Enddo

	Endif		! Function status is normal

	Return
	End

