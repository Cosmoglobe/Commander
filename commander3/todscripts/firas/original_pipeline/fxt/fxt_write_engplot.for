C-------------------------------------------------------------------------------

	Integer*4 Function FXT_Write_Engplot ( Plot_Table, 
	1	  Plot_Start, Plot_Stop, Plot, Report, Lun_Rpt)

C-------------------------------------------------------------------------------
C
C	Purpose: This routine accumulates build command lines for batch mode
C                runs to produce Engplots of engineering fields that 
C                exceeded extrema.
C
C	Author: Shirley M. Read
C		STX, November, 1988
C
C	Invocation: Status = FXT_Write_Engplot ( Plot_Table, Plot_Start,
C                   Plot_Stop, Plot, Report, Lun_Rpt)
C
CH	Change Log:
CH
CH		Version 4.2.3 05/08/89, SPR 3806, Shirley M. Read, STX
CH			Engplots has been modified to plot up to a maximum
CH			of 4 fields per plot so that the plot labels will be
CH			readable. FXT generates an Engplots command file. The
CH			maximum number of fields has been changed to 4.
CH
CH              SPR 6727, FXT must support talaris 1590t printer
CH               Version 6.1, H. Wang, STX, May 15, 1990
C	  ----------------------------------------------------------------------
C
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name		   Type    Description
C 	  ----------------------------------------------------------------------
C	  Plot_Table(118)  L*1     Plot table for Engplots.
C	  Plot_Start(2)    I*4     Plot start time.
C	  Plot_Stop(2)     I*4     Plot stop time.
C	  Plot             I*4     Flag for type of plot: 0=none,1=line,2=laser
C	  Report           L*1     Flag to write report.
C	  Lun_Rpt          I*4     Logical unit for report file.
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	
C	Subroutines Called:
C
C     	  Lib$Get_Lun
C         Lib$Signal
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C	  FXT_Msg.Txt
C	  FUT_Plot_Names
C         $SSDef
C
C	Processing Method:
C
C       Set the function status to FXT_Normal.
C	If first time, get a logical unit and open the report file. 
C	Convert the plot start and stop times to character and insert in 
C	   command line.
C	Do for all engineering fields which exceeded the previous extrema.
C	   Insert the fields to be plotted from the FUT_Plot_Names into the
C	   command line. When the number reaches 4 or the total characters
C	   equals 254, write this command line. Start the next command line
C	   with the new field.	
C	Enddo
C	If there is an error, write an error message and set the function 
C	   status to FXT_Aberr.
C	If a report is requested, write the plot information, including the 
C          start and stop plot times, in the report.
C	Return
C
C**************************************************************************************

	Implicit None

!	Input/Output Passed Parameters.

	Logical*1       Plot_Table(118)	  ! Plot flags for Engplots.

	Integer*4       Plot_Start(2) ! Binary quad start time for plot
	Integer*4       Plot_Stop(2)  ! Binary quad stop time for plot
	Integer*4       Plot          ! Flag for type of plot: 0 = none,
				      ! 1 = lineprinter, 2 = laserprinter
	Logical*1       Report        ! Flag to write report.
	Integer*4       Lun_Rpt       ! Logical unit for report file.

!	Functions.
	
	Integer*4	Lib$Get_Lun	! Get logical unit number

!	Include Files.

	Include         '(FXT_Msg)'
	Include	        '(FUT_Plot_Names)'
	Include		'($SSDef)'

!	Local Variables
        character*42    plot_file

	Character*254	Command_Line
	Data		Command_Line(1:40) /
	1    '$ENG/NOIN/CONV/NORE/JSTA=yydddhhmmssmmm/' /
	Data		Command_Line(41:66) /
	1    'JSTO=yydddhhmmssmmm/FIEL=('  /
	Character*11    Plot_Dev / '/PLO="/QMS"' /

	Character*230   Blank		! Blank char string
	Character*1	Blanks(188) 	/ 188 * ' ' /
	Equivalence     ( Blank, Blanks(1) )

	Character*15	Filename / 'FXT_Engplot.Com' /

	Integer*4	Zero    / 0 /   ! Zero value
	Integer*4	Status		! Dummy status variable

	Integer*4 	Lun		! Logical unit number to write text file
	Logical*1	First_Time /.True./	! First time in routine
	Character*14    Start, Stop     ! Plot start and stop times
	Integer*2	Max_Fields	! Max fields for 1 Engplot
	Parameter	( Max_Fields = 4 )
	Integer*2	Max		! Max characters on command line
	Parameter	( Max = 242 )   ! before the plot device
	Integer*2	Table_Max	! Max plot table element
	Parameter	( Table_Max = 118 )
	Integer*2	Fcount		! Counter for fields in 1 plot
	Integer*2	Bcount		! Byte or char counter for line
	Integer*2	Totcount	! Total line count
	Integer*2	Ix, Jx		! Indices

!	Set the function status to normal.

	FXT_Write_Engplot = %loc(FXT_Normal)

!	Initialization.

	If ( First_Time ) Then

	    Status = Lib$Get_Lun ( Lun )

 	    If ( Status .ne. SS$_Normal ) Then
		FXT_Write_Engplot= %loc(FXT_Aberr)
		Call Lib$Signal ( FXT_GetLunErr, %val(1), %val(Status))
	    Endif
	    If (FXT_Write_Engplot .Eq. %loc(FXT_Normal) ) Then

		Open ( Unit=Lun, File=Filename, Status= 'NEW', 
	1	  Recl=254, Carriagecontrol='LIST', Iostat=Status )

	        If ( Status .ne. Zero ) Then
		  FXT_Write_Engplot = %loc(FXT_Aberr)
		  Call Lib$Signal ( FXT_RMSOpenPlt, %val(1), %val(Status))
	        Endif
	    Endif

	    If ( Plot .Eq. 1 ) Then		! 1 = Lineprinter
	      Plot_Dev = '/PLO="/PRI"' 
	    Elseif ( Plot .Eq. 2 ) Then		! 2 = Laserprinter
	      Plot_Dev = '/PLO="/QMS"' 
	    Elseif ( Plot .Eq. 3 ) Then		! 2 = Laserprinter
	      Plot_Dev = '/PLO="/TF "' 
	    Endif
	    First_Time = .False.
	Endif

	Command_Line(67:254) = Blank(1:188) 
	Bcount = 66
	Fcount = 0

!	Convert plot times to ASCII characters and insert times on command line.

	Call CT_Binary_to_GMT (Plot_Start, Start)
	Call CT_Binary_to_GMT (Plot_Stop, Stop)

	Command_Line(26:39) = Start
	Command_Line(46:59) = Stop

!	Scan the Plot Table for fields to be plotted and insert the 
!	field names from FUT_Plot_Names.

	Do Ix = 1, Table_Max
	
	   If (( Plot_Table(Ix) ) .and. 
	1	( FXT_Write_Engplot .eq. %loc(FXT_Normal))) Then 
		Fcount = Fcount + 1
	        Totcount = Bcount + Name_Len(Ix)
	        If (( Totcount .gt. Max ) .or.
	1	   ( Fcount .gt. Max_Fields)) Then
		   Command_Line(Bcount:Bcount) = ')'
	           Command_Line(Bcount+1:Bcount+11) = Plot_Dev(1:11)
		   Write(Unit=Lun, Fmt=100, Iostat= Status) Command_Line
	           If ( Status .ne. Zero ) Then
		      FXT_Write_Engplot = %loc(FXT_Aberr)
		      Call Lib$Signal(FXT_RMSWritePlt, %val(1), %val(Status))
	           Endif
		   Command_Line(67:254) = Blank(1:188)
	           Bcount = 66
		   Fcount = 1
	        Endif

	        If ( FXT_Write_Engplot .Eq. %loc(FXT_Normal)) Then
		   Command_Line(Bcount+1:Bcount+Name_Len(Ix)) 
	1		= Plot_Name(Ix)
		   Bcount = Bcount + Name_Len(Ix) + 1
	           Command_Line(Bcount:Bcount) = ','
	        Endif
	   Endif	! Plot_Table(Ix) and function status is normal

	   If (( FXT_Write_Engplot .Eq. %loc(FXT_Normal) ) .and.
	1	(Ix .eq. Table_Max) .and. ( Fcount .gt. Zero )) Then
		Command_Line(Bcount:Bcount) = ')'
	        Command_Line(Bcount+1:Bcount+11) = Plot_Dev(1:11)
		Write(Unit=Lun, Fmt=100, Iostat= Status) Command_Line
	        If ( Status .ne. Zero ) Then
		   FXT_Write_Engplot = %loc(FXT_Aberr)
		   Call Lib$Signal(FXT_RMSWritePlt, %val(1), %val(Status))
	        Endif
	    Endif       ! End of loop and valid command line not yet written

	Enddo		! Ix = 1, Table_Max

!	Reset the Plot Table.

	Do Ix = 1, Table_Max
	  Plot_Table(Ix) = .False.
	Enddo

!	If report is requested, write information in report.

	If ((Report) .And. (FXT_Write_Engplot .Eq. %loc(FXT_Normal))) 
	1  Write (Lun_Rpt, 200) Start, Stop

	Return
 100	Format (A254)
 200    Format (1x/5x,'Engplot Command File written for time range: ',
	1       A, ' to ', A )

	End

