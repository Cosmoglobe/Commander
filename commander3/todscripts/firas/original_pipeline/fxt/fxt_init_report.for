C-------------------------------------------------------------------------------
	Integer*4 Function FXT_Init_Report ( Gmt_Start, Gmt_Stop,
	1	  File_Seg, Filter, Init_Xtrm, Plot, Min_Offtime,
	2	  Max_Offtime, Flags, Lun_Rpt, Report_File, CurGMT)
C-------------------------------------------------------------------------------
C
C	Purpose: To write initial information into the report file,
C		 primarily, the options selected from the command line.
C
C	Author: Shirley M. Read
C		STX, November, 1988
C
C	Invocation: Status = FXT_Init_Report ( Gmt_Start, Gmt_Stop,
C	  	    File_Seg, Filter, Init_Xtrm, Plot, Min_Offtime,
C		    Max_Offtime, Flags, Lun_Rpt, Report_File, CurGMT )
C
CH	Change Log:
CH        SPR 6727, FXT must support talaris 1590t printer
CH               H. Wang, STX, may 16, 1990
C	  SPR4171, Standardize report file name
C		Larry P. Rosen, STX, 29 August 1990
C	  ----------------------------------------------------------------------
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  GMT_Start     C*14            Selected start time
C	  GMT_Stop      C*14            Selected stop time
C         File_Seg      C*39            Selected FDQ_Eng segment
C         Filter        I*4             Filter value for extrema check
C         Init_Xtrm     L*1             Flag to initialize extrema record
C         Plot          I*4             Plot type
C         Min_Offtime   I*4             Minimum offset time for plot
C         Max_Offtime   I*4             Maximum offset time for plot
C         Flags         L*1             Indicator to use enable/disable flags
C                                       for extrema checking or not
C         Report_File   C*33		Report file name
C         CurGMT        C*14            Current GMT run time
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C         Lun_Rpt 	I*4             Logical unit for report file  
C	
C	Subroutines Called:
C
C	  Lib$Get_Lun
C         Sys$Gettim
C	  CT_Binary_to_GMT
C	  Lib$Signal
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C	  FXT_Msg.Txt
C         $SSDef
C
C	Processing Method: 
c	  Set the function status to FXT_Normal.
C	  Get a logical unit mnumber for the report file.
C         Get the system time.
C	  Build the report file name using the system time.
C         Open the report file.
C         If the return status is good, then
C           Write the selected options from the command line.
C         Else
C           Signal an error.
C           Set the function status to FXT_Aberr. 
C	  Endif
C         Return
C
CH        VERSION 4.4.1 SER 3306 Q. CHUNG STX 08/09/89
CH                      PROVIDE VERSION NUMBER TO TRACK SOFTWARE UPDATE.
C------------------------------------------------------------------------------
C	
	Implicit None

!	Passed Parameters.

	Character*14 GMT_Start		! Selected start time
	Character*14 GMT_Stop           ! Selected stop time
        Character*39 File_Seg           ! Selected FDQ_Eng segment
	Integer*4    Filter             ! Filter value for extrema check
	Logical*1    Init_Xtrm          ! Flag to initialize extrema record
	Integer*4    Plot               ! Plot type
	Integer*4    Min_Offtime        ! Minimum offset time for plot
	Integer*4    Max_Offtime        ! Maximum offset time for plot
	Logical*1    Flags              ! Indicator to use enable/disable flags
                                        ! for extrema checking or not
        Integer*4    Lun_Rpt            ! Logical unit for report file  
	Character*33 Report_File	! File name for report
	Character*14 CurGMT		! Current run time GMT

!	Include files.

	Include      '(FXT_Msg)'
	Include      '($SSDef)'

!	Functions.

	Integer*4    Lib$Get_Lun

!	Local Declarations.

	Character*12 Plot_Dev(3)  / 'Lineprinter ', 'Laserqms    ',
	1	'Lasertf     '/
	Character*39 Blank              ! Blank characters
	Character*1  Blanks(39) / 39 * ' ' /
	Equivalence  (Blanks(1),Blank)
	Integer*4    Status		! Return status
	Integer*4    Zero / 0 /         ! Zero value for FORTRAN I/O status
	Integer*4    One / 1 /          ! One value
C
 	Integer*4    Cut_Register_Version
       	Integer*4    Cut_Display_Banner
        Integer*4    rstatus
        Integer*4    num_vol/80/
        Character*5  version
        Parameter    (version='4.4.1')

!	Set the function status to Normal.

	FXT_Init_Report = %loc(FXT_Normal)

!	Get a unit number and open the report file.

	Status = Lib$Get_Lun (Lun_Rpt)

	If ( Status .Ne. SS$_Normal ) Then
	      FXT_Init_Report = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_GetLunErr, %val(1), %val(Status))
	Endif

!	Open the report file.

	If ( FXT_Init_Report .Eq. %loc(FXT_Normal) ) Then

	      Open( Unit=Lun_Rpt, File=Report_File, 
	1	    Status='NEW', Form='FORMATTED', Access='SEQUENTIAL',       
	2	    Organization='SEQUENTIAL', Iostat=Status )

	      Rstatus = Cut_Register_Version(version)
              Rstatus = Cut_Display_Banner(Lun_rpt,Num_vol,
	1	'FIRAS Facility FXT_Log_Extrema')
	      If ( Status .Ne. Zero ) Then
		    FXT_Init_Report = %loc(FXT_Aberr)
	            Call Lib$Signal(FXT_OpenErr,%val(1),%val(Status)) 
	      Endif
	Endif	! Function status is normal

!	Write the current time and selected command line options in 
!       the report file.

	If ( FXT_Init_Report .Eq. %loc(FXT_Normal) ) Then
	  Write (Unit=Lun_Rpt, FMT=100, Iostat=Status) Curgmt
	  If (Status .Ne. Zero) Then
	    FXT_Init_Report = %loc(FXT_Aberr)
	    Call Lib$Signal(FXT_WritErr, %val(1),%val(Status))
	  Endif
	Endif		! Function status is normal

	If ( FXT_Init_Report .Eq. %loc(FXT_Normal) ) Then
	  If ( File_Seg .Eq. Blank ) Then
	    Write (Unit=Lun_Rpt, FMT=200, Iostat=Status) 
	1	Gmt_Start, Gmt_Stop
	    If (Status .Ne. Zero) Then
	      FXT_Init_Report = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_WritErr, %val(1),%val(Status))
	    Endif
	  Else
	    Write (Unit=Lun_Rpt, FMT=250, Iostat=Status) File_Seg
	    If (Status .Ne. Zero) Then
	      FXT_Init_Report = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_WritErr, %val(1),%val(Status))
	    Endif
	  Endif
	Endif		! Function status is normal

	If ( FXT_Init_Report .Eq. %loc(FXT_Normal) ) Then
	  Write (Unit=Lun_Rpt, FMT=300, Iostat=Status) Filter
	  If (Status .Ne. Zero) Then
	    FXT_Init_Report = %loc(FXT_Aberr)
	    Call Lib$Signal(FXT_WritErr, %val(1),%val(Status))
	  Endif
	Endif		! Function status is normal

	If ( FXT_Init_Report .Eq. %loc(FXT_Normal) ) Then
	  If ( Init_Xtrm ) Then
	    Write (Unit=Lun_Rpt, FMT=400, Iostat=Status) 
	    If (Status .Ne. Zero) Then
	      FXT_Init_Report = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_WritErr, %val(1),%val(Status))
	    Endif
	  Else
	    Write (Unit=Lun_Rpt, FMT=450, Iostat=Status) 
	    If (Status .Ne. Zero) Then
	      FXT_Init_Report = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_WritErr, %val(1),%val(Status))
	    Endif
	  Endif
	Endif		! Function status is normal

	If ( FXT_Init_Report .Eq. %loc(FXT_Normal) ) Then
	  If ( Plot .Lt. One ) Then
	    Write (Unit=Lun_Rpt, FMT=500, Iostat=Status) 
	    If (Status .Ne. Zero) Then
	      FXT_Init_Report = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_WritErr, %val(1),%val(Status))
	    Endif
	  Else
	    Write (Unit=Lun_Rpt, FMT=550, Iostat=Status) 
	    If (Status .Ne. Zero) Then
	      FXT_Init_Report = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_WritErr, %val(1),%val(Status))
	    Endif
	    If ((Plot .Eq. 1) .Or. (Plot .Eq. 2) .or.
	1	  (plot .eq. 3)) Then
	      Write (Unit=Lun_Rpt, FMT=600, Iostat=Status)
	1	Plot_Dev(Plot), Min_Offtime, Max_Offtime
	    Else
	      Write (Unit=Lun_Rpt, FMT=650, Iostat=Status)
	1	Min_Offtime, Max_Offtime
	      FXT_Init_Report = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_InvPlot)
	    Endif
	    If (Status .Ne. Zero) Then
	      FXT_Init_Report = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_WritErr, %val(1),%val(Status))
	    Endif
	  Endif
	Endif		! Function status is normal

	If ( FXT_Init_Report .Eq. %loc(FXT_Normal) ) Then
	  If ( FLags ) Then
	    Write (Unit=Lun_Rpt, FMT=700, Iostat=Status) 
	    If (Status .Ne. Zero) Then
	      FXT_Init_Report = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_WritErr, %val(1),%val(Status))
	    Endif
	  Else
	    Write (Unit=Lun_Rpt, FMT=750, Iostat=Status) 
	    If (Status .Ne. Zero) Then
	      FXT_Init_Report = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_WritErr, %val(1),%val(Status))
	    Endif
	  Endif
	Endif		! Function status is normal

 100	Format (1x///25x,'FIRAS FXT_LOG_EXTREMA REPORT'//
	1	25x,'Run Time :    ',a/)
 200	Format (1x//5x,'FDQ_Eng Time Range: ' , A, ' to ', A)
 250    Format (1x/5x,'One FDQ_Eng Segment Specified: ', A)
 300    Format (1x/5x,'Filter value for number of extrema ',
	1	'exceeded =', I10 )
 400    Format (1x/5x,'The FXT_Eng_Xtrm Record will be initialized.'//
	1	10x,'The output record defines the minima and ',
	2	'maxima for the current run.')
 450    Format (1x/5x,'The FXT_Eng_Xtrm Record defines the cumulative',
	1	' minima and maxima.')
 500    Format (1x/5x,'No Engplots command file will be generated ',
	1	'for the run.')
 550    Format (1x/5x,'An Engplots command file will be generated ',
	1	'for the run.')
 600    Format (1x/10x,'Plot device: ', A //
	1	10x,'Minimum Offset Time in Hours for the Plots is ',I10//
	1	10x,'Maximum Offset Time in Hours for the Plots is ',I10)
 650    Format (1x/10x,'Invalid plot device selected. Run will abort.'//
	1	10x,'Minimum Offset Time in Hours for the Plots is ',I10//
	1	10x,'Maximum Offset Time in Hours for the Plots is ',I10)
 700    Format (1x/5x,'Enable/Disable Flags from FEX_Limflags will ',
	1	'be used in checking the extrema.')
 750    Format (1x/5x,'Enable/Disable Flags from FEX_Limflags will ',
	1	'not be used in checking the extrema.')

	Return
	End
