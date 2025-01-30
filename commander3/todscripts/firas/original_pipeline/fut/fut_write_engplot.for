	INTEGER*4 FUNCTION FUT_WRITE_ENGPLOT ( PLOT_TABLE, 
	1	  PLOT_START, PLOT_STOP, ENGLIM, PLOT_DEVICE )
C/
C/ PROGRAM NAME: 
C/	FUT_WRITE_ENGPLOT
C/
C/ PROGRAM DESCRIPTION:
C/	This routine accumulates build command lines for batch mode runs to
C/	produce Engplots of engineering fields that exceeded limits.
C/
C/ AUTHOR:
C/	Shirley M. Read ( STX )  May 1988.
C/
C/ Modified by:
C/      H. Wang, STX, 1/29/91
C/      New requirements for FDQ
C/      Write "command file for generating engplots" out only
C/      if engplots must be produced
C/ 
CH CHANGE LOG:
CH
CH	Version 4.2.1 02/09/89, SPR 2990, Shirley M. Read, STX
CH		FDQ needs a command line option for plot device. The batch mode
CH	        runs should be able to go to the lineprinter or laserprinter.
CH
CH		Version 4.2.3 05/08/89, SPR 3806, Shirley M. Read, STX
CH			Engplots has been modified to plot up to a maximum
CH			of 4 fields per plot so that the plot labels will be
CH			readable. FUT_Write_Engplot generates an Engplots 
CH			command file for FDQ. The maximum number of fields 
CH			has been changed to 4.
CH
CH               
CH      Version 6.1 5/15/90, SPR 6726, H. Wang, STX
CH              FDQ must support talaris 1590t printer
CH
CH      Version 6.7 8/20/90, SPR 67267148, H. Wang, STX
CH              FDQ plot command file need unique name. 
C/
C/ CALLING SEQUENCE:
C/	Status = FUT_Write_Engplot( Plot_Table, Plot_Start, Plot_Stop,
C/		 Englim, Plot_Device )
C/
C/ INPUT PARAMETERS:
C/	Name             Type    Description
C/	-----------------------------------------------------------------
C/	Plot_Table(118)  L*1     Plot table for Engplots.
C/	Plot_Start(2)    I*4     Plot start time.
C/	Plot_Stop(2)     I*4     Plot stop time.
C/	Plot_Device      I*4     Plot device: line or laser printer
C/
C/ OUTPUT PARAMETERS:
C/	EngLim		 I*4	 Flag indicating engineering limit violations.
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
C/	FUT_PLOT_NAMES
C/
C/ SUBROUTINES CALLED:
C/	Lib$Get_Lun
C/
C/ SHARED DATA AREAS:
C/
C/ METHOD USED:
C/	Convert the plot start and stop times to character and insert in 
C/	   command line.
C/	Do for all engineering fields which exceeded red and yellow limits
C/	   Insert the fields to be plotted from the FUT_Plot_Names into the
C/	   command line. When the number reaches 4 or the total characters
C/	   equals 254, write this command line. Start the next command line
C/	   with the new field.	
C/	Enddo
C/
C/ SPECIAL HANDLING:
C/	None
C/
C**************************************************************************************

	IMPLICIT NONE

!	External Parameters and Functions
	
	External   	FUT_NORMAL
	External	FUT_ABERR
	External	FUT_LUNGETERR
	External	FUT_WRITERR
	External	FUT_OPENERR
	External	FUT_CLOSERR

	Integer*4	Lib$Get_Lun	! Get logical unit number

!	Include Files

	Include	        '(FUT_PLOT_NAMES)'
	Include		'($SSDef)'
        Include         '($jpidef)'
	Include		'(fut_params)'

!	Input/Output Passed Parameters

	Logical*1       Plot_Table(118)	  ! Plot flags for Engplots.

	Integer*4       Plot_Start(2) ! Binary quad start time for plot
	Integer*4       Plot_Stop(2)  ! Binary quad stop time for plot
	Integer*4	EngLim
	Integer*4	Plot_Device   ! Plot device: line or laser printer

!	Local Variables
        Integer*2       Items(30)
        Integer*4       L_items(15)
        Integer*4       Pro_Id
        Character*8     PID
        Integer*4       IOSB(2)
        Equivalence (Items(1),L_Items(1)) 
        Integer*4       Sys$getJPI

	Character*254	Command_Line
	Data		Command_Line(1:40) /
	1    '$ENG/NOIN/CONV/NORE/JSTA=yydddhhmmssmmm/' /
	Data		Command_Line(41:66) /
	1    'JSTO=yydddhhmmssmmm/FIEL=('  /
	Character*11    Plot_Dev / '/PLO="/QMS"' /

	Character*230   Blank		! Blank char string
	Character*1	Blanks(188) 	/ 188 * ' ' /
	Equivalence     ( Blank, Blanks(1) )
        character*42    plot_file
	Character*30	Filename / 'Run_Engplot.Com' /

	Integer*4 	Retstat		! Return status
	Integer*4 	Success / 1 /, Error / 2 /  ! Values for status
	Integer*4	Zero    / 0 /   ! Zero value
	Integer*4	Status		! Dummy status variable

	Integer*4 	Lun		! Logical unit number to write text file
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
        Logical*1       OP_FILE/.false./           !
        Logical*1       AL_OPEN/.false./           !
!	Set return status to success.

	Retstat = Success

!	Initialization.
        Do IX = 1,118
          If (PLOT_table(ix)) OP_FILE = .true.
        ENDDO

	
        IF (OP_FILE .and. (.not. al_open)) then   
	    Status = Lib$Get_Lun ( Lun )

 	    If ( Status .ne. SS$_Normal ) Then
		Retstat = Error
		Call Lib$Signal ( FUT_LUNGETERR, %val(1), %val(Status))
	    Endif
	    If (Retstat .eq. Success ) Then
                Items(1) = 4
                Items(2) = JPI$_PID
                L_items(2) = %loc(Pro_ID)
                L_items(3) = 0
                Status = sys$getjpi(,,,%Ref(items),%ref(IOSB),,)
                write(PID(1:8),'(z8)') Pro_Id
                Filename= 'RUN_ENGPLOT_'//PID//'.COM'
		Open ( Unit=Lun, File=Filename, Status= 'NEW', 
	1	  Recl=254, Carriagecontrol='LIST', Iostat=Status )

	        If ( Status .ne. Zero ) Then
		  Retstat = Error
		  Call Lib$Signal ( FUT_OPENERR, %val(2), %val(Status),
	1				Filename )
	        Endif
	    Endif
	    If ( Plot_Device .eq. fac_lineprinter ) Then
		Plot_Dev = '/PLO="/PRI"' 
	    Elseif ( Plot_Device .eq. fac_laserqms ) Then
	        Plot_Dev = '/PLO="/QMS"'
	    Elseif ( Plot_Device .eq. fac_lasertf ) Then
	        Plot_Dev = '/PLO="/TF "'
	    Endif
	    OP_file = .false.
            AL_OPEN=.true.
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
!	field names from FUT_PLOT_NAMES.

	Do Ix = 1, Table_Max
	
	   If (( Plot_Table(Ix)) .and. (Retstat .eq. Success)) Then 
		EngLim = fac_present
		Fcount = Fcount + 1
	        Totcount = Bcount + Name_Len(Ix)
	        If (( Totcount .gt. Max ) .or.
	1	   ( Fcount .gt. Max_Fields)) Then
		   Command_Line(Bcount:Bcount) = ')'
	           Command_Line(Bcount+1:Bcount+11) = Plot_Dev(1:11)
		   Write(Unit=Lun, Fmt=100, Iostat= Status) Command_Line
	           If ( Status .ne. Zero ) Then
		      Retstat = Error
		      Call Lib$Signal( FUT_WRITERR, %val(1), %val(Status))
	           Endif
		   Command_Line(67:254) = Blank(1:188)
	           Bcount = 66
		   Fcount = 1
	        Endif

	        If ( Retstat .eq. Success ) Then
		   Command_Line(Bcount+1:Bcount+Name_Len(Ix)) 
	1		= Plot_Name(Ix)
		   Bcount = Bcount + Name_Len(Ix) + 1
	           Command_Line(Bcount:Bcount) = ','
	        Endif
	   Endif	! Plot_Table(Ix) and Retstat is success

	   If (( Retstat .eq. Success ) .and.
	1	(Ix .eq. Table_Max) .and. ( Fcount .gt. Zero )) Then
		Command_Line(Bcount:Bcount) = ')'
	        Command_Line(Bcount+1:Bcount+11) = Plot_Dev(1:11)
		Write(Unit=Lun, Fmt=100, Iostat= Status) Command_Line
	        If ( Status .ne. Zero ) Then
		   Retstat = Error
		   Call Lib$Signal( FUT_WRITERR, %val(1), %val(Status))
	        Endif
	    Endif       ! End of loop and valid command line not yet written

	Enddo		! Ix = 1, Table_Max

!	Reset the Plot Table.

	Do Ix = 1, Table_Max
	  Plot_Table(Ix) = .false.
	Enddo

!	Set function to return status

	If ( Retstat .eq. Success ) then
	  FUT_Write_Engplot = %loc(FUT_NORMAL)
	Else
	  FUT_Write_Engplot = %loc(FUT_ABERR)
	Endif

	Return
 100	Format (a254)

	End
