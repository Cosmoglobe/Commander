
C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Parse_Binfld ( Field, Start, Stop, 
	1		   Nbins, Bins, Type_Rpt)

C-------------------------------------------------------------------------------
C
C	Purpose: To parse the command line for the FTB_Binfld Program.
C
C	Author: Shirley M. Read
C		STX, January, 1990
C
C	Invocation: Status = FTB_Parse_Binfld (
C			     Field, Start, Stop, Nbins, Bins, Type_Rpt)
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
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C         Field         C*30            FIRAS Engineering Field name
C         Start         C*14            Start time for data access
C         Stop          C*14            Stop time for data access 
C	  Nbins         I*2             Number of bins for value range
C	  Bins(5)       R*4             Values for binning eng fields
C	  Type_Rpt      C*1             Type of report:Engineering,Table or both
C	
C	Subroutines Called:
C
C         UPM_Get_Value -- Returns a character string from the command line.
C         UPM_Get_Word  -- Returns an I*2 number from the command line.
C         UPM_Get_Float -- Returns a R*4 number from the command line.
C	  UPM_Present   -- True if a qualifier is on the command line.
C
C	Common Variables Used:
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C	
C         $SSDEF
C
C	Processing Method:
C
C         Set return status to normal.
C	  Invoke UPM_Get_Value to get the FIRAS engineering field name.
C	  Invoke UPM_Get_Value to get the start time.
C	  Invoke UPM_Get_Value to get the stop time.
C	  Invoke UPM_Get_Word to get the number of bins.
C	  Invoke UPM_Get_Float to get the values for bins.
C	  Invoke UPM_Get_Value to get the type of report.
C	  If any bad status is returned, set the return status to abort.
C	  Return with status.
C
C------------------------------------------------------------------------------
C	
	Implicit None

!	Passed Parameters.

	Character*30 Field	        ! Name of FIRAS engineering field
	Character*14 Start              ! Start time for data access.
	Character*14 Stop               ! Stop time for data access.
	Integer*2    Nbins   		! Number of bins
	Real*4       Bins(5)            ! Values for bins          
	Character*1  Type_Rpt           ! Type report: engineering, table, both
	
!	Include Files

	Include '($SSDef)'
	Include '(UPM_Stat_Msg)'

!       Functions.

	Integer*4 UPM_Get_Value
	Integer*4 UPM_Get_Word
	Integer*4 UPM_Get_Float
	Integer*4 UPM_Present
	Integer*4 Str$Upcase

!	Local Declarations.

	Integer*4 Status       	        ! Function return status
	Integer*2 Txtlen		! Length of string
	Integer*2 Five / 5 /            ! Maximum number of bins allowed

	Character*14 Txtval             ! Input string
	Integer*2 Count                 ! Counter for bin values
 
!	External Parameters.

	External FTB_Normal
	External FTB_Aberr
	External FTB_NumBins

	FTB_Parse_Binfld = %loc(FTB_Normal)

	Status = UPM_Get_Value ( 'FIELD', Field, Txtlen )
	If ( Status .NE. SS$_Normal ) Then
	  FTB_Parse_Binfld = %loc(FTB_Aberr)
	  Call Lib$Signal (%val(Status))
	Else
	  Status =  Str$Upcase(Field,Field)
	  If ( Status .NE. SS$_Normal ) Then
	    FTB_Parse_Binfld = %loc(FTB_Aberr)
	    Call Lib$Signal (%val(Status))
	  Endif
 	Endif
	
	Start = '85001000000000'
	Stop =  '99365595959900'
	If ( FTB_Parse_Binfld .EQ. %loc(FTB_Normal) ) Then
	  Status = UPM_Get_Value ( 'JSTART', Txtval, Txtlen )
	  If ( Status .NE. SS$_Normal ) Then
	    FTB_Parse_Binfld = %loc(FTB_Aberr)
	    Call Lib$Signal (%val(Status))
	  Else
	    Start(1:Txtlen) = Txtval(1:Txtlen)
	  Endif
	Endif

	If ( FTB_Parse_Binfld .EQ. %loc(FTB_Normal) ) Then
	  Status = UPM_Get_Value ( 'JSTOP', Txtval, Txtlen )
	  If ( Status .NE. SS$_Normal ) Then
	    FTB_Parse_Binfld = %loc(FTB_Aberr)
	    Call Lib$Signal (%val(Status))
	  Else
	    Stop(1:Txtlen) = Txtval(1:Txtlen)
	  Endif
	Endif

	Status = UPM_Get_Word ( 'Numbins', Nbins )
	If ( Status .NE. SS$_Normal ) Then
	  FTB_Parse_Binfld = %loc(FTB_Aberr)
	  Call Lib$Signal (%val(Status))
	Else
	  If (Nbins .GT. Five) Then
	    FTB_Parse_Binfld = %loc(FTB_Aberr)
	    Call Lib$Signal (FTB_NumBins)
	  Endif
 	Endif

	Count = 0
	Status = UPM_Comma
	If ( FTB_Parse_Binfld .EQ. %loc(FTB_Normal) ) Then
	  Do While (Status .EQ. UPM_Comma .AND. Count .LT. Nbins)
	    Count = Count + 1
	    Status = UPM_Get_Float ( 'BinVals', Bins(Count) )
 	  Enddo
	Endif

	Status = UPM_Get_Value ( 'Report', Type_Rpt, Txtlen )
	If ( Status .NE. SS$_Normal ) Then
	  FTB_Parse_Binfld = %loc(FTB_Aberr)
	  Call Lib$Signal (%val(Status))
	Else
	  Status =  Str$Upcase(Type_Rpt,Type_Rpt)
	  If ( Status .NE. SS$_Normal ) Then
	    FTB_Parse_Binfld = %loc(FTB_Aberr)
	    Call Lib$Signal (%val(Status))
	  Endif
 	Endif

	Return
	End
