
C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Parse_Checksum ( Channel, Start, Stop,
	1                            Type_Rpt, Checksum, Badtime )

C-------------------------------------------------------------------------------
C
C	Purpose: To parse the command line for the FTB_Checksum Program.
C
C	Written by: Nilo G. Gonzales
C		STX, April, 1990
C
C	Invocation: Status = FTB_Parse_Checksum ( Channel, Start, Stop,
C	                      Type_Rpt, Checksum, Badtime )
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
C         Channel       C*2             FIRAS Channel name
C         Start         C*14            Start time for data access
C         Stop          C*14            Stop time for data access 
C	  Type_Rpt      C*1             Type of report:Badtime, Goodtime or Both
C	  Checksum      L*1             Set data checksum flag
C	  Badtime       L*1             Set data badtime flag
C
C	Subroutines Called:
C
C         UPM_Get_Value -- Returns a character string from the command line.
C	  UPM_Present   -- True if a qualifier is on the command line.
C
C	Common Variables Used:
C
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
C	  Invoke UPM_Get_Value to get the FIRAS channel name.
C	  Invoke UPM_Get_Value to get the start time.
C	  Invoke UPM_Get_Value to get the stop time.
C	  Invoke UPM_Get_Value to get the type of report.
C         Invoke UPM_Present   to get the checksum flag, or badtime flag. 
C	  If any bad status is returned, set the return status to abort.
C	  Return with status.
C
C------------------------------------------------------------------------------
C	
	Implicit None

!	Passed Parameters.

	Character*2  Channel	        ! Name of FIRAS channel
	Character*14 Start              ! Start time for data access.
	Character*14 Stop               ! Stop time for data access.
	Character*1  Type_Rpt           ! Type report: goodtime, badtime, both
	Logical*1    Checksum           ! Set data checksum flag
	Logical*1    Badtime            ! Set data badtime flag
	
!	Include Files

	Include '($SSDef)'

!       Functions.

	Integer*4 UPM_Get_Value
	Integer*4 UPM_Present
	Integer*4 Str$Upcase

!	Local Declarations.

	Integer*4 Status       	        ! Function return status
	Integer*2 Txtlen		! Length of string
	Integer*2 Two / 2 /             ! Maximum char in channel name
	Character*14 Txtval             ! Input string
 
!	External Parameters.

	External FTB_Normal
	External FTB_Aberr
	External FTB_NFSNamLen

	FTB_Parse_Checksum = %loc(FTB_Normal)

	Status = UPM_Get_Value ( 'CHANNEL', Channel, Txtlen )
	If ( Status .NE. SS$_Normal ) Then
	  FTB_Parse_Checksum = %loc(FTB_Aberr)
	  Call Lib$Signal (%val(Status))
	Else
	  Status =  Str$Upcase(Channel,Channel)
	  If ( Status .NE. SS$_Normal ) Then
	    FTB_Parse_Checksum = %loc(FTB_Aberr)
	    Call Lib$Signal (%val(Status))
	  Endif
 	Endif
	
	If ( Txtlen .GT. Two ) Then
	  FTB_Parse_Checksum = %loc(FTB_Aberr)
	  Call Lib$Signal ( FTB_NFSNamLen, %val(1), Txtlen )
	Endif

	Start = '85001000000000'
	Stop =  '99365595959900'
	If ( FTB_Parse_Checksum .EQ. %loc(FTB_Normal) ) Then
	  Status = UPM_Get_Value ( 'JSTART', Txtval, Txtlen )
	  If ( Status .NE. SS$_Normal ) Then
	    FTB_Parse_Checksum = %loc(FTB_Aberr)
	    Call Lib$Signal (%val(Status))
	  Else
	    Start(1:Txtlen) = Txtval(1:Txtlen)
	  Endif
	Endif

	If ( FTB_Parse_Checksum .EQ. %loc(FTB_Normal) ) Then
	  Status = UPM_Get_Value ( 'JSTOP', Txtval, Txtlen )
	  If ( Status .NE. SS$_Normal ) Then
	    FTB_Parse_Checksum = %loc(FTB_Aberr)
	    Call Lib$Signal (%val(Status))
	  Else
	    Stop(1:Txtlen) = Txtval(1:Txtlen)
	  Endif
	Endif

	Status = UPM_Get_Value ( 'Report', Type_Rpt, Txtlen )
	If ( Status .NE. SS$_Normal ) Then
	  FTB_Parse_Checksum = %loc(FTB_Aberr)
	  Call Lib$Signal (%val(Status))
	Else
	  Status =  Str$Upcase(Type_Rpt,Type_Rpt)
	  If ( Status .NE. SS$_Normal ) Then
	    FTB_Parse_Checksum = %loc(FTB_Aberr)
	    Call Lib$Signal (%val(Status))
	  Endif
 	Endif

	If ( FTB_Parse_Checksum .EQ. %loc(FTB_Normal) ) Then
	  If ( UPM_Present ('CHECKSUM') ) Then
	    Checksum = .True.
	  Else
	    Checksum = .False.
	  Endif
	Endif

	If ( FTB_Parse_Checksum .EQ. %loc(FTB_Normal) ) Then
	  If ( UPM_Present ('BADTIME') ) Then
	    Badtime = .True.
	  Else
	    Badtime = .False.
	  Endif
	Endif

	Return
	End
