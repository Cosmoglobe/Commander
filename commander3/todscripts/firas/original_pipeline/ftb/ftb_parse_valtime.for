
C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Parse_Valtime ( Dataset_Name, Start, Stop,
	1	  Pre_IT, Collect, Listall, Dump)

C-------------------------------------------------------------------------------
C
C	Purpose: To parse the command line for the FTB_Valid_Time Program.
C
C	Author: Shirley M. Read
C		STX, November, 1988
C
C	Invocation: Status = FTB_Parse_Valtime(Dataset_Name, Start, Stop,
C		    Pre_IT, Collect, Listall, Dump)
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
C         Dataset_Name  C*10            COBETRIEVE dataset name
C         Start         C*14            Start time for data access
C         Stop          C*14            Stop time for data access 
C	  Pre_IT        L*1             Flag to check the science data run 
C			 		before the I&T problem was corrected
C	  Collect       L*1             Flag to check the midpoint of collect
C					time rather than CT header time
C	  Listall       L*1             Flag to list information for all records
C         Dump          L*1             On/Off switch to dump record.
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
C	  Invoke UPM_Get_Value to get the dataset name.
C	  Invoke UPM_Get_Value to get the start time.
C	  Invoke UPM_Get_Value to get the stop time.
C	  If Qualifier Pre_IT is present, set the Pre_IT flag
C	  If Qualifier Collect is present, set the collect time flag
C	  If Qualifier Listall is present, set the list all records flag
C         If any of the above calls have an error, signal the error and reset
C           the return status to abort.
C	  Return with status.
C
C------------------------------------------------------------------------------
C	
	Implicit None

!	Passed Parameters.

	Character*10 Dataset_Name	! Name of COBETRIEVE dataset.
	Character*14 Start              ! Start time for data access.
	Character*14 Stop               ! Stop time for data access.
	Logical*1    Pre_IT             ! Flag to check Pre IT data before
					! time tag problem was corrected
	Logical*1    Collect            ! Flag to check the midpoint of 
					! collect time
	Logical*1    Listall            ! Flag to list info for all records
        Logical*1    Dump               ! Dump Record On/Off Switch

!	Include Files

	Include '($SSDef)'

!       Functions.

	Integer*4 UPM_Get_Value
	Integer*4 UPM_Present

!	Local Declarations.

	Integer*4 Status       	        ! Function return status
	Integer*2 Txtlen		! Length of string
	Integer*2 Ten / 10 /            ! Maximum NFS dataset name length
	Character*14 Txtval             ! Input string
 
!	External Parameters.

	External FTB_Normal
	External FTB_Aberr
	External FTB_NFSNamLen

	FTB_Parse_Valtime = %loc(FTB_Normal)

	Status = UPM_Get_Value ( 'DATASET', Dataset_Name, Txtlen )
	If ( Status .NE. SS$_Normal ) Then
	  FTB_Parse_Valtime = %loc(FTB_Aberr)
	  Call Lib$Signal (%val(Status))
 	Endif
	
	If ( Txtlen .GT. Ten ) Then
	  FTB_Parse_Valtime = %loc(FTB_Aberr)
	  Call Lib$Signal ( FTB_NFSNamLen, %val(1), Txtlen )
	Endif

	Start = '85001000000000'
	Stop =  '99365595959900'
	If ( FTB_Parse_Valtime .EQ. %loc(FTB_Normal) ) Then
	  Status = UPM_Get_Value ( 'JSTART', Txtval, Txtlen )
	  If ( Status .NE. SS$_Normal ) Then
	    FTB_Parse_Valtime = %loc(FTB_Aberr)
	    Call Lib$Signal (%val(Status))
	  Else
	    Start(1:Txtlen) = Txtval(1:Txtlen)
	  Endif
	Endif

	If ( FTB_Parse_Valtime .EQ. %loc(FTB_Normal) ) Then
	  Status = UPM_Get_Value ( 'JSTOP', Txtval, Txtlen )
	  If ( Status .NE. SS$_Normal ) Then
	    FTB_Parse_Valtime = %loc(FTB_Aberr)
	    Call Lib$Signal (%val(Status))
	  Else
	    Stop(1:Txtlen) = Txtval(1:Txtlen)
	  Endif
	Endif

	If ( FTB_Parse_Valtime .EQ. %loc(FTB_Normal) ) Then
	  If ( UPM_Present ('PRE_IT') ) Then
	    Pre_IT = .True.
	  Else
	    Pre_IT = .False.
	  Endif
	Endif

	If ( FTB_Parse_Valtime .EQ. %loc(FTB_Normal) ) Then
	  If ( UPM_Present ('COLLECT') ) Then
	    Collect = .True.
	  Else
	    Collect = .False.
	  Endif
	Endif

	If ( FTB_Parse_Valtime .EQ. %loc(FTB_Normal) ) Then
	  If ( UPM_Present ('LISTALL') ) Then
	    Listall = .True.
	  Else
	    Listall = .False.
	  Endif
	Endif

	If ( FTB_Parse_Valtime .EQ. %loc(FTB_Normal) ) Then
	  If ( UPM_Present ('DUMP') ) Then
	    Dump = .True.
	  Else
	    Dump = .False.
	  Endif
	Endif

	Return
	End
