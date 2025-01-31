C-------------------------------------------------------------------------------

	Integer*4 Function FXT_Close_Xtrm ( Data_Start, 
	1	  Lun_Xtrm, Report, Lun_Rpt )

C-------------------------------------------------------------------------------
C
C	Purpose: To set the time tag on the FXT_Eng_Xtrm file and close the
C	         extrema archive.
C
C	Author: Shirley M. Read
C		STX, December, 1988
C
C	Invocation: Status = FXT_Close_Xtrm ( Data_Start, Lun_Xtrm, 
C		    Report, Lun_Rpt )
C
CH	Change Log:
CH
C	  ----------------------------------------------------------------------
C
C	Input Files:
C
C	Output Files:
C         FXT_Eng_Xtrm Archive file
C
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Data_Start(2) I*4		Binary data start time
C         Lun_Xtrm      I*4		Logical unit for extrema file
C         Report        L*1		Flag to write report
C         Lun_Rpt       I*4		logical unit for report file
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	
C	Subroutines Called:
C
C	  CCT_Set_Ttg_Time_Range
C	  CT_GMT_to_Binary
C 	  CT_Close_Arcv
C	  Lib$Signal
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C	  FXT_Msg.Txt
C	  CT$Library:CTUser.Inc
C
C	Processing Method:
C
C	  Set the function status to FXT_Normal.
C	  Call CCT_Set_Ttg_Time_Range to set the time tag on the 
C           FXT_Eng_Xtrm file.
C         If the status is not normal, 
C	    Signal the error.
C	    Set the function status to FXT_Aberr.
C	  Endif
C	  Close the FXT_Eng_Xtrm archive.
C         If the status is not normal, 
C	    Signal the error.
C	    Set the function status to FXT_Aberr.
C	  Endif
C	  If a report is requested, log the close of the FXT_Eng_Xtrm archive.
C	  Return.
C
C------------------------------------------------------------------------------
C	
	Implicit None

!	Passed parameters.

	Integer*4      Data_Start(2) 	! Binary data start time
	Integer*4      Lun_Xtrm    	! Logical unit for extrema file
	Logical*1      Report      	! Flag to write report
	Integer*4      Lun_Rpt          ! Logical unit for report file

!	Include files.

	Include '(FXT_Msg)'
	Include 'CT$Library:CTUser.Inc'

!	Functions.

	Integer*4      CCT_Set_Ttg_Time_Range

!	Local Declarations.

	Integer*4      Status           ! Return status
	Integer*4      Zero / 0 /     
	Integer*2      Ct_Status(20)    ! COBETRIEVE return status
	Character*14   GMT_Stop_Tag  / '99365235959990' /
	Integer*4      Stop_Tag(2)      ! CT final stop time tag 

!	Set function status to FXT_Normal.

	FXT_Close_Xtrm = %loc(FXT_Normal)

! 	Set the time tag for the FXT_Eng_Xtrm file.

	Call CT_GMT_to_Binary (GMT_Stop_Tag, Stop_Tag)

	Status = CCT_Set_Ttg_Time_Range (Lun_Xtrm, Data_Start, Stop_Tag)

	If ( .Not. Status ) Then
	  FXT_Close_Xtrm = %loc(FXT_Aberr)
	  Call Lib$Signal(FXT_CCTSetTtgErr, %val(1), %val(Status))
	Endif

!	If function status is normal, close the archive.

	If ( FXT_Close_Xtrm .Eq. %loc(FXT_Normal)) Then
	  Call CT_Close_Arcv (, Lun_Xtrm, CT_Status )
	  If ( CT_Status(1) .Ne. Ctp_Normal ) Then
	    Status = CT_Status(1)
	    FXT_Close_Xtrm = %loc(FXT_Aberr)
	    Call Lib$Signal(FXT_CTClosErr, %val(1), %val(Status))
	    If ( Report ) Then
	      Write (Unit=Lun_Rpt, FMT=100, Iostat=Status) CT_Status(1)
	      If ( Status .Ne. Zero ) Then
	        Call Lib$Signal(FXT_WritErr, %val(1), %val(Status))
	      Endif
	    Endif
	  Else		! CT_Status(1) .Eq. normal
	    If ( Report ) Then
	      Write (Unit=Lun_Rpt, FMT=200, Iostat=Status)
	      If ( Status .Ne. Zero ) Then
	        FXT_Close_Xtrm = %loc(FXT_Aberr)
	        Call Lib$Signal(FXT_WritErr, %val(1), %val(Status))
	      Endif
	    Endif
	  Endif		! CT_Status(1) .Ne. normal
	Endif		! Function status is normal

 100    Format(1x/5x,'FXT_Eng_Xtrm Archive was not closed successfully.',
	1  ' CT_Status = ', I8 )
 200    Format(1x/5x,'FXT_Eng_Xtrm Archive was closed successfully.')

	Return
	End
