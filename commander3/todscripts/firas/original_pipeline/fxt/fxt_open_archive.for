
C-------------------------------------------------------------------------------

	Integer*4 Function FXT_Open_Archive ( File_Seg, Time_Range, 
	1	  Num_Seg, Report, Lun_Rpt, Lun_Eng )

C-------------------------------------------------------------------------------
C
C	Purpose: To open the FDQ_Eng archive for read access.
C
C	Author: Shirley M. Read
C		STX, November, 1988
C
C	Invocation: Status = FXT (File_Seg, Time_Range, Num_Seg, 
C		    Report, Lun_Rpt, Lun_Eng)
C
CH	Change Log:
CH
CH
C	  ----------------------------------------------------------------------
C
C	Input Files:
C	  FDQ_Eng archive files for read access.
C
C	Output Files:
C
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  File_Seg      C*39            FDQ_Eng filename for open archive
C	  Time_Range    C*30            ASCII time range for open archive
C	  Num_Seg       I*4             Number of FDQ_Eng files found in 
C					selected time range
C	  Report        L*1             Flag indicating whether to write report
C	  Lun_Rpt       I*4             Logical unit number for report file
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Lun_Eng       I*4             Logical unit number for FDQ_Eng
C	
C	Subroutines Called:
C
C	  Lib$Get_Lun
C	  Lib$Signal
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C	  FXT_Msg.Txt -- Declarations of external message parameters
C	  SS$Def      -- System errors
C
C	Processing Method:
C
C	  Set the function return status to FXT_Normal.
C         Get a logical unit number to read the FDQ_Eng archive.
C	  If there is an error, Then 
C	      Set the function status to FXT_Aberr.
C	      Signal the error.
C	  Endif
C         If the status is OK, Then
C           Build the complete filename from the archive name and the 
C           filename or the dataset and time range.
C	    Open the FDQ_Eng archive file for read access.
C	    If there is an error, Then
C             Set the function status to FXT_Aberr.
C	      Signal the error.
C	    Endif
C	  Endif
C	  If a report is requested and the status is OK, Then
C	    Write the name of the file opened or the dataset and
C	      time range for the open in the report.
C	  Endif
C	  Return
C
C------------------------------------------------------------------------------
C	

	Implicit None

!	Passed Parameters

	Character*39    File_Seg        ! FDQ_Eng filename
	Character*30    Time_Range      ! Time range for open archive
	Integer*4       Num_Seg		! Number of FDQ_Eng files
	Logical*1       Report          ! Flag for report file
	Integer*4	Lun_Rpt         ! Logical unit number for report
	Integer*4       Lun_Eng         ! Logical unit number for FDQ_Eng

!	Include Files and External Parameters

	Include 'CT$Library:Ctuser.Inc'
	Include '(FXT_Msg)'
	Include '($SSDef)'

	Integer*4 	CT_Connect_Read
	External        CT_Connect_Read

!	Functions

	Integer*4 Lib$Get_Lun    !Get unit number (from system library)

!     Local variables     

	Integer*4       Status	        ! System status
	Integer*4       Iostatus        ! FORTRAN I/O status
	Integer*4       Zero  / 0 /
	Integer*4       One  / 1 /
	Character*15    Arc_Id /'CSDR$FIRAS_RAW:'/
	Character*23    Arc_Parm / 'CSDR$FIRAS_RAW:FDQ_ENG/' /
	Character*60    Arc_File	! Storage for complete filename

!	Set function return status to normal.

	FXT_Open_Archive = %loc(FXT_Normal)

!	Get logical unit number for the FDQ_Eng file.

	Status = Lib$Get_Lun ( Lun_Eng )
	If ( Status .Ne. SS$_Normal ) Then
	    FXT_Open_Archive = %loc(FXT_Aberr)
	    Call Lib$Signal (FXT_GetLunErr, %val(1), %val(Status))
	Endif		! Status is not normal

!       Open FDQ_Eng file for read access.                          

	If ( FXT_Open_Archive .Eq. %loc(FXT_Normal)) Then
	    If ( Num_Seg .Eq. One ) Then
	      Arc_File = Arc_Id // File_Seg
	    Else
	      Arc_File = Arc_Parm // Time_Range
	    Endif
	    
	    Open ( Unit=Lun_Eng,
	1	   File=Arc_File, 
	2	   Status='OLD', Iostat=Iostatus,
	3	   Useropen=CT_Connect_Read)
 
	    If (Iostatus .Ne. Zero) Then
	      FXT_Open_Archive = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_CTOpenErr,%val(1),%val(Iostatus))
	    Endif
	Endif		! Status is normal 

!	If a report is requested, write information in report.

	If ( Report .And. FXT_Open_Archive .Eq. %loc(FXT_Normal)) Then
	    Write (Unit=Lun_Rpt, FMT=100, Iostat=Iostatus) Arc_File
	    If (Iostatus .Ne. Zero) Then
	      Call Lib$Signal( FXT_Writerr, %val(1), %val(Iostatus))
	    Endif
	Endif	! Report

 100	Format(1x/5x,'File Opened for Read Access: ', A )

	Return
	End
