C-------------------------------------------------------------------------------
	Integer*4 Function FPP_Open_Archive ( Chan_Num, In_Files, Out_Files,
	1	Report, Lun_Rpt, Lun_Anc, Lun_In_Sci, Lun_Out_Sci, Day, Maxn ) 
C-------------------------------------------------------------------------------
C	Purpose: To open the science channel input and output file and
C		corresponding Word 31 ancillary file.
C
C	Author: Shirley M. Read
C		STX, February, 1989
C
C	Invocation: Status = FPP_Open_Archive ( Chan_Num, In_Files, Out_Files,
C		Report, Lun_Rpt, Lun_Anc, Lun_In_Sci, Lun_Out_Sci, Day, Maxn )
C
CH	Change Log:
CH
CH		Version 4.4.1 07/22/89, SER 4168, R. Kummerer, STX
CH		There have been problems during test operations for FIRAS
CH		processing due to the required clean-up of the archives after an
CH		FPP or FDQ abort. The FPR tracking system compounds the problems
CH		Files with non-matching version numbers seem often to result
CH		from improper clean-up. Bad record times cause SEGCTL to abort
CH		and mess up the tracking system. It was decided to change the
CH		modify of the science records in FPP and FDQ to a simple
CH		COBETRIEVE read of the existing records from a dataset and write
CH		a modifed dataset with the same information which was entered on
CH		the modify. Two new science data sets will be required: a science
CH		dataset of raw science data plus FPP input and a science dataset
CH		with FPP science data plus FDQ input. These datasets will be
CH		FPP_SDF_xx, where xx is the channel id (RH, RL, LH or LL) and
CH		FDQ_SDF_xx, where xx is the channel id. The new datasets must be
CH		opened and processed in FPP and FDQ.
CH
CH		Version 4.4.1 08/21/89 SER 4210, R. Kummerer, STX
CH			Prevent overlaps in raw science segments.
CH
CH		Version 4.4.4 11/29/89, SPR 5165, R. Kummerer STX
CH			Fails to skip missing segments.
CH
CH		New Version  March 1991, Larry Rosen, STX
CH		New requirements/design modifications: open all input and output
CH		files.  No tracking.
C	-----------------------------------------------------------------------
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Maxn		  I*2		Max number days
C	  Chan_Num        I*2 		Current channel number
C	  In_Files(maxn,5)  C*64	Input science and ancillary filenames
C	  Out_Files(maxn,4) C*64	Output science filenames
C	  Report          L*1           Flag for report
C	  Lun_Rpt 	  I*4           Logical unit for report file
C	  Day             I*2		Day for file open and close
C
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Lun_Anc         I*4           Unit for Word 31 ancillary data 
C	  Lun_In_Sci      I*4           Unit for input science channel data
C	  Lun_Out_Sci     I*4           Unit for output science channel data
C	
C	Subroutines Called:
C
C	  Lib$Get_Lun
C	  Lib$Locc
C	  Lib$Signal
C
C	Include Files:
C	  CT$Library:CTUser.Inc
C	  FPP_Msg.Txt
C         $SSDef
C
C	Processing Method: PDL for FPP_Open_Archive
C
C	If it is the first time in this routine, then
C	   Get a logical unit number for accessing the NFS_ANC archive.
C	   Get a logical unit number for accessing the NFS_SDF channel archive.
C	   Get a logical unit number for accessing the FPP_SDF channel archive.
C	   Set the first time flag to false.
C	Endif
C
C	Open the NFS_ANC file for read access.
C	If there is an error, then
C	   Signal an error and set an error status for the return.
C	Else
C	   If a report is requested, log the success for opening the file.
C	Endif
C	Open the NFS_SDF file for read access.
C	If there is an error, then
C	   Signal an error and set an error status for the return.
C	Else
C	   If a report is requested, log the success for opening the file.
C	Endif
C	Open the FPP_SDF file for read access.
C	If there is an error, then
C	   Signal an error and set an error status for the return.
C	Else
C	   If a report is requested, log the success for opening the file.
C	Endif
C
C	Return with normal or error status. 
C------------------------------------------------------------------------------
	Implicit None

C  Passed Parameters

	Integer*2     Chan_Num		! Current channel number
	Logical*1     Report		! Flag to enable writing a report
	Integer*4     Lun_Rpt		! Logical unit for report file  
	Integer*4     Lun_Anc		! Logical unit for Word 31 NFS_ANC data
	Integer*4     Lun_In_Sci	! Logical unit for input science data
	Integer*4     Lun_Out_Sci	! Logical unit for output science data
	Integer*2     Day		! Day of processing (file number)
	Integer*2     Maxn		! Max number days
	Character*64  In_Files(maxn,5)	! Input science and ancillary file names
	Character*64  Out_Files(maxn,4)	! Output science file names

C  Include files

	Include      'CT$Library:CTUser.Inc'
	Include      '(FPP_Msg)'
	Include      '($SSDef)'

C  Functions

	Integer*4    Lib$Get_Lun
	Integer*4    Lib$Locc

C  Externals

	Integer*4    CT_Connect_Read
	External     CT_Connect_Read	
	Integer*4    CT_Connect_Write
	External     CT_Connect_Write

C  Local Declarations

	Integer*4    Status, Iostatus	! Return status
	Integer*4    Zero / 0 /         ! Zero value for FORTRAN I/O status
	Logical*1    First_Time / .True. /  ! First time in routine

C  Set the function status to Normal.

	FPP_Open_Archive = %loc(FPP_Normal)

C  Get unit numbers for the NFS_ANC, NFS_SDF, and FPP_SDF files.

	If ( First_Time ) Then
	    Status = Lib$Get_Lun (Lun_Anc)
	    If ( Status .Ne. SS$_Normal ) Then
	        FPP_Open_Archive = %loc(FPP_Aberr)
	        Call Lib$Signal(FPP_GetLunErr, %val(1), %val(Status))
	    Endif
	    If ( FPP_Open_Archive .Eq. %loc(FPP_Normal) ) Then
	        Status = Lib$Get_Lun (Lun_In_Sci)
	        If ( Status .Ne. SS$_Normal ) Then
	            FPP_Open_Archive = %loc(FPP_Aberr)
	            Call Lib$Signal(FPP_GetLunErr, %val(1), %val(Status))
	        Endif
	    Endif
	    If ( FPP_Open_Archive .Eq. %loc(FPP_Normal) ) Then
	        Status = Lib$Get_Lun (Lun_Out_Sci)
	        If ( Status .Ne. SS$_Normal ) Then
	            FPP_Open_Archive = %loc(FPP_Aberr)
	            Call Lib$Signal(FPP_GetLunErr, %val(1), %val(Status))
	        Endif
	    Endif
	    First_Time = .False.
	Endif		! First time

C  Open the set of files.

	Open(Unit=Lun_Anc,File=In_Files(Day,5),Status='OLD',
	1	  Iostat=Iostatus,Useropen=CT_Connect_Read)
	If ( Iostatus .Ne. Zero ) Then
	    FPP_Open_Archive = %loc(FPP_Aberr)
	    Call Lib$Signal(FPP_CTOpenErr,%val(1), %val(Iostatus))
	Else
	    If (Report) Write(Unit=Lun_Rpt,Iostat=Iostatus,FMT=10) 
	1	    In_Files(Day,5)
	Endif
	If ( FPP_Open_Archive .Eq. %loc(FPP_Normal)) Then
	    Open(Unit=Lun_In_Sci,File=In_Files(Day,Chan_Num),
	1	    Status='OLD', Iostat=Iostatus, Useropen=CT_Connect_Read)
	    If ( Iostatus .Ne. Zero ) Then
	      FPP_Open_Archive = %loc(FPP_Aberr)
	      Call Lib$Signal(FPP_CTOpenErr,%val(1), %val(Iostatus))
	    Else
	      If (Report) Write(Unit=Lun_Rpt,Iostat=Iostatus,FMT=10) 
	1	      In_Files(Day,Chan_Num)
	    Endif
	Endif		! FPP_Open_Archive status is normal
	If ( FPP_Open_Archive .Eq. %loc(FPP_Normal)) Then
	    Open(Unit=Lun_Out_Sci,File=Out_Files(Day,Chan_Num),
	1	    Status='NEW',Iostat=Iostatus, Useropen=CT_Connect_Write)
	    If ( Iostatus .Ne. Zero ) Then
	      FPP_Open_Archive = %loc(FPP_Aberr)
	      Call Lib$Signal(FPP_CTOpenErr,%val(1), %val(Iostatus))
	    Else
	      If (Report) Write(Unit=Lun_Rpt,Iostat=Iostatus,FMT=10) 
	1	      Out_Files(Day,Chan_Num)
	    Endif
	Endif		! FPP_Open_Archive status is normal

  10	Format(5X,'Archive file successfully opened: ',/,15X,A)

 	Return
	End
