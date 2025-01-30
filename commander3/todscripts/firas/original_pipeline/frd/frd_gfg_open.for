C--------------------------------------------------------------------------
	Integer*4 Function FRD_Gfg_Open ( LUN_HKP, LUN_R_Gain, LUN_L_Gain,
	1     LUN_R_Fakeit,LUN_L_Fakeit, Report_File, LUN_Rpt, Report, CMD_Line,
	2     GMT_Start, GMT_Stop)
C--------------------------------------------------------------------------
C       Purpose: To open a logical unit number for each of the selected 
C                input and output.
C       
C       Programmer: Nilo G. Gonzales/STX
C                   March 25, 1991
C       Invocation:
C	Integer*4 Function FRD_Gfg_Open ( LUN_HKP, LUN_R_Gain, LUN_L_Gain,
C	1     LUN_R_Fakeit, LUN_L_Fakeit, Report_File, LUN_Rpt, Report,CMD_Line,
C	2     GMT_Start, GMT_Stop
CH      Change Log:
CH
C--------------------------------------------------------------------------
C       Output Parameters:
C         Name             Type        Description
C         -----------------------------------------------------------------
C 	  LUN_HKP  	    I*4 	! Unit numbers for NFS_HKP CT archive
C  	  LUN_R_Gain	    I*4 	! Unit numbers for FEX_GAIN CT archive
C  	  LUN_L_Gain	    I*4 	! Unit numbers for FEX_GAIN CT archive
C	  LUN_R_Fakeit      I*4 	! Unit numbers for FEX_FAKEIT archive 
C	  LUN_L_Fakeit      I*4 	! Unit numbers for FEX_FAKEIT archive 
C         LUN_Rpt           I*4		! Report unit number
C         Report_File       C*21	! Report file name
C         Report            L*1		! Flag to enable writing a report file
C	  CMD_Line(79)      C*79	! Command Line string with defaults
C	  GMT_Start	    C*14	! Data start time
C	  GMT_Stop	    C*14	! Data stop time
C
C       Subroutine Called:
C	  Lib$Get_Lun
C	  Lib$Locc
C	  Lib$Signal
C	  Lib$Establish
C
C       Include Files:
C      	  CT$Library:CTUser.Inc
C      	  $SSdef
C	  FUT_Error
C	  $jpidef
C
C  	Externals:
C	  CT_Connect_Read
C	  CT_Connect_Write
C	  FUT_Error
C	  Cut_translate_archive_id
C	  FRD_Normal
C	  FRD_Aberr
C	  FRD_GetLunErr
C	  FRD_CTOpenErr
C	  FRD_CTinit
C
C      Processing Method: PDL for FRD_Gfg_Open
C   
C      Begin
C           Get LUN for input archive dataset
C           Open housekeeping data for read
C           Get LUN for output retrieve dataset (FEX_FAKEIT.DAT)
C	    Get LUN for output retrieve dataset (FEX_GAIN.DAT)
C   	    Open FEX_FAKEIT.DAT for write
C	    Open FEX_GAIN.DAT for write
C	      
C	    If report is needed Then
C	       Get LUN for report file
C	       Open report file 
C              Write initial information to report file
C           Endif
C      End
C--------------------------------------------------------------------------
	Implicit None

C	Passed Parameters:

	Integer*4	LUN_HKP  	! Unit numbers for NFS_HKP CT archive
	Integer*4	LUN_R_Gain	! Unit numbers for FEX_GAIN CT archive
	Integer*4	LUN_L_Gain	! Unit numbers for FEX_GAIN CT archive
	Integer*4       LUN_R_Fakeit    ! Unit numbers for FEX_FAKEIT archive
	Integer*4       LUN_L_Fakeit    ! Unit numbers for FEX_FAKEIT archive
	Character*21    Report_File     ! Report file name 
	Integer*4	LUN_Rpt         ! Report unit number
	Logical*1       Report          ! Flag to enable writing a report file
	Character*79	CMD_Line(79)	! Command Line string with defaults
	Character*14	GMT_Start
	Character*14	GMT_Stop

C 	Include files:

      	Include 'CT$Library:CTUser.Inc'
       	Include '($SSdef)'
	Include '(FUT_Error)'
	INCLUDE '($jpidef)'

C  	Functions:

	Integer*4       Lib$Get_Lun
	Integer*4       Lib$Locc
	Integer*4       Lib$Signal
	Integer*4	Lib$Establish

C  	Externals:

	Integer*4       CT_Connect_Read
	External        CT_Connect_Read
	Integer*4       CT_Connect_Write
	External        CT_Connect_Write
	External	FUT_Error
	INTEGER*4	lib$getjpi
	INTEGER*4	cut_register_version
	INTEGER*4	cut_display_banner
	INTEGER*4	sys$asctim
	INTEGER*4 	cut_translate_archive_id
	EXTERNAL  	cut_translate_archive_id
	EXTERNAL        FRD_Normal
	EXTERNAL        FRD_Aberr
	EXTERNAL        FRD_GetLunErr
	EXTERNAL        FRD_CTOpenErr
	EXTERNAL        FRD_CTinit

C     	Local Variables:

	CHARACTER*6	version
	PARAMETER	(version='8.7')
	INTEGER*4	num_vol/80/
	Integer*4       Status, Iostatus   ! Return status
	Integer*4       Zero / 0 /         ! Zero value for FORTRAN status
	Integer*4       One / 1 /          ! One value
	Character*22    HKP_File / 'CSDR$FIRAS_RAW:NFS_HKP' /
	Character*31    R_FAKE_File / 'CSDR$FIRAS_REF:FEX_FAKEIT_R.DAT' /
	Character*29    R_GAIN_File / 'CSDR$FIRAS_REF:FEX_GAIN_R.DAT' /
	Character*31    L_FAKE_File / 'CSDR$FIRAS_REF:FEX_FAKEIT_L.DAT' /
	Character*29    L_GAIN_File / 'CSDR$FIRAS_REF:FEX_GAIN_L.DAT' /
	Character*64    HKP_TFile
	Integer*2	Ct_Stat(20)
	CHARACTER*8	owner			! invoking user name
	INTEGER*4	current_time(2)		! current system adt time
	INTEGER*2	time_len		! length of time string
	CHARACTER*32	time			! current system time string
	INTEGER*4	rstatus			! return status
	CHARACTER*72	logn, flogn, tlogn	! logical name and translted
	INTEGER*4	flen, tlen		! length of translated logicals
	
C  	Set the function status to Normal.

	FRD_Gfg_Open = %loc(FRD_Normal)

C      Error checking for the Logical Unit Number of the Input and Output.
 	
	Status = Lib$Get_Lun (LUN_HKP)
	If ( Status .Ne. SS$_Normal ) Then
	   FRD_Gfg_Open = %loc(FRD_Aberr)
	   Call Lib$Signal(FRD_GetLunErr, %val(1), %val(Status))
	Endif

	Status = Lib$Get_Lun (LUN_R_Gain)
	If ( Status .Ne. SS$_Normal ) Then
	   FRD_Gfg_Open = %loc(FRD_Aberr)
	   Call Lib$Signal(FRD_GetLunErr, %val(1), %val(Status))
	Endif

	Status = Lib$Get_Lun (LUN_L_Gain)
	If ( Status .Ne. SS$_Normal ) Then
	   FRD_Gfg_Open = %loc(FRD_Aberr)
	   Call Lib$Signal(FRD_GetLunErr, %val(1), %val(Status))
	Endif

	Status = Lib$Get_Lun (LUN_R_Fakeit)
	If ( Status .Ne. SS$_Normal ) Then
	   FRD_Gfg_Open = %loc(FRD_Aberr)
	   Call Lib$Signal(FRD_GetLunErr, %val(1), %val(Status))
	Endif

	Status = Lib$Get_Lun (LUN_L_Fakeit)
	If ( Status .Ne. SS$_Normal ) Then
	   FRD_Gfg_Open = %loc(FRD_Aberr)
	   Call Lib$Signal(FRD_GetLunErr, %val(1), %val(Status))
	Endif

	Status = Lib$Get_Lun (LUN_Rpt)
	If ( Status .Ne. SS$_Normal ) Then
	   FRD_Gfg_Open = %loc(FRD_Aberr)
	   Call Lib$Signal(FRD_GetLunErr, %val(1), %val(Status))
	Endif

C      Open a report file.

	If ( FRD_Gfg_Open .Eq. %loc(FRD_Normal) ) Then
	   IF (Report) Then
  	       Open(Unit=LUN_Rpt,File=Report_File,Status='NEW',
	1          Access='Sequential',Form='Formatted',Iostat=Iostatus)
	       FUT_Report_Lun = LUN_Rpt
	       Status = Lib$Establish (FUT_Error)
	       If ( Iostatus .Ne. Zero ) Then
	          FRD_Gfg_Open = %loc(FRD_Aberr)
	          Call Lib$Signal(FRD_CTOpenErr,%val(1), %val(Iostatus))
	       Else
	          status = cut_register_version(version)
                  status  = cut_display_banner(lun_rpt,num_vol,
	1             'FIRAS facility FRD_GFG')
	          WRITE(lun_rpt,10)
  10	          FORMAT(/)

C  Write user and current time to the report file.

	          status = lib$getjpi (jpi$_username,,,,owner,)
	          CALL sys$gettim ( current_time )
	          status = sys$asctim ( time_len, time, current_time, 0 )
	          WRITE (lun_rpt,20) owner, time
  20	          FORMAT (' Run by:   ', A, '   at  Time: ',A,/)
	          logn(1:15) = 'CSDR$FIRAS_RAW:'
	          rstatus = cut_translate_archive_id(logn,flogn,flen,tlogn,tlen)
	          WRITE (lun_rpt,40) logn(1:15), tlogn
  40	FORMAT (1X, 'Logical Translation for Input Archive:',/,4X,A,' = ',A)
	          logn(1:15) = 'CSDR$FIRAS_REF:'
	          rstatus = cut_translate_archive_id(logn,flogn,flen,tlogn,tlen)
	          WRITE (lun_rpt,50) logn(1:15), tlogn
  50	FORMAT(1X,'Logical Translation for Output Archive:',/,4X,A,' = ',A)
	          WRITE (lun_rpt,60) CMD_Line(1), CMD_Line(2)
  60    Format (1X,'Command Line with defaults: ',/,10X,A,/,10X,A,/)
	       Endif
	   Endif
C
	   CALL ct_init( ct_stat )
	   IF (ct_stat(1) .NE. ctp_normal) THEN
	      status = ct_stat(1)
	      CALL lib$signal(frd_ctinit,%val(1),%val(status))
	      FRD_GFG_OPEN = %loc(frd_aberr)
	   ENDIF
	EndIf

C      Open CT_Connect_Read for NFS_HKP dataset.

	If ( FRD_Gfg_Open .Eq. %loc(FRD_Normal) ) Then
	   HKP_TFile = Hkp_File // '/' // GMT_Start // ';'
	1   // GMT_Stop // ';'
	   Open(Unit=LUN_HKP,File=HKP_TFile,Status='OLD',
	1      Iostat=Iostatus,Useropen=CT_Connect_Read)
	    If ( Iostatus .Ne. Zero ) Then
	       FRD_Gfg_Open = %loc(FRD_Aberr)
	       Call Lib$Signal(FRD_CTOpenErr,%val(1), %val(Iostatus))
	    Else
	       If (Report) Write(Unit=Lun_Rpt,Iostat=Iostatus,FMT=200)
	1	    HKP_TFile
	    Endif
	Endif

C       Open CT_Connect_Write for FEX_Gain.Dat dataset.

	 If ( FRD_Gfg_Open .Eq. %loc(FRD_Normal)) Then
	    Open(Unit=LUN_R_Gain, File=R_Gain_File,
	1       Status='NEW',Iostat=Iostatus, Useropen=CT_Connect_Write)
	    If ( Iostatus .Ne. Zero ) Then
	        FRD_Gfg_Open = %loc(FRD_Aberr)
	        Call Lib$Signal(FRD_CTOpenErr,%val(1), %val(Iostatus))
	    Else
	        If (Report) Write(Unit=Lun_Rpt,Iostat=Iostatus,FMT=250)
	1	    R_Gain_File
	    Endif
	Endif

	 If ( FRD_Gfg_Open .Eq. %loc(FRD_Normal)) Then
	    Open(Unit=LUN_L_Gain, File=L_Gain_File,
	1       Status='NEW',Iostat=Iostatus, Useropen=CT_Connect_Write)
	    If ( Iostatus .Ne. Zero ) Then
	        FRD_Gfg_Open = %loc(FRD_Aberr)
	        Call Lib$Signal(FRD_CTOpenErr,%val(1), %val(Iostatus))
	    Else
	        If (Report) Write(Unit=Lun_Rpt,Iostat=Iostatus,FMT=250)
	1	    L_Gain_File
	    Endif
	Endif

C       Open CT_Connect_Write for FEX_Fakeit.Dat dataset.

	If ( FRD_Gfg_Open .Eq. %loc(FRD_Normal)) Then
	     Open(Unit=LUN_R_Fakeit,File=R_Fake_File,
	1        Status='NEW',Iostat=Iostatus, Useropen=CT_Connect_Write)
	     If ( Iostatus .Ne. Zero ) Then
	        FRD_Gfg_Open = %loc(FRD_Aberr)
	        Call Lib$Signal(FRD_CTOpenErr,%val(1), %val(Iostatus))
	     Else
	        If (Report) Write(Unit=Lun_Rpt,Iostat=Iostatus,FMT=250) 
	1	    R_Fake_File
	     Endif
	Endif

	If ( FRD_Gfg_Open .Eq. %loc(FRD_Normal)) Then
	     Open(Unit=LUN_L_Fakeit,File=L_Fake_File,
	1        Status='NEW',Iostat=Iostatus, Useropen=CT_Connect_Write)
	     If ( Iostatus .Ne. Zero ) Then
	        FRD_Gfg_Open = %loc(FRD_Aberr)
	        Call Lib$Signal(FRD_CTOpenErr,%val(1), %val(Iostatus))
	     Else
	        If (Report) Write(Unit=Lun_Rpt,Iostat=Iostatus,FMT=250) 
	1	    L_Fake_File
	     Endif
	Endif

 200	Format(1x,'Archive file successfully opened: ',/,10X,A)
 250	Format(1x,'Archive file successfully opened: ',A,/)

	  Return
	  End

