C--------------------------------------------------------------------------
	Integer*4 Function FRD_Gfg_Close ( LUN_HKP, LUN_R_Gain, LUN_L_Gain,
	1                  LUN_R_Fakeit, LUN_L_Fakeit, LUN_Rpt, Report )
C--------------------------------------------------------------------------
C       Purpose: To close the NFS_HKP file, FEX_GAIN.DAT file, 
C                FEX_FAKEIT.DAT file and Report file in the FIRAS archive.
C       
C       Programmer: Nilo G. Gonzales/STX
C                   March 18, 1991
C       Invocation: Status = FRD_Gfg_Close (LUN_HKP, LUN_R_Gain, LUN_L_Gain,
C                            LUN_R_Fakeit,LUN_L_Fakeit, LUN_Rpt, Report )
CH      Change Log:
CH
C--------------------------------------------------------------------------
C       Input Parameters:
C         Name             Type        Description
C         -----------------------------------------------------------------
C 	  LUN_HKP  	    I*4        ! Unit for NFS_HKP data
C  	  LUN_R_Gain	    I*4        ! Unit for output FEX_GAIN.DAT 
C  	  LUN_L_Gain	    I*4        ! Unit for output FEX_GAIN.DAT 
C	  LUN_R_Fakeit      I*4        ! Unit for output FEX_FAKEIT.DAT 
C	  LUN_L_Fakeit      I*4        ! Unit for output FEX_FAKEIT.DAT 
C         LUN_Rpt           I*4        ! Logical unit for report file
C         Report            L*1        ! Flag for report
C
C       Include Files:
C       
C      	  CT$Library:CTUser.Inc
C      	  $SSdef
C	  FUT_Params
C	  FUT_Error
C
C       Externals:
C
C         FRD_Normal
C         FRD_Aberr
C         FRD_CTClosErr
C
C      Processing Method: PDL for FRD_Gfg_CLose
C   
C      Begin
C           Call CT_Close_ARCV to close the NFS_HKP file.
C           
C           If there is an error, then
C              Signal an error and set an error status for the return.
C              If a report is requested, log the error.
C           Else
C              If a report is requested, log the success for closing
C              the file.
C           Endif
C
C           Call CT_Close_ARCV to close the FEX_GAIN.DAT file.
C           
C           If there is an error, then
C              Signal an error and set an error status for the return.
C              If a report is requested, log the error.
C           Else
C              If a report is requested, log the success for closing
C              the file.
C           Endif
C
C           Call CT_Close_ARCV to close the FEX_FAKEIT.DAT file.
C           
C           If there is an error, then
C              Signal an error and set an error status for the return.
C              If a report is requested, log the error.
C           Else
C              If a report is requested, log the success for closing
C              the file.
C           Endif
C   
C           Return with normal or error status
C      End
C--------------------------------------------------------------------------
	Implicit None

C    Passed Parameters:

	Integer*4	LUN_HKP  	! Unit for NFS_HKP data
	Integer*4	LUN_R_Gain	! Unit for output FEX_GAIN.DAT 
	Integer*4	LUN_L_Gain	! Unit for output FEX_GAIN.DAT 
	Integer*4       LUN_R_Fakeit    ! Unit for output FEX_FAKEIT.DAT  
	Integer*4       LUN_L_Fakeit    ! Unit for output FEX_FAKEIT.DAT  
	Integer*4	LUN_Rpt         ! Logical unit for report file
	Logical*1       Report          ! Flag for report 

C    Include files:

      	Include 'CT$Library:CTUser.Inc'
       	Include '($SSdef)'
	Include '(FUT_Params)'
	Include '(FUT_Error)'

C    Externals:

  	External        FRD_Normal
	External        FRD_Aberr
	External        FRD_CTClosErr

C    Functions:

	Integer*4       Lib$Signal  

C    Local Variables:

	Integer*4       Status, Iostatus   ! Return status
	Integer*4       Zero / 0 /         ! Zero value for FORTRAN status
	Integer*4       One / 1 /          ! One value
	Integer*2       CT_Stat(20)        ! COBETRIEVE return status
	
C   Set the function status to Normal.

	FRD_Gfg_Close = %loc(FRD_Normal)

C   Close the matching set of files.
 	
	Call CT_Close_Arcv (, Lun_Hkp, CT_Stat)  ! NFS_HKP data

	If (CT_Stat(1) .Ne. CTP_Normal) Then
	  Status = CT_Stat(1)
	  FRD_Gfg_Close = %loc(FRD_Aberr)
	  Call Lib$Signal (FRD_CTClosErr, %val(1), %val(Status))
	  If (Report) Write(Unit=Lun_Rpt, FMT=100, Iostat=Iostatus) Status
100	  Format(1x/5x,'Error on close of COBETRIEVE NFS_HKP file. ',
	1	 'Status = ',Z8.8 )
	Endif

	Call CT_Close_Arcv (, LUN_R_Gain, CT_Stat)  ! FEX_GAIN_R.DAT data

	If (CT_Stat(1) .Ne. CTP_Normal) Then
	  Status = CT_Stat(1)
	  FRD_Gfg_Close = %loc(FRD_Aberr)
	  Call Lib$Signal (FRD_CTClosErr, %val(1), %val(Status))
	  If (Report) Write(Unit=Lun_Rpt, FMT=200, Iostat=Iostatus) Status
200	  Format(1x/5x,'Error on close of COBETRIEVE FEX_GAIN_R.DAT file. ',
	1	 'Status = ',Z8.8 )
	Endif

	Call CT_Close_Arcv (, LUN_L_Gain, CT_Stat)  ! FEX_GAIN_L.DAT data

	If (CT_Stat(1) .Ne. CTP_Normal) Then
	  Status = CT_Stat(1)
	  FRD_Gfg_Close = %loc(FRD_Aberr)
	  Call Lib$Signal (FRD_CTClosErr, %val(1), %val(Status))
	  If (Report) Write(Unit=Lun_Rpt, FMT=250, Iostat=Iostatus) Status
250	  Format(1x/5x,'Error on close of COBETRIEVE FEX_GAIN_L.DAT file. ',
	1	 'Status = ',Z8.8 )
	Endif

	Call CT_Close_Arcv (, LUN_R_Fakeit, CT_Stat)  ! FEX_FAKEIT_R.DAT data

	If (CT_Stat(1) .Ne. CTP_Normal) Then
	  Status = CT_Stat(1)
	  FRD_Gfg_Close = %loc(FRD_Aberr)
	  Call Lib$Signal (FRD_CTClosErr, %val(1), %val(Status))
	  If (Report) Write(Unit=Lun_Rpt, FMT=300, Iostat=Iostatus) Status
300	  Format(1x/5x,'Error on close of COBETRIEVE FEX_FAKEIT_R.DAT file. ',
	1	 'Status = ',Z8.8 )
	Endif

	Call CT_Close_Arcv (, LUN_L_Fakeit, CT_Stat)  ! FEX_FAKEIT_L.DAT data

	If (CT_Stat(1) .Ne. CTP_Normal) Then
	  Status = CT_Stat(1)
	  FRD_Gfg_Close = %loc(FRD_Aberr)
	  Call Lib$Signal (FRD_CTClosErr, %val(1), %val(Status))
	  If (Report) Write(Unit=Lun_Rpt, FMT=350, Iostat=Iostatus) Status
350	  Format(1x/5x,'Error on close of COBETRIEVE FEX_FAKEIT_L.DAT file. ',
	1	 'Status = ',Z8.8 )
	Endif

	Return
	End
