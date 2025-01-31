	Integer*4  Function  FUT_FCS_REFERENCE  ( Arch_Ref, G_Intercept,
	1                                         G_Slope, Chan, Scan, Report,
	2                                         Lun_Rpt )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                            FUT_FMS_REFERENCE.FOR
C
C  Purpose:  Open, read, and close Glitch model file: FEX_GLTCHCOR.DAT
C     Used by FCS and FMS.
C
C  Author:  Larry P. Rosen, Hughes STX, July 1994.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Include files.

	Include		'($SSDef)'
	Include		'(CCT_Get_Config)'
	Include		'CSDR$Library:CTParams.inc'
	Include		'(CCT_Query_Catalog_Record)'

C  Passed parameters.

	Logical*1	Report			! Whether or not to report
	Integer*4	Lun_Rpt			! Report logical unit #
	Character*15	Arch_Ref		! C or D Vector archive
	Character*2	Chan, Scan		! channel, scan mode
	Real*8		G_Slope			! Glitch model slope
	Real*8		G_Intercept		! Glitch model intercept

C  Functions

	Integer*4	Lib$Get_Lun
	Integer*4	CCT_Query_Catalog
	Integer*4	CCT_Open_Config
	integer*4	CCT_Get_Config_TOD
	Integer*4	CCT_Close_Config

C  Local

	Character*16	Ref_File /'FEX_GLTCHCOR.DAT'/	! Reference file name
	Integer*4	Rstat				! Return status
	Integer*4	Lun			! Logical Unit Number
	Dictionary 'ccm_cme_catalog_entry'
	Record /query_catalog/		Query_Cat
	Record /CCM_CME_Catalog_Entry/	Cats (1)
	Integer*4	Access_Time (2)		! time for reference file
	Integer*2	Snum, Cnum		! scan number, chan number
	Integer*2	Index			! index to glitch values

C  Config stuff

	Character*14	Config_GMT_Start / '86001000001000' /
	Character*14	Config_GMT_Stop / '99365235958990' /
	Integer*4	Config_Start (2), Config_Stop (2)
	Character*38	Config_Name (1)
	Integer*4	Config_Size (1) / 256 /
	Integer*4	Config_Lun (1), Config_Index (1)
	Logical*1	Config_New_Segment (1)
	Record /Config_Status/ Config_Status(1)
	Integer*4	Config_Ref_Count
	Dictionary 'FEX_GLTCHCOR'
	Structure /Config_Record/
	   record /FEX_GLTCHCOR / GLTCH_REC
	EndStructure
	Record 	/Config_Record/ Config_Record

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FUT_FCS_Reference = 0

c Query catalog for time of input skymap for use in getting reference files.

	Query_Cat.Archive_ID = Arch_Ref
	Query_Cat.Filename = Ref_File
	Rstat = CCT_Query_Catalog ( Query_Cat, Cats(1) )
	If (.NOT. Rstat) Then
	   FUT_FCS_Reference = 1
	   Write (6,*) 'Error - Error querying catalog.  Status = ', Rstat
	   If (Report) Write (Lun_Rpt,*)
	1     'Error - Error querying catalog.  Status = ', Rstat
	Else
	   Access_Time (1) = Cats(1).Initial_Time (1)
	   Access_Time (2) = Cats(1).Initial_Time (2)

C  Get logical unit numbers, open, read, and close c and d vectors for each
C  input skymap.

	   Rstat = Lib$Get_Lun (Lun)
	   If (Rstat .NE. SS$_Normal) Then
	      FUT_FCS_Reference = 1
	      Write (6,*) 'Error - Getting Logical unit number.  Status = ',
	1        Rstat
	      If (Report) Write (Lun_Rpt,*)
	1        'Error - Getting Logical unit number.  Status = ', Rstat
	   Else
	      Call CT_GMT_To_Binary (Config_GMT_Start, Config_Start)
	      Call CT_GMT_To_Binary (Config_GMT_Stop, Config_Stop)
	      Config_Name (1) = Arch_Ref // Ref_File
	      Config_Lun (1) = Lun
	      Rstat = CCT_Open_Config ( Config_Start, Config_Stop, 1,
	1                               Config_Name, Config_Size, ' ', 1,
	2                               Config_Lun, Config_Index,
	3                               Config_Status, Config_Ref_Count )
	      If (.NOT. Rstat) Then
	         FUT_FCS_Reference = 1
	         Write (6,*) 'Error - opening reference file ', Ref_File,
	1           '  Status = ', Rstat
	         If (Report) Write (Lun_Rpt,*)
	1           'Error - opening reference file ', Ref_File,
	2           '  Status = ', Rstat
	      Else
	         If (Report) Then
	            Write (Lun_Rpt, *) 'Opened Reference file ', Ref_file
	         EndIf

	         Rstat = CCT_Get_Config_TOD ( Access_Time, 1, Config_Size,
	1                                     Config_Lun, Config_Index,
	2                                     Config_Record,
	3                                     Config_New_Segment,
	4                                     Config_Status )
	         If (.NOT. Rstat) Then
	            FUT_FCS_Reference = 1
	            Write (6,*) 'Error - reading reference file ', Ref_File,
	1              '  Status = ', Rstat
	            If (Report) Write (Lun_Rpt,*)
	1              'Error - reading reference file ', Ref_File,
	2              '  Status = ', Rstat
	         Else
	            If (Report) Then
	               Write (Lun_Rpt, *) 'Read Reference file ', Ref_file
	            EndIf

C Determine the index of the slope and intercept for the channel scan mode

	            If ( Chan .EQ. 'RH' ) Then
	               Cnum = 1
	            Else If ( Chan .EQ. 'RL' ) Then
	               Cnum = 2
	            Else If ( Chan .EQ. 'LH' ) Then
	               Cnum = 3
	            Else If ( Chan .EQ. 'LL' ) Then
	               Cnum = 4
	            Else
	               FUT_FCS_Reference = 1
	               Write (6,*) 'Error - Channel must be RH, RL, LH, or ',
	1                 'LL.  Not ', Chan
	               If (Report) Write (Lun_Rpt,*)
	1                 'Error - Channel must be RH, RL, LH, or LL.  Not ',
	2                 Chan
	            EndIf
	            If ( Scan .EQ. 'SS' ) Then
	               Snum = 1
	            Else If ( Scan .EQ. 'SF' ) Then
	               Snum = 2
	            Else If ( Scan .EQ. 'LF'.OR. Scan .EQ. 'FL' ) Then
	               Snum = 3
	            Else
	               FUT_FCS_Reference = 1
	               Write (6,*) 'Error - Scan mode must be SS, SF, LF, ',
	1                 'or FL. Not ', Scan
	               If (Report) Write (Lun_Rpt,*)
	1                 'Error - Scan mode must be SS, SF, LF, or FL.  Not ',
	2                 Scan
	            EndIf
	            If ( FUT_FCS_Reference .EQ. 0 ) Then
	               Index = 3 * (Cnum - 1) + Snum
	               G_Slope = Config_Record.GLTCH_REC.Slope (Index)
	               G_Intercept = Config_Record.GLTCH_REC.Intercept (Index)
	            EndIf

	            Rstat = CCT_Close_Config ( 1, Config_Lun, Config_Index )
	            If (.NOT. Rstat) Then
	               FUT_FCS_Reference = 1
	               Write (6,*) 'Error - Closing reference file ', Ref_File,
	1                 '  Status = ', Rstat
	               If (Report) Write (Lun_Rpt,*) 'Error - Closing ',
	1                 'reference file ', Ref_File, '  Status = ', Rstat
	            Else
	               If (Report) Then
	                  Write (Lun_Rpt, *) 'Closed Reference file ', Ref_File
	               EndIf
	               Call Lib$Free_Lun (Lun)
	            EndIf
	         EndIf
	      EndIf
	   EndIf
	EndIf
	Return
	End
