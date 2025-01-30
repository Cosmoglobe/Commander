	Integer*4  Function  FMS_READ_C_D_VECTORS  ( Model_Ext, In_Skymap,
	1                                            Specif, Report, Lun_Rep,
	2                                            Arch_Ref, Arch_In,
	3                                            C_Vector_Rec,
	4                                            D_Vector_Rec,
	5                                            Access_Time )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                            FMS_READ_C_D_VECTORS.FOR
C
C  Purpose:  Open, read, and close C and D vector files.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C  Modified:  10 / 94, LPR, moved Access_Time into a passed parameter so that
C     it can be used in the FMS_CORRECTION_SPECTRUM module.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Include files.

	Include		'($SSDef)'
	Include		'(FMS_Msg)'
	Include		'(CCT_Get_Config)'
	Include		'CSDR$Library:CTParams.inc'
	Include		'(CCT_Query_Catalog_Record)'

C  Passed parameters.

	Character*12	Model_Ext		! Cal model soln extension
	Character*33	In_Skymap (2)		! Input skymaps
	Character*4	Specif (2)		! "chanscan" of input skymaps
	Logical*1	Report			! Whether or not to report
	Integer*4	Lun_Rep			! Report logical unit #
	Character*15	Arch_Ref		! C or D Vector archive
	Character*14	Arch_In			! Input data archive
	Dictionary	'FEX_CVS'
	Record /FEX_CVS/	C_Vector_Rec (2)
	Dictionary	'FEX_VAR'
	Record /FEX_VAR/	D_Vector_Rec (2)
	Integer*4	Access_Time (2)		! time for reference file

C  Functions

	Integer*4	Lib$Get_Lun
	Integer*4	CCT_Query_Catalog
	Integer*4	CCT_Open_Config
	integer*4	CCT_Get_Config_TOD
	Integer*4	CCT_Close_Config

C  Local

	Integer*2	Map				! Input map number
	Character*48	Cfile, Dfile			! C & D vector files.
	Integer*4	Rstat				! Return status
	Integer*4	C_Lun, D_Lun			! Logical Unit Numbers
	Dictionary 'ccm_cme_catalog_entry'
	Record /query_catalog/		Query_Cat
	Record /CCM_CME_Catalog_Entry/	Cats (5)

C  Config stuff

	Character*14	Config_GMT_Start / '86001000001000' /
	Character*14	Config_GMT_Stop / '99365235958990' /
	Integer*4	Config_Start (2), Config_Stop (2)
	Character*48	Config_Name (2)
	Integer*4	Config_Size (2) / 2176, 2176 /
	Integer*4	Config_Lun (2), Config_Index (2)
	Logical*1	Config_New_Segment (2)
	Record /Config_Status/ Config_Status(2)
	Integer*4	Config_Ref_Count
	Structure /Config_Record/
	   record /FEX_CVS/ CVEC
	   record /FEX_VAR/ DVEC
	EndStructure
	Record 	/Config_Record/ Config_Record

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Read_C_D_Vectors = %Loc (FMS_Normal)

c Query catalog for time of input skymap for use in getting reference files.

	Query_Cat.Archive_ID = Arch_In
	Query_Cat.Filename = In_Skymap (1)
	Rstat = CCT_Query_Catalog ( Query_Cat, Cats(1) )
	If (.NOT. Rstat) Then
	   FMS_Read_C_D_Vectors = %Loc (FMS_Abort)
	   Call Lib$Signal ( fms_query_cat, %Val (1), %Val (rstat) )
	Else
	   Access_Time (1) = Cats(1).Initial_Time (1)
	   Access_Time (2) = Cats(1).Initial_Time (2)
	EndIf

C  Determine the C-Vector and D-Vector file names.  If input is FCS then use
C  FEX_VAR and FEX_CVS.  If input is FMS then use FMS_VAR and FMS_CVS.

	Map = 1					! Do for each input map
	Do While (Map .LE. 2 .AND. FMS_Read_C_D_Vectors .EQ. %Loc (FMS_Normal))

	   If (In_Skymap (Map)(1:3) .EQ. 'FCS') Then
	      Dfile = Arch_Ref // 'FEX_VAR_' // Specif (Map) // '.'// Model_Ext
	      Cfile = Arch_Ref // 'FEX_CVS_' // Specif (Map) // '.'// Model_Ext
	   Else
	      Dfile = Arch_In // 'FMS_VAR_' // Specif (Map) // '.'// Model_Ext
	      Cfile = Arch_In // 'FMS_CVS_' // Specif (Map) // '.'// Model_Ext
	   EndIf

C  Get logical unit numbers, open, read, and close c and d vectors for each
C  input skymap.

	   Rstat = Lib$Get_Lun (C_Lun)
	   If (Rstat .NE. SS$_Normal) Then
	      FMS_Read_C_D_Vectors = %Loc (FMS_Abort)
	      Call Lib$Signal (FMS_Lunerr, %Val(1), %Val(Rstat))
	   Else
	      Rstat = Lib$Get_Lun (D_Lun)
	      If (Rstat .NE. SS$_Normal) Then
	         FMS_Read_C_D_Vectors = %Loc (FMS_Abort)
	         Call Lib$Signal (FMS_Lunerr, %Val(1), %Val(Rstat))
	      Else
	      EndIf
	   EndIf
	   If (FMS_Read_C_D_Vectors .EQ. %Loc (FMS_Normal)) Then
	      Call CT_GMT_To_Binary (Config_GMT_Start, Config_Start)
	      Call CT_GMT_To_Binary (Config_GMT_Stop, Config_Stop)
	      Config_Name (1) = Cfile
	      Config_Name (2) = Dfile
	      Config_Lun (1) = C_Lun
	      Config_Lun (2) = D_Lun
	      Rstat = CCT_Open_Config ( Config_Start, Config_Stop, 2,
	1                               Config_Name, Config_Size, ' ', 1,
	2                               Config_Lun, Config_Index,
	3                               Config_Status, Config_Ref_Count )
	      If (.NOT. Rstat) Then
	         FMS_Read_C_D_Vectors = %Loc (FMS_Abort)
	         Call Lib$Signal (FMS_Openerr, %Val(2), Cfile, %Val(Rstat))
	         Call Lib$Signal (FMS_Openerr, %Val(2), Dfile, %Val(Rstat))
	      Else
	         Rstat = CCT_Get_Config_TOD ( Access_Time, 2, Config_Size,
	1                                     Config_Lun, Config_Index,
	2                                     Config_Record,
	3                                     Config_New_Segment,
	4                                     Config_Status )
	         If (.NOT. Rstat) Then
	            FMS_Read_C_D_Vectors = %Loc (FMS_Abort)
	            Call Lib$Signal (FMS_Readerr, %Val(2), Cfile, %Val(Rstat))
	            Call Lib$Signal (FMS_Readerr, %Val(2), Dfile, %Val(Rstat))
	         Else
	            C_Vector_Rec (Map) = Config_Record.CVEC
	            D_Vector_Rec (Map) = Config_Record.DVEC
	            Rstat = CCT_Close_Config ( 2, Config_Lun, Config_Index )
	            If (.NOT. Rstat) Then
	               FMS_Read_C_D_Vectors = %Loc (FMS_Abort)
	               Call Lib$Signal (FMS_Closerr, %Val(2),Cfile,%Val(Rstat))
	               Call Lib$Signal (FMS_Readerr, %Val(2),Dfile,%Val(Rstat))
	            EndIf
	         EndIf
	      EndIf
	   EndIf
	   Map = Map + 1
	EndDo
	If (FMS_Read_C_D_Vectors .EQ. %Loc (FMS_Normal)) Then
	   Call Lib$Free_Lun (C_Lun)
	   Call Lib$Free_Lun (D_Lun)
	EndIf
	Return
	End
