	Integer*4  Function  FMS_CLOSE  ( Lun_In, Lun_Out, Report, Lun_Rep,
	1                                 Numpix, In_Skymap, Specif_Out,
	2                                 File_Ext )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                            FMS_CLOSE.FOR
C
C  Purpose:  Close input and output skymap files, and report file.
C
C  Author:  Larry P. Rosen, Hughes STX, February 1993.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Include

	Include		'(FUT_Params)'
	Include		'(FUT_Error)'
	Include		'(FMS_MSG)'

C  Passed Parameters

	Integer*4	Lun_Out, Lun_In (2)	! Output and input logical unit
	Logical*1	Report			! Whether or not to report
	Integer*4	Lun_Rep			! Report logical unit #
	Integer*2	Numpix (3)		! Number of pixels with data
						! 1=map1, 2=map2, 3=merged map
	Character*33	In_Skymap (2)		! Input skymaps
	Character*4	Specif_Out		! Output file chanscan specify
	Character*20	File_Ext		! Output file extension

C  External

	External	CSA_Normal

C  Function

	Integer*4	CSA_Close_Skymap

C  Local
	Integer*4	Rstat
	Character*33	Outfile		! Output file name sans archive.

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Close = %Loc (FMS_Normal)

	If (Report) Then
	   Write (Lun_Rep,10) Numpix
   10	   Format ( 2X, 'The number of pixels with data in skymap 1 = ', I4,/,
	1           2X, 'The number of pixels with data in skymap 2 = ', I4,/,
	2           2X, 'The number of pixels with data in output skymap = ',
	3           I4, / )
	Else
	   Write (6,10) Numpix
	EndIf
	Rstat = CSA_Close_Skymap ( Lun_In (1), FAC_Skymap_No_Levels )
	If (Rstat .NE. %Loc (CSA_Normal)) Then
	   FMS_Close = %Loc (FMS_Abort)
	   Call Lib$Signal (FMS_CSAClose, %Val(2), In_Skymap(1), %Val(Rstat))
	Else
	   Rstat = CSA_Close_Skymap ( Lun_In (2), FAC_Skymap_No_Levels )
	   If (Rstat .NE. %Loc (CSA_Normal)) Then
	      FMS_Close = %Loc (FMS_Abort)
	      Call Lib$Signal (FMS_CSAClose, %Val(2), In_Skymap(2),%Val(Rstat))
	   Else
	      Rstat = CSA_Close_Skymap ( Lun_Out, FAC_Skymap_No_Levels )
	      If (Rstat .NE. %Loc (CSA_Normal)) Then
	         FMS_Close = %Loc (FMS_Abort)
	         Outfile = 'FMS_SKY_' // Specif_Out // '.' // File_Ext
	         Call Lib$Signal ( FMS_CSAClose, %Val(2), Outfile, %Val(Rstat))
	      EndIf
	   EndIf
	EndIf

C  Signal the processing status.

	If ( FMS_Close .EQ. %Loc (FMS_Normal) ) Then
	   Call Lib$Signal (FMS_Normal)
	Else
	   Call Lib$Signal (FMS_Abort)
	EndIf

C  Close report

	If (Report) Then
	   Close ( Unit=Lun_Rep, Iostat=Rstat )
	   FUT_Report_Lun = 0
	   If (Rstat .NE. 0) Then
	      FMS_Close = %Loc (FMS_Abort)
	      Call Lib$Signal ( FMS_Closerr, %Val(2), 'Report file',
	1                       %Val(Rstat) )
	   EndIf
	EndIf
	Return
	End
