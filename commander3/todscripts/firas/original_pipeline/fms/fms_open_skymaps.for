	Integer*4  Function  FMS_OPEN_SKYMAPS  ( Lun_In, Lun_Out, Arch_In,
	1                                        Arch_Out, In_Skymap,
	2                                        Specif_Out, File_Ext, Report,
	3                                        Lun_Rep )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                            FMS_OPEN_SKYMAPS.FOR
C
C  Purpose:  Open input and output skymap files.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Include files.

	Include		'($SSDef)'
	Include		'(FUT_Params)'
	Include		'(FMS_MSG)'

C  Passed parameters.

	Integer*4	Lun_In (2)			! Input logical units
	Integer*4	Lun_Out				! Output logical unit
	Character*14	Arch_In				! Input data archive
	Character*15	Arch_Out			! Output data archive
	Character*33	In_Skymap (2)			! Input skymaps
	Character*4	Specif_Out	! "chanscan" specification of output
	Character*20	File_Ext		! Output file extension
	Logical*1	Report			! Whether or not to report
	Integer*4	Lun_Rep			! Report logical unit #

C  Function

	Integer*4	Lib$Get_Lun
	Integer*4	CSA_Field_Offset_Values

C  External

	External	CSA_Normal
	External	CSA_Open_Skymap
	External	CSA_Field_Offset_Values

C  Local
	Character*47	In_File (2)
	Character*48	Out_File
	Integer*4	Rstat
	Integer*2	Skymaplen		! Skymap length in longwords
	Integer*4	F_Pix_Offset	! I*4 of fac_coad_spec_pix_offset
	Integer*2	Map			! Input map number
	Integer*2	Flen			! Length of filename string

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Open_Skymaps = %Loc (FMS_Normal)

C  Construct the input and output file names.

	In_File (1) = Arch_In // In_Skymap (1)
	In_File (2) = Arch_In // In_Skymap (2)
	Out_File = Arch_Out // 'FMS_SKY_' // Specif_Out // '.' // File_Ext

C  Get lun's.

	Rstat = Lib$Get_Lun (Lun_In (1))
	If (Rstat .NE. SS$_Normal) Then
	   FMS_Open_Skymaps = %Loc (FMS_Abort)
	   Call Lib$Signal (FMS_Lunerr, %Val(1), %Val(Rstat))
	Else
	   Rstat = Lib$Get_Lun (Lun_In (2))
	   If (Rstat .NE. SS$_Normal) Then
	      FMS_Open_Skymaps = %Loc (FMS_Abort)
	      Call Lib$Signal (FMS_Lunerr, %Val(1), %Val(Rstat))
	   Else
	      Rstat = Lib$Get_Lun (Lun_Out)
	      If (Rstat .NE. SS$_Normal) Then
	         FMS_Open_Skymaps = %Loc (FMS_Abort)
	         Call Lib$Signal (FMS_Lunerr, %Val(1), %Val(Rstat))
	      EndIf
	   EndIf
	EndIf

C  Open the output skymap.

	If (FMS_Open_Skymaps .EQ. %Loc (FMS_Normal)) Then
	   Skymaplen = FAC_Coad_Spec_Size / 4
	   F_Pix_Offset = FAC_Coad_Spec_Pix_Offset
	   Call Str$Trim ( Out_File, Out_File, Flen )
	   Open ( Unit=Lun_Out, File=Out_File (1:Flen), Status='New',
	1         Form='Unformatted', Recordtype='Fixed', Recl=Skymaplen,
	2         Useropen=CSA_Open_Skymap, Iostat=Rstat )
	   If (Rstat .NE. 0) Then
	      FMS_Open_Skymaps = %Loc (FMS_Abort)
	      Call Lib$Signal (FMS_CSAOpen, %Val(2), Out_File, %Val(Rstat))
	   Else
	      Rstat = CSA_Field_Offset_Values ( F_Pix_Offset, FAC_Time_Offset,
	1                                       -1, Lun_Out )
	      If (Rstat .NE. %Loc (CSA_Normal)) Then
	         FMS_Open_Skymaps = %Loc (FMS_Abort)
	         Call Lib$Signal (FMS_CSAOffset, %Val(1), %Val(Rstat))
	      EndIf
	   EndIf
	EndIf

C  Open the input skymaps.

	Do Map = 1, 2
	   If (FMS_Open_Skymaps .EQ. %Loc (FMS_Normal)) Then
	      Call Str$Trim ( In_File (Map), In_File (Map), Flen )
	      Open ( Unit=Lun_In (Map), File=In_File (Map) (1:Flen),
	1            Status='Old', Form='Unformatted', Recordtype='Fixed',
	2            Recl=Skymaplen, Readonly, Useropen=CSA_Open_Skymap,
	3            Iostat=Rstat )
	      If (Rstat .NE. 0) Then
	         FMS_Open_Skymaps = %Loc (FMS_Abort)
	         Call Lib$Signal(FMS_CSAOpen,%Val(2),In_File (Map),%Val(Rstat))
	      Else
	         Rstat = CSA_Field_Offset_Values ( F_Pix_Offset,
	1                                          FAC_Time_Offset, -1,
	2                                          Lun_In (Map) )
	         If (Rstat .NE. %Loc (CSA_Normal)) Then
	            FMS_Open_Skymaps = %Loc (FMS_Abort)
	            Call Lib$Signal (FMS_CSAOffset, %Val(1), %Val(Rstat))
	         EndIf
	      EndIf
	   EndIf
	EndDo
	Return
	End
