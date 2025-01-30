	Integer*4  Function  FMS_READ_SKYMAPS  ( Lun_In, Pixel, In_Recs,
	1                                        Rec_Found )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                             FMS_READ_SKYMAPS.FOR
C
C  Purpose:  Read 1 record from each input skymap for the given pixel.
C            Rec_Found will flag each input record as true if a record
C            was read successfully.  It is possible that neither skymap
C            has a record in this pixel.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Include files.

	Include		'(CSA_Pixel_Input_Rec)'
	Include		'(CSA_Pixel_Output_Rec)'
	Include		'(FUT_Params)'
	Include		'(FMS_Msg)'

C  Passed parameters.

	Integer*4	Lun_In (2)	! Logical unit numbers of open skymaps
	Integer*4	Pixel		! FIRAS pixel number
	Dictionary	'FMS_SKY'
	Record /FMS_SKY/	In_Recs (2)	! 1 record from each skymap
	Logical*1	Rec_Found (2)	! Flag record found in maps 1,2

C  External

	External	CSA_Normal
	
C  Functions:

	Integer*4	CSA_Read_Pixels

C  Local
	Integer*4	Rstat
	Record /Pixel_Input_List/	Inlist
	Record /Pixel_Output_List/	Outlist
	Integer*4	Dummy				! Junk data
	Integer*2	Map

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Read_Skymaps = %Loc (FMS_Normal)

	Inlist.Level_No = FAC_Skymap_Level
	Inlist.Pixel_No = Pixel
	Do Map = 1, 2
	   Rstat = CSA_Read_Pixels ( Lun_In (Map), Inlist, 1, In_Recs (Map),
	1                            1, Outlist, Dummy, 0 )

	   If ( Rstat .NE. %Loc (CSA_Normal) ) Then
	      FMS_Read_Skymaps = %Loc (FMS_Abort)
	      Call Lib$Signal ( FMS_CSARead, %Val(1), %Val (Rstat) )
	   Else
	      If (Outlist.No_Records .EQ. 1) Then
	         Rec_Found (Map) = .TRUE.
	      Else
	         Rec_Found (Map) = .FALSE.
	      EndIf
	   EndIf
	EndDo
	Return
	End
