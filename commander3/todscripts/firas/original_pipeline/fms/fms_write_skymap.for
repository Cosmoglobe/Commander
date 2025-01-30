	Integer*4  Function  FMS_WRITE_SKYMAP  ( Lun_Out, Out_Rec )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                            FMS_WRITE_SKYMAP.for
C
C  Purpose:  Write the record to the output skymap.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Include

	Include		'(FMS_MSG)'		! FMS messages

C  Passed parameters

	Integer*4		Lun_Out		! Output logical unit number
	Dictionary		'FMS_SKY'
	Record /FMS_SKY/	Out_Rec		! Output record for pixel

C  External

	External	CSA_Normal

C  Function

	Integer*4	CSA_Write_Pixels

C  Local

	Integer*4	Rstat			! Return status

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Write_Skymap = %Loc (FMS_Normal)

	Rstat = CSA_Write_Pixels ( Lun_Out, Out_Rec, 1, )

	If ( Rstat .NE. %Loc (CSA_Normal) ) Then
	   FMS_Write_Skymap = %Loc (FMS_Abort)
	   Call Lib$Signal ( FMS_CSAWrite, %Val (1), %Val (Rstat))
	EndIf
	Return
	End
