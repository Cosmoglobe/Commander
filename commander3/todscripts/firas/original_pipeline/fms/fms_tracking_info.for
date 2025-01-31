	Integer*4  Function FMS_TRACKING_INFO  ( In_Recs, Out_Rec, Specif,
	1                                        Nifgs, Order_Num, In_Skymap )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                              FMS_TRACKING_INFO.FOR
C
C  Purpose:  Fill in the tracking information for the merge in the COMBINATIONS
C            substructure, of the output record.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Passed Parameters

	Dictionary	'FMS_SKY'
	Record /FMS_SKY/	In_Recs (2)	! 1 record from each skymap
	Record /FMS_SKY/	Out_Rec		! Output record for pixel
	Character*4		Specif (2)	! "chanscan" of input skymaps
	Real*8			Nifgs (2)	! Floating point number of ifgs
	Integer*2		Order_Num	! Current order of combination
	Character*33		In_Skymap (2)	! Input skymaps

C  External

	External	FMS_Normal

C  Function

	Integer*4	FMS_TRACKING_COPY

C  Local

	Integer*2	INifgs (2)		! Integer version of Nifgs * 10
	Integer*4	Rstat			! Return status
	Integer*2	Map			! Input skymap number
	Integer*2	in, out			! Counters
	Character*4	Empty			! Empty character string

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Tracking_Info = %Loc (FMS_Normal)

	INifgs (1) = Nint (10. * Nifgs (1)) 
	INifgs (2) = Nint (10. * Nifgs (2))
	Out_Rec.Spec_Data.Combinations.Current_Order = Order_Num

	If ( Order_Num .EQ. 0 ) Then
	   Out_Rec.Spec_Data.Combinations.Order_0_Input (1) = Specif (1)
	   Out_Rec.Spec_Data.Combinations.Order_0_Input (2) = Specif (2)
	   Out_Rec.Spec_Data.Combinations.Order_0_Nifgs (1) = INifgs (1)
	   Out_Rec.Spec_Data.Combinations.Order_0_Nifgs (2) = INifgs (2)
	ElseIf ( Order_Num .EQ. 1 ) Then
	   Out_Rec.Spec_Data.Combinations.Order_1_Input (1) = Specif (1)
	   Out_Rec.Spec_Data.Combinations.Order_1_Input (2) = Specif (2)
	   Out_Rec.Spec_Data.Combinations.Order_1_Nifgs (1) = INifgs (1)
	   Out_Rec.Spec_Data.Combinations.Order_1_Nifgs (2) = INifgs (2)
	   Rstat = FMS_TRACKING_COPY (0, Out_Rec, In_Recs)
	ElseIf ( Order_Num .EQ. 2 ) Then
	   Out_Rec.Spec_Data.Combinations.Order_2_Input (1) = Specif (1)
	   Out_Rec.Spec_Data.Combinations.Order_2_Input (2) = Specif (2)
	   Out_Rec.Spec_Data.Combinations.Order_2_Nifgs (1) = INifgs (1)
	   Out_Rec.Spec_Data.Combinations.Order_2_Nifgs (2) = INifgs (2)
	   Rstat = FMS_TRACKING_COPY (0, Out_Rec, In_Recs)
	   Rstat = FMS_TRACKING_COPY (1, Out_Rec, In_Recs)
	Else
	   Out_Rec.Spec_Data.Combinations.Order_3_Input (1) = Specif (1)
	   Out_Rec.Spec_Data.Combinations.Order_3_Input (2) = Specif (2)
	   Out_Rec.Spec_Data.Combinations.Order_3_Nifgs (1) = INifgs (1)
	   Out_Rec.Spec_Data.Combinations.Order_3_Nifgs (2) = INifgs (2)
	   Rstat = FMS_TRACKING_COPY (0, Out_Rec, In_Recs)
	   Rstat = FMS_TRACKING_COPY (1, Out_Rec, In_Recs)
	   Rstat = FMS_TRACKING_COPY (2, Out_Rec, In_Recs)
	EndIf
	Return
	End
