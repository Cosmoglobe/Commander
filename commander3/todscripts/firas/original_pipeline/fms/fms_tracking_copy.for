	Integer*4  Function  FMS_TRACKING_COPY  ( Copy_Num, Out_Rec, In_Recs )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                              FMS_TRACKING_COPY.FOR
C
C  Purpose:  Copy tracking information from input maps to output maps.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Passed Parameters

	Dictionary	'FMS_SKY'
	Record /FMS_SKY/	In_Recs (2)	! 1 record from each skymap
	Record /FMS_SKY/	Out_Rec		! Output record for pixel
	Integer*2		Copy_Num	! Level number to copy

C  External

	External	FMS_Normal

C  Local

	Integer*2	Map			! Input skymap number
	Integer*2	in, out			! Counters

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Tracking_Copy = %Loc (FMS_Normal)

	If ( Copy_Num .EQ. 0 ) Then
	   out = 1
	   Do Map = 1, 2
	      in = 1
	      Do While (
	1           In_Recs(Map).Spec_Data.Combinations.Order_0_Input(in)(1:1)
	2           .GE. 'A' .AND.
	3           In_Recs(Map).Spec_Data.Combinations.Order_0_Input(in)(1:1)
	4           .LE. 'Z')
	         Out_Rec.Spec_Data.Combinations.Order_0_Input (out) =
	1           In_Recs (Map).Spec_Data.Combinations.Order_0_Input (in)
	         Out_Rec.Spec_Data.Combinations.Order_0_Nifgs (out) =
	1           In_Recs (Map).Spec_Data.Combinations.Order_0_Nifgs (in)
	         in = in + 1
	         out = out + 1
	      EndDo
	   EndDo
	ElseIf ( Copy_Num .EQ. 1 ) Then
	   out = 1
	   Do Map = 1, 2
	      in = 1
	      Do While (
	1           In_Recs(Map).Spec_Data.Combinations.Order_1_Input(in)(1:1)
	2           .GE. 'A' .AND.
	3           In_Recs(Map).Spec_Data.Combinations.Order_1_Input(in)(1:1)
	4           .LE. 'Z')
	         Out_Rec.Spec_Data.Combinations.Order_1_Input (out) =
	1           In_Recs (Map).Spec_Data.Combinations.Order_1_Input (in)
	         Out_Rec.Spec_Data.Combinations.Order_1_Nifgs (out) =
	1           In_Recs (Map).Spec_Data.Combinations.Order_1_Nifgs (in)
	         in = in + 1
	         out = out + 1
	      EndDo
	   EndDo
	Else
	   out = 1
	   Do Map = 1, 2
	      in = 1
	      Do While (
	1           In_Recs(Map).Spec_Data.Combinations.Order_2_Input(in)(1:1)
	2           .GE. 'A' .AND.
	3           In_Recs(Map).Spec_Data.Combinations.Order_2_Input(in)(1:1)
	4           .LE. 'Z')
	         Out_Rec.Spec_Data.Combinations.Order_2_Input (out) =
	1           In_Recs (Map).Spec_Data.Combinations.Order_2_Input (in)
	         Out_Rec.Spec_Data.Combinations.Order_2_Nifgs (out) =
	1           In_Recs (Map).Spec_Data.Combinations.Order_2_Nifgs (in)
	         in = in + 1
	         out = out + 1
	      EndDo
	   EndDo
	EndIf
	Return
	End
