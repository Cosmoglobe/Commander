	Integer*4  Function  FMS_COMBINE_OTHER  ( Rec_Found, In_Recs, Out_Rec,
	1                                         Specif, Specif_Out, Nifgs,
	2                                         No_Freq_Weights, Numpix )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                          FMS_COMBINE_OTHER.FOR
C
C  Purpose:  Combine all fields other than spectra and variances.  Forms sums
C            and computes weighted averages.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Include file.

	Include		'(FUT_FCS_Include)'
	Include		'(FUT_Params)'

C  Passed parameters.

	Logical*1	Rec_Found (2)		! Flag record found in maps 1,2
	Dictionary	'FMS_SKY'
	Record /FMS_SKY/	In_Recs (2)	! 1 record from each skymap
	Record /FMS_SKY/	Out_Rec		! Output record for pixel
	Character*4	Specif (2)		! "chanscan" of input skymaps
	Character*4	Specif_Out	! "chanscan" specification of output
	Real*8		Nifgs (2)		! Floating point number of ifgs
	Real*8		No_Freq_Weights (2)	! Freq indep. weight of skymaps
	Integer*2	Numpix (3)		! Number of pixels with data
						! 1=map1, 2=map2, 3=merged map

C  External

	External	FMS_Normal

C  Function

	Integer*4	FUT_FCS_NON_AVG
	Integer*4	FUT_FCS_MINMAX
	Integer*4	FUT_FCS_SUM
	Integer*4	FUT_FCS_AVG

C  Local

	Integer*2	Map			! Input skymap number
	Logical*1	First			! First map
	Character*3	Data /'FMS'/		! Data type.
	Real*8		Weight			! Weight for averaging
	Real*8		Total			! Total number of ifgs combined
	Real*8		Sum_Weight		! Total of weights
	Real*4		Moon_Phase, Orb_Phase, Scan_Angle	! some angles
	Record / Double_FCF / Sum_Rec, Rec1	! Record to store sums.
	Byte		Blank_Sum (Sum_Rec_Size) / Sum_Rec_Size * 0 /
	Integer*2	NWeight			! NINT of Weight, for angles.
	Integer*2	N			! Counter
	Integer*4	Anum			! Counter for each map
	Integer*2	Numrecs			! Number of records used
	Integer*4	Maxsum			! Max # for sum of nifgs
	Parameter	(Maxsum = 10000)	! Must be large to accomodate
					! high order merges of many ifgs.
	Real*4		Mphase (Maxsum), Ophase (Maxsum), Sangle (Maxsum)
	Integer*4	Time_Array (2, Maxsum)	! Time array to average
	Integer*4	Time_Weights (Maxsum)
	Character*1	Xcal_Pos
	Integer*4	Rstat			! Return status

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Combine_Other = %Loc (FMS_Normal)

C  Initialize 1st data record and sum record to zeros.

	Call Lib$MovC3 (Sum_Rec_Size, Blank_Sum, Rec1)

	Call Lib$MovC3 (Sum_Rec_Size, Blank_Sum, Sum_Rec)

C  Initialize some counters

	Anum = 0
	Numrecs = 0
	Total = 0.
	Sum_Weight = 0.
	First = .TRUE.

C  Record Loop: do for each input record, maximum of one from each map.

	Do Map = 1, 2

	   If (Rec_Found (Map) .AND. Nifgs (Map) .GT. 0.) Then

	      Numrecs = Numrecs + 1
	      Numpix (Map) = Numpix (Map) + 1

C  If first record, call FUT_FCS_NON_AVG to store non-averaged data in output
C  record.

	      If (First) Then
	         Rstat = FUT_FCS_NON_AVG ( Out_Rec, In_Recs (Map),
	1                                  Specif (1)(3:4), Specif (2)(3:4),
	2                                  Specif (1)(1:2), Specif (2)(1:2),
	3                                  Data )
	      EndIf

C  Increase total of floating point number of ifgs for this pixel.
C  Increase the sum of the spectrum weights for this pixel.
C  If number of ifgs is less than or equal to 1, we can't compute a sigma.

	      Total = Total + Nifgs (Map)
	      Weight = No_Freq_Weights (Map) * Nifgs (Map)
	      Sum_Weight = Sum_Weight + Weight

C  Store min and max of time and engineering quantities

	      Rstat = FUT_FCS_MINMAX ( First, Out_Rec, In_Recs (Map) )

C  Form sums of engineering quantities and attitude.

	      Rstat = FUT_FCS_SUM ( Weight, Out_Rec, In_Recs (Map),
	1                           Moon_Phase, Orb_Phase, Scan_Angle,
	2                           Sum_Rec, Rec1, First )

C  Store the attitude angles for averaging.
C  Note:  I have to use the record's integer number of ifgs for weighting the
C  angles and times because there are cases where the weights were less than
C  0.5 and the nearest integers were 0.

	      NWeight = Max ( In_Recs (Map).Coad_Spec_Head.Num_Ifgs, 1 )
	      Do N = 1, NWeight
	         Mphase (Anum+N) = Moon_Phase
	         Ophase (Anum+N) = Orb_Phase
	         Sangle (Anum+N) = Scan_Angle
	      EndDo
	      Time_Array (1,Numrecs) = In_recs (Map).CT_Head.Time (1)
	      Time_Array (2,Numrecs) = In_recs (Map).CT_Head.Time (2)
	      Time_Weights (Numrecs) = NWeight
	      Anum = Anum + NWeight
	      First = .FALSE.
	   EndIf					! If Nifgs GT 0
	EndDo						! Do for each map

	Out_Rec.Coad_Spec_Head.Num_Ifgs =
	1   Nint (Out_Rec.Coad_Spec_Head.Comb_Num_Ifgs)

	If (Total .GT. 0.) Then
	   Numpix (3) = Numpix (3) + 1

C  Call FUT_FCS_Avg to average engineering and attitude quantities.
C  Note that Anum is equal to the total number of Mphase, Ophase and Sangle.

	   Rstat = FUT_FCS_AVG  ( Maxsum, Time_Array, Numrecs, Time_Weights,
	1                         Out_Rec, Sum_Rec, Rec1, Sum_Weight,
	2                         Mphase, Ophase, Sangle, Anum )

c Use Encode to create text label.

	   If (Out_Rec.Coad_Spec_Data.Xcal_Pos .EQ. FAC_XcalOut) Then
	      Xcal_Pos = 'O'
	   Else
	      Xcal_Pos = 'I'
	   EndIf

	   If (Out_Rec.Attitude.Pixel_No .GE. 0) Then
	      Encode ( 60, 10, Out_Rec.Coad_Spec_Head.Label)
	1              Out_Rec.Coad_Spec_Head.First_GMT,
	2              Out_Rec.Coad_Spec_Head.Last_GMT,
	3              Out_Rec.Coad_Spec_Head.Num_IFGs,
	4              Specif_Out,
	5              Nint (Out_Rec.Coad_Spec_Data.Ical * 1000.0),
	6              Xcal_Pos, Out_Rec.Attitude.Pixel_No
  10	      Format (2(A,'_'),I3.3,'_',A,'_',I5.5,'_',A,'_',I4.4,9X)
	   Else
	      Encode ( 60, 10, Out_Rec.Coad_Spec_Head.Label)
	1              Out_Rec.Coad_Spec_Head.First_GMT,
	2              Out_Rec.Coad_Spec_Head.Last_GMT,
	3              Out_Rec.Coad_Spec_Head.Num_IFGs,
	4              Specif_Out,
	5              Nint (Out_Rec.Coad_Spec_Data.Ical * 1000.0),
	6              Xcal_Pos
  20	      Format (2(A,'_'),I3.3,'_',A,'_',I5.5,'_',A,14X)
	   EndIf
	EndIf
	Return
	End
