	Integer*4  Function  FMS_COMPUTE_WEIGHTS  ( Lun_In, Trans_Freq,
	1                       Specif, Weight, C_Vector_Rec, D_Vector_Rec,
	2                       Weights, Relative_Weights, Model, Nu_Last_Low,
	3                       No_Freq_Weights, C_Vector, D_Vector,
	4                       Combo_Type )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                            FMS_COMPUTE_WEIGHTS.FOR
C
C  Purpose:  Compute relative weights for combination.
C
C  Method:
C     1st check whether to use c-vector or d-vector.
C     Do for each map.
C        Other = other map number (ie, if map=1 then other=2 and vice versa).
C        Extract the c and d vectors from their records.
C        Set the Weights vector equal to c-vector or d-vector.
C        Compute relative weights. If either denominator is 0, set weight to 0.
C        Compute sum of relative weights for calculation of frequency
C           independent weight.
C        Extract the number of calibration ifgs.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C  Modified: April 1994, L. Rosen.
C     Frequency independent weight is supposed to be averaged over frequencies
C     2 - 21 cm-1 for low freq, and 2 - 97 cm-1 for high freq, whereas I
C     originally did average over all 257.  This correction means that I need
C     to read skymaps to get the nyquist frequency, so that I can check the
C     frequencies against the ranges.  Frequency of bin i is Nyq * (i-1) / 256.
C     Nu_Last_Low is the highest bin of the low frequency range.
C  Modified: 19 August 1994, L. Rosen.  Formula for Frequency Independent
C     Weight is different for zeroeth order combinations.
C  Modified: 16 November 1994, L. Rosen.  The 2-21 and 2-97 are not true for
C     high resolution maps which go down to 1 icm.  Also there is a bug in
C     the "total" when summing relative weights of value 0.  A much more
C     straightforward approach is to only sum weights when the C (or D)
C     vector is non-zero.  This will work for all data, high frequency,
C     low frequency, and high resolution.  The frequency ranges are still
C     used in zeroth order to be consistent with the code used for the
C     previous delivery (summer 1994) of the FMS zeroth order merges. 
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Passed parameters.

	Integer*4	Lun_In (2)
	Integer*2	Trans_Freq	! Transition Low-High Frequency in icm
	Character*4	Specif (2)		! "chanscan" of input skymaps
	Character*8	Weight			! C_VECTOR or D_VECTOR
	Dictionary	'FEX_CVS'
	Record /FEX_CVS/	C_Vector_Rec (2)
	Dictionary	'FEX_VAR'
	Record /FEX_VAR/	D_Vector_Rec (2)
	Real*8		Weights (257,2)		! C or D vector
	Real*8		Relative_Weights (257,2)
	Character*40	Model			! Calibration Solution
	Integer*2	Nu_Last_Low		! Bin # top of low freq range
	Real*8		No_Freq_Weights (2)	! Freq indep. weight skymap
	Real*8		C_Vector (257,2)	! C-vector for skymap 1 and 2
	Real*8		D_Vector (257,2)	! D-vector for skymap 1 and 2
	Character*1	Combo_Type		! Combination type.
			! Z=Zeroth, L=Low+Low, H=High+High, F=High+Low, M=Other

C  Include

	Include		'(FUT_Params)'		! Fac parameters
	Include		'(FMS_MSG)'
	Include		'(CSA_Pixel_Input_Rec)'
	Include		'(CSA_Pixel_Output_Rec)'

C  External

	External	CSA_Normal

C  Functions

	Integer*4	CSA_Read_Pixels

C  Local

	Integer*2	Map		! Input skymap number
	Integer*2	Other		! Other skymap number
	Integer*2	i		! Counter
	Real*8		RWeight_Sum	! Relative weight sum
	Real*8		Numer		! Numerator of sum
	Real*8		Denom		! Denominator of sum
	Logical*1	Rec_Found (2)	! Flag a record read from skymap
	Integer*4	Pixel
	Integer*4	Rstat		! Return status
	Dictionary	'FMS_SKY'
	Record /FMS_SKY/	In_Recs (2)	! 1 record from each skymap
	Real*4		Nyquist (2)	! Nyquist frequency (each map) / 256.
	Record /Pixel_Input_List/	Inlist
	Record /Pixel_Output_List/	Outlist
	Integer*4	Dummy				! Junk data for CSA
	Integer*2	Total		! Number of frequencies to avg over
	Real*4		Freq		! Frequency of bin
	Integer*2	LowF / 2 /	! Low frequency cut-off on sum (icm)
	Integer*2	HighF		! High frequency cut-off

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Compute_Weights = %Loc (FMS_Normal)

C  Loop through pixels until find first record in each skymap.  Get nyquist.

	Nu_Last_Low = 1
	Rec_Found (1) = .FALSE.
	Rec_Found (2) = .FALSE.
	Inlist.Level_No = FAC_Skymap_Level

	Do Map = 1,2
	   Pixel = 0
	   Inlist.Pixel_No = Pixel
	   Do While ( Pixel .LT. FAC_Num_Pixels .AND.
	1             (.NOT. Rec_Found (Map)) .AND.
	2             FMS_Compute_Weights .EQ. %Loc (FMS_Normal) )

	      Rstat = CSA_Read_Pixels ( Lun_In (Map), Inlist, 1, In_Recs (Map),
	1                               1, Outlist, Dummy, 0 )

	      If ( Rstat .NE. %Loc (CSA_Normal) ) Then
	         FMS_Compute_Weights = %Loc (FMS_Abort)
	         Call Lib$Signal ( FMS_CSARead, %Val(1), %Val (Rstat) )
	      Else
	         If (Outlist.No_Records .EQ. 1) Then
	            Rec_Found (Map) = .TRUE.
	         Else
	            Rec_Found (Map) = .FALSE.
	         EndIf
	      EndIf
	      Pixel = Pixel + 1
	      Inlist.Pixel_No = Pixel
	   EndDo
	   If (Rec_Found (Map)) Then
	      Nyquist (Map) = In_Recs (Map).Coad_Spec_Data.Nyquist_Icm / 256.
	   Else
	      FMS_Compute_Weights = %Loc (FMS_Abort)
	   EndIf
	EndDo						! for each map

C Calculate the last bin of the low frequencies less than the transition freq.

	Freq = Nyquist (1) * Float (Nu_Last_Low - 1)
	Do While (Freq .LT. Trans_Freq)
	   Freq = Nyquist (1) * Float (Nu_Last_Low)
	   Nu_Last_Low = Nu_Last_Low + 1
	EndDo
	Nu_Last_Low = Nu_Last_Low - 1

	If ( Combo_Type .EQ. 'Z' ) Then
	   If ( Specif (1) (2:2) .EQ. 'L' ) Then
	      HighF = Trans_Freq
	   Else
	      HighF = 97
	   EndIf
	EndIf

C  Assign weights based on C and D vectors.  Compute relative weights and
C  sum of relative weights in frequency band for frequency independent weight.

	If (Weight .EQ. 'CVECTOR') Then
	   Model = C_Vector_Rec (1).Model_Label
	   Do Map = 1,2
	      Other = Mod (Map, 2) + 1			! gives opposite map #
	      RWeight_Sum = 0.
	      Total = 0
	      Numer = 0.
	      Denom = 0.
	      Do i=1,257
	         C_Vector (i, Map) = C_Vector_Rec (Map).CVector (i)
	         D_Vector (i, Map) = D_Vector_Rec (Map).Variances (i)
	         Weights (i, Map) = C_Vector (i, Map)
	         If ( (C_Vector_Rec (Other).CVector (i) .NE. 0) .AND.
	1             (C_Vector_Rec (Other).CVector (i) .NE. -1 *
	2              C_Vector_Rec (Map).CVector (i) ) ) Then

	            Relative_Weights (i, Map) = 1. /
	1              ( 1. + ( C_Vector_Rec (Map).CVector (i) /
	2                       C_Vector_Rec (Other).CVector (i) ) )
	         Else
	            Relative_Weights (i, Map) = 0.
	         EndIf
	         If ( Combo_Type .EQ. 'Z' ) Then
	            Freq = Nyquist (Map) * (i-1)
	            If ( Freq .GT. LowF .AND. Freq .LT. HighF ) Then
	               If ( D_Vector_Rec (Map).Variances (i) .NE. 0. ) Then
	               Numer = Numer + 1.0 / D_Vector_Rec (Map).Variances (i)
	               EndIf
	               If ( D_Vector_Rec (Other).Variances (i) .NE. 0. ) Then
	               Denom = Denom + 1.0 / D_Vector_Rec (Other).Variances (i)
	               EndIf
	            EndIf
	         Else
	            If ( C_Vector_Rec (Map).CVector (i) .NE. 0. .AND.
	1                Relative_Weights (i, Map) .NE. 0. ) Then
	               RWeight_Sum = RWeight_Sum + Relative_Weights (i, Map)
	               Total = Total + 1
	            EndIf
	         EndIf
	      EndDo				! Frequency loop
	      If ( Combo_Type .EQ. 'Z' ) Then
	         If ( Denom .NE. 0. ) Then
	            No_Freq_Weights (Map) = SQRT ( Numer / Denom )
	         Else
	            No_Freq_Weights (Map) = 0.
	         EndIf
	      Else
	         If ( Total .NE. 0 ) Then
	            No_Freq_Weights (Map) = RWeight_Sum / DBLE (Total)
	         Else
	            No_Freq_Weights (Map) = 0.
	         EndIf
	      EndIf
	   EndDo					! Loop over maps
	Else			! Weight = 'D_VECTOR'
	   Model = D_Vector_Rec (1).Model_Label
	   Do Map = 1,2
	      Other = Mod (Map, 2) + 1			! gives opposite map #
	      RWeight_Sum = 0.
	      Total = 0
	      Numer = 0.
	      Denom = 0.
	      Do i=1,257
	         C_Vector (i, Map) = C_Vector_Rec (Map).CVector (i)
	         D_Vector (i, Map) = D_Vector_Rec (Map).Variances (i)
	         Weights (i, Map) = D_Vector (i, Map)
	         If ( (D_Vector_Rec (Other).Variances (i) .NE. 0) .AND.
	1             (D_Vector_Rec (Other).Variances (i) .NE. -1 *
	2              D_Vector_Rec (Map).Variances (i) ) ) Then

	            Relative_Weights (i, Map) = 1. /
	1              ( 1. + ( D_Vector_Rec (Map).Variances (i) /
	2                       D_Vector_Rec (Other).Variances (i) ) )
	         Else
	            Relative_Weights (i, Map) = 0.
	         EndIf
	         If ( Combo_Type .EQ. 'Z' ) Then
	            Freq = Nyquist (Map) * (i-1)
	            If ( Freq .GT. LowF .AND. Freq .LT. HighF ) Then
	               If ( D_Vector_Rec (Map).Variances (i) .NE. 0. ) Then
	               Numer = Numer + 1.0 / D_Vector_Rec (Map).Variances (i)
	               EndIf
	               If ( D_Vector_Rec (Other).Variances (i) .NE. 0. ) Then
	               Denom = Denom + 1.0 / D_Vector_Rec (Other).Variances (i)
	               EndIf
	            EndIf
	         Else
	            If ( D_Vector_Rec (Map).Variances (i) .NE. 0. .AND.
	1                Relative_Weights (i, Map) .NE. 0. ) Then
	               RWeight_Sum = RWeight_Sum + Relative_Weights (i, Map)
	               Total = Total + 1
	            EndIf
	         EndIf
	      EndDo
	      If ( Combo_Type .EQ. 'Z' ) Then
	         If ( Denom .NE. 0. ) Then
	            No_Freq_Weights (Map) = SQRT ( Numer / Denom )
	         Else
	            No_Freq_Weights (Map) = 0.
	         EndIf
	      Else
	         If ( Total .NE. 0 ) Then
	            No_Freq_Weights (Map) = RWeight_Sum / DBLE (Total)
	         Else
	            No_Freq_Weights (Map) = 0.
	         EndIf
	      EndIf
	   EndDo
	EndIf
	Return
	End
