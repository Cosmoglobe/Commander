 	Integer*4  Function  FMS_COMBINE_C_D_VECTORS
	1                       ( Weight, C_Vector, D_Vector, C_Vector_Rec,
	2                         D_Vector_Rec, Relative_Weights,
	3                         No_Freq_Weights, Specif_Out, Arch_Out,
	4                         Model_Ext, Report, Lun_Rep, Nu_Last_Low,
	5                         Specif, Combo_Type, Lun_In )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                            FMS_COMBINE_C_D_VECTORS.FOR
C
C  Purpose:  Combine the C and D vectors for the input skymaps and write
C            the combined vectors to new output files.  Also calculate
C            weighted floating point number of model ifgs.
C
C  Method:  D(nu) = RW1^2 * D1(nu) + RW2^2 * D2(nu) and
C           C(nu) = RW1^2 * C1(nu) + RW2^2 * C2(nu) where
C           RW1 is the relative weight for input skymap 1,...
C           N = N1 * T1 + N2 * T2 where
C           T = frequency independent weight.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C  Modification: L. Rosen, August 1994.  If zeroth order combination, then
C     just copy one of the C vectors to the output.  This is because
C     for zeroeth order the C vectors are the same in the two input maps.
C     For zeroth order, the D vectors are combined the way they are done
C     in FDS.  This uses the formula:
C               T1^2 * D1 * Sum_pix N1(p) + T2^2 * D2 * Sum_pix N2(p)
C        D =    -----------------------------------------------------
C                      T1 * Sum_pix N1(p) + T2 * Sum_pix N2(p)
C
C     Where Sum_pix N is the sum of the glitchrate weighted number of ifgs
C     over all pixels of the input skymap that are outside the galactic
C     latitude cut.  Order is zeroth if Combo_Type = 'Z'.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Include files.

	Include		'($SSDef)'
	Include		'(FMS_MSG)'
	Include		'CT$Library:CTUser.Inc'
	Include		'(FUT_Params)'		! Fac parameters

C  Passed parameters.

	Character*8	Weight			! C_VECTOR or D_VECTOR
	Real*8		C_Vector (257,2)	! C-vector for skymap 1 and 2
	Real*8		D_Vector (257,2)	! D-vector for skymap 1 and 2
	Dictionary	'FEX_CVS'
	Record /FEX_CVS/	C_Vector_Rec (2)  ! Input C-vector records
	Dictionary	'FEX_VAR'
	Record /FEX_VAR/	D_Vector_Rec (2)  ! Input D-vector records
	Real*8		Relative_Weights (257,2)
	Real*8		No_Freq_Weights (2)	! Freq indep. weight of skymaps
	Character*4	Specif_Out	! "chanscan" specification of output
	Character*15	Arch_Out		! Output data archive
	Character*12	Model_Ext		! Cal model soln extension
	Logical*1	Report			! Whether or not to report
	Integer*4	Lun_Rep			! Report logical unit #
	Integer*2	Nu_Last_Low	! Highest index for low frequency band
	Character*4	Specif (2)		! "chanscan" of input skymaps
	Character*1	Combo_Type		! Combination type.
			! Z=Zeroth, L=Low+Low, H=High+High, F=High+Low, M=Other
	Integer*4	Lun_In (2)		! Input skymaps

C  Function

	Integer*4	Lib$Get_Lun
	Integer*4	CCT_Set_TTG_Time_Range
	Integer*4	FMS_READ_SKYMAPS

C  External

	External	CT_Connect_Write

C  Local

	Character*14	Ref_Start_GMT / '89324000000000' /
	Character*14	Ref_Stop_GMT / '99365235959990' /
	Integer*4	Ref_Start_ADT (2), Ref_Stop_ADT (2)
	Record /FEX_CVS/ Out_C_Vector_Rec	! Output C-vector record
	Record /FEX_VAR/ Out_D_Vector_Rec	! Output D-vector record
	Integer*2	i			! Counter
	Integer*2	Map
	Integer*4	Rstat			! Return Status
	Integer*4	Lun_C, Lun_D	! Output luns for Cvector, Dvector
	Character*40	Cfile, Dfile	! Output filenames for Cvector, Dvector
	Integer*2	Ct_Stat (20)		! Cobetrieve status
	Integer*4	Pixel			! Pixel number
	Integer*4	Status, Err /2/, Ok /0/		! Processing status
	Logical*1	Rec_Found (2)		! Flag record found in maps 1,2
	Dictionary	'FMS_SKY'
	Record /FMS_SKY/	In_Recs (2)	! 1 record from each skymap
	Real*8		SumNum (2) / 0., 0. /	! sum of nifgs from each map
	Real*8		T2SumNum (2)		! no_freq_weight ^2 * SumNum
	Real*8		Denom			! Denominator of calculation

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Combine_C_D_Vectors  = %Loc (FMS_Normal)

C If combination is zeroth order (SF + LF or FL) then just copy C-vectors
C and merge D-vectors as stated in header notes above.

	If (Combo_Type .EQ. 'Z') Then

C  Do for each pixel.  FAC_Num_Pixels = 6144 (FUT_Params)

	   Status = Ok
	   Pixel = 0
	   Do While ((Pixel .LT. FAC_Num_Pixels) .AND. (Status .EQ. Ok))

C  Call FMS_READ_SKYMAPS to read a record from each skymap for this pixel.

	      Rstat = FMS_READ_SKYMAPS ( Lun_In, Pixel, In_Recs, Rec_Found )

	      If ( Rstat .NE. %Loc (FMS_Normal) ) Then
	         Status = Err
	      Else
	         If ( Rec_Found (1) .OR. Rec_Found (2) ) Then
	            If ( Rec_Found (1) ) Then
	               SumNum (1) = SumNum (1) +
	1                 In_Recs (1).Coad_Spec_Head.Comb_Num_Ifgs
	            EndIf
	            If ( Rec_Found (2) ) Then
	               SumNum (2) = SumNum (2) +
	1                 In_Recs (2).Coad_Spec_Head.Comb_Num_Ifgs
	            EndIf
	         EndIf
	      EndIf
	      Pixel = Pixel + 1
	   EndDo

C Rewind the file pointer to the beginning of the file

	   Call Sys$Rewind (Lun_In(1))
	   Call Sys$Rewind (Lun_In(2))

C Now copy C-vector and weight D-vectors.

	   T2SumNum (1) = SumNum (1) * No_Freq_Weights (1) ** 2
	   T2SumNum (2) = SumNum (2) * No_Freq_Weights (2) ** 2
	   Denom = T2SumNum (1) + T2SumNum (2)
	   Do i = 1, 257
	      Out_C_Vector_Rec.CVector (i) = C_Vector (i,1)
	      Out_D_Vector_Rec.Variances (i) =
	1        ( T2SumNum (1) * D_Vector (i,1) +
	2          T2SumNum (2) * D_Vector (i,2) ) / Denom
	   EndDo
	Else

C Else - not zeroth order.  Combine the vectors, low freqs then high.
C Low frequencies are always averaged after zeroth order.

	   Do i = 1, Nu_Last_Low
	      Out_C_Vector_Rec.CVector (i) =
	1        (Relative_Weights (i,1) ** 2) * C_Vector (i,1) +
	2        (Relative_Weights (i,2) ** 2) * C_Vector (i,2)
	      Out_D_Vector_Rec.Variances (i) =
	1        (Relative_Weights (i,1) ** 2) * D_Vector (i,1) +
	2        (Relative_Weights (i,2) ** 2) * D_Vector (i,2)
	   EndDo

C If combining High and Low frequency maps (combo_type=F) then just copy
C the high frequency vectors, otherwise do weighting.

	   If (Combo_Type .EQ. 'F') Then
	      If ( Specif (1)(2:2) .EQ. 'L' .OR. Specif (1)(1:2) .EQ. 'LO'
	1          .OR. Specif (1)(3:4) .EQ. 'LO' ) Then
	         Map = 2
	      Else
	         Map = 1
	      EndIf
	      Do i = Nu_Last_Low+1, 257
	         Out_C_Vector_Rec.CVector (i) = C_Vector (i,Map)
	         Out_D_Vector_Rec.Variances (i) = D_Vector (i,Map)
	      EndDo
	   Else
	      Do i = Nu_Last_Low+1, 257
	         Out_C_Vector_Rec.CVector (i) =
	1           (Relative_Weights (i,1) ** 2) * C_Vector (i,1) +
	2           (Relative_Weights (i,2) ** 2) * C_Vector (i,2)
	         Out_D_Vector_Rec.Variances (i) =
	1           (Relative_Weights (i,1) ** 2) * D_Vector (i,1) +
	2           (Relative_Weights (i,2) ** 2) * D_Vector (i,2)
	      EndDo
	   EndIf
	EndIf
	Out_C_Vector_Rec.GMT = C_Vector_Rec (1).GMT
	Out_D_Vector_Rec.GMT = D_Vector_Rec (1).GMT
	Out_C_Vector_Rec.Time (1) = C_Vector_Rec (1).Time (1)
	Out_C_Vector_Rec.Time (2) = C_Vector_Rec (1).Time (2)
	Out_D_Vector_Rec.Time (1) = D_Vector_Rec (1).Time (1)
	Out_D_Vector_Rec.Time (2) = D_Vector_Rec (1).Time (2)
	If (C_Vector_Rec (1).Channel .EQ. C_Vector_Rec (2).Channel) Then
	   Out_C_Vector_Rec.Channel = C_Vector_Rec (1).Channel
	   Out_D_Vector_Rec.Channel = D_Vector_Rec (1).Channel
	Else
	   Out_C_Vector_Rec.Channel = -1
	   Out_D_Vector_Rec.Channel = -1
	EndIf
	If (C_Vector_Rec(1).Scan_Length .EQ. C_Vector_Rec(2).Scan_Length) Then
	   Out_C_Vector_Rec.Scan_Length = C_Vector_Rec (1).Scan_Length
	   Out_D_Vector_Rec.Scan_Length = D_Vector_Rec (1).Scan_Length
	Else
	   Out_C_Vector_Rec.Scan_Length = -1
	   Out_D_Vector_Rec.Scan_Length = -1
	EndIf
	If (C_Vector_Rec(1).Scan_Speed .EQ. C_Vector_Rec(2).Scan_Speed) Then
	   Out_C_Vector_Rec.Scan_Speed = C_Vector_Rec (1).Scan_Speed
	   Out_D_Vector_Rec.Scan_Speed = D_Vector_Rec (1).Scan_Speed
	Else
	   Out_C_Vector_Rec.Scan_Speed = -1
	   Out_D_Vector_Rec.Scan_Speed = -1
	EndIf
	If (C_Vector_Rec(1).Model_Label .EQ. C_Vector_Rec(2).Model_Label) Then
	   Out_C_Vector_Rec.Model_Label = C_Vector_Rec (1).Model_Label
	   Out_D_Vector_Rec.Model_Label = D_Vector_Rec (1).Model_Label
	Else
	   Out_C_Vector_Rec.Model_Label = 'Mixed'
	   Out_D_Vector_Rec.Model_Label = 'Mixed'
	EndIf
	If (C_Vector_Rec(1).Galat_Exc .EQ. C_Vector_Rec(2).Galat_Exc) Then
	   Out_C_Vector_Rec.Galat_Exc = C_Vector_Rec (1).Galat_Exc
	   Out_D_Vector_Rec.Galat_Exc = D_Vector_Rec (1).Galat_Exc
	Else
	   Out_C_Vector_Rec.Galat_Exc = -1.
	   Out_D_Vector_Rec.Galat_Exc = -1.
	EndIf
	Out_C_Vector_Rec.NCal_Ifgs =
	1   No_Freq_Weights (1) * Dble (C_Vector_Rec (1).NCal_Ifgs) +
	2   No_Freq_Weights (2) * Dble (C_Vector_Rec (2).NCal_Ifgs)
	Out_D_Vector_Rec.NSky_Ifgs = 
	1   No_Freq_Weights (1) * Dble (D_Vector_Rec (1).NSky_Ifgs) +
	2   No_Freq_Weights (2) * Dble (D_Vector_Rec (2).NSky_Ifgs)

C  Construct the output file names.  Get Logical Unit Numbers.  Open the
C  output C & D Vector files.  Write the records.  Close the files.

	Cfile = Arch_Out // 'FMS_CVS_' // Specif_Out // '.' // Model_Ext
	Dfile = Arch_Out // 'FMS_VAR_' // Specif_Out // '.' // Model_Ext

	Rstat = Lib$Get_Lun (Lun_C)
	If (Rstat .NE. SS$_Normal) Then
	   FMS_Combine_C_D_Vectors = %Loc (FMS_Abort)
	   Call Lib$Signal (FMS_Lunerr, %Val(1), %Val(Rstat))
	Else
	   Open ( Unit=Lun_C, File=Cfile, Status='New', Iostat=Rstat,
	1         Useropen=CT_Connect_Write )

	   If (Rstat .NE. 0) Then
	      FMS_Combine_C_D_Vectors = %Loc (FMS_Abort)
	      Call Lib$Signal (FMS_Openerr, %Val(2), Cfile, %Val(Rstat))
	   Else
	      If (Report) Write (Lun_Rep, 10) 'Opened', Cfile
  10	      Format (2X, 'Successfully ', A, 1X, A)
	      Call CT_Write_Arcv ( , Lun_C, Out_C_Vector_Rec, CT_Stat )
	      If (CT_Stat(1) .NE. CTP_Normal) Then
	         FMS_Combine_C_D_Vectors = %Loc (FMS_Abort)
	         Call Lib$Signal (FMS_Writerr, %Val(2), Cfile, %Val(Rstat))
	      Else
	         Call CT_GMT_To_Binary ( Ref_Start_GMT, Ref_Start_ADT )
	         Call CT_GMT_To_Binary ( Ref_Stop_GMT, Ref_Stop_ADT )
	         Rstat = CCT_Set_TTG_Time_Range ( Lun_C, Ref_Start_ADT,
	1                                         Ref_Stop_ADT )
	         If (.NOT. Rstat ) Then
	            FMS_Combine_C_D_Vectors = %Loc (FMS_Abort)
	            Call Lib$Signal (FMS_TTGseterr, %Val(1), Cfile)
	            Call Lib$Signal (%Val (Rstat) )
	         Else
	            Call CT_Close_Arcv (, Lun_C, CT_Stat)
	            If (CT_Stat(1) .NE. CTP_Normal) Then
	               FMS_Combine_C_D_Vectors = %Loc (FMS_Abort)
	               Call Lib$Signal (FMS_Closerr,%Val(2),Cfile,%Val(Rstat))
	            Else
	               If (Report) Write (Lun_Rep, 10) 'Closed', Cfile
	            EndIf
	         EndIf
	      EndIf
	   EndIf
	   Call Lib$Free_Lun (Lun_C)
	EndIf
	If (FMS_Combine_C_D_Vectors .EQ. %Loc (FMS_Normal)) Then
	   Rstat = Lib$Get_Lun (Lun_D)
	   If (Rstat .NE. SS$_Normal) Then
	      FMS_Combine_C_D_Vectors = %Loc (FMS_Abort)
	      Call Lib$Signal (FMS_Lunerr, %Val(1), %Val(Rstat))
	   Else
	      Open ( Unit=Lun_D, File=Dfile, Status='New', Iostat=Rstat,
	1            Useropen=CT_Connect_Write )

	      If (Rstat .NE. 0) Then
	         FMS_Combine_C_D_Vectors = %Loc (FMS_Abort)
	         Call Lib$Signal (FMS_Openerr, %Val(2), Dfile, %Val(Rstat))
	      Else
	         If (Report) Write (Lun_Rep, 10) 'Opened', Dfile
	         Call CT_Write_Arcv ( , Lun_D, Out_D_Vector_Rec, CT_Stat )
	         If (CT_Stat(1) .NE. CTP_Normal) Then
	            FMS_Combine_C_D_Vectors = %Loc (FMS_Abort)
	            Call Lib$Signal (FMS_Writerr, %Val(2), Dfile, %Val(Rstat))
	         Else
	            Rstat = CCT_Set_TTG_Time_Range ( Lun_C, Ref_Start_ADT,
	1                                            Ref_Stop_ADT )
	            If ( .NOT. Rstat ) Then
	               FMS_Combine_C_D_Vectors = %Loc (FMS_Abort)
	               Call Lib$Signal (FMS_TTGseterr, %Val(1), Dfile )
	               Call Lib$Signal ( %Val (Rstat) )
	            Else
	               Call CT_Close_Arcv (, Lun_D, CT_Stat)
	               If (CT_Stat(1) .NE. CTP_Normal) Then
	                  FMS_Combine_C_D_Vectors = %Loc (FMS_Abort)
	                  Call Lib$Signal ( FMS_Closerr, %Val(2), Dfile,
	1                                   %Val(Rstat) )
	               Else
	                  If (Report) Write (Lun_Rep, 10) 'Closed', Dfile
	               EndIf
	            EndIf
	         EndIf
	      EndIf
	      Call Lib$Free_Lun (Lun_D)
	   EndIf
	EndIf
	Return
	End
