	Integer*4  Function  FMS_SPECTRA_VARIANCES
	1                  ( In_Skymap, Weights, Specif, Order, Offset, Gain,
	2                    Nu_Last_Low, Relative_Weights, Combo_Type, In_Rec,
	3                    Out_Rec, Nifgs, No_Freq_Weights, C_Vector )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                          FMS_SPECTRA_VARIANCES.FOR
C
C  Purpose:  Combine spectra and variances.  Compute chi squared for goodness
C            of combination.  Combine number of ifgs.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C  Modified: L. Rosen, HSTX, 16 August 1994.  If input is FCS then the field
C     SPEC_DATA.REAL_IMAG_VAR is replaced by the new COMB_SPEC_DATA.CVS_VAR
C     field.  This variance is calculated for each input spectrum as
C     cvs_var = C / N    where C is the C-Vector variance (c^2) and N is the
C     glitch rate weighted number of ifgs.  If the combination is zeroth order
C     then the N in the denominator gets multiplied by T, the No_Freq_Weights.
C     The combined CVS_VAR field is combined the same as other variance fields
C     and stored in the FMS record.
C     The variances are combined using Vout = (w1^2 V1 + w2^2 V2)/(w1+w2)^2
C     where w = TN for 0th order or RelativeWeight*N for higher order.
C  Modified: L. Rosen, HSTX, 24 August 1994.  Chi square is now calculated
C     for real and imaginary parts.
C  Modified: L. Rosen, HSTX, October & November 1994.  Correction spectra are
C     only applied for FIRST and THIRD order merges.  The offset correction
C     is applied to ALL spectra in the low frequencies, those below the
C     transition frequency specified by the command line.  The gain correction
C     is only applied to the high frequencies of maps from high frequency
C     detectors, and only in the first order.
C  18 November 1994, L. Rosen, HSTX:
C     HISL, LOSL, HRES merges use FCS input and are 1st order merges, and
C     should use frequency dependent weights for CVS calculation.
C     Variance calculations for High frequencies for 1st order need same
C     gain correction as the spectra.
C  5 December 1994, L. Rosen, HSTX:
C     Ken Jensen has demonstrated that we do need to apply a gain and offset
C     correction in second order merges.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  External

	External	FMS_Normal

C  Passed parameters.

	Character*33	In_Skymap (2)		! Input skymaps
	Real*8		Weights (257,2)		! C or D vector
	Character*4	Specif (2)		! "chanscan" of input skymaps
	Character*6	Order			! Command combination order
	Complex*8	Offset (257,2)		! Gain correction spectrum
	Real*4		Gain (257,2)		! Offset correction spectrum
	Integer*2	Nu_Last_Low	! Highest index for low frequency band
	Real*8		Relative_Weights (257,2)
	Character*1	Combo_Type		! Combination type.
			! Z=Zeroth, L=Low+Low, H=High+High, F=High+Low, M=Other
	Dictionary	'FMS_SKY'
	Record /FMS_SKY/	In_Rec (2)	! 1 record from each skymap
	Record /FMS_SKY/	Out_Rec		! Output record for pixel
	Real*8		Nifgs (2)	! Glitch rate weighted number of ifgs
	Real*8		No_Freq_Weights (2)	! Freq indep. weight of skymaps
	Real*8		C_Vector (257,2)	! C-vector for skymap 1 and 2

C  Local

	Integer*2	Imax		! Upper limit of loop
	Integer*2	Map
	Integer*2	HMap
	Integer*2	Freq
	Real*8		Denom		! real*8 denominator
	Real*8		NumerR		! real*8 numerator
	Complex*16	NumerC		! complex*16 numerator
	Complex*16	Ispec (2), Ospec	! Double complex of input and
						! output spectra
	Real*8		XR		! Parts of Real chi square computation
	Real*8		XI		! Parts of Imag chi square computation
	Real*8		TN (2)		! = No_Freq_Weights * Nifgs
	Real*8		TN2 (2)		! TN^2
	Real*8		Denom2		! Denom * Denom
	Real*8		CVS (2)		! If input maps are FCS, this = cvs_var
	Real*8		Part		! Part of a calculation
	Real*8		G		! real*8 Gain correction for frequency

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Spectra_Variances = %Loc (FMS_Normal)
	XR = 0.
	XI = 0.

C  Use glitch rate weighted number of ifgs.  Store the sum of the two.

	Do Map = 1, 2
	   Nifgs (Map) = Dble (In_Rec (Map).Coad_Spec_Head.Comb_Num_Ifgs)
	   TN (Map) = No_Freq_Weights (Map) * Nifgs (Map)
	   TN2 (Map) = TN (Map) * TN (Map)
	EndDo
	Out_Rec.Coad_Spec_Head.Comb_Num_Ifgs = TN (1) + TN (2)

C  If Zeroth Order Case : use frequency independent weighting

	If ( Combo_Type .EQ. 'Z' ) Then
	   Denom = TN (1) + TN (2)
	   Denom2 = Denom * Denom

C	   If Denom = 0 Then Spec, Real_Var, Imag_Var are all 0.

	   If (Denom .NE. 0.) Then
	      If (Specif(1)(2:2) .EQ. 'H') Then
	         Imax = 257
	      Else
	         Imax = Nu_Last_Low
	      EndIf
	      Do Freq = 1, Imax
	         Ispec (1) = DCMPLX (In_Rec (1).Spec_Data.Spec (Freq))
	         Ispec (2) = DCMPLX (In_Rec (2).Spec_Data.Spec (Freq))
	         Ospec = ( TN (1) * Ispec (1) + TN (2) * Ispec (2) ) / Denom
	         Out_Rec.Spec_Data.Spec (Freq) = CMPLX (Ospec)

	         Out_Rec.Spec_Data.Real_Var (Freq) = 
	1           ( TN2 (1) * In_Rec(1).Spec_Data.Real_Var (Freq) +
	2             TN2 (2) * In_Rec(2).Spec_Data.Real_Var (Freq) ) /
	3           Denom2
	         Out_Rec.Spec_Data.Imag_Var (Freq) = 
	1           ( TN2 (1) * In_Rec(1).Spec_Data.Imag_Var (Freq) +
	2             TN2 (2) * In_Rec(2).Spec_Data.Imag_Var (Freq) ) /
	3           Denom2

C If either map is missing then the chi square contribution is 0.
C Note that for zeroth order c-vectors are the same for map 1 and map 2.

	         If ( TN (1) .NE. 0. .AND. TN (2) .NE. 0. ) Then
	            If ( C_Vector (Freq, 1) .NE. 0. ) Then
	               XR = XR + ( TN (1) / C_Vector (Freq,1) ) *
	1                    ( DREAL (Ispec (1)) - DREAL (Ospec) )**2
	               XI = XI + ( TN (1) / C_Vector (Freq,1) ) *
	1                    ( DIMAG (Ispec (1)) - DIMAG (OSpec) )**2
	            EndIf
	            If ( C_Vector (Freq, 2) .NE. 0.) Then
	               XR = XR + ( TN (2) / C_Vector (Freq,2) ) *
	1                    ( DREAL (Ispec (2)) - DREAL (Ospec) )**2
	               XI = XI + ( TN (2) / C_Vector (Freq,2) ) *
	1                    ( DIMAG (Ispec (2)) - DIMAG (OSpec) )**2
	            EndIf
	         EndIf
	      EndDo			! Frequency loop
	   EndIf

	Else				! Else not zeroth order merge

C  Low frequencies.

	   Do Freq = 1, Nu_Last_Low
	      Denom = Relative_Weights (Freq,1) * Nifgs (1) +
	1             Relative_Weights (Freq,2) * Nifgs (2)

C	   If Denom = 0 Then Spec, Real_Var, Imag_Var are all 0.
C          Also, at all frequencies above Nu_Last_Low are all 0.

	      If (Denom .NE. 0.) Then
	         Denom2 = Denom * Denom
	         Ispec (1) = DCMPLX (In_Rec (1).Spec_Data.Spec (Freq))
	         Ispec (2) = DCMPLX (In_Rec (2).Spec_Data.Spec (Freq))

C  The following calculates the numerator of the combined spectrum.  The
C  same formula applies in all cases.
C  Of course this is only the low frequency part of the spectrum, where
C  only an offset correction is used.

	         NumerC = Relative_Weights (Freq,1) * Nifgs (1) *
	1           ( Ispec (1) - Offset (Freq,1) )
	2           +     Relative_Weights (Freq,2) * Nifgs (2) *
	3           ( Ispec (2) - Offset (Freq,2) )

	         Ospec = NumerC / Denom
	         Out_Rec.Spec_Data.Spec (Freq) = CMPLX (Ospec)

	         NumerR = (Relative_Weights(Freq,1) * Nifgs(1))**2 *
	1           DBLE (In_Rec(1).Spec_Data.Real_Var (Freq)) +
	2                (Relative_Weights(Freq,2) * Nifgs(2))**2 *
	3           DBLE (In_Rec(2).Spec_Data.Real_Var (Freq))

	         Out_Rec.Spec_Data.Real_Var (Freq) = SNGL (NumerR / Denom2)

	         NumerR = ( Relative_Weights (Freq,1) * Nifgs (1) )**2 *
	1           DBLE (In_Rec (1).Spec_Data.Imag_Var (Freq)) +
	2           ( Relative_Weights (Freq,2) * Nifgs (2) )**2 *
	3           DBLE (In_Rec (2).Spec_Data.Imag_Var (Freq))

	         Out_Rec.Spec_Data.Imag_Var (Freq) = SNGL (NumerR / Denom2)

	      EndIf

C If either map is missing then the chi square contribution is 0.
C Offset may be 0.  "Part" business is to protect from VMS underflow
C problems.

	      If ( Nifgs (1) .NE. 0. .AND. Nifgs (2) .NE. 0. ) Then
	         Do Map = 1, 2
	            If (C_Vector (Freq, Map) .NE. 0.) Then
	               Part = DREAL (Ispec (Map)) - DREAL (Ospec) -
	1                     DBLE (Offset (Freq,Map))
	               If (ABS (Part) .GT. 1.E-19) Then
	                  XR = XR + ( ( Part )**2 ) *
	1                        ( Nifgs (Map) / C_Vector (Freq, Map) )
	               EndIf
	               Part = DIMAG (Ispec (Map)) - DIMAG (Ospec) -
	1                     DBLE (AIMAG (Offset (Freq,Map)))
	               If (ABS (Part) .GT. 1.E-19) Then
	                  XI = XI + ( ( Part )**2 ) *
	1                           ( Nifgs (Map) / C_Vector (Freq, Map) )
	               EndIf
	            EndIf
	         EndDo			! for each map
	      EndIf
	   EndDo			! for all low frequencies

C  High frequencies.
C  If low + low (combo_type = L) then spectra and variances are 0 (skip loop).

	   If (Combo_Type .NE. 'L') Then

C  If High + Low or Low + High (combo_type='F') then copy high map data.
C  No additions to chi squared.
C  If Order is FIRST then the high frequency gets gain correction.

	      If (Combo_Type .EQ. 'F') Then
	         If ((Specif(1)(2:2) .EQ. 'H') .OR. (Specif (1)(1:2) .EQ. 'HI')
	1          .OR. (Specif (1)(3:4) .EQ. 'HI')) Then	! Map 1 is H
	            HMap = 1
	         Else
	            HMap = 2
	         EndIf
	         If (Order .EQ. 'FIRST' .OR. Order .EQ. 'SECOND') Then
	            Do Freq = Nu_Last_Low+1, 257
	               Out_Rec.Spec_Data.Spec (Freq) =
	1                 CMPLX ( DCMPLX (In_Rec(HMap).Spec_Data.Spec(Freq))
	2                         * DBLE (Gain (Freq,HMap)) )
	               Out_Rec.Spec_Data.Real_Var (Freq) =
	1                 In_Rec (HMap).Spec_Data.Real_Var (Freq) 
	2                 * (DBLE (Gain (Freq,HMap))**2)
	               Out_Rec.Spec_Data.Imag_Var (Freq) =
	1                 In_Rec (HMap).Spec_Data.Imag_Var (Freq)
	2                 * (DBLE (Gain (Freq,HMap))**2)
	            EndDo
	         Else
	            Do Freq = Nu_Last_Low+1, 257
	               Out_Rec.Spec_Data.Spec (Freq) =
	1                 In_Rec(HMap).Spec_Data.Spec(Freq)
	               Out_Rec.Spec_Data.Real_Var (Freq) =
	1                 In_Rec (HMap).Spec_Data.Real_Var (Freq) 
	               Out_Rec.Spec_Data.Imag_Var (Freq) =
	1                 In_Rec (HMap).Spec_Data.Imag_Var (Freq)
	            EndDo
	         EndIf

	      ElseIf (Combo_Type .EQ. 'H' .OR. Combo_Type .EQ. 'M') Then

C If first OR second order apply gain correction.

	         If (Order .EQ. 'FIRST' .OR. Order .EQ. 'SECOND') Then
	            Do Freq = Nu_Last_Low+1, 257
	               Denom = Relative_Weights (Freq,1) * Nifgs (1) +
	1                      Relative_Weights (Freq,2) * Nifgs (2)

	               If (Denom .EQ. 0.) Then
	                  Out_Rec.Spec_Data.Spec (Freq) = (0.,0.)
	                  Out_Rec.Spec_Data.Real_Var (Freq) = 0.
	                  Out_Rec.Spec_Data.Imag_Var (Freq) = 0.
	               Else
	                  Denom2 = Denom * Denom
	                  Ispec (1) = DCMPLX (In_Rec (1).Spec_Data.Spec (Freq))
	                  Ispec (2) = DCMPLX (In_Rec (2).Spec_Data.Spec (Freq))

	                  NumerC = Relative_Weights (Freq,1) * Nifgs (1) *
	1                    Ispec (1) * Gain (Freq,1)
	2                    + Relative_Weights (Freq,2) * Nifgs (2) *
	3                    Ispec (2) * Gain (Freq,2)

	                  Ospec = NumerC / Denom 
	                  Out_Rec.Spec_Data.Spec (Freq) = CMPLX (Ospec)

	                  NumerR = ( Relative_Weights (Freq,1) * Nifgs(1) *
	1                            Gain (Freq,1) )**2 *
	2                    DBLE (In_Rec(1).Spec_Data.Real_Var (Freq))
	3                    +     ( Relative_Weights (Freq,2) * Nifgs(2) *
	4                            Gain (Freq,2) )**2 *
	5                    DBLE (In_Rec(2).Spec_Data.Real_Var (Freq))

	                  Out_Rec.Spec_Data.Real_Var (Freq) =
	1                    SNGL (NumerR/Denom2)

	                  NumerR = ( Relative_Weights (Freq,1) * Nifgs (1) *
	1                            Gain (Freq,1) )**2 *
	2                    DBLE (In_Rec (1).Spec_Data.Imag_Var (Freq))
	3                    +     ( Relative_Weights (Freq,2) * Nifgs (2) *
	4                            Gain (Freq,2) )**2 *
	5                    DBLE (In_Rec (2).Spec_Data.Imag_Var (Freq))

	                  Out_Rec.Spec_Data.Imag_Var (Freq) =
	1                    SNGL (NumerR/Denom2)

C If either map is missing then the chi square contribution is 0.

	                  If ( Nifgs (1) .NE. 0. .AND. Nifgs (2) .NE. 0. ) Then
	                     Do Map = 1, 2
	                        If (C_Vector (Freq,Map) .NE. 0.) Then
	                           G = DBLE (Gain (Freq,Map))
	                           NumerR=(DREAL(Ispec(Map))*G-DREAL(OSpec))**2
	                           Part = Nifgs(Map)/(C_Vector(Freq,Map)*G**2)
	                           XR = XR + NumerR * Part
	                           NumerR=(DIMAG(Ispec(Map))*G-DIMAG(OSpec))**2
	                           XI = XI + NumerR * Part
	                        EndIf
	                     EndDo			! for each map
	                  EndIf
	               EndIf			! If Denom .EQ. 0
	            EndDo			! For high frequencies
	         Else				! Not first order
	            Do Freq = Nu_Last_Low+1, 257
	               Denom = Relative_Weights (Freq,1) * Nifgs (1) +
	1                      Relative_Weights (Freq,2) * Nifgs (2)

	               If (Denom .EQ. 0.) Then
	                  Out_Rec.Spec_Data.Spec (Freq) = (0.,0.)
	                  Out_Rec.Spec_Data.Real_Var (Freq) = 0.
	                  Out_Rec.Spec_Data.Imag_Var (Freq) = 0.
	               Else
	                  Denom2 = Denom * Denom
	                  Ispec (1) = DCMPLX (In_Rec (1).Spec_Data.Spec (Freq))
	                  Ispec (2) = DCMPLX (In_Rec (2).Spec_Data.Spec (Freq))

	                  NumerC = Relative_Weights (Freq,1) * Nifgs (1) *
	1                             Ispec (1)
	2                        + Relative_Weights (Freq,2) * Nifgs (2) *
	3                             Ispec (2)

	                  Ospec = NumerC / Denom 
	                  Out_Rec.Spec_Data.Spec (Freq) = CMPLX (Ospec)

	                  NumerR = (Relative_Weights (Freq,1) * Nifgs(1))**2 *
	1                    DBLE (In_Rec(1).Spec_Data.Real_Var (Freq)) +
	2                         (Relative_Weights (Freq,2) * Nifgs(2))**2 *
	3                    DBLE (In_Rec(2).Spec_Data.Real_Var (Freq))

	                  Out_Rec.Spec_Data.Real_Var (Freq) =
	1                    SNGL (NumerR/Denom2)

	                  NumerR = (Relative_Weights (Freq,1) * Nifgs (1))**2 *
	1                    DBLE (In_Rec (1).Spec_Data.Imag_Var (Freq)) +
	2                          (Relative_Weights (Freq,2) * Nifgs (2))**2 *
	3                    DBLE (In_Rec (2).Spec_Data.Imag_Var (Freq))

	                  Out_Rec.Spec_Data.Imag_Var (Freq) =
	1                    SNGL (NumerR/Denom2)

C If either map is missing then the chi square contribution is 0.

	                  If ( Nifgs (1) .NE. 0. .AND. Nifgs (2) .NE. 0. ) Then
	                     Do Map = 1, 2
	                        If (C_Vector (Freq,Map) .NE. 0.) Then
	                           XR = XR +(DREAL(Ispec(Map))-DREAL(OSpec))**2
	1                                 * ( Nifgs(Map) / C_Vector(Freq,Map) )
	                           XI = XI +(DIMAG(Ispec(Map))-DIMAG(OSpec))**2
	1                                 * ( Nifgs(Map) / C_Vector(Freq,Map) )
	                        EndIf
	                     EndDo		! for each map
	                  EndIf
	               EndIf			! If Denom .EQ. 0
	            EndDo			! For high frequencies
	         EndIf				! First or Third order
	      EndIf				! If combo_type
	   EndIf				! If not L combo
	EndIf					! If not Z combo

	Out_Rec.Coad_Spec_Head.Real_Comb_Chi_Square = SNGL (XR)
	Out_Rec.Coad_Spec_Head.Imag_Comb_Chi_Square = SNGL (XI)

C CVS Variances

	If ( In_Skymap (1)(1:3) .EQ. 'FCS' ) Then

	   If ( Combo_Type .EQ. 'Z' .AND. Denom2 .NE. 0. ) Then
	      Denom = TN (1) + TN (2)
	      Denom2 = Denom * Denom
	      Do Freq = 1, 257
	         Do Map = 1, 2
	            If ( TN (Map) .NE. 0.) Then
	               CVS (Map) = C_Vector (Freq,Map) / TN (Map)
	            EndIf
	         EndDo
	         Out_Rec.Spec_Data.Cvs_Var (Freq) =
	1           ( CVS (1) * TN2 (1) + CVS (2) * TN2 (2) ) / Denom2
	      EndDo
	   Else				! Not zeroth order, but FCS input
	      Do Freq = 1, 257
	         Denom = Relative_Weights (Freq,1) * Nifgs (1) +
	1                Relative_Weights (Freq,2) * Nifgs (2)
	         Denom2 = Denom * Denom
	         Do Map = 1, 2
	            If ( Nifgs (Map) .EQ. 0.) Then
	               CVS (Map) = 0.
	            Else
	               CVS (Map) = C_Vector (Freq,Map) / Nifgs (Map)
	            EndIf
	         EndDo
	         If ( Denom2 .NE. 0 ) Then
	            If ( ((Order .EQ. 'FIRST') .OR. (Order .EQ. 'SECOND'))
	1                .AND. (Freq .GT. Nu_Last_Low) ) Then
	               NumerR = ( Relative_Weights (Freq,1) * Nifgs (1) *
	1                         Gain (Freq,1) )**2 * CVS(1)
	2              +        ( Relative_Weights (Freq,2) * Nifgs (2) *
	3                         Gain (Freq,2) )**2 * CVS(2)
	            Else
	               NumerR = ( Relative_Weights (Freq,1) * Nifgs (1) )**2 *
	1                 CVS (1) +
	2                       ( Relative_Weights (Freq,2) * Nifgs (2) )**2 *
	3                 CVS(2)
	            EndIf
	            Out_Rec.Spec_Data.Cvs_Var (Freq) = SNGL (NumerR / Denom2)
	         EndIf
	      EndDo				! for each frequency
	   EndIf				! Zeroth order?
	Else					! Else not FCS input
	   If ((Order .EQ. 'FIRST') .OR. (Order .EQ. 'SECOND')) Then
	      Do Freq = 1, 257
	         Denom = Relative_Weights (Freq,1) * Nifgs (1) +
	1                Relative_Weights (Freq,2) * Nifgs (2)
	         Denom2 = Denom * Denom
	         If (Freq .GT. Nu_Last_Low .AND. Combo_Type .EQ. 'F') Then
	            Out_Rec.Spec_Data.Cvs_Var (Freq) = Gain (Freq,HMap)**2 *
	1              In_Rec (HMap).Spec_Data.Cvs_Var (Freq)
	         ElseIf (Denom .NE. 0.) Then
	            If (Freq .GT. Nu_Last_Low) Then
	               NumerR = ( Relative_Weights (Freq,1) * Nifgs (1) *
	1                         Gain (Freq,1) )**2 *
	2                 In_Rec (1).Spec_Data.Cvs_Var (Freq)
	3                 +     ( Relative_Weights (Freq,2) * Nifgs (2) *
	4                         Gain (Freq,2) )**2 *
	5                 In_Rec (2).Spec_Data.Cvs_Var (Freq)
	            Else
	               NumerR = ( Relative_Weights (Freq,1) * Nifgs (1) )**2 *
	1                 In_Rec (1).Spec_Data.Cvs_Var (Freq)
	2                 +     ( Relative_Weights (Freq,2) * Nifgs (2) )**2 *
	3                 In_Rec (2).Spec_Data.Cvs_Var (Freq)
	            EndIf
	            Out_Rec.Spec_Data.Cvs_Var (Freq) = SNGL(NumerR / Denom2)
	         EndIf
	      EndDo
	   Else					! Not 1st order, FMS input
	      Do Freq = 1, 257
	         Denom = Relative_Weights (Freq,1) * Nifgs (1) +
	1                Relative_Weights (Freq,2) * Nifgs (2)
	         Denom2 = Denom * Denom
	         If (Freq .GT. Nu_Last_Low .AND. Combo_Type .EQ. 'F') Then
	            Out_Rec.Spec_Data.Cvs_Var (Freq) =
	1              In_Rec (HMap).Spec_Data.Cvs_Var (Freq)

	         ElseIf (Denom .NE. 0.) Then
	            NumerR = ( Relative_Weights (Freq,1) * Nifgs (1) )**2 *
	1              In_Rec (1).Spec_Data.Cvs_Var (Freq) +
	2                    ( Relative_Weights (Freq,2) * Nifgs (2) )**2 *
	3              In_Rec (2).Spec_Data.Cvs_Var (Freq)

	            Out_Rec.Spec_Data.Cvs_Var (Freq) = SNGL (NumerR / Denom2)
	         EndIf
	      EndDo
	   EndIf
	EndIf
	Return
	End
