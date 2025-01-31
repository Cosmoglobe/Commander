	Integer * 4 Function FUT_Temp_List (rawtemps, rawsigs,
     &                GRTWT, GRTTRANS, combswitch, singlifg, outtemps, outsigs)
C------------------------------------------------------------------------------
C    PURPOSE:  (From "Requirements for FUT_Temperature_List, vsn 1991 Jan 15",
C		Rich Isaacman, GRC--revised 1991 Aug 26, R. E. Eplee)
C	Reduce the list of 64 GRT readings (A and B sides and high and low
C	    current readings for 12 GRTs + 4 cal resistors) down to a short
C	    list of ten "definitive" temperatures.  For each object (e.g., XCAL)
C	    form a weighted average of the high and low current readings
C	    for each GRT to get a reliable reading for that GRT, then form a
C	    weighted average of all the GRTs on each object to get a definitive
C	    temperature for the object.
C	The routine uses the coaddition sigmas as part of the averaging weights.
C	    It handles individual readings as well as multiple ones from coadds.
C	    It is also capable of combining the high/low currents alone
C	    without averaging the GRTs themselves.
C
C    AUTHOR:  Fred Shuman,  STX,  1991 Jan 31.  Supersedes FUT_Temperature_List
C	      by Wes K. Young,  SASC Technologies (now STX),  1985 October
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C    INPUT:
C	RawTemps,RawSigs  Structure  FUT_EngAnlg (*)
C	GRTWT		  Structure  FEX_GrtRawWt (or FEX_GrtCoaWt; same RDL)
C	GRTTRANS	  Structure  FEX_GrtTrans
C	combswitch	    I*4    0=Combine A&B; 1=A-side only; 2=B-side only
C	singlifg	    L*4    .True.= 1 IFG input; flag in- and out-sigmas
C                                  .False.= >1 IFG input; compute out-sigmas
C    OUTPUT:
C	outtemps(10), outsigs(10)  R*4    (**)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C   (*) FUT_EngAnlg orders the temps.     (**) Output temperatures and sigmas
C as follows within each side/current:            are in the following order:
C     1 XCal                                 1 XCal & (S5 or S6)
C     2 Sky Horn                             2 ICal
C     3 Ref Horn                             3 Sky Horn
C     4 ICal                                 4 Ref Horn
C     5 Dihedral                             5 Dihedral
C     6 RH Bolometer                         6 "Structure" (Mirror & Collimator)
C     7 RL Bolometer                         7 RH Bolometer
C     8 LH Bolometer                         8 RL Bolometer
C     9 LL Bolometer                         9 LH Bolometer
C    10 Mirror (Folding Flat)               10 LL Bolometer
C 11-14 Calibration Resistors
C    15 XCal Segment (S5 if A side; S6 if B side)
C    16 Collimator
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  NOTE:  The reference files FEX_GrtRawWt/FEX_GrtCoaWt contain the A & B GRT
C	  weights; the reference file FEX_GrtTrans contains the transition
C	  temperatures, and transition half-widths, each in FUT_EngAnlg order.
C
C    INCLUDE FILES:  FUT_Error.txt, $SSDef
C------------------------------------------------------------------------------
C    CHANGES:
C	SPRs 8835, 8899.  FEX_GRTSWT is being split into FEX_GrtTrans and
C	    FEX_GrtRawWt/FEX_GrtCoaWt.  Fred Shuman, STX, 1991 Aug 26.
C
C	SPR 9171.  In 'Combine A & B section', siga was being used to form the
C	    weight for B-side GRTs for k=2,10 (ICal, SHorn, RHorn, Dihed, and
C	    the 4 Boloms).  Fred Shuman, Hughes STX, 1991 Oct 18.
C
C	SPR 9529.  Contrary to our original assumption, FIC testing produced
C	    some coadd groups in which all the temperature values for a given
C	    GRT were identical, giving rise to a sigma of 0.0 for that GRT.
C	    This in turn causes division by zero when applying sigma-weights
C	    in this routine.  The remedy for this is a minimum, or 'floor value'
C	    for the input sigmas, which is found by noting that at T=2.53 K,
C	    delta_count=1 leads to delta_T=0.06 mK, so we take our minimum
C	    sigma to be (6e-5 K)*(T/2.53 K)^3/SQRT(12) = 1.070e-6 K *(T/K)^3.
C	    Fred Shuman, Hughes STX, 1992 Feb 27.
C	SPR 9613.  Requirements change.  The method of combining temps & sigmas
C	    from different GRTs was based on wrong assumptions and had to be
C	    changed.  Fred Shuman, Hughes STX, 1992 Mar 25.
C	SPR 9675.  Requirements change.  The method for combining hi/lo current
C	    temps & sigmas should be the same as that now being used for A/B
C	    combination (see SPR 9613).  Fred Shuman,  Hughes STX,  1992 Apr 27.
C------------------------------------------------------------------------------

	Implicit None

	Include '(FUT_Error)'
	Include '($ssdef)'

C  Calling Arguments

	Dictionary 'FUT_EngAnlg'
	Record /FUT_EngAnlg/ RawTemps	!Input temperatures
	Record /FUT_EngAnlg/ RawSigs	!Sigmas for input temperatures
	Record /FUT_EngAnlg/ ScratSigs	!'Scratch' sigmas--these can be altered
	Dictionary 'FEX_GrtRawWt'
	Record /FEX_GrtRawWt/ GrtWt	!Ref file containing GRT weights
	Dictionary 'FEX_GrtTrans'
	Record /FEX_GrtTrans/ GrtTrans	!Ref file of GRT transition temps,
					!    transition 1/2-widths
	Integer* 4  combswitch	  !Switch to control output; 0=Combine A&B
				  !    1=Output A only; 2=Output B only
	Logical* 4  singlifg	  !Flag: .True. = single IFG input; so
				  !        don't compute output sigmas;
				  !      .False.=compute output sigmas
	Real   * 4  outtemps(10)  !Combined temps to be output
	Real   * 4  outsigs(10)	  !Sigmas for combined temps to output

C  Local variables

	Integer* 4  i, j, k	  !counters
	Integer* 4  offgrts	  !lowest index number in the out-array
C			!  where all grt switches were off.  If none, then 0.
	Integer* 4  inout(10) /0,4,2,3,5,10,6,7,8,9/
C			! Which input temp corresponds to each output temp
	Real   * 4  grtwta(16), grtwtb(16)	!A- and B-side GRT weights
	Real   * 4  gwta(16), gwtb(16)	!GRT weights zeroed if temp=flag value
	Real   * 4  trantempa(16)	!A-, B-side hi/lo current temperature
	Real   * 4  trantempb(16)	!  transition values
	Real   * 4  tranhwa(16)		!A-, B-side hi/lo current temperature
	Real   * 4  tranhwb(16)		!  transition half-widths
	Real   * 4  tranlo		!transition temp minus half-width
	Real   * 4  tranhi		!transition temp plus half-width
	Real   * 4  tlo, thi		!temporaries for low, high GRT readings
	Real   * 4  slo, shi		!temporaries for low, high GRT sigmas
	Real   * 4  wlo, whi		!temporaries weights for low, high GRTs
	Real   * 4  wt			!hi/lo temperature combination weight
	Real   * 4  sumwts		!temporary sum of weights
	Real   * 4  tempa(16),tempb(16) !A-, B-side temps
	Real   * 4  siga(16), sigb(16)	!A-, B-side sigmas
	Real   * 4  temp(16)		!combined temps
	Real   * 4  infin /1.E+18/	!"infinity"  (We need to square this
	                                !  value without a Floating Overflow.)
	Real   * 4  coeff /1.070E-6/	!Coefficient for the 'floor' value for
	                                ! incoming sigmas.  Floor=coeff*(T/K)^3
C  From the message file...
	External       FUT_CombSwErr

	If (combswitch .Lt. 0 .Or. combswitch .Gt. 2) Then
	   FUT_Temp_List = %Loc(FUT_CombSwErr)
	   Call LIB$Signal(FUT_CombSwErr)
	   Return
	End If
c
c  Extract GRT weights from the input GRTWT record; extract GRT transition and
c    half-width temperatures from the input GRTTRANS record;
c    initialize the 'return status', offgrts.
c
	Do j=1,16
	   grtwta(j)    = GRTWT.GRT_A_Weight(j)
	   grtwtb(j)    = GRTWT.GRT_B_Weight(j)
	   gwta(j)      = grtwta(j)
	   gwtb(j)      = grtwtb(j)
	   trantempa(j) = GRTTRANS.GRT_A_Trans_Temp(j)
	   trantempb(j) = GRTTRANS.GRT_B_Trans_Temp(j)
	   tranhwa(j)   = GRTTRANS.GRT_A_Trans_HWid(j)
	   tranhwb(j)   = GRTTRANS.GRT_B_Trans_HWid(j)
	End Do
	offgrts = 0
c
c  Protect the input sigmas from the alterations to come
c
	Do j=1,64
	   ScratSigs.GRT(j) = RawSigs.GRT(j)
	End Do
c
c  Set the 8 hi-curr bols to "infinity"
c
	Do j=6,9
	   ScratSigs.a_hi_GRT(j) = infin
	   ScratSigs.b_hi_GRT(j) = infin
	End Do
c
c Set any positive sigmas that are too small to a minimum value
c
	Do j=1,64
	   If (ScratSigs.GRT(j) .Ge. 0.) Then
	      ScratSigs.GRT(j) = Max(ScratSigs.GRT(j), coeff*RawTemps.GRT(j)**3)
	   End If
	End Do
c
c A-side temps--Unless "B side only" was requested, combine Hi & Lo readings of
c                                                     each A-side GRT...
	If (combswitch .Ne. 2) Then
	   Do j=1,16
	      tlo = RawTemps.a_lo_grt(j)
	      thi = RawTemps.a_hi_grt(j)
	      slo = ScratSigs.a_lo_grt(j)
	      shi = ScratSigs.a_hi_grt(j)
	      If (grtwta(j) .Ne. 0 .And.
     &                  (tlo .Gt. -999. .Or. thi .Gt. -999.)) Then
	         tranlo = trantempa(j) - tranhwa(j)
	         tranhi = trantempa(j) + tranhwa(j)
	         If (thi .Lt. -999.) Then
	            wt = 0.
	         Else If (tlo .Lt. -999.) Then
	            wt = 1.
	         Else If (tlo .Lt. tranlo) Then
	            wt = 0.
	         Else If (tlo .Gt. tranhi) Then
	            wt = 1.
	         Else
	            wt = 0.5 + (tlo - trantempa(j))/tranhwa(j)/2
	         End If

	         tempa(j) = (1. - wt)*tlo + wt*thi
	         siga(j)  = Sqrt( ((1. - wt)*slo)**2 + (wt*shi)**2 )
	         If (tlo .Lt. -999.   .And.   thi .Lt. -999.) Then
	            tempa(j) = -9999.
	            siga(j)  = -9999.
	         End If
	      Else
	         tempa(j) = -9999.
	         siga(j)  = -9999.
	         gwta(j)  =     0.
	      End If	    ! (grtwta(j) .Ne. 0 .And. ...
	   End Do	  ! j=1,16
	End If		! (combswitch .Ne. 2)
c
c B-side temps--Unless "A side only" was requested, combine Hi & Lo readings of
c                                                     each B-side GRT...
	If (combswitch .Ne. 1) Then
	   Do j=1,16
	      tlo = RawTemps.b_lo_grt(j)
	      thi = RawTemps.b_hi_grt(j)
	      slo = ScratSigs.b_lo_grt(j)
	      shi = ScratSigs.b_hi_grt(j)
	      If (grtwtb(j) .Ne. 0 .And.
     &                  (tlo .Gt. -999. .Or. thi .Gt. -999.)) Then
	         tranlo = trantempb(j) - tranhwb(j)
	         tranhi = trantempb(j) + tranhwb(j)
	         If (thi .Lt. -999.) Then
	            wt = 0.
	         Else If (tlo .Lt. -999.) Then
	            wt = 1.
	         Else If (tlo .Lt. tranlo) Then
	            wt = 0.
	         Else If (tlo .Gt. tranhi) Then
	            wt = 1.
	         Else
	            wt = 0.5 + (tlo - trantempb(j))/tranhwb(j)/2
	         End If

	         tempb(j) = (1. - wt)*tlo + wt*thi
	         sigb(j)  = Sqrt( ((1. - wt)*slo)**2 + (wt*shi)**2 )
	         If (tlo .Lt. -999.   .And.   thi .Lt. -999.) Then
	            tempb(j) = -9999.
	            sigb(j)  = -9999.
	         End If
	      Else
	         tempb(j) = -9999.
	         sigb(j)  = -9999.
	         gwtb(j)  =     0.
	      End If	    ! (grtwtb(j) .Ne. 0 .And. ...
	   End Do	  ! j=1,16
	End If		! (combswitch .Ne. 1)
c
c  If "Combine A & B sides" was requested, take the merged current readings and
c  combine A and B sides to get a single definitive temperature for each object.
c  Otherwise, supply just the side requested--A or B.
c
	If (combswitch .Eq. 0) Then
c
c    XCal temp ( out#1 = in#1 + in#15 )
c
	   sumwts = gwta( 1) + gwtb( 1) + gwta(15) + gwtb(15)

	   If (sumwts .Gt. 0.) Then
	      outtemps(1) = (gwta( 1)*tempa( 1) + gwtb( 1)*tempb( 1) +
     &                       gwta(15)*tempa(15) + gwtb(15)*tempb(15))/sumwts
	      If (.Not. singlifg) Then
	         outsigs(1)  = Sqrt( (gwta( 1)*siga( 1))**2 +
     &                               (gwtb( 1)*sigb( 1))**2 +
     &                               (gwta(15)*siga(15))**2 +
     &                               (gwtb(15)*sigb(15))**2 )/sumwts
	      Else				! If 'single IFG' flag is on,
	         outsigs(1)  = -9999.		! set out-sigmas to -9999
	      End If
	   Else
	      outtemps(1) = -9999.
	      outsigs(1)  = -9999.
	      offgrts = 1
	   End If
c
c    Structure temp (out#6) = mirror mounts (in#10, aka folding flats)
c                                          + collimator (in#16)
c      (B-side collimator is unavailable)
c
	   sumwts = gwta(10) + gwtb(10) + gwta(16)

	   If (sumwts .Gt. 0.) Then
	      outtemps(6) = (gwta(10)*tempa(10) + gwtb(10)*tempb(10) +
     &                       gwta(16)*tempa(16))/sumwts
	      If (.Not. singlifg) Then
	         outsigs(6)  = Sqrt( (gwta(10)*siga(10))**2 +
     &                               (gwtb(10)*sigb(10))**2 +
     &                               (gwta(16)*siga(16))**2 ) /sumwts
	      Else				! If 'single IFG' flag is on,
	         outsigs(6)  = -9999.		! set out-sigmas to -9999
	      End If
	   Else
	      outtemps(6) = -9999.
	      outsigs(6)  = -9999.
	      If (offgrts .Eq. 0) Then
	         offgrts = 6
	      End If
	   End If
c
c    ICal, Sky Horn, Ref Horn, Dihedral, Bolometer RH, RL, LH, LL temps
c
	   Do k=2,10
	      j = inout(k)
	      If (j .Ne. 0 .And. k .Ne. 6) Then
	         sumwts = gwta(j) + gwtb(j)

	         If (sumwts .Gt. 0.) Then
	            outtemps(k) = (gwta(j)*tempa(j) + gwtb(j)*tempb(j))
     &                            /sumwts
	            If (.Not. singlifg) Then
	               outsigs(k)  = Sqrt( (gwta(j)*siga(j))**2 +
     &                                     (gwtb(j)*sigb(j))**2 )/sumwts
	            Else			! If 'single IFG' flag is on,
	               outsigs(k)  = -9999.	! set out-sigmas to -9999
	            End If
	         Else
	            outtemps(k) = -9999.
	            outsigs(k)  = -9999.
	            If (k .Lt. offgrts .Or. offgrts .Eq. 0) Then
	               offgrts = k
	            End If
	         End If		! (sumwts .Gt. 0.)
	      End If	    ! (j .Ne. 0 .And. k .Ne. 6)
	   End Do	! k=2,10
c
c If "A side only" was requested...
c
	Else If (combswitch .Eq. 1) Then
c
c    XCal temp ( out#1 = in#1 + in#15 )
c
	   sumwts = gwta( 1) + gwta(15)
	   If (sumwts .Gt. 0.) Then
	      outtemps(1) = (gwta( 1)*tempa( 1) + gwta(15)*tempa(15))/sumwts
	      If (.Not. singlifg) Then
	         outsigs(1)  = Sqrt( (gwta( 1)*siga( 1))**2 +
     &                               (gwta(15)*siga(15))**2 )/sumwts
	      Else				! If 'single IFG' flag is on,
	         outsigs(1)  = -9999.		! set out-sigmas to -9999
	      End If
	   Else
	      outtemps(1) = -9999.
	      outsigs(1)  = -9999.
	      offgrts = 1
	   End If
c
c    Structure temp (out#6) = mirror mounts (in#10, aka folding flats)
c                                          + collimator (in#16)
	   sumwts = gwta(10) + gwta(16)
	   If (sumwts .Gt. 0.) Then
	      outtemps(6) = (gwta(10)*tempa(10) + gwta(16)*tempa(16))/sumwts
	      If (.Not. singlifg) Then
	         outsigs(6)  = Sqrt( (gwta(10)*siga(10))**2 +
     &                               (gwta(16)*siga(16))**2 )/sumwts
	      Else				! If 'single IFG' flag is on,
	         outsigs(6)  = -9999.		! set out-sigmas to -9999
	      End If
	   Else
	      outtemps(6) = -9999.
	      outsigs(6)  = -9999.
	      If (offgrts .Eq. 0) Then
	         offgrts = 6
	      End If
	   End If
c
c    Sky Horn, Ref Horn, ICal, Dihedral, Bolometer RH, RL, LH, LL temps
c
	   Do k=2,10
	      j = inout(k)
	      If (j .Ne. 0 .And. k .Ne. 6) Then
	         outtemps(k) = tempa(j)
	         If (.Not. singlifg) Then
	            outsigs(k)  = siga(j)
	         Else				! If 'single IFG' flag is on,
	            outsigs(k)  = -9999.	! set out-sigmas to -9999
	         End If
	         If ( tempa(j) .Lt. -999.  .And.
     &                (k .Lt. offgrts .Or. offgrts .Eq. 0) ) Then
	            offgrts = k
	         End If
	      End If
	   End Do
c
c If "B side only" was requested...
c
	Else If (combswitch .Eq. 2) Then
c
c    XCal temp
c
	   sumwts = gwtb( 1) + gwtb(15)
	   If (sumwts .Gt. 0.) Then
	      outtemps(1) = (gwtb( 1)*tempb( 1) + gwtb(15)*tempb(15))/sumwts
	      If (.Not. singlifg) Then
	         outsigs(1)  = Sqrt( (gwtb( 1)*sigb( 1))**2 +
     &                               (gwtb(15)*sigb(15))**2 )/sumwts
	      Else				! If 'single IFG' flag is on,
	         outsigs(1)  = -9999.		! set out-sigmas to -9999
	      End If
	   Else
	      outtemps(1) = -9999.
	      outsigs(1)  = -9999.
	      offgrts = 1
	   End If
c
c    ICal, Sky Horn, Ref Horn, Dihedral, Structure (= mirror mounts, since
c       B-side collimator is unavailable), Bolometer RH, RL, LH, LL,
c
	   Do k=2,10
	      j = inout(k)
	      If (j .Ne. 0) Then
	         outtemps(k) = tempb(j)
	         If (.Not. singlifg) Then
	            outsigs(k)  = sigb(j)
	         Else				! If 'single IFG' flag is on,
	            outsigs(k)  = -9999.	! set out-sigmas to -9999
	         End If
	         If ( tempb(j) .Lt. -999.  .And.
     &                (k .Lt. offgrts .Or. offgrts .Eq. 0) ) Then
	            offgrts = k
	         End If
	      End If
	   End Do

	End If		! (combswitch .Eq. 0)
c
c   Return the lowest index number in the out-array where a 'flag' value was
c     returned because all the GRT switches in its make-up were off; if none, 0.
c
	FUT_Temp_List = offgrts

	Return
	End
