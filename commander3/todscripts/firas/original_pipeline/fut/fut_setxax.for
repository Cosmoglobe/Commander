	Integer*4 Function FUT_SetxAx(iflag, fake_it, mspeed, mlen, chan,
     &                                ngroup, scimode, firstsamp, nyquist,
     &                                startx, spacex)

C------------------------------------------------------------------------------
C      Function:
C         Subroutine to set x-axis parameters for FUT_DISPLAY.
C
C      Author: W.K. Young
C              SASC Technologies
C              August 1985
C
C      CALLING SEQUENCE:
C	status = FUT_SetxAx(XAXIS TYPE, FAKEIT MODE, MTM SPEED, MTM LENGTH,
C     &               CHANNEL, ADDS PER GROUP, SCIENCE MODE, FIRSTSAMP,
C     &               NYQUIST REFERENCE, XAXIS START, XAXIS SPACING)
C
C      INPUT:
C	IFLAG	I*4	TYPE OF X-AXIS SCALE [ 0=position;
C			   1=spatial freq (coadd); -1=spatial freq (noisetest) ]
C	FAKE_IT	I*4	FAKEIT MODE
C	MSPEED	I*4	MTM SPEED   [ 0=slow;  1=fast ]
C	MLEN	I*4	MTM LENGTH  [ 0=short; 1=long ]
C	CHAN	I*4	CHANNEL     [ 1; 2; 3; 4 ]
C	NGROUP	I*4	ADDS PER GROUP
C	SCIMODE	I*4	SCIENCE (MICROPROCESSOR) MODE  [ 0; 1; 2; 3; 4 ]
C	FIRSTSAMP I*4	FIRSTSAMP [ start sample no.; normally=1 ]
C	NYQUIST RECORD  Reference Dataset of Nyquist Frequencies
C
C      OUTPUT:
C	STARTX	R*4	X-AXIS START POINT (cm)
C	SPACEX	R*4	X-AXIS SPACE INTERVAL (cm)
C
C      SUBROUTINES CALLED: FUT_DEFAULT_PEAK
C
C      INCLUDE FILES: NONE
C---------------------------------------------------------------------------
C Changes:
C
C	R. Kummerer, December 5, 1986. STARTX position adjusted by -SPACEX
C	for coadd plots. SPACEX was calculated assuming noise data was
C	being plotted.
C
C	R. Kummerer, September 15, 1987. Use FUT_NYQUIST to calculate STARTX
C	and SPACEX.
C
C	R. Kummerer, October 28, 1987. Noisetest spectra are always in Hz.
C	FUT_NYQUIST was asking for wavenumbers.
C
C	R. Kummerer, December 29, 1987. Check when fake-it data is being
C	plotted and request Hz with FUT_NYQUIST.  SPR 1834.
C
C	F. Shuman, 1988 May 31.  Repair logic of checking both fake-it and
C	iflag to determine spacex and startx.  In collaboration with
C	R. Isaacman and R. Kummerer.
C
C	R. Kummerer, August 25, 1988. SPR 2434, Error status return.
C
C	R. Kummerer, March 3, 1989. SPR 3274, Initialize status return to
C	FUT_NORMAL.
C
C	R. Kummerer, August 14, 1989. SPR 4317, Correct offset term for
C	calculating MTM distance relative to the IFG peak.  Changed from
C	1.19 to 1.20.
C
C	F. Shuman, 1990 Oct 30.  SPR 7490, Include offset for microproc's
C	start-sampling counter when calculating the position of the zero-
C	path-difference point.  Inserted a new calling argument, OFFSET.
C
C	F. Shuman, 1990 Nov 09.  SPR 7494, Make MSPEED consistent with other
C	FIRAS s/w & document the meaning of calling arguments to FUT_NYQUIST.
C
C	F. Shuman, Hughes STX, 1992 Mar 04.  SPR 9514, Correct the zero-point
C	on plots and double SPACEX now that our spectra are 256 pts, not 512.
C	Renamed misleading OFFSET to FIRSTSAMP.
C
C	S. Alexander, Hughes STX, 1992 July 14.  SER 9790, Use new reference
C	dataset FEX_NYQUIST for Nyquist frequencies instead of FUT routine
C	FUT_NYQUIST.
C---------------------------------------------------------------------------

	Implicit None

C  Calling arguments:
	Integer*4	iflag		! coadd or spectrum flag
					!  -1=noisetest spectrum
					!   0=coadd
					!   1=spectrum
	Integer*4	fake_it		! fake-it flag:  0=real; 1=fake-it
	Integer*4	mspeed		! speed of mtm scan:   0=slow;  1=fast
	Integer*4	mlen		! length of mtm scan:  0=short; 1=long
	Integer*4	chan		! channel: 1=RH; 2=RL; 3=LH; 4=LL
	Integer*4	ngroup		! adds per group
	Integer*4	scimode		! science (microproc) mode  [0,1,2,3,4]
	Integer*4	linflag /0/	! (IFG is not linearized for plots)
	Integer*4	firstsamp	! 1st sample # for non-standard sweeps
	Real   *4	startx		! x-axis starting point
	Real   *4	spacex		! x-axis discretization

 	Dictionary 'fex_nyquist'
	Record /fex_nyquist/ nyquist    ! reference dataset record structure
C  Other:
	Integer*4	peak		! IFG pk sample # from FUT_Default_Peak
	Real   *4	fnyq_icm	! Nyquist frequency in icm
	Real   *4	fnyq_hz		! Nyquist frequency in hertz
        Integer*4	freq,scan_mode
        Integer*4       status
C  Routines:
	Integer*4	FUT_Default_Peak
C  Externals:
	External	FUT_Normal, FUT_PlotScale
        External	FUT_Fakeithz, FUT_Illadds
        External	FUT_Illchan, FUT_Illmode
        External	FUT_Illsmp

	status = %Loc(FUT_Normal)
C
C Retrieve Nyquist frequency in icm and hertz
C
	If (fake_it .eq. 1) Then
C
C Can only find Nyquist frequency in hertz for fake-it data.
C
	   If (ngroup .eq. 1) Then
              fnyq_hz = nyquist.hz(9)
           Else If (ngroup .eq. 2) Then
              fnyq_hz = nyquist.hz(10)
           Else If (ngroup .eq. 3) Then
              fnyq_hz = nyquist.hz(11)
           Else If (ngroup .eq. 8) Then
              fnyq_hz = nyquist.hz(12)
           Else If (ngroup .eq. 12) Then
              fnyq_hz = nyquist.hz(13)
           Else
C
C Illegal adds-per-group value input.
C
              fut_setxax = %loc(fut_illadds)
              call lib$signal(fut_illadds)
              Return
           End If
        Else If (fake_it .eq. 0) Then
C
C Determine high or low frequency and scan mode to get Nyquist frequency.
C
           If ((chan .eq. 1) .or. (chan .eq. 3)) Then
              freq = 0
           Else If ((chan .eq. 2) .or. (chan .eq. 4)) Then
              freq = 4
           Else
C
C Illegal channel value input.
C
              fut_setxax = %loc(fut_illchan)
              call lib$signal(fut_illchan)
              Return
           End If

           scan_mode = mlen*2 + mspeed + 1

           If ((scan_mode .lt. 1) .or. (scan_mode .gt. 4)) Then
C
C Illegal scan mode value input.
C
              fut_setxax = %loc(fut_illmode)
              call lib$signal(fut_illmode)
              Return
           End If

           fnyq_icm = nyquist.icm(freq + scan_mode)
           fnyq_hz =  nyquist.hz(freq + scan_mode)

        Else
C
C Illegal fake-it value input.
C
           fut_setxax = %loc(fut_illsmp)
           call lib$signal(fut_illsmp)
           Return
        End If         

	If (iflag .Eq. 0) Then
C
C Dealing with distance
C
	   If (fake_it .Eq. 1) Then
C
C ...Fake-it data
C
	      spacex = 1./512
	      startx = 0.
	   Else
C
C ...Real data; 
C
	      spacex = 0.5/fnyq_icm
	      status = FUT_Default_Peak(mspeed, mlen, chan, ngroup, scimode,
     &                                  linflag, peak)
	      startx = ((firstsamp-1.)/ngroup - peak + 1)*spacex
	   End If

	Else If (iflag .Eq. 1) Then
C
C Dealing in wave number space
C
	   If (fake_it .Eq. 1) Then
C
C ...Fake-it data; fails because Nyquist frequency can only be specified
C    in hertz for fake-it data.
C
              fut_setxax = %loc(fut_fakeithz)
              call lib$signal(fut_fakeithz)
              Return
	   Else
C
C ...Real data; 
C
	      spacex = fnyq_icm/256
	      startx = 0.
	   End If

	Else If (iflag .Eq. -1) Then
C
C Dealing with NOISETEST spectra; 
C
	   spacex = fnyq_hz/256
	   startx = 0.
	Else
C
C Unknown type flag
C
	   status = %Loc(FUT_PlotScale)
	   Call LIB$Signal(FUT_PlotScale)
	End If

	FUT_SetxAx = status

	Return
	End
