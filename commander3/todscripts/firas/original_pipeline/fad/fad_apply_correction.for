	Integer*4  Function  FAD_Apply_Correction  ( time, spec, variance,
	1                                            fex_ejv_rec )

c------------------------------------------------------------------------------
c   Purpose: Determine appropriate destriper correction and apply to spectrum.
c
c   Algorithm:
c      Correction (f) =  Sum_ufo_i A(f)_i * exp (-t/tau_i) +
c                        Sum_top_j B(f)_j * Top(t)_j +
c                        (V0 + V1*ty + V2*ty^2 + V3*ty^3 + V4*ty^4) * Vspec(f)
c      Corrected Spectra (f) = Spectra (f) - Correction (f)
c      where: f = frequency index, t is time in days since cover ejection
c      (893251118), ty = time in years since cover eject.  taus, As, and Bs
c      are in the reference file FEX_EJV_ccss.
c
c   Input Parameters:
c      integer*4  time(2)       -- time of the spectrum in adt
c      complex*8  spec(257)     -- spectrum to destripe
c      real*4     variance(257) -- real variance of spectrum
c      record /fex_ejv/    fex_ejv_rec    -- reference record
c
c   Output Parameters:
c      complex*8  spec(257)     -- destriped spectrum
c
c   Include Files:
c      FAD_msg       -- External message params needed for FAD
c      FUT_params    -- FIRAS parameters
c
c   Functions and Subroutines:
c      Lib$Signal
c      Time_GE
c      Lib$Subx
c      Lib$Ediv
c
c   Author:  Larry Paris Rosen, Hughes STX, 13 May 1993
c   Modified: L. Rosen, May 1994.  Remove all gain correction references.
c   Modified: L. Rosen, June 1994.  Correction now has several ufo's and tau's,
c      and vibration correction (quartic).  The corrections are summed first
c      so if they are small, then we will be summing small numbers and then
c      later subtracting from the spectrum value.  This is better numerically.
c   Modified: L. Rosen, July 29, 1994.  Spectrum should be zero at frequencies
c      where the variance is zero because data there is completely
c      undetermined.   SER  118847.
c------------------------------------------------------------------------------
	Implicit None

c Include files:

	Include		'(FAD_msg)'
	Include		'(FUT_params)'

c Passed Parameters:

	Integer*4	time (2)
	Complex*8	spec (257)
	Real*4		variance (257)
	Dictionary	'FEX_EJV'
	Record /FEX_EJV/	fex_ejv_rec

c Function
	Logical*1	Time_GE

c Local
	Integer*2	num_periods
	Parameter	(num_periods = 13)
	Integer*2	tim			! time period
	Integer*4	start_time (2)		! period start time
	Integer*4	stop_time (2)		! period stop time
	Integer*2	top (13)	! top hat = 1 if t in period else 0
	Integer*4	diff_adt (2)		! adt time - time_eject
	Integer*4	adt_to_secs		! diff_adt / adt_to_secs = secs
	Parameter	(adt_to_secs = 10000000)
	Integer*4	timesecs		! time in secs since eject
	Integer*4	remainder
	Real*8		timedays		! time in days since eject
	Real*8		timeyears		! time in years since eject
	Real*8		vib			! vibration correction
	Complex*16	correction (257)	! correction to spectrum
	Integer*2	specsize
	parameter	(specsize = 4112)
	Byte		zero_corr (specsize) / specsize * 0 /
						! initializes correction to 0
	Integer*2	f			! frequency index
	Integer*2	i			! counter
	Integer*2	tauindex		! tau number
	Integer*2	topindex		! top hat number
c------------------------------------------------------------------------------

c Begin

	FAD_Apply_Correction = %Loc (FAD_Normal)

c Set the top hats. 1 if time is in the period, else 0.

	Do tim = 1, num_periods
	   start_time (1) = fex_ejv_rec.Start_Times (1,tim)
	   start_time (2) = fex_ejv_rec.Start_Times (2,tim)
	   stop_time (1) = fex_ejv_rec.Stop_Times (1,tim)
	   stop_time (2) = fex_ejv_rec.Stop_Times (2,tim)
	   If (Time_GE (time, start_time) .AND. Time_GE (stop_time, time)) Then
	      top (tim) = 1
	   Else
	      top (tim) = 0
	   Endif
	Enddo

c Calculate time since apco cover eject in days and years.

	Call Lib$Subx (time, fex_ejv_rec.eject_time, diff_adt)
	Call Lib$Ediv (adt_to_secs, diff_adt, timesecs, remainder)
	timedays = timesecs / fac_day
	timeyears = timedays / 365.

c Calculate vibration correction - quartic in time in years.

	vib = fex_ejv_rec.vib_corr (1) +
	1     timeyears * (fex_ejv_rec.vib_corr (2) +
	2     (timeyears * (fex_ejv_rec.vib_corr (3) +
	3     (timeyears * (fex_ejv_rec.vib_corr (4) +
	4     (timeyears * fex_ejv_rec.vib_corr (5)))))))

c Initialize correction to 0

	Call Lib$MovC3 ( specsize, zero_corr, correction )

c Build up correction.  corr_index is used to determine correction type:
c 1 = exponential decay (ufo), 2 = vibration correction, 3 = top hat.

	tauindex = 1				! initialize tau number
	topindex = 1				! initialize tophat number
	i = 1
	Do While (i .LE. num_periods .AND. fex_ejv_rec.corr_index (i) .NE. 0)
	   If (fex_ejv_rec.corr_index (i) .EQ. 1) Then		! UFO
	      Do f = 1, 257
	         correction (f) = correction (f) + fex_ejv_rec.corr_spec (f,i)
	1                       * Exp (- timedays / fex_ejv_rec.Tau(tauindex))
	      Enddo
	      tauindex = tauindex + 1
	   Elseif (fex_ejv_rec.corr_index (i) .EQ. 2) Then	! Vibration
	      Do f = 1, 257
	         correction (f) = correction (f) + fex_ejv_rec.corr_spec (f,i)
	1                         * vib
	      Enddo
	   Elseif (fex_ejv_rec.corr_index (i) .EQ. 3) Then	! Tophat
	      If (top (topindex) .NE. 0) Then
	         Do f = 1, 257
	            correction (f) = correction (f) +
	1                            fex_ejv_rec.corr_spec (f,i)
	2                            * top (topindex)
	         Enddo
	      Endif
	      topindex = topindex + 1
	   Endif
	   i = i + 1
	Enddo
	Do f = 1, 257
           If (variance (f) .EQ. 0.) Then
	      spec (f) = cmplx (0.,0.)
	   Else
	      spec (f) = spec (f) - correction (f)
	   Endif
	Enddo
	Return
	End
