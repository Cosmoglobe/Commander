	integer * 4 function  fcf_autophase_correct (spec, pspec, variances,
     .						     phase_corr)

c-------------------------------------------------------------------------------
c
c	Function FCF_AUTOPHASE_CORRECT
c
c	This function computes a linear autophase corrector for the calibrated
c	spectrum.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  17 May 1993
c
c-------------------------------------------------------------------------------
c
c	Input:
c		spec		complex * 16	calibrated differential spectrum
c		pspec		complex * 16	reference spectrum spectrum
c		variances	real	*  8	variance estimate
c
c	Output:
c		phase_corr	real	*  8		linear phase corrector
c
c-------------------------------------------------------------------------------
c
c	Include files:
c		fcf_model.txt
c
c-------------------------------------------------------------------------------
c
c	Hard-coded constants:
c
c		ALPHA		1.0D+06			constraint on phase
c							corrector
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Changed weights from the calibrated variances CVAR to a more general
c	variance estimate VARIANCES.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11397
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fcf_model)'

	complex	* 16	pspec(257)	!  reference spectrum
	complex	* 16	spec(257)	!  calibrated differential spectrum

	integer *  4	j		!  a counter

	real	*  8	alpha		!  constraint on phase corrector
	parameter	(alpha = 1.0D+06)
	real	*  8	f		!  frequency bin
	real	*  8	phase_corr	!  phase corrector
	real	*  8	variances(257)	!  variance estimate
	real	*  8	x		!  real part of phase error
	real	*  8	y		!  imaginary part of phase error

	external	fcf_normal

c
c  Initialize the calculation.
c
	x = 0.0D0
	y = 0.0D0
	phase_corr = 0.0D0

c
c  Compute the autophase corrector.
c
	do j = 2,257
	   if (cdabs(emiss(1,j)) .gt. trlim) then
	      f = dble(j-1)
	      x = x + f*f*(dreal(spec(j)))**2/variances(j)
	      y = y + f*dimag(spec(j)-pspec(j))*dreal(spec(j))/variances(j)
	   endif
	enddo
	phase_corr = - y / (x + alpha)


	fcf_autophase_correct = %loc(fcf_normal)
	
	return
	end
