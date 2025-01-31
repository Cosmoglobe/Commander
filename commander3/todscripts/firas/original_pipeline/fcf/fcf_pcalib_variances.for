	integer * 4 function  fcf_pcalib_variances (vvar, cvar)

c-------------------------------------------------------------------------------
c
c	Function FCF_PCALIB_VARIANCES
c
c	This function partially applies the Firas calibration model solution to
c	the voltage variances, producing calibrated variances with the units of
c	(ergs/sec/cm**2/sr/icm)**2.  The input parameters are the voltage
c	variances, the detector responsivity and delay, and the instrument
c	transfer function.  The real-real, imaginary-imaginary, and
c	real-imaginary variances are all calibrated.
c
c	Author:  Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 23 September 1992
c
c-------------------------------------------------------------------------------
c
c	Input:
c		vvar(3,257)	real	*  8		voltage variances
c							 ( R*R, I*I, R*I )
c
c	Output:
c		cvar(3,257)	real * 8		calibrated variances
c							 ( R*R, I*I, R*I )
c
c	Include files:
c		fcf_model.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fcf_model)'

	complex * 16	pshift		!  phase shift of calibrated variances
	complex * 16	xfer	 	!  calibration function

	integer *  4	j		!  a counter

	real	*  8	cvar(3,257)	!  calibrated variances
	real	*  8	iphase		!  imaginary phase shift
	real	*  8	iphase2		!  square of imaginary phase shift
	real	*  8	rphase		!  real phase shift
	real	*  8	rphase2		!  square of real phase shift
	real	*  8	vvar(3,257)	!  voltage variances
	real	*  8	xfer2		!  square of calibration function

	external	fcf_normal

c
c  Initialize the arrays.
c
	call lib$movc5 (0,,0,6168,cvar)


C
C  Calibrate the variances.
C

c
c  Loop over frequencies, checking for zeroa in the transfer function.
c
	do j = 2, 257

	   if (cdabs(emiss(1,j)) .gt. trlim) then
c
c  Calculate the calibration normalization and phase shifts.
c
	      xfer    = B(j) / emiss(1,j) / S0
	      xfer2   = dreal(xfer * dconjg(xfer))
	      pshift  = xfer / cdabs(xfer)
	      rphase  = dreal(pshift)
	      rphase2 = rphase**2.0
	      iphase  = - dimag(pshift)
	      iphase2 = iphase**2.0

c
c  Calibrate the variances (real-real, imaginary-imaginary, and real-imaginary).
c
	      cvar(1,j) =  (rphase2*vvar(1,j) - 2.0*rphase*iphase*vvar(3,j)
     .			  + iphase2*vvar(2,j)) * xfer2
	      cvar(2,j) =  (rphase2*vvar(2,j) + 2.0*rphase*iphase*vvar(3,j)
     .			  + iphase2*vvar(1,j)) * xfer2
	      cvar(3,j) = ((iphase2-rphase2)*vvar(3,j)
     .			  - rphase*iphase*(vvar(1,j) - vvar(2,j))) * xfer2

	   endif	!  check for zeros in the transfer function

	enddo		!  loop over frequencies


	fcf_pcalib_variances = %loc(fcf_normal)

	return
	end
