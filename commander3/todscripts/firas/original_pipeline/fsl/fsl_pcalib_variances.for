	integer * 4 function  fsl_pcalib_variances (vvar, cvar)

c-------------------------------------------------------------------------------
c
c	Function FSL_PCALIB_VARIANCES
c
c	This function partially applies the Firas calibration model solution to
c	the voltage variances, producing calibrated variances with the units of
c	(ergs/sec/cm**2/sr/icm)**2.  The input parameters are the voltage
c	variances, the detector responsivity and delay, and the instrument
c	transfer function.  The real-real, imaginary-imaginary, and
c	real-imaginary variances are all calibrated.
c
c	Author:  
c                FCF_Pcalib_Variances
c                Gene Eplee
c		 General Sciences Corp.
c		 23 September 1992
c
c                FSL_Pcalib_Variances
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
c
c-------------------------------------------------------------------------------
c
c	Input:
c		vvar(3,361)	real	*  8		voltage variances
c							 ( R*R, I*I, R*I )
c
c	Output:
c		cvar(3,361)	real * 8		calibrated variances
c							 ( R*R, I*I, R*I )
c
c	Include files:
c		fsl_model.txt
c
c-------------------------------------------------------------------------------
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, August 9, 1995 
c       Modified FCF_Pcalib_Variances to FSL_Pcalib_Variances for the new 
c       FIRAS pipeline which will process long spectra to get improved 
c       frequency resolution.
c           1. Changed status, include file, and function names for FSL.
c           2. Changed array lengths to accomodate long spectra. 
c           3. Computed new number of bytes to initialize arrays to zero.
c           4. Used start and stop frequency indices corresponding to model
c              to process arrays.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fsl_model)'

	complex * 16	pshift		!  phase shift of calibrated variances
	complex * 16	xfer	 	!  calibration function

	integer *  4	j		!  a counter

	real	*  8	cvar(3,361)	!  calibrated variances
	real	*  8	iphase		!  imaginary phase shift
	real	*  8	iphase2		!  square of imaginary phase shift
	real	*  8	rphase		!  real phase shift
	real	*  8	rphase2		!  square of real phase shift
	real	*  8	vvar(3,361)	!  voltage variances
	real	*  8	xfer2		!  square of calibration function

	external	fsl_normal

c
c  Initialize the arrays.
c
	call lib$movc5 (0,,0,8664,cvar)

C
C  Calibrate the variances.
C

c
c  Loop over frequencies, checking for zeroa in the transfer function.
c
	do j = start_idx, stop_idx

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


	fsl_pcalib_variances = %loc(fsl_normal)

	return
	end
