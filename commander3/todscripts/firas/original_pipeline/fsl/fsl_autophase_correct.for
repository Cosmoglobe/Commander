	integer * 4 function  fsl_autophase_correct (spec, pspec, variances,
     .						     phase_corr)

c-------------------------------------------------------------------------------
c
c	Function FSL_AUTOPHASE_CORRECT
c
c	This function computes a linear autophase corrector for the calibrated
c	spectrum.
c
c	Author:   
c                FCF_Autophase_Correct
c                Gene Eplee
c		 General Sciences Corp.
c		 17 May 1993
c
c                FSL_Autophase_Correct
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
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
c		fsl_model.txt
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
c	Changes for FCF:
c
c	Changed weights from the calibrated variances CVAR to a more general
c	variance estimate VARIANCES.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11397
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, August 8, 1995 
c       Modified FCF_Autophase_Correct to FSL_Autophase_Correct for the new 
c       FIRAS pipeline which will process long spectra to get improved 
c       frequency resolution.
c           1. Changed status, include file, and function names for FSL.
c           2. Changed array lengths to accomodate long spectra. 
c           3. Used start and stop frequency indices corresponding to model
c              to process arrays.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fsl_model)'

	complex	* 16	pspec(361)	!  reference spectrum
	complex	* 16	spec(361)	!  calibrated differential spectrum

	integer *  4	j		!  a counter

	real	*  8	alpha		!  constraint on phase corrector
	parameter	(alpha = 1.0D+06)
	real	*  8	f		!  frequency bin
	real	*  8	phase_corr	!  phase corrector
	real	*  8	variances(361)	!  variance estimate
	real	*  8	x		!  real part of phase error
	real	*  8	y		!  imaginary part of phase error

	external	fsl_normal

c
c  Initialize the calculation.
c
	x = 0.0D0
	y = 0.0D0
	phase_corr = 0.0D0

c
c  Compute the autophase corrector.
c
	do j = start_idx,stop_idx
	   if (cdabs(emiss(1,j)) .gt. trlim) then
	      f = dble(j-1)
	      x = x + f*f*(dreal(spec(j)))**2/variances(j)
	      y = y + f*dimag(spec(j)-pspec(j))*dreal(spec(j))/variances(j)
	   endif
	enddo
	phase_corr = - y / (x + alpha)


	fsl_autophase_correct = %loc(fsl_normal)
	
	return
	end
