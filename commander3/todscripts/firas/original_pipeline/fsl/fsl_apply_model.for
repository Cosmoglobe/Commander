	integer * 4 function  fsl_apply_model (primary_vib, nifgs, adjnifgs, 
     .                        temp, tsig, vspec, cvar, power, phase_corr, cspec)

c-------------------------------------------------------------------------------
c
c	Function FSL_APPLY_MODEL
c
c	This function applies the Firas calibration model to a voltage spectrum
c	with the units of volts/cm**2/sr/icm, producing a calibrated spectrum
c	with the units of ergs/sec/cm**2/sr/icm.  The routine calibrates the
c	differential spectrum and produces a reference spectrum.  It then
c	computes a linear autophase corrector for the spectrum.  Finally, it
c	calibrates the spectrum.
c
c	The input parameters are the instrument temperatures, voltage spectrum,
c	the detector responsivity and delay, the primary vibration correction, 
c	and the instrument emissivities.  The output parameters are the total
c	IR power incident on the bolometer, the phase corrector, and the
c	calibrated spectrum.
c
c	Author:  
c                FCF_Apply_Model
c                Gene Eplee
c		 General Sciences Corp.
c		 11 May 1993
c
c                FSL_Apply_Model
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
c
c-------------------------------------------------------------------------------
c
c	Input:
c		vspec(361)	complex * 16		voltage spectrum
c               cvar(3,361)     real    *  8            calibrated variances
c		nifgs		integer *  4		number of ifgs in the
c							  spectrum
c               adjnifgs        real    *  4            glitch rate adjusted
c                                                         number of ifgs
c		primary_vib	real	*  8		primary (time-dependent)
c							  vibration corrrection
c		temp(10)	real 	*  8		instrument temperatures
c		tsig(10)	real	*  8		temperature sigmas
c
c	Output:
c		power		real	*  8		total IR power
c		phase_corr	real	*  8		linear phase corrector
c		cspec(361)	complex * 16		calibrated spectrum
c
c	Subroutines called:
c		dplanck_dist (uoe_dplanck_dist)
c			dplanck_dist calls dplanck (uoe_dplanck)
c		fsl_autophase_correct
c		lib$movc5
c
c	Include files:
c		fsl_config.txt
c		fsl_display.txt
c		fsl_invoc.txt
c		fsl_model.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Hard-coded constants (contained in UOE_DPLANCK):
c
c		HCK		1.438769 D+00  cm K		second radiation
c								   constant
c		THCC		1.1910439D-05  erg cm^2		first radiation
c								   constant
c
c	Hard-coded constants (contained in UOE_DPLANCK_DIST):
c
c		SQ2		1.414213562373095D+00		sqrt(2.0)
c
c-------------------------------------------------------------------------------
c
c	Changes for FCF:
c
c	Modified vibration correction to separate short fast and long fast
c	scan modes.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11395
c
c	Added D-Vector option to the determination of the variance estimate for
c	the autophase correction.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11397
c
c	Insert actual autophase-corrected spectrum into DISPLAY structure for
c	use by FCF_DISPLAY_SPECTRA.
c	Gene Eplee, GSC, 31 March 1994
c	SPR 11681
c
c	Put in optional differential spectrum output and add reference spectrum
c	to DISPLAY structure.
c	Gene Eplee, GSC,11 July 1994
c	SPR 11826
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, August 10, 1995 
c       Modified FCF_Apply_Model to FSL_Apply_Model for the new FIRAS pipeline 
c       which will process long spectra to get improved frequency resolution.
c           1. Changed status, include file, and function names for FSL.
c           2. Changed array lengths to accomodate long spectra. 
c           3. Computed new number of bytes to initialize arrays to zero.
c           4. Used start and stop frequency indices corresponding to model
c              to process arrays.
c           5. Check for flag to use FIL averaged variances to weight the 
c              autophase corrector computation. If it is set, allow the variance
c              vector to be set from the real calibrated variance vector 
c              regardless of number of IFGs, compute the autophase correctors,
c              and calibrate the spectrum.
c           6. Added glitch rate adjusted number of IFGs to calling sequence for
c              weighting the Dvector.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fsl_config)'
	include '(fsl_display)'
	include '(fsl_invoc)'
	include '(fsl_model)'

	complex * 16	Bvp		!  detector delay at ifp
	complex * 16	Bvs		!  detector delay at ifs
	complex * 16	cspec(361)	!  calibrated spectrum
	complex * 16	element(7)	!  optical model spectral components
	complex * 16	pspec(361)	!  reference spectrum (sum of
					!    emissivities * Planck function)
	complex * 16	spec(361)	!  harmonic, vibration corrected
					!    spectrum sample with detector
					!    model applied
	complex * 16	vspec(361)	!  uncalibrated voltage spectrum

	integer *  4	fp		!  frequency offset for primary
					!    vibration correction
	integer *  4	fs		!  frequency offset for secondary
					!    vibration correction
	integer *  4	ifp		!  frequency counter for primary
					!    vibration correction
	integer *  4	ifs		!  frequency counter for secondary
					!    vibration correction
	integer *  4	if2a		!  frequency counter for 2nd harmonic
	integer *  4	if2b		!  frequency counter for 2nd harmonic
	integer *  4	if3a		!  frequency counter for 3rd harmonic
	integer *  4	if3b		!  frequency counter for 3rd harmonic
	integer *  4	if3c		!  frequency counter for 3rd harmonic
	integer *  4	it		!  a temperature index
	integer *  4	iv		!  a vibration index
	integer *  4	j		!  a counter
	integer *  4	k		!  a counter
	integer *  4	nifgs		!  number of ifgs in the spectrum
	integer *  4	status		!  return status

	real    *  4    adjnifgs        !  glitch rate adjusted number of ifgs

	real	*  8	cvar(3,361)	!  calibrated spectrum variances
	real	*  8	phase_corr	!  linear phase corrector
	real	*  8	power		!  total IR power on the bolometer
	real	*  8	primary_vib	!  primary (time-dependent) vibration
					!    correction
	real	*  8	temp(10)	!  instrument temperatures
	real	*  8	tsig(10)	!  temperature sigmas
	real	*  8	variances(361)	!  variance estimate

	integer *  4	fsl_autophase_correct
	real	*  8	dplanck_dist

	external	fsl_normal

C
C  Initialize the arrays and the processing parameters.
C
	status = %loc(fsl_normal)
	call lib$movc5 (0,,0,5776,spec)
	call lib$movc5 (0,,0,5776,pspec)
	call lib$movc5 (0,,0,5776,cspec)
	call lib$movc5 (0,,0,2888,variances)
	it    = fcc_chan + 6
	iv    = 3*(fcc_chan-1) + (fcc_speed+1) + fcc_length
	power = 0.0D0
	if (fcc_sky .eq. fac_present) then
	   temp(1) = sky_temp
	   tsig(1) = sky_temp_sig
	endif


C
C  Calibrate the differential spectrum.
C  

c
c  Loop over frequencies, checking for zeros in the transfer function.
c
	do j = start_idx,stop_idx

	   if (cdabs(emiss(1,j)) .gt. trlim) then
c
c  Calculate the frequencies required for computing the vibration correction
c	detector delays.
c
	      fp = (j-1) - config.vibcorr.primary_offset(iv)
	      fs = (j-1) - config.vibcorr.secondary_offset(iv)
	      ifp = jmax0(jiabs(fp),1)
	      ifs = jmax0(jiabs(fs),1)
c
c  Calculate the vibration correction detector delays.
c
	      Bvp  = dcmplx(1.0D0, dble(jiabs(fp))*dw*tau)
	      Bvs  = dcmplx(1.0D0, dble(jiabs(fs))*dw*tau)

c
c  Calculate the frequencies required for the 2nd and 3rd harmonic corrections.
c
	      if2a = (j-1)/2
	      if2b = j/2
	      if3a = (j-1)/3
	      if3b = j/3
	      if3c = (j+1)/3
	      if2a = jmax0(if2a,1)
	      if2b = jmax0(if2b,1)
	      if3a = jmax0(if3a,1)
	      if3b = jmax0(if3b,1)
	      if3c = jmax0(if3c,1)

c
c  Convert the spectrum from volts/cm**2/sr/icm to ergs/sec/cm**2/sr/icm.
c  Correct the spectrum for the 3rd and 2nd harmonic terms and and for
c	the secondary and primary (time-dependent) vibration terms.
c
	      spec(j) = (B(j) * vspec(j)
     .			   - param(10) * (2.0D0 + B(j))/9.0D0 *
     .			    (vspec(if3a) + vspec(if3b) + vspec(if3c))
     .			   - param(11) * (1.0D0 + B(j))/4.0D0 *
     .			    (vspec(if2a) + vspec(if2b))
     .			   - param(15) * dble(fs) * Bvs * vspec(ifs)
     .			   - primary_vib * dble(fp) * Bvp * vspec(ifp))
     .			/ (S0 * emiss(1,j))
	      display.espec(j) = cmplx(spec(j))

	   endif	!  check for zeros in the transfer function
	enddo		!  loop over frequencies


C
C  Compute the reference spectrum (sum of emissivities * Planck function) and
C	the IR power incident on the bolometer.  Determine the variance estimate
C	for the autophase correction.
C

c
c  Loop over frequencies, checking for zeros in the transfer function.
c
	do j = start_idx,stop_idx

	   if (cdabs(emiss(1,j)) .gt. trlim) then
c
c  Compute the optical model contributions to the reference spectrum.
c
	      call lib$movc5 (0,,0,112,element)
	      do k = 1,6
	         element(k) = (emiss(k,j)/emiss(1,j)) *
     .			       dplanck_dist(freq(j),temp(k),tsig(k))
	      enddo
	      element(7)  = (emiss(7,j)/emiss(1,j)) *
     .			     dplanck_dist(freq(j),temp(it),tsig(it))

c
c  Compute the reference spectrum and the total IR power.
c
	      power = power + cdabs(element(1))
	      do k = 2,7
	         pspec(j) = pspec(j) + element(k)
	         power = power + cdabs(element(k))
	      enddo
	      display.pspec(j) = cmplx(pspec(j))

c
c  Determine the variance estimate. Weight Dvector by glitch rate adjusted
c  number of IFGs.
c
	      if (fcc_dvec .eq. fac_present) then
	         variances(j) = dvector(j)/dble(adjnifgs)
	      elseif ((nifgs .ge. 3) .or. (fcc_flv .eq. fac_present)) then
	         variances(j) = cvar(1,j)
	      endif

	   endif	!  check for zeros in the transfer function
	enddo		!  loop over frequencies
	power = rscale * power


	if ((fcc_dvec .eq. fac_present) .or. (nifgs .ge. 3) .or. 
     .     (fcc_flv .eq. fac_present)) then
c
c  Compute the linear autophase corrector.
c
	   status = fsl_autophase_correct (spec, pspec, variances, phase_corr)
	endif


C
C  Calibrate the spectrum.
C
	if (status .eq. %loc(fsl_normal)) then
	   if ((fcc_dvec .eq. fac_present) .or. (nifgs .ge. 3) .or.
     .        (fcc_flv .eq. fac_present)) then
c
c  Autophase correct and calibrate the spectrum.
c
	      do j=start_idx,stop_idx
	         if (cdabs(emiss(1,j)) .gt. trlim) then
	            if (fcc_diff .eq. fac_present) then
	               cspec(j) = spec(j) *
     .				      cdexp(dcmplx(0.0D0,dble(j-1)*phase_corr))
	            else
	               cspec(j) = spec(j) *
     .				      cdexp(dcmplx(0.0D0,dble(j-1)*phase_corr))
     .			        - pspec(j)
	            endif
	         endif
	      enddo
	   else
c
c  Calibrate the spectrum.
c
	      do j=start_idx,stop_idx
	         if (cdabs(emiss(1,j)) .gt. trlim) then
	            if (fcc_diff .eq. fac_present) then
	               cspec(j) = spec(j)
	            else
	               cspec(j) = spec(j) - pspec(j)
	            endif
	         endif
	      enddo
	   endif
	endif		!  status from autophase


	fsl_apply_model = status
	
	return
	end
