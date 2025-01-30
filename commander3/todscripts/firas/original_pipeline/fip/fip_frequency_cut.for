	integer * 4 function  fip_frequency_cut (fnyq_icm, fnyq_hz)

c-------------------------------------------------------------------------------
c
c	Function FIP_FREQUENCY_CUT
c
c	This function determines the initial and final frequency indices of the
c	FIRAS Initial Product datasets.  The routine uses the Nyquist frequency
c	of the scan mode in icm, the low frequency cutoff in icm, and the high 
c	frequency cutoff in icm to compute the initial and final frequency
c       indices and the number of frequency points.
c
c	For the optical frequencies in GHz, the routine computes the initial
c	frequency point and the frequency interval for entry into the FIP
c	datasets from the Nyquist frequency in icm.  For the data frequencies
c	in radians/sec, the routine computes the initial frequency point and
c	frequency interval in radians/sec for entry into the calibration model
c	initial products from the Nyquist freqeuncy in hz.
c
c	Note:  The only check that FIP_Frequency_Cut performs on the requested
c	frequency range is that the lower index is greater than or equal to one
c	and that the upper index is less than or equal to 180.
c
c	FIP_FREQUENCY_CUT DOES NOT VERIFY THAT THE SPECTRA OR EMISSIVITIES ARE
c	NON-ZERO OVER THE SPECIFIED FREQUENCY RANGE.  THIS DETERMINATION IS
c	LEFT TO THE USER BEFORE THE PROGRAM IS RUN.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  25 May 1993
c
c-------------------------------------------------------------------------------
c
c	Input:
c		fnyq_icm	real * 4		Nyquist frequency in icm
c		fnyq_hz		real * 4		Nyquist frequency in hz
c
c	Output:
c		none
c
c	Include files:
c		fip_frequency.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fip_frequency)'

	real	* 4	df		!  frequency interval in icm
	real	* 4	fnyq_hz		!  Nyquist frequency in hz
	real	* 4	fnyq_icm	!  Nyquist frequency in icm

	external	fip_normal

c
c  Compute the frequency interval in icm, radians/sec, and GHz.
c
	df      = fnyq_icm/256.0
	fcc_dw  = 2.0*fac_dpi*(fnyq_hz/256.0)
	fcc_dnu = df * fac_icm_ghz

	if (fcc_freq .eq. fac_present) then
c
c  Compute the lower and upper frequency indices
c
	   fcc_jlo = jint(fcc_lofreq/df) + 2
	   if (fcc_jlo .lt. 1) fcc_jlo = 1
	   fcc_jhi = jint(fcc_hifreq/df) + 1
	   if (fcc_jhi .gt. 180) fcc_jhi = 180
	else
	   fcc_jlo    = fcc_jlo_default
	   fcc_jhi    = fcc_jhi_default
	endif

c
c Compute the frequency index shift, the first frequency point, and the number
c	of frequency points.
c
	fcc_freq_offset = fcc_jlo - 1
	fcc_nu0         = float(fcc_freq_offset) * fcc_dnu
	fcc_w0          = float(fcc_freq_offset) * fcc_dw
	fcc_nfreq       = fcc_jhi - fcc_jlo + 1

	fip_frequency_cut = %loc(fip_normal)

	return
	end
