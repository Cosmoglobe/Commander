	integer * 4 function  ffp_frequency_cut (fnyq_icm, fnyq_hz)

c-------------------------------------------------------------------------------
c
c	Function FFP_FREQUENCY_CUT
c
c	For the optical frequencies in GHz, the routine computes the initial
c	frequency point and the frequency interval for entry into the FFP
c	datasets from the Nyquist frequency in icm.  For the data frequencies
c	in radians/sec, the routine computes the initial frequency point and
c	frequency interval in radians/sec for entry into the calibration model
c	products from the Nyquist freqeuncy in hz.
c
c	Note:  The only check that FFP_Frequency_Cut performs on the requested
c	frequency range is that the lower index is greater than or equal to one
c	and that the upper index is less than or equal to 210.
c
c	FFP_FREQUENCY_CUT DOES NOT VERIFY THAT THE SPECTRA OR EMISSIVITIES ARE
c	NON-ZERO OVER THE SPECIFIED FREQUENCY RANGE.  THIS DETERMINATION IS
c	LEFT TO THE USER BEFORE THE PROGRAM IS RUN.
c
c	Author:  S. Brodd, HSTX, 3/21/96
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
c		ffp_frequency.txt
c		fut_params.txt
c
c	Modifications:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(ffp_frequency)'

	real	* 4	df		!  frequency interval in icm
	real	* 4	fnyq_hz		!  Nyquist frequency in hz
	real	* 4	fnyq_icm	!  Nyquist frequency in icm

	external	ffp_normal
c
c  Compute the frequency interval in icm, radians/sec, and GHz.
c
	df      = (fnyq_icm * fac_nyq_correct) / (fac_spec_length(fcc_fsmode)-1)
	fcc_dw  = 2.0*fac_dpi*(fnyq_hz/(fac_spec_length(fcc_fsmode) - 1))
	fcc_dnu = df * fac_icm_ghz
c
c Compute the frequency index shift, the first frequency point, and the number
c	of frequency points.
c
	fcc_freq_offset = fcc_jlo - 1
	fcc_nu0         = float(fcc_freq_offset) * fcc_dnu
	fcc_w0          = float(fcc_freq_offset) * fcc_dw
	fcc_nfreq       = fcc_jhi - fcc_jlo + 1

	ffp_frequency_cut = %loc(ffp_normal)

	return
	end
