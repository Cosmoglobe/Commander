	integer * 4 function  fcf_compute_constants ()

c-------------------------------------------------------------------------------
c
c	Function FCF_COMPUTE_CONSTANTS
c
c	This function computes the constants required by FCF to process
c	spectra.  The constants computed are:
c
c	1)  the Nyquist frequency and frequency interval in icm and hz
c	2)  the phase shift for high frequency short fast spectra
c	3)  the normalization for the FFTed spectra
c
c	Author:   Gene Eplee
c		  General Sciences Corporation
c		  513-7768
c		  26 June 1992
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		none
c
c	Include files:
c		fcf_config.txt
c		fcf_invoc.txt
c		fcf_model.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c       SER 11395
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fcf_config)'
	include '(fcf_invoc)'
	include '(fcf_model)'

	integer * 4	frec		!  Nyquist frequency record number
	integer * 4	j		!  a counter

	external	fcf_normal


C
C  Compute the constants.
C

c
c  Get the Nyquist frequencies and frequency intervals in icm and hz.
c	Fill in the frequency arrays.
c
	if (fcc_xllf .eq. fac_present) then
	   frec     = 4*jmod((fcc_chan-1),2) + 2
	else
	   frec     = 4*jmod((fcc_chan-1),2) + fcc_smode
	endif
	fnyq_icm = dble(config.nyquist.icm(frec))
	df       = fnyq_icm / 256.0D0
	fnyq_hz  = dble(config.nyquist.hz(frec))
	dw       = fac_dpi * fnyq_hz / 128.0D0		! 2 * pi * f / 256.0

	do j = 1,257
	   freq(j)  = dble(j-1) * df
	   afreq(j) = dble(j-1) * dw
	enddo

c
c  Calculate the phase shift required to make the high frequency short fast
c	spectra compatable with the high frequency long fast spectra.
c
	if (fcc_xhsf .eq. fac_present) then
	   phase_shift = dcmplx(0.0D0,fac_dpi/512.0D0)
	else
	   phase_shift = dcmplx(0.0D0,0.0D0)
	endif

c
c  Calculate the normalization for the FFTed spectra.
c	The normalization factor is the product of the Nyquist frequency,
c	the instrument throughput, and the adc voltage scale factor.  The signs
c	of the spectra for the right side channels is flipped so the right side
c	transfer functions will be positive.
c
	if (fcc_chan .le. 2) then
           spec_norm = - fnyq_icm * dble(fac_etendu) * dble(fac_adc_scale)
	elseif (fcc_chan .ge. 3) then
           spec_norm = fnyq_icm * dble(fac_etendu) * dble(fac_adc_scale)
	endif


	fcf_compute_constants = %loc(fcf_normal)

	return
	end
