	integer * 4 function  fsl_compute_constants ()

c-------------------------------------------------------------------------------
c
c	Function FSL_COMPUTE_CONSTANTS
c
c	This function computes the constants required by FSL to process
c	spectra.  The constants computed are:
c
c	1)  the Nyquist frequency and frequency interval in icm and hz
c	2)  the phase shift for high frequency short fast spectra
c	3)  the normalization for the FFTed spectra
c
c	Author:   
c                FCF_Compute_Constants
c                Gene Eplee
c		 General Sciences Corporation
c		 26 June 1992
c      
c                FSL_Compute_Constants
c                Shirley M. Read
c                Hughes STX Corporation
c                July 1995
c
c-------------------------------------------------------------------------------
c
c	Input:	none
c
c	Output: none
c
c	Include files:
c		fsl_config.txt
c		fsl_invoc.txt
c		fsl_model.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes for FCF:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c       SER 11395
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, July 14, 1995 
c       Modified FCF_Compute_Constants to FSL_Compute_Constants for the new 
c       FIRAS pipeline which will process long spectra to get improved 
c       frequency resolution.
c         1. Change hard-coded spectra length and FFT length to fcc_spec_length
c            which varies for scan modes. Spectra and FFT lengths are set from 
c            the FUT_Params include file in the FSL_Parse routine.
c         2. Separate Nyquist frequencies for low frequency channels in LF 
c            scan mode have been created in the FEX_Nyquistl reference data set.
c            The use of the fcc_xllf flag is no longer needed and was thus 
c            removed from the code.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fsl_config)'
	include '(fsl_invoc)'
	include '(fsl_model)'

	integer * 4	frec		!  Nyquist frequency record number
	integer * 4	j		!  a counter
        real * 8        spec_len        !  spectra length in double precision
        real * 8        fft_len         !  FFT length in double precision

	external	fsl_normal


C
C  Compute the constants.
C

c
c  Get the Nyquist frequencies and frequency intervals in icm and hz.
c	Fill in the frequency arrays.
c
	frec     = 4*jmod((fcc_chan-1),2) + fcc_smode
	fnyq_icm = dble(config.nyquist.icm(frec))
        spec_len = dfloat(fcc_spec_length - 1)
	df       = fnyq_icm / spec_len
	fnyq_hz  = dble(config.nyquist.hz(frec))
	dw       = 2.0 * fac_dpi * fnyq_hz / spec_len  ! 2 * pi * f / spec_len

	do j = 1,fcc_spec_length
	   freq(j)  = dble(j-1) * df
	   afreq(j) = dble(j-1) * dw
	enddo

c
c  Calculate the phase shift required to make the high frequency short fast
c	spectra compatable with the high frequency long fast spectra.
c
	if (fcc_xhsf .eq. fac_present) then
           fft_len  = dfloat(fcc_fft_length)
	   phase_shift = dcmplx(0.0D0,fac_dpi/fft_len)
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


	fsl_compute_constants = %loc(fsl_normal)

	return
	end
