	integer * 4 function  ffi_compute_constants ()

c-------------------------------------------------------------------------------
c
c	Function FFI_COMPUTE_CONSTANTS
c
c	This function computes the constants required by FFI to process
c	spectra.  The constants computed are:
c
c	1)  the Nyquist frequency in icm
c	2)  the phase shift for high frequency short fast spectra
c	3)  the renormalization for the voltage variances
c	4)  the ifg peak position
c
c	Author:   Gene Eplee
c		  General Sciences Corporation
c		  513-7768
c		  9 March 1993
c		  SER 10763
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
c		ffi_config.txt
c		ffi_invoc.txt
c		ffi_spec.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11690
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(ffi_config)'
	include '(ffi_invoc)'
	include '(ffi_spec)'

	integer * 4	frec			!  Nyquist frequency record
						!    number
	integer * 4	j			!  a counter
	integer * 4	linearized/0/		!  linearization flag
	integer * 4	rstatus			!  return status
	integer * 4	upmode/4/		!  microprocessor mode

	real	*  8	fnyq_icm		!  Nyquist frequency in icm

	integer * 4	fut_default_peak

	external	ffi_normal


C
C  Compute the constants.
C

c
c  Get the Nyquist frequencies and frequency intervals in icm and hz.
c
	if (fcc_xllf .eq. fac_present) then
	   frec  = 4*jmod((fcc_chan-1),2) + 4
	else
	   frec  = 4*jmod((fcc_chan-1),2) + fcc_smode
	endif
	fnyq_icm = dble(config.nyquist.icm(frec))

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
c  Calculate the renormalization for the voltage variances.
c
	var_norm = fnyq_icm * fac_etendu

c
c  Determine the ifg peak position for hybrid dataset determination.
c
	if (fcc_xllf .eq. fac_present) then
	   rstatus = fut_default_peak (1, 1, fcc_chan, 8, upmode, linearized,
     .				       peak_pos)
	else
	   rstatus = fut_default_peak (fcc_speed, fcc_length, fcc_chan,
     .				       fcc_ngroup, upmode, linearized, peak_pos)
	endif


	ffi_compute_constants = %loc(ffi_normal)

	return
	end
