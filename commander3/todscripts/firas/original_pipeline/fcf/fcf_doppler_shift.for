	integer * 4 function  fcf_doppler_shift (vbary, cspec)

c-------------------------------------------------------------------------------
c
c	Function FCF_DOPPLER_SHIFT
c
c	This function Doppler shifts a calibrated spectrum from the spacecraft
c	frame of reference into the barycentric frame.  The frequency shift is
c	performed by a cubic spline interpolation, applied to the real and 
c	imaginary parts of the calibrated spectrum separately.  The frequency
c	scaling factor is (1 + v/c), while the photometric scaling factor is
c	(1 - v/c)**3.  Only first-order effects are considered since
c	|v/c| << 0.001.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  26 June 1992
c
c-------------------------------------------------------------------------------
c
c	Input:
c		vbary		real 	*  8		projection of the
c							spacecraft barycentric
c							velocity along the Firas
c							line of sight, in km/sec
c		cspec		complex * 16		calibrated spectrum
c
c	Output:
c		cspec		complex * 16		Doppler-shifted
c							calibrated spectrum
c
c	Subroutines called:
c		dcsint
c		dcsval
c		lib$movc5
c
c	Include files:
c		fcf_model.txt
c
c-------------------------------------------------------------------------------
c
c	Hard-coded constants:
c
c		VLIGHT		2.99792458D+05  km/sec		speed of light
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Restrict the spline fit for the Doppler shift to the frequency range
c	where the spectra are actually good.
c	Gene Eplee, GSC, 1 August 1994
c	SPR 11850
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fcf_model)'

	complex * 16	cspec(257)	!  calibrated spectrum

	integer *  4	j		!  a counter
	integer *  4	npts		!  number of good data points in
					!    input spectrum
	integer *  4	offset		!  good data spectrum offset

	logical *  1	first_time	!  status flag

	real	*  8	beta_minus	!  (1 - v/c)**3
	real	*  8	beta_plus	!  (1 + v/c)
	real	*  8	freqsp(257)	!  spline frequencies
	real	*  8	ibreak(257)	!  imaginary spline break points
	real	*  8	idopspec(257)	!  imaginary Doppler-shifted spectrum
	real	*  8	ispec(257)	!  imaginary part of calibrated spectrum
	real	*  8	ispline(4,257)	!  imaginary spline coefficients
	real	*  8	rbreak(257)	!  real spline break points
	real	*  8	rdopspec(257)	!  real Doppler-shifted spectrum
	real	*  8	rspec(257)	!  real part of calibrated spectrum
	real	*  8	rspline(4,257)	!  real spline coefficients
	real	*  8	vbary		!  projected barycentric velocity
	real	*  8	vlight		!  speed of light in km/sec
	parameter	(vlight = 2.99792458D+05)

	real	*  8	dcsval

	external	fcf_normal

C
C  Perform the Doppler shift
C

c
c  Extract the real and imaginary parts of the calibrated spectrum.
c     Loop over frequencies, checking for zeros in the transfer function.
c
	call lib$movc5 (0,,0,2056,freqsp)
	call lib$movc5 (0,,0,2056,rspec)
	call lib$movc5 (0,,0,2056,ispec)
	npts = 0
	first_time = .true.
	do j = 1, 257
	   if (cdabs(emiss(1,j)) .gt. trlim) then
	      if (first_time) then
	         first_time = .false.
	         offset = j - 1
	      endif
	      npts = npts + 1
	      freqsp(npts) = freq(j)
	      rspec(npts)  = dreal(cspec(j))
	      ispec(npts)  = dimag(cspec(j))
	   endif
	enddo

c
c  Calculate the spline interpolation coefficients.
c
	call lib$movc5 (0,,0,2056,rbreak)
	call lib$movc5 (0,,0,2056,ibreak)
	call lib$movc5 (0,,0,8224,rspline)
	call lib$movc5 (0,,0,8224,ispline)
	call dcsint (npts, freqsp, rspec, rbreak, rspline)
	call dcsint (npts, freqsp, ispec, ibreak, ispline)

c
c  Calculate the frequency points in the barycentric frame.  The Doppler-shifted
c	frequency array consists of the original frequencies scaled by
c	(1 + v/c).
c
	call lib$movc5 (0,,0,2056,freqsp)
	beta_plus  = 1.0D0 + vbary/vlight
	do j = 1, npts
	   freqsp(j) = freq(j+offset) * beta_plus
	   if (freqsp(j) .ge. freq(257)) freqsp(j) = freq(257)
	enddo

c
c  Interpolate the calibrated spectrum onto the barycentric frequency array.
c
	call lib$movc5 (0,,0,2056,rdopspec)
	call lib$movc5 (0,,0,2056,idopspec)
	do j = 1, npts
	   rdopspec(j) = dcsval (freqsp(j), npts-1, rbreak, rspline)
	   idopspec(j) = dcsval (freqsp(j), npts-1, ibreak, ispline)
	enddo

c
c  Rescale the Doppler-shifted spectrum by the photometric scaling factor of
c	(1 - v/c)**3.  Write the scaled, shifted spectrum back into the
c	calibrated spectrum array.
c
	call lib$movc5 (0,,0,2056,rspec)
	call lib$movc5 (0,,0,2056,ispec)
	call lib$movc5 (0,,0,4112,cspec)
	beta_minus = (1.0D0 - vbary/vlight)**3
	do j = 1, npts
	   rspec(j+offset) = rdopspec(j) * beta_minus
	   ispec(j+offset) = idopspec(j) * beta_minus
	   cspec(j+offset) = dcmplx(rspec(j+offset),ispec(j+offset))
	enddo


	fcf_doppler_shift = %loc(fcf_normal)

	return
	end
