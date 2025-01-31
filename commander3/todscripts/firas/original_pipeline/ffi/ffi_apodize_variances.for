	integer * 4 function  ffi_apodize_variances (vvar, vsigma)

c-------------------------------------------------------------------------------
c
c	Function FFI_APODIZE_VARIANCES
c
c	This function apodizes the real-real voltage variances, then computes
c	the real-real apodized voltage sigmas.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  24 September 1992
c		  SER 10763
c
c-------------------------------------------------------------------------------
c
c	Input:
c		vvar(3,257)	real * 8	voltage variances
c
c	Output:
c		vsigma(257)	real * 8	apodized real-real voltage
c						sigmas
c
c	Subroutines called:
c		lib$movc5
c
c	Include files:
c		ffi_config.txt
c		ffi_spec.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(ffi_config)'
	include '(ffi_spec)'

	integer * 4	j			!  a counter
	integer * 4	k			!  a counter

	real	* 8	cvar(257)		!  zapodized voltage variance
	real	* 8	rvar(512)		!  unpacked rr voltage variance
	real	* 8	vsigma(257)		!  apodized rr voltage sigmas
	real	* 8	vvar(3,257)		!  voltage variances

	external	ffi_normal

c
c  Initilize the variance arrays.
c
	call lib$movc5 (0,,0,2056,cvar)
	call lib$movc5 (0,,0,4096,rvar)

c
c  Unpack the real-real variance into the RVAR vector.
c
	rvar(1) = vvar(1,1)
	do j = 2,256
	   rvar(j) = vvar(1,j)
	   rvar(514-j) = vvar(1,j)
	enddo
	rvar(257) = vvar(1,257)

c
c  Convolve the RVAR vector with the zapodization function to produce the
c	zapodized variance.  Produce and renormalize the voltage sigma.
c
	do j = 1,257
	   do k = 1,512
	      cvar(j) = cvar(j) + rvar(k) * zapod_fcn(513+k-j)
	   enddo
	   vsigma(j) = dsqrt(cvar(j)) * var_norm
	enddo


	ffi_apodize_variances = %loc(ffi_normal)

	return
	end
