	integer * 4 function ffi_apodize_and_rotate (aifg, vvar, rifg, vsigma)

c-------------------------------------------------------------------------------
c
c	Function FFI_APODIZE_AND_ROTATE
c
c	This function apodizes and rotates coadded IFGs, using apodization
c	functions read from the binary reference dataset FEX_APOD.  It also
c	apodizes the real-real voltage variances by convolving them with
c	FFT'ed apodization functions, then computes the real-real voltage
c	simgas.
c
c	Author: Gene Eplee
c		General Sciences Corp.
c		513-7768
c		2 November 1992
c		SER 10763
c
c-------------------------------------------------------------------------------
c
c	Input:
c		aifg(512)	real	* 8		input ifg
c		vvar(3,257)	real	* 8		voltage variances
c
c	Output:
c		rifg(512)	real	* 8		apodized, rotated ifg
c		vsigma(257)	real	* 8		apodized real-real
c							voltage sigmas
c
c	Subroutines called:
c		lib$movc5
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

	integer * 4	j		!  a counter
	integer * 4	k		!  a counter

	real 	* 8	aifg(512)	!  coadded ifg
	real	* 8	cvar(257)	!  zapodized voltage variance
	real	* 8	rifg(512)	!  apodized, rotated ifg
	real	* 8	rvar(512)	!  unpacked rr voltage variance
	real	* 8	vsigma(257)	!  apodized rr voltage sigmas
	real	* 8	vvar(3,257)	!  voltage variances

	external	ffi_normal

c
c  Initialize the apodized, rotated ifg and the variance arrays.
c
	call lib$movc5 (0,,0,2056,cvar)
	call lib$movc5 (0,,0,4096,rifg)
	call lib$movc5 (0,,0,4096,rvar)

	if ((fcc_xlsf .eq. fac_present)  .or.
     .	    (fcc_xllf .eq. fac_present)) then
c
c  Apodize and rotate the 128-point IFG.
c
	   do j = 1,128
	      k = j + 1 - ictr
	      if (k .le. 0) k = k + 128
	      rifg(k) = aifg(j) * apod_fcn(j)
	   enddo
	else
c
c  Apodize and rotate the 512-point IFG.
c
	   do j = 1,512
	      k = j + 1 - ictr
	      if (k .le. 0) k = k + 512
	      rifg(k) = aifg(j) * apod_fcn(j)
	   enddo
	endif


	if ((fcc_xlsf .eq. fac_present)  .or.
     .	    (fcc_xllf .eq. fac_present)) then

c
c  Unpack the 64-point real-real variances into the RVAR vector.
c
	   do j = 2,64
	      rvar(j) = vvar(1,j)
	      rvar(130-j) = vvar(1,j)
	   enddo
	   rvar(65) = vvar(1,65)

c
c  Convolve the real-real variances with the FFT'ed apodization function to
c	apodize the variances.  Produce and renormalize the voltage sigmas.
c
	   do j = 1,65
	      do k = 1,128
	         cvar(j) = cvar(j) + rvar(k) * zapod_fcn(129+k-j)
	      enddo
	      vsigma(j) = dsqrt(cvar(j)) * var_norm
	   enddo

	else
c
c  Unpack the 256-point real-real variances into the RVAR vector.
c
	   do j = 2,256
	      rvar(j) = vvar(1,j)
	      rvar(514-j) = vvar(1,j)
	   enddo
	   rvar(257) = vvar(1,257)

c
c  Convolve the real-real variances with the FFT'ed apodization function to
c	apodize the variances.  Produce and renormalize the voltage sigmas.
c
	   do j = 1,257
	      do k = 1,512
	         cvar(j) = cvar(j) + rvar(k) * zapod_fcn(513+k-j)
	      enddo
	      vsigma(j) = dsqrt(cvar(j)) * var_norm
	   enddo

	endif


	ffi_apodize_and_rotate = %loc(ffi_normal)

	return
	end
