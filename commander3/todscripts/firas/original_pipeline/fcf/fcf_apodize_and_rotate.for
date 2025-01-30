	integer * 4 function fcf_apodize_and_rotate (aifg, ivar, rifg, vvar)

c-------------------------------------------------------------------------------
c
c	Function FCF_APODIZE_AND_ROTATE
c
c	This function apodizes and rotates coadded IFGs, using apodization
c	functions read from the binary reference dataset FEX_APOD.  It also
c	apodizes voltage variances by convolving them with FFT'ed apodization
c	functions.
c
c	Author: Gene Eplee
c		General Sciences Corp.
c		513-7768
c		2 November 1992
c
c-------------------------------------------------------------------------------
c
c	Input:
c		aifg(512)	real* 8		input ifg
c		ivar(3,257)	real * 8	unapodized voltage variances
c
c	Output:
c		rifg(512)	real* 8		apodized, rotated ifg
c		vvar(3,257)	real * 8	apodized voltage variances
c
c	Subroutines called:
c		lib$movc5
c
c	Include files:
c		fcf_config.txt
c		fcf_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11395
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fcf_config)'
	include '(fcf_invoc)'

	integer * 4	i		!  a counter
	integer * 4	j		!  a counter
	integer * 4	k		!  a counter

	real 	* 8	aifg(512)	!  coadded ifg
	real	* 8	ivar(3,257)	!  unapdoized voltage variances
	real	* 8	rifg(512)	!  apodized, rotated ifg
	real	* 8	uvar(3,512)	!  unpacked voltage variance
	real	* 8	vvar(3,257)	!  apodized voltage variances

	external	fcf_normal

c
c  Initialize the apodized, rotated ifg and the variance arrays.
c
	call lib$movc5 (0,,0,4096,rifg)
	call lib$movc5 (0,,0,12288,uvar)
	call lib$movc5 (0,,0,6168,vvar)

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

	if (fcc_single .eq. fac_present) then
c
c  Set the voltage variances to flag values.
c
	   do j = 1,3
	      do k = 1,257
	         vvar(j,k) = fcc_varflag
	      enddo
	   enddo

	else
	   if ((fcc_xlsf .eq. fac_present)  .or.
     .	       (fcc_xllf .eq. fac_present)) then
c
c  Unpack the 64-point voltage variances.
c
	      do i = 1,3
	         do j = 2,64
	            uvar(i,j) = ivar(i,j)
	            uvar(i,130-j) = ivar(i,j)
	         enddo
	         uvar(i,65) = ivar(i,65)
	      enddo

c
c  Convolve the unpacked variances with the FFT'ed apodization function to
c	apodize the variances.
c
	      do i = 1,3
	         do j = 1,65
	            do k = 1,128
	               vvar(i,j) = vvar(i,j) + uvar(i,k) * zapod_fcn(129+k-j)
	            enddo
	         enddo
	      enddo

	   else
c
c  Unpack the 256-point voltage variances.
c
	      do i = 1,3
	         do j = 2,256
	            uvar(i,j) = ivar(i,j)
	            uvar(i,514-j) = ivar(i,j)
	         enddo
	         uvar(i,257) = ivar(i,257)
	      enddo

c
c  Convolve the unpacked variances with the FFT'ed apodization function to
c	apodize the variances.
c
	      do i = 1,3
	         do j = 1,257
	            do k = 1,512
	               vvar(i,j) = vvar(i,j) + uvar(i,k) * zapod_fcn(513+k-j)
	            enddo
	         enddo
	      enddo

	   endif	!	(xllf or xlsf

	endif		! 	(fcc_single


	fcf_apodize_and_rotate = %loc(fcf_normal)

	return
	end
