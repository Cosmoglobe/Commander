	integer * 4 function  ffp_galactic_cut (pixel_no)

c-------------------------------------------------------------------------------
c
c	Function FFP_GALACTIC_CUT
c
c	This function filters pixels by galactic latitude.  For each input
c	pixel number it determines the latitude of the galactic center and
c	compares that to the cutoff latitude.  If the pixel latitude is greater
c	than the cutoff latitude, the pixel number is advanced by one and the
c	check is done again.
c
c	Author:  S. Brodd, HSTX, 3/21/96
c
c	Input:
c		pixel_no	integer * 4		input pixel number
c
c	Output:
c		pixel_no	integer * 4		output pixel number
c
c	Subroutines called:
c		firas_cenpix
c		xcc_e_to_g
c
c	Include files:
c		ffp_invoc_sky.txt
c		fut_params.txt
c
c	Modifications:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(ffp_invoc_sky)'

	integer * 4	pixel_no
	integer * 4	status

	real	* 4	glat
	real	* 4	evec(3)
	real	* 4	gvec(3)

	external	ffp_eof
	external	ffp_normal
c
c  Determine the galactic latitude of the input pixel number.
c
	call firas_cenpix (pixel_no, evec)
	call xcc_e_to_g (evec, fac_epoch, gvec)
	glat = atan2d(gvec(3),sqrt(gvec(1)**2 + gvec(2)**2))
c
c  Filter excluded galactic latitudes.
c
	do while ((abs(glat) .gt. fcc_glat)  .and.
     .		  (pixel_no .lt. 6143))
	   pixel_no = pixel_no + 1
	   call firas_cenpix (pixel_no, evec)
	   call xcc_e_to_g (evec, fac_epoch, gvec)
	   glat = atan2d(gvec(3),sqrt(gvec(1)**2 + gvec(2)**2))
	enddo
c
c  Check the status of the output pixel number.
c
	if ((pixel_no .lt. 6143)  .or.
     .	    ((pixel_no .eq. 6143) .and. (abs(glat) .le. fcc_glat))) then
	   status = %loc(ffp_normal)
	else
	   status = %loc(ffp_eof)
	endif

	ffp_galactic_cut = status

	return
	end
