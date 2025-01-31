	integer * 4 function  fip_galactic_cut (pixel_no)

c-------------------------------------------------------------------------------
c
c	Function FIP_GALACTIC_CUT
c
c	This function filters pixels by galactic latitude.  For each input
c	pixel number it determines the latitude of the galactic center and
c	compares that to the cutoff latitude.  If the pixel latitude is greater
c	than the cutoff latitude, the pixel number is advanced by one and the
c	check is done again.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  10 May 1993
c
c-------------------------------------------------------------------------------
c
c	Input:
c		pixel_no	integer * 4		input pixel number
c
c	Output:
c		pixel_no	integer * 4		output pixel number
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		firas_cenpix
c		xcc_e_to_g
c
c	Include files:
c		fip_invoc_sky.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fip_invoc_sky)'

	integer * 4	pixel_no
	integer * 4	status

	real	* 4	glat
	real	* 4	evec(3)
	real	* 4	gvec(3)

	external	fip_eof
	external	fip_normal

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
	   status = %loc(fip_normal)
	else
	   status = %loc(fip_eof)
	endif


	fip_galactic_cut = status

	return
	end
