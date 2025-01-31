	integer * 4 function  frd_apod_fcn (ictr, good_peak, apod_fcn)

c-------------------------------------------------------------------------------
c
c	Function FRD_APOD_FCN
c
c	This routine generates a real*8, 512-point apodization function whose
c	shape depends only on the position of the IFG peak.  The routine first
c	forms a function in which a symmetric region about the peak, equal in
c	length to the short side of the stroke, is "averaged", whereas the
c	remaining long part of the stroke is given full weight.  The routine
c	then multiplies this function by a smoothing function, which is an
c	inverted octic gentered at the peak position and tapering to zero at the
c	end of the long side of the stroke.
c
c	The routine zeros out the first two points of the IFG to compensate for
c	the effects of the digital transient response subtraction.
c
c	This routine is based on the routine FRD_SHAPE written by Rich Isaacman,
c	but it has been modified to reflect the changes made in the algorithm
c	made by Dale Fixsen.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  24 January 1992
c		  SER 8836
c
c-------------------------------------------------------------------------------
c
c	Input:
c		ictr		integer * 4		ifg peak position
c		good_peak	integer * 4		peak position flag
c
c	Output:
c		apod_fcn(512)	real	* 8		apodization function
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	integer * 4	good_peak	!  peak position flag
	integer * 4	ictr		!  ifg peak position
	integer * 4	j		!  a counter
	integer * 4	k		!  a counter

	real	* 8	apod_fcn(512)	!  apodization function
	real	* 8	pi		!  pi
	parameter	(pi = 3.141592653589793)
	real	* 8	x		!  apodization function scale factor

	external	frd_normal

c
c  Initialize the apodization function.
c
	do j = 1, 2
	   apod_fcn(j) = 0.0D0
	enddo
	do j = 3, 512
	   apod_fcn(j) = 1.0D0
	enddo


c
c  Average the apodization function according to the peak position.
c

	if (.not. good_peak) then
c
c  If the peak position is bad, set it to 512 and do not average the function.
c
	   if (ictr .ge. 513) then
	      ictr = 512
	   endif
	   x = dble(ictr - 2)

	else
c
c  Average the function.
c
	   if (ictr .lt. 257) then
	      x = dble(513 - ictr)
	      do j = 1, 30
	         apod_fcn(2+j) = 0.5D0 * (1.0D0 - dcos(dble(j)*pi/30.0D0))
	      enddo
	      do j = 1, ictr-2
	         apod_fcn(ictr+j) = 2.0D0 - apod_fcn(ictr-j)
	      enddo
	      do j = 2*ictr-1, 512
	         apod_fcn(j) = 2.0D0
	      enddo
	   else
	      x = dble(ictr - 1)
	      do j = 1, 30
	         apod_fcn(513-j) = 0.5D0 * (1.0D0 - dcos(dble(j)*pi/30.0D0))
	      enddo
	      do j = 1, 512-ictr
	         apod_fcn(ictr-j) = 2.0D0 - apod_fcn(ictr+j)
	      enddo
	      do j = 3, 2*ictr-513
	         apod_fcn(j) = 2.0D0
	      enddo
	   endif

	endif		!  good_peak


C
C  Smooth the apodization function.
C
	do j = 1, 512
	   apod_fcn(j) = apod_fcn(j) * (1.0D0 - (dble(j-ictr)/x)**4)**2
	enddo


	frd_apod_fcn = %loc(frd_normal)

	return
	end
