	integer * 4 function  frd_apod_fcnl (ictr, good_peak, short, xhsf,
     .			xhlf, apod_fcn)
c-------------------------------------------------------------------------------
c
c	Function FRD_APOD_FCNL
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
c-------------------------------------------------------------------------------
c
c	Changes:  Added "short" flag and formation of length 128 apodization
c       functions to aid in processing of special low fast data in which
c       interferograms have either been truncated or smoothed and decimated.
c       Changed flag for bad peak if "short" and peak bigger than 128.
c       Alice Trenholme, GSC, 7 February 1995
c       
c          Changed passed peak position to real * 4 to accomodate HSF data
c
c          Changed apodization for HSF data so that
c              1.  Passed peak position is shifted DOWN 1/2 bin if 
c                  nonlinearized peak position is used and microprocessor 
c                  is on
c              2.  Passed peak position is shifted UP 1/2 bin otherwise
c              3.  Apodization is symmetric about the fractional peak position
c              4.  First four points are zeroed
c          Changed apodization for HLF data so that
c              1.  Points 1, 510, 511, and 512 are zeroed
c              2.  Point 2 is not zeroed
c                        Alice Trenholme, GSC, 15 February 1995, SER 12244
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

	implicit none

	integer * 4	good_peak	!  peak position flag
	real    * 4	ictr		!  ifg peak position
	integer * 4	j		!  a counter
	integer * 4	k		!  a counter
        integer * 4     short           !  0 = make a length 512 function
                                        !  1 = make a length 128 function
        integer * 4     xhsf            !  High short fast indicator
        integer * 4     xhlf            !  High long fast indicator
        integer * 4     upperlim        !  last nonzero data point
        integer * 4     halflim         !  65 if short=1, 257 otherwise
        integer * 4     lowerlim        !  first nonzero data point
        integer * 4     twidth         !  transition length; 
                                        !  30 if short = 0, 8 if short = 1
        real    * 8     dtwidth        !  double precision twidth
	real	* 8	apod_fcn(512)	!  apodization function
	real	* 8	pi		!  pi
	parameter	(pi = 3.141592653589793)
	real	* 8	x		!  apodization function scale factor

	external	frd_normal

c
c  Initialize the apodization function.
c
        halflim = 257
        twidth = 30
        upperlim = 512
        lowerlim = 3

        if (short) then
		upperlim = 128
		lowerlim = 3
                halflim = 65
                twidth = 8
                if (ictr .gt. 128) good_peak = 0
	else if (xhsf) then
		lowerlim = 5
	else if (xhlf) then
		upperlim = 509
		lowerlim = 2
	endif

	do j = 1, lowerlim - 1
	   apod_fcn(j) = 0.0D0
	enddo
	do j = lowerlim, upperlim
	   apod_fcn(j) = 1.0D0
	enddo
        do j = 1 + upperlim, 512
	   apod_fcn(j) = 0.0D0
        enddo

        dtwidth = dble(twidth)
c
c  Average the apodization function according to the peak position.
c
	if (.not. good_peak) then
c
c  If the peak position is bad, set it to upperlim and do not average 
c      the function.
c
	   if (ictr .ge. 1 + upperlim) then
	      ictr = upperlim
	   endif
	   x = dble(ictr - 2)

	else
c
c  Average the function.
c
	   if (ictr .lt. halflim) then
	      x = dble(1 + upperlim - ictr)
	      do j = 1, twidth - 1
	         apod_fcn(j + lowerlim - 1) = 
     .                  0.5D0 * (1.0D0 - dcos(dble(j)*pi/DTWIDTH))
                 apod_fcn(j + 2*ictr - lowerlim - twidth + 1) =
     .                  0.5D0 * (3.0D0 - dcos(dble(j)*pi/DTWIDTH))
	      enddo

              do j = 2*ictr - lowerlim + 1, upperlim 
		 apod_fcn(j) = 2.d0
              enddo

	   else

	      x = dble(ictr + 1 - lowerlim)

	      do j = 1, twidth - 1
		apod_fcn(j + 2*ictr - upperlim - 1) =
     .                  0.5D0 * (3.0D0 + dcos(dble(j)*pi/DTWIDTH))
		apod_fcn(j + upperlim - twidth + 1) = 
     .                  0.5D0 * (1.0D0 + dcos(dble(j)*pi/DTWIDTH))
	      enddo

	      do j = lowerlim, 2*ictr - upperlim - 1
	         apod_fcn(j) = 2.0D0
	      enddo

	   endif

	endif		!  good_peak

C
C  Smooth the apodization function.
C
	do j = lowerlim, upperlim
	   apod_fcn(j) = apod_fcn(j) * (1.0D0 - (dble(j-ictr)/x)**4)**2
	enddo


	frd_apod_fcnl = %loc(frd_normal)

	return
	end
