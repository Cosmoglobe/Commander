	integer * 4 function  frd_dgtl_transient_fcn (flag, response)

c-------------------------------------------------------------------------------
c
c	Function FRD_DGTL_TRANSIENT_FCN
c
c	This function calculates the transient response functions for the Firas
c	digital filters.  Either the high or low frequency response function is
c	calculated.
c
c	Adapted from an algorithm developed by Rich Isaacman.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  11 September 1991, SER 7985
c
c-------------------------------------------------------------------------------
c
c	Input:
c		flag		integer * 4	hi/lo frequency flag
c
c	Output:
c		response(512)	real 	* 4	digital filter transient 
c						response function
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	integer * 4	flag		!  hi/lo freuqency flag
	integer * 4	j		!  a counter

	real	* 4	response(512)	!  transient response function
	real	* 4	X(516)		!  filter input signal
	real	* 4	Y(516)		!  filter response

	external	frd_normal


c
c  Initialize the frequency response arrays.
c
	call lib$movc5(0,,0,2048,response)
	call lib$movc5(0,,0,2064,X)
	call lib$movc5(0,,0,2064,Y)
	do j = 4,516
	   X(j) = 1.0
	enddo


	if (flag .eq. 1) then
c
c  Calculate the transient high frequency response function.
c
	   do j = 5,516
	      Y(j) = (X(j) + 2.0*X(j-1) + X(j-2) +
     .		      11.0*Y(j-1) - 5.0*Y(j-2)) / 8.0
	      response(j-4) = 1.0 - Y(j)/2.0
	   enddo

	elseif (flag .eq. 2) then
c
c  Calculate the transient low frequency response function.
c
	   do j = 5,516
	      Y(j) = (X(j) + 2.0*X(j-1) + 2.0*X(j-2) + 2.0*X(j-3) +
     .		      X(j-4) + 8.0*Y(j-2) - Y(j-4)) / 8.0
	      response(j-4) = 1.0 - Y(j)/8.0
	   enddo

	endif


	frd_dgtl_transient_fcn = %loc(frd_normal)

	return
	end
