	subroutine fut_ttinterp (ttag1, ttag2, target, coefs)
c----------------------------------------------------------------------
c  This routine returns the interpolation coefficients at time 
c  tag "target" for data associated with time tags "ttag1" and "ttag2".
c  In other words, the weight to be applied to the data from time "ttag1"
c  is:
c
c			  ttag2 - target
c	       coef(1) =  --------------
c			  ttag2 - ttag1
c
c  and the weight for the data from time ttag2 is coef(2) = 1 - coef(1).
c  Obviously, if the "target" time tag is midway between ttag1 and ttag2 
c  then coef(1) = coef(2) = 0.5. 
c
c  Time tags ttag1, ttag2, and target are supplied as ADT binary times.
c  Coefficients are returned as real*4. It doesn't matter whether ttag1 
c  or ttag2 is the later time, but "target" must lie in between them; 
c  otherwise coef(1) = coef(2) = 0. is returned.
c-------------------------------------------------------------------------
c  Calling sequence:
c
c	call ttinterp 
c			(ttag1,		! integer*4 array: length 2 (ADT)
c			 ttag2,		! integer*4 array: length 2 (ADT)
c			 target,	! integer*4 array: length 2 (ADT)
c			 coefs)		! real*4 array: length 2 (interp coeffs)
c-----------------------------------------------------------------------------
c  Subroutines called:	lib$day		! converts ADT into day number + units
c					!   of 0.01-second units since midnight
c-----------------------------------------------------------------------------
c  Written by:		Rich Isaacman
c			Applied Research Corp.
c			21 Dec 1987
c-----------------------------------------------------------------------------
	implicit none
	integer*4 	ttag1(2)
	integer*4	ttag2(2)
	integer*4	target(2)
	integer*4	numday1, numday2, num_tgt
	integer*4	hunsec1, hunsec2, hunsec_tgt
	real*4		coefs(2)
	real*8		diff21
	real*8		difft1
	real*8		hsec_per_day

	parameter	(hsec_per_day = 86400.d+02)	!0.01-sec units in a day

	call lib$day (numday1, ttag1, hunsec1)
	call lib$day (numday2, ttag2, hunsec2)
	call lib$day (num_tgt, target, hunsec_tgt)
c
c  Calculate the intervals ttag2-ttag1 (= diff21) and target-ttag1 
c  (=difft1) as real*8 numbers in units of days, to a precision
c  of 0.01 sec. (0.01 second  = 1/8640000 of a day = 1/hsec_per_day.)
c
	diff21 = (numday2 - numday1) + (hunsec2 - hunsec1)/hsec_per_day
	difft1 = (num_tgt - numday1) + 
     .				    (hunsec_tgt - hunsec1)/hsec_per_day
	coefs(2) = difft1/diff21
	coefs(1) = 1. - coefs(2)
c
c  One coefficient will be negative if target lies outside [ttag1, ttag2].
c
	if (coefs(1)*coefs(2) .lt. 0.) then
	   coefs(2) = 0.
	   coefs(1) = 0.
	endif
	return
	end
