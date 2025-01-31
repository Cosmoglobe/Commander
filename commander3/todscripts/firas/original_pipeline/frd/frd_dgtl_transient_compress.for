	integer * 4 function  frd_dgtl_transient_compress (ngroup, offset,
     .							   response, trf)

c-------------------------------------------------------------------------------
c
c	Function FRD_DGTL_TRANSIENT_COMPRESS
c
c	This function compresses the digital filter transient response function
c	by the adds per group and the microprocessor offset for a given scan 
c	mode.  It then normalizes the response function to a magnitude of unity.
c	A 128-point response function is returned, since this is all that is
c	needed for fitting transient responses in the ifgs.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  11 September 1991, SER 7985
c
c-------------------------------------------------------------------------------
c
c	Input:
c		ngroup		integer * 4		adds per group
c		offset		integer * 4		microprocesor offset
c		response(512)	real	* 4		uncompressed digital
c							filter transient 
c							response function
c
c	Output:
c		trf(128)	real	* 4		normalized, compressed
c							digital filter
c							transient response 
c							function
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	integer * 4	j		!  a counter
	integer * 4	k		!  a counter
	integer * 4	lstart		!  first index for compression
	integer * 4	lstop		!  last index for compression
	integer * 4	ngroup		!  adds per group
	integer * 4	offset		!  microprocessor offset
	integer * 4	ulim		!  number of compressed data points

	real	* 4	adds		!  adds per group for compression
	real	* 4	buffer		!  compression data buffer
	real	* 4	cresp(512)	!  normal compressed transient function
	real	* 4	norm		!  normalization factor for trf
	real	* 4	response(512)	!  uncompressed transient function
	real	* 4	sum_sq		!  sum squares of transient function
	real	* 4	trf(128)	!  normal compressed transient function

	external	frd_normal


c
c  Initialize the normalized, compressed transient response function and
c	compute the constants for the compression.
c
	call lib$movc5(0,,0,2048,cresp)
	call lib$movc5(0,,0,512,trf)
	adds   = floatj(ngroup)
	sum_sq = 0.0
	ulim   = jnint(floatj(513-offset)/floatj(ngroup))

c
c  Compress the response function.
c
	do j = 1,ulim

	   lstart = (j-1)*ngroup + offset
	   lstop  = lstart + ngroup - 1
	   if (lstop .gt. 512) lstop = 512

	   buffer = 0.0
	   do k = lstart,lstop
	      buffer = buffer + response(k)
	   enddo
	   cresp(j) = buffer/adds
	   sum_sq = sum_sq + cresp(j)**2

	enddo

c
c  Normalize the response function
c
	norm = sqrt(sum_sq)
	do j = 1,128
	   trf(j) = cresp(j)/norm
	enddo


	frd_dgtl_transient_compress = %loc(frd_normal)

	return
	end
