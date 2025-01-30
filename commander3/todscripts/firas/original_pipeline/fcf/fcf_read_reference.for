	integer * 4 function  fcf_read_reference (fakeit, upmode)

c-------------------------------------------------------------------------------
c
c	Function FCF_READ_REFERENCE
c
c	This function determines the IFG peak position for the appropriate 
c	channel, scan mode, and microprocessor mode.  It reads and FFTs the
c	appropriate apodization function from the reference dataset.  It then
c	reads the appropriate electronics transfer function from the
c	reference dataset.
c
c	Author:  Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 9 October 1992
c
c-------------------------------------------------------------------------------
c
c	Input:
c		fakeit		integer * 4		fakeit flag
c		upmode		integer * 4		microprocessor mode
c
c	Output:
c		none
c
c	Subroutines called:
c		dfftrf
c		fut_apod_recnum
c		fut_default_peak
c		fut_get_recnum
c		lib$movc5
c		lib$signal
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

	integer * 4	arecno		!  apodization function record number
	integer * 4	erecno		!  etf record number
	integer * 4	fakeit		!  fakeit flag
	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	k		!  a counter
	integer * 4	linearized/1/	!  linearization flag
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status
	integer * 4	upmode 		!  microprocessor mode

	real	* 8	apod(512)	!  apodization function
	real	* 8	norm		!  normalization factor for
	parameter	(norm = 512.0D0**2)  !  FFT'ed apodization function
	real	* 8	norm_short	!  normalization factor for short
	parameter	(norm_short = 128.0D0**2)  ! FFT'ed apodization function
	real	* 8	rapod(512)	!  FFT'ed apodization function
	real	* 8	zapod(257)	!  unpacked FFT'ed apodization
					!    function

	integer * 4	fut_apod_recnum
	integer * 4	fut_default_peak

	external	fcf_normal
	external	fcf_readapod
	external	fcf_readetf

	status = %loc(fcf_normal)

C
C  Determine the IFG peak position.
C
	if ((fcc_xlsf .eq. fac_present)  .or. (fcc_xllf .eq. fac_present)) then
	   rstatus = fut_default_peak (1, 1, fcc_chan, 8, upmode, linearized,
     .				       ictr)
	else
	   rstatus = fut_default_peak (fcc_speed, fcc_length, fcc_chan,
     .				       fcc_ngroup, upmode, linearized, ictr)
	endif


C
C  Read and FFT the appropriate apodization function.
C

c
c  Get the record number.
c
	rstatus = fut_apod_recnum (fcc_speed, fcc_ngroup, fakeit, fcc_length,
     .				   fcc_chan, upmode, linearized, arecno)

c
c  Read the apodization function.
c
	call lib$movc5 (0,,0,4096,apod_fcn)
	read (dir_lun(1)'arecno, iostat=io_stat) apod
	if ((fcc_xlsf .eq. fac_present)  .or.  (fcc_xllf .eq. fac_present)) then
	   do j=1,127
	      apod_fcn(j) = (apod(4*j+1) + apod(4*j+2)) / 2.0D0
	   enddo
	else
	   do j=1,512
	      apod_fcn(j) = apod(j)
	   enddo
	endif

	if (io_stat .ne. 0) then
	   status = %loc(fcf_readapod)
	   call lib$signal (fcf_readapod, %val(2), %val(arecno), %val(io_stat))
	endif

	if (status .eq. %loc(fcf_normal)) then
c
c  Produce the zapodization function by FFT'ing the apodization function.
c
	   call lib$movc5 (0,,0,4096,rapod)
	   call lib$movc5 (0,,0,2056,zapod)
	   call lib$movc5 (0,,0,8192,zapod_fcn)

	   if ((fcc_xlsf .eq. fac_present)  .or.
     .	       (fcc_xllf .eq. fac_present)) then
c
c  Do a 128-point fft for the low frequency fast ifgs.
c
	      call dfftrf (128, apod_fcn, rapod)
	      zapod(1) = rapod(1)**2/norm_short
	      do j = 2,126,2
	         k = j/2 + 1
	         zapod(k) = (rapod(j)**2 + rapod(j+1)**2)/norm_short
	      enddo
	      zapod(65) = rapod(128)**2/norm_short

	      zapod_fcn(1)   = zapod(1)
	      zapod_fcn(129) = zapod(1)
	      do j = 2,64
	         zapod_fcn(j)      = zapod(j)
	         zapod_fcn(130-j)  = zapod(j)
	         zapod_fcn(128+j)  = zapod(j)
	         zapod_fcn(258-j) = zapod(j)
	      enddo
	      zapod_fcn(65) = zapod(65)
	      zapod_fcn(193) = zapod(65)

	   else
c
c  Do a 512-point fft for the rest of the data.
c
	      call dfftrf (512, apod_fcn, rapod)
	      zapod(1) = rapod(1)**2/norm
	      do j = 2,510,2
	         k = j/2 + 1
	         zapod(k) = (rapod(j)**2 + rapod(j+1)**2)/norm
	      enddo
	      zapod(257) = rapod(512)**2/norm

	      zapod_fcn(1)   = zapod(1)
	      zapod_fcn(513) = zapod(1)
	      do j = 2,256
	         zapod_fcn(j)      = zapod(j)
	         zapod_fcn(514-j)  = zapod(j)
	         zapod_fcn(512+j)  = zapod(j)
	         zapod_fcn(1026-j) = zapod(j)
	      enddo
	      zapod_fcn(257) = zapod(257)
	      zapod_fcn(769) = zapod(257)
	   endif


C
C  Read the appropriate electronics transfer function.
C

c
c  Get the record number.
c
	   call fut_get_recnum (fakeit, fcc_speed, fcc_chan, upmode, fcc_ngroup,
     .			        erecno)

c
c  Read the ETF.
c
	   read (dir_lun(2)'erecno, iostat=io_stat) ztrans

	   if (io_stat .ne. 0) then
	      status = %loc(fcf_readetf)
	      call lib$signal (fcf_readetf, %val(1), %val(erecno),
     .					    %val(io_stat))
	   endif

	endif	!  status from apodization function read


	fcf_read_reference = status

	return
	end
