	Program frd_dgtl_transient

c-------------------------------------------------------------------------------
c
c	Program FRD_DGTL_TRANSIENT
c
c	This program generates a direct-access file containing transient 
c	response functions for the Firas digital filters.  The program generates
c	two response functions, one for the high frequency channels and one for
c	the low frequency channels.  These functions are then compressed 
c	by the adds per group and microprocessor offset for each scan mode into 
c	128-point, real*4 arrays.  The eight resulting response functions are
c	written into the binary file FEX_DTRF.DAT.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  11 September 1991, SER 7985
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		frd_dgtl_transient_compress
c		frd_dgtl_transient_fcn
c		fut_free_lun
c		fut_get_lun
c
c	Include files:
c		$ssdef
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '($ssdef)'

	character * 12  trf_file	!  transient response function file name

	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	k		!  a counter
	integer * 4	ngroup(4,2)	!  adds per group
	data  ngroup	/ 3, 2, 3, 2, 3, 2, 12, 8 /
 	integer * 4	nrec		!  transient response function rec num
	integer * 4	offset(4,2)	!  microprocessor offset
	data  offset	/ 4, 3, 10, 10, 4, 3, 10, 9 /
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status
	integer * 4	t_lun		!  transient response function file lun

	real	* 4	response(512)	!  transient response function
	real	* 4	trf(128)	!  compressed transient response fcn

	integer * 4	frd_dgtl_transient_compress
	integer * 4	frd_dgtl_transient_fcn
	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	external	frd_normal
	external	frd_rmsclose
	external	frd_rmsopen
	external	frd_rmswrite
	external	fut_normal

c
c  Open the transient response function file.
c
	status = fut_get_lun(t_lun)

	if (status .eq. %loc(fut_normal)) then
	   
	   trf_file = 'fex_dtrf.dat'
	   open (unit=t_lun, file=trf_file, status='new', recl=128,
     .		 form = 'unformatted', access='direct', iostat = io_stat)


	   if (io_stat .eq. 0) then
	      status = %loc(frd_normal)
	      nrec = 0
c
c  Loop over high and low frequencies.
c
	      do j = 1,2		!  j=1 ==> high;  j=2 ==> low

c
c  Generate the digital transient response functions.
c
	         status = frd_dgtl_transient_fcn (j, response)

c
c  Compress the transient response function by the adds per group and the
c	microprocessor offset for each scan mode.  Normalize the compressed
c	function to a value of unity.
c
	         do k = 1,4
	            if (io_stat .eq. 0) then
	               status = frd_dgtl_transient_compress (ngroup(k,j),
     .							     offset(k,j),
     .							     response, trf)
c
c	Write the normalized, compressed functions to the output file.
c
	               nrec = nrec + 1
	               write (t_lun'nrec, iostat = io_stat) trf
	            endif
	         enddo
	      enddo

c
c  Check the status of the writes to the transient response function file and
c	close the file.
c
	      if (io_stat .ne. 0) then
	         status = %loc(frd_rmswrite)
	         call lib$signal (frd_rmswrite,%val(2),trf_file,%val(io_stat))
	      endif

	      close (t_lun, iostat=io_stat)
	      if (io_stat .ne. 0) then
	         status = %loc(frd_rmsclose)
	         call lib$signal (frd_rmsclose,%val(2),trf_file,%val(io_stat))
	      endif

	   else			!  File open status
	      status = %loc(frd_rmsopen)
	      call lib$signal (frd_rmsopen,%val(2),trf_file,%val(io_stat))
	   endif

	   rstatus = fut_free_lun(t_lun)

	else			!  FUT
	   call lib$signal (%val(status))
	endif


c
c  Terminate the program
c
	if (status .eq. %loc(frd_normal)) then
	   call lib$signal(frd_normal)
	else
	   call lib$signal(%val(ss$_abort))
	endif


	end
