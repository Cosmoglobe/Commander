	Program FRD_APOD

c-------------------------------------------------------------------------------
c
c	FRD_APOD
c
c	This program generates a direct-access file containing apodization
c	functions for coadded IFGs.  The apodization functions are 512-point,
c	real*8 arrays that have not been rotated to put the peak in the first
c	position.  They are written into the binary file FEX_APOD.DAT
c
c	This program is based on the program FRD_RPAD_APOD written by Rich
c	Isaacman.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  28 January 1992
c		  SER 8836
c	
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		frd_apod_fcn
c		fut_default_peak
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
c	Shift low frequency short fast peak position by -3 to allow data
c	recovery.  Gene Eplee, GSC, 6 August 1993, SER 11394.
c
c-------------------------------------------------------------------------------

	implicit none

	include '($ssdef)'

	character * 12  apod_file	!  apodization function file name

	integer * 4	a_lun		!  apodization function file lun
	integer * 4	chan		!  channel number
	integer * 4	ictr		!  ifg peak position
	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	linearized	!  ifg linearization flag
	integer * 4	mtm_length	!  scan length
	integer * 4	mtm_speed	!  scan speed
	integer * 4	ngroup		!  adds per group
	integer * 4	nrec		!  apodization function record number
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status
	integer * 4	upmode		!  microprocessor mode

	real	* 8	apod_fcn(512)	!  apodization function

	integer * 4	frd_apod_fcn
	integer * 4	fut_default_peak
	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	external	frd_normal
	external	frd_rmsclose
	external	frd_rmsopen
	external	frd_rmswrite
	external	fut_normal

c
c  Open the apodization function file.
c
	status = fut_get_lun(a_lun)

	if (status .eq. %loc(fut_normal)) then
	   
	   apod_file = 'fex_apod.dat'
	   open (unit=a_lun, file=apod_file, status='new', recl=1024,
     .		 form = 'unformatted', access='direct', iostat = io_stat)


	   if (io_stat .eq. 0) then
	      status = %loc(frd_normal)
	      nrec = 0
c
c  Cycle through the instrument states.
c
	      do ngroup = 1,12			!  adds per group
	       do mtm_speed = 0,1		!  scan speed (slow/fast)
	        do mtm_length = 0,1		!  scan length (short/long)
	         do chan = 2,3			!  channel (low/high)
	          do upmode = 1,2		!  digital filters off/on
	           do linearized = 0,1		!  unlinearized/linearized

	              if (io_stat .eq. 0) then
c
c  For each instrument state, determine the IFG peak position, then generate
c  the apodization function.
c
	                 rstatus = fut_default_peak (mtm_speed, mtm_length,
     .					 chan, ngroup, upmode, linearized, ictr)
	                 if ((mtm_speed .eq. 1)  .and.  (mtm_length .eq. 0)
     .			     .and.  (chan .eq. 2)) ictr = ictr - 3
	                 status = frd_apod_fcn (ictr, rstatus, apod_fcn)

c
c  Write the apodization function to the output file
c
	                 nrec = nrec + 1
	                 write (a_lun'nrec, iostat = io_stat) apod_fcn

	              endif

	           enddo
	          enddo
	         enddo
	        enddo
	       enddo
	      enddo

c
c  Generate an apodization function full of 1's for fakeit data.
c
	      if (io_stat .eq. 0) then
	         do j = 1,512
	            apod_fcn(j) = 1.0D0
	         enddo
	         write (a_lun'nrec+1, iostat = io_stat) apod_fcn
	      endif


c
c  Check the status of the writes to the apodization function file and close the
c  file.
c
	      if (io_stat .ne. 0) then
	         status = %loc(frd_rmswrite)
	         call lib$signal (frd_rmswrite,%val(2),apod_file,%val(io_stat))
	      endif

	      close (a_lun, iostat=io_stat)
	      if (io_stat .ne. 0) then
	         status = %loc(frd_rmsclose)
	         call lib$signal (frd_rmsclose,%val(2),apod_file,%val(io_stat))
	      endif

	   else			!  File open status
	      status = %loc(frd_rmsopen)
	      call lib$signal (frd_rmsopen,%val(2),apod_file,%val(io_stat))
	   endif

	   rstatus = fut_free_lun(a_lun)

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
