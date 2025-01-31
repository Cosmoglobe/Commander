	Program FRD_APODL

c-------------------------------------------------------------------------------
c
c	FRD_APODL
c
c	This program generates a direct-access file containing apodization
c	functions for coadded IFGs.  The apodization functions are 512-point,
c	real*8 arrays that have not been rotated to put the peak in the first
c	position.  They are written into the binary file FEX_APODL.DAT
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
c	Changes:
c
c	Shift low frequency short fast peak position by -3 to allow data
c	recovery.  Gene Eplee, GSC, 6 August 1993.
c
c          Removed shift in short fast peak position since calls to 
c       FUT_DEFAULT_PEAK are now geared to recover correct peak position 
c       for low short fast smoothed and decimated interferograms.  Added 
c       apodization functons of length 128 for use with low short fast 
c       smoothed or truncated interferograms.
c                        Alice Trenholme, GSC, 7 February 1995
c       
c          Changed peak position for HSF data so that
c              1.  Peak position is shifted DOWN 1/2 bin if nonlinearized 
c                  peak position is used and microprocessor is on
c              2.  Peak position is shifted UP 1/2 bin otherwise
c          Passed flags to FRD_APOD_FCNL if HSF or HLF functions are to be 
c          generated
c                        Alice Trenholme, GSC, 15 February 1995, SER 12244
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		frd_apod_fcnl
c		fut_default_peak
c		fut_free_lun
c		fut_get_lun
c
c	Include files:
c		$ssdef
c
c-------------------------------------------------------------------------------

	implicit none

	include '($ssdef)'

	character * 13  apod_file	!  apodization function file name

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
        integer * 4     short           !  0=write length 512 function
                                        !  1=write length 128 function
        integer * 4     xhlf            !  HLF flag
        integer * 4     xhsf            !  HSF flag

        real    * 4     peak            !  non-integer peak for HSF data
	real	* 8	apod_fcn(512)	!  apodization function

	integer * 4	frd_apod_fcnl
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
	   
	   apod_file = 'fex_apodl.dat'
	   open (unit=a_lun, file=apod_file, status='new', recl=1024,
     .		 form = 'unformatted', access='direct', iostat = io_stat)


	   if (io_stat .eq. 0) then
	      status = %loc(frd_normal)
	      nrec = 0
c
c  Cycle through the instrument states.  This is for the "regular" length
c  512 apodization functions, so set short = 0.
c
              short = 0
	      do ngroup = 1,12			!  adds per group
	       do mtm_speed = 0,1		!  scan speed (slow/fast)
	        do mtm_length = 0,1		!  scan length (short/long)
	         do chan = 2,3			!  channel (low/high)
                  xhsf =  (chan - 2)*  (1 - mtm_length) * (mtm_speed)
                  xhlf =  (chan - 2)*  (mtm_length) * (mtm_speed)
	          do upmode = 1,2		!  digital filters off/on
	           do linearized = 0,1		!  unlinearized/linearized

	              if (io_stat .eq. 0) then
c
c  For each instrument state, determine the IFG peak position, then generate
c  the apodization function.
c
	                 rstatus = fut_default_peak (mtm_speed, mtm_length,
     .					 chan, ngroup, upmode, linearized, ictr)

c
c    Peak shift for High Short Fast data.
c
                         if (xhsf) then
                           if ((linearized .eq. 0) .and. (upmode .eq. 2)) then
                               peak = ictr + .5
                           else
                               peak = ictr - .5
                           endif
                         else
                           peak = ictr
                         endif
c
c
c
	                 status = frd_apod_fcnl (peak, rstatus, short, 
     .                     xhsf, xhlf, apod_fcn)
c
c  Write the apodization function to the output file
c
	                 nrec = nrec + 1
	                 write (a_lun''nrec, iostat = io_stat) apod_fcn

	              endif

	           enddo
	          enddo
	         enddo
	        enddo
	       enddo
	      enddo

c  Cycle through the instrument states.  This is for the short length
c  128 apodization functions, so set short = 1.  Manufacture these only for 
c  special long fast cases (special short fast data will have been decimated
c  and smoothed to match special long fast).
c
              xhsf = 0
              xhlf = 0
              short = 1
	      do ngroup = 1,12			!  adds per group
	       mtm_speed = 1                    !  fast scan speed 
	        mtm_length = 1                  !  long scan length (correct
                                                !   for peak retrieval)
	         chan = 2			!  low channel 
	          do upmode = 1,2		!  digital filters off/on
	           do linearized = 0,1		!  unlinearized/linearized

	              if (io_stat .eq. 0) then
c
c  For each instrument state, determine the IFG peak position, then generate
c  the apodization function.
c
	                 rstatus = fut_default_peak (mtm_speed, mtm_length,
     .					 chan, ngroup, upmode, linearized, ictr)
                         peak = ictr
	                 status = frd_apod_fcnl (peak, rstatus, short, 
     .                     xhsf, xhlf, apod_fcn)
c
c  Write the apodization function to the output file
c
	                 nrec = nrec + 1
	                 write (a_lun'nrec, iostat = io_stat) apod_fcn

	              endif

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
