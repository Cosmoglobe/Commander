	subroutine ftb_mtm_sumplot (histo, spec, avg_time, rms, 
     1				    nrec, ybuff, xbuff, plotmode,
     2				    plt_com, plt_com_file, spec2, avg_time2)
c
c Module Name: FTB_Mtm_Sumplot
c
c Purpose: To
c----------------------------------------------------------------------------
c
c Changes:
c	Jody Caldwell, STX, September, 1988 : Added the FTB_Labels include 
c            file and a call to FTB_Vecplt which interfaces with PLT.(SPR 2161)
c
c       Shirley M. Read, STX, September, 1988 : Added the capability of
c	     automatic hardcopy. (SPR 2430)
c
c       Shirley M. Read, STX, October, 1988 : Corrected the inclusion of 
c	     the FTB_Labels text file in the subroutine. (SPR 2558)
c
c       R. Kummerer, STX, May, 1989 : Use FUT_VECPLT instead of FTB_VECPLT.
c
c	S. Alexander, STX, March, 1990 : Add capability to pass in a
c		PLT command file.
c	L. Rosen, HSTX, June 1992 : Make additional plot of spec2, the fft of
c		the rms average, instead of the average of the fft. SER 2865,
c               SPR 9765.
c----------------------------------------------------------------------------

        implicit none

c	External Messages

	external FUT_Normal
	external FTB_Normal
	external FTB_Aberr

c	Functions

	integer*4 FUT_Vecplt

c	Local Declarations

        real*4 histo(1), spec(1), ybuff(1), xbuff(1)
	real*4 freqs(512), spec2(1)
        real*4 avg_time, df, rms, avg_time2
	integer*4 plotmode	 !Screen display or hardcopy
	integer*4 plt_com
	character*64 plt_com_file
        integer*4 j, k, npts, nrec
	integer*4 status
	logical*1 proceed
	include '(fut_vecplt_labels)'

	proceed = .true.
c
c  First plot the grand average spectrum with correct frequency scale
c
	xlabl = 'Frequency (Hz)'
	ylabl = 'Spectral Amplitude (uSec/Hz)'
	title = 'Average Jitter Spectrum '//time_range
	df = 1.e06/avg_time/1024.
	do k=1,512
	   ybuff(k) = spec(k)
	   freqs(k) = (k - 1) * df
	   xbuff(k) = freqs(k)
	enddo
	npts = 512
	status = fut_vecplt (xbuff,ybuff,npts,plotmode,0,0,0,0,0,0,0,0,
     1			     plt_com, plt_com_file)
	if ( status .ne. %loc(FUT_Normal)) then
	  proceed = .false.
	endif
c
c  Second plot the (non-magnitude averaged) average spectra.
c
        if ( proceed ) then
	   xlabl = 'Frequency (Hz)'
	   ylabl = 'Spectral Amplitude (uSec/Hz)'
	   title = '(non-magnitude) Average Jitter Spectrum '//time_range
	   df = 1.e06/avg_time2/1024.
	   do k=1,512
	      ybuff(k) = spec2(k)
	      freqs(k) = (k - 1) * df
	      xbuff(k) = freqs(k)
	   enddo
	   npts = 512
	   status = fut_vecplt (xbuff,ybuff,npts,plotmode,0,0,0,0,0,0,0,0,
     1			     plt_com, plt_com_file)
	   if ( status .ne. %loc(FUT_Normal)) then
	      proceed = .false.
	   endif
	endif
c
c  Now plot the pulse interval histogram
c
	if ( proceed ) then
	  xlabl = 'Sample Pulse Interval (\guSec)'
	  ylabl = 'Number of Samples'
	  title = 'Pulse Interval Histogram '//time_range
	  write(oxlabl,'(a,f6.1,a,f6.1,a)')
     &	     'All Data Average Interval: ',avg_time,' +/- ',rms,
     &	     ' (rms) \guSec'
	  do k=1,256
	     ybuff(k) = histo(k)
	     freqs(k) = (k - 1)*10. + 45.	!Actually times in uSec
	     xbuff(k) = freqs(k)
	  enddo
c
c  Find the non-zero range
c
	  k = 1
	  do while( (histo(k).eq. 0) .and. (k.lt.256))
	     k = k + 1
	  enddo
	  k = max ( 1, k-5)  !  Move back 4 points from the first non-zero
	  npts = 256
	  do while( (histo(npts).eq. 0) .and. (npts.gt.k))
	     npts = npts - 1
	  end do
	  npts = min (256, npts+5)  ! Move forward 4 points from last non-zero
	  npts = npts - k + 1
	  status = fut_vecplt ( xbuff(k), ybuff(k),
     &	     npts,plotmode,0,0,0,0,0,0,0,0, plt_com, plt_com_file)
c
c  Finally, display the values of the grand averages
c
	  call lib$erase_page(1,1)
	  type 100, nrec, avg_time, 1.e06/avg_time, rms
100	  format (///' Number of records: ',i3,//,
     1		' Average sample interval in the time range: ',
     2			f6.1,' microsec.'//' Corresponding',
     3			' sampling frequency: ',f6.2,' Hz'//
     4			' RMS sampling jitter: ',f6.1,' microsec.'//)
	endif
	return
	end
