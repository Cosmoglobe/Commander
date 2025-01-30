	subroutine ftb_plotspec (time, avg_time, timebuff,
     1				 freqbuff, pflag, plotmode,
     2				 plt_com, plt_com_file)

c----------------------------------------------------------------------------
c
c Changes:
c
c	SPR 1351, Convert to new IMSL. R. Kummerer, April 14, 1988.
c
c	Added the include file, FTB_Labels, and a call to FTB_Vecplt which
c	interfaces with PLT.           J. Caldwell, STX, September, 1988.
c
c	Added the capability of obtaining automatic hardcopy plots.
c	                               Shirley M Read, STX, September 1988.
c
c       Shirley M. Read, STX, October, 1988 : Corrected the inclusion of
c	     the FTB_Labels text file in the subroutine. (SPR 2558)
c
c       R. Kummerer, STX, May, 1989 : Use FUT_VECPLT instead of FTB_VECPLT.
c
c	R. Isaacman, GSC, Sept 1989: Apodize time series before FFT
c            SPR 4586 - Incorrect spectra in MTM_Jitter. 09/21/89
c	     Put under configuration. SMR
c
c	S. Alexander, STX, March, 1990:  SER 5726 -- Add capability to pass
c	     in a PLT command file.
c
c	F. Shuman, STX, 1990 Mar 09:  SPRs 5396, 5410 -- recent Jitter runs
c	     disagree with old ones.
c
c----------------------------------------------------------------------------

	implicit none

C	Functions

	integer*4 FUT_Vecplt

C	Local Declarations

	logical*1 pflag
	real*4 time(1), timebuff(1), freqbuff(1), freqs(512)
	real*4 avg_time, df, apod
	integer*4 status
	integer*4 j, k, npts
	integer*4 plotmode	!CRT screen or auto hardcopy
	integer*4 plt_com
	character*64 plt_com_file
	complex*8 zspec(1024)

	include '(fut_vecplt_labels)'

	xlabl = 'Frequency (Hz)'
	ylabl = 'Spectral Amplitude (uSec/Hz)'
	df = 1.e06/avg_time/1024.

	do j=1,1024
	   apod = 2./3.*( 1. + cos((j-512.5)*3.1415927/511.5) )
	   zspec(j) = cmplx(time(j)-avg_time,0.) * apod
	enddo

	call fftcf (1024, zspec, zspec)

	do k=1,512
	   time(k) = sqrt (real (zspec(k) * conjg(zspec(k))))/512.
	   timebuff(k) = time(k)
	   freqs(k) = (k - 1) * df
	   freqbuff(k) = freqs(k)
	enddo

	npts = 512

	if(pflag) then
	  status = fut_vecplt (freqbuff,timebuff,npts,plotmode,
	1	   0,0,0,0,0,0,0,0,plt_com,plt_com_file)
	endif

	return
	end
