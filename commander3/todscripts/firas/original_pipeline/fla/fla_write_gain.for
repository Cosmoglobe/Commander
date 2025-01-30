	integer * 4 function  fla_write_gain ()

c-------------------------------------------------------------------------------
c
c	Function FLA_WRITE_GAIN
c
c	This function computes the FIRAS gain function for a given channel and
c	scan mode.  It then writes the gain function to the ASCII file
c	CSDR$FIRAS_OUT:FEX_GAIN_CCSS.VVV_XXXXXXXXXX, where CCSS is channel and
c	scan mode, and VVV_XXXXXXXXXX is the calibration model solution file
c	name extension.
c
c	Author:  Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 1 July 1993
c
c------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		none
c
c	Subroutines called:
c		fut_free_lun
c		fut_get_lun
c		str$trim
c
c	Include files:
c		fla_config_gain.txt
c		fla_invoc_gain.txt
c		fla_model.txt
c
c-------------------------------------------------------------------------------
c
c	Hard-Coded Constants:
c
c		TRLIM		1.0D-10		lower limit for a non-zero
c						transfer function
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Changed LOLIM and UPLIM indices to accomodate low frequency short fast
c	calibration model solutions.
c	Gene Eplee, GSC, 19 November 1993
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fla_config_gain)'
	include '(fla_invoc_gain)'
	include '(fla_model)'

	complex * 16	gain(257)	!  gain function

	integer *  4	g_lun		!  gain function lun
	integer *  4	io_stat		!  I/O return status
	integer * 4	is		!  scan mode pointer
	integer *  4	j		!  a counter
	integer *  4	rstatus		!  return status
	integer *  4	status		!  return status

	real	*  8	trlim		!  lower limit for transfer function
	parameter	(trlim = 1.0D-10)

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	external	fla_normal
	external	fla_rmsclose
	external	fla_rmsopen
	external	fla_rmswrite
	external	fut_normal

C
C  Generate the gain function.
C

	status = %loc(fla_normal)
	call lib$movc5 (0,,0,4112,gain)

c
c  Loop over frequencies, checking for zeros in the transfer function.
c
	do j = 2,257

	   if (cdabs(fex_model.transfer(j)) .gt. trlim) then
	      gain(j) = - dcmplx(1.0D0,afreq(j)*tau)
     .			 / S0 / ztrans(j) / fex_model.transfer(j) / spec_norm
	   endif
	enddo


C
C  Write the gain function to an ASCII file.
C

c
c  Get the gain function file name.
c
	fcc_gain_file = 'CSDR$FIRAS_OUT:FEX_GAIN_' // fcc_scan_mode // '.' //
     .			 fcc_model_ext
	call str$trim (fcc_gain_file, fcc_gain_file, fcc_glen)

c
c  Open the gain function file.
c
	rstatus = fut_get_lun(g_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=g_lun, file=fcc_gain_file, status='new', form='formatted',
     .	      access='sequential', iostat=io_stat)


	if (io_stat .eq. 0) then
c
c  Write the gain function.
c
	   write (g_lun,10,iostat=io_stat) fex_model.mod_head.label
  10	   format (x,a)
	   is = 3*(fcc_chan-1) + (fcc_speed+1) + fcc_length
	   do j=lolim(is),uplim(is)
	      write (g_lun,20,iostat=io_stat) freq(j), dreal(gain(j)),
     .						       dimag(gain(j))
	   enddo
  20	   format (2x,f10.7,2x,e20.13,2x,e20.13)

	   if (io_stat .ne. 0) then
	      status = %loc(fla_rmswrite)
	      call lib$signal (fla_rmswrite, %val(2), fcc_gain_file(1:fcc_glen),
     .					     %val(io_stat))
	   endif

c
c  Close the gain function file.
c
	   close (g_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fla_rmsclose)
	      call lib$signal (fla_rmsclose, %val(2), fcc_gain_file(1:fcc_glen),
     .					     %val(io_stat))
	   endif
	   rstatus = fut_free_lun(g_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal(rstatus)
	   endif

	else
	   status = %loc(fla_rmsopen)
	   call lib$signal (fla_rmsopen, %val(2), fcc_gain_file(1:fcc_glen),
     .					 %val(io_stat))
	endif		!  status from open


	fla_write_gain = status

	return
	end
