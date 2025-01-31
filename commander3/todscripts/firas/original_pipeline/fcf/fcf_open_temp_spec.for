	integer * 4 function fcf_open_temp_spec (flag, lun)

c-------------------------------------------------------------------------------
c
c	Function FCF_OPEN_TEMP_SPEC
c
c		This function opens the temporary RMS spectra output file for
c	sequential write access.
c
c	Author:		Gene Eplee
c			General Sciences Corp.
c			513-7768
c			5 June 1992
c
c-------------------------------------------------------------------------------
c
c	Inputs:
c		flag		integer * 4		spectrum type flag
c							1 = voltage
c							2 = differential
c							3 = calibrated
c
c	Outputs:
c		lun		integer * 4		temporary file lun
c
c	Subroutines called:
c		lib$signal
c
c	Include files:
c		fcf_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Put in optional differential spectrum output.
c	Gene Eplee, GSC, 11 July 1994
c	SPR 11826
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fcf_invoc)'

	character * 52	temp_file	!  temporary file name

	integer * 2	temp_len	!  file name buffer length

	integer * 4	flag		!  spectrum type flag
	integer * 4	io_stat		!  open return status
	integer * 4	lun		!  RMS file logical unit number
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	fut_get_lun

	external	fcf_normal
	external	fcf_tempopen
	external	fut_normal


	status = %loc(fcf_normal)

c
c  Determine the file name.
c
	if (flag .eq. 1) then
	   temp_file = fcc_temp_outfile_vs
	   temp_len  = fcc_toutlenvs
	elseif (flag .eq. 2) then
	   temp_file = fcc_temp_outfile_ds
	   temp_len  = fcc_toutlends
	elseif (flag .eq. 3) then
	   temp_file = fcc_temp_outfile_cs
	   temp_len  = fcc_toutlencs
	endif

c
c  Open the temporary output file for direct write access.
c
	rstatus = fut_get_lun(lun)

	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=lun, file=temp_file, access='sequential', status='new',
     .	      form='unformatted', recordtype='fixed', recl=fcc_record_size,
     .	      iostat=io_stat)
	if (io_stat .ne. 0) then
           status = %loc(fcf_tempopen)
	   call lib$signal (fcf_tempopen, %val(2), temp_file(1:temp_len),
     .					  %val(io_stat))
	end if


	fcf_open_temp_spec = status

	return
	end
