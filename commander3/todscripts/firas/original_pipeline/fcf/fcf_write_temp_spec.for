	integer * 4 function fcf_write_temp_spec (flag, lun, num, spec_recs)

c-------------------------------------------------------------------------------
c
c	Function FCF_WRITE_TEMP_SPEC
c
c	This function writes spectrum records to the temporary RMS output file.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  13 March 1992
c
c-------------------------------------------------------------------------------
c
c	Inputs:
c		flag		integer * 4		spectrum type flag
c							1 = voltage
c							2 = differential
c							3 = calibrated
c		lun		integer * 4		logical unit number
c		num		integer * 4		number of spectra in
c							current ensemble
c		spec_recs(num)				spectrum records
c
c	Outputs:
c		none
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
	integer * 4	io_stat		!  write return status
	integer * 4	j		!  a counter
	integer * 4	lun		!  RMS file logical unit number
	integer * 4	num		!  number of spectra to be written
	integer * 4	status		!  return status

	dictionary 'fcf_sky'
	record /fcf_sky/ spec_recs(num)

	external	fcf_normal
	external	fcf_tempwrite


c
c  Write calibrated spectra to the temporary output file.
c
	status = %loc(fcf_normal)
	j = 1

	do while ((status .eq. %loc(fcf_normal))  .and.  (j .le. num))

	   write (unit=lun, iostat=io_stat) spec_recs(j)

	   if (io_stat .ne. 0) then
	      status = %loc(fcf_tempwrite)
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
	      call lib$signal (fcf_tempwrite, %val(2), temp_file(1:temp_len),
     .					      %val(io_stat))
	   end if

	   j = j + 1

	end do


	fcf_write_temp_spec = status

	return
	end
