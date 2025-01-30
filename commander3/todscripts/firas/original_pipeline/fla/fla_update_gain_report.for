	integer * 4 function  fla_update_gain_report ()

c-------------------------------------------------------------------------------
c
c	Function FLA_UPDATE_GAIN_REPORT
c
c	This function writes the input and output file names and the data
c	selection cuts to the FLA_GAIN processing report.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  30 June 1993
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		none
c
c	Subroutines called:
c		lib$signal
c		str$trim
c
c	Include files:
c		fla_invoc_gain.txt
c		fut_error.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_error)'
	include '(fut_params)'
	include '(fla_invoc_gain)'

	integer * 4	io_stat		!  I/O return status

	external	fla_normal
	external	fla_repwrite

c
c  Write out the files read and written.
c
	write (fut_report_lun,10,iostat=io_stat)
	write (fut_report_lun,20,iostat=io_stat) fcc_fex_file(1:fcc_fexlen)
	write (fut_report_lun,30,iostat=io_stat) fcc_gain_file(1:fcc_glen)
  10	format (//, x, 'FLA_GAIN Processing Summary:', /)
  20	format (4x, 'FEX Model Solution File:     ', a)
  30	format (4x, 'Gain Function File:          ', a, /)

c
c  Check the status of the writes:
c
	if (io_stat .ne. 0) then
	   fla_update_gain_report = %loc(fla_repwrite)
	   call lib$signal (fla_repwrite, %val(2),
     .			    fcc_report_file(1:fcc_replen), %val(io_stat))
	else
	   fla_update_gain_report = %loc(fla_normal)
	endif


	return
	end
