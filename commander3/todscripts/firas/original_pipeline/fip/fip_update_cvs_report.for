	integer * 4 function  fip_update_cvs_report ()

c-------------------------------------------------------------------------------
c
c	Function FIP_UPDATE_CVS_REPORT
c
c	This function writes the input and output file names and the data
c	selection cuts to the FIP_CVS processing report.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  17 November 1994
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
c		fip_frequency.txt
c		fip_invoc_cvs.txt
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
	include '(fip_invoc_cvs)'
	include '(fip_frequency)'

	integer * 4	io_stat		!  I/O return status

	external	fip_normal
	external	fip_repwrite

c
c  Write out the files read and written.
c
	write (fut_report_lun,10,iostat=io_stat)
	write (fut_report_lun,20,iostat=io_stat) fcc_fms_cfile(1:fcc_fmslen)
	write (fut_report_lun,30,iostat=io_stat) fcc_fip_cfile(1:fcc_fiplen)
  10	format (//, x, 'FIP_CVS Processing Summary:', /)
  20	format (4x, 'FMS C-Vector File:  ', a)
  30	format (4x, 'FIP C-Vector File:  ', a)

c
c  Write out the data selection cuts:
c
	if (fcc_freq .eq. fac_present) then
	   write (fut_report_lun,40,iostat=io_stat) fcc_lofreq, fcc_hifreq
	   write (fut_report_lun,50,iostat=io_stat) fcc_jlo, fcc_jhi
	   write (fut_report_lun,60,iostat=io_stat)fcc_nu0, fcc_dnu, fcc_nfreq
	else
	   write (fut_report_lun,70,iostat=io_stat)
	endif
  40	format (4x, 'Frequency Range:             ', I3, ' - ',I3, '  icm')
  50	format (4x, 'Frequency Indices:           ', I3, ' - ',I3)
  60	format (4x, 'Initial Optical Frequency:   ', F7.3, '    GHz', /,
     .		4x, 'Optical Frequency Interval:  ', F7.3, '    GHz', /,
     .		4x, 'Number of Frequency Points:  ', I3, /)
  70    format (4x, 'All Frequencies Included.')

c
c  Check the status of the writes:
c
	if (io_stat .ne. 0) then
	   fip_update_cvs_report = %loc(fip_repwrite)
	   call lib$signal (fip_repwrite, %val(2),
     .			    fcc_report_file(1:fcc_replen), %val(io_stat))
	else
	   fip_update_cvs_report = %loc(fip_normal)
	endif


	return
	end
