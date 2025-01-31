	integer * 4 function  fip_update_sky_report (nspec, npix)

c-------------------------------------------------------------------------------
c
c	Function FIP_UPDATE_SKY_REPORT
c
c	This function writes the input and output file names, the data
c	selection cuts,  and the number of spectra and pixels processed to the
c	FIP_SKY processing report.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  26 May 1993
c
c-------------------------------------------------------------------------------
c
c	Input:
c		nspec		integer * 4		number of spectra
c							processed
c
c		npix		integer * 4		number of pixels
c							containing spectra
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
c		fip_invoc_sky.txt
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
	include '(fip_invoc_sky)'
	include '(fip_frequency)'

	integer * 4	io_stat		!  I/O return status
	integer * 4	npix		!  number of pixels containing spectra
	integer * 4	nspec		!  number of spectra processed

	external	fip_normal
	external	fip_repwrite

c
c  Write out the files read and written.
c
	write (fut_report_lun,10,iostat=io_stat)
	write (fut_report_lun,20,iostat=io_stat) fcc_infile(1:fcc_inlen)
	write (fut_report_lun,30,iostat=io_stat) fcc_outfile(1:fcc_outlen)
	if (fcc_destriped .eq. fac_present) then
	   write (fut_report_lun,40,iostat=io_stat)
	endif
  10	format (//, x, 'FIP_SKY Processing Summary:', /)
  20	format (4x, 'Input Skymap File:       ', a)
  30	format (4x, 'Output Skymap File:      ', a, /)
  40    format (4x, 'Spectra have been destriped.' /)

c
c  Write out the data selection cuts:
c
	if (fcc_galexc .eq. fac_present) then
	   write (fut_report_lun,50,iostat=io_stat) fcc_glat
	else
	   write (fut_report_lun,60,iostat=io_stat)
	endif
	if (fcc_freq .eq. fac_present) then
	   write (fut_report_lun,70,iostat=io_stat) fcc_lofreq, fcc_hifreq
	   write (fut_report_lun,80,iostat=io_stat) fcc_jlo, fcc_jhi
	   write (fut_report_lun,90,iostat=io_stat) fcc_nu0, fcc_dnu, fcc_nfreq
	else
	   write (fut_report_lun,100,iostat=io_stat)
	endif
  50	format (4x, 'Galactic Latitude Cutoff:    ', 11x, F6.2, '     degrees')
  60    format (4x, 'All Galactic Latitudes Included.')
  70	format (4x, 'Frequency Range:             ', 11x, I3, ' - ',I3, '  icm')
  80	format (4x, 'Frequency Indices:           ', 11x, I3, ' - ',I3)
  90	format (4x, 'Initial Optical Frequency:   ', 11x, F7.3, '    GHz', /,
     .		4x, 'Optical Frequency Interval:  ', 11x, F7.3, '    GHz', /,
     .		4x, 'Number of Frequency Points:  ', 11x, I3, /)
 100    format (4x, 'All Frequencies Included.')

c
c  Write the number of spectra processed:
c
	write (fut_report_lun,110,iostat=io_stat) nspec, npix
 110	format (4x, 'Number of sky spectra processed:      ', I5, /,
     .		4x, 'Number of pixels containing spectra:  ', I5, //)

c
c  Check the status of the writes:
c
	if (io_stat .ne. 0) then
	   fip_update_sky_report = %loc(fip_repwrite)
	   call lib$signal (fip_repwrite, %val(2),
     .			    fcc_report_file(1:fcc_replen), %val(io_stat))
	else
	   fip_update_sky_report = %loc(fip_normal)
	endif


	return
	end
