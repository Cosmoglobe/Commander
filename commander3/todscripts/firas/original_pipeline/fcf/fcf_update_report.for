	integer * 4 function  fcf_update_report ()

c-------------------------------------------------------------------------------
c
c	Function FCF_UPDATE_REPORT
c
c	This function writes the model solution label and model solution file
c	name, the input and output file names, and the number of spectra and
c	pixels processed to the processing report.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  16 June 1992
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
c		fcf_invoc.txt
c		fcf_model.txt
c		fut_error.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Added input file FEX_VAR.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11397
c
c	Put in optional differential spectrum output.
c	Gene Eplee, GSC, 11 July 1994
c	SPR 11826
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_error)'
	include '(fut_params)'
	include '(fcf_invoc)'
	include '(fcf_model)'

	integer * 4	io_stat		!  I/O return status

	external	fcf_normal
	external	fcf_repwrite

c
c  Write out the files read and written.
c
	write (fut_report_lun,10,iostat=io_stat)
	if (fcc_calibrate .eq. fac_present) then
	   write (fut_report_lun,20,iostat=io_stat) model_label(1:mod_lablen)
	   write (fut_report_lun,30,iostat=io_stat) fcc_model_file(1:fcc_modlen)
	   if (fcc_dvec .eq. fac_present) then
	     write (fut_report_lun,35,iostat=io_stat) fcc_var_file(1:fcc_varlen)
	   endif
	endif
	if (fcc_sky .eq. fac_present) then
	   write (fut_report_lun,40,iostat=io_stat) fcc_infile(1:fcc_inlen)
	else
	   write (fut_report_lun,50,iostat=io_stat) fcc_infile(1:27),
     .						    fcc_infile(28:fcc_inlen)
	endif
	if (fcc_write_vs .eq. fac_present) then
	   write (fut_report_lun,60,iostat=io_stat)
     .						  fcc_outfile_vs(1:fcc_outlenvs)
	endif
	if (fcc_write_ds .eq. fac_present) then
	   write (fut_report_lun,70,iostat=io_stat)
     .						  fcc_outfile_ds(1:fcc_outlends)
	endif
	if (fcc_write_cs .eq. fac_present) then
	   write (fut_report_lun,80,iostat=io_stat)
     .						  fcc_outfile_cs(1:fcc_outlencs)
	endif
  10	format (//, x, 'FCF Processing Summary:', /)
  20	format ( 4x, 'Calibration Model', /,
     .		 4x, '   Solution Label:       ', a) 
  30	format ( 4x, 'Calibration Model', /,
     .		 4x, '    Solution File:       ', a) 
  35	format ( 4x, 'Sky Variance File:       ', a)
  40	format ( 4x, 'Input Coadded IFG', /,
     .		 4x, '             File:       ', a) 
  50	format ( 4x, 'Input Coadded IFG', /,
     .		 4x, '             File:       ', a, /,
     .           4x, '                                       ', a) 
  60	format ( 4x, 'Output Voltage', /,
     .		 4x, '    Spectrum File:       ', a) 
  70	format ( 4x, 'Output Differential', /,
     .		 4x, '    Spectrum File:       ', a) 
  80	format ( 4x, 'Output Calibrated', /,
     .		 4x, '    Spectrum File:       ', a) 

c
c  Write the number of spectra processed:
c
	if (fcc_sky .eq. fac_present) then
	   write (fut_report_lun,90,iostat=io_stat) fcc_nspec, fcc_npix
	else
	   write (fut_report_lun,100,iostat=io_stat) fcc_nspec
	endif
  90	format (/, 4x, 'Number of sky spectra processed:      ', I5,
     .		/, 4x, 'Number of pixels containing spectra:  ', I5, //)
 100	format (/, 4x, 'Number of calibration spectra processed:  ', I5, //)

c
c  Check the status of the writes:
c
	if (io_stat .ne. 0) then
	   fcf_update_report = %loc(fcf_repwrite)
	   call lib$signal (fcf_repwrite, %val(2),
     .			    fcc_report_file(1:fcc_replen), %val(io_stat))
	else
	   fcf_update_report = %loc(fcf_normal)
	endif


	return
	end
