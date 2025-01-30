	integer * 4 function  ffi_update_report ()

c-------------------------------------------------------------------------------
c
c	Function FFI_UPDATE_REPORT
c
c	This function writes the model solution label and model solution file
c	name, the input and output file names, and the number of spectra
c	processed to the processing report.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  9 March 1993
c		  SER 10763
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
c		ffi_invoc.txt
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
	include '(ffi_invoc)'

	integer * 4	io_stat		!  I/O return status

	external	ffi_normal
	external	ffi_repwrite

c
c  Write out the files read and written.
c
	write (fut_report_lun,10,iostat=io_stat)
	write (fut_report_lun,21,iostat=io_stat) fcc_infile1(1:28),
     .						 fcc_infile1(29:fcc_inlen1)
	if (fcc_hybrid .eq. fac_present) then
	   write (fut_report_lun,22,iostat=io_stat) fcc_infile2(1:28),
     .						    fcc_infile2(29:fcc_inlen2)
	endif
	write (fut_report_lun,30,iostat=io_stat) fcc_outfile(1:fcc_outlen)
  10	format (//, x, 'FFI Processing Summary:', /)
  21	format ( 4x, 'Input Coadded IFG        ', a, /,
     .		 4x, '           File 1:                     ', a) 
  22	format ( 4x, 'Input Coadded IFG        ', a, /,
     .		 4x, '           File 2:                     ', a) 
  30	format ( 4x, 'Output Voltage', /,
     .		 4x, '    Spectrum File:       ', a) 

c
c  Write the number of spectra processed:
c
	if (fcc_hybrid .eq. fac_present) then
	   write (fut_report_lun,40,iostat=io_stat) fcc_nspec, fcc_nspec1,
     .						    fcc_nspec2
	else
	   write (fut_report_lun,50,iostat=io_stat) fcc_nspec
	endif
  40	format (/, 4x, 'Number of calibration spectra processed:  ', I5,
     .		/, 4x, '    Number from input 1:                  ', I5,
     .		/, 4x, '    Number from input 2:                  ', I5, //)
  50	format (/, 4x, 'Number of calibration spectra processed:  ', I5, //)

c
c  Check the status of the writes:
c
	if (io_stat .ne. 0) then
	   ffi_update_report = %loc(ffi_repwrite)
	   call lib$signal (ffi_repwrite, %val(2),
     .			    fcc_report_file(1:fcc_replen), %val(io_stat))
	else
	   ffi_update_report = %loc(ffi_normal)
	endif


	return
	end
