	Integer * 4 Function  ffl_update_report (ffli)

c-------------------------------------------------------------------------------
c
c	Function FFL_UPDATE_REPORT
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
c	    ffli record		invoc structure (defined in ffl_invoc.txt)
c
c	Output:
c	    none
c
c	Subroutines called:
c	    LIB$Signal
c
c	Include files:
c	    ffl_invoc.txt
c	    fut_error.txt
c	    fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Name of facility changed from FFI to FFL.  Fred Shuman, HSTX
c	1995 May 19.
c
c	In order to remove some of the obscurity arising from Include files that
c	   harbor hidden Common blocks, converted these Commons to Structures.
c	   This necessitates adding their Record names to the calling lists of
c	   functions that use their variables.
c	Fred Shuman, HSTX, 1995 June 14.
c
c-------------------------------------------------------------------------------

	Implicit None

	Include '(fut_error)'
	Include '(fut_params)'
	Include '(ffl_invoc)'

	Integer * 4	io_stat		!  I/O return status

	External	ffl_normal
	External	ffl_repwrite

c
c  Write out the files read and written.
c

	Write (fut_report_lun,10,iostat=io_stat)
	Write (fut_report_lun,21,iostat=io_stat) ffli.infile1(1:28),
     &						 ffli.infile1(29:ffli.inlen1)
	If (ffli.hybrid .Eq. fac_present) Then
	   Write (fut_report_lun,22,iostat=io_stat) ffli.infile2(1:28),
     &						    ffli.infile2(29:ffli.inlen2)
	EndIf
	Write (fut_report_lun,30,iostat=io_stat) ffli.outfile(1:ffli.outlen)
  10	Format (//, x, 'FFL Processing Summary:', /)
  21	Format ( 4x, 'Input Coadded IFG        ', a, /,
     &		 4x, '           File 1:                     ', a)
  22	Format ( 4x, 'Input Coadded IFG        ', a, /,
     &		 4x, '           File 2:                     ', a)
  30	Format ( 4x, 'Output Voltage', /,
     &		 4x, '    Spectrum File:       ', a)

c
c  Write the number of spectra processed:
c
	If (ffli.hybrid .Eq. fac_present) Then
	   Write (fut_report_lun,40,iostat=io_stat) ffli.nspec, ffli.nspec1,
     &						    ffli.nspec2
	Else
	   Write (fut_report_lun,50,iostat=io_stat) ffli.nspec
	EndIf
  40	Format (/, 4x, 'Number of calibration spectra processed:  ', I5,
     &		/, 4x, '    Number from input 1:                  ', I5,
     &		/, 4x, '    Number from input 2:                  ', I5, //)
  50	Format (/, 4x, 'Number of calibration spectra processed:  ', I5, //)

c
c  Check the status of the writes:
c
	If (io_stat .Ne. 0) Then
	   ffl_update_report = %Loc(ffl_repwrite)
	   Call LIB$Signal (ffl_repwrite, %Val(2),
     &			    ffli.report_file(1:ffli.replen), %Val(io_stat))
	Else
	   ffl_update_report = %Loc(ffl_normal)
	EndIf


	Return
	End
