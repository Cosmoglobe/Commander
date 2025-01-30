	integer * 4 function  fsl_read_dvector ()

c-------------------------------------------------------------------------------
c
c	Function FSL_READ_DVECTOR
c
c	This function opens, reads, and closes the FIRAS D-Vector reference
c	dataset FEX_VAR_CCSS.VVV_XXXXXXXXXX from the directory CSDR$FIRAS_CAL 
c	(CC is channel, SS is scan mode, VVV is the version of the calibration
c	program that generated the solution, and XXXXXXXXXX is the model
c	solution label).  The file extension is specified by the command line
c	qualifier /MODEL_EXT.  The D-Vector is written into the include file
c	FSL_MODEL for use by fsl.
c
c	Author:   
c		  FCF_Read_Dvector
c                 Gene Eplee
c		  General Sciences Corp.
c		  25 October 1993
c		  SER 11397
c       
c                FSL_Read_Dvector
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
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
c		fut_free_lun
c		fut_get_lun
c		lib$movc5
c		lib$signal
c		str$trim
c
c	Include files:
c		fsl_invoc.txt
c		fsl_model.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, August 7, 1995 
c       Modified FCF_Read_Dvector to FSL_Read_Dvector for the new FIRAS 
c       pipeline which will process long spectra to get improved frequency 
c       resolution. 
c           1. An array of 32 real*8 bolometer parameters and arrays of 5 
c              257 point complex*16 instrument emissivities comprise the 
c              calibration model solution. Since model array sizes for the 
c              new pipeline are unchanged, the array sizes for the spectral 
c              variances are limited to 257 by the restrictions to the record 
c              length of the model. Thus the length of the dvector is 
c              unchanged for FSL.
c           2. Changed status, include file, and function names for FSL.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fsl_invoc)'
	include '(fsl_model)'

	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status
	integer * 4	var_lun		!  D-Vector file lun

	logical * 1	var_ver		!  D-Vector verification flag

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fex_var'
	record /fex_var/ fex_var

	external	fsl_varclose
	external	fsl_varinval
	external	fsl_varopen
	external	fsl_varread
	external	fsl_normal
	external	fut_normal

	status = %loc(fsl_normal)

C
C  Open, read, and close the D-Vector file.
C

c
c  Find the D-Vector filename.
c
	fcc_var_file = 'CSDR$FIRAS_CAL:FEX_VAR_' // fcc_scan_mode // '.' //
     .			fcc_model_ext
	call str$trim (fcc_var_file, fcc_var_file, fcc_varlen)

c
c  Open the D-Vector file.
c
	rstatus = fut_get_lun(var_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=var_lun, file=fcc_var_file, status='old',
     .	      form='unformatted', recordtype='fixed', recl=fcc_var_size,
     .	      access='sequential', readonly, iostat=io_stat)

	if (io_stat .eq. 0) then
c
c  Read the D-Vector.
c
	   read (var_lun, iostat=io_stat) fex_var
	   if (io_stat .ne. 0) then
	      status = %loc(fsl_varread)
	      call lib$signal (fsl_varread, %val(2),
     .			       fcc_var_file(1:fcc_varlen), %val(io_stat))
	   endif

c
c  Close the D-Vector file.
c
	   close (unit=var_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fsl_varclose)
	      call lib$signal (fsl_varclose, %val(2),
     .			       fcc_var_file(1:fcc_varlen), %val(io_stat))
	   endif

	   rstatus = fut_free_lun(var_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	else
	   status = %loc(fsl_varopen)
	   call lib$signal (fsl_varopen, %val(2),
     .			    fcc_var_file(1:fcc_varlen), %val(io_stat))
	endif	!  (open status


	if (status .eq. %loc(fsl_normal)) then
C
C  Process the D-Vector.
C

c
c  Verify that the correct D-Vector has been read.
c
	   var_ver = .true.
	   if (fcc_chan .ne. fex_var.channel) var_ver = .false.
	   if (fcc_length .ne. fex_var.scan_length) var_ver = .false.
	   if (fcc_speed .ne. fex_var.scan_speed) var_ver = .false.
	   if (var_ver .eq. .false.) then
	      status = %loc(fsl_varinval)
	      call lib$signal (fsl_varinval, %val(1),
     .			       fcc_var_file(1:fcc_varlen))
	   endif

	endif

	if (status .eq. %loc(fsl_normal)) then
c
c  Insert the D-Vector into the common block.
c
	   call lib$movc5 (0,,0,2056,dvector)
	   do j = 1,257
	      dvector(j) = fex_var.variances(j)
	   enddo 
	endif


	fsl_read_dvector = status

	return
	end
