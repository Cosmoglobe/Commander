	integer * 4 function  fip_read_model ()

c-------------------------------------------------------------------------------
c
c	Function FIP_READ_MODEL
c
c	This function opens, reads, and closes the FIRAS calibration model
c	solution FEX_MOD_CCSS.VVV_XXXXXXXXXX from the directory CSDR$FIRAS_CAL 
c	(CC is channel, SS is scan mode, VVV is the version of the calibration
c	program that generated the solution, and XXXXXXXXXX is the model
c	solution label).  The file extension is specified by the command line
c	qualifier /MODEL_EXT.
c
c	The function also reads the C-Vector^2 from the reference dataset
c	CSDR$FIRAS_IN:FEX_CVS_CCSS.VVV_XXXXXXXXXX.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  19 May 1993
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
c		lib$signal
c		str$trim
c
c	Include files:
c		fip_invoc_model.txt
c		fip_model.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Read C-Vector^2 from FEX_CVS_CCSS reference dataset.
c	   Gene Eplee, GSC, 29 September 1993.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fip_invoc_model)'
	include '(fip_model)'

	integer * 4	cvs_lun		!  C-vector file lun
	integer * 4	io_stat		!  I/O return status
	integer * 4	mod_lun		!  model solution file lun
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	external	fip_normal
	external	fip_rmsclose
	external	fip_rmsopen
	external	fip_rmsread
	external	fut_normal

	status = %loc(fip_normal)

C
C  Open, read, and close the calibration model solution file.
C

c
c  Find the model solution file name.
c
	fcc_mod_file = 'CSDR$FIRAS_CAL:FEX_MOD_' // fcc_scan_mode // '.' //
     .			  fcc_model_ext
	call str$trim (fcc_mod_file, fcc_mod_file, fcc_modlen)

c
c  Open the model solution file.
c
	rstatus = fut_get_lun(mod_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=mod_lun, file=fcc_mod_file, status='old',
     .	      form='unformatted', recordtype='fixed', recl=fcc_mod_size,
     .	      access='sequential', readonly, iostat=io_stat)

	if (io_stat .eq. 0) then
c
c  Read the model solution.
c
	   read (mod_lun, iostat=io_stat) fex_model
	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmsread)
	      call lib$signal (fip_rmsread, %val(2),
     .			       fcc_mod_file(1:fcc_modlen), %val(io_stat))
	   endif

c
c  Close the model solution file.
c
	   close (unit=mod_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmsclose)
	      call lib$signal (fip_rmsclose, %val(2),
     .			       fcc_mod_file(1:fcc_modlen), %val(io_stat))
	   endif

	   rstatus = fut_free_lun(mod_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	else
	   status = %loc(fip_rmsopen)
	   call lib$signal (fip_rmsopen, %val(2),
     .			    fcc_mod_file(1:fcc_modlen), %val(io_stat))
	endif	!  (open status


	if (status .eq. %loc(fip_normal)) then
C
C  Read the C-vector^2 from the reference dataset.
C

c
c  Find the FEX_CVS file name.
c
	   fcc_cvs_file = 'CSDR$FIRAS_IN:FEX_CVS_' // fcc_scan_mode // '.'
     .			  // fcc_model_ext
	   call str$trim (fcc_cvs_file, fcc_cvs_file, fcc_cvslen)

c
c  Open the FEX_CVS file.
c
	   rstatus = fut_get_lun(cvs_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	   open (unit=mod_lun, file=fcc_cvs_file, status='old',
     .		 form='unformatted', recordtype='fixed', recl=fcc_cvs_size,
     .		 access='sequential', readonly, iostat=io_stat)

	   if (io_stat .eq. 0) then
c
c  Read the C-Vector^2.
c
	      read (cvs_lun, iostat=io_stat) fex_cvs
	      if (io_stat .ne. 0) then
	         status = %loc(fip_rmsread)
	         call lib$signal (fip_rmsread, %val(2),
     .				  fcc_cvs_file(1:fcc_cvslen), %val(io_stat))
	      endif

c
c  Close the FEX_CVS file.
c
	      close (unit=cvs_lun, iostat=io_stat)
	      if (io_stat .ne. 0) then
	         status = %loc(fip_rmsclose)
	         call lib$signal (fip_rmsclose, %val(2),
     .				  fcc_cvs_file(1:fcc_cvslen), %val(io_stat))
	      endif

	      rstatus = fut_free_lun(cvs_lun)
	      if (rstatus .ne. %loc(fut_normal)) then
	         call lib$signal (%val(rstatus))
	      endif

	   else
	      status = %loc(fip_rmsopen)
	      call lib$signal (fip_rmsopen, %val(2),
     .			       fcc_cvs_file(1:fcc_cvslen), %val(io_stat))
	   endif	!  (open status

	endif		!   status from model read


	fip_read_model = status

	return
	end
