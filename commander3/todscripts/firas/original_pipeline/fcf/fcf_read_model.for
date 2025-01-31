	integer * 4 function  fcf_read_model ()

c-------------------------------------------------------------------------------
c
c	Function FCF_READ_MODEL
c
c	This function opens, reads, and closes the FIRAS calibration model
c	solution FEX_MOD_CCSS.VVV_XXXXXXXXXX from the directory CSDR$FIRAS_CAL 
c	(CC is channel, SS is scan mode, VVV is the version of the calibration
c	program that generated the solution, and XXXXXXXXXX is the model
c	solution label).  The file extension is specified by the command line
c	qualifier /MODEL_EXT.  The model solution is written into the include
c	file FCF_MODEL for use by FCF.  The routine optionally calls
c	FCF_DISPLAY_MODEL to plot the model solution.
c
c	The routine switches the order of the indices of the emissivity array
c	to improve memory management for computational purposes.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  17 June 1992
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
c		fcf_display_model
c		fut_free_lun
c		fut_get_lun
c		lib$movc5
c		lib$signal
c		str$trim
c
c	Include files:
c		fcf_invoc.txt
c		fcf_model.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Removed undefined variable MODEL_TYPE from error trapping.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11395
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fcf_invoc)'
	include '(fcf_model)'

	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	k		!  a counter
	integer * 4	mod_lun		!  model solution file lun
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	logical * 1	mod_ver		!  model verification flag

	integer * 4	fcf_display_model
	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fex_mod'
	record /fex_mod/ fex_model

	external	fcf_modelclose
	external	fcf_modelinval
	external	fcf_modelopen
	external	fcf_modelread
	external	fcf_normal
	external	fut_normal

	status = %loc(fcf_normal)

C
C  Open, read, and close the calibration model solution file.
C

c
c  Find the model solution file name.
c
	fcc_model_file = 'CSDR$FIRAS_CAL:FEX_MOD_' // fcc_scan_mode // '.' //
     .			  fcc_model_ext
	call str$trim (fcc_model_file, fcc_model_file, fcc_modlen)

c
c  Open the model solution file.
c
	rstatus = fut_get_lun(mod_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=mod_lun, file=fcc_model_file, status='old',
     .	      form='unformatted', recordtype='fixed', recl=fcc_model_size,
     .	      access='sequential', readonly, iostat=io_stat)

	if (io_stat .eq. 0) then
c
c  Read the model solution.
c
	   read (mod_lun, iostat=io_stat) fex_model
	   if (io_stat .ne. 0) then
	      status = %loc(fcf_modelread)
	      call lib$signal (fcf_modelread, %val(2),
     .			       fcc_model_file(1:fcc_modlen), %val(io_stat))
	   endif

c
c  Close the model solution file.
c
	   close (unit=mod_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fcf_modelclose)
	      call lib$signal (fcf_modelclose, %val(2),
     .			       fcc_model_file(1:fcc_modlen), %val(io_stat))
	   endif

	   rstatus = fut_free_lun(mod_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	else
	   status = %loc(fcf_modelopen)
	   call lib$signal (fcf_modelopen, %val(2),
     .			    fcc_model_file(1:fcc_modlen), %val(io_stat))
	endif	!  (open status


	if (status .eq. %loc(fcf_normal)) then
C
C  Verify that the correct model solution has been read.
C
	   mod_ver = .true.
	   if (fcc_chan .ne. fex_model.mod_head.chan_id) mod_ver = .false.
	   if (fcc_length .ne. fex_model.mod_head.mtm_length) mod_ver = .false.
	   if (fcc_speed .ne. fex_model.mod_head.mtm_speed) mod_ver = .false.
	   if (fcc_ngroup.ne.fex_model.mod_head.adds_per_group) mod_ver=.false.
	   if (mod_ver .eq. .false.) then
	      status = %loc(fcf_modelinval)
	      call lib$signal (fcf_modelinval, %val(1),
     .			       fcc_model_file(1:fcc_modlen))
	   endif

	endif


	if (status .eq. %loc(fcf_normal)) then
C
C  Insert the model solution into the common block.
C

c
c  Initialize the common block.
c
	   call lib$movc5 (0,,0,8*npar,param)
	   call lib$movc5 (0,,0,80,bolparm)
	   call lib$movc5 (0,,0,28784,emiss)

c
c  Insert the model solution id.
c
	   model_label = fex_model.mod_head.label
	   call str$trim (model_label, model_label, mod_lablen)
	   model_ttag  = fex_model.mod_head.gmt

c
c  Insert the bolometer parameters.
c
	   do j = 1,npar
	      param(j) = fex_model.bolparm(j)
	   enddo
	   do j = 1,9
	      bolparm(j) = fex_model.bolparm(j)
	   enddo
	   bolparm(10) = load_resistance

c
c  Insert the emissivities, switching the order of the indices for improved 
c	memory management.
c
	   do k = 1,257
	      do j = 1,7
	         emiss(j,k) = fex_model.emissivity(k,j)
	      enddo
	   enddo

	endif	!	(status from reading model


C
C  Optionally, plot the model.
C
	if ((status .eq. %loc(fcf_normal))  .and.
     .	    (fcc_plot .eq. fac_present)) then
	   status = fcf_display_model ()
	endif


	if (status .eq. %loc(fcf_normal)) type 10, model_label(1:mod_lablen)
 10	format (x, 'Retrieved model solution  ', a, /)


	fcf_read_model = status

	return
	end
