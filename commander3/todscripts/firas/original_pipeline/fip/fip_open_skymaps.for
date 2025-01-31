	integer * 4 function  fip_open_skymaps (flun, alun)

c-------------------------------------------------------------------------------
c
c	Function FIP_OPEN_SKYMAPS
c
c	This function opens the FIRAS input spectrum skymap and the ADB output
c	spectrum skymap.
c
c	Author:  Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 10 May 1993
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		alun		integer * 4		ADB input skymap lun
c		flun		integer * 4		FIRAS input skymap lun
c
c	Subroutines called:
c		csa_field_offset_values
c		csa_open_skymap
c		fut_get_lun
c		lib$signal
c
c	Include files:
c		fip_invoc_sky.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Move filename determination to the parsing routine.
c	Gene Eplee, GSC, 10 February 1994.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fip_invoc_sky)'

	integer * 4	alun			!  ADB skymap lun
	integer * 4	cstatus			!  CT or CSA return status
	integer * 4	flun			!  FIRAS skymap lun
	integer * 4	io_stat			!  I/O return status
	integer * 4	rstatus			!  return status
	integer * 4	status			!  return status

	integer * 4	csa_field_offset_values
	integer * 4	csa_open_skymap
	external	csa_open_skymap
	integer * 4	fut_get_lun

	external	csa_normal
	external	fip_csafld
	external	fip_csaopen
	external	fip_normal
	external	fut_normal


C
C  Initialize the open.
C
	status = %loc(fip_normal)

c
c  Get the logical unit numbers for the files.
c
	rstatus = fut_get_lun (flun)
	rstatus = fut_get_lun (alun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal(rstatus)
	endif


C
C  Open the files.
C

c
c  Open the input FIRAS skymap file.
c
	open (unit=flun, file=fcc_infile, iostat=io_stat, status='old',
     .	      form='unformatted', recordtype='fixed', readonly,
     .	      useropen=csa_open_skymap)

	if (io_stat .eq. 0) then
c
c  Open the output ADB skymap file.
c
	   open (unit=alun, file=fcc_outfile, iostat=io_stat, status='new',
     .		 form='unformatted', recordtype='fixed', recl=fcc_alen,
     .		 useropen=csa_open_skymap)

	   if (io_stat .eq. 0) then
c
c  Set up the skymap offset values.
c
	      cstatus = csa_field_offset_values (fcc_pix_offset, -1, -1, alun)
	      if (cstatus .ne. %loc(csa_normal)) then
	         status = %loc(fip_csafld)
	         call lib$signal (fip_csafld, %val(2),
     .				  fcc_outfile(1:fcc_outlen), %val(cstatus))
	      endif
	   else
 	      status = %loc(fip_csaopen)
	      call lib$signal (fip_csaopen, %val(2), fcc_outfile(1:fcc_outlen),
     .					    %val(io_stat))
	   endif

	else
 	   status = %loc(fip_csaopen)
	   call lib$signal (fip_csaopen, %val(2), fcc_infile(1:fcc_inlen),
     .					 %val(io_stat))
	endif


	fip_open_skymaps = status

	return
	end
