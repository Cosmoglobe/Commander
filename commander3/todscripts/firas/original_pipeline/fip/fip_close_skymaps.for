	integer * 4 function  fip_close_skymaps (flun, alun)

c-------------------------------------------------------------------------------
c
c	Function FIP_CLOSE_SKYMAPS
c
c	This function closes the FIRAS input spectrum skymap and the ADB output
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
c		alun		integer * 4		ADB output skymap lun
c		flun		integer * 4		FIRAS input skymap lun
c
c	Output:
c		none
c
c	Subroutines called:
c		csa_close_skymap
c		fut_free_lun
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

	integer * 4	csa_close_skymap
	integer * 4	fut_free_lun

	external	csa_normal
	external	fip_csaclose
	external	fip_normal
	external	fut_normal


	status = %loc(fip_normal)

C
C  Close the files.
C

c
c  Close the input FIRAS skymap file.
c
	cstatus = csa_close_skymap (flun, fac_skymap_no_levels)
	if (cstatus .ne. %loc(csa_normal)) then
	   status = %loc(fip_csaclose)
	   call lib$signal (fip_csaclose, %val(2), fcc_infile(1:fcc_inlen),
     .					  %val(cstatus))
	endif

c
c  Close the output ADB skymap file.
c
	cstatus = csa_close_skymap (alun, fac_skymap_no_levels)
	if (cstatus .ne. %loc(csa_normal)) then
	   status = %loc(fip_csaclose)
	   call lib$signal (fip_csaclose, %val(2), fcc_outfile(1:fcc_outlen),
     .					  %val(cstatus))
	endif

c
c  Free the logical unit numbers.
c
	rstatus = fut_free_lun(flun)
	rstatus = fut_free_lun(alun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif


	fip_close_skymaps = status

	return
	end
