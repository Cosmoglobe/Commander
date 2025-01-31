	integer * 4 function  fcf_read_coadd (ct_lun, num, coadd_recs)

c-------------------------------------------------------------------------------
c
c	Function FCF_READ_COADD
c
c	This function reads in coadded IFGs from the Cobetrieve sky or 
c	calibration archive.  Each time through, the function reads one
c	populated pixel from the sky archive or up to FCC_MAX_COADD coadd
c	records from the cal archive.  The last time through, the function
c	closes the Cobetrieve archive.
c
c	If requested, the function reads data for sky pixels specified in the
c	array FCC_PLIST.
c
c	Author:	 Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 5 June 1992
c
c-------------------------------------------------------------------------------
c
c	Input:
c		ct_lun		integer * 4		coadd file lun	
c
c	Output:
c		num		integer * 4		number of coadds read
c		coadd_recs	coadd records		input coadd records
c
c	Subroutines called:
c		csa_close_skymap
c		ct_close_arcv
c		fcf_read_cal_coadd
c		fcf_read_sky_coadd
c		fcf_read_sky_coadd_list
c		fut_free_lun
c		lib$signal
c
c	Include files:
c		ct$library:ctuser.inc
c		fcf_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include 'ct$library:ctuser.inc'
	include '(fut_params)'
	include '(fcf_invoc)'

	integer * 2	ct_stat(20)	!  CT return status

	integer * 4	cstatus		!  return status
	integer * 4	ct_lun		!  CT logical unit number
	integer * 4	num		!  number of coadd records
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	csa_close_skymap
	integer * 4	fcf_read_cal_coadd
	integer * 4	fcf_read_sky_coadd
	integer * 4	fcf_read_sky_coadd_list
	integer * 4	fut_free_lun

	dictionary 'fic_sky'
	record /fic_sky/ coadd_recs(fcc_max_coadd)

	external	csa_normal
	external	fcf_csaifgclose
	external	fcf_ctifgclose
	external	fcf_eof
	external	fcf_normal
	external	fut_normal

C
C  Read in the coadded IFGs.
C

	if (fcc_sky .eq. fac_present) then
c
c  Read the sky coadds.
c
	   if (fcc_pixel .eq. fac_present) then
	      status = fcf_read_sky_coadd_list (ct_lun, num, coadd_recs)
	   else
	      status = fcf_read_sky_coadd (ct_lun, num, coadd_recs)
	   endif

	else
c
c  Read the calibration coadds.
c
	   status = fcf_read_cal_coadd (ct_lun, num, coadd_recs)

	endif


	if (status .eq. %loc(fcf_eof)) then
C
C  Close the archive after the final read.
C

	   if (fcc_sky .eq. fac_present) then
c
c  Close the sky coadd archive.
c
	      cstatus = csa_close_skymap (ct_lun, fac_skymap_no_levels)
	      if (cstatus .ne. %loc(csa_normal)) then
	         status = %loc(fcf_csaifgclose)
	         call lib$signal (fcf_csaifgclose, %val(2),
     .				  fcc_infile(1:fcc_inlen), %val(cstatus))
	      end if

	   else
c
c  Close the calibration coadd archive.
c
	      call ct_close_arcv (, ct_lun, ct_stat)
	      if (ct_stat(1) .ne. ctp_normal) then
	         status = %loc(fcf_ctifgclose)
	         call lib$signal (fcf_ctifgclose, %val(2),
     .				  fcc_infile(1:fcc_inlen), %val(cstatus))
	      end if

	   endif

c
c  Free the Cobetrieve logical unit number.
c
	   rstatus = fut_free_lun (ct_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	endif		!  (end of file


	fcf_read_coadd = status

	return
	end
