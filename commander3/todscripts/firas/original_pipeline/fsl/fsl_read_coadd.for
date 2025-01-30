	integer * 4 function  fsl_read_coadd (ct_lun, num, coadd_recs)

c-------------------------------------------------------------------------------
c
c	Function FSL_READ_COADD
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
c	Author:	 
c                FCF_Read_Coadd
c		 Gene Eplee
c		 General Sciences Corp.
c		 5 June 1992
c
c       
c                FSL_Read_Coadd
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
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
c		fsl_read_cal_coadd
c		fsl_read_sky_coadd
c		fsl_read_sky_coadd_list
c		fut_free_lun
c		lib$signal
c
c	Include files:
c		ct$library:ctuser.inc
c		fsl_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, August 7, 1995 
c       Modified FCF_Read_Coadd to FSL_Read_Coadd for the new FIRAS pipeline 
c       which will process long spectra to get improved frequency resolution.
c           1. Changed status, include file, and function names for FSL.
c           2. Changed dictionary and record from FIC to FIL.
c
c-------------------------------------------------------------------------------

	implicit none

	include 'ct$library:ctuser.inc'
	include '(fut_params)'
	include '(fsl_invoc)'

	integer * 2	ct_stat(20)	!  CT return status

	integer * 4	cstatus		!  return status
	integer * 4	ct_lun		!  CT logical unit number
	integer * 4	num		!  number of coadd records
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	csa_close_skymap
	integer * 4	fsl_read_cal_coadd
	integer * 4	fsl_read_sky_coadd
	integer * 4	fsl_read_sky_coadd_list
	integer * 4	fut_free_lun

	dictionary 'fil_sky'
	record /fil_sky/ coadd_recs(fcc_max_coadd)

	external	csa_normal
	external	fsl_csaifgclose
	external	fsl_ctifgclose
	external	fsl_eof
	external	fsl_normal
	external	fut_normal

C
C  Read in the coadded IFGs.
C

	if (fcc_sky .eq. fac_present) then
c
c  Read the sky coadds.
c
	   if (fcc_pixel .eq. fac_present) then
	      status = fsl_read_sky_coadd_list (ct_lun, num, coadd_recs)
	   else
	      status = fsl_read_sky_coadd (ct_lun, num, coadd_recs)
	   endif

	else
c
c  Read the calibration coadds.
c
	   status = fsl_read_cal_coadd (ct_lun, num, coadd_recs)

	endif


	if (status .eq. %loc(fsl_eof)) then
C
C  Close the archive after the final read.
C

	   if (fcc_sky .eq. fac_present) then
c
c  Close the sky coadd archive.
c
	      cstatus = csa_close_skymap (ct_lun, fac_skymap_no_levels)
	      if (cstatus .ne. %loc(csa_normal)) then
	         status = %loc(fsl_csaifgclose)
	         call lib$signal (fsl_csaifgclose, %val(2),
     .				  fcc_infile(1:fcc_inlen), %val(cstatus))
	      end if

	   else
c
c  Close the calibration coadd archive.
c
	      call ct_close_arcv (, ct_lun, ct_stat)
	      if (ct_stat(1) .ne. ctp_normal) then
	         status = %loc(fsl_ctifgclose)
	         call lib$signal (fsl_ctifgclose, %val(2),
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


	fsl_read_coadd = status

	return
	end
