	integer * 4 function  fsl_open_coadd (ct_lun)

c-------------------------------------------------------------------------------
c
c	Function FSL_OPEN_COADD
c
c	This function opens the coadd archive CSDR$FIRAS_IN for either sky or
c	cal data.  It gets the reference datasets from the reference archive
c	CSDR$FIRAS_REF using CCT_Get_Config routines, calls an FSL routine to
c       compute constants required for processing the spectra, reads the 
c       calibration model solution from the reference archive CSDR$FIRAS_CAL. 
c       and optionally reads the D_Vector or averaged FIL variances from the 
c       reference archive CSDR$FIRAS_UREF. 
c       
c
c	Author:	 
c                FCF_Open_Coadd
c                Gene Eplee
c		 General Sciences Corp.
c		 14 July 1992
c         
c                FSL_Open_Coadd
c                Shirley M. Read
c                Hughes STX Corporation
c                July 1995
c
c-------------------------------------------------------------------------------
c
c	Input:
c		No parameters passed. Input is obtained from include files.
c
c	Output:
c		ct_lun		integer * 4		coadd CT file lun
c
c	Subroutines called:
c		cct_get_config_idx_tod
c		cct_get_config_tod
c		fsl_compute_constants
c		fsl_open_cal_coadd
c		fsl_open_sky_coadd
c		fsl_read_dvector
c		fsl_read_flv
c		fsl_read_model
c		fut_get_lun
c		lib$signal
c
c	Include files:
c		fsl_config.txt
c		fsl_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes for FCF:
c
c	Added call to FCF_READ_DVECTOR.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11397
c
c	Put in optional differential spectrum output.
c	Gene Eplee, GSC, 11 July 1994
c	SPR 11826
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, July 25, 1995 
c       Modified FCF_Open_Coadd to FSL_Open_Coadd for the new FIRAS pipeline 
c       which will process long spectra to get improved frequency resolution.
c           1. Changed status and function names.
c           2. Removed the open of temporary RMS files.
c           3. Added call to FSL_Read_FLV.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fsl_config)'
	include '(fsl_invoc)'

	integer * 4	cstatus		!  return status
	integer * 4	ct_lun		!  CT logical unit number
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	cct_get_config_idx_tod
	integer * 4	cct_get_config_tod
	integer * 4	fsl_compute_constants
	integer * 4	fsl_open_cal_coadd
	integer * 4	fsl_open_sky_coadd
	integer * 4	fsl_read_dvector
	integer * 4	fsl_read_flv
	integer * 4	fsl_read_model
	integer * 4	fut_get_lun

	external	cct_normal
	external	fsl_cfgdirget
	external	fsl_cfgseqget
	external	fsl_normal
	external	fut_normal

C
C  Open the coadd file.
C
c
c  Get the logical unit number for Cobetrieve.
c
	rstatus = fut_get_lun (ct_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

c
c  Open the coadd file for direct file access.
c
	if (fcc_sky .eq. fac_present) then
	   status = fsl_open_sky_coadd (ct_lun)
	else
	   status = fsl_open_cal_coadd (ct_lun)
	endif


	if (status .eq. %loc(fsl_normal)) then
C
C  Get the reference datasets.
C

	   cstatus = cct_get_config_tod (fcc_jstart, ndset_tod, size_tod,
     .					 tod_lun, tod_index, config,
     .				         new_tod_segment, tod_stat)

	   if (cstatus .ne. %loc(cct_normal)) then
	      status = %loc(fsl_cfgseqget)
	      call lib$signal (fsl_cfgseqget, %val(1), %val(cstatus))
	   else
	      cstatus = cct_get_config_idx_tod (fcc_jstart, ndset_dir,
     .						dir_lun, dir_index,
     .						new_dir_segment, dir_stat)
	      if (cstatus .ne. %loc(cct_normal)) then
	         status = %loc(fsl_cfgdirget)
	         call lib$signal (fsl_cfgdirget, %val(1), %val(cstatus))
	      endif
	   endif

	   if (status .eq. %loc(fsl_normal)) then
c
c  Compute the constants required for processing the spectra.
c
	      status = fsl_compute_constants ()
	   endif


	   if ((status .eq. %loc(fsl_normal))  .and.
     .	       (fcc_calibrate .eq. fac_present)) then
C
C  Read the calibration model solution and the D-Vector file or the FIL
C  averaged variances file.
C
	      status = fsl_read_model ()
	      if ((status .eq. %loc(fsl_normal))  .and.
     .		  (fcc_dvec .eq. fac_present)) then
	         status = fsl_read_dvector ()
	      endif

	      if ((status .eq. %loc(fsl_normal))  .and.
     .		  (fcc_flv .eq. fac_present)) then
	         status = fsl_read_flv ()
	      endif

	   endif    !   (calibrate flag is set

	endif	!	(status from coadd file open


	fsl_open_coadd = status

	return
	end
