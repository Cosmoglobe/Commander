	integer * 4 function  fcf_open_coadd (ct_lun, vs_lun, cs_lun)

c-------------------------------------------------------------------------------
c
c	Function FCF_OPEN_COADD
c
c	This function opens the coadd archive CSDR$FIRAS_IN for either sky or
c	cal data.  It gets the reference datasets from the reference archive
c	CSDR$FIRAS_REF, then reads the calibration model solution from the
c	reference archive CSDR$FIRAS_CAL.  If required, it opens the temporary
c	spectrum RMS files.
c
c	Author:	 Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 14 July 1992
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none	
c
c	Output:
c		cs_lun		integer * 4		calibrated spectra RMS
c							   file lun
c		ct_lun		integer * 4		coadd CT file lun
c		vs_lun		integer * 4		voltage spectra RMS
c							   file lun
c
c	Subroutines called:
c		cct_get_config_idx_tod
c		cct_get_config_tod
c		fcf_compute_constants
c		fcf_open_cal_coadd
c		fcf_open_sky_coadd
c		fcf_open_temp_spec
c		fcf_read_dvector
c		fcf_read_model
c		fut_get_lun
c		lib$signal
c
c	Include files:
c		fcf_config.txt
c		fcf_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Added call to FCF_READ_DVECTOR.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11397
c
c	Put in optional differential spectrum output.
c	Gene Eplee, GSC, 11 July 1994
c	SPR 11826
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fcf_config)'
	include '(fcf_invoc)'

	integer * 4	cs_lun		!  calibrated spectra RMS file lun
	integer * 4	cstatus		!  return status
	integer * 4	ct_lun		!  CT logical unit number
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status
	integer * 4	vs_lun		!  voltage spectra RMS file lun

	integer * 4	cct_get_config_idx_tod
	integer * 4	cct_get_config_tod
	integer * 4	fcf_compute_constants
	integer * 4	fcf_open_cal_coadd
	integer * 4	fcf_open_sky_coadd
	integer * 4	fcf_open_temp_spec
	integer * 4	fcf_read_dvector
	integer * 4	fcf_read_model
	integer * 4	fut_get_lun

	external	cct_normal
	external	fcf_cfgdirget
	external	fcf_cfgseqget
	external	fcf_normal
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
	   status = fcf_open_sky_coadd (ct_lun)
	else
	   status = fcf_open_cal_coadd (ct_lun)
	endif


	if (status .eq. %loc(fcf_normal)) then
C
C  Get the reference datasets.
C

	   cstatus = cct_get_config_tod (fcc_jstart, ndset_tod, size_tod,
     .					 tod_lun, tod_index, config,
     .				         new_tod_segment, tod_stat)

	   if (cstatus .ne. %loc(cct_normal)) then
	      status = %loc(fcf_cfgseqget)
	      call lib$signal (fcf_cfgseqget, %val(1), %val(cstatus))
	   else
	      cstatus = cct_get_config_idx_tod (fcc_jstart, ndset_dir,
     .						dir_lun, dir_index,
     .						new_dir_segment, dir_stat)
	      if (cstatus .ne. %loc(cct_normal)) then
	         status = %loc(fcf_cfgdirget)
	         call lib$signal (fcf_cfgdirget, %val(1), %val(cstatus))
	      endif
	   endif

	   if (status .eq. %loc(fcf_normal)) then
c
c  Compute the constants required for processing the spectra.
c
	      status = fcf_compute_constants ()
	   endif


	   if ((status .eq. %loc(fcf_normal))  .and.
     .	       (fcc_calibrate .eq. fac_present)) then
C
C  Read the calibration model solution and the D-Vector file.
C
	      status = fcf_read_model ()
	      if ((status .eq. %loc(fcf_normal))  .and.
     .		  (fcc_dvec .eq. fac_present)) then
	         status = fcf_read_dvector ()
	      endif

	   endif

	endif	!	(status from coadd file open


C
C  Open the temporary RMS spectrum files.
C

c
c  The voltage spectrum file.
c
	vs_lun = 0
	if ((status .eq. %loc(fcf_normal))  .and.
     .	    (fcc_write_vs .eq. fac_present)) then
	   status = fcf_open_temp_spec (1, vs_lun)
	endif

c
c  The calibrated or differential spectrum file.
c
	cs_lun = 0
	if (status .eq. %loc(fcf_normal)) then
	   if (fcc_write_ds .eq. fac_present) then
	      status = fcf_open_temp_spec (2, cs_lun)
	   elseif (fcc_write_cs .eq. fac_present) then
	      status = fcf_open_temp_spec (3, cs_lun)
	   endif
	endif


	fcf_open_coadd = status

	return
	end
