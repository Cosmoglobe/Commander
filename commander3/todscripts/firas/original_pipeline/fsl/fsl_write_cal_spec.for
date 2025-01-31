	integer * 4 function  fsl_write_cal_spec (flag, cs_lun, num, spec_recs)

c-------------------------------------------------------------------------------
c
c	Function FSL_WRITE_CAL_SPEC
c
c	This function writes calibration voltage spectrum records and calibrated
c	spectrum records out to the Cobetrieve time-ordered archives. The 
c       spectra are written out by the Cobetrieve routine CT_Connect_Write.
c
c	Author:	  
c                FCF_Write_Cal_Spec
c                Gene Eplee
c	         General Sciences Corp.
c	         21 April 1992
c         
c                FSL_Write_Cal_Spec
c                Shirley M. Read
c                Hughes STX Corporation
c                July 1995
c
c-------------------------------------------------------------------------------
c
c       Input:
c               flag            integer * 4             spectrum type flag
c                                                       1 = voltage
c                                                       2 = differential
c                                                       3 = calibrated
c		cs_lun		integer * 4		spectra file lun
c               num             integer * 4             number of spectra in
c                                                       current ensemble 
c               spec_recs(num)                          spectrum records

c
c	Output:
c		none
c
c	Subroutines called:
c		ct_write_arcv
c		lib$signal
c
c	Include files:
c		ct$library:ctuser.inc
c               fut_params
c               fsl_invoc
c
c-------------------------------------------------------------------------------
c
c	Changes for FCF:
c
c	Put in optional differential spectrum output.
c	Gene Eplee, GSC, 11 July 1994
c	SPR 11826
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, July 26, 1995 
c       Modified FCF_Write_Cal_Spec to FSL_Write_Cal_Spec for the new FIRAS 
c       pipeline which will process long spectra to get improved frequency 
c       resolution.
c           1. Changed status and function names.
c           2. Changed calling sequence to input spectra file logical unit,
c              number of records, and spectrum records, 
c           3. Removed open of Cobetrieve archive and RMS spectrum files.
c           4. Removed close of Cobetrieve archive and RMS spectrum files.
c           5. Removed read of RMS files and related temporary file parameters.
c           
c-------------------------------------------------------------------------------

	implicit none

	include 'ct$library:ctuser.inc'
	include '(fut_params)'
	include '(fsl_invoc)'

	integer * 4     flag            !  spectrum type flag
	integer * 4	cs_lun		!  CT logical unit number for cal spec
	integer * 4     num             !  number of spectra in ensemble

	character * 48  out_file        !  Cobetrieve file name

	integer * 2     out_len         !  length of skymap file name

	integer * 4	numout		!  number of spectra written out
	integer * 4	status		!  return status
	integer * 2	ct_stat(20)	!  CT return status

	dictionary 'fsl_sky'
	record /fsl_sky/ spec_recs(num)

	external	fsl_ctspecwrite
	external	fsl_normal

C
C  Initialize function return.
C

	status = %loc(fsl_normal)

C
C  Determine the file name for message.
C
	if (flag .eq. 1) then
	   out_file  = fcc_outfile_vs
	   out_len   = fcc_outlenvs
	elseif (flag .eq. 2) then
	   out_file  = fcc_outfile_ds
	   out_len   = fcc_outlends
	elseif (flag .eq. 3) then
	   out_file  = fcc_outfile_cs
	   out_len   = fcc_outlencs
	endif

C
C  Write the spectrum records to the Cobetrieve file.
C

	ct_stat(1) = ctp_normal
	numout = 1

	do while ((numout .le. num)  .and.  (ct_stat(1) .eq. ctp_normal))

	      spec_recs(numout).coad_spec_head.coadd_no = numout
	      call ct_write_arcv (, cs_lun, spec_recs(numout), ct_stat)

	      numout = numout + 1

	end do

	if (ct_stat(1) .ne. ctp_normal) then
	      status = %loc(fsl_ctspecwrite)
	      call lib$signal (fsl_ctspecwrite, %val(2), out_file(1:out_len),
     .						   %val(ct_stat(1)))
	end if


	fsl_write_cal_spec = status

	return
	end
