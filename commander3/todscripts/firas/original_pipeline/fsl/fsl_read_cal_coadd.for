	integer * 4 function   fsl_read_cal_coadd (ct_lun, num, coadd_recs)

c-------------------------------------------------------------------------------
c
c	Function FSL_READ_CAL_COADD
c
c	This function reads the calibration coadd records from the Cobetrieve
c	archive and filters them by the data quality summary flags. 
c	FCC_MAX_COAD records are returned to the calling routine.
c
c	Author:	  
c                FCF_Read_Cal_Coadd
c                Gene Eplee
c		 General Sciences Corp.
c		 21 April 1992
c       
c                FSL_Read_Cal_Coadd
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
c
c-------------------------------------------------------------------------------
c
c	Input:
c		ct_lun		integer * 4		CT logical unit number
c
c	Output:
c		num		integer * 4		number of coadds read
c		coadd_recs	coadd records		input coadd records
c
c	Subroutines called:
c		ct_read_arcv
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
c       Modified FCF_Read_Cal_Coadd to FSL_Read_Cal_Coadd for the new FIRAS 
c       pipeline which will process long spectra to get improved frequency 
c       resolution.
c           1. Changed status, include file, and function names for FSL.
c           2. Changed dictionary and record from FIC to FIL.
c
c-------------------------------------------------------------------------------

	implicit none

	include 'ct$library:ctuser.inc'
	include '(fut_params)'
	include '(fsl_invoc)'

	integer * 2	ct_stat(20)	!  ct return status

	integer * 4	ct_lun		!  ct logical unit number
	integer * 4	num		!  number of coadds read
	integer * 4	status		!  return status

	dictionary 'fil_sky'
	record /fil_sky/ coadd_rec
	record /fil_sky/ coadd_recs(fcc_max_coadd)

	external	fsl_ctifgread
	external	fsl_eof
	external	fsl_normal

	num = 0
	ct_stat(1) = ctp_normal

C
C  Read and filter the data.
C
	do while ((ct_stat(1) .eq. ctp_normal)  .and.
     .	          (num .lt. fcc_max_coadd))

c
c  Read the coadd records.
c
	   call ct_read_arcv (, ct_lun, coadd_rec, ct_stat)

	   if ((ct_stat(1) .eq. ctp_normal)  .and.
     .	       (coadd_rec.coad_spec_data.dq_summary_flag .le.
     .						    fcc_instr_qual)  .and.
     .	       (coadd_rec.coad_spec_data.att_summary_flag .le.
     .						     fcc_attit_qual)) then

c
c  Put the coadds into the coadd buffer.
c
	      num = num + 1
	      coadd_recs(num) = coadd_rec

	   endif

	enddo	!  while (ct_stat


C
C  Check the status of the read.
C

	if (ct_stat(1) .eq. ctp_normal) then
	   status = %loc(fsl_normal)
	elseif (ct_stat(1) .eq. ctp_endoffile) then
	   status = %loc(fsl_eof)
	else
	   status = %loc(fsl_ctifgread)
	   call lib$signal (fsl_ctifgread, %val(2), fcc_infile(1:fcc_inlen),
     .					   %val(ct_stat(1)))
	endif
	          

	fsl_read_cal_coadd = status

	return
	end
