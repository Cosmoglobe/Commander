	integer * 4 function  fsl_write_sky_spec (flag, cs_lun, num, spec_recs)

c-------------------------------------------------------------------------------
c
c	Function FSL_WRITE_SKY_SPEC
c
c	This function writes sky voltage spectrum records and calibrated
c	spectrum records out to the Cobetrieve skymap archives.  
c	The spectra are written out by the CSA routine CSA_Write_Pixels.
c
c	Author:	  
c                FCF_Write_Sky_Spec
c                Gene Eplee
c		 General Sciences Corp.
c		 13 March 1992
c       
c                FSL_Write_Sky_Spec
c                Shirley M. Read
c                Hughes STX Corporation
c                July 1995
c
c-------------------------------------------------------------------------------
c
c	Input:
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
c		csa_write_pixels
c		lib$signal
c
c	Include files:
c		fsl_invoc.txt
c		fut_params.txt
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
c       Modified FCF_Write_Sky_Spec to FSL_Write_Sky_Spec for the new FIRAS 
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

	include '(fut_params)'
	include '(fsl_invoc)'

	integer * 4	flag		!  spectrum type flag
	integer * 4	cs_lun		!  CT logical unit number for sky spec
	integer * 4	num		!  number of spectra in ensemble

	character * 48  out_file        !  Cobetrieve file name

	integer * 2	out_len		!  length of skymap file name

	integer * 4	blocks		!  block count for csa_write_pixels
	integer * 4	cstatus		!  CSA return status
	integer * 4	status		!  return status
	integer * 4     ix              !  index

	integer * 4	csa_write_pixels

	dictionary 'fsl_sky'
	record /fsl_sky/ spec_recs(num)

	external	csa_normal
	external	fsl_csaspecwrite
	external	fsl_normal

C
C  Initialize function return.
C

	status = %loc(fsl_normal)

C
C  Determine the file names for message.
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

	do ix = 1, num

	      spec_recs(ix).coad_spec_head.coadd_no = ix

	end do

	blocks = 0

	cstatus = csa_write_pixels (cs_lun, spec_recs, num, blocks)

	if (cstatus .ne. %loc(csa_normal)) then

	      status = %loc(fsl_csaspecwrite)

	      call lib$signal (fsl_csaspecwrite, %val(2),
     .                         out_file(1:out_len), %val(cstatus))
	endif


	fsl_write_sky_spec = status

	return
	end
