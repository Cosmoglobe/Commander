	integer * 4 function  fsl_read_sky_coadd_list (ct_lun, num, coadd_recs)

c-------------------------------------------------------------------------------
c
c	Function FSL_READ_SKY_COADD_LIST
c
c	This function reads the sky coadd records from the Cobetrieve archive
c	for the fcc_num_pix pixels specified in the array fcc_plist via CSA and
c	filters them by the data quality summary flags.  The function reads one
c	populated pixel per invocation and returns up to FCC_MAX_COADD records
c	to the calling routine.
c
c	Author:   
c 		 FCF_Read_Sky_Coadd_List
c                Gene Eplee
c		 General Sciences Corp.
c		 21 April 1992
c       
c                FSL_Read_Sky_Coadd_List
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
c		csa_read_pixels
c		lib$signal
c
c	Include files:
c		csa_pixel_input_rec.txt
c		csa_pixel_output_rec.txt
c		fsl_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, August 7, 1995 
c       Modified FCF_Read_Sky_Coadd_List to FSL_Read_Sky_Coadd_List for the 
c       new FIRAS pipeline which will process long spectra to get improved 
c       frequency resolution.
c           1. Changed status, include file, and function names for FSL.
c           2. Changed dictionary and record from FIC to FIL.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(csa_pixel_input_rec)'
	include '(csa_pixel_output_rec)'
	include '(fut_params)'
	include '(fsl_invoc)'

	integer * 4	blocks		!  block count for csa_read_pixels
	integer * 4	cstatus		!  CSA return status
	integer * 4	ct_lun		!  Cobetrieve logical unit number
	integer * 4	j		!  a counter
	integer * 4	max_in		!  maximun number of coadds to read
					!    in one call to csa
	integer * 4	num_pix		!  number of pixel in fcc_plist
					!    to be read by csa
	integer * 4	num		!  number of coadd records read
	integer * 4	num_in		!  number of pixels to be read
					!    in one call to csa
	integer * 4	num_out		!  number of pixels read by csa
	integer * 4	status		!  return status

	logical * 1	first_time /.true./	!  status flag

	integer * 4	csa_read_pixels

	dictionary 'fil_sky'
	record /fil_sky/ in_recs(fcc_max_coadd)
	record /fil_sky/ coadd_recs(fcc_max_coadd)

	record /pixel_input_list/  inlist
	record /pixel_output_list/ outlist

	external	csa_normal
	external	fsl_csaifgread
	external	fsl_eof
	external	fsl_csapixnoinval
	external	fsl_maxcoaddexc
	external	fsl_normal

C
C  Initialize the read.
C

	if (first_time) then
	   first_time = .false.
	   blocks = 0
	   fcc_npix = 0
	   inlist.level_no = fac_skymap_level
	   max_in = fcc_max_coadd
	   num_in = 1
	   num_pix = 0
	endif

	num = 0
	cstatus = %loc(csa_normal)
	status = %loc(fsl_normal)


C
C  Read and filter the data, one pixel at a time.
C

	do while ((cstatus .eq. %loc(csa_normal))  .and.
     .		  (status .eq. %loc(fsl_normal))  .and.
     .		  (num_pix .lt. fcc_num_pix)  .and.
     .		  (num .eq. 0))

c
c  Read the specified pixels.
c
	   num_pix = num_pix + 1
	   inlist.pixel_no = fcc_plist(num_pix)
	   cstatus = csa_read_pixels (ct_lun, inlist, num_in, in_recs,
     .				      max_in, outlist, num_out, blocks)

	   if (outlist.no_records .gt. fcc_max_coadd) then
	      status = %loc(fsl_maxcoaddexc)
	      call lib$signal (fsl_maxcoaddexc, %val(2),
     .			     %val(outlist.no_records), %val(fcc_plist(num_pix)))
	   elseif ((outlist.no_records .ne. 0)  .and.
     .		   (cstatus .eq. %loc(csa_normal))) then
	      fcc_npix = fcc_npix + 1
c
c  Filter the coadds.
c
	      do j = 1, outlist.no_records
	         if ((in_recs(j).coad_spec_data.dq_summary_flag .le.
     .						       fcc_instr_qual)  .and.
     .		     (in_recs(j).coad_spec_data.att_summary_flag .le.
     .						       fcc_attit_qual)) then
c
c  Put the coadds into the coadd buffer.
c
	            num = num + 1
	            coadd_recs(num) = in_recs(j)
	         endif
	      enddo
	   endif	!  (outlist

	enddo		!  while (cstatus


C
C  Check the status of the read.
C

	if (status .eq. %loc(fsl_normal)) then
	   if (cstatus .eq. %loc(csa_normal)) then
	      if (num_pix .lt. fcc_num_pix) then
	         status = %loc(fsl_normal)
	      elseif (num_pix .eq. fcc_num_pix) then
	         status = %loc(fsl_eof)
	      else
	         status = %loc(fsl_csapixnoinval)
	         call lib$signal (fsl_csapixnoinval, %val(2),
     .				  %val(fcc_plist(num_pix)),
     .				  fcc_infile(1:fcc_inlen))
	      endif
	   else
	      status = %loc(fsl_csaifgread)
	      call lib$signal (fsl_csaifgread, %val(2), fcc_infile(1:fcc_inlen),
     .					       %val(cstatus))
	   endif
	endif


	fsl_read_sky_coadd_list = status

	return
	end
