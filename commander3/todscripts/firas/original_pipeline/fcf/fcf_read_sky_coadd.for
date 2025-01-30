	integer * 4 function  fcf_read_sky_coadd (ct_lun, num, coadd_recs)

c-------------------------------------------------------------------------------
c
c	Function FCF_READ_SKY_COADD
c
c	This function reads the sky coadd records from the Cobetrieve archive
c	via CSA and filters them by the data quality summary flags.  The
c	function reads one populated pixel per invocation and returns up to 
c	FCC_MAX_COADD records to the calling routine.
c	
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  21 April 1992
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
c		fcf_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(csa_pixel_input_rec)'
	include '(csa_pixel_output_rec)'
	include '(fut_params)'
	include '(fcf_invoc)'

	integer * 4	blocks		!  block count for csa_read_pixels
	integer * 4	cstatus		!  CSA return status
	integer * 4	ct_lun		!  Cobetrieve logical unit number
	integer * 4	j		!  a counter
	integer * 4	max_in		!  maximun number of coadds to read
					!    in one call to csa
	integer * 4	num		!  number of coadd records read
	integer * 4	num_in		!  number of pixels to be read
					!    in one call to csa
	integer * 4	num_out		!  number of pixels read by csa
	integer * 4	pixel_no	!  pixel number to be read by csa
	integer * 4	status		!  return status

	logical * 1	first_time /.true./	!  status flag

	integer * 4	csa_read_pixels

	dictionary 'fic_sky'
	record /fic_sky/ in_recs(fcc_max_coadd)
	record /fic_sky/ coadd_recs(fcc_max_coadd)

	record /pixel_input_list/  inlist
	record /pixel_output_list/ outlist

	external	csa_normal
	external	fcf_csaifgread
	external	fcf_eof
	external	fcf_csapixnoinval
	external	fcf_maxcoaddexc
	external	fcf_normal

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
	   pixel_no = -1
	endif

	num = 0
	cstatus = %loc(csa_normal)
	status = %loc(fcf_normal)


C
C  Read and filter the data, one pixel at a time.
C

	do while ((cstatus .eq. %loc(csa_normal))  .and.
     .		  (status .eq. %loc(fcf_normal))  .and.
     .		  (pixel_no .lt. 6143)  .and.
     .		  (num .eq. 0))

c
c  Read the pixel.
c
	   pixel_no = pixel_no + 1
	   inlist.pixel_no = pixel_no
	   cstatus = csa_read_pixels (ct_lun, inlist, num_in, in_recs,
     .				      max_in, outlist, num_out, blocks)

	   if (outlist.no_records .gt. fcc_max_coadd) then
	      status = %loc(fcf_maxcoaddexc)
	      call lib$signal (fcf_maxcoaddexc, %val(2),
     .			       %val(outlist.no_records), %val(pixel_no))
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

	if (status .eq. %loc(fcf_normal)) then
	   if (cstatus .eq. %loc(csa_normal)) then
	      if (pixel_no .lt. 6143) then
	         status = %loc(fcf_normal)
	      elseif (pixel_no .eq. 6143) then
	         status = %loc(fcf_eof)
	      else
	         status = %loc(fcf_csapixnoinval)
	         call lib$signal (fcf_csapixnoinval, %val(2), %val(pixel_no),
     .						     fcc_infile(1:fcc_inlen))
	      endif
	   else
	      status = %loc(fcf_csaifgread)
	      call lib$signal (fcf_csaifgread, %val(2), fcc_infile(1:fcc_inlen),
     .					       %val(cstatus))
	   endif
	endif


	fcf_read_sky_coadd = status

	return
	end
