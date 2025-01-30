	integer * 4 function  ffi_read_hybrid_coadd (ct_lun1, ct_lun2, num,
     .					       coadd_recs)

c-------------------------------------------------------------------------------
c
c	Function FFI_READ_HYBRID_COADD
c
c	This function reads in coadded IFGs from the Cobetrieve calibration 
c	archives.  Each time through, the function reads one up to FAC_MAX_COAD
c	coadd records from the cal archive.  The routine decides which input
c	stream coadd to place in the input buffer.  The last time through, the 
c	function closes the Cobetrieve archive.
c
c	Author:	 Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 9 March 1993
c		 SER 10763
c
c-------------------------------------------------------------------------------
c
c	Input:
c		ct_lun1		integer * 4		coadd file lun	
c		ct_lun2		integer * 4		coadd file lun	
c
c	Output:
c		num		integer * 4		number of coadds read
c		coadd_recs	coadd records		input coadd records
c
c	Subroutines called:
c		ct_close_arcv
c		ffi_bracket_coadd
c		fut_free_lun
c		lib$signal
c
c	Include files:
c		ct$library:ctuser.inc
c		ffi_invoc.txt
c		ffi_config.txt
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
	include '(ffi_config)'
	include '(ffi_invoc)'

	integer * 2	ct_stat1(20)	!  CT return status
	integer * 2	ct_stat2(20)	!  CT return status

	integer * 4	ct_lun1		!  CT logical unit number
	integer * 4	ct_lun2		!  CT logical unit number
	integer * 4	num		!  number of coadd records
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	logical * 1	bracket /.true./ !  bracketing coadd flag

	real	* 4	peak_height	!  ifg peak height

	integer * 4	ffi_bracket_coadd
	integer * 4	fut_free_lun

	dictionary 'fic_sky'
	record /fic_sky/ coadd_rec1
	record /fic_sky/ coadd_rec2
	record /fic_sky/ coadd_recs(fac_max_coad)

	external	ffi_ctifgclose
	external	ffi_ctifgread
	external	ffi_eof
	external	ffi_normal
	external	fut_normal

C
C  Read and filter the coadded IFGs.
C
	num = 0
	status = %loc(ffi_normal)

	do while ((status .eq. %loc(ffi_normal))  .and.
     .		  (num .lt. fac_max_coad))

c
c  Read the coadd records.
c
	   status = ffi_bracket_coadd (ct_lun1, ct_lun2, coadd_rec1,
     .				       coadd_rec2, bracket)

	   if ((status .eq. %loc(ffi_normal))  .and.
     .	       (coadd_rec1.coad_spec_data.dq_summary_flag .le.
     .						    fcc_instr_qual)  .and.
     .	       (coadd_rec1.coad_spec_data.att_summary_flag .le.
     .						    fcc_attit_qual)) then

c
c  Put the coadds (from the correct stream) into the coadd buffer.
c
	      num = num + 1
	      if ((bracket .eq. .true.)  .and.
     .	          (coadd_rec2.coad_spec_data.dq_summary_flag .le.
     .					     fcc_instr_qual)  .and.
     .	          (coadd_rec2.coad_spec_data.att_summary_flag .le.
     .					     fcc_attit_qual)) then
	         peak_height = abs(coadd_rec1.coad_data.ifg(peak_pos))
	         if (peak_height .gt. fcc_threshold) then
	            coadd_recs(num) = coadd_rec2
	            fcc_nspec2 = fcc_nspec2 + 1
	         else
	            coadd_recs(num) = coadd_rec1
	            fcc_nspec1 = fcc_nspec1 + 1
	         endif
	      else
	         coadd_recs(num) = coadd_rec1
	         fcc_nspec1 = fcc_nspec1 + 1
	      endif

	   endif

	enddo	!  while (status from read


	if (status .eq. %loc(ffi_eof)) then
C
C  Close the archive after the final read.
C
	   call ct_close_arcv (, ct_lun1, ct_stat1)
	   if (ct_stat1(1) .ne. ctp_normal) then
	      status = %loc(ffi_ctifgclose)
	      call lib$signal (ffi_ctifgclose, %val(2),
     .			       fcc_infile1(1:fcc_inlen1), %val(ct_stat1(1)))
	   endif

	   call ct_close_arcv (, ct_lun2, ct_stat2)
	   if (ct_stat2(1) .ne. ctp_normal) then
	      status = %loc(ffi_ctifgclose)
	      call lib$signal (ffi_ctifgclose, %val(2),
     .			       fcc_infile2(1:fcc_inlen2), %val(ct_stat2(1)))
	   endif

c
c  Free the Cobetrieve logical unit number.
c
	   rstatus = fut_free_lun (ct_lun1)
	   rstatus = fut_free_lun (ct_lun2)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	endif		!  (end of file


	ffi_read_hybrid_coadd = status

	return
	end
