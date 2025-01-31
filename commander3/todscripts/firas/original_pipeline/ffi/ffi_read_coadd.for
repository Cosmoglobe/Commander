	integer * 4 function  ffi_read_coadd (ct_lun1, num, coadd_recs)

c-------------------------------------------------------------------------------
c
c	Function FFI_READ_COADD
c
c	This function reads in coadded IFGs from the Cobetrieve calibration 
c	archive.  Each time through, the function reads one up to FCC_MAX_COAD
c	coadd records from the cal archive.  The last time through, the function
c	closes the Cobetrieve archive.
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
c
c	Output:
c		num		integer * 4		number of coadds read
c		coadd_recs	coadd records		input coadd records
c
c	Subroutines called:
c		ct_close_arcv
c		ct_read_arcv
c		fut_free_lun
c		lib$signal
c
c	Include files:
c		ct$library:ctuser.inc
c		ffi_invoc.txt
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
	include '(ffi_invoc)'

	integer * 2	ct_stat(20)	!  CT return status

	integer * 4	cstatus		!  return status
	integer * 4	ct_lun1		!  CT logical unit number
	integer * 4	num		!  number of coadd records
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status

	integer * 4	fut_free_lun

	dictionary 'fic_sky'
	record /fic_sky/ coadd_rec
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
	ct_stat(1) = ctp_normal

	do while ((ct_stat(1) .eq. ctp_normal)  .and.
     .	          (num .lt. fac_max_coad))

c
c  Read the coadd records.
c
	   call ct_read_arcv (, ct_lun1, coadd_rec, ct_stat)

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
	   status = %loc(ffi_normal)
	elseif (ct_stat(1) .eq. ctp_endoffile) then
	   status = %loc(ffi_eof)
	else
	   status = %loc(ffi_ctifgread)
	   call lib$signal (ffi_ctifgread, %val(2), fcc_infile1(1:fcc_inlen1),
     .					   %val(ct_stat(1)))
	endif
	          

	if (status .eq. %loc(ffi_eof)) then
C
C  Close the archive after the final read.
C
	   call ct_close_arcv (, ct_lun1, ct_stat)
	   if (ct_stat(1) .ne. ctp_normal) then
	      status = %loc(ffi_ctifgclose)
	      call lib$signal (ffi_ctifgclose, %val(2),
     .			       fcc_infile1(1:fcc_inlen1), %val(cstatus))
	   endif

c
c  Free the Cobetrieve logical unit number.
c
	   rstatus = fut_free_lun (ct_lun1)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	endif		!  (end of file


	ffi_read_coadd = status

	return
	end
