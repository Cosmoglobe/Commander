	integer * 4 function  ffi_bracket_coadd (ct_lun1, ct_lun2, coadd_rec1,
     .						 coadd_rec2, bracket)

c-------------------------------------------------------------------------------
c
c	Function FFI_BRACKET_COADD
c
c	This function reads in coadded IFGs from the Cobetrieve calibration 
c	archives.  The function matches coadds from the two input streams, if
c	possible.
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
c		coadd_rec1	coadd record		stream 1 input coadd
c		coadd_rec2	coadd record		stream 2 input coadd
c		bracket		logical * 1		bracketing coadd flag
c
c	Subroutines called:
c		ct_read_arcv
c		lib$signal
c		time_gt
c		time_lt
c
c	Include files:
c		ct$library:ctuser.inc
c		ffi_invoc.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include 'ct$library:ctuser.inc'
	include '(ffi_invoc)'

	integer * 2	ct_stat1(20)	!  CT return status
	integer * 2	ct_stat2(20)	!  CT return status

	integer * 4	ct_lun1		!  CT logical unit number
	integer * 4	ct_lun2		!  CT logical unit number
	integer * 4	j		!  a counter
	integer * 4	status		!  return status
	integer * 4	time1a(2)	!  first ifg time of coadd 1
	integer * 4	time1b(2)	!  last ifg time of coadd 1
	integer * 4	time2(2)	!  coadd time of coadd 2

	logical * 1	bracket
	logical * 1	stream /.true./

	logical * 1	time_gt
	logical * 1	time_lt

	dictionary 'fic_sky'
	record /fic_sky/ coadd_rec1
	record /fic_sky/ coadd_rec2
	record /fic_sky/ buffer

	external	ffi_ctifgread
	external	ffi_eof
	external	ffi_normal

C
C  Read the coadd records.
C

	ct_stat1(1) = ctp_normal
	ct_stat2(1) = ctp_normal
	status = %loc(ffi_normal)
c
c  Read the records from the two input streams
c
	call ct_read_arcv (, ct_lun1, coadd_rec1, ct_stat1)
	if (stream .eq. .true.) then
	   if (bracket .eq. .false.) then
	      bracket = .true.
	      coadd_rec2 = buffer
	   else
	      call ct_read_arcv (, ct_lun2, coadd_rec2, ct_stat2)
	   endif
	
	   if ((ct_stat1(1) .eq. ctp_normal)  .and.
     .	       (ct_stat2(1) .eq. ctp_normal)) then
c
c  Check timetags for bracketing records.
c
	      do j=1,2
	         time1a(j) = coadd_rec1.coad_spec_head.first_time(j)
	         time1b(j) = coadd_rec1.coad_spec_head.last_time(j)
	         time2(j)  = coadd_rec2.ct_head.time(j)
	      enddo

	      if (time_lt(time2,time1a)) then
	         do while ((time_lt(time2,time1a))  .and.
     .			   (ct_stat2(1) .eq. ctp_normal))
	            call ct_read_arcv (, ct_lun2, coadd_rec2, ct_stat2)
	            do j=1,2
	               time2(j)  = coadd_rec2.ct_head.time(j)
	            enddo
	         enddo
	      endif

	      if (ct_stat2(1) .eq. ctp_normal) then
	         if (time_gt(time2,time1b)) then
	            bracket = .false.
	            buffer = coadd_rec2
	         else
	            bracket = .true.
	         endif
	      endif
	   endif	!  (ct_stat from initial read

	endif		!  (stream


C
C  Check the status of the read.
C

	if (ct_stat1(1) .eq. ctp_normal) then
	   status = %loc(ffi_normal)
	elseif (ct_stat1(1) .eq. ctp_endoffile) then
	   status = %loc(ffi_eof)
	else
	   status = %loc(ffi_ctifgread)
	   call lib$signal (ffi_ctifgread, %val(2), fcc_infile1(1:fcc_inlen1),
     .					   %val(ct_stat1(1)))
	endif

	if (stream .eq. .true.) then
	   if (ct_stat2(1) .eq. ctp_normal) then
	      stream = .true.
	   elseif (ct_stat2(1) .eq. ctp_endoffile) then
	      stream = .false.
	      bracket = .false.
	   else
	      status = %loc(ffi_ctifgread)
	      call lib$signal (ffi_ctifgread, %val(2),
     .			       fcc_infile2(1:fcc_inlen2), %val(ct_stat2(1)))
	   endif
	endif


	ffi_bracket_coadd = status

	return
	end
