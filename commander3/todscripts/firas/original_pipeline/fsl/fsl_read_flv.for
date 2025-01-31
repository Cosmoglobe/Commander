	integer * 4 function  fsl_read_flv ()

c-------------------------------------------------------------------------------
c
c	Function FSL_READ_FLV
c
c	This function opens, reads, and closes the FIRAS all-sky, averaged FIL
c       variance vectors reference dataset FEX_FLV_CCSS.DAT from the directory 
c       CSDR$FIRAS_UREF (CC is channel, SS is scan mode), The FIL variances are
c       written into the include file FSL_MODEL for use by FSL.
c
c	Author:   
c                FSL_Read_FLV
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		none
c
c	Subroutines called:
c		fut_free_lun
c		fut_get_lun
c		lib$movc5
c		lib$signal
c		str$trim
c
c	Include files:
c		fsl_invoc.txt
c		fsl_model.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c       Changes for FSL:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fsl_invoc)'
	include '(fsl_model)'

	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status
	integer * 4	flv_lun		!  FIL variances file lun

	logical * 1	flv_ver		!  FIL variances verification flag

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fex_flv'
	record /fex_flv/ fex_flv

	external	fsl_flvclose
	external	fsl_flvinval
	external	fsl_flvopen
	external	fsl_flvread
	external	fsl_normal
	external	fut_normal

	status = %loc(fsl_normal)

C
C  Open, read, and close the FIL variances file.
C

c
c  Find the FIL variances filename.
c
	fcc_flv_file = 'CSDR$FIRAS_UREF:FEX_FLV_' // fcc_scan_mode // '.DAT' 
	call str$trim (fcc_flv_file, fcc_flv_file, fcc_flvlen)

c
c  Open the FIL variances file.
c
	rstatus = fut_get_lun(flv_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=flv_lun, file=fcc_flv_file, status='old',
     .	      form='unformatted', recordtype='fixed', recl=fcc_flv_size,
     .	      access='sequential', readonly, iostat=io_stat)

	if (io_stat .eq. 0) then
c
c  Read the FIL variances.
c
	   read (flv_lun, iostat=io_stat) fex_flv
	   if (io_stat .ne. 0) then
	      status = %loc(fsl_flvread)
	      call lib$signal (fsl_flvread, %val(2),
     .			       fcc_flv_file(1:fcc_flvlen), %val(io_stat))
	   endif

c
c  Close the FIL variances file.
c
	   close (unit=flv_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fsl_flvclose)
	      call lib$signal (fsl_flvclose, %val(2),
     .			       fcc_flv_file(1:fcc_flvlen), %val(io_stat))
	   endif

	   rstatus = fut_free_lun(flv_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	else
	   status = %loc(fsl_flvopen)
	   call lib$signal (fsl_flvopen, %val(2),
     .			    fcc_flv_file(1:fcc_flvlen), %val(io_stat))
	endif	!  (open status


	if (status .eq. %loc(fsl_normal)) then
C
C  Process the FIL variances.
C

c
c  Verify that the correct FIL variances has been read.
c
	   flv_ver = .true.
	   if (fcc_chan .ne. fex_flv.channel) flv_ver = .false.
	   if (fcc_smode .ne. fex_flv.scan_mode) flv_ver = .false.
	   if (flv_ver .eq. .false.) then
	      status = %loc(fsl_flvinval)
	      call lib$signal (fsl_flvinval, %val(1),
     .			       fcc_flv_file(1:fcc_flvlen))
	   endif

	endif

	if (status .eq. %loc(fsl_normal)) then
c
c  Insert the FIL variances into the common block.
c
	   call lib$movc5 (0,,0,8664,fil_var)
	   do j = 1,361
	      fil_var(1,j) = fex_flv.rr_variances(j)
	      fil_var(2,j) = fex_flv.ii_variances(j)
	      fil_var(3,j) = fex_flv.ri_variances(j)
	   enddo 
	endif


	fsl_read_flv = status

	return
	end
