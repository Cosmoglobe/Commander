	program frd_variances

c-------------------------------------------------------------------------------
c
c	Program FRD_VARIANCES
c
c	This program produces all-sky average variances for calibrated FIRAS
c	spectra.  The program reads the spectrum variances of the skymaps
c	specified in an input script file.  These variances are averaged to
c	produce an all-sky average variance.  This average variance is written
c	to the reference dataset FEX_VAR_CCSS.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  27 October 1993
c		  SER 11409
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		ct_init
c		cut_display_banner
c		cut_register_version
c		frd_parse_variances
c		frd_read_variances
c		fut_free_lun
c		fut_get_lun
c		lib$signal
c
c	Include files:
c		csdr$library:ctuser.inc
c		frd_invoc_variances.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include 'csdr$library:ctuser.inc'
	include '(fut_params)'
	include '(frd_invoc_variances)'

	character * 14	current_gmt		!  GMT time of invocation
	character * 40  label			!  cal model solution label

	integer * 2	ct_stat(20)		!  CT return status

	integer * 4	current_time(2)		!  ADT time of invocation
	integer * 4	io_stat			!  I/O return status
	integer * 4	j			!  a counter
	integer * 4 	lun_out /6/		!  terminal lun
	integer * 4	nifgs			!  number of ifgs in spectra
	integer * 4	nspec			!  number of spectra read
	integer * 4	rstatus			!  return status
	integer * 4	out_lun			!  output file lun
	integer * 4	status			!  return status

	real	* 8	sum_nk			!  sum of (n-k)
	real	* 8	sum_var(257)		!  sum of variances

	integer * 4	cut_display_banner
	integer * 4	cut_register_version
	integer * 4	frd_parse_variances
	integer * 4	frd_read_variances
	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fex_var'
	record /fex_var/ out_rec

	external	frd_ctinit
	external	frd_failure
	external	frd_normal
	external	frd_rmsclose
	external	frd_rmsopen
	external	frd_rmswrite
	external	fut_normal

c
c  Parse the command line.
c
	status = frd_parse_variances (current_time, current_gmt)

c
c  Print the banner.
c
	rstatus = cut_register_version (fcc_version)
	rstatus = cut_display_banner (lun_out, 80,
     .				     'FIRAS Facility FRD_VARIANCES')
	write(lun_out,10)
 10	format (/)

	if (status .eq. %loc(frd_normal)) then
c
c  Initialize Cobetrieve.
c
	   call ct_init(ct_stat)


	   if (ct_stat(1) .eq. ctp_normal) then
c
c  Read the variances from the spectrum records.
c
	      status = frd_read_variances (label, nifgs, nspec, sum_var, sum_nk)

	      if (status .eq. %loc(frd_normal)) then
c
c  Compute the average variances.
c
	         do j = 1,257
	            out_rec.variances(j) = sum_var(j) / sum_nk
	         enddo
	         out_rec.nsky_ifgs = floatj(nifgs)

c
c  Put the scan mode information into the RDL.
c
	         out_rec.gmt         = current_gmt
	         out_rec.time(1)     = current_time(1)
	         out_rec.time(2)     = current_time(2)
	         out_rec.channel     = fcc_chan
	         out_rec.scan_length = fcc_length
	         out_rec.scan_speed  = fcc_speed
	         out_rec.model_label = label
	         out_rec.galat_exc   = fcc_glat

c
c  Write out the average variances to the binary file.
c
	         rstatus = fut_get_lun (out_lun)
	         if (rstatus .ne. %loc(fut_normal)) then
	            call lib$signal (%val(rstatus))
	         endif

	         open (unit=out_lun, file=fcc_outfile, status='new',
     .		       form='unformatted', recordtype='fixed',
     .		       recl=fcc_rec_size, access='sequential', iostat=io_stat)

	         if (io_stat .eq. 0) then
	            write (out_lun, iostat=io_stat) out_rec
	            if (io_stat .ne. 0) then
	               status = %loc(frd_rmswrite)
	               call lib$signal(frd_rmswrite, %val(2),
     .				       fcc_outfile(1:fcc_outlen), %val(io_stat))
	            endif

c
c  Close the binary file.
c
	            close (out_lun, iostat = io_stat)
	            if (io_stat .ne. 0) then
	               status = %loc(frd_rmsclose)
	               call lib$signal (frd_rmsclose, %val(2),
     .				       fcc_outfile(1:fcc_outlen), %val(io_stat))
	            endif

	         else
	            status = %loc(frd_rmsopen)
	            call lib$signal (frd_rmsopen, %val(2),
     .				     fcc_outfile(1:fcc_outlen), %val(io_stat))
	         endif	!  status from RMS open

c
c  Free the output file logical unit number.
c
	         rstatus = fut_free_lun (out_lun)
	         if (rstatus .ne. %loc(fut_normal)) then
	            call lib$signal (%val(rstatus))
	         endif

	      endif	!  status from read variances

	   else
	      status = %loc(frd_ctinit)
	      call lib$signal (frd_ctinit, %val(1), %val(ct_stat(1)))
	   endif	!  ct_stat from ct_init

	endif		!  status from parse

c
c  Write out number of ifgs, spectra, and skymaps processed.
c
	if (status .eq. %loc(frd_normal)) then
	   write (lun_out, 20, iostat=io_stat) fcc_mapnum, nspec, nifgs
	endif
  20	format (/, 4x, 'Number of skymaps processed:     ', I6,
     .		/, 4x, 'Number of sky spectra read:      ', I6,
     .		/, 4x, 'Number of sky ifgs represented:  ', I6, //)

c
c  Signal the program completion status.
c
	if (status .eq. %loc(frd_normal)) then
	   call lib$signal (frd_normal)
	else
	   call lib$signal (frd_failure)
	endif


	end
