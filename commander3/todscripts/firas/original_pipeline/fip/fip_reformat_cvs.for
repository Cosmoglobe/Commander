	integer * 4 function  fip_reformat_cvs ()

c-------------------------------------------------------------------------------
c
c	Function FIP_REFORMAT_CVS
c
c	This function reads the FIRAS C-Vector from the binary file
c	CSDR$FIRAS_IN:FMS_CVS_CCSS.VVV_XXXXXXXXXX.  It performs the frequency
c	cut, renormalizes the C-Vector to MJy/sr, then writes the C-Vector to
c	the FIRAS project dataset file
c	CSDR$FIRAS_OUT:FIP_CVS_CCSS.VVV_XXXXXXXXXX.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  17 November 1994
c
c-------------------------------------------------------------------------------
c
c	Subroutines called:
c		fut_free_lun
c		fut_get_lun
c		lib$movc5
c		lib$signal
c		str$trim
c
c	Include files:
c		fip_frequency.txt
c		fip_invoc_cvs.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fip_invoc_cvs)'
	include '(fip_frequency)'

	integer * 4	fip_lun			!  fip file logical unit number
	integer * 4	fms_lun			!  fms file logical unit number
	integer * 4	io_stat			!  I/O return status
	integer * 4	j			!  a counter
	integer * 4	rstatus			!  return status
	integer * 4	status			!  return status

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fms_cvs'
	record /fms_cvs/ fms_cvs
	dictionary 'fip_cvs'
	record /fip_cvs/ fip_cvs

	external	fip_normal
	external	fip_rmsclose
	external	fip_rmsopen
	external	fip_rmsread
	external	fip_rmswrite
	external	fut_normal

c
c  Initialize the routine.
c
	status = %loc(fip_normal)


C
C  Read the C-Vector from the FMS file.
C

c
c  Find the C-Vector file name.
c
	fcc_fms_cfile = 'CSDR$FIRAS_IN:FMS_CVS_' // fcc_scan_mode // '.' //
     .			fcc_file_ext
	call str$trim (fcc_fms_cfile, fcc_fms_cfile, fcc_fmslen)

c
c  Open the C-Vector file.
c
	rstatus = fut_get_lun(fms_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=fms_lun, file=fcc_fms_cfile, status='old',
     .	      form='unformatted', access='sequential', readonly, iostat=io_stat)

	if (io_stat .eq. 0) then
c
c  Read the C-Vector file.
c
	   read (fms_lun, iostat=io_stat) fms_cvs
	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmsread)
	      call lib$signal (fip_rmsread, %val(2),
     .			       fcc_fms_cfile(1:fcc_fmslen), %val(io_stat))
	   endif

c
c  Close the C-Vector file.
c
	   close (unit=fms_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmsclose)
	      call lib$signal (fip_rmsclose, %val(2),
     .			       fcc_fms_cfile(1:fcc_fmslen), %val(io_stat))
	   endif

	   rstatus = fut_free_lun(fms_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	else
	   status = %loc(fip_rmsopen)
	   call lib$signal (fip_rmsopen, %val(2), fcc_fms_cfile(1:fcc_fmslen),
     .					 %val(io_stat))
	endif	!  (open status


	if (status .eq. %loc(fip_normal)) then
C
C  Reformat the C-Vector.
C

c
c  Fill in the C-Vector identification fields.
c
	   fip_cvs.model_label = fms_cvs.model_label
	   fip_cvs.chanscan    = fcc_scan_mode
	   fip_cvs.galatexc    = fms_cvs.galat_exc
	   fip_cvs.nu_zero     = fcc_nu0
	   fip_cvs.delta_nu    = fcc_dnu
	   fip_cvs.num_freq    = fcc_nfreq

c
c  Write the C-Vector to the FIP RDL, performing the frequency cut
c	and the renormalization.
c
	   do j = fcc_jlo,fcc_jhi
	      fip_cvs.c_vector(j-fcc_freq_offset) =
     .				sqrt(fms_cvs.cvector(j)) * fac_erg_to_mjy
	   enddo


C
C  Write the FIP RDL to the project dataset file.
C

c
c  Get the project dataset file name.
c
	   fcc_fip_cfile = 'CSDR$FIRAS_OUT:FIP_CVS_' // fcc_scan_mode // '.' //
     .			    fcc_file_ext
	   call str$trim (fcc_fip_cfile, fcc_fip_cfile, fcc_fiplen)

c
c  Open the project dataset file.
c
	   rstatus = fut_get_lun(fip_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	   open (unit=fip_lun, file=fcc_fip_cfile, status='new',
     .		 form='unformatted', recordtype='fixed', recl=fcc_fip_csize,
     .		 access='sequential', iostat=io_stat)

	   if (io_stat .eq. 0) then
c
c  Write the C-Vector to the project dataset file.
c
	      write (fip_lun, iostat=io_stat) fip_cvs
	      if (io_stat .ne. 0) then
	         status = %loc(fip_rmswrite)
	         call lib$signal (fip_rmswrite, %val(2),
     .				  fcc_fip_cfile(1:fcc_fiplen), %val(io_stat))
	      endif

c
c  Close the project dataset file.
c
	      close (unit=fip_lun, iostat=io_stat)
	      if (io_stat .ne. 0) then
	         status = %loc(fip_rmsclose)
	         call lib$signal (fip_rmsclose, %val(2),
     .				  fcc_fip_cfile(1:fcc_fiplen), %val(io_stat))
	      endif

	      rstatus = fut_free_lun(fip_lun)
	      if (rstatus .ne. %loc(fut_normal)) then
	         call lib$signal (%val(rstatus))
	      endif

	   else
	      status = %loc(fip_rmsopen)
	      call lib$signal (fip_rmsopen, %val(2),
     .			       fcc_fip_cfile(1:fcc_fiplen), %val(io_stat))
	   endif	!  (open status

	endif		!  (status from FMS file I/O


	fip_reformat_cvs = status

	return
	end
