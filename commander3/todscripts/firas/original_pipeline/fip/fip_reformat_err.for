	integer * 4 function  fip_reformat_err ()

c-------------------------------------------------------------------------------
c
c	Function FIP_REFORMAT_ERR
c
c	This function reads the FIRAS error terms from the binary file
c	CSDR$FIRAS_IN:FER_ERR_CCSS.VVV_XXXXXXXXXX.  It performs the frequency
c	cut, renormalizes the spectra to MJy/sr, then writes the spectra to the
c	FIRAS project dataset file CSDR$FIRAS_OUT:FIP_ERR_CCSS.VVV_XXXXXXXXXX.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  11 August 1994
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
c		fip_invoc_err.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fip_invoc_err)'
	include '(fip_frequency)'

	integer * 4	fer_lun			!  fer file logical unit number
	integer * 4	fip_lun			!  fip file logical unit number
	integer * 4	io_stat			!  I/O return status
	integer * 4	j			!  a counter
	integer * 4	k			!  a counter
	integer * 4	rstatus			!  return status
	integer * 4	status			!  return status

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fer_err'
	record /fer_err/ fer_err
	dictionary 'fip_err'
	record /fip_err/ fip_err

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
C  Read the errors from the FER file.
C

c
c  Find the error file name.
c
	fcc_fer_file = 'CSDR$FIRAS_IN:FER_ERR_' // fcc_scan_mode // '.' //
     .			fcc_file_ext
	call str$trim (fcc_fer_file, fcc_fer_file, fcc_ferlen)

c
c  Open the error file.
c
	rstatus = fut_get_lun(fer_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=fer_lun, file=fcc_fer_file, status='old',
     .	      form='unformatted', access='sequential', readonly, iostat=io_stat)

	if (io_stat .eq. 0) then
c
c  Read the error file.
c
	   read (fer_lun, iostat=io_stat) fer_err
	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmsread)
	      call lib$signal (fip_rmsread, %val(2),
     .			       fcc_fer_file(1:fcc_ferlen), %val(io_stat))
	   endif

c
c  Close the error file.
c
	   close (unit=fer_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmsclose)
	      call lib$signal (fip_rmsclose, %val(2),
     .			       fcc_fer_file(1:fcc_ferlen), %val(io_stat))
	   endif

	   rstatus = fut_free_lun(fer_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	else
	   status = %loc(fip_rmsopen)
	   call lib$signal (fip_rmsopen, %val(2), fcc_fer_file(1:fcc_ferlen),
     .					 %val(io_stat))
	endif	!  (open status


	if (status .eq. %loc(fip_normal)) then
C
C  Reformat the errors.
C

c
c  Fill in the error identification fields.
c
	   fip_err.model_label = fer_err.model_label
	   fip_err.chanscan    = fcc_scan_mode
	   fip_err.nu_zero     = fcc_nu0
	   fip_err.delta_nu    = fcc_dnu
	   fip_err.num_freq    = fcc_nfreq
	   fip_err.galatexc    = fer_err.galatexc

c
c  Write the error terms to the FIP RDL, performing the frequency cut
c	and the renormalization.
c
	   do j = fcc_jlo,fcc_jhi
	      fip_err.d_vector(j-fcc_freq_offset)   = fer_err.d_vector(j) *
     .							 fac_erg_to_mjy
	      fip_err.c_vector(j-fcc_freq_offset)   = fer_err.c_vector(j) *
     .							 fac_erg_to_mjy
	      fip_err.pep_gain(j-fcc_freq_offset)   = fer_err.pep_gain(j)
	      fip_err.pep_offset(j-fcc_freq_offset) = fer_err.pep_offset(j) *
     .							 fac_erg_to_mjy
	      do k=1,14
	         fip_err.jcj_gain(k,j-fcc_freq_offset)   = fer_err.jcj_gain(k,j)
	         fip_err.jcj_offset(k,j-fcc_freq_offset) =
     .					fer_err.jcj_offset(k,j) * fac_erg_to_mjy
	      enddo
	      fip_err.pup_temp			    = fer_err.pup_temp
	      fip_err.pup_spec(j-fcc_freq_offset)   = fer_err.pup_spec(j) *
     .							 fac_erg_to_mjy
	      fip_err.ptp_temp			    = fer_err.ptp_temp
	      fip_err.ptp_spec(j-fcc_freq_offset)   = fer_err.ptp_spec(j) *
     .							 fac_erg_to_mjy
	   enddo


C
C  Write the FIP RDL to the project dataset file.
C

c
c  Get the project dataset file name.
c
	   fcc_fip_file = 'CSDR$FIRAS_OUT:FIP_ERR_' // fcc_scan_mode // '.' //
     .			    fcc_file_ext
	   call str$trim (fcc_fip_file, fcc_fip_file, fcc_fiplen)

c
c  Open the project dataset file.
c
	   rstatus = fut_get_lun(fip_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	   open (unit=fip_lun, file=fcc_fip_file, status='new',
     .		 form='unformatted', recordtype='fixed', recl=fcc_fip_size,
     .		 access='sequential', iostat=io_stat)

	   if (io_stat .eq. 0) then
c
c  Write the correction spectra to the project dataset file.
c
	      write (fip_lun, iostat=io_stat) fip_err
	      if (io_stat .ne. 0) then
	         status = %loc(fip_rmswrite)
	         call lib$signal (fip_rmswrite, %val(2),
     .				  fcc_fip_file(1:fcc_fiplen), %val(io_stat))
	      endif

c
c  Close the project dataset file.
c
	      close (unit=fip_lun, iostat=io_stat)
	      if (io_stat .ne. 0) then
	         status = %loc(fip_rmsclose)
	         call lib$signal (fip_rmsclose, %val(2),
     .				  fcc_fip_file(1:fcc_fiplen), %val(io_stat))
	      endif

	      rstatus = fut_free_lun(fip_lun)
	      if (rstatus .ne. %loc(fut_normal)) then
	         call lib$signal (%val(rstatus))
	      endif

	   else
	      status = %loc(fip_rmsopen)
	      call lib$signal (fip_rmsopen, %val(2),
     .			       fcc_fip_file(1:fcc_fiplen), %val(io_stat))
	   endif	!  (open status

	endif		!  (status from FER file I/O


	fip_reformat_err = status

	return
	end
