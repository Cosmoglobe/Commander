	integer * 4 function  fip_reformat_ejv ()

c-------------------------------------------------------------------------------
c
c	Function FIP_REFORMAT_EJV
c
c	This function reads the FIRAS destriper correction spectra from the
c	binary file CSDR$FIRAS_IN:FEX_EJV_CCSS.VVV_XXXXXXXXXX.  It performs the
c	frequency cut, renormalizes the spectra to MJy/sr, then writes the
c	spectra to the FIRAS project dataset file
c	CSDR$FIRAS_OUT:FIP_EJV_CCSS.VVV_XXXXXXXXXX.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  8 July 1994
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
c		fip_invoc_ejv.txt
c		fut_params.txt
c		uoe_constants.txt
c
c-------------------------------------------------------------------------------
c
c	Changes
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fip_invoc_ejv)'
	include '(fip_frequency)'
	include '(uoe_constants)'

	integer * 4	atime(2)		!  ADT time
	integer * 4	btime(2)		!  ADT time
	integer * 4	fex_lun			!  fex file logical unit number
	integer * 4	fip_lun			!  fip file logical unit number
	integer * 4	io_stat			!  I/O return status
	integer * 4	j			!  a counter
	integer * 4	k			!  a counter
	integer * 4	rstatus			!  return status
	integer * 4	status			!  return status

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	real	* 8	uoe_adt2t68

	dictionary 'fex_ejv'
	record /fex_ejv/ fex_ejv
	dictionary 'fip_ejv'
	record /fip_ejv/ fip_ejv

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
C  Read the offset correction spectra from the FEX file.
C

c
c  Find the correction spectra file name.
c
	fcc_fex_ofile = 'CSDR$FIRAS_IN:FEX_EJV_' // fcc_scan_mode // '.' //
     .			fcc_file_ext
	call str$trim (fcc_fex_ofile, fcc_fex_ofile, fcc_fexlen)

c
c  Open the correction spectra file.
c
	rstatus = fut_get_lun(fex_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=fex_lun, file=fcc_fex_ofile, status='old',
     .	      form='unformatted', access='sequential', readonly, iostat=io_stat)

	if (io_stat .eq. 0) then
c
c  Read the correction spectra file.
c
	   read (fex_lun, iostat=io_stat) fex_ejv
	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmsread)
	      call lib$signal (fip_rmsread, %val(2),
     .			       fcc_fex_ofile(1:fcc_fexlen), %val(io_stat))
	   endif

c
c  Close the correction spectra file.
c
	   close (unit=fex_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmsclose)
	      call lib$signal (fip_rmsclose, %val(2),
     .			       fcc_fex_ofile(1:fcc_fexlen), %val(io_stat))
	   endif

	   rstatus = fut_free_lun(fex_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	else
	   status = %loc(fip_rmsopen)
	   call lib$signal (fip_rmsopen, %val(2), fcc_fex_ofile(1:fcc_fexlen),
     .					 %val(io_stat))
	endif	!  (open status


	if (status .eq. %loc(fip_normal)) then
C
C  Reformat the correction spectra.
C

c
c  Fill in the correcton spectra identification fields.
c
	   fip_ejv.model_label = fex_ejv.label
	   fip_ejv.chanscan    = fcc_scan_mode
	   fip_ejv.nu_zero     = fcc_nu0
	   fip_ejv.delta_nu    = fcc_dnu
	   fip_ejv.num_freq    = fcc_nfreq
	   fip_ejv.galatexc    = fex_ejv.gallat

c
c  Fill in the mission period time fields.
c
	   fip_ejv.eject_time = uoe_adt2t68(fex_ejv.eject_time) - T68IRASBase
	   do j=1,5
	      do k=1,2
	         atime(k) = fex_ejv.start_times(k,j)
	         btime(k) = fex_ejv.stop_times(k,j)
	      enddo
	      fip_ejv.start_times(j) = uoe_adt2t68(atime) - T68IRASBase
	      fip_ejv.stop_times(j)  = uoe_adt2t68(btime) - T68IRASBase
	   enddo

c
c  Write the correction processing parameters to the FIP RDL.
c
	   do j=1,5
	      fip_ejv.corr_index(j) = fex_ejv.corr_index(j)
	      fip_ejv.tau(j)        = fex_ejv.tau(j)
	      fip_ejv.vib_corr(j)   = fex_ejv.vib_corr(j)
	   enddo

c
c  Write the correction spectra to the FIP RDL, performing the frequency cut
c	and the renormalization.
c
	   do j = fcc_jlo,fcc_jhi
	      do k = 1,5
	         fip_ejv.rcorr_spec(j-fcc_freq_offset,k) =
     .				real(fex_ejv.corr_spec(j,k)) * fac_erg_to_mjy
	         fip_ejv.icorr_spec(j-fcc_freq_offset,k) =
     .				aimag(fex_ejv.corr_spec(j,k)) * fac_erg_to_mjy
	      enddo
	   enddo


C
C  Write the FIP RDL to the project dataset file.
C

c
c  Get the project dataset file name.
c
	   fcc_fip_ofile = 'CSDR$FIRAS_OUT:FIP_EJV_' // fcc_scan_mode // '.' //
     .			    fcc_file_ext
	   call str$trim (fcc_fip_ofile, fcc_fip_ofile, fcc_fiplen)

c
c  Open the project dataset file.
c
	   rstatus = fut_get_lun(fip_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	   open (unit=fip_lun, file=fcc_fip_ofile, status='new',
     .		 form='unformatted', recordtype='fixed', recl=fcc_fip_osize,
     .		 access='sequential', iostat=io_stat)

	   if (io_stat .eq. 0) then
c
c  Write the correction spectra to the project dataset file.
c
	      write (fip_lun, iostat=io_stat) fip_ejv
	      if (io_stat .ne. 0) then
	         status = %loc(fip_rmswrite)
	         call lib$signal (fip_rmswrite, %val(2),
     .				  fcc_fip_ofile(1:fcc_fiplen), %val(io_stat))
	      endif

c
c  Close the project dataset file.
c
	      close (unit=fip_lun, iostat=io_stat)
	      if (io_stat .ne. 0) then
	         status = %loc(fip_rmsclose)
	         call lib$signal (fip_rmsclose, %val(2),
     .				  fcc_fip_ofile(1:fcc_fiplen), %val(io_stat))
	      endif

	      rstatus = fut_free_lun(fip_lun)
	      if (rstatus .ne. %loc(fut_normal)) then
	         call lib$signal (%val(rstatus))
	      endif

	   else
	      status = %loc(fip_rmsopen)
	      call lib$signal (fip_rmsopen, %val(2),
     .			       fcc_fip_ofile(1:fcc_fiplen), %val(io_stat))
	   endif	!  (open status

	endif		!  (status from FEX file I/O


	fip_reformat_ejv = status

	return
	end
