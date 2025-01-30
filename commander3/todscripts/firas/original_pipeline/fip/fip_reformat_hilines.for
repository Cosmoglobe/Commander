	integer * 4 function  fip_reformat_hilines ()

c-------------------------------------------------------------------------------
c
c	Function FIP_REFORMAT_HILINES
c
c	This function reads the symthetic FIRAS line profiles from the binary
c	file CSDR$FIRAS_IN:FLA_HLP_CCSS.VVV_XXXXXXXXXX for the high frequency
c	channels.  It performs the frequency cut, renormalizes the profiles to
c	MJr/sr, and then writes the profiles to the FIRAS Initial Product file
c	CSDR$FIRAS_OUT:FIP_HLP_CCSS.VVV_XXXXXXXXXX.
c
c	The input line profiles have the units of 1/icm.  Consequently, the
c	output line profiles are converted to units of 1/GHz.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  22 November 1993
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
c		fip_invoc_lines.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fip_invoc_lines)'
	include '(fip_frequency)'

	complex * 8	buffer(167)		!  line profile buffer

	integer * 4	fex_lun			!  fex file logical unit number
	integer * 4	fip_lun			!  fip file logical unit number
	integer * 4	io_stat			!  I/O return status
	integer * 4	is			!  scan mode pointer
	integer * 4	j			!  a counter
	integer * 4	rstatus			!  return status
	integer * 4	start			!  lower profile index
	integer * 4	status			!  return status

	real	* 4	rco_115(180)		!  real  3.845 icm  CO    line
	real	* 4	ico_115(180)		!  imag  3.845 icm  CO    line
	real	* 4	rco_230(180)		!  real  7.690 icm  CO    line
	real	* 4	ico_230(180)		!  imag  7.690 icm  CO    line
	real	* 4	rco_345(180)		!  real 11.535 icm  CO    line
	real	* 4	ico_345(180)		!  imag 11.535 icm  CO    line
	real	* 4	ro2_424(180)		!  real 14.168 icm  O2    line
	real	* 4	io2_424(180)		!  imag 14.168 icm  O2    line
	real	* 4	rco_461(180)		!  real 15.379 icm  CO    line
	real	* 4	ico_461(180)		!  imag 15.379 icm  CO    line
	real	* 4	rc_i_492(180)		!  real 15.379 icm [C I]  line
	real	* 4	ic_i_492(180)		!  imag 15.379 icm [C I]  line
	real	* 4	rh2o_556(180)		!  real 15.379 icm  H2O   line
	real	* 4	ih2o_556(180)		!  imag 15.379 icm  H2O   line
	real	* 4	rco_576(180)		!  real 19.222 icm  CO    line
	real	* 4	ico_576(180)		!  imag 19.222 icm  CO    line

	real	* 4	rco_691(180)		!  real 23.065 icm  CO    line
	real	* 4	ico_691(180)		!  imag 23.065 icm  CO    line
	real	* 4	rc_i_809(180)		!  real 27.00  icm [C I]  line
	real	* 4	ic_i_809(180)		!  imag 27.00  icm [C I]  line
	real	* 4	rh2o_1113(180)		!  real 37.136 icm  H2O   line
	real	* 4	ih2o_1113(180)		!  imag 37.136 icm  H2O   line
	real	* 4	rn_ii_1461(180)		!  real 48.738 icm [N II] line
	real	* 4	in_ii_1461(180)		!  imag 48.738 icm [N II] line
	real	* 4	rc_ii_1900(180)		!  real 63.395 icm [C II] line
	real	* 4	ic_ii_1900(180)		!  imag 63.395 icm [C II] line
	real	* 4	ro_i_2060(180)		!  real 68.716 icm [O I]  line
	real	* 4	io_i_2060(180)		!  imag 68.716 icm [O I]  line
	real	* 4	rsi_i_2311(180)		!  real 77.11  icm [Si I] line
	real	* 4	isi_i_2311(180)		!  imag 77.11  icm [Si I] line
	real	* 4	rn_ii_2459(180)		!  real 82.036 icm [N II] line
	real	* 4	in_ii_2459(180)		!  imag 82.036 icm [N II] line

	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fip_hlp'
	record /fip_hlp/ fip_lin

	external	fip_normal
	external	fip_rmsclose
	external	fip_rmsopen
	external	fip_rmsread
	external	fip_rmswrite
	external	fut_normal

c
c  Initialize the routine.
c
	is = 3*(fcc_chan-1) + (fcc_speed+1) + fcc_length
	start = lolim(is)-1
	call lib$movc5 (0,,0,720,rco_115)
	call lib$movc5 (0,,0,720,ico_115)
	call lib$movc5 (0,,0,720,rco_230)
	call lib$movc5 (0,,0,720,ico_230)
	call lib$movc5 (0,,0,720,rco_345)
	call lib$movc5 (0,,0,720,ico_345)
	call lib$movc5 (0,,0,720,ro2_424)
	call lib$movc5 (0,,0,720,io2_424)
	call lib$movc5 (0,,0,720,rco_461)
	call lib$movc5 (0,,0,720,ico_461)
	call lib$movc5 (0,,0,720,rc_i_492)
	call lib$movc5 (0,,0,720,ic_i_492)
	call lib$movc5 (0,,0,720,rh2o_556)
	call lib$movc5 (0,,0,720,ih2o_556)
	call lib$movc5 (0,,0,720,rco_576)
	call lib$movc5 (0,,0,720,ico_576)
	call lib$movc5 (0,,0,720,rco_691)
	call lib$movc5 (0,,0,720,ico_691)
	call lib$movc5 (0,,0,720,rc_i_809)
	call lib$movc5 (0,,0,720,ic_i_809)
	call lib$movc5 (0,,0,720,rh2o_1113)
	call lib$movc5 (0,,0,720,ih2o_1113)
	call lib$movc5 (0,,0,720,rn_ii_1461)
	call lib$movc5 (0,,0,720,in_ii_1461)
	call lib$movc5 (0,,0,720,rc_ii_1900)
	call lib$movc5 (0,,0,720,ic_ii_1900)
	call lib$movc5 (0,,0,720,ro_i_2060)
	call lib$movc5 (0,,0,720,io_i_2060)
	call lib$movc5 (0,,0,720,rsi_i_2311)
	call lib$movc5 (0,,0,720,isi_i_2311)
	call lib$movc5 (0,,0,720,rn_ii_2459)
	call lib$movc5 (0,,0,720,in_ii_2459)
	call lib$movc5 (0,,0,11576,fip_lin)
	status = %loc(fip_normal)


C
C  Read the line profiles from the binary file.
C

c
c  Find the line profile file name.
c
	fcc_fex_file = 'CSDR$FIRAS_IN:FLA_HLP_' // fcc_scan_mode // '.' //
     .			fcc_lines_ext
	call str$trim (fcc_fex_file, fcc_fex_file, fcc_fexlen)

c
c  Open the line profiles file.
c
	rstatus = fut_get_lun(fex_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

	open (unit=fex_lun, file=fcc_fex_file, status='old', form='unformatted',
     .	      access='sequential', readonly, iostat=io_stat)

	if (io_stat .eq. 0) then
c
c  Read the line profiles from the file.
c
	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rco_115(j) = real(buffer(j-start))
	      ico_115(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rco_230(j) = real(buffer(j-start))
	      ico_230(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rco_345(j) = real(buffer(j-start))
	      ico_345(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      ro2_424(j) = real(buffer(j-start))
	      io2_424(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rco_461(j) = real(buffer(j-start))
	      ico_461(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rc_i_492(j) = real(buffer(j-start))
	      ic_i_492(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rh2o_556(j) = real(buffer(j-start))
	      ih2o_556(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rco_576(j) = real(buffer(j-start))
	      ico_576(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rco_691(j) = real(buffer(j-start))
	      ico_691(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rc_i_809(j) = real(buffer(j-start))
	      ic_i_809(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rh2o_1113(j) = real(buffer(j-start))
	      ih2o_1113(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rn_ii_1461(j) = real(buffer(j-start))
	      in_ii_1461(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rc_ii_1900(j) = real(buffer(j-start))
	      ic_ii_1900(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      ro_i_2060(j) = real(buffer(j-start))
	      io_i_2060(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rsi_i_2311(j) = real(buffer(j-start))
	      isi_i_2311(j) = aimag(buffer(j-start))
	   enddo

	   call lib$movc5 (0,,0,1336,buffer)
	   read (fex_lun, iostat=io_stat) buffer
	   do j=lolim(is),uplim(is)
	      rn_ii_2459(j) = real(buffer(j-start))
	      in_ii_2459(j) = aimag(buffer(j-start))
	   enddo

	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmsread)
	      call lib$signal (fip_rmsread, %val(2), fcc_fex_file(1:fcc_fexlen),
     .					    %val(io_stat))
	   endif

c
c  Close the line profiles file.
c
	   close (unit=fex_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(fip_rmsclose)
	      call lib$signal (fip_rmsclose, %val(2),
     .			       fcc_fex_file(1:fcc_fexlen), %val(io_stat))
	   endif

	   rstatus = fut_free_lun(fex_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	else
	   status = %loc(fip_rmsopen)
	   call lib$signal (fip_rmsopen, %val(2), fcc_fex_file(1:fcc_fexlen),
     .					 %val(io_stat))
	endif	!  (open status


	if (status .eq. %loc(fip_normal)) then
C
C  Reformat the line profiles.
C

c
c  Fill in the line profile identification fields.
c
	   fip_lin.model_label(1:fcc_extlen) = fcc_lines_ext
	   fip_lin.chanscan                  = fcc_scan_mode
	   fip_lin.nu_zero                   = fcc_nu0
	   fip_lin.delta_nu                  = fcc_dnu
	   fip_lin.num_freq                  = fcc_nfreq

c
c  Write the line profiles to the RDL, performing the frequency cut and the
c	renormalization.
c
	   do j = fcc_jlo,fcc_jhi
c
c  3.845 icm CO J 1-0 line.
c
	      fip_lin.rco_115(j-fcc_freq_offset) = rco_115(j)/fac_icm_ghz
	      fip_lin.ico_115(j-fcc_freq_offset) = ico_115(j)/fac_icm_ghz
c
c  7.690 icm CO J 2-1 line.
c
	      fip_lin.rco_230(j-fcc_freq_offset) = rco_230(j)/fac_icm_ghz
	      fip_lin.ico_230(j-fcc_freq_offset) = ico_230(j)/fac_icm_ghz
c
c  11.535 icm CO J 3-2 line.
c
	      fip_lin.rco_345(j-fcc_freq_offset) = rco_345(j)/fac_icm_ghz
	      fip_lin.ico_345(j-fcc_freq_offset) = ico_345(j)/fac_icm_ghz
c
c  14.168 icm O2 line.
c
	      fip_lin.ro2_424(j-fcc_freq_offset) = ro2_424(j)/fac_icm_ghz
	      fip_lin.io2_424(j-fcc_freq_offset) = io2_424(j)/fac_icm_ghz
c
c  15.378 icm CO J 4-3 line.
c
	      fip_lin.rco_461(j-fcc_freq_offset) = rco_461(j)/fac_icm_ghz
	      fip_lin.ico_461(j-fcc_freq_offset) = ico_461(j)/fac_icm_ghz
c
c  16.419 icm [C I] line.
c
	      fip_lin.rc_i_492(j-fcc_freq_offset) = rc_i_492(j)/fac_icm_ghz
	      fip_lin.ic_i_492(j-fcc_freq_offset) = ic_i_492(j)/fac_icm_ghz
c
c  18.576 icm H2O line.
c
	      fip_lin.rh2o_556(j-fcc_freq_offset) = rh2o_556(j)/fac_icm_ghz
	      fip_lin.ih2o_556(j-fcc_freq_offset) = ih2o_556(j)/fac_icm_ghz
c
c  19.222 icm CO J 6-5 line.
c
	      fip_lin.rco_576(j-fcc_freq_offset) = rco_576(j)/fac_icm_ghz
	      fip_lin.ico_576(j-fcc_freq_offset) = ico_576(j)/fac_icm_ghz
c
c  23.065 icm CO J 6-5 line.
c
	      fip_lin.rco_691(j-fcc_freq_offset) = rco_691(j)/fac_icm_ghz
	      fip_lin.ico_691(j-fcc_freq_offset) = ico_691(j)/fac_icm_ghz
c
c  27.000 icm [C I] line.
c
	      fip_lin.rc_i_809(j-fcc_freq_offset) = rc_i_809(j)/fac_icm_ghz
	      fip_lin.ic_i_809(j-fcc_freq_offset) = ic_i_809(j)/fac_icm_ghz
c
c  37.136 icm H2O line.
c
	      fip_lin.rh2o_1113(j-fcc_freq_offset) = rh2o_1113(j)/fac_icm_ghz
	      fip_lin.ih2o_1113(j-fcc_freq_offset) = ih2o_1113(j)/fac_icm_ghz
c
c  48.738 icm [N II] line.
c
	      fip_lin.rn_ii_1461(j-fcc_freq_offset) = rn_ii_1461(j)/fac_icm_ghz
	      fip_lin.in_ii_1461(j-fcc_freq_offset) = in_ii_1461(j)/fac_icm_ghz
c
c  63.395 icm [C II] line.
c
	      fip_lin.rc_ii_1900(j-fcc_freq_offset) = rc_ii_1900(j)/fac_icm_ghz
	      fip_lin.ic_ii_1900(j-fcc_freq_offset) = ic_ii_1900(j)/fac_icm_ghz
c
c  68.716 icm [O I] line.
c
	      fip_lin.ro_i_2060(j-fcc_freq_offset) = ro_i_2060(j)/fac_icm_ghz
	      fip_lin.io_i_2060(j-fcc_freq_offset) = io_i_2060(j)/fac_icm_ghz
c
c  77.110 icm [Si I] line.
c
	      fip_lin.rsi_i_2311(j-fcc_freq_offset) = rsi_i_2311(j)/fac_icm_ghz
	      fip_lin.isi_i_2311(j-fcc_freq_offset) = isi_i_2311(j)/fac_icm_ghz
c
c  82.036 icm [N II] line.
c
	      fip_lin.rn_ii_2459(j-fcc_freq_offset) = rn_ii_2459(j)/fac_icm_ghz
	      fip_lin.in_ii_2459(j-fcc_freq_offset) = in_ii_2459(j)/fac_icm_ghz
	   enddo


C
C  Write the line profile RDL to the initial product file.
C

c
c  Get the initial product file name.
c
	   fcc_fip_file = 'CSDR$FIRAS_OUT:FIP_HLP_' // fcc_scan_mode // '.' //
     .			   fcc_lines_ext
	   call str$trim (fcc_fip_file, fcc_fip_file, fcc_fiplen)

c
c  Open the initial product file.
c
	   rstatus = fut_get_lun(fip_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	   open (unit=fip_lun, file=fcc_fip_file, status='new',
     .		 form='unformatted', recordtype='fixed', recl=fcc_hifip_size,
     .		 access='sequential', iostat=io_stat)

	   if (io_stat .eq. 0) then
c
c  Write the line profiles to the initial product file.
c
	      write (fip_lun, iostat=io_stat) fip_lin
	      if (io_stat .ne. 0) then
	         status = %loc(fip_rmswrite)
	         call lib$signal (fip_rmswrite, %val(2),
     .				  fcc_fip_file(1:fcc_fiplen), %val(io_stat))
	      endif

c
c  Close the initial product file.
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

	endif		!  (status from FEX file I/O


	fip_reformat_hilines = status

	return
	end
