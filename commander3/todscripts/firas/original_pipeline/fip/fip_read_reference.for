	integer * 4 function  fip_read_reference ()

c-------------------------------------------------------------------------------
c
c	Function FIP_READ_REFERENCE
c
c	This function reads the reference datasets for FIP_MODEL.  The routine
c	initializes Cobetrieve.  It then uses the open config, get config, and
c	close config routines to get the Nyquist frequency, the apodization
c	function, and the electronics transfer function for FIP_MODEL.  The
c	routine also rotates the apodization function to that the ifg peak
c	position is in the first sample.
c
c
c	Author:  Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 25 May 1993
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
c		cct_close_config
c		cct_get_config_idx_tod
c		cct_get_config_tod
c		cct_open_config
c		ct_init
c		fut_apod_recnum
c		fut_get_recnum
c		lib$signal
c
c	Include files:
c		ct$library:ctuser.inc
c		fip_config_model.txt
c		fip_invoc_model.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include 'ct$library:ctuser.inc'
	include '(fip_config_model)'
	include '(fip_invoc_model)'
	include '(fut_params)'

	integer * 2	ct_stat(20)	!  CT_Init return status

	integer * 4	cstatus		!  CT return status
	integer * 4	arecno		!  apodization function record number
	integer * 4	erecno		!  ETF record number
	integer * 4	fakeit /0/	!  fakeit flag
	integer * 4	frecno		!  Nyquist frequency record number
	integer * 4	ictr		!  ifg peak position
	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	k		!  a counter
	integer * 4	linearized /1/	!  linearization status flag
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status
	integer * 4	upmode /4/	!  microprocessor mode flag

	real	* 8	apod(512)	!  unrotated apodization function

	integer * 4	cct_close_config
	integer * 4	cct_get_config_idx_tod
	integer * 4	cct_get_config_tod
	integer * 4	cct_open_config
	integer * 4	fut_apod_recnum
	integer * 4	fut_default_peak
	integer * 4	fut_get_recnum

	external	cct_normal
	external	fip_cfgseqclose
	external	fip_cfgseqget
	external	fip_cfgseqopen
	external	fip_cfgdirclose
	external	fip_cfgdirget
	external	fip_cfgdiropen
	external	fip_ctinit
	external	fip_normal
	external	fip_readapod
	external	fip_readetf

c
c  Initialize Cobetrieve.
c
	status = %loc(fip_normal)
	call ct_init(ct_stat)

	if (ct_stat(1) .eq. ctp_normal) then
C
C  Get the sequential-access reference dataset.
C

	   cstatus = cct_open_config (ref_start, ref_stop, ndset_tod,
     .				      dset_tod, size_tod, seq_access, ncache,
     .				      tod_lun, tod_index, tod_stat, ref_count)
	   if (cstatus .ne. %loc(cct_normal)) then
	      status = %loc(fip_cfgseqopen)
	      call lib$signal (fip_cfgseqopen, %val(1), %val(cstatus))
	   else
	      cstatus = cct_get_config_tod (ref_time, ndset_tod, size_tod,
     .					    tod_lun, tod_index, nyquist,
     .				            new_tod_segment, tod_stat)
	      if (cstatus .ne. %loc(cct_normal)) then
	         status = %loc(fip_cfgseqget)
	         call lib$signal (fip_cfgseqget, %val(1), %val(cstatus))
	      else
c
c  Get the Nyquist frequency
c
	         frecno = 4*jmod((fcc_chan-1),2) + fcc_smode
	         fnyq_icm  = nyquist.icm(frecno)
	         fnyq_hz   = nyquist.hz(frecno)
	         spec_norm = fac_watt_to_mjy / 
     .                       (fnyq_icm * fac_etendu * fac_adc_scale)
	      endif
	      cstatus = cct_close_config (ndset_tod, tod_lun, tod_index)
	      if (cstatus .ne. %loc(cct_normal)) then
	         status = %loc(fip_cfgseqclose)
	         call lib$signal (fip_cfgseqclose, %val(1), %val(cstatus))
	      endif
	   endif


	   if (status .eq. %loc(fip_normal)) then
C
C  Get the direct-access reference dataset.
C

	      cstatus = cct_open_config (ref_start, ref_stop, ndset_dir,
     .                                   dset_dir, size_dir, dir_access, 
     .                                   ncache, dir_lun, dir_index,
     .                                   dir_stat, ref_count)
	      if (cstatus .ne. %loc(cct_normal)) then
	         status = %loc(fip_cfgdiropen)
	         call lib$signal (fip_cfgdiropen, %val(1), %val(cstatus))
	      else
	         cstatus = cct_get_config_idx_tod (ref_time, ndset_dir,
     .						   dir_lun, dir_index,
     .						   new_dir_segment, dir_stat)
	         if (cstatus .ne. %loc(cct_normal)) then
	            status = %loc(fip_cfgdirget)
	            call lib$signal (fip_cfgdirget, %val(1), %val(cstatus))
	         else

c
c  Read the apodization function.
c
	            rstatus = fut_apod_recnum (fcc_speed, fcc_ngroup, fakeit,
     .				               fcc_length, fcc_chan, upmode,
     .					       linearized, arecno)
	            read (dir_lun(1)'arecno, iostat=io_stat) apod
	            if (io_stat .ne. 0) then
	               status = %loc(fip_readapod)
	               call lib$signal (fip_readapod, %val(2), %val(arecno),
     .						      %val(io_stat))
	            else
c
c  Rotate the apodization function.
c
	               rstatus = fut_default_peak (fcc_speed, fcc_length,
     .						   fcc_chan, fcc_ngroup, upmode,
     .						   linearized, ictr)
	               call lib$movc5 (0,,0,4096,apod_fcn)
	               do j = 1,512
	                  k = j + 1 - ictr
	                  if (k .le. 0) k = k + 512
	                  apod_fcn(k) = apod(j)
	               enddo

c
c  Read the ETF.
c
	               call fut_get_recnum (fakeit, fcc_speed, fcc_chan, upmode,
     .					    fcc_ngroup, erecno)
	               read (dir_lun(2)'erecno, iostat=io_stat) ztrans
	               if (io_stat .ne. 0) then
	                  status = %loc(fip_readetf)
	                  call lib$signal (fip_readetf, %val(1), %val(erecno),
     .							%val(io_stat))
	               endif
	            endif
	         endif
	         cstatus = cct_close_config (ndset_dir, dir_lun, dir_index)
	         if (cstatus .ne. %loc(cct_normal)) then
	            status = %loc(fip_cfgdirclose)
	            call lib$signal (fip_cfgdirclose, %val(1), %val(cstatus))
	         endif
	      endif

	   endif	!  status from tod reference dataset

	else
	   status = %loc(fip_ctinit)
	   call lib$signal (fip_ctinit, %val(1), %val(ct_stat(1)))
	endif	!  status from ct_init


	fip_read_reference = status

	return
	end
