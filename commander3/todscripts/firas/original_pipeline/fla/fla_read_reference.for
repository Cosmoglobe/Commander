	integer * 4 function  fla_read_reference ()

c-------------------------------------------------------------------------------
c
c	Function FLA_READ_REFERENCE
c
c	This function reads the reference datasets for FLA_GAIN.  The routine
c	initializes Cobetrieve.  It then uses the open config, get config, and
c	close config routines to get the Nyquist frequency and the electronics
c	transfer function for FLA_GAIN.
c
c	Author:  Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 1 July 1993
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
c		fut_get_recnum
c		lib$signal
c
c	Include files:
c		ct$library:ctuser.inc
c		fla_config_gain.txt
c		fla_invoc_gain.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include 'ct$library:ctuser.inc'
	include '(fut_params)'
	include '(fla_config_gain)'
	include '(fla_invoc_gain)'

	integer * 2	ct_stat(20)	!  CT_Init return status

	integer * 4	cstatus		!  CT return status
	integer * 4	erecno		!  ETF record number
	integer * 4	fakeit /0/	!  fakeit flag
	integer * 4	frecno		!  Nyquist frequency record number
	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	status		!  return status
	integer * 4	upmode /4/	!  microprocessor mode flag

	real	* 8	df		!  frequency interval in icm
	real	* 8	dw		!  frequency interval in radians/sec
	real	* 8	fnyq_hz		!  Nyquist frequency in hz
	real	* 8	fnyq_icm	!  Nyquist frequency in icm

	integer * 4	cct_close_config
	integer * 4	cct_get_config_idx_tod
	integer * 4	cct_get_config_tod
	integer * 4	cct_open_config
	integer * 4	fut_get_recnum

	external	cct_normal
	external	fla_cfgseqclose
	external	fla_cfgseqget
	external	fla_cfgseqopen
	external	fla_cfgdirclose
	external	fla_cfgdirget
	external	fla_cfgdiropen
	external	fla_ctinit
	external	fla_normal
	external	fla_readetf

c
c  Initialize Cobetrieve.
c
	status = %loc(fla_normal)
	call ct_init(ct_stat)

	if (ct_stat(1) .eq. ctp_normal) then
C
C  Get the sequential-access reference dataset.
C

	   cstatus = cct_open_config (ref_start, ref_stop, ndset_tod,
     .				      dset_tod, size_tod, seq_access, ncache,
     .				      tod_lun, tod_index, tod_stat, ref_count)
	   if (cstatus .ne. %loc(cct_normal)) then
	      status = %loc(fla_cfgseqopen)
	      call lib$signal (fla_cfgseqopen, %val(1), %val(cstatus))
	   else
	      cstatus = cct_get_config_tod (ref_time, ndset_tod, size_tod,
     .					    tod_lun, tod_index, nyquist,
     .				            new_tod_segment, tod_stat)
	      if (cstatus .ne. %loc(cct_normal)) then
	         status = %loc(fla_cfgseqget)
	         call lib$signal (fla_cfgseqget, %val(1), %val(cstatus))
	      else
c
c  Get the Nyquist frequency.
c
	         frecno = 4*jmod((fcc_chan-1),2) + fcc_smode
	         fnyq_icm  = dble(nyquist.icm(frecno))
	         fnyq_hz   = dble(nyquist.hz(frecno))
	      endif
	      cstatus = cct_close_config (ndset_tod, tod_lun, tod_index)
	      if (cstatus .ne. %loc(cct_normal)) then
	         status = %loc(fla_cfgseqclose)
	         call lib$signal (fla_cfgseqclose, %val(1), %val(cstatus))
	      endif
	   endif
c
c  Compute the frequency constants.
c
	   df = fnyq_icm / 256.0D0
	   dw = fac_dpi * fnyq_hz / 128.0D0
	   do j = 1,257
	      freq(j)  = dble(j-1) * df
	      afreq(j) = dble(j-1) * dw
	   enddo
	   spec_norm = fnyq_icm * dble(fac_etendu) * dble(fac_adc_scale)



	   if (status .eq. %loc(fla_normal)) then
C
C  Get the direct-access reference dataset.
C

	      cstatus = cct_open_config (ref_start, ref_stop, ndset_dir,
     .                                   dset_dir, size_dir, dir_access, 
     .                                   ncache, dir_lun, dir_index,
     .                                   dir_stat, ref_count)
	      if (cstatus .ne. %loc(cct_normal)) then
	         status = %loc(fla_cfgdiropen)
	         call lib$signal (fla_cfgdiropen, %val(1), %val(cstatus))
	      else
	         cstatus = cct_get_config_idx_tod (ref_time, ndset_dir,
     .						   dir_lun, dir_index,
     .						   new_dir_segment, dir_stat)
	         if (cstatus .ne. %loc(cct_normal)) then
	            status = %loc(fla_cfgdirget)
	            call lib$signal (fla_cfgdirget, %val(1), %val(cstatus))
	         else
c
c  Read the ETF.
c
	            call fut_get_recnum (fakeit, fcc_speed, fcc_chan, upmode,
     .					 fcc_ngroup, erecno)
	            read (dir_lun'erecno, iostat=io_stat) ztrans
	            if (io_stat .ne. 0) then
	               status = %loc(fla_readetf)
	               call lib$signal (fla_readetf, %val(1), %val(erecno),
     .						     %val(io_stat))
	            endif
	         endif
	         cstatus = cct_close_config (ndset_dir, dir_lun, dir_index)
	         if (cstatus .ne. %loc(cct_normal)) then
	            status = %loc(fla_cfgdirclose)
	            call lib$signal (fla_cfgdirclose, %val(1), %val(cstatus))
	         endif
	      endif

	   endif	!  status from tod reference dataset

	else
	   status = %loc(fla_ctinit)
	   call lib$signal (fla_ctinit, %val(1), %val(ct_stat(1)))
	endif	!  status from ct_init


	fla_read_reference = status

	return
	end
