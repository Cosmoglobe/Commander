	integer * 4 function  fip_read_nyquist ()

c-------------------------------------------------------------------------------
c
c	Function FIP_READ_NYQUIST
c
c	This function reads the Nyquist frequency reference dataset for FIP
c	The routine initializes Cobetrieve.  It then uses the open config, get
c	config, and close config routines to get the Nyquist frequency dataset.
c
c	Author:  Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 1 June 1993
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
c		cct_get_config_tod
c		cct_open_config
c		ct_init
c		lib$signal
c
c	Include files:
c		ct$library:ctuser.inc
c		fip_config_freq.txt
c		fip_frequency.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Changed FIP_CONFIG_LINES.TXT to FIP_CONFIG_FREQ.TXT and moved fcc_chan
c	and fcc_smode information to FIP_FREQUENCY.TXT.
c	Gene Eplee, GSC, 10 February 1994.
c
c-------------------------------------------------------------------------------

	implicit none

	include 'ct$library:ctuser.inc'
	include '(fip_config_freq)'
	include '(fip_frequency)'

	integer * 2	ct_stat(20)	!  CT_Init return status

	integer * 4	cstatus		!  CT return status
	integer * 4	frecno		!  Nyquist frequency record number
	integer * 4	status		!  return status

	integer * 4	cct_close_config
	integer * 4	cct_get_config_tod
	integer * 4	cct_open_config

	external	cct_normal
	external	fip_cfgseqclose
	external	fip_cfgseqget
	external	fip_cfgseqopen
	external	fip_ctinit
	external	fip_normal

c
c  Initialize Cobetrieve.
c
	status = %loc(fip_normal)
	call ct_init(ct_stat)

	if (ct_stat(1) .eq. ctp_normal) then
c
c  Get the sequential-access reference dataset.
c
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
	         frecno = 4*jmod((fcc_fchan-1),2) + fcc_fsmode
	         fnyq_icm  = nyquist.icm(frecno)
	         fnyq_hz   = nyquist.hz(frecno)
	      endif
	      cstatus = cct_close_config (ndset_tod, tod_lun, tod_index)
	      if (cstatus .ne. %loc(cct_normal)) then
	         status = %loc(fip_cfgseqclose)
	         call lib$signal (fip_cfgseqclose, %val(1), %val(cstatus))
	      endif
	   endif

	else
	   status = %loc(fip_ctinit)
	   call lib$signal (fip_ctinit, %val(1), %val(ct_stat(1)))
	endif	!  status from ct_init


	fip_read_nyquist = status

	return
	end
