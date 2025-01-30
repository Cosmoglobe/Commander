	integer * 4 function  fsl_open_config (ref_open_seq, ref_open_dir)

c-------------------------------------------------------------------------------
c
c	Function FSL_OPEN_CONFIG
c
c	This function opens the configuration files for FIRAS reference data
c       sets in CSDR$FIRAS_REF: FEX_GRTCOAWT, FEX_GRTTRANS, FEX_NYQUISTL,
c       FEX_VIBCORRL, FEX_APODL, and FEX_ETFL using the COBEtrieve routine,
c       CCT_OPEN_CONFIG.
c
c	Author:	 
c                FSL_Open_Config
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
c
c-------------------------------------------------------------------------------
c
c	Input:
c		No input parameters. Data required for call to CCT_Open_Config
c               is pre-stored in the FSL_Config include file.
c
c	Output:
c		ref_open_seq	logical * 1	sequential access reference
c						   data set open flag
c		ref_open_dir	logical * 1	direct access reference
c						   data set open flag
c
c		Data returned from CCT_Open_Config is put into the FSL_Config
c               include file. 
c
c	Subroutines called:
c		cct_open_config
c		lib$signal
c
c	Include files:
c		fsl_config.txt
c
c-------------------------------------------------------------------------------
c
c       Changes for FSL:
c
c-------------------------------------------------------------------------------


	implicit none

	include '(fsl_config)'

	logical *  1	ref_open_dir		!  reference dataset open flag
	logical *  1	ref_open_seq		!  reference dataset open flag

	integer *  4    status                  !  processing status
	integer *  4    cstatus                 !  CCT_Open_Config return status

	integer *  4	cct_open_config

	external	fsl_cfgdiropen
	external	fsl_cfgseqopen
	external        fsl_normal
	external        cct_normal
c
c  Initialize function.
c
	status = %loc(fsl_normal)
c
c  Open the sequential access configuration reference data sets.
c
	call ct_gmt_to_binary(ref_gmt_start, ref_start)
	call ct_gmt_to_binary(ref_gmt_stop, ref_stop)

	cstatus = cct_open_config (ref_start, ref_stop, ndset_tod, dset_tod, 
     .                             size_tod, seq_access, ncache, tod_lun, 
     .                             tod_index, tod_stat, ref_count)
	if (cstatus .ne. %loc(cct_normal)) then
	     status = %loc(fsl_cfgseqopen)
	     call lib$signal (fsl_cfgseqopen, %val(1), %val(cstatus))
	else
	     ref_open_seq = .true.
c
c  Open the direct access configuration reference data sets.
c

	     cstatus = cct_open_config (ref_start, ref_stop, ndset_dir,
     .                                  dset_dir, size_dir, dir_access, ncache,
     .                                  dir_lun, dir_index, dir_stat, ref_count)
	     if (cstatus .ne. %loc(cct_normal)) then
	         status = %loc(fsl_cfgdiropen)
	         call lib$signal (fsl_cfgdiropen, %val(1), %val(cstatus))
	     else
	         ref_open_dir = .true.
	     endif
	endif

	fsl_open_config = status

	return
	end	
