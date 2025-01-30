	integer * 4 function  fsl_close_config (ref_open_seq, ref_open_dir)

c-------------------------------------------------------------------------------
c
c	Function FSL_CLOSE_CONFIG
c
c	This function closes the configuration files for FIRAS reference data
c       sets in CSDR$FIRAS_REF: FEX_GRTCOAWT, FEX_GRTTRANS, FEX_NYQUISTL,
c       FEX_VIBCORRL, FEX_APODL, and FEX_ETFL. using the COBEtrieve routine,
c       CCT_CLOSE_CONFIG.
c
c	Author:	 
c                FSL_Close_Config
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
c
c-------------------------------------------------------------------------------
c
c	Input:
c		ref_open_seq	logical * 1	sequential access reference
c						   data set open flag
c		ref_open_dir	logical * 1	direct access reference
c						   data set open flag
c	Output:
c		ref_open_seq	logical * 1	sequential access reference
c						   data set open flag
c		ref_open_dir	logical * 1	direct access reference
c						   data set open flag
c
c	Subroutines called:
c		cct_close_config
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
	integer *  4    cstatus                 !  CCT_Close_Config return 
						!  status
	integer *  4	cct_close_config

	external	fsl_cfgdirclose
	external	fsl_cfgseqclose
	external        fsl_normal
	external        cct_normal
c
c  Initialize function.
c
	status = %loc(fsl_normal)
c
c  Close the sequential access configuration reference data sets.
c
	if (ref_open_seq) then
	    cstatus = cct_close_config (ndset_tod, tod_lun, tod_index)
	    if (cstatus .ne. %loc(cct_normal)) then
	        status = %loc(fsl_cfgseqclose)
	        call lib$signal (fsl_cfgseqclose, %val(1), %val(cstatus))
	    else
		ref_open_seq = .false.
	    endif
	endif
c
c  Close the direct access configuration reference data sets.
c
	if (ref_open_dir) then
	    cstatus = cct_close_config (ndset_dir, dir_lun, dir_index)
	    if (cstatus .ne. %loc(cct_normal)) then
	        status = %loc(fsl_cfgdirclose)
	        call lib$signal (fsl_cfgdirclose, %val(1), %val(cstatus))
	    else
		ref_open_dir = .false.
	    endif
	endif

	fsl_close_config = status

	return
	end	
