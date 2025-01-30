	integer * 4 function  ffi_get_reference (ref_open_dir, ref_open_seq)

c-------------------------------------------------------------------------------
c
c	Function FFI_GET_REFERENCE
c
c	This function opens and gets the reference datasets for FFI via the
c	Get_Config routines.  The function then computes the constants required
c	to process the data.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  22 September 1992
c		  SER 10763
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		ref_open_dir		logical * 1		open status flag
c		ref_open_seq		logical * 1		open status flag
c
c	Subroutines called:
c		cct_get_config_idx_tod
c		cct_get_config_tod
c		cct_open_config
c		ct_gmt_to_binary
c		ffi_compute_constants
c		lib$signal
c
c	Include files:
c		ffi_config.txt
c		ffi_invoc.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c-------------------------------------------------------------------------------

	implicit none

	include '(ffi_config)'
	include '(ffi_invoc)'

	integer * 4	cstatus			!  return status
	integer * 4	status			!  return status

	logical * 1	ref_open_dir		!  reference file open flag
	logical * 1	ref_open_seq		!  reference file open flag

	integer * 4	cct_get_config_idx_tod
	integer * 4	cct_get_config_tod
	integer * 4	cct_open_config
	integer * 4	ffi_compute_constants

	external	cct_normal
	external	ffi_cfgdirget
	external	ffi_cfgdiropen
	external	ffi_cfgseqget
	external	ffi_cfgseqopen
	external	ffi_normal

	status = %loc(ffi_normal)
c
c  Open the configuration files for sequential, and direct access.
c
	ref_open_dir = .false.
	ref_open_seq = .false.
	call ct_gmt_to_binary(ref_gmt_start,ref_start)
	call ct_gmt_to_binary(ref_gmt_stop,ref_stop)

	cstatus = cct_open_config (ref_start, ref_stop, ndset_tod,
     .                             dset_tod, size_tod, seq_access,
     .                             ncache, tod_lun, tod_index,
     .                             tod_stat, ref_count)
	if (cstatus .ne. %loc(cct_normal)) then
	   status = %loc(ffi_cfgseqopen)
	   call lib$signal (ffi_cfgseqopen, %val(1), %val(cstatus))
	else
	   ref_open_seq = .true.

	   cstatus = cct_open_config (ref_start, ref_stop, ndset_dir,
     .                                dset_dir, size_dir, dir_access, 
     .                                ncache, dir_lun, dir_index,
     .                                dir_stat, ref_count)
	   if (cstatus .ne. %loc(cct_normal)) then
	      status = %loc(ffi_cfgdiropen)
	      call lib$signal (ffi_cfgdiropen, %val(1), %val(cstatus))
	   else
	      ref_open_dir = .true.
	   endif
	endif


	if (status .eq. %loc(ffi_normal)) then
c
c  Get the reference datasets.
c

	   cstatus = cct_get_config_tod (fcc_jstart, ndset_tod, size_tod,
     .					 tod_lun, tod_index, config,
     .				         new_tod_segment, tod_stat)
	   if (cstatus .ne. %loc(cct_normal)) then
	      status = %loc(ffi_cfgseqget)
	      call lib$signal (ffi_cfgseqget, %val(1), %val(cstatus))
	   else

	      cstatus = cct_get_config_idx_tod (fcc_jstart, ndset_dir,
     .						dir_lun, dir_index,
     .						new_dir_segment, dir_stat)
	      if (cstatus .ne. %loc(cct_normal)) then
	         status = %loc(ffi_cfgdirget)
	         call lib$signal (ffi_cfgdirget, %val(1), %val(cstatus))
	      endif
	   endif


	   if (status .eq. %loc(ffi_normal)) then
c
c  Compute the constants required for processing the spectra.
c
	      status = ffi_compute_constants ()
	   endif

	endif	!   (status from Open_Config


	ffi_get_reference = status

	return
	end
