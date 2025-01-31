	integer * 4 function  fsl_open_spectra (vs_lun, cs_lun)

c-------------------------------------------------------------------------------
c
c	Function FSL_OPEN_SPECTRA
c
c	This function opens the spectra archives CSDR$FIRAS_OUT for either sky
c	or calibration data. The command line qualifiers allow selections of
c       writing voltage spectra and either calibrated spectra or differential
c       spectra to the archives. This function opens voltage spectra and
c       either calibrated or differential spectra if requested.
c
c	Author:
c                FSL_Open_Spectra
c                Shirley M. Read
c                Hughes STX Corporation
c                July 27, 1995
c
c-------------------------------------------------------------------------------
c
c	Input:
c		No parameters passed. Input is obtained from include files.
c
c	Output:
c		vs_lun		integer * 4		file lun for voltage
c                                                       spectra
c		cs_lun		integer * 4		file lun for calibrated
c                                                       or differential spectra
c
c	Subroutines called:
c               FORTRAN open with ct_connect_write and csa_open_skymap
c               csa_field_offset_values
c		fut_get_lun
c		lib$signal
c
c	Include files:
c		fsl_config.txt
c		fsl_invoc.txt
c		fut_params.txt
c-------------------------------------------------------------------------------
c	Changes for FSL:
c
c	Fred Shuman, Hughes STX, 1995 Sep 26
c	    1. "fac_coad_spec_pix_offset" replaced with
c	       "fac_coad_spec_pix_offsetL".
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fsl_config)'
	include '(fsl_invoc)'

	integer * 4	vs_lun		!  voltage spectra logical unit
	integer * 4	cs_lun		!  calibrated spectra logical unit

	character * 48	out_file	!  Cobetrieve file name

	integer * 4	io_stat		!  open return status
	integer * 4	rstatus		!  return status
	integer * 4	cstatus		!  CSA return status
	integer * 4	status		!  return status
	integer * 2	ct_stat(20)	!  CT return status
	integer * 2     outlen		!  length of file name

	integer * 4	ct_connect_write
	external	ct_connect_write
	integer * 4	csa_open_skymap
	external	csa_open_skymap

	integer * 4	fut_get_lun
	integer * 4	csa_field_offset_values

	external	fsl_ctspecnum
	external	fsl_ctspecopen
	external	fsl_csafldoffset
	external	fsl_csaspecopen
	external	fsl_normal
	external	csa_normal
	external	fut_normal
c
c  Initialize function.
c
	status = %loc(fsl_normal)
c
c  Determine if voltage spectra are to be archived. If so,
c  open the Cobetrieve spectrum file.
c
	vs_lun = 0
	if ( fcc_write_vs .eq. fac_present) then
	   rstatus = fut_get_lun(vs_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif

	   out_file(1:fcc_outlenvs) = fcc_outfile_vs(1:fcc_outlenvs)

	   if ((fcc_cal .eq. fac_present) .and.        ! open the CT TOD file
     .         (rstatus .eq. %loc(fut_normal))) then    ! for calibration spectra
	      open (unit=vs_lun, file=out_file, iostat=io_stat,
     .		    status='new', useropen=ct_connect_write)
	      if (io_stat .ne. 0 ) then
	         status = %loc(fsl_ctspecopen)
	         call lib$signal (fsl_ctspecopen, %val(2),
     .                            out_file(1:fcc_outlenvs), %val(io_stat))
	      endif
	   elseif (rstatus .eq. %loc(fut_normal)) then	! open spectra skymap

	      out_file(1:fcc_outlenvs) = fcc_outfile_vs(1:fcc_outlenvs)
	      open (unit=vs_lun, file=out_file, iostat=io_stat, status='new',
     .		 form='unformatted', recordtype='fixed', recl=fcc_record_size,
     .		 useropen=csa_open_skymap)
	      if (io_stat .ne. 0 ) then
		 status = %loc(fsl_csaspecopen)
		 call lib$signal (fsl_csaspecopen, %val(2),
     .                            out_file(1:fcc_outlenvs), %val(io_stat))
	      else
c
c  Set up the skymap offset values.
c
	         cstatus = csa_field_offset_values (fac_coad_spec_pix_offsetL,
     .						 fac_time_offset, -1, vs_lun)
	         if (cstatus .ne. %loc(csa_normal)) then
	            status = %loc(fsl_csafldoffset)
	            call lib$signal (fsl_csafldoffset, %val(2),
     .				     out_file(1:fcc_outlenvs), %val(cstatus))
	         endif
	      endif     ! open sky voltage spectra file
	   endif	! open cal or sky voltage spectra file
	endif	! fcc_write_vs is requested
c
c  Determine if differential or calibrated spectra are to be archuived. If so,
c  open the Cobetrieve spectrum file.
c
	cs_lun = 0
	if ((fcc_write_ds .eq. fac_present). or.
     .       (fcc_write_cs .eq. fac_present)) then
	   rstatus = fut_get_lun(cs_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   endif
	   if (fcc_write_ds .eq. fac_present) then
	      out_file(1:fcc_outlends) = fcc_outfile_ds(1:fcc_outlends)
	      outlen = fcc_outlends
	   elseif (fcc_write_cs .eq. fac_present) then
	      out_file(1:fcc_outlencs) = fcc_outfile_cs(1:fcc_outlencs)
	      outlen = fcc_outlencs
	   endif

	   if ((fcc_cal .eq. fac_present) .and.        ! open the CT TOD file
     .         (rstatus .eq. %loc(fut_normal))) then
	      open (unit=cs_lun, file=out_file, iostat=io_stat,
     .		    status='new', useropen=ct_connect_write)
	      if (io_stat .ne. 0 ) then
	         status = %loc(fsl_ctspecopen)
	         call lib$signal (fsl_ctspecopen, %val(2),
     .                            out_file(1:outlen), %val(io_stat))
	      endif

	   elseif (rstatus .eq. %loc(fut_normal)) then	! open spectra skymap

	      open (unit=cs_lun, file=out_file, iostat=io_stat, status='new',
     .		 form='unformatted', recordtype='fixed', recl=fcc_record_size,
     .		 useropen=csa_open_skymap)
	      if (io_stat .ne. 0 ) then
		 status = %loc(fsl_csaspecopen)
		 call lib$signal (fsl_csaspecopen, %val(2),
     .                            out_file(1:outlen), %val(io_stat))
	      else
c
c  Set up the skymap offset values.
c
	         cstatus = csa_field_offset_values (fac_coad_spec_pix_offsetL,
     .						 fac_time_offset, -1, cs_lun)
	         if (cstatus .ne. %loc(csa_normal)) then
	            status = %loc(fsl_csafldoffset)
	            call lib$signal (fsl_csafldoffset, %val(2),
     .				  out_file(1:outlen), %val(cstatus))
	         endif
	      endif     ! open sky calibrated or differential spectra file
	   endif   ! open cal or sky calibrated or differential spectra file
	endif	! fcc_write_ds or fcc_write_cs  is requested


	fsl_open_spectra = status

	return
	end
