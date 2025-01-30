	integer * 4 function  fsl_close_spectra (vs_lun, cs_lun)

c-------------------------------------------------------------------------------
c
c	Function FSL_CLOSE_SPECTRA
c
c	This function closes the spectra archives in CSDR$FIRAS_OUT for either 
c       sky or calibration data, It closes archives for both voltage spectra 
c       and either calibration or differential spectra.
c
c	Author:	 
c               FSL_Close_Spectra
c               Shirley M. Read
c               Hughes STX Corporation
c               August 11 1995
c
c-------------------------------------------------------------------------------
c
c	Input:
c		vs_lun          integer * 4             voltage spectra lun
c		cs_lun          integer * 4             calibrated or 
c                                                       differential spectra lun
c
c	Output:
c               none
c
c	Subroutines called:
c               ct_close_arcv
c               csa_close_skymap
c               fut_free_lun
c		lib$signal
c
c	Include files:
c		fsl_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c       Changes for FSL:
c
c-------------------------------------------------------------------------------

	implicit none
	
	include '(fut_params)'
	include '(fsl_invoc)'
	include 'ct$library:ctuser.inc'

	integer * 4     vs_lun			! voltage spectra lun
	integer * 4     cs_lun			! calibrated or differential 
						!   spectra lun
	integer * 4	status			! function status
	integer * 4	cstatus			! csa status return
	integer * 4	rstatus			! fut_free_lun status return

	integer * 2	ct_stat(20)		! Cobetrieve status return
	integer * 2     outlen                  ! Cobetrieve file lemgth

	character * 48  out_file                ! Cobetrieve output file name

	integer * 4     csa_close_skymap
	integer * 4     fut_free_lun

	external	fsl_ctspecclose
	external        fsl_csaspecclose
	external        fsl_normal
	external	fut_normal
	external        csa_normal

c
c  Initialize function.
c
	status = %loc(fsl_normal)
c
c  Close the Cobetrieve voltage spectrum file.
c
	if ( fcc_write_vs .eq. fac_present ) then   ! voltage spectra written

	   if ( fcc_cal .eq. fac_present ) then     ! calibration data

	      call ct_close_arcv (,vs_lun, ct_stat)
	      if (ct_stat(1) .ne. ctp_normal) then
	         status = %loc(fsl_ctspecclose)
	         call lib$signal (fsl_ctspecclose, %val(2), 
     .                fcc_outfile_vs(1:fcc_outlenvs), %val(ct_stat(1)))
	      end if

	   else					    ! sky data

	      cstatus = csa_close_skymap (vs_lun, fac_skymap_no_levels)
	      if (cstatus .ne. %loc(csa_normal)) then
	         status = %loc(fsl_csaspecclose)
	         call lib$signal (fsl_csaspecclose, %val(2),
     .		      fcc_outfile_vs(1:fcc_outlenvs), %val(cstatus))
	      endif
	     
	   end if	! calibration or sky data written
c
c  Free the logical unit number.
c
	   rstatus = fut_free_lun(vs_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   end if
	end if		! voltage spectra written

c
c  Close the Cobetrieve calibrated or differential spectra file.
c
	if (( fcc_write_cs .eq. fac_present ) .or.  ! calibrated or differential
     .      ( fcc_write_ds .eq. fac_present )) then   ! spectra written
	   if ( fcc_write_cs .eq. fac_present ) then  ! calibrated spectra
              outlen = fcc_outlencs
	      out_file(1:outlen) = fcc_outfile_cs(1:outlen)
	   else					      ! differential spectra
	      outlen = fcc_outlends
	      out_file(1:outlen) = fcc_outfile_ds(1:outlen)
	   endif
	
	   if ( fcc_cal .eq. fac_present ) then     ! calibration data
	      
	      call ct_close_arcv (,cs_lun, ct_stat)
	      if (ct_stat(1) .ne. ctp_normal) then
	         status = %loc(fsl_ctspecclose)
	         call lib$signal (fsl_ctspecclose, %val(2), 
     .                out_file(1:outlen), %val(ct_stat(1)))
	      end if

	   else					    ! sky data

	      cstatus = csa_close_skymap (cs_lun, fac_skymap_no_levels)
	      if (cstatus .ne. %loc(csa_normal)) then
	         status = %loc(fsl_csaspecclose)
	         call lib$signal (fsl_csaspecclose, %val(2),
     .		      out_file(1:outlen), %val(cstatus))
	      endif
	     
	   end if	! calibration or sky data written
c
c  Free the logical unit number.
c
	   rstatus = fut_free_lun(cs_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal (%val(rstatus))
	   end if
	end if		! calibration or differential spectra written

	fsl_close_spectra = status

	return
	end
