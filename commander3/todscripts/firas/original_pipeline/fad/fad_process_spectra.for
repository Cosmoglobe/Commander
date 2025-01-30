	Integer*4  Function  FAD_Process_Spectra  ( lun_in, lun_out, lun_rpt,
	1                                           report, fex_ejv_rec )

c------------------------------------------------------------------------------
c   Purpose: Process all the spectra.  Read spectra, apply correction,
c            and write spectra.
c
c   Input Parameters:
c      integer*4     lun_in        -- logical unit number for input file
c      integer*4     lun_out       -- logical unit number for output file
c      logical*1     report        -- true if report is to be written (default)
c      integer*4     lun_rpt       -- logical unit number for report file
c      record /fex_ejv/    fex_ejv_rec    -- EJV (FDS) reference record
c
c   Output Parameters:
c      none
c
c   Include Files:
c      FAD_msg       -- External message params needed for FAD
c      FUT_params    -- FIRAS global parameters
c      csa_pixel_input_rec
c      csa_pixel_output_rec
c
c   Functions and Subroutines:
c      CT_GMT_To_Binary
c      Lib$Signal
c      FAD_Apply_Correction
c      CSA_Read_Pixels, CSA_Write_Pixels
c      Lib$MovC3
c      Firas_Cenpix - return ecliptic unit vector of pixel center
c      Xcc_E_to_G - convert ecliptic unit vector to galactic unit vector
c
c   Author:  Larry Paris Rosen, Hughes STX, 20 April 1993
c   Modified: L. Rosen, May 1994.  Remove all references to gain correction.
c      FDS will not be generating a gain correction.
c   Modified: L. Rosen, June 1994.  Correction now includes vibration
c      correction polynomial (quartic) and several ufo terms with different
c      tau's.  Eject time is now in record as ADT.
c   Modified: L. Rosen, July 29, 1994.  Pass real variance to function
c      fad_apply_correction, so it can zero spectra where variance is zero.
c      SER 11847.
c   Modified: L. Rosen, August 10, 1994.  FAD now must check if each record's
c      pixel center is within or outside the galactic latitude cut used by
c      FDS.  FAD sets a new flag in the attitude record accordingly. SER 11869
c------------------------------------------------------------------------------
	Implicit None

c Include files

	Include		'(fad_msg)'
	Include		'(fut_params)'
	Include		'(csa_pixel_input_rec)'
	Include		'(csa_pixel_output_rec)'

c Passed Parameters:

	Integer*4	lun_in
	Integer*4	lun_out
	Logical*1	report
	Integer*4	lun_rpt
	Dictionary	'FEX_EJV'
	Record /FEX_EJV/ fex_ejv_rec

c External

	External	csa_normal

c Functions

	Integer*4	CSA_Read_Pixels, CSA_Write_Pixels
	Integer*4	FAD_Apply_Correction

c Local
	Integer*4	rstat
	Integer*2	npix			! Number of FIRAS pixels
	Integer*4	pix			! Pixel number
	Integer*2	numpix			! Number of pixels with data.
	Integer*4	numspec			! Number of spectra processed.
	Record /pixel_input_list/  inlist
	Record /pixel_output_list/ outlist
	Dictionary	'FCF_SKY'
	Record / FCF_SKY / fcf_rec (fac_max_skymap_recs) ! record read by CSA
	Dictionary	'FAD_SKY'
	Record / FAD_SKY / fad_rec (fac_max_skymap_recs) ! record written by CSA
	Integer*4	dummy				! Junk data
	Integer*4	rec			! record number (within pixel)
	Logical*1	first /.True./		! only check different model
						! name once.
	Real*4		glat		! galactic latitude of pixel in degrees
	Real*4		evec (3)		! unit vector in ecliptic
	Real*4		gvec (3)		! unit vector in galactic
c------------------------------------------------------------------------------

c Begin

	FAD_Process_Spectra = %Loc (FAD_Normal)

c Write FISH model name to report file (if report).

	If (report) Then
	   Write (lun_rpt, 10) fex_ejv_rec.Label
  10	   Format (1X,'FISH model solution name: ',A)
	Endif

c Pixel Loop: do for each FIRAS pixel.

	numpix = 0
	numspec = 0
	npix = ((4 ** fac_skymap_level) * 6) - 1
	inlist.level_no = fac_skymap_level
	pix = 0
	Do While ( (pix .LE. npix) .AND.
	1          (FAD_Process_Spectra .EQ. %Loc (FAD_Normal) ) )

	   inlist.pixel_no = pix

	   rstat = CSA_Read_Pixels (lun_in, inlist, 1, fcf_rec,
	1                           fac_max_skymap_recs, outlist, dummy, 0)

	   if ( rstat .NE. %Loc (csa_normal) ) Then
	      FAD_Process_Spectra = %Loc (FAD_Abort)
	      Call Lib$Signal (fad_csaread, %Val(1), %Val(rstat))
	   Else

	      If (outlist.no_records .GT. 0) Then
	         numpix = numpix + 1

c Record Loop: do for each input record.

	         rec = 1
	         Do While ( (rec .LE. outlist.no_records) .AND.
	1                   (FAD_Process_Spectra .EQ. %Loc (FAD_Normal)) )

c Check that FISH model labels match.  If not, give informational message.

	            If (first .AND. 
	1		(fcf_rec(rec).spec_data.model_label(6:) .NE.
	2                           fex_ejv_rec.Label(6:))) Then

	               Call Lib$Signal (fad_model_mismatch, %Val(2),
	1                               fex_ejv_rec.Label,
	2                               fcf_rec(rec).spec_data.model_label)
	               first = .False.
	            Endif

c Copy all record information to the fad record, and set destriped flag.

	            Call Lib$MovC3(fac_coad_spec_size,fcf_rec(rec),fad_rec(rec))

	            fad_rec (rec).spec_data.destriped = 1

c Call FAD_Apply_Correction

	            rstat = FAD_Apply_Correction
	1                      ( fad_rec (rec).ct_head.time,
	2                        fad_rec (rec).spec_data.spec,
	3                        fad_rec (rec).spec_data.real_var,
	4                        fex_ejv_rec )

	            If (rstat .NE. %Loc (FAD_Normal)) Then
	               FAD_Process_Spectra = %Loc (FAD_Abort)
	            Else

c Check if each records pixel center is within or outside the galactic
c latitude cut used by FDS.  FAD sets a new flag in the attitude record
c accordingly:      OUTSIDE_GALAXY_CUT is a byte
c                   0 = not set yet (pre-FAD data)
c                   1 = record included (outside galactic cut; abs(lat) > cut)
c                   2 = record excluded (within galactic cut; abs(lat) <= cut)

	               Call Firas_Cenpix ( pix, evec )
	               Call Xcc_E_to_G ( evec, fac_epoch, gvec )
	               glat = Atan2d ( gvec(3), sqrt (gvec(1)**2 + gvec(2)**2))
	               If ( Abs (glat) .LE. fex_ejv_rec.gallat ) Then
	                  fad_rec (rec).attitude.outside_galaxy_cut = 2
	               Else
	                  fad_rec (rec).attitude.outside_galaxy_cut = 1
	               Endif
	               numspec = numspec + 1
	               rec = rec + 1
	            Endif
	         Enddo                                     ! do for each record

c Write the FAD records for this pixel.

	         If (FAD_Process_Spectra .EQ. %Loc (FAD_Normal)) Then
	            rstat = CSA_Write_Pixels (lun_out, fad_rec,
	1                                     outlist.no_records, )

	            If ( rstat .NE. %Loc (CSA_Normal) ) Then
	               FAD_Process_Spectra = %Loc (FAD_Abort)
	               Call Lib$Signal (fad_csawrite, %Val(1), %Val(rstat))
	            Endif
	         Endif
	      Endif                                  ! if number of records > 0
	   Endif                               ! csa_read_pixels was successful
	   pix = pix + 1
	Enddo                                               ! do for each pixel

	If (report .AND. (FAD_Process_Spectra .EQ. %Loc (FAD_Normal))) Then
	   Write (lun_rpt, 20) numspec
  20	   Format (1X, 'Number of spectra processed: ', I5)
	   Write (lun_rpt, 30) numpix
  30	   Format (1X, 'Number of pixels with data: ', I5)
	Endif
	Return
	End
