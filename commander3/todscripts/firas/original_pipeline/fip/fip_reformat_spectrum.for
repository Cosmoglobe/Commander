	integer * 4 function  fip_reformat_spectrum (in_rec, out_rec)

c-------------------------------------------------------------------------------
c	Function FIP_REFORMAT_SPECTRUM
c
c	This function reformats a FIRAS Pipeline spectrum into an FIP_SKY
c	spectrum.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  10 June 1993
c-------------------------------------------------------------------------------
c	Input:
c		in_rec		FCF_SKY			input spectrum record
c
c	Output:
c		out_rec		FIP_SKY			output spectrum record
c
c	Subroutines called:
c		lib$movc5
c		xcc_e_to_g
c		xcc_q_to_e
c		UPX_Pixel_Number
c
c	Include files:
c		fip_frequency.txt
c		fip_invoc_sky.txt
c		fut_params.txt
c-------------------------------------------------------------------------------
c	Changes:
c
c	Convert num_ifgs fields to floating point for use in merged scan modes.
c	Gene Eplee, GSC, 24 May 1994.
c
c	Convert FIP_SKY RDL for use with FCF spectrum records as well.
c	Gene Eplee, GSC, 15 July 1994.
c
c	Implement FCC_NUM_TYPE flag for the number of ifgs field.
c	Gene Eplee, GSC, 27 July 1994.
c
c	For merged skymaps (FMS skymaps), write the propagated spectrum
c	variance to the sigmas field instead of the real-to-real variance.
c       NOTE:  The field SPEC_DATA.REAL_IMAG_VAR in the FCF_SKY RDL corresponds
c	to the field SPEC_DATA.CVS_VAR in the FMS_SKY RDL.
c	Gene Eplee, GSC, 17 November 1994.  SPR 11989.
c
c	For calibration coadds (for which the sky fields really are irrelevant,
c	because, unlike sky coadds, the constituent IFGs can and do have
c	different pixel numbers), FIP was producing sky coordinates (equatorial,
c	ecliptic, galactic) which disagreed with pixel number.  Patch is to
c	compute all of the sky fields from the equatorial unit vector.
c	Fred Shuman, HSTX, 1995 Apr 10, 12174.
c
c       Steve Brodd, HSTX, 6/26/95, SPR 12208.  Remove the moon_angle and
c       galactic attitude fields; remove the destriped and calibrated flags.
c       
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fip_invoc_sky)'
	include '(fip_frequency)'
	include '(uoe_constants)'

	integer * 2	res	!Pixel resolution for UPX_Pixel_Number (FIRAS=6)

	integer * 2	cmd_bias		!  commanded bolometer bias
	integer * 2	mrads			!  angles in 10e-4 radians

	integer * 4	j			!  a counter

					!  Attitude unit vectors in:
	real	* 4	equatorial(3)		!    equatorial coordinates
	real	* 4	ecliptic(3)		!    ecliptic coordinates
	real	* 4	galactic(3)		!    galactic coordinates

	real	* 8	time			!  collect time in seconds
						!    since T68 reference time
	dictionary 'fcf_sky'
	record /fcf_sky/ in_rec
	dictionary 'fip_sky'
	record /fip_sky/ out_rec

	real	* 4	rad_deg			! radians-to-degrees conversion
						!    function
	real	* 8	uoe_adt2t68

	Parameter	(res = 6)

	external	fip_adttime
	external	fip_normal

C
C  Initialize the routine.
C

c
c  Define the radians to degrees conversion statement function.
c
	rad_deg(mrads) = floati(mrads) * 0.018 / fac_pi

c
c  Initialize the output spectrum record.
c
	call lib$movc5 (0,,0,2304,out_rec)

C
C  Fill in the primary data fields.
C

c
c  Convert the pointing unit vectors from equatorial coordinates to galactic
c	and ecliptic coordinates and fill the attitude data fields.
c
	Do j = 1,3
	   equatorial(j) = in_rec.attitude.equatorial(j)
	End Do
	call xcc_q_to_e (equatorial, fac_epoch, ecliptic)
	call xcc_e_to_g (ecliptic, fac_epoch, galactic)

	Call UPX_Pixel_Number(ecliptic, res, out_rec.pixel)

	If (equatorial(1) .Eq. 0. .And. equatorial(2) .Eq. 0.) Then
	   out_rec.ra = 0.
	Else
	   out_rec.ra = AMod(ATan2D(equatorial(2), equatorial(1)) + 360., 360.)
	Endif
	out_rec.dec   = ASinD(equatorial(3))

	If (ecliptic(1) .Eq. 0. .And. ecliptic(2) .Eq. 0.) Then
	   out_rec.eclon = 0.
	Else
	   out_rec.eclon = AMod(ATan2D(ecliptic(2), ecliptic(1)) + 360., 360.)
	Endif
	out_rec.eclat = ASinD(ecliptic(3))

	If (galactic(1) .Eq. 0. .And. galactic(2) .Eq. 0.) Then
	   out_rec.galon = 0.
	Else
	   out_rec.galon = AMod(ATan2D(galactic(2), galactic(1)) + 360., 360.)
	Endif
	out_rec.galat = ASinD(galactic(3))

	out_rec.galatexc   = fcc_glat

c
c  The spectrum field.
c	Convert from ergs/sec/cm^2/sr/icm to MJy/sr.
c
	do j = fcc_jlo,fcc_jhi
	   out_rec.real_spectrum(j-fcc_freq_offset) =
     .			real(in_rec.spec_data.spec(j)) * fac_erg_to_mjy
	   out_rec.imag_spectrum(j-fcc_freq_offset) =
     .			aimag(in_rec.spec_data.spec(j)) * fac_erg_to_mjy
	   out_rec.sigmas(j-fcc_freq_offset) =
     .			sqrt(in_rec.spec_data.real_var(j)) * fac_erg_to_mjy
	enddo

c
c  The sigmas field.
c	Convert from ergs/sec/cm^2/sr/icm to MJy/sr.
c       Get the propagated spectrum variance from the FMS files by reading the
c	SPEC_DATA.REAL_IMAG_VAR field of the FCF_SKY RDL, which corresponds to
c	the SPEC_DATA.CVS_VAR field of the FMS_SKY RDL.
c
	if (fcc_data_type .ge. 3) then
	   do j = fcc_jlo,fcc_jhi
	      out_rec.sigmas(j-fcc_freq_offset) =
     .			   sqrt(in_rec.spec_data.real_var(j)) * fac_erg_to_mjy
	   enddo
	else
	   do j = fcc_jlo,fcc_jhi
	      out_rec.sigmas(j-fcc_freq_offset) =
     .		sqrt(in_rec.spec_data.real_imag_var(j)) * fac_erg_to_mjy
	   enddo
	endif

c
c  The spectrum-specific data fields.
c
	out_rec.glitch_rate = in_rec.coad_spec_data.glitch_rate
	if (fcc_num_type .eq. fac_present) then
	   out_rec.num_ifgs = in_rec.coad_spec_head.comb_num_ifgs
	else
	   out_rec.num_ifgs = floati(in_rec.coad_spec_head.num_ifgs)
	endif
	out_rec.num_templates = in_rec.coad_spec_data.prim_template.subtracted +
     .				in_rec.coad_spec_data.sec_template.subtracted
	time = uoe_adt2t68 (in_rec.ct_head.time)
	if (time .ge. 0.0D0) then
	   out_rec.time = time - T68IRASBase
	else
	   fip_reformat_spectrum = %loc(fip_adttime)
	   call lib$signal (fip_adttime, %val(2), in_rec.ct_head.gmt,
     .					 %val(in_rec.attitude.pixel_no))
	   return
	endif


C
C  Fill in the secondary data fields.
C

c
c  The data and scan mode identification fields.
c
	out_rec.chanscan  = fcc_scan_mode
	out_rec.datatype  = fcc_data_type
	out_rec.nu_zero   = fcc_nu0
	out_rec.delta_nu  = fcc_dnu
	out_rec.num_freq  = fcc_nfreq

c
c  The bolometer response function fields.
c
	cmd_bias = in_rec.coad_spec_data.bol_cmd_bias
	if (cmd_bias .lt. 0) cmd_bias = cmd_bias + 256  !  Two's complement
	out_rec.bolom_bias    = floati(cmd_bias)/25.5
	out_rec.bolom_voltage = in_rec.coad_spec_data.bol_volt
	out_rec.dc_response   = in_rec.spec_data.responsivity
	out_rec.time_constant = in_rec.spec_data.time_constant
	out_rec.phase_corr    = in_rec.spec_data.phase_corr

c
c  The instrument temperature fields.
c
	out_rec.xcal_temp     =  in_rec.coad_spec_data.xcal
	out_rec.ical_temp     =  in_rec.coad_spec_data.ical
	out_rec.skyhorn_temp  =  in_rec.coad_spec_data.skyhorn
	out_rec.refhorn_temp  =  in_rec.coad_spec_data.refhorn
	out_rec.dihedral_temp =  in_rec.coad_spec_data.dihedral
	out_rec.mirror_temp   =  in_rec.coad_spec_data.collimator_mirror
	out_rec.bolom_temp    =  in_rec.coad_spec_data.bol_assem(1)
	out_rec.bath_temp     = (in_rec.en_analog.a_lo_bol_assem(1) +
     .				 in_rec.en_analog.b_lo_bol_assem(1))/2.0

	fip_reformat_spectrum = %loc(fip_normal)

	return
	end
