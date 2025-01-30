	integer * 4 function  ffp_reformat_spectrum (in_rec, out_rec)

c-------------------------------------------------------------------------------
c	Function FFP_REFORMAT_SPECTRUM
c
c	This function reformats a FIRAS Pipeline spectrum into an FFP_SKY
c	spectrum.
c
c	Author:  S. Brodd, HSTX, 3/21/96
c
c	Input:
c		in_rec		FSL_SKY			input spectrum record
c
c	Output:
c		out_rec		FFP_SKY			output spectrum record
c
c	Subroutines called:
c		lib$movc5
c		xcc_e_to_g
c		xcc_q_to_e
c		UPX_Pixel_Number
c
c	Include files:
c		ffp_frequency.txt
c		ffp_invoc_sky.txt
c		fut_params.txt
c
c	Modifications:
c               S. Brodd, HSTX, 1/23/97, SPR 12339. Now delivering only
c               merged channel and scan mode destriped spectra.
c       
c               S. Brodd, HSTX, 2/27/97, SPR 12344. Do not take square root
c               of flag value -9999.0.
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(ffp_invoc_sky)'
	include '(ffp_frequency)'
	include '(uoe_constants)'

	integer * 2	res	!Pixel resolution for UPX_Pixel_Number (FIRAS=6)
	integer * 2	cmd_bias		!  commanded bolometer bias
	integer * 2	mrads			!  angles in 10e-4 radians

	integer * 4	j,nifg,tnifg		!  a counter; number of ifgs;
                                                !  number of template ifgs
					!  Attitude unit vectors in:
	real	* 4	equatorial(3)		!    equatorial coordinates
	real	* 4	ecliptic(3)		!    ecliptic coordinates
	real	* 4	galactic(3)		!    galactic coordinates

	real	* 8	time			!  collect time in seconds
						!    since T68 reference time
	dictionary 'fsl_sky'
	record /fsl_sky/ in_rec
	dictionary 'ffp_spc'
	record /ffp_spc/ out_rec

	real	* 4	rad_deg			! radians-to-degrees conversion
						!    function
	real	* 8	uoe_adt2t68

	Parameter	(res = 6)

	external	ffp_adttime
	external	ffp_normal
c
c  Define the radians to degrees conversion statement function.
c
	rad_deg(mrads) = floati(mrads) * 0.018 / fac_pi
c
c  Initialize the output spectrum record.
c
	call lib$movc5 (0,,0,4352,out_rec)
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
	enddo
c
c  The sigmas field.
c	Convert from ergs/sec/cm^2/sr/icm to MJy/sr.
c       Get the propagated spectrum variance from the FMS files by reading the
c	SPEC_DATA.REAL_IMAG_VAR field of the FSL_SKY RDL, which corresponds to
c	the SPEC_DATA.CVS_VAR field of the FMS_SKY RDL.
c
	if (fcc_data_type .ge. 2) then
	   do j = fcc_jlo,fcc_jhi
              if (in_rec.spec_data.real_var(j) .le. 0.0) then
	         out_rec.sigmas(j-fcc_freq_offset) =
     .		         in_rec.spec_data.real_var(j)
              else
	         out_rec.sigmas(j-fcc_freq_offset) =
     .	                 sqrt(in_rec.spec_data.real_var(j)) * fac_erg_to_mjy
              endif
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
        nifg = in_rec.coad_spec_head.num_ifgs
c
c  In the future, may need to do this for adj_num_ifgs (replace with
c  coad_spec_head.adj_comb_num_ifgs)
c
	if (fcc_num_type .eq. fac_present) then
	   out_rec.num_ifgs = in_rec.coad_spec_head.comb_num_ifgs
	else
	   out_rec.num_ifgs = floati(in_rec.coad_spec_head.num_ifgs)
	endif
	out_rec.adj_num_ifgs = in_rec.coad_spec_head.adj_num_ifgs
	time = uoe_adt2t68 (in_rec.ct_head.time)
	if (time .ge. 0.0D0) then
	   out_rec.time = time - T68IRASBase
	else
	   ffp_reformat_spectrum = %loc(ffp_adttime)
	   call lib$signal (ffp_adttime, %val(2), in_rec.ct_head.gmt,
     .					 %val(in_rec.attitude.pixel_no))
	   return
	endif
        do j = 1,nifg
           out_rec.ifg_times(j) = in_rec.coad_spec_head.times(j)
        enddo
        out_rec.orphans = in_rec.coad_spec_data.orphans
	out_rec.num_templates = in_rec.coad_spec_data.prim_template.subtracted +
     .				in_rec.coad_spec_data.sec_template.subtracted
        tnifg = in_rec.coad_spec_data.template.num_ifgs
        out_rec.tpl_num_ifgs = tnifg
        out_rec.neighbors = in_rec.coad_spec_data.template.neighbors
        out_rec.nbr_num_ifgs = 
     .                   in_rec.coad_spec_data.template.neighbor_num_ifgs
        do j = 1,tnifg
   	   out_rec.tpl_times(j) = 
     .             uoe_adt2t68(in_rec.coad_spec_data.template.ifgs(j).time) - 
     .             T68IRASBase
           out_rec.tpl_pixels(j) = 
     .             in_rec.coad_spec_data.template.ifgs(j).pixel_no
        enddo
c
c  The data and scan mode identification fields.
c
	out_rec.chanscan  = fcc_scan_mode
	out_rec.datatype  = fcc_data_type
        out_rec.fft_length = in_rec.spec_data.fft_length
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
	out_rec.bolom_temp    =  in_rec.coad_spec_data.bol_assem(fcc_fchan)
	out_rec.bath_temp     = (in_rec.en_analog.a_lo_bol_assem(fcc_fchan) +
     .				 in_rec.en_analog.b_lo_bol_assem(fcc_fchan))/2.0

	ffp_reformat_spectrum = %loc(ffp_normal)

	return
	end
