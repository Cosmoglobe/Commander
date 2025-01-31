	integer * 4 function  ffp_reformat_ifg (in_rec, out_rec)

c------------------------------------------------------------------------------
c	Function FFP_REFORMAT_IFG
c	This function reformats a FIRAS Pipeline coadded ifg into an FFP_IFG
c	coadded ifg.
c
c	Author:  S. Brodd, HSTX, 3/21/96
c
c	Input:
c		in_rec		FIL_SKY			input coadd record
c
c	Output:
c		out_rec		FFP_IFG			output coadd record
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
c
c------------------------------------------------------------------------------

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
        integer * 4     adt_time(2)
					!  Attitude unit vectors in:
	real	* 4	equatorial(3)		!    equatorial coordinates
	real	* 4	ecliptic(3)		!    ecliptic coordinates
	real	* 4	galactic(3)		!    galactic coordinates

	real	* 8	time			!  collect time in seconds
						!    since T68 reference time
	dictionary 'fil_sky'
	record /fil_sky/ in_rec
	dictionary 'ffp_ifg'
	record /ffp_ifg/ out_rec

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
	call lib$movc5 (0,,0,6400,out_rec)
c
c  Convert the pointing unit vectors from equatorial coordinates to galactic
c	and ecliptic coordinates and fill the attitude data fields.
c
	Do j=1,3
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
c  The ifg and variance fields.
c	Convert the variance from (ergs/sec/cm^2/sr/icm_^2 to (MJy/sr)^2.
c
	do j = 1,512
	   out_rec.coadded_ifg(j) = in_rec.coad_data.ifg(j)
	enddo
	do j = fcc_jlo,fcc_jhi
	   out_rec.real_variance(j-fcc_freq_offset) =
     .			in_rec.coad_data.real_var(j) * fac_erg_to_mjy**2.0
	   out_rec.imag_variance(j-fcc_freq_offset) =
     .			in_rec.coad_data.imag_var(j) * fac_erg_to_mjy**2.0
	   out_rec.real_imag_variance(j-fcc_freq_offset) =
     .			in_rec.coad_data.real_imag_var(j) * fac_erg_to_mjy**2.0
	enddo
c
c  The coadd-specific data fields.
c
	out_rec.glitch_rate = in_rec.coad_spec_data.glitch_rate
        nifg = in_rec.coad_spec_head.num_ifgs
	out_rec.num_ifgs = nifg
	out_rec.adj_num_ifgs = in_rec.coad_spec_head.adj_num_ifgs
	time = uoe_adt2t68 (in_rec.ct_head.time)
	if (time .ge. 0.0D0) then
	   out_rec.time = time - T68IRASBase
	else
	   ffp_reformat_ifg = %loc(ffp_adttime)
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
	out_rec.chanscan      = fcc_scan_mode
        out_rec.datatype      = fcc_data_type
	out_rec.peak_position = in_rec.coad_spec_data.peak_pos
        out_rec.fft_length    = in_rec.coad_data.fft_length
	cmd_bias              = in_rec.coad_spec_data.bol_cmd_bias
	if (cmd_bias .lt. 0) cmd_bias = cmd_bias + 256  !  Two's complement
	out_rec.bolom_bias    = floati(cmd_bias)/25.5
	out_rec.bolom_voltage = in_rec.coad_spec_data.bol_volt
	out_rec.nu_zero       = fcc_nu0
	out_rec.delta_nu      = fcc_dnu
	out_rec.num_freq      = fcc_nfreq
c
c  The instrument temperature fields.
c
	out_rec.xcal_temp     =  in_rec.coad_spec_data.xcal
	out_rec.ical_temp     =  in_rec.coad_spec_data.ical
	out_rec.skyhorn_temp  =  in_rec.coad_spec_data.skyhorn
	out_rec.refhorn_temp  =  in_rec.coad_spec_data.refhorn
	out_rec.dihedral_temp =  in_rec.coad_spec_data.dihedral
	out_rec.mirror_temp   =  in_rec.coad_spec_data.collimator_mirror
	out_rec.bath_temp     = (in_rec.en_analog.a_lo_bol_assem(fcc_fchan) +
     .				 in_rec.en_analog.b_lo_bol_assem(fcc_fchan))/2.0

	ffp_reformat_ifg = %loc(ffp_normal)

	return
	end
