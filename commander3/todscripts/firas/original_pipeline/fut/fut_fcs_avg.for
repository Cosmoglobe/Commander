	Integer*4  Function  FUT_FCS_Avg  ( maxsum, time_array, numrecs,
	1                                   time_weights, out_rec, sum, rec1,
	2                                   sum_weight, mphase, ophase, sangle,
	3                                   total )

c-----------------------------------------------------------------------------
c	Purpose:  To compute final averages of all values, except spectra and
c          variances, during sky map combination.  Used by FCS and FMS.
c	Author: Larry P. Rosen, Hughes STX, October 1993.
c
c	Modified: LPR, HSTX, August 1994, SPR 11878
c	Certain en_analog fields need to be checked for flag values instead of
c	just blindly averaging.  FUT_FCS_SUM puts -9999 in the rec1 if a flag
c	was found.
c-----------------------------------------------------------------------------

	Implicit None

c Include.

	Include		'(fut_fcs_include)'
	Include		'(fut_params)'

c Passed parameters.

	Integer*4	maxsum			! max sum of nifgs per record
	Integer*4	time_array (2, maxsum)	! Time array to average
	Integer*2	numrecs			! Number of records to average
	Integer*4	time_weights (maxsum)	! Time average weights = nifgs
	Dictionary	'FCS_SKY'		! Same as record type for FMS
	Record / fcs_sky / out_rec		! Output averaged map record
	Record / double_fcf / sum		! Record holding sums.
	Record / double_fcf / rec1		! Record storing 1st data points
	Real*8		sum_weight		! Sum of weights for average
	Real*4		mphase (maxsum)		! Moon phase array
	Real*4		ophase (maxsum)		! Orbital phase array
	Real*4		sangle (maxsum)		! Sun angle array
	Integer*4	total			! Total # of IFGs combined

c Functions

	Integer*4	FUT_Ave_Bin_Times
	Integer*4	FUT_Average_Angles

c Local

	Integer*4	rstat				! return status
	Integer*2	i, j				! counter
	Real*8		variance
	Real*8		mean				! r*8 for calculations
	Real*8 		norm_equat, norm_terr, norm_elimb
	Real*8		norm_mnang, norm_gal, norm_ecl
	Real*8		avg_ecl(3), avg_gal(3), avg_elimb(3)
	Real*8		avg_mnang(3), avg_terr(3)
	Real*8		ra, dec, terr_lat, terr_long, elimb, el_az
	Real*8		moon_angle, moon_az, gal_lat, gal_long
	Real*8		ecl_lat, ecl_long
	Real*4		avg_angle		! real*4 angle from fut_avg_ang
	Real*4		rterr(3)		! real*4 version of avg_terr

c-----------------------------------------------------------------------------
c Begin

	FUT_FCS_Avg = 0

	rstat = FUT_Ave_Bin_Times ( time_array, numrecs, time_weights,
	1                           out_rec.ct_head.time )

	Call CT_Binary_To_GMT ( out_rec.ct_head.time, out_rec.ct_head.gmt )

c coad_spec_data

	Do i = 1, 10
	   mean = sum.coad_spec_data.temp (i) / sum_weight +
	1            rec1.coad_spec_data.temp (i)
	   out_rec.coad_spec_data.temp (i) = Sngl ( mean )
	   variance = (sum.coad_spec_data.temp_sigma (i) - sum_weight *
	1     (mean - rec1.coad_spec_data.temp (i))**2 ) / (sum_weight - 1)
	   If (variance .LE. 0) Then
	      out_rec.coad_spec_data.temp_sigma (i) = 1.E15
	   Else
	      out_rec.coad_spec_data.temp_sigma (i) = Sngl (Sqrt (variance))
	   Endif
	Enddo

	out_rec.coad_spec_data.glitch_rate =
	1   Sngl ( sum.coad_spec_data.glitch_rate / sum_weight +
	2          rec1.coad_spec_data.glitch_rate )
	out_rec.coad_spec_data.adds_per_group =
	1   Nint ( sum.coad_spec_data.adds_per_group / sum_weight +
	2   rec1.coad_spec_data.adds_per_group )
	out_rec.coad_spec_data.sweeps = 
	1   Nint ( sum.coad_spec_data.sweeps / sum_weight +
	2   rec1.coad_spec_data.sweeps )
	out_rec.coad_spec_data.peak_pos =
	1   Nint ( sum.coad_spec_data.peak_pos / sum_weight +
	2   rec1.coad_spec_data.peak_pos )
	out_rec.coad_spec_data.nyquist_hertz =
	1   Sngl ( sum.coad_spec_data.nyquist_hertz / sum_weight +
	2          rec1.coad_spec_data.nyquist_hertz )
	out_rec.coad_spec_data.nyquist_icm =
	1   Sngl ( sum.coad_spec_data.nyquist_icm / sum_weight +
	2          rec1.coad_spec_data.nyquist_icm )
	out_rec.coad_spec_data.prim_template.amplitude =
	1   Sngl ( sum.coad_spec_data.prim_template.amplitude / sum_weight +
	2          rec1.coad_spec_data.prim_template.amplitude )
	out_rec.coad_spec_data.prim_template.snr =
	1   Sngl ( sum.coad_spec_data.prim_template.snr / sum_weight +
	2          rec1.coad_spec_data.prim_template.snr )
	out_rec.coad_spec_data.sec_template.amplitude =
	1   Sngl ( sum.coad_spec_data.sec_template.amplitude / sum_weight +
	2          rec1.coad_spec_data.sec_template.amplitude )
	out_rec.coad_spec_data.sec_template.snr =
	1   Sngl ( sum.coad_spec_data.sec_template.snr / sum_weight +
	2          rec1.coad_spec_data.sec_template.snr )
	out_rec.coad_spec_data.sec_template.variance =
	1   Sngl ( ( 511. * sum.coad_spec_data.sec_template.variance ) /
	2          ( (512. * sum_weight) - 1 ) )

	mean = sum.coad_spec_data.sec_template.b_average / sum_weight +
	1          rec1.coad_spec_data.sec_template.b_average
	out_rec.coad_spec_data.sec_template.b_average = Sngl ( mean )
	out_rec.coad_spec_data.sec_template.b_variance = Sngl (
	1  (sum.coad_spec_data.sec_template.b_variance - sum_weight *
	2   (mean - rec1.coad_spec_data.sec_template.b_average)**2 ) /
	3   (sum_weight-1) )

	out_rec.coad_spec_data.noise =
	1   Sngl ( sum.coad_spec_data.noise / sum_weight +
	2          rec1.coad_spec_data.noise )
	Do i=1,5
	   out_rec.coad_spec_data.bl_coeffs (i) =
	1     Sngl ( sum.coad_spec_data.bl_coeffs (i) / sum_weight +
	2            rec1.coad_spec_data.bl_coeffs (i) )
	Enddo

	mean = sum.coad_spec_data.transient.c_average / sum_weight +
	1          rec1.coad_spec_data.transient.c_average
	out_rec.coad_spec_data.transient.c_average = Sngl ( mean )
	out_rec.coad_spec_data.transient.c_variance = Sngl (
	1  (sum.coad_spec_data.transient.c_variance - sum_weight *
	2   (mean - rec1.coad_spec_data.transient.c_average)**2 ) /
	3   (sum_weight-1) )

	out_rec.coad_spec_data.transient.bl_trans_coeff =
	1   Sngl ( sum.coad_spec_data.transient.bl_trans_coeff / sum_weight +
	2          rec1.coad_spec_data.transient.bl_trans_coeff )
	out_rec.coad_spec_data.deglitch.glitch_iter =
	1   Nint ( sum.coad_spec_data.deglitch.glitch_iter / sum_weight ) +
	2   rec1.coad_spec_data.deglitch.glitch_iter
	out_rec.coad_spec_data.deglitch.glitch_signal =
	1   Sngl ( sum.coad_spec_data.deglitch.glitch_signal / sum_weight +
	2          rec1.coad_spec_data.deglitch.glitch_signal )
	out_rec.coad_spec_data.bol_cmd_bias =
	1   Nint ( sum.coad_spec_data.bol_cmd_bias / sum_weight +
	2   rec1.coad_spec_data.bol_cmd_bias )
	out_rec.coad_spec_data.bol_volt =
	1   Sngl ( sum.coad_spec_data.bol_volt / sum_weight +
	2          rec1.coad_spec_data.bol_volt )

C spec_data

	mean = sum.spec_data.responsivity / sum_weight +
	1         rec1.spec_data.responsivity
	out_rec.spec_data.responsivity = Sngl ( mean )
	variance = (sum.spec_data.resp_sigma - sum_weight *
	1  (mean - rec1.spec_data.responsivity)**2) / (sum_weight-1)
	If (variance .LE. 0) Then
	   out_rec.spec_data.resp_sigma = 1.E15
	Else
	   out_rec.spec_data.resp_sigma = Sngl (Sqrt (variance))
	Endif

	mean = sum.spec_data.time_constant / sum_weight +
	1          rec1.spec_data.time_constant
	out_rec.spec_data.time_constant = Sngl ( mean )
	variance = (sum.spec_data.tc_sigma - sum_weight *
	1  (mean - rec1.spec_data.time_constant)**2 ) / (sum_weight-1)
	If (variance .LE. 0) Then
	   out_rec.spec_data.tc_sigma = 1.E15
	Else
	   out_rec.spec_data.tc_sigma = Sngl (Sqrt (variance))
	Endif

	mean = sum.spec_data.phase_corr / sum_weight +
	1          rec1.spec_data.phase_corr
	out_rec.spec_data.phase_corr = Sngl ( mean )
	variance = ( sum.spec_data.phase_corr_sqrd - sum_weight *
	1   (mean - rec1.spec_data.phase_corr)**2 ) / (sum_weight-1)
	If (variance .LE. 0) Then
	   out_rec.spec_data.pc_sigma = 1.E15         ! flag value
	Else
	   out_rec.spec_data.pc_sigma = Sngl (Sqrt (variance))
	Endif

	mean = sum.spec_data.qrad / sum_weight + rec1.spec_data.qrad
	out_rec.spec_data.qrad = Sngl ( mean )
	variance = ( sum.spec_data.qrad_sqrd - sum_weight *
	1   (mean - rec1.spec_data.qrad)**2 ) / (sum_weight-1)
	If (variance .LE. 0) Then
	   out_rec.spec_data.qrad_sigma = 1.E15         ! flag value
	Else
	   out_rec.spec_data.qrad_sigma = Sngl (Sqrt (variance))
	Endif

	mean = sum.spec_data.ir_power / sum_weight + rec1.spec_data.ir_power
	out_rec.spec_data.ir_power = Sngl ( mean )
	variance = ( sum.spec_data.ir_power_sqrd - sum_weight *
	1   (mean - rec1.spec_data.ir_power)**2 ) / (sum_weight-1)
	If (variance .LE. 0) Then
	   out_rec.spec_data.ir_power_sigma = 1.E15        ! flag value
	Else
	   out_rec.spec_data.ir_power_sigma = Sngl (Sqrt (variance))
	Endif

c eng_status

	Do i = 1, 14
	   out_rec.en_stat.group1 (i) =
	1     Nint ( sum.en_stat.group1 (i) / sum_weight +
	2            rec1.en_stat.group1 (i) )
	   out_rec.en_stat.group2 (i) =
	1     Nint ( sum.en_stat.group2 (i) / sum_weight +
	2            rec1.en_stat.group2 (i) )
	Enddo
	out_rec.en_stat.group1(15) =
	1   Nint ( sum.en_stat.group1(15) / sum_weight +
	2          rec1.en_stat.group1(15) )
	out_rec.en_stat.group1(16) =
	1   Nint ( sum.en_stat.group1(16) / sum_weight +
	2          rec1.en_stat.group1(16) )
	out_rec.en_stat.hot_spot_cmd(1) =
	1   Nint ( sum.en_stat.hot_spot_cmd(1) / sum_weight +
	2          rec1.en_stat.hot_spot_cmd(1) )
	out_rec.en_stat.hot_spot_cmd(2) =
	1   Nint ( sum.en_stat.hot_spot_cmd(2) / sum_weight +
	2          rec1.en_stat.hot_spot_cmd(2) )
	out_rec.en_stat.power_a_status(1) =
	1   Nint ( sum.en_stat.power_a_status(1) / sum_weight +
	2          rec1.en_stat.power_a_status(1) )
	out_rec.en_stat.power_a_status(2) =
	1   Nint ( sum.en_stat.power_a_status(2) / sum_weight +
	2          rec1.en_stat.power_a_status(2) )
	out_rec.en_stat.power_b_status(1) =
	1   Nint ( sum.en_stat.power_b_status(1) / sum_weight +
	2          rec1.en_stat.power_b_status(1) )
	out_rec.en_stat.power_b_status(2) =
	1   Nint ( sum.en_stat.power_b_status(2) / sum_weight +
	2          rec1.en_stat.power_b_status(2) )

c eng_analog & eng_sigma

	Do i = 1,62
	   If (rec1.en_analog.group1 (i) .EQ. -9999.) Then
	      out_rec.en_analog.group1 (i) = -9999.
	      variance = -9999.
	   Else
	      mean = sum.en_analog.group1 (i) / sum_weight +
	1               rec1.en_analog.group1 (i)
	      out_rec.en_analog.group1 (i) = Sngl ( mean )
	      variance = (sum.en_sigma.group2 (i) - sum_weight *
	1        (mean - rec1.en_analog.group1 (i))**2 ) / (sum_weight-1)
	   EndIf
	   If (variance .LE. 0) Then
	      out_rec.en_sigma.group2 (i) = 1.E15
	   Else
	      out_rec.en_sigma.group2 (i) = Sngl (Sqrt (variance))
	   Endif

	   mean = sum.en_analog.grt (i) / sum_weight +
	1            rec1.en_analog.grt (i)
	   out_rec.en_analog.grt (i) = Sngl ( mean )
	   variance = (sum.en_sigma.sig_grt (i) - sum_weight *
	1     (mean - rec1.en_analog.grt (i))**2) / (sum_weight-1)
	   If (variance .LE. 0) Then
	      out_rec.en_sigma.sig_grt (i) = 1.E15
	   Else
	      out_rec.en_sigma.sig_grt (i) = Sngl (Sqrt (variance))
	   Endif
	Enddo
	Do i = 63,64
	   mean = sum.en_analog.grt (i) / sum_weight + rec1.en_analog.grt (i)
	   out_rec.en_analog.grt (i) = Sngl ( mean )
	   variance = (sum.en_sigma.sig_grt (i) - sum_weight *
	1     (mean - rec1.en_analog.grt (i))**2) / (sum_weight-1)
	   If (variance .LE. 0) Then
	      out_rec.en_sigma.sig_grt (i) = 1.E15
	   Else
	      out_rec.en_sigma.sig_grt (i) = Sngl (Sqrt (variance))
	   Endif
	Enddo

c eng_tempdiff

	Do i = 1,2
	   out_rec.en_tempdiff(i).xcal =
	1     Nint ( sum.en_tempdiff(i).xcal / sum_weight +
	2            rec1.en_tempdiff(i).xcal )
	   out_rec.en_tempdiff(i).ical =
	1     Nint ( sum.en_tempdiff(i).ical / sum_weight +
	2            rec1.en_tempdiff(i).ical )
	   out_rec.en_tempdiff(i).skyhorn =
	1     Nint ( sum.en_tempdiff(i).skyhorn / sum_weight +
	2            rec1.en_tempdiff(i).skyhorn )
	   out_rec.en_tempdiff(i).refhorn =
	1     Nint ( sum.en_tempdiff(i).refhorn / sum_weight +
	2            rec1.en_tempdiff(i).refhorn )
	   out_rec.en_tempdiff(i).dihedral =
	1     Nint ( sum.en_tempdiff(i).dihedral / sum_weight +
	2            rec1.en_tempdiff(i).dihedral )
	   out_rec.en_tempdiff(i).collimator_mirror =
	1     Nint ( sum.en_tempdiff(i).collimator_mirror / sum_weight +
	2            rec1.en_tempdiff(i).collimator_mirror )
	   Do j = 1,4
	      out_rec.en_tempdiff(i).bol_assem(j) =
	1        Nint ( sum.en_tempdiff(i).bol_assem(j) / sum_weight +
	2               rec1.en_tempdiff(i).bol_assem(j) )
	   Enddo
	Enddo

C attitude

	out_rec.attitude.sun_moon_dist =
	1   Sngl ( sum.attitude.sun_moon_dist / sum_weight +
	2          rec1.attitude.sun_moon_dist )
	out_rec.attitude.cobe_moon_dist =
	1   Sngl ( sum.attitude.cobe_moon_dist / sum_weight +
	2          rec1.attitude.cobe_moon_dist )
	out_rec.attitude.altitude =
	1   Nint ( sum.attitude.altitude / sum_weight +
	2   rec1.attitude.altitude )
	out_rec.attitude.projected_barycentric_velocity =
	1   Nint ( sum.attitude.projected_barycentric_velocity / sum_weight +
	2   rec1.attitude.projected_barycentric_velocity )
	out_rec.attitude.mcilwain_l_param =
	1   Nint ( sum.attitude.mcilwain_l_param / sum_weight +
	2   rec1.attitude.mcilwain_l_param )
	norm_equat = Sqrt ( sum.attitude.equat(1)**2 +
	1                   sum.attitude.equat(2)**2 +
	2                   sum.attitude.equat(3)**2 )
	norm_terr  = Sqrt ( sum.attitude.terr(1)**2 +
	1                   sum.attitude.terr(2)**2 +
	2                   sum.attitude.terr(3)**2 )
	norm_elimb = Sqrt ( sum.attitude.elimb(1)**2 +
	1                   sum.attitude.elimb(2)**2 +
	2                   sum.attitude.elimb(3)**2 )
	norm_mnang = Sqrt ( sum.attitude.mnang(1)**2 +
	1                   sum.attitude.mnang(2)**2 +
	2                   sum.attitude.mnang(3)**2 )
	norm_gal   = Sqrt ( sum.attitude.gal(1)**2 +
	1                   sum.attitude.gal(2)**2 +
	2                   sum.attitude.gal(3)**2 )
	norm_ecl   = Sqrt ( sum.attitude.ecl(1)**2 +
	1                   sum.attitude.ecl(2)**2 +
	2                   sum.attitude.ecl(3)**2 )
	Do i = 1,3
	   out_rec.attitude.equatorial(i) =
	1     Sngl ( sum.attitude.equat(i) / norm_equat )
	   avg_ecl(i) = sum.attitude.ecl(i) / norm_ecl
	   avg_gal(i) = sum.attitude.gal(i) / norm_gal
	   avg_elimb(i) = sum.attitude.elimb(i) / norm_elimb
	   avg_mnang(i) = sum.attitude.mnang(i) / norm_mnang
	   avg_terr(i) = sum.attitude.terr(i) / norm_terr
	EndDo
	ra = 0.0
	If (out_rec.attitude.equatorial(1) .NE. 0.0) Then
	   ra = Atan2 ( (sum.attitude.equat(2) / norm_equat),
	1               (sum.attitude.equat(1) / norm_equat) )
	Endif
	out_rec.attitude.ra = Nint (ra * 10000.0)
	dec = asin (sum.attitude.equat(3) / norm_equat)
	out_rec.attitude.dec = Nint (dec * 10000.0)
	terr_lat = Atan2d ( avg_terr(3),
	1                   Sqrt (avg_terr(1)**2 + avg_terr(2)**2) )
	terr_long = 0.0
	If (avg_terr(1) .NE. 0.0) Then
	    terr_long = Atan2d (avg_terr(2), avg_terr(1))
	Endif
	out_rec.attitude.terr_latitude = Nint (terr_lat / fac_att_conv)
	out_rec.attitude.terr_longitude = Nint (terr_long / fac_att_conv)

	Do i=1,3
	   rterr (i) = SNGL ( avg_terr (i) )
	Enddo
	Call Firas_pixno ( rterr, out_rec.attitude.terr_pixel_no )

	elimb = atan2d ( Sqrt (avg_elimb(1)**2 + avg_elimb(2)**2),
	1                avg_elimb(3) )
	el_az = 0.0
	If (avg_elimb(1) .NE. 0.0) Then
	   el_az = Atan2d (avg_elimb(2), avg_elimb(1))
	Endif
	out_rec.attitude.earth_limb = Nint ((elimb + 90.0) / fac_att_conv)
	out_rec.attitude.earth_limb_azimuth = Nint (el_az / fac_att_conv)
	moon_angle = Atan2d ( Sqrt (avg_mnang(1)**2 + avg_mnang(2)**2),
	1                     avg_mnang(3) )
	moon_az = 0.0
	If (avg_mnang(1) .NE. 0.0) Then
	    moon_az = Atan2d (avg_mnang(2), avg_mnang(1))
	Endif
	out_rec.attitude.moon_angle = Nint ((moon_angle + 90.0)/fac_att_conv)
	out_rec.attitude.moon_az_angle = Nint (moon_az / fac_att_conv)
	gal_lat = Atan2d (avg_gal(3), Sqrt (avg_gal(1)**2 + avg_gal(2)**2))
	gal_long = 0.0
	If (avg_gal(1) .NE. 0.0) Then
	    gal_long = atan2d (avg_gal(2), avg_gal(1))
	Endif
	out_rec.attitude.galactic_latitude = Nint (gal_lat / fac_att_conv)
	out_rec.attitude.galactic_longitude = Nint (gal_long / fac_att_conv)
	ecl_lat = Atan2d (avg_ecl(3), Sqrt (avg_ecl(1)**2 + avg_ecl(2)**2))
	ecl_long = 0.0
	If (avg_ecl(1) .NE. 0.0) Then
	   ecl_long = Atan2d (avg_ecl(2), avg_ecl(1))
	Endif
	out_rec.attitude.ecliptic_latitude = Nint (ecl_lat / fac_att_conv)
	out_rec.attitude.ecliptic_longitude = Nint (ecl_long / fac_att_conv)
	out_rec.attitude.sun_angle =
	1  Nint ( sum.attitude.sun_angle / sum_weight )
	out_rec.attitude.projected_geocentric_velocity =
	1  Nint ( sum.attitude.geovel / sum_weight )
	out_rec.attitude.sc_rotation_angle =
	1  Nint ( sum.attitude.sc_rot_angle / sum_weight )
	rstat = FUT_Average_Angles (mphase, total, avg_angle)
	out_rec.attitude.moon_phase = Nint (avg_angle / fac_att_conv_rad)
	rstat = FUT_Average_Angles (ophase, total, avg_angle)
	out_rec.attitude.orbital_phase = Nint (avg_angle / fac_att_conv_rad)
	rstat = FUT_Average_Angles (sangle, total, avg_angle)
	out_rec.attitude.scan_angle = Nint (avg_angle / fac_att_conv_rad)
	out_rec.attitude.exc_galactic_lat =
	1  Nint (sum.attitude.exc_gal_lat / sum_weight)
	Return
	End
