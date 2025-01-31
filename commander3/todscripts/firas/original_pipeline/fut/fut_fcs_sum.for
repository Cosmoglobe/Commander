	Integer*4  Function  FUT_FCS_Sum  ( weight, out_rec, in_rec,
	1                                   moon_phase, orb_phase, scan_angle,
	2                                   sum, rec1, first )
c-----------------------------------------------------------------------------
c	Purpose: Store sum of values from FCF, FCS, or FMS record in "sum"
c	   record, to compute averages and sigmas.
c	Author: Larry P. Rosen, Hughes STX, March-April 1992.
c	This uses the trick of subtracting the first datum from each of the
c	data to be summed, which is then added to the computed average.  This
c	reduces errors due to round off.  The first datum record is stored in
c	"rec1".  See requirements and design for details.
c	Modification: LPR, October 1993.
c	   Routine is now an FUT routine.  Weighting is now all done with
c	   whatever is passed in as "weight".  When run in FCS the weight is
c	   the number of ifgs.  When run in FMS, the weight is different.
c      Modification: LPR, HSTX, August 1994  SPR 11878
c          Sweeps are not used anymore to weight b-variance or other
c          combinations.  The glitch rate weighted number of ifgs is fine
c          since it gets the d-vector ratio weight during zeroth order
c          combination.
c      Modification: LPR, HSTX, August 1994  SPR
c          Some en_analog fields need to be checked for flag values rather
c          than blindly averaging.
c-----------------------------------------------------------------------------
	Implicit None

	Include		'(fut_params)'
	Include		'(fut_fcs_include)'

c Passed parameters.

	Real*8		weight
	Dictionary	'FCF_SKY'
	Record / fcf_sky / in_rec
	Dictionary	'FCS_SKY'
	Record / fcs_sky / out_rec
	Real*4		moon_phase, orb_phase, scan_angle
	Record / double_fcf / sum			! Record to store sums.
	Record / double_fcf / rec1	! Record to store first data record.
	Logical*1	first

c Local

	Integer*4	m, n			! counters
	Real*8		terr_lat, terr_long, elimb, el_az, moon_angle, moon_az
	Real*8		gal_lat, gal_long, ecl_lat, ecl_long
	Real*8		wm1		! weight - 1
c-----------------------------------------------------------------------------
c Begin

	FUT_FCS_Sum = 0
	wm1 = weight - 1.
	If (first) Then				! store 1st data in rec1 record
c------------------------------------------------------------------------------

c coad_spec_data

	   out_rec.coad_spec_data.gain_sum =
	1     DBLE (in_rec.coad_spec_data.gain_sum)

	   Do m = 1, 10
	      rec1.coad_spec_data.temp(m) =
	1        DBLE (in_rec.coad_spec_data.temp(m))
	   Enddo
	   rec1.coad_spec_data.glitch_rate =
	1     DBLE (in_rec.coad_spec_data.glitch_rate)
	   rec1.coad_spec_data.adds_per_group =
	1     DBLE (Zext (in_rec.coad_spec_data.adds_per_group))
	   rec1.coad_spec_data.sweeps =
	1     DBLE (Zext (in_rec.coad_spec_data.sweeps))
	   rec1.coad_spec_data.peak_pos =
	1     DBLE (in_rec.coad_spec_data.peak_pos)
	   rec1.coad_spec_data.nyquist_hertz =
	1     DBLE (in_rec.coad_spec_data.nyquist_hertz)
	   rec1.coad_spec_data.nyquist_icm =
	1     DBLE (in_rec.coad_spec_data.nyquist_icm)
	   rec1.coad_spec_data.prim_template.amplitude =
	1     DBLE (in_rec.coad_spec_data.prim_template.amplitude)
	   rec1.coad_spec_data.prim_template.snr =
	1     DBLE (in_rec.coad_spec_data.prim_template.snr)
	   rec1.coad_spec_data.sec_template.amplitude =
	1     DBLE (in_rec.coad_spec_data.sec_template.amplitude)
	   sum.coad_spec_data.sec_template.variance =
	1     DBLE (in_rec.coad_spec_data.sec_template.variance * weight)
	   rec1.coad_spec_data.sec_template.snr =
	1     DBLE (in_rec.coad_spec_data.sec_template.snr)
	   rec1.coad_spec_data.sec_template.b_average =
	1     DBLE (in_rec.coad_spec_data.sec_template.b_average)
	   rec1.coad_spec_data.transient.c_average =
	1     DBLE (in_rec.coad_spec_data.transient.c_average)
	   rec1.coad_spec_data.transient.bl_trans_coeff =
	1     DBLE (in_rec.coad_spec_data.transient.bl_trans_coeff)
	   rec1.coad_spec_data.deglitch.glitch_iter =
	1     DBLE (in_rec.coad_spec_data.deglitch.glitch_iter)
	   rec1.coad_spec_data.deglitch.glitch_signal =
	1     DBLE (in_rec.coad_spec_data.deglitch.glitch_signal)
	   rec1.coad_spec_data.noise = DBLE (in_rec.coad_spec_data.noise)
	   Do m=1,5
	      rec1.coad_spec_data.bl_coeffs(m) =
	1        DBLE (in_rec.coad_spec_data.bl_coeffs(m))
	   Enddo
	   rec1.coad_spec_data.bol_cmd_bias =
	1     DBLE (in_rec.coad_spec_data.bol_cmd_bias)
	   rec1.coad_spec_data.bol_volt =
	1     DBLE (in_rec.coad_spec_data.bol_volt)

c spec_data

	   rec1.spec_data.responsivity = DBLE (in_rec.spec_data.responsivity)
	   rec1.spec_data.time_constant = DBLE (in_rec.spec_data.time_constant)
	   rec1.spec_data.phase_corr = DBLE (in_rec.spec_data.phase_corr)
	   rec1.spec_data.qrad = DBLE (in_rec.spec_data.qrad)
	   rec1.spec_data.ir_power = DBLE (in_rec.spec_data.ir_power)

c eng_status

	   Do m = 1,14
	      rec1.en_stat.group1(m) = DBLE (in_rec.en_stat.group1(m))
	      rec1.en_stat.group2(m) = DBLE (Zext (in_rec.en_stat.group2(m)))
	   Enddo
	   rec1.en_stat.group1(15) = DBLE (in_rec.en_stat.group1(15))
	   rec1.en_stat.group1(16) = DBLE (in_rec.en_stat.group1(16))
	   rec1.en_stat.hot_spot_cmd(1) =
	1     DBLE (Zext (in_rec.en_stat.hot_spot_cmd(1)))
	   rec1.en_stat.hot_spot_cmd(2) =
	1     DBLE (Zext (in_rec.en_stat.hot_spot_cmd(2)))
	   rec1.en_stat.power_a_status(1) =
	1     DBLE (Zext (in_rec.en_stat.power_a_status(1)))
	   rec1.en_stat.power_a_status(2) =
	1     DBLE (Zext (in_rec.en_stat.power_a_status(2)))
	   rec1.en_stat.power_b_status(1) =
	1     DBLE (Zext (in_rec.en_stat.power_b_status(1)))
	   rec1.en_stat.power_b_status(2) =
	1     DBLE (Zext (in_rec.en_stat.power_b_status(2)))

c eng_analog

	   Do m = 1,62
	      rec1.en_analog.grt(m) = DBLE (in_rec.en_analog.grt(m))
	      rec1.en_analog.group1(m) = DBLE (in_rec.en_analog.group1(m))
	   Enddo
	   Do m = 63,64
	      rec1.en_analog.grt(m) = DBLE (in_rec.en_analog.grt(m))
	   Enddo

c eng_tempdiff

	   Do m = 1,2
	      rec1.en_tempdiff(m).xcal = DBLE (in_rec.en_tempdiff(m).xcal)
	      rec1.en_tempdiff(m).ical = DBLE (in_rec.en_tempdiff(m).ical)
	      rec1.en_tempdiff(m).skyhorn =
	1        DBLE (in_rec.en_tempdiff(m).skyhorn)
	      rec1.en_tempdiff(m).refhorn =
	1        DBLE (in_rec.en_tempdiff(m).refhorn)
	      rec1.en_tempdiff(m).dihedral =
	1        DBLE (in_rec.en_tempdiff(m).dihedral)
	      rec1.en_tempdiff(m).collimator_mirror =
	1        DBLE (in_rec.en_tempdiff(m).collimator_mirror)
	      Do n = 1,4
	         rec1.en_tempdiff(m).bol_assem(n) =
	1           DBLE (in_rec.en_tempdiff(m).bol_assem(n))
	      Enddo
	   Enddo

c attitude

	   rec1.attitude.sun_moon_dist = DBLE (in_rec.attitude.sun_moon_dist)
	   rec1.attitude.cobe_moon_dist = DBLE (in_rec.attitude.cobe_moon_dist)
	   rec1.attitude.altitude = DBLE (in_rec.attitude.altitude)
	   rec1.attitude.projected_barycentric_velocity =
	1     DBLE (in_rec.attitude.projected_barycentric_velocity)
	   rec1.attitude.mcilwain_l_param =
	1     DBLE (in_rec.attitude.mcilwain_l_param)

c Calculate sigmas.

	   Do m = 1,10
	      sum.coad_spec_data.temp_sigma(m) =
	1        Dble (in_rec.coad_spec_data.temp_sigma(m))**2 * wm1
	   Enddo
	   sum.spec_data.resp_sigma =
	1     Dble (in_rec.spec_data.resp_sigma)**2 * wm1
	   sum.spec_data.tc_sigma = 
	1     Dble (in_rec.spec_data.tc_sigma)**2 * wm1
	   Do m = 1,62
	      sum.en_sigma.group2(m) =
	1        Dble (in_rec.en_sigma.group2(m))**2 * wm1
	      sum.en_sigma.sig_grt(m) =
	1        Dble (in_rec.en_sigma.sig_grt(m))**2 * wm1
	   Enddo
	   Do m = 63,64
	      sum.en_sigma.sig_grt(m) =
	1        Dble (in_rec.en_sigma.sig_grt(m))**2 * wm1
	   Enddo
c new
	   sum.coad_spec_data.sec_template.b_variance =
	1     Dble (in_rec.coad_spec_data.sec_template.b_variance) * wm1
	   sum.coad_spec_data.transient.c_variance =
	1     Dble (in_rec.coad_spec_data.transient.c_variance) * wm1

c------------------------------------------------------------------------------
	Else		! not first record; general case.

c coad_spec_data

	   Do m = 1,10
	      sum.coad_spec_data.temp(m) = sum.coad_spec_data.temp(m) +
	1        weight * ( Dble (in_rec.coad_spec_data.temp(m)) -
	2                   rec1.coad_spec_data.temp(m) )
	   Enddo
	   sum.coad_spec_data.glitch_rate = sum.coad_spec_data.glitch_rate +
	1     weight * ( Dble (in_rec.coad_spec_data.glitch_rate) -
	2                rec1.coad_spec_data.glitch_rate )
	   sum.coad_spec_data.adds_per_group =
	1     sum.coad_spec_data.adds_per_group + weight *
	2     ( Zext (in_rec.coad_spec_data.adds_per_group) -
	3       Nint (rec1.coad_spec_data.adds_per_group) )
	   sum.coad_spec_data.sweeps = sum.coad_spec_data.sweeps + weight *
	1     ( Zext (in_rec.coad_spec_data.sweeps) -
	2       Nint (rec1.coad_spec_data.sweeps) )
	   sum.coad_spec_data.peak_pos = sum.coad_spec_data.peak_pos + weight
	1     * ( Zext (in_rec.coad_spec_data.peak_pos) -
	2         Nint (rec1.coad_spec_data.peak_pos) )
	   sum.coad_spec_data.nyquist_hertz = sum.coad_spec_data.nyquist_hertz
	1     + weight * ( Dble (in_rec.coad_spec_data.nyquist_hertz) -
	2                  rec1.coad_spec_data.nyquist_hertz )
	   sum.coad_spec_data.nyquist_icm = sum.coad_spec_data.nyquist_icm +
	1     weight * ( Dble (in_rec.coad_spec_data.nyquist_icm) -
	2                rec1.coad_spec_data.nyquist_icm )

c The following block of quantities are weighted by the number of sweeps *
c the number of ifgs, as per the requirements of 27 May 1993.
c Not Anymore! LPR. August 1994.

	   sum.coad_spec_data.prim_template.amplitude =
	1     sum.coad_spec_data.prim_template.amplitude +
	2     weight * ( Dble (in_rec.coad_spec_data.prim_template.amplitude) -
	3                rec1.coad_spec_data.prim_template.amplitude )
	   sum.coad_spec_data.prim_template.snr =
	1     sum.coad_spec_data.prim_template.snr +
	2     weight * ( Dble (in_rec.coad_spec_data.prim_template.snr) -
	3                rec1.coad_spec_data.prim_template.snr )
	   sum.coad_spec_data.sec_template.amplitude =
	1     sum.coad_spec_data.sec_template.amplitude +
	2     weight * ( Dble (in_rec.coad_spec_data.sec_template.amplitude) -
	3                rec1.coad_spec_data.sec_template.amplitude )
	   sum.coad_spec_data.sec_template.snr =
	1     sum.coad_spec_data.sec_template.snr +
	2     weight * ( Dble (in_rec.coad_spec_data.sec_template.snr) -
	3                rec1.coad_spec_data.sec_template.snr )
	   sum.coad_spec_data.sec_template.b_average =
	1     sum.coad_spec_data.sec_template.b_average + weight *
	2     ( Dble (in_rec.coad_spec_data.sec_template.b_average) -
	3       rec1.coad_spec_data.sec_template.b_average )
	   sum.coad_spec_data.transient.c_average =
	1     sum.coad_spec_data.transient.c_average + weight *
	2     ( Dble (in_rec.coad_spec_data.transient.c_average) -
	3       rec1.coad_spec_data.transient.c_average )
	   sum.coad_spec_data.transient.bl_trans_coeff =
	1     sum.coad_spec_data.transient.bl_trans_coeff + weight *
	2     ( Dble (in_rec.coad_spec_data.transient.bl_trans_coeff) -
	3       rec1.coad_spec_data.transient.bl_trans_coeff )
	   sum.coad_spec_data.deglitch.glitch_iter =
	1     sum.coad_spec_data.deglitch.glitch_iter + weight *
	2     ( in_rec.coad_spec_data.deglitch.glitch_iter -
	3       rec1.coad_spec_data.deglitch.glitch_iter )
	   sum.coad_spec_data.deglitch.glitch_signal =
	1     sum.coad_spec_data.deglitch.glitch_signal + weight *
	2     ( Dble (in_rec.coad_spec_data.deglitch.glitch_signal) -
	3       rec1.coad_spec_data.deglitch.glitch_signal )
	   sum.coad_spec_data.noise = sum.coad_spec_data.noise + weight *
	1     ( Dble (in_rec.coad_spec_data.noise) -
	3       rec1.coad_spec_data.noise)
	   Do m=1,5
	      sum.coad_spec_data.bl_coeffs(m) = sum.coad_spec_data.bl_coeffs(m)
	1        + weight * ( Dble (in_rec.coad_spec_data.bl_coeffs(m)) -
	2                     rec1.coad_spec_data.bl_coeffs(m) )
	   Enddo
c
	   sum.coad_spec_data.sec_template.variance =
	1     sum.coad_spec_data.sec_template.variance +
	2     in_rec.coad_spec_data.sec_template.variance * weight
	   out_rec.coad_spec_data.gain_sum = out_rec.coad_spec_data.gain_sum +
	1         in_rec.coad_spec_data.gain_sum
	   sum.coad_spec_data.bol_cmd_bias =
	1     sum.coad_spec_data.bol_cmd_bias + weight *
	2     ( in_rec.coad_spec_data.bol_cmd_bias -
	3       rec1.coad_spec_data.bol_cmd_bias )
	   sum.coad_spec_data.bol_volt =
	1     sum.coad_spec_data.bol_volt + weight *
	2     ( Dble (in_rec.coad_spec_data.bol_volt) -
	3       rec1.coad_spec_data.bol_volt )

c spec_data

	   sum.spec_data.responsivity = sum.spec_data.responsivity + weight *
	1     ( Dble (in_rec.spec_data.responsivity) -
	2       rec1.spec_data.responsivity )
	   sum.spec_data.time_constant = sum.spec_data.time_constant + weight *
	1     ( Dble (in_rec.spec_data.time_constant) -
	2       rec1.spec_data.time_constant )
	   sum.spec_data.phase_corr = sum.spec_data.phase_corr + weight *
	1     ( Dble (in_rec.spec_data.phase_corr) - rec1.spec_data.phase_corr)
	   sum.spec_data.qrad = sum.spec_data.qrad + weight *
	1     ( Dble (in_rec.spec_data.qrad) - rec1.spec_data.qrad )
	   sum.spec_data.ir_power = sum.spec_data.ir_power + weight *
	1     ( Dble (in_rec.spec_data.ir_power) - rec1.spec_data.ir_power )

c phase_corr, qrad, and ir_power don't have sigmas to average.
c Instead, calculate a sigma for the FCS average:
c   sigma = sqrt [ sum {w*(p(i)-po)^2}  -  N * (<p>-po)^2 ] / (N-1)
c     where w = weight, N = sum w, p(i) = ith phase_corr, po = constant =
c     first phase_corr.  phase_corr_sqrd = sum {w*(p(i)-po)^2}.

	   sum.spec_data.phase_corr_sqrd = sum.spec_data.phase_corr_sqrd +
	1     weight * ( Dble (in_rec.spec_data.phase_corr) -
	2                rec1.spec_data.phase_corr )**2
	   sum.spec_data.qrad_sqrd = sum.spec_data.qrad_sqrd + weight *
	1     ( Dble (in_rec.spec_data.qrad) - rec1.spec_data.qrad )**2
	   sum.spec_data.ir_power_sqrd = sum.spec_data.ir_power_sqrd + weight *
	1     ( Dble (in_rec.spec_data.ir_power) - rec1.spec_data.ir_power )**2

c eng_status

	   Do m = 1,14
	      sum.en_stat.group1(m) = sum.en_stat.group1(m) + weight *
	1        ( in_rec.en_stat.group1(m) - Nint (rec1.en_stat.group1(m)) )
	      sum.en_stat.group2(m) = sum.en_stat.group2(m) + weight *
	1        ( Zext (in_rec.en_stat.group2(m)) -
	2          Nint (rec1.en_stat.group2(m)) )
	   Enddo
	   sum.en_stat.group1(15) = sum.en_stat.group1(15) + weight *
	1     ( in_rec.en_stat.group1(15) -
	2     Nint (rec1.en_stat.group1(15)) )
	   sum.en_stat.group1(16) = sum.en_stat.group1(16) + weight *
	1     ( in_rec.en_stat.group1(16) -
	2     Nint (rec1.en_stat.group1(16)) )
	   sum.en_stat.hot_spot_cmd(1) = sum.en_stat.hot_spot_cmd(1) + weight *
	1     ( Zext (in_rec.en_stat.hot_spot_cmd(1)) -
	2     Nint (rec1.en_stat.hot_spot_cmd(1)) )
	   sum.en_stat.hot_spot_cmd(2) = sum.en_stat.hot_spot_cmd(2) + weight *
	1     ( Zext (in_rec.en_stat.hot_spot_cmd(2)) -
	2     Nint (rec1.en_stat.hot_spot_cmd(2)) )
	   sum.en_stat.power_a_status(1) = sum.en_stat.power_a_status(1) +
	1     weight * ( Zext (in_rec.en_stat.power_a_status(1)) -
	2                Nint (rec1.en_stat.power_a_status(1)) )
	   sum.en_stat.power_a_status(2) = sum.en_stat.power_a_status(2) +
	1     weight * ( Zext (in_rec.en_stat.power_a_status(2)) -
	2                Nint (rec1.en_stat.power_a_status(2)) )
	   sum.en_stat.power_b_status(1) = sum.en_stat.power_b_status(1) +
	1     weight * ( Zext (in_rec.en_stat.power_b_status(1)) -
	2                Nint (rec1.en_stat.power_b_status(1)) )
	   sum.en_stat.power_b_status(2) = sum.en_stat.power_b_status(2) +
	1     weight * ( Zext (in_rec.en_stat.power_b_status(2)) -
	2                Nint (rec1.en_stat.power_b_status(2)) )

c eng_analog

	   Do m = 1,62
	      sum.en_analog.grt(m) = sum.en_analog.grt(m) + weight *
	1        ( Dble (in_rec.en_analog.grt(m)) -
	2          rec1.en_analog.grt(m) )
	      If (in_rec.en_analog.group1(m) .EQ. -9999.) Then
	         rec1.en_analog.group1(m) = -9999.
	      Else
	         sum.en_analog.group1(m) = sum.en_analog.group1(m) + weight *
	1           ( Dble (in_rec.en_analog.group1(m)) -
	2             rec1.en_analog.group1(m) )
	      EndIf
	   Enddo
	   Do m = 63,64
	      sum.en_analog.grt(m) = sum.en_analog.grt(m) + weight *
	1        ( Dble (in_rec.en_analog.grt(m)) -
	2          rec1.en_analog.grt(m) )
	   Enddo

c eng_tempdiff

	   Do m = 1,2
	      sum.en_tempdiff(m).xcal = sum.en_tempdiff(m).xcal + weight *
	1        ( Dble (in_rec.en_tempdiff(m).xcal) -
	2          rec1.en_tempdiff(m).xcal )
	      sum.en_tempdiff(m).ical = sum.en_tempdiff(m).ical + weight *
	1        ( Dble (in_rec.en_tempdiff(m).ical) -
	2          rec1.en_tempdiff(m).ical )
	      sum.en_tempdiff(m).skyhorn = sum.en_tempdiff(m).skyhorn + weight
	1        * ( Dble (in_rec.en_tempdiff(m).skyhorn) -
	2            rec1.en_tempdiff(m).skyhorn )
	      sum.en_tempdiff(m).refhorn = sum.en_tempdiff(m).refhorn + weight
	1        * ( Dble (in_rec.en_tempdiff(m).refhorn) -
	2            rec1.en_tempdiff(m).refhorn )
	      sum.en_tempdiff(m).dihedral = sum.en_tempdiff(m).dihedral +
	1        weight * ( Dble (in_rec.en_tempdiff(m).dihedral) -
	2                   rec1.en_tempdiff(m).dihedral )
	      sum.en_tempdiff(m).collimator_mirror =
	1        sum.en_tempdiff(m).collimator_mirror + weight *
	2        ( Dble (in_rec.en_tempdiff(m).collimator_mirror) -
	3          rec1.en_tempdiff(m).collimator_mirror )
	      Do n = 1,4
	         sum.en_tempdiff(m).bol_assem(n) =
	1           sum.en_tempdiff(m).bol_assem(n)
	2           + weight * ( Dble (in_rec.en_tempdiff(m).bol_assem(n)) -
	3                        rec1.en_tempdiff(m).bol_assem(n) )
	      Enddo
	   Enddo

c attitude

	   sum.attitude.sun_moon_dist = sum.attitude.sun_moon_dist + weight *
	1      ( Dble (in_rec.attitude.sun_moon_dist) -
	2        rec1.attitude.sun_moon_dist )
	   sum.attitude.cobe_moon_dist = sum.attitude.cobe_moon_dist + weight *
	1      ( Dble (in_rec.attitude.cobe_moon_dist) -
	2        rec1.attitude.cobe_moon_dist )
	   sum.attitude.altitude = sum.attitude.altitude + weight *
	1      ( DBLE (in_rec.attitude.altitude) - rec1.attitude.altitude )
	   sum.attitude.projected_barycentric_velocity =
	1      sum.attitude.projected_barycentric_velocity + weight *
	2      ( DBLE (in_rec.attitude.projected_barycentric_velocity) -
	3        rec1.attitude.projected_barycentric_velocity )
	   sum.attitude.mcilwain_l_param = sum.attitude.mcilwain_l_param +
	1     weight * ( DBLE (in_rec.attitude.mcilwain_l_param) -
	2                rec1.attitude.mcilwain_l_param )

c Calculate sigmas.

	   Do m = 1,10
	      sum.coad_spec_data.temp_sigma(m) =
	1        sum.coad_spec_data.temp_sigma(m) +
	2        Dble (in_rec.coad_spec_data.temp_sigma(m))**2 * wm1
	3        + weight * (in_rec.coad_spec_data.temp(m) -
	4                    rec1.coad_spec_data.temp(m))**2
	   Enddo

	   sum.spec_data.resp_sigma = sum.spec_data.resp_sigma +
	1     Dble (in_rec.spec_data.resp_sigma)**2 * wm1 +
	2     weight * (in_rec.spec_data.responsivity -
	3               rec1.spec_data.responsivity)**2

	   sum.spec_data.tc_sigma = sum.spec_data.tc_sigma +
	1     Dble (in_rec.spec_data.tc_sigma)**2 * wm1 +
	2     weight * (in_rec.spec_data.time_constant - 
	3               rec1.spec_data.time_constant)**2

c variances
	   sum.coad_spec_data.sec_template.b_variance =
	1     sum.coad_spec_data.sec_template.b_variance +
	2     Dble (in_rec.coad_spec_data.sec_template.b_variance) * wm1 +
	3     weight * (in_rec.coad_spec_data.sec_template.b_average -
	4               rec1.coad_spec_data.sec_template.b_average)**2

	   sum.coad_spec_data.transient.c_variance =
	1     sum.coad_spec_data.transient.c_variance +
	2     Dble (in_rec.coad_spec_data.transient.c_variance) * wm1 +
	3     weight * (in_rec.coad_spec_data.transient.c_average -
	4               rec1.coad_spec_data.transient.c_average)**2

	   Do m = 1,62
	      If (in_rec.en_analog.group1(m) .EQ. -9999.) Then
	         sum.en_sigma.group2(m) = 1.E15       ! Flag value
	      Else
	         sum.en_sigma.group2(m) = sum.en_sigma.group2(m) +
	1           Dble (in_rec.en_sigma.group2(m))**2 * wm1 +
	2           weight * (in_rec.en_analog.group1(m) -
	3                     rec1.en_analog.group1(m))**2
	      EndIf
	      sum.en_sigma.sig_grt(m) = sum.en_sigma.sig_grt(m) +
	1        Dble (in_rec.en_sigma.sig_grt(m))**2 * wm1 +
	2        weight * (in_rec.en_analog.grt(m) -
	3                  rec1.en_analog.grt(m))**2
	   Enddo
	   Do m = 63,64
	      sum.en_sigma.sig_grt(m) = sum.en_sigma.sig_grt(m) +
	1        Dble (in_rec.en_sigma.sig_grt(m))**2 * wm1 +
	2        weight * (in_rec.en_analog.grt(m) -
	3                  rec1.en_analog.grt(m))**2
	   Enddo
	Endif				! If first record (2nd line of code)
c------------------------------------------------------------------------------

c  Theses attitude sums are the same for both 1st record and otherwise.

	Do m = 1,3
	   sum.attitude.equat(m) = sum.attitude.equat(m) +
	1     in_rec.attitude.equatorial(m) * weight
	Enddo
	terr_lat  = in_rec.attitude.terr_latitude * fac_att_conv
	terr_long = in_rec.attitude.terr_longitude * fac_att_conv
	sum.attitude.terr(1) = sum.attitude.terr(1) + Cosd (terr_lat) *
	1   Cosd (terr_long) * weight
	sum.attitude.terr(2) = sum.attitude.terr(2) + Cosd (terr_lat) *
	1   Sind (terr_long) * weight
	sum.attitude.terr(3) = sum.attitude.terr(3) + Sind (terr_lat) * weight
	elimb  = in_rec.attitude.earth_limb * fac_att_conv - 90.0
	el_az  = in_rec.attitude.earth_limb_azimuth * fac_att_conv
	sum.attitude.elimb(1) = sum.attitude.elimb(1) +
	1   Sind (elimb) * Cosd (el_az) * weight
	sum.attitude.elimb(2) = sum.attitude.elimb(2) +
	1   Sind (elimb) * Sind (el_az) * weight
	sum.attitude.elimb(3) = sum.attitude.elimb(3) + Cosd (elimb) * weight
	sum.attitude.sun_angle = sum.attitude.sun_angle +
	1   in_rec.attitude.sun_angle * weight
	moon_angle  = in_rec.attitude.moon_angle * fac_att_conv - 90.0
	moon_az     = in_rec.attitude.moon_az_angle * fac_att_conv
	sum.attitude.mnang(1) = sum.attitude.mnang(1) +
	1   Sind (moon_angle) * Cosd (moon_az) * weight
	sum.attitude.mnang(2) = sum.attitude.mnang(2) +
	1   Sind (moon_angle) * Sind (moon_az) * weight
	sum.attitude.mnang(3) = sum.attitude.mnang(3) + Cosd (moon_angle) *
	1   weight
	gal_lat     = in_rec.attitude.galactic_latitude * fac_att_conv
	gal_long    = in_rec.attitude.galactic_longitude * fac_att_conv
	sum.attitude.gal(1) = sum.attitude.gal(1) +
	1   Cosd (gal_lat) * Cosd (gal_long) * weight
	sum.attitude.gal(2) = sum.attitude.gal(2) +
	1   Cosd (gal_lat) * Sind (gal_long) * weight
	sum.attitude.gal(3) = sum.attitude.gal(3) + Sind (gal_lat) * weight
	ecl_lat     = in_rec.attitude.ecliptic_latitude * fac_att_conv
	ecl_long    = in_rec.attitude.ecliptic_longitude * fac_att_conv
	sum.attitude.ecl(1) = sum.attitude.ecl(1) + 
	1   Cosd (ecl_lat) * Cosd (ecl_long) * weight
	sum.attitude.ecl(2) = sum.attitude.ecl(2) +
	1   Cosd (ecl_lat) * Sind (ecl_long) * weight
	sum.attitude.ecl(3) = sum.attitude.ecl(3) + Sind (ecl_lat) * weight
	sum.attitude.geovel = sum.attitude.geovel +
	1   in_rec.attitude.projected_geocentric_velocity * weight
	sum.attitude.sc_rot_angle = sum.attitude.sc_rot_angle +
	1   in_rec.attitude.sc_rotation_angle * weight
	sum.attitude.exc_gal_lat = sum.attitude.exc_gal_lat +
	1   in_rec.attitude.exc_galactic_lat * weight
	moon_phase = in_rec.attitude.moon_phase * fac_att_conv_rad
	orb_phase  = in_rec.attitude.orbital_phase * fac_att_conv_rad
	scan_angle = in_rec.attitude.scan_angle * fac_att_conv_rad
	Return
	End
