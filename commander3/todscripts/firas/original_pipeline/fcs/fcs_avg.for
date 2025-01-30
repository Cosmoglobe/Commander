	Integer*4  Function  FCS_Avg  ( total, sum_weight, fcs_rec, maxsum,
	1                               mphase, ophase, sangle,
	2                               sum, rec1, time_array, time_weights,
	3                               numrecs, sum_dof, scan )

c-----------------------------------------------------------------------------
c Function to compute final averages of values, spectra and variances.
c Author: Larry P. Rosen, Hughes STX, March 1992
c 23 February 1993. Corrected mixed mode in encode.
c See requirements and design for details.
c Modified:  Larry Rosen, October 1993.  All the fields except for spectra
c    and variance have been moved into subroutine FUT_FCS_Avg.
c-----------------------------------------------------------------------------
	Implicit None

c Include.

	Include		'(fut_params)'
	Include		'(fcs_msg)'
	Include		'(fut_fcs_include)'

c Passed parameters.

	Integer*4	total			! total number of ifgs
	Integer*4	maxsum			! max sum of nifgs per record
	Real*8		sum_weight	! Sum of weights for average (= total)
	Real*4		mphase (maxsum)		! Moon phase array
	Real*4		ophase (maxsum)		! Orbital phase array
	Real*4		sangle (maxsum)		! Sun angle array

	Dictionary	'FCS_SKY'
	Record / fcs_sky / fcs_rec		! Output averaged map record
	Record / double_fcf / sum		! Record holding sums.
	Record / double_fcf / rec1		! Record storing 1st data points
	Integer*4	time_array (2, maxsum)	! Time array to average
	Integer*2	numrecs			! Number of records to average
	Integer*4	time_weights (maxsum)	! Time average weights = nifgs
	Real*8		sum_dof		! sum of degrees of freedom (dof) where
					!   dof = weight - k  where k=2 if 2nd
					!   template subtracted, 1 if NOT.
	Character*2	scan		! scan mode from command line

c Function
	Integer*4	FUT_FCS_Avg

c Local
	Real*8		wt_dof		! sum of weights times sum of dof
	Integer*4	rstat		! return status
	Integer*4	i		! counter
	Real*8		variance
	Character*1	xcal_pos
	Integer*2	channel
	Real*8		real_part, imag_part
	Real*8		real_part2, imag_part2
c-----------------------------------------------------------------------------
c Begin

	FCS_Avg = %loc (fcs_normal)

	fcs_rec.coad_spec_head.num_ifgs = total

	channel = fcs_rec.coad_spec_data.chan_id
	fcs_rec.coad_spec_data.bol_cmd_bias =
	1   fcs_rec.en_stat.bol_cmd_bias (channel)
	fcs_rec.coad_spec_data.bol_volt = fcs_rec.en_analog.bol_volt (channel)

C Spectra and its variances.  If values for spectra are too small use complex 0.

	Do i = 1, 257
	   real_part = Dreal (sum.spec_data.spec (i))
	   imag_part = Dimag (sum.spec_data.spec (i))
	   If ( Abs (real_part) .LE. 1.0e-25) Then
	      real_part = 0.0
	   Else
	      real_part = real_part / sum_weight
	   Endif
	   If ( Abs (imag_part) .LE. 1.0e-25) Then
	      imag_part = 0.0
	   Else
	      imag_part = imag_part / sum_weight
	   Endif
	   real_part2 = Dreal (rec1.spec_data.spec(i))
	   If ( (Abs (real_part) .LE. 1.0e-36) .AND.
	1       (Abs (real_part2) .LE. 1.0e-36) ) Then
	      real_part = 0.0
	   Else
	      real_part = real_part + real_part2
	   Endif
	   imag_part2 = Dimag (rec1.spec_data.spec(i))
	   If ( (Abs (imag_part) .LE. 1.0e-36) .AND.
	1       (Abs (imag_part2) .LE. 1.0e-36) ) Then
	      imag_part = 0.0
	   Else
	      imag_part = imag_part + imag_part2
	   Endif
	   fcs_rec.spec_data.spec(i) =  Cmplx (real_part, imag_part)

	   wt_dof = sum_weight * sum_dof

	   fcs_rec.spec_data.real_var(i) = sum.spec_data.real_var(i) / wt_dof
	   fcs_rec.spec_data.imag_var(i) = sum.spec_data.imag_var(i) / wt_dof
	   fcs_rec.spec_data.real_imag_var(i) =
	1     sum.spec_data.real_imag_var(i) / wt_dof
	Enddo

c Call FUT_FCS_Avg to average engineering and attitude quantities.

	rstat = FUT_FCS_Avg  ( maxsum, time_array, numrecs, time_weights,
	1                      fcs_rec, sum, rec1, sum_weight,
	2                      mphase, ophase, sangle, total )

c Use Encode to create text label.

	If (fcs_rec.coad_spec_data.xcal_pos .EQ. fac_xcalout) Then
	   xcal_pos = 'O'
	Else
	   xcal_pos = 'I'
	Endif

	If (fcs_rec.attitude.pixel_no .GE. 0) Then
	   Encode (60, 10, fcs_rec.coad_spec_head.label)
	1     fcs_rec.coad_spec_head.first_gmt,
	2     fcs_rec.coad_spec_head.last_gmt,
	3     fcs_rec.coad_spec_head.num_ifgs,
	4     scan,
	5     Nint (fcs_rec.coad_spec_data.ical * 1000.0),
	6     xcal_pos, fcs_rec.attitude.pixel_no
  10	   Format (2(A,'_'),I3.3,'_',A,'_',I5.5,'_',A,'_',I4.4,11X)
	Else
	   Encode (60, 20, fcs_rec.coad_spec_head.label)
	1     fcs_rec.coad_spec_head.first_gmt,
	2     fcs_rec.coad_spec_head.last_gmt,
	3     fcs_rec.coad_spec_head.num_ifgs,
	4     scan,
	5     Nint (fcs_rec.coad_spec_data.ical * 1000.0),
	6     xcal_pos
  20	   Format (2(A,'_'),I3.3,'_',A,'_',I5.5,'_',A,16X)
	Endif

	Return
	End
