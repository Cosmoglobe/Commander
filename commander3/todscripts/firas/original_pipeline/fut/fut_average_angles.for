	Integer*4 Function FUT_Average_Angles ( angles, num, avg_angle )

C------------------------------------------------------------------------
C    PURPOSE: Calculates the average angle from a list of input angles.
C	      All calculations are performed in radians.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            ST Systems Corporation
C            September 15, 1988
C
C    INVOCATION: STATUS = FUT_Average_Angles ( ANGLES, NUM, AVG_ANGLE )
C
C    INPUT PARAMETERS:
C	ANGLES(NUM)	R*4	List in angles to be averaged.
C	NUM		I*4	Number of angles in list.
C
C    OUTPUT PARAMETERS:
C	AVG_ANGLE	R*4	Average angle.
C	STATUS		I*4	Success status.
C
C    SUBROUTINES CALLED:
C	FUT_Params
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: None
C
C------------------------------------------------------------------------------
C
C    REVISIONS:
C
C    1990 Feb 14,  Fred Shuman,  STX.  SPR 5699, Wrong ave sometimes returned,
C      e.g., inputs = (pi-.1, -pi+.3).  Note that the proposed solution given
C      there is also sometimes wrong, e.g., inputs = (-.1, .1).  This routine
C      assumes the input NUMERICAL VALUES all lie within some interval of length
C      2 pi.  If the input ANGLES are not covered by any semicircle, the result
C      returned by this routine may be untrustworthy, but 'average angle' begins
C      to lose its meaning as the angles being averaged spread out.  If an
C      average MUST be found in such a case, alter this routine to average the
C      unit vectors.  The length of the average vector is then a good indicator
C      of the lack of 'dispersion' of the inputs, and thus of the reliability of
C      the resulting average angle. 
C------------------------------------------------------------------------------

	Implicit None

	Include		'(FUT_Params)'

	Integer		*4	num
	Real		*4	angles(num)
	Real		*4	avg_angle

	Integer		*4	i
	Real		*4	sum_angle
	Real		*4	corr_angle

	External	FUT_Normal
	External	FUT_InvAngle

C
	sum_angle = 0.0

	Do i=1,num
C
C  Subtract the first value.  This will put all the values between +/- 2 pi.  If
C    they are not between +/- pi, get them there by adding or subtracting 2 pi.
C
	  corr_angle = angles(i) - angles(1)

	  If (corr_angle .Gt. fac_pi) Then
	    corr_angle = corr_angle - 2.0 * fac_pi
	  Else If (corr_angle .Lt. -fac_pi) Then
	    corr_angle = corr_angle + 2.0 * fac_pi
	  End If

	  sum_angle = sum_angle + corr_angle

	End Do
C
C   If the ave is not between +/- pi, add or subtract 2 pi.
C
	avg_angle = (sum_angle / num) + angles(1)
	If (avg_angle .Gt. fac_pi) Then
	  avg_angle = avg_angle - 2.0 * fac_pi
	Else If (avg_angle .Lt. -fac_pi) Then
	  avg_angle = avg_angle + 2.0 * fac_pi
	End If

	Return
	End
