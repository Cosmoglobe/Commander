	INTEGER*4 FUNCTION FUT_GET_PLOT_PARMS
	1				(S_PARAMS
	2				,S_VAR
	3				,S_COVAR
	4				,S_PLOT
	5				,S_PLOT_SIG)

C-FUT_GET_PLOT_PARAMS-Compute Plot Values from Calibration Model Solution
C
C The calibration model solution routine 'FFC_SOLVE_WEIGHTED_MODEL'
C returns the product H(f) x Ex(f) for all the 'emissivity' type parameters.
C This is because we do a simple weighted linear regression, and the natural
C parametrization for the model has this product as a model parameter.
C For plotting and display of these parameters, it is desirable to
C have the actual emissivities (and errors for the emissivities).
C
C Given the calibration model parameters, variances and covariance matrices,
C this subroutine returns the emissivities and the errors for the emissivities.
C Note that the covariances are not converted, only the errors (i.e., we
C do not compute the correlations for the actual emissivities with
C each other or with any other parameters).
C
C Calling Sequence:	ISTAT = FUT_GET_PLOT_PARMS
C					(S_PARAMS
C					,S_VAR
C					,S_COVAR
C					,S_PLOT
C					,S_PLOT_SIG)
C
C Parameters:
C
C	ISTAT		i*4	Return status.
C	S_PARAMS	r*4	Array of 6 x 512 elements returning the model
C				solution parameters for each frequency.  The
C				first index definitions are:
C
C				    1 -	Transfer function.
C				    2 -	ICAL emissivity.
C				    3 -	Skyhorn emissivity.
C				    4 -	Reference horn emissivity.
C				    5 - Structure emissivity.
C				    6 -	Constant term.
C
C				Parameters not fit should have value
C				0. and variance/covariance -1.0e+10. 
C				This array is field <SOLUTION.MODEL_SOLN>
C				in file 'FFC_MODEL'.
C	S_VAR		r*4	Array of 6 x 512 elements returning the
C				variances for the parameters in <S_PARAMS>.
C				This array corresponds to field <SOLUTION.
C				PARAM_VAR> in file 'FFC_MODEL'.
C	S_COVAR		r*4	Array of the off-diagonal terms of the 512
C				covariance matrices, one for each frequency,
C				each of which is derived from a symmetrical
C				6 x 6 matrix (where 6 is the maximum number of
C				parameters in the model solution) but without
C				the diagonal terms which are stored in
C				<S_VAR>.  Each individual matrix is stored
C				in 'modified symmetrical storage mode'.
C				See statement function <MSSMI>. This array
C				is field <FIT_COVAR.MODEL_COVAR> in file
C				'FFC_MODEL_COVAR'.
C	S_PLOT		r*4	Output 512 x 6 array. The first index is the
C				frequency index, the second index has the
C				same interpretation as the first index in
C				<S_PARAMS>.  Note that the transfer function
C				and constant parameters are copied from the
C				input arrays with no changes (it is easier
C				to plot from this array than from the input
C				array).  Undefined parameters are set to 0.
C	S_PLOT_SIG	r*4	Output 512 x 6 array specifying the statistical
C				errors for the numbers in <S_PLOT>. Undefined
C				parameters are flagged with errors of -1.
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C MODIFICATIONS:
C
C   NAME	  DATE				CHANGE
C E. Cheng	08-Jul-87	Original version.
C		02-Aug-87	Change record structure to match current
C				thinking.
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C
C
C COBE CSDR requirements
C
	implicit NONE
C
C Maximum number of parameters (must have the same value as the parameter in 
C 'FFC_SOLVE_WEIGHTED_MODEL')
C
	include '(FUT_PARAMS)'
C
C Calling parameters
C
	real*4 S_PARAMS(fac_max_struct,512)
	real*4 S_VAR(fac_max_struct,512)
	real*4 S_COVAR(((fac_max_struct-1)*fac_max_struct)/2,512)
	real*4 S_PLOT(512,fac_max_struct)
	real*4 S_PLOT_SIG(512,fac_max_struct)
C
C Local variables
C
	real B1INV			!inverse of transfer function
	real B1INV2,B1INV3,B1INV4	!various powers of above
	integer I,J			!small change
	integer IFREQ			!frequency loop index
	integer IP			!parameter loop index
	real S				!temporary parameter value
	real VII,VI1,V11		!temporary covariances
C
C External routines
C
	external FUT_NORMAL
C
C 'Modified symmetrical storage mode' indexing function (IMSL symmetrical
C storage mode without the diagonal elements).
C
C This function is valid for (I .gt. J) and (I, J .gt. 1).  The first row
C and column are labeled 1 (not 0).
C
	integer MSSMI
	MSSMI(I,J) = ((I-1)*(I-2))/2 + J
C
C Loop over all frequency bins
C
	do IFREQ = 1, 512
C
C Transfer function is valid if the number is positive, will not
C over/underflow, and the variance is positive.
C
		if (S_PARAMS(1,IFREQ).gt.1.0e-08 .and.
	1	    S_PARAMS(1,IFREQ).lt.1.0e+08 .and.
	2	    S_VAR(1,IFREQ).gt.0.) then
C
C Copy transfer function
C
			S_PLOT(IFREQ,1) = S_PARAMS(1,IFREQ)
			S_PLOT_SIG(IFREQ,1) = sqrt(S_VAR(1,IFREQ))
C
C Save inverse powers of transfer function for later
C
			B1INV = 1. / S_PARAMS(1,IFREQ)
			B1INV2 = B1INV * B1INV
			B1INV3 = B1INV2 * B1INV
			B1INV4 = B1INV2 * B1INV2
C
C Copy constant if necessary
C
			if (S_VAR(fac_max_struct,IFREQ) .gt. 0.) then
				S_PLOT(IFREQ,fac_max_struct) = S_PARAMS(fac_max_struct,IFREQ)
				S_PLOT_SIG(IFREQ,fac_max_struct) = sqrt(S_VAR(fac_max_struct,IFREQ))
C
C Constant not defined
C
			else
				S_PLOT(IFREQ,fac_max_struct) = 0.
				S_PLOT_SIG(IFREQ,fac_max_struct) = 0.
			endif
C
C Loop on emissivity type parameters
C
			do IP = 2, fac_max_struct-1
			    VII = S_VAR(IP,IFREQ)
			    VI1 = S_COVAR(MSSMI(IP,1),IFREQ)
			    V11 = S_VAR(1,IFREQ)
			    if (VII .gt. 0.) then
				S = S_PARAMS(IP,IFREQ)
				S_PLOT(IFREQ,IP) = S * B1INV
				S_PLOT_SIG(IFREQ,IP) = sqrt(abs(
	1						  VII*B1INV2
	2						+ S*S*V11/B1INV4
	3						- S*VI1/B1INV3
	4								))
			    else
				S_PLOT(IFREQ,IP) = 0.	
				S_PLOT_SIG(IFREQ,IP) = 0.
			    endif
			enddo
C
C If the transfer function is bad, then all other parameters are assumed bad
C
		else
			do IP = 1, fac_max_struct
				S_PLOT(IFREQ,IP) = 0.
				S_PLOT_SIG(IFREQ,IP) = 0.
			enddo
		endif
C
C End of loop on frequency bins
C
	enddo
C
C Normal completion
C
	FUT_GET_PLOT_PARMS = %loc(FUT_NORMAL)
C
	return
	end
