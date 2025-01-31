	Integer*4 Function FUT_AttOffset(ain, aout)
C------------------------------------------------------------------------------
C    PROGRAM NAME:
C	FUT_AttOffset
C
C    PURPOSE:  Take the attitude tensor (a proper orthogonal matrix) and apply
C	a known, small offset rotation to bring the spacecraft negative spin
C	axis to the FIRAS beam direction.
C
C    AUTHOR: Fred Shuman,  ST Systems Corporation,  1991 April 25
C
C    INVOCATION: status = FUT_AttOffset(ain, aout)
C
C    INPUT ARGUMENTS:
C	AIN(3,3)	R*4	Attitude matrix from UAX_Get_Attitude
C
C    OUTPUT ARGUMENT:
C	AOUT(3,3)	R*4	Attitude matrix, offset to the FIRAS beam
C
C    SUBROUTINES CALLED:  None
C
C    INCLUDE FILES:  None
C
C    PROCESSING METHOD:
C------------------------------------------------------------------------------

	Implicit	None

	Real	* 4	ain(3,3), aout(3,3)

	Integer	* 4	j, k

	External fut_normal

	Do j=1,3
	   Do k=1,3
	      aout(j,k) = ain(j,k)
	   End Do
	End Do

!	Set return status to success.

	FUT_AttOffset = %Loc(FUT_Normal)

	Return

	End
