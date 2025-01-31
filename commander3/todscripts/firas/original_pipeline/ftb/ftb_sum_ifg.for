
C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Sum_Ifg ( Ifg )

C-------------------------------------------------------------------------------
C
C	Purpose: To compute the sum of the absolute values of the 512 points
C	         of a FIRAS interferogram and return this sum as the function
c		 value.
C
C	Author: Shirley M. Read
C		STX, November 8, 1988
C
C	Invocation: Status = FTB_Sum_Ifg ( Ifg )
C
CH	Change Log:
CH
C	  ----------------------------------------------------------------------
C
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Ifg (512)     I*2             A FIRAS interferogram of 512 points
C	
C	Output Function Value:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Sum_Ifg       I*4             Sum of absolute values of Ifg points
C	
C	Subroutines Called:
C
C	  ABS -- FORTRAN absolute value function
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C
C	Processing Method:
C
C	  Initialize the sum accumulator to zero.
C	  Do for 512 Interferogram points
C	    Add the absolute value of the point to the sum.
C	  Enddo
C	  Set the function value to the sum.
C	  Return.
C
C------------------------------------------------------------------------------
C
	Implicit None
	
C	Passed Parameters.

	Integer*2 	Ifg (512)	! FIRAS Interferogram of 512 points

C	Local Declarations.	

	Integer*4	Sum		! Sum accumulator
	Integer*4	Zero / 0 /      ! Initial value
	Integer*4	Ix		! Index
	Integer*4       Ifgi4           ! I*4 Ifg for function Abs

	Sum = Zero

	Do Ix = 1, 512

	    Ifgi4 = Ifg(Ix)
	    Sum = Sum + Abs ( Ifgi4 )

	Enddo
	
	FTB_Sum_Ifg = Sum

	Return 
	End
