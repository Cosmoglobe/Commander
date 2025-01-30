	INTEGER*4 FUNCTION FEP_INVERT (PCOEFS,XX,YY)

C--------------------------------------------------------------------------
C
C  Solves for the 2nd order polynomial (parabola) that goes through 3
C  consecutve entries in the lookup tables.
C
C  Input Parameters:  X = Array of 3 consecutive independent vars from table
C		      Y = 3 corresponding dependent variables from table
C
C  Output Parameter:  PCOEFS = Array of 3 parabola coefficients P0,P1,P2 such
C			       that Y(n) = P0 + P1*X(n) + P2*X(n)**2
C
C  This subroutine uses the method of inverting a 3 x 3 matrix by expansion
C  in 2 x 2 minors. This technique uses only 16 adds and 26 multiplies.
C
C  First calculate frequently-used quantities to save time
C
C	Modified by Don Stevens-Rayburn on 10-Jul-1989 to use REAL * 8 
C		arithmetic in computing the coefficients of the 
C		interpolation polynomical.
C
C---------------------------------------------------------------------------

	IMPLICIT	NONE

	INTEGER		*4	I

	REAL		*4	XX(1)
	REAL		*4	YY(1)
	REAL		*4	PCOEFS(1)

	REAL		*8	X(3)
	REAL		*8	Y(3)
	REAL		*8	X0X1
	REAL		*8	X0X2
	REAL		*8	X1X2
	REAL		*8	X1MX0
	REAL		*8	X2MX1
	REAL		*8	X2MX0
	REAL		*8	DET

	EXTERNAL	FEP_NORMAL
	EXTERNAL	FEP_INVTMAT

	DO I = 1,3
	   X(I)=XX(I)
	   Y(I)=YY(I)
	ENDDO

	X0X1 = X(1) * X(2)
	X0X2 = X(1) * X(3)
	X1X2 = X(2) * X(3)
	X1MX0 = X(2) - X(1)
	X2MX1 = X(3) - X(2)
	X2MX0 = X(3) - X(1)
C
C  Calculate determinant of 3 x 3 matrix for normalization. Return error 
C  if matrix is singular.
C
	DET = X1X2*X2MX1 - X(1)*X2MX1*(X(2)+X(3)) + X(1)*X(1)*X2MX1
	IF (ABS(DET) .LT. 1.E-25) THEN
	   FEP_INVERT = %LOC(FEP_INVTMAT)
	   CALL LIB$SIGNAL(FEP_INVTMAT)
	   RETURN
	ENDIF
C
C  Now calculate coefficients
C
	PCOEFS(1) = (X1X2*X2MX1*Y(1) - X0X2*X2MX0*Y(2)
     .					 + X0X1*X1MX0*Y(3))/DET
	PCOEFS(2) = (-X2MX1*(X(2)+X(3))*Y(1)
     .		     +X2MX0*(X(1)+X(3))*Y(2)
     .		     -X1MX0*(X(1)+X(2))*Y(3))/DET	
	PCOEFS(3) = (X2MX1*Y(1)-X2MX0*Y(2)+X1MX0*Y(3))/DET

	FEP_INVERT = %LOC(FEP_NORMAL)

	RETURN
	END
