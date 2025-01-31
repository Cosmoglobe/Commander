	SUBROUTINE FUT_SYMCENTER (DATA,ICENTR,SYM)

C------------------------------------------------------------------------
C    PURPOSE:
C
C  This subroutine takes the data array DATA and calculates a symmetry
C  parameter ISYM in a 21-point wide window centered on index number ICENTR.
C  The algorithm is to calculate (1) a weighted sum of the absolute values of
C  the data points in the window, and (2) the sum of the absolute differences
C  between points on the left side of the window and their reflections on the
C  right side. The difference bewteen these two quantities is the symmetry
C  parameter SYM; the larger SYM, the more symmetrical the window about ICENTR
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rich Isaacman
C            ARC
C
C    INVOCATION: CALL FUT_SYMCENTER(IFG DATA, IFG CENTER, SYMMETRY PARAMETER)
C
C    INPUT PARAMETERS:
C	     IFG
C	     Point about which symmetry parameter is to be calculated
C
C    OUTPUT PARAMETERS: 
C	     Symmetry parameter
C
C    SUBROUTINES CALLED: None
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: None
C
C----------------------------------------------------------------------

	IMPLICIT NONE
	REAL*4 SYM
	REAL*4 SUM
	REAL*4 DIF
	REAL*4 DATA(1)
	INTEGER*4 J
	INTEGER*4 ICENTR
	SUM = 2. * ABS(DATA(ICENTR))
	DIF = 0.
	DO J=1,10
	   SUM = SUM + ABS(DATA(ICENTR-J)) + ABS(DATA(ICENTR+J))
	   DIF = DIF + ABS(DATA(ICENTR-J) - DATA(ICENTR+J))
	ENDDO
	SYM = SUM - DIF
	RETURN
	END
