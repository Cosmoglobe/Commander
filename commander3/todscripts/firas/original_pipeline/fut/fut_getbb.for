	SUBROUTINE FUT_GETBB (TEMPK, FNYQ, NPTS, SPECTRUM)

C-------------------------------------------------------------------------
C
C  Function
C  	Returns NPT blackbody spectrum at temperature TEMPK up
C	to frequency FNYQ.
C
C  Written by Rich Isaacman (Applied Research Corp.)
C
C  CALLING SEQUENCE:
C	SPECTRUM = FUT_GETBB(TEMPK,FNYQ,NPTS,SPECTRUM)
C
C  INPUT:
C	TEMPERATURE IN DEGREES KELVIN
C	NYQUIST FREQUENCY
C	NUMBER OF POINTS
C
C  OUTPUT:
C	BLACKBODY SPECTRUM
C
C  SUBROUTINES CALLED:
C	PLANCK
C
C  INCLUDE FILES: NONE
C
C-------------------------------------------------------------------------
C
C  Changes:
C
C	R. Kummerer, January 12, 1987. Remove multiplication by DF factor.
C
C	R. Kummerer, March 30, 1987. Move from FFC to FUT facility.
C
C	R. Kummerer, January 11, 1989. SPR 2695, Use UCB_PLANCK instead of
C		FUT_PLANCK.
C
C	R. Isaacman, R. Kummerer, November 17, 1989.  SPR 5082, Correct
C		frequency interval calculation.  Divide by NPTS.
C
C-------------------------------------------------------------------------

	IMPLICIT NONE
	INTEGER*4	NPTS
	REAL*4		SPECTRUM(NPTS)
	REAL*4		DF
	REAL*4		FNYQ
	REAL*4		FREQ
	REAL*4		TEMPK
	INTEGER*4	K
	REAL*4		PLANCK
	DF = FNYQ/NPTS
	SPECTRUM(1) = 0.
	DO K=2,NPTS
	   FREQ = (K - 1) * DF
	   SPECTRUM(K) = 1.5 * PLANCK(FREQ,TEMPK)
	ENDDO
	RETURN
	END
