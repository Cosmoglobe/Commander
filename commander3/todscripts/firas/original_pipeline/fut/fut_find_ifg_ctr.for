	INTEGER*4 FUNCTION FUT_FIND_IFG_CTR(COAD_IFG, NPTS, NGROUP, MTM_SPEED,
     .						ICTR)

c-----------------------------------------------------------------------------
c
c  This subroutine is to be used in conjunction with subroutine SYMCENTER to
c  find the nominal zero-path difference point of an interferogram. It returns
c  an error flag when the interferogram center cannot be located unambiguously
c  due to low signal-to-noise ratio.
c
c  Written August 1985 by Rich Isaacman (Applied Research Corp.)
c
c  Input parameters:	COAD_IFG  = Input interferogram array (R*4)
c			NPTS      = Number of points in interferogram (I*4)
c			NGROUP    = Adds per group (I*4)
c			MTM_SPEED = MTM speed (1 = slow, 2 = fast) (I*4)
c
c  Output parameters:	ICTR   =   Index of interferogram peak (I*4)
c			STATUS =   FUT_NORMAL if all went well
c				   FUT_DEFPEAK if normal peak finding failed but
c					a default was successfully determined
c				   FUT_INVDEFPEAK if no peak was successfully
c					determined
c
c  Restrictions:	This routine will fail if the interferogram peak falls
c			at ICTR < 11 or ICTR > NPTS-10. In normal scan modes
c			this should not occur.
c
c-----------------------------------------------------------------------------
c  Changes:
c
c   9/19/86    J.DURACHTA     ERROR MESSAGES ADDED.ERROR TRAPPING ALTERED.
c
c   8/18/87    R.KUMMERER     MOVE FROM FCI INTO FUT.
c
c   10/23/87   R.ISACCMAN     CHANGE EPSILON TO PREVENT RMS FROM BEING 0
c                             FROM 0.01 TO 1.0E-10.
c
c   1/9/89     R. Kummerer    SPR 3109, Return a default IFG peak if normal
c				peak finding fails.
c
c-----------------------------------------------------------------------------

        IMPLICIT NONE

        INTEGER*4 J           !A COUNTER
        INTEGER*4 NPTS
        INTEGER*4 ICTR
        INTEGER*4 R_STATUS    !RETURN STATUS
        INTEGER*4 MTM_SPEED
	INTEGER*4 NGROUP

        INTEGER*4 FUT_DEFAULT_PEAK

        REAL*4 COAD_IFG
	DIMENSION COAD_IFG(1)
        REAL*4 RMS, SN
        REAL*4 SYMMAX, PARMAV, PARMSQ, SYMPARM

        EXTERNAL FUT_NORMAL
        EXTERNAL FUT_ICTRBDA
        EXTERNAL FUT_ICTRBDB
        EXTERNAL FUT_NOSYMCTR
        EXTERNAL FUT_DEFPEAK
        EXTERNAL FUT_INVDEFPEAK

	SYMMAX = -10000.
	PARMAV = 0.
	PARMSQ = 0.

        R_STATUS = %LOC(FUT_NORMAL)
c
c  Get symmetry parameter at every point (from FUT_SYMCENTER),
c  and accumulate its mean and mean-square.
c
	J=11

	DO WHILE (J .LE. NPTS-10 .AND. R_STATUS .EQ. %LOC(FUT_NORMAL))

	   CALL FUT_SYMCENTER (COAD_IFG,J,SYMPARM)

	   IF (SYMPARM .GT. SYMMAX) THEN

	      SYMMAX = SYMPARM
	      ICTR = J

              IF (ICTR .LT. 11 .OR. ICTR .GT. NPTS-10) THEN

c _FIND_IFG_CTR ERROR: INDEX OF IFG PEAK IS OUT OF BOUNDS
c ICTR='!UL'.NOMINAL VALUES ARE 11 .LT. ICTR .LT. NPTS - 10.

                CALL LIB$SIGNAL(FUT_ICTRBDA,%VAL(0),
     .                           FUT_ICTRBDB,%VAL(1),%VAL(ICTR))
                R_STATUS = %LOC(FUT_ICTRBDA)
              END IF

	   ENDIF

	   PARMAV = PARMAV + SYMPARM
	   PARMSQ = PARMSQ + SYMPARM*SYMPARM

           J=J+1

	ENDDO
c
c  Calculate rms of symmetry parameter SYMPARM; using this, calculate the 
c  signal-to-noise ratio (SN) at the peak value of SYMPARM. If SN < 5.5 the
c  interferogram is too noisy to clearly locate the peak.
c
	IF (R_STATUS .EQ. %LOC(FUT_NORMAL)) THEN

	   PARMAV = PARMAV/(NPTS - 20)
	   PARMSQ = PARMSQ/(NPTS - 20)
	   RMS = SQRT (PARMSQ - PARMAV*PARMAV) + 1.0E-10	!Avoid RMS=0

	   SN = (SYMMAX - PARMAV)/RMS

	   IF(SN .LE. 5.5)THEN
             R_STATUS = %LOC(FUT_NOSYMCTR)

c _FIND_IFG_CTR: IFG TOO NOISY TO LOCATE PEAK.

             CALL LIB$SIGNAL(FUT_NOSYMCTR)
           END IF

	END IF

c
c Supply the default IFG peak position if normal processing fails.
c
	IF (R_STATUS .NE. %LOC(FUT_NORMAL)) THEN

	   CALL LIB$SIGNAL(FUT_DEFPEAK)

	   R_STATUS = FUT_DEFAULT_PEAK(MTM_SPEED-1,NGROUP,ICTR)

	   IF (R_STATUS .EQ. %LOC(FUT_NORMAL)) THEN
	      R_STATUS = %LOC(FUT_DEFPEAK)
	   ELSE
	      CALL LIB$SIGNAL(FUT_INVDEFPEAK)
	      R_STATUS = %LOC(FUT_INVDEFPEAK)
	   END IF

	END IF

	FUT_FIND_IFG_CTR = R_STATUS

	RETURN
	END
