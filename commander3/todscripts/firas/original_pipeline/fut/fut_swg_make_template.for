	INTEGER*4 FUNCTION FUT_SWG_MAKE_TEMPLATE (AIFGS,NIFG, 
     .                                 TEMPLATE,SIGMAS,AVG_SIGMA)
c
c  This routine takes an array of interferograms (e.g. the ensemble retrieved
c  by Category Select) and forms the the template used in the Consistency
c  Checker and Deglitcher. The template is essentially the "median 
c  interferogram" of the ensemble. In addition to making the template, it also
c  returns (1) an estimate of the rms noise in each interferogram, (2) the 
c  median noise of the ensemble, and (3) the original ifgs with dither 
c  and template subtracted. The original interferograms -- less dither --
c  can be reconstructed by adding the template back in.
c
c  Input parameters:	AIFGS	   = Real*4 array of interferograms, 512, x NIFG
c				     where NIFG is the number of ifgs in the 
c				     group. On output, AIFGS contains the same
c				     interferograms, but with the dither and
c				     template subtracted off.
c
c			NIFG      =  Integer*4 number of interferograms in the
c				     group.
c
c  Output parameters:   TEMPLATE   = Real*4 array of length 512 holding the 
c				     median ifg of the ensemble.
c
c			SIGMAS     = Real*4 array of length NIFG holding the
c				     rms noise (in counts) for each ifg.
c
c			AVG_SIGMA  = Real*4 median of the elements in array 
c				     SIGMAS.
c
c  Internal arrays:	BUFFER(512)  used for buffering sorts.
c
c  Subroutines called:  VSRTA (IMSL algebraic sort)
c
c  Written 29 January 1986 by Rich Isaacman (Applied Research Corp.)
c
C  CHANGES:
C
C   9/19/86   J.DURACHTA    ERROR MESSAGES ADDED. ERROR TRAPPING ALTERED.
C

        IMPLICIT NONE

C       ERROR MESSAGES

        external fut_sig2sm

        INTEGER*4 NIFG
        INTEGER*4 NROW,IPT,IMID      !COUNTERS
        INTEGER*4 NCENTER
        INTEGER*4 NQUART
        INTEGER*4 NMIDAV

        REAL*4 DITHER
        REAL*4 AVG
        REAL*4 AVGSQ
        REAL*4 DATAPT
        REAL*4 SIGMA
        REAL*4 AIFGS(512,NIFG)
	REAL*4 TEMPLATE(512)
	REAL*4 SIGMAS(NIFG)
	REAL*4 AVG_SIGMA
	REAL*4 BUFFER(512)
c
c  Subtract dither by subtracting median value from each ifg in group.
c
	DO NROW=1,NIFG
	  DO IPT=1,512
	    BUFFER(IPT) = AIFGS(IPT,NROW)
	  ENDDO
	  CALL SVRGN(512,BUFFER,BUFFER)
	  DITHER = (BUFFER(256) + BUFFER(257))/2.
	  DO IPT=1,512
	    AIFGS(IPT,NROW) = AIFGS(IPT,NROW) - DITHER
	  ENDDO
	ENDDO
c
c  Form template by sorting along columns to find "midaveraged ifg" at each pt
c
	NQUART = NIFG/4
	NMIDAV = NIFG - 2*NQUART
	DO IPT=1,512
	   DO NROW=1,NIFG
	      BUFFER(NROW) = AIFGS(IPT,NROW)
	   ENDDO
	   CALL SVRGN(NIFG,BUFFER,BUFFER)
	   TEMPLATE(IPT) = 0.
	   DO IMID = NQUART+1,NQUART+NMIDAV
	      TEMPLATE(IPT) = TEMPLATE(IPT) + BUFFER(IMID)
	   ENDDO
	   TEMPLATE(IPT) = TEMPLATE(IPT)/NMIDAV
	ENDDO
c
c  Estimate noise of each ifg to be the rms value after template subtraction.
c
	DO NROW=1,NIFG
	  AVG = 0.
	  AVGSQ = 0.
	  DO IPT=1,512
	    DATAPT = AIFGS(IPT,NROW) - TEMPLATE(IPT)
	    AIFGS(IPT,NROW) = DATAPT
	    AVG = AVG + DATAPT
	    AVGSQ = AVGSQ + DATAPT*DATAPT
	  ENDDO
	  AVG = AVG/512.
	  SIGMA = SQRT (AVGSQ/512. - AVG*AVG)
	  SIGMAS(NROW) = AMAX1(SIGMA, 0.5)			! Set sigma > .5
	  BUFFER(NROW) = SIGMAS(NROW)

c     _SWG_MAKE_TEMPLATE: SIGMAS('!UL') IS < .5 . IT HAS BEEN SET TO .5 

          IF(SIGMA .LT. .5)CALL LIB$SIGNAL(FUT_SIG2SM,%VAL(1),%VAL(NROW))

	ENDDO
c
c  Take median value of all sigmas to be characteristic noise of the ensemble
c
	CALL SVRGN(NIFG,BUFFER,BUFFER)
	NCENTER = NIFG/2 + 1
	AVG_SIGMA = BUFFER(NCENTER)
	RETURN
	END
