	INTEGER*4 FUNCTION FEP_COUNTS_TO_OHMS ( CAL_LUN,first,
	1	                                KOUNTS_CAL, KOUNTS_GRT,
	1					NR,  ICHAN,  RESIS_GRT,
	2                                       IREC,		!***
	3                                       CO_COEFFS   )	!***

C---------------------------------------------------------------------------
C
C  WARNING:	The calling sequence for this version of the routine has
C		changed from that of Build 4.2.  This necessitates chnages
C		to the routine FEP_Convert_To_Eng as well.
C
C               Conversion of GRT telemetry counts to ohms
C	Written 23 Sept 1985 by Rich Isaacman (Applied Research Corp.)
C	Modified 28 Oct 1985 to account for ADC offset and nonlinearity
C	    "	 10 Jul 1986 to get cal resistors from database
C	    "    19 May 1989 by RBI to do double-precision solution
C
C	Modified 23-Jun-1989 by Don Stevens-Rayburn to "remember" all coef-
C		ficients that have been previously calculated.  All modifi-
C		cations are flagged with the "comment" !*** in the right
C		margin.  In addition, some redundant code was removed.
C
C	Modified 14-Jul-1989 by Don Stevens-Rayburn.  Added the "pre-computed"
C		counts to ohms coefficients 'CO_COEFFS' to the calling sequence.
C		If these are non-zero, this routine uses them to perform the
C		conversion; otherwise it does the conversion in the usual way.
C		This allows us to use an average over the entire timerange for
C		historical data while maintaining the computation of "instant-
C		aneous" conversion coefficients for real-time data.
C
C	Modified 17-Jul-1989 by Don Stevens-Rayburn to do a correct double
C		precision least squares solution.
C 
C        SPR 3921, Add the capabalities to read reference data from the
C                  reference archive.
C                  STX, Harte Wang, Feb. 10, 1990 
C----------------------------------------------------------------------------
C Subroutine name and calling sequence:
C  
C	INTEGER*4 FUNCTION FEP_COUNTS_TO_OHMS ( CAL_LUN,First,
C                                               KOUNTS_CAL, KOUNTS_GRT,
C						NR,  ICHAN,  RESIS_GRT,
C	                                        IREC,		
C	                                        CO_COEFFS   )
C
C Purpose:
C  This routine takes the telemetry counts associated with all GRTs and the
C  counts associated with all four onboard calibration resistors and returns
C  the resistance of the GRTs in ohms for a given current range. There is an
C  implicit assumption that the calibration and GRT readouts are taken in
C  the same current range.
C
C  The derivation is done via the IMSL least-squares routine LLSQF. The columns
C  of the coefficient matrix Q are arranged such that when only one cal PRT 
C  falls in the middle ADC range then the "fit" yields a simple ratio of
C  ohms per count. When 2 cal GRTs are good, the fit solves a straight line,
C  ax + b. When 3 are good it solves a parabola, and when all 4 are good it
C  derives a least-squares solution for a parabola.
C
C  Input Parameters:	KOUNTS_CAL  =  4-element integer array holding the
C				       TLM counts of the four cal resistors
C				       used for this current range.
C
C			KOUNTS_GRT  =  NR-element array of GRT TLM counts
C
C			NR	    =  Number of GRTs to convert (usually 23)
C
C                       IREC        = the Nth record.
C			ICHAN	    =  1 for A side low-current data 
C				    =  2  "  A side high-  "     "   
C				    =  3  "  B side low-  "     "   
C				    =  2  "  B side high-  "     "   
C
C			CO_COEFFS (3, 4) = real * 4 array containing the
C					conversion coefficients used to 
C					convert from counts to ohms.  These
C					coefficients are compute by the 
C					new routine FEP_FIX_COUNTS_TO_OHMS
C					and are used if we are not running
C					in real time.
C
C  Output Parameter:	RESIS_GRT   =  Array of the GRT resistances (ohms)
C
C
C
C  Internal Arrays:	RESIS_CAL   =  4 x 2 array holding the resistances
C				       of the four cal resistors for each
C				       current range. Data supplied by
C				       Dave Huff and John Sutton (728).
C
C			Q(4,3)	    =  Least-squares coefficient matrix
C
C			RHS	    =  array of "good" cal resistances for fit
C
C			COEFFS	    =  Linear or parabolic coefficients of fit
C
C			IUSE        =  array of indexes of usable cal PRTs
C
C			S1, S2	    =  scratch arrays for IMSL routine
C
C
C---------------------------------------------------------------------------
C  PDL:
C    Set the Return status to good
C 
C	If we have "pre-computed" coefficients, 
C       Then	
C         Do for each Grt Points
C          use pre-computed coefficients to
C           compute the GRT resistance 
C         Enddo 							
C	  Return						
C	End If		
C       Call  Cct_Get_Config_Idx_Tod to get
C         index and new_reference_data_flag of ref data base file
C       If new_reference_data_flag is true
C       Then
C        Call FEP_GET_Calres to get resistor calibration coefficients
C        form reference archive.
C        If bad status return
C        Then
C          Set return statusfor this function to bad
C          Return
C        endif
C       ENdif
C      
C      Determine which calibration resistors are usable and fill the least-square 
C      coefficient matrix with counts from the good GRTs only. Use the number 
C      of good GRTs to determine how many columns of the Q array to use as 
C      bases for the fit (i.e. ratio, linear, or parabolic).
C
C	RETURN
C	END
C---------------------------------------------------------------------------

	IMPLICIT	NONE

        Include		'(FUT_Params)'					!***
        Include		'(fep_data)'					!***
        Include         '(cct_get_config)'
        Include         '(cct_status_record)'
C
        record /config_status/ stat1
C
        INTEGER         *4              CAL_LUN
        Logical         *1              First
	INTEGER		*4		KOUNTS_GRT(1)
	INTEGER		*4		KOUNTS_CAL(1)
	REAL		*8		Q(4,3)				!***
	REAL		*8		RHS(4)				!***
	REAL		*8		COEFFS(3,fac_max_hskp_recs,4)	!***
	REAL		*4		CO_COEFFS ( 3, 4 )		!***
	REAL		*4		RESIS_GRT(1)
	REAL		*4		RESIS_CAL(4,4)
	REAL		*8		RES(4)				!***
	REAL		*4		S1(3)
	REAL		*4		S2(3)
C	INTEGER		*4		IUSE(4)				!***
	Integer		*4		IREC				!***
	INTEGER		*4		NCALL
	INTEGER		*4		KAL_OK
	INTEGER		*4		KK
	INTEGER		*4		KALCOUNTS
C	INTEGER		*4		NCAL				!***
	REAL		*4		TLM
	INTEGER		*4		KBASIS
	REAL		*8		TOL				!***
	INTEGER		*4		STATUS
	INTEGER		*4		ICHAN
	INTEGER		*4		IER
	INTEGER		*4		Index
	INTEGER		*4		J
	INTEGER		*4		M
	INTEGER		*4		NR
	INTEGER		*4		FEP_GET_CALRES
        INTEGER         *4              CCT_GET_CONFIG_IDX_TOD
        LOGICAL         *1              NEW_CAL
	EXTERNAL	FEP_LLSQFERR
	EXTERNAL	FEP_NORMAL
	EXTERNAL        FEP_GETCONFIGERR
        DATA NCALL /0/


	FEP_COUNTS_TO_OHMS = %LOC(FEP_NORMAL)

C	If we have "pre-computed" coefficients, use them!		!***
									!***
	If ( co_coeffs ( 1, 1 ) .ne. 0.0 ) then				!***
	   Do J = 1, nr							!***
	      tlm = kounts_grt ( j )					!***
	      resis_grt( j ) = co_coeffs ( 2, ichan )			!***
	1                    + co_coeffs ( 1, ichan ) * tlm		!***
	1                    + co_coeffs ( 3, ichan ) * tlm * tlm	!***
	   End Do							!***
	   Return							!***
	End If								!***

C
C
C  Get Cal resistances from current database on first pass only
C
C	IF (NCALL .EQ. 0) THEN
C	   NCALL = 1
           status= cct_get_config_idx_tod(timetags(irec).time,
	1                1,cal_lun,index,new_cal,stat1)
           If (.not. status) then
             Call lib$signal(fep_getconfigerr,%val(1),%val(Status))
             Fep_Counts_to_ohms = status
             return
           Endif
           if (new_cal .or. first) then
	     STATUS = FEP_GET_CALRES (cal_lun,RESIS_CAL)
	     IF (STATUS .NE. %LOC(FEP_NORMAL)) THEN
	      FEP_COUNTS_TO_OHMS = STATUS
	      RETURN
	     END IF
           ENDIF
C	ENDIF
C
C  Determine which cal resistors are usable and fill the least-square 
C  coefficient matrix with counts from the good GRTs only. Use the number 
C  of good GRTs to determine how many columns of the Q array to use as bases 
C  for the fit (i.e. ratio, linear, or parabolic).
C
C

	IF ( COEFFS ( 1, IREC, ICHAN ) .EQ. 0. ) THEN
	   KAL_OK = 0
	   DO J=1,4
	      KALCOUNTS = KOUNTS_CAL(J)
	      IF (KALCOUNTS.GT.100 .AND. KALCOUNTS.LT.16000) THEN
	         KAL_OK = KAL_OK + 1
	         TLM = KOUNTS_CAL(j)					!***
	         Q(KAL_OK,1) = TLM					!***
	         Q(KAL_OK,2) = 1.					!***
	         Q(KAL_OK,3) = TLM * TLM				!***
	         RHS(KAL_OK) = RESIS_CAL(j,ICHAN)			!***
	      ENDIF
	   ENDDO
	   IF (KAL_OK .EQ. 0) THEN
	      DO M=1,NR
	         RESIS_GRT(M) = 0.
	      ENDDO
	      RETURN
	   ENDIF

	   TOL = 0.							  !***	
	   CALL DLSQRR (KAL_OK, 3, Q, 4, RHS, TOL,			  !***
	1		COEFFS ( 1, IREC, ICHAN ) , RES, KBASIS)	  !***

	ENDIF

C
C  Apply fit to convert all GRT counts to ohms
C
	DO J=1,NR
	   TLM = KOUNTS_GRT(J)
	   RESIS_GRT( J ) = COEFFS( 2, IREC, ICHAN )			!***
	2                 + COEFFS( 1, IREC, ICHAN )*TLM		!***
	3                 + COEFFS( 3, IREC, ICHAN )*TLM*TLM		!***
	ENDDO

	RETURN
	END
