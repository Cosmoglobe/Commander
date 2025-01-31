	INTEGER*4 FUNCTION FNT_DISPLAY_NOISE ( CHAN,SPEC_REC,FAKE_IT,DIF_FLAG,
	1                              NUMBER, SIZE, LUN, CINDEX, CONFIG,
	2                              NEW_SEGMENT, STAT)

C------------------------------------------------------------------------
C    PURPOSE: 
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            STX
C            May 7, 1987
C
C	     Adpated from code written by W.K.Young
C
C    INVOCATION: 
C
C    INPUT PARAMETERS:
C
C    OUTPUT PARAMETERS: 
C
C    SUBROUTINES CALLED: 
C
C    COMMON VARIABLES USED: 
C
C    INCLUDE FILES: 
C
C----------------------------------------------------------------------
C
C  Changes:
C
C	SER 3493, R. Kummerer, August 18, 1989.  Use new PLT graphics
C		to display IFGs and spectra.
C
C	SER 5779, S. Alexander, March 5, 1990.  Add capability to pass in
C		a PLT command file.
C       SER 6626, H. Wang. 4/12/90, FNT difference in spectra needed for 
C                noise monitoring.            
C
C	SPR 7729, S. Alexander, November 27, 1990.  Add OFFSET argument
C		to FUT_SETXAX call.
C
C	SPR 9583, S. Alexander, August 3, 1992.  Change calling list to
C		  FUT_SETXAX.
C
c       SPR 9790, Update FNT to use new FEX_Nyquist and FUT_Setxax.
c                 Nilo G. Gonzales, Hughes STX, 1992 August 4.
C----------------------------------------------------------------------

	IMPLICIT NONE

	INCLUDE '(FUT_PARAMS)'
	INCLUDE '(FNT_INVOC)'
	INCLUDE '(CCT_GET_CONFIG)'
C 
C Functions and Subroutines
C 
	INTEGER   *  4 FUT_SETXAX, FUT_PLOT_TITLE, FUT_PLOT
	INTEGER   *  4 CCT_GET_CONFIG_TOD
C
C Records
C
	DICTIONARY 'FNT_NOISE'
	RECORD /FNT_NOISE/ SPEC_REC

	DICTIONARY 'FEX_NYQUIST'
	STRUCTURE /CONFIG_DATA/
	   RECORD /FEX_NYQUIST/FEX_NYQUIST
	ENDSTRUCTURE

	RECORD /CONFIG_DATA/ CONFIG
	RECORD /CONFIG_STATUS/ STAT(1)

	EXTERNAL FNT_GETCONFIGERR

	INTEGER   *  4 STATUS
	INTEGER   *  4 I
	CHARACTER * 16 ANSWER

	INTEGER   *  4 NGROUP
	INTEGER   *  4 SWEEPS
	INTEGER   *  4 MTM_SPEED
	INTEGER   *  4 MTM_LENGTH
	INTEGER   *  4 CHAN
	INTEGER   *  4 UPMODE
	INTEGER   *  4 XCAL_POS
	INTEGER   *  4 FAKE_IT
	INTEGER   *  4 GAIN
        Logical   *  1 Dif_flag
	CHARACTER * 60 LABEL
	CHARACTER * 60 PLOT_LABEL
	CHARACTER *100 TITLE(3)
        INTEGER   *  4 FIRSTSAMP /1/

	INTEGER   *  4 NUMBER   	 ! number of data sets
	INTEGER   *  4 SIZE              ! size of data sets in bytes
	INTEGER   *  4 LUN(1)		 ! logical unit numbers
	INTEGER   *  4 CINDEX(1)         ! initial cache pointers
	LOGICAL   *  1 NEW_SEGMENT(1)	 ! flag for new segments
	INTEGER   *  4 JTIME(2)
	CHARACTER * 14 GMTFIRST   
	INTEGER   *  4 ZP
	REAL      *  4 STARTX
	REAL      *  4 SPACEX
	INTEGER   *  4 INTERACTIVE
	CHARACTER * 32 PLOT_DEVICE

	COMPLEX   *  8 DSPEC(1024)

C
C Initialize for plotting.
C
	ZP = FAC_PRESENT
	INTERACTIVE = FTC_INTERACTIVE
	PLOT_DEVICE = FTC_PLOTS_DEVICE

	CHAN = SPEC_REC.CHAN.CHAN_ID
	NGROUP = SPEC_REC.CHAN.NGROUP
        SWEEPS = SPEC_REC.CHAN.SWEEPS
	MTM_SPEED = SPEC_REC.CHAN.MTM_SPEED
	MTM_LENGTH = SPEC_REC.CHAN.MTM_LENGTH
        UPMODE = SPEC_REC.CHAN.SCI_MODE
	XCAL_POS = SPEC_REC.COA_HEAD.XCAL_POS
	GAIN = FAC_GAINS(SPEC_REC.CHAN.GAIN)
	LABEL = SPEC_REC.COA_HEAD.LABEL
	GMTFIRST = SPEC_REC.COA_HEAD.FIRST_GMT

	CALL CT_GMT_TO_BINARY ( GMTFIRST, JTIME)
c
c Access FEX_Nquist reference dataset using CCT_Get_Config. 
c
	STATUS = CCT_GET_CONFIG_TOD ( JTIME, NUMBER, SIZE, LUN, CINDEX,
	1                           CONFIG, NEW_SEGMENT, STAT )

	IF (.NOT. STATUS) THEN
	    CALL LIB$SIGNAL(FNT_GETCONFIGERR,%VAL(1),%VAL(STATUS))
	ENDIF

	CALL FUT_SETXAX(-1,FAKE_IT,MTM_SPEED,MTM_LENGTH,CHAN,
     .                  NGROUP,UPMODE,FIRSTSAMP,CONFIG,STARTX,SPACEX)

C
C Query user for displaying noise data.
C
	IF (INTERACTIVE .EQ. FAC_PRESENT) THEN

	   ANSWER = 'Y'

	   DO WHILE (ANSWER(1:1) .NE. 'Q')

              TYPE 10
10            FORMAT(/,3X,'Choose one of the following:',/,
	1		5x,'Noise spectrum',/,
	2		5x,'Skip n spectra',/,
	3		5x,'Quit displaying plots',/,
	4		5x,'Hit return for next channel => ',$)
	      ACCEPT 20, ANSWER
20	      FORMAT(A)

	      CALL STR$UPCASE ( ANSWER, ANSWER )

	      IF (ANSWER(1:1) .EQ. 'N') THEN
C
C Plot spectra.
C
                  DO I=1,257
                        DSPEC(I) = SPEC_REC.CHAN.SPEC(I)
                  END DO
                  If (Dif_flag) then
 		     PLOT_LABEL = 'Noise Spectrum (Difference)'
                  else
		    PLOT_LABEL = 'Noise Spectrum'
                  Endif
		  CALL FUT_PLOT_TITLE(LABEL,NGROUP,SWEEPS,MTM_SPEED,MTM_LENGTH,
	1			      CHAN,UPMODE,XCAL_POS,FAKE_IT,GAIN,
	2			      PLOT_LABEL,TITLE)
		  CALL FUT_PLOT(DSPEC,257,STARTX,SPACEX,TITLE,
	1			'Frequency (Hz)','Volt/Hz**0.5',
	2			ZP,INTERACTIVE,PLOT_DEVICE,
	3			FTC_PLT_COM,FTC_PLT_COM_FILE)

	      ELSE IF (ANSWER(1:4) .EQ. '    ') THEN
C
C On to next plot.
C
                  ANSWER = 'Q'

	      ELSE IF (ANSWER(1:1) .EQ. 'Q') THEN
C
C Set flag not to display any more plots.
C
		  FTC_DISPLAY = .FALSE.

	      ELSE IF (ANSWER(1:1) .EQ. 'S') THEN
C
C Jump over N spectra.
C
		  TYPE 35
35                FORMAT(3X,' Enter number of spectra to jump => ',$)
	          ACCEPT *, FTC_JUMP
		  ANSWER = 'Q'

	      ELSE
C
C Reply not recognized; try again.
C
                  TYPE 60
60      	  FORMAT(3x,'Unrecognized input -- please try again')

	      END IF

	   END DO

	ELSE

C
C Plot spectra.
C
	   DO I=1,257
              DSPEC(I) = SPEC_REC.CHAN.SPEC(I)
           END DO

	   PLOT_LABEL = 'Noise Spectrum'

	   CALL FUT_PLOT_TITLE(LABEL,NGROUP,SWEEPS,MTM_SPEED,MTM_LENGTH,
	1		       CHAN,UPMODE,XCAL_POS,FAKE_IT,GAIN,
	2		       PLOT_LABEL,TITLE)
	   CALL FUT_PLOT(DSPEC,257,STARTX,SPACEX,TITLE,
	1		 'Frequency (Hz)','Volt/Hz**0.5',
	2		 ZP,INTERACTIVE,PLOT_DEVICE,
	3		 FTC_PLT_COM,FTC_PLT_COM_FILE)

	END IF

	FNT_DISPLAY_NOISE = FAC_NORMAL

	RETURN
	END
