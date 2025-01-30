
 	INTEGER*4 FUNCTION fep_get_calres( res_lun,resis_cal)
C/----------------------------------------------------------------------
C/	PROGRAM NAME:
C/	  Fep_Get_Calres
C/
C/	PROGRAM DESCRIPTION:
C/	  This subroutine reads the GRT calibration resistor
C/          reference data from the archive and loads the
C/	  resistor information to the buffer
C/
C/
C/	AUTHOR:
C/	  Harte Wang
C/	  STX
C/	  Feb. 9, 1989
C/
C/	MODIFIED BY:
C/
C/
C/      INPUT PARAMETERS:
C/        RES_LUN
C/	OUTPUT PARAMETERS:
C/	  Resis_cal
C/
C/
C/
C/	INPUT FILES:
C/
C/	OUTPUT FILES:
C/	 
C/	 
C/
C/	INCLUDE FILES USED:
C/	  SYS$LIBRARY:FORIOSDEF
C/
C/	SUBROUTINES USED:
C/	  LIB$MOVC3 (from system library)
C/
C/	ERROR HANDLING:
C/
C/	METHOD USED:
C/	PDL --
C/       Set this function return status to success
C/        Do for all 4 calibration sets
C/	    Do for all 4 calibration resistors
C/	      Read configured file FDB_GRTCAL using calibration resistor's 
C/              database name as the key
C/            If no error 
C/            Then
C/   	        Get the conversion coefficients from the reference archive 
C/                        dataset. 
C/            Else
C/              Set this function return status to error
C/              Write Error message to screen
C/              Return
C/            Endif
C/	    Enddo;
C/	  Enddo;
C/	  Return
C/	  End
C/-------------------------------------------------------------------------

	IMPLICIT	NONE
	INCLUDE		'SYS$LIBRARY:FORIOSDEF/NOLIST'
	INCLUDE		'($SSDEF)'

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         !
!	Passed Parameters !
! 			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

	INTEGER*4	res_lun        		!Units for configuration files
        REAL*4          resis_cal(4,4)
C
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			  !
!     Local Variables     !
!			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!


	INTEGER*2	I, IST

	LOGICAL*1	BUFF20(20)  	! Read buffer

	CHARACTER*8	CALRES_NAME(4, 4)
	1		/'FACR1TLO', 'FACR2TLO',
	1		'FACR3TLO', 'FACR4TLO',
	1		'FACR1THI', 'FACR2THI',
	1		'FACR3THI', 'FACR4THI',
	1		'FBCR1TLO', 'FBCR2TLO',
	1		'FBCR3TLO', 'FBCR4TLO',
	1		'FBCR1THI', 'FBCR2THI',
	1		'FBCR3THI', 'FBCR4THI' /

	INTEGER*2	CALSET

	EXTERNAL   	FEP_NORMAL
	EXTERNAL	FEP_Cresread
        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!	Read the calibrator resistor values from the configuration dataset.  !
!	 
!	  CALRES(4,4)           R*4     Calibrator resistor values           ! 
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	DO CALSET = 1, 4	 
	  DO	I = 1, 4
	    READ (res_lun, KEY=CALRES_NAME(I, CALSET), IOSTAT=IST)
	1	BUFF20
	    IF (IST .NE. 0) THEN
	      CALL LIB$SIGNAL(FEP_Cresread,%VAL(2),%VAL(IST),
	1		'FDB_GRTCAL')
	      Fep_Get_Calres=%loc(Fep_Cresread)
              Return
	    ELSE
	      CALL LIB$MOVC3(4, BUFF20(15), Resis_cal(I, CALSET))
	    ENDIF
	  ENDDO
	ENDDO		! CALSET = 1, 4


C	Set function to return status

	fep_get_calres = %loc(FEP_NORMAL)
	RETURN
	END

