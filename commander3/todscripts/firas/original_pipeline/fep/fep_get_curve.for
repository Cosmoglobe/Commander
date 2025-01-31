
 	subroutine FEP_GET_CURVE(
	1	  poly_Lun, Db_word,ndeg,poly_coeffs,istat )
C/-----------------------------------------------------------------------
C/	PROGRAM NAME:
C/	  FEP_GET_CURVE
C/
C/	PROGRAM DESCRIPTION:
C/	  This subroutine reads the polynomial coeffs. reference data 
C/        from the archive  and loads the
C/	  polynomial coefficients into output buffer  
C/
C/	AUTHOR:
C/	  Harte Wang
C/	  STX
C/	  OCTOBER 25, 1989
C/
C/	MODIFIED BY:
C/
C/
C/	CALLING SEQUENCE:
C/	  Fep_Get_Curve(Poly_Lun, Db_Word,Ndeg,Poly_Coeffs,Istat)
C/
C/	INPUT PARAMETERS:
C/	  Poly_Lun        I*4  Unit numbers for reading conversion information
C/			       : Conversion coefficients
C/	  Db_Word         C*8  Data Base word
C/
C/	OUTPUT PARAMETERS:
C/        Poly_Coeffs(6)       R*4  polynomial coffecients          
C/        Ndeg                 i*2  number of degree in polynomial
C/	  Istat                I*2  error return code
C/                               1 = no error
C/                               2 = error 
C/
C/
C/	INPUT FILES:
C/
C/	OUTPUT FILES:
C/	   
C/	INCLUDE FILES USED:
C/	
C/	
C/	
C/	  FUT_LCKEYPARS.TXT
C/	  SYS$LIBRARY:FORIOSDEF
C/
C/	SUBROUTINES USED:
C/
C/	ERROR HANDLING:
C/	  Passed back in output parameter ISTAT
C/
C/-----------------------------------------------------------------------------
C/	METHOD USED:
C/	 PDL --
C/
C/	Set return status to success.
C/       Load conversions coefficients:		 
C/          for HSKP quantity which has polynomial-type conversions     
C/	    Read FDB_LIMCUVKY with database name as key;
C/          If any error then
C/            Set return status to ERROR
C/            Write error message to screen      
C/          ELSE		
C/            Move coefficients to Poly_coef buffer
C/            **note**
C/            The coefficients in FDB_LIMCUVKY  are in ascending order;
C/            but in POLY_COEF it is
C/            in descending order because
C/            later on LIB$POLY is expecting
C/	      coefficients in descending order
C/          ENDIF
C/	RETURN
C/	END
C/------------------------------------------------------------------------

	IMPLICIT	NONE

	INCLUDE		'(FUT_LCKEYPARS)'
	INCLUDE		'SYS$LIBRARY:FORIOSDEF/NOLIST'
	INCLUDE		'($SSDEF)'

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         !
!	Passed Parameters !
! 			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

	INTEGER*4	poly_lun		!Units for configuration files
        CHARACTER*8     Db_word                 !Data base word
        Integer*2       ndeg
        Real*4          PolY_Coeffs(6)
        INTEGER*2       Istat 

        External        FEP_NOPOLYDEF
        External        FEP_LIMread
C
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			  !
!     Local Variables     !
!			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!


	byte		LCK_BUFF(LCK_RECSZ)

	INTEGER*2	 II,  NOC_EQV, IST

	REAL*4		CFF_EQV(6)
	EQUIVALENCE	(CFF_EQV,	LCK_BUFF(LCK_CFF))	! Poly. Coef.
	EQUIVALENCE	(NOC_EQV,	LCK_BUFF(LCK_NOC))	! Number of Coef.



!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set return status to success.

	istat = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									 !
!    The following  will load conversions coefficients		 !
!    for all HSKP quantities' which has polynomial-type conversions     !
!									 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	  READ (poly_lun, KEY=Db_word, IOSTAT=IST) LCK_BUFF
	  IF (IST .NE. 0) THEN
            Istat = 2 
	    IF (IST .EQ. FOR$IOS_ATTACCNON) THEN
	      Call Lib$Signal(Fep_NOPOLYDEF,%val(1),db_word)
	    ELSE
	        call lib$signal(Fep_Limread,%val(2),%val(ist),
	1	'FDB_LIMCUVKY')
	    ENDIF
	  ELSE						! GOOD READ
	    ndeg = NOC_EQV - 1
	    DO	II=1, NOC_EQV				! Move coef.
	      POLY_COEFfs(II) = CFF_EQV(NOC_EQV-II+1)	! *** NOTE ***
	    ENDDO        				! The coef. in DB file
							! are in ascending order;
							! but in POLY_COEF it is
							! in descending order because
							! later on LIB$POLY is expecting
							! coef. in descending order
	  ENDIF		! IST .NE. 0

	RETURN
	END

