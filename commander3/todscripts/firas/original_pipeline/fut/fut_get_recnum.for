	INTEGER*4 FUNCTION FUT_GET_RECNUM (FAKEIT, MTMSPEED, ICHAN,
     .					MICRO_MODE, IADDS_GRP, IRECNO)

c------------------------------------------------------------------------------
c
c  Function:
c
c  This routine takes the numerical values associated with the appropriate
c  instrument status bits and returns the corresponding record number in 
c  the direct-access Electronics Transfer Function file. That record holds
c  the transfer function associated with the given instrument configuration.
c
c  Written 29 April 1984 by Rich Isaacman (Applied Research Corp.)
c
c  Input parameters:	FAKEIT		: Fakeit pulse status. 0=off, 1=on
c
c			MTMSPEED	: 0=slow scan, 1=fast scan
c
c			ICHAN		: Channel number (1=RH, 2=RL, etc.)
c
c			MICRO_MODE	: uProc science mode ("block type" 0-4)
c
c			IADDS_GRP	: uP data compression (adds per group)
c
c  Output parameter:	IRECNO		: Desired record number in xfcn file
c
c  Subroutines called: None
c
c  Include files: None
c
c------------------------------------------------------------------------------
c
c Changes:
c
c	R. Isaacman, 21 Jan 1987 to account for micro science mode 4.
c
c	R. Kummerer, 20 Apr 1987 to move from FFC to FUT.
c
c	R. Isaacman, 28 Oct 1987 to reflect change in transfer function 
c		file structure, with different records for all 4 channels
c
c------------------------------------------------------------------------------

	IMPLICIT NONE
	INCLUDE  '(FUT_ERROR)'
	INTEGER * 4	MTMSPEED
	INTEGER * 4	FAKEIT
	INTEGER * 4	IRECNO
	INTEGER * 4	ICHAN
	INTEGER * 4	IDIG_FLTR
	INTEGER * 4	MICRO_MODE
	INTEGER * 4	IADDS_GRP

	EXTERNAL	FUT_NORMAL
	EXTERNAL	FUT_ETFRECNO

	IRECNO = 96 * (1 - MTMSPEED)
	IF (FAKEIT .EQ. 1) IRECNO = 192
	IRECNO = IRECNO + 24*(ICHAN - 1)
	IF (MICRO_MODE.EQ.2 .OR. MICRO_MODE.EQ.4) THEN
	   IDIG_FLTR = 0				!Digital filter on
	ELSE
	   IDIG_FLTR = 1				!Digital filter off
	ENDIF
	IRECNO = IRECNO + 12*IDIG_FLTR
	IRECNO = IRECNO + MAX0 (IADDS_GRP,1)
	IF(IRECNO .LE. 0)THEN
	   FUT_GET_RECNUM = %LOC(FUT_ETFRECNO)
	   CALL LIB$SIGNAL(FUT_ETFRECNO,%VAL(1),%VAL(IRECNO))
	ELSE IF(IRECNO .GT. 288)THEN	
	   FUT_GET_RECNUM = %LOC(FUT_ETFRECNO)
	   CALL LIB$SIGNAL(FUT_ETFRECNO,%VAL(1),%VAL(IRECNO))
	ELSE
	   FUT_GET_RECNUM = %LOC(FUT_NORMAL)
	END IF
	RETURN
	END
