
	SUBROUTINE FTB_GET_RAI_TIM (WRITE_LUN, NRECS, BUFF)
C/
C/	PROGRAM NAME:
C/	  FTB_GET_RAI_TIM
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine will print the GMT time tags for a given
C/	  Index buffer.
C/
C/	AUTHOR:
C/	  E.FUNG
C/	  GSFC
C/	  NOV. 2, 1985
C/
C/	MODIFIED BY:
C/	  N/A
C/
C/	CALLING SEQUENCE:
C/	  CALL GET_RAI_TIM (WRITE_LUN, NRECS, BUFF)
C/
C/	INPUT PARAMETERS:
C/	  WRITE_LUN	I*2	Logical unit number for the list file.
C/	  NRECS		I*4	The count of the record processed so far (counted
C/				in the calling routine.)
C/	  BUFF(5632)	BYTE	Buffer holding the Science record.
C/
C/	OUTPUT PARAMETERS:
C/	  NONE
C/
C/	INPUT FILE:
C/	  NONE
C/
C/	OUTPUT FILE:
C/	  File for the time tags (at logical unit WRITE_LUN).
C/
C/	INCLUDE FILES USED:
C/	  NONE
C/
C/	SUBROUTINES CALLED:
C/	  TBD
C/
C/	ERROR HANDLING
C/	  NONE
C/
C/	METHOD USED:
C/	 Trivial
C/
	IMPLICIT	NONE

	BYTE		BUFF(512)

	INTEGER*2	WRITE_LUN

	INTEGER*4	NRECS

!
!     Local variables
!

	BYTE		END_GMT(14), START_GMT(14)

	INTEGER*2	I, J, K, L

	INTEGER*4	END_BNTM(2), START_BNTM(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                          !
!     Code begins here     !
!                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	CALL LIB$MOVC3 (8, BUFF(29), START_BNTM)
	CALL LIB$MOVC3 (8, BUFF(37), END_BNTM)

	CALL BINARY_TO_GMT (START_BNTM, START_GMT)
	CALL BINARY_TO_GMT (END_BNTM, END_GMT)

	WRITE (WRITE_LUN, 6001)		nrecs,(BUFF(I),I=1,14),
	1				START_BNTM(2), START_BNTM(1),
	1				(START_GMT(J),J=1,14),
	1				(BUFF(K),K=15, 28),
	1				END_BNTM(2), END_BNTM(1),
	1				(END_GMT(L),L=1,14)

6001	FORMAT (/' # RAI records read in: ', I7 /
	1	'        Start GMT of Index record:   ', 
	1	2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1,
	1	2X, '[ ', Z8.8, 1X, Z8.8, ' = ',
	1	2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, ' ]'/
	1	'        End GMT of Index record:     ', 
	1	2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1,
	1	2X, '[ ', Z8.8, 1X, Z8.8, ' = ',
	1	2A1, '-', 3A1, '-', 3(2A1, '-'), 3A1, ' ]')

	RETURN
	END 
