	INTEGER*4 FUNCTION FDQ_GET_HSK_FLDS (HSK_ANALOG, 
	1	IPDU_RELAY, IDX_REC)
C/
C/	PROGRAM NAME:
C/	  FDQ_GET_HSK_FLDS 
C/
C/	PROGRAM DESCRIPTION:
C/	  This routines fill the index fields that can only be picked up
C/	  from the HSK record (analog counts and the IPDU individual
C/	  statuses).
C/
C/	AUTHOR:
C/	  EDWIN FUNG
C/	  GSFC
C/	  MARCH 30, 1986
C/
C/	MODIFIED BY:
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 14, 1988
C/	  REASON: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal.
C/
C/	CALLING SEQUENCE:
C/	  Status = FDQ_GET_HSK_FLDS  (HSK_ANALOG, IPDU_RELAY, IDX_REC)
C/
C/	INPUT PARAMETERS:
C/	  HSK_ANALOG(52)	BYTE	HSK analog values in counts
C/	  IPDU_RELAY(8)		BYTE	The IPDU power relay status.  The
C/					individual status bits are taken
C/					from these 8 bytes as well.
C/
C/	OUTPUT PARAMETERS:
C/	  IDX_REC(512)	BYTE	Current Index record buffer
C/
C/	INPUT/OUTPUT FILES:
C/	  NONE
C/
C/	SUBROUTINES CALLED:
C/	  LIB$MOVC3
C/
C/	INCLUDE FILES USED:
C/	  FIR_IDX_PARMS.INC
C/
C/	ERROR HANDLING:
C/	  None
C/
C/	METHOD USED:
C/	  Trivial
C/

	IMPLICIT	NONE

        dictionary 'fdq_idx'
        record /fdq_idx/   idx_rec

	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR

	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status

	BYTE		HSK_ANALOG(52), IPDU_RELAY(0:7)

	BYTE		TEMP_BYTE_EQV

	INTEGER*2	I, IPDU_RELAY_I2(0:7), TEMP

	EQUIVALENCE	(TEMP_BYTE_EQV,	TEMP)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               !
!       Code Starts Here        !
!                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set return status to success.

	RETSTAT = SUCCESS

        do i = 1,14
          idx_rec.idx_analog.temp(i) = hsk_analog(i)
        end do

        do i = 1,32
          idx_rec.idx_analog.v_and_i(i) = hsk_analog(i+14)
        end do

        idx_rec.idx_analog.mtm_cal_motor(1) = hsk_analog(47)
        idx_rec.idx_analog.mtm_cal_motor(2) = hsk_analog(48)

!
!     The last 4 bytes in HSK_ANALOG are the bolometer bias voltages.
!     They are in the channel specific fields of the IDX record
!
        do i = 1,4
          idx_rec.chan(i).bol_readout_volt = hsk_analog(i+48)
        end do

!
!     Move the IPDU status words and pick out the individual IPDU
!     status bits
!

        call lib$movc3(8,ipdu_relay,idx_rec.ipdu_relay.group4)

	DO	I = 0,7 		          ! It is necessary to move
	  IPDU_RELAY_I2(I) = ipdu_relay(i)        ! the statuses to INTEGER
                                                  ! in order to use instrinsic
	ENDDO				          ! routines IBITS and MVBITS

        do i = 1,2                               !loop over side a & b

          TEMP = 0
          CALL MVBITS(IPDU_RELAY_I2(i-1), 6, 1, TEMP, 1)
          CALL MVBITS(IPDU_RELAY_I2(i+3), 6, 1, TEMP, 0)
          IDX_REC.ipdu_stat(i).group5(1) = TEMP_BYTE_EQV

          TEMP = 0
          CALL MVBITS(IPDU_RELAY_I2(i-1), 4, 1, TEMP, 1)
          CALL MVBITS(IPDU_RELAY_I2(i+1), 4, 1, TEMP, 0)
          IDX_REC.ipdu_stat(i).group5(2) = TEMP_BYTE_EQV

          TEMP = 0
          CALL MVBITS(IPDU_RELAY_I2(i-1), 5, 1, TEMP, 1)
          CALL MVBITS(IPDU_RELAY_I2(i+1), 5, 1, TEMP, 0)
          IDX_REC.ipdu_stat(i).group5(3) = TEMP_BYTE_EQV

          IDX_REC.ipdu_stat(i).group5(4) = IBITS(IPDU_RELAY_I2(i-1), 1, 1)
          IDX_REC.ipdu_stat(i).group5(5) = IBITS(IPDU_RELAY_I2(i-1), 0, 1)

          TEMP = 0
          CALL MVBITS(IPDU_RELAY_I2(i+1), 3, 1, TEMP, 1)
          CALL MVBITS(IPDU_RELAY_I2(i+3), 3, 1, TEMP, 0)
          IDX_REC.ipdu_stat(i).group5(6) = TEMP_BYTE_EQV

          TEMP = 0
          CALL MVBITS(IPDU_RELAY_I2(i+1), 2, 1, TEMP, 1)
          CALL MVBITS(IPDU_RELAY_I2(i+3), 2, 1, TEMP, 0)
          IDX_REC.ipdu_stat(i).group5(7) = TEMP_BYTE_EQV

          TEMP = 0
          CALL MVBITS(IPDU_RELAY_I2(i+1), 1, 1, TEMP, 1)
          CALL MVBITS(IPDU_RELAY_I2(i+3), 1, 1, TEMP, 0)
          IDX_REC.ipdu_stat(i).group5(8) = TEMP_BYTE_EQV

          TEMP = 0
          CALL MVBITS(IPDU_RELAY_I2(i+1), 0, 1, TEMP, 1)
          CALL MVBITS(IPDU_RELAY_I2(i+3), 0, 1, TEMP, 0)
          IDX_REC.ipdu_stat(i).group5(9) = TEMP_BYTE_EQV

          IDX_REC.ipdu_stat(i).group5(10) = IBITS(IPDU_RELAY_I2(i+5), 3, 1)
          IDX_REC.ipdu_stat(i).group5(11) = IBITS(IPDU_RELAY_I2(i+5), 2, 1)
          IDX_REC.ipdu_stat(i).group5(12) = IBITS(IPDU_RELAY_I2(i+5), 1, 1)
          IDX_REC.ipdu_stat(i).group5(13) = IBITS(IPDU_RELAY_I2(i+5), 0, 1)
          IDX_REC.ipdu_stat(i).group5(14) = IBITS(IPDU_RELAY_I2(i+5), 4, 1)
          IDX_REC.ipdu_stat(i).group5(15) = IBITS(IPDU_RELAY_I2(i+5), 5, 1)

        end do

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_GET_HSK_FLDS = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_GET_HSK_FLDS = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END

