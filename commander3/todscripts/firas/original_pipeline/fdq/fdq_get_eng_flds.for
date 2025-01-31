	INTEGER*4 FUNCTION FDQ_GET_ENG_FLDS (ENG_REC, IDX_REC)
C/
C/	PROGRAM NAME:
C/	  FDQ_GET_ENG_FLDS 
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine picks up the pertinent fields in the Engineering
C/	  record buffer which are related to the Index record.
C/
C/	AUTHOR:
C/	  Edwin Fung
C/	  GSFC
C/	  March 30, 1986
C/
C/	MODIFIED:
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 14, 1988
C/	  REASON: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal.
C/
C/	  Shirley M. Read
C/	  STX
C/	  June 3, 1988
C/	  REASON: 	Due to changes in telemetry and thus to NFS_HKP and 
C/			FDQ_ENG records, the storing algorithm for Sys_Status
C/			in the index record was modified.
C/
C/	  H. Wang
C/	  STX
C/	  Jan. 29, 1991
C/	  REASON:   New requirements for FDQ
C/                1)  Get the LVDT status from the FDQ_ENG
C/                2)  Get Fakeit data from the FDQ_ENG
C/
CH	CHANGE LOG: 		New Format for Build 4.2  STX, 10/15/88
CH
CH      Version 4.1.1 10/15/88, SPR 2622, Shirley M. Read, STX
CH              Flag values for missing conversions need changes in FDQ.
CH	        The flag values for converted fields which are out of range
CH		need the same changes. GRTs in dwell mode need a flag value
CH              also. The SWG has selected a flag value of -9999 for all cases.
CH              The limit checking algorithms need to bypass setting the quality
CH              flags if the GRT or engineering analog has the flag value.
CH		Interpolation and averaging algorithms need to be modified to
CH              include the flag value. 
CH
CH      Spr 7892, Correct the Power Status data in FDQ_IDX record
CH                H. Wang, STX, Dec 5, 1990
C-----------------------------------------------------------------------------
C/
C/	CALLING SEQUENCE:
C/	  Status = FDQ_GET_ENG_FLDS  (CUR_ENG_REC, IDX_REC)
C/
C/	INPUT PARAMETERS:
C/	  CUR_ENG_REC(1024)		BYTE	Buffer holding the current Engineering record
C/
C/	OUTPUT PARAMETERS:
C/	  IDX_REC(512)	BYTE	Buffer holding the current Index record
C/
C/	INPUT/OUTPUT FILES:
C/	  NONE
C/
C/	SUBROUTINES CALLED:
C/	  TBD
C/
C/	INCLUDE FILES USED:
C/	  FIR_IDX_PARMS.INC
C/
C/	ERROR HANDLING:
C/	  TBD
C/
C/	METHOD USED:
C/	  To be filled in later
C/

	IMPLICIT	NONE

        dictionary 'fdq_eng'
        record /fdq_eng/   eng_rec

        dictionary 'fdq_idx'
        record /fdq_idx/   idx_rec

	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR

	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status

	INTEGER*2	eng_stmon(2,0:7),
	1		GRT_I2, TC_I2, 
	1		I, J, K, L, TWO / 2 /

	INTEGER*2       IDXF, ENGF  ! Temp storage for byte fields
	BYTE            IDXB(2), ENGB(2) 
	EQUIVALENCE     (IDXF, IDXB(1))
	EQUIVALENCE     (ENGF, ENGB(1))

        integer*4       cur_grt           !counter
	REAL*4		ENG_TC, grt_temp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


C	Set return status to success.

	RETSTAT = SUCCESS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									  !
!     Move GRT temperatures after converting to milli-Kelvin integers     !
!     (skipping the cal resistors!)
!									  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do i = 1,4          ! loop over side a & b and current modes
          do cur_grt = 1,16     ! loop over the grts

            grt_temp = eng_rec.en_analog.grt(cur_grt + (i-1)*16)

            if(cur_grt.ge.11 .and. cur_grt.le.14)then       ! put cal resist counts in index rec
              grt_i2 = nint(grt_temp)
            elseif(grt_temp .eq. -9999.0)then
              grt_i2 = -9999
            elseif(grt_temp .lt. 32.767)then
              grt_i2 = nint(1000.*grt_temp)
            else
              grt_i2 = nint(-100.*grt_temp)
            endif

          idx_rec.grts(i).group1(cur_grt) = grt_i2

          end do                ! cur_grt = 1,16
        end do                  ! i = 1,4

!
!    Temperature controllers
!    Since Build 4.0 there have been no conversions for the commanded 
!    readbacks which were stored in the temperature controller fields.
!    The counts are put into the engineering record in the floating point
!    fields. There will be no flag value. Thus the current algorithm is
!    modified to the following code. 10/18/88
!
	do	i = 1, 8
          eng_tc = eng_rec.en_analog.temp_ctrl(i)
	  TC_I2 = NINT (ENG_TC)
          idx_rec.temp_ctrl.group2(i) = tc_i2
	end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!					   !
!     The status monitor fields            !
!					   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do i = 1,8
          eng_stmon (1,i-1) = eng_rec.en_stat.group1(i)
          eng_stmon (2,i-1) = eng_rec.en_stat.group1(i+8)
        end do

c        DO I = 1,2
c          DO J = 1,2
c            IDX_REC.CHAN((I-1)*2 + J).FAKEIT = 
c     &                          IBITS (eng_stmon(i,0), (10-J), 1)   ! Fake it 
c          END DO
c        END DO

         DO I = 1,4
          IDX_REC.chan(i).fakeit = eng_rec.chan(i).fakeit
         enddo   
        DO I = 1,2           !LOOP OVER SIDES A & B

          IDX_REC.XCAL.POS(I) = IBITS (eng_stmon(i,0), 0, 2)   ! Bits 0-1

          idx_rec.stat_mon(i).group3(1) = IBITS (eng_stmon(i,0), 15, 1)! TC Integrator Circuit Statuses
          idx_rec.stat_mon(i).group3(2) = IBITS (eng_stmon(i,0), 13, 1)
          idx_rec.stat_mon(i).group3(3) = IBITS (eng_stmon(i,0), 12, 1)
          idx_rec.stat_mon(i).group3(4) = IBITS (eng_stmon(i,4), 13, 1)
          idx_rec.stat_mon(i).group3(5) = IBITS (eng_stmon(i,4), 12, 1)

          idx_rec.stat_mon(i).group3(6) = IBITS (eng_stmon(i,0), 11, 1)	! Dither
          idx_rec.stat_mon(i).group3(7) = IBITS (eng_stmon(i,0), 10, 1)

          idx_rec.stat_mon(i).group3(8) = IBITS (eng_stmon(i,0), 7, 1)	! ROM
          idx_rec.stat_mon(i).group3(9) = IBITS (eng_stmon(i,0), 6, 1)

          idx_rec.stat_mon(i).group3(10) = IBITS (eng_stmon(i,0), 4, 2)	! Bits 4-5
          idx_rec.stat_mon(i).group3(11) = IBITS (eng_stmon(i,0), 2, 2)	! Bits 2-3

          idx_rec.stat_mon(i).group3(12) = IBITS (eng_stmon(i,3), 14, 2)
          idx_rec.stat_mon(i).group3(13) = IBITS (eng_stmon(i,3), 12, 2)
          idx_rec.stat_mon(i).group3(14) = IBITS (eng_stmon(i,7), 14, 2)
          idx_rec.stat_mon(i).group3(15) = IBITS (eng_stmon(i,7), 12, 2)

          idx_rec.stat_mon(i).group3(16) = IBITS (eng_stmon(i,3), 9, 3)
          idx_rec.stat_mon(i).group3(17) = IBITS (eng_stmon(i,3), 3, 3)
          idx_rec.stat_mon(i).group3(18) = IBITS (eng_stmon(i,7), 9, 3)
          idx_rec.stat_mon(i).group3(19) = IBITS (eng_stmon(i,7), 3, 3)

          idx_rec.stat_mon(i).group3(20) = IBITS (eng_stmon(i,3), 6, 3)
          idx_rec.stat_mon(i).group3(21) = IBITS (eng_stmon(i,3), 0, 3)
          idx_rec.stat_mon(i).group3(22) = IBITS (eng_stmon(i,7), 6, 3)
          idx_rec.stat_mon(i).group3(23) = IBITS (eng_stmon(i,7), 0, 3)

          idx_rec.stat_mon(i).group3(24) = IBITS (eng_stmon(i,4), 15, 1)
          idx_rec.stat_mon(i).group3(25) = IBITS (eng_stmon(i,4), 5, 1)
          idx_rec.stat_mon(i).group3(26) = IBITS (eng_stmon(i,4), 4, 1)
          idx_rec.stat_mon(i).group3(27) = IBITS (eng_stmon(i,4), 3, 1)
          idx_rec.stat_mon(i).group3(28) = IBITS (eng_stmon(i,4), 2, 1)
          idx_rec.stat_mon(i).group3(29) = IBITS (eng_stmon(i,4), 1, 1)

!                            30-32 are spares

          idx_rec.stat_mon(i).group3(33) = IBITS (eng_stmon(i,4), 0, 0)

        END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!				       !
!     Dwell statuses and addresses     !
!				       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        IDX_REC.MISC_STAT.DWELL(1) = eng_rec.en_stat.dwell_stat(1)
        IDX_REC.MISC_STAT.DWELL(3) = eng_rec.en_stat.dwell_stat(2)
        IDX_REC.MISC_STAT.DWELL(2) = eng_rec.en_stat.grt_addr(1)
        IDX_REC.MISC_STAT.DWELL(4) = eng_rec.en_stat.grt_addr(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!				       !
!     LVDT statuses                    !
!				       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IDX_REC.MISC_STAT.LVDT_STAT(1) = eng_rec.en_stat.LVDT_STAT(1)
        IDX_REC.MISC_STAT.LVDT_STAT(2) = eng_rec.en_stat.LVDT_STAT(2)

!
!     micro-processor readouts
!

        DO I = 1,4
          IDX_REC.MISC_STAT.STAT_RD_BUS(I) = ENG_REC.EN_STAT.MICRO_STAT_BUS(I)
        END DO
!
!     Bolometer bias commanded status
!
        DO I = 1,4
          IDX_REC.CHAN(I).BOL_CMD_BIAS = ENG_REC.EN_STAT.BOL_CMD_BIAS(I)
        END DO
!
!     POWER statuses
!
	K = 0
	DO I = 1, 4
	  L = 0
	  DO J = 1, 2
	    ENGB(1) = ENG_REC.EN_STAT.POWER_A_STATUS(J)
	    CALL MVBITS( ENGF, K, TWO, IDXF, L)
	    IDX_REC.MISC_STAT.POWER_A_STATUS(I) = IDXB(1)
	    L = L + 2
	  ENDDO
	  K = K + 2
	ENDDO
	K = 0
	DO I = 1, 4
	  L = 0
	  DO J = 1, 2
	    ENGB(1) = ENG_REC.EN_STAT.POWER_B_STATUS(J)
	    CALL MVBITS( ENGF, K, TWO, IDXF, L)
            IDX_REC.MISC_STAT.POWER_B_STATUS(I) = IDXB(1)
	    L = L + 2
	  ENDDO
	  K = K + 2
	ENDDO
        Idx_rec.misc_stat.lvdt_stat(1)=eng_rec.en_stat.lvdt_stat(1)
        Idx_rec.misc_stat.lvdt_stat(2)=eng_rec.en_stat.lvdt_stat(2)

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_GET_ENG_FLDS = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_GET_ENG_FLDS = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END
