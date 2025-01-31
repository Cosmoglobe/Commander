	INTEGER*4 FUNCTION FUT_MJF_CHANGE ( HSK_RECS, LIMFLAGS, 
	1	FIRST_MJ, CHANNEL, QUALFLAGS )
C/
C/	PROGRAM NAME:
C/	  FUT_MJF_CHANGE
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine takes as input the 2-record Housekeeping buffer,
C/	  the pointer indicating which of the 4 major frames is first in the 
C/	  buffer and the current channel number. It compares indicated status 
C/	  bits or words from the two bracketing housekeeping records and sets a
C/	  bit in one of the four quality flags if the status has changed across
C/	  the major frame boundary.
C/
C/	AUTHOR:
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 22, 1988 
C/
C/	MODIFIED BY:
C/	  For Build 3.2 only 3 HKP status bits or words are being checked.
C/	  Shirley M. Read, STX, May 1988.
C/	    Added Limflags to calling sequence so that the bit mask for
C/	    FLG_STCHG_MJ_ST will be avaliable for the AND of the limit check.
C/
C/	CALLING SEQUENCE:
C/
C/	  Status = FUT_MJF_CHANGE  (HSK_RECS, LIMFLAGS, FIRST_MJ, 
C/		   CHANNEL, QUALFLAGS )
C/
C/	INPUT PARAMETERS:
C/	  HSK_RECS(1024)	BYTE	Buffer holding 2 Housekeeping records
C?	  LIMFLAGS              Record  Limit check enable/disable flags
C/	  FIRST_MJ		I*2	Pointer to which of the 4 major frames
C/					in HSK_RECS.
C/						1 = 1st maj. fr., 1st record
C/						2 = 2nd maj. fr., 1st record
C/						3 = 1st maj. fr., 2nd record
C/						4 = 2nd maj. fr., 2nd record
C/	  CHANNEL               I*2     Current channel being processed.
C/	
C/	OUTPUT PARAMETERS:
C/	  QUALFLAGS(4)          Byte    Four byte flag corresponding to changes
C/					across the major frame boundary.
C/	
C/	INPUT/OUTPUT FILES:
C/	  NONE
C/
C/	INCLUDE FILES USED:
C/
C/	SUBROUTINES CALLED:
C/
C/	ERROR HANDLING:
C/	  Calls to Lib$Signal and its interface to Fut_Error, the condition
C/	  handler.
C/
C/	METHOD USED:
C/	  The following is the PDL --
C/        set the return status to success.
C/	  Determine which housekeeping records contain the two bracketing
C/	  major frames of data and which major frames within each record
C/	  are involved.
C/	  For each of the status words or bits to be checked, compare the
C/	  status value for the two major frames. 
C/	  If the status has changed, set the corresponding bit in one of the
C/	  four Qualflags.
C/	  AND the result with the Limflags.
C/ 	  Set the function value to the processing status.
C/	  Return;
C/	  End.
C/

	IMPLICIT	NONE

	include		'($ssdef)'

	EXTERNAL   	FUT_NORMAL
	EXTERNAL	FUT_ABERR

!	Input/Output Parameters 

        dictionary 'nfs_hkp'
        record /nfs_hkp/    hsk_recs(2)

        dictionary 'fex_limflags'
        record /fex_limflags/ limflags

	integer*2	FIRST_MJ

	integer*2	CHANNEL

        byte		QUALFLAGS(4)

!	Local variables

	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status

        integer*4       i,j             ! counters
        integer*4       mode            ! counter to loop over
                                        ! side-current configurations
        integer*4       act_mj_a(2)     ! actual a side major frame
        integer*4       act_mj_b(2)     ! actual b side major frame

	byte		dwell_1_bytes(2) !two bytes storage for 1 dwell status 
	integer*2	dwell_1_i2	  !I*2 word storage for 1 dwell status
	equivalence     (dwell_1_bytes(1),dwell_1_i2)
	byte		dwell_2_bytes(2) !two bytes storage for 2 dwell status 
	integer*2	dwell_2_i2	  !I*2 word storage for 2 dwell status
	equivalence     (dwell_2_bytes(1),dwell_2_i2)
	integer*2       dwell_1_mode, dwell_2_mode  ! Dwell mode for mjf 1 and 2

	INTEGER*2	bn(2),mj(2),CUR,
	1		ENG_A(2), ENG_B(2),
	1		FIRAS_MJ_A(2),
	1		FIRAS_MJ_B(2),
	1		I1, I2, MJFR, 
	1		sys_stat_i2, STAT_WD,
	1		zeroes(4) / 4*0 /,
	1               hsk_a_dwell,
	1               hsk_b_dwell,
	1		curside

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set return status to success.

	RETSTAT = SUCCESS
!
!     Get the Bracketing Housekeeping major frames
!     from HSK_RECS and pointer FIRST_MJ
!
        if(first_mj .eq. 1)then       ! First major frame in buffer 1 of HSK_RECS
          bn(1) = 1
          bn(2) = 1
          mj(1) = 1
          mj(2) = 2
        elseif(first_mj .eq. 2)then   ! First major frame in buffer 2 of HSK_RECS
          bn(1) = 1
          bn(2) = 2
          mj(1) = 2
          mj(2) = 1
        elseif(first_mj .eq. 3)then   ! First major frame in buffer 3 of HSK_RECS
          bn(1) = 2
          bn(2) = 2
          mj(1) = 1
          mj(2) = 2
        else                          ! First major frame in buffer 4 of HSK_RECS
          bn(1) = 2
          bn(2) = 1
          mj(1) = 2
          mj(2) = 1
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !
!     The following code may or may not be used in the final FDQ.   !
!     Check which of the 2 buffers is FIRAS Major Frame 1 and	    !
!     which is FIRAS Major Frame 2.  It is possible that within	    !
!     the same S/C major frame, Side A could be at FIRAS Major      !
!     Frame 1 while Side B is at FIRAS Major Frame 2!		    !
!								    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       stat_wd = hsk_recs(bn(1)).frame(mj(1)).hskp_head.stat_monitor_cmd(1)
!	IF (BTEST(STAT_WD, 14)) THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									  !
!     The (FIRAS) Major Frame ID is in the first command status word,	  !
!     bit 14, but since HSK_RECS is defined as a byte array, we have	  !
!     to move it to an I*2 variable, STAT_WD, in order to use the BTEST   !
!     intrinsic function.  If this bit = 0, then it is FIRAS Major 	  !
!     Frame 1, if it is 1, then it is FIRAS Major Frame 2.		  !
!									  !
!     Control is here if the BTEST is true, which means the bit is 1,     !
!     equivalent to having FIRAS Major Frame 2 here for side A.		  !
!									  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	  FIRAS_MJ_A(1) = bn(2) 	! The 2nd (S/C) major frame is
!	  FIRAS_MJ_A(2) = bn(1)  	! actually FIRAS major frame 1
!         act_mj_a(1) = mj(2)             ! adjust associated frame
!         act_mj_a(2) = mj(1)
!
!     FIRAS_MJ_A(1) is the buff num in the HSK_RECS buffer for FIRAS major frame 1
!     FIRAS_MJ_A(2) is the buff num in the HSK_RECS buffer for FIRAS major frame 2
!
!	  ENG_A(1) = 2		! ENG_A(1) indicates which of ENG_BUFFS is FIRAS major frame 1
!	  ENG_A(2) = 1		! ENG_A(2) indicates which of ENG_BUFFS is FIRAS major frame 2

!	ELSE

!	  FIRAS_MJ_A(1) = bn(1) 	! The first S/C major frame is
!	  FIRAS_MJ_A(2) = bn(2) 	! also FIRAS major frame 1
!         act_mj_a(1) = mj(1)             ! adjust associated frame
!         act_mj_a(2) = mj(2)
!	  ENG_A(1) = 1
!	  ENG_A(2) = 2

!	ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!					   !
!     Now do the same check for Side B     !
!					   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       stat_wd = hsk_recs(bn(1)).frame(mj(1)).hskp_head.stat_monitor_cmd(5)
!	IF (BTEST(STAT_WD, 14)) THEN

!	  FIRAS_MJ_B(1) = bn(2) 	! First S/C major frame
!	  FIRAS_MJ_B(2) = bn(1) 	! is FIRAS major frame 2
!         act_mj_b(1) = mj(2)
!         act_mj_b(2) = mj(1)
!	  ENG_B(1) = 2
!	  ENG_B(2) = 1

!	ELSE

!	  FIRAS_MJ_B(1) = bn(1) 	! First S/C major frame
!	  FIRAS_MJ_B(2) = bn(2) 	! is FIRAS major frame 1
!         act_mj_b(1) = mj(1)
!         act_mj_b(2) = mj(2)
!	  ENG_B(1) = 1
!	  ENG_B(2) = 2

!	ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     									  !
!     Now the status words; start with the command status words, they     !
!     also depend on whether it is FIRAS major frame 1 or 2		  !
!     									  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         do i = 1,2
!         do j = 1,4

! Move A-side status words, 8 bytes from FIRAS MJ 1 and 8 from FIRAS MJ 2

!            eng_rec.en_stat.group1(j + 4*(i-1)) = 
!    &      hsk_recs(firas_mj_a(i)).frame(act_mj_a(i))
!    &                                      .hskp_head.stat_monitor_cmd(j)

! Move B-side status words, 8 bytes from FIRAS MJ 1 and 8 from FIRAS MJ 2

!            eng_rec.en_stat.group1(j + 4*(i+1)) = 
!    &      hsk_recs(firas_mj_b(i)).frame(act_mj_b(i))
!    &                                      .hskp_head.stat_monitor_cmd(j+4)

!          end do
!        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								            !
!     Check status words:					            !
!     Dwell Mode Status, Micro Bus Readouts and Bolometer Bias statuses.    !
!				  				            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	Determine Current A or B side for channel IFG.

	If (channel .eq. 1 .or. channel .eq. 2 ) then
	  curside = 1		! A side
	Else
	  curside = 2		! B side
	Endif 

!	hsk_a_dwell = HSK_RECS(bn(1)).frame(mj(1)).hskp_head.dwell_stat(1)
                                                            ! Dwell on a side
!	hsk_b_dwell = HSK_RECS(bn(1)).frame(mj(1)).hskp_head.dwell_stat(2)
							    ! Dwell on b side
!	Temporarily for Build 3.2 the dwell status will be checked only for
!	the A or B side for a channel IFG. The above algorithm assuming A
!       corresponds to the right channels and B to the left channels will be
!	used. It is assumed that in the HKP record dwell_stat(1) = A side and
!	dwell_stat(2) = B side.

	dwell_1_bytes(1) = 
	1  HSK_RECS(bn(1)).frame(mj(1)).hskp_head.dwell_stat(curside)

	dwell_2_bytes(1) = 
	1  HSK_RECS(bn(2)).frame(mj(2)).hskp_head.dwell_stat(curside)

	call mvbits(dwell_1_i2,7,1,dwell_1_mode,0)
	call mvbits(dwell_2_i2,7,1,dwell_2_mode,0)

	IF ( dwell_1_mode .NE. dwell_2_mode ) 
     &     Qualflags(1) = Qualflags(1) .OR. '01'x

	IF(HSK_RECS(bn(1)).frame(mj(1)).hskp_head.u_proc_stat(channel).NE.
     &     HSK_RECS(bn(2)).frame(mj(2)).hskp_head.u_proc_stat(channel))
     &     Qualflags(1) = Qualflags(1) .OR. '02'x 

	IF(HSK_RECS(bn(1)).frame(mj(1)).hskp_head.bol_cmd_bias(channel).NE.
     &     HSK_RECS(bn(2)).frame(mj(2)).hskp_head.bol_cmd_bias(channel))
     &     Qualflags(1) = Qualflags(1) .OR. '04'x

!	AND the results for the four Qualflags with the Limflags.

	DO i = 1, 4
	
	   Qualflags(i) = Qualflags(i) .AND. 
	1	Limflags.Lim_Flags.FLG_STCHG_MJ_ST(i)

	ENDDO
!	Set function to return status

	if (retstat.eq.success) then
	  FUT_MJF_CHANGE = %loc(FUT_NORMAL)
	else
	  FUT_MJF_CHANGE = %loc(FUT_ABERR)
	endif

	RETURN
	END
