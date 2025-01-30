	INTEGER*4 FUNCTION FDQ_NEW_IDX ( new_segments, IDX_REC, LAST_IDX_REC,
	1	eng_rec, idxtols, idxflags,telm_flag)
C/
C/	PROGRAM NAME:
C/	  FDQ_NEW_IDX
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine does the comparision of the "current" Index record buffer
C/	  with the last index record written and return "true" if a new Index
C/	  record ought to be written, "false" otherwise.
C/
C/	AUTHOR:
C/	  J. Durachta
C/	  STX
C/	  March 11, 1987
C/        Based on the work of Edwin Fung
C/
C/	MODIFIED BY:
C/	  E. FUNG
C/	  GSFC
C/	  JUNE 7, 1987
C/	  REASON:	To add a new parameter NEW_SEGMENTS.  This is passed
C/			from the calling routine.  If this is true it signals
C/			to FDQ_NEW_IDX that it is the start of a  new
C/			set of segments so it should start with a new index
C/			record.
C/
C/        J. Durachta
C/        ARC
C/        July 28, 1987
C/        REASON:       Altered to prevent index creation for analog values
C/                      which toggle between 2 values due to digitization.
C/
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 6, 1988
C/	  REASON: 	Changed to an I*4 function, added error checking and 
C/			inserted calls to Lib$Signal for interface with 
C/			the Fut_Error condition handler. 
C/
C/	   Shirley M. Read
C/	   STX
C/	   June 3, 1988
C/	   REASON:      Several FDQ text files were moved to the FUT Facility.
C/			Changed include file names to FUT. This routine includes
C/			FUT_Text_Item.Txt. Added idxtols and idxflags to the
C/			input calling sequence. These records are now obtained
C/			from CCT_Get_Config.
C/
C/	   H. Wang
C/	   STX
C/	   Jan. 29, 1991
C/	   REASON:      New requirements for FDQ
C/                New quantities will be tested for the generation of
C/                Idex records
C/                * Bolometer CMD Bias
C/                * IPDU relays
C/                * LVDT status
C/                * MTM CAL MOTOR
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
CH	Version XXX, 30 August 1990, SPR 4171, Larry P. Rosen, STX
CH		Standardize the IDX report file name to by of form:
CH		FDQ_IDX_yydddhh_yydddhh.REP_yydddhhmm  , where the times are the
CH		start, stop, and run.  REPORT logical flag whether to make file.
C-----------------------------------------------------------------------------
C/
C/	CALLING SEQUENCE:
C/	  Return status = FDQ_NEW_IDX (new_segments, IDX_REC, LAST_IDX_REC,
C/			eng_rec, idxtols, idxflags,telm_flag)
C/
C/	INPUT PARAMETERS:
C/	  new_segments		L*1	Flag to indicate whether we are
C/					starting with a new set of segments
C/	  IDX_REC(512)		BYTE	Current Index record buffer
C/	  LAST_IDX_REC(512)	BYTE	Buffer containing the last written
C/					Index record
C/	  eng_rec	        Rec	Current ENG record
C/	  IDXTOLS               Rec     Current index tolerances.
C/ 	  IDXFLAGS              Rec     Current index check enable flags.
C/        TELM_FLAG             I*2     Telemetry quality of Eng. record
C/                                      0 = good
C/                                      1 = bad
C/	OUTPUT PARAMETERS:
C/	  None other than the function value of FDQ_NEW_IDX itself.
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
C/	  TBD
C/	METHOD USED:
C/	  PDL for this routine --
C/
C/	  Make the IDX report file name.
C/	  If (first time through) then
C/	    Set first time flag to false;
C/	    Set NEW_SEGMENTS to false;
C/	    Move current index buffer to last written index buffer;
C/	    Set FDQ_NEW_IDX to false;
C/	    Return;
C/	  Endif;
C/
C/	  If NEW_SEGMENTS (but not first time through) then
C/	    Set NEW_SEGMENTS to false;
C/	    Move current index buffer to last written index buffer;
C/	    Set FDQ_NEW_IDX to false;
C/	    Return;
C/	  Endif;
C/	  Compare channel-specific fields;
C/	  Compare engineering status fields;
C/	  Compare analog fields;
C/
C/	  If #fields that have changed > 0 then
C/	    Move tolerances to Index buffer;
C/	    Set FDQ_NEW_IDX to true;
C/	  Else
C/	    Set FDQ_NEW_IDX to false;
C/	  Endif;
C/
C/	  Return;
C/	  End.

	IMPLICIT	NONE


        include         '(fdq_new_idx_pars)'

        dictionary 'fdq_idx'
        record /fdq_idx/      idx_rec
        record /fdq_idx/      last_idx_rec

        dictionary 'fdq_eng'
        record/fdq_eng/  eng_rec
        dictionary 'fex_idx_tols'
        record/fex_idx_tols/  idxtols
        dictionary 'fex_idx_flag'
        record/fex_idx_flag/  idxflags

	EXTERNAL	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_WRITEIDX
	EXTERNAL	FDQ_LUNGETERR
	EXTERNAL	FDQ_LUNFREER
	EXTERNAL	FDQ_OPENERR
	EXTERNAL	FDQ_WRITERR
	EXTERNAL	FDQ_READERR
	EXTERNAL	FDQ_CLOSERR
	EXTERNAL	FDQ_IDXFLGTRUE
	EXTERNAL        FDQ_IDXTOLONE
	EXTERNAL	FDQ_NOIDXRPT
	INTEGER*4	Lib$Get_Lun
	INTEGER*4	Lib$Free_Lun

	INCLUDE 	'($SSDEF)'

	LOGICAL*1	FIRST /.TRUE./, new_segments,
	1		tols_byte(90), idx_flags(284),
	1		new_idx			! Now a local flag

	INTEGER*2	class, cur, I,
	1		ITEMS_CHANGED, J, K,
	1		place,
	1		side, TC(8),
	1		temp1, temp2,
	1		TOLS(256), TELM_FLAG
	INTEGER*4      	NEW_I4, OLD_I4

	REAL*4		TOLS_R4(90)
	integer*4	STATUS		! Status variable
	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status
	integer*4	ZERO / 0 /      ! Good return status for Fortran I/O
	integer*2       FV / -9999 /    ! Flag value for missing or bad GRT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set return status to success.

	RETSTAT = SUCCESS

	IF (FIRST .and. (telm_flag .eq. 0)) THEN
	    FIRST = .FALSE.
	    new_segments = .false.
	    DO i = 1, nclass		! tols is in percent units
	      tols_r4(i) = float(idxtols.idx_tols.tols(i)) / 100.
	    End Do				! tols_r4 is a fraction

	    idx_rec.header.gmt_start_time_ascii = eng_rec.ct_head.gmt
	    idx_rec.header.binary_start_time(1) = eng_rec.ct_head.time(1)
	    idx_rec.header.binary_start_time(2) = eng_rec.ct_head.time(2)

	    do i = 1, nclass
	      idx_rec.idx_tail.tols(i) = idxtols.idx_tols.tols(i)
	    enddo
	    call lib$movc3 (512, idx_rec, last_idx_rec)
	    new_idx = .false.

	    If ((retstat .eq. success) .and. (new_idx)) Then
		FDQ_NEW_IDX = %loc(FDQ_WRITEIDX)
	    elseif ( retstat .eq. success ) then
		FDQ_NEW_IDX = %loc(FDQ_NORMAL)
 	    else
		FDQ_NEW_IDX = %loc(FDQ_ABERR)
	    endif

	    return
	endif 	! IF (FIRST)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									   !
!     The following IF-block is added June 7, 1987 for DQ version 2.0,	   !
!     to handle multiple set of SCI segments (and creates separate RSI	   !
!     and ENG segments instead of one big segment).  Because of this	   !
!     new requirement a new flag NEW_SEGMENTS has to be passed down to     !
!     indicate to FDQ_NEW_IDX that it should react differently	   !
!     (almost like the first time).  The following IF-block is the	   !
!     code when we just start a new set of SCI segments.  Note that	   !
!     it does "RETURN" at the end of the IF-block.			   !
!									   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (new_segments .and.(TELM_FLAG .eq. 0)) then
	    new_segments = .false.
!
!

	    IDX_REC.HEADER.GMT_START_TIME_ASCII = ENG_REC.CT_HEAD.GMT
	    IDX_REC.HEADER.BINARY_START_TIME(1) = ENG_REC.CT_HEAD.TIME(1)
	    IDX_REC.HEADER.BINARY_START_TIME(2) = ENG_REC.CT_HEAD.TIME(2)
	     call lib$movc3 (512, idx_rec, last_idx_rec)
	     new_idx = .false.
	     if ((retstat .eq. success) .and. (new_idx)) then
		FDQ_NEW_IDX = %loc(FDQ_WRITEIDX)
	    elseif ( retstat .eq. success ) then
		FDQ_NEW_IDX = %loc(FDQ_NORMAL)
	    else
		FDQ_NEW_IDX = %loc(FDQ_ABERR)
	    endif

	    return
	Endif		! if (new_segments)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!     Begin block for determining whether new index record should be   !
!     written on calls subsequent to start of new segments.            !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
	Do i = 1, 284
	  idx_flags(i) = idxflags.idx_flags.flags(i)
	Enddo
	Do i = 1, nclass
	  tols_r4(i) = float(idxtols.idx_tols.tols(i)) / 100.	! tols is in percent units
	Enddo					! tols_r4 is a fraction
!
!     External Calibrator Positions
!
	ITEMS_CHANGED = 0
        if(idx_flags(1))then
          if (last_idx_rec.XCAL.POS(1) .ne. idx_rec.XCAL.POS(1))then
	    ITEMS_CHANGED = ITEMS_CHANGED + 1
	  endif
        endif

        if(idx_flags(2))then
	  IF (LAST_IDX_REC.xcal.pos(2) .NE. IDX_REC.xcal.pos(2))THEN
	    ITEMS_CHANGED = ITEMS_CHANGED + 1
	  ENDIF
        endif

!
!     Compare channel-specific fields
!

	DO k = 1,4                         ! k represents the channels
	  do j = 1,6
            if(idx_flags(2 + j + (k-1)*9))then
	      IF (LAST_IDX_REC.chan(k).group6(j) .NE.
     &               IDX_REC.chan(k).group6(j))then
	        ITEMS_CHANGED = ITEMS_CHANGED + 1
	      ENDIF
            endif
	  enddo		! j = 1,6

          if(idx_flags(9 + (k-1)*9))then
	    IF (idx_rec.chan(k).science_gain .NE.
     &                        last_idx_rec.chan(k).science_gain) THEN
	      ITEMS_CHANGED = ITEMS_CHANGED + 1
	    ENDIF
          endif

!         Interval check added 7/24/87

          if(idx_flags(10 + (k-1)*9))then
	    IF (IDX_REC.chan(k).bol_cmd_bias .NE.
     &               LAST_IDX_REC.chan(k).bol_cmd_bias) then

	      ITEMS_CHANGED = ITEMS_CHANGED + 1
	    ENDIF
          endif
!
!     Bolometer Voltage Readouts (Analog quantities)
!
          if(idx_flags(11 + (k-1)*9))then
	    temp2 = ZEXT(last_idx_rec.chan(k).bol_readout_volt)
	    temp1 = ZEXT(idx_rec.chan(k).bol_readout_volt)
	    temp1 = temp1 - temp2
            if(iiabs(temp1) .gt. 1)then
	      IF (TEMP2 .EQ. 0) THEN
	        items_changed = items_changed + 1
	      ELSE
	        if((float(abs(temp1)) / float(temp2)) .ge. 
     &              tols_r4(bol_rdout_tol) ) then
	          items_changed = items_changed + 1
	        endif
	      ENDIF
            endif     !if(iiabs
          endif       !if(idx_
    	ENDDO		! Do K (each of the 4 channels)

!
!     Compare GRT's
!

        DO side = 0,1                   ! side = a, b
          do cur = 1,2                  ! current = low, high
            do j = 1,16                 ! loop over the grts

              if(j .le. 4)then
                class = grt_controlled_lo_tol
              elseif(j .eq. 5)then
                class = grt_dihedral_lo_tol
              elseif(j .le. 9)then
                class = grt_bol_lo_tol
              elseif(j .eq. 10)then
                class = grt_mir_lo_tol
              elseif(j .le. 14)then
                class = grt_controlled_lo_tol
              else
                class = grt_mir_lo_tol
              endif

              if(idx_flags(38 + j + (cur-1)*16 + side*32))then
                old_i4 = last_idx_rec.grts(side*2 + cur).group1(j) 
                new_i4 = idx_rec.grts(side*2 + cur).group1(j) 
                if ((( new_i4 .eq. FV ) .and. ( old_i4 .ne. FV )) .OR.
	1	    (( old_i4 .eq. FV ) .and. ( new_i4 .ne. FV ))) then
     	          new_i4 = FV			! Check for flag 10/17/88
     	        else
                  new_i4 = new_i4 - old_i4
     		endif
                if((jiabs(new_i4) .gt. 1) .or. (new_i4 .eq. FV))then  ! 10/17/88
                  if(old_i4 .eq. 0)then

                    items_changed = items_changed + 1
                  elseif((float (abs(new_i4)) / (float(old_i4)))
     &                    .ge. tols_r4(class + cur-1))then
                    items_changed = items_changed + 1
                  endif    !if(old_i4
                endif      !if(jiabs
              endif        !if(idx_

            end do         ! j = 1,16
          end do           ! cur = 0,1
        END DO             ! side = 0,1

!
!     Temperature Controllers
!
	DO j = 1,8
          if(idx_flags(102 + j))then
            OLD_I4 =  last_idx_rec.temp_ctrl.group2(j)
            NEW_I4 =  idx_rec.temp_ctrl.group2(j)
            new_i4 = new_i4 - old_i4
            if(jiabs(new_i4) .gt. 1)then
              if (old_i4 .eq. 0) then
                ITEMS_CHANGED = ITEMS_CHANGED + 1
	      elseif((FLOAT(ABS(NEW_I4)) / FLOAT(OLD_I4)) 
     &                                   .GE. TOLS_R4(TC_TOL)) THEN
	        ITEMS_CHANGED = ITEMS_CHANGED + 1
	      endif    !if(old_i4
	    endif      !if(jiabs
          endif        !if(idx_

	END DO
!
!     Status Monitor Bits (A & B sides)
!
	DO k = 0,1	! Side A to Side B
	  do j = 0,28	! 29 items on each side
            if(idx_flags(111 + j + k*29))then
              IF (IDX_REC.stat_mon(k+1).group3(j+1) .NE.
     &                  LAST_IDX_REC.stat_mon(k+1).group3(j+1))then
	        ITEMS_CHANGED = ITEMS_CHANGED + 1
	      ENDIF
            endif       !if(idx_
	  end do			! DO J (29 items on each side)

!     Status Monitor Sync (Major Frame ID Toggle)

          if(idx_flags(169 + k))then
	    if (idx_rec.idx_sync.side(k+1) .ne.
     &                  last_idx_rec.idx_sync.side(k+1)) then
	      items_changed = items_changed + 1
	    endif
          endif            !if(idx_

	END DO			! DO K (Side A to Side B)

!
!     IPDU Individual Statuses Plus Dwell (Side A & B)
!

	DO	k = 0,1 	! IPDU: Side A to Side B

	  do	j = 0,14	! 15 items on each side

            if(idx_flags(171 + j + k*15))then
	      IF (IDX_REC.ipdu_stat(k+1).group5(j+1) .NE.
     &                  LAST_IDX_REC.ipdu_stat(k+1).group5(j+1))then
	        ITEMS_CHANGED = ITEMS_CHANGED + 1
	      ENDIF
            endif        ! if(idx_

	  end do			! DO J (15 IPDU items on each side)

          if(idx_flags(201 + k))then
	    if(idx_rec.misc_stat.dwell(2*k+1) .ne. 
     &                  last_idx_rec.misc_stat.dwell(2*k+1)) then
	      items_changed = items_changed + 1
	    endif
          endif         ! if(idx_

          if(idx_flags(203 + k))then
	    if(idx_rec.misc_stat.dwell(2*(k+1)) .ne. 
     &             last_idx_rec.misc_stat.dwell(2*(k+1)))then
	      items_changed = items_changed + 1
	    endif
          endif         ! if(idx_

	END DO			! DO K (IPDU: SIDE A TO SIDE B)
!
!     Microprocessor Bus Readouts
!
	do j = 0,3	! 4 channels

          if(idx_flags(205 + j))then
	    IF(IDX_REC.misc_stat.stat_rd_bus(j+1) .NE. 
     &             LAST_IDX_REC.misc_stat.stat_rd_bus(j+1))THEN
	      ITEMS_CHANGED = ITEMS_CHANGED + 1
	    ENDIF
          endif         ! if(idx_

	end do			! DO J (4 channels of Microprocessor readouts)
!
!     SIS Statuses
!
	do j = 0,3 	! 4 Power A statuses

          if(idx_flags(209 + j))then
	    IF (IDX_REC.misc_stat.power_a_status(j+1) .NE. 
     &             LAST_IDX_REC.misc_stat.power_a_status(j+1))then
	      ITEMS_CHANGED = ITEMS_CHANGED + 1
	    ENDIF
          endif         ! if(idx_

	end do			! do j(4 Power A statuses

	do j = 0,3 	! 4 Power B statuses

          if(idx_flags(213 + j))then
	    IF (IDX_REC.misc_stat.power_b_status(j+1) .NE. 
     &             LAST_IDX_REC.misc_stat.power_b_status(j+1))then
	      ITEMS_CHANGED = ITEMS_CHANGED + 1
	    ENDIF
          endif         ! if(idx_

	end do			! do j(4 Power B statuses

!
!     Byte-Size Analogs
!

	do	j = 0,13

          if(idx_flags(217 + j))then
	    IF(j .EQ. 10 .OR. j .EQ. 11)THEN	! The 2 hot spot heater currents
	      if(idx_rec.idx_analog.temp(j+1) .lt.
     &            last_idx_rec.idx_analog.temp(j+1)-1 .or.
     &            idx_rec.idx_analog.temp(j+1) .gt.
     &            last_idx_rec.idx_analog.temp(j+1)+1)then
	        items_changed = items_changed + 1
	      endif
	    ELSE				! "Other temperatures"
	      TEMP2 = ZEXT(LAST_IDX_REC.idx_analog.temp(j+1))
	      TEMP1 = ZEXT(IDX_REC.idx_analog.temp(j+1))
	      TEMP1 = TEMP1 - TEMP2
              if(iiabs(temp1) .gt. 1)then
                if(temp2 .eq. 0) then
	          ITEMS_CHANGED = ITEMS_CHANGED + 1
                elseif( (FLOAT(ABS(TEMP1)) / FLOAT(TEMP2)) .GE.
     &                 TOLS_R4(OTHER_TEMP_TOL) ) THEN
	          ITEMS_CHANGED = ITEMS_CHANGED + 1
	        ENDIF       !if(temp2
	      endif         !if(iiabs
	    ENDIF           !if(j
          endif             !if(idx_

	end do
!
!     IPDU voltages (20 of them) and IPDU currents (12 of them)
!
	do j = 0, 31

          if(idx_flags(231 + j))then
	    IF (J .LE. 5) THEN			! IPDU VOL 1 - 6
	      CLASS = IPDU_DIG_VOL_TOL
	    ELSEIF (J .LE. 9) THEN		! IPDU VOL 7 - 10
	      CLASS = IPDU_ANA_VOL_TOL
	    ELSEIF (J .LE. 11) THEN		! IPDU VOL 11 - 12
	      CLASS = IPDU_PREREG_VOL_TOL
	    ELSEIF (J .LE. 13) THEN		! IPDU VOL 13 - 14
	      CLASS = IPDU_28V_VOL_TOL
	    ELSEIF (J .LE. 17) THEN		! IPDU VOL 15 - 18
	      CLASS = IPDU_15V_VOL_TOL
	    ELSEIF (J .LE. 19) THEN		! IPDU VOL 19 - 20
	      CLASS = IPDU_5V_VOL_TOL
	    ELSEIF (J .LE. 21) THEN		! IPDU CUR 1 - 2
	      CLASS = IPDU_PREREG_CUR_TOL
	    ELSEIF (J .LE. 23) THEN		! IPDU CUR 3 - 4
	      CLASS = IPDU_ANA_CUR_TOL
	    ELSEIF (J .LE. 25) THEN		! IPDU CUR 5 - 6
	      CLASS = IPDU_DIG_CUR_TOL
	    ELSEIF (J .LE. 29) THEN		! IPDU CUR 7 - 10
	      CLASS = IPDU_CONST_CUR_TOL
	    ELSE				! IPDU CUR 11 - 12
	      CLASS = IPDU_INT_CUR_TOL
	    ENDIF

	    temp2 = ZEXT(last_idx_rec.idx_analog.v_and_i(j+1))
	    temp1 = ZEXT(idx_rec.idx_analog.v_and_i(j+1))
	    temp1 = temp1 - temp2
            if(iiabs(temp1) .gt. 1)then
              IF(TEMP2 .EQ. 0)THEN
	        items_changed = items_changed + 1
	      elseif((float(abs(temp1))/float(temp2)) .ge. tols_r4(CLASS))then
	        items_changed = items_changed + 1
	      endif  !if(temp2
	    ENDIF    !if(iiabs
          endif      !if(idx_

	enddo

!
!     MTM CAL MOTOR
!
        if(idx_flags(263))then
	  if (idx_rec.idx_analog.mtm_cal_motor(1) .NE. 
     &         last_idx_rec.idx_analog.mtm_cal_motor(1)) then
	    items_changed = items_changed + 1
	  endif
        endif

        if(idx_flags(264))then
	  if (idx_rec.idx_analog.mtm_cal_motor(2) .NE. 
     &         last_idx_rec.idx_analog.mtm_cal_motor(2))then
	    items_changed = items_changed + 1
	  endif
        endif

!
!     Notch Filter Statuses
!
	do j = 1,10        ! 10 Statuses, 5A side & 5 B side

          if(idx_flags(264 + j))then
	    IF (IDX_REC.misc_stat.ntch_fltr_stat(j) .NE. 
     &             LAST_IDX_REC.misc_stat.ntch_fltr_stat(j))then
	      ITEMS_CHANGED = ITEMS_CHANGED + 1
	    ENDIF
          endif         ! if(idx_
        end do

!
!     IPDU_relay
!
	do j = 1,8        ! 8 IPDU relay

          if(idx_flags(274 + j))then
	    IF (IDX_REC.IPDU_relay.group4(j) .NE. 
     &             LAST_IDX_REC.IPDU_relay.group4(j)) then
	      ITEMS_CHANGED = ITEMS_CHANGED + 1
	    ENDIF
          endif         ! if(idx_
        end do
!
!     LVDT Statuses
!
	do j = 1,2        ! 2 LVDT STATUS

          if(idx_flags(282 + j))then
	    IF (IDX_REC.misc_stat.LVDT_STAT(j) .NE. 
     &             LAST_IDX_REC.misc_stat.LVDT_STAT(j))then
	      ITEMS_CHANGED = ITEMS_CHANGED + 1
	    ENDIF
          endif         ! if(idx_
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							      !
!     Now see if any changes occur and set function value     !
!							      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if ((items_changed .gt. 0) .and. (telm_flag .eq. 0)) then
!
!     Move tolerances to index buffer
!
	  do i = 1, nclass
	    idx_rec.idx_tail.tols(i) = idxtols.idx_tols.tols(i)
	  enddo
	  new_idx = .true.
	else
	  new_idx = .false.
	endif

!	Set function to return status

	if ((retstat .eq. success) .and. (new_idx)) then
	  FDQ_NEW_IDX = %loc(FDQ_WRITEIDX)
	elseif ( retstat .eq. success ) then
	  FDQ_NEW_IDX = %loc(FDQ_NORMAL)
	else
	  FDQ_NEW_IDX = %loc(FDQ_ABERR)
	endif

	RETURN
	END

