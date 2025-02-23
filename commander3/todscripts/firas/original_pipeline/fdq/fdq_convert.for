        INTEGER*4 FUNCTION FDQ_CONVERT (HSK_RECS, FIRST_MJ, CAL_LUN,
        1       Time_range,ENG_BUFFS)
C/
C/      PROGRAM NAME:
C/        FDQ_CONVERT
C/
C/      PROGRAM DESCRIPTION:
C/        This routine takes as input the housekeeping buffer for 2
C/        records as well as the pointer to the first HSKP major
C/        frame and converts the counts into engineering units.  It
C/        makes use of the common CONVBUFF which contains polynomial
C/        coefficients and GRT look-up tables, which have already
C/        been loaded (off from disk files) by subroutine FDQ_LOAD_CONV.
C/
C/      AUTHOR:
C/        EDWIN H. FUNG
C/        GSFC
C/        OCTOBER 24, 1985
C/      
C/      MODIFIED BY:
C/        E.FUNG
C/        GSFC
C/        NOV. 6, 1985
C/        REASON:       To add FDQ_FIRNAMES.TXT include file, so that when FDQ_OL_GRT_CVT
C/                      return an abnormal status it will be able to report to
C/                      SYS$OUTPUT and the log file which GRT has the problem.
C/
C/        E.FUNG
C/        GSFC
C/        NOV. 14, 1985
C/        REASON:       To add logic to handle the case when there is no
C/                      GRT lookup tables defined for a GRT (manifested by
C/                      a calibration group ID of 0).
C/      
C/        E.FUNG
C/        GSFC
C/        DEC. 7, 1985
C/        REASON:       To change the writing to the log file to be D-lines.
C/                      Also there is a new input parameter, to indicate
C/                      whether any logging to LUN 6 is desired at all.
C/
C/        E. FUNG
C/        GSFC
C/        FEB. 26, 1986
C/        REASON:       To use include file NARCVSIZE.INC instead of ARCVSIZE.INC
C/                      to correspond to the new archive formats for the 2nd ITD.
C/
C/        D. WARD
C/        GSFC/STX
C/        OCT. 13, 1987
C/        REASN:        To delete NARCVSIZE.INC from the include list.
C/
C/        D. WARD
C/        GSFC/STX
C/        OCT. 20, 1987
C/        REASN:        To delete HSKINDX AND ENGINDX from the include list.
C/
C/        Shirley M. Read
C/        STX
C/        January 6, 1988
C/        REASON:       Converted from subroutine to function for interface
C/                      with Fut_Error condition handler. Added error checking
C/                      and calls to Lib$Signal. Add poly coefficient index
C/                      and GRT index to call sequences for signal of database
C/                      names from FDQ_FIRNAMES.TXT when errors occur.
C/        Shirley M. Read
C/        STX
C/        May 17, 1988
C/        REASON:       To store the values of the calibrator resistors from
C/                      the housekeeping record into the engineering record.
C/                      Since there are no poly coefficients for the command
C/                      readback bits of the temperature controllers, the
C/                      code was modified to store the counts instead of the
C/                      converted values. (Cur_Poly 1 to 8 in FDQ) 
C/                      FDQ_Firnames.Txt  and FDQ_Convbuff.Txt are moved to 
C/                      the FUT Facility; thus the include files are now 
C/                      FUT_Firnames and FUT_Convbuff. 
C/                      The logging option was never fully implemented. Other
C/                      reporting functions replaced logging.  
C/
C/      H. Wang
C/      STX
C/      Jan. 29, 1991
C/      Reason:
C/                New requirements for FDQ
C/                The calibrator resistor counts used for the conversion
C/                of temperature GRTs will be read from a FEX_AV_CALRS
C/                Reference file
C/  
CH      CHANGE LOG:             New Format for Build 4.2  STX, 10/15/88
CH
CH      Version 4.1.1 10/15/88, SPR 2622, Shirley M. Read, STX
CH              Flag values for missing conversions need changes in FDQ.
CH              The flag values for converted fields which are out of range
CH              need the same changes. GRTs in dwell mode need a flag value
CH              also. The SWG has selected a flag value of -9999 for all cases.
CH              The limit checking algorithms need to bypass setting the quality
CH              flags if the GRT or engineering analog has the flag value.
CH              Interpolation and averaging algorithms need to be modified to
CH              include the flag value. 
CH      Version 4.1.1 10/19/88, SPR 2657, Shirley M. Read, STX
CH              FDQ Convert needs to convert the new LMACs. The conversion
CH              coefficients have now been defined and the include text file
CH              has been expanded to include the LMACs.
CH      Version 4.2.1 01/12/89, SPR 3133, Shirley M. Read, STX
CH              FDQ fails to flag all GRTs in dwell mode. FDQ_Convert needs a
CH              loop over all 16 GRTs for each of the four groups: A High, 
CH              A Low, B High and B Low to set the engineering buffer to the
CH              flag value when the side is in dwell mode.
CH
CH      Version 7.3 SPRs 7852, Glitch Channel PREAMP
CH              H. Wang, STX, 12/11/90
C-------------------------------------------------------------------------
C/
C/      CALLING SEQUENCE:
C/        STATUS = FDQ_CONVERT (HSK_RECS, FIRST_MJ, CAL_LUN,time_range,
C/                 ENG_BUFFS)
C/      
C/      INPUT PARAMETERS:
C/        HSK_RECS(1024)        BYTE    Buffer holding 2 records of HSKP records
C/        FIRST_MJ              I*2     Pointer to one of the 4 major frames held
C/                                      in buffer HSK_RECS
C/        CAL_LUN               I*2     Ave. cal. resistor File log. Unit num. 
C/       *** Also: ***          The common CONVBUFF, which contains the polynomial
C/                              coefficients, GRT lookup tables and calres information.
C/
C/      OUTPUT PARAMETERS:
C/        ENG_BUFFS(1024, 2)    BYTE    Buffer holding 2 records of engineering
C/                                      records encompassing the science record
C/                                      in time; these will later be interpolated
C/
C/      INPUT/OUTPUT FILES:
C/        NONE
C/
C/      INCLUDE FILES USED:
C/        FUT_CONVBUFF.TXT
C/        FUT_FIRNAMES.TXT
C/
C/      SUBROUTINES CALLED:
C/        FDQ_OL_GRT_CVT
C/        FDQ_OL_POLY_CVT
C/      
C/      ERROR HANDLING:
C/        Call Lib$Signal with FDQ defined parameters. A condition handler
C/        was established in FDQ to signal the errors and return control
C/        to the FDQ Program.
C/      
C/      METHOD USED:
C/        The following is the PDL --
C/        If first time then
C/          Read cal. resistor counts from fex_av_calrs
C/        Endif  
C/        Compare the hkp gmt and FEX_av_calrs gmt 
C/        If Hkp GMT not match  the FEX_AV_CALRS GMT
C/        Then
C/          Read cal. resistor counts from fex_av_calrs
C/        Endif  
C/
C/        Move the 2 pertinent (known through FIRST_MJ) mj. frames of
C/               HSKP_RECS into local buffers (so that equivalencing is possible);
C/        Do for both major frames
C/
C/          Do for all poly-type quantities
C/            Using info from CONVBUFF call POLYCVT to convert to engr. units;
C/          Enddo
C/        
C/          Do for all GRTs
C/            Using info from CONVBUFF and fex_av_clars ave. cal. resistor
C/               counts  call GRTCVT to convert to engr. units (temps);
C/          Enddo
c/
C/        Enddo
C/        
        IMPLICIT        NONE

        INCLUDE         '(FUT_CONVBUFF)'

        INCLUDE         '(FUT_FIRNAMES)'
        INCLUDE         'ct$library:ctuser.inc'

        dictionary 'fdq_eng'
        record /fdq_eng/ eng_buffs(2)

        dictionary 'nfs_hkp'
        record /nfs_hkp/  hsk_recs(2)

        dictionary 'fex_av_calrs'
        record /fex_av_calrs/ calrs

        EXTERNAL        FDQ_NORMAL
        EXTERNAL        FDQ_ABERR
        EXTERNAL        FDQ_MJFERR
        EXTERNAL        FDQ_OLPOLYCVER
        EXTERNAL        FDQ_OLGRTCVTER
        EXTERNAL        FDQ_NCALREAD
          
        INTEGER*4       FDQ_OL_POLY_CVT
        INTEGER*4       FDQ_OL_GRT_CVT
        INTEGER*4       Ct_read_arcv
        INTEGER*4       Ct_connect_read
        External        Ct_connect_read 
        INTEGER*4       Ct_close_arcv
        INTEGER*4       Time_gt

        integer*4       RETSTAT         ! Return status
        integer*4       STATUS          ! Status from function
        integer*4       SUCCESS / 1 /, ERROR / 2 /  ! Values for status
        integer*4       ZERO / 0 /
        character*14    fr1_gmt, fr2_gmt, temp_gmt
        Integer*4       temp_time(2) 
        BYTE            i1              ! Local variable to hold the 
                                        ! unconverted byte size counts

        integer*2       i2              ! Local variable to hold the 
                                        ! unconverted int*2  size counts
        INTEGER*2       FIRST_MJ, ISTAT, ct_stat(20)
        INTEGER*2       CAL_LUN 
        INTEGER*2       CALGRP, CNTS, CUR_GRT,
        1               CUR_POLY, MJ,
        1               NBYTES, RES_CNTS(4)
        integer*2       curpoly                 ! Index of array
        integer*2       mjf                     ! Index of array
        Real*4          AVE_RES_CNTS(4)
        INTEGER*2       CAL_HSK_INDX(4)         ! These are the indices into the
        1               / 57, 89, 121, 153 /    ! 216-byte HSK buffer for the start
                                                ! of the 4 cal resistors group
                                                ! CAL_HSK_INDX(1) : for A, low
                                                ! CAL_HSK_INDX(2) : for A, high
                                                ! CAL_HSK_INDX(3) : for B, low
                                                ! CAL_HSK_INDX(4) : for B, high
        integer*4       loc(8)
     &                  /2,3,2,3,6,7,6,7/       ! Current location of temp ctrl
        integer*4       mf(2)                   ! major frame of hsk rec under
                                                ! consideration
        integer*4       rec_num(2)              ! hsk_rec under consideration
                                                ! command words
        integer*4       i,j                     ! counters
        integer*4       place                   ! a counter
        Character*30    time_range
        Character*60    Arc_file
        Integer*4       iostatus,time_lt,ct_binary_to_gmt
        Logical*1       FIRST_CAL/.true./
        Logical*1       Telm_Flag 
        REAL*4          CVALUE, CORR_TEMP
        REAL*4          FV / -9999.0 /          ! Flag value, in this case 
                                                ! dwell mode or no cal. group
        PARAMETER       BAD = 2                 ! Return status
        PARAMETER       GOOD = 1                ! Return status

        integer*4 offset /180/
        integer*4 delta_time(2)
        Integer*4 qoffset(2)
        Integer*4 new_time(2)
        Character*14 new_gmt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                          !
!     Code begins here     !
!                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C       Set return status to success.
        RETSTAT = SUCCESS
        IF (first_cal) then
          Status = ct_read_arcv(, cal_lun,calrs,ct_stat)
          If (ct_stat(1) .ne. ctp_normal) then
            call lib$signal(fdq_Ncalread)
            retstat = error
          else
            If(time_lt(hsk_recs(1).ct_head.time, calrs.ct_head.time)) then
             status= ct_close_arcv(,cal_lun,ct_stat)
             Call Lib$subx(calrs.data_stop_time,calrs.ct_head.time, delta_time)
             Call lib$emul(offset,10000000,0,qoffset)
             Call Lib$addx(delta_time, qoffset,delta_time)
             Call lib$subx(calrs.ct_head.time, delta_time,new_time)
             status = ct_binary_to_gmt(new_time, new_gmt) 
             time_range(1:14)=new_gmt
             arc_file='csdr$firas_ref:fex_av_calrs'//'/'//time_range
             open(unit=cal_lun,file=arc_file,status='old',iostat=iostatus,
        1       useropen=ct_connect_read)
             If (iostatus .ne. 0) then
               retstat=error
               call lib$signal(fdq_Ncalread)
             else
               Status = ct_read_arcv(, cal_lun,calrs,ct_stat)
               If (ct_stat(1) .ne. ctp_normal) then
                 call lib$signal(fdq_Ncalread)
                 retstat = error
               Endif
             Endif
            endif
            first_cal = .false.
          endif     
        Endif       

        IF (first_mj.lt.1 .OR. first_mj.gt.4) THEN      ! ERROR!

          RETSTAT = ERROR
          CALL LIB$SIGNAL(FDQ_MJFERR,%VAL(1),%VAL(first_mj))
        ENDIF

        if(first_mj .eq. 1)then
          rec_num(1) = 1
          mf(1) = 1
          rec_num(2) = 1
          mf(2) = 2
          fr1_gmt=HSk_recs(rec_num(1)).ct_head.gmt
          fr2_gmt=HSk_recs(rec_num(2)).hskp_tail.gmt_mjf2
        elseif(first_mj .eq. 2)then
          rec_num(1) = 1
          mf(1) = 2
          rec_num(2) = 2
          mf(2) = 1
          fr2_gmt=HSk_recs(rec_num(2)).ct_head.gmt
          fr1_gmt=HSk_recs(rec_num(1)).hskp_tail.gmt_mjf2
        elseif(first_mj .eq. 3)then
          rec_num(1) = 2
          mf(1) = 1
          rec_num(2) = 2
          mf(2) = 2
          fr1_gmt=HSk_recs(rec_num(1)).ct_head.gmt
          fr2_gmt=HSk_recs(rec_num(2)).hskp_tail.gmt_mjf2
        else
          rec_num(1) = 2
          mf(1) = 2
          rec_num(2) = 1
          mf(2) = 1
          fr2_gmt=HSk_recs(rec_num(2)).ct_head.gmt
          fr1_gmt=HSk_recs(rec_num(1)).hskp_tail.gmt_mjf2
        endif
        DO      MJ = 1, 2
         IF ( RETSTAT .EQ. SUCCESS ) THEN
          If (mj .eq. 1) temp_gmt = fr1_gmt 
          If (mj .eq. 2) temp_gmt = fr2_gmt 
          Call ct_gmt_to_binary(temp_gmt, temp_time)   
          If (time_gt(temp_time,calrs.data_stop_time)) then
            Status = ct_read_arcv(, cal_lun,calrs,ct_stat)
            If (ct_stat(1) .ne. ctp_normal) then
              call lib$signal(fdq_Ncalread)
              retstat = error
            endif
          endif           
         Endif
         IF ( RETSTAT .EQ. SUCCESS ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !
!    Handle all polynomial conversions in the following DO loop     !
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
          nbytes = 2

          DO    CUR_POLY = 1,8
           IF ( RETSTAT .EQ. SUCCESS ) THEN
            i2 = hsk_recs(rec_num(mj)).frame(mf(mj)).
     &                                hskp_head.stat_monitor_cmd(loc(cur_poly))

!   There are no poly coefficients for the command readback of the temperature
!   controllers. Until a conversion table is available, the code is removed.
!   10/15/88.

            eng_buffs(mj).en_analog.temp_ctrl(cur_poly) = i2
           ENDIF        ! return status is success
          ENDDO         ! CUR_POLY = 1, 8

          nbytes = 1

          do cur_poly = 9,10
           IF ( RETSTAT .EQ. SUCCESS ) THEN
            i1 = hsk_recs(rec_num(mj)).frame(mf(mj)).temps.ipdu(cur_poly-8)
            mjf = MJ
            curpoly = CUR_POLY

            STATUS = FDQ_OL_POLY_CVT ( NBYTES,  ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

            IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
            ENDIF

            eng_buffs(mj).en_analog.ipdu_temp(cur_poly-8) = cvalue
           ENDIF        ! return status is success
          end do

          do cur_poly = 11,14
           IF ( RETSTAT .EQ. SUCCESS ) THEN
            i1 = hsk_recs(rec_num(mj)).frame(mf(mj)).temps.
     &                                                 chan_temp(cur_poly-10)
            mjf = MJ
            curpoly = CUR_POLY

            STATUS = FDQ_OL_POLY_CVT ( NBYTES,  ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

            IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
            ENDIF

            eng_buffs(mj).en_analog.cna_temp(cur_poly-10) = cvalue
           ENDIF        ! return status is success
          end do

          IF ( RETSTAT .EQ. SUCCESS ) THEN
          cur_poly = 15
          i1 = hsk_recs(rec_num(mj)).frame(mf(mj)).temps.drive_box_a
          mjf = MJ
          curpoly = CUR_POLY

          STATUS = FDQ_OL_POLY_CVT ( NBYTES,    ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

            IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
            ENDIF

          eng_buffs(mj).en_analog.dbx_temp(1) = cvalue
          ENDIF ! return status is success

          IF ( RETSTAT .EQ. SUCCESS ) THEN
          cur_poly = 16
          i1 = hsk_recs(rec_num(mj)).frame(mf(mj)).temps.drive_box_b
          mjf = MJ
          curpoly = CUR_POLY

          STATUS = FDQ_OL_POLY_CVT ( NBYTES,    ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

            IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
            ENDIF

          eng_buffs(mj).en_analog.dbx_temp(2) = cvalue
          ENDIF ! return status is success

          do cur_poly = 17,18
           IF ( RETSTAT .EQ. SUCCESS ) THEN
            i1 = hsk_recs(rec_num(mj)).frame(mf(mj)).temps.
     &                                                   stat_mon(cur_poly-16)
            mjf = MJ
            curpoly = CUR_POLY

            STATUS = FDQ_OL_POLY_CVT ( NBYTES,  ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

            IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
            ENDIF

            eng_buffs(mj).en_analog.stat_mon_temp(cur_poly-16) = cvalue
           ENDIF        ! return status is success
          end do

          IF ( RETSTAT .EQ. SUCCESS ) THEN
          cur_poly = 19
          i1 = hsk_recs(rec_num(mj)).frame(mf(mj)).temps.chan_pre_amp
          mjf = MJ
          curpoly = CUR_POLY

          STATUS = FDQ_OL_POLY_CVT ( NBYTES,    ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

          IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
          ENDIF

          eng_buffs(mj).en_analog.pamp_chan = cvalue
          ENDIF ! return status is success

          IF ( RETSTAT .EQ. SUCCESS ) THEN
          cur_poly = 20
          i1 = hsk_recs(rec_num(mj)).frame(mf(mj)).temps.optical_preamp
          mjf = MJ
          curpoly = CUR_POLY

          STATUS = FDQ_OL_POLY_CVT ( NBYTES,    ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

           IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
           ENDIF

          eng_buffs(mj).en_analog.pamp_op = cvalue
          ENDIF ! return status is success

          do cur_poly = 21,22
           IF ( RETSTAT .EQ. SUCCESS ) THEN
            i1 = hsk_recs(rec_num(mj)).frame(mf(mj)).temps.
     &                                         hot_spot_current(cur_poly-20)
            mjf = MJ
            curpoly = CUR_POLY

            STATUS = FDQ_OL_POLY_CVT ( NBYTES,  ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

            IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
            ENDIF

            eng_buffs(mj).en_analog.hot_spot(cur_poly-20) = cvalue
           ENDIF        ! return status is success
          end do

          do cur_poly = 23,24
           IF ( RETSTAT .EQ. SUCCESS ) THEN
            i1 = hsk_recs(rec_num(mj)).frame(mf(mj)).v_and_i.
     &                                          mtm_cal_motor(cur_poly-22)
            mjf = MJ
            curpoly = CUR_POLY

            STATUS = FDQ_OL_POLY_CVT ( NBYTES,  ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

            IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
            ENDIF

            eng_buffs(mj).en_analog.mtm_cal_mtr(cur_poly-22) = cvalue
           ENDIF        ! return status is success
          end do

          do cur_poly = 25,26
           IF ( RETSTAT .EQ. SUCCESS ) THEN
            i1 = hsk_recs(rec_num(mj)).frame(mf(mj)).hskp_head.
     &                                            lvdt_stat(cur_poly-24)
            mjf = MJ
            curpoly = CUR_POLY

            STATUS = FDQ_OL_POLY_CVT ( NBYTES,  ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

            IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
            ENDIF

            eng_buffs(mj).en_analog.mtm_pos(cur_poly-24) = cvalue
           ENDIF        ! return status is success
          end do

          do cur_poly = 27,30
           IF ( RETSTAT .EQ. SUCCESS ) THEN
            i1 = hsk_recs(rec_num(mj)).frame(mf(mj)).v_and_i.
     &                                             bol_bias_volt(cur_poly-26)
            mjf = MJ
            curpoly = CUR_POLY

            STATUS = FDQ_OL_POLY_CVT ( NBYTES,  ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

            IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
            ENDIF

            eng_buffs(mj).en_analog.bol_volt(cur_poly-26) = cvalue
           ENDIF        ! return status is success
          end do

          do cur_poly = 31,50
           IF ( RETSTAT .EQ. SUCCESS ) THEN
            i1 = hsk_recs(rec_num(mj)).frame(mf(mj)).v_and_i.
     &                                            ipdu_volt(cur_poly-30)
            mjf = MJ
            curpoly = CUR_POLY

            STATUS = FDQ_OL_POLY_CVT ( NBYTES,  ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

            IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
            ENDIF

            eng_buffs(mj).en_analog.ipdu_volt(cur_poly-30) = cvalue
           ENDIF        ! return status is success
          end do

          do cur_poly = 51,62
           IF ( RETSTAT .EQ. SUCCESS ) THEN
            i1 = hsk_recs(rec_num(mj)).frame(mf(mj)).v_and_i.
     &                                                     ipdu_amp(cur_poly-50)
            mjf = MJ
            curpoly = CUR_POLY

            STATUS = FDQ_OL_POLY_CVT ( NBYTES,  ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

            IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
            ENDIF

            eng_buffs(mj).en_analog.ipdu_amp(cur_poly-50) = cvalue
           ENDIF        ! return status is success
          end do

C       LMAC temperatures added 10/19/88.

          IF ( RETSTAT .EQ. SUCCESS ) THEN
            cur_poly = 63  
            i1 = hsk_recs(rec_num(mj)).mj_frm(mf(mj)).lmac_analog_temp
            mjf = MJ
            curpoly = CUR_POLY

            STATUS = FDQ_OL_POLY_CVT ( NBYTES,  ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

            IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
            ENDIF

            eng_buffs(mj).en_tail.lmac_analog_temp = cvalue
          ENDIF ! return status is success

          IF ( RETSTAT .EQ. SUCCESS ) THEN
            cur_poly = 64  
            i1 = hsk_recs(rec_num(mj)).mj_frm(mf(mj)).lmac_digital_temp
            mjf = MJ
            curpoly = CUR_POLY

            STATUS = FDQ_OL_POLY_CVT ( NBYTES,  ! # bytes in count
        1               i1,                     ! Value in counts to be converted
        1               POLY_NC(CUR_POLY),      ! # coef.
        1               POLY_COEF(1, CUR_POLY), ! coefficients
        1               CVALUE,                 ! Converted value
        1               mjf,                    ! Major frame ( 1 or 2 )
        1               curpoly )               ! Polynomial array index

            IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
              IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
              CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
              RETSTAT = ERROR
            ENDIF

            eng_buffs(mj).en_tail.lmac_digital_temp = cvalue
          ENDIF ! return status is success

          place = 0

          do i = 1,4               ! loop over side a & b and the current modes

           If (((i .eq. 1 .or. i .eq. 2) .and.   ! Check for dwell mode 10/17/88
        1       (eng_buffs(mj).en_stat.dwell_stat(1) .eq. 0)) .OR.
        2      ((i .eq. 3 .or. i .eq. 4) .and.                ! Not dwell mode
        3       (eng_buffs(mj).en_stat.dwell_stat(2) .eq. 0))) Then

           IF ( RETSTAT .EQ. SUCCESS ) THEN
            do cur_grt = 1,16
             IF ( RETSTAT .EQ. SUCCESS ) THEN
              if(cur_grt .eq. 1)then
                do j = 1,4
                  res_cnts(j) = hsk_recs(rec_num(mj)).frame(mf(mj)).temps.
     &                                             side_amp(i).cal_resist(j)
                  eng_buffs(mj).en_analog.grt(j + 10 + 16*(i-1)) =
     &                  res_cnts(j)     ! Cal. resistors are positions 11 - 14
                                        ! in FDQ_ENG for each group of GRTs
                end do
              endif
              IF (I .eq. 1) then
                DO j=1,4
                  AVE_RES_CNTS(j) = CALRS.calres_ave_a_lo(j)
                ENDDO 
              ENDIF
              IF (I .eq. 2) then
                DO j=1,4
                  AVE_RES_CNTS(j) = CALRS.calres_ave_a_Hi(j)
                ENDDO 
              ENDIF
              IF (I .eq. 3) then
                DO j=1,4
                  AVE_RES_CNTS(j) = CALRS.calres_ave_B_lo(j)
                ENDDO 
              ENDIF
              IF (I .eq. 4) then
                DO j=1,4
                  AVE_RES_CNTS(j) = CALRS.calres_ave_B_Hi(j)
                ENDDO 
              ENDIF
              if(cur_grt.le.10 .or. cur_grt.ge.15)then    !exclude cal resistors

!               Place takes on values from 1 to 48

                place = place + 1

!      Grt_calgrp is determined in load conversions and passed through convbuff

                calgrp = grt_calgrp(place)

                if(calgrp .ne. 0)then

                  cnts = hsk_recs(rec_num(mj)).frame(mf(mj)).temps.
     &                                              side_amp(i).grt(cur_grt)
                  If (hsk_recs(rec_num(mj)).mj_frm(mf(mj)).tlm_qual_maj_frm 
        1             .ne. 0) then
                    telm_flag =.true.
                  else
                    telm_flag = .false.
                  Endif 
                  mjf = MJ
                  If (((i .eq. 3 .or. i .eq. 4) .and. (cur_grt .eq. 16)) .or.
        1           ((i .eq.2 .or. i .eq. 4) .and. (cur_grt .ge. 6 .and. cur_grt 
        1           .le. 9)) .or.((i .eq. 3 .or. i .eq. 4) .and. ( 
        1            cur_grt .eq. 1))) 
        1             then
                    corr_temp = FV
                  Else 
                   STATUS = FDQ_OL_GRT_CVT (CNTS,       ! GRT_CNT
        1                       AVE_RES_CNTS,           ! Count of cal. resistors
        1                       CALGRP,                 ! Cal group
        1                       GRT_NPTS(place),        ! # points in lookup table
        1                       GRT_XTAB(1, place),     ! Resistance values in table
        1                       GRT_YTAB(1, place),     ! Temperature values in lookup table
        1                       GRT_COEF(1, 1, place),  ! Quad. coef. for fits
        1                       CORR_TEMP,              ! Corrected GRT temperature
        1                       ISTAT,                  ! Return status for data
        1                       mjf,                    ! Major frame ( 1 or 2 )
        1                       telm_flag,      
        1                       place)                  ! Array index of GRT
!
!     Check for bad status from FDQ_OL_GRT_CVT and report which GRT (Added 11-6-85)
!
                  IF ( STATUS .NE. %loc(FDQ_NORMAL)) THEN
                    IF ( STATUS .EQ. %loc(FDQ_ABERR)) STATUS = ZERO
                    CALL LIB$SIGNAL(FDQ_OLPOLYCVER,%VAL(1),%VAL(STATUS))
                    RETSTAT = ERROR
                  ELSEIF ( ISTAT .NE. GOOD ) THEN
                  ENDIF
                 Endif  
                else

                  corr_temp = FV

                endif
!
!     Now move corrected temperature to the converted buffer ENG_BUFFS
!
                eng_buffs(mj).en_analog.grt(cur_grt + 16*(i-1)) = corr_temp


              endif
            ENDIF       ! return status is success
            end do       !cur_grt
           ENDIF        ! return status is success
           Else          ! Dwell Mode 10/17/88
            do cur_grt = 1, 16
             eng_buffs(mj).en_analog.grt(cur_grt + 16*(i-1)) = FV
            enddo
           Endif
          end do         !i = 1,4
         ENDIF           !Return status is success
        end do           !mj = 1,2


!       Set function to return status

        IF (RETSTAT.EQ.SUCCESS) THEN
          FDQ_CONVERT = %loc(FDQ_NORMAL)
        ELSE
          FDQ_CONVERT = %loc(FDQ_ABERR)
        ENDIF

        RETURN
        END

