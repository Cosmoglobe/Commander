	INTEGER*4 FUNCTION FDQ_INTERPOLATE (HSK_RECS, FIRST_MJ, 
	1  IFG_TIME, ENG_BUFFS, ENG_REC, HSK_ANALOG, XCAL_SET, TELM_FLAG,
	1	grtwt,grttrans)
C/
C/	PROGRAM NAME:
C/	  FDQ_INTERPOLATE
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine takes as input the 2-record Housekeeping buffer,
C/	  the pointer indicating which of the 4 major frames in the buffer,
C/	  the converted Engineering buffers emanating from the 2 Housekeeping
C/	  major frames, as well as the Science record through named common
C/	  SCIREC, and produce as output the interpolated Engineering record.
C/	  It also copy the permanent status words from the Housekeeping
C/	  buffer to the Engineering record.
C/
C/	AUTHOR:
C/	  Edwin H. Fung
C/	  GSFC
C/	  October 25, 1985 (2:40 am!)
C/
C/	MODIFIED BY:
C/	  Edwin H. Fung
C/	  GSFC
C/	  October 26, 1985
C/	  REASON:	To basically complete the subroutine:  Add logic
C/			to handle Temperature Controllers and move the
C/			housekeeping status words.
C/
C/	  E. FUNG
C/	  GSFC
C/	  November 14, 1985
C/	  REASON:	To correct bug in converting 100-nanoseconds units
C/			to minor frames.
C/
C/	  E. FUNG
C/	  GSFC
C/	  November 16, 1985
C/	  REASON:	To calculate interpolated counts for calibrator
C/			resistors.
C/
C/	  E. FUNG
C/	  GSFC
C/	  November 18, 1985
C/	  REASON:	To correct bug introduced in the last modification:
C/			Have to convert cal resistors counts from integer
C/			to floating point data type!
C/
C/	  E.FUNG
C/	  GSFC
C/	  DEC. 7, 1985
C/	  REASON:	To change the writing to the log file to D-lines.
C/
C/	  E. FUNG
C/	  GSFC
C/	  FEB. 25, 1986
C/	  REASON:	To change the indices into fields into the ENG record
C/			because of the new format for the ENG archive.
C/
C/	  E. FUNG
C/	  GSFC
C/	  MARCH 30, 1986
C/	  REASON:	To make changes to allow DQ to write the Index.  Among
C/			additions are the moving of the first major frame of
C/			HSK to the new passed parameter HSK_ANALOG.  The decision
C/			NOT to interpolate the 2 bracketing HSK major frame
C/			analog values is made because of the additional overhead
C/			involved in the data conversion (from integer to real
C/			and then to integer again) as well as the way the VAX
C/			handle byte integers as signed rather than unsigned,
C/			so there would be more overhead to handle this by software.
C/
C/        J. Durachta
C/        ARC
C/        July 27, 1987
C/        REASON:       Alter NFS_HDF to NFS_HKP for notch filters.
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 6, 1988
C/	  REASON: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal. Removed quality flag checking
C/			from this routine to FDQ_Get_Qualflags. Added user
C/	                selected option of setting Xcal or not.
C/
C/	  Shirley M. Read
C/	  STX
C/	  June 3, 1988
C/	  REASON: 	Changed the algorithm for storing the A and B power
C/			status in the engineering record. The telemetry 
C/			Subcom and consequently the housekeeeping record have
C/			been modified. There are now 2 bits for each power
C/			status component instead of 1.
C/
C/       H. Wang
C/       STX
C/       Jan. 29, 1991
C/       Reason:        New requirements for FDQ
C/                      1) The instrument temperature changes between bracketing
C/                         major frames will be calculated and stored in the
C/                         engineering file(FDQ_ENG)
C/                      2) Gain and fakeit will no longer be figured out by
C/                         this routine
C/ 
CH	CHANGE LOG:		New Format for Build 4.2  STX, 10/15/88
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
CH      Version 4.1.1 10/19/88, SPR 2657, Shirley M. Read, STX
CH	   	FDQ Convert needs to convert the new LMACs. The conversion
CH	        coefficients have now been defined and the include text file
CH	        has been expanded to include the LMACs.
CH
CH	Version 4.2.1 02/09/89, SPR 2700, Shirley M. Read, STX
CH		FDQ falls over for out of order and bad time ranges.The FIRAS
CH	 	Stripper now writes the science records using the telemetry
CH		minor frame time at transmit for the primary key in order to
CH		ensure a valid time tag. A new FIRAS preprocessor reads the raw
CH		science records, computes the midpoint of collect time and flags
CH		any records with bad midpoint times. FDQ must now access the
CH		computed collect time from the new field in the science record
CH	        in order to group the four science channels with the bracketing
CH		housekeeping frames and set the data quality summary flag when
Ch		the badtime flag is set. The bad record is written back to the
CH		archive immediately so that the next record may be considered
CH	        for the group. The average midpoint of collect times is still
CH		used as the time tag of the FDQ engineering record and is thus
CH		stored in the science record as the cross reference. However,
CH		the transmit times of the science records are stored in the 
CH		FDQ engineering record. FDQ must now enter the science gain into
CH		each science record from the bracketing housekeeping frames. The
CH		two LS bits must now be checked in the certification mask when
CH		running FDQ. The bit must be set for FPP but not for FDQ.
CH
CH	Version 4.4.1 05/30/89, SPR 3916, Shirley M. Read, STX
CH		The two telemetry bits indicating the XCAL position have been
CH		reversed from their original meaning. The MSB "XCAL Position
CH	        In" is set to 0 when the XCAL is put in the horn and the LSB
CH		"XCAL Position Out" is set to 1. The two bits produce a result
CH		of 1 for CAL Position. The reverse is true for stow position
CH	        (sky data) with a resulting value of 2. The command line set
CH	        with XCAL_SET has been updated to set the bits for the new
CH		convention. This command line option will probably never be
CH		used again since the telemetry readout is now available. The
CH		simplest fix is to modify the setting rather than remove the
CH		option from the command line.
CH
C-----------------------------------------------------------------------------
C/
C/	CALLING SEQUENCE:
C/	  Status = FDQ_INTERPOLATE (HSK_RECS, FIRST_MJ, IFG_TIME, ENG_BUFFS,
C/	           ENG_REC, HSK_ANALOG, XCAL_SET,TELM_FLAG,grtwt,grttrans)
C/
C/	INPUT PARAMETERS:
C/	  HSK_RECS(1024)	BYTE	Buffer holding 2 Housekeeping records
C/	  FIRST_MJ		I*2	Pointer to which of the 4 major frames
C/					in HSK_RECS.
C/						1 = 1st maj. fr., 1st record
C/						2 = 2nd maj. fr., 1st record
C/						3 = 1st maj. fr., 2nd record
C/						4 = 2nd maj. fr., 2nd record
C/	  IFG_TIME(2)		I*4	Time of the midpoint of the collect
C/					cycle in VAX binary format.
C/	  ENG_BUFFS(1024,2)	BYTE	Buffers containing the converted values
C/					of the 2 major frames pointed to by FIRST_MJ
C/	  XCAL_SET              L*1     Flag indicating whether to set Xcal or
C/					not. User selected option at runtime.
C/        GRTWT                 record  FEX_GRTRAWWT       
C/        GRTTRANS              record  FEX_GRTtrans       
C/       
C/	  *** ALSO: Science record SCI_REC in named common SCIREC ***
C/	
C/	OUTPUT PARAMETERS:
C/	  ENG_REC(1024)		BYTE	The interpolated Engineering record,
C/					ready to go out to the archive (with
C/					the exceptions of the quality flags)
C/	  HSK_ANALOG(52)	BYTE	The HSK analog values that will remain
C/					as counts in the Index that DQ will
C/					write (uninterpolated; from 1st major frame)
C/	  TELM_FLAG             I*2     Telemetry flag of Engineering record
C/                                      0 = telemetry quality good
C/                                      1 = telemetry quality bad 
C/	INPUT/OUTPUT FILES:
C/	  NONE
C/
C/	INCLUDE FILES USED:
C/
C/	SUBROUTINES CALLED:
C/	  GMT_TO_BINARY (COBETRIEVE routine)
C/	  LIB$MOVC3
C/	  LIB$SUBX
C/
C/	ERROR HANDLING:
C/	  TBD
C/
C/	METHOD USED:
C/	  The following is the PDL --
C/	  
C/
C/        Calculate the instrument temperature changes between
C/         bracketing major frames 
C/
C/	  Find the difference between the Science hskp_head time and the
C/	  	GMT time of the first major frame of Housekeeping data;
C/	  Convert the difference from seconds to minor frames;
C/	  Adjust the difference (in minor frames) according to the
C/	  	information given in the channel hskp_head in the Science record;
C/	  Use the resulting difference as the "weight" in interpolating the
C/	  	converted Engineering values;
C/	  Move other Housekeeping data (2-level subcoms and statuses) to
C/		the interpolated Engineering record buffer;
C/	  Return;
C/	  End.
C/

	IMPLICIT	NONE

        include 	'ct$library:ctuser.inc'
	include		'($ssdef)'
        dictionary 'fdq_eng'
        record /fdq_eng/   eng_buffs(2)
        record /fdq_eng/   eng_rec

        dictionary 'nfs_hkp'
        record /nfs_hkp/    hsk_recs(2)

        dictionary 'fut_enganlg'
        record /fut_enganlg/    rawsigmas
        dictionary 'fex_grtrawwt'
        dictionary 'fex_grttrans'
	record /fex_grtrawwt/    grtwt
	record /fex_grttrans/    grttrans

        byte		HSK_ANALOG(52)

	INTEGER*2	FIRST_MJ, TELM_FLAG

	INTEGER*4	IFG_TIME(2)

        logical*1       XCAL_SET


	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_ERSUBX
	INTEGER*4       LIB$SUBX	! System function - subtract quadword
	
	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status

        integer*4       i,j             ! counters
        integer*4       mode            ! counter to loop over
                                        ! side-current configurations
        integer*4       act_mj_a(2)     ! actual a side major frame
        integer*4       act_mj_b(2)     ! actual b side major frame

	INTEGER*2	bn(2),mj(2),CUR,
	1		ENG_A(2), ENG_B(2),
	1		FIRAS_MJ_A(2),
	1		FIRAS_MJ_B(2),
	1		I1, I2, MJFR, 
	1		STAT_WD,
	1		zeroes(4) / 4*0 /,
	1               hsk_a_dwell,
	1               hsk_b_dwell


	INTEGER*4	DIFF(2), DIFF12(2), COMBSWITCH,
	1		HSK_TIM1(2), HSK_TIM2(2)

	REAL*4		REAL1, REAL2,
	1		RESULT, WEIGHT

	REAL*4		FV / -9999.0 /
        REAL*4          outtemps_a1(10),outtemps_a2(10),outsigmas(10)
        REAL*4          outtemps_b1(10),outtemps_b2(10),outsigs(10)
        INTEGER*2       IOUTTEMPS_A(10), IOUTTEMPS_B(10)	
        Byte            Hsk_A_LVDT,HSK_B_LVDT       

        logical*4       singlifg/.true./

        character*1     ans

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


C	Set return status to success.

	RETSTAT = SUCCESS

        TELM_FLAG = 0
!
!     Get the start time of the first Housekeeping major frame
!     from HSK_RECS and pointer FIRST_MJ
!

        if(first_mj .eq. 1)then       ! First major frame in buffer 1 of HSK_RECS
          bn(1) = 1
          bn(2) = 1
          mj(1) = 1
          mj(2) = 2
          hsk_tim1(1) = hsk_recs(bn(1)).ct_head.time(1)
          hsk_tim1(2) = hsk_recs(bn(1)).ct_head.time(2)
          call gmt_to_binary(%ref(hsk_recs(bn(2)).hskp_tail.gmt_mjf2),hsk_tim2)
        elseif(first_mj .eq. 2)then   ! First major frame in buffer 2 of HSK_RECS
          bn(1) = 1
          bn(2) = 2
          mj(1) = 2
          mj(2) = 1
          hsk_tim2(1) = hsk_recs(bn(2)).ct_head.time(1)
          hsk_tim2(2) = hsk_recs(bn(2)).ct_head.time(2)
          call gmt_to_binary(%ref(hsk_recs(bn(1)).hskp_tail.gmt_mjf2),hsk_tim1)
        elseif(first_mj .eq. 3)then   ! First major frame in buffer 3 of HSK_RECS
          bn(1) = 2
          bn(2) = 2
          mj(1) = 1
          mj(2) = 2
          hsk_tim1(1) = hsk_recs(bn(1)).ct_head.time(1)
          hsk_tim1(2) = hsk_recs(bn(1)).ct_head.time(2)
          call gmt_to_binary(%ref(hsk_recs(bn(2)).hskp_tail.gmt_mjf2),hsk_tim2)
        else                          ! First major frame in buffer 4 of HSK_RECS
          bn(1) = 2
          bn(2) = 1
          mj(1) = 2
          mj(2) = 1
          hsk_tim2(1) = hsk_recs(bn(2)).ct_head.time(1)
          hsk_tim2(2) = hsk_recs(bn(2)).ct_head.time(2)
          call gmt_to_binary(%ref(hsk_recs(bn(1)).hskp_tail.gmt_mjf2),hsk_tim1)
        endif
        IF ((Hsk_recs(bn(1)).mj_frm(mj(1)).tlm_qual_maj_frm .eq. 1) .or.
	1            (Hsk_recs(bn(2)).mj_frm(mj(2)).tlm_qual_maj_frm .eq. 1)) then
           telm_flag = 1
        endif
!
!     calculate the instrument temperature changes between
!     bracketing major frames 
!
         combswitch = 1   
         Call fut_temp_list(ENG_buffs(1).en_analog,rawsigmas,grtwt,grttrans,
	1	combswitch,singlifg,outtemps_a1,outsigmas,outsigs)
         combswitch = 2   
         Call fut_temp_list(ENG_buffs(1).en_analog,rawsigmas,grtwt,grttrans,
	1	combswitch,singlifg,outtemps_b1,outsigmas,outsigs)
         combswitch = 1   
         Call fut_temp_list(ENG_buffs(2).en_analog,rawsigmas,grtwt,grttrans,
	1	combswitch,singlifg,outtemps_a2,outsigmas,outsigs)
         combswitch = 2   
         Call fut_temp_list(ENG_buffs(2).en_analog,rawsigmas,grtwt,grttrans,
	1	combswitch,singlifg,outtemps_b2,outsigmas,outsigs)
       
         DO i=1, 10
           If ((outtemps_a2(i) .eq. fv) .or. (outtemps_a1(i) .eq. fv)) then
             Iouttemps_a(i) = -9999
           else
             Iouttemps_a(i)=ININT(10000*(outtemps_a2(i) - outtemps_a1(i)))
           Endif  
           If ((outtemps_b2(i) .eq. fv) .or. (outtemps_b1(i) .eq. fv)) then
             Iouttemps_b(i) = -9999
           else
             Iouttemps_b(i)=ININT(10000*(outtemps_b2(i) - outtemps_b1(i)))
           Endif  
         ENDDO
C A side         
         Eng_rec.en_tempdiff(1).xcal = iouttemps_a(1)               
         Eng_rec.en_tempdiff(1).ical = iouttemps_a(2)               
         Eng_rec.en_tempdiff(1).skyhorn = iouttemps_a(3)               
         Eng_rec.en_tempdiff(1).refhorn = iouttemps_a(4)               
         Eng_rec.en_tempdiff(1).dihedral = iouttemps_a(5)               
         Eng_rec.en_tempdiff(1).collimator_mirror= iouttemps_a(6)
         Eng_rec.en_tempdiff(1).bol_assem(1) = iouttemps_a(7)       
         Eng_rec.en_tempdiff(1).bol_assem(2) = iouttemps_a(8)       
         Eng_rec.en_tempdiff(1).bol_assem(3) = iouttemps_a(9)       
         Eng_rec.en_tempdiff(1).bol_assem(4) = iouttemps_a(10)       
C  B side
         Eng_rec.en_tempdiff(2).xcal = iouttemps_b(1)               
         Eng_rec.en_tempdiff(2).ical = iouttemps_b(2)               
         Eng_rec.en_tempdiff(2).skyhorn = iouttemps_b(3)               
         Eng_rec.en_tempdiff(2).refhorn = iouttemps_b(4)               
         Eng_rec.en_tempdiff(2).dihedral = iouttemps_b(5)               
         Eng_rec.en_tempdiff(2).collimator_mirror= iouttemps_b(6)
         Eng_rec.en_tempdiff(2).bol_assem(1) = iouttemps_b(7)       
         Eng_rec.en_tempdiff(2).bol_assem(2) = iouttemps_b(8)       
         Eng_rec.en_tempdiff(2).bol_assem(3) = iouttemps_b(9)       
         Eng_rec.en_tempdiff(2).bol_assem(4) = iouttemps_b(10)       
               
!
!     Move the HSK analog values (except for the GRT's and the TC's) to
!     passed parameter HSK_ANALOG so as to put these into the index 
!     record buffer when DQ write the Index File.
!
!     There are 52 of these analog values:
!	IPDU temps		 2
!	Drive box temps		 2
!	Pre-Amps temps		 2
!	Channel temps		 4
!	Status Monitor temps	 2
!	Hot spot heater curr.	 2
!	IPDU voltages		20
!	IPDU currents		12
!	MTM cal motor		 2
!	Bolometer bias voltages	 4
!

	call lib$movc3 (52, hsk_recs(bn(1)).frame(mj(1)).temps.ipdu(1),
     &                                                           hsk_analog)

!	Store the A and B power status.

	do	i =1, 2
          eng_rec.en_stat.power_A_status(i) =
	1   hsk_recs(bn(i)).mj_frm(mj(i)).power_A_status
	enddo
	do	i =1, 2
          eng_rec.en_stat.power_B_status(i) =
	1   hsk_recs(bn(i)).mj_frm(mj(i)).power_B_status
	enddo

!
!     Put the difference of HSK_TIM1 and HSK_TIM2 in DIFF12
!     and the difference of IFG_TIME (at midpoint of collect cycle)
!     and HSK_TIM1 (midpoint of 1st housekeeping major frame)
!     in DIFF and calculate the relative weight of the 2 major frames
!

	RETSTAT = LIB$SUBX (HSK_TIM2, HSK_TIM1, DIFF12)
	IF (RETSTAT .NE. SS$_NORMAL) THEN
	  CALL LIB$SIGNAL(FDQ_ERSUBX,%VAL(1),%VAL(RETSTAT))
	  RETSTAT = ERROR
	ENDIF

	IF (RETSTAT .EQ. SUCCESS ) THEN
	  RETSTAT = LIB$SUBX (IFG_TIME, HSK_TIM1, DIFF)
	  IF (RETSTAT .NE. SS$_NORMAL) THEN
	    CALL LIB$SIGNAL(FDQ_ERSUBX,%VAL(1),%VAL(RETSTAT))
	    RETSTAT = ERROR
	  ENDIF
	ENDIF

C	Check the status before proceeding.

	IF ( RETSTAT .EQ. SUCCESS ) THEN

	WEIGHT = FLOAT(DIFF(1)) / FLOAT(DIFF12(1))

        
	do i = 1,64                                ! The 16 cal resistors are included
	  If (( eng_buffs(1).en_analog.grt(i) .eq. FV ) .or. 
	1     ( eng_buffs(2).en_analog.grt(i) .eq. FV )) Then
	    eng_rec.en_analog.grt(i) = FV                  
          Else                                        ! to make the logic simple.They
            real1 = eng_buffs(1).en_analog.grt(i)  ! should all be 0 in eng_bufs.
            real2 = eng_buffs(2).en_analog.grt(i)  ! Also, the 8 TC's are 2 level subcoms 
                                                      ! and should not be intepollated.
            result = real1 + weight*(real2 - real1)   ! Again, their inclusion makes the logic
                                                      ! simple; they will be corrected in the 
            eng_rec.en_analog.grt(i) = result      ! next do loop.
          Endif
        end do
	do i = 1,62                                ! The 16 cal resistors are included
	  If (( eng_buffs(1).en_analog.group1(i) .eq. FV ) .or. 
	1     ( eng_buffs(2).en_analog.group1(i) .eq. FV )) Then
	    eng_rec.en_analog.group1(i) = FV                  
          Else                                        ! to make the logic simple.They
            real1 = eng_buffs(1).en_analog.group1(i)  ! should all be 0 in eng_bufs.
            real2 = eng_buffs(2).en_analog.group1(i)  ! Also, the 8 TC's are 2 level subcoms 
                                                      ! and should not be intepollated.
            result = real1 + weight*(real2 - real1)   ! Again, their inclusion makes the logic
                                                      ! simple; they will be corrected in the 
            eng_rec.en_analog.group1(i) = result      ! next do loop.
          Endif
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								!
!     Now put the interpolated calibrator resistors' counts     !
!     into the engineering record buffer			!
!								!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
        do mode = 1,4       ! 2 sides, 2 current modes
          do i = 1,4       ! 4 cal resistors per group

            i1 = hsk_recs(bn(1)).frame(mj(1)).temps.side_amp(mode).cal_resist(i)
            i2 = hsk_recs(bn(2)).frame(mj(2)).temps.side_amp(mode).cal_resist(i)

            result = float(i1) + weight*(float(i2) - float(i1))

            eng_rec.en_analog.grt(10+i+16*(mode-1)) = result

          end do
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                !
!     Interpolate the converted LMACs for the Eng_Rec. 10/20/88  !
!                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  If (( eng_buffs(1).en_tail.lmac_analog_temp .eq. FV ) .or. 
	1     ( eng_buffs(2).en_tail.lmac_analog_temp .eq. FV )) Then
	    eng_rec.en_tail.lmac_analog_temp = FV                  
          Else                              
            real1 = eng_buffs(1).en_tail.lmac_analog_temp
            real2 = eng_buffs(2).en_tail.lmac_analog_temp
            result = real1 + weight*(real2 - real1) 
            eng_rec.en_tail.lmac_analog_temp = result   
          Endif
	  If (( eng_buffs(1).en_tail.lmac_digital_temp .eq. FV ) .or. 
	1     ( eng_buffs(2).en_tail.lmac_digital_temp .eq. FV )) Then
	    eng_rec.en_tail.lmac_digital_temp = FV                  
          Else                              
            real1 = eng_buffs(1).en_tail.lmac_digital_temp
            real2 = eng_buffs(2).en_tail.lmac_digital_temp
            result = real1 + weight*(real2 - real1) 
            eng_rec.en_tail.lmac_digital_temp = result   
          Endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								    !
!     Now the temperature controllers:				    !
!     They have already been converted and are in ENG_BUFFS,	    !
!     but which of the 2 buffers depend on whether which one is	    !
!     "FIRAS major frame 1" and which is "FIRAS major frame 2";	    !
!     They may or may not correspond with the S/C major frames!	    !
!     Check which of the 2 buffers is FIRAS Major Frame 1 and	    !
!     which is FIRAS Major Frame 2.  It is possible that within	    !
!     the same S/C major frame, Side A could be at FIRAS Major      !
!     Frame 1 while Side B is at FIRAS Major Frame 2!		    !
!								    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        stat_wd = hsk_recs(bn(1)).frame(mj(1)).hskp_head.stat_monitor_cmd(1)
	IF (BTEST(STAT_WD, 14)) THEN

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

	  FIRAS_MJ_A(1) = bn(2) 	! The 2nd (S/C) major frame is
	  FIRAS_MJ_A(2) = bn(1)  	! actually FIRAS major frame 1
          act_mj_a(1) = mj(2)             ! adjust associated frame
          act_mj_a(2) = mj(1)
!
!     FIRAS_MJ_A(1) is the buff num in the HSK_RECS buffer for FIRAS major frame 1
!     FIRAS_MJ_A(2) is the buff num in the HSK_RECS buffer for FIRAS major frame 2
!
	  ENG_A(1) = 2		! ENG_A(1) indicates which of ENG_BUFFS is FIRAS major frame 1
	  ENG_A(2) = 1		! ENG_A(2) indicates which of ENG_BUFFS is FIRAS major frame 2

	ELSE

	  FIRAS_MJ_A(1) = bn(1) 	! The first S/C major frame is
	  FIRAS_MJ_A(2) = bn(2) 	! also FIRAS major frame 1
          act_mj_a(1) = mj(1)             ! adjust associated frame
          act_mj_a(2) = mj(2)
	  ENG_A(1) = 1
	  ENG_A(2) = 2

	ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!					   !
!     Now do the same check for Side B     !
!					   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        stat_wd = hsk_recs(bn(1)).frame(mj(1)).hskp_head.stat_monitor_cmd(5)
	IF (BTEST(STAT_WD, 14)) THEN

	  FIRAS_MJ_B(1) = bn(2) 	! First S/C major frame
	  FIRAS_MJ_B(2) = bn(1) 	! is FIRAS major frame 2
          act_mj_b(1) = mj(2)
          act_mj_b(2) = mj(1)
	  ENG_B(1) = 2
	  ENG_B(2) = 1

	ELSE

	  FIRAS_MJ_B(1) = bn(1) 	! First S/C major frame
	  FIRAS_MJ_B(2) = bn(2) 	! is FIRAS major frame 1
          act_mj_b(1) = mj(1)
          act_mj_b(2) = mj(2)
	  ENG_B(1) = 1
	  ENG_B(2) = 2

	ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									    !
!     Now that we have figured out where the FIRAS major frames are,	    !
!     we can intelligently move the converted T.C. values and the 	    !
!     Command Status Words in the correct positions in the Engr.	    !
!     buffer.								    !
!									    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do i = 1,2
          eng_rec.en_analog.temp_ctrl(i) = 
     &                          eng_buffs(eng_a(1)).en_analog.temp_ctrl(i)
          eng_rec.en_analog.temp_ctrl(i+2) = 
     &                          eng_buffs(eng_a(2)).en_analog.temp_ctrl(i+2)
          eng_rec.en_analog.temp_ctrl(i+4) = 
     &                          eng_buffs(eng_b(1)).en_analog.temp_ctrl(i+4)
          eng_rec.en_analog.temp_ctrl(i+6) = 
     &                          eng_buffs(eng_b(2)).en_analog.temp_ctrl(i+6)
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     									  !
!     Now the status words; start with the command status words, they     !
!     also depend on whether it is FIRAS major frame 1 or 2		  !
!     									  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do i = 1,2
          do j = 1,4

! Move A-side status words, 8 bytes from FIRAS MJ 1 and 8 from FIRAS MJ 2

            eng_rec.en_stat.group1(j + 4*(i-1)) = 
     &      hsk_recs(firas_mj_a(i)).frame(act_mj_a(i))
     &                                      .hskp_head.stat_monitor_cmd(j)

! Move B-side status words, 8 bytes from FIRAS MJ 1 and 8 from FIRAS MJ 2

            eng_rec.en_stat.group1(j + 4*(i+1)) = 
     &      hsk_recs(firas_mj_b(i)).frame(act_mj_b(i))
     &                                      .hskp_head.stat_monitor_cmd(j+4)

          end do
        end do

!  Extract the science gains for the 4 channels.


        if (XCAL_SET) then		! SMR 05/30/89 reverse the setting
          eng_rec.en_stat.group1(1) = iibset(eng_rec.en_stat.group1(1),0)
          eng_rec.en_stat.group1(9) = iibset(eng_rec.en_stat.group1(9),0)
          eng_rec.en_stat.group1(1) = iibclr(eng_rec.en_stat.group1(1),1)
          eng_rec.en_stat.group1(9) = iibclr(eng_rec.en_stat.group1(9),1)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								!
!     Move other status words:					!
!     Dwell, Micro Bus Readouts and Bolometer Bias statuses.    !
!				  				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        hsk_a_dwell = HSK_RECS(bn(1)).frame(mj(1)).hskp_head.dwell_stat(1)
        hsk_a_lvdt = HSK_RECS(bn(1)).frame(mj(1)).hskp_head.lvdt_stat(1)
       
	ENG_REC.en_stat.LVDT_stat(1) = hsk_a_lvdt
	ENG_REC.en_stat.dwell_stat(1) = ibits(hsk_a_dwell,7,1)
	ENG_REC.en_stat.grt_addr(1) = ibits(hsk_a_dwell,0,5)  
                                                            ! Dwell on a side

	ENG_REC.chan(1).dither = ibits(hsk_a_dwell,5,1)
	ENG_REC.chan(2).dither = ibits(hsk_a_dwell,6,1)
                                                            ! Dither chan 1,2

!	IF (HSK_RECS(bn(1)).frame(mj(1)).hskp_head.dwell_stat(1) .NE.
!    &	              HSK_RECS(bn(2)).frame(mj(2)).hskp_head.dwell_stat(1))
!    &                QUALFLAGS(FLG_STCHG_MJ_ST) = 1


        hsk_b_dwell = HSK_RECS(bn(1)).frame(mj(1)).hskp_head.dwell_stat(2)
        hsk_b_lvdt = HSK_RECS(bn(1)).frame(mj(1)).hskp_head.lvdt_stat(2)
       
	ENG_REC.en_stat.LVDT_stat(2) = hsk_b_lvdt
	ENG_REC.en_stat.dwell_stat(2) = ibits(hsk_b_dwell,7,1)
	ENG_REC.en_stat.grt_addr(2) = ibits(hsk_b_dwell,0,5)  
                                                            ! Dwell on a side


	ENG_REC.chan(3).dither = ibits(hsk_b_dwell,5,1)
	ENG_REC.chan(4).dither = ibits(hsk_b_dwell,6,1)
                                                            ! Dither chan 3,4

	do i = 1,4

	  ENG_REC.en_stat.micro_stat_bus(i) = 		    ! Micro bus readouts
     &                HSK_RECS(bn(1)).frame(mj(1)).hskp_head.u_proc_stat(i)
	  ENG_REC.en_stat.bol_cmd_bias(i) =		    ! Bolometer biases
     &                HSK_RECS(bn(1)).frame(mj(1)).hskp_head.bol_cmd_bias(i)


	ENDDO
	ENDIF	! Return status is 'Success'.

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_INTERPOLATE = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_INTERPOLATE = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END

