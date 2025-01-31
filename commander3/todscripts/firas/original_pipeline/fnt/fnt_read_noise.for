      integer * 4 function fnt_read_noise(chan_id,archive_type,first_IFG,
     .					  last_IFG,num,sci_data,fake_it)

c--------------------------------------------------------------------------
c	Author: W. K. Young
c		STI Inc.
c		July 1986
c
c	Input:	Channel number
c
c	Output:	number of records found
c		the data
c		fake-it flag  (0 = false;  1 = true)
c
c	Include files:
c		fut_params
c		fnt_invoc
c		ct$library:ctuser.inc
c
c	CALLing Sequence
c		status = fnt_read_noise(channel number,
c					number of records found,
c					the time series,
c					fake-it flag)
c	Gist of routine
c		IF channel 1 THEN
c		   QUERY for start and stop times
c                  PAD with zeros IF necessary
c		END IF
c		CALL CT_OPEN to open science data for reading
c		IF return status is an error THEN
c		   SET return status to an error
c		ELSE
c		   CALL CT_READ_ARCV to read first record
c		   IF return status is an error THEN
c		      SET return status to an ERROR
c		   ELSE
c		      DO UNTIL an error occurs or end of file
c		         INCREMENT record counter
c			 CALL CT_READ_ARCV to read science records
c		      END DO
c		      IF an error occured THEN
c		         SET return status to ERROR
c		      ELSE
c	                 SET return status to normal
c		      END IF
c                  END IF
c               END IF
c	        RETURN
c	        END 
c

c ********************************************************
c * Modification: 1989 01 13 (Friday the Thirteenth)     *
c * Purpose: To address SPR 2900. FNT needs to write     *
c *          several parameters to the noise spectrum    *
c *          records which are not currently written     *
c *          there.                                      *
c * Solution: Take the following steps:                  *
c *        ----- In FNT_power_density -----              *
c *           1. Assign Start_time & Stop_time, which are*
c *              already determined, to the proper       *
c *              Spec_rec variable.                      *
c *           2. Assign Num, which is already determined,*
c *              to the proper Spec_rec variable.        *
c *           3. Assign Bin_time_ave, which is already   *
c *              determined, to the proper Spec_rec vari-*
c *              able.                                   *
c *           4. Handle MTM length in a manner similar to*
c *              MTM speed. Then assign MTM_length to the*
c *              proper Spec_rec variable.               *
c *    ----- In FNT_power_density & FNT_read_noise ----- *
c *           5. Pass the commanded bolometer biases from*
c *              FNT_read_noise to FNT_power_density and *
c *              then assign them to the proper Spec_rec *
c *              variable.                               *
c * Programmer: Nick Iascone, STX/COBE                   *
c ********************************************************
c *******************************************************
c * Mod 1988 10 19                                      *
c * Purpose: To make enhancements per SPR 2649:         *
c *          NOISETEST gives incorrect spectral averages*
c *          over time ranges containing more than a    *
c *          single instrument configuration.           *
c * Solution: Perform a consistency check (compare first*
c *           value with ensuing values) of the follow- *
c *           ing variables - adds_per_group, mtm_speed,*
c *           mtm_length, gain, sci_mode, fakeit,       *
c *           bol_cmd_status and sweeps.                *
c * Programmer: Nick Iascone, STX/COBE                  *
c *******************************************************
c *******************************************************
c * Mod 1988 10 17                                      *
c * Purpose: To make corrections per SPR 2610: Noisetest*
c *          examines only the fake-it bits from side A *
c *          of the commAND status word 0. Noisetest    *
c *          should set the fake-it flag for each chan- *
c *          nel by examining the proper status bit.    *
c * Programmer: Nick Iascone, STX/COBE                  *
c *******************************************************
c *******************************************************
c * Mod 1988 10 10                                      *
c * Purpose: To incorporate function FUT_GOODSCI to     *
c *          avoid processing of questionable IFG (see  *
c *          SPR 2573).                                 *
c * Programmer: Nick Iascone, STX/COBE                  *
c *******************************************************
c *******************************************************
c * Mod: 1988 10 08                                     *
c * Purpose: To correct problem noted in SPR 2569. The  *
c *          NOISETEST software crashes when the number *
c *          of IFGs ingested is apparently 0.          *
c * Solution: Do not attempt to process data when there *
c *           is none!                                  *
c * Programmer: Nick Iascone, STX/COBE                  *
c *******************************************************
c *******************************************************
c * Mods: 1988 07 28                                    *
c *       Per SPR 2158 "Need to use GET_CONFIG for      *
c *       reference data sets"                          *
c * Programmer: Nick Iascone                            *
c *******************************************************
c--------------------------------------------------------------------------
c  Changes:
c	Trap for reaching IFG read limit -- 1987 Jan 27  Fred Shuman
c
c--------------------------------------------------------------------------
c **********************************************************************
c * Mods: 1988 07 27                                                   *
c *                                                                    *
c *       General: These items address SOR of 26-Jul-1988              *
c *                SPR 2158 - Need to use GET_CONFIG for reference     *
c *                           data sets                                *
c *       ModIFier: Nick Iascone                                       *
c **********************************************************************
c                                                                      
c         SPR 4174 
c         The implementation of the access to the new science data set
c         FDQ_SDF_xx, xx=RH, RL, LH, LL
c         August 21, 1989, Jon Smid
c
c	Version 4.4.2 SPR 5041, 5057, R. Kummerer, Nov 15, 1989.  Inappropriate
c		FIRAS logical for fetching raw data; change CSDR$FIRAS_ARCHIVE 
c		to CSDR$FIRAS_RAW.  Correct interpretation of /INPUT.
c	Version 4.4.02, SPR 5063, Nov 16, 1989, R. Kummerer, STX
c		Direct output with CSDR$FIRAS_OUT.
c
c	Version 5.9, SPR 6583,5129, Apr 12, 1990, H. Wang, STX
c		correction of input default of FNT command line.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IMPLICIT NONE

	include '(fut_params)'
	include '(fnt_invoc)'
	include '(cct_filespec_fields_record)'
	include 'ct$library:ctuser.inc'

	integer * 4	chan_id		!channel id
	integer * 4	start_chan	!first channel to be run
	integer * 4	i		!a counter
	integer * 4	num		!number of records found
	integer * 2	ct_stat(20)	!CT return status
	integer * 4	fake_it		!fake it flag
	integer * 4	fak_tim(2)	!time to find fake it bit
	integer * 4	ct_unit		!CT logical unit
	integer * 2	check_fake	!temporary variable
        integer * 4     len1 
        integer * 4     len2
        integer * 4     ios        

	record /filespec_fields/ file_spec
	character * 2	archive_type
	integer * 4	first_IFG(2)
	integer * 4	last_IFG(2)

	Integer		*4	status

	Integer		*4	start_time(2)
	Integer		*4	end_time(2)
        Common /Times/ Start_time, End_time
	Character	*14	starting_time
	Character	*14	ending_time
        Character       *3      input
        Character       *2      chan
        Character       *64     infile 

	Integer		*4	STR$UpCase
	Integer		*4	FUT_Timerange
        Integer		*4      LIB$Get_LUN 
        Integer		*4      CLI$Get_Value 
        Integer		*4	CT_Connect_Read
        External		CT_Connect_Read
	Integer		*4	CCT_Get_Filespec_Fields

	logical * 4 	first_time/.true./	!flag to query for times

	character * 30	time_span	!combined start-stop times

	EXTERNAL FNT_MAXRECCOU
        EXTERNAL FNT_RNOISE_NODATA    !per SPR 2569
        EXTERNAL FNT_RNOISE_INCONSIST !per SPR 2649

c ***********************************************
c * Items related to mod 1988 10 10 (SPR 2573)  *
c *                                             *
      Logical * 4 FUT_goodsci
      Logical * 4 alltests
      Data alltests /.TRUE./
      Logical * 1 tests(2)
      Logical * 1 flags(2)
c *                                             *
c ***********************************************
c ***********************************************
c * Items related to mod 1988 10 17 (SPR 2610)  *
c *                                             *
      Logical * 4 MJF1, MJF2, STILL_need
      Integer * 4 Try_num
c *                                             *
c ***********************************************
c ***********************************************
c * Items related to mod 1988 10 19: SPR 2649   *
c *                                             *
c * Miscellaneous flags, etc.                   *
c *                                             *
      Logical * 4 OK_fine, FIRST_val
c *                                             *
c * The following arrays hold the first oc-     *
c * curence of a value in element 1 and the last*
c * occurence in element 2.                     *
c *                                             *
c * Science record items                        *
c *                                             *
      Integer * 2 Gains,       MTM_speeds
      Integer * 2 MTM_lengths, SC_head9s
      Integer * 2 SC_head11s
      Byte        SC_head1as
c *                                             *
c * Housekeeping record items                   *
c *                                             *
      Byte Bol_cmd_biases, Bol_cmd_biases_first
      Integer * 4 Fake_it_first
c *                                             *
c ***********************************************
c * Item related to Mod 1989 01 13/SPR2900      *
c * For the purpose of communicating the value  *
c * of Bol_cmd_biases to FNT_power_density      *
c *                                             *
      COMMON /BOL_BLOCK/ Bol_cmd_biases
c *                                             *
c ***********************************************
c
	dictionary 'NFS_SDF'
	record /NFS_SDF/ sci_data(fac_max_num)	!science records
	record /NFS_SDF/ dum_rec	!dummy records
	dictionary 'NFS_HKP'
	record /NFS_HKP/ hskp_rec

C Initialize.

	fake_it = 0		!set fake it mode off until later

	fnt_read_noise = fac_normal

	IF (first_time) THEN

	   IF (ftc_interactive .EQ. fac_present .AND.
     .	        (ftc_jstart_select .EQ. fac_not_present .AND.
     .		 ftc_jstop_select  .EQ. fac_not_present)) THEN

	       status = FUT_Timerange ( starting_time,
     .				   start_time, ending_time, end_time )

            ELSE

	       start_time(1) = ftc_jstart(1)
	       start_time(2) = ftc_jstart(2)
	       end_time(1) = ftc_jstop(1)
	       end_time(2) = ftc_jstop(2)

	    END IF

	    CALL CT_Binary_To_GMT(start_time,starting_time)
            CALL CT_Binary_To_GMT(end_time,ending_time)

	    time_span = starting_time // ';' // ending_time // ';'
	    first_time = .false.

            CALL ct_init(ct_stat)

            If (ftc_interactive .eq. fac_present) then
               Type 51
51             format (' Enter Data Type (Raw, FPP, or FDQ)/[FDQ]: ',$ )
               Accept 55, input
55             format ( A )      
	       status = STR$UpCase(input,input)
	    else
               input = ftc_input
            end if

	END IF

C Read the raw science data.

	Type 5, fac_channel_ids(chan_id)
5	format (/,' Fetching data for channel ', a2)

        status = lib$get_lun ( ct_unit )

        chan = fac_channel_ids(chan_id)  
        
        If ( input(1:1) .eq. 'R' ) then
             infile = 'CSDR$FIRAS_RAW:NFS_SDF_' // chan //'/'// Time_span
        else if (input(1:2) .eq. 'FP' ) then
             infile = 'CSDR$FIRAS_RAW:FPP_SDF_' // chan //'/'// Time_span
        else
             infile = 'CSDR$FIRAS_RAW:FDQ_SDF_' // chan //'/'// Time_span
        end if

        open(unit = ct_unit,
     .       file = infile,
     .       status = 'old',
     .       iostat = ios,
     .       useropen = ct_connect_read )

        If ( ios .ne. 0 ) then
	   fnt_read_noise = ct_stat(1)
           type *,'Failed to open COBETRIEVE. CT error = ',ct_stat(1)
	ELSE
	   CALL ct_read_arcv(,ct_unit,dum_rec,ct_stat)
           First_time = .TRUE.     !Mod 1988 10 19: SPR 2649

	   num = 0
           OK_fine = .TRUE.

	   DO WHILE(ct_stat(1) .EQ. ctp_normal .AND. num .LT. fac_max_num
     .              .AND. OK_fine)
              IF (FUT_goodsci(dum_rec,alltests,tests,flags)) THEN !SPR 2573
                 num = num + 1
                 sci_data(num) = dum_rec
                 IF (First_time) THEN !SPR 2649: First occurence of values
                    Gains       = dum_rec.SCI_HEAD.GAIN
                    MTM_speeds  = dum_rec.SCI_HEAD.MTM_SPEED
                    MTM_lengths = dum_rec.SCI_HEAD.MTM_LENGTH
                    SC_head1as  = dum_rec.SCI_HEAD.SC_HEAD1A !Science mode
                    SC_HEAD9S   = dum_rec.SCI_HEAD.SC_HEAD9  !Adds per group
                    SC_HEAD11S  = dum_rec.SCI_HEAD.SC_HEAD11 !Sweeps
		    first_IFG(1) = dum_rec.ct_head.time(1)
		    first_IFG(2) = dum_rec.ct_head.time(2)
		    last_IFG(1) = dum_rec.ct_head.time(1)
		    last_IFG(2) = dum_rec.ct_head.time(2)
		    status = cct_get_filespec_fields(ct_unit,file_spec)
		    archive_type = file_spec.filename_extension(1:2)
                    First_time = .FALSE.
                 ELSE !Not the first occurence of values
                    IF (Gains .NE. dum_rec.SCI_HEAD.GAIN)
     .                 OK_fine= .FALSE.
                    IF (MTM_speeds  .NE. dum_rec.SCI_HEAD.MTM_SPEED)
     .                 OK_fine = .FALSE.
                    IF (MTM_lengths .NE. dum_rec.SCI_HEAD.MTM_LENGTH)
     .                 OK_fine = .FALSE.
                    IF (SC_head1as  .NE. dum_rec.SCI_HEAD.SC_HEAD1A)
     .                 OK_fine = .FALSE.      !Science mode
                    IF (SC_HEAD9S   .NE. dum_rec.SCI_HEAD.SC_HEAD9 )
     .                 OK_fine = .FALSE.      !Adds per group
                    IF (SC_HEAD11S  .NE. dum_rec.SCI_HEAD.SC_HEAD11)
     .                 OK_fine = .FALSE.      !Sweeps
		    last_IFG(1) = dum_rec.ct_head.time(1)
		    last_IFG(2) = dum_rec.ct_head.time(2)
                 END IF

              ELSE !SPR 2573
                 WRITE(*,'(2a)')' **** Bad IFG at: ',dum_rec.ct_head.gmt
                 IF(.NOT.flags(1))
     .           WRITE(*,*)' Failed Data Ready Flag test'
                 IF(.NOT.flags(2))
     .           WRITE(*,*)' Failed Data Checksum test'
              END IF
              CALL ct_read_arcv(,ct_unit,dum_rec,ct_stat)
	   END DO
c
c ********************************************************
c * Change required by SPR 2569                          *
c *                                                      *
           IF (num .EQ. 0) THEN
              FNT_READ_NOISE = 10 !No data: SPR 2569
              CALL LIB$SIGNAL(FNT_RNOISE_NODATA) !SPR 2569
              RETURN
           END IF
c *                                                      *
c ********************************************************
c
c *****************************************************************
c * Mod 1988 10 19 (SPR 2649) items. See if the first value       *
c * has changed.                                                  *
c *                                                               *
        IF (.NOT. OK_fine) THEN ! Data inconsistent. Bye.
           FNT_READ_NOISE = 15
           CALL LIB$signal(FNT_RNOISE_INCONSIST)
           RETURN
        END IF
c *                                                               *
c *****************************************************************

	   IF (ct_stat(1) .EQ. ctp_endoffile .OR.
     .         ct_stat(1) .EQ. ctp_normal) THEN

c IF max SCI records were read without encountering
c end-of-file, issue Info-level message ...

	      IF (ct_stat(1) .EQ. ctp_normal) THEN
	         CALL LIB$SIGNAL(FNT_MAXRECCOU, %Val(1), %Val(num))
	      END IF

	      fnt_read_noise = fac_normal
	      CALL ct_close_arcv(,ct_unit,ct_stat)
	   ELSE
	      type *,'Did''t read in data. CT code is ',ct_stat(1)
	      fnt_read_noise = ct_stat(1)
              RETURN
	   END IF
	END IF

      type 10, fac_channel_ids(chan_id)
10    format (' Processing data for channel ', a)

c Find a fake-it byte from housekeeping file ...

      infile = 'CSDR$FIRAS_RAW:NFS_HKP' // '/' // Time_span
      Open(unit = ct_unit,
     .     file = infile,
     .     status = 'old',
     .     iostat = ios,
     .     useropen = ct_connect_read )

      IF (ios .NE. 0) THEN
         fnt_read_noise = ios
         TYPE *,'While attempting HSKP_REC OPEN,'
         TYPE *,'Failed to open COBETRIEVE. CT error = ',
     .          ios
         fake_it = 1 !Set the fake it because no house keeping
         RETURN
      ELSE
         CALL ct_read_arcv(,ct_unit,hskp_rec,ct_stat)
         Try_num = 0
         OK_fine = .TRUE.
         FIRST_val = .TRUE.
         DO WHILE (ct_stat(1) .EQ. ctp_normal .AND.
     .             Try_num .lt. fac_max_num .AND. OK_fine)
             Try_num = Try_num + 1
c
c ********************************************
c * Channel 1 case                           *
c *                                          *
            IF (Chan_id .EQ. 1) THEN
               MJF1 = .FALSE.
               MJF2 = .FALSE.
               check_fake =
     .         hskp_rec.frame(1).hskp_head.stat_monitor_cmd(1)
               IF (.NOT. btest(check_fake,14)) MJF1 = .TRUE.
               check_fake =
     .         hskp_rec.frame(2).hskp_head.stat_monitor_cmd(1)
               IF (.NOT. btest(check_fake,14)) MJF2 = .TRUE.
               IF ((.NOT. MJF1) .AND. (.NOT. MJF2)) THEN
               ELSE
                  IF (MJF1) THEN
                     check_fake =
     .               hskp_rec.frame(1).hskp_head.stat_monitor_cmd(1)
                     Bol_cmd_biases =
     .               hskp_rec.frame(1).hskp_head.bol_cmd_bias(Chan_id)
                  END IF
                  IF (MJF2) THEN
                     check_fake =
     .               hskp_rec.frame(2).hskp_head.stat_monitor_cmd(1)
                     Bol_cmd_biases =
     .               hskp_rec.frame(2).hskp_head.bol_cmd_bias(Chan_id)
                  END IF

                  fake_it = 0
                  IF (btest(check_fake,9)) fake_it = 1 !Channel 1 bit
                  IF (FIRST_val) Fake_it_first = Fake_it
                  IF (FIRST_val) Bol_cmd_biases_first = Bol_cmd_biases
                  FIRST_val = .FALSE.
                  IF (Fake_it .NE. Fake_it_first) OK_fine = .FALSE.
                  IF (Bol_cmd_biases .NE. Bol_cmd_biases_first)
     .            OK_fine = .FALSE.
               END IF
	    END IF
c *                                          *
c ********************************************
c
c ********************************************
c * Channel 2 case                           *
c *                                          *
            IF (Chan_id .EQ. 2) THEN
               MJF1 = .FALSE.
               MJF2 = .FALSE.
               check_fake =
     .         hskp_rec.frame(1).hskp_head.stat_monitor_cmd(1)
               IF (.NOT. btest(check_fake,14)) MJF1 = .TRUE.
               check_fake =
     .         hskp_rec.frame(2).hskp_head.stat_monitor_cmd(1)
               IF (.NOT. btest(check_fake,14)) MJF2 = .TRUE.
               IF ((.NOT. MJF1) .AND. (.NOT. MJF2)) THEN
               ELSE
                  IF (MJF1) THEN
                     check_fake =
     .               hskp_rec.frame(1).hskp_head.stat_monitor_cmd(1)
                     Bol_cmd_biases =
     .               hskp_rec.frame(1).hskp_head.bol_cmd_bias(Chan_id)
                  END IF
                  IF (MJF2) THEN
                     check_fake =
     .               hskp_rec.frame(2).hskp_head.stat_monitor_cmd(1)
                     Bol_cmd_biases =
     .               hskp_rec.frame(2).hskp_head.bol_cmd_bias(Chan_id)
                  END IF

                  fake_it = 0
                  IF (btest(check_fake,8)) fake_it = 1 !Channel 2 bit
                  IF (FIRST_val) Fake_it_first = Fake_it
                  IF (FIRST_val) Bol_cmd_biases_first = Bol_cmd_biases
                  FIRST_val = .FALSE.
                  IF (Fake_it .NE. Fake_it_first) OK_fine = .FALSE.
                  IF (Bol_cmd_biases .NE. Bol_cmd_biases_first)
     .            OK_fine = .FALSE.
               END IF
	    END IF
c *                                          *
c ********************************************
c
c ********************************************
c * Channel 3 case                           *
c *                                          *
            IF (Chan_id .EQ. 3) THEN
               MJF1 = .FALSE.
               MJF2 = .FALSE.
               check_fake =
     .         hskp_rec.frame(1).hskp_head.stat_monitor_cmd(5)
               IF (.NOT. btest(check_fake,14)) MJF1 = .TRUE.
               check_fake =
     .         hskp_rec.frame(2).hskp_head.stat_monitor_cmd(5)
               IF (.NOT. btest(check_fake,14)) MJF2 = .TRUE.
               IF ((.NOT. MJF1) .AND. (.NOT. MJF2)) THEN
               ELSE
                  IF (MJF1) THEN
                     check_fake =
     .               hskp_rec.frame(1).hskp_head.stat_monitor_cmd(5)
                     Bol_cmd_biases =
     .               hskp_rec.frame(1).hskp_head.bol_cmd_bias(Chan_id)
                  END IF
                  IF (MJF2) THEN
                     check_fake =
     .               hskp_rec.frame(2).hskp_head.stat_monitor_cmd(5)
                     Bol_cmd_biases =
     .               hskp_rec.frame(2).hskp_head.bol_cmd_bias(Chan_id)
                  END IF

                  fake_it = 0
                  IF (btest(check_fake,9)) fake_it = 1 !Channel 3 bit
                  IF (FIRST_val) Fake_it_first = Fake_it
                  IF (FIRST_val) Bol_cmd_biases_first = Bol_cmd_biases
                  FIRST_val = .FALSE.
                  IF (Fake_it .NE. Fake_it_first) OK_fine = .FALSE.
                  IF (Bol_cmd_biases .NE. Bol_cmd_biases_first)
     .            OK_fine = .FALSE.
               END IF
	    END IF
c *                                          *
c ********************************************
c ********************************************
c * Channel 4 case                           *
c *                                          *
            IF (Chan_id .EQ. 4) THEN
               MJF1 = .FALSE.
               MJF2 = .FALSE.
               check_fake =
     .         hskp_rec.frame(1).hskp_head.stat_monitor_cmd(5)
               IF (.NOT. btest(check_fake,14)) MJF1 = .TRUE.
               check_fake =
     .         hskp_rec.frame(2).hskp_head.stat_monitor_cmd(5)
               IF (.NOT. btest(check_fake,14)) MJF2 = .TRUE.
               IF ((.NOT. MJF1) .AND. (.NOT. MJF2)) THEN
               ELSE
                  IF (MJF1) THEN
                     check_fake =
     .               hskp_rec.frame(1).hskp_head.stat_monitor_cmd(5)
                     Bol_cmd_biases =
     .               hskp_rec.frame(1).hskp_head.bol_cmd_bias(Chan_id)
                  END IF
                  IF (MJF2) THEN
                     check_fake =
     .               hskp_rec.frame(2).hskp_head.stat_monitor_cmd(5)
                     Bol_cmd_biases =
     .               hskp_rec.frame(2).hskp_head.bol_cmd_bias(Chan_id)
                  END IF

                  fake_it = 0
                  IF (btest(check_fake,8)) fake_it = 1 !Channel 4 bit
                  IF (FIRST_val) Fake_it_first = Fake_it
                  IF (FIRST_val) Bol_cmd_biases_first = Bol_cmd_biases
                  FIRST_val = .FALSE.
                  IF (Fake_it .NE. Fake_it_first) OK_fine = .FALSE.
                  IF (Bol_cmd_biases .NE. Bol_cmd_biases_first)
     .            OK_fine = .FALSE.
               END IF
	    END IF
c *                                          *
c ********************************************

            IF (OK_fine) CALL ct_read_arcv(,ct_unit,hskp_rec,ct_stat)
         END DO

c *****************************************************************
c * Mod 1988 10 19 (SPR 2649) items. See if the first value       *
c * has changed.                                                  *
c *                                                               *
        IF (.NOT. OK_fine) THEN ! Data inconsistent. Bye.
           FNT_READ_NOISE = 15
           CALL LIB$signal(FNT_RNOISE_INCONSIST)
           RETURN
        END IF
c *                                                               *
c *****************************************************************
      END IF

      CALL ct_close_arcv(,ct_unit,ct_stat)

      RETURN
      END
