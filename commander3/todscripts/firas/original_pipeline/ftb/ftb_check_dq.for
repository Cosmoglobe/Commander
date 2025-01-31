	Program FTB_CHECK_DQ

C------------------------------------------------------------------------
C
C    PURPOSE: FTB Check data quality
C 
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Shirely READ, H. WANG
C            STX           
C            OCT, 3, 1989
C
C    INVOCATION: FTB_CHECK_DQ
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS:
C
C    SUBROUTINES CALLED:
C     FTB_LIST_SCI
C     FTB_LIST_ENGSTAT
C    COMMON VARIABLES USED:
C
C    INCLUDE FILES:
C	$SSDEF
C	CT Lib
C
C----------------------------------------------------------------------
C    CHANGE LOG:
C
CC----------------------------------------------------------------------
C----------------------------------------------------------------------
C    PDL 
C  Establish condition handler
C  Set status to success
C  Get the FIRAS Channel from the command line
C  If error then
C    Set status to error
C    write message and abort
C  Endif
C  Get File extension from the command line
C  If error then
C    Set status to error
C    write message and abort
C  Endif
C  Build the complete filenames for FDQ Science and FDQ Engineering data sets
C  Open the report file for write access
C  Write the run date, Channel and filenames in the report
C  Open the Science and Engineering files for CobeTrieve read access
C  If error then
C    Set status to error
C    write message and abort
C  Endif
C  Set total_rec to Zero, set badqual_rec to zero
C  Set badtime_rec to zero, set counts for each type of failure to zero
C  Set EOF to false
C  Read one Science record and one engineering record 
C  Increment total_rec
C  If error then
C    Set status to error
C    write message and abort
C  Endif
C  Do while not EOF and status is success
C   If (Science.Collect_time.badtime_flag .GT. Zero) Then
C     Increment badtime_rec
C     Read Science record and increment total_rec
C   ElseIf (Science.Ct_Head.time .GT. Engineering.EN_head.Sci_time(Channel)
C      .Bin_time) 
C     Then
C       Read Engineering Record
C   Elseif (Science.Ct_Head.time .LT. Engineering.EN_head.Sci_time(Channel).
C          Bin_time) Then
C         Write Warning Message to Sys$out and report
C           (include the record times)
C         Read science record
C         Increment total_rec
C    Else Check the data quality
C      If (Science.DQ_data.data_quality(110)
C           .GT. FAC_good_DQ) Then
C        Increment badqual_rec
C        Write total_rec in report
C        If any of the following Data_quality flag are set, 
C           Write this information in the report, set countes for each
C            Saturated_sample, Checksum, telemetry_ready,
C            Data_ready, Micro_Status, Status_change_across_major_frame
C        Call FTB_List_SCI
C        Call Ftb_List_ENGSTAT
C      Endif
C      Read Science and Engineering records
C      Increment total_rec
C    Endif
C   Enddo
C   Write the total_rec, badqual_rec, badtime_rec to
C    the report file
C   Write the countes for types of failure
C   Report the processing status
C   Stop      
C     
C    
C     
C       
C
CC----------------------------------------------------------------------

	implicit none

C	FTB messages

	external	ftb_normal
	external	ftb_aberr
	external	ftb_ctinit
	external        fut_normal
	external        cli$_present
	external        cli$_defaulted
        external        ct_connect_read
C	Condition Handler

	external        fut_error
	external        ftb_get_chan
	external        ftb_get_chan_file
        external        ftb_chkwarn
	external        ftb_get_file
	external        ftb_open
	external        ftb_ctopen
	external        ftb_ctread

C	Include Files

	include '($ssdef)'
        Include '(FUT_params)'
        Include '(FUT_Qualflags)'
	Include 'CT$Library:CTUser.Inc'

C       Dictionaries
        
        Dictionary 'NFS_SDF'
        Dictionary 'FDQ_ENG'

C       Records

        Record/nfs_sdf/ science
        Record/fdq_eng/ Engineering

C	Local Declarations


	integer         *4      cli$present
	integer         *4      cli$get_value
	integer         *4      LUN_RPT
        INTEGER         *4      LUN_SDF
        INTEGER         *4      LUN_ENG
        INTEGER         *4      Total_rec      
        INTEGER         *4      badqual_rec
        INTEGER         *4      badtime_rec    
        INTEGER         *4      CHK_sum_count     
        INTEGER         *4      SatU_SAM_count      
        INTEGER         *4      TELM_QUAL_COUNT      
        INTEGER         *4      DATA_READ_COUNT
        INTEGER         *4      MICR_STAT_COUNT      
        INTEGER         *4      STAU_CHANG_count
        INTEGER         *4      LEN
        INTEGER         *2      CT_return_sta(20)            
	integer		*2	chan_num
	integer		*4	iostatus
        integer         *2      zero/0/
	logical		*1	EOF 
        logical         *1      status
        CHARACTER       *19     Arch_id/'CSDR$FIRAS_ARCHIVE:'/
        CHARACTER       *16     Rpt_file
	character	*24     EXT
	character	*64     SDF_FILE
	integer  	*2      IRUN_TiME(2)
        integer         *2      nrec
        integer         *4      cat_entry/0/
        logical         *1      write_ifg/.false./       
        Character       *14     RUN_time
        CHARACTER       *29     IN_FILE
        Character       *64     ENG_file
        character       *1      blk(64)/64*' '/
        character       *1      blk1(64)/64*' '/
	character       *2      chan
        logical         *2      Time_GT
        logical         *2      Time_LT
        Character       *80     line/'  FIRAS SCIENCE  DATA LISTING '/

        equivalence (blk(1),sdf_file(1:1))   
        equivalence (blk1(1),ENG_file(1:1))   
C	Functions
C----------------------------------------------------------------------
C             BODY
C
C----------------------------------------------------------------------
C
C
C     Establish condition handler.
C
	call lib$establish ( fut_error )
C 
C     Set the status to success
C
         status = .true.
C 
c Get the channel
c
C
C
C
	If ( cli$present('channel').eq.%loc(cli$_present)) Then
	   status = CLI$Get_Value('Channel',Chan,Len)
	   If ( chan .EQ.'RH' .or. chan .EQ. 'RL' .or. chan .EQ. 'LH' .or.
     *           chan .EQ. 'LL') then	  
             status = .true. 
           else         
             call lib$signal(ftb_aberr)
             call exit(ss$_abort)
	   endif
         else
           status = .false.
           call lib$signal(ftb_get_chan)
           call lib$signal(ftb_aberr)
           call exit(ss$_abort)
    
	endif

          IF (Chan .eq. 'RH') Chan_num = 1
          IF (Chan .eq. 'RL') Chan_num = 2
          IF (Chan .eq. 'LH') Chan_num = 3
          IF (Chan .eq. 'LL') Chan_num = 4
c
c  get the file extension from the command line
c
	If ( cli$present('filename').eq.%loc(cli$_present)) Then
	   status = CLI$Get_Value('filename',in_file,Len)
           If (in_file(9:10) .ne. chan) then
            status = .false.
            call lib$signal(ftb_get_chan_file)
            call lib$signal(ftb_aberr)
            call exit(ss$_abort)
          Endif 
       Else 
           status = .false.
           call lib$signal(ftb_get_File)
           call lib$signal(ftb_aberr)
           call exit(ss$_abort)
       Endif
       Ext = In_file(12:27) 
       sdf_file =ARCH_id//'FDQ_SDF_'//Chan//'.'//ext
       Eng_file =Arch_id//'FDQ_ENG'//'.'//ext
C  open the report file
       Call Lib$Get_lun(Lun_rpt)
       Rpt_file = 'FTB_CHECK_DQ.RPT'
       Open(unit=Lun_rpt,file=Rpt_File,status='NEW',Iostat=Iostatus)    
       If (Iostatus .gt. 0) then
           call lib$signal(ftb_OPen,%val(1),%val(iostatus))
           call lib$signal(ftb_aberr)
           call exit(ss$_abort)
       Endif
C    Write The RUN date , channel and filenames in the report
C   
       Write(Lun_rpt,500) 
       Call SyS$Gettim(IRUN_time)
       CALL CT_BINARY_TO_GMT(IRUN_TIME,RUN_TIME) 
       Write(lun_rpt,501) Run_time
       Write(lun_rpt,502) SDf_file
       Write(lun_rpt,503) Eng_file
C
C    Open the Sciences and Engineering Files for cobetrieve read access 
c
c   Initialize COBETRIEVE
c
	call ct_init(ct_return_sta)
	if (ct_return_sta(1) .ne. ctp_normal) then
	   call lib$signal(ftb_ctinit, %val(1), %val(ct_return_sta(1)))
	   call lib$signal(ftb_aberr)
	   call exit(ss$_abort)
	end if
C       
       Call Lib$Get_lun(Lun_SDF)
       Call Lib$Get_lun(Lun_ENG)
       Open(Unit=Lun_sdf, File=SDF_file,status='old', iostat=iostatus,
     *           useropen=CT_CONNECT_READ)
        IF (iostatus .ne. zero) then
         Call Lib$signal(FTB_CTopen,%Val(1),%val(iostatus))
         call lib$signal(ftb_aberr)
         call exit(ss$_abort)
       Endif
       
       Open(Unit=Lun_eng, File=Eng_file,status='old', iostat=iostatus,
     *      useropen=CT_CONNECT_READ)
       IF (iostatus .ne. 0) then
         Call Lib$signal(FTB_CTopen,%Val(1),%val(iostatus))
         call lib$signal(ftb_aberr)
         call exit(ss$_abort)
       Endif
C   Set the total_rec to zero

       Total_rec = 0
       Badqual_rec= 0
       badtime_rec = 0
C   set countes for each type of failure to zero
       CHK_SUM_COUNT=0
       SATU_SAM_COUNT=0
       TELM_QUAL_COUNT= 0
       DATA_read_COUNT=0
       MICR_STAT_COUNT=0
       StAU_chang_count=0
       EOF = .false.
       Call CT_READ_ARCV(,LUN_SDF, Science,CT_return_sta)
       IF (CT_return_sta(1) .ne. CTP_normal) then            
         Call Lib$signal(FTB_CTread,%Val(1),%val(Ct_return_Sta(1)))
         call lib$signal(ftb_aberr)
         call exit(ss$_abort)
       Endif
       Call CT_READ_ARCV(,LUN_eng, ENGineering,Ct_return_sta)
       IF (Ct_return_sta(1) .ne. CTP_normal) then            
         Call Lib$signal(FTB_CTread,%Val(1),%val(Ct_return_Sta(1)))
         call lib$signal(ftb_aberr)
         call exit(ss$_abort)
       Endif
C    
       Total_rec = Total_rec + 1
       EOF = .false.
       Status = .true.
       DO While (.not. EOF .and. Status)
        If (Science.Collect_time.badtime_flag .Gt. 0) then
           badtime_rec= badtime_rec + 1
         Call CT_READ_ARCV(,LUN_SDF, Science,CT_return_sta)
         IF (CT_return_sta(1) .ne. CTP_normal) then            
           IF (Ct_return_sta(1) .EQ. CTP_endoffile) then
              EOF = .true. 
           Else        
            Call Lib$signal(FTB_CTread,%Val(1),%val(Ct_return_Sta(1)))
            call lib$signal(ftb_aberr)
            call exit(ss$_abort)
           Endif
         Endif 
         Total_rec = total_rec + 1
        ElseIf (TIME_GT(science.ct_head.time , Engineering.En_head.sci_time
     *          (Chan_num).Bin_time))
     *    Then
           Call CT_READ_ARCV(,LUN_eng, ENGineering,Ct_return_sta)
           IF (Ct_return_sta(1) .ne. CTP_normal) then
             IF (Ct_return_sta(1) .EQ. CTP_endoffile) then
               EOF = .true. 
             Else        
              Call Lib$signal(FTB_CTread,%Val(1),%val(Ct_return_Sta(1)))
              call lib$signal(ftb_aberr)
              call exit(ss$_abort)
            endif   
           Endif
        ElseIf (TIME_LT(Science.Ct_head.time, Engineering.EN_Head.sci_time
     *       (chan_num).Bin_time)) Then
           write(lun_rpt,504)science.ct_head.gmt, total_rec 
             Write(lun_rpt,*) ' '
            Call Lib$signal(FTB_ChkWarn)
           Call CT_READ_ARCV(,LUN_SDF, Science,CT_return_sta)
           IF (CT_return_sta(1) .ne. CTP_normal) then            
             IF (Ct_return_sta(1) .EQ. CTP_endoffile) then
               EOF = .true. 
             Else        
              Call Lib$signal(FTB_CTread,%Val(1),%val(Ct_return_Sta(1)))
              call lib$signal(ftb_aberr)
              call exit(ss$_abort)
             Endif
           Endif
            Total_rec = Total_rec + 1
        Else
          If (Science.Dq_data.data_quality(Flg_Summary) .Gt. Fac_good_DQ ) Then
           badqual_rec = badqual_rec + 1
           Write(lun_rpt,505) total_rec
           If (Science.dq_data.data_quality(Flg_Badsci) .ne. 0) then
              DATA_read_COUNT= DATA_read_count + 1
              Write(lun_rpt,506)
           ENDIF
           If (Science.dq_data.data_quality(Flg_Badhkp) .ne. 0) then
              Telm_Qual_COUNT= Telm_qual_count + 1
              Write(lun_rpt,507)
           ENDIF
           If (Science.dq_data.data_quality(Flg_Saturate) .ne. 0) then
              satu_sam_COUNT= satu_sam_count + 1
              Write(lun_rpt,508)
           ENDIF
              If (Science.dq_data.data_quality(Flg_Cksm_Err_St) .ne. 0) then
              CHK_SUM_COUNT= CHK_SUM_count + 1
              Write(lun_rpt,509)
           ENDIF
           If (Science.dq_data.data_quality(Flg_Stchg_mj_St) .ne. 0) then
              STAU_chang_COUNT= Stau_chang_count + 1
              Write(lun_rpt,510)
           ENDIF
           If (Science.dq_data.data_quality(Flg_Micro) .ne. 0) then
              Micr_stat_count= Micr_stat_count + 1
              Write(lun_rpt,511)
           ENDIF
             nrec=total_rec
             Write(lun_rpt,622)
             Call FTB_list_sci(LUN_rpt, science,nrec,line,cat_entry,write_ifg)
             Write(lun_rpt, 623) 
             CALL FTB_list_ENGStat(LUN_rpt,ENGINEERING,engineering.EN_stat)
           WRITE(LUN_RPT, 621) Total_rec
         Endif
         
         Call CT_READ_ARCV(,LUN_SDF, Science,CT_return_sta)
         IF (CT_return_sta(1) .ne. CTP_normal) then            
           IF (Ct_return_sta(1) .EQ. CTP_endoffile) then
              EOF = .true. 
           Else        
            Call Lib$signal(FTB_CTread,%Val(1),%val(Ct_return_Sta(1)))
            call lib$signal(ftb_aberr)
            call exit(ss$_abort)
           Endif
         Endif 
        
        
         Call CT_READ_ARCV(,LUN_eng, ENgineering,Ct_return_sta)
         IF (Ct_return_sta(1) .ne. CTP_normal) then            
           IF (Ct_return_sta(1) .EQ. CTP_endoffile) then
              EOF = .true. 
           Else        
            Call Lib$signal(FTB_CTread,%Val(1),%val(Ct_return_Sta(1)))
            call lib$signal(ftb_aberr)
            call exit(ss$_abort)
          Endif
        ENDIF
         TOtal_rec= total_rec + 1
        ENDIF
       ENDDO
        write(lun_rpt,620) 
        Write(lun_rpt,600)Total_rec,Badqual_rec, badtime_rec
        Write(lun_rpt,601)chk_sum_count,satu_sam_count,telm_qual_count,
     *        data_read_count, Micr_stat_count, Stau_chang_count        
500     format(20x,' FIRAS DATA QUALITY CHECKING REPORT'///)
501     format(20x,' Run Time : ', a14//)
502     FORmat(20x, ' Input Science file name : ', a64//)
503     format(20x, ' Input Engineering file name: ', a64/////)
504     format(10x, '--- WARNING :'/
     *        15x,' Science record header time is less than',
     *            ' Engineering record header time , Time = ', a14
     *            , ' record # = ', i5 )
505     format(1x, '--------  BAD QUALITY ERROR CHECKING  at record # = ',
     *          i5,' ---------------------------------------------------'/)
506     format(5x,' **   Missing data ready bits in science data flag is set:', 
     *             ' Any bit is off ')
507     Format(5x, ' **   Missing minor frames in housekeeping record')
511     Format(5x, ' **   Microprocessor Status word flag is set ')
508     format(5x, ' **   Saturated sample count flag is set: Too high ')
509     Format(5x, ' **   Checksum errors flag is set : Does not agree with ',
     *          'checksum in science header ')
510     Format(5x, ' **   Status changes across Major frame boundary flag is ',
     *              'set ')
600     Format(///20x,  ' The total of records : ', i5// 
     *         20x, ' The total of bad quality records : ', I5// 20x, 
     *         ' The total of bad time records : ', I5)
601     Format(//20x, ' The number of failures of checksum errors : ', i4/
     *         20x, ' The number of failures of Saturated sample count : ',
     *         i4/ 20x, ' The number of failures of Telemetry quality : ',
     *            i4/ 20x, ' The number of failures of Missing data ready :',
     *             i4/ 20x, ' The number of failures of micro status : ',
     *          i4/ 20x, ' The number of failures of status change across ',
     *         'major frame : ', i4) 
620    Format(//15x, '     ********    SUMMARY REPORT    ********************') 
621    format(/1x, '---------- END of RECORD ', i5, '-----------------------',
     *        ' ------------------------------------------------------ '///// )
622    format(/2x, ' ------  SCIENCE RECORD --------------------------------')
623    Format(//2x, ' ------  ENGINEERING RECORD ---------------------------'/) 
        call lib$signal(FTB_normal)  
        CAll EXIT( SS$_normal)
        stop
	end
