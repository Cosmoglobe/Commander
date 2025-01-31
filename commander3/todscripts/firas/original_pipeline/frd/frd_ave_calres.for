C-----------------------------------------------------------------------------
	Program FRD_AVE_CALRES
C-----------------------------------------------------------------------------
C 	Program Description:
C	    The purpose of this facility is to create a time ordered
C	  reference file called FEX_AV_CALRS.DAT. This file will contain 
C	  the average values of the calibrator resistors for FDQ to use.
C         The values of the calibrator resistors found in the housekeeping
C         major frames will be averaged over a desired period of time.
C
C 	Calling Sequence:
C 	  This is a main program.
C 
C       Programmer: Nilo G. Gonzales, Harte Wang/STX, April 5, 1991  
C
C       Change Log:
CH
CH
C 	Input Files:
C 	  FIRAS Raw Housekeeping Data Archive Files : NFS_HKP
C       
C 	Output Files:
C 	  FIRAS Reference Archive File: FEX_AV_CALRS.DAT
C 	  Report file for run: FRD_FACR.REP_yydddhhmm 
C 
C 	Include Files Used:
C
C      	  CT$Library:CTUser.Inc
C      	  $SSdef
C         FUT_Params
C	  FUT_Error
C	  $jpidef
C 
C 	Subroutines and Functions Called:
C
C 	  FRD_FACR_Get_Options
C	  FRD_Read_Sum_Avg
C	  FUT_GMT_CONV
C	  CT_Connect_Read
C	  CT_Connect_Read
C	  CT_Connect_Write
C	  CT_Connect_Write
C	  CT_Read_Arcv
C	  CT_Write_Arcv
C	  CT_Close_Arcv
C         Lib$getjpi
C	  Sys$Gettim
C	  Lib$Get_Lun
C	  Lib$Signal
C	  Lib$Locc
C	  Lib$Establish
C	  Cut_register_version
C	  Cut_display_banner
C	  Sys$asctim
C	  Cut_translate_archive_id
C	  Cut_translate_archive_id
C	  FUT_Error       
C 
C 	Method Used:  PDL for FRD_AVE_CALRES
C
C    Begin  
C 	  Establish the condition handler.
C
C	  Call FRD_FACR_GET_OPTIONS to parse user options from the command line.
C         If (Status .Eq. Success .And. Report) Then
C            Open a report file
C            Write user and current time to the report file
C         Endif
C         Call CT_Init (CT_Stat)
C
C  Convert GMT times to time in seconds. Compute time difference, num_recs which
C  are part of the equation for solving Averaging time period.
C    
C         Call FUT_GMT_CONV      ! Convert GMT_Start time to time unit seconds
C         Call FUT_GMT_CONV      ! Convert GMT_Stop time to time unit seconds
C         Time difference = Time_Sec_Stop - Time_Sec_Start
C         Number of records = NINT((Time difference)/Average time*24*3600))
C         Ave_Period = ((Time difference)/(Avetime*24*3600))/Number or records
C         Average_P = Ave_Period * Avetime
C         Do I = 1, Number of records - 1
C            Jstart_rec(I) = Time_sec_start 
C            Jstop_rec(I) = Time_sec_start + (Ave_Period*avetime*24*3600)
C            Time_sec_start = Jstop_rec(I) + 1
C         Enddo
C  Covert times in unit seconds to GMT
C         Do I = 1, Number of records
C            Call FUT_GMT_CONV      ! Convert GMT_Start from unit seconds to GMT
C            Call FUT_GMT_CONV      ! Convert GMT_Stop from unit seconds to GMT
C         Enddo
C  Get LUN for input and out reference files.
C         Call Lib$Get_Lun (Lun_Hkp)
C         Call Lib$Get_Lun (Lun_out_calrs)
C         If (Return status .Eq. Normal) Then
C             Open CT_connect_write for FEX_AV_Calrs dataset
C         Endif
C  Open CT_connect read of the housekeeping records.
C         Set I flag to 1
C         Do while (Status .EQ. Success .and. I .Le. Number of records)
C            Set reference fields structure of average period, previous GMT,
C            ADT times
C            Call CT_READ_Arcv to read housekeeping record from the archive
C            Set reference fields structure previous GMT and ADT times
C            Do while (More data .and. Status .Eq. Success)
C               Call FRD_READ_Sum_Avg to extract calibrator resistors values
C               Set First time flag to false
C               Call CT_READ_Arcv to read housekeeping record from the archive
C               If (CT_Status(1) .EQ. End of file) Then
C                  Set End segment to true
C                  Set reference fields structure of data stop time(2)
C                  Call FRD_READ_Sum_Avg to extract calibrator resistors values
C                  Call CT_Write_arcv to write reference file
C                  Call CT_Close_Arcv to close housekeeping logical unit
C                  Set More data flag to false
C               Endif
C            Enddo           !More data and status equal success
C            Add 1 to I flag
C            Call CT_Close_Arcv to close FEX logical unit
C            If (Status .Eq. Success) Then
C                Call EXit (SS$_Normal)
C            Else
C                Call EXit (SS$_Abort)
C            Endif
C       Enddo                !Status .EQ. Success .and. I .Le. Number of records
C       Turn off FUT_Report_Lun.
C    End
C----------------------------------------------------------------------------- 
	Implicit	None

C	Include Files

      	Include 'CT$Library:CTUser.Inc'
       	Include '($SSdef)'
	Include '(FUT_Params)'
	Include '(FUT_Error)'
	Include '($jpidef)'

C	Functions

 	Integer*4       FRD_FACR_Get_Options
	Integer*4       FRD_Read_Sum_Avg
	Integer*4       FUT_GMT_CONV
	Integer*4       CT_Connect_Read
	External        CT_Connect_Read
	Integer*4       CT_Connect_Write
	External        CT_Connect_Write
	Integer*4       CT_Read_Arcv
	Integer*4       CT_Write_Arcv
	Integer*4       CT_Close_Arcv

	Integer*4	Lib$getjpi
	Integer*4       Sys$Gettim
	Integer*4       Lib$Get_Lun
	Integer*4       Lib$Signal
	Integer*4       Lib$Locc
	Integer*4       Lib$Establish
	Integer*4	Cut_register_version
	Integer*4	Cut_display_banner
	Integer*4	Sys$asctim
	Integer*4 	Cut_translate_archive_id
	External  	Cut_translate_archive_id
	External        FUT_Error       ! Condition handler
	External        FRD_Normal       
	External        FRD_CTInit       
	External        FRD_CTOpenErr    
	External        FRD_CTReadErr    
	External        FRD_CTWritErr    
	External        FRD_CTClosErr    

C	Local Variables
	Character*14    gmt_start       ! gmt start time
	Character*14    gmt_stop        ! gmt start time
        Character*14    Prev_gmt        ! previous gmt time
	CHARACTER*6	version
	PARAMETER	(version='8.7')
	INTEGER*4	num_vol/80/
	CHARACTER*8	owner		! invoking user name
	INTEGER*4	current_time(2)	! current system adt time
	INTEGER*2	time_len	! length of time string
	CHARACTER*32	time		! current system time string
	INTEGER*4	rstatus		! return status
	CHARACTER*72	logn, flogn, tlogn ! logical name and translted
	INTEGER*4	flen, tlen	! length of translated logicals
	Integer*4 	Status		! Status of processing
	Integer*4 	Retstat		! Return status from function call
	Integer*4       Iostatus        ! Return status
	Integer*4 	Success / 1 /, Err / 2 /  ! Values for status
	Integer*4       Zero / 0 /      ! Zero value for Fortran status
	Integer*4       One  /0 /       ! One value
	Real*4          Avetime         ! Averaging time	
	Real*4          Average_P       ! Computed averaging period
	Integer*4	LUN_HKP  	! Unit numbers for NFS_HKP CT archive
	Integer*4	LUN_Out_CALRS	! Unit numbers for FEX_AV_CALRS.DAT
	Integer*4	LUN_Rpt         ! Report unit number
	Character*34    Report_File     ! Report file name
	Logical*1       Report          ! Flag to enable writing a report file
	Character*79	CMD_Line(2)	! Command line string with defaults
	Integer*2       FLag            ! Flag = 1 or 2 to convert time 
	Integer*2       I               ! Index
	Real*4          Ave_period      ! Computed average period      
	Integer*2       num_recs        ! Number of records to be processed

	Character*22    HKP_File / 'CSDR$FIRAS_RAW:NFS_HKP' /
	Character*31    CALRS_File / 'CSDR$FIRAS_REF:FEX_AV_CALRS.DAT' /
	Character*64    HKP_TFile       ! Housekeeping file

        Real*8          jstart_rec(700) ,jstop_rec(700)
        Character*14    Jstart_gmt(700) ,Jstop_gmt(700) 
        Real*8          Time_diff, time_sec_stop, time_sec_start
	Logical*1       More_data /.True./      ! Flag for more data
	Integer*2	CT_Stat(20)             ! Cobetrieve status
        Logical*1       First_time, End_segment ! Flags
        Integer*4       Bi_Time(2)              ! Binary time
        Integer*4       Cpbi_Time(2)              ! Binary time
        Integer*4       Time_gt 
        Dictionary      'NFS_HKP'
        Record          /NFS_HKP/  HKP_REC,LHKP_REC

	Dictionary      'FEX_AV_CALRS'
	Record          /FEX_AV_CALRS/ Fex_av_calrs_Rec

C	Set status for FRD processing to Success.

	Status = Success

C       Establish condition handler.       

	Call Lib$Establish ( FUT_ERROR )

C       Get processing options from command line.

	Retstat = FRD_FACR_Get_Options ( gmt_start, gmt_stop, Avetime,
	1                               Report, Report_File, Cmd_Line )
	If ( Retstat .ne. %loc(FRD_Normal) ) Then
	  Status = Err
	Endif
	If ((Status .Eq. Success) .And. Report) Then
C    Open a report file.
            Retstat = Lib$Get_Lun (LUN_Rpt)
	    If ( Retstat .ne. SS$_Normal ) Then
	        Status = Err
	    Else
                Open (Unit=LUN_Rpt, File=Report_File, Status='NEW',
	1       Access='Sequential',Form='Formatted',Iostat=Iostatus)
	        If ( Iostatus .Ne. Zero ) Then
	           Call Lib$Signal(FRD_CTOpenErr,%val(1), %val(Iostatus))
                   Status = Err
	        Else
	            Retstat = cut_register_version(version)
                    Retstat = cut_display_banner(lun_rpt,num_vol,
	1                  'FIRAS facility FRD_Ave_Calres')
	            WRITE(lun_rpt,10)
  10	            FORMAT(/)
	            Fut_Report_Lun = Lun_Rpt

C  Write user and current time to the report file.

	            Retstat = lib$getjpi (jpi$_username,,,,owner,)
	            CALL sys$gettim ( current_time )
	            Retstat = sys$asctim ( time_len, time, current_time, 0 )
	            WRITE (lun_rpt,20) owner, time
  20	            FORMAT (' Run by:   ', A, '   at  Time: ',A,/)
	            logn(1:15) = 'CSDR$FIRAS_RAW:'
	            rstatus = cut_translate_archive_id (logn, flogn,
	1                       flen,tlogn,tlen)
	            WRITE (lun_rpt,40) logn(1:15), tlogn
  40	            FORMAT (1X, 'Logical Translation for Input Archive:',
	1                       /,4X,A,' = ',A)
	            logn(1:15) = 'CSDR$FIRAS_REF:'
	            rstatus = cut_translate_archive_id(logn,flogn,
	1                      flen,tlogn,tlen)
	            WRITE (lun_rpt,50) logn(1:15), tlogn
  50                FORMAT(1X,'Logical Translation for Output Archive:',
	1                       /,4X,A,' = ',A)
	            WRITE (lun_rpt,60) CMD_Line(1)
	            WRITE (lun_rpt,61) CMD_Line(2)
  60                Format (1X,'Command Line with defaults: ',/,10X,A)
  61                Format(10x,A)
  70      	    Format (1x,'Calculated Averaging Period: ',F5.3,/)
	        Endif
	    Endif
	    CALL ct_init( ct_stat )
	    IF (ct_stat(1) .NE. ctp_normal) THEN
	       Retstat = ct_stat(1)
	       CALL lib$signal(frd_ctinit,%val(1),)
	       Status = Err
	    ENDIF
          Endif
C  Convert GMT times to time in seconds. Compute time difference, num_recs which
C  are part of the equation for solving Averaging time period.

        Flag = 2
        Retstat = FUT_GMT_CONV(time_sec_start,gmt_start,flag)
        Retstat = FUT_GMT_CONV(time_sec_stop,gmt_stop,flag)
            Time_diff = time_sec_stop - time_sec_start
            Num_recs = NINT((time_diff)/(avetime*24*3600))
            Ave_period = ((time_diff)/(avetime*24*3600))/num_recs
	            Average_P = Ave_Period * Avetime
	            WRITE (lun_rpt,70) Average_P
            Do I=1, num_recs-1
               Jstart_rec(I) = time_sec_start 
               jstop_rec(I) = time_sec_start + (ave_period*avetime*24*3600 )
               Time_sec_start = jstop_rec(I) + 1
            Enddo 
            Jstart_rec(num_recs) = Time_sec_start
            Jstop_rec(num_recs) = Time_sec_stop

C  Covert times in seconds to GMT.

            Do I = 1, num_recs
               Flag = 1
               Retstat= FUT_GMT_CONV(jstart_rec(I),jstart_gmt(I),flag)
               Retstat= FUT_GMT_CONV(jstop_rec(I),jstop_gmt(I),flag)
            Enddo

C  Get LUN for input and out reference files.

        Retstat = Lib$get_lun(Lun_hkp)
        Retstat = Lib$get_lun(Lun_out_calrs)
	If ( Retstat .ne. SS$_Normal ) then
	    Status = Err

	Else       ! Open CT_Connect_Write for FEX_Av_Calrs dataset.

           Open (Unit=LUN_Out_Calrs, File=Calrs_File, Status='NEW',
	1        Iostat=Iostatus, Useropen=CT_Connect_Write)
	   If ( Iostatus .Ne. Zero ) Then
	     Call Lib$Signal(FRD_CTOpenErr,%val(1), %val(Iostatus))
             Status = Err
	   Else
	       If (Report) Then
	           Write (Unit=Lun_Rpt,Iostat=Iostatus,FMT=200) Calrs_File
	       Endif      
	   Endif
        Endif
C  Open CT_Connect read of HKP records.      
        Hkp_Tfile = Hkp_file//'/'//jstart_gmt(1)//';'//jstop_gmt(num_recs)//';'
	Open (Unit=LUN_HKP, File=HKP_TFile, Status='OLD',
	1        Iostat=Iostatus ,Useropen=CT_Connect_Read)
	If (Iostatus .Ne. Zero) Then
	    Call Lib$Signal(FRD_CTOpenErr,%val(1), %val(Iostatus))
            Status = err
	Else
	    If (Report) Then
	      Write (Unit=Lun_Rpt,Iostat=Iostatus,FMT=200) HKP_TFile
	    Endif      
	Endif
        I = 1
        Do While ( Status .eq. success .and. I .le. num_recs)

C  Set reference fields structure of average period, previous GMT,
C  and ADT times. 
           call ct_gmt_to_binary(jstop_gmt(i),cpbi_time)
           Fex_av_calrs_rec.ave_period = ave_period*avetime
           If ( I .eq. 1) then
             Fex_av_calrs_rec.prev_data_start = '              '
             Retstat = CT_Read_Arcv (, Lun_Hkp, Hkp_Rec, CT_Stat )
   	     If (CT_Stat(1) .Ne. CTP_Normal) Then
	      Call Lib$Signal(FRD_CTReaderr)
              Status = err
             else
              prev_gmt = Hkp_rec.ct_head.gmt  
              Lhkp_rec = hkp_rec
             Endif
           else 
             Fex_av_Calrs_rec.prev_data_start= prev_gmt
	   endif
           If (status .eq. success) then
               Fex_av_calrs_rec.ct_head.gmt=Hkp_rec.ct_head.gmt
               call ct_gmt_to_binary(hkp_rec.ct_head.gmt,bi_time)
               fex_av_calrs_rec.ct_head.time(1)=bi_time(1)
               fex_av_calrs_rec.ct_head.time(2)=bi_time(2)
               More_data = .true.
               First_time = .true.
               End_segment = .false.
               prev_gmt = Hkp_rec.ct_head.gmt  
           endif 

C  Read housekeeping record from Cobetrieve archive.


C  Set reference fields structure of previous GMT, ADT times
	   Do While (More_Data .and. Status .Eq. Success)
C  Extract Calibrator resistors values from housekeeping records.
                  retstat=FRD_Read_Sum_AVG (Hkp_rec,first_time,End_segment,
	1	                         fex_av_calrs_rec)
                  First_time = .false.
C  Read housekeeping record from Cobetrieve archive.
	          retstat= CT_Read_Arcv (, Lun_Hkp, Hkp_Rec, CT_Stat )
	          If (CT_Stat(1) .Ne. CTP_Normal ) then
                     If (CT_Stat(1) .eq. CTP_Endoffile) Then
                        End_segment = .True.
                        Fex_av_calrs_rec.data_stop=Lhkp_rec.hskp_tail.gmt_mjf2
                        call ct_gmt_to_binary(Lhkp_rec.hskp_tail.gmt_mjf2,
	1	             bi_time)
                        fex_av_calrs_rec.data_stop_time(1)=bi_time(1)
                        fex_av_calrs_rec.data_stop_time(2)=bi_time(2)
                        retstat= FRD_Read_Sum_AVG (Hkp_rec, First_time,
	1	        End_segment, Fex_av_calrs_rec)
C  Write reference file into Cobetrieve archive.
         	        retstat=CT_Write_arcv (,LUN_Out_Calrs,
	1                           Fex_av_calrs_rec, Ct_stat)
                        If (Ct_Stat(1) .Ne. CTP_Normal) Then
                            Status = Err
        	            Call Lib$Signal (FRD_CTWRITERR)
                        else
	                 If (Report) Then
	                  Write (Unit=Lun_Rpt,Iostat=Iostatus,FMT=201)
	1	                 jstart_gmt(i), jstop_gmt(i)
	                 Endif      
                        Endif
               	        More_data = .False.
                     Else
                         status = err
        	         Call Lib$Signal(FRD_CTReaderr)
                     Endif
                  else
                     If (time_gt(hkp_rec.ct_head.time,cpbi_time)) then
                        More_data = .false.
                        End_segment = .True.
                        Fex_av_calrs_rec.data_stop=Lhkp_rec.hskp_tail.gmt_mjf2
                        call ct_gmt_to_binary(Lhkp_rec.hskp_tail.gmt_mjf2,
	1	             bi_time)
                        fex_av_calrs_rec.data_stop_time(1)=bi_time(1)
                        fex_av_calrs_rec.data_stop_time(2)=bi_time(2)
                        retstat= FRD_Read_Sum_AVG (Hkp_rec, First_time,
	1	        End_segment, Fex_av_calrs_rec)
C  Write reference file into Cobetrieve archive.
         	        retstat=CT_Write_arcv (,LUN_Out_Calrs,
	1                           Fex_av_calrs_rec, Ct_stat)
                        If (Ct_Stat(1) .Ne. CTP_Normal) Then
                            Status = Err
        	            Call Lib$Signal (FRD_CTWRITERR)
                        else
	                 If (Report) Then
	                  Write (Unit=Lun_Rpt,Iostat=Iostatus,FMT=201)
	1	                 jstart_gmt(i), jstop_gmt(i)
	                 Endif      
                        Endif
                      Else
                        Lhkp_rec = hkp_rec
                      Endif
                  Endif
               Enddo        ! More_data and Status .Eq. Success
           I = I + 1        ! Add 1 to I flag
        Enddo               ! Status .eq. success .and. I .le. num_recs 

C  Close Cobetrieve archive.

        Retstat = CT_CLOSE_ARCV(,LUN_OUT_CALRS,CT_STAT)
	If (Ct_Stat(1) .Ne. CTP_Normal) Then
           Status = Err
           Call Lib$Signal (FRD_CTClosErr)
        Endif
C  Close Cobetrieve archive if end of file.
        retstat=CT_CLOSE_ARCV(,LUN_HKP,CT_STAT)
        If (Ct_Stat(1) .Ne. CTP_Normal) Then
           Status = Err
           Call Lib$Signal (FRD_CTClosErr)
        Endif
 200    Format(1x,'Archive file successfully opened: ',/10x,a) 
 201	Format(1x,'Reference record written : ',a14,';', a14)

        If (status .eq. success) then
           Call Lib$Signal (FRD_Normal)
           CALL EXIT(ss$_normal)
        Else
           CALL EXIT(ss$_Abort)
        Endif
C   Turn FUT_Report_Lun to off.
	If (Report) Then
	   Close (Lun_rpt)
	   FUT_Report_Lun = 0
	Endif
       END        
