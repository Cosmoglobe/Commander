	Program FXT_Log_Extrema
C
C-----------------------------------------------------------------------------
C
C 	Program Name:
C 	  FXT_Log_Extrema
C 
C 	Program Description:
C 	  This program writes the FIRAS Engineering Extrema records to the
C         COBETRIEVE archive. The latest extrema record and FDQ_Eng records
C	  for a selected segment or time range are read from the archive.
C	  If the Init option is selected, a new FXT Extrema record is 
C	  initialized with values which will trigger a new extrema for every
C	  engineering analog for the run. This option acts as a reset tool.
C	  The engineering field values from each FDQ_Eng record are compared
C	  to the minimum and maximum values on the current FXT Extrema record.
C	  If either the current mimimum or maximum is exceeded for the selected
C	  filter number of records, the new extrema value along with the data 
C	  time is recorded the extrema record buffer. If any extrema has been 
C	  exceeded at the end of the user specified time period, the new FXT 
C	  Extrema record is written to the archive. The program will also 
C	  generate a command file for engineering plots if the user requests 
C	  this option. A minimum and maximum time offset may be selected for
C	  the start and stop times of the plots.
C
C	Author:
C	  Shirley M. Read, STX, September 1988.
C 
C       Change Log:
C--------------------------------------------------------------------------
CH      VERSION 4.4.1 SER 3306 Q. CHUNG STX 08/09/89
CH                    PROVIDE VERSION NUMBER TO TRACK SOFTWARE UPDATE.
C        SER 6413, Version 5.8, Modify FXT to initilize if no FXT_ENG_XTRM
C                  record is in the archive.
C                  H. Wang, STX, Mar. 8, 1990
C
C        SER 6602, Version 5.9, Increase the array size of catalog_rec(*)
C                   from 20 to 500
C                  H. Wang, STX, Mar. 30, 1990
C
C        SER 6727, Version 6.1, FXT must support talaris 1590t printer
C                  H. Wang, STX, May. 15, 1990
C	 SER 4171, Standardize report file name.
C		   Larry P. Rosen, STX, 29 August 1990.
C---------------------------------------------------------------------------
C 	Calling Sequence:
C 	  This is a main program.
C 
C 	Input/Output Parameters:
C 	  None
C 
C 	Input Files:
C 	  FIRAS Engineering Data Archive File : FDQ_Eng
C         FIRAS Extrema Data Archive File : FXT_Eng_Xtrm
C       
C 	Output Files:
C         FIRAS Extrema Data Archive File : FXT_Eng_Xtrm
C	  Command file to run Engplots : FXT_Engplot.Com
C 	  Report file for run : FXT_Report.Lis_yydddhhmm 
C 
C 	Include Files Used:
C 	  CT$Library:CTUser.Inc (COBETRIEVE return status definitions)
C 
C 	Subroutines and Functions Called:
C 	  FXT_Parse_Command
C	  FXT_Cat_Info
C	  FXT_Init_Report
C	  FXT_Get_Minmax
C 	  FXT_Open_Archive
C	  FXT_Get_Flags
C 	  FXT_Compare_Eng
C 	  FXT_Trigger_Plot
C	  FXT_Write_Extrema
C	  FXT_Write_Engplot
C	  FXT_Close_Xtrm
C	  CT_Close_Arcv
C	  CT_Read_Arcv
C	  Lib$Establish
C	  Lib$Signal
C 
C 	Error Handling:
C	  Establish the condition handler. Signal all errors via Lib$Signal.
C 	  Error messages to SYS$OUTPUT (terminal if run interactively, log
C 	  file if run batch) 
C 
C 	Method Used: PDL for FXT
C 	  
C	  Establish the condition handler.
C
C 	  Call FXT_Parse_Command to parse user options from the command line.
C
C	  If a report is requested,
C	       Call FXT_Init_Report to open the FXT Report
C	            File and write initial information for the run.
C
C	  Call FXT_Cat_Info to get the catalog entries for all FDQ_Eng
C	       segments for the run if selected by time range and check for
C	       valid FDQ_Eng file if selected by filename. The complete data 
C	       time range for the run is determined from the catalog 
C	       information if selected by filename or from the user specified
C	       start and stop if selected by time range.
C
C         Call FXT_Get_Minmax to get the latest time tagged FXT_Eng_Xtrm 
C	       record of engineering extrema and open a new file segment
C	       with write access for the current run.
C
C	  If enable/disable flags are requested, 
C	     Call CCT_Open_Config to open the time tagged FEX_Limflags file.
C
C 	  Call FXT_Open_Archive to open the selected FDQ_Eng archive for 
C            read access.
C
C	  Set Eof to False.
C
C 	  Do While ( not Eof )
C 
C 	     Call CT_Read_Arcv to read the FDQ_Eng record.
C
C            If (CT status is not end of file ) Then
C	        If enable/disable flags requested, 
C	           Call FXT_Get_Flags to get the corresponding enable flags.
C		Call FXT_Compare_Eng to compare engineering fields of
C		     the FDQ_Eng record to the current minima and maxima.
C               If (Engplots are requested and any extrema are exceeded
C		     for this engineering record) Then
C 	             Call FXT_Trigger_Plot to update the plot table.
C		Endif
C	     Else
C		Set Eof to True
C	     Endif
C	     
C	  Enddo While not Eof
C
C	  Call CT_Close_Arcv to close the FDQ_Eng file.
C
C	  If ( Any extrema are exceeded for the current run ) Then 
C
C	     Call FXT_Write_Extrema to write the new extrema record in
C		     the archive.
C            
C	     If ( Engplots are requested ) Then
C	        Call FXT_Write_Engplot to write the Engplot 
C		     command file.
C	     Endif
C	  Endif
C
C         Call FXT_Close_Xtrm to invoke COBETRIEVE to set the time tag for
C	       the FXT_Eng_Xtrm file and to close the FXT_Eng_Xtrm archive.
C
C	  If enable/disable flags were requested, call CCT_Close_Config.
C
C	  If a report was requested, close the report file.
C
C	  Signal the processing status and Exit.
C
C----------------------------------------------------------------------------- 

	Implicit	None

	External	FUT_Error	! Condition handler

C	Include Files

	Include		'CT$Library:CTUser.Inc'
	Include		'($SSDef)'
	Include		'(FXT_Msg)'
	Include         '(CCT_Get_Config)'

C 	Record Structures     

	Dictionary 'FDQ_Eng'
	Record /FDQ_Eng/ Eng_Rec

	Dictionary 'FXT_Eng_Xtrm'
	Record /FXT_Eng_Xtrm/ Xtrm_Rec

C	Functions

 	Integer*4 FXT_Parse_Command
 	Integer*4 FXT_Cat_Info
	Integer*4 FXT_Init_Report
	Integer*4 FXT_Get_Minmax
 	Integer*4 FXT_Open_Archive
	Integer*4 FXT_Get_Flags
 	Integer*4 FXT_Compare_Eng
 	Integer*4 FXT_Trigger_Plot
	Integer*4 FXT_Write_Extrema
	Integer*4 FXT_Write_Engplot
	Integer*4 FXT_Close_Xtrm
	Integer*4 CCT_Open_Config
	Integer*4 CCT_Close_Config

C	Local Variables

	Integer*4 	Status		! Status of processing
	Integer*4 	Retstat		! Return status from function call
	Integer*4 	Success / 1 /, Error / 2 /  ! Values for status
	Integer*4       Zero / 0 / 

	Integer*4       Data_Start(2), Data_Stop(2)  ! Data start and stop 
	Integer*4       Plot_Start(2), Plot_Stop(2)  ! Plot start and stop 

	Character*14    GMT_Start       ! Requested start GMT for data 
	Character*14    GMT_Stop        ! Requested stop GMT for data 
	Character*30	Time_Range      ! Time range of data
	Character*39    File_Seg	! User specified file-segment
	Character*2     Data_Type       ! Source type for data
	Character*33	Report_File	! File name for report
	Character*14	CurGMT		! Current run time GMT

	Integer*4	Lun_Eng, Lun_Xtrm  ! Unit numbers for CT archive
	Integer*4	Lun_Rpt         ! Report unit number
	Integer*4       Num_Seg         ! Number of FDQ_Eng segments found
					! in data time range
	Integer*4       Filter          ! Number of points to filter new extrema
	Logical*1       Init_Xtrm       ! Flag to initialize extrema record
	Logical*1       Flags           ! Enable flags for extrema checking
	Integer*4       Plot		! Plot flag ( 0=no plot,1=line,2=laser )
	Integer*4 	Min_Offtime     ! Hours to offset plot start
	Integer*4       Max_Offtime     ! Hours to offset plot stop
	Logical*1       Report          ! Flag to enable writing a report file
	Integer*2       CT_Stat(20)	! COBETRIEVE status
	Logical*1	Eof		! Flag for eof on FDQ_Eng file
	Logical*1       New_Xtrm        ! Flag to write new extrema record
	Logical*1	Exceed          ! Flag indicating extrema exceeded for
					! current record
	Logical*1       Xtrm_Table(118,2) ! Table of flags indicting extremes
					  ! are exceeded for the current run
	Logical*1       Plot_Table(118) ! Table of flags to trigger Engplots
	Logical*1       Cur_Table(118)  ! Table of flags indicting extrema
					! are exceeded for the current record
	Logical*1       Enable_Flags(118) ! Enable flags for extrema checking
	Integer*2	Ix, Jx          ! Indices

	Record  / Config_Status / Stat	        ! Structure for Get Config
	Integer*4 	Ndset / 1 /		! Number of datasets
	Character*27    Dset			! Names of datasets
	Data dset / 'CSDR$FIRAS_REF:FEX_LIMFLAGS' /
	Integer*4 	Size / 512 /            ! Dataset size
	Integer*4       Con_Lun   		! Logical unit numbers
	Integer*4       Ncache / 1 /		! Number of caches- not used
	Integer*4       Index    		! Initial cache pointers
	Integer*4 	Ref_Count		! Reference counter for cache
	Logical*1       New_Segment		! Flag for new segments accessed
	Character*1     Access_Mode / ' ' /     ! Access mode sequential-default
	Character*19    Dummy_space		! Dummy space pad
	Character*14    Con_Time		! Data time for report

	Integer*4  Cut_Register_Version
	Integer*4  Cut_Display_Banner
	Integer*4  num_vol/80/
	Integer*4  rstatus
	Integer*4  lun_out/6/
	Character*5 Version
	Parameter   (version='6.1')
C  
C	Set status for FXT processing to Success.

	Status = Success

C       Establish condition handler.       

	Call Lib$Establish ( FUT_ERROR )

	Rstatus = Cut_Register_Version(version)
	Rstatus = Cut_Display_Banner(lun_out,Num_vol,
	1 ' FIRAS Facility FXT_Log_Extrema')
	write(lun_out,61)
  61	format(//)

C       Get processing options from command line.

	Retstat = FXT_Parse_Command ( GMT_Start, GMT_Stop, 
	1	  File_Seg, Filter, Init_Xtrm, Plot, Min_offtime,
	2	  Max_Offtime, Flags, Report, Report_File, CurGMT )
	If ( Retstat .ne. %loc(FXT_Normal) ) then
	  Status = Error
	Endif

C	Open a report file for the run and write initial information.

	If ( (Status .eq. Success) .and. (Report) ) Then
	  Retstat = FXT_Init_Report ( GMT_Start, GMT_Stop, File_Seg,
	1	  Filter, Init_Xtrm, Plot, Min_Offtime, Max_Offtime,
	2         Flags, Lun_Rpt, Report_File, CurGMT )
	  If ( Retstat .ne. %loc(FXT_Normal) ) Then
	    Status = Error
	  Endif
	Endif

C	Get the catalog information for the FDQ_Eng records to be processed.

	If ( Status .eq. Success) Then
	  Retstat = FXT_Cat_Info ( GMT_Start, GMT_Stop, File_Seg,
	1	Report, Lun_Rpt, Time_Range, Data_Start, Data_Stop, 
	2       Num_Seg, Data_Type )
	  If ( retstat .ne. %loc(FXT_Normal) ) then
	    Status = Error
	  Endif
	Endif
	If ( Num_Seg .Eq. Zero ) Then
	  Status = Error
	Endif

C	Get the latest time tagged FXT_Eng_Xtrm or initialize a new 
C	FXT_Eng_Xtrm record and open a new segment for write access.

	If ( Status .eq. Success ) Then
	  Retstat = FXT_Get_Minmax ( Data_Start, Data_Stop, Time_Range, 
	1   Init_Xtrm, Data_Type, Report, Lun_Rpt, Xtrm_Rec, Lun_Xtrm )
	  If ( Retstat .ne. %loc(FXT_Normal) ) Then
	    Status = Error
	  Endif
	Endif

C	If enable/disable flags requested, open the time tagged 
C	configuration file.

	If ((Status .Eq. Success) .And. (Flags)) Then
	  Retstat = CCT_Open_Config ( Data_Start, Data_Stop, 
	1	    Ndset, Dset, Size, Access_Mode, Ncache, Con_Lun,
	2	    Index, Stat, Ref_Count )
	  If ( .Not. Retstat ) Then
	    Status = Error
	    Call Lib$Signal(FXT_CCTOpnConfig,%val(1),%val(Retstat)) 
	  Endif
	Endif

C	Initialize flags.

	If ( Status .eq. Success ) Then
	  Eof = .False.
	  Exceed = .False.
	  New_Xtrm = .False.
	  Do Ix = 1, 118
	    Do Jx = 1, 2
              Xtrm_Table(Ix,Jx) = .False.
	    Enddo
	    Cur_Table(Ix) = .False.
	    Plot_Table(Ix) = .False.
	    Enable_Flags(Ix) = .True.
	  Enddo
	Endif

C	Open the FDQ_ENG archive for read access.
 
	If ( Status .eq. Success ) Then
	    Retstat = FXT_Open_Archive ( File_Seg, Time_Range,
	1	Num_Seg, Report, Lun_Rpt, Lun_Eng )
	    If ( Retstat .ne. %loc(FXT_Normal) ) Then
	      Status = Error
	    Endif
	    If (Num_Seg .eq. Zero ) Then
	      Eof = .True.
            Endif
	Endif

C	Process FDQ_Eng records one by one until all are done.

	Do While ((.not. Eof) .and. (Status .eq. Success))

	   Call CT_Read_Arcv ( , Lun_Eng, Eng_Rec, CT_Stat )

	   If (CT_Stat(1) .ne. CTP_Normal) Then

	      If ( CT_Stat(1) .eq. CTP_Endoffile ) Then
	         Eof = .True.
	      Else		! Read error
		 Status = CT_Stat(1)
		 Call Lib$Signal(FXT_CTReadErr,%val(2),%val(Status),
	1	      'FDQ_Eng')
		 Status = Error
	      Endif		! Read error
	   Else		! Normal return from CT_Read_Arcv
	      If ( Flags ) Then
	   	 Retstat = FXT_Get_Flags ( Eng_Rec.CT_Head.Time,
	1	    Con_Lun, Report, Lun_Rpt, Enable_Flags )
	         If ( Retstat .Ne. %loc(FXT_Normal)) Then
		    Status = Error
	      	 Endif
	      Endif
	      Exceed = .False.
	      Do ix = 1,118
	         Cur_Table(ix) = .False.
	      Enddo
	      If (Status .Eq. Success) Then
		 Retstat = FXT_Compare_Eng ( Eng_Rec, Xtrm_Rec, Filter,
	1          Xtrm_Table, Cur_Table, Exceed, New_Xtrm, Enable_Flags )
		 If (Retstat .ne. %loc(FXT_Normal)) Then
	           Status = Error
		 Else		! Good return from FXT_Compare_Eng
	           If ((Plot .Gt. Zero) .And. (Exceed)) Then  
		      Retstat = FXT_Trigger_Plot ( Eng_Rec.CT_Head.Time, 
	1		Cur_Table, Min_Offtime, Max_Offtime, Plot,
	2		Report, Lun_Rpt, Plot_Table, Plot_Start, Plot_Stop)
		      If (Retstat .ne. %loc(FXT_Normal)) Then
	                Status = Error
	              Endif   ! Bad return from FXT_Trigger_Plot
		   Endif     ! Plot Gt zero and exceed
	  	 Endif       ! Good return from FXT_Compare_Eng
	      Endif         ! Good return from FXT_Get_Flags
	   Endif            ! Normal return from CT_Read_Arcv
	Enddo 		! While ((.not. eof).and.(status .eq. success))

C	Close the FDQ_Engineering Archive File.

	Call CT_Close_Arcv (, Lun_Eng, CT_Stat )
	If (CT_Stat(1) .ne. CTP_Normal) Then
	   Retstat = CT_Stat(1)
	   Call Lib$Signal(FXT_CTClosErr,%val(1),%val(Retstat))
	   Status = Error
	Endif

C	If extrema are exceeded, write the Xtrm_Rec to the archive.

	If ( (New_Xtrm) .and. (Status .eq. Success)) Then
	   Retstat = FXT_Write_Extrema ( Xtrm_Rec, Lun_Xtrm, Report,
	1     Lun_Rpt, Xtrm_Table )
	   If ( Retstat .ne. %loc(FXT_Normal)) Then
	      Status = Error
           Endif
        Endif 

C	If plot Gt zero and extrema are exceeded, 
C	write the Engplots Command file.

	If ( (New_Xtrm) .and. (Plot .Gt. Zero) .and. 
	1  (Status .eq. Success)) Then
	   Retstat = FXT_Write_Engplot ( Plot_Table, Plot_Start,
	1     Plot_Stop, Plot, Report, Lun_Rpt)
	   If ( Retstat .ne. %loc(FXT_Normal)) Then
	      Status = Error
           Endif
        Endif 

C	Close the FXT_Eng_Xtrm Archive File.

	Retstat = FXT_Close_Xtrm (Data_Start, Lun_Xtrm, Report, Lun_Rpt)
	If ( Retstat .ne. %loc(FXT_Normal)) Then
	   Status = Error
        Endif

!	If enable/disable flags were requested, close the configuration file.

	If ((Status .Eq. Success) .And. (Flags)) Then
	  Retstat = CCT_Close_Config ( Ndset, Con_Lun, Index )
	  If ( .Not. Retstat ) Then
	    Status = Error
	    Call Lib$Signal(FXT_CCtCloConfig,%val(1),%val(Retstat)) 
	  Endif
	Endif

C	Signal the processing status and exit.

	If (Status .eq. Success) Then
	  Call Lib$Signal(FXT_Normal)
	  Call Exit(SS$_Normal)
	Else
	  Call Lib$Signal(FXT_Aberr)
	  Call Exit(SS$_Abort)
	Endif
	End
