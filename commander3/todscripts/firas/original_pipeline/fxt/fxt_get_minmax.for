C-------------------------------------------------------------------------------

	Integer*4 Function FXT_Get_Minmax ( Data_Start, Data_Stop,
	1	  Time_Range, Init_Xtrm, Data_Type, Report, Lun_Rpt, 
	2	  Xtrm_Rec, Lun_Xtrm )

C-------------------------------------------------------------------------------
C
C	Purpose: To get the FXT_Eng_Xtrm record for the run. If the Xtrm_Init
C	         flag is set, a new initialized record will be prepared for
C	         the current run. If the flag is not set, the data time tag 
C	         will be used to get an FXT_Eng_Xtrm record from the archive.
C
C	Author: Shirley M. Read
C		STX, December, 1988
C
C	Invocation: Status = FXT_Get_Minmax ( Data_Start, Data_Stop,
C                   Time_Range, Init_Xtrm, Data_Type, Report, Lun_Rpt, 
C		    Xtrm_Rec, Lun_Xtrm )
C
CH	Change Log:
CH        SER 6413, Modify FXT to initialize if no fxt_eng_xtrm record is
CH                  in the archive.
CH               H. Wang, STX, Mar. 8, 1990
Ch
CH        SPR 6602, Increase the array number of catalog_rec(*) from 20 to
CH                  500.
CH               H. Wang, STX, Mar. 30, 1990
C	  ----------------------------------------------------------------------
C
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name		  Type		Description
C 	  ----------------------------------------------------------------------
C         Data_Start(2)   I*4		Start time of data	    
C	  Data_Stop(2)    I*4		Stop time of data
C	  Time_Range      C*30		ASCII data time range
C	  Init_Xtrm       L*1           Flag to initialize a new record
C	  Data_Type       C*2           Source type for data
C	  Report          L*1		Flag to write report
C	  Lun_Rpt         I*4		Logical unit for report
C	
C	Output Parameters:
C	  Name		  Type		Description
C 	  ----------------------------------------------------------------------
C	  Xtrm_Rec	  FXT_Eng_Xtrm  Extrema archive record
C	  Lun_Xtrm 	  I*4		Logical unit for FXT_Eng_Xtrm file
C	
C	Subroutines Called:
C
C	  CCT_Open_Config 
C	  CCT_Get_Config_TOD
C	  CCT_Close_Config
C	  Lib$Get_Lun  
C         Lib$Signal
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C	  FXT_Msg.Txt
C         CCT_Get_Config
C	  CT$Library:CTUser.Inc
C	  $SSDef
C
C	Processing Method:
C
C	  Set the function status to FXT_Normal.
C	  If the Init_Xtrm flag is set,
C	    Initialize a new FXT_Eng_Xtrm record for the run.
C	  Else
C   	    Call CCT_Open_Config using the FXT_Eng_Xtrm dataset name and the
C	       data start and stop times.
C	    Call CCT_Get_Config_Tod using the data start time to get the latest
C              time tagged FXT_Eng_Xtrm record.
C	    Call CCT_Close_Config for the FXT_Eng_Xtrm dataset.
C	    If no engineering extrema record was found, then 
C 	      Signal an error and set the function status to FXT_Aberr.
C	  Endif
C	  If an engineering extrema record was found or initialized, then
C	     Get a logical unit number for opening the new FXT_Eng_Xtrm file
C	       to contain the data for the current run.
C	     Load the times in the CT Header.
C	     Build the new filename and extension using the FXT_Eng_Xtrm 
C	       datset name, the data type and the data start time.
C 	     Open the FXT_Eng_Xtrm file with write access.
C	     If a report is requested, log the filename into the report file.
C	  Endif
C  	  Return with normal or error status.
C
C--------------------------------------------------------------------------------
C
	Implicit  None

!	Passed Parameters.

	Integer*4         Data_Start(2)     ! Start time of data	    
	Integer*4         Data_Stop(2)      ! Stop time of data
	Character*30	  Time_Range        ! ASCII data time range
	Logical*1	  Init_Xtrm         ! Flag to initialize a new record
	Character*2       Data_Type         ! Source type for data
	Logical*1	  Report            ! Flag to write report
	Integer*4	  Lun_Rpt           ! Logical unit for report
	Integer*4	  Lun_Xtrm          ! Logical unit for FXT_Eng_Xtrm file

!	Include Files.

	Include		'($SSDef)'
	Include		'(FXT_Msg)'
C
	Dictionary        'FXT_Eng_Xtrm'
	Record	/FXT_Eng_Xtrm/ Xtrm_Rec	    ! Extrema archive record
C

!	Structure for CCT_Get_Config_TOD

	Include	        '(CCT_Get_Config)'

	Record /Config_Status/ Stat

!	Local Declarations.

	Integer*4 	CT_Connect_Write
	External 	CT_Connect_Write
	Real*4          Initval(2) / 1.7E+38, - 1.7E+38 /
	Integer*4 	Ndset / 1 /		! Number of datasets
	Character*26    Dset			! Name of dataset for qcat
	Data Dset / 'CSDR$FIRAS_IN:FXT_ENG_XTRM' /
	Character*27    Out_Dset		! Name of dataset for write
	Data Out_Dset / 'CSDR$FIRAS_OUT:FXT_ENG_XTRM' /
	Integer*4 	Size / 3072 /           ! Dataset size
	Integer*4       Con_Lun   		! Logical unit numbers
	Integer*4       Ncache / 1 /		! Number of caches- not used
	Integer*4       Index    		! Initial cache pointer
	Integer*4 	Ref_Count		! Reference counter for cache
        Logical*1       Reinit_Fxt
	Logical*1       New_Segment		! Flag for new segments accessed
	Character*1     Access_Mode / ' ' /     ! Access mode sequential-default
	Character*19    Dummy_space		! Dummy space pad
	Character*14    Con_Time		! Data time converted for report
	Integer*4	Iostatus                ! Return status from I/O
	Integer*4       Zero / 0 /              ! Zero value
	Character*7     CT_End / '9936523' /    ! Final CT file extension time
	Character*46    Filename                ! Complete filename
	Character*46    Blank                   ! Blank characters
	Character*1     Blanks(46) / 46 * ' ' / 
	Equivalence     ( Blank, Blanks(1) )
	Integer*2       Ix, Jx, Kx              ! Indices
!	Functions.
	Include	        '(CCT_Query_ttg_catalog_record)'
	Include		'CT$Library:CTUser.Inc'
        Record/query_ttg_catalog/query_ttg_catalog_rec
	Dictionary        'CCM_CME_CATALOG_ENTRY'
	Record	/CCM_CME_CATALOG_ENTRY/ CATALOG_REC(500)	    ! 
	Integer*4	CCT_query_ttg_catalog
	Integer*4 	Status		        ! System return status 
        Integer*2       no_catalog_records
        integer*2       returned_status_block(20)
        character*30    primary_key_value 
        Character*14    Quer_Start
        Character*14    Quer_Stop

	Integer*4 	CCT_Open_Config     
	Integer*4       CCT_Get_Config_Tod
	Integer*4	CCT_Close_Config
	Integer*4       Lib$Get_Lun     ! Get logical unit number


!	Set function status to Normal
	
	FXT_Get_MinMax = %loc(FXT_Normal)
        quer_start =  '80001000000000'
        quer_stop  ='99001000000000'

        Reinit_Fxt = .false.
        Query_ttg_catalog_rec.archive_id='CSDR$FIRAS_IN' 
        Query_ttg_catalog_rec.dataset_name='FXT_ENG_XTRM' 
        Call CT_GMT_TO_BINARY(Quer_start,query_ttg_catalog_rec.start_time)
        Call CT_GMT_TO_BINARY(Quer_stop,query_ttg_catalog_rec.stop_time)
C
        Status = CCT_QUERY_TTG_CATALOG(Query_ttg_catalog_rec,
	1             Catalog_rec,no_catalog_records)
	  If ( .Not. Status ) Then
	    FXT_Get_MinMax = %loc(FXT_Aberr)
	    Call Lib$Signal(FXT_CCTOpnConfig,%val(1),%val(Status)) 
	  Endif
        If (status .and. (no_catalog_records .eq. 0)) then
          Reinit_fxt = .true.
        Endif                                
!	If the Init_Xtrm flag is set, initialize a new FXT_Eng_Xtrm 
!       record for the run.

	If ( Init_Xtrm .or. reinit_fxt) Then

	  Do Ix = 1, 2		! Loop over minima and maxima

	    Do Jx = 1, 64	! Loop over GRTs
	      Xtrm_Rec.MinMax(Ix).Grt(Jx) = Initval(Ix)
	    Enddo
	    Do Jx = 1, 64	! Loop over GRT times
	      Do Kx = 1, 2
	        Xtrm_Rec.MinMax(Ix).Grt_Time(Kx,Jx) = Data_Start(Kx)
	      Enddo
	    Enddo

	    Do Jx = 1, 16	! Loop over temperatures and currents
	      Xtrm_Rec.MinMax(Ix).T_and_I(Jx) = Initval(Ix)
	    Enddo
	    Do Jx = 1, 16	! Loop over temperature and current times
	      Do Kx = 1, 2
	        Xtrm_Rec.MinMax(Ix).T_and_I_Time(Kx,Jx) = Data_Start(Kx)
	      Enddo
	    Enddo
	   
	    Do Jx = 1, 36	! Loop over volts and currents
	      Xtrm_Rec.MinMax(Ix).V_and_I(Jx) = Initval(Ix)
	    Enddo
	    Do Jx = 1, 36	! Loop over volt and current times
	      Do Kx = 1, 2
	        Xtrm_Rec.MinMax(Ix).V_and_I_Time(Kx,Jx) = Data_Start(Kx)
	      Enddo
	    Enddo

	    Do Jx = 1, 2	! Loop over LMACs
	      Xtrm_Rec.MinMax(Ix).LMACs(Jx) = Initval(Ix)
	    Enddo
	    Do Jx = 1, 2	! Loop over LMAC times
	      Do Kx = 1, 2
	        Xtrm_Rec.MinMax(Ix).LMAC_Time(Kx,Jx) = Data_Start(Kx)
	      Enddo
	    Enddo

	Enddo			! Loop for minima and maxima

	Else			! Not Init_Xtrm .or. not reinit

!       Open the time tagged configuration file.

	  Status = CCT_Open_Config ( Data_Start, Data_Stop, 
	1	    Ndset, Dset, Size, Access_Mode, Ncache, Con_Lun,
	2	    Index, Stat, Ref_Count )
	  If ( .Not. Status ) Then
	    FXT_Get_MinMax = %loc(FXT_Aberr)
	    Call Lib$Signal(FXT_CCTOpnConfig,%val(1),%val(Status)) 
	  Endif

!	Get the FXT_Eng_Xtrm time tagged, configuration file.

	  If ( FXT_Get_MinMax .Eq. %loc(FXT_Normal) ) Then
     	    Call CT_Binary_To_Gmt(Data_Start,Con_Time)

	    Status = CCT_Get_Config_Tod( Data_Start, Ndset,
	1     Size, Con_Lun, Index, Xtrm_Rec, New_Segment, Stat)
	    If ( .Not. Status ) Then
	      FXT_Get_MinMax = %loc(FXT_Aberr)
	      Call Lib$Signal( FXT_CCTGetConfig, %val(1), %val(Status))
	      If ( Report ) Then
     		Write(Unit=Lun_Rpt,FMT=100,Iostat=Iostatus) Con_Time,
	1	     Time_Range, Status
 100 	        Format(1x/5x,'Error on CCT_Get_Config_TOD for Extrema ',
	1	     'Record. Input Time= ',a//
	2	     10x,'Open_Config Time: ',a,' Status(Hex)= ',z8.8)
	      Endif
	    Else
	      If ( Report .And. New_Segment ) Then
     		Write(Unit=Lun_Rpt,FMT=200,Iostat=Iostatus) 
	1	  Con_Time, Dset
 200		Format(1x/5x,'Extrema Segment Accessed at Time: ',a//
	1         10x,' for Configuration Dataset: ',a )
	      Endif
	    Endif
	  Endif		! FXT_Get_MinMax Eq FXT_Normal	

!	Close the Configuration File

	  If ( FXT_Get_MinMax .Eq. %loc(FXT_Normal) ) Then
	    Status = CCT_Close_Config ( Ndset, Con_Lun, Index )
	    If ( .Not. Status ) Then
	      FXT_Get_MinMax = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_CCTCloConfig,%val(1),%val(Status)) 
	    Endif
	  Endif
	Endif		! Init_Xtrm or not

!	If an engineering extrema record was found or initialized, proceed
!	with the opening of the output FXT_Eng_Xtrm file.

	If ( FXT_Get_MinMax .Eq. %loc(FXT_Normal)) Then

!	Set the times in the CT Header for the output FXT_Eng_Xtrm record.

	  Xtrm_Rec.CT_Head.Gmt(1:14) = Time_Range(1:14)

	  Do Ix = 1, 2
	    Xtrm_Rec.CT_Head.Time(Ix) = Data_Start(Ix)
	  Enddo

!	Build the filename for the output FXT_Eng_Xtrm record and open the
!	file with COBETRIEVE write access.

	  Filename(1:46) = Blank(1:46)
	  Filename = Out_Dset // '.' // Data_Type // '_' // 
	1	     Time_Range(1:7) // '_' // CT_End

!	Get a logical unit number for opening the new FXT_Eng_Xtrm file
!	to contain the data for the current run.

	  Status = Lib$Get_Lun (Lun_Xtrm)
	  If ( Status .Ne. SS$_Normal ) Then
	    FXT_Get_MinMax = %loc(FXT_Aberr)
	    Call Lib$Signal(FXT_GetLunErr, %val(1), %val(Status))
	  Endif
	Endif		! Function status is normal

	If ( FXT_Get_MinMax .Eq. %loc(FXT_Normal)) Then

!	Open the file for write access.

	  Open ( Unit=Lun_Xtrm, File=Filename, Status='NEW', 
	1      Iostat=Iostatus, Useropen=CT_Connect_Write )

	  If (Iostatus .Ne. Zero ) Then
	    FXT_Get_Minmax = %loc(FXT_Aberr)
	    Call Lib$Signal (FXT_CTOpenErr, %val(1), %val(Iostatus))
	    If (Report) Then
	       Write(Unit=Lun_Rpt,FMT=300,Iostat=Status) 
	1      Filename, Iostatus
	       If (Status .Ne. Zero ) Then
	 	 FXT_Get_MinMax = %loc(FXT_Aberr)
	         Call Lib$Signal(FXT_WritErr, %val(1), %val(Status))
	       Endif
	    Endif	! Report
	  Elseif (Report) Then
	    Write(Unit=Lun_Rpt,FMT=400,Iostat=Status) Filename    
	    If (Status .Ne. Zero ) Then
	      FXT_Get_MinMax = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_WritErr, %val(1), %val(Status))
	    Endif
	  Endif		! Iostatus Ne zero or Eq zero
	Endif		! Function status is normal

 300	Format (1x/5x,'Failure on open of ',A//10x,'Status(Hex) = ',Z8.8)
 400    Format (1x/5x,'Successful open of ',A,' for write.')

	Return
	End	            
