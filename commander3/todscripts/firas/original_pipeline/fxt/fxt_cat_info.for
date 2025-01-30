C-------------------------------------------------------------------------------

	Integer*4 Function FXT_Cat_Info ( Gmt_Start, Gmt_Stop,
	1	  File_Seg, Report, Lun_Rpt, Time_Range, Data_Start,
	2	  Data_Stop, Num_Seg, Data_Type )

C-------------------------------------------------------------------------------
C
C	Purpose: To convert engineering analogs from counts to engineering
C		 units. 
C
C	Author: Shirley M. Read
C		STX, November, 1988
C
C	Invocation: Status = FXT_Cat_Info ( Gmt_Start, Gmt_Stop,
C		  File_Seg, Report, Lun_Rpt, Time_Range, Data_Start,
C		  Data_Stop, Num_Seg, Data_Type )
C
CH	Change Log:
CH
C	  ----------------------------------------------------------------------
C
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Gmt_Start     C*14            Selected Gmt start time
C	  Gmt_Stop      C*14            Selected Gmt stop time
C	  File_Seg      C*39            Selected FDQ_Eng filename
C	  Report        L*1             Flag to write report
C	  Lun_Rpt       I*4             Logical unit for report
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C         Time_Range    C*30            Complete data time range for segments
C 	  Data_Start(2) I*4             Binary data start time
C 	  Data_Stop(2)  I*4             Binary data stop time
C 	  Num_Seg       I*4             Actual number of segments in time range
C	  Data_Type     C*2		Source type for data
C	
C	Subroutines Called:
C
C	  CT_Init
C	  CCT_Query_Catalog
C	  CCT_Query_TOD_Catalog
C	  CT_Gmt_to_Binary
C	  CT_Binary_to_Gmt
C	  Time_GT
C 	  SYS$Asctim
C	  Lib$Locc
C	  Lib$Signal
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C	  FXT_Msg.Txt
C	  CT$Library:CTUser.Inc
C
C	Processing Method:
C
C	  This routine will query the COBETRIEVE catalog to find out what 
C	  segment or segments are to be processed. If more than one segment is
C	  selected, it will build the filenames for the FDQ_Eng segments to be
C	  processed from the returned Catalog_Rec and write it to the report 
C         file. The actual start and stop times of the data are determined from
C         the Catalog_Rec if one file has been selected on the command line.
C	  Otherwise, the user selected start and stop times will be used.
C	  If the status is OK, the number of segments for the dataset is stored 
C	  in Num_Seg. If the number of segments exceeds the maximum allowed 
C	  number of entries, a message is signalled, but no further action is
C	  taken.
C
C	  < PDL for this routine > :
C
C         Set the function status to normal.
C	  Call Ct_Init to initialize to COBETRIEVE.
C	  If a specific FDQ_Eng segment requested, then
C	    Call CCT_Query_Catalog_Record.
C           If return status is error or no catalog record was found, then
C	      Signal error and set function status to abort.
C	    Elseif OK, then
C	      Determine the data start and stop times from the catalog record.
C	      Convert the data times to an ASCII time range for the data.
C	      Set Num_Seg to one.
C	      If a report is requested, log the filename and time range.
C	    Endif
C	  Elseif a time range is used, then 
C	    Call CCT_Query_TOD_Catalog with time range to get all FDQ_Eng
C             segments from the catalog entries.
C	    If the number of entries = 0 or an error is returned, then
C	      Signal error and set the function status to abort.
C	    Elseif there are too many segments for the catalog array, then
C	      Signal an information message.
C	    Endif
C	    If status is OK, then
C	      Store the file extension for the first FDQ_Eng segment in 
C	        the catalog record and extract the data type to be used in
C	        the file extension name of the output FXT_Eng_Xtrm file.
C	      If a report is requested, Then
C 		Build the FDQ_Eng filename for each undeleted segment from 
C		  the catalog and log it to the report file.
C	      Endif
C	      Convert the ASCII GMT start and stop to binary start and stop.
C	      Store the number of segments to be processed.	
C	      If a report is requested, Then
C               Log the names of the undeleted files in the time range into 
C                 the report.
C	        Log the data time range and number of segments.
C	      Endif
C           Endif
C	    If the number of segments is one, set it to Num_for_Time_Range.
C	  Endif
C	  Return.
C
C------------------------------------------------------------------------------
C	
	Implicit None

!	Passed Parameters

	Character*14	Gmt_Start       ! Selected Gmt start time
	Character*14	Gmt_Stop        ! Selected Gmt stop time
	Character*39	File_Seg        ! Selected FDQ_Eng filename
	Logical*1	Report          ! Flag to write report
	Integer*4 	Lun_Rpt         ! Logical unit for report
	Character*30	Time_Range    	! Complete data time range for segments
	Integer*4	Data_Start(2)   ! Binary data start time
	Integer*4      	Data_Stop(2)    ! Binary data stop time
	Integer*4 	Num_Seg         ! Number of segments for processing
	Character*2     Data_Type       ! Source type for data

!	Include files and external parameters.

	Include    '(FXT_Msg)'
	Include	   '($SSDef)'
	Include    'CT$Library:CTUser.Inc'

	External   Cut_No_Cat_Records_Fnd
	external   cct_q_no_cat_entry
	integer*4  cct_q_no_cat_entry
	Integer*4  	Max_Entries	! Maximum # of segments allowed
	Parameter  	(Max_entries = 50)	

	Dictionary 'CCM_CME_Catalog_Entry'	!Info retrieved from CT
	Record     /CCM_CME_Catalog_Entry/Catalog_Rec(Max_Entries)

	Include '(CCT_Query_Catalog_Record)'
	Record  /Query_Catalog/Query_Catalog_Rec

	Include '(CCT_Query_TOD_Catalog_Record)'
	Record  /Query_TOD_Catalog/Query_TOD_Catalog_Rec

!	Functions

	Integer*4	CCT_Query_Catalog
	Integer*4	CCT_Query_TOD_Catalog
	Integer*4       Time_GT
	Integer*4	SYS$Asctim
	Integer*4       Lib$Locc

!	Local Declarations

	Character*14    Arc_Id /'CSDR$FIRAS_RAW'/
	Character*7     Dataset / 'FDQ_ENG' /
	Character*2     DT              ! Data type
	Character*39	Filename        ! Name of FDQ_Eng file in time range
	Character*39    Blank           ! Blanks
	Character*1     Blanks(39) / 39 * ' ' /
	Equivalence     ( Blanks(1), Blank )
	Integer*2       CT_Status(20)   ! COBETRIEVE return status
	Integer*4	End_Bntm(2)
	Integer*4 	Start_Bntm(2)
	Integer*4       Num_Cat_Rec     ! Number of catalog records from CT
	Integer*4       Ipos            ! Character position in string
	Integer*4	Iext		!Number of characters in file extension
	Integer*4	Cur_Info	! Index for catalog records
	Integer*4       File_Ptr        ! Index for files
	Integer*4       First_File      ! First file to process
	Integer*4       Count		! COunter for number of segments 
	Integer*2       Ix              ! Index
	Integer*4       Zero / 0 /      ! Zero value
	Integer*4       One / 1 /       ! One value
	Integer*4	Status		! Return status variable
	Logical*1       Good_Set(50)    ! Flag for good data set
	Data            Good_Set / 50 * .False. /
	Logical*1       Found / .False. /
	Integer*4       Num_for_Time_Range / 1000 /  ! Flag for number of 
						     ! segments if one segment

!	Set function status to normal.

	FXT_Cat_Info = %loc(FXT_Normal)

!	Initialize to COBETRIEVE.

	Call CT_Init(CT_Status)
	If ( CT_Status(1). Ne. CTP_Normal ) Then
	  FXT_Cat_Info = %loc(FXT_Aberr)
	  Status = Ct_Status(1)
	  Call Lib$Signal(FXT_CTInit, %val(1), %val(Status))
	Endif

!     A specific FDQ_Eng segment has been specified by the user.      
!     Call the Query COBETRIEVE Catalog routine for the information.

	IF (( File_Seg .Ne. Blank ) .And. 
	1   ( FXT_Cat_Info .Eq. %loc(FXT_Normal) )) Then
	
!	  Get catalog information for the FDQ_Eng file.

	  Query_Catalog_Rec.Archive_Id = Arc_Id
	  Query_Catalog_Rec.Filename(1:39) = File_Seg(1:39)
	  Status = CCT_Query_Catalog( Query_Catalog_Rec,
	1            Catalog_Rec)
	  If (.Not. Status ) Then		! Failed to get cat rec
	    FXT_Cat_Info = %loc(FXT_Aberr)
	    If (Status .Eq. %loc(cct_q_no_cat_entry)) Then
	      Call Lib$Signal( FXT_NoCatRec, %val(1), File_Seg)
	    Else
	      Call Lib$Signal(FXT_QCatErr,%val(2),%val(Status))
	    Endif
	  Elseif (Status .Eq. %loc(cct_q_no_cat_entry)) Then
	      Call Lib$Signal( FXT_NoCatRec, %val(1), File_Seg)
	      FXT_Cat_Info = %loc(FXT_Aberr)
	  Else			! Call to CCT_query_catalog is OK

!	Save the time range for the data in the segment.

	    Do Ix = 1, 2
	      Data_Start(Ix) = Catalog_Rec(1).Initial_Time(Ix)
	    Enddo
	    Do Ix = 1, 2
	      Data_Stop(Ix) = Catalog_Rec(1).Final_Time(Ix)
	    Enddo
	    Call CT_Binary_To_GMT ( Data_Start, Time_Range(1:14))
	    Time_Range(15:15) = ';'
	    Time_Range(16:29) = Catalog_Rec(1).Final_Time_Key(1:14)
	    Time_Range(30:30) = ';'

!	Save the number of segments.

	    Num_Seg = 1

!	Save the data type.

	    Data_Type(1:2) = 
	1     Catalog_Rec(1).Filename_Extension(1:2)

	    If (Report) Then
	      Write(Unit=Lun_Rpt, FMT=100, Iostat=Status)
	1	File_Seg, Time_Range(1:14), Time_Range(16:29), Data_Type
	      If ( Status .Ne. Zero ) Then
		FXT_Cat_Info = %loc(FXT_Aberr)
	        Call Lib$Signal(FXT_WritErr, %val(1), %val(Status))
              Endif
	    Endif
	  Endif

!	End of one segment selected.

	ELSEIF (FXT_Cat_Info .Eq. %loc(FXT_Normal)) Then

!     A time range was selected by the user.                           
!     Convert GMT start and stop to start and stop VAX quadwords.               
!     Query COBETRIEVE TOD catalog to get files in time range.          

	  Call CT_Gmt_to_Binary (Gmt_Start, Start_Bntm)
	  Call CT_Gmt_to_Binary (Gmt_Stop, End_Bntm)
	  Query_Tod_Catalog_Rec.Archive_Id = Arc_Id 
	  Query_Tod_Catalog_Rec.Dataset_Name = Dataset
	  Do Ix = 1, 2
	    Query_Tod_Catalog_Rec.Start_Time(Ix)= Start_Bntm(Ix)
            Query_Tod_Catalog_Rec.Stop_Time(Ix) = End_Bntm(Ix) 
	  Enddo

!	  Get catalog information for the FDQ_Eng files.

	  Status = CCT_Query_Tod_Catalog( Query_Tod_Catalog_Rec,
	1	   Catalog_Rec, Num_Cat_Rec)

!	  Check return status. If bad, signal error and set function to abort.

	  If (.Not. Status ) Then  ! Failed to get cat rec
            FXT_Cat_Info = %loc(FXT_Aberr)
	    If (Status .Eq. %loc(Cut_No_Cat_Records_Fnd)) Then
	      Call Lib$Signal(FXT_NoCatRecTOD, %val(1), %val(Status)) 
	    Else		! Status is bad from CCT_query_catalog
	      Call Lib$Signal(FXT_QTODCatErr,%val(1),%val(Status))
	    Endif		! Not status from CCT_query_catalog

	  Elseif (Status .Eq. %loc(Cut_No_Cat_Records_Fnd)) Then
            FXT_Cat_Info = %loc(FXT_Aberr)
	    Call Lib$Signal(FXT_NoCatRecTOD, %val(1), %val(Status)) 

	  Elseif ( Num_Cat_Rec .Gt. Max_Entries ) Then

!	Signal message if there are too many catalog entries for query.
!       This should not affect the reading of the FDQ_ENG records.

	    Call Lib$Signal(FXT_TooManyCat,%val(1),%val(Num_Cat_Rec))
	  Endif     ! Status from CCT_query_TOD_catalog

	  If ( FXT_Cat_Info .Eq. %loc(FXT_Normal)) Then	  ! Good query status

!	Print the names of the files to be processed. Find the first segment
!       in the time range which does not have a delete flag set. Store this
!	data type for use in the file extension name for the FXT_Eng_Xtrm file.

	    File_Ptr = 0
	    Do While ((.Not. Found) .And. (File_Ptr .Lt. Num_Cat_Rec))
	      File_Ptr = File_Ptr + 1
	      Ipos = Lib$Locc ( ' ',
	1	     Catalog_Rec(File_Ptr).Filename_Extension)
	      If ( Ipos .Eq. Zero ) Then
	        Iext = 20
	      Else
	        Iext = Ipos - 1
	      Endif
	      Filename(1:39) = Blank(1:39)
	      If (Catalog_Rec(File_Ptr).Deletion_Flag .Eq. Zero) Then
		Found = .True.
	        First_File = File_Ptr
		Good_Set(File_Ptr) = .True.
	          Data_Type(1:2) = 
	1	  Catalog_Rec(File_Ptr).Filename_Extension(1:2)
	        Filename = Dataset // '.' //
	1	  Catalog_Rec(File_Ptr).Filename_Extension(1:Iext) 
	2	  // ';' // Catalog_Rec(File_Ptr).Version_Number
		If (Report) Then
		  Write(Unit=Lun_Rpt,FMT=200,Iostat=Status) Filename
 	        Endif
	      Elseif (Report) Then
		Filename = Dataset//'.'//
	1	  Catalog_Rec(File_Ptr).Filename_Extension//';'
	2	  //Catalog_Rec(File_Ptr).Version_Number
	        Write(Unit=Lun_Rpt,FMT=300,Iostat=Status) Filename
	      Endif
	      If (.Not. Found) File_Ptr = File_Ptr + 1
            Enddo
	    If (.Not. Found ) Then
	      FXT_Cat_Info = %loc(Fxt_Aberr)
	      Call Lib$Signal(FXT_AllCatRecDel)
	    Endif
	    If ((Report) .And. (Status .Ne. Zero)) Then
	      FXT_Cat_Info = %loc(FXT_Aberr)
	      Call Lib$Signal(FXT_WritErr, %val(1), %val(Status))
	      Found = .False.
            Endif

!	If found, convert the user selected time range for the FDQ_ENG files
!	to binary start and stop times. Set the time range to GMT times.

	    If ((Found) .And. (FXT_Cat_Info .Eq. %loc(FXT_Normal))) Then

	        Call CT_GMT_To_Binary (GMT_Start, Data_Start)
	        Call CT_GMT_To_Binary (GMT_Stop, Data_Stop)

	        Time_Range(1:14) = GMT_Start(1:14)
	        Time_Range(15:15) = ';'
	        Time_Range(16:29) = GMT_Stop(1:14)
	        Time_Range(30:30) = ';'
	        Num_Seg = 0

		Cur_Info = File_Ptr 

	        Do While ((Cur_Info .Lt. Num_Cat_Rec) .And. 
	1	  (Cur_Info .Lt. Max_Entries))
		  Cur_Info = Cur_Info + 1
	          Filename(1:39) = Blank(1:39)
	          Ipos = Lib$Locc ( ' ',
	1	     Catalog_Rec(Cur_Info).Filename_Extension)
	          If ( Ipos .Eq. Zero ) Then
	             Iext = 20
	          Else
	             Iext = Ipos - 1
	          Endif
	          Filename = Dataset // '.' //
	1	    Catalog_Rec(Cur_Info).Filename_Extension(1:Iext) 
	2	    // ';' // Catalog_Rec(File_Ptr).Version_Number
	          If (Catalog_Rec(Cur_Info).Deletion_Flag .Eq. Zero) Then
		    Good_Set(Cur_Info) = .True.
		    If (Report) Then
		      Write(Unit=Lun_Rpt,FMT=200,Iostat=Status) Filename
 	            Endif
	          Elseif (Report) Then
	            Write(Unit=Lun_Rpt,FMT=300,Iostat=Status) Filename
	          Endif
	        Enddo

	        Do While ((Count .Lt. Num_Cat_Rec) .And. 
	1	  (Count .Lt. Max_Entries))
		  Count = Count + 1
	          If (Good_Set(Count)) Then
	  	    Num_Seg = Num_Seg + 1
	          Endif
	        Enddo

	        If ( Report ) Then
	          Write(Unit=Lun_Rpt, FMT=400, Iostat=Status) 
	1	    Time_Range(1:14), Time_Range(16:29), Data_Type, Num_Seg
	          If ( Status .Ne. Zero ) Then
	            FXT_Cat_Info = %loc(FXT_Aberr)
	            Call Lib$Signal(FXT_WritErr, %val(1), %val(Status))
                  Endif
		Endif   ! Report
	        If (Num_Seg .Eq. One) Num_Seg = Num_for_Time_Range
	    Endif 	! Found and function status is good
	  Endif		! Good status from CCT_query_catalog
	ENDIF		! Time range specified.

 100    Format(1x/5x,'Segment from Catalog to be processed: ',
	1     A //10x,'Time Range of Data in Segment: ',A,' to ',A,
	2     ' Data_Type: ',A)
 200    Format(1x/5x,'Segment from Catalog to be processed: ',A )
 300    Format(1x/5x,'Segment from Catalog to be omitted: ', A //
	1     10x,'The deletion flag is set in the catalog record.')
 400    Format(1x/5x,'Time Range for Data in Current Run: ',A,' to ',A//
	1     10x,' Data Type: ', A, '. Number of catalog segments: ',I10)
	
	Return
	End
 
