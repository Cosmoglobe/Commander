
        Program FTB_Binfld

c-----------------------------------------------------------------------------
c
c	Program Name: FTB_Binfld
c
c	Program Description:
c               This program reads records from the FDQ_ENG dataset for
c               a specified time range to generate statistics for selected
c	        engineering fields. The binning of temperatures is selected
c               on the command line, up to 5 bins. The number of records is
c               accumulated for each temperature range by MTM scan mode.
c		A table is generated with these statistics and a histogram
c		is plotted.
c
c	Author: Shirley M. Read
c		STX
c		January 1990
c
c	Calling sequence:
c		invoked from DCL
c
c	Include files:
c		$SSdef
c		Ct$Library:Ctuser.Inc
c
c	Input:
c		None
c
c	Output:
c		None
c
c-----------------------------------------------------------------------------
c
ch Change Log:
ch
c-----------------------------------------------------------------------------

	Implicit None

c	Include Files.

	Include 'CT$Library:CTUser.Inc'
	Include '($SSDef)'

c	Functions

	Integer*4 FTB_Parse_Binfld
	Integer*4 FTB_Fld_Idx
	Integer*4 FTB_Get_Fldval
	Integer*4 FTB_Get_MTM_Mode
	Integer*4 FTB_Print_Eng_Info
	Integer*4 FTB_Accum_Eng_Table
	Integer*4 FTB_Print_Eng_Table
	Integer*4 Lib$Locc
	Integer*4 Lib$Get_Lun	
	Integer*4 Str$Upcase

c	External Parameters.

	External FTB_Normal
	External FTB_Aberr
	External FTB_CTInit
	External FTB_CTOpen
	External FTB_CTRead
	External FTB_CTClos

	Integer*4 CT_Connect_Read
	External  CT_Connect_Read
	External  FUT_Error

c	Data Dictionary and Records.

	Dictionary 'FDQ_ENG'
	Record /FDQ_ENG / Eng_Rec	! FIRAS Engineering records

c	Local Declarations.

	Integer*4  Status		! General program status
	Integer*4  Rstatus              ! Function return status
	Integer*4  Zero / 0 /
	Character*30 Field              ! Selected engineering field
	Character*64 Filename           ! String to hold filename
	Character*10 Dataset_Name 
	Character*7  EDataset / 'FDQ_ENG'/   ! CT Dataset Name
	Character*14 Start, Stop        ! Data start and stop times
	Character*19 Archive / 'CSDR$FIRAS_ARCHIVE:' /
	Integer*2 CT_Status(20)         ! COBETRIEVE status
	Integer*4 Elun, SLun, ERpt_Lun, TRpt_Lun    ! Logical units
	Character*7 Eng_File_1 / 'LatLon_' / 
 	Character*10 Table_File_1 / 'Eng_Table_' /         
	Character*64 Eng_File 
 	Character*64 Table_File 
	Logical*1 Eof / .False. /       ! End of file
	Logical*1 Goodval               ! Good field value
	Integer*4 Pos, Pos1		! Line position
	Integer*2 Index			! Index for field in FUT_Plot_Names
	Integer*2 Scan_Mode(4), Dom_Mode  ! MTM scan modes for 4 chan, dominant
	Integer*2 Nbins                 ! Number of bins for statistics
	Real*4	Bins(5)                 ! Bin ranges, max 5, ascending order
	Real*4  Fieldval	        ! Engineering field value
	Integer*4  Eng_Table(6,4)       ! Count of Eng recs for each bin & mode
	Integer*4  Ix                   ! Index
	Character*1 Type_Rpt            ! Type of timetag file output
	logical*1 Erpt / .False. /, Trpt / .False. /
	Character*1 Semi / ';' /
	Character*1 Slash / '/' /

c Establish the error handler and initialize.

	Call Lib$Establish ( Fut_Error )

	Status = SS$_Normal

c Parse command line.

	Rstatus = FTB_Parse_Binfld ( Field, Start, Stop, Nbins, Bins, Type_Rpt)
	If ( Rstatus .NE. %loc(FTB_Normal) ) Then
	  Status = SS$_Abort
	Else
	  If (Type_Rpt .EQ. 'T') Then
	    Trpt = .True.
	  Endif
	  If (Type_Rpt .EQ. 'E') Then
	    Erpt = .True.
	  Endif
	  If (Type_Rpt .EQ. 'B') Then
	    Erpt = .True.
	    Trpt = .True.
	  Endif
	Endif

c Check if field name is valid and get the index of the field.

	Rstatus = FTB_Fld_Idx ( Field, Index )
	If ( Rstatus .NE. %loc(FTB_Normal) ) Then
	  Status = SS$_Abort
	Endif

	Pos = Lib$Locc(' ',Field)
	Pos1 = Pos - 1

c Open report file.

	If ((Erpt) .and. (Status .eq. SS$_Normal)) Then
	  Rstatus = Lib$Get_Lun (ERpt_Lun)
	  If ( Rstatus .NE. SS$_Normal ) Then
	    Status = SS$_Abort
	    Call Lib$Signal (%val(Rstatus))
	  Endif
	Endif
	If ((Trpt) .and. (Status .eq. SS$_Normal)) Then
	  Rstatus = Lib$Get_Lun (TRpt_Lun)
	  If ( Rstatus .NE. SS$_Normal ) Then
	    Status = SS$_Abort
	    Call Lib$Signal (%val(Rstatus))
	  Endif
	Endif

	Pos = Lib$Locc(' ',Field)
	Pos1 = Pos - 1

	If ( Status .EQ. SS$_Normal ) Then
	  If ( ERpt ) Then
	    Eng_File = Eng_File_1//Field(1:Pos1)//'.'
	1	//Start(1:7)//'_'//Stop(1:7)
	    Open (Unit=ERpt_Lun, File=Eng_File, Status='NEW', 
	1	Iostat=Rstatus)
	    If ( Rstatus .NE. Zero ) Then
	      Status = SS$_Abort
	      Call Lib$Signal (%val(Rstatus))
	    Endif
	  Endif
	  If ( TRpt ) Then
	    Table_File = Table_File_1//Field(1:Pos1)//'.'
	1	//Start(1:7)//'_'//Stop(1:7)
	    Open (Unit=TRpt_Lun, File=Table_File, Status='NEW', 
	1	Iostat=Rstatus)
	    If ( Rstatus .NE. Zero ) Then
	      Status = SS$_Abort
	      Call Lib$Signal (%val(Rstatus))
	    Endif
	  Endif
	Endif	  

c Initialize to COBETRIEVE.
	
	If ( Status .EQ. SS$_Normal ) Then
	  Call CT_Init ( CT_Status )
	  If (Ct_Status(1) .NE. 1 ) Then
	    Rstatus = Ct_Status(1)
	    Call Lib$Signal (FTB_CTInit, %val(1), %val(Rstatus))
	    Status = SS$_Abort
	  Endif
 	Endif

c Get logical unit.

	Rstatus = Lib$Get_Lun ( ELun )
	If ( Rstatus .NE. SS$_Normal ) Then
	  Status = SS$_Abort
	Endif

c Open the COBETRIEVE File.

	If ( Status .EQ. SS$_Normal ) Then
	  Filename = Archive//EDataset//Slash//Start//Semi//Stop//Semi
	  Open (Unit=ELun, File=Filename, Useropen=CT_Connect_Read,
	1	Status='OLD', Iostat=Rstatus )
	  If ( Rstatus .NE. Zero ) Then
	    Call Lib$Signal (FTB_CTOpen, %val(1), %val(Rstatus))
	    Status = SS$_Abort
	  Endif
	Endif

c Read and process the records.

	If ( Status .EQ. SS$_Normal ) Then
	  Eof = .False.
	  Do While ( .NOT. Eof ) 
	      Call CT_Read_Arcv ( ,ELun, Eng_Rec, CT_Status )
	      If (Ct_Status(1) .Eq. CTP_Endoffile) Then
	          Eof = .True.
	      ElseIf ( CT_Status(1) .NE. CTP_Normal ) Then
	          Status = SS$_Abort
	          Rstatus = CT_Status(1)
	          Call Lib$Signal (FTB_CTRead, %val(1), %val(Rstatus))
	          Eof = .True.
	      Endif
	      If ((Status .EQ. SS$_Normal) .and. (.not. Eof)) Then  

		Rstatus = FTB_Get_Fldval( Eng_Rec, Index, Fieldval, Goodval)
		If ( Rstatus .NE. %loc(FTB_Normal)) Status = SS$_Abort

		Rstatus = FTB_Get_MTM_Mode( Eng_Rec, Scan_Mode, Dom_Mode)
		If ( Rstatus .NE. %loc(FTB_Normal)) Status = SS$_Abort

	        If (Erpt .AND. Goodval .AND. Status .EQ. SS$_Normal) Then
		  Rstatus = FTB_Print_Eng_Info  !Print engineering record info
	1	    ( Field, Eng_Rec, Fieldval, Scan_Mode, Dom_Mode, 
	2	     ERpt_Lun, Start, Stop )
	          If ( Rstatus .NE. %loc(FTB_Normal)) Status = SS$_Abort
	        Endif
	        If (Trpt .AND. Goodval .AND. Status .EQ. SS$_Normal) Then
		  Rstatus = FTB_Accum_Eng_Table
	1	    (Fieldval, Dom_Mode, Nbins, Bins, Eng_TabLe)
	          If ( Rstatus .NE. %loc(FTB_Normal)) Status = SS$_Abort
	        Endif
	      Endif
	      If ( Status .Eq. SS$_Abort ) Eof = .True.
	  Enddo 	!Do While ( .NOT. Eof ) 
	Endif		! Status .EQ. SS$_Normal 

c If statistics table requested, print the table.

	If ( Trpt .AND. Status .Eq.SS$_Normal ) Then
		Rstatus = FTB_Print_Eng_Table
	1	  ( Field, Start, Stop, TRpt_Lun, Nbins, Bins,
	2          Eng_TabLe )
	        If ( Rstatus .NE. %loc(FTB_Normal)) Status = SS$_Abort
	Endif

c Close the Archive.

	Call Ct_Close_Arcv( , ELun, Ct_Status)

c Exit with status.

	If ( Status .EQ. SS$_Normal ) Then
	  Call Lib$Signal(FTB_Normal)
	  Call Exit(Status)
	Else
	  Call Lib$Signal(FTB_Aberr)
	  Call Exit(Status)
	Endif	  
	End
