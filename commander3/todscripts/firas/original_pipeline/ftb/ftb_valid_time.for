        Program FTB_Valid_Time

c-----------------------------------------------------------------------------
c
c	Program Name: FTB_Valid_Time
c
c	Program Description:
c               This program scans a selected NFS dataset for a specified time
c               range to verify that the record times are monotonically 
c	        increasing. If a repeated time or a reverse time is detected,
c               the record time and related information is logged to an output
c               file. The bracketing normal records are also logged to give a
c	        sequential event overview of the bad time periods. The collect
c	        option applies only to the science records which have a computed
c		midpoint of collect time. If the collect flag is set, the
c	        midpoint time field is used in the checking algorithm instead
c		of the COBETRIEVE header time. If the listall flag is set,
c		information from every record is listed.
c
c	Author: Shirley M. Read
c		STX
c		November 1988
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

	Integer*4 FTB_Parse_Valtime
	Integer*4 FTB_Init_Report
	Integer*4 FTB_Check_Sci_Time
	Integer*4 FTB_Check_Pre_IT_Sci_Time
	Integer*4 FTB_List_Sci_Time
	Integer*4 FTB_Check_Eng_Time
	Integer*4 FTB_Check_Hkp_Time
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

	Dictionary 'NFS_SDF'
	Record / NFS_SDF / Sci_Rec(3)	! FIRAS Raw Science records
	Dictionary 'NFS_EMF'
	Record / NFS_EMF / Eng_Rec(3)   ! FIRAS Engineering records
	Dictionary 'NFS_HKP'
	Record / NFS_HKP / Hkp_Rec(3)   ! FIRAS Housekeeping records

c	Local Declarations.

	Integer*4  Status		! General program status
	Integer*4  Rstatus              ! Function return status
	Integer*4  Zero / 0 /
	Character*64 Filename           ! String to hold filename
	Character*10 Dataset_Name       ! CT Dataset Name
	Character*2 Channel             ! FIRAS channel -- RH, RL, LH, LL
	Character*14 Start, Stop        ! Data start and stop times
	Character*19 Archive / 'CSDR$FIRAS_ARCHIVE:' /
	Integer*2 CT_Status(20)         ! COBETRIEVE status
	Integer*4 Lun, Rpt_Lun		! Logical units
	Logical*1 First / .True. /      ! First time
	Logical*1 Eof / .False. /       ! End of file
	Logical*1  Pre_IT		! Flag to check Pre I&T science data
					! before time tag problem was fixed
	Logical*1  Collect              ! Flag to check the midpoint of
					! collect time
	Logical*1 Listall               ! Flag to list info for all records
        Logical*1 Dump                  ! Dump record On/Off switch
	Integer*4  Ix                   ! Index
	Character*1 Type                ! S = sci, E = eng, H = hkp
	Character*1 Semi / ';' /
	Character*1 Slash / '/' /

c Establish the error handler and initialize.

	Call Lib$Establish ( Fut_Error )

	Status = SS$_Normal

c Parse command line.

	Rstatus = FTB_Parse_Valtime ( Dataset_Name, Start, Stop, 
	1	  Pre_IT, Collect, Listall, Dump )
	If ( Rstatus .NE. %loc(FTB_Normal) ) Then
	  Status = SS$_Abort
	Else
	  Rstatus = Str$Upcase(Dataset_Name, Dataset_Name)
	Endif

c Open report file and write initial information.

	IF (Status .EQ. SS$_Normal) Then

        Rstatus = FTB_Init_Report (Dataset_Name, Start, Stop, 
	1         Status, Rpt_Lun, Listall, Collect, Pre_It) 

	Else

       	    If (Rstatus .NE. %loc(FTB_Normal)) Then
                Status = SS$_Abort
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

c Get a logical unit.

	Rstatus = Lib$Get_Lun ( Lun )
	If ( Rstatus .NE. SS$_Normal ) Then
	  Status = SS$_Abort
	Endif

c Open the File.

	If ( Status .EQ. SS$_Normal ) Then
	  Filename = Archive//Dataset_Name//Slash//Start//Semi//Stop//Semi
	  Open (Unit=Lun, File=Filename, Useropen=CT_Connect_Read,
	1	Status='OLD', Iostat=Rstatus )
	  If ( Rstatus .NE. Zero ) Then
	    Call Lib$Signal (FTB_CTOpen, %val(1), %val(Rstatus))
	    Status = SS$_Abort
	  Endif
	Endif

c Read and process the records.

	If ( Dataset_Name(1:7) .EQ. 'NFS_HKP' ) Then
	  Type = 'H'
	Else
	  Channel(1:2) = Dataset_Name(9:10)
	  If ( Dataset_Name(1:7) .EQ. 'NFS_SDF'
	1      .OR. Dataset_name(1:7) .EQ. 'FPP_SDF'
	2      .OR. Dataset_name(1:7) .EQ. 'FDQ_SDF'  ) Then
	    Type = 'S'
	    If ( .Not. Listall ) Then
	      If ( Collect ) Then
	        Write (Unit=Rpt_Lun, FMT=200, Iostat=Rstatus)
 200            Format (1x//20x,'Time Check is Performed on the Computed ',
	1	'Midpoint of Collect Times for the Interferogram'//)

	        Write (Unit=Rpt_Lun, FMT=210, Iostat=Rstatus)
 210            Format (1x//20x,'Note !!! **** The Original Record is',
	1       ' included within the Duplicate Records ****'//) 

	      Else
	        Write (Unit=Rpt_Lun, FMT=250, Iostat=Rstatus)
 250            Format (1x//20x,'Time Check is Performed on the Telemetry',
	1	' Minor Frame Transmit Times for the Interferogram'//)
	      Endif
	    Endif
	  Elseif ( Dataset_Name(1:7) .EQ. 'NFS_EMF' ) Then
	    Type = 'E'
	  Endif
	Endif
	If ( Status .EQ. SS$_Normal ) Then
	  Eof = .False.
	  Do While ( .NOT. Eof ) 
	    If ( First .AND. (.Not. Listall)) Then
              Do Ix = 1, 2
	        If ( Type .EQ. 'S' ) Then
	          Call CT_Read_Arcv ( ,Lun, Sci_Rec(Ix), CT_Status )
	        Elseif ( Type .EQ. 'E' ) Then
	          Call CT_Read_Arcv ( ,Lun, Eng_Rec(Ix), CT_Status )
	        Elseif ( Type .EQ. 'H' ) Then
	          Call CT_Read_Arcv ( ,Lun, Hkp_Rec(Ix), CT_Status )
	        Endif
	        If ( CT_Status(1) .NE. 1 ) Then
	          Status = SS$_Abort
	          Rstatus = CT_Status(1)
	          Call Lib$Signal (FTB_CTRead, %val(1), %val(Rstatus))
	          Eof = .True.
	        Endif
	      Enddo
	      If ( Status .EQ. SS$_Normal ) Then
	        If ( Type .EQ. 'S' ) Then
	          If ( Pre_IT ) Then
	            Rstatus = FTB_Check_Pre_IT_Sci_Time 
	1		( Channel, Sci_Rec, First, Eof, Rpt_Lun, Dump)
	          Else
	            Rstatus = FTB_Check_Sci_Time 
	1		( Channel, Sci_Rec, First, Eof, Rpt_Lun,
	1		  Collect, Dump)
	          Endif
	        Elseif ( Type .EQ. 'E' ) Then
	          Rstatus = FTB_Check_Eng_Time 
	1		( Channel, Eng_Rec, First, Eof, Rpt_Lun , Dump)
	        Else
	          Rstatus = FTB_Check_Hkp_Time 
	1		( Hkp_Rec, First, Eof, Rpt_Lun , Dump)
	    	Endif
	        If ( Rstatus .NE. %loc(FTB_Normal)) Status = SS$_Abort
	      Endif
	      First = .False.
	    Else		! Listall or Not First Time
	      If ( Type .EQ. 'S' ) Then
	        Call CT_Read_Arcv ( ,Lun, Sci_Rec(3), CT_Status )
	      Elseif ( Type .EQ. 'E' ) Then
	        Call CT_Read_Arcv ( ,Lun, Eng_Rec(3), CT_Status )
              Elseif ( Type .EQ. 'H' ) Then
	        Call CT_Read_Arcv ( ,Lun, Hkp_Rec(3), CT_Status )
	      Endif
	      If ( CT_Status(1) .EQ. CTP_Endoffile ) Eof = .True.
	      If (( CT_Status(1) .EQ. 1 ) .OR. 
	1	   (CT_Status(1) .EQ. CTP_Endoffile)) Then
	        If ( Type .EQ. 'S' ) Then
		  If ( Listall ) Then
	            Rstatus = FTB_List_Sci_Time 
	1		( Channel, Sci_Rec(3), First, Eof, Rpt_Lun )
	          Elseif ( Pre_IT ) Then
	            Rstatus = FTB_Check_Pre_IT_Sci_Time 
	1		( Channel, Sci_Rec, First, Eof, Rpt_Lun, Dump)
	          Else
	            Rstatus = FTB_Check_Sci_Time 
	1		( Channel, Sci_Rec, First, Eof, Rpt_Lun,
	1		  Collect, Dump )
	          Endif
	        Elseif ( Type .EQ. 'E' ) Then
	          Rstatus = FTB_Check_Eng_Time 
	1		( Channel, Eng_Rec, First, Eof, Rpt_Lun, Dump)
	        Else
	          Rstatus = FTB_Check_Hkp_Time 
	1		( Hkp_Rec, First, Eof, Rpt_Lun, Dump)
	    	Endif
	        If ( Rstatus .NE. %loc(FTB_Normal)) Status = SS$_Abort
	      Else  
	        Status = SS$_Abort
	        Rstatus = CT_Status(1)
	        Call Lib$Signal (FTB_CTRead, %val(1), %val(Rstatus))
	        Eof = .True.
	      Endif
	    Endif		! Not First Time
	  Enddo 	!Do While ( .NOT. Eof ) 
	Endif		! Status .EQ. SS$_Normal 

c Close the Archive.

	Call Ct_Close_Arcv( , Lun, Ct_Status)

c Exit with status.

	If ( Status .EQ. SS$_Normal ) Then
	  Call Lib$Signal(FTB_Normal)
	  Call Exit(Status)
	Else
	  Call Lib$Signal(FTB_Aberr)
	  Call Exit(Status)
	Endif	  
	End
