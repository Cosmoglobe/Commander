
	Integer*4 Function FTB_Init_Report (Dataset_name, Start, 
	1         Stop, Status, Rpt_Lun, Listall, Collect, Pre_It )
c
c
c-----------------------------------------------------------------------------
c
c       Purpose: To Open a report file and write initial information 
c                about the output file.
c                    
c	written by:  N. Gonzales
c		     STX
c		     March 1990
c
c	Invocation: Status = FTB_Init_Report (Dataset_name, Start, 
c                         Stop, Status, Rpt_Lun, Listall, Collect, Pre_It )
c
ch Change Log:
ch
c-----------------------------------------------------------------------------
c
c
c 	Input Files:
c
c	Output Files:
c         FTB_Valid_Time Report
c
c	Input Parameters:
c 	Name	  Type      Use           Discription
c-----------------------------------------------------------------------------
c 
c	Character*10 Dataset_Name       ! CT Dataset Name
c	Character*14 Start, Stop        ! Data start and stop times
c	Integer*4  Rpt_Lun		! Logical units
c	Logical*1  Pre_IT		! Flag to check Pre I&T science data
c					! before time tag problem was fixed
c	Logical*1  Collect              ! Flag to check the midpoint of
c					! collect time
c	Logical*1 Listall               ! Flag to list info for all records
c
c	Otput Parameters:
c 	Name	  Type      Use           Discription
c-----------------------------------------------------------------------------
c
c
c	Subroutine Called:
c
c       Common Variables Used:
c
c 	Name	  Type      Use           Discription
c-----------------------------------------------------------------------------
c
c       Include Files:
c       $SSdef
c
c       Processing Method:
c       Set function return
c        Open a report file and write initial 
c	information about the output file.
c       Return
c-----------------------------------------------------------------------------
	Implicit None

c       Passed Parameters:

	Character*10 Dataset_Name       ! CT Dataset Name
	Character*14 Start, Stop        ! Data start and stop times
	Integer*4  Status		! General program status
	Integer*4 Rpt_Lun		! Logical Report unit
	Logical*1  Pre_IT		! Flag to check Pre I&T science data
					! before time tag problem was fixed
	Logical*1  Collect              ! Flag to check the midpoint of
					! collect time
	Logical*1 Listall               ! Flag to list info for all records

c	Include Files.

	Include '($SSDef)'

c	Functions

	Integer*4 Lib$Get_Lun

c	External Parameters.

	External FTB_Normal
	External FTB_Aberr

c	Local Declarations.

        Integer*4  Zero / 0 /
	Integer*4  Rstatus              ! Function return status
	Character*12 Mid_Rpt /'Midpoint.Rpt'/         
	Character*10 Tlm_Rpt /'Tlm_Mf.Rpt'/         
	Character*12 Rel_Rpt /'Relevant.Rpt'/         

c Open report file and write initial information.

	Rstatus = Lib$Get_Lun (Rpt_Lun)
	If ( Rstatus .NE. SS$_Normal ) Then
	  Status = SS$_Abort
	  Call Lib$Signal (%val(Rstatus))
	Endif
	If ( Status .EQ. SS$_Normal ) Then
          If ( .not. Listall ) Then
             If ( Collect ) then 
	         Open (Unit=Rpt_Lun, File=Mid_Rpt, Status='NEW',
	1        Iostat=Rstatus)
             Else
	         Open (Unit=Rpt_Lun, File=Tlm_Rpt, Status='NEW',
	2        Iostat=Rstatus)
             Endif
         
          Else
	         Open (Unit=Rpt_Lun, File=Rel_Rpt, Status='NEW',
	3        Iostat=Rstatus)
          Endif
       
	  If ( Rstatus .NE. Zero ) Then
	    Status = SS$_Abort
	    Call Lib$Signal (%val(Rstatus))
	  Endif
	Endif	  
	If ( Status .EQ. SS$_Normal ) Then
	  If ( .not. Listall ) Then
	    Write (Unit=Rpt_Lun, FMT=100, Iostat=Rstatus)
	1	 Dataset_Name, Start, Stop
 100        Format (1x//20x,'INVALID RECORD START TIMES for DATASET ',a//
	1	 20x,'COBETRIEVE Time Range Selected Is ',a,' to ',a//)
	  Else
	    Write (Unit=Rpt_Lun, FMT=150, Iostat=Rstatus)
	1	 Dataset_Name, Start, Stop
 150        Format (1x//20x,'LISTING of RECORD START TIMES and '
	1	 'RELEVANT INFORMATION for DATASET ',a//
	2	 20x,'COBETRIEVE Time Range Selected Is ',a,' to ',a//)
	  Endif
	  If ( Rstatus .NE. Zero ) Then
	    Status = SS$_Abort
	    Call Lib$Signal (%val(Rstatus))
  	  Endif
	Endif

	FTB_INIT_REPORT = STATUS        	
	Return
        End
