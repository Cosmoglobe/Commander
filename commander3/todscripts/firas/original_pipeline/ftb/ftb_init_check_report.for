
	Integer*4 Function FTB_Init_Check_Report ( Dataset_name,
	1   	     Channel, Start, Stop, Type_Rpt, Checksum, 
	2            Badtime,Status, GRpt, BRpt, GRpt_Lun, BRpt_Lun )
c
c
c-----------------------------------------------------------------------------
c
c       Purpose: To open a report file and write initial information
c                of the output file. 
c                    
c	Written:  N. Gonzales
c		STX
c		March 1990
c
c	Invocation: Status = FTB_Init_Check_Report ( Dataset_name,
c	1   	     Channel, Start, Stop, Type_Rpt, Checksum, 
c	2            Badtime,Status, GRpt, BRpt, GRpt_Lun, BRpt_Lun )
ch Change Log:
ch
c-----------------------------------------------------------------------------
c
c
c 	Input Files:
c
c	Output Files:
c       Good_Checksum_Flag.Rpt, or Bad_Checksum_Flag.Rpt  
c
c	Input Parameters:
c 	Name	  Type      Use           Discription
c-----------------------------------------------------------------------------
c 
c	Character*10 Dataset_Name       ! CT Dataset Name
c       Character*2  Channel            ! FIRAS channel - RH, RL, LH, LL
c	Character*14 Start, Stop        ! Data start and stop times
c	Integer*4    GRpt_Lun, BRpt_Lun ! Logical units
c	Character*1  Type_Rpt           ! Type of checksum file output
c	Logical*1    GRpt  		! Flag to check good checksum data
c	Logical*1    BRpt               ! Flag to check bad checksum data
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
c       Open a report file and write initial information.
c       Return
c-----------------------------------------------------------------------------
	Implicit None

c       Passed Parameters:

	Character*10 Dataset_Name       ! CT Dataset Name
        Character*2  Channel            ! FIRAS channel - RH, RL, LH, LL
	Character*14 Start, Stop        ! Data start and stop times
	Integer*4    GRpt_Lun, BRpt_Lun ! Logical units
	Character*1  Type_Rpt           ! Type of checksum file output
	Logical*1    GRpt  		! Flag to check good checksum data
	Logical*1    BRpt               ! Flag to check bad checksum data
        Logical*1    Checksum           ! set checksum qualifier
       	Logical*1    Badtime            ! set badtime  qualifier

c	Include Files.

	Include '($SSDef)'

c	Functions

	Integer*4 Lib$Get_Lun

c	External Parameters.

	External FTB_Normal
	External FTB_Aberr

c	Local Declarations.

        Integer*4  Zero / 0 /
	Integer*4 Status           !General program status
	Integer*4 Rstatus          !Function return status
	Character*22 Good_Rpt_File /'Good_Checksum_Flag.Rpt'/      
	Character*21 Bad_Rpt_File /'Bad_Checksum_Flag.Rpt'/         
	Character*18 Good_Time_Flag /'Good_Time_Flag.Dat'/      
	Character*17 Bad_Time_Flag  /'Bad_Time_Flag.Dat'/      

c Open report file and write initial information.

	If (Checksum) Then
            If ((Grpt) .and. (Status .eq. SS$_Normal)) Then   
	         Rstatus = Lib$Get_Lun (GRpt_Lun)
	         If ( Rstatus .NE. SS$_Normal ) Then
	              Status = SS$_Abort
	              Call Lib$Signal (%val(Rstatus))
	         Endif  
             Endif
            If ((Brpt) .and. (Status .eq. SS$_Normal)) Then   
	         Rstatus = Lib$Get_Lun (BRpt_Lun)
	         If ( Rstatus .NE. SS$_Normal ) Then
	              Status = SS$_Abort
	              Call Lib$Signal (%val(Rstatus))
	         Endif  
             Endif
            If ((Grpt) .and. (Status .eq. SS$_Normal)) Then   
	         Open (Unit=GRpt_Lun, File=Good_Rpt_File, 
	1        Status='NEW', Iostat=Rstatus)
                 If ( Rstatus .NE. Zero ) Then
	              Status = SS$_Abort
	              Call Lib$Signal (%val(Rstatus))
  	         Endif
             Endif
            If ((Brpt) .and. (Status .eq. SS$_Normal)) Then   
	         Open (Unit=BRpt_Lun, File=Bad_Rpt_File, 
	1        Status='NEW', Iostat=Rstatus)
                 If ( Rstatus .NE. Zero ) Then
	              Status = SS$_Abort
	              Call Lib$Signal (%val(Rstatus))
  	         Endif
             Endif
             If ((Status .eq. SS$_Normal) .and. (Grpt)) Then
	         Write (Unit=GRpt_Lun, FMT=100, Iostat=Rstatus)
	1	       Dataset_Name, Start, Stop
 100             Format (1x//5x,'SCIENCE INTERFEROGRAMS with GOOD ',
	1	'CHECKSUMS for DATASET ',a//
	1	 5x,'COBETRIEVE Time Range Selected Is ',a,' to ',a//)
                 If ( Rstatus .NE. Zero ) Then
	              Status = SS$_Abort
	              Call Lib$Signal (%val(Rstatus))
  	         Endif
              Endif
             If ((Status .eq. SS$_Normal) .and. (Brpt)) Then
 	         Write (Unit=BRpt_Lun, FMT=150, Iostat=Rstatus)
	1	        Dataset_Name, Start, Stop
 150             Format (1x//5x,'SCIENCE INTERFEROGRAMS with BAD ',
	1	 'CHECKSUMS for DATASET ',a//
	1	 5x,'COBETRIEVE Time Range Selected Is ',a,' to ',a//)
                 If ( Rstatus .NE. Zero ) Then
	              Status = SS$_Abort
	              Call Lib$Signal (%val(Rstatus))
  	         Endif
              Endif
	Endif
        If (Badtime) Then
            If ((Grpt) .and. (Status .eq. SS$_Normal)) Then   
	         Rstatus = Lib$Get_Lun (GRpt_Lun)
	         If ( Rstatus .NE. SS$_Normal ) Then
	              Status = SS$_Abort
	              Call Lib$Signal (%val(Rstatus))
	         Endif  
             Endif
            If ((Brpt) .and. (Status .eq. SS$_Normal)) Then   
	         Rstatus = Lib$Get_Lun (BRpt_Lun)
	         If ( Rstatus .NE. SS$_Normal ) Then
	              Status = SS$_Abort
	              Call Lib$Signal (%val(Rstatus))
	         Endif  
             Endif
            If ((Grpt) .and. (Status .eq. SS$_Normal)) Then   
	         Open (Unit=GRpt_Lun, File=Good_Time_Flag, 
	1        Status='NEW', Iostat=Rstatus)
                 If ( Rstatus .NE. Zero ) Then
	              Status = SS$_Abort
	              Call Lib$Signal (%val(Rstatus))
  	         Endif
             Endif
            If ((Brpt) .and. (Status .eq. SS$_Normal)) Then   
	         Open (Unit=BRpt_Lun, File=Bad_Time_Flag, 
	1        Status='NEW', Iostat=Rstatus)
                 If ( Rstatus .NE. Zero ) Then
	              Status = SS$_Abort
	              Call Lib$Signal (%val(Rstatus))
  	         Endif
             Endif
             If ((Status .eq. SS$_Normal) .and. (Grpt)) Then
	         Write (Unit=GRpt_Lun, FMT=200, Iostat=Rstatus)
	1	       Dataset_Name, Start, Stop
 200             Format (1x//5x,'SCIENCE INTERFEROGRAMS with GOOD ',
	1	'TIME FLAGS for DATASET ',a//
	1	 5x,'COBETRIEVE Time Range Selected Is ',a,' to ',a//)
                 If ( Rstatus .NE. Zero ) Then
	              Status = SS$_Abort
	              Call Lib$Signal (%val(Rstatus))
  	         Endif
              Endif
             If ((Status .eq. SS$_Normal) .and. (Brpt)) Then
 	         Write (Unit=BRpt_Lun, FMT=250, Iostat=Rstatus)
	1	        Dataset_Name, Start, Stop
 250             Format (1x//5x,'SCIENCE INTERFEROGRAMS with BAD ',
	1	 'TIME FLAGS for DATASET ',a//
	1	 5x,'COBETRIEVE Time Range Selected Is ',a,' to ',a//)
                 If ( Rstatus .NE. Zero ) Then
	              Status = SS$_Abort
	              Call Lib$Signal (%val(Rstatus))
  	         Endif
              Endif
	Endif
        
       	FTB_Init_Check_Report = Status
        Return
	End
