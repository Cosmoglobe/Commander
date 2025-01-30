C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Print_Eng_Table ( Field, Start, Stop, 
	1	  Trpt_Lun, Nbins, Bins, Eng_Table )

C-------------------------------------------------------------------------------
C
C	Purpose: To print Eng_Table and Percent Table of record 
C	         counts for each engineering bin and scan mode.
C
C	Author: Shirley M. Read
C		STX, January, 1990
C
C	Invocation: Status = FTB_Print_Eng_Table ( Field, Start, Stop,
C		             Trpt_Lun, Nbins, Bins, Eng_Table )
C
CH	Change Log:
CH
C	  ----------------------------------------------------------------------
C
C	Input Files:
C
C	Output Files:
C	  Engineering Table report
C
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Field          C*30           FIRAS engineering field
C	  Start          C*14           Data start time
C	  Stop           C*14           Data stop time
C	  Trpt_Lun       I*4            Logical report unit
C	  Nbins          I*2            Number of bins for engineering values
C	  Bins(5)        R*4            Bin values
C	  Eng_Table(6,4) I*4            Record counts: bins versus scan mode
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	
C	Subroutines Called:
C
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C
C	Processing Method:
C	  Set function return     
C	  Write title, time range and processing information in report file.
C	  Write the engineering table and the percent table in one report file.
C	  Return.
C
C------------------------------------------------------------------------------
C	
	Implicit None

C	Passed Parameters.

	Character*30  Field           ! Engineering field name
	Character*14  Start           ! Data start time
	Character*14  Stop            ! Data stop time
	Integer*4     Trpt_Lun        ! Logical report unit
	Integer*2     Nbins           ! Number of bins for engineering values
        Real*4        Last_Nbins      ! Last Number of bins for engrg. values
	Real*4        Bins(5)         ! Bin values
	Integer*4     Eng_Table(6,4)  ! Record counts: bins versus scan mode

C	Functions
	
	Integer*4     Lib$Locc 

C	Local Declarations.

	Integer*4     Pos, Pos1	      ! Position on line
	Integer*4     Rstatus         ! Write status
        Real*4        Eng_Percent (6,4) ! Engineering Percent Table  
	Integer*4     Zero / 0 /      ! Status value
	Integer*2     Ix, Jx
	Character*71  Dashead
        Character*1   Dasheader(68) / 68 * '_' /
        Equivalence   ( Dashead, Dasheader(1) )
	Character*67  Dash
	Character*1   Dashes(63) / 63 * '_' /
	Equivalence   ( Dash, Dashes(1) )
        Integer*2     Eng_Total (4)  ! Total number of records in each scan mode
        Integer*2     Ind            ! Index

C	External Parameters.

	External FTB_Normal
	External FTB_Aberr

	FTB_Print_Eng_Table = %loc(FTB_Normal)

C	Write title, time range and processing information.

	Pos = Lib$Locc(' ',Field)
	Pos1 = Pos - 1
	    Write (Unit=TRpt_Lun, FMT=100, Iostat=Rstatus)
	1	 Field(1:pos1), Start, Stop
 100        Format (1x//7x,'TABLE of ',a,' FIELD VALUES BINNED BY ',
	1	'SELECTED RANGES '//7x,
	1	'FROM FIRAS ENGINEERING RECORDS SORTED BY MTM SCAN MODE'//
	1	 7x,'COBETRIEVE Time Range Selected Is ',a,' to ',a//)
	    If ( Rstatus .NE. Zero ) Then
	      FTB_Print_Eng_Table = %loc(FTB_Aberr)
	      Call Lib$Signal (%val(Rstatus))
  	    Endif

C	Write Engineering Table subtitle.

	IF ( FTB_Print_Eng_Table .EQ. %loc(FTB_Normal) ) Then
	    Write (Unit=TRpt_Lun, FMT=150, Iostat=Rstatus)
 150        Format (7x,'Eng_Val_Bins    Scan Mode 0  Scan Mode 1  Scan Mode 2', 
	1   '    Scan Mode 3') 
	    Write (Unit=TRpt_Lun, FMT=175, Iostat=Rstatus)
 175        Format (7x,'                    SS           SF           LS ', 
	1   '            LF  ') 
	    Write (Unit=TRpt_Lun, FMT=200, Iostat=Rstatus)

	1     Dashead
 200	    Format (7x,a)
	Endif

C	Write the contents of the Engineering Table.

	IF ( FTB_Print_Eng_Table .EQ. %loc(FTB_Normal) ) Then
	    Do Ix = 1, Nbins
	      Write (Unit=TRpt_Lun, FMT=250, Iostat=Rstatus)
	1       Bins(Ix), (Eng_Table(Ix,Jx),Jx=1,4)
 250 	      Format (7x,'<',2x,F10.3,1x,4(' |',I11),1x,'|'/)
	      Write (Unit=TRpt_Lun, FMT=300, Iostat=Rstatus)
	1       Dash
 300          Format (12x,a)
	    Enddo

C	Write the Last row and its contents of the Engineering Table.

            Last_Nbins=Bins(Nbins)
	    Ix = Nbins+1
            Write (Unit=TRpt_Lun, FMT=350, Iostat=Rstatus)
	1       Last_Nbins, (Eng_Table(Ix,Jx),Jx=1,4)
 350 	    Format (7x,'>',2x,F10.3,1x,4(' |',I11),1x,'|'/)
	      Write (Unit=TRpt_Lun, FMT=400, Iostat=Rstatus)
	1       Dash
 400          Format (12x,a//)
	Endif

C	Write Title and subtitle of the Percent Table.

	IF ( FTB_Print_Eng_Table .EQ. %loc(FTB_Normal) ) Then
	    Write (Unit=TRpt_Lun, FMT=500, Iostat=Rstatus)
 500        Format (27x,'*** PERCENT % TABLE ***'/) 
        Endif
	IF ( FTB_Print_Eng_Table .EQ. %loc(FTB_Normal) ) Then
	    Write (Unit=TRpt_Lun, FMT=550, Iostat=Rstatus)
 550        Format (7x,'Eng_Val_Bins    Scan Mode 0  Scan Mode 1  Scan Mode 2', 
	1   '    Scan Mode 3') 
	    Write (Unit=TRpt_Lun, FMT=575, Iostat=Rstatus)
 575        Format (7x,'                    SS           SF           LS ', 
	1   '            LF  ') 
	    Write (Unit=TRpt_Lun, FMT=600, Iostat=Rstatus)

	1     Dashead
 600	    Format (7x,a)
	Endif

C	Compute Eng_Total from the Engineering table.

        Do Ind = 1,4
           Eng_Total(Ind) = 0
           Do Ix = 1, Nbins + 1
              Eng_Total(Ind) = Eng_Total(Ind) + Eng_Table (Ix,Ind)
           Enddo
        Enddo

C	Compute Eng_Percent from the Engineering table.
        
        Do Ind = 1,4
           Do Ix = 1, Nbins + 1
              If (Eng_Total(Ind) .NE. 0) Then
                  Eng_Percent (Ix,Ind) = ((Float(Eng_Table (Ix,Ind)) 
	1                           / Float(Eng_Total(Ind))) * 100.)
              Else
		  Eng_Total(Ind) = 0
              Endif              
           Enddo
        Enddo

C	Write the contents of the Percent table.

	IF ( FTB_Print_Eng_Table .EQ. %loc(FTB_Normal) ) Then
	    Do Ix = 1, Nbins
	      Write (Unit=TRpt_Lun, FMT=650, Iostat=Rstatus)
	1       Bins(Ix), (Eng_Percent(Ix,Jx),Jx=1,4)
 650 	      Format (7x,'<',2x,F10.3,2x,'|',F9.1,1x,'% |',
	1          F9.1,1x,'% |',F9.1,1x,'% |',F9.1,1x,'% |'/)
	      Write (Unit=TRpt_Lun, FMT=700, Iostat=Rstatus)
	1       Dash
 700          Format (12x,a)
	    Enddo

C	Write the last row and its contents of the Percent Table.

            Last_Nbins=Bins(Nbins)
	    Ix = Nbins+1
            Write (Unit=TRpt_Lun, FMT=750, Iostat=Rstatus)
	1       Last_Nbins, (Eng_Percent(Ix,Jx),Jx=1,4)
 750 	    Format (7x,'>',2x,F10.3,2x,'|',F9.1,1x,'% |',
	1          F9.1,1x,'% |',F9.1,1x,'% |',F9.1,1x,'% |'/)
	      Write (Unit=TRpt_Lun, FMT=800, Iostat=Rstatus)
	1       Dash
 800          Format (12x,a)
	Endif

	Return
	End
