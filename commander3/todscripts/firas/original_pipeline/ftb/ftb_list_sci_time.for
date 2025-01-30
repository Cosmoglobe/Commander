
C-------------------------------------------------------------------------------

	Integer*4 Function FTB_List_Sci_Time 
	1	  ( Channel, Sci_Rec, First, Eof, Rpt_Lun )

C-------------------------------------------------------------------------------
C
C	Purpose: To list the time sequence for NFS raw science records and
C	         other relevant microprocessor information in a compact, one
C	         line per record format. Both transmit time and collect time
c		 will be printed.
C
C	Author: Shirley M. Read
C		STX, March, 1988
C
C	Invocation: Status = FTB_List_Sci_Time 
C		    ( Channel, Sci_Rec, First, Eof, Rpt_Lun )
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
C	  Channel       C*2             FIRAS Channel -- RH, RL, LH or LL
C	  Sci_Rec       NFS_SDF record  FIRAS Raw Science record
C	  First         L*1             First set of records for run
C	  Eof           L*1             Last set of records for run	   
C	  Rpt_Lun       I*4             Report unit number
C
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	
C	Subroutines Called:
C
C	  FTB_Dump_Sci_Data  -- Print routine for FIRAS raw science 
C			        record data
C	  FTB_Sum_Ifg        -- Sum of absolute values for IFG points
C         Time_LT	     -- COBETRIEVE time less than function
C         Time_GT            -- COBETRIEVE time greater than function
C 	  CT_Binary_to_GMT   -- CT binary to GMT conversion
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C
C	Processing Method:
C
C         Set function return status to success.
C         Extract the primary time key from the COBETRIEVE header for the
C	     science record and the midpoint of collect time.
C	  Extract other information needed for the print file from the record.
C	  Check the telemetry quality.
C	  Set the Badtime flag in ASCII.
C	  Compute the sum of absolute values for IFG points.
C	  Print the information for the record.
C	  If a bad status is detected from a call to the dump routine, setgmt
C	     the function status to abort.
C	  Return. 
C
C------------------------------------------------------------------------------
C	
	Implicit None

! 	Passed Parameters.

	Character*2 Channel		! FIRAS channel -- RH, RL, LH, LL
	Dictionary 'NFS_SDF'         
	Record / NFS_SDF / Sci_Rec      ! FIRAS raw science record
	Logical*1 First, Eof            ! First and last calls
	Integer*4 Rpt_Lun               ! Report unit number

! 	Functions.

	Integer*4 FTB_Dump_Sci_Data
	Integer*4 FTB_Sum_IFG
           
!	External Parameters.

	External FTB_Normal
	External FTB_Aberr

!	Local Declarations.

	Integer*4 Time(2)
	Integer*4 Status
	Integer*4 Rstatus
	Integer*4 Start, Stop, Ix
	Integer*4 Absum                 ! Sum of absolute values for IFG Points
	Integer*2 Status_Bits           ! Microprocessor status bits
	Integer*2 Gain                  ! Science gain     
	Integer*2 Speed, Length         ! MTM scan spedd and length
	Integer*2 Sweeps_per_Ifg        ! MTM sweeps per IFG
	Integer*2 Sampcol               !ADC samples collected for collect cycle
	Integer*2 Points                !ADC samples processed for collect cycle
	Integer*2 Sampergr              ! Samples averaged for group

	Integer*4  Transmit_Frame       ! Microprocessor header transmit frame
	Integer*2  Transmit_Word(2)     ! counter
	Equivalence ( Transmit_Frame, Transmit_Word(1) )
	Integer*4  Collect_Frame        ! Microprocessor header collect frame
	Integer*2  Collect_Word(2)      ! counter
	Equivalence ( Collect_Frame, Collect_Word(1) )
	Integer*4  Collect_Time(2)      ! Midpoint of collect time
	Character*14 Collect_Gmt        ! Corresponding Ascii collect time
	Character*14 Gmt                ! Record time in CT Header.
	Byte Qual                       ! Value of telemetry quality
	Byte Badmp                      ! Value of badtime flag
	Character*1 Badtime             ! Value of badtime flag in ASCII
	Character*2 Tlm_Flag / ' 0' /	! Flag indicating telemetry quality:
					! -1 is bad or missing, 0 is good,
					! 1 is questionable
	Character*2 Badval(10) /'0','1','2','3','4','5','6','7','8','9'/
        Character*1 Dash(130)/ 130 * '_'/
	Byte Zero / 0 /
	Byte One / 1 /
	Byte Nine / 9 /
	Byte Xv  / 15 /
	Byte Xvi / 16 /
	Byte Xvii / 17 /
        Logical*4  First_Gmt/.True./
        Integer*4  Nrec/0/        

	FTB_List_Sci_Time = %loc(FTB_Normal)

        If (First_Gmt) Then
          Nrec = 0
          First_Gmt =.False.
        Endif
	If (Eof) Write ( Unit=Rpt_Lun, FMT=100, Iostat=Rstatus)
 100    Format(1x/1x,'End of Science File')

	If ( .Not. Eof ) Then
          Nrec = Nrec + 1
c
	  Do Ix = 1 ,2 
	    Time(Ix) = Sci_Rec.CT_Head.Time(Ix)
	  Enddo	 

	  Gmt = Sci_Rec.CT_Head.Gmt
	  Do Ix = 1, 2
	    Collect_Time(Ix) = Sci_Rec.Collect_Time.Midpoint_Time(Ix)
	  Enddo
	  Call CT_Binary_to_GMT ( Collect_Time , Collect_Gmt )

	  Status_Bits = Sci_Rec.Sci_Head.SC_Head3
	  Gain = Sci_Rec.Sci_Head.Gain
	  Speed = Sci_Rec.Sci_Head.MTM_Speed
	  Length = Sci_Rec.Sci_Head.MTM_Length
	  Sweeps_per_IFG = Sci_Rec.Sci_Head.SC_Head11
	  Sampcol = Sci_Rec.Sci_Head.SC_Head7
	  Points = Sci_Rec.Sci_Head.SC_Head8
	  Sampergr = Sci_Rec.Sci_head.SC_Head9

	  Transmit_Word(1) = Sci_Rec.Sci_Head.SC_Head4  !LSW
	  Transmit_Word(2) = Sci_Rec.Sci_Head.SC_Head5  !MSW
	  Collect_Word(1) = Sci_Rec.Sci_Head.SC_Head12  !LSW
	  Collect_Word(2) = Sci_Rec.Sci_Head.SC_Head13  !MSW

!       Check the telemetry quality.

	  Tlm_Flag = ' 0'
	  Do Ix = 1, 60
	    Qual = Sci_Rec.Sci_Head.Data_Qual(Ix) 	    
	    If ((Qual .Lt. Zero) .Or. (Qual .Eq. Xv) .Or.
	1	(Qual .Gt. Xvii)) Then
	      Tlm_Flag = '-1'
	    Elseif ((Tlm_Flag .Eq. ' 0') .And. ((Qual .Eq. One) .Or.
	1     (Qual .Eq. Xvi) .Or. (Qual .Eq. Xvii))) Then
	      Tlm_Flag = ' 1'
	    Endif
	  Enddo

!	Set the value of the badtime flag.

	  If (( Sci_Rec.Collect_Time.Badtime_Flag .Ge. Zero) .And.
	1    ( Sci_Rec.Collect_Time.Badtime_Flag .Le. Nine)) Then
	    Badmp = Sci_Rec.Collect_Time.Badtime_Flag
	    Badtime = Badval(Badmp + 1)
	  Else
	    Badtime = 'X'
	  Endif

!	Compute the sum of the absolute values of the IFG Points.

	  Absum = FTB_Sum_Ifg ( Sci_Rec.IFG_Data.Ifg )

!	Dump the record information.

	  Status = FTB_Dump_Sci_Data ( Channel, Gmt, Status_Bits, 
	1	Absum, Gain, Speed, Length, Transmit_Frame, Collect_Frame, 
	2	Sweeps_per_IFg, Sampcol, Points, Sampergr, Collect_Gmt, 
	3	Nrec, Tlm_Flag, Badtime, Rpt_Lun)

	Endif 	! Not Eof

	Return
	End
