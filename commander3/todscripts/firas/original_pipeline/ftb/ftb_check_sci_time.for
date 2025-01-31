
C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Check_Sci_Time 
	1	  ( Channel, Sci_Rec, First, Eof, Rpt_Lun,
	2	    Collect, Dump )

C-------------------------------------------------------------------------------
C
C	Purpose: To check the time sequence for NFS raw science records and
C	         log all records with repeating time keys or time keys which
C		 occur in reverse chronological order. If the collect flag is
C		 set, the midpoint time fields will be used for the testing.
C
C	Author: Shirley M. Read
C		STX, November, 1988
C
C	Invocation: Status = FTB_Check_Sci_Time 
C		    ( Channel, Sci_Rec, First, Eof, Rpt_Lun, Collect, Dump )
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
C	  Sci_Rec(3)    NFS_SDF record  FIRAS Raw Science records
C	  First         L*1             First set of records for run
C	  Eof           L*1             Last set of records for run	   
C	  Rpt_Lun       I*4             Report unit number
C	  Collect       L*1             Flag to check the midpoint of collect
C         Dump          L*1             Dump Record On/Off switch
C
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	
C	Subroutines Called:
C
C	  FTB_Dump_Sci_Time  -- Print routine for FIRAS raw science record times
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
C	  Set the start and stop pointers for the raw science records in core
C            according to whether it is the first call or not.
C         Extract the primary time key from the COBETRIEVE header for the
C	     start and stop science records.
C	  Extract other information needed for the print file from the record.
C	  Compare the time keys to see if the records are in ascending order.
C	     If the time is in ascending order, dump the start record if it's
C	        time was bad on a previous call or it marks the end of a block
C               with invalid time sequencing, and set flags for the next call.
C	     If the time is in reverse order, dump the start record with 
C                relevant information and set flags for the next call.
C	     If the time is repeated, dump the start record only if it is the 
C	         last record, increment the repeated record counter and set
C	         flags for the next call.
C	  If a bad status is detected from a call to the dump routine, set
C	     the function status to abort.
C	  Return. 
C
C------------------------------------------------------------------------------
C	
	Implicit None

! 	Passed Parameters.

	Character*2 Channel		! FIRAS channel -- RH, RL, LH, LL
	Dictionary 'NFS_SDF'         
	Record / NFS_SDF / Sci_Rec(3)   ! FIRAS raw science records
	Logical*1 First, Eof            ! First and last calls
	Integer*4 Rpt_Lun               ! Report unit number
	Logical*1 Collect               ! Flag to check the collect time
        Logical*1 Dump                  ! Dump Record On/Off Switch

! 	Functions.

	Logical*1 Time_LT		! CT Time less than
	Logical*1 Time_GT               ! CT Time grater then
	Integer*4 FTB_Dump_Sci_Time
	Integer*4 FTB_Sum_IFG
	Integer*4 FTB_If_Midpoint_Time

!	External Parameters.

	External FTB_Normal
	External FTB_Aberr

!	Local Declarations.

	Integer*4 Time1(2)
	Integer*4 Time2(2)
	Integer*4 Status
	Integer*4 Rstatus
	Integer*4 Start, Stop, Ix
	Integer*4 Numrec / 1 /
	Character*2 Time_Flag / 'BB' /	! Flag indicating reason for printing
					! the record. BB = begin block, 
					! EB = end block, EQ = time equal,
 					! LT = time less than.
	Logical*1 Fault / .False. /     ! Bad time sequence detected
	Logical*1 Prtnxt / .False. /    ! Log next record as end of block
	Integer*4 Absum                 ! Sum of absolute values for IFG Points
	Integer*2 Status_Bits           ! Microprocessor status bits
	Integer*2 Speed, Length         ! MTM scan spedd and length
	Integer*2 Sweeps_per_Ifg        ! MTM sweeps per IFG
	Integer*2 Neg_Sweeps            ! Negative of sweeps per IFG
	Integer*2 Points                ! Points process for collect cycle
	Integer*4  Transmit_Frame       ! Microprocessor header transmit frame
	Integer*2  Transmit_Word(2)     ! counter
	Equivalence ( Transmit_Frame, Transmit_Word(1) )
	Integer*4  Collect_Frame        ! Microprocessor header collect frame
	Integer*2  Collect_Word(2)      ! counter
	Equivalence ( Collect_Frame, Collect_Word(1) )
	Integer*4  Collect_Time(2)      ! Computed collect time
	Character*14 Collect_Gmt        ! Corresponding Ascii collect time
	Character*14 Gmt                ! Record time in CT Header.
	Character*14 Contime            ! GMT converted header time
	Logical*1 Check / .False. /
        Character*1 Dash(130)/ 130 * '_'/
c
        Logical*4  First_Gmt/.True./
        Integer*4  Nrec/0/        
	Integer*4  iflag

	FTB_Check_Sci_Time = %loc(FTB_Normal)
	Check = .False.

	If ( First ) Then
	  Start = 1
	  Stop = 2
c	  iflag = sci_rec(start).collect_time.badtime_flag
	Else
	  Start = 2
	  Stop = 3
c	  iflag = sci_rec(start).collect_time.badtime_flag
	Endif

	Do Ix = 1 ,2 
	  If (Collect) Then
	    Time1(Ix) = Sci_Rec(Start).Collect_Time.Midpoint_Time(Ix)
	    Time2(Ix) = Sci_Rec(Stop).Collect_Time.Midpoint_Time(Ix)
	  Else
	    Time1(Ix) = Sci_Rec(Start).CT_Head.Time(Ix)
	    Time2(Ix) = Sci_Rec(Stop).CT_Head.Time(Ix)
	  Endif
	Enddo	 

c
c *** initial the first_gmt record number
c *** QCC/STX 12/22/88
C
        if (first_gmt) then
        nrec = 0
        first_gmt =.false.
        endif
c
	Gmt = Sci_Rec(Start).CT_Head.Gmt
	Do Ix = 1, 2
	   Collect_Time(Ix) = Sci_Rec(Start).Collect_Time.Midpoint_Time(Ix)
	Enddo
	Call CT_Binary_to_GMT ( Collect_Time , Collect_Gmt )
	Status_Bits = Sci_Rec(Start).Sci_Head.SC_Head3
	Speed = Sci_Rec(Start).Sci_Head.MTM_Speed
	Length = Sci_Rec(Start).Sci_Head.MTM_Length
	Sweeps_per_IFG = Sci_Rec(Start).Sci_Head.SC_Head11
	Points = Sci_Rec(Start).Sci_Head.SC_Head8

	Transmit_Word(1) = Sci_Rec(Start).Sci_Head.SC_Head4  !LSW
	Transmit_Word(2) = Sci_Rec(Start).Sci_Head.SC_Head5  !MSW
	Collect_Word(1) = Sci_Rec(Start).Sci_Head.SC_Head12  !LSW
	Collect_Word(2) = Sci_Rec(Start).Sci_Head.SC_Head13  !MSW

!	Compute the sum of the absolute values of the IFG Points.

	Absum = FTB_Sum_Ifg ( Sci_Rec(Start).IFG_Data.Ifg )

	If ( Time_GT( Time2, Time1 )) Then	! Time OK on current records 
	  If (Eof) Time_Flag = 'EF'
	  If ( Fault .OR. Eof  ) THEN ! Previous time seq. bad
	    Status = FTB_Dump_Sci_Time ( Channel, Collect_Gmt, Status_Bits, 
	1	Absum, Speed, Length, Transmit_Frame, Collect_Frame, 
	2	Sweeps_per_IFg, Points, Gmt, Numrec, 
	3       Time_Flag, Rpt_Lun)
	    Check = .True.
	    Numrec = 1
	    Fault = .False.
	    Time_Flag = 'EB'
	    Prtnxt = .True.
	  Else
	    If ( Prtnxt ) Then	! Log end of block
 	      Status = FTB_Dump_Sci_Time ( Channel, Collect_Gmt, Status_Bits, 
	1	Absum, Speed, Length, Transmit_Frame, Collect_Frame, 
	2	Sweeps_per_Ifg, Points, Gmt, Numrec, 
	3	Time_Flag, Rpt_Lun)
	      Check = .True.   
	      Prtnxt = .False.
             Write( Unit=Rpt_lun, Fmt=100, Iostat=Rstatus) Dash
 100         Format(1x, 130a)
	    Endif
	    Numrec = 1
	    Time_Flag = 'BB'
	  Endif
	Elseif ( Time_LT ( Time2, Time1 )) Then ! Dump start record regardless
	                                        ! of how flags are set, i.e.,
						! First and Fault
	        If (Eof) Time_Flag = 'EF'
	  Status = FTB_Dump_Sci_Time ( Channel, Collect_Gmt, Status_Bits, 
	1	Absum, Speed, Length, Transmit_Frame, Collect_Frame, 
	2	Sweeps_per_IFg, Points, Gmt, Numrec,
	3       Time_Flag, Rpt_Lun)
	  Check = .True.   
	  Time_Flag = 'LT'
	  Fault = .True.
	  Numrec = 1
	    
	Else					! Times are equal 
          If ( (EOF) .OR. (Time_flag .EQ.'BB') ) Then ! If new block or Eof, 
				                      ! dump it
	      If (Eof) Time_Flag = 'EF'           ! End of file only
	    Status = FTB_Dump_Sci_Time ( Channel, Collect_Gmt, Status_Bits, 
	1	Absum, Speed, Length, Transmit_Frame, Collect_Frame, 
	2	Sweeps_per_IFg, Points, Gmt, Numrec,
	3       Time_Flag, Rpt_Lun)
	    Check = .True.   
	  Endif
	  Numrec = Numrec + 1
	  Time_Flag = 'EQ'
	  Fault = .True.
	Endif

	If ((Check) .AND. (Status .NE. %loc(FTB_Normal))) Then	    
	  FTB_Check_Sci_Time = %loc(FTB_Aberr)
	Endif

	If (.Not. First) Then
	  Sci_Rec(1) = Sci_Rec(2)
	  Sci_Rec(2) = Sci_Rec(3)
	Endif
c
c *** add codes here to write out gmt,eof,time_flag,fault,time_gt,time_lt
c *** QCC/STX 12/22/88
c
        If ( Dump ) Then
          Nrec = Nrec + 1
          Print *,'rec#: ',nrec,' gmt: ',gmt,' eof: ',eof,
     1         ' time_flag: ',time_flag, ' fault: ',fault,
     2         ' time_gt: ',time_gt(time2,time1),
     3         ' time_lt: ',time_lt(time2,time1), 
     4         ' time_converted: ', contime,iflag
        Endif
	Return
	End
