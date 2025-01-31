
C-------------------------------------------------------------------------------

	Integer*4 Function FTB_Check_Hkp_Time 
	1	  ( Hkp_Rec, First, Eof, Rpt_Lun , Dump )

C-------------------------------------------------------------------------------
C
C	Purpose: To check the time sequence for NFS housekeeping records and
C	         log all records with repeating time keys or time keys which
C		 occur in reverse chronological order.
C
C	Author: Shirley M. Read
C		STX, November, 1988
C
C	Invocation: Status = FTB_Check_Hkp_Time 
C			     ( Channel, Hkp_Rec, First, Eof, Rpt_Lun, Dump )
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
C	  Hkp_Rec(3)    NFS_HKP record  FIRAS Housewkeeping records
C	  First         L*1             First set of records for run
C	  Eof           L*1             Last set of records for run	   
C	  Rpt_Lun       I*4             Report unit number
C         Dump          L*1             Dump Record On/Off Switch
C
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	
C	Subroutines Called:
C
C	  FTB_Dump_Hkp_Time  -- Print routine for FIRAS housekeeping record times
C         Time_LT	     -- COBETRIEVE time less than function
C         Time_GT            -- COBETRIEVE time greater than function
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
C	  Set the start and stop pointers for the housekeeping records in core
C            according to whether it is the first call or not.
C         Extract the primary time key from the COBETRIEVE header for the
C	     start and stop records.
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

	Dictionary 'NFS_HKP'         
	Record / NFS_HKP / Hkp_Rec(3)   ! FIRAS housekeeping records
	Logical*1 First, Eof            ! First and last calls
	Integer*4 Rpt_Lun               ! Report unit number
        Logical*1 Dump                  ! Dump record On/Off Switch

! 	Functions.

	Logical*1 Time_LT		! CT Time less than
	Logical*1 Time_GT               ! CT Time grater then
	Integer*4 FTB_Dump_Hkp_Time

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
	Character*14 Gmt                ! Record time in CT Header.
	Logical*1 Check / .False. /
        Character*1  Dash(130)/ 130 * '_'/
        Character*14 Contime
        Logical*4    First_gmt/.True./          
        Integer*4    Nrec               ! dump record number

	FTB_Check_Hkp_Time = %loc(FTB_Normal)

	Check = .False.

	If ( First ) Then
	  Start = 1
	  Stop = 2
	Else
	  Start = 2
	  Stop = 3
	Endif

	Do Ix = 1 ,2 
	  Time1(Ix) = Hkp_Rec(Start).CT_Head.Time(Ix)
	  Time2(Ix) = Hkp_Rec(Stop).CT_Head.Time(Ix)
	Enddo	 

        call ct_binary_to_gmt(hkp_rec(start).ct_head.time,contime)
c
c **** initial the first gmt record number
c
        if ( first_gmt) then
         nrec = 0
         first_gmt = .False.
        endif

	Gmt = Hkp_Rec(Start).CT_Head.Gmt

	If ( Time_GT( Time2, Time1 )) Then	! Time OK on current records 
	  If (Eof) Time_Flag = 'EF'  		! End of file only
	  If ( Fault .OR. Eof ) Then            ! Previous time seq. bad or Eof
	    Status = FTB_Dump_Hkp_Time (Gmt, Numrec, Time_Flag, Rpt_Lun)
	    Check = .True.
	    Numrec = 1
	    Fault = .False.
	    Time_Flag = 'EB'
	    Prtnxt = .True.
	  Else
	    If ( Prtnxt ) Then			! Log end of block
	      Status = FTB_Dump_Hkp_Time (Gmt, Numrec, Time_Flag, Rpt_Lun)
	      Check = .True.
	      Prtnxt = .False.
              Write(Unit=Rpt_lun, Fmt=100, Iostat=Rstatus) Dash
 100          Format(1x,130a)
	    Endif
	    Numrec = 1
	    Time_Flag = 'BB'
	  Endif
	Elseif ( Time_LT ( Time2, Time1 )) Then ! Dump start record regardless
	                                        ! of how flags are set, i.e.,
                                                ! First and Fault
	  If (Eof) Time_Flag = 'EF'  		! End of file only
	  Status = FTB_Dump_Hkp_Time (Gmt, Numrec, Time_Flag, Rpt_Lun)	    
	  Check = .True.
	  Time_Flag = 'LT'
	  Fault = .True.
	  Numrec = 1
	    
	Else					! Times are equal 
	  If ( (Eof) .OR. (Time_Flag .EQ. 'BB')) Then ! If last record, dump it
	    If (Eof) Time_Flag = 'EF'  		! End of file only
	    Status = FTB_Dump_Hkp_Time (Gmt, Numrec, Time_Flag, Rpt_Lun)    
	    Check = .True.
	  Endif
	  Numrec = Numrec + 1
	  Time_Flag = 'EQ'
	  Fault = .True.
	Endif

	If ((Check) .AND. (Status .NE. %loc(FTB_Normal))) Then	    
	  FTB_Check_Hkp_Time = %loc(FTB_Aberr)
	Endif
	If (.Not. First) Then
	  Hkp_Rec(1) = Hkp_Rec(2)
	  Hkp_Rec(2) = Hkp_Rec(3)
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
     4         ' time_converted: ', contime
        Endif
	Return
	End
