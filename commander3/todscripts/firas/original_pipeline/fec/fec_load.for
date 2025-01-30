      integer*4 function FEC_Load (CT_Unit,
     .                             XCal_Pos,
     .                             Query_Flag,
     .                             RSE_Flag,
     .                             Fakeit_Flag,
     .                             Data_Array,
     .                             Data_Array_Size,
     .                             Plat_Array,
     .                             Plat_Array_Size)
c
c     By Ken Jensen, STX, 25-FEB-1991
c
c     Revision of FEC_LOAD_ARRAY_FROM_ARCHIVE.FOR
c     By Reid Wilson, SASC Technologies Inc,  8-AUG-1986
c
c     Revision made to distinguish final requirements version of FEC
c     from older versions of FEC.
c
c     Description:
c         FEC_Load uses the COBETRIEVE utilities
c         to read information from an archive file into an array of
c         engineering records.  The records in which certain variables
c         (see routine FEC_Command_Change for the particulars) change
c         are noted and the index saved.
c     
c     Format of call:
c         Ret_Code = FEC_Load ( CT_Unit,
c                               XCal_Pos,
c                               Query_Flag,
c                               RSE_Flag,
c                               Fakeit_Flag,
c                               Data_Array,
c                               Data_Array_Size,
c                               Plat_Array,
c                               Plat_Array_Size)
c
c     Input Parameters:
c
c         CT_Unit            = logical unit connected to archive file
c                              integer*2
c	  XCal_Pos	     = acceptable XCal positions
c
c         Query_Flag        = Select by Query List
c 
c         RSE_Flag        = Select by RSE List
c 
c         Fakeit_Flag        = Select FAKEIT or No FAKEIT
c 
c     Output Parameters:
c
c         Data_Array         = array of engineering records
c                              RECORD/FDQ_ENG_DATA/ Data_Array(*)
c         Data_Array_Size    = the number of engineering records actually read
c                              integer*4 Data_Array_Size
c         Plat_Array         = array containing the indices of engineering
c                              records in Data_Array that precede a change in
c                              selected variables (see FEC_Command_Change for
c                              more information).  This corresponds to the
c                              rightmost point of the plateau.
c                              integer*4 Plat_Array(0:MAX_PLATEAUS)
c         Plat_Array_Size    = the number of plateau points found
c                              integer*4 Plat_Array_Size
c         Ret_Code           = status of operation
c                              integer*4 Ret_Code
c
c     Modules called:
c         FEC_Command_Change
c         FEC_Hot_Spot_Command_Change
c         FEC_Hot_Flag_Plateau
c         FEC_Sort
c
c     Include files:
c         CT$Library:CTUSER.INC      ! COBETRIEVE routines
c         FEC_INC                    ! FEC constants
c         FEC_MSG                    ! FEC message definition
c
c-------------------------------------------------------------------------
c Changes
c	SPR 2531, R. Kummerer, User specified acceptable XCal positons.
c	October 11, 1988.
c
c	SER 4210, R. Kummerer, Prevent overlaps in raw science related
c	changes. Track with most recent segment catalog/IFG matrix.
c	August 24, 1989.
c
c       SER 6851, H. Wang, improvements to fec
c                 add the capability to check hot spot command
c                 June 13, 1990, STX. 
c
c       SPR 9517, K. Jensen, STX, Properly handle hot spot plateaus
c                 when A side and B side plateaus overlap.
c                 February 13, 1992.
c
c-------------------------------------------------------------------------
c
      implicit none
c
c     Include files:
c
      include 'CT$Library:CTUSER.INC'    !contains COBETRIEVE information
      include '(FEC_INC)'                !contains constants
      include '(FEC_MSG)'                !contains FEC messages
      include '(FUT_PARAMS)'             !contains XCAL in values
      include '(CCT_FileSpec_Fields_Record)'  !contains filespec information
c
c     Formal Parameters:
c
      integer * 4                 CT_Unit
      integer * 4                 Ret_Code
      DICTIONARY 'FDQ_ENG'  
      RECORD /FDQ_ENG/            Data_Array(*)
      integer * 4                 Data_Array_Size
      integer * 4                 Plat_Array(0:MAX_PLATEAUS)
      integer * 4                 Plat_Array_Size
      integer * 4                 Plat_Hot_a_Command_Array_Size
      integer * 4                 Plat_Hot_a_Command_start_Array(50)
      integer * 4                 Plat_Hot_a_Command_end_Array(50)
      integer * 4                 Plat_Hot_b_Command_start_Array(50)
      integer * 4                 Plat_Hot_b_Command_end_Array(50)
      integer * 4                 Plat_Hot_b_Command_Array_Size
      integer * 4                 Plat_Hot_Command_start_Array(50)
      integer * 4                 Plat_Hot_Command_end_Array(50)
      integer * 4                 Plat_Hot_Command_Array_Size
      integer * 4                 Start_Array_Size
      integer * 4                 End_Array_Size
      integer * 4                 Start_Index_A(50)
      integer * 4                 Start_Index_B(50)
      integer * 4                 End_Index_A(50)
      integer * 4                 End_Index_B(50)
      integer * 4                 StartA
      integer * 4                 StartB
      integer * 4                 EndA
      integer * 4                 EndB
      integer*2			  XCal_Pos
      integer*4                   Query_Flag
      integer*4                   RSE_Flag
      integer*4                   Fakeit_Flag
      integer*4                   status
      integer*4                   i,k,j
      Logical*1                   Hot_Spot_a_Flag
      Logical*1                   Hot_Spot_b_Flag
      Logical*1                   Hot_Spot_start(2)
      Logical*1                   Hot_Spot_end(2)
c
c     Called routines:
c
      integer*4                   FEC_Command_Change
      integer*4                   FEC_Hot_Spot_Command_Change
      integer*4                   FEC_hot_flag_Plateau
      integer*4                   FEC_Sort
c
c     Local Variables:
c
      integer*4                   side
      logical*1                   EOF/.FALSE./
      integer*2                   CT_Status(20)
c
c     ** BEGIN FEC_Load **
c
c     Initialize counters to 0
c
      Data_Array_Size = 0
      Plat_Array_Size = 0
      Plat_Hot_a_command_Array_Size = 0
      Hot_Spot_a_flag = .false.
      Plat_Hot_b_command_Array_Size = 0
      Hot_Spot_b_flag = .false.
      hot_spot_start(1)=.false.
      hot_spot_end(1) = .false.
      hot_spot_start(2)=.false.
      hot_spot_end(2) = .false.
      Plat_Hot_command_Array_Size = 0
      Start_Array_Size = 0
      End_Array_Size = 0
c
c     Repeat loop until EOF on archive reached.  Each iteration read
c     a complete record into the next spot in the array.  Whenever the
c     change command occurs, the index of the previous record is pushed
c     onto the array Plat_Array (it contains the index of the rightmost 
c     point of each plateau.)
c
      do while ((.not. EOF) .and. (Data_Array_Size .lt. MAX_DATA)) !loop until EOF
          Data_Array_Size = Data_Array_Size + 1  !Push new record onto array,
          if( Data_Array_Size .gt. MAX_DATA) then
              FEC_Load = %loc(FEC_DATAOVERFLOW)
              call lib$signal (FEC_DATAOVERFLOW)
          else
              if( query_flag .eq. 1 .or. rse_flag .eq. 1 )then
                call CT_Query_Get (,CT_Unit,Data_Array(Data_Array_Size),CT_Status)
              else
                call CT_Read_Arcv (,CT_Unit,Data_Array(Data_Array_Size),CT_Status)
              endif
              if(CT_Status(1) .eq. CTP_ENDOFFILE) then
                  EOF = .TRUE.                          !Set to drop out of loop
                  Data_Array_Size = Data_Array_Size - 1 !Decrement to correct value
                  Plat_Array_Size = Plat_Array_Size + 1
                  Plat_Array(Plat_Array_Size) = Data_Array_Size
	      else if (CT_Status(1) .ne. CTP_NORMAL) then
		  Call LIB$Signal(FEC_CTQGETERROR,%Val(1),%Val(ct_status(1)))
		  FEC_Load = %Loc(FEC_CTQGETERROR)
		  Return
              else if(Fakeit_Flag .Eq. 0 .And.
     .                (Data_Array(Data_Array_Size).Chan(1).Fakeit .ne. 0 .or.
     .                 Data_Array(Data_Array_Size).Chan(2).Fakeit .ne. 0 .or.
     .                 Data_Array(Data_Array_Size).Chan(3).Fakeit .ne. 0 .or.
     .                 Data_Array(Data_Array_Size).Chan(4).Fakeit .ne. 0 ))then
                   !if /NOFAKEIT and Fakeit Data then
                   !do not add to the array
                  Data_Array_Size = Data_Array_Size - 1
              else if(Fakeit_Flag .Eq. 1 .And.
     .                (Data_Array(Data_Array_Size).Chan(1).Fakeit .ne. 1 .or.
     .                 Data_Array(Data_Array_Size).Chan(2).Fakeit .ne. 1 .or.
     .                 Data_Array(Data_Array_Size).Chan(3).Fakeit .ne. 1 .or.
     .                 Data_Array(Data_Array_Size).Chan(4).Fakeit .ne. 1 ))then
                   !if /FAKEIT and Not Fakeit Data then
                   !do not add to the array
                  Data_Array_Size = Data_Array_Size - 1
              else if(XCal_Pos .Eq. fac_xcalin .And.
     .		      (Data_Array(Data_Array_Size).En_XCal.Pos(1) .ne. fac_xcalin .or.
     .                 Data_Array(Data_Array_Size).En_XCal.Pos(2) .ne. fac_xcalin)) then
                   !if the XCAL ("POS") is not IN for side A or B then
                   !do not add to the array
                  Data_Array_Size = Data_Array_Size - 1
              else if(XCal_Pos .Eq. fac_xcalout .And.
     .		      (Data_Array(Data_Array_Size).En_XCal.Pos(1) .ne. fac_xcalout .or.
     .                 Data_Array(Data_Array_Size).En_XCal.Pos(2) .ne. fac_xcalout)) then
                   !if the XCAL ("POS") is not OUT for side A or B then
                   !do not add to the array
                  Data_Array_Size = Data_Array_Size - 1
              else if(XCal_Pos .Eq. fac_both .And.
     .		      ((Data_Array(Data_Array_Size).En_XCal.Pos(1) .ne. fac_xcalin .and.
     .		        Data_Array(Data_Array_Size).En_XCal.Pos(1) .ne. fac_xcalout) .or.
     .                 (Data_Array(Data_Array_Size).En_XCal.Pos(2) .ne. fac_xcalin .and.
     .                  Data_Array(Data_Array_Size).En_XCal.Pos(2) .ne. fac_xcalout))) then
                   !if the XCAL ("POS") is not IN or OUT for side A or B then
                   !do not add to the array
                  Data_Array_Size = Data_Array_Size - 1
             else if ((data_array(data_array_size).en_head.dq_summary_flag(1)
     1    	       .eq. 127) .or.
     1   	       (data_array(data_array_size).en_head.dq_summary_flag(2)
     1   	       .eq. 127) .or.
     1  	       (data_array(data_array_size).en_head.dq_summary_flag(3)
     1  	       .eq. 127) .or.
     1  	       (data_array(data_array_size).en_head.dq_summary_flag(4)
     1  	       .eq. 127)) then
                   data_array_size = data_array_size  - 1
              else 
	        status=FEC_Command_Change(Data_Array_Size,Data_Array)
                  !if change command has been given, then previous point is the
                  !rightmost point of the previous plateau.  Add index to array
                  !that contains indices of plateau points.
                if(status) then
                  Plat_Array_Size = Plat_Array_Size + 1
                  if( Plat_Array_Size .gt. MAX_PLATEAUS) then
                      FEC_Load = %loc(FEC_PLATOVERFLOW)
                      call lib$signal(FEC_PLATOVERFLOW)
                  else
                    Plat_Array(Plat_Array_Size) = 
     1  	    Data_Array_Size - 1
                  end if
                endif
              
C
C        call fec_hot_spot_command to find hot spot a side plateau
c
               side = 1
	       status=FEC_Hot_Spot_Command_Change(side,hot_spot_start,
     1  	        hot_spot_end,Data_Array_Size,Data_Array)
                  !if change command has been given, then previous point is the
                  !rightmost point of the previous plateau.  Add index to array
                  !that contains indices of plateau points.
               if(status) then
                 if (hot_spot_start(side)) then
                  Plat_Hot_a_command_Array_Size = Plat_Hot_a_command_Array_Size+1
                  if( Plat_Hot_a_command_Array_Size .gt. 50) then
                      FEC_Load = %loc(FEC_PLATOVERFLOW)
                      call lib$signal(FEC_PLATOVERFLOW)
                  else
                   hot_spot_start(side) = .false.
                   Hot_Spot_a_Flag = .true.
                   Plat_hot_a_command_start_Array(Plat_Hot_a_command_Array_Size) = 
     1   	    Data_Array_Size - 1 
                  end if
                 endif
                 if (hot_spot_end(side)) then
                   Hot_Spot_a_Flag = .true.
                   Plat_hot_a_command_end_Array(Plat_Hot_a_command_Array_Size) = 
     1  	      Data_Array_Size - 1
                   hot_spot_end(side)=.false.
                 endif
               endif
               side = 2
	       status=FEC_Hot_Spot_Command_Change(side,hot_spot_start,
     1  	               hot_spot_end,Data_Array_Size,
     1  	                                  Data_Array)
                  !if change command has been given, then previous point is the
                  !rightmost point of the previous plateau.  Add index to array
                  !that contains indices of plateau points.
                if(status) then
                 if (hot_spot_start(side)) then
                 Plat_Hot_b_command_Array_Size = Plat_Hot_b_command_Array_Size+1 
                 if( Plat_Hot_b_command_Array_Size .gt. 50) then
                      FEC_Load = %loc(FEC_PLATOVERFLOW)
                      call lib$signal(FEC_PLATOVERFLOW)
                 else
                   hot_spot_start(side) = .false.
                   Hot_Spot_b_Flag = .true.
                   Plat_hot_b_command_start_Array(Plat_Hot_b_command_Array_Size) = 
     1                 Data_Array_Size - 1
                 end if
                endif
                if (hot_spot_end(side)) then
                 Hot_Spot_b_Flag = .true.
                 Plat_hot_b_command_end_Array(Plat_Hot_b_command_Array_Size) = 
     1   	    Data_Array_Size - 1
                 hot_spot_end(side)=.false.
                endif
              endif
          end if
        endif
      end do
c
c  if hot spot plateau exist then
c    call fec_hot_flag_plateau to  flag hot spot into  plateau array
c  endif   
c
      If (Hot_Spot_a_Flag) Then
        if (plat_hot_a_command_array_size .eq. 0) hot_spot_a_flag = .false.
      endif           
      If (Hot_Spot_b_Flag) Then
        if (plat_hot_b_command_array_size .eq. 0) hot_spot_b_flag = .false.
      endif 
      if ( (hot_spot_a_flag) .and. (.not. hot_spot_b_flag) ) then
        do i = 1, plat_hot_a_command_array_size
          plat_hot_command_start_array(i) = plat_hot_a_command_start_array(i)
          plat_hot_command_end_array(i) = plat_hot_a_command_end_array(i)
        enddo
        start_array_size = plat_hot_a_command_array_size
        end_array_size = start_array_size
      endif
      if ( (hot_spot_b_flag) .and. (.not. hot_spot_a_flag) ) then
        do i = 1, plat_hot_b_command_array_size
          plat_hot_command_start_array(i) = plat_hot_b_command_start_array(i)
          plat_hot_command_end_array(i) = plat_hot_b_command_end_array(i)
        enddo
        start_array_size = plat_hot_b_command_array_size
        end_array_size = start_array_size
      endif
      if ( (hot_spot_b_flag) .and. (hot_spot_a_flag) ) then
c
c      Eliminate A side indices which are overlapped by B side plateaus
c
        do i = 1, plat_hot_a_command_array_size
          starta = plat_hot_a_command_start_array(i)          
          enda = plat_hot_a_command_end_array(i)
          start_index_a(i) = 1          
          end_index_a(i) = 1          
          do j=1,plat_hot_b_command_array_size
            startb = plat_hot_b_command_start_array(j)
            endb = plat_hot_b_command_end_array(j)
            if((startb .lt. starta).and.(endb .ge. starta))start_index_a(i) = 0
            if((startb .le. enda).and.(endb .gt. enda))end_index_a(i) = 0
          enddo
        enddo
c
c      Eliminate B side indices which are overlapped by A side plateaus
c
        do i = 1, plat_hot_b_command_array_size
          startb = plat_hot_b_command_start_array(i)          
          endb = plat_hot_b_command_end_array(i)
          start_index_b(i) = 1          
          end_index_b(i) = 1          
          do j=1,plat_hot_a_command_array_size
            starta = plat_hot_a_command_start_array(j)
            enda = plat_hot_a_command_end_array(j)
            if((starta .le. startb).and.(enda .ge. startb))start_index_b(i) = 0
            if((starta .le. endb).and.(enda .ge. endb))end_index_b(i) = 0
          enddo
        enddo
c
        start_array_size = 0
        end_array_size = 0
c
        do i=1,plat_hot_a_command_array_size
          if (start_index_a(i) .eq. 1) then
            start_array_size = start_array_size + 1
            plat_hot_command_start_array(start_array_size)=plat_hot_a_command_start_array(i)
          endif
          if (end_index_a(i) .eq. 1) then
            end_array_size = end_array_size + 1
            plat_hot_command_end_array(end_array_size)=plat_hot_a_command_end_array(i)
          endif
        enddo
c
        do i=1,plat_hot_b_command_array_size
          if (start_index_b(i) .eq. 1) then
            start_array_size = start_array_size + 1
            plat_hot_command_start_array(start_array_size)=plat_hot_b_command_start_array(i)
          endif
          if (end_index_b(i) .eq. 1) then
            end_array_size = end_array_size + 1
            plat_hot_command_end_array(end_array_size)=plat_hot_b_command_end_array(i)
          endif
        enddo
c
        If ( start_array_size .ne. end_array_size ) Then
          FEC_Load = %loc(FEC_HOTMERGE)
          call lib$signal(FEC_HOTMERGE)
          return
        Endif
        Ret_Code = FEC_Sort (plat_hot_command_start_array,start_array_size)
        If (.not. Ret_Code) Then
          call lib$signal(%val(Ret_Code))
          FEC_Load = Ret_Code
          return
        Endif
        Ret_Code = FEC_Sort (plat_hot_command_end_array,end_array_size)
        If (.not. Ret_Code) Then
          call lib$signal(%val(Ret_Code))
          FEC_Load = Ret_Code
          return
        Endif
c
      endif
c
      if ( start_array_size .gt. 0 ) then
       do i=1,start_array_size
         Status=Fec_hot_flag_plateau(plat_hot_command_start_array(i),
     1                     plat_hot_command_end_array(i),
     1                     plat_array_size, plat_array)
       enddo
      endif
c
      FEC_Load = %loc(FEC_NORMAL)
      return 
      end !FEC_Load
