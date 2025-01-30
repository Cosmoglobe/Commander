	Integer*4 Function FEP_Fix_Counts_to_Ohms (Fr_lun,cal_lun,
	1	  co_coeffs )

C------------------------------------------------------------------------
C    PURPOSE: Scans the housekeeping data and averages the counts from
C	      each of the calibration resistors and computes the standard
C	      deviation as well.  Rescans, throwing out all of the counts
C	      that differ by more than three standard deviations, and
C	      recomputes the averages.  Finally, computes the coefficients
C	      required to convert counts to ohms for each spacecraft side and
C	      current range.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Don Stevens-Rayburn
C            Applied Research Corporation
C            July 13, 1989
C
C    INVOCATION: STATUS = FEP_Fix_Counts_to_Ohms ( FR_LUN, Cal_Lun,co_coeffs )
C
C    INPUT PARAMETERS:
C	Fr_Lun			I*4		Fortran unit number
C       Cal_lun                 I*4             fortran unit number for 
C                                               reference file
C
C    OUTPUT PARAMETERS:
C	co_coeffs ( 3, 4 )	R*4		Coefficients required to
C						convert counts to ohms.
C
C    SUBROUTINES CALLED:
C
C	FEP_Get_Calres
C	DRline
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: FEP_Data
C
C----------------------------------------------------------------------
C    Changes:
C
C	SPR 4463.  Engplots aborts:  divide by zero in FEP_Fix_Counts_to_Ohms. 
C	When all values are equal, their sample std dev will be 0.  Test for
C	non-outliers was "x .gt. avg-3s .and. x .lt. avg+3s", so when s=0,
C	nothing passed.  Changed .gt. and .lt. to .ge. and .le.  Also put in a
C	check for n=0 or 1, so as not to divide by 0 when computing avg and s. 
C	Normally, with lots of data points, we needn't worry about s=0, but when
C	the chosen timerange has only two points, it becomes pretty likely. 
C	Qfix 422. Fred Shuman,  STX,  1989 Aug 31. 
C
C       SPR 4852, add the capability to handle the dwell data right.
C       Qfix #, Harte Wang, STX, 1989 Nov. 3
C
C       SPR 5068, avoid plotting data points when telm. qual. is bad
C       Qfix #, Harte Wang, STX, 1989 DEC. 2
C
C       SPR 5540, Confusion over DWELL MODE.
C                 Harte Wang, STX, 1990 Jan. 11.
C
C       SPR 3921, Add the capabilities to read reference data from the
C                 reference archive.   
C                 Harte Wang, STX, 1990 Feb. 9.
C
C----------------------------------------------------------------------
C PDL:
C  Write  String :'   Computing Counts to Ohms Conversion Coefficients'
C    to Screen
C  Initialized local variables
C  Do for each record
C   Do for each side/current (AL,AH,BL,BH)
C       Get the sum of calibration resistors 
C       Get the sum of square of calibration resistors
C      Enddo
C    ENDDO
C  ENDDO    
C   Get the average of calibration resistor values
C   Compute the standard deviation  
C   call Fut_Ave_Bin_Times to get the average time 
C   If return status for the function is bad
C   Then
C     Set return status for the function to bad
C     Write error message to screen
C     Return
C   Endif
C   Call CCT_Get_Config_Idx_Tod to get the index and new_reference_data_flag 
C   If return status is bad
C   Then
C     Set return status to bad
C     Write error message to screen
C     Return
C   Endif
C   If new_reference_data_flag is true
C   Then
C    Call Fep_Get_Calres to get 
C          the conversion coefficients of resistor calibration.
C   Endif
C   If return status is bad
C   Then
C     set return status to bad
C     Write error message to screen
C     Return
C   Endif
C   Do for each Side/current
C    compute coefficients required to convert counts to ohms for each 
C      spacecraft side and current range.
C   Enddo
C   Set the return status to good
C   Return
C   End 

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FUT_Invoc)'
	Include		'(FUT_Params)'
	Include		'(FEP_Data)'
	Include		'(CCT_GET_CONFIG)'
        Record /config_status/ stat1
        Integer         *4      Bin_Time(2,2)
        Integer         *4      Ave_time(2)
        Integer         *4      Weight(2)  
	Integer		*4	caloff ( 4 ) / 56, 88, 120, 152 /
	Integer		*4	cal_counts ( 4, 4 )
	Integer		*4	cal_used   ( 4, 4 )
	Integer		*4	iside
        Integer         *4      index
	Integer		*4	i
	Integer		*4	j
        Integer         *4      fr_lun
        Character       *64     dwell_field
        Character       *64     Telm_field
        Integer         *2      dwell_offset_a
        Byte                    TELM_stat
        Integer         *2      telm_offset
        Integer         *2      dwell_offset_b
        Integer         *2      TELM_length
        Integer         *2      dwell_length
        Byte                    dwell_stat_a
        byte                    dwell_stat_b
        Byte                    dwell_byte
        integer         *2      equiv,dwecont/0/ 
        Byte                    null(2)/0,0/
        equivalence (null(1),equiv )        
	Integer		*4	k
	Integer		*4	kbasis
	Integer		*4	kal_ok
        Integer         *4      Cal_lun
        Logical         *1      New_Cal
	Integer		*4	status

	Real		*4	co_coeffs ( 3, 4 )
	Real		*8	tmp_coeffs ( 3 )
	Real		*8	sum   ( 4, 4 )
	Real		*8	sumsq ( 4, 4 )
	Real		*8	aver  ( 4, 4 )
	Real		*8	three_sigma ( 4, 4 )
	Real		*8	q ( 4, 3 )
	Real		*8	rhs ( 4 )
	Real		*8	tlm
	Real		*8	tol
	Real		*8	res ( 4 )

	Real		*4	resis_cal ( 4, 4 )

	Integer		*4	FEP_Get_Calres
	Integer		*4	FUT_field_ATTRIBUTES
        Integer         *4      CCT_Get_Config_Idx_Tod
        Integer         *4      Fut_Ave_Bin_Times
	External		FEP_Normal
	External		FUT_Normal
        EXternal                FEP_GETCONFIGERR
C-----------------------------------------------------------------------------
C         CODE Start
C-----------------------------------------------------------------------------
	Type *, '   Computing Counts to Ohms Conversion Coefficients'

C
C Get the dwell status mode offsets, side a and side b.
C	

        TELM_field='FEP_HKP.CURR_MENU.TELEMETRY_QUALITY'
        STATUS = FUT_FIELD_Attributes(FR_LUN,TELM_field,TELM_length,
	1        TELM_offset)
        If (status .Ne. %loc(FUT_Normal)) Then
          Fep_Fix_counts_to_Ohms = Status
          Return
        Endif
        Dwell_field='FEP_HKP.HSKP_MENU.DWELL_STATUS_A'
        STATUS = FUT_FIELD_Attributes(FR_LUN,dwell_field,dwell_length,
	1        dwell_offset_a)

        If (status .Ne. %loc(FUT_Normal)) Then
          Fep_Fix_counts_to_Ohms = Status
          Return
        Endif
C
        Dwell_field='FEP_HKP.HSKP_MENU.DWELL_STATUS_B'
        STATUS = FUT_FIELD_Attributes(FR_LUN,dwell_field,dwell_length,
	1        dwell_offset_b)

        If (status .Ne. %loc(FUT_Normal)) Then
          Fep_Fix_counts_to_Ohms = Status
          Return
        Endif



        Do i = 1,4
	   Do j = 1,4
	      sum      ( i, j ) = 0.0
	      sumsq    ( i, j ) = 0.0
	      cal_used ( i, j ) = 0
	   End Do
	End Do

	Do i = 1, nrecs
    	   dwell_stat_a=0
           dwell_stat_b=0
           telm_stat = 0 
           If (hskp_data(telm_offset+1,i) .ne. 0) telm_stat = 1
           Dwell_byte = hskp_data(dwell_offset_a+1,i)
           null(1) = dwell_byte
           If ( dwell_byte .lt. 0) then
            call MVbits(equiv,0,5,dwecont,0)
            if (dwecont .lt. 31) dwell_stat_a = 1
           Endif
           Dwell_byte = hskp_data(dwell_offset_b+1,i)
           null(1) = dwell_byte
           If ( dwell_byte .lt. 0) then
            call MVbits(equiv,0,5,dwecont,0)
            if (dwecont .lt. 31) dwell_stat_b = 1
           Endif
           IF (telm_stat .eq. 0) then
           
            Do iside = 1, 4
              IF (((iside .eq. 1 .or. iside .eq. 2) .and. dwell_stat_a .eq. 0)
	1        .or. ((iside .eq. 3 .or. iside .eq. 4) .and. dwell_stat_b .eq. 
	2         0) ) then               
              Do j = 1, 4
	         k = 0
		 Call LIB$Movc3 ( 2, hskp_data( caloff( iside )+j*2-1,i), k )
	         If ( ( k .gt. 100 ) .and. ( k .lt. 16000 ) ) Then
	            sum   ( j, iside ) = sum   ( j, iside ) + dfloat ( k )
	            sumsq ( j, iside ) = sumsq ( j, iside ) + dfloat ( k * k )
	            cal_used ( j, iside ) = cal_used ( j, iside ) + 1
	         End If
	      End Do
             Endif
	   End Do
         endif
	End Do

	Do i = 1, 4
	   Do j = 1, 4
	      If ( cal_used ( i, j ) .eq. 0 ) Then
	         aver  ( i, j ) = 0.0
	         three_sigma ( i, j ) = 0.0
	      Else If ( cal_used ( i, j ) .eq. 1 ) Then
	         aver  ( i, j ) = sum ( i, j )
	         three_sigma ( i, j ) = 0.0
	      Else
	         aver  ( i, j ) = sum ( i, j ) / float ( cal_used ( i, j ) )
	         three_sigma ( i, j ) = 3.0 * sqrt ( ( sumsq ( i, j )
	1                                   - sum ( i, j )**2
	2                                   / float ( cal_used ( i, j ) ) )
	3                                   / ( cal_used ( i, j ) - 1 ) )
	      End If
	   End Do
	End Do

	Do i = 1,4
	   Do j = 1,4
	      sum      ( i, j ) = 0.0
	      cal_used ( i, j ) = 0
	   End Do
	End Do

	Do i = 1, nrecs
    	   dwell_stat_a=0
           dwell_stat_b=0
           Telm_stat=0
           If (hskp_data(telm_offset+1,i) .ne. 0) telm_stat = 1
           Dwell_byte = hskp_data(dwell_offset_a+1,i)
           null(1) = dwell_byte
           If ( dwell_byte .lt. 0) then
            call MVbits(equiv,0,5,dwecont,0)
            if (dwecont .lt. 31) dwell_stat_a = 1
           Endif
           Dwell_byte = hskp_data(dwell_offset_b+1,i)
           null(1) = dwell_byte
           If ( dwell_byte .lt. 0) then
            call MVbits(equiv,0,5,dwecont,0)
            if (dwecont .lt. 31) dwell_stat_b = 1
           Endif
           if (telm_stat .eq. 0) then
            Do iside = 1, 4
              IF (((iside .eq. 1 .or. iside .eq. 2) .and. dwell_stat_a .eq. 0)
	1        .or. ((iside .eq. 3 .or. iside .eq. 4) .and. dwell_stat_b .eq. 
	2         0)  ) then               
               Do j = 1, 4
	         k = 0
		 Call LIB$Movc3 ( 2, hskp_data( caloff( iside )+j*2-1,i), k )
	         If ( ( k .ge. aver ( j, iside )
	1                    - three_sigma ( j, iside ) ) .and.
	2             ( k .le. aver ( j, iside )
	3                    + three_sigma ( j, iside ) ) ) Then
	            sum   ( j, iside ) = sum   ( j, iside ) + dfloat ( k )
	            cal_used ( j, iside ) = cal_used ( j, iside ) + 1
	         End If
	       End Do
              Endif
	   End Do
          endif
	End Do

	Do i = 1, 4
	   Do j = 1, 4
	      aver  ( i, j ) = sum ( i, j ) / float ( cal_used ( i, j ) )
	   End Do
	End Do
           Weight(1) = 1
           Bin_time(1,1) = timetags(1).time(1)
           Bin_time(2,1) = timetags(1).time(2)  
           Weight(2) = 1
           Bin_time(1,2) = timetags(nrecs).time(1)
           Bin_time(2,2) = timetags(nrecs).time(2)  
        Status = Fut_Ave_Bin_Times(Bin_time,2,Weight,Ave_time)
        If (.not. Status) Then
         Fep_Fix_Counts_to_Ohms=Status
         Call Lib$Signal(%val(Status))
         Return
        Endif
        Status = CCT_Get_Config_Idx_Tod(Ave_time,1,Cal_Lun,index,New_cal,
	1	Stat1)
        If (.not. Status) Then
         Fep_Fix_Counts_to_Ohms=Status
         Call Lib$Signal(FEP_GETCONFIGERR,%val(1),%val(Status))
         Return
        Endif
  	status = FEP_Get_Calres (cal_lun, resis_cal )
	If ( .not. status ) Then
	   Call LIB$Signal ( %val ( status ) )
	End If

	Do iside = 1, 4
	   kal_ok = 0
	   Do j = 1, 4
              tlm = aver ( j, iside )
	      If ( tlm .GT. 100 .AND. tlm .LT. 16000 ) Then
	         kal_ok = kal_ok + 1
	         q   ( kal_ok, 1 ) = tlm
	         q   ( kal_ok, 2 ) = 1
	         q   ( kal_ok, 3 ) = tlm * tlm
	         rhs ( kal_ok    ) = resis_cal ( j, iside )
	      End If
	   End Do
	   If (kal_ok .Eq. 0) Then
	      Return
	   End If

	   tol = 0.
	   Call DLsqrr ( kal_ok, 3, q, 4, rhs, tol, tmp_coeffs, res, kbasis )
	   Do j = 1, 3
	      co_coeffs ( j, iside ) = tmp_coeffs ( j )
	   End Do
	End Do

	FEP_Fix_Counts_to_Ohms = %Loc(FEP_Normal)

	Return
	End
