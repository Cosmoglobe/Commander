C-------------------------------------------------------------------------------

	Integer*4 Function FXT_Compare_Eng ( Eng_Rec, Xtrm_Rec, Filter,
	1	  Xtrm_Table, Cur_Table, Exceed, New_Xtrm, Enable_Flags)

C-------------------------------------------------------------------------------
C
C	Purpose: To compare engineering field values of the FDQ_Eng record
C	         to the current minima and maxima.
C
C	Author: Shirley M. Read
C		STX, December, 1988
C
C	Invocation: Status = FXT_Compare_Eng ( Eng_Rec, Xtrm_Rec, Filter,
C		    Xtrm_Table, Cur_Table, Exceed, New_Xtrm, Enable_Flags )
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
C	  Name		    Type		 Description
C 	  ----------------------------------------------------------------------
C	  Eng_Rec           FDQ_Eng record       FDQ_Eng archive record
C	  Xtrm_Rec          FXT_Eng_Xtrm record  FXT_Eng_Xtrm archive record
C	  Filter            I*4		         Filter value for checking
C	  Enable_Flags(118) L*1                  Flags to enable checking
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C         Xtrm_Table(118,2) L*1                  Flags indicating extrema 
C	                                         exceeded for current run
C	  Cur_Table(118)    L*1                  Flags indicating extrema
C						 exceeded for current record
C	  Exceed            L*1                  Flag indicating at least one
C					         exceeded on current record
C	  New_Xtrm          L*1                  Flag indicating at least one
C						 exceeded for current run
C	Subroutines Called:
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C	  FXT_Msg.Txt
C
C	Processing Method: PDL for FXT_Compare_Eng
C
C	  Set function status to normal.
C	  Do for each of the GRTs 
C	    If the GRT is not in Dwell Mode, Then
C	      If the FDQ_Eng value is less than the FXT_Eng_Xtrm minimum value 
C               and the Enable Flag is set, Then
C		Increment the corresponding GRT counter.
C	        If the counter = 1 or the new Eng value is less than the 
c	          old Eng value, Then
C	          Set the Eng value to the new value.
C	          Store the corresponding time.
c		Endif
C               If the counter is greater or equal to the Filter value, Then 
C		  Store the new Eng value in the FXT_Eng_Xtrm record along 
C		    with the FDQ_Eng time.
C	          Set all the corresponding extrema flags.
C	          Reset the counter.
C	        Endif
C	      Else
C	        Reset the counter.
C	      Endif
C             If the FDQ_Eng value is greater than FXT_Eng_Xtrm maximum value
C	        and the Enable flag is set, Then
C	        Perform the same sequence as for the minimum.
C	      Endif
C           Else in Dwell Mode
C	      Reset the counter.
C	    Endif Dwell Mode or not
C	  Enddo for the GRTs
C    	  Do for all Engineering Fields (Temperatures and Currents, Voltages 
C           and Currents, LMACs) in the FXT_Eng_Xtrm record:
C	    If the FDQ_Eng engineering analog value is less than the
C	      FXT_Eng_Xtrm minimum and the Enable Flag is set, Then
C	      Perform the same sequence as for the GRTs.
C	    Endif
C	    If the FDQ_Eng engineering analog value is greater than the
C	      FXT_Eng_Xtrm maximum and the Enable Flag is set, Then
C	      Perform the same sequence as for the GRTs.
C	    Endif
C	  Enddo for all engineering fields
C	  Return with normal or error status.
C
C------------------------------------------------------------------------------
C	
	Implicit	None

!	Passed Parameters.

	Dictionary 'FDQ_Eng'
	Record /FDQ_Eng/ Eng_Rec          ! FDQ_Eng archive record
	Dictionary 'FXT_Eng_Xtrm'
	Record /FXT_Eng_Xtrm/ Xtrm_Rec    ! FXT_Eng_Xtrm archive record
	Integer*4  Filter                 ! Filter value for checking
	Logical*1  Enable_Flags(118)      ! Flags to enable extrema checking
	Logical*1  Xtrm_Table(118,2)      ! Flags indicating extrema 
	                                  ! exceeded for current run
	Logical*1  Cur_Table(118)         ! Flags indicating extrema
					  ! exceeded for current record
	logical*1  Exceed                 ! Flag indicating at least one
					  ! exceeded on current record
	Logical*1  New_Xtrm               ! Flag indicating at least one
					  ! exceeded for current run
!	Include Files.

	Include    '(FXT_Msg)'

!	Local variables

	Integer*4  Status		  ! Return status variable
	Integer*4  Numex(118,2)           ! Number of consecutive records for 
					  ! which extrema are exceeded
	Real*4     Engval(118,2)          ! Exceeded engineering values
	Integer*4  Timex(2,118,2)         ! Times at which extrema are exceeded 
	Real*4     Tempval(36)            ! Temporary Eng_Rec storage
	Integer*4  One / 1 /, Zero / 0 /
	Integer*2  Min / 1 /, Max / 2 /   ! Indices for Xtrm record
	Integer*4  FV / - 9999 /          ! Integer Engineering record flag
	Integer*4  Neval                  ! Nearest Integer value for Eng field
	Integer*2  Ptr                    ! Index pointer for arrays
	Integer*2  Ix, Jx                 ! Indices
	Byte       Not_Dwell / 0 /        ! Status value for not in dwell mode

!	Set the function return status to Normal.

	FXT_Compare_Eng = %loc(FXT_Normal)

!	Check the GRTs for exceeded minima and maxima.
	
	Do Ix = 1, 64
	  If(((Eng_Rec.En_Stat.Dwell_Stat(1) .Eq. Not_Dwell) . And.
	1     (Ix .Le. 32)) .Or.
	2     ((Eng_Rec.En_Stat.Dwell_Stat(2) .Eq. Not_Dwell) . And.
	3     (Ix .Ge. 33))) Then
	    Neval = Nint(Eng_Rec.En_Analog.GRT(Ix))

!       Check the mimima.

	    If((Eng_Rec.En_Analog.GRT(Ix) .Lt. Xtrm_Rec.MinMax(Min).GRT(Ix))
	1	 .And. Enable_Flags(Ix) .And. (Neval .Ne. FV)) Then
	      Numex(Ix,Min) = Numex(Ix,Min) + 1
	      If((Numex(Ix,Min) .Eq. One) .Or. 
	1	(Eng_Rec.En_Analog.GRT(Ix) .Lt. Engval(Ix,Min))) Then
	        Engval(Ix,Min) = Eng_Rec.En_Analog.GRT(Ix)
	        Do Jx = 1, 2
	          Timex(Jx,Ix,Min) = Eng_Rec.Ct_Head.Time(Jx)
		Enddo
	      Endif
	      If(Numex(Ix,Min) .Ge. Filter) Then
		Xtrm_Rec.MinMax(Min).GRT(Ix) = Engval(Ix,Min)
		Do Jx = 1, 2
		  Xtrm_Rec.MinMax(Min).GRT_Time(Jx,Ix) = Timex(Jx,Ix,Min)
		Enddo
	        Xtrm_Table(Ix,Min) = .True.
	        Cur_Table(Ix) = .True.
		Exceed = .True.
		New_Xtrm = .True.
	        Numex(Ix,Min) = 0
	      Endif
	    Else 	! Eng_Rec value not < Xtrm_Rec or Enable_Flag not set
	      Numex(Ix,Min) = 0
	    Endif

!       Check the maxima.

	    If((Eng_Rec.En_Analog.GRT(Ix) .Gt. Xtrm_Rec.MinMax(Max).GRT(Ix))
	1	 .And. Enable_Flags(Ix) .And. (Neval .Ne. FV)) Then
	      Numex(Ix,Max) = Numex(Ix,Max) + 1
	      If((Numex(Ix,Max) .Eq. One) .Or. 
	1	(Eng_Rec.En_Analog.GRT(Ix) .Gt. Engval(Ix,Max))) Then
	        Engval(Ix,Max) = Eng_Rec.En_Analog.GRT(Ix)
	        Do Jx = 1, 2
	          Timex(Jx,Ix,Max) = Eng_Rec.Ct_Head.Time(Jx)
		Enddo
	      Endif
	      If(Numex(Ix,Max) .Ge. Filter) Then
		Xtrm_Rec.MinMax(Max).GRT(Ix) = Engval(Ix,Max)
		Do Jx = 1, 2
		  Xtrm_Rec.MinMax(Max).GRT_Time(Jx,Ix) = Timex(Jx,Ix,Max)
		Enddo
	        Xtrm_Table(Ix,Max) = .True.
	        Cur_Table(Ix) = .True.
		Exceed = .True.
		New_Xtrm = .True.
	        Numex(Ix,Max) = 0
	      Endif
	    Else 	! Eng_Rec value not > Xtrm_Rec or Enable_Flag not set
	      Numex(Ix,Max) = 0
	    Endif
	  Else		! Side of GRT in dwell mode
 	    Numex(Ix,Min) = 0
	    Numex(Ix,Max) = 0
	  Endif
	Enddo		! Ix = 1, 64 for GRTs		    

!	Check the Temperatures and Currents for exceeded minima and maxima.

	Call Lib$Movc3(64, Eng_Rec.En_Analog.IPDU_Temp(1), Tempval(1))
	Ptr = 64

	Do Ix = 1, 16
	    Neval = Nint(Tempval(Ix))
	    Ptr = Ptr + 1

!       Check the mimima.

	    If((Tempval(Ix) .Lt. Xtrm_Rec.MinMax(Min).T_and_I(Ix))
	1	 .And. Enable_Flags(Ptr) .And. (Neval .Ne. FV)) Then
	      Numex(Ptr,Min) = Numex(Ptr,Min) + 1
	      If((Numex(Ptr,Min) .Eq. One) .Or. 
	1	(Tempval(Ix) .Lt. Engval(Ptr,Min))) Then
	        Engval(Ptr,Min) = Tempval(Ix)
	        Do Jx = 1, 2
	          Timex(Jx,Ptr,Min) = Eng_Rec.Ct_Head.Time(Jx)
		Enddo
	      Endif
	      If(Numex(Ptr,Min) .Ge. Filter) Then
		Xtrm_Rec.MinMax(Min).T_and_I(Ix) = Engval(Ptr,Min)
		Do Jx = 1, 2
		  Xtrm_Rec.MinMax(Min).T_and_I_Time(Jx,Ix) = Timex(Jx,Ptr,Min)
		Enddo
	        Xtrm_Table(Ptr,Min) = .True.
	        Cur_Table(Ptr) = .True.
		Exceed = .True.
		New_Xtrm = .True.
	        Numex(Ptr,Min) = 0
	      Endif
	    Else 	! Eng_Rec value not < Xtrm_Rec or Enable_Flag not set
	      Numex(Ptr,Min) = 0
	    Endif

!       Check the maxima.

	    If((Tempval(Ix) .Gt. Xtrm_Rec.MinMax(Max).T_and_I(Ix))
	1	 .And. Enable_Flags(Ptr) .And. (Neval .Ne. FV)) Then
	      Numex(Ptr,Max) = Numex(Ptr,Max) + 1
	      If((Numex(Ptr,Max) .Eq. One) .Or. 
	1	(Tempval(Ix) .Gt. Engval(Ptr,Max))) Then
	        Engval(Ptr,Max) = Tempval(Ix)
	        Do Jx = 1, 2
	          Timex(Jx,Ptr,Max) = Eng_Rec.Ct_Head.Time(Jx)
		Enddo
	      Endif
	      If(Numex(Ptr,Max) .Ge. Filter) Then
		Xtrm_Rec.MinMax(Max).T_and_I(Ix) = Engval(Ptr,Max)
		Do Jx = 1, 2
		  Xtrm_Rec.MinMax(Max).T_and_I_Time(Jx,Ix) = Timex(Jx,Ptr,Max)
		Enddo
	        Xtrm_Table(Ptr,Max) = .True.
	        Cur_Table(Ptr) = .True.
		Exceed = .True.
		New_Xtrm = .True.
	        Numex(Ptr,Max) = 0
	      Endif
	    Else 	! Eng_Rec value not > Xtrm_Rec or Enable_Flag not set
	      Numex(Ptr,Max) = 0
	    Endif
	Enddo		! Check the T_and_I for extrema

!	Check the Voltages and Currents for exceeded minima and maxima.

	Call Lib$Movc3(144, Eng_Rec.En_Analog.bol_volt(1), Tempval(1))
	Ptr = 80

	Do Ix = 1, 36
	    Neval = Nint(Tempval(Ix))
	    Ptr = Ptr + 1

!       Check the mimima.

	    If((Tempval(Ix) .Lt. Xtrm_Rec.MinMax(Min).V_and_I(Ix))
	1	 .And. Enable_Flags(Ptr) .And. (Neval .Ne. FV)) Then
	      Numex(Ptr,Min) = Numex(Ptr,Min) + 1
	      If((Numex(Ptr,Min) .Eq. One) .Or. 
	1	(Tempval(Ix) .Lt. Engval(Ptr,Min))) Then
	        Engval(Ptr,Min) = Tempval(Ix)
	        Do Jx = 1, 2
	          Timex(Jx,Ptr,Min) = Eng_Rec.Ct_Head.Time(Jx)
		Enddo
	      Endif
	      If(Numex(Ptr,Min) .Ge. Filter) Then
		Xtrm_Rec.MinMax(Min).V_and_I(Ix) = Engval(Ptr,Min)
		Do Jx = 1, 2
		  Xtrm_Rec.MinMax(Min).V_and_I_Time(Jx,Ix) = Timex(Jx,Ptr,Min)
		Enddo
	        Xtrm_Table(Ptr,Min) = .True.
	        Cur_Table(Ptr) = .True.
		Exceed = .True.
		New_Xtrm = .True.
	        Numex(Ptr,Min) = 0
	      Endif
	    Else 	! Eng_Rec value not < Xtrm_Rec or Enable_Flag not set
	      Numex(Ptr,Min) = 0
	    Endif

!       Check the maxima.

	    If((Tempval(Ix) .Gt. Xtrm_Rec.MinMax(Max).V_and_I(Ix))
	1	 .And. Enable_Flags(Ptr) .And. (Neval .Ne. FV)) Then
	      Numex(Ptr,Max) = Numex(Ptr,Max) + 1
	      If((Numex(Ptr,Max) .Eq. One) .Or. 
	1	(Tempval(Ix) .Gt. Engval(Ptr,Max))) Then
	        Engval(Ptr,Max) = Tempval(Ix)
	        Do Jx = 1, 2
	          Timex(Jx,Ptr,Max) = Eng_Rec.Ct_Head.Time(Jx)
		Enddo
	      Endif
	      If(Numex(Ptr,Max) .Ge. Filter) Then
		Xtrm_Rec.MinMax(Max).V_and_I(Ix) = Engval(Ptr,Max)
		Do Jx = 1, 2
		  Xtrm_Rec.MinMax(Max).V_and_I_Time(Jx,Ix) = Timex(Jx,Ptr,Max)
		Enddo
	        Xtrm_Table(Ptr,Max) = .True.
	        Cur_Table(Ptr) = .True.
		Exceed = .True.
		New_Xtrm = .True.
	        Numex(Ptr,Max) = 0
	      Endif
	    Else 	! Eng_Rec value not > Xtrm_Rec or Enable_Flag not set
	      Numex(Ptr,Max) = 0
	    Endif
	Enddo		! Check the V_and_I for extrema

!       Check the LMAC temperatures for analog and digital converters for 
!	exceeded mimima and maxima.

	Call Lib$Movc3(8, Eng_Rec.En_Tail.LMAC_Analog_Temp, Tempval(1))
	Ptr = 116

	Do Ix = 1, 2
	    Neval = Nint(Tempval(Ix))
	    Ptr = Ptr + 1

!       Check the mimima.

	    If((Tempval(Ix) .Lt. Xtrm_Rec.MinMax(Min).LMACs(Ix))
	1	 .And. Enable_Flags(Ptr) .And. (Neval .Ne. FV)) Then
	      Numex(Ptr,Min) = Numex(Ptr,Min) + 1
	      If((Numex(Ptr,Min) .Eq. One) .Or. 
	1	(Tempval(Ix) .Lt. Engval(Ptr,Min))) Then
	        Engval(Ptr,Min) = Tempval(Ix)
	        Do Jx = 1, 2
	          Timex(Jx,Ptr,Min) = Eng_Rec.Ct_Head.Time(Jx)
		Enddo
	      Endif
	      If(Numex(Ptr,Min) .Ge. Filter) Then
		Xtrm_Rec.MinMax(Min).LMACs(Ix) = Engval(Ptr,Min)
		Do Jx = 1, 2
		  Xtrm_Rec.MinMax(Min).LMAC_Time(Jx,Ix) = Timex(Jx,Ptr,Min)
		Enddo
	        Xtrm_Table(Ptr,Min) = .True.
	        Cur_Table(Ptr) = .True.
		Exceed = .True.
		New_Xtrm = .True.
	        Numex(Ptr,Min) = 0
	      Endif
	    Else 	! Eng_Rec value not < Xtrm_Rec or Enable_Flag not set
	      Numex(Ptr,Min) = 0
	    Endif

!       Check the maxima.

	    If((Tempval(Ix) .Gt. Xtrm_Rec.MinMax(Max).LMACs(Ix))
	1	 .And. Enable_Flags(Ptr) .And. (Neval .Ne. FV)) Then
	      Numex(Ptr,Max) = Numex(Ptr,Max) + 1
	      If((Numex(Ptr,Max) .Eq. One) .Or. 
	1	(Tempval(Ix) .Gt. Engval(Ptr,Max))) Then
	        Engval(Ptr,Max) = Tempval(Ix)
	        Do Jx = 1, 2
	          Timex(Jx,Ptr,Max) = Eng_Rec.Ct_Head.Time(Jx)
		Enddo
	      Endif
	      If(Numex(Ptr,Max) .Ge. Filter) Then
		Xtrm_Rec.MinMax(Max).LMACs(Ix) = Engval(Ptr,Max)
		Do Jx = 1, 2
		  Xtrm_Rec.MinMax(Max).LMAC_Time(Jx,Ix) = Timex(Jx,Ptr,Max)
		Enddo
	        Xtrm_Table(Ptr,Max) = .True.
	        Cur_Table(Ptr) = .True.
		Exceed = .True.
		New_Xtrm = .True.
	        Numex(Ptr,Max) = 0
	      Endif
	    Else 	! Eng_Rec value not > Xtrm_Rec or Enable_Flag not set
	      Numex(Ptr,Max) = 0
	    Endif
	Enddo		! Check the LMACs for extrema

	Return
	End
