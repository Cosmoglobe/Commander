C-------------------------------------------------------------------------------
	Integer*4 Function FRD_Read_Sum_AVG ( Hkp_rec, First_time, End_segment,
	1                                     Fex_av_calrs_rec )
C-------------------------------------------------------------------------------

C  Purpose:  This program compute the average as well as standard deviation
C           of the calibrator resistors high and low current. 
C
C  Programmer: Nilo G. Gonzales/Harte Wang, STX
C              April 5, 1991
C
C  Invocation:   Status = FRD_Read_Sum_AVG ( Hkp_rec, First_time,
C                              End_segment, Fex_av_calrs_rec )
C
C  Include files:
C
C	CT$Library:CTUser.Inc
C       $SSDef
C	FUT_CALRS
C	Fut_Error
C
C  Method Used:  PDL for FRD_READ_SUM_AVG
C
C    Begin  
C	Do i =1, 4
C          Set sidemax, sidemin arrays with nominal 
C          Calibrator resistors values in unit counts. 
C	Enddo
C   Initialize counters.
C	If (First_time) Then
C	   Set good_rec counter to 0
C	   Set bad_rec counter to 0
C	   Do for each (I,J) 1 through 4
C	      Set Resistor sum to 0
C	      Set Resistor counter to 0
C             Set bad points to 0 
C	      Set Sigma to 0.0
C	      Set Sum of square to 0.0
C	      Set Real resistor average to 0.0
C	   Enddo
C	Endif
C	If (.Not. End_segment) Then
C          Do Major Frame 1 and 2
C             If (Telemetry is good) Then
C                 Add 1 to Good frame counter
C                 Do side 1 through 4
C                    Do I 1 throught 4
C                       If (HKP calres data is within nominal counts) Then
C                          Add 1 to Resistor counter (side,I) 
C                          Resistor Sum (side,i) = Resistor Sum (side,i) +
C                          HKP calres data (I) - Offset(side,i)
C                          Sum_Square(side,I) = Sum_Square(side,I) +
C                                               (HKP calres data (I) - 
C                                               Offset(side,i)) ** 2
C                       Else
C                           Add 1 to Bad points (side,I)
C                       Endif
C                    Enddo
C                  Enddo
C             Else
C                Add 1 to Bad frame counter
C             Endif
C       Else         ! End segment
C   Check if Resistor counter(side,I) is greater than 0, Then compute offset
C   averages, resistor averages in real numbers, sum of averages squares and
C   Sigmas.
C          Do side 1 through 4
C             Do I 1 through 4
C                If (Resistor counter is greater than 0) Then
C                    offset average = resistor sum / resistor counter
C                    Real Resistor average = offset average + offset
C                    If (Resistor counter is greater than 1) Then
C                       Sigma = Sqrt(Sum Square - Resistor counter *
C                               offset average **2)/resistor counter/(resistor
C                               counter - 1)
C                    Else
C                        Sigma = - 1.0
C                    Endif
C                Else
C                    Real Resistor average = 0.0
C                    Sigma = - 1.0
C                Endif
C             Enddo
C          Enddo
C          Set reference fields structure of good and bad frames
C          Do for each (side,I) 1 through 4
C             set reference fields structure of calres averages, Calres
C             deviation averages and bad points
C          Enddo 
C       Endif
C       Return
C    End PDL
C-------------------------------------------------------------------------------
	Implicit  None

C  Passed Parameters

	Dictionary   'NFS_HKP'
	Record       /NFS_HKP/ Hkp_Rec

	Dictionary   'FEX_AV_CALRS'
	Record       /FEX_AV_CALRS/ FEX_AV_CALRS_REC
	Logical*1    First_Time, End_Segment     ! Flags

C  Include Files

	Include	     'CT$Library:CTUser.Inc'
	Include	     '($SSDef)'
	Include      '(FUT_CALRS)'
	Include	     '(Fut_Error)'

C  Externals

  	External	frd_normal
	External	frd_aberr

C  Local variables

	Integer*4     Res_Ctr (4,4)    ! Resistors counter
	Real*4        Res_Sum (4,4)    ! Resistors Sum counter
	Integer*4     Good_rec /0/     ! Good record counter
	Integer*4     Bad_rec  /0/     ! Bad record counter
	Integer*2     Bad_Pnts (4,4)   ! Bad record counter
	Real*4        RAve_Res(4,4)    ! Average resistor counter (real)
	Real*4        Sigma(4,4)       ! Sigma
	Real*4        Sum_sq(4,4)      ! Sum of the Square
	Real*4        offsetave        ! Offset average
	Integer*2     I, j, Side, mj   ! Indexs
	Integer*2     sidemin(4,4), sidemax(4,4)
	real*4         offset(4,4) /10242.0,  4925.0, 10252.0,  4923.0,
	1                          14314.0, 10267.0, 14325.0, 10271.0,
	2                          15648.0, 14322.0, 15659.0, 14332.0,
	3                          15983.0, 15344.0, 15995.0, 15355.0/
	Integer*4     n, Status
	real*4        x1

	FRD_Read_Sum_AVG = %loc (FRD_Normal)

C   Set sidemax, sidemin arrays with high and low calibrator variables.

	Do i =1, 4
	  sidemax(1,i) = A_lo_max(i)
	  sidemax(2,i) = A_hi_max(i)
	  sidemax(3,i) = B_lo_max(i)
	  sidemax(4,i) = B_hi_max(i)
	  sidemin(1,i) = A_lo_min(i)
	  sidemin(2,i) = A_hi_min(i)
	  sidemin(3,i) = B_lo_min(i)
	  sidemin(4,i) = B_hi_min(i)
	Enddo

C   Initialize counters.

	If (First_time) then
	  Good_rec = 0
	  BAD_rec = 0
	  Do i=1,4
	    Do j=1,4
	      Res_sum(i,j)  = 0
	      Res_ctr(i,j)  = 0
              bad_pnts(i,j) = 0 
	      Sigma(i,j)    = 0.0
	      Sum_sq(i,j)   = 0.0
	      rave_res(i,j) = 0.0
 
	    Enddo
	  Enddo
	Endif
	If (.Not. End_segment) then
C   Checks good telemetry quality, nominal cal resistor ranges in counts
C   and increment flags.
	  Do mj=1,2
	    If (hkp_rec.mj_frm(mj).tlm_qual_maj_frm .eq. 0) Then
	      GOOD_REC = GOOD_REC + 1
	      Do side=1,4
	        Do i=1,4
	          If (hkp_rec.frame(mj).temps.side_amp(side).cal_resist(i)
	1                    .Le. sidemax(side,i)
	2           .And. hkp_rec.frame(mj).temps.side_amp(side).cal_resist(i)
	3                    .Ge. sidemin(side,i)) Then
	            Res_Ctr(side,i) = Res_Ctr(side,i) + 1
	            Res_sum(side,i) = Res_sum(side,i) +
	1              hkp_rec.frame(mj).temps.side_amp(side).cal_resist(i) -
	2              offset(side,i)
	            Sum_sq(side,i)= sum_sq(side,i) +
	1              (hkp_rec.frame(mj).temps.side_amp(side).cal_resist(i) -
	2               offset(side,i)) ** 2
                  Else
                    Bad_pnts(side,i) = Bad_pnts(side,i) + 1
	          Endif
	        Enddo
	      Enddo
	    Else
	      Bad_rec = Bad_rec + 1
C              Do side=1,4
C                Do i=1,4
C                 Bad_pnts(side,i) = bad_pnts(side,i) + 1
C                Enddo
C              Enddo 
	    Endif
	  Enddo
	Else ! End segment

C   Check if Res_ctr(I) is greater than 0, then computes offset averages, 
C   resistors averages in real numbers, sum of averages square and Sigmas.
           Do side=1,4
	    Do i=1,4
	     If (Res_ctr(side,i) .Gt. 0) then
	        n = res_ctr(side,i)
	        offsetave = res_sum(side,i)/n
	        Rave_res(side,i) = offsetave + offset(side,i)
	        If (n .Gt. 1) Then
	           x1 = (SUM_SQ(side,i) - n* offsetave**2)/n/(n-1)
	           Sigma(side,i) = Sqrt(x1)
	        Else
	           Sigma(side,i) = - 1.0
	           FRD_Read_Sum_AVG = %loc (FRD_Aberr)
	           Call lib$signal (%val(status))
	        Endif
	     Else
	       RAve_res(side,i) = 0.0
	       Sigma(side,i) = - 1.0
	       FRD_Read_Sum_AVG = %loc (FRD_Aberr)
	       Call lib$signal (%val(status))
	     Endif
	    Enddo
	  Enddo

C    Set reference structure of good and bad fields.
 
	  Fex_av_calrs_rec.num_good_record = good_rec
	  Fex_av_calrs_rec.num_bad_record = bad_rec
	  Do side = 1,4

C    Set reference structure fields of calres averages, Calres deviation 
C    averages and bad points. 

	    If (side .eq. 1) then
	      do i=1,4
	        fex_av_calrs_rec.calres_Bad_PNTS_a_lo(i) = bad_pnts(side,i)
	        fex_av_calrs_rec.calres_ave_a_lo(i) = Rave_res(side,i)
	        fex_av_calrs_rec.calres_ave_a_lo_dev(i) = Sigma(side,i)
	      Enddo
	    Endif
	    If (side .eq. 2) then
	      do i=1,4
	        fex_av_calrs_rec.calres_Bad_PNTS_a_hi(i) = bad_pnts(side,i)
	        fex_av_calrs_rec.calres_ave_a_Hi(i) = Rave_res(side,i)
	        fex_av_calrs_rec.calres_ave_a_Hi_dev(i) = Sigma(side,i)
	      Enddo
	    Endif
	    If (side .eq. 3) then
	      do i=1,4
	        fex_av_calrs_rec.calres_Bad_PNTS_b_lo(i) = bad_pnts(side,i)
	        fex_av_calrs_rec.calres_ave_B_lo(i) = Rave_res(side,i)
	        fex_av_calrs_rec.calres_ave_B_lo_dev(i) = Sigma(side,i)
	      Enddo
	    Endif
	    If (side .eq. 4) then
	      do i=1,4
	        fex_av_calrs_rec.calres_Bad_PNTS_B_Hi(i) = bad_pnts(side,i)
	        fex_av_calrs_rec.calres_ave_B_Hi(i) = Rave_res(side,i)
	        fex_av_calrs_rec.calres_ave_B_Hi_dev(i) = Sigma(side,i)
	      Enddo
	    Endif
	  Enddo
	Endif

	Return

	End
