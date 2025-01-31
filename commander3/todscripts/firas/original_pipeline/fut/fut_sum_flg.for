       	INTEGER*4 FUNCTION FUT_SUM_FLG ( CHANNEL, QUALFLAGS )
C/
C/ PROGRAM NAME: 
C/	FUT_SUM_FLG
C/
C/ PROGRAM DESCRIPTION:
C/	This routine will count RED and YELLOW flags and determine a value 
C/	for the overall quality summary flag.
C/
C/ AUTHOR:
C/	Shirley M. Read ( STX )  January 1988.
C/
C/ MODIFIED BY:
C/	Shirley M. Read ( STX )  February 1988. Comment out signal of bad
C/		data quality. GRTs for test data were not functioning and
C/		too many messages were being output.
C/
CH CHANGE LOG:
CH
CH      Version 4.1.1 10/19/88, SPR 2657, Shirley M. Read, STX
CH	   	FDQ Convert needs to convert the new LMACs. The conversion
CH	        coefficients have now been defined and the include text file
CH	        has been expanded to include the LMACs. The LMACs flags must
CH		noe be counted.
CH      Version 4.1.1 10/20/88, SPR 2667, Shirley M. Read, STX
CH		FDQ has an error in counting the number of quality flags set
CH		for exceeding limits or encounters other bad data. The bit
CH		flags can have a negative value due to the msb being set.
CH	        In addition, this routine counted the Flg_Badsci twice.
CH		Note: The counting algorithms for Build 4.1.1 and ultimately
CH	        for Build 4.2 are based on the current algorithms for setting
CH		the data quality flags in FDQ. There are still some questions
CH		about setting the flags according to the FIRAS side which
CH		is powering the particular instrument component. The post 
CH		I and T analysis will determine the final algorithms. The
CH		analogs most likely to be affected are the IPDU temperatures,
CH		the drive box temperatures and the MTM cal motor current.
CH		These are the primary items changed for this build.
CH	Version 4.2.1 02/09/89, SPR 3062, Shirley M. Read, STX
CH		Moon contaminated IFGs failed by data quality checking. FES and
CH		FCI will interpret the data as being unusable. It is needed for
CH		moon modeling. The moon flag bit will not be counted in the
CH		data quality summary flag.
CH	Version 4.2.1 02/08/89, SPR 3297, Shirley M. Read, STX
CH		An attitude data quality summary flag should maintained apart
CH		from the existing over-all data quality summary flag. It will
CH		be assigned, just as the existing summary flag, a quality
CH		level dependent on the number of yellow and red limit 
CH		violations of the attitude quantities checked. FDQ will fill
CH		the attitude summary flag.
CH
C----------------------------------------------------------------------------
C/
C/ CALLING SEQUENCE:
C/	Status = FUT_Sum_Flg ( Channel, Qualflags )
C/
C/ INPUT PARAMETERS:
C/	Channel    -- Current channel being processed
C/	Qualflags  -- 110 quality flags of current science record
C/
C/ OUTPUT PARAMETERS:
C/	The summary flag in Qualflags
C/
C/ ENTRY/EXIT:
C/	Normal function entry and return.
C/
C/ ERROR HANDLING
C/	Calls to Lib$Signal and interface with Fut_Error condition handler.
C/	
C/ INPUT/OUTPUT DISK FILES: 
C/	None
C/
C/ PROCS/MACROS/INCLUDE FILES:
C/	FUT_QUALFLAGS
C/
C/ SUBROUTINES CALLED:
C/
C/ SHARED DATA AREAS:
C/
C/ METHOD USED: 
C/	Count the red and yellow flags. Determine the value of the overall
C/	quality summary flag according to the counts. Also count the red and
C/	yellow attitude flags. If the summary flag exceeds an acceptable value,
C/	send a message through Lib$Signal.
C/      Set function value to success for return.
C/
C/ PDL for FUT_SUM_FLG:
C/ 
C/	BEGIN
C/
C/	Initialize red and yellow counts to zero.
C/	Count the red and yellow flags for quality flags having only a red
C/	or yellow value.
C/	Count the red and yellow flags for attitude quantities.
C/	Count red and yellow flags from quality flags having values of one
C/	or two for yellow and red, respectively.
C/	Determine the value of the summary quality flag from the red and 
C/	yellow counts.
C/	Determine the value of the attitude summary quality flag from the red 
C/	and yellow counts.
C/	If the telemetry is bad, set the summary flag to telemetry failure.
C/	If the value of the summary flag indicates unacceptable data, send
C/	a message via Lib$Signal.
C/
C/	SET the return status for the function to success.
C/
C/	RETURN
C/	END
C/
C/ SPECIAL HANDLING:
C/	None
C/
C**************************************************************************************

	IMPLICIT NONE
	
	EXTERNAL   	FUT_NORMAL
	EXTERNAL	FUT_ABERR
	EXTERNAL	FUT_SUMQUALBAD
	Include         '(FUT_QUALFLAGS)'
	Include	        '(FUT_PARAMS)'

!	Input/Output Passed Parameters

	integer*2	CHANNEL		! Micro-processor channel
	byte		QUALFLAGS(110)	! Science record data quality flags

!	Local Variables
	
	integer*4 	retstat		! Return status
	integer*4 	success / 1 /, error / 2 /  ! Values for status
	integer*4	status		! Dummy status variable
	integer*4       chani4, flag
	integer*2 	red, yellow     ! Counters for red and yellow flags
	integer*2 	att_red, att_yellow  ! Counters for red and yellow flags
	byte            zero	        ! Zero value
	parameter       ( zero = 0 )
	byte		one             ! Value one
	parameter       ( one = 1 )	
	byte	        two             ! Value two
	parameter       ( two = 2 )
	integer*2	no_red / 0 /    ! No red flags	
	integer*2       no_yel / 0 /	! No yellow flags
	integer*2       some_yel  / 2 / ! Some yellow flags
        integer*2       some_red  / 2 / ! Some red flags
	integer*2	ix, jx
	integer*2       attitude_flag   ! Temporary attitude quality flag
!	integer*2	moon_mask / '00FB'x /  ! Mask out moon bit in byte flag

C**	Set return status to success.

	Retstat = Success

!	Initialize red, yellow and engineering limits violation
!	counts to zero.

	red = 0
	yellow = 0
	att_red = 0
	att_yellow = 0

!	Count the red and yellow flags for quality flags having only a red
!	or yellow value.

	DO Ix = 1, 5					! Red Group 1, 2 & 3
	  If ( Qualflags(Ix) .ne. zero ) red = red + 1
	ENDDO

	attitude_flag = Qualflags(FLG_ATT_ALT_RED) 
	DO Ix = 0, 7
	  If (BTest(attitude_flag,Ix))			! Red attitude Gr 3
	1	att_red = att_red + one 
	ENDDO

	attitude_flag = Qualflags(FLG_ATT_ALT_YEL) 
	DO Ix = 0, 7
	  If ( BTest(attitude_flag,Ix))			! Yellow attitude Gr 3
	1	att_yellow = att_yellow + one  
	ENDDO

	DO Ix = 8, 10					! Red Group 2
	  If ( Qualflags(Ix) .ne. zero ) red = red + 1
	ENDDO

	DO Ix = 11, 18					! Red Group 4
	  If ( Qualflags(Ix) .ne. zero ) Then
	    red = red + 1
	  End If
	ENDDO

	DO Ix =19, 26					! Yellow Group 4
	  If ( Qualflags(Ix) .ne. zero ) Then
	    yellow = yellow + 1
	  End If
	ENDDO

	If ( Qualflags(FLG_TC_Red) .ne. zero ) Then	! Group 4
	  red = red + 1
	End If
	If ( Qualflags(FLG_TC_Yel) .ne. zero ) Then	! Group 4
	  yellow = yellow + 1
	End If
	 	 
!	Count red and yellow flags from quality flags having values of one
!	or two for yellow and red, respectively.

	If ( Qualflags(FLG_IPDU_TMP) .eq. one ) then
	  yellow = yellow + 1
	Elseif ( Qualflags(FLG_IPDU_TMP) .eq. two ) then
	  red = red + 1
	Endif
	If ( Qualflags(FLG_IPDU_TMP + 1) .eq. one) then
	  yellow = yellow + 1
	Elseif ( Qualflags(FLG_IPDU_TMP + 1) .eq. two ) then
	  red = red + 1
	endif	

	If ( Qualflags(FLG_CHAN_TMP_ST) .eq. one ) then
	  yellow = yellow + 1
	Elseif ( Qualflags(FLG_CHAN_TMP_ST) .eq. two ) then
	  red = red + 1
	Endif

	If ( Qualflags(FLG_DBOX_TMP) .eq. one ) then
	  yellow = yellow + 1
	Elseif ( Qualflags(FLG_DBOX_TMP) .eq. two ) then
	  red = red + 1
	Endif
	If ( Qualflags(FLG_DBOX_TMP + 1) .eq. one ) then
	  yellow = yellow + 1
	Elseif ( Qualflags(FLG_DBOX_TMP + 1) .eq. two ) then
	  red = red + 1
	Endif	

	DO Ix = 34, 39					! Group 4
	  If ( Qualflags(Ix) .eq. one ) then
	    yellow = yellow + 1
	  Elseif ( Qualflags(Ix) .eq. two ) then
	    red = red + 1
	  Endif
	ENDDO

	If ( Qualflags(FLG_MTMCAL_MTR + 1) .gt. zero ) Then	! Group 4
	  yellow = yellow + 1
	Elseif ( Qualflags(FLG_MTMCAL_MTR + 1) .eq. two ) then
	  red = red + 1
	End If

	DO Ix = 41, 57					! Group 4
	  If ( Qualflags(Ix) .eq. one ) then
	    yellow = yellow + 1
	  Elseif ( Qualflags(Ix) .eq. two ) then
	    red = red + 1
	  Endif
	ENDDO

	DO Ix = 1, 4					! Group 6
	  Jx = FLG_STCHG_MJ_ST + Ix - 1
	  If ( Qualflags(Jx) .ne. zero ) red = red + 1
	ENDDO

	DO Ix = 84, 85					! LMACs
	  If ( Qualflags(Ix) .eq. one ) then
	    yellow = yellow + 1
	  Elseif ( Qualflags(Ix) .eq. two ) then
	    red = red + 1
	  Endif
	ENDDO


!	Group 7 will be added at a later build.
	
!	Determine the value of the overall summary quality flag from the red 
!	and yellow counts.

	If (( yellow .eq. no_yel ) .and. ( red .eq. no_red )) then
	   Qualflags(FLG_SUMMARY) = fac_good_dq
	Elseif ( yellow .le. some_yel .and. red .eq. no_red ) then
	   Qualflags(FLG_SUMMARY) = fac_some_yellow
	Elseif ( yellow .gt. some_yel .and. red .eq. no_red ) then
	   Qualflags(FLG_SUMMARY) = fac_many_yellow
	Elseif ( yellow .eq. no_yel .and. red .le. some_red ) then
	   Qualflags(FLG_SUMMARY) = fac_no_y_some_r
	Elseif ( yellow .le. some_yel .and. red .le. some_red ) then
	   Qualflags(FLG_SUMMARY) = fac_some_y_some_r
	Elseif ( yellow .gt. some_yel .and. red .le. some_red ) then
	   Qualflags(FLG_SUMMARY) = fac_many_y_some_r
	Elseif ( red .gt. some_red ) then
	   Qualflags(FLG_SUMMARY) = fac_many_red
	Endif

!	Determine the value of the attitude summary quality flag from 
!	the red and yellow counts.

        If (qualflags(flg_no_attitude) .eq. 1) then
           Qualflags(flg_att_sum) = fac_no_att_soln
	Elseif (( att_yellow .eq. no_yel ) .and. ( att_red .eq. no_red )) then
	   Qualflags(FLG_ATT_SUM) = fac_good_dq
	Elseif ( att_yellow .le. some_yel .and. att_red .eq. no_red ) then
	   Qualflags(FLG_ATT_SUM) = fac_some_yellow
	Elseif ( att_yellow .gt. some_yel .and. att_red .eq. no_red ) then
	   Qualflags(FLG_ATT_SUM) = fac_many_yellow
	Elseif ( att_yellow .eq. no_yel .and. att_red .le. some_red ) then
	   Qualflags(FLG_ATT_SUM) = fac_no_y_some_r
	Elseif ( att_yellow .le. some_yel .and. att_red .le. some_red ) then
	   Qualflags(FLG_ATT_SUM) = fac_some_y_some_r
	Elseif ( att_yellow .gt. some_yel .and. att_red .le. some_red ) then
	   Qualflags(FLG_ATT_SUM) = fac_many_y_some_r
	Elseif ( att_red .gt. some_red ) then
	   Qualflags(FLG_ATT_SUM) = fac_many_red
	Endif

!	If the telemetry is bad, set summary flag to telemetry failure.

	If (( Qualflags(FLG_BADHKP) .gt. zero ) .or. 
	1    ( Qualflags(FLG_BADSCI) .gt. zero )) then
	   Qualflags(FLG_SUMMARY) = fac_tlm_error
	Endif

!	If the value of the summary flag indicates unacceptable data, send
!	a message via Lib$Signal.

	If ( Qualflags(FLG_SUMMARY) .gt. two ) then
	    flag = Qualflags(FLG_SUMMARY)
	    chani4 = channel
!	    call Lib$Signal( FUT_SUMQUALBAD, %val(2),
!	1	 %val(flag), %val(chani4))   
	Endif

!	Set function to return status

	IF ( Retstat .eq. Success) Then
	  FUT_SUM_FLG = %loc(FUT_NORMAL)
	Else
	  FUT_SUM_FLG = %loc(FUT_ABERR)
	Endif

	Return
	End

