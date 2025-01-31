	INTEGER*4 FUNCTION FDQ_SET_LIMFLAGS ( LIMFLAGS, 
	1	SC_SIDE, ATTITUDE_TYPE, LIMIT_ON ) 
C/
C/ PROGRAM NAME: 
C/	FDQ_SET_LIMFLAGS
C/
C/ PROGRAM DESCRIPTION:
C/	This routine sets the 110 flags in the internal array, Limit_On, 
C/	to enable or disable limit checking. The criteria for setting the
C/	flags are obtained from the time_tagged archive dataset, Limflags,
C/	and in some cases the side of the spacecraft powering FIRAS.
C/
C/ AUTHOR:
C/	Shirley M. Read ( STX )  May 1988.
C/
C/ MODIFIED BY:
C/        H. Wang
C/        STX
C/        Jan. 29, 1991
C/        Reason:   New requirements for FDQ
C/             Flg_moon has been replaced by Flg_NO_ATTITUDE
C/
C/      STX for Observatory I and T, October 1988. Started Change Log Format.
C----------------------------------------------------------------------------
CH CHANGE LOG:
CH
CH	Version 4.1.1 10/15/88, SPR 2620, Shirley M. Read, STX
CH		Disable limit checking of attitude fields for /ATT=NONE.
CH
CH	Version 4.2.1 02/09/89, SPR 3062, Shirley M. Read, STX
CH		Moon contaminated IFGs failed by data quality checking. FES and
CH		FCI will interpret the data as being unusable. It is needed for
CH		moon modeling. The moon flag bit will not be counted in the
CH		data quality summary flag.
CH
CH	Version 4.2.2 03/29/89, SPR 3495, Shirley M. Read, STX
CH		Disable limit checking of Moon Angle to set Flg_Moon and 
CH	        Flg_Att_Sum for command line option, /ATT=NONE.
CH
C----------------------------------------------------------------------------
C/
C/ CALLING SEQUENCE:
C/	Status = FDQ_Set_Limflags( Limflags, SC_Side, Attitude_Type, Limit_On )
C/
C/ INPUT PARAMETERS:
C/	Name          Type                Description
C/	----------------------------------------------------------------- 
C/	Limflags      Record Fex_Limflags Enable/disable flags for limit
C/			   		  checking
C/	SC_Side       I*2                 Spacecraft side powering FIRAS
C/	Attitude_Type I*4                 Type of attitude solution: 0 = none,
C/					  1 = sim, 2 - 6 real SC types 
C/
C/ OUTPUT PARAMETERS:
C/	Name          Type                Description
C/	----------------------------------------------------------------- 
C/	Limit_On(110) L*1		  Limit check flags corresponding
C/					  to the 110 data quality flags
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
C/	FUT_Qualflags
C/      
C/
C/ SUBROUTINES CALLED:
C/
C/ SHARED DATA AREAS:
C/
C/ METHOD USED: 
C/	This routine will set the flags of the output Limit_On array for the 
C/	next science set to true. Then by examining the current Limflags
C/	dataset and when applicable, the spacecraft side and attitude solution
C/	type, the corresponding flag in the Limit_On array will be set to 
C/	false if indicated.
C/
C/ SPECIAL HANDLING:
C/	None
C/
C**************************************************************************************

	IMPLICIT NONE
	
	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	Include         '(FUT_QUALFLAGS)'
	Include         '(FUT_PARAMS)'		!Added for attitude type
	
!	Input/Output Passed Parameters

	Dictionary 'Fex_Limflags'
	Record /Fex_Limflags/ Limflags

	Integer*2	SC_Side 	! Spacecraft side 1 or 2
	Integer*4       Attitude_Type   ! Type of attitude solution
	Logical*1	Limit_On(110)   ! Limits on/off flags

!	Local Variables
	
	Integer*4 	retstat		! Return status
	Integer*4 	success / 1 /, error / 2 /  ! Values for status
	Integer*2       zero / 0 /, one / 1 /, two / 2 /
	Integer*4	status		! Dummy status variable
	Integer*2	ix, jx
	Logical*1       offred, offyel  ! Flags for setting to false
	Integer*2       Amask(10) / 1,3,5,7,9,11,13,15,17,19 /
	Integer*2       Bmask(10) / 2,4,6,8,10,12,14,16,18,20 /

!	Set return status to success.

	retstat = success

!       Initialize
 
	Do ix = 1, 110
	   Limit_On(ix) = .true.
	Enddo

!	Set the science attitude limits. 
!       The red and yellow attitude flags must be checked for "attitude =
!       none" from the command line input. The original code loop is 
!       therefore split into two loops. SMR 10/15/88.

	Do ix = 1, 3
	  If ( Limflags.Lim_Flags.Sci_Att(ix) .eq. zero ) Then
	     Limit_On(ix) = .false.
	  Endif
	Enddo
	If (( Limflags.Lim_Flags.Flg_Micro(1) .eq. zero ) .and.
	1  ( Limflags.Lim_Flags.Flg_Micro(2) .eq. zero )) Then
	     Limit_On(4) = .false.
	     Limit_On(5) = .false.
	Endif
	Do ix = 6, 7
	  If (( Limflags.Lim_Flags.Sci_Att(ix) .eq. zero ) .or.
	1     ( Attitude_Type .eq. fac_none )) Then
	     Limit_On(ix) = .false.
	  Endif
	Enddo
	Do ix = 8, 10
	  If ( Limflags.Lim_Flags.Sci_Att(ix) .eq. zero ) Then
	     Limit_On(ix) = .false.
	  Endif
	Enddo

!       Set the GRT flags. Ignore the TC flags (no limits for Build 4.0).

	offred = .true.
	offyel = .true.
	Do ix = 1, 8
	  If ( Limflags.Lim_Flags.Flg_Grtred_St(ix) .ne. zero ) Then
	     offred = .false.
	  Endif
	  If ( Limflags.Lim_Flags.Flg_Grtyel_St(ix) .ne. zero ) Then
	     offyel = .false.
	  Endif
	Enddo
	If ( offred ) Then
	jx = Flg_Grtred_St -1
	  Do ix = 1, 8
	    jx = jx + 1
	    Limit_On(jx) = .false.
	  Enddo
	Endif
	If ( offyel ) Then
	jx = Flg_Grtyel_St -1
	  Do ix = 1, 8
	    jx = jx + 1
	    Limit_On(jx) = .false.
	  Enddo
	Endif

!	Set the Engineering Analog Flags. Leave channel specific alone.

	If ( SC_Side .eq. one ) Then  ! Analog Conv. A, Dig.Conv. B
	   If ( Limflags.Lim_Flags.Flg_Ipdu_Tmp(1) .eq. zero ) 
	1       Limit_On(Flg_Ipdu_Tmp) = .false.
	   If ( Limflags.Lim_Flags.Flg_Ipdu_Tmp(2) .eq. zero ) 
	1       Limit_On(Flg_Ipdu_Tmp + 1) = .false.
	Elseif ( Sc_Side .eq. two ) Then  ! Analog Conv. B, Dig. Conv. A
	   If ( Limflags.Lim_Flags.Flg_Ipdu_Tmp(3) .eq. zero ) 
	1       Limit_On(Flg_Ipdu_Tmp) = .false.
	   If ( Limflags.Lim_Flags.Flg_Ipdu_Tmp(4) .eq. zero ) 
	1       Limit_On(Flg_Ipdu_Tmp + 1) = .false.
	Endif
	If ( Sc_Side .eq. one ) Then  ! Drive box temp. SC 1 reads A, 2 reads B
	   Limit_On(Flg_Dbox_Tmp+1) = .false.
	   If ( Limflags.Lim_Flags.Flg_Dbox_Tmp(1) .eq. zero )
	1     Limit_On(Flg_Dbox_Tmp) = .false.
	Elseif ( Sc_Side .eq. two ) Then
	   Limit_On(Flg_Dbox_Tmp) = .false.
	   If ( Limflags.Lim_Flags.Flg_Dbox_Tmp(2) .eq. zero )
	1     Limit_On(Flg_Dbox_Tmp+1) = .false.
	Endif

	If (( Limflags.Lim_Flags.Flg_Stmon_Tmp(1) .eq. zero) .and.
	1   ( Limflags.Lim_Flags.Flg_Stmon_Tmp(2) .eq. zero)) Then
           Limit_On(Flg_Stmon_Tmp) = .false.	! Status monitor temperatures.
	Endif

	If ( SC_Side .eq. one ) Then	! Optical and channel preamps
	   If ( Limflags.Lim_Flags.Flg_Chpa_Tmp(1) .eq. zero )
	1	Limit_On(Flg_Chpa_Tmp) = .false.
	   If ( Limflags.Lim_Flags.Flg_Oppa_Tmp(1) .eq. zero )
	1	Limit_On(Flg_Oppa_Tmp) = .false.
	Elseif ( Sc_Side .eq. two ) Then
	   If ( Limflags.Lim_Flags.Flg_Chpa_Tmp(2) .eq. zero )
	1	Limit_On(Flg_Chpa_Tmp) = .false.
	   If ( Limflags.Lim_Flags.Flg_Oppa_Tmp(2) .eq. zero )
	1	Limit_On(Flg_Oppa_Tmp) = .false.
	Endif 

	If ( Limflags.Lim_Flags.Flg_HS_Heat_A(1) .eq. zero)  !Hot spot heaters
	1	Limit_On(Flg_HS_Heat_A) = .false.
	If ( Limflags.Lim_Flags.Flg_HS_Heat_B(1) .eq. zero)  !Hot spot heaters
	1	Limit_On(Flg_HS_Heat_B) = .false.

	If ( Limflags.Lim_Flags.Flg_Mtmcal_Mtr(1) .eq. zero )  
  	1    Limit_On(Flg_Mtmcal_Mtr) = .false.  !Mtm/Cal latch & motor current 
	If ( Limflags.Lim_Flags.Flg_Mtmcal_Mtr(2) .eq. zero )  
  	1    Limit_On(Flg_Mtmcal_Mtr+1) = .false.  

	jx = Flg_IPDU_Vol_St - 1	
	Do ix = 1, 10	! IPDU voltages, 10 flags can be set on A or B side
	  jx = jx + 1
	  If ((Limflags.Lim_Flags.Flg_IPDU_Vol_St(Amask(ix)) .eq. zero) .and.
	1   (Limflags.Lim_Flags.Flg_IPDU_Vol_St(Bmask(ix)) .eq. zero)) 
	1	Limit_On(jx) = .false.
	Enddo	  
	jx = Flg_IPDU_Cur_St - 1	
	Do ix = 1, 6	! IPDU currents, 6 flags can be set on A or B side
	  jx = jx + 1
	  If ((Limflags.Lim_Flags.Flg_IPDU_Cur_St(Amask(ix)) .eq. zero) .and.
	1   (Limflags.Lim_Flags.Flg_IPDU_Cur_St(Bmask(ix)) .eq. zero))
	1	Limit_On(jx) = .false.
	Enddo	  

!	Set the status changes flags and miscellaneous flags.

	jx = Flg_Stchg_Mj_St - 1
	Do ix = 1, 4
	   jx = jx + 1
	   If ( Limflags.Lim_Flags.Flg_Stchg_Mj_St(ix) .eq. zero )
	1     Limit_On(jx) = .false.
	Enddo

	jx = Flg_Dirbe_St - 1
	Do ix = 1, 22
	   jx = jx + 1
	   If ( Limflags.Lim_Flags.External(ix) .eq. zero )
	1     Limit_On(jx) = .false.	
	Enddo

	If ( Limflags.Lim_Flags.Flg_ACLmac_Tmp(1) .eq. zero )
	1     Limit_On(Flg_ACLmac_Tmp) = .false.	

	If ( Limflags.Lim_Flags.Flg_DCLmac_Tmp(1) .eq. zero )
	1     Limit_On(Flg_DCLmac_Tmp) = .false.	

!	Set the non attitude flag.

	If (( Limflags.Lim_Flags.FLG_NO_ATTITUDE(1) .eq. zero ) .or.
	1     ( Attitude_Type .eq. fac_none ))
	2     Limit_On(Flg_NO_ATTITUDE) = .false.	

!	Set the summary flags.

	If (( Limflags.Lim_Flags.Flg_Att_Sum(1) .eq. zero ) .or.
	1     ( Attitude_Type .eq. fac_none ))
	2     Limit_On(Flg_Att_Sum) = .false.	

	If ( Limflags.Lim_Flags.Flg_Summary(1) .eq. zero )
	1     Limit_On(Flg_Summary) = .false.	

!	Set function to return status

	If (retstat.eq.success) then
	  FDQ_Set_Limflags = %loc(FDQ_NORMAL)
	Else
	  FDQ_Set_Limflags = %loc(FDQ_ABERR)
	Endif

	Return
	End

