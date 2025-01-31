C-------------------------------------------------------------------------------

	Integer*4 Function FXT_Get_Flags ( Eng_Time, Con_Lun, Report,
	1	  Lun_Rpt, Enable_Flags ) 

C-------------------------------------------------------------------------------
C
C	Purpose: To get the time tagged configuration record, FEX_Limflags,
C	         corresponding to the current FDQ_Eng record and set the
C	         Enable_Flags for the extrema checking.
C
C	Author: Shirley M. Read
C		STX, December, 1988
C
C	Invocation: Status = FDQ_Get_Flags ( Eng_Time, Con_Lun,
C	            Report, Lun_Rpt, Enable_Flags )
C
CH	Change Log:
CH
C	  ----------------------------------------------------------------------
C
C	Input Files:
C	  Fex_Limflags -- Time tagged limit flags, archived configuration file
C
C	Output Files:
C
C	Input Parameters:
C	  Name		      Type	Description
C 	  ----------------------------------------------------------------------
C	  Eng_Time(2)         I*4       Binary data engineering record time
C	  Con_Lun             I*4       Logical unit for Get Config
C	  Report              L*1       Flag to write report
C	  Lun_Report          I*4       Logical unit for writing report     
C	
C	Output Parameters:
C	  Name		      Type	Description
C 	  ----------------------------------------------------------------------
C	  Enable_Flags(118)   L*1       Flags to enable extrema checking
C	
C	Subroutines Called:
C
C	  CCT_Get_Config_TOD  
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C	  FXT_Msg.Txt
C	  FUT_Qualflags.Txt
C         CCT_Get_Config
C	  CT$Library:CTUser.Inc
C	  $SSDef
C
C	Processing Method:
C
C	This routine calls the CCT_Get_Config_Tod to get the FEX_Limflags
C	time tagged, configuration record corresponding to the FDQ_Eng record.
C	Then it sets the 118 flags in the internal array, Enable_Flags, 
C	to enable or disable extrema checking. The criteria for setting the
C	flags are obtained from the FEX_Limflags. There is no requirement, as
C	there is for FDQ, to check the Spacecraft Side to see which FIRAS side
C	is being sampled in the telemetry. Thus the spacecraft archive will
C	not be read to check the Spacecraft Side for the Build 4.2 version of
C	FXT_Log_Extrema.
C
C------------------------------------------------------------------------------

	Implicit None
	
!	Input/Output Passed Parameters

	Integer*4	Eng_Time(2)         ! Binary engineering record time
	Integer*4       Con_lun             ! Logical unit for Get Config
	Logical*1       Report              ! Flag to write report
	Integer*4       Lun_Rpt		    ! Logical unit for report file 	
	Logical*1	Enable_Flags(118)   ! Enable flags to check extrema

!	Include Files.

	Include		'CT$Library:CTUser.Inc'
	Include		'($SSDef)'
	Include 	'(FUT_Qualflags)'
	Include		'(FXT_Msg)'
	Include	        '(CCT_Get_Config)'

!	Structure for CCT_Get_Config_TOD

	Record /Config_Status/ Stat

!	Local Variables
	
	Dictionary 'Fex_Limflags'
	Record /Fex_Limflags/ Limflags

	Integer*4 	Status		      ! Return status
	Integer*4	Iostatus              ! Return status from I/O
	Logical*1       Limit_On_Grtred(16)   ! Red GRTs Limits enable flags
	Logical*1       Limit_On_Grtyel(16)   ! yellow GRTs Limits enable flags
	Integer*2       Grtflags
	Integer*2       Flag                    ! Flag pointer
	Integer*2       Zero / 0 /, One / 1 /, Two / 2 /
	Integer*2	Ix, Jx, Ptr, Bp, Ep   ! Indices
	Integer*2       Amask(10) / 1,3,5,7,9,11,13,15,17,19 /
	Integer*2       Bmask(10) / 2,4,6,8,10,12,14,16,18,20 /

	Integer*4 	Ndset / 1 /		! Number of datasets
	Character*27    Dset			! Names of datasets
	Data dset / 'CSDR$FIRAS_REF:FEX_LIMFLAGS' /
	Integer*4 	Size / 512 /            ! Dataset size
	Integer*4       Ncache / 1 /		! Number of caches- not used
	Integer*4       Index    		! Initial cache pointers
	Integer*4 	Ref_Count		! Reference counter for cache
	Logical*1       New_Segment		! Flag for new segments accessed
	Character*1     Access_Mode / ' ' /     ! Access mode sequential-default
	Character*19    Dummy_space		! Dummy space pad
	Character*14    Con_Time		! Data time for report

!	Functions.

	Integer*4       CCT_Get_Config_Tod

!	Set function return status to normal.

	FXT_Get_FLags = %loc(FXT_Normal)

!       Initialize
 
	Do Ix = 1, 118
	   Enable_Flags(Ix) = .True.
	Enddo
     	Call CT_Binary_To_Gmt(Eng_Time,Con_Time)

!	Get the FEX_Limflags time tagged, configuration file.

	Status = CCT_Get_Config_Tod( Eng_Time, Ndset,
	1     Size, Con_Lun, Index, Limflags, New_Segment, Stat)

	If ( .Not. Status ) Then
	  FXT_Get_Flags = %loc(FXT_Aberr)
	  Call Lib$Signal( FXT_CCTGetConfig, %val(1), %val(Status))
	  If ( Report ) Then
     	    Write(Unit=Lun_Rpt,FMT=100,Iostat=Iostatus) 
	1	  Dset, Con_Time, Status
 100 	    Format(1x/5x,'Error CCT_Get_Config_TOD for ',a,
	1	   ' Input Time= ',a, ' Status= ',z8.8)
	  Endif		! Report
	Else		! Status is OK
	  If ( Report .And. New_Segment ) Then
            Write(Unit=Lun_Rpt,FMT=200,Iostat=Iostatus) Con_Time, Dset
 200	    Format(1x/5x,'New Segment Accessed at time: ',a /
	1          5x,' for Configuration Dataset: ',a )
	  Endif 	! Report and new segment
	Endif		! Status from CCT_Get_Config

!	If function status is normal, set the Enable_Flags.

!       Set the GRT flags. First check the Limflags bit mask for individual
!       GRT fields. Using the bit masks, set the Limit_On_Grtred values and
!       Limit_On_Grtyel values. Then set the corresponding Enable_Flags for
!       the individual GRT fields to false (disable) only if both red and
!       yellow limit checking is disabled. The four passes through the 'Do Ptr'
!       loop correspond to the four GRT group fields: A low GRT, A high GRT,
!       B low GRT and B high GRT.
	
	If (FXT_Get_Flags .Eq. %loc(FXT_Normal)) Then
	  Ep = 0
	  Do Ptr = 0, 6, 2
	    Bp = Ep + 1
	    Ep = Bp + 15 
	    Do Ix = 1, 2
	      Grtflags = Limflags.Lim_Flags.FLG_GRTRED_ST(Ix+Ptr)
	      Do Jx = 1, 8
	        Limit_On_Grtred((Ix-1)*8 + Jx) = BITest( Grtflags, Jx-1)
	      Enddo
	    Enddo
	    Do Ix = 1, 2
	      Grtflags = Limflags.Lim_Flags.FLG_GRTYEL_ST(Ix+Ptr)
	      Do Jx = 1, 8
	        Limit_On_Grtyel((Ix-1)*8 + Jx) = BITest( Grtflags, Jx-1)
	      Enddo
	    Enddo
	    Jx = 0
	    Do Ix = Bp, Ep
	      Jx = Jx + 1
	      If (( .Not. Limit_On_Grtred(Jx)) .And.
	1	( .Not. Limit_On_Grtyel(Jx))) Enable_Flags(Ix) = .False.   
	    Enddo
	  Enddo		! Ptr = 0, 6, 2

!	Set the Engineering Analog Flags. 

	  Ptr = 65		! Beginning of IPDU Temperatures. Ptr 65 and 66.

!	If the Spacecraft Side equals one, then Analog Conv. A, Dig.Conv. B 
!	Elseif the Spacecraft Side equals two then Analog Conv. B, Dig. Conv. A.
!	FXT will not check the SC Side. The extrema checking will not be
!	turned off for these converters unless all flags in FEX_Limflags are
!	set to zero for both A and B converters. The assumption is made that
!	the converters not being sampled in the telemetry will de detected
!       in the conversions routine of FDQ and get the flag value of -9999.0.

	  If (( Limflags.Lim_Flags.Flg_Ipdu_Tmp(1) .eq. zero ) .And.
	1     ( Limflags.Lim_Flags.Flg_Ipdu_Tmp(3) .eq. zero )) Then
	    Enable_Flags(Ptr) = .False.
	  Endif
	  If (( Limflags.Lim_Flags.Flg_Ipdu_Tmp(2) .eq. zero ) .And.
	1     ( Limflags.Lim_Flags.Flg_Ipdu_Tmp(4) .eq. zero )) Then
	    Enable_Flags(Ptr+1) = .False.
	  Endif

	  Ptr = 66		! Set index for channel temperatures loop.
				! Ptr 67 to 70.

	  Do Ix = 1, 4
	    Ptr = Ptr + 1
	    If ( Limflags.Lim_Flags.Flg_Chan_Tmp_St(Ix) .Eq. Zero ) 
	1	 Enable_Flags(Ptr) = .False.
	  Enddo
	
! 	Drive box temperatures. SC Side 1 reads A, SC Side 2 reads B.
! 	Extrema checking will be turned off only if both Limflags are off.
!	Ptr 71 to 72.

	  Ptr = 71
	  If (( Limflags.Lim_Flags.Flg_Dbox_Tmp(1) .Eq. Zero ) .And.
	1     ( Limflags.Lim_Flags.Flg_Dbox_Tmp(2) .eq. zero )) Then
	    Enable_Flags(Ptr) = .False.
	    Enable_Flags(Ptr+1) = .False.
	  Endif

	  Ptr = 72		! Set index for status monitor temperatures.
				! Ptr 73 to 74.
	  Do Ix = 1, 2
	    Ptr = Ptr + 1
	    If ( Limflags.Lim_Flags.Flg_Stmon_Tmp(Ix) .Eq. Zero) 
	1     Enable_FLags(Ptr) = .False.
	  Enddo

!	For the channel pre_amp temperatures and optical pre-amp temperatures
!       Spacecraft Side 1 reads FIRAS A and Spacecraft Side 2 reads FIRAS B.
!	Extrema checking will be turned off only if the Limflags for both SC 
!	sides are zero. Ptr 75 and 76.

	  Ptr = 75
	  If (( Limflags.Lim_Flags.Flg_Chpa_Tmp(1) .Eq. Zero ) .And.
	1     ( Limflags.Lim_Flags.Flg_Chpa_Tmp(2) .Eq. Zero )) Then
	    Enable_Flags(Ptr) = .False.
	  Endif

	  If (( Limflags.Lim_Flags.Flg_Oppa_Tmp(1) .Eq. Zero ) .And.
	1     ( Limflags.Lim_Flags.Flg_Oppa_Tmp(2) .Eq. Zero )) Then
            Enable_Flags(Ptr+1) = .False.
	  Endif 

!	Set the hot spot heaters and MTM cal motor drive currents for 
!	FIRAS A and B sides. Ptr 77 to 80.

	  Ptr = 77
	  If ( Limflags.Lim_Flags.Flg_HS_Heat_A(1) .Eq. Zero)  !Hot spot heaters
	1   Enable_Flags(Ptr) = .False.
	  Ptr = Ptr + 1
	  If ( Limflags.Lim_Flags.Flg_HS_Heat_B(1) .Eq. Zero)  !Hot spot heaters
	1   Enable_Flags(Ptr) = .False.
	  Ptr = Ptr + 1
	  If ( Limflags.Lim_Flags.Flg_Mtmcal_Mtr(1) .Eq. Zero )  
  	1   Enable_Flags(Ptr) = .False.  !Mtm/Cal latch & motor current 
	  Ptr = Ptr + 1
	  If ( Limflags.Lim_Flags.Flg_Mtmcal_Mtr(2) .Eq. Zero )  
  	1   Enable_Flags(Ptr) = .False.  

!	Set the bolometer bias voltages for the four channels. Ptr 81 to 84.

	  Do Ix = 1, 4
	    Ptr = Ptr + 1
	    If ( Limflags.Lim_Flags.Flg_Bol_Vol_St(Ix) .Eq. Zero ) Then
	      Enable_Flags(Ptr) = .False.
	    Endif
	  Enddo

!	Set the IPDU voltages. The FIRAS sides alternate between A and B.
!	Ptr 85 to 104.

	  Do Ix = 1, 10	! IPDU voltages, 10 flags can be set on A or B side
	    Ptr = Ptr + 1
	    If (Limflags.Lim_Flags.Flg_IPDU_Vol_St(Amask(Ix)) .Eq. Zero) 
	1     Enable_Flags(Ptr) = .False.
	    Ptr = Ptr + 1
	    If (Limflags.Lim_Flags.Flg_IPDU_Vol_St(Bmask(Ix)) .Eq. Zero)
	1     Enable_FLags(Ptr)  = .False.
	  Enddo	  

!	Set the IPDU currents. The FIRAS sides alternate between A and B.
!	Ptr 105 to 116.

	  Do Ix = 1, 6	! IPDU currents, 6 flags can be set on A or B side
	    Ptr = Ptr + 1
	    If (Limflags.Lim_Flags.Flg_IPDU_Cur_St(Amask(Ix)) .Eq. Zero) 
	1     Enable_Flags(Ptr) = .False.
	    Ptr = Ptr + 1
	    If (Limflags.Lim_Flags.Flg_IPDU_Cur_St(Bmask(Ix)) .Eq. Zero)
	1     Enable_Flags(Ptr) = .False.
	  Enddo	  

!	Set the LMAC temperatures. Ptr 117 to 118.

	  Ptr = 117
	  If ( Limflags.Lim_Flags.Flg_ACLmac_Tmp(1) .Eq. Zero )
	1   Enable_Flags(Ptr) = .False.	
	  Ptr = 118
	  If ( Limflags.Lim_Flags.Flg_DCLmac_Tmp(1) .Eq. Zero )
	1   Enable_Flags(Ptr) = .False.	

	Endif	! Function status is normal

	Return
	End
