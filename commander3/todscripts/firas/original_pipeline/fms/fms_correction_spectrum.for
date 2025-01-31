	Integer*4  Function  FMS_CORRECTION_SPECTRUM
	1                    ( Arch_Cal, Specif, Model_Ext, Access_Time,
	2                      Corr_Offset_Spec, Corr_Gain_Spec )

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                          FMS_CORRECTION_SPECTRUM.FOR
C
C  Purpose:  Read the correction spectrum.
C
C  The following paragraphs are left here for history sake.  Jump right to
C  WARNING modification below.
C
C  Equations, et al:  (See "Combining FIRAS Skymaps - A General Approach")
C
C     Use integer number of ifgs if input is FCS, else floating point for FMS.
C     Compute the pixel weight WP (f,p) = RW1(f)*N1(p)*RW2(f)*N2(p) /
C                                         (RW1(f)*N1(p) + RW2(f)*N2(p))
C        where N is the glitch rate weighted number of ifgs, RW is relative
C        weights, p is pixel, f is frequency bin, and 1,2 refer to the input
C        skymaps.
C
C  Calculation of correction spectrum depends on combination:
C     Low frequency + Low frequency: offset correction spectrum
C     High frequency + High frequency: gain and offset correction spectrum
C     Low frequency + High frequency: offset correction spectrum
C     Others: gain and offset correction spectrum
C
C     For combining low frequency channels and/or scan modes calculate an
C     offset error correction spectrum:
C        CS(f) = SUM(p) {(Spec1(f,p) - Spec2(f,p)) * WP(f,p)} / SUM(p) WP(f,p)
C     Weights for applying correction spectrum:
C        WC1(f) = 1 / { 1 + [{NC1*D2} / {NC2*D1}] }  and
C        WC2 is gotten by swapping 1's and 2's.  NC1 is num ifgs in model
C        solution 1...
C
C     IMPORTANT!  Note that the formula in the "Combining FIRAS Skymaps -
C     A General Approach" are NOT symmetric in Left and Right.  In particular
C     the correction spectrum components.  Since the user may have entered
C     input skymaps in either order, we must determine and specify which side
C     "map" refers to in arrays with dimension (freq,map).  This is determined
C     by "in_skymap(1)(9:9) which is 'R' or 'L', and is used to set the
C     variable "Side".  Similarly "Fchan" gets set depending on which map is
C     High or Low.
C
C  Author:  Larry P. Rosen, Hughes STX, January 1994.  Important note and
C  relevant modifications added April 1994.
C  Modification: Larry P. Rosen, Hughes STX, July 1994.  Use glitch rate
C     weighted number of ifgs (which is stored in comb_num_ifgs.
C
C  * * * WARNING! * * * Major changes * * *
C     Modification: Larry P. Rosen, Hughes STX, October 1994.
C     Combining FIRAS Skymaps document has changed, and so must FMS.  Now, gain
C     and offset correction spectra are read from an FRD reference data set,
C     FEX_MCS_<chanscan>.<model>.  These spectra are created by FEF.  The
C     correction spectra are applied during the first and third order merges
C     only. Note, this routine is only called for first and third order merges.
C     FMS no longer writes out correction spectra files.
C
C     5 December 1994, Larry P. Rosen, Hughes STX.
C     This routine is now also called for second order merges.  FMS now
C     applies a gain and offset correction for second order merges.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Implicit None

C  Include files.

	Include		'($ssdef)'
C	Include		'CT$Library:CTUser.Inc'
	Include		'(FUT_Params)'		! Fac parameters
	Include		'(FMS_Msg)'
	Include		'(CCT_Get_Config)'
	Include		'CSDR$Library:CTParams.inc'

C  Passed parameters.

	Character*15	Arch_Cal		! Correction spectra archive
	Character*4	Specif (2)	! "chanscan" specification of input
	Character*12	Model_Ext		! Cal model soln extension
	Integer*4	Access_Time (2)		! Time for reference file
	Complex*8	Corr_Offset_Spec (257,2)  ! Gain correction spectrum
	Real*4		Corr_Gain_Spec (257,2)    ! Offset correction spectrum

C  Functions

	Integer*4	Lib$Get_Lun
	Integer*4	FMS_READ_SKYMAPS
	Integer*4	CCT_Open_Config
	integer*4	CCT_Get_Config_TOD
	Integer*4	CCT_Close_Config

C  Config stuff

	Character*14	Config_GMT_Start / '89324000001000' /
	Character*14	Config_GMT_Stop / '99365235958990' /
	Integer*4	Config_Start (2), Config_Stop (2)
	Character*48	Config_Name (1)
	Integer*4	Config_Size (1) / 772 /
	Integer*4	Config_Lun (1), Config_Index (1)
	Logical*1	Config_New_Segment (1)
	Record /Config_Status/ Config_Status(1)
	Integer*4	Config_Ref_Count
	Dictionary 'fex_mcs'
	Structure /Config_Record/
	   record /FEX_MCS/ CorrSpec_Rec
	EndStructure
	Record 	/Config_Record/ Config_Record

C  External

	External	CT_Connect_Write

C  Local

	Integer*4	Rstat			! Return status
	Integer*2	Freq			! Frequency bin
	Integer*2	Map			! Map number 1, 2

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FMS_Correction_Spectrum = %Loc (FMS_Normal)

C  Open, Read, and Close the fex_mcs files.

	Rstat = Lib$Get_Lun (Config_Lun(1))
	If (Rstat .NE. SS$_Normal) Then
	   FMS_Correction_Spectrum = %Loc (FMS_Abort)
	   Call Lib$Signal (FMS_Lunerr, %Val(1), %Val(Rstat))
	Else
	   Call CT_GMT_To_Binary (Config_GMT_Start, Config_Start)
	   Call CT_GMT_To_Binary (Config_GMT_Stop, Config_Stop)
	   Do Map = 1, 2
	      Config_Name(1)=Arch_Cal//'FEX_MCS_'//Specif(Map)//'.'//Model_Ext
	      Rstat = CCT_Open_Config ( Config_Start, Config_Stop, 1,
	1                               Config_Name, Config_Size, ' ', 1,
	2                               Config_Lun, Config_Index,
	3                               Config_Status, Config_Ref_Count )
	      If (.NOT. Rstat) Then
	         FMS_Correction_Spectrum = %Loc (FMS_Abort)
	         Call Lib$Signal ( FMS_Openerr, %Val(2), Config_Name(1),
	1                          %Val(Rstat) )
	      Else
	         Rstat = CCT_Get_Config_TOD ( Access_Time, 1, Config_Size,
	1                                     Config_Lun, Config_Index,
	2                                     Config_Record,Config_New_Segment,
	4                                     Config_Status )
	         If (.NOT. Rstat) Then
	            FMS_Correction_Spectrum = %Loc (FMS_Abort)
	            Call Lib$Signal (FMS_Readerr, %Val(2), Config_Name(1),
	1                            %Val(Rstat))
	         Else
	            Do Freq = 1, 257
	               Corr_Offset_Spec (Freq,Map) =
	1                 Config_Record.CorrSpec_Rec.Offset_Spec (Freq)
	               Corr_Gain_Spec (Freq,Map) =
	1                 Config_Record.CorrSpec_Rec.Gain_Spec (Freq)
	            EndDo
	            Rstat = CCT_Close_Config ( 1, Config_Lun, Config_Index )
	            If (.NOT. Rstat) Then
	               FMS_Correction_Spectrum = %Loc (FMS_Abort)
	               Call Lib$Signal (FMS_Closerr, %Val(2),Config_Name(1),
	1                               %Val(Rstat))
	            EndIf
	         EndIf
	      EndIf
	   EndDo
	EndIf
	Return
	End
