C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                          FMS - FIRAS  Merge  Skymaps
C
C  Purpose and  Top  Level  Description.
C
C  Combine the skymaps produced by the FIRAS pipeline.  The pipeline produces
C  a skymap for 12 channel and scan mode combinations of the instrument with
C  one spectrum per pixel.  FMS combines skymaps across channel and scan mode,
C  correcting the spectra for calibration gain and offset errors in the
C  process.  FMS is run once for each combination, initially taking FCS
C  skymaps as input, and later taking FMS products as input.
C
C  This file contains the top level, main calling program of FMS.
C
C  References (these documents should be in FIRUSER:[FIRDOC.FMS]):
C  1. "Combining FIRAS Skymaps - A General Approach" by Gene Eplee and
C         Dale Fixsen for explanations for algorithms used.
C  2. "FMS Requirements" by Gene Eplee.
C  3. "FMS Design Walkthrough" by Larry Rosen.
C
C  SER 11418    : Need Program to Merge Skymaps...
C
C  Author:  Larry P. Rosen, Hughes STX, January 1993.
C
C    Modified: July 1994, Larry P. Rosen, Hughes STX.
C       NIFGs for weighting is now glitch rate weighted.
C    10/25/94, Larry P. Rosen.  Combining FIRAS Skymaps document has changed,
C       and so must FMS.  Now, gain and offset correction spectra are read
C       from an FRD reference data set, FRD_MCS_<chanscan>.<model>.
C       These spectra are created by FEF.  The correction spectra are applied
C       during the first and third order merges only.  At second order only the
C       usual glitch rate weighted number of ifgs and c-vectors are used to
C       weight the spectra. SER 11977.
C    12/5/94, Larry P. Rosen.  Nix the previous comment about second order.
C       Ken Jensen demonstrated that we do need an offset and gain correction
C       in second order.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Program  FMS

	Implicit None

C  Include files.

	Include		'($SSDef)'
	Include		'(FUT_Error)'
	Include		'(FUT_Params)'		! Fac parameters
	Include		'(FMS_MSG)'
	Include		'CSDR$Library:CTParams.Inc'

C  External

	External	FUT_Error

C  Functions

	Integer*4	CUT_Register_Version,  CUT_Display_Banner
	Integer*4	FMS_PARSE
	Integer*4	FMS_INIT_REPORT
	Integer*4	FMS_VERIFY_MERGE
	Integer*4	FMS_READ_C_D_VECTORS
	Integer*4	FMS_OPEN_SKYMAPS
	Integer*4	FMS_COMPUTE_WEIGHTS
	Integer*4	FMS_COMBINE_C_D_VECTORS
	Integer*4	FMS_CORRECTION_SPECTRUM
	Integer*4	FMS_READ_SKYMAPS
	Integer*4	FMS_SPECTRA_VARIANCES
	Integer*4	FMS_COMBINE_OTHER
	Integer*4	FMS_TRACKING_INFO
	Integer*4	FMS_WRITE_SKYMAP
	Integer*4	FMS_CLOSE

C  Local  -  note: All fields dimensioned 2 refer to each input skymap.

	Character*6	Version
	Parameter	(Version = '12.6  ')
	Integer*4	Rstat				! Return status
	Integer*4	Status, Err /2/, Ok /0/		! Processing status
	Character*80	Script			! Script file name
	Integer*2	Slen			! Length of script name string
	Character*9	Combo			! Combination type
	Character*6	Order			! Command combination order
	Character*12	Model_Ext		! Cal model soln extension
	Character*8	Weight			! C_VECTOR or D_VECTOR
	Integer*2	Trans_Freq	! Transition Low-High Frequency in icm
	Logical*1	Report			! Whether or not to report
	Character*45	Reportfile		! Filename for report
	Integer*2	Rlen			! Length of report file name
	Logical*1	Rep_Default		! Flag for default report name
	Character*4	Specif (2)		! "chanscan" of input skymaps
	Character*33	In_Skymap (2)		! Input skymaps
	Character*20	File_Ext		! Output file extension
	Character*79	Cmdline (3)		! Command line with defaults
	Integer*2	Ct_Stat (20)		! Cobetrieve status
	Integer*4	Lun_Rep			! Report logical unit #
	Character*14	Arch_In	/'CSDR$FIRAS_IN:'/	! Input data archive
	Character*15	Arch_Out /'CSDR$FIRAS_OUT:'/	! Output data archive
	Character*15	Arch_Ref /'CSDR$FIRAS_REF:'/	! C or D Vector archive
	Character*15	Arch_Cal /'CSDR$FIRAS_CAL:'/	! Correction spect arch
	Character*14	Tstart			! Start time in GMT
	Character*14	Tstop			! Stop time in GMT
	Integer*4	Current_Time(2)		! Current system adt time
	Character*4	Specif_Out	! "chanscan" specification of output
	Integer*2	Order_Num		! Command order as an integer
	Dictionary	'FEX_CVS'
	Record /FEX_CVS/	C_Vector_Rec (2)
	Dictionary	'FEX_VAR'
	Record /FEX_VAR/	D_Vector_Rec (2)
	Real*8		Weights (257,2)		! C or D vector
	Real*8		Relative_Weights (257,2)
	Character*40	Model			! Calibration Solution
	Real*8		No_Freq_Weights (2)	! Freq indep. weight of skymaps
	Real*8		C_Vector (257,2)	! C-vector for skymap 1 and 2
	Real*8		D_Vector (257,2)	! D-vector for skymap 1 and 2
	Integer*4	Access_Time (2)		! For reference data
	Integer*4	Lun_Out, Lun_In (2)	! Output and input logical unit
	Integer*4	Pixel			! Pixel number
	Dictionary	'FMS_SKY'
	Record /FMS_SKY/	In_Recs (2)	! 1 record from each skymap
	Record /FMS_SKY/	Out_Rec		! Output record for pixel
	Byte	Blank_Record (Fac_Coad_Spec_Size) / Fac_Coad_Spec_Size * 0 /
		! A blank FMS record.  fac_coad_spec_size = 7680 (FUT_Params)
	
	Complex*8	Corr_Offset_Spec (257,2) / 514 * (0.,0.) /
			! Offset correction spectrum
	Real*4		Corr_Gain_Spec (257,2) / 514 * 1. /
			! Gain correction spectrum
	Integer*2	Nu_Last_Low	! Highest index for low frequency band
	Character*1	Combo_Type		! Combination type.
			! Z=Zeroth, L=Low+Low, H=High+High, F=High+Low, M=Other
	Logical*1	Rec_Found (2)		! Flag record found in maps 1,2
	Real*8		Nifgs (2)	! Glitch rate weighted number of ifgs
	Integer*2	Numpix (3) /0,0,0/	! Number of pixels with data
						! 1=map1, 2=map2, 3=merged map

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	Status = Ok

C  Print software version number and display banner.

	Rstat = CUT_Register_Version (Version)
	Rstat = CUT_Display_Banner (6, 80, 'FIRAS  FMS_Merge_Skymaps')

C  Call FMS_PARSE to parse the command line, to open, read, and close script
C  file.  Returns qualifiers and input skymap file names.

	Rstat = FMS_PARSE ( Script, Slen, Combo, Order, Model_Ext, Weight,
	1                   Trans_Freq, Report, Reportfile, Rlen, Rep_Default,
	2                   Specif, In_Skymap, File_Ext, Cmdline )

	If ( Rstat .NE. %Loc (FMS_Normal) ) Status = Err

C  Call CT_INIT to initialize Cobetrieve.

	If (Status .EQ. Ok) Then
	   Call CT_Init ( Ct_Stat )
	   If (Ct_Stat(1) .NE. Ctp_Normal) Then
	      Call Lib$Signal ( FMS_Ctinit, %Val(1), %Val (Ct_Stat(1)) )
	      Status = Err
	   EndIf
	EndIf

C  If (Report) Call FMS_INIT_REPORT to open and initialize the report file.

	If ((Status .EQ. Ok) .AND. (Report)) Then

	   Rstat = FMS_INIT_REPORT ( Lun_Rep, Reportfile, Rep_Default, Rlen,
	1                            Version, Cmdline, Script, Slen,
	2                            Trans_Freq, In_Skymap, File_Ext, Arch_In,
	3                            Arch_Out, Arch_Ref, Tstart, Tstop,
	4                            Current_Time )

	   If ( Rstat .NE. %Loc (FMS_Normal) ) Then
	      Status = Err
	   Else
	      Call Lib$Establish ( FUT_Error )
	   EndIf
	EndIf

C  Call FMS_VERIFY_MERGE to check whether the combination is valid.

	If (Status .EQ. Ok) Then

	   Rstat = FMS_VERIFY_MERGE ( Combo, Order, Weight, Report, Lun_Rep,
	1                             Rep_Default, Reportfile, Specif,
	2                             Specif_Out, Tstart, Tstop, Current_Time,
	3                             Combo_Type, Order_Num )

	   If ( Rstat .NE. %Loc (FMS_Normal) ) Status = Err
	EndIf

C  Call FMS_READ_C_D_VECTORS to open, read, and close C and D vector files.

	If (Status .EQ. Ok) Then

	   Rstat = FMS_READ_C_D_VECTORS ( Model_Ext, In_Skymap, Specif,
	1                                 Report, Lun_Rep, Arch_Ref, Arch_In,
	2                                 C_Vector_Rec, D_Vector_Rec,
	3                                 Access_Time )

	   If ( Rstat .NE. %Loc (FMS_Normal) ) Status = Err
	EndIf

C  Call FMS_OPEN_SKYMAPS to open input and output skymap files.

	If (Status .EQ. Ok) Then

	   Rstat = FMS_OPEN_SKYMAPS ( Lun_In, Lun_Out, Arch_In, Arch_Out,
	1                             In_Skymap, Specif_Out, File_Ext, Report,
	2                             Lun_Rep )

	   If ( Rstat .NE. %Loc (FMS_Normal) ) Status = Err
	EndIf

C  Call FMS_COMPUTE_WEIGHTS to compute relative weights for combination.

	If (Status .EQ. Ok) Then

	   Rstat = FMS_COMPUTE_WEIGHTS ( Lun_In, Trans_Freq, Specif, Weight,
	1                                C_Vector_Rec, D_Vector_Rec, Weights,
	2                                Relative_Weights, Model, Nu_Last_Low,
	3                                No_Freq_Weights, C_Vector, D_Vector,
	4                                Combo_Type )

C  Call FMS_COMBINE_C_D_VECTORS to combine C and D vectors, and write to
C  output file.

	   Rstat = FMS_COMBINE_C_D_VECTORS ( Weight, C_Vector, D_Vector,
	1                                    C_Vector_Rec, D_Vector_Rec,
	2                                    Relative_Weights,
	3                                    No_Freq_Weights, Specif_Out,
	4                                    Arch_Out, Model_Ext, Report,
	5                                    Lun_Rep, Nu_Last_Low, Specif,
	6                                    Combo_Type, Lun_In )

	   If ( Rstat .NE. %Loc (FMS_Normal) ) Status = Err
	EndIf

C  If Order is not ZEROTH (ie. FIRST, SECOND, or THIRD) call
C  FMS_CORRECTION_SPECTRUM to read the correction spectra.

	If ( (Order .NE. 'ZEROTH') .AND. (Status .EQ. Ok) ) Then

	   Rstat = FMS_CORRECTION_SPECTRUM ( Arch_Cal, Specif, Model_Ext,
	1                                    Access_Time, Corr_Offset_Spec,
	2                                    Corr_Gain_Spec )

	   If ( Rstat .NE. %Loc (FMS_Normal) ) Status = Err
	EndIf

C  Do for each pixel.  FAC_Num_Pixels = 6144 (FUT_Params)

	Pixel = 0
	Do While ((Pixel .LT. FAC_Num_Pixels) .AND. (Status .EQ. Ok))

C  Initialize an empty input and output records.  (copy Blank_Record)

	   Call Lib$MovC3 ( Fac_Coad_Spec_Size, Blank_Record, Out_Rec )
	   Call Lib$MovC3 ( Fac_Coad_Spec_Size, Blank_Record, In_Recs (1) )
	   Call Lib$MovC3 ( Fac_Coad_Spec_Size, Blank_Record, In_Recs (2) )

C  Call FMS_READ_SKYMAPS to read a record from each skymap for this pixel.

	   Rstat = FMS_READ_SKYMAPS ( Lun_In, Pixel, In_Recs, Rec_Found )

	   If ( Rstat .NE. %Loc (FMS_Normal) ) Then
	      Status = Err
	   ElseIf ( Rec_Found (1) .OR. Rec_Found (2) ) Then

C  Call FMS_SPECTRA_VARIANCES to combine spectra and variances.

	      Rstat = FMS_SPECTRA_VARIANCES ( In_Skymap, Weights, Specif,
	1                                     Order, Corr_Offset_Spec,
	2                                     Corr_Gain_Spec, Nu_Last_Low,
	3                                     Relative_Weights, Combo_Type,
	4                                     In_Recs, Out_Rec, Nifgs,
	5                                     No_Freq_Weights, C_Vector )

C  Call FMS_COMBINE_OTHER to combine all fields except spectra and variances.

	      Rstat = FMS_COMBINE_OTHER ( Rec_Found, In_Recs, Out_Rec, Specif,
	1                                 Specif_Out, Nifgs, No_Freq_Weights,
	2                                 Numpix )

C  Call FMS_TRACKING_INFO to store tracking information in the output record.

	      Rstat = FMS_TRACKING_INFO ( In_Recs, Out_Rec, Specif, Nifgs,
	1                                 Order_Num, In_Skymap )

C  Call FMS_WRITE_SKYMAP to write the merged spectral record.

	      Rstat = FMS_WRITE_SKYMAP ( Lun_Out, Out_Rec )

	      If ( Rstat .NE. %Loc (FMS_Normal) ) Status = Err
	   EndIf

	   Pixel = Pixel + 1
	EndDo

C  Call FMS_CLOSE to close input and output skymap files, write statuses to
C  and close the report.

	If (Status .EQ. Ok) Then

	   Rstat = FMS_CLOSE ( Lun_In, Lun_Out, Report, Lun_Rep, Numpix,
	1                      In_Skymap, Specif_Out, File_Ext )

	   If ( Rstat .NE. %Loc (FMS_Normal) ) Status = Err
	EndIf

C Exit

	If (Status .EQ. Ok) Then
	   Call Exit (SS$_Normal)
	Else
	   Call Exit (SS$_Abort)
	EndIf
	Stop
	End
