C******************************************************************************
	Integer*4 Function FPP_Midpoint_Time (Scan_Length, Scan_Speed,
	1   Sweeps_per_IFG, Transmit_Time, Sync, Start_Collect, End_Collect, 
	2   IFG_time, AncFilename)
C******************************************************************************
C FUNCTION: This routine computes the midpoint time of an interferogram.
C	It returns an odd number if successful and an even number if the
C	computation could not be made.
C
C CREATED BY: S Hilinski			DATE: 2-4-85
C
C argument		type	input/output	description
C --------              ----    ------------    -----------
C Scan_Length		I*2	input		1=long, 0=short.
C Scan_Speed		I*2	input		1=fast, 0=slow.
C Sweeps_per_IFG	I*2	input		Mirror sweeps in collection
C Transmit_Time(2)	I*4	input		Spacecraft time associated with
C						   Transmit_frame.
C Sync			L*1	input		MTM in Sync flag. True = in sync
C Start_Collect(2)	I*4(2)	in or out	Start of Collection time
C End_Collect(2)	I*4(2)	in or out	Start of Collection time
C IFG_Time(2)		I*4(2)	output		Midpoint time.
C AncFilename		C*64    input		Name of Anc file for open times.
C
C Include Files:  FPP_Msg.Txt, CCT_GET_CONFIG
C Reference File: FEX_MTMSweep.DAT
C******************************************************************************
C CHANGE LOG:
C
C	Version 4.2.0 01/21/89, SPR 2700, Shirley M. Read, STX
C	   During the Observatory I & T in the fall of 1988, a critical time tag
C          problem developed for the the FIRAS science archive. The computed
C	   midpoint of collect time caused the archiving of some of the science
C	   records with a bad primary key. FDQ and many subsequent Firas
C	   pipeline programs could not handle bad time tags.  Two SPRs were
C	   filed against the stripper and the above SPR against FDQ. The Firas
C	   Task wrote CCR 148 which was passed by the CCB on January 6, 1989.
C	   The FIRAS Stripper will now archive the science data by telemetry
C	   minor frame of transmit time. A new FIRAS Preprocessor will modify
C	   the science records with the computed midpoint of collect time, the
C	   MTM scan speed and length and a badtime flag for the computed time.
C	   The preprocessor will use this Midpoint_Time function instead of the
C	   stripper. The name has been updated to the FPP, in accord with the
C	   FIRAS software prefix convention. The values of the half sweep
C	   flyback times had changed since the FIRAS Stripper was written. The
C	   new values were put in an include file, FPP_MTM_Sweep.Txt, instead of
C	   being hard coded. The return status values were updated to reflect
C	   the types of error in the microprocessor information. The numbers
C	   correspond to the array postion in the fail-counter statistics. The
C	   extended precision arithmetic routines are now invoked as functions,
C	   instead of subroutines so that errors may be trapped and the
C	   computed IFG time may be flagged.
C
C	New Version, 27 March 1991, Larry P. Rosen, STX
C	   New Design and requirements lead to the following changes: frame
C	   counter checks are done in FPP_Collect_Time. If MTM is in sync,
C	   Start_Collect is known and input, so calculate End_Collect and
C	   midpoint of collect. If MTM is out of sync, then End_Collect is
C	   known and input. Include file MTM_SWEEP.TXT has been made into a
C	   reference data file since the originally used values have been found
C	   to be not quite right.  Ref file is FEX_MTMSWEEP.DAT and contains
C	   the total sweep + turnaround + flyback times and flyback times for
C	   scan modes: SS, LS, SF, LF respectively.  AncFilename contains times
C	   suitable for cct_open_config for reference data.
C******************************************************************************
	Implicit None

C  Include File:
	Include		'(FPP_Msg)'
	Include		'(CCT_GET_CONFIG)'

C  Passed Params:
	Integer*2	Scan_Length
	Integer*2	Scan_Speed
	Integer*2	Sweeps_per_IFG
	Integer*4	Transmit_Time(2)	! Transmission time VAX quadword
	Logical*1	Sync			! MTM sync
	Integer*4	Start_Collect(2)	!Start of IFG time in quadword
	Integer*4	End_Collect(2)	!End of IFG time in quadword
	Integer*4	IFG_time(2)	!Midpoint of IFG time in VAX quadword
	Character*64	AncFilename	!Name of anc file contains ref times.

C  Local Variables:
	Logical*1	First_Time /.True./	! first call to this routine.
	Integer*4	spos			! position in string
	Character*14	Ref_GMT_Start	! start gmt time for open ref data
	Character*14	Ref_GMT_Stop	! stop gmt time for open ref data
	Integer*4	Ref_Time_Start(2)	! start time for open ref data
	Integer*4	Ref_Time_Stop(2) 	! stop time for open ref data
	Integer*4	Status			! Function return status

C  Junk for CCT config routines:
	Integer*4	Lib$Get_Lun
	Integer*4	CCT_OPEN_CONFIG
	Integer*4       CCT_GET_CONFIG_TOD
	Integer*4	CCT_CLOSE_CONFIG
	Integer*4	Ndset / 1 /		! Number of datasets
	Character*32	Dset			! Names of datasets
	Data dset / 'CSDR$FIRAS_REF:FEX_MTMSWEEP' /
	Integer*4	Size / 128 /		! Dataset sizes
	character*1     access_mode / ' ' /     ! Access mode sequential-default
	Integer*4	Ncache / 1 /		! Number of caches- not used
	Integer*4	Con_lun			! Logical unit number for config
	Integer*4	Indx			! Initial cache pointers
	Record /Config_Status/ CStat
	Integer*4	Ref_count		! Reference counter for cache
	Logical*1       New_segment		! Flag for new segment accessed

C  More local variables:
	Integer*4	Sweeps			! I*4 version of sweeps per IFG
	Integer*4	Lastfly(2)		! time of last flyback to sub.
	Integer*2	Ti			! Time index
	Integer*4	TSweep			! Total time for 1 sweep f+r
	Integer*4	Total(2)		! Tsweep * # sweeps
	Integer*4	TrimTotal(2)		! Total - last flyback
	Integer*4	Two /2/
	Integer*4	Half			! Half of trim total as Longword
	Integer*4	Addend / 0 /		! Space for arithmetic
	Integer*4	Remainder
	Integer*4	Half_Duration(2)	! 1/2 IFG duration time.
	Dictionary	'FEX_MTMSWEEP'
	Record /FEX_MTMSWEEP/	FEX_Rec

C  Functions for extended precision arithmetic.

	Integer*4 Lib$Emul, Lib$Subx, Lib$Addx, Lib$Ediv

	FPP_Midpoint_Time = 1

C Initialize Midpoint time to transmit time.

	IFG_Time(1)=Transmit_Time(1)
	IFG_Time(2)=Transmit_Time(2)

C  First time, extract times from AncFilename for calling CCT_Open_Config.
C  Open and extract reference data FEX_MTMSWEEP.DAT.  This gets record
C  containing Total times ( = sweep + turnaround + flyback time) and
C  flyback times.

	If (First_Time) Then
	   First_Time = .False.
	   spos = Index (AncFilename,'/')
	   Ref_GMT_Start = AncFilename(spos+1:spos+14)
	   Ref_GMT_Stop = AncFilename(spos+16:spos+29)
	   Call CT_GMT_TO_BINARY (Ref_GMT_Start,Ref_Time_Start)
	   Call CT_GMT_TO_BINARY (Ref_GMT_Stop,Ref_Time_Stop)
	   Status = Lib$Get_Lun (Con_Lun)
	   If (.Not. Status) Then
	      Fpp_Midpoint_Time = 0
	      Call Lib$Signal(FPP_GetLunErr, %val(1), %val(Status))
	   Else
	      Status = CCT_OPEN_CONFIG ( Ref_Time_Start, Ref_Time_Stop, ndset,
	1        dset, size, access_mode, ncache, con_lun,indx,cstat,ref_count)
	      If ( .Not. Status ) Then
	         Call Lib$Signal (FPP_OPNCONFIGERR, %val(1), %val(status))
	         Fpp_Midpoint_Time = 0
	      Else
	         Status = CCT_GET_CONFIG_TOD ( Transmit_Time, ndset, size,
	1           con_lun, indx, Fex_Rec, new_segment, cstat)
	         If ( .Not. Status ) Then
	            Call Lib$Signal (FPP_GETCONFIGERR, %val(1), %val(status))
	            Fpp_Midpoint_Time = 0
	         EndIf
	         Status = CCT_CLOSE_CONFIG ( ndset, con_lun, indx )
	         If ( .Not. Status ) Then
	            Call Lib$Signal (FPP_CLSCONFIGERR, %val(1), %val(status))
	            Fpp_Midpoint_Time = 0
	         EndIf
	      EndIf
	   EndIf
	   If (.Not. FPP_Midpoint_Time) Return
	EndIf

C Compute 1/2 the IFG duration from the times in FEX_MTMSWEEP.DAT.
C 1/2 duration = (total time * # sweeps - 1 flyback time) /2

	Sweeps = Sweeps_per_IFG
	Ti = Scan_Length + 2 * Scan_Speed + 1
	Lastfly(1) = Fex_Rec.Flyback(Ti)
	TSweep = Fex_Rec.Total_Sweep_Flyback(Ti)
	Status = Lib$Emul (TSweep, Sweeps, Addend, Total)
	If (.Not. Status) Then
	   Call Lib$Signal (FPP_EmulErr,%val(1),%val(Status))
	   FPP_Midpoint_Time=0
	Else
	   Status = Lib$Subx (Total, LastFly, TrimTotal)
	   If (.Not. Status) Then
	      Call Lib$Signal (FPP_SubxErr,%val(1),%val(Status))
	      FPP_Midpoint_Time=0
	   Else
	      Status = Lib$EDiv ( Two, TrimTotal, Half, Remainder)
	      If (.Not. Status) Then
	         Call Lib$Signal (FPP_EDivErr,%val(1),%val(Status))
	         FPP_Midpoint_Time=0
	      Else
	         Half_Duration(1) = Half
	      EndIf
	   EndIf
	EndIf
	If (FPP_Midpoint_Time) Then
	   If (Sync) Then			! MTM in Sync

C The IFG midpoint is the start of collection time plus the half duration time.
C The End of Collect time is the midpoint time plus the half duration time.

	      Status = Lib$Addx (Start_Collect, Half_Duration, IFG_Time)
	      If (.Not. Status) Then
	         Call Lib$Signal (FPP_AddxErr,%val(1),%val(Status))
	         FPP_Midpoint_Time=0
	      Else
	         Status = Lib$Addx (IFG_Time, Half_Duration, End_Collect)
	         If (.Not. Status) Then
	            Call Lib$Signal (FPP_AddxErr,%val(1),%val(Status))
	            FPP_Midpoint_Time=0
	         Else
	            FPP_Midpoint_Time=1
	         EndIf
	      EndIf
	   Else				! MTM out of Sync

C The IFG midpoint is the End of Collection time minus the half duration time.
C The Start of Collect time is the midpoint time minus the half duration time.

	      Status = Lib$Subx (End_Collect, Half_Duration, IFG_Time)
	      If (.Not. Status) Then
	         Call Lib$Signal (FPP_SubxErr,%val(1),%val(Status))
	         FPP_Midpoint_Time=0
	      Else
	         Status = Lib$Subx (IFG_Time, Half_Duration, Start_Collect)
	         If (.Not. Status) Then
	            Call Lib$Signal (FPP_SubxErr,%val(1),%val(Status))
	            FPP_Midpoint_Time = 0
	         Else
	            FPP_Midpoint_Time = 1
	         EndIf
	      EndIf
	   EndIf
	EndIf

	Return
	End
