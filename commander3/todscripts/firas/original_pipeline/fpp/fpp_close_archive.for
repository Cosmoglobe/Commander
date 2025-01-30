C-------------------------------------------------------------------------------
	Integer*4 Function FPP_Close_Archive ( Chan_Num, Report, Lun_Rpt,
	1	  Lun_Anc, Lun_In_Sci, Lun_Out_Sci, Lun_Fake, Lun_Gain )
C-------------------------------------------------------------------------------
C	Purpose: To close the Word 31 ancillary file, the science channel file,
C		and the FEX_Gain and FEX_Fakeit reference files.
C
C	Author: Shirley M. Read
C		STX, February, 1989
C
CH	Change Log:
CH
CH	SPR 4030, Fill FPP run flag into segment catalog.
CH	    R. Kummerer, June 16, 1989.
CH
CH	Version 4.4.1 07/22/89, SER 4168, R. Kummerer, STX
CH	    There have been problems during test operations for FIRAS processing
CH	    due to the required clean-up of the archives after an FPP or FDQ
CH	    abort. The FPR tracking system compounds the problems. Files with
CH	    non-matching version numbers seem often to result from improper
CH	    clean-up. Bad record times cause SEGCTL to abort and mess up the
CH	    tracking system. It was decided to change the modify of the science
CH	    records in FPP and FDQ to a simple COBETRIEVE read of the existing
CH	    records from a dataset and write a modifed dataset with the same
CH	    information which was entered on the modify. Two new science data
CH	    sets will be required: a science dataset of raw science data plus
CH	    FPP input and a science dataset with FPP science data plus FDQ
CH	    input. These datasets will be FPP_SDF_xx, where xx is the channel
CH	    id (RH, RL, LH or LL) and FDQ_SDF_xx, where xx is the channel id.
CH	    The new datasets must be opened and processed in FPP and FDQ. 
CH
CH	Version 4.4.1 08/21/89 SER 4210, R. Kummerer, STX
CH		Prevent overlaps in raw science segments.
CH
CH	New Version, New requirements/design, STX, Larry P. Rosen, 5 April 1991
CH		Remove tracking. Close FEX reference files.
C------------------------------------------------------------------------------
C	Input Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	  Chan_Num        I*2 		Current channel number
C	  Report          L*1           Flag for report
C         Lun_Rpt 	  I*4           Logical unit for report file  
C	  Lun_Anc         I*4           Unit for Word 31 ancillary data 
C	  Lun_In_Sci      I*4           Unit for input science channel data
C	  Lun_Out_Sci     I*4           Unit for output science channel data
C	  Lun_Fake        I*4           Unit for fakeit reference file
C	  Lun_Gain        I*4           Unit for gain reference file
C
C	Subroutines Called:
C	  Lib$Signal
C	  Cct_Get_FileSpec_Fields
C	  Cct_Query_Catalog
C
C	Include Files:
C	  CT$Library:CTUser.Inc
C	  FPP_Msg.Txt
C         $SSDef
C
C	Processing Method: PDL for FPP_Close_Archive
C
C	Call CT_Close_Arcv to close the NFS_ANC file.
C	If there is an error, then
C	   Signal an error and set an error status for the return.
C	Endif
C
C	Call CT_Close_Arcv to close the NFS_SDF channel file.
C	If there is an error, then
C	   Signal an error and set an error status for the return.
C	Endif
C
C	Call CT_Close_Arcv to close the FPP_SDF channel file.
C	If there is an error, then
C	   Signal an error and set an error status for the return.
C	Endif
C
C	Call CT_Close_Arcv to close the FEX_Fakeit file.
C	If there is an error, then
C	   Signal an error and set an error status for the return.
C	Endif
C
C	Call CT_Close_Arcv to close the FEX_Gain file.
C	If there is an error, then
C	   Signal an error and set an error status for the return.
C	Endif
C
C	Return with normal or error status. 
C
C------------------------------------------------------------------------------
	Implicit None

C  Passed Parameters.
	Integer*2     Chan_Num          ! Current channel number
	Logical*1     Report            ! Flag to enable writing a report
	Integer*4     Lun_Rpt           ! Logical unit for report file  
	Integer*4     Lun_Anc           ! Logical unit for Word 31 NFS_ANC data
	Integer*4     Lun_In_Sci        ! Logical unit for input science data
	Integer*4     Lun_Out_Sci       ! Logical unit for output science data
	Integer*4     Lun_Fake		! Unit for fakeit reference file
	Integer*4     Lun_Gain		! Unit for gain reference file

C  Include files.
	Include      'CT$Library:CTUser.Inc'
	Include      '(FPP_Msg)'

C  Local Declarations.
	Integer*2    CT_Stat(20)        ! COBETRIEVE return status

C  Set the function status to Normal.

	FPP_Close_Archive = %loc(FPP_Normal)

C  Close the matching set of files.

	Call CT_Close_Arcv (, Lun_Anc, CT_Stat)  ! NFS_ANC ancillary data
	If (CT_Stat(1) .NE. CTP_Normal) Then
	   FPP_CLose_Archive = %loc(FPP_Aberr)
	   Call Lib$Signal (FPP_CTClosErr, %val(1), %val(CT_Stat(1)))
	EndIf

	Call CT_Close_Arcv (, Lun_In_Sci, CT_Stat)  ! NFS_SDF science data
	If (CT_Stat(1) .NE. CTP_Normal) Then
	   FPP_CLose_Archive = %loc(FPP_Aberr)
	   Call Lib$Signal (FPP_CTClosErr, %val(1), %val(CT_Stat(1)))
	EndIf

	Call CT_Close_Arcv (, Lun_Out_Sci, CT_Stat)  ! FPP_SDF science data
	If (CT_Stat(1) .NE. CTP_Normal) Then
	   FPP_CLose_Archive = %loc(FPP_Aberr)
	   Call Lib$Signal (FPP_CTClosErr, %val(1), %val(CT_Stat(1)))
	EndIf

	Call CT_Close_Arcv (, Lun_Fake, CT_Stat)  ! Fex_Fakeit data
	If (CT_Stat(1) .NE. CTP_Normal) Then
	   FPP_CLose_Archive = %loc(FPP_Aberr)
	   Call Lib$Signal (FPP_CTClosErr, %val(1), %val(CT_Stat(1)))
	EndIf

	Call CT_Close_Arcv (, Lun_Gain, CT_Stat)  ! Fex_Gain data
	If (CT_Stat(1) .NE. CTP_Normal) Then
	   FPP_CLose_Archive = %loc(FPP_Aberr)
	   Call Lib$Signal (FPP_CTClosErr, %val(1), %val(CT_Stat(1)))
	EndIf

 	Return
	End
