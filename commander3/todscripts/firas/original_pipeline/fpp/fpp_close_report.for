C-------------------------------------------------------------------------------
	Integer*4 Function FPP_Close_Report ( Maxn, Num_Chan_Tbl, Proc_Chan_Tbl,
	1  Num_Days, Out_Files, Num_Rec, Num_Good, Num_Bad, BadFlg, Lun_Rpt,
	2  Status )
C-------------------------------------------------------------------------------
C	Purpose: To write summary information into the report file,
C		 primarily, the statistics for the run.
C
C	Author: Shirley M. Read
C		STX, January, 1989
C
C	Invocation: Status = FPP_Close_Report ( Maxn, Num_Chan_Tbl,
C		    Proc_Chan_Tbl, Num_Days, Out_Files, Num_Rec, Num_Good,
C		    Num_Bad, BadFlg, Lun_Rpt, Status )
C
CH	Change Log:
CH
CH	Version 4.4.1 07/22/89, SER 4168, R. Kummerer, STX
CH	There have been problems during test operations for FIRAS processing
CH	due to the required clean-up of the archives after an FPP or FDQ abort.
CH	The FPR tracking system compounds the problems. Files with non-matching
CH	version numbers seem often to result from improper clean-up. Bad record
CH	times cause SEGCTL to abort and mess up the tracking system. It was
CH	decided to change the modify of the science records in FPP and FDQ to a
CH	simple COBETRIEVE read of the existing records from a dataset and write
CH	a modifed dataset with the same information which was entered on the
CH	modify. Two new science data sets will be required: a science dataset
CH	of raw science data plus FPP input and a science dataset with FPP
CH	science data plus FDQ input. These datasets will be FPP_SDF_xx, where
CH	xx is the channel id (RH, RL, LH or LL) and FDQ_SDF_xx, where xx is the
CH	channel id. The new datasets must be opened and processed in FPP and FDQ
CH
CH	Version 4.4.1 08/21/89 SER 4210, R. Kummerer, STX
CH		Prevent overlaps in raw science segments.
CH
CH	Version 4.4.4 11/29/89, SPR 5165, R. Kummerer STX
CH		Fails to skip missing segments.
CH 
CH	New Version, New Requirements, Larry P. Rosen, STX, 5 April 1991
CH	   New flags for BadTime are listed in new chart.
CH
CH	Add BadFlag #18, for bad number of sweeps. Larry P. Rosen, 
CH	Hughes STX, 17 October 1991, SER 9167.
C-------------------------------------------------------------------------------
C	Input Parameters:
C	  Name		  Type		  Description
C 	  ----------------------------------------------------------------------
C	  Maxn		   I*2		  Number of files
C	  Num_Chan_Tbl(Maxn) I*2            Number of channel segments processed
C	  Proc_Chan(Maxn,4)  I*2            Channel segmetns processed
C	  Num_Days         I*4            Number of files for each channel
C	  Out_Files(Maxn,4) C*64            Names of files that were created
C	  Num_Rec(Maxn,4)   I*4             Number of records in each file
C	  Num_Good(Maxn,4)  I*4             Number of good records in each file
C	  Num_Bad(Maxn,4)   I*4             Number of bad records in each file
C	  BadFlg(18,Maxn,4)  I*4        Counters for failures on time tags
C         Lun_Rpt 	  I*4             Logical unit for report file  
C	  Status          I*4             Processing status for the run
C
C	Subroutines Called:
C	  Lib$Signal
C	  Lib$Locc
C
C	Include Files:
C	  FPP_Msg.Txt
C         $SSDef
C
C	Processing Method: PDL for FPP_Close_Report
C
C	If the processing status is normal, then
C
C	   Initialize a summary table for the production run.
C
C	   Do for Each Channel Processed for the Run
C
C	      Write the dataset name, the channel number and the number 
C                of files processed in the report file.
C
C	      Write the filenames, the number of catalog records for the file
C                and the number of records actually processed.
C
C	      If the number processed does not agree with the number of records
C                in the catalog, write a warning message.
C
C	      Write the counter array for failures on time tags.
C
C	   Enddo
C
C	   Write a success message for the completion of the run.
C	
C	Elseif the processing status is abort due to error, then
C
C	   Write an abort message in the report file.
C
C	Endif
C
C	Close the report file.
C
C	Return
C------------------------------------------------------------------------------
	Implicit None

!  Passed Parameters.

	Integer*2	Maxn			! Number of files
	Integer*2       Num_Chan_Tbl(Maxn)        ! Number of channel segments
						! processed
	Integer*2	Proc_Chan_Tbl(Maxn,4)     ! Channel segments processed
	Integer*4	Num_Days         ! Number of files for each channel
	Character*64	Out_Files(Maxn,4) ! Names of files to be processed
	Integer*4	Num_Rec(Maxn,4)   ! Number of records in each file
	Integer*4	Num_Good(Maxn,4)  ! Number of good records in each file
	Integer*4	Num_Bad(Maxn,4)   ! Number of bad records in each file
	Integer*4       BadFlg(18,Maxn,4)  ! Counters for failures on time tag
        Integer*4       Lun_Rpt         ! Logical unit for report file  
	Integer*4	Status          ! Processing status for the run

!	Include files.

	INCLUDE		'(fut_error)'
	Include      	'(FPP_Msg)'
	Include      	'($SSDef)'

!	Functions

	Integer*4 	Lib$Locc        ! Locate char in string

!	Local Declarations.

	Integer*4    	Iostatus	! Return IO status
	Character*1	Dashes(79)  / 79 * '-' /
	Character*79    Dash 
	Equivalence     ( Dashes(1), Dash )
	Integer*2	Ix, Jx		! Indices
	Integer*2       Chan		! Channel index
	Integer*2       Channel         ! Channel number
	Integer*4       Pos1, Pos2      ! String positions
	Integer*4	Tot		! Total record count
	Integer*4	BadFlgTab(18,4)	! table of bad flags for 4 chan
	Integer*4	success / 1 /, err / 2 /  ! values for status
! Bad Time Flags meanings
	Character*29	Flag(18) / 'Ancillary GAP Before IFG     ' ,
	1		           'Transmit Time Before Collect ' ,
	1		           'Transmit Count < 0           ' ,
	1		           'Collect Count < 0            ' ,
	1		           'Transmit = Collect Count     ' ,
	1		           'Midpoint Time > Transmit time' ,
	1		           'Midpoint < Previous Transmit ' ,
	1		           'Mode Consistency Check Failed' ,
	1		           'Ancillary GAP During IFG     ' ,
	1		           'Computational Error          ' ,
	1		           'Gain Unknown                 ' ,
	1		           'FakeIt Unknown               ' ,
	1		           'Science GAP                  ' ,
	1		           'Housekeeping Data GAP        ' ,
	1		           'HKP Telemetry Quality Bad    ' ,
	1		           'Non-Science Telemetry Format ' ,
	1		           'Bad Number of Adds per Group ' ,
	1		           'Bad Number of MTM Sweeps' /

!  Set the function status to Normal.

	FPP_Close_Report = %loc(FPP_Normal)

!  Write the summary title for the run statistics in the report file.

	Write (Unit=Lun_Rpt, FMT=100, Iostat=Iostatus) Dash(1:79)
	If (Iostatus .NE. 0) Then
	   FPP_Close_Report = %loc(FPP_Aberr)
	   Call Lib$Signal(FPP_WritErr, %val(1),%val(Iostatus))
	   status = err
	EndIf

	Do Ix = 1, Num_Days

	   Write (Unit=Lun_Rpt, FMT=200,Iostat=Iostatus) Ix, Num_Days,
	1     Num_Chan_Tbl(Ix)
	   If (Iostatus .NE. 0) Then
	      FPP_Close_Report = %loc(FPP_Aberr)
	      Call Lib$Signal(FPP_WritErr, %val(1),%val(Iostatus))
	      status = err
	   EndIf
	   Write (Unit=Lun_Rpt, FMT=300, Iostat=Iostatus)
	   If (Iostatus .NE. 0) Then
	      FPP_Close_Report = %loc(FPP_Aberr)
	      Call Lib$Signal(FPP_WritErr, %val(1),%val(Iostatus))
	      status = err
	   EndIf

	   Do Chan = 1, Num_Chan_Tbl(Ix)

	      Channel = Proc_Chan_Tbl(Ix,Chan)
	      Pos1 = Lib$Locc(':',Out_Files(Ix,Channel))
	      Pos1 = Pos1 + 1
	      Pos2 = Pos1 + 23
	      Write (Unit=Lun_Rpt, FMT=400, Iostat=Iostatus) 
	1	Out_Files(Ix,Channel)(Pos1:Pos2), Num_Rec(Ix,Channel), 
	2	Num_Good(Ix, Channel), Num_Bad(Ix,Channel)
	      If (Iostatus .NE. 0) Then
	         FPP_Close_Report = %loc(FPP_Aberr)
	         Call Lib$Signal(FPP_WritErr, %val(1),%val(Iostatus))
	         status = err
	      EndIf
	      Tot = Num_Good(Ix,Channel) + Num_Bad(Ix,Channel)
	      If ( Num_Rec(Ix,Channel) .NE. Tot ) Then
	         Write ( Unit=Lun_Rpt, FMT=500, Iostat=Iostatus)
	         If (Iostatus .NE. 0) Then
	            FPP_Close_Report = %loc(FPP_Aberr)
	            Call Lib$Signal(FPP_WritErr, %val(1),%val(Iostatus))
	            status = err
	         EndIf
	      EndIf
	      Do Jx = 1,18
	         BadFlgTab(Jx,Channel) = BadFlg(Jx,Ix,Channel)
	      EndDo
	   Enddo 			! 1 to Num_Chan_Tbl
	   Write (Unit=Lun_Rpt, FMT=600, Iostat=Iostatus) Dash(1:74)
	   Do Jx = 1,18
	      Write (Unit=Lun_Rpt, FMT=700) Jx, Flag(Jx),
	1        (BadFlgTab(Jx,Chan),Chan=1,4)
	   EndDo				! each flag
	EndDo					! each day of data

C  Signal the processing status.

	IF (status .EQ. success) THEN
	  CALL lib$signal (fpp_normal)
	ELSE
	  CALL lib$signal (fpp_aberr)
	ENDIF

	Close (Unit=Lun_Rpt, Iostat=Iostatus)
	FUT_Report_Lun = 0

	If (Iostatus .NE. 0) Then
	  FPP_Close_Report = %loc(FPP_Aberr)
	  Call LIB$Signal(FPP_ClosErr, %val(1),%val(Iostatus))
	Endif

 100	Format (1x,/,20x,'Run Statistics for FPP_PRE_PROCESSOR',/,A)
 200	Format (1x,/,5x,'Segment ',I5,2x,'out of ',I5,6x,'Number of Channels:',
	1	I2)
 300    Format (1x,/,4x,'Science Filename',10x,'Total Cat. Rec.',3x,
	1	'Valid Time Rec.',3x,'Bad Time Rec.',/)
 400    Format (2x, A23, 9x, I5, 14x, I5, 12x, I5)
 500    Format (1x,/,5x,'Warning: Number of Valid Time Records + ',
	1       'Number of Bad Time Records',/,
	2	10x,'is not equal to Total Number of Records in Catalog.')
 600	Format (1x,/,4x,'Table of Number of Records with Each Bad Time Flag ',
	1	'for Each Channel:',/,/,8x,'Bad Time Flag',21x,'RH',6x,'RL',
	2	6x,'LH',6x,'LL',/,5x,A)
 700	Format (2x,I2,1x,A29,4(3x,I5))
	Return
	End
