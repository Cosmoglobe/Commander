C----------------------------------------------------------------------------
	Integer*4 Function Fpp_Cat_Info (Time_Present, Gmt_Start, Gmt_Stop,
	1	File_Present, File_Name, Proc_Chan, Num_Chan, Report, Lun_Rpt,
	2	Maxn, Num_Chan_Tbl, Proc_Chan_Tbl, Num_Days, In_Files,Out_Files,
	3	Num_Rec, Check_Jstart, Check_Jstop, Arch_In, Arch_Out  )
C-----------------------------------------------------------------------------
C	Purpose: To query and to save The information from the catolog of
C	         NFS_SDF AND NFS_ANC.
C
C	Author:  Quoc C Chung / STX, February, 1989
C
C	Invocation:  Rstatus = Fpp_Cat_Info(Time_Present, GMT_start, GMT_stop,
C	1	File_Present, File_Name, Proc_Chan, Num_Chan, Report, Lun_rpt,
C	2	Maxn, Num_Chan_Tbl, Proc_Chan_Tbl, Num_Days, In_Files,Out_Files,
C	3	Num_rec, Check_JStart, Check_JStop, Arch_In, Arch_Out  )
C
CH	Change Log:
CH		SPR 3993, Fill FILENAMES using file counter variable to ensure
CH		raw science and ancillary are grouped together properly when
CH		some raw science segments are marked as modified by a previous
CH		run of FPP.  S. Read, R. Kummerer, June 9, 1989. STX
CH
CH		SPR 4006, Fails to ensure that raw science always matches
CH		ancillary segments.  R. Kummerer, June 13, 1989. STX
CH
CH		Version 4.4.1 07/22/89, SER 4168, R. Kummerer, STX
CH		There have been problems during test operations for FIRAS
CH		processing due to the required clean-up of the archives after an
CH		FPP or FDQ abort. The FPR tracking system compounds the problems
CH		Files with non-matching version numbers seem often to result
CH		from improper clean-up. Bad record times cause SEGCTL to abort
CH		and mess up the tracking system. It was decided to change the
CH		modify of the science records in FPP and FDQ to a simple
CH		COBETRIEVE read of the existing records from a dataset and write
CH		a modifed dataset with the same information which was entered on
CH		the modify. Two new science data sets will be required: a science
CH		dataset of raw science data plus FPP input and a science dataset
CH		with FPP science data plus FDQ input. These datasets will be
CH		FPP_SDF_xx, where xx is the channel id (RH, RL, LH or LL) and
CH		FDQ_SDF_xx, where xx is the channel id. The new datasets must be
CH		opened and processed in FPP and FDQ.
CH
CH		Version 4.4.1 08/16/89 SPR 4268, R. Kummerer, STX
CH			Correct CCT_Query_TOD_Catalog NUM_SCI_SEG argument
CH			declaration.  Changed from I*4 to I*2.
CH
CH		Version 4.4.1 08/21/89 SER 4210, R. Kummerer, STX
CH			Prevent overlaps in raw science segments.
CH
CH		Version 4.4.1 09/04/89 SPR 4480, R. Kummerer, STX
CH			Restrict overlap and gap check to a four week
CH			window about the segment being processed.
CH
CH		Version 4.4.1 10/03/89 SPR 4642, R. Kummerer, STX
CH			Remove FDQ_INIT source code.
CH
CH		Version 4.4.2 10/12/89, SPR 4711, R. Kummerer STX
CH			Correct record counts displayed in report file
CH			when processing by timerange.
CH
CH		Version 4.4.3 11/14/89, SPR 5034, R. Kummerer STX
CH			FPP must note when a channel is missing a raw science
CH			segment.
CH
CH		version 4.4.4 11/28/89, spr 5165, r. kummerer stx
CH			Fails to skip missing channel segments.  Maintain
CH			Proc_Chan_Tbl array for each IN_FILES segment set.
CH
CH		New Version 3/1/91, New requirements, Larry Rosen, STX.
CH		Cleaning up this routine to work with new design.  Major changes
CH		include: using cobetrieve open ct_connect_read means that we
CH		don't need each input file name if opening by timerange.
CH		One output file for each day overlapped unless user chooses
CH		filename, which limits open to one file and write one file
CH		(for each channel).  Since using timerange to read, we don't
CH		need to test all input NFS_ANC extensions against NFS_SDF.
CH		ANC records need to be opened +- an offset (64 seconds) so that
CH		they bracket the science collect time.
C----------------------------------------------------------------------------
C	Input Parameters:
C	  Name		Type		Description
C         ---------------------------------------------------------------------
C	  Time_Present	I*4		Presence of /JSTART,/JSTOP on command line
C	  GMT_start	C*14		Selected start time
C	  GMT_stop	C*14		Selected stop time
C	  File_Present	I*4		Presence of /FILENAME on command line
C	  File_Name	C*39		Selected data segment file
C	  Proc_Chan(4)	I*2		Channels to be processed
C	  Num_Chan	I*2		Number of channels to be processed
C	  Report	L*1		Flag for printing a report file or not
C	  Lu_rpt	I*4		logical unit for report file
C	  Maxn		I*2		Max number of output files per channel
C	  Arch_In	C*(*)		Input archive
C	  Arch_Out	C*(*)		Output archive
C------------------------------------------------------------------------------
C	Output Parameters:
C	  Name			Type	Description
C         ---------------------------------------------------------------------
C	  Num_Days	I*4		Number of output files for each channel
C	  Num_Chan_Tbl(maxn)	I*2	Number of channels to be processed for
C					  each segment
C	  Proc_Chan_Tbl(maxn,4)	I*2	Channel numbers to be processed
C	  In_Files(maxn,5)	C*64	Names of input files
C	  Out_Files(maxn,4)	C*64	Names of files to be created
C	  Num_rec(maxn,4)	I*4	Number of records in each file
C	  Check_JStart(2)	I*4	Base start time from which overlap check
C					  time is computed.
C	  Check_JStop(2)	I*4	Base stop time from which overlap check
C					  time is computed.
C-------------------------------------------------------------------------------
C	Subroutines And Functions Used :
C	------------------------------
C	  CT_INIT
C	  LIB$LOCC
C	  CCT_QUERY_CATALOG
C	  CCT_QUERY_TOD_CATALOG
C	  LIB$SIGNAL
C
C	Include Files :
C	---------------
C	  FPP_MSG			External message params needed for FPP
C	  CT$LIBRARY:CTUSER:INC	        COBETRIEVE return status definitions
C	  $SSDEF				System status values
C	  FUT_PARAMS
C	  CCT_QUERY_CATALOG_RECORD
C	  CCT_QUERY_TOD_CATALOG_RECORD
C-------------------------------------------------------------------------------
C	Processing Method : PDL for FPP_CAT_INFO
C
C	Call CT_Init to initialize to COBETRIEVE.
C	If the NFS_SDF data was selected by filename, Then
C
C	  Determine the channel ID from the filename.
C	  Set the reference channel to this ID.
C	  IF one channel was selected on the command line, Then
C	    Compare the channel ID from the filename, to the selected
C             channel: RH, RL LH or LL.
C	    If the channel ID's match, Then
C	      Set the corresponding channel number in the channel
C               process array.
C	      Set the number of channels to 0ne.
C           Else
C	      Signal an error and set an error status for the return.
C	    Endif
C	  Else if all channels are requested, then
C	    Set all channel numbers in the channel process array.
C	    Set the number of channels to Four.
C	  Endif
C	   Do for each NFS_SDF channel to be processed
C	     If this channel is not the reference channel, then
C	        Build the corresponding filename (including extension
C                 and version number since FDQ requires a matching set)
C	     Endif
C	     Call CCT_Query_Catalog to get the entry in the catalog record.
C            If the entry is valid, then
C	        Call FPP_Save_Cat to save the information from the catalog.
C            Elseif no catalog record is available for this channel, then
C		Signal this information.
C	        If a report is requested, log the missing channel.
C             Else
C		Signal an error and set an error status for the return.
C	        If a report is requested, log the filename and error.
C	     Endif
C	  Enddo for each NFS_SDF channel
C	  Build the corresponding filename (including extension and
C          version number since FDQ requires a matching set) for the NFS_ANC.
C	  Call CCT_Query_Catalog to get the NFS_ANC catalog record.
C          If the entry is valid, then
C	        Call FPP_Save_Cat to save the information from the catalog.
C	  Elseif no catalog record is available, then
C 	      Signal an error and set an error status for the return.
C	      If a report is requested, log the missing filename.
C          Else
C	      Signal an error and set an error status for the return.
C	      If a report is requested, log the filename and error.
C	  Endif
C
C	Elseif the NFS_SDF data was selected by time range, then
C
C	   Call CCT_Query_TOD_Catalog to get the names of all NFS_ANC files
C		in the selected time range.
C           If the return status is good, then
C	      Call FPP_Save_Cat to save the information from the catalog.
C	   Elseif no catalog records are available, then
C	      Signal an error and set an error status for the return.
C	      If a report is requested, log the missing file.
C           Else
C	      Signal an error and set an error status for the return.
C	      If a report is requested, log the filename and error.
C	   Endif
C	   If the return status is still good, then
C            Mark which channels are to be processed in the processing array.
C             Do for each NFS_SDF channel to be processed
C	       Call CCT_Query_TOD_Catalog to get the names of all files
C		 in the selected time range.
C	       If the return status is good, then
C	         Call FPP_Save_Cat to save the information from the catalog.
C	       Elseif no catalog records are available, then
C	         Signal an error and set an error status for the return.
C	         If a report is requested, log the missing channel.
C               Else
C	         Signal an error and set an error status for the return.
C	         If a report is requested, log the dataset and error.
C	       Endif
C             Enddo for each NFS_SDF channel to be processed
C	     If there is no error yet, then
C	       Do for each channel with catalog entries
C	         Store the filenames of existing NFS_SDF segments in the
C                   processing array with the matching NFS_ANC.
C	         If there is no matching science file for an NFS_ANC file, then
C                   Remove the NFS_ANC from the processing array.
C		   If a report is requested, log the filename and message.
C	         Endif
C	         If there is no matching NFS_ANC, then
C	           Signal an error and set an error status for the return.
C	           If a report is requested, log the error.
C                 Endif
C	       Enddo
C	     Endif
C  	     Determine number of daily output files to write.
C	     Construct input and output file names.
C	  Endif return status is still good
C	Endif for method of data selection
C	Return with normal or error status.
C------------------------------------------------------------------------------
	Implicit None

C  Include Files

	Include 'Ct$Library:Ctuser.Inc'
	Include '($Ssdef)'
	Include '(Cct_Query_Catalog_Record)'
	Include '(Cct_Query_Tod_Catalog_Record)'
	Include '(Fut_Params)'
	Include '(Fpp_Msg)'

C  Passed In parameters.

	Integer*4	Time_Present	! presence of /jstart,/jstop on command line
	Character*14	Gmt_Start	! selected start time
	Character*14	Gmt_Stop	! selected stop time
	Integer*4	File_Present	! presence of /filename on command line
	Character*39	File_Name	! user entered file name
	Character*4	Channel		! channel selected: all, rh, rl, lh, ll
	Logical*1	Report		! flag for printing output report
	Integer*4	Lun_Rpt		! logical unit for output report
	Integer*2	Maxn		! maximum # of output files per channel
	Integer*2	Proc_Chan(4)	! channel numbers to be processed
	Integer*2	Num_Chan	! number of channels to be processed
	Character*(*)	Arch_In		! input archive
	Character*(*)	Arch_Out	! output archive

C  Out Parameters.

	Character*64	In_Files(Maxn,5)	! names of input files
	Character*64	Out_Files(Maxn,4)	! names of output files to be created
	Character*30	Time_Range		! time tag for file segment
	Integer*2	Num_Chan_Tbl(Maxn)	! # channel segments for processing
	Integer*2	Proc_Chan_Tbl(Maxn,4)	! channel # to be processed
	Integer*4	Num_Days		! number of output science files, 1 per day
	Integer*4	Num_Rec(Maxn,4)		! number of record in each file
	Integer*4	Cct_Q_No_Cat_Entry	! cat rec not found
	Integer*2	Num_Sci_Seg(4)		! no. rec found
	Integer*4	Check_Jstart(2)		! times from which overlap
	Integer*4	Check_Jstop(2)		! timerange are computed

C  Other Parameters
	Integer*4	Nfiles			! should = maxn
	Parameter	(Nfiles = 50)
	Real*8		Offset			! offset to open anc records
	Parameter	(Offset = 128.)		! to bracket science collection.
	Character*9	Eod /'235959999'/	! end of day
	Character*9	Sod /'000000000'/	! start of day
	Character*9	Eodo /'000104000'/	! end of day + offset
	Character*9	Sodo /'235855999'/	! start of day - offset

C  Dictionaries
	Dictionary	'NFS_SDF'
	Dictionary	'NFS_ANC'
	Dictionary	'Ccm_Cme_Catalog_Entry'
C  Records
	Record/NFS_SDF/Sci_Rec
	Record/NFS_ANC/Anc_Rec
	Record/Query_Catalog/Query_Catalog_Rec
	Record/Query_Tod_Catalog/Query_Tod_Catalog_Rec
	Record/Ccm_Cme_Catalog_Entry/Catalog_Rec(Nfiles,5)

C  Used as Function Values

	Integer*4	Sys$Bintim	! binary time conversion
	Integer*4	Lib$Addx	! extended precision add
	Integer*4	Lib$Subx	! extended precision subtract
	Integer*4	Lib$Locc	! system fn to locate selected string
	Integer*4	Index		! get byte length of the located char
	Integer*4	Cct_Query_Catalog	! cobetrieve catalog routine
	Integer*4	Cct_Query_Tod_Catalog	! cobetrieve catalog routine
	Logical*2	Time_Lt		! ct time comparsion routine
	Logical*2	Time_Gt		! ct time comparsion routine
	Integer*4	Fut_Clean_Catalog	! cobetrieve catalog routine
	Integer*4	FUT_Gmt_Conv	! converts gmt to seconds and back

C  Local Variables

	Integer*4    Rstatus		! function return status
	Integer*4    Char_Pos		! character string location
	Integer*4    Chan		! index for channel number
	Integer*4    Start_Bntime(2)	!start binary time
	Integer*4    Stop_Bntime(2)	!stop binary time
	Integer*4    Ref_Chan		! reference channel number
	Integer*4    Sci_Cur_Seg(4)/4 * 1 / ! current science files pointer
	Integer*4    I,J,Ix,Jx		! index
	Integer*4    Cntr		!file counter
	Integer*4    Num_Seg/0/		! number of files found in archive
	Integer*2    Id			! channel id
	Integer*2    Blank_Len		! positon of blank in a string
	Integer*2    Anc_Counter /0/, Num_Anc_Seg
	Integer*4    File_Cntr /1/	! filenames processing groups counter
	Integer*2    Ct_Stat(20)	! cobetrieve status
	Integer*4    Two_Weeks(2)	! two weeks overlap check window
	Integer*4    Chan_Pres(4)/4*0/  ! check for data from each channel
	Integer*4    Iext
	Character*64 Tempfile			! temporary file name
	Character*20 File_Ext			! temporary file extension
	Character*20 Fext,Fext_Sav		! temporary file extension
	Character*20 Anc_Filext			! final anc file extension
	Character*20 Cert_Ancfile_Ext(Nfiles)	! file extension for anc
	Character*20 File_Fext(Nfiles)		! file extension per file
	Character*20 Files_Ext(Nfiles,5)	! set of file extension
	Character*20 Cert_Scifile_Ext(Nfiles,4)	! file extension for raw
	Character*2  Ref_Chanid 		! file segment channel id
	Logical*1    Zero / 0 /
	Character*1  Blank64(64)/ 64 * ' '/
	Character*1  Blank39(39)/ 39 * ' '/
	Character*1  Blank/' '/
	Character*1  Blank14(14)/ 14 * ' '/
	Character*1  Blank5(5)/ 5 * ' '/
	Logical*1    Match /.False./
	Logical*1    Deleted /.False./
	Logical*1    Del_Flags(Nfiles,5)	! deleted flag for catalog
	Equivalence (Blank,Blank64(1))
	Equivalence (Blank,Blank39(1))
	Equivalence (Blank,Blank14(1))
	Equivalence (Blank,Blank5(1))
	Integer*2    Iyr, Iday, Ihr, Imin, Jyr, Jday	! gmt date/time parts
	Integer*4    Msec			! gmt date/time parts
	Integer*4    Day			! day for segment
	Integer*4    Numdays1			! num_days - 1
	Character*5  Seg			! gmt(1:5) for daily file
	Character*14 Gstarto			! offset start time for anc
	Character*14 Gstopo			! offset stop time for anc
	Character*5  Starto			! offset start time for anc
	Character*5  Stopo			! offset stop time for anc
	Real*8       Tsec			! time in seconds since 1989

C  External used

	External    Cct_Q_No_Cat_Entry
	External    Fut_Normal

C  Initialization

	Fpp_Cat_Info = %Loc(Fpp_Normal)
	If (Nfiles .Ne. Maxn) Write (6,10) Maxn, Nfiles
  10	Format (1X,'*** Warning: Max Number Files In Fpp: ',I,
	1   ' Does Not Match Number Of Files In Fpp_Cat_Info:',I)
	Call Ct_Init( Ct_Stat )
	If (Ct_Stat(1) .Ne. Ctp_Normal) Then
	  Rstatus = Ct_Stat(1)
	  Call Lib$Signal(Fpp_Ctinit,%Val(1),%Val(Rstatus))
	  Fpp_Cat_Info = %Loc(Fpp_Aberr)
	Endif

C  User selects specific filename; not time range.

	If (File_Present .Eq. Fac_Present)  Then
	  If ( File_Name(1:7) .Eq. 'NFS_SDF') Then
	    Char_Pos = Lib$Locc('.' , File_Name)
	    If ( Char_Pos .Ne. 0 ) Then
	      Ref_Chanid = File_Name(Char_Pos-2:Char_Pos-1)
	      If ( Ref_Chanid .Eq. 'RH') Then
	        Ref_Chan = 1
	      Elseif ( Ref_Chanid .Eq. 'RL') Then
	        Ref_Chan = 2
	      Elseif (Ref_Chanid .Eq. 'LH' ) Then
	        Ref_Chan = 3
	      Elseif (Ref_Chanid .Eq. 'LL' ) Then
	        Ref_Chan = 4
	      Else
	        Call Lib$Signal(Fpp_Fileseg)
	        Fpp_Cat_Info = %Loc(Fpp_Aberr)
	      Endif  ! ref_chanid
	    Else
	      Call Lib$Signal(Fpp_Fileseg)
	      Fpp_Cat_Info = %Loc(Fpp_Aberr)
	    Endif ! for char_pos
	  Else
	    Call Lib$Signal(Fpp_Fileseg)
	    Fpp_Cat_Info = %Loc(Fpp_Aberr)
	  Endif				! for file_name is NFS_SDF...

C  User also selects individual channel

	  If ( Num_Chan .Ne. 4 ) Then
	    If ( Ref_Chan .Ne. Proc_Chan(1) ) Then
	      Fpp_Cat_Info = %Loc(Fpp_Aberr)
	      Call Lib$Signal(Fpp_Brefchan,%Val(1),%Val(Ref_Chan))
	    Endif
	  Endif

C  Get file extension

	  Blank_Len = Index(File_Name,Blank)
	  File_Ext = File_Name(Char_Pos:Blank_Len)
	  Ix = 0
	  Num_Days = 1

C  Set query catalog record archive id and filename and query the catalog.

	  Query_Catalog_Rec.Archive_Id = Arch_In
	  Do Chan = 1 , Num_Chan
	    Id = Proc_Chan(Chan)
	    Query_Catalog_Rec.Filename = 'NFS_SDF_' // Fac_Channel_Ids(Id)
	1       // File_Ext
	    Rstatus = Cct_Query_Catalog(Query_Catalog_Rec,Catalog_Rec(1,Id))
	    If (Rstatus .Eq. %Loc(Cct_Q_No_Cat_Entry) ) Then
	      Call Lib$Signal(Fpp_Nocatrec,%Val(1),Query_Catalog_Rec.Filename)
	    Elseif ( Rstatus .NE. 0) Then
	      Num_Sci_Seg(Id) = 1
	      Num_Rec(1,Id) = Catalog_Rec(1,Id).No_Of_Records
	      Iext = Index(Catalog_Rec(1,Id).Filename_Extension, ' ') - 1
	      Fext = Catalog_Rec(1,Id).Filename_Extension(1:Iext)
	      If ( .Not. Catalog_Rec(1,Id).Deletion_Flag ) Then
	        In_Files(1,Id) = Arch_In// 'NFS_SDF_' // Fac_Channel_Ids(Id)
	1           // '.' // Fext
	        Out_Files(1,Id) = Arch_Out// 'FPP_SDF_' // Fac_Channel_Ids(Id)
	1           // '.' // Fext
	        Ix = Ix + 1
	        Num_Chan_Tbl(1) = Num_Chan_Tbl(1) + 1
	        Proc_Chan_Tbl(1,Ix) = Id
	      Endif
	    Else
	      Fpp_Cat_Info = %Loc(Fpp_Aberr)
	      Call Lib$Signal(Fpp_Qcaterr,%Val(1),%Val(Rstatus))
	    Endif				! for catalog query rstatus
	  Enddo					! for each NFS_SDF channel

C  Build the file name for NFS_ANC

	  Query_Catalog_Rec.Filename = 'NFS_ANC' // File_Ext
	  Id = 5
	  Rstatus = Cct_Query_Catalog(Query_Catalog_Rec,Catalog_Rec(1,Id))
	  If (Rstatus .Eq. %Loc(Cct_Q_No_Cat_Entry) ) Then
	    Fpp_Cat_Info = %Loc(Fpp_Aberr)
	    Call Lib$Signal(Fpp_Nocatrec,%Val(1),Query_Catalog_Rec.Filename)
	  Elseif ( Rstatus .NE. 0) Then
	    Iext = Index(Catalog_Rec(1,Id).Filename_Extension, ' ') - 1
	    Fext = Catalog_Rec(1,Id).Filename_Extension(1:Iext)
	    If (Rstatus .NE. 0 .And. ( .Not. Catalog_Rec(1,Id).Deletion_Flag ) )
	1       In_Files(1,Id) = Arch_In // 'NFS_ANC' // '.' // Fext
	    Check_Jstart(1) = Catalog_Rec(1,Id).Initial_Time(1)
	    Check_Jstart(2) = Catalog_Rec(1,Id).Initial_Time(2)
	    Check_Jstop(1) = Catalog_Rec(1,Id).Final_Time(1)
	    Check_Jstop(2) = Catalog_Rec(1,Id).Final_Time(2)
	  Else
	    Fpp_Cat_Info = %Loc(Fpp_Aberr)
	    Call Lib$Signal(Fpp_Qcaterr,%Val(1),%Val(Rstatus))
	    If (Report) Write(Lun_Rpt,20) Query_Catalog_Rec.Filename
  20	    Format (10X,'File: ',A39,' Not Found?')
	  Endif					! for catalog entries rstatus
	  If ( Num_Chan_Tbl(1) .Eq. 0 ) Then
	    Fpp_Cat_Info = %Loc(Fpp_Aberr)
	    Call Lib$Signal(Fpp_Fnumchan)
	  Endif

C  Time range JSTART AND JSTOP are selected by the user.
C  Construct input/output files.

	Elseif ( Time_Present .Eq. Fac_Present ) Then

C  Set query catalog record archive id

	  Query_Tod_Catalog_Rec.Archive_Id = Arch_In
	  Call Ct_Gmt_To_Binary(Gmt_Start,Start_Bntime)
	  Call Ct_Gmt_To_Binary(Gmt_Stop,Stop_Bntime)
	  Do I = 1 , 2
	    Query_Tod_Catalog_Rec.Start_Time(I) = Start_Bntime(I)
	    Query_Tod_Catalog_Rec.Stop_Time(I)  = Stop_Bntime(I)
	  Enddo
	  Query_Tod_Catalog_Rec.Dataset_Name = 'NFS_ANC'
	  Id = 5
	  Rstatus = Cct_Query_Tod_Catalog(Query_Tod_Catalog_Rec,
	1     Catalog_Rec(1,Id),Num_Anc_Seg)
	  If (Rstatus .Eq. 0 .Or. Num_Anc_Seg .Eq. 0) Then
	    Call Lib$Signal(Fpp_Nocatrectod,%Val(1),
	1       Query_Tod_Catalog_Rec.Dataset_Name)
	    Fpp_Cat_Info = %Loc(Fpp_Aberr)
	  Elseif ( Rstatus .NE.0 ) Then
	    Rstatus = Fut_Clean_Catalog(Catalog_Rec(1,Id),Num_Anc_Seg,
	1       Fac_Present)
	    If (Rstatus .Ne. %Loc(Fut_Normal)) Then
	      Call Lib$Signal(Fpp_Clncat,%Val(1),%Val(Rstatus))
	      Fpp_Cat_Info = %Loc(Fpp_Aberr)
	    Else
	      Cntr = 0
	      Call Ct_Gmt_To_Binary(Fac_Jstop_Default,Check_Jstart)
	      Call Ct_Gmt_To_Binary(Fac_Jstart_Default,Check_Jstop)
C  store file extension for out files.
	      Iext = Index(Catalog_Rec(1,Id).Filename_Extension, ' ') - 1
	      Fext_Sav = Catalog_Rec(1,Id).Filename_Extension(1:Iext)
	      Do While (Cntr .Lt. Num_Anc_Seg)
	        Cntr = Cntr + 1
	        Iext = Index(Catalog_Rec(Cntr,Id).Filename_Extension, ' ') - 1
	        File_Ext = Catalog_Rec(Cntr,Id).Filename_Extension(1:Iext)
	        Cert_Ancfile_Ext(Cntr) = File_Ext
	        If(Time_Lt(Catalog_Rec(Cntr,Id).Initial_Time,Check_Jstart))Then
	          Check_Jstart(1) = Catalog_Rec(Cntr,Id).Initial_Time(1)
	          Check_Jstart(2) = Catalog_Rec(Cntr,Id).Initial_Time(2)
	        Endif
	        If (Time_Gt(Catalog_Rec(Cntr,Id).Final_Time,Check_Jstop)) Then
	          Check_Jstop(1) = Catalog_Rec(Cntr,Id).Final_Time(1)
	          Check_Jstop(2) = Catalog_Rec(Cntr,Id).Final_Time(2)
	        Endif
	      Enddo
	    Endif
	  Else
	    Fpp_Cat_Info = %Loc(Fpp_Aberr)
	    Call Lib$Signal(Fpp_Qtodcaterr,%Val(1),%Val(Rstatus))
	  Endif  ! for catalog entries rstatus

C  Query Science files If Return status of Query NFS_ANC still good.

	  Do Chan = 1 , Num_Chan
	    Id = Proc_Chan(Chan)
	    Query_Tod_Catalog_Rec.Dataset_Name = 'NFS_SDF_' //
	1       Fac_Channel_Ids(Id)
	    Rstatus = Cct_Query_Tod_Catalog(Query_Tod_Catalog_Rec,
	1       Catalog_Rec(1,Id), Num_Sci_Seg(Id))
	    If (Rstatus .Eq. 0 .Or. Num_Sci_Seg(Id) .Eq. 0 ) Then
	      Call Lib$Signal(Fpp_Nocatrectod,%Val(1),
	1         Query_Tod_Catalog_Rec.Dataset_Name)
	    Elseif ( Rstatus .NE. 0) Then

C  Extract other pertinent information.

	      Rstatus = Fut_Clean_Catalog(Catalog_Rec(1,Id),Num_Sci_Seg(Id),
	1         Fac_Present)
	      If (Rstatus .Ne. %Loc(Fut_Normal)) Then
	        Call Lib$Signal(Fpp_Clncat,%Val(1),%Val(Rstatus))
	        Fpp_Cat_Info = %Loc(Fpp_Aberr)
	      Else
 	        Cntr = 0
	        Do While (Cntr .Lt. Num_Sci_Seg(Id) .And. Rstatus .NE. 0)
	          Cntr = Cntr + 1
	          Iext = Index(Catalog_Rec(Cntr,Id).Filename_Extension, ' ')-1
	          File_Ext = Catalog_Rec(Cntr,Id).Filename_Extension(1:Iext)
	          Cert_Scifile_Ext(Cntr,Id) = File_Ext
	          Del_Flags(Cntr,Id) = Catalog_Rec(Cntr,Id).Deletion_Flag
	        Enddo
	      Endif
	    Else
	      Fpp_Cat_Info = %Loc(Fpp_Aberr)
	      Call Lib$Signal(Fpp_Qtodcaterr,%Val(1),%Val(Rstatus))
	    Endif ! for catalog query rstatus
	  Enddo   ! for each NFS_SDF channel to be processed

C  If status still good, store the science file names with the matching ANC file

	  If ( Rstatus .NE. 0) Then

C  Loop to get next group of certified file extension.

	    Do While ( (.Not. Match) .And. (Anc_Counter .Lt. Num_Anc_Seg ) )
	      Anc_Counter = Anc_Counter + 1
	      Anc_Filext =  Cert_Ancfile_Ext(Anc_Counter)

C  Compare ancillary file extension to one set of science file extension

	      Ix = 0
	      Do Chan = 1 , Num_Chan
	        Id = Proc_Chan(Chan)
	        Do Jx = Sci_Cur_Seg(Id) , Num_Sci_Seg(Id)
	          If ( Cert_Scifile_Ext(Jx,Id) .Eq. Anc_Filext) Then
	            Sci_Cur_Seg(Id) = Jx      ! update sci current segment
	            If ( Del_Flags(Jx,Id) ) Then
	              Deleted = .True.
	              Tempfile = Arch_In // 'NFS_SDF_' // Fac_Channel_Ids(Id) //
	1                 '.' // Cert_Scifile_Ext(Jx,Id)
	              Call Lib$Signal(Fpp_Segdelt,%Val(1),Tempfile)
	            Else
	              Match = .True.
	              Ix = Ix + 1
	              Proc_Chan_Tbl(File_Cntr,Ix) = Id
	              Num_Chan_Tbl(File_Cntr) = Num_Chan_Tbl(File_Cntr) + 1
	              Chan_Pres(Id) = Chan_Pres(Id) + 1
	              File_Fext(File_Cntr) = '.' // Cert_Scifile_Ext(Jx,Id)
	            Endif    ! del_flag
	          Endif      ! scifile_ext .eq. anc_filext
	        Enddo        ! sci_cur_seg to num_sci_seg
	      Enddo          ! channel 1 to 4

	      If ( Match ) Then
	        File_Cntr = File_Cntr + 1
	        Num_Seg = Num_Seg + 1
	        Match = .False.
	      Elseif ( Deleted ) Then
	        Continue
	      Else
	        Tempfile = Arch_In // 'NFS_ANC' // '.' // Anc_Filext
	        Call Lib$Signal(Fpp_Ancnomat,%Val(1),Tempfile)
	      Endif    ! found
            Enddo      ! do while

C  Note missing raw science segments.

	    Do Ix = 1 , Num_Seg
	      If (Num_Chan_Tbl(Ix) .Lt. Num_Chan) Then
	        File_Ext = File_Fext(Ix)
	        Do Chan = 1 , Num_Chan
	          Id = Proc_Chan(Chan)
	          Match = .False.
	          Do Jx = 1 , Num_Chan_Tbl(Ix)
	            If (Id .Eq. Proc_Chan_Tbl(Ix,Jx)) Then
	              Match = .True.
	            Endif
	          Enddo
		  If (.Not. Match) Then
	            Tempfile = 'NFS_SDF_' // Fac_Channel_Ids(Chan) // File_Ext
	            Call Lib$Signal(Fpp_Nocatrec,%Val(1),Tempfile)
	          Endif
	        Enddo
	      Endif
	    End Do		! for all segments

C  Check for raw science segments that do not match an ancillary file.

	    Do Ix = 1 , Num_Chan
	      Id = Proc_Chan(Ix)
	      Do Jx = Sci_Cur_Seg(Id) , Num_Sci_Seg(Id)
	        Match = .False.
	        Anc_Counter = 0
	        Do While ( (.Not. Match) .And. (Anc_Counter .Lt. Num_Anc_Seg ) )
	          Anc_Counter = Anc_Counter + 1
	          Anc_Filext =  Cert_Ancfile_Ext(Anc_Counter)
	          If ( Cert_Scifile_Ext(Jx,Id) .Eq. Anc_Filext) Then
	            Match = .True.
	          Endif      ! scifile_ext .eq. anc_filext
	        Enddo        ! anc_counter to num_anc_seg
	        If ( .Not. Match ) Then
	          Tempfile = Arch_In // 'NFS_SDF_' // Fac_Channel_Ids(Id) // '.'
	1	      // Cert_Scifile_Ext(Jx,Id)
	          Call Lib$Signal(Fpp_Scinomat,%Val(1),Tempfile)
	        Endif    ! found
	      Enddo      ! sci_cur_seg to num_sci_seg
	    Enddo        ! channel 1 to 4

C  Blank out the filenames buffer

	    Do Ix = 1,5
	      Do Jx = 1 , Maxn
	        In_Files(Jx,Ix) = Blank
	      Enddo
	    Enddo
	    Do Ix = 1,4
	      Do Jx = 1 , Maxn
	        Out_Files(Jx,Ix) = Blank
	      Enddo
	    Enddo

C  Get offset times for ANC files.

	    Read (Gmt_Start,30) Iyr,Iday,Ihr,Imin,Msec
	    Read (Gmt_Stop,30) Jyr,Jday,Ihr,Imin,Msec
  30	    Format (I2,I3,I2,I2,I5)
	    Rstatus = FUT_Gmt_Conv (Tsec,Gmt_Start,2)
	    Rstatus = FUT_Gmt_Conv (Tsec-Offset,Gstarto,1)
	    Rstatus = FUT_Gmt_Conv (Tsec,Gmt_Stop,2)
	    Rstatus = FUT_Gmt_Conv (Tsec+Offset,Gstopo,1)

C  Determine number of daily output files to write.

	    Num_Days = (Jyr-Iyr)*365 + (Jday-Iday) + 1
	    If (Num_Days .Le. 0) Then
	      Fpp_Cat_Info = %Loc(Fpp_Aberr)
	      Return
	    Else

C  Fill in processing channel table.

	       Ix = 0
	       Do Chan = 1 , Num_Chan
	          Id = Proc_Chan(Chan)
	          If (Chan_Pres(Id) .Gt. 0) Then
	             Ix = Ix + 1
	             Do Jx = 1,Num_Days
	                Proc_Chan_Tbl(Jx,Ix) = Id
	             EndDo
	          Endif
	       Enddo
	       Do Jx = 1,Num_Days
	          Num_Chan_Tbl(Jx) = Ix
	       EndDo
	    Endif

C  Create the file names.

	    If (Num_Days .Eq. 1) Then
	      In_Files(1,5) = Arch_In // 'NFS_ANC' // '/' // Gstarto // ';' //
	1         Gstopo // ';'
	      Do Ix = 1 , Num_Chan
	        Id = Proc_Chan(Ix)
	        In_Files(1,Id) = Arch_In // 'NFS_SDF_' // Fac_Channel_Ids(Id)
	1	    // '/' // Gmt_Start // ';' // Gmt_Stop // ';'
	        Out_Files(1,Id) = Arch_Out // 'FPP_SDF_' // Fac_Channel_Ids(Id)
	1	    // '.' // Fext_Sav(1:3) // Gmt_Start(1:9) // ';'
	      Enddo

	    Else
C first day
	      If (Iday.Eq.365) Then
	        Stopo(1:2) = '90'
	        Iday = 0
	      Else
	        Stopo(1:2) = Gmt_Start(1:2)
	      Endif
	      Write (Stopo(3:5),40) Iday + 1
  40	      Format (I3.3)
	      In_Files(1,5) = Arch_In // 'NFS_ANC' // '/' // Gstarto // ';' //
	1         Stopo // Eodo // ';'
	      Do Ix = 1 , Num_Chan
	        Id = Proc_Chan(Ix)
	        In_Files(1,Id) = Arch_In // 'NFS_SDF_' // Fac_Channel_Ids(Id)
	1           // '/' // Gmt_Start // ';' // Gmt_Start(1:5) // Eod // ';'
	        Out_Files(1,Id) = Arch_Out // 'FPP_SDF_' // Fac_Channel_Ids(Id)
	1           // '.' // Fext_Sav(1:3) // Gmt_Start(1:9) // ';'
	      Enddo
C last day
	      If (Jday.Eq.1) Then
	        Starto(1:2) = '89'
	        Jday = 366
	      Else
	        Starto(1:2) = Gmt_Stop(1:2)
	      Endif
	      Write (Starto(3:5),40) Jday - 1
	      In_Files(Num_Days,5) = Arch_In // 'NFS_ANC' // '/' //
	1	  Starto // Sodo // ';' // Gstopo // ';'
	      Do Ix = 1 , Num_Chan
	        Id = Proc_Chan(Ix)
	        In_Files(Num_Days,Id) = Arch_In // 'NFS_SDF_' //
	1           Fac_Channel_Ids(Id) // '/' // Gmt_Stop(1:5) // Sod // ';'
	2           // Gmt_Stop // ';'
	        Out_Files(Num_Days,Id) = Arch_Out // 'FPP_SDF_' //
	1           Fac_Channel_Ids(Id) // '.' // Fext_Sav(1:3) //
	2           Gmt_Stop(1:5) // Sod(1:4) // ';'
	      Enddo
C In Between Days (The Cure)
	      Numdays1 = Num_Days - 1
	      Do Day = 2,Numdays1
	        Jday = Iday + Day - 1
	        If (Jday.Gt.365) Then
	          Jday = Jday - 365
	          Iyr = 90
	        Endif
	        Write(Seg,50) Iyr,Jday
  50	        Format(I2,I3.3)
	        If (Jday.Eq.1) Then
	          Starto = '89365'
	          Stopo = '90002'
	        Elseif (Jday.Eq.365) Then
	          Stopo = '90001'
	          Starto = '89364'
	        Else
	          Write (Starto(3:5),40) Jday-1
	          Write (Starto(1:2),60) Iyr
  60	          Format (I2.2)
	          Write (Stopo(3:5),40) Jday+1
	          Stopo(1:2) = Starto(1:2)
	        Endif
	        In_Files(Day,5) = Arch_In // 'NFS_ANC' // '/' // Starto // Sodo
	1           // ';' // Stopo // Eodo // ';'
	        Do Ix = 1 , Num_Chan
	          Id = Proc_Chan(Ix)
	          In_Files(Day,Id) = Arch_In // 'NFS_SDF_' //
	1             Fac_Channel_Ids(Id) // '/' // Seg // Sod // ';' // Seg //
	2             Eod // ';'
	          Out_Files(Day,Id) = Arch_Out // 'FPP_SDF_' //
	1             Fac_Channel_Ids(Id) // '.' // Fext_Sav(1:3) // Seg //
	2             Sod(1:4) // ';'
	        Enddo
	      Enddo
	    Endif

	  Endif  ! number of days

C  Unrecognized data selection.

        Else
          Call Lib$Signal(Fpp_Fileseg)
          Fpp_Cat_Info = %Loc(Fpp_Aberr)
        Endif  ! for method of data selection

C  Compute overlap and gap check timerange by subtracting two weeks from
C  the segment start time and adding two weeks to the segment stop time.

	Rstatus = Sys$Bintim('14 0:',Two_Weeks)
	Rstatus = Lib$Addx(Check_Jstart,Two_Weeks,Check_Jstart)
	Rstatus = Lib$Subx(Check_Jstop,Two_Weeks,Check_Jstop)

        Return
        End
