	Integer*4 Function FSD_Astroplots_Catinfo(File_Seg,scitype,
	2                     GMT_Start,GMT_Stop,numchans,sel_chan,extension)
C-----------------------------------------------------------------------------
C
C	This routine is called only if the user has selected a segment 
C	(filename).  It receives the name of the file requested and returns
C	start and stop times in GMT (Char*14) format of that segment.
C
C	AUTHOR:  Fred Shuman, STX, 1989 Jan 25
C	         adapted from FXT_Cat_Info by Shirley Read
C
C	INVOCATION:
C	  status = FSD_Astroplots_Catinfo(File_Seg,GMT_Start,GMT_Stop,
C	                                  numchans,sel_chan,extension)
C
C	INPUT LIST:
C	    Ch* 39     File_Seg      ! Selected NFS_SDF filename
C	    Ch*  3     scitype       ! Type of science data: NFS, FPP, or FDQ
C	    Ch* 14     GMT_Start     ! Start time of data
C	    Ch* 14     GMT_Stop      ! Stop time of data
C	     I*  4     numchans      ! Number of channels selected by user
C	     I*  4     sel_chan(4)   ! Channels selected by user
C
C	OUTPUT LIST:
C	    Ch* 14     GMT_Start     ! Start time of data
C	    Ch* 14     GMT_Stop      ! Stop time of data
C	     I*  4     numchans      ! Number of channels selected by user
C	     I*  4     sel_chan(4)   ! Channels selected by user
C	    Ch* 20     extension     ! Filename extension of skymap
C	                               to be written later
C	SUBROUTINES CALLED:
C
C	COMMON VARIABLES USED:
C	    None
C
C	INCLUDE FILES:
C	    CT$LIBRARY:CTUSER.INC
C	    $SSDEF
C-----------------------------------------------------------------------------
C  Changes:
C
C	v4.4 Access violation following call to CCT_Query_TOD_Catalog.
C	    SPR 4007.  Fred Shuman, STX.  1989 Jun 15.
C
C	SPR 3945, 4178, Allow Astroplots to work on the new FPP_SDF and FDQ_SDF
C	    files, as well as the raw science, NFS_SDF, so that it can run
C	    whether or not FPP or FDQ have.  FSD_Astroplots, _Init, _Catinfo,
C	    and FSD.CLD.  Fred Shuman, STX / 1989 Aug 29.
C
C	SPR 6430, Access violation from CCT_QUERY_CATALOG_TOD; increase limit
C	    on maximum number of raw science segments that can be queried.
C	    R. Kummerer, STX / 1990 Mar 8.
C-----------------------------------------------------------------------------

	Implicit None

	Include   'ct$library:ctuser.inc'
	Include   '($SSDEF)'
	Include   '(CUT_Params)'

	Integer   *  4      status

	Character * 39      File_Seg      ! Selected NFS_SDF filename
	Character *  3      scitype       ! Type of sci data: NFS, FPP, or FDQ
	Character * 14      GMT_Start     ! Start time of segment
	Character * 14      GMT_Stop      ! Stop time of segment
	Integer   *  4      numchans      ! Number of channels selected by user
	Integer   *  4      sel_chan(4)   ! Channels selected by user
	Character * 20      extension     ! Skymap filename extension

	Integer   *  4      i
	Integer   *  4      j
	Integer   *  2      Ix              ! Index
	Integer   *  4      Data_Start(2) ! Binary start time of data
	Integer   *  4      Data_Stop(2)  ! Binary stop time of data
	Integer   *  4      loc
	Integer   *  2      numrecs
	Character *  2      channame
	Character *  2      chanlist(4) /'RH', 'RL', 'LH', 'LL'/
	Character * 39      Blank		! Blanks
	Character *  1      Blanks(39) / 39 * ' ' /
	Equivalence  ( Blank, Blanks(1) )

	Dictionary 'CCM_CME_Catalog_Entry'	!Info retrieved from CT
	Record     /CCM_CME_Catalog_Entry/CatRec
	Record     /CCM_CME_Catalog_Entry/CatRecs(CUT$_Max_Segments)

	Include '(CCT_Query_Catalog_Record)'
	Record  /Query_Catalog/QCatRec

	Include '(CCT_Query_TOD_Catalog_Record)'
	Record  /Query_TOD_Catalog/QTODCatRec

	Integer   *  4      CCT_Query_Catalog
	Integer   *  4      CCT_Query_TOD_Catalog
	Integer   *  4      Lib$Locc

	External            FSD_Normal
	External            FSD_Aberr1
	External            FSD_NoCatRec
	External            FSD_QCatErr
	External            CCT_Q_No_Cat_Entry

	FSD_Astroplots_CatInfo = %loc(FSD_Normal)
c
c  A specific NFS_SDF segment has been specified by the user.
c  Call the COBETRIEVE Query Catalog routine for the information.
c
c   Get catalog information for the NFS_SDF file.
c

	If (file_seg .Eq. blank) Then
c
c  User has supplied a timerange, so we have GMT_Start & Stop, but no file_seg.
c  We need the first two chars of the input filename extension (the data type).
c
	   QTODCatRec.Archive_Id = 'CSDR$FIRAS_ARCHIVE'
	   QTODCatRec.Dataset_Name = scitype // '_SDF_' // chanlist(sel_chan(1))
	   Call CT_GMT_to_Binary(GMT_Start, Data_Start)
	   Call CT_GMT_to_Binary(GMT_Stop, Data_Stop)
	   Do i=1,2
	      QTODCatRec.Start_Time(i) = Data_Start(i)
	      QTODCatRec.Stop_Time(i)  = Data_Stop(i)
	   End Do

	   Status = CCT_Query_TOD_Catalog( QTODCatRec, CatRecs, numrecs )
	   
	   extension( 1: 2) = CatRecs(1).Filename_Extension(1:2)
	   extension( 3:18) = '_' // GMT_Start(1:7) // '_' // GMT_Stop(1:7)

	Else
c
c  User has supplied a file_seg, so we need to find GMT_Start & Stop.
c
	   QCatRec.Archive_Id = 'CSDR$FIRAS_ARCHIVE'
	   QCatRec.Filename(1:39) = File_Seg(1:39)
	   Status = CCT_Query_Catalog( QCatRec, CatRec)
	   If (.Not. Status ) Then
C
C   Failed to get catalog record
C
	      FSD_Astroplots_CatInfo = %loc(FSD_Aberr1)
	      If (Status .Eq. %loc(CCT_Q_No_Cat_Entry)) Then
	         Call Lib$Signal( FSD_NoCatRec, %val(1), File_Seg)
	      Else
	         Call Lib$Signal(FSD_QCatErr,%val(2),%val(Status))
	      Endif
	   Else If (Status .Eq. %loc(CCT_Q_No_Cat_Entry)) Then
	      Call Lib$Signal( FSD_NoCatRec, %val(1), File_Seg)
	      FSD_Astroplots_CatInfo = %loc(FSD_Aberr1)
	   Else
C
C Call to CCT_query_catalog is OK.
C   Save the time range for the data in the segment.
C
	      Do Ix = 1, 2
	         Data_Start(Ix) = CatRec.Initial_Time(Ix)
	      Enddo
	      Call CT_Binary_To_GMT ( Data_Start, GMT_Start )
	      GMT_Stop = CatRec.Final_Time_Key(1:14)
	      loc = Lib$Locc('.', file_seg)
	      channame = file_seg(loc-2:loc-1)
	      j = 1
	      Do While (channame .Ne. chanlist(j) .And. j .Le. 5)
	         j = j + 1
	      End Do
	      numchans = 1
	      sel_chan(1) = j
	   Endif

	Endif

	End
