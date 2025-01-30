C-----------------------------------------------------------------------------

	Integer*4 Function FTB_Overlap_Set_Stop ( Cat_Rec, Report, Report_LUN )

C-----------------------------------------------------------------------------
C
C	Purpose: To refine an overlapping segment stop time by using the
C		 timetag of the last raw science record in the overlapping
C		 timerange.
C		 
C
C       Author:  R. Kummerer / STX, October 31, 1989
C
C   	Invocation:  Status = FTB_Overlap_Set_Stop ( Cat_Rec, Report,
C						     Report_LUN )
C
CH	Change Log:
CH         SPR 7020, FTB_OVERLAP CT read error during period:90175-90179
CH               H.Wang, STX, 7/2/90
CH
C-----------------------------------------------------------------------------
C
C       Input Files :
C
C	Output Files :
C
C  	Input Parameters:
C     
C	  Name           Type		Description
C         ---------------------------------------------------------------------
C         Cat_Rec	 Record		Raw science catalog record with an
C					approximate non-overlapping start time
C	  Report	 L*1		Generate a report file.
C	  Report_LUN	 I*4		Report file LUN.
C-----------------------------------------------------------------------------
C
C	Output Parameters:
C	  Name            Type		Description
C         ---------------------------------------------------------------------
C         Cat_Rec	 Record		Raw science catalog record with a
C					precise non-overlapping start time
C-----------------------------------------------------------------------------
C
C	 Subroutines And Functions Used :
C	 ------------------------------
C	 LIB$GET_LUN
C	 LIB$FREE_LUN
C        LIB$SIGNAL
C
C	 Common Variables Used :
C
C	 Name		Type	 Use		Description
C-----------------------------------------------------------------------------
C	 Include Files :
C	 ---------------
C-----------------------------------------------------------------------------
C
C	 Processing Method : PDL for FTB_Overlap_Set_Stop
C
C			PDL for FTB_Overlap_Set_Stop
C			--------------------------------
C
C	 CALL LIB$GET_LUN to get a logical unit number on which to open
C	      the raw science data file
C	 OPEN the raw sceince data file using CT_CONNECT_READ and the data
C	      timerange
C	 CALL CT_READ_ARCV_REV to read the last raw science record
C	 UPDATE the catalog record initial time with the timetag of the
C	      raw science data record
C	 CALL CT_CLOSE_ARCV to close the raw science data file
C	 CALL LIB$FREE_LUN to free the logical unit number on which
C	      the raw science data file was open
C
C-----------------------------------------------------------------------------
C

	Implicit None

C	Include Files  

	Include	'($SSDef)'
	Include 'CT$Library:CTUser.Inc'

C       Dictionaries

	Dictionary 'NFS_SDF'
	Dictionary 'CCM_CME_Catalog_Entry'

C       Records

	Record /NFS_SDF/ SCI_REC

C       Passed parameters.

	Record /CCM_CME_CATALOG_ENTRY/ Cat_Rec
	Logical*1	Report
	Integer*4	Report_LUN

C	Out Parameters.

C       Local Variables

	Character*64	Catalog_File
	Character*14	GMT_Start
	Character*14	GMT_Stop
	Integer*4	LUN
	Integer*4	Status
	Integer*2	CT_Status(20)
	Integer*4	j
	Integer*4	Retstat
	Integer*4	Success / 1 /, Error / 2 /

C  	Used as Function Values

	Integer*4	LIB$Get_LUN
	Integer*4	LIB$Free_LUN
	Integer*4	CT_Connect_Read

C       Externals used

	External	CT_Connect_Read
	External	FTB_Normal
	External	FTB_Aberr
	External	FTB_OpenSDF
	External	FTB_ReadSDF
	External	FTB_CloseSDF

C Fetch a LUN and open the raw science segment catalog archive file.

	Retstat = Success

	Status = LIB$Get_LUN ( LUN )

	If (Status .Eq. SS$_Normal) Then

	   Call CT_Binary_To_GMT ( Cat_Rec.Initial_Time, GMT_Start )
	   Call CT_Binary_To_GMT ( Cat_Rec.Final_Time, GMT_Stop )

	   Catalog_File = 'CSDR$FIRAS_RAW:' // Cat_Rec.Dataset_Name // '/' //
	1		  GMT_Start // ';' // GMT_Stop // ';'

	   Open ( Unit=LUN, File=Catalog_File, Status='Old', 
	1		  IOStat=Status, Useropen=CT_Connect_Read )

	   If (Status .Eq. 0) Then

C Read the entry from the archive.

	      Call CT_Read_Arcv_Rev ( , LUN, SCI_REC, CT_Status )

	      If (CT_Status(1) .Ne. CTP_Normal) Then
		 Retstat = Error
	 	 Call LIB$Signal ( FTB_ReadSDF, %Val(2),
	1			Cat_Rec.Dataset_Name, %Val(CT_Status(1)) )
		 If (Report) Then
		    Write (Report_LUN, 10) CT_Status(1)
		 Endif
	      Else
		 Cat_Rec.Final_Time(1) = SCI_REC.CT_Head.Time(1)
		 Cat_Rec.Final_Time(2) = SCI_REC.CT_Head.Time(2)
	      Endif

C Close the file.

	      Call CT_Close_Arcv ( , LUN, CT_Status )

	      If (CT_Status(1) .Ne. CTP_Normal) Then
                 Retstat = Error
	         Call LIB$Signal ( FTB_CloseSDF, %Val(2),
	1			Cat_Rec.Dataset_Name, %Val(CT_Status(1)) )
		 If (Report) Then
		    Write (Report_LUN, 20) CT_Status(1)
		 Endif
	      Endif

	   Else

	      Retstat = Error
	      Call LIB$Signal ( FTB_OpenSDF, %Val(2),
	1			Cat_Rec.Dataset_Name, %Val(Status) )
	      If (Report) Then
		 Write (Report_LUN, 40) Status
	      Endif

	   Endif

	   Status = LIB$Free_LUN ( LUN )

	   If (Status .Ne. SS$_Normal) Then
	      Call LIB$Signal ( %Val(Status) )
	      If (Report) Then
	         Write (Report_LUN, 30) Status
	      Endif
	   Endif

	Else

	   Retstat = Error
	   Call LIB$Signal ( %Val(Status) )
	   If (Report) Then
	      Write (Report_LUN, 50) Status
	   Endif

	Endif

	If (Retstat .Eq. Success) Then
	   FTB_Overlap_Set_Stop = %Loc(FTB_Normal)
	Else
	   FTB_Overlap_Set_Stop = %Loc(FTB_Aberr)
	Endif

C	Formats statements.

10	Format (//, 5x, '*** FTB_Overlap_Set_Stop CT READ error = ', i)
20	Format (//, 5x, '*** FTB_Overlap_Set_Stop CT CLOSE error = ', i)
30	Format (//, 5x, '*** FTB_Overlap_Set_Stop LIB$FREE_LUN error = ', i)
40	Format (//, 5x, '*** FTB_Overlap_Set_Stop CT OPEN error = ', i)
50	Format (//, 5x, '*** FTB_Overlap_Set_Stop LIB$GET_LUN error = ', i)

	Return
	End
