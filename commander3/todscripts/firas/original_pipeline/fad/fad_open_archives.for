	Integer*4  Function  FAD_Open_Archives  ( lun_in, lun_out, lun_rpt,
	1                                         report, chan, scan, filext,
	2                                         arch_in, arch_out, filein,
	3                                         fileout, access_time )

c------------------------------------------------------------------------------
c   Purpose: Open the input, reference, and output files.
c
c   Input Parameters:
c      integer*4     lun_rpt       -- logical unit number for report file
c      logical*1     report        -- true if report is to be written (default)
c      character*2   chan          -- channel (RH, RL, LH, or LL)
c      character*2   scan          -- scan mode (SS, SF, LS, or LF)
c      character*22  filext        -- file extension of FCF_SKY file
c      character*14  arch_in       -- input archive name
c      character*15  arch_out      -- output archive name
c
c   Output Parameters:
c      integer*4     lun_in        -- logical unit number for input file
c      integer*4     lun_out       -- logical unit number for output file
c      character*50  filein, fileout  -- file names
c      integer*4     access_time   -- time for reference data
c
c   Include Files:
c      FAD_msg       -- External message params needed for FAD
c      $SSdef        -- System status values
c      FUT_params    -- FIRAS global parameters
c      csdr$library:ctparams.inc  -- Cobetrieve parameters
c      cct_query_catalog_record
c
c   Functions and Subroutines:
c      Lib$Get_Lun
c      Lib$Signal
c      CT_Init
c      CSA_Field_Offset_Values
c      CCT_Query_Catalog
c
c   Author:  Larry Paris Rosen, Hughes STX, 19 April 1993
c   Modifications:  Larry Rosen, HSTX, June 7, 1993, SPR 11022.  Need to add
c                   readonly to open statement for FEX_EJ reference file.
c                   L. Rosen, February 1994.  Move reference data opening to
c                   FAD_READ_REFERENCE_DATA file.
c                   L. Rosen, April 1994.  Use catalog query to get time of
c                   skymap to pass to reference file getting routine.
c------------------------------------------------------------------------------
	Implicit None

c Include files

	Include		'($ssdef)'
	Include		'(fut_params)'
	Include		'(fad_msg)'
	Include		'csdr$library:ctparams.inc'
	Include		'(cct_query_catalog_record)'

c Passed parameters

	Integer*4	lun_rpt				! unit number of report
	Logical*1	report
	Character*2	chan, scan     ! FIRAS channel and scan mode (2 letter)
	Character*22	filext                 ! File extension of FCF_SKY file
	Character*14	arch_in				! input data archive
	Character*15	arch_out			! output data archive
	Integer*4	lun_in				! unit number for input
	Integer*4	lun_out				! unit number for out
	Character*50	filein, fileout			! file names
	Integer*4	access_time (2)		! time for reference file

c Functions

	Integer*4	Lib$Get_Lun
	Integer*4	CSA_Field_Offset_Values
	Integer*4	CCT_Query_Catalog

c External

	External	csa_normal
	External	csa_open_skymap
	External	csa_field_offset_values

c Local

	Integer*4	rstat		! Return status
	Integer*2	ct_stat(20)	! Cobetrieve status
	Integer*2	skymaplen	! skymap length in longwords
	Integer*4	f_pix_offset	! I*4 of fac_coad_spec_pix_offset
	Integer*4	f_time_offset	! I*4 of fac_time_offset
	Dictionary 'ccm_cme_catalog_entry'
	Record /query_catalog/		query_cat
	Record /ccm_cme_catalog_entry/	cats(5)

c------------------------------------------------------------------------------

c Begin

	FAD_Open_Archives = %Loc (FAD_Normal)

c Get logical unit numbers.

	rstat = Lib$Get_Lun (lun_in)
	If (rstat .NE. SS$_Normal) Then
	   FAD_Open_Archives = %Loc (FAD_Abort)
	   Call Lib$Signal (fad_lunerr, %Val(1), %Val(rstat))
	Else
	   rstat = Lib$Get_Lun (lun_out)
	   If (rstat .NE. SS$_Normal) Then
	      FAD_Open_archives = %Loc (FAD_Abort)
	      Call Lib$Signal (fad_lunerr, %Val(1), %Val(rstat))
	   Endif
	Endif

c Call CT_INIT to initialize Cobetrieve.

	If (FAD_Open_Archives .EQ. %Loc (FAD_Normal)) Then
	   Call CT_INIT ( ct_stat )
	   If (ct_stat(1) .NE. ctp_normal) Then
	      Call LIB$SIGNAL (fad_ctinit, %val(1), %val(ct_stat(1)))
	      FAD_Open_Archives = %Loc (FAD_Abort)
	   Endif
	Endif

c Open FCF input file.

	If (FAD_Open_Archives .EQ. %Loc (FAD_Normal)) Then
	   filein = arch_in // 'FCF_SKY_' // chan // scan // '.' // filext
	   skymaplen = fac_coad_spec_size / 4
	   f_pix_offset = fac_coad_spec_pix_offset
	   f_time_offset = fac_time_offset
	   Open ( Unit=lun_in, File=filein, Status='OLD', Form='UNFORMATTED',
	1         Recordtype='FIXED', READONLY, Recl=skymaplen,
	2         Useropen=CSA_Open_Skymap, Iostat=rstat )

	   If (rstat .NE. 0) Then
	      FAD_Open_Archives = %Loc (FAD_Abort)
	      Call Lib$Signal (fad_csaopen, %Val(2), filein, %Val(rstat))
	   Else
	      If (report) Then
	         Write (lun_rpt,*)
	         Write (lun_rpt, 10) filein
  10	         Format (1X, 'Successfully opened: ',A)
	      Endif
	      rstat = CSA_Field_Offset_Values (f_pix_offset, f_time_offset, -1,
	1                                      lun_in)
	      If (rstat .NE. %Loc (CSA_Normal)) Then
	         FAD_Open_Archives = %Loc (FAD_Abort)
	         Call Lib$Signal (fad_csaoffset, %Val(1), %Val(rstat))
	      Endif
	   Endif
	Endif

c Query catalog for time of input skymap for later use in getting reference
c files.

	If (FAD_Open_Archives .EQ. %Loc (FAD_Normal)) Then
	   query_cat.archive_id = arch_in
	   query_cat.filename = 'FCF_SKY_' // chan // scan // '.' // filext
	   rstat = cct_query_catalog (query_cat, cats(1))
	   If (.NOT. rstat) Then
	      FAD_Open_Archives = %Loc (FAD_Abort)
	      Call Lib$Signal ( fad_query_cat, %Val (1), %Val (rstat) )
	   Else
	      access_time (1) = cats(1).initial_time (1)
	      access_time (2) = cats(1).initial_time (2)
	   Endif
	EndIf

c Open FAD output file.

	If (FAD_Open_Archives .EQ. %Loc (FAD_Normal)) Then
	   fileout = arch_out // 'FAD_SKY_' // chan // scan // '.' // filext
	   Open ( Unit=lun_out, File=fileout, Status='NEW', Form='UNFORMATTED',
	1         Recordtype='FIXED', Recl=skymaplen,
	2         Useropen=CSA_Open_Skymap, Iostat=rstat )

	   If (rstat .NE. 0) Then
	      FAD_Open_Archives = %Loc (FAD_Abort)
	      Call Lib$Signal (fad_csaopen, %Val(2), fileout, %Val(rstat))
	   Else
	      If (report) Then
	         Write (lun_rpt, 10) fileout
	      Endif
	      rstat = CSA_Field_Offset_Values (f_pix_offset, f_time_offset, -1,
	1                                      lun_out)
	      If (rstat .NE. %Loc (CSA_Normal)) Then
	         FAD_Open_Archives = %Loc (FAD_Abort)
	         Call Lib$Signal (fad_csaoffset, %Val(1), %Val(rstat))
	      Endif
	   Endif
	Endif
	Return
	End
