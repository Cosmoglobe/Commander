	Integer*4  Function  FAD_Read_Reference_Data  ( arch_cal, modelext,
	1                                               chan, scan, lun_rpt,
	2                                               report, fex_ejv_rec,
	3                                               file_ejv, access_time )

c------------------------------------------------------------------------------
c   Purpose: Read the FEX_EJV reference file.
c
c   Input Parameters:
c     character*15    arch_cal    -- reference archive
c     character*22    modelext    -- FISH model name - FEX_EJV file extension
c     character*2     chan, scan  -- FIRAS channel and scan mode (2 letter)
c     integer*4       lun_rpt     -- report file logical unit number
c     logical*1       report      -- flag whether to write report
c     integer*4       access_time -- time for reference data
c
c   Output Parameters:
c      record /fex_ejv/ fex_ejv_rec  -- ejv reference record
c      character*50    file_ejv     -- ejv file name
c
c   Include Files:
c      FAD_msg       -- External message params needed for FAD
c
c   Functions:
c      Lib$Signal
c
c   Author:  Larry Paris Rosen, Hughes STX, 20 April 1993
c   Modified: L. Rosen, February 1994.  Read FEX_GN file in addition to FEX_EJ
c      file.  Now uses cct_get_config.  Add 20 seconds to config open time to
c      get a start time for the reference records.
c   Modified: L. Rosen, April 1994.  Change get config time to be skymap time,
c      passed in as "access_time".  It is found in FAD_open_archives.for.
c   Modified: L. Rosen, May 1994.  Eliminate all gain references.  No more
c      fex_gn file or record.
c   Modified: L. Rosen, June 1994.  EJ is now EJV and has a vibration
c      correction, and more than 1 time constant.  Corrected spectrum is now
c      CSP = SP - Sum_ufos Ai exp (-t/taui) - Sum_tophat Bj Topj (t)
c               - P4 (t) Vib,
c      where P4 is the quartic vibration correction.
c------------------------------------------------------------------------------
	Implicit None

c Include

	Include		'(fad_msg)'
	Include		'(cct_get_config)'	! Needed for config_status

c Passed Parameters

	Character*15	arch_cal			! reference archive
	Character*22	modelext      ! FISH model name - FEX_EJV file extension
	Character*2	chan, scan     ! FIRAS channel and scan mode (2 letter)
	Integer*4	lun_rpt				! unit number of report
	Logical*1	report
	Dictionary	'FEX_EJV'
	Record /FEX_EJV/ fex_ejv_rec
	Character*51	file_ejv
	Integer*4	access_time (2)

c Function

	Integer*4	CCT_Close_Config
	Integer*4	CCT_Open_Config
	Integer*4	CCT_Get_Config_Tod

c Local

	Integer*4	rstat

C  Config stuff

	Character*14	config_GMT_start / '86001000001000' /
	Character*14	config_GMT_stop / '99365235958990' /
	Integer*4	config_start (2)
	Integer*4	config_stop (2)
	Integer*4	config_size_ejv (1) / 27648 /
	Character*51	config_name_ejv (1)
	Integer*4	config_ref_count
	Logical*1	config_new_segment (1)
	Integer*4	config_nrecs (1) / 1 /
	Integer*4	config_index (1)
	Record /Config_Status/ config_status(1)
	Integer*4	config_lun (1)

c------------------------------------------------------------------------------
c Begin

	FAD_Read_Reference_Data = %Loc (FAD_Normal)

c Use cct_config methods of opening reference file.

	Call CT_GMT_To_Binary (config_gmt_start, config_start)
	Call CT_GMT_To_Binary (config_gmt_stop, config_stop)
	file_ejv = arch_cal // 'FEX_EJV_' // chan // scan // '.' // modelext
	config_name_ejv (1) = arch_cal // 'FEX_EJV_' // chan // scan

C Get ej record

c Open the sequential access reference data sets.

	rstat = CCT_Open_Config ( config_start, config_stop, 1,
	1                         config_name_ejv, config_size_ejv,' ', 1,
	2                         config_lun, config_index, config_status,
	3                         config_ref_count )

	If (.NOT. rstat) Then
	   FAD_Read_Reference_Data = %Loc (FAD_Abort)
	   Call Lib$Signal (fad_openerr, %Val(2), file_ejv, %Val(rstat))
	Else
	   if (report) Then
	      Write (lun_rpt, 10) file_ejv
  10	      Format (1X, 'Successfully opened: ',A)
	   Endif

	   rstat = CCT_Get_Config_Tod ( access_time, 1, config_size_ejv,
	1                               config_lun, config_index, fex_ejv_rec,
	2                               config_new_segment, config_status )

	   If (.NOT. rstat) Then
	      FAD_Read_Reference_Data = %Loc (FAD_Abort)
	      Call Lib$Signal (fad_readerr, %Val(2), file_ejv, %Val(rstat))
	   Else

c Close reference file.

	      rstat = CCT_Close_Config ( 1, config_lun, config_index )
	      If (.NOT. rstat) Then
	         FAD_Read_Reference_Data = %Loc (FAD_Abort)
	         Call Lib$Signal (fad_closerr, %Val(2), file_ejv, %Val(rstat))
	      Elseif (report) Then
	         Write (lun_rpt, 20) file_ejv
  20	         Format (1X, 'Successfully closed: ',A)
	      Endif
	   Endif
	Endif

	Return
	End
