        PROGRAM FEC
c
c      By Ken Jensen, STX, April 24, 1991
c      FIRAS SWG Final Requirements version.  
c      
c      Revision of original FEC
c      By Connie Lau, SASC Technologies Inc., May 1986
c      By Reid Wilson, SASC Technologies Inc., July 1986
c
c      Description:
c            FEC analyzes a file of engineering records in order to extract
c            from it all stable plateaus.  A stable plateau is defined as a
c            consecutive series of records starting at the last record with
c            a given commanded change value and proceeding temporally back-
c            wards for at least a given N records. The value of N is located
c            in a reference data set. The plateau ends whenever a different
c            commanded change value is found or a temperature in an engineering
c            record falls outside of the tolerance of the average for that
c            temperature.  The input file is the engineering data archive
c            FDQ_ENG files created by Firas_Data_Qualify (FDQ).
c            FEC products include a report file summarizing the condition of
c            each stable plateau, if requested, and a file for each of the
c            requested science channels which contains short science records  
c            culled from the engineering archive for all of the stable plateau
c            records. The short science files (FEC_SSCAL_ch) serve as pointers,
c            identifying coaddable calibration IFGs for the next facility in
c            the Firas pipeline, Firas_Interferogram_Coadd (FIC).
c
c       Modification History:
c            12/07/87 Reid Wilson
c                     Add check to see whether any good data was found.  If no
c                     good data is found, set return code and signal condition.
c                     Also, change IF-THEN-ENDIF structures with only one state-
c                     ment in them to simple "IF( cond ) CALL X" statements.
c
c	     10/11/88 R. Kummerer
c		      SPR 2531, User specified acceptable XCal positons.
c
c	     12/01/88 J.T.Bonnell, STX@GSFC, version 4.1.1
c		      SER 2379, This program was modified to refer to the
c		      new firas archive logical names in response
c		      to SER 2379 (csdr$firas_in, _out, _raw,
c		      _ref, _uref, and _cal).
c
c	     01/06/89 R. Kummerer
c		      SPR 3098, File closes done only on success status
c		      because FEC_ANALYZE_STABILITY always returned a
c		      success.
c
c            01/20/89 D. Bouler
c                     SPR 2944, Added report_flag to parameter list of 
c                     get_qualifiers,open_files, analyze_stability and
c                     close_files to pass presence/absence of REPORT qualifier.
c
c	     03/08/89 R. Kummerer
c		      SER 3299, Select data by commandable qualities.
c
c	     08/28/89 Q. Chung
c                     SER 3306, Provide version number to track software update.
c
c            01/91-02/91 K. Jensen
c                        FEC Enhancements to satisfy new Requirements.
c                        Subroutine names and argument lists have been
c                        changed to distinguish final requirements version
c                        of FEC from older versions.
c
c            01/17/92  K. Jensen
c                      SPR 9374, Add a flag to enable program to abort
c                                gracefully if there is no data.
c
c            01/23/92     K. Jensen
c                      SPR 9428, Add a  command line qualifier, MAX_DIHED,
c                                to allow for a dihedral temperature quality
c                                check
c            10/20/95     K. Jensen
c                      SPR 12270, Add a  command line qualifier, MIN_DIHED,
c                                 to allow for a dihedral temperature quality
c                                 check
        implicit none
c
c       Include Files:
c
        include '(FEC_INC)'
        include '(FEC_MSG)'
        include '($SSDEF)'
c
c       External Routines:
c
        integer*4       FEC_Load,
     .                  FEC_Analyze,
     .                  FEC_Get_Quals,
     .                  FEC_Open,
     .                  FEC_Close,
     .                  CUT_Register_Version,
     .                  CUT_Display_Banner,
     .                  SYS$GetTim,
     .                  SYS$ASCTim,
     .                  FUT_ERROR
        external        FUT_ERROR
c
        integer*4       report_flag
        integer*4       report_all_flag
        integer*4       rse_flag
        integer*4       query_flag
        integer*4       fakeit_flag
        integer*4       nodata_flag
        integer*4       ct_unit
        integer*4       LU_Report
        integer*4       FUT_Report_LUN
        integer*4       LU_Short_Sci(4)
        integer*4       xcal_error_ctr
        integer*4       temp_change_ctr
        integer*4       temp_change_record
        integer*4       Ret_Code    
        integer*4       Data_Array_Size,
     .                  Plat_Array(0:MAX_PLATEAUS),
     .                  Plat_Array_Size,
     .                  Window_Size
        integer*2       XCal_Pos
        integer*4       Start_Chan
        integer*4       Stop_Chan
        integer*4       Skip_Chan
        character*14    JStart
        character*14    JStop
        character*255   RSE_Filename
        integer*4       Write_flag
        integer*4       Instr_Qual
        integer*4       Attit_Qual
        real*4          Temp_size(Num_Temps_Used)
        real*4          Min_DiHed
        real*4          Max_DiHed
        integer*4       Tol_Type
        integer*4       tstatus
        character*79    command_line(6)
        character*80    report_filename
        integer*4       current_time(2)
        integer*2       time_len
        character*32    runtime
        character*14    rungmt
c
        integer*4       num_vol/80/
        integer*4       lun_out/6/
        character*5     version
        parameter       (version='13.6')
c
        dictionary      'FDQ_ENG'
        record          /FDQ_ENG/ Data_Array(MAX_DATA)
c
c       Initialize error handler:
c
        call lib$establish (FUT_ERROR)
c
c       Obtain current time for run
c

        Call SYS$GetTim( current_time )
        Call CT_Binary_To_GMT( current_time, rungmt )
        Ret_Code = SYS$ASCTim ( time_len, runtime, current_time, 0 )

c
c       Display which version and facility are being run
c

        Ret_Code = CUT_Register_Version(version)
        Ret_Code = CUT_Display_Banner(lun_out,num_vol,
     .			'FIRAS Facility FEC_Extract_Calibration')

        
c       Acquire command qualifiers and keywords from the Command Line.
c
        Ret_Code = FEC_Get_Quals (Start_Chan,
     .                            Stop_Chan,
     .                            Skip_Chan,
     .                            Jstart,
     .                            Jstop,
     .                            Instr_Qual,
     .                            Attit_Qual,
     .                            Query_Flag,
     .                            Fakeit_Flag,
     .                            Report_Flag,
     .                            Report_All_Flag,
     .                            Runtime,
     .                            RunGMT,
     .                            Report_Filename,
     .                            RSE_Flag,
     .                            RSE_Filename,
     .                            Write_Flag,
     .                            XCal_Pos,
     .                            Tol_Type,
     .                            Window_Size,
     .                            Temp_Size,
     .                            Min_DiHed,
     .                            Max_DiHed,
     .                            Command_Line)

c
c       Open the engineering archive ( CSDR$FIRAS_RAW ), open a report
c          file, if report_flag is set.
c
        if( Ret_Code ) 
     -     Ret_Code = FEC_Open (Report_flag,
     .                          CT_unit,
     .                          LU_Report,
     .                          Version,
     .                          Runtime,
     .                          Report_Filename,
     .                          Start_Chan,
     .                          Stop_Chan,
     .                          Skip_Chan,
     .                          Jstart,
     .                          Jstop,
     .                          Write_Flag,
     .                          Query_flag,
     .                          RSE_Flag,
     .                          RSE_Filename,
     .                          Command_Line)

c
c       Read data from the engineering archive into an array of engineering
c          records. Identify calibration plateaus where the commandable
c          temperatures are unchanged and where hot spot power is ON.
c
        if( Ret_Code ) 
     -     Ret_Code = FEC_Load (CT_Unit,
     .                          XCal_Pos,
     .                          Query_Flag,
     .                          RSE_Flag,
     .                          Fakeit_Flag,
     .                          Data_Array,
     .                          Data_Array_Size,
     .                          Plat_Array,
     .                          Plat_Array_Size,)


        nodata_flag = 0
c*      IF there is no data THEN
        if( Ret_Code .and. ( Data_Array_Size .eq. 0 )) then
c*         SET the return code 
           Ret_Code = %loc(FEC_NODATA)
c*         SIGNAL THE ERROR
           call lib$signal( FEC_NODATA )
c*         SET THE NODATA FLAG for FEC_CLOSE
           nodata_flag = 1
c*      END IF
        end if

c
c       Check each calibration plateau for temperature stability. Form
c          short science records from each engineering record belonging
c          to a stable plateau. Write the short science records to an
c          output archive file, if write_flag is set. Write a report
c          file, if report_flag is set.
c

        If ( nodata_flag .eq. 0) Then
         if( Ret_Code ) 
     -     Ret_Code = FEC_Analyze (Report_flag,
     .                             Report_All_flag,
     .                             LU_Report,
     .                             CT_Unit,
     .                             LU_Short_Sci,
     .                             Data_Array,
     .                             Data_Array_Size,
     .                             Plat_Array,
     .                             Plat_Array_Size,
     .                             Window_Size,
     .                             Tol_Type,
     .                             Temp_size,
     .                             Min_DiHed,
     .                             Max_DiHed,
     .                             Start_Chan,
     .                             Stop_Chan,
     .                             Skip_Chan,
     .                             Write_flag,
     .                             Instr_Qual,
     .                             Attit_Qual)
        Endif
c
c
c       Close the files opened by FEC_Open.
c

        tstatus = FEC_Close (report_flag,
     .                       ct_unit, 
     .                       LU_Report, 
     .                       LU_Short_Sci,
     .                       Start_Chan,
     .                       Stop_Chan,
     .                       Skip_Chan,
     .                       nodata_flag,
     .                       Write_flag)

	If (Ret_Code .and. tstatus) Then
	  Call LIB$Signal ( FEC_Normal )
	Else
	  Call LIB$Signal ( %Val(SS$_Abort) )
	End If

	END
