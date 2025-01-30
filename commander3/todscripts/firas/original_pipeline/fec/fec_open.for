        integer*4 function FEC_OPEN( Report_flag,
     .                               CT_Unit,
     .				     LU_Report,
     .                               Version,
     .                               Runtime,
     .                               Report_Filename,
     .                               Start_Chan,
     .                               Stop_Chan,
     .                               Skip_Chan,
     .                               JStart,
     .				     JStop,
     .				     Write_Flag,
     .				     Query_Flag,
     .                               RSE_Flag,
     .				     RSE_Filename,
     .                               Command_Line)
c
c
c       By Ken Jensen, STX, 25-FEB-1991
c
c       Revision of FEC_OPEN_FILES_A.FOR
c	Author:  C. Lau
c		 SASC Technologies
c		 May 1986
c
c       Revision made to distinguish between final requirements FEC
c       and older versions of FEC.
c
c	version 4.1.1 12/01/88, ser 2379, J.T.Bonnell, STX@GSFC
c		This program was modified to refer to the
c		new firas archive logical names in response
c		to SER 2379 (csdr$firas_in, _out, _raw,
c		_ref, _uref, and _cal).
c       version 4.2.? 01/20/89 spr 2944, D. Bouler, stx gsfc
c               added report_flag to parameter list so that
c               presence/absence of REPORT qualifier is passed.
c       version 4.4 05/17/89 spr 3805, R. Kummerer stx gsfc
c		Fetch RSE from archive.
c       version 4.4 06/08/89 spr 3972, R. Kummerer stx gsfc
c       	Replace LIB$GET_LUN / LIB$FREE_LUN with analogous FUT
c		routines.
c	version 4.4.1 08/28/89 ser 3306, Q. Chung stx
c               Provide version number to track software update.
c
c       Jan-Feb, 1991 : K. Jensen (STX)
c                       FEC enhancements to satisfy new SWG requirements.
c
c	Include files
c		 FEC_INC
c                FEC_MSG
c                FUT_ERROR
c                FUT_PARAMS
c	
c	Subroutines called
c		 
       implicit none
c
       include  'CT$LIBRARY:ctuser.inc'      ! Access COBETRIEVE routines
       include '(FEC_INC)'
       include '(FEC_MSG)'
c       external FEC_MSG
       include '($SSDEF)'
       include '($JPIDEF)'
       include '(FUT_PARAMS)'                ! record size definitions
       include '(FUT_ERROR)'
       dictionary 'FDQ_ENG'                  ! structure definition of eng rec
c
       integer * 4     report_flag       !flag to indicate REPORT presence    
       integer * 2     ct_status(20)     !I/O status for COBETRIEVE read
       integer * 4     ct_unit           !unit number selected by 
                                         !COBETRIEVE
       integer * 4     Ret_Code
       integer * 4     LU_Report
       integer * 4     io_stat           !I/O status for GRT_SWITCH read
       integer * 4     status

       character * 32  runtime
       character * 80  report_filename
       character * 30  time_range
       character * 128 rse(16)/16*' '/
       character * 128 rse0
       logical   * 1   oldrse/.false./
       character * 14  jstart
       character * 14  jstop
       integer * 4     cjstart
       integer * 4     cjstop
       character * 255 RSE_Filename
       integer * 4     Query_flag
       integer * 4     RSE_flag
       integer * 4     Write_flag
       integer * 4     Start_Chan
       integer * 4     Stop_Chan
       integer * 4     Skip_Chan
       integer * 4     LU_RSE
       integer * 4     RSE_Index/1/
       integer * 4     nrse
       integer * 4     FUT_Get_LUN
       integer * 4     FUT_Free_LUN
       character * 64  in_file
       integer * 4     FUT_Query_RSE
       integer * 4     ct_connect_query
       external        ct_connect_query
       integer * 4     ct_connect_read
       external        ct_connect_read
       integer * 4     cut_translate_archive_id
       external        cut_translate_archive_id
       integer * 4     FUT_Error
       external        FUT_Error
       integer * 4     cli$get_value
       integer * 4     comm_len
       character * 512 comm_str
       character * 79  command_line(6)
       record /FDQ_ENG/ dum_eng
       integer * 4     counter

       integer * 4     CUT_Register_Version
       integer * 4     CUT_Display_Banner
       character * 8   owner
c       integer * 4     LIB$GetJPI
       character * 72  line
       character * 1   lin(72) / 72 * ' ' /
       equivalence     ( line, lin(1) )
       character * 72  log
       character * 15  log15
       character * 72  flog
       integer * 4     flen
       character * 72  tlog
       integer * 4     tlen
       integer * 4     num_vol/80/
       integer * 4     rstatus
       character * 5   version
c
c      ******* BEGIN *******
c
       io_stat = 0
       CALL ct_init(ct_status)


c Open report file. First check to see if report flag is set.
c
       if (report_flag .eq. 1) then
           If( (io_stat .eq. 0) .and. (ct_status(1) .eq. ctp_normal)) then
                status = fut_get_lun(LU_Report)
                If (.Not. status) Then
	          FEC_Open = status
	          Call LIB$Signal(%Val(status))
                  Return
	        End If
                OPEN (unit=LU_Report, file=report_filename,
     -                status='new',iostat=io_stat)
                IF (io_stat .ne. 0) THEN
                   call lib$signal(fec_openrep,%val(1),%val(io_stat))
                   FEC_Open=%loc(FEC_OPENREP)
                   return
                ELSE

                    rstatus = CUT_Register_Version(version)
                    rstatus = CUT_Display_Banner(lu_report,num_vol,
     -                        'FIRAS Facility FEC_Extract_Calibration')

                    fut_report_lun=lu_report
                    call lib$establish(fut_error)

                ENDIF

                Call LIB$GetJPI (JPI$_UserName,,,,owner,)
                Write (LU_Report,55)
                Write (LU_Report,60) owner
	        Write (LU_Report,65) runtime


c     Write Command Line Invocation to Report File

                Write (LU_Report,70)
                Write (LU_Report,75) command_line(1)
                Write (LU_Report,75) command_line(2)
                Write (LU_Report,75) command_line(3)
                Write (LU_Report,75) command_line(4)
                Write (LU_Report,75) command_line(5)
                Write (LU_Report,75) command_line(6)

c     Write Logical Translations to Report File

                Write (LU_Report,80)
                log(1:72) = line
                log15 = 'CSDR$FIRAS_RAW '
                log(1:15)=log15
                call str$upcase(log,log)
                status = cut_translate_archive_id(log, flog, flen, tlog, tlen)
                Write (LU_Report,85)log15,tlog
                log15 = 'CSDR$FIRAS_REF '
                log(1:15)=log15
                call str$upcase(log,log)
                status = cut_translate_archive_id(log, flog, flen, tlog, tlen)
                Write (LU_Report,85)log15,tlog
                log15 = 'CSDR$FIRAS_OUT '
                log(1:15)=log15
                call str$upcase(log,log)
                status = cut_translate_archive_id(log, flog, flen, tlog, tlen)
                Write (LU_Report,85)log15,tlog

           Endif
       Endif

       If ( Query_flag .Eq. 1) Then
         If( ct_status(1) .eq. ctp_normal) 
     -      Call ct_query_bld (ctu_$firas, time_range, rse, ct_status, oldrse)
       Else
         time_range = JStart // ';' // JStop // ';'
         ct_Status(1) = ctp_normal
         if( rse_flag .eq. 1 )then
	   If (rse_filename .Eq. ' ') Then
               Status = %loc(FEC_NORSEFILE)
               Call LIB$Signal(FEC_NORSEFILE)
	       FEC_Open = status
	       Return
	   Else
             status = FUT_Get_LUN(LU_RSE)
	     If (.Not. status) Then
	       FEC_Open = status
	       Call LIB$Signal(%Val(status))
	       Return
	     End If
             Open (LU_RSE, status='OLD',file=RSE_Filename,iostat=io_stat,readonly)
             Do While ((io_stat .eq. 0) .and. (RSE_Index .lt. FEC_MAX_RSE))
               Read (LU_RSE, '(A128)', iostat=io_stat) RSE(RSE_Index)
               RSE_Index = RSE_Index + 1
             End Do
             Close (LU_RSE)
	     If (io_stat .Ne. 0) Then
	       Call errsns (io_stat,,,,io_stat)
	       FEC_Open = io_stat
	       Call LIB$Signal(%Val(io_stat))
	       Return
	     End If
             status = FUT_Free_LUN(LU_RSE)
	     If (.Not. status) Then
               FEC_Open = status
	       Call LIB$Signal(%Val(status))
               Return
	     End If
           End If
         end If
       End If

c
c  Expansion of all RSEs to report file.
c
       If ( query_flag .eq. 1 ) Then
          Write (LU_Report,90)
          Write (LU_Report,95)
          Write (LU_Report,97) Time_Range
       Endif
       If ( rse_flag .eq. 1 ) Then
          Write (LU_Report,92) RSE_Filename(1:70)
          Write (LU_Report,95)
          nrse=1
          Do While (nrse .le. 16)
             rse0=rse(nrse)
             Write (LU_Report,100) rse0(1:79)
             nrse = nrse + 1
          End do
       Endif

       Write (LU_Report,105)


c Open engineering archive.
c
       If( ct_status(1) .Ne. ctp_normal) Then
          FEC_Open = %loc(FEC_CTQBLDRERROR)
          Call lib$signal(FEC_CTQBLDRERROR)
          Return
       Else

	 status = fut_get_lun(ct_unit)
	 if(.not. status)then
	    fec_open = status
	    call lib$signal(%val(status))
	    return
	 end if

	 in_file = 'CSDR$FIRAS_RAW:fdq_eng/' // time_range

         IF ( query_flag .eq. 1 .Or. rse_flag .eq. 1 ) THEN

  	    open (unit=ct_unit, file=in_file, status='old',
     .		iostat=io_stat, useropen=ct_connect_query)
            If( io_stat .eq. 0) Then
                Call CT_Query_Arcv ( , ct_unit, rse, ct_status)
                If ( ct_Status(1) .ne. ctp_normal) Then
                   FEC_Open = %loc(FEC_QARCVERR)
                   Call lib$signal (FEC_QARCVERR)
                   Return
                End If
            Else
                FEC_Open = %loc(FEC_OPENARCVERR)
                Call lib$signal (FEC_OPENARCVERR)
                Return
            End If

         ELSE
  	    open (unit=ct_unit, file=in_file, status='old',
     .		iostat=io_stat, useropen=ct_connect_read)
            if ( io_stat .ne. 0 ) then
                FEC_Open = %loc(FEC_OPENARCVERR)
                Call lib$signal (FEC_OPENARCVERR)
                Return
            endif

         ENDIF
c
c
c Check condition status.
c
         If( (io_stat .eq. 0) .and. (ct_status(1) .eq. ctp_normal)) then
             FEC_Open = %loc(FEC_NORMAL)
         Else
             If (ct_status(1) .NE. ctp_normal) Then
                If (ct_status(1) .Eq. ctp_incomplete) Then
                   FEC_Open = %loc(FEC_USERUNDFND)
                   call lib$signal (FEC_USERUNDFND)
                   Return
                Else If (ct_status(1) .EQ. ctp_lunquota_exceeded) Then
                   FEC_Open = %loc(FEC_NOOPENLUN)
                   call lib$signal (FEC_NOOPENLUN)
                   Return
                Else
                   FEC_Open = %loc(FEC_BADLUN)
                   call lib$signal (FEC_BADLUN, %val(1),%val(ct_status(1)))
                   Return
                End If
             Else If( io_stat .ne. 0) Then
                call errsns (io_stat,,,,io_stat)
                FEC_Open = io_stat
                call lib$signal (%val(io_stat))
                Return
             End If
         End If
       End If


55    Format (//,' ',20x,'FEC_EXTRACT_CALIBRATION Calibration Plateau Report',//)
60    Format (' Run by :     ',a,/)
65    Format (' Run Time :   ', a,/)
70    Format (//,' Command Line Invocation : ',/)
75    Format (' ',a)
80    Format (///,' Logical Translations :',/)
85    Format ('  ',a,' = ',a)
90    Format (////,' Cobetrieve Query Builder Used :',/)
92    Format (////,' RSE File Used : ',/,'     ',a,/)
95    Format (' Record Selection Expressions Applied to Archive :',/)
97    Format ('  TIME_RANGE EQ ',a)
100   Format ('  ',a)
105   Format (//,' ',//)

       FEC_Open = %loc(FEC_NORMAL)
       Return
       End
