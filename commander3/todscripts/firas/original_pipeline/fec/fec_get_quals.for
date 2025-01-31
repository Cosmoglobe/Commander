        integer*4 function FEC_GET_QUALS(Start_Chan,
     .                                   Stop_Chan,
     .                                   Skip_Chan,
     .                                   JStart,
     .                                   JStop,
     .                                   Instr_Qual,
     .                                   Attit_Qual,
     .                                   Query_Flag,
     .                                   Fakeit_Flag,
     .                                   Report_Flag,
     .                                   Report_All_Flag,
     .                                   Runtime,
     .                                   RunGMT,
     .                                   Report_filename,
     .                                   RSE_flag,
     .                                   RSE_filename,
     .                                   Write_flag,
     .                                   Xcal_pos,
     .                                   Tol_Type,
     .                                   Window_Size,
     .                                   Temp_Size,
     .                                   Min_DiHed,
     .                                   Max_DiHed,
     .                                   Command_Line)
c
c      By Ken Jensen, STX, 25-FEB-1991
c
c      Revision of FEC_GET_QUALIFIERS
c      Written by Connie Lau, SASC Technologies Inc., May 1986
c      Commented by Reid Wilson, SASC Technologies Inc.,  22-AUG-1986
c
c      Revision to satisfy SWG final requirements.
c
c      Modification History:
c            
c      20-DEC-1987 Reid Wilson, STX
c            Modified to write short science records into split files.
c
c      20-jan-1989 D. Bouler, stx
c            added report flag to calling parameter list in order
c            to report whether the REPORT qualifier was present.
c
c      23-Jan-1992 K. Jensen, STX
c            SPR 9428, added max_dihed qualifier.
c
c      20-Oct-1995 K. Jensen, Hughes STX
c            SPR 12270, added min_dihed qualifier.
c
c
c      Description:
c            This function gets the information supplied by the user through
c            the command line.
c
c      Format:
c            Ret_Code = FEC_GET_QUALS(Start_Chan,
c                                     Stop_Chan,
c                                     Skip_Chan,
c                                     JStart,
c                                     JStop,
c                                     Instr_Qual,
c                                     Attit_Qual,
c                                     Query_Flag,
c                                     Fakeit_Flag,
c                                     Report_Flag,
c                                     Report_All_Flag,
c                                     Runtime,
c                                     RunGMT,
c                                     Report_filename,
c                                     RSE_flag,
c                                     RSE_filename,
c                                     Write_flag,
c                                     Xcal_Pos,
c                                     Tol_Type,
c                                     Window_Size,
c                                     Temp_size,
c                                     Min_DiHed,
c                                     Max_DiHed,
c                                     Command_Line)
c
c      Input Parameters:
c            None
c
c      Output Parameters:
c
c            Start_Chan           = first channel to process
c                                   integer*4
c            Stop_Chan            = last channel to process
c                                   integer*4
c            Skip_Chan            = channel skip increment
c                                   integer*4
c            JStart               = Start time of input file data set
c                                   character*14
c            JStop                = Stop time of input file data set
c                                   character*14
c            Instr_Qual           = Instrument data quality threshold
c                                   integer*4
c            Attit_Qual           = Attitude data quality threshold
c                                   integer*4
c            Query_Flag           = Indicates whether the Query Builder is to
c                                   be used to build the RSE list : integer*4
c            Fakeit_Flag          = Indicates whether to select Fakeit or
c                                   No Fakeit records : integer*4
c            Report_Flag          = Indicates whether report is to be produced
c                                   integer*4
c            Report_All_Flag      = Indicates whether report file will include
c                                   null and hot spot plateaus : integer*4
c            Runtime              = Real Time of Facility Invocation
c                                   character*32
c            RunGMT               = Runtime in GMT units
c                                   character*14
c            Report_Filename      = Report filename
c                                   character*80
c            RSE_Flag             = Indicates whether RSE list will be used
c                                   integer*4
c            RSE_Filename         = Name of file containing RSE list
c                                   character*255
c            Write_Flag           = Indicates whether short science records will
c                                   be written to the output archive : integer*4
c            Xcal_Pos             = acceptable Xcal Positions
c                                   integer*2
c            Tol_Type             = 1=Absolute Tolerance : 2=Relative Tolerance
c                                   integer*4
c            Window_Size          = Number of points in window look ahead
c                                   integer*4
c            Temp_Size            = Amount temperature may differ from average
c                                   real*4 Temp_Size(Num_Temps_Used)
c            Min_DiHed            = Minimum allowable Dihedral Temperature
c                                   real*4
c            Max_DiHed            = Maximum allowable Dihedral Temperature
c                                   real*4
c            Command_Line         = Command line
c                                   character*79 (4)
c
c      Modules Called:
c            UPM_Present          = tests for presence of qualifier
c            UPM_Get_Value        = gets string value of qualifier
c            UPM_Get_Float        = gets a real*4 number
c            UPM_Get_LongWord     = gets an integer*4 number
c                  
c      Include Files:            
c            $SSDEF             
c            UPM$USER:UPM_STAT_MSG.PAR = UPM routines and parameters
c            FEC_INC                   = FEC parameters
c            FEC_MSG                   = FEC message definition
c
c--------------------------------------------------------------------------
c Changes
c	SPR 2531, R. Kummerer, User specified acceptable XCal positons.
c	October 11, 1988.
c	SER 3299, R. Kummerer, Select data by commandable quality.
c	March 8, 1989.
c	SER 6851, H. Wang, Select Tolerance method by commandable quality.
c	June 13, 1990.
c       SER 7567, K. Jensen (STX), Add Xcal_position keyword "Error".
c       Nov. 2, 1990.
c       SER 7816, K. Jensen (STX), FEC Enhancements
c
c--------------------------------------------------------------------------

        implicit none
c
        INCLUDE '($SSDEF)'
        INCLUDE '(UPM_STAT_MSG)'      ! UPM routines and parameters
        include '(FEC_INC)'           ! FEC parameters
        include '(FEC_MSG)'           ! FEC message definition
        include '(FUT_Params)'
c
        integer*4       report_flag      !flag to indicate REPORT presence 
        integer*4       report_all_flag  !flag to indicate REPORT of null
                                         !and hot spot plateaus
        integer*4       ret_status       !return status
        integer*4       status           !status from upm subroutines
        integer*4       status2          !status from upm subroutines
        integer*4       upm_present      !The extreme bozonic variable
                                         !indicates the level of linguini
                                         !in the algorithm.
        integer*4       upm_get_float    !routine to get real*4 number
        integer*4       upm_get_longword !routine to get integer*4 number
        integer*4       upm_get_value    !routine to get a string
        integer*4       num_qual/16/     !number of qualifiers
        integer*4       num_key(16)      !number of keywords
        integer*4       i
        integer*4       j
        integer*4       nj
        Integer*4       Tol_Type
        integer*4       upm_status
        integer*4       current_time(2)
        integer*2       time_len
        character*32    runtime
        character*30    text
        character*1     qtext
        character*2     qtext2
        character*4     trsize
        character*11    temp_text(Num_Temps_Used)
        character*11    ttext
        character*11    text_default
        real*4          temp_default
        real*4          dihed_default(2)
        integer*4       temp_len(Num_Temps_Used)
        real*4          Temp_Size(Num_Temps_Used)
        real*4          Tsize(Num_Temps_Used)
        real*4          Min_DiHed
        real*4          Max_DiHed
        integer*4       Window_Size
        integer*4       Query_Flag
        integer*4       Fakeit_Flag
        integer*4       RSE_Flag
        integer*2       XCal_Pos
        Integer*4       start_chan
        Integer*4       stop_chan
        Integer*4       skip_chan
        character*14    JStart
        character*14    JStart_Default/'85001000000000'/
        character*14    JStop
        character*14    JStop_Default/'99365235959990'/
        character*14    RunGMT
        integer*4       Cut_Translate_Archive_ID
        external        Cut_Translate_Archive_ID
        character*72    log
        character*15    log15
        character*72    flog
        character*72    tlog
        integer*4       flen
        integer*4       ttlen
        integer*4       tlen
        character*255   RSE_Filename
        integer*4       Write_Flag
        integer*4       Instr_Qual
        integer*4       Attit_Qual
        integer*4       quality
        integer*4       rsize
        real*4          frsize
        integer*4       clen
        integer*4       clen0
        character*256   Invoke_Line
        character*79    Command_line(6)
        character*79    Command_line1
        character*79    Command_line2
        character*79    Command_line3
        character*79    Command_line4
        character*79    Command_line5
        character*79    Command_line6
        Integer*2       len
        Character*32    keyword(16,9)

        character*16    qualifier(16)    !command line qualifiers
        character*80    def_files(2)
	character*80    report_filename
c
c	   Invocation command line
c
        status = UPM_Get_Value ( '$LINE', Invoke_Line, len )

        qualifier(1) = 'channel'
        qualifier(2) = 'jstart'
        qualifier(3) = 'jstop'
        qualifier(4) = 'quality'
        qualifier(5) = 'query'
        qualifier(6) = 'fakeit'
        qualifier(7) = 'write'
        qualifier(8) = 'report'
        qualifier(9) = 'reportall'
        qualifier(10)= 'rse'
        qualifier(11)= 'xcal_position'
        qualifier(12)= 'tolerance'
        qualifier(13)= 'window_size'
        qualifier(14)= 'temp_size'
        qualifier(15)= 'min_dihed'
        qualifier(16)= 'max_dihed'

c   Assign keywords

        do i=1,num_qual
           num_key(i)=1
        enddo
        num_key(1)=9
        num_key(4)=2
        num_key(11)=4
        num_key(12)=2


        keyword(1,1) = 'channel.rh'
        keyword(1,2) = 'channel.rl'
        keyword(1,3) = 'channel.lh'
        keyword(1,4) = 'channel.ll'
        keyword(1,5) = 'channel.right'
        keyword(1,6) = 'channel.left'
        keyword(1,7) = 'channel.high'
        keyword(1,8) = 'channel.low'
        keyword(1,9) = 'channel.all'

        keyword(4,1) = 'quality.instrument'
        keyword(4,2) = 'quality.attitude'

        keyword(11,1) = 'xcal_position.in'
        keyword(11,2) = 'xcal_position.out'
        keyword(11,3) = 'xcal_position.transit'
        keyword(11,4) = 'xcal_position.error'

        keyword(12,1) = 'tolerance.absolute'
        keyword(12,2) = 'tolerance.relative'

        command_line1(1:13)=' FEC/CHANNEL='
        clen=13

c       set up default qualifiers
c
        start_chan = 1
        stop_chan = 4
        skip_chan = 1
        jstart = jstart_default
        jstop = jstop_default
        instr_qual = fac_many_yellow
        attit_qual = 33
        query_flag = 0
        fakeit_flag = 0
        report_flag = 1
        report_all_flag = 0
        rse_flag = 0
        write_flag = 0
        xcal_pos = fac_xcalin
        tol_type = 2        !
        window_size = 7     !
        temp_default=0.005   !
        text_default(1:11)='     0.0050'
        dihed_default(1) = 0.0
        dihed_default(2) = 3.5
        
        do i = 1, num_temps_used
           temp_size(i)=temp_default
           temp_text(i)=text_default
           temp_len(i)=6
           tsize(i)=-9999.
        end do
c
c          DO FOR ALL QUALIFIERS
c
        ret_status = %loc(FEC_NORMAL)
	DO i = 1, num_qual

	   status = upm_present(qualifier(i))
           IF (i .eq. 1) THEN

              IF (status .EQ. upm_absent .or. status .EQ. upm_defaulted) THEN
          
                 command_line1(clen+1:clen+3)='ALL'
                 clen=clen+3

              ENDIF

              IF (status .EQ. upm_pres) THEN

c             << Get Start_Chan, Stop_Chan, Skip_Chan >>

                 nj = 0
                 do j=1,num_key(i)
                    status = upm_present(keyword(i,j))
                    if(status .eq. upm_pres)then
                       if(j .eq. 1)then
                          start_chan = 1
                          stop_chan = 1
                          skip_chan = 1
                          command_line1(clen+1:clen+2)='RH'
                          clen=clen+2
                          nj=1
                       else if(j .eq. 2)then
                          start_chan = 2
                          stop_chan = 2
                          skip_chan = 1
                          command_line1(clen+1:clen+2)='RL'
                          clen=clen+2
                          nj=1
                       else if(j .eq. 3)then
                          start_chan = 3
                          stop_chan = 3
                          skip_chan = 1
                          command_line1(clen+1:clen+2)='LH'
                          clen=clen+2
                          nj=1
                       else if(j .eq. 4)then
                          start_chan = 4
                          stop_chan = 4
                          skip_chan = 1
                          command_line1(clen+1:clen+2)='LL'
                          clen=clen+2
                          nj=1
                       else if(j .eq. 5)then
                          start_chan = 1
                          stop_chan = 2
                          skip_chan = 1
                          command_line1(clen+1:clen+5)='RIGHT'
                          clen=clen+5
                          nj=1
                       else if(j .eq. 6)then
                          start_chan = 3
                          stop_chan = 4
                          skip_chan = 1
                          command_line1(clen+1:clen+4)='LEFT'
                          clen=clen+4
                          nj=1
                       else if(j .eq. 7)then
                          start_chan = 1
                          stop_chan = 3
                          skip_chan = 2
                          command_line1(clen+1:clen+4)='HIGH'
                          clen=clen+4
                          nj=1
                       else if(j .eq. 8)then
                          start_chan = 2
                          stop_chan = 4
                          skip_chan = 2
                          command_line1(clen+1:clen+3)='LOW'
                          clen=clen+3
                          nj=1
                       else if(j .eq. 9)then
                          start_chan = 1
                          stop_chan = 4
                          skip_chan = 1
                          command_line1(clen+1:clen+3)='ALL'
                          clen=clen+3
                          nj=1
                       endif
                    endif
                 end do

                 if ( nj .eq. 0 ) then
                       start_chan = 1
                       stop_chan = 4
                       skip_chan = 1
                       command_line1(clen+1:clen+3)='ALL'
                       clen=clen+3
                 endif
              ENDIF
           
           command_line1(clen+1:clen+8)='/JSTART='
           clen=clen+8

	   ELSE if(i .eq. 2)then

c             << Get Jstart >>

              IF (status .EQ. upm_absent .or. status .EQ. upm_defaulted) THEN
          
                 command_line1(clen+1:clen+14)=jstart_default(1:14)
                 clen=clen+14

              ENDIF

              IF (status .EQ. upm_pres) THEN
                 status = upm_get_value(qualifier(i),jstart,len)
                 if(status .eq. upm_absent)then
                    jstart = jstart_default
                 else
                    if (len .lt. 14) jstart(len+1:14)=jstart_default(len+1:14)
                 endif

                 command_line1(clen+1:clen+14)=jstart(1:14)
                 clen=clen+14

              ENDIF

           command_line1(clen+1:clen+7)='/JSTOP='
           clen=clen+7

	   ELSE if(i .eq. 3)then

c             << Get Jstop >>

              IF (status .EQ. upm_absent .or. status .EQ. upm_defaulted) THEN
          
                 command_line1(clen+1:clen+14)=jstop_default(1:14)
                 clen=clen+14

              ENDIF

              IF (status .EQ. upm_pres) THEN
                 status = upm_get_value(qualifier(i),jstop,len)
                 if(status .eq. upm_absent)then
                    jstop = jstop_default
                 else
                    if (len .lt. 14) jstop(len+1:14)=jstop_default(len+1:14)
                 endif

                 command_line1(clen+1:clen+14)=jstop(1:14)
                 clen=clen+14

              ENDIF

              command_line1(clen+1:clen+2)=' -'

           ELSE if(i .eq. 4)then

c             << Get Instr_Qual, Attit_Qual >>


              command_line2(1:26)='    /QUALITY=(INSTRUMENT=3'
              command_line2(27:38)=',ATTITUDE=33'
              command_line2(39:39)=')'
              clen = 39

              IF (status .EQ. upm_pres) THEN

                 do j=1,num_key(i)
                    status = upm_get_longword(keyword(i,j),quality)
                    if(status .eq. ss$_normal)then
                       if(j .eq. 1)then
                          if(quality .gt. 7)then
                             quality=7
                             qtext='7'
                          else
                             write(qtext,10)quality
 10                          Format(I1)
                          endif
                          instr_qual = quality
                          command_line2(26:26)=qtext(1:1)
                       else if(j .eq. 2)then
                          if(quality .le. 9)then
                             write(qtext,10)quality
                             command_line2(37:37)=qtext(1:1)
                             command_line2(38:38)=')'
                             clen = 38
                          else if(quality .gt. 9)then
                             write(qtext2,101)quality
 101                         Format(I2)
                             command_line2(37:38)=qtext2(1:2)
                             command_line2(39:39)=')'
                             clen = 39
                          endif
                          attit_qual = quality
                       endif
                    endif
                 end do

              ENDIF


           ELSE if(i .eq. 5)then

              IF (status .EQ. upm_absent .OR. status .EQ. upm_negated) THEN
          
                 command_line2(clen+1:clen+8)='/NOQUERY'
                 clen=clen+8

              ENDIF

              IF (status .EQ. upm_pres) THEN

c             <<  Query enabled  >>

                 query_flag = 1
                 command_line2(clen+1:clen+6)='/QUERY'
                 clen=clen+6

              ENDIF


           ELSE if(i .eq. 6)then

              IF (status .EQ. upm_absent .OR. status .EQ. upm_negated) THEN
          
                 command_line2(clen+1:clen+9)='/NOFAKEIT'
                 clen=clen+9

              ENDIF

              IF (status .EQ. upm_pres) THEN

c             <<  Select Fakeit Records  >>

                 fakeit_flag = 1
                 command_line2(clen+1:clen+7)='/FAKEIT'
                 clen=clen+7

              ENDIF

           ELSE if(i .eq. 7)then

              IF (status .EQ. upm_absent .OR. status .EQ. upm_negated) THEN

c             << Disable write of short science records to output archive >>

                 write_flag = 0
                 command_line2(clen+1:clen+10)='/NOWRITE -'
                 clen=clen+10

              ENDIF

              IF (status .EQ. upm_pres) THEN

c             << Enable write of short science records to output archive >>

                 write_flag = 1
                 command_line2(clen+1:clen+8)='/WRITE -'
                 clen=clen+8

              ENDIF


           ELSE if(i .eq. 8)then

              IF (status .EQ. upm_pres .OR. status .EQ. upm_defaulted) THEN
          
c             << Report enabled. Get Report_Filename >>

                 report_flag = 1
                 command_line3(1:12)='    /REPORT='
                 clen=12

                 IF (status .EQ. upm_pres) THEN
                   status2 = upm_get_value(qualifier(i),report_filename,tlen)
                   if (status2 .EQ. upm_absent) then
		      report_filename = 'FEC_' // jstart(1:7) // '_' 
     .                                  // jstop(1:7) // '.REP_' // rungmt(1:9)
                      command_line3(clen+1:clen+33)=report_filename(1:33)
                      clen=clen+33
                   else
                     command_line3(clen+1:clen+tlen)=report_filename(1:tlen)
                     clen=clen+tlen
	           end if
                 ELSE
                   report_filename = 'FEC_' // jstart(1:7) // '_' 
     .                                  // jstop(1:7) // '.REP_' // rungmt(1:9)
                   command_line3(clen+1:clen+33)=report_filename(1:33)
                   clen=clen+33
                 ENDIF

              ENDIF

              IF (status .EQ. upm_negated) THEN

c             <<  Report disabled  >>

                 report_flag = 0
                 command_line3(1:15)='    /NOREPORT -'
                 clen=13

              ENDIF


           ELSE if(i .eq. 9)then

              IF (status .EQ. upm_absent .OR. status .EQ. upm_negated) THEN
          
                 command_line3(clen+1:clen+14)='/NOREPORTALL -'
                 clen=clen+14

              ENDIF

              IF (status .EQ. upm_pres) THEN

c             <<  Report null and hot spot plateaus  >>

                 report_all_flag = 1
                 command_line3(clen+1:clen+12)='/REPORTALL -'
                 clen=clen+12

              ENDIF


           ELSE if(i .eq. 10)then

              IF (status .EQ. upm_negated .OR. status .EQ. upm_absent) THEN

c             <<  /NORSE  >>

                 rse_flag = 0
                 command_line4(1:10)='    /NORSE'
                 clen=10

              ENDIF

              IF (status .EQ. upm_pres ) THEN
          
c             << RSE enabled. Get RSE_Filename >>

                 rse_flag = 1
                 command_line6(1:9)='    /RSE='

                 status2 = upm_get_value(qualifier(i),rse_filename,tlen)
                 if (status2 .eq. upm_absent) then
                    Ret_Status = %loc(FEC_NORSEFILE)
                    call lib$signal (FEC_NORSEFILE)
                    FEC_Get_Quals = Ret_Status
                    return
                 else
                    command_line6(10:9+tlen)=rse_filename(1:tlen)
                 endif

              ENDIF


           ELSE if(i .eq. 11)then

              if(rse_flag .eq. 0)then
                 command_line4(clen+1:clen+15)='/XCAL_POSITION='
                 clen=clen+15
              endif
              if(rse_flag .eq. 1)then
                 command_line4(1:19)='    /XCAL_POSITION='
                 clen=19
              endif

c             << Get Xcal_pos >>

              IF (status .EQ. upm_absent .or. status .EQ. upm_defaulted) THEN

                 xcal_pos=fac_xcalin
                 command_line4(clen+1:clen+2)='IN'
                 clen=clen+2

              ENDIF

              IF (status .EQ. upm_pres) THEN

                 do j=1,num_key(i)
                    status = upm_present(keyword(i,j))
                    if(status .eq. upm_pres)then
                       if(j .eq.1)then
                          xcal_pos = fac_xcalin
                          command_line4(clen+1:clen+2)='IN'
                          clen=clen+2
                       else if(j .eq. 2)then
                          xcal_pos = fac_xcalout
                          command_line4(clen+1:clen+3)='OUT'
                          clen=clen+3
                       else if(j .eq. 3)then
                          xcal_pos = fac_xcaltrans
                          command_line4(clen+1:clen+4)='TRANSIT'
                          clen=clen+4
                       else if(j .eq. 4)then
                          xcal_pos=fac_xcalposerr
                          command_line4(clen+1:clen+5)='ERROR'
                          clen=clen+5
                       endif
                    endif
                 end do

                 if (clen .eq. clen0) then
                    xcal_pos = fac_xcalin
                    command_line4(clen+1:clen+2)='IN'
                    clen=clen+2
                 endif

              ENDIF


           ELSE if(i .eq. 12)then

              command_line4(clen+1:clen+11)='/TOLERANCE='
              clen=clen+11

c             << Get Tol_Type >>

              IF (status .EQ. upm_absent .or. status .EQ. upm_defaulted) THEN

                 tol_type = 2
                 command_line4(clen+1:clen+8)='RELATIVE'
                 clen=clen+8

              ENDIF

              IF (status .EQ. upm_pres) THEN

                 status = upm_get_value(qualifier(i),text,tlen)
                 if (text(1:3) .eq. 'ABS') then
                    tol_type = 1
                    command_line4(clen+1:clen+8)='ABSOLUTE'
                    clen=clen+8
                 else if(text(1:3) .eq. 'REL') then
                    tol_type = 2
                    command_line4(clen+1:clen+8)='RELATIVE'
                    clen=clen+8
                 else
                    tol_type = 2
                    command_line4(clen+1:clen+8)='RELATIVE'
                    clen=clen+8
                 end if

              ENDIF


           ELSE if(i .eq. 13)then

              command_line4(clen+1:clen+13)='/WINDOW_SIZE='
              clen=clen+13

c             << Get Window_Size >>

              IF (status .EQ. upm_absent .or. status .EQ. upm_defaulted) THEN

                 window_size = 7
                 command_line4(clen+1:clen+1)='7'
                 clen=clen+1
   
              ENDIF

              IF (status .EQ. upm_pres) THEN

                 status = upm_get_longword(qualifier(i), rsize, len)
                 window_size = rsize
                 write(trsize,40)rsize
 40              Format(I4)
                 tlen=1
                 if(rsize .ge. 10)tlen=2
                 if(rsize .ge. 100)tlen=3
                 if(rsize .ge. 1000)tlen=4                    
                 command_line4(clen+1:clen+tlen)=trsize(5-tlen:4)
                 clen=clen+tlen
                 command_line4(clen+1:clen+2)=' -'

              ENDIF

           ELSE if(i .eq. 14)then
c
c             << Get Temp_Sizes >>
c
                 
              IF (status .EQ. upm_pres) THEN

                 j = 0
                 do while ( (j .lt. Num_Temps_Used) .and. Ret_Status)
                    j = j + 1
                    UPM_Status = UPM_Get_Float ( qualifier(i), Tsize(j) )
                    if ( UPM_Status .eq. UPM_INVNUM) then
                       Ret_Status = %loc(FEC_INVTOL)
                       call lib$signal (FEC_INVTOL)
                       FEC_Get_Quals = Ret_Status
                       return
                    end if
                    if (tsize(j) .gt. 0.0) then
                       temp_size(j)=tsize(j)
                       write(temp_text(j),99)tsize(j)
  99                   Format(F10.4)
                       temp_len(j) = 6
                       if(tsize(j) .ge. 10.)temp_len(j)=7
                       if(tsize(j) .ge. 100.)temp_len(j)=8
                       if(tsize(j) .ge. 1000.)temp_len(j)=9
                       if(tsize(j) .ge. 10000.)temp_len(j)=10
                    endif
                 enddo
                 if( Ret_Status ) then
                    if((UPM_Status .eq. UPM_Concat).or.(UPM_Status .eq.
     .              UPM_Comma)) then
                       Ret_Status = %loc(FEC_TOOMANYTOL)
                       call lib$signal (FEC_TOOMANYTOL,%val(Num_Temps_Used))
                       FEC_Get_Quals = Ret_Status
                       return
                    end if
                 end if

             END IF

             command_line5(1:16)='    /TEMP_SIZE=('
             clen=16
             j = 0
             do while ( (j .lt. Num_Temps_Used) .and. Ret_Status)
                j = j + 1
                ttext=temp_text(j)
                command_line5(clen+1:clen+temp_len(j))=ttext(12-temp_len(j):11)
                clen=clen+temp_len(j)
                if ( j .lt. Num_Temps_Used ) then
                   command_line5(clen+1:clen+1)=','
                   clen=clen+1
                else
                   command_line5(clen+1:clen+1)=')'
                   clen=clen+1
                end if
             enddo
              
           ELSE if(i .eq. 15)then

              command_line5(clen+1:clen+11)='/MIN_DIHED='
              clen=clen+11

c             << Min_DiHed >>

              IF (status .EQ. upm_absent .or. status .EQ. upm_defaulted) THEN

                 min_dihed = 0.0
                 command_line5(clen+1:clen+4)='0.00'
                 clen=clen+4
   
              ENDIF

              IF (status .EQ. upm_pres) THEN

                 status = upm_get_float(qualifier(i), frsize )
                 min_dihed = frsize
                 write(trsize,50)frsize
 50              Format(F4.2)
                 command_line5(clen+1:clen+4)=trsize
                 clen=clen+4

              ENDIF

           ELSE if(i .eq. 16)then

              command_line5(clen+1:clen+11)='/MAX_DIHED='
              clen=clen+11

c             << Max_DiHed >>

              IF (status .EQ. upm_absent .or. status .EQ. upm_defaulted) THEN

                 max_dihed = 3.5
                 command_line5(clen+1:clen+4)='3.50'
                 clen=clen+4
   
              ENDIF

              IF (status .EQ. upm_pres) THEN

                 status = upm_get_float(qualifier(i), frsize )
                 max_dihed = frsize
                 write(trsize,50)frsize
                 command_line5(clen+1:clen+4)=trsize
                 clen=clen+4

              ENDIF

           END IF

	END DO
c
c

        Command_Line(1)=Command_Line1
        Command_Line(2)=Command_Line2
        Command_Line(3)=Command_Line3
        Command_Line(4)=Command_Line4
        Command_Line(5)=Command_Line5
        Command_Line(6)=Command_Line6

c
c

        FEC_Get_Quals = Ret_Status
	RETURN
	END
