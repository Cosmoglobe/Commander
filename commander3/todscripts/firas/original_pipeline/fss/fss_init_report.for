      integer*4 function fss_init_report(repfile, lu_rep, c_array,
     .                                   qual_stat, cnum, rse, rsefile)

C------------------------------------------------------------------------
C    PURPOSE: This function initializes the report. It obtains the logical
C             unit number for the report and the report file name, 
C             opens the report file, and writes the header information.
C
C    AUTHOR: D. Bouler, STX, Jan, 1990
C
C    INVOCATION: status = fss_init_report(repfile,lu_rep,c_array,
C                                         qual_stat,cnum,rse, rsefile)
C
C    INPUT PARAMETERS:   C*128 repfile         report file name
C                        C*79  c_array(30)     command line
C                        I*4   qual_stat       status from fss_get_quals 
C                        I*4   cnum            number of command line elements
C                        c*128 rse(16)         RSE lines
C                        c*128 rsefile         RSE file name
C
C    OUTPUT PARAMETERS:  I*4   lu_rep          logical unit for report
C                        I*4   fut_report_lun  Logical unit for fut_error
C
C    SUBROUTINES CALLED:    cut_register_version
C                           cut_display_banner
C                           str$upcase
C
C    COMMON VARIABLES USED: fut_report_lun
C 
C    INCLUDE FILES:         $jpidef
C                           fss_include
C                           fut_params
C                           fut_error
C
C----------------------------------------------------------------------
C
C Changes:
C	SPR 9943 - Report file includes CSDR$ archive references which
C			are not accessed by FSS.
C
C
C----------------------------------------------------------------------
C                  PDL FOR fss_init_report
C
C                  D. BOULER, JAN, 1991
C
C fss_init_report  = %loc(fss_normal)         
C status           = %loc(fss_normal)
C
C get today's date and time
C get process id of person running FSS
C
C if (status) then
C     get LUN for report
C     if (.not. status) signal error
C endif
C
C if (status) then
C     open report file
C     if (.not. status) then
C          signal error
C     else
C          set fut_report_lun to report file lun
C          call lib$establish to set up FUT_ERROR error handler
C     endif
C endif
C
C if (status) then
C     write report banner and user ID to report and log file 
C     if (.not. status) signal error
C endif
C
C if (status) then 
C     translate FIRAS archive logicals and write to report and log file
C     write command line to report and log file
C endif
c
c if (status) then
c     if rse file name has a logical in it, translate and report.
c endif
c
C if (status) then 
C     write RSEs to report and log file
C endif
c
c if (.not. qual_stat) then
c      signal previous error from fss_get_quals so it is in report file
C      status = qual_stat
C endif
C
C fss_init_report = status 
C
C return
C end (pdl)
C*****************************************************************************
C
      implicit none
c
c Include files:
c
      include  '($jpidef)'
      include  '(fss_include)'
      include  '(fut_params)'
      include  '(fut_error)'
c
c Declare externals
c
      external fss_normal
      external fut_error
      external fss_badrepopen
c
c Function declarations:
c
      integer*4 cut_translate_archive_id
      integer*4 cut_register_version
      integer*4 cut_display_banner
      integer*4 lib$get_lun
C
C Variable declarations
C
      integer*4 i,j,k,l                ! Loop counters
      integer*4 qual_stat              ! Status from fss_get_quals
      integer*4 status                 ! General status return
      integer*4 t_stat                 ! Status return for translation
      integer*4 ios                    ! i/o status return
      integer*4 lu_rep                 ! Logical unit for report
      integer*4 len                    ! length of report file name
      integer*4 today                  ! Current date and time (binary)
      integer*4 name_given             ! Was report name specified ?
      integer*4 report                 ! Was report requested ?
      integer*4 temp_len               ! Temporary translated logical length
      integer*4 trans_len(4)           ! Length of translated logicals
      integer*4 cnum                   ! Number of command line array elements
      integer*4 colon                  ! Position of ":" in file name

      character*20   acct              ! Account of person running FSS

      character*24   atoday            ! Ascii date and time

      character*79   in_logs(4)        ! Firas logical name s to translate
      character*79   temp_log          ! Temporary logical name
      character*79   trans_logs(4)     ! Translated logical names
      character*79   c_array(30)       ! Command line
      character*79   rse_log           ! RSE logical in file name 

      character*128  repfile           ! Report file name
      character*128  inskymap          ! input skymap name
      character*128  outskymap         ! output skymap name
      character*128  rse(16)           ! RSEs
      character*128  rsefile           ! RSE file name
c*
c******  Begin code *******************************************
c*
      fss_init_report = %loc(fss_normal)         
      status          = %loc(fss_normal)         
c
c Get today's date and the user running FSS.
c
      call sys$asctim(, atoday, today)
      call lib$getjpi(jpi$_username,,,,acct,)
c
c Get logical unit number for report
c
      if (status) then
          status = lib$get_lun(lu_rep)
          if (.not. status) call lib$signal(%val(status))
      endif
c
c Open report file 
c
      if (status) then
          open(lu_rep,file=repfile,iostat=ios,status='new')
          if (ios .ne. 0) then
              status = %loc(fss_badrepopen)
              call lib$signal(fss_badrepopen, %val(1), %val(ios) )
          else
              fut_report_lun = lu_rep
              call lib$establish(fut_error)
          endif
      endif
c
c Write banner and account.
c
      if (status) then
          call Cut_Register_Version(fss_version)
          call Cut_Display_Banner(lu_rep,80,'FIRAS Facility FSS_Sort_Sky')
          Write(lu_rep,fmt='(/)')
          write(lu_rep,fmt='(x,2a)') 'Run by: ',acct
      endif     
c
c Set up names of logicals to be translated
c
      if (status) then
          in_logs(1) = 'CSDR$FIRAS_IN'
          in_logs(2) = 'CSDR$FIRAS_OUT'
          in_logs(3) = 'CSDR$FIRAS_RAW'
          in_logs(4) = 'CSDR$FIRAS_REF'
      endif
c
c Get translations for all software logicals
c
      if (status) then
          do i = 1,4
             t_stat = cut_translate_archive_id(in_logs(i), temp_log, temp_len,
     .                                         trans_logs(i), trans_len(i))
             if (.not. t_stat) call lib$signal(%val(t_stat))
          enddo
      endif
c
c Write translated logicals to report and log file
c
      if (status) then
          do j = 1,4
              call str$trim(trans_logs(j), trans_logs(j), trans_len(j) )
              k = min(50, trans_len(j))
              write(lu_rep,fmt='(x,a,t25,a)') in_logs(j)(:20),trans_logs(j)(:k)
          enddo
      endif
c
c If RSE file name has a logical, translate and report.
c
      if (status) then
          colon = index(rsefile, ':')
          if (colon .gt. 0) then
              rse_log = rsefile(:colon-1)  
              t_stat = cut_translate_archive_id(rse_log, temp_log, temp_len,
     .                                          trans_logs(1), trans_len(1))
              if (t_stat) then
                  call str$trim(rsefile, rsefile, len)
                  rsefile = trans_logs(1)(:trans_len(1)) // rsefile(colon:len)
                  write(lu_rep,fmt='(x,/,x,a,t25,a,/)') 'RSE file:',rsefile(:50)
              endif
           endif
      endif
c
c Write RSEs to report and log file
c
      if (status) then
          write(lu_rep,fmt='(x,/,x,a)') 'Record selection expressions:'
          do i = 1,16
             if (rse(i)(1:1) .ne. ' ') write(lu_rep,*) rse(i)(:75) 
          enddo
      endif
c
c Write command line to report and log file
c
      if (status) then
          write(6     ,fmt='(x,/x,a)') 'Command line:'
          write(lu_rep,fmt='(x,/x,a)') 'Command line:'
          do j = 1,cnum
             write(6     ,fmt='(x,a)') c_array(j)
             write(lu_rep,fmt='(x,a)') c_array(j)
          enddo
      endif
      write(6,*) ' '
c
c If status from fss_get_quals is bad, signal the error.
c Write lines to indicate this error came from fss_get_quals.
c Set return status for fss_init_report to the bad return from get_quals.
c
      if (.not. qual_stat) then
           write(6     ,*) ' '
           write(lu_rep,*) ' '

           write(6,     10) 'Error from FSS_GET_QUALS'
           write(lu_rep,10) 'Error from FSS_GET_QUALS'
10         format(x,25('*'),a,25('*')) 

           call lib$signal(%val(qual_stat))

           write(6,     20)
           write(lu_rep,20)
20         format(x,74('*')) 

           status = qual_stat
      endif

      fss_init_report = status 

      return
      end
