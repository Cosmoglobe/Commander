      integer*4 function fic_open_summary(file,file_ext,jstart,jstop,rse_input,
     .                                    rse_name,report,report_name,
     .                                    command_line,parse_status,account,
     .                                    current_gmt)
c-------------------------------------------------------------------------------
c
c     Purpose: Open and initialize summary report.
c
c     Author: S. Alexander, STX, 9/91, SER 7985
c
c     Input: file               i*4  Input filename specified on command line.
c            file_ext           ch*72  Value of input filename extension on
c                                      command line.
c            jstart             ch*14  Value of start time.
c            jstop              ch*14  Value of stop time.
c            rse_input          i*4  Rse specified on command line.
c            rse_name           ch*72  Rse filename.
c            report             i*4  Write consistency and deglitch reports.
c            report_name        ch*72  Consistency and deglitch report name on
c                                      command line.
c            command_line       ch*79(9)  Fully defaulted command line.
c            parse_status       i*4  Return status from fic_parse.
c
c     Output: report_name        ch*72  Consistency and deglitch report name on
c                                       command line or default.
c             account            ch*12  Account owner name.
c             current_gmt        ch*14  Current time in Julian format.
c
c     Modifications:
c    
c-------------------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_error)'
      include '(fut_params)'
      include '($jpidef)'
c
c     Return statuses.
c
      external fic_normal
      external fic_repopen
      external fic_repwrite
      external fic_invtime
c
c     External references.
c
      external fut_error
c
c     Functions.
c
      integer*4 fut_get_lun
      integer*4 cut_display_banner
      integer*4 cut_translate_archive_id
c
c     Input parameters.
c
      integer*4 file,rse_input,report,parse_status
      character*72 file_ext,rse_name
      character*14 jstart,jstop
      character*79 command_line(9)
c
c     Input/output parameters.
c
      character*72 report_name
c
c     Output parameters.
c
      character*12 account
      character*14 current_gmt
c
c     Local variables.
c
      integer*4 status,num_logs,pos,file_ext_len,period
      integer*4 current_time(2)
      integer*4 rep_len
      character*72 logs(7),int_log,trans_logs(7)
      integer*4 rse_colon,rep_colon
      integer*4 int_len,trans_len(7)

      fic_open_summary = %loc(fic_normal)         
c
c     Get current time.
c
      call sys$gettim(current_time)
      call ct_binary_to_gmt(current_time,current_gmt)
c
c     If no input filename extension was specified on the command line, create
c     a default for the report filenames.
c
      if (file .eq. fac_not_present) then
         file_ext = jstart(1:7)//'_'//jstop(1:7)
         file_ext_len = 15
      else
         call str$trim(file_ext,file_ext,file_ext_len)
      end if
c
c     If a report filename was specified on the command line, use it.  
c     Otherwise, construct a default.
c
      call str$trim(report_name,report_name,rep_len)

      if (rep_len .ne. 0) then
         period = index(report_name,'.')
c
c     Create a report filename extension if none was specified.
c
         if (period .eq. 0) then
            report_name = report_name(:rep_len)//'.REP_'//current_gmt(1:9)
         end if
      else
         report_name = 'FIC_'//file_ext(:file_ext_len)//
     .                  '.REP_'//current_gmt(1:9)
      end if

      call str$trim(report_name,report_name,rep_len)
c
c     Open the summary report file.
c
      status = fut_get_lun(fut_report_lun)
      if (.not. status) then
         fic_open_summary = status
         return
      end if

      open(unit=fut_report_lun,file=report_name,status='new',iostat=status)
      if (status .ne. 0) then
         fic_open_summary = %loc(fic_repopen)
         call lib$signal(fic_repopen,%val(2),report_name(:rep_len),
     .                   %val(status))
         return
      else
         call lib$establish(fut_error)
      endif
c
c     Write banner, summary name, account, and current time.
c
      status = cut_display_banner(fut_report_lun,80,
     .                            'FIRAS Facility FIC_Interferogram_Coadd')

      write(fut_report_lun,10,iostat=status) 'Summary Report File: ',
     .                                       report_name(:rep_len)
10    format(/x,2a)
      if (status .ne. 0) then
         fic_open_summary = %loc(fic_repwrite)
         call lib$signal(fic_repwrite,%val(2),report_name(:rep_len),
     .                   %val(status))
         return
      end if

      call lib$getjpi(jpi$_username,,,,account,)

      write(fut_report_lun,20,iostat=status) 'Account: ',account,
     .                                       '    Time: ',current_gmt(1:9)
20    format(/x,4a)
      if (status .ne. 0) then
         fic_open_summary = %loc(fic_repwrite)
         call lib$signal(fic_repwrite,%val(2),report_name(:rep_len),
     .                   %val(status))
         return
      end if
c
c     If jstop was less than jstart in fic_parse, signal the error now
c     so that it is written to the summary report file.
c
      if (parse_status .eq. %loc(fic_invtime)) then
         fic_open_summary = %loc(fic_invtime)
         call lib$signal(fic_invtime,%val(2),jstop,jstart)
         return
      endif
c
c     Assign names of logicals to be translated.
c
      num_logs = 5
      logs(1) = 'CSDR$FIRAS_IN'
      logs(2) = 'CSDR$FIRAS_RAW'
      logs(3) = 'CSDR$FIRAS_OUT'
      logs(4) = 'CSDR$FIRAS_STAT'
      logs(5) = 'CSDR$FIRAS_REF'
c
c     Translate logical names in rse file name and report file name, if they
c     exist.
c
      if (rse_input .eq. fac_present) then
         rse_colon = index(rse_name,':')
         if (rse_colon .gt. 0) then
            num_logs = 6
            logs(6) = rse_name(:rse_colon-1)  
         end if
      end if

      if (report .eq. fac_present) then
         rep_colon = index(report_name,':')
         if (rep_colon .gt. 0) then
            num_logs = 7
            logs(7) = report_name(:rep_colon-1)  
         end if
      end if
c
c     Get translations for logical names and write them to summary report.
c
      write(fut_report_lun,30,iostat=status) 'Logical Name Translations:'
30    format(/x,a)
      if (status .ne. 0) then
         fic_open_summary = %loc(fic_repwrite)
         call lib$signal(fic_repwrite,%val(2),report_name(:rep_len),
     .                   %val(status))
         return
      end if

      do pos = 1,num_logs
         status = cut_translate_archive_id(logs(pos),int_log,int_len,
     .                                     trans_logs(pos),trans_len(pos))
         call str$trim(trans_logs(pos),trans_logs(pos),trans_len(pos))

         if (pos .le. 5) then
            write(fut_report_lun,40,iostat=status) logs(pos)(:15),'   ',
     .                              trans_logs(pos)(:trans_len(pos))
40          format(x,3a)
            if (status .ne. 0) then
               fic_open_summary = %loc(fic_repwrite)
               call lib$signal(fic_repwrite,%val(2),report_name(:rep_len),
     .                         %val(status))
               return
            end if
         end if
      end do
c
c     Write translations of report file logical and rse file if they exist.
c
      if ((report .eq. fac_present) .and. (rep_colon .gt. 0)) then
         write(fut_report_lun,50,iostat=status) 'Report File: ',
     .             trans_logs(7)(:trans_len(7))//report_name(rep_colon:rep_len)
50    format(x,2a)
         if (status .ne. 0) then
            fic_open_summary = %loc(fic_repwrite)
            call lib$signal(fic_repwrite,%val(2),report_name(:rep_len),
     .                      %val(status))
            return
         end if
      end if

      if (rse_input .eq. fac_present) then
         call str$trim(rse_name,rse_name,int_len)
         if (rse_colon .gt. 0) then
            write(fut_report_lun,10,iostat=status) 'RSE File: ',
     .            trans_logs(6)(:trans_len(6))//rse_name(rse_colon:int_len)
         else
            write(fut_report_lun,10,iostat=status) 'RSE File: ',
     .                                             rse_name(:int_len)
         end if
         if (status .ne. 0) then
            fic_open_summary = %loc(fic_repwrite)
            call lib$signal(fic_repwrite,%val(2),report_name(:rep_len),
     .                      %val(status))
            return
         end if
      end if
c
c     Write fully defaulted command line to summary report.
c
      write(fut_report_lun,30,iostat=status) 'Command Line:'
      if (status .ne. 0) then
         fic_open_summary = %loc(fic_repwrite)
         call lib$signal(fic_repwrite,%val(2),report_name(:rep_len),
     .                   %val(status))
         return
      end if

      do pos = 1,9
         write(fut_report_lun,60,iostat=status) command_line(pos)
60       format(x,a)
         if (status .ne. 0) then
            fic_open_summary = %loc(fic_repwrite)
            call lib$signal(fic_repwrite,%val(2),report_name(:rep_len),
     .                      %val(status))
            return
         end if
      end do

      return
      end
