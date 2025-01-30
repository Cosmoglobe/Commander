      integer*4 function fcl_report(channel,scan_mode,cov_ext,label,
     .                              command_line,parse_status,in_cov_name,
     .                              out_cov_name)
c-------------------------------------------------------------------------------
c
c     Purpose: Open and write processing report.
c
c     Author: S. Brodd, HSTX, 12/95, SPR 12291
c
c     Input: channel       i*4  Value of channel.
c            scan_mode     i*4  Value of scan mode.
c            cov_ext       ch*72  Value of input filename extension.
c            label         l*1  Label present or not.
c            command_line  ch*79(7)  Fully defaulted command line.
c            parse_status  i*4  Return status from fcl_parse.
c
c     Output: in_cov_name   ch*72  Filename for input covariance matrix.
c             out_cov_name  ch*72  Filename for output covariance matrix.
c
c     Modifications:
c    
c-------------------------------------------------------------------------------
      implicit none
c
c     Include files.
c
      include '(fut_error)'
      include '($jpidef)'
      include '(fut_params)'
c
c     Return statuses.
c
      external fcl_normal
      external fcl_repopen,fcl_repwrite
c
c     External references.
c
      external fut_error,fcl_parseerr
c
c     Functions.
c
      integer*4 fut_get_lun
      integer*4 cut_display_banner
      integer*4 cut_translate_archive_id
c
c     Input parameters.
c
      integer*4 channel,scan_mode,parse_status
      character*72 cov_ext
      logical*1 label
      character*79 command_line(7)
c
c     Output parameters.
c
      character*72 in_cov_name,out_cov_name
c
c     Local variables.
c
      integer*4 current_time(2),len,status,pos,int_len,trans_len(3)
      character*14 current_gmt
      character*72 report,logs(3),int_log,trans_logs(3)
      character*12 account

      fcl_report = %loc(fcl_normal)         
c
c     Get current time.
c
      call sys$gettim(current_time)
      call ct_binary_to_gmt(current_time,current_gmt)

      call str$trim(cov_ext,cov_ext,len)
c
c     Create input and output covariance matrix filenames.
c
      in_cov_name = 'CSDR$FIRAS_IN:FIL_COV_'//fac_channel_ids(channel)//
     .               fac_scan_mode_idsl(scan_mode)//'.'//cov_ext(:len)
      out_cov_name = 'CSDR$FIRAS_OUT:FCL_COV_'//fac_channel_ids(channel)//
     .               fac_scan_mode_idsl(scan_mode)//'.'//cov_ext(:len)
c
c     Create report filename.
c
      report = 'FCL_'//fac_channel_ids(channel)//fac_scan_mode_idsl(scan_mode)
     .          //'_'//cov_ext(4:len)//'.REP_'//current_gmt(1:9)

      call str$trim(report,report,len)
c
c     Open the report file.
c
      status = fut_get_lun(fut_report_lun)
      if (.not. status) then
         fcl_report = status
         return
      end if

      open(unit=fut_report_lun,file=report,status='new',iostat=status)
      if (status .ne. 0) then
         fcl_report = %loc(fcl_repopen)
         call lib$signal(fcl_repopen,%val(2),report(:len),%val(status))
         return
      else
         call lib$establish(fut_error)
      end if
c
c     Write banner, report name, account, and current time.
c
      status = cut_display_banner(fut_report_lun,80,
     .                     'FIRAS Facility FCL_CalibrateCovariances_Long')

      write(fut_report_lun,10,iostat=status) 'Report File: ',
     .                                       report(:len)
10    format(/x,2a)
      if (status .ne. 0) then
         fcl_report = %loc(fcl_repwrite)
         call lib$signal(fcl_repwrite,%val(2),report(:len),
     .                   %val(status))
         return
      end if

      call lib$getjpi(jpi$_username,,,,account,)

      write(fut_report_lun,20,iostat=status) 'Account: ',account,
     .                                       '    Time: ',current_gmt(1:9)
20    format(/x,4a)
      if (status .ne. 0) then
         fcl_report = %loc(fcl_repwrite)
         call lib$signal(fcl_repwrite,%val(2),report(:len),
     .                   %val(status))
         return
      end if
c
c     If fcl_parse returned an error, signal it error now so that it is 
c     written to the report file.
c
      if (parse_status .eq. %loc(fcl_parseerr)) then
         fcl_report = %loc(fcl_parseerr)
         call lib$signal(fcl_parseerr)
         return
      end if
c
c     Assign names of logicals to be translated.
c
      logs(1) = 'CSDR$FIRAS_IN'
      logs(2) = 'CSDR$FIRAS_OUT'
      logs(3) = 'CSDR$FIRAS_REF'
c
c     Get translations for logical names and write them to the report.
c
      write(fut_report_lun,30,iostat=status) 'Logical Name Translations:'
30    format(/x,a)
      if (status .ne. 0) then
         fcl_report = %loc(fcl_repwrite)
         call lib$signal(fcl_repwrite,%val(2),report(:len),
     .                   %val(status))
         return
      end if

      do pos = 1,3
         status = cut_translate_archive_id(logs(pos),int_log,int_len,
     .                                     trans_logs(pos),trans_len(pos))
         call str$trim(trans_logs(pos),trans_logs(pos),trans_len(pos))

         write(fut_report_lun,40,iostat=status) logs(pos)(:15),'   ',
     .                              trans_logs(pos)(:trans_len(pos))
40       format(x,3a)
         if (status .ne. 0) then
            fcl_report = %loc(fcl_repwrite)
            call lib$signal(fcl_repwrite,%val(2),report(:len),
     .                      %val(status))
            return
         end if
      end do
c
c     Write fully defaulted command line to the report.
c
      write(fut_report_lun,30,iostat=status) 'Command Line:'
      if (status .ne. 0) then
         fcl_report = %loc(fcl_repwrite)
         call lib$signal(fcl_repwrite,%val(2),report(:len),
     .                   %val(status))
         return
      end if

      do pos = 1,6
         write(fut_report_lun,60,iostat=status) command_line(pos)
60       format(x,a)
         if (status .ne. 0) then
            fcl_report = %loc(fcl_repwrite)
            call lib$signal(fcl_repwrite,%val(2),report(:len),
     .                      %val(status))
            return
         end if
      end do

      if (label) then
         write(fut_report_lun,60,iostat=status) command_line(7)
         if (status .ne. 0) then
            fcl_report = %loc(fcl_repwrite)
            call lib$signal(fcl_repwrite,%val(2),report(:len),
     .                      %val(status))
            return
         end if
      end if
c
c     Write input and output covariance matrix filenames to the report.
c
      write(fut_report_lun,30,iostat=status) 'Input Filename:'
      if (status .ne. 0) then
         fcl_report = %loc(fcl_repwrite)
         call lib$signal(fcl_repwrite,%val(2),report(:len),
     .                   %val(status))
         return
      end if

      write(fut_report_lun,60,iostat=status) in_cov_name
      if (status .ne. 0) then
         fcl_report = %loc(fcl_repwrite)
         call lib$signal(fcl_repwrite,%val(2),report(:len),
     .                   %val(status))
         return
      end if

      write(fut_report_lun,30,iostat=status) 'Output Filename:'
      if (status .ne. 0) then
         fcl_report = %loc(fcl_repwrite)
         call lib$signal(fcl_repwrite,%val(2),report(:len),
     .                   %val(status))
         return
      end if

      write(fut_report_lun,60,iostat=status) out_cov_name
      if (status .ne. 0) then
         fcl_report = %loc(fcl_repwrite)
         call lib$signal(fcl_repwrite,%val(2),report(:len),
     .                   %val(status))
         return
      end if

      return
      end
