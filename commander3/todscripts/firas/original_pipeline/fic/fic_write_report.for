      integer*4 function fic_write_report(channel,scan_mode,deglitch,
     .                                    report_detail,report_name,account,
     .                                    current_gmt,mincoadd,cth,num_recs,
     .                                    fake_it,sci_mode,adds_per_group,
     .                                    con_recs,first_sci_rec,last_sci_rec,
     .                                    avg_noise,glitch_pos,glitch_amp,
     .                                    report_luns,report_filenames)
c-------------------------------------------------------------------------------
c
c     Purpose: Write the consistency/deglitcher report section for a coadd 
c              group.
c
c     Author: S. Alexander, STX, 9/91, SER 7985
c
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode          i*4  Value of scan mode, 1-4=SS-LF.
c            deglitch           i*4  Perform deglitching.
c            report_detail      i*4  Detailed report option.
c            report_name        ch*72  Consistency and deglitch report name on
c                                      command line.
c            account            ch*12  Account owner name.
c            current_gmt        ch*14  Current time in Julian format.
c            mincoadd           rec  Minimum interferograms to coadd.
c            cth                rec  Consistency check parameters.
c            num_recs           i*4  Number of input science and engineering
c                                    records.
c            fake_it            i*4  Value of fake-it bit, 0-1.
c            sci_mode           i*4  Value of science mode, 0-4.
c            adds_per_group     i*4  Value of adds per group, 1-12.
c            con_recs           rec(num_recs)  Consistency check records.
c            first_sci_rec      i*4  Pointer to earliest science record in
c                                    coadd group.
c            last_sci_rec       i*4  Pointer to latest science record in 
c                                    coadd group.
c            avg_noise          r*4  Average noise of coadd group.
c            glitch_pos         i*4(1000,num_recs)  Glitch positions.
c            glitch_amp         r*4(1000,num_recs)  Glitch amplitudes.
c            report_luns        i*4(4)  Logical units for reports by scan mode.
c            report_filenames   ch*72(4,4)  Report filenames by channel and
c                                           scan mode.
c
c     Output: report_luns        i*4(4)  Logical units for reports by scan mode.
c             report_filenames   ch*72(4,4)  Report filenames by channel and
c                                            scan mode.
c
c     Modifications:
c
c-----------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_params)'
c
c     Return statuses.
c
      external fic_normal
      external fic_repwrite
c
c     Functions.
c
      integer*4 fic_open_report
c
c     Input parameters.
c
      integer*4 channel,scan_mode,deglitch,report_detail
      character*72 report_name
      character*12 account
      character*14 current_gmt

      dictionary 'fex_mincoadd'
      dictionary 'fex_cth'
      record /fex_mincoadd/ mincoadd
      record /fex_cth/ cth

      integer*4 num_recs,first_sci_rec,last_sci_rec
      integer*4 fake_it,sci_mode,adds_per_group

      dictionary 'fic_scc'
      record /fic_scc/ con_recs(num_recs)

      real*4 avg_noise
      integer*4 glitch_pos(fac_max_deglitch,num_recs)
      real*4 glitch_amp(fac_max_deglitch,num_recs)
c
c     Input/output parameters.
c
      integer*4 report_luns(4)
      character*72 report_filenames(4,4)
c
c     Local variables.
c
      integer*4 lun,status,pos,rec
      character*72 filename
      integer*4 file_len
      integer*4 num_good,num_fail(32)
      logical*1 first
      integer*4 fail,bit,iter
      integer*4 percent_good,percent_fail(32)

      fic_write_report = %loc(fic_normal)
c
c     If a report has not yet been opened for this channel and scan mode, call
c     fic_open_report to do so.
c
      lun = report_luns(scan_mode)

      if (lun .eq. 0) then

         status = fic_open_report(channel,scan_mode,report_name,account,
     .                            current_gmt,mincoadd,cth,lun,report_filenames)
         if (status .ne. %loc(fic_normal)) then
            fic_write_report = status
            return
         end if

         report_luns(scan_mode) = lun
      end if
c
c     Write information about the coadd group.
c
      write(lun,10,iostat=status) con_recs(1).coadd_gmt,num_recs,
     .                            con_recs(first_sci_rec).gmt,
     .                            con_recs(last_sci_rec).gmt,
     .                            fake_it,sci_mode,adds_per_group,avg_noise
10    format(/x,'Coadd Time    Num Recs First Sci Time  Last Sci Time Fake ',
     .      'SciMode Adds Noise',/x,a,3x,i3,3x,a,2x,a,2x,i2,3x,i2,4x,i2,2x,g9.3)
      if (status .ne. 0) then
         fic_write_report = %loc(fic_repwrite)
         inquire(lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fic_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if
c
c     Accumulate number of good and failing interferograms.
c
      num_good = 0
      do pos = 1,32
         num_fail(pos) = 0
      end do
      first = .true.

      do rec = 1,num_recs
         if (con_recs(rec).con_check .eq. 0) then
            fail = 0
            num_good = num_good + 1
         else
            fail = 0
            bit = 0
            do while ((fail .eq. 0) .and. (bit .le. 31))
               if (btest(con_recs(rec).con_check,bit)) then
                  fail = bit + 1
                  num_fail(bit+1) = num_fail(bit+1) + 1
               end if
               bit = bit + 1
            end do
         end if
c
c     If deglitching was performed, write summary deglitcher information
c     for the interferogram.
c
         if (deglitch .eq. fac_present) then
            if (first) then
               first = .false.
               write(lun,20,iostat=status)
20             format(45x,'First Glitch     Total Glitch',
     .                /x,'Record Time    Pixel  Fail Reason  Noise     ',
     .                'Pos  Amp        Iter  Signal')
               if (status .ne. 0) then
                  fic_write_report = %loc(fic_repwrite)
                  inquire(lun,name=filename)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fic_repwrite,%val(2),filename(:file_len),
     .                            %val(status))
                  return
               end if
            end if
            write(lun,30,iostat=status) con_recs(rec).gmt,
     .                                  con_recs(rec).pixel_no,fail,
     .                                  con_recs(rec).noise,
     .                                  con_recs(rec).glitch_pos,
     .                                  con_recs(rec).glitch_amp,
     .                                  con_recs(rec).glitch_iter,
     .                                  con_recs(rec).glitch_signal
30          format(x,a,2x,i4,7x,i2,4x,g9.3,3x,i3,2x,g10.3,x,i4,2x,g10.3)
            if (status .ne. 0) then
               fic_write_report = %loc(fic_repwrite)
               inquire(lun,name=filename)
               call str$trim(filename,filename,file_len)
               call lib$signal(fic_repwrite,%val(2),filename(:file_len),
     .                         %val(status))
               return
            end if
c
c     If a detailed report was requested, give information for each glitch.
c
            if (report_detail .eq. fac_present) then
               do iter = 2,con_recs(rec).glitch_iter
                  write(lun,40,iostat=status) iter,glitch_pos(iter,rec),
     .                                        glitch_amp(iter,rec)
40                format(x,'Glitch ',i4,'; Pos ',i3,'; Amp ',g10.3)
                  if (status .ne. 0) then
                     fic_write_report = %loc(fic_repwrite)
                     inquire(lun,name=filename)
                     call str$trim(filename,filename,file_len)
                     call lib$signal(fic_repwrite,%val(2),filename(:file_len),
     .                               %val(status))
                     return
                  end if
               end do
            end if
         else if (con_recs(rec).con_check .ne. 0) then
c
c     If deglitching was not performed, write information for each 
c     interferogram without deglitching data.
c
            if (first) then
               first = .false.
               write(lun,50,iostat=status)
50             format(/x,'Record Time    Pixel  Fail Reason  Noise')
               if (status .ne. 0) then
                  fic_write_report = %loc(fic_repwrite)
                  inquire(lun,name=filename)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fic_repwrite,%val(2),filename(:file_len),
     .                            %val(status))
                  return
               end if
            end if
            write(lun,60,iostat=status) con_recs(rec).gmt,
     .                                  con_recs(rec).pixel_no,fail,
     .                                  con_recs(rec).noise
60          format(x,a,2x,i4,7x,i2,4x,g9.3)
            if (status .ne. 0) then
               fic_write_report = %loc(fic_repwrite)
               inquire(lun,name=filename)
               call str$trim(filename,filename,file_len)
               call lib$signal(fic_repwrite,%val(2),filename(:file_len),
     .                         %val(status))
               return
            end if
         end if
      end do
c
c     Calculate and write statistics for passing and failing interferograms.
c
      percent_good = nint(num_good * 100.0 / num_recs)

      do pos = 1,32
         percent_fail(pos) = nint(num_fail(pos) * 100.0 / num_recs)
      end do

      write(lun,70,iostat=status) (num_fail(pos), pos = 1,16),
     .                            (percent_fail(pos), pos = 1,16),
     .                            (num_fail(pos), pos = 17,32),
     .                            (percent_fail(pos), pos = 17,32)
70    format(/x,'Number and Percent of Records Failing for Each Reason:',/x,
     .       16(i4,x),/x,16(i3,'% '),/x,16(i4,x),/x,16(i3,'% '))
      if (status .ne. 0) then
         fic_write_report = %loc(fic_repwrite)
         inquire(lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fic_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if

      write(lun,80,iostat=status) num_good,percent_good
80    format(/x,'Number and Percent of Good Records: ',i3,'; ',i3,'%')
      if (status .ne. 0) then
         fic_write_report = %loc(fic_repwrite)
         inquire(lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fic_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if

      return
      end
