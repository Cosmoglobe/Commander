      integer*4 function fil_close_summary(covar,deglitch,output,report,
     .                                     mincoadd,cth,num_sci_recs,num_groups,
     .                                     num_sci_good,num_sci_fail,
     .                                     num_sci_noise,tot_noise,
     .                                     num_sci_deglitch,glitch_tot_iter,
     .                                     glitch_tot_signal,num_coadds,
     .                                     output_filenames,report_filenames)
c-------------------------------------------------------------------------------
c
c     Purpose: Complete summary report.
c
c     Author: S. Brodd, HSTX, 4/95
c
c     Input: covar              i*4  Covariance matrix to be written.
c            deglitch           i*4  Perform deglitching.
c            output             i*4(6) Selection of output types. 
c            report             i*4  Write consistency and deglitch reports.
c            mincoadd           rec  Minimum interferograms to coadd.
c            cth                rec  Consistency check parameters.
c            num_sci_recs       i*4(4,4)  Number of input science and
c                                         engineering records by channel and 
c                                         scan mode.
c            num_groups         i*4(4,4)  Number of input coadd groups by
c                                         channel and scan mode.
c            num_sci_good       i*4(4,4,20)  Number of good interferogrmas
c                                            contributing to an output coadd 
c                                            record by channel and scan mode; 
c                                            grouped as 1-5, ..., 96-100.
c            num_sci_fail       i*4(4,4,32)  Number of science or engineering
c                                            records failing by channel and 
c                                            scan mode for each of 32 
c                                            consistency checks.
c            num_sci_noise      i*4(4,4)  Number of interferograms for which a
c                                         noise was calculated by channel and 
c                                         scan mode.
c            tot_noise          r*4(4,4)  Total noise by channel and scan mode.
c            num_sci_deglitch   i*4(4,4)  Number of deglitched interferograms
c                                         by channel and scan mode.
c            glitch_tot_iter    i*4(4,4)  Total deglitcher iterations by
c                                         channel and scan mode.
c            glitch_tot_signal  r*4(4,4)  Total deglitcher signal by channel 
c                                         and scan mode.
c            num_coadds         i*4(4,6)  Number of output coadd records by
c                                         channel and scan mode, including
c                                         FS and FL.
c            output_filenames   ch*72(4,6,7)  Output filenames by channel,
c                                             scan mode, including FS and FL,
c                                             and output type, including 
c                                             covariance matrices.
c            report_filenames   ch*72(4,4)  Report filenames by channel and
c                                           scan mode.
c
c     Output:
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
c
c     Return statuses.
c
      external fil_normal
      external fil_repwrite
c
c     Input parameters.
c
      integer*4 covar,deglitch,output(6),report

      dictionary 'fex_mincoadd'
      dictionary 'fex_cth'
      record /fex_mincoadd/ mincoadd
      record /fex_cth/ cth

      integer*4 num_sci_recs(4,4),num_groups(4,4)
      integer*4 num_sci_good(4,4,20),num_sci_fail(4,4,32)
      integer*4 num_sci_noise(4,4)
      real*4 tot_noise(4,4)
      integer*4 num_sci_deglitch(4,4),glitch_tot_iter(4,4)
      real*4 glitch_tot_signal(4,4)
      integer*4 num_coadds(4,6)

      character*72 output_filenames(4,6,7),report_filenames(4,4)
c
c     Local variables.
c
      integer*4 status,pos,channel,scan_mode,end_scan_mode,temp_scan_mode
      character*72 filename
      integer*4 file_len
      integer*4 tot_sci_recs,tot_groups,tot_sci_good(20),tot_sci_fail(32)
      integer*4 tot_sci_noise,tot_sci_deglitch,all_iter,iter
      real*4 all_noise,all_signal,noise,signal
      integer*4 tot_coadds,num_recs,num_fail,percent_fail(4,4,32)
      integer*4 num_good,percent_good,percent_tot_fail(32)
      logical*1 first

      fil_close_summary = %loc(fil_normal)
c
c     Write list of failure reasons.
c
      write(fut_report_lun,10,iostat=status)
10    format(/x,'Failure Reasons: ',/x,
     .       '1 - Data or Attitude Quality, Flagged Temperatures, ',
     .       'Unknown Gain, Zero Sweeps;',/x,
     .       '2 - Channel; 3 - MTM Speed; 4 - MTM Length; 5 - Fake-it; ',
     .       '6 - Science Mode;',/x,
     .       '7 - Adds-per-group; 8 - Xcal Position; 9 - Commanded Bias; ',
     .       '10 - Bolo Voltage;',/x,
     .       'Temperatures: 11 - Xcal; 12 - Ical; 13 - Skyhorn; 14 - Refhorn; ',
     .       '15 - Bolometer;',/x,
     .       'Spatial Gradients: 16 - Xcal S5-S6; 17 - Xcal Tip-Cone; ',
     .       '18 - Ical A-B;',/20x,
     .       '19 - Skyhorn A-B; 20 - Refhorn A-B;',/x,
     .       'Temporal Gradients: 21 - Xcal; 22 - Ical; 23 - Skyhorn; ',
     .       '24 - Refhorn;',/x,
     .       '25 - Temp Controller Integral and Proportional Gain Status Bits;',
     .       /x,'26 - Pixel Number; 27 - Pixel Definition, Skymap Index; ',
     .       '28 - Deglitch Failure;',/x,
     .       '29 - Low Noise; 30 - High Noise; 31 - Outlier Points; ',
     .       '32 - Too Few IFGs Left;')
      if (status .ne. 0) then
         fil_close_summary = %loc(fil_repwrite)
         inquire(fut_report_lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fil_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if
c
c     Write consistency check parameters from mincoadd and cth records.
c
      write(fut_report_lun,20,iostat=status) 
     .                            (mincoadd.min_ifg_coadd(pos), pos = 1,4),
     .                            (cth.min_ifg_shape_coadd(pos), pos = 1,4),
     .                            (cth.bolometer_voltage_tolerances(pos), 
     .                                                   pos = 1,4),
     .                            (cth.grt_tolerances(pos), pos = 1,4),
     .                            (cth.grt_tolerances(pos), pos = 7,10),
     .                            (cth.spatial_gradients(1,pos), pos = 1,3),
     .                            (cth.spatial_gradients(2,pos), pos = 1,3),
     .                            (cth.spatial_gradients(3,pos), pos = 1,3),
     .                            (cth.spatial_gradients(4,pos), pos = 1,3),
     .                            (cth.spatial_gradients(5,pos), pos = 1,3),
     .                            (cth.temporal_gradients(pos), pos = 1,4),
     .                            (cth.temporal_gradients(pos), pos = 11,14),
     .                            (cth.min_ifg_noise(pos), pos = 1,4),
     .                            (cth.max_ifg_noise(pos), pos = 1,4),
     .                            (cth.max_point_deviation(pos), pos = 1,4),
     .                            (cth.max_bad_points(pos), pos = 1,4)
20    format(/x,'Consistency Thresholds:',/13x,
     .       'Minimum Number IFGs Before Shape Check by Channel ',4(i2,'; '),
     .       /13x,
     .       'Minimum Number IFGs After Shape Check by Channel  ',4(i2,'; '),/x,
     .       'Tolerances: Bolometer Voltage by Channel ',4(f4.2,'; '),/13x,
     .       'Xcal ',f5.3,'; Ical ',f5.3,'; ',
     .       'Skyhorn ',f5.3,'; Refhorn ',f5.3,';',/13x,
     .       'Bolometer by Channel ',4(f5.3,'; '),/x,
     .       'Spatial Gradients:',/9x,
     .       'Xcal S5-S6    ',3(g15.9,3x),/9x,
     .       'Xcal Tip-Cone ',3(g15.9,3x),/9x,
     .       'Ical A-B      ',3(g15.9,3x),/9x,
     .       'Skyhorn A-B   ',3(g15.9,3x),/9x,
     .       'Refhorn A-B   ',3(g15.9,3x),/x,
     .       'Temporal Gradients:',/10x,
     .       'Xcal A ',i4,'; Ical A ',i4,'; ',
     .       'Skyhorn A ',i4,'; Refhorn A ',i4,';',/10x,  
     .       'Xcal B ',i4,'; Ical B ',i4,'; ',
     .       'Skyhorn B ',i4,'; Refhorn B ',i4,';',/x,  
     .       'Minimum IFG Noise Ratio by Channel ',4(f3.1,'; '),/x,
     .       'Maximum IFG Noise Ratio by Channel ',4(f3.1,'; '),/x,
     .       'Maximum Point Deviation by Channel ',4(f3.1,'; '),/x,
     .       'Maximum Bad Points by Channel ',4(i2,'; '))
      if (status .ne. 0) then
         fil_close_summary = %loc(fil_repwrite)
         inquire(fut_report_lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fil_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if
c
c     Loop over channels and scan modes.
c
      do channel = 1,4

         if ((channel .eq. 1) .or. (channel .eq. 3)) then
            end_scan_mode = 4
         else
            end_scan_mode = 6
         end if

         do scan_mode = 1,end_scan_mode
            if (scan_mode .le. 4) then
               temp_scan_mode = scan_mode
               num_recs = num_sci_recs(channel,scan_mode)
c
c     If any records were read for this channel and scan mode, calculate the
c     number of input records, number of input groups, number of good input
c     records, and number that failed.
c
               if (num_recs .gt. 0) then

                  tot_sci_recs = tot_sci_recs + num_sci_recs(channel,scan_mode)
                  tot_groups = tot_groups + num_groups(channel,scan_mode)
                  do pos = 1,20
                     tot_sci_good(pos) = tot_sci_good(pos) + 
     .                                   num_sci_good(channel,scan_mode,pos)
                  end do
                  do pos = 1,32
                     tot_sci_fail(pos) = tot_sci_fail(pos) + 
     .                                   num_sci_fail(channel,scan_mode,pos)
                  end do
c
c     Calculate summary information for noise and deglitching.
c
                  tot_sci_noise = tot_sci_noise + 
     .                            num_sci_noise(channel,scan_mode)
                  all_noise = all_noise + tot_noise(channel,scan_mode)
                  tot_sci_deglitch = tot_sci_deglitch + 
     .                               num_sci_deglitch(channel,scan_mode)
                  all_iter = all_iter + glitch_tot_iter(channel,scan_mode)
                  all_signal = all_signal + glitch_tot_signal(channel,scan_mode)
               end if
            else if (scan_mode .eq. 5) then
               temp_scan_mode = 2
               num_recs = num_sci_recs(channel,temp_scan_mode)
            else if (scan_mode .eq. 6) then
               temp_scan_mode = 4
               num_recs = num_sci_recs(channel,temp_scan_mode)
            end if

            tot_coadds = tot_coadds + num_coadds(channel,scan_mode)

            if (num_recs .gt. 0) then

               write(fut_report_lun,30,iostat=status) 
     .                                        fac_channel_ids(channel),
     .                                        fac_scan_mode_idsl(scan_mode)
30             format(/x,79('*')//x,'Summary for Channel ',a,'; ',
     .                'Scan Mode ',a,':')
               if (status .ne. 0) then
                  fil_close_summary = %loc(fil_repwrite)
                  inquire(fut_report_lun,name=filename)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fil_repwrite,%val(2),filename(:file_len),
     .                            %val(status))
                  return
               end if
c
c     Write summary statistics.
c
               write(fut_report_lun,40,iostat=status) num_recs,
     .                                 num_groups(channel,temp_scan_mode),
     .                                 num_coadds(channel,scan_mode)
40             format(/x,'Sci Records Read ',i6,'; ',
     .                'Coadd Groups Read ',i5,'; ',
     .                'Coadd Records Written ',i5,';')
               if (status .ne. 0) then
                  fil_close_summary = %loc(fil_repwrite)
                  inquire(fut_report_lun,name=filename)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fil_repwrite,%val(2),filename(:file_len),
     .                            %val(status))
                  return
               end if
c
c     Write pass and fail statistics.
c
               num_fail = 0
               do pos = 1,32
                  num_fail = num_fail + num_sci_fail(channel,temp_scan_mode,pos)
                  percent_fail(channel,temp_scan_mode,pos) = 
     .             nint(num_sci_fail(channel,temp_scan_mode,pos)*100.0/num_recs)
               end do
               num_good = num_recs - num_fail
               percent_good = nint(num_good * 100.0 / num_recs)
               if (num_sci_noise(channel,temp_scan_mode) .gt. 0) then
                  noise = tot_noise(channel,temp_scan_mode) / 
     .                    num_sci_noise(channel,temp_scan_mode)
               else
                  noise = 0.0
               end if

               write(fut_report_lun,50,iostat=status)
     .                  (num_sci_fail(channel,temp_scan_mode,pos), pos =  1,16),
     .                  (percent_fail(channel,temp_scan_mode,pos), pos =  1,16),
     .                  (num_sci_fail(channel,temp_scan_mode,pos), pos = 17,32),
     .                  (percent_fail(channel,temp_scan_mode,pos), pos = 17,32)
50    format(/x,'Number and Percent of Records Failing for Each Reason:',/,
     .       16(i5),/x,16(i3,'% '),/,16(i5),/x,16(i3,'% '))
               if (status .ne. 0) then
                  fil_close_summary = %loc(fil_repwrite)
                  inquire(fut_report_lun,name=filename)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fil_repwrite,%val(2),filename(:file_len),
     .                            %val(status))
                  return
               end if
    
               write(fut_report_lun,60,iostat=status) 
     .                  (num_sci_good(channel,temp_scan_mode,pos), pos =  1,10),
     .                  (num_sci_good(channel,temp_scan_mode,pos), pos = 11,20)
60             format(/x,'Number of Good Science Records Per Coadd ',
     .                '(1-5, 6-10, ..., 96-100):',/x,10(i5,' '),/x,10(i5,' '))
               if (status .ne. 0) then
                  fil_close_summary = %loc(fil_repwrite)
                  inquire(fut_report_lun,name=filename)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fil_repwrite,%val(2),filename(:file_len),
     .                            %val(status))
                  return
               end if

               write(fut_report_lun,70,iostat=status) num_good,percent_good,
     .                                                noise
70             format(/x,'Number and Percent of Good Records ',i6,'; ',
     .                i3,'%; Average Noise ',g9.3)
               if (status .ne. 0) then
                  fil_close_summary = %loc(fil_repwrite)
                  inquire(fut_report_lun,name=filename)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fil_repwrite,%val(2),filename(:file_len),
     .                            %val(status))
                  return
               end if
c
c     Write deglitcher summary statistics.
c
               if ((deglitch .eq. fac_present) .and.
     .             (num_sci_deglitch(channel,temp_scan_mode) .gt. 0)) then
                  iter = nint(float(glitch_tot_iter(channel,temp_scan_mode)) / 
     .                        num_sci_deglitch(channel,temp_scan_mode))
                  signal = glitch_tot_signal(channel,temp_scan_mode) / 
     .                     num_sci_deglitch(channel,temp_scan_mode)
                  write(fut_report_lun,80,iostat=status) iter,signal
80                format(/x,'Average Deglitcher Iterations ',i4,
     .                   '; Average Total Glitch Signal ',g10.3)
                  if (status .ne. 0) then
                     fil_close_summary = %loc(fil_repwrite)
                     inquire(fut_report_lun,name=filename)
                     call str$trim(filename,filename,file_len)
                     call lib$signal(fil_repwrite,%val(2),filename(:file_len),
     .                               %val(status))
                     return
                  end if
               end if
c
c     Write output filenames.
c
               first = .true.
               do pos = 1,6
                  if (output(pos) .eq. fac_present) then
                     filename = output_filenames(channel,scan_mode,pos)
                     call str$trim(filename,filename,file_len)
                     if (file_len .gt. 0) then
                        if (first) then
                           first = .false.
                           write(fut_report_lun,90,iostat=status)
90                         format(/x,'Output Filenames:')
                           if (status .ne. 0) then
                              fil_close_summary = %loc(fil_repwrite)
                              inquire(fut_report_lun,name=filename)
                              call str$trim(filename,filename,file_len)
                              call lib$signal(fil_repwrite,%val(2),
     .                                        filename(:file_len),%val(status))
                              return
                           end if
                        end if
                        write(fut_report_lun,100,iostat=status) 
     .                        filename(:file_len)
100                     format(x,a)
                        if (status .ne. 0) then
                           fil_close_summary = %loc(fil_repwrite)
                           inquire(fut_report_lun,name=filename)
                           call str$trim(filename,filename,file_len)
                           call lib$signal(fil_repwrite,%val(2),
     .                                     filename(:file_len),%val(status))
                           return
                        end if
                     end if
                  end if
               end do
c
c     Write covariance matrix names if they were calculated.
c
               if ((covar .eq. fac_present) .and. 
     .             (num_coadds(channel,scan_mode) .gt. 0)) then
                  filename = output_filenames(channel,scan_mode,7)
                  call str$trim(filename,filename,file_len)
                  if (file_len .gt. 0) then
                     write(fut_report_lun,110,iostat=status) filename(:file_len)
110                  format(/x,'Covariance Matrix: ',a)
                     if (status .ne. 0) then
                        fil_close_summary = %loc(fil_repwrite)
                        inquire(fut_report_lun,name=filename)
                        call str$trim(filename,filename,file_len)
                        call lib$signal(fil_repwrite,%val(2),
     .                                  filename(:file_len),%val(status))
                        return
                     end if
                  end if
               end if
c
c     Write report names if they were produced.
c
               if ((report .eq. fac_present) .and. (scan_mode .le. 4)) then
                  filename = report_filenames(channel,scan_mode)
                  call str$trim(filename,filename,file_len)
                  if (file_len .gt. 0) then
                     write(fut_report_lun,120,iostat=status) filename(:file_len)
120                  format(/x,'Report Filename: ',a)
                     if (status .ne. 0) then
                        fil_close_summary = %loc(fil_repwrite)
                        inquire(fut_report_lun,name=filename)
                        call str$trim(filename,filename,file_len)
                        call lib$signal(fil_repwrite,%val(2),
     .                                  filename(:file_len),%val(status))
                        return
                     end if
                  end if
               end if
            end if
         end do
      end do
c
c     Write overall summary statistics.
c
      write(fut_report_lun,130,iostat=status) 
130   format(/x,79('*')//x,'Summary Over All Channels and Scan Modes:')
      if (status .ne. 0) then
         fil_close_summary = %loc(fil_repwrite)
         inquire(fut_report_lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fil_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if

      write(fut_report_lun,40,iostat=status) tot_sci_recs,tot_groups,tot_coadds
      if (status .ne. 0) then
         fil_close_summary = %loc(fil_repwrite)
         inquire(fut_report_lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fil_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if

      num_fail = 0
      do pos = 1,32
         num_fail = num_fail + tot_sci_fail(pos)
         if (tot_sci_recs .gt. 0) then
            percent_tot_fail(pos) = nint(tot_sci_fail(pos)*100.0 / tot_sci_recs)
         else
            percent_tot_fail(pos) = 0
         end if
      end do

      num_good = tot_sci_recs - num_fail
      if (tot_sci_recs .gt. 0) then
         percent_good = nint(num_good * 100.0 / tot_sci_recs)
      else
         percent_good = 0
      end if

      if (tot_sci_noise .gt. 0) then
         noise = all_noise/tot_sci_noise
      else
         noise = 0.0
      end if

      write(fut_report_lun,50,iostat=status) (tot_sci_fail(pos), pos = 1,16),
     .                                       (percent_tot_fail(pos), pos=1,16),
     .                                       (tot_sci_fail(pos), pos = 17,32),
     .                                       (percent_tot_fail(pos), pos=17,32)
      if (status .ne. 0) then
         fil_close_summary = %loc(fil_repwrite)
         inquire(fut_report_lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fil_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if

      write(fut_report_lun,60,iostat=status) (tot_sci_good(pos), pos =  1,10),
     .                                       (tot_sci_good(pos), pos = 11,20)
      if (status .ne. 0) then
         fil_close_summary = %loc(fil_repwrite)
         inquire(fut_report_lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fil_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if

      write(fut_report_lun,70,iostat=status) num_good,percent_good,noise
      if (status .ne. 0) then
         fil_close_summary = %loc(fil_repwrite)
         inquire(fut_report_lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fil_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if
c
c     Write overall deglitcher statistics.
c
      if ((deglitch .eq. fac_present) .and. (tot_sci_deglitch .gt. 0)) then
         iter = nint(float(all_iter) / tot_sci_deglitch)
         signal = all_signal / tot_sci_deglitch
         write(fut_report_lun,80,iostat=status) iter,signal
      end if
      if (status .ne. 0) then
         fil_close_summary = %loc(fil_repwrite)
         inquire(fut_report_lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fil_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if

      return
      end
