      integer*4 function fil_open_report(channel,scan_mode,report_name,account,
     .                                   current_gmt,mincoadd,cth,lun,
     .                                   report_filenames)
c----------------------------------------------------------------------------
c
c     Purpose: Open and initialize a report file for a channel and scan mode.
c
c     Author: S. Brodd, HSTX, 4/95
c
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode          i*4  Value of scan mode, 1-4 = SS-LF.
c            report_name        ch*72  Consistency and deglitch report name on
c                                      command line or from summary report name.
c            account            ch*12  Account owner name.
c            current_gmt        ch*14  Current time in Julian format.
c            mincoadd           rec  Minimum interferograms to coadd.
c            cth                rec  Consistency check parameters.
c            report_filenames   ch*72(4,4)  Report filenames by channel and
c                                           scan mode.
c
c     Output: lun                i*4  Logical unit for output.
c             report_filenames   ch*72(4,4)  Report filenames by channel and
c                                            scan mode.
c
c     Modifications:
c
c----------------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_params)'
c
c     Return statuses.
c
      external fil_normal
      external fil_repopen
      external fil_repwrite
c
c     Functions.
c
      integer*4 fut_get_lun
      integer*4 cut_display_banner
c
c     Input parameters.
c
      integer*4 channel,scan_mode
      character*72 report_name
      character*12 account
      character*14 current_gmt

      dictionary 'fex_mincoadd'
      dictionary 'fex_cth'
      record /fex_mincoadd/ mincoadd
      record /fex_cth/ cth
c
c     Input/output parameters.
c
      character*72 report_filenames(4,4)
c
c     Output parameters.
c
      integer*4 lun
c    
c     Local variables.
c
      integer*4 period,status,pos
      integer*4 report_name_len
      character*72 filename
      integer*4 file_len

      fil_open_report = %loc(fil_normal)
c
c     Assign the report filename for this channel and scan mode based on the
c     report name specified on the command line or the default summary report
c     name if none was specified.
c
      call str$trim(report_name,report_name,report_name_len)
      period = index(report_name,'.')

      filename = report_name(:period-1)//'_'//fac_channel_ids(channel)//
     .           fac_scan_mode_idsl(scan_mode)//
     .           report_name(period:report_name_len)

      call str$trim(filename,filename,file_len)
c
c     Save the report filename for the summary report.
c
      report_filenames(channel,scan_mode) = filename
c
c     Open report file.
c
      status = fut_get_lun(lun)
      if (.not. status) then
         fil_open_report = status
         return
      end if

      open(unit=lun,file=filename,status='new',iostat=status)
      if (status .ne. 0) then
         fil_open_report = %loc(fil_repopen)
         call lib$signal(fil_repopen,%val(2),filename(:file_len),%val(status))
         return
      end if
c
c     Write banner, report name, account, and current time.
c
      status = cut_display_banner(lun,80,
     .'FIRAS Facility FIL_Interferogram_Long Consistency Check Deglitch Report'
     .)

      write(lun,10,iostat=status) 'Report File: ',filename(:file_len)
10    format(/x,2a)
      if (status .ne. 0) then
         fil_open_report = %loc(fil_repwrite)
         call lib$signal(fil_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if

      write(lun,20,iostat=status) 'Account: ',account,
     .                            '     Time: ',current_gmt(1:9)
20    format(/x,4a)
      if (status .ne. 0) then
         fil_open_report = %loc(fil_repwrite)
         call lib$signal(fil_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if
c
c     Write consistency check failure reasons.
c
      write(lun,30,iostat=status)
30    format(/x,'Failure Reasons: ',/x,
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
         fil_open_report = %loc(fil_repwrite)
         call lib$signal(fil_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if
c
c     Write consistency check parameters from the mincoadd and cth files.
c
      write(lun,40,iostat=status) mincoadd.min_ifg_coadd(channel),
     .                            cth.min_ifg_shape_coadd(channel),
     .                            cth.bolometer_voltage_tolerances(channel),
     .                            (cth.grt_tolerances(pos), pos = 1,4),
     .                            cth.grt_tolerances(channel+6),
     .                            (cth.spatial_gradients(1,pos), pos = 1,3),
     .                            (cth.spatial_gradients(2,pos), pos = 1,3),
     .                            (cth.spatial_gradients(3,pos), pos = 1,3),
     .                            (cth.spatial_gradients(4,pos), pos = 1,3),
     .                            (cth.spatial_gradients(5,pos), pos = 1,3),
     .                            (cth.temporal_gradients(pos), pos = 1,4),
     .                            (cth.temporal_gradients(pos), pos = 11,14),
     .                            cth.min_ifg_noise(channel),
     .                            cth.max_ifg_noise(channel),
     .                            cth.max_point_deviation(channel),
     .                            cth.max_bad_points(channel)
40    format(/x,'Consistency Thresholds:',/13x,
     .       'Minimum Number IFGs Before Shape Check ',i2,'; After ',i2,';',/x,
     .       'Tolerances: Bolometer Voltage ',f4.2,'; Xcal ',f5.3,'; ',
     .       'Ical ',f5.3,';',/13x,
     .       'Skyhorn ',f5.3,'; Refhorn ',f5.3,'; Bolometer ',f5.3,';',/x,
     .       'Spatial Gradients:',/9x,
     .       'Xcal S5-S6    ',3(g15.9,3x),/9x,
     .       'Xcal Tip-Cone ',3(g15.9,3x),/9x,
     .       'Ical A-B      ',3(g15.9,3x),/9x,
     .       'Skyhorn A-B   ',3(g15.9,3x),/9x,
     .       'Refhorn A-B   ',3(g15.9,3x),/x,
     .       'Temporal Gradients:',/10x,
     .       'Xcal A ',i4,'; Ical A ',i4,'; Skyhorn A ',i4,'; ',
     .       'Refhorn A ',i4,';',/10x,  
     .       'Xcal B ',i4,'; Ical B ',i4,'; Skyhorn B ',i4,'; ',
     .       'Refhorn B ',i4,';',/x,  
     .       'Minimum IFG Noise Ratio ',f3.1,'; ',
     .       'Maximum IFG Noise Ratio ',f3.1,';',/x,
     .       'Maximum Point Deviation ',f3.1,'; Maximum Bad Points ',i2,';')
      if (status .ne. 0) then
         fil_open_report = %loc(fil_repwrite)
         call lib$signal(fil_repwrite,%val(2),filename(:file_len),%val(status))
         return
      end if

      return
      end
