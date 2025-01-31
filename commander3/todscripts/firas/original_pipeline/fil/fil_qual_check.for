      integer*4 function fil_qual_check(channel,scan_mode_input,scan_mode,
     .                                  qual_inst,qual_att,single,deglitch,
     .                                  mincoadd,grttrans,grtcoawt,cmdgain,
     .                                  samprate,nyquistl,num_recs,num_cgr_recs,
     .                                  sci_recs,eng_recs,num_sci_recs,
     .                                  first_sci_rec,last_sci_rec,num_groups,
     .                                  num_good,num_cgr_good,min_good,
     .                                  num_sci_fail,fake_it,sci_mode,
     .                                  adds_per_group,xcal_pos,sweeps,temps,
     .                                  con_recs,coa_rec)
c-------------------------------------------------------------------------------
c
c     Purpose: Check science and engineering records for data and attitude 
c              quality, and for consistency in various instrument parameters.
c
c     Author: S. Brodd, HSTX, 4/95
c
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode_input    i*4  Input scan mode, 1-6 = SS-FL; 0 if none.
c            qual_inst          i*4  Maximum acceptable data quality.
c            qual_att           i*4  Maximum acceptable attitude quality.
c            single             i*4  Process single ifgs.
c            deglitch           i*4  Perform deglitching.
c            mincoadd           rec  Minimum interferograms to coadd.
c            grttrans           rec  Transition temperatures.
c            grtcoawt           rec  Coadd temperature weights.
c            cmdgain            rec  Real values of commanded gains.
c            samprate           rec  Sampling rate.
c            nyquistl           rec  Nyquist frequencies.
c            num_recs           i*4  Number of input science and engineering
c                                    records, including neighbors.
c            num_cgr_recs       i*4  Number of input science and engineering
c                                    records, excluding neighbors.
c            sci_recs           rec(num_recs)  Science records.
c            eng_recs           rec(num_recs)  Engineering records.
c            num_sci_recs       i*4(4,4)  Number of input science and
c                                         engineering records by channel and 
c                                         scan mode.
c            num_groups         i*4(4,4)  Number of input coadd groups by
c                                         channel and scan mode.
c            num_sci_fail       i*4(4,4,32)  Number of science or engineering
c                                            records failing by channel and 
c                                            scan mode for each of 32 
c                                            consistency checks.
c
c     Output: scan_mode          i*4  Value of scan mode, 1-4 = SS-LF.
c             num_sci_recs       i*4(4,4)  Number of input science and
c                                          engineering records by channel and 
c                                          scan mode.
c             first_sci_rec      i*4  Pointer to earliest science record in
c                                     coadd group.
c             last_sci_rec       i*4  Pointer to latest science record in 
c                                     coadd group.
c             num_groups         i*4(4,4)  Number of input coadd groups by
c                                          channel and scan mode.
c             num_good           i*4  Current number of good interferograms
c                                     including neighbors.
c             num_cgr_good       i*4  Current number of good interferograms
c                                     excluding neighbors.
c             min_good           i*4  Minimum number of good interferograms
c                                     required to coadd.
c             num_sci_fail       i*4(4,4,32)  Number of science or engineering
c                                             records failing by channel and 
c                                             scan mode for each of 32 
c                                             consistency checks.
c             fake_it            i*4  Value of fake-it bit, 0-1.
c             sci_mode           i*4  Value of science mode, 0-4.
c             adds_per_group     i*4  Value of adds per group, 1-12.
c             xcal_pos           i*4  External calibrator position,
c                                     1=in=calibration data, 2=out=sky data.
c             sweeps             i*4  Value of sweeps, 1-16.
c             temps              r*4(10,num_recs)  Combined engineering
c                                                  temperatures.
c             con_recs           rec(num_recs)  Consistency check records.
c             coa_rec            rec  Coadd record.
c
c     Modifications:
c    
c-------------------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_params)'
c
c     Return statuses.
c 
      external fil_normal
c
c     Functions.
c
      logical*1 time_lt,time_gt
      integer*4 fut_ave_bin_times
      integer*4 fut_temp_list
      integer*4 fut_default_peak
c
c     Input parameters.
c
      integer*4 channel,scan_mode_input,qual_inst,qual_att
      integer*4 single,deglitch,num_recs,num_cgr_recs

      dictionary 'fex_mincoadd'
      dictionary 'fex_grttrans'
      dictionary 'fex_grtcoawt'
      dictionary 'fex_cmdgain'
      dictionary 'fex_samprate'
      dictionary 'fex_nyquistl'
      dictionary 'fdq_sdf'
      dictionary 'fdq_eng'

      record /fex_mincoadd/ mincoadd
      record /fex_grttrans/ grttrans
      record /fex_grtcoawt/ grtcoawt
      record /fex_cmdgain/ cmdgain
      record /fex_samprate/ samprate
      record /fex_nyquistl/ nyquistl
      record /fdq_sdf/ sci_recs(num_recs)
      record /fdq_eng/ eng_recs(num_recs)
c
c     Input/output parameters.
c
      integer*4 num_sci_recs(4,4),num_groups(4,4),num_sci_fail(4,4,32) 
c
c     Output parameters.
c
      integer*4 scan_mode,first_sci_rec,last_sci_rec
      integer*4 num_good,num_cgr_good,min_good
      integer*4 fake_it,sci_mode,adds_per_group,xcal_pos,sweeps
      real*4 temps(10,num_recs)

      dictionary 'fil_scc'
      dictionary 'fil_sky'

      record /fil_scc/ con_recs(num_recs)
      record /fil_sky/ coa_rec
c
c     Local variables.
c
      integer*4 first_time(2),last_time(2),rec
      integer*4 time_array(2,fac_max_num),time_weights(fac_max_num)
      integer*4 status,mtm_speed,mtm_length
      logical*1 singlifg
      integer*4 pos,rec_scan_mode,val,bit,con_most
      integer*4 con_array(0:16),index
      real*4 pts_per_sweep

      fil_qual_check = %loc(fil_normal)
c
c     Establish default first and last record times.
c
      call ct_gmt_to_binary(fac_jstop_default,first_time)
      call ct_gmt_to_binary(fac_jstart_default,last_time)
c
c     Put information into consistency check records so they correspond with
c     science records.
c
      do rec = 1,num_recs
         con_recs(rec).gmt = sci_recs(rec).ct_head.gmt
         con_recs(rec).time(1) = sci_recs(rec).ct_head.time(1)
         con_recs(rec).time(2) = sci_recs(rec).ct_head.time(2)
         con_recs(rec).pixel_no = sci_recs(rec).attitude.pixel_no
         con_recs(rec).earth_limb = sci_recs(rec).attitude.earth_limb
         con_recs(rec).earth_limb_azimuth = 
     .                 sci_recs(rec).attitude.earth_limb_azimuth
         con_recs(rec).sun_angle = sci_recs(rec).attitude.sun_angle
         con_recs(rec).moon_angle = sci_recs(rec).attitude.moon_angle
         con_recs(rec).moon_az_angle = sci_recs(rec).attitude.moon_az_angle
         con_recs(rec).moon_phase = sci_recs(rec).attitude.moon_phase
         con_recs(rec).sun_moon_dist = sci_recs(rec).attitude.sun_moon_dist
         con_recs(rec).cobe_moon_dist = sci_recs(rec).attitude.cobe_moon_dist
c
c     Determine earliest and latest science records, excluding neighbors.
c
         if (rec .le. num_cgr_recs) then
            if (time_lt(con_recs(rec).time,first_time)) then
               first_time(1) = con_recs(rec).time(1)
               first_time(2) = con_recs(rec).time(2)
               first_sci_rec = rec
            end if
            if (time_gt(con_recs(rec).time,last_time)) then
               last_time(1) = con_recs(rec).time(1)
               last_time(2) = con_recs(rec).time(2)
               last_sci_rec = rec
            end if
            time_array(1,rec) = con_recs(rec).time(1)
            time_array(2,rec) = con_recs(rec).time(2)
            time_weights(rec) = 1
         end if
      end do
c
c     Record times for single interferograms, otherwise determine an initial
c     coadd time in case the group is never coadded.
c
      if (single .eq. fac_present) then
         coa_rec.ct_head.gmt = sci_recs(1).ct_head.gmt
         coa_rec.ct_head.time(1) = sci_recs(1).ct_head.time(1)
         coa_rec.ct_head.time(2) = sci_recs(1).ct_head.time(2)
      else
         status = fut_ave_bin_times(time_array,num_cgr_recs,time_weights,
     .                              coa_rec.ct_head.time)
         call ct_binary_to_gmt(coa_rec.ct_head.time,coa_rec.ct_head.gmt)
      end if
c
c     Put times in the consistency check and coadd records.
c
      do rec = 1,num_recs
         con_recs(rec).coadd_gmt = coa_rec.ct_head.gmt
         con_recs(rec).coadd_time(1) = coa_rec.ct_head.time(1)
         con_recs(rec).coadd_time(2) = coa_rec.ct_head.time(2)
      end do

      coa_rec.coad_spec_head.first_gmt = sci_recs(first_sci_rec).ct_head.gmt
      coa_rec.coad_spec_head.first_time(1) = 
     .        sci_recs(first_sci_rec).ct_head.time(1)
      coa_rec.coad_spec_head.first_time(2) = 
     .        sci_recs(first_sci_rec).ct_head.time(2)
      do pos = 1,6
         coa_rec.coad_spec_head.first_space_time(pos) = 
     .           sci_recs(first_sci_rec).ct_head.space_time(pos)
      end do
      coa_rec.coad_spec_head.first_mjr_frm_no = 
     .        sci_recs(first_sci_rec).ct_head.mjr_frm_no

      coa_rec.coad_spec_head.last_gmt = sci_recs(last_sci_rec).ct_head.gmt
      coa_rec.coad_spec_head.last_time(1) = 
     .        sci_recs(last_sci_rec).ct_head.time(1)
      coa_rec.coad_spec_head.last_time(2) = 
     .        sci_recs(last_sci_rec).ct_head.time(2)
      do pos = 1,6
         coa_rec.coad_spec_head.last_space_time(pos) = 
     .           sci_recs(last_sci_rec).ct_head.space_time(pos)
      end do
      coa_rec.coad_spec_head.last_mjr_frm_no = 
     .        sci_recs(last_sci_rec).ct_head.mjr_frm_no

      num_good = num_recs
      num_cgr_good = num_cgr_recs
      coa_rec.coad_spec_head.num_ifgs = num_cgr_recs
c
c     Determine minimum number of good interferograms required to coadd.
c
      if (single .eq. fac_present) then
         min_good = 1
      else
         min_good = mincoadd.min_ifg_coadd(channel)
         if (num_cgr_recs .lt. min_good) then
            coa_rec.coad_spec_data.orphans = fac_present
         else
            coa_rec.coad_spec_data.orphans = fac_not_present
         end if 
      end if
c
c     Set defaults for instrument parameters.
c
      if (scan_mode_input .eq. 5) then
         scan_mode = 2
      else if (scan_mode_input .eq. 6) then
         scan_mode = 4
      else
         scan_mode = scan_mode_input
      end if

      mtm_speed = -1
      mtm_length = -1
      fake_it = -1
      sci_mode = -1
      adds_per_group = -1
      xcal_pos = -1
      sweeps = -1
c
c     Get combined engineering temperatures by calling fut_temp_list.  Check
c     each record for bad temperatures, instrument and attitude quality, gain
c     outside of allowed values, and zero sweeps.  Use the coadd record as a
c     temporary storage place for the temperatures.
c
      rec = 1
      singlifg = .true.
      do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .          (rec .le. num_recs))
         status = fut_temp_list(eng_recs(rec).en_analog,coa_rec.en_sigma,
     .                          grtcoawt,grttrans,0,singlifg,
     .                          coa_rec.coad_spec_data.temp,
     .                          coa_rec.coad_spec_data.temp_sigma)
c
c     If record fails these checks, mark the consistency check record and
c     decrement the number of good records remaining.
c
         if ((status .ne. 0) .or.
     .       (sci_recs(rec).dq_data.data_quality(110) .gt. qual_inst) .or.
     .       (sci_recs(rec).dq_data.data_quality(109) .gt. qual_att) .or.
     .       (sci_recs(rec).sci_head.gain .lt. fac_min_gain_set) .or.
     .       (sci_recs(rec).sci_head.gain .gt. fac_max_gain_set) .or.
     .       (sci_recs(rec).sci_head.sc_head11 .eq. 0)) then
            con_recs(rec).con_check = ibset(con_recs(rec).con_check,0)
            num_good = num_good - 1
            if (rec .le. num_cgr_recs) then
               num_cgr_good = num_cgr_good - 1
            end if
         else
            do pos = 1,10
               temps(pos,rec) = coa_rec.coad_spec_data.temp(pos)
            end do
         end if
         rec = rec + 1
      end do
c
c     Consistency check the channel value.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
         rec = 1
         do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .             (rec .le. num_recs))
            if (con_recs(rec).con_check .eq. 0) then
c
c     If record fails this check, mark the consistency check record and
c     decrement the number of good records remaining.
c
               if (sci_recs(rec).sci_head.chan_id .ne. channel) then
                  con_recs(rec).con_check = ibset(con_recs(rec).con_check,1)
                  num_good = num_good - 1
                  if (rec .le. num_cgr_recs) then
                     num_cgr_good = num_cgr_good - 1
                  end if
               end if
            end if
            rec = rec + 1
         end do
      end if
c
c     Consistency check the input scan mode, if it was specified.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
         if (scan_mode .ne. 0) then
            mtm_speed = mod(scan_mode-1,2)
            mtm_length = int((scan_mode-1)/2)
            rec = 1
            do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .                (rec .le. num_recs))
               if (con_recs(rec).con_check .eq. 0) then
                  rec_scan_mode = sci_recs(rec).sci_head.mtm_length*2 +
     .                            sci_recs(rec).sci_head.mtm_speed + 1
c
c     If record fail this check, mark the consistency check record and
c     decrement the number of good records remaining.
c
                  if (rec_scan_mode .ne. scan_mode) then
                     con_recs(rec).con_check = ibset(con_recs(rec).con_check,2)
                     num_good = num_good - 1
                     if (rec .le. num_cgr_recs) then
                        num_cgr_good = num_cgr_good - 1
                     end if
                  end if
               end if
               rec = rec + 1
            end do
         else
c
c     If no scan mode was specified on the command line, consistency check 
c     the plurality of the mtm speed; legal values are 0 and 1.
c 
            if (single .eq. fac_present) then
               mtm_speed = sci_recs(1).sci_head.mtm_speed
               mtm_length = sci_recs(1).sci_head.mtm_length
            else
               do val = 0,1
                  con_array(val) = 0
               end do
               do rec = 1,num_recs
                  if (con_recs(rec).con_check .eq. 0) then
                     do val = 0,1
                        if (sci_recs(rec).sci_head.mtm_speed .eq. val) then
                           con_array(val) = con_array(val) + 1
                        end if
                     end do
                  end if
               end do
               if (con_array(0) .ge. con_array(1)) then
                  mtm_speed = 0
               else
                  mtm_speed = 1
               end if

               rec = 1
               do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)
     .                   .and. (rec .le. num_recs))
                  if (con_recs(rec).con_check .eq. 0) then
c
c     If record fail this check, mark the consistency check record and
c     decrement the number of good records remaining.
c
                     if (sci_recs(rec).sci_head.mtm_speed .ne. mtm_speed) then
                        con_recs(rec).con_check = 
     .                           ibset(con_recs(rec).con_check,2)
                        num_good = num_good - 1
                        if (rec .le. num_cgr_recs) then
                           num_cgr_good = num_cgr_good - 1
                        end if
                     end if                                                    
                  end if
                  rec = rec + 1
               end do
c
c     Consistency check the plurality of the mtm length; legal values 
c     are 0 and 1.
c 
               if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
                  do val = 0,1
                     con_array(val) = 0
                  end do
    
                  do rec = 1,num_recs
                     if (con_recs(rec).con_check .eq. 0) then
                        do val = 0,1
                           if (sci_recs(rec).sci_head.mtm_length .eq. val) then
                              con_array(val) = con_array(val) + 1
                           end if
                        end do
                     end if
                  end do
                  if (con_array(0) .ge. con_array(1)) then
                     mtm_length = 0
                  else
                     mtm_length = 1
                  end if

                  rec = 1
                  do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)
     .                      .and. (rec .le. num_recs))
                     if (con_recs(rec).con_check .eq. 0) then
c
c     If record fail this check, mark the consistency check record and
c     decrement the number of good records remaining.
c
                        if (sci_recs(rec).sci_head.mtm_length 
     .                      .ne. mtm_length) then
                           con_recs(rec).con_check = 
     .                              ibset(con_recs(rec).con_check,3)
                           num_good = num_good - 1
                           if (rec .le. num_cgr_recs) then
                              num_cgr_good = num_cgr_good - 1
                           end if
                        end if
                     end if
                     rec = rec + 1
                  end do
               end if
            end if
         end if
      end if
c
c     Determine scan mode if none was specified on the command line.
c
      if (scan_mode .eq. 0) then
         if ((mtm_speed .ge. 0) .and. (mtm_length .ge.0)) then
            scan_mode = mtm_length*2 + mtm_speed + 1
         else
            scan_mode = sci_recs(1).sci_head.mtm_length*2 + 
     .                  sci_recs(1).sci_head.mtm_speed + 1
         end if
      end if
c
c     Update statistics arrays for reports.
c
      num_sci_recs(channel,scan_mode) = num_sci_recs(channel,scan_mode) + 
     .                                  num_cgr_recs
      num_groups(channel,scan_mode) = num_groups(channel,scan_mode) + 1

      do rec = 1,num_cgr_recs
         do bit = 0,3
            if (btest(con_recs(rec).con_check,bit)) then
               num_sci_fail(channel,scan_mode,bit+1) =
     .                      num_sci_fail(channel,scan_mode,bit+1) + 1
            end if
         end do
      end do
c
c     Consistency check the plurality of the fake-it bit; legal values
c     are 0 and 1.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
         if (single .eq. fac_present) then
            fake_it = sci_recs(1).dq_data.fake
         else
            do val = 0,1
               con_array(val) = 0
            end do
            
            do rec = 1,num_recs
               if (con_recs(rec).con_check .eq. 0) then
                  do val = 0,1
                     if (sci_recs(rec).dq_data.fake .eq. val) then
                        con_array(val) = con_array(val) + 1
                     end if
                  end do
               end if
            end do

            if (con_array(0) .ge. con_array(1)) then
               fake_it = 0
            else
               fake_it = 1
            end if

            rec = 1
            do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .                (rec .le. num_recs))
               if (con_recs(rec).con_check .eq. 0) then
c
c     If record fail this check, mark the consistency check record and
c     decrement the number of good records remaining.
c
                  if (sci_recs(rec).dq_data.fake .ne. fake_it) then
                     con_recs(rec).con_check = ibset(con_recs(rec).con_check,4)
                     num_good = num_good - 1
                     if (rec .le. num_cgr_recs) then
                        num_cgr_good = num_cgr_good - 1
                        num_sci_fail(channel,scan_mode,5) =
     .                               num_sci_fail(channel,scan_mode,5) + 1 
                     end if
                  end if
               end if
               rec = rec + 1
            end do
         end if
      end if
c
c     Consistency check the plurality of the science mode; legal values 
c     are 0 through 4.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
         if (single .eq. fac_present) then
            sci_mode = sci_recs(1).sci_head.sc_head1a
         else

            do val = 0,4
               con_array(val) = 0
            end do

            do rec = 1,num_recs
               if (con_recs(rec).con_check .eq. 0) then
                  do val = 0,4
                     if (sci_recs(rec).sci_head.sc_head1a .eq. val) then
                        con_array(val) = con_array(val) + 1
                     end if
                  end do
               end if
            end do
            con_most = -1
            do val = 0,4
               if (con_array(val) .gt. con_most) then
                  con_most = con_array(val)
                  sci_mode = val
               end if
            end do
            rec = 1
            do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .                (rec .le. num_recs))
               if (con_recs(rec).con_check .eq. 0) then
c
c     If record fail this check, mark the consistency check record and
c     decrement the number of good records remaining.
c
                  if (sci_recs(rec).sci_head.sc_head1a .ne. sci_mode) then
                     con_recs(rec).con_check = ibset(con_recs(rec).con_check,5)
                     num_good = num_good - 1
                     if (rec .le. num_cgr_recs) then
                        num_cgr_good = num_cgr_good - 1
                        num_sci_fail(channel,scan_mode,6) = 
     .                               num_sci_fail(channel,scan_mode,6) + 1 
                     end if
                  end if
               end if
               rec = rec + 1
            end do
         end if
      end if
c
c     Consistency check the pluality of the adds per group value; legal
c     values are 1 through 12.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
         if (single .eq. fac_present) then
            adds_per_group = sci_recs(1).sci_head.sc_head9
         else
            do val = 1,12
               con_array(val) = 0
            end do
            rec = 1
            do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .                (rec .le. num_recs))
               if (con_recs(rec).con_check .eq. 0) then
                  val = sci_recs(rec).sci_head.sc_head9
c
c     Fail the record if the adds per group value is outside the legal limits,
c     or if the deglitcher will be used or this is fake-it data and the adds
c     per group value is other than 1,2,3,8, or 12.  The deglitcher glitch
c     profiles are defined only for these values of adds per group, and the
c     Nyquist frequency for fake-it data is defined only for these values.
c
                  if ((val .lt. fac_min_adds) .or. (val .gt. fac_max_adds) .or.
     .                (((deglitch .eq. fac_present) .or. (fake_it .eq. 1)) .and.
     .                (((val .ge. 4) .and. (val .le. 7)) .or.
     .                ((val .ge. 9) .and. (val .le. 11))))) then
                     con_recs(rec).con_check = ibset(con_recs(rec).con_check,6)
                     num_good = num_good - 1
                     if (rec .le. num_cgr_recs) then
                        num_cgr_good = num_cgr_good - 1
                        num_sci_fail(channel,scan_mode,7) = 
     .                               num_sci_fail(channel,scan_mode,7) + 1 
                     end if
                  else
                     con_array(val) = con_array(val) + 1
                  end if
               end if
               rec = rec + 1
            end do
            if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
               con_most = -1
               do val = 1,12
                  if (con_array(val) .gt. con_most) then
                      con_most = con_array(val)
                      adds_per_group = val
                  end if
               end do
               rec = 1
               do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) 
     .                   .and. (rec .le. num_recs))
                  if (con_recs(rec).con_check .eq. 0) then
c
c     If record fail this check, mark the consistency check record and
c     decrement the number of good records remaining.
c
                     if (sci_recs(rec).sci_head.sc_head9 
     .                   .ne. adds_per_group) then
                        con_recs(rec).con_check = 
     .                                ibset(con_recs(rec).con_check,6)
                        num_good = num_good - 1
                        if (rec .le. num_cgr_recs) then
                           num_cgr_good = num_cgr_good - 1
                           num_sci_fail(channel,scan_mode,7) = 
     .                                  num_sci_fail(channel,scan_mode,7) + 1 
                        end if
                     end if
                  end if
                  rec = rec + 1
               end do
            end if
         end if
      end if
c
c     Determine the majority external calibrator position for processing
c     raw data.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
         if (single .eq. fac_present) then
            xcal_pos = sci_recs(1).dq_data.xcal_pos
         else
            do val = 1,2
               con_array(val) = 0
            end do

            do rec = 1,num_recs
               con_recs(rec).xcal_pos = sci_recs(rec).dq_data.xcal_pos
               if (con_recs(rec).con_check .eq. 0) then
                  do val = 1,2
                     if (con_recs(rec).xcal_pos .eq. val) then
                        con_array(val) = con_array(val) + 1
                     end if
                  end do
               end if
            end do

            if (con_array(1) .ge. con_array(2)) then
               xcal_pos = 1
            else
               xcal_pos = 2
            end if
         end if
      end if
c
c     Consistency check the sweeps value; legal values are 1 through 16.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
         if (single .eq. fac_present) then
            sweeps = sci_recs(1).sci_head.sc_head11
         else

            do val = 1,16
               con_array(val) = 0
            end do

            do rec = 1,num_recs
               if (con_recs(rec).con_check .eq. 0) then
                  do val = 1,16
                     if (sci_recs(rec).sci_head.sc_head11 .eq. val) then
                        con_array(val) = con_array(val) + 1
                     end if
                  end do
               end if
            end do
            con_most = -1
            do val = 1,16
               if (con_array(val) .gt. con_most) then
                  con_most = con_array(val)
                  sweeps = val
               end if
            end do
            rec = 1
            do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .                (rec .le. num_recs))
               if (con_recs(rec).con_check .eq. 0) then
c
c     If record fail this check, mark the consistency check record and
c     decrement the number of good records remaining.
c
                  if (sci_recs(rec).sci_head.sc_head11 .ne. sweeps) then
                     con_recs(rec).con_check = ibset(con_recs(rec).con_check,6)
                     num_good = num_good - 1
                     if (rec .le. num_cgr_recs) then
                        num_cgr_good = num_cgr_good - 1
                        num_sci_fail(channel,scan_mode,7) = 
     .                               num_sci_fail(channel,scan_mode,7) + 1 
                     end if
                  end if
               end if
               rec = rec + 1
            end do
         end if
      end if
c
c     If too few good interferograms remain to coadd, then mark the rest of the
c     group as failing.
c
      if ((num_good .lt. min_good) .or. (num_cgr_good .lt. 1)) then
         num_good = 0
         num_cgr_good = 0
         do rec = 1,num_recs
            if (con_recs(rec).con_check .eq. 0) then
               con_recs(rec).con_check = ibset(con_recs(rec).con_check,31)
               if (rec .le. num_cgr_recs) then
                  num_sci_fail(channel,scan_mode,32) = 
     .                         num_sci_fail(channel,scan_mode,32) + 1 
               end if
            end if
         end do
      else
c
c     Find the real-valued gain for each record and put it in the consistency
c     check record.  Also calculate the glitch rate and the threshold bit noise
c     and store them in the consistency check record.
c
         do rec = 1,num_recs
            if (con_recs(rec).con_check .eq. 0) then
               con_recs(rec).gain =
     .             cmdgain.chan(channel).cmdgain(sci_recs(rec).sci_head.gain+1)
               pts_per_sweep = sci_recs(rec).sci_head.sc_head10
               if (pts_per_sweep .eq. 0.0) then
                  pts_per_sweep = 1.0
               end if
               con_recs(rec).glitch_rate = (samprate.sampling_rate * 
     .                       sci_recs(rec).sci_head.sc_head21) /
     .                       (pts_per_sweep * sweeps * adds_per_group)
               con_recs(rec).noise = 1.0 / (con_recs(rec).gain * 
     .                                      sqrt(12.0 * adds_per_group))
            end if
         end do
c
c     Get the interferogram peak position from fut_default_peak; set it to a
c     flag value if fut_default_peak would return an error condition.  In the
c     function call, the zero indicates that the interferogram is not
c     linearized.
c
         if ((adds_per_group .eq. 1) .or. ((adds_per_group .eq. 2) .and.
     .       (mtm_speed .eq. 0))) then
            coa_rec.coad_spec_data.peak_pos = -9999
         else
            status = fut_default_peak(mtm_speed,mtm_length,channel,
     .                                adds_per_group,sci_mode,0,
     .                                coa_rec.coad_spec_data.peak_pos)
         end if
c
c     Assign the Nyquist frequencies.
c
         if (fake_it .eq. 1) then
            if (adds_per_group .eq. 1) then
               coa_rec.coad_spec_data.nyquist_hertz = nyquistl.hz(11)
            else if (adds_per_group .eq. 2) then
               coa_rec.coad_spec_data.nyquist_hertz = nyquistl.hz(12)
            else if (adds_per_group .eq. 3) then
               coa_rec.coad_spec_data.nyquist_hertz = nyquistl.hz(13)
            else if (adds_per_group .eq. 8) then
               coa_rec.coad_spec_data.nyquist_hertz = nyquistl.hz(14)
            else if (adds_per_group .eq. 12) then
               coa_rec.coad_spec_data.nyquist_hertz = nyquistl.hz(15)
            end if
         else
            index = 4*(mod(channel-1,2)) + scan_mode
            coa_rec.coad_spec_data.nyquist_icm = nyquistl.icm(index)
            coa_rec.coad_spec_data.nyquist_hertz = nyquistl.hz(index)
         end if
      end if
c
c     Put determined instrument parameters into the consistency check records
c     and the coadd record.
c
      do rec = 1,num_recs
         con_recs(rec).chan_id = channel
         con_recs(rec).mtm_speed = mtm_speed
         con_recs(rec).mtm_length = mtm_length
         con_recs(rec).fakeit = fake_it
         con_recs(rec).sci_mode = sci_mode
         con_recs(rec).adds_per_group = adds_per_group
         con_recs(rec).sweeps = sweeps
      end do

      coa_rec.coad_spec_data.chan_id = channel
      coa_rec.coad_spec_data.mtm_speed = mtm_speed
      coa_rec.coad_spec_data.mtm_length = mtm_length
      coa_rec.coad_spec_data.fakeit = fake_it
      coa_rec.coad_spec_data.sci_mode = sci_mode
      coa_rec.coad_spec_data.adds_per_group = adds_per_group
      coa_rec.coad_spec_data.sweeps = sweeps

      return
      end
