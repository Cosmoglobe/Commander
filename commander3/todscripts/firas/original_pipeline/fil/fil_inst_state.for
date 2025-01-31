      integer*4 function fil_inst_state(input,channel,scan_mode,pixel_no,
     .                                  mincoadd,grttrans,grtrawwt,cth,
     .                                  num_recs,num_cgr_recs,sci_recs,
     .                                  eng_recs,num_good,num_cgr_good,
     .                                  min_good,num_sci_fail,xcal_pos,temps,
     .                                  con_recs)
c-------------------------------------------------------------------------------
c
c     Purpose: Check science and engineering values for consistency in 
c              instrument state.
c
c     Authors: J. Durachta, 4/87
c              R. Kummerer, STX, 1/89
c              S. Brodd, HSTX, 4/95
c
c     Input: input              i*4  Type of input: sky, calibration, or raw.
c            channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode          i*4  Value of scan mode, 1-4=SS-LF.
c            pixel_no           i*4  Consistency check the pixel number.
c            mincoadd           rec  Minimum interferograms to coadd.
c            grttrans           rec  Transition temperatures.
c            grtrawwt           rec  Raw temperature weights.
c            cth                rec  Consistency check parameters.
c            num_recs           i*4  Number of input science and engineering
c                                    records, including neighbors.
c            num_cgr_recs       i*4  Number of input science and engineering
c                                    records, excluding neighbors.
c            sci_recs           rec(num_recs)  Science records.
c            eng_recs           rec(num_recs)  Engineering records.
c            num_good           i*4  Current number of good interferograms,
c                                    including neighbors.
c            num_cgr_good       i*4  Current number of good interferograms,
c                                    excluding neighbors.
c            num_sci_fail       i*4(4,4,32)  Number of science or engineering
c                                            records failing by channel and 
c                                            scan mode for each of 32 
c                                            consistency checks.
c            temps              r*4(10,num_recs)  Combined engineering
c                                                 temperatures.
c            con_recs           rec(num_recs)  Consistency check records.
c
c
c     Output: num_good           i*4  Current number of good interferograms,
c                                     including neighbors.
c             num_cgr_good       i*4  Current number of good interferograms,
c                                     excluding neighbors.
c             min_good           i*4  Minimum number of good interferograms
c                                     required to coadd.
c             num_sci_fail       i*4(4,4,32)  Number of science or engineering
c                                             records failing by channel and 
c                                             scan mode for each of 32 
c                                             consistency checks.
c             xcal_pos           i*4  External calibrator position,
c                                     1=in=calibration data, 2=out=sky data.
c             con_recs           rec(num_recs)  Consistency check records.
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
      integer*4 fut_temp_list
      integer*4 fut_combine_hilo
c
c     Input parameters.
c
      integer*4 input,channel,scan_mode,pixel_no,num_recs,num_cgr_recs

      dictionary 'fex_mincoadd'
      dictionary 'fex_grttrans'
      dictionary 'fex_grtrawwt'
      dictionary 'fex_cth'

      record /fex_mincoadd/ mincoadd
      record /fex_grttrans/ grttrans
      record /fex_grtrawwt/ grtrawwt
      record /fex_cth/ cth

      dictionary 'fdq_sdf'
      dictionary 'fdq_eng'

      record /fdq_sdf/ sci_recs(num_recs)
      record /fdq_eng/ eng_recs(num_recs)

      real*4 temps(10,num_recs)
c
c     Input/output parameters.
c
      integer*4 num_good,num_cgr_good,num_sci_fail(4,4,32)

      dictionary 'fil_scc'
      record /fil_scc/ con_recs(num_recs)
c
c     Output parameters.
c
      integer*4 min_good,xcal_pos
c
c     Local variables.
c
      integer*4 val,con_array(-3600:6143),rec
      integer*4 con_most,con_val,rec_good
      real*4 sort(fac_max_num),sum
      integer*4 first,last,num_midav_pts,mid
      integer*4 temp_pos,con_pos,status
      logical*1 singlifg

      dictionary 'fut_enganlg'
      record /fut_enganlg/ loc_sigmas

      real*4 temps_a(10,fac_max_num),temps_b(10,fac_max_num)
      real*4 loc_temps(10),loc_temp_sigmas(10),weight

      integer*4 stat_pos,group_num,bit_num
      character*1 char_val

      fil_inst_state = %loc(fil_normal)
c
c     Establish the minimum number of good interferograms required to coadd.
c
      min_good = mincoadd.min_ifg_coadd(channel)
c
c     Check the plurality of the external calibrator position.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
         do val = 1,2
            con_array(val) = 0
         end do
         do rec = 1,num_recs
            if (con_recs(rec).con_check .eq. 0) then
               do val = 1,2
                  if (sci_recs(rec).dq_data.xcal_pos .eq. val) then
                     con_array(val) = con_array(val) + 1
                  end if
               end do
            end if
         end do
         con_most = -1
         do val = 1,2
            if (con_array(val) .gt. con_most) then
                con_most = con_array(val)
                xcal_pos = val
            end if
         end do
         rec = 1
         do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .             (rec .le. num_recs))
            if (con_recs(rec).con_check .eq. 0) then
c
c     If external calibrator position does not match the plurality, mark
c     the record as bad and decrement the number of good records.
c
               if (sci_recs(rec).dq_data.xcal_pos .ne. xcal_pos) then
                  con_recs(rec).con_check = ibset(con_recs(rec).con_check,7)
                  num_good = num_good - 1
                  if (rec .le. num_cgr_recs) then
                     num_cgr_good = num_cgr_good - 1
                     num_sci_fail(channel,scan_mode,8) = 
     .                            num_sci_fail(channel,scan_mode,8) + 1
                  end if
               end if
            end if
            rec = rec + 1
         end do
      end if
c
c     Check the plurality of the commanded bias.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
         do val = 0,255
            con_array(val) = 0
         end do
         do rec = 1,num_recs
            if (con_recs(rec).con_check .eq. 0) then
               val = zext(eng_recs(rec).en_stat.bol_cmd_bias(channel))
               con_array(val) = con_array(val) + 1
            end if
         end do
         con_val = -1
         con_most = -1
         do val = 0,255
            if (con_array(val) .gt. con_most) then
               con_most = con_array(val)
               con_val = val
            end if
         end do
         rec = 1
         do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .             (rec .le. num_recs))
            if (con_recs(rec).con_check .eq. 0) then
c
c     If the commanded bias does not match the plurality, mark the record as 
c     bad and decrement the number of good records.
c
               val = zext(eng_recs(rec).en_stat.bol_cmd_bias(channel))
               if (val .ne. con_val) then
                  con_recs(rec).con_check = ibset(con_recs(rec).con_check,8)
                  num_good = num_good - 1
                  if (rec .le. num_cgr_recs) then
                     num_cgr_good = num_cgr_good - 1
                     num_sci_fail(channel,scan_mode,9) = 
     .                            num_sci_fail(channel,scan_mode,9) + 1
                  end if
               end if
            end if
            rec = rec + 1
         end do
      end if
c
c     Check that bolometer voltages fall within a tolerated range around the
c     midaverage.  The midaverage is calculated without any neighbors.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
         if (num_cgr_good .eq. 1) then
            sum = eng_recs(1).en_analog.bol_volt(channel)
         else if (num_cgr_good .eq. 2) then
            sum = (eng_recs(1).en_analog.bol_volt(channel) +
     .             eng_recs(2).en_analog.bol_volt(channel)) / 2.0    
         else
            rec_good = 0
            do rec = 1,num_cgr_recs
               if (con_recs(rec).con_check .eq. 0) then
                  rec_good = rec_good + 1
                  sort(rec_good) = eng_recs(rec).en_analog.bol_volt(channel)
               end if
            end do
            call svrgn(num_cgr_good,sort,sort)
c
c     Calculate the midaverage.
c
            sum = 0.0
            first = nint(num_cgr_good/4.0)
            last = nint((3.0*num_cgr_good)/4.0)
            num_midav_pts = last - first

            do mid = first+1,last
               sum = sum + sort(mid)
            end do
            sum = sum/num_midav_pts
         end if
         
         rec = 1
         do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .             (rec .le. num_recs))
            if (con_recs(rec).con_check .eq. 0) then
c
c     If the bolometer voltage does not fall within the allowed range, mark
c     the record as bad and decrement the number of good records.
c
               if ((eng_recs(rec).en_analog.bol_volt(channel) .lt. 
     .             (sum*(1.0 - cth.bolometer_voltage_tolerances(channel)))) .or.
     .             (eng_recs(rec).en_analog.bol_volt(channel) .gt. 
     .             (sum*(1.0+ cth.bolometer_voltage_tolerances(channel))))) then
                  con_recs(rec).con_check = ibset(con_recs(rec).con_check,9)
                  num_good = num_good - 1
                  if (rec .le. num_cgr_recs) then
                     num_cgr_good = num_cgr_good - 1
                     num_sci_fail(channel,scan_mode,10) = 
     .                            num_sci_fail(channel,scan_mode,10) + 1
                  end if
               end if
            end if
            rec = rec + 1
         end do
      end if
c
c     Check that the combined temperatures for the external calibrator, the
c     internal calibrator, the skyhorn, the reference horn, and the bolometer
c     for the current channel fall within tolerated ranges of the midaverage.
c
      temp_pos = 1
      do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .          (temp_pos .le. 5))
         con_pos = temp_pos + 9
         if (temp_pos .eq. 5) then
            temp_pos = temp_pos + channel + 1
         end if
c
c     Check the external calibrator temperature only for calibration data.
c
         if ((temp_pos .gt. 1) .or. ((input .eq. fac_cal) .or. 
     .       ((input .eq. fac_raw) .and. (xcal_pos .eq. fac_xcalin)))) then
c
c     Calculate the midaverage, excluding any neighbors.
c
            if (num_cgr_good .eq. 1) then
               sum = temps(temp_pos,1)
            else if (num_cgr_good .eq. 2) then
               sum = (temps(temp_pos,1) + temps(temp_pos,2)) / 2.0
            else
               first = nint(num_cgr_good/4.0)
               last = nint((3.0*num_cgr_good)/4.0)
               num_midav_pts = last - first

               rec_good = 0

               do rec = 1,num_cgr_recs
                  if (con_recs(rec).con_check .eq. 0) then
                     rec_good = rec_good + 1
                     sort(rec_good) = temps(temp_pos,rec)
                  end if
               end do
               call svrgn(num_cgr_good,sort,sort)

               sum = 0.0
               do mid = first+1,last
                  sum = sum + sort(mid)
               end do
               sum = sum/num_midav_pts
            end if
   
            rec = 1

            do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .                (rec .le. num_recs))
               if (con_recs(rec).con_check .eq. 0) then
c   
c     If the temperature is outside the allowed range, mark the record as bad
c     and decrement the number of good interferograms.
c
                  if ((temps(temp_pos,rec) .lt. 
     .                (sum*(1.0 - cth.grt_tolerances(temp_pos)))) .or.
     .                (temps(temp_pos,rec) .gt. 
     .                (sum*(1.0 + cth.grt_tolerances(temp_pos))))) then
                     con_recs(rec).con_check = 
     .                             ibset(con_recs(rec).con_check,con_pos)
                     num_good = num_good - 1
                     if (rec .le. num_cgr_recs) then
                        num_cgr_good = num_cgr_good - 1
                        num_sci_fail(channel,scan_mode,con_pos+1) = 
     .                      num_sci_fail(channel,scan_mode,con_pos+1) + 1
                     end if
                  end if
               end if
               rec = rec + 1
            end do
         end if
         temp_pos = temp_pos + 1
      end do
c
c     Use fut_temp_list to get the uncombined a and b side temperatures in
c     order to check the spatial gradients of the temperatures.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
         singlifg = .true.
         do rec = 1,num_recs
            if (con_recs(rec).con_check .eq. 0) then
               status = fut_temp_list(eng_recs(rec).en_analog,loc_sigmas,
     .                                grtrawwt,grttrans,1,singlifg,loc_temps,
     .                                loc_temp_sigmas)
               do temp_pos = 1,10
                  temps_a(temp_pos,rec) = loc_temps(temp_pos)
               end do
               status = fut_temp_list(eng_recs(rec).en_analog,loc_sigmas,
     .                                grtrawwt,grttrans,2,singlifg,loc_temps,
     .                                loc_temp_sigmas)
               do temp_pos = 1,10
                  temps_b(temp_pos,rec) = loc_temps(temp_pos)
               end do
            end if
         end do
      end if
c
c     For calibration data only, check the spatial gradient of the external
c     calibrator s5 and s6 grts, and the spatial gradient of the external
c     calibrator tip and cone grts.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .    ((input .eq. fac_cal) .or. ((input .eq. fac_raw) .and. 
     .    (xcal_pos .eq. fac_xcalin)))) then
         rec = 1
         do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .             (rec .le. num_recs))
            if (con_recs(rec).con_check .eq. 0) then
c
c     Use fut_combine_hilo to get the s5 and s6 grts.
c
               loc_temps(1) = eng_recs(rec).en_analog.a_hi_xcal_cone
               loc_temps(2) = eng_recs(rec).en_analog.a_lo_xcal_cone
               status = fut_combine_hilo(loc_temps(1),loc_temps(2),
     .                                   grttrans.grt_a_trans_temp(15),
     .                                   grttrans.grt_a_trans_hwid(15),weight)
               temps_a(1,rec) = (1.0-weight)*loc_temps(2) + weight*loc_temps(1)
               loc_temps(1) = eng_recs(rec).en_analog.b_hi_xcal_cone
               loc_temps(2) = eng_recs(rec).en_analog.b_lo_xcal_cone
               status = fut_combine_hilo(loc_temps(1),loc_temps(2),
     .                                   grttrans.grt_b_trans_temp(15),
     .                                   grttrans.grt_b_trans_hwid(15),weight)
               temps_b(1,rec) = (1.0-weight)*loc_temps(2) + weight*loc_temps(1)
c
c     If the spatial gradient is too great, mark the record as bad and 
c     decrement the number of good records.
c
               if (abs(temps_b(1,rec) - temps_a(1,rec)) .gt.
     .             (cth.spatial_gradients(1,1) + cth.spatial_gradients(1,2) *
     .             temps_a(1,rec) + cth.spatial_gradients(1,3) * 
     .             temps_b(1,rec))) then
                  con_recs(rec).con_check = ibset(con_recs(rec).con_check,15)
                  num_good = num_good - 1
                  if (rec .le. num_cgr_recs) then
                     num_cgr_good = num_cgr_good - 1
                     num_sci_fail(channel,scan_mode,16) = 
     .                            num_sci_fail(channel,scan_mode,16) + 1
                  end if
               else
c
c     Use the average of s5 and s6 as the external calibrator cone temperature.
c
                  temps_b(1,rec) = (temps_a(1,rec) + temps_b(1,rec))/2.0
c
c     Use fut_combine_hilo to get the external calibrator tip temperature.
c
                  loc_temps(1) = eng_recs(rec).en_analog.a_hi_xcal_tip
                  loc_temps(2) = eng_recs(rec).en_analog.a_lo_xcal_tip
                  status = fut_combine_hilo(loc_temps(1),loc_temps(2),
     .                                      grttrans.grt_a_trans_temp(1),
     .                                      grttrans.grt_a_trans_hwid(1),weight)
                  temps_a(1,rec)=(1.0-weight)*loc_temps(2) + weight*loc_temps(1)
c
c     If the spatial gradient between the tip and cone is too great, mark the 
c     record as bad and decrement the number of good records.
c
                  if (abs(temps_b(1,rec) - temps_a(1,rec)) .gt.
     .                (cth.spatial_gradients(2,1) + cth.spatial_gradients(2,2) *
     .                temps_a(1,rec) + cth.spatial_gradients(2,3) * 
     .                temps_b(1,rec))) then
                     con_recs(rec).con_check = ibset(con_recs(rec).con_check,16)
                     num_good = num_good - 1
                     if (rec .le. num_cgr_recs) then
                        num_cgr_good = num_cgr_good - 1
                        num_sci_fail(channel,scan_mode,17) = 
     .                               num_sci_fail(channel,scan_mode,17) + 1
                     end if
                  end if
               end if
            end if
            rec = rec + 1
         end do
      end if
c
c     Check the spatial gradients of the internal calibrator, skyhorn, and
c     reference horn.
c
      temp_pos = 2
      do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .          (temp_pos .le. 4))
         rec = 1
         do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .             (rec .le. num_recs))
            if (con_recs(rec).con_check .eq. 0) then
c
c     If the spatial gradient is too great, mark the record as bad and 
c     decrement the number of good records.
c
               if (abs(temps_b(temp_pos,rec) - temps_a(temp_pos,rec)) .gt.
     .             (cth.spatial_gradients(temp_pos+1,1) + 
     .              cth.spatial_gradients(temp_pos+1,2)*temps_a(temp_pos,rec) + 
     .              cth.spatial_gradients(temp_pos+1,3)*temps_b(temp_pos,rec)))
     .             then
                  con_recs(rec).con_check = 
     .                          ibset(con_recs(rec).con_check,temp_pos+15)
                  num_good = num_good - 1
                  if (rec .le. num_cgr_recs) then
                     num_cgr_good = num_cgr_good - 1
                     num_sci_fail(channel,scan_mode,temp_pos+16) = 
     .                   num_sci_fail(channel,scan_mode,temp_pos+16) + 1
                  end if
               end if
            end if
            rec = rec + 1
         end do
         temp_pos = temp_pos + 1
      end do
c
c     Check the temporal gradient of the external calibrator for calibration
c     data only.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .    ((input .eq. fac_cal) .or. ((input .eq. fac_raw) .and. 
     .    (xcal_pos .eq. fac_xcalin)))) then
         rec = 1
         do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .             (rec .le. num_recs))
            if (con_recs(rec).con_check .eq. 0) then
c
c     If the temporal gradient is too great, mark the record as bad and 
c     decrement the number of good records.
c
               if ((abs(eng_recs(rec).en_tempdiff(1).xcal) .gt. 
     .              cth.temporal_gradients(1)) .or.
     .             (abs(eng_recs(rec).en_tempdiff(2).xcal) .gt.
     .              cth.temporal_gradients(11))) then
                  con_recs(rec).con_check = ibset(con_recs(rec).con_check,20)
                  num_good = num_good - 1
                  if (rec .le. num_cgr_recs) then
                     num_cgr_good = num_cgr_good - 1
                     num_sci_fail(channel,scan_mode,21) = 
     .                            num_sci_fail(channel,scan_mode,21) + 1
                  end if
               end if
            end if
            rec = rec + 1
         end do
      end if
c
c     Check the temporal gradient of the internal calibrator.
c
      rec = 1
      do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .          (rec .le. num_recs))
         if (con_recs(rec).con_check .eq. 0) then
c
c     If the temporal gradient is too great, mark the record as bad and 
c     decrement the number of good records.
c
            if ((abs(eng_recs(rec).en_tempdiff(1).ical) .gt. 
     .           cth.temporal_gradients(2)) .or.
     .          (abs(eng_recs(rec).en_tempdiff(2).ical) .gt.
     .           cth.temporal_gradients(12))) then
               con_recs(rec).con_check = ibset(con_recs(rec).con_check,21)
               num_good = num_good - 1
                  if (rec .le. num_cgr_recs) then
                     num_cgr_good = num_cgr_good - 1
                     num_sci_fail(channel,scan_mode,22) = 
     .                            num_sci_fail(channel,scan_mode,22) + 1
                  end if
            end if
         end if
         rec = rec + 1
      end do
c
c     Check the temporal gradient of the skyhorn.
c
      rec = 1
      do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .          (rec .le. num_recs))
         if (con_recs(rec).con_check .eq. 0) then
c
c     If the temporal gradient is too great, mark the record as bad and 
c     decrement the number of good records.
c
            if ((abs(eng_recs(rec).en_tempdiff(1).skyhorn) .gt. 
     .           cth.temporal_gradients(3)) .or.
     .          (abs(eng_recs(rec).en_tempdiff(2).skyhorn) .gt.
     .           cth.temporal_gradients(13))) then
               con_recs(rec).con_check = ibset(con_recs(rec).con_check,22)
               num_good = num_good - 1
                  if (rec .le. num_cgr_recs) then
                     num_cgr_good = num_cgr_good - 1
                     num_sci_fail(channel,scan_mode,23) = 
     .                            num_sci_fail(channel,scan_mode,23) + 1
                  end if
            end if
         end if
         rec = rec + 1
      end do
c
c     Check the temporal gradient of the reference horn.
c
      rec = 1
      do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .          (rec .le. num_recs))
         if (con_recs(rec).con_check .eq. 0) then
c
c     If the temporal gradient is too great, mark the record as bad and 
c     decrement the number of good records.
c
            if ((abs(eng_recs(rec).en_tempdiff(1).refhorn) .gt. 
     .           cth.temporal_gradients(4)) .or.
     .          (abs(eng_recs(rec).en_tempdiff(2).refhorn) .gt.
     .           cth.temporal_gradients(14))) then
               con_recs(rec).con_check = ibset(con_recs(rec).con_check,23)
               num_good = num_good - 1
                  if (rec .le. num_cgr_recs) then
                     num_cgr_good = num_cgr_good - 1
                     num_sci_fail(channel,scan_mode,24) = 
     .                            num_sci_fail(channel,scan_mode,24) + 1
                  end if
            end if
         end if
         rec = rec + 1
      end do
c
c     Check the plurality of the temperature controller integral and
c     proportional gain status bits.  These are found in the engineering
c     record as follows:
c          field      bits  controller  int/prop  side
c        stat_word_4   0:2   refhorn      prop      a
c        stat_word_4   3:5   refhorn      int       a
c        stat_word_4   6:8   ical         prop      a
c        stat_word_4   9:11  ical         int       a
c        stat_word_8   0:2   xcal         prop      a
c        stat_word_8   3:5   xcal         int       a
c        stat_word_8   6:8   skyhorn      prop      a
c        stat_word_8   9:11  skyhorn      int       a
c        stat_word_12  0:2   refhorn      prop      b
c        stat_word_12  3:5   refhorn      int       b
c        stat_word_12  6:8   ical         prop      b
c        stat_word_12  9:11  ical         int       b
c        stat_word_16  0:2   xcal         prop      b
c        stat_word_16  3:5   xcal         int       b
c        stat_word_16  6:8   skyhorn      prop      b
c        stat_word_16  9:11  skyhorn      int       b
c     The bits are checked in eight groups of six bits each.
c
      stat_pos = 1
      do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .          (stat_pos .le. 8))
         group_num = (((stat_pos-1)/2)+1)*4
         bit_num = mod(stat_pos-1,2)*6
c
c     Check the external calibrator bits only for calibration data.
c
         if ((mod(stat_pos,4) .ne. 3) .or. ((input .eq. fac_cal) .or. 
     .       ((input .eq. fac_raw) .and. (xcal_pos .eq. fac_xcalin)))) then
            do val = 0,63
               con_array(val) = 0
            end do
            do rec = 1,num_recs
               if (con_recs(rec).con_check .eq. 0) then
                  val = ibits(eng_recs(rec).en_stat.group1(group_num),bit_num,6)
                  con_array(val) = con_array(val) + 1
               end if
            end do
            con_val = -1
            con_most = -1
            do val = 0,63
               if (con_array(val) .gt. con_most) then
                  con_most = con_array(val)
                  con_val = val
               end if
            end do
            rec = 1
            do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .                (rec .le. num_recs))
               if (con_recs(rec).con_check .eq. 0) then
c
c     If the records bits do not match the plurality, mark the record as
c     bad and decrement the number of good records.
c
                  if (ibits(eng_recs(rec).en_stat.group1(group_num),bit_num,6) 
     .                .ne. con_val) then
                     con_recs(rec).con_check = ibset(con_recs(rec).con_check,24)
                     num_good = num_good - 1
                     if (rec .le. num_cgr_recs) then
                        num_cgr_good = num_cgr_good - 1
                        num_sci_fail(channel,scan_mode,25) = 
     .                               num_sci_fail(channel,scan_mode,25) + 1
                     end if
                  end if
               end if
               rec = rec + 1
            end do
         end if
         stat_pos = stat_pos + 1
      end do
c
c     Check the plurality of the pixel number if processing sky data and the
c     pixel number consistency check is enabled.  Exclude neighbors since they
c     have different pixel numbers.
c
      if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .    ((input .eq. fac_sky) .or. ((input .eq. fac_raw) .and. 
     .    (xcal_pos .eq. fac_xcalout)))) then
         if (pixel_no .eq. fac_present) then
            do val = -1,6143
               con_array(val) = 0
            end do
            rec = 1
            do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .                (rec .le. num_recs))
               if (con_recs(rec).con_check .eq. 0) then
                  val = sci_recs(rec).attitude.pixel_no
c
c     If the records pixel number is outside the allowed range, mark the 
c     record as bad and decrement the number of good records.
c
                  if ((val .lt. -1) .or. (val .gt. 6143)) then
                     con_recs(rec).con_check = ibset(con_recs(rec).con_check,25)
                     num_good = num_good - 1
                     if (rec .le. num_cgr_recs) then
                        num_cgr_good = num_cgr_good - 1
                        num_sci_fail(channel,scan_mode,26) = 
     .                               num_sci_fail(channel,scan_mode,26) + 1
                     end if
                  else if (rec .le. num_cgr_recs) then
                     con_array(val) = con_array(val) + 1
                  end if
               end if
               rec = rec + 1
            end do
            if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
               con_val = -1
               con_most = -1
               do val = -1,6143
                  if (con_array(val) .gt. con_val) then
                      con_most = con_array(val)
                      con_val = val
                  end if
               end do
               rec = 1
               do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)
     .                   .and. (rec .le. num_cgr_recs))
                  if (con_recs(rec).con_check .eq. 0) then
c
c     If the records pixel number does not match the plurality, mark the
c     record as bad and decrement the number of good records.
c
                     if (sci_recs(rec).attitude.pixel_no .ne. con_val) then
                        con_recs(rec).con_check = 
     .                                ibset(con_recs(rec).con_check,25)
                        num_good = num_good - 1
                        num_cgr_good = num_cgr_good - 1
                        num_sci_fail(channel,scan_mode,26) = 
     .                               num_sci_fail(channel,scan_mode,26) + 1
                     end if
                  end if
                  rec = rec + 1
               end do
            end if
         end if
c          
c     Check the plurality of the pixel definition.
c          
         if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
            do val = 0,127
               con_array(val) = 0
            end do
            rec = 1
            do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .                (rec .le. num_recs))
               if (con_recs(rec).con_check .eq. 0) then
                  char_val = sci_recs(rec).attitude.pixel_definition
                  if (char_val .eq. ' ') then
                     char_val = 'q'
                  end if
c
c     If the records pixel definition does not match one of the allowed 
c     values, mark the record as bad and decrement the number of good records.
c
                  if ((char_val .ne. 'q') .and. (char_val .ne. 'O') .and.
     .                (char_val .ne. 'S') .and. (char_val .ne. 'E')) then
                     con_recs(rec).con_check = ibset(con_recs(rec).con_check,26)
                     num_good = num_good - 1
                     if (rec .le. num_cgr_recs) then
                        num_cgr_good = num_cgr_good - 1
                        num_sci_fail(channel,scan_mode,27) = 
     .                               num_sci_fail(channel,scan_mode,27) + 1
                     end if
                  else
                     val = ichar(char_val)
                     con_array(val) = con_array(val) + 1
                  end if
               end if
               rec = rec + 1
            end do
            if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
               con_val = -1
               con_most = -1
               do val = 0,127
                  if (con_array(val) .gt. con_val) then
                      con_most = con_array(val)
                      con_val = val
                  end if
               end do
               rec = 1
               do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)
     .                   .and. (rec .le. num_recs))
                  if (con_recs(rec).con_check .eq. 0) then
c
c     If the records pixel definition does not match the plurality, mark the
c     record as bad and decrement the number of good records.
c
                     if (ichar(sci_recs(rec).attitude.pixel_definition) .ne. 
     .                   con_val) then
                        con_recs(rec).con_check = 
     .                                ibset(con_recs(rec).con_check,26)
                        num_good = num_good - 1
                        if (rec .le. num_cgr_recs) then
                           num_cgr_good = num_cgr_good - 1
                           num_sci_fail(channel,scan_mode,27) = 
     .                                  num_sci_fail(channel,scan_mode,27) + 1
                        end if
                     end if
                  end if
                  rec = rec + 1
               end do
            end if
         end if
c
c     Check the plurality of the skymap index.
c
         if ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1)) then
            do val = -3600,3600
               con_array(val) = 0
            end do
            do rec = 1,num_recs
               if (con_recs(rec).con_check .eq. 0) then
                  val = sci_recs(rec).attitude.skymap_index
                  con_array(val) = con_array(val) + 1
               end if
            end do

            con_most = -1
            con_val = -1
            do val = -3600,3600
               if (con_array(val) .gt. con_most) then
                  con_most = con_array(val)
                  con_val = val
               end if
            end do
            rec = 1
            do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .                (rec .le. num_recs))
               if (con_recs(rec).con_check .eq. 0) then
c
c     If the records skymap index does not match the plurality, mark the
c     record as bad and decrement the number of good records.
c
                  if (sci_recs(rec).attitude.skymap_index .ne. con_val) then
                     con_recs(rec).con_check = ibset(con_recs(rec).con_check,26)
                     num_good = num_good - 1
                     if (rec .le. num_cgr_recs) then
                        num_cgr_good = num_cgr_good - 1
                        num_sci_fail(channel,scan_mode,27) = 
     .                               num_sci_fail(channel,scan_mode,27) + 1
                     end if
                  end if
               end if
               rec = rec + 1
            end do
         end if
      end if
c
c     If too few good interferograms remain to coadd, mark the remaining good
c     interferograms as bad.
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
      end if

      return
      end
