      integer*4 function fil_coadd(channel,scan_mode_input,scan_mode,covar,
     .                             single,grttrans,grtcoawt,nyquistl,cth,
     .                             gltchcor,apodl_lun,etfl_lun,num_cgr_recs,
     .                             num_cgr_good,num_sci_good,sci_recs,eng_recs,
     .                             fake_it,sci_mode,adds_per_group,temps,
     .                             con_recs,num_coadds,coa_rec,aifgs,template,
     .                             sec_template,sec_scan_mode,start_repeat,
     .                             stop_repeat,var_vecs,cov_recs,
     .                             first_covar_times,last_covar_times)
c-------------------------------------------------------------------------------
c
c     Purpose: Coadd science and engineering data and form variance vectors
c              from spectra; call fil_covar to calculate other statistical
c              quantities.
c
c     Author: S. Brodd, HSTX, 4/95
c
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode_input    i*4  Value of scan mode specified on command
c                                    line, 1-6 = SS-FL.
c            scan_mode          i*4  Value of scan mode, 1-4 = SS-LF.
c            covar              i*4  Covariance matrix to be written.
c            single             i*4  Process single ifgs.
c            grttrans           rec  Transition temperatures.
c            grtcoawt           rec  Coadd temperature weights.
c            nyquistl           rec  Nyquist frequencies.
c            cth                rec  Consistency check parameters.
c            gltchcor           rec  Glitch rate weighting parameters.
c            apodl_lun          i*4  Logical unit for apodization functions.
c            etfl_lun           i*4  Logical unit for electronics transfer
c                                    functions.
c            num_cgr_recs       i*4  Number of input science and engineering
c                                    records, excluding neighbors.
c            num_cgr_good       i*4  Current number of good interferograms,
c                                    excluding neighbors.
c            num_sci_good       i*4(4,4,20)  Number of good interferogrmas
c                                            contributing to an output coadd 
c                                            record by channel and scan mode; 
c                                            grouped as 1-5, ..., 96-100.
c            sci_recs           rec(num_cgr_recs)  Science records.
c            eng_recs           rec(num_cgr_recs)  Engineering records.
c            fake_it            i*4  Value of fake-it bit, 0-1.
c            sci_mode           i*4  Value of science mode, 0-4.
c            adds_per_group     i*4  Value of adds per group, 1-12.
c            temps              r*4(10,num_cgr_recs)  Combined engineering
c                                                 temperatures.
c            con_recs           rec(num_cgr_recs)  Consistency check records.
c            num_coadds         i*4(4,6)  Number of output coadd records by
c                                         channel and scan mode.
c            coa_rec            rec  Coadd record.
c            aifgs              r*4(512,num_cgr_recs)  Real-valued interferograms.
c            template           r*4(512)  Primary template.
c            sec_template       r*4(512)  Secondary template.
c            cov_recs           rec(1285)  Covariance matrix records.
c            first_covar_times  i*4(2,5)  Binary times of first records
c                                         contributing to covariance matrices.
c            last_covar_times   i*4(2,5)  Binary times of last records
c                                         contributing to covariance matrices.
c
c     Output: num_sci_good       i*4(4,4,20)  Number of good interferogrmas
c                                             contributing to an output coadd 
c                                             record by channel and scan mode; 
c                                             grouped as 1-5, ..., 96-100.
c             num_coadds         i*4(4,6)  Number of output coadd records by
c                                          channel and scan mode.
c             coa_rec            rec  Coadd record.
c             sec_scan_mode      i*4  Secondary scan mode, 5 = FS, 6 = FL.
c             start_repeat       i*4  Start flag for writing FS or FL coadds.
c             stop_repeat        i*4  Stop flag for writing FS or FL coadds.
c             var_vecs           r*4(3,361)  Variance vectors for FS or FL scan
c                                            modes.
c             cov_recs           rec(1285)  Covariance matrix records.
c             first_covar_times  i*4(2,5)  Binary times of first records
c                                          contributing to covariance matrices.
c             last_covar_times   i*4(2,5)  Binary times of last records
c                                          contributing to covariance matrices.
c
c     Modifications: S. Brodd, HSTX, 8/30/95, SPR 12246. Correctly fill in
c                    the fft length field in the coadd record in the case of
c                    a coadd with only one interferogram.
c                    S. Brodd, HSTX, 9/20/95, SPR 12256. Add array of glitch 
c                    rates to pass to fil_covar for glitch rate statistical 
c                    information.
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
      external fil_cfdirread
c
c     Functions.
c
      logical*1 time_lt,time_gt,time_ge
      integer*4 fut_apod_recnuml
      integer*4 fut_apod_rotl
      integer*4 fut_ave_bin_times
      integer*4 fut_get_recnum
      integer*4 fut_temp_list
      integer*4 fil_covar
      integer*4 fil_short
c
c     Input parameters.
c
      integer*4 channel,scan_mode_input,scan_mode
      integer*4 covar,single

      dictionary 'fex_grttrans'
      dictionary 'fex_grtcoawt'
      dictionary 'fex_nyquistl'
      dictionary 'fex_cth'
      dictionary 'fex_gltchcor'
      record /fex_grttrans/ grttrans
      record /fex_grtcoawt/ grtcoawt
      record /fex_nyquistl/ nyquistl
      record /fex_cth/ cth
      record /fex_gltchcor/ gltchcor

      integer*4 apodl_lun,etfl_lun,num_cgr_recs,num_cgr_good

      dictionary 'fdq_sdf'
      dictionary 'fdq_eng'
      dictionary 'fil_scc'
      record /fdq_sdf/ sci_recs(num_cgr_recs)
      record /fdq_eng/ eng_recs(num_cgr_recs)
      record /fil_scc/ con_recs(num_cgr_recs)

      integer*4 fake_it,sci_mode,adds_per_group
      real*4 temps(10,num_cgr_recs),aifgs(512,num_cgr_recs)
      real*4 template(512),sec_template(512)
c
c     Input/output parameters.
c
      integer*4 num_sci_good(4,4,20),num_coadds(4,6)

      dictionary 'fil_sky'
      dictionary 'fil_cov'
      record /fil_sky/ coa_rec
      record /fil_cov/ cov_recs(1285)

      integer*4 first_covar_times(2,5),last_covar_times(2,5)
c
c     Output parameters.
c
      integer*4 sec_scan_mode,start_repeat,stop_repeat
      real*4 var_vecs(3,361)
c
c     Local variables.
c
      integer*4 pos,rec_good,temp_scan_mode,rec,index
      real*8 glitch_intercept,glitch_slope,glwt(fac_max_num),glwt_tot
      integer*4 first_time(2),last_time(2)

      real*4 xcal_sum,glitch_rate_sum,glrt(fac_max_num)
      real*4 b_sum,c_sum,b_ave,c_ave,int_most

      integer*4 bins(fac_max_num),first_sci_rec,last_sci_rec
      integer*4 time_array(2,fac_max_num),time_weights(fac_max_num)
      integer*4 ifg_pos,status,int_val

      real*8 sum_equat(3),sum_terr(3),sum_elimb(3)
      real*8 sum_mnang(3),sum_gal(3),sum_ecl(3)

      real*4 terr_lat,terr_long,elimb,el_az
      real*4 moon_angle,moon_az
      real*4 gal_lat,gal_long,ecl_lat,ecl_long

      real*4 sum_sun_angle,sum_moon_phase,moon_phase(fac_max_num)
      real*4 sum_sm_dist,sum_cm_dist
      real*4 sum_altitude,sum_baryvel,sum_lparam
      real*4 sum_orb_phase,orb_phase(fac_max_num),sum_geovel
      real*4 sum_scan_angle,scan_angle(fac_max_num),sum_sc_rot_angle

      real*4 sum_att_sol(0:6),sum_pixel_no(4),sum_pix_def(4)
      real*4 sum_sky_idx(4),sum_exc_gal_lat

      real*8 norm_equat,norm_terr,norm_elimb
      real*8 norm_mnang,norm_gal,norm_ecl
      real*4 avg_terr(3),avg_elimb(3),avg_mnang(3)
      real*4 avg_gal(3),avg_ecl(3),conv_ecl(3)

      real*4 ra,dec,corr_angle,avg_angle

      integer*4 temp_pos,eng_pos
      real*4 sum,val,sum_eng(6)
      real*8 dsum,dval,rec_glwt
      logical*1 singlifg
      character*1 xcal_pos

      integer*4 repeat,fft_len,spec_len,peak
      integer*4 mtm_speed,temp_adds_per_group,etfl_rec,apodl_rec

      record /fil_sky/ stat_coa_rec

      dictionary 'fex_etfl'
      dictionary 'fex_apodl'
      record /fex_etfl/ etfl
      record /fex_apodl/ apodl

      character*72 filename
      integer*4 file_len
      integer*4 vec_pos,check_time(2)

      complex*16 dis_vecs(369,fac_max_num),coa_vec(369)
      complex*16 phase,scale(361)
      real*4 divisor,aifg(512)
      real*8 difg(720)

      fil_coadd = %loc(fil_normal)
c
c     Accumulate good interferogram and number of coadd statistics.
c
      pos = (num_cgr_good - 1)/5 + 1
      num_sci_good(channel,scan_mode,pos) = num_sci_good(channel,scan_mode,pos) 
     .                                      + 1
c
c     Determine whether statistics calculation and subsequent coadd record
c     writing for scan modes FS or FL will be required.  Set appropriate flags
c     if so.
c
      if (scan_mode_input .ge. 5) then
         num_coadds(channel,scan_mode_input) =
     .                      num_coadds(channel,scan_mode_input) + 1
         coa_rec.coad_spec_head.coadd_no = num_coadds(channel,scan_mode_input)
         sec_scan_mode = scan_mode_input
         start_repeat = 2
         stop_repeat = 2
      else if (scan_mode_input .ne. 0) then
         num_coadds(channel,scan_mode) = num_coadds(channel,scan_mode) + 1
         coa_rec.coad_spec_head.coadd_no = num_coadds(channel,scan_mode)
         sec_scan_mode = 0
         start_repeat = 1
         stop_repeat = 1
      else if ((mod(channel,2) .eq. 0) .and. (scan_mode .eq. 2)) then
         num_coadds(channel,scan_mode) = num_coadds(channel,scan_mode) + 1
         num_coadds(channel,5) = num_coadds(channel,5) + 1
         coa_rec.coad_spec_head.coadd_no = num_coadds(channel,scan_mode)
         sec_scan_mode = 5
         start_repeat = 1
         stop_repeat = 2
      else if ((mod(channel,2) .eq. 0) .and. (scan_mode .eq. 4)) then
         num_coadds(channel,scan_mode) = num_coadds(channel,scan_mode) + 1
         num_coadds(channel,6) = num_coadds(channel,6) + 1
         coa_rec.coad_spec_head.coadd_no = num_coadds(channel,scan_mode)
         sec_scan_mode = 6
         start_repeat = 1
         stop_repeat = 2
      else
         num_coadds(channel,scan_mode) = num_coadds(channel,scan_mode) + 1
         coa_rec.coad_spec_head.coadd_no = num_coadds(channel,scan_mode)
         sec_scan_mode = 0
         start_repeat = 1
         stop_repeat = 1
      end if

      coa_rec.coad_spec_head.num_ifgs = num_cgr_good

      call ct_gmt_to_binary(fac_jstop_default,first_time)
      call ct_gmt_to_binary(fac_jstart_default,last_time)
c
c     Initialize local variables.
c
      rec_good = 0
      glwt_tot = 0.0
      xcal_sum = 0.0
      glitch_rate_sum = 0.0
      b_sum = 0.0
      c_sum = 0.0
      b_ave = 0.0
      c_ave = 0.0
c
c     Determine glitch rate weight constants.
c
      if (scan_mode .eq. 3) then
         glitch_intercept = 1.0
         glitch_slope = 0.0
      else
         temp_scan_mode = scan_mode
         if (temp_scan_mode .eq. 4) then
            temp_scan_mode = 3
         end if
         index = (channel-1)*3 + temp_scan_mode
         glitch_intercept = gltchcor.intercept(index)
         glitch_slope = gltchcor.slope(index)
      end if

      do rec = 1,num_cgr_recs
c
c     Coadd science data which passed consistency checking.
c
         if (con_recs(rec).con_check .eq. 0) then
            rec_good = rec_good + 1
c
c     Initialize statistical bins.
c
            bins(rec_good) = 0
c
c     Record the time of the interferogram in the coadd record as days since
c     aperture cover ejection times 100.
c
            coa_rec.coad_spec_head.times(rec_good) = 
     .          sngl(dble(sci_recs(rec).ct_head.time(2) - fac_apco_date) /
     .          fac_vax_year_len) * 36525.0
c
c     Find the earliest and latest science records.
c
            if (single .ne. fac_present) then
               if (time_lt(sci_recs(rec).ct_head.time,first_time)) then
                  first_time(1) = sci_recs(rec).ct_head.time(1)
                  first_time(2) = sci_recs(rec).ct_head.time(2)
                  first_sci_rec = rec
               end if
               if (time_gt(sci_recs(rec).ct_head.time,last_time)) then
                  last_time(1) = sci_recs(rec).ct_head.time(1)
                  last_time(2) = sci_recs(rec).ct_head.time(2)
                  last_sci_rec = rec
               end if

               time_array(1,rec_good) = 
     .                    sci_recs(rec).collect_time.midpoint_time(1)
               time_array(2,rec_good) = 
     .                    sci_recs(rec).collect_time.midpoint_time(2)
               time_weights(rec_good) = 1
            end if
c
c     Determine the glitch rate weight.
c
            glwt(rec) = 1.0 / (glitch_intercept + 
     .                        (glitch_slope * con_recs(rec).glitch_rate))
            glwt_tot = glwt_tot + glwt(rec)
c
c     Accumulate gain, xcal position, and data and attitude quality.
c
            coa_rec.coad_spec_data.gain_sum = coa_rec.coad_spec_data.gain_sum +
     .              fac_gains(sci_recs(rec).sci_head.gain)
            xcal_sum = xcal_sum + (sci_recs(rec).dq_data.xcal_pos * glwt(rec))
            coa_rec.coad_spec_data.dq_summary_flag =
     .              max(zext(coa_rec.coad_spec_data.dq_summary_flag),
     .                  zext(sci_recs(rec).dq_data.data_quality(110)))
            coa_rec.coad_spec_data.att_summary_flag =
     .              max(zext(coa_rec.coad_spec_data.att_summary_flag),
     .                  zext(sci_recs(rec).dq_data.data_quality(109)))
c
c     Sum the interferogram.
c
            do ifg_pos = 1,512
               coa_rec.coad_data.ifg(ifg_pos) = coa_rec.coad_data.ifg(ifg_pos) +
     .                                          (aifgs(ifg_pos,rec) * glwt(rec))
            end do
c
c     Accumulate glitch rate and "b" and "c" coefficient values.
c
            glitch_rate_sum = glitch_rate_sum + (con_recs(rec).glitch_rate *
     .                        glwt(rec))
c
c     Use the glitch rate for binning and put it in the dispersion vectors.
c
            if (con_recs(rec).glitch_rate .ge. cth.glitch_rate(channel)) then
               bins(rec_good) = ibset(bins(rec_good),4)
            end if
            dis_vecs(369,rec_good) = dcmplx(con_recs(rec).glitch_rate)
            glrt(rec_good) = con_recs(rec).glitch_rate
            b_sum = b_sum + (con_recs(rec).b * glwt(rec))
            c_sum = c_sum + (con_recs(rec).c * glwt(rec))
c
c     Bin by galactic latitude and by interferogram ct_head time.
c
            if (abs(sci_recs(rec).attitude.galactic_latitude) .ge.
     .          cth.abs_galactic_lat) then
               bins(rec_good) = ibset(bins(rec_good),5)
            end if
            if (time_ge(sci_recs(rec).ct_head.time,cth.time_saa)) then
               bins(rec_good) = ibset(bins(rec_good),6)
            end if
         end if
      end do
c
c     Calculate the average time.
c
      if (single .ne. fac_present) then
         status = fut_ave_bin_times(time_array,num_cgr_good,time_weights,
     .            coa_rec.ct_head.time)
         call ct_binary_to_gmt(coa_rec.ct_head.time,coa_rec.ct_head.gmt)
c
c     Put the coadd time into the consistency check records.
c
         do rec = 1,num_cgr_recs
            con_recs(rec).coadd_gmt = coa_rec.ct_head.gmt
            con_recs(rec).coadd_time(1) = coa_rec.ct_head.time(1)
            con_recs(rec).coadd_time(2) = coa_rec.ct_head.time(2)
        end do
c
c     Record the various times of the earliest and latest good science records.
c
         coa_rec.coad_spec_head.first_gmt = sci_recs(first_sci_rec).ct_head.gmt
         coa_rec.coad_spec_head.first_time(1) = 
     .           sci_recs(first_sci_rec).ct_head.time(1)
         coa_rec.coad_spec_head.first_time(2) = 
     .           sci_recs(first_sci_rec).ct_head.time(2)
         do pos = 1,6
            coa_rec.coad_spec_head.first_space_time(pos) =
     .              sci_recs(first_sci_rec).ct_head.space_time(pos)
         end do
         coa_rec.coad_spec_head.first_mjr_frm_no = 
     .           sci_recs(first_sci_rec).ct_head.mjr_frm_no

         coa_rec.coad_spec_head.last_gmt = sci_recs(last_sci_rec).ct_head.gmt
         coa_rec.coad_spec_head.last_time(1) = 
     .           sci_recs(last_sci_rec).ct_head.time(1)
         coa_rec.coad_spec_head.last_time(2) = 
     .           sci_recs(last_sci_rec).ct_head.time(2)
         do pos = 1,6
            coa_rec.coad_spec_head.last_space_time(pos) =
     .              sci_recs(last_sci_rec).ct_head.space_time(pos)
         end do
         coa_rec.coad_spec_head.last_mjr_frm_no = 
     .           sci_recs(last_sci_rec).ct_head.mjr_frm_no
      end if
c
c     Record the glitch rate weight total.
c
      coa_rec.coad_spec_head.adj_num_ifgs = glwt_tot
c
c     Coadd the interferogram, xcal position, and glitch rate.
c
      do ifg_pos=1,512
         coa_rec.coad_data.ifg(ifg_pos) =
     .           coa_rec.coad_data.ifg(ifg_pos)/glwt_tot
      end do

      coa_rec.coad_spec_data.xcal_pos = nint(xcal_sum/glwt_tot)
      coa_rec.coad_spec_data.glitch_rate = glitch_rate_sum/glwt_tot
      coa_vec(369) = dcmplx(coa_rec.coad_spec_data.glitch_rate)
c
c     Calculate the averages and variances of the "b" and "c" coefficients.
c
      b_ave = b_sum/glwt_tot
      coa_rec.coad_spec_data.sec_template.b_average = b_ave
      c_ave = c_sum/glwt_tot
      coa_rec.coad_spec_data.transient.c_average = c_ave
      b_sum = 0.0
      c_sum = 0.0
      do rec = 1,num_cgr_recs
         if (con_recs(rec).con_check .eq. 0) then
            b_sum = b_sum + ((con_recs(rec).b - b_ave)**2 * glwt(rec))
            c_sum = c_sum + ((con_recs(rec).c - c_ave)**2 * glwt(rec))
         end if
      end do
      coa_rec.coad_spec_data.sec_template.b_variance = b_sum/glwt_tot
      coa_rec.coad_spec_data.transient.c_variance = c_sum/glwt_tot
c
c     If processing single interferograms, record the attitude directly.
c
      if ((single .eq. fac_present) .or. (num_cgr_good .eq. 1)) then
         coa_rec.attitude.pixel_no = sci_recs(1).attitude.pixel_no
         do pos = 1,3
            coa_rec.attitude.equatorial(pos) = 
     .              sci_recs(1).attitude.equatorial(pos)
         end do
         coa_rec.attitude.ra = sci_recs(1).attitude.ra
         coa_rec.attitude.dec = sci_recs(1).attitude.dec
         coa_rec.attitude.terr_pixel_no = sci_recs(1).attitude.terr_pixel_no
         coa_rec.attitude.terr_latitude = sci_recs(1).attitude.terr_latitude
         coa_rec.attitude.terr_longitude = sci_recs(1).attitude.terr_longitude
         coa_rec.attitude.earth_limb = sci_recs(1).attitude.earth_limb
         coa_rec.attitude.earth_limb_azimuth = 
     .           sci_recs(1).attitude.earth_limb_azimuth
         coa_rec.attitude.sun_angle= sci_recs(1).attitude.sun_angle
         coa_rec.attitude.moon_angle = sci_recs(1).attitude.moon_angle
         coa_rec.attitude.moon_az_angle = sci_recs(1).attitude.moon_az_angle
         coa_rec.attitude.moon_phase = sci_recs(1).attitude.moon_phase
         coa_rec.attitude.sun_moon_dist = sci_recs(1).attitude.sun_moon_dist
         coa_rec.attitude.cobe_moon_dist = sci_recs(1).attitude.cobe_moon_dist
         coa_rec.attitude.altitude = sci_recs(1).attitude.altitude
         coa_rec.attitude.projected_barycentric_velocity = 
     .           sci_recs(1).attitude.projected_barycentric_velocity
         coa_rec.attitude.mcilwain_l_param = 
     .           sci_recs(1).attitude.mcilwain_l_param
         coa_rec.attitude.galactic_latitude = 
     .           sci_recs(1).attitude.galactic_latitude
         coa_rec.attitude.galactic_longitude = 
     .           sci_recs(1).attitude.galactic_longitude
         coa_rec.attitude.ecliptic_latitude = 
     .           sci_recs(1).attitude.ecliptic_latitude
         coa_rec.attitude.ecliptic_longitude = 
     .           sci_recs(1).attitude.ecliptic_longitude
         coa_rec.attitude.orbital_phase = sci_recs(1).attitude.orbital_phase
         coa_rec.attitude.projected_geocentric_velocity = 
     .           sci_recs(1).attitude.projected_geocentric_velocity
         coa_rec.attitude.scan_angle = sci_recs(1).attitude.scan_angle
         coa_rec.attitude.sc_rotation_angle = 
     .           sci_recs(1).attitude.sc_rotation_angle
         coa_rec.attitude.solution= sci_recs(1).attitude.solution
         coa_rec.attitude.pixel_definition = 
     .           sci_recs(1).attitude.pixel_definition
         coa_rec.attitude.skymap_index = sci_recs(1).attitude.skymap_index
         coa_rec.attitude.exc_galactic_lat = 
     .           sci_recs(1).attitude.exc_galactic_lat
         coa_rec.attitude.terr_rad_byte = sci_recs(1).attitude.terr_rad_byte

      else
c
c     For a coadd group, initialize the local attitude variables.
c
         do pos = 1,3
            sum_equat(pos) = 0.0
            sum_terr(pos) = 0.0
            sum_elimb(pos) = 0.0
            sum_mnang(pos) = 0.0
            sum_gal(pos) = 0.0
            sum_ecl(pos) = 0.0
         end do
         sum_sun_angle = 0.0
         sum_moon_phase = 0.0
         sum_sm_dist = 0.0
         sum_cm_dist = 0.0
         sum_altitude = 0.0
         sum_baryvel = 0.0
         sum_lparam = 0.0
         sum_orb_phase = 0.0
         sum_geovel = 0.0
         sum_scan_angle = 0.0
         sum_sc_rot_angle = 0.0
         sum_exc_gal_lat = 0.0
         do pos = 0,6
            sum_att_sol(pos) = 0.0
         end do
         do pos = 1,4
            sum_pixel_no(pos) = 0.0
            sum_pix_def(pos) = 0.0
            sum_sky_idx(pos) = 0.0
         end do
c
c     Accumulate attitude summations over all good interferograms in the group.
c
         rec_good = 0
         do rec = 1,num_cgr_recs
            if (con_recs(rec).con_check .eq. 0) then
               rec_good = rec_good + 1
c
c     Accumulate equatorial coordinates.
c
               do pos = 1,3
                  sum_equat(pos) = sum_equat(pos) + 
     .                (sci_recs(rec).attitude.equatorial(pos) * glwt(rec))
               end do
c
c     Accumulate terrestrial coordinates.
c
               terr_lat = sci_recs(rec).attitude.terr_latitude * fac_att_conv
               terr_long = sci_recs(rec).attitude.terr_longitude * fac_att_conv
               sum_terr(1) = sum_terr(1) + 
     .                       (cosd(terr_lat) * cosd(terr_long) * glwt(rec))
               sum_terr(2) = sum_terr(2) + 
     .                       (cosd(terr_lat) * sind(terr_long) * glwt(rec))
               sum_terr(3) = sum_terr(3) + (sind(terr_lat) * glwt(rec))
c
c     Accumulate earth limb information.
c
               elimb = sci_recs(rec).attitude.earth_limb * fac_att_conv
               el_az = sci_recs(rec).attitude.earth_limb_azimuth * fac_att_conv
               sum_elimb(1) = sum_elimb(1) + 
     .                        (sind(elimb) * cosd(el_az) * glwt(rec))
               sum_elimb(2) = sum_elimb(2) + 
     .                        (sind(elimb) * sind(el_az) * glwt(rec))
               sum_elimb(3) = sum_elimb(3) + (cosd(elimb) * glwt(rec))
c
c     Accumulate sun and moon information.
c
               sum_sun_angle = sum_sun_angle + 
     .                         (sci_recs(rec).attitude.sun_angle * glwt(rec))

               moon_angle = sci_recs(rec).attitude.moon_angle*fac_att_conv
               moon_az = sci_recs(rec).attitude.moon_az_angle * fac_att_conv
               sum_mnang(1) = sum_mnang(1) + 
     .                        (sind(moon_angle) * cosd(moon_az) * glwt(rec))
               sum_mnang(2) = sum_mnang(2) + 
     .                        (sind(moon_angle) * sind(moon_az) * glwt(rec))
               sum_mnang(3) = sum_mnang(3) + (cosd(moon_angle) * glwt(rec))

               moon_phase(rec_good) = sci_recs(rec).attitude.moon_phase * 
     .                                fac_att_conv_rad
               corr_angle = moon_phase(rec) - moon_phase(1)
               if (corr_angle .gt. fac_pi) then
                  corr_angle = corr_angle - 2.0 * fac_pi
               else if (corr_angle .lt. -fac_pi) then
                  corr_angle = corr_angle + 2.0 * fac_pi
               end if
               sum_moon_phase = sum_moon_phase + (corr_angle * glwt(rec))

               sum_sm_dist = sum_sm_dist + 
     .                       (sci_recs(rec).attitude.sun_moon_dist * glwt(rec))
               sum_cm_dist = sum_cm_dist + 
     .                       (sci_recs(rec).attitude.cobe_moon_dist * glwt(rec))
c
c     Accumulate miscellaneous attitude parameters.
c
               sum_altitude = sum_altitude + 
     .                        (sci_recs(rec).attitude.altitude * glwt(rec))
               sum_baryvel = sum_baryvel + 
     .         (sci_recs(rec).attitude.projected_barycentric_velocity*glwt(rec))
               sum_lparam = sum_lparam + 
     .             (sci_recs(rec).attitude.mcilwain_l_param * glwt(rec))
c
c     Accumulate galactic and ecliptic latitude and longitude.
c
               gal_lat = sci_recs(rec).attitude.galactic_latitude * fac_att_conv
               gal_long = sci_recs(rec).attitude.galactic_longitude*fac_att_conv
               sum_gal(1) = sum_gal(1) + 
     .                      (cosd(gal_lat) * cosd(gal_long) * glwt(rec))
               sum_gal(2) = sum_gal(2) + 
     .                      (cosd(gal_lat) * sind(gal_long) * glwt(rec))
               sum_gal(3) = sum_gal(3) + (sind(gal_lat) * glwt(rec))

               ecl_lat = sci_recs(rec).attitude.ecliptic_latitude * fac_att_conv
               ecl_long = sci_recs(rec).attitude.ecliptic_longitude*fac_att_conv
               sum_ecl(1) = sum_ecl(1) + 
     .                      (cosd(ecl_lat) * cosd(ecl_long) * glwt(rec))
               sum_ecl(2) = sum_ecl(2) + 
     .                      (cosd(ecl_lat) * sind(ecl_long) * glwt(rec))
               sum_ecl(3) = sum_ecl(3) + (sind(ecl_lat) * glwt(rec))
c
c     Accumulate miscellaneous attitude parameters.
c
               orb_phase(rec_good) = sci_recs(rec).attitude.orbital_phase * 
     .                               fac_att_conv_rad
               corr_angle = orb_phase(rec) - orb_phase(1)
               if (corr_angle .gt. fac_pi) then
                  corr_angle = corr_angle - 2.0 * fac_pi
               else if (corr_angle .lt. -fac_pi) then
                  corr_angle = corr_angle + 2.0 * fac_pi
               end if
               sum_orb_phase = sum_orb_phase + (corr_angle * glwt(rec))

               sum_geovel = sum_geovel + 
     .         (sci_recs(rec).attitude.projected_geocentric_velocity*glwt(rec))

               scan_angle(rec_good) = sci_recs(rec).attitude.scan_angle * 
     .                                fac_att_conv_rad
               corr_angle = scan_angle(rec) - scan_angle(1)
               if (corr_angle .gt. fac_pi) then
                  corr_angle = corr_angle - 2.0 * fac_pi
               else if (corr_angle .lt. -fac_pi) then
                  corr_angle = corr_angle + 2.0 * fac_pi
               end if
               sum_scan_angle = sum_scan_angle + (corr_angle * glwt(rec))

               sum_sc_rot_angle = sum_sc_rot_angle + 
     .             (sci_recs(rec).attitude.sc_rotation_angle * glwt(rec))
c
c     Accumulate attitude solution type.
c
               int_val = sci_recs(rec).attitude.solution
               if ((int_val .lt. 0) .or. (int_val .gt. 6)) then
                  int_val = 0
               end if
               sum_att_sol(int_val) = sum_att_sol(int_val) + glwt(rec)
c
c     Accumulate pixelization information.
c
               if (sci_recs(rec).attitude.pixel_definition .eq. 'O') then
                  sum_pixel_no(2) = sum_pixel_no(2) + 
     .                      (sci_recs(rec).attitude.pixel_no * glwt(rec))
                  sum_pix_def(2) = sum_pix_def(2) + glwt(rec)
                  sum_sky_idx(2) = sum_sky_idx(2) + 
     .                    (sci_recs(rec).attitude.skymap_index * glwt(rec))
               else if (sci_recs(rec).attitude.pixel_definition .eq. 'S') then
                  sum_pixel_no(3) = sum_pixel_no(3) + 
     .                      (sci_recs(rec).attitude.pixel_no * glwt(rec))
                  sum_pix_def(3) = sum_pix_def(3) + glwt(rec)
                  sum_sky_idx(3) = sum_sky_idx(3) + 
     .                    (sci_recs(rec).attitude.skymap_index * glwt(rec))
               else if (sci_recs(rec).attitude.pixel_definition .eq. 'E') then
                  sum_pix_def(4) = sum_pix_def(4) + glwt(rec)
                  sum_sky_idx(4) = sum_sky_idx(4) + 
     .                    (sci_recs(rec).attitude.skymap_index * glwt(rec))
               else
                  sum_pixel_no(1) = sum_pixel_no(1) + 
     .                      (sci_recs(rec).attitude.pixel_no * glwt(rec))
                  sum_pix_def(1) = sum_pix_def(1) + glwt(rec)
                  sum_sky_idx(1) = sum_sky_idx(1) + 
     .                    (sci_recs(rec).attitude.skymap_index * glwt(rec))
               end if
c
c     Accumulate miscellaneous attitude parameters.
c
               sum_exc_gal_lat = sum_exc_gal_lat +
     .                 (sci_recs(rec).attitude.exc_galactic_lat * glwt(rec))
               coa_rec.attitude.terr_rad_byte = 
     .                 max(zext(coa_rec.attitude.terr_rad_byte),
     .                 zext(sci_recs(rec).attitude.terr_rad_byte))
            end if
         end do
c
c     Calculate coadded attitude if possible.
c
         if ((sum_pixel_no(1) .lt. 0.0) .or. ((sum_equat(1) .eq. 0.0) .and. 
     .       (sum_equat(2) .eq. 0.0) .and. (sum_equat(3) .eq. 0.0))) then
            coa_rec.attitude.pixel_no = -1
            coa_rec.attitude.terr_pixel_no = -1
         else
            norm_equat = sqrt(sum_equat(1)**2+sum_equat(2)**2+sum_equat(3)**2)
            norm_terr = sqrt(sum_terr(1)**2 + sum_terr(2)**2 + sum_terr(3)**2)
            norm_elimb = sqrt(sum_elimb(1)**2+sum_elimb(2)**2+sum_elimb(3)**2)
            norm_mnang = sqrt(sum_mnang(1)**2+sum_mnang(2)**2+sum_mnang(3)**2)
            norm_gal = sqrt(sum_gal(1)**2 + sum_gal(2)**2 + sum_gal(3)**2)
            norm_ecl = sqrt(sum_ecl(1)**2 + sum_ecl(2)**2 + sum_ecl(3)**2)
c
c     Coadd equatorial, ecliptic, galactic, earth limb, moon angle, and
c     terrestrial information.
c
            do pos = 1,3
               if (norm_equat .ne. 0.0) then
                  coa_rec.attitude.equatorial(pos) = sum_equat(pos) / norm_equat
               end if
               if (norm_ecl .ne. 0.0) then
                  avg_ecl(pos) = sum_ecl(pos) / norm_ecl
               end if
               if (norm_gal .ne. 0.0) then
                  avg_gal(pos) = sum_gal(pos) / norm_gal
               end if
               if (norm_elimb .ne. 0.0) then
                  avg_elimb(pos) = sum_elimb(pos) / norm_elimb
               end if
               if (norm_mnang .ne. 0.0) then
                  avg_mnang(pos) = sum_mnang(pos) / norm_mnang
                end if
               if (norm_terr .ne. 0.0) then
                  avg_terr(pos) = sum_terr(pos) / norm_terr
               end if
            end do
c
c     Calculate pixel number.
c
            call xcc_q_to_e(coa_rec.attitude.equatorial,fac_epoch,conv_ecl)
            call firas_pixno(conv_ecl,coa_rec.attitude.pixel_no)
c
c     Calculate ra and dec.
c
            ra = 0.0
            if (coa_rec.attitude.equatorial(1) .ne. 0.0) then
               ra = atan2(coa_rec.attitude.equatorial(2),
     .                    coa_rec.attitude.equatorial(1))
            end if
            coa_rec.attitude.ra = nint(ra*10000.0)
            dec = asin(coa_rec.attitude.equatorial(3))
            coa_rec.attitude.dec = nint(dec*10000.0)
c
c     Coadd terrestrial latitude and longitude.
c
            terr_lat = atan2d(avg_terr(3),sqrt(avg_terr(1)**2 + avg_terr(2)**2))
            terr_long = 0.0
            if (avg_terr(1) .ne. 0.0) then
               terr_long = atan2d(avg_terr(2),avg_terr(1))
            end if
            coa_rec.attitude.terr_latitude = terr_lat / fac_att_conv
            coa_rec.attitude.terr_longitude = terr_long / fac_att_conv
            call firas_pixno(avg_terr,coa_rec.attitude.terr_pixel_no)
c
c     Coadd earth limb information.
c
            elimb = atan2d(sqrt(avg_elimb(1)**2+avg_elimb(2)**2),avg_elimb(3))
            el_az = 0.0
            if (avg_elimb(1) .ne. 0.0) then
               el_az = atan2d(avg_elimb(2),avg_elimb(1))
            end if
            coa_rec.attitude.earth_limb =  elimb / fac_att_conv
            coa_rec.attitude.earth_limb_azimuth = el_az / fac_att_conv
c
c     Coadd sun and moon information.
c
            coa_rec.attitude.sun_angle = sum_sun_angle / glwt_tot

            moon_angle = atan2d(sqrt(avg_mnang(1)**2+avg_mnang(2)**2),
     .                          avg_mnang(3))
            moon_az = 0.0
            if (avg_mnang(1) .ne. 0.0) then
               moon_az = atan2d(avg_mnang(2),avg_mnang(1))
            end if
            coa_rec.attitude.moon_angle = moon_angle / fac_att_conv
            coa_rec.attitude.moon_az_angle = moon_az / fac_att_conv

            avg_angle = (sum_moon_phase / glwt_tot) + moon_phase(1)
            if (avg_angle .gt. fac_pi) then
               avg_angle = avg_angle - 2.0 * fac_pi
            else if (avg_angle .lt. -fac_pi) then
               avg_angle = avg_angle + 2.0 * fac_pi
            end if

            coa_rec.attitude.moon_phase = nint(avg_angle / fac_att_conv_rad)

            coa_rec.attitude.sun_moon_dist = sum_sm_dist / glwt_tot
            coa_rec.attitude.cobe_moon_dist = sum_cm_dist / glwt_tot
c
c     Coadd miscellaneous attitude parameters.
c
            coa_rec.attitude.altitude = sum_altitude / glwt_tot
            coa_rec.attitude.projected_barycentric_velocity = 
     .                    sum_baryvel / glwt_tot
            coa_rec.attitude.mcilwain_l_param = sum_lparam / glwt_tot
c
c     Coadd galactic and ecliptic latitude and longitude.
c
            gal_lat = atan2d(avg_gal(3),sqrt(avg_gal(1)**2 + avg_gal(2)**2))
            gal_long = 0.0
            if (avg_gal(1) .ne. 0.0) then
               gal_long = atan2d(avg_gal(2),avg_gal(1))
            end if
            coa_rec.attitude.galactic_latitude = gal_lat / fac_att_conv
            coa_rec.attitude.galactic_longitude = gal_long / fac_att_conv

            ecl_lat = atan2d(avg_ecl(3),sqrt(avg_ecl(1)**2 + avg_ecl(2)**2))
            ecl_long = 0.0
            if (avg_ecl(1) .ne. 0.0) then
               ecl_long = atan2d(avg_ecl(2),avg_ecl(1))
            end if
            coa_rec.attitude.ecliptic_latitude = ecl_lat / fac_att_conv
            coa_rec.attitude.ecliptic_longitude = ecl_long / fac_att_conv
c
c     Coadd miscellaneous attitude parameters.
c
            avg_angle = (sum_orb_phase / glwt_tot) + orb_phase(1)
            if (avg_angle .gt. fac_pi) then
               avg_angle = avg_angle - 2.0 * fac_pi
            else if (avg_angle .lt. -fac_pi) then
               avg_angle = avg_angle + 2.0 * fac_pi
            end if
            coa_rec.attitude.orbital_phase = nint(avg_angle/fac_att_conv_rad)

            coa_rec.attitude.projected_geocentric_velocity = sum_geovel/glwt_tot

            avg_angle = (sum_scan_angle / glwt_tot) + scan_angle(1)
            if (avg_angle .gt. fac_pi) then
               avg_angle = avg_angle - 2.0 * fac_pi
            else if (avg_angle .lt. -fac_pi) then
               avg_angle = avg_angle + 2.0 * fac_pi
            end if
            coa_rec.attitude.scan_angle = nint(avg_angle / fac_att_conv_rad)

            coa_rec.attitude.sc_rotation_angle = sum_sc_rot_angle / glwt_tot
c
c     Coadd attitude solution type.
c
            int_most = -1.0
            do pos = 0,6
               if (sum_att_sol(pos) .gt. int_most) then
                  int_most = sum_att_sol(pos)
                  int_val = pos
               end if
            end do

            coa_rec.attitude.solution = int_val
c
c     Coadd pixelization information.
c
            int_most = -1.0
            do pos = 1,4
               if (sum_pix_def(pos) .gt. int_most) then
                  int_most = sum_pix_def(pos)
                  int_val = pos
               end if
            end do

            if (int_val .eq. 1) then
               coa_rec.attitude.pixel_definition = 'q'
               coa_rec.attitude.skymap_index = 
     .                 nint(sum_sky_idx(1)/sum_pix_def(1))
               if (coa_rec.attitude.skymap_index .lt. 5) then
                  coa_rec.attitude.pixel_no = 
     .                    nint(sum_pixel_no(1)/sum_pix_def(1))
               end if
            else if (int_val .eq. 2) then
               coa_rec.attitude.pixel_definition = 'O'
               coa_rec.attitude.pixel_no = 
     .                 nint(sum_pixel_no(2)/sum_pix_def(2))
               coa_rec.attitude.skymap_index = 
     .                 nint(sum_sky_idx(2)/sum_pix_def(2))
            else if (int_val .eq. 3) then
               coa_rec.attitude.pixel_definition = 'S'
               coa_rec.attitude.pixel_no = 
     .                 nint(sum_pixel_no(3)/sum_pix_def(3))
               coa_rec.attitude.skymap_index = 
     .                 nint(sum_sky_idx(3)/sum_pix_def(3))
            else if (int_val .eq. 4) then
               coa_rec.attitude.pixel_definition = 'E'
               coa_rec.attitude.skymap_index = 
     .                       nint(sum_sky_idx(4)/sum_pix_def(4))
               if (coa_rec.attitude.skymap_index .lt. 5) then
                  coa_rec.attitude.pixel_no = 
     .                    nint(sum_pixel_no(4)/sum_pix_def(4))
               end if
            end if
c
c     Coadd miscellaneous attitude parameters.
c
            coa_rec.attitude.exc_galactic_lat = sum_exc_gal_lat/glwt_tot
         end if
      end if
c
c     If processing single interferograms, record engineering data and set
c     minima and maxima to flag values.
c
      if ((single .eq. fac_present) .or. (num_cgr_good .eq. 1)) then
         do temp_pos = 1,4
            coa_rec.coad_spec_data.bol_volt_min(temp_pos) = -9999.0
            coa_rec.coad_spec_data.bol_volt_max(temp_pos) = -9999.0
         end do
         do temp_pos = 1,2
            coa_rec.en_tempdiff(temp_pos).xcal = 
     .              eng_recs(1).en_tempdiff(temp_pos).xcal
            coa_rec.en_tempdiff(temp_pos).ical = 
     .              eng_recs(1).en_tempdiff(temp_pos).ical
            coa_rec.en_tempdiff(temp_pos).skyhorn = 
     .              eng_recs(1).en_tempdiff(temp_pos).skyhorn
            coa_rec.en_tempdiff(temp_pos).refhorn = 
     .              eng_recs(1).en_tempdiff(temp_pos).refhorn
            coa_rec.en_tempdiff(temp_pos).dihedral = 
     .              eng_recs(1).en_tempdiff(temp_pos).dihedral
            coa_rec.en_tempdiff(temp_pos).collimator_mirror = 
     .              eng_recs(1).en_tempdiff(temp_pos).collimator_mirror
            do eng_pos = 1,4
               coa_rec.en_tempdiff(temp_pos).bol_assem(eng_pos) = 
     .                 eng_recs(1).en_tempdiff(temp_pos).bol_assem(eng_pos)
            end do
         end do
         do temp_pos = 1,10
            coa_rec.coad_spec_data.temp_min(temp_pos) = -9999.0
            coa_rec.coad_spec_data.temp_max(temp_pos) = -9999.0
         end do
      else
c
c     If processing a coadd group, find the minima and maxima of engineering
c     quantities.
c
         do eng_pos = 1,4
            coa_rec.coad_spec_data.bol_volt_min(eng_pos) = 1000.0
         end do
         do temp_pos = 1,10
            coa_rec.coad_spec_data.temp_min(temp_pos) = 1000.0
         end do

         rec_good = 0
         do rec = 1,num_cgr_recs
            if (con_recs(rec).con_check .eq. 0) then
               rec_good = rec_good + 1
               do eng_pos = 1,4
                   coa_rec.coad_spec_data.bol_volt_min(eng_pos) = 
     .                     min(coa_rec.coad_spec_data.bol_volt_min(eng_pos),
     .                     eng_recs(rec).en_analog.bol_volt(eng_pos))
                   coa_rec.coad_spec_data.bol_volt_max(eng_pos) = 
     .                     max(coa_rec.coad_spec_data.bol_volt_max(eng_pos),
     .                     eng_recs(rec).en_analog.bol_volt(eng_pos))
               end do
               do eng_pos = 1,2
                  coa_rec.en_tempdiff(eng_pos).xcal = 
     .                    max(coa_rec.en_tempdiff(eng_pos).xcal,
     .                    abs(eng_recs(rec).en_tempdiff(eng_pos).xcal))
                  coa_rec.en_tempdiff(eng_pos).ical = 
     .                    max(coa_rec.en_tempdiff(eng_pos).ical,
     .                    abs(eng_recs(rec).en_tempdiff(eng_pos).ical))
                  coa_rec.en_tempdiff(eng_pos).skyhorn = 
     .                    max(coa_rec.en_tempdiff(eng_pos).skyhorn,
     .                    abs(eng_recs(rec).en_tempdiff(eng_pos).skyhorn))
                  coa_rec.en_tempdiff(eng_pos).refhorn = 
     .                    max(coa_rec.en_tempdiff(eng_pos).refhorn,
     .                    abs(eng_recs(rec).en_tempdiff(eng_pos).refhorn))
                  coa_rec.en_tempdiff(eng_pos).collimator_mirror = 
     .                    max(coa_rec.en_tempdiff(eng_pos).collimator_mirror,
     .                   abs(eng_recs(rec).en_tempdiff(eng_pos).collimator_mirror))
                  coa_rec.en_tempdiff(eng_pos).dihedral = 
     .                    max(coa_rec.en_tempdiff(eng_pos).dihedral,
     .                    abs(eng_recs(rec).en_tempdiff(eng_pos).dihedral))
                  do temp_pos = 1,4
                     coa_rec.en_tempdiff(eng_pos).bol_assem(temp_pos) = 
     .                       max(coa_rec.en_tempdiff(eng_pos).bol_assem(temp_pos),
     .                  abs(eng_recs(rec).en_tempdiff(eng_pos).bol_assem(temp_pos)))
                  end do
               end do
               do temp_pos = 1,10
                  coa_rec.coad_spec_data.temp_min(temp_pos) = 
     .                    min(coa_rec.coad_spec_data.temp_min(temp_pos),
     .                    temps(temp_pos,rec))
                  coa_rec.coad_spec_data.temp_max(temp_pos) = 
     .                    max(coa_rec.coad_spec_data.temp_max(temp_pos),
     .                    temps(temp_pos,rec))
               end do
c
c     Use combined engineering temperatures for statistical binning, and 
c     put them into the dispersion vectors.
c
               dis_vecs(362,rec_good) = dcmplx(temps(1,rec))
               if (temps(1,rec) .ge. cth.xcal_temp) then
                  bins(rec_good) = ibset(bins(rec_good),0)
               end if
               dis_vecs(363,rec_good) = dcmplx(temps(2,rec))
               if (temps(2,rec) .ge. cth.ical_temp) then
                   bins(rec_good) = ibset(bins(rec_good),1)
               end if
               val = (temps(3,rec) + temps(4,rec)) / 2.0
               dis_vecs(364,rec_good) = dcmplx(val)
               if (val .ge. cth.skyhorn_refhorn_avg) then
                  bins(rec_good) = ibset(bins(rec_good),2)
               end if
               if (temps(5,rec) .ge. cth.dihedral_temp) then
                  bins(rec_good) = ibset(bins(rec_good),7)
               end if
               val = (temps(3,rec) - temps(4,rec)) / 2.0
               dis_vecs(365,rec_good) = dcmplx(val)
               dis_vecs(366,rec_good) = dcmplx(temps(5,rec))
               dis_vecs(367,rec_good) = dcmplx(temps(6,rec))
               dis_vecs(368,rec_good) = dcmplx(temps(6+channel,rec))
               if (temps(6+channel,rec) .ge. cth.bolometer_temp(channel)) then
                  bins(rec_good) = ibset(bins(rec_good),3)
               end if
c
c     Increment each statistical bin value to set the range of values to 1-256.
c
               bins(rec_good) = bins(rec_good) + 1
               coa_rec.coad_spec_data.bin_info(bins(rec_good)) = 
     .                 coa_rec.coad_spec_data.bin_info(bins(rec_good)) + 1   
            end if
         end do
      end if
c
c     Coadd the engineering status bits, hot spot commands, and power statuses;
c     record them if processing single interferograms.
c
      if ((single .eq. fac_present) .or. (num_cgr_good .eq. 1)) then
         do eng_pos=1,16
            coa_rec.en_stat.group1(eng_pos)=eng_recs(1).en_stat.group1(eng_pos)
         end do
         do eng_pos=1,14
            coa_rec.en_stat.group2(eng_pos)=eng_recs(1).en_stat.group2(eng_pos)
         end do
         do eng_pos=1,2
            coa_rec.en_stat.hot_spot_cmd(eng_pos) = 
     .              eng_recs(1).en_stat.hot_spot_cmd(eng_pos)
            coa_rec.en_stat.power_a_status(eng_pos) = 
     .              eng_recs(1).en_stat.power_a_status(eng_pos)
            coa_rec.en_stat.power_b_status(eng_pos) = 
     .              eng_recs(1).en_stat.power_b_status(eng_pos)
         end do
      else
c
c     Coadd engineering status bits.
c 
        do eng_pos=1,16
            sum = 0.0
            do rec = 1,num_cgr_recs
               if (con_recs(rec).con_check .eq. 0) then
                  sum = sum + (eng_recs(rec).en_stat.group1(eng_pos)*glwt(rec))
               end if
            end do
            coa_rec.en_stat.group1(eng_pos) = nint(sum / glwt_tot)
         end do

         do eng_pos=1,14
            sum = 0.0
            do rec = 1,num_cgr_recs
               if (con_recs(rec).con_check .eq. 0) then
                  sum = sum + (eng_recs(rec).en_stat.group2(eng_pos)*glwt(rec))
               end if
            end do
            coa_rec.en_stat.group2(eng_pos) = nint(sum / glwt_tot)
         end do
c
c     Coadd hot spot command and power statuses.
c
         do eng_pos = 1,6
            sum_eng(eng_pos) = 0.0
         end do
         do rec = 1,num_cgr_recs
            if (con_recs(rec).con_check .eq. 0) then
               sum_eng(1) = sum_eng(1) + 
     .             (eng_recs(rec).en_stat.hot_spot_cmd(1) * glwt(rec))
               sum_eng(2) = sum_eng(2) + 
     .             (eng_recs(rec).en_stat.hot_spot_cmd(2) * glwt(rec))
               sum_eng(3) = sum_eng(3) + 
     .             (eng_recs(rec).en_stat.power_a_status(1) * glwt(rec))
               sum_eng(4) = sum_eng(4) + 
     .             (eng_recs(rec).en_stat.power_a_status(2) * glwt(rec))
               sum_eng(5) = sum_eng(5) + 
     .             (eng_recs(rec).en_stat.power_b_status(1) * glwt(rec))
               sum_eng(6) = sum_eng(6) + 
     .             (eng_recs(rec).en_stat.power_b_status(2) * glwt(rec))
            end if
         end do
         coa_rec.en_stat.hot_spot_cmd(1) = nint(sum_eng(1) / glwt_tot)
         coa_rec.en_stat.hot_spot_cmd(2) = nint(sum_eng(2) / glwt_tot)
         coa_rec.en_stat.power_a_status(1) = nint(sum_eng(3) / glwt_tot)
         coa_rec.en_stat.power_a_status(2) = nint(sum_eng(4) / glwt_tot)
         coa_rec.en_stat.power_b_status(1) = nint(sum_eng(5) / glwt_tot)
         coa_rec.en_stat.power_b_status(2) = nint(sum_eng(6) / glwt_tot)
      end if

      coa_rec.coad_spec_data.bol_cmd_bias=coa_rec.en_stat.bol_cmd_bias(channel)
c
c     Coadd grt temperatures and remaining engineering quantities and 
c     calculate sigmas; record them and set flag values if processing single 
c     interferograms.
c
      if ((single .eq. fac_present) .or. (num_cgr_good .eq. 1)) then
         do eng_pos = 1,64
            coa_rec.en_analog.grt(eng_pos) = eng_recs(1).en_analog.grt(eng_pos)
         end do
         do eng_pos = 1,62
            coa_rec.en_analog.group1(eng_pos) = 
     .              eng_recs(1).en_analog.group1(eng_pos)
         end do
         do eng_pos = 1,64
            coa_rec.en_sigma.sig_grt(eng_pos) = -9999.0
         end do
         do eng_pos = 1,62
            coa_rec.en_sigma.group2(eng_pos) = -9999.0
         end do
      else
c
c     Coadd grt temperatures.
c
         do eng_pos = 1,64
            rec_good = 0
            sum = 0.0
            rec_glwt = 0.0
            do rec = 1,num_cgr_recs
               if (con_recs(rec).con_check .eq. 0) then
                  val = eng_recs(rec).en_analog.grt(eng_pos) * glwt(rec)
                  if (val .gt. 0.0) then
                     rec_good = rec_good + 1
                     sum = sum + val
                     rec_glwt = rec_glwt + glwt(rec)
                  end if
               end if
            end do
            if (rec_good .ne. 0) then
               coa_rec.en_analog.grt(eng_pos) = sum/rec_glwt
               dsum = 0.0
               do rec = 1,num_cgr_recs
                  if (con_recs(rec).con_check .eq. 0) then
                     dval = eng_recs(rec).en_analog.grt(eng_pos)
                     if (dval .gt. 0.0) then
                        dsum = dsum + (dval - 
     .                         coa_rec.en_analog.grt(eng_pos))**2 * glwt(rec)
                     end if
                  end if
               end do
               coa_rec.en_sigma.sig_grt(eng_pos) = sqrt(dsum/rec_glwt)
            else
               coa_rec.en_analog.grt(eng_pos) = -9999.0
               coa_rec.en_sigma.sig_grt(eng_pos) = 1.0E15
            end if
         end do
c
c     Coadd remaining engineering analog quantities.
c
         do eng_pos = 1,62
            rec_good = 0
            sum = 0.0
            rec_glwt = 0.0
            do rec = 1,num_cgr_recs
               if (con_recs(rec).con_check .eq. 0) then
                  val = eng_recs(rec).en_analog.group1(eng_pos) * glwt(rec)
                  if (val .gt. 0.0) then
                     rec_good = rec_good + 1
                     sum = sum + val
                     rec_glwt = rec_glwt + glwt(rec)
                  end if
               end if
            end do
            if (rec_good .ne. 0) then
               coa_rec.en_analog.group1(eng_pos) = sum/rec_glwt
               dsum = 0.0
               do rec = 1,num_cgr_recs
                  if (con_recs(rec).con_check .eq. 0) then
                     dval = eng_recs(rec).en_analog.group1(eng_pos)
                     if (dval .gt. 0.0) then
                        dsum = dsum + (dval - 
     .                    coa_rec.en_analog.group1(eng_pos))**2 * glwt(rec)
                     end if
                  end if
               end do
               coa_rec.en_sigma.group2(eng_pos) = sqrt(dsum/rec_glwt)
            else
               coa_rec.en_analog.group1(eng_pos) = -9999.0
               coa_rec.en_sigma.group2(eng_pos) = 1.0E15
            end if
         end do
      end if

      coa_rec.coad_spec_data.bol_volt = coa_rec.en_analog.bol_volt(channel)
c
c     If processing single interferograms, set combined temperature sigmas to
c     a flag value, otherwise, call fut_temp_list to calculate the sigmas.
c
      if ((single .eq. fac_present) .or. (num_cgr_good .eq. 1)) then
         singlifg = .true.
         do temp_pos = 1,10
            coa_rec.coad_spec_data.temp_sigma(temp_pos) = -9999.0
         end do
      else
         singlifg = .false.
      end if

      status = fut_temp_list(coa_rec.en_analog,coa_rec.en_sigma,
     .                       grtcoawt,grttrans,0,singlifg,
     .                       coa_rec.coad_spec_data.temp,
     .                       coa_rec.coad_spec_data.temp_sigma)
c
c     Construct the coadd label.
c
      if (coa_rec.coad_spec_data.xcal_pos .eq. fac_xcalout) then
         xcal_pos = 'O'
      else
         xcal_pos = 'I'
      end if

      if (coa_rec.attitude.pixel_no .ge. 0) then
         encode (60,10,coa_rec.coad_spec_head.label)
     .           coa_rec.coad_spec_head.first_gmt,
     .           coa_rec.coad_spec_head.last_gmt,
     .           coa_rec.coad_spec_head.num_ifgs,
     .           fac_scan_mode_idsl(scan_mode),
     .           nint(coa_rec.coad_spec_data.ical*1000.0),
     .           xcal_pos,coa_rec.attitude.pixel_no
10       format (2(a,'_'),i3.3,'_',a,'_',i5.5,'_',a,'_',i4.4,11x)
      else
         encode (60,20,coa_rec.coad_spec_head.label)
     .           coa_rec.coad_spec_head.first_gmt,
     .           coa_rec.coad_spec_head.last_gmt,
     .           coa_rec.coad_spec_head.num_ifgs,
     .           fac_scan_mode_idsl(scan_mode),
     .           nint(coa_rec.coad_spec_data.ical*1000.0),
     .           xcal_pos
20       format (2(a,'_'),i3.3,'_',a,'_',i5.5,'_',a,16x)
      end if
c
c     If processing single interferograms, set the three statistical vectors
c     in the coadd record to flag values.
c
      if ((single .eq. fac_present) .or. (num_cgr_good .eq. 1)) then
         do temp_pos = 1,361
            coa_rec.coad_data.real_var(temp_pos) = -9999.0
            coa_rec.coad_data.imag_var(temp_pos) = -9999.0
            coa_rec.coad_data.real_imag_var(temp_pos) = -9999.0
            do pos = 1,3
               var_vecs(pos,temp_pos) = -9999.0
            end do
         end do
         coa_rec.coad_data.fft_length = fac_fft_length(scan_mode)
      else
c
c     When processing a coadd group, add the engineering data to the end
c     of the coadd dispersion vector.
c
         coa_vec(362) = dcmplx(coa_rec.coad_spec_data.xcal)
         coa_vec(363) = dcmplx(coa_rec.coad_spec_data.ical)
         coa_vec(364) = dcmplx((coa_rec.coad_spec_data.skyhorn +
     .                   coa_rec.coad_spec_data.refhorn) / 2.0)
         coa_vec(365) = dcmplx((coa_rec.coad_spec_data.skyhorn - 
     .                   coa_rec.coad_spec_data.refhorn) / 2.0)
         coa_vec(366) = dcmplx(coa_rec.coad_spec_data.dihedral)
         coa_vec(367) = dcmplx(coa_rec.coad_spec_data.collimator_mirror)
         coa_vec(368) = dcmplx(coa_rec.coad_spec_data.bol_assem(channel))
c
c     Calculate a divisor value for the statistics based on the template
c     subtraction information.
c         
         if (coa_rec.coad_spec_data.sec_template.subtracted .eq. 
     .       fac_present) then
            divisor = glwt_tot * (num_cgr_good - 2)
         else if (coa_rec.coad_spec_data.prim_template.subtracted .eq.
     .            fac_present) then
            divisor = glwt_tot * (num_cgr_good - 1)
         else
            divisor = glwt_tot
         end if
c
c     Repeat the statistics calculations as needed for FS and FL scan modes.
c
         do repeat = start_repeat,stop_repeat
c
c     If not FS or FL, get fft length.
c
            if (repeat .eq. 1) then
               stat_coa_rec = coa_rec
               fft_len = fac_fft_length(scan_mode)
               coa_rec.coad_data.fft_length = fft_len
               temp_scan_mode = scan_mode
               mtm_speed = mod(temp_scan_mode-1,2)
               temp_adds_per_group = adds_per_group
            else
c
c     Scan mode is FS or FL; decimate or truncate the coadded interferogram
c     and retrieve correct Nyquist frequency and fft length.  Also set adds
c     per group so that the electronics transfer function for low channel SF
c     scan mode is retrieved.
c
               status = fil_short(channel,sec_scan_mode,nyquistl,sci_mode,
     .                            coa_rec,var_vecs,stat_coa_rec)
               fft_len = stat_coa_rec.coad_data.fft_length
               temp_scan_mode = sec_scan_mode
               mtm_speed = 1
               temp_adds_per_group = 2
            end if
c
c     The spectrum length is half the fft length plus one.
c 
            spec_len = (fft_len / 2) + 1
c
c     Find record number of the correct electronics transfer function.
c
            status = fut_get_recnum(fake_it,mtm_speed,channel,sci_mode,
     .                              temp_adds_per_group,etfl_rec)
c
c     Read the electronics transfer function.
c
            read(etfl_lun,rec=etfl_rec,iostat=status) etfl
            if (status .ne. 0) then
               fil_coadd = %loc(fil_cfdirread)
               inquire(etfl_lun,name=filename)
               call str$trim(filename,filename,file_len)
               call lib$signal(fil_cfdirread,%val(2),filename(:file_len),
     .                         %val(status))
               return
            end if
c
c     Find record number of the correct apodization function.  The zero
c     indicates an unlinearized interferogram.
c
            if (repeat .eq. 2) then
               temp_adds_per_group = 8
            end if
            status = fut_apod_recnuml(channel,temp_scan_mode,fake_it,sci_mode,
     .                                temp_adds_per_group,0,apodl_rec)
c
c     Read the apodization function.
c
            read(apodl_lun,rec=apodl_rec,iostat=status) apodl
            if (status .ne. 0) then
               fil_coadd = %loc(fil_cfdirread)
               inquire(apodl_lun,name=filename)
               call str$trim(filename,filename,file_len)
               call lib$signal(fil_cfdirread,%val(2),filename(:file_len),
     .                         %val(status))
               return
            end if
c
c     Construct the scaling vector from the electronics transfer function, the
c     etendu, the adc scale, and the nyquist frequency in icm.  Add a phase
c     shift for high channels and sf scan mode.
c
            do vec_pos = 1,spec_len
               if ((mod(channel,2) .eq. 1) .and. (scan_mode .eq. 2)) then
                  phase = dcmplx(0.d0,fac_dpi/fft_len)
                  if (etfl.ztransfcn(vec_pos) .eq. 0.0) then
                     scale(vec_pos) = exp((spec_len-1) * phase) / (fac_etendu *
     .                                fac_adc_scale * 
     .                                stat_coa_rec.coad_spec_data.nyquist_icm)
                  else
                     scale(vec_pos) = exp((spec_len-1) * phase) / 
     .                                (etfl.ztransfcn(vec_pos) * fac_etendu * 
     .                                fac_adc_scale * 
     .                                stat_coa_rec.coad_spec_data.nyquist_icm)
                  end if
               else
                  if (etfl.ztransfcn(vec_pos) .eq. 0.0) then
                     scale(vec_pos) = 1.0 / (fac_etendu * fac_adc_scale * 
     .                                stat_coa_rec.coad_spec_data.nyquist_icm)
                  else
                     scale(vec_pos) = 1.0 / (etfl.ztransfcn(vec_pos) * 
     .                                fac_etendu * fac_adc_scale * 
     .                                stat_coa_rec.coad_spec_data.nyquist_icm)
                  end if
               end if
            end do
c
c     Call the routine fut_apod_rotl to rotate, pad and apodize the coadded
c     interferogram.  Then apply the Fast Fourier transform.
c
            peak = stat_coa_rec.coad_spec_data.peak_pos
            status = fut_apod_rotl(stat_coa_rec.coad_data.ifg,apodl,fft_len,
     .                             peak,difg)
            call dfftrf(fft_len,difg,difg)
c
c     Set the first point of the coadd dispersion vector to zero.
c
            coa_vec(1) = dcmplx(0.0,0.0)
c
c     Unpack the Fourier transform and multiply by the scaling vector.
c
            do vec_pos = 2,spec_len-1
               coa_vec(vec_pos) = dcmplx(difg(vec_pos*2-2),difg(vec_pos*2-1)) * 
     .                            scale(vec_pos)
            end do
            coa_vec(spec_len) = dcmplx(difg(fft_len)) * scale(spec_len)
         
            rec_good = 0
c
c     Loop over the good interferograms to form the dispersion vectors.
c
            do rec = 1,num_cgr_recs
               if (con_recs(rec).con_check .eq. 0) then
                  rec_good = rec_good + 1
                  do ifg_pos = 1,512
                     aifg(ifg_pos) = aifgs(ifg_pos,rec)
                  end do
c
c     Decimate or truncate the interferogram as necessary for scan mode FS, FL.
c
                  if (repeat .eq. 2) then
                     if (temp_scan_mode .eq. 5) then
                        do ifg_pos = 1,127
                           pos = ifg_pos * 4
                           aifg(ifg_pos) = 
     .                         (aifgs(pos,rec) + 
     .                          aifgs(pos+1,rec) +
     .                          aifgs(pos+2,rec) + 
     .                          aifgs(pos+3,rec)) / 4.0
                        end do
                        aifg(128) = aifgs(512,rec)
                     end if
                     do ifg_pos = 129,512
                        aifg(ifg_pos) = 0.0
                     end do
                  end if
c
c     Call the routine fut_apod_rotl to rotate, pad and apodize the
c     interferogram.  Then apply the Fast Fourier transform.
c
                  status = fut_apod_rotl(aifg,apodl,fft_len,peak,difg)
                  call dfftrf(fft_len,difg,difg)
c        
c     Set the first point of the dispersion vector to zero.
c
                  dis_vecs(1,rec_good) = dcmplx(0.0,0.0)
c
c     Unpack the Fourier transform, multiply by the scaling vector, and
c     subtract the coadd dispersion vector.
c
                  do vec_pos = 2,spec_len-1
                     dis_vecs(vec_pos,rec_good) = 
     .                        dcmplx(difg(vec_pos*2-2),difg(vec_pos*2-1)) *
     .                        scale(vec_pos) - coa_vec(vec_pos)
                  end do
                  dis_vecs(spec_len,rec_good) = dcmplx(difg(fft_len)) * 
     .                     scale(spec_len) - coa_vec(spec_len)
                  do vec_pos = 362,369
                     dis_vecs(vec_pos,rec_good) = dis_vecs(vec_pos,rec_good) -
     .                                            coa_vec(vec_pos)
                  end do
c
c     Calculate the three statistical vectors:  real-real variance, 
c     imaginary-imaginary variance, and real-imaginary variance; leave them
c     set at zero if the divisor value is zero.
c
                  if (divisor .ne. 0.0) then
                     if (repeat .eq. 1) then
                        do vec_pos = 1,spec_len
                           coa_rec.coad_data.real_var(vec_pos) =
     .                             coa_rec.coad_data.real_var(vec_pos) + 
     .                             (dreal(dis_vecs(vec_pos,rec_good))**2 *
     .                             glwt(rec)) / divisor
                           coa_rec.coad_data.imag_var(vec_pos) =
     .                             coa_rec.coad_data.imag_var(vec_pos) + 
     .                             (dimag(dis_vecs(vec_pos,rec_good))**2 * 
     .                             glwt(rec)) / divisor
                           coa_rec.coad_data.real_imag_var(vec_pos) =
     .                             coa_rec.coad_data.real_imag_var(vec_pos) + 
     .                             (dreal(dis_vecs(vec_pos,rec_good)) *
     .                             dimag(dis_vecs(vec_pos,rec_good)) *
     .                             glwt(rec)) / divisor
                        end do
                     else
                        do vec_pos = 1,spec_len	
                           var_vecs(1,vec_pos) = var_vecs(1,vec_pos) +
     .                             (dreal(dis_vecs(vec_pos,rec_good))**2 *
     .                             glwt(rec)) / divisor
                           var_vecs(2,vec_pos) = var_vecs(2,vec_pos) +
     .                             (dimag(dis_vecs(vec_pos,rec_good))**2 * 
     .                             glwt(rec)) / divisor
                           var_vecs(3,vec_pos) = var_vecs(3,vec_pos) +
     .                             (dreal(dis_vecs(vec_pos,rec_good)) *
     .                             dimag(dis_vecs(vec_pos,rec_good)) *
     .                             glwt(rec)) / divisor
                        end do
                     end if
                  end if
               end if
            end do
c
c     Update covariance matrix if required.  Only coadd groups for which no
c     secondary template was subtracted contribute to the covariance matrix.
c     Temporarily reset the scan mode value if it is 4.
c
            if ((covar .eq. fac_present) .and. (temp_scan_mode .ne. 3) .and.
     .          (coa_rec.coad_spec_data.sec_template.subtracted .eq. 
     .           fac_not_present)) then
               if (temp_scan_mode .ge. 4) then
                  temp_scan_mode = temp_scan_mode - 1
               end if
c
c     Compare the earliest and lastest science record times to the current
c     earliest and latest covariance matrix times to see if the latter need
c     to be updated.
c
               check_time(1) = first_covar_times(1,temp_scan_mode)
               check_time(2) = first_covar_times(2,temp_scan_mode)
               if (time_lt(first_time,check_time)) then
                  first_covar_times(1,temp_scan_mode) = first_time(1)
                  first_covar_times(2,temp_scan_mode) = first_time(2)
               end if
         
               check_time(1) = last_covar_times(1,temp_scan_mode)
               check_time(2) = last_covar_times(2,temp_scan_mode)
               if (time_gt(last_time,check_time)) then
                  last_covar_times(1,temp_scan_mode) = last_time(1)
                  last_covar_times(2,temp_scan_mode) = last_time(2)
               end if
c
c     Perform covariance matrix calculations.
c
               status = fil_covar(temp_scan_mode,num_cgr_good,glwt_tot,glrt,
     .                            dis_vecs,sci_mode,temp_adds_per_group,
     .                            coa_vec,bins,cov_recs)
            end if
         end do
      end if
c
c     If the primary or secondary templates were subtracted, add them back to
c     the coadded interferogram.
c
      if (single .ne. fac_present) then
         if (coa_rec.coad_spec_data.prim_template.subtracted .eq. 
     .       fac_present) then
            do ifg_pos = 1,512
               coa_rec.coad_data.ifg(ifg_pos) = coa_rec.coad_data.ifg(ifg_pos) +
     .                                          template(ifg_pos)
            end do
         end if
         if (coa_rec.coad_spec_data.sec_template.subtracted .eq. 
     .       fac_present) then
            do ifg_pos = 1,512
               coa_rec.coad_data.ifg(ifg_pos) = coa_rec.coad_data.ifg(ifg_pos) +
     .                 sec_template(ifg_pos) * 
     .                 coa_rec.coad_spec_data.sec_template.b_average
            end do
         end if
      end if

      return
      end
