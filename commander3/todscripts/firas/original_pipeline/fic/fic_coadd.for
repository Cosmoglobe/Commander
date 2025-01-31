      integer*4 function fic_coadd(input,channel,scan_mode,covar,single,
     .                             grttrans,grtcoawt,cth,etf_lun,num_recs,
     .                             num_good,num_sci_good,sci_recs,eng_recs,
     .                             fake_it,sci_mode,adds_per_group,temps,
     .                             con_recs,num_coadds,coa_rec,aifgs,template,
     .                             sec_template,cov_recs,first_covar_times,
     .                             last_covar_times)
c-------------------------------------------------------------------------------
c
c     Purpose: Coadd science and engineering data and form variance vectors
c              from spectra; call fic_covar to calculate other statistical
c              quantities.
c
c     Authors: R. Kummerer, STI, 6/86
c              S. Alexander, STX, 9/91, SER 7985, 7564
c
c     Input: input              i*4  Type of input: sky, calibration, or raw.
c            channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode          i*4  Value of scan mode, 1-4=SS-LF.
c            covar              i*4  Covariance matrix to be written.
c            single             i*4  Process single ifgs.
c            grttrans           rec  Transition temperatures.
c            grtcoawt           rec  Coadd temperature weights.
c            cth                rec  Consistency check parameters.
c            etf_lun            i*4  Logical unit for electronics transfer
c                                    functions.
c            num_recs           i*4  Number of input science and engineering
c                                    records.
c            num_good           i*4  Current number of good interferograms.
c            num_sci_good       i*4(4,4,20)  Number of good interferogrmas
c                                            contributing to an output coadd 
c                                            record by channel and scan mode; 
c                                            grouped as 1-5, ..., 96-100.
c            sci_recs           rec(num_recs)  Science records.
c            eng_recs           rec(num_recs)  Engineering records.
c            fake_it            i*4  Value of fake-it bit, 0-1.
c            sci_mode           i*4  Value of science mode, 0-4.
c            adds_per_group     i*4  Value of adds per group, 1-12.
c            temps              r*4(10,num_recs)  Combined engineering
c                                                 temperatures.
c            con_recs           rec(num_recs)  Consistency check records.
c            num_coadds         i*4(4,4)  Number of output coadd records by
c                                         channel and scan mode.
c            coa_rec            rec  Coadd record.
c            aifgs              r*4(512,num_recs)  Real-valued interferograms.
c            template           r*4(512)  Primary template.
c            sec_template       r*4(512)  Secondary template.
c            cov_recs           rec(771)  Covariance matrix records.
c            first_covar_times  i*4(2,3)  Binary times of first records
c                                         contributing to covariance matrices.
c            last_covar_times   i*4(2,3)  Binary times of last records
c                                         contributing to covariance matrices.
c
c     Output: num_sci_good       i*4(4,4,20)  Number of good interferogrmas
c                                             contributing to an output coadd 
c                                             record by channel and scan mode; 
c                                             grouped as 1-5, ..., 96-100.
c             num_coadds         i*4(4,4)  Number of output coadd records by
c                                          channel and scan mode.
c             coa_rec            rec  Coadd record.
c             cov_recs           rec(771)  Covariance matrix records.
c             first_covar_times  i*4(2,3)  Binary times of first records
c                                          contributing to covariance matrices.
c             last_covar_times   i*4(2,3)  Binary times of last records
c                                          contributing to covariance matrices.
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
      external fic_normal
      external fic_cfdirread
c
c     Functions.
c
      logical*1 time_lt,time_gt,time_ge
      integer*4 fut_ave_bin_times
      integer*4 fut_average_angles
      integer*4 fut_temp_list
      integer*4 fic_covar
c
c     Input parameters.
c
      integer*4 input,channel,scan_mode
      integer*4 covar,single

      dictionary 'fex_grttrans'
      dictionary 'fex_grtcoawt'
      dictionary 'fex_cth'
      record /fex_grttrans/ grttrans
      record /fex_grtcoawt/ grtcoawt
      record /fex_cth/ cth

      integer*4 etf_lun,num_recs,num_good

      dictionary 'fdq_sdf'
      dictionary 'fdq_eng'
      dictionary 'fic_scc'
      record /fdq_sdf/ sci_recs(num_recs)
      record /fdq_eng/ eng_recs(num_recs)
      record /fic_scc/ con_recs(num_recs)

      integer*4 fake_it,sci_mode,adds_per_group
      real*4 temps(10,num_recs),aifgs(512,num_recs)
      real*4 template(512),sec_template(512)
c
c     Input/output parameters.
c
      integer*4 num_sci_good(4,4,20),num_coadds(4,4)

      dictionary 'fic_sky'
      dictionary 'fic_cov'
      record /fic_sky/ coa_rec
      record /fic_cov/ cov_recs(771)

      integer*4 first_covar_times(2,3),last_covar_times(2,3)
c
c     Local variables.
c
      integer*4 pos,rec_good,rec
      integer*4 first_time(2),last_time(2)

      real*4 xcal_sum,glitch_rate_sum
      real*4 b_sum,c_sum,b_ave,c_ave

      integer*4 bins(fac_max_num),first_sci_rec,last_sci_rec
      integer*4 time_array(2,fac_max_num),time_weights(fac_max_num)
      integer*4 ifg_pos,status,int_val,int_most

      real*8 sum_equat(3),sum_terr(3),sum_elimb(3)
      real*8 sum_mnang(3),sum_gal(3),sum_ecl(3)

      real*4 terr_lat,terr_long,elimb,el_az
      real*4 moon_angle,moon_az
      real*4 gal_lat,gal_long,ecl_lat,ecl_long

      real*4 sum_sun_angle,moon_phase(fac_max_num)
      real*4 sum_sm_dist,sum_cm_dist
      real*4 sum_altitude,sum_baryvel,sum_lparam
      real*4 orb_phase(fac_max_num),sum_geovel
      real*4 scan_angle(fac_max_num),sum_sc_rot_angle

      integer*4 sum_att_sol(0:6)
      integer*4 sum_pixel_no(4),sum_pix_def(4)
      real*4 sum_sky_idx(4),sum_exc_gal_lat

      real*8 norm_equat,norm_terr,norm_elimb
      real*8 norm_mnang,norm_gal,norm_ecl
      real*4 avg_terr(3),avg_elimb(3),avg_mnang(3)
      real*4 avg_gal(3),avg_ecl(3),conv_ecl(3)

      real*4 ra,dec,avg_angle

      integer*4 temp_pos,eng_pos
      real*4 sum,val,sum_eng(6)
      real*8 dsum,dval
      logical*1 singlifg
      character*1 xcal_pos

      integer*4 speed,mode,etf_rec

      dictionary 'fex_etf'
      record /fex_etf/ etf

      character*72 filename
      integer*4 file_len
      integer*4 vec_pos,temp_scan_mode,check_time(2)

      complex*16 dis_vecs(265,fac_max_num),coa_vec(265)
      complex*16 scale(257)
      real*8 difg(512)
      real*4 divisor

      fic_coadd = %loc(fic_normal)
c
c     Accumulate good interferogram and number of coadd statistics.
c
      pos = (num_good - 1)/5 + 1
      num_sci_good(channel,scan_mode,pos) = num_sci_good(channel,scan_mode,pos) 
     .                                      + 1
      num_coadds(channel,scan_mode) = num_coadds(channel,scan_mode) + 1

      coa_rec.coad_spec_head.coadd_no = num_coadds(channel,scan_mode)
      coa_rec.coad_spec_head.num_ifgs = num_good

      call ct_gmt_to_binary(fac_jstop_default,first_time)
      call ct_gmt_to_binary(fac_jstart_default,last_time)
c
c     Initialize local variables.
c
      rec_good = 0
      xcal_sum = 0.0
      glitch_rate_sum = 0.0
      b_sum = 0.0
      c_sum = 0.0
      b_ave = 0.0
      c_ave = 0.0

      do rec = 1,num_recs
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
c     Accumulate gain, xcal position, and data and attitude quality.
c
            coa_rec.coad_spec_data.gain_sum = coa_rec.coad_spec_data.gain_sum +
     .              fac_gains(sci_recs(rec).sci_head.gain)
            xcal_sum = xcal_sum + sci_recs(rec).dq_data.xcal_pos
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
     .                                          aifgs(ifg_pos,rec)
            end do
c
c     Accumulate glitch rate and "b" and "c" coefficient values.
c
            glitch_rate_sum = glitch_rate_sum + con_recs(rec).glitch_rate
c
c     Use the glitch rate for binning and put it in the dispersion vectors.
c
            if (con_recs(rec).glitch_rate .ge. cth.glitch_rate(channel)) then
               bins(rec_good) = ibset(bins(rec_good),4)
            end if
            dis_vecs(265,rec_good) = dcmplx(con_recs(rec).glitch_rate)
            b_sum = b_sum + con_recs(rec).b
            c_sum = c_sum + con_recs(rec).c
c
c     Bin by galactic latitude, interferogram ct_head time, and secondary
c     template subtraction.
c
            if (abs(sci_recs(rec).attitude.galactic_latitude) .ge.
     .          cth.abs_galactic_lat) then
               bins(rec_good) = ibset(bins(rec_good),5)
            end if
            if (time_ge(sci_recs(rec).ct_head.time,cth.time_saa)) then
               bins(rec_good) = ibset(bins(rec_good),6)
            end if
            if (coa_rec.coad_spec_data.sec_template.subtracted .eq.
     .          fac_present) then
               bins(rec_good) = ibset(bins(rec_good),7)
            end if
         end if
      end do
c
c     Calculate the average time.
c
      if (single .ne. fac_present) then
         status = fut_ave_bin_times(time_array,num_good,time_weights,
     .            coa_rec.ct_head.time)
         call ct_binary_to_gmt(coa_rec.ct_head.time,coa_rec.ct_head.gmt)
c
c     Put the coadd time into the consistency check records.
c
         do rec = 1,num_recs
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
c     Coadd the interferogram, xcal position, and glitch rate.
c
      do ifg_pos=1,512
         coa_rec.coad_data.ifg(ifg_pos) =
     .           coa_rec.coad_data.ifg(ifg_pos)/num_good
      end do

      coa_rec.coad_spec_data.xcal_pos = nint(xcal_sum/num_good)
      coa_rec.coad_spec_data.glitch_rate = glitch_rate_sum/num_good
      coa_vec(265) = dcmplx(coa_rec.coad_spec_data.glitch_rate)
c
c     Calculate the averages and variances of the "b" and "c" coefficients.
c
      b_ave = b_sum/num_good
      coa_rec.coad_spec_data.sec_template.b_average = b_ave
      c_ave = c_sum/num_good
      coa_rec.coad_spec_data.transient.c_average = c_ave
      b_sum = 0.0
      c_sum = 0.0
      do rec = 1,num_recs
         if (con_recs(rec).con_check .eq. 0) then
            b_sum = b_sum + (con_recs(rec).b - b_ave)**2
            c_sum = c_sum + (con_recs(rec).c - c_ave)**2
         end if
      end do
      if (num_good .eq. 1) then
         coa_rec.coad_spec_data.sec_template.b_variance = -9999.0
         coa_rec.coad_spec_data.transient.c_variance = -9999.0
      else
         coa_rec.coad_spec_data.sec_template.b_variance = b_sum/(num_good-1)
         coa_rec.coad_spec_data.transient.c_variance = c_sum/(num_good-1)
      end if
c
c     If processing single interferograms, record the attitude directly.
c
      if (single .eq. fac_present) then
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
         sum_sm_dist = 0.0
         sum_cm_dist = 0.0
         sum_altitude = 0.0
         sum_baryvel = 0.0
         sum_lparam = 0.0
         sum_geovel = 0.0
         sum_sc_rot_angle = 0.0
         sum_exc_gal_lat = 0.0
         do pos = 0,6
            sum_att_sol(pos) = 0
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
         do rec = 1,num_recs
            if (con_recs(rec).con_check .eq. 0) then
               rec_good = rec_good + 1
c
c     Accumulate equatorial coordinates.
c
               do pos = 1,3
                  sum_equat(pos) = sum_equat(pos) + 
     .                             sci_recs(rec).attitude.equatorial(pos)
               end do
c
c     Accumulate terrestrial coordinates.
c
               terr_lat = sci_recs(rec).attitude.terr_latitude * fac_att_conv
               terr_long = sci_recs(rec).attitude.terr_longitude * fac_att_conv
               sum_terr(1) = sum_terr(1) + cosd(terr_lat) * cosd(terr_long)
               sum_terr(2) = sum_terr(2) + cosd(terr_lat) * sind(terr_long)
               sum_terr(3) = sum_terr(3) + sind(terr_lat)
c
c     Accumulate earth limb information.
c
               elimb = sci_recs(rec).attitude.earth_limb * fac_att_conv - 90.0
               el_az = sci_recs(rec).attitude.earth_limb_azimuth * fac_att_conv
               sum_elimb(1) = sum_elimb(1) + sind(elimb) * cosd(el_az)
               sum_elimb(2) = sum_elimb(2) + sind(elimb) * sind(el_az)
               sum_elimb(3) = sum_elimb(3) + cosd(elimb)
c
c     Accumulate sun and moon information.
c
               sum_sun_angle = sum_sun_angle + sci_recs(rec).attitude.sun_angle

               moon_angle = sci_recs(rec).attitude.moon_angle*fac_att_conv-90.0
               moon_az = sci_recs(rec).attitude.moon_az_angle * fac_att_conv
               sum_mnang(1) = sum_mnang(1) + sind(moon_angle) * cosd(moon_az)
               sum_mnang(2) = sum_mnang(2) + sind(moon_angle) * sind(moon_az)
               sum_mnang(3) = sum_mnang(3) + cosd(moon_angle)

               moon_phase(rec_good) = sci_recs(rec).attitude.moon_phase * 
     .                                fac_att_conv_rad

               sum_sm_dist = sum_sm_dist + sci_recs(rec).attitude.sun_moon_dist
               sum_cm_dist = sum_cm_dist + sci_recs(rec).attitude.cobe_moon_dist
c
c     Accumulate miscellaneous attitude parameters.
c
               sum_altitude = sum_altitude + sci_recs(rec).attitude.altitude
               sum_baryvel = sum_baryvel + 
     .                     sci_recs(rec).attitude.projected_barycentric_velocity
               sum_lparam = sum_lparam + sci_recs(rec).attitude.mcilwain_l_param
c
c     Accumulate galactic and ecliptic latitude and longitude.
c
               gal_lat = sci_recs(rec).attitude.galactic_latitude * fac_att_conv
               gal_long = sci_recs(rec).attitude.galactic_longitude*fac_att_conv
               sum_gal(1) = sum_gal(1) + cosd(gal_lat) * cosd(gal_long)
               sum_gal(2) = sum_gal(2) + cosd(gal_lat) * sind(gal_long)
               sum_gal(3) = sum_gal(3) + sind(gal_lat)

               ecl_lat = sci_recs(rec).attitude.ecliptic_latitude * fac_att_conv
               ecl_long = sci_recs(rec).attitude.ecliptic_longitude*fac_att_conv
               sum_ecl(1) = sum_ecl(1) + cosd(ecl_lat) * cosd(ecl_long)
               sum_ecl(2) = sum_ecl(2) + cosd(ecl_lat) * sind(ecl_long)
               sum_ecl(3) = sum_ecl(3) + sind(ecl_lat)
c
c     Accumulate miscellaneous attitude parameters.
c
               orb_phase(rec_good) = sci_recs(rec).attitude.orbital_phase * 
     .                               fac_att_conv_rad
               sum_geovel = sum_geovel + 
     .                      sci_recs(rec).attitude.projected_geocentric_velocity
               scan_angle(rec_good) = sci_recs(rec).attitude.scan_angle * 
     .                                fac_att_conv_rad
               sum_sc_rot_angle = sum_sc_rot_angle + 
     .                            sci_recs(rec).attitude.sc_rotation_angle
c
c     Accumulate attitude solution type.
c
               int_val = sci_recs(rec).attitude.solution
               if ((int_val .lt. 0) .or. (int_val .gt. 6)) then
                  int_val = 0
               end if
               sum_att_sol(int_val) = sum_att_sol(int_val) + 1
c
c     Accumulate pixelization information.
c
               if (sci_recs(rec).attitude.pixel_definition .eq. 'O') then
                  sum_pixel_no(2) = sum_pixel_no(2) + 
     .                              sci_recs(rec).attitude.pixel_no
                  sum_pix_def(2) = sum_pix_def(2) + 1
                  sum_sky_idx(2) = sum_sky_idx(2) + 
     .                             sci_recs(rec).attitude.skymap_index
               else if (sci_recs(rec).attitude.pixel_definition .eq. 'S') then
                  sum_pixel_no(3) = sum_pixel_no(3) + 
     .                              sci_recs(rec).attitude.pixel_no
                  sum_pix_def(3) = sum_pix_def(3) + 1
                  sum_sky_idx(3) = sum_sky_idx(3) + 
     .                             sci_recs(rec).attitude.skymap_index
               else if (sci_recs(rec).attitude.pixel_definition .eq. 'E') then
                  sum_pix_def(4) = sum_pix_def(4) + 1
                  sum_sky_idx(4) = sum_sky_idx(4) + 
     .                             sci_recs(rec).attitude.skymap_index
               else
                  sum_pixel_no(1) = sum_pixel_no(1) + 
     .                              sci_recs(rec).attitude.pixel_no
                  sum_pix_def(1) = sum_pix_def(1) + 1
                  sum_sky_idx(1) = sum_sky_idx(1) + 
     .                             sci_recs(rec).attitude.skymap_index
               end if
c
c     Accumulate miscellaneous attitude parameters.
c
               sum_exc_gal_lat = sum_exc_gal_lat +
     .                           sci_recs(rec).attitude.exc_galactic_lat
               coa_rec.attitude.terr_rad_byte = 
     .                 max(zext(coa_rec.attitude.terr_rad_byte),
     .                 zext(sci_recs(rec).attitude.terr_rad_byte))
            end if
         end do
c
c     Calculate coadded attitude if possible.
c
         if ((sum_pixel_no(1) .lt. 0) .or. ((sum_equat(1) .eq. 0.) .and. 
     .       (sum_equat(2) .eq. 0.) .and. (sum_equat(3) .eq. 0.))) then
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
            coa_rec.attitude.earth_limb =  (elimb + 90.0) / fac_att_conv
            coa_rec.attitude.earth_limb_azimuth = el_az / fac_att_conv
c
c     Coadd sun and moon information.
c
            coa_rec.attitude.sun_angle = sum_sun_angle / num_good

            moon_angle = atan2d(sqrt(avg_mnang(1)**2+avg_mnang(2)**2),
     .                          avg_mnang(3))
            moon_az = 0.0
            if (avg_mnang(1) .ne. 0.0) then
               moon_az = atan2d(avg_mnang(2),avg_mnang(1))
            end if
            coa_rec.attitude.moon_angle = (moon_angle + 90.0) / fac_att_conv
            coa_rec.attitude.moon_az_angle = moon_az / fac_att_conv

            status = fut_average_angles(moon_phase,num_good,avg_angle)
            coa_rec.attitude.moon_phase = nint(avg_angle / fac_att_conv_rad)

            coa_rec.attitude.sun_moon_dist = sum_sm_dist / num_good
            coa_rec.attitude.cobe_moon_dist = sum_cm_dist / num_good
c
c     Coadd miscellaneous attitude parameters.
c
            coa_rec.attitude.altitude = sum_altitude / num_good
            coa_rec.attitude.projected_barycentric_velocity = 
     .                    sum_baryvel / num_good
            coa_rec.attitude.mcilwain_l_param = sum_lparam / num_good
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
            status = fut_average_angles(orb_phase,num_good,avg_angle)
            coa_rec.attitude.orbital_phase = nint(avg_angle/fac_att_conv_rad)

            coa_rec.attitude.projected_geocentric_velocity = sum_geovel/num_good

            status = fut_average_angles(scan_angle,num_good,avg_angle)
            coa_rec.attitude.scan_angle = nint(avg_angle / fac_att_conv_rad)

            coa_rec.attitude.sc_rotation_angle = sum_sc_rot_angle / num_good
c
c     Coadd attitude solution type.
c
            int_most = -1
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
            int_most = -1
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
               if (coa_rec.attitude.skymap_index .lt. 6) then
                  coa_rec.attitude.pixel_no = sum_pixel_no(1)/sum_pix_def(1)
               end if
            else if (int_val .eq. 2) then
               coa_rec.attitude.pixel_definition = 'O'
               coa_rec.attitude.pixel_no = sum_pixel_no(2)/sum_pix_def(2)
               coa_rec.attitude.skymap_index = 
     .                 nint(sum_sky_idx(2)/sum_pix_def(2))
            else if (int_val .eq. 3) then
               coa_rec.attitude.pixel_definition = 'S'
               coa_rec.attitude.pixel_no = sum_pixel_no(3)/sum_pix_def(3)
               coa_rec.attitude.skymap_index = 
     .                 nint(sum_sky_idx(3)/sum_pix_def(3))
            else if (int_val .eq. 4) then
               coa_rec.attitude.pixel_definition = 'E'
               coa_rec.attitude.skymap_index = 
     .                       nint(sum_sky_idx(4)/sum_pix_def(4))
               if (coa_rec.attitude.skymap_index .lt. 6) then
                  coa_rec.attitude.pixel_no = sum_pixel_no(4)/sum_pix_def(4)
               end if
            end if
c
c     Coadd miscellaneous attitude parameters.
c
            coa_rec.attitude.exc_galactic_lat = sum_exc_gal_lat/num_good
         end if
      end if
c
c     If processing single interferograms, record engineering data and set
c     minima and maxima to flag values.
c
      if (single .eq. fac_present) then
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
         do rec = 1,num_recs
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
               dis_vecs(258,rec_good) = dcmplx(temps(1,rec))
               if (temps(1,rec) .ge. cth.xcal_temp) then
                  bins(rec_good) = ibset(bins(rec_good),0)
               end if
               dis_vecs(259,rec_good) = dcmplx(temps(2,rec))
               if (temps(2,rec) .ge. cth.ical_temp) then
                   bins(rec_good) = ibset(bins(rec_good),1)
               end if
               val = (temps(3,rec) + temps(4,rec)) / 2.0
               dis_vecs(260,rec_good) = dcmplx(val)
               if (val .ge. cth.skyhorn_refhorn_avg) then
                  bins(rec_good) = ibset(bins(rec_good),2)
               end if
               val = (temps(3,rec) - temps(4,rec)) / 2.0
               dis_vecs(261,rec_good) = dcmplx(val)
               dis_vecs(262,rec_good) = dcmplx(temps(5,rec))
               dis_vecs(263,rec_good) = dcmplx(temps(6,rec))
               dis_vecs(264,rec_good) = dcmplx(temps(6+channel,rec))
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
      if (single .eq. fac_present) then
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
            do rec = 1,num_recs
               if (con_recs(rec).con_check .eq. 0) then
                  sum = sum + eng_recs(rec).en_stat.group1(eng_pos)
               end if
            end do
            coa_rec.en_stat.group1(eng_pos) = nint(sum / num_good)
         end do

         do eng_pos=1,14
            sum = 0.0
            do rec = 1,num_recs
               if (con_recs(rec).con_check .eq. 0) then
                  sum = sum + eng_recs(rec).en_stat.group2(eng_pos)
               end if
            end do
            coa_rec.en_stat.group2(eng_pos) = nint(sum / num_good)
         end do
c
c     Coadd hot spot command and power statuses.
c
         do eng_pos = 1,6
            sum_eng(eng_pos) = 0.0
         end do
         do rec = 1,num_recs
            if (con_recs(rec).con_check .eq. 0) then
               sum_eng(1) = sum_eng(1) + eng_recs(rec).en_stat.hot_spot_cmd(1)
               sum_eng(2) = sum_eng(2) + eng_recs(rec).en_stat.hot_spot_cmd(2)
               sum_eng(3) = sum_eng(3) + eng_recs(rec).en_stat.power_a_status(1)
               sum_eng(4) = sum_eng(4) + eng_recs(rec).en_stat.power_a_status(2)
               sum_eng(5) = sum_eng(5) + eng_recs(rec).en_stat.power_b_status(1)
               sum_eng(6) = sum_eng(6) + eng_recs(rec).en_stat.power_b_status(2)
            end if
         end do
         coa_rec.en_stat.hot_spot_cmd(1) = nint(sum_eng(1) / num_good)
         coa_rec.en_stat.hot_spot_cmd(2) = nint(sum_eng(2) / num_good)
         coa_rec.en_stat.power_a_status(1) = nint(sum_eng(3) / num_good)
         coa_rec.en_stat.power_a_status(2) = nint(sum_eng(4) / num_good)
         coa_rec.en_stat.power_b_status(1) = nint(sum_eng(5) / num_good)
         coa_rec.en_stat.power_b_status(2) = nint(sum_eng(6) / num_good)
      end if

      coa_rec.coad_spec_data.bol_cmd_bias=coa_rec.en_stat.bol_cmd_bias(channel)
c
c     Coadd grt temperatures and remaining engineering quantities and 
c     calculate sigmas; record them and set flag values if processing single 
c     interferograms.
c
      if (single .eq. fac_present) then
         do eng_pos = 1,64
            coa_rec.en_analog.grt(eng_pos) = eng_recs(1).en_analog.grt(eng_pos)
         end do
         do eng_pos = 1,62
            coa_rec.en_analog.group1(eng_pos) = 
     .              eng_recs(1).en_analog.group1(eng_pos)
         end do
         do eng_pos = 1,64
            coa_rec.en_sigma.sig_grt(eng_pos) = -9999.0
            coa_rec.en_sigma.group2(eng_pos) = -9999.0
         end do
      else
c
c     Coadd grt temperatures.
c
         do eng_pos = 1,64
            rec_good = 0
            sum = 0.0
            do rec = 1,num_recs
               if (con_recs(rec).con_check .eq. 0) then
                  val = eng_recs(rec).en_analog.grt(eng_pos)
                  if (val .gt. 0.0) then
                     rec_good = rec_good + 1
                     sum = sum + val
                  end if
               end if
            end do
            if (rec_good .ne. 0) then
               coa_rec.en_analog.grt(eng_pos) = sum/rec_good
               dsum = 0.0
               do rec = 1,num_recs
                  if (con_recs(rec).con_check .eq. 0) then
                     dval = eng_recs(rec).en_analog.grt(eng_pos)
                     if (dval .gt. 0.0) then
                        dsum = dsum + (dval - coa_rec.en_analog.grt(eng_pos))**2
                        coa_rec.en_sigma.sig_grt(eng_pos) = sqrt(dsum/rec_good)
                     end if
                  end if
               end do
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
            do rec = 1,num_recs
               if (con_recs(rec).con_check .eq. 0) then
                  val = eng_recs(rec).en_analog.group1(eng_pos)
                  if (val .gt. 0.0) then
                     rec_good = rec_good + 1
                     sum = sum + val
                  end if
               end if
            end do
            if (rec_good .ne. 0) then
               coa_rec.en_analog.group1(eng_pos) = sum/rec_good
               dsum = 0.0
               do rec = 1,num_recs
                  if (con_recs(rec).con_check .eq. 0) then
                     dval = eng_recs(rec).en_analog.group1(eng_pos)
                     if (dval .gt. 0.0) then
                        dsum = dsum+(dval-coa_rec.en_analog.group1(eng_pos))**2
                        coa_rec.en_sigma.group2(eng_pos) = sqrt(dsum/rec_good)
                     end if
                  end if
               end do
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
      if (single .eq. fac_present) then
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
     .           fac_scan_mode_ids(scan_mode),
     .           nint(coa_rec.coad_spec_data.ical*1000.0),
     .           xcal_pos,coa_rec.attitude.pixel_no
10       format (2(a,'_'),i3.3,'_',a,'_',i5.5,'_',a,'_',i4.4,11x)
      else
         encode (60,20,coa_rec.coad_spec_head.label)
     .           coa_rec.coad_spec_head.first_gmt,
     .           coa_rec.coad_spec_head.last_gmt,
     .           coa_rec.coad_spec_head.num_ifgs,
     .           fac_scan_mode_ids(scan_mode),
     .           nint(coa_rec.coad_spec_data.ical*1000.0),
     .           xcal_pos
20       format (2(a,'_'),i3.3,'_',a,'_',i5.5,'_',a,16x)
      end if
c
c     If processing single interferograms, set the three statistical vectors
c     in the coadd record to flag values.
c
      if (single .eq. fac_present) then
         do temp_pos = 1,257
            coa_rec.coad_data.real_var(temp_pos) = -9999.0
            coa_rec.coad_data.imag_var(temp_pos) = -9999.0
            coa_rec.coad_data.real_imag_var(temp_pos) = -9999.0
         end do
      else
c
c     When processing a coadd group, add the engineering data to the end
c     of the coadd dispersion vector.
c
         coa_vec(258) = dcmplx(coa_rec.coad_spec_data.xcal)
         coa_vec(259) = dcmplx(coa_rec.coad_spec_data.ical)
         coa_vec(260) = dcmplx((coa_rec.coad_spec_data.skyhorn +
     .                   coa_rec.coad_spec_data.refhorn) / 2.0)
         coa_vec(261) = dcmplx((coa_rec.coad_spec_data.skyhorn - 
     .                   coa_rec.coad_spec_data.refhorn) / 2.0)
         coa_vec(262) = dcmplx(coa_rec.coad_spec_data.dihedral)
         coa_vec(263) = dcmplx(coa_rec.coad_spec_data.collimator_mirror)
         coa_vec(264) = dcmplx(coa_rec.coad_spec_data.bol_assem(channel))
c
c     Determine the correct record for the electronics transfer function.
c
         if (fake_it .eq. 1) then
            speed = 2
         else
            speed = coa_rec.coad_spec_data.mtm_speed
         end if
         
         if ((sci_mode .eq. 2) .or. (sci_mode .eq. 4)) then
            mode = 0
         else
            mode = 1
         end if

         etf_rec = 96*speed + 24*(channel-1) + 12*mode + adds_per_group
c
c     Read the electronics transfer function.
c
         read(etf_lun,rec=etf_rec,iostat=status) etf
         if (status .ne. 0) then
            fic_coadd = %loc(fic_cfdirread)
            inquire(etf_lun,name=filename)
            call str$trim(filename,filename,file_len)
            call lib$signal(fic_cfdirread,%val(2),filename(:file_len),
     .                      %val(status))
            return
         end if
c
c     Construct the scaling vector from the electronics transfer function, the
c     etendu, the adc scale, and the nyquist frequency in icm.
c
         do vec_pos = 1,257
            if (etf.ztransfcn(vec_pos) .eq. 0.0) then
               scale(vec_pos) = 1.0 / (fac_etendu * fac_adc_scale * 
     .                          coa_rec.coad_spec_data.nyquist_icm)
            else
               scale(vec_pos) = 1.0 / (etf.ztransfcn(vec_pos) * 
     .                          fac_etendu * fac_adc_scale * 
     .                          coa_rec.coad_spec_data.nyquist_icm)
            end if
         end do
c
c     Fast Fourier transform the coadded interferogram.
c
         do ifg_pos = 1,512
            difg(ifg_pos) = coa_rec.coad_data.ifg(ifg_pos)
         end do
         call dfftrf(512,difg,difg)
c
c     Set the first point of the coadd dispersion vector to zero.
c
         coa_vec(1) = dcmplx(0.0,0.0)
c
c     Unpack the Fourier transform and multiply by the scaling vector.
c
         do vec_pos = 2,256
            coa_vec(vec_pos) = dcmplx(difg(vec_pos*2-2),difg(vec_pos*2-1)) * 
     .                         scale(vec_pos)
         end do
         coa_vec(257) = dcmplx(difg(512)) * scale(257)
c
c     Calculate a divisor value for the statistics based on the template
c     subtraction information.
c         
         if (coa_rec.coad_spec_data.sec_template.subtracted .eq. 
     .       fac_present) then
            divisor = num_good * (num_good - 2)
         else if (coa_rec.coad_spec_data.prim_template.subtracted .eq.
     .            fac_present) then
            divisor = num_good * (num_good - 1)
         else
            divisor = num_good
         end if

         rec_good = 0
c
c     Loop over the good interferograms to form the dispersion vectors.
c
         do rec = 1,num_recs
            if (con_recs(rec).con_check .eq. 0) then
               rec_good = rec_good + 1
c
c     Fast Fourier transform the interferogram.
c
               do ifg_pos = 1,512
                  difg(ifg_pos) = aifgs(ifg_pos,rec)
               end do
               call dfftrf(512,difg,difg)
c
c     Set the first point of the dispersion vector to zero.
c
               dis_vecs(1,rec_good) = dcmplx(0.0,0.0)
c
c     Unpack the Fourier transform, multiply by the scaling vector, and
c     subtract the coadd dispersion vector.
c
               do vec_pos = 2,256
                  dis_vecs(vec_pos,rec_good) = 
     .                     dcmplx(difg(vec_pos*2-2),difg(vec_pos*2-1)) *
     .                     scale(vec_pos) - coa_vec(vec_pos)
               end do
               dis_vecs(257,rec_good) = dcmplx(difg(512)) * scale(257) -
     .                                  coa_vec(257)
               do vec_pos = 258,265
                  dis_vecs(vec_pos,rec_good) = dis_vecs(vec_pos,rec_good) -
     .                                         coa_vec(vec_pos)
               end do
c
c     Calculte the three statistical vectors:  real-real variance, 
c     imaginary-imaginary variance, and real-imaginary variance; leave them
c     set at zero if the divisor value is zero.
c
               if (divisor .ne. 0.0) then
                  do vec_pos = 1,257
                     coa_rec.coad_data.real_var(vec_pos) =
     .                       coa_rec.coad_data.real_var(vec_pos) + 
     .                       dreal(dis_vecs(vec_pos,rec_good))**2 / divisor
                     coa_rec.coad_data.imag_var(vec_pos) =
     .                       coa_rec.coad_data.imag_var(vec_pos) + 
     .                       dimag(dis_vecs(vec_pos,rec_good))**2 / divisor
                     coa_rec.coad_data.real_imag_var(vec_pos) =
     .                       coa_rec.coad_data.real_imag_var(vec_pos) + 
     .                       (dreal(dis_vecs(vec_pos,rec_good)) *
     .                       dimag(dis_vecs(vec_pos,rec_good))) / divisor
                  end do
               end if
            end if
         end do
c
c     Temporarily reset the scan mode value if it is 4.
c
         if ((covar .eq. fac_present) .and. (scan_mode .ne. 3)) then
            temp_scan_mode = scan_mode
            if (temp_scan_mode .eq. 4) then
               temp_scan_mode = 3
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
            status = fic_covar(temp_scan_mode,num_good,dis_vecs,sci_mode,
     .                         adds_per_group,coa_vec,bins,cov_recs)
         end if
c
c     If the primary or secondary templates were subtracted, add them back to
c     the coadded interferogram.
c
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
