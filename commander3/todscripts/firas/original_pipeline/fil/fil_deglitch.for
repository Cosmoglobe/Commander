      integer*4 function fil_deglitch(channel,scan_mode,mincoadd,gltchpro_lun,
     .                                num_recs,num_cgr_recs,num_good,
     .                                num_cgr_good,min_good,num_sci_deglitch,
     .                                num_sci_fail,fake_it,adds_per_group,
     .                                con_recs,coa_rec,aifgs,glitch_pos,
     .                                glitch_amp,glitch_tot_iter,
     .                                glitch_tot_signal)
c-------------------------------------------------------------------------------
c
c     Purpose: Deglitch interferograms using the "clean" algorithm from radio 
c              astronomy: search for peaks well above the noise, and subtract 
c              some fraction of a glitch profile whenever one is found. The 
c              glitch profile used is found by calculating the "true" position 
c              of the peak and seeing which profile has its peak closest to 
c              the calculated position.
c
c     Author: R. Isaacman, GSC, 1/91
c  
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode          i*4  Value of scan mode, 1-4 = SS-LF.
c            mincoadd           rec  Minimum interferograms to coadd.
c            gltchpro_lun       i*4  Logical unit for glitch profiles.
c            num_recs           i*4  Number of input science and engineering
c                                    records, including neighbors.
c            num_cgr_recs       i*4  Number of input science and engineering
c                                    records, excluding neighbors.
c            num_good           i*4  Current number of good interferograms,
c                                    including neighbors.
c            num_cgr_good       i*4  Current number of good interferograms,
c                                    excluding neighbors.
c            num_sci_deglitch   i*4(4,4)  Number of deglitched interferograms
c                                         by channel and scan mode.
c            num_sci_fail       i*4(4,4,32)  Number of science or engineering
c                                            records failing by channel and 
c                                            scan mode for each of 32 
c                                            consistency checks.
c            fake_it            i*4  Value of fake-it bit, 0-1.
c            adds_per_group     i*4  Value of adds per group, 1-12.
c            con_recs           rec(num_recs)  Consistency check records.
c            coa_rec            rec  Coadd record.
c            aifgs              r*4(512,num_recs)  Real-valued interferograms.
c            glitch_tot_iter    i*4(4,4)  Total deglitcher iterations by
c                                         channel and scan mode.
c            glitch_tot_signal  r*4(4,4)  Total deglitcher signal by channel 
c                                         and scan mode.
c
c     Output: num_good           i*4  Current number of good interferograms,
c                                     including neighbors.
c             num_cgr_good       i*4  Current number of good interferograms,
c                                     excluding neighbors.
c             min_good           i*4  Minimum number of good interferograms
c                                     required to coadd.
c             num_sci_deglitch   i*4(4,4)  Number of deglitched interferograms
c                                          by channel and scan mode.
c             num_sci_fail       i*4(4,4,32)  Number of science or engineering
c                                             records failing by channel and 
c                                             scan mode for each of 32 
c                                             consistency checks.
c             con_recs           rec(num_recs)  Consistency check records.
c             coa_rec            rec  Coadd record.
c             aifgs              r*4(512,num_recs)  Real-valued interferograms.
c             glitch_pos         i*4(1000,num_recs)  Glitch positions.
c             glitch_amp         r*4(1000,num_recs)  Glitch amplitudes.
c             glitch_tot_iter    i*4(4,4)  Total deglitcher iterations by
c                                          channel and scan mode.
c             glitch_tot_signal  r*4(4,4)  Total deglitcher signal by channel 
c                                          and scan mode.
c
c     Modifications:
c
c------------------------------------------------------------------------

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
      integer*4 fil_find_glitch
c
c     Input parameters.
c
      integer*4 channel,scan_mode,gltchpro_lun

      dictionary 'fex_mincoadd'
      record /fex_mincoadd/ mincoadd

      integer*4 num_recs,num_cgr_recs,fake_it,adds_per_group
c
c     Input/output parameters.
c
      integer*4 num_good,num_cgr_good,num_sci_deglitch(4,4),num_sci_fail(4,4,32)

      dictionary 'fil_scc'
      dictionary 'fil_sky'

      record /fil_scc/ con_recs(num_recs)
      record /fil_sky/ coa_rec

      real*4 aifgs(512,num_recs)
      integer*4 glitch_tot_iter(4,4)
      real*4 glitch_tot_signal(4,4)
c
c     Output parameters.
c
      integer*4 min_good,glitch_pos(fac_max_deglitch,num_recs)
      real*4 glitch_amp(fac_max_deglitch,num_recs)
c
c     Local variables.
c
      integer*4 speed,skip,start,adds_rec,status

      dictionary 'fex_gltchpro'
      record /fex_gltchpro/ gltchpro(12)

      character*72 filename
      integer*4 file_len
      integer*4 ifg_pos,iter,gl_pos
      real*4 sort(512),noise,median
      real*4 aifg(512),gl_amp
      real*4 sum,offset,delta,min_delta
      integer*4 rec,rec_best,act_pos
      integer*4 first,last,index

      fil_deglitch = %loc(fil_normal)

      min_good = mincoadd.min_ifg_coadd(channel)
c
c     Determine correct glitch profiles to retrieve based on fake-it status,
c     mtm speed, and adds per group.
c
      if (fake_it .eq. 1) then
         speed = 2
      else
         speed = coa_rec.coad_spec_data.mtm_speed
      end if

      if (adds_per_group .eq. 1) then
         skip = 0
      else if (adds_per_group .eq. 2) then
         skip = 1
      else if (adds_per_group .eq. 3) then
         skip = 3
      else if (adds_per_group .eq. 8) then
         skip = 6
      else if (adds_per_group .eq. 12) then
         skip = 14
      end if

      start = 104*speed + 26*(channel-1) + skip + 1
c
c     Read glitch profiles, one for each values of adds per group.
c
      do adds_rec = 1,adds_per_group
         read(gltchpro_lun,rec=start+adds_rec-1,iostat=status) 
     .        gltchpro(adds_rec)
         if (status .ne. 0) then
            fil_deglitch = %loc(fil_cfdirread)
            inquire(gltchpro_lun,name=filename)
            call str$trim(filename,filename,file_len)
            call lib$signal(fil_cfdirread,%val(2),filename(:file_len),
     .                      %val(status))
            return
         end if
      end do
c
c     For each good interferogram, calculate the noise as 1.25 * median 
c     absolute value in order to scale the deglitching process.
c
      rec = 1
      do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .          (rec .le. num_recs))
         if (con_recs(rec).con_check .eq. 0) then
            do ifg_pos = 1,512
               aifg(ifg_pos) = aifgs(ifg_pos,rec)
            end do
            call svrgn(512,aifg,sort)
            median = (sort(256) + sort(257))/2.0
            con_recs(rec).pre_dg_median = median
            call svrbn(512,aifg,sort)
            noise = 1.25 * ((abs(sort(256)) + abs(sort(257))) / 2.0)
c
c     Scale the noise according to whether the secondary template has been
c     subtracted or not.
c
            if (coa_rec.coad_spec_data.sec_template.subtracted .eq. 
     .          fac_present) then
               noise = (noise * num_good) / (num_good - 1 - 
     .                                       min(abs(con_recs(rec).b),1.0))
            else
               noise = (noise * num_good) / (num_good - 1)
            end if
            con_recs(rec).pre_dg_noise = noise
c
c     Set the noise to the maximum of the interferogram noise and the 
c     threshold bit noise from fil_qual_check.
c
            noise = max(noise,con_recs(rec).noise)
c
c     Find the first peak (potential glitch) in the interferogram to initialize
c     the deglitching loop. Every time fil_find_glitch returns a peak with 
c     significant height, subtract the glitch profile centered on that peaks
c     position. Continue the process as long as fil_find_glitch returns peaks 
c     with nonzero scaled heights, or until a maximum number of iterations is 
c     reached.
c
            iter = 1

            status = fil_find_glitch(aifg,noise,gl_pos,gl_amp)

            if (gl_amp .ne. 0.0) then
               glitch_pos(iter,rec) = gl_pos
               glitch_amp(iter,rec) = gl_amp
               con_recs(rec).glitch_pos = gl_pos
               con_recs(rec).glitch_amp = gl_amp
               sum = gl_amp
            else
               sum = 0.0
            end if
            do while ((gl_amp. ne. 0.0) .and. (iter .lt. fac_max_deglitch))
c
c     Find the offset to the "true" peak position by fitting a parabola
c     to the three points defining the peak of the glitch.  If the glitch
c     position is 1 or 512, arbitrarily assign the offset to zero. 
c
               if ((gl_pos .ne. 1) .and. (gl_pos .ne. 512)) then
                  offset = (aifgs(gl_pos-1,rec) - aifgs(gl_pos+1,rec)) / 
     .                     (aifgs(gl_pos-1,rec) - 2.0*aifgs(gl_pos,rec) + 
     .                      aifgs(gl_pos+1,rec))/2.0
               else
                  offset = 0.0
               end if
c
c     Step through the glitch profiles and find the one with its peak closest 
c     to that of the "true" peak. 
c
               min_delta = 1000.
               do adds_rec = 1,adds_per_group
                  delta = abs(gltchpro(adds_rec).gltchpro(512) - offset) 
                  if (delta .lt. min_delta) then
                     min_delta = delta
                     rec_best = adds_rec
                  end if
               end do
               act_pos = gltchpro(rec_best).gltchpro(511)
c
c     Deglitch with the profile found above.
c
               first = max(1,gl_pos-act_pos+1)
               last = min(512,gl_pos+510-act_pos)
               do ifg_pos = first,last
                  index = ifg_pos - gl_pos + act_pos
                  aifgs(ifg_pos,rec) = aifgs(ifg_pos,rec) -
     .                                 gl_amp*gltchpro(rec_best).gltchpro(index)
               end do
               do ifg_pos = 1,512
                  aifg(ifg_pos) = aifgs(ifg_pos,rec)
               end do

               status = fil_find_glitch(aifg,noise,gl_pos,gl_amp)

               iter = iter + 1
               if (gl_amp .ne. 0.0) then
                  glitch_pos(iter,rec) = gl_pos
                  glitch_amp(iter,rec) = gl_amp
                  sum = sum + gl_amp
               end if
            end do
c
c     If the maximum number of deglitcher iterations was reached, mark the
c     interferogram as bad in the consistency check record and decrement the
c     number of good interferograms.
c
            if (iter .eq. fac_max_deglitch) then
               con_recs(rec).glitch_iter = iter
               con_recs(rec).glitch_signal = sum
               con_recs(rec).con_check = ibset(con_recs(rec).con_check,27)
               num_good = num_good - 1
               if (rec .le. num_cgr_recs) then
                  num_cgr_good = num_cgr_good - 1
                  num_sci_fail(channel,scan_mode,28) = 
     .                         num_sci_fail(channel,scan_mode,28) + 1
               end if
            else
c
c     Accumulate the deglitcher iterations and signal removed.
c
               con_recs(rec).glitch_iter = iter - 1
               con_recs(rec).glitch_signal = sum
               if (rec .le. num_cgr_recs) then
                  coa_rec.coad_spec_data.deglitch.glitch_iter = 
     .                    coa_rec.coad_spec_data.deglitch.glitch_iter + iter - 1
                  coa_rec.coad_spec_data.deglitch.glitch_signal = 
     .                    coa_rec.coad_spec_data.deglitch.glitch_signal + sum
               end if
            end if
         end if
         rec = rec + 1
      end do
c
c     Update the deglitcher statistics arrays.
c
      if (num_cgr_good .gt. 0) then
         num_sci_deglitch(channel,scan_mode)=num_sci_deglitch(channel,scan_mode)
     .                                         + num_cgr_good
         glitch_tot_iter(channel,scan_mode)=glitch_tot_iter(channel,scan_mode) +
     .                   coa_rec.coad_spec_data.deglitch.glitch_iter
         glitch_tot_signal(channel,scan_mode) = 
     .                     glitch_tot_signal(channel,scan_mode) + 
     .                     coa_rec.coad_spec_data.deglitch.glitch_signal
c
c     Store deglitcher information in the coadd record.
c
         coa_rec.coad_spec_data.deglitch.glitch_iter = 
     .           coa_rec.coad_spec_data.deglitch.glitch_iter/num_cgr_good
         coa_rec.coad_spec_data.deglitch.glitch_signal = 
     .           coa_rec.coad_spec_data.deglitch.glitch_signal/num_cgr_good
      end if
c
c     If too few interferograms remain in the group to coadd, mark the
c     remaining good interferograms as failing.
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
