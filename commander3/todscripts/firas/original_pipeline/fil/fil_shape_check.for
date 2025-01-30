      integer*4 function fil_shape_check(channel,scan_mode,cth,num_recs,
     .                                   num_cgr_recs,num_good,num_cgr_good,
     .                                   min_good,num_sci_fail,con_recs,aifgs,
     .                                   avg_noise)
c------------------------------------------------------------------------------
c
c     Purpose: Check interferograms for consistency in shape.
c
c     Author: R. Isaacman, ARC, 2/86
c
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode          i*4  Value of scan mode, 1-4 = SS-LF.
c            cth                rec  Consistency check parameters.
c            num_recs           i*4  Number of input science and engineering
c                                    records, including neighbors.
c            num_cgr_recs       i*4  Number of input science and engineering
c                                    records, excluding neighbors.
c            num_good           i*4  Current number of good interferograms,
c                                    including neighbors.
c            num_cgr_good       i*4  Current number of good interferograms,
c                                    excluding neighbors.
c            num_sci_fail       i*4(4,4,32)  Number of science or engineering
c                                            records failing by channel and 
c                                            scan mode for each of 32 
c                                            consistency checks.
c            con_recs           rec(num_recs)  Consistency check records.
c            aifgs              r*4(512,num_recs)  Real-valued interferograms.
c            avg_noise          r*4  Average noise of coadd group.
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
c             con_recs           rec(num_recs)  Consistency check records.
c
c     Modifications:
c
c------------------------------------------------------------------------------

      implicit none
c
c     Return statuses.
c
      external fil_normal
c
c     Input parameters.
c
      integer*4 channel,scan_mode,num_recs,num_cgr_recs

      dictionary 'fex_cth'
      record /fex_cth/ cth

      real*4 aifgs(512,num_recs),avg_noise
c
c     Input/output parameters.
c
      integer*4 num_good,num_cgr_good,num_sci_fail(4,4,32)

      dictionary 'fil_scc'
      record /fil_scc/ con_recs(num_recs)
c
c     Output parameters.
c
      integer*4 min_good
c
c     Local variables.
c
      integer*4 freq,rec,num_dev,ifg_pos
      real*4 noise,dev

      fil_shape_check = %loc(fil_normal)
c
c     Determine number of good interferograms required to coadd.
c
      min_good = cth.min_ifg_shape_coadd(channel)
c
c     Set frequency to high or low.
c
      if ((channel .eq. 1) .or. (channel .eq. 3)) then
         freq = 1
      else
         freq = 2
      end if
c
c     Loop through good interferograms.
c
      rec = 1
      do while ((num_good .ge. min_good) .and. (num_cgr_good .ge. 1) .and.
     .          (rec .le. num_recs))
         if (con_recs(rec).con_check .eq. 0) then
c
c     Find the ratio of the interferogram noise to the coadd group noise.
c
            noise = con_recs(rec).noise/avg_noise
c
c     Check that the ratio is within tolerance specified in the cth reference
c     data set.  If not, fail the record and decrement the number of good
c     interferograms remaining.
c
            if (noise .lt. cth.min_ifg_noise(channel)) then
               con_recs(rec).con_check = ibset(con_recs(rec).con_check,28)
               num_good = num_good - 1
               if (rec .le. num_cgr_recs) then
                  num_cgr_good = num_cgr_good - 1
                  num_sci_fail(channel,scan_mode,29) = 
     .                         num_sci_fail(channel,scan_mode,29) + 1
               end if
            else if (noise .gt. cth.max_ifg_noise(channel)) then
               con_recs(rec).con_check = ibset(con_recs(rec).con_check,29)
               num_good = num_good - 1
               if (rec .le. num_cgr_recs) then
                  num_cgr_good = num_cgr_good - 1
                  num_sci_fail(channel,scan_mode,30) = 
     .                         num_sci_fail(channel,scan_mode,30) + 1
               end if
            else
c
c     Determine at each point of the interferogram, from a cth mask value to 
c     512, the absolute value of the ratio of the point to the coadd group
c     noise.  If this ratio is larger than a maximum in the cth file, count
c     it as a bad point.  If the number of bad points exceed another cth
c     maximum, mark the record as bad.
c
               num_dev = 0

               do ifg_pos = cth.mask_ifg(freq,scan_mode),512
                  dev = abs(aifgs(ifg_pos,rec)/avg_noise)
                  if (dev .gt. cth.max_point_deviation(channel)) then
                     num_dev = num_dev + 1
                  end if
               end do

               if (num_dev .gt. cth.max_bad_points(channel)) then
                  con_recs(rec).con_check = ibset(con_recs(rec).con_check,30)
                  num_good = num_good - 1
                  if (rec .le. num_cgr_recs) then
                     num_cgr_good = num_cgr_good - 1
                     num_sci_fail(channel,scan_mode,31) = 
     .                            num_sci_fail(channel,scan_mode,31) + 1
                  end if
               end if
            end if
         end if
         rec = rec + 1
      end do
c
c     If too few good interferograms remain to coadd, mark the remaining 
c     records.
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
