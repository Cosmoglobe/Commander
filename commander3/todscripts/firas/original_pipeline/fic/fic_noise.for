      integer*4 function fic_noise(channel,scan_mode,num_recs,num_good,
     .                             num_sci_noise,adds_per_group,con_recs,
     .                             coa_rec,aifgs,avg_noise,tot_noise)
c-------------------------------------------------------------------------------
c
c     Purpose: Calculate the noise amplitude for each interferogram after
c              template subtraction and deglitching, and calculate the noise
c              of the coadd group.
c
c     Author: S. Read, STX, 6/92, SER 7985
c
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode          i*4  Value of scan mode, 1-4=SS-LF.
c            num_recs           i*4  Number of input science and engineering
c                                    records.
c            num_good           i*4  Current number of good interferograms.
c            num_sci_noise      i*4(4,4)  Number of interferograms for which a
c                                         noise was calculated by channel and 
c                                         scan mode.
c            adds_per_group     i*4  Value of adds per group, 1-12.
c            con_recs           rec(num_recs)  Consistency check records.
c            coa_rec            rec  Coadd record.
c            aifgs              r*4(512,num_recs)  Real-valued interferograms.
c            tot_noise          r*4(4,4)  Total noise by channel and scan mode.
c
c     Output: num_sci_noise      i*4(4,4)  Number of interferograms for which a
c                                          noise was calculated by channel and 
c                                          scan mode.
c             con_recs           rec(num_recs)  Consistency check records.
c             coa_rec            rec  Coadd record.
c             avg_noise          r*4  Average noise of coadd group.
c             tot_noise          r*4(4,4)  Total noise by channel and scan mode.
c
c     Modifications:
c
c-------------------------------------------------------------------------------

      implicit none
c
c     Return statuses.
c
      external fic_normal
c
c     Input parameters.
c
      integer*4 channel,scan_mode,num_recs,num_good,adds_per_group

      real*4 aifgs(512,num_recs)
c
c     Input/output parameters.
c
      integer*4 num_sci_noise(4,4)

      dictionary 'fic_scc'
      record /fic_scc/ con_recs(num_recs)

      dictionary 'fic_sky'
      record /fic_sky/ coa_rec

      real*4 tot_noise(4,4)
c
c     Output parameters.
c
      real*4 avg_noise
c
c     Local variables.
c
      integer*4 rec,ifg_pos,rec_good
      real*4 sort(512),median,sum_gain

      fic_noise = %loc(fic_normal)
c
c     Update the noise statistics array.
c
      num_sci_noise(channel,scan_mode) = num_sci_noise(channel,scan_mode) +
     .                                   num_good
c
c     Calculate the noise for each interferogram as 1.25 times the median of
c     the absolute values of the interferogram minus the median of the
c     interferogram.
c
      do rec = 1,num_recs
         if (con_recs(rec).con_check .eq. 0) then
            do ifg_pos = 1,512
               sort(ifg_pos) = aifgs(ifg_pos,rec)
            end do
            call svrgn(512,sort,sort)
            median = (sort(256) + sort(257))/2.0
            do ifg_pos = 1,512
               sort(ifg_pos) = abs(aifgs(ifg_pos,rec) - median)
            end do
            call svrgn(512,sort,sort)
            con_recs(rec).noise = 1.25*(sort(256) + sort(257))/2.0
            tot_noise(channel,scan_mode) = tot_noise(channel,scan_mode) +
     .                                     con_recs(rec).noise
         end if
      end do
c
c     Calculate the noise of the coadd group as the median of the noises of
c     the individual interferograms.  Accumulate the reciprocals of the gain
c     values squared for threshold bit noise calculation.
c
      sum_gain = 0.0

      rec_good = 0
      do rec = 1,num_recs
         if (con_recs(rec).con_check .eq. 0) then
            rec_good = rec_good + 1
            sort(rec_good) = con_recs(rec).noise
            sum_gain = sum_gain + (1.0 / con_recs(rec).gain**2)
         end if
      end do
      call svrgn(num_good,sort,sort)

      if (mod(num_good,2) .eq. 1) then
         median = sort(num_good/2 + 1)
      else
         median = (sort(num_good/2) + sort(num_good/2 + 1))/2.0
      end if
c
c     Determine the threshold bit noise of the coadd group.
c
      avg_noise = sqrt(sum_gain / (12.0 * num_good * adds_per_group))
c
c     Take the maximum of the computed median noise and the threshold bit 
c     noise as the coadd group noise.
c
      avg_noise = max(avg_noise,median)

      coa_rec.coad_spec_data.noise = avg_noise

      return
      end
