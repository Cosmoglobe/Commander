      integer*4 function fil_convert(channel,num_recs,sci_recs,sweeps,con_recs,
     .                               aifgs)
c-------------------------------------------------------------------------------
c
c     Purpose: To convert the interferograms to real numbers, divide by the
c              real-valued gain, divide by the number of sweeps, and calculate
c              and subtract the dither.
c
c     Author: S. Brodd, HSTX, 4/95
c
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            num_recs           i*4  Number of input science and engineering
c                                    records, including neighbors.
c            sci_recs           rec(num_recs)  Science records.
c            sweeps             i*4  Value of sweeps, 1-16.
c            con_recs           rec(num_recs)  Consistency check records.
c
c     Output: aifgs              r*4(512,num_recs)  Real-valued interferograms.
c
c     Modifications:
c
c-----------------------------------------------------------------------------

      implicit none
c
c     Return statuses.
c
      external fil_normal
c
c     Input parameters.
c
      integer*4 channel,num_recs,sweeps

      dictionary 'fdq_sdf'
      dictionary 'fil_scc'

      record /fdq_sdf/ sci_recs(num_recs)
      record /fil_scc/ con_recs(num_recs)
c
c     Output parameters.
c
      real*4 aifgs(512,num_recs)
c
c     Local variables.
c
      integer*4 rec,ifg_pos

      real*4 divisor,dither,sort(492)

      fil_convert = %loc(fil_normal)
c
c     Divide the good interferograms by the gain and the number of sweeps.
c
      do rec = 1,num_recs
         if (con_recs(rec).con_check .eq. 0) then
            divisor = con_recs(rec).gain*sweeps

            do ifg_pos = 1,512
               aifgs(ifg_pos,rec) = sci_recs(rec).ifg_data.ifg(ifg_pos)/divisor
            end do
c
c     Calculate the dither by finding the median value at each point in the
c     group, and subtract it.  Omit the first twenty points of the
c     interferogram to minimize possible digital transient effects.
c
            do ifg_pos = 1,492
               sort(ifg_pos) = aifgs(ifg_pos+20,rec)
            end do
            call svrgn(492,sort,sort)
            dither = (sort(246) + sort(247))/2.0
            do ifg_pos = 1,512
               aifgs(ifg_pos,rec) = aifgs(ifg_pos,rec) - dither
            end do
         end if
      end do

      return
      end
