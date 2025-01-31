      integer*4 function fic_template(num_recs,num_good,con_recs,coa_rec,aifgs,
     .                                template)
c-----------------------------------------------------------------------------
c
c     Purpose: Calculate and subtract the primary template.
c
c     Authors: R. Isaacman, ARC, 1/29/86
c              W. Young, SASC, 2/86
c              S. Alexander, STX, 9/91, SER 7985
c
c     Input: num_recs           i*4  Number of input science and engineering
c                                    records.
c            num_good           i*4  Current number of good interferograms.
c            con_recs           rec(num_recs)  Consistency check records.
c            coa_rec            rec  Coadd record.
c            aifgs              r*4(512,num_recs)  Real-valued interferograms.
c
c     Output: coa_rec            rec  Coadd record.
c             aifgs              r*4(512,num_recs)  Real-valued interferograms.
c             template           r*4(512)  Primary template.
c
c     Modifications:
c
c-----------------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_params)'
c
c     Return statuses.
c
      external fic_normal
      external fic_tempzero
c
c     Input parameters.
c
      integer*4 num_recs,num_good

      dictionary 'fic_scc'
      record /fic_scc/ con_recs(num_recs)
c
c     Input/output parameters.
c
      dictionary 'fic_sky'
      record /fic_sky/ coa_rec

      real*4 aifgs(512,num_recs)
c
c     Output parameters.
c
      real*4 template(512)
c
c     Local variables.
c
      integer*4 first,last,num_midav_pts
      real*4 sum,sort(512)
      integer*4 rec_good,ifg_pos,rec,mid

      fic_template = %loc(fic_normal)
c
c     Form the primary template by sorting the coadd group and finding the
c     midaverage at each of the 512 interferogram points.
c
      first = nint(num_good/4.0)
      last = nint((3.0*num_good)/4.0)
      num_midav_pts = last - first

      sum = 0.0

      do ifg_pos=1,512
         rec_good = 0
         do rec=1,num_recs
            if (con_recs(rec).con_check .eq. 0) then
               rec_good = rec_good + 1
               sort(rec_good) = aifgs(ifg_pos,rec)
            end if
         end do
         call svrgn(num_good,sort,sort)
         template(ifg_pos) = 0.0
         do mid = first+1,last
            template(ifg_pos) = template(ifg_pos) + sort(mid)
         end do
         template(ifg_pos) = template(ifg_pos)/num_midav_pts
         sum = sum + abs(template(ifg_pos))          
c
c      Subtract the primary template.
c
         do rec = 1,num_recs
            if (con_recs(rec).con_check .eq. 0) then
               aifgs(ifg_pos,rec) = aifgs(ifg_pos,rec) - template(ifg_pos)
            end if
         end do
      end do

      coa_rec.coad_spec_data.prim_template.subtracted = fac_present
c
c     If the primary template sum is zero, an error has occurred.
c
      if (sum .eq. 0.0) then
         fic_template = %loc(fic_tempzero)
         call lib$signal(fic_tempzero,%val(2),
     .                   fac_channel_ids(con_recs(1).chan_id),
     .                   con_recs(1).coadd_gmt)
      end if

      return
      end
