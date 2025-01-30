      integer*4 function fil_template(num_recs,num_good,num_cgr_good,con_recs,
     .                                coa_rec,aifgs,template)
c-----------------------------------------------------------------------------
c
c     Purpose: Calculate and subtract the primary template.
c
c     Authors: R. Isaacman, ARC, 1/86
c              W. Young, SASC, 2/86
c              S. Brodd, HSTX, 4/95
c
c     Input: num_recs           i*4  Number of input science and engineering
c                                    records, including neighbors.
c            num_good           i*4  Current number of good interferograms,
c                                    including neighbors.
c            num_cgr_good       i*4  Current number of good interferograms,
c                                    excluding neighbors.
c            con_recs           rec(num_recs)  Consistency check records.
c            coa_rec            rec  Coadd record.
c            aifgs              r*4(512,num_recs)  Real-valued interferograms.
c
c     Output: coa_rec            rec  Coadd record.
c             aifgs              r*4(512,num_recs)  Real-valued interferograms.
c             template           r*4(512)  Primary template.
c
c     Modifications: S. Brodd, HSTX, 9/95, SPR 12256. Change algorithm for 
c                    cases of num_good = 6, 10, 14, ... to prevent bias.
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
      external fil_normal
      external fil_tempzero
c
c     Input parameters.
c
      integer*4 num_recs,num_good,num_cgr_good

      dictionary 'fil_scc'
      record /fil_scc/ con_recs(num_recs)
c
c     Input/output parameters.
c
      dictionary 'fil_sky'
      record /fil_sky/ coa_rec

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

      fil_template = %loc(fil_normal)
c
c     Form the primary template by sorting the coadd group and finding the
c     midaverage at each of the 512 interferogram points.
c
      first = nint(num_good/4.0)
      last = nint((3.0*num_good)/4.0)
      num_midav_pts = last - first
c    
c     For the case where num_good is 6, 10, 14, ..., decrement the last pointer
c     so that the central 2, 4, 6, ..., points are initially included in the
c     midaverage.  The num_midav_pts parameter remains set at 3, 5, 7, ...,
c     because half of each of points (2,5), (3,8), (4,11) ..., will be added.
c
      if (mod(num_good,4) .eq. 2) then
         last = last - 1
      end if

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
c
c     Add in half of each of points (2,5), (3,8), (4,11) ... to avoid biasing
c     the midaverage.
c
         if (mod(num_good,4) .eq. 2) then
            template(ifg_pos) = template(ifg_pos) + 0.5 * sort(first) + 
     .                                              0.5 * sort(last+1)
         end if
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
c     Put information about interferograms used in template formation including
c     neighbors, if any, in the coadd record.
c
      coa_rec.coad_spec_data.template.num_ifgs = num_good
      if (num_good .gt. num_cgr_good) then
         coa_rec.coad_spec_data.template.neighbors = fac_present
         coa_rec.coad_spec_data.template.neighbor_num_ifgs = 
     .                                   num_good - num_cgr_good
      else
         coa_rec.coad_spec_data.template.neighbors = fac_not_present
         coa_rec.coad_spec_data.template.neighbor_num_ifgs = 0
      end if

      rec_good = 0
      do rec=1,num_recs
         if (con_recs(rec).con_check .eq. 0) then
            rec_good = rec_good + 1
            coa_rec.coad_spec_data.template.ifgs(rec_good).time(1) =
     .                                      con_recs(rec).time(1)
            coa_rec.coad_spec_data.template.ifgs(rec_good).time(2) =
     .                                      con_recs(rec).time(2)
            coa_rec.coad_spec_data.template.ifgs(rec_good).pixel_no =
     .                                      con_recs(rec).pixel_no
         end if
      end do
c
c     If the primary template sum is zero, an error has occurred.
c
      if (sum .eq. 0.0) then
         fil_template = %loc(fil_tempzero)
         call lib$signal(fil_tempzero,%val(2),
     .                   fac_channel_ids(con_recs(1).chan_id),
     .                   con_recs(1).coadd_gmt)
      end if

      return
      end
