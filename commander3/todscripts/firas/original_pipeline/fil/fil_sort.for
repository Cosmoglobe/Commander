      integer*4 function fil_sort(channel,neighbors,ifgs_max,group_max,cth,
     .                            num_recs,num_cgr_recs,sci_recs,num_good,
     .                            num_cgr_good,con_recs)
c-------------------------------------------------------------------------------
c
c     Purpose: Sort the neighbor interferograms and add the appropriate number
c              of the best ones to the coadd group for template formation.
c
c     Author: S. Brodd, HSTX, 4/95
c
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            neighbors          i*4  Use neighbors in template formation.
c            ifgs_max           i*4  Value of neighbors qualifier.
c            group_max          i*4  Number of interferograms, including
c                                    neighbors, to use in template formation.
c            cth                rec  Consistency check parameters.
c            num_recs           i*4  Number of input science and engineering
c                                    records, including neighbors.
c            num_cgr_recs       i*4  Number of input science and engineering
c                                    records, excluding neighbors.
c            sci_recs           rec(num_recs)  Science records.
c            num_good           i*4  Current number of good interferograms
c                                    including neighbors.
c            num_cgr_good       i*4  Current number of good interferograms
c                                    excluding neighbors.
c            con_recs           rec(num_recs)  Consistency check records.
c
c     Output: num_good           i*4  Current number of good interferograms
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
c     Input parameters.
c
      integer*4 channel,neighbors,ifgs_max,group_max
      integer*4 num_recs,num_cgr_recs,num_cgr_good

      dictionary 'fex_cth'
      dictionary 'fdq_sdf'

      record /fex_cth/ cth
      record /fdq_sdf/ sci_recs(num_recs)
c
c     Input/output parameters.
c
      integer*4 num_good

      dictionary 'fil_scc'

      record /fil_scc/ con_recs(num_recs)
c
c     Local variables.
c
      integer*4 rec,rec_del
      real*4 ave,val,val_most

      fil_sort = %loc(fil_normal)
c
c     If the number of good interferograms excluding neighbors is greater than 
c     or equal to group_max or is greater than ifgs_max from the command line, 
c     then use no neighbors by marking them all as having bad data quality.
c
      if ((num_cgr_good .ge. group_max) .or. ((neighbors .eq. fac_present) 
     .    .and. (num_cgr_good .gt. ifgs_max))) then
         num_good = num_cgr_good
         do rec = num_cgr_recs+1,num_recs
            if (con_recs(rec).con_check .eq. 0) then
               con_recs(rec).con_check = ibset(con_recs(rec).con_check,0)
            end if
         end do
c
c     Otherwise, if there are more good interferograms including neighbors than
c     group_max from the command line, some of the neighbors must be
c     eliminated.
c
      else if (num_good .gt. group_max) then
c
c     Determine the average galactic latitude of the group.
c
         ave = 0.0
         do rec = 1,num_recs   
            if (con_recs(rec).con_check .eq. 0) then
               ave = ave + sci_recs(rec).attitude.galactic_latitude*fac_att_conv
            end if
         end do
         ave = ave / num_good
c
c     If the average galactic latitude is less than the cutoff value from the
c     cth file, then use galactic latitude as the criterion for deciding among
c     the neighbors.  Otherwise, use the glitch rate.
c
         if (abs(ave) .lt. cth.gal_lat_cutoff(channel)) then
c
c     Until the number of good interferograms including neighbors is less than
c     or equal to group_max, find the neighbor interferogram with the greatest
c     absolute value of the difference between that interferograms galactic 
c     latitude and the average and mark it as bad.
c
            do while (num_good .gt. group_max)
               val_most = 0.0
               do rec = num_cgr_recs+1,num_recs   
                  if (con_recs(rec).con_check .eq. 0) then
                     val = abs(ave - sci_recs(rec).attitude.galactic_latitude * 
     .                     fac_att_conv)
                     if (val .gt. val_most) then
                        rec_del = rec 
                        val_most = val
                     end if
                  end if
               end do
               num_good = num_good - 1
               con_recs(rec_del).con_check=ibset(con_recs(rec_del).con_check,0)
            end do
         else
c
c     Until the number of good interferograms including neighbors is less than
c     or equal to group_max, find the neighbor interferogram with the greatest
c     glitch rate and mark it as bad.
c
            do while (num_good .gt. group_max)
               val_most = 0.0
               do rec = num_cgr_recs+1,num_recs   
                  if (con_recs(rec).con_check .eq. 0) then
                     val = con_recs(rec).glitch_rate
                     if (val .gt. val_most) then
                        rec_del = rec 
                        val_most = val
                     end if
                  end if
               end do
               num_good = num_good - 1
               con_recs(rec_del).con_check=ibset(con_recs(rec_del).con_check,0)
            end do
         end if
      end if

      return
      end
