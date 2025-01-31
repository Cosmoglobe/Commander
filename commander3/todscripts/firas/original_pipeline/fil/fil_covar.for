      integer*4 function fil_covar(temp_scan_mode,num_vecs,adj_num_vecs,glrt,
     .                             dis_vecs,sci_mode,adds_per_group,coa_vec,
     .                             bins,cov_recs)
c-------------------------------------------------------------------------------
c
c     Purpose: Calculate covariance matrix and associated statistical
c              quantities.
c
c     Author: S. Brodd, HSTX, 4/95
c
c     Input: temp_scan_mode     i*4  Converts scan mode 4 to 3, 5 to 4, and
c                                    6 to 5 for covariance matrix calculation.
c            num_vecs           i*4  Number of dispersion vectors.
c            adj_num_vecs       r*4  Glitch rate weighted number of dispersion 
c                                    vectors.
c            glrt               r*4(num_vecs)  Glitch rates of interferograms.
c            dis_vecs           c*16(369,num_vecs)  Interferogram dispersion
c                                                   vectors.
c            sci_mode           i*4  Value of science mode, 0-4.
c            adds_per_group     i*4  Value of adds per group, 1-12.
c            coa_vec            c*16(369)  Coadd dispersion vector.
c            bins               i*4(num_vecs)  Statistical bins.
c            cov_recs           rec(1285)  Covariance matrix records.
c
c     Output: cov_recs           rec(1285)  Covariance matrix records.
c
c     Modifications: S. Brodd, HSTX, 9/95, SPR 12256. Add statistical data 
c                    associated with glitch rate.
c    
c-------------------------------------------------------------------------------

      implicit none
c
c     Return statuses.
c
      external fil_normal
c
c     Input parameters.
c
      integer*4 temp_scan_mode,num_vecs
      real*4 adj_num_vecs,glrt(num_vecs)
      complex*16 dis_vecs(369,num_vecs),coa_vec(369)
      integer*4 sci_mode,adds_per_group
      integer*4 bins(num_vecs)
c
c     Input/output parameters.
c
      dictionary 'fil_cov'
      record /fil_cov/ cov_recs(1285)
c
c     Local variables.
c
      integer*4 offset,rec,loc,vec,vec_pos,row,col
      real*8 quot,glrt_sum,norm,col_vec(369),val

      fil_covar = %loc(fil_normal)
c
c     Calculate the offset to the correct group of 257 covariance matrix
c     records associated with the current scan mode; SS=1-257, SF=258-514,
c     LF=515-771, FS=772-1028, FL=1029-1285.
c
      offset = (temp_scan_mode - 1)*257
c
c     Put the science mode and adds per group values in each covariance matrix
c     record.
c
      do rec = 1,257
         loc = rec + offset
         if (cov_recs(loc).ident.sci_mode .eq. 0) then
            cov_recs(loc).ident.sci_mode = sci_mode
         else if (cov_recs(loc).ident.sci_mode .ne. sci_mode) then
            cov_recs(loc).ident.sci_mode = -1
         end if

         if (cov_recs(loc).ident.adds_per_group .eq. 0) then
            cov_recs(loc).ident.adds_per_group = adds_per_group
         else if (cov_recs(loc).ident.adds_per_group .ne. 
     .            adds_per_group) then
            cov_recs(loc).ident.adds_per_group = -1
         end if
      end do
c
c     Increment the number of interferograms and coadds contributing to the
c     covariance matrix in the first record.
c
      rec = offset + 1
      cov_recs(rec).num_ifgs = cov_recs(rec).num_ifgs + num_vecs
      cov_recs(rec).adj_num_ifgs = cov_recs(rec).adj_num_ifgs + adj_num_vecs
      cov_recs(rec).num_coadds = cov_recs(rec).num_coadds + 1
c
c     Accumulate the statistical data from the coadded voltage spectrum in
c     the first covariance matrix record.
c
      do vec_pos = 1,361
         cov_recs(rec).avg_rvoltspec(vec_pos) = 
     .                 cov_recs(rec).avg_rvoltspec(vec_pos) +
     .                 num_vecs * dreal(coa_vec(vec_pos))
         cov_recs(rec).avg_ivoltspec(vec_pos) = 
     .                 cov_recs(rec).avg_ivoltspec(vec_pos) +
     .                 num_vecs * dimag(coa_vec(vec_pos))
         cov_recs(rec).volt_rsqsum(vec_pos) = 
     .                 cov_recs(rec).volt_rsqsum(vec_pos) +
     .                 num_vecs * dreal(coa_vec(vec_pos))**2
         cov_recs(rec).volt_isqsum(vec_pos) = 
     .                 cov_recs(rec).volt_isqsum(vec_pos) +
     .                 num_vecs * dimag(coa_vec(vec_pos))**2
      end do
c
c     Accumulate the statistical data associated with glitch rates in the
c     first covariance matrix record.
c
      quot = (num_vecs - 1.0) / num_vecs
      glrt_sum = 0.0
      do vec = 1,num_vecs
         glrt_sum = glrt_sum + glrt(vec)
         do vec_pos = 1,334
            cov_recs(rec).wtd_disp_tot(vec_pos) = 
     .                    cov_recs(rec).wtd_disp_tot(vec_pos) + 
     .                    (abs(dis_vecs(vec_pos,vec)))**2 * glrt(vec)
         end do
      end do
      cov_recs(rec).sum_quots = cov_recs(rec).sum_quots + quot
      cov_recs(rec).sum_wts = cov_recs(rec).sum_wts + (quot * glrt_sum)
      cov_recs(rec).sum_wts_sq = cov_recs(rec).sum_wts_sq + 
     .                           (quot * glrt_sum**2)
c
c     Loop over the dispersion vectors formed from good interferograms.
c
      do vec = 1,num_vecs
         norm = 0.0
c
c     Locate the appropriate record based on the statistical bin of the
c     dispersion vector.
c
         rec = bins(vec) + offset + 1
c
c     Increment the bin total and the weighted bin total.
c
         cov_recs(rec).bin_total = cov_recs(rec).bin_total + 1
         cov_recs(rec).wtd_bin_total = cov_recs(rec).wtd_bin_total +
     .                                 1.0 - 1.0/num_vecs
c
c     Accumulate the dispersion vectors.
c
         do vec_pos = 1,361
            cov_recs(rec).disp(vec_pos) = cov_recs(rec).disp(vec_pos) +
     .                                    dis_vecs(vec_pos,vec)
            norm = norm + dis_vecs(vec_pos,vec) * conjg(dis_vecs(vec_pos,vec))
         end do
c
c     Calculate the squares and fourth powers of the norms of the dispersion
c     vector totals.
c
         cov_recs(rec).norms(1) = cov_recs(rec).norms(1) + norm
         cov_recs(rec).norms(2) = cov_recs(rec).norms(2) + norm**2
c
c     Accumulate the temperature and glitch rate information in the 
c     dispersion vectors.
c
         do vec_pos = 1,8
            cov_recs(rec).temp(vec_pos) = cov_recs(rec).temp(vec_pos) +
     .                                    dreal(dis_vecs(vec_pos+361,vec))
         end do
c
c     Extract the real parts of the dispersion vector for subsequent use.
c
         do col = 1,369
            col_vec(col) = dreal(dis_vecs(col,vec))
         end do
c
c     Calculate the real-real portion of the covariance matrix and store it in 
c     records 2-67 plus the scan mode offset with 1050 values in each record;
c     there will be 15 values in the last record for a total of 68265 values.
c
         loc = 0
         rec = offset + 2
         do row = 1,369
            val = dreal(dis_vecs(row,vec))
            do col = row,369
               loc = loc + 1
               if (loc .gt. 1050) then
                  rec = rec + 1
                  loc = 1
               end if
               cov_recs(rec).covar(loc) = cov_recs(rec).covar(loc) + 
     .                                    val * col_vec(col)
            end do
         end do
c
c     Extract the imaginary parts of the dispersion vector for subsequent use.
c
         do col = 1,361
            col_vec(col) = dimag(dis_vecs(col,vec))
         end do
c
c     Calculate the real-imaginary portion of the covariance matrix and store 
c     it in records 68-194 plus the scan mode offset with 1050 values in each 
c     record; there will be 909 values in the last record for a total of
c     133209 values.
c
         loc = 0
         rec = offset + 68
         do row = 1,369
            val = dreal(dis_vecs(row,vec))
            do col = 1,361
               loc = loc + 1
               if (loc .gt. 1050) then
                  rec = rec + 1
                  loc = 1
               end if
               cov_recs(rec).covar(loc) = cov_recs(rec).covar(loc) +
     .                                    val * col_vec(col)
            end do
         end do
c
c     Calculate the imaginary-imaginary portion of the covariance matrix and 
c     store it in records 195-257 plus the scan mode offset with 1050 values in 
c     each record; there will be 241 values in the last record for a total of
c     65341 values.
c
         loc = 0
         rec = offset + 195
         do row = 1,361
            val = dimag(dis_vecs(row,vec))
            do col = row,361
               loc = loc + 1
               if (loc .gt. 1050) then
                  rec = rec + 1
                  loc = 1
               end if
               cov_recs(rec).covar(loc) = cov_recs(rec).covar(loc) +
     .                                    val * col_vec(col)
            end do
         end do
      end do

      return
      end
