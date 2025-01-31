      integer*4 function fic_covar(temp_scan_mode,num_vecs,dis_vecs,sci_mode,
     .                             adds_per_group,coa_vec,bins,cov_recs)
c-------------------------------------------------------------------------------
c
c     Purpose: Calculate covariance matrix and associated statistical
c              quantities.
c
c     Author: S. Alexander, STX, 9/91, SER 7985
c
c     Input: temp_scan_mode     i*4  Converts scan mode 4 to 3 for covariance
c                                    matrix calculation.
c            num_vecs           i*4  Number of dispersion vectors.
c            dis_vecs           c*16(265,num_recs)  Interferogram dispersion
c                                                   vectors.
c            sci_mode           i*4  Value of science mode, 0-4.
c            adds_per_group     i*4  Value of adds per group, 1-12.
c            coa_vec            c*16(265)  Coadd dispersion vector.
c            bins               i*4(100)  Statistical bins for interferograms.
c            cov_recs           rec(771)  Covariance matrix records.
c
c     Output: cov_recs           rec(771)  Covariance matrix records.
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
      integer*4 temp_scan_mode,num_vecs
      complex*16 dis_vecs(265,num_vecs),coa_vec(265)
      integer*4 sci_mode,adds_per_group
      integer*4 bins(num_vecs)
c
c     Input/output parameters.
c
      dictionary 'fic_cov'
      record /fic_cov/ cov_recs(771)
c
c     Local variables.
c
      integer*4 offset,rec,loc,vec,vec_pos,row,col
      real*8 norm,col_vec(265),val

      fic_covar = %loc(fic_normal)
c
c     Calculate the offset to the correct group of 257 covariance matrix
c     records associated with the current scan mode; SS=1-257, SF=258-514,
c     LF=515-771.
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
      cov_recs(rec).num_coadds = cov_recs(rec).num_coadds + 1
c
c     Accumulate the statistical data from the coadded voltage spectrum in
c     the first covariance matrix record.
c
      do vec_pos = 1,257
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
     .                      cov_recs(rec).volt_isqsum(vec_pos) +
     .                      num_vecs * dimag(coa_vec(vec_pos))**2
      end do
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
         do vec_pos = 1,257
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
     .                                    dreal(dis_vecs(vec_pos+257,vec))
         end do
c
c     Extract the real parts of the dispersion vector for subsequent use.
c
         do col = 1,265
            col_vec(col) = dreal(dis_vecs(col,vec))
         end do
c
c     Calculate the real-real portion of the covariance matrix and store it in 
c     records 2-66 plus the scan mode offset with 550 values in each record;
c     there will be 45 values in the last record for a total of 35245 values.
c
         loc = 0
         rec = offset + 2
         do row = 1,265
            val = dreal(dis_vecs(row,vec))
            do col = row,265
               loc = loc + 1
               if (loc .gt. 550) then
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
         do col = 1,257
            col_vec(col) = dimag(dis_vecs(col,vec))
         end do
c
c     Calculate the real-imaginary portion of the covariance matrix and store 
c     it in records 67-190 plus the scan mode offset with 550 values in each 
c     record; there will be 455 values in the last record for a total of
c     68105 values.
c
         loc = 0
         rec = offset + 67
         do row = 1,265
            val = dreal(dis_vecs(row,vec))
            do col = 1,257
               loc = loc + 1
               if (loc .gt. 550) then
                  rec = rec + 1
                  loc = 1
               end if
               cov_recs(rec).covar(loc) = cov_recs(rec).covar(loc) +
     .                                    val * col_vec(col)
            end do
         end do
c
c     Calculate the imaginary-imaginary portion of the covariance matrix and 
c     store it in records 191-251 plus the scan mode offset with 550 values in 
c     each record; there will be 153 values in the last record for a total of
c     33153 values.
c
         loc = 0
         rec = offset + 191
         do row = 1,257
            val = dimag(dis_vecs(row,vec))
            do col = row,257
               loc = loc + 1
               if (loc .gt. 550) then
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
