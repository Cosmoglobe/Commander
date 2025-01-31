      integer*4 function fcl_insert(c,n,n2,eq,da,out_cov)
c-------------------------------------------------------------------------------
c
c     Purpose: Calculate and store the output covariance matrix from the four 
c              calibrated intermediate matrices; store the calibrated coadd 
c              and binned dispersion vectors.
c
c     Author: S. Brodd, HSTX, 12/95, SPR 12291
c
c     Input: c            r*8(730,730)   Output covariance matrix
c            n            c*16(361,361)  Intermediate matrix N (calibrated)
c            n2           c*16(361,361)  Intermediate matrix N2 (calibrated)
c            eq           c*16(8,361)    Intermediate matrix EQ (calibrated)
c            da           c*16(257,361)  Intermediate matrix DA (calibrated)
c            out_cov      rec(257)       Output covariance matrix records
c
c     Output: out_cov      rec(257)       Output covariance matrix records
c
c     Modifications:
c
c-------------------------------------------------------------------------------
      implicit none
c
c     Return status.
c
      external fcl_normal
c
c     Input parameters.
c
      real*8 c(730,730)
      complex*16 n(361,361),n2(361,361),eq(8,361),da(257,361)
c
c     Input/output parameters.
c
      dictionary 'fcl_cov'
      record /fcl_cov/ out_cov(257)
c
c     Local variables.
c
      integer*4 deg_freedom,num_ifgs,rec,pos,row,col

      fcl_insert = %loc(fcl_normal)
c
c     Retrieve the degrees of freedom and the number of interferograms from
c     the first output record.
c
      deg_freedom = out_cov(1).deg_freedom
      num_ifgs = out_cov(1).num_ifgs
c
c     Put the real and imaginary parts of the first 256 rows of the
c     intermediate matrix DA into the output records as the calibrated binned 
c     dispersion totals.
c
      do rec = 2,257
         do pos = 1,361
            out_cov(rec).rdisp(pos) = dreal(da(rec-1,pos))
            out_cov(rec).idisp(pos) = dimag(da(rec-1,pos))
         end do
      end do
c
c     Put the real and imaginary parts of the 257th row of DA into the first 
c     output record as the calibrated coadd voltage spectra.
c
      do pos = 1,361
         out_cov(1).avg_rcalspec(pos) = dreal(da(257,pos)) / num_ifgs
         out_cov(1).avg_icalspec(pos) = dimag(da(257,pos)) / num_ifgs
      end do
c
c     Construct the upper right half of the 730 by 730 FCL symmetric calibrated 
c     output covariance matrix C from the intermediate matrices N, N2, and EQ.
c
      do row = 1,361
         do col = row,361
            c(row,col) = (dreal(n(row,col)) + dreal(n2(row,col))) / 
     .                   (2.0 * deg_freedom)
         end do
         do col = 362,369
            c(row,col) = dreal(eq(col-361,row)) / deg_freedom
         end do
         do col = 370,730
            c(row,col) = (dimag(n(col-369,row)) + dimag(n2(row,col-369))) / 
     .                   (2.0 * deg_freedom)
         end do
      end do
      do row = 362,369
         do col = row,369	
            c(row,col) = c(row,col) / deg_freedom
         end do
         do col = 370,730
            c(row,col) = dimag(eq(row-361,col-369)) / deg_freedom
         end do
      end do
      do row = 370,730
         do col = row,730
            c(row,col) = (dreal(n(row-369,col-369)) - 
     .                    dreal(n2(row-369,col-369))) / (2.0 * deg_freedom)
         end do
      end do
c
c     Store the calibrated output covariance matrix in the output records.
c     The real-real portion is stored in records 2-67 with 1050 values in
c     records 2-66 and 15 values in record 67 for a total of 68265 values.
c
      pos = 0
      rec = 2
      do row = 1,369
         do col = row,369
            pos = pos + 1
            if (pos .gt. 1050) then
               rec = rec + 1
               pos = 1
            end if
            out_cov(rec).covar(pos) = c(row,col)
         end do
      end do
c
c     The real-imaginary portion is stored in records 68-194 with 1050 values in
c     records 68-193 and 909 values in record 190 for a total of 133209 values.
c
      pos = 0
      rec = 68
      do row = 1,369
         do col = 370,730
            pos = pos + 1
            if (pos .gt. 1050) then
               rec = rec + 1
               pos = 1
            end if
            out_cov(rec).covar(pos) = c(row,col)
         end do
      end do
c
c     The imaginary-imaginary portion is stored in records 195-257 with 1050 
c     values in records 195-256 and 241 values in record 257 for a total of 
c     65341 values.
c
      pos = 0
      rec = 195
      do row = 370,730
         do col = row,730
            pos = pos + 1
            if (pos .gt. 1050) then
               rec = rec + 1
               pos = 1
            end if
            out_cov(rec).covar(pos) = c(row,col)
         end do
      end do

      return
      end
