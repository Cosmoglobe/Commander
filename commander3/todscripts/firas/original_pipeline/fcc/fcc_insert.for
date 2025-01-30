      integer*4 function fcc_insert(c,n,n2,eq,da,out_cov)
c-------------------------------------------------------------------------------
c
c     Purpose: Calculate and store the output covariance matrix from the four 
c              calibrated intermediate matrices; store the calibrated coadd 
c              and binned dispersion vectors.
c
c     Author: S. Alexander, HSTX, 7/93, SER 11189
c
c     Input: c            r*8(522,522)   Output covariance matrix
c            n            c*16(512,512)  Intermediate matrix N (calibrated)
c            n2           c*16(512,512)  Intermediate matrix N2 (calibrated)
c            eq           c*16(8,512)    Intermediate matrix EQ (calibrated)
c            da           c*16(257,512)  Intermediate matrix DA (calibrated)
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
      external fcc_normal
c
c     Input parameters.
c
      real*8 c(522,522)
      complex*16 n(512,512),n2(512,512),eq(8,512),da(257,512)
c
c     Input/output parameters.
c
      dictionary 'fcc_cov'
      record /fcc_cov/ out_cov(257)
c
c     Local variables.
c
      integer*4 deg_freedom,num_ifgs,rec,pos,row,col

      fcc_insert = %loc(fcc_normal)
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
         do pos = 1,257
            out_cov(rec).rdisp(pos) = dreal(da(rec-1,pos))
            out_cov(rec).idisp(pos) = dimag(da(rec-1,pos))
         end do
      end do
c
c     Put the real and imaginary parts of the 257th row of DA into the first 
c     output record as the calibrated coadd voltage spectra.
c
      do pos = 1,257
         out_cov(1).avg_rcalspec(pos) = dreal(da(257,pos)) / num_ifgs
         out_cov(1).avg_icalspec(pos) = dimag(da(257,pos)) / num_ifgs
      end do
c
c     Construct the upper right half of the 522 by 522 FCC symmetric calibrated 
c     output covariance matrix C from the intermediate matrices N, N2, and EQ.
c
      do row = 1,257
         do col = row,257
            c(row,col) = (dreal(n(row,col)) + dreal(n2(row,col))) / 
     .                   (2*deg_freedom)
         end do
         do col = 258,265
            c(row,col) = dreal(eq(col-257,row)) / deg_freedom
         end do
         do col = 266,522
            c(row,col) = (dimag(n(col-265,row)) + dimag(n2(row,col-265))) / 
     .                   (2*deg_freedom)
         end do
      end do
      do row = 258,265
         do col = row,265	
            c(row,col) = c(row,col) / deg_freedom
         end do
         do col = 266,522
            c(row,col) = dimag(eq(row-257,col-265)) / deg_freedom
         end do
      end do
      do row = 266,522
         do col = row,522
            c(row,col) = (dreal(n(row-265,col-265)) - 
     .                    dreal(n2(row-265,col-265))) / (2*deg_freedom)
         end do
      end do
c
c     Store the calibrated output covariance matrix in the output records.
c     The real-real portion is stored in records 2-66 with 550 values in
c     records 2-65 and 45 values in record 66 for a total of 35245 values.
c
      pos = 0
      rec = 2
      do row = 1,265
         do col = row,265
            pos = pos + 1
            if (pos .gt. 550) then
               rec = rec + 1
               pos = 1
            end if
            out_cov(rec).covar(pos) = c(row,col)
         end do
      end do
c
c     The real-imaginary portion is stored in records 67-190 with 550 values in
c     records 67-189 and 455 values in record 190 for a total of 68105 values.
c
      pos = 0
      rec = 67
      do row = 1,265
         do col = 266,522
            pos = pos + 1
            if (pos .gt. 550) then
               rec = rec + 1
               pos = 1
            end if
            out_cov(rec).covar(pos) = c(row,col)
         end do
      end do
c
c     The imaginary-imaginary portion is stored in records 191-251 with 550 
c     values in records 191-250 and 153 values in record 251 for a total of 
c     33153 values.
c
      pos = 0
      rec = 191
      do row = 266,522
         do col = row,522
            pos = pos + 1
            if (pos .gt. 550) then
               rec = rec + 1
               pos = 1
            end if
            out_cov(rec).covar(pos) = c(row,col)
         end do
      end do

      return
      end
