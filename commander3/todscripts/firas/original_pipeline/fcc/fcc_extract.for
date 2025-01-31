      integer*4 function fcc_extract(label,label_val,bol_avg,cbias_avg,
     .                               volt_avg,model,s0,tau,tbol,in_cov,out_cov,
     .                               c,n,n2,eq,da)
c-------------------------------------------------------------------------------
c
c     Purpose: Copy unchanged quantities from the input covariance matrix to
c              the output, calculate the degrees of freedom and effective
c              weights and create the four intermediate matrices N, N2, EQ,
c              and DA.
c
c     Author: S. Alexander, HSTX, 7/93, SER 11189
c
c     Input: label      l*1  Label present or not.
c            label_val  ch*40  Optional label for output covariance matrix.
c            bol_avg    r*4  Average bolometer temperature in degrees K.
c            cbias_avg  r*4  Average commanded bias in counts.
c            volt_avg   r*4  Average readout voltage in volts.
c            model      rec  Calibration model.
c            s0         r*8  Detector responsivity.
c            tau        r*8  Detector time constant.
c            tbol       r*8  Derived bolometer temperature.
c            in_cov     rec(257)  Input covariance matrix records.
c
c     Output: out_cov  rec(257)  Output covariance matrix records.
c             c        r*8(522,522)   Reconstructed covariance matrix.
c             n        c*16(512,512)  Intermediate matrix N.
c             n2       c*16(512,512)  Intermediate matrix N2.
c             eq       c*16(8,512)    Intermediate matrix EQ.
c             da       c*16(257,512)  Intermediate matrix DA.
c
c     Modifications:
c
c-------------------------------------------------------------------------------
      implicit none
c
c     Include files.
c
      include '(fut_error)'
c
c     Return statuses.
c
      external fcc_normal
      external fcc_repwrite
c
c     External references.
c
      external fut_error
c
c     Input parameters.
c
      logical*1 label
      character*40 label_val

      real*4 bol_avg,cbias_avg,volt_avg

      dictionary 'fex_mod'
      record /fex_mod/ model

      real*8 s0,tau,tbol

      dictionary 'fic_cov'
      record /fic_cov/ in_cov(257)
c
c     Output parameters.
c
      dictionary 'fcc_cov'
      record /fcc_cov/ out_cov(257)

      real*8 c(522,522)
      complex*16 n(512,512),n2(512,512),eq(8,512),da(257,512)
c
c     Local variables.
c
      integer*4 rec,num_ifgs,num_coadds,pos,sum_st,deg_freedom
      real*4 prop_st
      integer*4 status,len,row,col
      character*72 filename
      real*8 sum_norm,avg_norm

      fcc_extract = %loc(fcc_normal)
c
c     Put the model label and optional input label into the output records;
c     copy the instrument parameters from input to output.
c
      do rec = 1,257
         out_cov(rec).model_label = model.mod_head.label
         if (label) then
            out_cov(rec).cov_label = label_val
         end if
         out_cov(rec).ident.chan_id = in_cov(rec).ident.chan_id
         out_cov(rec).ident.mtm_speed = in_cov(rec).ident.mtm_speed
         out_cov(rec).ident.mtm_length = in_cov(rec).ident.mtm_length
         out_cov(rec).ident.sci_mode = in_cov(rec).ident.sci_mode
         out_cov(rec).ident.adds_per_group = in_cov(rec).ident.adds_per_group
      end do
c
c     Copy the command line instrument values and the calculated detector
c     values into the first output record.
c
      out_cov(1).bol_avg = bol_avg
      out_cov(1).cbias_avg = cbias_avg
      out_cov(1).volt_avg = volt_avg
      out_cov(1).s0 = s0
      out_cov(1).tau = tau
      out_cov(1).tbol = tbol
c
c     Copy the number of interferograms and number of coadds from input to
c     output.
c
      num_ifgs = in_cov(1).num_ifgs
      out_cov(1).num_ifgs = num_ifgs
      num_coadds = in_cov(1).num_coadds
      out_cov(1).num_coadds = num_coadds
c
c     Copy the bin totals, weighted bin totals, and temperature and glitch 
c     rate dispersion totals from input to output.
c
      do rec = 2,257
         out_cov(rec).bin_total = in_cov(rec).bin_total
         out_cov(rec).wtd_bin_total = in_cov(rec).wtd_bin_total
         do pos = 1,8
            out_cov(rec).temp(pos) = in_cov(rec).temp(pos)
         end do
      end do
c
c     Sum the bin totals over bins 129 through 256; this is the number of
c     interferograms from which the secondary template was subtracted.
c
      sum_st = 0
      do rec = 130,257
         sum_st = sum_st + in_cov(rec).bin_total
      end do
c
c     Divide this sum by the number of interferograms; this is the proportion
c     of interferograms from which the secondary template was subtracted.
c      
      prop_st = float(sum_st)/num_ifgs
c
c     Determine the number of degrees of freedom and write it to the output
c     covariance matrix.
c
      deg_freedom = nint(num_ifgs - (1.0 + prop_st) * num_coadds)
      out_cov(1).deg_freedom = deg_freedom
c
c     Write the three calculated quantities to the report file.
c
      write(fut_report_lun,10,iostat=status) 
     .      'Number of Interferograms Contributing to Matrix:',num_ifgs,
     .      'Number of Coadds Contributing to Matrix:',num_coadds,
     .      'Number of Interferograms with Secondary Template Subtracted:',
     .       sum_st,
     .      'Proportion of Interferograms with Secondary Template Subtracted:',
     .       prop_st,
     .      'Degrees of Freedom:',deg_freedom
10    format(3(/x,a,t73,i8,/),/x,a,t73,f8.6,/,/x,a,t73,i8,/)
      if (status .ne. 0) then
         fcc_extract = %loc(fcc_repwrite)
         inquire(fut_report_lun,name=filename)
         call str$trim(filename,filename,len)
         call lib$signal(fcc_repwrite,%val(2),filename(:len),%val(status))
         return
      end if
c
c     Sum the norms array over all bins.
c
      sum_norm = 0.0
      do rec = 2,257
         sum_norm = sum_norm + in_cov(rec).norms(1)
      end do
c
c     Define the average norm to be this sum divided by the number of ifgs.
c
      avg_norm = sum_norm/num_ifgs
c
c     Define the effective weight at a bin to be the norms at that bin divided
c     by the bin total at that bin times the average norm.
c
      do rec = 2,257
         if (in_cov(rec).bin_total .eq. 0) then
            out_cov(rec).eff_wt = 0.0
         else
            out_cov(rec).eff_wt = in_cov(rec).norms(1) / 
     .                            (in_cov(rec).bin_total * avg_norm)
         end if
      end do
c
c     Reconstruct the 522 by 522 FIC covariance matrix from the input records;
c     only the upper right half of the symmetric matrix is stored.  The 
c     real-real portion of the covariance matrix is in records 2-66 with 550 
c     values in records 2-65 and 45 values in record 66 for a total of 35245 
c     values.
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
            c(row,col) = in_cov(rec).covar(pos)
         end do
      end do
c
c     The real-imaginary portion of the covariance matrix is in records 67-190 
c     with 550 values in records 67-189 and 455 values in record 190 for a 
c     total of 68105 values.
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
            c(row,col) = in_cov(rec).covar(pos)
         end do
      end do
c
c     The imaginary-imaginary portion of the covariance matrix is in records 
c     191-251 with 550 values in records 191-250 and 153 values in record 251
c     for a total of 33153 values.
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
            c(row,col) = in_cov(rec).covar(pos)
         end do
      end do
c
c     Construct the lower left half of C by symmetry from the retrieved
c     upper right half.
c
      do row = 2,522
         do col = 1,row-1
            c(row,col) = c(col,row)
         end do
      end do
c
c     Construct the intermediate matrices N and N2 from the reconstructed
c     FIC covariance matrix C.
c
      do row = 1,257
         do col = 1,257
            n(row,col) = dcmplx(c(row,col)+c(row+265,col+265),
     .                          c(row+265,col)-c(col+265,row))
            n2(row,col) = dcmplx(c(row,col)-c(row+265,col+265),
     .                           c(row+265,col)+c(col+265,row))
         end do
         do col = 258,512
            n(row,col) = n2(row,514-col)
            n2(row,col) = n(row,514-col)
         end do
      end do
      do row = 258,512
         do col = 1,257
            n(row,col) = conjg(n(col,row))
            n2(row,col) = n2(col,row)
         end do
         do col = 258,512
            n(row,col) = n(514-col,514-row)
            n2(row,col) = conjg(n2(514-row,514-col))
         end do
      end do
c
c     Construct the intermediate matrix EQ from the reconstructed FIC 
c     covariance matrix C.
c
      do row = 1,8
         do col = 1,257
            eq(row,col) = dcmplx(c(row+257,col),-c(row+257,col+265))
         end do
         do col = 258,512
            eq(row,col) = conjg(eq(row,514-col))
         end do
      end do
c
c     Construct the intermediate matrix DA from the binned totals of voltage
c     spectra dispersions and the real and imaginary parts of the coadd 
c     voltage spectra.
c
      do row = 1,256
         do col = 1,257
            da(row,col) = in_cov(row+1).disp(col)
         end do
         do col = 258,512
            da(row,col) = conjg(da(row,514-col))
         end do
      end do
      do col = 1,257
         da(257,col) = dcmplx(in_cov(1).avg_rvoltspec(col), 
     .                        in_cov(1).avg_ivoltspec(col))
      end do
      do col = 258,512
         da(257,col) = conjg(da(257,514-col))
      end do

      return
      end
