      integer*4 function fic_sec_template(input,channel,scan_mode,
     .                                    sec_template_input,cth,reftemps_lun,
     .                                    num_recs,num_good,xcal_pos,con_recs,
     .                                    coa_rec,aifgs,template,sec_template,
     .                                    num_peak_pts,x,b)
c-------------------------------------------------------------------------------
c
c     Purpose: Calculate and subtract the secondary template to remove signal
c              contributed by dust and gas in the galactic plane.
c
c     Authors: S. Read, STX, 6/90
c              S. Alexander, STX, 9/91, SER 7985
c
c     Input: input              i*4  Type of input: sky, calibration, or raw.
c            channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode          i*4  Value of scan mode, 1-4=SS-LF.
c            sec_template_input i*4  Perform secondary template subtraction.
c            cth                rec  Consistency check parameters.
c            reftemps_lun       i*4  Logical unit for reference templates.
c            num_recs           i*4  Number of input science and engineering
c                                    records.
c            num_good           i*4  Current number of good interferograms.
c            xcal_pos           i*4  External calibrator position,
c                                    1=in=calibration data, 2=out=sky data.
c            con_recs           rec(num_recs)  Consistency check records.
c            coa_rec            rec  Coadd record.
c            aifgs              r*4(512,num_recs)  Real-valued interferograms.
c            template           r*4(512)  Primary template.
c            num_peak_pts       i*4  Number of points around interferogram
c                                    peak for secondary template fitting.
c            x                  r*4(num_peak_pts,2)  Fitting matrix for
c                                                    secondary template 
c                                                    calculation.
c
c     Output: con_recs           rec(num_recs)  Consistency check records.
c             coa_rec            rec  Coadd record.
c             aifgs              r*4(512,num_recs)  Real-valued interferograms.
c             sec_template       r*4(512)  Secondary template.
c             b                  r*4(num_recs)  Secondary template coefficients.
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
      external fic_normal
      external fic_nosectemp
      external fic_cfdirread
      external fic_primrms
      external fic_secrms
c
c     Input parameters.
c
      integer*4 input,channel,scan_mode,sec_template_input

      dictionary 'fex_cth'
      record /fex_cth/ cth

      integer*4 reftemps_lun,num_recs,num_good,xcal_pos,num_peak_pts
      real*4 x(num_peak_pts,2),template(512)
c
c     Input/output parameters.
c
      dictionary 'fic_scc'
      record /fic_scc/ con_recs(num_recs)

      dictionary 'fic_sky'
      record /fic_sky/ coa_rec

      real*4 aifgs(512,num_recs)
c
c     Output parameters.
c
      real*4 sec_template(512),b(num_recs)
c
c     Local variables.
c
      integer*4 peak_pos,rec,posneg(fac_max_num)
      integer*4 num_neg,ifg_pos
      integer*4 first,last,num_midav_pts
      real*4 sum,ave,sort(512)
      integer*4 rec_good,mid
      integer*4 offset,reftemps_rec,status

      dictionary 'fex_reftemps'
      record /fex_reftemps/ reftemps

      character*72 filename
      integer*4 file_len,index
      real*4 amp,prim_sum,sec_sum,prim_ave,sec_ave,sae
      integer*4 rank,iter,nrmiss

      fic_sec_template = %loc(fic_normal)

      peak_pos = coa_rec.coad_spec_data.peak_pos
c
c     Check the peak position before proceeding; it may have been set to a
c     flag value in fic_qual_check or may be too small or large to compute
c     template statistics.
c
      if ((peak_pos .lt. 51) .or. (peak_pos .gt. 462)) then
         call lib$signal(fic_nosectemp,%val(2),fac_channel_ids(channel),
     .                   coa_rec.ct_head.gmt)
         return
      end if
c
c     Determine whether each peak is positive or negative.
c
      do rec = 1,num_recs
         posneg(rec) = 1
         if (con_recs(rec).con_check .eq. 0) then
            num_neg = 0
            do ifg_pos = peak_pos-1,peak_pos+1
               if (aifgs(ifg_pos,rec) .lt. 0.0) num_neg = num_neg + 1
            end do
            if (num_neg .ge. 2) posneg(rec) = -1
         end if
      end do
c
c     Form the secondary template by sorting the coadd group to find the
c     midaverage of the group at each of the 512 interferogram points.  If the
c     peak position was negative, multiply by -1 to reverse the interferogram.
c     Retrieve the upper and lower bounds for the midaverage from the cth
c     reference data set.
c
      first = nint(num_good*cth.lower_bounds(channel))
      last = nint(num_good*cth.upper_bounds(channel))
      num_midav_pts = last - first

      sum = 0.0

      do ifg_pos = 1,512
         rec_good = 0
         do rec = 1,num_recs
            if (con_recs(rec).con_check .eq. 0) then
               rec_good = rec_good + 1
               sort(rec_good) = aifgs(ifg_pos,rec)*posneg(rec)
            end if
         end do
         call svrgn(num_good,sort,sort)
         sec_template(ifg_pos) = 0.0
         do mid = first+1,last
            sec_template(ifg_pos) = sec_template(ifg_pos) + sort(mid)
         end do
         sec_template(ifg_pos) = sec_template(ifg_pos)/num_midav_pts
         sum = sum + sec_template(ifg_pos)
      end do
c
c     Calculate the variance of the secondary template.
c
      ave = sum/512
      sum = 0.0

      do ifg_pos = 1,512
         sum = sum + (sec_template(ifg_pos) - ave)**2
      end do

      coa_rec.coad_spec_data.sec_template.variance = sum/511
c
c     Retrieve the appropriate reference primary and secondary templates from 
c     the fex_reftemps reference data set.
c
      if ((input .eq. fac_cal) .or. 
     .    ((input .eq. fac_raw) .and. (xcal_pos .eq. fac_xcalin))) then
         offset = 4
      else
         offset = 0
      end if

      reftemps_rec = 8*(channel-1) + scan_mode + offset

      read(reftemps_lun,rec=reftemps_rec,iostat=status) reftemps
      if (status .ne. 0) then
         fic_sec_template = %loc(fic_cfdirread)
         inquire(reftemps_lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fic_cfdirread,%val(2),filename(:file_len),%val(status))
         return
      end if
c
c     Compute the primary and secondary template amplitudes by finding the inner
c     product of each template and its corresponding reference template.
c
      amp = 0.0
      index = 0
c
c     The 41 points around the peak of the primary template are used.
c
      do ifg_pos = peak_pos-20,peak_pos+20
         index = index + 1
         amp = amp + template(ifg_pos)*reftemps.prim_temp(index)
      end do
 
      coa_rec.coad_spec_data.prim_template.amplitude = abs(amp)

      amp = 0.0
      index = 0
c
c     The three points around the peak of the secondary template are used.
c
      do ifg_pos = peak_pos-1,peak_pos+1
         index = index + 1
         amp = amp + sec_template(ifg_pos)*reftemps.sec_temp(index)
      end do

      coa_rec.coad_spec_data.sec_template.amplitude = abs(amp)
c
c     Compute the off-peak root-mean-square of the templates.
c
      prim_sum = 0.0
      sec_sum = 0.0

      do ifg_pos = 21,peak_pos-50
         prim_sum = prim_sum + template(ifg_pos)
         sec_sum = sec_sum + sec_template(ifg_pos)
      end do
      do ifg_pos = peak_pos+50,512
         prim_sum = prim_sum + template(ifg_pos)
         sec_sum = sec_sum + sec_template(ifg_pos)
      end do

      prim_ave = prim_sum/393
      sec_ave = sec_sum/393

      prim_sum = 0.0
      sec_sum = 0.0

      do ifg_pos = 21,peak_pos-50
         prim_sum = prim_sum + (template(ifg_pos) - prim_ave)**2      
         sec_sum = sec_sum + (sec_template(ifg_pos) - sec_ave)**2      
      end do
      do ifg_pos = peak_pos+50,512
         prim_sum = prim_sum + (template(ifg_pos) - prim_ave)**2      
         sec_sum = sec_sum + (sec_template(ifg_pos) - sec_ave)**2      
      end do
c
c     If the root-mean-square is zero, signal this information.  Otherwise, set
c     the signal-to-noise ratio to be the amplitude divided by the
c     root-mean-square value.
c
      if (prim_sum .eq. 0.0) then
         call lib$signal(fic_primrms,%val(2),fac_channel_ids(channel),
     .                   coa_rec.ct_head.gmt)
         coa_rec.coad_spec_data.prim_template.snr = 0.0
      else 
         coa_rec.coad_spec_data.prim_template.snr = 
     .           coa_rec.coad_spec_data.prim_template.amplitude / 
     .           sqrt(prim_sum/392)
      end if

      if (sec_sum .eq. 0.0) then
         call lib$signal(fic_secrms,%val(2),fac_channel_ids(channel),
     .                   coa_rec.ct_head.gmt)
         coa_rec.coad_spec_data.sec_template.snr = 0.0
      else 
         coa_rec.coad_spec_data.sec_template.snr =
     .           coa_rec.coad_spec_data.sec_template.amplitude /
     .           sqrt(sec_sum/392)
      end if
c
c     Decide whether or not the secondary template should be subtracted based
c     on the amplitudes, the root-mean-squares, and values in the cth reference
c     data set.
c
      if ((sec_template_input .eq. fac_present) .and.
     .    (coa_rec.coad_spec_data.prim_template.amplitude .gt.
     .     cth.prim_temp_amp(channel,scan_mode)) .and.
     .    (coa_rec.coad_spec_data.prim_template.snr .gt.
     .     cth.prim_temp_snr(channel,scan_mode)) .and.
     .    (coa_rec.coad_spec_data.sec_template.amplitude .gt.
     .     cth.sec_temp_amp(channel,scan_mode)) .and.
     .    (coa_rec.coad_spec_data.sec_template.snr .gt.
     .     cth.sec_temp_snr(channel,scan_mode))) then
         coa_rec.coad_spec_data.sec_template.subtracted = fac_present
      else
         coa_rec.coad_spec_data.sec_template.subtracted = fac_not_present
      end if
c
c     Compute a coefficient for each interferogram which minimizes the 
c     difference between the interferogram and the secondary template around 
c     the interferogram peak position.
c                                     
      sum = 0.0
      do rec = 1, num_recs
         if (con_recs(rec).con_check .eq. 0) then
            ifg_pos = peak_pos - (num_peak_pts/2 + 1)
            do index = 1, num_peak_pts
               ifg_pos = ifg_pos + 1
               x(index,1) = sec_template(ifg_pos)   
               x(index,2) = aifgs(ifg_pos,rec)
            end do
c
c     Parameters to the imsl routine rlav are: number of points around peak 
c     (number of observations); number of columns (2); matrix of data with 
c     independent variable in first column and dependent variable in second; 
c     leading dimension of matrix; no intercept term (0); negative of number of
c     independent variables (-1); index vector of columns of independent 
c     variable (not used); column number of dependent variable (2); estimated 
c     coefficients; rank of matrix of regressors; sum of absolute values of 
c     errors; number of iterations required for solution; number of rows of 
c     data containing invalid numbers.
c
            call rlav(num_peak_pts,2,x,num_peak_pts,0,-1,index,2,b(rec),rank,
     .                 sae,iter,nrmiss)
            sum = sum + b(rec)
            con_recs(rec).b = b(rec)
c
c     Unless disabled on the command line or by prior calculations, subtract 
c     the secondary template times the calculated coefficient from the
c     interferogram.
c
            if (coa_rec.coad_spec_data.sec_template.subtracted .eq. 
     .          fac_present) then
               do ifg_pos = 1, 512
                  aifgs(ifg_pos,rec) = aifgs(ifg_pos,rec) - 
     .                                 b(rec)*sec_template(ifg_pos)
               end do
            end if
         end if
      end do
c
c     Calculate the average and variance of the "b" secondary template
c     coefficients.
c
      ave = sum/num_good
      coa_rec.coad_spec_data.sec_template.b_average = ave
      sum = 0.0
      do rec = 1,num_recs
         if (con_recs(rec).con_check .eq. 0) then
            sum = sum + (b(rec) - ave)**2
         end if
      end do
      coa_rec.coad_spec_data.sec_template.b_variance = sum/(num_good-1)

      return
      end
