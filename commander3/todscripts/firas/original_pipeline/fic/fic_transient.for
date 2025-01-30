      integer*4 function fic_transient(channel,scan_mode,dtrf_lun,num_recs,
     .                                 num_good,con_recs,coa_rec,aifgs,
     .                                 transient_coadd,c)
c-----------------------------------------------------------------------------
c
c     Purpose: Retrieve and subtract the least squares projection of the
c              digital transient response function.
c
c     Authors: S. Alexander, STX, 9/91, SER 7985
c
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode          i*4  Value of scan mode, 1-4=SS-LF.
c            dtrf_lun           i*4  Logical unit for digital transient
c                                    response functions.
c            num_recs           i*4  Number of input science and engineering
c                                    records.
c            num_good           i*4  Current number of good interferograms.
c            con_recs           rec(num_recs)  Consistency check records.
c            coa_rec            rec  Coadd record.
c            aifgs              r*4(512,num_recs)  Real-valued interferograms.
c            transient_coadd    l*1  Flag for digital transient subtraction
c                                    from a coadded interferogram.
c
c     Output: coa_rec            rec  Coadd record.
c             aifgs              r*4(512,num_recs)  Real-valued interferograms.
c             c                  r*4(num_recs)  Digital transient coefficients.
c   
c     Modifications:
c
c-----------------------------------------------------------------------------

      implicit none
c
c     Return statuses.
c
      external fic_normal
      external fic_cfdirread
c
c     Input parameters.
c
      integer*4 channel,scan_mode,dtrf_lun
      integer*4 num_recs,num_good

      dictionary 'fic_scc'

      record /fic_scc/ con_recs(num_recs)

      logical*1 transient_coadd
c
c     Input/output parameters.
c
      dictionary 'fic_sky'
      record /fic_sky/ coa_rec

      real*4 aifgs(512,num_recs)
c
c     Output parameters.
c
      real*4 c(num_recs)
c
c     Local variables.
c
      integer*4 dtrf_rec,status

      dictionary 'fex_dtrf'
      record /fex_dtrf/ dtrf

      character*72 filename
      integer*4 file_len
      integer*4 ifg_pos,rec
      real*4 sum,ave

      fic_transient = %loc(fic_normal)
c
c     Retrieve the appropriate digital transient response function.
c
      dtrf_rec = 4*mod((channel-1),2) + scan_mode

      read(dtrf_lun,rec=dtrf_rec,iostat=status) dtrf
      if (status .ne. 0) then
         fic_transient = %loc(fic_cfdirread)
         inquire(dtrf_lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fic_cfdirread,%val(2),filename(:file_len),%val(status))
         return
      end if
c
c     If removing the transient from the coadded interferogram, put the 
c     interferogram into the aifgs array.
c
      if (transient_coadd) then
         do ifg_pos = 1,128
            aifgs(ifg_pos,1) = coa_rec.coad_data.ifg(ifg_pos)
         end do
      end if
c
c     Determine the coefficient for the least squares projection of the
c     interferogram onto the response function.  The norm of the digital 
c     transient response function vector is unity.
c
      sum = 0.0
      do rec = 1,num_recs
         if ((con_recs(rec).con_check .eq. 0) .or. transient_coadd) then
            c(rec) = 0.0
            do ifg_pos = 1,128
               c(rec) = c(rec) + aifgs(ifg_pos,rec) * dtrf.trans(ifg_pos)
            end do
            sum = sum + c(rec)
c
c     If processing a coadd group, store the "c" coefficients in the
c     consistency check records.
c
            if (.not. transient_coadd) then
               con_recs(rec).c = c(rec)
            end if
c
c     Subtract the 128 points of the digital transient response function times
c     the coefficient from the interferogram.
c
            do ifg_pos = 1,128
               aifgs(ifg_pos,rec) = aifgs(ifg_pos,rec) - 
     .                              c(rec) * dtrf.trans(ifg_pos)
            end do
         end if
      end do
c
c     If processing the coadded interferogram, put the interferogram back into
c     the coadd record.
c
      if (transient_coadd) then
         do ifg_pos = 1,128
            coa_rec.coad_data.ifg(ifg_pos) = aifgs(ifg_pos,1)
         end do
         coa_rec.coad_spec_data.transient.bl_trans_coeff = sum
      else
c
c     If processing a coadd group, calculate the average and variance of the
c     "c" coefficients.
c
         ave = sum/num_good
         coa_rec.coad_spec_data.transient.c_average = ave
         sum = 0.0
         do rec = 1,num_recs
            if (con_recs(rec).con_check .eq. 0) then
                  sum = sum + (c(rec) - ave)**2
            end if
         end do
c
c     If processing single interferograms, set a flag value for the variance.
c
         if (num_good .eq. 1) then
            coa_rec.coad_spec_data.transient.c_variance = -9999.0
         else
            coa_rec.coad_spec_data.transient.c_variance = sum/(num_good-1)
         end if
      end if

      return
      end
