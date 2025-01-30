      integer*4 function fil_short(channel,scan_mode,nyquistl,sci_mode,
     .                             in_coa_rec,var_vecs,out_coa_rec)
c----------------------------------------------------------------------------
c
c     Purpose: Change coadd record fields as necessary to convert SF or LF
c              coadds into corresponding FS or FL coadds.
c
c     Author: S. Brodd, HSTX, 4/95
c
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode          i*4  Short scan mode, 5 = FS, 6 = FL.
c            nyquistl           rec  Nyquist frequencies.
c            sci_mode           i*4  Value of science mode, 0-4.
c            in_coa_rec         rec  Input coadd record.
c            var_vecs           r*4(3,361)  Variance vectors for FS or FL.
c
c     Output: out_coa_rec        rec  Modified output coadd record.
c
c     Modifications:
c
c----------------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_params)'
c
c     Return status.
c
      external fil_normal
c
c     Functions.
c
      integer*4 fut_default_peak
c
c     Input parameters.
c
      integer*4 channel,scan_mode,sci_mode

      dictionary 'fex_nyquistl'
      dictionary 'fil_sky'

      record /fex_nyquistl/ nyquistl
      record /fil_sky/ in_coa_rec

      real*4 var_vecs(3,361)
c
c     Output parameters.
c
      record /fil_sky/ out_coa_rec
c
c     Local variables.
c
      integer*4 ifg_pos,pos,index,status
      character*1 xcal_pos

      fil_short = %loc(fil_normal)
c
c     Begin by copying input coadd into output coadd.
c
      out_coa_rec = in_coa_rec
c
c     For FS, change the adds per group to 8, find a new peak position, and
c     decimate the coadded interferogram.
c
      if (scan_mode .eq. 5) then
         out_coa_rec.coad_spec_data.adds_per_group = 8
         status = fut_default_peak(1,1,channel,8,sci_mode,0,
     .                             out_coa_rec.coad_spec_data.peak_pos)
         do ifg_pos = 1,127
            pos = ifg_pos * 4
            out_coa_rec.coad_data.ifg(ifg_pos) = 
     .          (in_coa_rec.coad_data.ifg(pos) + 
     .           in_coa_rec.coad_data.ifg(pos+1) +
     .           in_coa_rec.coad_data.ifg(pos+2) + 
     .           in_coa_rec.coad_data.ifg(pos+3)) / 4.0
         end do
         out_coa_rec.coad_data.ifg(128) = in_coa_rec.coad_data.ifg(512)
      else
c
c     For FL, change the mtm length to 0.
c
         out_coa_rec.coad_spec_data.mtm_length = 0
      end if
c
c     For both FS and FL, change the coadd label, find the new nyquist
c     frequencies in hertz and icm, set the points of the coadded interferogram
c     from 129 to 512 to zero, find the new fft length, and replace the
c     variance vectors with the correct values.
c
      if (in_coa_rec.coad_spec_data.xcal_pos .eq. fac_xcalout) then
         xcal_pos = 'O'
      else
         xcal_pos = 'I'
      end if

      if (in_coa_rec.attitude.pixel_no .ge. 0) then
         encode (60,10,out_coa_rec.coad_spec_head.label)
     .           in_coa_rec.coad_spec_head.first_gmt,
     .           in_coa_rec.coad_spec_head.last_gmt,
     .           in_coa_rec.coad_spec_head.num_ifgs,
     .           fac_scan_mode_idsl(scan_mode),
     .           nint(in_coa_rec.coad_spec_data.ical*1000.0),
     .           xcal_pos,in_coa_rec.attitude.pixel_no
10       format (2(a,'_'),i3.3,'_',a,'_',i5.5,'_',a,'_',i4.4,11x)
      else
         encode (60,20,out_coa_rec.coad_spec_head.label)
     .           in_coa_rec.coad_spec_head.first_gmt,
     .           in_coa_rec.coad_spec_head.last_gmt,
     .           in_coa_rec.coad_spec_head.num_ifgs,
     .           fac_scan_mode_idsl(scan_mode),
     .           nint(in_coa_rec.coad_spec_data.ical*1000.0),
     .           xcal_pos
20       format (2(a,'_'),i3.3,'_',a,'_',i5.5,'_',a,16x)
      end if

      do ifg_pos = 129,512
         out_coa_rec.coad_data.ifg(ifg_pos) = 0.0
      end do

      index = 4*(mod(channel-1,2)) + scan_mode
      out_coa_rec.coad_spec_data.nyquist_hertz = nyquistl.hz(index)
      out_coa_rec.coad_spec_data.nyquist_icm = nyquistl.icm(index)

      out_coa_rec.coad_data.fft_length = fac_fft_length(scan_mode)

      do pos = 1,361
         out_coa_rec.coad_data.real_var(pos) = var_vecs(1,pos)
         out_coa_rec.coad_data.imag_var(pos) = var_vecs(2,pos)
         out_coa_rec.coad_data.real_imag_var(pos) = var_vecs(3,pos)
      end do

      return
      end
