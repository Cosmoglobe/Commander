      integer*4 function fut_setxaxl(channel,scan_mode,fake_it,sci_mode,
     .                               adds_per_group,nyquistl,scale,first_samp,
     .                               startx,spacex)
c------------------------------------------------------------------------------
c
c     Purpose: Set x-axis parameters for plotting.
c
c     Authors: W.K. Young, SASC, 8/85
c              A.R. Trenholme, GSC, 2/95
c
c     Input: channel         i*4  Value of channel (1-4 = RH-LL).
c            scan_mode       i*4  Value of scan mode (1-6 = SS-FL).
c            fake_it         i*4  Value of fake-it bit.
c            sci_mode        i*4  Value of science mode (0-4).
c            adds_per_group  i*4  Value of adds per group (1-12).
c            nyquistl        rec  Nyquist frequencies.
c            scale           i*4  Type of x-axis scale (0=coadd, 1=spectrum,
c                                  -1=noisetest spectrum)
c            first_samp      i*4  Start sample number, usually 1.
c
c     Output: startx          r*4  X-axis start point in centimeters.
c             spacex          r*4  X-axis space interval in centimeters.
c                         
c     Modifications: S. Brodd, HSTX, SPR 12288.  Version for FSL.
c
c---------------------------------------------------------------------------
      implicit none
c
c     Include files.
c
      include '(fut_params)'
c
c     Return status.
c
      external fut_normal
c
c     Functions called.
c
      integer*4 fut_default_peak
c
c     Input parameters.
c
      integer*4 channel,scan_mode,fake_it,sci_mode
      integer*4 adds_per_group,scale,first_samp

      dictionary 'fex_nyquistl'
      record /fex_nyquistl/ nyquistl
c
c     Output parameters.
c
      real*4 startx,spacex
c
c     Local parameters.
c
      integer*4 status,peak_pos,mtm_speed,mtm_length,index
      real*4 nyq_icm,nyq_hz

      fut_setxaxl = %loc(fut_normal)
c
c     Retrieve the peak position, using the unlinearized position for plotting.
c
      if (scan_mode .ge. 5) then
         status = fut_default_peak(1,1,channel,8,sci_mode,0,peak_pos)
      else
         mtm_speed = mod(scan_mode-1,2)
         mtm_length = int((scan_mode-1)/2)
         status = fut_default_peak(mtm_speed,mtm_length,channel,adds_per_group,
     .                             sci_mode,0,peak_pos)
      end if
c
c     Retrieve the nyquist frequency in icm and hertz.
c
      if (fake_it .eq. 1) then
         if (adds_per_group .eq. 1) then
            nyq_hz = nyquistl.hz(11)
         else if (adds_per_group .eq. 2) then
            nyq_hz = nyquistl.hz(12)
         else if (adds_per_group .eq. 3) then
            nyq_hz = nyquistl.hz(13)
         else if (adds_per_group .eq. 8) then
            nyq_hz = nyquistl.hz(14)
         else if (adds_per_group .eq. 12) then
            nyq_hz = nyquistl.hz(15)
         end if
      else
         index = 4*(mod(channel-1,2)) + scan_mode
         nyq_icm = nyquistl.icm(index)
         nyq_hz =  nyquistl.hz(index)
      end if         

      if (scale .eq. 0) then
c
c     Find spacex and startx for coadd plotting.
c
         if (fake_it .eq. 1) then
            spacex = 1./512
            startx = 0.
         else
            spacex = 0.5/nyq_icm
            startx = ((first_samp-1.)/adds_per_group - peak_pos + 1)*spacex
         end if
      else if (scale .eq. 1) then
c
c     Find spacex and startx for spectrum plotting.
c
         spacex = nyq_icm/(fac_fft_length(scan_mode)/2)
         startx = 0.
      else
c
c     Find spacex and startx for noisetest spectrum plotting.
c
         spacex = nyq_hz/(fac_fft_length(scan_mode)/2)
         startx = 0.
      end if

      return
      end
