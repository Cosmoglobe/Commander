      integer*4 function fcc_reference(channel,scan_mode,sci_mode,
     .                                 adds_per_group,start,nyquist_hz,model,
     .                                 apod,etf,peak)
c-------------------------------------------------------------------------------
c
c     Purpose: Retrieve all reference data.
c
c     Author: S. Alexander, HSTX, 7/93, SER 11189
c
c     Input: channel         i*4  Channel value, 1-4.
c            scan_mode       i*4  Scan mode, 1,2, or 4.
c            sci_mode        i*4  Science mode, 0-4.
c            adds_per_group  i*4  Adds per group value, 1-12.
c            start           i*4(2)  Calibration matrix start time tag.
c
c     Output: nyquist_hz      r*4  Nyquist frequency in Hertz.
c             model           rec  Calibration model.
c             apod            r*8(512)  Apodization function.
c             etf             c*16(512) Extended electronics transfer function.
c             peak            i*4  Peak location.
c      
c     Modifications:
c    
c-------------------------------------------------------------------------------
      implicit none
c
c     Include files.
c
      include '(fut_params)'
      include '(cct_get_config)'
c
c     Return statuses.
c 
      external fcc_normal
      external fcc_cfseqopen,fcc_cfseqget,fcc_cfseqclose
      external fcc_cfdiropen,fcc_cfdirget,fcc_cfdirread,fcc_cfdirclose
c
c     Functions.
c
      integer*4 cct_open_config
      integer*4 cct_get_config_tod
      integer*4 cct_get_config_idx_tod
      integer*4 cct_close_config      
      integer*4 fut_default_peak
c
c     Input parameters.
c
      integer*4 channel,scan_mode,sci_mode,adds_per_group,start(2)
c
c     Output parameters.
c
      real*4 nyquist_hz

      dictionary 'fex_mod'
      record /fex_mod/ model

      real*8 apod(512)
      complex*16 etf(512)
      integer*4 peak
c
c     Local variables.
c
      integer*4 status

      character*14 config_gmt_start,config_gmt_stop,jstart
      integer*4 config_start(2),config_stop(2)

      character*28 config_name(2),filename
      integer*4 config_size(2),config_lun(2),config_index(2)

      record /config_status/ config_status(2)

      integer*4 config_ref_count

      dictionary 'fex_nyquist'
      structure /config_record/
         record /fex_nyquist/ nyquist
         record /fex_mod/ model
      endstructure        

      record /config_record/ config_record

      logical*1 config_new_segment(2)

      integer*4 freq,speed,length,mode,rec,len,pos

      dictionary 'fex_apod'
      record /fex_apod/ apod_rec

      dictionary 'fex_etf'
      record /fex_etf/ etf_rec

      fcc_reference = %loc(fcc_normal)
c     
c     Assign the sequential read configuration files (those with only one
c     record).  The start time must be advanced by one second and the stop
c     time decreased, or the open will fail.
c
      config_gmt_start = '86001000001000'
      call ct_gmt_to_binary(config_gmt_start,config_start)
      config_gmt_stop  = '99365235958990'
      call ct_gmt_to_binary(config_gmt_stop,config_stop)

      config_name(1) = 'CSDR$FIRAS_REF:FEX_NYQUIST'
      config_name(2) = 'CSDR$FIRAS_REF:FEX_MOD_'//fac_channel_ids(channel)//
     .                                            fac_scan_mode_ids(scan_mode)
     
      config_size(1) = 128
      config_size(2) = 29696 
c
c     Open the sequential access reference data sets.
c
      status = cct_open_config(config_start,config_stop,2,config_name,
     .                         config_size,' ',1,config_lun,config_index,
     .                         config_status,config_ref_count)  
      if (.not. status) then
         fcc_reference = %loc(fcc_cfseqopen)
         call lib$signal(fcc_cfseqopen,%val(1),%val(status))
         return
      end if
c
c     Get the reference data sets at the specified start time.
c
      status = cct_get_config_tod(start,2,config_size,config_lun,config_index,
     .                            config_record,config_new_segment,
     .                            config_status)
      if (.not. status) then
         fcc_reference = %loc(fcc_cfseqget)
         call ct_binary_to_gmt(start,jstart)
         call lib$signal(fcc_cfseqget,%val(2),jstart,%val(status))
         return
      end if
c
c     Get nyquist frequency and calibration model.
c
      if ((channel .eq. 1) .or. (channel .eq. 3)) then
         nyquist_hz = config_record.nyquist.hz(scan_mode)
      else
         nyquist_hz = config_record.nyquist.hz(scan_mode+4)
      end if

      model = config_record.model
c
c     Close sequential access configuration files.
c
      status = cct_close_config(2,config_lun,config_index)
      if (.not. status) then
         fcc_reference = %loc(fcc_cfseqclose)
         call lib$signal(fcc_cfseqclose,%val(1),%val(status))
         return
      end if
c
c     Assign direct read configuration files (those with multiple records).
c 
      config_name(1) = 'CSDR$FIRAS_REF:FEX_APOD'
      config_name(2) = 'CSDR$FIRAS_REF:FEX_ETF'
                     
      config_size(1) = 4096
      config_size(2) = 4112
c
c     Open the direct access reference data sets.
c
      status = cct_open_config(config_start,config_stop,2,config_name,
     .                         config_size,'direct',1,config_lun,config_index,
     .                         config_status,config_ref_count)
      if (.not. status) then
         fcc_reference = %loc(fcc_cfdiropen)
         call lib$signal(fcc_cfdiropen,%val(1),%val(status))
         return
      end if
c
c     Get the direct access reference data sets at the specified start time.
c
      status = cct_get_config_idx_tod(start,2,config_lun,config_index,
     .                                config_new_segment,config_status)
      if (.not. status) then
         fcc_reference = %loc(fcc_cfdirget)
         call ct_binary_to_gmt(start,jstart)
         call lib$signal(fcc_cfdirget,%val(2),jstart,%val(status))
         return
      end if
c
c     Determine the correct record for the apodization function.
c
      if ((channel .eq. 2) .or. (channel .eq. 4)) then
         freq = 0
      else
         freq = 1
      end if
      if (scan_mode .eq. 1) then
         speed = 0
      else
         speed = 1
      end if
      if (scan_mode .eq. 4) then
         length = 1
      else
         length = 0
      end if
      if ((sci_mode .eq. 2) .or. (sci_mode .eq. 4)) then
         mode = 0
      else
         mode = 1
      end if

      rec = 32*(adds_per_group-1) + 16*speed + 8*length + 4*freq +
     .      2*abs(mode-1) + 2
c
c     Read the apodization function.
c
      read(config_lun(1),rec=rec,iostat=status) apod_rec
      if (status .ne. 0) then
         fcc_reference = %loc(fcc_cfdirread)
         call str$trim(config_name(1),filename,len)
         call lib$signal(fcc_cfdirread,%val(2),filename(:len),%val(status))
         return
      end if

      do pos = 1,512
         apod(pos) = apod_rec.apodfcn(pos)
      end do
c
c     Determine the correct record for the electronics transfer function.
c
      rec = 96*speed + 24*(channel-1) + 12*mode + adds_per_group
c
c     Read the electronics transfer function.
c
      read(config_lun(2),rec=rec,iostat=status) etf_rec
      if (status .ne. 0) then
         fcc_reference = %loc(fcc_cfdirread)
         call str$trim(config_name(2),filename,len)
         call lib$signal(fcc_cfdirread,%val(2),filename(:len),%val(status))
         return
      end if

      do pos = 1,257
         etf(pos) = etf_rec.ztransfcn(pos)
      end do
      do pos = 258,512
         etf(pos) = conjg(etf(514-pos))
      end do
c
c     Close the direct access reference data sets.
c
      status = cct_close_config(2,config_lun,config_index)
      if (.not. status) then
         fcc_reference = %loc(fcc_cfdirclose)
         call lib$signal(fcc_cfdirclose,%val(1),%val(status))
         return
      end if
c
c     Get the peak position.
c
      status = fut_default_peak(speed,length,channel,adds_per_group,sci_mode,1,
     .                          peak)

      return
      end
