      integer*4 function fcl_reference(channel,scan_mode,start,nyquistl_hz,
     .                                 model)
c-------------------------------------------------------------------------------
c
c     Purpose: Retrieve all reference data.
c
c     Author: S. Brodd, HSTX, 12/95, SPR 12291
c
c     Input: channel         i*4  Channel value, 1-4.
c            scan_mode       i*4  Scan mode, 1,2,4,5, or 6.
c            start           i*4(2)  Calibration matrix start time tag.
c
c     Output: nyquistl_hz     r*4  Nyquist frequency in Hertz.
c             model           rec  Calibration model.
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
      external fcl_normal
      external fcl_cfseqopen,fcl_cfseqget,fcl_cfseqclose
      external fcl_cfdiropen,fcl_cfdirget,fcl_cfdirread,fcl_cfdirclose
c
c     Functions.
c
      integer*4 cct_open_config
      integer*4 cct_get_config_tod
      integer*4 cct_get_config_idx_tod
      integer*4 cct_close_config      
c
c     Input parameters.
c
      integer*4 channel,scan_mode,start(2)
c
c     Output parameters.
c
      real*4 nyquistl_hz

      dictionary 'fex_mod'
      record /fex_mod/ model
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

      dictionary 'fex_nyquistl'
      structure /config_record/
         record /fex_nyquistl/ nyquistl
         record /fex_mod/ model
      endstructure        

      record /config_record/ config_record

      logical*1 config_new_segment(2)

      fcl_reference = %loc(fcl_normal)
c     
c     Assign the sequential read configuration files (those with only one
c     record).  The start time must be advanced by one second and the stop
c     time decreased, or the open will fail.
c
      config_gmt_start = '86001000001000'
      call ct_gmt_to_binary(config_gmt_start,config_start)
      config_gmt_stop  = '99365235958990'
      call ct_gmt_to_binary(config_gmt_stop,config_stop)

      config_name(1) = 'CSDR$FIRAS_REF:FEX_NYQUISTL'
      config_name(2) = 'CSDR$FIRAS_REF:FEX_MOD_'//fac_channel_ids(channel)//
     .                                            fac_scan_mode_idsl(scan_mode)
     
      config_size(1) = 128
      config_size(2) = 29696 
c
c     Open the sequential access reference data sets.
c
      status = cct_open_config(config_start,config_stop,2,config_name,
     .                         config_size,' ',1,config_lun,config_index,
     .                         config_status,config_ref_count)  
      if (.not. status) then
         fcl_reference = %loc(fcl_cfseqopen)
         call lib$signal(fcl_cfseqopen,%val(1),%val(status))
         return
      end if
c
c     Get the reference data sets at the specified start time.
c
      status = cct_get_config_tod(start,2,config_size,config_lun,config_index,
     .                            config_record,config_new_segment,
     .                            config_status)
      if (.not. status) then
         fcl_reference = %loc(fcl_cfseqget)
         call ct_binary_to_gmt(start,jstart)
         call lib$signal(fcl_cfseqget,%val(2),jstart,%val(status))
         return
      end if
c
c     Get nyquist frequency and calibration model.
c
      if ((channel .eq. 1) .or. (channel .eq. 3)) then
         nyquistl_hz = config_record.nyquistl.hz(scan_mode)
      else
         nyquistl_hz = config_record.nyquistl.hz(scan_mode+4)
      end if

      model = config_record.model
c
c     Close sequential access configuration files.
c
      status = cct_close_config(2,config_lun,config_index)
      if (.not. status) then
         fcl_reference = %loc(fcl_cfseqclose)
         call lib$signal(fcl_cfseqclose,%val(1),%val(status))
         return
      end if

      return
      end
