      integer*4 function fil_reference(jstart,rse_input,rse_name,init,rse,
     .                                 mincoadd,grttrans,grtrawwt,grtcoawt,
     .                                 cmdgain,samprate,nyquistl,cth,gltchcor,
     .                                 basis,reftemps_lun,dtrf_lun,gltchpro_lun,
     .                                 apodl_lun,etfl_lun)
c-------------------------------------------------------------------------------
c
c     Purpose: Retrieve all reference data.
c
c     Author: S. Brodd, HSTX, 4/95
c
c     Input: jstart             ch*14  Value of start time.
c            rse_input          i*4  Rse specified on command line.
c            rse_name           ch*72  Rse filename.
c            init               l*1  Open or close reference data sets.
c
c     Output: rse                ch*128(16)  Record selection expressions.
c             mincoadd           rec  Minimum interferograms to coadd.
c             grttrans           rec  Transition temperatures.
c             grtrawwt           rec  Raw temperature weights.
c             grtcoawt           rec  Coadd temperature weights.
c             cmdgain            rec  Real values of commanded gains.
c             samprate           rec  Sampling rate.
c             nyquistl           rec  Nyquist frequencies.
c             cth                rec  Consistency check parameters.
c             gltchcor           rec  Glitch rate weighting parameters.
c             basis              rec  Legendre polynomial basis vectors.
c             reftemps_lun       i*4  Logical unit for reference templates.
c             dtrf_lun           i*4  Logical unit for digital transient
c                                     response functions.
c             gltchpro_lun       i*4  Logical unit for glitch profiles.
c             apodl_lun          i*4  Logical unit for apodization functions.
c             etfl_lun           i*4  Logical unit for electronics transfer
c                                     functions.
c      
c     Modifications:
c    
c-------------------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_params)'
      include '(fut_error)'
      include '(cct_get_config)'
c
c     Return statuses.
c 
      external fil_normal
      external fil_rseopen,fil_rseread,fil_rseclose
      external fil_repwrite
      external fil_cfseqopen,fil_cfseqget,fil_cfseqclose
      external fil_cfdiropen,fil_cfdirget,fil_cfdirclose
c
c     Functions.
c
      integer*4 fut_get_lun
      integer*4 fut_free_lun
      integer*4 cct_open_config
      integer*4 cct_get_config_tod
      integer*4 cct_get_config_idx_tod
      integer*4 cct_close_config      
c
c     Input parameters.
c
      character*14 jstart
      integer*4 rse_input
      character*72 rse_name
      logical*1 init
c
c     Output parameters.
c
      character*128 rse(16)

      dictionary 'fex_mincoadd'
      dictionary 'fex_grttrans'
      dictionary 'fex_grtrawwt'
      dictionary 'fex_grtcoawt'
      dictionary 'fex_cmdgain'
      dictionary 'fex_samprate'
      dictionary 'fex_nyquistl'
      dictionary 'fex_cth'
      dictionary 'fex_gltchcor'
      dictionary 'fex_basis'

      record /fex_mincoadd/ mincoadd
      record /fex_grttrans/ grttrans
      record /fex_grtrawwt/ grtrawwt
      record /fex_grtcoawt/ grtcoawt
      record /fex_cmdgain/ cmdgain
      record /fex_samprate/ samprate
      record /fex_nyquistl/ nyquistl
      record /fex_cth/ cth
      record /fex_gltchcor/ gltchcor
      record /fex_basis/ basis

      integer*4 reftemps_lun,dtrf_lun,gltchpro_lun,apodl_lun,etfl_lun
c
c     Local variables.
c
      integer*4 rse_pos,status,rse_lun
      character*72 summary_name
      integer*4 file_len

      character*14 config_gmt_start,config_gmt_stop
      integer*4 config_start(2),config_stop(2)

      character*28 config_name(10)
      integer*4 config_size(10),config_lun(10),config_index(10)

      record /config_status/ config_status(10)

      integer*4 config_ref_count,start(2)

      structure /config_record/
         record /fex_mincoadd/ mincoadd
         record /fex_grttrans/ grttrans
         record /fex_grtrawwt/ grtrawwt
         record /fex_grtcoawt/ grtcoawt
         record /fex_cmdgain/ cmdgain
         record /fex_samprate/ samprate
         record /fex_nyquistl/ nyquistl
         record /fex_cth/ cth
         record /fex_gltchcor/ gltchcor
         record /fex_basis/ basis
      endstructure        

      record /config_record/ config_record

      logical*1 config_new_segment(10)

      fil_reference = %loc(fil_normal)
c
c     If this is the beginning of processing, retrieve all reference data.
c
      if (init) then
c
c     If an rse has been specified, then retrieve the rse.
c
         if (rse_input .eq. fac_present) then

            do rse_pos = 1,16
               rse(rse_pos) = ' '
            end do
c
c     Open the rse file.
c
            status = fut_get_lun(rse_lun)
            if (.not. status) then
               fil_reference = status
               return
            end if 

            open(unit=rse_lun,file=rse_name,status='old',iostat=status,shared,
     .           readonly) 
            if (status .ne. 0) then
               fil_reference = %loc(fil_rseopen)
               call str$trim(rse_name,rse_name,file_len)
               call lib$signal(fil_rseopen,%val(2),rse_name(:file_len),
     .                         %val(status))
               return
            end if
c
c     Read the record selection expressions.
c
            read(rse_lun,10,iostat=status) (rse(rse_pos),rse_pos=1,16)
10          format (a128)
            if (status .ne. 0) then
               fil_reference = %loc(fil_rseread)
               call str$trim(rse_name,rse_name,file_len)
               call lib$signal(fil_rseread,%val(2),rse_name(:file_len),
     .                         %val(status))
               return
            end if
c
c     Close the rse file.
c
            close(rse_lun,iostat=status)
            if (status .ne. 0) then
               fil_reference = %loc(fil_rseclose)
               call str$trim(rse_name,rse_name,file_len)
               call lib$signal(fil_rseclose,%val(2),rse_name(:file_len),
     .                         %val(status))
               return
            end if

            status = fut_free_lun(rse_lun)
            if (.not. status) then
               fil_reference = status
               return
            end if
c
c     Write record selection expressions to the summary report.
c
            write(fut_report_lun,20,iostat=status) 
     .            'Record Selection Expressions:'
20          format(/x,a)
            if (status .ne. 0) then
               fil_reference = %loc(fil_repwrite)
               inquire(fut_report_lun,name=summary_name)
               call str$trim(summary_name,summary_name,file_len)
               call lib$signal(fil_repwrite,%val(2),summary_name(:file_len),
     .                         %val(status))
               return
            end if

            do rse_pos = 1,16
               if (rse(rse_pos)(1:1) .ne. ' ') then
                  write(fut_report_lun,30,iostat=status) rse(rse_pos)(:79) 
30                format(x,a)
                  if (status .ne. 0) then
                     fil_reference = %loc(fil_repwrite)
                     inquire(fut_report_lun,name=summary_name)
                     call str$trim(summary_name,summary_name,file_len)
                     call lib$signal(fil_repwrite,%val(2),
     .                               summary_name(:file_len),%val(status))
                     return
                  end if
               end if
            end do
         end if
c     
c     Assign the sequential read configuration files (those with only one
c     record).  The start time must be advanced by one second and the stop
c     time decreased, or the open will fail.
c
         config_gmt_start = '86001000001000'
         call ct_gmt_to_binary(config_gmt_start,config_start)
         config_gmt_stop  = '99365235958990'
         call ct_gmt_to_binary(config_gmt_stop,config_stop)

         config_name(1) = 'CSDR$FIRAS_REF:FEX_MINCOADD'
         config_name(2) = 'CSDR$FIRAS_REF:FEX_GRTTRANS'
         config_name(3) = 'CSDR$FIRAS_REF:FEX_GRTRAWWT'
         config_name(4) = 'CSDR$FIRAS_REF:FEX_GRTCOAWT'
         config_name(5) = 'CSDR$FIRAS_REF:FEX_CMDGAIN'
         config_name(6) = 'CSDR$FIRAS_REF:FEX_SAMPRATE'
         config_name(7) = 'CSDR$FIRAS_REF:FEX_NYQUISTL'
         config_name(8) = 'CSDR$FIRAS_REF:FEX_CTH'
         config_name(9) = 'CSDR$FIRAS_REF:FEX_GLTCHCOR'
         config_name(10) = 'CSDR$FIRAS_REF:FEX_BASIS'

         config_size(1) = 128
         config_size(2) = 384
         config_size(3) = 256
         config_size(4) = 256
         config_size(5) = 256
         config_size(6) = 64
         config_size(7) = 128
         config_size(8) = 1024
         config_size(9) = 256
         config_size(10) = 20480
c
c     Open the sequential access reference data sets.
c
         status = cct_open_config(config_start,config_stop,10,config_name,
     .                            config_size,' ',1,config_lun,config_index,
     .                            config_status,config_ref_count)  
         if (.not. status) then
            fil_reference = %loc(fil_cfseqopen)
            call lib$signal(fil_cfseqopen,%val(1),%val(status))
            return
         end if
c
c     Get the reference data sets at the specified start time.
c
         call ct_gmt_to_binary(jstart,start)

         status = cct_get_config_tod(start,10,config_size,config_lun,
     .                               config_index,config_record,
     .                               config_new_segment,config_status)
         if (.not. status) then
            fil_reference = %loc(fil_cfseqget)
            call lib$signal(fil_cfseqget,%val(2),jstart,%val(status))
            return
         end if
c
c     Assign the reference data sets to the appropriate records.
c
         mincoadd = config_record.mincoadd
         grttrans = config_record.grttrans
         grtrawwt = config_record.grtrawwt
         grtcoawt = config_record.grtcoawt
         cmdgain = config_record.cmdgain
         samprate = config_record.samprate
         nyquistl = config_record.nyquistl
         cth = config_record.cth
         gltchcor = config_record.gltchcor
         basis = config_record.basis
c
c     Close sequential access configuration files.
c
         status = cct_close_config(10,config_lun,config_index)
         if (.not. status) then
            fil_reference = %loc(fil_cfseqclose)
            call lib$signal(fil_cfseqclose,%val(1),%val(status))
            return
         end if
c
c     Assign direct read configuration files (those with multiple records).
c 
         config_name(1) = 'CSDR$FIRAS_REF:FEX_REFTEMPS'
         config_name(2) = 'CSDR$FIRAS_REF:FEX_DTRF'
         config_name(3) = 'CSDR$FIRAS_REF:FEX_GLTCHPRO'
         config_name(4) = 'CSDR$FIRAS_REF:FEX_APODL'
         config_name(5) = 'CSDR$FIRAS_REF:FEX_ETFL'
                     
         config_size(1) = 256
         config_size(2) = 512
         config_size(3) = 2048
         config_size(4) = 4096
         config_size(5) = 6144
c
c     Open the direct access reference data sets.
c
         status = cct_open_config(config_start,config_stop,5,config_name,
     .                            config_size,'direct',1,config_lun,
     .                            config_index,config_status,config_ref_count)
         if (.not. status) then
            fil_reference = %loc(fil_cfdiropen)
            call lib$signal(fil_cfdiropen,%val(1),%val(status))
            return
         end if
c
c     Get the direct access reference data sets at the specified start time.
c
         status = cct_get_config_idx_tod(start,5,config_lun,config_index,
     .                                   config_new_segment,config_status)
         if (.not. status) then
            fil_reference = %loc(fil_cfdirget)
            call lib$signal(fil_cfdirget,%val(2),jstart,%val(status))
            return
         end if
c
c     Assign logical units for the direct access reference data sets.
c
         reftemps_lun = config_lun(1)
         dtrf_lun = config_lun(2)
         gltchpro_lun = config_lun(3)
         apodl_lun = config_lun(4)
         etfl_lun = config_lun(5)

      else
c
c     If this is the end of processing, close the reference data sets that have
c     logical units.
c
         status = cct_close_config(5,config_lun,config_index)
         if (.not. status) then
            fil_reference = %loc(fil_cfdirclose)
            call lib$signal(fil_cfdirclose,%val(1),%val(status))
            return
         end if
      end if

      return
      end
