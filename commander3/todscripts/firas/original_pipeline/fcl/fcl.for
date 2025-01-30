      program fcl
c----------------------------------------------------------------------------
c
c     Purpose: Calibrate covariance matrices.
c
c     Author: S. Brodd, HSTX, 12/95, SPR 12291
c
c     Modifications:
c
c----------------------------------------------------------------------------
      implicit none
c
c     Include files.
c
      include '(fut_error)'
      include 'csdr$library:ctparams.inc'
c
c     Return statuses.
c
      external fcl_normal
      external fcl_ctinit
      external fcl_failure
c
c     External references.
c
      external fut_error
c
c     Functions.
c
      integer*4 cut_register_version,cut_display_banner
      integer*4 fcl_parse,fcl_report,fcl_read
      integer*4 fcl_reference,fcl_gain,fcl_extract
      integer*4 fcl_calibrate,fcl_insert,fcl_write
c
c     Local variables for fcl_parse.
c
      integer*4 channel,scan_mode
      character*72 cov_ext,mod_ext
      real*4 bol_avg,cbias_avg,volt_avg
      logical*1 label
      character*40 label_val
      character*79 command_line(7)
c
c     Local variables for fcl_report.
c
      character*72 in_cov_name,out_cov_name
c
c     Local variables for fcl_read.
c     
      dictionary 'fil_cov'
      record /fil_cov/ in_cov(257)

      integer*4 start(2),stop(2)
c
c     Local variables for fcl_reference.
c
      real*4 nyquistl_hz

      dictionary 'fex_mod'
      record /fex_mod/ model
c
c     Local variables for fcl_gain.
c
      real*8 s0,tau,tbol
      complex*16 gain(361)
c
c     Local variables for fcl_extract.
c
      dictionary 'fcl_cov'
      record /fcl_cov/ out_cov(257)

      real*8 c(730,730)
      complex*16 n(361,361),n2(361,361),eq(8,361),da(257,361)
c
c     Local variables.
c
      integer*4 status
      character*6 version
      parameter(version='13.9')
      integer*2 ct_status(20)
c
c     Register version and display banner.
c
      status = cut_register_version(version)
      status = cut_display_banner(6,80,
     .             'FIRAS Facility FCL_CalibrateCovariances_Long')
c
c     Parse command line.
c
      status = fcl_parse(channel,scan_mode,cov_ext,mod_ext,bol_avg,cbias_avg,
     .                   volt_avg,label,label_val,command_line)
c
c     Initialize report.
c
      status = fcl_report(channel,scan_mode,cov_ext,label,command_line,
     .                    status,in_cov_name,out_cov_name)
c
c     Initialize Cobetrieve and error handler.
c
      if (status .eq. %loc(fcl_normal)) then
         call lib$establish(fut_error)
         call ct_init(ct_status)
         if (ct_status(1) .ne. ctp_normal) then
            status = %loc(fcl_ctinit)  
            call lib$signal(fcl_ctinit,%val(1),%val(ct_status(1)))
         end if
      end if
c
c     Read input covariance matrix.
c
      if (status .eq. %loc(fcl_normal)) then
         status = fcl_read(in_cov_name,in_cov,start,stop)
      end if
c
c     Read reference data sets.
c
      if (status .eq. %loc(fcl_normal)) then
         status = fcl_reference(channel,scan_mode,start,nyquistl_hz,model)
      end if
c
c     Calculate gain function.
c
      if (status .eq. %loc(fcl_normal)) then
         status = fcl_gain(scan_mode,bol_avg,cbias_avg,volt_avg,nyquistl_hz,
     .                     model,s0,tau,tbol,gain)
      end if
c
c     Create four intermediate matrices.
c
      if (status .eq. %loc(fcl_normal)) then
         status = fcl_extract(label,label_val,bol_avg,cbias_avg,volt_avg,model,
     .                        s0,tau,tbol,in_cov,out_cov,c,n,n2,eq,da)
      end if
c
c     Calibrate four intermediate matrices.
c
      if (status .eq. %loc(fcl_normal)) then
         status = fcl_calibrate(gain,n,n2,eq,da)
      end if
c
c     Calculate output covariance matrix.
c
      if (status .eq. %loc(fcl_normal)) then
         status = fcl_insert(c,n,n2,eq,da,out_cov)
      end if
c
c     Write output covariance matrix records.
c
      if (status .eq. %loc(fcl_normal)) then
         status = fcl_write(out_cov_name,out_cov,start,stop)
      end if

      if (status .eq. %loc(fcl_normal)) then
         call lib$signal(fcl_normal)
      else
         call lib$signal(fcl_failure)
      end if

      close(fut_report_lun)

      end
