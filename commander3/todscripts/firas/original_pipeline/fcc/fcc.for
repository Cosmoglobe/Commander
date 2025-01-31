      program fcc
c----------------------------------------------------------------------------
c
c     Purpose: Calibrate covariance matrices.
c
c     Author: S. Alexander, HSTX, 7/93, SER 11189
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
      external fcc_normal
      external fcc_ctinit
      external fcc_failure
c
c     External references.
c
      external fut_error
c
c     Functions.
c
      integer*4 cut_register_version,cut_display_banner
      integer*4 fcc_parse,fcc_report,fcc_read
      integer*4 fcc_reference,fcc_gain,fcc_extract
      integer*4 fcc_calibrate,fcc_insert,fcc_write
c
c     Local variables for fcc_parse.
c
      integer*4 channel,scan_mode
      character*72 cov_ext,mod_ext
      real*4 bol_avg,cbias_avg,volt_avg
      logical*1 label
      character*40 label_val
      character*79 command_line(7)
c
c     Local variables for fcc_report.
c
      character*72 in_cov_name,out_cov_name
c
c     Local variables for fcc_read.
c     
      dictionary 'fic_cov'
      record /fic_cov/ in_cov(257)

      integer*4 start(2),stop(2),sci_mode,adds_per_group
c
c     Local variables for fcc_reference.
c
      real*4 nyquist_hz

      dictionary 'fex_mod'
      record /fex_mod/ model

      real*8 apod(512)
      complex*16 etf(512)
      integer*4 peak
c
c     Local variables for fcc_gain.
c
      real*8 s0,tau,tbol
      complex*16 gain(257)
c
c     Local variables for fcc_extract.
c
      dictionary 'fcc_cov'
      record /fcc_cov/ out_cov(257)

      real*8 c(522,522)
      complex*16 n(512,512),n2(512,512),eq(8,512),da(257,512)
c
c     Local variables.
c
      integer*4 status
      character*6 version
      parameter(version='11.4')
      integer*2 ct_status(20)
c
c     Register version and display banner.
c
      status = cut_register_version(version)
      status = cut_display_banner(6,80,
     .                            'FIRAS Facility FCC_Calibrate_Covariances')
c
c     Parse command line.
c
      status = fcc_parse(channel,scan_mode,cov_ext,mod_ext,bol_avg,cbias_avg,
     .                   volt_avg,label,label_val,command_line)
c
c     Initialize report.
c
      status = fcc_report(channel,scan_mode,cov_ext,label,command_line,
     .                    status,in_cov_name,out_cov_name)
c
c     Initialize Cobetrieve and error handler.
c
      if (status .eq. %loc(fcc_normal)) then
         call lib$establish(fut_error)
         call ct_init(ct_status)
         if (ct_status(1) .ne. ctp_normal) then
            status = %loc(fcc_ctinit)  
            call lib$signal(fcc_ctinit,%val(1),%val(ct_status(1)))
         end if
      end if
c
c     Read input covariance matrix.
c
      if (status .eq. %loc(fcc_normal)) then
         status = fcc_read(in_cov_name,in_cov,start,stop,sci_mode,
     .                     adds_per_group)
      end if
c
c     Read reference data sets.
c
      if (status .eq. %loc(fcc_normal)) then
         status = fcc_reference(channel,scan_mode,sci_mode,adds_per_group,
     .                          start,nyquist_hz,model,apod,etf,peak)
      end if
c
c     Calculate gain function.
c
      if (status .eq. %loc(fcc_normal)) then
         status = fcc_gain(bol_avg,cbias_avg,volt_avg,nyquist_hz,model,s0,tau,
     .                     tbol,gain)
      end if
c
c     Create four intermediate matrices.
c
      if (status .eq. %loc(fcc_normal)) then
         status = fcc_extract(label,label_val,bol_avg,cbias_avg,volt_avg,model,
     .                        s0,tau,tbol,in_cov,out_cov,c,n,n2,eq,da)
      end if
c
c     Calibrate four intermediate matrices.
c
      if (status .eq. %loc(fcc_normal)) then
         status = fcc_calibrate(apod,etf,gain,peak,n,n2,eq,da)
      end if
c
c     Calculate output covariance matrix.
c
      if (status .eq. %loc(fcc_normal)) then
         status = fcc_insert(c,n,n2,eq,da,out_cov)
      end if
c
c     Write output covariance matrix records.
c
      if (status .eq. %loc(fcc_normal)) then
         status = fcc_write(out_cov_name,out_cov,start,stop)
      end if

      if (status .eq. %loc(fcc_normal)) then
         call lib$signal(fcc_normal)
      else
         call lib$signal(fcc_failure)
      end if

      close(fut_report_lun)

      end
