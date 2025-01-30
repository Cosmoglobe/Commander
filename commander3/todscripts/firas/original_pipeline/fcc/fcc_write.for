      integer*4 function fcc_write(out_cov_name,out_cov,start,stop)
c-------------------------------------------------------------------------------
c
c     Purpose: Write output covariance matrix.
c
c     Author: S. Alexander, HSTX, 7/93, SER 11189
c
c     Input: out_cov_name  ch*72  Output covariance matrix filename.
c            out_cov       rec(257)  Output covariance matrix records.
c            start         i*4(2)  Binary start time of covariance matrix.
c            stop          i*4(2)  Binary stop time of covariance matrix.
c
c     Output: 
c
c     Modifications:
c    
c-------------------------------------------------------------------------------
      implicit none
c
c     Include files.
c
      include 'csdr$library:ctparams.inc'
c
c     Return statuses.
c
      external fcc_normal
      external fcc_ctopen
      external fcc_ctwrite
      external fcc_ctclose
c
c     External references.
c
      external ct_connect_write
c
c     Functions.
c
      integer*4 fut_get_lun
      integer*4 cct_set_ttg_time_range
      integer*4 fut_free_lun
c
c     Input parameters.
c
      character*72 out_cov_name

      dictionary 'fcc_cov'
      record /fcc_cov/ out_cov(257)

      integer*4 start(2),stop(2)
c
c     Local variables.
c
      integer*4 status,lun,len,rec
      integer*2 ct_status(20)

      fcc_write = %loc(fcc_normal)
c
c     Open output covariance matrix.
c
      status = fut_get_lun(lun)
      if (.not. status) then
         fcc_write = status
         return
      end if

      open(unit=lun,file=out_cov_name,status='new',iostat=status,
     .     useropen=ct_connect_write)
      if (status .ne. 0) then
         fcc_write = %loc(fcc_ctopen)
         call str$trim(out_cov_name,out_cov_name,len)
         call lib$signal(fcc_ctopen,%val(2),out_cov_name(:len),%val(status))
         return
      end if
c
c     Write covariance matrix records.
c
      do rec = 1,257
         call ct_write_arcv(,lun,out_cov(rec),ct_status)
         if (ct_status(1) .ne. ctp_normal) then
            fcc_write = %loc(fcc_ctwrite)
            call str$trim(out_cov_name,out_cov_name,len)
            call lib$signal(fcc_ctwrite,%val(2),out_cov_name(:len),
     .                      %val(ct_status(1)))
            return
         end if
      end do
c
c     Set the time tags for the covariance matrix in the catalog.
c
      status = cct_set_ttg_time_range(lun,start,stop)
c
c     Close covariance matrix.
c
      call ct_close_arcv(,lun,ct_status)
      if (ct_status(1) .ne. ctp_normal) then
         fcc_write = %loc(fcc_ctclose)
         call str$trim(out_cov_name,out_cov_name,len)
         call lib$signal(fcc_ctclose,%val(2),out_cov_name(:len),
     .                   %val(ct_status(1)))
         return
      end if

      status = fut_free_lun(lun)
      if (.not. status) then
         fcc_write = status
         return
      end if

      return
      end
