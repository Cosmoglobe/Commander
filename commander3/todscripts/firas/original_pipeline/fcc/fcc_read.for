      integer*4 function fcc_read(in_cov_name,in_cov,start,stop,sci_mode,
     .                            adds_per_group)
c-------------------------------------------------------------------------------
c
c     Purpose: Read input covariance matrix.
c
c     Author: S. Alexander, HSTX, 7/93, SER 11189
c
c     Input: in_cov_name  ch*72  Input covariance matrix filename.
c
c     Output: in_cov          rec(257)  Input covariance matrix records.
c             start           i*4(2)  Start time tag in catalog.
c             stop            i*4(2)  Stop time tag in catalog.
c             sci_mode        i*4  Value of science mode.
c             adds_per_group  i*4  Value of adds per group.
c
c     Modifications:
c    
c-------------------------------------------------------------------------------
      implicit none
c
c     Include files.
c
      include 'csdr$library:ctparams.inc'
      include '(cct_query_catalog_record)'
c
c     Return statuses.
c
      external fcc_normal
      external fcc_ctopen,fcc_ctread,fcc_ctclose
c
c     External references.
c
      external ct_connect_read
c
c     Functions.
c
      integer*4 fut_get_lun
      integer*4 fut_free_lun
      integer*4 cct_query_catalog
c
c     Input parameters.
c
      character*72 in_cov_name
c
c     Output parameters.
c
      dictionary 'fic_cov'
      record /fic_cov/ in_cov(257)

      integer*4 start(2),stop(2),sci_mode,adds_per_group
c
c     Local variables.
c
      integer*4 status,lun,len,rec
      integer*2 ct_status(20)
   
      record /query_catalog/ query_cat

      dictionary 'ccm_cme_catalog_entry'
      record /ccm_cme_catalog_entry/ cat

      fcc_read = %loc(fcc_normal)
c
c     Open the input covariance matrix.
c
      status = fut_get_lun(lun)
      if (.not. status) then
         fcc_read = status
         return
      end if

      open(unit=lun,file=in_cov_name,status='old',iostat=status,
     .      useropen=ct_connect_read)
      if (status .ne. 0) then
         fcc_read = %loc(fcc_ctopen)
         call str$trim(in_cov_name,in_cov_name,len)
         call lib$signal(fcc_ctopen,%val(2),in_cov_name(:len),%val(status))
         return
      end if
c
c     Read the covariance matrix records.
c         
      do rec = 1,257
         call ct_read_arcv(,lun,in_cov(rec),ct_status)
         if (ct_status(1) .ne. ctp_normal) then
            fcc_read = %loc(fcc_ctread)
            call str$trim(in_cov_name,in_cov_name,len)
            call lib$signal(fcc_ctread,%val(2),in_cov_name(:len),
     .                      %val(ct_status(1)))
            return
         end if
      end do
c
c     Close the covariance matrix.
c
      call ct_close_arcv(,lun,ct_status)
      if (ct_status(1) .ne. ctp_normal) then
         fcc_read = %loc(fcc_ctclose)
         call str$trim(in_cov_name,in_cov_name,len)
         call lib$signal(fcc_ctclose,%val(2),in_cov_name(:len),
     .                   %val(ct_status(1)))
         return
      end if

      status = fut_free_lun(lun)
      if (.not. status) then
         fcc_read = status
         return
      end if
c
c     Query the catalog for the start and stop times so that they may be
c     copied to the output covariance matrix.
c
      query_cat.archive_id = 'CSDR$FIRAS_IN'
      query_cat.filename = in_cov_name(15:)

      status = cct_query_catalog(query_cat,cat)
      if (.not. status) then
         fcc_read = status
         return
      end if

      start(1) = cat.initial_time(1)
      start(2) = cat.initial_time(2)
      stop(1) = cat.final_time(1)
      stop(2) = cat.final_time(2)
c
c     Return the science mode and adds-per-group values from the matrix.
c
      sci_mode = in_cov(1).ident.sci_mode
      adds_per_group = in_cov(1).ident.adds_per_group

      return
      end
