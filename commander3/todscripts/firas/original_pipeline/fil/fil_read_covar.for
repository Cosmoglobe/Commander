      integer*4 function fil_read_covar(channel,scan_mode_input,covar_new,
     .                                  current_gmt,cov_recs,first_covar_times,
     .                                  last_covar_times)
c-------------------------------------------------------------------------------
c
c     Purpose: Initialize covariance matrix records to zero, or read them if
c              they are to be updated and they exist.
c
c     Author: S. Brodd, HSTX, 4/95
c
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode_input    i*4  Input scan mode, 1-6 = SS-FL; 0 if none.
c            covar_new          i*4  Covariance matrix to be initialized.
c            current_gmt        ch*14  Current time in Julian format.
c
c     Output: cov_recs           rec(1285)  Covariance matrix records.
c             first_covar_times  i*4(2,5)  Binary times of first records
c                                          contributing to covariance matrices.
c             last_covar_times   i*4(2,5)  Binary times of last records
c                                          contributing to covariance matrices.
c
c     Modifications:
c    
c-------------------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_params)'
      include '(cct_query_ttg_catalog_record)'
      include 'csdr$library:ctparams.inc'
c
c     Return statuses.
c
      external fil_normal
      external fil_ctopen,fil_ctread,fil_ctclose
c
c     External references.
c
      external ct_connect_read
c
c     Functions.
c
      integer*4 cct_query_ttg_catalog
      logical*1 time_gt
      integer*4 fut_get_lun
      integer*4 fut_free_lun
c
c     Input parameters.
c
      integer*4 channel,scan_mode_input,covar_new
      character*14 current_gmt
c
c     Output parameters.
c
      dictionary 'fil_cov'
      record /fil_cov/ cov_recs(1285)
      integer*4 first_covar_times(2,5),last_covar_times(2,5)
c
c     Local variables.
c
      integer*4 rec,time(2)
      record /query_ttg_catalog/ query_ttg_cat
      integer*4 archive_len,dataset_len,file_ext_len

      integer*4 first_mode,last_mode,mode
      integer*4 offset,scan_mode,status

      dictionary 'ccm_cme_catalog_entry'
      record /ccm_cme_catalog_entry/ cats(50)

      integer*4 num_cats,activity(2),latest,lun
      character*72 filename
      integer*4 file_len
      integer*2 ct_status(20)

      fil_read_covar = %loc(fil_normal)
c
c     Initialize 1285 covariance matrix records, 257 for each of the scan modes
c     SS, SF, LF, FS, FL.
c
      do rec = 1,257
         call lib$movc5(0,,0,fac_covar_sizel,cov_recs(rec))
         cov_recs(rec).ident.chan_id = channel
         cov_recs(rec).ident.mtm_speed = 0
         cov_recs(rec).ident.mtm_length = 0
      end do
      do rec = 258,514
         call lib$movc5(0,,0,fac_covar_sizel,cov_recs(rec))
         cov_recs(rec).ident.chan_id = channel
         cov_recs(rec).ident.mtm_speed = 1
         cov_recs(rec).ident.mtm_length = 0
      end do
      do rec = 515,771
         call lib$movc5(0,,0,fac_covar_sizel,cov_recs(rec))
         cov_recs(rec).ident.chan_id = channel
         cov_recs(rec).ident.mtm_speed = 1
         cov_recs(rec).ident.mtm_length = 1
      end do
      do rec = 772,1285
         call lib$movc5(0,,0,fac_covar_sizel,cov_recs(rec))
         cov_recs(rec).ident.chan_id = channel
         cov_recs(rec).ident.mtm_speed = 1
         cov_recs(rec).ident.mtm_length = 0
      end do
c
c     Initialize the earliest and latest times of records in the covariance
c     matrices; these will be used as time tags in the archive.
c
      call ct_gmt_to_binary(fac_jstop_default,time)
      do rec = 1,5
         first_covar_times(1,rec) = time(1)
         first_covar_times(2,rec) = time(2)
      end do
      call ct_gmt_to_binary(fac_jstart_default,time)
      do rec = 1,5
         last_covar_times(1,rec) = time(1)
         last_covar_times(2,rec) = time(2)
      end do
c
c     If requested on command line, query catalog for covariance matrices from
c     a previous processing run.
c
      if (covar_new .ne. fac_present) then

         query_ttg_cat.archive_id = 'CSDR$FIRAS_STAT'
         archive_len = 15
         call ct_gmt_to_binary(fac_jstart_default,query_ttg_cat.start_time)
         call ct_gmt_to_binary(fac_jstop_default,query_ttg_cat.stop_time)

         if (scan_mode_input .eq. 0) then
            first_mode = 1
            last_mode = 5
         else if (scan_mode_input .ge. 4) then
            first_mode = scan_mode_input - 1
            last_mode = scan_mode_input - 1
         else 
            first_mode = scan_mode_input
            last_mode = scan_mode_input
         end if

         do mode = first_mode,last_mode

            offset = (mode-1)*257
            scan_mode = mode
            if (scan_mode .ge. 3) then
               scan_mode = scan_mode + 1
            end if

            query_ttg_cat.dataset_name = 'FIL_COV_'//fac_channel_ids(channel)//
     .                                    fac_scan_mode_idsl(scan_mode)
            dataset_len = 12

            status = cct_query_ttg_catalog(query_ttg_cat,cats,num_cats)

            if (num_cats .gt. 0) then
               activity(1) = 0
               activity(2) = 0
c
c     Choose from among multiple covariance matrices by the latest activity
c     time.
c
               do rec = 1,num_cats
                  if (time_gt(cats(rec).activity_time,activity)) then
                     activity(1) = cats(rec).activity_time(1)
                     activity(2) = cats(rec).activity_time(2)
                     latest = rec
                  end if
               end do

               first_covar_times(1,mode) = cats(latest).initial_time(1)
               first_covar_times(2,mode) = cats(latest).initial_time(2)
               last_covar_times(1,mode) = cats(latest).final_time(1)
               last_covar_times(2,mode) = cats(latest).final_time(2)
c
c     Determine the covariance matrix filename and open it.
c
               call str$trim(cats(latest).filename_extension,
     .                       cats(latest).filename_extension,file_ext_len)
               filename = query_ttg_cat.archive_id(:archive_len)//':'//
     .                    query_ttg_cat.dataset_name(:dataset_len)//'.'//
     .                    cats(latest).filename_extension(:file_ext_len)

               status = fut_get_lun(lun)
               if (.not. status) then
                  fil_read_covar = status
                  return
               end if

               open (unit=lun,file=filename,status='old',iostat=status,
     .               useropen=ct_connect_read)
               if (status .ne. 0) then
                  fil_read_covar = %loc(fil_ctopen)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fil_ctopen,%val(2),filename(:file_len),
     .                            %val(status))
                  return
               end if
c
c     Read the covariance matrix records into the appropriate subset of the
c     1285 records.
c         
               do rec = 1,257
                  call ct_read_arcv(,lun,cov_recs(rec+offset),ct_status)
                  if (ct_status(1) .ne. ctp_normal) then
                     fil_read_covar = %loc(fil_ctread)
                     call str$trim(filename,filename,file_len)
                     call lib$signal(fil_ctread,%val(2),filename(:file_len),
     .                               %val(ct_status(1)))
                     return
                  end if
               end do
c
c     Close the covariance matrix.
c
               call ct_close_arcv(,lun,ct_status)
               if (ct_status(1) .ne. ctp_normal) then
                  fil_read_covar = %loc(fil_ctclose)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fil_ctclose,%val(2),filename(:file_len),
     .                            %val(ct_status(1)))
                  return
               end if

               status = fut_free_lun(lun)
               if (.not. status) then
                  fil_read_covar = status
                  return
               end if

            end if
         end do
      end if
c
c     Put the current time into the ct_head of the covariance matrix records.
c
      call ct_gmt_to_binary(current_gmt,time)

      do rec = 1,1285
         cov_recs(rec).ct_head.gmt = current_gmt
         cov_recs(rec).ct_head.time(1) = time(1)
         cov_recs(rec).ct_head.time(2) = time(2)
      end do

      return
      end
