      integer*4 function fil_read(input,channel,last_channel,next_channel,
     .                            scan_mode_input,rse_input,rse,neighbors,
     .                            orphans,pool_max,single,short_lun,sci_lun,
     .                            eng_lun,num_recs,num_cgr_recs,sci_recs,
     .                            eng_recs)
c------------------------------------------------------------------------------
c
c     Purpose: Read coaddable groups of science and engineering records.
c
c     Author: S. Brodd, HSTX, 4/95
c
c     Input: input              i*4  Type of input: sky, calibration, or raw.
c            channel            i*4  Value of channel, 1-4 = RH-LL.
c            last_channel       i*4  Value of channel of last coadd group, 1-4.
c            next_channel       l*1  Process the next channel.
c            scan_mode_input    i*4  Input scan mode, 1-6 = SS-FL; 0 if none.
c            rse_input          i*4  Rse specified on command line.
c            rse                ch*128(16)  Record selection expressions.
c            neighbors          i*4  Use neighbor interferograms.
c            orphans            i*4  Read orphan interferograms.
c            pool_max           i*4  Maximum number of interferograms to read
c                                    if neighbors are required.
c            single             i*4  Process single ifgs.
c            short_lun          i*4  Logical unit for short science file(s).
c            sci_lun            i*4  Logical unit for science file(s).
c            eng_lun            i*4  Logical unit for engineering file(s).
c
c     Output: last_channel       i*4  Value of channel of last coadd group, 1-4.
c             next_channel       l*1  Process the next channel.
c             num_recs           i*4  Number of input science and engineering
c                                     records, including neighbors.
c             num_cgr_recs       i*4  Number of input science and engineering
c                                     records, excluding neighbors.
c             sci_recs           rec(num_recs)  Science records.
c             eng_recs           rec(num_recs)  Engineering records.
c
c     Modifications:
c
c------------------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_params)'
      include 'csdr$library:ctparams.inc'
      include '(csa_pixel_input_rec)'
      include '(csa_pixel_output_rec)'
c
c     Return statuses.
c
      external fil_normal
      external fil_csaread
      external fil_ctread
      external fil_ctkeyread
      external fil_ctquery
c
c     Functions.
c
      integer*4 csa_read_pixels
c
c     Input parameters.
c
      integer*4 input,channel,scan_mode_input,rse_input
      integer*4 neighbors,orphans,pool_max,single
      character*128 rse(16)
      integer*4 short_lun,sci_lun,eng_lun
c
c     Input/output parameters.
c
      integer*4 last_channel
      logical*1 next_channel
c
c     Output parameters.
c
      integer*4 num_recs,num_cgr_recs

      dictionary 'fdq_sdf'
      dictionary 'fdq_eng'
      record /fdq_sdf/ sci_recs(fac_max_num)
      record /fdq_eng/ eng_recs(fac_max_num)
c
c     Local variables.
c
      record /pixel_input_list/ in_pixel
      record /pixel_output_list/ out_pixel

      logical*1 next_pixel,trans_flag
      integer*4 scan_mode,status,num_pixels,num_pixel_recs
      integer*4 check_rec,rec_scan_mode,pixel_rec,rec

      dictionary 'fss_sssky'
      record /fss_sssky/ short_recs(fac_max_short)

      character*72 filename
      integer*4 file_len
      integer*2 ct_status(20)
      character*14 jtime

      fil_read = %loc(fil_normal)
c
c     If requested scan mode is FS, read SF coadd groups; if FL, read LF
c     coadd groups.
c
      scan_mode = scan_mode_input
      if (scan_mode .eq. 5) then
         scan_mode = 2
      else if (scan_mode .eq. 6) then
         scan_mode = 4
      end if
c
c     Read short science data.
c
      if (input .ne. fac_raw) then
c
c     Read short science skymap.
c
         if (input .eq. fac_sky) then

            in_pixel.level_no = fac_skymap_level
c
c     If reading a new channel, set pixel number to -1 and flag for next pixel.
c
            if (channel .ne. last_channel) then
               last_channel = channel
               in_pixel.pixel_no = -1
               next_pixel = .true.
            end if
c
c     Loop to find a pixel with data, and if an input scan mode was specified
c     on the command line, with data in that scan mode.
c
            do while (next_pixel)
c
c     Increment pixel number.  If it is greater than 6143, stop and return flag
c     to go to next channel.
c
               in_pixel.pixel_no = in_pixel.pixel_no + 1
               if (in_pixel.pixel_no .gt. 6143) then 
                  next_channel = .true.
                  next_pixel = .false.
                  num_recs = 0
                  num_cgr_recs = 0
               else
c
c     Read data from this pixel.
c
                  status = csa_read_pixels(short_lun,in_pixel,1,short_recs,
     .                                     fac_max_short,out_pixel,
     .                                     num_pixels,0)
                  if (.not. status) then
                     fil_read = %loc(fil_csaread)
                     inquire(short_lun,name=filename)
                     call str$trim(filename,filename,file_len)
                     call lib$signal(fil_csaread,%val(2),filename(:file_len),
     .                               %val(status))
                     return
                  end if
                  num_pixel_recs = out_pixel.no_records
c
c     If there are records for this pixel and a scan mode was specified on the
c     command line, check each record.
c
                  if ((scan_mode .gt. 0) .and. 
     .                (num_pixel_recs .gt. 0)) then
                     check_rec = 0
                     do rec = 1,num_pixel_recs
                        rec_scan_mode = short_recs(rec).mtm_scan_length*2 +
     .                                  short_recs(rec).mtm_scan_speed + 1
                        if (rec_scan_mode .eq. scan_mode) then
                           check_rec = check_rec + 1
                           short_recs(check_rec) = short_recs(rec)
                        end if
                     end do
                     num_pixel_recs = check_rec
                  end if
c
c     If there are usable records left for this pixel, go on.  Otherwise, read
c     the next pixel.
c
                  if (num_pixel_recs .gt. 0) then 
                     next_pixel = .false.
                     pixel_rec = 0
                  end if
               end if
            end do
c
c     If a pixel with usable data was found in the current channel, then
c     return a single record if processing single interferograms or determine
c     the coadd group marked in the short science skymap.
c
            if (.not. next_channel) then
               num_recs = 1
               num_cgr_recs = 1
               if (single .eq. fac_present) then
                  pixel_rec = pixel_rec + 1
                  short_recs(num_recs) = short_recs(pixel_rec)
c
c     When all records in this pixel have been returned, go on to the next 
c     pixel.
c
                  if (pixel_rec .ge. num_pixel_recs) then
                     next_pixel = .true.
                  end if
               else
c
c     If orphans should not be used, find the beginning of a coadd group by 
c     finding the first transition flag equal to 2.
c
                  trans_flag = 0
                  if (orphans .eq. fac_not_present) then
                     do while ((.not. next_pixel) .and. (trans_flag .ne. 2))
                        pixel_rec = pixel_rec + 1
                        short_recs(num_recs) = short_recs(pixel_rec)
                        trans_flag = short_recs(num_recs).transition_flag
c
c     If there are no records with transition flag 2, there are no coaddable
c     groups.  Go on to the next pixel.
c
                        if (pixel_rec .ge. num_pixel_recs) then
                           num_recs = 0
                           num_cgr_recs = 0
                           next_pixel = .true.
                        end if
                     end do
                  else
c
c     If orphans should be used, find the beginning of a coadd group by 
c     finding the first transition flag equal to 1 or 2.
c
                     do while ((.not. next_pixel) .and. ((trans_flag .eq. 0)
     .                          .or. (trans_flag .eq. 3)))
                        pixel_rec = pixel_rec + 1
                        short_recs(num_recs) = short_recs(pixel_rec)
                        trans_flag = short_recs(num_recs).transition_flag
c
c     If there are no records with transition flags 1 or 2, there are no 
c     coaddable or orphan groups.  Go on to the next pixel.
c
                        if (pixel_rec .ge. num_pixel_recs) then
                           num_recs = 0
                           num_cgr_recs = 0
                           next_pixel = .true.
                        end if
                     end do
                  end if
                  trans_flag = 0
c
c     Determine the coaddable group by reading records with transition flag
c     of zero.  The next transition flag of 2 will indicate the next coaddable
c     group.  A transition flag of 1 will indicate an orphan group.  A
c     transition flag or 3 indicates the beginning of a neighbor group.  A
c     maximum of fac_max_num records can be coadded.
c
                  do while ((.not. next_pixel) .and. (trans_flag .eq. 0) .and.
     .                      (num_recs .lt. fac_max_num))
                     num_recs = num_recs + 1
                     num_cgr_recs = num_cgr_recs + 1
                     pixel_rec = pixel_rec + 1
                     short_recs(num_recs) = short_recs(pixel_rec)
                     trans_flag = short_recs(num_recs).transition_flag
c
c     If this is the last group in the pixel, return the coadd group and go on 
c     to the next pixel for the next read.
c
                     if (pixel_rec .ge. num_pixel_recs) then
                        next_pixel = .true.
                     end if
                  end do
c
c     If the last transition flag read was not zero, then the current record
c     is the beginning of another group.  Decrement pointers appropriately.
c
                  if (trans_flag .ne. 0) then
                     num_cgr_recs = num_cgr_recs - 1
                  end if
c
c     Read neighbor records if required and if they exist.
c
                  if ((.not. next_pixel) .and. ((neighbors .eq. fac_present) 
     .                .or. (orphans .eq. fac_present)) .and. 
     .                (num_cgr_recs .lt. pool_max) .and. 
     .                (trans_flag .eq. 3)) then
                     trans_flag = 0
c
c     Read neighbor records until the end of records for this pixel is reached,
c     the next coadd, orphan, or neighbor group is encountered, or until
c     pool_max records have been read.
c
                     do while ((.not. next_pixel) .and. (trans_flag .eq. 0) 
     .                         .and. (num_recs .lt. pool_max))
                        num_recs = num_recs + 1
                        pixel_rec = pixel_rec + 1
                        short_recs(num_recs) = short_recs(pixel_rec)
                        trans_flag = short_recs(num_recs).transition_flag
c
c     If this is the last group in the pixel, return the coadd group and go on 
c     to the next pixel for the next read.
c
                        if (pixel_rec .ge. num_pixel_recs) then
                           next_pixel = .true.
                        end if
                     end do
                  end if
c
c     If the last transition flag read was not zero, then the current record
c     is the beginning of another group.  Decrement pointers appropriately.
c
                  if (trans_flag .ne. 0) then
                     num_recs = num_recs - 1
                     pixel_rec = pixel_rec - 1
                  end if
c
c     For the special case in which a lone neighbor interferogram is the last
c     record for the pixel, and if it is needed, then the num_rec pointer
c     should not be decremented.  
c
                  if (next_pixel .and. (trans_flag .eq. 3) .and.
     .                ((neighbors .eq. fac_present) .or.
     .                (orphans .eq. fac_present)) .and.
     .                (num_recs .ge. 1) .and.
     .                (num_recs .lt. pool_max)) then
                     num_recs = num_recs + 1
                  end if
               end if
            end if
         else
c
c     Read short science records for calibration data.  In this case, only
c     coaddable groups will be found in the short science file(s), and the last
c     record of each coadd group will be marked with a transition flag of true.
c     Only a maximum of fac_max_num records can be coadded.
c
            num_recs = 0
            trans_flag = .false.
            do while ((.not. trans_flag) .and. (num_recs .lt. fac_max_num))
               num_recs = num_recs + 1
               call ct_read_arcv(,short_lun,short_recs(num_recs),ct_status)
c
c     When the end of the short science file(s) is reached, go on to the 
c     next channel.
c
               if (ct_status(1) .eq. ctp_endoffile) then
                  next_channel = .true.
                  trans_flag = .true.
                  num_recs = num_recs - 1
               else if (ct_status(1) .ne. ctp_normal) then
                  fil_read = %loc(fil_ctread)
                  inquire(short_lun,name=filename)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fil_ctread,%val(2),filename(:file_len),
     .                            %val(ct_status(1)))
                  return
c
c     If a scan mode was specified on the command line, check the record.
c
               else if ((scan_mode .gt. 0) .and.
     .                  ((short_recs(num_recs).mtm_scan_length*2 + 
     .                    short_recs(num_recs).mtm_scan_speed + 1) .ne.
     .                  scan_mode)) then
                  num_recs = num_recs - 1
               else
                  trans_flag = short_recs(num_recs).transition_flag
               end if
c
c     Return the short science record when processing single interferograms.
c
               if (single .eq. fac_present) then
                  trans_flag = .true.
               end if
            end do
c
c     The number of records is always the same as the number of records in
c     the coadd group for calibration data, as no neighbors are used.
c
            num_cgr_recs = num_recs
         end if
c
c     Retrieve associated full science record.
c
         if (num_recs .gt. 0) then
            do rec=1,num_recs

               call ct_keyedread_arcv(,sci_lun,sci_recs(rec),
     .                                short_recs(rec).time(1),ct_status)
               if (ct_status(1) .ne. ctp_normal) then
                  fil_read = %loc(fil_ctkeyread)
                  call ct_binary_to_gmt(short_recs(rec).time,jtime)
                  inquire(sci_lun,name=filename)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fil_ctkeyread,%val(3),filename(:file_len),
     .                            jtime,%val(ct_status(1)))
                  return
               end if
c
c     Write short science pixelization information into the full science 
c     records if processing sky data.
c
               if (input .eq. fac_sky) then
c
c     Copy pixel number only for non-neighbor records.
c
                  if (rec .le. num_cgr_recs) then
                     sci_recs(rec).attitude.pixel_no = short_recs(rec).pixel_no
                  end if
                  sci_recs(rec).attitude.pixel_definition = 
     .                          short_recs(rec).pixel_definition
                  sci_recs(rec).attitude.skymap_index = 
     .                          short_recs(rec).skymap_index
                  sci_recs(rec).attitude.exc_galactic_lat =
     .                          short_recs(rec).exc_galactic_lat
               end if
c
c     Retrieve the engineering record.
c
               call ct_keyedread_arcv(,eng_lun,eng_recs(rec),
     .                                sci_recs(rec).dq_data.eng_time(1),ct_status)
               if (ct_status(1) .ne. ctp_normal) then
                  fil_read = %loc(fil_ctkeyread)
                  call ct_binary_to_gmt(sci_recs(rec).ct_head.time,jtime)
                  inquire(eng_lun,name=filename)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fil_ctkeyread,%val(3),filename(:file_len),
     .                            jtime,%val(ct_status(1)))
                  return
               end if
            end do
         end if
      else 
c
c     When processing raw data and an rse file has been specified, query
c     the archive on the first read for the current channel.
c
         if ((rse_input .eq. fac_present) .and. 
     .       (channel .ne. last_channel)) then
            last_channel = channel
            call ct_query_arcv(,eng_lun,rse,ct_status)
            if (ct_status(1) .ne. ctp_normal) then
               fil_read = %loc(fil_ctquery)
               inquire(eng_lun,name=filename)
               call str$trim(filename,filename,file_len)
               call lib$signal(fil_ctquery,%val(2),filename(:file_len),
     .                         %val(ct_status(1)))
               return
            end if
         end if

         num_recs = 0
         trans_flag = .false.
c
c     Read engineering records until the end of the file, or until fac_max_num
c     records have been read.
c
         do while ((num_recs .lt. fac_max_num) .and. (.not. trans_flag))

            num_recs = num_recs + 1
c
c     Read by query if an rse is being used.
c
            if (rse_input .eq. fac_present) then
               call ct_query_get(,eng_lun,eng_recs(num_recs),ct_status)
            else
               call ct_read_arcv(,eng_lun,eng_recs(num_recs),ct_status)
            end if

            if (ct_status(1) .eq. ctp_endoffile) then
               next_channel = .true.
               trans_flag = .true.
               num_recs = num_recs - 1
            else if (ct_status(1) .ne. ctp_normal) then
               fil_read = %loc(fil_ctread)
               inquire(eng_lun,name=filename)
               call str$trim(filename,filename,file_len)
               call lib$signal(fil_ctread,%val(2),filename(:file_len),
     .                         %val(ct_status(1)))
               return
c
c     If the associated science records time field is zero, skip this 
c     engineering record.
c
            else if ((eng_recs(num_recs).en_head.sci_time(channel).bin_time(1) 
     .                .eq. 0) .and.
     .               (eng_recs(num_recs).en_head.sci_time(channel).bin_time(2) 
     .                .eq. 0)) then
                  num_recs = num_recs - 1
            else
               call ct_keyedread_arcv(,sci_lun,sci_recs(num_recs),
     .              eng_recs(num_recs).en_head.sci_time(channel).bin_time(1),
     .              ct_status)
               if (ct_status(1) .ne. ctp_normal) then
                  fil_read = %loc(fil_ctkeyread)
                  call ct_binary_to_gmt(eng_recs(num_recs).ct_head.time,jtime)
                  inquire(sci_lun,name=filename)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fil_ctkeyread,%val(3),filename(:file_len),
     .                            jtime,%val(ct_status(1)))
                  return
               end if
               if (single .eq. fac_present) then
                  trans_flag = .true.
               end if
            end if
         end do
c
c     The number of records is always the same as the number of records in
c     the coadd group for raw data, as no neighbors are used.
c
         num_cgr_recs = num_recs
      end if

      return
      end
