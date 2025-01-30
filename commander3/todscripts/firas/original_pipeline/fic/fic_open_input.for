      integer*4 function fic_open_input(input,channel,file,file_ext,jstart,
     .                                  jstop,rse_input,outfile,outfile_ext,
     .                                  short_lun,sci_lun,eng_lun)
c-------------------------------------------------------------------------------
c
c     Purpose: Open input files for a channel.
c
c     Author: S. Alexander, STX, 9/91, SER 7985
c
c     Input: input              i*4  Type of input: sky, calibration, or raw.
c            channel            i*4  Value of channel, 1-4 = RH-LL.
c            file               i*4  Input filename specified on command line.
c            file_ext           ch*72  Value of input filename extension on
c                                      command line.
c            jstart             ch*14  Value of start time.
c            jstop              ch*14  Value of stop time.
c            rse_input          i*4  Rse specified on command line.
c            outfile            i*4  Output filename specified on command line.
c
c     Output: jstart             ch*14  Value of start time.
c             jstop              ch*14  Value of stop time.
c             outfile_ext        ch*72  Value of output filename extension on
c                                       command line.
c             short_lun          i*4  Logical unit for short science file(s).
c             sci_lun            i*4  Logical unit for science file(s).
c             eng_lun            i*4  Logical unit for engineering file(s).
c
c     Modifications:
c    
c-------------------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_params)'
      include 'csdr$library:ctparams.inc'
      include '(cct_query_catalog_record)'
      include '(cct_query_ttg_catalog_record)'
c
c     Return statuses.
c
      external fic_normal
      external fic_nocatrecs
      external fic_csaopen
      external fic_ctopen
c
c     External references.
c
      external csa_open_skymap
      external ct_connect_read
      external ct_connect_query
c
c     Functions.
c
      integer*4 cct_query_catalog
      integer*4 cct_query_ttg_catalog
      logical*1 time_gt
      integer*4 fut_get_lun
c
c     Input parameters.
c
      integer*4 input,channel,file,rse_input,outfile
      character*72 file_ext
c
c     Input/output parameters.
c
      character*14 jstart,jstop
c
c     Output parameters.
c
      character*72 outfile_ext
      integer*4 short_lun,sci_lun,eng_lun
c
c     Local variables.
c
      character*14 archive
      character*12 dataset
      integer*4 archive_len,dataset_len,file_ext_len
      character*72 filename
      integer*4 file_len

      dictionary 'ccm_cme_catalog_entry'
      record /query_catalog/ query_cat
      record /ccm_cme_catalog_entry/ cats(50)
      record /query_ttg_catalog/ query_ttg_cat

      integer*4 status, num_cats, activity(2), rec
      character*14 sky_jstart, sky_jstop
      integer*2 ct_status(20)

      fic_open_input = %loc(fic_normal)
c
c     Assign archive and input data set names.
c
      archive = 'CSDR$FIRAS_IN'
      archive_len = 13
c
c     Short science file(s) for sky or calibration input.
c
      if (input .eq. fac_sky) then
         dataset = 'FSS_SSSKY_'//fac_channel_ids(channel)
         dataset_len = 12
      else if (input .eq. fac_cal) then
         dataset = 'FEC_SSCAL_'//fac_channel_ids(channel)
         dataset_len = 12
      else 
c
c     Engineering file(s) for raw input.
c
         archive = 'CSDR$FIRAS_RAW'
         archive_len = 14
         dataset = 'FDQ_ENG'
         dataset_len = 7
      end if

      call str$trim(file_ext,file_ext,file_ext_len)
c
c     If input file is specified on command line, query the catalog using the
c     input filename extension.
c
      if (file .eq. fac_present) then

         filename = archive(:archive_len)//':'//dataset(:dataset_len)//
     .              '.'//file_ext(:file_ext_len)

         query_cat.archive_id = archive(:archive_len)
         query_cat.filename = dataset(:dataset_len)//
     .                        '.'//file_ext(:file_ext_len)

         status = cct_query_catalog(query_cat,cats(1))
c
c     Check the status of the catalog query.
c
         if (.not. status) then
            fic_open_input = %loc(fic_nocatrecs)
            call lib$signal(fic_nocatrecs,%val(2),archive(:archive_len),
     .                       dataset(:dataset_len))
            return
         end if
c
c     Determine jstart and jstop from the catalog entry.
c
         call ct_binary_to_gmt(cats(1).initial_time,jstart)
         call ct_binary_to_gmt(cats(1).final_time,jstop)
c
c     With sky input and jstart and jstop specified on the command line, query
c     the catalog with these time tags to find the short science skymap whose
c     time range overlaps the specified time range and whose activity time is 
c     the latest among this set.
c     
      else if (input .eq. fac_sky) then

         query_ttg_cat.archive_id = archive(:archive_len)
         query_ttg_cat.dataset_name = dataset(:dataset_len)

         call ct_gmt_to_binary(jstart,query_ttg_cat.start_time)
         call ct_gmt_to_binary(jstop,query_ttg_cat.stop_time)

         status = cct_query_ttg_catalog(query_ttg_cat,cats,num_cats)

         if (num_cats .eq. 0) then
            fic_open_input = %loc(fic_nocatrecs)
            call lib$signal(fic_nocatrecs,%val(2),archive(:archive_len),
     .                      dataset(:dataset_len))
            return
         end if

         activity(1) = 0
         activity(2) = 0

         do rec = 1,num_cats
            if (time_gt(cats(rec).activity_time,activity)) then
               activity(1) = cats(rec).activity_time(1)
               activity(2) = cats(rec).activity_time(2)
               file_ext = cats(rec).filename_extension
            end if
         end do
c
c     Assign the correct filename for the short science skymap.
c
         call str$trim(file_ext,file_ext,file_ext_len)
         filename = archive(:archive_len)//':'//dataset(:dataset_len)//
     .              '.'//file_ext(:file_ext_len)

      else
c
c     For calibration or raw input, append the jstart and jstop times to the
c     end of the filename for Cobetrieve to handle.
c
         filename = archive(:archive_len)//':'//dataset(:dataset_len)//
     .              '/'//jstart//';'//jstop//';'

      end if
c
c     A separate pair of times must be established for opening the science and
c     engineering files for two reasons.  A short science skymap can refer to
c     data outside the time range established by its catalog start and stop 
c     times, and the last engineering record in a time range can refer to a 
c     science record whose time is after the jstop of that time range.
c
      sky_jstart = jstart
      sky_jstop = jstop
c
c     Open short science file(s) if input is sky or calibration.
c
      if (input .ne. fac_raw) then

         status = fut_get_lun(short_lun)
         if (.not. status) then
            fic_open_input = status
            return
         end if
c
c     For sky input, open the short science skymap and establish the jstart and
c     jstop times for opening the science and engineering file(s) to be the
c     entire mission time period.
c
         if (input .eq. fac_sky) then

            open (unit=short_lun,file=filename,status='old',iostat=status,
     .            form='unformatted',recordtype='fixed',readonly,
     .            useropen=csa_open_skymap)
            if (status .ne. 0) then
               fic_open_input = %loc(fic_csaopen)
               call str$trim(filename,filename,file_len)
               call lib$signal(fic_csaopen,%val(2),filename(:file_len),
     .                         %val(status))
               return
            end if

            sky_jstart = '89324000000000'
            sky_jstop = '90264235959990'

         else
c
c     Open the short science file(s) for calibration data.
c
            open (unit=short_lun,file=filename,status='old',iostat=status,
     .            useropen=ct_connect_read)
            if (status .ne. 0) then
               fic_open_input = %loc(fic_ctopen)
               call str$trim(filename,filename,file_len)
               call lib$signal(fic_ctopen,%val(2),filename(:file_len),
     .                         %val(status))
               return
            end if
         end if
      end if
c
c     Open the engineering file(s).
c
      status = fut_get_lun(eng_lun)
      if (.not. status) then
         fic_open_input = status
         return
      end if
c
c     For raw data input, open the engineering file(s) for a query or a simple
c     read, depending on the presence of a specified rse.  Set the jstop time
c     for opening the science file(s) to be the end of the mission time period.
c
      if (input .eq. fac_raw) then

         sky_jstop = '90264235959990'

         if (rse_input .eq. fac_present) then

             open(unit=eng_lun,file=filename,status='old',iostat=status,
     .            useropen=ct_connect_query)
            if (status .ne. 0) then
               fic_open_input = %loc(fic_ctopen)
               call str$trim(filename,filename,file_len)
               call lib$signal(fic_ctopen,%val(2),filename(:file_len),
     .                         %val(status))
               return
            end if

         else

            open (unit=eng_lun,file=filename,status='old',iostat=status,
     .            useropen=ct_connect_read)
            if (status .ne. 0) then
               fic_open_input = %loc(fic_ctopen)
               call str$trim(filename,filename,file_len)
               call lib$signal(fic_ctopen,%val(2),filename(:file_len),
     .                         %val(status))
               return
            end if

         end if

      else
c
c     Open the engineering file(s) for sky or calibration input.
c
         filename = 'CSDR$FIRAS_RAW:FDQ_ENG'//'/'//sky_jstart//';'//
     .              sky_jstop//';'

         open (unit=eng_lun,file=filename,status='old',iostat=status,
     .         useropen=ct_connect_read)
         if (status .ne. 0) then
            fic_open_input = %loc(fic_ctopen)
            call str$trim(filename,filename,file_len)
            call lib$signal(fic_ctopen,%val(2),filename(:file_len),%val(status))
            return
         end if
      end if
c
c     Open the full science file(s) using the jstart and jstop determined
c     earlier.
c
      status = fut_get_lun(sci_lun)
      if (.not. status) then
         fic_open_input = status
         return
      end if

      filename = 'CSDR$FIRAS_RAW:FDQ_SDF_'//fac_channel_ids(channel)//
     .            '/'//sky_jstart//';'//sky_jstop//';'

      open (unit=sci_lun,file=filename,status='old',iostat=status,
     .      useropen=ct_connect_read)
      if (status .ne. 0) then
         fic_open_input = %loc(fic_ctopen)
         call str$trim(filename,filename,file_len)
         call lib$signal(fic_ctopen,%val(2),filename(:file_len),%val(status))
         return
      end if
c
c     Assign the default output filename extension if none was specified on the
c     command line.
c
      if (outfile .eq. fac_not_present) then
         outfile_ext = 'ED_'//jstart(1:7)//'_'//jstop(1:7)
      end if

      return
      end
