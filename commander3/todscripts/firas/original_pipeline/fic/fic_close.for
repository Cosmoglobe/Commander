      integer*4 function fic_close(input,short_lun,sci_lun,eng_lun,output_luns,
     .                             report_luns)
c----------------------------------------------------------------------------
c
c     Purpose: Close all opened files for a channel.
c
c     Author: S. Alexander, STX, 9/91, SER 7985
c
c     Input: input              i*4  Type of input: sky, calibration, or raw.
c            short_lun          i*4  Logical unit for short science file(s).
c            sci_lun            i*4  Logical unit for science file(s).
c            eng_lun            i*4  Logical unit for engineering file(s).
c            output_luns        i*4(7,4)  Logical units for output by output
c                                         type, including covariance matrices,
c                                         and scan mode.
c            report_luns        i*4(4)  Logical units for reports by scan mode.
c
c     Output:
c
c     Modifications:
c
c----------------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_params)'
      include 'csdr$library:ctparams.inc'
c
c     Return statuses.
c
      external fic_normal
      external fic_csaclose
      external fic_ctclose
      external fic_repclose
c
c     Functions.
c
      integer*4 csa_close_skymap
      integer*4 fut_free_lun
c
c     Input parameters.
c
      integer*4 input,short_lun,sci_lun,eng_lun
      integer*4 output_luns(7,4),report_luns(4)
c
c     Local variables.
c
      integer*4 status,output_type,mode
      character*72 filename
      integer*4 file_len
      integer*2 ct_status(20)

      fic_close = %loc(fic_normal)
c
c     Close files for sky data.
c
      if (input .eq. fac_sky) then
c
c     Close input short science skymap.
c
         status = csa_close_skymap(short_lun,fac_skymap_no_levels)
         if (.not. status) then
            fic_close = %loc(fic_csaclose)
            inquire(short_lun,name=filename)
            call str$trim(filename,filename,file_len)
            call lib$signal(fic_csaclose,%val(2),filename(:file_len),
     .                      %val(status))
            return
         end if

         status = fut_free_lun(short_lun)
         if (.not. status) then
            fic_close = status
            return
         end if
         short_lun = 0
c
c     Close any opened output skymaps for all scan modes.
c
         do output_type = 1,6
            do mode = 1,4
               if (output_luns(output_type,mode) .ne. 0) then
                  status = csa_close_skymap(output_luns(output_type,mode),
     .                                      fac_skymap_no_levels)
                  if (.not. status) then
                     fic_close = %loc(fic_csaclose)
                     inquire(output_luns(output_type,mode),name=filename)
                     call str$trim(filename,filename,file_len)
                     call lib$signal(fic_csaclose,%val(2),filename(:file_len),
     .                               %val(status))
                     return
                  end if
                  status = fut_free_lun(output_luns(output_type,mode))
                  if (.not. status) then
                     fic_close = status
                     return
                  end if
                  output_luns(output_type,mode) = 0
               end if
            end do
         end do
c
c     Close covariance matrices if they have been opened.
c
         output_type = 7
         do mode = 1,4
            if (output_luns(output_type,mode) .ne. 0) then
               call ct_close_arcv(,output_luns(output_type,mode),ct_status)
               if (ct_status(1) .ne. ctp_normal) then
                  fic_close = %loc(fic_ctclose)
                  inquire(output_luns(output_type,mode),name=filename)
                  call str$trim(filename,filename,file_len)
                  call lib$signal(fic_ctclose,%val(2),filename(:file_len),
     .                            %val(ct_status(1)))
                  return
               end if
               status = fut_free_lun(output_luns(output_type,mode))
               if (.not. status) then
                  fic_close = status
                  return
               end if
               output_luns(output_type,mode) = 0
            end if
         end do

      else
c
c     Close files for calibration data.
c
         if (input .eq. fac_cal) then
c
c     Close short science file.
c
            call ct_close_arcv(,short_lun,ct_status)
            if (ct_status(1) .ne. ctp_normal) then
               fic_close = %loc(fic_ctclose)
               inquire(short_lun,name=filename)
               call str$trim(filename,filename,file_len)
               call lib$signal(fic_ctclose,%val(2),filename(:file_len),
     .                         %val(ct_status(1)))
               return
            end if
            status = fut_free_lun(short_lun)
            if (.not. status) then
               fic_close = status
               return
            end if
            short_lun = 0
         end if
c
c     Close any opened output files, including covariance matrices, for all
c     scan modes.
c
         do output_type = 1,7
            do mode = 1,4
               if (output_luns(output_type,mode) .ne. 0) then
                  call ct_close_arcv(,output_luns(output_type,mode),ct_status)
                  if (ct_status(1) .ne. ctp_normal) then
                     fic_close = %loc(fic_ctclose)
                     inquire(output_luns(output_type,mode),name=filename)
                     call str$trim(filename,filename,file_len)
                     call lib$signal(fic_ctclose,%val(2),filename(:file_len),
     .                               %val(ct_status(1)))
                     return
                  end if
                  status = fut_free_lun(output_luns(output_type,mode))
                  if (.not. status) then
                     fic_close = status
                     return
                  end if
                  output_luns(output_type,mode) = 0
               end if
            end do
         end do
      end if
c
c     Close the full science file(s).
c
      call ct_close_arcv(,sci_lun,ct_status)
      if (ct_status(1) .ne. ctp_normal) then
         fic_close = %loc(fic_ctclose)
         inquire(sci_lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fic_ctclose,%val(2),filename(:file_len),
     .                   %val(ct_status(1)))
         return
      end if
      status = fut_free_lun(sci_lun)
      if (.not. status) then
         fic_close = status
         return
      end if
      sci_lun = 0
c
c     Close the enginering file(s).
c
      call ct_close_arcv(,eng_lun,ct_status)
      if (ct_status(1) .ne. ctp_normal) then
         fic_close = %loc(fic_ctclose)
         inquire(eng_lun,name=filename)
         call str$trim(filename,filename,file_len)
         call lib$signal(fic_ctclose,%val(2),filename(:file_len),
     .                   %val(ct_status(1)))
         return
      end if
      status = fut_free_lun(eng_lun)
      if (.not. status) then
         fic_close = status
         return
      end if
      eng_lun = 0
c
c     Close any opened report files for all scan modes.
c
      do mode = 1,4
         if (report_luns(mode) .ne. 0) then
            close(report_luns(mode),iostat=status)
            if (status .ne. 0) then
               fic_close = %loc(fic_repclose)
               inquire(report_luns(mode),name=filename)
               call str$trim(filename,filename,file_len)
               call lib$signal(fic_repclose,%val(2),filename(:file_len),
     .                         %val(status))
               return
            end if
            status = fut_free_lun(report_luns(mode))
            if (.not. status) then
               fic_close = status
               return
            end if
            report_luns(mode) = 0
         end if
      end do

      return
      end
