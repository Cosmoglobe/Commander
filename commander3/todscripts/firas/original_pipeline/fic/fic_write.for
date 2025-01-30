      integer*4 function fic_write(input,channel,scan_mode,outfile_ext,con_rec,
     .                             coa_rec,cov_rec,output_type,output_luns,
     .                             output_filenames)
c----------------------------------------------------------------------------
c
c     Purpose: Write a coadd, consistency check, or covariance matrix record
c              to an output file.
c
c     Author: S. Alexander, STX, 9/91, SER 7985
c
c     Input: input              i*4  Type of input: sky, calibration, or raw.
c            channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode          i*4  Value of scan mode, 1-4=SS-LF.
c            outfile_ext        ch*72  Value of output filename extension on
c                                      command line or default.
c            con_rec            rec  Consistency check record.
c            coa_rec            rec  Coadd record.
c            cov_rec            rec  Covariance matrix record.
c            output_type        i*4  Type of output data, 1-6; or 7 for
c                                    covariance matrices.
c            output_luns        i*4(7,4)  Logical units for output by output
c                                         type, including covariance matrices,
c                                         and scan mode.
c            output_filenames   ch*72(4,4,7)  Output filenames by channel,
c                                             scan mode, and output type, 
c                                             including covariance matrices.
c
c     Output: output_luns        i*4(7,4)  Logical units for output by output
c                                          type, including covariance matrices,
c                                          and scan mode.
c             output_filenames   ch*72(4,4,7)  Output filenames by channel,
c                                              scan mode, and output type, 
c                                              including covariance matrices.
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
      external fic_csawrite,fic_ctwrite
c
c     Functions.
c
      integer*4 fic_open_output
      integer*4 csa_write_pixels
c
c     Input parameters.
c
      integer*4 input,channel,scan_mode
      character*72 outfile_ext

      dictionary 'fic_scc'
      dictionary 'fic_sky'
      dictionary 'fic_cov'

      record /fic_scc/ con_rec
      record /fic_sky/ coa_rec
      record /fic_cov/ cov_rec
c
c     Input/output parameters.
c
      integer*4 output_type,output_luns(7,4)
      character*72 output_filenames(4,4,7)
c
c     Local variables.
c
      integer*4 lun,status
      integer*2 ct_status(20)
      character*72 filename
      integer*4 file_len

      fic_write = %loc(fic_normal)
c
c     If an output file for this output type, channel, and scan mode has not
c     yet been opened, call fic_open_output to do so.
c
      lun = output_luns(output_type,scan_mode)

      if (lun .eq. 0) then

         status = fic_open_output(input,channel,scan_mode,outfile_ext,
     .                            output_type,lun,output_filenames)
         if (status .ne. %loc(fic_normal)) then
            fic_write = status
            return
         end if

         output_luns(output_type,scan_mode) = lun
      end if
c
c     Write record to the output file.
c
      if (input .eq. fac_sky) then
c
c     Write consistency check record to output skymap.
c
         if (output_type .eq. fac_out_con_check) then
            status = csa_write_pixels(lun,con_rec,1,0)
            if (.not. status) then
               fic_write = %loc(fic_csawrite)
               filename = output_filenames(channel,scan_mode,output_type)
               call str$trim(filename,filename,file_len)
               call lib$signal(fic_csawrite,%val(2),filename(:file_len),
     .                         %val(status))
               return
            end if
         else
c
c     Write coadd record to output skymap.
c
            status = csa_write_pixels(lun,coa_rec,1,0)
            if (.not. status) then
               fic_write = %loc(fic_csawrite)
               filename = output_filenames(channel,scan_mode,output_type)
               call str$trim(filename,filename,file_len)
               call lib$signal(fic_csawrite,%val(2),filename(:file_len),
     .                         %val(status))
               return
            end if
         end if
      else
c
c     Write consistency check record to output time-ordered file.
c
         if (output_type .eq. fac_out_con_check) then
            call ct_write_arcv(,lun,con_rec,ct_status)
            if (ct_status(1) .ne. ctp_normal) then
               fic_write = %loc(fic_ctwrite)
               filename = output_filenames(channel,scan_mode,output_type)
               call str$trim(filename,filename,file_len)
               call lib$signal(fic_ctwrite,%val(2),filename(:file_len),
     .                         %val(ct_status(1)))
               return
            end if
         else if (output_type .eq. fac_out_covar) then
c
c     Write covariance matrix record.
c
            call ct_write_arcv(,lun,cov_rec,ct_status)
            if (ct_status(1) .ne. ctp_normal) then
               fic_write = %loc(fic_ctwrite)
               filename = output_filenames(channel,scan_mode,output_type)
               call str$trim(filename,filename,file_len)
               call lib$signal(fic_ctwrite,%val(2),filename(:file_len),
     .                         %val(ct_status(1)))
               return
            end if
         else
c
c     Write coadd record to output time-ordered file.
c
            call ct_write_arcv(,lun,coa_rec,ct_status)
            if (ct_status(1) .ne. ctp_normal) then
               fic_write = %loc(fic_ctwrite)
               filename = output_filenames(channel,scan_mode,output_type)
               call str$trim(filename,filename,file_len)
               call lib$signal(fic_ctwrite,%val(2),filename(:file_len),
     .                         %val(ct_status(1)))
               return
            end if
         end if
      end if

      return
      end
