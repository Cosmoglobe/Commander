      integer*4 function fic_open_output(input,channel,scan_mode,outfile_ext,
     .                                   output_type,lun,output_filenames)
c----------------------------------------------------------------------------
c
c     Purpose: Open an output file.
c
c     Author: S. Alexander, STX, 9/91, SER 7985
c
c     Input: input              i*4  Type of input: sky, calibration, or raw.
c            channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode          i*4  Value of scan mode, 1-4=SS-LF.
c            outfile_ext        ch*72  Value of output filename extension on
c                                      command line.
c            output_type        i*4  Type of output data, 1-6; or 7 for
c                                    covariance matrices.
c            output_filenames   ch*72(4,4,7)  Output filenames by channel,
c                                             scan mode, and output type, 
c                                             including covariance matrices.
c
c     Output: lun                i*4  Logical unit for output.
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
c
c     Return statuses.
c
      external fic_normal
      external fic_csaopen,fic_ctopen
c
c     External references.
c
      external csa_open_skymap
      external ct_connect_write
c
c     Functions.
c
      integer*4 fut_get_lun
      integer*4 csa_field_offset_values
c
c     Input parameters.
c
      integer*4 input,channel,scan_mode,output_type
      character*72 outfile_ext
c
c     Input/output parameters.
c
      character*72 output_filenames(4,4,7)
c
c     Output parameters.
c
      integer*4 lun
c    
c     Local variables.
c
      character*16 archive
      integer*4 archive_len,outfile_ext_len
      character*3 output_ext
      character*72 filename
      integer*4 file_len
      integer*4 status,rec_size,pix_offset

      fic_open_output = %loc(fic_normal)
c
c     Assign the output archive and filename based on the input and output
c     types.
c
      archive = 'CSDR$FIRAS_OUT'
      archive_len = 14

      if (input .eq. fac_sky) then
         output_ext = 'SKY'
      else
         output_ext = 'CAL'
      end if

      if (output_type .eq. fac_out_template) then
         output_ext(2:3) = 'PT'
      else if (output_type .eq. fac_out_sec_template) then
         output_ext(2:3) = 'ST'
      else if (output_type .eq. fac_out_deglitch) then
         output_ext(2:3) = 'DG'
      else if (output_type .eq. fac_out_con_check) then
         output_ext(2:3) = 'CC'
      else if (output_type .eq. fac_out_coadd) then
         output_ext(2:3) = 'CO'
c
c     Assign archive and output filename for covariance matrices.
c
      else if (output_type .eq. fac_out_covar) then
         archive = 'CSDR$FIRAS_STAT'
         archive_len = 15
         output_ext = 'COV'
      end if

      call str$trim(outfile_ext,outfile_ext,outfile_ext_len)

      filename = archive(:archive_len)//':'//'FIC_'//output_ext//'_'//
     .           fac_channel_ids(channel)//fac_scan_mode_ids(scan_mode)//
     .           '.'//outfile_ext(:outfile_ext_len)
c
c     Record the filename for use in reports.
c
      output_filenames(channel,scan_mode,output_type) = filename
c
c     Open output file.
c
      status = fut_get_lun(lun)
      if (.not. status) then
         fic_open_output = status
         return
      end if
c
c     For output skymaps, determine the record size and pixel offset.  All
c     output skymaps have the same record size and pixel offset except for
c     consistency check skymaps.
c   
      if (input .eq. fac_sky) then
  
         if (output_type .eq. fac_out_con_check) then
            rec_size = fac_con_check_size/4
            pix_offset = fac_con_check_pix_offset
         else
            rec_size = fac_coad_spec_size/4
            pix_offset = fac_coad_spec_pix_offset
         end if

         open(unit=lun,file=filename,status='new',iostat=status,recl=rec_size,
     .        form='unformatted',recordtype='fixed',useropen=csa_open_skymap)
         if (status .ne. 0) then
            fic_open_output = %loc(fic_csaopen)
            call str$trim(filename,filename,file_len)
            call lib$signal(fic_csaopen,%val(2),filename(:file_len),
     .                      %val(status))
            return
         end if
c
c     Assign the pixel and time offsets in the output skymap.
c
         status = csa_field_offset_values(pix_offset,fac_time_offset,-1,lun)

      else
c
c     Open output file for calibration or raw time-ordered data.
c   
         open(unit=lun,file=filename,status='new',iostat=status,
     .        useropen=ct_connect_write)
         if (status .ne. 0) then
            fic_open_output = %loc(fic_ctopen)
            call str$trim(filename,filename,file_len)
            call lib$signal(fic_ctopen,%val(2),filename(:file_len),%val(status))
            return
         end if
      end if

      return
      end
