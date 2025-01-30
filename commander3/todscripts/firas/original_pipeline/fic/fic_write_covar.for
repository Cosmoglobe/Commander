      integer*4 function fic_write_covar(channel,outfile_ext,cov_recs,
     .                                   first_covar_times,last_covar_times,
     .                                   output_luns,num_coadds,
     .                                   output_filenames)
c-------------------------------------------------------------------------------
c
c     Purpose: Write covariance matrices.
c
c     Author: S. Alexander, STX, 9/91, SER 7985
c
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            outfile_ext        ch*72  Value of output filename extension on
c                                      command line or default.
c            cov_recs           rec(771)  Covariance matrix records.
c            first_covar_times  i*4(2,3)  Binary times of first records
c                                         contributing to covariance matrices.
c            last_covar_times   i*4(2,3)  Binary times of last records
c                                         contributing to covariance matrices.
c            output_luns        i*4(7,4)  Logical units for output by output
c                                         type, including covariance matrices,
c                                         and scan mode.
c            num_coadds         i*4(4,4)  Number of output coadd records by
c                                         channel and scan mode.
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
c-------------------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_params)'
c
c     Return statuses.
c
      external fic_normal
c
c     Functions.
c
      integer*4 fic_write
      integer*4 cct_set_ttg_time_range
c
c     Input parameters.
c
      integer*4 channel,num_coadds(4,4)
      character*72 outfile_ext

      dictionary 'fic_cov'
      record /fic_cov/ cov_recs(771)

      integer*4 first_covar_times(2,3),last_covar_times(2,3)
c
c     Input/output parameters.
c
      integer*4 output_luns(7,4)
      character*72 output_filenames(4,4,7)
c
c     Local variables.
c
      integer*4 mode,temp_mode,offset,rec,status

      dictionary 'fic_scc'
      dictionary 'fic_sky'
      record /fic_scc/ con_rec
      record /fic_sky/ coa_rec

      integer*4 first_time(2),last_time(2)

      fic_write_covar = %loc(fic_normal)
c
c     Loop over scan modes SS, SF, and LF.
c
      do mode = 1,3
         temp_mode = mode
         if (temp_mode .eq. 3) then
            temp_mode = 4
         end if
         offset = (mode-1)*257
c
c     If any coadds were produced for this channel and scan mode, then the
c     covariance matrix was updated and needs to be written.
c
         if (num_coadds(channel,temp_mode) .ne. 0) then
            do rec = 1+offset,257+offset
               status = fic_write(fac_raw,channel,temp_mode,outfile_ext,con_rec,
     .                            coa_rec,cov_recs(rec),fac_out_covar,
     .                            output_luns,output_filenames)
               if (status .ne. %loc(fic_normal)) then
                  fic_write_covar = status
                  return
               end if
            end do
c
c     Set the time tags for the covariance matrices in the catalog.
c
            first_time(1) = first_covar_times(1,mode)
            first_time(2) = first_covar_times(2,mode)
            last_time(1) = last_covar_times(1,mode)
            last_time(2) = last_covar_times(2,mode)

            status=cct_set_ttg_time_range(output_luns(fac_out_covar,temp_mode),
     .                                    first_time,last_time)
         end if
      end do

      return
      end
