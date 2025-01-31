      integer*4 function frd_read_flv(galexc,galexc_val,min_ifgs,infile,
     .                                infile_num,num_ifgs,adj_num_ifgs,
     .                                deg_freedom,var)
c-------------------------------------------------------------------------------
c
c     Purpose: Read in all coadds from specified FIL skymaps and accumulate
c              three weighted average variance vectors from those coadds
c              passing galactic latitude exclusion and minimum interferograms
c              criteria.
c
c     Author: S. Brodd, HSTX, 10/95.
c
c     Input: galexc        integer*4         Galactic exclusion requested.
c            galexc_val    real*4            Value of galactic exclusion.
c            min_ifgs      integer*4         Minimum number of interferograms.
c            infile        character*48(20)  Input FIL filenames.
c            infile_num    integer*4         Number of input FIL filenames.
c
c     Output: num_ifgs      integer*4         Number of interferograms
c                                             contributing to average variances.
c             adj_num_ifgs  real*4            Glitch rate weighted number of
c                                             contributing interferograms.
c             deg_freedom   integer*4         Degrees of freedom.
c             var           real*8(3,361)     Three average variance vectors.
c
c-------------------------------------------------------------------------------
      implicit none
c
c     Include files.
c
      include '(fut_params)'
      include '(csa_pixel_input_rec)'
      include '(csa_pixel_output_rec)'
c
c     Return statuses.
c
      external frd_normal,fut_normal,csa_normal
      external frd_csaopen,frd_csaread,frd_csaclose
c
c     External references.
c
      external csa_open_skymap
c
c     Functions.
c
      integer*4 frd_gal_cut_flv,fut_get_lun,fut_free_lun
      integer*4 csa_open_skymap,csa_read_pixels,csa_close_skymap
c
c     Input parameters.
c
      integer*4 galexc,min_ifgs,infile_num
      real*4 galexc_val
      character*48 infile(infile_num)
c
c     Output parameters.
c
      integer*4 num_ifgs,deg_freedom
      real*4 adj_num_ifgs
      real*8 var(3,361)
c
c     Local variables.
c
      integer*4 status,lun,file_num,length,pixel,num_out
      integer*4 rec_num,rec_num_ifgs,rec_mult,var_num,vec_num
      real*4 rec_adj_num_ifgs

      record /pixel_input_list/  inlist
      record /pixel_output_list/ outlist

      dictionary 'fil_sky'
      record /fil_sky/ recs(fac_max_coad)

      frd_read_flv = %loc(frd_normal)
c
c     Get logical unit number.
c
      status = fut_get_lun(lun)
      if (status .ne. %loc(fut_normal)) then
         frd_read_flv = status
         call lib$signal(%val(status))
         return
      end if
c
c     Read the records for each FIL skymap.
c
      num_ifgs = 0
      adj_num_ifgs = 0.0
      deg_freedom = 0
      call lib$movc5(0,,0,8664,var)
      file_num = 0

      do while (file_num .lt. infile_num)
         file_num = file_num + 1
         open(unit=lun,file=infile(file_num),status='old',iostat=status,
     .        form='unformatted',recordtype='fixed',readonly,
     .        useropen=csa_open_skymap)
         if (status .ne. 0) then
            frd_read_flv = %loc(frd_csaopen)
            call str$trim(infile(file_num),infile(file_num),length)
            call lib$signal(frd_csaopen,%val(2),infile(file_num)(1:length),
     .                      %val(status))
            return
         end if
c
c     Read the pixels.
c
         inlist.level_no = fac_skymap_level
         pixel = -1

         do while (pixel .lt. 6143)
            pixel = pixel + 1
            if (galexc .eq. fac_present) then
               status = frd_gal_cut_flv(pixel,galexc_val)
               if (status .ne. %loc(frd_normal)) then
                  frd_read_flv = status
                  call lib$signal(%val(status))
                  return
               end if
            end if

            inlist.pixel_no = pixel
            status = csa_read_pixels(lun,inlist,1,recs,fac_max_coad,outlist,
     .                               num_out,0)
            if (status .ne. %loc(csa_normal)) then
               frd_read_flv = %loc(frd_csaread)
               call str$trim(infile(file_num),infile(file_num),length)
               call lib$signal(frd_csaread,%val(2),infile(file_num)(1:length),
     .                         %val(status))
               return
            end if

            if (outlist.no_records .gt. 0) then
c
c     Sum the number of interferograms, glitch rate weighted number, degrees of
c     freedom, and three variance vectors for those coadds with at least
c     min_ifgs in them.
c
               do rec_num = 1,outlist.no_records
                  rec_num_ifgs = recs(rec_num).coad_spec_head.num_ifgs
                  if (rec_num_ifgs .ge. min_ifgs) then
                     num_ifgs = num_ifgs + rec_num_ifgs
                     rec_adj_num_ifgs = 
     .                   recs(rec_num).coad_spec_head.adj_num_ifgs
                     adj_num_ifgs = adj_num_ifgs + rec_adj_num_ifgs
                     deg_freedom = deg_freedom + (rec_num_ifgs - 1)
                     rec_mult = rec_num_ifgs - 1 - 
     .                   recs(rec_num).coad_spec_data.sec_template.subtracted
                     do var_num = 1,361
                        var(1,var_num) = var(1,var_num) +
     .                     dble(recs(rec_num).coad_data.real_var(var_num)) *
     .                     rec_mult * rec_adj_num_ifgs
                     end do
                     do var_num = 1,361
                        var(2,var_num) = var(2,var_num) +
     .                     dble(recs(rec_num).coad_data.imag_var(var_num)) *
     .                     rec_mult * rec_adj_num_ifgs
                     end do
                     do var_num = 1,361
                        var(3,var_num) = var(3,var_num) +
     .                     dble(recs(rec_num).coad_data.real_imag_var(var_num))*
     .                     rec_mult * rec_adj_num_ifgs
                     end do
                  end if
               end do
            end if
         end do
c
c     Close the input skymap file.
c
         status = csa_close_skymap(lun,fac_skymap_no_levels)
         if (status .ne. %loc(csa_normal)) then
            frd_read_flv = %loc(frd_csaclose)
            call str$trim(infile(file_num),infile(file_num),length)
            call lib$signal(frd_csaclose,%val(2),infile(file_num)(1:length),
     .                      %val(status))
            return
         end if
      end do
c
c     Divide the variances vectors by the degrees of freedom.
c
      do vec_num = 1,3
         do var_num = 1,361
            var(vec_num,var_num) = var(vec_num,var_num) / deg_freedom
         end do
      end do
c
c     Free logical unit number.
c
      status = fut_free_lun(lun)
      if (status .ne. %loc(fut_normal)) then
         frd_read_flv = status
         call lib$signal (%val(status))
         return
      end if

      return
      end
