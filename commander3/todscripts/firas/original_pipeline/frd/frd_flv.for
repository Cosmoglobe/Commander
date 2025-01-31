      program frd_flv
c-------------------------------------------------------------------------------
c
c     Purpose: Produce all-sky average variance vectors from coadded
c              interferograms for use in calibrating spectra.
c
c     Author: S. Brodd, HSTX, 10/95, SPR 12286.
c
c     Modifications:
c
c-------------------------------------------------------------------------------
      implicit none
c
c     Include files.
c
      include 'csdr$library:ctuser.inc'
      include '(fut_params)'
c
c     Return statuses.
c
      external frd_normal,frd_failure,frd_ctinit
      external frd_rmsopen,frd_rmswrite,frd_rmsclose,fut_normal
c
c     Functions.
c
      integer*4 cut_register_version,cut_display_banner
      integer*4 frd_parse_flv,frd_read_flv,fut_get_lun,fut_free_lun
c
c     Local variables for frd_parse_flv.
c
      integer*4 channel,scan_mode,galexc,min_ifgs,label,infile_num
      real*4 galexc_val
      character*40 label_val
      character*48 infile(20),outfile
c
c     Local variables for frd_read_flv.
c
      integer*4 num_ifgs,deg_freedom
      real*4 adj_num_ifgs
      real*8 var(3,361)
c
c     Local variables.
c
      integer*4 current_time(2),status,sec_status,var_num,lun,length
      character*14 current_gmt
      character*6 version
      parameter(version='13.8')
      integer*2 ct_status(20)

      dictionary 'fex_flv'
      record /fex_flv/ out_rec
c
c     Get the invocation time.
c
      call sys$gettim(current_time)
      call ct_binary_to_gmt(current_time,current_gmt)
c
c     Register version and display banner.
c
      status = cut_register_version(version)
      status = cut_display_banner(6,80,'FIRAS Facility FRD_FLV')
c
c     Parse the command line.
c
      status = frd_parse_flv(channel,scan_mode,galexc,galexc_val,min_ifgs,
     .                       label,label_val,infile,infile_num,outfile)
      if (status .eq. %loc(frd_normal)) then
c
c     Initialize Cobetrieve.
c
         call ct_init(ct_status)
         if (ct_status(1) .ne. ctp_normal) then
            status = %loc(frd_ctinit)
            call lib$signal(frd_ctinit,%val(1),%val(ct_status(1)))
         else
c
c     Read the variances from the FIL records.
c
            status = frd_read_flv(galexc,galexc_val,min_ifgs,infile,infile_num,
     .                            num_ifgs,adj_num_ifgs,deg_freedom,var)
            if (status .eq. %loc(frd_normal)) then
c
c     Put information into the output record.
c
               out_rec.gmt = current_gmt
               out_rec.time(1) = current_time(1)
               out_rec.time(2) = current_time(2)
               out_rec.channel = channel
               out_rec.scan_mode = scan_mode
               if (label .eq. fac_present) then
                  out_rec.label = label_val
               end if
               out_rec.galat_exc = galexc_val
               out_rec.min_ifgs = min_ifgs
               out_rec.nsky_ifgs = float(num_ifgs)
               out_rec.adj_nsky_ifgs = adj_num_ifgs
               out_rec.deg_freedom = deg_freedom
               do var_num = 1,361
                  out_rec.rr_variances(var_num) = var(1,var_num)
                  out_rec.ii_variances(var_num) = var(2,var_num)
                  out_rec.ri_variances(var_num) = var(3,var_num)
               end do
c
c     Open the FEX_FLV file and write the record to it.
c
               sec_status = fut_get_lun(lun)
               if (sec_status .ne. %loc(fut_normal)) then
                  call lib$signal(%val(sec_status))
               end if

               open(unit=lun,file=outfile,status='new',form='unformatted',
     .              recordtype='fixed',recl=2240,access='sequential',
     .              iostat=sec_status)

               if (sec_status .eq. 0) then
                  write(lun,iostat=sec_status) out_rec
                  if (sec_status .ne. 0) then
                     sec_status = %loc(frd_rmswrite)
                     call str$trim(outfile,outfile,length)
                     call lib$signal(frd_rmswrite,%val(2),outfile(1:length),
     .               %val(sec_status))
                  end if
c
c     Close the file.
c
                  close(lun,iostat = sec_status)
                  if (sec_status .ne. 0) then
                     sec_status = %loc(frd_rmsclose)
                     call str$trim(outfile,outfile,length)
                     call lib$signal(frd_rmsclose,%val(2),outfile(1:length),
     .                               %val(sec_status))
                  end if
               else
                  sec_status = %loc(frd_rmsopen)
                  call str$trim(outfile,outfile,length)
                  call lib$signal(frd_rmsopen,%val(2),outfile(1:length),
     .                            %val(sec_status))
               end if
c
c     Free logical unit number.
c
               sec_status = fut_free_lun(lun)
               if (sec_status .ne. %loc(fut_normal)) then
                  call lib$signal(%val(sec_status))
               end if
            end if
         end if
      end if
c
c     Write out the number of skymaps and interferograms processed.
c
      if (status .eq. %loc(frd_normal)) then
         write(6,10,iostat=sec_status) infile_num,num_ifgs
      end if
 10   format(/,4x,'Number of skymaps processed:          ',i6,
     .       /,4x,'Number of interferograms represented: ',i6,//)
c
c     Signal the program completion status.
c
      if (status .eq. %loc(frd_normal)) then
         call lib$signal(frd_normal)
      else
         call lib$signal(frd_failure)
      end if

      end
