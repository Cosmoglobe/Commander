      program fic

c----------------------------------------------------------------------------
c
c     Purpose: Perform interferogram coaddition and related functions:
c              instrument state and shape consistency checking, primary and
c              secondary template formation, deglitching, baseline subtraction,
c              and computation of covariance matrices and other statistical
c              information.
c
c     Author: S. Alexander, STX, 9/91, SER 7985
c
c     Modifications:
c
c----------------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_params)'
      include '(fut_error)'
      include 'csdr$library:ctparams.inc'
c
c     Return statuses.
c
      external fic_normal
      external fic_ctinit
      external fic_failure
c
c     External references.
c
      external fut_error
c
c     Functions.
c
      integer*4 cut_register_version,cut_display_banner
      integer*4 fic_parse,fic_open_summary,fic_reference
      integer*4 fic_read_covar,fic_open_input,fic_read
      integer*4 fic_qual_check,fic_inst_state,fic_convert
      integer*4 fic_template,fic_sec_template
      integer*4 fic_deglitch,fic_noise
      integer*4 fic_shape_check,fic_coadd,fic_baseline_sub
      integer*4 fic_transient,fic_write,fic_write_report
      integer*4 fic_close,fic_write_covar,fic_close_summary
c
c     Local variables for fic_parse.
c
      integer*4 input
      integer*4 start_channel,stop_channel,skip_channel,scan_mode_input
      integer*4 file,rse_input
      integer*4 covar,covar_new,qual_inst,qual_att,single
      integer*4 inst_state,pixel_no,sec_template_input
      integer*4 deglitch,shape_check,coadd,baseline_sub
      integer*4 output(6)
      integer*4 outfile,report,report_detail
      character*72 file_ext,rse_name,outfile_ext,report_name
      character*14 jstart,jstop
      character*79 command_line(9)
c
c     Local variables for fic_open_summary.
c
      character*12 account
      character*14 current_gmt
c
c     Local variables for fic_read_covar.
c     
      dictionary 'fic_cov'
      record /fic_cov/ cov_recs(771)

      integer*4 channel,first_covar_times(2,3),last_covar_times(2,3)
c
c     Local variables for fic_open_input.
c     
      integer*4 short_lun,sci_lun,eng_lun
c
c     Local variables for fic_reference.
c
      logical*1 init
      character*128 rse(16)

      dictionary 'fex_mincoadd'
      dictionary 'fex_grttrans'
      dictionary 'fex_grtrawwt'
      dictionary 'fex_grtcoawt'
      dictionary 'fex_cmdgain'
      dictionary 'fex_samprate'
      dictionary 'fex_nyquist'
      dictionary 'fex_cth'
      dictionary 'fex_basis'

      record /fex_mincoadd/ mincoadd
      record /fex_grttrans/ grttrans
      record /fex_grtrawwt/ grtrawwt
      record /fex_grtcoawt/ grtcoawt
      record /fex_cmdgain/ cmdgain
      record /fex_samprate/ samprate
      record /fex_nyquist/ nyquist
      record /fex_cth/ cth
      record /fex_basis/ basis

      integer*4 reftemps_lun,dtrf_lun,gltchpro_lun,etf_lun
c
c     Local variables for fic_read.
c
      integer*4 last_channel,num_recs
      logical*1 next_channel

      dictionary 'fdq_sdf'
      dictionary 'fdq_eng'

      record /fdq_sdf/ sci_recs(fac_max_num)
      record /fdq_eng/ eng_recs(fac_max_num)
c
c     Local variables for fic_qual_check.
c
      integer*4 scan_mode,num_sci_recs(4,4),first_sci_rec,last_sci_rec
      integer*4 num_groups(4,4)
      integer*4 num_good,min_good,num_sci_fail(4,4,32) 
      integer*4 fake_it,sci_mode,adds_per_group,xcal_pos,sweeps
      real*4 temps(10,fac_max_num)

      dictionary 'fic_scc'
      dictionary 'fic_sky'

      record /fic_scc/ con_recs(fac_max_num)
      record /fic_sky/ coa_rec
c
c     Local variables for fic_convert.
c
      real*4 aifgs(512,fac_max_num)
c
c     Local variables for fic_template.
c
      real*4 template(512)
c
c     Local variables for fic_sec_template.
c
      real*4 sec_template(512),x(512,2),b(fac_max_num)
      integer*4 num_peak_pts
c
c     Local variables for fic_transient.
c
      logical*1 transient_coadd
      real*4 c(fac_max_num)
c
c     Local variables for fic_deglitch.
c
      integer*4 num_sci_deglitch(4,4)
      integer*4 glitch_pos(fac_max_deglitch,fac_max_num)
      real*4 glitch_amp(fac_max_deglitch,fac_max_num)
      integer*4 glitch_tot_iter(4,4)
      real*4 glitch_tot_signal(4,4)
c
c     Local variables for fic_noise.
c
      integer*4 num_sci_noise(4,4)
      real*4 avg_noise
      real*4 tot_noise(4,4)
c
c     Local variables for fic_coadd.
c
      integer*4 num_sci_good(4,4,20),num_coadds(4,4)
c
c     Local variables for fic_write.
c
      record /fic_scc/ con_rec
      record /fic_cov/ cov_rec
      integer*4 output_type,output_luns(7,4)
      character*72 output_filenames(4,4,7)
c
c     Local variables for fic_write_report.
c
      integer*4 report_luns(4)
      character*72 report_filenames(4,4)
c
c     Local variables.
c
      integer*4 status
      character*6 version
      parameter(version='10.3')
      integer*2 ct_status(20)
      logical*1 write_template,write_sec_template
      integer*4 rec,temp_single,temp_num_recs,temp_num_good,ifg_pos

      record /fic_sky/ temp_coa_rec
   
      real*4 temp_aifgs(512,1)
      integer*4 temp_num_sci_good(4,4,20),temp_num_coadds(4,4)
c
c     Register version and display banner.
c
      status = cut_register_version(version)
      status = cut_display_banner(6,80,'FIRAS Facility FIC_Interferogram_Coadd')
c
c     Parse command line.
c
      status = fic_parse(input,start_channel,stop_channel,skip_channel,
     .                   scan_mode_input,file,file_ext,jstart,jstop,rse_input,
     .                   rse_name,covar,covar_new,qual_inst,qual_att,single,
     .                   inst_state,pixel_no,sec_template_input,deglitch,
     .                   shape_check,coadd,baseline_sub,output,outfile,
     .                   outfile_ext,report,report_detail,report_name,
     .                   command_line)
c
c     Set local variables for writing deglitched interferograms.
c
      temp_single = fac_present
      temp_num_recs = 1
      temp_num_good = 1

      init = .true.
c
c     Initialize summary report.
c
      status = fic_open_summary(file,file_ext,jstart,jstop,rse_input,rse_name,
     .                          report,report_name,command_line,status,account,
     .                          current_gmt)
c
c     Initialize Cobetrieve and error handler.
c
      if (status .eq. %loc(fic_normal)) then
         call lib$establish(fut_error)
         call ct_init(ct_status)
         if (ct_status(1) .ne. ctp_normal) then
            status = %loc(fic_ctinit)  
            call lib$signal(fic_ctinit,%val(1),%val(ct_status(1)))
         end if
      end if
c
c     Process each requested channel.
c
      last_channel = 0
      channel = start_channel
      do while ((status .eq. %loc(fic_normal)) .and. 
     .          (channel .le. stop_channel))
c
c     Read current covariance matrices if they exist; initialize if not.
c
         if (covar .eq. fac_present) then 
            status = fic_read_covar(channel,scan_mode_input,covar_new,
     .                              current_gmt,cov_recs,first_covar_times,
     .                              last_covar_times)
         end if
c
c     Open all input files.
c     
         if (status .eq. %loc(fic_normal)) then
            status = fic_open_input(input,channel,file,file_ext,jstart,jstop,
     .                              rse_input,outfile,outfile_ext,short_lun,
     .                              sci_lun,eng_lun)
         end if
c
c     Open reference data sets.
c
         if ((status .eq. %loc(fic_normal)) .and. init) then
            status = fic_reference(jstart,rse_input,rse_name,init,rse,mincoadd,
     .                             grttrans,grtrawwt,grtcoawt,cmdgain,samprate,
     .                             nyquist,cth,basis,reftemps_lun,dtrf_lun,
     .                             gltchpro_lun,etf_lun)
            init = .false.
         end if
c
c     Process current channel.
c
         next_channel = .false.
         do while ((status .eq. %loc(fic_normal)) .and. (.not. next_channel))
c
c     Initialize consistency check and coadd records.
c
            do rec = 1,fac_max_num
               call lib$movc5(0,,0,fac_con_check_size,con_recs(rec))
            end do
            call lib$movc5(0,,0,fac_coad_spec_size,coa_rec)

            write_template = .false.
            write_sec_template = .false.
            num_good = 0
            min_good = 1
c
c     Read a coaddable group of interferograms.
c 
           status = fic_read(input,channel,last_channel,next_channel,
     .                       scan_mode_input,rse_input,rse,single,short_lun,
     .                       sci_lun,eng_lun,num_recs,sci_recs,eng_recs)
c
c     Perform quality checking.
c
            if ((status .eq. %loc(fic_normal)) .and. (num_recs .gt. 0)) then
               status = fic_qual_check(channel,scan_mode_input,scan_mode,
     .                                 qual_inst,qual_att,single,deglitch,
     .                                 mincoadd,grttrans,grtcoawt,cmdgain,
     .                                 samprate,nyquist,num_recs,sci_recs,
     .                                 eng_recs,num_sci_recs,first_sci_rec,
     .                                 last_sci_rec,num_groups,num_good,
     .                                 min_good,num_sci_fail,fake_it,sci_mode,
     .                                 adds_per_group,xcal_pos,sweeps,temps,
     .                                 con_recs,coa_rec)
            end if
c
c     Perform instrument state consistency checking.
c
            if ((status .eq. %loc(fic_normal)) .and. (num_good .ge. min_good)
     .          .and. (inst_state .eq. fac_present)) then
               status = fic_inst_state(input,channel,scan_mode,pixel_no,
     .                                 mincoadd,grttrans,grtrawwt,cth,num_recs,
     .                                 sci_recs,eng_recs,num_good,min_good,
     .                                 num_sci_fail,xcal_pos,temps,con_recs)
            end if
c
c     Convert the interferograms by dividing by gain and sweeps and subtracting
c     the dither.
c
            if ((status .eq. %loc(fic_normal)) .and. 
     .          (num_good .ge. min_good)) then
               status = fic_convert(channel,num_recs,sci_recs,sweeps,con_recs,
     .                              aifgs)
            end if
c
c     Calculate and subtract primary template.
c
            if ((status .eq. %loc(fic_normal)) .and. (num_good .ge. min_good)
     .          .and. ((deglitch .eq. fac_present) .or.
     .          (shape_check .eq. fac_present))) then
               status = fic_template(num_recs,num_good,con_recs,coa_rec,aifgs,
     .                               template)
               if (status .eq. %loc(fic_normal)) then
                  if (output(fac_out_template) .eq. fac_present) then
                     write_template = .true.
                  end if
c
c     Calculate and subtract secondary template.
c
                  num_peak_pts = cth.peak_points(channel)
                  status = fic_sec_template(input,channel,scan_mode,
     .                                      sec_template_input,cth,reftemps_lun,
     .                                      num_recs,num_good,xcal_pos,con_recs,
     .                                      coa_rec,aifgs,template,sec_template,
     .                                      num_peak_pts,x,b)
                  if (status .eq. %loc(fic_normal)) then
                     if (output(fac_out_sec_template) .eq. fac_present) then
                        write_sec_template = .true.
                     end if
                  end if
               end if
            end if
c
c     Perform digital transient function subtraction on interferograms.
c
            if ((status .eq. %loc(fic_normal)) .and. 
     .          (num_good .ge. min_good)) then
               transient_coadd = .false.
               status = fic_transient(channel,scan_mode,dtrf_lun,num_recs,
     .                                num_good,con_recs,coa_rec,aifgs,
     .                                transient_coadd,c)
            end if
c
c     Perform deglitching.
c
            if ((status .eq. %loc(fic_normal)) .and. (num_good .ge. min_good)
     .          .and. (deglitch .eq. fac_present)) then
               status = fic_deglitch(channel,scan_mode,mincoadd,gltchpro_lun,
     .                               num_recs,num_good,min_good,
     .                               num_sci_deglitch,num_sci_fail,fake_it,
     .                               adds_per_group,con_recs,coa_rec,aifgs,
     .                               glitch_pos,glitch_amp,glitch_tot_iter,
     .                               glitch_tot_signal)
            end if
c
c     Calculate noise of interferograms and coadd.
c
            if ((status .eq. %loc(fic_normal)) .and. (num_good .ge. min_good)
     .          .and. ((deglitch .eq. fac_present) .or.
     .          (shape_check .eq. fac_present))) then
               status = fic_noise(channel,scan_mode,num_recs,num_good,
     .                            num_sci_noise,adds_per_group,con_recs,
     .                            coa_rec,aifgs,avg_noise,tot_noise)
            end if
c
c     Write deglitched interferograms.
c
            if ((status .eq. %loc(fic_normal)) .and. 
     .          (output(fac_out_deglitch) .eq. fac_present)) then
               rec = 1
               do while ((status .eq. %loc(fic_normal)) .and. 
     .                   (rec .le. num_recs))
                  if (con_recs(rec).con_check .eq. 0) then
                     temp_coa_rec = coa_rec
                     do ifg_pos = 1,512
                        temp_aifgs(ifg_pos,1) = aifgs(ifg_pos,rec)
                     end do
c
c     Invoke fic_coadd to translate science and engineering record into a
c     coadd record.
c
                     status = fic_coadd(input,channel,scan_mode,covar,
     .                                  temp_single,grttrans,grtcoawt,cth,
     .                                  etf_lun,temp_num_recs,temp_num_good,
     .                                  temp_num_sci_good,sci_recs(rec),
     .                                  eng_recs(rec),fake_it,sci_mode,
     .                                  adds_per_group,temps,con_recs(rec),
     .                                  temp_num_coadds,temp_coa_rec,
     .                                  temp_aifgs,template,sec_template,
     .                                  cov_recs,first_covar_times,
     .                                  last_covar_times)
                     if (status .eq. %loc(fic_normal)) then
                        temp_coa_rec.ct_head.gmt = sci_recs(rec).ct_head.gmt
                        temp_coa_rec.ct_head.time(1) = 
     .                               sci_recs(rec).ct_head.time(1) 
                        temp_coa_rec.ct_head.time(2) = 
     .                               sci_recs(rec).ct_head.time(2) 
c
c     Write individual b and c values into output coadd record.
c
                        temp_coa_rec.coad_spec_data.sec_template.b_average = 
     .                               b(rec)
                        temp_coa_rec.coad_spec_data.transient.c_average = 
     .                               c(rec)
                        output_type = fac_out_deglitch
                        status = fic_write(input,channel,scan_mode,outfile_ext,
     .                                     con_rec,temp_coa_rec,cov_rec,
     .                                     output_type,output_luns,
     .                                     output_filenames)
                     end if
                  end if
                  rec = rec + 1
               end do
            end if
c
c     Perform shape consistency checking.
c
            if ((status .eq. %loc(fic_normal)) .and. (num_good .ge. min_good)
     .          .and. (shape_check .eq. fac_present)) then
               status = fic_shape_check(channel,scan_mode,cth,num_recs,
     .                                  num_good,min_good,num_sci_fail,
     .                                  con_recs,aifgs,avg_noise)
            end if
c
c     Perform coaddition.
c 
           if ((status .eq. %loc(fic_normal)) .and. (num_good .ge. min_good)
     .          .and. (coadd .eq. fac_present)) then
               status = fic_coadd(input,channel,scan_mode,covar,single,
     .                            grttrans,grtcoawt,cth,etf_lun,num_recs,
     .                            num_good,num_sci_good,sci_recs,eng_recs,
     .                            fake_it,sci_mode,adds_per_group,temps,
     .                            con_recs,num_coadds,coa_rec,aifgs,template,
     .                            sec_template,cov_recs,first_covar_times,
     .                            last_covar_times)
c
c     Write coadd record.
c
               if ((status .eq. %loc(fic_normal)) .and. 
     .             (output(fac_out_coadd) .eq. fac_present)) then
                  output_type = fac_out_coadd
                  status = fic_write(input,channel,scan_mode,outfile_ext,
     .                               con_rec,coa_rec,cov_rec,output_type,
     .                               output_luns,output_filenames)
               end if
            end if
c
c     Write primary template.
c
            if ((status .eq. %loc(fic_normal)) .and. (num_good .ge. min_good)
     .           .and. write_template) then
               temp_coa_rec = coa_rec
               do ifg_pos = 1,512
                  temp_coa_rec.coad_data.ifg(ifg_pos) = template(ifg_pos)
               end do
               output_type = fac_out_template
               status = fic_write(input,channel,scan_mode,outfile_ext,con_rec,
     .                            temp_coa_rec,cov_rec,output_type,output_luns,
     .                            output_filenames)
            end if
c
c     Write secondary template.
c
            if ((status .eq. %loc(fic_normal)) .and. (num_good .ge. min_good)
     .           .and. write_sec_template) then
               temp_coa_rec = coa_rec
               do ifg_pos = 1,512
                  temp_coa_rec.coad_data.ifg(ifg_pos) = sec_template(ifg_pos)
               end do
               output_type = fac_out_sec_template
               status = fic_write(input,channel,scan_mode,outfile_ext,con_rec,
     .                            temp_coa_rec,cov_rec,output_type,output_luns,
     .                            output_filenames)
            end if
c
c     Perform baseline subtraction.
c                  
            if ((status .eq. %loc(fic_normal)) .and. (num_good .ge. min_good)
     .          .and. (baseline_sub .eq. fac_present)) then
               status = fic_baseline_sub(basis,coa_rec)
            end if
c
c     Perform digital transient function subtraction on the coadded
c     interferogram.
c
            if ((status .eq. %loc(fic_normal)) .and. 
     .          (num_good .ge. min_good)) then
               transient_coadd = .true.
               status = fic_transient(channel,scan_mode,dtrf_lun,temp_num_recs,
     .                                temp_num_good,con_recs,coa_rec,aifgs,
     .                                transient_coadd,c)
            end if
c
c     Write baseline-subtracted coadd record.
c
            if ((status .eq. %loc(fic_normal)) .and. (num_good .ge. min_good)
     .          .and. (output(fac_out_baseline_sub) .eq. fac_present)) then
               output_type = fac_out_baseline_sub
               status = fic_write(input,channel,scan_mode,outfile_ext,con_rec,
     .                            coa_rec,cov_rec,output_type,output_luns,
     .                            output_filenames)
            end if
c
c     Write consistency check records.
c
            if ((status .eq. %loc(fic_normal)) .and. (num_recs .gt. 0)
     .          .and. (output(fac_out_con_check) .eq. fac_present)) then
               output_type = fac_out_con_check
               rec = 1
               do while ((status .eq. %loc(fic_normal)) .and. 
     .                   (rec .le. num_recs))
                  status = fic_write(input,channel,scan_mode,outfile_ext,
     .                               con_recs(rec),coa_rec,cov_rec,output_type,
     .                               output_luns,output_filenames)
                  rec = rec + 1
               end do
            end if
c
c     Write report for channel and scan mode.
c
            if ((status .eq. %loc(fic_normal)) .and. (num_recs .gt. 0)
     .          .and. (report .eq. fac_present)) then
               status = fic_write_report(channel,scan_mode,deglitch,
     .                                   report_detail,report_name,account,
     .                                   current_gmt,mincoadd,cth,num_recs,
     .                                   fake_it,sci_mode,adds_per_group,
     .                                   con_recs,first_sci_rec,last_sci_rec,
     .                                   avg_noise,glitch_pos,glitch_amp,
     .                                   report_luns,report_filenames)
            end if
         end do
c
c     Write covariance records.
c
         if ((status .eq. %loc(fic_normal)) .and. (covar .eq. fac_present)) then
            status = fic_write_covar(channel,outfile_ext,cov_recs,
     .                               first_covar_times,last_covar_times,
     .                               output_luns,num_coadds,output_filenames)
         end if
c
c     Close all files for the channel.
c
         if (status .eq. %loc(fic_normal)) then
            status = fic_close(input,short_lun,sci_lun,eng_lun,output_luns,
     .                         report_luns)
         end if

         channel = channel + skip_channel
      end do
c
c     Close reference data sets.
c
      if (status .eq. %loc(fic_normal)) then
         status = fic_reference(jstart,rse_input,rse_name,init,rse,mincoadd,
     .                          grttrans,grtrawwt,grtcoawt,cmdgain,samprate,
     .                          nyquist,cth,basis,reftemps_lun,dtrf_lun,
     .                          gltchpro_lun,etf_lun)
      end if
c
c     Complete summary report.
c
      if (status .eq. %loc(fic_normal)) then
         status = fic_close_summary(covar,deglitch,output,report,mincoadd,cth,
     .                              num_sci_recs,num_groups,num_sci_good,
     .                              num_sci_fail,num_sci_noise,tot_noise,
     .                              num_sci_deglitch,glitch_tot_iter,
     .                              glitch_tot_signal,num_coadds,
     .                              output_filenames,report_filenames)
      end if

      if (status .eq. %loc(fic_normal)) then
         call lib$signal(fic_normal)
      else
         call lib$signal(fic_failure)
      end if

      close(fut_report_lun)

      end
