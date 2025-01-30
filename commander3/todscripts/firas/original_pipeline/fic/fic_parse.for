      integer*4 function fic_parse(input,start_channel,stop_channel,
     .                             skip_channel,scan_mode_input,file,file_ext,
     .                             jstart,jstop,rse_input,rse_name,covar,
     .                             covar_new,qual_inst,qual_att,single,
     .                             inst_state,pixel_no,sec_template_input,
     .                             deglitch,shape_check,coadd,baseline_sub,
     .                             output,outfile,outfile_ext,report,
     .                             report_detail,report_name,command_line)
c----------------------------------------------------------------------------
c
c     Purpose: Parse the command line.
c
c     Author: S. Alexander, STX, 9/91, SER 7985
c
c     Input:
c
c     Output: input              i*4  Type of input: sky, calibration, or raw.
c             start_channel      i*4  Value of first channel to process, 1-4.
c             stop_channel       i*4  Value of last channel to process, 1-4.
c             skip_channel       i*4  Skip nonsequential channels.
c             scan_mode_input    i*4  Input scan mode, 1-4=SS-LF; 0 if none.
c             file               i*4  Input filename specified on command line.
c             file_ext           ch*72  Value of input filename extension on
c                                       command line.
c             jstart             ch*14  Value of start time.
c             jstop              ch*14  Value of stop time.
c             rse_input          i*4  Rse specified on command line.
c             rse_name           ch*72  Rse filename.
c             covar              i*4  Covariance matrix to be written.
c             covar_new          i*4  Covariance matrix to be initialized.
c             qual_inst          i*4  Maximum acceptable data quality.
c             qual_att           i*4  Maximum acceptable attitude quality.
c             single             i*4  Process single ifgs.
c             inst_state         i*4  Perform state consistency checking.
c             pixel_no           i*4  Consistency check the pixel number.
c             sec_template_input i*4  Perform secondary template subtraction.
c             deglitch           i*4  Perform deglitching.
c             shape_check        i*4  Perform shape consistency checking.
c             coadd              i*4  Perform coaddition.
c             baseline_sub       i*4  Perform baseline subtraction.
c             output             i*4(6) Selection of output types. 
c             outfile            i*4  Output filename specified on command line.
c             outfile_ext        ch*72  Value of output filename extension on
c                                       command line.
c             report             i*4  Write consistency and deglitch reports.
c             report_detail      i*4  Detailed report option.
c             report_name        ch*72  Consistency and deglitch report name on
c                                       command line.
c             command_line       ch*79(9)  Fully defaulted command line.
c
c     Modifications:
c
c----------------------------------------------------------------------------

      implicit none
c
c     Include files.
c
      include '(fut_params)'
      include '(upm_stat_msg)'
c
c     Return statuses.
c
      external fic_normal
      external fic_invtime
c
c     Functions.
c
      integer*4 upm_present
      integer*4 upm_get_value
      logical*2 time_lt
      integer*4 upm_get_longword
c
c     Output parameters.
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
c     Local variables.
c
      integer*4 line_pos,status,length
      integer*4 start(2),stop(2)
      character*2 qual_val(2)
      logical*1 other_output

      fic_parse = %loc(fic_normal)
c
c     Parse input qualifier.
c
      command_line(1)(1:10) = 'FIC/INPUT='

      if (upm_present('input.calibration')) then
         input = fac_cal
         command_line(1)(11:13) = 'CAL'
      else if (upm_present('input.raw')) then
         input = fac_raw
         command_line(1)(11:13) = 'RAW'
      else
         input = fac_sky
         command_line(1)(11:13) = 'SKY'
      end if
c
c     Parse channel qualifier.
c
      command_line(1)(14:22) = '/CHANNEL='

      if (upm_present('channel.rh')) then
         start_channel = 1
         stop_channel = 1
         skip_channel = 1
         command_line(1)(23:24) = 'RH'
         line_pos = 25
      else if (upm_present('channel.rl')) then
         start_channel = 2
         stop_channel = 2
         skip_channel = 1
         command_line(1)(23:24) = 'RL'
         line_pos = 25
      else if (upm_present('channel.lh')) then
         start_channel = 3
         stop_channel = 3
         skip_channel = 1
         command_line(1)(23:24) = 'LH'
         line_pos = 25
      else if (upm_present('channel.ll')) then
         start_channel = 4
         stop_channel = 4
         skip_channel = 1
         command_line(1)(23:24) = 'LL'
         line_pos = 25
      else if (upm_present('channel.right')) then
         start_channel = 1
         stop_channel = 2
         skip_channel = 1
         command_line(1)(23:27) = 'RIGHT'
         line_pos = 28
      else if (upm_present('channel.left')) then
         start_channel = 3
         stop_channel = 4
         skip_channel = 1
         command_line(1)(23:26) = 'LEFT'
         line_pos = 27
      else if (upm_present('channel.high')) then
         start_channel = 1
         stop_channel = 3
         skip_channel = 2
         command_line(1)(23:26) = 'HIGH'
         line_pos = 27
      else if (upm_present('channel.low')) then
         start_channel = 2
         stop_channel = 4
         skip_channel = 2
         command_line(1)(23:25) = 'LOW'
         line_pos = 26
      else
         start_channel = 1
         stop_channel = 4
         skip_channel = 1
         command_line(1)(23:25) = 'ALL'
         line_pos = 26
      end if
c
c     Parse scan_mode qualifier.
c
      command_line(1)(line_pos:line_pos+10) = '/SCAN_MODE='
      line_pos = line_pos+11

      if (upm_present('scan_mode.ss')) then
         scan_mode_input = 1
         command_line(1)(line_pos:line_pos+1) = 'SS'
      else if (upm_present('scan_mode.sf')) then
         scan_mode_input = 2
         command_line(1)(line_pos:line_pos+1) = 'SF'
      else if (upm_present('scan_mode.ls')) then
         scan_mode_input = 3
         command_line(1)(line_pos:line_pos+1) = 'LS'
      else if (upm_present('scan_mode.lf')) then
         scan_mode_input = 4
         command_line(1)(line_pos:line_pos+1) = 'LF'
      else
         scan_mode_input = 0
         command_line(1)(line_pos:line_pos+2) = 'ALL'
      end if
c
c     Parse file qualifier and retrieve input filename if present.
c
      if (upm_present('file') .eq. upm_pres) then
         file = fac_present
         status = upm_get_value('file',file_ext,length)
         call str$upcase(file_ext,file_ext)
         command_line(2)(1:6+length) = '/FILE='//file_ext(1:length)
c
c     Return only the file extension for use in opening input files.
c
         length = index(file_ext,'.')
         file_ext = file_ext(length+1:)
      else
c
c     Get input or defaulted jstart and jstop.
c
         file = fac_not_present
         status = upm_get_value('jstart',jstart,length)
         if (status .eq. upm_absent) then
            jstart = fac_jstart_default
         else
            jstart = jstart(1:length)//fac_jstart_default(length+1:)
         end if

         status = upm_get_value('jstop',jstop,length)
         if (status .eq. upm_absent) then
            jstop = fac_jstop_default
         else
            jstop = jstop(1:length)//fac_jstop_default(length+1:)
         end if

         command_line(2)(1:50) = '/NOFILE/JSTART='//jstart//'/JSTOP='//jstop
c
c     Test that jstart is less than jstop; return if error.
c           
         call ct_gmt_to_binary(jstart,start)
         call ct_gmt_to_binary(jstop,stop)

         if (time_lt(stop,start)) then
            fic_parse = %loc(fic_invtime)
            return
         end if
      end if
c
c     Parse rse qualifer and retrieve rse filename if present.
c
      if (upm_present('rse') .eq. upm_pres) then
         rse_input = fac_present
         status = upm_get_value('rse',rse_name,length)
         call str$upcase(rse_name,rse_name)
         command_line(3)(1:5+length) = '/RSE='//rse_name(1:length)
      else
         rse_input = fac_not_present
         command_line(3)(1:6) = '/NORSE'
      end if
c
c     Parse covar qualifer.
c
      if (upm_present('covar.new') .eq. upm_pres) then
         covar = fac_present
         covar_new = fac_present
         command_line(4)(1:10) = '/COVAR=NEW'
         line_pos = 11
      else if (upm_present('covar') .eq. upm_pres) then
         covar = fac_present
         covar_new = fac_not_present
         command_line(4)(1:13) = '/COVAR=UPDATE'
         line_pos = 14
      else
         covar = fac_not_present
         covar_new = fac_not_present
         command_line(4)(1:8) = '/NOCOVAR'
         line_pos = 9
      end if
c
c     Parse quality qualifier, returning instrument and attitude quality.
c
      qual_inst = 3
      qual_att = 3
      status = upm_get_longword('quality.instrument',qual_inst)
      status = upm_get_longword('quality.attitude',qual_att)
      write(qual_val(1),10) qual_inst
      write(qual_val(2),10) qual_att
 10   format(i2.2)
      command_line(4)(line_pos:line_pos+35) = '/QUALITY=(INSTRUMENT='//
     .             qual_val(1)//',ATTITUDE='//qual_val(2)//')'
      line_pos = line_pos + 36
c
c     Parse single qualifier.
c
      if (upm_present('single') .eq. upm_pres) then
         single = fac_present
         command_line(4)(line_pos:line_pos+6) = '/SINGLE'
      else
         single = fac_not_present
         command_line(4)(line_pos:line_pos+8) = '/NOSINGLE'
      end if
c
c     Parse inst_state qualifier.
c
      if ((upm_present('inst_state') .eq. upm_negated) .or. 
     .    (single .eq. fac_present)) then
         inst_state = fac_not_present
         command_line(5)(1:13) = '/NOINST_STATE'
         line_pos = 14
      else
         inst_state = fac_present
         command_line(5)(1:11) = '/INST_STATE'
         line_pos = 12
      end if
c
c     Parse pixel_no qualifier.
c
      if ((upm_present('pixel_no') .eq. upm_negated) .or.
     .    (single .eq. fac_present)) then
         pixel_no = fac_not_present
         command_line(5)(line_pos:line_pos+10) = '/NOPIXEL_NO'
         line_pos = line_pos+11
      else
         pixel_no = fac_present
         command_line(5)(line_pos:line_pos+8) = '/PIXEL_NO'
         line_pos = line_pos+9
      end if
c
c     Parse sec_template qualifier.
c
      if ((upm_present('sec_template') .eq. upm_negated) .or.
     .    (single .eq. fac_present)) then
         sec_template_input = fac_not_present
         command_line(5)(line_pos:line_pos+14) = '/NOSEC_TEMPLATE'
         line_pos = line_pos+15
      else
         sec_template_input = fac_present
         command_line(5)(line_pos:line_pos+12) = '/SEC_TEMPLATE'
         line_pos = line_pos+13
      end if
c
c     Parse deglitch qualifier.
c
      if ((upm_present('deglitch') .eq. upm_negated) .or.
     .    (single .eq. fac_present)) then
         deglitch = fac_not_present
         command_line(5)(line_pos:line_pos+10) = '/NODEGLITCH'
         line_pos = line_pos+11
      else
         deglitch = fac_present
         command_line(5)(line_pos:line_pos+8) = '/DEGLITCH'
         line_pos = line_pos+9
      end if
c
c     Parse shape_check qualifier.
c
      if ((upm_present('shape_check') .eq. upm_negated) .or.
     .    (single .eq. fac_present)) then
         shape_check = fac_not_present
         command_line(6)(1:14) = '/NOSHAPE_CHECK'
         line_pos = 15
      else
         shape_check = fac_present
         command_line(6)(1:12) = '/SHAPE_CHECK'
         line_pos = 13
      end if
c
c     Parse coadd qualifier.
c
      if (upm_present('coadd') .eq. upm_negated) then
         coadd = fac_not_present
         command_line(6)(line_pos:line_pos+7) = '/NOCOADD'
         line_pos = line_pos+8
      else
         coadd = fac_present
         command_line(6)(line_pos:line_pos+5) = '/COADD'
         line_pos = line_pos+6
      end if
c
c     Parse baseline_sub qualifier.
c
      if ((upm_present('baseline_sub') .eq. upm_negated) .or.
     .    (coadd .eq. fac_not_present)) then
         baseline_sub = fac_not_present
         command_line(6)(line_pos:line_pos+14) = '/NOBASELINE_SUB'
      else
         baseline_sub = fac_present
         command_line(6)(line_pos:line_pos+12) = '/BASELINE_SUB'
      end if
c
c     Parse output qualifier.  Initially set all outputs to not present.
c
      output(fac_out_template) = fac_not_present
      output(fac_out_sec_template) = fac_not_present
      output(fac_out_deglitch) = fac_not_present
      output(fac_out_con_check) = fac_not_present
      output(fac_out_coadd) = fac_not_present
      output(fac_out_baseline_sub) = fac_not_present

      if (upm_present('output') .eq. upm_negated) then
         command_line(7)(1:9) = '/NOOUTPUT'
      else
c
c     Check for presence of each output individually, noting if any are found.
c
         other_output = .false.
         command_line(7)(1:9) = '/OUTPUT=('
         line_pos = 10

         if (upm_present('output.template')) then
            output(fac_out_template) = fac_present
            other_output = .true.
            command_line(7)(line_pos:line_pos+8) = 'TEMPLATE,'
            line_pos = line_pos+9
         end if
    
         if (upm_present('output.sec_template')) then
            output(fac_out_sec_template) = fac_present
            other_output = .true.
            command_line(7)(line_pos:line_pos+12) = 'SEC_TEMPLATE,'
            line_pos = line_pos+13
         end if
   
         if (upm_present('output.deglitch')) then
            output(fac_out_deglitch) = fac_present
            other_output = .true.
            command_line(7)(line_pos:line_pos+8) = 'DEGLITCH,'
            line_pos = line_pos+9
         end if
   
         if (upm_present('output.con_check')) then
            output(fac_out_con_check) = fac_present
            other_output = .true.
            command_line(7)(line_pos:line_pos+9) = 'CON_CHECK,'
            line_pos = line_pos+10
         end if
   
         if (upm_present('output.coadd')) then
            output(fac_out_coadd) = fac_present
            other_output = .true.
            command_line(7)(line_pos:line_pos+5) = 'COADD,'
            line_pos = line_pos+6
         end if
c
c     Baseline_sub output is the default.
c
         if ((upm_present('output.baseline_sub') .or. (.not. other_output))
     .        .and. (baseline_sub .eq. fac_present)) then
            output(fac_out_baseline_sub) = fac_present
            other_output = .true.
            command_line(7)(line_pos:line_pos+12) = 'BASELINE_SUB,'
            line_pos = line_pos+13
         else if (.not. other_output) then
            command_line(7)(1:9) = '/NOOUTPUT'
         end if

         if (other_output) then
            command_line(7)(line_pos-1:line_pos-1) = ')'
         end if
      end if
c
c     Parse outfile qualifier and retrieve output filename if present.
c
      if (upm_present('outfile') .eq. upm_pres) then
         outfile = fac_present
         status = upm_get_value('outfile',outfile_ext,length)
         call str$upcase(outfile_ext,outfile_ext)
         command_line(8)(1:9+length) = '/OUTFILE='//outfile_ext(1:length)
c
c     Return only the filename extension for use in opening output files.
c
         length = index(outfile_ext,'.')
         outfile_ext = outfile_ext(length+1:)
      else
         outfile = fac_not_present
         command_line(8)(1:10) = '/NOOUTFILE'
      end if
c
c     Parse report qualifier and retrieve report file name if present.
c
      if (upm_present('report.detail') .eq. upm_pres) then
         report = fac_present
         report_detail = fac_present
         status = upm_get_value('report.detail',report_name,length)
         if (status .eq. upm_absent) then 
            command_line(9)(1:14) = '/REPORT=DETAIL'
         else
            call str$upcase(report_name,report_name)
            command_line(9)(1:15+length) = 
     .                     '/REPORT=DETAIL='//report_name(1:length)
         end if
      else if (upm_present('report') .eq. upm_pres) then
         report = fac_present
         report_detail = fac_not_present
         status = upm_get_value('report.summary',report_name,length)
         if (status .eq. upm_absent) then 
            command_line(9)(1:15) = '/REPORT=SUMMARY'
         else
            call str$upcase(report_name,report_name)
            command_line(9)(1:16+length) =
     .         '/REPORT=SUMMARY='//report_name(1:length)
         end if
      else
         report = fac_not_present
         report_detail = fac_not_present
         report_name = ' '
         command_line(9)(1:9) = '/NOREPORT'
      end if
 
      return
      end
