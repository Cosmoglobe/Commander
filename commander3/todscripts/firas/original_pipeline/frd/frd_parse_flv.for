      integer*4 function frd_parse_flv(channel,scan_mode,galexc,galexc_val,
     .                                 min_ifgs,label,label_val,infile,
     .                                 infile_num,outfile)
c-------------------------------------------------------------------------------
c
c     Purpose: Parse the command line for FRD_FLV.
c
c     Author: S. Brodd, HSTX, 10/95.
c
c     Input: none
c
c     Output: channel       integer*4         Channel, 1-4 = RH-LL.
c             scan_mode     integer*4         Scan mode, 1-6 = SS-FL.
c             galexc        integer*4         Galactic exclusion requested.
c             galexc_val    real*4            Value of galactic exclusion.
c             min_ifgs      integer*4         Minimum number of interferograms.
c             label         integer*4         Label requested.
c             label_val     character*40      Value of label.
c             infile        character*48(20)  Input FIL filenames.
c             infile_num    integer*4         Number of input FIL filenames.
c             outfile       character*48      Output FEX_FLV filename.
c
c     Modifications:
c
c-------------------------------------------------------------------------------
      implicit none
c
c     Include files.
c
      include '(upm_stat_msg)'
      include '(fut_params)'
c
c     Return statuses.
c
      external frd_normal,frd_galexc,frd_noscript,frd_nofilext
      external frd_rmsopen,frd_rmsread,frd_rmsclose,fut_normal
c
c     Functions.
c
      integer*4 upm_present,upm_get_value,upm_get_float,upm_get_longword
      integer*4 fut_get_lun,fut_free_lun
c
c     Output parameters.
c
      integer*4 channel,scan_mode,galexc,min_ifgs,label,infile_num
      real*4 galexc_val
      character*40 label_val
      character*48 infile(20),outfile
c
c     Local variables.
c
      integer*4 status,length,script_file_len,lun,file_num
      character*56 script_file
      namelist / filexts / infile_ext,outfile_ext
      character*20 infile_ext(20)
      character*20 outfile_ext

      frd_parse_flv = %loc(frd_normal)
c
c     Parse channel qualifier.
c
      if (upm_present('channel.rh')) then
         channel = 1
      else if (upm_present('channel.rl')) then
         channel = 2
      else if (upm_present('channel.lh')) then
         channel = 3
      else if (upm_present('channel.ll')) then
         channel = 4
      end if
c
c     Parse scan_mode qualifier.
c
      if (upm_present('scan_mode.ss')) then
         scan_mode = 1
      else if (upm_present('scan_mode.sf')) then
         scan_mode = 2
      else if (upm_present('scan_mode.lf')) then
         scan_mode = 4
      else if (upm_present('scan_mode.fs')) then
         scan_mode = 5
      else if (upm_present('scan_mode.fl')) then
         scan_mode = 6
      end if
c
c     Parse exc_glat qualifier.
c
      status = upm_get_float('exc_glat',galexc_val)
      if (status .eq. upm_absent) then 
         galexc = fac_not_present
         galexc_val = 0.0
      else if ((galexc_val .ge. 0.0)  .or. (galexc_val .le. 90.0)) then
         galexc = fac_present
      else
         frd_parse_flv = %loc(frd_galexc)
         call lib$signal(frd_galexc,%val(1),%val(status))
         return
      end if
c
c     Parse min_ifgs qualifier.
c
      min_ifgs = 3
      status = upm_get_longword('min_ifgs',min_ifgs)
c
c     Parse label qualifier and retrieve value if present.
c
      status = upm_get_value('label',label_val,length)
      if (status .ne. upm_absent) then 
         label = fac_present
         call str$upcase(label_val,label_val)
      else
         label = fac_not_present
      end if
c
c     Parse script qualifier.
c
      status = upm_get_value('script',script_file,script_file_len)
      if (status .eq. upm_absent) then 
         frd_parse_flv = %loc(frd_noscript)
         call lib$signal(frd_noscript)
         return
      else
         call str$upcase(script_file,script_file)
      end if
c
c     Open the script file.
c
      status = fut_get_lun(lun)
      if (status .ne. %loc(fut_normal)) then
         frd_parse_flv = status
         call lib$signal(%val(status))
         return
      end if

      open(unit=lun,file=script_file,access='sequential',status='old',
     .     readonly,iostat=status)

      if (status .ne. 0) then
         frd_parse_flv = %loc(frd_rmsopen)
         call lib$signal(frd_rmsopen,%val(2),script_file(1:script_file_len),
     .                   %val(status))
         return
      end if
c
c     Read the filename extensions from the file and determine input file names.
c
      do file_num = 1,20
         infile_ext(file_num)(1:1) = ' '
      end do

      read(lun,nml=filexts,iostat=status)
      if (status .ne. 0) then
         frd_parse_flv = %loc(frd_rmsread)
         call lib$signal(frd_rmsread,%val(2),script_file(1:script_file_len),
     .                   %val(status))
         return
      end if

      file_num = 1
      do while (file_num .le. 20)
         if (infile_ext(file_num)(1:1) .ne. ' ') then
            infile_num = infile_num + 1
            infile(file_num) = 'CSDR$FIRAS_IN:FIL_SKY_'//
     .                         fac_channel_ids(channel)//
     .                         fac_scan_mode_idsl(scan_mode)//'.'//
     .                         infile_ext(file_num)
            call str$upcase(infile(file_num),infile(file_num))
            call str$trim(infile(file_num),infile(file_num),length)
         end if
         file_num = file_num + 1
      end do
c
c     Determine the output file name.
c
      if (outfile_ext(1:1) .eq. ' ') then
         frd_parse_flv = %loc(frd_nofilext)
         call lib$signal(frd_nofilext,%val(1),script_file(1:script_file_len))
         return
      end if

      outfile = 'CSDR$FIRAS_OUT:FEX_FLV_'//fac_channel_ids(channel)//
     .          fac_scan_mode_idsl(scan_mode)//'.'//outfile_ext
      call str$upcase(outfile,outfile)
      call str$trim(outfile,outfile,length)
c
c     Close the script file.
c
      close(lun,iostat=status)
      if (status .ne. 0) then
         frd_parse_flv = %loc(frd_rmsclose)
         call lib$signal(frd_rmsclose,%val(2),script_file(1:script_file_len),
     .                   %val(status))
         return
      end if

      status = fut_free_lun(lun)
      if (status .ne. %loc(fut_normal)) then
         frd_parse_flv = status
         call lib$signal(%val(status))
         return
      end if

      return
      end
