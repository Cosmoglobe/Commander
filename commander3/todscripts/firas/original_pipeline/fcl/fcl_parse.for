      integer*4 function fcl_parse(channel,scan_mode,cov_ext,mod_ext,bol_avg,
     .                             cbias_avg,volt_avg,label,label_val,
     .                             command_line)
c----------------------------------------------------------------------------
c
c     Purpose: Parse the command line.
c
c     Author: S. Brodd, HSTX, 12/95, SPR 12291
c
c     Input:
c
c     Output: channel       i*4    Value of channel, 1-4.
c             scan_mode     i*4    Value of scan mode, 1-6=SS-FL (3 excluded).
c             cov_ext       ch*72  Value of covariance matrix file extension.
c             mod_ext       ch*72  Value of model file extension.
c             bol_avg       r*4    Average bolometer temperature in degrees K.
c             cbias_avg     r*4    Average commanded bias in counts.
c             volt_avg      r*4    Average readout voltage in volts.
c             label         l*1    Label present or not.
c             label_val     ch*40  Optional label for output covariance matrix.
c             command_line  ch*79(7)  Fully defaulted command line.
c
c     Modifications:
c
c----------------------------------------------------------------------------
      implicit none
c
c     Include files.
c
      include '(upm_stat_msg)'
c
c     Return statuses.
c
      external fcl_normal
      external fcl_parseerr
c
c     Functions.
c
      integer*4 upm_present
      integer*4 upm_get_value
      integer*4 upm_get_float
c
c     Output parameters.
c
      integer*4 channel,scan_mode
      character*72 cov_ext,mod_ext
      real*4 bol_avg,cbias_avg,volt_avg
      logical*1 label
      character*40 label_val
      character*79 command_line(7)
c
c     Local variables.
c
      integer*4 status,length
      character*12 val

      fcl_parse = %loc(fcl_normal)
c
c     Parse channel qualifier.
c
      command_line(1)(1:12) = 'FCL/CHANNEL='

      if (upm_present('channel.rh')) then
         channel = 1
         command_line(1)(13:14) = 'RH'
      else if (upm_present('channel.rl')) then
         channel = 2
         command_line(1)(13:14) = 'RL'
      else if (upm_present('channel.lh')) then
         channel = 3
         command_line(1)(13:14) = 'LH'
      else if (upm_present('channel.ll')) then
         channel = 4
         command_line(1)(13:14) = 'LL'
      else
         fcl_parse = %loc(fcl_parseerr)
         return
      end if
c
c     Parse scan_mode qualifier.
c
      command_line(1)(15:25) = '/SCAN_MODE='

      if (upm_present('scan_mode.ss')) then
         scan_mode = 1
         command_line(1)(26:27) = 'SS'
      else if (upm_present('scan_mode.sf')) then
         scan_mode = 2
         command_line(1)(26:27) = 'SF'
      else if (upm_present('scan_mode.lf')) then
         scan_mode = 4
         command_line(1)(26:27) = 'LF'
      else if (upm_present('scan_mode.fs')) then
         scan_mode = 5
         command_line(1)(26:27) = 'FS'
      else if (upm_present('scan_mode.fl')) then
         scan_mode = 6
         command_line(1)(26:27) = 'FL'
      else
         fcl_parse = %loc(fcl_parseerr)
         return
      end if
c
c     Parse cov_ext qualifier.
c
      status = upm_get_value('cov_ext',cov_ext,length)
      if (status .eq. upm_absent) then 
         fcl_parse = %loc(fcl_parseerr)
         return
      else
         call str$upcase(cov_ext,cov_ext)
         command_line(2)(1:9+length) = '/COV_EXT='//cov_ext(1:length)
c
c     Return only the file extension for use in opening input files.
c
         length = index(cov_ext,'.')
         cov_ext = cov_ext(length+1:)
      end if
c
c     Parse mod_ext qualifier.
c
      status = upm_get_value('mod_ext',mod_ext,length)
      if (status .eq. upm_absent) then 
         fcl_parse = %loc(fcl_parseerr)
         return
      else
         call str$upcase(mod_ext,mod_ext)
         command_line(3)(1:9+length) = '/MOD_EXT='//mod_ext(1:length)
c
c     Return only the file extension for use in opening input files.
c
         length = index(mod_ext,'.')
         mod_ext = mod_ext(length+1:)
      end if
c
c     Parse bol_avg qualifier.
c
      status = upm_get_float('bol_avg',bol_avg)
      if (status .eq. upm_absent) then 
         fcl_parse = %loc(fcl_parseerr)
         return
      else
         write(val,10) bol_avg
         command_line(4)(1:21) = '/BOL_AVG='//val
      end if
c
c     Parse cbias_avg qualifier.
c
      status = upm_get_float('cbias_avg',cbias_avg)
      if (status .eq. upm_absent) then 
         fcl_parse = %loc(fcl_parseerr)
         return
      else
         write(val,10) cbias_avg
         command_line(5)(1:23) = '/CBIAS_AVG='//val
      end if
c
c     Parse volt_avg qualifier.
c
      status = upm_get_float('volt_avg',volt_avg)
      if (status .eq. upm_absent) then 
         fcl_parse = %loc(fcl_parseerr)
         return
      else
         write(val,10) volt_avg
         command_line(6)(1:22) = '/VOLT_AVG='//val
      end if
c
c     Parse label qualifier and retrieve value if present.
c
      label = .false.
      if (upm_present('label') .eq. upm_pres) then
         status = upm_get_value('label',label_val,length)
         if (status .ne. upm_absent) then 
            label = .true.
            call str$upcase(label_val,label_val)
            command_line(7)(1:7+length) = '/LABEL='//label_val(1:length)
         end if
      end if

 10   format(f12.8)
 
      return
      end
