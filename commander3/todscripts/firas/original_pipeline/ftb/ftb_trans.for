      program ftb_trans

c------------------------------------------------------------------------------
c
c     Purpose:  To rewrite output of FFI (FIRAS_Fish_Input) to ASCII files 
c               so that they may be transferred via FTP to Ultrix machines 
c               in the most error-free manner.
c
c     Author:  David Cottingham, USRA, 4/5/93, SER 10764
c
c     Modifications:
c
c------------------------------------------------------------------------------

      implicit none
c
c     Include file.
c
      include '(upm_stat_msg)'
c
c     External declarations.
c
      external fut_error,ftb_normal,ftb_aberr
      external ftb_rmsopen,ftb_rmsread,ftb_rmswrite
c
c     Functions.
c
      integer*4 fut_error
      integer*4 cut_register_version,cut_display_banner
      integer*4 upm_get_value,lib$get_lun
c
c     Local variables.
c 
      integer*4  rstatus,num_vol/80/,lun_out/6/,iostatus
      character*6 version
      parameter (version='10.6')

      logical*1 status,eof

      integer*4 inlun,outlun,len
      character*80 infile,outfile

      integer*4 nifgs,gain
      real*4 time,temps(7),temp_sigs(7)
      real*8 volt,bias,sigma(175)
      complex*16 spec(175)
c
c     Establish the condition handler and display the program banner.
c
      call lib$establish(fut_error)
      rstatus = cut_register_version(version)
      rstatus = cut_display_banner(lun_out,num_vol,'FIRAS Facility FTB_Trans')
c
c     Get the file names and open the files.
c
      rstatus = upm_get_value('INFILE',infile,len)
      rstatus = lib$get_lun(inlun)

      status = .true.

      open(unit=inlun,file=infile,status='old',readonly,form='unformatted',
     .     recordtype='fixed',recl=1071,iostat=iostatus)
      if (iostatus .ne. 0) then
         status = .false.
         call str$trim(infile,infile,len)
         call lib$signal(ftb_rmsopen,%val(2),infile(:len),%val(iostatus))
      end if

      if (status) then
         rstatus = upm_get_value('OUTFILE',outfile,len)
         rstatus = lib$get_lun(outlun)
         open(unit=outlun,file=outfile,status='unknown',iostat=iostatus)
         if (iostatus .ne. 0) then
            status = .false.
            call str$trim(outfile,outfile,len)
            call lib$signal(ftb_rmsopen,%val(2),outfile(:len),%val(iostatus))
         end if
      end if

      eof = .false.
      do while (status .and. (.not. eof))
         read (inlun,iostat=iostatus) 
     .         nifgs,time,temps,temp_sigs,volt,bias,gain,spec,sigma
         if (iostatus .lt. 0) then
            eof = .true.
         else if (iostatus .gt. 0) then
            status = .false.
            call str$trim(infile,infile,len)
            call lib$signal(ftb_rmsread,%val(2),infile(:len),%val(iostatus))
         else
            write(outlun,*,iostat=iostatus) 
     .            nifgs,time,temps,temp_sigs,volt,bias,gain,spec,sigma
            if (iostatus .ne. 0) then
               status = .false.
               call str$trim(outfile,outfile,len)
               call lib$signal(ftb_rmswrite,%val(2),outfile(:len),%val(iostatus))
            end if
         end if
      end do

      close(inlun)
      close(outlun)

      if (status) then
        call lib$signal(ftb_normal)
      else
        call lib$signal(ftb_aberr)
      endif

      end
