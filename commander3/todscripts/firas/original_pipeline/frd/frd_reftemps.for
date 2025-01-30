      program frd_reftemps

c------------------------------------------------------------------------
c
c     Purpose: Put data from the FEX_REFTEMPS.TXT reference text file into 
c              the FEX_REFTEMPS.DAT binary file.
c
c     Input Data: FEX_REFTEMPS.TXT
c
c     Output Data: FEX_REFTEMPS.DAT
c
c     Author: S. Alexander
c             Hughes STX
c             August 6, 1992, SER 7985
c
c     Modifications:
c
c----------------------------------------------------------------------

      implicit none

      include '($ssdef)'

      integer*4 lib$get_lun,lib$free_lun

      integer*4 i,j,k,ios,status,in_lun,out_lun

      character*16 in_file /'FEX_REFTEMPS.TXT'/
      character*16 out_file /'FEX_REFTEMPS.DAT'/

      real*4 prim_temp(41),sec_temp(3)

      dictionary 'FEX_REFTEMPS'
      record /FEX_REFTEMPS/ reftemps

      external      frd_normal
      external      frd_rmsopen
      external      frd_rmsread
      external      frd_rmswrite
      external      frd_rmsclose
c
c     Open the reference template text file.
c
      status = lib$get_lun(in_lun)
      if (status .ne. ss$_normal) then
         call lib$signal(%val(status))
      end if

      open (unit=in_lun,name=in_file,status='old',iostat=ios,readonly,shared)

      if (ios .eq. 0) then
c
c     Open the reference template data file.
c
         status = lib$get_lun(out_lun)
         if (status .ne. ss$_normal) then
            call lib$signal(%val(status))
         end if
       
         open (unit=out_lun,name=out_file,status='new',form='unformatted',
     .         access='direct',recl=64,iostat=ios)

         if (ios .eq. 0) then
c
c     Read reference templates.
c
            do i = 1,32
               do j = 1,10
                  read (in_lun,*) (reftemps.prim_temp((j-1)*4+k),k=1,4)
               end do
               read (in_lun,*) reftemps.prim_temp(41)
               read (in_lun,*) (reftemps.sec_temp(k),k=1,3)
c
c     Write reference templates.
c
               write (out_lun,rec=i,iostat=ios) reftemps
               if (ios .ne. 0) then
                  status = %loc(frd_rmswrite)
                  call lib$signal(frd_rmswrite,%val(2),out_file,%val(ios))
               end if
            end do
c
c     Close input and output files.
c
            close (in_lun,iostat=ios)
            if (ios .ne. 0) then
               status = %loc(frd_rmsclose)
               call lib$signal(frd_rmsclose,%val(2),in_file,%val(ios))
            end if

            close (out_lun,iostat=ios)
            if (ios .ne. 0) then
               status = %loc(frd_rmsclose)
               call lib$signal(frd_rmsclose,%val(2),out_file,%val(ios))
            end if

         else
            status = %loc(frd_rmsopen)
            call lib$signal(frd_rmsopen,%val(2),out_file,%val(ios))
         end if

      else
         status = %loc(frd_rmsopen)
         call lib$signal(frd_rmsopen,%val(2),in_file,%val(ios))
      end if

      status = lib$free_lun(in_lun)
      if (status .ne. ss$_normal) then
         call lib$signal(%val(status))
      end if

      status = lib$free_lun(out_lun)
      if (status .ne. ss$_normal) then
         call lib$signal(%val(status))
      end if


      if (status .eq. ss$_normal) then
         call lib$signal(frd_normal)
      else
         call lib$signal(%val(ss$_abort))
      end if

      stop
      end
