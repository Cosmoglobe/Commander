      integer*4 function FSS_get_rse ( rsefile, RSE )

C------------------------------------------------------------------------
C    PURPOSE: Return the contents of an RSE file. 
C
C    AUTHOR: D. Bouler, STX, Jan, 1990
C
C    INVOCATION:  status = fss_get_rse (rsefile, RSE)
C
C    INPUT PARAMETERS:         c*14       jstart      start time
C                              C*14       jstop       stop  time
C                              c*128      rsefile     rse file name
C
C    OUTPUT PARAMETERS:        c*128(16)  RSE         Record selection exps
C
C    SUBROUTINES CALLED:       None.
C
C    COMMON VARIABLES USED:    None.
C 
C    INCLUDE FILES:            None.
C
C----------------------------------------------------------------------
c
c                  PDL for FSS_GET_RSE
c                  D. Bouler, Jan, 1990
c
c   Set return status to fss_normal
c
c   call lib$get_lun to get lu_rse
c   if (.not. status) then 
c        call lib$signal(%val(status))
c        set status to error
c   endif
c
c   if (status) then
c       open RSE file
c       if (ios .ne. 0) then
c           signal error
c           set status to error
c       endif
c   endif        
c
c   if (status) then
c       rse_index = 1
c       do while ((ios .eq. 0) .and. (rse_index .le. 16) )
c          read rse(rse_index) from RSE file
c          if (ios .ne. 0) then
c              signal error
c              set status to error
c          else
c              rse_index = rse_index + 1
c          endif
c       end do
c   endif       
c
c   close (lu_rse)
c   call lib$free_lun(lu_rse)
c            
c   return
c   end (pdl)
c
c-----------------------------------------------------------------------------
C
C Changes:
C
C
C----------------------------------------------------------------------
c

      implicit none
c
c Functions
c
      integer*4 lib$get_lun
      integer*4 lib$free_lun
c
c Externals
c
      external fss_normal                  !Everything OK signal
c
c Variables:
c
      integer*4 i,j,k,l                    !Loop counters
      integer*4 ios                        !Input/Output status
      integer*4 status                     !general status
      integer*4 lu_rse                     !LU for rse file
      integer*4 rse_index                  !Which rse is being read from file
      integer*4 start(2)                   !binary start time
      integer*4 stop(2)                    !binary stop  time
      integer*4 len                        !Length of rse file name

      character*14 jstart                  !Start time
      character*14 jstop                   !Stop  time

      character*128 rsefile                !Rse file name
      character*128 rse(16)                !Record selection expression

c
c***  BEGIN CODE  ******************************************************
c
      fss_get_rse = %loc(fss_normal)
      status      = %loc(fss_normal)

      do i = 1,16
         rse(i) = ' '
      enddo

      call str$trim(rsefile, rsefile, len)

      status = lib$get_lun(lu_rse)
      if (.not. status) then 
           call lib$signal(%val(status))
      endif

      if (status) then
          open (lu_rse,status='old',file=rsefile,iostat=ios,readonly)
          if (ios .ne. 0) then
              call errsns(ios,,,,status)
              write(6,*) 'Could not open RSE file ',rsefile(:len)
          endif
       endif

       if (status) then
           rse_index = 1
           do while ((ios .eq. 0) .and. (rse_index .le. 16) )
              read (lu_rse, '(A128)', iostat=ios, end=20) rse(rse_index)
              if (ios .ne. 0) then
                  call errsns(ios,,,,status)
                  write(6,*) 'Error reading from RSE file ',rsefile(:len)
              else
                  rse_index = rse_index + 1
              endif
            end do
20          continue
       endif

       close (lu_rse)
       call lib$free_lun(lu_rse)

       fss_get_rse = status

       return
       end
