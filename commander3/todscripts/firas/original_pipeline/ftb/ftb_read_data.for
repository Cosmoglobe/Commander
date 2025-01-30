       integer*4 function ftb_read_data(word31, time_flags, num_recs, num_pts, lun)
c
c      Purpose:
c
c          To read the desired data and form the time flags array for the 
c          purpose of identifying time gaps.
c
c      Written by:
c
c          J. W. Durachta
c              ARC
c          August, 1988
c
c      Modified by:
c
c	   Shirley M. Read
c	       STX
c	   August, 1988
c	   Reason:  Output messages via Lib$Signal.
c
c---------------------------------------------------------------------------
c

       implicit none

       include 'ct$library:ctuser.inc'

       dictionary 'nfs_anc'
       record/nfs_anc/ word31(1)

       external  ftb_ctread

       integer*2 ct_stat(20)

       integer*4 status
       integer*4 num_recs
       integer*4 lun
       integer*4 sec80/800000000/
       integer*4 lib$subx
       integer*4 lib$ediv
       integer*4 lib$emul
       integer*4 ptr, num_pts
       integer*4 diff(2)
       integer*4 diff_mj(2)
       integer*4 diff_rec
       integer*4 quo, rem
       integer*4 i

       logical*4 time_flags(1)

       common /data_info/ diff_mj, diff_rec

       ftb_read_data = 1
       ct_stat(1) = ctp_normal
       num_recs = 0

!   Read the data.

       do while(ct_stat(1) .eq. ctp_normal)
         num_recs = num_recs + 1
         call ct_read_arcv(, lun, word31(num_recs), ct_stat)
       end do

       num_recs = num_recs - 1

       if(ct_stat(1) .ne. ctp_endoffile)then
	  call lib$signal(ftb_ctread,%val(1),%val(ct_stat(1)))

!         type*, 'ERROR: Read error; CT status = ', ct_stat(1), '.
!    & WORD31 Aborting.'
         ftb_read_data = 0

       else

!   Find the time difference between rec(i+1) and rec(i). If > 1 rec, calculate
!   the time gap. Set time flags accordingly.

         num_pts = 1
         ptr = 1
         time_flags(1) = .true.
         status = lib$subx(word31(1).gmt_mjf2, word31(1).ct_head.time, diff_mj)
         diff_rec = 2 * diff_mj(1)

         do while(ptr .lt. num_recs)

           status = lib$subx(word31(ptr+1).ct_head.time,
     &                                word31(ptr).ct_head.time, diff)

           if(diff(1) .gt. sec80 .or. diff(1) .lt. 0
     &           .or. diff(2) .gt. 0)then

             status = lib$ediv(diff_rec, diff, quo, rem)
             if(rem .le. diff_mj(1))then

               do i = 1, quo - 1
                 time_flags(num_pts+i) = .false.
               end do

             else

               do i = 1, quo - 2
                 time_flags(num_pts+i) = .false.
               end do

             endif

             num_pts = num_pts + i
             time_flags(num_pts) = .true.

           else

             num_pts = num_pts + 1
             time_flags(num_pts) = .true.

           endif

           ptr = ptr + 1

         end do

       endif

       return
       end
