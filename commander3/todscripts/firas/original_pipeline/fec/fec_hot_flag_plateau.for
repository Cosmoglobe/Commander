      integer*4 function Fec_hot_flag_plateau ( start_rec, end_rec,
     r                                       index,
     R                                       Array)
c
c     By H. Wang, STX,  22-May-1990
c
c     Description:
c         fec_hot_flag_plateau flag the hot spot plateau within the
C          grt command  plateau array.
c
c     Format of Call:
c         Ret_Code = fec_hot_flag_plateau (start_rec,end_rec,
c                                       Index,array)
c     Input Parameters:
c         start_rec        = the start point of the hot spot plateau 
c         end_rec          = the end   point of the hot spot plateau 
c         index            = the size of the plateau array
c         Array            = the plateau array
c
c     Output Parameters:
c         index            = the size of the plateau array
c         Array            = the merged plateau array
c         Ret_Code         = TRUE if a change command occurred, otherwise FALSE
c                            integer*4
c
c     Modules Called:
c         None
c
c     Include Files: 
c         FEC.PAR          = parameters used by FEC (FEC_TRUE, FEC_FALSE)
c         FIRAS_REC_DEF    = definition of engineering record
c
c     Changes:
c----------------------------------------------------------------------------
c     PDL:
C
c           
C         do while (hot spot starting point not found in the grt command plateau
c                   array)
c           if the grt command plateau point(record number) is .gt.
c              hot spot starting point(record number)
C           then
c              hot spot starting point has been found in the grt
c              command plateau array
c           endif
c         endo        
c           
C         do while (hot spot ending point not found in the grt command plateau
c                   array)
c           if the grt command plateau point(record number) is .gt.
c              hot spot ending point(record number)
C           then
c              hot spot ending point has been found in the grt
c              command plateau array
c           endif
c         endo        
c          flag the hot spot starting point in the grt commanded plateau 
c          array
c          flag the hot spot ending point in the grt commanded plateau 
c          array
c         return
c         end 
c------------------------------------------------------------------------
c     ******* BEGIN *******
c
      implicit none
c
c     Include files:
c
      include '(FEC_INC)'             ! contains TRUE/FALSE values to return
      include '(FEC_MSG)'
c
c     Formal Parameters:
c
      integer*4 start_rec
      integer*4 end_rec  
      integer*4 Index
      integer*4 Array(0:max_plateaus)
      integer*4 Array_temp(0:max_plateaus)
      integer*4 i,k, j, start_temp, end_temp
      logical*1 found
      logical*1 start_dup
      logical*1 end_dup
c
c     *** BEGIN ***
c
          start_dup = .false.
          end_dup = .false.
          k = 1 
          found = .false.
         do while ((.not. found) .and. (k .le. index))
           if ((array(k) .ge. start_rec) .and. (array(k) .ne. 0)) then
             found = .true.
             start_temp = k
           else
             k = k + 1
           endif
          enddo
          k = 1 
          found = .false.
         do while ((.not. found) .and. (k .le. index)) 
           if ((array(k) .ge. end_rec) .and. (array(k) .ne. 0)) then
             found = .true.
             end_temp = k
           else
             k = k + 1
           endif
         enddo
         if (array(start_temp) .eq. start_rec) start_dup=.true. 
         if (array(end_temp) .eq. end_rec) end_dup=.true.
         do i = 0, start_temp - 1
           array_temp(i) = array(i) 
         enddo
	   array_temp(start_temp)= start_rec
           array_temp(start_temp +1) = 0
           array_temp(start_temp + 2) = end_rec
           k = end_temp
           if (end_dup) k = k + 1
           j = start_temp + 2 + 1  
           do while ( k .le. index)
             array_temp(j) = array(k)
             j = j + 1
             k = k + 1
           enddo
            index = j - 1
           do i = 0, index
            array(i) = array_temp(i)
           enddo
       fec_hot_flag_plateau = %loc(FEC_NORMAL)
      return 
      end !fec_hot_flag_plateau
