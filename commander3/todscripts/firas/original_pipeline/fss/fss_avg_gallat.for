      integer*4  function fss_avg_gallat  (sky_buff, rec_num, num_recs, avg)

C------------------------------------------------------------------------
C    PURPOSE: Calculate average galactic latitude of the coadd group.
C
C    AUTHOR: L. Rosen, HSTX, 13 February 1995
C
C    INVOCATION:  status = fss_avg_gallat  (sky_buff, rec_num, num_recs, avg)
C
C    INPUT PARAMETERS:
C     integer*4      rec_num                   ! Record number of 1st of group
C     record/fss_sssky/ sky_buff(fss_max_buff) ! Array of short sky recs
C     integer*4      num_recs                  ! Number of records in array
C    OUTPUT PARAMETERS:
C     real*4      avg                          ! Average galactic latitude
C
C    INCLUDE FILES:  fss_include
C
C    METHOD:
C       rec_num indicates the first record of the coadd group in the sky_buff
C       array.  Sum the galactic latitudes of each record minus the first
C       record, after the first.  Divide the sum by the number of records
C       in the group - 1.  Add that value to the first one.  Galactic
C       latitudes are stored as 2 byte integers of latitude in 1e-4 radians.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
c
c Include files
c
      include       '(fss_include)'
c
c Passed parameters
c
      integer*4      rec_num                   ! Record number of 1st of group
      dictionary 'fss_sssky'
      record/fss_sssky/ sky_buff(fss_max_buff)       ! Array of short sky recs
      integer*4      num_recs                  ! Number of records in array
      real*4         avg
c
c External
c
      external       fss_normal
c
c Local
c
      integer*2      first
      integer*4      sum
      logical*1      more
      integer*4      n
c*
c******  Begin code *******************************************
c*
      fss_avg_gallat = %loc(fss_normal)

      first = sky_buff (rec_num).galactic_latitude
      sum = 0
      more = .true.
      n = rec_num + 1
      do while ( n .le. num_recs .and. more )
         if (sky_buff(n).transition_flag .eq. 0) then
            sum = sum + (sky_buff(n).galactic_latitude - first)
            n = n + 1
         else
            more = .false.
         endif
      enddo
      avg = real (first) + real (sum) / real (n - rec_num)
      return
      end
